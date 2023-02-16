#
# (C) Copyright IBM 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

from copy import copy, deepcopy
from dataclasses import dataclass
from typing import Dict, List, Set
from rustworkx import PyGraph
from qiskit_qec.circuits.repetition_code import ArcCircuit, RepetitionCodeCircuit
from qiskit_qec.circuits import SurfaceCodeCircuit

from qiskit_qec.decoders.decoding_graph import DecodingGraph


@dataclass
class SpanningForest:
    vertices: Dict[int, List[int]]
    edges: List[int]


@dataclass
class BoundaryEdge:
    index: int
    cluster_vertex: int
    neighbour_vertex: int
    data: Dict[str, object]

    def reverse(self):
        return BoundaryEdge(
            index=self.index,
            cluster_vertex=self.neighbour_vertex,
            neighbour_vertex=self.cluster_vertex,
            data=self.data
        )


@dataclass
class UnionFindDecoderCluster:
    boundary: List[BoundaryEdge]
    atypical_nodes: Set[int]
    fully_grown_edges: Set[int]
    size: int


@dataclass
class FusionEntry:
    u: int
    v: int
    connecting_edge: BoundaryEdge


class UnionFindDecoder:
    """
    Decoder based on growing clusters around syndrome errors to 
    "convert" them into erasure errors, which can be corrected easily,
    by the peeling decoder in case of the surface code, or by checking for
    interference with the boundary in case of an abritrary ARC.

    TODO: Add weights to edges of graph according to Huang et al (see. arXiv:2004.04693, section III)

    See arXiv:1709.06218v3 for more details.
    """

    def __init__(
        self,
        code_circuit: SurfaceCodeCircuit | RepetitionCodeCircuit | ArcCircuit,
        logical: str
    ) -> None:
        self.code_circuit = code_circuit
        self.decoding_graph = DecodingGraph(
            code_circuit,
        )
        self.logical = logical

    def process(self, string: str):
        """
        Process an output string and return corrected final outcomes.

        Args:
            string (str): Output string of the code.
        Returns:
            corrected_z_logicals (list): A list of integers that are 0 or 1.
        These are the corrected values of the final transversal
        measurement, corresponding to the logical operators of
        self.z_logicals.
        """
        self.graph = deepcopy(self.decoding_graph.graph)
        string = "".join([str(c) for c in string[::-1]])
        output = [int(bit) for bit in list(string.split(" ")[0])][::-1]
        highlighted_nodes = self.code_circuit.string2nodes(
            string, logical=self.logical)
        if not highlighted_nodes:
            return output  # There's nothing for us to do here
        clusters = self.cluster(highlighted_nodes)

        for cluster in clusters:
            erasure = self.graph.subgraph(cluster)
            if isinstance(self.code_circuit, ArcCircuit):
                # NOTE: it just corrects for final logical readout
                for node in erasure.nodes():
                    if node["is_boundary"]:
                        qubit_to_be_corrected = node["qubits"][0]
                        output[qubit_to_be_corrected] = (
                            output[qubit_to_be_corrected]+1) % 2
                continue

            flipped_qubits = self.peeling(erasure)
            for qubit_to_be_corrected in flipped_qubits:
                output[qubit_to_be_corrected] = (
                    output[qubit_to_be_corrected]+1) % 2

        return output

    def cluster(self, nodes):
        """
        Create clusters using the union-find algorithm.

        Args:
            nodes (List): List of non-typical nodes in the syndrome graph, 
            of the type produced by `string2nodes`.

        Returns:
            FIXME: Make this more expressive.
            clusters (List[List[int]]): List of Lists of indices of nodes in clusters
        """
        node_indices = [self.graph.nodes().index(node) for node in nodes]
        for node_index, _ in enumerate(self.graph.nodes()):
            self.graph[node_index]["syndrome"] = node_index in node_indices
            self.graph[node_index]["root"] = node_index

        for edge in self.graph.edges():
            edge["growth"] = 0
            edge["fully_grown"] = False

        self.clusters: Dict[int, UnionFindDecoderCluster] = {}
        self.odd_cluster_roots = set(node_indices)
        for node_index in self.graph.node_indices():
            boundary_edges = []
            for edge_index, (_, neighbour, data) in dict(self.graph.incident_edge_index_map(node_index)).items():
                boundary_edges.append(
                    BoundaryEdge(
                        edge_index,
                        node_index,
                        neighbour,
                        data
                    )
                )
            self.clusters[node_index] = UnionFindDecoderCluster(
                boundary=boundary_edges,
                fully_grown_edges=set(),
                atypical_nodes=set(
                    [node_index]) if node_index in self.odd_cluster_roots else set([]),
                size=1
            )

        while self.odd_cluster_roots:
            fusion_edge_list = self._grow_clusters()
            self._merge_clusters(fusion_edge_list)

        cluster_nodes = []
        for _, cluster in self.clusters.items():
            if not cluster.atypical_nodes:
                continue
            nodes = set()
            for edge in cluster.fully_grown_edges:
                nodes |= set(self.graph.get_edge_endpoints_by_index(edge))
            cluster_nodes.append(list(nodes))
        return cluster_nodes

    def find(self, u: int) -> int:
        """
        Find() function as described in the paper that returns the root 
        of the cluster of a node, including path compression.

        Args:
            u (int): The index of the node in the decoding graph.

        Returns:
            root (int): The root of the cluster of node u.
        """
        if self.graph[u]["root"] == u:
            return self.graph[u]["root"]

        self.graph[u]["root"] = self.find(self.graph[u]["root"])
        return self.graph[u]["root"]

    def _grow_clusters(self) -> List[FusionEntry]:
        """
        Grow every "odd" cluster by half an edge.

        Returns:
            fusion_edge_list (List[FusionEntry]): List of edges that connect two 
            clusters that will be merged in the next step.
        """
        fusion_edge_list: List[FusionEntry] = []
        for root in self.odd_cluster_roots:
            cluster = self.clusters[root]
            for edge in cluster.boundary:
                edge.data["growth"] += 0.5
                if edge.data["growth"] >= edge.data["weight"] and not edge.data["fully_grown"]:
                    edge.data["fully_grown"] = True
                    cluster.fully_grown_edges.add(edge.index)
                    fusion_entry = FusionEntry(
                        u=edge.cluster_vertex, v=edge.neighbour_vertex, connecting_edge=edge)
                    fusion_edge_list.append(fusion_entry)
        return fusion_edge_list

    def _merge_clusters(self, fusion_edge_list: List[FusionEntry]) -> None:
        """
        Merges the clusters based on the fusion_edge_list computed in _grow_clusters().
        Updates the odd_clusters list by recomputing the neutrality of the newly merged clusters.

        Args:
            fusion_edge_list (List[FusionEntry]): List of edges that connect two 
            clusters that was computed in _grow_clusters().
        """
        for entry in fusion_edge_list:
            root_u, root_v = self.find(entry.u), self.find(entry.v)
            if root_u == root_v:
                continue
            new_root = root_v if self.clusters[root_v].size > self.clusters[root_u].size else root_u
            root_to_update = root_v if new_root == root_u else root_u

            cluster = self.clusters[new_root]
            other_cluster = self.clusters.pop(root_to_update)

            # Merge boundaries
            cluster.boundary += other_cluster.boundary
            cluster.boundary.remove(entry.connecting_edge)
            cluster.boundary.remove(entry.connecting_edge.reverse())

            cluster.atypical_nodes |= other_cluster.atypical_nodes
            cluster.fully_grown_edges |= other_cluster.fully_grown_edges
            cluster.size += other_cluster.size

            # update odd_cluster_roots
            if not self.code_circuit.is_cluster_even([self.graph[node] for node in cluster.atypical_nodes]):
                self.odd_cluster_roots.add(new_root)
            else:
                self.odd_cluster_roots.discard(new_root)
            self.odd_cluster_roots.discard(root_to_update)
            self.graph[root_to_update]["root"] = new_root

    def peeling(self, erasure: PyGraph) -> List[int]:
        """"
        Runs the peeling decoder on the erasure provided.
        Assumes that the erasure is one connected component, if not it will run in an 
        infinite loop in the tree construction.
        It works by first producing a spanning forest of the erasure and then 
        going backwards through the edges of the tree computing the error based on the syndrome.
        Based on arXiv:1703.01517.

        TODO: Extract to a separate decoder.

        Args:
            erasure (PyGraph): subgraph of the syndrome graph that represents the erasure.

        Returns:
            errors (List[int]): List of qubit indices on which Pauli errors occurred. 
        """
        tree = SpanningForest(vertices={}, edges=[])

        # Construct spanning forest
        # Pick starting vertex
        for vertex in erasure.node_indices():
            if erasure[vertex]["is_boundary"]:
                tree.vertices[vertex] = []
                break
        if not tree.vertices:
            tree.vertices[erasure.node_indices()[0]] = []

        # Expand forest |V| - 1 times, constructing it
        while len(tree.edges) < len(erasure.nodes()) - 1:
            vertices = copy(tree.vertices)
            for node in vertices.keys():
                for edge in erasure.incident_edges(node):
                    neighbour = list(
                        set(erasure.get_edge_endpoints_by_index(edge)) - set([node]))[0]
                    if not neighbour in tree.vertices.keys() and not erasure[neighbour]["is_boundary"]:
                        tree.edges.append(edge)
                        tree.vertices[neighbour] = []
                        tree.vertices[node].append(edge)
                        break

        edges = set()
        for edge in tree.edges[::-1]:
            endpoints = erasure.get_edge_endpoints_by_index(edge)
            pendant_vertex = endpoints[0] if not tree.vertices[endpoints[0]
                                                               ] else endpoints[1]
            if erasure[pendant_vertex]["is_boundary"]:
                pendant_vertex = endpoints[0] if pendant_vertex == endpoints[1] else endpoints[1]
            tree_vertex = endpoints[0] if pendant_vertex == endpoints[1] else endpoints[1]
            tree.vertices[tree_vertex].remove(edge)
            if erasure[pendant_vertex]["syndrome"] and not erasure[pendant_vertex]["is_boundary"]:
                edges.add(edge)
                erasure[tree_vertex]["syndrome"] = not erasure[tree_vertex]["syndrome"]
                erasure[pendant_vertex]["syndrome"] = False

        return [erasure.edges()[edge]["qubits"][0] for edge in edges if erasure.edges()[edge]["qubits"]]
