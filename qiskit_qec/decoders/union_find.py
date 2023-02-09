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
from typing import Dict, List, Set, Tuple
from rustworkx import PyGraph, visualization

from qiskit_qec.decoders.decoding_graph import DecodingGraph
from qiskit_qec.circuits import SurfaceCodeCircuit

from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class UnionFindDecoder:
    """
    Decoder based on growing clusters around syndrome errors to 
    "convert" them into erasure errors, which can be corrected easily, by using a weighted growth version of 
    the original algorithm to achieve almost-linear time decoding.

    See arXiv:1709.06218v3 for more details.
    """

    def __init__(
        self,
        code_circuit: SurfaceCodeCircuit
    ) -> None:
        self.code_circuit = code_circuit
        self.decoding_graph = DecodingGraph(
            code_circuit
        )

    def process(self, string: str):
        self.graph = deepcopy(self.decoding_graph.graph)
        string = "".join([str(c) for c in string[::-1]])
        output = [int(bit) for bit in list(string.split(" ")[0])][::-1]
        highlighted_nodes = self.code_circuit.string2nodes(string)
        if not highlighted_nodes: return output # There's nothing for us to do here
        highlighted_nodes_indices = [self.graph.nodes().index(
            node) for node in highlighted_nodes]
        for index, _ in enumerate(self.graph.nodes()):
            self.graph[index]["syndrome"] = index in highlighted_nodes_indices

        for edge in self.graph.edges():
            edge["growth"] = 0

        self.odd_clusters: Set[UnionFindDecoderCluster] = set(
            [UnionFindDecoderCluster(self.graph, self.graph[index])
             for index in highlighted_nodes_indices]
        )
        self.clusters: Set[UnionFindDecoderCluster] = set()
        while self.odd_clusters:
            self._grow_clusters()
            self._merge_clusters()
            self._update_clusters()
        
        erasure_vertices = set()
        for i, edge in enumerate(self.graph.edges()):
            if edge["weight"] <= edge["growth"]:
                endpoints = self.graph.get_edge_endpoints_by_index(i)
                for endpoint in endpoints: erasure_vertices.add(endpoint)
        
        erasure = self.graph.subgraph(list(erasure_vertices))
        
        flipped_qubits = self.peeling(erasure)

        for qubit_to_be_corrected in flipped_qubits:
            output[qubit_to_be_corrected] = (output[qubit_to_be_corrected]+1) % 2

        return output

    def _grow_clusters(self) -> None:
        for cluster in self.odd_clusters:
            outer_edges = copy(cluster.outer_edges)
            for index in outer_edges:
                edge = self.graph.edges()[index]
                edge["growth"] += 0.5
                self.graph.update_edge_by_index(index, edge)
                if edge["growth"] >= edge["weight"]:
                    nodes = list(self.graph.get_edge_endpoints_by_index(index))
                    cluster.add_nodes(nodes)
                    cluster.remove_outer_edge(index)

    def _merge_clusters(self) -> None:
        odd_clusters = copy(self.odd_clusters)
        for cluster_x in odd_clusters:
            for cluster_y in odd_clusters:
                if cluster_x == cluster_y:
                    continue
                if not cluster_x in self.odd_clusters or not cluster_y in self.odd_clusters:
                    continue
                if not (cluster_x.nodes & cluster_y.nodes):
                    continue  # They have no overlap, ignore them

                cluster_x.merge(cluster_y)
                self.odd_clusters.remove(cluster_y)

    def _update_clusters(self) -> None:
        odd_clusters = copy(self.odd_clusters)
        for cluster in odd_clusters:
            if not cluster.is_odd():
                self.clusters.add(cluster)
                self.odd_clusters.remove(cluster)

    def peeling(self, erasure: PyGraph) -> List[int]:
        """"
        Peeling decoder based on arXiv:1703.01517.
        """

        tree = SpanningForest(vertices={}, edges=[])

        ## Construct spanning forest
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
                    neighbour = list(set(erasure.get_edge_endpoints_by_index(edge)) - set([node]))[0]
                    if not neighbour in tree.vertices.keys():
                        tree.edges.append(edge)
                        tree.vertices[neighbour] = []
                        tree.vertices[node].append(edge)
                        break
        
        edges = set()
        for edge in tree.edges[::-1]:
            endpoints = erasure.get_edge_endpoints_by_index(edge)
            pendant_vertex = endpoints[0] if not tree.vertices[endpoints[0]] else endpoints[1]
            tree_vertex = endpoints[0] if pendant_vertex == endpoints[1] else endpoints[1]
            tree.vertices[tree_vertex].remove(edge)
            if erasure[pendant_vertex]["syndrome"] and not erasure[pendant_vertex]["is_boundary"]:
                edges.add(edge)
                erasure[tree_vertex]["syndrome"] = not erasure[tree_vertex]["syndrome"]
                erasure[pendant_vertex]["syndrome"] = False
            
        return [erasure.edges()[edge]["qubits"][0] for edge in edges if erasure.edges()[edge]["qubits"]]


@dataclass
class SpanningForest:
    vertices: Dict[int, List[int]]
    edges: List[int]


class UnionFindDecoderCluster:
    def __init__(self, graph: PyGraph, initial_node=None) -> None:
        self.graph = graph
        initial_node_index = self.graph.nodes().index(initial_node)
        self.nodes = set([initial_node_index])
        self.outer_edges = set(
            [edge_index for edge_index in self.graph.incident_edges(initial_node_index)])
    
    def add_nodes(self, nodes_indices: List[int]) -> None:
        for node_index in nodes_indices:
            self.nodes.add(node_index)
            for edge_index in self.graph.incident_edges(node_index):
                self.outer_edges.add(edge_index)

    def remove_outer_edge(self, edge_index):
        self.outer_edges.remove(edge_index)

    def is_odd(self):
        atypical_node_count = 0
        for index in self.nodes:
            if self.graph[index]["syndrome"]:
                atypical_node_count += 1
        return atypical_node_count % 2 == 1
    
    def merge(self, rhs):
        self.nodes |= rhs.nodes
        self.outer_edges |= rhs.outer_edges
