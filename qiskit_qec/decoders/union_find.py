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
from typing import Dict, List
from rustworkx import PyGraph

from qiskit_qec.decoders.decoding_graph import DecodingGraph

@dataclass
class SpanningForest:
    vertices: Dict[int, List[int]]
    edges: List[int]

@dataclass
class BoundaryEdge:
    cluster_vertex: int
    neighbour_vertex: int
    data: Dict[str, object]

@dataclass
class UnionFindDecoderCluster:
    boundary: List[BoundaryEdge]
    is_odd: bool
    size: int 

@dataclass
class FusionEntry:
    u: int
    v: int
    connecting_edge: BoundaryEdge

class UnionFindDecoder:
    """
    Decoder based on growing clusters around syndrome errors to 
    "convert" them into erasure errors, which can be corrected easily, by using a weighted growth version of 
    the original algorithm to achieve almost-linear time decoding.

    See arXiv:1709.06218v3 for more details.
    """

    def __init__(
        self,
        code_circuit,
        logical: str
    ) -> None:
        self.code_circuit = code_circuit
        self.decoding_graph = DecodingGraph(
            code_circuit,
        )
        self.logical = logical

    def process(self, string: str):
        self.graph = deepcopy(self.decoding_graph.graph)
        string = "".join([str(c) for c in string[::-1]])
        output = [int(bit) for bit in list(string.split(" ")[0])][::-1]
        highlighted_nodes = self.code_circuit.string2nodes(string, logical=self.logical)
        if not highlighted_nodes: return output # There's nothing for us to do here
        highlighted_nodes_indices = [self.graph.nodes().index(
            node) for node in highlighted_nodes]
        for index, _ in enumerate(self.graph.nodes()):
            self.graph[index]["syndrome"] = index in highlighted_nodes_indices
            self.graph[index]["root"] = index

        for edge in self.graph.edges():
            edge["growth"] = 0
            edge["fully_grown"] = False

        self.clusters: Dict[int, UnionFindDecoderCluster] = {}
        self.odd_cluster_roots = highlighted_nodes_indices
        for index in self.graph.node_indices():
            boundary_edges = []
            for _, (_, neighbour, data) in dict(self.graph.incident_edge_index_map(index)).items():
                boundary_edges.append(
                    BoundaryEdge(
                        index,
                        neighbour,
                        data
                    )
                )
            self.clusters[index] = UnionFindDecoderCluster(
                boundary=boundary_edges,
                is_odd=index in highlighted_nodes,
                size=1
            )
        
        while self.odd_cluster_roots:
            fusion_edge_list = self._grow_clusters()
            self._merge_clusters(fusion_edge_list)
        
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

    def find(self, u: int) -> int:
        if self.graph[u]["root"] == u:
            return self.graph[u]["root"]
        
        self.graph[u]["root"] = self.find(self.graph[u]["root"])
        return self.graph[u]["root"]
    
    def _grow_clusters(self) -> List[FusionEntry]:
        fusion_edge_list: List[FusionEntry] = []
        for root in self.odd_cluster_roots:
            cluster = self.clusters[root]
            for edge in cluster.boundary:
                edge.data["growth"] += 0.5
                if edge.data["growth"] >= edge.data["weight"] and not edge.data["fully_grown"]:
                    edge.data["fully_grown"] = True
                    fusion_entry = FusionEntry(u=edge.cluster_vertex, v=edge.neighbour_vertex, connecting_edge=edge)
                    fusion_edge_list.append(fusion_entry)
        return fusion_edge_list


    def _merge_clusters(self, fusion_edge_list: List[FusionEntry]) -> None:
        for entry in fusion_edge_list:
            root_u, root_v = self.find(entry.u), self.find(entry.v)
            if root_u == root_v:
                continue
            new_root = root_v if self.clusters[root_v].size > self.clusters[root_u].size else root_u
            root_to_update = root_v if new_root == root_u else root_u

            cluster = self.clusters[new_root]
            # Merge boundaries
            cluster.boundary += self.clusters[root_to_update].boundary
            reverse_edge = BoundaryEdge(
                cluster_vertex=entry.connecting_edge.neighbour_vertex,
                neighbour_vertex=entry.connecting_edge.cluster_vertex,
                data=entry.connecting_edge.data
            )
            cluster.boundary.remove(entry.connecting_edge)
            cluster.boundary.remove(reverse_edge)
            
            # update size
            cluster.size += self.clusters[root_to_update].size
            # update parity
            cluster.is_odd ^= self.clusters[root_to_update].is_odd
            # update root
            self.graph[root_to_update]["root"] = new_root
            
            if root_to_update in self.odd_cluster_roots:
                self.odd_cluster_roots.remove(root_to_update)

            # update odd_cluster_roots
            if not cluster.is_odd and new_root in self.odd_cluster_roots:
                self.odd_cluster_roots.remove(new_root)    
            if cluster.is_odd and not new_root in self.odd_cluster_roots:
                self.odd_cluster_roots.append(new_root)

    def _update_clusters(self) -> None:
        odd_clusters = copy(self.odd_clusters)
        for root, cluster in odd_clusters.items():
            if not cluster.is_odd:
                self.clusters[root] = self.odd_clusters.pop(root, None)

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
