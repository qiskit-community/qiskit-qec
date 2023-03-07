# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2022.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# pylint: disable=invalid-name

"""Hard decision renormalization group decoders."""

from copy import copy, deepcopy
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple, Tuple
from rustworkx import connected_components, distance_matrix, PyGraph

from qiskit_qec.circuits.repetition_code import ArcCircuit, RepetitionCodeCircuit
from qiskit_qec.analysis.decoding_graph import DecodingGraph, Node, Edge
from qiskit_qec.exceptions import QiskitQECError


class ClusteringDecoder:
    """
    Generic base class for clustering decoders.
    """

    def __init__(
        self,
        code_circuit,
        decoding_graph: DecodingGraph = None,
    ):
        self.code = code_circuit
        if decoding_graph:
            self.decoding_graph = decoding_graph
        else:
            self.decoding_graph = DecodingGraph(self.code)


class BravyiHaahDecoder(ClusteringDecoder):
    """Decoder based on finding connected components within the decoding graph."""

    def __init__(
        self,
        code_circuit,
        decoding_graph: DecodingGraph = None,
    ):
        if not isinstance(code_circuit, (ArcCircuit, RepetitionCodeCircuit)):
            raise QiskitQECError("Error: code_circuit not supported.")

        super().__init__(code_circuit, decoding_graph)

        if isinstance(self.code, ArcCircuit):
            self.z_logicals = self.code.z_logicals
        elif isinstance(self.code, RepetitionCodeCircuit):
            if self.code._xbasis:
                self.z_logicals = self.code.css_x_logical[0]
            else:
                self.z_logicals = self.code.css_z_logical[0]
        if isinstance(self.code, ArcCircuit):
            self.code_index = self.code.code_index
        elif isinstance(self.code, RepetitionCodeCircuit):
            self.code_index = {2 * j: j for j in range(self.code.d)}

    def _cluster(self, ns, dist_max):
        """
        Finds connected components in the given nodes, for nodes connected by at most the given distance
        in the given decoding graph.
        """

        # calculate distance for the graph
        dg = self.decoding_graph.graph
        distance = distance_matrix(dg)

        # create empty `DecodingGraph`
        cluster_graph = DecodingGraph(None)
        cg = cluster_graph.graph
        # add all the given nodes to cg
        d2c = {}
        c2g = {}
        for n in ns:
            node = dg.nodes()[n]
            d2c[n] = cg.add_node(node)
            c2g[d2c[n]] = n
        # add an edge between a pair of the given nodes if their distance is small enough
        for n0 in ns:
            for n1 in ns:
                if n0 < n1:
                    dist = distance[n0, n1]
                    if dist <= dist_max:
                        cg.add_edge(d2c[n0], d2c[n1], {"distance": dist})
        # find the connected components of cg
        con_comps = connected_components(cg)

        # use these to define clusters
        clusters = {}
        con_comp_dict = {}
        for c, con_comp in enumerate(con_comps):
            con_comp_dict[c] = []

            # check the neutrality of each connected component
            con_nodes = [cg[n] for n in con_comp]
            neutral, logicals, num_errors = self.code.check_nodes(
                con_nodes, ignore_extra_boundary=True
            )

            # it's fully neutral if no extra logicals are needed
            # and if the error num is less than the max dist
            fully_neutral = neutral and logicals == [] and num_errors < dist_max

            # if a cluster is neutral, all nodes are labelled with c
            # otherwise, it gets a None
            for n in con_comp:
                if fully_neutral:
                    clusters[c2g[n]] = c
                else:
                    clusters[c2g[n]] = None
                con_comp_dict[c].append(c2g[n])

        return clusters, con_comp_dict

    def _get_boundary_nodes(self):
        boundary_nodes = []
        for element, z_logical in enumerate(self.z_logicals):
            node = Node(
                is_boundary=True,
                qubits=[z_logical],
                index=element
            )
            if isinstance(self.code, ArcCircuit):
                node.properties["link qubit"] = None
            boundary_nodes.append(node)
        return boundary_nodes

    def cluster(self, nodes):
        """

        Args:
            nodes (list): List of nodes, of the type produced by `string2nodes`.
        Returns:
            final_clusters (dict): Dictionary with the indices of the given node
            as keys and an integer specifying their cluster as the corresponding
            value.
        """

        # get indices for nodes and boundary nodes
        dg = self.decoding_graph.graph
        ns = set(dg.nodes().index(node) for node in nodes)
        bns = set(dg.nodes().index(node) for node in self._get_boundary_nodes())

        dist_max = 0
        final_clusters = {}
        con_comps = []
        clusterss = []
        while ns and dist_max <= self.code.d:
            dist_max += 1
            # add boundary nodes to unpaired nodes
            ns = set(ns).union(bns)

            # cluster nodes and contract decoding graph given the current distance
            clusters, con_comp = self._cluster(ns, dist_max)
            # record the clustered and unclustered nodes
            ns = []
            for n, c in clusters.items():
                if c is not None:
                    final_clusters[n] = c
                else:
                    if not dg[n].is_boundary:
                        ns.append(n)
            con_comps.append(con_comp)
            clusterss.append(clusters)

        return final_clusters

    def process(self, string):
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
        code = self.code
        decoding_graph = self.decoding_graph

        # turn string into nodes and cluster
        nodes = code.string2nodes(string)
        clusters = self.cluster(nodes)

        # get the list of bulk nodes for each cluster
        cluster_nodes = {c: [] for c in clusters.values()}
        for n, c in clusters.items():
            node = decoding_graph.graph[n]
            if not node.is_boundary:
                cluster_nodes[c].append(node)

        # get the list of required logicals for each cluster
        cluster_logicals = {}
        for c, nodes in cluster_nodes.items():
            _, logical_nodes, _ = code.check_nodes(nodes)
            z_logicals = [node.qubits[0] for node in logical_nodes]
            cluster_logicals[c] = z_logicals

        # get the net effect on each logical
        net_z_logicals = {z_logical: 0 for z_logical in self.z_logicals}
        for c, z_logicals in cluster_logicals.items():
            for z_logical in self.z_logicals:
                if z_logical in z_logicals:
                    net_z_logicals[z_logical] += 1
        for z_logical, num in net_z_logicals.items():
            net_z_logicals[z_logical] = num % 2

        corrected_z_logicals = []
        string = string.split(" ")[0]
        for z_logical in self.z_logicals:
            raw_logical = int(string[-1 - self.code_index[z_logical]])
            corrected_logical = (raw_logical + net_z_logicals[z_logical]) % 2
            corrected_z_logicals.append(corrected_logical)

        return corrected_z_logicals


@dataclass
class SpanningForest:
    """
    Spanning forest for the peeling decoder.
    """

    vertices: Dict[int, List[int]]
    edges: List[int]


@dataclass
class BoundaryEdge:
    """
    Boundary edge for the boundary of a UnionFindDecoderCluster.
    """

    index: int
    cluster_vertex: int
    neighbour_vertex: int
    data: Edge

    def reverse(self):
        """
        Returns a reversed version of the boundary edge (cluster and neighbour vertex flipped)
        """
        return BoundaryEdge(
            index=self.index,
            cluster_vertex=self.neighbour_vertex,
            neighbour_vertex=self.cluster_vertex,
            data=self.data,
        )


@dataclass
class UnionFindDecoderCluster:
    """
    Cluster for the UnionFindDecoder
    """

    boundary: List[BoundaryEdge]
    atypical_nodes: Set[int]
    fully_grown_edges: Set[int]
    size: int


@dataclass
class FusionEntry:
    """
    Entry for the fusion list between the growing and merging of the union find decoder.
    """

    u: int
    v: int
    connecting_edge: BoundaryEdge


class UnionFindDecoder(ClusteringDecoder):
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
        code,
        logical: str,
        decoding_graph: DecodingGraph = None,
    ) -> None:
        super().__init__(code, decoding_graph)
        self.logical = logical
        self.graph = deepcopy(self.decoding_graph.graph)
        self.clusters: List[List[int]] = []
        self.odd_cluster_roots: Set[int] = []

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
        output = [int(bit) for bit in list(string.split(" ", maxsplit=self.code.d)[0])][::-1]
        highlighted_nodes = self.code.string2nodes(string, logical=self.logical)
        if not highlighted_nodes:
            return output  # There's nothing for us to do here
        clusters = self.cluster(highlighted_nodes)

        for cluster in clusters:
            erasure = self.graph.subgraph(cluster)
            if isinstance(self.code, ArcCircuit):
                # NOTE: it just corrects for final logical readout
                for node in erasure.nodes():
                    if node.is_boundary:
                        # FIXME: Find a general way to go from physical qubit
                        # index to code qubit index
                        qubit_to_be_corrected = int(node.qubits[0] / 2)
                        output[qubit_to_be_corrected] = (output[qubit_to_be_corrected] + 1) % 2
                continue

            flipped_qubits = self.peeling(erasure)
            for qubit_to_be_corrected in flipped_qubits:
                output[qubit_to_be_corrected] = (output[qubit_to_be_corrected] + 1) % 2

        return output

    def cluster(self, nodes) -> List[List[int]]:
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
            self.graph[node_index].properties["syndrome"] = node_index in node_indices
            self.graph[node_index].properties["root"] = node_index

        for edge in self.graph.edges():
            edge.properties["growth"] = 0
            edge.properties["fully_grown"] = False

        self.clusters: Dict[int, UnionFindDecoderCluster] = {}
        self.odd_cluster_roots = set(node_indices)
        for node_index in self.graph.node_indices():
            boundary_edges = []
            for edge_index, (_, neighbour, data) in dict(
                self.graph.incident_edge_index_map(node_index)
            ).items():
                boundary_edges.append(BoundaryEdge(edge_index, node_index, neighbour, data))
            self.clusters[node_index] = UnionFindDecoderCluster(
                boundary=boundary_edges,
                fully_grown_edges=set(),
                atypical_nodes=set([node_index])
                if node_index in self.odd_cluster_roots
                else set([]),
                size=1,
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
        if self.graph[u].properties["root"] == u:
            return self.graph[u].properties["root"]

        self.graph[u].properties["root"] = self.find(self.graph[u].properties["root"])
        return self.graph[u].properties["root"]

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
                edge.data.properties["growth"] += 0.5
                if edge.data.properties["growth"] >= edge.data.weight and not edge.data.properties["fully_grown"]:
                    edge.data.properties["fully_grown"] = True
                    cluster.fully_grown_edges.add(edge.index)
                    fusion_entry = FusionEntry(
                        u=edge.cluster_vertex, v=edge.neighbour_vertex, connecting_edge=edge
                    )
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
            if not self.code.is_cluster_neutral(
                [self.graph[node] for node in cluster.atypical_nodes]
            ):
                self.odd_cluster_roots.add(new_root)
            else:
                self.odd_cluster_roots.discard(new_root)
            self.odd_cluster_roots.discard(root_to_update)
            self.graph[root_to_update].properties["root"] = new_root

    def peeling(self, erasure: PyGraph) -> List[int]:
        """ "
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
            if erasure[vertex].is_boundary:
                tree.vertices[vertex] = []
                break
        if not tree.vertices:
            tree.vertices[erasure.node_indices()[0]] = []

        # Expand forest |V| - 1 times, constructing it
        while len(tree.edges) < len(erasure.nodes()) - 1:
            vertices = copy(tree.vertices)
            for node in vertices.keys():
                if len(tree.edges) >= len(erasure.nodes()) - 1:
                    break
                for edge, (_, neighbour, _) in dict(erasure.incident_edge_index_map(node)).items():
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
            if erasure[pendant_vertex].properties["syndrome"] and not erasure[pendant_vertex].is_boundary:
                edges.add(edge)
                erasure[tree_vertex].properties["syndrome"] = not erasure[tree_vertex].properties["syndrome"]
                erasure[pendant_vertex].properties["syndrome"] = False

        return [
            erasure.edges()[edge].qubits[0] for edge in edges if erasure.edges()[edge].qubits
        ]

@dataclass
class Cluster:
    """
    Cluster class for the ClAYG decoder. 
    FIXME: Remove, when ClAYG decoder uses Union Find infrastructure.
    """
    boundary: List[Tuple[int, int]] # List[(edge, neighbour)]
    fully_grown_edges: List[int] # List[edge]
    nodes: List[int] # List[node_index]
    atypical_nodes: List[int] # List[node_index]

class ClAYGDecoder(UnionFindDecoder):
    """
    Decoder that is very similar to the Union Find decoder, but instead of adding clusters all at once, 
    adds them separated by syndrome round with a growth and merge phase in between.
    Then it just proceeds like the Union Find decoder.

    FIXME: Use the Union Find infrastructure and just change the self.cluster() method. Problem is that
    the peeling decoder needs a modified version the graph with the syndrome nodes marked, which is done 
    in the process method. For now it is mostly its separate thing, but merging them shouldn't be 
    too big of a hassle.
    Merge method should also be modified, as boundary clusters are not marked as odd clusters.
    """
    def __init__(self, code, logical: str, decoding_graph: DecodingGraph = None) -> None:
        super().__init__(code, logical, decoding_graph)
        self.graph = deepcopy(self.decoding_graph.graph)
    
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
        nodes_at_time_zero = []
        for index, node in enumerate(self.graph.nodes()):
            if node.time == 0 or node.is_boundary:
                nodes_at_time_zero.append(index)
        self.graph = self.graph.subgraph(nodes_at_time_zero)

        for edge in self.graph.edges():
            edge.weight = 1
            edge.properties["growth"] = 0
        
        string = "".join([str(c) for c in string[::-1]])
        output = [int(bit) for bit in list(string.split(" ", maxsplit=self.code.d)[0])][::-1]
        nodes = self.code.string2nodes(string, logical=self.logical)

        clusters = self.cluster(nodes)

        for cluster in clusters:
            erasure_graph = deepcopy(self.graph)
            for node in cluster.nodes:
                erasure_graph[node].properties["syndrome"] = False
            for node in cluster.atypical_nodes:
                erasure_graph[node].properties["syndrome"] = True
            erasure = erasure_graph.subgraph(cluster.nodes + cluster.atypical_nodes)
            qubits_to_be_corrected = self.peeling(erasure)
            for idx in qubits_to_be_corrected:
                output[idx] = (output[idx] + 1) % 2
        
        return output
        

    def cluster(self, nodes) -> List[List[int]]:
        """
        Create clusters using the union-find algorithm.

        Args:
            nodes (List): List of non-typical nodes in the syndrome graph,
            of the type produced by `string2nodes`.

        Returns:
            FIXME: Make this more expressive. Maybe return a list of separate PyGraphs? Would fix the infrastructure-sharing-issue mentioned above.
            clusters (List[List[int]]): List of Lists of indices of nodes in clusters.
        """
        self.roots = {}
        self.odd = {}
        for i, node in enumerate(self.graph.nodes()):
            self.roots[i] = i
            self.odd[i] = False
        
        self.clusters: Dict[int, Cluster] = dict([(node_index, None) for node_index in self.graph.node_indices()])
        self.odd_cluster_roots: List[int] = []
        times = [[] for _ in range(self.code.T+1)]
        boundaries = []
        for node in deepcopy(nodes):
            if nodes.count(node) > 1: continue
            if node.is_boundary:
                boundaries.append(node)
            else:
                node.time = 0
                times[node.time].append(node)
        
        neutral_clusters = []

        for time in times:
            if not time: continue
            for node in time:
                neutral_clusters += self.add_atypical_node_to_decoding_graph(node, True)
            for root in self.odd_cluster_roots:
                neutral_clusters += self.grow_cluster_and_merge(root)
        
        for node in boundaries:
            neutral_clusters += self.add_atypical_node_to_decoding_graph(node, False)

        while self.odd_cluster_roots:
            for root in self.odd_cluster_roots:
                neutral_clusters += self.grow_cluster_and_merge(root)
        
        return neutral_clusters

    def add_atypical_node_to_decoding_graph(self, node, add_odd_cluster: bool) -> List[Cluster]:
        """
        Adds non-typical syndrome nodes to the graph and neutralize/create clusters around them if necessary

        Args:
            node: dictionary with node data in the form produced by string2nodes.
            add_odd_cluster (bool): specifices whether the newly created cluster is going to be added to the 
            odd_clusters_list.
        """
        node_index = self.graph.nodes().index(node)
        current_cluster_root = self.find(node_index)
        cluster = self.clusters[current_cluster_root]
        current_cluster_odd = self.odd[current_cluster_root]
        neutral_clusters = []
        # If cluster that it's in is odd set it to even and add it to the error log
        if current_cluster_odd:
            self.odd[current_cluster_root] = False
            self.odd_cluster_roots.remove(current_cluster_root)
            cluster.atypical_nodes.append(node_index)
            if not node_index == current_cluster_root:
                # Simple measurement error, don't add it to the error log
                # FIXME: Make peeling decoder prestage handle this
                neutral_clusters.append(cluster)
            for edge, _ in cluster.boundary:
                self.graph.edges()[edge].properties["growth"] = 0
            for node in cluster.nodes + cluster.atypical_nodes:
                self.roots[node] = node
            self.clusters[current_cluster_root] = None
        # Else create a new cluster around it and set it to odd
        else: 
            self.roots[node_index] = node_index
            self.odd[node_index] = True
            if add_odd_cluster:
                self.odd_cluster_roots.append(node_index)
            boundary: List[Tuple[int, int]] = []
            for edge, (_, neighbour, _) in dict(self.graph.incident_edge_index_map(node_index)).items():
                boundary.append((edge, neighbour))
            self.clusters[node_index] = Cluster(
                boundary=boundary,
                fully_grown_edges=[],
                nodes=[],
                atypical_nodes=[node_index]
            )

        return neutral_clusters
    
    def grow_cluster_and_merge(self, root: int):
        """
        Grows the cluster specified by root by half an edge and merges them if necessary.

        Args:
            root (int): index of the root node of the cluster.
        """
        cluster = self.clusters[root]
        if not cluster: return
        for edge, neighbour in copy(cluster.boundary):
            self.graph.edges()[edge].properties["growth"] += 0.5
            if self.graph.edges()[edge].properties["growth"] < self.graph.edges()[edge].weight:
                continue
            cluster.boundary.remove((edge, neighbour))
            cluster.fully_grown_edges.append(edge)
            self.graph.edges()[edge].properties["growth"] = 0
            neighbour_root = self.find(neighbour)
            if neighbour_root == root: continue
            neighbour_odd = self.odd[neighbour_root]
            if neighbour_odd:
                # It is odd, so there has to be a cluster
                neighbour_cluster = self.clusters[neighbour_root]
                cluster.boundary += neighbour_cluster.boundary
                cluster.fully_grown_edges += neighbour_cluster.fully_grown_edges
                cluster.nodes += neighbour_cluster.nodes
                cluster.atypical_nodes += neighbour_cluster.atypical_nodes
                for edge, _ in cluster.boundary:
                    self.graph.edges()[edge].properties["growth"] = 0
                for node in cluster.nodes + cluster.atypical_nodes:
                    self.roots[node] = node
                for root in [root, neighbour_root]:
                    if self.graph[root].is_boundary: continue
                    self.odd[root] = False
                    self.clusters[root] = None
                    self.odd_cluster_roots.remove(root)
                return [cluster]
            else:
                cluster.nodes += [neighbour]
                self.roots[neighbour] = root
                for edge, (_, neighbour_neighbour, _) in dict(self.graph.incident_edge_index_map(neighbour)).items():
                    if neighbour_neighbour == neighbour: continue
                    cluster.boundary.append((edge, neighbour_neighbour))
        return []

    def find(self, node_index):
        """
        Returns the root of the cluster the node belongs to.

        Args:
            node_index (int): index of the node in self.graph
        """
        if self.roots[node_index] == node_index:
            return node_index
        self.roots[node_index] = self.find(self.roots[node_index])
        return self.roots[node_index]
        