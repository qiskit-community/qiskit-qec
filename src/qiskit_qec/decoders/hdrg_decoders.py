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

from abc import ABC
from copy import copy, deepcopy
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple

import numpy as np
from rustworkx import PyGraph, connected_components, distance_matrix, adjacency_matrix
from rustworkx.visualization import graphviz_draw

from qiskit_qec.codes.bb_code import BBCode
from qiskit_qec.decoders.decoding_graph import DecodingGraph
from qiskit_qec.linear.matrix import solve2, LinAlgError
from qiskit_qec.utils import DecodingGraphEdge


class ClusteringDecoder(ABC):
    """
    Generic base class for clustering decoders.
    """

    def __init__(
        self,
        code_circuit,
        decoding_graph: DecodingGraph = None,
    ):
        self.code = code_circuit

        if hasattr(self.code, "code_index"):
            self.code_index = self.code.code_index
        else:
            self.code_index = {j: j for j in range(self.code.n)}

        if decoding_graph:
            self.decoding_graph = decoding_graph
        else:
            self.decoding_graph = DecodingGraph(self.code)

    def get_corrections(self, string, clusters):
        """
        Turn a set of neutral clusters into corrections.

        Args:
            string (str): Output string of the code
            clusters (dict): Dictionary with the indices of the given node
            as keys and an integer specifying their cluster as the corresponding
            value.
        Returns:
            corrected_logicals (list): A list of integers that are 0 or 1.
        These are the corrected values of the final transversal
        measurement, in the same form as given by the code's `string2raw_logicals`.
        """

        # get the list of bulk nodes for each cluster
        cluster_nodes = {c: [] for c in clusters.values()}
        for n, c in clusters.items():
            node = self.decoding_graph.graph[n]
            if not node.is_logical:
                cluster_nodes[c].append(node)

        # get the list of required logicals for each cluster
        cluster_logicals = {}
        for c, nodes in cluster_nodes.items():
            _, logical_nodes, _ = self.code.check_nodes(nodes, minimal=True)
            log_indexes = [node.index for node in logical_nodes]
            cluster_logicals[c] = log_indexes

        # get the net effect on each logical
        net_logicals = {node.index: 0 for node in self.decoding_graph.logical_nodes}
        for c, log_indexes in cluster_logicals.items():
            for log_index in log_indexes:
                net_logicals[log_index] += 1
        for log_index, num in net_logicals.items():
            net_logicals[log_index] = num % 2

        corrected_logicals = self.code.string2raw_logicals(string)
        for log_index, log_value in enumerate(corrected_logicals):
            corrected_logicals[log_index] = (net_logicals[log_index] + int(log_value)) % 2

        return corrected_logicals


class BravyiHaahDecoder(ClusteringDecoder):
    """Decoder based on finding connected components within the decoding graph."""

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
                con_nodes, ignore_extra_logical=True
            )

            # it's fully neutral if no extra logicals are needed
            # and if the error num is less than the max dist
            fully_neutral = neutral and logicals == []
            if num_errors:
                fully_neutral = fully_neutral and num_errors < dist_max

            # if a cluster is neutral, all nodes are labelled with c
            # otherwise, it gets a None
            for n in con_comp:
                if fully_neutral:
                    clusters[c2g[n]] = c
                else:
                    clusters[c2g[n]] = None
                con_comp_dict[c].append(c2g[n])

        return clusters, con_comp_dict

    def cluster(self, nodes):
        """

        Args:
            nodes (list): List of nodes, of the type produced by `string2nodes`.
        Returns:
            final_clusters (dict): Dictionary with the indices of the given node
            as keys and an integer specifying their cluster as the corresponding
            value.
        """

        # get indices for nodes and logical nodes
        dg = self.decoding_graph.graph
        ns = set(dg.nodes().index(node) for node in nodes)
        lns = set(dg.nodes().index(node) for node in self.decoding_graph.logical_nodes)

        dist_max = 0
        final_clusters = {}
        con_comps = []
        clusterss = []
        while ns and dist_max <= self.code.d:
            dist_max += 1
            # add logical nodes to unpaired nodes
            ns = set(ns).union(lns)

            # cluster nodes and contract decoding graph given the current distance
            clusters, con_comp = self._cluster(ns, dist_max)
            # record the clustered and unclustered nodes
            ns = []
            for n, c in clusters.items():
                if c is not None:
                    final_clusters[n] = c
                else:
                    if not dg[n].is_logical:
                        ns.append(n)
            con_comps.append(con_comp)
            clusterss.append(clusters)

        return final_clusters

    def process(self, string, predecoder=None):
        """
        Process an output string and return corrected final outcomes.

        Args:
            string (str): Output string of the code.
            predecoder (callable): Function that takes in and returns
            a list of nodes. Used to do preprocessing on the nodes
            corresponding to the input string.

        Returns:
            corrected_logicals (list): A list of integers that are 0 or 1.
        These are the corrected values of the final transversal
        measurement, in the same form as given by the code's `string2raw_logicals`.
        """

        # turn string into nodes and cluster
        nodes = self.code.string2nodes(string, all_logicals=True)
        # apply predecoder if one is given
        if predecoder:
            nodes = predecoder(nodes)
        # then cluster
        clusters = self.cluster(nodes)

        return self.get_corrections(string, clusters)


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
    data: DecodingGraphEdge

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
    boundary_nodes: Set[int]
    nodes: Set[int]
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
    by the peeling decoder for compatible codes or by the standard HDRG
    method in general.

    TODO: Add weights to edges of graph according to Huang et al (see. arXiv:2004.04693, section III)

    See arXiv:1709.06218v3 for more details.
    """

    def __init__(self, code, decoding_graph: DecodingGraph = None, use_peeling=True) -> None:
        super().__init__(code, decoding_graph=decoding_graph)
        self.graph = deepcopy(self.decoding_graph.graph)
        self.clusters: Dict[int, UnionFindDecoderCluster] = {}
        self.odd_cluster_roots: List[int] = []
        self.use_peeling = use_peeling
        self._clusters4peeling = []

    def process(self, string: str, predecoder=None):
        """
        Process an output string and return corrected final outcomes.

        Args:
            string (str): Output string of the code.
            predecoder (callable): Function that takes in and returns
            a list of nodes. Used to do preprocessing on the nodes
            corresponding to the input string.
        Returns:
            corrected_logicals (list): A list of integers that are 0 or 1.
        These are the corrected values of the final logical measurement.
        """

        if self.use_peeling:
            self.graph = deepcopy(self.decoding_graph.graph)
            highlighted_nodes = self.code.string2nodes(string, all_logicals=True)
            if predecoder:
                highlighted_nodes = predecoder(highlighted_nodes)

            # call cluster to do the clustering, but actually use the peeling form
            self.cluster(highlighted_nodes)
            clusters = self._clusters4peeling

            # determine the net logical z
            measured_logicals = {}
            for node in self.decoding_graph.logical_nodes:
                measured_logicals[node.index] = node.qubits
            net_z_logicals = {tuple(z_logical): 0 for z_logical in measured_logicals.values()}
            for cluster_nodes, _ in clusters:
                erasure = self.graph.subgraph(cluster_nodes)
                flipped_qubits = self.peeling(erasure)
                for qubit_to_be_corrected in flipped_qubits:
                    for z_logical in net_z_logicals:
                        if qubit_to_be_corrected in z_logical:
                            net_z_logicals[z_logical] += 1
            for z_logical, num in net_z_logicals.items():
                net_z_logicals[z_logical] = num % 2

            # apply this to the raw readout
            corrected_z_logicals = []
            raw_logicals = self.code.string2raw_logicals(string)
            for j, z_logical in measured_logicals.items():
                raw_logical = int(raw_logicals[j])
                corrected_logical = (raw_logical + net_z_logicals[tuple(z_logical)]) % 2
                corrected_z_logicals.append(corrected_logical)
            return corrected_z_logicals
        else:
            # turn string into nodes and cluster
            nodes = self.code.string2nodes(string, all_logicals=True)
            if predecoder:
                nodes = predecoder(nodes)
            clusters = self.cluster(nodes)
            return self.get_corrections(string, clusters)

    def cluster(self, nodes: List):
        """
        Create clusters using the union-find algorithm.

        Args:
            nodes (List): List of non-typical nodes in the syndrome graph,
            of the type produced by `string2nodes`.

        Returns:
            clusters (dict): Dictionary with the indices of
            the given node as keys and an integer specifying their cluster as the corresponding
            value.
        """
        node_indices = [self.decoding_graph.node_index(node) for node in nodes]
        for node_index, _ in enumerate(self.graph.nodes()):
            self.graph[node_index].properties["syndrome"] = node_index in node_indices
            self.graph[node_index].properties["root"] = node_index

        for edge in self.graph.edges():
            edge.properties["growth"] = 0
            edge.properties["fully_grown"] = False

        self.clusters: Dict[int, UnionFindDecoderCluster] = {}
        self.odd_cluster_roots = []
        for node_index in node_indices:
            self._create_new_cluster(node_index)

        j = 0
        while self.odd_cluster_roots and j < 2 * self.code.d * (self.code.T + 1):
            self._grow_and_merge_clusters()
            j += 1

        # compile info into standard clusters dict
        clusters = {}
        for c, cluster in self.clusters.items():
            # determine which nodes exactly are in the neutral cluster
            neutral_nodes = list(cluster.atypical_nodes | cluster.boundary_nodes)
            # put them in the required dict
            for n in neutral_nodes:
                clusters[n] = c

        # also compile into form required for peeling
        self._clusters4peeling = []
        for _, cluster in self.clusters.items():
            if not cluster.atypical_nodes:
                continue
            self._clusters4peeling.append(
                (list(cluster.nodes), list(cluster.atypical_nodes | cluster.boundary_nodes))
            )

        return clusters

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

    def _create_new_cluster(self, node_index):
        node = self.graph[node_index]
        if not node.is_logical:
            self.odd_cluster_roots.insert(0, node_index)
        boundary_edges = []
        for edge_index, neighbour, data in self.neighbouring_edges(node_index):
            boundary_edges.append(BoundaryEdge(edge_index, node_index, neighbour, data))
        self.clusters[node_index] = UnionFindDecoderCluster(
            boundary=boundary_edges,
            fully_grown_edges=set(),
            atypical_nodes=set([node_index]) if not node.is_logical else set([]),
            boundary_nodes=set([node_index]) if node.is_logical else set([]),
            nodes=set([node_index]),
            size=1,
        )

    def _grow_and_merge_clusters(self) -> Set[int]:
        fusion_edge_list = self._grow_clusters()
        return self._merge_clusters(fusion_edge_list)

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
                if (
                    edge.data.properties["growth"] >= edge.data.weight
                    and not edge.data.properties["fully_grown"]
                ):
                    neighbour_root = self.find(edge.neighbour_vertex)
                    if not neighbour_root in self.clusters:
                        boundary_edges = []
                        for edge_index, neighbour_neighbour, data in self.neighbouring_edges(
                            edge.neighbour_vertex
                        ):
                            boundary_edges.append(
                                BoundaryEdge(
                                    edge_index, edge.neighbour_vertex, neighbour_neighbour, data
                                )
                            )
                        self.graph[edge.neighbour_vertex].properties["root"] = edge.neighbour_vertex
                        self.clusters[edge.neighbour_vertex] = UnionFindDecoderCluster(
                            boundary=boundary_edges,
                            fully_grown_edges=set(),
                            atypical_nodes=set(),
                            boundary_nodes=set([edge.neighbour_vertex])
                            if self.graph[edge.neighbour_vertex].is_logical
                            else set([]),
                            nodes=set([edge.neighbour_vertex]),
                            size=1,
                        )
                    fusion_entry = FusionEntry(
                        u=edge.cluster_vertex, v=edge.neighbour_vertex, connecting_edge=edge
                    )
                    fusion_edge_list.append(fusion_entry)
        return fusion_edge_list

    def _merge_clusters(self, fusion_edge_list: List[FusionEntry]):
        """
        Merges the clusters based on the fusion_edge_list computed in _grow_clusters().
        Updates the odd_clusters list by recomputing the neutrality of the newly merged clusters.

        Args:
            fusion_edge_list (List[FusionEntry]): List of edges that connect two
            clusters that was computed in _grow_clusters().
        Returns:
            new_neutral_cluster_roots (List[int]): List of roots of newly neutral clusters
        """
        new_neutral_clusters = []
        for entry in fusion_edge_list:
            root_u, root_v = self.find(entry.u), self.find(entry.v)
            if root_u == root_v:
                continue
            new_root = root_v if self.clusters[root_v].size > self.clusters[root_u].size else root_u
            root_to_update = root_v if new_root == root_u else root_u

            if new_root in new_neutral_clusters or root_to_update in new_neutral_clusters:
                continue

            cluster = self.clusters[new_root]
            other_cluster = self.clusters.pop(root_to_update)

            entry.connecting_edge.data.properties["growth"] = 0
            entry.connecting_edge.data.properties["fully_grown"] = True
            cluster.fully_grown_edges.add(entry.connecting_edge.index)

            # Merge boundaries
            cluster.boundary += other_cluster.boundary
            cluster.boundary.remove(entry.connecting_edge)
            cluster.boundary.remove(entry.connecting_edge.reverse())

            cluster.nodes |= other_cluster.nodes
            cluster.atypical_nodes |= other_cluster.atypical_nodes
            cluster.boundary_nodes |= other_cluster.boundary_nodes
            cluster.fully_grown_edges |= other_cluster.fully_grown_edges
            cluster.size += other_cluster.size

            # update odd_cluster_roots
            if self.code.is_cluster_neutral(
                [self.graph[node] for node in cluster.atypical_nodes]
            ) or self.code.is_cluster_neutral(
                [
                    self.graph[node]
                    for node in cluster.atypical_nodes
                    | (set(list(cluster.boundary_nodes)[:1]) if cluster.boundary_nodes else set())
                ]
            ):
                if new_root in self.odd_cluster_roots:
                    self.odd_cluster_roots.remove(new_root)
                    new_neutral_clusters.append(new_root)
            else:
                if not new_root in self.odd_cluster_roots:
                    self.odd_cluster_roots.append(new_root)

            if root_to_update in self.odd_cluster_roots:
                self.odd_cluster_roots.remove(root_to_update)
            self.graph[root_to_update].properties["root"] = new_root
            self.odd_cluster_roots = sorted(
                self.odd_cluster_roots, key=lambda c: self.clusters[c].size
            )

    def peeling(self, erasure: PyGraph) -> List[int]:
        """ "
        Runs the peeling decoder on the erasure provided.
        Assumes that the erasure is one connected component, if not it will run in an
        infinite loop in the tree construction.
        It works by first producing a spanning forest of the erasure and then
        going backwards through the edges of the tree computing the error based on the syndrome.
        Based on arXiv:1703.01517.

        Args:
            erasure (PyGraph): subgraph of the syndrome graph that represents the erasure.

        Returns:
            errors (List[int]): List of qubit indices on which Pauli errors occurred.
        """
        tree = SpanningForest(vertices={}, edges=[])

        # Construct spanning forest
        # Pick starting vertex
        for vertex in erasure.node_indices():
            if erasure[vertex].is_logical and erasure[vertex].properties["syndrome"]:
                tree.vertices[vertex] = []
                break

        if not tree.vertices:
            for vertex in erasure.node_indices():
                if erasure[vertex].properties["syndrome"]:
                    tree.vertices[vertex] = []
                    break

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
            if erasure[pendant_vertex].properties["syndrome"]:
                edges.add(edge)
                erasure[tree_vertex].properties["syndrome"] = not erasure[tree_vertex].properties[
                    "syndrome"
                ]
                erasure[pendant_vertex].properties["syndrome"] = False

        return [erasure.edges()[edge].qubits[0] for edge in edges if erasure.edges()[edge].qubits]

    def neighbouring_edges(self, node_index) -> List[Tuple[int, int, DecodingGraphEdge]]:
        """Returns all of the neighbouring edges of a node in the decoding graph.

        Args:
            node_index (int): The index of the node in the graph.

        Returns:
            neighbouring_edges (List[Tuple[int, int, DecodingGraphEdge]]): List of neighbouring edges

            In following format::

                {
                    index of edge in graph,
                    index of neighbour node in graph,
                    data payload of the edge
                }

        """
        return [
            (edge, neighbour, data)
            for edge, (_, neighbour, data) in dict(
                self.graph.incident_edge_index_map(node_index)
            ).items()
        ]
    
class TannerUnionFind:
    def __init__(self, adj_mat: np.ndarray) -> None:
        self.adj_mat = adj_mat

    def decode(self, syndrome: np.ndarray):
        checks = syndrome.astype(bool)
        data = np.zeros(self.adj_mat.shape[1], dtype=bool)
        cluster = (checks, data)

        return self.decode_recursive(cluster, syndrome)
    
    def decode_recursive(self, cluster, syndrome):
        ccs = self.connected_components(cluster)
        all_valid = True
        decoded_error = np.zeros_like(cluster[1])
        for cc in ccs:
            valid, x = self.is_valid(cc, syndrome)
            if not valid:
                all_valid = False
                break
            decoded_error |= x
        
        if all_valid:
            return decoded_error

        cluster = self.grow(cluster)
        return self.decode_recursive(cluster, syndrome)     

    def is_valid(self, cluster, syndrome):
        checks, data = cluster
        soe, interior = self.relevant_soe(cluster)
        b = syndrome[checks]
        try:
            x, _, _ = solve2(soe,b)
            tmp_x = np.zeros_like(data)
            tmp_x[interior] = x.astype(bool)
            return True, tmp_x
        except LinAlgError:
            return False, None
        

    def relevant_soe(self, cluster):
        """ajd_mat is the reduced (bipartite) adjaceny matrix (n//2 x n) of the whole tanner (X or Z) graph.
        cluster = (check_selector, data_selector)
        check_selector is a binary array of size n//2, where a 1 means the corresponing check node is in E.
        data_selector is a binary array of size n, where a 1 means the corresponing data node is in E."""
        
        check_selector, data_selector = cluster

        # First we want to find the interior of E. But we only care about the data qubits in Int(E)
        no_outside = self.adj_mat[~check_selector].sum(axis=0) == 0 # selector for all data qubits that do not have connection to any check outside of E
        #print(no_outside)
        interior = data_selector & no_outside
        #print(data_selector)
        #print(interior)
        return self.adj_mat[check_selector][:,interior], interior
    
    def grow(self, cluster):
        checks, data = cluster
        new_checks = self.adj_mat[:, data].sum(axis=1) > 0
        new_data = self.adj_mat[checks].sum(axis=0) > 0
        return (checks|new_checks, data|new_data)
    
    def connected_components(self, cluster):
        checks, data = cluster

        subgraph = self.adj_mat[checks][:, data] # or np._ix
        ms, ns = subgraph.shape
        clusters = []
        checks_accounted = np.zeros(ms, dtype=bool)
        while not np.all(checks_accounted):
            cluster_checks = np.zeros(ms, dtype=bool)
            cluster_data = np.zeros(ns, dtype=bool)
            new_checks = np.zeros_like(cluster_checks)
            new_checks[np.argmax(~checks_accounted)] = True
            while True:
                new_data = (subgraph[new_checks].sum(axis=0) > 0) & ~cluster_data
                cluster_checks |= new_checks
                if not np.any(new_data):
                    break
                new_checks = (subgraph[:, new_data].sum(axis=1) > 0) & ~cluster_checks
                cluster_data |= new_data
                if not np.any(new_checks):
                    break
            checks_accounted |= cluster_checks
            
            tmp_checks = np.zeros_like(cluster[0])
            tmp_checks[checks] = cluster_checks
            tmp_data = np.zeros_like(cluster[1])
            tmp_data[data] = cluster_data
            clusters.append((tmp_checks, tmp_data))
        return clusters
    

class TannerUnionFind_GraphLegacy:
    def __init__(self, code: BBCode) -> None:
        self.code = code

    def extract_syndrome_naive(self, output):
        """
        output is a substring (over two SM cycles) of an actual output of the syndrome measurement circuit
        with 'zx' SM order.
        Returns the X syndrome sx and the Y syndrome sy
        """
        z0, x0, z1, x1 = output[::-1].split(' ')
        sz = (np.array(list(z0)).astype(int) + np.array(list(z1)).astype(int)) % 2
        sx = (np.array(list(x0)).astype(int) + np.array(list(x1)).astype(int)) % 2
        return sx, sz
    
    def extract_syndromes(self, output, round_schedule='zx'):
        """
        output is the whole output string of the syndrome measurement with round schedule 'zx' (otherwise not implemented yet)
        """
    
    def expand(self, tanner, subnodes):
        nodes = []
        for node in subnodes:
            nodes += list(tanner.adj(node).keys())
        nodes += list(subnodes)
        return nodes
    
    def is_valid(self, graph, subnodes):
        subgraph = graph.subgraph(subnodes)
        data = subgraph.filter_nodes(lambda node: node['type']=='data')
        checks = subgraph.filter_nodes(lambda node: node['type']=='check')

        A = adjacency_matrix(subgraph)[np.ix_(checks, data)].astype(int)
        b = np.ones(A.shape[0])

        try:
            x, Ap, bp = solve2(A,b)
            # is a binary vector with indicating errors on this subgraph
            # painfully convert this to indices of the original graph
            x = [subgraph.get_node_data(i)['index'] for i in np.where(x)[0]]
            return True, x, Ap, bp
        except:
            return False, None, None, None

    def interior(self, graph, subnodes):
        _int = []
        for node in subnodes:
            is_in_int = True
            for neighbour in graph.adj(node).keys():
                if neighbour not in subnodes:
                    is_in_int = False
                    break
            if is_in_int:
                _int.append(node)
        return _int
    
    def get_connected_components(self, tanner, E):
        """ Takes as input the full tanner graph (X or Z) and a subset of its nodes.
        Returns the connected components as a list of sets of indices of the original graph."""
        subgraph = tanner.subgraph(E)
        ccs = connected_components(subgraph) # css is a list of sets of indices of the subgraph
        return [[subgraph.get_node_data(sub_index)['index'] for sub_index in cc] for cc in ccs]
    
    def decode(self, output):
        sx, sy = self.extract_syndrome_naive(output)

        E_x = np.where(sx)[0] + self.code.n
        E_z = np.where(sy)[0] + self.code.n

        ex = self.decode_tanner_recursive(self.code.tanner_graph_X, E_x)
        ez = self.decode_tanner_recursive(self.code.tanner_graph_Z, E_z)

        from qiskit_qec.operators.pauli import Pauli

        return Pauli(np.hstack([ex, ez]))
    
    def decode_x(self, output, visual=False):
        sx, sz = self.extract_syndrome_naive(output)
        tanner = self.code.tanner_graph_X
        #initial E
        E = np.where(sx)[0] + self.code.n
        if visual:
            return self.decode_tanner_recursive_visual(tanner, E)
        return self.decode_tanner_recursive(tanner, E)
    
    def decode_z(self, output, visual=False):
        sx, sz = self.extract_syndrome_naive(output)
        tanner = self.code.tanner_graph_Z
        #initial E
        E = np.where(sz)[0] + self.code.n
        if visual:
            return self.decode_tanner_recursive_visual(tanner, E)
        return self.decode_tanner_recursive(tanner, E)

    def decode_tanner_recursive(self, tanner, E):
        """
        Performs one step of the decoding. Get all connected components of E
        Checks if there are invalid ones, if so, calls itself again
        """
        ccs = self.get_connected_components(tanner, E)
        all_valid = True
        e = np.zeros(self.code.n) # initial error
        for cc in ccs:
            #subgraph = tanner.subgraph(cc)
            #fig = graphviz_draw(subgraph, node_attr_fn=lambda node: node['node_attr'])
            interior_ = self.interior(tanner, cc)
            #plot_info_round[2].append(graphviz_draw(tanner.subgraph(interior_), node_attr_fn=lambda node: node['node_attr']))
            valid, x, _, _ = self.is_valid(tanner, interior_)
            # x contains error indices
            e[x] = 1
            #plot_info_round[1].append((fig, valid, x))
            if not valid:
                all_valid = False
                break
        
        #plot_info.append(plot_info_round)

        if all_valid:
            return e#, plot_info
        
        E = self.expand(tanner, E)

        return self.decode_tanner_recursive(tanner, E)
    
    def decode_tanner_recursive_visual(self, tanner, E, plot_info=[]):
        """
        Performs one step of the decoding. Get all connected components of E
        Checks if there are invalid ones, if so, calls itself again
        """
        ccs = self.get_connected_components(tanner, E)
        plot_info_round = (len(ccs), [], [])
        all_valid = True
        e = np.zeros(self.code.n) # initial error
        for cc in ccs:
            subgraph = tanner.subgraph(cc)
            fig = graphviz_draw(subgraph, node_attr_fn=lambda node: node['node_attr'])
            interior_ = self.interior(tanner, cc)
            plot_info_round[2].append(graphviz_draw(tanner.subgraph(interior_), node_attr_fn=lambda node: node['node_attr']))
            valid, x, _, _ = self.is_valid(tanner, interior_)
            # x contains error indices
            e[x] = 1
            plot_info_round[1].append((fig, valid, x))
            if not valid:
                all_valid = False
                break
        
        plot_info.append(plot_info_round)

        if all_valid:
            return e, plot_info
        
        E = self.expand(tanner, E)

        return self.decode_tanner_recursive_visual(tanner, E, plot_info=plot_info)
    