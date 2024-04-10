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
from time import perf_counter
from typing import Dict, List, Set, Tuple

import numpy as np
from rustworkx import PyGraph, connected_components, distance_matrix, adjacency_matrix
from rustworkx.visualization import graphviz_draw

from qiskit_qec.codes.bb_code import BBCode
from qiskit_qec.decoders.decoding_graph import DecodingGraph
#from qiskit_qec.linear.matrix import solve2, LinAlgError
from qiskit_qec.analysis.lse_solvers import solve
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
    
class TannerUnionFindOptimized:
    def __init__(self, adj_mat: np.ndarray) -> None:
        """ Initializes the Tanner graph (simply with reduced adjacency matrix) """
        self.adj_mat = adj_mat
        self.decoders = {}

    class Cluster:
        def __init__(self, checks_interior, data_interior, cluster_surface, cluster_surface_type, valid, internal_error) -> None:
            self.checks_interior = checks_interior
            self.data_interior = data_interior
            self.cluster_surface = cluster_surface
            self.cluster_surface_type = cluster_surface_type
            self.valid = valid
            self.internal_error = internal_error

        @property
        def checks(self):
            """ return all checks """
            if self.cluster_surface_type == 'c':
                return self.checks_interior | self.cluster_surface
            return self.checks_interior
        
        @property
        def data(self):
            """ return all data qubits """
            if self.cluster_surface_type == 'd':
                return self.data_interior | self.cluster_surface
            return self.data_interior

        def touches(self, other: "TannerUnionFindOptimized.Cluster") -> bool:
            return self.cluster_surface_type == other.cluster_surface_type and np.any(self.cluster_surface & other.cluster_surface)
        
        def merge(self, other: "TannerUnionFindOptimized.Cluster") -> None:
            """ merging happens in place. This assumes that touches is true """
            self.checks_interior |= other.checks_interior
            self.data_interior |= other.data_interior
            self.cluster_surface |= other.cluster_surface
            # self.cluster_surface_type not affected
            self.valid = False
            self.internal_error = None

        def size(self):
            return self.checks.sum(), self.data.sum()
        
        def __repr__(self) -> str:
            string = f'checks: {list(np.where(self.checks_interior)[0])}\n'
            string += f'data: {list(np.where(self.data_interior)[0])}\n'
            string += f'surface: {list(np.where(self.cluster_surface)[0])}, type: {self.cluster_surface_type}\n'
            if self.valid:
                string += f'internal error: {list(np.where(self.internal_error)[0])}\n'
            else:
                string += f'Not valid'
            return string
        
    class Decoder:
        def __init__(self, adj_mat: np.ndarray, T: int) -> None:
            m, n = adj_mat.shape
            if T > 1: # in this case (more than two consecutive SM rounds) assume with measurement error
                self.adj_mat_T = np.zeros((T*m, (T+1)*m + T*n), dtype=bool) # reduced adjacency matrix of T-foliated Tanner graph
                for i in range(T):
                    self.adj_mat_T[i*m:i*m+m, i*(m+n): i*(m+n) + m] = np.eye(m, dtype=bool)
                    self.adj_mat_T[i*m:i*m+m, i*(m+n) + m: i*(m+n) + m+n] = adj_mat
                    self.adj_mat_T[i*m:i*m+m, i*(m+n) + m+n: i*(m+n) + 2*m+n] = np.eye(m, dtype=bool)
            else:
                self.adj_mat_T = adj_mat # otherwise assume no measurement errors -> decode on unfoliated (in particular no boundary virtual qubits) Tanner graph
            self.m = m
            self.n = n
            self.T = T
            #self.adj_mat_T = adj_mat

        def predecode(self, syndrome):
            pre_error = self.adj_mat_T.T @ syndrome == 3
            pre_syndrome = self.adj_mat_T @ pre_error % 2 ^ syndrome
            return pre_error, pre_syndrome

        def decode(self, clusters: List["TannerUnionFindOptimized.Cluster"], syndrome, history=None, lse_solver=None):
            if history is None:# benchmarking
                history = []# benchmarking

            step_history = {}# benchmarking
            step_history['num_clust'] = len(clusters)# benchmarking

            cluster_figs = []  # benchmarking
            all_valid = True
            decoded_error = np.zeros(self.adj_mat_T.shape[1], dtype=np.uint8)
            #new_clusters = []
            for cluster in clusters:
                if cluster.valid:
                    decoded_error |= cluster.internal_error
                    continue

                valid, x, figs = self.is_valid(cluster, syndrome, lse_solver=lse_solver)  # benchmarking
                cluster_figs.append(figs)  # benchmarking
                if not valid:
                    all_valid = False
                else:
                    decoded_error |= x
                    cluster.valid = True
                    cluster.internal_error = x

            step_history['clusters'] =  cluster_figs# benchmarking
            
            if all_valid:
                step_history['t_grow'] = 0
                #figs = (size, t_con, num_clust, cluster_figs, None) # benchmarking
                history.append(step_history)
                #history.append(figs) # benchmarking
                return decoded_error, history # benchmarking
            t0 = perf_counter() # benchmarking
            clusters = self.grow(clusters)
            t_grow = perf_counter() - t0 # benchmarking
            #figs = (size, t_con, num_clust, cluster_figs, t_grow) # benchmarking
            step_history['t_grow'] = t_grow# benchmarking
            history.append(step_history) # benchmarking
            return self.decode(clusters, syndrome, history=history, lse_solver=lse_solver)      # benchmarking

        def is_valid(self, cluster: "TannerUnionFindOptimized.Cluster", syndrome, lse_solver=None):
            # checks_interior, data_interior, cluster_surface, cluster_surface_type, valid, internal_error = cluster
            # if cluster_surface_type == 'c':
            #     checks = checks_interior | cluster_surface
            #     data = data_interior
            # else:
            #     checks = checks_interior
            #     data = data_interior | cluster_surface

            clust_size = cluster.size() # benchmarking

            t0 = perf_counter() # benchmarking
            a, true_interior_data = self.get_cluster_LSE(cluster) # benchmarking
            b = syndrome[cluster.checks]
            t_int = perf_counter() - t0 # t_int really doesn't say much, just means how long to copy the array

            int_size = true_interior_data.sum() # benchmarking

            #t0 = perf_counter() # benchmarking

            valid, x, solve_stats = solve(a, b, minimize_weight=True, stats=True, lse_solver=lse_solver)
            if valid:
            # x is sparse, contains true indices of error within interior
                int_x = np.zeros_like(cluster.data)
                int_x[true_interior_data] = x
            else: int_x = None
            # now int_x contains the interior error in dense format with respect to whole graph
            
            #t_valid = perf_counter() - t0 # benchmarking

            stats = {'valid': valid, 'clust_size': clust_size, 'int_size': int_size, 't_int': t_int}
            for solve_stat in solve_stats:
                stats[solve_stat] = solve_stats[solve_stat]

            # if figs is not None: # benchmarking
            #     t_valid -= (figs[0] + figs[2]) # benchmarking
            #     figs = (clust_size, t_int, int_size, t_valid, *figs) # benchmarking
            # else:
            #     figs = (clust_size, t_int, int_size, t_valid, None, None, None)
            
            return valid, int_x, stats
            

        def get_cluster_LSE(self, cluster: "TannerUnionFindOptimized.Cluster"):
            """ajd_mat is the reduced (bipartite) adjaceny matrix (n//2 x n) of the whole tanner (X or Z) graph.
            cluster = (check_selector, data_selector)
            check_selector is a binary array of size n//2, where a 1 means the corresponing check node is in E.
            data_selector is a binary array of size n, where a 1 means the corresponing data node is in E."""
            # checks_interior, data_interior, cluster_surface, cluster_surface_type, valid, internal_error = cluster
            # if cluster_surface_type == 'c':
            #     check_selector = checks_interior | cluster_surface
            #     data_selector = data_interior
            # else:
            #     check_selector = checks_interior
            #     data_selector = data_interior | cluster_surface

            # First we want to find the interior of E. But we only care about the data qubits in Int(E)
            no_outside = self.adj_mat_T[~cluster.checks].sum(axis=0) == 0 # selector for all data qubits that do not have connection to any check outside of E
            interior = cluster.data & no_outside
            return self.adj_mat_T[cluster.checks][:,interior], interior
        
        def grow(self, clusters: List["TannerUnionFindOptimized.Cluster"], grow_valid=False):
            grown_clusters = []
            for cluster in clusters:
                #checks_interior, data_interior, cluster_surface, cluster_surface_type, valid, internal_error = cluster
                if not cluster.valid or grow_valid:
                    if cluster.cluster_surface_type == 'c':
                        cluster.checks_interior |= cluster.cluster_surface
                        cluster.cluster_surface = np.any(self.adj_mat_T[cluster.cluster_surface], axis=0) &~cluster.data_interior
                        cluster.cluster_surface_type = 'd'
                    elif cluster.cluster_surface_type == 'd':
                        cluster.data_interior |= cluster.cluster_surface
                        cluster.cluster_surface =  np.any(self.adj_mat_T[:, cluster.cluster_surface], axis=1) &~cluster.checks_interior
                        cluster.cluster_surface_type = 'c'

                # merge directly
                merged_away = []
                for i, other in enumerate(grown_clusters):
                    if cluster.touches(other): # these 2 clusters are now connected
                        # keep track of merged
                        merged_away.append(i)
                        cluster.merge(other) # merges other into cluster in place

                for i in merged_away[::-1]: # delete the merged clusters in reverse order so indices keep making sense
                    del grown_clusters[i]
                # append the merged cluster
                grown_clusters.append(cluster)
            return grown_clusters

    def get_decoder(self, T: int) -> "TannerUnionFindOptimized.Decoder":
        if T not in self.decoders:
            self.decoders[T] = TannerUnionFindOptimized.Decoder(self.adj_mat, T)
        return self.decoders[T]

    def decode(self, syndrome: np.ndarray, predecode=False, lse_solver=None) -> Tuple[np.ndarray, List[Tuple[Tuple[int,int], float, int, List[Tuple[Tuple[int,int], float, int, float, float, int, float]], float]]]:
        """ Main entry point for decoding, takes syndrome and initialized clusters. In optimized version actually creates the clusters. """

        "Syndrome can be from multiple measurement rounds"
        m, n = self.adj_mat.shape
        l = len(syndrome)
        # syndrome length should be T*m where T+1 is number of SM rounds (i.e T number of parity comparisons between them)
        if l % m != 0:
            raise ValueError("Syndrome not of valid length (valid <=> multiple of number of checks)")
        
        T = l // m

        decoder = self.get_decoder(T)

        if predecode:
            pre_error, syndrome = decoder.predecode(syndrome=syndrome)

        clusters = []
        for check in np.where(syndrome)[0]:
            checks = np.zeros(T*m, dtype=bool)
            if T > 1:
                data = np.zeros((T+1)*m + T*n, dtype=bool)
            else:
                data = np.zeros(n, dtype=bool)
            cluster_surface = np.zeros(T*m, dtype=bool)
            cluster_surface[check] = True
            cluster_surface_type = 'c'
            valid = False
            internal_error = None
            cluster = TannerUnionFindOptimized.Cluster(checks, data, cluster_surface, cluster_surface_type, valid, internal_error)
            clusters.append(cluster)


        decoded_error, history = decoder.decode(clusters, syndrome, lse_solver=lse_solver)
        if predecode:
            decoded_error ^= pre_error

        return decoded_error, history
