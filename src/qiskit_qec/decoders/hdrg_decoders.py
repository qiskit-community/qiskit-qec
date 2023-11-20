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

from rustworkx import PyGraph, connected_components, distance_matrix

from qiskit_qec.circuits.repetition_code import ArcCircuit
from qiskit_qec.decoders.decoding_graph import DecodingGraph
from qiskit_qec.utils import DecodingGraphEdge, DecodingGraphNode


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

        self.measured_logicals = self.code.measured_logicals()

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
        measurement, corresponding to the logical operators of
        self.measured_logicals.
        """

        # get the list of bulk nodes for each cluster
        cluster_nodes = {c: [] for c in clusters.values()}
        for n, c in clusters.items():
            node = self.decoding_graph.graph[n]
            if not node.is_boundary:
                cluster_nodes[c].append(node)

        # get the list of required logicals for each cluster
        cluster_logicals = {}
        for c, nodes in cluster_nodes.items():
            _, logical_nodes, _ = self.code.check_nodes(nodes, minimal=True)
            z_logicals = [node.qubits[0] for node in logical_nodes]
            cluster_logicals[c] = z_logicals

        # get the net effect on each logical
        net_z_logicals = {z_logical[0]: 0 for z_logical in self.measured_logicals}
        for c, z_logicals in cluster_logicals.items():
            for z_logical in self.measured_logicals:
                if z_logical[0] in z_logicals:
                    net_z_logicals[z_logical[0]] += 1
        for z_logical, num in net_z_logicals.items():
            net_z_logicals[z_logical] = num % 2

        corrected_z_logicals = []
        string = string.split(" ")[0]
        for z_logical in self.measured_logicals:
            raw_logical = int(string[-1 - self.code_index[z_logical[0]]])
            corrected_logical = (raw_logical + net_z_logicals[z_logical[0]]) % 2
            corrected_z_logicals.append(corrected_logical)

        return corrected_z_logicals


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
                con_nodes, ignore_extra_boundary=True
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

    def _get_boundary_nodes(self):
        boundary_nodes = []
        for element, z_logical in enumerate(self.measured_logicals):
            node = DecodingGraphNode(is_boundary=True, qubits=z_logical, index=element)
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
            corrected_logicals (list): A list of integers that are 0 or 1.
        These are the corrected values of the final transversal
        measurement, corresponding to the logical operators of
        self.measured_logicals.
        """

        # turn string into nodes and cluster
        nodes = self.code.string2nodes(string, all_logicals=True)
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
    by the peeling decoder in case of the surface code, or by checking for
    interference with the boundary in case of an abritrary ARC.

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

        if self.use_peeling:
            self.graph = deepcopy(self.decoding_graph.graph)
            highlighted_nodes = self.code.string2nodes(string, all_logicals=True)

            # call cluster to do the clustering, but actually use the peeling form
            self.cluster(highlighted_nodes)
            clusters = self._clusters4peeling

            # determine the net logical z
            net_z_logicals = {tuple(z_logical): 0 for z_logical in self.measured_logicals}
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
            for j, z_logical in enumerate(self.measured_logicals):
                raw_logical = int(raw_logicals[j])
                corrected_logical = (raw_logical + net_z_logicals[tuple(z_logical)]) % 2
                corrected_z_logicals.append(corrected_logical)
            return corrected_z_logicals
        else:
            # turn string into nodes and cluster
            nodes = self.code.string2nodes(string, all_logicals=True)
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
        node_indices = [self.graph.nodes().index(node) for node in nodes]
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

        while self.odd_cluster_roots:
            self._grow_and_merge_clusters()

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
        if not node.is_boundary:
            self.odd_cluster_roots.insert(0, node_index)
        boundary_edges = []
        for edge_index, neighbour, data in self.neighbouring_edges(node_index):
            boundary_edges.append(BoundaryEdge(edge_index, node_index, neighbour, data))
        self.clusters[node_index] = UnionFindDecoderCluster(
            boundary=boundary_edges,
            fully_grown_edges=set(),
            atypical_nodes=set([node_index]) if not node.is_boundary else set([]),
            boundary_nodes=set([node_index]) if node.is_boundary else set([]),
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
                            if self.graph[edge.neighbour_vertex].is_boundary
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
            if erasure[vertex].is_boundary and erasure[vertex].properties["syndrome"]:
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

    def __init__(self, code, decoding_graph: DecodingGraph = None) -> None:
        super().__init__(code, decoding_graph)
        self.graph = deepcopy(self.decoding_graph.graph)
        self.r = 1
        self._clusters4peeling = []

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

        nodes_at_time_zero = []
        for index, node in zip(
            self.decoding_graph.graph.node_indices(), self.decoding_graph.graph.nodes()
        ):
            if node.time == 0 or node.is_boundary:
                nodes_at_time_zero.append(index)
        self.graph = self.decoding_graph.graph.subgraph(nodes_at_time_zero)
        for index, node in zip(self.graph.node_indices(), self.graph.nodes()):
            node.properties["root"] = index
        for edge in self.graph.edges():
            edge.properties["growth"] = 0
            edge.properties["fully_grown"] = False

        string = "".join([str(c) for c in string[::-1]])
        output = [int(bit) for bit in list(string.split(" ", maxsplit=self.code.d)[0])][::-1]
        highlighted_nodes = self.code.string2nodes(string, all_logicals=True)
        if not highlighted_nodes:
            return output  # There's nothing for us to do here

        self.cluster(highlighted_nodes)
        clusters = self._clusters4peeling

        flattened_highlighted_nodes: List[DecodingGraphNode] = []
        for highlighted_node in highlighted_nodes:
            highlighted_node.time = 0
            flattened_highlighted_nodes.append(self.graph.nodes().index(highlighted_node))

        for cluster_nodes, cluster_atypical_nodes in clusters:
            if not cluster_nodes:
                continue
            erasure_graph = deepcopy(self.graph)
            for node in cluster_nodes:
                erasure_graph[node].properties["syndrome"] = node in cluster_atypical_nodes
            erasure = erasure_graph.subgraph(cluster_nodes)
            qubits_to_be_corrected = self.peeling(erasure)
            for idx in qubits_to_be_corrected:
                output[idx] = (output[idx] + 1) % 2

        return output

    def cluster(self, nodes):
        """
        Args:
            nodes (List): List of non-typical nodes in the syndrome graph,
            of the type produced by `string2nodes`.

        Returns:
            clusters (dict): Ddictionary with the indices of
            the given node as keys and an integer specifying their cluster as the corresponding
            value.
        """
        self.clusters: Dict[int, UnionFindDecoderCluster] = {}
        self.odd_cluster_roots = []

        times: List[List[DecodingGraphNode]] = [[] for _ in range(self.code.T + 1)]
        boundaries = []
        for node in deepcopy(nodes):
            if node.is_boundary:
                boundaries.append(node)
            else:
                times[node.time].append(node)
                node.time = 0
        # FIXME: I am not sure when the optimal time to add the boundaries is. Maybe the middle?
        # for node in boundaries:
        times.insert(len(times) // 2, boundaries)

        neutral_clusters = []
        for time in times:
            if not time:
                continue
            for node in time:
                self._add_node(node)
            neutral_clusters += self._collect_neutral_clusters()
            for _ in range(self.r):
                self._grow_and_merge_clusters()
            neutral_clusters += self._collect_neutral_clusters()

        while self.odd_cluster_roots:
            self._grow_and_merge_clusters()

        neutral_clusters += self._collect_neutral_clusters()

        # compile info into standard clusters dict
        clusters = {}
        for c, cluster in enumerate(neutral_clusters):
            # determine which nodes exactly are in the neutral cluster
            neutral_nodes = list(cluster.atypical_nodes | cluster.boundary_nodes)
            # put them in the required dict
            for n in neutral_nodes:
                clusters[n] = c

        neutral_cluster_nodes: List[List[int]] = []
        for cluster in neutral_clusters:
            neutral_cluster_nodes.append(
                (list(cluster.nodes), list(cluster.atypical_nodes | cluster.boundary_nodes))
            )

        self._clusters4peeling = neutral_cluster_nodes

        return neutral_cluster_nodes

    def _add_node(self, node):
        node_index = self.graph.nodes().index(node)
        root = self.find(node_index)
        cluster = self.clusters.get(root)
        if cluster and not node.is_boundary:
            # Add the node to the cluster or remove it if it's already present
            if node_index in cluster.atypical_nodes:
                cluster.atypical_nodes.remove(node_index)
            else:
                cluster.atypical_nodes.add(node_index)
        else:
            self.graph[node_index].properties["root"] = node_index
            self._create_new_cluster(node_index)

    def _collect_neutral_clusters(self):
        neutral_clusters = []
        for root, cluster in self.clusters.copy().items():
            if self.code.is_cluster_neutral(
                [
                    self.graph[node]
                    for node in cluster.atypical_nodes
                    | (set([list(cluster.boundary_nodes)[0]]) if cluster.boundary_nodes else set())
                ]
            ):
                if root in self.odd_cluster_roots:
                    self.odd_cluster_roots.remove(root)
                cluster = self.clusters.pop(root)
                if cluster.atypical_nodes:
                    neutral_clusters.append(cluster)
                for edge in cluster.fully_grown_edges:
                    self.graph.edges()[edge].properties["fully_grown"] = False
                for edge in cluster.boundary:
                    self.graph.edges()[edge.index].properties["growth"] = 0
                for node in cluster.nodes:
                    if self.graph[node].is_boundary:
                        self._create_new_cluster(node)
                    self.graph[node].properties["root"] = node
        return neutral_clusters
