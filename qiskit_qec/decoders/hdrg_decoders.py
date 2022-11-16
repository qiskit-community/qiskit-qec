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

from rustworkx import distance_matrix
from rustworkx import connected_components

from qiskit_qec.decoders.decoding_graph import DecodingGraph


class ClusteringDecoder:
    """Decoder based on finding connected components within the decoding graph."""

    def __init__(
        self,
        code_circuit,
        logical: str,
        decoding_graph: DecodingGraph = None,
    ):
        self.code = code_circuit
        self.logical = logical
        if decoding_graph:
            self.decoding_graph = decoding_graph
        else:
            self.decoding_graph = DecodingGraph(self.code)

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
        for element, z_logical in enumerate(self.code.z_logicals):
            node = {"time": 0, "link qubit": None, "is_boundary": True}
            node["qubits"] = [z_logical]
            node["element"] = element
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
        while ns and dist_max <= len(self.code.links):

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
                    if not dg[n]["is_boundary"]:
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
        self.code.z_logicals.
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
            if not node["is_boundary"]:
                cluster_nodes[c].append(node)

        # get the list of required logicals for each cluster
        cluster_logicals = {}
        for c, nodes in cluster_nodes.items():
            _, z_logicals, _ = code.check_nodes(nodes)
            cluster_logicals[c] = z_logicals

        # get the net effect on each logical
        net_z_logicals = {z_logical: 0 for z_logical in code.z_logicals}
        for c, z_logicals in cluster_logicals.items():
            for z_logical in code.z_logicals:
                if z_logical in z_logicals:
                    net_z_logicals[z_logical] += 1
        for z_logical, num in net_z_logicals.items():
            net_z_logicals[z_logical] = num % 2

        corrected_z_logicals = []
        for z_logical in code.z_logicals:
            raw_logical = int(string[-1 - code.code_index[z_logical]])
            corrected_logical = (raw_logical + net_z_logicals[z_logical]) % 2
            corrected_z_logicals.append(corrected_logical)

        return corrected_z_logicals
