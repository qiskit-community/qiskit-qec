# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2019.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# pylint: disable=invalid-name

"""
Decoders for quantum error correction codes, with a focus on those that can be
expressed as solving a graph-theoretic problem.
"""


import warnings
import retworkx as rx


class GraphDecoder:
    """
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction code, and then run suitable decoders.
    """

    def __init__(self, graph):
        """
        Args:
            graph (GraphDecoder): The decoder graph on which decoding is performed.
        """

        self.graph = graph

    def matching(self, string):
        """
        Args:
            string (str): A string describing the output from the code.

        Returns:
            str: Corrected values for logical operators that correspond to nodes.

        Additional information:
            This function can be run directly, or used indirectly to
            calculate a logical error probability with `get_logical_prob`
        """

        # this matching algorithm is designed for a single graph
        E = self.graph.make_error_graph(string)

        # set up graph that is like E, but each syndrome node is connected to a
        # separate copy of the nearest logical node
        E_matching = rx.PyGraph(multigraph=False)
        syndrome_nodes = []
        logical_nodes = []
        logical_neighbours = []
        for node in E.nodes():
            E_matching.add_node(node)
            if node["is_boundary"]:
                logical_nodes.append(node)
            else:
                syndrome_nodes.append(node)
        for source in syndrome_nodes:
            n0 = E.nodes().index(source)
            for target in syndrome_nodes:
                n1 = E.nodes().index(target)
                if target != source:
                    E_matching.add_edge(
                        n0,
                        n1,
                        E.get_edge_data(n0, n1),
                    )

            potential_logical = {}
            for target in logical_nodes:
                n1 = E.nodes().index(target)
                potential_logical[n1] = E.get_edge_data(n0, n1)
            nearest_logical = max(potential_logical, key=potential_logical.get)
            nl_target = E[nearest_logical].copy()
            nl_target["connected_to"] = source.copy()
            if nl_target not in E_matching.nodes():
                E_matching.add_node(nl_target)
            E_matching.add_edge(
                n0,
                E_matching.nodes().index(nl_target),
                potential_logical[nearest_logical],
            )
            logical_neighbours.append(nl_target)
        for source in logical_neighbours:
            for target in logical_neighbours:
                if target != source:
                    n0 = E_matching.nodes().index(source)
                    n1 = E_matching.nodes().index(target)
                    E_matching.add_edge(n0, n1, 0)
        # do the matching on this
        matches = rx.max_weight_matching(E_matching, max_cardinality=True, weight_fn=lambda x: x)
        # use it to construct and return a corrected logical string
        logicals = self.graph.code.string2raw_logicals(string)
        for (n0, n1) in matches:
            source = E_matching[n0]
            target = E_matching[n1]
            sil = E_matching[n0]["is_boundary"]
            til = E_matching[n1]["is_boundary"]
            if sil and not til:
                elem = E_matching[n0]["element"]
                logicals[elem] = str((int(logicals[elem]) + 1) % 2)
            if til and not sil:
                elem = E_matching[n1]["element"]
                logicals[elem] = str((int(logicals[elem]) + 1) % 2)

        return logicals

    def get_logical_prob(self, results, logical="0", algorithm="matching"):
        """
        Args:
            results (dict): A results dictionary from running a circuit
            of the code.
            logical (str): Encoded logical value at readout.
            algorithm (str): Choice of which decoder to use.

        Returns:
            dict: Dictionary of logical error probabilities for
            each of the encoded logical states whose results were given in
            the input.
        """

        shots = 0
        incorrect_shots = 0

        corrected_results = {}
        if algorithm == "matching":
            for string in results:
                corr_str = self.matching(string)[0]
                if corr_str in corrected_results:
                    corrected_results[corr_str] += results[string]
                else:
                    corrected_results[corr_str] = results[string]
        else:
            warnings.warn(
                "The requested algorithm " + str(algorithm) + " is not known.",
                Warning,
            )

        for string, samples in corrected_results.items():
            shots += samples
            if string != str(logical):
                incorrect_shots += samples

        logical_prob = incorrect_shots / shots

        return logical_prob
