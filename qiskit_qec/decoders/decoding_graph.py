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
Graph used as the basis of decoders.
"""

import itertools
import retworkx as rx
import numpy as np

from qiskit_qec.analysis.faultenumerator import FaultEnumerator


class DecodingGraph:
    """
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction code, and then run suitable decoders.
    """

    def __init__(self, code):
        """
        Args:
            code (CodeCircuit): The QEC code circuit object for which this decoder
                will be used.
        """

        self.code = code

        self.S = self._make_syndrome_graph()

    def _make_syndrome_graph(self):

        S = rx.PyGraph(multigraph=False)

        qc = self.code.circuit["0"]
        fe = FaultEnumerator(qc, method="stabilizer")
        blocks = list(fe.generate_blocks())
        fault_paths = list(itertools.chain(*blocks))

        for _, _, _, output in fault_paths:
            string = "".join([str(c) for c in output[::-1]])
            nodes = self.code.string2nodes(string)
            for node in nodes:
                if node not in S.nodes():
                    S.add_node(node)
            for source in nodes:
                for target in nodes:
                    if target != source:
                        n0 = S.nodes().index(source)
                        n1 = S.nodes().index(target)
                        S.add_edge(n0, n1, 1)

        return S

    def get_error_probs(self, results, logical="0"):
        """
        Generate probabilities of single error events from result counts.
        Args:
            results (dict): A results dictionary.
            logical (string): Logical value whose results are used.
        Returns:
            dict: Keys are the edges for specific error
            events, and values are the calculated probabilities
        Additional information:
            Uses `results` to estimate the probability of the errors that
            create the pairs of nodes specified by the edge.
            Default calculation method is that of Spitz, et al.
            https://doi.org/10.1002/qute.201800012
        """

        shots = sum(results.values())

        neighbours = {}
        av_v = {}
        for n in self.S.node_indexes():
            av_v[n] = 0
            neighbours[n] = []

        av_vv = {}
        av_xor = {}
        for n0, n1 in self.S.edge_list():
            av_vv[n0, n1] = 0
            av_xor[n0, n1] = 0
            neighbours[n0].append(n1)
            neighbours[n1].append(n0)

        for string in results:

            # list of i for which v_i=1
            error_nodes = self.code.string2nodes(string, logical=logical)

            for node0 in error_nodes:
                n0 = self.S.nodes().index(node0)
                av_v[n0] += results[string]
                for n1 in neighbours[n0]:
                    node1 = self.S[n1]
                    if node1 in error_nodes and (n0, n1) in av_vv:
                        av_vv[n0, n1] += results[string]
                    if node1 not in error_nodes:
                        if (n0, n1) in av_xor:
                            av_xor[n0, n1] += results[string]
                        else:
                            av_xor[n1, n0] += results[string]

        for n in self.S.node_indexes():
            av_v[n] /= shots
        for n0, n1 in self.S.edge_list():
            av_vv[n0, n1] /= shots
            av_xor[n0, n1] /= shots

        boundary = []
        error_probs = {}
        for n0, n1 in self.S.edge_list():

            if self.S[n0]["is_logical"]:
                boundary.append(n1)
            elif self.S[n1]["is_logical"]:
                boundary.append(n0)
            else:
                if (1 - 2 * av_xor[n0, n1]) != 0:
                    x = (av_vv[n0, n1] - av_v[n0] * av_v[n1]) / (1 - 2 * av_xor[n0, n1])
                    error_probs[n0, n1] = max(0, 0.5 - np.sqrt(0.25 - x))
                else:
                    error_probs[n0, n1] = np.nan

        prod = {}
        for n0 in boundary:
            for n1 in self.S.node_indexes():
                if n0 != n1:
                    if n0 not in prod:
                        prod[n0] = 1
                    if (n0, n1) in error_probs:
                        prod[n0] *= 1 - 2 * error_probs[n0, n1]
                    elif (n1, n0) in error_probs:
                        prod[n0] *= 1 - 2 * error_probs[n1, n0]

        for n0 in boundary:
            error_probs[n0, n0] = 0.5 + (av_v[n0] - 0.5) / prod[n0]

        return error_probs

    def weight_syndrome_graph(self, results):
        """Generate weighted syndrome graph from result counts.

        Args:
            results (dict): A results dictionary, as produced by the
            `process_results` method of the code.

        Additional information:
            Uses `results` to estimate the probability of the errors that
            create the pairs of nodes in S. The edge weights are then
            replaced with the corresponding -log(p/(1-p).
        """

        error_probs = self.get_error_probs(results)

        for edge in self.S.edge_list():
            p = error_probs[self.S[edge[0]], self.S[edge[1]]]
            if p == 0:
                w = np.inf
            elif 1 - p == 1:
                w = -np.inf
            else:
                w = -np.log(p / (1 - p))
            self.S.update_edge(edge[0], edge[1], w)

    def make_error_graph(self, string):
        """
        Args:
            string (str): A string describing the output from the code.

        Returns:
            E: The subgraph of S which corresponds to the non-trivial
            syndrome elements in the given string.
        """

        E = rx.PyGraph(multigraph=False)
        nodes = self.code.string2nodes(string, all_logicals=True)
        for node in nodes:
            if node not in E.nodes():
                E.add_node(node)

        # for each pair of nodes in error create an edge and weight with the
        # distance
        distance_matrix = rx.graph_floyd_warshall_numpy(self.S, weight_fn=float)

        for source_index in E.node_indexes():
            for target_index in E.node_indexes():
                source = E[source_index]
                target = E[target_index]
                if target != source:
                    distance = int(
                        distance_matrix[self.S.nodes().index(source)][self.S.nodes().index(target)]
                    )
                    E.add_edge(source_index, target_index, -distance)
        return E
