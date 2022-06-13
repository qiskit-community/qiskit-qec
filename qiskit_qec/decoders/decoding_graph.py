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
from typing import List, Tuple

import itertools
import logging
import numpy as np
import retworkx as rx

from qiskit_qec.analysis.faultenumerator import FaultEnumerator


class DecodingGraph:
    """
    Class to construct the decoding graph for the code given by a CodeCircuit object,
    for use in a decoder.
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

        self.graph = self._make_syndrome_graph()

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
                        qubits = []
                        if not (source["is_boundary"] and target["is_boundary"]):
                            qubits = list(set(source["qubits"]).intersection(target["qubits"]))
                        if source["time"] != target["time"] and len(qubits) > 1:
                            qubits = []
                        edge = {"qubits": qubits, "weight": 1}
                        S.add_edge(n0, n1, edge)

        return S

    def get_error_probs(self, results, logical="0"):
        """Generate probabilities of single error events from result counts.

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
        for n in self.graph.node_indexes():
            av_v[n] = 0
            neighbours[n] = []

        av_vv = {}
        av_xor = {}
        for n0, n1 in self.graph.edge_list():
            av_vv[n0, n1] = 0
            av_xor[n0, n1] = 0
            neighbours[n0].append(n1)
            neighbours[n1].append(n0)

        for string in results:

            # list of i for which v_i=1
            error_nodes = self.code.string2nodes(string, logical=logical)

            for node0 in error_nodes:
                n0 = self.graph.nodes().index(node0)
                av_v[n0] += results[string]
                for n1 in neighbours[n0]:
                    node1 = self.graph[n1]
                    if node1 in error_nodes and (n0, n1) in av_vv:
                        av_vv[n0, n1] += results[string]
                    if node1 not in error_nodes:
                        if (n0, n1) in av_xor:
                            av_xor[n0, n1] += results[string]
                        else:
                            av_xor[n1, n0] += results[string]

        for n in self.graph.node_indexes():
            av_v[n] /= shots
        for n0, n1 in self.graph.edge_list():
            av_vv[n0, n1] /= shots
            av_xor[n0, n1] /= shots

        boundary = []
        error_probs = {}
        for n0, n1 in self.graph.edge_list():

            if self.graph[n0]["is_boundary"]:
                boundary.append(n1)
            elif self.graph[n1]["is_boundary"]:
                boundary.append(n0)
            else:
                if (1 - 2 * av_xor[n0, n1]) != 0:
                    x = (av_vv[n0, n1] - av_v[n0] * av_v[n1]) / (1 - 2 * av_xor[n0, n1])
                    error_probs[n0, n1] = max(0, 0.5 - np.sqrt(0.25 - x))
                else:
                    error_probs[n0, n1] = np.nan

        prod = {}
        for n0 in boundary:
            for n1 in self.graph.node_indexes():
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

        for edge in self.graph.edge_list():
            p = error_probs[self.graph[edge[0]], self.graph[edge[1]]]
            if p == 0:
                w = np.inf
            elif 1 - p == 1:
                w = -np.inf
            else:
                w = -np.log(p / (1 - p))
            self.graph.update_edge(edge[0], edge[1], w)

    def make_error_graph(self, string: str):
        """Returns error graph.

        Args:
            string (str): A string describing the output from the code.

        Returns:
            The subgraph of S which corresponds to the non-trivial
            syndrome elements in the given string.
        """

        E = rx.PyGraph(multigraph=False)
        nodes = self.code.string2nodes(string, all_logicals=True)
        for node in nodes:
            if node not in E.nodes():
                E.add_node(node)

        # for each pair of nodes in error create an edge and weight with the
        # distance
        def weight_fn(edge):
            return int(edge["weight"])

        distance_matrix = rx.graph_floyd_warshall_numpy(self.graph, weight_fn=weight_fn)

        for source_index in E.node_indexes():
            for target_index in E.node_indexes():
                source = E[source_index]
                target = E[target_index]
                if target != source:
                    distance = int(
                        distance_matrix[self.graph.nodes().index(source)][
                            self.graph.nodes().index(target)
                        ]
                    )
                    E.add_edge(source_index, target_index, -distance)
        return E


class CSSDecodingGraph:
    """
    Class to construct the decoding graph required for the CircuitModelMatchingDecoder
    for a generic CSS code.
    """

    def __init__(
        self,
        css_x_gauge_ops: List[Tuple[int]],
        css_x_stabilizer_ops: List[Tuple[int]],
        css_x_boundary: List[Tuple[int]],
        css_z_gauge_ops: List[Tuple[int]],
        css_z_stabilizer_ops: List[Tuple[int]],
        css_z_boundary: List[Tuple[int]],
        blocks: int,
        round_schedule: str,
        basis: str,
    ):

        self.css_x_gauge_ops = css_x_gauge_ops
        self.css_x_stabilizer_ops = css_x_stabilizer_ops
        self.css_x_boundary = css_x_boundary
        self.css_z_gauge_ops = css_z_gauge_ops
        self.css_z_stabilizer_ops = css_z_stabilizer_ops
        self.css_z_boundary = css_z_boundary
        self.blocks = blocks
        self.round_schedule = round_schedule
        self.basis = basis

        self.layer_types = self._layer_types(self.blocks, self.round_schedule, self.basis)

        self._decoding_graph()

    def _layer_types(self, blocks: int, round_schedule: str, basis: str) -> List[str]:
        """Return a list of decoding graph layer types.

        The entries are 'g' for gauge and 's' for stabilizer.
        """
        layer_types = []
        last_step = basis
        for _ in range(blocks):
            for step in round_schedule:
                if basis == "z" and step == "z" and last_step == "z":
                    layer_types.append("g")
                elif basis == "z" and step == "z" and last_step == "x":
                    layer_types.append("s")
                elif basis == "x" and step == "x" and last_step == "x":
                    layer_types.append("g")
                elif basis == "x" and step == "x" and last_step == "z":
                    layer_types.append("s")
                last_step = step
        if last_step == basis:
            layer_types.append("g")
        else:
            layer_types.append("s")
        return layer_types

    def _decoding_graph(self):
        """Construct the decoding graph for the given basis.

        This method sets edge weights all to 1 and is based on
        computing intersections of operator supports.

        Returns a tuple (idxmap, node_layers, G)
        where idxmap is a dict
        mapping tuples (t, qubit_set) to integer vertex indices in the
        decoding graph G. The list node_layers contains lists of nodes
        for each time step.
        """
        graph = rx.PyGraph(multigraph=False)
        gauges = []
        stabilizers = []
        boundary = []
        if self.basis == "z":
            gauges = self.css_z_gauge_ops
            stabilizers = self.css_z_stabilizer_ops
            boundary = self.css_z_boundary
        elif self.basis == "x":
            gauges = self.css_x_gauge_ops
            stabilizers = self.css_x_stabilizer_ops
            boundary = self.css_x_boundary

        # Construct the decoding graph
        idx = 0  # vertex index counter
        idxmap = {}  # map from vertex data (t, qubits) to vertex index
        node_layers = []
        for time, layer in enumerate(self.layer_types):
            # Add vertices at time t
            node_layer = []
            if layer == "g":
                all_z = gauges
            elif layer == "s":
                all_z = stabilizers
            for supp in all_z:
                node = {"time": time, "qubits": supp, "highlighted": False}
                graph.add_node(node)
                logging.debug("node %d t=%d %s", idx, time, supp)
                idxmap[(time, tuple(supp))] = idx
                node_layer.append(idx)
                idx += 1
            for supp in boundary:
                # Add optional is_boundary property for pymatching
                node = {"time": time, "qubits": supp, "highlighted": False, "is_boundary": True}
                graph.add_node(node)
                logging.debug("boundary %d t=%d %s", idx, time, supp)
                idxmap[(time, tuple(supp))] = idx
                node_layer.append(idx)
                idx += 1
            node_layers.append(node_layer)
            if layer == "g":
                all_z = gauges + boundary
            elif layer == "s":
                all_z = stabilizers + boundary
            # Add space-like edges at time t
            # The qubit sets of any pair of vertices at time
            # t can intersect on multiple qubits.
            # If they intersect, we add an edge and label it by
            # one of the common qubits. This makes an assumption
            # that the intersection operator is equivalent to a single
            # qubit operator modulo the gauge group.
            # Space-like edges do not correspond to syndrome errors, so the
            # syndrome property is an empty list.
            for i, op_g in enumerate(all_z):
                for j in range(i + 1, len(all_z)):
                    op_h = all_z[j]
                    com = list(set(op_g).intersection(set(op_h)))
                    if -1 in com:
                        com.remove(-1)
                    if len(com) > 0:
                        # Include properties for use with pymatching:
                        # qubit_id is an integer or set of integers
                        # weight is a floating point number
                        # error_probability is a floating point number
                        edge = {
                            "qubits": [com[0]],
                            "measurement_error": 0,
                            "weight": 1,
                            "highlighted": False,
                        }
                        graph.add_edge(
                            idxmap[(time, tuple(op_g))], idxmap[(time, tuple(op_h))], edge
                        )
                        logging.debug("spacelike t=%d (%s, %s)", time, op_g, op_h)
                        logging.debug(
                            " qubits %s",
                            [com[0]],
                        )

            # Add boundary space-like edges
            for i in range(len(boundary) - 1):
                bound_g = boundary[i]
                bound_h = boundary[i + 1]
                # Include properties for use with pymatching:
                # qubit_id is an integer or set of integers
                # weight is a floating point number
                # error_probability is a floating point number
                edge = {
                    "qubits": [],
                    "measurement_error": 0,
                    "weight": 0,
                    "highlighted": False,
                }
                graph.add_edge(idxmap[(time, tuple(bound_g))], idxmap[(time, tuple(bound_h))], edge)
                logging.debug("spacelike boundary t=%d (%s, %s)", time, bound_g, bound_h)

            # Add (space)time-like edges from t to t-1
            # By construction, the qubit sets of pairs of vertices at S and T
            # at times t-1 and t respectively
            # either (a) contain each other (S subset T or T subset S) and
            # |S|,|T|>1,
            # (b) intersect on one or more qubits, or (c) are disjoint.
            # In case (a), we add an edge that corresponds to a syndrome bit
            # error at time t-1.
            # In case (b), we add an edge that corresponds to a spacetime hook
            # error, i.e. a syndrome bit error at time t-1
            # together with an error on one of the common qubits. Again
            # this makes an assumption that all such errors are equivalent.
            # In case (c), we do not add an edge.
            # Important: some space-like hooks are not accounted for.
            # They can have longer paths between non-intersecting operators.
            # We will account for these in _revise_decoding_graph if needed.
            if time > 0:
                current_sets = gauges
                prior_sets = gauges
                if self.layer_types[time] == "s":
                    current_sets = stabilizers
                if self.layer_types[time - 1] == "s":
                    prior_sets = stabilizers
                for op_g in current_sets:
                    for op_h in prior_sets:
                        com = list(set(op_g).intersection(set(op_h)))
                        if -1 in com:
                            com.remove(-1)
                        if len(com) > 0:  # not Case (c)
                            # Include properties for use with pymatching:
                            # qubit_id is an integer or set of integers
                            # weight is a floating point number
                            # error_probability is a floating point number
                            # Case (a)
                            if set(com) == set(op_h) or set(com) == set(op_g):
                                edge = {
                                    "qubits": [],
                                    "measurement_error": 1,
                                    "weight": 1,
                                    "highlighted": False,
                                }
                                graph.add_edge(
                                    idxmap[(time - 1, tuple(op_h))],
                                    idxmap[(time, tuple(op_g))],
                                    edge,
                                )
                                logging.debug("timelike t=%d (%s, %s)", time, op_g, op_h)
                            else:  # Case (b)
                                edge = {
                                    "qubits": [com[0]],
                                    "measurement_error": 1,
                                    "weight": 1,
                                    "highlighted": False,
                                }
                                graph.add_edge(
                                    idxmap[(time - 1, tuple(op_h))],
                                    idxmap[(time, tuple(op_g))],
                                    edge,
                                )
                                logging.debug("spacetime hook t=%d (%s, %s)", time, op_g, op_h)
                                logging.debug(" qubits %s", [com[0]])
                # Add a single time-like edge between boundary vertices at
                # time t-1 and t
                edge = {
                    "qubits": [],
                    "measurement_error": 0,
                    "weight": 0,
                    "highlighted": False,
                }
                graph.add_edge(
                    idxmap[(time - 1, tuple(boundary[0]))], idxmap[(time, tuple(boundary[0]))], edge
                )
                logging.debug("boundarylink t=%d", time)

        self.idxmap = idxmap
        self.node_layers = node_layers
        self.graph = graph
