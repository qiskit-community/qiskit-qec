# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2019-2023.
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
import logging
from typing import List, Tuple

import numpy as np
import rustworkx as rx
from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.exceptions import QiskitQECError
from qiskit_qec.utils import DecodingGraphNode, DecodingGraphEdge


class DecodingGraph:
    """
    Class to construct the decoding graph for the code given by a CodeCircuit object,
    for use in a decoder.
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction code, and then run suitable decoders.
    """

    METHOD_SPITZ: str = "spitz"
    METHOD_NAIVE: str = "naive"
    AVAILABLE_METHODS = {METHOD_SPITZ, METHOD_NAIVE}

    def __init__(self, code, brute=False, graph=None):
        """
        Args:
            code (CodeCircuit): The QEC code circuit object for which this decoding
                graph will be created. If None, graph will initialized as empty.
            brute (bool): Whether to create the graph by analysing the circuits,
            or to use a helper method from the code class (if available).
        """

        self.code = code
        self.brute = brute

        if graph:
            self.graph = graph
        else:
            self._make_syndrome_graph()

    def _make_syndrome_graph(self):
        if not self.brute and hasattr(self.code, "_make_syndrome_graph"):
            self.graph, self.hyperedges = self.code._make_syndrome_graph()
        else:
            graph = rx.PyGraph(multigraph=False)
            self.hyperedges = []

            if self.code is not None:
                # get the circuit used as the base case
                if isinstance(self.code.circuit, dict):
                    if "base" not in dir(self.code):
                        base = "0"
                    else:
                        base = self.code.base
                    qc = self.code.circuit[base]
                else:
                    qc = self.code.circuit

                fe = FaultEnumerator(qc, method="stabilizer")
                blocks = list(fe.generate_blocks())
                fault_paths = list(itertools.chain(*blocks))

                for _, _, _, output in fault_paths:
                    string = "".join([str(c) for c in output[::-1]])
                    nodes = self.code.string2nodes(string)
                    for node in nodes:
                        if node not in graph.nodes():
                            graph.add_node(node)
                    hyperedge = {}
                    for source in nodes:
                        for target in nodes:
                            if target != source:
                                n0 = graph.nodes().index(source)
                                n1 = graph.nodes().index(target)
                                qubits = []
                                if not (source.is_boundary and target.is_boundary):
                                    qubits = list(set(source.qubits).intersection(target.qubits))
                                    if not qubits:
                                        continue
                                if (
                                    source.time != target.time
                                    and len(qubits) > 1
                                    and not source.is_boundary
                                    and not target.is_boundary
                                ):
                                    qubits = []
                                edge = DecodingGraphEdge(qubits, 1)
                                graph.add_edge(n0, n1, edge)
                                if (n1, n0) not in hyperedge:
                                    hyperedge[n0, n1] = edge
                    if hyperedge and hyperedge not in self.hyperedges:
                        self.hyperedges.append(hyperedge)

            self.graph = graph

    def get_error_probs(
        self, counts, logical: str = "0", method: str = METHOD_SPITZ
    ) -> List[Tuple[Tuple[int, int], float]]:
        """
        Generate probabilities of single error events from result counts.

        Args:
            counts (dict): Counts dictionary of the results to be analyzed.
            logical (string): Logical value whose results are used.
            method (string): Method to used for calculation. Supported
            methods are 'spitz' (default) and 'naive'.
        Returns:
            dict: Keys are the edges for specific error
            events, and values are the calculated probabilities.
        Additional information:
            Uses `counts` to estimate the probability of the errors that
            create the pairs of nodes specified by the edge.
            Default calculation method is that of Spitz, et al.
            https://doi.org/10.1002/qute.201800012
        """

        shots = sum(counts.values())

        if method not in self.AVAILABLE_METHODS:
            raise QiskitQECError("fmethod {method} is not supported.")

        # method for edges
        if method == self.METHOD_SPITZ:
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

            for string in counts:
                # list of i for which v_i=1
                error_nodes = self.code.string2nodes(string, logical=logical)

                for node0 in error_nodes:
                    n0 = self.graph.nodes().index(node0)
                    av_v[n0] += counts[string]
                    for n1 in neighbours[n0]:
                        node1 = self.graph[n1]
                        if node1 in error_nodes and (n0, n1) in av_vv:
                            av_vv[n0, n1] += counts[string]
                        if node1 not in error_nodes:
                            if (n0, n1) in av_xor:
                                av_xor[n0, n1] += counts[string]
                            else:
                                av_xor[n1, n0] += counts[string]

            for n in self.graph.node_indexes():
                av_v[n] /= shots
            for n0, n1 in self.graph.edge_list():
                av_vv[n0, n1] /= shots
                av_xor[n0, n1] /= shots

            boundary = []
            error_probs = {}
            for n0, n1 in self.graph.edge_list():
                if self.graph[n0].is_boundary:
                    boundary.append(n1)
                elif self.graph[n1].is_boundary:
                    boundary.append(n0)
                else:
                    if (1 - 2 * av_xor[n0, n1]) != 0:
                        x = (av_vv[n0, n1] - av_v[n0] * av_v[n1]) / (1 - 2 * av_xor[n0, n1])
                        if x < 0.25:
                            error_probs[n0, n1] = max(0, 0.5 - np.sqrt(0.25 - x))
                        else:
                            error_probs[n0, n1] = np.nan
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

        # generally applicable but approximate method
        elif method == self.METHOD_NAIVE:
            # for every edge in the graph, we'll determine the histogram
            # of whether their nodes are in the error nodes
            count = {
                edge: {element: 0 for element in ["00", "01", "10", "11"]}
                for edge in self.graph.edge_list()
            }
            for string in counts:
                error_nodes = self.code.string2nodes(string, logical=logical)
                for edge in self.graph.edge_list():
                    element = ""
                    for j in range(2):
                        if self.graph[edge[j]] in error_nodes:
                            element += "1"
                        else:
                            element += "0"
                    count[edge][element] += counts[string]

            # ratio of error on both to error on neither is, to first order,
            # p/(1-p) where p is the prob to be determined.
            error_probs = {}
            for n0, n1 in self.graph.edge_list():
                if count[n0, n1]["00"] > 0:
                    ratio = count[n0, n1]["11"] / count[n0, n1]["00"]
                else:
                    ratio = np.nan
                p = ratio / (1 + ratio)
                if self.graph[n0].is_boundary and not self.graph[n1].is_boundary:
                    edge = (n1, n1)
                elif not self.graph[n0].is_boundary and self.graph[n1].is_boundary:
                    edge = (n0, n0)
                else:
                    edge = (n0, n1)
                error_probs[edge] = p

        return error_probs

    def weight_syndrome_graph(self, counts, method: str = METHOD_SPITZ):
        """Generate weighted syndrome graph from result counts.

        Args:
            counts (dict): Counts dictionary of the results used to calculate
            the weights.
            method (string): Method to used for calculation. Supported
            methods are 'spitz' (default) and 'naive'.

        Additional information:
            Uses `counts` to estimate the probability of the errors that
            create the pairs of nodes in graph. The edge weights are then
            replaced with the corresponding -log(p/(1-p).
        """

        error_probs = self.get_error_probs(counts, method=method)

        boundary_nodes = []
        for n, node in enumerate(self.graph.nodes()):
            if node.is_boundary:
                boundary_nodes.append(n)

        for edge in self.graph.edge_list():
            if edge not in error_probs:
                # these are associated with the boundary, and appear as loops in error_probs
                # they need to be converted to this standard form
                bulk_n = list(set(edge).difference(set(boundary_nodes)))[0]
                p = error_probs[bulk_n, bulk_n]
            else:
                p = error_probs[edge]
            if 0 < p < 1:
                w = -np.log(p / (1 - p))
            elif p <= 0:  # negative values are assumed 0
                w = np.inf
            elif p == 1:
                w = -np.inf
            else:  # nan values are assumed maximally random
                w = 0
            edge_data = self.graph.get_edge_data(edge[0], edge[1])
            edge_data["weight"] = w
            self.graph.update_edge(edge[0], edge[1], edge_data)

    def make_error_graph(self, string: str):
        """Returns error graph.

        Args:
            string (str): A string describing the output from the code.

        Returns:
            The subgraph of graph which corresponds to the non-trivial
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

    @staticmethod
    def _layer_types(blocks: int, round_schedule: str, basis: str) -> List[str]:
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
            for index, supp in enumerate(all_z):
                node = DecodingGraphNode(time=time, qubits=supp, index=index)
                node.properties["highlighted"] = True
                graph.add_node(node)
                logging.debug("node %d t=%d %s", idx, time, supp)
                idxmap[(time, tuple(supp))] = idx
                node_layer.append(idx)
                idx += 1
            for index, supp in enumerate(boundary):
                # Add optional is_boundary property for pymatching
                node = DecodingGraphNode(is_boundary=True, qubits=supp, index=index)
                node.properties["highlighted"] = False
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
                        edge = DecodingGraphEdge(qubits=[com[0]], weight=1)
                        edge.properties["highlighted"] = False
                        edge.properties["measurement_error"] = 0
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
                edge = DecodingGraphEdge(qubits=[], weight=0)
                edge.properties["highlighted"] = False
                edge.properties["measurement_error"] = 0
                graph.add_edge(idxmap[(time, tuple(bound_g))], idxmap[(time, tuple(bound_h))], edge)
                logging.debug("spacelike boundary t=%d (%s, %s)", time, bound_g, bound_h)

            # Add (space)time-like edges from t to t-1
            # By construction, the qubit sets of pairs of vertices at graph and T
            # at times t-1 and t respectively
            # either (a) contain each other (graph subset T or T subset graph) and
            # |graph|,|T|>1,
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
                                edge = DecodingGraphEdge(qubits=[], weight=1)
                                edge.properties["highlighted"] = False
                                edge.properties["measurement_error"] = 1
                                graph.add_edge(
                                    idxmap[(time - 1, tuple(op_h))],
                                    idxmap[(time, tuple(op_g))],
                                    edge,
                                )
                                logging.debug("timelike t=%d (%s, %s)", time, op_g, op_h)
                            else:  # Case (b)
                                edge = DecodingGraphEdge(qubits=[com[0]], weight=1)
                                edge.properties["highlighted"] = False
                                edge.properties["measurement_error"] = 1
                                graph.add_edge(
                                    idxmap[(time - 1, tuple(op_h))],
                                    idxmap[(time, tuple(op_g))],
                                    edge,
                                )
                                logging.debug("spacetime hook t=%d (%s, %s)", time, op_g, op_h)
                                logging.debug(" qubits %s", [com[0]])
                # Add a single time-like edge between boundary vertices at
                # time t-1 and t
                edge = DecodingGraphEdge(qubits=[], weight=0)
                edge.properties["highlighted"] = False
                edge.properties["measurement_error"] = 0
                graph.add_edge(
                    idxmap[(time - 1, tuple(boundary[0]))], idxmap[(time, tuple(boundary[0]))], edge
                )
                logging.debug("boundarylink t=%d", time)

        self.idxmap = idxmap
        self.node_layers = node_layers
        self.graph = graph
