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
import copy
from typing import List, Tuple, Union

import numpy as np
import rustworkx as rx
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.exceptions import QiskitQECError
from qiskit_qec.utils import DecodingGraphEdge, DecodingGraphNode


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

    def __init__(self, code, brute=False, graph=None, hyperedges=None):
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
            self.hyperedges = hyperedges
        else:
            self._make_syndrome_graph()

        self.logical_nodes = []
        for node in self.graph.nodes():
            if node.is_logical:
                self.logical_nodes.append(node)

        self.update_attributes()

    def update_attributes(self):
        """
        Calculates properties of the graph used by `node_index` and `edge_in_graph`.
        If `graph` is updated this method should called to update these properties.
        """
        self._edge_set = set(self.graph.edge_list())
        self._node_index = {}
        for n, node in enumerate(self.graph.nodes()):
            clean_node = copy.deepcopy(node)
            clean_node.properties = {}
            self._node_index[clean_node] = n

    def node_index(self, node):
        """
        Given a node of `graph`, returns the corrsponding index.

        Args:
            node (DecodingGraphNode): Node of the graph.

        Returns:
            n (int): Index corresponding to the node within the graph.
        """
        clean_node = copy.deepcopy(node)
        clean_node.properties = {}
        return self._node_index[clean_node]

    def edge_in_graph(self, edge):
        """
        Given a pair of node indices for `graph`, determines whether
        the edge exists within the graph.

        Args:
            edge (tuple): Pair of node indices for the graph.

        Returns:
            in_graph (bool): Whether the edge is within the graph.
        """
        return edge in self._edge_set

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
                                if not (source.is_logical and target.is_logical):
                                    qubits = list(set(source.qubits).intersection(target.qubits))
                                    if not qubits:
                                        continue
                                if (
                                    source.time != target.time
                                    and len(qubits) > 1
                                    and not source.is_logical
                                    and not target.is_logical
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
        self, counts, logical: str = "0", method: str = METHOD_SPITZ, return_samples=False
    ) -> List[Tuple[Tuple[int, int], float]]:
        """
        Generate probabilities of single error events from result counts.

        Args:
            counts (dict): Counts dictionary of the results to be analyzed.
            logical (string): Logical value whose results are used.
            method (string): Method to used for calculation. Supported
            methods are 'spitz' (default) and 'naive'.
            return_samples (bool): Whether to also return the number of
            samples used to calculated each probability.
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
                error_nodes = set(self.code.string2nodes(string, logical=logical))

                for node0 in error_nodes:
                    n0 = self.node_index(node0)
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
                if self.graph[n0].is_logical:
                    boundary.append(n1)
                elif self.graph[n1].is_logical:
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
                error_nodes = set(self.code.string2nodes(string, logical=logical))
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
            samples = {}
            for n0, n1 in self.graph.edge_list():
                if count[n0, n1]["00"] > 0:
                    ratio = count[n0, n1]["11"] / count[n0, n1]["00"]
                else:
                    ratio = np.nan
                p = ratio / (1 + ratio)
                if self.graph[n0].is_logical and not self.graph[n1].is_logical:
                    edge = (n1, n1)
                elif not self.graph[n0].is_logical and self.graph[n1].is_logical:
                    edge = (n0, n0)
                else:
                    edge = (n0, n1)
                error_probs[edge] = p
                samples[edge] = count[n0, n1]["11"] + count[n0, n1]["00"]

        if return_samples:
            return error_probs, samples
        else:
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
            if node.is_logical:
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
            edge_data.weight = w
            self.graph.update_edge(edge[0], edge[1], edge_data)

    def make_error_graph(self, data: Union[str, List], all_logicals=True):
        """Returns error graph.

        Args:
            data: Either an ouput string of the code, or a list of
            nodes for the code.
            all_logicals(bool): Whether to do all logicals

        Returns:
            The subgraph of graph which corresponds to the non-trivial
            syndrome elements in the given string.
        """

        E = rx.PyGraph(multigraph=False)
        if isinstance(data, str):
            nodes = self.code.string2nodes(data, all_logicals=all_logicals)
        else:
            if all_logicals:
                nodes = list(set(data).union(set(self.logical_nodes)))
            else:
                nodes = data
        for node in nodes:
            if node not in E.nodes():
                E.add_node(node)

        # for each pair of nodes in error create an edge and weight with the
        # distance
        def weight_fn(edge):
            return int(edge.weight)

        distance_matrix = rx.graph_floyd_warshall_numpy(self.graph, weight_fn=weight_fn)

        for source_index in E.node_indexes():
            for target_index in E.node_indexes():
                source = E[source_index]
                target = E[target_index]
                if target != source:
                    ns = self.node_index(source)
                    nt = self.node_index(target)
                    distance = distance_matrix[ns][nt]
                    if np.isfinite(distance):
                        qubits = list(set(source.qubits).intersection(target.qubits))
                        distance = int(distance)
                        E.add_edge(source_index, target_index, DecodingGraphEdge(qubits, distance))
        return E

    def clean_measurements(self, nodes: List):
        """
        Removes pairs of nodes that obviously correspond to measurement errors
        from a list of nodes.

        Args:
            nodes: A list of nodes.
        Returns:
            nodes: The input list of nodes, with pairs removed if they obviously
            correspond to a measurement error.

        """

        # order the nodes by where and when
        node_pos = {}
        for node in nodes:
            if not node.is_boundary:
                if node.index not in node_pos:
                    node_pos[node.index] = {}
                node_pos[node.index][node.time] = self.node_index(node)
        # find pairs corresponding to time-like edges
        all_pairs = set()
        for node_times in node_pos.values():
            ts = list(node_times.keys())
            ts.sort()
            for j in range(len(ts) - 1):
                if ts[j + 1] - ts[j] <= 2:
                    n0 = node_times[ts[j]]
                    n1 = node_times[ts[j + 1]]
                    if self.edge_in_graph((n0, n1)) or self.edge_in_graph((n1, n0)):
                        all_pairs.add((n0, n1))
        # filter out those that share nodes
        all_nodes = set()
        common_nodes = set()
        for pair in all_pairs:
            for n in pair:
                if n in all_nodes:
                    common_nodes.add(n)
                all_nodes.add(n)
        paired_ns = set()
        for pair in all_pairs:
            if pair[0] not in common_nodes:
                if pair[1] not in common_nodes:
                    for n in pair:
                        paired_ns.add(n)
        # return the nodes that were not paired
        ns = set(self.node_index(node) for node in nodes)
        unpaired_ns = ns.difference(paired_ns)
        return [self.graph.nodes()[n] for n in unpaired_ns]

    def get_edge_graph(self):
        """
        Returns a copy of the graph that uses edges to store information
        about the effects of errors on logical operators. This is done
        via the `'fault_ids'` of the edges. No logical nodes are present
        in such a graph.

        Returns:
            edge_graph (rx.PyGraph): The edge graph.
        """

        nodes = self.graph.nodes()
        # get a list of boundary nodes
        bns = []
        for n, node in enumerate(nodes):
            if node.is_logical:
                bns.append(n)
        # find pairs of bulk edges that have overlap with a boundary
        bedge = {}
        # and their edges connecting to the boundary, that we'll discard
        spares = set()
        for edge, (n0, n1) in zip(self.graph.edges(), self.graph.edge_list()):
            if not nodes[n0].is_logical and not nodes[n1].is_logical:
                for n2 in bns:
                    adj = set(edge.qubits).intersection(set(nodes[n2].qubits))
                    if adj:
                        if (n0, n1) not in bedge:
                            bedge[n0, n1] = {nodes[n2].index}
                        else:
                            bedge[n0, n1].add(nodes[n2].index)
                        for n in (n0, n1):
                            spares.add((n, n2))
                            spares.add((n2, n))
        # find bulk-boundary pairs not covered by the above
        for (n0, n1) in self.graph.edge_list():
            n2 = None
            for n in (n0, n1):
                if nodes[n].is_logical:
                    n2 = n
            if n2 is not None:
                if (n0, n1) not in spares:
                    adj = set(nodes[n2].qubits)
                    for n in (n0, n1):
                        adj = adj.intersection(set(nodes[n].qubits))
                    if (n0, n1) not in bedge:
                        bedge[n0, n1] = {nodes[n2].index}
                    else:
                        bedge[n0, n1].add(nodes[n2].index)
        # make a new graph with fault_ids on boundary edges, and ignoring the spare edges
        edge_graph = rx.PyGraph(multigraph=False)
        for node in nodes:
            edge_graph.add_node(copy.copy(node))
        for edge, (n0, n1) in zip(self.graph.edges(), self.graph.edge_list()):
            if (n0, n1) in bedge:
                edge.fault_ids = bedge[n0, n1]
                edge_graph.add_edge(n0, n1, edge)
            elif (n0, n1) not in spares and (n1, n0) not in spares:
                edge.fault_ids = set()
                edge_graph.add_edge(n0, n1, edge)
        # turn logical nodes into boundary nodes
        for node in edge_graph.nodes():
            if node.is_logical:
                node.is_boundary = True
                node.is_logical = False
        return edge_graph

    def get_node_graph(self):
        """
        Returns a copy of the graph that uses logical nodes to store information
        about the effects of errors on logical operators. No non-trivial `'fault_ids'`
        are present in such a graph.

        Returns:
            node_graph (rx.PyGraph): The node graph.
        """
        node_graph = self.graph.copy()
        for edge, (n0, n1) in zip(self.graph.edges(), self.graph.edge_list()):
            if edge.fault_ids:
                # is the edge has fault ids, make corresponding logical nodes
                # and connect them to these edges
                for index in edge.fault_ids:
                    node2 = DecodingGraphNode(is_logical=True, index=index)
                    n2 = node_graph.add_node(node2)
                    node_graph.add_edge(n0, n2, copy.copy(edge))
                    node_graph.add_edge(n1, n2, copy.copy(edge))
        for edge in self.graph.edges():
            edge.fault_ids = set()
        return node_graph


def make_syndrome_graph_from_aer(code, shots=1):
    """
    Generates a graph and list of hyperedges for a given code by inserting Pauli errors
    around the gates of the base circuit for that code. Also supplied information regarding which
    edges where generated by which Pauli insertions.

    Args:
        code (CodeCircuit): Code for which the graph is to be made
        shots (int): Number of shots used in simulations.
    Returns:
        graph (PyGraph): Graph of the form used in `DecodingGraph`.
        hyperedges (list): List of hyperedges of the form used in `DecodingGraph`.
        hyperedge_errors (list): A list of the Pauli insertions that causes each hyperedge.
        These are specified by a tuple that specifies the type of Pauli and where it was
        inserted: (circuit depth, qreg index, qubit index, Pauli type).
        error_circuit (dict): Keys are the above tuples, and values are the corresponding
        circuits.
    """

    graph = rx.PyGraph(multigraph=False)
    hyperedges = []
    hyperedge_errors = []

    qc = code.circuit[code.base]
    blank_qc = QuantumCircuit()
    for qreg in qc.qregs:
        blank_qc.add_register(qreg)
    for creg in qc.cregs:
        blank_qc.add_register(creg)

    error_circuit = {}
    depth = len(qc)
    for j in range(depth):
        gate = qc.data[j][0].name
        qubits = qc.data[j][1]
        if gate not in ["measure", "reset", "barrier"] and len(qubits) != 2:
            for error in ["x", "y", "z"]:
                for qubit in qubits:
                    temp_qc = copy.deepcopy(blank_qc)
                    for qreg in qc.qregs:
                        if qubit in qreg:
                            break
                    temp_qc_name = (j, qc.qregs.index(qreg), qreg.index(qubit), error)
                    temp_qc.data = qc.data[0:j]
                    getattr(temp_qc, error)(qubit)
                    temp_qc.data += qc.data[j : depth + 1]
                    error_circuit[temp_qc_name] = temp_qc
        elif len(qubits) == 2:
            qregs = []
            for qubit in qubits:
                for qreg in qc.qregs:
                    if qubit in qreg:
                        qregs.append(qreg)
                        break
            for pauli_0 in ["id", "x", "y", "z"]:
                for pauli_1 in ["id", "x", "y", "z"]:
                    if not pauli_0 == pauli_1 == "id":
                        temp_qc = copy.deepcopy(blank_qc)
                        temp_qc_name = (
                            j,
                            (qc.qregs.index(qregs[0]), qc.qregs.index(qregs[1])),
                            (qregs[0].index(qubits[0]), qregs[1].index(qubits[1])),
                            pauli_0 + "," + pauli_1,
                        )
                        temp_qc.data = qc.data[0:j]
                        getattr(temp_qc, pauli_0)(qubits[0])
                        getattr(temp_qc, pauli_1)(qubits[1])
                        temp_qc.data += qc.data[j : depth + 1]
                        error_circuit[temp_qc_name] = temp_qc
        elif gate == "measure":
            pre_error = "x"
            for post_error in ["id", "x"]:
                for qubit in qubits:
                    temp_qc = copy.deepcopy(blank_qc)
                    for qreg in qc.qregs:
                        if qubit in qreg:
                            break
                    temp_qc_name = (
                        j,
                        qc.qregs.index(qreg),
                        qreg.index(qubit),
                        pre_error + "_" + post_error,
                    )
                    temp_qc.data = qc.data[0:j]
                    getattr(temp_qc, pre_error)(qubit)
                    temp_qc.data.append(qc.data[j])
                    getattr(temp_qc, post_error)(qubit)
                    temp_qc.data += qc.data[j + 1 : depth + 1]
                    error_circuit[temp_qc_name] = temp_qc

    errors = []
    circuits = []
    for error, circuit in error_circuit.items():
        errors.append(error)
        circuits.append(circuit)

    result = AerSimulator().run(circuits, shots=shots).result()
    no_nodes = []
    for j, circuit in enumerate(circuits):
        for string in result.get_counts(j):
            nodes = code.string2nodes(string)
            for node in nodes:
                if node not in graph.nodes():
                    graph.add_node(node)
            hyperedge = {}
            for source in nodes:
                for target in nodes:
                    if target != source or (len(nodes) == 1):
                        n0 = graph.nodes().index(source)
                        n1 = graph.nodes().index(target)
                        if not (source.is_logical and target.is_logical):
                            qubits = list(set(source.qubits).intersection(target.qubits))
                        if source.time != target.time and len(qubits) > 1:
                            qubits = []
                        edge = DecodingGraphEdge(qubits, 1)
                        graph.add_edge(n0, n1, edge)
                        if (n1, n0) not in hyperedge:
                            hyperedge[n0, n1] = edge
            if hyperedge:
                if hyperedge not in hyperedges:
                    hyperedges.append(hyperedge)
                    hyperedge_errors.append([])
                k = hyperedges.index(hyperedge)
                if errors[j] not in hyperedge_errors[k]:
                    hyperedge_errors[k].append(errors[j])
            else:
                if errors[j] not in no_nodes:
                    no_nodes.append(errors[j])
    hyperedges.append({})
    hyperedge_errors.append(no_nodes)

    return graph, hyperedges, hyperedge_errors, error_circuit
