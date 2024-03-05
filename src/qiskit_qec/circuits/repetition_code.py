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

"""Generates circuits based on repetition codes."""
from copy import deepcopy
from typing import List, Optional, Tuple

import numpy as np
import rustworkx as rx
import networkx as nx

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, transpile
from qiskit.circuit.library import XGate, RZGate
from qiskit.transpiler import PassManager, InstructionDurations
from qiskit_ibm_provider.transpiler.passes.scheduling import DynamicCircuitInstructionDurations
from qiskit_ibm_provider.transpiler.passes.scheduling import PadDynamicalDecoupling
from qiskit_ibm_provider.transpiler.passes.scheduling import ALAPScheduleAnalysis

from qiskit_qec.circuits.code_circuit import CodeCircuit
from qiskit_qec.utils import DecodingGraphEdge, DecodingGraphNode
from qiskit_qec.utils.decoding_graph_attributes import _nodes2cpp
from qiskit_qec.circuits._c_circuits import _c_check_nodes, _c_is_cluster_neutral


def _separate_string(string):
    separated_string = []
    for syndrome_type_string in string.split("  "):
        separated_string.append(syndrome_type_string.split(" "))
    return separated_string


class RepetitionCodeCircuit(CodeCircuit):
    """RepetitionCodeCircuit class."""

    def __init__(
        self,
        d: int,
        T: int,
        xbasis: bool = False,
        resets: bool = False,
        delay: Optional[int] = None,
        barriers: bool = False,
    ):
        """
        Creates the circuits corresponding to a logical 0 and 1 encoded
        using a repetition code.

        Implementation of a distance d repetition code, implemented over
        T syndrome measurement rounds.

        Args:
            d (int): Number of code qubits (and hence repetitions) used.
            T (int): Number of rounds of ancilla-assisted syndrome measurement.
            xbasis (bool): Whether to use the X basis to use for encoding (Z basis used by default).
            resets (bool): Whether to include a reset gate after mid-circuit measurements.
            delay (float): Time (in dt) to delay after mid-circuit measurements (and delay).
            barriers (bool): Whether to include barriers between different sections of the code.


        Additional information:
            No measurements are added to the circuit if `T=0`. Otherwise
            `T` rounds are added, followed by measurement of the code
            qubits (corresponding to a logical measurement and final
            syndrome measurement round).
        """

        super().__init__()

        self.n = d
        self.d = d
        self.T = 0

        self.code_qubit = QuantumRegister(d, "code_qubit")
        self.link_qubit = QuantumRegister((d - 1), "link_qubit")
        self.qubit_registers = {"code_qubit", "link_qubit"}

        self.link_bits = []
        self.code_bit = ClassicalRegister(d, "code_bit")

        self.circuit = {}
        for log in ["0", "1"]:
            self.circuit[log] = QuantumCircuit(self.link_qubit, self.code_qubit, name=log)

        self._xbasis = xbasis
        self._resets = resets
        self._barriers = barriers

        self._preparation()

        delay = delay or 0

        for _ in range(T - 1):
            self.syndrome_measurement(delay=delay)

        if T != 0:
            self.syndrome_measurement(final=True)
            self.readout()

        self.gauge_ops = [[j, j + 1] for j in range(self.d - 1)]
        self.measured_logical = [[0], [self.d - 1]]
        self.basis = "x"

        self.resets = resets
        self.delay = delay
        self.base = "0"

    def get_circuit_list(self) -> List[QuantumCircuit]:
        """Returns circuit list.

        circuit_list: self.circuit as a list, with
        circuit_list[0] = circuit['0']
        circuit_list[1] = circuit['1']
        """
        circuit_list = [self.circuit[log] for log in ["0", "1"]]
        return circuit_list

    def x(self, logs=("0", "1"), barrier=False):
        """Applies a logical x to the circuits for the given logical values.

        Args:
            logs (list or tuple): List or tuple of logical values expressed as
                strings.
            barrier (bool): Boolean denoting whether to include a barrier at
                the start.
        """
        barrier = barrier or self._barriers
        for log in logs:
            if barrier and (log == "1" or self._xbasis):
                self.circuit[log].barrier()
            if self._xbasis:
                self.circuit[log].z(self.code_qubit)
            else:
                self.circuit[log].x(self.code_qubit)

    def _preparation(self):
        """Prepares logical bit states by applying an x to the circuit that will
        encode a 1.
        """
        for log in ["0", "1"]:
            if self._xbasis:
                self.circuit[log].h(self.code_qubit)
        self.x(["1"])

    def syndrome_measurement(self, final: bool = False, barrier: bool = False, delay: int = 0):
        """Application of a syndrome measurement round.

        Args:
            final (bool): Whether to disregard the reset (if applicable) due to this
            being the final syndrome measurement round.
            barrier (bool): Boolean denoting whether to include a barrier at the start.
            delay (float): Time (in dt) to delay after mid-circuit measurements (and delay).
        """
        barrier = barrier or self._barriers

        self.link_bits.append(ClassicalRegister((self.d - 1), "round_" + str(self.T) + "_link_bit"))

        for log in ["0", "1"]:
            self.circuit[log].add_register(self.link_bits[-1])

            # entangling gates
            if barrier:
                self.circuit[log].barrier()
            if self._xbasis:
                self.circuit[log].h(self.link_qubit)
            for j in range(self.d - 1):
                if self._xbasis:
                    self.circuit[log].cx(self.link_qubit[j], self.code_qubit[j])
                else:
                    self.circuit[log].cx(self.code_qubit[j], self.link_qubit[j])
            for j in range(self.d - 1):
                if self._xbasis:
                    self.circuit[log].cx(self.link_qubit[j], self.code_qubit[j + 1])
                else:
                    self.circuit[log].cx(self.code_qubit[j + 1], self.link_qubit[j])
            if self._xbasis:
                self.circuit[log].h(self.link_qubit)

            # measurement
            if barrier:
                self.circuit[log].barrier()
            for j in range(self.d - 1):
                self.circuit[log].measure(self.link_qubit[j], self.link_bits[self.T][j])

            # resets
            if self._resets and not final:
                if barrier:
                    self.circuit[log].barrier()
                for j in range(self.d - 1):
                    self.circuit[log].reset(self.link_qubit[j])

            # delay
            if delay > 0 and not final:
                if barrier:
                    self.circuit[log].barrier()
                for j in range(self.d - 1):
                    self.circuit[log].delay(delay, self.link_qubit[j])

        self.T += 1

    def readout(self):
        """
        Readout of all code qubits, which corresponds to a logical measurement
        as well as allowing for a measurement of the syndrome to be inferred.
        """
        for log in ["0", "1"]:
            if self._xbasis:
                self.circuit[log].h(self.code_qubit)
            self.circuit[log].add_register(self.code_bit)
            self.circuit[log].measure(self.code_qubit, self.code_bit)

    def _process_string(self, string):
        # logical readout taken from
        measured_log = string[0] + " " + string[self.d - 1]

        if self._resets:
            syndrome = string[self.d :]
        else:
            # if there are no resets, results are cumulative and need to be separated
            cumsyn_list = string[self.d :].split(" ")
            syndrome_list = []
            for tt, cum_syn in enumerate(cumsyn_list[0:-1]):
                syn = ""
                for j in range(len(cum_syn)):
                    syn += str(int(cumsyn_list[tt][j] != cumsyn_list[tt + 1][j]))
                syndrome_list.append(syn)
            syndrome_list.append(cumsyn_list[-1])
            syndrome = " ".join(syndrome_list)

        # final syndrome deduced from final code qubit readout
        full_syndrome = ""
        for j in range(self.d - 1):
            full_syndrome += "0" * (string[j] == string[j + 1]) + "1" * (string[j] != string[j + 1])
        # results from all other syndrome measurements then added
        full_syndrome = full_syndrome + syndrome

        # changes between one syndrome and the next then calculated
        syndrome_list = full_syndrome.split(" ")
        syndrome_changes = ""
        for t in range(self.T + 1):
            for j in range(self.d - 1):
                if t == 0:
                    change = syndrome_list[-1][j] != "0"
                else:
                    change = syndrome_list[-t - 1][j] != syndrome_list[-t][j]
                syndrome_changes += "0" * (not change) + "1" * change
            syndrome_changes += " "

        # the space separated string of syndrome changes then gets a
        # double space separated logical value on the end
        new_string = measured_log + "  " + syndrome_changes[:-1]

        return new_string

    def string2nodes(self, string, **kwargs):
        """Convert output string from circuits into a set of nodes.

        Args:
            string (string): Results string to convert.
            kwargs (dict): Additional keyword arguments.
                    logical (str): Logical value whose results are used ('0' as default).
                    all_logicals (bool): Whether to include logical nodes
                    irrespective of value. (False as default).

        Returns:
            dict: List of nodes corresponding to to the non-trivial
            elements in the string.

        Additional information:
        Strings are read right to left, but lists*
        are read left to right. So, we have some ugly indexing
        code whenever we're dealing with both strings and lists.
        """

        # set kwargs
        all_logicals = kwargs.get("all_logicals")
        logical = kwargs.get("logical")
        if logical is None:
            logical = "0"

        string = self._process_string(string)
        # [ <logical>, <syn>, <syn>,...]
        separated_string = _separate_string(string)
        nodes = []

        # logical/boundary nodes
        boundary = separated_string[0]  # [<last_elem>, <init_elem>]
        for bqec_index, belement in enumerate(boundary[::-1]):
            if all_logicals or belement != logical:
                bqubits = self.measured_logical[bqec_index]
                bnode = DecodingGraphNode(
                    is_logical=True, is_boundary=True, qubits=bqubits, index=bqec_index
                )
                nodes.append(bnode)

        # bulk nodes
        for syn_type in range(1, len(separated_string)):
            for syn_round in range(len(separated_string[syn_type])):
                elements = separated_string[syn_type][syn_round]
                for qec_index, element in enumerate(elements[::-1]):
                    if element == "1":
                        qubits = self.gauge_ops[qec_index]
                        node = DecodingGraphNode(time=syn_round, qubits=qubits, index=qec_index)
                        nodes.append(node)
        return nodes

    def string2raw_logicals(self, string):
        """
        Extracts raw logicals from output string.
        Args:
            string (string): Results string from which to extract logicals
        Returns:
            list: Raw values for logical operators that correspond to nodes.
        """
        return [string.split(" ", maxsplit=1)[0][-1]]

    def check_nodes(self, nodes, ignore_extras=False, minimal=False):
        """
        Determines whether a given set of nodes are neutral. If so, also
        determines any additional logical readout qubits that would be
        flipped by the errors creating such a cluster and how many errors
        would be required to make the cluster.
        Args:
            nodes (list): List of nodes, of the type produced by `string2nodes`.
            ignore_extras (bool): If `True`, undeeded boundary nodes are
            ignored.
            minimal (bool): Whether output should only reflect the minimal error
            case.
        Returns:
            neutral (bool): Whether the nodes independently correspond to a valid
            set of errors.
            flipped_logical_nodes (list): List of qubits nodes for logical
            operators that are flipped by the errors, that were not included
            in the original nodes.
            num_errors (int): Minimum number of errors required to create nodes.
        """

        # see which qubits for logical zs are given and collect bulk nodes
        given_logicals = []
        for node in nodes:
            if node.is_logical:
                given_logicals += node.qubits
        given_logicals = set(given_logicals)

        # bicolour code qubits according to the domain walls
        walls = []
        for node in nodes:
            if not node.is_logical:
                walls.append(node.qubits[1])
        walls.sort()
        c = 0
        colors = ""
        for j in range(self.d):
            if walls:
                if walls[0] == j:
                    c = (c + 1) % 2
                    walls.remove(j)
            colors += str(c)
        colors = colors[::-1]

        # determine which were in the minority
        error_c_min = str(int(colors.count("1") < self.d / 2))
        # and majority
        error_c_max = str((int(error_c_min) + 1) % 2)

        # list the colours with the max error one first
        # (unless we do min only)
        error_cs = []
        if not minimal:
            error_cs.append(error_c_max)
        error_cs.append(error_c_min)

        # see what happens for both colours
        # if neutral for maximal, it's neutral
        # otherwise, it is whatever it is for the minimal
        for error_c in error_cs:
            num_errors = colors.count(error_c)

            # determine the corresponding flipped logicals
            flipped_logicals = []
            for j in [0, self.d - 1]:
                if colors[-1 - j] == error_c:
                    flipped_logicals.append(j)
            flipped_logicals = set(flipped_logicals)

            # if unneeded logical zs are given, cluster is not neutral
            # (unless this is ignored)
            if (not ignore_extras) and given_logicals.difference(flipped_logicals):
                neutral = False
            # otherwise, report only needed logicals that aren't given
            else:
                neutral = True
                flipped_logicals = flipped_logicals.difference(given_logicals)

            flipped_logical_nodes = []
            for flipped_logical in flipped_logicals:
                qubits = [flipped_logical]
                elem = self.measured_logical.index(qubits)
                node = DecodingGraphNode(
                    is_logical=True, is_boundary=True, qubits=qubits, index=elem
                )
                flipped_logical_nodes.append(node)

            if neutral and not flipped_logical_nodes:
                break

        return neutral, flipped_logical_nodes, num_errors

    def is_cluster_neutral(self, nodes):
        """
        Determines whether or not the cluster is neutral, meaning that one or more
        errors could have caused the set of nodes (syndrome changes) passed
        to the method.
        Args:
            nodes (list of nodes)
        """
        return not bool(len(nodes) % 2)

    def partition_outcomes(
        self, round_schedule: str, outcome: List[int]
    ) -> Tuple[List[List[int]], List[List[int]], List[int]]:
        """Extract measurement outcomes."""
        # split into gauge and final outcomes
        outcome = "".join([str(c) for c in outcome])
        outcome = outcome.split(" ")
        gs = outcome[0:-1]
        gauge_outcomes = [[int(c) for c in r] for r in gs]
        finals = outcome[-1]
        # if circuit did not use resets, construct standard output
        if not self.resets:
            for i, layer in enumerate(gauge_outcomes):
                for j, gauge_op in enumerate(layer):
                    if i > 0:
                        gauge_outcomes[i][j] = (gauge_op + gauge_outcomes[i - 1][j]) % 2
        # assign outcomes to the correct gauge ops
        if round_schedule == "z":
            x_gauge_outcomes = []
            z_gauge_outcomes = gauge_outcomes
        else:
            x_gauge_outcomes = gauge_outcomes
            z_gauge_outcomes = []
        final_outcomes = [int(c) for c in finals]

        return x_gauge_outcomes, z_gauge_outcomes, final_outcomes


def add_edge(graph, pair, edge=None):
    """
    Adds an edge correspoding to the given pair of nodes to the given graph,
    adding also the nodes themselves if not already present.
    """
    ns = []
    for node in pair:
        if node not in graph.nodes():
            ns.append(graph.add_node(node))
        else:
            ns.append(graph.nodes().index(node))
    graph.add_edge(ns[0], ns[1], edge)
    return ns


class ArcCircuit(CodeCircuit):
    """Anisotropic repetition code class."""

    METHOD_SPITZ: str = "spitz"
    METHOD_NAIVE: str = "naive"
    AVAILABLE_METHODS = {METHOD_SPITZ, METHOD_NAIVE}

    def __init__(
        self,
        links: list,
        T: int,
        basis: str = "xy",
        logical: str = "0",
        resets: bool = True,
        delay: Optional[int] = None,
        barriers: bool = True,
        color: Optional[dict] = None,
        max_dist: int = 2,
        schedule: Optional[list] = None,
        run_202: bool = True,
        rounds_per_202: int = 9,
        conditional_reset: bool = False,
    ):
        """
        Creates circuits corresponding to an anisotropic repetition code implemented over T syndrome
        measurement rounds, with the syndrome measurements provided.
        Args:
            links (list): List of tuples (c0, a, c1), where c0 and c1 are the two code qubits in each
            syndrome measurement, and a is the auxiliary qubit used.
            T (int): Number of rounds of syndrome measurement.
            basis (string): Pair of `'x'`, `'y'` and `'z'`, specifying the pair of local bases to be
            used.
            logical (string): Logical value to store (`'0'` or `'1'`).
            resets (bool): Whether to include a reset gate after mid-circuit measurements.
            ff (bool): Whether to correct the effects of [[2,0,2]] sequences via feed forward.
            delay (float): Time (in dt) to delay after mid-circuit measurements (and delay).
            barriers (bool): Whether to include barriers between different sections of the code.
            color (dict): Dictionary with code qubits as keys and 0 or 1 for each value, to specify
            a predetermined bicoloring. If not provided, a bicoloring is found on initialization.
            max_dist (int): Maximum edge distance used when determining the bicoloring of code qubits.
            schedule (list): Specifies order in which entangling gates are applied in each syndrome
            measurement round. Each element is a list of lists [c, a] for entangling gates to be
            applied simultaneously.
            run_202 (bool): Whether to run [[2,0,2]] sequences. This will be overwritten if T is not high
            enough (at least rounds_per_202xlen(links)).
            rounds_per_202 (int): Number of rounds that are part of the 202, including the typical link
            measurements at the beginning and end. At least 9 are required to get an event dedicated to
            conjugate errors.
            conditional_reset: Whether to apply conditional resets (an x conditioned on the result of the
            previous measurement), rather than a reset gate.
        """

        super().__init__()

        self.links = links
        self.basis = basis
        self.logical = logical
        self._barriers = barriers
        self._max_dist = max_dist
        self.delay = delay or 0
        self.conditional_reset = conditional_reset

        # calculate coloring and schedule, etc
        if color is None:
            self._coloring()
        else:
            self.color = color
        if schedule is None:
            self._scheduling()
        else:
            self.schedule = schedule
        self._get_cycles()
        self._preparation()

        # determine the placement of [2,0,2] rounds
        self.rounds_per_202 = rounds_per_202
        if run_202:
            self.links_202 = []
            for link in self.links:
                logical_overlap = {link[0], link[2]}.intersection(set(self.z_logicals))
                if not logical_overlap:
                    self.links_202.append(link)
            num_links = len(self.links_202)
            if num_links > 0:
                self.rounds_per_link = int(np.floor(T / num_links))
                self.metabuffer = np.ceil((T - num_links * self.rounds_per_link) / 2)
                self.roundbuffer = np.ceil((self.rounds_per_link - self.rounds_per_202) / 2)
                if self.roundbuffer > 0:
                    self.roundbuffer -= 1
                self.run_202 = self.rounds_per_link >= self.rounds_per_202
            else:
                self.run_202 = False
        else:
            self.run_202 = False
        self.resets = resets or self.run_202
        if not self.run_202:
            self.rounds_per_link = np.inf

        # create the circuit
        self.base = basis
        self.T = 0
        for _ in range(T - 1):
            self._syndrome_measurement()
        if T != 0:
            self._syndrome_measurement(final=True)
        self._readout()

        self._cpp_link_graph, self._cpp_link_neighbors = self._links2cpp()

    def _get_link_graph(self, max_dist=1):
        graph = rx.PyGraph()
        for link in self.links:
            add_edge(graph, (link[0], link[2]), {"distance": 1, "link qubit": link[1]})
        distance = rx.distance_matrix(graph)
        edges = graph.edge_list()
        for n0, node0 in enumerate(graph.nodes()):
            for n1, node1 in enumerate(graph.nodes()):
                if n0 < n1:
                    if (n0, n1) not in edges:
                        dist = distance[n0, n1]
                        if dist < max_dist:
                            add_edge(graph, (node0, node1), {"distance": dist})
        return graph

    def _get_cycles(self):
        """
        For each edge in the link graph (expressed in terms of the pair of qubits), the
        set of qubits around adjacent cycles is found.
        """

        self.link_graph = self._get_link_graph()
        self.degree = {}
        for n, q in enumerate(self.link_graph.nodes()):
            self.degree[q] = self.link_graph.degree(n)
        degrees = list(self.degree.values())
        self._linear = degrees.count(1) == 2 and degrees.count(2) == len(degrees) - 2
        lg_edges = set(self.link_graph.edge_list())
        lg_nodes = self.link_graph.nodes()
        ng = nx.Graph()
        for n0, n1 in self.link_graph.edge_list():
            ng.add_edge(n0, n1)
        # express the cycles in terms of the ns of the link graph
        self.cycles = nx.minimum_cycle_basis(ng)
        # and for each pair of data qubits, list the cycles it is a part of
        self.cycle_dict = {(lg_nodes[edge[0]], lg_nodes[edge[1]]): set() for edge in lg_edges}
        for c, cycle in enumerate(self.cycles):
            for n0 in cycle:
                for n1 in cycle:
                    for edge in [(n0, n1), (n1, n0)]:
                        if edge in lg_edges:
                            self.cycle_dict[lg_nodes[edge[0]], lg_nodes[edge[1]]].add(c)

    def _coloring(self):
        """
        Creates a graph with a weight=1 edge for each link, and additional edges up to `max_weight`.
        Then performs a matching of edges in this graph. The code qubits of each pair are then
        bicolored. All unmatched code qubits are alternately bicolored.
        """

        graph = self._get_link_graph(self._max_dist)
        matching = rx.max_weight_matching(
            graph, max_cardinality=True, weight_fn=lambda edge: -int(edge["distance"])
        )
        self.color = {}
        unmatched = list(graph.node_indices())
        nodes = graph.nodes()
        for pair in matching:
            for j, n in enumerate(pair):
                self.color[nodes[n]] = j
                unmatched.remove(n)
        for j, n in enumerate(unmatched):
            # color opposite to a colored neighbor
            neighbors = graph.neighbors(n)
            if neighbors:
                for nn in neighbors:
                    if nodes[nn] in self.color:
                        self.color[nodes[n]] = (self.color[nodes[nn]] + 1) % 2
            else:
                # otherwise color arbitrarily
                self.color[nodes[n]] = j % 2

    def _get_coupling_graph(self, aux=None):
        """
        Returns a graph for pairs of nodes on which entangling gates are applied in the code.
        """

        if aux is None:
            aux = [link[1] for link in self.links]

        graph = rx.PyGraph()
        for link in self.links:
            if link[1] in aux:
                add_edge(graph, (link[0], link[1]), {})
                add_edge(graph, (link[1], link[2]), {})
        # we use max degree of the nodes as the edge weight, to delay bottlenecks
        for e, (n0, n1) in enumerate(graph.edge_list()):
            graph.edges()[e]["weight"] = max(graph.degree(n0), graph.degree(n1))

        return graph

    def _scheduling(self):
        """
        Determines the order in which entangling gates should be applied in each round.
        """

        link_dict = {link[1]: link for link in self.links}
        aux = set(link_dict.keys())

        def weight_fn(edge):
            return -int(edge["weight"])

        schedule = []
        while aux:
            # construct coupling graph for as yet unpaired auxiliaries (i.e. link qubits)
            graph = self._get_coupling_graph(aux)

            # find a min weight matching, and then another that exlcudes the pairs from the first
            matching = [rx.max_weight_matching(graph, max_cardinality=True, weight_fn=weight_fn)]
            cut_graph = deepcopy(graph)
            for n0, n1 in matching[0]:
                cut_graph.remove_edge(n0, n1)
            matching.append(
                rx.max_weight_matching(cut_graph, max_cardinality=True, weight_fn=weight_fn)
            )

            # rewrite the matchings to use nodes instead of indices, and to always place
            # the auxilliary second
            nodes = [graph.nodes(), cut_graph.nodes()]
            for j in range(2):
                matching[j] = list(matching[j])
                for p, pair in enumerate(matching[j]):
                    node_pair = [None, None]
                    for n in pair:
                        node = nodes[j][n]
                        node_pair[node in aux] = node
                    matching[j][p] = node_pair

            # determine which links are covered by the conjuction of these matchings
            matched_aux = [set(), set()]
            for j in range(2):
                for pair in matching[j]:
                    matched_aux[j].add(pair[1])
            completed = matched_aux[0].intersection(matched_aux[1])

            # add these matched pairs to the schedule
            for j in range(2):
                schedule.append([pair for pair in matching[j] if pair[1] in completed])

            # update the list of auxilliaries for links yet to be paired
            aux = aux.difference(completed)

        self.schedule = schedule

    def _rotate(self, basis, c, regqubit, inverse):
        """
        Rotates the given qubit to (or from) the basis specified by the color and the pair of
        bases used for the code.
        """
        if not inverse:
            if basis[c] in ["x", "y"]:
                self.circuit[basis].h(regqubit)
            if basis[c] == "y":
                self.circuit[basis].s(regqubit)
        else:
            if basis[c] == "y":
                self.circuit[basis].sdg(regqubit)
            if basis[c] in ["x", "y"]:
                self.circuit[basis].h(regqubit)

    def _basis_change(self, basis, inverse=False):
        """
        Rotates all code qubits to (or from) the bases required for the code.
        """
        for qubit, q in self.code_index.items():
            c = self.color[qubit]
            self._rotate(basis, c, self.code_qubit[q], inverse)

    def _preparation(self):
        """
        Creates the circuits and their registers.
        """

        # get a list of all code qubits (qubits[0]) and link qubits (qubits[1])
        qubits = [[], []]
        for link in self.links:
            for code_qubit in [link[0], link[2]]:
                if code_qubit not in qubits[0]:
                    qubits[0].append(code_qubit)
            qubits[1].append(link[1])
        self.qubits = [qubits[j] for j in range(2)]
        self.num_qubits = [len(qubits[j]) for j in range(2)]
        self.d = self.num_qubits[0]

        # define the quantum egisters
        self.code_qubit = QuantumRegister(self.num_qubits[0], "code_qubit")
        self.link_qubit = QuantumRegister(self.num_qubits[1], "link_qubit")
        self.qubit_registers = {"code_qubit", "link_qubit"}

        # for looking up where each qubit lives on the quantum registers
        self.code_index = {qubit: q for q, qubit in enumerate(list(qubits[0]))}
        self.link_index = {qubit: q for q, qubit in enumerate(list(qubits[1]))}

        # set up the classical registers
        self.link_bits = []
        self.code_bit = ClassicalRegister(len(qubits[0]), "code_bit")

        # create the circuits and initialize the code qubits
        self.circuit = {}
        for basis in list({self.basis, self.basis[::-1]}):
            self.circuit[basis] = QuantumCircuit(self.link_qubit, self.code_qubit, name=basis)
            if self.logical == "1":
                self.circuit[basis].x(self.code_qubit)
            self._basis_change(basis)

        # use degree 1 code qubits for logical z readouts (and boundary)
        graph = self._get_coupling_graph()
        self._leaves = False
        z_logicals = []
        self.boundary = []
        for n, node in enumerate(graph.nodes()):
            if graph.degree(n) == 1:
                z_logicals.append(node)
                self.boundary.append(node)
                self._leaves = True
        # if there are none, just use the first (not boundary)
        if not z_logicals:
            z_logicals = [min(self.code_index.keys())]
        self.z_logicals = z_logicals

    def _get_202(self, t):
        """
        Returns the position within a 202 sequence for the current measurement round:
        * `False` means not part of a 202 sequence;
        * Even taus use the standard coloring;
        * Odd taus use the flipped coloring.
        Also returns the link qubits for the link for which the 202 sequence is run and its neigbours.
        """
        # null values in case no 202 done during this round
        tau, qubit_l_202, qubit_l_nghbrs = None, None, [[], []]
        if self.run_202 and int(t / self.rounds_per_link) < len(self.links_202):
            # determine the link qubit for which the 202 sequence is run
            link = self.links_202[int(t / self.rounds_per_link)]
            # set the 202 link
            qubit_l_202 = link[1]
            #  determine where we are in the sequence
            tau = int(t % self.rounds_per_link - self.roundbuffer)
            if t < 0 or tau not in range(self.rounds_per_202):
                tau = False
            # determine the neighbouring link qubits that are suppressed
            graph = self._get_link_graph(0)
            nodes = graph.nodes()
            edges = graph.edge_list()
            ns = [nodes.index(link[j]) for j in [0, 2]]
            qubit_l_nghbrs = []
            for n in ns:
                neighbors = list(graph.incident_edges(n))
                neighbors.remove(list(edges).index(tuple(ns)))
                qubit_l_nghbrs.append(
                    [graph.get_edge_data_by_index(ngbhr)["link qubit"] for ngbhr in neighbors]
                )
        return tau, qubit_l_202, qubit_l_nghbrs

    def _syndrome_measurement(self, final: bool = False):
        """
        Adds a syndrome measurement round.
        Args:
            final (bool): Whether or not this is the final round.
        """

        self.link_bits.append(
            ClassicalRegister(self.num_qubits[1], "round_" + str(self.T) + "_link_bit")
        )

        tau, qubit_l_202, qubit_l_nghbrs = self._get_202(self.T)
        links_to_measure = set()
        links_to_reset = set()
        for basis, qc in self.circuit.items():
            if self._barriers:
                # pre round barrier
                qc.barrier()
            for pairs in self.schedule:
                for qubit_c, qubit_l in pairs:
                    q_c = self.code_index[qubit_c]
                    q_l = self.link_index[qubit_l]
                    neighbor = qubit_l in qubit_l_nghbrs[0] + qubit_l_nghbrs[1]
                    if not (tau in range(1, self.rounds_per_202 - 1) and neighbor):
                        c = self.color[qubit_c]
                        if qubit_l == qubit_l_202:
                            c = (c + tau) % 2
                        self._rotate(basis, c, self.code_qubit[q_c], True)
                        qc.cx(self.code_qubit[q_c], self.link_qubit[q_l])
                        self._rotate(basis, c, self.code_qubit[q_c], False)
                        links_to_measure.add(q_l)
                        if not (not isinstance(tau, bool) and tau == 0 and neighbor):
                            links_to_reset.add(q_l)

        # measurement and resets
        for basis, qc in self.circuit.items():
            # measurement
            if final:
                # already prep code qubits for readout
                self._basis_change(basis, inverse=True)
            if self._barriers:
                # post-round, pre-measurement barrier
                qc.barrier()
            qc.add_register(self.link_bits[self.T])
            for q_l in links_to_measure:
                qc.measure(self.link_qubit[q_l], self.link_bits[self.T][q_l])

            # resets
            if self.resets and not final:
                for q_l in links_to_reset:
                    if self.conditional_reset:
                        qc.x(self.link_qubit[q_l]).c_if(self.link_bits[self.T][q_l], 1)
                    else:
                        qc.reset(self.link_qubit[q_l])

            # correct
            if self.run_202:
                if tau == (self.rounds_per_202 - 1):
                    target_link = self.links[self.link_index[qubit_l_202]]
                    # for neighbouring links on both sides of the 202 link
                    for j in range(2):
                        if qubit_l_nghbrs[j]:
                            # get the first listed neighbouring link
                            control_link = self.links[self.link_index[qubit_l_nghbrs[j][0]]]
                            # find the qubit on which it overlaps with the 202
                            qubit_t = list(set(target_link).intersection(set(control_link)))[0]
                            # and the qubit whose result controls the feedforward
                            qubit_c = control_link[1]
                            # get their indices
                            q_t = self.code_index[qubit_t]
                            q_c = self.link_index[qubit_c]
                            # and the colour of the targeted qubit
                            c = self.color[qubit_t]
                            self._rotate(basis, c, self.code_qubit[q_t], True)
                            qc.x(self.code_qubit[q_t]).c_if(self.link_bits[self.T][q_c], 1)
                            self._rotate(basis, c, self.code_qubit[q_t], False)

            # delay
            if self.delay > 0 and not final:
                if self._barriers:
                    # post-reset, pre-delay barrier
                    qc.barrier()
                qc.delay(self.delay, self.link_qubit)

        self.T += 1

    def _readout(self):
        """
        Readout of all code qubits, which corresponds to a logical measurement
        as well as allowing for a measurement of the syndrome to be inferred.
        """

        for basis, qc in self.circuit.items():
            if self.T == 0:
                self._basis_change(basis, inverse=True)
                if self._barriers:
                    qc.barrier()
                # otherwise, code qubits are already prepped
            qc.add_register(self.code_bit)
            qc.measure(self.code_qubit, self.code_bit)

    def _process_string(self, string):
        # logical readout taken from assigned qubits
        measured_log = ""
        for qubit in self.z_logicals:
            j = self.code_index[qubit]
            measured_log += string[self.num_qubits[0] - j - 1] + " "

        if self.resets:
            syndrome = string[self.num_qubits[0] :]
        else:
            # if there are no resets, results are cumulative and need to be separated
            cumsyn_list = string[self.num_qubits[0] :].split(" ")
            syndrome_list = []
            for tt, cum_syn in enumerate(cumsyn_list[0:-1]):
                syn = ""
                for j in range(len(cum_syn)):
                    syn += str(int(cumsyn_list[tt][j] != cumsyn_list[tt + 1][j]))
                syndrome_list.append(syn)
            syndrome_list.append(cumsyn_list[-1])
            syndrome = " ".join(syndrome_list)

        # final syndrome deduced from final code qubit readout
        full_syndrome = ""
        for link in self.links:
            q = [self.num_qubits[0] - 1 - self.code_index[link[j]] for j in [0, -1]]
            full_syndrome += "0" * (string[q[0]] == string[q[1]]) + "1" * (
                string[q[0]] != string[q[1]]
            )
        full_syndrome = full_syndrome[::-1]
        # results from all other syndrome measurements then added
        full_syndrome = full_syndrome + syndrome

        # changes between appropriate results are then calculated
        syndrome_list = full_syndrome.split(" ")
        syndrome_changes = ""
        last_neighbors = []
        just_finished = False
        for t in range(self.T + 1):
            tau, qubit_l_202, qubit_l_nghbrs = self._get_202(t)
            all_neighbors = qubit_l_nghbrs[0] + qubit_l_nghbrs[1]
            for j in range(self.num_qubits[1]):
                dt = None
                q_l = self.num_qubits[1] - 1 - j
                qubit_l = self.links[q_l][1]
                # the first results are themselves the changes
                if t == 0:
                    change = syndrome_list[-1][j] != "0"
                # if the link was involved in a just finished 202...
                elif just_finished:
                    # skip back self.rounds_per_202 for a neighbouring link
                    if qubit_l in last_neighbors:
                        dt = self.rounds_per_202
                    # and just 1 for all others (as normal)
                    else:
                        dt = 1
                # otherwise, everything not during a 202 just compares
                # results with the previous round (as normal)
                elif tau not in range(1, self.rounds_per_202):
                    dt = 1
                # and all others depend on the placement of the link
                # within the current 202
                else:
                    # if this link is the 202 link
                    if qubit_l == qubit_l_202:
                        if tau == 1:
                            change = False
                        else:
                            dt = 2
                    # if this link neighbours the 202 link
                    elif qubit_l in all_neighbors:
                        change = False
                    # if the link is not near the 202 link (or there are no 202s)
                    else:
                        dt = 1
                # for those where we now have a dt, calculate the change
                if dt:
                    change = syndrome_list[-t - 1][j] != syndrome_list[-t - 1 + dt][j]
                syndrome_changes += "0" * (not change) + "1" * change
            syndrome_changes += " "
            last_neighbors = all_neighbors.copy()
            just_finished = tau == (self.rounds_per_202 - 1)

        # the space separated string of syndrome changes then gets a
        # double space separated logical value on the end
        new_string = measured_log + " " + syndrome_changes[:-1]

        return new_string

    def string2raw_logicals(self, string):
        """
        Extracts raw logicals from output string.
        Args:
            string (string): Results string from which to extract logicals
        Returns:
            list: Raw values for logical operators that correspond to nodes.
        """
        return _separate_string(self._process_string(string))[0]

    def string2nodes(self, string, **kwargs) -> List[DecodingGraphNode]:
        """Convert output string from circuits into a set of nodes.

        Args:
            string (string): Results string to convert.
            kwargs (dict): Additional keyword arguments. See below.

        kwargs:
            all_logicals (bool): Whether to include logical nodes irrespective
                of value. (False as default).

        Returns:
            dict: List of nodes corresponding to to the non-trivial
            elements in the string.
        """

        all_logicals = kwargs.get("all_logicals")
        string = self._process_string(string)
        separated_string = _separate_string(string)
        nodes = []
        for syn_type, _ in enumerate(separated_string):
            for syn_round in range(len(separated_string[syn_type])):
                elements = separated_string[syn_type][syn_round]
                for elem_num, element in enumerate(elements):
                    if (syn_type == 0 and (all_logicals or element != self.logical)) or (
                        syn_type != 0 and element == "1"
                    ):
                        is_logical = syn_type == 0
                        if is_logical:
                            elem_num = syn_round
                            syn_round = 0
                            code_qubits = [self.z_logicals[elem_num]]
                            link_qubit = None
                        else:
                            link = self.links[-elem_num - 1]
                            code_qubits = [link[0], link[2]]
                            link_qubit = link[1]
                        tau, _, _ = self._get_202(syn_round)
                        if not tau:
                            tau = 0
                        node = DecodingGraphNode(
                            is_logical=is_logical,
                            is_boundary=(is_logical and self._leaves),
                            time=syn_round if not is_logical else None,
                            qubits=code_qubits,
                            index=elem_num,
                        )
                        node.properties["conjugate"] = ((tau % 2) == 1) and tau > 1
                        node.properties["link qubit"] = link_qubit
                        nodes.append(node)
        return nodes

    @staticmethod
    def flatten_nodes(nodes: List[DecodingGraphNode]):
        """
        Removes time information from a set of nodes, and consolidates those on
        the same position at different times. Also removes nodes corresponding
        to the conjugate error from [[2,0,2]]s.
        Args:
            nodes (list): List of nodes, of the type produced by `string2nodes`, to be flattened.
        Returns:
            flat_nodes (list): List of flattened nodes.
        """
        # strip out conjugate nodes
        non_conj_nodes = []
        for node in nodes:
            if not node.properties["conjugate"]:
                non_conj_nodes.append(node)
        nodes = non_conj_nodes
        # remove time info
        nodes_per_link = {}
        for node in nodes:
            link_qubit = node.properties["link qubit"]
            if link_qubit in nodes_per_link:
                nodes_per_link[link_qubit] += 1
            else:
                nodes_per_link[link_qubit] = 1
        flat_nodes = []
        for node in nodes:
            if node.is_logical or node.is_boundary:
                flat_nodes.append(node)
            elif nodes_per_link[node.properties["link qubit"]] % 2:
                flat_node = deepcopy(node)
                flat_node.time = None
                if flat_node not in flat_nodes:
                    flat_nodes.append(flat_node)

        return flat_nodes

    def _links2cpp(self):
        """
        Convert data about the link graph to the form required by C++ functions.
        """
        nodes = self.link_graph.nodes()
        link_graph = []
        for edge in self.link_graph.edge_list():
            link_graph.append((nodes[edge[0]], nodes[edge[1]]))
        link_neighbors = {}
        for n, node in enumerate(self.link_graph.nodes()):
            link_neighbors[node] = []
            for j in self.link_graph.neighbors(n):
                link_neighbors[node].append(nodes[j])
        return link_graph, link_neighbors

    def _extras2cpp(self):
        """
        Returns logical and boundary nodes as tuples. First value is the qubit,
        second is 1 for logical only, 2 for boundary only and 3 for both.
        """
        extras = {}
        for q in self.z_logicals:
            extras[q] = 1 + 2 * (q in self.boundary)
        for q in self.boundary:
            if q not in self.z_logicals:
                extras[q] = 2
        return extras

    def check_nodes(self, nodes, ignore_extras=False, minimal=False):
        """
        Determines whether a given set of nodes are neutral. If so, also
        determines any additional logical readout qubits that would be
        flipped by the errors creating such a cluster and how many errors
        would be required to make the cluster.
        Args:
            nodes (list): List of nodes, of the type produced by `string2nodes`.
            ignore_extras (bool): If `True`, undeeded boundary and logical nodes are
            ignored.
            minimal (bool): Whether output should only reflect the minimal error
            case.
        Returns:
            neutral (bool): Whether the nodes independently correspond to a valid
            set of errors.
            flipped_logical_nodes (list): List of qubits nodes for logical
            operators that are flipped by the errors, that were not included
            in the original nodes.
            num_errors (int): Minimum number of errors required to create nodes.
        """

        nodes = _nodes2cpp(nodes)

        cpp_output = _c_check_nodes(
            nodes,
            ignore_extras,
            minimal,
            self.cycle_dict,
            self._cpp_link_graph,
            self._cpp_link_neighbors,
            self._extras2cpp(),
        )

        neutral = bool(cpp_output[0])
        num_errors = cpp_output[1]
        flipped_extra_nodes = []
        for flipped_extra in cpp_output[2::]:
            is_logical = flipped_extra in self.z_logicals
            is_boundary = flipped_extra in self.boundary
            if is_logical:
                index = self.z_logicals.index(flipped_extra)
            else:
                index = self.boundary.index(flipped_extra)
            node = DecodingGraphNode(
                is_logical=is_logical,
                is_boundary=is_boundary,
                qubits=[flipped_extra],
                index=index,
            )
            flipped_extra_nodes.append(node)

        return neutral, flipped_extra_nodes, num_errors

    def is_cluster_neutral(self, nodes: dict):
        """
        Determines whether or not the cluster is neutral, meaning that one or more
        errors could have caused the set of nodes (syndrome changes) passed
        to the method.
        Args:
            nodes: dictionary in the form of the return value of string2nodes
        """
        nodes = _nodes2cpp(nodes)
        return _c_is_cluster_neutral(
            nodes,
            False,
            False,
            self.cycle_dict,
            self._cpp_link_graph,
            self._cpp_link_neighbors,
            self._extras2cpp(),
            self._linear,
        )

    def transpile(self, backend, echo=("X", "X"), echo_num=(2, 0)):
        """
        Args:
            backend (qiskit.providers.ibmq.IBMQBackend): Backend to transpile and schedule the
            circuits for. The numbering of the qubits in this backend should correspond to
            the numbering used in `self.links`.
            echo (tuple): List of gate sequences (expressed as strings) to be used on code qubits and
            link qubits, respectively. Valid strings are `'X'` and `'XZX'`.
            echo_num (tuple): Number of times to repeat the sequences for code qubits and
            link qubits, respectively.
        Returns:
            transpiled_circuit: As `self.circuit`, but with the circuits scheduled, transpiled and
            with dynamical decoupling added.
        """

        bases = list(self.circuit.keys())
        circuits = [self.circuit[basis] for basis in bases]

        initial_layout = []
        for qreg in circuits[0].qregs:
            qreg_index = int("link" in qreg.name)
            initial_layout += [
                self.qubits[qreg_index][q] for q in range(self.num_qubits[qreg_index])
            ]

        # transpile to backend
        circuits = transpile(circuits, backend, initial_layout=initial_layout)

        # then dynamical decoupling if needed
        if any(echo_num):
            if self.run_202:
                durations = DynamicCircuitInstructionDurations().from_backend(backend)
            else:
                durations = InstructionDurations().from_backend(backend)

            # set up the dd sequences
            dd_sequences = []
            spacings = []
            for j in range(2):
                if echo[j] == "X":
                    dd_sequences.append([XGate()] * echo_num[j])
                    spacings.append(None)
                elif echo[j] == "XZX":
                    dd_sequences.append([XGate(), RZGate(np.pi), XGate()] * echo_num[j])
                    d = 1.0 / (2 * echo_num[j] - 1 + 1)
                    spacing = [d / 2] + ([0, d, d] * echo_num[j])[:-1] + [d / 2]
                    for _ in range(2):
                        spacing[0] += 1 - sum(spacing)
                    spacings.append(spacing)
                else:
                    dd_sequences.append(None)
                    spacings.append(None)

            # add in the dd sequences
            for j, dd_sequence in enumerate(dd_sequences):
                if dd_sequence:
                    if echo_num[j]:
                        qubits = self.qubits[j]
                    else:
                        qubits = None
                    pm = PassManager(
                        [
                            ALAPScheduleAnalysis(durations),
                            PadDynamicalDecoupling(
                                durations, dd_sequence, qubits=qubits, spacings=spacings[j]
                            ),
                        ]
                    )
                    circuits = pm.run(circuits)
            if not isinstance(circuits, list):
                circuits = [circuits]

        return {basis: circuits[j] for j, basis in enumerate(bases)}

    def _make_syndrome_graph(self):
        # get the list of nodes
        string = (
            "1" * len(self.code_qubit)
            + " "
            + ("0" * len(self.links) + " ") * (self.T - 1)
            + "1" * len(self.links)
        )
        nodes: List[DecodingGraphNode] = []
        for node in self.string2nodes(string, all_logicals=True):
            if not node.is_logical:
                for t in range(self.T + 1):
                    new_node = deepcopy(node)
                    new_node.time = t
                    if new_node not in nodes:
                        nodes.append(new_node)
            else:
                node.time = None
                nodes.append(node)

        # find pairs that should be connected
        edges: List[Tuple[int, int]] = []
        for n0, node0 in enumerate(nodes):
            for n1, node1 in enumerate(nodes):
                if n0 < n1:
                    # just record all possible edges for now
                    dt = abs((node1.time or 0) - (node0.time or 0))
                    adj = set(node0.qubits).intersection(set(node1.qubits))
                    if adj:
                        if (node0.is_logical ^ node1.is_logical) or dt <= 1:
                            edges.append((n0, n1))
                        elif not self.resets:
                            if node0.qubits == node1.qubits and dt == 2:
                                edges.append((n0, n1))

        # put it all in a graph
        S = rx.PyGraph(multigraph=False)
        for node in nodes:
            S.add_node(node)
        for n0, n1 in edges:
            source = nodes[n0]
            target = nodes[n1]
            qubits = []
            if not (source.is_logical and target.is_logical):
                qubits = list(set(source.qubits).intersection(target.qubits))
            if source.time != target.time and len(qubits) > 1:
                qubits = []
            edge = DecodingGraphEdge(qubits=qubits, weight=1)
            S.add_edge(n0, n1, edge)

        # remove invalid edges
        self.get_error_coords(None, S, remove_invalid_edges=True)

        # just record edges as hyperedges for now (should be improved later)
        hyperedges = []
        for e, n0n1 in enumerate(S.edge_list()):
            hyperedges.append({n0n1: S.edges()[e]})

        return S, hyperedges

    def get_error_coords(
        self,
        counts,
        decoding_graph,
        method="spitz",
        remove_invalid_edges=False,
        return_samples=False,
    ):
        """
        Uses the `get_error_probs` method of the given decoding graph to generate probabilities
        of single error events from given counts. The location and time of each error is
        also calculated.

        Args:
            counts (dict): Counts dictionary of the results to be analyzed.
            decoding_graph (DecodingGraph): Decoding graph object constructed
            from this code.
            method (string): Method to used for calculation. Supported
            methods are 'spitz' (default) and 'naive'.
            remove_invalid_edges (string): Whether to delete edges from the graph if
            they are found to be invalid.
            return_samples (bool): Whether to also return the number of
            samples used to calculated each probability.
        Returns:
            dict: Keys are the coordinates (qubit, start_time, end_time) for specific error
            events. Time refers to measurement rounds. Values are a dictionary whose keys are
            the edges that detected the event, and whose keys are the calculated probabilities.
        Additional information:
            Time calculation does not take into account get lengths. It assumes that the
            subrounds within the schedule and the measurement all take the same time. Time
            is in units of rounds.
        """

        # though the documented use case requires a decoding graph and a counts dict, there is also an
        # undocumented internal use case, where just the bare graph is provided and no counts. This is
        # to find and delete invalid edges
        if isinstance(decoding_graph, rx.PyGraph):
            graph = decoding_graph
        else:
            graph = decoding_graph.graph
        nodes = graph.nodes()
        if counts:
            if return_samples:
                error_probs, samples = decoding_graph.get_error_probs(
                    counts, method=method, return_samples=True
                )
            else:
                error_probs = decoding_graph.get_error_probs(counts, method=method)
        else:
            error_probs = {}
            for n0, n1 in graph.edge_list():
                if nodes[n0].is_logical:
                    edge = (n1, n1)
                elif nodes[n1].is_logical:
                    edge = (n0, n0)
                else:
                    edge = (n0, n1)
                error_probs[edge] = np.nan

        if hasattr(self, "z_logicals"):
            z_logicals = set(self.z_logicals)
        elif hasattr(self, "z_logical"):
            z_logicals = {self.z_logical}
        else:
            print("No qubits for z logicals found. Proceeding without.")
            z_logicals = set()

        round_length = len(self.schedule) + 1

        error_coords = {}
        sample_coords = {}
        for (n0, n1), prob in error_probs.items():
            node0 = nodes[n0]
            node1 = nodes[n1]
            if n0 != n1:
                qubits = graph.get_edge_data(n0, n1).qubits
                if qubits:
                    # error on a code qubit between rounds, or during a round
                    assert (node0.time == node1.time and node0.qubits != node1.qubits) or (
                        node0.time != node1.time and node0.qubits != node1.qubits
                    )
                    qubit = qubits[0]
                    # error between rounds
                    if node0.time == node1.time:
                        dts = []
                        for node in [node0, node1]:
                            pair = [qubit, node.properties["link qubit"]]
                            for dt, pairs in enumerate(self.schedule):
                                if pair in pairs or tuple(pair) in pairs:
                                    dts.append(dt)
                        time = [max(0, node0.time - 1 + (max(dts) + 1) / round_length)]
                        time.append(min(self.T, node0.time + min(dts) / round_length))
                    # error during a round
                    else:
                        # put nodes in descending time order
                        if node0.time < node1.time:
                            node_pair = [node1, node0]
                        else:
                            node_pair = [node0, node1]
                        # see when in the schedule each node measures the qubit
                        dts = []
                        for node in node_pair:
                            pair = [qubit, node.properties["link qubit"]]
                            for dt, pairs in enumerate(self.schedule):
                                if pair in pairs or tuple(pair) in pairs:
                                    dts.append(dt)
                        # use to define fractional time
                        if dts[0] < dts[1]:
                            time = [node_pair[1].time + (dts[0] + 1) / round_length]
                            time.append(node_pair[1].time + dts[1] / round_length)
                        else:
                            # impossible cases get no valid time
                            time = []
                            if remove_invalid_edges:
                                graph.remove_edge(n0, n1)
                else:
                    # measurement error
                    assert node0.time != node1.time and node0.qubits == node1.qubits
                    qubit = node0.properties["link qubit"]
                    t0 = min(node0.time, node1.time)
                    if abs(node0.time - node1.time) == 1:
                        if self.resets:
                            time = [t0, t0 + 1]
                        else:
                            time = [t0, t0 + (round_length - 1) / round_length]
                    else:
                        time = [t0 + (round_length - 1) / round_length, t0 + 1]
            else:
                # detected only by one stabilizer
                boundary_qubits = list(set(node0.qubits).intersection(z_logicals))
                # for the case of boundary stabilizers
                if boundary_qubits:
                    qubit = boundary_qubits[0]
                    pair = [qubit, node0.properties["link qubit"]]
                    for dt, pairs in enumerate(self.schedule):
                        if pair in pairs or tuple(pair) in pairs:
                            time = [max(0, node0.time - 1 + (dt + 1) / round_length)]
                            time.append(min(self.T, node0.time + dt / round_length))

                else:
                    qubit = tuple(node0.qubits + [node0.properties["link qubit"]])
                    time = [node0.time, node0.time + (round_length - 1) / round_length]

            if time:  # only record if not nan
                if (qubit, time[0], time[1]) not in error_coords:
                    error_coords[qubit, time[0], time[1]] = {}
                    sample_coords[qubit, time[0], time[1]] = {}
                error_coords[qubit, time[0], time[1]][n0, n1] = prob
                if return_samples:
                    sample_coords[qubit, time[0], time[1]][n0, n1] = samples[n0, n1]

        if return_samples:
            return error_coords, sample_coords
        else:
            return error_coords

    def clean_code(self, string):
        """
        Given an output string of the code, obvious code qubit errors are identified and their effects
        are removed.

        Args:
            string (str): Output string of the code.

        Returns:
            string (str): Modifed output string of the code.
        """

        # get the parities for the rounds and turn them into lists of integers
        # (also turn them the right way around)
        parities = []
        for rstring in string.split(" ")[1:]:
            parities.append([int(p) for p in rstring][::-1])
        parities = parities[::-1]

        # calculate the final parities from the final readout and add them on
        final = string.split(" ")[0]
        final_parities = [0] * self.num_qubits[1]
        for c0, a, c1 in self.links:
            final_parities[-self.link_index[a] - 1] = (
                int(final[-self.code_index[c0] - 1]) + int(final[-self.code_index[c1] - 1])
            ) % 2
        parities.append(final_parities[::-1])

        flips = {c: 0 for c in self.code_index}
        for rparities in parities:
            # see how many links around each code qubit detect a flip
            link_count = {c: 0 for c in self.code_index}
            for c0, a, c1 in self.links:
                # we'll need to determine whether the as yet uncorrected parity
                # checks from this round should be flipped, based on results
                # from previous rounds
                flip = (flips[c0] + flips[c1]) % 2
                b = self.link_index[a]
                for c in [c0, c1]:
                    link_count[c] += (rparities[b] + flip) % 2
            # if it's all of them, assume a flip
            for c in link_count:
                if link_count[c] == self.degree[c]:
                    flips[c] = (flips[c] + 1) % 2
            # modify the parities to remove the effect
            for c0, a, c1 in self.links:
                flip = (flips[c0] + flips[c1]) % 2
                b = self.link_index[a]
                rparities[b] = (rparities[b] + flip) % 2
        # turn the results back into a string
        new_string = ""
        for rparities in parities[:-1][::-1]:
            new_string += " " + "".join([str(p) for p in rparities][::-1])
        final_string = [int(p) for p in string.split(" ", maxsplit=1)[0]]
        for c, flip in flips.items():
            b = self.code_index[c]
            final_string[-b - 1] = (final_string[-b - 1] + flip) % 2
        final_string = "".join([str(p) for p in final_string])

        return final_string + new_string
