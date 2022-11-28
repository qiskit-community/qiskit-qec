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
from typing import List, Optional, Tuple

import numpy as np
import rustworkx as rx

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, transpile
from qiskit.circuit.library import XGate, RZGate
from qiskit.transpiler import PassManager, InstructionDurations
from qiskit.transpiler.passes import DynamicalDecoupling


class RepetitionCodeCircuit:
    """RepetitionCodeCircuit class."""

    def __init__(
        self,
        d: int,
        T: Optional[int] = None,
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

        gauge_ops = [[j, j + 1] for j in range(self.d - 1)]
        measured_logical = [[0]]
        flip_logical = list(range(self.d))
        boundary = [[0], [self.d - 1]]

        if xbasis:
            self.css_x_gauge_ops = gauge_ops
            self.css_x_stabilizer_ops = gauge_ops
            self.css_x_logical = measured_logical
            self.css_x_boundary = boundary
            self.css_z_gauge_ops = []
            self.css_z_stabilizer_ops = []
            self.css_z_logical = flip_logical
            self.css_z_boundary = []
            self.basis = "x"
        else:
            self.css_x_gauge_ops = []
            self.css_x_stabilizer_ops = []
            self.css_x_logical = flip_logical
            self.css_x_boundary = []
            self.css_z_gauge_ops = gauge_ops
            self.css_z_stabilizer_ops = gauge_ops
            self.css_z_logical = measured_logical
            self.css_z_boundary = boundary
            self.basis = "z"
        self.round_schedule = self.basis
        self.blocks = T

        self.resets = resets
        self.delay = delay

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

    def _separate_string(self, string):

        separated_string = []
        for syndrome_type_string in string.split("  "):
            separated_string.append(syndrome_type_string.split(" "))
        return separated_string

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

    def string2nodes(self, string, logical="0", all_logicals=False):
        """
        Convert output string from circuits into a set of nodes.
        Args:
            string (string): Results string to convert.
            logical (string): Logical value whose results are used.
            all_logicals (bool): Whether to include logical nodes
            irrespective of value.
        Returns:
            dict: List of nodes corresponding to to the non-trivial
            elements in the string.

        Additional information:
        Strings are read right to left, but lists*
        are read left to right. So, we have some ugly indexing
        code whenever we're dealing with both strings and lists.
        """

        string = self._process_string(string)
        separated_string = self._separate_string(string)  # [ <boundary>, <syn>, <syn>,...]
        nodes = []

        # boundary nodes
        boundary = separated_string[0]  # [<last_elem>, <init_elem>]
        for bqec_index, belement in enumerate(boundary[::-1]):
            if all_logicals or belement != logical:
                bnode = {"time": 0}
                i = [0, -1][bqec_index]
                if self.basis == "z":
                    bqubits = [self.css_x_logical[i]]
                else:
                    bqubits = [self.css_z_logical[i]]
                bnode["qubits"] = bqubits
                bnode["is_boundary"] = True
                bnode["element"] = bqec_index
                nodes.append(bnode)

        # bulk nodes
        for syn_type in range(1, len(separated_string)):
            for syn_round in range(len(separated_string[syn_type])):
                elements = separated_string[syn_type][syn_round]
                for qec_index, element in enumerate(elements[::-1]):
                    if element == "1":
                        node = {"time": syn_round}
                        if self.basis == "z":
                            qubits = self.css_z_gauge_ops[qec_index]
                        else:
                            qubits = self.css_x_gauge_ops[qec_index]
                        node["qubits"] = qubits
                        node["is_boundary"] = False
                        node["element"] = qec_index
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
        return self._separate_string(self._process_string(string))[0]

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


class ArcCircuit:
    """Anisotropic repetition code class."""

    def __init__(
        self,
        links: list,
        T: int,
        basis: str = "xy",
        resets: bool = True,
        ff: bool = True,
        delay: Optional[int] = None,
        barriers: bool = False,
        color: Optional[dict] = None,
        max_dist: int = 2,
        schedule: Optional[list] = None,
        run_202: bool = True,
    ):
        """
        Creates circuits corresponding to an anisotropic repetition code implemented over T syndrome
        measurement rounds, with the syndrome measurements provided.
        Args:
            links (list): List of tuples (c0, a, c1), where c0 and c1 are the two code qubits in each
            syndrome measurement, and a is the auxiliary qubit used.
            T (int): Number of rounds of syndrome measurement.
            basis (list): Pair of `'x'`, `'y'` and `'z'`, specifying the pair of local bases to be
            used.
            resets (bool): Whether to include a reset gate after mid-circuit measurements.
            ff (bool): Whether to correct the effects of [[2,0,2]] sequences via feed forward.
            delay (float): Time (in dt) to delay after mid-circuit measurements (and delay).
            barriers (bool): Whether to include barriers between different sections of the code.
            color (dict): Dictionary with code qubits as keys and 0 or 1 for each value, to specify
            a predetermined bicoloring. If not provided, a bicoloring is found on initialization.
            max_dist (int): Maximum edge distance used when determining the bicoloring of code qubits.
            schedule(list): Specifies order in which entangling gates are applied in each syndrome
            measurement round. Each element is a list of lists [c, a] for entangling gates to be
            applied simultaneously.
            run_202 (bool): Whether to run [[2,0,2]] sequences. This will be overwritten if T is not high
            enough (at least 5xlen(links)).
        """

        self.links = links
        self.basis = basis
        self._resets = resets
        self._barriers = barriers
        self._max_dist = max_dist
        self.delay = delay or 0

        # calculate coloring and schedule, etc
        if color is None:
            self._coloring()
        else:
            self.color = color
        if schedule is None:
            self._scheduling()
        else:
            self.schedule = schedule
        self._preparation()

        # determine the placement of [2,0,2] rounds
        if run_202:
            self.links_202 = []
            for link in self.links:
                logical_overlap = {link[0], link[2]}.intersection(set(self.z_logicals))
                if not logical_overlap:
                    self.links_202.append(link)
            num_links = len(self.links_202)
            if num_links > 0:
                self.rounds_per_link = np.floor(T / num_links)
                self.metabuffer = np.ceil((T - num_links * self.rounds_per_link) / 2)
                self.roundbuffer = np.ceil((self.rounds_per_link - 5) / 2)
                self.run_202 = self.rounds_per_link >= 5
            else:
                self.run_202 = False
        else:
            self.run_202 = False
        self._ff = ff and self.run_202

        # create the circuit
        self.base = basis
        self.T = 0
        for _ in range(T - 1):
            self._syndrome_measurement()
        if T != 0:
            self._syndrome_measurement(final=True)
        self._readout()

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

        weight_fn = lambda edge: -int(edge["weight"])

        schedule = []
        while aux:

            # construct coupling graph for as yet unpaired auxiliaries (i.e. link qubits)
            graph = self._get_coupling_graph(aux)

            # find a min weight matching, and then another that exlcudes the pairs from the first
            matching = [rx.max_weight_matching(graph, max_cardinality=True, weight_fn=weight_fn)]
            cut_graph = graph.copy()
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
        qubits = [set(), set()]
        for link in self.links:
            qubits[0] = qubits[0].union({link[0], link[2]})
            qubits[1].add(link[1])
        self.qubits = [list(qubits[j]) for j in range(2)]
        self.num_qubits = [len(qubits[j]) for j in range(2)]

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
        for basis in [self.basis, self.basis[::-1]]:
            self.circuit[basis] = QuantumCircuit(self.link_qubit, self.code_qubit, name=basis)
            self._basis_change(basis)

        # use degree 1 code qubits for logical z readouts
        graph = self._get_coupling_graph()
        z_logicals = []
        for n, node in enumerate(graph.nodes()):
            if graph.degree(n) == 1:
                z_logicals.append(node)
        # if there are none, just use the first
        if z_logicals == []:
            z_logicals = [min(self.code_index.keys())]
        self.z_logicals = z_logicals

    def _get_202(self, t):
        """
        Returns the position within a 202 sequence for the current measurement round:
        * `False` means not part of a 202 sequence;
        * 0, 2 and 4 use the standard coloring;
        * 1 and 3 use the flipped coloring.
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
            if t < 0 or tau not in range(5):
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
        qubits_to_correct = set()
        for basis, qc in self.circuit.items():
            if self._barriers:
                qc.barrier()
            for pairs in self.schedule:
                for qubit_c, qubit_l in pairs:
                    q_c = self.code_index[qubit_c]
                    q_l = self.link_index[qubit_l]
                    neighbor = qubit_l in qubit_l_nghbrs[0] + qubit_l_nghbrs[1]
                    if not (tau in [1, 2, 3] and neighbor):
                        c = self.color[qubit_c]
                        if qubit_l == qubit_l_202:
                            c = (c + tau) % 2
                        self._rotate(basis, c, self.code_qubit[q_c], True)
                        qc.cx(self.code_qubit[q_c], self.link_qubit[q_l])
                        self._rotate(basis, c, self.code_qubit[q_c], False)
                        links_to_measure.add(q_l)
                        if type(tau) == int and tau == 0 and neighbor:
                            if not self._ff:
                                links_to_reset.add(q_l)
                        else:
                            links_to_reset.add(q_l)
                        if tau == 4 and neighbor:
                            qubits_to_correct.add(q_l)

        # measurement, etc
        for basis, qc in self.circuit.items():

            # measurement
            if self._barriers:
                if final:
                    qc.barrier(self.link_qubit)
                else:
                    qc.barrier()
            qc.add_register(self.link_bits[self.T])
            for q_l in links_to_measure:
                qc.measure(self.link_qubit[q_l], self.link_bits[self.T][q_l])

            # resets
            if self._resets and not final:
                for q_l in links_to_reset:
                    qc.reset(self.link_qubit[q_l])
                    # might at some point add an option for reset via
                    # qc.x(self.link_qubit[q_l]).c_if(self.link_bits[self.T][q_l], 1)

            # correct
            if self._ff:
                for q_l in qubits_to_correct:
                    link = set(self.links[q_l])
                    link_202 = set(self.links[self.link_index[qubit_l_202]])
                    qubit_c = list(link.intersection(link_202))[0]
                    q_c = self.code_index[qubit_c]
                    c = self.color[qubit_c]
                    self._rotate(basis, c, self.code_qubit[q_c], True)
                    qc.x(self.code_qubit[q_c]).c_if(self.link_bits[self.T][q_l], 1)
                    self._rotate(basis, c, self.code_qubit[q_c], False)

            # delay
            if self.delay > 0 and not final:
                if self._barriers:
                    qc.barrier()
                qc.delay(self.delay, self.link_qubit)

        self.T += 1

    def _readout(self):
        """
        Readout of all code qubits, which corresponds to a logical measurement
        as well as allowing for a measurement of the syndrome to be inferred.
        """

        for basis, qc in self.circuit.items():
            self._basis_change(basis, inverse=True)
            if self._barriers:
                qc.barrier(self.code_qubit)
            qc.add_register(self.code_bit)
            qc.measure(self.code_qubit, self.code_bit)

    def _separate_string(self, string):
        separated_string = []
        for syndrome_type_string in string.split("  "):
            separated_string.append(syndrome_type_string.split(" "))
        return separated_string

    def _process_string(self, string):

        # logical readout taken from assigned qubits
        measured_log = ""
        for qubit in self.z_logicals:
            j = self.code_index[qubit]
            measured_log += string[self.num_qubits[0] - j - 1] + " "

        if self._resets:
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

        # changes between one syndrome and the next then calculated
        syndrome_list = full_syndrome.split(" ")
        syndrome_changes = ""
        last_neighbors = []
        just_finished = False
        for t in range(self.T + 1):
            tau, qubit_l_202, qubit_l_nghbrs = self._get_202(t)
            all_neighbors = qubit_l_nghbrs[0] + qubit_l_nghbrs[1]
            for j in range(self.num_qubits[1]):
                q_l = self.num_qubits[1] - 1 - j
                qubit_l = self.links[q_l][1]
                if qubit_l in all_neighbors and tau in [1, 2, 3]:
                    # don't calculate changes for neighbours of the 202
                    change = False
                elif qubit_l in all_neighbors and tau == 4:
                    # index for neighbours on the other side of the 202 link
                    opp_index = int(qubit_l in qubit_l_nghbrs[0])
                    # first listed link on the other side
                    qubit_l_opp = qubit_l_nghbrs[opp_index][0]
                    # position in register for this link
                    k = self.num_qubits[1] - 1 - self.link_index[qubit_l_opp]
                    # determine change for product of these two over the past four rounds
                    changes = 0
                    for jj in [j, k]:
                        dt = 4
                        changes += syndrome_list[-t - 1][jj] != syndrome_list[-t - 1 + dt][jj]
                    change = changes % 2 == 1
                else:
                    if t == 0:
                        change = syndrome_list[-1][j] != "0"
                    else:
                        if qubit_l == qubit_l_202:
                            if tau == 1:
                                change = False
                            else:
                                if tau in [2, 3, 4]:
                                    dt = 2
                                else:
                                    if self._ff and qubit_l in last_neighbors and just_finished:
                                        dt = 5
                                    else:
                                        dt = 1
                                change = syndrome_list[-t - 1][j] != syndrome_list[-t - 1 + dt][j]
                        elif qubit_l in last_neighbors and just_finished:
                            if self._ff:
                                dt = 5
                            else:
                                dt = 1
                            change = syndrome_list[-t - 1][j] != syndrome_list[-t - 1 + dt][j]
                        else:
                            change = syndrome_list[-t - 1][j] != syndrome_list[-t][j]

                syndrome_changes += "0" * (not change) + "1" * change
            syndrome_changes += " "
            last_neighbors = all_neighbors.copy()
            just_finished = tau == 4

        # the space separated string of syndrome changes then gets a
        # double space separated logical value on the end
        new_string = measured_log + " " + syndrome_changes[:-1]

        return new_string

    def string2nodes(self, string, logical="0", all_logicals=False):
        """
        Convert output string from circuits into a set of nodes.
        Args:
            string (string): Results string to convert.
            logical (string): Logical value whose results are used.
            all_logicals (bool): Whether to include logical nodes
            irrespective of value.
        Returns:
            dict: List of nodes corresponding to to the non-trivial
            elements in the string.
        """

        string = self._process_string(string)
        separated_string = self._separate_string(string)
        nodes = []
        for syn_type, _ in enumerate(separated_string):
            for syn_round in range(len(separated_string[syn_type])):
                elements = separated_string[syn_type][syn_round]
                for elem_num, element in enumerate(elements):
                    if (syn_type == 0 and (all_logicals or element != logical)) or (
                        syn_type != 0 and element == "1"
                    ):
                        is_boundary = syn_type == 0
                        if is_boundary:
                            elem_num = syn_round
                            syn_round = 0
                            code_qubits = [self.z_logicals[-elem_num]]
                            link_qubit = None
                        else:
                            link = self.links[-elem_num - 1]
                            code_qubits = [link[0], link[2]]
                            link_qubit = link[1]
                        node = {"time": syn_round}
                        node["qubits"] = code_qubits
                        node["link qubit"] = link_qubit
                        node["is_boundary"] = is_boundary
                        node["element"] = elem_num
                        nodes.append(node)
        return nodes

    def flatten_nodes(self, nodes):
        """
        Removes time information from a set of nodes, and consolidates those on
        the same position at different times.
        Args:
            nodes (dict): List of nodes, of the type produced by `string2nodes`, to be flattened.
        Returns:
            flat_nodes (dict): List of flattened nodes.
        """
        nodes_per_link = {}
        for node in nodes:
            link_qubit = node["link qubit"]
            if link_qubit in nodes_per_link:
                nodes_per_link[link_qubit] += 1
            else:
                nodes_per_link[link_qubit] = 1
        flat_nodes = []
        for node in nodes:
            if nodes_per_link[node["link qubit"]] % 2:
                flat_node = node.copy()
                if "time" in flat_node:
                    flat_node.pop("time")
                flat_nodes.append(flat_node)
        return flat_nodes

    def check_nodes(self, nodes):
        """
        Determines whether a given set of nodes are neutral. If so, also
        determines the logical readout qubits they contain.
        Args:
            nodes (list): List of nodes, of the type produced by `string2nodes`.
        Returns:
            neutral (bool): Whether the nodes independently correspond to a valid
            set of errors.
            flipped_logicals (list): List of qubits within `z_logicals`
            enclosed by the nodes.
        """
        nodes = self.flatten_nodes(nodes)
        link_qubits = set(node["link qubit"] for node in nodes)
        node_color = {0: 0}
        neutral = True
        link_graph = self._get_link_graph()
        ns_to_do = set(n for n in range(1, len(link_graph.nodes())))
        n = 0
        while ns_to_do and neutral:
            for n in link_graph.neighbors(n):
                if n in node_color:
                    incident_es = link_graph.incident_edges(n)
                    for e in incident_es:
                        edge = link_graph.edges()[e]
                        n0, n1 = link_graph.edge_list()[e]
                        if n0 == n:
                            nn = n1
                        else:
                            nn = n0
                        dc = edge["link qubit"] in link_qubits
                        if nn not in node_color:
                            node_color[nn] = (node_color[n] + dc) % 2
                            ns_to_do.remove(nn)
                        else:
                            neutral = neutral and (node_color[nn] == (node_color[n] + dc) % 2)

        flipped_logicals = []
        if neutral:
            inside_c = int(sum(node_color.values()) < len(node_color) / 2)
            for n, c in node_color.items():
                node = link_graph.nodes()[n]
                if node in self.z_logicals and c == inside_c:
                    flipped_logicals.append(node)

        return neutral, flipped_logicals

    def transpile(self, backend, echo=("X", "X"), echo_num=(2, 0)):
        """
        Args:
            backend (qiskit.providers.ibmq.IBMQBackend): Backend to transpile and schedule the
            circuits for. The numbering of the qubits in this backend should correspond to
            the numbering used in `self.links`.
            echo (tuple): List of gate sequences (expressed as strings) to be used on code qubits and
            link qubits, respectively. Valid strings are `'X'` and `'XZX'`.
            echo_num(tuple): Number of times to repeat the sequences (as a list) for code qubits and
            link qubits, respectively.
        Returns:
            transpiled_circuit: As `self.circuit`, but with the circuits scheduled, transpiled and
            with dynamical decoupling added.
        """

        circuits = [self.circuit[basis] for basis in [self.basis, self.basis[::-1]]]

        initial_layout = []
        for qreg in circuits[0].qregs:
            qreg_index = int("link" in qreg.name)
            initial_layout += [
                self.qubits[qreg_index][q] for q in range(self.num_qubits[qreg_index])
            ]

        # transpile to backend and schedule
        circuits = transpile(
            circuits, backend, initial_layout=initial_layout, scheduling_method="alap"
        )

        # then dynamical decoupling if needed
        if any(echo_num):

            # set up the dd sequences
            dd_sequences = []
            spacings = []
            for j in range(2):
                if echo[j] == "X":
                    dd_sequences.append([XGate()] * echo_num[j])
                    spacings.append(None)
                elif echo[j] == "XZX":
                    dd_sequences.append([XGate(), RZGate(np.pi), XGate()] * echo_num)
                    d = 1.0 / (2 * echo_num - 1 + 1)
                    spacing = [d / 2] + ([0, d, d] * echo_num[j])[:-1] + [d / 2]
                    for _ in range(2):
                        spacing[0] += 1 - sum(spacing)
                    spacings.append(spacing)
                else:
                    dd_sequences.append(None)
                    spacings.append(None)

            # add in the dd sequences
            durations = InstructionDurations().from_backend(backend)
            for j, dd_sequence in enumerate(dd_sequences):
                if dd_sequence:
                    if echo_num[j]:
                        qubits = self.qubits[j]
                    else:
                        qubits = None
                    pm = PassManager(
                        [
                            DynamicalDecoupling(
                                durations, dd_sequence, qubits=qubits, spacing=spacings[j]
                            )
                        ]
                    )
                    circuits = pm.run(circuits)

            # make sure delays are a multiple of 16 samples, while keeping the barriers
            # as aligned as possible
            for qc in circuits:
                total_delay = [{q: 0 for q in qc.qubits} for _ in range(2)]
                for gate in qc.data:
                    if gate[0].name == "delay":
                        q = gate[1][0]
                        t = gate[0].params[0]
                        total_delay[0][q] += t
                        new_t = 16 * np.ceil((total_delay[0][q] - total_delay[1][q]) / 16)
                        total_delay[1][q] += new_t
                        gate[0].params[0] = new_t

            # transpile to backend and schedule again
            circuits = transpile(circuits, backend, scheduling_method="alap")

        return {basis: circuits[j] for j, basis in enumerate([self.basis, self.basis[::-1]])}

    def _make_syndrome_graph(self):

        # get the list of nodes
        string = (
            "1" * len(self.code_qubit)
            + " "
            + ("0" * len(self.links) + " ") * (self.T - 1)
            + "1" * len(self.links)
        )
        nodes = []
        for node in self.string2nodes(string):
            if not node["is_boundary"]:
                for t in range(self.T + 1):
                    new_node = node.copy()
                    new_node["time"] = t
                    if new_node not in nodes:
                        nodes.append(new_node)
            else:
                node["time"] = 0
                nodes.append(node)

        # find pairs that should be connected
        edges = []
        for n0, node0 in enumerate(nodes):
            for n1, node1 in enumerate(nodes):
                if n0 < n1:
                    # just record all possible edges for now (should be improved later)
                    dt = abs(node1["time"] - node0["time"])
                    adj = set(node0["qubits"]).intersection(set(node1["qubits"]))
                    if adj:
                        if (node0["is_boundary"] ^ node1["is_boundary"]) or dt <= 1:
                            edges.append((n0, n1))

        # put it all in a graph
        S = rx.PyGraph(multigraph=False)
        hyperedges = []
        for node in nodes:
            S.add_node(node)
        for n0, n1 in edges:
            source = nodes[n0]
            target = nodes[n1]
            qubits = []
            if not (source["is_boundary"] and target["is_boundary"]):
                qubits = list(set(source["qubits"]).intersection(target["qubits"]))
            if source["time"] != target["time"] and len(qubits) > 1:
                qubits = []
            edge = {"qubits": qubits, "weight": 1}
            S.add_edge(n0, n1, edge)
            # just record edges as hyperedges for now (should be improved later)
            hyperedges.append({(n0, n1): edge})

        return S, hyperedges
