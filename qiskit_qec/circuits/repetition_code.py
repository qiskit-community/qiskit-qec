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

"""Generates circuits for quantum error correction."""
from typing import Optional, List, Tuple

from qiskit import QuantumRegister, ClassicalRegister
from qiskit import QuantumCircuit


class RepetitionCodeCircuit:
    """RepetitionCodeCircuit class."""

    def __init__(
        self,
        d: int,
        T: Optional[int] = None,
        xbasis: bool = False,
        resets: bool = False,
        delay: Optional[int] = None,
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
                the end.
        """
        for log in logs:
            if self._xbasis:
                self.circuit[log].z(self.code_qubit)
            else:
                self.circuit[log].x(self.code_qubit)
            if barrier:
                self.circuit[log].barrier()

    def _preparation(self, barrier=False):
        """Prepares logical bit states by applying an x to the circuit that will
        encode a 1.
        """

        for log in ["0", "1"]:
            if self._xbasis:
                self.circuit[log].h(self.code_qubit)
            if barrier:
                self.circuit[log].barrier()

        self.x(["1"])

    def syndrome_measurement(self, final: bool = False, barrier: bool = False, delay: int = 0):
        """Application of a syndrome measurement round.

        Args:
            final (bool): Whether to disregard the reset (if applicable) due to this
            being the final syndrome measurement round.
            barrier (bool): Boolean denoting whether to include a barrier at the end.
            delay (float): Time (in dt) to delay after mid-circuit measurements (and delay).
        """
        self.link_bits.append(ClassicalRegister((self.d - 1), "round_" + str(self.T) + "_link_bit"))

        for log in ["0", "1"]:

            self.circuit[log].add_register(self.link_bits[-1])

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

            for j in range(self.d - 1):
                self.circuit[log].measure(self.link_qubit[j], self.link_bits[self.T][j])
                if self._resets and not final:
                    self.circuit[log].reset(self.link_qubit[j])
                if delay > 0 and not final:
                    self.circuit[log].delay(delay, self.link_qubit[j])

            if barrier:
                self.circuit[log].barrier()

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

        # final syndrome deduced from final code qubit readout
        full_syndrome = ""
        for j in range(self.d - 1):
            full_syndrome += "0" * (string[j] == string[j + 1]) + "1" * (string[j] != string[j + 1])
        # results from all other syndrome measurements then added
        full_syndrome = full_syndrome + string[self.d :]

        # changes between one syndrome and the next then calculated
        syndrome_list = full_syndrome.split(" ")
        syndrome_changes = ""
        for t in range(self.T + 1):
            for j in range(self.d - 1):
                if self._resets:
                    if t == 0:
                        change = syndrome_list[-1][j] != "0"
                    else:
                        change = syndrome_list[-t - 1][j] != syndrome_list[-t][j]
                else:
                    if t <= 1:
                        if t != self.T:
                            change = syndrome_list[-t - 1][j] != "0"
                        else:
                            change = syndrome_list[-t - 1][j] != syndrome_list[-t][j]
                    elif t == self.T:
                        last3 = ""
                        for dt in range(3):
                            last3 += syndrome_list[-t - 1 + dt][j]
                        change = last3.count("1") % 2 == 1
                    else:
                        change = syndrome_list[-t - 1][j] != syndrome_list[-t + 1][j]
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
                        if syn_type == 0:
                            elem_num = syn_round
                            syn_round = 0
                        node = {"time": syn_round}
                        is_boundary = syn_type == 0
                        if is_boundary:
                            i = [0, -1][elem_num]
                            if self.basis == "z":
                                qubits = [self.css_x_logical[i]]
                            else:
                                qubits = [self.css_z_logical[i]]
                        else:
                            if self.basis == "z":
                                qubits = self.css_z_gauge_ops[elem_num]
                            else:
                                qubits = self.css_x_gauge_ops[elem_num]
                        node["qubits"] = qubits
                        node["is_boundary"] = is_boundary
                        node["element"] = elem_num
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
