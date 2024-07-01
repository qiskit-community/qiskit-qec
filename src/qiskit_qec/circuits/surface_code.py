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

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

from qiskit_qec.utils import DecodingGraphNode
from qiskit_qec.circuits.code_circuit import CodeCircuit


class SurfaceCodeCircuit(CodeCircuit):

    """Distance d rotated surface code with  T syndrome measurement rounds."""

    def __init__(self, d: int, T: int, basis: str = "z", resets=True):
        """Creates the circuits corresponding to logical basis states.

        Creates the circuits corresponding to logical basis states encoded
        using a rotated surface code.

        Args:
            d (int): Number of code qubits (and hence repetitions) used.
            T (int): Number of rounds of ancilla-assisted syndrome measurement.
            basis (str): Basis used to initialize qubit.
            resets (bool): Whether to include a reset gate after mid-circuit measurements.


        Additional information:
            No measurements are added to the circuit if `T=0`. Otherwise
            `T` rounds are added, followed by measurement of the code
            qubits (corresponding to a logical measurement and final
            syndrome measurement round).
        """
        super().__init__()

        self.d = d
        self.n = d**2
        self.T = 0
        self.basis = basis
        self._resets = resets

        # get layout of plaquettes
        self.zplaqs, self.xplaqs, self._zplaq_coords, self._xplaq_coords = self._get_plaquettes()

        self._logicals = {"x": [], "z": []}
        # X logicals for left and right sides
        self._logicals["x"].append([j * self.d for j in range(self.d)])
        self._logicals["x"].append([(j + 1) * self.d - 1 for j in range(self.d)])
        # Z logicals for top and bottom rows
        self._logicals["z"].append(list(range(self.d)))
        self._logicals["z"].append([self.d**2 - 1 - j for j in range(self.d)])

        # set gauge and stabilizer info
        self.x_gauge_ops = [[q for q in plaq if q is not None] for plaq in self.xplaqs]
        self.x_stabilizer_ops = self.x_gauge_ops
        self.x_logical = [self._logicals["x"][0]]
        self.x_boundary = [self._logicals["x"][0] + self._logicals["x"][1]]
        self.z_gauge_ops = [[q for q in plaq if q is not None] for plaq in self.zplaqs]
        self.z_stabilizer_ops = self.z_gauge_ops
        self.z_logical = [self._logicals["z"][0]]
        self.z_boundary = [self._logicals["z"][0] + self._logicals["z"][1]]

        # quantum registers
        self._num_xy = int((d**2 - 1) / 2)
        self.code_qubit = QuantumRegister(d**2, "code_qubit")
        self.zplaq_qubit = QuantumRegister(self._num_xy, "zplaq_qubit")
        self.xplaq_qubit = QuantumRegister(self._num_xy, "xplaq_qubit")
        self.qubit_registers = [self.code_qubit, self.zplaq_qubit, self.xplaq_qubit]

        # classical registers
        self.xplaq_bits = []
        self.zplaq_bits = []
        self.code_bit = ClassicalRegister(d**2, "code_bit")

        # create the circuits
        self.circuit = {}
        for log in ["0", "1"]:
            self.circuit[log] = QuantumCircuit(
                self.code_qubit, self.zplaq_qubit, self.xplaq_qubit, name=log
            )
        self.base = "0"

        # apply initial logical paulis for encoded states
        self._preparation()

        # add the gates required for syndrome measurements
        for _ in range(T - 1):
            self.syndrome_measurement()
        if T != 0:
            self.syndrome_measurement(final=True)
            self.readout()

    def _get_plaquettes(self):
        """
        Returns `zplaqs` and `xplaqs`, which are lists of the Z and X type
        stabilizers. Each plaquettes is specified as a list of four qubits,
        in the order in which entangling gates are applied.
        """

        d = self.d

        zplaqs = []
        xplaqs = []
        zplaq_coords = []
        xplaq_coords = []
        for y in range(-1, d):
            for x in range(-1, d):
                bulk = x in range(d - 1) and y in range(d - 1)
                ztab = (x == -1 and y % 2 == 0) or (x == d - 1 and y % 2 == 1)
                xtab = (y == -1 and x % 2 == 1) or (y == d - 1 and x % 2 == 0)

                if (x in range(d - 1) or y in range(d - 1)) and (bulk or ztab or xtab):
                    plaq = []
                    for dy in range(2):
                        for dx in range(2):
                            if x + dx in range(d) and y + dy in range(d):
                                plaq.append(x + dx + d * (y + dy))
                            else:
                                plaq.append(None)

                    if (x + y) % 2 == 0:
                        xplaqs.append([plaq[0], plaq[1], plaq[2], plaq[3]])
                        xplaq_coords.append((x, y))
                    else:
                        zplaqs.append([plaq[0], plaq[2], plaq[1], plaq[3]])
                        zplaq_coords.append((x, y))

        return zplaqs, xplaqs, zplaq_coords, xplaq_coords

    def _preparation(self):
        """
        Prepares logical bit states by applying an x to the circuit that will
        encode a 1.
        """
        if self.basis == "z":
            self.x(["1"])
        else:
            for qc in self.circuit.values():
                qc.h(self.code_qubit)
            self.z(["1"])

    def get_circuit_list(self):
        """
        Returns:
            circuit_list: self.circuit as a list, with
            circuit_list[0] = circuit['0']
            circuit_list[1] = circuit['1']
        """
        circuit_list = [self.circuit[log] for j, log in enumerate(["0", "1"])]
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
            for j in range(self.d):
                self.circuit[log].x(self.code_qubit[j * self.d])
            if barrier:
                self.circuit[log].barrier()

    def z(self, logs=("0", "1"), barrier=False):
        """Applies a logical z to the circuits for the given logical values.

        Args:
            logs (list or tuple): List or tuple of logical values expressed as
                strings.
            barrier (bool): Boolean denoting whether to include a barrier at
                the end.
        """
        for log in logs:
            for j in range(self.d):
                self.circuit[log].z(self.code_qubit[j])
            if barrier:
                self.circuit[log].barrier()

    def syndrome_measurement(self, final=False, barrier=False):
        """Application of a syndrome measurement round.

        Args:
            final (bool): Whether to disregard the reset (if applicable) due to this
                being the final syndrome measurement round.
            barrier (bool): Boolean denoting whether to include a barrier at the end.

        """

        zplaqs, xplaqs = self.zplaqs, self.xplaqs

        # classical registers for this round
        self.zplaq_bits.append(
            ClassicalRegister(self._num_xy, "round_" + str(self.T) + "_zplaq_bit")
        )
        self.xplaq_bits.append(
            ClassicalRegister(self._num_xy, "round_" + str(self.T) + "_xplaq_bit")
        )

        for log in ["0", "1"]:
            self.circuit[log].add_register(self.zplaq_bits[-1])
            self.circuit[log].add_register(self.xplaq_bits[-1])

            self.circuit[log].h(self.xplaq_qubit)

            for j in range(4):
                for p, plaq in enumerate(zplaqs):
                    c = plaq[j]
                    if c is not None:
                        self.circuit[log].cx(self.code_qubit[c], self.zplaq_qubit[p])
                for p, plaq in enumerate(xplaqs):
                    c = plaq[j]
                    if c is not None:
                        self.circuit[log].cx(self.xplaq_qubit[p], self.code_qubit[c])

            self.circuit[log].h(self.xplaq_qubit)

            for j in range(self._num_xy):
                self.circuit[log].measure(self.xplaq_qubit[j], self.xplaq_bits[self.T][j])
                self.circuit[log].measure(self.zplaq_qubit[j], self.zplaq_bits[self.T][j])
                if self._resets and not final:
                    self.circuit[log].reset(self.xplaq_qubit[j])
                    self.circuit[log].reset(self.zplaq_qubit[j])

            if barrier:
                self.circuit[log].barrier()

        self.T += 1

    def readout(self):
        """
        Readout of all code qubits, which corresponds to a logical measurement
        as well as allowing for a measurement of the syndrome to be inferred.
        """

        for log in ["0", "1"]:
            if self.basis == "x":
                self.circuit[log].h(self.code_qubit)
            self.circuit[log].add_register(self.code_bit)
            self.circuit[log].measure(self.code_qubit, self.code_bit)

    def _string2changes(self, string):
        basis = self.basis

        # final syndrome for plaquettes deduced from final code qubit readout
        final_readout = string.split(" ")[0][::-1]
        if basis == "z":
            plaqs = self.zplaqs
        else:
            plaqs = self.xplaqs
        full_syndrome = ""
        for plaq in plaqs:
            parity = 0
            for q in plaq:
                if q is not None:
                    parity += int(final_readout[q])
            full_syndrome = str(parity % 2) + full_syndrome

        # results from all other plaquette syndrome measurements then added
        if basis == "z":
            full_syndrome = full_syndrome + " " + " ".join(string.split(" ")[2::2])
        else:
            full_syndrome = full_syndrome + " " + " ".join(string.split(" ")[1::2])

        # changes between one syndrome and the next then calculated
        syndrome_list = full_syndrome.split(" ")
        height = len(syndrome_list)
        width = len(syndrome_list[0])
        syndrome_changes = ""
        for t in range(height):
            for j in range(width):
                if self._resets:
                    if t == 0:
                        change = syndrome_list[-1][j] != "0"
                    else:
                        change = syndrome_list[-t][j] != syndrome_list[-t - 1][j]
                    syndrome_changes += "0" * (not change) + "1" * change
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
        syndrome_changes = syndrome_changes[0:-1]

        if basis != self.basis:
            # trim the noisy nonsense (first and last rounds)
            syndrome_changes = " ".join(syndrome_changes.split(" ")[1:-1])

        return syndrome_changes

    def string2raw_logicals(self, string):
        """
        Extracts raw logicals from output string.
        Args:
            string (string): Results string from which to extract logicals
        Returns:
            list: Raw values for logical operators that correspond to nodes.
        """
        final_readout = string.split(" ")[0][::-1]
        # get logical readout
        # (though it's called Z, it actually depends on the basis)
        Z = [0, 0]
        for j in range(self.d):
            if self.basis == "z":
                # evaluated using top row
                Z[0] = (Z[0] + int(final_readout[j])) % 2
                # evaluated using bottom row
                Z[1] = (Z[1] + int(final_readout[self.d**2 - 1 - j])) % 2
            else:
                # evaluated using left side
                Z[0] = (Z[0] + int(final_readout[j * self.d])) % 2
                # evaluated using right side
                Z[1] = (Z[1] + int(final_readout[(j + 1) * self.d - 1])) % 2
        return [str(Z[0]), str(Z[1])]

    def _process_string(self, string):
        # get logical readout
        measured_Z = self.string2raw_logicals(string)

        # then get syndrome changes
        syndrome_changes = self._string2changes(string)

        # the space separated string of syndrome changes then gets a
        # double space separated logical value on the end
        new_string = " ".join(measured_Z) + "  " + syndrome_changes

        return new_string

    def _separate_string(self, string):
        separated_string = []
        for syndrome_type_string in string.split("  "):
            separated_string.append(syndrome_type_string.split(" "))
        return separated_string

    def string2nodes(self, string, **kwargs):
        """Convert output string from circuits into a set of nodes.

        Args:
            string (string): Results string to convert.
            kwargs (dict): Additional keyword arguments. See below:

        kwargs:
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

        all_logicals = kwargs.get("all_logicals")
        logical = kwargs.get("logical")
        if logical is None:
            logical = "0"

        string = self._process_string(string)
        separated_string = self._separate_string(string)
        nodes = []

        # boundary nodes
        boundary = separated_string[0]  # [<last_elem>, <init_elem>]
        for bqec_index, belement in enumerate(boundary[::-1]):
            if all_logicals or belement != logical:
                node = DecodingGraphNode(
                    is_logical=True,
                    is_boundary=True,
                    qubits=self._logicals[self.basis][-bqec_index - 1],
                    index=1 - bqec_index,
                )
                nodes.append(node)

        # bulk nodes
        for syn_type in range(1, len(separated_string)):
            for syn_round in range(len(separated_string[syn_type])):
                elements = separated_string[syn_type][syn_round]
                for qec_index, element in enumerate(elements[::-1]):
                    if element == "1":
                        if self.basis == "x":
                            qubits = self.x_stabilizer_ops[qec_index]
                        else:
                            qubits = self.z_stabilizer_ops[qec_index]
                        node = DecodingGraphNode(time=syn_round, qubits=qubits, index=qec_index)
                        nodes.append(node)
        return nodes

    def check_nodes(self, nodes, ignore_extras=False, minimal=False):
        """
        Determines whether a given set of nodes are neutral. If so, also
        determines any additional logical readout qubits that would be
        flipped by the errors creating such a cluster and how many errors
        would be required to make the cluster.
        Args:
            nodes (list): List of nodes, of the type produced by `string2nodes`.
            ignore_extras (bool): If `True`, undeeded logical nodes are
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

        bulk_nodes = [node for node in nodes if not node.is_logical]
        logical_nodes = [node for node in nodes if node.is_logical]
        given_logicals = set(node.index for node in logical_nodes)

        if self.basis == "z":
            coords = self._zplaq_coords
        else:
            coords = self._xplaq_coords

        if (len(bulk_nodes) % 2) == 0:
            if (len(logical_nodes) % 2) == 0 or ignore_extras:
                neutral = True
                flipped_logicals = set()
                # estimate num_errors from size
                if bulk_nodes:
                    xs = []
                    ys = []
                    for node in bulk_nodes:
                        x, y = coords[node.index]
                        xs.append(x)
                        ys.append(y)
                    dx = max(xs) - min(xs)
                    dy = max(ys) - min(ys)
                    num_errors = dx + dy
                    if dx > 0 and dy > 0:
                        num_errors -= 1
                else:
                    num_errors = 0
            else:
                neutral = False
                flipped_logicals = set()
                num_errors = 0
        else:
            # find nearest boundary
            num_errors = (self.d - 1) / 2
            for node in bulk_nodes:
                x, y = coords[node.index]
                if self.basis == "z":
                    p = y
                else:
                    p = x
                num_errors = min(num_errors, p + 1, self.d - p)
            flipped_logicals = {1 - int(p < (self.d - 1) / 2)}

        # if unneeded logical zs are given, cluster is not neutral
        # (unless this is ignored)
        if (not ignore_extras) and given_logicals.difference(flipped_logicals):
            neutral = False
        # otherwise, report only needed logicals that aren't given
        else:
            neutral = True
            flipped_logicals = flipped_logicals.difference(given_logicals)

        # get the required boundary nodes
        flipped_logical_nodes = []
        for elem in flipped_logicals:
            node = DecodingGraphNode(
                is_logical=True,
                is_boundary=True,
                qubits=self._logicals[self.basis][elem],
                index=elem,
            )
            flipped_logical_nodes.append(node)

        return neutral, flipped_logical_nodes, num_errors

    def is_cluster_neutral(self, nodes):
        """
        Determines whether or not the cluster is neutral, meaning that one or more
        errors could have caused the set of nodes (syndrome changes) passed
        to the method.
        Args:
            nodes (dictionary in the form of the return value of string2nodes)
        """
        return not bool(len(nodes) % 2)
