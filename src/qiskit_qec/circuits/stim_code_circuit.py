# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# pylint: disable=invalid-name, disable=no-name-in-module, disable=no-member

"""Generates CodeCircuits from stim circuits"""
import warnings

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

from qiskit.circuit.library.standard_gates import (
    IGate,
    XGate,
    YGate,
    ZGate,
    HGate,
    SGate,
    SdgGate,
    CXGate,
    CYGate,
    CZGate,
    SwapGate,
)

from stim import CircuitInstruction, CircuitRepeatBlock, target_inv
from stim import Circuit as StimCircuit

from qiskit_qec.circuits.code_circuit import CodeCircuit
from qiskit_qec.utils.stim_tools import (
    detector_error_model_to_rx_graph,
    string2nodes_with_detectors,
    string2logical_meas,
)


class StimCodeCircuit(CodeCircuit):
    """
    Prepares a CodeCircuit class based on the supplied stim circuit.
    """

    def __init__(
        self,
        stim_circuit: StimCircuit,
        barriers: bool = True,
    ):
        """
        Args:
            stim_circuit: stim circuit to be coverted
            barriers (optional): whether to include barriers (coeersponding to stim TICK instructions)
                in the qiskit circuits. Default is True
        Examples:
            Prepare and measure a Bell-state, checking the parity with a measurement
            comparison (DETECTOR)
            >>> from qiskit_qec.circuits.stim_code_circuit import StimCodeCircuit
            >>> stim_ex1 = stim.Circuit('''
            >>>     H 0
            >>>     TICK
            >>>     CX 0 1
            >>>     X_ERROR(0.2) 0 1
            >>>     TICK
            >>>     M 0 1
            >>>     DETECTOR rec[-1] rec[-2]
            >>> ''')
            >>> stim_code = StimCodeCircuit(stim_circuit = stim_ex1)
            >>> stim_code.qc
        """
        super().__init__()
        self.stim_circuit = stim_circuit

        self.measurement_data = []
        self.detectors = []
        self.logicals = []

        self.qc = QuantumCircuit()
        self.qc.add_register(QuantumRegister(self.stim_circuit.num_qubits))

        single_qubit_gate_dict = {
            "I": IGate(),
            "X": XGate(),
            "Y": YGate(),
            "Z": ZGate(),
            "H": HGate(),
            "S": SGate(),
            "S_DAG": SdgGate(),
        }

        two_qubit_gate_dict = {
            "CX": CXGate(),
            "CY": CYGate(),
            "CZ": CZGate(),
            "SWAP": SwapGate(),
        }

        def _helper(stim_circuit: StimCircuit, reps: int):
            nonlocal rep_block_count
            nonlocal block_count
            for rep_ind in range(reps):
                meas_count = 0
                for instruction in stim_circuit:
                    if isinstance(instruction, CircuitRepeatBlock):
                        _helper(instruction.body_copy(), instruction.repeat_count)
                        rep_block_count += 1
                    elif isinstance(instruction, CircuitInstruction):
                        inst_name = instruction.name
                        if inst_name in single_qubit_gate_dict:
                            qubits = [target.value for target in instruction.targets_copy()]
                            for q in qubits:
                                self.qc.append(single_qubit_gate_dict[inst_name], qargs=[q])
                        elif inst_name in two_qubit_gate_dict:
                            stim_targets = instruction.targets_copy()
                            for t1, t2 in zip(stim_targets[::2], stim_targets[1::2]):
                                if t1.is_qubit_target:
                                    q1, q2 = t1.value, t2.value
                                    self.qc.append(two_qubit_gate_dict[inst_name], qargs=[q1, q2])
                                elif t1.is_measurement_record_target:
                                    c1, q2 = t1.value, t2.value
                                    inst_name = inst_name[1:]  # remove the C from CX,CY,CZ
                                    self.qc.append(
                                        single_qubit_gate_dict[inst_name], qargs=[q2]
                                    ).c_if(c1, 1)

                        elif inst_name == "M":
                            qubits = [target.value for target in instruction.targets_copy()]
                            invert_result = [
                                target.is_inverted_result_target
                                for target in instruction.targets_copy()
                            ]
                            if reps > 1:
                                cl_reg_name = (
                                    "rep_block_"
                                    + str(rep_block_count)
                                    + "_rep_"
                                    + str(rep_ind)
                                    + "_meas_block_"
                                    + str(meas_count)
                                )
                            else:
                                cl_reg_name = (
                                    "block_"
                                    + str(rep_block_count)
                                    + "_meas_block_"
                                    + str(meas_count)
                                )
                            creg = ClassicalRegister(len(qubits), name=cl_reg_name)
                            self.qc.add_register(creg)
                            flip_qubits = [q for q, inv in zip(qubits, invert_result) if inv]
                            if flip_qubits != []:
                                self.qc.x(flip_qubits)
                            self.qc.measure(qubits, creg)
                            for i, q in enumerate(qubits):
                                self.measurement_data.append((cl_reg_name, i))
                            if flip_qubits != []:
                                self.qc.x(flip_qubits)
                            meas_count += 1
                        elif inst_name == "R":
                            qubits = [target.value for target in instruction.targets_copy()]
                            self.qc.reset(qubits)

                        elif inst_name == "TICK" and barriers:
                            self.qc.barrier()

                        elif inst_name == "DETECTOR":
                            meas_targets = [t.value for t in instruction.targets_copy()]
                            self.detectors.append(
                                {
                                    "clbits": [self.measurement_data[ind] for ind in meas_targets],
                                    "time": [],
                                    "qubits": [],
                                    "stim_coords": instruction.gate_args_copy(),
                                }
                            )
                        elif inst_name == "OBSERVABLE_INCLUDE":
                            meas_targets = [t.value for t in instruction.targets_copy()]
                            self.logicals.append(
                                {"clbits": [self.measurement_data[ind] for ind in meas_targets]}
                            )

        rep_block_count = 0
        block_count = 0
        self.decomp_stim_circuit = self.decompose_stim_circuit(self.stim_circuit)
        _helper(self.decomp_stim_circuit, 1)

        self.circuit = {"": self.qc}
        self.base = ""

        # if a set of measurement comparisons is deterministically 1 in the absence of errors,
        # the set of syndromes is compared to that
        noiseless_measurements = self.decomp_stim_circuit.compile_sampler().sample(1)[0]
        clbit_dict = {
            (clbit._register.name, clbit._index): clind
            for clind, clbit in enumerate(self.qc.clbits)
        }
        detector_meas_indices = [
            [clbit_dict[clbit] for clbit in det["clbits"]] for det in self.detectors
        ]
        self.det_ref_values = [
            sum(noiseless_measurements[meas_list]) % 2 for meas_list in detector_meas_indices
        ]

        # further code parameters
        try:
            self.d = len(self.stim_circuit.shortest_graphlike_error())  # code distance
        except ValueError:
            self.d = 0
        self.n = stim_circuit.num_qubits
        # the number of rounds is not necessarily well-defined (Floquet codes etc.)

    def decompose_stim_circuit(self, stim_circuit):
        """
        Decompose gates in the stim circuit into Clifford gates that have an equivalent qiskit gate.
        Errors are not included.
        """
        decompose_1_dict = {
            # single-qubit gates
            "C_XYZ": [("S_DAG", [0]), ("H", [0])],
            "C_ZYX": [("H", [0]), ("S", [0])],
            "H_XY": [("X", [0]), ("S", [0])],
            "H_XZ": [("H", [0])],
            "H_YZ": [("H", [0]), ("S", [0]), ("H", [0]), ("Z", [0])],
            "SQRT_X": [("H", [0]), ("S", [0]), ("H", [0])],
            "SQRT_X_DAG": [("S", [0]), ("H", [0]), ("S", [0])],
            "SQRT_Y": [("Z", [0]), ("H", [0])],
            "SQRT_Y_DAG": [("H", [0]), ("Z", [0])],
            "SQRT_Z": [("S", [0])],
            "SQRT_Z_DAG": [("S_DAG", [0])],
        }

        decompose_2_dict = {
            # two-qubit gates
            "CNOT": [("CX", [0, 1])],
            "CXSWAP": [("CX", [0, 1]), ("SWAP", [0, 1])],
            "SWAPCX": [("SWAP", [0, 1]), ("CX", [0, 1])],
            "ISWAP": [
                ("H", [0]),
                ("SWAP", [0, 1]),
                ("CX", [0, 1]),
                ("S", [0]),
                ("H", [1]),
                ("S", [1]),
            ],
            "ISWAP_DAG": [
                ("S_DAG", [0]),
                ("S_DAG", [1]),
                ("H", [1]),
                ("CX", [0, 1]),
                ("SWAP", [0, 1]),
                ("H", [0]),
            ],
            "SQRT_XX": [
                ("H", [0]),
                ("CX", [0, 1]),
                ("S", [0]),
                ("H", [0]),
                ("H", [1]),
                ("S", [1]),
                ("H", [1]),
            ],
            "SQRT_XX_DAG": [
                ("H", [0]),
                ("CX", [0, 1]),
                ("S_DAG", [0]),
                ("H", [0]),
                ("H", [1]),
                ("S_DAG", [1]),
                ("H", [1]),
            ],
            "SQRT_YY": [
                ("S_DAG", [0]),
                ("H", [0]),
                ("S_DAG", [1]),
                ("CX", [0, 1]),
                ("S", [0]),
                ("H", [0]),
                ("S", [0]),
                ("H", [1]),
                ("S", [1]),
                ("H", [1]),
                ("S", [1]),
            ],
            "SQRT_YY_DAG": [
                ("S_DAG", [0]),
                ("H", [0]),
                ("S", [1]),
                ("CX", [0, 1]),
                ("S", [0]),
                ("H", [0]),
                ("S", [0]),
                ("H", [1]),
                ("S", [1]),
                ("H", [1]),
                ("S_DAG", [1]),
            ],
            "SQRT_ZZ": [("H", [1]), ("CX", [0, 1]), ("S", [0]), ("H", [1]), ("S", [1])],
            "SQRT_ZZ_DAG": [
                ("H", [1]),
                ("CX", [0, 1]),
                ("S_DAG", [0]),
                ("H", [1]),
                ("S_DAG", [1]),
            ],
            "XCX": [("H", [0]), ("CX", [0, 1]), ("H", [0])],
            "XCY": [("H", [0]), ("CY", [0, 1]), ("H", [0])],
            "XCZ": [("H", [0]), ("CZ", [0, 1]), ("S", [0])],
            "YCX": [("S_DAG", [0]), ("H", [0]), ("CX", [0, 1]), ("H", [0]), ("S", [0])],
            "YCY": [("S_DAG", [0]), ("H", [0]), ("CY", [0, 1]), ("H", [0]), ("S", [0])],
            "YCZ": [("S_DAG", [0]), ("H", [0]), ("CZ", [0, 1]), ("H", [0]), ("S", [0])],
            "ZCX": [("CX", [0, 1])],
            "ZCY": [("CY", [0, 1])],
            "ZCZ": [("CZ", [0, 1])],
        }
        decompose_m_dict = {
            # measurements
            "MX": [("H", [0]), ("M", [0]), ("H", [0])],
            "MY": [
                ("S_DAG", [0]),
                ("H", [0]),
                ("M", [0]),
                ("H", [0]),
                ("S", [0]),
            ],
            "MZ": [("M", [0])],
            "RX": [("H", [0]), ("R", [0]), ("H", [0])],
            "RY": [
                ("S_DAG", [0]),
                ("H", [0]),
                ("R", [0]),
                ("H", [0]),
                ("S", [0]),
            ],
            "RZ": [("R", [0])],
            "MRX": [("H", [0]), ("M", [0]), ("R", [0]), ("H", [0])],
            "MRY": [
                ("S_DAG", [0]),
                ("H", [0]),
                ("M", [0]),
                ("R", [0]),
                ("H", [0]),
                ("S", [0]),
            ],
            "MRZ": [("M", [0]), ("R", [0])],
            "MR": [("M", [0]), ("R", [0])],
            "MXX": [("CX", [0, 1]), ("H", [0]), ("M", [0]), ("H", [0]), ("CX", [0, 1])],
            "MYY": [
                ("S_DAG", [0]),
                ("S_DAG", [1]),
                ("CX", [0, 1]),
                ("H", [0]),
                ("M", [0]),
                ("H", [0]),
                ("CX", [0, 1]),
                ("S", [0]),
                ("S", [1]),
            ],
            "MZZ": [("CX", [0, 1]), ("M", [1]), ("CX", [0, 1])],
        }

        stim_error_list = [
            "CORRELATED_ERROR",
            "DEPOLARIZE1",
            "DEPOLARIZE2",
            "E",
            "ELSE_CORRELATED_ERROR",
            "HERALDED_ERASE",
            "HERALDED_PAULI_CHANNEL_1",
            "PAULI_CHANNEL_1",
            "PAULI_CHANNEL_2",
            "X_ERROR",
            "Y_ERROR",
            "Z_ERROR",
        ]

        decomposed_stim_circuit = StimCircuit()
        for instruction in stim_circuit:
            if isinstance(instruction, CircuitInstruction):
                if instruction.name in decompose_1_dict:
                    targets = instruction.targets_copy()
                    decomp_gate = StimCircuit()
                    basis_gates = decompose_1_dict[instruction.name]
                    for target in targets:
                        for gate, _ in basis_gates:
                            decomp_gate.append(gate, [target], [])
                    decomposed_stim_circuit += decomp_gate
                elif instruction.name in decompose_2_dict:
                    targets = instruction.targets_copy()
                    decomp_gate = StimCircuit()
                    basis_gates = decompose_2_dict[instruction.name]
                    for target_pair in zip(targets[::2], targets[1::2]):
                        for gate, qubit_ind in basis_gates:
                            qubits = [target_pair[qubit_ind[0]]]
                            if len(qubit_ind) == 2:
                                qubits.append(target_pair[qubit_ind[1]])
                            decomp_gate.append(gate, qubits, [])
                    decomposed_stim_circuit += decomp_gate
                elif instruction.name in decompose_m_dict:
                    arg = instruction.gate_args_copy()
                    target_list = instruction.targets_copy()
                    decomp_gate = StimCircuit()
                    basis_gates = decompose_m_dict[instruction.name]
                    if instruction.name in ("MXX", "MYY", "MZZ"):
                        if len(set(target_list)) < len(
                            target_list
                        ):  # avoiding overlaps that result from stim broadcating
                            target_list = [
                                target_list[i : i + 2] for i in range(0, len(target_list), 2)
                            ]
                        else:
                            target_list = [target_list]
                        for targets in target_list:
                            for gate, qubit_ind in basis_gates:
                                if len(qubit_ind) == 2:
                                    qubit_list = [
                                        target_pair[q_ind].value
                                        for target_pair in zip(targets[::2], targets[1::2])
                                        for q_ind in qubit_ind
                                    ]
                                else:
                                    qubit_list = [
                                        target_pair[qubit_ind[0]].value
                                        for target_pair in zip(targets[::2], targets[1::2])
                                    ]

                                if gate == "M":
                                    inv_list = [
                                        (
                                            target_pair[0].is_inverted_result_target
                                            + target_pair[1].is_inverted_result_target
                                        )
                                        % 2
                                        for target_pair in zip(targets[::2], targets[1::2])
                                    ]
                                    for i, inv in enumerate(inv_list):
                                        if inv:
                                            qubit_list[i] = target_inv(qubit_list[i])
                                    decomp_gate.append(gate, qubit_list, arg)
                                else:
                                    decomp_gate.append(gate, qubit_list, [])

                    else:
                        for gate, _ in basis_gates:
                            qubit_list = [target.value for target in target_list]
                            if gate == "M":
                                inv_list = [
                                    target.is_inverted_result_target for target in target_list
                                ]
                                for i, inv in enumerate(inv_list):
                                    if inv:
                                        qubit_list[i] = target_inv(qubit_list[i])
                                decomp_gate.append(gate, qubit_list, arg)
                            else:
                                decomp_gate.append(gate, qubit_list, [])

                    decomposed_stim_circuit += decomp_gate

                elif instruction.name == "MPP":
                    decomposed_stim_circuit += self.MPP_circuit(instruction)
                elif instruction.name == "MPAD":
                    # MPAD does affect the stim measurement sample output, but not the detector sample.
                    # This is due to the fact that detectors are sensitive to changes
                    # wrt the perfect output,
                    # R 0
                    # MPAD 1
                    # M 0
                    # DETECTOR rec[-1] rec[-2]
                    # The above example results in a deterministinc detector samples False, not True!
                    #
                    # MPAD can have an effect on measurements, when a CX is controlled on an MPAD bit...
                    # but what is the point of that anyway?
                    warnings.warn(
                        "The circuit contains MPAD instructions that are ignored in the conversion."
                        "This can affect the measurement outcomes, but not the detectors."
                    )
                    pass
                elif instruction.name in stim_error_list:
                    # do not include errors
                    pass
                else:
                    gate = instruction.name
                    arg = instruction.gate_args_copy()
                    targets = instruction.targets_copy()
                    decomposed_stim_circuit.append(gate, targets, arg)
            elif isinstance(instruction, CircuitRepeatBlock):
                decomposed_stim_circuit.append(
                    CircuitRepeatBlock(
                        instruction.repeat_count,
                        self.decompose_stim_circuit(instruction.body_copy()),
                    )
                )

        return decomposed_stim_circuit

    def MPP_circuit(self, stim_instruction):
        """Handle MPP measurements."""
        MPP_stim_circuit = StimCircuit()
        arg = stim_instruction.gate_args_copy()

        # break it down into individual Pauli products
        target_lists = []
        target_list = []
        prev_is_combiner = True
        for target in stim_instruction.targets_copy():
            if prev_is_combiner:
                target_list.append(target)
                prev_is_combiner = False
            elif target.is_combiner:
                prev_is_combiner = True
            else:
                target_lists.append(target_list)
                target_list = [target]
                prev_is_combiner = False
        target_lists.append(target_list)

        for target_list in target_lists:
            invert = 0
            first_target_qubit = target_list[0].value
            for target in target_list:
                if target.is_x_target:
                    MPP_stim_circuit.append("H", target.value)
                elif target.is_y_target:
                    MPP_stim_circuit.append("S_DAG", target.value)
                    MPP_stim_circuit.append("H", target.value)
                invert = (invert + target.is_inverted_result_target) % 2
                if target.value != first_target_qubit:
                    MPP_stim_circuit.append("CX", [target.value, first_target_qubit])
            if invert:
                MPP_stim_circuit.append("M", target_inv(first_target_qubit), arg)
            else:
                MPP_stim_circuit.append("M", first_target_qubit, arg)

            for target in target_list[::-1]:
                if target.value != first_target_qubit:
                    MPP_stim_circuit.append("CX", [target.value, first_target_qubit])
                if target.is_x_target:
                    MPP_stim_circuit.append("H", target.value)
                elif target.is_y_target:
                    MPP_stim_circuit.append("H", target.value)
                    MPP_stim_circuit.append("S", target.value)

        return MPP_stim_circuit

    def string2nodes(self, string, **kwargs):
        """
        Convert output string from circuits into a set of nodes for `DecodingGraph`.
        Args:
            string (string): Results string to convert.
            kwargs (dict): Any additional keyword arguments.
                logical (str): Logical value whose results are used ('0' as default).
                all_logicals (bool): Whether to include logical nodes
                irrespective of value. (False as default).
        """

        nodes = string2nodes_with_detectors(
            string=string,
            detectors=self.detectors,
            logicals=self.logicals,
            clbits=self.qc.clbits,
            det_ref_values=self.det_ref_values,
            **kwargs,
        )
        return nodes

    def string2raw_logicals(self, string):
        """
        Converts output string into a list of logical measurement outcomes
        Logicals are the logical measurements produced by self.stim_detectors()
        """
        _, self.logicals = self.stim_detectors()

        log_outs = string2logical_meas(string, self.logicals, self.circuit.clbits)

        return log_outs

    def _make_syndrome_graph(self):
        e = self.stim_circuit.detector_error_model(
            decompose_errors=True, approximate_disjoint_errors=True
        )
        graph, hyperedges = detector_error_model_to_rx_graph(e, detectors=self.detectors)
        return graph, hyperedges

    def check_nodes(self, nodes, ignore_extras=False, minimal=False):
        raise NotImplementedError
