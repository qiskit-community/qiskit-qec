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

# pylint: disable=invalid-name, disable=no-name-in-module

"""Generates circuits for CSS codes."""
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer.noise import depolarizing_error, pauli_error

from qiskit_qec.circuits.code_circuit import CodeCircuit
from qiskit_qec.utils.stim_tools import (
    noisify_circuit,
    get_stim_circuits,
    detector_error_model_to_rx_graph,
    string2nodes_with_detectors,
    string2logical_meas,
)
from qiskit_qec.codes import StabSubSystemCode
from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.linear.symplectic import normalizer
from qiskit_qec.exceptions import QiskitQECError


class CSSCodeCircuit(CodeCircuit):
    """
    CodeCircuit class for generic CSS codes.
    """

    def __init__(
        self,
        code,
        T: int,
        basis: str = "z",
        round_schedule: str = "zx",
        noise_model=None,
    ):
        """
        Args:
            code: A CSS code class which is either
                a) StabSubSystemCode
                b) a class with the following methods:
                    'x_gauges' (as a list of list of qubit indices),
                    'z_gauges',
                    'x_stabilizers',
                    'z_stabilizers',
                    'logical_x',
                    'logical_z',
                    'n' (number of qubits),
            T: Number of syndrome measurement rounds
            basis: basis for encoding ('x' or 'z')
            round_schedule: Order in which to measureme gauge operators ('zx' or 'xz')
            noise_model: Pauli noise model used in the construction of noisy circuits.
            If a tuple, a pnenomological noise model is used with the entries being
            probabity of depolarizing noise on code qubits between rounds and
            probability of measurement errors, respectively.
        Examples:
            The QuantumCircuit of a memory experiment for the distance-3 HeavyHEX code
            >>> from qiskit_qec.codes.hhc import HHC
            >>> from qiskit_qec.circuits.css_code import CSSCodeCircuit
            >>> code = CSSCodeCircuit(HHC(3),T=3,basis='x',noise_model=(0.01,0.01),round_schedule='xz')
            >>> code.circuit['0']
        """

        super().__init__()

        self.code = code
        self._get_code_properties()
        self.T = T
        self.basis = basis
        self.base = "0"
        self.round_schedule = round_schedule
        self.noise_model = noise_model
        self._phenom = isinstance(noise_model, tuple)
        if self._phenom:
            p_depol, p_meas = self.noise_model
            self._depol_error = depolarizing_error(p_depol, 1)
            self._meas_error = pauli_error([("X", p_meas), ("I", 1 - p_meas)])

        circuit = {}
        states = ["0", "1"]
        if self._phenom:
            states += ["0n", "1n"]
        for state in states:
            qc = QuantumCircuit()
            qregs = []
            qregs.append(QuantumRegister(code.n, name="code qubits"))
            qregs.append(QuantumRegister(len(self.z_gauges), name="z auxs"))
            qregs.append(QuantumRegister(len(self.x_gauges), name="x auxs"))
            for qreg in qregs:
                qc.add_register(qreg)
            self._prepare_initial_state(qc, qregs, state)
            self._perform_syndrome_measurements(qc, qregs, state)
            creg = ClassicalRegister(code.n, name="final_readout")
            qc.add_register(creg)
            self._final_readout(qc, qregs, creg, state)
            circuit[state] = qc

        self.circuit = {}
        self.noisy_circuit = {}
        for state, qc in circuit.items():
            if state[-1] == "n" and self._phenom:
                self.noisy_circuit[state[0]] = qc
            else:
                self.circuit[state] = qc
        if noise_model and not self._phenom:
            for state, qc in circuit.items():
                self.noisy_circuit[state] = noisify_circuit(qc, noise_model)

        self._gauges4stabilizers = []
        self._stabilizers = [self.x_stabilizers, self.z_stabilizers]
        self._gauges = [self.x_gauges, self.z_gauges]
        for j in range(2):
            self._gauges4stabilizers.append([])
            for stabilizer in self._stabilizers[j]:
                gauges = []
                for g, gauge in enumerate(self._gauges[j]):
                    if set(stabilizer).intersection(set(gauge)) == set(gauge):
                        gauges.append(g)
                self._gauges4stabilizers[j].append(gauges)
        self.detectors, self.logicals = self.stim_detectors()

    def _get_code_properties(self):
        if isinstance(self.code, StabSubSystemCode):
            is_css = True

            raw_gauges = self.code.gauge_group.generators
            center, log, conj_log = normalizer(self.code.generators.matrix)
            raw_stabilizers = PauliList(center)
            raw_logicals = PauliList(log) + PauliList(conj_log)

            gauges = [[], []]
            stabilizers = [[], []]
            logicals = [[], []]

            for (raw_ops, ops,) in zip(
                [raw_gauges, raw_stabilizers, raw_logicals],
                [gauges, stabilizers, logicals],
            ):
                for op in raw_ops:
                    op = str(op)
                    for j, pauli in enumerate(["X", "Z"]):
                        if (op.count(pauli) + op.count("I")) == self.code.n:
                            ops[j].append([k for k, p in enumerate(op[::-1]) if p == pauli])
                is_css = is_css and (len(ops[0]) + len(ops[1])) == len(raw_ops)

            # extra stabilizers: the product of all others
            for j in range(2):
                combined = []
                for stabilizer in stabilizers[j]:
                    combined += stabilizer
                stabilizers[j].append([])
                for q in combined:
                    if combined.count(q) % 2:
                        stabilizers[j][-1].append(q)

            if is_css:
                self.x_gauges = gauges[0]
                self.z_gauges = gauges[1]
                self.x_stabilizers = stabilizers[0]
                self.z_stabilizers = stabilizers[1]
                self.logical_x = logicals[0]
                self.logical_z = logicals[1]
            else:
                raise QiskitQECError("Code is not obviously CSS.")

        else:
            # otherwise assume it has the info
            self.x_gauges = self.code.x_gauges
            self.z_gauges = self.code.z_gauges
            self.x_stabilizers = self.code.x_stabilizers
            self.z_stabilizers = self.code.z_stabilizers
            self.logical_x = self.code.logical_x
            self.logical_z = self.code.logical_z
        # for the unionfind decoder
        self.css_x_logical = self.logical_x
        self.css_z_logical = self.logical_z

    def measured_logicals(self):
        if self.basis == "x":
            measured_logicals = self.logical_x
        else:
            measured_logicals = self.logical_z
        return measured_logicals

    def _prepare_initial_state(self, qc, qregs, state):
        if state[0] == "1":
            if self.basis == "z":
                qc.x(self.logical_x[0])
            else:
                qc.x(self.logical_z[0])
        if self.basis == "x":
            qc.h(qregs[0])

    def _perform_syndrome_measurements(self, qc, qregs, state):
        for t in range(self.T):
            if state[-1] == "n" and self._phenom:
                for q in qregs[0]:
                    qc.append(self._depol_error, [q])
            # gauge measurements
            if self.round_schedule == "zx":
                self._z_gauge_measurements(qc, t, state)
                self._x_gauge_measurements(qc, t, state)
            elif self.round_schedule == "xz":
                self._x_gauge_measurements(qc, t, state)
                self._z_gauge_measurements(qc, t, state)
            else:
                raise NotImplementedError(
                    "Round schedule " + self.round_schedule + " not supported."
                )

    def _final_readout(self, qc, qregs, creg, state):
        if self.basis == "x":
            qc.h(qregs[0])
        if state[-1] == "n" and self._phenom:
            for q in qregs[0]:
                qc.append(self._meas_error, [q])
        qc.measure(qregs[0], creg)

    def _z_gauge_measurements(self, qc, t, state):
        creg = ClassicalRegister(len(self.z_gauges), name="round_" + str(t) + "_z_bits")
        qc.add_register(creg)
        for g, z_gauge in enumerate(self.z_gauges):
            for q in z_gauge:
                qc.cx(qc.qregs[0][q], qc.qregs[1][g])
            if state[-1] == "n" and self._phenom:
                qc.append(self._meas_error, [qc.qregs[1][g]])
            qc.measure(qc.qregs[1][g], creg[g])
            qc.reset(qc.qregs[1][g])

    def _x_gauge_measurements(self, qc, t, state):
        creg = ClassicalRegister(len(self.x_gauges), name="round_" + str(t) + "_x_bits")
        qc.add_register(creg)
        for g, x_gauge in enumerate(self.x_gauges):
            for q in x_gauge:
                qc.h(qc.qregs[0][q])
                qc.cx(qc.qregs[0][q], qc.qregs[2][g])
                qc.h(qc.qregs[0][q])
            if state[-1] == "n" and self._phenom:
                qc.append(self._meas_error, [qc.qregs[2][g]])
            qc.measure(qc.qregs[2][g], creg[g])
            qc.reset(qc.qregs[2][g])

    def string2nodes(self, string, **kwargs):
        """
        Convert output string from circuits into a set of nodes for
        `DecodingGraph`.
        Args:
            string (string): Results string to convert.
            kwargs (dict): Any additional keyword arguments.
                logical (str): Logical value whose results are used ('0' as default).
                all_logicals (bool): Whether to include logical nodes
                irrespective of value. (False as default).
        Returns:
            a list of 'DecodingGraphNode()'-s corresponding to the triggered detectors
        """

        nodes = string2nodes_with_detectors(
            string=string,
            detectors=self.detectors,
            logicals=self.logicals,
            clbits=self.circuit["0"].clbits,
            det_ref_values=0,
            **kwargs,
        )
        return nodes

    def string2raw_logicals(self, string):
        """
        Converts output string into a list of logical measurement outcomes
        Logicals are the logical measurements produced by self.stim_detectors()
        """
        _, self.logicals = self.stim_detectors()

        log_outs = string2logical_meas(string, self.logicals, self.circuit["0"].clbits)
        return log_outs

    def check_nodes(self, nodes, ignore_extra_boundary=False, minimal=False):
        raise NotImplementedError

    def is_cluster_neutral(self, atypical_nodes):
        raise NotImplementedError

    def stim_detectors(self):
        """
        Returns:
            detectors (list[dict]) are dictionaries containing
                a) 'clbits', the classical bits (register, index) included in the measurement comparisons
                b) 'qubits', the qubits (list of indices) participating in the stabilizer measurements
                c) 'time', measurement round (int) of the earlier measurements in the detector
                d) 'basis', the pauli basis ('x' or 'z') of the stabilizers
            logicals (list[dict]) are dictionaries containing
                a) 'clbits', the classical bits (register, index) included in the logical measurement
                b) 'basis', the pauli basis ('x' or 'z') of the logical
        """

        detectors = []
        logicals = []

        ## 0th round of measurements
        if self.basis == "x":
            reg = "round_0_x_bits"
            for stabind, stabilizer in enumerate(self.x_stabilizers):
                det = {"clbits": []}
                for gauge_ind in self._gauges4stabilizers[0][stabind]:
                    det["clbits"].append((reg, gauge_ind))
                det["qubits"] = stabilizer.copy()
                det["time"] = 0
                det["basis"] = "x"
                detectors.append(det)

        else:
            reg = "round_0_z_bits"
            for stabind, stabilizer in enumerate(self.z_stabilizers):
                det = {"clbits": []}
                for gauge_ind in self._gauges4stabilizers[1][stabind]:
                    det["clbits"].append((reg, gauge_ind))
                det["qubits"] = stabilizer.copy()
                det["time"] = 0
                det["basis"] = "z"
                detectors.append(det)

        # adding first x and then z stabilizer comparisons
        for j, basis in enumerate(["x", "z"]):
            for t in range(
                1, self.T
            ):  # compare stabilizer measurements with previous in each round
                reg_prev = "round_" + str(t - 1) + "_" + basis + "_bits"
                reg_t = "round_" + str(t) + "_" + basis + "_bits"
                for gind, gs in enumerate(self._gauges4stabilizers[j]):
                    det = {"clbits": []}
                    for gauge_ind in gs:
                        det["clbits"].append((reg_t, gauge_ind))
                        det["clbits"].append((reg_prev, gauge_ind))
                    det["qubits"] = self._stabilizers[j][gind].copy()
                    det["time"] = t
                    det["basis"] = basis
                    detectors.append(det)

        ## final measurements
        if self.basis == "x":
            reg_prev = "round_" + str(self.T - 1) + "_x_bits"
            reg_T = "final_readout"
            for stabind, stabilizer in enumerate(self.x_stabilizers):
                det = {"clbits": []}
                for q in stabilizer:
                    det["clbits"].append((reg_T, q))
                for gauge_ind in self._gauges4stabilizers[0][stabind]:
                    det["clbits"].append((reg_prev, gauge_ind))
                det["qubits"] = stabilizer.copy()
                det["time"] = self.T
                det["basis"] = "x"
                detectors.append(det)
            logicals.append(
                {
                    "clbits": [(reg_T, q) for q in sorted(self.logical_x[0])],
                    "basis": "z",
                }
            )
        else:
            reg_prev = "round_" + str(self.T - 1) + "_z_bits"
            reg_T = "final_readout"
            for stabind, stabilizer in enumerate(self.z_stabilizers):
                det = {"clbits": []}
                for q in stabilizer:
                    det["clbits"].append((reg_T, q))
                for gauge_ind in self._gauges4stabilizers[1][stabind]:
                    det["clbits"].append((reg_prev, gauge_ind))
                det["qubits"] = stabilizer.copy()
                det["time"] = self.T
                det["basis"] = "x"
                detectors.append(det)
            logicals.append(
                {
                    "clbits": [(reg_T, q) for q in sorted(self.logical_z[0])],
                    "basis": "z",
                }
            )

        return detectors, logicals

    def _make_syndrome_graph(self):
        """
        Used by the DecodingGraph class to build the decoding graph and the obtain hyperedges
        """
        detectors, logicals = self.stim_detectors()
        stim_circuit = get_stim_circuits(
            self.noisy_circuit["0"], detectors=detectors, logicals=logicals
        )[0][0]
        e = stim_circuit.detector_error_model(
            decompose_errors=True, approximate_disjoint_errors=True
        )
        graph, hyperedges = detector_error_model_to_rx_graph(e, detectors=detectors)
        return graph, hyperedges
