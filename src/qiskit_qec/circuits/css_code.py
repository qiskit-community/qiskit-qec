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

from stim import Circuit as StimCircuit
from stim import target_rec as StimTarget_rec

from qiskit_qec.utils import DecodingGraphNode
from qiskit_qec.circuits.code_circuit import CodeCircuit
from qiskit_qec.utils.stim_tools import (
    noisify_circuit,
    get_stim_circuits,
    detector_error_model_to_rx_graph,
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
        self, code, T: int, basis: str = "z", round_schedule: str = "zx", noise_model=None
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

            for (
                raw_ops,
                ops,
            ) in zip([raw_gauges, raw_stabilizers, raw_logicals], [gauges, stabilizers, logicals]):
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
        """

        all_logicals = kwargs.get("all_logicals")
        logical = kwargs.get("logical")
        if logical is None:
            logical = "0"

        output = string.split(" ")[::-1]
        gauge_outs = [[], []]
        for t in range(self.T):
            gauge_outs[0].append(
                [int(b) for b in output[2 * t + self.round_schedule.find("x")]][::-1]
            )
            gauge_outs[1].append(
                [int(b) for b in output[2 * t + self.round_schedule.find("z")]][::-1]
            )
        final_outs = [int(b) for b in output[-1]]

        stabilizer_outs = []
        for j in range(2):
            stabilizer_outs.append([])
            for t in range(self.T):
                round_outs = []
                for gs in self._gauges4stabilizers[j]:
                    out = 0
                    for g in gs:
                        out += gauge_outs[j][t][g]
                    out = out % 2
                    round_outs.append(out)
                stabilizer_outs[j].append(round_outs)

        bases = ["x", "z"]
        j = bases.index(self.basis)
        final_gauges = []
        for gauge in self._gauges[j]:
            out = 0
            for q in gauge:
                out += final_outs[-q - 1]
            out = out % 2
            final_gauges.append(out)
        final_stabilizers = []
        for gs in self._gauges4stabilizers[j]:
            out = 0
            for g in gs:
                out += final_gauges[g]
            out = out % 2
            final_stabilizers.append(out)
        stabilizer_outs[j].append(final_stabilizers)

        stabilizer_changes = []
        for j in range(2):
            stabilizer_changes.append([])
            for t in range(self.T + (bases[j] == self.basis)):
                stabilizer_changes[j].append([])
                for e in range(len(stabilizer_outs[j][t])):
                    if t == 0 and j == bases.index(self.basis):
                        stabilizer_changes[j][t].append(stabilizer_outs[j][t][e])
                    else:
                        stabilizer_changes[j][t].append(
                            (stabilizer_outs[j][t][e] + stabilizer_outs[j][t - 1][e]) % 2
                        )

        nodes = []
        for j in range(2):
            for t, round_changes in enumerate(stabilizer_changes[j]):
                for e, change in enumerate(round_changes):
                    if change == 1:
                        node = DecodingGraphNode(time=t, qubits=self._stabilizers[j][e], index=e)
                        node.properties["basis"] = bases[j]
                        nodes.append(node)

        if self.basis == "x":
            logicals = self.logical_x
        else:
            logicals = self.logical_z

        for index, logical_op in enumerate(logicals):
            logical_out = 0
            for q in logical_op:
                logical_out += final_outs[-q - 1]
            logical_out = logical_out % 2

            if all_logicals or str(logical_out) != logical:
                node = DecodingGraphNode(
                    is_boundary=True,
                    qubits=logical,
                    index=index,
                )
                node.properties["basis"] = self.basis
                nodes.append(node)

        return nodes

    def check_nodes(self, nodes, ignore_extra_boundary=False, minimal=False):
        raise NotImplementedError

    def is_cluster_neutral(self, atypical_nodes):
        raise NotImplementedError

    def stim_circuit_with_detectors(self):
        """Converts the qiskit circuits into stim ciruits and add detectors.
        This is required for the stim-based construction of the DecodingGraph.
        """
        stim_circuits, _ = get_stim_circuits(self.noisy_circuit)
        measurements_per_cycle = len(self.x_gauges) + len(self.z_gauges)

        if self.round_schedule[0] == "x":
            measurement_round_offset = [0, len(self.x_gauges)]
        else:
            measurement_round_offset = [len(self.z_gauges), 0]

        ## 0th round of measurements
        if self.basis == "x":
            for stabind, stabilizer in enumerate(self.x_stabilizers):
                record_targets = []
                for gauge_ind in self._gauges4stabilizers[0][stabind]:
                    record_targets.append(
                        StimTarget_rec(
                            measurement_round_offset[0]
                            + gauge_ind
                            - (self.T * measurements_per_cycle + self.code.n)
                        )
                    )
                qubits_and_time = stabilizer.copy()
                qubits_and_time.extend([0])
                stim_circuits["0"].append("DETECTOR", record_targets, qubits_and_time)
                stim_circuits["1"].append("DETECTOR", record_targets, qubits_and_time)
        else:
            for stabind, stabilizer in enumerate(self.z_stabilizers):
                record_targets = []
                for gauge_ind in self._gauges4stabilizers[1][stabind]:
                    record_targets.append(
                        StimTarget_rec(
                            measurement_round_offset[1]
                            + gauge_ind
                            - (self.T * measurements_per_cycle + self.code.n)
                        )
                    )
                qubits_and_time = stabilizer.copy()
                qubits_and_time.extend([0])
                stim_circuits["0"].append("DETECTOR", record_targets, qubits_and_time)
                stim_circuits["1"].append("DETECTOR", record_targets, qubits_and_time)

        # adding first x and then z stabilizer comparisons
        for j in range(2):
            circuit = StimCircuit()
            for t in range(
                1, self.T
            ):  # compare stabilizer measurements with previous in each round
                for gind, gs in enumerate(self._gauges4stabilizers[j]):
                    record_targets = []
                    for gauge_ind in gs:
                        record_targets.append(
                            StimTarget_rec(
                                t * measurements_per_cycle
                                + measurement_round_offset[j]
                                + gauge_ind
                                - (self.T * measurements_per_cycle + self.code.n)
                            )
                        )
                        record_targets.append(
                            StimTarget_rec(
                                (t - 1) * measurements_per_cycle
                                + measurement_round_offset[j]
                                + gauge_ind
                                - (self.T * measurements_per_cycle + self.code.n)
                            )
                        )
                    qubits_and_time = self._stabilizers[j][gind].copy()
                    qubits_and_time.extend([t])
                    circuit.append("DETECTOR", record_targets, qubits_and_time)
            stim_circuits["0"] += circuit
            stim_circuits["1"] += circuit

        ## final measurements
        if self.basis == "x":
            for stabind, stabilizer in enumerate(self.x_stabilizers):
                record_targets = []
                for q in stabilizer:
                    record_targets.append(StimTarget_rec(q - self.code.n))
                for gauge_ind in self._gauges4stabilizers[0][stabind]:
                    record_targets.append(
                        StimTarget_rec(
                            measurement_round_offset[0]
                            + gauge_ind
                            - self.code.n
                            - measurements_per_cycle
                        )
                    )
                qubits_and_time = stabilizer.copy()
                qubits_and_time.extend([self.T])
                stim_circuits["0"].append("DETECTOR", record_targets, qubits_and_time)
                stim_circuits["1"].append("DETECTOR", record_targets, qubits_and_time)
            stim_circuits["0"].append(
                "OBSERVABLE_INCLUDE",
                [StimTarget_rec(q - self.code.n) for q in sorted(self.logical_x[0])],
                0,
            )
            stim_circuits["1"].append(
                "OBSERVABLE_INCLUDE",
                [StimTarget_rec(q - self.code.n) for q in sorted(self.logical_x[0])],
                0,
            )
        else:
            for stabind, stabilizer in enumerate(self.z_stabilizers):
                record_targets = []
                for q in stabilizer:
                    record_targets.append(StimTarget_rec(q - self.code.n))
                for gauge_ind in self._gauges4stabilizers[1][stabind]:
                    record_targets.append(
                        StimTarget_rec(
                            measurement_round_offset[1]
                            + gauge_ind
                            - self.code.n
                            - measurements_per_cycle
                        )
                    )
                qubits_and_time = stabilizer.copy()
                qubits_and_time.extend([self.T])
                stim_circuits["0"].append("DETECTOR", record_targets, qubits_and_time)
                stim_circuits["1"].append("DETECTOR", record_targets, qubits_and_time)
            stim_circuits["0"].append(
                "OBSERVABLE_INCLUDE",
                [StimTarget_rec(q - self.code.n) for q in sorted(self.logical_z[0])],
                0,
            )
            stim_circuits["1"].append(
                "OBSERVABLE_INCLUDE",
                [StimTarget_rec(q - self.code.n) for q in sorted(self.logical_z[0])],
                0,
            )

        return stim_circuits

    def _make_syndrome_graph(self):
        stim_circuit = self.stim_circuit_with_detectors()["0"]
        e = stim_circuit.detector_error_model(
            decompose_errors=True, approximate_disjoint_errors=True
        )
        graph, hyperedges = detector_error_model_to_rx_graph(e)
        return graph, hyperedges
