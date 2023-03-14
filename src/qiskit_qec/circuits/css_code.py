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

"""Generates circuits for CSS codes."""
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer.noise import depolarizing_error, pauli_error

from qiskit_qec.circuits.code_circuit import CodeCircuit
from qiskit_qec.utils.stim_tools import noisify_circuit
from qiskit_qec.codes import StabSubSystemCode
from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.linear.symplectic import normalizer


class CssCodeCircuit(CodeCircuit):
    """
    CodeCircuit class for generic CSS codes.
    """

    def __init__(self, code, T, basis="z", round_schedule="zx", noise_model=None):
        """
        Args:
            code: A CSS code
            T: Number of syndrome measurement rounds
            basis: basis for encoding ('x' or 'z')
            round_schedule: Order in which to measureme gauge operators ('zx' or 'xz')
            noise_model: Pauli noise model used in the construction of noisy circuits.
            If a tuple, a pnenomological noise model is used with the entries being
            probabity of depolarizing noise on code qubits between rounds and
            probability of measurement errors, respectively.
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
            # prepare initial state
            if state[0] == "1":
                if basis == "z":
                    qc.x(self.logical_x[0])
                else:
                    qc.x(self.logical_z[0])
            if basis == "x":
                qc.h(qregs[0])
            # peform syndrome measurements
            for t in range(T):
                if state[-1] == "n" and self._phenom:
                    for q in qregs[0]:
                        qc.append(self._depol_error, [q])
                # gauge measurements
                if round_schedule == "zx":
                    self._z_gauge_measurements(qc, t, state)
                    self._x_gauge_measurements(qc, t, state)
                elif round_schedule == "xz":
                    self._x_gauge_measurements(qc, t, state)
                    self._z_gauge_measurements(qc, t, state)
                else:
                    print("Round schedule " + round_schedule + " not supported.")
            # final readout
            creg = ClassicalRegister(code.n, name="final_readout")
            qc.add_register(creg)
            if basis == "x":
                qc.h(qregs[0])
            if state[-1] == "n" and self._phenom:
                for q in qregs[0]:
                    qc.append(self._meas_error, [q])
            qc.measure(qregs[0], creg)
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

            raw_gauges = code.gauge_group.generators
            center, log, conj_log = normalizer(code.generators.matrix)
            raw_stabilizers = PauliList(center)
            raw_logicals = PauliList(log) + PauliList(conj_log)

            gauges = [[],[]]
            stabilizers = [[],[]]
            logicals = [[],[]]

            for raw_ops, ops, in zip(
                [raw_gauges,raw_stabilizers,raw_logicals],
                [gauges, stabilizers, logicals]
                ):

                for op in raw_ops:
                    op = str(op)
                    for j, pauli in enumerate(['X','Z']):
                        if (op.count(pauli) + op.count('I')) == code.n:
                            ops[j].append([k for k,p in enumerate(op[::-1]) if p==pauli])
                is_css = is_css and (len(ops[0]) + len(ops[1])) == len(raw_ops)

            if is_css:
                self.x_gauges = gauges[0]
                self.z_gauges = gauges[1]
                self.x_stabilizers = stabilizers[0]
                self.z_stabilizers = stabilizers[1]
                self.logical_x = logicals[0]
                self.logical_z = logicals[1]
            else:
                raise

        else:
            # otherwise assume it has the info
            self.x_gauges = self.code.x_gauges
            self.z_gauges = self.code.z_gauges
            self.x_stabilizers = self.code.x_stabilizers
            self.z_stabilizers = self.code.z_stabilizers
            self.logical_x = self.code.logical_x
            self.logical_z = self.code.logical_z

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
        """
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
        round_outs = []
        for stabilizer in self._stabilizers[j]:
            out = 0
            for q in stabilizer:
                out += final_outs[q]
            out = out % 2
            round_outs.append(out)
        stabilizer_outs[j].append(round_outs)

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
                        node = {
                            "time": t,
                            "basis": bases[j],
                            "qubits": self._stabilizers[j][e],
                            "element": e,
                            "is_boundary": False,
                        }
                        nodes.append(node)

        if self.basis == 'x':
            logicals = self.logical_x
        else:
            logicals = self.logical_z

        for index, logical in enumerate(logicals):
            logical_out = 0
            for q in logical:
                logical_out += final_outs[q]
            logical_out = logical_out % 2

            if logical_out == 1:
                node = {
                    "time": 0,
                    "basis": self.basis,
                    "qubits": logical,
                    "element": index,
                    "is_boundary": True,
                }
                nodes.append(node)

        return nodes

    def check_nodes(self, nodes, ignore_extra_boundary=False):
        pass

    def is_cluster_neutral(self, atypical_nodes):
        pass
