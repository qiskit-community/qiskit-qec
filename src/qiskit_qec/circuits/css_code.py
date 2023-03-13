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

import stim
from qiskit_aer.noise.errors.quantum_error import QuantumChannelInstruction

class CssCodeCircuit():
    """
    CodeCircuit class for generic CSS codes.
    """

    def __init__(self, code, T, basis="z", round_schedule="zx", p_depol=0, p_meas=0):
        """
        Args:
            code: A CSS code
            T: Number of syndrome measurement rounds
            basis: basis for encoding ('x' or 'z')
            round_schedule: Order in which to measureme gauge operators ('zx' or 'xz')
            p_depol: Probabity of depolarizing noise on code qubits between rounds
            p_meas: Probability of measurement errors
        """

        self.code = code
        self.T = T
        self.basis = basis
        self.base = "0"
        self.round_schedule = round_schedule
        self.p_depol = p_depol
        self.p_meas = p_meas
        self._noise = p_depol > 0 or p_meas > 0

        self._depol_error = depolarizing_error(p_depol, 1)
        self._meas_error = pauli_error([("X", p_meas), ("I", 1 - p_meas)])

        circuit = {}
        states = ["0", "1"]
        if self._noise:
            states += ["0n", "1n"]
        for state in states:
            qc = QuantumCircuit()
            qregs = []
            qregs.append(QuantumRegister(code.n, name="code qubits"))
            qregs.append(QuantumRegister(len(code.z_gauges), name="z auxs"))
            qregs.append(QuantumRegister(len(code.x_gauges), name="x auxs"))
            for qreg in qregs:
                qc.add_register(qreg)
            # prepare initial state
            if state[0] == "1":
                if basis == "z":
                    qc.x(code.logical_x[0])
                else:
                    qc.x(code.logical_z[0])
            if basis == "x":
                qc.h(qregs[0])
            # peform syndrome measurements
            for t in range(T):
                if state[-1] == 'n':
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
            if state[-1] == 'n':
                for q in qregs[0]:
                    qc.append(self._meas_error, [q])
            qc.measure(qregs[0], creg)
            circuit[state] = qc

        self.circuit = {}
        self.noisy_circuit = {}
        for state, qc in circuit.items():
            if state[-1] == 'n':
                self.noisy_circuit[state[0]] = qc
            else:
                self.circuit[state] = qc

        self._gauges4stabilizers = []
        self._stabilizers = [code.x_stabilizers, code.z_stabilizers]
        self._gauges = [code.x_gauges, code.z_gauges]
        for j in range(2):
            self._gauges4stabilizers.append([])
            for stabilizer in self._stabilizers[j]:
                gauges = []
                for g, gauge in enumerate(self._gauges[j]):
                    if set(stabilizer).intersection(set(gauge)) == set(gauge):
                        gauges.append(g)
                self._gauges4stabilizers[j].append(gauges)

    def _z_gauge_measurements(self, qc, t, state):
        creg = ClassicalRegister(len(self.code.z_gauges), name="round_" + str(t) + "_z_bits")
        qc.add_register(creg)
        for g, z_gauge in enumerate(self.code.z_gauges):
            for q in z_gauge:
                qc.cx(qc.qregs[0][q], qc.qregs[1][g])
            if state[-1] == 'n':
                qc.append(self._meas_error, [qc.qregs[1][g]])
            qc.measure(qc.qregs[1][g], creg[g])
            qc.reset(qc.qregs[1][g])

    def _x_gauge_measurements(self, qc, t, state):
        creg = ClassicalRegister(len(self.code.x_gauges), name="round_" + str(t) + "_x_bits")
        qc.add_register(creg)
        for g, x_gauge in enumerate(self.code.x_gauges):
            for q in x_gauge:
                qc.h(qc.qregs[0][q])
                qc.cx(qc.qregs[0][q], qc.qregs[2][g])
                qc.h(qc.qregs[0][q])
            if state[-1] == 'n':
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

        j = bases.index(self.basis)
        boundary = [self.code.x_boundary, self.code.z_boundary]
        boundary_qubits = [q[0] for q in boundary[j]]

        boundary_out = 0
        for q in boundary_qubits:
            boundary_out += final_outs[q]
        boundary_out = boundary_out % 2

        if boundary_out == 1:
            node = {
                "time": 0,
                "basis": self.basis,
                "qubits": boundary_qubits,
                "element": 0,
                "is_boundary": True,
            }
            nodes.append(node)

        return nodes
    
