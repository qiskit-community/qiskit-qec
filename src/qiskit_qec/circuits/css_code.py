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


class CssCodeCircuit(CodeCircuit):
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

        super().__init__()

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

        self.circuit = {}
        for state in ["0", "1"]:
            qc = QuantumCircuit()
            qregs = []
            qregs.append(QuantumRegister(code.n, name="code qubits"))
            qregs.append(QuantumRegister(len(code.z_gauges), name="z auxs"))
            qregs.append(QuantumRegister(len(code.x_gauges), name="x auxs"))
            for qreg in qregs:
                qc.add_register(qreg)
            # prepare initial state
            if state == "1":
                if basis == "z":
                    qc.x(code.logical_x[0])
                else:
                    qc.x(code.logical_z[0])
            if basis == "x":
                qc.h(qregs[0])
            # peform syndrome measurements
            for t in range(T):
                if self._noise:
                    for q in qregs[0]:
                        qc.append(self._depol_error, [q])
                # gauge measurements
                if round_schedule == "zx":
                    self._z_gauge_measurements(qc, t)
                    self._x_gauge_measurements(qc, t)
                elif round_schedule == "xz":
                    self._x_gauge_measurements(qc, t)
                    self._z_gauge_measurements(qc, t)
                else:
                    print("Round schedule " + round_schedule + " not supported.")
            # final readout
            creg = ClassicalRegister(code.n, name="final_readout")
            qc.add_register(creg)
            if basis == "x":
                qc.h(qregs[0])
            if self._noise:
                for q in qregs[0]:
                    qc.append(self._meas_error, [q])
            qc.measure(qregs[0], creg)
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

    def _z_gauge_measurements(self, qc, t):
        creg = ClassicalRegister(len(self.code.z_gauges), name="round_" + str(t) + "_z_bits")
        qc.add_register(creg)
        for g, z_gauge in enumerate(self.code.z_gauges):
            for q in z_gauge:
                qc.cx(qc.qregs[0][q], qc.qregs[1][g])
            if self._noise:
                qc.append(self._meas_error, [qc.qregs[1][g]])
            qc.measure(qc.qregs[1][g], creg[g])
            qc.reset(qc.qregs[1][g])

    def _x_gauge_measurements(self, qc, t):
        creg = ClassicalRegister(len(self.code.x_gauges), name="round_" + str(t) + "_x_bits")
        qc.add_register(creg)
        for g, x_gauge in enumerate(self.code.x_gauges):
            for q in x_gauge:
                qc.h(qc.qregs[0][q])
                qc.cx(qc.qregs[0][q], qc.qregs[2][g])
                qc.h(qc.qregs[0][q])
            if self._noise:
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
    
    def to_stim_circuit(self):
        '''Returns a list of dictionaries. The first one contains the stim circuits for the two logicals,
           the second contatins the data how the stim measurement outcomes are ordered
        '''
        stim_circuits = {}
        stim_measurement_data = {}
        for circ_label, circuit in self.circuit.items():
            stim_circuit = stim.Circuit()

            '''Dictionaries are not complete. For the stim definitions see: 
               https://github.com/quantumlib/Stim/blob/main/doc/gates.md
            ''' 
            qiskit_to_stim_dict = {'id':'I','x':'X','y':'Y','z':'Z','h':'H','s':'S','sdg':'S_DAG',
                                'cx':'CNOT','cy':'CY','cz':'CZ','swap':'SWAP',
                                'reset':'R','measure':'M','barrier':'TICK'}
            pauli_error_1_stim_order = {'id':0,'I':0,'X':1,'x':1,'Y':2,'y':2,'Z':3,'z':3}
            pauli_error_2_stim_order = {'II':0,'IX':1,'IY':2,'IZ':3,
                                        'XI':4,'XX':5,'XY':6,'XZ':7,
                                        'YI':8,'YX':9,'YY':10,'YZ':11,
                                        'ZI':12,'ZX':13,'ZY':14,'ZZ':15} 

            measurement_data = []
            register_offset = {}
            previous_offset = 0
            for inst, qargs, cargs in circuit.data:
                for qubit in qargs:
                    if qubit._register.name not in register_offset:
                        register_offset[qubit._register.name] = previous_offset
                        previous_offset += qubit._register.size

                qubit_indices = [qargs[i]._index+register_offset[qargs[i]._register.name] for i in range(len(qargs))]

                if isinstance(inst, QuantumChannelInstruction):
                    '''Errors: it would be good to add exeptions for delpoarizinng noise, Z error etc,
                    because the stim detctor error model relies on the disjoint error assumption 
                    for general Pauli channel
                    '''
                    qerror = inst._quantum_error
                    pauli_errors_types=qerror.to_dict()['instructions']
                    pauli_probs = qerror.to_dict()['probabilities']
                    if pauli_errors_types[0][0]['name'] in pauli_error_1_stim_order.keys():
                        probs = 4*[0.]
                        for pind,ptype in enumerate(pauli_errors_types):
                            probs[pauli_error_1_stim_order[ptype[0]['name']]] = pauli_probs[pind]
                        stim_circuit.append("PAULI_CHANNEL_1",qubit_indices,probs[1:])
                    elif pauli_errors_types[0][0]['params'][0] in pauli_error_2_stim_order.keys(): 
                        #here the name is always 'pauli' and the params gives the Pauli type
                        probs = 16*[0.]
                        for pind,ptype in enumerate(pauli_errors_types):
                            probs[pauli_error_2_stim_order[ptype[0]['params'][0]]] = pauli_probs[pind]
                        stim_circuit.append("PAULI_CHANNEL_2",qubit_indices,probs[1:])
                    else:
                        error_types = qerror.to_dict()['instructions']
                        print("I didn't see this comming: "+str([inst, qargs, cargs]))+str(qerror)
                else:
                    '''Gates and measurements'''
                    if inst.name in qiskit_to_stim_dict:
                        if len(cargs)>0: #keeping track of measurement indices in stim
                            measurement_data.append([cargs[0]._index+register_offset[qargs[0]._register.name],qargs[0]._register.name])
                        if qiskit_to_stim_dict[inst.name] == "TICK": #barrier
                            stim_circuit.append("TICK")
                        else: #gates/measurements acting on qubits
                            stim_circuit.append(qiskit_to_stim_dict[inst.name],qubit_indices)
                    else: 
                        print("I didn't see this comming: "+str([inst, qargs, cargs]))

            stim_circuits[circ_label] = stim_circuit
            stim_measurement_data[circ_label] = measurement_data
        return [stim_circuits, stim_measurement_data]
    
    def get_counts_via_stim(self, logical: str = '0', shots: int = 1024):
        '''Returns a qiskit compatible dictionary of measurement outcomes'''
        stim_circuits, stim_measurement_data = self.to_stim_circuit()
        stim_circuit = stim_circuits[logical]
        measurement_data = stim_measurement_data[logical]

        stim_samples = stim_circuit.compile_sampler().sample(shots=shots)
        qiskit_counts = {}
        for stim_sample in stim_samples:
            prev_reg = measurement_data[-1][1]
            qiskit_count = ''
            for idx,meas in enumerate(measurement_data[::-1]):
                measind,reg = meas
                if reg != prev_reg:
                    qiskit_count += ' '
                qiskit_count += str(int(stim_sample[-idx-1]))
                prev_reg = reg
            if qiskit_count in qiskit_counts.keys():
                qiskit_counts[qiskit_count] = qiskit_counts[qiskit_count]+1
            else:
                qiskit_counts[qiskit_count] = 1
        
        return qiskit_counts
