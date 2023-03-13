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

# pylint: disable=invalid-name

"""Tools to use functionality from Stim."""

from stim import Circuit as StimCircuit
from qiskit_aer.noise.errors.quantum_error import QuantumChannelInstruction


def get_stim_circuits(circuit_dict):
    """Returns a list of dictionaries. The first one contains the stim circuits for the two logicals,
    the second contatins the data how the stim measurement outcomes are ordered.

    Args:
        circuit_dict: Dictionary with Qiskit circuits as values. Compatible gates are paulis,
        controlled paulis, h, s, and sdg, swap, reset, measure and barrier. Compatible noise
        operators correspond to a single or two qubit pauli channel.

    Returns:
        stim_circuits: TODO
        stim_measurement_data: TODO
    """
    stim_circuits = {}
    stim_measurement_data = {}
    for circ_label, circuit in circuit_dict.items():
        stim_circuit = StimCircuit()

        # Dictionaries are not complete. For the stim definitions see:
        #   https://github.com/quantumlib/Stim/blob/main/doc/gates.md
        qiskit_to_stim_dict = {
            "id": "I",
            "x": "X",
            "y": "Y",
            "z": "Z",
            "h": "H",
            "s": "S",
            "sdg": "S_DAG",
            "cx": "CNOT",
            "cy": "CY",
            "cz": "CZ",
            "swap": "SWAP",
            "reset": "R",
            "measure": "M",
            "barrier": "TICK",
        }
        pauli_error_1_stim_order = {"id": 0, "I": 0, "X": 1, "x": 1, "Y": 2, "y": 2, "Z": 3, "z": 3}
        pauli_error_2_stim_order = {
            "II": 0,
            "IX": 1,
            "IY": 2,
            "IZ": 3,
            "XI": 4,
            "XX": 5,
            "XY": 6,
            "XZ": 7,
            "YI": 8,
            "YX": 9,
            "YY": 10,
            "YZ": 11,
            "ZI": 12,
            "ZX": 13,
            "ZY": 14,
            "ZZ": 15,
        }

        measurement_data = []
        register_offset = {}
        previous_offset = 0
        for inst, qargs, cargs in circuit.data:
            for qubit in qargs:
                if qubit._register.name not in register_offset:
                    register_offset[qubit._register.name] = previous_offset
                    previous_offset += qubit._register.size

            qubit_indices = [
                qargs[i]._index + register_offset[qargs[i]._register.name]
                for i in range(len(qargs))
            ]

            if isinstance(inst, QuantumChannelInstruction):
                # Errors: it would be good to add exeptions for delpoarizinng noise, Z error etc,
                # because the stim detctor error model relies on the disjoint error assumption
                # for general Pauli channel
                qerror = inst._quantum_error
                pauli_errors_types = qerror.to_dict()["instructions"]
                pauli_probs = qerror.to_dict()["probabilities"]
                if pauli_errors_types[0][0]["name"] in pauli_error_1_stim_order:
                    probs = 4 * [0.0]
                    for pind, ptype in enumerate(pauli_errors_types):
                        probs[pauli_error_1_stim_order[ptype[0]["name"]]] = pauli_probs[pind]
                    stim_circuit.append("PAULI_CHANNEL_1", qubit_indices, probs[1:])
                elif pauli_errors_types[0][0]["params"][0] in pauli_error_2_stim_order:
                    # here the name is always 'pauli' and the params gives the Pauli type
                    probs = 16 * [0.0]
                    for pind, ptype in enumerate(pauli_errors_types):
                        probs[pauli_error_2_stim_order[ptype[0]["params"][0]]] = pauli_probs[pind]
                    stim_circuit.append("PAULI_CHANNEL_2", qubit_indices, probs[1:])
                else:
                    raise Exception("Unexpected operations: " + str([inst, qargs, cargs]))
            else:
                # Gates and measurements
                if inst.name in qiskit_to_stim_dict:
                    if len(cargs) > 0:  # keeping track of measurement indices in stim
                        measurement_data.append(
                            [
                                cargs[0]._index + register_offset[qargs[0]._register.name],
                                qargs[0]._register.name,
                            ]
                        )
                    if qiskit_to_stim_dict[inst.name] == "TICK":  # barrier
                        stim_circuit.append("TICK")
                    else:  # gates/measurements acting on qubits
                        stim_circuit.append(qiskit_to_stim_dict[inst.name], qubit_indices)
                else:
                    raise Exception("Unexpected operations: " + str([inst, qargs, cargs]))

        stim_circuits[circ_label] = stim_circuit
        stim_measurement_data[circ_label] = measurement_data
    return stim_circuits, stim_measurement_data


def get_counts_via_stim(circuit, shots: int = 4000):
    """Returns a qiskit compatible dictionary of measurement outcomes

    Args:
        circuit: Qiskit circuit compatible with `get_stim_circuits`.
        shots: Number of samples to be generated.

    Returns:
        counts: Counts dictionary in standard Qiskit form.
    """

    stim_circuits, stim_measurement_data = get_stim_circuits({"": circuit})
    stim_circuit = stim_circuits[""]
    measurement_data = stim_measurement_data[""]

    stim_samples = stim_circuit.compile_sampler().sample(shots=shots)
    qiskit_counts = {}
    for stim_sample in stim_samples:
        prev_reg = measurement_data[-1][1]
        qiskit_count = ""
        for idx, meas in enumerate(measurement_data[::-1]):
            _, reg = meas
            if reg != prev_reg:
                qiskit_count += " "
            qiskit_count += str(int(stim_sample[-idx - 1]))
            prev_reg = reg
        if qiskit_count in qiskit_counts:
            qiskit_counts[qiskit_count] += 1
        else:
            qiskit_counts[qiskit_count] = 1

    return qiskit_counts
