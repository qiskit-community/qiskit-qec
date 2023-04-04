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

# pylint: disable=invalid-name, disable=no-name-in-module

"""Tools to use functionality from Stim."""

from stim import Circuit as StimCircuit
from stim import DetectorErrorModel as StimDetectorErrorModel
from stim import DemInstruction as StimDemInstruction
from stim import DemTarget as StimDemTarget

from qiskit import QuantumCircuit
from qiskit_aer.noise.errors.quantum_error import QuantumChannelInstruction
from qiskit_aer.noise import pauli_error
from qiskit_qec.utils.decoding_graph_attributes import DecodingGraphNode, DecodingGraphEdge

import rustworkx as rx
from typing import List, Dict
import numpy as np
from math import log

def get_stim_circuits(circuit_dict):
    """Converts compatible qiskit circuits to stim circuits.

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


def get_counts_via_stim(circuits, shots: int = 4000, noise_model=None):
    """Returns a qiskit compatible dictionary of measurement outcomes

    Args:
        circuit: Qiskit circuit compatible with `get_stim_circuits` or list thereof.
        shots: Number of samples to be generated.
        noise_model: Pauli noise model for any additional noise to be applied.

    Returns:
        counts: Counts dictionary in standard Qiskit form or list thereof.
    """

    if noise_model:
        circuits = noisify_circuit(circuits, noise_model)

    single_circuit = isinstance(circuits, QuantumCircuit)
    if single_circuit:
        circuits = [circuits]

    counts = []
    for circuit in circuits:

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
        counts.append(qiskit_counts)

    if single_circuit:
        counts = counts[0]

    return counts

def detector_error_model_to_rx_graph(model: StimDetectorErrorModel) -> rx.PyGraph:
    """Convert a stim error model into a RustworkX graph.
       It assumes that the stim circuit does not contain repeat blocks.
       Later on repeat blocks should be handled to make this function compatible with
       user-defined stim circuits.
    """

    g = rx.PyGraph(multigraph=False)

    index_to_DecodingGraphNode = {}

    for instruction in model:
        if instruction.type == "detector":
            a = np.array(instruction.args_copy())
            time = a[-1]
            qubits = [int(qubit_ind) for qubit_ind in a[:-1]]
            for t in instruction.targets_copy():
                node = DecodingGraphNode(index=t.val, time=time, qubits=qubits)
                index_to_DecodingGraphNode[t.val]=node
                g.add_node(node)
    
    trivial_boundary_node = DecodingGraphNode(index=model.num_detectors, time=0, is_boundary=True)
    g.add_node(trivial_boundary_node)
    index_to_DecodingGraphNode[model.num_detectors]=trivial_boundary_node

    def handle_error(p: float, dets: List[int], frame_changes: List[int], hyperedge: Dict):
        if p == 0:
            return
        if len(dets) == 0:
            return
        if len(dets) == 1:
            dets = [dets[0], model.num_detectors]
            # if frame_changes == []:
            #     dets = [dets[0], model.num_detectors]
            # else:
            #     dets = [dets[0], model.num_detectors+1]
        if len(dets) > 2:
            raise NotImplementedError(
                f"Error with more than 2 symptoms can't become an edge or boundary edge: {dets!r}.")
        if g.has_edge(dets[0],dets[1]):
            edge_ind = [dets for dets in g.edge_list()].index((dets[0],dets[1]))
            edge_data = g.edges()[edge_ind].properties
            old_p = edge_data["error_probability"]
            old_frame_changes = edge_data["fault_ids"]
            # If frame changes differ, the code has distance 2; just keep whichever was first.
            if set(old_frame_changes) == set(frame_changes):
                p = p * (1 - old_p) + old_p * (1 - p)
                g.remove_edge(dets[0],dets[1])
        if p > 0.5:
            p = 1 - p
        if p > 0:
            qubits = list(set(index_to_DecodingGraphNode[dets[0]].qubits).intersection(index_to_DecodingGraphNode[dets[1]].qubits))
            edge = DecodingGraphEdge(qubits=qubits,weight=log((1 - p) / p), properties={'fault_ids': set(frame_changes), 'error_probability': p})
            g.add_edge(dets[0],dets[1], edge)
            hyperedge[dets[0],dets[1]] = edge

    hyperedges = []

    for instruction in model:
        if isinstance(instruction, StimDemInstruction):
            if instruction.type == "error":
                dets: List[int] = []
                frames: List[int] = []
                t: StimDemTarget
                p = instruction.args_copy()[0]
                hyperedge = {}
                for t in instruction.targets_copy():
                    if t.is_relative_detector_id():
                        dets.append(t.val)
                    elif t.is_logical_observable_id():
                        frames.append(t.val)
                    elif t.is_separator():
                        # Treat each component of a decomposed error as an independent error.
                        handle_error(p, dets, frames, hyperedge)
                        frames = []
                        dets = []
                # Handle last component.
                handle_error(p, dets, frames, hyperedge)
                if len(hyperedge)>1:
                    hyperedges.append(hyperedge)
            elif instruction.type == "detector":
                pass
            elif instruction.type == "logical_observable":
                pass
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()

    return g, hyperedges

def noisify_circuit(circuits, noise_model):
    """
    Inserts error operations into a circuit according to a pauli noise model.

    Args:
        circuits: Circuit or list thereof.
        noise_model: Pauli noise model.

    Returns:
        noisy_circuits: Corresponding circuit or list thereof.
    """

    single_circuit = isinstance(circuits, QuantumCircuit)
    if single_circuit:
        circuits = [circuits]

    # create pauli errors for all errors in noise model
    errors = {}
    for g, noise in noise_model.to_dict().items():
        errors[g] = []
        for pauli, prob in noise["chan"].items():
            pauli = pauli.upper()
            errors[g].append(pauli_error([(pauli, prob), ("I" * len(pauli), 1 - prob)]))

    noisy_circuits = []
    for qc in circuits:

        noisy_qc = QuantumCircuit()
        for qreg in qc.qregs:
            noisy_qc.add_register(qreg)
        for creg in qc.cregs:
            noisy_qc.add_register(creg)

        for gate in qc:
            g = gate[0].name
            qubits = gate[1]
            pre_error = g == "reset"
            # add gate if it needs to go before the error
            if pre_error:
                noisy_qc.append(gate)
            # then the error
            if g in errors:
                for error_op in errors[g]:
                    noisy_qc.append(error_op, qubits)
            # add gate if it needs to go after the error
            if not pre_error:
                noisy_qc.append(gate)

        noisy_circuits.append(noisy_qc)

    if single_circuit:
        noisy_circuits = noisy_circuits[0]

    return noisy_circuits
