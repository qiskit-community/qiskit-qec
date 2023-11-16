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

# pylint: disable=invalid-name, disable=no-name-in-module, disable=unused-argument

"""Tools to use functionality from Stim."""
from typing import Union, List, Dict, Callable
from math import log as loga
from stim import Circuit as StimCircuit
from stim import DetectorErrorModel as StimDetectorErrorModel
from stim import DemInstruction as StimDemInstruction
from stim import DemRepeatBlock as StimDemRepeatBlock
from stim import DemTarget as StimDemTarget
from stim import target_rec as StimTarget_rec

import numpy as np
import rustworkx as rx

from qiskit import QuantumCircuit
from qiskit_aer.noise.errors.quantum_error import QuantumChannelInstruction
from qiskit_aer.noise import pauli_error
from qiskit_qec.utils.decoding_graph_attributes import (
    DecodingGraphNode,
    DecodingGraphEdge,
)
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


def get_stim_circuits(
    circuit: Union[QuantumCircuit, List],
    detectors: List[Dict] = None,
    logicals: List[Dict] = None,
):
    """Converts compatible qiskit circuits to stim circuits.
       Dictionaries are not complete. For the stim definitions see:
       https://github.com/quantumlib/Stim/blob/main/doc/gates.md
    Args:
        circuit: Compatible gates are Paulis, controlled Paulis, h, s,
        and sdg, swap, reset, measure and barrier. Compatible noise operators
        correspond to a single or two qubit pauli channel.
        detectors: A list of measurement comparisons. A measurement comparison
        (detector) is either a list of measurements given by a the name and index
        of the classical bit or a list of dictionaries, with a mandatory clbits
        key containing the classical bits. A dictionary can contain keys like
        'qubits', 'time', 'basis' etc.
        logicals: A list of logical measurements. A logical measurement is a
        list of classical bits whose total parity is the logical eigenvalue.
        Again it can be a list of dictionaries.

    Returns:
        stim_circuits, stim_measurement_data
    """

    if detectors is None:
        detectors = [{}]
    if logicals is None:
        logicals = [{}]

    if len(detectors) > 0 and isinstance(detectors[0], List):
        detectors = [{"clbits": det, "qubits": [], "time": 0} for det in detectors]

    if len(logicals) > 0 and isinstance(logicals[0], List):
        logicals = [{"clbits": log} for log in logicals]

    stim_circuits = []
    stim_measurement_data = []
    if isinstance(circuit, QuantumCircuit):
        circuit = [circuit]
    for circ in circuit:
        stim_circuit = StimCircuit()

        qiskit_to_stim_dict = {
            "id": "I",
            "x": "X",
            "y": "Y",
            "z": "Z",
            "h": "H",
            "s": "S",
            "sdg": "S_DAG",
            "cx": "CX",
            "cy": "CY",
            "cz": "CZ",
            "swap": "SWAP",
            "reset": "R",
            "measure": "M",
            "barrier": "TICK",
        }
        pauli_error_1_stim_order = {
            "id": 0,
            "I": 0,
            "X": 1,
            "x": 1,
            "Y": 2,
            "y": 2,
            "Z": 3,
            "z": 3,
        }
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
        qreg_offset = {}
        creg_offset = {}
        prevq_offset = 0
        prevc_offset = 0
        for inst, qargs, cargs in circ.data:
            for qubit in qargs:
                if qubit._register.name not in qreg_offset:
                    qreg_offset[qubit._register.name] = prevq_offset
                    prevq_offset += qubit._register.size
            for bit in cargs:
                if bit._register.name not in creg_offset:
                    creg_offset[bit._register.name] = prevc_offset
                    prevc_offset += bit._register.size

            qubit_indices = [
                qargs[i]._index + qreg_offset[qargs[i]._register.name] for i in range(len(qargs))
            ]

            if isinstance(inst, QuantumChannelInstruction):
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
                        measurement_data.append([cargs[0]._register.name, cargs[0]._index])

                    if qiskit_to_stim_dict[inst.name] == "TICK":  # barrier
                        stim_circuit.append("TICK")
                    elif inst.condition is not None:  # handle c_ifs
                        if inst.name in "xyz":
                            if inst.condition[1] == 1:
                                clbit = inst.condition[0]
                                stim_circuit.append(
                                    qiskit_to_stim_dict["c" + inst.name],
                                    [
                                        StimTarget_rec(
                                            measurement_data.index(
                                                [clbit._register.name, clbit._index]
                                            )
                                            - len(measurement_data)
                                        ),
                                        qubit_indices[0],
                                    ],
                                )
                            else:
                                raise Exception(
                                    "Classically controlled gate must be conditioned on bit value 1"
                                )
                        else:
                            raise Exception(
                                "Classically controlled " + inst.name + " gate is not supported"
                            )
                    else:  # gates/measurements acting on qubits
                        stim_circuit.append(qiskit_to_stim_dict[inst.name], qubit_indices)
                else:
                    raise Exception("Unexpected operations: " + str([inst, qargs, cargs]))

        if detectors != [{}]:
            for det in detectors:
                stim_record_targets = []
                for reg, ind in det["clbits"]:
                    stim_record_targets.append(
                        StimTarget_rec(measurement_data.index([reg, ind]) - len(measurement_data))
                    )
                if det["time"] != []:
                    stim_circuit.append(
                        "DETECTOR", stim_record_targets, det["qubits"] + [det["time"]]
                    )
                else:
                    stim_circuit.append("DETECTOR", stim_record_targets, [])
        if logicals != [{}]:
            for log_ind, log in enumerate(logicals):
                stim_record_targets = []
                for reg, ind in log["clbits"]:
                    stim_record_targets.append(
                        StimTarget_rec(measurement_data.index([reg, ind]) - len(measurement_data))
                    )
                stim_circuit.append("OBSERVABLE_INCLUDE", stim_record_targets, log_ind)

        stim_circuits.append(stim_circuit)
        stim_measurement_data.append(measurement_data)

    return stim_circuits, stim_measurement_data


def get_counts_via_stim(
    circuits: Union[List, QuantumCircuit],
    shots: int = 4000,
    noise_model: PauliNoiseModel = None,
):
    """Returns a qiskit compatible dictionary of measurement outcomes

    Args:
        circuits: Qiskit circuit compatible with `get_stim_circuits` or list thereof.
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
        stim_circuits, stim_measurement_data = get_stim_circuits(circuit)
        stim_circuit = stim_circuits[0]
        measurement_data = stim_measurement_data[0]

        stim_samples = stim_circuit.compile_sampler().sample(shots=shots)
        qiskit_counts = {}
        for stim_sample in stim_samples:
            prev_reg = measurement_data[-1][0]
            qiskit_count = ""
            for idx, meas in enumerate(measurement_data[::-1]):
                reg, _ = meas
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


def iter_flatten_model(
    model: StimDetectorErrorModel,
    handle_error: Callable[[float, List[int], List[int]], None],
    handle_detector_coords: Callable[[int, np.ndarray], None],
    detectors: List[Dict],
    hyperedges: List[Dict],
):
    """
    This function have been copied from the built-in method of
    stim: stim.Circuit.generated("surface_code:rotated_memory_z",...)
    """

    det_offset = 0

    def _helper(m: StimDetectorErrorModel, reps: int):
        nonlocal det_offset
        for _ in range(reps):
            for instruction in m:
                if isinstance(instruction, StimDemRepeatBlock):
                    _helper(instruction.body_copy(), instruction.repeat_count)
                elif isinstance(instruction, StimDemInstruction):
                    if instruction.type == "error":
                        dets: List[int] = []
                        frames: List[int] = []
                        t: StimDemTarget
                        p = instruction.args_copy()[0]
                        hyperedge = {}
                        for t in instruction.targets_copy():
                            if t.is_relative_detector_id():
                                dets.append(t.val + det_offset)
                            elif t.is_logical_observable_id():
                                frames.append(t.val)
                            elif t.is_separator():
                                # Treat each component of a decomposed error as an independent error.
                                handle_error(p, dets, frames, hyperedge)
                                frames = []
                                dets = []
                        # Handle last component.
                        handle_error(p, dets, frames, hyperedge)
                        if len(hyperedge) > 1:
                            hyperedges.append(hyperedge)
                    elif instruction.type == "shift_detectors":
                        det_offset += instruction.targets_copy()[0]
                    elif instruction.type == "detector":
                        t = instruction.targets_copy()[0]
                        det_ind = t.val + det_offset
                        if detectors == [{}]:
                            a = np.array(instruction.args_copy())
                            time = a[-1]
                            qubits = [int(qubit_ind) for qubit_ind in a[:-1]]
                            det = {}
                        else:
                            det = detectors[det_ind].copy()
                            time = det.pop("time")
                            qubits = det.pop("qubits")
                            del det["clbits"]
                        for t in instruction.targets_copy():
                            handle_detector_coords(
                                detector_index=det_ind,
                                time=time,
                                qubits=qubits,
                                det_props=det,
                            )
                    elif instruction.type == "logical_observable":
                        pass
                    else:
                        raise NotImplementedError()
                else:
                    raise NotImplementedError()

    _helper(model, 1)


def detector_error_model_to_rx_graph(
    model: StimDetectorErrorModel, detectors: List[Dict] = None
) -> rx.PyGraph:
    """Convert a stim error model into a RustworkX graph.
    It assumes that the stim circuit does not contain repeat blocks.
    Later on repeat blocks should be handled to make this function compatible with
    user-defined stim circuits.

    Args:
        detectors:
        coordinate included as the last element for every detector in the stim detector error model
    """

    if detectors is None:
        detectors = [{}]

    g = rx.PyGraph(multigraph=False)

    index_to_DecodingGraphNode = {}

    def skip_error(p: float, dets: List[int], frame_changes: List[int], hyperedge: Dict):
        pass

    def handle_error(p: float, dets: List[int], frame_changes: List[int], hyperedge: Dict):
        if p == 0:
            return
        if len(dets) == 0:
            return
        if len(dets) == 1:
            dets = [dets[0], model.num_detectors]
        if len(dets) > 2:
            raise NotImplementedError(
                f"Error with more than 2 symptoms can't become an edge or boundary edge: {dets!r}."
            )
        if g.has_edge(dets[0], dets[1]):
            edge_ind = list(g.edge_list()).index((dets[0], dets[1]))
            edge_data = g.edges()[edge_ind].properties
            old_p = edge_data["error_probability"]
            old_frame_changes = edge_data["fault_ids"]
            # If frame changes differ, the code has distance 2; just keep whichever was first.
            if set(old_frame_changes) == set(frame_changes):
                p = p * (1 - old_p) + old_p * (1 - p)
                g.remove_edge(dets[0], dets[1])
        if p > 0.5:
            p = 1 - p
        if p > 0:
            qubits = list(
                set(index_to_DecodingGraphNode[dets[0]].qubits).intersection(
                    index_to_DecodingGraphNode[dets[1]].qubits
                )
            )
            edge = DecodingGraphEdge(
                qubits=qubits,
                weight=loga((1 - p) / p),
                properties={"fault_ids": set(frame_changes), "error_probability": p},
            )
            g.add_edge(dets[0], dets[1], edge)
            hyperedge[dets[0], dets[1]] = edge

    def skip_detector_coords(detector_index: int, time, qubits, det_props):
        pass

    def handle_detector_coords(detector_index: int, time, qubits, det_props):
        node = DecodingGraphNode(index=detector_index, time=time, qubits=qubits)
        node.properties = det_props
        index_to_DecodingGraphNode[detector_index] = node
        g.add_node(node)

    hyperedges = []

    iter_flatten_model(
        model,
        handle_error=skip_error,
        handle_detector_coords=handle_detector_coords,
        detectors=detectors,
        hyperedges=hyperedges,
    )

    trivial_boundary_node = DecodingGraphNode(index=model.num_detectors, time=0, is_boundary=True)
    g.add_node(trivial_boundary_node)
    index_to_DecodingGraphNode[model.num_detectors] = trivial_boundary_node

    iter_flatten_model(
        model,
        handle_error=handle_error,
        handle_detector_coords=skip_detector_coords,
        detectors=detectors,
        hyperedges=hyperedges,
    )

    return g, hyperedges


def string2nodes_with_detectors(
    string: str,
    detectors: List[Dict],
    logicals: List[Dict],
    clbits: QuantumCircuit.clbits,
    det_ref_values: Union[List, int] = 0,
    **kwargs,
):
    """
    Convert output string from circuits into a set of nodes for
    `DecodingGraph`.
    Args:
        string (string): Results string to convert.
        detectors: A list of measurement comparisons. A measurement comparison
                (detector) is either a list of measurements given by a the name and index
                of the classical bit or a list of dictionaries, with a mandatory clbits
                key containing the classical bits. A dictionary can contain keys like
                'qubits', 'time', 'basis' etc.
        logicals: A list of logical measurements. A logical measurement is a
                list of classical bits whose total parity is the logical eigenvalue.
                Again it can be a list of dictionaries.
        clbits: classical bits of the qiskit circuit, needed to identify
                measurements in the output string
        det_ref_values: Reference value for the detector outcomes, 0 by default

        kwargs (dict): Any additional keyword arguments.
            logical (str): Logical value whose results are used ('0' as default).
            all_logicals (bool): Whether to include logical nodes
            irrespective of value. (False as default).
    """

    output_bits = np.array([int(char) for char in string.replace(" ", "")[::-1]])

    clbit_dict = {(clbit._register.name, clbit._index): clind for clind, clbit in enumerate(clbits)}

    if isinstance(det_ref_values, int):
        det_ref_values = [det_ref_values] * len(detectors)

    nodes = []
    for ind, det in enumerate(detectors):
        det = det.copy()
        outcomes = [clbit_dict[clbit_key] for clbit_key in det.pop("clbits")]
        if sum(output_bits[outcomes]) % 2 != det_ref_values[ind]:
            node = DecodingGraphNode(time=det.pop("time"), qubits=det.pop("qubits"), index=ind)
            node.properties = det
            nodes.append(node)

    log_nodes = string2rawlogicals_with_detectors(
        string=string, logicals=logicals, clbits=clbits, start_ind=len(detectors), **kwargs
    )

    for node in log_nodes:
        nodes.append(node)

    return nodes


def string2rawlogicals_with_detectors(
    string: str,
    logicals: List[Dict],
    clbits: QuantumCircuit.clbits,
    start_ind: int = 0,
    **kwargs,
):
    """
    Convert output string from circuits into raw logical values.
    """

    all_logicals = kwargs.get("all_logicals")
    logical = kwargs.get("logical")
    if logical is None:
        logical = "0"

    output_bits = np.array([int(char) for char in string.replace(" ", "")[::-1]])

    clbit_dict = {(clbit._register.name, clbit._index): clind for clind, clbit in enumerate(clbits)}

    nodes = []
    for index, logical_op in enumerate(logicals, start=start_ind):
        logical_out = 0
        for q in logical_op["clbits"]:
            qind = clbit_dict[q]
            logical_out += output_bits[qind]
        logical_out = logical_out % 2

        if all_logicals or str(logical_out) != logical:
            node = DecodingGraphNode(
                is_boundary=True,
                qubits=[],
                index=index,
            )
            nodes.append(node)

    return nodes


def string2logical_meas(
    string: str,
    outcomes_in_logical: List[Dict],
    clbits: QuantumCircuit.clbits,
):
    """
    Args:
        string (string): Results string from qiskit circuit
        outcomes_in_logical: the detector-style logical outcome
        clbits: classical bits of the qiskit circuit, needed to identify
        measurements in the output string
    """

    output_bits = np.array([int(char) for char in string.replace(" ", "")[::-1]])

    clbit_dict = {(clbit._register.name, clbit._index): clind for clind, clbit in enumerate(clbits)}

    log_outs = []
    for logical_op in outcomes_in_logical:
        logical_out = 0
        for q in logical_op["clbits"]:
            qind = clbit_dict[q]
            logical_out += output_bits[qind]
        logical_out = logical_out % 2
        log_outs.append(logical_out)

    return log_outs


def noisify_circuit(circuits: Union[List, QuantumCircuit], noise_model: PauliNoiseModel):
    """
    Inserts error operations into a circuit according to a pauli noise model.
    Handles idling errors in the form of custom gates "idle_#" which are assumed to
    encode the identity gate only.
    qc = QuantumCircuit(1, name='idle_1')
    qc.i(0)
    idle_1 = qc.to_instruction()

    Args:
        circuits: Circuit or list thereof to which noise is added.
        noise_model: Pauli noise model used to define types of errors to add to circuit.

    Returns:
        noisy_circuits: Corresponding circuit or list thereof.
    """

    single_circuit = isinstance(circuits, QuantumCircuit)
    if single_circuit:
        circuits = [circuits]

    # create pauli errors for all errors in noise model
    errors = {}
    for g, noise in noise_model.to_dict().items():
        paulis = [pauli.upper() for pauli in noise["chan"].keys()]
        probs = list(noise["chan"].values())
        errors[g] = pauli_error(list(zip(paulis, probs)))

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
                noisy_qc.append(errors[g], qubits)
            # add gate if it needs to go after the error
            if not pre_error:
                if not g.startswith("idle_"):
                    noisy_qc.append(gate)

        noisy_circuits.append(noisy_qc)

    if single_circuit:
        noisy_circuits = noisy_circuits[0]

    return noisy_circuits
