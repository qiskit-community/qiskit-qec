"""Module of utilities for analyzing fault-tolerant circuits."""

from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag
from qiskit.circuit.library import XGate, YGate, ZGate
from qiskit import execute, Aer
from itertools import combinations, combinations_with_replacement

import cpp_helper  # lives in the .gitignore. Andrew Cross custom code


def pauli_fault_events(circ, order, location_types=[], pauli_types="all"):
    """Generator for iterating over faulty circuits.

    Insert all paths of exactly "order" faults into "circ".
    Operations whose names are in "location_types" are considered
    faulty. If location_types is empty, then all operations are faulty.
    A faulty operation has Pauli errors inserted after it (or
    before in the case of measurement). The pauli_types string
    indicates the types of Pauli errors that are inserted.

    circ = QuantumCircuit to insert faults into
    order = number of faults, positive integer
    location_types = list of operation names
    pauli_types = string that indicates the type of Pauli errors

    Return a tuple with data about the combination and the faulty circuit.
    """
    dag = circuit_to_dag(circ)
    counter = 0
    # TODO: take location_types into account
    for comb in combinations(
        filter(lambda x: x.name != "barrier", dag.topological_op_nodes()), order
    ):
        # Count the total number of qubits involved in this configuration
        total_qubits = 0
        for node in comb:
            total_qubits += len(node.qargs)
        # Iterate over non-identity Pauli errors on total_qubits
        # TODO: take pauli_types into account
        for pauli_tup in combinations_with_replacement("IXYZ", total_qubits):
            pauli = "".join(pauli_tup)
            if pauli == "I" * total_qubits:
                continue
            # Iterate over nodes in the combination
            start = 0
            circuit = QuantumCircuit(
                *dag.qregs.values(),
                *dag.cregs.values(),
                name=dag.name,
                global_phase=dag.global_phase
            )
            circuit.calibrations = dag.calibrations
            # Iterate over the original DAGCircuit
            for orig_node in dag.topological_op_nodes():
                # If the node is in the combination, replace it
                if orig_node in comb:
                    m = len(orig_node.qargs)
                    if orig_node.name == "measure":
                        for i in range(m):
                            if pauli[start + i] == "X":
                                circuit._append(
                                    XGate(), [orig_node.qargs[i]], orig_node.cargs
                                )
                            elif pauli[start + i] == "Y":
                                circuit._append(
                                    YGate(), [orig_node.qargs[i]], orig_node.cargs
                                )
                            elif pauli[start + i] == "Z":
                                circuit._append(
                                    ZGate(), [orig_node.qargs[i]], orig_node.cargs
                                )
                    inst = orig_node.op.copy()
                    inst.condition = orig_node.condition
                    circuit._append(inst, orig_node.qargs, orig_node.cargs)
                    if orig_node.name != "measure":
                        for i in range(m):
                            if pauli[start + i] == "X":
                                circuit._append(
                                    XGate(), [orig_node.qargs[i]], orig_node.cargs
                                )
                            elif pauli[start + i] == "Y":
                                circuit._append(
                                    YGate(), [orig_node.qargs[i]], orig_node.cargs
                                )
                            elif pauli[start + i] == "Z":
                                circuit._append(
                                    ZGate(), [orig_node.qargs[i]], orig_node.cargs
                                )
                    start += m
                # Otherwise append it to the circuit
                else:
                    inst = orig_node.op.copy()
                    inst.condition = orig_node.condition
                    circuit._append(inst, orig_node.qargs, orig_node.cargs)
            circuit.duration = dag.duration
            circuit.unit = dag.unit
            combstr = tuple(map(lambda x: x.name, comb))
            result = execute(
                circuit,
                Aer.get_backend("qasm_simulator"),
                method="stabilizer",
                shots=1,
                optimization_level=0,
                seed_simulator=0,
            ).result()
            outcomes = result.get_counts(circuit)
            raw_outcome = list(outcomes.keys())[0]
            outcome = list(map(int, raw_outcome[::-1]))
            yield (counter, combstr, pauli, outcome)
            counter += 1


def pauli_fault_events_error_propagation(
    circ, order, location_types=[], pauli_types="all"
):
    """Generator for iterating over faulty circuits.

    Insert all paths of exactly "order" faults into "circ".
    Operations whose names are in "location_types" are considered
    faulty. If location_types is empty, then all operations are faulty.
    A faulty operation has Pauli errors inserted after it (or
    before in the case of measurement). The pauli_types string
    indicates the types of Pauli errors that are inserted.

    This is a couple thousand times faster than pauli_fault_events.

    circ = QuantumCircuit to insert faults into
    order = number of faults, positive integer
    location_types = list of operation names
    pauli_types = string that indicates the type of Pauli errors

    Return a tuple with data about the combination and the error propagation
    outcomes from the faulty circuits.
    """
    # TODO generalize to more than one reg later
    if circ.qubits[0].register.size != len(circ.qubits):
        raise Exception("expected only one QuantumRegister")
    if circ.clbits[0].register.size != len(circ.clbits):
        raise Exception("expected only one ClassicalRegister")
    # TODO: enlarge this set of operations
    stabilizer_op_names = ["cx", "id", "h", "reset", "measure", "barrier", "x", "z"]
    name_to_index = {"cx": 0, "id": 1, "h": 2, "reset": 3, "measure": 4, "x": 6, "z": 7}
    # Count non-barrier stabilizer ops and record their data
    total_ops = 0
    non_barrier_ops = []
    for inst in circ.data:
        iname = inst[0].name
        iqubits = inst[1]
        iclbits = inst[2]
        if iname != "barrier":
            if iname not in stabilizer_op_names:
                raise Exception('op "%s" not recognized stabilizer op' % iname)
            total_ops += 1
            # Re-encode in a trivial way
            inum = name_to_index[iname]  # instruction number
            num_q = len(iqubits)  # num. qubits
            q_idx = [iqubits[j].index for j in range(len(iqubits))]  # indices
            c_idx = [iclbits[j].index for j in range(len(iclbits))]  # indices
            non_barrier_ops.append((inum, num_q, q_idx, c_idx))
    # Set up sparse matrix for storing measurement outcomes
    qreg_size = len(circ.qubits)
    creg_size = len(circ.clbits)

    # Store qreg_size X errors followed by qreg_size Z errors
    qubit_shape = (15 * total_ops, 2 * qreg_size)  # over-estimate rows
    clbit_shape = (15 * total_ops, creg_size)  # over-estimate rows
    qubit_mat = [0] * (qubit_shape[0] * qubit_shape[1])
    clbit_mat = [0] * (clbit_shape[0] * clbit_shape[1])

    # TODO: take location types into account
    # Iterate over combinations of instructions
    counter = 0
    for comb in combinations(range(len(non_barrier_ops)), order):
        # Count the total number of qubits involved in this combination
        total_qubits = sum([non_barrier_ops[i][1] for i in comb])
        # Iterate over non-identity Pauli errors on total_qubits
        # TODO: take pauli_types into account
        for pauli_tup in combinations_with_replacement("IXYZ", total_qubits):
            pauli = "".join(pauli_tup)
            if pauli == "I" * total_qubits:
                continue
            start = 0

            # Iterate over the circuit
            for j in range(len(non_barrier_ops)):
                dat = non_barrier_ops[j]
                # Simulate error propagation through operation
                # TODO ignored conditionals for now
                if dat[0] == 0:  # cx
                    q0 = dat[2][0]  # control
                    q1 = dat[2][1]  # target
                    # bit flips from control to target
                    qubit_mat[counter * 2 * qreg_size + q1] ^= qubit_mat[
                        counter * 2 * qreg_size + q0
                    ]
                    # phase flips from target to control
                    qubit_mat[counter * 2 * qreg_size + q0 + qreg_size] ^= qubit_mat[
                        counter * 2 * qreg_size + q1 + qreg_size
                    ]
                elif dat[0] == 1:  # id
                    pass
                elif dat[0] == 2:  # h
                    q = dat[2][0]
                    z = qubit_mat[counter * 2 * qreg_size + q + qreg_size]
                    qubit_mat[counter * 2 * qreg_size + q + qreg_size] = qubit_mat[
                        counter * 2 * qreg_size + q
                    ]
                    qubit_mat[counter * 2 * qreg_size + q] = z
                elif dat[0] == 3:  # reset
                    q = dat[2][0]
                    qubit_mat[counter * 2 * qreg_size + q] = 0
                    qubit_mat[counter * 2 * qreg_size + q + qreg_size] = 0
                elif dat[0] == 4:  # measure
                    q = dat[2][0]
                    c = dat[3][0]
                    # If the inst is in the combination and it's a measurement
                    # apply the error before the measurement
                    if j in comb:
                        if pauli[start] == "X" or pauli[start] == "Y":
                            qubit_mat[counter * 2 * qreg_size + q] ^= 1
                        start += 1
                    clbit_mat[counter * creg_size + c] = qubit_mat[
                        counter * 2 * qreg_size + q
                    ]
                elif dat[0] == 6:  # x (nop)
                    pass
                elif dat[0] == 7:  # z (nop)
                    pass
                # If the inst is in the combination and it's not a
                # measurement, apply the error after the operation
                if j in comb and dat[0] != 4:
                    for i in range(dat[1]):
                        q = dat[2][i]  # ith qubit
                        if pauli[start + i] == "X" or pauli[start + i] == "Y":
                            qubit_mat[counter * 2 * qreg_size + q] ^= 1
                        if pauli[start + i] == "Z" or pauli[start + i] == "Y":
                            qubit_mat[counter * 2 * qreg_size + q + qreg_size] ^= 1
                    start += dat[1]
            if counter == 0:
                raw_outcome = "".join(
                    map(
                        lambda x: str(x),
                        list(clbit_mat[counter * creg_size + creg_size - 1 :: -1]),
                    )
                )
            else:
                raw_outcome = "".join(
                    map(
                        lambda x: str(x),
                        list(
                            clbit_mat[
                                counter * creg_size
                                + creg_size
                                - 1 : counter * creg_size
                                - 1 : -1
                            ]
                        ),
                    )
                )
            outcome = list(map(int, raw_outcome[::-1]))

            combstr = tuple(
                map(lambda x: stabilizer_op_names[non_barrier_ops[x][0]], comb)
            )
            yield (counter, combstr, pauli, outcome)
            counter += 1


def pauli_fault_events_error_propagation_cpp(
    circ, order, location_types=[], pauli_types="all"
):
    """Generator for iterating over faulty circuits.

    Insert all paths of exactly "order" faults into "circ".
    Operations whose names are in "location_types" are considered
    faulty. If location_types is empty, then all operations are faulty.
    A faulty operation has Pauli errors inserted after it (or
    before in the case of measurement). The pauli_types string
    indicates the types of Pauli errors that are inserted.

    This implementation calls a cpp function to propagate errors.

    circ = QuantumCircuit to insert faults into
    order = number of faults, positive integer
    location_types = list of operation names
    pauli_types = string that indicates the type of Pauli errors

    Return a tuple with data about the combination and the error propagation
    outcomes from the faulty circuits.
    """
    # TODO generalize to more than one reg later
    if circ.qubits[0].register.size != len(circ.qubits):
        raise Exception("expected only one QuantumRegister")
    if circ.clbits[0].register.size != len(circ.clbits):
        raise Exception("expected only one ClassicalRegister")
    # TODO: enlarge this set of operations
    stabilizer_op_names = ["cx", "id", "h", "reset", "measure", "barrier", "x", "z"]
    name_to_index = {"cx": 0, "id": 1, "h": 2, "reset": 3, "measure": 4, "x": 6, "z": 7}
    # Count non-barrier stabilizer ops and record their data
    # in two different ways for later use
    total_ops = 0
    op_data_list = []
    non_barrier_ops = []
    for inst in circ.data:
        iname = inst[0].name
        iqubits = inst[1]
        iclbits = inst[2]
        if iname != "barrier":
            if iname not in stabilizer_op_names:
                raise Exception('op "%s" not recognized stabilizer op' % iname)
            total_ops += 1
            # Re-encode to a flat array of integers
            inum = name_to_index[iname]  # instruction number
            num_q = len(iqubits)  # num. qubits
            q_idx = [iqubits[j].index for j in range(len(iqubits))]  # indices
            c_idx = [iclbits[j].index for j in range(len(iclbits))]  # indices
            op_data_list.append(inum)
            for j in range(num_q):
                op_data_list.append(iqubits[j].index)
            if iname == "measure":
                op_data_list.append(iclbits[0].index)
            non_barrier_ops.append((inum, num_q, q_idx, c_idx))
    # Set up sparse matrix for storing measurement outcomes
    qreg_size = len(circ.qubits)
    creg_size = len(circ.clbits)

    # TODO: take location types into account
    # Iterate over combinations of instructions
    counter = 0
    for comb in combinations(range(total_ops), order):
        # Count the total number of qubits involved in this combination
        total_qubits = sum([non_barrier_ops[i][1] for i in comb])
        # Iterate over non-identity Pauli errors on total_qubits
        # TODO: take pauli_types into account
        for pauli_tup in combinations_with_replacement("IXYZ", total_qubits):
            pauli = "".join(pauli_tup)
            if pauli == "I" * total_qubits:
                continue
            outcome = cpp_helper.cpp_error_propagation(
                op_data_list, comb, pauli, total_ops, qreg_size, creg_size
            )
            combstr = tuple(
                map(lambda x: stabilizer_op_names[non_barrier_ops[x][0]], comb)
            )
            yield (counter, combstr, pauli, outcome)
            counter += 1
