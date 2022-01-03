Module errorpropagator
======================
Pauli error propagator.

Classes
-------

`ErrorPropagator(qreg_size=1, creg_size=1)`
:   Circuit error propagator interface.
    
    Create new error propagator.

    ### Ancestors (in MRO)

    * abc.ABC

    ### Methods

    `apply_error(self, q_idx, err_str)`
    :   Apply a single-qubit Pauli error during error propagation.
        
        q_idx = list of qubits the gate acts on
        err_str = string of "ixyz" characters describing Pauli error
        
        Method acts on qubit_array reference.

    `cx(self, qc, qt)`
    :   Apply CX gate.

    `get_bit_array(self)`
    :   Return the classical bit array.

    `get_error(self)`
    :   Return the qubit error state as a lowercase string.

    `get_qubit_array(self)`
    :   Return the qubit array.

    `h(self, q)`
    :   Apply Hadamard gate.

    `load_circuit(self, circ)`
    :   Express (stabilizer) circuit operations as a list of opcodes.
        
        circ = QuantumCircuit
        
        Encoded circuit is a list of tuples (opcode, qubits, clbits, label).
        The operations are visited in topologically sorted order.
        Return tuple: encoded circuit, qreg size, creg size.

    `measure(self, q, c)`
    :   Apply measure operation.
        
        Returns the outcome bit.

    `propagate_faults(self, icomb, error)`
    :   Insert a set of faults and propagate through a circuit.
        
        icomb = integer tuple of failed operations' indices
        error = tuple of pauli strings
        
        Return: measurement outcome discrepancies.

    `reset(self, q)`
    :   Apply reset operation.

    `s(self, q)`
    :   Apply Phase gate.