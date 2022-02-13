"""Python Pauli error propagator."""

from qiskit.converters import circuit_to_dag

from qiskit_qec.analysis.errorpropagator import ErrorPropagator


class PyErrorPropagator(ErrorPropagator):
    """ErrorPropagator object in pure Python."""

    def __init__(self, qreg_size=1, creg_size=1):
        """Create new error propagator."""
        self.stabilizer_op_names = [
            "h",
            "s",
            "x",
            "y",
            "z",
            "cx",
            "id",
            "reset",
            "measure",
            "barrier",
        ]
        self.stabilizer_op_codes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        self.name_to_code = {
            k: v for k, v in zip(self.stabilizer_op_names, self.stabilizer_op_codes)
        }
        self.gate_dispatch = {
            0: self._faulty_h,
            1: self._faulty_s,
            2: self._faulty_id,
            3: self._faulty_id,
            4: self._faulty_id,
            5: self._faulty_cx,
            6: self._faulty_id,
            7: self._faulty_reset,
            8: self._faulty_measure,
            9: self._faulty_barrier,
        }
        self.encoded_circ = None
        self.qreg_size = qreg_size
        self.creg_size = creg_size
        # Store qreg_size X errors followed by qreg_size Z errors
        self.qubit_array = [0] * (2 * qreg_size)
        self.clbit_array = [0] * creg_size

    def _range_check(self, q):
        """Check that qubit index is in bounds.

        Not used, so be careful :)
        """
        if q >= self.qreg_size or q < 0:
            raise IndexError("qubit index out of range")

    def apply_error(self, q_idx, err_str):
        """Apply a single-qubit Pauli error during error propagation.

        q_idx = list of qubits the gate acts on
        err_str = string of "ixyz" characters describing Pauli error

        Method acts on qubit_array reference.
        """
        assert len(q_idx) == len(err_str), "length mismatch"
        assert set(err_str) <= set("ixyz"), "bad error string"
        for i in range(len(q_idx)):
            q = q_idx[i]
            # self._range_check(q)
            if err_str[i] == "x" or err_str[i] == "y":
                self.qubit_array[q] ^= 1
            if err_str[i] == "z" or err_str[i] == "y":
                self.qubit_array[q + self.qreg_size] ^= 1

    def _node_name_label(self, node):
        """Form an identifier string for a node's operation.

        Use node.op._label if it exists. Otherwise use node.name.
        Return a string.
        """
        if "_label" in node.op.__dict__ and node.op._label is not None:
            name_label = node.op._label
        else:
            name_label = node.name
        return name_label

    def load_circuit(self, circ):
        """Express (stabilizer) circuit operations as a list of opcodes.

        circ = QuantumCircuit

        Encoded circuit is a list of tuples (opcode, qubits, clbits, label).
        The operations are visited in topologically sorted order.
        Return tuple: encoded circuit, qreg size, creg size.
        """
        # We expect a single quantum and classical register
        if circ.qubits[0].register.size != len(circ.qubits):
            raise Exception("expected only one QuantumRegister")
        if circ.clbits[0].register.size != len(circ.clbits):
            raise Exception("expected only one ClassicalRegister")
        dag = circuit_to_dag(circ)
        self.encoded_circ = []
        for node in dag.topological_op_nodes():
            name = node.name
            label = self._node_name_label(node)
            qubits = node.qargs
            clbits = node.cargs
            if name not in self.stabilizer_op_names:
                raise Exception('op "%s" not recognized' % name)
            opcode = self.name_to_code[name]
            q_idx = [qubits[j].index for j in range(len(qubits))]
            c_idx = [clbits[j].index for j in range(len(clbits))]
            # TODO: conditionals currently ignored, fix later
            self.encoded_circ.append((opcode, q_idx, c_idx, label))
        self.qreg_size = len(circ.qubits)
        self.creg_size = len(circ.clbits)
        self.qubit_array = [0] * (2 * self.qreg_size)
        self.clbit_array = [0] * self.creg_size

    def cx(self, qc, qt):
        """Apply CX gate."""
        # self._range_check(qc)
        # self._range_check(qt)
        self.qubit_array[qt] ^= self.qubit_array[qc]
        self.qubit_array[qc + self.qreg_size] ^= self.qubit_array[qt + self.qreg_size]

    def h(self, q):
        """Apply Hadamard gate."""
        # self._range_check(q)
        z = self.qubit_array[q + self.qreg_size]
        self.qubit_array[q + self.qreg_size] = self.qubit_array[q]
        self.qubit_array[q] = z

    def s(self, q):
        """Apply Phase gate."""
        # self._range_check(q)
        if self.qubit_array[q]:
            self.qubit_array[q + self.qreg_size] ^= 1

    def reset(self, q):
        """Apply reset operation."""
        # self._range_check(q)
        self.qubit_array[q] = 0
        self.qubit_array[q + self.qreg_size] = 0

    def measure(self, q, c):
        """Apply measure operation.

        Returns the outcome bit.
        """
        # self._range_check(q)
        # if c >= self.creg_size or c < 0:
        #     raise IndexError("bit index out of range")
        self.clbit_array[c] = self.qubit_array[q]
        return self.clbit_array[c]

    def _faulty_cx(self, op_idx, q_idx, c_idx, icomb, error):
        """Apply faulty CX gate."""
        self.cx(q_idx[0], q_idx[1])
        if op_idx in icomb:  # i.e., this operation failed
            j = icomb.index(op_idx)
            self.apply_error(q_idx, error[j])

    def _faulty_id(self, op_idx, q_idx, c_idx, icomb, error):
        """Apply faulty Identity gate."""
        if op_idx in icomb:  # i.e., this operation failed
            j = icomb.index(op_idx)
            self.apply_error(q_idx, error[j])

    def _faulty_h(self, op_idx, q_idx, c_idx, icomb, error):
        """Apply faulty H gate."""
        self.h(q_idx[0])
        if op_idx in icomb:  # i.e., this operation failed
            j = icomb.index(op_idx)
            self.apply_error(q_idx, error[j])

    def _faulty_s(self, op_idx, q_idx, c_idx, icomb, error):
        """Apply faulty S gate."""
        self.s(q_idx[0])
        if op_idx in icomb:  # i.e., this operation failed
            j = icomb.index(op_idx)
            self.apply_error(q_idx, error[j])

    def _faulty_reset(self, op_idx, q_idx, c_idx, icomb, error):
        """Apply faulty reset operation."""
        self.reset(q_idx[0])
        if op_idx in icomb:  # i.e., this operation failed
            j = icomb.index(op_idx)
            self.apply_error(q_idx, error[j])

    def _faulty_measure(self, op_idx, q_idx, c_idx, icomb, error):
        """Apply faulty measure operation."""
        if op_idx in icomb:  # i.e., this operation failed
            j = icomb.index(op_idx)
            self.apply_error(q_idx, error[j])
        self.measure(q_idx[0], c_idx[0])  # don't return the value

    def _faulty_barrier(self, op_idx, q_idx, c_idx, icomb, error):
        """Apply faulty barrier operation."""
        pass  # ignore barriers

    def propagate_faults(self, icomb, error):
        """Insert a set of faults and propagate through a circuit.

        icomb = integer tuple of failed operations' indices
        error = tuple of pauli strings

        Return: measurement outcome discrepancies.
        """
        if self.encoded_circ is None:
            raise Exception("no circuit loaded")
        self.qubit_array = [0] * (2 * self.qreg_size)
        self.clbit_array = [0] * self.creg_size
        for j in range(len(self.encoded_circ)):
            opcode, q_idx, c_idx, label = self.encoded_circ[j]
            self.gate_dispatch[opcode](j, q_idx, c_idx, icomb, error)
        return self.clbit_array

    def get_qubit_array(self):
        """Return the qubit array."""
        return self.qubit_array

    def get_bit_array(self):
        """Return the classical bit array."""
        return self.clbit_array

    def get_error(self):
        """Return the qubit error state as a lowercase string."""
        error = ""
        for j in range(self.qreg_size):
            if self.qubit_array[j] == 0 and self.qubit_array[j + self.qreg_size] == 0:
                error += "i"
            if self.qubit_array[j] == 0 and self.qubit_array[j + self.qreg_size] == 1:
                error += "z"
            if self.qubit_array[j] == 1 and self.qubit_array[j + self.qreg_size] == 0:
                error += "x"
            if self.qubit_array[j] == 1 and self.qubit_array[j + self.qreg_size] == 1:
                error += "y"
        return error
