"""Compiled Pauli error propagator."""

from typing import List, Tuple

from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag

from qiskit_qec.analysis.extensions import C_ERROR_PROPAGATOR

if C_ERROR_PROPAGATOR:
    from qiskit_qec.analysis.baseerrorpropagator import BaseErrorPropagator

    from qiskit_qec.analysis.extensions import _CErrorPropagator

    class CErrorPropagator(BaseErrorPropagator):
        """ErrorPropagator object using compiledextensions."""

        def __init__(self, qreg_size: int = 1, creg_size: int = 1):
            """Create new error propagator."""
            # pylint: disable=super-init-not-called
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
            self.name_to_code = dict(zip(self.stabilizer_op_names, self.stabilizer_op_codes))
            self.encoded_circ = None
            # pylint: disable=c-extension-no-member
            self.cep = _CErrorPropagator(qreg_size, creg_size)

        def apply_error(self, q_idx: List[int], err_str: str):
            """Apply a single-qubit Pauli error during error propagation.

            q_idx = list of qubits the gate acts on
            err_str = string of "ixyz" characters describing Pauli error

            Method acts on qubit_array reference.
            """
            assert len(q_idx) == len(err_str), "length mismatch"
            assert set(err_str) <= set("ixyz"), "bad error string"
            self.cep.apply_error(q_idx, err_str)

        def load_circuit(self, circ: QuantumCircuit):
            """Express (stabilizer) circuit operations as a list of opcodes.

            circ = QuantumCircuit

            Encoded circuit is a list of tuples (opcode, qubits, clbits, label).
            The operations are visited in topologically sorted order.

            Return tuple: encoded circuit, qreg size, creg size.
            """
            dag = circuit_to_dag(circ)
            self.encoded_circ = []
            qubit_indices = {bit: index for index, bit in enumerate(circ.qubits)}
            clbit_indices = {bit: index for index, bit in enumerate(circ.clbits)}
            for node in dag.topological_op_nodes():
                name = node.name
                if name not in self.stabilizer_op_names:
                    raise Exception(f'op "{name}" not recognized')
                opcode = self.name_to_code[name]
                q_idx = [qubit_indices[qarg] for qarg in node.qargs]
                c_idx = [clbit_indices[carg] for carg in node.cargs]
                # TODO: conditionals currently ignored, fix later
                if name in ["h", "s", "x", "y", "z", "id", "reset"]:
                    self.encoded_circ.append([opcode, q_idx[0]])
                elif name == "barrier":
                    self.encoded_circ.append([opcode])
                elif name == "cx":
                    self.encoded_circ.append([opcode, q_idx[0], q_idx[1]])
                elif name == "measure":
                    self.encoded_circ.append([opcode, q_idx[0], c_idx[0]])
                else:
                    raise Exception("bad opcode")
            self.cep.load_circuit(len(circ.qubits), len(circ.clbits), self.encoded_circ)

        def cx(self, qc: int, qt: int):
            """Apply CX gate."""
            self.cep.cx(qc, qt)

        def h(self, q: int):
            """Apply Hadamard gate."""
            self.cep.h(q)

        def s(self, q: int):
            """Apply Phase gate."""
            self.cep.s(q)

        def reset(self, q: int):
            """Apply reset operation."""
            self.cep.reset(q)

        def measure(self, q: int, c: int):
            """Apply measure operation.

            Returns the outcome bit.
            """
            return self.cep.measure(q, c)

        def propagate_faults(self, icomb: Tuple[int], error: Tuple[str]):
            """Insert a set of faults and propagate through a circuit.

            icomb = integer tuple of failed operations' indices
            error = tuple of pauli strings

            Return: measurement outcome discrepancies.
            """
            if self.encoded_circ is None:
                raise Exception("no circuit loaded")
            return self.cep.propagate(icomb, error)

        def get_qubit_array(self):
            """Return the qubit array."""
            return self.cep.get_qubit_array()

        def get_bit_array(self):
            """Return the classical bit array."""
            return self.cep.get_cbits()

        def get_error(self):
            """Return the qubit error state as a lowercase string."""
            return self.cep.get_qubits()
