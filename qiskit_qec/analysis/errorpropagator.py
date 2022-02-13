"""Pauli error propagator."""

from abc import ABC, abstractmethod


class ErrorPropagator(ABC):
    """Circuit error propagator interface."""

    @abstractmethod
    def __init__(self, qreg_size=1, creg_size=1):
        """Create new error propagator."""
        pass

    @abstractmethod
    def apply_error(self, q_idx, err_str):
        """Apply a single-qubit Pauli error during error propagation.

        q_idx = list of qubits the gate acts on
        err_str = string of "ixyz" characters describing Pauli error

        Method acts on qubit_array reference.
        """
        pass

    @abstractmethod
    def load_circuit(self, circ):
        """Express (stabilizer) circuit operations as a list of opcodes.

        circ = QuantumCircuit

        Encoded circuit is a list of tuples (opcode, qubits, clbits, label).
        The operations are visited in topologically sorted order.
        Return tuple: encoded circuit, qreg size, creg size.
        """
        pass

    @abstractmethod
    def cx(self, qc, qt):
        """Apply CX gate."""
        pass

    @abstractmethod
    def h(self, q):
        """Apply Hadamard gate."""
        pass

    @abstractmethod
    def s(self, q):
        """Apply Phase gate."""
        pass

    @abstractmethod
    def reset(self, q):
        """Apply reset operation."""
        pass

    @abstractmethod
    def measure(self, q, c):
        """Apply measure operation.

        Returns the outcome bit.
        """
        pass

    @abstractmethod
    def propagate_faults(self, icomb, error):
        """Insert a set of faults and propagate through a circuit.

        icomb = integer tuple of failed operations' indices
        error = tuple of pauli strings

        Return: measurement outcome discrepancies.
        """
        pass

    @abstractmethod
    def get_qubit_array(self):
        """Return the qubit array."""
        pass

    @abstractmethod
    def get_bit_array(self):
        """Return the classical bit array."""
        pass

    @abstractmethod
    def get_error(self):
        """Return the qubit error state as a lowercase string."""
        pass
