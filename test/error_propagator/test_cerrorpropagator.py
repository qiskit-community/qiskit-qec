"""Test C error propagator."""

import unittest
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_qec.analysis.cerrorpropagator import CErrorPropagator


class TestCErrorPropagator(unittest.TestCase):
    """Unit tests for C error propagator."""

    def test_error_propagation(self):
        """Test basic error propagator methods."""
        ep = CErrorPropagator(3, 3)
        ep.apply_error([0], "y")
        ep.h(0)
        ep.cx(0, 1)
        ep.cx(1, 2)
        ep.measure(0, 0)
        ep.measure(1, 1)
        ep.measure(2, 2)
        self.assertEqual(ep.get_qubit_array(), [1, 1, 1, 1, 0, 0])
        self.assertEqual(ep.get_bit_array(), [1, 1, 1])
        self.assertEqual(ep.get_error(), "yxx")

    def test_load_circuit(self):
        """Test loading a circuit and propagating errors."""
        qc = QuantumCircuit(1, 1)
        qc.h(0)
        qc.measure(0, 0)
        ep = CErrorPropagator()
        ep.load_circuit(qc)
        outcome = ep.propagate_faults([0], ["x"])
        self.assertEqual(outcome, [1])
        outcome = ep.propagate_faults([0], ["y"])
        self.assertEqual(outcome, [1])
        outcome = ep.propagate_faults([0], ["z"])
        self.assertEqual(outcome, [0])
        outcome = ep.propagate_faults([1], ["x"])
        self.assertEqual(outcome, [1])

    def test_load_another_circuit(self):
        """Test loading another circuit."""
        qc = QuantumCircuit(3, 3)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(1, 2)
        qc.id(0)
        qc.measure(0, 0)
        qc.measure(1, 1)
        qc.measure(2, 2)
        ep = CErrorPropagator()
        ep.load_circuit(qc)
        outcome = ep.propagate_faults([0], ["x"])
        self.assertEqual(outcome, [1, 1, 1])
        outcome = ep.propagate_faults([0], ["y"])
        self.assertEqual(outcome, [1, 1, 1])
        outcome = ep.propagate_faults([0], ["z"])
        self.assertEqual(outcome, [0, 0, 0])
        outcome = ep.propagate_faults([1], ["ix"])
        self.assertEqual(outcome, [0, 1, 1])

    def test_load_multiregister_circuit(self):
        """Test loading another circuit."""
        qrega = QuantumRegister(2, "qra")
        qregb = QuantumRegister(1, "qrb")
        crega = ClassicalRegister(3)
        qc = QuantumCircuit(qrega, qregb, crega)
        qc.h(qrega[0])
        qc.cx(qrega[0], qrega[1])
        qc.cx(qrega[1], qregb[0])
        qc.id(qrega[0])
        qc.measure(qrega[0], crega[0])
        qc.measure(qrega[1], crega[1])
        qc.measure(qregb[0], crega[2])
        ep = CErrorPropagator()
        ep.load_circuit(qc)
        outcome = ep.propagate_faults([0], ["x"])
        self.assertEqual(outcome, [1, 1, 1])
        outcome = ep.propagate_faults([0], ["y"])
        self.assertEqual(outcome, [1, 1, 1])
        outcome = ep.propagate_faults([0], ["z"])
        self.assertEqual(outcome, [0, 0, 0])
        outcome = ep.propagate_faults([1], ["ix"])
        self.assertEqual(outcome, [0, 1, 1])


if __name__ == "__main__":
    unittest.main()
