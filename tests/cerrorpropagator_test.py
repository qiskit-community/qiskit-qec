import unittest
from qiskit import QuantumCircuit
from qiskit_qec.analysis.cerrorpropagator import CErrorPropagator


class TestCErrorPropagator(unittest.TestCase):
    def test_error_propagation(self):
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


if __name__ == "__main__":
    unittest.main()
