"""Test python error propagator selector."""
import unittest
from qiskit_qec.analysis.epselector import EPSelector


class TestErrorPropagatorSelector(unittest.TestCase):
    """Pure python error propagator selection test."""

    def test_get_error_propagator_auto(self):
        """Test interaction with auto error propagator."""
        eps = EPSelector()
        ep = eps.get_error_propagator("auto", 3, 3)
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

    def test_get_error_propagator_c(self):
        """Test interaction with c error propagator."""
        eps = EPSelector()
        ep = eps.get_error_propagator("c", 3, 3)
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

    def test_get_error_propagator_py(self):
        """Test interaction with py error propagator."""
        eps = EPSelector()
        ep = eps.get_error_propagator("py", 3, 3)
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


if __name__ == "__main__":
    unittest.main()
