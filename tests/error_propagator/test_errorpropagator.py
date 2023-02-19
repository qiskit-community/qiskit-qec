"""Test python error propagator."""
import unittest
from qiskit_qec.analysis.errorpropagator import ErrorPropagator


class TestErrorPropagator(unittest.TestCase):
    """Pure python error propagator selection test."""

    def test_get_error_propagator_auto(self):
        """Test interaction with auto error propagator."""
        ep = ErrorPropagator("auto", 3, 3)  # pylint: disable=abstract-class-instantiated
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
        ep = ErrorPropagator("c", 3, 3)  # pylint: disable=abstract-class-instantiated
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
        ep = ErrorPropagator("c", 3, 3)  # pylint: disable=abstract-class-instantiated
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
