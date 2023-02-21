"""Test the compiled error propagator directly."""
import unittest
from qiskit_qec.analysis.extensions import _CErrorPropagator


class Testcompilederrorpropagator(unittest.TestCase):
    """Tests for the compiled error propagator."""

    def test_errorpropagator_interact(self):
        """Test interaction with error propagator."""
        # pylint: disable=c-extension-no-member
        ep = _CErrorPropagator(3, 3)
        ep.apply_error([0, 2], "xy")
        self.assertEqual(ep.get_qubits(), "xiy")
        self.assertEqual(ep.get_qubit_array(), [1, 0, 1, 0, 0, 1])
        ep.cx(1, 2)
        self.assertEqual(ep.get_qubits(), "xzy")
        ep.measure(0, 0)
        ep.measure(1, 1)
        ep.measure(2, 2)
        self.assertEqual(ep.get_cbits(), [1, 0, 1])

    def test_errorpropagator_propagate(self):
        """Test loading circuit and propagating."""
        # pylint: disable=c-extension-no-member
        ep = _CErrorPropagator(3, 4)
        # H 0, CX 0 1, CX 1 2, M 0 0, M 1 1, M 2 2
        circ = [[0, 0], [5, 0, 1], [5, 1, 2], [8, 0, 0], [8, 1, 1], [8, 2, 2]]
        ep.load_circuit(3, 4, circ)
        self.assertEqual(ep.get_qreg_size(), 3)
        self.assertEqual(ep.get_creg_size(), 4)
        self.assertEqual(ep.get_circuit_size(), 6)
        # H fails with 'x'
        self.assertEqual(ep.propagate([0], ["x"]), [1, 1, 1, 0])


if __name__ == "__main__":
    unittest.main()
