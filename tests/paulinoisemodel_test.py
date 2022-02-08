"""Pauli noise model tests."""
import unittest
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class TestPauliNoiseModel(unittest.TestCase):
    """Test Pauli noise model."""

    def test_from_dict(self):
        """Test creating from a dictionary."""
        def_as_dict = {
            "cx": {
                "p": 0.001,
                "chan": {
                    "ix": 1,
                    "iy": 1,
                    "iz": 1,
                    "xi": 1,
                    "xx": 1,
                    "xy": 1,
                    "xz": 1,
                    "yi": 1,
                    "yx": 1,
                    "yy": 1,
                    "yz": 1,
                    "zi": 1,
                    "zx": 1,
                    "zy": 1,
                    "zz": 1,
                },
            },
            "h": {"p": 0.0001, "chan": {"x": 1, "y": 1, "z": 1}},
            "id": {"p": 0.0001, "chan": {"x": 1, "y": 1, "z": 1}},
            "idm": {"p": 0.001, "chan": {"x": 1, "y": 1, "z": 1}},
            "reset": {"p": 0.01, "chan": {"x": 1}},
            "measure": {"p": 0.01, "chan": {"x": 1}},
            "x": {"p": 0.0001, "chan": {"x": 1, "y": 1, "z": 1}},
            "y": {"p": 0.0001, "chan": {"x": 1, "y": 1, "z": 1}},
            "z": {"p": 0.0001, "chan": {"x": 1, "y": 1, "z": 1}},
        }
        pnm = PauliNoiseModel(def_as_dict)
        self.assertEqual(pnm.get_error_probability("h"), 0.0001)
        self.assertAlmostEqual(pnm.get_pauli_weight("cx", "xx"), 1.0 / 15.0)

    def test_scale_factor(self):
        """Test setting a scale factor."""
        pnm = PauliNoiseModel()
        pnm.add_operation("h", {"x": 1, "y": 1, "z": 1})
        pnm.add_operation("s", {"x": 1, "y": 1, "z": 1})
        pnm.add_operation("y", {"x": 1, "y": 1, "z": 1})
        pnm.set_scale_factor("h", 1)
        pnm.set_scale_factor("s", 1 / 10)
        pnm.set_error_probability("y", 0.01)
        pnm.set_scaled_error_probabilities(0.001)
        self.assertAlmostEqual(pnm.get_error_probability("h"), 0.001)
        self.assertAlmostEqual(pnm.get_error_probability("s"), 0.0001)
        self.assertAlmostEqual(pnm.get_error_probability("y"), 0.01)
        pnm.set_scaled_error_probabilities(0.1)
        self.assertAlmostEqual(pnm.get_error_probability("h"), 0.1)
        self.assertAlmostEqual(pnm.get_error_probability("s"), 0.01)
        self.assertAlmostEqual(pnm.get_error_probability("y"), 0.01)

    def test_add_operation(self):
        """Test adding an operation."""
        pnm = PauliNoiseModel()
        pnm.add_operation("cx", {"ix": 1, "xi": 1, "xx": 1})
        pnm.set_error_probability("cx", 0.01)
        self.assertEqual(pnm.get_operations(), ["cx"])
        self.assertEqual(pnm.get_error_probability("cx"), 0.01)
        self.assertEqual(pnm.get_pauli_weight("cx", "iz"), 0)
        self.assertAlmostEqual(pnm.get_pauli_weight("cx", "xx"), 1.0 / 3.0)
        self.assertAlmostEqual(pnm.get_pauli_error_probability("cx", "xx"), 0.01 / 3.0)
        self.assertEqual(pnm.get_pauli_error_types(), {"cx": ["ix", "xi", "xx"]})

    def test_get_before_add(self):
        """Get a gate before adding it."""
        with self.assertRaises(Exception):
            pnm = PauliNoiseModel()
            pnm.get_error_probability("cx")

    def test_no_error_probability(self):
        """Get a probability before setting."""
        with self.assertRaises(Exception):
            pnm = PauliNoiseModel()
            pnm.add_operation("h", {"x": 1})
            pnm.get_error_probability("h")

    def test_bad_input_1(self):
        """Pass out of range value."""
        with self.assertRaises(Exception):
            pnm = PauliNoiseModel()
            pnm.add_operation("h", {"x": 1})
            pnm.set_error_probability("h", 5)

    def test_bad_input_2(self):
        """Pass out of range value."""
        with self.assertRaises(Exception):
            pnm = PauliNoiseModel()
            pnm.add_operation("h", {"x": 1})
            pnm.set_error_probability("h", -0.1)

    def test_bad_input_3(self):
        """Get something that doesn't exist."""
        with self.assertRaises(Exception):
            pnm = PauliNoiseModel()
            pnm.add_operation("h", {"x": 1})
            pnm.set_error_probability("h", 0.01)
            pnm.get_pauli_weight("cx", "ix")

    def test_bad_input_4(self):
        """Get the wrong term."""
        with self.assertRaises(Exception):
            pnm = PauliNoiseModel()
            pnm.add_operation("h", {"x": 1})
            pnm.set_error_probability("h", 0.01)
            pnm.get_pauli_weight("h", "ix")

    def test_bad_input_5(self):
        """Get a malformed term."""
        with self.assertRaises(Exception):
            pnm = PauliNoiseModel()
            pnm.add_operation("h", {"x": 1})
            pnm.set_error_probability("h", 0.01)
            pnm.get_pauli_weight("h", "fx")


if __name__ == "__main__":
    unittest.main()
