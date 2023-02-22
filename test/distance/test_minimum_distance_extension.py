"""Test minimum distance computations."""
import unittest
import numpy as np

from qiskit_qec.linear.symplectic import normalizer
from qiskit_qec.analysis.extensions import _c_minimum_distance


def strarray(label: str) -> np.ndarray:
    """Convert a Pauli label string into a symplectic vector."""
    n = len(label)
    vec = np.zeros((2 * n,))
    for i, pauli in enumerate(label):
        if pauli in ("x", "y"):
            vec[i] = 1
        if pauli in ("z", "y"):
            vec[i + n] = 1
    return vec


# pylint: disable=c-extension-no-member
class TestMinimumDistanceCompiled(unittest.TestCase):
    """Tests of minimum distance property in compiled extensions."""

    def test_minimum_distance_1(self):
        """Test Steane code."""
        paulis = ["iiixxxx", "ixxiixx", "xixixix", "iiizzzz", "izziizz", "ziziziz"]
        stabilizer = np.asarray(list(map(strarray, paulis)))
        inputform = stabilizer.astype(np.int32).tolist()
        d = _c_minimum_distance(inputform, inputform)
        self.assertEqual(d, 3)

    def test_minimum_distance_2(self):
        """Test four qubit code."""
        paulis = ["xxxx", "zzzz"]
        stabilizer = np.asarray(list(map(strarray, paulis)))
        inputform = stabilizer.astype(np.int32).tolist()
        d = _c_minimum_distance(inputform, inputform)
        self.assertEqual(d, 2)

    def test_minimum_distance_3(self):
        """Test five qubit code."""
        paulis = ["xzzxi", "ixzzx", "xixzz", "zxixz"]
        stabilizer = np.asarray(list(map(strarray, paulis)))
        inputform = stabilizer.astype(np.int32).tolist()
        d = _c_minimum_distance(inputform, inputform)
        self.assertEqual(d, 3)

    def test_minimum_distance_4(self):
        """Test [[10,2,4]] code from codetables.de."""
        paulis = [
            "xzizixizzi",
            "iyizzyiizz",
            "izxiiyiyzx",
            "izzyiyziiz",
            "iizzyyzxzx",
            "izzizzxxzz",
            "iziiizzzyy",
            "zzzzzziiii",
        ]
        stabilizer = np.asarray(list(map(strarray, paulis)))
        inputform = stabilizer.astype(np.int32).tolist()
        d = _c_minimum_distance(inputform, inputform)
        self.assertEqual(d, 4)

    def test_minimum_distance_4_overcomplete(self):
        """Test [[10,2,4]] code from codetables.de, overcomplete version."""
        paulis = [
            "xzizixizzi",
            "iyizzyiizz",
            "izxiiyiyzx",
            "izzyiyziiz",
            "iizzyyzxzx",
            "izzizzxxzz",
            "iziiizzzyy",
            "zzzzzziiii",
            "zizzzizzyy",
            "ziiziixxzz",
        ]
        stabilizer = np.asarray(list(map(strarray, paulis)))
        inputform = stabilizer.astype(np.int32).tolist()
        d = _c_minimum_distance(inputform, inputform)
        self.assertEqual(d, 4)

    def test_minimum_distance_subsystem(self):
        """Test [[9,1,4,3]] Bacon-Shor code."""
        # Overcomplete generating set of gauge group
        paulis = [
            "xiixiiiii",
            "ixiixiiii",
            "iixiixiii",
            "iiixiixii",
            "iiiixiixi",
            "iiiiixiix",
            "zziiiiiii",
            "izziiiiii",
            "iiizziiii",
            "iiiizziii",
            "iiiiiizzi",
            "iiiiiiizz",
        ]
        gauge = np.asarray(list(map(strarray, paulis)), dtype=bool)
        stabilizer, _, _ = normalizer(gauge)
        inputform1 = stabilizer.astype(np.int32).tolist()
        inputform2 = gauge.astype(np.int32).tolist()
        d = _c_minimum_distance(inputform1, inputform2)
        self.assertEqual(d, 3)

    def test_minimum_distance_subsystem_2(self):
        """Test [[9,0,4,2]] Bacon-Shor code state."""
        # Overcomplete generating set of gauge group
        paulis = [
            "xiixiiiii",
            "ixiixiiii",
            "iixiixiii",
            "iiixiixii",
            "iiiixiixi",
            "iiiiixiix",
            "zziiiiiii",
            "izziiiiii",
            "iiizziiii",
            "iiiizziii",
            "iiiiiizzi",
            "iiiiiiizz",
            "zzzzzzzzz",
        ]
        gauge = np.asarray(list(map(strarray, paulis)), dtype=bool)
        stabilizer, _, _ = normalizer(gauge)
        inputform1 = stabilizer.astype(np.int32).tolist()
        inputform2 = gauge.astype(np.int32).tolist()
        d = _c_minimum_distance(inputform1, inputform2)
        self.assertEqual(d, 2)


if __name__ == "__main__":
    unittest.main()
