"""Test the compiled fault sampler."""
import unittest
from qiskit_qec.analysis.extensions import _CFaultSampler


class TestCompiledFaultSampler(unittest.TestCase):
    """Tests of the compiled fault sampler."""

    def test_faultsampler(self):
        """Test sampling multiple faults in a circuit."""
        # R 0, R 1, R 2, H 0, CX 0 1, I 2, CX 1 2, I 0, M 0 0, M 1 1, M 2 2
        circ = [
            [7, 0],
            [7, 1],
            [7, 2],
            [0, 0],
            [5, 0, 1],
            [6, 2],
            [5, 1, 2],
            [6, 0],
            [8, 0, 0],
            [8, 1, 1],
            [8, 2, 2],
        ]
        # H to second I
        faulty_ops_indices = [3, 4, 5, 6, 7]
        faulty_ops_labels = ["h", "cx", "id", "cx2", "id"]
        label_to_pauli_weight = {
            "h": [("x", 1), ("y", 1), ("z", 1)],
            "cx": [("ix", 1), ("xi", 1), ("xx", 1)],
            "id": [("x", 1), ("y", 1)],
            "cx2": [("ix", 1), ("xi", 1), ("xx", 1), ("iz", 1), ("zi", 1), ("zz", 1)],
        }
        label_to_error_probability = {"h": 0.1, "cx": 0.1, "id": 0.1, "cx2": 0.1}
        # pylint: disable=c-extension-no-member
        fe = _CFaultSampler(
            3,
            3,
            circ,
            faulty_ops_indices,
            faulty_ops_labels,
            label_to_pauli_weight,
            label_to_error_probability,
            10,
        )
        block = fe.sample(100000)
        num_faults = {}
        for x in block:
            failures = len(x[1])
            if failures not in num_faults:
                num_faults[failures] = 1
            else:
                num_faults[failures] += 1
        for key in num_faults:
            num_faults[key] /= 100000
        limit = {
            0: (1 - 0.1) ** 5,
            1: 5 * 0.1 * (1 - 0.1) ** 4,
            2: 10 * 0.1**2 * (1 - 0.1) ** 3,
            3: 10 * 0.1**3 * (1 - 0.1) ** 2,
            4: 5 * 0.1**4 * (1 - 0.1),
            5: 0.1**5,
        }
        for key, value in num_faults.items():
            self.assertAlmostEqual(value, limit[key], 2)

    def test_faultsampler_2(self):
        """Test sampling in a circuit with one faulty location."""
        # R 0, H 0, M 0 0
        circ = [[7, 0], [0, 0], [8, 0, 0]]
        faulty_ops_indices = [1]
        faulty_ops_labels = ["h"]
        label_to_pauli_weight = {"h": [("x", 1), ("y", 1), ("z", 1)]}
        label_to_error_probability = {"h": 0.2}
        # pylint: disable=c-extension-no-member
        fe = _CFaultSampler(
            1,
            1,
            circ,
            faulty_ops_indices,
            faulty_ops_labels,
            label_to_pauli_weight,
            label_to_error_probability,
            10,
        )
        block = fe.sample(100000)
        failures = 0
        for x in block:
            if x[3][0] == 1:
                failures += 1
        failures /= 100000
        self.assertAlmostEqual(failures, 0.2 * 2.0 / 3.0, 2)


if __name__ == "__main__":
    unittest.main()
