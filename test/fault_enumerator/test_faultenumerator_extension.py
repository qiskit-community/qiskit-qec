"""Test the compiled fault enumerator."""
import unittest
from qiskit_qec.analysis.extensions import _CFaultEnumerator


class TestFaultEnumeratorExtension(unittest.TestCase):
    """Tests for compiled fault enumerator."""

    def test_faultenumerator(self):
        """Test example inputs/outputs for compiled fault enumerator."""
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
        faulty_ops_pauli_errors = [
            ["x", "y", "x"],
            ["ix", "xi", "xx"],
            ["x", "y"],
            ["ix", "xi", "xx", "iz", "zi", "zz"],
            ["x"],
        ]
        # pylint: disable=c-extension-no-member
        fe = _CFaultEnumerator(
            1,
            3,
            3,
            circ,
            faulty_ops_indices,
            faulty_ops_labels,
            faulty_ops_pauli_errors,
        )
        output = [
            (
                3,
                [0],
                False,
                [
                    (0, ["h"], ["x"], [1, 1, 1]),
                    (1, ["h"], ["y"], [1, 1, 1]),
                    (2, ["h"], ["x"], [1, 1, 1]),
                ],
            ),
            (
                6,
                [1],
                False,
                [
                    (3, ["cx"], ["ix"], [0, 1, 1]),
                    (4, ["cx"], ["xi"], [1, 0, 0]),
                    (5, ["cx"], ["xx"], [1, 1, 1]),
                ],
            ),
            (
                14,
                [3],
                False,
                [
                    (6, ["id"], ["x"], [0, 0, 1]),
                    (7, ["id"], ["y"], [0, 0, 1]),
                    (8, ["cx2"], ["ix"], [0, 0, 1]),
                    (9, ["cx2"], ["xi"], [0, 1, 0]),
                    (10, ["cx2"], ["xx"], [0, 1, 1]),
                    (11, ["cx2"], ["iz"], [0, 0, 0]),
                    (12, ["cx2"], ["zi"], [0, 0, 0]),
                    (13, ["cx2"], ["zz"], [0, 0, 0]),
                ],
            ),
            (15, [4], True, [(14, ["id"], ["x"], [1, 0, 0])]),
        ]

        step = 0
        while not fe.done():
            block = fe.enumerate(3)
            self.assertEqual(fe.get_index(), output[step][0])
            self.assertEqual(fe.get_state(), output[step][1])
            self.assertEqual(fe.done(), output[step][2])
            self.assertEqual(block, output[step][3])
            step += 1

        fe.reset()
        step = 0
        while not fe.done():
            block = fe.enumerate(3)
            self.assertEqual(fe.get_index(), output[step][0])
            self.assertEqual(fe.get_state(), output[step][1])
            self.assertEqual(fe.done(), output[step][2])
            self.assertEqual(block, output[step][3])
            step += 1


if __name__ == "__main__":
    unittest.main()
