"""Test heavy-hexagon code family definition."""
import unittest
from qiskit_qec.codes.hhc import HHC


class TestHHC(unittest.TestCase):
    """Test heavy-hexagon code family."""

    def test_even_distance(self):
        """Only odd distance implemented."""
        with self.assertRaises(Exception):
            HHC(4)

    def test_d3(self):
        """Check d=3 output."""
        c = HHC(3)
        self.assertEqual(c.n, 9)
        self.assertEqual(c.k, 1)
        self.assertEqual(c.d, 3)
        self.assertEqual(c.logical_x, [[0, 1, 2, 3, 4, 5, 6, 7, 8]])
        self.assertEqual(c.logical_z, [[0, 1, 2, 3, 4, 5, 6, 7, 8]])
        self.assertEqual(c.x_boundary, [[0], [1], [2], [6], [7], [8]])
        self.assertEqual(c.z_boundary, [[0], [3], [6], [2], [5], [8]])
        self.assertEqual(c.x_gauges, [[0, 3], [1, 4], [2, 5], [3, 6], [4, 7], [5, 8]])
        self.assertEqual(c.z_gauges, [[0, 1], [1, 2, 4, 5], [3, 4, 6, 7], [7, 8]])
        self.assertEqual(c.x_stabilizers, [[0, 1, 3, 4], [2, 5], [3, 6], [4, 5, 7, 8]])
        self.assertEqual(c.z_stabilizers, [[0, 1, 3, 4, 6, 7], [1, 2, 4, 5, 7, 8]])

    def test_d5(self):
        """Check d=5 output."""
        c = HHC(5)
        self.assertEqual(c.n, 25)
        self.assertEqual(c.k, 1)
        self.assertEqual(c.d, 5)
        self.assertEqual(
            c.logical_x,
            [
                [
                    0,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                ]
            ],
        )
        self.assertEqual(
            c.logical_z,
            [
                [
                    0,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                ]
            ],
        )
        self.assertEqual(c.x_boundary, [[0], [1], [2], [3], [4], [20], [21], [22], [23], [24]])
        self.assertEqual(c.z_boundary, [[0], [5], [10], [15], [20], [4], [9], [14], [19], [24]])
        self.assertEqual(
            c.x_gauges,
            [
                [0, 5],
                [1, 6],
                [2, 7],
                [3, 8],
                [4, 9],
                [5, 10],
                [6, 11],
                [7, 12],
                [8, 13],
                [9, 14],
                [10, 15],
                [11, 16],
                [12, 17],
                [13, 18],
                [14, 19],
                [15, 20],
                [16, 21],
                [17, 22],
                [18, 23],
                [19, 24],
            ],
        )
        self.assertEqual(
            c.z_gauges,
            [
                [0, 1],
                [2, 3],
                [1, 2, 6, 7],
                [3, 4, 8, 9],
                [5, 6, 10, 11],
                [7, 8, 12, 13],
                [11, 12, 16, 17],
                [13, 14, 18, 19],
                [15, 16, 20, 21],
                [17, 18, 22, 23],
                [21, 22],
                [23, 24],
            ],
        )
        self.assertEqual(
            c.x_stabilizers,
            [
                [0, 1, 5, 6],
                [2, 3, 7, 8],
                [4, 9],
                [5, 10],
                [6, 7, 11, 12],
                [8, 9, 13, 14],
                [10, 11, 15, 16],
                [12, 13, 17, 18],
                [14, 19],
                [15, 20],
                [16, 17, 21, 22],
                [18, 19, 23, 24],
            ],
        )
        self.assertEqual(
            c.z_stabilizers,
            [
                [0, 1, 5, 6, 10, 11, 15, 16, 20, 21],
                [1, 2, 6, 7, 11, 12, 16, 17, 21, 22],
                [2, 3, 7, 8, 12, 13, 17, 18, 22, 23],
                [3, 4, 8, 9, 13, 14, 18, 19, 23, 24],
            ],
        )


if __name__ == "__main__":
    unittest.main()
