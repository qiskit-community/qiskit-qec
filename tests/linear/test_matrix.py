# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""Tests symplectic."""

from unittest import TestCase

import numpy as np
from qiskit import QiskitError

from qiskit_qec.linear.matrix import (
    _augment_mat,
    _rref_complete,
    create_lambda_matrix,
    augment_mat,
    rref_complete,
    rank,
)


class TestLinearMatrix(TestCase):
    """Tests matrix."""

    def test_create_lambda_matrix(self):
        """Tests creation of lambda matrix."""
        lmatrix = create_lambda_matrix(2)
        correct_lmatrix = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
        self.assertTrue(isinstance(lmatrix, np.ndarray))
        self.assertTrue(np.equal(lmatrix, correct_lmatrix).all())

    def test_invalid_create_lambda_matrix(self):
        """Invalid creation."""
        self.assertRaises(QiskitError, create_lambda_matrix, -2)
        self.assertRaises(QiskitError, create_lambda_matrix, 1.1)

    def test_augment_mat(self):
        """Tests augment matrix."""
        matrix = np.array(
            [
                [False, False, False, True, False, False],
                [False, False, False, False, True, False],
                [False, False, False, False, False, True],
                [True, False, False, False, False, False],
                [False, True, False, False, False, False],
                [False, False, True, False, False, False],
            ]
        )

        correct_left_mat = np.array(
            [
                [
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    True,
                    False,
                    False,
                ],
                [
                    False,
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    True,
                    False,
                ],
                [
                    False,
                    False,
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    True,
                ],
                [
                    False,
                    False,
                    False,
                    True,
                    False,
                    False,
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                ],
                [
                    False,
                    False,
                    False,
                    False,
                    True,
                    False,
                    False,
                    True,
                    False,
                    False,
                    False,
                    False,
                ],
                [
                    False,
                    False,
                    False,
                    False,
                    False,
                    True,
                    False,
                    False,
                    True,
                    False,
                    False,
                    False,
                ],
            ]
        )
        left_mat = _augment_mat(matrix, "left")
        self.assertTrue(np.equal(left_mat, correct_left_mat).all())

        correct_right_mat = np.array(
            [
                [
                    False,
                    False,
                    False,
                    True,
                    False,
                    False,
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                ],
                [
                    False,
                    False,
                    False,
                    False,
                    True,
                    False,
                    False,
                    True,
                    False,
                    False,
                    False,
                    False,
                ],
                [
                    False,
                    False,
                    False,
                    False,
                    False,
                    True,
                    False,
                    False,
                    True,
                    False,
                    False,
                    False,
                ],
                [
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    True,
                    False,
                    False,
                ],
                [
                    False,
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    True,
                    False,
                ],
                [
                    False,
                    False,
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    True,
                ],
            ]
        )

        right_mat = _augment_mat(matrix, "right")
        self.assertTrue(np.equal(right_mat, correct_right_mat).all())

        correct_top_mat = np.array(
            [
                [True, False, False, False, False, False],
                [False, True, False, False, False, False],
                [False, False, True, False, False, False],
                [False, False, False, True, False, False],
                [False, False, False, False, True, False],
                [False, False, False, False, False, True],
                [False, False, False, True, False, False],
                [False, False, False, False, True, False],
                [False, False, False, False, False, True],
                [True, False, False, False, False, False],
                [False, True, False, False, False, False],
                [False, False, True, False, False, False],
            ]
        )

        top_mat = _augment_mat(matrix, "top")
        self.assertTrue(np.equal(top_mat, correct_top_mat).all())

        correct_bottom_mat = np.array(
            [
                [False, False, False, True, False, False],
                [False, False, False, False, True, False],
                [False, False, False, False, False, True],
                [True, False, False, False, False, False],
                [False, True, False, False, False, False],
                [False, False, True, False, False, False],
                [True, False, False, False, False, False],
                [False, True, False, False, False, False],
                [False, False, True, False, False, False],
                [False, False, False, True, False, False],
                [False, False, False, False, True, False],
                [False, False, False, False, False, True],
            ]
        )

        bottom_mat = _augment_mat(matrix, "bottom")
        self.assertTrue(np.equal(bottom_mat, correct_bottom_mat).all())

    def test_invalid_augment_mat(self):
        """Tests invalid augment matrix."""
        matrix = np.array([[[1, 1], [1, 2]]])
        self.assertRaises(QiskitError, augment_mat, {"matrix": matrix, "pos": "left"})
        matrix = np.array([[1, 2, 3], [4, 45, 6]])
        self.assertRaises(QiskitError, augment_mat, {"matrix": matrix, "pos": "middle"})

    def test_rref_complete(self):
        """Tests rref complete."""
        matrix = np.array(
            [
                [1, 0, 0, 1, 0, 0, 1, 0],
                [0, 1, 1, 1, 0, 0, 0, 1],
                [1, 1, 1, 0, 1, 0, 0, 0],
                [1, 0, 0, 1, 0, 1, 0, 1],
            ],
            dtype=np.bool_,
        )
        heads, rref_mat, transform_mat, rank_ = _rref_complete(matrix)
        expected_heads = [1, 1, 0, 0, 1, 1, 0, 0]
        expected_rref_mat = np.array(
            [
                [1, 0, 0, 1, 0, 0, 1, 0],
                [0, 1, 1, 1, 0, 0, 0, 1],
                [0, 0, 0, 0, 1, 0, 1, 1],
                [0, 0, 0, 0, 0, 1, 1, 1],
            ]
        )
        expected_transform_mat = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 1, 0], [1, 0, 0, 1]])
        expected_rank = 4
        self.assertTrue(np.array_equal(heads, expected_heads))
        self.assertTrue(np.array_equal(expected_rref_mat, rref_mat.astype(int)))
        self.assertTrue(np.array_equal(expected_transform_mat, transform_mat.astype(int)))
        self.assertEqual(rank_, expected_rank)

    def test_invalid_rref_complete(self):
        """Tests invalid rref complete."""
        not_array_obj = np.array
        self.assertRaises(QiskitError, rref_complete, not_array_obj)
        invalid_array = np.array([[[1, 1], [2, 1]]])
        self.assertRaises(QiskitError, rref_complete, invalid_array)

    def test_rank(self):
        """Tests rank."""
        matrix = np.array(
            [
                [1, 0, 0, 1, 0, 0, 1, 0],
                [0, 1, 1, 1, 0, 0, 0, 1],
                [1, 1, 1, 0, 1, 0, 0, 0],
                [1, 0, 0, 1, 0, 1, 0, 1],
            ],
            dtype=np.bool_,
        )
        self.assertEqual(rank(matrix), 4)
