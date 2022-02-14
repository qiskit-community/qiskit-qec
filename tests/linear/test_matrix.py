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

from typing import List
from unittest import TestCase

import numpy as np

from qiskit import QiskitError

from qiskit_qec.linear.matrix import *


class TestLinearMatrix(TestCase):
    def test__create_lambda_matrix(self):
        lmatrix = _create_lambda_matrix(4)
        correct_lmatrix = np.array(
            [
                [False, False, False, True, False, False],
                [False, False, False, False, True, False],
                [False, False, False, False, False, True],
                [True, False, False, False, False, False],
                [False, True, False, False, False, False],
                [False, False, True, False, False, False],
            ]
        )
        self.assertTrue(isinstance(lmatrix, np.ndarray))
        self.assertEqual(lmatrix, correct_lmatrix)

    def test_invalid_create_lambda_matrix(self):
        self.assertRaises(QiskitError, create_lambda_matrix(-2))
        self.assertRaises(QiskitError, create_lambda_matrix(1.1))

    def test__augment_mat(self):
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
                [True, False, False, False, False, False, False, False, False, True, False, False],
                [False, True, False, False, False, False, False, False, False, False, True, False],
                [False, False, True, False, False, False, False, False, False, False, False, True],
                [False, False, False, True, False, False, True, False, False, False, False, False],
                [False, False, False, False, True, False, False, True, False, False, False, False],
                [False, False, False, False, False, True, False, False, True, False, False, False],
            ]
        )
        left_mat = _augment_mat(matrix, "left")
        self.assertEqual(left_mat, correct_left_mat)

        correct_right_mat = np.array(
            [
                [False, False, False, True, False, False, True, False, False, False, False, False],
                [False, False, False, False, True, False, False, True, False, False, False, False],
                [False, False, False, False, False, True, False, False, True, False, False, False],
                [True, False, False, False, False, False, False, False, False, True, False, False],
                [False, True, False, False, False, False, False, False, False, False, True, False],
                [False, False, True, False, False, False, False, False, False, False, False, True],
            ]
        )

        right_mat = _augment_mat(matrix, "right")
        self.assertEqual(right_mat, correct_right_mat)

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
        self.assertEqual(top_mat, correct_top_mat)

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
        self.assertEqual(bottom_mat, correct_bottom_mat)

    def test_invalid_augment_mat(self):
        matrix = np.array([[[1, 1], [1, 2]]])
        self.assertRaises(QiskitError, augment_mat(matrix, pos="left"))
        matrix = np.array([[1, 2, 3], [4, 45, 6]])
        self.assertRaises(QiskitError, augment_mat(matrix, pos="middle"))

    def test__rref_complete(self):
        matrix = np.array(
            [[1, 0, 0, 1, 1, 0, 1, 1], [1, 1, 1, 0, 0, 1, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0]],
            dtype=bool,
        )
        correct_heads = [1, 1, 0, 1, 0, 0, 0, 0]
        correct_rref_mat = np.array(
            [
                [True, False, False, False, False, False, False, False],
                [False, True, True, False, False, True, True, False],
                [False, False, False, True, True, False, True, True],
            ]
        )
        correct_transform_mat = np.array(
            [[False, False, True], [False, True, True], [True, False, True]]
        )
        correct_rank = 3
        heads, rref_mat, transform_mat, rank = _rref_complete(matrix)
        self.assertEqual(heads, correct_heads)
        self.assertEqual(rref_mat, correct_rref_mat)
        self.assertEqual(transform_mat, correct_transform_mat)
        self.assertEqual(rank, correct_rank)

    def test_invalid_rref_complete(self):
        not_array_obj = np.array
        self.assertRaises(QiskitError, rref_complete(not_array_obj))
        invalid_array = np.array([[[1, 1], [2, 1]]])
        self.assertRaises(QiskitError, rref_complete(invalid_array))
