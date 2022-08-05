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

from qiskit_qec.linear.symplectic import (
    all_commute,
    symplectic_product,
    locate_hyper_partner,
    make_commute_hyper,
    build_hyper_partner,
    symplectic_gram_schmidt,
)


class TestSymplectic(TestCase):
    """Tests simplectic."""

    def test_all_commute(self):
        """Tests all commute."""
        matrix = np.array(
            [
                [1, 0, 0, 1, 0, 0, 1, 0],
                [0, 1, 1, 1, 0, 0, 0, 1],
                [1, 1, 1, 0, 1, 0, 0, 0],
                [1, 0, 0, 1, 0, 1, 0, 1],
            ],
            dtype=np.bool_,
        )
        self.assertFalse(all_commute(matrix))

        matrix = np.array(
            [
                [1, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0, 0, 1],
            ],
            dtype=np.bool_,
        )
        self.assertTrue(all_commute(matrix))

    def test_symplectic_product(self):
        """Tests symplectic product."""
        test_mata = [[1, 0, 1, 1], [0, 1, 1, 0]]
        test_matb = [[0, 1, 1, 1], [0, 0, 1, 0]]

        # input matrices as lists
        mat1 = test_mata
        mat2 = test_matb
        result = symplectic_product(mat1, mat2)
        answer = np.array([[0, 1], [1, 0]], dtype=np.int8)
        self.assertTrue(np.equal(result, answer).all())

        # input matrices as numpy arrays

        mat1 = np.array(test_mata)
        mat2 = np.array(test_matb)
        result = symplectic_product(mat1, mat2)
        answer = np.array([[0, 1], [1, 0]], dtype=np.int8)
        self.assertTrue(np.equal(result, answer).all())

        # input vectors as lists

        mat1 = test_mata[0]
        mat2 = test_matb[0]
        result = symplectic_product(mat1, mat2)
        answer = 0
        self.assertTrue(np.equal(result, answer).all())

        # input vectors as numpy array

        mat1 = np.array(test_mata)[0]
        mat2 = np.array(test_matb)[0]
        result = symplectic_product(mat1, mat2)
        answer = 0
        self.assertTrue(np.equal(result, answer).all())

        # input a single Pauli

        mat = np.array([1, 1, 0, 0], dtype=np.bool_)
        result = symplectic_product(mat, mat)
        answer = 0
        self.assertTrue(answer == result)

    def test_invalid_symplectic_product(self):
        """Tests invalid symplectic product."""
        test_mata = [[2, 0, 1, 1], [0, 1, 1, 0]]
        test_matb = [[0, 1, 1, 1], [0, 0, 1, 0]]
        self.assertRaises(QiskitError, symplectic_product, test_mata, test_matb)

        test_mata = [[1, 0, 1, 1], [0, 1, 1, 0]]
        test_matb = [[0, 2, 1, 1], [0, 0, 1, 0]]
        self.assertRaises(QiskitError, symplectic_product, test_mata, test_matb)

        test_mata = [[2, 0, 1, 1], [0, 1, 1, 0], [0, 0, 0, 0]]
        test_matb = [[0, 1, 1, 1], [0, 0, 1, 0]]
        self.assertRaises(QiskitError, symplectic_product, test_mata, test_matb)

        # dim = 1 test cases

        test_mata = True
        test_matb = False
        self.assertRaises(QiskitError, symplectic_product, test_mata, test_matb)

        test_mata = [2, 0, 1, 1]
        test_matb = [0, 1, 1, 1, 1, 1]
        self.assertRaises(QiskitError, symplectic_product, test_mata, test_matb)

        test_mata = [2, 0, 1, 1, 1]
        test_matb = [0, 1, 1, 1, 1]
        self.assertRaises(QiskitError, symplectic_product, test_mata, test_matb)

        # dim = 2 test cases

        test_mata = [[1, 0, 1, 1, 1], [0, 1, 1, 0, 1]]
        test_matb = [[0, 1, 1, 1, 1], [0, 0, 1, 0, 1]]
        self.assertRaises(QiskitError, symplectic_product, test_mata, test_matb)

        # dim > 2

        test_mata = [[[1, 0, 1, 1], [0, 1, 1, 0]]]
        test_matb = [[[0, 1, 1, 1], [0, 0, 1, 0]]]
        self.assertRaises(QiskitError, symplectic_product, test_mata, test_matb)

    def test_make_commute_hyper(self):
        """Tests make_commute_hyper."""
        a = np.array([1, 1, 1, 0, 0, 0], dtype=np.bool_)
        x = np.array([0, 0, 1, 0, 0, 0], dtype=np.bool_)
        z = np.array([0, 0, 0, 0, 0, 1], dtype=np.bool_)
        a = make_commute_hyper(a, x, z).astype(int)
        a_expected = [1, 1, 0, 0, 0, 0]
        self.assertTrue(np.array_equal(a, a_expected))

        a = np.array([1, 1, 1, 0, 0, 0, 0, 0], dtype=np.bool_)
        x = np.array([[0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0]], dtype=np.bool_)
        z = np.array([[0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0]], dtype=np.bool_)
        xrange = [0, 1]
        zrange = [0, 1]
        a = make_commute_hyper(a, x, z, xrange=xrange, zrange=zrange).astype(int)
        a_expected = [1, 0, 0, 0, 0, 0, 0, 0]
        self.assertTrue(np.array_equal(a, a_expected))

        a = np.array(
            [[1, 1, 1, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0]], dtype=np.bool_
        )  # X1X2X3, X2X3
        x = np.array([0, 1, 0, 0, 0, 0, 0, 0], dtype=np.bool_)  # X2
        z = np.array([0, 0, 0, 0, 0, 1, 0, 0], dtype=np.bool_)  # Z2
        arange = [0, 1]
        a = make_commute_hyper(a, x, z, arange).astype(int)
        a_expected = [[1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0]]
        self.assertTrue(np.array_equal(a, a_expected))

        a = np.array([[1, 1, 1, 0, 0, 0, 0, 0], [0, 1, 1, 1, 0, 0, 0, 0]], dtype=np.bool_)
        x = np.array([[0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0]], dtype=np.bool_)
        z = np.array([[0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0]], dtype=np.bool_)
        arange = [0, 1]
        a = make_commute_hyper(a, x, z, arange).astype(int)
        a_expected = [[1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0]]
        self.assertTrue(np.array_equal(a, a_expected))

    def test_invalid_element_commute_with_hyper_pair(self):
        """Tests make element commute invalid."""
        # assert vector.ndim == 1 and hyper1.ndim == 1 and hyper2.ndim == 1
        test_matrix_bad = np.array([[0, 1, 0, 0, 1], [0, 0, 0, 0, 1]])
        test_matrix_good = np.array([0, 1, 0, 0])
        test_hyper1_bad = np.array([[0, 1, 0, 0, 1], [0, 0, 0, 0, 1]])
        test_hyper1_good = np.array([0, 1, 0, 0])
        test_hyper2_bad = np.array([[0, 1, 0, 0, 1], [0, 0, 0, 0, 1]])
        test_hyper2_good = np.array([0, 1, 0, 0])

        for vector in [test_matrix_bad, test_matrix_good]:
            for hyper1 in [test_hyper1_bad, test_hyper1_good]:
                for hyper2 in [test_hyper2_bad, test_hyper2_good]:
                    if (
                        not np.array_equal(vector, test_matrix_good)
                        and not np.array_equal(hyper1, test_hyper1_good)
                        and not np.array_equal(hyper2, test_hyper2_good)
                    ):
                        self.assertRaises(
                            QiskitError,
                            make_commute_hyper,
                            vector,
                            hyper1,
                            hyper2,
                        )

        # assert not (vector.shape[0]%2 or hyper1.shape[0]%2 or hyper2.shape[0]%2)
        test_matrix_bad = np.array([0, 1, 0, 0, 1])
        test_hyper1_bad = np.array([0, 1, 0, 0, 1])
        test_hyper2_bad = np.array([0, 1, 0, 0, 1])
        for vector in [test_matrix_bad, test_matrix_good]:
            for hyper1 in [test_hyper1_bad, test_hyper1_good]:
                for hyper2 in [test_hyper2_bad, test_hyper2_good]:
                    if (
                        not np.array_equal(vector, test_matrix_good)
                        and not np.array_equal(hyper1, test_hyper1_good)
                        and not np.array_equal(hyper2, test_hyper2_good)
                    ):
                        self.assertRaises(
                            QiskitError,
                            make_commute_hyper,
                            vector,
                            hyper1,
                            hyper2,
                        )

    def test_find_noncommutative_partner(self):
        """Tests find noncommutive partner."""
        matrix = np.array([[1, 0, 1, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0]], dtype=np.bool_)
        vector = np.array([0, 0, 0, 0, 0, 1, 0, 0], dtype=np.bool_)
        av, index = locate_hyper_partner(matrix, vector)
        av.astype(int)
        self.assertTrue(np.array_equal(av.astype(int), [0, 1, 1, 0, 0, 0, 0, 0]))
        self.assertEqual(index, 1)

    def test_build_hyper_partner(self):
        """Tests build_hyper_partner."""
        matrix = np.array(
            [
                [1, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0],
            ],
            dtype=np.bool_,
        )
        av = build_hyper_partner(matrix, 0).astype(int)
        self.assertTrue(np.array_equal([0, 0, 0, 0, 1, 0, 0, 0], av))

    def test_symplectic_gram_schmidt(self):
        """Tests gram schmidt."""
        a = np.array(
            [
                [0, 1, 0, 0, 1, 0, 1, 0],
                [0, 0, 0, 0, 1, 1, 0, 1],
                [1, 1, 1, 0, 0, 1, 0, 0],
                [1, 1, 0, 1, 0, 0, 0, 0],
            ],
            dtype=np.bool_,
        )
        center, x, z = symplectic_gram_schmidt(a)
        self.assertTrue(
            np.array_equal(center.astype(int), [[1, 1, 1, 0, 1, 0, 0, 1], [1, 0, 0, 1, 0, 1, 1, 1]])
        )
        self.assertTrue(np.array_equal(x.astype(int), [[0, 1, 0, 0, 1, 0, 1, 0]]))
        self.assertTrue(np.array_equal(z.astype(int), [[0, 0, 0, 0, 1, 1, 0, 1]]))
