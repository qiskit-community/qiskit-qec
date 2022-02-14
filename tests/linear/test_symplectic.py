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
    make_element_commute_with_hyper_pair,
    make_element_commute_with_hyper_pairs,
    make_elements_commute_with_hyper_pair,
    symplectic_gram_schmidt,
    locate_hyper_partner,
)


class TestSymplectic(TestCase):
    """Tests simplectic."""

    def test_all_commute(self):
        """Tests all commute."""
        test_mat = np.array(
            [[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1], [1, 1, 0, 0, 0, 0]]
        )
        self.assertTrue(all_commute(test_mat))
        test_mat = np.array(
            [[1, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1], [1, 1, 0, 0, 0, 0]]
        )
        self.assertFalse(all_commute(test_mat))

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

    def test_make_element_commute_with_hyper_pair(self):
        """Tests make element commute."""
        test_matrix = np.array(
            [[1, 0, 0, 0, 0, 0], [1, 1, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1]],
            dtype=np.bool_,
        )
        target_answer = np.array([1, 1, 0, 0, 0, 0], dtype=np.bool_)
        answer = make_element_commute_with_hyper_pair(
            test_matrix[1], test_matrix[2], test_matrix[3]
        )
        self.assertTrue(np.array_equal(target_answer, answer))

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
                            make_element_commute_with_hyper_pair,
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
                            make_element_commute_with_hyper_pair,
                            vector,
                            hyper1,
                            hyper2,
                        )

    def test_make_element_commute_with_hyper_pairs(self):
        """Tests make element with pair commute."""
        elem = np.array([1, 1, 1, 0, 0, 0, 0, 0], dtype=np.bool_)  # X1X2X3
        hyper1 = np.array(
            [[0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0]], dtype=np.bool_
        )  # X2, X3
        hyper2 = np.array(
            [[0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0]], dtype=np.bool_
        )  # Z2, Z3
        range1 = [0, 1]
        range2 = [0, 1]
        target_answer = np.array([1, 0, 0, 0, 0, 0, 0, 0])
        answer = make_element_commute_with_hyper_pairs(elem, hyper1, hyper2, range1, range2)
        self.assertTrue(np.array_equal(answer, target_answer))

    def test_invalid_make_element_commute_with_hyper_pairs(self):
        """Tests make element commute with pairs invalid."""
        # assert vector.ndim == 1 and hyper1.ndim > 1 and hyper2.ndim > 1
        test_matrix_bad = np.array([[0, 1, 0, 0], [1, 1, 1, 1]])
        test_matrix_good = np.array([0, 1, 0, 0])
        test_hyper1_bad = np.array([0, 1, 0, 0])
        test_hyper1_good = np.array([[0, 1, 0, 0], [0, 0, 0, 1]])
        test_hyper2_bad = np.array([0, 1, 0, 0])
        test_hyper2_good = np.array([[0, 1, 0, 0], [0, 0, 0, 1]])
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
                            make_element_commute_with_hyper_pair,
                            vector,
                            hyper1,
                            hyper2,
                        )

        # assert vector.shape[0] == hyper1.shape[1] == hyper2.shape[1]
        vector = np.array([0, 1, 0, 0])
        test_hyper1_bad = np.array([[0, 1, 0, 0, 1], [0, 0, 0, 1, 1]])
        test_hyper1_good = np.array([[0, 1, 0, 0], [0, 0, 0, 1]])
        test_hyper2_bad = np.array([[0, 1, 0, 0, 1], [0, 0, 0, 1, 1]])
        test_hyper2_good = np.array([[0, 1, 0, 0], [0, 0, 0, 1]])
        for hyper1 in [test_hyper1_bad, test_hyper1_good]:
            for hyper2 in [test_hyper2_bad, test_hyper2_good]:
                if not (
                    np.array_equal(hyper1, test_hyper1_good)
                    and np.array_equal(hyper2, test_hyper2_good)
                ):
                    self.assertRaises(
                        QiskitError, make_element_commute_with_hyper_pair, vector, hyper1, hyper2
                    )

        # assert not (vector.shape[0]%2 or hyper1.shape[1]%2 or hyper2.shape[1]%2)
        test_matrix_bad = np.array([0, 1, 0, 0, 1])
        test_matrix_good = np.array([0, 1, 0, 0])
        test_hyper1_bad = np.array([[0, 1, 0, 0, 1], [0, 0, 0, 1, 1]])
        test_hyper1_good = np.array([[0, 1, 0, 0], [0, 0, 0, 1]])
        test_hyper2_bad = np.array([[0, 1, 0, 0, 1], [0, 0, 0, 1, 1]])
        test_hyper2_good = np.array([[0, 1, 0, 0], [0, 0, 0, 1]])

        self.assertRaises(
            QiskitError,
            make_element_commute_with_hyper_pair,
            test_matrix_bad,
            test_hyper1_good,
            test_hyper2_good,
        )

        self.assertRaises(
            QiskitError,
            make_element_commute_with_hyper_pair,
            test_matrix_good,
            test_hyper1_bad,
            test_hyper2_good,
        )

        self.assertRaises(
            QiskitError,
            make_element_commute_with_hyper_pair,
            test_matrix_good,
            test_hyper1_good,
            test_hyper2_bad,
        )

        range1 = False  # Non iterable object
        range2 = [0, 1]
        self.assertRaises(
            QiskitError,
            make_element_commute_with_hyper_pair,
            test_matrix_good,
            test_hyper1_good,
            test_hyper2_good,
        )

        self.assertRaises(
            QiskitError,
            make_element_commute_with_hyper_pair,
            test_matrix_good,
            test_hyper1_good,
            test_hyper2_good,
        )

        range1_good = [0, 1]
        range1_bad = [8, 10]
        range2_good = [0, 1]
        range2_bad = [17, 29]

        for range1 in [range1_good, range1_bad]:
            for range2 in [range2_good, range2_bad]:
                if not (range1 == range1_good and range2 == range2_good):
                    self.assertRaises(
                        QiskitError,
                        make_element_commute_with_hyper_pair,
                        test_matrix_good,
                        test_hyper1_good,
                        test_hyper2_good,
                    )

    def test_make_elements_commute_with_hyperbolic_pair(self):
        """Tests elements commute with hyperbolic pair."""
        matrix = np.array(
            [[1, 1, 1, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0]], dtype=np.bool_
        )  # X1X2X3, X2X3
        hyper1 = np.array([0, 1, 0, 0, 0, 0, 0, 0], dtype=np.bool_)  # X2
        hyper2 = np.array([0, 0, 0, 0, 0, 1, 0, 0], dtype=np.bool_)  # Z2
        mrange = [0, 1]
        target_result = np.array(
            [[1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0]], dtype=np.bool_
        )
        result = make_elements_commute_with_hyper_pair(matrix, mrange, hyper1, hyper2)
        self.assertTrue(np.array_equal(result, target_result))

    def test_find_noncommutative_partner(self):
        """Tests find noncommutive partner."""
        matrix = np.array(
            [[1, 0, 1, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0]], dtype=np.bool_
        )  # X1X2X3, X2X3
        vector = np.array([0, 0, 0, 0, 0, 1, 0, 0], dtype=np.bool_)  # Z2
        target_answer = (np.array([0, 1, 1, 0, 0, 0, 0, 0], dtype=np.bool_), 1)
        result = locate_hyper_partner(matrix, vector)
        self.assertTrue(np.array_equal(result[0], target_answer[0]))
        self.assertEqual(result[1], target_answer[1])

    def test_symplectic_gram_schmidt(self):
        """Tests gram schmidt."""
        matrix = np.array(
            [
                [0, 1, 0, 0, 1, 0, 1, 0],
                [0, 0, 0, 0, 1, 1, 0, 1],
                [1, 1, 1, 0, 0, 1, 0, 0],
                [1, 1, 0, 1, 0, 0, 0, 0],
            ],
            dtype=np.bool_,
        )
        target_center = np.array(
            [[1, 1, 1, 0, 1, 0, 0, 1], [1, 0, 0, 1, 0, 1, 1, 1]], dtype=np.bool_
        )
        target_hyper1 = np.array([[0, 1, 0, 0, 1, 0, 1, 0]], dtype=np.bool_)
        target_hyper2 = np.array([[0, 0, 0, 0, 1, 1, 0, 1]], dtype=np.bool_)
        center, hyper1, hyper2 = symplectic_gram_schmidt(matrix)
        self.assertTrue(
            np.array_equal(center, target_center)
            and np.array_equal(hyper1, target_hyper1)
            and np.array_equal(hyper2, target_hyper2)
        )
