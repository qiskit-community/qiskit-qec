# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""Tests for XPPauliList class."""

import unittest
import numpy as np
from ddt import ddt
from qiskit.test import QiskitTestCase
from qiskit_qec.operators.xp_pauli_list import XPPauliList


class TestXPPauliListInit(QiskitTestCase):
    """Tests for XPPauliList initialization."""

    def test_array_init(self):
        """Test array initialization."""
        # Matrix array initialization
        precision = 8
        matrix1 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0], dtype=np.int64)
        phase_exp1 = 12
        matrix2 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0], dtype=np.int64)
        phase_exp2 = 2
        matrix = np.array([matrix1, matrix2])
        phase_exp = np.array([phase_exp1, phase_exp2])
        xppaulilist = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
        np.testing.assert_equal(xppaulilist.matrix, matrix)
        np.testing.assert_equal(xppaulilist._phase_exp, phase_exp)
        np.testing.assert_equal(xppaulilist.precision, precision)


@ddt
class TestXPPauliListOperator(QiskitTestCase):
    """Tests for XPPauliList base operator methods."""

    def test_precision_rescale(self):
        """Test precision rescaling method."""
        matrix1 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0], dtype=np.int64)
        phase_exp1 = 12
        matrix2 = np.array([1, 1, 0, 1, 0, 0, 0, 0, 0, 4, 0, 0, 0, 6], dtype=np.int64)
        phase_exp2 = 8
        precision = 8
        new_precision = 4
        matrix = np.array([matrix1, matrix2])
        phase_exp = np.array([phase_exp1, phase_exp2])

        xppaulilist = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
        rescaled_xppaulilist = xppaulilist.rescale_precision(new_precision=new_precision)
        target_matrix1 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0], dtype=np.int64)
        target_phase_exp1 = 6
        target_matrix2 = np.array([1, 1, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3], dtype=np.int64)
        target_phase_exp2 = 4
        target_matrix = np.array([target_matrix1, target_matrix2])
        target_phase_exp = np.array([target_phase_exp1, target_phase_exp2])
        target_xppaulilist = XPPauliList(
            data=target_matrix, phase_exp=target_phase_exp, precision=new_precision
        )
        np.testing.assert_equal(target_xppaulilist.matrix, rescaled_xppaulilist.matrix)
        np.testing.assert_equal(target_xppaulilist._phase_exp, rescaled_xppaulilist._phase_exp)
        np.testing.assert_equal(target_xppaulilist.precision, rescaled_xppaulilist.precision)

    def test_weight(self):
        """Test weight method."""
        precision = 8
        matrix1 = np.array([1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 1], dtype=np.int64)
        phase_exp1 = 12
        matrix2 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 4], dtype=np.int64)
        phase_exp2 = 3
        matrix = np.array([matrix1, matrix2])
        phase_exp = np.array([phase_exp1, phase_exp2])
        xppaulilist = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
        value = xppaulilist.weight()
        target = np.array([5, 4])
        np.testing.assert_equal(target, value)

    def test_diagonal(self):
        """Test is_diagonal method."""
        precision = 8
        matrix1 = np.array([0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3], dtype=np.int64)
        phase_exp1 = 0

        matrix2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 3, 3, 3, 3], dtype=np.int64)
        phase_exp2 = 0

        matrix3 = np.array([0, 1, 0, 1, 0, 0, 0, 0, 1, 3, 3, 3, 3, 3], dtype=np.int64)
        phase_exp3 = 12
        matrix = np.array([matrix1, matrix2, matrix3])
        phase_exp = np.array([phase_exp1, phase_exp2, phase_exp3])
        xppaulilist = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)

        value = xppaulilist.is_diagonal()
        target = np.array([True, True, False])
        np.testing.assert_equal(target, value)

    def test_antisymmetric_op(self):
        """Test antisymmetric_op method."""

        matrix = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 3, 3, 3],
                [0, 0, 0, 0, 0, 0, 0, 3, 1, 2, 3, 7, 6, 3],
            ],
            dtype=np.int64,
        )
        phase_exp = np.array([0, 0])
        precision = 8
        xppauli_list = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
        dinput = xppauli_list.z
        value = xppauli_list.antisymmetric_op(dinput)

        target_matrix = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0, 0, -1, -2, -3, -3, -3, -3],
                [0, 0, 0, 0, 0, 0, 0, -3, -1, -2, -3, -7, -6, -3],
            ],
            dtype=np.int64,
        )
        target_phase_exp = np.array([15, 25])
        target_precision = 8
        target = XPPauliList(
            data=target_matrix, phase_exp=target_phase_exp, precision=target_precision
        )
        np.testing.assert_equal(target.matrix, value.matrix)
        np.testing.assert_equal(target._phase_exp, value._phase_exp)
        np.testing.assert_equal(target.precision, value.precision)

    def test_multiplication(self):
        """Test multiplication method."""
        a_matrix = np.array([[0, 1, 0, 0, 2, 0], [0, 1, 0, 0, 2, 0]], dtype=np.int64)
        a_phase_exp = np.array([6, 6])
        a_precision = 4
        a = XPPauliList(data=a_matrix, phase_exp=a_phase_exp, precision=a_precision)
        b_matrix = np.array([[1, 1, 1, 3, 3, 0], [1, 1, 1, 3, 3, 0]], dtype=np.int64)
        b_phase_exp = np.array([2, 2])
        b_precision = 4
        b = XPPauliList(data=b_matrix, phase_exp=b_phase_exp, precision=b_precision)
        value = a.compose(b)

        target_matrix = np.array([[1, 0, 1, 3, 3, 0], [1, 0, 1, 3, 3, 0]], dtype=np.int64)
        target_phase_exp = np.array([6, 6])
        target_precision = 4
        target = XPPauliList(
            data=target_matrix, phase_exp=target_phase_exp, precision=target_precision
        )
        np.testing.assert_equal(target.matrix, value.matrix)
        np.testing.assert_equal(target._phase_exp, value._phase_exp)
        np.testing.assert_equal(target.precision, value.precision)

    def test_degree(self):
        """Test degree method."""
        matrix = np.array([[0, 0, 0, 2, 1, 0], [0, 0, 0, 2, 1, 0]], dtype=np.int64)
        phase_exp = np.array([2, 2])
        precision = 4
        xppauli_list = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
        value = xppauli_list.degree()

        target = np.array([4, 4])
        np.testing.assert_equal(target, value)

    def test_power(self):
        """Test power method."""
        matrix = np.array(
            [
                [1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 1],
                [1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 1],
            ],
            dtype=np.int64,
        )
        phase_exp = np.array([12, 12])
        precision = 8
        n = 5
        xppauli_list = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
        value = xppauli_list.power(n=n)

        target_matrix = np.array(
            [
                [1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 5],
                [1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 5],
            ],
            dtype=np.int64,
        )
        target_phase_exp = np.array([8, 8])
        target_precision = 8
        target = XPPauliList(
            data=target_matrix, phase_exp=target_phase_exp, precision=target_precision
        )
        np.testing.assert_equal(target.matrix, value.matrix)
        np.testing.assert_equal(target._phase_exp, value._phase_exp)
        np.testing.assert_equal(target.precision, value.precision)


if __name__ == "__main__":
    unittest.main()
