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
        matrix1 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0], dtype=np.int64)
        phase_exp1 = 12
        precision1 = 8
        matrix2 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0], dtype=np.int64)
        phase_exp2 = 2
        precision2 = 4
        matrix = np.array([matrix1, matrix2])
        phase_exp = np.array([phase_exp1, phase_exp2])
        precision = np.array([precision1, precision2])
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
        precision1 = 8
        new_precision1 = 2
        matrix2 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 8], dtype=np.int64)
        phase_exp2 = 6
        precision2 = 10
        new_precision2 = 5
        matrix = np.array([matrix1, matrix2])
        phase_exp = np.array([phase_exp1, phase_exp2])
        precision = np.array([precision1, precision2])
        new_precision = np.array([new_precision1, new_precision2])

        xppaulilist = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
        rescaled_xppaulilist = xppaulilist.rescale_precision(new_precision=new_precision)
        target_matrix1 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype=np.int64)
        target_phase_exp1 = 3
        target_matrix2 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 4], dtype=np.int64)
        target_phase_exp2 = 3
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
        matrix1 = np.array([1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 1], dtype=np.int64)
        phase_exp1 = 12
        precision1 = 8
        matrix2 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 4], dtype=np.int64)
        phase_exp2 = 3
        precision2 = 5
        matrix = np.array([matrix1, matrix2])
        phase_exp = np.array([phase_exp1, phase_exp2])
        precision = np.array([precision1, precision2])
        xppaulilist = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
        value = xppaulilist.weight()
        target = np.array([5, 4])
        np.testing.assert_equal(target, value)

    def test_diagonal(self):
        """Test is_diagonal method."""
        matrix1 = np.array([0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3], dtype=np.int64)
        phase_exp1 = 0
        precision1 = 8

        matrix2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 3, 3, 3, 3], dtype=np.int64)
        phase_exp2 = 0
        precision2 = 5

        matrix3 = np.array([0, 1, 0, 1, 0, 0, 0, 0, 1, 3, 3, 3, 3, 3], dtype=np.int64)
        phase_exp3 = 12
        precision3 = 8
        matrix = np.array([matrix1, matrix2, matrix3])
        phase_exp = np.array([phase_exp1, phase_exp2, phase_exp3])
        precision = np.array([precision1, precision2, precision3])
        xppaulilist = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)

        value = xppaulilist.is_diagonal()
        target = np.array([True, True, False])
        np.testing.assert_equal(target, value)


if __name__ == "__main__":
    unittest.main()
