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

# pylint: disable=invalid-name

"""Tests for Pauli operator class."""

import unittest
import itertools as it
from functools import lru_cache

import numpy as np
from ddt import ddt, data, unpack

from qiskit import QuantumCircuit
from qiskit.exceptions import QiskitError
from qiskit.test import QiskitTestCase

from qiskit.quantum_info.operators import Operator

from qiskit_qec.operators.random import random_clifford, random_pauli, random_xppauli
from qiskit_qec.operators.base_xp_pauli import BaseXPPauli
from qiskit_qec.operators.xp_pauli import XPPauli

# TODO from qiskit_qec.utils.pauli_rep import split_pauli, cpxstr2exp

# from qiskit.quantum_info.operators.symplectic.pauli import _split_pauli_label, _phase_from_label


@ddt 
class TestXPPauliInit(QiskitTestCase):
    """Tests for XPPauli initialization."""

    def test_array_init(self):
        """Test array initialization."""
        # Test case taken from Mark's paper on XPF
        matrix = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0], dtype=np.int64)
        phase_exp = 12
        precision = 8
        xppauli = XPPauli(data=matrix, phase_exp=phase_exp, precision=precision)
        np.testing.assert_equal(xppauli.matrix, np.atleast_2d(matrix))
        np.testing.assert_equal(xppauli._phase_exp, phase_exp)
        np.testing.assert_equal(xppauli.precision, precision)


@ddt
class TestXPPauli(QiskitTestCase):
    """Tests for XPPauli operator class."""

    def test_precision_rescale(self):
        """Test precision rescaling method."""
        # Test case taken from Mark's paper on XPF
        matrix = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0], dtype=np.int64)
        phase_exp = 12
        precision = 8
        new_precision = 2
        xppauli = XPPauli(data=matrix, phase_exp=phase_exp, precision=precision)
        rescaled_xppauli = xppauli.rescale_precision(new_precision=new_precision)
        target_matrix = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype=np.int64)
        target_phase_exp = 3
        target_xppauli = XPPauli(data=target_matrix, phase_exp=target_phase_exp, precision=new_precision)
        np.testing.assert_equal(target_xppauli.matrix, rescaled_xppauli.matrix)
        np.testing.assert_equal(target_xppauli._phase_exp, rescaled_xppauli._phase_exp)
        np.testing.assert_equal(target_xppauli.precision, rescaled_xppauli.precision)

    def test_weight(self):
        """Test weight method."""
        matrix = np.array([1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 1], dtype=np.int64)
        phase_exp = 12
        precision = 8
        xppauli = XPPauli(data=matrix, phase_exp=phase_exp, precision=precision)
        value = xppauli.weight()
        target = 5
        self.assertEqual(target, value)

    def test_diagonal(self):
        """Test is_diagonal method."""
        # Test case taken from Mark's paper, Table 5.
        matrix = np.array([0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3], dtype=np.int64)
        phase_exp = 0
        precision = 8
        xppauli = XPPauli(data=matrix, phase_exp=phase_exp, precision=precision)
        value = xppauli.is_diagonal()
        target = np.array([True])
        self.assertEqual(target, value)

        # Test case taken from Mark's paper, Table 5.
        matrix = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 3, 3, 3, 3], dtype=np.int64)
        phase_exp = 0
        precision = 8
        xppauli = XPPauli(data=matrix, phase_exp=phase_exp, precision=precision)
        value = xppauli.is_diagonal()
        target = np.array([True])
        self.assertEqual(target, value)

        matrix = np.array([0, 1, 0, 1, 0, 0, 0, 0, 1, 3, 3, 3, 3, 3], dtype=np.int64)
        phase_exp = 12
        precision = 8
        xppauli = XPPauli(data=matrix, phase_exp=phase_exp, precision=precision)
        value = xppauli.is_diagonal()
        target = np.array([False])
        self.assertEqual(target, value)

    def test_antisymmetric_op(self):
        """Test antisymmetric_op method."""

        matrix = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 3, 3, 3], dtype=np.int64)
        phase_exp = 0
        precision = 8
        xppauli = XPPauli(data=matrix, phase_exp=phase_exp, precision=precision)
        value = xppauli.antisymmetric_op()

        target_matrix = np.array([0, 0, 0, 0, 0, 0, 0, 0, -1, -2, -3, -3, -3, -3], dtype=np.int64)
        target_phase_exp = 15
        target_precision = 8
        target = XPPauli(data=target_matrix, phase_exp=target_phase_exp, precision=target_precision)
        self.assertEqual(target, value)

    def test_power(self):
        """Test power method."""
        matrix = np.array([1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 1], dtype=np.int64)
        phase_exp = 12
        precision = 8
        n = 5
        xppauli = XPPauli(data=matrix, phase_exp=phase_exp, precision=precision)
        value = xppauli.power(n=n)

        target_matrix = np.array([1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 5], dtype=np.int64)
        target_phase_exp = 15
        target_precision = 8
        target = XPPauli(data=target_matrix, phase_exp=target_phase_exp, precision=target_precision)
        self.assertEqual(target, value)

    def test_multiplication(self):
        """Test multiplication method."""
        # Test case taken from Mark's code.
        a_matrix = np.array([0,1,0,0,2,0], dtype=np.int64)
        a_phase_exp = 6
        a_precision = 4
        a = XPPauli(data=a_matrix, phase_exp=a_phase_exp, precision=a_precision)
        b_matrix = np.array([1,1,1,3,3,0], dtype=np.int64)
        b_phase_exp = 2
        b_precision = 4
        b = XPPauli(data=b_matrix, phase_exp=b_phase_exp, precision=b_precision)
        value = XPPauli.compose(a, b)

        target_matrix = np.array([1,0,1,3,3,0], dtype=np.int64)
        target_phase_exp = 6
        target_precision = 4
        target = XPPauli(data=target_matrix, phase_exp=target_phase_exp, precision=target_precision)

        self.assertEqual(target, value)

if __name__ == "__main__":
    unittest.main()
