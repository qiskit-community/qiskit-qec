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
from matplotlib.contour import QuadContourSet

import numpy as np

from qiskit import QiskitError

from qiskit_qec.linear.symplectic import *

class TestSymplectic(TestCase):

    def test_all_commute(self):
        test_mat = np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,0,0,1],[1,1,0,0,0,0]])
        self.assertTrue(all_commute(test_mat))
        test_mat = np.array([[1,0,0,0,0,0],[1,1,0,0,0,0],[0,0,0,0,0,1],[1,1,0,0,0,0]])
        self.assertFalse(all_commute(test_mat))

    def test_symplectic_product(self):
        test_mata=[[1,0,1,1],[0,1,1,0]]
        test_matb=[[0,1,1,1],[0,0,1,0]]

        # input matrices as lists
        mat1 = test_mata
        mat2 = test_matb
        result = symplectic_product(mat1, mat2)
        answer = np.array([[0, 1], [1, 0]], dtype=np.int8)
        self.assertEqual(result, answer)

        # input matrices as numpy arrays

        mat1 = np.array(test_mata)
        mat2 = np.array(test_matb)
        result = symplectic_product(mat1, mat2)
        answer = np.array([[0, 1], [1, 0]], dtype=np.int8)
        self.assertEqual(result, answer)

        # input vectors as lists
        
        mat1 = test_mata[0]
        mat2 = test_matb[0]
        result = symplectic_product(mat1, mat2)
        answer = 0
        self.assertEqual(result, answer)

        # input vectors as numpy array

        mat1 = np.array(test_mata)[0]
        mat2 = np.array(test_matb)[0]
        result = symplectic_product(mat1, mat2)
        answer = 0
        self.assertEqual(result, answer)

    def test_invalid_symplectic_product(self):
        test_mata=[[2,0,1,1],[0,1,1,0]]
        test_matb=[[0,1,1,1],[0,0,1,0]]
        self.assertRaises(QiskitError, symplectic_product(test_mata, test_matb))

        test_mata=[[1,0,1,1],[0,1,1,0]]
        test_matb=[[0,2,1,1],[0,0,1,0]]
        self.assertRaises(QiskitError, symplectic_product(test_mata, test_matb))

        test_mata=[[2,0,1,1],[0,1,1,0],[0,0,0,0]]
        test_matb=[[0,1,1,1],[0,0,1,0]]
        self.assertRaises(QiskitError, symplectic_product(test_mata, test_matb))

        # dim = 1 test cases

        test_mata=True
        test_matb=False
        self.assertRaises(QiskitError, symplectic_product(test_mata, test_matb))

        test_mata=[2,0,1,1]
        test_matb=[0,1,1,1,1,1]
        self.assertRaises(QiskitError, symplectic_product(test_mata, test_matb))

        test_mata=[2,0,1,1,1]
        test_matb=[0,1,1,1,1]
        self.assertRaises(QiskitError, symplectic_product(test_mata, test_matb))

        # dim = 2 test cases

        test_mata=[[1,0,1,1,1],[0,1,1,0,1]]
        test_matb=[[0,1,1,1,1],[0,0,1,0,1]]
        self.assertRaises(QiskitError, symplectic_product(test_mata, test_matb))

        # dim > 2

        test_mata=[[[1,0,1,1],[0,1,1,0]]]
        test_matb=[[[0,1,1,1],[0,0,1,0]]]
        self.assertRaises(QiskitError, symplectic_product(test_mata, test_matb))

