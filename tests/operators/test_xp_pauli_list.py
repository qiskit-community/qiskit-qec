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


import itertools as it
import numpy as np
from ddt import ddt

from qiskit import QiskitError
from qiskit.quantum_info.operators import (
    Clifford,
    Operator,
    PauliTable,
    StabilizerTable,
)
from qiskit.quantum_info.random import random_clifford, random_pauli_list
from qiskit.test import QiskitTestCase

# Note: In Qiskit tests below is just test

from qiskit_qec.operators.base_xp_pauli import BaseXPPauli
from qiskit_qec.operators.xp_pauli import XPPauli
from qiskit_qec.operators.xp_pauli_list import XPPauliList

from tests import combine


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
    pass


if __name__ == "__main__":
    unittest.main()
