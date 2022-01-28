# surface code to colour code
# "I'm a leave my autograph" -- Queen LISA

# This code is part of Qiskit.
#
# (C) Copyright IBM 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""Tests for geometry folder."""

from unittest import TestCase

import numpy as np

from qiskit_qec.geometry.model.qubit_count import QubitCount
from qiskit_qec.geometry.model.qubit_data import QubitData
from qiskit_qec.geometry.tiles.tiling import Tiling
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.geometry.shape import Shape
from qiskit_qec.geometry.tiles.RSSC import RSSC
from qiskit_qec.geometry.lattice import Lattice


class TestGeometry(TestCase):
    """Test for geometry/"""

    def test_basic_regional_intersection(self):
        """Test basic intersection"""
        # rssc = RSSCode(Params)
        #