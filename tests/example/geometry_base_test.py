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
    """Test for geometry"""

    def test_basic_regional_intersection(self):
        """Test basic intersection"""
        # Set parameters for RSSC code on the plane
        d = 3
        manifold = Plane()

        # Create the bounding or cutter shape
        cutter = Shape.square(
            origin=(0, -1), direction=(1, 1), length=d - 1, manifold=manifold, dtype=int
        )

        # Set up the qubit counter and aux data structures
        qubit_count = QubitCount()
        qubit_data = QubitData()

        # Choose a tile
        tile = RSSC

        # Basis for Tiling or get from Tile directly
        u_vec = np.array([2, 0])
        v_vec = np.array([0, 2])

        # Create Tiling lattice
        lattice = Lattice(u_vec, v_vec)

        lattice = lattice.restrict_for_tiling(cutter, tile=tile, expand_value=np.array([0.5, 0.5]))

        # Tile the restriced lattice L_l with RSSC tiles
        tiling = Tiling(
            tile_type=tile,
            lattice=lattice,
            qubit_count=qubit_count,
            qubit_data=qubit_data,
        )

        # Determine the intersection of faces/operators with defined region/cutter
        intersection = cutter.intersection(tiling)

        # TODO: fix tests
        self.assertTrue(intersection)
