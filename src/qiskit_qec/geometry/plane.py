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
# Part of the QEC framework
"""Module for Plane"""
from math import cos, sin

import numpy as np
from qiskit_qec.geometry.two_manifold import TwoManifold


class Plane(TwoManifold):
    """Represents a R2 plane two manifold"""

    # pylint: disable=useless-super-delegation
    def __init__(self):
        super().__init__()

    @staticmethod
    def ison(point):
        """Check that point is on the plane"""
        # Qick check that point is on the
        return point.shape == (2,) and len(point) == 2

    @staticmethod
    def rotate(theta, vector):
        """Rotate `vector` around `theta`"""
        theta = np.deg2rad(theta)
        rot = Plane.rot_matrix(theta)
        return np.dot(rot, vector)

    @staticmethod
    def rot_matrix(theta):
        """Create rotation matrix that rotates by theta"""
        return np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
