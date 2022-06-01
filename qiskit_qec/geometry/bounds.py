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
"""Module for Geometry Bounds"""
import copy

import numpy as np
from qiskit.exceptions import QiskitError


class GeometryBounds:
    """A simple Bounding Box class for an AABB in Rn"""

    def __init__(self, center=None, size=None, dim: int = None) -> None:
        """Initialize Geometry Bounds

        Args:
            center (nd.ndarry, optional): Center of bounding box. Defaults to None.
            size (nd.ndarry, optional): Size of boundeing box. Defaults to None.
            dim (int, optional): Dimensions of geometric space. Defaults to None.

        Raises:
            QiskitError: [description]
        """

        if center is None:
            if size is None:
                if dim is None:
                    center = np.array([0, 0, 0])
                    size = np.array([0, 0, 0])
                else:
                    center = np.zeros(dim)
                    size = np.zeros(dim)
            else:
                raise QiskitError("Both center and size must be spefified togther or not at all")

        self.center = np.array(center)
        self.size = np.array(size)

        assert self.center.shape == (self.center.size,), "Center not correct shape"
        assert self.size.shape == (self.size.size,), "Size not correct shape"

        self.min = self.center - self.size / 2
        self.max = self.center + self.size / 2

    def __str__(self):
        retstr = "GeometryBounds\n"
        retstr += "center ={self.center}\n"
        retstr += f"size   ={self.size}\n"
        retstr += f"min    = {self.min}\n"
        retstr += f"max    = {self.max}"
        return retstr

    def copy(self):
        """Copying self"""
        bounds = copy.copy(self)  # TODO shouldn't we do bounds = GeometryBounds()?
        bounds.center = self.center.copy()
        bounds.size = self.size.copy()
        bounds.min = self.min.copy()
        bounds.max = self.max.copy()
        return bounds

    def contains(self, point) -> bool:
        """Returns true is point is within bounds

        Args:
            point (nd.ndarry): Point to test

        Returns:
            bool: True is point in within bounds. Else, false.
        """
        return np.all((point - self.min) >= 0) and np.all((self.max - point) >= 0)

    def set_min_max(self, min_point, max_point):
        """Set the bound max to be max_point and bound min to min_point

        Args:
            min_point (nd.ndarry):  min point
            max_point (nd.ndarry): max
        """
        min_point = np.array(min_point)
        max_point = np.array(max_point)

        assert min_point.shape == (min_point.size,), "Min point not correct shape"
        assert max_point.shape == (max_point.size,), "Max point not correct shape"
        assert min_point.size == max_point.size, "min and max vectors must be the same dimension"

        self.min = min_point
        self.max = max_point
        self.center = (self.min + self.max) / 2
        self.size = 2 * (self.center - self.min)

    @classmethod
    def combine(cls, bounds1, bounds2):
        """Create the smallest AABB bounding box that includes bounds1 and bounds2

        Args:
            bounds1 (GeometryBounds): AABB
            bounds2 (GeometryBounds): AABB

        Returns:
            GeometryBounds: Smallest AABB bounding box that includes bounds1 and bounds2
        """
        new_min = np.minimum(bounds1.min, bounds2.min)
        new_max = np.maximum(bounds1.max, bounds2.max)
        bounds = GeometryBounds()
        bounds.set_min_max(new_min, new_max)

        return bounds

    def expand(self, amount: np.array):
        """Expand the bounds of the AABB by increasing its size by amount along each axis

        Args:
            amount (nd.ndarry): Vector to increase bounds of AABB along each axis
        """
        assert (
            self.center.size == amount.size
        ), "Amount vector not the same dimension as bounding box"
        self.set_min_max(self.min - amount, self.max + amount)

    @staticmethod  # TODO should probably put this in init or make factories?
    def bounding_box_from_line(point1, point2):
        """Create bounding box from points1 and points2"""

        new_min = np.minimum(point1, point2)
        new_max = np.maximum(point1, point2)

        bounds = GeometryBounds()
        bounds.set_min_max(new_min, new_max)

        return bounds
