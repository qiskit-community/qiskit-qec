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
from typing import List
import copy

import numpy as np
from qiskit.exceptions import QiskitError


class GeometryBounds:
    """A simple Bounding Box class for an AABB in R^n"""

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
        retstr += f"center ={self.center}\n"
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

    def intercepts(self, line: List[float]) -> List:
        """Returns the intercepts of the input line with the aabb.

        The line is inputed as a vector [a,b,c] which represents the
        line given by ax+by=c

        Args:
            line (List[float, float, float]): _description_

        Returns:
            List: _description_
        """

        a, b, c = line

        xmin, ymin = self.min
        xmax, ymax = self.max

        intercepts = []

        if abs(b) < 0.0000001:
            if abs(c - xmin) < 0.0000001:
                intercepts.append([xmin, ymin])
                intercepts.append([xmin, ymax])
            elif xmin <= c <= xmax:
                intercepts.append([c, ymin])
                intercepts.append([c, ymax])
            elif abs(c - xmax) < 0.0000001:
                intercepts.append([xmax, ymin])
                intercepts.append([xmax, ymax])
        else:
            # |.
            y = (c - a * xmin) / b
            if ymin <= y <= ymax:
                intercepts.append([xmin, y])
            # .|
            y = (c - a * xmax) / b
            if ymin <= y <= ymax:
                intercepts.append([xmax, y])
            # ‾
            if abs(a) < 0.0000001 and abs(c / b - ymax) < 0.0000001:
                intercepts.append([xmin, ymax])
                intercepts.append([xmax, ymax])
            elif abs(a) > 0.0000001:
                x = (c - b * ymax) / a
                if xmin <= x <= xmax:
                    intercepts.append([x, ymax])
            # _
            if abs(a) < 0.0000001 and abs(c / b - ymin) < 0.0000001:
                intercepts.append([xmin, ymin])
                intercepts.append([xmax, ymin])
            elif abs(a) > 0.0000001:
                x = (c - b * ymin) / a
                if xmin <= x <= xmax:
                    intercepts.append([x, ymin])
        result = []
        if len(intercepts) > 0:
            for item in intercepts:
                flag = True
                for value in result:
                    if abs(item[0] - value[0]) < 0.0000001 and abs(item[1] - value[1]) < 0.0000001:
                        flag = False
                        continue
                if flag:
                    result.append(item)
        return result
