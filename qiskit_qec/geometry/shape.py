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
"""Module for shape"""
from typing import List

from math import copysign

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.bounds import GeometryBounds
from qiskit_qec.geometry.manifold import Manifold
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.geometry.model.shell import Shell
from qiskit_qec.geometry.face_partition import FacePartition


# pylint: disable=anomalous-backslash-in-string
class Shape:
    """This class is used to store a boundary shape on a given manifold that is used to
    select sublattices.

    """

    OUTSIDE = 0
    INSIDE = 1

    def __init__(self, points: List, lines: List[List], indices: List = None):
        """Init shape

        Args:
            points (List): points on the shape
            lines (List[List]): lines on the shape
            indices (List, optional): indices of the shape Defaults to None.
        """
        self.points = points
        self.lines = lines
        if indices is not None:
            self.indices = indices
        else:
            self.indices = list(range(len(self.points)))

        self.bounds = self.bounding_box_from_lines()

    def bounding_box_from_lines(self):
        """Generate bounding box from lines on the shape"""
        line = self.lines[0]
        start = self.points[line[0]]
        finish = self.points[line[1]]
        bounds = GeometryBounds.bounding_box_from_line(start, finish)
        for line in self.lines:
            bounds2 = GeometryBounds.bounding_box_from_line(
                self.points[line[0]], self.points[line[1]]
            )
            bounds = GeometryBounds.combine(bounds, bounds2)
        return bounds

    @classmethod
    def square(cls, origin: List, direction: List, length, manifold: Manifold, dtype=int):

        """Create a square (2d) shape on the 2-manifold.

        The square has one corner placed at the origin (origin must be on the plane). The the
        square is defined as the path following heading in the direction <direction> towards the
        next corner of length <length>. The proceeding to the next corner and so on.

        Args:
            origin ([coordinate]): a coordinate on the plane
            direction ([vector]): vector describing the direction of the second corner
            length ([real]): length of the side of the square
            manifold ([TwoManifold]): 2-manifold
            dtype ([dtype], optional): Date type of coordinate entries. Defaults to int.

        Returns:
            [Shape]: Shape obj describing a square of length <length> in the manifold <manifold>
        """

        origin = np.asarray(origin)
        direction = np.asarray(direction)

        return cls.rect(
            origin=origin,
            direction=direction,
            length1=length,
            length2=length,
            manifold=manifold,
            dtype=dtype,
        )

    @classmethod
    def rect(cls, origin, direction, length1, length2, manifold, dtype=int):

        """Create a rectangle on a manifold

                    r2
                    o
                /      \
         r3  o          \
              \          o r1
               \       /
                  o
              r0=origin

        direction (r0 to r1)
        length1 = length(r0,r1)
        length2 = length(r1,r2)

        Args:
            origin ([type]): [description]
            direction ([type]): [description]
            length1 ([type]): [description]
            length2 ([type]): [description]
            dtype ([type], optional): [description]. Defaults to int.

        Returns:
            (Shape): rectangle

        Raises:
            QiskitError: [description]
        """

        assert isinstance(origin, np.ndarray)
        assert isinstance(direction, np.ndarray)

        assert length1 > 0 and length2 > 0, "Rectangle lengths must be positive"
        assert manifold.ison(origin), f"{origin} must be on the surface if the manifold"

        if isinstance(manifold, Plane):
            r0 = origin.astype(dtype)
            r1 = origin + length1 * direction
            r1 = r1.astype(int)
            direction = Plane.rotate(theta=90, vector=direction)
            r2 = r1 + length2 * direction
            r2 = r2.astype(int)
            direction = Plane.rotate(theta=90, vector=direction)
            r3 = r2 + length1 * direction
            r3 = r3.astype(int)

            points = [r0, r1, r2, r3]
            lines = [[0, 1], [1, 2], [2, 3], [3, 0]]

            return cls(points, lines)

        else:
            raise QiskitError(f"Manifold {manifold} not yet supported")

    def contains(self, point, on_boundary=True, epsilon=0.001):
        """Check if inside the bounded region using an infinite horizontal line from the point
        to +infinity

        This implementation assumes that the path is given by a sequence of
        straight lines.

        A more sophisticated verison is required if the path is otherwise.

        Args:
            point (point): point to check if in region (including on boundary)
            on_boundary (bool): True if boundary to be include in "inside". Default is True
            epsilon (real): Tolerance distance between two points to be equal. Default is 0.001

        Returns: (bool) Check if inside the bounded region using an infinite
                        horizontal line from the point to +infinity
        """
        x = point[0]
        y = point[1]
        count = 0
        k = len(self.lines)

        for index, line in enumerate(self.lines):
            [x0, y0] = self.points[line[0]]
            [x1, y1] = self.points[line[1]]
            [_, y2] = self.points[self.lines[(index + 1) % k][1]]
            m = (y0 - y1) / (x0 - x1)
            sign = copysign(1, m)
            c = y0 - m * x0
            yp = m * x + c

            # Point on line
            if y == yp and Shape.is_between(x, x0, x1):
                return True and on_boundary

            # Does the ray pass through the middle of a line
            if Shape.is_between(y, y0, y1, strict=True):
                if x <= min(x0, x1) or sign * y < sign * yp:
                    count += 1
                else:
                    continue

            # Does ray pass though a vertex

            if abs(y - y1) < epsilon:
                if x < x1:
                    if y1 > max(y0, y2) or y1 < min(y0, y2):
                        # ray passes through vertex y1 /\ or \/
                        count += 2
                    else:
                        # ray passes through vertex y1 |
                        count += 1

        return bool(count % 2)

    @staticmethod
    def is_between(p, a, b, strict=False):
        """is p between a and b: a <= p < =b

        Args:
            p (real): point
            a (real): end point
            b (real): end point
            strict (bool): is strict comparison

        Returns:
            bool: True is a is between a and b, False if not
        """
        # pylint: disable=chained-comparison
        if strict:
            return (((a <= b) * a + (b < a) * b) < p) and (p < ((a > b) * a + (b >= a) * b))

        return (((a <= b) * a + (b < a) * b) <= p) and (p <= ((a > b) * a + (b >= a) * b))

    def intersection(self, tiling: Shell, levels: int = 4) -> FacePartition:
        """Find the vertex intersection of faces/operators with region defined by
        the cutter/shape (self) on a tiling (Shell)

        Args:
            tiling (Shell): A Shell representation a set of faces/operators
            levels (int, optional): Bound on the number of intersection per face to include.
                So if levels is 4 then all levels < 4 will be included. Defaults to 4

        Returns:
            FacePartition: Details of faces partitioned by the cutter on the tiling (Shell)
        """
        face_partition = FacePartition()

        for face in tiling.faces:
            face_partition.add(face)
            for vertex in face.vertices:
                if self.contains(vertex.pos):
                    face_partition.faces[face][Shape.INSIDE].append(vertex)
                else:
                    face_partition.faces[face][Shape.OUTSIDE].append(vertex)
            # Add procedure here to detect and store information about
            # edges that cross the cutter/region boundary. This is needed
            # to allow the tiling to be mapped to other non-plane surfaces
            # such as the torus etc.

        # Split faces/operators into intersection level categories
        parts = [[] for i in range(levels)]
        for key, vert_lists in face_partition.faces.items():
            in_vert_count = len(vert_lists[Shape.INSIDE])
            if in_vert_count < levels:
                parts[in_vert_count].append(key)

        face_partition.set_parts(parts)

        return face_partition
