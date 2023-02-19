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

from numbers import Real
from typing import List, Dict, Any, Tuple, Optional, Union
from math import sqrt, copysign

import numpy as np

from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.bounds import GeometryBounds
from qiskit_qec.geometry.manifold import Manifold
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.geometry.model.shell import Shell


class Shape:
    """This class is used to store a boundary shape on a given manifold that is used to
    select sublattices.
    """

    _DEBUG = False

    def __init__(self, points: List, lines: Optional[List[List]] = None, indices: List = None):
        """Init shape

        Args:
            points (List): points on the shape
            lines (List[List], optional): lines on the shape. Defaults to None
            indices (List, optional): indices of the shape Defaults to None.

        Raise:
            QiskitError: If <indices> are defined then lines must also be defined
        """
        self.points = points
        if lines is None:
            if lines is None and indices is not None:
                raise QiskitError("If <indices> are defined then lines must also be defined")
            lines = self.create_lines(points)
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

    @staticmethod
    def create_lines(points: List):
        """Creates Lines from a set of points"""
        num_points = len(points)
        lines = [[index, (index + 1) % num_points] for index in range(num_points)]
        return lines

    @classmethod
    def square(
        cls, origin: List, direction: List, scale, manifold: Manifold, delta: float = 0, dtype=int
    ):
        """Create a square (2d) shape on the 2-manifold.

        The square has one corner placed at the origin (origin must be on the plane). The the
        square is defined as the path following heading in the direction <direction> towards the
        next corner of length <length>. The proceeding to the next corner and so on.

        Args:
            origin ([coordinate]): a coordinate on the plane
            direction ([vector]): vector describing the direction of the second corner
            scale ([real]): side of square = scale * ||direction||
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
            scale1=scale,
            scale2=scale,
            manifold=manifold,
            delta=delta,
            dtype=dtype,
        )

    @classmethod
    def rect(
        cls,
        origin: Union[List, Tuple, np.ndarray],
        direction: Union[List, Tuple, np.ndarray],
        scale1: Union[int, float],
        scale2: Union[int, float],
        manifold: Manifold = Plane(),
        delta: float = 0,
        dtype: type = float,
    ) -> "Shape":
        r"""Create a rectangle on a manifold

                      r2
                      o
                  /       \
            r3  o          \
                 \          o r1
                  \       /
                      o
                     r0=origin
            direction (r0 to r1)
            scale1 = scale(r0,r1)
            scale2 = scale(r1,r2)

        Args:
            origin : Origin (r0) of rectangle. The origin is one of the
                corners of the rectangle.
            direction: Direction vector from r0 to r1
            scale1: how large to scale up (r0,r1) side
            scale2: how large to scale up (r1,r2) side
            manifold: Manifold that rectangle lives on. Defaults to Plane().
            delta (optional): How large to increase the rectangle by. Defaults to 0.
            dtype (type, optional): type of entries for definiting points to use.
                Defaults to float.

        Raises:
            QiskitError: Rectangle scales must be positive
            QiskitError: Origin must be on the surface if the manifold
            QiskitError: Manifold not yet supported

        Returns:
            rect: rectangle of given parameters
        """
        origin = np.asarray(origin)
        direction = np.asarray(direction)

        if scale1 < 0 or scale2 < 0:
            raise QiskitError(f"Rectangle scales must be positive: scale1:{scale1} scale2:{scale2}")
        if not manifold.ison(origin):
            raise QiskitError(f"{origin} must be on the surface if the manifold")

        if isinstance(manifold, Plane):
            r0 = origin.astype(dtype)
            r1 = origin + scale1 * direction
            r1 = r1.astype(dtype)
            direction = Plane.rotate(theta=90, vector=direction)
            r2 = r1 + scale2 * direction
            r2 = r2.astype(dtype)
            direction = Plane.rotate(theta=90, vector=direction)
            r3 = r2 + scale1 * direction
            r3 = r3.astype(dtype)

            # scale up the rectangle to have diagonal lenths + 2 * delta
            center = (r1 + r3) / 2
            scale = 1 + delta / np.linalg.norm(r1)
            r0 = scale * r0 + (1 - scale) * center
            r1 = scale * r1 + (1 - scale) * center
            r2 = scale * r2 + (1 - scale) * center
            r3 = scale * r3 + (1 - scale) * center

            points = [r0, r1, r2, r3]
            lines = [[0, 1], [1, 2], [2, 3], [3, 0]]

            return cls(points, lines)

        else:
            raise QiskitError(f"Manifold {manifold} not yet supported")

    @staticmethod
    def _l2distance(p: Tuple, q: Tuple) -> Real:
        return sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2)

    @staticmethod
    def is_between(p: Real, a: Real, b: Real, strict: bool = False, epsilon: Real = 0.0):
        """is p between a and b: a <= p < =b

        Args:
            p: point
            a: end point
            b: end point
            strict: is strict comparison. Default is False
            epsilon: Will squeeze inequality by epsilon

        Returns:
            bool: True is a is between a and b, False if not
        """

        if strict:
            if a == b:
                return False
            if a < b:
                return a + epsilon < p < b - epsilon
            else:
                return b + epsilon < p < a - epsilon
        else:
            if a == b:
                return p == a
            if a < b:
                return a + epsilon <= p <= b - epsilon
            else:
                return b + epsilon <= p <= a - epsilon

    def contains(self, point, on_boundary=True, epsilon=0.01, method="winding") -> bool:
        """

        Args:
            point (_type_): _description_
            on_boundary (bool, optional): _description_. Defaults to True.
            epsilon (float, optional): _description_. Defaults to 0.01.
            method (str, optional): _description_. Defaults to "winding".

        Returns:
            bool: _description_
        """
        if method == "winding":
            return self.contains_quad_winding_number(
                point, on_boundary=on_boundary, epsilon=epsilon
            )
        elif method == "raytracing":
            return self.contains_ray_trace(point, on_boundary=on_boundary, epsilon=epsilon)
        else:
            raise QiskitError(f"Unknown contains method: {method}")

    def contains_quad_winding_number(
        self, point: Union[list, np.ndarray], on_boundary: bool = False, epsilon: float = 0.01
    ):
        """Deterine if a point is inside a polygon (PIP Problem)

        On or near boundary method uses l2 distances and interior method uses
        algorithm from Hornmann and Agathos, Computational Geometry 20 (2001) 131-144

        With optional on_boundary set to False (default) the boundary is not included
        on the inside. If on_boundary is False then points positioned on the boundary
        may of may not be included. Increase the size of the shape to avoid this problem
        or set on_boundary to True. Setting on_boundary to True will be slower as it checks
        l2 distances from boundary.

        Currently only a simple version of the algorithm is implemented. This should be corrected
        at a later date or in next version.

        Args:
            point: Point to check if inside the given shape
            on_boundary: Set True to include the boundary. Slow for larger boundaries.
                Default is False. Points will be included if with epsilon of the boundary.
            epsilon: Value used to indicate if a point is close enough to the boundary to be
                included. Default is 0.01

        Returns:
            bool: True is point in on the polygon and False otherwise.
        """
        x = point[0]
        y = point[1]

        if on_boundary:
            # Determine if the given point is within (l2 distance) epsilon of the shape boundary.
            d = 3 * epsilon
            for index, line in enumerate(self.lines):
                [x0, y0] = self.points[line[0]]
                [x1, y1] = self.points[line[1]]

                if (
                    abs(x0 - x1) < 0.000000001
                ):  # It is assumed that the line sigement is not enormous!
                    # m = float('inf')
                    if Shape.is_between(y, y0, y1):
                        d = min(d, abs(x - x0))
                    else:
                        if y - y0 > y - y1:
                            d = min(d, Shape._l2distance((x, y), (x1, y1)))
                        else:
                            d = min(d, Shape._l2distance((x, y), (x0, y0)))
                else:
                    m = (y0 - y1) / (x0 - x1)
                    c = y0 - m * x0
                    c0 = y0 + m * x0
                    c1 = y1 + m * x1
                    # yp = m * x + c
                    yp0 = -m * x + c0
                    yp1 = -m * x + c1
                    if abs(m) < 0.000000001:
                        if Shape.is_between(x, x0, x1):
                            d = min(d, abs(y - y0))
                        else:
                            if x - x0 > x - x1:
                                d = min(d, Shape._l2distance((x, y), (x1, y1)))
                            else:
                                d = min(d, Shape._l2distance((x, y), (x0, y0)))
                    else:
                        if self.is_between(y, yp0, yp1):
                            d = min(d, abs(m * x - y + c) / sqrt(m**2 + 1))
                        else:
                            if y > yp0 and y > yp1:
                                d = min(d, Shape._l2distance((x, y), (x0, y0)))
                            else:
                                d = min(d, Shape._l2distance((x, y), (x1, y1)))
            if bool(d < epsilon):
                return True

        def det(point, start, end):
            return (start[0] - point[0]) * (end[1] - point[1]) - (end[0] - point[0]) * (
                start[1] - point[1]
            )

        # Determine which quandrant each defining path point is in
        quadrants = [0] * len(self.points)
        for index, line_point in enumerate(self.points):
            if line_point[0] > x and line_point[1] >= y:
                quadrants[index] = 0
            elif line_point[0] <= x and line_point[1] > y:
                quadrants[index] = 1
            elif line_point[0] < x and line_point[1] <= y:
                quadrants[index] = 2
            elif line_point[0] >= x and line_point[1] < y:
                quadrants[index] = 3
        omega = 0
        n = len(quadrants)

        # Calculate winding number
        for index, _ in enumerate(quadrants):
            val = quadrants[(index + 1) % n] - quadrants[index]
            if val in [1, -3]:
                omega = omega + 1
            elif val in [-1, 3]:
                omega = omega - 1
            elif val in [2, -2]:
                val_det = det(point, self.points[index], self.points[(index + 1) % n])
                omega = omega + copysign(2, val_det)
        return bool(omega / 4)

    def contains_ray_trace(self, point, on_boundary=True, epsilon=0.01):
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
        epsilon = 0.01
        x = point[0]
        y = point[1]
        count = 0
        k = len(self.lines)

        # First check if point lines on the boundary. Here a point is defined to lie
        # on the boundary if dist(p,P)<epsilon where p is the point and P is the path
        # of the boundary. When P consists of a collection of line segments then
        # dist(p,P) = min dist(p,l) for l in P where l is a line segement of P.

        def dprint(string):
            if Shape._DEBUG:
                print(string)

        d = 3 * epsilon
        dprint(f"Processing point [{x},{y}]")
        for index, line in enumerate(self.lines):
            dprint(f"Processing line {index}")
            [x0, y0] = self.points[line[0]]
            [x1, y1] = self.points[line[1]]
            dprint(f"{[x0, y0]} to {[x1, y1]}")
            if abs(x0 - x1) < 0.000000001:  # It is assumed that the line sigement is not enormous!
                m = float("inf")
                if on_boundary:
                    if Shape.is_between(y, y0, y1):
                        d = min(d, abs(x - x0))
                    else:
                        if y - y0 > y - y1:
                            d = min(d, Shape._l2distance((x, y), (x1, y1)))
                        else:
                            d = min(d, Shape._l2distance((x, y), (x0, y0)))
            else:
                m = (y0 - y1) / (x0 - x1)
                c = y0 - m * x0
                c0 = y0 + m * x0
                c1 = y1 + m * x1
                yp = m * x + c
                yp0 = -m * x + c0
                yp1 = -m * x + c1
                if on_boundary:
                    if abs(m) < 0.000000001:
                        if Shape.is_between(x, x0, x1):
                            d = min(d, abs(y - y0))
                        else:
                            if x - x0 > x - x1:
                                d = min(d, Shape._l2distance((x, y), (x1, y1)))
                            else:
                                d = min(d, Shape._l2distance((x, y), (x0, y0)))
                    else:
                        if self.is_between(y, yp0, yp1):
                            d = min(d, abs(m * x - y + c) / sqrt(m**2 + 1))
                        else:
                            if y > yp0 and y > yp1:
                                d = min(d, Shape._l2distance((x, y), (x0, y0)))
                            else:
                                d = min(d, Shape._l2distance((x, y), (x1, y1)))

            if m == float("inf"):
                if Shape.is_between(y, y0, y1, epsilon=0.000000001, strict=True) and x < x0:
                    count += 1
                    dprint(f"m=inf inside count = {count}")
            elif abs(m) < 0.000001:
                pass
                dprint(f"m~0 inside count = {count}")
            elif m > 0:
                if Shape.is_between(y, y0, y1, epsilon=0.000000001, strict=True) and y > yp:
                    count += 1
                    dprint(f"m>0 inside count = {count}")
            else:
                if self.is_between(y, y0, y1, epsilon=0.000000001, strict=True) and y < yp:
                    count += 1
                    dprint(f"m<0 inside count = {count}")

            [x2, y2] = self.points[self.lines[(index + 1) % k][1]]
            [_, y3] = self.points[self.lines[(index + 2) % k][1]]

            # Does ray pass though a vertex

            if abs(y - y1) < epsilon:
                if y1 > max(y0, y2) or y1 < min(y0, y2):
                    dprint(r"/\ or \/ vertex count")
                    pass
                elif self.is_between(y1, y0, y2, epsilon=0, strict=True):
                    if x <= x1:
                        count += 1
                        dprint(f"| vertex count = {count}")
                elif abs(y1 - y2) < 0.000001:
                    if y1 > max(y0, y3) or y1 < min(y0, y3):
                        dprint(r"/_\ or \_/ vertex count")
                        pass
                    elif self.is_between(y1, y0, y3, epsilon=0, strict=True):
                        if (x <= x1 < x2) or (x1 > x2 >= x):
                            count += 1
                            dprint(r"/_/ or \_\ vertex count =" + f"{count}")

        dprint(f"d<epsilon = {d < epsilon} and count %2 = {count %2}")
        return bool(d < epsilon) or bool(int(count) % 2)

    def inside(
        self, tiling: Shell, on_boundary: bool = False, epsilon: Real = 0.1
    ) -> Dict[Any, bool]:
        """Returns the in-out data for a shape and a tiling

        Args:
            tiling: Input Tiling (shell)
            on_boundary: If True then the in-out data will include those vertices within epsilon of
                the boundary. This is more expensive. It is better to use a shape which is slightly
                larger. Default is False.
            epsilon: Real number use when on_boundary is True. Default is 0.1

        Returns:
            Dict[bool]: vertex -> bool dictionary. True for inside the shape and False otherwise
        """
        inout = {}
        for vertex in tiling.vertices:
            if Shape._DEBUG:
                print(
                    f"Vertex id: {vertex.id} of face :{vertex.parents[0].parents[0].parents[0].id}"
                )
            inout[vertex] = self.contains(vertex.pos, on_boundary=on_boundary, epsilon=epsilon)
        return inout
