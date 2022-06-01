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
from typing import List, Dict, Any, Tuple

from math import sqrt

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.bounds import GeometryBounds
from qiskit_qec.geometry.manifold import Manifold
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.geometry.model.shell import Shell


# pylint: disable=anomalous-backslash-in-string
class Shape:
    """This class is used to store a boundary shape on a given manifold that is used to
    select sublattices.
    """

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
    def square(cls, origin: List, direction: List, scale, manifold: Manifold, dtype=int):
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
            dtype=dtype,
        )

    @classmethod
    def rect(cls, origin, direction, scale1, scale2, manifold, dtype=int):
        r"""Create a rectangle on a manifold

        .. parsed-literal::
            ```
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
            ```

        Args:
            origin ([]): [description]
            direction ([type]): [description]
            scale1 ([type]): [description]
            scale2 ([type]): [description]
            dtype ([type], optional): [description]. Defaults to int.

        Returns:
            (Shape): rectangle

        Raises:
            QiskitError: qiskit error
        """

        assert isinstance(origin, np.ndarray)
        assert isinstance(direction, np.ndarray)

        assert scale1 > 0 and scale2 > 0, "Rectangle lengths must be positive"
        assert manifold.ison(origin), f"{origin} must be on the surface if the manifold"

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

            points = [r0, r1, r2, r3]
            lines = [[0, 1], [1, 2], [2, 3], [3, 0]]

            return cls(points, lines)

        else:
            raise QiskitError(f"Manifold {manifold} not yet supported")
    
    @staticmethod
    def _l2distance(p:Tuple, q:Tuple)->Real:
        return sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)

    @staticmethod
    def is_between(p, a, b, strict=False, epsilon:Real = 0.0):
        """is p between a and b: a <= p < =b

        Args:
            p (real): point
            a (real): end point
            b (real): end point
            strict (bool): is strict comparison

        Returns:
            bool: True is a is between a and b, False if not
        """
        #print(f"y={p}, y0={a}, y1={b}, epsilon={epsilon}")

        if strict:
            if a == b:
                return False
            if a < b:
                return a+epsilon <= p and p <= b-epsilon
            else:
                return b+epsilon <= p and p <= a-epsilon
        else:
            if a == b:
                return (p ==a)
            if a < b:
                #print("a<b")
                return a+epsilon <= p and p <= b-epsilon
            else:
                #print("b<a")
                return b+epsilon <= p and p <= a-epsilon


        ## pylint: disable=chained-comparison
        #if strict:
        #    return (((a <= b) * a + (b < a) * b) < p) and (p < ((a > b) * a + (b >= a) * b))
        #
        #return (((a <= b) * a + (b < a) * b) <= p) and (p <= ((a > b) * a + (b >= a) * b))

    

    def contains(self, point, on_boundary=True, epsilon=0.1):
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

        # First check if point lines on the boundary. Here a point is defined to lie
        # on the boundary if dist(p,P)<epsilon where p is the point and P is the path 
        # of the boundary. When P consists of a collection of line segments then
        # dist(p,P) = min dist(p,l) for l in P where l is a line segement of P.

        d = 3*epsilon
        for index, line in enumerate(self.lines):
            #print(f"Processing line {index}")
            [x0, y0] = self.points[line[0]]
            [x1, y1] = self.points[line[1]]
            if abs(x0-x1) < 0.0000000000001: # It is assumed that the line sigement is not enormous!
                m = float('inf')
                if on_boundary:
                    if Shape.is_between(y, y0, y1):
                        d = min(d, abs(x-x0))
                    else:
                        if y-y0 > y-y1:
                            d = min(d, Shape._l2distance((x,y),(x1,y1)))
                        else:
                            d = min(d, Shape._l2distance((x,y),(x0,y0)))
            else:
                m = (y0 - y1) / (x0 - x1)
                c  = y0 - m * x0
                c0 = y0 + m * x0
                c1 = y1 + m * x1
                yp = m * x + c
                yp0 = -m * x + c0
                yp1 = -m * x + c1
                if on_boundary:
                    if m == 0:
                        if Shape.is_between(x, x0, x1):
                            d = min(d, abs(y-y0))
                        else:
                            if x - x0 > x - x1:
                                d = min(d, Shape._l2distance((x,y),(x1,y1)))
                            else:
                                d = min(d, Shape._l2distance((x,y),(x0,y0)))
                    else:
                        if self.is_between(y, yp0, yp1):
                            d = min(d, abs(m * x - y + c)/sqrt(m**2+1))
                        else:
                            if y > yp0 and y > yp1:
                                d=min(d,Shape._l2distance((x,y),(x0,y0)))
                            else:
                                d=min(d,Shape._l2distance((x,y),(x1,y1)))

            if m == float('inf'):
                if Shape.is_between(y, y0, y1, epsilon=epsilon) and x < x0:
                    count += 1
                    print(f"0 count={count}")
            elif m == 0:
                pass
            elif m > 0:
                if Shape.is_between(y, y0, y1, epsilon=epsilon) and y > yp:
                    count += 1
            else:
                if self.is_between(y, y0, y1, epsilon=epsilon) and y < yp:
                    count += 1

            [_, y2] = self.points[self.lines[(index + 1) % k][1]]

            # Does ray pass though a vertex

            if abs(y - y1) < epsilon:
                if x < x1:
                    if y1 > max(y0, y2) or y1 < min(y0, y2):
                        # ray passes through vertex y1 /\ or \/
                        pass
                    else:
                        # ray passes through vertex y1 |
                        count += 1

        return bool(d < epsilon) or bool(count % 2)

    
    def intersection(self, tiling:Shell, on_boundary=False) -> Dict[Any, bool]:
        """Returns the inout data for a shape and a tiling

        Args:
            tiling (Shell): _description_

        Returns:
            Dict[bool]: _description_
        """
        inout = {}
        for vertex in tiling.vertices:
            inout[vertex] = self.contains(vertex.pos, on_boundary=on_boundary)
        return inout

