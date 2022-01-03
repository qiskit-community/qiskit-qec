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

import numpy as np

from qiskit_qec.geometry.plane import Plane
from qiskit_qec.geometry.bounds import GeometryBounds

from qiskit.exceptions import QiskitError


class Shape:
    """This class is used to store a boundary shape on a given manifold that is used to
    select sublattices. 
    
    """
    def __init__(self, points, lines, indices=None):
        self.points = points
        self.lines = lines
        if indices is not None:
            self.indices = indices
        else:
            self.indices = list(range(len(self.points)))

        self.bounds = self.bounding_box_from_lines()

    def bounding_box_from_lines(self):
        line = self.lines[0]
        start = self.points[line[0]]
        finish = self.points[line[1]]
        bounds = GeometryBounds.bounding_box_from_line(start, finish)
        for line in self.lines:
            bounds2 = GeometryBounds.bounding_box_from_line(self.points[line[0]], self.points[line[1]])
            bounds = GeometryBounds.combine(bounds,bounds2)
        return bounds

    @classmethod
    def square(cls, 
        origin, 
        direction, 
        length, 
        manifold,
        dtype=int):

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
            [Shape]: Shape object describing a square of length <length> in the manifold <manifold>
        """

        origin = np.asarray(origin)
        direction = np.asarray(direction)

        return cls.rect(
            origin=origin, 
            direction=direction, 
            length1=length, 
            length2=length,
            manifold=manifold,
            dtype=dtype)

    @classmethod
    def rect(cls, 
        origin, 
        direction, 
        length1, 
        length2, 
        manifold,
        dtype=int):

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

        Raises:
            QiskitError: [description]
        """         

        assert(isinstance(origin, np.ndarray))
        assert(isinstance(direction, np.ndarray))

        assert length1>0 and length2>0, "Rectangle lengths must be positive"
        assert manifold.ison(origin), f"{origin} must be on the surface if the manifold"

        if isinstance(manifold, Plane):
            zero_origin = np.zeros(2,dtype=dtype)
            r0 = origin.astype(dtype)
            r1 = origin + length1 * direction
            r1 = r1.astype(int)
            direction = Plane.rotate(theta=90, vector=direction)
            r2 = r1 + length2 * direction
            r2 = r2.astype(int)
            direction = Plane.rotate(theta=90, vector=direction)
            r3 = r2 + length1 * direction
            r3 = r3.astype(int)
            
            points = [r0,r1,r2,r3]
            lines = [[0,1],[1,2],[2,3],[3,0]]

            return cls(points, lines)

        else:
            raise QiskitError(f"Manifold {manifold} not yet supported")

    
    def contains(self, point):
        """Simple check if inside region using an infinite horizontal line from point
        point to +infinity

        This simple implementation assumes that the path is given by a sequence of 
        straight lines.

        A more sophisticated verison is required if the path is otherwise.

        Args:
            point (point): point to check if in region (including on boundary)
        """
        x = point[0]
        y = point[1]
        print(f"Check {[x,y]}")
        count = 0
        for line in self.lines:
            [x0, y0] = self.points[line[0]]
            [x1, y1] = self.points[line[1]]
            if  x <= min(x0, x1) and Shape.is_between(y, y0, y1):
                count += 1
                print(f"{[x0, y0]} to {[x1, y1]} : Intersect")
            else:
                print(f"{[x0, y0]} to {[x1, y1]} : Not intersect")

        return bool(count%2)

    def is_between(p, a, b):
        """is p between a and b: a <= p < =b

        Args:
            p (real): point
            a (real): end point
            b (real): end point

        Returns:
            bool: True is a is between a and b, False if not
        """
        return (((a <= b)*a +(b<a)*b) <= p) and (p <= ((a>b)*a + (b>=a)*b))