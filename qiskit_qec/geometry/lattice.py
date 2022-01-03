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

from math import ceil
import numpy as np

from qiskit.exceptions import QiskitError

class Lattice:
    def __init__(self, u=np.array([1,0]), v=np.array([0,1]), size=(np.inf, np.inf), points=None) -> None:
        self.u = u
        self.v = v
        self.transform = Lattice.make_transform(u,v)
        self.unorm = np.linalg.norm(u)
        self.vnorm = np.linalg.norm(v)

        if size != (np.inf, np.inf):
            if (size[0] == np.inf) ^ (size[1]==np.inf):
                raise QiskitError("Half infinite lattices not yet supported")
            else:
                if points is None:
                    self.points = self.generate_points(np.array(size))
        else:
            self.points = points

    def make_transform(u,v):
        return np.vstack((u,v))

    def find_shear_length(self, size):
        xynorm = np.linalg.norm(size)
        dem1 = min(self.unorm, self.vnorm)
        dem2 = np.sqrt(1-np.dot(self.u, self.v)/(self.unorm * self.vnorm))
        l = ceil(xynorm/(dem1*dem2)) 
        l += (l%2)
        return l

    def generate_points(self, size):
        """lattice of size points centered at (0,0)

        Args:
            size: size vector for lattice of points
        """
        points = []
        bounds = np.ceil(size/2).astype(int) 
        for i in range(-bounds[0],bounds[0]+1):
            for j in range(-bounds[1], bounds[1]+1):
                points.append(i * self.u + j * self.v)

        return points

    def apply_transform_from(self, lattice):
        points = []
        if self.points is None:
            raise QiskitError("Lattice points must first be generate")
        for point in self.points:
            points.append(np.matmul(point,lattice.transform))

        return Lattice(lattice.u, lattice.v, points=points)


    def restrict(self, bounding_box):
        if self.points is None:
            raise QiskitError("Points must first be generated")

        points = []
        for point in self.points:
            if bounding_box.contains(point):
                points.append(point)

        return Lattice(self.u, self.v, points=points)



