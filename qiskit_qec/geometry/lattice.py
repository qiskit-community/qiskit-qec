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
"""Module for lattice"""

from typing import Union, List, Optional

from math import ceil, floor

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.bounds import GeometryBounds
from qiskit_qec.geometry.shape import Shape
from qiskit_qec.geometry.tiles.tile import Tile


class Lattice:
    """Lattice on which tiles are tiled"""

    def __init__(self, u_vec=None, v_vec=None, size=None, points=None) -> None:
        """Lattice

        Args:
            u_vec ([type], optional): Zeroth basis vector. Defaults to None.
            v_vec ([type], optional): First basis vector. Defaults to None.
            size ([type], optional): Width/Height of lattice. Defaults to None.
            points ([type], optional): Points on generated lattice. Defaults to None.


        Raises:
            QiskitError: Something went wrong.
        """

        if u_vec is None:
            u_vec = np.array([1, 0])
        if v_vec is None:
            v_vec = np.array([0, 1])
        if size is None:
            size = (np.inf, np.inf)

        self.u_vec = u_vec
        self.v_vec = v_vec
        self.transform = Lattice.make_transform(u_vec, v_vec)
        self.unorm = np.linalg.norm(u_vec)
        self.vnorm = np.linalg.norm(v_vec)

        if size != (np.inf, np.inf):
            if (size[0] == np.inf) ^ (size[1] == np.inf):
                raise QiskitError("Half infinite lattices not yet supported")
            if points is None:
                self.points = self.generate_points(np.array(size))
        else:
            self.points = points

    def __str__(self) -> str:
        if self.points is None:
            return f"Lattice(u_vec={np.array2string(self.u_vec, separator=', ')},\
                 v_vec={np.array2string(self.v_vec, separator=', ')})"
        outstr = "Lattive["
        for point in self.points[:-1]:
            outstr += np.array2string(point, separator=", ")
            outstr += ", "
        outstr += np.array2string(self.points[-1], separator=", ")
        outstr += "]"
        return outstr

    @classmethod
    def make_transform(cls, u_vec, v_vec):
        """Generate a transformation to be used by other lattices"""
        return np.vstack((u_vec, v_vec))

    def find_pre_transform_length(self, size):
        """Shear length.

        Args:
            size (int|float): size
        """
        xynorm = np.linalg.norm(size)
        dem1 = min(self.unorm, self.vnorm)
        dem2 = np.sqrt(1 - np.dot(self.u_vec, self.v_vec) / (self.unorm * self.vnorm))
        length = ceil(xynorm / (dem1 * dem2))
        length += length % 2
        return length

    def generate_points(self, size) -> List:
        """lattice of size points centered at (0,0)

        Args:
            size (int|float): size vector for lattice of points
        """
        points = []
        bounds = np.ceil(size / 2).astype(int)
        for i in range(-bounds[0], bounds[0] + 1):
            for j in range(-bounds[1], bounds[1] + 1):
                points.append(i * self.u_vec + j * self.v_vec)

        return points

    def apply_transform_from(self, lattice: "Lattice") -> "Lattice":
        """Apply transformation to self from lattice"""
        points = []
        if self.points is None:
            raise QiskitError("Lattice points must first be generated")
        for point in self.points:
            points.append(np.matmul(point, lattice.transform))

        return Lattice(lattice.u_vec, lattice.v_vec, points=points)

    def restrict(self, region: GeometryBounds, *, in_place: bool = False) -> "Lattice":
        """Create a new lattice with same basis but that fits inside the supplied region

        Args:
            region (GeometryBounds): region to which the lattice will be restricted.
            in_place (bool): Modify lattice in place. Default is False

        Raises:
            QiskitError: something went wrong

        Returns:
            Lattice : A new lattice with same basis vectors as self but
                restricted to bounding_box
        """

        if isinstance(region, GeometryBounds):
            if self.points is None:
                raise QiskitError("Points must first be generated")
            points = []
            for point in self.points:
                if region.contains(point):
                    points.append(point)

            if in_place:
                self.points = points
                return self

            return Lattice(self.u_vec, self.v_vec, points=points)
        else:
            raise QiskitError("Region should be an instance of GeometryBounds")

    def restrict_for_tiling(
        self,
        region: Shape,
        *,
        tile: Optional[Tile] = None,
        size: List[Union[float, int]] = None,
        expand_value: np.array = None,
        in_place: bool = False,
        alpha: float = 1,
    ) -> List:
        """Given a Shape to tile based on the lattice (self), restrict lattice (self) to
        the provided shape such that a tiling of that shape with the given tile will
        completely fill that shape. A size can be provided instead of a tile. If both are
        provided the size attribute will be used. If no Tile or size is provided then
        it will be assumed that a tiling will use tiles of the size of the lattices (self)
        given basis.

        Args:
            region (Shape): region that needs to be tiled
            tile (Tile, optional): Tile to be used for tiling. Defaults to None.
            size (List[Union[float, int]], optional): Size of tile. Defaults to None.
            expand_value (np.array, optional): Amount to expand the region AABB by
            to ensure that the entire region is filled. If no expand_value is provied then
            an approximate value will be determined. Defaults to None.
            in_place (bool, optional): Perform in place on lattice if set to True. Defaults to False.
            alpha: tile_size = alpha * tile.size - Used for optimization of a factory

        Returns:
            Union[Lattice, None]: Lattice for tiling region, or None if done in place.
        """

        # Create an extended bounding box AABB
        extended_aabb = region.bounds.copy()

        # Determine the size of tile expected

        if tile is not None:
            tile_size = alpha * tile.size
        elif size is not None:
            tile_size = np.array(size)
        else:
            tile_size = np.array((self.unorm, self.vnorm))

        # Determine how much to expand the region's AABB by

        if expand_value is None:
            expand_value = tile_size

        extended_aabb.expand(expand_value)

        # Choose vert_vec so that it has some y component and that horz_vec has some x component
        if abs(self.v_vec[1]) > 0.1:
            if abs(self.u_vec[0]) > 0.1:
                vert_vec = self.v_vec
                horz_vec = self.u_vec
            else:
                vert_vec = self.u_vec
                horz_vec = self.v_vec
        elif abs(self.u_vec[1]) > 0.1:
            if abs(self.v_vec[0]) > 0.1:
                vert_vec = self.u_vec
                horz_vec = self.v_vec
            else:
                raise QiskitError(
                    "Lattice basis not independent or too skewed or vector lengths \
            to small. Resize and/or run LLL to improve basis before continuing"
                )

        assert horz_vec[0] != 0, "Division by zero - error in lattice basis selection"

        def _find_points(aabb, u_vec, v_vec, point, pts, direction="up"):
            while True:
                if direction == "up":
                    point = point + v_vec
                else:
                    point = point - v_vec
                m = u_vec[1] / u_vec[0]
                c = point[1] - m * point[0]
                line = [-m, 1, c]
                intersects = aabb.intercepts(line)
                if len(intersects) != 2:
                    break
                x0, y0 = intersects[0]
                x1, y1 = intersects[1]

                min_x = min(x0, x1)
                max_x = max(x0, x1)

                ta_min = ceil((min_x - point[0]) / u_vec[0])
                ta_max = floor((max_x - point[0]) / u_vec[0])

                if u_vec[1] != 0:
                    min_y = min(y0, y1)
                    max_y = max(y0, y1)

                    tb_min = ceil((min_y - point[1]) / u_vec[1])
                    tb_max = floor((max_y - point[1]) / u_vec[1])

                    t_min = max(ta_min, tb_min)
                    t_max = min(ta_max, tb_max)
                else:
                    t_min = ta_min
                    t_max = ta_max

                t_range = range(t_min, t_max + 1)
                for t_val in t_range:
                    new_point = point + t_val * u_vec
                    pts.append(new_point.copy())

            return pts

        # find points on lattice on vert_vec line that are in AABB (and one just outside)
        origin = np.array((0, 0))
        point = origin - vert_vec
        points = []
        points = _find_points(extended_aabb, horz_vec, vert_vec, point, points, direction="up")
        point = origin.copy()
        points = _find_points(extended_aabb, horz_vec, vert_vec, point, points, direction="down")

        # Create a lattice with points generated

        lattice_l = Lattice(u_vec=horz_vec, v_vec=vert_vec, points=points)

        if in_place:
            self.points = lattice_l.points
            return self

        return lattice_l
