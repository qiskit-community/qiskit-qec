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

from typing import Union, List

from math import ceil

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

    @classmethod
    def make_transform(cls, u_vec, v_vec):
        """Generate a transformation to be used by other lattices"""
        return np.vstack((u_vec, v_vec))

    def find_pre_transform_length(self, size):
        """Shear length.

        Args:
            size (int|float): size
        """
        # TODO I don't know why we need size
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
        # TODO: Not sure what this actually does
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
        tile: Tile = None,
        size: List[Union[float, int]] = None,
        expand_value: np.array = None,
        in_place: bool = False,
    ) -> "Lattice":
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

        Returns:
            Union[Lattice, None]: Lattice for tiling region, or None if done in place.
        """

        # Create an extended bounding box AABB
        extended_aabb = region.bounds.copy()

        # Determine the size of tile expetected

        if tile is not None:
            tile_size = tile.size
        elif size is not None:
            tile_size = np.array(size)
        else:
            tile_size = np.array((self.unorm, self.vnorm))

        # Determine how much to expand the region's AABB by

        if expand_value is None:
            expand_value = tile_size

        extended_aabb.expand(expand_value)

        # Calculate pre-transform length
        _len = self.find_pre_transform_length(extended_aabb.size)

        # Create standard _len x _len lattice to be transformed
        standard_lattice_l = Lattice(size=(_len, _len))

        # Transform the standard lattice into a sublattice of self
        lattice_l = standard_lattice_l.apply_transform_from(self)

        # Mask the lattice with the expanded AABB
        lattice_l.restrict(extended_aabb, in_place=True)

        if in_place:
            self.points = lattice_l.points
            return self

        return lattice_l
