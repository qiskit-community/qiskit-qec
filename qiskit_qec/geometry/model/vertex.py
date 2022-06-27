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
"""Module for Vertex"""
from typing import List, Union

from qiskit_qec.geometry.model.shape_object import ShapeObject


class Vertex(ShapeObject):
    """ "Vertex" inherits from "ShapeObject" """

    def __init__(self, pos: List[Union[float, int]]) -> None:
        """Inits Vertex

        Args:
            pos (List[Union[float, int]]): position of vertex
        """
        self.pos = pos
        super().__init__()

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        return self.pos.__str__()

    def set_position(self, pos: List[Union[float, int]]):
        """Sets global position of Vertex

        Args:
            pos (List[Union[float, int]]): global position
        """
        self.pos = pos

    @property
    def position(self) -> List[Union[float, int]]:
        """Return position of vertex"""
        return self.pos

    def shallowcopy(self) -> "Vertex":
        """Returns a shallow copy of the Vertex.

        A shallow copy creates a new Vertex with only the position being copied

        Returns:
            Vertex:  Shallow copy of vertex
        """
        return Vertex(self.pos)
