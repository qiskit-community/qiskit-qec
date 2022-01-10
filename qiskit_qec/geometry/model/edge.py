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
"""Edge Module"""
from typing import List

from qiskit_qec.geometry.model.shape_object import ShapeObject
from qiskit_qec.geometry.model.vertex import Vertex


# pylint: disable=invalid-name
class Edge(ShapeObject):
    """Edge"""

    def __init__(self, vertices: List[Vertex]) -> None:
        """Inits Edge

        Args:
            vertices (Tuple[Vertex, Vertex]): Endpoints of the edge.
        """
        super().__init__()
        self.vertices = vertices
        vertices[0].add_parent(self)
        vertices[1].add_parent(self)

    def set_v0(self, v0: Vertex):
        """Sets 0th vertex"""

        self.vertices[0] = v0

    def set_v1(self, v1: Vertex):
        """Sets 1st vertex"""
        self.vertices[1] = v1

    @property
    def v0(self) -> Vertex:
        """Return 0th vertex

        Returns:
            Vertex: 0th vertex
        """
        return self.vertices[0]

    @property
    def v1(self) -> Vertex:
        """Return 1st vertex

        Returns:
            Vertex: 1st vertex
        """
        return self.vertices[1]
