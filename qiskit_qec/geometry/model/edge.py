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
from typing import List, Optional

from qiskit_qec.geometry.model.shape_object import ShapeObject
from qiskit_qec.geometry.model.vertex import Vertex


# pylint: disable=invalid-name
class Edge(ShapeObject):
    """Edge"""

    def __init__(
        self,
        vertices: List[Vertex],
        next_edge: Optional["Edge"] = None,
        previous_edge: Optional["Edge"] = None,
    ) -> None:
        """Inits Edge

        Args:
            vertices (Tuple[Vertex, Vertex]): Endpoints of the edge.
            next_edge: Next edge
            previous_edge: Previous edge
        """
        super().__init__()
        self.vertices = vertices
        self.coincident_edges = []
        self.next = next_edge
        self.previous = previous_edge

        for item in self.vertices:
            item.add_parent(self)

    def __repr__(self) -> str:
        string = "Edge" + self.__str__()
        return string

    def __str__(self):
        string = "["
        for vertex in self.vertices[:-1]:
            string += vertex.__str__()
            string += ","
        string += self.vertices[-1].__str__()
        string += "]"
        return string
