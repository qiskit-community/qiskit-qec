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
"""Module of Wireframe"""
from typing import List

from qiskit_qec.geometry.model.edge import Edge
from qiskit_qec.geometry.model.shape_object import ShapeObject


class WireFrame(ShapeObject):
    """`WireFrame` inherits from `ShapeObject`"""

    def __init__(self, edges: List[Edge]) -> None:
        """Inits WireFrame

        A wireframe should have only a single connected component (but are not
        linited to having a single component)

        Args:
            edges (List[Edge]): Edges that comprise the WireFrame
        """
        super().__init__()

        self.edges = edges
        self.vertices = []

        # Recording edges and vertices for this wireframe
        for edge in self.edges:
            edge.add_parent(self)
            for vertex in edge.vertices:
                if vertex not in self.vertices:
                    self.vertices.append(vertex)
