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
"""Module for Face"""
from typing import List

from qiskit_qec.geometry.model.shape_object import ShapeObject
from qiskit_qec.geometry.model.wireframe import WireFrame


class Face(ShapeObject):
    """`Face` inherits from `ShapeObject`"""

    def __init__(self, wireframes: List[WireFrame]) -> None:
        """Inits Face

        Args:
            wireframe (WireFrame): Wireframe for face
        """
        super().__init__()
        self.wireframes = wireframes
        self.edges = []
        self.vertices = []
        for wf in self.wireframes:
            self.edges += wf.edges
            self.vertices += wf.vertices
            wf.add_parent(self)
