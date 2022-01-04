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
import numpy as np

from qiskit_qec.geometry.model.shape_object import ShapeObject


class Shell(ShapeObject):
    def __init__(self, faces) -> None:
        self.faces = faces
        self.wireframes = []
        self.edges = []
        self.vertices = []
        super().__init__(stype=Shell, child=self)

        for face in self.faces:
            self.wireframes.append(face.wireframe)
            self.edges += face.edges
            self.vertices += face.vertices

    def union(self, shell):
        self.faces += shell.faces
        self.wireframes += shell.wireframes
        self.edges += shell.edges
        self.vertices += shell.vertices

    def delete_subtree(self, roots):
        pass
            

