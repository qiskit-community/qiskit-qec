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

from qiskit_qec.geometry.model.shape_object import ShapeObject

class Edge(ShapeObject):
    def __init__(self, vertices) -> None:
        self.vertices = vertices
        super().__init__(stype=Edge, child=self)
        vertices[0].add_parent(self)
        vertices[1].add_parent(self)

    def set_v0(self, v0):
        self.vertices[0] = v0

    def set_v1(self, v1):
        self.vertices[1] = v1

    @property
    def v0(self):
        return self.vertices[0]
    
    @property
    def v1(self):
        return self.vertices[1]
    

