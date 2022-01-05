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

class Vertex(ShapeObject):
    def __init__(self, pos) -> None:
        self.pos = pos
        super().__init__(stype=Vertex, child=self)

    def set_position(self, pos):
        self.pos = pos

    @property
    def position(self):
        return self.pos

