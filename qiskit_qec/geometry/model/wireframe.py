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

class WireFrame(ShapeObject):
    def __init__(self, edges) -> None:
        self.edges = edges
        self.vertices = []
        super().__init__(stype=WireFrame, child=self)

        # Recording edges and vertices for this wireframe
        for edge in self.edges:
            edge.add_parent(self)
            for vertex in edge.vertices:
                if vertex not in self.vertices:
                    self.vertices.append(vertex)