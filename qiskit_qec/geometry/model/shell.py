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
"""Module for Shell"""

from typing import List

# import numpy as np

from qiskit_qec.geometry.model.shape_object import ShapeObject
from qiskit_qec.geometry.model.face import Face


class Shell(ShapeObject):
    """`Shell` inherits from `ShapeObject`"""

    def __init__(self, faces: List[Face]) -> None:
        """Inits Shell.
        Each shell keeps track internally of several lists for
        each type of subcomponent (ex: faces, wireframes, edges, vertices)
        Args:
            faces (List[Face]): Faces that make up the `Shell` instance
        """
        super().__init__()

        self.faces = faces
        self.wireframes = []
        self.edges = []
        self.vertices = []

        for face in self.faces:
            self.wireframes.append(face.wireframe)
            self.edges += face.edges
            self.vertices += face.vertices

    def union(self, other_shell: "Shell"):
        # TODO control for overlapping subcomponents of other_shell and self
        """Add all subcomponents of other_shell to self

        Args:
            other_shell (Shell): A different shell
        """
        self.faces += other_shell.faces
        self.wireframes += other_shell.wireframes
        self.edges += other_shell.edges
        self.vertices += other_shell.vertices

    def delete_subtree(self, roots: ShapeObject):
        """Delete subcomponents of self"""
        pass

    # @property
    # def size(self):
    #    """Get the size set of vertex
    #
    #    Returns:
    #        [type]: [description]
    #    """
    #    vert = np.asarray(self.vertices)
    #    return np.amax(vert)
