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
# Part of the QEC framework
"""Module for FacePartition"""

from typing import List

from qiskit_qec.geometry.model.face import Face


class FacePartition:
    """FacePartition class
    
        This is a crappy little class to store data on how the faces
        that intersect with the region/cutter intersect
    """

    def __init__(self) -> None:
        """A class to hold the date for partitioning faces of a shell
        relative to a give shape
        """
        self.faces = {}

    def add(self, face: Face):
        """Add a face to the faces dictionary

        Face dictionary : entries [[inside vertices], [outside vertices]]

        Shape.OUTSIDE = 0 -> list of vertices that are outside the region
        Shape.INSIDE = 1 -> list of vertices that are in the region
        Shape.INOUT = 2 -> bool list of which vertices in PATH are In=True or Out=False
        Shape.PATH = 3 ->  list of vertices around face

        self.faces[face][Shape.var] = var list

        Args:
            face (Face): Face to add to dictionary

        """
        self.faces[face] = [[], [], [], []]
