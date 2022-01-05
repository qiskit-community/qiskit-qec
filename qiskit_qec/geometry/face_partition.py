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
    """FacePartition class"""

    def __init__(self) -> None:
        """A class to hold the date for partitioning faces of a shell
        relative to a give shape
        """
        self.faces = {}
        self.parts = None

    def add(self, face: Face):
        """Add a face to the faces dictionary

        Face dictionary : entries [[inside vertices], [outside vertices]]

        Args:
            face (Face): Face to add to dictionary
        """
        self.faces[face] = [[], []]

    def set_parts(self, parts: List[List]):
        """Set parts

        Args:
            parts (List[List]): List of lists of faces with increasing vertex count
        """
        self.parts = parts
