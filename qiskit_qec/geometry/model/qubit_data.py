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
"""Module for Qubit Data"""


class QubitData:
    """Class for containing qubit information"""

    def __init__(self) -> None:
        # TODO do QubitData and QubitCount share a qubit ID?
        """Init Qubit Data"""
        self.operator = {} # Vertex_id to 
        self.qubit = {} # vertex_id to qubit_id
        self.orientation = {} # face_id -> orientation
