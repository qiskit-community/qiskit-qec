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
from typing import Union, List, Tuple
from numpy import ndarray
from qiskit import QiskitError


class QubitData:
    """Class for containing qubit information"""

    def __init__(self) -> None:
        """Init Qubit Data"""

        # Vertex Data
        self.operator = {}  # Vertex_id to operator list for given vertex
        self.qubit = {}  # vertex_id to qubit_id

        self.index = {}  # vertex_id to PauliList index
        # Edge Data

        # Wireframe Data

        # Face Data
        self.face_colors = {}  # Face id to color str

        # Other

        self.qubit_to_index = {}  # qubit_id to PauliList index
        self.index_to_qubit = {}  # PauliList index to qubit_id

        self.data_arrays = {}

    def add_data_array(self, data_array: Union[List, Tuple, ndarray], name: str) -> None:
        """Adds a data array to the QubitData class instance

        Args:
            data_array (Union[List, Tuple, ndarray]): _description_
            name: string name of array
        """
        self.data_arrays[name] = data_array

    def del_data_array(self, name: str):
        """Deletes a given data_array from the QubitData class instance

        Args:
            name (str): Name of data array to delete
        """
        try:
            del self.data_arrays[name]
        except KeyError as keyerror:
            raise QiskitError(f"Data array {name} does not exist") from keyerror
