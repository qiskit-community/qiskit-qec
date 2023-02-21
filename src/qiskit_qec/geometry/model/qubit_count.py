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
"""Module for keeping track of qubits in a given geometry"""


class QubitCount:
    """Each geometry will have a QubitCount class to maintain pointers to all qubits currently in use"""

    def __init__(self) -> None:
        """QubitCount inits with 0 qubits and an empty qubit:references dictionary"""
        self.num_qubits = 0
        self.qubits_count = {}

    def new_qubit(self) -> int:
        """Creates a new qubits_count dictionary entry.
        The key is ID of qubit. The value is the reference.

        Returns:
            int: Qubit ID
        """
        self.qubits_count[self.num_qubits] = 0
        self.num_qubits += 1
        return self.num_qubits - 1

    def increment_qubit(self, key: int) -> None:
        """Increment number of references to qubit with ID key

        Args:
            key (int): Unique ID of qubit
        """

        self.qubits_count[key] += 1

    def decrement_qubit(self, key: int) -> None:
        """Decrement number of references to qubit with ID key

        Args:
            key (int): Unique ID of qubit
        """

        self.qubits_count[key] -= 1
