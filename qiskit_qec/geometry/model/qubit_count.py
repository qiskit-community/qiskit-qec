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

class QubitCount:
    def __init__(self) -> None:
        self.num_qubits = 0
        self.qubits_count = dict()

    def new_qubit(self):
        self.qubits_count[self.num_qubits] = 0
        self.num_qubits += 1
        return self.num_qubits - 1

    def increment_qubit(self, key):
        self.qubits_count[key] += 1

    def decrement_qubit(self, key):
        self.qubits_count[key] -= 1

    @property
    def qubits(self):
        keys = self.num_qubits.keys()
        qubits = [item for item in keys if item>0]
        return qubits

