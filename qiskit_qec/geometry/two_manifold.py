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

from abc import abstractclassmethod
from qiskit.exceptions import QiskitError

from qiskit_qec.geometry.manifold import Manifold

class TwoManifold(Manifold):
    def __init__(self):
        dim=2
        super().__init__(dim=2)
    
    @abstractclassmethod
    def ison(self, point):
        pass
    