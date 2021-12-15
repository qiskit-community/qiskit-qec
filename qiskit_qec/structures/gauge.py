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

import copy

from qiskit.exceptions import QiskitError

from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.operators.pauli import Pauli
from qiskit_qec.info.properties import Properties
from qiskit_qec.structures.group import Group

class GaugeGroup(Group):
    def __init__(self, 
                generators=None, 
                gauge_group=None, 
                code_id=None, 
                code_family=None, 
                code_parameters=None,
                deep_copy=False) -> None:
        if generators is not None:
            try:
                self.generators = PauliList(generators)
            except:
                raise QiskitError("Error: generators cannot generate PauliList")
        elif gauge_group is not None:
            if isinstance(gauge_group, GaugeGroup):
                if deep_copy == False:
                    self.generators = gauge_group.generators
                else:
                    self.generators = gauge_group.generators.copy()
                self.properties = Properties.duplicate(gauge_group.properties)
        elif code_id is not None:
            raise QiskitError("Note yet implemeneted")
        elif code_family is not None:
            raise QiskitError("Note yet implemeneted")
        super().__init__()

    def from_id(code_id=None):
        raise QiskitError("from_id not yet implemented")

    def from_family(code_family=None, code_parameters=None):
        raise QiskitError("from_family not yet implemented")