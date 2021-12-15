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

from qiskit.exceptions import QiskitError

from qiskit_qec.structures.gauge import GaugeGroup
from qiskit_qec.info.properties import Properties
from qiskit_qec.codes.code import Code

class SubSystemCode(Code):
    def __init__(self, 
                gauge_group=None, 
                code_id=None, 
                code_family=None, 
                code_parameters=None, 
                parameters=None) -> None:
        """SystemSystemCode class

        """

        if gauge_group is not None:
            self.gauge_group = GaugeGroup(gauge_group=gauge_group)
        elif code_id is not None:
            self.gauge_group = GaugeGroup.from_id(code_id=code_id)
        elif code_family is not None:
            self.gauge_group = GaugeGroup.from_family(
                code_family=code_family, 
                code_parameters=code_parameters)
        else:
            raise QiskitError("Error: Input format not supported")


        self.properties = Properties(__class__.__name__, parameters)

        self.model_array = dict()

        super().__init__()

    def add_model(self, model, key):
        self.model_array[key] = model

    def get_model(self, key):
        return self.model_array[key]







    