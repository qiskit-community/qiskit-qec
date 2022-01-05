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
""" This is the main module that defines a Subsystem Code. """

from qiskit.exceptions import QiskitError
from qiskit_qec.codes.code import Code
from qiskit_qec.structures.gauge import GaugeGroup


class SubSystemCode(Code):
    """The "SubSystemCode" inherits the "Code" class"""

    def __init__(
        self,
        gauge_group: GaugeGroup = None,
        parameters: dict = None,
    ) -> None:
        """SubSystemCode

        Args:
            gauge_group (GaugeGroup, optional): Defaults to None.
            parameters (dict, optional): Additional parameters user may create
            on the fly when designing codes. Defaults to None.

        Raises:
            QiskitError: Error
        """
        if gauge_group is not None:
            self.gauge_group = GaugeGroup(gauge_group=gauge_group)
        else:
            raise QiskitError("Error: Input format not supported")

        self.parameters = parameters

        super().__init__()
