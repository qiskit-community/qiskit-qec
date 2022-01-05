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
"""Module for Gauge Groups"""

from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.structures.group import Group


class GaugeGroup(Group):
    """`GaugeGroup` inherits from `Group`"""

    def __init__(
        self,
        generators: PauliList = None,
        gauge_group=None,
        deep_copy=False,
    ) -> None:
        """[summary]

        Args:
            generators (PauliList, optional): Generator for the GaugeGroup. Defaults to None.
            gauge_group (GaugeGroup, optional): GaugeGroup to either copy or deep copy depending
            on `deep_copy` param. Defaults to None.
            deep_copy (bool, optional): Whether to deep or shallow copy gauge_group. Defaults to False.

        Raises:
            QiskitError: Something went wrong.
        """
        super().__init__()
        if generators is not None:
            self.generators = PauliList(generators)
        elif gauge_group is not None:
            if isinstance(gauge_group, GaugeGroup):
                if deep_copy is False:
                    self.generators = gauge_group.generators
                else:
                    self.generators = gauge_group.generators.copy()
