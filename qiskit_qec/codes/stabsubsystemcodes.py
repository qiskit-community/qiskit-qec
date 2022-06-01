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
""" This is the main module that defines a StabSubsystem Code. """

from typing import Union, Tuple, List

from qiskit_qec.codes.code import Code
from qiskit_qec.structures.gauge import GaugeGroup


class StabSubSystemCode(Code):
    """The "StabSubSystemCode" inherits the "Code" class"""

    def __init__(self, gauge_group: GaugeGroup) -> None:
        """SubSystemCode

        Args:
            gauge_group (GaugeGroup)

        Raises:
            QiskitError: Error
        """

        self.gauge_group = gauge_group
        self._n = self.gauge_group.generators.num_qubits  # pylint: disable=invalid-name
        super().__init__()

    @property
    def n(self) -> int:
        """_summary_

        Returns:
            int: _description_
        """
        return self._n

    def k(self) -> int:
        """_summary_

        Returns:
            int: _description_
        """
        pass

    def gauge_degree(self) -> int:
        """_summary_

        Returns:
            int: _description_
        """
        pass

    def id(self) -> Union[Tuple[int, int, int], None]:  # pylint: disable=invalid-name
        """_summary_

        Returns:
            Union[Tuple[int,int,int],None]: _description_
        """
        pass

    def d(self) -> int:
        """_summary_

        Returns:
            int: _description_
        """
        pass

    def dx(self) -> int:  # pylint: disable=invalid-name
        """_summary_

        Returns:
            int: _description_
        """
        pass

    def dz(self) -> int:  # pylint: disable=invalid-name
        """_summary_

        Returns:
            int: _description_
        """
        pass

    def dy(self) -> int:  # pylint: disable=invalid-name
        """_summary_

        Returns:
            int: _description_
        """
        pass

    def is_css(self) -> bool:
        """_summary_

        Returns:
            bool: _description_
        """
        pass

    def is_decomposable(self) -> bool:
        """_summary_

        Returns:
            bool: _description_
        """
        pass

    def is_degenerate(self) -> bool:
        """_summary_

        Returns:
            bool: _description_
        """
        pass

    def is_gf4linear(self) -> bool:
        """_summary_

        Returns:
            bool: _description_
        """
        pass

    def is_triorthogonal(self) -> bool:
        """_summary_

        Returns:
            bool: _description_
        """
        pass

    def weight_enumerator(self) -> List[int]:
        """_summary_

        Returns:
            List[int]: _description_
        """
        pass
