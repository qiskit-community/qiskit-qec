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

from typing import Optional

from qiskit import QiskitError

from qiskit_qec.codes.code import Code
from qiskit_qec.geometry.model.qubit_count import QubitCount
from qiskit_qec.geometry.model.qubit_data import QubitData
from qiskit_qec.geometry.model.shell import Shell

from qiskit_qec.structures.gauge import GaugeGroup
from qiskit_qec.operators.pauli_list import PauliList


class StabSubSystemCode(Code):
    """The "StabSubSystemCode" inherits the "Code" class"""

    def __init__(
        self,
        gauge_group: Optional[GaugeGroup] = None,
        *,
        shell: Optional[Shell] = None,
        qubit_data: Optional[QubitData] = None,
        qubit_count: Optional[QubitCount] = None,
    ) -> None:
        """Stabilizer Subsystem System Code Class

        A stabilizer subsystem code can be initialized by providing either a gauge group
        of a shell and associated data but not both.

        Args:
            gauge_group (optional): Gauge group defining the subsystem. Defaults to None.
            shell (optional): Shell definging the subsystem. Defaults to None.
            qubit_data (optional): Qubit data associated with defining shell. Defaults to None.
            qubit_count (optional): Qubit count data associated with defining shell. Defaults to None.

        Examples:
            >>> generators = PauliList(['X1X2','Z3Z4', 'X3Z9'])
            >>> gauge_group = GaugeGroup(generators)
            >>> code = StabSubSystemCode(gauge_group)

            >>> qubit_data = QubitData()
            >>> qubit_count = QubitCount()
            >>> origin = numpy.array([0,0])
            >>> shell = CheckerBoardTile(origin=origin,
                                         qubit_data=qubit_data,
                                         qubit_count=qubit_count)
            >>> code = StabSubSystemCode(shell)
        """

        if gauge_group is None and shell is None:
            raise QiskitError("A Gauge group or Shell and associated data is required")

        if gauge_group is not None:
            self.gauge_group = gauge_group
            self._n = self.gauge_group.generators.num_qubits  # pylint: disable=invalid-name

        if gauge_group is not None and shell is not None:
            raise QiskitError("Only one of gauge_group or shell should be provided")

        if shell is not None:
            if qubit_data is None or qubit_count is None:
                raise QiskitError("Both qubit_data and qubit_count inputs are needed for a shell")

            self.shell = shell
            self.qubit_data = qubit_data
            self.qubit_count = qubit_count

            generators, self.qubit_data = Shell.shell2symplectic(
                self.shell, self.qubit_data, self.qubit_count
            )

            gauge_group = GaugeGroup(generators)
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

    @property
    def generators(self) -> PauliList:
        """Returns the Pauli generators for the subsystem"""
        return self.gauge_group.generators

    def draw(self, **kwargs) -> None:
        """Draws the subsytem code if a shell exists"""
        if self.shell is not None:
            self.shell.draw(qubit_data=self.qubit_data, qubit_count=self.qubit_count, **kwargs)
