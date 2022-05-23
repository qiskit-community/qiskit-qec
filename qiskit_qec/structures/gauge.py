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

from typing import Union, List, Any
import warnings

import numpy as np
from qiskit import QiskitError

import qiskit_qec.linear.symplectic as sym
from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.operators.base_pauli import BasePauli
from qiskit_qec.structures.group import Group


class GaugeGroup(Group):
    """`GaugeGroup` inherits from `Group`"""

    def __init__(
        self,
        generators: Union[str, List[str], BasePauli, None] = None,
        name: Any = None,
        *,
        isotropic_generators: Union[str, List[str], np.ndarray, None] = None,
        hyperbolic_generators: Union[str, List[str], np.ndarray, None] = None,
        with_generators=True,
    ) -> None:
        """_summary_

        Args:
            generators (Union[str,List[str], BasePauli]): _description_
            name (Any, optional): _description_. Defaults to None.
            isotropic_generators (Union[str,List[str],np.ndarray,None], optional): _description_. Defaults to None.
            hyperbolic_generators (Union[str,List[str],np.ndarray,None], optional): _description_. Defaults to None.
            with_generators (bool, optional): _description_. Defaults to True.
        """

        if generators is None and isotropic_generators is None and hyperbolic_generators is None:
            raise QiskitError("Some generators are required to create GaugeGroup")

        generators = (
            PauliList(generators)
            + PauliList(isotropic_generators)
            + PauliList(hyperbolic_generators)
        )

        if with_generators:
            generators = PauliList(sym.min_generating(generators.matrix))

        self.generators = generators
        self._with_generators = with_generators
        if name is None:
            self.name = ""
        else:
            self.name = name

        super().__init__()

    @property
    def n(self):
        """Returns the number of qubits"""
        return self.generators.num_qubits

    @property
    def k(self):
        """Returns the number of logical operators"""
        warnings.warn("Number of logical operators - Property not yet implemented")
        return None

    @property
    def num_gen(self):
        """Returns the number of generators"""
        warnings.warn("Number of generators - Property not yet implemented")
        return None

    def __str__(self) -> str:
        return f"GaugeGroup[ [[{self.n},{self.k}]] with {self.num_gen} min generators]"


class GaugeGroupWithGenerators:
    """Returns a GaugeGroup instance generated by generators. object.generators
    does not necessarily return the input generators"""

    def __new__(cls, generators: Union[str, List[str], BasePauli], name: Any = None) -> None:
        """_summary_

        Args:
            generators (Union[str,List[str], BasePauli]): _description_
            name (Any, optional): _description_. Defaults to None.
        """
        return GaugeGroup(generators, name=name, with_generators=True)


class GaugeGroupByGenerators:
    """Returns a GaugeGroup instance generated by generators. object.generators
    will return the input generators"""

    def __new__(cls, generators: Union[str, List[str], BasePauli], name: Any = None) -> None:
        """_summary_

        Args:
            generators (Union[str,List[str], BasePauli]): _description_
            name (Any, optional): _description_. Defaults to None.
        """
        return GaugeGroup(generators, name=name, with_generators=False)
