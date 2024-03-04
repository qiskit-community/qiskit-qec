# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2022
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Subsystem Surface code builder example"""

from typing import Optional, Sequence

import numpy as np

import qiskit_qec.utils.pauli_rep as pauli_rep
from qiskit_qec.codes.codebuilders.builder import Builder
from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode
from qiskit_qec.operators.base_pauli import BasePauli
from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.structures.gauge import GaugeGroup


class BivariateBicycleCodeBuilder(Builder):
    """Bivariate Bicycle Code (arXiv:2308.07915) Builder Class"""

    # pylint: disable=anomalous-backslash-in-string
    def __init__(
        self,
        l: int,
        m: int,
        p1: Optional[Sequence[int]] = None,
        p2: Optional[Sequence[int]] = None
    ) -> None:
        """Initializes a bivariate bicycle code builder

        If p1 is not specified, then A has to be specified manually by set_A.
        If p2 is not specified, then B has to be specified manually by set_B.

        Args:
            l: dimension of first space of generators x and y,
            m: dimension of second space of generators x and y,
            p1 (optional): p1 = (a,b,c) => A = x^a + y^b + y^c. Defaults to None.
            p2 (optional): p2 = (d,e,f) => B = y^d + x^e + x^f. Defaults to None.

        Examples:
            Example 1:
            >>> code = BivariateBicycleCodeBuilder(l=6, m=6, p1=(3,1,2), p2=(3,1,2)).build()

            Example 2:
            >>> builder = BivariateBicycleCodeBuilder(l=6, m=6)
            >>> builder.set_A(builder.gen_x(3), builder.gen_y(1), builder.gen_y(2))
            >>> builder.set_B(builder.gen_y(3), builder.gen_x(1), builder.gen_x(2))
            >>> code = builder.build()

            Examples 1 and 2 produce the same code. The example 1 syntax is easier when creating codes in the standard
            xyy, yxx form. The example 2 syntax allows for more flexible code building with non-standard forms.

        """

        self.l = l
        self.m = m
        self.n = 2*l*m

        if p1 is not None:
            a, b, c = p1
            self.set_A(self.gen_x(a), self.gen_y(b), self.gen_y(c))
        else: self.A = None

        if p2 is not None:
            d, e, f = p2
            self.set_B(self.gen_y(d), self.gen_x(e), self.gen_x(f))
        else: self.B = None

    def build(self) -> StabSubSystemCode:
        """
        Buils the Stabilizer surface code and returns it.

        Returns
        -------
        StabSubSystemCode
            The Stabilizer code defined by the parity check matrices
            hx = [A|B] and hz = [B^T, A^T].

        Raises
        ------
        AttributeError
            When A or B are None (have not been defined yet).
        
        """

        gauge_group = self._create_gauge_group()
        return StabSubSystemCode(gauge_group=gauge_group)
    
    def set_A(self, A1, A2, A3) -> None:
        """
        Set the matrix A as the sum of A1, A2 and A3 modulo 2.

        Parameters
        ----------
        A1 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.
        A2 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.
        A3 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.

        Examples
        --------

        >>> builder.set_A(builder.gen_x(3), builder.gen_y(1), builder.gen_y(2))
        Sets A = x^3 + y + y^2
        """

        self.A = ((A1 + A2 + A3) % 2).astype(np.uint8)

    def set_B(self, B1, B2, B3) -> None:
        """
        Set the matrix B as the sum of B1, B2 and B3 modulo 2.

        Parameters
        ----------
        B1 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.
        B2 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.
        B3 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.

        Examples
        --------

        >>> builder.set_B(builder.gen_y(3), builder.gen_x(1), builder.gen_x(2))
        Sets B = y^3 + x + x^2
        """

        self.B = ((B1 + B2 + B3) % 2).astype(np.uint8)

    def gen_x(self, power: int = 1):
        """
        Calculates and returns x^power

        Parameters
        ----------
        power : int
            Power to which x is raised. Defaults to 1.

        Returns
        --------
        np.array
            x^power

        Examples
        --------

        >>> builder.gen_x(3) # x^3
        >>> builder.gen_x(1) # x
        >>> builder.gen_x() # x
        """

        return np.kron(BivariateBicycleCodeBuilder.cs_pow(self.l, power=power), np.eye(self.m)).astype(np.uint8)
    
    def gen_y(self, power=1):
        """
        Calculates and returns y^power

        Parameters
        ----------
        power : int
            Power to which y is raised. Defaults to 1.

        Returns
        --------
        np.array
            y^power

        Examples
        --------

        >>> builder.gen_y(3) # y^3
        >>> builder.gen_y(1) # y
        >>> builder.gen_y() # y
        """

        return np.kron(np.eye(self.l), BivariateBicycleCodeBuilder.cs_pow(self.m, power=power)).astype(np.uint8)
    
    @staticmethod
    def cs_pow(l, power=1):
        """
        Calculates and returns a power of the cyclic shift matrix of size lxl (C_l)

        Parameters
        ----------
        l     : int
            size of cyclic shift matrix.
        power : int
            Power to which Cl is raised. Defaults to 1.

        Returns
        --------
        np.array
            C_l^power

        Examples
        --------

        >>> BivariateBicycleCodeBuilder.cs_pow(3, 2) # (C_3)^2
        """

        return np.roll(np.eye(l), shift=power, axis=1).astype(np.uint8)

    def _create_gauge_group(self):
        """
        Creates the gauge group of this code.

        Returns
        -------
        GaugeGroup
            The Gauge Group (Stabilizer Group) defined by the parity check matrices
            hx = [A|B] and hz = [B^T, A^T].

        Raises
        ------
        AttributeError
            When A or B are None (have not been defined yet).
        
        """
        if self.A is None or self.B is None:
            raise AttributeError(f'A or B is undefined, first set them via set_A and set_B')
        
        # create parity check matrices, padded with zeros (for conversion to PauliList from sympletic representation)
        hx_s = np.hstack([self.A, self.B, np.zeros((self.n//2, self.n))])
        hz_s = np.hstack([np.zeros((self.n//2, self.n)), self.B.T, self.A.T])

        s = np.vstack([hx_s, hz_s]) # total sympletic representation of all stabilizers of this code

        generators = PauliList(s)
        gauge_group = GaugeGroup(generators)
        return gauge_group


def main():
    """ test """
    code1 = BivariateBicycleCodeBuilder(l=6, m=6, p1=(3,1,2), p2=(3,1,2)).build()

    builder = BivariateBicycleCodeBuilder(l=6, m=6)
    builder.set_A(builder.gen_x(3), builder.gen_y(1), builder.gen_y(2))
    builder.set_B(builder.gen_y(3), builder.gen_x(1), builder.gen_x(2))
    code2 = builder.build()

    print(code1._n)
    BasePauli.set_syntax(pauli_rep.INDEX_SYNTAX)
    print(code1.gauge_group.generators)
    print(code1.gauge_group.generators[0])
    BasePauli.set_syntax(pauli_rep.PRODUCT_SYNTAX)
    print(code1.gauge_group.generators[0])
    BasePauli.set_syntax(pauli_rep.LATEX_SYNTAX)
    print(code1.gauge_group.generators[0])
    print(code1.gauge_group.generators==code2.gauge_group.generators)

if __name__ == "__main__":
    main()