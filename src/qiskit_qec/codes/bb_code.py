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
"""Define the bivariate bicylce code."""

from typing import Optional, Sequence

import numpy as np

from qiskit_qec.linear.symplectic import normalizer
import qiskit_qec.utils.pauli_rep as pauli_rep
from qiskit_qec.codes.codebuilders.builder import Builder
from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode
from qiskit_qec.operators.base_pauli import BasePauli
from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.structures.gauge import GaugeGroup


class BBCode:
    """Bivariate Bicycle code data.

    The X and Z gauge operator lists are given as lists of supports.
    There is a consistent qubit ordering of these lists so that
    we can construct gate schedules for circuits.
    """

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

        self._hx = None
        self._hz = None
        self._x_stabilizers = None
        self._z_stabilizers = None
        self._logical_z = None
        self._logical_x = None

    @property
    def hx(self) -> np.array:
        if not self.is_defined():
            raise AttributeError(f'A or B undefined, first set them via set_A and set_B')
        if self._hx is None: self._hx = np.hstack([self.A, self.B])
        return self._hx

    @property
    def hz(self) -> np.array:
        if not self.is_defined():
            raise AttributeError(f'A or B undefined, first set them via set_A and set_B')
        if self._hz is None: self._hz = np.hstack([self.B.T, self.A.T])
        return self._hz
    
    @property
    def x_stabilizers(self) -> PauliList:
        if self._x_stabilizers is None: self._x_stabilizers = PauliList(np.hstack([self.hx, np.zeros((self.n//2, self.n))]))
        return self._x_stabilizers
    
    @property
    def z_stabilizers(self) -> PauliList:
        if self._z_stabilizers is None: self._z_stabilizers = PauliList(np.hstack([np.zeros((self.n//2, self.n)), self.hz]))
        return self._z_stabilizers
    
    @property
    def x_gauges(self) -> PauliList:
        return self.x_stabilizers
    
    @property
    def z_gauges(self) -> PauliList:
        return self.z_stabilizers
    
    @property
    def logical_z(self):
        if self._logical_z is None:
            full_sym = np.vstack([np.hstack([self.hx, np.zeros((self.n//2, self.n))]),
                                  np.hstack([np.zeros((self.n//2, self.n)), self.hz])])
            center_, x_new, z_new = normalizer(full_sym)
            self._logical_z = PauliList(z_new)
            self._logical_x = PauliList(x_new)

        #return self._logical_z.matrix.astype(int).tolist()
        return self._logical_z
    
    @property
    def logical_x(self):
        if self._logical_x is None:
            full_sym = np.vstack([np.hstack([self.hx, np.zeros((self.n//2, self.n))]),
                                  np.hstack([np.zeros((self.n//2, self.n)), self.hz])])
            center_, x_new, z_new = normalizer(full_sym.astype(np.bool_))
            self._logical_z = PauliList(z_new)
            self._logical_x = PauliList(x_new)

        return self._logical_x
    
    @property
    def x_boundary(self):
        raise NotImplementedError()
    
    @property
    def z_boundary(self):
        raise NotImplementedError()
    
    @property
    def k(self):
        raise NotImplementedError()

    @property
    def d(self):
        raise NotImplementedError()

    def __str__(self) -> str:
        """Formatted string."""
        return f"(l={self.l},m={self.m}) bivariate bicycle code"
        #return f"[[{self.n}, {self.k}, {self.d}]] heavy-hexagon compass code"

    def __repr__(self) -> str:
        """String representation."""
        val = str(self)
        val += f"\nx_gauges = {self.x_gauges}"
        val += f"\nz_gauges = {self.z_gauges}"
        val += f"\nx_stabilizers = {self.x_stabilizers}"
        val += f"\nz_stabilizers = {self.z_stabilizers}"
        val += f"\nlogical_x = {self.logical_x}"
        val += f"\nlogical_z = {self.logical_z}"
        val += f"\nx_boundary = {self.x_boundary}"
        val += f"\nz_boundary = {self.z_boundary}"
        return val
    
    def is_defined(self):
        return self.A is not None and self.B is not None
    
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

        return np.kron(BBCode.cs_pow(self.l, power=power), np.eye(self.m)).astype(np.uint8)
    
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

        return np.kron(np.eye(self.l), BBCode.cs_pow(self.m, power=power)).astype(np.uint8)
    
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
    
    @staticmethod
    def symplectic_to_indices(arr):
        """ Returns """