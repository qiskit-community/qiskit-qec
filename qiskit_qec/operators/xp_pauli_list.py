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
#
# This code is based on the paper: "The XP Stabiliser Formalism: a
# Generalisation of the Pauli Stabiliser Formalism with Arbitrary Phases", Mark
# A. Webster, Benjamin J. Brown, and Stephen D. Bartlett. Quantum 6, 815
# (2022).
"""Module for XPPauli List"""
import numbers
from typing import Iterable, List, Tuple, Union, Optional

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit.quantum_info.operators.mixins import GroupMixin, LinearMixin
from qiskit_qec.operators.base_xp_pauli import BaseXPPauli
from qiskit_qec.utils import xp_pauli_rep


# pylint: disable=unused-argument
# pylint: disable=no-member
class XPPauliList(BaseXPPauli, LinearMixin, GroupMixin):
    """`XPPauliList` inherits from `BaseXPPauli`"""

    # Set the max number of qubits * xppaulis before string truncation
    _truncate__ = 2000

    def __init__(
        self,
        data: Union[BaseXPPauli, np.ndarray, Tuple[np.ndarray], Iterable, None] = None,
        phase_exp: Union[int, np.ndarray, None] = None,
        precision: Union[int, np.ndarray] = None,
        *,
        input_pauli_encoding: str = BaseXPPauli.EXTERNAL_XP_PAULI_ENCODING,
        input_qubit_order: str = "right-to-left",
        tuple_order: str = "zx",
    ) -> None:
        """Inits a XPPauliList

        Args:
            data (str): List of XPPauli Operators.
            phase_exp (int, optional): i**phase_exp. Defaults to 0.
            input_qubit_order (str, optional): Order to read pdata. Defaults to "right-to-left".
            precision: Precision of XP operators. Must be an integer/array of
            integers greater than or equal to 2.

        Raises:
            QiskitError: Something went wrong.
        """
        if data is None:
            matrix = np.empty(shape=(0, 0), dtype=np.bool_)
            phase_exp = np.empty(shape=(0,), dtype=np.int8)
        elif isinstance(data, BaseXPPauli):
            matrix = data.matrix
            phase_exp = data._phase_exp
            precision = data.precision
        # TODO elif isinstance(data, StabilizerTable), elif isinstance(data, PauliTable)
        elif isinstance(data, np.ndarray):
            if data.size == 0:
                matrix = np.empty(shape=(0, 0), dtype=np.bool_)
                phase_exp = np.empty(shape=(0,), dtype=np.int8)
            # TODO elif isinstance(data[0], str):
            else:
                if phase_exp is None:
                    phase_exp = 0
                matrix, phase_exp, precision = xp_pauli_rep.from_array(
                    data, phase_exp, precision, input_pauli_encoding=input_pauli_encoding
                )
        # TODO elif isinstance(data, tuple)
        else:
            # TODO Conversion as iterable of Paulis
            pass

        super().__init__(matrix, phase_exp, precision)

        # TODO
        # self.paulis = [
        #     Pauli(self.matrix[i], phase_exp=self._phase_exp[i]) for i in range(self.matrix.shape[0])
        # ]

    # ---------------------------------------------------------------------
    # Init Methods
    # ---------------------------------------------------------------------

    # ---------------------------------------------------------------------
    # Property Methods
    # ---------------------------------------------------------------------

    @property
    def phase(self):
        """Return the phase vector of the XPPauliList.

        Note: This is different from the quantum_info phase property which
        instead returns the phase_exp
        """
        # TODO
        pass

    @phase.setter
    def phase(self, phase):
        """Set the phase vector of the XPPauliList

        Args:
            phase (numpy.ndarray or complex numbers): Array of phases,
                phases must be one of [1,-1, 1j, -1j]
        """
        # TODO
        pass

    @property
    def shape(self):
        """The full shape of the :meth:`array`"""
        return self._num_xppaulis, self.num_qubits

    @property
    def size(self):
        """The number of XPPauli rows in the table."""
        return self._num_xppaulis

    @property
    def num_xppaulis(self) -> int:
        """Returns the number of XPPauli's in List"""
        return self._num_xppaulis

    @property
    def phase_exp(self):
        """Return the phase exponent vector of the XPPauliList"""
        # TODO
        pass

    @phase_exp.setter
    def phase_exp(self, phase_exp, input_phase_encoding=BaseXPPauli.EXTERNAL_PHASE_ENCODING):
        """Set the phase exponent vector of the XPPauliList. Note that this method
        converts the phase exponents directly and does not take into account the
        number of Y paulis in the representation.

        Args:
            phase_exp (_type_): _description_
            input_phase_encoding (_type_, optional): _description_. Defaults to
                BaseXPPauli.EXTERNAL_PHASE_ENCODING.
        """
        # TODO
        pass

    @property
    def settings(self):
        """Return settings."""
        return {"data": self.to_labels()}

    # ---------------------------------------------------------------------
    # Magic Methods and related methods
    # ---------------------------------------------------------------------

    def __getitem__(self, index):
        """Return a view of the XPPauliList."""
        # Returns a view of specified rows of the XPPauliList
        # This supports all slicing operations the underlying array supports.
        # TODO
        pass

    def getaslist(self, slc: Union[numbers.Integral, slice]) -> List["XPPauli"]:
        """_summary_

        Returns:
            _type_: _description_
        """
        # TODO
        pass

    def __setitem__(self, index, value):
        """Update XPPauliList."""
        # TODO
        pass

    def __repr__(self):
        """Display representation."""
        return self._truncated_str(True)

    def __str__(self):
        """Print representation."""
        return self._truncated_str(False)

    def _truncated_str(self, show_class):
        # TODO
        pass

    def __array__(self, dtype=None):
        """Convert to numpy array"""
        # pylint: disable=unused-argument
        shape = (len(self),) + 2 * (2**self.num_qubits,)
        ret = np.zeros(shape, dtype=complex)
        for i, mat in enumerate(self.matrix_iter()):
            ret[i] = mat
        return ret

    def __eq__(self, other):
        """Entrywise comparison of XPPauli equality."""
        if not isinstance(other, XPPauliList):
            other = XPPauliList(other)
        if not isinstance(other, BaseXPPauli):
            return False
        return self._eq(other)

    def __len__(self):
        """Return the number of XPPauli rows in the table."""
        return self._num_xppaulis

    def _add(self, other, qargs=None):
        """summary"""
        pass

    # ----
    #
    # ----

    # ---------------------------------------------------------------------
    # BaseOperator methods
    # ---------------------------------------------------------------------

    def tensor(self, other):
        """Return the tensor product with each XPPauli in the list.

        Args:
            other (XPPauliList): another XPPauliList.

        Returns:
            XPPauliList: the list of tensor product XPPaulis.

        Raises:
            QiskitError: if other cannot be converted to a XPPauliList, does
                         not have either 1 or the same number of XPPaulis as
                         the current list.
        """
        # TODO
        pass

    def compose(
        self,
        other: "XPPauliList",
        qargs: Optional[list] = None,
        front: bool = False,
        inplace: bool = False,
    ) -> "XPPauliList":
        """Return the composition self∘other for each XPPauli in the list.

        Note:
            This method is adapted from the method XPMul from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            other (XPPauliList): another XPPauliList.
            qargs (None or list): qubits to apply dot product on (Default: None).
            front (bool): If True use `dot` composition method [default: False].
            inplace (bool): If True update in-place (default: False).

        Returns:
            XPPauliList: the list of composed XPPaulis.

        Raises:
            QiskitError: if other cannot be converted to a XPPauliList, does
                         not have either 1 or the same number of XPPaulis as
                         the current list, has the wrong number of qubits
                         for the specified qargs, or if precision of other
                         does not match precision of self.

        Examples:
            >>> a = XPPauliList(
            ... data=np.array([[0, 1, 0, 0, 2, 0], [0, 1, 0, 0, 2, 0]], dtype=np.int64),
            ... phase_exp=np.array([6, 6]), precision=4)
            >>> b = XPPauliList(
            ... data=np.array([[1, 1, 1, 3, 3, 0], [1, 1, 1, 3, 3, 0]], dtype=np.int64),
            ... phase_exp=np.array([2, 2]), precision=4)
            >>> value = a.compose(b)
            >>> value.matrix
            array([[1, 0, 1, 3, 3, 0], [1, 0, 1, 3, 3, 0]], dtype=int64)
            >>> value._phase_exp
            array([6, 6])

        See also:
            _compose
        """
        if qargs is None:
            qargs = getattr(other, "qargs", None)
        if not isinstance(other, XPPauliList):
            other = XPPauliList(other)
        if len(other) not in [1, len(self)]:
            raise QiskitError(
                "Incompatible XPPauliLists. Other list must "
                "have either 1 or the same number of XPPaulis."
            )
        if self.precision != other.precision:
            raise QiskitError(
                "Precision of the two XPPauliLists to be multiplied must be the same."
            )

        return XPPauliList(super().compose(other, qargs=qargs, front=front, inplace=inplace))

    def rescale_precision(self, new_precision: int) -> "XPPauliList":
        """Rescale the generalized symplectic vector components
        of XPPauli operator to the new precision. Returns the rescaled XPPauli object.

        Note:
            This method is adapted from method XPSetN from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            new_precision: The target precision in which XPPauli is to be expressed

        Returns:
            XPPauliList: Resultant of rescaling the precision of XPPauliList

        Raises:
            QiskitError: If it is not possible to express XPPauliList in new_precision

        Examples:
            >>> matrix1 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0], dtype=np.int64)
            >>> phase_exp1 = 12
            >>> matrix2 = np.array([1, 1, 0, 1, 0, 0, 0, 0, 0, 4, 0, 0, 0, 6], dtype=np.int64)
            >>> phase_exp2 = 8
            >>> precision = 8
            >>> new_precision = 4
            >>> matrix = np.array([matrix1, matrix2])
            >>> phase_exp = np.array([phase_exp1, phase_exp2])
            >>> xppauli_list = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
            >>> rescaled_xppauli_list = xppaulilist.rescale_precision(new_precision=new_precision)
            >>> rescaled_xppauli_list.matrix
            array([[1, 1, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],
                [1, 1, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3]] dtype=np.int64)
            >>> rescaled_xppauli_list._phase_exp
            array([6, 4])

        See also:
            _rescale_precision
        """
        return XPPauliList(super().rescale_precision(new_precision))

    def antisymmetric_op(self, int_vec: np.ndarray) -> "XPPauliList":
        """Return the antisymmetric operators corresponding to the list of
        integer vectors, with precision specified by BaseXPPauli.

        Note:
            This method is adapted from method XPD from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            int_vec (np.ndarray): Array containing integer vectors

        Returns:
            XPPauliList: The antisymmetric operators corresponding to the input vectors

        Examples:
            >>> matrix = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 3, 3, 3],
            ... [0, 0, 0, 0, 0, 0, 0, 3, 1, 2, 3, 7, 6, 3]], dtype=np.int64)
            >>> phase_exp = np.array([0, 0])
            >>> precision = 8
            >>> xppauli_list = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
            >>> value = xppauli_list.antisymmetric_op(xppauli_list.z)
            >>> value.matrix
            array([[0, 0, 0, 0, 0, 0, 0, 0, -1, -2, -3, -3, -3, -3],
                [0, 0, 0, 0, 0, 0, 0, -3, -1, -2, -3, -7, -6, -3]], dtype=np.int64)
            >>> value._phase_exp
            array([15, 25])

        See also:
            _antisymmetric_op
        """
        return XPPauliList(super().antisymmetric_op(int_vec))

    def inverse(self) -> "XPPauli":
        """Return the inverse of the list of XP operators.

        Note:
            This method is adapted from method XPInverse from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            XPPauliList: Inverse of XPPauliList
        
        Examples:
            >>> matrix = np.array([[1, 1, 0, 1, 1, 0, 1, 2, 4, 4, 3, 1, 6, 1],
            ... [0, 1, 0, 0, 1, 0, 1, 7, 7, 3, 4, 6, 2, 7]], dtype=np.int64)
            >>> phase_exp = np.array([1, 0])
            >>> precision = 8
            >>> xppauli_list = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
            >>> value = xppauli_list.inverse()
            >>> value.matrix
            array([[1, 1, 0, 1, 1, 0, 1, 2, 4, 4, 3, 1, 2, 1],
                [0, 1, 0, 0, 1, 0, 1, 1, 7, 5, 4, 6, 6, 7]], dtype=np.int64)
            >>> value._phase_exp
            np.array([9, 8])
        """
        return XPPauliList(super().inverse())

    def power(self, n: int) -> "XPPauliList":
        """Return the XP operators of specified precision raised to the power n.

        Note:
            This method is adapted from method XPPower from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            n: The power to which XPPauliList is to be raised

        Returns:
            XPPauliList: XPPauliList raised to the power n

        Examples:
            >>> matrix = np.array([[1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 1],
            ... [1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 1]], dtype=np.int64)
            >>> phase_exp = np.array([12, 12])
            >>> precision = 8
            >>> n = 5
            >>> xppauli_list = XPPauliList(data=matrix, phase_exp=phase_exp, precision=precision)
            >>> value = xppauli_list.power(n=n)
            >>> value.matrix
            array([[1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 5],
                [1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 5]], dtype=np.int64)
            >>> value._phase_exp
            np.array([8, 8])

        See also:
            _power
        """
        return XPPauliList(super().power(n))

    def conjugate(self, other: "XPPauliList", front: bool = True, inplace: bool = False) -> "XPPauliList":
        """Return the conjugation of two XP operators.

        For single XP operators, this means

        A.conjugate(B, front=True) = A . B . A^{-1},

        where . is the XP Pauli multiplication and A^{-1} is the inverse of A.

        Likewise, 

        A.conjugate(B, front=False) = B . A . B^{-1}.

        For a list of XP operators, conjugation is performed element-wise:

        [A_1, ..., A_k].conjugate([B_1, ..., B_k]) = [A_1.conjugate(B_1), ..., A_k.conjugate(B_k)].

        TODO: This method currently only supports conjugation of two XP operator
        lists of the same length.

        Note:
            This method is adapted from method XPConjugate from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            other: List of XP operators to be conjugated with self
            front (bool, optional): Whether to conjugate in front (True) or
            behind (False), defaults to True
            inplace (bool, optional): Whether to perform the conjugation in
            place (True) or to return a new BaseXPPauli (False), defaults to
            False

        Returns:
            XPPauliList: List of conjugated XP operator

        Raises:
            QiskitError: Other list must have either 1 or the same number of
            XPPaulis

        Examples:
            >>> a = XPPauliList(data=np.array([[1, 0, 1, 1, 5, 3, 5, 4], 
            ... [1, 0, 1, 0, 1, 5, 2, 0]], dtype=np.int64),
            ... phase_exp=np.array([4, 7]), precision=6)
            >>> b = XPPauliList(data=np.array([[1, 0, 0, 1, 4, 1, 0, 1], 
            ... [0, 1, 1, 0, 1, 3, 0, 5]], dtype=np.int64),
            ... phase_exp=np.array([11, 2]), precision=6)
            >>> value = a.conjugate(b)
            >>> value.matrix
            array([[1, 0, 0, 1, 0, 1, 0, 1], 
            [0, 1, 1, 0, 5, 5, 4, 5]], dtype=int64)
            >>> value._phase_exp
            array([3, 10])
        """
        if not isinstance(other, XPPauliList):
            other = XPPauliList(other)
        if len(other) not in [1, len(self)]:
            raise QiskitError(
                "Incompatible XPPauliLists. Other list must "
                "have either 1 or the same number of XPPaulis."
            )

        return XPPauliList(super().conjugate(other, front=front, inplace=inplace))

    def commutator(self, other: "XPPauliList", front: bool = True, inplace: bool = False) -> "XPPauliList":
        """Return the commutator of two XP operators.

        For single XP operators, this means

        A.commutator(B, front=True) = [A, B] = A . B . A^{-1} . B^{-1},

        where . is the XP Pauli multiplication and A^{-1} is the inverse of A.

        Likewise, 

        A.commutator(B, front=False) = [B, A] = B . A . B^{-1}. A^{-1}.

        For a list of XP operators, commutator is computed element-wise:

        [A_1, ..., A_k].commutator([B_1, ..., B_k]) = [A_1.commutator(B_1), ..., A_k.commutator(B_k)].

        TODO: This method currently only supports commutator of two XP operator
        lists of the same length.

        Note:
            This method is adapted from method XPCommutator from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            other: List of XP operators to be in the commutator with self
            front (bool, optional): Whether self is the first element in the
            commutator (True) or second (False), defaults to True
            inplace (bool, optional): Whether to compute the commutator in
            place (True) or to return a new BaseXPPauli (False), defaults to
            False

        Returns:
            XPPauliList: List of commutators of XP operators

        Raises:
            QiskitError: Other list must have either 1 or the same number of
            XPPaulis

        Examples:
            >>> a = XPPauliList(data=np.array([[1, 0, 1, 1, 5, 3, 5, 4], 
            ... [1, 0, 1, 0, 1, 5, 2, 0]], dtype=np.int64),
            ... phase_exp=np.array([4, 7]), precision=6)
            >>> b = XPPauliList(data=np.array([[1, 0, 0, 1, 4, 1, 0, 1], 
            ... [0, 1, 1, 0, 1, 3, 0, 5]], dtype=np.int64),
            ... phase_exp=np.array([11, 2]), precision=6)
            >>> value = a.commutator(b)
            >>> value.matrix
            array([[0, 0, 0, 0, 4, 0, 0, 0], 
            [0, 0, 0, 0, 4, 4, 2, 0]], dtype=int64)
            >>> value._phase_exp
            array([8, 8])
        """
        if not isinstance(other, XPPauliList):
            other = XPPauliList(other)
        if len(other) not in [1, len(self)]:
            raise QiskitError(
                "Incompatible XPPauliLists. Other list must "
                "have either 1 or the same number of XPPaulis."
            )

        return XPPauliList(super().commutator(other, front=front, inplace=inplace))

    # def conjugate(self):
    #     """Return the conjugate of each XPPauli in the list."""
    #     # TODO
    #     pass

    # def transpose(self):
    #     """Return the transpose of each XPPauli in the list."""
    #     # TODO
    #     pass

    # def adjoint(self):
    #     """Return the adjoint of each XPPauli in the list."""
    #     # TODO
    #     pass

    # def inverse(self):
    #     """Return the inverse of each XPPauli in the list."""
    #     # TODO
    #     pass

    # ---------------------------------------------------------------------
    # Utility methods
    # ---------------------------------------------------------------------

    # ---------------------------------------------------------------------
    # Custom Iterators
    # ---------------------------------------------------------------------

    # ---------------------------------------------------------------------
    # Class methods
    # ---------------------------------------------------------------------

    @classmethod
    def from_symplectic(cls, z, x, phase_exp=0):
        """Construct a XPPauliList from a symplectic data.

        Args:
            z (np.ndarray): 2D boolean Numpy array.
            x (np.ndarray): 2D boolean Numpy array.
            phase_exp (np.ndarray or None): Optional, 1D integer array from Z_4.

        Returns:
            XPPauliList: the constructed XPPauliList.

        Note: Initialization this way will copy matrices and not reference them.

        TODO: Fix this method to be more general and not in old form only
            (i.e. include matrix inputs ...)
        """
        # TODO
        pass
