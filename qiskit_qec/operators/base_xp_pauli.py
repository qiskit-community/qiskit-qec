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
"""Module for base XP pauli"""

import numbers
from typing import List, Optional, Union

import numpy as np

# Must be imported as follows to avoid circular import errors
import qiskit.quantum_info.operators.symplectic.clifford
from qiskit import QiskitError
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.barrier import Barrier
from qiskit.circuit.delay import Delay
from qiskit.circuit.instruction import Instruction
from qiskit.quantum_info.operators.base_operator import BaseOperator
from qiskit.quantum_info.operators.mixins import AdjointMixin, MultiplyMixin
from qiskit_qec.utils import xp_pauli_rep


# pylint: disable=no-member
class BaseXPPauli(BaseOperator, AdjointMixin, MultiplyMixin):
    r"""Symplectic representation of a list of N-qubit XP operators with phases
    using numpy arrays for generalized symplectic matrices and phase vectors.

    Base class for XPPauli and XPPauliList.
    """

    # External string formats used when displaying Pauli's as strings
    EXTERNAL_TENSOR_ENCODING = xp_pauli_rep.DEFAULT_EXTERNAL_TENSOR_ENCODING
    EXTERNAL_PHASE_ENCODING = xp_pauli_rep.DEFAULT_EXTERNAL_PHASE_ENCODING
    EXTERNAL_XP_PAULI_ENCODING = EXTERNAL_PHASE_ENCODING + EXTERNAL_TENSOR_ENCODING

    EXTERNAL_SYNTAX = xp_pauli_rep.PRODUCT_SYNTAX
    EXTERNAL_QUBIT_ORDER = xp_pauli_rep.DEFAULT_QUBIT_ORDER

    PRINT_PHASE_ENCODING = None

    # pylint: disable=unused-argument
    def __init__(
        self,
        matrix: Union[np.ndarray, None] = None,
        phase_exp: Union[None, np.ndarray, np.integer] = None,
        precision: int = None,
        order: str = "xz",
    ) -> None:
        # TODO: Need to update this docstring once the representation to be
        # used for XP operators has been decided. In addition,
        # (-i)^phase_exp Z^z X^x and GF(2) need to be updated to XP operators.
        """A BaseXPPauli object represents a list N-qubit XP Pauli operators with phases.
        Numpy arrays are used to represent the generalized symplectic matrix represention of these
        XP operators. The phases of the XP operators are stored encoded. The phases of the XP operators
        operators are internally encoded in the '-iZX' Pauli encoding (See the xp_pauli_rep
        module for more details). That is a XP operator is represented as generalized symplectic
        vector V and a phase exponent phase_exp such that:

        (-i)^phase_exp Z^z X^x

        where V = [x, z] and phase_exp is a vector of Z_2N elements (0,1,2,...,2N-1). A list
        of XP operators is represented as a generalized symplectic matrix S and a phase exponent
        vector phase_exp such that the rows or S are the generalized symplectic vector representations
        of the XP operators and the phase_exp vector store the phase exponent of each
        associated XP Operator.

        Args:
            matrix: Input GF(2) symplectic matrix
            phase_exp (optional): Phase exponent vector for input matrix. A value of None will
                result in an a complex coefficients of 1 for each XP operator. Defaults to None.
            precision: Precision of XP operators. Must be an integer greater than or equal to 2.
            order: Set to 'xz' or 'zx'. Defines which side the x and z parts of the input matrix

        Raises:
            QiskitError: matrix and phase_exp sizes are not compatible, or if precision is less than 2.

        Examples:
            >>> matrix = numpy.array([[1,1,0,0],[0,1,0,1]])
            >>> base_xp_pauli = BaseXPPauli(matrix)

        See Also:
            XPPauli, XPPauliList
        """

        if not (isinstance(precision, int) and (precision > 1)):
            raise QiskitError(
                "Precision of XP operators must be an integer greater than or equal to 2."
            )

        if matrix is None or matrix.size == 0:
            matrix = np.empty(shape=(0, 0), dtype=np.int64)
            phase_exp = np.empty(shape=(0,), dtype=np.int64)
        matrix = np.atleast_2d(matrix)
        num_qubits = matrix.shape[1] >> 1

        if order == "zx":
            nmatrix = np.empty(shape=matrix.shape, dtype=matrix.dtype)
            nmatrix[:, :num_qubits] = matrix[:, num_qubits:]
            nmatrix[:, num_qubits:] = matrix[:, :num_qubits]
            matrix = nmatrix

        self.matrix = matrix
        self.precision = precision
        self._num_xppaulis = self.matrix.shape[0]
        if phase_exp is None:
            self._phase_exp = np.zeros(shape=(self.matrix.shape[0],), dtype=np.int64)
        else:
            self._phase_exp = np.atleast_1d(phase_exp)
        if self._phase_exp.shape[0] != self.matrix.shape[0]:
            raise QiskitError("matrix and phase_exp sizes are not compatible")

        super().__init__(num_qubits=num_qubits)

    # ---------------------------------------------------------------------
    # Properties
    # ---------------------------------------------------------------------

    @property
    def x(self):
        """_summary_"""
        return self.matrix[:, : self.num_qubits]

    @x.setter
    def x(self, val: np.ndarray):
        """_summary_"""
        self.matrix[:, : self.num_qubits] = val

    # @final Add when python >= 3.8
    @property
    def _x(self):  # pylint: disable=invalid-name
        """_summary_"""
        return self.matrix[:, : self.num_qubits]

    # @final Add when python >= 3.8
    @_x.setter
    def _x(self, val):  # pylint: disable=invalid-name
        """_summary_"""
        self.matrix[:, : self.num_qubits] = val

    @property
    def z(self):
        """_summary_"""
        return self.matrix[:, self.num_qubits :]

    @z.setter
    def z(self, val):
        """_summary_"""
        self.matrix[:, self.num_qubits :] = val

    # @final Add when python >= 3.8
    @property
    def _z(self):  # pylint: disable=invalid-name
        """_summary_"""
        return self.matrix[:, self.num_qubits :]

    # @final Add when python >= 3.8
    @_z.setter
    def _z(self, val):  # pylint: disable=invalid-name
        """_summary_"""
        self.matrix[:, self.num_qubits :] = val

    @property
    def num_y(self):
        """_summary_"""
        pass

    @property
    def tensor_encoding(self):
        """_summary_"""
        pass

    @classmethod
    def set_tensor_encoding(cls, encoding: str = xp_pauli_rep.DEFAULT_EXTERNAL_TENSOR_ENCODING):
        """_summary_"""
        pass

    @property
    def phase_encoding(self):
        """_summary_"""
        pass

    @classmethod
    def set_phase_encoding(cls, encoding: str = xp_pauli_rep.DEFAULT_EXTERNAL_PHASE_ENCODING):
        """_summary_"""
        pass

    @property
    def pauli_encoding(self):
        """_summary_"""
        pass

    @classmethod
    def set_pauli_encoding(cls, encoding: str = xp_pauli_rep.DEFAULT_EXTERNAL_XP_PAULI_ENCODING):
        """_summary_"""
        pass

    @property
    def syntax(self):
        """_summary_"""
        pass

    @classmethod
    def set_syntax(cls, syntax_code: Optional[int] = None, syntax_str: Optional[str] = "Product"):
        """_summary_"""
        pass

    @property
    def print_phase_encoding(self):
        """_summary_"""
        pass

    @classmethod
    def set_print_phase_encoding(cls, phase_encoding: Optional[str] = None):
        """_summary_"""
        pass

    @property
    def qubit_order(self):
        """_summary_"""
        pass

    @classmethod
    def set_qubit_order(cls, qubit_order: Optional[str] = None):
        """_summary_"""
        pass

    # ---------------------------------------------------------------------
    # Magic Methods
    # ---------------------------------------------------------------------

    def __imul__(self, other: "BaseXPPauli") -> "BaseXPPauli":
        """_summary_"""
        pass

    def __neg__(self) -> "BaseXPPauli":
        """_summary_"""
        pass

    # ---------------------------------------------------------------------

    def copy(self) -> "BaseXPPauli":
        """_summary_"""
        pass

    # ---------------------------------------------------------------------
    # BaseOperator methods
    # ---------------------------------------------------------------------
    # Needed by GroupMixin class from BaseOperator class

    # pylint: disable=arguments-differ
    def compose(
        self,
        other: "BaseXPPauli",
        qargs: Optional[list] = None,
        front: bool = False,
        inplace: bool = False,
    ) -> "BaseXPPauli":
        r"""Return the composition of XPPaulis lists

        To be consistent with other compose functions in Qiskit, composition is defined via
        left multiplication. That is

        A.compose(B) = B.A = B.dot(A) = A.compose(B, front=False)

        where . is the Pauli group multiplication and so B is applied after A. Likewise

        A.compose(B, front=True) = A.B = A.dot(B)

        That is B is applied first or at the front.

        This compose is:

        [A_1,A_2,...,A_k].compose([B_1,B_2,...,B_k]) = [A_1.compose(B_1),...,A_k.compose(B_k)]

        or

        [A].compose([B_1,B_2,...,B_k])) = [A.compose(B_1),...,A.compose(B_k)]

        Note:
            This method does compose coordinate wise (which is different from the PauliTable compose
            which should be corrected at some point).
            This method is adapted from method XPMul from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            other: BaseXPPauli
            front (bool): (default: False)
            qargs (list or None): Optional, qubits to apply compose on
                                  on (default: None->All).
            inplace (bool): If True update in-place (default: False).

        Returns:
            BaseXPPauli : Compositon of self and other

        Raises:
            QiskitError: if number of qubits of other does not match qargs, or
                         if precision of other does not match precision of self.

        Examples:
            >>> a = BaseXPPauli(matrix=np.array([0, 1, 0, 0, 2, 0], dtype=np.int64),
            ... phase_exp=6, precision=4)
            >>> b = BaseXPPauli(matrix=np.array([1, 1, 1, 3, 3, 0], dtype=np.int64),
            ... phase_exp=2, precision=4)
            >>> value = BaseXPPauli.compose(a, b)
            >>> value.matrix
            array([[1, 0, 1, 3, 3, 0]], dtype=int64)
            >>> value._phase_exp
            array([6])

        See also:
            _compose
        """

        # Validation
        if qargs is None and other.num_qubits != self.num_qubits:
            raise QiskitError(f"other {type(self).__name__} must be on the same number of qubits.")

        if qargs and other.num_qubits != len(qargs):
            raise QiskitError(
                f"Number of qubits of the other {type(self).__name__} does not match qargs."
            )

        if other._num_xppaulis not in [1, self._num_xppaulis]:
            raise QiskitError(
                "Incompatible BaseXPPaulis. Second list must "
                "either have 1 or the same number of XPPaulis."
            )

        if self.precision != other.precision:
            raise QiskitError(
                "Precision of the two BaseXPPaulis to be multiplied must be the same."
            )

        return self._compose(self, other, qargs=qargs, front=front, inplace=inplace)

    @staticmethod
    def _compose(
        a: "BaseXPPauli",
        b: "BaseXPPauli",
        qargs: Optional[list] = None,
        front: bool = False,
        inplace: bool = False,
    ) -> "BaseXPPauli":
        """Returns the composition of two BaseXPPauli objects.

        Note:
            This method is adapted from method XPMul from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            a : BaseXPPauli object
            b : BaseXPPauli object
            qargs (Optional[list], optional): _description_. Defaults to None.
            front (bool, optional): _description_. Defaults to False.
            inplace (bool, optional): _description_. Defaults to False.

        Returns:
            BaseXPPauli: _description_

        See also:
            unique_vector_rep, _unique_vector_rep
        """

        if qargs is not None:
            qargs = list(qargs) + [item + a.num_qubits for item in qargs]
            amat = a.matrix[:, qargs]
        else:
            amat = a.matrix
        bmat = b.matrix

        # Calculate the sum of generalized symplectic matrix for the composition, excluding D
        x = np.logical_xor(amat[:, : a.num_qubits], bmat[:, : b.num_qubits])
        z = amat[:, a.num_qubits :] + bmat[:, b.num_qubits :]
        mat = np.concatenate((x, z), axis=-1)

        # Calculate the phase of the composition, excluding D
        phase_exp = a._phase_exp + b._phase_exp
        # Calculate antisymmetric operator, i.e. D
        if front:
            dx = np.zeros(np.shape(a.x))
            dz = 2 * np.multiply(b.x, a.z)
            dmat = np.concatenate((dx, dz), axis=-1)
            d = BaseXPPauli(matrix=dmat, precision=a.precision)._antisymmetric_op()
        else:
            dx = np.zeros(np.shape(a.x))
            dz = 2 * np.multiply(a.x, b.z)
            dmat = np.concatenate((dx, dz), axis=-1)
            d = BaseXPPauli(matrix=dmat, precision=a.precision)._antisymmetric_op()

        if qargs is None:
            if not inplace:
                result_x = np.logical_xor(x, d.x)
                result_z = z + d.z
                result_phase_exp = phase_exp + d._phase_exp
                result_mat = np.concatenate((result_x, result_z), axis=-1)
                return BaseXPPauli(
                    matrix=result_mat, phase_exp=result_phase_exp, precision=a.precision
                )._unique_vector_rep()
            # Inplace update
            a.x = np.logical_xor(x, d.x)
            a.z = z + d.z
            a._phase_exp = phase_exp + d._phase_exp
            return a._unique_vector_rep()

        # Qargs update
        ret = a if inplace else a.copy()
        ret.matrix[:, qargs] = mat
        ret._phase_exp = phase_exp + d._phase_exp
        ret = ret._unique_vector_rep()
        return ret

    # ---------------------------------------------------------------------

    def tensor(self, other):
        """_summary_"""
        pass

    @staticmethod
    def _tensor(a: "BaseXPPauli", b: "BaseXPPauli"):
        """_summary_"""
        pass

    # ---------------------------------------------------------------------

    def expand(self, other):
        """_summary_"""
        pass

    # ---------------------------------------------------------------------
    # Needed by MultiplyMixin class
    # ---------------------------------------------------------------------

    def _multiply(self, phase, roundit=True) -> "BaseXPPauli":  # pylint: disable=arguments-renamed
        """_summary_"""
        pass

    # Needed by AdjointMixin class

    def conjugate(self, inplace=False) -> "BaseXPPauli":
        """_summary_"""
        pass

    def transpose(self, inplace: bool = False) -> "BaseXPPauli":
        """_summary_"""
        pass

    def commutes(self, other: "BaseXPPauli", qargs: List = None) -> np.ndarray:
        """_summary_"""
        pass

    def _commutes(self, other: "BaseXPPauli", qargs: List = None) -> np.ndarray:
        """_summary_"""
        pass

    # ---------------------------------------------------------------------
    # Extra all_commutes method
    # ---------------------------------------------------------------------

    def all_commutes(self, other: "BaseXPPauli") -> np.ndarray:
        """_summary_"""
        pass

    # ---------------------------------------------------------------------
    # Evolve Methods
    # ---------------------------------------------------------------------

    def evolve(self, other: "BaseXPPauli", qargs: Union[None, List, int] = None, frame: str = "h"):
        """_summary_"""
        pass

    def _evolve_clifford(
        self,
        other: qiskit.quantum_info.operators.symplectic.clifford.Clifford,
        qargs: Union[None, List, int] = None,
        frame: str = "h",
    ):
        """_summary_"""
        pass

    # ---------------------------------------------------------------------
    # Helper Metods
    # ---------------------------------------------------------------------

    def _eq(self, other):
        """_summary_"""
        pass

    # ---------------------------------------------------------------------
    # Class methods for creating labels -> Uses xp_pauli_rep suite of methods
    # ---------------------------------------------------------------------

    def to_label(
        self,
        output_pauli_encoding: Optional[str] = None,
        no_phase: bool = False,
        return_phase: bool = False,
        syntax: Optional[int] = None,
        qubit_order: Optional[str] = None,
        index_start: int = 0,
        squeeze: bool = True,
        index_str: str = "",
    ) -> Union[str, List[str]]:
        """_summary_"""
        pass

    # ---------------------------------------------------------------------
    # The methods below should be deprecated eventually as they simply refer
    # to the newer versions that are now elsewhere.
    # ---------------------------------------------------------------------

    def _count_y(self) -> np.ndarray:
        """_summary_"""
        pass

    @staticmethod
    def _stack(array: np.ndarray, size: int, vertical: bool = True) -> np.ndarray:
        """_summary_"""
        pass

    @staticmethod
    def _phase_from_complex(coeff: numbers.Complex) -> np.ndarray:
        """_summary_"""
        pass

    # ---------------------------------------------------------------------
    # Apply Clifford to BaseXPPauli
    # ---------------------------------------------------------------------

    def _append_circuit(
        self,
        circuit: Union[Barrier, Delay, QuantumCircuit, Instruction],
        qargs: Optional[List] = None,
    ) -> "BaseXPPauli":
        """_summary_"""
        pass

    # ---------------------------------------------------------------------
    # BaseXPPauli methods for XP arithmetic
    # ---------------------------------------------------------------------

    def unique_vector_rep(self) -> "BaseXPPauli":
        """Convert the BaseXPPauli operator into unique vector form, ie
        phase_exp in Z modulo 2*precision, x in Z_2, z in Z modulo
        precision.

        Note:
            This method is adapted from method XPRound from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            BaseXPPauli: Unique vector representation of self

        Examples:
            >>> a = BaseXPPauli(matrix=np.array([0, 3, 1, 6, 4, 3], dtype=np.int64),
            ... phase_exp=11, precision=4)
            >>> a = a.unique_vector_rep()
            >>> a.matrix
            np.array([[0, 1, 1, 2, 0, 3]], dtype=int64)
            >>> a._phase_exp
            array([3], dtype=int32)

        See also:
            _unique_vector_rep
        """
        return self._unique_vector_rep()

    def _unique_vector_rep(self) -> "BaseXPPauli":
        """Convert the BaseXPPauli operator into unique vector form, ie
        phase_exp in Z modulo 2*precision, x in Z_2, z in Z modulo
        precision.

        Note:
            This method is adapted from method XPRound from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            BaseXPPauli: Unique vector representation of self
        """
        matrix = np.empty(shape=np.shape(self.matrix), dtype=np.int64)

        phase_exp = np.mod(self._phase_exp, 2 * self.precision)
        matrix[:, : self.num_qubits] = np.mod(self.x, 2)
        matrix[:, self.num_qubits :] = np.mod(self.z, np.expand_dims(self.precision, axis=-1))

        return BaseXPPauli(matrix, phase_exp, self.precision)

    def rescale_precision(self, new_precision: int) -> "BaseXPPauli":
        """Rescale the generalized symplectic vector components of BaseXPPauli
        operator to the new precision.

        Note:
            This method is adapted from method XPSetN from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            new_precision: The target precision in which BaseXPPauli is to be expressed

        Returns:
            BaseXPPauli: Resultant of rescaling the precision of BaseXPPauli

        Raises:
            QiskitError: If it is not possible to express BaseXPPauli in new_precision

        Examples:
            >>> a = BaseXPPauli(
            ... matrix=np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0], dtype=np.int64),
            ... phase_exp=12, precision=8)
            >>> a = a.rescale_precision(new_precision=2)
            >>> a.matrix
            array([[1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]], dtype=int64)
            >>> a._phase_exp
            array([3, dtype=int32])

        See also:
            _rescale_precision
        """
        unique_xp_op = self.unique_vector_rep()
        old_precision = self.precision
        if new_precision < old_precision:
            scale_factor = old_precision // new_precision
        if (new_precision > old_precision) and (new_precision % old_precision > 0):
            raise QiskitError("XP Operator can not be expressed in new_precision.")
        if (new_precision < old_precision) and (
            (old_precision % new_precision > 0)
            or (np.sum(np.mod(unique_xp_op._phase_exp, scale_factor)) > 0)
            or (np.sum(np.mod(unique_xp_op.z, scale_factor)) > 0)
        ):
            raise QiskitError("XP Operator can not be expressed in new_precision.")

        return self._rescale_precision(new_precision)

    def _rescale_precision(self, new_precision: int) -> "BaseXPPauli":
        """Rescale the generalized symplectic vector components
        of BaseXPPauli operator to the new precision. Returns None if the
        rescaling is not possible, else returns the rescaled BaseXPPauli object.

        Note:
            This method is adapted from method XPSetN from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            new_precision: The target precision in which BaseXPPauli is to be expressed

        Returns:
            BaseXPPauli: Resultant of rescaling the precision of BaseXPPauli

        See also:
            unique_vector_rep
        """

        # TODO Currently, if any operator in an XPPauliList can not be
        # rescaled, this function will return None.
        unique_xp_op = self.unique_vector_rep()
        old_precision = unique_xp_op.precision
        matrix = np.empty(shape=np.shape(unique_xp_op.matrix), dtype=np.int64)
        phase_exp = np.empty(shape=np.shape(unique_xp_op._phase_exp))

        if new_precision > old_precision:
            if np.mod(new_precision, old_precision > 0):
                return None
            scale_factor = new_precision // old_precision
            phase_exp = scale_factor * unique_xp_op._phase_exp
            matrix[:, unique_xp_op.num_qubits :] = scale_factor * np.atleast_2d(unique_xp_op.z)

        else:
            scale_factor = old_precision // new_precision
            if (
                (old_precision % new_precision > 0)
                or (np.sum(np.mod(unique_xp_op._phase_exp, scale_factor)) > 0)
                or (np.sum(np.mod(unique_xp_op.z, scale_factor)) > 0)
            ):
                return None
            phase_exp = unique_xp_op._phase_exp // scale_factor
            matrix[:, unique_xp_op.num_qubits :] = np.atleast_2d(unique_xp_op.z) // scale_factor

        matrix[:, 0 : unique_xp_op.num_qubits] = unique_xp_op.x

        return BaseXPPauli(matrix, phase_exp, new_precision)

    def weight(self) -> Union[int, np.ndarray]:
        """Return the weight, i.e. count of qubits where either z or x component is nonzero.

        Note:
            This method is adapted from method XPDistance from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            Union[int, np.ndarray]: Weight of BaseXPPauli

        Examples:
            >>> a = BaseXPPauli(
            ... matrix=np.array([1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 1], dtype=np.int64),
            ... phase_exp = 12, precision = 8)
            >>> a.weight()
            array([5])

        See also:
            _weight
        """
        return self._weight()

    def _weight(self) -> Union[int, np.ndarray]:
        """Return the weight, i.e. count of qubits where either z or x component is nonzero.

        Note:
            This method is adapted from method XPDistance from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            Union[int, np.ndarray]: Weight of BaseXPPauli
        """
        return np.sum(np.logical_or(self.x, self.z), axis=-1)

    def is_diagonal(self) -> np.ndarray:
        """Return True if the XP operator is diagonal.

        Note:
            This method is adapted from method XPisDiag from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            np.ndarray: True, if BaseXPPauli is diagonal

        Examples:
            >>> a = BaseXPPauli(
            ... matrix=np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 3, 3, 3, 3], dtype=np.int64),
            ... phase_exp=0, precision=8)
            >>> a.is_diagonal()
            array([True])

            >>> b = BaseXPPauli(
            ... matrix=np.array([0, 1, 0, 1, 0, 0, 0, 0, 1, 3, 3, 3, 3, 3], dtype=np.int64),
            ... phase_exp=12, precision=8)
            >>> b.is_diagonal()
            array([False])

        See also:
            _is_diagonal
        """
        return self._is_diagonal()

    def _is_diagonal(self) -> np.ndarray:
        """Return True if the XP operator is diagonal.

        Note:
            This method is adapted from method XPisDiag from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            np.ndarray: True, if BaseXPPauli is diagonal
        """
        return np.where(np.sum(self.x, axis=-1) == 0, True, False)

    def antisymmetric_op(self) -> "BaseXPPauli":
        """Return the antisymmetric operator corresponding to the
        z component of XP operator, only if x component is 0.

        Note:
            This method is adapted from method XPD from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            BaseXPPauli: Antisymmetric operator corresponding to BaseXPPauli, if x is 0

        Examples:
            >>> a = BaseXPPauli(
            ... matrix=np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 3, 3, 3], dtype=np.int64),
            ... phase_exp=0, precision=8)
            >>> value = a.antisymmetric_op()
            >>> value.matrix
            array([0, 0, 0, 0, 0, 0, 0, 0, -1, -2, -3, -3, -3, -3], dtype=int64)
            >>> value._phase_exp
            array([15])

        See also:
            _antisymmetric_op
        """
        return self._antisymmetric_op()

    def _antisymmetric_op(self) -> "BaseXPPauli":
        """Return the antisymmetric operator corresponding to the
        z component of XP operator, only if x component is 0, else it returns
        None.

        Note:
            This method is adapted from method XPD from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            BaseXPPauli: Antisymmetric operator corresponding to BaseXPPauli, if x is 0
        """

        if np.any(self.x):
            # TODO should there be an assertion here?
            return None

        phase_exp = np.sum(self.z, axis=-1)
        x = np.zeros(np.shape(self.z))
        matrix = np.concatenate((x, -self.z), axis=-1)

        return BaseXPPauli(matrix=matrix, phase_exp=phase_exp, precision=self.precision)

    def power(self, n: int) -> "BaseXPPauli":
        """Return the XP operator of specified precision raised to the power n.

        Note:
            This method is adapted from method XPPower from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            n: The power to which BaseXPPauli is to be raised

        Returns:
            BaseXPPauli: BaseXPPauli raised to the power n

        Examples:
            >>> a = BaseXPPauli(
            ... matrix=np.array([1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 1], dtype=np.int64),
            ... phase_exp=12, precision=8)
            >>> value = a.power(n=5)
            >>> value.matrix
            array([1, 1, 1, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 5], dtype=int64)
            >>> value._phase_exp
            array([8])

        See also:
            _power
        """
        return self._power(n)

    def _power(self, n: int) -> "BaseXPPauli":
        """Return the XP operator of specified precision raised to the power n.

        Note:
            This method is adapted from method XPPower from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Args:
            n: The power to which BaseXPPauli is to be raised

        Returns:
            BaseXPPauli: BaseXPPauli raised to the power n

        See also:
            _unique_vector_rep
        """
        # TODO at present, this function only handles positive powers. If it is
        # supposed to calculate inverses as well, that functionality needs to
        # be coded.

        a = np.mod(n, 2)

        x = np.multiply(self.x, a)
        z = np.multiply(self.z, n)
        phase_exp = np.multiply(self._phase_exp, n)
        matrix = np.concatenate((x, z), axis=-1)
        first = BaseXPPauli(matrix=matrix, phase_exp=phase_exp, precision=self.precision)

        x = np.zeros(np.shape(self.z))
        z = np.multiply((n - a), np.multiply(self.x, self.z))
        matrix = np.concatenate((x, z), axis=-1)
        second = BaseXPPauli(matrix=matrix, precision=self.precision).antisymmetric_op()

        product = BaseXPPauli(
            matrix=first.matrix + second.matrix,
            phase_exp=first._phase_exp + second._phase_exp,
            precision=self.precision,
        )

        return product._unique_vector_rep()

    def degree(self) -> np.ndarray:
        """Return the degree of XP operator.

        Note:
            This method is adapted from method XPDegree from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            np.ndarray: Degree of BaseXPPauli

        Examples:
        >>> a = BaseXPPauli(matrix=np.array([0, 0, 0, 2, 1, 0], dtype=np.int64),
        ... phase_exp=2, precision=4)
        >>> a.degree()
        array([4], dtype=int64)

        See also:
            _degree
        """
        return self._degree()

    def _degree(self) -> np.ndarray:
        """Return the degree of XP operator.

        Note:
            This method is adapted from method XPDegree from XPFpackage:
            https://github.com/m-webster/XPFpackage, originally developed by
            Mark Webster. The original code is licensed under the GNU General
            Public License v3.0 and Mark Webster has given permission to use
            the code under the Apache License v2.0.

        Returns:
            np.ndarray: Degree of BaseXPPauli

        See also:
            is_diagonal
        """

        gcd = np.gcd(self.z, self.precision)
        precision_by_gcd = np.atleast_2d(np.floor_divide(self.precision, gcd))
        lcm = np.atleast_2d(precision_by_gcd)[:, 0]
        for i, val in enumerate(precision_by_gcd):
            for j in val:
                lcm[i] = np.lcm(lcm[i], j)

        square = BaseXPPauli.compose(self, self)
        if not isinstance(square, type(self)):
            square = type(self)(square)
        gcd_square = np.gcd(square.z, square.precision)
        precision_by_gcd_square = np.atleast_2d(np.floor_divide(square.precision, gcd_square))
        lcm_square = np.atleast_2d(precision_by_gcd_square)[:, 0]
        for i, val in enumerate(precision_by_gcd_square):
            for j in val:
                lcm_square[i] = np.lcm(lcm_square[i], j)

        lcm_square = 2 * lcm_square

        # Do not modify the logic of this function. Naively, it looks like the
        # algorithm that is used when the XP operator is non-diagonal (which
        # involves squaring) gives the correct output for diagonal XP operators
        # as well. However, that is not true. Counter example given by Mark
        # Webster is the operator -I, where the faulty method would give the
        # degree 2, while the actual degree is 1.
        return np.where(self.is_diagonal(), lcm, lcm_square)


# ---------------------------------------------------------------------
# Evolution by Clifford gates
# ---------------------------------------------------------------------


# pylint: disable=unused-argument
def _evolve_h(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _evolve_s(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _evolve_sdg(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _evolve_i(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    """_summary_"""
    pass


def _evolve_x(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    """_summary_"""
    pass


def _evolve_y(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    """_summary_"""
    pass


def _evolve_z(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    """_summary_"""
    pass


def _evolve_cx(base_xp_pauli: "BaseXPPauli", qctrl: int, qtrgt: int) -> "BaseXPPauli":
    """_summary_"""
    pass


def _evolve_cz(  # pylint: disable=invalid-name
    base_xp_pauli: "BaseXPPauli", q1: int, q2: int  # pylint: disable=invalid-name
) -> "BaseXPPauli":
    """_summary_"""
    pass


def _evolve_cy(base_xp_pauli: "BaseXPPauli", qctrl: int, qtrgt: int) -> "BaseXPPauli":
    """_summary_"""
    pass


def _evolve_swap(  # pylint: disable=invalid-name
    base_xp_pauli: "BaseXPPauli", q1: int, q2: int  # pylint: disable=invalid-name
) -> "BaseXPPauli":
    """_summary_"""
    pass
