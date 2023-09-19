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
"""Module for base pauli"""

import copy
import numbers
from typing import List, Optional, Union

import numpy as np

# Must be imported as follows to avoid circular import errors
from qiskit import QiskitError
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.barrier import Barrier
from qiskit.circuit.delay import Delay
from qiskit.circuit.instruction import Instruction
from qiskit.quantum_info.operators.base_operator import BaseOperator
from qiskit.quantum_info.operators.mixins import AdjointMixin, MultiplyMixin

from qiskit_qec.linear import matrix as mt
from qiskit_qec.linear.symplectic import symplectic_product
from qiskit_qec.utils import pauli_rep


# pylint: disable=no-member
class BasePauli(BaseOperator, AdjointMixin, MultiplyMixin):
    r"""Symplectic representation of a list of N-qubit Paulis with phases using
    numpy arrays for symplectic matrices and phase vectors.

    Base class for Pauli and PauliList.
    """

    # External string formats used when displaying Pauli's as strings
    EXTERNAL_TENSOR_ENCODING = pauli_rep.DEFAULT_EXTERNAL_TENSOR_ENCODING
    EXTERNAL_PHASE_ENCODING = pauli_rep.DEFAULT_EXTERNAL_PHASE_ENCODING
    EXTERNAL_PAULI_ENCODING = EXTERNAL_PHASE_ENCODING + EXTERNAL_TENSOR_ENCODING

    EXTERNAL_SYNTAX = pauli_rep.PRODUCT_SYNTAX
    EXTERNAL_QUBIT_ORDER = pauli_rep.DEFAULT_QUBIT_ORDER

    PRINT_PHASE_ENCODING = None

    def __init__(
        self,
        matrix: Union[np.ndarray, None] = None,
        phase_exp: Union[None, np.ndarray, np.integer] = None,
        order: str = "xz",
    ) -> None:
        """A BasePauli object represents a list N-qubit Pauli operators with phases.
        Numpy arrays are used to represent the symplectic matrix represention of these
        Paulis. The phases of the Paulis are stored encoded. The phases of the Pauli
        operators are internally encoded in the '-iZX' Pauli encoding (See the pauli_rep
        module for more details). That is a Pauli operator is represented as symplectic
        vector V and a phase exponent phase_exp such that:

        (-i)^phase_exp Z^z X^x

        where V = [x, z] and phase_exp is a vector of Z_4 elements (0,1,2,3). A list
        of Pauli operators is represented as a symplectic matrix S and a phase exponent
        vector phase_exp such that the rows or S are the symplectic vector representations
        of the Paulis and the phase_exp vector store the phase exponent of each
        associated Pauli Operator.

        Args:
            matrix: Input GF(2) symplectic matrix
            phase_exp (optional): Phase exponent vector for imput matrix. A value of None will
                result in an a complex coefficients of 1 for each Pauli operator. Defaults to None.
            order: Set to 'xz' or 'zx'. Defines which side the x and z parts of the input matrix

        Raises: QiskitError: matrix and phase_exp sizes are not compatible

        Examples:
            >>> matrix = numpy.array([[1,1,0,0],[0,1,0,1]])
            >>> base_pauli = BasePauli(matrix)

        See Also:
            Pauli, PauliList
        """

        if matrix is None or matrix.size == 0:
            matrix = np.empty(shape=(0, 0), dtype=np.bool_)
            phase_exp = np.empty(shape=(0,), dtype=np.int8)
        matrix = np.atleast_2d(matrix)
        num_qubits = matrix.shape[1] >> 1

        if order == "zx":
            nmatrix = np.empty(shape=matrix.shape, dtype=matrix.dtype)
            nmatrix[:, :num_qubits] = matrix[:, num_qubits:]
            nmatrix[:, num_qubits:] = matrix[:, :num_qubits]
            matrix = nmatrix

        self.matrix = matrix
        self._num_paulis = self.matrix.shape[0]
        if phase_exp is None:
            self._phase_exp = np.zeros(shape=(self.matrix.shape[0],), dtype=np.int8)
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
        """Returns the X part of symplectic representation as a 2d matrix.

        Note: The Pauli class over writes this method to return
        a 1d array instead of a 2d array. Use the self._x method
        if a 2d array is needed as _x method is markeded as @final

        Examples:
            >>> matrix = numpy.array([[1,0,0,0],[0,1,1,1]], dtype=numpy.bool_)
            >>> phase_exp = numpy.array([0,1])
            >>> base_pauli = BasePauli(matrix, phase_exp)
            >>> base_pauli.x.astype(int)
            array([[1, 0],
                   [0, 1]])

        See Also:
            _x, z, _z

        """
        return self.matrix[:, : self.num_qubits]

    @x.setter
    def x(self, val: np.ndarray):
        """Sets the X part of symplectic representation

        Args:
            val: GF(2) matrix used to set the X part of
                the symplectic representation

        Examples:
            >>> matrix = numpy.array([[1,0,0,0],[0,1,1,1]], dtype=numpy.bool_)
            >>> phase_exp = numpy.array([0,1])
            >>> base_pauli = BasePauli(matrix, phase_exp)
            >>> base_pauli.x = numpy.array([[1,1],[0,0]], dtype=numpy.bool_)
            >>> base_pauli.x.astype(int)
            array([[1, 1],
                   [0, 0]])

        See Also:
            x, z, _z
        """
        self.matrix[:, : self.num_qubits] = val

    # @final Add when python >= 3.8
    @property
    def _x(self):  # pylint: disable=invalid-name
        """Returns the X part of symplectic representation as a 2d matrix.

        Note: The Pauli class over writes this method to return
        a 1d array instead of a 2d array. Use the self._x method
        if a 2d array is needed as _x method is markeded as @final

        Examples:
            >>> matrix = numpy.array([[1,0,0,0],[0,1,1,1]], dtype=numpy.bool_)
            >>> phase_exp = numpy.array([0,1])
            >>> base_pauli = BasePauli(matrix, phase_exp)
            >>> base_pauli._x.astype(int)
            array([[1, 0],
                   [0, 1]])

        See Also:
            x, z, _z
        """
        return self.matrix[:, : self.num_qubits]

    # @final Add when python >= 3.8
    @_x.setter
    def _x(self, val):  # pylint: disable=invalid-name
        self.matrix[:, : self.num_qubits] = val

    @property
    def z(self):
        """The z array for the symplectic representation."""
        return self.matrix[:, self.num_qubits :]

    @z.setter
    def z(self, val):
        self.matrix[:, self.num_qubits :] = val

    # @final Add when python >= 3.8
    @property
    def _z(self):  # pylint: disable=invalid-name
        """The z array for the symplectic representation."""
        return self.matrix[:, self.num_qubits :]

    # @final Add when python >= 3.8
    @_z.setter
    def _z(self, val):  # pylint: disable=invalid-name
        self.matrix[:, self.num_qubits :] = val

    @property
    def num_y(self):
        """Return the number of Y for each operator"""
        return np.sum(np.logical_and(self.x, self.z), axis=1)

    @property
    def tensor_encoding(self):
        """Return the external symplectic matrix encoding"""
        return BasePauli.EXTERNAL_TENSOR_ENCODING

    @classmethod
    def set_tensor_encoding(cls, encoding: str = pauli_rep.DEFAULT_EXTERNAL_TENSOR_ENCODING):
        """Set the external symplectic matrix format

        Args:
            encoding (optional): Symplectic matrix tensor encoding.
            Defaults to pauli_rep.DEFAULT_EXTERNAL_TENSOR_ENCODING.
        """
        assert encoding in pauli_rep.get_tensor_encodings(), QiskitError(
            f"Invalid symplectic matrix encoding: {encoding}. Must be one \
                of {pauli_rep.get_tensor_encodings()}"
        )
        BasePauli.EXTERNAL_TENSOR_ENCODING = encoding

    @property
    def phase_encoding(self):
        """Return the phase encoding"""
        return BasePauli.EXTERNAL_PHASE_ENCODING

    @classmethod
    def set_phase_encoding(cls, encoding: str = pauli_rep.DEFAULT_EXTERNAL_PHASE_ENCODING):
        """Set the phase encoding

        Args:
            encoding (optional): phase encoding.
            Defaults to pauli_rep.DEFAULT_EXTERNAL_PHASE_ENCODING.
        """
        assert encoding in pauli_rep.get_phase_encodings(), QiskitError(
            f"Invalid phase encoding: {encoding}. Must be one of {pauli_rep.get_phase_encodings()}"
        )
        BasePauli.EXTERNAL_PHASE_ENCODING = encoding

    @property
    def pauli_encoding(self):
        """Pauli format."""
        return BasePauli.EXTERNAL_PAULI_ENCODING

    @classmethod
    def set_pauli_encoding(cls, encoding: str = pauli_rep.DEFAULT_EXTERNAL_PAULI_ENCODING):
        """Set the Pauli encoding

        Args:
            encoding (optional): Pauli encoding.
            Defaults to pauli_rep.DEFAULT_EXTERNAL_PAULI_REP_FORMAT.
        """
        assert encoding in pauli_rep.get_pauli_encodings(), QiskitError(
            f"Invalid pauli encoding: {encoding}. Must be one of {pauli_rep.get_pauli_encodings()}"
        )
        phase_encoding, tensor_encoding = pauli_rep._split_pauli_enc(encoding)
        BasePauli.EXTERNAL_PHASE_ENCODING = phase_encoding
        BasePauli.EXTERNAL_TENSOR_ENCODING = tensor_encoding
        BasePauli.EXTERNAL_PAULI_ENCODING = phase_encoding + tensor_encoding

    @property
    def syntax(self):
        """_summary_"""
        return BasePauli.EXTERNAL_SYNTAX

    @classmethod
    def set_syntax(cls, syntax_code: Optional[int] = None, syntax_str: Optional[str] = "Product"):
        """Sets the global input and output format

        Args:
            syntax_code (Optional[int], optional): sets the syntax of Pauli tensors. Possible inputs
                                                   are 0 for product syntax, 1 for index syntax and 2
                                                   for latex syntax. Defaults to None.
            syntax_str (Optional[str], optional): sets the syntax of Pauli tensors. Possible inputs are
                                                  Product or Latex, if another input is given the syntax
                                                   is set to Order. Defaults to "Product".

        Raises:
            QiskitError: Unknown syntax: {syntax_code}. See pauli_rep for options.
        """
        if syntax_code is None:
            if syntax_str == "Product":
                BasePauli.EXTERNAL_SYNTAX = pauli_rep.PRODUCT_SYNTAX
            elif syntax_str == "Latex":
                BasePauli.EXTERNAL_SYNTAX = pauli_rep.LATEX_SYNTAX
            else:
                BasePauli.EXTERNAL_SYNTAX = pauli_rep.INDEX_SYNTAX
        else:
            if syntax_code not in [0, 1, 2]:
                raise QiskitError("Unknown syntax: {syntax_code}. See pauli_rep for options.")
            BasePauli.EXTERNAL_SYNTAX = syntax_code

    @property
    def print_phase_encoding(self):
        """Prints how the phase will be displayed in when printing."""
        return BasePauli.PRINT_PHASE_ENCODING

    @classmethod
    def set_print_phase_encoding(cls, phase_encoding: Optional[str] = None):
        """_summary_

        Args:
            phase_encoding (Optional[str], optional): _description_. Defaults to None.

        Raises:
            QiskitError: _description_
        """
        if phase_encoding is None or phase_encoding in pauli_rep.PHASE_ENCODINGS:
            BasePauli.PRINT_PHASE_ENCODING = phase_encoding
        else:
            raise QiskitError(
                f"Unknown print phase encoding {phase_encoding}. Encoding \
                must be None or one of {pauli_rep.get_phase_encodings}"
            )

    @property
    def qubit_order(self):
        """Get external qubit order"""
        return BasePauli.EXTERNAL_QUBIT_ORDER

    @classmethod
    def set_qubit_order(cls, qubit_order: Optional[str] = None):
        """Set external qubit order

        Args:
            qubit_order (Optional[str], optional): _description_. Defaults to None.

        Raises:
            QiskitError: _description_
        """
        if qubit_order is None:
            BasePauli.EXTERNAL_QUBIT_ORDER = pauli_rep.DEFAULT_QUBIT_ORDER
        else:
            if qubit_order not in pauli_rep.QUBIT_ORDERS:
                raise QiskitError(f"Unknown qubit order: {qubit_order}")

            BasePauli.EXTERNAL_QUBIT_ORDER = qubit_order

    # ---------------------------------------------------------------------
    # Magic Methods
    # ---------------------------------------------------------------------

    def __imul__(self, other: "BasePauli") -> "BasePauli":
        return self.compose(other, front=True, inplace=True)

    def __neg__(self) -> "BasePauli":
        ret = copy.copy(self)
        ret._phase_exp = np.mod(self._phase_exp + 2, 4)
        return ret

    # ---------------------------------------------------------------------

    def copy(self) -> "BasePauli":
        """Make a deep copy of current operator."""
        # Deepcopy has terrible performance on objects with Numpy arrays
        # attributes so we make a shallow copy and then manually copy the
        # Numpy arrays to efficiently mimic a deepcopy
        ret = copy.copy(self)
        ret.matrix = self.matrix.copy()
        ret._phase_exp = self._phase_exp.copy()
        return ret

    # ---------------------------------------------------------------------
    # BaseOperator methods
    # ---------------------------------------------------------------------
    # Needed by GroupMixin class from BaseOperator class

    # pylint: disable=arguments-differ
    def compose(
        self,
        other: "BasePauli",
        qargs: Optional[list] = None,
        front: bool = False,
        inplace: bool = False,
    ) -> "BasePauli":
        r"""Return the composition of Paulis lists

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

        Args:
            other: BasePauli
            front (bool): (default: False)
            qargs (list or None): Optional, qubits to apply compose on
                                  on (default: None->All).
            inplace (bool): If True update in-place (default: False).

        Returns:
            BasePauli : Compositon of self and other

        Raises:
            QiskitError: if number of qubits of other does not match qargs.
        """

        # Validation
        if qargs is None and other.num_qubits != self.num_qubits:
            raise QiskitError(f"other {type(self).__name__} must be on the same number of qubits.")

        if qargs and other.num_qubits != len(qargs):
            raise QiskitError(
                f"Number of qubits of the other {type(self).__name__} does not match qargs."
            )

        if other._num_paulis not in [1, self._num_paulis]:
            raise QiskitError(
                "Incompatible BasePaulis. Second list must "
                "either have 1 or the same number of Paulis."
            )

        return self._compose(self, other, qargs=qargs, front=front, inplace=inplace)

    @staticmethod
    def _compose(
        a: "BasePauli",
        b: "BasePauli",
        qargs: Optional[list] = None,
        front: bool = False,
        inplace: bool = False,
    ) -> "BasePauli":
        """Returns the composition of two BasePauli objects.

        Args:
            a : BasePauli object
            b : BasePauli object
            qargs (Optional[list], optional): _description_. Defaults to None.
            front (bool, optional): _description_. Defaults to False.
            inplace (bool, optional): _description_. Defaults to False.

        Returns:
            BasePauli: _description_
        """
        if qargs is not None:
            qargs = list(qargs) + [item + a.num_qubits for item in qargs]
            amat = a.matrix[:, qargs]
        else:
            amat = a.matrix
        bmat = b.matrix

        # Calculate the symplectic matrix for the composition
        mat = np.logical_xor(amat, bmat)

        # Calculate the phase of the composition
        phase_exp = a._phase_exp + b._phase_exp
        if front:
            phase_exp += 2 * np.sum(
                np.logical_and(amat[:, : b.num_qubits], bmat[:, b.num_qubits :]), axis=1
            )
        else:
            phase_exp += 2 * np.sum(
                np.logical_and(bmat[:, : b.num_qubits], amat[:, b.num_qubits :]), axis=1
            )
        phase_exp = np.mod(phase_exp, 4)

        if qargs is None:
            if not inplace:
                return BasePauli(mat, phase_exp)
            # Inplace update
            a.matrix = mat
            a._phase_exp = phase_exp
            return a

        # Qargs update
        ret = a if inplace else a.copy()
        ret.matrix[:, qargs] = mat
        ret._phase_exp = np.mod(phase_exp, 4)
        return ret

    # ---------------------------------------------------------------------

    def tensor(self, other):
        return self._tensor(self, other)

    @staticmethod
    def _tensor(a: "BasePauli", b: "BasePauli"):
        """Returns the generalized tensor product of Paulis that are represented by
        symplectic matrice and a phase exponent

        Args:
            a (BasePauli): List of Pauli operators
            b (BasePauli): List of Pauli operators

        if a=[A_1,A_2,...,A_k] and b=[B_1,B_2,...,B_v] then

        a tensor b = [A_1 tensor B_1, A_1 tensor B_2, ..., A_k tensor B_v]

        Returns:
            [PaulisBase]: a gtensor b
        """

        x1 = mt.istack(a.x, b._num_paulis, True)
        x2 = mt.istack(b.x, a._num_paulis, False)
        z1 = mt.istack(a.z, b._num_paulis, True)
        z2 = mt.istack(b.z, a._num_paulis, False)
        phase1_exp = (
            np.vstack(b._num_paulis * [a._phase_exp])
            .transpose(1, 0)
            .reshape(a._num_paulis * b._num_paulis)
        )
        phase2_exp = mt.istack(b._phase_exp, a._num_paulis)

        xz_mat = np.hstack((x2, x1, z2, z1))
        phase_exp = np.mod(phase1_exp + phase2_exp, 4)

        return BasePauli(xz_mat, phase_exp)

    # ---------------------------------------------------------------------

    def expand(self, other):
        return self._tensor(other, self)

    # ---------------------------------------------------------------------
    # Needed by MultiplyMixin class
    # ---------------------------------------------------------------------

    def _multiply(self, phase, roundit=True) -> "BasePauli":  # pylint: disable=arguments-renamed
        """Return the {cls} phase * self where phase is in ``[1, -1j, -1, 1j]``.

        Args:
            phase (complex): a complex number(s) in ``[1, -1j, -1, 1j]``.
            roundit (bool): Set True to round components of other. Default=True
        Returns:
            {cls}: the {cls} phase * self.

        Raises:
            QiskitError: if the phase is not in the set ``[1, -1j, -1, 1j]``.
        """.format(
            cls=type(self).__name__
        )
        phase_exp = pauli_rep.cpx2exp(
            phase, output_encoding=pauli_rep.INTERNAL_PHASE_ENCODING, roundit=roundit
        )
        return BasePauli(self.matrix, np.mod(self._phase_exp + phase_exp, 4))

    # Needed by AdjointMixin class

    def conjugate(self, inplace=False) -> "BasePauli":
        """Return the conjugate of each Pauli in the list.

        Args:
            inplace (boolean) : If True will modify inplace. Default: False,

        Returns:
            {cls} : a new {cls} which has phases conjugates (if replace=False) or will change
            the phase of the clasing instance if replace=True
        """
        new_phase_exp = pauli_rep.exp2exp(self._phase_exp, output_encoding="i")
        if not inplace:
            return BasePauli(self.matrix, new_phase_exp).copy()
        elif new_phase_exp == self._phase_exp:
            return self
        else:
            self._phase_exp = new_phase_exp
            return self

    def transpose(self, inplace: bool = False) -> "BasePauli":
        """Return the transpose of each Pauli in the list."""
        new_phase_exp = np.mod(self._phase_exp + 2 * (self.num_y % 2), 4)
        if not inplace:
            return BasePauli(self.matrix, new_phase_exp).copy()
        elif new_phase_exp == self._phase_exp:
            return self
        else:
            self._phase_exp = new_phase_exp
            return self

    def commutes(self, other: "BasePauli", qargs: List = None) -> np.ndarray:
        """Return True if Pauli that commutes with other.

        Args:
            other (PaulisBase): another PaulisBase operator.
            qargs (list): qubits to apply dot product on (default: None).

        Returns:
            np.array: Boolean array of True if Pauli's commute, False if
                      they anti-commute.

        Raises:
            QiskitError: if number of qubits of other does not match qargs.
        """
        if qargs is not None and len(qargs) != other.num_qubits:
            raise QiskitError(
                f"Number of qubits of other Pauli does not match number of \
                qargs ({other.num_qubits} != {len(qargs)})."
            )
        if qargs is None and self.num_qubits != other.num_qubits:
            raise QiskitError(
                "Number of qubits of other Pauli does not match the current \
                Pauli ({other.num_qubits} != {self.num_qubits})."
            )
        return self._commutes(other, qargs=qargs)

    def _commutes(self, other: "BasePauli", qargs: List = None) -> np.ndarray:
        if qargs is not None:
            inds = list(qargs)
            sinds = [index + self.num_qubits for index in inds]
            # x1, z1 = self.x[:, inds], self.z[:, inds]
            x1, z1 = self.matrix[:, inds], self.matrix[:, sinds]
        else:
            # x1, z1 = self.x, self.z
            x1, z1 = self.matrix[:, : self.num_qubits], self.matrix[:, self.num_qubits :]

        return np.logical_not(np.sum((x1 & other.z) ^ (z1 & other.x), axis=1) % 2)

    # ---------------------------------------------------------------------
    # Extra all_commutes method
    # ---------------------------------------------------------------------

    def all_commutes(self, other: "BasePauli") -> np.ndarray:
        """_summary_

        Args:
            other (BasePauli): _description_

        Returns:
            np.ndarray: _description_
        """
        return np.logical_not(symplectic_product(self.matrix, other.matrix))

    # ---------------------------------------------------------------------
    # Evolve Methods
    # ---------------------------------------------------------------------

    def evolve(self, other: "BasePauli", qargs: Union[None, List, int] = None, frame: str = "h"):
        r"""Heisenberg picture evolution of a Pauli by a Clifford.

        This returns the Pauli :math:`P^\prime = C^\dagger.P.C`.

        By choosing the parameter frame='s', this function returns the Schrödinger evolution of the Pauli
        :math:`P^\prime = C.P.C^\dagger`. This option yields a faster calculation.

        Args:
            other (BasePauli or QuantumCircuit): The Clifford circuit to evolve by.
            qargs (list): a list of qubits to apply the Clifford to.
            frame (string): 'h' for Heisenberg or 's' for Schrödinger framework.

        Returns:
            BasePauli: the Pauli :math:`C^\dagger.P.C`.

        Raises:
            QiskitError: if the Clifford number of qubits and qargs don't match.
        """
        # Check dimension
        if qargs is not None and len(qargs) != other.num_qubits:
            raise QiskitError(
                f"Incorrect number of qubits for Clifford circuit \
                    ({other.num_qubits} != {len(qargs)})."
            )
        if qargs is None and self.num_qubits != other.num_qubits:
            raise QiskitError(
                f"Incorrect number of qubits for Clifford circuit \
                    ({other.num_qubits} != {self.num_qubits})."
            )

        # Evolve via Pauli
        if isinstance(other, BasePauli):
            if frame == "s":
                ret = self.compose(other, qargs=qargs)
                ret = ret.compose(other.adjoint(), front=True, qargs=qargs)
            else:
                ret = self.compose(other.adjoint(), qargs=qargs)
                ret = ret.compose(other, front=True, qargs=qargs)
            return ret

        # Evolve by the inverse circuit to compute C^dg.P.C
        if frame == "s":
            return self.copy()._append_circuit(other, qargs=qargs)
        return self.copy()._append_circuit(other.inverse(), qargs=qargs)

    # ---------------------------------------------------------------------
    # Helper Metods
    # ---------------------------------------------------------------------

    def _eq(self, other):
        """Entrywise comparison of Pauli equality."""
        return (
            self.num_qubits == other.num_qubits
            and np.all(np.mod(self._phase_exp, 4) == np.mod(other._phase_exp, 4))
            and np.all(self.matrix == other.matrix)
        )

    # ---------------------------------------------------------------------
    # Class methods for creating labels -> Uses pauli_rep suite of methods
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
        """Returns the string representatiojn for a Pauli or Paulis.

        Args:
            output_pauli_encoding (optional): Encoding used to represent phases.
                A value of None will result in complex phases notation. Defaults
                to None which will in turn use BasePauli.EXTERNAL_PAULI_ENCODING.
            no_phase (optional): When set to True, no phase will appear no matter
                what encoding is selected. So the symplectic matrix [1, 1] will produce
                the operator Y in 'XZY' encoding but also (XZ) in the 'XZ' encoding which
                are different operators if phases are considered. Defaults to False.
            return_phase (optional): If True return the adjusted phase for the coefficient
                of the returned Pauli label. This can be used even if
                ``full_group=False``.
            syntax (optional): Syntax of pauli tensor. Values are
                PRODUCT_SYNTAX = 0 and INDEX_SYNTAX=1. Defaults to INDEX_SYNTAX.
            qubit_order (optional): Order in which qubits are read. options aree
                "right-to-left" and "left-to-right". Defaults to "right-to-left".
            index_start (optional): Lowest value for index in index syntax tensors.
                Defaults to 0
            squeeze (optional): Squeezes the list of reults to a scalar if the number
                of Paulis is one. Defaults to True.
            index_str (optional): String that get inserted between operator and numbers in
                index format. Default is "".


        Returns:
            str: the Pauli label(string) from the full Pauli group (if ``no_phase=False``) or
                from the unsigned Pauli group (if ``no_phase=True``).
            Tuple[str or List[str], Any or List[Any]]: if ``return_phase=True`` returns a
                tuple of the Pauli
                label (from either the full or unsigned Pauli group) and
                the phase ``q`` for the coefficient :math:`(-i)^(q + x.z)`
                for the label from the full Pauli group.
        """
        if output_pauli_encoding is None:
            output_phase_encoding = BasePauli.PRINT_PHASE_ENCODING
            output_tensor_encoding = BasePauli.EXTERNAL_TENSOR_ENCODING
        else:
            output_phase_encoding, output_tensor_encoding = pauli_rep.split_pauli_enc(
                output_pauli_encoding
            )
        if syntax is None:
            syntax = BasePauli.EXTERNAL_SYNTAX

        if qubit_order is None:
            qubit_order = BasePauli.EXTERNAL_QUBIT_ORDER

        pauli_str = pauli_rep.symplectic2str(
            self.matrix,
            self._phase_exp,
            output_phase_encoding=output_phase_encoding,
            no_phase=no_phase,
            output_tensor_encoding=output_tensor_encoding,
            syntax=syntax,
            qubit_order=qubit_order,
            index_start=index_start,
            same_type=squeeze,
            index_str=index_str,
        )

        if return_phase:
            phase_exp = pauli_rep.change_pauli_encoding(
                self._phase_exp,
                self.num_y,
                input_pauli_encoding=pauli_rep.INTERNAL_PAULI_ENCODING,
                output_pauli_encoding=output_pauli_encoding,
            )
            return pauli_str, phase_exp
        else:
            return pauli_str

    # ---------------------------------------------------------------------
    # The methods below should be deprecated eventually as they simply refer
    # to the newer versions that are now elsewhere.
    # ---------------------------------------------------------------------

    def _count_y(self) -> np.ndarray:
        """Count the number of Y Pauli's"""
        return self.num_y

    @staticmethod
    def _stack(array: np.ndarray, size: int, vertical: bool = True) -> np.ndarray:
        return mt.istack(array, size, interleave=vertical)

    @staticmethod
    def _phase_from_complex(coeff: numbers.Complex) -> np.ndarray:
        return pauli_rep.cpxstr2exp(coeff, encoding=pauli_rep.INTERNAL_PHASE_ENCODING)

    # ---------------------------------------------------------------------
    # Apply Clifford to BasePauli
    # ---------------------------------------------------------------------

    def _append_circuit(
        self,
        circuit: Union[Barrier, Delay, QuantumCircuit, Instruction],
        qargs: Optional[List] = None,
    ) -> "BasePauli":
        """Update BasePauli inplace by applying a Clifford circuit.

        Args:
            circuit (QuantumCircuit or Instruction): the gate or composite gate to apply.
            qargs (list or None): The qubits to apply gate to.

        Returns:
            BasePauli: the updated Pauli.

        Raises:
            QiskitError: if input gate cannot be decomposed into Clifford gates.
        """
        if isinstance(circuit, (Barrier, Delay)):
            return self

        if qargs is None:
            qargs = list(range(self.num_qubits))

        if isinstance(circuit, QuantumCircuit):
            gate = circuit.to_instruction()
        else:
            gate = circuit

        # Basis Clifford Gates
        basis_1q = {
            "i": _evolve_i,
            "id": _evolve_i,
            "iden": _evolve_i,
            "x": _evolve_x,
            "y": _evolve_y,
            "z": _evolve_z,
            "h": _evolve_h,
            "s": _evolve_s,
            "sdg": _evolve_sdg,
            "sinv": _evolve_sdg,
        }
        basis_2q = {"cx": _evolve_cx, "cz": _evolve_cz, "cy": _evolve_cy, "swap": _evolve_swap}

        # Non-Clifford gates
        non_clifford = ["t", "tdg", "ccx", "ccz"]

        if isinstance(gate, str):
            # Check if gate is a valid Clifford basis gate string
            if gate not in basis_1q and gate not in basis_2q:
                raise QiskitError(f"Invalid Clifford gate name string {gate}")
            name = gate
        else:
            # Assume gate is an Instruction
            name = gate.name

        # Apply gate if it is a Clifford basis gate
        if name in non_clifford:
            raise QiskitError(f"Cannot update Pauli with non-Clifford gate {name}")
        if name in basis_1q:
            if len(qargs) != 1:
                raise QiskitError("Invalid qubits for 1-qubit gate.")
            return basis_1q[name](self, qargs[0])
        if name in basis_2q:
            if len(qargs) != 2:
                raise QiskitError("Invalid qubits for 2-qubit gate.")
            return basis_2q[name](self, qargs[0], qargs[1])

        # If not a Clifford basis gate we try to unroll the gate and
        # raise an exception if unrolling reaches a non-Clifford gate.
        if gate.definition is None:
            raise QiskitError(f"Cannot apply Instruction: {gate.name}")
        if not isinstance(gate.definition, QuantumCircuit):
            raise QiskitError(
                f"{gate.name} instruction definition is {type(gate.definition)}; \
                    expected QuantumCircuit"
            )

        flat_instr = gate.definition
        bit_indices = {
            bit: index
            for bits in [flat_instr.qubits, flat_instr.clbits]
            for index, bit in enumerate(bits)
        }

        for instr, qregs, cregs in flat_instr:
            if cregs:
                raise QiskitError(
                    f"Cannot apply Instruction with classical registers: {instr.name}"
                )
            # Get the integer position of the flat register
            new_qubits = [qargs[bit_indices[tup]] for tup in qregs]
            self._append_circuit(instr, new_qubits)

        # Since the individual gate evolution functions don't take mod
        # of phase we update it at the end
        self._phase_exp %= 4
        return self


# ---------------------------------------------------------------------
# Evolution by Clifford gates
# ---------------------------------------------------------------------


def _evolve_h(base_pauli: "BasePauli", qubit: int) -> "BasePauli":
    """Update P -> H.P.H"""
    x = base_pauli.matrix[:, qubit].copy()
    z = base_pauli.matrix[:, qubit + base_pauli.num_qubits].copy()
    base_pauli.matrix[:, qubit] = z
    base_pauli.matrix[:, qubit + base_pauli.num_qubits] = x
    base_pauli._phase_exp += 2 * np.logical_and(x, z).T
    return base_pauli


def _evolve_s(base_pauli: "BasePauli", qubit: int) -> "BasePauli":
    """Update P -> S.P.Sdg"""
    x = base_pauli.matrix[:, qubit]
    base_pauli.matrix[:, qubit + base_pauli.num_qubits] ^= x
    base_pauli._phase_exp += x.T
    return base_pauli


def _evolve_sdg(base_pauli: "BasePauli", qubit: int) -> "BasePauli":
    """Update P -> Sdg.P.S"""
    x = base_pauli.matrix[:, qubit]
    base_pauli.matrix[:, qubit + base_pauli.num_qubits] ^= x
    base_pauli._phase_exp -= x.T
    return base_pauli


# pylint: disable=unused-argument
def _evolve_i(base_pauli: "BasePauli", qubit: int) -> "BasePauli":
    """Update P -> P"""
    return base_pauli


def _evolve_x(base_pauli: "BasePauli", qubit: int) -> "BasePauli":
    """Update P -> X.P.X"""
    base_pauli._phase_exp += 2 * base_pauli.matrix[:, qubit + base_pauli.num_qubits].T
    return base_pauli


def _evolve_y(base_pauli: "BasePauli", qubit: int) -> "BasePauli":
    """Update P -> Y.P.Y"""
    base_pauli._phase_exp += (
        2 * base_pauli.matrix[:, qubit].T
        + 2 * base_pauli.matrix[:, qubit + base_pauli.num_qubits].T
    )
    return base_pauli


def _evolve_z(base_pauli: "BasePauli", qubit: int) -> "BasePauli":
    """Update P -> Z.P.Z"""
    base_pauli._phase_exp += 2 * base_pauli.matrix[:, qubit].T
    return base_pauli


def _evolve_cx(base_pauli: "BasePauli", qctrl: int, qtrgt: int) -> "BasePauli":
    """Update P -> CX.P.CX"""
    base_pauli.matrix[:, qtrgt] ^= base_pauli.matrix[:, qctrl]
    base_pauli.matrix[:, qctrl + base_pauli.num_qubits] ^= base_pauli.matrix[
        :, qtrgt + base_pauli.num_qubits
    ]
    return base_pauli


def _evolve_cz(  # pylint: disable=invalid-name
    base_pauli: "BasePauli", q1: int, q2: int  # pylint: disable=invalid-name
) -> "BasePauli":
    """Update P -> CZ.P.CZ"""
    x1 = base_pauli.matrix[:, q1].copy()
    x2 = base_pauli.matrix[:, q2].copy()
    base_pauli.matrix[:, q1 + base_pauli.num_qubits] ^= x2
    base_pauli.matrix[:, q2 + base_pauli.num_qubits] ^= x1
    base_pauli._phase_exp += 2 * np.logical_and(x1, x2).T
    return base_pauli


def _evolve_cy(base_pauli: "BasePauli", qctrl: int, qtrgt: int) -> "BasePauli":
    """Update P -> CY.P.CY"""
    x1 = base_pauli.matrix[:, qctrl].copy()
    x2 = base_pauli.matrix[:, qtrgt].copy()
    z2 = base_pauli.matrix[:, qtrgt + base_pauli.num_qubits].copy()
    base_pauli.matrix[:, qtrgt] ^= x1
    base_pauli.matrix[:, qtrgt + base_pauli.num_qubits] ^= x1
    base_pauli.matrix[:, qctrl + base_pauli.num_qubits] ^= np.logical_xor(x2, z2)
    base_pauli._phase_exp += x1 + 2 * np.logical_and(x1, x2).T
    return base_pauli


def _evolve_swap(  # pylint: disable=invalid-name
    base_pauli: "BasePauli", q1: int, q2: int  # pylint: disable=invalid-name
) -> "BasePauli":
    """Update P -> SWAP.P.SWAP"""
    x1 = base_pauli.matrix[:, q1].copy()
    z1 = base_pauli.matrix[:, q1 + base_pauli.num_qubits].copy()
    base_pauli.matrix[:, q1] = base_pauli.matrix[:, q2]
    base_pauli.matrix[:, q1 + base_pauli.num_qubits] = base_pauli.matrix[
        :, q2 + base_pauli.num_qubits
    ]
    base_pauli.matrix[:, q2] = x1
    base_pauli.matrix[:, q2 + base_pauli.num_qubits] = z1
    return base_pauli
