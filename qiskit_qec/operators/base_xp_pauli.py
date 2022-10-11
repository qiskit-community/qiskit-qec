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
"""Module for base XP pauli"""

import copy
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
from qiskit_qec.linear import matrix as mt
from qiskit_qec.linear.symplectic import symplectic_product
from qiskit_qec.utils import xp_pauli_rep


# pylint: disable=no-member
class BaseXPPauli(BaseOperator, AdjointMixin, MultiplyMixin):
    r"""Symplectic representation of a list of N-qubit XP operators with phases using
    numpy arrays for symplectic matrices and phase vectors.

    Base class for XPPauli and XPPauliList.
    """

    # External string formats used when displaying Pauli's as strings
    EXTERNAL_TENSOR_ENCODING = xp_pauli_rep.DEFAULT_EXTERNAL_TENSOR_ENCODING
    EXTERNAL_PHASE_ENCODING = xp_pauli_rep.DEFAULT_EXTERNAL_PHASE_ENCODING
    EXTERNAL_XP_PAULI_ENCODING = EXTERNAL_PHASE_ENCODING + EXTERNAL_TENSOR_ENCODING

    EXTERNAL_SYNTAX = xp_pauli_rep.PRODUCT_SYNTAX
    EXTERNAL_QUBIT_ORDER = xp_pauli_rep.DEFAULT_QUBIT_ORDER

    PRINT_PHASE_ENCODING = None

    def __init__(
        self,
        matrix: Union[np.ndarray, None] = None,
        phase_exp: Union[None, np.ndarray, np.integer] = None,
        order: str = "xz",
    ) -> None:
        """A BaseXPPauli object represents a list N-qubit Pauli operators with phases.
        Numpy arrays are used to represent the symplectic matrix represention of these
        Paulis. The phases of the Paulis are stored encoded. The phases of the Pauli
        operators are internally encoded in the '-iZX' Pauli encoding (See the xp_pauli_rep
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
            >>> base_xp_pauli = BaseXPPauli(matrix)

        See Also:
            Pauli, PauliList
        """

        if matrix is None or matrix.size == 0:
            matrix = np.empty(shape=(0, 0), dtype=np.int8)
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
        pass

    @x.setter
    def x(self, val: np.ndarray):
        pass

    # @final Add when python >= 3.8
    @property
    def _x(self):  # pylint: disable=invalid-name
        pass

    # @final Add when python >= 3.8
    @_x.setter
    def _x(self, val):  # pylint: disable=invalid-name
        pass

    @property
    def z(self):
        pass

    @z.setter
    def z(self, val):
        pass

    # @final Add when python >= 3.8
    @property
    def _z(self):  # pylint: disable=invalid-name
        pass

    # @final Add when python >= 3.8
    @_z.setter
    def _z(self, val):  # pylint: disable=invalid-name
        pass

    @property
    def num_y(self):
        pass

    @property
    def tensor_encoding(self):
        pass

    @classmethod
    def set_tensor_encoding(cls, encoding: str = xp_pauli_rep.DEFAULT_EXTERNAL_TENSOR_ENCODING):
        pass

    @property
    def phase_encoding(self):
        pass

    @classmethod
    def set_phase_encoding(cls, encoding: str = xp_pauli_rep.DEFAULT_EXTERNAL_PHASE_ENCODING):
        pass

    @property
    def pauli_encoding(self):
        pass

    @classmethod
    def set_pauli_encoding(cls, encoding: str = xp_pauli_rep.DEFAULT_EXTERNAL_XP_PAULI_ENCODING):
        pass

    @property
    def syntax(self):
        pass

    @classmethod
    def set_syntax(cls, syntax_code: Optional[int] = None, syntax_str: Optional[str] = "Product"):
        pass

    @property
    def print_phase_encoding(self):
        pass

    @classmethod
    def set_print_phase_encoding(cls, phase_encoding: Optional[str] = None):
        pass

    @property
    def qubit_order(self):
        pass

    @classmethod
    def set_qubit_order(cls, qubit_order: Optional[str] = None):
        pass

    # ---------------------------------------------------------------------
    # Magic Methods
    # ---------------------------------------------------------------------

    def __imul__(self, other: "BaseXPPauli") -> "BaseXPPauli":
        pass

    def __neg__(self) -> "BaseXPPauli":
        pass

    # ---------------------------------------------------------------------

    def copy(self) -> "BaseXPPauli":
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
        pass

    @staticmethod
    def _compose(
        a: "BaseXPPauli",
        b: "BaseXPPauli",
        qargs: Optional[list] = None,
        front: bool = False,
        inplace: bool = False,
    ) -> "BaseXPPauli":
        pass

    # ---------------------------------------------------------------------

    def tensor(self, other):
        pass

    @staticmethod
    def _tensor(a: "BaseXPPauli", b: "BaseXPPauli"):
        pass

    # ---------------------------------------------------------------------

    def expand(self, other):
        pass

    # ---------------------------------------------------------------------
    # Needed by MultiplyMixin class
    # ---------------------------------------------------------------------

    def _multiply(self, phase, roundit=True) -> "BaseXPPauli":  # pylint: disable=arguments-renamed
        pass

    # Needed by AdjointMixin class

    def conjugate(self, inplace=False) -> "BaseXPPauli":
        pass

    def transpose(self, inplace: bool = False) -> "BaseXPPauli":
        pass

    def commutes(self, other: "BaseXPPauli", qargs: List = None) -> np.ndarray:
        pass

    def _commutes(self, other: "BaseXPPauli", qargs: List = None) -> np.ndarray:
        pass

    # ---------------------------------------------------------------------
    # Extra all_commutes method
    # ---------------------------------------------------------------------

    def all_commutes(self, other: "BaseXPPauli") -> np.ndarray:
        pass

    # ---------------------------------------------------------------------
    # Evolve Methods
    # ---------------------------------------------------------------------

    def evolve(self, other: "BaseXPPauli", qargs: Union[None, List, int] = None, frame: str = "h"):
        pass

    def _evolve_clifford(
        self,
        other: qiskit.quantum_info.operators.symplectic.clifford.Clifford,
        qargs: Union[None, List, int] = None,
        frame: str = "h",
    ):
        pass

    # ---------------------------------------------------------------------
    # Helper Metods
    # ---------------------------------------------------------------------

    def _eq(self, other):
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
        pass

    # ---------------------------------------------------------------------
    # The methods below should be deprecated eventually as they simply refer
    # to the newer versions that are now elsewhere.
    # ---------------------------------------------------------------------

    def _count_y(self) -> np.ndarray:
        pass

    @staticmethod
    def _stack(array: np.ndarray, size: int, vertical: bool = True) -> np.ndarray:
        pass

    @staticmethod
    def _phase_from_complex(coeff: numbers.Complex) -> np.ndarray:
        pass

    # ---------------------------------------------------------------------
    # Apply Clifford to BaseXPPauli
    # ---------------------------------------------------------------------

    def _append_circuit(
        self,
        circuit: Union[Barrier, Delay, QuantumCircuit, Instruction],
        qargs: Optional[List] = None,
    ) -> "BaseXPPauli":
        pass

# ---------------------------------------------------------------------
# Evolution by Clifford gates
# ---------------------------------------------------------------------


def _evolve_h(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    pass


def _evolve_s(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    pass

def _evolve_sdg(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    pass

# pylint: disable=unused-argument
def _evolve_i(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    pass


def _evolve_x(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    pass

def _evolve_y(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    pass

def _evolve_z(base_xp_pauli: "BaseXPPauli", qubit: int) -> "BaseXPPauli":
    pass

def _evolve_cx(base_xp_pauli: "BaseXPPauli", qctrl: int, qtrgt: int) -> "BaseXPPauli":
    pass

def _evolve_cz(  # pylint: disable=invalid-name
    base_xp_pauli: "BaseXPPauli", q1: int, q2: int  # pylint: disable=invalid-name
) -> "BaseXPPauli":
    pass

def _evolve_cy(base_xp_pauli: "BaseXPPauli", qctrl: int, qtrgt: int) -> "BaseXPPauli":
    pass

def _evolve_swap(  # pylint: disable=invalid-name
    base_xp_pauli: "BaseXPPauli", q1: int, q2: int  # pylint: disable=invalid-name
) -> "BaseXPPauli":
    pass
