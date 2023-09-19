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
"""Module for Pauli"""
from typing import Any, Dict, List, Optional, Union

import numpy as np
from qiskit.circuit import Instruction, QuantumCircuit
from qiskit.circuit.barrier import Barrier
from qiskit.circuit.delay import Delay
from qiskit.circuit.library.generalized_gates import PauliGate
from qiskit.circuit.library.standard_gates import IGate, XGate, YGate, ZGate
from qiskit.exceptions import QiskitError
from qiskit.quantum_info.operators.mixins import generate_apidocs
from qiskit.quantum_info.operators.scalar_op import ScalarOp
from qiskit.utils.deprecation import deprecate_function

from qiskit_qec.operators.base_pauli import BasePauli
from qiskit_qec.utils import pauli_rep


class Pauli(BasePauli):
    """`Pauli` inherits from `BasePauli`"""

    # Set the max Pauli string size before truncation
    _truncate__ = 50

    # Pauli (x,z) lookup table to string
    pltb_str = {(0, 0): "I", (1, 0): "X", (0, 1): "Z", (1, 1): "Y"}
    # Pauli (x,z) lookup table to int
    pltb_int = {(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 3}

    def __init__(
        self,
        data: Any,
        *,
        x: Union[List, np.ndarray, None] = None,
        z: Union[List, np.ndarray, None] = None,
        phase_exp: Union[int, np.ndarray, None] = None,
        input_pauli_encoding: str = BasePauli.EXTERNAL_PAULI_ENCODING,
        label: Optional[str] = None,
        input_qubit_order: str = "right-to-left",
        tuple_order: str = "zx",
    ):
        """Pauli Init

        Args:
            data (str): Still in progress
            x ([type], optional): [description]. Defaults to None.
            z ([type], optional): [description]. Defaults to None.
            phase_exponent ([type], optional): [description]. Defaults to None.
            stype (str, optional): [description]. Defaults to "numpy".
            label ([type], optional): [description]. Defaults to None.
            input_qubit_order (str, optional): [description]. Defaults to "right-to-left".
            tuple_prder (optional): Defaults to 'zx'

        Raises:
            QiskitError: Something went wrong.
        """
        if isinstance(data, np.ndarray):
            matrix = np.atleast_2d(data)
            if phase_exp is None:
                phase_exp = 0
        elif isinstance(data, BasePauli):
            matrix = data.matrix[:, :]
            phase_exp = data._phase_exp[:]
        elif isinstance(data, (tuple, list)):
            if len(data) not in [2, 3]:
                raise QiskitError(
                    "Invalid input tuple for Pauli, input tuple must be `(z, x, phase)` or `(z, x)`"
                )
            if tuple_order not in ["zx", "xz"]:
                raise QiskitError(f"`tuple_order` {tuple_order} not valid")
            if len(data) == 3:
                if tuple_order == "zx":
                    matrix, phase_exp = pauli_rep.from_split_array(
                        data[1], data[0], data[2], input_pauli_encoding=input_pauli_encoding
                    )
                else:
                    matrix, phase_exp = pauli_rep.from_split_array(
                        *data, input_pauli_encoding=input_pauli_encoding
                    )
            elif tuple_order == "zx":
                matrix, phase_exp = pauli_rep.from_split_array(
                    data[1], data[0], 0, input_pauli_encoding=input_pauli_encoding
                )
            else:
                matrix, phase_exp = pauli_rep.from_array(
                    *data, 0, input_pauli_encoding=input_pauli_encoding
                )

        elif isinstance(data, str):
            matrix, phase_exp = pauli_rep.str2symplectic(data, qubit_order=input_qubit_order)
        elif isinstance(data, ScalarOp):
            matrix, phase_exp = pauli_rep.scalar_op2symplectic(
                data, output_encoding=pauli_rep.INTERNAL_PHASE_ENCODING
            )
        elif isinstance(data, (QuantumCircuit, Instruction)):
            matrix, phase_exp = self.instrs2symplectic(data)
        elif x is not None:  # DEPRECATED
            if z is None:
                # Using old Pauli initialization with positional args instead of kwargs
                z = data
            matrix, phase_exp = self._from_split_array_deprecated(z, x, phase_exp)

        elif label is not None:  # DEPRECATED
            matrix, phase_exp = self._from_label_deprecated(label, qubit_order=input_qubit_order)
        else:
            raise QiskitError("Invalid input data for Pauli.")

        # Initialize BasePauli
        if matrix.shape[0] != 1:
            raise QiskitError("Input is not a single Pauli")

        super().__init__(matrix, phase_exp)
        self.vlist = self.matrix[0].tolist()

    # ---------------------------------------------------------------------
    # Property Methods
    # ---------------------------------------------------------------------

    @property
    def name(self):
        """Unique string identifier for operation type."""
        return "pauli"

    @property
    def num_clbits(self):
        """Number of classical bits."""
        return 0

    @property
    def settings(self) -> Dict:
        """Return settings."""
        return {"data": self.to_label()}

    @property
    def num_y(self):
        """Return the number of Y for each operator"""
        return np.sum(np.logical_and(self.x, self.z), axis=0)

    @classmethod
    def set_truncation(cls, val):
        """Set the max number of Pauli characters to display before truncation/

        Args:
            val (int): the number of characters.

        .. note::

            Truncation will be disabled if the truncation value is set to 0.
        """
        cls._truncate__ = int(val)

    @property
    def phase_exp(self):
        """Return the group phase exponent for the Pauli."""
        # Convert internal Pauli encoding to external Pauli encoding
        return pauli_rep.change_pauli_encoding(
            self._phase_exp, self.num_y, output_pauli_encoding=BasePauli.EXTERNAL_PAULI_ENCODING
        )

    @phase_exp.setter
    def phase_exp(self, value):
        # Convert external Pauli encoding to the internal Pauli Encoding
        self._phase_exp[:] = pauli_rep.change_pauli_encoding(
            value,
            self.num_y,
            input_pauli_encoding=BasePauli.EXTERNAL_PAULI_ENCODING,
            output_pauli_encoding=pauli_rep.INTERNAL_PAULI_ENCODING,
            same_type=False,
        )

    @property
    def phase(self):
        """Return the complex phase of the Pauli"""
        return pauli_rep.exp2cpx(self.phase_exp, input_encoding=BasePauli.EXTERNAL_PHASE_ENCODING)

    @property
    def x(self):
        """The x vector for the Pauli."""
        return self.matrix[:, : self.num_qubits][0]

    @x.setter
    def x(self, val):
        self.matrix[:, : self.num_qubits][0] = val

    @property
    def z(self):
        """The z vector for the Pauli."""
        return self.matrix[:, self.num_qubits :][0]

    @z.setter
    def z(self, val):
        self.matrix[:, self.num_qubits :][0] = val

    # ---------------------------------------------------------------------
    # Magic Methods
    # ---------------------------------------------------------------------

    def __repr__(self):
        """Display representation."""
        return f"Pauli('{self.__str__()}')"

    def __str__(self):
        """Print representation."""
        if (
            self._truncate__
            and self.num_qubits > self._truncate__
            and BasePauli.EXTERNAL_SYNTAX == pauli_rep.PRODUCT_SYNTAX
        ):
            front = self[-self._truncate__ :].to_label()
            return front + "..."
        return self.to_label()

    def __array__(self, dtype=None):
        if dtype:
            return np.asarray(self.to_cpx_matrix(), dtype=dtype)
        return self.to_cpx_matrix()

    def __eq__(self, other):
        """Test if two Paulis are equal."""
        if not isinstance(other, BasePauli):
            return False
        return self._eq(other)

    def __len__(self):
        """Return the number of qubits in the Pauli."""
        return self.num_qubits

    def __hash__(self):
        """Make hashable based on string representation."""
        return hash(self.to_label())

    def __setitem__(self, qubits, value):
        """Update the Pauli for a subset of qubits."""
        if not isinstance(value, Pauli):
            value = Pauli(value)
        self._z[0, qubits] = value._z
        self._x[0, qubits] = value._x
        # Add extra phase from new Pauli to current
        self._phase_exp = (self._phase_exp + value._phase_exp) % 4

    def __getitem__(
        self, qubits: Union[int, np.integer, slice, List[int], List[np.integer]]
    ) -> "Pauli":
        """Get no phase Pauli for specific qubits

        Args:
            qubits: index of qubit

        Returns:
            Pauli: Pauli acting on qubit i
        """
        # This is expensive ~6us compared to ~100ns for a Pauli get from a PauliList. Better
        # than the ~25mu that the existing Pauli takes in quantum_info
        # Does allow referencing but limited to single indexing.
        #
        # There are pros and cons associated with how we represent the symplectic
        # matrix. If you represent this matrix as a X and Z matrix then more slicing options
        # are possible here that would allow referening but only to a limited extent.
        #
        # Linear slicing with steps could be down with full vector representations. Linear slicing
        # with steps is chosing a contigous block of qubits to select and steping in that block.
        # for expample say [l:k:s]. To do this you can reshape the full vector from
        # (1, 2*num_qubits) to (2, num_qubits) and then slicing that matrix with [l:k:s]. This
        # gives some matrix of size (2, selected_num_qubits) with the 0th row representing the
        # X part and the 1st row representing the Z part. Yoiu cannot reshape this as it will distroy
        # the refering.
        #
        # This could be used by a second form of the vector (a matrix) would need to be stored and
        # interpreted appropriately.
        #
        # If refercing is not necessary then it is easy to extract the necessary value in what
        # ever form is desired.
        #
        # Suggest that __getitem__ be as general as possible to allow all sorts of possible selctions
        # some of which will and some will not provide referencing. Then have some specific
        # fast functions
        # that check less and can do less.
        #
        # With relative speeds
        #
        # __getitem__  6us (will get slower as more checks and capability added)
        # getitem 6us
        # fast_getitem_int 1.5us
        # fast_getitem_str 1.5us
        # etc
        # If you just reture the tuple (x,z) value then you get
        # getxz
        #
        # One can get faster results if for Paulis we also store a list version
        #
        # vlist_getitem_raw 530ns
        # vlist_getitem_int 580ns
        #
        # These results are fore small dimensions (and using nump arrays). For different situtations
        # different methods may should be used.

        if isinstance(qubits, (int, np.integer)):
            qubits = [qubits]
        return Pauli((self.z[qubits], self.x[qubits]))

    def _getitem(self, i):
        """Get Pauli for qubit i

        Args:
            i (int): index of qubit

        Returns:
            Pauli: Pauli acting on qubit i
        """
        return Pauli(self.matrix[0][i : self.matrix.shape[1] : self.num_qubits])

    def _fast_getitem_str(self, i):
        """Get Pauli for qubit i

        Args:
            i (int): index of qubit

        Returns:
            str: Streing representing the Pauli acting on qubit i,
            (0,0):"I", (1,0):"X", (0,1):"Z", (1,1):"Y"
        """
        return Pauli.pltb_str[(self.matrix[0][i], self.matrix[0][i + self.num_qubits])]

    def _fast_getitem_int(self, i):
        """Get Pauli for qubit i

        Args:
            i (int): index of qubit

        Returns:
            int: Integer representing the Pauli acting on qubit i, (I:0, X:1, Z:2, Y:4)
        """

        return Pauli.pltb_int[(self.matrix[0][i], self.matrix[0][i + self.num_qubits])]

    def _getitem_raw(self, i):
        """Get Pauli symplectic (x,z)-tuple for qubit i

        Args:
            i (int): index of qubit

        Returns:
            (list): Symplectic (x,z)-tuple for the Pauli on qubit i
        """
        return self.matrix[0][i], self.matrix[0][i + self.num_qubits]

    def _vlist_getitem_raw(self, i):
        """Get Pauli symplectic (x,z)-tuple for qubit i. Requires extra storage
        for the list representation of the Pauli

        Args:
            i (int): index of qubit

        Returns:
            [tuple]: Symplectic (x,z)-tuple for the Pauli on qubit i
        """
        return self.vlist[i], self.vlist[i + self.num_qubits]

    def _vlist_getitem_int(self, i):
        """Get Pauli for qubit i using a stored list of the Pauli.

        Args:
            i (int): index of qubit

        Returns:
            int: Integer representing the Pauli acting on qubit i, (I:0, X:1, Z:2, Y:4)
        """
        return Pauli.pltb_int[(self.vlist[i], self.vlist[i + self.num_qubits])]

    # ---------------------------------------------------------------------
    # Insert/Delete Methods
    # ---------------------------------------------------------------------

    def delete(self, qubits: Union[int, List[int]]) -> "Pauli":
        """Return a Pauli with qubits deleted.

        Args:
            qubits (int or list): qubits to delete from Pauli.

        Returns:
            Pauli: the resulting Pauli with the specified qubits removed.

        Raises:
            QiskitError: if ind is out of bounds for the array size or
                         number of qubits.
        """
        if isinstance(qubits, (int, np.integer)):
            qubits = [qubits]
        if max(qubits) > self.num_qubits - 1:
            raise QiskitError(
                f"Qubit index is larger than the number of qubits \
                ({max(qubits)}>{self.num_qubits - 1})."
            )
        if len(qubits) == self.num_qubits:
            raise QiskitError("Cannot delete all qubits of Pauli")
        z = np.delete(self._z, qubits, axis=1)
        x = np.delete(self._x, qubits, axis=1)
        matrix = np.hstack((x, z))
        return Pauli(matrix, phase_exp=self._phase_exp)

    def insert(self, qubits: Union[int, List[int]], value: "Pauli") -> "Pauli":
        """Insert a Pauli at specific qubit value.

        Args:
            qubits (int or list): qubits index to insert at.
            value (Pauli): value to insert.

        Returns:
            Pauli: the resulting Pauli with the entries inserted.

        Raises:
            QiskitError: if the insertion qubits are invalid.
        """
        if not isinstance(value, Pauli):
            value = Pauli(value)

        # Initialize empty operator
        ret_qubits = self.num_qubits + value.num_qubits
        ret = Pauli(np.zeros(2 * ret_qubits, dtype=bool))
        if isinstance(qubits, (int, np.integer)):
            if value.num_qubits == 1:
                qubits = [qubits]
            else:
                qubits = list(range(qubits, qubits + value.num_qubits))
        if len(qubits) != value.num_qubits:
            raise QiskitError(
                f"Number of indices does not match number of qubits for \
                 the inserted Pauli ({len(qubits)}!={value.num_qubits})"
            )
        if max(qubits) > ret.num_qubits - 1:
            raise QiskitError(
                f"Index is too larger for combined Pauli number of qubits \
                 ({max(qubits)}>{ret.num_qubits - 1})."
            )
        # Qubit positions for original op
        self_qubits = [i for i in range(ret.num_qubits) if i not in qubits]
        ret[self_qubits] = self
        ret[qubits] = value
        return ret

    def equiv(self, other):
        """Return True if Pauli's are equivalent up to group phase.

        Args:
            other (Pauli): an operator object.

        Returns:
            bool: True if the Pauli's are equivalent up to group phase.
        """
        if not isinstance(other, Pauli):
            try:
                other = Pauli(other)
            except QiskitError:
                return False
        return np.all(self._z == other._z) and np.all(self._x == other._x)

    @classmethod
    def _from_scalar_op(cls, op):
        return pauli_rep.scalar_op2symplectic(op)

    @classmethod
    def _from_pauli_instruction(cls, instr):
        return pauli_rep.gate2symplectic(instr)

    @classmethod
    def _from_circuit(cls, instr):
        return cls.instrs2symplectic(instr)

    # ---------------------------------------------------------------------
    # Output representation methods: Can also use the pauli_rep methods
    # ---------------------------------------------------------------------

    # Should be deprecated in favour of to_cpx_matrix
    def to_matrix(self, sparse: bool = False) -> np.ndarray:
        """_summary_

        Args:
            sparse (bool, optional): _description_. Defaults to False.

        Returns:
            np.ndarray: _description_
        """
        return self.to_cpx_matrix(sparse=sparse)

    def to_cpx_matrix(self, sparse: bool = False) -> np.ndarray:
        """_summary_

        Args:
            sparse (bool, optional): _description_. Defaults to False.

        Returns:
            np.ndarray: _description_
        """
        return pauli_rep.to_cpx_matrix(self.matrix, self._phase_exp, sparse=sparse)

    def to_instruction(self):
        """Convert to Pauli circuit instruction."""
        from math import pi

        pauli, phase_exp = self.to_label(
            output_pauli_encoding=pauli_rep.DEFAULT_EXTERNAL_PAULI_ENCODING,
            no_phase=True,
            return_phase=True,
            syntax=pauli_rep.PRODUCT_SYNTAX,
            qubit_order="right-to-left",
        )

        if len(pauli) == 1:
            gate = {"I": IGate(), "X": XGate(), "Y": YGate(), "Z": ZGate()}[pauli]
        else:
            gate = PauliGate(pauli)
        if not phase_exp[0]:
            return gate
        # Add global phase
        circuit = QuantumCircuit(self.num_qubits, name=str(self))
        circuit.global_phase = -phase_exp[0] * pi / 2
        circuit.append(gate, range(self.num_qubits))
        return circuit.to_instruction()

    # ---------------------------------------------------------------------
    # BaseOperator methods
    # ---------------------------------------------------------------------

    def compose(
        self,
        other: Union["Pauli", BasePauli],
        qargs: Optional[List[int]] = None,
        front: bool = False,
        inplace: bool = False,
    ) -> "Pauli":
        """Return the operator composition with another Pauli.

        Args:
            other (Pauli): a Pauli object.
            qargs (list or None): Optional, qubits to apply dot product
                                  on (default: None).
            front (bool): If True compose using right operator multiplication,
                          instead of left multiplication [default: False].
            inplace (bool): If True update in-place (default: False).

        Returns:
            Pauli: The composed Pauli.

        Raises:
            QiskitError: if other cannot be converted to an operator, or has
                         incompatible dimensions for specified subsystems.

        .. note::
            Composition (``&``) by default is defined as `left` matrix multiplication for
            matrix operators, while :meth:`dot` is defined as `right` matrix
            multiplication. That is that ``A & B == A.compose(B)`` is equivalent to
            ``B.dot(A)`` when ``A`` and ``B`` are of the same type.

            Setting the ``front=True`` kwarg changes this to `right` matrix
            multiplication and is equivalent to the :meth:`dot` method
            ``A.dot(B) == A.compose(B, front=True)``.
        """
        if qargs is None:
            qargs = getattr(other, "qargs", None)
        if not isinstance(other, Pauli):
            other = Pauli(other)
        return Pauli(super().compose(other, qargs=qargs, front=front, inplace=inplace))

    # pylint: disable=arguments-differ
    def dot(
        self,
        other: Union["Pauli", BasePauli],
        qargs: Optional[List[int]] = None,
        inplace: bool = False,
    ) -> "Pauli":
        """Return the right multiplied operator self * other.

        Args:
            other (Pauli): an operator object.
            qargs (list or None): Optional, qubits to apply dot product
                                  on (default: None).
            inplace (bool): If True update in-place (default: False).

        Returns:
            Pauli: The operator self * other.
        """
        return self.compose(other, qargs=qargs, front=True, inplace=inplace)

    def tensor(self, other) -> "Pauli":
        if not isinstance(other, Pauli):
            other = Pauli(other)
        return Pauli(super().tensor(other))

    def expand(self, other: Union["Pauli", BasePauli]) -> "Pauli":
        if not isinstance(other, Pauli):
            other = Pauli(other)
        return Pauli(super().expand(other))

    def _multiply(self, phase: Union["Pauli", BasePauli]) -> "Pauli":
        return Pauli(super()._multiply(phase))

    def conjugate(self):
        return Pauli(super().conjugate())

    def transpose(self) -> "Pauli":
        return Pauli(super().transpose())

    def adjoint(self) -> "Pauli":
        return Pauli(super().adjoint())

    def inverse(self) -> "Pauli":
        """Return the inverse of the Pauli."""
        return Pauli(super().adjoint())

    # ---------------------------------------------------------------------
    # Utility methods
    # ---------------------------------------------------------------------

    def commutes(self, other, qargs=None):
        """Return True if the Pauli commutes with other.

        Args:
            other (Pauli or PauliList): another Pauli operator.
            qargs (list): qubits to apply dot product on (default: None).

        Returns:
            bool: True if Pauli's commute, False if they anti-commute.
        """
        if qargs is None:
            qargs = getattr(other, "qargs", None)
        if not isinstance(other, BasePauli):
            other = Pauli(other)
        ret = super().commutes(other, qargs=qargs)
        if len(ret) == 1:
            return ret[0]
        return ret

    def anticommutes(self, other, qargs=None):
        """Return True if other Pauli anticommutes with self.

        Args:
            other (Pauli): another Pauli operator.
            qargs (list): qubits to apply dot product on (default: None).

        Returns:
            bool: True if Pauli's anticommute, False if they commute.
        """
        return np.logical_not(self.commutes(other, qargs=qargs))

    def evolve(self, other, qargs=None, frame="h"):
        r"""Heisenberg picture evolution of a Pauli by a Clifford.

        This returns the Pauli :math:`P^\prime = C^\dagger.P.C`.

        By choosing the parameter frame='s', this function returns the Schrödinger evolution of the Pauli
        :math:`P^\prime = C.P.C^\dagger`. This option yields a faster calculation.

        Args:
            other (Pauli or Clifford or QuantumCircuit): The Clifford operator to evolve by.
            qargs (list): a list of qubits to apply the Clifford to.
            frame (string): 'h' for Heisenberg or 's' for Schrödinger framework.

        Returns:
            Pauli: the Pauli :math:`C^\dagger.P.C`.

        Raises:
            QiskitError: if the Clifford number of qubits and qargs don't match.
        """
        if qargs is None:
            qargs = getattr(other, "qargs", None)

        if not isinstance(other, (Pauli, Instruction, QuantumCircuit)):
            # Convert to a Pauli
            other = Pauli(other)

        return Pauli(super().evolve(other, qargs=qargs, frame=frame))

    @staticmethod
    def instrs2symplectic(instr: Union[Instruction, QuantumCircuit]):
        """Convert a Pauli circuit to BasePauli data."""
        # Try and convert single instruction
        if isinstance(instr, (PauliGate, IGate, XGate, YGate, ZGate)):
            return pauli_rep.gate2symplectic(instr)
        if isinstance(instr, Instruction):
            # Convert other instructions to circuit definition
            if instr.definition is None:
                raise QiskitError(f"Cannot apply Instruction: {instr.name}")
            # Convert to circuit
            instr = instr.definition

        # Initialize identity Pauli
        ret = Pauli(np.zeros((1, 2 * instr.num_qubits), dtype=bool), phase_exp=0)
        # Add circuit global phase if specified
        if instr.global_phase:
            ret._phase_exp = pauli_rep.cpx2exp(
                np.exp(1j * float(instr.global_phase)),
                output_encoding=pauli_rep.INTERNAL_PHASE_ENCODING,
            )
        # Recursively apply instructions
        for dinstr, qregs, cregs in instr.data:
            if cregs:
                raise QiskitError(
                    f"Cannot apply instruction with classical registers: {dinstr.name}"
                )
            if not isinstance(dinstr, (Barrier, Delay)):
                next_instr = BasePauli(*Pauli.instrs2symplectic(dinstr))
                if next_instr is not None:
                    qargs = [instr.find_bit(tup)[0] for tup in qregs]
                    ret = ret.compose(next_instr, qargs=qargs)
        return ret.matrix, ret._phase_exp

    # ---------------------------------------------------------------------
    # DEPRECATED methods from old Pauli class
    # ---------------------------------------------------------------------

    @classmethod
    @deprecate_function(
        "Initializing Pauli from `Pauli(label=l)` kwarg is deprecated as of \
         version 0.17.0 and will be removed no earlier than 3 months after \
         the release date. Use `Pauli(l)` instead."
    )
    def _from_label_deprecated(cls, label, qubit_order):
        # Deprecated wrapper of `_from_label` so that a deprecation warning
        # can be displaced during initialization with deprecated kwarg
        return pauli_rep.str2symplectic(label, qubit_order=qubit_order)

    @classmethod
    @deprecate_function(
        "Initializing Pauli from `Pauli(z=z, x=x)` kwargs is deprecated as of \
         version 0.17.0 and will be removed no earlier than 3 months after \
         the release date. Use tuple initialization `Pauli((z, x))` instead."
    )
    def _from_split_array_deprecated(cls, z, x, phase_exp: int = 0):
        # Deprecated wrapper of `_from_array` so that a deprecation warning
        # can be displaced during initialization with deprecated kwarg
        return pauli_rep.from_split_array(x, z, phase_exp)

    @staticmethod
    def _make_np_bool(arr):
        if not isinstance(arr, (list, np.ndarray, tuple)):
            arr = [arr]
        arr = np.asarray(arr).astype(bool)
        return arr

    @staticmethod
    @deprecate_function(
        "`from_label` is deprecated and will be removed no earlier than \
         3 months after the release date. Use Pauli(label) instead."
    )
    def from_label(label):
        """DEPRECATED: Construct a Pauli from a string label.

        This function is deprecated use ``Pauli(label)`` instead.

        Args:
            label (str): Pauli string label.

        Returns:
            Pauli: the constructed Pauli.

        Raises:
            QiskitError: If the input list is empty or contains invalid
            Pauli strings.
        """
        if isinstance(label, tuple):
            # Legacy usage from aqua
            label = "".join(label)
        return Pauli(label)

    @staticmethod
    @deprecate_function(
        "sgn_prod is deprecated and will be removed no earlier than \
         3 months after the release date. Use `dot` instead."
    )  # pylint: disable=invalid-name
    def sgn_prod(p1, p2):  # pylint: disable=invalid-name
        r"""
        DEPRECATED: Multiply two Paulis and track the phase.

        This function is deprecated. The Pauli class now handles full
        Pauli group multiplication using :meth:`compose` or :meth:`dot`.

        $P_3 = P_1 \otimes P_2$: X*Y

        Args:
            p1 (Pauli): pauli 1
            p2 (Pauli): pauli 2

        Returns:
            Pauli: the multiplied pauli (without phase)
            complex: the sign of the multiplication, 1, -1, 1j or -1j
        """
        pauli = p1.dot(p2)
        return pauli[:], (-1j) ** pauli.phase

    @deprecate_function(
        "`to_spmatrix` is deprecated and will be removed no earlier than \
         3 months after the release date. Use `to_matrix(sparse=True)` instead."
    )
    def to_spmatrix(self):
        r"""
        DEPRECATED Convert Pauli to a sparse matrix representation (CSR format).

        This function is deprecated. Use :meth:`to_matrix` with kwarg
        ``sparse=True`` instead.

        Returns:
            scipy.sparse.csr_matrix: a sparse matrix with CSR format that
            represents the pauli.
        """
        return pauli_rep.to_cpx_matrix(self.matrix, self._phase_exp, sparse=True)

    @deprecate_function(
        "`kron` is deprecated and will be removed no earlier than \
        3 months after the release date of Qiskit Terra 0.17.0. \
        Use `expand` instead, but note this does not change \
        the operator in-place."
    )
    def kron(self, other):
        r"""DEPRECATED: Kronecker product of two paulis.

        This function is deprecated. Use :meth:`expand` instead.

        Order is $P_2 (other) \otimes P_1 (self)$

        Args:
            other (Pauli): P2

        Returns:
            Pauli: self
        """
        pauli = self.expand(other)
        self.matrix = pauli.matrix
        self._phase_exp = pauli._phase_exp
        self._op_shape = self._op_shape.expand(other._op_shape)
        return self

    @deprecate_function(
        "`update_z` is deprecated and will be removed no earlier than \
        3 months after the release date. Use `Pauli.z = val` or \
        `Pauli.z[indices] = val` instead."
    )
    def update_z(self, z, indices=None):
        """
        DEPRECATED: Update partial or entire z.

        This function is deprecated. Use the setter for :attr:`z` instead.

        Args:
            z (numpy.ndarray or list): to-be-updated z
            indices (numpy.ndarray or list or optional): to-be-updated qubit indices

        Returns:
            Pauli: self

        Raises:
            QiskitError: when updating whole z, the number of qubits must be the same.
        """
        phase_exp = self._phase_exp
        z = self._make_np_bool(z)
        if indices is None:
            if len(self.z) != len(z):
                raise QiskitError(
                    "During updating whole z, you can not change the number of qubits."
                )
            self.z = z
        else:
            if not isinstance(indices, list) and not isinstance(indices, np.ndarray):
                indices = [indices]
            for p, idx in enumerate(indices):
                self.z[idx] = z[p]
        self._phase_exp = phase_exp
        return self

    @deprecate_function(
        "`update_z` is deprecated and will be removed no earlier than \
         3 months after the release date. Use `Pauli.x = val` or \
         `Pauli.x[indices] = val` instead."
    )
    def update_x(self, x, indices=None):
        """
        DEPRECATED: Update partial or entire x.

        This function is deprecated. Use the setter for :attr:`x` instead.

        Args:
            x (numpy.ndarray or list): to-be-updated x
            indices (numpy.ndarray or list or optional): to-be-updated qubit indices

        Returns:
            Pauli: self

        Raises:
            QiskitError: when updating whole x, the number of qubits must be the same.
        """
        phase_exp = self._phase_exp
        x = self._make_np_bool(x)
        if indices is None:
            if len(self.x) != len(x):
                raise QiskitError(
                    "During updating whole x, you can not change the number of qubits."
                )
            self.x = x
        else:
            if not isinstance(indices, list) and not isinstance(indices, np.ndarray):
                indices = [indices]
            for p, idx in enumerate(indices):
                self.x[idx] = x[p]
        self._phase_exp = phase_exp
        return self

    @deprecate_function(
        "`insert_paulis` is deprecated and will be removed no earlier than \
         3 months after the release date. For similar functionality use \
         `Pauli.insert` instead."
    )
    def insert_paulis(self, indices=None, paulis=None, pauli_labels=None):
        """
        DEPRECATED: Insert or append pauli to the targeted indices.

        This function is deprecated. Similar functionality can be obtained
        using the :meth:`insert` method.

        If indices is None, it means append at the end.

        Args:
            indices (list[int]): the qubit indices to be inserted
            paulis (Pauli): the to-be-inserted or appended pauli
            pauli_labels (list[str]): the to-be-inserted or appended pauli label

        Note:
            the indices refers to the location of original paulis,
            e.g. if indices = [0, 2], pauli_labels = ['Z', 'I'] and original pauli = 'ZYXI'
            the pauli will be updated to ZY'I'XI'Z'
            'Z' and 'I' are inserted before the qubit at 0 and 2.

        Returns:
            Pauli: self

        Raises:
            QiskitError: provide both `paulis` and `pauli_labels` at the same time
        """
        if pauli_labels is not None:
            if paulis is not None:
                raise QiskitError("Please only provide either `paulis` or `pauli_labels`")
            if isinstance(pauli_labels, str):
                pauli_labels = list(pauli_labels)
            # since pauli label is in reversed order.
            label = "".join(pauli_labels[::-1])
            paulis = self.from_label(label)

        # Insert and update self
        if indices is None:  # append
            z = np.concatenate((self.z, paulis.z))
            x = np.concatenate((self.x, paulis.x))
        else:
            if not isinstance(indices, list):
                indices = [indices]
            z = np.insert(self.z, indices, paulis.z)
            x = np.insert(self.x, indices, paulis.x)
        matrix = np.hstack((x, z))
        pauli = Pauli(matrix, phase_exp=self._phase_exp + paulis._phase_exp)
        self.z = pauli.z
        self.x = pauli.x
        self._phase_exp = pauli._phase_exp
        self._op_shape = pauli._op_shape
        return self

    @deprecate_function(
        "`append_paulis` is deprecated and will be removed no earlier than \
         3 months after the release date. Use `Pauli.expand` instead."
    )
    def append_paulis(self, paulis=None, pauli_labels=None):
        """
        DEPRECATED: Append pauli at the end.

        Args:
            paulis (Pauli): the to-be-inserted or appended pauli
            pauli_labels (list[str]): the to-be-inserted or appended pauli label

        Returns:
            Pauli: self
        """
        return self.insert_paulis(None, paulis=paulis, pauli_labels=pauli_labels)

    @deprecate_function(
        "`append_paulis` is deprecated and will be removed no earlier than \
         3 months after the release date. For equivalent functionality \
         use `Pauli.delete` instead."
    )
    def delete_qubits(self, indices):
        """
        DEPRECATED: Delete pauli at the indices.

        This function is deprecated. Equivalent functionality can be obtained
        using the :meth:`delete` method.

        Args:
            indices(list[int]): the indices of to-be-deleted paulis

        Returns:
            Pauli: self
        """
        pauli = self.delete(indices)
        self._z = pauli._z  # pylint: disable=invalid-name
        self._x = pauli._x  # pylint: disable=invalid-name
        self._phase_exp = pauli._phase_exp
        self._op_shape = pauli._op_shape
        return self

    @classmethod
    @deprecate_function(
        "`pauli_single` is deprecated and will be removed no earlier than \
         3 months after the release date."
    )
    def pauli_single(cls, num_qubits, index, pauli_label):
        """
        DEPRECATED: Generate single qubit pauli at index with pauli_label with length num_qubits.

        Args:
            num_qubits (int): the length of pauli
            index (int): the qubit index to insert the single qubit
            pauli_label (str): pauli

        Returns:
            Pauli: single qubit pauli
        """
        tmp = Pauli(pauli_label)
        ret = Pauli((np.zeros(num_qubits, dtype=bool), np.zeros(num_qubits, dtype=bool)))
        ret.x[index] = tmp.x[0]
        ret.z[index] = tmp.z[0]
        ret.phase = tmp.phase
        return ret

    @classmethod
    @deprecate_function(
        "`random` is deprecated and will be removed no earlier than \
         3 months after the release date. \
         Use `qiskit.quantum_info.random_pauli` instead"
    )
    def random(cls, num_qubits, seed=None):
        """DEPRECATED: Return a random Pauli on number of qubits.

        This function is deprecated use
        :func:`~qiskit.quantum_info.random_pauli` instead.

        Args:
            num_qubits (int): the number of qubits
            seed (int): Optional. To set a random seed.
        Returns:
            Pauli: the random pauli
        """
        # pylint: disable=cyclic-import
        from qiskit.quantum_info.operators.symplectic.random import random_pauli

        return random_pauli(num_qubits, group_phase=False, seed=seed)


# Update docstrings for API docs
generate_apidocs(Pauli)
