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
"""
Optimized list of Pauli operators
"""

"""
Change Log:

Jul 27, 2021 Drew Vandeth

    Added PauliRep and a new base class to handle different representations and
    conversion utilities

    Moved _to_label, _to_matrix, _from_array, _phase_from_label to PauliRep class.

August 16, 2021 Drew Vandeth

    Improved commutes performance by simply returning
    
    np.logical_not(np.sum((x1&other._z)^(z1&other._x),axis=1) % 2) after x1 and
    z1 are set. Since is is a direct description of the symplectic product it
    is still readable and faster at the same time.

December 7th, 2021 Drew Vandeth

Notes:

    Inplace versus not inplace. The variuos methods need to be gone over to be
    consistent with inplace versus non-inplace options or behaviours. For operations
    such a conjugate, or transpose, or _multiply, etc. we should consider splitting the
    two variations up: one for inplace and one for not inplace. So for example having
    both a iconjugate for inplace conjugation and conjugate for non-inplace conjugation.
    Or something like that.

    I like the idea of having a <name> method that has the safety checks etc that calls
    the <_name> method that does not perform checks or input conversions. The second one
    can be used for speed up when appropriate. Perhaps a better name than <_name> should be
    provided such as <f_name> for fast_name variant.

"""
# pylint: disable=invalid-name

import copy
import sys

import numpy as np
import qiskit_qec.utils.pauli_rep as pauli_rep
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.barrier import Barrier
from qiskit.exceptions import QiskitError
from qiskit.quantum_info.operators.base_operator import BaseOperator
from qiskit.quantum_info.operators.mixins import AdjointMixin, MultiplyMixin
from qiskit_qec.structures.symplectic_matrix_dense_numpy_old import *
from qiskit_qec.structures.symplectic_matrix_dense_numpy_old import (
    SymplecticMatrixBase,
    SymplecticMatrixDenseNumPy,
)


class PaulisBase(BaseOperator, AdjointMixin, MultiplyMixin):
    r"""Base representation of a list of N-qubit Paulis


    Base class for Pauli and PauliList.
    """
    DENSENUMPY = SymplecticMatrixBase.DENSENUMPY

    external_symp_rep_format = pauli_rep.DEFAULT_EXTERNAL_SYMP_FORMAT
    external_phase_rep_format = pauli_rep.DEFAULT_EXTERNAL_PHASE_REP_FORMAT
    external_pauli_rep_format = external_phase_rep_format + external_symp_rep_format

    # -------------------------------------------------------------------------------
    # Initilization Method
    # -------------------------------------------------------------------------------

    def __init__(self, xz_matrix, phase_exponent=None, stype=DENSENUMPY):
        """Initialize the PaulisBase.

        This is an array of k N-qubit Paulis represented in the '-iZX'
        frmt:
        P = (-i)^phase Z^z X^x.
        This representation is chosen to simplify various operations and to be consistent
        with other aspects within Qiskit. This frmt only changes
        the interpretation of the phase which is recorded as an exponent of -i.

        Args:
            xz_matrix: input symplectic matrix
            phase_exponent (np.ndarray): input phase exponent vector in '-i' representation frmt
            stype (int):  Specifies which interal representation type to use for symplectic matrices. The
            possibles types are returned by SymplecticMatrixBase.stypes(). By default the Dense NumPy
            sympletic matrix module is used. If using a different symplectic matrix module from the default
            then this module must be imported first.
        """

        # Note: The code below to check and then set self.symplectic_class should probably be
        # replaced with something better like defined types etc. But this works and is simple
        # enough for now.

        # Check if the requested symplectic matrix module is supported
        if stype not in SymplecticMatrixBase.stypes():
            print("Error: Unknown symplectic matric type")
            exit

        self.symplectic_class_type = stype
        self.symplectic_class_name = SymplecticMatrixBase.classnames(stype)

        try:
            self.symplectic_class = getattr(sys.modules[__name__], self.symplectic_class_name)
        except:
            print(
                "Error: To use a non-default symplectic representation first import the needed module"
            )
            exit

        self._matrix = self.symplectic_class(xz_matrix)

        self._phase_exponent = np.asarray(phase_exponent)

        super().__init__(num_qubits=self._matrix.num_qubits)

    @classmethod
    def _set_formats(cls, *, phase_format=None, symp_format=None):
        """Update the external phase and symplectic represenations

        Calling with no parameters will reset to the default external
        representations
        Args:
            phase_format (str) : Phase frmt string from pauli_rep.PHASE_REP_FORMATS. Default = '-i'
            symp_format (str): Symplectic frmt string from pauli_rep.SYMP_REP_FORMATS. Default = 'YZX'
        Raiss:
            QiskitError: If formats are not implemented
        """

        if phase_format is not None:
            if phase_format in pauli_rep.PHASE_REP_FORMATS:
                cls.external_phase_rep_format = phase_format
            else:
                raise QiskitError("Invalid phase frmt")

        if symp_format is not None:
            if symp_format in pauli_rep.SYMP_REP_FORMATS:
                cls.external_symp_rep_format = symp_format
            else:
                raise QiskitError("Invalid symplectic frmt")

        if phase_format is None:
            phase_format = pauli_rep.DEFAULT_EXTERNAL_PHASE_REP_FORMAT
        if symp_format is None:
            symp_format = pauli_rep.DEFAULT_EXTERNAL_SYMP_FORMAT

        cls.external_phase_rep_format = phase_format
        cls.external_symp_rep_format = symp_format
        cls.external_pauli_rep_format = phase_format + symp_format

    # ---------------------------------------------------------------------
    # Property methods
    # ---------------------------------------------------------------------

    @staticmethod
    def set_formats(*, phase_format=None, symp_format=None):
        """Update the external phase and symplectic represenations

        Calling with no parameters will reset to the default external
        representations
        Args:
            phase_format (str) : Phase frmt string from PHASE_REP_FORMATS. Default = '-i'
            symp_format (str): Symplectic frmt string from SYMP_REP_FORMATS. Default = 'YZX'
        Raiss:
            QiskitError: If formats are not implemented
        """
        PaulisBase._set_formats(phase_format=phase_format, symp_format=symp_format)

    # External formats

    @staticmethod
    def external_pauli_format():
        """Display the external representation frmt for Pauli's
        Returns:
            str: Representation frmt for Pauli operator
        """
        return PaulisBase.external_pauli_rep_format

    @staticmethod
    def set_external_pauli_format(format):
        """[summary]
        Args:
            format ([type]): [description]
        Raises:
            QiskitError: [description]
        """
        if format in pauli_rep.pauli_formats():
            phase_format, symp_format = pauli_rep._split_rep(format)

            PaulisBase.set_formats(phase_format=phase_format, symp_format=symp_format)
        else:
            raise QiskitError("Invalid Pauli representation frmt or unsupported frmt")

    @staticmethod
    def external_phase_format():
        """Display the external phase representation frmt for Pauli's
        Returns:
            str: Phase representation frmt for Pauli operator
        """
        return PaulisBase.external_phase_rep_format

    @staticmethod
    def set_external_phase_format(phase_format=None):
        """[summary]
        Args:
            phase_format ([type]): [description]
        """
        PaulisBase.set_formats(phase_format=phase_format)

    @staticmethod
    def external_symp_format():
        """Display the external symplectic representation frmt for Pauli's
        Returns:
            str: Symplectic representation frmt for Pauli operator
        """
        return PaulisBase.external_symp_rep_format

    @staticmethod
    def set_external_symp_format(symp_format=None):
        """[summary]
        Args:
            symp_format ([type], optional): [description]. Defaults to None.
        """
        PaulisBase.set_formats(symp_format=symp_format)

    # Internal formats
    @staticmethod
    def internal_phase_format():
        """[summary]
        Returns:
            [type]: [description]
        """
        return pauli_rep.INTERNAL_PHASE_REP_FORMAT

    @staticmethod
    def internal_symp_format():
        """[summary]
        Returns:
            [type]: [description]
        """
        return pauli_rep.INTERNAL_SYMP_REP_FORMAT

    @staticmethod
    def internal_pauli_format():
        """[summary]
        Returns:
            [type]: [description]
        """
        return pauli_rep.INTERNAL_PAULI_REP_FORMAT

    @property
    def num_paulis(self):
        """Number of Paulis represented in list/symplectic matrix

        Returns:
            int: Number of Paulis
        """
        return self._matrix.num_paulis

    # Note this replaces the _count_y etc
    @property
    def num_y(self):
        """Count the number of Y Pauli's"""
        return self._matrix.num_y()

    def copy(self):
        """Make a deep copy of current operator."""
        # Deepcopy has terrible performance on objects with Numpy arrays
        # attributes so we make a shallow copy and then manually copy the
        # Numpy arrays to efficiently mimic a deepcopy
        ret = copy.copy(self)
        ret._matrix = self._matrix.copy()
        ret._phase_exponent = self._phase_exponent.copy()
        return ret

    # ---------------------------------------------------------------------
    # BaseOperator methods
    # ---------------------------------------------------------------------

    def tensor(self, other):
        return self._tensor(self, other)

    def expand(self, other):
        return self._tensor(other, self)

    @classmethod
    def _tensor(cls, a, b):
        """Generalized tensor product of Paulis (gtensor)

        Args:
            a (PauliList or Pauli): List of Pauli operators
            b (PauliList of Pauli): List of Pauli operators

        if a=[A_1,A_2,...,A_k] and b=[B_1,B_2,...,B_v] then

        a gtensor b = [A_1 tensor B_1, A_1 tensor B_2, ..., A_k tensor B_v]

        Returns:
            [PaulisBase]: a gtensor b
        """

        SympMatrix = a.symplectic_class

        a_matrix = a._matrix

        if SympMatrix.__name__ != b.symplectic_class.__name__:
            c_matrix = SympMatrix.convert(b._matrix)
        else:
            c_matrix = b._matrix

        x1 = SympMatrix.istack(a_matrix.x, c_matrix.num_paulis, True)
        x2 = SympMatrix.istack(c_matrix.x, a_matrix.num_paulis, False)
        z1 = SympMatrix.istack(a_matrix.z, c_matrix.num_paulis, True)
        z2 = SympMatrix.istack(c_matrix.z, a_matrix.num_paulis, False)
        phase1 = (
            np.vstack(c_matrix.num_paulis * [a._phase_exponent])
            .transpose(1, 0)
            .reshape(a_matrix.num_paulis * c_matrix.num_paulis)
        )
        phase2 = SympMatrix.istack(b._phase_exponent, a_matrix.num_paulis)

        xz_matrix = SympMatrix.hstack((x2, x1, z2, z1))
        phase_exponent = np.mod(phase1 + phase2, 4)

        return PaulisBase(xz_matrix=xz_matrix, phase_exponent=phase_exponent)

    # pylint: disable=arguments-differ
    def compose(self, other, qargs=None, front=False, inplace=False):
        """Return the composition of Paulis lists

        To be consistent with other compose functions in Qiskit, composition is defined via 
        left matrix multiplication. That is

        A.compose(B) = B.A =B.dot(A) = A.compose(B, front=False)
        
        and so B is applied after A. Likewise

        A.compose(B, front=True) = A.B = A.dot(B)

        That is B is applied first or at the front.

        This method does compose coordinate wise (which is different from the PauliTable complose
        whish should be corrected at some point). That is 

        [A_1,A_2,...,A_k].compose([B_1,B_2,...,B_k]) = [A_1.compose(B_1),...,A_k.compose(B_k)]

        or

        [A].compose([B_1,B_2,...,B_k])) = [A.compose(B_1),...,A.compose(B_k)]

        etc.

        Args:
            front (bool): (default: False)
            qargs (list or None): Optional, qubits to apply compose on
                                  on (default: None->All).
            inplace (bool): If True update in-place (default: False).

        Returns:
            {cls}: The operator self.compose(other, qargs=qargs, front=front, inplace=inplace)

        Raises:
            QiskitError: if number of qubits of other does not match qargs.
        """.format(
            cls=type(self).__name__
        )
        # Validation
        if qargs is None and other.num_qubits != self.num_qubits:
            raise QiskitError(f"other {type(self).__name__} must be on the same number of qubits.")

        if qargs and other.num_qubits != len(qargs):
            raise QiskitError(
                f"Number of qubits of the other {type(self).__name__} does not match qargs."
            )

        if other.num_paulis not in [1, self.num_paulis]:
            raise QiskitError(
                "Incompatible PaulisBases. Second list must "
                "either have 1 or the same number of Paulis."
            )

        SympMatrix = self.symplectic_class

        # Compute phase shift
        if qargs is not None:
            x1, z1 = self._matrix.x[:, qargs], self._matrix.z[:, qargs]
        else:
            x1, z1 = self._matrix.x, self._matrix.z

        if SympMatrix.__name__ != other.symplectic_class.__name__:
            new_other = PaulisBase(
                xz_matrix=SympMatrix.convert(other._matrix), phase_exponent=other._phase_exponent
            )
        else:
            new_other = other

        x2, z2 = new_other._matrix.x, new_other._matrix.z

        # Get phase shift from reordering X and Z
        phase_exponent = self._phase_exponent + new_other._phase_exponent
        if front:
            phase_exponent += 2 * SympMatrix.sum(SympMatrix.logical_and(x1, z2), axis=1)
        else:
            phase_exponent += 2 * SympMatrix.sum(SympMatrix.logical_and(z1, x2), axis=1)

        phase_exponent = np.mod(phase_exponent, 4)

        # Since we already have the seprated objects just use them and recombine after
        x = SympMatrix.logical_xor(x1, x2)
        z = SympMatrix.logical_xor(z1, z2)
        xz_matrix = SympMatrix.combine(x, z)

        if qargs is None:
            if not inplace:
                return PaulisBase(xz_matrix=xz_matrix, phase_exponent=phase_exponent)
            # Inplace update
            self._matrix = xz_matrix
            self._phase = phase_exponent
            return self

        # Qargs update
        ret = self if inplace else self.copy()
        ret._matrix.x[:, qargs] = xz_matrix.x
        ret._matrix.z[:, qargs] = xz_matrix.z
        # TODO(Drew Vandeth): Possible incorrect calculation on next line. Takes entire new phase and not a modified phase
        # Need to come back and correct if needed.
        ret._phase_exponent = phase_exponent
        return ret

    def _multiply(self, phase, roundit=True):
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
        phase_exponent = pauli_rep.phase_to_exponent(
            phase, output_phase_format=pauli_rep.INTERNAL_PHASE_REP_FORMAT, roundit=roundit
        )
        return PaulisBase(self._matrix, np.mod(self._phase_exponent + phase_exponent, 4))

    def conjugate(self, inplace=False):
        """Return the conjugate of each Pauli in the list.

        Converts internal '-i' frmt into 'i' frmt

        Args:
            inplace (boolean) : Default: False, If True will modify inplace

        Returns:
            {cls} : a new {cls} which has phases conjugates (if replace=False) or will change
            the phase of the clasing instance if replace=True
        """

        # Note: The original checked if the input set of Paulis would be invariant under
        # conjugation and return itself but it is changed it would return a new Pauli
        # that referenced the old Pauli with a new phase. The problem I have with this
        # can be seen in the following example:
        #
        # A=PauliList(['XX','ZZ'])
        # A_conjugate=A.conjugate()
        # A[0]='YI'
        # A_conjugate==A.conjugate() -> False
        #
        # Then A_conjugate does not equal to A.conjugate() by A_conjugate no longer equal the
        # orginal A.conjugate.
        # For this reason the code below removes this behaviour

        new_phase_exp = pauli_rep.convert_phase_exp(self._phase_exponent, output_format="i")
        if inplace == False:
            return PaulisBase(self._matrix, new_phase_exp).copy()
        elif new_phase_exp == self._phase_exponent:
            return self
        else:
            self._phase_exponent = new_phase_exp
            return self

    def transpose(self, inplace=False):
        """Return the transpose of each Pauli in the list.
        Transpose sets X->X, Z->Z, and Y -> -Y.

        """
        new_phase_exp = np.mod(self._phase_exponent + 2 * (self.num_y % 2), 4)
        if inplace == False:
            return PaulisBase(self._matrix, new_phase_exp).copy()
        elif new_phase_exp == self._phase_exponent:
            return self
        else:
            self._phase_exponent = new_phase_exp
            return self

    def _commutes(self, other, qargs=None):
        if qargs is not None:
            inds = list(qargs)
            x1, z1 = self.x[:, inds], self.z[:, inds]
        else:
            x1, z1 = self.x, self.z

        SympMatrix = self.symplectic_class
        return np.logical_not(SympMatrix.sum((x1 & other.z) ^ (z1 & other.x), axis=1) % 2)

    def commutes(self, other, qargs=None):
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
                "Number of qubits of other Pauli does not match number of "
                "qargs ({} != {}).".format(other.num_qubits, len(qargs))
            )
        if qargs is None and self.num_qubits != other.num_qubits:
            raise QiskitError(
                "Number of qubits of other Pauli does not match the current "
                "Pauli ({} != {}).".format(other.num_qubits, self.num_qubits)
            )
        return self._commutes(other, qargs=qargs)

    def evolve(self, other, qargs=None):
        r"""Heisenberg picture evolution of a Pauli by a Clifford.

        This returns the Pauli :math:`P^\prime = C^\dagger.P.C`.

        Args:
            other (PaulisBase or QuantumCircuit): The Clifford circuit to evolve by.
            qargs (list): a list of qubits to apply the Clifford to.

        Returns:
            PaulisBase: the Pauli :math:`C^\dagger.P.C`.

        Raises:
            QiskitError: if the Clifford number of qubits and qargs don't match.
        """
        # Check dimension
        if qargs is not None and len(qargs) != other.num_qubits:
            raise QiskitError(
                "Incorrect number of qubits for Clifford circuit ({} != {}).".format(
                    other.num_qubits, len(qargs)
                )
            )
        if qargs is None and self.num_qubits != other.num_qubits:
            raise QiskitError(
                "Incorrect number of qubits for Clifford circuit ({} != {}).".format(
                    other.num_qubits, self.num_qubits
                )
            )

        # Evolve via Pauli
        if isinstance(other, PaulisBase):
            ret = self.compose(other.adjoint(), qargs=qargs)
            ret = ret.compose(other, front=True, qargs=qargs)
            return ret

        # pylint: disable=cyclic-import
        from qiskit.quantum_info.operators.symplectic.clifford import Clifford

        # Convert Clifford to quantum circuits
        if isinstance(other, Clifford):
            return self._evolve_clifford(other, qargs=qargs)

        # Otherwise evolve by the inverse circuit to compute C^dg.P.C
        return self.copy()._append_circuit(other.inverse(), qargs=qargs)

    # To be refactored
    def _evolve_clifford(self, other, qargs=None):
        """Heisenberg picture evolution of a Pauli by a Clifford."""
        if qargs is None:
            idx = slice(None)
            num_act = self.num_qubits
        else:
            idx = list(qargs)
            num_act = len(idx)

        # Set return to I on qargs
        ret = self.copy()
        ret.x[:, idx] = False
        ret.z[:, idx] = False

        # pylint: disable=cyclic-import
        from qiskit.quantum_info.operators.symplectic.pauli import Pauli
        from qiskit.quantum_info.operators.symplectic.pauli_list import PauliList

        # Get action of Pauli's from Clifford
        adj = other.adjoint()
        pauli_list = []
        for z in self._z[:, idx]:
            pauli = Pauli("I" * num_act)
            for row in adj.stabilizer[z]:
                pauli.compose(Pauli((row.Z[0], row.X[0], 2 * row.phase[0])), inplace=True)
            pauli_list.append(pauli)
        ret.dot(PauliList(pauli_list), qargs=qargs, inplace=True)

        pauli_list = []
        for x in self._x[:, idx]:
            pauli = Pauli("I" * num_act)
            for row in adj.destabilizer[x]:
                pauli.compose(Pauli((row.Z[0], row.X[0], 2 * row.phase[0])), inplace=True)
            pauli_list.append(pauli)
        ret.dot(PauliList(pauli_list), qargs=qargs, inplace=True)
        return ret

    def _eq(self, other):
        """Entrywise comparison of Pauli equality."""
        return (
            self.num_qubits == other.num_qubits
            and np.all(np.mod(self._phase, 4) == np.mod(other._phase, 4))
            and np.all(self._z == other._z)
            and np.all(self._x == other._x)
        )

    # ---------------------------------------------------------------------
    # Helper Methods
    # ---------------------------------------------------------------------
    # To be refactored
    def __imul__(self, other):
        return self.compose(other, front=True, inplace=True)

    # To be refactored
    def __neg__(self):
        ret = copy.copy(self)
        ret._phase = np.mod(self._phase + 2, 4)
        return ret

    # Replaced by specific methoids istack
    @staticmethod
    def _stack(array, size, vertical=True):
        """Stack array."""
        if size == 1:
            return array
        if vertical:
            return np.vstack(size * [array]).reshape((size * len(array),) + array.shape[1:])
        return np.hstack(size * [array]).reshape((size * len(array),) + array.shape[1:])

    # To be refactored
    def _append_circuit(self, circuit, qargs=None):
        """Update PaulisBase inplace by applying a Clifford circuit.

        Args:
            circuit (QuantumCircuit or Instruction): the gate or composite gate to apply.
            qargs (list or None): The qubits to apply gate to.

        Returns:
            PaulisBase: the updated Pauli.

        Raises:
            QiskitError: if input gate cannot be decomposed into Clifford gates.

        #TODO: Update phase conversion to general method
        """
        if isinstance(circuit, Barrier):
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
                "{} instruction definition is {}; expected QuantumCircuit".format(
                    gate.name, type(gate.definition)
                )
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
        self._phase %= 4
        return self


# ---------------------------------------------------------------------
# Evolution by Clifford gates
# ---------------------------------------------------------------------

# To be refactored
def _evolve_h(base_pauli, qubit):
    """Update P -> H.P.H"""
    x = base_pauli.x[:, qubit].copy()
    z = base_pauli.z[:, qubit].copy()
    base_pauli.x[:, qubit] = z
    base_pauli.z[:, qubit] = x
    base_pauli.phase_exponent += 2 * np.logical_and(x, z).T
    return base_pauli


# To be refactored
def _evolve_s(base_pauli, qubit):
    """Update P -> S.P.Sdg"""
    x = base_pauli._x[:, qubit]
    base_pauli._z[:, qubit] ^= x
    base_pauli._phase += x.T
    return base_pauli


# To be refactored
def _evolve_sdg(base_pauli, qubit):
    """Update P -> Sdg.P.S"""
    x = base_pauli._x[:, qubit]
    base_pauli._z[:, qubit] ^= x
    base_pauli._phase -= x.T
    return base_pauli


# To be refactored
# pylint: disable=unused-argument
def _evolve_i(base_pauli, qubit):
    """Update P -> P"""
    return base_pauli


# To be refactored
def _evolve_x(base_pauli, qubit):
    """Update P -> X.P.X"""
    base_pauli._phase += 2 * base_pauli._z[:, qubit].T
    return base_pauli


# To be refactored
def _evolve_y(base_pauli, qubit):
    """Update P -> Y.P.Y"""
    base_pauli._phase += 2 * base_pauli._x[:, qubit].T + 2 * base_pauli._z[:, qubit].T
    return base_pauli


# To be refactored
def _evolve_z(base_pauli, qubit):
    """Update P -> Z.P.Z"""
    base_pauli._phase += 2 * base_pauli._x[:, qubit].T
    return base_pauli


# To be refactored
def _evolve_cx(base_pauli, qctrl, qtrgt):
    """Update P -> CX.P.CX"""
    base_pauli._x[:, qtrgt] ^= base_pauli._x[:, qctrl]
    base_pauli._z[:, qctrl] ^= base_pauli._z[:, qtrgt]
    return base_pauli


# To be refactored
def _evolve_cz(base_pauli, q1, q2):
    """Update P -> CZ.P.CZ"""
    x1 = base_pauli._x[:, q1].copy()
    x2 = base_pauli._x[:, q2].copy()
    base_pauli._z[:, q1] ^= x2
    base_pauli._z[:, q2] ^= x1
    base_pauli._phase += 2 * np.logical_and(x1, x2).T
    return base_pauli


# To be refactored
def _evolve_cy(base_pauli, qctrl, qtrgt):
    """Update P -> CY.P.CY"""
    x1 = base_pauli._x[:, qctrl].copy()
    x2 = base_pauli._x[:, qtrgt].copy()
    z2 = base_pauli._z[:, qtrgt].copy()
    base_pauli._x[:, qtrgt] ^= x1
    base_pauli._z[:, qtrgt] ^= x1
    base_pauli._z[:, qctrl] ^= np.logical_xor(x2, z2)
    base_pauli._phase += x1 + 2 * np.logical_and(x1, x2).T
    return base_pauli


# To be refactored
def _evolve_swap(base_pauli, q1, q2):
    """Update P -> SWAP.P.SWAP"""
    x1 = base_pauli._x[:, q1].copy()
    z1 = base_pauli._z[:, q1].copy()
    base_pauli._x[:, q1] = base_pauli._x[:, q2]
    base_pauli._z[:, q1] = base_pauli._z[:, q2]
    base_pauli._x[:, q2] = x1
    base_pauli._z[:, q2] = z1
    return base_pauli
