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
# pylint: disable=invalid-name

import copy
import numpy as np

from qiskit.exceptions import QiskitError
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.barrier import Barrier
from qiskit.quantum_info.operators.base_operator import BaseOperator
from qiskit.quantum_info.operators.mixins import AdjointMixin, MultiplyMixin
from qiskit.quantum_info.operators.symplectic.pauli_rep import PauliRep

class BasePauli(BaseOperator, AdjointMixin, MultiplyMixin, PauliRep):
    r"""Phase symplectic representation of a list of N-qubit Paulis.

    An array of k N-qubit Paulis represented in the '-iZX'
        format:

    Base class for Pauli and PauliList.
    """

    def __init__(self, z, x, phase):
        """Initialize the BasePauli.

        An array of k N-qubit Paulis represented in the '-iZX'
        format:
        math::
            P = (-i)^phase Z^z X^x.

        Note that the <phase> variable is not the actual phase of the Pauli 
        operator but the exponent of the when the true phase is encoded as
        a power of :math:`-i`.

        Args:
            z (np.ndarray): input Z matrix.
            x (np.ndarray): input X matrix.
            phase (np.ndarray): input phase exponent vector.
        """
        self._z = z
        self._x = x
        self._phase = phase
        self._num_paulis, num_qubits = self._z.shape
        super().__init__(num_qubits=num_qubits)

    def copy(self):
        """Make a deep copy of current operator."""
        # Deepcopy has terrible performance on objects with Numpy arrays
        # attributes so we make a shallow copy and then manually copy the
        # Numpy arrays to efficiently mimic a deepcopy
        ret = copy.copy(self)
        ret._z = self._z.copy()
        ret._x = self._x.copy()
        ret._phase = self._phase.copy()
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
        x1 = cls._stack(a._x, b._num_paulis, False)
        x2 = cls._stack(b._x, a._num_paulis)
        z1 = cls._stack(a._z, b._num_paulis, False)
        z2 = cls._stack(b._z, a._num_paulis)
        phase1 = (
            np.vstack(b._num_paulis * [a._phase])
            .transpose(1, 0)
            .reshape(a._num_paulis * b._num_paulis)
        )
        phase2 = cls._stack(b._phase, a._num_paulis)
        z = np.hstack([z2, z1])
        x = np.hstack([x2, x1])
        phase = np.mod(phase1 + phase2, 4)
        return BasePauli(z, x, phase)

    # pylint: disable=arguments-differ
    def compose(self, other, qargs=None, front=False, inplace=False):
        """Return the composition of Paulis.

        Args:
            a ({cls}): an operator object.
            b ({cls}): an operator object.
            qargs (list or None): Optional, qubits to apply dot product
                                  on (default: None).
            inplace (bool): If True update in-place (default: False).

        Returns:
            {cls}: The operator a.compose(b)

        Raises:
            QiskitError: if number of qubits of other does not match qargs.
        """.format(
            cls=type(self).__name__
        )
        # Validation
        if qargs is None and other.num_qubits != self.num_qubits:
            raise QiskitError(
                "other {} must be on the same number of qubits.".format(type(self).__name__)
            )

        if qargs and other.num_qubits != len(qargs):
            raise QiskitError(
                "Number of qubits of the other {} does not match qargs.".format(type(self).__name__)
            )

        if other._num_paulis not in [1, self._num_paulis]:
            raise QiskitError(
                "Incompatible BasePaulis. Second list must "
                "either have 1 or the same number of Paulis."
            )

        # Compute phase shift
        if qargs is not None:
            x1, z1 = self._x[:, qargs], self._z[:, qargs]
        else:
            x1, z1 = self._x, self._z
        x2, z2 = other._x, other._z

        # Get phase shift
        phase = self._phase + other._phase
        if front:
            phase += 2 * np.sum(np.logical_and(x1, z2), axis=1)
        else:
            phase += 2 * np.sum(np.logical_and(z1, x2), axis=1)

        # Update Pauli
        x = np.logical_xor(x1, x2)
        z = np.logical_xor(z1, z2)

        if qargs is None:
            if not inplace:
                return BasePauli(z, x, phase)
            # Inplace update
            self._x = x
            self._z = z
            self._phase = phase
            return self

        # Qargs update
        ret = self if inplace else self.copy()
        ret._x[:, qargs] = x
        ret._z[:, qargs] = z
        ret._phase = np.mod(phase, 4)
        return ret

    def _multiply(self, other, roundit=True):
        """Return the {cls} other * self.

        Args:
            other (complex): a complex number(s) in ``[1, -1j, -1, 1j]``.
            round (bool): Set True to round components of other. Default=True
        Returns:
            {cls}: the {cls} other * self.

        Raises:
            QiskitError: if the phase is not in the set ``[1, -1j, -1, 1j]``.
        """.format(cls=type(self).__name__)
        phase = self.coeff_to_exponent(
                other,
                output_phase_format=PauliRep.__INTERNAL_PHASE_REP_FORMAT__,
                roundit=roundit
            )
        return BasePauli(self._z, self._x, np.mod(self._phase + phase, 4))

    def conjugate(self):
        """Return the conjugate of each Pauli in the list."""
        complex_phase = np.mod(self._phase, 2)
        if np.all(complex_phase == 0):
            return self
        return BasePauli(self._z, self._x, np.mod(self._phase + 2 * complex_phase, 4))

    def transpose(self):
        """Return the transpose of each Pauli in the list."""
        # Transpose sets Y -> -Y. This has effect on changing the phase
        parity_y = self.count_y() % 2
        if np.all(parity_y == 0):
            return self
        return BasePauli(self._z, self._x, np.mod(self._phase + 2 * parity_y, 4))

    def commutes(self, other, qargs=None):
        """Return True if Pauli that commutes with other.

        Args:
            other (BasePauli): another BasePauli operator.
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
        if qargs is not None:
            inds = list(qargs)
            x1, z1 = self._x[:, inds], self._z[:, inds]
        else:
            x1, z1 = self._x, self._z
        a_dot_b = np.mod(np.sum(np.logical_and(x1, other._z), axis=1), 2)
        b_dot_a = np.mod(np.sum(np.logical_and(z1, other._x), axis=1), 2)
        return a_dot_b == b_dot_a

    def evolve(self, other, qargs=None):
        r"""Heisenberg picture evolution of a Pauli by a Clifford.

        This returns the Pauli :math:`P^\prime = C^\dagger.P.C`.

        Args:
            other (BasePauli or QuantumCircuit): The Clifford circuit to evolve by.
            qargs (list): a list of qubits to apply the Clifford to.

        Returns:
            BasePauli: the Pauli :math:`C^\dagger.P.C`.

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
        if isinstance(other, BasePauli):
            ret = self.compose(other.adjoint(), qargs=qargs)
            ret = ret.compose(other, front=True, qargs=qargs)
            return ret

        # Otherwise evolve by the inverse circuit to compute C^dg.P.C
        return self.copy()._append_circuit(other.inverse(), qargs=qargs)

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

    def __imul__(self, other):
        return self.compose(other, front=True, inplace=True)

    def __neg__(self):
        ret = copy.copy(self)
        ret._phase = np.mod(self._phase + 2, 4)
        return ret

    def count_y(self):
        """Count the number of I Pauli's"""
        return self._count_y(self._x, self._z)

    @staticmethod
    def _stack(array, size, vertical=True):
        """Stack array."""
        if size == 1:
            return array
        if vertical:
            return np.vstack(size * [array]).reshape((size * len(array),) + array.shape[1:])
        return np.hstack(size * [array]).reshape((size * len(array),) + array.shape[1:])

    def _append_circuit(self, circuit, qargs=None):
        """Update BasePauli inplace by applying a Clifford circuit.

        Args:
            circuit (QuantumCircuit or Instruction): the gate or composite gate to apply.
            qargs (list or None): The qubits to apply gate to.

        Returns:
            BasePauli: the updated Pauli.

        Raises:
            QiskitError: if input gate cannot be decomposed into Clifford gates.
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
                raise QiskitError("Invalid Clifford gate name string {}".format(gate))
            name = gate
        else:
            # Assume gate is an Instruction
            name = gate.name

        # Apply gate if it is a Clifford basis gate
        if name in non_clifford:
            raise QiskitError("Cannot update Pauli with non-Clifford gate {}".format(name))
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
            raise QiskitError("Cannot apply Instruction: {}".format(gate.name))
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
                    "Cannot apply Instruction with classical registers: {}".format(instr.name)
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


def _evolve_h(base_pauli, qubit):
    """Update P -> H.P.H"""
    x = base_pauli._x[:, qubit].copy()
    z = base_pauli._z[:, qubit].copy()
    base_pauli._x[:, qubit] = z
    base_pauli._z[:, qubit] = x
    base_pauli._phase += 2 * np.logical_and(x, z).T
    return base_pauli


def _evolve_s(base_pauli, qubit):
    """Update P -> S.P.Sdg"""
    x = base_pauli._x[:, qubit]
    base_pauli._z[:, qubit] ^= x
    base_pauli._phase += x.T
    return base_pauli


def _evolve_sdg(base_pauli, qubit):
    """Update P -> Sdg.P.S"""
    x = base_pauli._x[:, qubit]
    base_pauli._z[:, qubit] ^= x
    base_pauli._phase -= x.T
    return base_pauli


# pylint: disable=unused-argument
def _evolve_i(base_pauli, qubit):
    """Update P -> P"""
    return base_pauli


def _evolve_x(base_pauli, qubit):
    """Update P -> X.P.X"""
    base_pauli._phase += 2 * base_pauli._z[:, qubit].T
    return base_pauli


def _evolve_y(base_pauli, qubit):
    """Update P -> Y.P.Y"""
    base_pauli._phase += 2 * base_pauli._x[:, qubit].T + 2 * base_pauli._z[:, qubit].T
    return base_pauli


def _evolve_z(base_pauli, qubit):
    """Update P -> Z.P.Z"""
    base_pauli._phase += 2 * base_pauli._x[:, qubit].T
    return base_pauli


def _evolve_cx(base_pauli, qctrl, qtrgt):
    """Update P -> CX.P.CX"""
    base_pauli._x[:, qtrgt] ^= base_pauli._x[:, qctrl]
    base_pauli._z[:, qctrl] ^= base_pauli._z[:, qtrgt]
    return base_pauli


def _evolve_cz(base_pauli, q1, q2):
    """Update P -> CZ.P.CZ"""
    x1 = base_pauli._x[:, q1].copy()
    x2 = base_pauli._x[:, q2].copy()
    base_pauli._z[:, q1] ^= x2
    base_pauli._z[:, q2] ^= x1
    base_pauli._phase += 2 * np.logical_and(x1, x2).T
    return base_pauli


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


def _evolve_swap(base_pauli, q1, q2):
    """Update P -> SWAP.P.SWAP"""
    x1 = base_pauli._x[:, q1].copy()
    z1 = base_pauli._z[:, q1].copy()
    base_pauli._x[:, q1] = base_pauli._x[:, q2]
    base_pauli._z[:, q1] = base_pauli._z[:, q2]
    base_pauli._x[:, q2] = x1
    base_pauli._z[:, q2] = z1
    return base_pauli
