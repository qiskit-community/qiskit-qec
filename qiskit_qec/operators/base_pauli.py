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

from qiskit import QiskitError
from qiskit.quantum_info.operators.base_operator import BaseOperator
from qiskit.quantum_info.operators.mixins import AdjointMixin, MultiplyMixin
from qiskit_qec.linear.smatrix_api.smatrix import SMatrix
from qiskit_qec.utils import pauli_rep


# pylint: disable=no-member
class BasePauli(BaseOperator, AdjointMixin, MultiplyMixin):
    """Base Pauli"""

    # External string formats used when displaying Pauli's as strings
    EXTERNAL_SYMP_FORMAT = pauli_rep.DEFAULT_EXTERNAL_SYMP_FORMAT
    EXTERNAL_PHASE_REP_FORMAT = pauli_rep.DEFAULT_EXTERNAL_PHASE_REP_FORMAT
    EXTERNAL_PAULI_REP_FORMAT = EXTERNAL_PHASE_REP_FORMAT + EXTERNAL_SYMP_FORMAT

    def __init__(self, matrix, phase_exponent=0, stype="numpy"):
        """Generate Base Pauli

        Args:
            matrix ([type]): Corresponding pauli matrix
            phase_exponent (int, optional): i**phase_exponent Defaults to 0.
            stype (str, optional): Type of matrix being used Defaults to "numpy".
        """
        self.smatrix = SMatrix.get_methods(stype=stype)
        self.matrix = SMatrix(self.smatrix.atleast_2d(matrix), stype=stype)
        self.stype = stype
        self.phase_exponent = phase_exponent
        super().__init__(num_qubits=self.matrix.shape[1] >> 1)

    @property
    def num_paulis(self):
        """Return pauli shape"""
        return self.matrix.shape[0]

    @property
    def num_y(self):
        """Return the number of Y for each operator"""
        return self.smatrix.sum(self.smatrix.logical_and(self.x, self.z), axis=1)

    @property
    def symp_format(self):
        """Return the external symplectic matrix frmt"""
        return BasePauli.EXTERNAL_SYMP_FORMAT

    @symp_format.setter
    def symp_format(self, frmt=pauli_rep.DEFAULT_EXTERNAL_SYMP_FORMAT):
        """Set the external symplectic matrix format

        Args:
            frmt (str, optional): Symplectic matrix format.
            Defaults to pauli_rep.DEFAULT_EXTERNAL_SYMP_FORMAT.
        """
        assert frmt in pauli_rep.symp_formats(), QiskitError(
            f"Invalid symplectic matrix frmt: {frmt}. Must be one of {pauli_rep.symp_formats()}"
        )
        BasePauli.EXTERNAL_SYMP_FORMAT = frmt

    @property
    def phase_format(self):
        """Return the phase frmt"""
        return BasePauli.EXTERNAL_PHASE_REP_FORMAT

    @phase_format.setter
    def phase_format(self, frmt=pauli_rep.DEFAULT_EXTERNAL_PHASE_REP_FORMAT):
        """Set the phase frmt

        Args:
            frmt (str, optional): phase frmt.
            Defaults to pauli_rep.DEFAULT_EXTERNAL_PHASE_REP_FORMAT.
        """
        assert frmt in pauli_rep.phase_formats(), QiskitError(
            f"Invalid phase frmt: {frmt}. Must be one of {pauli_rep.phase_formats()}"
        )
        BasePauli.EXTERNAL_PHASE_REP_FORMAT = frmt

    @property
    def pauli_format(self):
        """Pauli format."""
        return BasePauli.EXTERNAL_PAULI_REP_FORMAT

    @pauli_format.setter
    def pauli_format(self, frmt=pauli_rep.DEFAULT_EXTERNAL_PAULI_REP_FORMAT):
        """Set the Pauli frmt

        Args:
            frmt (str, optional): Pauli frmt.
            Defaults to pauli_rep.DEFAULT_EXTERNAL_PAULI_REP_FORMAT.
        """
        assert frmt in pauli_rep.pauli_formats(), QiskitError(
            f"Invalid pauli frmt: {frmt}. Must be one of {pauli_rep.pauli_formats()}"
        )
        phase_format, symp_format = pauli_rep._split_rep(frmt)
        BasePauli.EXTERNAL_PHASE_REP_FORMAT = phase_format
        BasePauli.EXTERNAL_SYMP_FORMAT = symp_format
        BasePauli.EXTERNAL_PAULI_REP_FORMAT = phase_format + symp_format

    def compose(self, other, qargs=None, front=False):
        pass

    def tensor(self, other):
        pass

    def expand(self, other):
        pass

    def _multiply(self, other):
        pass

    def conjugate(self):
        pass

    def transpose(self):
        pass
