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
"""Module for XPPauli"""
from typing import Any, List, Optional, Union

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit.quantum_info.operators.mixins import generate_apidocs
from qiskit_qec.operators.base_xp_pauli import BaseXPPauli

# from qiskit_qec.utils import xp_pauli_rep


class XPPauli(BaseXPPauli):
    """`XPPauli` inherits from `BaseXPPauli`"""

    # Set the max XPPauli string size before truncation
    _truncate__ = 50

    # pylint: disable=unused-argument
    def __init__(
        self,
        data: Any,
        *,
        x: Union[List, np.ndarray, None] = None,
        z: Union[List, np.ndarray, None] = None,
        phase_exp: Union[int, np.ndarray, None] = None,
        precision: int = None,
        input_xppauli_encoding: str = BaseXPPauli.EXTERNAL_XP_PAULI_ENCODING,
    ):
        """XPPauli Init

        Args:
            data (str): Still in progress
            x ([type], optional): [description]. Defaults to None.
            z ([type], optional): [description]. Defaults to None.
            phase_exponent ([type], optional): [description]. Defaults to None.
            precision: Precision of XP operators. Must be an integer greater than or equal to two.

        Raises:
            QiskitError: Something went wrong.
        """
        if isinstance(data, np.ndarray):
            matrix = np.atleast_2d(data)
            if phase_exp is None:
                phase_exp = 0
        elif isinstance(data, BaseXPPauli):
            matrix = data.matrix[:, :]
            phase_exp = data._phase_exp[:]
            precision = data.precision
        # TODO: elif isinstance(data, (tuple, list)), isinstance(data, str),
        # isinstance(data, ScalarOp) may be implemented later, like Pauli class
        else:
            raise QiskitError("Invalid input data for XPPauli.")

        # Initialize BaseXPPauli
        if matrix.shape[0] != 1:
            raise QiskitError("Input is not a single XPPauli")

        super().__init__(matrix, phase_exp, precision)
        # TODO check if this is needed
        # self.vlist = self.matrix[0].tolist()

    # ---------------------------------------------------------------------
    # Property Methods
    # ---------------------------------------------------------------------

    @property
    def name(self):
        """Unique string identifier for operation type."""
        return "xppauli"

    # TODO: several @property methods exist in Pauli class, analogous methods
    # may be added here later.

    @property
    def phase_exp(self):
        """Return the group phase exponent for the Pauli."""
        # Convert internal Pauli encoding to external Pauli encoding
        # return pauli_rep.change_pauli_encoding(
        #     self._phase_exp, self.num_y, output_pauli_encoding=BasePauli.EXTERNAL_PAULI_ENCODING
        # )
        pass

    @phase_exp.setter
    def phase_exp(self, value):
        # Convert external Pauli encoding to the internal Pauli Encoding
        # self._phase_exp[:] = pauli_rep.change_pauli_encoding(
        #     value,
        #     self.num_y,
        #     input_pauli_encoding=BasePauli.EXTERNAL_PAULI_ENCODING,
        #     output_pauli_encoding=pauli_rep.INTERNAL_PAULI_ENCODING,
        #     same_type=False,
        # )
        pass

    @property
    def phase(self):
        """Return the complex phase of the Pauli"""
        # return pauli_rep.exp2cpx(self.phase_exp, input_encoding=BasePauli.EXTERNAL_PHASE_ENCODING)
        pass

    @property
    def x(self):
        """The x vector for the XPPauli."""
        return self.matrix[:, : self.num_qubits][0]

    @x.setter
    def x(self, val):
        self.matrix[:, : self.num_qubits][0] = val

    @property
    def z(self):
        """The z vector for the XPPauli."""
        return self.matrix[:, self.num_qubits :][0]

    @z.setter
    def z(self, val):
        self.matrix[:, self.num_qubits :][0] = val

    # ---------------------------------------------------------------------
    # BaseOperator methods
    # ---------------------------------------------------------------------

    def compose(
        self,
        other: Union["XPPauli", BaseXPPauli],
        qargs: Optional[List[int]] = None,
        front: bool = False,
        inplace: bool = False,
    ) -> "XPPauli":
        """Return the operator composition with another XPPauli.

        Args:
            other (XPPauli): a XPPauli object.
            qargs (list or None): Optional, qubits to apply dot product
                                  on (default: None).
            front (bool): If True compose using right operator multiplication,
                          instead of left multiplication [default: False].
            inplace (bool): If True update in-place (default: False).

        Returns:
            XPPauli: The composed XPPauli.

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
        if not isinstance(other, XPPauli):
            other = XPPauli(other)
        return XPPauli(super().compose(other, qargs=qargs, front=front, inplace=inplace))

    # ---------------------------------------------------------------------

    def unique_vector_rep(self) -> "XPPauli":
        return XPPauli(super().unique_vector_rep())

    def rescale_precision(self, new_precision) -> "XPPauli":
        return XPPauli(super().rescale_precision(new_precision))

    def antisymmetric_op(self) -> "XPPauli":
        return XPPauli(super().antisymmetric_op())

    def power(self, n) -> "XPPauli":
        return XPPauli(super().power(n))


# Update docstrings for API docs
generate_apidocs(XPPauli)
