# pylint: disable=missing-module-docstring
# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020, 2022
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# Part of the QEC framework


from qiskit_qec.analysis.pyerrorpropagator import PyErrorPropagator
from qiskit_qec.analysis.baseerrorpropagator import BaseErrorPropagator

from qiskit_qec.analysis.extensions import C_ERROR_PROPAGATOR

if C_ERROR_PROPAGATOR:
    from qiskit_qec.analysis.cerrorpropagator import CErrorPropagator


class ErrorPropagator(BaseErrorPropagator):
    """General ErrorPropagator class

    Return an error propagator.

    auto = try to automatically detect extensions
    c = force the compiled version
    py = force the python version

    """

    def __new__(cls, eptype: str = "auto", qreg_size: int = 1, creg_size: int = 1):
        """_summary_

        Args:
            qreg_size (int, optional): _description_. Defaults to 1.
            creg_size (int, optional): _description_. Defaults to 1.
            eptype (str, optional): _description_. Defaults to "auto".

        Raises:
            Exception: _description_
            Exception: _description_

        Returns:
            BaseErrorPropagator: _description_
        """
        if eptype == "auto":
            if C_ERROR_PROPAGATOR is True:
                return CErrorPropagator(qreg_size, creg_size)
            return PyErrorPropagator(qreg_size, creg_size)
        if eptype == "c":
            if C_ERROR_PROPAGATOR:
                return CErrorPropagator(qreg_size, creg_size)
            raise ImportError("C/C++ ErrorPropagator extension not loaded.")
        if eptype == "py":
            return PyErrorPropagator(qreg_size, creg_size)
        raise NotImplementedError("Unsupported error propagator type")
