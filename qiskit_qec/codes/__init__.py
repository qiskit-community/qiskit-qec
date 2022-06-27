# This code is part of Qiskit.
#
# (C) Copyright IBM 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
===============================
Codes (:mod:`qiskit_qec.codes`)
===============================

.. currentmodule:: qiskit_qec.codes

Codes module classes and functions
==================================

.. autosummary::
    :toctree: ../stubs/

    QECCode
    StabSubSystemCode
    CodeLibrary
    CodeLibrarian
    small_code
    all_small_codes
"""

from .qec_code import QECCode
from .stabsubsystemcodes import StabSubSystemCode
from .codebase import CodeLibrary, CodeLibrarian, small_code, all_small_codes
