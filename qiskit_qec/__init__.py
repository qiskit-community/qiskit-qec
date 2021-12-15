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

r"""
=============================================
Qiskit QEC module (:mod:`qiskit_qec`)
=============================================
.. currentmodule:: qiskit_qec
Qiskit Framework for Quantum Error Correction is an open-source framework for developers,
 experimentalist and theorists of Quantum Error Correction (QEC).
The top-level classes and submodules of qiskit_nature are:
.. autosummary::
   :toctree: ../stubs/
   QiskitQECError
Submodules
==========
.. autosummary::
   :toctree:
   codes
   info
   linear
   models
   modelers
   operators
   structures
   utils
"""

from .exceptions import QiskitQECError

__version__ = "0.0.1"
__all__ = ["__version__", "QiskitQECError"]
