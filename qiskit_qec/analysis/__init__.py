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
=====================================
Analysis (:mod:`qiskit_qec.analysis`)
=====================================

.. currentmodule:: qiskit_qec.analysis

This module contains an :class:`ErrorPropagator` a circuit error propagator interface.

Analysis module classes and functions
=====================================

.. autosummary::
    :toctree: ../stubs/

    ErrorPropagator
    CErrorPropagator
    PyErrorPropagator
    EPSelector
    FaultEnumerator
    minimum_distance
"""

from .pyerrorpropagator import PyErrorPropagator
from .cerrorpropagator import CErrorPropagator
from .epselector import EPSelector
from .errorpropagator import ErrorPropagator
from .properties import minimum_distance
from .faultenumerator import FaultEnumerator
