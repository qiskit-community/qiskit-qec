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
============================================
Analysis module (:mod:`qiskit_qec.analysis`)
============================================

.. currentmodule:: qiskit_qec.analysis

This module contains an :class:`ErrorPropagator` a circuit error propagator interface.

Analysis Class
==============

.. autosummary::
    :toctree: ../stubs/

    PyErrorPropagator
"""

from . import cerrorpropagator
from . import epselector
from . import errorpropagator
from . import pyerrorpropagator

from .pyerrorpropagator import PyErrorPropagator
