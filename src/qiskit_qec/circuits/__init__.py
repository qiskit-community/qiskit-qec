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
Circuits (:mod:`qiskit_qec.circuits`)
=====================================

.. currentmodule:: qiskit_qec.circuits

Circuits module classes and functions
=====================================

.. autosummary::
    :toctree: ../stubs/

    CodeCircuit
    RepetitionCodeCircuit
    ArcCircuit
    SurfaceCodeCircuit
    CSSCodeCircuit
"""

from .code_circuit import CodeCircuit
from .repetition_code import RepetitionCodeCircuit, ArcCircuit
from .surface_code import SurfaceCodeCircuit
from .css_code import CSSCodeCircuit
from .stim_code_circuit import StimCodeCircuit
