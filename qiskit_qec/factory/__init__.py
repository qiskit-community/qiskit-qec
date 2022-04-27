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
===================================
Factory (:mod:`qiskit_qec.factory`)
===================================

.. currentmodule:: qiskit_qec.factory

Factory module classes and functions
====================================

.. autosummary::
    :toctree: ../stubs/

    Factory
    SubSystemFactory
    CSSSubSystemFactory
    HeavyHexFactory
"""

from .subsystem.subsystem_factory import SubSystemFactory
from .subsystem.css.css_subsystem_factory import CSSSubSystemFactory
from .subsystem.css.heavy_hex_factory import HeavyHexFactory
from .factory import Factory
