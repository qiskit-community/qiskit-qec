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
""" This is the main module that defines what a general QEC Code. """


class QECCode:
    """`QECCode` is the core class for all QEC codes and is the
    central construct from which all QEC codes are derived.

    A QEC Code consists of a Code, an Error Model, and a Recovery (Decoder + stuff)
    """

    def __init__(self) -> None:
        """Init function for QECCode"""
        pass
