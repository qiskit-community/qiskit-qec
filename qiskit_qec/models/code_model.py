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

"""Code models."""


class CodeModel:
    """Code model class.
    Args:
        parameter (str): Parameter for code model.
    Attributes:
        parameter (str): Parameter for code model.
    """

    def __init__(self, parameter: str):
        self.parameter = parameter

    def convert(self) -> bool:
        """Converts code model."""
        return True
