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
"""Spanning Error Generator for Fibonacci Code"""
import logging
from typing import Optional

import numpy as np

from qiskit_qec.codes.classic.fibonacci_code import ClassicFibonacciCode
from qiskit_qec.exceptions import QiskitQECError

logger = logging.getLogger(__name__)


def generate_spanning_error(
    code: ClassicFibonacciCode,
    offset: int = 0,
    probability_of_error: Optional[float] = 1e-7,
    width: Optional[int] = None,
    is_vertical: Optional[bool] = True,
) -> (np.array, np.array):
    """Generates spanning errors for fibonacci code (ClassicFibonacciCode)
    that can be fed to the fibonacci decoder (ClassicFibonacciSpanningErrorDecoder)
    See section III subsection D figure 7 on arXiv:2002.11738 for more details

    Args:
        code (ClassicFibonacciCode): Fibonacci Code, has information about the codeword that
        generate_spanning_error will use as a base for it's error
        offset (int, optional): Offset for where spanning error begins. Offset will be a
        horizontal offset from the left if is_vertical=False and, otherwise,
        a vertical offsetfrom the top if is_vertical=True). Defaults to 0.
        probability_of_error (float, optional): Error probability. Defaults to 1e-7.
        width (int, optional): How wide the error should be. If is_vertical=True, width
        willbe how wide the spanning error is across from left to right. Otherwrise, if
        is_vertical=False width will be the height of the spanning error from top to bottom.
        Defaults to None.
        is_vertical (bool, optional): If true, the error will span top to bottom. Otherwise,
        the error will span left to right. Defaults to True.

    Raises:
        QiskitQECError: probability_of_error must be <= 1 and >=0
    Returns:
        error_board (np.array): element-wise "bitwise XOR" of error_mask and codeword
        error_mask (np.array): an np.array indicating which bits are flipped by errors
    """

    codeword = code.code_word
    if not width:
        width = code.num_logical_bits
    if probability_of_error < 0 or probability_of_error > 1:
        raise QiskitQECError(
            f"probability_of_error must be > 0 and < 1 but is: {probability_of_error}"
        )
    error_mask = np.zeros(codeword.shape, dtype=int)

    num_rows = len(codeword)
    num_col = len(codeword[0])
    if is_vertical:
        if probability_of_error == 1:
            error_mask[:, offset : width + offset] = 1
        else:
            for i in range(num_rows):
                for j in range(width):
                    error_mask[i][(j + offset) % num_col] = np.random.choice(
                        [0, 1], p=[1 - probability_of_error, probability_of_error]
                    )

    else:
        height = width
        if probability_of_error == 1:
            error_mask[offset : offset + height, :] = 1
        else:
            for i in range(width):
                for j in range(num_col):
                    error_mask[(i + offset) % num_rows][j] = np.random.choice(
                        [0, 1], p=[1 - probability_of_error, probability_of_error]
                    )

    error_board = error_mask ^ codeword
    return error_board, error_mask
