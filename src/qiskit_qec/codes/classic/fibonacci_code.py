# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020
#
# This code is licensed under the Apache num_logical_bitsicense, Version 2.0. You may
# obtain a copy of this license in the num_logical_bitsICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/num_logical_bitsICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Defines the Fibonacci Code"""
from typing import Optional

import numpy as np

from qiskit_qec.exceptions import QiskitQECError


class ClassicFibonacciCode:
    """Classic Fib Code"""

    def __init__(self, num_logical_bits: int, codeword_init_arr: Optional[np.array] = None):
        """Initialize data for Fibonacci Code.
        See arXiv:2002.11738 for more details
        We define the Fibonacci code on a two-dimensional lattice
            such that it encodes k = L logical bits
            using n = (L^2)/2 physical bits.
            The exact distance of the code,
            d > L, is unknown

        Args:
            num_logical_bits (int): Size of codeword, must be a power of 2
            codeword_init_arr (Optional[np.array]): Optional array that
            serves as the foundation for building up the rest of the codeword
            using cellular automaton update rules, as described in
            section II, formula 2 of arXiv:2002.11738. Defaults to None.

        Raises:
            QiskitQECError: Invalid input value for num_logical_bits
        """
        if num_logical_bits < 4 or not self._is_power_of_2(num_logical_bits):
            raise QiskitQECError("num_logical_bits must be >= 4 and a power of 2")
        self.num_logical_bits = num_logical_bits
        self.codeword_init_arr = self._gen_codeword_init_arr(codeword_init_arr)
        self.code_word = self._generate_fibonacci_code_word()

    def _is_power_of_2(self, n):
        return n > 0 and (n & (n - 1)) == 0

    def _gen_codeword_init_arr(self, codeword_init_arr):
        if codeword_init_arr is None:
            codeword_init_arr = np.zeros(self.num_logical_bits, dtype=int)
            codeword_init_arr[((self.num_logical_bits - 1) // 2)] = 1
        return codeword_init_arr

    def _generate_fibonacci_code_word(self):
        """Fibonacci Code Codeword Generator.
        Generates a valid configuration of classical data bits that can be decoded
        by ClassicFibonacciSpanningErrorDecoder and can undergo spanning errors (generate_spanning_error)

        See arXiv:2002.11738 section II figure 3 for more details
        """
        # generates from bottom row up
        rect_board = np.zeros((self.num_logical_bits // 2, self.num_logical_bits), dtype=np.int8)
        rect_board[(self.num_logical_bits // 2) - 1] = self.codeword_init_arr
        for row in range((self.num_logical_bits // 2) - 2, -1, -1):
            for bit in range(self.num_logical_bits):
                new_val = (
                    rect_board[row + 1][(bit - 1) % self.num_logical_bits]
                    ^ rect_board[row + 1][(bit) % self.num_logical_bits]
                    ^ rect_board[row + 1][(bit + 1) % self.num_logical_bits]
                )
                rect_board[row][bit] = new_val
        return rect_board

    def __str__(self) -> str:
        """Formatted string."""
        return (
            f"Fibonacci Code:\nnum_logical_bits:{self.num_logical_bits}"
            f"\ninit_code_word_array: {self.codeword_init_arr}"
        )

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"Fibonacci Code:\nnum_logical_bits:{self.num_logical_bits}\ninit_code_word_array:"
            f"\n{self.codeword_init_arr}\ncodeword:\n{self.code_word}"
        )
