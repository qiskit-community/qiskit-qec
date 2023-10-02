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
import unittest

from qiskit_qec.codes.classic.fibonacci_code import ClassicFibonacciCode
from qiskit_qec.decoders.classic.fibonacci_decoder import (
    ClassicFibonacciSpanningErrorDecoder,
)
from qiskit_qec.noise.classic.fibonacci_spanning_noise_generator import (
    generate_spanning_error,
)


class TestClassicFibCode(unittest.TestCase):
    """Test Classic Fib Code"""

    def test_full_fib_decoding(self):
        test_cases = [
            [4, [(0.01, 1.0), (0.105, 0.8), (0.2, 0.51)]],
            [8, [(0.01, 1.0), (0.105, 0.84), (0.2, 0.43)]],
        ]
        num_shots = 100

        for test in test_cases:
            L = test[0]
            for p, accuracy_sol in test[1]:
                success_no = 0
                for round in range(num_shots):
                    code = ClassicFibonacciCode(
                        L
                    )  # generate an initial codeword. The default one bit at bottom center and cellular automata rules upward

                    error_board, error_mask = generate_spanning_error(
                        code, probability_of_error=p
                    )  # setting width to L and vertical=True makes iid noise
                    f = ClassicFibonacciSpanningErrorDecoder(
                        error_board
                    )  # give this class the errored codeword
                    f.decode_fib_code()
                    f.board.shape = (L // 2, L)
                    if (f.board == code.code_word).all():
                        success_no += 1
                p_success = success_no / num_shots
                self.assertAlmostEqual(
                    p_success,
                    accuracy_sol,
                    delta=0.1,  # small decoders are pretty vulnerable to variations in success rates
                ), f"L={L} with p={p} gave: {p_success} instead of expected: {accuracy_sol}"
