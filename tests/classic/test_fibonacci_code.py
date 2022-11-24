import unittest

from qiskit_qec.codes.classic.fibonacci_code import generate_fibonacci_code_word
from qiskit_qec.decoders.classic.fibonacci_decoding.fibonacci_decoder import ClassicFibDecoder
from qiskit_qec.noise.classic.swath_error_generator import generate_swath_error


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
                    codeword = generate_fibonacci_code_word(
                        L
                    )  # generate an initial codeword. The default one bit at bottom center and cellular automata rules upward
                    error_board, error_mask = generate_swath_error(
                        codeword, L, probability_of_error=p
                    )  # setting width to L and vertical=True makes iid noise
                    f = ClassicFibDecoder(error_board)  # give this class the errored codeword
                    f.decode_fib_code()
                    f.board.shape = (L // 2, L)
                    if (f.board == codeword).all():
                        success_no += 1
                p_success = success_no / num_shots
                self.assertAlmostEqual(
                    p_success,
                    accuracy_sol,
                    delta=0.1,  # small decoders are pretty vulnerable to variations in success rates
                ), f"L={L} with p={p} gave: {p_success} instead of expected: {accuracy_sol}"
