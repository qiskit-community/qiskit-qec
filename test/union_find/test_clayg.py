# This code is part of Qiskit.
#
# (C) Copyright IBM 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""Tests for template."""

import math
import random
import unittest

from qiskit_qec.circuits import RepetitionCodeCircuit
from qiskit_qec.decoders import ClAYGDecoder
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


def flip_with_probability(p, val):
    """
    Flips parity of val with probability p.
    """
    if random.random() <= p:
        val = (val + 1) % 2
    return val


def noisy_surface_code_outcome(d, p):
    """
    Generates outcome for surface code with phenomenological noise built in.
    """
    string = ""
    qubits = [0 for _ in range(d**2)]
    for _ in range(d):
        for qubit in qubits:
            qubit = flip_with_probability(p, qubit)
        # Top ancillas
        for i in [2 * i for i in range((d - 1) // 2)]:
            ancilla = (qubits[i] + qubits[i + 1]) % 2
            ancilla = flip_with_probability(p, ancilla)
            string += str(ancilla)
        for row in range(d - 1):
            offset = (row + 1) % 2
            for topleft in [offset + row * d + 2 * i for i in range((d - 1) // 2)]:
                ancilla = (
                    qubits[topleft]
                    + qubits[topleft + 1]
                    + qubits[topleft + d]
                    + qubits[topleft + d + 1]
                ) % 2
                ancilla = flip_with_probability(p, ancilla)
                string += str(ancilla)
        for i in [d * (d - 1) + 1 + 2 * i for i in range((d - 1) // 2)]:
            ancilla = (qubits[i] + qubits[i + 1]) % 2
            ancilla = flip_with_probability(p, ancilla)
            string += str(ancilla)
        string += " "
    for qubit in qubits:
        qubit = flip_with_probability(p, qubit)
        string += str(qubit)
    return string


class ClAYGDecoderTest(unittest.TestCase):
    """Tests will be here."""

    def setUp(self) -> None:
        # Bit-flip circuit noise model
        p = 0.05
        noise_model = PauliNoiseModel()
        noise_model.add_operation("cx", {"ix": 1, "xi": 1, "xx": 1})
        noise_model.add_operation("id", {"x": 1})
        noise_model.add_operation("reset", {"x": 1})
        noise_model.add_operation("measure", {"x": 1})
        noise_model.add_operation("x", {"x": 1, "y": 1, "z": 1})
        noise_model.set_error_probability("cx", p)
        noise_model.set_error_probability("x", p)
        noise_model.set_error_probability("id", p)
        noise_model.set_error_probability("reset", p)
        noise_model.set_error_probability("measure", p)
        self.noise_model = noise_model

        self.fault_enumeration_method = "stabilizer"

        return super().setUp()

    def test_surface_code_d3(self):
        """
        Test the ClAYG decoder on a surface code with d=3 and T=3
        with faults inserted by FaultEnumerator by checking if the syndromes
        have even parity (if it's a valid code state) and if the logical value measured
        is the one encoded by the circuit.

        This test won't complete atm, as the ClAYG decoder isn't able to decode some errors produced
        by this error model, because of the placement of the boundary nodes and clusters neutralizing
        when they shouldn't.
        """
        return

    def test_repetition_code_d5(self):
        """
        Test the ClAYG decoder on a repetition code with d=5 and T=5
        with faults inserted by FaultEnumerator by checking if the syndromes
        have even parity (if it's a valid code state) and if the logical value measured
        is the one encoded by the circuit.

        This test won't complete atm, as the ClAYG decoder isn't able to decode some errors produced
        by this error model, because of the placement of the boundary nodes and clusters neutralizing
        when they shouldn't.
        """
        return

    def test_error_rates(self):
        """
        Test the error rates using some repetition codes.
        """
        d = 8
        p = 0.01
        samples = 1000

        testcases = []
        testcases = [
            "".join([random.choices(["0", "1"], [1 - p, p])[0] for _ in range(d)])
            for _ in range(samples)
        ]
        codes = self.construct_codes(d)

        # now run them all and check it works
        for code in codes:
            decoder = ClAYGDecoder(code)
            z_logicals = code.css_z_logical[0]

            logical_errors = 0
            min_flips_for_logical = code.d
            for sample in range(samples):
                # generate random string
                string = ""
                for _ in range(code.T):
                    string += "0" * (d - 1) + " "
                string += testcases[sample]
                # get and check corrected_z_logicals
                outcome = decoder.process(string)
                # pylint: disable=consider-using-generator
                logical_outcome = sum([outcome[int(z_logical / 2)] for z_logical in z_logicals]) % 2
                if not logical_outcome == 0:
                    logical_errors += 1
                    min_flips_for_logical = min(min_flips_for_logical, string.count("1"))

            # check that error rates are at least <p^/2
            # and that min num errors to cause logical errors >d/3
            self.assertTrue(
                logical_errors / samples
                < (math.factorial(d)) / (math.factorial(int(d / 2)) ** 2) * p**4,
                "Logical error rate shouldn't exceed d!/((d/2)!^2)*p^(d/2).",
            )
            self.assertTrue(
                min_flips_for_logical >= d / 2,
                "Minimum amount of errors that also causes logical errors shouldn't be lower than d/2.",
            )

    def construct_codes(self, d):
        """
        Construct codes for the logical error rate test.
        """
        codes = []
        # TODO: Add more codes
        codes.append(RepetitionCodeCircuit(d=d, T=1))
        return codes


if __name__ == "__main__":
    unittest.main()
