# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2019.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# pylint: disable=invalid-name

"""Run codes and decoders."""

import itertools
import unittest
from random import choices

from qiskit_qec.circuits import ArcCircuit, RepetitionCodeCircuit
from qiskit_qec.decoders import PyMatchingDecoder
from qiskit_qec.codes.hhc import HHC
from qiskit_qec.circuits.css_code import CSSCodeCircuit
from qiskit_qec.utils.stim_tools import get_counts_via_stim


class TestMatching(unittest.TestCase):
    """Test the PyMatching decoder."""

    def test_repetition_codes(self, code):
        """
        Test on repetition codes.
        """

        d = 8
        p = 0.1
        N = 1000

        codes = []
        codes.append(RepetitionCodeCircuit(8, 1))
        codes.append(ArcCircuit([(2 * j, 2 * j + 1, 2 * (j + 1)) for j in range(d - 1)], 1))
        codes.append(
            ArcCircuit([(2 * j, 2 * j + 1, (2 * (j + 1)) % (2 * d - 2)) for j in range(d - 1)], 1)
        )
        for c, code in enumerate(codes):
            matcher = PyMatchingDecoder(code)
            min_error_num = code.d
            min_error_string = ""
            for _ in range(N):
                string = "".join([choices(["1", "0"], [1 - p, p])[0] for _ in range(d)])
                for _ in range(code.T):
                    string = string + " " + "0" * (d - 1)
                    # get and check corrected_z_logicals
                    corrected_z_logicals = matcher.process(string)
                    for node in matcher.decoding_graph.logical_nodes:
                        if node.index < len(corrected_z_logicals):
                            error = corrected_z_logicals[node.index] != 1
                            if error:
                                error_num = string.split(" ", maxsplit=1)[0].count("0")
                                if error_num < min_error_num:
                                    min_error_num = error_num
                                    min_error_string = string
            # check that min num errors to cause logical errors >d/2
            self.assertTrue(
                min_error_num >= d / 2,
                str(min_error_num)
                + " errors cause logical error despite d="
                + str(code.d)
                + " for code "
                + str(c)
                + " with "
                + min_error_string
                + ".",
            )
            print(c, min_error_num, min_error_string)

    def test_css_codes(self, code):
        """
        Test on CSS codes.
        """
        d = 5
        p = 0.01
        N = 1000

        codes = [CSSCodeCircuit(HHC(d), T=1, basis="x", noise_model=(p, p))]
        for c, code in enumerate(codes):
            matcher = PyMatchingDecoder(code)
            stim_counts = get_counts_via_stim(code.noisy_circuit["1"], N)
            errors = 0
            for string in stim_counts:
                corrected_z_logicals = matcher.process(string)
                if corrected_z_logicals:
                    if corrected_z_logicals[0] != 1:
                        errors += 1
        self.assertTrue(
            errors / N < p,
            "Logical error rate of"
            + str(errors / N)
            + " despite p="
            + str(p)
            + " for code "
            + str(c)
            + ".",
        )


if __name__ == "__main__":
    unittest.main()
