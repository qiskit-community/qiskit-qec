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

from random import choices
import unittest
import math
from unittest import TestCase
from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.decoders import UnionFindDecoder
from qiskit_qec.circuits import SurfaceCodeCircuit, RepetitionCodeCircuit, ArcCircuit
from qiskit_qec.decoders.temp_code_util import temp_syndrome
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class UnionFindDecoderTest(TestCase):
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
        Test the union find decoder on a surface code with d=3 and T=3
        with faults inserted by FaultEnumerator by checking if the syndromes
        have even parity (if it's a valid code state) and if the logical value measured
        is the one encoded by the circuit.
        """
        for logical in ["0", "1"]:
            code = SurfaceCodeCircuit(d=3, T=3)
            decoder = UnionFindDecoder(code)
            fault_enumerator = FaultEnumerator(
                code.circuit[logical], method=self.fault_enumeration_method, model=self.noise_model
            )
            for fault in fault_enumerator.generate():
                outcome = "".join([str(x) for x in fault[3]])
                corrected_outcome = decoder.process(outcome[::-1])
                stabilizers = temp_syndrome(corrected_outcome, code.css_z_stabilizer_ops)
                for syndrome in stabilizers:
                    self.assertEqual(syndrome, 0)
                logical_measurement = temp_syndrome(corrected_outcome, [code.css_z_logical])[0]
                self.assertEqual(str(logical_measurement), logical)

    def test_repetition_code_d5(self):
        """
        Test the union find decoder on a repetition code with d=3 and T=3
        with faults inserted by FaultEnumerator by checking if the syndromes
        have even parity (if it's a valid code state) and if the logical value measured
        is the one encoded by the circuit.
        """
        for logical in ["0", "1"]:
            code = RepetitionCodeCircuit(d=5, T=5)
            decoder = UnionFindDecoder(code)
            fault_enumerator = FaultEnumerator(
                code.circuit[logical], method=self.fault_enumeration_method, model=self.noise_model
            )
            for fault in fault_enumerator.generate():
                outcome = "".join([str(x) for x in fault[3]])
                corrected_outcome = decoder.process(outcome[::-1])
                stabilizers = temp_syndrome(corrected_outcome, code.css_z_stabilizer_ops)
                for syndrome in stabilizers:
                    self.assertEqual(syndrome, 0)
                logical_measurement = temp_syndrome(corrected_outcome, code.css_z_logical)[0]
                self.assertEqual(str(logical_measurement), logical)

    def test_circular_arc_code(self):
        """
        Test the union find decoder on a circular ARC code with faults inserted
        by FaultEnumerator by checking if the syndromes have even parity
        (if it's a valid code state) and if the logical value measured
        is the one encoded by the circuit (only logical 0 for the ARC circuit,
        see issue #309).
        """
        links = [(0, 1, 2), (2, 3, 4), (4, 5, 6), (6, 7, 0)]
        code = ArcCircuit(links=links, T=len(links), resets=False)
        decoder = UnionFindDecoder(code)
        fault_enumerator = FaultEnumerator(
            code.circuit[code.base], method=self.fault_enumeration_method, model=self.noise_model
        )
        for fault in fault_enumerator.generate():
            outcome = "".join([str(x) for x in fault[3]])
            corrected_outcome = decoder.process(outcome[::-1])
            logical_measurement = temp_syndrome(
                corrected_outcome, [[int(q / 2) for q in code.z_logicals]]
            )[0]
            self.assertEqual(str(logical_measurement), "0")

    def test_error_rates(self):
        """
        Test the error rates using some ARCs.
        """
        d = 8
        p = 0.01
        samples = 1000

        testcases = []
        testcases = [
            "".join([choices(["0", "1"], [1 - p, p])[0] for _ in range(d)]) for _ in range(samples)
        ]
        codes = self.construct_codes(d)

        # now run them all and check it works
        for code in codes:
            decoder = UnionFindDecoder(code)
            if isinstance(code, ArcCircuit):
                z_logicals = code.z_logicals
            else:
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
                outcome = decoder.process(string[::-1])
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
        # parameters for test
        codes = []
        # add them to the code list
        # TODO: Add ARCs to tests as soon as a better general alternative to the peeling is found,
        # instead of just looking at what logicals are affected by checking for boundary nodes.
        # for links in [links_ladder, links_line, links_cross]:
        #      codes.append(ArcCircuit(links, 0))
        codes.append(RepetitionCodeCircuit(d=d, T=1))
        return codes


if __name__ == "__main__":
    unittest.main()
