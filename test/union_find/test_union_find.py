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

import unittest
from unittest import TestCase
from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.decoders import UnionFindDecoder
from qiskit_qec.circuits import SurfaceCodeCircuit
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class UnionFindDecoderTest(TestCase):
    """Tests for the UnionFind decoder not covered elsewhere"""

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
                corrected_outcome = decoder.process(outcome)
                self.assertTrue(corrected_outcome[0] == 0, "Correction for surface code failed.")


if __name__ == "__main__":
    unittest.main()
