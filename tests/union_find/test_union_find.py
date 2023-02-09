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

"""Tests for template."""

from unittest import TestCase
from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.decoders import UnionFindDecoder
from qiskit_qec.circuits import SurfaceCodeCircuit
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

    def test_d3(self):
        for logical in ["0", "1"]:
            code = SurfaceCodeCircuit(d=3, T=3)
            decoder = UnionFindDecoder(code)
            fault_enumerator = FaultEnumerator(code.circuit[logical], method=self.fault_enumeration_method, model=self.noise_model)
            for i, fault in enumerate(fault_enumerator.generate()):
                outcome = "".join([str(x) for x in fault[3]])
                corrected_outcome = decoder.process(outcome)
                syndromes = temp_syndrome(corrected_outcome, code.css_z_gauge_ops)
                for syndrome in syndromes:
                    self.assertEqual(syndrome, 0)
