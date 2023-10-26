# -*- coding: utf-8 -*-

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

"""Test for the CSSCodeCircuit class for Heavy-HEX code with pymatching"""
import unittest
import pymatching

from qiskit_qec.codes.hhc import HHC
from qiskit_qec.circuits.css_code import CSSCodeCircuit
from qiskit_qec.decoders.decoding_graph import DecodingGraph
from qiskit_qec.utils.stim_tools import get_stim_circuits


class TestCircuitMatcher(unittest.TestCase):
    """Test for the CSSCodeCircuit class for Heavy-HEX code with pymatching"""

    def log_failure_dists(self, error_rate: float):
        """Constructs the stim circuit and the decoding graph (via a stim DetectorErrorModel)
        for the heavy-hex code and tests it on 10_000 samples for distance 3 and 5.
        Returns the logical failure for distance 3 and 5 at the specified error rate.
        Below ~100_000 shotss, the runtime is limited by the decoding graph construction,
        not the sampling."""
        dist_list = [3, 5]
        num_shots = 10_000
        log_fail_d = []
        for d in dist_list:
            code = HHC(d)
            css_code = CSSCodeCircuit(code, T=d, basis="x", noise_model=(error_rate, error_rate))
            graph = DecodingGraph(css_code).graph
            m = pymatching.Matching(graph)
            detectors, logicals = css_code.stim_detectors()
            stim_circuit = get_stim_circuits(
                css_code.noisy_circuit["0"], detectors=detectors, logicals=logicals
            )[0][0]
            stim_sampler = stim_circuit.compile_detector_sampler()
            num_correct = 0
            stim_samples = stim_sampler.sample(num_shots, append_observables=True)
            for sample in stim_samples:
                actual_observable = sample[-1]
                detectors_only = sample.copy()
                detectors_only[-1] = 0
                predicted_observable = m.decode(detectors_only)[0]
                num_correct += actual_observable == predicted_observable
            log_fail_d.append((num_shots - num_correct) / num_shots)
        return log_fail_d

    def test_HHC(self):
        """Tests the order of logical failure rates for two distances.
        One test is below threshold, one is above. (Threshold ~3.5% for the test code)
        Runtime is approx 5 seconds."""
        error_rate = 0.01  # this should be below threshold
        log_fail_dist = self.log_failure_dists(error_rate)
        self.assertTrue(log_fail_dist[0] > 2 * log_fail_dist[1])

        error_rate = 0.1  # this should be above threshold
        log_fail_dist = self.log_failure_dists(error_rate)
        self.assertTrue(log_fail_dist[0] < log_fail_dist[1])


if __name__ == "__main__":
    unittest.main()
