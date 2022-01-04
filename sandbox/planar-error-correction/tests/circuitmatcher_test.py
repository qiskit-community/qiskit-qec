import unittest

from qiskit import execute, QuantumCircuit, Aer

from plerco.pec_python.qec_code.css_code import CSSCode
from plerco.pec_python.qec_noise_model.paulinoisemodel import PauliNoiseModel
from plerco.pec_python.qec_decoder.decoder_utils.faultenumerator import FaultEnumerator
from plerco.pec_python.qec_decoder.circuit_matching_decoder import (
    CircuitModelMatchingDecoder,
)
from plerco.pec_python.configuration.config import Config


# Bit-flip circuit noise model
p = 0.05
pnm = PauliNoiseModel()
pnm.add_operation("cx", {"ix": 1, "xi": 1, "xx": 1})
pnm.add_operation("id", {"x": 1})
pnm.add_operation("reset", {"x": 1})
pnm.add_operation("measure", {"x": 1})
pnm.add_operation("x", {"x": 1, "y": 1, "z": 1})
pnm.set_error_probability("cx", p)
pnm.set_error_probability("x", p)
pnm.set_error_probability("id", p)
pnm.set_error_probability("reset", p)
pnm.set_error_probability("measure", p)

# 3-bit repetition code
X_stabilizers = []
Z_stabilizers = [[0, 1], [1, 2]]
X_logical = [[0, 1, 2]]
Z_logical = [[0]]
code = CSSCode(
    X_stabilizers,
    Z_stabilizers,
    3,
    1,
    1,
    X_stabilizers,
    Z_stabilizers,
    Z_logical,
    X_logical,
)
code.z_boundary = [[0], [2]]  # extra info

# Construct the 3-bit syndrome measurement circuit
qc = QuantumCircuit(4, 7)
qc.reset(0)
qc.reset(1)
qc.reset(2)
for i in range(2):
    qc.reset(3)
    qc.cx(0, 3)
    qc.cx(1, 3)
    qc.measure(3, 0 + 2 * i)
    qc.reset(3)
    qc.cx(1, 3)
    qc.cx(2, 3)
    qc.measure(3, 1 + 2 * i)
qc.measure(0, 4)
qc.measure(1, 5)
qc.measure(2, 6)


class ThreeBitDecoder(CircuitModelMatchingDecoder):
    """Simple 3-bit code matching decoder."""

    def __init__(self, code, circuit, model, config):
        """Constructor."""
        config["circuit"]["rounds"] = 2
        config["circuit"]["round_schedule"] = "z"
        config["circuit"]["basis"] = "z"
        self.bits_per_round = 2
        super().__init__(code, circuit, model, config)

    def _partition_outcomes(self, rounds, round_schedule, outcome):
        """Extract measurement outcomes."""
        assert rounds == 2
        assert round_schedule == "z"
        x_gauge_outcomes = []
        z_gauge_outcomes = [outcome[0:2], outcome[2:4]]
        final_outcomes = outcome[4:7]
        return x_gauge_outcomes, z_gauge_outcomes, final_outcomes


class TestCircuitMatcher(unittest.TestCase):
    def test_no_errors(self):
        shots = 100
        seed = 100
        config = Config()
        config["decoder"]["method"] = "matching_networkx"
        dec = ThreeBitDecoder(code, qc, pnm, config)
        result = execute(
            qc,
            Aer.get_backend("qasm_simulator"),
            method="stabilizer",
            shots=shots,
            optimization_level=0,
            seed_simulator=seed,
        ).result()
        counts = result.get_counts(qc)
        dec.update_edge_weights(pnm)
        failures = 0
        for (outcome, num) in counts.items():
            reversed_outcome = list(map(int, outcome[::-1]))
            corrected_outcomes = dec.process(reversed_outcome)
            fail = code.logical_x_error(corrected_outcomes)
            failures += fail[0]
        self.assertEqual(failures, 0)

    def test_no_errors_pymatching(self):
        shots = 100
        seed = 100
        config = Config()
        config["decoder"]["method"] = "matching_pymatching"
        dec = ThreeBitDecoder(code, qc, pnm, config)
        result = execute(
            qc,
            Aer.get_backend("qasm_simulator"),
            method="stabilizer",
            shots=shots,
            optimization_level=0,
            seed_simulator=seed,
        ).result()
        counts = result.get_counts(qc)
        dec.update_edge_weights(pnm)
        failures = 0
        for (outcome, num) in counts.items():
            reversed_outcome = list(map(int, outcome[::-1]))
            corrected_outcomes = dec.process(reversed_outcome)
            fail = code.logical_x_error(corrected_outcomes)
            failures += fail[0]
        self.assertEqual(failures, 0)

    def test_correct_single_errors(self):
        config = Config()
        config["decoder"]["method"] = "matching_networkx"
        dec = ThreeBitDecoder(code, qc, pnm, config)
        dec.update_edge_weights(pnm)
        fe = FaultEnumerator(qc, method="stabilizer", model=pnm)
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = code.logical_x_error(corrected_outcomes)
            self.assertEqual(fail[0], 0)

    def test_correct_single_errors_pymatching(self):
        config = Config()
        config["decoder"]["method"] = "matching_pymatching"
        dec = ThreeBitDecoder(code, qc, pnm, config)
        dec.update_edge_weights(pnm)
        fe = FaultEnumerator(qc, method="stabilizer", model=pnm)
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = code.logical_x_error(corrected_outcomes)
            self.assertEqual(fail[0], 0)

    def test_error_pairs(self):
        config = Config()
        config["decoder"]["method"] = "matching_networkx"
        dec = ThreeBitDecoder(code, qc, pnm, config)
        dec.update_edge_weights(pnm)
        fe = FaultEnumerator(qc, order=2, method="stabilizer", model=pnm)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = code.logical_x_error(corrected_outcomes)
            failures += fail[0]
        self.assertEqual(failures, 128)

    def test_error_pairs_propagator_pymatching(self):
        config = Config()
        config["decoder"]["method"] = "matching_pymatching"
        dec = ThreeBitDecoder(code, qc, pnm, config)
        dec.update_edge_weights(pnm)
        fe = FaultEnumerator(qc, order=2, method="propagator", model=pnm)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = code.logical_x_error(corrected_outcomes)
            failures += fail[0]
        self.assertEqual(failures, 128)


if __name__ == "__main__":
    unittest.main()
