"""Tests for the subsystem CSS circuit-level matching decoder."""
import unittest

from qiskit import execute, QuantumCircuit, Aer

from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.decoders.circuit_matching_decoder import temp_syndrome
from qiskit_qec.decoders.three_bit_decoder import ThreeBitDecoder
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class TestCircuitMatcher(unittest.TestCase):
    """Tests for the three bit decoder example."""

    def setUp(self) -> None:
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
        x_stabilizers = []
        z_stabilizers = [[0, 1], [1, 2]]
        x_logical = [[0, 1, 2]]
        z_logical = [[0]]
        x_boundary = []
        z_boundary = [[0], [2]]

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

        self.pnm = pnm
        self.qc = qc
        self.x_stabilizers = x_stabilizers
        self.x_boundary = x_boundary
        self.z_stabilizers = z_stabilizers
        self.z_boundary = z_boundary
        self.z_logical = z_logical
        self.x_logical = x_logical

    def test_no_errors(self):
        """Test the case with no errors using rustworkx."""
        shots = 100
        seed = 100
        dec = ThreeBitDecoder(
            3,
            self.x_stabilizers,
            self.x_stabilizers,
            self.x_boundary,
            self.z_stabilizers,
            self.z_stabilizers,
            self.z_boundary,
            self.qc,
            self.pnm,
            "z",
            "z",
            2,
            "rustworkx",
            False,
        )
        result = execute(
            self.qc,
            Aer.get_backend("aer_simulator"),
            method="stabilizer",
            shots=shots,
            optimization_level=0,
            seed_simulator=seed,
        ).result()
        counts = result.get_counts(self.qc)
        dec.update_edge_weights(self.pnm)
        failures = 0
        for outcome, _ in counts.items():
            reversed_outcome = list(map(int, outcome[::-1]))
            corrected_outcomes = dec.process(reversed_outcome)
            fail = temp_syndrome(corrected_outcomes, self.z_logical)
            failures += fail[0]
        self.assertEqual(failures, 0)

    def test_no_errors_pymatching(self):
        """Test the case with no errors using pymatching."""
        shots = 100
        seed = 100
        dec = ThreeBitDecoder(
            3,
            self.x_stabilizers,
            self.x_stabilizers,
            self.x_boundary,
            self.z_stabilizers,
            self.z_stabilizers,
            self.z_boundary,
            self.qc,
            self.pnm,
            "z",
            "z",
            2,
            "pymatching",
            False,
        )
        result = execute(
            self.qc,
            Aer.get_backend("aer_simulator"),
            method="stabilizer",
            shots=shots,
            optimization_level=0,
            seed_simulator=seed,
        ).result()
        counts = result.get_counts(self.qc)
        dec.update_edge_weights(self.pnm)
        failures = 0
        for outcome, _ in counts.items():
            reversed_outcome = list(map(int, outcome[::-1]))
            corrected_outcomes = dec.process(reversed_outcome)
            fail = temp_syndrome(corrected_outcomes, self.z_logical)
            failures += fail[0]
        self.assertEqual(failures, 0)

    def test_correct_single_errors(self):
        """Test the case with single faults using rustworkx."""
        dec = ThreeBitDecoder(
            3,
            self.x_stabilizers,
            self.x_stabilizers,
            self.x_boundary,
            self.z_stabilizers,
            self.z_stabilizers,
            self.z_boundary,
            self.qc,
            self.pnm,
            "z",
            "z",
            2,
            "rustworkx",
            False,
        )
        dec.update_edge_weights(self.pnm)
        fe = FaultEnumerator(self.qc, method="stabilizer", model=self.pnm)
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = temp_syndrome(corrected_outcomes, self.z_logical)
            self.assertEqual(fail[0], 0)

    def test_correct_single_errors_uniform(self):
        """Test the case with single faults using rustworkx."""
        dec = ThreeBitDecoder(
            3,
            self.x_stabilizers,
            self.x_stabilizers,
            self.x_boundary,
            self.z_stabilizers,
            self.z_stabilizers,
            self.z_boundary,
            self.qc,
            self.pnm,
            "z",
            "z",
            2,
            "rustworkx",
            True,
        )
        dec.update_edge_weights(self.pnm)
        fe = FaultEnumerator(self.qc, method="stabilizer", model=self.pnm)
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = temp_syndrome(corrected_outcomes, self.z_logical)
            self.assertEqual(fail[0], 0)

    def test_correct_single_errors_pymatching(self):
        """Test the case with single faults using pymatching."""
        dec = ThreeBitDecoder(
            3,
            self.x_stabilizers,
            self.x_stabilizers,
            self.x_boundary,
            self.z_stabilizers,
            self.z_stabilizers,
            self.z_boundary,
            self.qc,
            self.pnm,
            "z",
            "z",
            2,
            "pymatching",
            False,
        )
        dec.update_edge_weights(self.pnm)
        fe = FaultEnumerator(self.qc, method="stabilizer", model=self.pnm)
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = temp_syndrome(corrected_outcomes, self.z_logical)
            self.assertEqual(fail[0], 0)

    def test_correct_single_errors_pymatching_uniform(self):
        """Test the case with single faults using pymatching."""
        dec = ThreeBitDecoder(
            3,
            self.x_stabilizers,
            self.x_stabilizers,
            self.x_boundary,
            self.z_stabilizers,
            self.z_stabilizers,
            self.z_boundary,
            self.qc,
            self.pnm,
            "z",
            "z",
            2,
            "pymatching",
            True,
        )
        dec.update_edge_weights(self.pnm)
        fe = FaultEnumerator(self.qc, method="stabilizer", model=self.pnm)
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = temp_syndrome(corrected_outcomes, self.z_logical)
            self.assertEqual(fail[0], 0)

    def test_error_pairs(self):
        """Test the case with two faults using rustworkx."""
        dec = ThreeBitDecoder(
            3,
            self.x_stabilizers,
            self.x_stabilizers,
            self.x_boundary,
            self.z_stabilizers,
            self.z_stabilizers,
            self.z_boundary,
            self.qc,
            self.pnm,
            "z",
            "z",
            2,
            "rustworkx",
            False,
        )
        dec.update_edge_weights(self.pnm)
        fe = FaultEnumerator(self.qc, order=2, method="stabilizer", model=self.pnm)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = temp_syndrome(corrected_outcomes, self.z_logical)
            failures += fail[0]
        self.assertEqual(failures, 140)

    def test_error_pairs_uniform(self):
        """Test the case with two faults using rustworkx."""
        dec = ThreeBitDecoder(
            3,
            self.x_stabilizers,
            self.x_stabilizers,
            self.x_boundary,
            self.z_stabilizers,
            self.z_stabilizers,
            self.z_boundary,
            self.qc,
            self.pnm,
            "z",
            "z",
            2,
            "rustworkx",
            True,
        )
        dec.update_edge_weights(self.pnm)
        fe = FaultEnumerator(self.qc, order=2, method="stabilizer", model=self.pnm)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = temp_syndrome(corrected_outcomes, self.z_logical)
            failures += fail[0]
        self.assertEqual(failures, 152)

    def test_error_pairs_propagator_pymatching(self):
        """Test the case with two faults using error propagator and pymatching."""
        dec = ThreeBitDecoder(
            3,
            self.x_stabilizers,
            self.x_stabilizers,
            self.x_boundary,
            self.z_stabilizers,
            self.z_stabilizers,
            self.z_boundary,
            self.qc,
            self.pnm,
            "z",
            "z",
            2,
            "pymatching",
            False,
        )
        dec.update_edge_weights(self.pnm)
        fe = FaultEnumerator(self.qc, order=2, method="propagator", model=self.pnm)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = temp_syndrome(corrected_outcomes, self.z_logical)
            failures += fail[0]
        self.assertEqual(failures, 128)

    def test_error_pairs_propagator_pymatching_uniform(self):
        """Test the case with two faults using error propagator and pymatching."""
        dec = ThreeBitDecoder(
            3,
            self.x_stabilizers,
            self.x_stabilizers,
            self.x_boundary,
            self.z_stabilizers,
            self.z_stabilizers,
            self.z_boundary,
            self.qc,
            self.pnm,
            "z",
            "z",
            2,
            "pymatching",
            True,
        )
        dec.update_edge_weights(self.pnm)
        fe = FaultEnumerator(self.qc, order=2, method="propagator", model=self.pnm)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = temp_syndrome(corrected_outcomes, self.z_logical)
            failures += fail[0]
        # For pymatching v0.x, there are 168 failures, whereas for pymatching >v2.0.0 there are 156.
        # The reason for the difference is that many of these test cases have degenerate solutions
        # (Both versions of pymatching are giving valid minimum-weight perfect matching solutions, but
        # the predictions they make are not always the same when there is more than one valid
        # minimum-weight solution.)
        self.assertTrue(failures in {156, 168})


if __name__ == "__main__":
    unittest.main()
