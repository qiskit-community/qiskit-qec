"""Tests for the subsystem CSS circuit-level matching decoder."""
import unittest

from qiskit_aer import Aer

from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.decoders.circuit_matching_decoder import temp_syndrome
from qiskit_qec.decoders.repetition_decoder import RepetitionDecoder
from qiskit_qec.circuits.repetition_code import RepetitionCodeCircuit
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class TestRepetitionCircuitMatcher(unittest.TestCase):
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
        self.pnm = pnm

        # 3-bit, 2 round repetition code
        self.code_circuit = RepetitionCodeCircuit(3, 2)
        self.z_logical = self.code_circuit.css_z_logical

        # 5-bit, 2 round repetition code
        self.code_circuit_5 = RepetitionCodeCircuit(5, 2)
        self.z_logical_5 = self.code_circuit.css_z_logical

    def test_no_errors(self, method="rustworkx"):
        """Test the case with no errors using rustworkx."""

        def gint(c):
            """Casts to int if possible"""
            if c.isnumeric():
                return int(c)
            else:
                return c

        shots = 100
        seed = 100
        for logical in ["0", "1"]:
            dec = RepetitionDecoder(self.code_circuit, self.pnm, method, False, logical)
            qc = self.code_circuit.circuit[logical]
            backend = Aer.get_backend("aer_simulator")
            options = {"method": "stabilizer", "shots": shots, "seed_simulator": seed}
            result = backend.run(qc, **options).result()
            counts = result.get_counts(qc)
            dec.update_edge_weights(self.pnm)
            failures = 0
            for outcome, _ in counts.items():
                reversed_outcome = list(map(gint, outcome[::-1]))
                corrected_outcomes = dec.process(reversed_outcome)
                fail = temp_syndrome(corrected_outcomes, self.z_logical)
                failures += str(fail[0]) != logical
            self.assertEqual(failures, 0)

    def test_no_errors_pymatching(self):
        """Test the case with no errors using pymatching."""
        self.test_no_errors(method="pymatching")

    def test_correct_single_errors(self, method="rustworkx"):
        """Test the case with single faults using rustworkx."""
        for logical in ["0", "1"]:
            dec = RepetitionDecoder(self.code_circuit, self.pnm, method, False, logical)
            qc = self.code_circuit.circuit[logical]
            dec.update_edge_weights(self.pnm)
            fe = FaultEnumerator(qc, method="stabilizer", model=self.pnm)
            for faultpath in fe.generate():
                outcome = faultpath[3]
                corrected_outcomes = dec.process(outcome)
                fail = temp_syndrome(corrected_outcomes, self.z_logical)
                self.assertEqual(str(fail[0]), logical)

    def test_correct_single_errors_pymatching(self):
        """Test the case with two faults using pymatching."""
        self.test_correct_single_errors(method="pymatching")

    def test_error_pairs(self, dec_method="rustworkx", fe_method="stabilizer"):
        """Test the case with two faults on a d=5 code using rustworkx."""
        expected_failures = {"0": 0, "1": 0}
        for logical in ["0", "1"]:
            dec = RepetitionDecoder(self.code_circuit_5, self.pnm, dec_method, False, logical)
            dec.update_edge_weights(self.pnm)
            qc = self.code_circuit_5.circuit[logical]
            fe = FaultEnumerator(qc, order=2, method=fe_method, model=self.pnm)
            failures = 0
            for faultpath in fe.generate():
                outcome = faultpath[3]
                corrected_outcomes = dec.process(outcome)
                fail = temp_syndrome(corrected_outcomes, self.z_logical_5)
                failures += str(fail[0]) != logical
            self.assertEqual(failures, expected_failures[logical])

    def test_error_pairs_pymatching(self):
        """Test the case with two faults on a d=5 code using pymatching."""
        self.test_error_pairs(dec_method="pymatching", fe_method="stabilizer")


if __name__ == "__main__":
    unittest.main()
