import unittest
from typing import List, Dict, Tuple

from qiskit import execute, QuantumCircuit, Aer

from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel
from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.decoders.circuit_matching_decoder import (
    CircuitModelMatchingDecoder,
)


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
X_boundary = []
Z_boundary = [[0], [2]]

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

    def __init__(
        self,
        css_x_gauge_ops: List[Tuple[int]],
        css_x_stabilizer_ops: List[Tuple[int]],
        css_x_boundary: List[Tuple[int]],
        css_z_gauge_ops: List[Tuple[int]],
        css_z_stabilizer_ops: List[Tuple[int]],
        css_z_boundary: List[Tuple[int]],
        circuit: QuantumCircuit,
        model: PauliNoiseModel,
        basis: str,
        round_schedule: str,
        blocks: int,
        method: str,
        uniform: bool,
    ):
        """Constructor."""
        self.bits_per_round = 2
        super().__init__(
            css_x_gauge_ops,
            css_x_stabilizer_ops,
            css_x_boundary,
            css_z_gauge_ops,
            css_z_stabilizer_ops,
            css_z_boundary,
            circuit,
            model,
            basis,
            round_schedule,
            blocks,
            method,
            uniform,
        )

    def _partition_outcomes(
        self, blocks: int, round_schedule: str, outcome: List[int]
    ) -> Tuple[List[List[int]], List[List[int]], List[int]]:
        """Extract measurement outcomes."""
        assert blocks == 2
        assert round_schedule == "z"
        x_gauge_outcomes = []
        z_gauge_outcomes = [outcome[0:2], outcome[2:4]]
        final_outcomes = outcome[4:7]
        return x_gauge_outcomes, z_gauge_outcomes, final_outcomes


class TestCircuitMatcher(unittest.TestCase):
    def test_no_errors(self):
        shots = 100
        seed = 100
        dec = ThreeBitDecoder(
            X_stabilizers,
            X_stabilizers,
            X_boundary,
            Z_stabilizers,
            Z_stabilizers,
            Z_boundary,
            qc,
            pnm,
            "z",
            "z",
            2,
            "networkx",
            False,
        )
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
        dec = ThreeBitDecoder(
            X_stabilizers,
            X_stabilizers,
            X_boundary,
            Z_stabilizers,
            Z_stabilizers,
            Z_boundary,
            qc,
            pnm,
            "z",
            "z",
            2,
            "pymatching",
            False,
        )
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
        dec = ThreeBitDecoder(
            X_stabilizers,
            X_stabilizers,
            X_boundary,
            Z_stabilizers,
            Z_stabilizers,
            Z_boundary,
            qc,
            pnm,
            "z",
            "z",
            2,
            "networkx",
            False,
        )
        dec.update_edge_weights(pnm)
        fe = FaultEnumerator(qc, method="stabilizer", model=pnm)
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = code.logical_x_error(corrected_outcomes)  # TODO fix
            self.assertEqual(fail[0], 0)

    def test_correct_single_errors_pymatching(self):
        dec = ThreeBitDecoder(
            X_stabilizers,
            X_stabilizers,
            X_boundary,
            Z_stabilizers,
            Z_stabilizers,
            Z_boundary,
            qc,
            pnm,
            "z",
            "z",
            2,
            "pymatching",
            False,
        )
        dec.update_edge_weights(pnm)
        fe = FaultEnumerator(qc, method="stabilizer", model=pnm)
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = code.logical_x_error(corrected_outcomes)  # TODO fix
            self.assertEqual(fail[0], 0)

    def test_error_pairs(self):
        dec = ThreeBitDecoder(
            X_stabilizers,
            X_stabilizers,
            X_boundary,
            Z_stabilizers,
            Z_stabilizers,
            Z_boundary,
            qc,
            pnm,
            "z",
            "z",
            2,
            "networkx",
            False,
        )
        dec.update_edge_weights(pnm)
        fe = FaultEnumerator(qc, order=2, method="stabilizer", model=pnm)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = code.logical_x_error(corrected_outcomes)  # TODO fix
            failures += fail[0]
        self.assertEqual(failures, 128)

    def test_error_pairs_propagator_pymatching(self):
        dec = ThreeBitDecoder(
            X_stabilizers,
            X_stabilizers,
            X_boundary,
            Z_stabilizers,
            Z_stabilizers,
            Z_boundary,
            qc,
            pnm,
            "z",
            "z",
            2,
            "pymatching",
            False,
        )
        dec.update_edge_weights(pnm)
        fe = FaultEnumerator(qc, order=2, method="propagator", model=pnm)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            fail = code.logical_x_error(corrected_outcomes)  # TODO fix
            failures += fail[0]
        self.assertEqual(failures, 128)


if __name__ == "__main__":
    unittest.main()
