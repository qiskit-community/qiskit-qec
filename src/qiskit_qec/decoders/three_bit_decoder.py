"""Three bit decoder."""
from typing import Tuple, List

from qiskit import QuantumCircuit

from qiskit_qec.decoders.circuit_matching_decoder import CircuitModelMatchingDecoder
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class ThreeBitDecoder(CircuitModelMatchingDecoder):
    """Simple 3-bit code matching decoder."""

    def __init__(
        self,
        n: int,
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
            n,
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
