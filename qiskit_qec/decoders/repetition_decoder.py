"""Matching Decoder for Repetition Codes."""
from typing import Tuple, List

from qiskit_qec.decoders.circuit_matching_decoder import CircuitModelMatchingDecoder
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class RepetitionDecoder(CircuitModelMatchingDecoder):
    """Instance of CircuitModelMatchingDecoder for use with
    circuits from RepetitionCodeCircuit."""

    def __init__(
        self,
        code_circuit,
        model: PauliNoiseModel,
        method: str,
        uniform: bool,
        logical: str,
    ):
        """Constructor."""
        self.resets = code_circuit.resets
        super().__init__(
            code_circuit.css_x_gauge_ops,
            code_circuit.css_x_stabilizer_ops,
            code_circuit.css_x_boundary,
            code_circuit.css_z_gauge_ops,
            code_circuit.css_z_stabilizer_ops,
            code_circuit.css_z_boundary,
            code_circuit.circuit[logical],
            model,
            code_circuit.basis,
            code_circuit.round_schedule,
            code_circuit.blocks,
            method,
            uniform,
        )

    def _partition_outcomes(
        self, blocks: int, round_schedule: str, outcome: List[int]
    ) -> Tuple[List[List[int]], List[List[int]], List[int]]:
        """Extract measurement outcomes."""
        # split into gauge and final outcomes
        outcome = "".join([str(c) for c in outcome])
        outcome = outcome.split(" ")
        gs = outcome[0:-1]
        gauge_outcomes = [[int(c) for c in r] for r in gs]
        finals = outcome[-1]
        # if circuit did not use resets, construct standard output
        if not self.resets:
            for i, layer in enumerate(gauge_outcomes):
                for j, gauge_op in enumerate(layer):
                    if i > 0:
                        gauge_outcomes[i][j] = (gauge_op + gauge_outcomes[i - 1][j]) % 2
        # assign outcomes to the correct gauge ops
        if round_schedule == "z":
            x_gauge_outcomes = []
            z_gauge_outcomes = gauge_outcomes
        else:
            x_gauge_outcomes = gauge_outcomes
            z_gauge_outcomes = []
        final_outcomes = [int(c) for c in finals]

        return x_gauge_outcomes, z_gauge_outcomes, final_outcomes
