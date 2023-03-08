"""Matching Decoder for Repetition Codes."""
from typing import Tuple, List

from qiskit_qec.decoders.circuit_matching_decoder import CircuitModelMatchingDecoder
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel
from qiskit_qec.decoders.decoding_graph import DecodingGraph


class RepetitionDecoder(CircuitModelMatchingDecoder):
    """Instance of CircuitModelMatchingDecoder for use with
    circuits from RepetitionCodeCircuit.

    Args:
        code_circuit: The QEC code circuit object for which this decoder
                will be used.
        model: Noise model used to generate syndrome graph.
        uniform: Whether to use uniform weights for the syndrome graph.
        logical: Logical value for the circuit to be used.
    """

    def __init__(
        self,
        code_circuit,
        model: PauliNoiseModel,
        method: str,
        uniform: bool,
        logical: str,
    ):
        """Constructor."""
        self.code_circuit = code_circuit
        dg = DecodingGraph(code_circuit)
        super().__init__(
            code_circuit.n,
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
            dg,
        )

    def _partition_outcomes(
        self, blocks: int, round_schedule: str, outcome: List[int]
    ) -> Tuple[List[List[int]], List[List[int]], List[int]]:
        """Extract measurement outcomes."""
        return self.code_circuit.partition_outcomes(round_schedule, outcome)
