"""Object to construct decoders and decode the HHC."""

from typing import List, Tuple
import logging

from qiskit import QuantumCircuit

from qiskit_qec.decoders.decoding_graph import DecodingGraph
from qiskit_qec.decoders.circuit_matching_decoder import CircuitModelMatchingDecoder
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel
from qiskit_qec.decoders.temp_code_util import temp_syndrome


class HHCDecoder(CircuitModelMatchingDecoder):
    """Decoder for the heavy-hexagon compass code."""

    def __init__(
        self,
        n: int,
        css_x_gauge_ops: List[Tuple[int]],
        css_x_stabilizer_ops: List[Tuple[int]],
        css_x_boundary: List[int],
        css_z_gauge_ops: List[Tuple[int]],
        css_z_stabilizer_ops: List[Tuple[int]],
        css_z_boundary: List[int],
        circuit: QuantumCircuit,
        model: PauliNoiseModel,
        basis: str,
        round_schedule: str,
        blocks: int,
        method: str,
        uniform: bool,
        decoding_graph: DecodingGraph = None,
        annotate: bool = False,
    ):
        """Create a decoder object."""
        # Sum the total number of bits per round
        self.bits_per_round = 0
        self.round_schedule = round_schedule
        if not isinstance(self.round_schedule, str):
            raise Exception("expected round_schedule to be a string")
        if set(self.round_schedule) > set("xz"):
            raise Exception("expected round schedule of 'x', 'z' chars")
        for step in self.round_schedule:
            if step == "z":
                # need to include 2 flags per z-gauge outcome
                self.bits_per_round += 3 * len(css_z_gauge_ops)
            elif step == "x":
                self.bits_per_round += len(css_x_gauge_ops)
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
            decoding_graph,
            annotate,
        )

    def _partition_outcomes(self, blocks: int, round_schedule: str, outcome: List[int]):
        """Process the raw outcome and return results.

        blocks = number of blocks
        round_schedule = string of z and x characters
        outcome = list of 0, 1 outcomes

        Return lists x_gauge_outcomes, z_gauge_outcomes, final_outcomes.
        """
        # partition the outcome list by outcome type
        x_gauge_outcomes = []
        z_gauge_outcomes = []
        left_flag_outcomes = []
        right_flag_outcomes = []
        final_outcomes = []
        for r in range(blocks):
            bits_into_round = 0
            for rs in round_schedule:
                if rs == "z":
                    z_gauge_outcomes.append(
                        outcome[
                            r * self.bits_per_round
                            + bits_into_round : r * self.bits_per_round
                            + bits_into_round
                            + len(self.css_z_gauge_ops)
                        ]
                    )
                    bits_into_round += len(self.css_z_gauge_ops)
                    left_flag_outcomes.append(
                        outcome[
                            r * self.bits_per_round
                            + bits_into_round : r * self.bits_per_round
                            + bits_into_round
                            + len(self.css_z_gauge_ops)
                        ]
                    )
                    bits_into_round += len(self.css_z_gauge_ops)
                    right_flag_outcomes.append(
                        outcome[
                            r * self.bits_per_round
                            + bits_into_round : r * self.bits_per_round
                            + bits_into_round
                            + len(self.css_z_gauge_ops)
                        ]
                    )
                    bits_into_round += len(self.css_z_gauge_ops)
                if rs == "x":
                    x_gauge_outcomes.append(
                        outcome[
                            r * self.bits_per_round
                            + bits_into_round : r * self.bits_per_round
                            + bits_into_round
                            + len(self.css_x_gauge_ops)
                        ]
                    )
                    bits_into_round += len(self.css_x_gauge_ops)
        final_outcomes = outcome[-self.n :]
        # Process the flags
        logging.debug("left_flag_outcomes = %s", left_flag_outcomes)
        logging.debug("right_flag_outcomes = %s", right_flag_outcomes)
        logging.debug("x_gauge_outcomes (before deflag) = %s", x_gauge_outcomes)
        x_gauge_outcomes, final_outcomes = self._process_flags(
            self.blocks,
            self.round_schedule,
            x_gauge_outcomes,
            left_flag_outcomes,
            right_flag_outcomes,
            final_outcomes,
        )
        logging.debug("x_gauge_outcomes (after deflag) = %s", x_gauge_outcomes)
        return x_gauge_outcomes, z_gauge_outcomes, final_outcomes

    def _process_flags(
        self,
        blocks: int,
        round_schedule: str,
        x_gauge_outcomes: List[Tuple[int]],
        left_flag_outcomes: List[Tuple[int]],
        right_flag_outcomes: List[Tuple[int]],
        final_outcomes: List[int],
    ):
        """Process the flag data for a set of outcomes.

        The outcomes are 0-1 lists.

        For each Z gauge measurement, we look at the left and right flag
        outcomes. If only the left flag is raised, we apply (in software)
        a Z to the qubit in the upper left corner of the plaquette, after
        that round of Z gauge measurements. If only the right flag is
        raised, we instead apply Z to the qubit in the lower right corner.
        These errors are propagate and change the X gauge outcomes
        and the final outcomes if measured in the X basis.

        We return a new list of X gauge outcomes and final outcomes.
        """
        frame = self.n * [0]  # Z error frame
        zidx = 0  # index z measurement cycles
        xidx = 0  # index x measurement cycles
        for _ in range(blocks):
            for rs in round_schedule:
                if rs == "z":
                    # Examine the left/right flags and update frame
                    for j, zg in enumerate(self.css_z_gauge_ops):
                        # Only consider weight 4 operators
                        if len(zg) == 4:
                            if (
                                left_flag_outcomes[zidx][j] == 1
                                and right_flag_outcomes[zidx][j] == 0
                            ):
                                # upper left qubit
                                qubit = zg[0]
                                frame[qubit] ^= 1
                                logging.debug("frame, cycle %d -> qubit %d", zidx, qubit)
                            if (
                                left_flag_outcomes[zidx][j] == 0
                                and right_flag_outcomes[zidx][j] == 1
                            ):
                                # lower right qubit
                                qubit = zg[3]
                                frame[qubit] ^= 1
                                logging.debug("frame, cycle %d -> qubit %d", zidx, qubit)
                    zidx += 1
                if rs == "x":
                    # Update the X gauge syndromes
                    syn = temp_syndrome(frame, self.css_x_gauge_ops)
                    logging.debug("frame syndrome, cycle %d -> %s", xidx, syn)
                    block = x_gauge_outcomes[xidx]
                    block = list(u ^ v for u, v in zip(block, syn))
                    x_gauge_outcomes[xidx] = block
                    logging.debug("x gauge update, cycle %d -> %s", xidx, block)
                    xidx += 1
        # Update the final outcomes if X basis measurement
        if self.basis == "x":
            final_outcomes = list(u ^ v for u, v in zip(final_outcomes, frame))
        return x_gauge_outcomes, final_outcomes
