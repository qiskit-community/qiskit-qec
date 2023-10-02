# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Decoder for Spanning Errors on the the ClassicFibonacciCode"""
import copy
import logging
import math
from typing import Dict, Optional, Set, Tuple

import numpy as np
import pymatching as pm
import rustworkx as rx

from qiskit_qec.exceptions import QiskitQECError

logger = logging.getLogger(__name__)


class ClassicFibonacciSpanningErrorDecoder:
    """
    self.board =
    [L*(L//2 - 1)....                   ((L**2)//2) - 1]
     .
     .
     .
     2L
     L
     0 1 2 3 ....                                L - 1]
    """

    def __init__(
        self,
        original_errorword: np.array,
        halt: int = 9,
    ):
        """Decoder for Fibonacci Code (ClassicFibonacciCode)

        This decoder is a Minimum Weight Perfect Matching (MWPM)
        Decoder implemented in Pymatching, modified to work for
        spanning errors on the Fibonacci Code.
        See arXiv:2105.13082 for more information about PyMatching.
        See section III of arXiv:2002.11738 for more information
        about modifying MWPM for spanning errors on the Fibonacci Code.

        L is the number of encoded logical bits

        The code is envisioned as being on a ((L//2, L)) codeboard,
        which is represented by an (L//2, L) numpy array.

        self.board is where the results of the decoder
        are applied during self.decode_fib_code().

        Args:
            original_errorword (np.array): An (L/bl/2, L) np.array
            representing an errored Fibonnaci Code
            halt (int, optional): Number of iterations to run
            before terminating unsuccessfully. Defaults to 9.
        """
        self.num_logical_bits = len(original_errorword[0])  # len
        assert (
            math.log2(self.num_logical_bits) % 1 == 0
        ), "L must be some 2**n where n is an int >= 1"
        self.no_cols = self.num_logical_bits
        self.no_rows = self.num_logical_bits // 2
        self.no_bits = (self.num_logical_bits**2) // 2  # no bits
        self.halt = halt

        # fund_sym
        self.original_errorword = original_errorword

        self.board = copy.deepcopy(self.original_errorword)
        self.board.shape = self.no_bits
        self.fund_symmetry = self._generate_init_symmetry()
        self.fund_symmetry.shape = (self.no_rows, self.no_cols)
        (
            self.fund_stab_parity_check_matrix,
            self.fund_parity_rows_to_faces,
        ) = self.generate_check_matrix_from_faces(self.fund_symmetry)
        self.fund_symmetry.shape = self.no_bits

        self.fund_single_error_syndromes = self.fib_code_fault_enumeration(
            self.fund_stab_parity_check_matrix
        )
        self.hx_mat = self._generate_plus_x_trans_matrix()
        self.hy_mat = self._generate_plus_y_trans_matrix()

        self.fund_stab_faces = copy.deepcopy(self.fund_symmetry)
        self.fund_hori_probe_indx = self.no_bits - (self.num_logical_bits // 2) - 1
        self.fund_verti_probe_indx = self.no_bits - 1
        self.fund_stab_faces.shape = (
            self.num_logical_bits // 2,
            self.num_logical_bits,
        )  # TODO should work
        (
            self.fund_check_matrix,
            self.board2stab,
        ) = self.generate_check_matrix_from_faces(self.fund_stab_faces)
        fund_error_pairs = self.fib_code_fault_enumeration(self.fund_check_matrix)
        self.fund_matching_graph, self.fundstab2node = self.enumerated_errors2matching_graph(
            fund_error_pairs
        )

        self.all_stab_faces = np.ones(
            (self.num_logical_bits // 2, self.num_logical_bits), dtype=int
        )
        (
            self.all_stabs_check_mat,
            self.all_stabs_parity_rows_to_faces,
        ) = self.generate_check_matrix_from_faces(self.all_stab_faces)

        # pymatching
        self.matching_decoder = pm.Matching(self.fund_matching_graph)

        self.has_decoder_run = False  # lets user know whether decoder has already run

        logger.info(" original_errorword is  %s", self.original_errorword)
        logger.info(" error board is code %s", self.board)
        logger.info(" initial symmetry is: %s", self.fund_symmetry)
        logger.info(
            "fund_stab_parity_check_matrix is: %s",
            self.fund_stab_parity_check_matrix,
        )
        logger.info("fund_single_error_syndromes is : %s", self.fund_single_error_syndromes)
        logger.info(" hx %s", self.hx_mat)
        logger.info(" hy is code %s", self.hy_mat)

    def decode_prob(self, syndrome: np.array):
        """Returns whether the horizontal probe edge
        (aka the bottom middle bit of the probe triangle)
        exists an even (0 aka no flip) or odd (1 aka flip)
        number of times in the matching graph && the vertical probe
        """
        res = self.matching_decoder.decode(syndrome)
        return res[self.fund_hori_probe_indx], res[self.fund_verti_probe_indx], res

    def bit_to_rc(self, bit: int):
        """Maps maps from bit index (0 to   (L**2)//2) - 1) to (row, column) indexing

        In bit notation we think of all the bits of the codeword as being in a line:
        So we have bit 0, bit 1, ... all the way until the last bit ((L**2)//2) - 1
        Example:
        [0, 1, 2,  ...................  ((L**2)//2) - 1 ]

        however, we can picture these as being on the
        L//2 by L board.

        In numpy, if we reshape the  ((L**2)//2) - 1  to a L//2 by L array,
        then bit 0 will get mapped to the
        (0th row, 0th column).
        This is called row, column notation
        Bit 1 will get mapped to the 0th row, 1st column.
        ...
        L//2 gets mapped to the 1st row, 0th column
        etc.

        I show here:

        [  [0 1 2  3 ...............................    (L//2) -1]
            L//2
                .
                .
                .
                [(L - 1) * (L//2)  ....................    ((L**2)//2) - 1 ]
            ]

        Args:
            bit (int): Index of the bit in bit notation
        """
        row_len = self.num_logical_bits

        rindx = bit // row_len
        cindx = bit % row_len
        return (rindx, cindx)

    def rc_to_bit(self, row: int, col: int):
        """Maps maps from (row, column) index to bit index (0 to   (L**2)//2) - 1)

        In bit notation we think of all the bits of the codeword as being in a line:
        So we have bit 0, bit 1, ... all the way until the last bit ((L**2)//2) - 1
        Example:
        [0, 1, 2,  ...................  ((L**2)//2) - 1 ]

        however, we can picture these as being on the
        L//2 by L board.

        In numpy, if we reshape the  ((L**2)//2) - 1  to a L//2 by L array,
        then bit 0 will get mapped to the
        (0th row, 0th column).
        This is called row, column notation
        Bit 1 will get mapped to the 0th row, 1st column.
        ...
        L//2 gets mapped to the 1st row, 0th column
        etc.

        I show here:

        [  [0 1 2  3 ...............................    (L//2) -1]
            L//2
                .
                .
                .
                [(L - 1) * (L//2)  ....................    ((L**2)//2) - 1 ]
            ]

        Args:
            row (int): Row index
            col (int): Column index
        """
        bit = (row * self.num_logical_bits) + col
        return bit

    def _generate_plus_x_trans_matrix(self):
        """Takes bit to bit + 1 mod rownum aka shifts bit
        to the right but wraps around its current row"""
        h_mat = np.zeros((self.no_bits, self.no_bits), dtype=int)
        for b in range(self.no_bits):
            new_bit = self.shift_by_x_scalar(b)
            h_mat[new_bit][b] = 1
        return h_mat

    def _generate_plus_y_trans_matrix(self):
        """takes bit to (bit + L) mod (L//2) aka shifts bit the row
        below but very bottom row shifts to be the 0th row"""
        h_mat = np.zeros((self.no_bits, self.no_bits), dtype=int)
        for b in range(self.no_bits):
            new_bit = self.shift_by_y_scalar(b)
            h_mat[new_bit][b] = 1
        return h_mat

    def shift_by_x(self, bitarr: np.array, power: int = 1) -> np.array:
        """shifts every entry in array matrix right by power w/ wrap around

        Args:
            bitarr (np.array): original code array that serves as a starting point for the shift
            Please note, bitarr itself is *not* shifted
            power (int, optional): Number of times to shift. Defaults to 1.

        Returns:
            np.array: shifted array
        """

        power = power % self.num_logical_bits
        hx = np.linalg.matrix_power(self.hx_mat, power)
        sol = np.matmul(hx, bitarr)
        sol = sol.astype(int)
        return sol

    def shift_by_y(self, bitarr: np.array, power: int = 1):
        """shifts every entry in array matrix down by power w/ wrap around

        Args:
            bitarr (np.array): original code array that serves as a starting point for the shift
            Please note, bitarr itself is *not* shifted
            power (int, optional): Number of times to shift. Defaults to 1.

        Returns:
            np.array: shifted array
        """

        power = power % (self.num_logical_bits // 2)
        hy = np.linalg.matrix_power(self.hy_mat, power)
        sol = np.matmul(hy, bitarr)
        sol = sol.astype(int)
        return sol

    def shift_by_y_scalar(self, bit: int, shift_no: int = 1):
        """shifts entry in array matrix down by 1 w/ wrap around

        Args:
            bit (int): bit to be shifted
            shift_no (int, optional): Number of times to shift. Defaults to 1.

        Returns:
            int: new location of the bit after shift
        """

        new_bit = bit
        for _ in range(shift_no):
            new_bit = (new_bit + self.num_logical_bits) % self.no_bits
        return new_bit

    def shift_by_x_scalar(self, bit: int, shift_no: int = 1):
        """shifts entry in array matrix right by 1 w/ wrap around
        Args:
            bit (int): bit to be shifted
            shift_no (int, optional): Number of times to shift. Defaults to 1.

        Returns:
            int: new location of the bit after shift
        """
        new_bit = bit
        for _ in range(shift_no):
            new_bit = ((new_bit + 1) % self.num_logical_bits) + (
                (new_bit // self.num_logical_bits) * (self.num_logical_bits)
            )
        return new_bit

    def shift_parity_mat_by_y(self, parity_mat: np.array, power: int = 1):
        """Shifts parity matrix vertically by power.
        Please note this *modifies* the parity_mat.

        Args:
            parity_mat (np.array): parity matrix to be shifted
            power (int, optional): Number of times to shift. Defaults to 1.

        Returns:
            np.array: original parity_mat
        """

        for row, _ in enumerate(parity_mat):
            parity_mat[row] = self.shift_by_y(parity_mat[row], power=power)
        return parity_mat

    def shift_parity_mat_by_x(self, parity_mat: np.array, power: int = 1):
        """Shifts parity matrix horizontall by power.
        Please note this *modifies* the parity_mat.

        Args:
            parity_mat (np.array): parity matrix to be shifted
            power (int, optional): Number of times to shift. Defaults to 1.

        Returns:
            np.array: original parity_mat
        """
        for row, _ in enumerate(parity_mat):
            parity_mat[row] = self.shift_by_x(parity_mat[row], power=power)
        return parity_mat

    def _calc_syndrome(self, check_matr: np.array, board: np.array = None) -> np.array:
        """Calculates the syndrome of the check_matrix on the board.

        Args:
            check_matr (np.array): check_matrix for calculating the syndrome
            board (np.array, optional): board for calculating the syndrome.
            If None, defaults to self.board. Defaults to None

        Returns:
            np.array: syndrome results
        """
        if board is None:
            board = self.board
        sol = np.matmul(check_matr, board) % 2
        sol = sol.astype(int)
        return sol

    def _generate_init_symmetry(self, start_arr: Optional[np.array] = None) -> np.array:
        """Generates initial symmetry

        Args:
            start_arr (np.array, optional): Optional array that
            serves as the foundation for building up the rest of the codeword
            using cellular automaton update rules, as described in
            section II, formula 2 of arXiv:2002.11738. Defaults to None.

        Raises:
            ArgumentError: Validate start_arr

        Returns:
            np.array: Initial symmetry of the code
        """
        if start_arr and sum(start_arr) != 1:
            raise QiskitQECError(
                "Can only have a single 1 in start_arr."
                + f"All else should be 0 but you have: {start_arr}"
            )
        # fund symmetries start from the top instead of the bottom because numpy
        rect_board = np.zeros((self.num_logical_bits // 2, self.num_logical_bits), dtype=int)
        if start_arr is None:
            start_arr = np.zeros(self.num_logical_bits, dtype=int)
            start_arr[(self.num_logical_bits // 2) - 1] = 1
        rect_board[0] = start_arr
        for row in range(1, self.num_logical_bits // 2):
            for bit in range(self.num_logical_bits):
                new_val = (
                    rect_board[row - 1][(bit - 1) % self.num_logical_bits]
                    ^ rect_board[row - 1][(bit) % self.num_logical_bits]
                    ^ rect_board[row - 1][(bit + 1) % self.num_logical_bits]
                )
                rect_board[row][bit] = new_val
        return rect_board

    def generate_check_matrix_from_faces(self, stab_faces: np.array):
        """Use stabilizer faces to make check matrix

        Args:
            stab_faces (np.array): Stabilizer faces

        Returns:
            np.array: parity matrix
            Dict[int, int]: mapping from
            parity_check_faces to bit notation,
                A dictionary, {x:y}, meaning
                that row x of the parity check matrix
                has face centered at bit b

        """

        stab_row2bit_face = {}
        parity_mat = np.empty((0, self.no_bits), int)

        stab_row = 0
        for row in range(self.no_rows):
            for col in range(self.no_cols):
                if stab_faces[row][col] == 1:
                    a = self.rc_to_bit(row, col)
                    b = self.rc_to_bit((row) % self.no_rows, (col - 1) % self.no_cols)
                    c = self.rc_to_bit((row) % self.no_rows, (col + 1) % self.no_cols)
                    d = self.rc_to_bit(
                        (row - 1) % self.no_rows, col % self.no_cols
                    )  # changed to point the other direction
                    new_stab = np.array([0] * self.no_bits)
                    new_stab[a] = 1
                    new_stab[b] = 1
                    new_stab[c] = 1
                    new_stab[d] = 1
                    parity_mat = np.append(parity_mat, [new_stab], axis=0)
                    stab_row2bit_face[stab_row] = b
                    stab_row += 1

        return parity_mat, stab_row2bit_face

    def fib_code_fault_enumeration(
        self, parity_check_matrix: np.array, no_bits: int = None
    ) -> Set[Tuple]:
        """Return a set containing every possible single
        bit flip error and what stabilizers they light up

        Args:
            parity_check_matrix (np.array): parity check matrix,
            no_bits (int): Number of bits in code array

        Returns:
            Set[Tuple]: A set containing every possible single
            bit flip error and what stabilizers they light up.
            Of the form: {(stab_face1, stab_face2, errorbit0),
            (stab_face1, stab_face2, errorbit1), ...}

        """
        if no_bits is None:
            no_bits = self.no_bits
        error_pairs = set()  # (stab_face1, stab_face2, errorbit)
        single_error = np.zeros(no_bits, dtype=int)

        for b in range(no_bits):
            if no_bits > 10 and b % (no_bits // 10) == 0:
                logger.debug("on bit: %s and error set looks like: %s", b, error_pairs)

            ## set up new single error
            # clear prev_bit
            prev_bit = (b - 1) % no_bits
            single_error[prev_bit] = 0
            # set new error
            single_error[b] = 1

            ## what do it light?
            lighted = self._calc_syndrome(parity_check_matrix, single_error)
            stabs = (lighted == 1).nonzero()[0]

            if len(stabs) % 2 != 0:
                emsg = (
                    f"Minor panic. Error on  bit {b}"
                    + f"causes a BAD syndrome: {stabs} for lighted: {lighted}"
                )
                logger.error(
                    emsg
                )  # TODO just do this via inspection on 1s per column in stab parity check matrix
                raise Exception(emsg)

            if len(stabs) > 0:
                for indx in range(0, len(stabs), 2):
                    error_pairs.add((stabs[indx], stabs[indx + 1], b))

        return error_pairs

    def enumerated_errors2matching_graph(
        self, error_graphs: Set[Tuple]
    ) -> (rx.PyGraph, Dict[int, int]):
        """Turn stabilizer pairs lit by a single error
        into edges in the matching graph

        Args:
            error_graphs (Set[Tuple]): A set containing
            every possible single bit flip error and
            what stabilizers they light up.
            Of the form: {(stab_face1, stab_face2, errorbit0),
            (stab_face1, stab_face2, errorbit1), ...}

        Returns:
            rx.PyGraph: Matching Graph for PyMatching
            Dict[int, int]: Mapping from stabilizer id to graph node id

        """
        stab2node = {}
        graph = rx.PyGraph()

        def add_to_graph(stabid):
            if stabid not in stab2node:
                nodeid = graph.add_node({"element": stabid})
                stab2node[stabid] = nodeid
            return stab2node[stabid]

        for stab0, stab1, fund_e in error_graphs:
            graph_node_n0 = add_to_graph(stab0)
            graph_node_n1 = add_to_graph(stab1)
            graph.add_edge(graph_node_n0, graph_node_n1, {"fault_ids": {fund_e}})

        return graph, stab2node

    def decode_fib_code(self, force: Optional[bool] = False):
        """Run decoder on self.board

        Args:
            force (bool, optional): If true, will run
            the decoder, even if it has already been
            run on self.board
            Defaults to False.

        Raises:
            QiskitQECError: Raises an exception if
            decoder has already been run and force=False
        """
        if not force and self.has_decoder_run:
            raise QiskitQECError("In order to rerun decoder, you must pass in force=True")

        self.has_decoder_run = True
        # generate graphs and mappings
        h_correction = np.zeros(self.no_bits, dtype=int)
        v_correction = np.zeros(self.no_bits, dtype=int)
        parity_check_matrix = copy.deepcopy(self.fund_check_matrix)
        hori_probe_indx = self.fund_hori_probe_indx
        verti_probe_indx = self.fund_verti_probe_indx

        cur_all_syndrome = prev_all_syndrome = (
            (self._calc_syndrome(self.all_stabs_check_mat, self.board) == 1).sum(),
        )
        start_flag = True
        meta_round_count = 0
        round_count = 0
        self.fund_stab_faces.shape = self.no_bits
        while (
            (cur_all_syndrome < prev_all_syndrome or start_flag)
            and cur_all_syndrome != 0
            and meta_round_count < self.halt
        ):
            start_flag = False
            prev_all_syndrome = cur_all_syndrome

            for _ in range(self.num_logical_bits // 2):  # will wrap around to all bits
                parity_check_matrix = self.shift_parity_mat_by_y(parity_check_matrix)
                self.fund_stab_faces = self.shift_by_y(self.fund_stab_faces)
                hori_probe_indx = self.shift_by_y_scalar(hori_probe_indx)
                verti_probe_indx = self.shift_by_y_scalar(verti_probe_indx)

                for _ in range(self.num_logical_bits):
                    parity_check_matrix = self.shift_parity_mat_by_x(parity_check_matrix)
                    self.fund_stab_faces = self.shift_by_x(self.fund_stab_faces)
                    hori_probe_indx = self.shift_by_x_scalar(hori_probe_indx)
                    verti_probe_indx = self.shift_by_x_scalar(verti_probe_indx)

                    self.fund_stab_faces.shape = (
                        self.num_logical_bits // 2,
                        self.num_logical_bits,
                    )
                    self.fund_stab_faces.shape = self.no_bits

                    cur_syndrome = self._calc_syndrome(parity_check_matrix)
                    # convert syndrome to node
                    cur_node_syndrome = [0] * len(cur_syndrome)
                    for stabindx, value in enumerate(cur_syndrome):
                        nodeindx = self.fundstab2node[stabindx]
                        cur_node_syndrome[nodeindx] = value  # TODO is right?
                    hcorval, vcorval, res = self.decode_prob(cur_node_syndrome)

                    h_correction[hori_probe_indx] = hcorval
                    v_correction[verti_probe_indx] = vcorval

                    round_count += 1

                    logger.debug("current fund stabilizer faces:\n%s", self.fund_stab_faces)
                    logger.debug("current_parity_check_mat:\n%s", parity_check_matrix)
                    logger.debug("cur-syndrome-symm: %s", cur_syndrome)
                    logger.debug("res                             is: %s", res)
                    logger.debug("hcorval: %s\nvcorval:%s", hcorval, vcorval)
                    logger.debug(
                        "hori probd inex: %s, verti_probe_inx: %s",
                        hori_probe_indx,
                        verti_probe_indx,
                    )
                    logger.debug("h_corr: %s\nv_corr:%s", h_correction, v_correction)

            logger.info("Meta-Round %s:", meta_round_count)
            logger.info("h_correction: %s\nv_correction:%s", h_correction, v_correction)

            meta_round_count += 1
            d_correction = h_correction * v_correction
            hboard = self.board ^ h_correction  # apply correction
            vboard = self.board ^ v_correction
            dboard = self.board ^ d_correction

            hcorsynd = [
                (self._calc_syndrome(self.all_stabs_check_mat, hboard) == 1).sum(),
                hboard,
                "hori",
            ]
            vcorsynd = [
                (self._calc_syndrome(self.all_stabs_check_mat, vboard) == 1).sum(),
                vboard,
                "verti",
            ]
            dcorsynd = [
                (self._calc_syndrome(self.all_stabs_check_mat, dboard) == 1).sum(),
                dboard,
                "dcor",
            ]

            opts = [hcorsynd, vcorsynd, dcorsynd]

            winner = min(opts, key=lambda x: x[0])
            cur_all_syndrome = winner[0]
            self.board = winner[1]  # update board to best one
            logger.info("Updated board is: \n%s", self.board)

        logger.info("decoding finished")
