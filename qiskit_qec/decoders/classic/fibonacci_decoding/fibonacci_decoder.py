import copy
import logging
import math
from argparse import ArgumentError

import numpy as np
import rustworkx as rx

from qiskit_qec.decoders.classic.fibonacci_decoding.fibonacci_decoding_graph import (
    FibonacciDecoderGraph,
)

logger = logging.getLogger(__name__)


class ClassicFibDecoder:
    """
    self.board =
    [L*(L//2 - 1)....          ((L**2)//2) - 1]
     .
     .
     .
     2L
     L
     0 1 2 3 ....                                L - 1]
    """

    def __init__(
        self,
        original_errorword,
        halt=9,
        name="",
    ):

        self.L = len(original_errorword[0])  # len
        assert math.log2(self.L) % 1 == 0, "L must be some 2**n where n is an int >= 1"
        self.no_cols = self.L
        self.no_rows = self.L // 2
        self.no_bits = (self.L**2) // 2  # no bits
        self.halt = halt
        self.name = name

        # fund_sym
        self.original_errorword = original_errorword

        self.board = copy.deepcopy(self.original_errorword)
        self.board.shape = self.no_bits
        self.fundamental_symmetry = self._generate_init_symmetry()
        self.fundamental_symmetry.shape = (self.no_rows, self.no_cols)
        (
            self.fundamental_stabilizer_parity_check_matrix,
            self.fundamental_parity_rows_to_faces,
        ) = self.generate_check_matrix_from_faces(self.fundamental_symmetry)
        self.fundamental_symmetry.shape = self.no_bits

        self.fundamental_single_error_syndromes = self.generate_all_possible_error_syndromes(
            self.fundamental_stabilizer_parity_check_matrix
        )
        self.Hx = self._generate_plus_x_trans_matrix()
        self.Hy = self._generate_plus_y_trans_matrix()

        self.fundamental_stab_faces = copy.deepcopy(self.fundamental_symmetry)
        self.fundamental_hori_probe_indx = self.no_bits - (self.L // 2) - 1
        self.fundamental_verti_probe_indx = self.no_bits - 1
        self.fundamental_stab_faces.shape = (self.L // 2, self.L)  # TODO should work
        (
            self.fundamental_check_matrix,
            self.board2stab,
        ) = self.generate_check_matrix_from_faces(self.fundamental_stab_faces)
        fund_error_pairs = self.generate_all_possible_error_syndromes(self.fundamental_check_matrix)
        self.fund_matching_graph, self.fundstab2node = self.error_pairs2graph(fund_error_pairs)

        # Develop a generic decoder w/ fundamental_hori_probe_indx & fundamental_verti_probe_indx as anchors for indexing
        self.decoder = FibonacciDecoderGraph(
            self.fund_matching_graph,
            self.fundamental_hori_probe_indx,
            self.fundamental_verti_probe_indx,
            self.fundstab2node,
        )
        self.all_stab_faces = np.ones((self.L // 2, self.L), dtype=int)
        (
            self.all_stabs_check_mat,
            self.all_stabs_parity_rows_to_faces,
        ) = self.generate_check_matrix_from_faces(self.all_stab_faces)

        logger.info(f" original_errorword is  {self.original_errorword}")
        logger.info(f" error board is code {self.board}")
        logger.info(f" initial symmetry is: {self.fundamental_symmetry}")
        logger.info(
            f"fundamental_stabilizer_parity_check_matrix is : {self.fundamental_stabilizer_parity_check_matrix}"
        )
        logger.info(
            f"fundamental_single_error_syndromes is : {self.fundamental_single_error_syndromes}"
        )
        logger.info(f" Hx {self.Hx}")
        logger.info(f" Hy is code {self.Hy}")

    def bit_to_rc(self, bit):
        """
         This functions maps from bit index (0 to   (L**2)//2) - 1) to row, column mapping
        # DOESNT yet handle bits that are too big


        in bit notation we think of all the bits being in a line:
        So we have bit 0, bit 1, ... all the way until the last bit ((L**2)//2) - 1
        [0, 1, 2,  ...................  ((L**2)//2) - 1 ]

        However, we can picture these as being on the
        L//2 by L  board.

        In numpy, if we reshape the  ((L**2)//2) - 1  to a L//2 by L array,
        then bit 0 will get mapped to the
        0th row, 0th column.
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

        """
        row_len = self.L

        rindx = bit // row_len
        cindx = bit % row_len
        return (rindx, cindx)

    def rc_to_bit(self, row, col):
        """
        # DOESNT yet handle row/columns that are too big
        This functions maps from (row, column) indexing to bit  (0 to   (L**2)//2) - 1) indexing

        in bit notation we think of all the bits being in a line:
        So we have bit 0, bit 1, ... all the way until the last bit ((L**2)//2) - 1
        [0, 1, 2,  ...................  ((L**2)//2) - 1 ]

        However, we can picture these as being on the
        L//2 by L  board.

        In numpy, if we reshape the  ((L**2)//2) - 1  to a L//2 by L array,
        then bit 0 will get mapped to the
        0th row, 0th column.
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

        """
        bit = (row * self.L) + col
        return bit

    def _generate_plus_x_trans_matrix(self):
        "Takes bit to bit + 1 mod rownum aka shifts bit to the right but wraps around its current row"
        H = np.zeros((self.no_bits, self.no_bits), dtype=int)
        for b in range(self.no_bits):
            new_bit = self.shift_by_x_scalar(b)
            H[new_bit][b] = 1
        return H

    def _generate_plus_y_trans_matrix(self):
        "takes bit to bit + L mod L//2 aka shifts bit the row below but very bottom row shifts to be the 0th row"
        H = np.zeros((self.no_bits, self.no_bits), dtype=int)
        for b in range(self.no_bits):
            new_bit = self.shift_by_y_scalar(b)
            H[new_bit][b] = 1
        return H

    def shift_by_x(self, bitarr, power=1):
        # shifts every entry in board matrix right by 1 w/ wrap around
        power = power % self.L
        Hx = np.linalg.matrix_power(self.Hx, power)
        sol = np.matmul(Hx, bitarr)
        sol = sol.astype(int)
        return sol

    def shift_by_y(self, bitarr, power=1):
        # shifts every entry in board matrix down by 1 w/ wrap around
        power = power % (self.L // 2)
        Hy = np.linalg.matrix_power(self.Hy, power)
        sol = np.matmul(Hy, bitarr)
        sol = sol.astype(int)
        return sol

    def shift_by_y_scalar(self, bit, shift_no=1):
        # shifts entry in board matrix down by 1 w/ wrap around
        new_bit = bit
        for _ in range(shift_no):
            new_bit = (new_bit + self.L) % self.no_bits
        return new_bit

    def shift_by_x_scalar(self, bit, shift_no=1):
        # shifts entry in board matrix right by 1 w/ wrap around
        new_bit = bit
        for _ in range(shift_no):
            new_bit = ((new_bit + 1) % self.L) + ((new_bit // self.L) * (self.L))
        return new_bit

    def shift_parity_mat_by_y(self, parity_mat, power=1):
        "EDITS PARITY MAT"
        for row in range(len(parity_mat)):
            parity_mat[row] = self.shift_by_y(parity_mat[row], power=power)
        return parity_mat

    def shift_parity_mat_by_x(self, parity_mat, power=1):
        "EDITS PARITY MAT"
        for row in range(len(parity_mat)):
            parity_mat[row] = self.shift_by_x(parity_mat[row], power=power)
        return parity_mat

    def _calc_syndrome(self, check_matr, board=None):
        if board is None:
            board = self.board
        #  % 2 # TODO numpy almost certainly has a way of efficiently dealing w binary matrices -- figure that out
        sol = np.matmul(check_matr, board) % 2
        sol = sol.astype(int)
        return sol

    def _generate_init_symmetry(self, start_arr=None):
        if start_arr and sum(start_arr) != 1:
            raise ArgumentError(
                f"Can only have a single 1 in start_arr. All else should be 0 but you have: {start_arr}"
            )
        # fundamental symmetries start from the top instead of the bottom because numpy
        rect_board = np.zeros((self.L // 2, self.L), dtype=int)
        if start_arr is None:
            start_arr = np.zeros(self.L, dtype=int)
            start_arr[(self.L // 2) - 1] = 1
        rect_board[0] = start_arr
        for row in range(1, self.L // 2):
            for bit in range(self.L):
                new_val = (
                    rect_board[row - 1][(bit - 1) % self.L]
                    ^ rect_board[row - 1][(bit) % self.L]
                    ^ rect_board[row - 1][(bit + 1) % self.L]
                )
                rect_board[row][bit] = new_val
        return rect_board

    def generate_check_matrix_from_faces(self, stab_faces):
        """x
        STABS:           xxx
        """
        # TODO slow because of all the shape changing
        """REQUIRES L//2 x L """
        # AHHHH
        # create parity check matrix for each
        # np.append([[1, 2, 3], [4, 5, 6]], [[7, 8, 9]], axis=0)

        faces_to_stabs_rows = (
            {}
        )  # {x:y} where row x of the parity check matrix has face centered at bit b
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
                    faces_to_stabs_rows[stab_row] = b
                    stab_row += 1

        return parity_mat, faces_to_stabs_rows

    def generate_all_possible_error_syndromes(self, parity_check_matrix, no_bits=None):

        # "there's gotta be a smarter way to do this "
        if no_bits is None:
            no_bits = self.no_bits
        error_pairs = set()  # (stab_face1, stab_face2, errorbit)
        single_error = np.zeros(no_bits, dtype=int)

        for b in range(no_bits):
            if no_bits > 10 and b % (no_bits // 10) == 0:
                logger.debug(f"on bit: {b} and error set looks like: {error_pairs}")

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
                emsg = f"Minor panic. Error on  bit {b} causes a BAD syndrome: {stabs} for lighted: {lighted}"
                logger.error(
                    emsg
                )  # TODO just do this via inspection on 1s per column in stab parity check matrix
                raise Exception(emsg)

            if len(stabs) > 0:
                for indx in range(0, len(stabs), 2):
                    error_pairs.add((stabs[indx], stabs[indx + 1], b))

        return error_pairs

    def error_pairs2graph(self, error_graphs):
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

    def decode_fib_code(self):
        # generate graphs and mappings

        h_correction = np.zeros(self.no_bits, dtype=int)
        v_correction = np.zeros(self.no_bits, dtype=int)
        parity_check_matrix = copy.deepcopy(self.fundamental_check_matrix)
        hori_probe_indx = self.fundamental_hori_probe_indx
        verti_probe_indx = self.fundamental_verti_probe_indx

        cur_all_syndrome = prev_all_syndrome = (
            (self._calc_syndrome(self.all_stabs_check_mat, self.board) == 1).sum(),
        )
        start_flag = True
        meta_round_count = 0
        round_count = 0
        self.fundamental_stab_faces.shape = self.no_bits
        while (
            (cur_all_syndrome < prev_all_syndrome or start_flag)
            and cur_all_syndrome != 0
            and meta_round_count < self.halt
        ):
            start_flag = False
            prev_all_syndrome = cur_all_syndrome

            for y_offset in range(self.L // 2):  # will wrap around to all bits
                parity_check_matrix = self.shift_parity_mat_by_y(parity_check_matrix)
                self.fundamental_stab_faces = self.shift_by_y(self.fundamental_stab_faces)
                hori_probe_indx = self.shift_by_y_scalar(hori_probe_indx)
                verti_probe_indx = self.shift_by_y_scalar(verti_probe_indx)

                for x_offset in range(self.L):
                    parity_check_matrix = self.shift_parity_mat_by_x(parity_check_matrix)
                    self.fundamental_stab_faces = self.shift_by_x(self.fundamental_stab_faces)
                    hori_probe_indx = self.shift_by_x_scalar(hori_probe_indx)
                    verti_probe_indx = self.shift_by_x_scalar(verti_probe_indx)

                    self.fundamental_stab_faces.shape = (self.L // 2, self.L)
                    self.fundamental_stab_faces.shape = self.no_bits

                    cur_syndrome = self._calc_syndrome(parity_check_matrix)
                    # convert syndrome to node
                    cur_node_syndrome = [0] * len(cur_syndrome)
                    for stabindx, value in enumerate(cur_syndrome):
                        nodeindx = self.fundstab2node[stabindx]
                        cur_node_syndrome[nodeindx] = value  # TODO is right?
                    hcorval, vcorval, res = self.decoder.decode_prob(cur_node_syndrome)

                    h_correction[hori_probe_indx] = hcorval
                    v_correction[verti_probe_indx] = vcorval

                    round_count += 1

                    logger.debug(f"PROBs:current fundy:\n{self.fundamental_stab_faces}")
                    logger.debug(f"current_parity_check_mat:\n{parity_check_matrix}")
                    logger.debug(f"cur-syndrome-symm: {cur_syndrome}")
                    logger.debug(f"res                             is: {res}")
                    logger.debug(f"hcorval: {hcorval}\nvcorval:{vcorval}")
                    logger.debug(
                        f"hori probd inex: {hori_probe_indx}, verti_probe_inx: {verti_probe_indx}"
                    )
                    logger.debug(f"h_corr: {h_correction}\nv_corr:{v_correction}")

            logger.info(f"Meta-Round {meta_round_count}:")
            logger.info(f"h_correction: {h_correction}\nv_correction:{v_correction}")

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
            logger.info(f"Updated board is: \n{self.board}")

        logger.info("decoding finished")
