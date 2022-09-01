"""Object to construct quantum circuits for the RSSC."""

import copy
import logging
from tokenize import group
from typing import List, Tuple

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.circuit.library import IGate
from qiskit_qec.codes.rotated_surface_code import RSSC


class RSSCCircuit:
    """Create quantum circuits for syndrome measurements.

    Specialized to rotated subsystem surface code.
    """

    def __init__(
        self,
        rssc,
        idles:bool,
        basis:str,
        rounds:int,
        schedule,
        barriers:bool,
        group_meas,
        initial_state:str,
        round_schedule:str,
        logical_paulis:str,
        distinct_measurement_idle:bool,
        num_initialize:int=None,
    ):
        """Create an object associated to a RSSC.

        rssc = RSSC code object
        config = Config object
        """
        self.code = rssc
        self.schedule = schedule
        if self.schedule == "heavy-hex":
            self.total_ancilla, self.z_ancilla_indices = self._hex_ancillas(rssc)
        elif self.schedule == "higgott-breuckmann":
            self.total_ancilla = 2 * (len(rssc.x_gauges) + len(rssc.z_gauges))
        self.total_qubits = rssc.n + self.total_ancilla
        self.qreg = QuantumRegister(self.total_qubits)
        self.imlabel = "idm"  # label for idle during measurement

        self.barriers = barriers
        self.idles = idles 
        self.distinct_measurement_idle = distinct_measurement_idle
        if rounds < 1:
            raise Exception("expected positive integer rounds")
        self.rounds= rounds  
        if set(round_schedule) > set("xz"):
            raise Exception("expected round schedule of 'x', 'z' chars")
        self.round_schedule = round_schedule 
        if not (basis == "x" or basis == "z"):
            raise Exception("expected basis to be 'x' or 'z'")
        self.basis = basis 
        if not (initial_state == "+" or initial_state == "-"):
            raise Exception("expected initial state '+' or '-'")
        self.initial_state = initial_state
        if set(logical_paulis) > set("ixyz"):
            raise Exception("expected 'i', 'x', 'y', 'z' logical paulis")
        if len(logical_paulis) > 0 and len(logical_paulis) != len(round_schedule):
            raise Exception("len(logical paulis) != len(round schedule)")
        self.logical_paulis = logical_paulis
        self.num_initialize = num_initialize
        self.group_meas = group_meas

    def _hex_ancillas(self, code):
        """Get ancilla count and indices for Z gauge measurements.

        code = rssc code object
        Return:
            ancilla_count with total number of ancillas.
            ancilla_indices list whose jth entry is the list of
            indices for the ancillas of the jth z gauge operator.
            This list has length 2 for bulk Z gauges and length
            3 for boundary Z gauges.
        """
        ancilla_count = 0  # total ancillas
        ancilla_indices = []  # shared Z gauge ancilla indices
        for j in range(len(code.z_gauges)):
            if -1 in code.z_gauges[j]:  # boundary Z gauge
                # Hex layout has 3 extra ancillas per boundary Z gauge
                start = ancilla_count + len(code.x_gauges)
                ancilla_indices.append([start, start + 1, start + 2])
                ancilla_count += 3
                logging.debug(
                    "%s: boundary, ancillas %d, %d, %d"
                    % (code.z_gauges[j], start, start + 1, start + 2)
                )
            else:  # bulk Z gauge
                # These borrow ancillas from the adjacent X gauges
                if code.z_orientations[j] == 0:  # left triangle
                    for k in range(len(code.x_gauges)):
                        # qubits at indices 0 and 2 are shared
                        # with adjacent down X triangle
                        if (
                            code.z_gauges[j][0] in code.x_gauges[k]
                            and code.z_gauges[j][2] in code.x_gauges[k]
                        ):
                            down = k
                        # qubits at indices 1 and 2 are shared
                        # with adjacent up X triangle
                        if (
                            code.z_gauges[j][1] in code.x_gauges[k]
                            and code.z_gauges[j][2] in code.x_gauges[k]
                        ):
                            up = k
                    logging.debug("%s: left, bulk, ancillas %d, %d" % (code.z_gauges[j], down, up))
                elif code.z_orientations[j] == 1:  # right triangle
                    for k in range(len(code.x_gauges)):
                        # qubits at indices 0 and 1 are shared
                        # with adjacent down X triangle
                        if (
                            code.z_gauges[j][0] in code.x_gauges[k]
                            and code.z_gauges[j][1] in code.x_gauges[k]
                        ):
                            down = k
                        # qubits at indices 1 and 2 are shared
                        # with adjacent up X triangle
                        if (
                            code.z_gauges[j][1] in code.x_gauges[k]
                            and code.z_gauges[j][2] in code.x_gauges[k]
                        ):
                            up = k
                    logging.debug("%s: right, bulk, ancillas %d, %d" % (code.z_gauges[j], down, up))
                ancilla_indices.append([down, up])
        ancilla_count += len(code.x_gauges)
        logging.info("ancilla_count = %d" % ancilla_count)
        logging.info("ancilla_indices = %s" % ancilla_indices)
        return ancilla_count, ancilla_indices

    def _hb_round(self, leaving, entering, circ, creg, cbits):
        """Measure the requested gauge operators.
        Use the circuits described by Higgott-Breuckmann Fig 10ab.
        Append the gates to the input circuit circ. Measure results
        into the bits at the indices in cbits (list) in creg.
        The leaving and entering strings belong to the set
        {"", "zx", "zz", "xx"} and are used to determine what
        circuit to apply this round. For example, if leaving=""
        and entering="zx" then we do the ancilla qubit
        initialization and subsequent CNOT gates in Fig 10a,
        but do not apply any operations before the corresponding
        ancilla is initialized.
        The cbits are ordered as given in the leaving string.
        For example, if leaving="" then the cbit list is empty,
        but if leaving="zx" then the cbit list contains the
        Z gauge outcomes followed by X gauge outcomes.
        X gauge operators are lists of qubit indices in clockwise
        order around the face starting with the left-most qubit.
        Z gauge operators are lists of qubits indices in clockwise
        order around the face starting with the top-most qubit.
        There are three leaving schedules given as lists of times
        when corresponding qubits interact with ancillary qubits.
        In these schedules we get the outcomes of the indicated
        gauge operators.
        Leaving ZX:
        X down {1, 0, -1}, measured in step 2 -- ancilla 0
        X up {1, 0, 2}, measured in step 3 -- ancilla 0
        Z right {-1, -1, -1}, measured in step 0 -- ancilla 0
        Z left {-1, 0, -1}, measured in step 1 -- ancilla 0
        0: A0_X_down -> X_down[1]
           A0_X_up -> X_up[1]
           measure A0_Z_right
           Z_left[1] -> A0_Z_left
        1: A0_X_down -> X_down[0]
           A0_X_up -> X_up[0]
           measure A0_Z_left
        2: measure A0_X_down
           A0_X_up -> X_up[2]
        3: measure A0_X_up
        Leaving ZZ:
        Z right {1, -1, 0}, measured in step 2 -- ancilla 1
                {-1, -1, -1}, measured in step 0 -- ancilla 0
        Z left {1, 2, 0}, measured in step 3 -- ancilla 1
               {-1, 0, -1}, measured in step 1 -- ancilla 0
        0: Z_right[2] -> A1_Z_right
           measure A0_Z_right
           Z_left[2] -> A1_Z_left
           Z_left[1] -> A0_Z_left
        1: Z_right[0] -> A1_Z_right
           Z_left[0] -> A1_Z_left
           measure A0_Z_left
        2: measure A1_Z_right
           Z_left[1] -> A1_Z_left
        3: measure A1_Z_left
        Leaving XX:
        X down {1, 0, -1}, measured in step 2 -- ancilla 0
               {-1, -1, -1}, measured in step 0 -- ancilla 1
        X up {1, 0, 2}, measured in step 3 -- ancilla 0
             {-1, -1, 0}, measured in step 1 -- ancilla 1
        0: A0_X_down -> X_down[1]
           measure A1_X_down
           A0_X_up -> X_up[1]
           A1_X_up -> X_up[2]
        1: A0_X_down -> X_down[0]
           A0_X_up -> X_up[0]
           measure A1_X_up
        2: measure A0_X_down
           A0_X_up -> X_up[2]
        3: measure A0_X_up
        There are three entering schedules given in the same way.
        These schedules do part of the syndrome measurement circuit
        that initializes the ancillas and starts the round.
        Entering ZX:
        X down {-1, -1, 3}, prepared in step 2 -- ancilla 0
        X up {-1, -1, -1}, prepared in step 3 -- ancilla 0
        Z right {3, 1, 2}, prepared in step 0 -- ancilla 0
        Z left {3, -1, 2}, prepared in step 1 -- ancilla 0
        0: prepare A0_Z_right
        1: Z_right[1] -> A0_Z_right
           prepare A0_Z_left
        2: prepare A0_X_down
           Z_right[2] -> A0_Z_right
           Z_left[2] -> A0_Z_left
        3: A0_X_down -> X_down[2]
           prepare A0_X_up
           Z_right[0] -> A0_Z_right
           Z_left[0] -> A0_Z_left
        Entering ZZ:
        Z right {3, 1, 2}, prepared in step 0 -- ancilla 0
                {-1, 3, -1}, prepared in step 2 -- ancilla 1
        Z left {3, -1, 2}, prepared in step 1 -- ancilla 0
               {-1, -1, -1}, prepared in step 3 -- ancilla 1
        0: prepare A0_Z_right
        1: Z_right[1] -> A0_Z_right
           prepare A0_Z_left
        2: Z_right[2] -> A0_Z_right
           prepare A1_Z_right
           Z_left[2] -> A0_Z_left
        3: Z_right[0] -> A0_Z_right
           Z_right[1] -> A1_Z_right
           Z_left[0] -> A0_Z_left
           prepare A1_Z_left
        Entering XX:
        X down {-1, -1, 3}, prepared in step 2 -- ancilla 0
               {3, 2, 1}, prepared in step 0 -- ancilla 1
        X up {-1, -1, -1}, prepared in step 3 -- ancilla 0
             {3, 2, -1}, prepared in step 1 -- ancilla 1
        0: prepare A1_X_down
        1: A1_X_down -> X_down[2]
           prepare A1_X_up
        2: prepare A0_X_down
           A1_X_down -> X_down[1]
           A1_X_up -> X_up[1]
        3: A0_X_down -> X_down[2]
           A1_X_down -> X_down[0]
           prepare A0_X_up
           A1_X_up -> X_up[0]
        The leaving and entering circuits can be paired in any way
        to give a schedule for a round.
        The ancilla are prepared and measured as needed and are
        not idle.
        We ignore idle errors on the data qubits since these do not
        appear to have been included in HB's analysis. They only
        occur on the boundaries.
        """
        assert self.schedule == "higgott-breuckmann"
        add_barriers = self.barriers
        nxg = len(self.code.x_gauges)
        nzg = len(self.code.z_gauges)
        # Iterate over rounds 0, 1, 2, 3
        for r in range(4):
            # Iterate over X gauge operators
            for i in range(nxg):
                up = True if self.code.x_orientations[i] == 0 else False
                if leaving == "zx":
                    if not up:
                        if r == 0:
                            # A0_X_down -> X_down[1]
                            if self.code.x_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][1]]
                                )
                        if r == 1:
                            # A0_X_down -> X_down[0]
                            if self.code.x_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][0]]
                                )
                        if r == 2:
                            # measure A0_X_down
                            circ.h(self.qreg[self.code.n + i])
                            circ.measure(self.qreg[self.code.n + i], creg[cbits[nzg + i]])
                    if up:
                        if r == 0:
                            # A0_X_up -> X_up[1]
                            if self.code.x_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][1]]
                                )
                        if r == 1:
                            # A0_X_up -> X_up[0]
                            if self.code.x_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][0]]
                                )
                        if r == 2:
                            # A0_X_up -> X_up[2]
                            if self.code.x_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][2]]
                                )
                        if r == 3:
                            # measure A0_X_up
                            circ.h(self.qreg[self.code.n + i])
                            circ.measure(self.qreg[self.code.n + i], creg[cbits[nzg + i]])
                if leaving == "xx":
                    if not up:
                        if r == 0:
                            # A0_X_down -> X_down[1]
                            if self.code.x_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][1]]
                                )
                            # measure A1_X_down
                            circ.h(self.qreg[self.code.n + nxg + i])
                            circ.measure(self.qreg[self.code.n + nxg + i], creg[cbits[i]])
                        if r == 1:
                            # A0_X_down -> X_down[0]
                            if self.code.x_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][0]]
                                )
                        if r == 2:
                            # measure A0_X_down
                            circ.h(self.qreg[self.code.n + i])
                            circ.measure(self.qreg[self.code.n + i], creg[cbits[nxg + i]])
                    if up:
                        if r == 0:
                            # A0_X_up -> X_up[1]
                            if self.code.x_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][1]]
                                )
                            # A1_X_up -> X_up[2]
                            if self.code.x_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + nxg + i],
                                    self.qreg[self.code.x_gauges[i][2]],
                                )
                        if r == 1:
                            # A0_X_up -> X_up[0]
                            if self.code.x_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][0]]
                                )
                            # measure A1_X_up
                            circ.h(self.qreg[self.code.n + nxg + i])
                            circ.measure(self.qreg[self.code.n + nxg + i], creg[cbits[i]])
                        if r == 2:
                            # A0_X_up -> X_up[2]
                            if self.code.x_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][2]]
                                )
                        if r == 3:
                            # measure A0_X_up
                            circ.h(self.qreg[self.code.n + i])
                            circ.measure(self.qreg[self.code.n + i], creg[cbits[nxg + i]])
                if entering == "zx":
                    if not up:
                        if r == 2:
                            # prepare A0_X_down
                            circ.reset(self.qreg[self.code.n + i])
                            circ.h(self.qreg[self.code.n + i])
                        if r == 3:
                            # A0_X_down -> X_down[2]
                            if self.code.x_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][2]]
                                )
                    if up:
                        if r == 3:
                            # prepare A0_X_up
                            circ.reset(self.qreg[self.code.n + i])
                            circ.h(self.qreg[self.code.n + i])
                elif entering == "xx":
                    if not up:
                        if r == 0:
                            # prepare A1_X_down
                            circ.reset(self.qreg[self.code.n + nxg + i])
                            circ.h(self.qreg[self.code.n + nxg + i])
                        if r == 1:
                            # A1_X_down -> X_down[2]
                            if self.code.x_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + nxg + i],
                                    self.qreg[self.code.x_gauges[i][2]],
                                )
                        if r == 2:
                            # prepare A0_X_down
                            circ.reset(self.qreg[self.code.n + i])
                            circ.h(self.qreg[self.code.n + i])
                            # A1_X_down -> X_down[1]
                            if self.code.x_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + nxg + i],
                                    self.qreg[self.code.x_gauges[i][1]],
                                )
                        if r == 3:
                            # A0_X_down -> X_down[2]
                            if self.code.x_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][2]]
                                )
                            # A1_X_down -> X_down[0]
                            if self.code.x_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + nxg + i],
                                    self.qreg[self.code.x_gauges[i][0]],
                                )
                    if up:
                        if r == 1:
                            # prepare A1_X_up
                            circ.reset(self.qreg[self.code.n + nxg + i])
                            circ.h(self.qreg[self.code.n + nxg + i])
                        if r == 2:
                            # A1_X_up -> X_up[1]
                            if self.code.x_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + nxg + i],
                                    self.qreg[self.code.x_gauges[i][1]],
                                )
                        if r == 3:
                            # prepare A0_X_up
                            circ.reset(self.qreg[self.code.n + i])
                            circ.h(self.qreg[self.code.n + i])
                            # A1_X_up -> X_up[0]
                            if self.code.x_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.n + nxg + i],
                                    self.qreg[self.code.x_gauges[i][0]],
                                )
            # Iterate over Z gauge operators
            for i in range(nzg):
                left = True if self.code.z_orientations[i] == 0 else False
                if leaving == "zx":
                    if left:
                        if r == 0:
                            # Z_left[1] -> A0_Z_left
                            if self.code.z_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][1]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                        if r == 1:
                            # measure A0_Z_left
                            circ.measure(self.qreg[self.code.n + 2 * nxg + i], creg[cbits[i]])
                    if not left:
                        if r == 0:
                            # measure A0_Z_right
                            circ.measure(self.qreg[self.code.n + 2 * nxg + i], creg[cbits[i]])
                if leaving == "zz":
                    if left:
                        if r == 0:
                            # Z_left[2] -> A1_Z_left
                            if self.code.z_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][2]],
                                    self.qreg[self.code.n + 2 * nxg + nzg + i],
                                )
                            # Z_left[1] -> A0_Z_left
                            if self.code.z_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][1]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                        if r == 1:
                            # Z_left[0] -> A1_Z_left
                            if self.code.z_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][0]],
                                    self.qreg[self.code.n + 2 * nxg + nzg + i],
                                )
                            # measure A0_Z_left
                            circ.measure(self.qreg[self.code.n + 2 * nxg + i], creg[cbits[i]])
                        if r == 2:
                            # Z_left[1] -> A1_Z_left
                            if self.code.z_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][1]],
                                    self.qreg[self.code.n + 2 * nxg + nzg + i],
                                )
                        if r == 3:
                            # measure A1_Z_left
                            circ.measure(
                                self.qreg[self.code.n + 2 * nxg + nzg + i], creg[cbits[nzg + i]]
                            )
                    if not left:
                        if r == 0:
                            # Z_right[2] -> A1_Z_right
                            if self.code.z_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][2]],
                                    self.qreg[self.code.n + 2 * nxg + nzg + i],
                                )
                            # measure A0_Z_right
                            circ.measure(self.qreg[self.code.n + 2 * nxg + i], creg[cbits[i]])
                        if r == 1:
                            # Z_right[0] -> A1_Z_right
                            if self.code.z_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][0]],
                                    self.qreg[self.code.n + 2 * nxg + nzg + i],
                                )
                        if r == 2:
                            # measure A1_Z_right
                            circ.measure(
                                self.qreg[self.code.n + 2 * nxg + nzg + i], creg[cbits[nzg + i]]
                            )
                if entering == "zx":
                    if left:
                        if r == 1:
                            # prepare A0_Z_left
                            circ.reset(self.qreg[self.code.n + 2 * nxg + i])
                        if r == 2:
                            # Z_left[2] -> A0_Z_left
                            if self.code.z_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][2]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                        if r == 3:
                            # Z_left[0] -> A0_Z_left
                            if self.code.z_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][0]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                    if not left:
                        if r == 0:
                            # prepare A0_Z_right
                            circ.reset(self.qreg[self.code.n + 2 * nxg + i])
                        if r == 1:
                            # Z_right[1] -> A0_Z_right
                            if self.code.z_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][1]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                        if r == 2:
                            # Z_right[2] -> A0_Z_right
                            if self.code.z_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][2]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                        if r == 3:
                            # Z_right[0] -> A0_Z_right
                            if self.code.z_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][0]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                if entering == "zz":
                    if left:
                        if r == 1:
                            # prepare A0_Z_left
                            circ.reset(self.qreg[self.code.n + 2 * nxg + i])
                        if r == 2:
                            # Z_left[2] -> A0_Z_left
                            if self.code.z_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][2]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                        if r == 3:
                            # Z_left[0] -> A0_Z_left
                            if self.code.z_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][0]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                            # prepare A1_Z_left
                            circ.reset(self.qreg[self.code.n + 2 * nxg + nzg + i])
                    if not left:
                        if r == 0:
                            # prepare A0_Z_right
                            circ.reset(self.qreg[self.code.n + 2 * nxg + i])
                        if r == 1:
                            # Z_right[1] -> A0_Z_right
                            if self.code.z_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][1]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                        if r == 2:
                            # Z_right[2] -> A0_Z_right
                            if self.code.z_gauges[i][2] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][2]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                            # prepare A1_Z_right
                            circ.reset(self.qreg[self.code.n + 2 * nxg + nzg + i])
                        if r == 3:
                            # Z_right[0] -> A0_Z_right
                            if self.code.z_gauges[i][0] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][0]],
                                    self.qreg[self.code.n + 2 * nxg + i],
                                )
                            # Z_right[1] -> A1_Z_right
                            if self.code.z_gauges[i][1] != -1:
                                circ.cx(
                                    self.qreg[self.code.z_gauges[i][1]],
                                    self.qreg[self.code.n + 2 * nxg + nzg + i],
                                )
            if add_barriers:
                circ.barrier(self.qreg)

    def _x_gauge_one_round_hex(self, circ, creg, cbits, finalRound=False):
        """Measure all of the X gauge operators.

        Uses a single ancilla per X gauge operator.
        Append the gates to the input circuit circ. Measure results into
        the bits at the indices in cbits (list) in creg.
        """
        assert self.schedule == "heavy-hex"
        # Get options from self.config
        add_barriers = self.barriers
        add_idles = self.idles
        distinct_measurement_idle = self.distinct_measurement_idle
        num_initialize = self.num_initialize
        group_meas = self.group_meas
        # Initialize the ancilla, idle locations added at measurement
        for i in range(len(self.code.x_gauges)):
            if num_initialize is None:
                circ.reset(self.qreg[self.code.n + i])
            circ.h(self.qreg[self.code.n + i])
        if add_barriers:
            circ.barrier(self.qreg)
        # Interact with the ancilla
        for step in range(3):
            idle = list(range(self.code.n))
            for i in range(len(self.code.x_gauges)):
                if self.code.x_gauges[i][step] >= 0:
                    circ.cx(self.qreg[self.code.n + i], self.qreg[self.code.x_gauges[i][step]])
                    # data is involved in a gate this step
                    idle.remove(self.code.x_gauges[i][step])
                else:
                    # ancilla is idle this step
                    idle.append(self.code.n + i)
            if add_idles:
                for i in idle:
                    circ.i(i)
            if add_barriers:
                circ.barrier(self.qreg)
        # Measure the ancilla
        for i in range(len(self.code.x_gauges)):
            circ.h(self.qreg[self.code.n + i])
            if not group_meas:
                circ.measure(self.qreg[self.code.n + i], creg[cbits[i]])
                if not finalRound and num_initialize is not None:
                    circ.reset(self.qreg[self.code.n + i])

        if group_meas:
            if not add_barriers:
                circ.barrier(self.qreg)
            for i in range(len(self.code.x_gauges)):
                circ.measure(self.qreg[self.code.n + i], creg[cbits[i]])
                if not finalRound and num_initialize is not None:
                    circ.reset(self.qreg[self.code.n + i])
        # Data is idle during measurement/reset
        if add_idles:
            for i in range(self.code.n):
                if distinct_measurement_idle:
                    circ.append(IGate(label=self.imlabel), [self.qreg[i]])
                else:
                    circ.i(self.qreg[i])
        if add_barriers:
            circ.barrier(self.qreg)

    def _z_gauge_one_round_hex(self, circ, creg, cbits, finalRound=False):
        """Measure all of the Z gauge operators.

        Reuses the ancillas for X gauge operator.
        Append the gates to the input circuit circ. Measure results into
        the bits at the indices in cbits (list) in creg.
        """
        assert self.schedule == "heavy-hex"
        # Get options from self.config
        add_barriers = self.barriers
        add_idles = self.idles
        distinct_measurement_idle = self.distinct_measurement_idle
        num_initialize = self.num_initialize
        group_meas = self.group_meas
        # Measure 'left' triangles first
        for j in range(len(self.code.z_gauges)):
            if self.code.z_orientations[j] == 0:
                if len(self.z_ancilla_indices[j]) == 3:
                    # left triangle boundary circuit
                    data = copy.copy(self.code.z_gauges[j])
                    data.remove(-1)
                    a0 = self.code.n + self.z_ancilla_indices[j][0]
                    a1 = self.code.n + self.z_ancilla_indices[j][1]
                    a2 = self.code.n + self.z_ancilla_indices[j][2]
                    # timestep 1
                    if add_idles and not distinct_measurement_idle:
                        circ.i(data[0])
                        circ.i(data[1])
                    if num_initialize is None:
                        circ.reset(a2)
                    # timestep 2
                    if add_idles:
                        circ.i(data[0])
                    if num_initialize is None:
                        circ.reset(a0)
                        circ.reset(a1)
                    circ.cx(data[1], a2)
                    # timestep 3
                    circ.cx(data[0], a0)
                    circ.cx(a2, a1)
                    if add_idles:
                        circ.i(data[1])
                    # timestep 4
                    if add_idles:
                        circ.i(data[0])
                    circ.cx(a0, a1)
                    circ.cx(data[1], a2)
                    # timestep 5
                    circ.cx(data[0], a0)
                    circ.measure(a1, creg[cbits[j]])
                    if not finalRound and num_initialize is not None:
                        circ.reset(a1)
                    if add_idles:
                        if distinct_measurement_idle:
                            circ.append(IGate(label=self.imlabel), [data[0]])
                            circ.append(IGate(label=self.imlabel), [data[1]])
                        else:
                            circ.i(data[1])
                else:
                    # left triangle bulk circuit
                    # (data clockwise from top of triangle, so 0, 2, 1)
                    da = self.code.z_gauges[j][0]
                    db = self.code.z_gauges[j][2]
                    dc = self.code.z_gauges[j][1]
                    # down X ancilla
                    a0 = self.code.n + self.z_ancilla_indices[j][0]
                    # up X ancilla
                    a1 = self.code.n + self.z_ancilla_indices[j][1]
                    # timestep - (idle accounted elsewhere)
                    if num_initialize is None:
                        circ.reset(a0)
                    # timestep 1
                    circ.cx(da, a0)
                    if add_idles:
                        circ.i(db)
                        circ.i(dc)
                    if num_initialize is None:
                        circ.reset(a1)
                    # timestep 2
                    if add_idles:
                        circ.i(da)
                    circ.cx(a0, db)
                    circ.cx(dc, a1)
                    # timestep 3
                    if add_idles:
                        circ.i(da)
                        circ.i(a0)
                        circ.i(dc)
                    circ.cx(db, a1)
                    # timestep 4
                    if add_idles:
                        circ.i(da)
                        if not distinct_measurement_idle:
                            circ.i(dc)
                    circ.cx(a0, db)
                    circ.measure(a1, creg[cbits[j]])
                    if not finalRound and num_initialize is not None:
                        circ.reset(a1)
                    # timestep 5
                    circ.cx(da, a0)
                    if add_idles:
                        if distinct_measurement_idle:
                            circ.append(IGate(label=self.imlabel), [da])
                            circ.append(IGate(label=self.imlabel), [db])
                            circ.append(IGate(label=self.imlabel), [dc])
                        else:
                            circ.i(db)
                            circ.i(dc)
        if add_barriers:
            circ.barrier(self.qreg)
        # Measure 'right' triangles second
        for j in range(len(self.code.z_gauges)):
            if self.code.z_orientations[j] == 1:
                if len(self.z_ancilla_indices[j]) == 3:
                    # right triangle boundary circuit
                    # (same as left triangle, currently)
                    data = copy.copy(self.code.z_gauges[j])
                    data.remove(-1)
                    a0 = self.code.n + self.z_ancilla_indices[j][0]
                    a1 = self.code.n + self.z_ancilla_indices[j][1]
                    a2 = self.code.n + self.z_ancilla_indices[j][2]
                    # timestep 1
                    if add_idles and not distinct_measurement_idle:
                        circ.i(data[0])
                        circ.i(data[1])
                    if num_initialize is None:
                        circ.reset(a2)
                    # timestep 2
                    if add_idles:
                        circ.i(data[0])
                    if num_initialize is None:
                        circ.reset(a0)
                        circ.reset(a1)
                    circ.cx(data[1], a2)
                    # timestep 3
                    circ.cx(data[0], a0)
                    circ.cx(a2, a1)
                    if add_idles:
                        circ.i(data[1])
                    # timestep 4
                    if add_idles:
                        circ.i(data[0])
                    circ.cx(a0, a1)
                    circ.cx(data[1], a2)
                    # timestep 5
                    circ.cx(data[0], a0)
                    circ.measure(a1, creg[cbits[j]])
                    if not finalRound and num_initialize is not None:
                        circ.reset(a1)
                    if add_idles:
                        if distinct_measurement_idle:
                            circ.append(IGate(label=self.imlabel), [data[0]])
                            circ.append(IGate(label=self.imlabel), [data[1]])
                        else:
                            circ.i(data[1])
                else:
                    # right triangle bulk circuit
                    da = self.code.z_gauges[j][0]
                    db = self.code.z_gauges[j][1]
                    dc = self.code.z_gauges[j][2]
                    # down X ancilla
                    a0 = self.code.n + self.z_ancilla_indices[j][0]
                    # up X ancilla
                    a1 = self.code.n + self.z_ancilla_indices[j][1]
                    # timestep - (idle accounted elsewhere)
                    if num_initialize is None:
                        circ.reset(a1)
                    # timestep 1
                    if add_idles:
                        circ.i(da)
                        circ.i(db)
                    if num_initialize is None:
                        circ.reset(a0)
                    circ.cx(dc, a1)
                    # timestep 2
                    circ.cx(da, a0)
                    circ.cx(a1, db)
                    if add_idles:
                        circ.i(dc)
                    # timestep 3
                    if add_idles:
                        circ.i(da)
                        circ.i(a1)
                        circ.i(dc)
                    circ.cx(db, a0)
                    # timestep 4
                    if add_idles:
                        if not distinct_measurement_idle:
                            circ.i(da)
                        circ.i(dc)
                    circ.measure(a0, creg[cbits[j]])
                    if not finalRound and num_initialize is not None:
                        circ.reset(a0)
                    circ.cx(a1, db)
                    # timestep 5
                    circ.cx(dc, a1)
                    if add_idles:
                        if distinct_measurement_idle:
                            circ.append(IGate(label=self.imlabel), [da])
                            circ.append(IGate(label=self.imlabel), [db])
                            circ.append(IGate(label=self.imlabel), [dc])
                        else:
                            circ.i(da)
                            circ.i(db)
        if add_barriers:
            circ.barrier(self.qreg)

    def syndrome_measurement(
        self, rounds=None, round_schedule=None, basis=None, initial_state=None, logical_paulis=None
    ):
        """Construct repeated syndrome measurement circuit.

        Method parameters will override the configuration that
        was provided to the constructor.

        rounds = number of repetitions of round_schedule (int)
        round_schedule = schedule of x/z gauge rounds (str)
        basis = initialization and measurement basis, x or z (str)
        initial_state = eigenvalue + or - of basis (str)
        logical_paulis = what logical Pauli to apply after each
          X or Z gauge round in round_schedule (str)

        Additional options from self.config:

        self.schedule: (str) If this equals
        "higgott-breuckmann", we will the schedule in
        Phys. Rev. X 11, 031039 (2021). If this equals
        "heavy-hex", we will use our own circuits.

        self.barriers: (bool) If this is True,
        insert barrier commands between steps of the syndrome circuit.

        self.idles: (bool) If this is True,
        insert identity gates at time steps where qubits are idle.
        Use a timing model where two-qubit gates and measurements
        have the same duration and single-qubit gates have negligible
        duration.

        self.distinct_measurement_idle: (bool) If
        this is True, insert an 'idm' labeled identity gate on the data
        qubits while the ancillas are measured, so the idle errors during
        measurement can be changed independently from idle errors elsewhere
        in the circuit.
        """
        # Get configuration from self.config
        add_barriers = self.barriers
        distinct_measurement_idle = self.distinct_measurement_idle
        rounds = self.rounds if rounds is None else rounds
        round_schedule = self.round_schedule if round_schedule is None else round_schedule
        basis = self.basis if basis is None else basis
        initial_state = self.initial_state if initial_state is None else initial_state
        logical_paulis = self.logical_paulis if logical_paulis is None else logical_paulis
        num_initialize = self.num_initialize
        # Compute the total number of classical bits
        xg = len(self.code.x_gauges)
        zg = len(self.code.z_gauges)
        total_cbits = 0
        for i in round_schedule:
            if i == "x":
                total_cbits += xg
            elif i == "z":
                total_cbits += zg
            else:
                raise Exception("round_schedule should contain x or z only")
        if len(logical_paulis) > 0 and len(logical_paulis) != len(round_schedule):
            raise Exception("len(logical_paulis) != len(round_schedule)")
        if self.schedule == "higgott-breuckmann":
            if len(logical_paulis) > 0 or distinct_measurement_idle:
                raise Exception("incompatible with higgott-breuckmann")
        total_cbits *= rounds
        total_cbits += self.code.n
        creg = ClassicalRegister(total_cbits)
        circ = QuantumCircuit(self.qreg, creg)
        # Initialize the data qubits
        if num_initialize is None:
            for i in range(self.code.n):
                circ.reset(self.qreg[i])
        else:
            for _ in range(num_initialize):
                circ.reset(self.qreg)
                circ.barrier(self.qreg)
        # Apply bit-wise Hadamard gate if X basis state preparation
        if basis == "x":
            for i in range(self.code.n):
                circ.h(self.qreg[i])
        # Apply a logical X to the first (only) logical qubit
        # if Z basis state preparation, otherwise apply
        # a logical Z if X basis state preparation
        if initial_state == "-":
            if basis == "z":
                for i in self.code.logical_x[0]:
                    circ.x(self.qreg[i])
            elif basis == "x":
                for i in self.code.logical_z[0]:
                    circ.z(self.qreg[i])
        if add_barriers:
            circ.barrier(self.qreg)
        # Construct the syndrome measurement circuit
        if self.schedule == "heavy-hex":
            start = 0
            finalRound = False
            for i in range(rounds):
                for j in range(len(round_schedule)):
                    if (i == rounds - 1) and (j == len(round_schedule) - 1):
                        finalRound = True
                    if round_schedule[j] == "x":
                        self._x_gauge_one_round_hex(
                            circ, creg, list(range(start, start + xg)), finalRound=finalRound
                        )
                        start += xg
                    if round_schedule[j] == "z":
                        self._z_gauge_one_round_hex(
                            circ, creg, list(range(start, start + zg)), finalRound=finalRound
                        )
                        start += zg
                    # Apply logical Pauli operator, if any
                    # (omit idle locations)
                    if len(logical_paulis) > 0:
                        if logical_paulis[j] == "x" or logical_paulis[j] == "y":
                            for k in self.code.logical_x[0]:
                                circ.x(self.qreg[k])
                        if logical_paulis[j] == "z" or logical_paulis[j] == "y":
                            for k in self.code.logical_z[0]:
                                circ.z(self.qreg[k])
                        if add_barriers:
                            circ.barrier(self.qreg)
        elif self.schedule == "higgott-breuckmann":
            start = 0
            full_schedule = round_schedule * rounds
            pairs = [full_schedule[i : i + 2] for i in range(0, len(full_schedule), 2)]
            leaving = ""
            for j in range(len(pairs) + 1):
                if j == len(pairs):
                    entering = ""
                else:
                    entering = pairs[j]
                num_bits_map = {"xx": 2 * xg, "zz": 2 * zg, "zx": zg + xg, "": 0}
                if leaving not in num_bits_map:
                    raise Exception("unknown round schedule %s" % leaving)
                numbits = num_bits_map[leaving]
                self._hb_round(leaving, entering, circ, creg, list(range(start, start + numbits)))
                start += numbits
                if j != len(pairs):
                    leaving = pairs[j]
        # Measure the data qubits
        for i in range(self.code.n):
            if basis == "x":
                circ.h(self.qreg[i])
            circ.measure(self.qreg[i], creg[start])
            start += 1

        return circ
