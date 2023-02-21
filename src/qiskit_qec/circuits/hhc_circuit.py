"""Object to construct quantum circuits for the HHC."""

from typing import List, Tuple
import logging

from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.circuit.library import IGate
from qiskit.circuit import Gate
from qiskit_qec.codes.hhc import HHC


def _hex_ancillas(code: HHC) -> Tuple[int, List[List[int]], List[List[int]]]:
    """Assign ancilla indices for X and Z gauge measurements.
    code = HHC object
    Return:
        total_ancillas is number of total ancillas.
        x_ancilla_indices list whose jth entry is the list of
        indices for the ancillas of the jth x gauge operator.
        z_ancilla_indices list defined likewise.
        The ancilla indices are absolute indices into self.qreg.
    """
    total_ancilla = 0
    # Append ancillas for X gauge operators
    x_ancilla_indices = []
    for j in range(len(code.x_gauges)):
        x_ancilla_indices.append([code.n + j])
    total_ancilla += len(code.x_gauges)
    # Append additional ancillas for Z gauge operators
    z_ancilla_indices = []
    for j, zg in enumerate(code.z_gauges):
        if len(zg) == 2:  # boundary Z gauge
            # Hex layout has 3 extra ancillas per boundary Z gauge
            start = total_ancilla
            z_ancilla_indices.append([code.n + start, code.n + start + 1, code.n + start + 2])
            total_ancilla += 3
        else:  # bulk Z gauge
            # These borrow ancillas from adjacent X gauges
            # and use one additional ancilla in the center of the face
            for k, xg in enumerate(code.x_gauges):
                # qubits at indices 0 and 2 are shared
                # with adjacent X gauge on left
                if zg[0] in xg and zg[2] in xg:
                    left = code.n + k
                # qubits at indices 1 and 3 are shared
                # with adjacent X gauge on right
                if zg[1] in xg and zg[3] in xg:
                    right = code.n + k
            z_ancilla_indices.append([left, code.n + total_ancilla, right])
            total_ancilla += 1
    logging.info("total_ancilla = %s", total_ancilla)
    logging.info("x_ancilla_indices = %s", x_ancilla_indices)
    logging.info("z_ancilla_indices = %s", z_ancilla_indices)
    return total_ancilla, x_ancilla_indices, z_ancilla_indices


class HHCCircuit:
    """Create quantum circuits for syndrome measurements.

    Specialized to heavy hexagon subsystem code.
    """

    def __init__(
        self,
        hhc: HHC,
        barriers: bool,
        idles: bool,
        distinct_measurement_idle: bool,
        init_error: bool,
        group_meas: bool,
        xprs: bool,
        blocks: int,
        round_schedule: str,
        basis: str,
        initial_state: str,
        logical_paulis: str,
        num_initialize: int,
        idle_before_measure: bool,
    ):
        """Create an object associated to a HHC."""
        self.code = hhc
        self._validate(
            barriers,
            idles,
            distinct_measurement_idle,
            init_error,
            group_meas,
            xprs,
            blocks,
            round_schedule,
            basis,
            initial_state,
            logical_paulis,
            num_initialize,
            idle_before_measure,
        )
        self.total_ancilla, self.x_ancilla_indices, self.z_ancilla_indices = _hex_ancillas(hhc)
        self.total_qubits = hhc.n + self.total_ancilla
        self.qreg = QuantumRegister(self.total_qubits)
        self.imlabel = "idm"  # label for idle during measurement
        self.ibmlabel = "idbm"  # label for idle before measurement
        self.iinitlabel = (
            "idinit"  # label for idle in the beginning of circuit to add initialization error
        )

    def _validate(
        self,
        barriers: bool,
        idles: bool,
        distinct_measurement_idle: bool,
        init_error: bool,
        group_meas: bool,
        xprs: bool,
        blocks: int,
        round_schedule: str,
        basis: str,
        initial_state: str,
        logical_paulis: str,
        num_initialize: int,
        idle_before_measurement: bool,
    ):
        """Validate parameters."""
        self.barriers = barriers
        self.idles = idles
        self.distinct_measurement_idle = distinct_measurement_idle
        self.init_error = init_error
        self.group_meas = group_meas
        self.xprs = xprs
        self.blocks = blocks
        if self.blocks < 1:
            raise Exception("expected positive integer rounds")
        self.round_schedule = round_schedule
        if set(self.round_schedule) > set("xz"):
            raise Exception("expected round schedule of 'x', 'z' chars")
        self.basis = basis
        if not (self.basis in ["x", "z"]):
            raise Exception("expected basis to be 'x' or 'z'")
        self.initial_state = initial_state
        if not (self.initial_state in ["+", "-"]):
            raise Exception("expected initial state '+' or '-'")
        self.logical_paulis = logical_paulis
        if set(self.logical_paulis) > set("ixyz"):
            raise Exception("expected 'i', 'x', 'y', 'z' logical paulis")
        if len(self.logical_paulis) > 0 and len(self.logical_paulis) != len(round_schedule):
            raise Exception("len(logical paulis) != len(round schedule)")
        self.num_initialize = num_initialize
        if self.num_initialize is not None:
            if self.num_initialize < 0:
                raise Exception("expected zero or positive integer number of initialization")
        self.idle_before_measurement = idle_before_measurement

    def _x_gauge_one_round_hex(
        self, circ, creg, cbits, basis=None, finalRound=False, logical_pauli=None
    ):
        """Measure all of the X gauge operators.

        Uses a single ancilla per X gauge operator.
        Append the gates to the input circuit circ. Measure results into
        the bits at the indices in cbits (list) in creg.
        """
        # Get options from configuration
        add_barriers = self.barriers
        add_idles = self.idles
        distinct_measurement_idle = self.distinct_measurement_idle
        num_initialize = self.num_initialize
        group_meas = self.group_meas
        with_xprs = self.xprs
        basis = self.basis if basis is None else basis
        ibm = self.idle_before_measurement
        # Add xprs to use as reset right after measurement
        xprs = Gate(name="xprs", num_qubits=1, params=[])

        # Initialize the ancilla (assume occured in previous step)
        # i.e., no idle locations here
        for i in range(len(self.code.x_gauges)):
            for j in self.x_ancilla_indices[i]:
                if num_initialize is None:
                    circ.reset(self.qreg[j])
                circ.h(self.qreg[j])
        if add_barriers and num_initialize is None:
            circ.barrier(self.qreg)
        # Interact with the ancilla
        for step in range(2):
            idle = list(range(self.code.n))
            for i, xg in enumerate(self.code.x_gauges):
                circ.cx(self.qreg[self.x_ancilla_indices[i][0]], self.qreg[xg[step]])
                idle.remove(xg[step])
            if add_idles:
                for i in idle:
                    circ.i(self.qreg[i])
            if add_barriers:
                circ.barrier(self.qreg)
        # Measure the ancilla
        for i in range(len(self.code.x_gauges)):
            for j in self.x_ancilla_indices[i]:
                circ.h(self.qreg[j])
                if not group_meas:
                    if ibm:
                        circ.append(IGate(label=self.ibmlabel + str(cbits[i])), [self.qreg[j]])
                    circ.measure(self.qreg[j], creg[cbits[i]])
                    if not finalRound and num_initialize is not None:
                        if with_xprs:
                            circ.append(xprs, [self.qreg[j]])
                        else:
                            circ.reset(self.qreg[j])
        if group_meas:
            if logical_pauli in ["x", "y"]:
                for k in self.code.logical_x[0]:
                    circ.x(self.qreg[k])
            if logical_pauli in ["z", "y"]:
                for k in self.code.logical_z[0]:
                    circ.z(self.qreg[k])
            if add_barriers and logical_pauli is not None:
                circ.barrier(self.qreg)
            if basis == "x" and finalRound:
                for i in range(self.code.n):
                    circ.h(self.qreg[i])
            if not add_barriers:
                circ.barrier(self.qreg)
            for i in range(len(self.code.x_gauges)):
                for j in self.x_ancilla_indices[i]:
                    if ibm:
                        circ.append(IGate(label=self.ibmlabel + str(cbits[i])), [self.qreg[j]])
                    circ.measure(self.qreg[j], creg[cbits[i]])
                    if not finalRound and num_initialize is not None:
                        if with_xprs:
                            circ.append(xprs, [self.qreg[j]])
                        else:
                            circ.reset(self.qreg[j])
        # Data is idle during measurement/reset
        if add_idles:
            for i in range(self.code.n):
                if distinct_measurement_idle:
                    if group_meas:
                        if not finalRound:
                            circ.append(IGate(label=self.imlabel), [self.qreg[i]])
                    else:
                        circ.append(IGate(label=self.imlabel), [self.qreg[i]])
                else:
                    circ.i(self.qreg[i])
        if add_barriers or (group_meas and not finalRound):
            circ.barrier(self.qreg)

    def _z_gauge_one_round_hex(
        self,
        circ: QuantumCircuit,
        creg,
        cbits,
        lflags,
        rflags,
        basis=None,
        finalRound=False,
        logical_pauli=None,
    ):
        """Measure all of the Z gauge operators.

        Reuses the ancillas for X gauge operator and one additional per gauge.
        Append the gates to the input circuit circ. Measure results into
        the bits at the indices in cbits (list) in creg. Measure the left flags
        into the bits at indices in lflags (list) and right flags into rflags.
        """
        # Get options from configuration
        add_barriers = self.barriers
        add_idles = self.idles
        distinct_measurement_idle = self.distinct_measurement_idle
        num_initialize = self.num_initialize
        group_meas = self.group_meas
        with_xprs = self.xprs
        basis = self.basis if basis is None else basis
        ibm = self.idle_before_measurement
        # Add xprs to use as reset right after measurement
        xprs = Gate(name="xprs", num_qubits=1, params=[])
        for step in range(7):
            idle = list(range(self.code.n))
            for j, zg in enumerate(self.code.z_gauges):
                # Label the ancillas
                left = self.qreg[self.z_ancilla_indices[j][0]]
                middle = self.qreg[self.z_ancilla_indices[j][1]]
                right = self.qreg[self.z_ancilla_indices[j][2]]
                # Identify the type of gauge operator
                is_boundary = len(zg) == 2
                is_top_row = zg[0] < self.code.d
                # Follow the schedule in Fig 2b of arXiv:1907.09528,
                # except conjugate the circuit by Hadamard.
                # (i.e. switch +/0, X/Z, and CNOT direction in Fig 4)
                if step == 0:  # done in previous steps, no idles
                    # (although we added them anyway)
                    if num_initialize is None:
                        circ.reset(middle)
                        circ.reset(right)
                    circ.h(right)
                elif step == 1:
                    circ.cx(right, middle)
                    if num_initialize is None:
                        circ.reset(left)
                    circ.h(left)
                elif step == 2:
                    right_idle = True
                    circ.cx(left, middle)
                    if is_boundary and not is_top_row:
                        circ.cx(self.qreg[zg[1]], right)
                        idle.remove(zg[1])
                        right_idle = False
                    elif not is_boundary:
                        circ.cx(self.qreg[zg[1]], right)
                        idle.remove(zg[1])
                        right_idle = False
                    if right_idle and add_idles:
                        circ.i(right)
                elif step == 3:
                    ancilla_idle = True
                    if is_boundary and is_top_row:
                        circ.cx(zg[0], left)
                        circ.cx(zg[1], right)
                        idle.remove(zg[0])
                        idle.remove(zg[1])
                        ancilla_idle = False
                    elif not is_boundary:
                        circ.cx(zg[2], left)
                        circ.cx(zg[3], right)
                        idle.remove(zg[2])
                        idle.remove(zg[3])
                        ancilla_idle = False
                    if ancilla_idle and add_idles:
                        circ.i(left)
                        circ.i(right)
                elif step == 4:
                    left_idle = True
                    circ.cx(right, middle)
                    if is_boundary and not is_top_row:
                        circ.cx(zg[0], left)
                        idle.remove(zg[0])
                        left_idle = False
                    elif not is_boundary:
                        circ.cx(zg[0], left)
                        idle.remove(zg[0])
                        left_idle = False
                    if left_idle and add_idles:
                        circ.i(left)
                elif step == 5:
                    circ.cx(left, middle)
                    if group_meas:
                        if add_idles:
                            circ.i(right)
                    else:
                        circ.h(right)
                        if ibm:
                            circ.append(IGate(label=self.ibmlabel + str(rflags[j])), [right])
                        circ.measure(right, creg[rflags[j]])
                        if not finalRound and num_initialize is not None:
                            if with_xprs:
                                circ.append(xprs, [right])
                            else:
                                circ.reset(right)
                elif step == 6:
                    circ.h(left)
                    if group_meas:
                        circ.h(right)
                    else:
                        if ibm:
                            circ.append(IGate(label=self.ibmlabel + str(cbits[j])), [middle])
                            circ.append(IGate(label=self.ibmlabel + str(lflags[j])), [left])
                        circ.measure(middle, creg[cbits[j]])
                        circ.measure(left, creg[lflags[j]])
                        if not finalRound and num_initialize is not None:
                            if with_xprs:
                                circ.append(xprs, [middle])
                                circ.append(xprs, [left])
                            else:
                                circ.reset(middle)
                                circ.reset(left)

            if step == 6 and group_meas:
                if logical_pauli in ["x", "y"]:
                    for k in self.code.logical_x[0]:
                        circ.x(self.qreg[k])
                if logical_pauli in ["z", "y"]:
                    for k in self.code.logical_z[0]:
                        circ.z(self.qreg[k])
                if add_barriers and logical_pauli is not None:
                    circ.barrier(self.qreg)
                if basis == "x" and finalRound:
                    for i in range(self.code.n):
                        circ.h(self.qreg[i])
                if not add_barriers:
                    circ.barrier(self.qreg)
                for j in range(len(self.code.z_gauges)):
                    # Label the ancillas
                    left = self.qreg[self.z_ancilla_indices[j][0]]
                    middle = self.qreg[self.z_ancilla_indices[j][1]]
                    right = self.qreg[self.z_ancilla_indices[j][2]]
                    if ibm:
                        circ.append(IGate(label=self.ibmlabel + str(rflags[j])), [right])
                        circ.append(IGate(label=self.ibmlabel + str(cbits[j])), [middle])
                        circ.append(IGate(label=self.ibmlabel + str(lflags[j])), [left])
                    circ.measure(right, creg[rflags[j]])
                    circ.measure(middle, creg[cbits[j]])
                    circ.measure(left, creg[lflags[j]])
                    if not finalRound and num_initialize is not None:
                        if with_xprs:
                            circ.append(xprs, [middle])
                            circ.append(xprs, [left])
                            circ.append(xprs, [right])
                        else:
                            circ.reset(middle)
                            circ.reset(left)
                            circ.reset(right)
            if add_idles:
                for i in idle:
                    if distinct_measurement_idle and step == 6:
                        if group_meas:
                            if not finalRound:
                                circ.append(IGate(label=self.imlabel), [self.qreg[i]])
                        else:
                            circ.append(IGate(label=self.imlabel), [self.qreg[i]])
                    elif step == 0:
                        if num_initialize is None:
                            circ.i(self.qreg[i])
                    else:
                        circ.i(self.qreg[i])
            if add_barriers or (step == 6 and group_meas and not finalRound):
                if step != 0:
                    circ.barrier(self.qreg)
                elif num_initialize is None:
                    circ.barrier(self.qreg)

    def syndrome_measurement(
        self,
        rounds=None,
        round_schedule=None,
        basis=None,
        initial_state=None,
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

        Additional options from configuration:

        barriers: (bool) If this is True,
        insert barrier commands between steps of the syndrome circuit.

        idles: (bool) If this is True,
        insert identity gates at time steps where qubits are idle.
        Use a timing model where two-qubit gates and measurements
        have the same duration and single-qubit gates have negligible
        duration.

        distinct_measurement_idle: (bool) If
        this is True, insert an 'idm' labeled identity gate on the data
        qubits while the ancillas are measured, so the idle errors during
        measurement can be changed independently from idle errors elsewhere
        in the circuit.
        """
        # Get configuration
        add_barriers = self.barriers
        rounds = self.blocks if rounds is None else rounds
        round_schedule = self.round_schedule if round_schedule is None else round_schedule
        basis = self.basis if basis is None else basis
        initial_state = self.initial_state if initial_state is None else initial_state
        logical_paulis = self.logical_paulis
        num_initialize = self.num_initialize
        init_error = self.init_error
        group_meas = self.group_meas
        ibm = self.idle_before_measurement
        # Compute the total number of classical bits
        xg = len(self.code.x_gauges)
        zg = len(self.code.z_gauges)
        total_cbits = 0
        for i in round_schedule:
            if i == "x":
                total_cbits += xg
            elif i == "z":
                total_cbits += 3 * zg  # include bits for flag outcomes
            else:
                raise Exception("round_schedule should contain x or z only")
        if len(logical_paulis) > 0 and len(logical_paulis) != len(round_schedule):
            raise Exception("logical_paulis not same length as round_schedule")
        total_cbits *= rounds
        total_cbits += self.code.n
        creg = ClassicalRegister(total_cbits)
        circ = QuantumCircuit(self.qreg, creg)
        # Initialize the data qubits
        if num_initialize is None:
            for i in range(self.code.n):
                circ.reset(self.qreg[i])
        else:
            if num_initialize == 0 and init_error:
                circ.append(IGate(label=self.iinitlabel), [self.qreg])
                circ.barrier(self.qreg)
            for _ in range(num_initialize):
                circ.reset(self.qreg)
            if num_initialize > 0:
                circ.barrier(self.qreg)

        # Apply bit-wise Hadamard if X basis state preparation
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
            if num_initialize is None or basis == "x" or (basis == "z" and initial_state == "-"):
                circ.barrier(self.qreg)

        # Construct the syndrome measurement circuit
        start = 0
        final_round = False
        for i in range(rounds):
            for j, rs in enumerate(round_schedule):
                if (i == rounds - 1) and (j == len(round_schedule) - 1):
                    final_round = True
                log_paul = "" if len(logical_paulis) == 0 else logical_paulis[j]
                if rs == "x":
                    self._x_gauge_one_round_hex(
                        circ,
                        creg,
                        list(range(start, start + xg)),
                        basis=basis,
                        finalRound=final_round,
                        logical_pauli=log_paul,
                    )
                    start += xg
                if rs == "z":
                    # Outcomes go into first zg bits
                    # Left flags into second zg bits
                    # Right flags into third zg bits
                    self._z_gauge_one_round_hex(
                        circ,
                        creg,
                        list(range(start, start + zg)),
                        list(range(start + zg, start + 2 * zg)),
                        list(range(start + 2 * zg, start + 3 * zg)),
                        basis=basis,
                        finalRound=final_round,
                        logical_pauli=log_paul,
                    )
                    start += 3 * zg
                # Apply logical Pauli operator, if any
                # (omit idle locations)
                if len(logical_paulis) > 0 and not group_meas:
                    if logical_paulis[j] == "x" or logical_paulis[j] == "y":
                        for k in self.code.logical_x[0]:
                            circ.x(self.qreg[k])
                    if logical_paulis[j] == "z" or logical_paulis[j] == "y":
                        for k in self.code.logical_z[0]:
                            circ.z(self.qreg[k])
                    if add_barriers:
                        circ.barrier(self.qreg)
        # Measure the data qubits
        for i in range(self.code.n):
            if basis == "x" and not group_meas:
                circ.h(self.qreg[i])
            if ibm:
                circ.append(IGate(label=self.ibmlabel + str(start)), [self.qreg[i]])
            circ.measure(self.qreg[i], creg[start])
            start += 1

        return circ
