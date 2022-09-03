import unittest

from qiskit import Aer, execute
from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.circuits.rotated_subsystem_surface_code_circuit import RSSCCircuit
from qiskit_qec.codes.rotated_surface_code import RSSC
from qiskit_qec.decoders.rssc_decoder import RSSCDecoder
from qiskit_qec.noise import PauliNoiseModel


def make_model(p):
    """Make a Pauli model for depolarizing noise."""
    pnm = PauliNoiseModel()
    pnm.add_operation(
        "cx",
        {
            "ix": 1,
            "iy": 1,
            "iz": 1,
            "xi": 1,
            "xx": 1,
            "xy": 1,
            "xz": 1,
            "yi": 1,
            "yx": 1,
            "yy": 1,
            "yz": 1,
            "zi": 1,
            "zx": 1,
            "zy": 1,
            "zz": 1,
        },
    )
    pnm.add_operation("id", {"x": 1, "y": 1, "z": 1})
    pnm.add_operation("reset", {"x": 1})
    pnm.add_operation("measure", {"x": 1})
    pnm.add_operation("h", {"x": 1, "y": 1, "z": 1})
    pnm.add_operation("x", {"x": 1, "y": 1, "z": 1})
    pnm.add_operation("y", {"x": 1, "y": 1, "z": 1})
    pnm.add_operation("z", {"x": 1, "y": 1, "z": 1})
    pnm.add_operation("idm", {"x": 1, "y": 1, "z": 1})
    pnm.set_error_probability("cx", p)
    pnm.set_error_probability("id", p)
    pnm.set_error_probability("reset", p)
    pnm.set_error_probability("measure", p)
    pnm.set_error_probability("h", p)
    pnm.set_error_probability("x", p)
    pnm.set_error_probability("y", p)
    pnm.set_error_probability("z", p)
    pnm.set_error_probability("idm", p)
    return pnm


class TestRSSCDecoder(unittest.TestCase):
    def correct_all_1(self, code, circ, dec, model, good=0, xbasis=False, method="propagator"):
        dec.update_edge_weights(model)
        fe = FaultEnumerator(circ, method=method, model=model)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            if xbasis:
                fail = code.logical_z_error(corrected_outcomes)
            else:
                fail = code.logical_x_error(corrected_outcomes)
            if fail[0] != good:
                failures += 1
        self.assertEqual(failures, 0)

    def count_all_2(self, code, circ, dec, model, good=0, xbasis=False, method="propagator"):
        dec.update_edge_weights(model)
        # This is too slow to use in our tests
        fe = FaultEnumerator(circ, order=2, method=method, model=model)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            if xbasis:
                fail = code.logical_z_error(corrected_outcomes)
            else:
                fail = code.logical_x_error(corrected_outcomes)
            if fail[0] != good:
                failures += 1
        return failures

    def no_faults_success(self, code, circ, dec, model, good=0, xbasis=False):
        shots = 10
        seed = 100
        result = execute(
            circ,
            Aer.get_backend("qasm_simulator"),
            method="stabilizer",
            shots=shots,
            optimization_level=0,
            seed_simulator=seed,
        ).result()
        counts = result.get_counts(circ)
        dec.update_edge_weights(model)
        failures = 0
        for (outcome, num) in counts.items():
            reversed_outcome = list(map(int, outcome[::-1]))
            corrected_outcomes = dec.process(reversed_outcome)
            if xbasis:
                fail = code.logical_z_error(corrected_outcomes)
            else:
                fail = code.logical_x_error(corrected_outcomes)
            if fail[0] != good:
                failures += 1
        self.assertEqual(failures, 0)

    def test_d3_1(self):  # 3 zx z pymatching
        c = RSSC(3)

        schedule = "higgott-breuckmann"
        rounds = 3
        round_schedule = "zx"
        basis = "z"
        method = "matching_pymatching"
        distinct_measurement_idle = False
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            distinct_measurement_idle=distinct_measurement_idle,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(
            c,
            circ,
            model,
            method=method,
        )
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_2(self):  # 3 zx z pymatching hex-layout
        c = RSSC(3)

        schedule = "heavy-hex"
        rounds = 3
        round_schedule = "zx"
        basis = "z"
        method = "matching_pymatching"
        gen = RSSCCircuit(
            rssc=c, schedule=schedule, rounds=rounds, round_schedule=round_schedule, basis=basis
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_3(self):  # 3 zx z pymatching hex-layout logical=xx
        c = RSSC(3)

        schedule = "heavy-hex"
        rounds = 3
        round_schedule = "zx"
        basis = "z"
        method = "matching_pymatching"
        logical_paulis = "xx"
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            logical_paulis=logical_paulis,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_4(self):  # 1 zxzxzx z pymatching hex-layout logical=xyzxyz
        c = RSSC(3)

        schedule = "heavy-hex"
        rounds = 1
        round_schedule = "zxzxzx"
        basis = "z"
        method = "matching_pymatching"
        logical_paulis = "xyzxyz"
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            logical_paulis=logical_paulis,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_5(self):  # 3 zx z pymatching -1 eigenvalue
        c = RSSC(3)

        schedule = "higgott-breuckmann"
        rounds = 3
        round_schedule = "zx"
        basis = "z"
        method = "matching_pymatching"
        initial_state = "-"
        distinct_measurement_idle = False
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            distinct_measurement_idle=distinct_measurement_idle,
            initial_state=initial_state,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model, good=1)
        self.correct_all_1(c, circ, dec, model, good=0, xbasis=False, method="propagator")

    def test_d3_6(self):  # 3 zx z pymatching uniform
        c = RSSC(3)

        schedule = "higgott-breuckmann"
        rounds = 3
        round_schedule = "zx"
        basis = "z"
        method = "matching_pymatching"
        uniform = True
        distinct_measurement_idle = False
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            distinct_measurement_idle=distinct_measurement_idle,
            uniform=uniform,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_7(self):  # 3 zx z networkx
        c = RSSC(3)

        schedule = "higgott-breuckmann"
        rounds = 3
        round_schedule = "zx"
        basis = "z"
        method = "matching_networkx"
        distinct_measurement_idle = False
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            distinct_measurement_idle=distinct_measurement_idle,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_8(self):  # 3 zx z networkx uniform
        c = RSSC(3)

        schedule = "higgott-breuckmann"
        rounds = 3
        round_schedule = "zx"
        basis = "z"
        method = "matching_networkx"
        uniform = True
        distinct_measurement_idle = False
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            distinct_measurement_idle=distinct_measurement_idle,
            uniform=uniform,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_9(self):  # 3 zx x pymatching
        c = RSSC(3)

        schedule = "higgott-breuckmann"
        rounds = 3
        round_schedule = "zx"
        basis = "x"
        method = "matching_pymatching"
        distinct_measurement_idle = False
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            distinct_measurement_idle=distinct_measurement_idle,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model, good=0, xbasis=True)
        self.correct_all_1(c, circ, dec, model, good=0, xbasis=True, method="propagator")

    def test_d3_10(self):  # 3 zx x pymatching -1 eigenvalue
        c = RSSC(3)

        schedule = "higgott-breuckmann"
        rounds = 3
        round_schedule = "zx"
        basis = "x"
        method = "matching_pymatching"
        initial_state = "-"
        distinct_measurement_idle = False
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            distinct_measurement_idle=distinct_measurement_idle,
            initial_state=initial_state,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model, good=1, xbasis=True)
        self.correct_all_1(c, circ, dec, model, good=0, xbasis=True, method="propagator")

    def test_d3_11(self):  # 3 zx x networkx
        c = RSSC(3)

        schedule = "higgott-breuckmann"
        rounds = 3
        round_schedule = "zx"
        basis = "x"
        method = "matching_networkx"
        distinct_measurement_idle = False
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            distinct_measurement_idle=distinct_measurement_idle,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model, good=0, xbasis=True)
        self.correct_all_1(c, circ, dec, model, good=0, xbasis=True, method="propagator")

    def test_d5_1(self):  # 5 zx z pymatching
        c = RSSC(5)

        schedule = "higgott-breuckmann"
        rounds = 5
        round_schedule = "zx"
        basis = "z"
        method = "matching_pymatching"
        distinct_measurement_idle = False
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            distinct_measurement_idle=distinct_measurement_idle,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d5_2(self):  # 3 zxzx x pymatching
        c = RSSC(5)

        schedule = "higgott-breuckmann"
        rounds = 3
        round_schedule = "zx"
        basis = "x"
        method = "matching_pymatching"
        distinct_measurement_idle = False
        gen = RSSCCircuit(
            rssc=c,
            schedule=schedule,
            rounds=rounds,
            round_schedule=round_schedule,
            basis=basis,
            distinct_measurement_idle=distinct_measurement_idle,
        )
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, method=method)
        self.no_faults_success(c, circ, dec, model, good=0, xbasis=True)
        self.correct_all_1(c, circ, dec, model, good=0, xbasis=True, method="propagator")


if __name__ == "__main__":
    unittest.main()
