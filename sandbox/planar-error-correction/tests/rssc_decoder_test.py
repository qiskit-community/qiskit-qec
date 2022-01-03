import unittest

from qiskit import execute, Aer

from plerco.pec_python.qec_code.rssc import RSSC
from plerco.pec_python.qec_circuit.rssc_circuit import RSSCCircuit
from plerco.pec_python.qec_decoder.rssc_decoder import RSSCDecoder

from plerco.pec_python.qec_noise_model.paulinoisemodel import PauliNoiseModel
from plerco.pec_python.qec_decoder.decoder_utils.faultenumerator import FaultEnumerator
from plerco.pec_python.configuration.config import Config


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
    def correct_all_1(
        self, code, circ, dec, model, good=0, xbasis=False, method="propagator"
    ):
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

    def count_all_2(
        self, code, circ, dec, model, good=0, xbasis=False, method="propagator"
    ):
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
        config = Config()
        config["circuit"]["schedule"] = "higgott-breuckmann"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["distinct_measurement_idle"] = False
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_2(self):  # 3 zx z pymatching hex-layout
        c = RSSC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_3(self):  # 3 zx z pymatching hex-layout logical=xx
        c = RSSC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["logical_paulis"] = "xx"
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_4(self):  # 1 zxzxzx z pymatching hex-layout logical=xyzxyz
        c = RSSC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 1
        config["circuit"]["round_schedule"] = "zxzxzx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["logical_paulis"] = "xyzxyz"
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_5(self):  # 3 zx z pymatching -1 eigenvalue
        c = RSSC(3)
        config = Config()
        config["circuit"]["schedule"] = "higgott-breuckmann"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["initial_state"] = "-"
        config["circuit"]["distinct_measurement_idle"] = False
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model, good=1)
        self.correct_all_1(
            c, circ, dec, model, good=0, xbasis=False, method="propagator"
        )

    def test_d3_6(self):  # 3 zx z pymatching uniform
        c = RSSC(3)
        config = Config()
        config["circuit"]["schedule"] = "higgott-breuckmann"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["decoder"]["uniform"] = True
        config["circuit"]["distinct_measurement_idle"] = False
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_7(self):  # 3 zx z networkx
        c = RSSC(3)
        config = Config()
        config["circuit"]["schedule"] = "higgott-breuckmann"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_networkx"
        config["circuit"]["distinct_measurement_idle"] = False
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_8(self):  # 3 zx z networkx uniform
        c = RSSC(3)
        config = Config()
        config["circuit"]["schedule"] = "higgott-breuckmann"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_networkx"
        config["decoder"]["uniform"] = True
        config["circuit"]["distinct_measurement_idle"] = False
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_9(self):  # 3 zx x pymatching
        c = RSSC(3)
        config = Config()
        config["circuit"]["schedule"] = "higgott-breuckmann"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "x"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["distinct_measurement_idle"] = False
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model, good=0, xbasis=True)
        self.correct_all_1(
            c, circ, dec, model, good=0, xbasis=True, method="propagator"
        )

    def test_d3_10(self):  # 3 zx x pymatching -1 eigenvalue
        c = RSSC(3)
        config = Config()
        config["circuit"]["schedule"] = "higgott-breuckmann"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "x"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["initial_state"] = "-"
        config["circuit"]["distinct_measurement_idle"] = False
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model, good=1, xbasis=True)
        self.correct_all_1(
            c, circ, dec, model, good=0, xbasis=True, method="propagator"
        )

    def test_d3_11(self):  # 3 zx x networkx
        c = RSSC(3)
        config = Config()
        config["circuit"]["schedule"] = "higgott-breuckmann"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "x"
        config["decoder"]["method"] = "matching_networkx"
        config["circuit"]["distinct_measurement_idle"] = False
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model, good=0, xbasis=True)
        self.correct_all_1(
            c, circ, dec, model, good=0, xbasis=True, method="propagator"
        )

    def test_d5_1(self):  # 5 zx z pymatching
        c = RSSC(5)
        config = Config()
        config["circuit"]["schedule"] = "higgott-breuckmann"
        config["circuit"]["rounds"] = 5
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["distinct_measurement_idle"] = False
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d5_2(self):  # 3 zxzx x pymatching
        c = RSSC(5)
        config = Config()
        config["circuit"]["schedule"] = "higgott-breuckmann"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "x"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["distinct_measurement_idle"] = False
        gen = RSSCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = RSSCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model, good=0, xbasis=True)
        self.correct_all_1(
            c, circ, dec, model, good=0, xbasis=True, method="propagator"
        )


if __name__ == "__main__":
    unittest.main()
