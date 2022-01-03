import unittest

from qiskit import execute, Aer

from plerco.pec_python.qec_code.hhc import HHC
from plerco.pec_python.qec_circuit.hhc_circuit import HHCCircuit
from plerco.pec_python.qec_decoder.hhc_decoder import HHCDecoder

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


class TestHHCDecoder(unittest.TestCase):
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

    def test_d3_1(self):  # 3 x z pymatching
        # Case that threw pymatching's contiguous qubit_id exception
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "x"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_2(self):  # 3 zx z pymatching
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_3(self):  # 3 zx z networkx
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_networkx"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_4(self):  # 3 zx z pymatching logical=xx
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["logical_paulis"] = "xx"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_5(self):  # 1 zxzxzx z pymatching logical=xyzxyz
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 1
        config["circuit"]["round_schedule"] = "zxzxzx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["logical_paulis"] = "xyzxyz"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_6(self):  # 3 zx z pymatching ted's schedule
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex-yoder"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_7(self):  # 3 zx z pymatching -1 eigenstate
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["initial_state"] = "-"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model, 1)
        # The stabilizer method is pretty slow (2-3 minutes)
        # self.correct_all_1(c, circ, dec, model, 1, False, "stabilizer")
        # The propagator method does not treat Paulis in circ as
        # errors, so the observed outcomes are not flipped (0 in arguments)
        self.correct_all_1(c, circ, dec, model, 0, False, "propagator")

    def test_d3_8(self):  # 3 zx z pymatching uniform
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["decoder"]["uniform"] = True
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_9(self):  # 3 zx z networkx uniform
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_networkx"
        config["decoder"]["uniform"] = True
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d3_10(self):  # 3 zx x pymatching
        # Some of these test fails if you make p = 0.01 since the edge
        # weights will lead to some 1st order errors not being
        # correctable. This is one of them.
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "x"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model, 0, True)
        self.correct_all_1(c, circ, dec, model, 0, True)

    def test_d3_11(self):  # 3 zx x networkx
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "x"
        config["decoder"]["method"] = "matching_networkx"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model, 0, True)
        self.correct_all_1(c, circ, dec, model, 0, True)

    def test_d3_12(self):  # 3 zx x pymatching -1 eigenstate
        c = HHC(3)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 3
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "x"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["initial_state"] = "-"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.0001)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model, 1, True)
        # The stabilizer method is pretty slow (2-3 minutes)
        # self.correct_all_1(c, circ, dec, model, 1, True, "stabilizer")
        # The propagator method does not treat Paulis in circ as
        # errors, so the observed outcomes are not flipped (0 in arguments)
        self.correct_all_1(c, circ, dec, model, 0, True, "propagator")

    def test_d5_1(self):  # 5 zx z pymatching
        c = HHC(5)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 5
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "z"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.005)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model)
        self.correct_all_1(c, circ, dec, model)

    def test_d5_2(self):  # 5 zx x pymatching
        c = HHC(5)
        config = Config()
        config["circuit"]["schedule"] = "heavy-hex"
        config["circuit"]["rounds"] = 5
        config["circuit"]["round_schedule"] = "zx"
        config["circuit"]["basis"] = "x"
        config["decoder"]["method"] = "matching_pymatching"
        config["circuit"]["distinct_measurement_idle"] = True
        config["circuit"]["barriers"] = False
        config["circuit"]["idles"] = True
        config["circuit"]["num_initialize"] = 1
        config["circuit"]["group_meas"] = True
        gen = HHCCircuit(c, config)
        circ = gen.syndrome_measurement()
        model = make_model(0.01)
        dec = HHCDecoder(c, circ, model, config)
        self.no_faults_success(c, circ, dec, model, 0, True)
        self.correct_all_1(c, circ, dec, model, 0, True)


if __name__ == "__main__":
    unittest.main()
