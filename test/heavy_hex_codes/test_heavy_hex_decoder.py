"""Test a heavy-hexagon code decoder"""
import unittest

from qiskit_aer import Aer

from qiskit_qec.codes.hhc import HHC
from qiskit_qec.circuits.hhc_circuit import HHCCircuit
from qiskit_qec.decoders.hhc_decoder import HHCDecoder
from qiskit_qec.decoders.temp_code_util import temp_syndrome
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel
from qiskit_qec.analysis.faultenumerator import FaultEnumerator


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
    """Tests for a heavy-hexagon code decoder."""

    def setUp(self) -> None:
        """Work we can do once."""
        self.model = make_model(0.0001)

    def correct_all_1(self, code, circ, dec, model, good=0, xbasis=False, method="propagator"):
        """Test if we can correct all single-location faults."""
        dec.update_edge_weights(model)
        fe = FaultEnumerator(circ, method=method, model=model)
        failures = 0
        for faultpath in fe.generate():
            outcome = faultpath[3]
            corrected_outcomes = dec.process(outcome)
            if xbasis:
                fail = temp_syndrome(corrected_outcomes, [range(code.n)])
            else:
                fail = temp_syndrome(corrected_outcomes, [range(code.n)])
            if fail[0] != good:
                failures += 1
                print(good, fail, faultpath, corrected_outcomes)
        self.assertEqual(failures, 0)

    def no_faults_success(self, code, circ, dec, model, good=0, xbasis=False):
        """Test for correct behavior without faults."""
        shots = 10
        seed = 100
        backend = Aer.get_backend("aer_simulator")
        options = {"method": "stabilizer", "shots": shots, "seed_simulator": seed}
        result = backend.run(circ, **options).result()
        counts = result.get_counts(circ)
        dec.update_edge_weights(model)
        failures = 0
        for outcome, _ in counts.items():
            reversed_outcome = list(map(int, outcome[::-1]))
            corrected_outcomes = dec.process(reversed_outcome)
            if xbasis:
                fail = temp_syndrome(corrected_outcomes, [range(code.n)])
            else:
                fail = temp_syndrome(corrected_outcomes, [range(code.n)])
            if fail[0] != good:
                print(good, fail, reversed_outcome, corrected_outcomes)
                failures += 1
        self.assertEqual(failures, 0)

    def test_d3_2(self):
        """Check 3, zx, z, pymatching."""
        blocks = 3
        round_schedule = "zx"
        basis = "z"
        logical_paulis = "ii"
        c = HHC(3)
        gen = HHCCircuit(
            c,
            barriers=True,
            idles=True,
            distinct_measurement_idle=True,
            init_error=True,
            group_meas=False,
            xprs=False,
            blocks=blocks,
            round_schedule=round_schedule,
            basis=basis,
            initial_state="+",
            logical_paulis=logical_paulis,
            num_initialize=1,
            idle_before_measure=False,
        )
        circ = gen.syndrome_measurement()
        dec = HHCDecoder(
            n=c.n,
            css_x_gauge_ops=c.x_gauges,
            css_x_stabilizer_ops=c.x_stabilizers,
            css_x_boundary=c.x_boundary,
            css_z_gauge_ops=c.z_gauges,
            css_z_stabilizer_ops=c.z_stabilizers,
            css_z_boundary=c.z_boundary,
            circuit=circ,
            model=self.model,
            basis=basis,
            round_schedule=round_schedule,
            blocks=blocks,
            method="pymatching",
            uniform=False,
        )
        self.no_faults_success(c, circ, dec, self.model)
        self.correct_all_1(c, circ, dec, self.model)

    def test_d3_3(self):
        """Check 3, zx, z, rustworkx."""
        blocks = 3
        round_schedule = "zx"
        basis = "z"
        logical_paulis = "ii"
        c = HHC(3)
        gen = HHCCircuit(
            c,
            barriers=True,
            idles=True,
            distinct_measurement_idle=True,
            init_error=True,
            group_meas=False,
            xprs=False,
            blocks=blocks,
            round_schedule=round_schedule,
            basis=basis,
            initial_state="+",
            logical_paulis=logical_paulis,
            num_initialize=1,
            idle_before_measure=False,
        )
        circ = gen.syndrome_measurement()
        dec = HHCDecoder(
            n=c.n,
            css_x_gauge_ops=c.x_gauges,
            css_x_stabilizer_ops=c.x_stabilizers,
            css_x_boundary=c.x_boundary,
            css_z_gauge_ops=c.z_gauges,
            css_z_stabilizer_ops=c.z_stabilizers,
            css_z_boundary=c.z_boundary,
            circuit=circ,
            model=self.model,
            basis=basis,
            round_schedule=round_schedule,
            blocks=blocks,
            method="rustworkx",
            uniform=False,
        )
        # self.no_faults_success(c, circ, dec, self.model)
        self.correct_all_1(c, circ, dec, self.model)

    def test_d3_5(self):
        """Check 1, zxzxzx, z, pymatching, logical=xyzxyz."""
        blocks = 1
        round_schedule = "zxzxzx"
        basis = "z"
        logical_paulis = "xyzxyz"
        c = HHC(3)
        gen = HHCCircuit(
            c,
            barriers=True,
            idles=True,
            distinct_measurement_idle=True,
            init_error=True,
            group_meas=False,
            xprs=False,
            blocks=blocks,
            round_schedule=round_schedule,
            basis=basis,
            initial_state="+",
            logical_paulis=logical_paulis,
            num_initialize=1,
            idle_before_measure=False,
        )
        circ = gen.syndrome_measurement()
        dec = HHCDecoder(
            n=c.n,
            css_x_gauge_ops=c.x_gauges,
            css_x_stabilizer_ops=c.x_stabilizers,
            css_x_boundary=c.x_boundary,
            css_z_gauge_ops=c.z_gauges,
            css_z_stabilizer_ops=c.z_stabilizers,
            css_z_boundary=c.z_boundary,
            circuit=circ,
            model=self.model,
            basis=basis,
            round_schedule=round_schedule,
            blocks=blocks,
            method="pymatching",
            uniform=False,
        )
        self.no_faults_success(c, circ, dec, self.model)
        self.correct_all_1(c, circ, dec, self.model)

    def test_d3_7(self):
        """Check 3, zx, z, pymatching, -1 eigenstate."""
        blocks = 3
        round_schedule = "zx"
        basis = "z"
        logical_paulis = "ii"
        c = HHC(3)
        gen = HHCCircuit(
            c,
            barriers=True,
            idles=True,
            distinct_measurement_idle=True,
            init_error=True,
            group_meas=False,
            xprs=False,
            blocks=blocks,
            round_schedule=round_schedule,
            basis=basis,
            initial_state="-",
            logical_paulis=logical_paulis,
            num_initialize=1,
            idle_before_measure=False,
        )
        circ = gen.syndrome_measurement()
        dec = HHCDecoder(
            n=c.n,
            css_x_gauge_ops=c.x_gauges,
            css_x_stabilizer_ops=c.x_stabilizers,
            css_x_boundary=c.x_boundary,
            css_z_gauge_ops=c.z_gauges,
            css_z_stabilizer_ops=c.z_stabilizers,
            css_z_boundary=c.z_boundary,
            circuit=circ,
            model=self.model,
            basis=basis,
            round_schedule=round_schedule,
            blocks=blocks,
            method="pymatching",
            uniform=False,
        )
        self.no_faults_success(c, circ, dec, self.model, 1)
        # The propagator method does not treat Paulis in circ as
        # errors, so the observed outcomes are not flipped (0 in arguments)
        self.correct_all_1(c, circ, dec, self.model, 0, False, "propagator")

    def test_d3_10(self):
        """Check 3, zx, x, pymatching."""
        blocks = 3
        round_schedule = "zx"
        basis = "x"
        logical_paulis = "ii"
        c = HHC(3)
        gen = HHCCircuit(
            c,
            barriers=True,
            idles=True,
            distinct_measurement_idle=True,
            init_error=True,
            group_meas=False,
            xprs=False,
            blocks=blocks,
            round_schedule=round_schedule,
            basis=basis,
            initial_state="+",
            logical_paulis=logical_paulis,
            num_initialize=1,
            idle_before_measure=False,
        )
        circ = gen.syndrome_measurement()
        dec = HHCDecoder(
            n=c.n,
            css_x_gauge_ops=c.x_gauges,
            css_x_stabilizer_ops=c.x_stabilizers,
            css_x_boundary=c.x_boundary,
            css_z_gauge_ops=c.z_gauges,
            css_z_stabilizer_ops=c.z_stabilizers,
            css_z_boundary=c.z_boundary,
            circuit=circ,
            model=self.model,
            basis=basis,
            round_schedule=round_schedule,
            blocks=blocks,
            method="pymatching",
            uniform=True,
        )
        self.no_faults_success(c, circ, dec, self.model, 0, True)
        self.correct_all_1(c, circ, dec, self.model, 0, True)


if __name__ == "__main__":
    unittest.main()
