"""Test fault enumerator."""
import unittest
from qiskit import QuantumCircuit
from qiskit.circuit.library import IGate
from qiskit_qec.analysis.faultsampler import FaultSampler
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class TestFaultSampler(unittest.TestCase):
    """Tests of fault sampler class."""

    def test_fault_sampler_stabilizer(self):
        """Test sampling multiple faults using the stabilizer simulator."""
        qc = QuantumCircuit(2, 2)
        qc.reset(0)
        qc.reset(1)
        qc.h(0)
        qc.cx(0, 1)
        qc.append(IGate(label="special"), [0])
        qc.measure(0, 0)
        qc.measure(1, 1)
        pnm = PauliNoiseModel()
        p = 0.2
        pnm.add_operation("h", {"x": 1, "y": 1, "z": 1})
        pnm.add_operation("cx", {"ix": 1, "xi": 1, "xx": 1})
        pnm.add_operation("special", {"x": 1, "y": 1})
        pnm.add_operation("measure", {"x": 1})
        pnm.add_operation("reset", {"x": 1})
        pnm.set_error_probability("h", p)
        pnm.set_error_probability("cx", p)
        pnm.set_error_probability("special", p)
        pnm.set_error_probability("measure", 0.0)
        pnm.set_error_probability("reset", 0.0)
        fs = FaultSampler(qc, model=pnm, method="stabilizer", sim_seed=0)
        num_samples = 1000
        fault_paths = list(fs.sample(num_samples))
        # Test probability circuit output is incorrect
        failures = 0
        for event in fault_paths:
            if event[3] == [0, 1] or event[3] == [1, 0]:
                failures += 1
        # cx fails with ix or xi OR special fails OR cx fails with xx and special fails
        pfail = 5.0 / 3.0 * p * (1 - p) + 1.0 / 3.0 * p**2
        self.assertAlmostEqual(failures / num_samples, pfail, 1)
        # Test probability of each number of faults
        num_faults = {}
        for x in fault_paths:
            failures = len(x[1])
            if failures not in num_faults:
                num_faults[failures] = 1
            else:
                num_faults[failures] += 1
        for key in num_faults:
            num_faults[key] /= num_samples
        limit = {0: (1 - p) ** 3, 1: 3 * p * (1 - p) ** 2, 2: 3 * p**2 * (1 - p), 3: p**3}
        for key, value in num_faults.items():
            self.assertAlmostEqual(value, limit[key], 1)

    def test_fault_sampler_propagator(self):
        """Test sampling multiple faults using the error propagator."""
        qc = QuantumCircuit(2, 2)
        qc.reset(0)
        qc.reset(1)
        qc.h(0)
        qc.cx(0, 1)
        qc.append(IGate(label="special"), [0])
        qc.measure(0, 0)
        qc.measure(1, 1)
        pnm = PauliNoiseModel()
        p = 0.2
        pnm.add_operation("h", {"x": 1, "y": 1, "z": 1})
        pnm.add_operation("cx", {"ix": 1, "xi": 1, "xx": 1})
        pnm.add_operation("special", {"x": 1, "y": 1})
        pnm.add_operation("measure", {"x": 1})
        pnm.add_operation("reset", {"x": 1})
        pnm.set_error_probability("h", p)
        pnm.set_error_probability("cx", p)
        pnm.set_error_probability("special", p)
        pnm.set_error_probability("measure", 0.0)
        pnm.set_error_probability("reset", 0.0)
        fs = FaultSampler(qc, model=pnm, method="propagator", sim_seed=0)
        num_samples = 10000
        fault_paths = list(fs.sample(num_samples))
        # Test probability circuit output is incorrect
        failures = 0
        for event in fault_paths:
            if event[3] == [0, 1] or event[3] == [1, 0]:
                failures += 1
        # cx fails with ix or xi OR special fails OR cx fails with xx and special fails
        pfail = 5.0 / 3.0 * p * (1 - p) + 1.0 / 3.0 * p**2
        self.assertAlmostEqual(failures / num_samples, pfail, 1)
        # Test probability of each number of faults
        num_faults = {}
        for x in fault_paths:
            failures = len(x[1])
            if failures not in num_faults:
                num_faults[failures] = 1
            else:
                num_faults[failures] += 1
        for key in num_faults:
            num_faults[key] /= num_samples
        limit = {0: (1 - p) ** 3, 1: 3 * p * (1 - p) ** 2, 2: 3 * p**2 * (1 - p), 3: p**3}
        for key, value in num_faults.items():
            self.assertAlmostEqual(value, limit[key], 1)

    def test_fault_sampler_stabilizer_2(self):
        """Test sampling in a circuit with one faulty location."""
        qc = QuantumCircuit(1, 1)
        qc.reset(0)
        qc.id(0)
        qc.measure(0, 0)
        model = PauliNoiseModel()
        model.add_operation("id", {"x": 1, "y": 1, "z": 1})
        model.set_error_probability("id", 0.2)
        fs = FaultSampler(qc, model=model, method="stabilizer", sim_seed=0)
        # pylint: disable=c-extension-no-member
        block = fs.sample(1000)
        failures = 0
        for x in block:
            if x[3][0] == 1:
                failures += 1
        failures /= 1000
        self.assertAlmostEqual(failures, 0.2 * 2.0 / 3.0, 1)

    def test_fault_sampler_propagator_2(self):
        """Test sampling in a circuit with one faulty location."""
        qc = QuantumCircuit(1, 1)
        qc.reset(0)
        qc.h(0)
        qc.measure(0, 0)
        model = PauliNoiseModel()
        model.add_operation("h", {"x": 1, "y": 1, "z": 1})
        model.set_error_probability("h", 0.2)
        fs = FaultSampler(qc, model=model, method="propagator", sim_seed=0)
        # pylint: disable=c-extension-no-member
        block = fs.sample(10000)
        failures = 0
        for x in block:
            if x[3][0] == 1:
                failures += 1
        failures /= 10000
        self.assertAlmostEqual(failures, 0.2 * 2.0 / 3.0, 1)


if __name__ == "__main__":
    unittest.main()
