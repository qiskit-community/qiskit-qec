"""Test fault enumerator."""
import unittest
import itertools
from qiskit import QuantumCircuit
from qiskit.circuit.library import IGate
from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class TestFaultEnumerator(unittest.TestCase):
    """Tests of fault enumerator class."""

    def test_fault_enumerator_stabilizer(self):
        """Test an example using the stabilizer simulator."""
        qc = QuantumCircuit(3, 3)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(1, 2)
        qc.append(IGate(label="special"), [0])
        qc.measure(0, 0)
        qc.measure(1, 1)
        qc.measure(2, 2)
        pnm = PauliNoiseModel()
        pnm.add_operation("h", {"x": 1, "y": 1, "z": 1})
        pnm.add_operation("cx", {"ix": 1, "xi": 1, "xx": 1})
        pnm.add_operation("special", {"y": 1})
        pnm.add_operation("measure", {"x": 1})
        fe = FaultEnumerator(qc, order=1, method="stabilizer", model=pnm, sim_seed=0)
        fault_paths = list(fe.generate())
        self.assertEqual(
            fault_paths,
            [
                (0, ["h"], ["x"], [0, 0, 0]),
                (1, ["h"], ["y"], [0, 0, 0]),
                (2, ["h"], ["z"], [0, 0, 0]),
                (3, ["cx"], ["ix"], [0, 1, 1]),
                (4, ["cx"], ["xi"], [0, 1, 1]),
                (5, ["cx"], ["xx"], [0, 0, 0]),
                (6, ["special"], ["y"], [0, 1, 1]),
                (7, ["cx"], ["ix"], [0, 0, 1]),
                (8, ["cx"], ["xi"], [0, 1, 0]),
                (9, ["cx"], ["xx"], [0, 1, 1]),
                (10, ["measure"], ["x"], [0, 1, 1]),
                (11, ["measure"], ["x"], [0, 1, 0]),
                (12, ["measure"], ["x"], [0, 0, 1]),
            ],
        )

    def test_fault_enumerator_propagator(self):
        """Test an example using the error propagator."""
        qc = QuantumCircuit(3, 3)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(1, 2)
        qc.append(IGate(label="special"), [0])
        qc.measure(0, 0)
        qc.measure(1, 1)
        qc.measure(2, 2)
        pnm = PauliNoiseModel()
        pnm.add_operation("h", {"x": 1, "y": 1, "z": 1})
        pnm.add_operation("cx", {"ix": 1, "xi": 1, "xx": 1})
        pnm.add_operation("special", {"y": 1})
        pnm.add_operation("measure", {"x": 1})
        fe = FaultEnumerator(qc, order=1, method="propagator", model=pnm, sim_seed=0)
        fault_paths = list(fe.generate())
        self.assertEqual(
            fault_paths,
            [
                (0, ["h"], ["x"], [1, 1, 1]),
                (1, ["h"], ["y"], [1, 1, 1]),
                (2, ["h"], ["z"], [0, 0, 0]),
                (3, ["cx"], ["ix"], [0, 1, 1]),
                (4, ["cx"], ["xi"], [1, 0, 0]),
                (5, ["cx"], ["xx"], [1, 1, 1]),
                (6, ["special"], ["y"], [1, 0, 0]),
                (7, ["cx"], ["ix"], [0, 0, 1]),
                (8, ["cx"], ["xi"], [0, 1, 0]),
                (9, ["cx"], ["xx"], [0, 1, 1]),
                (10, ["measure"], ["x"], [1, 0, 0]),
                (11, ["measure"], ["x"], [0, 1, 0]),
                (12, ["measure"], ["x"], [0, 0, 1]),
            ],
        )

    def test_fault_enumerator_stabilizer_blocks(self):
        """Test an example using stabilizer simulator and returning blocks."""
        qc = QuantumCircuit(3, 3)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(1, 2)
        qc.append(IGate(label="special"), [0])
        qc.measure(0, 0)
        qc.measure(1, 1)
        qc.measure(2, 2)
        pnm = PauliNoiseModel()
        pnm.add_operation("h", {"x": 1, "y": 1, "z": 1})
        pnm.add_operation("cx", {"ix": 1, "xi": 1, "xx": 1})
        pnm.add_operation("special", {"y": 1})
        pnm.add_operation("measure", {"x": 1})
        fe = FaultEnumerator(qc, order=1, method="stabilizer", model=pnm, sim_seed=0)
        blocks = list(fe.generate_blocks(3))
        fault_paths = list(itertools.chain(*blocks))
        self.assertEqual(
            fault_paths,
            [
                (0, ["h"], ["x"], [0, 0, 0]),
                (1, ["h"], ["y"], [0, 0, 0]),
                (2, ["h"], ["z"], [0, 0, 0]),
                (3, ["cx"], ["ix"], [0, 1, 1]),
                (4, ["cx"], ["xi"], [0, 1, 1]),
                (5, ["cx"], ["xx"], [0, 0, 0]),
                (6, ["special"], ["y"], [0, 1, 1]),
                (7, ["cx"], ["ix"], [0, 0, 1]),
                (8, ["cx"], ["xi"], [0, 1, 0]),
                (9, ["cx"], ["xx"], [0, 1, 1]),
                (10, ["measure"], ["x"], [0, 1, 1]),
                (11, ["measure"], ["x"], [0, 1, 0]),
                (12, ["measure"], ["x"], [0, 0, 1]),
            ],
        )

    def test_fault_enumerator_propagator_blocks(self):
        """Test an example with the error propagator and returning blocks."""
        qc = QuantumCircuit(3, 3)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(1, 2)
        qc.append(IGate(label="special"), [0])
        qc.measure(0, 0)
        qc.measure(1, 1)
        qc.measure(2, 2)
        pnm = PauliNoiseModel()
        pnm.add_operation("h", {"x": 1, "y": 1, "z": 1})
        pnm.add_operation("cx", {"ix": 1, "xi": 1, "xx": 1})
        pnm.add_operation("special", {"y": 1})
        pnm.add_operation("measure", {"x": 1})
        fe = FaultEnumerator(qc, order=1, method="propagator", model=pnm, sim_seed=0)
        blocks = list(fe.generate_blocks(3))
        fault_paths = list(itertools.chain(*blocks))
        self.assertEqual(
            fault_paths,
            [
                (0, ["h"], ["x"], [1, 1, 1]),
                (1, ["h"], ["y"], [1, 1, 1]),
                (2, ["h"], ["z"], [0, 0, 0]),
                (3, ["cx"], ["ix"], [0, 1, 1]),
                (4, ["cx"], ["xi"], [1, 0, 0]),
                (5, ["cx"], ["xx"], [1, 1, 1]),
                (6, ["special"], ["y"], [1, 0, 0]),
                (7, ["cx"], ["ix"], [0, 0, 1]),
                (8, ["cx"], ["xi"], [0, 1, 0]),
                (9, ["cx"], ["xx"], [0, 1, 1]),
                (10, ["measure"], ["x"], [1, 0, 0]),
                (11, ["measure"], ["x"], [0, 1, 0]),
                (12, ["measure"], ["x"], [0, 0, 1]),
            ],
        )


if __name__ == "__main__":
    unittest.main()
