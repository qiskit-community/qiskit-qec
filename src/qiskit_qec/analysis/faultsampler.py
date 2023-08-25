"""Quantum circuit fault path sampler."""
from typing import Tuple, List

import numpy as np

from qiskit import QuantumCircuit
from qiskit.dagcircuit.dagnode import DAGNode
from qiskit.converters import circuit_to_dag
from qiskit.circuit.library import IGate, XGate, YGate, ZGate
from qiskit_aer import Aer

from qiskit_qec.utils.dag import node_name_label

from qiskit_qec.analysis.extensions import C_FAULT_SAMPLER

from qiskit_qec.analysis.errorpropagator import ErrorPropagator
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel
from qiskit_qec.exceptions import QiskitQECError

if C_FAULT_SAMPLER:
    from qiskit_qec.analysis.extensions import _CFaultSampler


class FaultSampler:
    """Samples faults in a circuit according to a noise model."""

    METHOD_STABILIZER: str = "stabilizer"
    METHOD_PROPAGATOR: str = "propagator"
    AVAILABLE_METHODS = {METHOD_STABILIZER, METHOD_PROPAGATOR}

    def __init__(
        self,
        circ: QuantumCircuit,
        model: PauliNoiseModel,
        method: str = METHOD_PROPAGATOR,
        sim_seed: int = 0,
    ):
        """Construct a fault sampler object.

        circ = QuantumCircuit object
        model = PauliNoiseModel object
        method = simulation method

        Faulty operations are identified by a string
        containing the Qiskit gate name or label
        if the latter exists.
        """
        self.model = model
        if method not in self.AVAILABLE_METHODS:
            raise QiskitQECError("fmethod {methid} is not supported.")
        self.method = method

        self.location_types = self.model.get_operations()
        self.pauli_error_types = self.model.get_pauli_error_types()
        self.sim_seed = sim_seed
        self._process_circuit(circ)
        self.faulty_nodes = list(filter(lambda x: x[1], self.tagged_nodes))
        self.faulty_ops_indices = [x[2] for x in self.faulty_nodes]
        self.faulty_ops_labels = [node_name_label(x[0]) for x in self.faulty_nodes]
        self.label_to_pauli_weight = {
            label: {
                pauli: self.model.get_pauli_weight(label, pauli)
                for pauli in self.pauli_error_types[label]
            }
            for label in self.faulty_ops_labels
        }
        self.label_to_pauli_weight_tuple = {
            label: list(pdict.items()) for label, pdict in self.label_to_pauli_weight.items()
        }
        self.label_to_error_probability = {
            label: self.model.get_error_probability(label) for label in self.faulty_ops_labels
        }
        self.propagator = None
        self.use_compiled = False
        self.faultsamp = None
        if method == "propagator":
            self.propagator = ErrorPropagator()  # pylint: disable=abstract-class-instantiated
            self.propagator.load_circuit(circ)
            self.reg_sizes = [len(reg) for reg in circ.cregs]
            if C_FAULT_SAMPLER:
                self.use_compiled = True
                self.faultsamp = _CFaultSampler(
                    len(self.propagator.get_error()),
                    len(self.propagator.get_bit_array()),
                    self.propagator.encoded_circ,  # pylint: disable=no-member
                    self.faulty_ops_indices,
                    self.faulty_ops_labels,
                    self.label_to_pauli_weight_tuple,
                    self.label_to_error_probability,
                    self.sim_seed,
                )

    def _process_circuit(self, circ: QuantumCircuit):
        """Precompute some data about a circuit."""
        self.dag = circuit_to_dag(circ)
        # Construct list of tuples of topologically sorted nodes
        # (node, faulty?, index)
        self.tagged_nodes = []
        index = 0
        for node in self.dag.topological_op_nodes():
            if node_name_label(node) in self.location_types:
                self.tagged_nodes.append((node, True, index))
            else:
                self.tagged_nodes.append((node, False, index))
            index += 1

    def _faulty_circuit(self, comb: Tuple[DAGNode], error: Tuple[str]) -> QuantumCircuit:
        """Construct faulty QuantumCircuit with the given faults.

        comb = tuple of DAG nodes
        error = tuple of Pauli strings
        Return QuantumCircuit
        """
        error_to_gate = {"i": IGate(), "x": XGate(), "y": YGate(), "z": ZGate()}
        circ = QuantumCircuit(
            *self.dag.qregs.values(),
            *self.dag.cregs.values(),
            name=self.dag.name,
            global_phase=self.dag.global_phase,
        )
        circ.calibrations = self.dag.calibrations
        for orig_node in self.dag.topological_op_nodes():
            if orig_node in comb:
                index = comb.index(orig_node)
                errs = [error_to_gate[x] for x in error[index]]
                if orig_node.name == "measure":
                    for i, err in enumerate(errs):
                        circ._append(err, [orig_node.qargs[i]], orig_node.cargs)
            inst = orig_node.op.copy()
            circ._append(inst, orig_node.qargs, orig_node.cargs)
            if orig_node in comb:
                if orig_node.name != "measure":
                    for i, err in enumerate(errs):
                        circ._append(err, [orig_node.qargs[i]], orig_node.cargs)
        circ.duration = self.dag.duration
        circ.unit = self.dag.unit
        return circ

    def _stabilizer_simulation(self, circ: QuantumCircuit) -> List[int]:
        """Run stabilizer simulation.

        circ = input QuantumCircuit
        Return REVERSED outcome.
        """

        def gint(c):
            # Casts to int if possible
            if c.isnumeric():
                return int(c)
            else:
                return c

        backend = Aer.get_backend("aer_simulator")
        options = {"method": "stabilizer", "shots": 1, "seed_simulator": self.sim_seed}
        result = backend.run(circ, **options).result()
        outcomes = result.get_counts(circ)
        raw_outcome = list(outcomes.keys())[0]
        outcome = list(map(gint, raw_outcome[::-1]))
        return outcome

    def sample_one(self):
        """Sample faults in the quantum circuit.

        Note: The bit order of outcome is reversed from Qiskit.

        Returns (index, labels, error, outcome).
        """
        # sample faulty nodes according to model
        failed_nodes = []
        errors = []
        indices = []
        uniform_samples = np.random.uniform(low=0.0, high=1.0, size=len(self.faulty_nodes))
        for i, node_dat in enumerate(self.faulty_nodes):
            label = node_name_label(node_dat[0])
            if uniform_samples[i] < self.label_to_error_probability[label]:
                failed_nodes.append(node_dat[0])
                indices.append(node_dat[2])
                paulis = []
                probs = []
                for pauli, prob in self.label_to_pauli_weight[label].items():
                    paulis.append(pauli)
                    probs.append(prob)
                errors.append("".join(np.random.choice(paulis, 1, p=probs)))
        labels = [node_name_label(x) for x in failed_nodes]
        if self.method == "stabilizer":
            fcirc = self._faulty_circuit(failed_nodes, errors)
            outcome = self._stabilizer_simulation(fcirc)
        elif self.method == "propagator":
            raw_outcome = self.propagator.propagate_faults(indices, errors)
            if len(self.reg_sizes) == 1:
                outcome = raw_outcome
            else:
                outcome = []
                j = 0
                for reg_size in self.reg_sizes:
                    for _ in range(reg_size):
                        outcome.append(raw_outcome[j])
                        j += 1
                    outcome.append(" ")
                outcome.pop(-1)
        return (-1, labels, errors, outcome)

    def sample(self, blocksize: int = 10000):
        """Sample faults in a quantum circuit.

        blocksize = number of samples to produce per call

        Note: The bit order of outcome is reversed from Qiskit.

        Returns [(index, labels, error, outcome), ...] with approximately
        blocksize elements. May be slightly more per call, or less on the
        final call.
        """
        if self.use_compiled:
            return self.faultsamp.sample(blocksize)
        else:
            block = []
            count = 0
            while count < blocksize:
                block.append(self.sample_one())
                count += 1
            return block
