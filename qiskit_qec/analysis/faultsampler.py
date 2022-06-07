"""Quantum circuit fault path sampler."""
from itertools import combinations, product
from typing import Tuple

from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag
from qiskit.circuit.library import IGate, XGate, YGate, ZGate
from qiskit import execute, Aer

from qiskit_qec.analysis.epselector import EPSelector
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel
from qiskit_qec.exceptions import QiskitQECError


class FaultSampler:
    """Samples faults in a circuit according to a noise model."""

    METHOD_STABILIZER: str = "stabilizer"
    METHOD_PROPAGATOR: str = "propagator"
    AVAILABLE_METHODS = {METHOD_STABILIZER, METHOD_PROPAGATOR}

    def __init__(
        self,
        circ,
        model,
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
        self.propagator = None
        self.use_compiled = False
        self.faultsamp = None
        if method == "propagator":
            eps = EPSelector()
            self.propagator = eps.get_error_propagator()
            self.propagator.load_circuit(circ)
            self.reg_sizes = [len(reg) for reg in circ.cregs]
            try:
                from qiskit_qec.extensions.compiledextension import (
                    FaultSampler as compiledFaultSampler,
                )
            except ImportError:
                self.use_compiled = False
            else:
                self.use_compiled = True
                #### TODO: get the right data into the constructor
                faulty_nodes = filter(lambda x: x[1], self.tagged_nodes)
                faulty_ops_indices = [x[2] for x in faulty_nodes]
                faulty_nodes = filter(lambda x: x[1], self.tagged_nodes)
                faulty_ops_labels = [self._node_name_label(x[0]) for x in faulty_nodes]
                faulty_ops_pauli_errors = [self.pauli_error_types[x] for x in faulty_ops_labels]

                self.faultsamp = compiledFaultSampler(
                    len(self.propagator.get_error()),
                    len(self.propagator.get_bit_array()),
                    self.propagator.encoded_circ,
                    faulty_ops_indices,
                    faulty_ops_labels,
                    faulty_ops_pauli_errors,
                )

    def _node_name_label(self, node):
        """Form an identifier string for a node's operation.

        Use node.op._label if it exists. Otherwise use node.name.
        Return a string.
        """
        if "_label" in node.op.__dict__ and node.op._label is not None:
            name_label = node.op._label
        else:
            name_label = node.name
        return name_label

    def _process_circuit(self, circ):
        """Precompute some data about a circuit."""
        self.dag = circuit_to_dag(circ)
        # Construct list of tuples of topologically sorted nodes
        # (node, faulty?, index)
        self.tagged_nodes = []
        index = 0
        for node in self.dag.topological_op_nodes():
            if self._node_name_label(node) in self.location_types:
                self.tagged_nodes.append((node, True, index))
            else:
                self.tagged_nodes.append((node, False, index))
            index += 1

    def _faulty_circuit(self, comb, error: Tuple[str]):
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

    def _stabilizer_simulation(self, circ):
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

        result = execute(
            circ,
            Aer.get_backend("aer_simulator"),
            method="stabilizer",
            shots=1,
            optimization_level=0,
            seed_simulator=self.sim_seed,
        ).result()
        outcomes = result.get_counts(circ)
        raw_outcome = list(outcomes.keys())[0]
        outcome = list(map(gint, raw_outcome[::-1]))
        return outcome

    def generate(self):
        """Generator to iterate over faults in a quantum circuit.

        Note: The bit order of outcome is reversed from Qiskit.

        Yields (index, labels, error, outcome).
        """
        index = 0
        if self.method == "stabilizer":
            faulty_nodes = filter(lambda x: x[1], self.tagged_nodes)
            for comb in combinations(faulty_nodes, self.order):
                nodes = [x[0] for x in comb]
                labels = [self._node_name_label(x) for x in nodes]
                iterable = [self.pauli_error_types[x] for x in labels]
                for error in product(*iterable):
                    fcirc = self._faulty_circuit(nodes, error)
                    outcome = self._stabilizer_simulation(fcirc)
                    yield (index, labels, list(error), outcome)
                    index += 1
        elif self.method == "propagator":
            faulty_nodes = filter(lambda x: x[1], self.tagged_nodes)
            for comb in combinations(faulty_nodes, self.order):
                nodes = [x[0] for x in comb]
                labels = [self._node_name_label(x) for x in nodes]
                indices = [x[2] for x in comb]
                iterable = [self.pauli_error_types[x] for x in labels]
                for error in product(*iterable):
                    raw_outcome = self.propagator.propagate_faults(indices, error)
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
                    yield (index, labels, list(error), outcome)
                    index += 1

    def generate_blocks(self, blocksize: int = 10000):
        """Generator to iterate over sequences of faults in a quantum circuit.

        blocksize = number of faults to process per call

        Note: The bit order of outcome is reversed from Qiskit.

        Yields [(index, labels, error, outcome), ...] with approximately
        blocksize elements. May be slightly more per call, or less on the
        final call.
        """
        if self.use_compiled:
            while not self.faultenum.done():
                block = self.faultenum.enumerate(blocksize)
                yield block
        else:
            # Fall back to calling the generate() method repeatedly
            block = []
            count = 0
            for x in self.generate():
                block.append(x)
                count += 1
                if count >= blocksize:
                    yield block
                    block = []
                    count = 0
            yield block
