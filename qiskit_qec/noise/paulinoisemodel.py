"""Pauli circuit-level noise model."""
import copy
from typing import Union, Dict
from qiskit_aer import noise


class PauliNoiseModel:
    """Encapsulate a circuit-level Pauli noise model."""

    def __init__(self, fromdict=None):
        """Create a new Pauli noise model.

        Optionally provide a dictionary for the model.
        """
        self.definition = {}
        self.error_probabilities = {}
        self.scale_factors = {}
        self.num_qubits = {}
        self.chan_label = "chan"
        self.p_label = "p"
        self.scale_label = "scale"
        if fromdict is not None:
            for name, propdict in fromdict.items():
                assert self.chan_label in propdict
                self.add_operation(name, propdict[self.chan_label])
                if self.p_label in propdict and propdict[self.p_label] is not None:
                    self.set_error_probability(name, propdict[self.p_label])
                if self.scale_label in propdict and propdict[self.scale_label] is not None:
                    self.set_scale_factor(name, propdict[self.scale_label])

    def to_dict(self):
        """Return corresponding dictionary."""
        asdict = {}
        for name, channel in self.definition.items():
            asdict[name] = {}
            chan = copy.deepcopy(channel)
            asdict[name][self.chan_label] = chan
            if name in self.error_probabilities:
                asdict[name][self.p_label] = self.error_probabilities[name]
            if name in self.scale_factors:
                asdict[name][self.scale_label] = self.scale_factors[name]
        return asdict

    def get_operations(self):
        """Return the list of defined names."""
        return list(self.definition.keys())

    def get_pauli_error_types(self):
        """Return a dict of error types for each operation."""
        error_types = {}
        for k, v in self.definition.items():
            error_types[k] = list(v.keys())
        return error_types

    def get_error_probability(self, name: str):
        """Get the error probability of an operation."""
        if name not in self.error_probabilities:
            raise Exception(f'no error probability for "{name}"')
        return self.error_probabilities[name]

    def get_pauli_weight(self, name: str, paulistring: str):
        """Get the weight for a particular term.

        name = string label for operation
        paulistring = string containing only "i", "x", "y", and "z".
        """
        if name not in self.definition:
            raise Exception(f'"{name}" is not an operation')
        if set(paulistring) > set("ixyz"):
            raise Exception(f'bad paulistring "{paulistring}"')
        if len(paulistring) != self.num_qubits[name]:
            raise Exception(f"paulistring not on {self.num_qubits[name]} qubits")
        if paulistring not in self.definition[name]:
            return 0
        return self.definition[name][paulistring]

    def get_pauli_error_probability(self, name: str, paulistring: str):
        """Get the error probability for a particular term.

        name = string label for operation
        paulistring = string containing only "i", "x", "y", and "z".
        """
        p = self.get_error_probability(name)
        return p * self.get_pauli_weight(name, paulistring)

    def set_scale_factor(self, name: str, factor: Union[int, float]):
        """Assign a scaling factor to an operation."""
        if isinstance(factor, int):
            factor = float(factor)
        if name not in self.definition:
            raise Exception(f'"{name}" is not an operation')
        if factor < 0:
            raise Exception("expected non-negative factor")
        self.scale_factors[name] = factor

    def set_scaled_error_probabilities(self, p: float):
        """Scale and assign error probabilities to operations.

        Only operations with a scale factor property will be changed.
        """
        for name in self.get_operations():
            if name in self.scale_factors:
                self.set_error_probability(name, self.scale_factors[name] * p)

    def set_error_probability(self, name: str, p: float):
        """Assign an error probability to an operation."""
        if name not in self.definition:
            raise Exception(f'"{name}" is not an operation')
        if p < 0 or p > 1:
            raise Exception("expected a probability")
        self.error_probabilities[name] = p

    def add_operation(self, name: str, paulichanneldict: Dict[str, float]):
        """Add a new faulty operation.

        name = string label for operation
        paulichanneldict = dictionary {paulistring: weight, ...}

        Each paulistring contains "i", "x", "y", and "z".
        The weights do not need to be normalized.
        """
        if len(paulichanneldict) == 0:
            raise Exception("expected non-empty dictionary")
        # Check input and compute total weight
        total_weight = 0
        num_qubits = -1
        for k, v in paulichanneldict.items():
            if set(k) > set("ixyz"):
                raise Exception(f'bad paulistring "{k}"')
            if num_qubits == -1:
                num_qubits = len(k)
            elif len(k) != num_qubits:
                raise Exception(f"expected paulistring on {num_qubits} qubits")
            total_weight += v
        # Add definition
        self.definition[name] = copy.deepcopy(paulichanneldict)
        self.num_qubits[name] = num_qubits
        # Normalize weights
        for k in paulichanneldict.keys():
            self.definition[name][k] /= total_weight

    def as_aer_noise_model(self):
        """Return corresponding Aer noise model."""
        model = noise.noise_model.NoiseModel()
        for name, paulichanneldict in self.definition.items():
            if name not in self.error_probabilities:
                raise Exception(f'no error probability for "{name}"')
            p = self.error_probabilities[name]
            terms = [("I" * self.num_qubits[name], 1.0 - p)]
            for pauli, weight in paulichanneldict.items():
                terms.append((pauli.upper(), weight * p))
            # Create the noise operator and add to the model
            if name == "measure":
                # Measurement is a special case
                mat = [[1.0 - p, p], [p, 1.0 - p]]
                error = noise.errors.readout_error.ReadoutError(mat)
                model.add_all_qubit_readout_error(error)
            else:
                # Everything else handled here
                error = noise.errors.pauli_error(terms)
                model.add_all_qubit_quantum_error(error, [name])
        return model
