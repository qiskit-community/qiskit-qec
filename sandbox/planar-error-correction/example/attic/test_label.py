from qiskit import execute, Aer, QuantumCircuit
from qiskit.providers.aer import noise
from qiskit.providers.aer.noise.errors import *

from qiskit.circuit.library import IGate

import logging

logging.basicConfig(filename="output.log", filemode="w", level=logging.INFO)


def make_depolarizing_noise_model(p1, p2, q):
    """Create an Aer noise model for circuit depolarizing noise.

    p1 = error probability for h, id
    p2 = error probability for cx
    q = error probability for measure, reset

    Return the noise model object.
    """
    model = noise.noise_model.NoiseModel()
    perror1 = noise.errors.pauli_error(
        [("I", 1.0 - p1), ("X", p1 / 3.0), ("Y", p1 / 3.0), ("Z", p1 / 3.0)]
    )
    perrorX = noise.errors.pauli_error([("I", 1.0 - q), ("X", q)])
    perror2 = noise.errors.pauli_error(
        [
            ("II", 1.0 - p2),
            ("IX", p2 / 15.0),
            ("IY", p2 / 15.0),
            ("IZ", p2 / 15.0),
            ("XI", p2 / 15.0),
            ("XX", p2 / 15.0),
            ("XY", p2 / 15.0),
            ("XZ", p2 / 15.0),
            ("YI", p2 / 15.0),
            ("YX", p2 / 15.0),
            ("YY", p2 / 15.0),
            ("YZ", p2 / 15.0),
            ("ZI", p2 / 15.0),
            ("ZX", p2 / 15.0),
            ("ZY", p2 / 15.0),
            ("ZZ", p2 / 15.0),
        ]
    )
    rerror = noise.errors.readout_error.ReadoutError([[1.0 - q, q], [q, 1.0 - q]])
    model.add_all_qubit_quantum_error(perror1, ["h", "im", "x", "y", "z"])
    model.add_all_qubit_quantum_error(perror2, ["cx"])
    model.add_all_qubit_quantum_error(perrorX, ["reset"])
    model.add_all_qubit_readout_error(rerror)
    return model


# This seems to define an identity gate that carries a Qiskit label,
# which can be used to single that gate out in a noise model.
# When the circuit is drawn, the label is used for the gate name.
# Q: How do I change the QASM representation? It is still "i".
qc = QuantumCircuit(1, 1)
qc.i(0)
qc.i(0)
qc.i(0)
qc.i(0)
qc.i(0)
IGateM = IGate(label="im")
qc.append(IGateM, [0])
qc.measure(0, 0)
print(qc.qasm())
print(qc.draw())
model = make_depolarizing_noise_model(0.2, 0.2, 0.2)
result_noise = execute(
    qc,
    Aer.get_backend("qasm_simulator"),
    noise_model=model,
    basis_gates=model.basis_gates,
    method="stabilizer",
    shots=10,
    optimization_level=0,
    seed_simulator=19,
).result()
print(result_noise.to_dict())
counts_noise = result_noise.get_counts(qc)
print(counts_noise)
