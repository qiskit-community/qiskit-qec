from qiskit import execute, Aer
from qiskit.providers.aer import noise

from plerco.pec_python.qec_code.hhc import HHC
from plerco.pec_python.qec_circuit.hhc_circuit import HHCCircuit
from plerco.pec_python.qec_decoder.hhc_decoder import HHCDecoder

from qiskit.converters import circuit_to_dag

import numpy as np

import time


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
    perrorx = noise.errors.pauli_error([("I", 1.0 - q), ("X", q)])
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
    model.add_all_qubit_quantum_error(perror1, ["h", "id"])
    model.add_all_qubit_quantum_error(perror2, ["cx"])
    model.add_all_qubit_quantum_error(perrorx, ["reset"])
    model.add_all_qubit_readout_error(rerror)
    return model


def run_circuit(qc, shots, model=None, seed=100):
    """Simulate a circuit with the stabilizer simulator.

    qc = QuantumCircuit to simulate
    shots = number of shots to run
    model = noise model to use, None if no noise
    seed = seed value for random number generator

    Return the dictionary of counts.
    """
    qasmsim = Aer.get_backend("qasm_simulator")
    if model is None:
        result = execute(
            qc,
            qasmsim,
            method="stabilizer",
            shots=shots,
            optimization_level=0,
            seed_simulator=seed,
        ).result()
    else:
        result = execute(
            qc,
            qasmsim,
            noise_model=model,
            basis_gates=model.basis_gates,
            method="stabilizer",
            shots=shots,
            optimization_level=0,
            seed_simulator=seed,
        ).result()
    return result.get_counts(qc)


def sigma(p, shots):
    """Compute the sample standard deviation.

    p = Bernoulli parameter
    shots = number of samples

    Return the value.
    """
    return np.sqrt(p * (1 - p) / float(shots - 1))


def decode_outcome(outcome, num, decoder):
    """Decode an outcome to decide failure.

    outcome = outcome string from simulator
    num = number of occurrences of this outcome
    decoder = decoder object

    Return number of failures associated to this outcome.
    """
    # convert to list of integers
    reversed_outcome = list(map(int, outcome[::-1]))
    corrected_outcomes = decoder.process(reversed_outcome)
    if basis == "z":
        fail = decoder.code.logical_x_error(corrected_outcomes)
    elif basis == "x":
        fail = decoder.code.logical_z_error(corrected_outcomes)
    if bool(fail[0]):
        return num
    else:
        return 0


def sim(circ, shots, p, decoder, code, param, basis):
    """Simulate one circuit at error p.

    circ = quantum circuit
    shots = number of shots
    p = physical error rate
    decoder = decoder object
    code = code object
    param = list of gate error parameter values
    basis = 'x' or 'z' prep/meas basis

    Return pfail, sigma
    """
    tic = time.perf_counter()
    model = make_depolarizing_noise_model(p, p, p)
    counts_n = run_circuit(circ, shots, model)
    toc = time.perf_counter()
    print(f"simulate circuits in {toc - tic:0.4f} seconds", flush=True)
    decoder.update_edge_weights(param)
    tic = time.perf_counter()
    failures = sum(map(lambda x: decode_outcome(x[0], x[1], decoder), counts_n.items()))
    toc = time.perf_counter()
    print(f"decode outcomes in {toc - tic:0.4f} seconds", flush=True)
    pfail = float(failures) / float(shots)
    sig = sigma(pfail, shots)
    print("pfail = %f +- %f" % (pfail, sig), flush=True)
    return pfail, sig


desc = "hhc"
d = 3
rounds = 1
round_schedule = "xzxzx"
basis = "z"
shots = 2

pin = [0.001]

print("desc = %s" % desc, flush=True)
print("d = %d" % d, flush=True)
print("rounds = %d" % rounds, flush=True)
print("round_schedule = %s" % round_schedule, flush=True)
print("basis = %s" % basis, flush=True)
print("shots = %d" % shots, flush=True)

code = HHC(d)
gen = HHCCircuit(code)
circ = gen.syndrome_measurement(rounds, round_schedule, basis)
print("generated circuit")

dag = circuit_to_dag(circ)
qr = list(dag.qregs.values())[0]
cr = list(dag.cregs.values())[0]
print("q = QuantumRegister(%d, '%s')" % (qr.size, qr.name))
print("c = ClassicalRegister(%d, '%s')" % (cr.size, cr.name))
print('circ = QuantumCircuit(q, c, name="hhc3")')
for orig_node in dag.topological_op_nodes():
    if orig_node.name == "h":
        print("circ.h(q[%d])" % orig_node.qargs[0].index)
    elif orig_node.name == "cx":
        print(
            "circ.cx(q[%d], q[%d])"
            % (orig_node.qargs[0].index, orig_node.qargs[1].index)
        )
    elif orig_node.name == "reset":
        print("circ.reset(q[%d])" % orig_node.qargs[0].index)
    elif orig_node.name == "measure":
        print(
            "circ.measure(q[%d], c[%d])"
            % (orig_node.qargs[0].index, orig_node.cargs[0].index)
        )
    elif orig_node.name == "id":
        print("circ.i(q[%d])" % orig_node.qargs[0].index)
    elif orig_node.name == "barrier":
        # We know that for this circuit all barriers are barrier(q)
        print("circ.barrier(q)")
    else:
        raise Exception("unknown %s" % orig_node.name)
    if not orig_node.condition is None:
        raise Exception("orig_node.condition is not None")

tic = time.perf_counter()
decoder = HHCDecoder(
    code,
    rounds,
    round_schedule,
    basis,
    circuit=circ,
    parameters=["cx", "id", "reset", "measure", "h"],
    usepymatching=True,
)
toc = time.perf_counter()
print(f"construct decoder in {toc - tic:0.4f} seconds", flush=True)

print(decoder.g.nodes(data=True))
print(decoder.g.edges(data=True))

print("cx = Symbol('cx')")
print("id = Symbol('id')")
print("reset = Symbol('reset')")
print("h = Symbol('h')")
print("measure = Symbol('measure')")
qubit_map = {}
for ed in decoder.g.edges(data=True):
    if len(ed[2]["qubits"]) > 0:
        qubit_map[ed[2]["qubits"][0]] = ed[2]["qubit_id"]
print("qubit_map = %s" % qubit_map)
print("g = nx.Graph()")
for ed in decoder.g.edges(data=True):
    if "weight_poly" in ed[2]:
        print(
            "g.add_edge(%d, %d, qubit_id=%d, weight=%d, weight_poly=%s)"
            % (ed[0], ed[1], ed[2]["qubit_id"], ed[2]["weight"], ed[2]["weight_poly"])
        )
    else:
        print(
            "g.add_edge(%d, %d, qubit_id=%d, weight=%d)"
            % (ed[0], ed[1], ed[2]["qubit_id"], ed[2]["weight"])
        )
for no in decoder.g.nodes(data=True):
    if "is_boundary" in no[1] and no[1]["is_boundary"]:
        print("g.nodes[%d]['is_boundary'] = True" % no[0])

print(decoder.idxmap)
print(decoder.layer_types)

print("xg = ", decoder.code.x_gauges)
print("zg = ", decoder.code.z_gauges)
print("xs = ", decoder.code.x_stabilizers)
print("zs = ", decoder.code.z_stabilizers)
print("lz = ", decoder.code.logical_z)
print("lx = ", decoder.code.logical_x)

print("xboundary = ", decoder.code.x_boundary)
print("zboundary = ", decoder.code.z_boundary)

print("zsmat = ", decoder.code._to_matrix(decoder.code.z_stabilizers))
print("xsmat = ", decoder.code._to_matrix(decoder.code.x_stabilizers))
print("zlmat = ", decoder.code._to_matrix(decoder.code.logical_z))
print("xlmat = ", decoder.code._to_matrix(decoder.code.logical_x))

pfails = []
sigs = []
for p in pin:
    print("running %d %f" % (shots, p), flush=True)
    pf, s = sim(circ, shots, p, decoder, code, [p, p, p, p, p], basis)
    pfails.append(pf)
    sigs.append(s)

print("pfails = %s" % pfails, flush=True)
print("sigs = %s" % sigs, flush=True)
