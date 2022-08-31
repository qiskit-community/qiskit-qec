from qiskit import Aer

from qiskit_qec.circuits.repetition_code import RepetitionCodeCircuit
from qiskit_qec.decoders import DecodingGraph
from qiskit_qec.decoders.repetition_decoder import RepetitionDecoder

backend_sim = Aer.get_backend("qasm_simulator")

d = 3
dx = 2
d = 3
dx = 2
code = RepetitionCodeCircuit(d, dx=dx, xbasis=False, resets=True)


for logical, qc in code.circuit.items():
    job = backend_sim.run(qc)
    counts = job.result().get_counts()
    output = list(counts.keys())[0]
    print("Output for a stored", logical, "is", output)
