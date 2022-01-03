from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit import execute, Aer
from qiskit.providers.aer import noise

from sympy import Symbol, Poly

import networkx as nx

from pymatching import Matching

from math import log
from copy import copy

import numpy as np

p = 0.01


def make_noise_model(p1, p2, q):
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


model = make_noise_model(p, p, p)


xgauges = [[0, 3], [1, 4], [2, 5], [3, 6], [4, 7], [5, 8]]
zgauges = [[0, 1], [1, 2, 4, 5], [3, 4, 6, 7], [7, 8]]
xstab = [[0, 1, 3, 4], [2, 5], [3, 6], [4, 5, 7, 8]]
xgaugeproduct = [[0, 1], [2], [3], [4, 5]]
zstab = [[0, 1, 3, 4, 6, 7], [1, 2, 4, 5, 7, 8]]
zgaugeproduct = [[0, 2], [1, 3]]
xboundary = [[0], [1], [2], [6], [7], [8]]
zboundary = [[0], [3], [6], [2], [5], [8]]

zstabmat = np.array([[1, 1, 0, 1, 1, 0, 1, 1, 0], [0, 1, 1, 0, 1, 1, 0, 1, 1]])
xstabmat = np.array(
    [
        [1, 1, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 1, 0, 1, 1],
    ]
)
allonesmat = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1]])


q = QuantumRegister(23, "q0")
c = ClassicalRegister(51, "c0")
circ = QuantumCircuit(q, c, name="hhc3z")
circ.reset(q[0])
circ.reset(q[1])
circ.reset(q[2])
circ.reset(q[3])
circ.reset(q[4])
circ.reset(q[5])
circ.reset(q[6])
circ.reset(q[7])
circ.reset(q[8])
circ.barrier(q)
circ.reset(q[10])
circ.h(q[10])
circ.reset(q[11])
circ.h(q[11])
circ.reset(q[12])
circ.h(q[12])
circ.reset(q[13])
circ.h(q[13])
circ.reset(q[14])
circ.h(q[14])
circ.reset(q[9])
circ.h(q[9])
circ.barrier(q)
circ.cx(q[10], q[1])
circ.cx(q[11], q[2])
circ.cx(q[12], q[3])
circ.cx(q[13], q[4])
circ.cx(q[14], q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.cx(q[9], q[0])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.cx(q[10], q[4])
circ.cx(q[11], q[5])
circ.cx(q[12], q[6])
circ.cx(q[13], q[7])
circ.cx(q[14], q[8])
circ.i(q[2])
circ.cx(q[9], q[3])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.h(q[10])
circ.h(q[11])
circ.h(q[12])
circ.h(q[13])
circ.h(q[14])
circ.i(q[2])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.h(q[9])
circ.measure(q[9], c[0])
circ.measure(q[10], c[1])
circ.measure(q[11], c[2])
circ.measure(q[12], c[3])
circ.measure(q[13], c[4])
circ.measure(q[14], c[5])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.reset(q[11])
circ.h(q[11])
circ.reset(q[13])
circ.h(q[13])
circ.reset(q[16])
circ.reset(q[17])
circ.h(q[17])
circ.reset(q[18])
circ.reset(q[19])
circ.i(q[2])
circ.reset(q[21])
circ.reset(q[22])
circ.h(q[22])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.reset(q[10])
circ.h(q[10])
circ.cx(q[11], q[18])
circ.reset(q[12])
circ.h(q[12])
circ.cx(q[13], q[19])
circ.reset(q[15])
circ.h(q[15])
circ.cx(q[17], q[16])
circ.i(q[2])
circ.reset(q[20])
circ.h(q[20])
circ.cx(q[22], q[21])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.cx(q[10], q[18])
circ.cx(q[12], q[19])
circ.cx(q[15], q[16])
circ.i(q[17])
circ.cx(q[2], q[11])
circ.cx(q[20], q[21])
circ.i(q[3])
circ.cx(q[4], q[13])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.cx(q[8], q[22])
circ.barrier(q)
circ.cx(q[0], q[15])
circ.cx(q[1], q[17])
circ.i(q[2])
circ.i(q[20])
circ.i(q[22])
circ.i(q[3])
circ.cx(q[4], q[10])
circ.cx(q[5], q[11])
circ.cx(q[6], q[12])
circ.cx(q[7], q[13])
circ.i(q[8])
circ.barrier(q)
circ.i(q[0])
circ.cx(q[1], q[10])
circ.cx(q[11], q[18])
circ.cx(q[13], q[19])
circ.i(q[15])
circ.cx(q[17], q[16])
circ.i(q[2])
circ.cx(q[22], q[21])
circ.cx(q[3], q[12])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.cx(q[7], q[20])
circ.i(q[8])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.cx(q[10], q[18])
circ.h(q[11])
circ.cx(q[12], q[19])
circ.h(q[13])
circ.cx(q[15], q[16])
circ.h(q[17])
circ.i(q[2])
circ.cx(q[20], q[21])
circ.h(q[22])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.measure(q[17], c[14])
circ.measure(q[11], c[15])
circ.measure(q[13], c[16])
circ.measure(q[22], c[17])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.h(q[10])
circ.measure(q[10], c[11])
circ.h(q[12])
circ.measure(q[12], c[12])
circ.h(q[15])
circ.measure(q[15], c[10])
circ.measure(q[16], c[6])
circ.measure(q[18], c[7])
circ.measure(q[19], c[8])
circ.i(q[2])
circ.h(q[20])
circ.measure(q[20], c[13])
circ.measure(q[21], c[9])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.barrier(q)
circ.reset(q[10])
circ.h(q[10])
circ.reset(q[11])
circ.h(q[11])
circ.reset(q[12])
circ.h(q[12])
circ.reset(q[13])
circ.h(q[13])
circ.reset(q[14])
circ.h(q[14])
circ.reset(q[9])
circ.h(q[9])
circ.barrier(q)
circ.cx(q[10], q[1])
circ.cx(q[11], q[2])
circ.cx(q[12], q[3])
circ.cx(q[13], q[4])
circ.cx(q[14], q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.cx(q[9], q[0])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.cx(q[10], q[4])
circ.cx(q[11], q[5])
circ.cx(q[12], q[6])
circ.cx(q[13], q[7])
circ.cx(q[14], q[8])
circ.i(q[2])
circ.cx(q[9], q[3])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.h(q[10])
circ.h(q[11])
circ.h(q[12])
circ.h(q[13])
circ.h(q[14])
circ.i(q[2])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.h(q[9])
circ.measure(q[9], c[18])
circ.measure(q[10], c[19])
circ.measure(q[11], c[20])
circ.measure(q[12], c[21])
circ.measure(q[13], c[22])
circ.measure(q[14], c[23])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.reset(q[11])
circ.h(q[11])
circ.reset(q[13])
circ.h(q[13])
circ.reset(q[16])
circ.reset(q[17])
circ.h(q[17])
circ.reset(q[18])
circ.reset(q[19])
circ.i(q[2])
circ.reset(q[21])
circ.reset(q[22])
circ.h(q[22])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.reset(q[10])
circ.h(q[10])
circ.cx(q[11], q[18])
circ.reset(q[12])
circ.h(q[12])
circ.cx(q[13], q[19])
circ.reset(q[15])
circ.h(q[15])
circ.cx(q[17], q[16])
circ.i(q[2])
circ.reset(q[20])
circ.h(q[20])
circ.cx(q[22], q[21])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.cx(q[10], q[18])
circ.cx(q[12], q[19])
circ.cx(q[15], q[16])
circ.i(q[17])
circ.cx(q[2], q[11])
circ.cx(q[20], q[21])
circ.i(q[3])
circ.cx(q[4], q[13])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.cx(q[8], q[22])
circ.barrier(q)
circ.cx(q[0], q[15])
circ.cx(q[1], q[17])
circ.i(q[2])
circ.i(q[20])
circ.i(q[22])
circ.i(q[3])
circ.cx(q[4], q[10])
circ.cx(q[5], q[11])
circ.cx(q[6], q[12])
circ.cx(q[7], q[13])
circ.i(q[8])
circ.barrier(q)
circ.i(q[0])
circ.cx(q[1], q[10])
circ.cx(q[11], q[18])
circ.cx(q[13], q[19])
circ.i(q[15])
circ.cx(q[17], q[16])
circ.i(q[2])
circ.cx(q[22], q[21])
circ.cx(q[3], q[12])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.cx(q[7], q[20])
circ.i(q[8])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.cx(q[10], q[18])
circ.h(q[11])
circ.cx(q[12], q[19])
circ.h(q[13])
circ.cx(q[15], q[16])
circ.h(q[17])
circ.i(q[2])
circ.cx(q[20], q[21])
circ.h(q[22])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.measure(q[17], c[32])
circ.measure(q[11], c[33])
circ.measure(q[13], c[34])
circ.measure(q[22], c[35])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.h(q[10])
circ.measure(q[10], c[29])
circ.h(q[12])
circ.measure(q[12], c[30])
circ.h(q[15])
circ.measure(q[15], c[28])
circ.measure(q[16], c[24])
circ.measure(q[18], c[25])
circ.measure(q[19], c[26])
circ.i(q[2])
circ.h(q[20])
circ.measure(q[20], c[31])
circ.measure(q[21], c[27])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.barrier(q)
circ.reset(q[10])
circ.h(q[10])
circ.reset(q[11])
circ.h(q[11])
circ.reset(q[12])
circ.h(q[12])
circ.reset(q[13])
circ.h(q[13])
circ.reset(q[14])
circ.h(q[14])
circ.reset(q[9])
circ.h(q[9])
circ.barrier(q)
circ.cx(q[10], q[1])
circ.cx(q[11], q[2])
circ.cx(q[12], q[3])
circ.cx(q[13], q[4])
circ.cx(q[14], q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.cx(q[9], q[0])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.cx(q[10], q[4])
circ.cx(q[11], q[5])
circ.cx(q[12], q[6])
circ.cx(q[13], q[7])
circ.cx(q[14], q[8])
circ.i(q[2])
circ.cx(q[9], q[3])
circ.barrier(q)
circ.i(q[0])
circ.i(q[1])
circ.h(q[10])
circ.h(q[11])
circ.h(q[12])
circ.h(q[13])
circ.h(q[14])
circ.i(q[2])
circ.i(q[3])
circ.i(q[4])
circ.i(q[5])
circ.i(q[6])
circ.i(q[7])
circ.i(q[8])
circ.h(q[9])
circ.measure(q[9], c[36])
circ.measure(q[10], c[37])
circ.measure(q[11], c[38])
circ.measure(q[12], c[39])
circ.measure(q[13], c[40])
circ.measure(q[14], c[41])
circ.barrier(q)
circ.measure(q[0], c[42])
circ.measure(q[1], c[43])
circ.measure(q[2], c[44])
circ.measure(q[3], c[45])
circ.measure(q[4], c[46])
circ.measure(q[5], c[47])
circ.measure(q[6], c[48])
circ.measure(q[7], c[49])
circ.measure(q[8], c[50])

print(circ.qasm())


cx = Symbol("cx")
id = Symbol("id")
reset = Symbol("reset")
h = Symbol("h")
measure = Symbol("measure")
symbols = {"cx": cx, "id": id, "reset": reset, "h": h, "measure": measure}
qubit_map = {1: 0, 0: 1, 3: 2, 6: 3, 2: 4, 5: 5, 8: 6}
reverse_qubit_map = {0: 1, 1: 0, 2: 3, 3: 6, 4: 2, 5: 5, 6: 8}
stabilizer_to_vertex_map = {
    (0, (0, 1, 3, 4, 6, 7)): 0,
    (0, (1, 2, 4, 5, 7, 8)): 1,
    (0, (0,)): 2,
    (0, (3,)): 3,
    (0, (6,)): 4,
    (0, (2,)): 5,
    (0, (5,)): 6,
    (0, (8,)): 7,
    (1, (0, 1, 3, 4, 6, 7)): 8,
    (1, (1, 2, 4, 5, 7, 8)): 9,
    (1, (0,)): 10,
    (1, (3,)): 11,
    (1, (6,)): 12,
    (1, (2,)): 13,
    (1, (5,)): 14,
    (1, (8,)): 15,
    (2, (0, 1, 3, 4, 6, 7)): 16,
    (2, (1, 2, 4, 5, 7, 8)): 17,
    (2, (0,)): 18,
    (2, (3,)): 19,
    (2, (6,)): 20,
    (2, (2,)): 21,
    (2, (5,)): 22,
    (2, (8,)): 23,
}
g = nx.Graph()
g.add_edge(
    0,
    1,
    qubit_id=0,
    weight=1,
    weight_poly=Poly(9 / 5 * cx + 26 / 3 * id + 3 * reset, cx, id, reset, domain="QQ"),
)
g.add_edge(
    0,
    2,
    qubit_id=1,
    weight=1,
    weight_poly=Poly(9 / 5 * cx + 10 * id + 3 * reset, cx, id, reset, domain="QQ"),
)
g.add_edge(0, 3, qubit_id=2, weight=1000000)
g.add_edge(0, 4, qubit_id=3, weight=1000000)
g.add_edge(
    0,
    8,
    qubit_id=-1,
    weight=1,
    weight_poly=Poly(
        16 / 5 * cx + 4 / 3 * id + 2 * measure + 2 * reset,
        cx,
        id,
        measure,
        reset,
        domain="QQ",
    ),
)
g.add_edge(
    1,
    5,
    qubit_id=4,
    weight=1,
    weight_poly=Poly(9 / 5 * cx + 8 * id + 3 * reset, cx, id, reset, domain="QQ"),
)
g.add_edge(1, 6, qubit_id=5, weight=1000000)
g.add_edge(1, 7, qubit_id=6, weight=1000000)
g.add_edge(1, 8, qubit_id=0, weight=1, weight_poly=Poly(cx, cx, domain="ZZ"))
g.add_edge(
    1,
    9,
    qubit_id=-1,
    weight=1,
    weight_poly=Poly(
        16 / 5 * cx + 4 / 3 * id + 2 * measure + 2 * reset,
        cx,
        id,
        measure,
        reset,
        domain="QQ",
    ),
)
g.add_edge(2, 3, qubit_id=-1, weight=0)
g.add_edge(2, 10, qubit_id=-1, weight=0)
g.add_edge(3, 4, qubit_id=-1, weight=0)
g.add_edge(4, 5, qubit_id=-1, weight=0)
g.add_edge(5, 6, qubit_id=-1, weight=0)
g.add_edge(6, 7, qubit_id=-1, weight=0)
g.add_edge(
    8,
    9,
    qubit_id=0,
    weight=1,
    weight_poly=Poly(11 / 5 * cx + 40 / 3 * id, cx, id, domain="QQ"),
)
g.add_edge(
    8,
    10,
    qubit_id=1,
    weight=1,
    weight_poly=Poly(11 / 5 * cx + 46 / 3 * id, cx, id, domain="QQ"),
)
g.add_edge(8, 11, qubit_id=2, weight=1000000)
g.add_edge(8, 12, qubit_id=3, weight=1000000)
g.add_edge(
    8,
    16,
    qubit_id=-1,
    weight=1,
    weight_poly=Poly(
        16 / 5 * cx + 4 / 3 * id + 2 * measure + 2 * reset,
        cx,
        id,
        measure,
        reset,
        domain="QQ",
    ),
)
g.add_edge(
    9,
    13,
    qubit_id=4,
    weight=1,
    weight_poly=Poly(11 / 5 * cx + 46 / 3 * id, cx, id, domain="QQ"),
)
g.add_edge(9, 14, qubit_id=5, weight=1000000)
g.add_edge(9, 15, qubit_id=6, weight=1000000)
g.add_edge(9, 16, qubit_id=0, weight=1, weight_poly=Poly(cx, cx, domain="ZZ"))
g.add_edge(
    9,
    17,
    qubit_id=-1,
    weight=1,
    weight_poly=Poly(
        16 / 5 * cx + 4 / 3 * id + 2 * measure + 2 * reset,
        cx,
        id,
        measure,
        reset,
        domain="QQ",
    ),
)
g.add_edge(10, 11, qubit_id=-1, weight=0)
g.add_edge(10, 18, qubit_id=-1, weight=0)
g.add_edge(11, 12, qubit_id=-1, weight=0)
g.add_edge(12, 13, qubit_id=-1, weight=0)
g.add_edge(13, 14, qubit_id=-1, weight=0)
g.add_edge(14, 15, qubit_id=-1, weight=0)
g.add_edge(
    16,
    17,
    qubit_id=0,
    weight=1,
    weight_poly=Poly(8 / 5 * cx + 8 * id + 3 * measure, cx, id, measure, domain="QQ"),
)
g.add_edge(
    16,
    18,
    qubit_id=1,
    weight=1,
    weight_poly=Poly(
        8 / 5 * cx + 26 / 3 * id + 3 * measure, cx, id, measure, domain="QQ"
    ),
)
g.add_edge(16, 19, qubit_id=2, weight=1000000)
g.add_edge(16, 20, qubit_id=3, weight=1000000)
g.add_edge(
    17,
    21,
    qubit_id=4,
    weight=1,
    weight_poly=Poly(
        8 / 5 * cx + 32 / 3 * id + 3 * measure, cx, id, measure, domain="QQ"
    ),
)
g.add_edge(17, 22, qubit_id=5, weight=1000000)
g.add_edge(17, 23, qubit_id=6, weight=1000000)
g.add_edge(18, 19, qubit_id=-1, weight=0)
g.add_edge(19, 20, qubit_id=-1, weight=0)
g.add_edge(20, 21, qubit_id=-1, weight=0)
g.add_edge(21, 22, qubit_id=-1, weight=0)
g.add_edge(22, 23, qubit_id=-1, weight=0)
g.nodes[2]["is_boundary"] = True
g.nodes[3]["is_boundary"] = True
g.nodes[4]["is_boundary"] = True
g.nodes[5]["is_boundary"] = True
g.nodes[6]["is_boundary"] = True
g.nodes[7]["is_boundary"] = True
g.nodes[10]["is_boundary"] = True
g.nodes[11]["is_boundary"] = True
g.nodes[12]["is_boundary"] = True
g.nodes[13]["is_boundary"] = True
g.nodes[14]["is_boundary"] = True
g.nodes[15]["is_boundary"] = True
g.nodes[18]["is_boundary"] = True
g.nodes[19]["is_boundary"] = True
g.nodes[20]["is_boundary"] = True
g.nodes[21]["is_boundary"] = True
g.nodes[22]["is_boundary"] = True
g.nodes[23]["is_boundary"] = True

print(g.nodes(data=True))
print(g.edges(data=True))

seed = 1337
sim = Aer.get_backend("qasm_simulator")
ideal_result = execute(
    circ, sim, method="stabilizer", shots=100, optimization_level=0, seed_simulator=seed
).result()
ideal_counts = ideal_result.get_counts(circ)
print(ideal_counts)

noisy_result = execute(
    circ,
    sim,
    noise_model=model,
    basis_gates=model.basis_gates,
    method="stabilizer",
    shots=100,
    optimization_level=0,
    seed_simulator=seed,
).result()
noisy_counts = noisy_result.get_counts(circ)
print(noisy_counts)


def set_weights_and_create_matcher(p):
    param_names = ["cx", "h", "measure", "reset", "id"]
    param_values = [p, p, p, p, p]
    symbol_list = [symbols[s] for s in param_names]
    assignment = {k: v for (k, v) in zip(symbol_list, param_values)}
    for e in g.edges(data=True):
        if "weight_poly" in e[2]:
            restriction = {x: assignment[x] for x in e[2]["weight_poly"].gens}
            pflip = e[2]["weight_poly"].eval(restriction).evalf()
            e[2]["weight"] = log((1 - p) / p)
    matcher = Matching(g)
    return matcher


matcher = set_weights_and_create_matcher(p)
print(matcher)


def partition_z(outcome):
    outcome = list(map(int, outcome[::-1]))
    x_out = [outcome[0:6], outcome[18:24], outcome[36:42]]
    z_out = [outcome[6:10], outcome[24:28]]
    lflag_out = [outcome[10:14], outcome[28:32]]
    rflag_out = [outcome[14:18], outcome[32:36]]
    final_out = outcome[42:51]
    return z_out, final_out


def highlighted_vertices_z(z_outcome, final_outcome):
    final_gauges = []
    for g in zgauges:
        final_gauges.append(sum([final_outcome[i] for i in g]) % 2)
    gauge_outcomes = z_outcome
    gauge_outcomes.append(final_gauges)
    highlighted = []
    for i in range(3):
        for j in range(len(zstab)):
            prior = 0
            if i > 0:
                prior = sum([gauge_outcomes[i - 1][k] for k in zgaugeproduct[j]]) % 2
            current = sum([gauge_outcomes[i][k] for k in zgaugeproduct[j]]) % 2
            if current != prior:
                highlighted.append((i, tuple(zstab[j])))
    # If odd number, add vertex at boundary
    if len(highlighted) % 2 == 1:
        highlighted.append((0, tuple(zboundary[0])))
    return highlighted


def process_outcome_z(outcome, num, matcher):
    z_out, final_out = partition_z(outcome)
    highlighted = highlighted_vertices_z(z_out, final_out)
    syndrome = [0] * len(stabilizer_to_vertex_map)
    for stab in highlighted:
        syndrome[stabilizer_to_vertex_map[stab]] = 1
    correction = matcher.decode(syndrome)
    corrected_out = copy(final_out)
    for i in range(len(correction)):
        if correction[i] == 1:
            j = reverse_qubit_map[i]
            corrected_out[j] = (corrected_out[j] + 1) % 2
    zsynd = np.dot(zstabmat, np.array(corrected_out)) % 2
    if zsynd.any():
        raise Exception("non-trivial syndrome after QEC -- bug!")
    logical = np.dot(allonesmat, np.array(corrected_out)) % 2
    return logical.any()


failures = sum(
    map(lambda x: process_outcome_z(x[0], x[1], matcher), ideal_counts.items())
)

print(failures)

failures = sum(
    map(lambda x: process_outcome_z(x[0], x[1], matcher), noisy_counts.items())
)

print(failures)

pin = [0.0001, 0.0003, 0.0005, 0.0007, 0.0009]
pout = []
samples = [100000, 40000, 10000, 1000, 1000]
for i in range(5):
    p = pin[i]
    print(p, flush=True)
    model = make_noise_model(p, p, p)
    matcher = set_weights_and_create_matcher(p)
    noisy_result = execute(
        circ,
        sim,
        noise_model=model,
        basis_gates=model.basis_gates,
        method="stabilizer",
        shots=samples[i],
        optimization_level=0,
        seed_simulator=seed,
    ).result()
    noisy_counts = noisy_result.get_counts(circ)
    failures = sum(
        map(lambda x: process_outcome_z(x[0], x[1], matcher), noisy_counts.items())
    )
    print(failures, flush=True)
    pout.append(float(failures) / float(samples[i]))

print(pin)
print(pout)
