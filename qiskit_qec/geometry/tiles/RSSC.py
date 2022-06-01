# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""RSSC."""

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.model.edge import Edge
from qiskit_qec.geometry.model.face import Face
from qiskit_qec.geometry.model.shell import Shell
from qiskit_qec.geometry.model.vertex import Vertex
from qiskit_qec.geometry.model.wireframe import WireFrame
from qiskit_qec.geometry.tiles.tile import Tile
from qiskit_qec.operators.pauli_list import PauliList


# pylint: disable=anomalous-backslash-in-string)
class RSSC(Tile):
    """
        (Square or non-Rotated orientation)

         q0    q1  q1   q2
         v0    v1  v0   v1
           o--o     o--o
           |0/       \1|
           o/         \o
         v2             v2
         q3             q4

         q3            q4
         v0            v0
           o           o
           |\         /|
           |2\       /3|
           o--o     o--o
         v1    v2 v1    v2
         q5    q6 q6    q7

    Orientations:

    Face: 0 [-1,1]  -> after 45 deg ac rot [-1,0] left
    Face: 1 [1,1]
    Face: 2 [-1,-1]
    Face: 4 [1,-1]

    """

    op_dict = {
        "XZXZ": PauliList(["XXX", "ZZZ", "ZZZ", "XXX"]),
        "ZXZX": PauliList(["ZZZ", "XXX", "XXX", "ZZZ"]),
        "XXXX": PauliList(["XXX", "XXX", "XXX", "XXX"]),
        "ZZZZ": PauliList(["ZZZ", "ZZZ", "ZZZ", "ZZZ"]),
    }

    q_indices = [[0, 1, 3], [1, 2, 4], [3, 5, 6], [4, 6, 7]]

    op_indices = [
        [[-1, 1], [0, 1], [-1, 0]],
        [[0, 1], [1, 1], [1, 0]],
        [[-1, 0], [-1, -1], [0, -1]],
        [[1, 0], [0, -1], [1, -1]],
    ]

    size = np.array([2, 2])

    def __new__(
        cls,
        center: np.array,
        qubit_count=None,
        qubit_data=None,
        operators=None,
        optype="XZXZ",
    ) -> None:

        if operators is None:
            try:
                operators = RSSC.op_dict[optype]
            except KeyError as key_error:
                raise QiskitError(f"Unsupported operator type: {optype}") from key_error

        # Qubits

        if qubit_count is not None and qubit_data is not None:
            qubits = [qubit_count.new_qubit() for i in range(8)]

        faces = []

        for i in range(4):
            faces.append(cls._make_operator(vertices=RSSC.op_indices[i], center=center))

            # This should be made much better to allow general data etc
            if qubit_count is not None and qubit_data is not None:
                cls._assign_qubit_data(
                    faces[i].vertices,
                    RSSC.q_indices[i],
                    qubit_data,
                    qubit_count,
                    qubits,
                    operators[i],
                )

        # Shell

        shell = Shell(faces)

        return shell

    @staticmethod
    def _make_operator(vertices, center):
        # Vertices
        v = [Vertex(np.asarray(vertex) + center) for vertex in vertices]

        # Edges
        edges = [Edge([v[i], v[(i + 1) % 3]]) for i in range(len(vertices))]

        # Wireframe
        wf = WireFrame(edges)

        return Face(wf)

    @staticmethod
    def _assign_qubit_data(vertices, q_indices, data, count, qubits, op):

        for index, vertex in enumerate(vertices):
            # Assign a qubit id to vertices
            data.qubit[vertex.id] = qubits[q_indices[index]]
            # Update the qubit count for the qubit
            count.increment_qubit(qubits[q_indices[index]])
            # Assign a Pauli operator to vertices
            data.operator[vertex.id] = op[index]
