# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2022
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Module for Tile Factory"""
from typing import List, Optional

import numpy as np
from qiskit_qec.geometry.model.edge import Edge
from qiskit_qec.geometry.model.face import Face
from qiskit_qec.geometry.model.shell import Shell
from qiskit_qec.geometry.model.vertex import Vertex
from qiskit_qec.geometry.model.wireframe import WireFrame
from qiskit_qec.geometry.model.qubit_count import QubitCount
from qiskit_qec.geometry.model.qubit_data import QubitData


class TileFactory:
    """Base class for all geometric tiles"""

    def __new__(
        cls,
        *,
        origin: np.ndarray,
        wf_coordinates: List,
        wf_q_indices: List,
        wf_loop_indicator: List[bool],
        faces_wf_components: List[List],
        num_qubits: int,
        face_colors: List,
        qubit_count: Optional[QubitCount] = None,
        qubit_data: Optional[QubitData] = None,
        operators: Optional[List] = None,
        **kwargs,
    ) -> Shell:
        """Create tile"""

        # Qubits

        if qubit_count is not None and qubit_data is not None:
            qubits = [qubit_count.new_qubit() for i in range(num_qubits)]

        # Create Shell
        wf_id_to_index = {}

        faces = []

        # print(f"faces_wf_components={faces_wf_components}")
        for wf_list in faces_wf_components:
            # print(f"wf_list={wf_list}")
            wfs = []
            for wf_index in wf_list:
                # print(f"wf_index={wf_index}")
                # print(f"operators[wf_index]={operators[wf_index]}")
                if operators[wf_index] is not None:
                    # print(f"Making wf for index = {wf_index}")
                    wf = cls._make_wireframe(
                        vertices=wf_coordinates[wf_index],
                        origin=origin,
                        loop_indicator=wf_loop_indicator[wf_index],
                    )
                    wfs.append(wf)
                    wf_id_to_index[wf.id] = wf_index
            if len(wfs) != 0:
                faces.append(Face(wfs))

        shell = Shell(faces)

        # Set the date for the shell
        # Vertex Data
        for face in shell.faces:
            for wf in face.wireframes:
                q_indices = wf_q_indices[wf_id_to_index[wf.id]]
                op_list = operators[wf_id_to_index[wf.id]]
                for index, vertex in enumerate(wf.vertices):
                    # Assign a qubit id to vertices
                    qubit_data.qubit[vertex.id] = qubits[q_indices[index]]
                    # Update the qubit count for the qubit
                    qubit_count.increment_qubit(qubits[q_indices[index]])
                    # Assign a Pauli operator(s) to vertices
                    qubit_data.operator[vertex.id] = [op[index] for op in op_list]

        # TODO: Add code to add extra vertex data

        # Edge Data
        # TODO:  Add code to add edge data

        # Wireframe Data
        prefix = "wf_"
        object_arrays = {key: array for key, array in kwargs.items() if key.startswith(prefix)}
        for name in object_arrays.keys():
            qubit_data.add_data_array(data_array={}, name=name)
        for wf in shell.wireframes:
            for name, array in object_arrays.items():
                qubit_data.data_arrays[name][wf.id] = array[wf_id_to_index[wf.id]]

        # Face Data

        for index, face in enumerate(shell.faces):
            qubit_data.face_colors[face.id] = face_colors[index]

        # TODO: Add code to add extra face data (if exists)

        # Shell Data
        # TODO:  Add code to add shell data

        return shell

    @staticmethod
    def _make_wireframe(vertices, origin, loop_indicator):
        # Vertices
        v = [Vertex(np.asarray(vertex) + origin) for vertex in vertices]

        # Edges
        num_vertices = len(vertices)
        if loop_indicator:
            edges = [Edge([v[i], v[(i + 1) % num_vertices]]) for i in range(num_vertices)]
        else:
            edges = [Edge([v[i], v[(i + 1)]]) for i in range(num_vertices - 1)]

        # Wireframe
        return WireFrame(edges)
