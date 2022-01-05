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

import numpy as np

from qiskit.exceptions import QiskitError

from qiskit_qec.operators.pauli_list import PauliList

from qiskit_qec.geometry.model.vertex import Vertex
from qiskit_qec.geometry.model.edge import Edge
from qiskit_qec.geometry.model.wireframe import WireFrame
from qiskit_qec.geometry.model.face import Face
from qiskit_qec.geometry.model.shell import Shell
from qiskit_qec.geometry.model.tile import Tile
from qiskit_qec.geometry.lattice import Lattice

class Tiling:
    def __new__(
        cls,
        *,
        tile: Tile,
        lattice: Lattice,
        qubit_count=None, 
        qubit_data=None,
        epsilon=0.01
    ):
        if lattice.points is None:
            raise QiskitError("Infinite lattices not yet supported")

        point = lattice.points[0]
        M = tile(center=point, qubit_count=qubit_count, qubit_data=qubit_data)

        for point in lattice.points[1:]:
            T = tile(center=point, qubit_count=qubit_count, qubit_data=qubit_data)
            i_vertices_num = len(M.vertices)
            M.union(T)
        
            for i_vertex in M.vertices[:i_vertices_num]:
                for u_vertex in M.vertices[i_vertices_num:]:
                    if cls.distance(i_vertex.pos, u_vertex.pos) < epsilon:
                        #print("\n")
                        #print(f"u_vertex.pos={u_vertex.pos} i_vertex.pos={i_vertex.pos}")

                        # This is a very simple form of snapping
                        # Works only due to the orientations all being
                        # the same for each tile.
                        u_vertex.pos = i_vertex.pos

                        # Align the qubits: set the two vertices to use the same qubit

                        u_vertex_qubit_id = qubit_data.qubit[u_vertex.id]
                        old_u_vertex_qubit_id = u_vertex_qubit_id
                        i_vertex_qubit_id = qubit_data.qubit[i_vertex.id]
                        
                        #print(f"u_vertex_qubit_id = {u_vertex_qubit_id}")
                        #print(f"i_vertex_qubit_id= {i_vertex_qubit_id}")

                        qubit_data.qubit[u_vertex.id] = i_vertex_qubit_id
                        u_vertex_qubit_id = i_vertex_qubit_id

                        #print("After merge")
                        #print(f"u_vertex_qubit_id = {u_vertex_qubit_id}")
                        #print(f"i_vertex_qubit_id = {i_vertex_qubit_id}")

                        #print("qubit counters")
                        #print(f"Old: {qubit_count.qubits_count[old_u_vertex_qubit_id]}")
                        #print(f"New: {qubit_count.qubits_count[i_vertex_qubit_id]}")

                        # Adjust the qubit_counters
                        qubit_count.decrement_qubit(old_u_vertex_qubit_id)
                        qubit_count.increment_qubit(i_vertex_qubit_id)
                        
                        #print("After update")
                        #print(f"Old: {qubit_count.qubits_count[old_u_vertex_qubit_id]}")
                        #print(f"New: {qubit_count.qubits_count[i_vertex_qubit_id]}")
        
        return M

            
    def distance(point1, point2):
        return np.linalg.norm(point1-point2)
    
