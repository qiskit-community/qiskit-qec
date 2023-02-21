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
"""Module for Tile Generation"""
from typing import Type

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.lattice import Lattice
from qiskit_qec.geometry.model.qubit_count import QubitCount
from qiskit_qec.geometry.model.qubit_data import QubitData
from qiskit_qec.geometry.tiles.tile import Tile


class Tiling:
    """Create a tiling[Shell] from a tile and lattice

    Args:
        tile_type (Type[Tile]): Type of tile to create. Ex. RSSC
        lattice (Lattice): Lattice on which to place tiles
        qubit_count (QubitCount, optional): Keeps track of qubits
        qubit_data (QubitData, optional): Keeps track of data associated code associated to tiling
        epsilon (float, optional): Minimum distance between two vertices to be considered separate
        relative to being assigned different qubits Defaults to 0.01.
    """

    def __new__(
        cls,
        *,
        tile_type: Type[Tile],
        tile_optype: str,
        lattice: Lattice,
        qubit_count: QubitCount = None,
        qubit_data: QubitData = None,
        epsilon: float = 0.01,
    ):
        if lattice.points is None:
            raise QiskitError("Infinite lattices not yet supported")

        point = lattice.points[0]
        shell = tile_type(
            origin=point, optype=tile_optype, qubit_count=qubit_count, qubit_data=qubit_data
        )

        for point in lattice.points[1:]:
            tile = tile_type(
                origin=point, optype=tile_optype, qubit_count=qubit_count, qubit_data=qubit_data
            )
            i_vertices_num = len(shell.vertices)
            shell.union(tile)

            for i_vertex in shell.vertices[:i_vertices_num]:
                for u_vertex in shell.vertices[i_vertices_num:]:
                    if cls.distance(i_vertex.pos, u_vertex.pos) < epsilon:
                        # This is a very simple form of snapping
                        # Works only due to the orientations all being
                        # the same for each tile.
                        u_vertex.pos = i_vertex.pos

                        # Align the qubits: set the two vertices to use the same qubit

                        u_vertex_qubit_id = qubit_data.qubit[u_vertex.id]
                        old_u_vertex_qubit_id = u_vertex_qubit_id
                        i_vertex_qubit_id = qubit_data.qubit[i_vertex.id]

                        qubit_data.qubit[u_vertex.id] = i_vertex_qubit_id
                        u_vertex_qubit_id = i_vertex_qubit_id

                        # Adjust the qubit_counters
                        qubit_count.decrement_qubit(old_u_vertex_qubit_id)
                        qubit_count.increment_qubit(i_vertex_qubit_id)

        return shell

    @classmethod
    def distance(cls, point1, point2):
        """Distance between two points."""
        return np.linalg.norm(point1 - point2)
