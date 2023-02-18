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

"""SquareDiamnondTile class."""

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.model.shell import Shell
from qiskit_qec.geometry.tiles.tile import Tile
from qiskit_qec.geometry.tiles.tilefactory import TileFactory
from qiskit_qec.operators.pauli_list import PauliList


# pylint: disable=anomalous-backslash-in-string)
class SquareDiamondTile(Tile):
    """Square Diamond Tile
    (Square or non-Rotated orientation)

     q0    q1  q1   q2
     v0    v1  v0   v1
       o--o     o--o
       |0/       \1|
       |/         \|
       o           o
     v2             v2
     q3             q4
            .(0,0)
     q3            q4
     v0            v0
       o           o
       |\         /|
       |2\       /3|
       o--o     o--o
     v1    v2 v1    v2
     q5    q6 q6    q7

       o--o--o
       |0/ \1|
       |/   \|
       o     o
       |\   /|
       |2\ /3|
       o--o--o

    """

    wf_operator_dict = {
        "pXZXZ": [PauliList(["XXX"]), PauliList(["ZZZ"]), PauliList(["ZZZ"]), PauliList(["XXX"])],
        "pZXZX": [PauliList(["ZZZ"]), PauliList(["XXX"]), PauliList(["XXX"]), PauliList(["ZZZ"])],
        "cXXXX": [PauliList(["XXX"]), PauliList(["XXX"]), PauliList(["XXX"]), PauliList(["XXX"])],
        "cZZZZ": [PauliList(["ZZZ"]), PauliList(["ZZZ"]), PauliList(["ZZZ"]), PauliList(["ZZZ"])],
    }

    # Descriptions of wireframes
    # qubit indices for each wireframe
    wf_q_indices = [[0, 1, 3], [1, 2, 4], [3, 5, 6], [4, 6, 7]]

    # coordinates for each wireframe vertex in path list form to enable the
    # creation of the associate edges
    wf_coordinates = [
        [[-1, 1], [0, 1], [-1, 0]],
        [[0, 1], [1, 1], [1, 0]],
        [[-1, 0], [-1, -1], [0, -1]],
        [[1, 0], [0, -1], [1, -1]],
    ]

    # If the wf's are closed loops or not
    wf_loop_indicator = [True, True, True, True]

    # Wireframe orientations
    wf_orientation = [[-1, 1], [1, 1], [-1, -1], [1, -1]]

    # How wireframes components combine into faces/operators
    faces_wf_components = [[0], [1], [2], [3]]

    # Face colors (wf's inherit colors from faces)
    face_colors = ["yellowgreen", "tomato", "tomato", "yellowgreen"]

    num_faces = 4

    size = np.array([2, 2])

    num_qubits = 8

    u_vec = np.array([0, 2])
    v_vec = np.array([2, 0])

    def __new__(
        cls,
        origin: np.array,
        qubit_count=None,
        qubit_data=None,
        operators=None,
        optype="pXZXZ",
    ) -> Shell:
        """Square Diamond Tile

         q0    q1  q1   q2
         v0    v1  v0   v1
           o--o     o--o
           |0/       \1|
           |/         \|
           o           o
         v2             v2
         q3             q4
                .(0,0)
         q3            q4
         v0            v0
           o           o
           |\         /|
           |2\       /3|
           o--o     o--o
         v1    v2 v1    v2
         q5    q6 q6    q7

           o--o--o
           |0/ \1|
           |/   \|
           o     o
           |\   /|
           |2\ /3|
           o--o--o


        Face colors for faces [0,1,2] are ["yellowgreen","tomato","steelblue"]

        Preformatted operators are stored in HexagonTile.op_dict. Keys for op_dict are of the
        form [p|c|d]PPPP... where p = pattern and c = copy, d=double and P is a Pauli opertor X, Z or Y.

        "pXZXZ" -> face #0 operator is PauliList(["XXX"]),
                   face #1 operator is PauliList(["ZZZ"]),
                   face #2 operator is PauliList(["ZZZ"]), and
                   face #3 operator is PauliList(["XXX"])]

        Available precomputed operator layouts are:

        "pXZXZ", "pZXZX", "cXXXX", "cZZZZ"

        The operator variable may be used to define the operators specifically. The operator must be
        a list of PauliList objects where each PauliList describes the opertors to be built for the
        faces as indexed above 0,1,2,3, ... If the PauliList contains k Paulis then k operators will
        be created for the given face.

        Args:
            origin (np.array): Coordinates of origin of tile (shell)
            qubit_count: Qubit counter. Defaults to None.
            qubit_data: Qubit data. Defaults to None.
            operators: Operators for tile faces. Defaults to None.
            optype (optional): Which of the listed opertor mapppings to used. Defaults to "pXZXZ".

        A appropriately scalled lattice basis for tiling with the HexagonTile can be accessed
        from [HexagonTile.u_vec, HexagonTile.v_vec]

        Raises:
            QiskitError: Unsupported operator type

        Returns:
            Shell: Returns a Square Diamond tile (shell) with provided origin
        """

        if operators is None:
            try:
                operators = SquareDiamondTile.wf_operator_dict[optype]
            except KeyError as key_error:
                raise QiskitError(f"Unsupported operator type: {optype}") from key_error

        return TileFactory(
            origin=origin,
            wf_coordinates=SquareDiamondTile.wf_coordinates,
            wf_q_indices=SquareDiamondTile.wf_q_indices,
            wf_loop_indicator=SquareDiamondTile.wf_loop_indicator,
            faces_wf_components=SquareDiamondTile.faces_wf_components,
            num_qubits=SquareDiamondTile.num_qubits,
            face_colors=SquareDiamondTile.face_colors,
            qubit_count=qubit_count,
            qubit_data=qubit_data,
            operators=operators,
            wf_orientation=SquareDiamondTile.wf_orientation,
        )
