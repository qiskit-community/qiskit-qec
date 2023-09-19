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

"""CheckerBoardTile class."""

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.model.shell import Shell
from qiskit_qec.geometry.tiles.tile import Tile
from qiskit_qec.geometry.tiles.tilefactory import TileFactory
from qiskit_qec.operators.pauli_list import PauliList


# pylint: disable=anomalous-backslash-in-string)
class CheckerBoardTile(Tile):
    """Checker Board Tile

    q0       q1  q1       q2
    v0       v1  v0       v1
      o-----o      o-----o
      |  0  |      |  1  |
      o-----o      o-----o
    v2       v3  v2       v3
    q3       q4  q4       q5

                .(0,0)

    q3       q4  q4       q5
    v0       v1  v0       v1
      o-----o      o-----o
      |  2  |      |  3  |
      o-----o      o-----o
    v2       v3  v2       v3
    q6       q7  q7       q8


      o-----o-----o
      |  0  |  1  |
      o-----o-----o
      |  2  |  3  |
      o-----o-----o

    """

    wf_operator_dict = {
        "pXZXZ": [
            PauliList(["XXXX"]),
            PauliList(["ZZZZ"]),
            PauliList(["ZZZZ"]),
            PauliList(["XXXX"]),
        ],
        "pZXZX": [
            PauliList(["ZZZZ"]),
            PauliList(["XXXX"]),
            PauliList(["XXXX"]),
            PauliList(["ZZZZ"]),
        ],
        "cXXXX": [
            PauliList(["XXXX"]),
            PauliList(["XXXX"]),
            PauliList(["XXXX"]),
            PauliList(["XXXX"]),
        ],
        "cZZZZ": [
            PauliList(["ZZZZ"]),
            PauliList(["ZZZZ"]),
            PauliList(["ZZZZ"]),
            PauliList(["ZZZZ"]),
        ],
        "cXZZX": [
            PauliList(["XZZX"]),
            PauliList(["XZZX"]),
            PauliList(["XZZX"]),
            PauliList(["XZZX"]),
        ],
        "cZXXZ": [
            PauliList(["ZXXZ"]),
            PauliList(["ZXXZ"]),
            PauliList(["ZXXZ"]),
            PauliList(["ZXXZ"]),
        ],
    }

    # Descriptions of wireframes
    # qubit indices for each wireframe
    wf_q_indices = [[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]]

    # coordinates for each wireframe vertex in path list form to enable the
    # creation of the associate edges
    wf_coordinates = [
        [[-1, 1], [0, 1], [0, 0], [-1, 0]],
        [[0, 1], [1, 1], [1, 0], [0, 0]],
        [[-1, 0], [0, 0], [0, -1], [-1, -1]],
        [[0, 0], [1, 0], [1, -1], [0, -1]],
    ]

    # If the wf's are closed loops or not
    wf_loop_indicator = [True, True, True, True]

    # How wireframes components combine into faces/operators
    faces_wf_components = [[0], [1], [2], [3]]

    # Face colors (wf's inherit colors from faces)
    face_colors = ["yellowgreen", "tomato", "tomato", "yellowgreen"]

    num_faces = 4

    size = np.array([2, 2])

    num_qubits = 9

    u_vec = np.array([2, 0])
    v_vec = np.array([0, 2])

    def __new__(
        cls,
        origin: np.array,
        qubit_count=None,
        qubit_data=None,
        operators=None,
        optype="pXZXZ",
    ) -> Shell:
        """Creates a Checker Board Tile (Shell)

            q0       q1  q1       q2
            v0       v1  v0       v1
              o-----o      o-----o
              |  0  |      |  1  |
              o-----o      o-----o
            v2       v3  v2       v3
            q3       q4  q4       q5

                        .(0,0)

            q3       q4  q4       q5
            v0       v1  v0       v1
              o-----o      o-----o
              |  2  |      |  3  |
              o-----o      o-----o
            v2       v3  v2       v3
            q6       q7  q7       q8


              o-----o-----o
              |  0  |  1  |
              o-----o-----o
              |  2  |  3  |
              o-----o-----o

        Face colors for faces [0,1,2,3] are ["yellowgreen", "tomato", "tomato", "yellowgreen"]

        Preformatted operators are stored in CheckerBoardTile.op_dict. Keys for op_dict are of the
        form [p|c]PPPP where p = pattern and c = copy and P is a Pauli opertor X, Z or Y.

        "pXZXZ" -> #0 face is Pauli("XXXX") operator, #1 face is Pauli("ZZZZ"), etc.
        "cXZZX" -> #0 face is Pauli("XZZX") operator, #1 face is Pauli("XZZX"), etc.

        Available precomputed operator layouts are:

        "pXZXZ", "pZXZX"
        "cXXXX", "cZZZZ", "cXZZX", "cZXXZ"

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

        A appropriately scalled lattice basis for tiling with the CheckerBoardTile can be accessed
        from [CheckerBoardTile.u_vec, CheckerBoardTile.v_vec]

        Raises:
            QiskitError: Unsupported operator type

        Returns:
            Shell: Returns a Checker Board tile (shell) with provided origin
        """

        if operators is None:
            try:
                operators = CheckerBoardTile.wf_operator_dict[optype]
            except KeyError as key_error:
                raise QiskitError(f"Unsupported operator type: {optype}") from key_error

        return TileFactory(
            origin=origin,
            wf_coordinates=CheckerBoardTile.wf_coordinates,
            wf_q_indices=CheckerBoardTile.wf_q_indices,
            wf_loop_indicator=CheckerBoardTile.wf_loop_indicator,
            faces_wf_components=CheckerBoardTile.faces_wf_components,
            num_qubits=CheckerBoardTile.num_qubits,
            face_colors=CheckerBoardTile.face_colors,
            qubit_count=qubit_count,
            qubit_data=qubit_data,
            operators=operators,
        )
