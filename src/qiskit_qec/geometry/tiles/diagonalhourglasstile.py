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

"""DiagonalHourGlassTile class."""

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.model.shell import Shell
from qiskit_qec.geometry.tiles.tile import Tile
from qiskit_qec.geometry.tiles.tilefactory import TileFactory
from qiskit_qec.operators.pauli_list import PauliList


# pylint: disable=anomalous-backslash-in-string)
class DiagonalHourGlassTile(Tile):
    r"""Diagonal Hour Glass Tile

    q0           q1     q1           q2
    v0           v0     v0           v1
     o- - - - - -o       o- - - - - -o
     |\   q3    /|       |           |
     |  \ v1  /  |       |           |
     | 0 o   o 1 |       |     2     |
     |  /   v2\  |       |           |
     |/     q3  \|       |           |
     o- - - - - -o       o- - - - - -o
    v2           v1     v3           v2
    q4           q5     q5           q6


    q4           q5     q5           q6
    v0           v1     v0           v0
     o- - - - - -o      o- - - - - -o
     |           |      |\   q7    /|
     |           |      |  \ v1  /  |
     |     3     |      | 4 o   o 5 |
     |           |      |  /   v2\  |
     |           |      |/     q7  \|
     o- - - - - -o      o- - - - - -o
    v3           v2    v2           v1
    v8           v9    q9           q10
    """

    wf_operator_dict = {
        "pXXZZXX": [
            PauliList(["XXX"]),
            PauliList(["XXX"]),
            PauliList(["ZZZZ"]),
            PauliList(["ZZZZ"]),
            PauliList(["XXX"]),
            PauliList(["XXX"]),
        ],
        "pZZXXZZ": [
            PauliList(["ZZZ"]),
            PauliList(["ZZZ"]),
            PauliList(["XXXX"]),
            PauliList(["XXXX"]),
            PauliList(["ZZZ"]),
            PauliList(["ZZZ"]),
        ],
    }

    # Descriptions of wireframes
    # qubit indices for each wireframe
    wf_q_indices = [[0, 3, 4], [1, 5, 3], [1, 2, 6, 5], [4, 5, 9, 8], [5, 7, 9], [6, 10, 7]]

    # coordinates for each wireframe vertex in path list form to enable the
    # creation of the associate edges
    wf_coordinates = [
        [[-2, 2], [-1, 1], [-2, 0]],
        [[0, 2], [0, 0], [-1, 1]],
        [[0, 2], [2, 2], [2, 0], [0, 0]],
        [[-2, 0], [0, 0], [0, -2], [-2, -2]],
        [[0, 0], [1, -1], [0, -2]],
        [[2, 0], [2, -2], [1, -1]],
    ]

    # If the wf's are closed loops or not
    wf_loop_indicator = [True, True, True, True, True, True]

    # How wireframes components combine into faces/operators
    faces_wf_components = [[0], [1], [2], [3], [4], [5]]

    # Face colors (wf's inherit colors from faces)
    face_colors = ["yellowgreen", "yellowgreen", "tomato", "tomato", "yellowgreen", "yellowgreen"]

    num_faces = 6

    size = np.array([4, 4])

    num_qubits = 11

    u_vec = np.array([4, 0])
    v_vec = np.array([0, 4])

    def __new__(
        cls,
        origin: np.array,
        qubit_count=None,
        qubit_data=None,
        operators=None,
        optype="pXXZZXX",
    ) -> Shell:
        r"""Diagonal Hour Glass Tile


            o- - - - - -o- - - - - -o
            |\         /|           |
            |  \     /  |           |
            | 0 o   o 1 |     2     |
            |  /     \  |           |
            |/         \|           |
            o- - - - - -o- - - - - -o
            |           |\         /|
            |           |  \     /  |
            |     3     | 4 o   o 5 |
            |           |  /     \  |
            |           |/         \|
            o- - - - - -o- - - - - -o


           q0           q1     q1           q2
           v0           v0     v0           v1
            o- - - - - -o       o- - - - - -o
            |\   q3    /|       |           |
            |  \ v1  /  |       |           |
            | 0 o   o 1 |       |     2     |
            |  /   v2\  |       |           |
            |/     q3  \|       |           |
            o- - - - - -o       o- - - - - -o
           v2           v1     v3           v2
           q4           q5     q5           q6


           q4           q5     q5           q6
           v0           v1     v0           v0
            o- - - - - -o      o- - - - - -o
            |           |      |\   q7    /|
            |           |      |  \ v1  /  |
            |     3     |      | 4 o   o 5 |
            |           |      |  /   v2\  |
            |           |      |/     q7  \|
            o- - - - - -o      o- - - - - -o
           v3           v2    v2           v1
           v8           v9    q9           q10


        Face colors for faces [0,1,2,3,4,5] are ["yellowgreen", "yellowgreen",
        "tomato", "tomato", "yellowgreen","yellowgreen"]

        Preformatted operators are stored in DiagonalHourGlassTile.op_dict.
        Keys for op_dict are of the form [p|c]PPPP... where p = pattern and
        c = copy and P is a Pauli opertor X, Z, Y.

        "pXXZZXX" -> #0 face is Pauli('XXX') operator,
                     #1 face is Pauli('XXX'),
                     #3 face is Pauli('ZZZZ') etc.

        Available precomputed operator layouts are:

        "pXXZZXX", "pZZXXZZ"

        The operator variable may be used to define the operators specifically.
        The operator must be a list of PauliList objects where each PauliList
        describes the opertors to be built for the faces as indexed above 0,1,2,3, ...
        If the PauliList contains k Paulis then k operators will be created for
        the given face.

        Args:
            origin (np.array): Coordinates of origin of tile (shell)
            qubit_count: Qubit counter. Defaults to None.
            qubit_data: Qubit data. Defaults to None.
            operators: Operators for tile faces. Defaults to None.
            optype (optional): Which of the listed opertor mapppings to used.
                Defaults to "pXXZZXX".

        Raises:
            QiskitError: Unsupported operator type

        Returns:
            Shell: Returns a Diagonal Hour Glass tile (shell) with provided origin
        """

        if operators is None:
            try:
                operators = DiagonalHourGlassTile.wf_operator_dict[optype]
            except KeyError as key_error:
                raise QiskitError(f"Unsupported operator type: {optype}") from key_error

        return TileFactory(
            origin=origin,
            wf_coordinates=DiagonalHourGlassTile.wf_coordinates,
            wf_q_indices=DiagonalHourGlassTile.wf_q_indices,
            wf_loop_indicator=DiagonalHourGlassTile.wf_loop_indicator,
            faces_wf_components=DiagonalHourGlassTile.faces_wf_components,
            num_qubits=DiagonalHourGlassTile.num_qubits,
            face_colors=DiagonalHourGlassTile.face_colors,
            qubit_count=qubit_count,
            qubit_data=qubit_data,
            operators=operators,
        )
