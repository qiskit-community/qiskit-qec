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

"""OctaSquareTile class."""
from math import sqrt
import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.model.shell import Shell
from qiskit_qec.geometry.tiles.tile import Tile
from qiskit_qec.geometry.tiles.tilefactory import TileFactory
from qiskit_qec.operators.pauli_list import PauliList


# pylint: disable=anomalous-backslash-in-string)
class OctaSquareTile(Tile):
    r"""Octa-Square Tile
                            q0        q1
                            v0        v1
      q2     q3               o------o
      v0     v1         q3   /        \   q4
       o-----o          v7  o          o  v2
       |     |              |          |
       |  0  |          q6  |    1     |  q7
       o-----o          v6  o          o  v3
      v3     v4              \        /
      q5     q6          q9   o-----o  q10
                         v5            v4
     q5        q6
     v0        v1
       o-----o               q9      q10
 q8   /       \   q9         v0      v1
 v7  o         o  v2          o-----o
     |    2    |              |     |
q11  |         |  q12         |  3  |
 v6  o         o  v3          o-----o
`     \       /              q12     q13
       o-----o               v3      v2
    q14        q15
    v5         v4



               o------o
              /        \
       o-----o          o
       |     |          |
       |  0  |    1     |
       o-----o          o
      /       \        /
     o  (0,0)  o-----o
     |    .    |     |
     |    2    |  3  |
     o         o-----o
`     \       /
       o-----o

    """

    wf_operator_dict = {
        "cXZZX": [
            PauliList(["XXXX"], input_qubit_order="left-to-right"),
            PauliList(["ZZZZZZZZ"], input_qubit_order="left-to-right"),
            PauliList(["ZZZZZZZZ"], input_qubit_order="left-to-right"),
            PauliList(["XXXX"], input_qubit_order="left-to-right"),
        ],
        "cZXXZ": [
            PauliList(["ZZZZ"], input_qubit_order="left-to-right"),
            PauliList(["XXXXXXXX"], input_qubit_order="left-to-right"),
            PauliList(["XXXXXXXX"], input_qubit_order="left-to-right"),
            PauliList(["ZZZZ"], input_qubit_order="left-to-right"),
        ],
    }

    # Descriptions of wireframes
    # qubit indices for each wireframe
    wf_q_indices = [
        [2, 3, 6, 5],
        [0, 1, 4, 7, 10, 9, 6, 3],
        [5, 6, 9, 12, 15, 14, 11, 8],
        [9, 10, 13, 12],
    ]

    # pylint: disable=invalid-name
    s8 = sqrt(2 - sqrt(2)) / 2
    c8 = sqrt(2 + sqrt(2)) / 2

    # coordinates for each wireframe vertex in path list form to enable the
    # creation of the associate edges
    wf_coordinates = [
        [[-s8, c8 + 2 * s8], [s8, c8 + 2 * s8], [s8, c8], [-s8, c8]],
        [
            [c8, 2 * c8 + s8],
            [c8 + 2 * s8, 2 * c8 + s8],
            [2 * c8 + s8, c8 + 2 * s8],
            [2 * c8 + s8, c8],
            [c8 + 2 * s8, s8],
            [c8, s8],
            [s8, c8],
            [s8, c8 + 2 * s8],
        ],
        [[-s8, c8], [s8, c8], [c8, s8], [c8, -s8], [s8, -c8], [-s8, -c8], [-c8, -s8], [-c8, s8]],
        [[c8, s8], [c8 + 2 * s8, s8], [c8 + 2 * s8, -s8], [c8, -s8]],
    ]

    # If the wf's are closed loops or not
    wf_loop_indicator = [True, True, True, True]

    # How wireframes components combine into faces/operators
    faces_wf_components = [[0], [1], [2], [3]]

    # Face colors (wf's inherit colors from faces)
    face_colors = ["yellowgreen", "tomato", "tomato", "yellowgreen"]

    num_faces = 4

    size = np.array([3 * c8 + s8, 3 * c8 + s8])

    num_qubits = 16

    u_vec = np.array([2 * (c8 + s8), 0])
    v_vec = np.array([0, 2 * (c8 + s8)])

    def __new__(
        cls,
        origin: np.array,
        qubit_count=None,
        qubit_data=None,
        operators=None,
        optype="cXZZX",
    ) -> Shell:
        r"""Octa Square Tile
                                    q0        q1
                                    v0        v1
              q2     q3               o------o
              v0     v1         q3   /        \   q4
               o-----o          v7  o          o  v2
               |     |              |          |
               |  0  |          q6  |    1     |  q7
               o-----o          v6  o          o  v3
              v3     v4              \        /
              q5     q6          q9   o-----o  q10
                                 v5            v4
             q5        q6
             v0        v1
               o-----o               q9      q10
         q8   /       \   q9         v0      v1
         v7  o         o  v2          o-----o
             |    2    |              |     |
        q11  |         |  q12         |  3  |
         v6  o         o  v3          o-----o
        `     \       /              q12     q13
               o-----o               v3      v2
            q14        q15
            v5         v4



                       o------o
                      /        \
               o-----o          o
               |     |          |
               |  0  |    1     |
               o-----o          o
              /       \        /
             o  (0,0)  o-----o
             |    .    |     |
             |    2    |  3  |
             o         o-----o
        `     \       /
               o-----o

        Face colors for faces [0,1,2,3] are  ["yellowgreen","tomato", "tomato", "yellowgreen"]

        Preformatted operators are stored in HexagonTile.op_dict. Keys for op_dict are of the
        form [p|c|d]PPPP... where p = pattern and c = copy, d=double and P is a Pauli opertor X, Z or Y.

        "cX" -> face #0 operator is Pauli("XXXXXX"), face #1 operator is Pauli("ZZZZZZ"), ...
        "dXZ" -> face #0 operators are Pauli("XXXXXX") and Pauli("ZZZZZZ"), face #1 operators are
        Pauli("ZZZZZZ") and  ...

        Available precomputed operator layouts are:

        cXZZX, cZXXZ

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

        A appropriately scalled lattice basis for tiling with the OctaSquareTile can be accessed
        from [OctaSquareTile.u_vec, OctaSquareTile.v_vec]

        Raises:
            QiskitError: Unsupported operator type

        Returns:
            Shell: Returns a Octa-Square tile (shell) with provided origin
        """

        if operators is None:
            try:
                operators = OctaSquareTile.wf_operator_dict[optype]
            except KeyError as key_error:
                raise QiskitError(f"Unsupported operator type: {optype}") from key_error

        return TileFactory(
            origin=origin,
            wf_coordinates=OctaSquareTile.wf_coordinates,
            wf_q_indices=OctaSquareTile.wf_q_indices,
            wf_loop_indicator=OctaSquareTile.wf_loop_indicator,
            faces_wf_components=OctaSquareTile.faces_wf_components,
            num_qubits=OctaSquareTile.num_qubits,
            face_colors=OctaSquareTile.face_colors,
            qubit_count=qubit_count,
            qubit_data=qubit_data,
            operators=operators,
        )
