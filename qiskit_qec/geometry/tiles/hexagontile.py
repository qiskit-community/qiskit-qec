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

"""HexagonTile class."""
from math import sqrt
import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.model.shell import Shell
from qiskit_qec.geometry.tiles.tile import Tile
from qiskit_qec.geometry.tiles.tilefactory import TileFactory
from qiskit_qec.operators.pauli_list import PauliList


# pylint: disable=anomalous-backslash-in-string)
class HexagonTile(Tile):
    """Hexagon Tile

                          q0      q1
                          v0      v1
                           o-----o
                      q3  /       \  q4
                      v5 o    0    o v2
        q2      q3        \       /
        v0      v1      q6 o-----o q7
         o-----o        v4         v3
    q5  /       \  q6  .(0,0)
    v5 o    1    o v2    q6       q7
        \       /        v0       v1
         o-----o           o-----o
        v4      v3    q9  /       \ q10
        q8      q9    v5 o    2    o v2
                          \       /
                           o-----o
                         v4       v3
                        q11       q12


                      o-----o
                     /       \
              o-----o    0    o
             /       \       /
            o    1    o-----o
             \       /       \
              o-----o    2    o
                     \       /
                      o-----o

    """

    wf_operator_dict = {
        "cX" : [PauliList(["XXXXXX"], input_qubit_order="left-to-right"),
                PauliList(["XXXXXX"], input_qubit_order="left-to-right"),
                PauliList(["XXXXXX"], input_qubit_order="left-to-right")],
        "cZ" : [PauliList(["ZZZZZZ"], input_qubit_order="left-to-right"),
                PauliList(["ZZZZZZ"], input_qubit_order="left-to-right"),
                PauliList(["ZZZZZZ"], input_qubit_order="left-to-right")],
        "dXZ": [PauliList(["XXXXXX", "ZZZZZZ"], input_qubit_order="left-to-right"),
                PauliList(["XXXXXX", "ZZZZZZ"], input_qubit_order="left-to-right"),
                PauliList(["XXXXXX", "ZZZZZZ"], input_qubit_order="left-to-right")],
        "dZX": [PauliList(["ZZZZZZ", "XXXXXX"], input_qubit_order="left-to-right"),
                PauliList(["ZZZZZZ", "XXXXXX"], input_qubit_order="left-to-right"),
                PauliList(["ZZZZZZ", "XXXXXX"], input_qubit_order="left-to-right")]     
    }

    # Descriptions of wireframes
    # qubit indices for each wireframe
    wf_q_indices = [[0,1,4,7,6,3],[2,3,6,9,8,5],[6,7,10,12,11,9]]

    # pylint: disable=invalid-name
    r = sqrt(3)/2
    h = 1/2

    # coordinates for each wireframe vertex in path list form to enable the 
    # creation of the associate edges
    wf_coordinates = [
        [[0 , 2 * r], [2 * h, 2 * r], [1 + h, r], [1, 0], [0, 0], [-h, r]],
        [[-1 - h, r],[-h, r], [0, 0],[-h,-r],[-1-h,-r],[-2,0]],
        [[0,0],[1,0],[1+h,-r],[2*h,-2*r],[0,-2*r],[-h,-r]]
    ]

    # If the wf's are closed loops or not
    wf_loop_indicator = [True, True, True]

    # How wireframes components combine into faces/operators
    faces_wf_components = [[0], [1], [2]]

    # Face colors (wf's inherit colors from faces)
    face_colors = ["yellowgreen","tomato","steelblue"]

    num_faces = 3

    size = np.array([2 * h + 2, 4 * r])

    num_qubits = 13

    u_vec = np.array([3, 0])
    v_vec = np.array([3/2,3*r])

    def __new__(
        cls,
        origin: np.array,
        qubit_count=None,
        qubit_data=None,
        operators=None,
        optype="dXZ",
    ) -> Shell:
        """Hexagon Tile
                                     q0      q1
                                     v0      v1
                                       o-----o
                                  q3  /       \  q4
                                  v5 o    0    o v2
                    q2      q3        \       /
                    v0      v1      q6 o-----o q7
                     o-----o        v4         v3
                q5  /       \  q6  .(0,0)
                v5 o    1    o v2    q6       q7
                    \       /        v0       v1
                     o-----o           o-----o
                    v4      v3    q9  /       \ q10
                    q8      q9    v5 o    2    o v2
                                      \       /
                                       o-----o
                                     v4       v3
                                    q11       q12


                                 o-----o
                                /       \
                         o-----o    0    o
                        /       \       /
                       o    1    o-----o
                        \       /       \
                         o-----o    2    o
                                \       /
                                 o-----o

        Face colors for faces [0,1,2] are ["yellowgreen","tomato","steelblue"]

        Preformatted operators are stored in HexagonTile.op_dict. Keys for op_dict are of the
        form [p|c|d]PPPP... where p = pattern and c = copy, d=double and P is a Pauli opertor X, Z or Y.

        "cX" -> face #0 operator is Pauli("XXXXXX"), face #1 operator is Pauli("ZZZZZZ"), ...
        "dXZ" -> face #0 operators are Pauli("XXXXXX") and Pauli("ZZZZZZ"), face #1 operators are
        Pauli("ZZZZZZ") and  ...

        Available precomputed operator layouts are:

        "cX", "cZ", "dXZ"

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
            Shell: Returns a Checker Board tile (shell) with provided origin
        """

        if operators is None:
            try:
                operators = HexagonTile.wf_operator_dict[optype]
            except KeyError as key_error:
                raise QiskitError(f"Unsupported operator type: {optype}") from key_error

        return TileFactory(origin=origin,
                           wf_coordinates=HexagonTile.wf_coordinates,
                           wf_q_indices=HexagonTile.wf_q_indices,
                           wf_loop_indicator=HexagonTile.wf_loop_indicator,
                           faces_wf_components=HexagonTile.faces_wf_components,
                           num_qubits=HexagonTile.num_qubits,
                           face_colors=HexagonTile.face_colors,
                           qubit_count=qubit_count,
                           qubit_data=qubit_data,
                           operators=operators,
                           )