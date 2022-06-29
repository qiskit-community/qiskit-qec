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
"""Heavy Hex code builder example"""

from typing import List, Optional
import numpy as np

from qiskit import QiskitError
from qiskit_qec.codes.codebuilders.builder import Builder
from qiskit_qec.codes.codefactory.tilecodefactory import TileCodeFactory
from qiskit_qec.geometry.model.qubit_data import QubitData
from qiskit_qec.geometry.shape import Shape
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.geometry.model.vertex import Vertex
from qiskit_qec.geometry.tiles.diagonalbartile import DiagonalBarTile
from qiskit_qec.geometry.lattice import Lattice
from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode
from qiskit_qec.operators.pauli import Pauli


class HeavyHexCodeBuilder(Builder):
    """Heavy Hex Code Builder Class"""

    # pylint: disable=anomalous-backslash-in-string
    def __init__(
        self,
        d: Optional[int] = None,
        *,
        dx: Optional[int] = None,
        dz: Optional[int] = None,
        w2_op: Optional[Pauli] = Pauli("Z"),
    ) -> None:
        """Initializes a heavy hex code builder

        If d is specified then dx and dz are ignored.

        Args:
            d: distance of code
            dx (optional): X distance of code. Default is d.
            dz (optional): Z distance of code. Default is d.
            w2_op (optional): operator type for weight two gauge operators. Resulting code
            is a CSS code and do only a single Pauli is necessary to define the gauge operator.
            Default is Pauli('Z'). Weight four operators will be Pauli('Z') or Pauli('X') choosen
            to be whatever w2_op are not.

        Examples:
            >>> code = HeavyHexCodeBuilder(d=5).build()

        """
        # Create cutter
        if d is not None:
            if not bool(d % 2) or d < 3:
                raise QiskitError(f"Distance d={d} must be an odd positive integer ≥ 3")
            dx = d
            dz = d
        elif dx is None or dz is None:
            raise QiskitError(f"Both dx:{dx} and dz:{dz} must be speicifed.")
        elif not bool(dx % 2) or dx < 3 or not bool(dz % 2) or dz < 3:
            raise QiskitError(f"dx:{dx} and dz:{dz} must be odd positive integers ≥ 3")

        if w2_op == Pauli("Z"):
            self.w2_op = Pauli("Z")
            self.w4_op = Pauli("X")
            self.optype = "pZZXXZZ"
        else:
            self.w2_op = Pauli("X")
            self.w4_op = Pauli("Z")
            self.optype = "pXXZZXX"

        self.cutter = Shape.rect(
            origin=(-1, -1),
            direction=(1, 0),
            scale1=dx - 1,
            scale2=dz - 1,
            manifold=Plane(),
            delta=0.25,
            dtype=float,
        )

        # Cutter with delta=0 used by the exclude method for boundary selection
        self.cutter_ex = Shape.rect(
            origin=(-1, -1),
            direction=(1, 0),
            scale1=dx - 1,
            scale2=dz - 1,
            manifold=Plane(),
            delta=0,
            dtype=float,
        )

        def exclude(vertex_paths: List[List[Vertex]], qubit_data: QubitData) -> bool:
            def _weight_len(path: List) -> int:
                length = len(path)
                if path[0] == path[-1] and length > 1:
                    return length - 1
                return length

            weights = [_weight_len(path) for path in vertex_paths]
            weight = sum(weights)
            if weight != 2:
                return False
            else:
                v0_pos = vertex_paths[0][0].pos
                v1_pos = vertex_paths[0][1].pos
                if abs(v0_pos[0] - v1_pos[0]) < 0.01:
                    if (
                        abs(v0_pos[0] - self.cutter_ex.bounds.min[0]) < 0.01
                        or abs(v0_pos[0] - self.cutter_ex.bounds.max[0]) < 0.01
                    ):
                        if qubit_data.operator[vertex_paths[0][0].id][0] == self.w4_op:
                            return True
            return False

        self.exclude = exclude

    def build(self) -> StabSubSystemCode:
        """Builds a triangular code code"""
        # Create a code factory
        heaveyhex_code_factory = TileCodeFactory()

        # Configure the code factory
        heaveyhex_code_factory.set_parameters(
            manifold=Plane(),
            tile=DiagonalBarTile,
            tile_optype=self.optype,
            lattice=Lattice(u_vec=DiagonalBarTile.u_vec, v_vec=DiagonalBarTile.v_vec),
            cutter=self.cutter,
            expand_value=np.array([2, 2]),
            on_boundary=False,
            boundary_strategy="combine",
            levels=[2, 4],
            integer_snap=True,
            exclude=self.exclude,
            lattice_view=False,
            precut_tiling_view=False,
            show_qubit_indices=False,
        )

        # Update the factory is_configure check. This is used since we
        # directly updated the TileCodeFactory configuration instead of
        # using the individual TileCodeFactory configuration methods.
        heaveyhex_code_factory.update_is_configure()

        # Create the base triangular color code
        return heaveyhex_code_factory.make_code()
