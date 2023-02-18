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
"""Rotated Surface code builder example"""

from typing import List, Optional
import numpy as np

from qiskit import QiskitError
from qiskit_qec.codes.codebuilders.builder import Builder
from qiskit_qec.codes.codefactory.tilecodefactory import TileCodeFactory
from qiskit_qec.geometry.model.qubit_data import QubitData
from qiskit_qec.geometry.shape import Shape
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.geometry.model.vertex import Vertex
from qiskit_qec.geometry.tiles.checkerboardtile import CheckerBoardTile
from qiskit_qec.geometry.lattice import Lattice
from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode
from qiskit_qec.operators.pauli import Pauli


class RotatedSurfaceCodeBuilder(Builder):
    """Rotated Surface Code Builder Class"""

    # pylint: disable=anomalous-backslash-in-string
    def __init__(
        self,
        d: Optional[int] = None,
        *,
        dx: Optional[int] = None,
        dz: Optional[int] = None,
        ul_op: Optional[Pauli] = Pauli("Z"),
    ) -> None:
        """Initializes a heavy hex code builder

        If d is specified then dx and dz are ignored.

        Args:
            d: distance of code
            dx (optional): X distance of code. Default is d.
            dz (optional): Z distance of code. Default is d.
            ul_op (optional): operator type for upper left corder weight
                four stabilizer operator.

        Examples:
            >>> code = RotatedSurfaceCodeBuilder(d=9, ul_op=Pauli('X')).build()

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

        if ul_op == Pauli("Z"):
            self.ul_op = Pauli("Z")
            self.nul_op = Pauli("X")
            self.optype = "pZXZX"
        else:
            self.ul_op = Pauli("X")
            self.nul_op = Pauli("Z")
            self.optype = "pXZXZ"

        self.cutter = Shape.rect(
            origin=(-1, -1),
            direction=(1, 0),
            scale1=dx - 1,
            scale2=dz - 1,
            manifold=Plane(),
            delta=0.25,
            dtype=float,
        )

        # Cutter with delta=0 used for boundary selection exclude method
        self.cutter_ex = Shape.rect(
            origin=(-1, -1),
            direction=(1, 0),
            scale1=dx - 1,
            scale2=dz - 1,
            manifold=Plane(),
            delta=0,
            dtype=float,
        )

        # Exclude method used to select the required boundary
        def exclude(vertex_paths: List[List[Vertex]], qubit_data: QubitData) -> bool:
            def _weight_len(path: List) -> int:
                """Find the weight of the operator from the vertex path listing"""
                length = len(path)
                if path[0] == path[-1] and length > 1:
                    return length - 1
                return length

            weights = [_weight_len(path) for path in vertex_paths]
            weight = sum(weights)
            # Exclude any operator that is not of weight 2
            if weight != 2:
                return False
            else:
                v0_pos = vertex_paths[0][0].pos
                v1_pos = vertex_paths[0][1].pos
                if abs(v0_pos[0] - v1_pos[0]) < 0.01:  # vertical line
                    if (
                        abs(v0_pos[0] - self.cutter_ex.bounds.min[0]) < 0.01
                        or abs(v0_pos[0] - self.cutter_ex.bounds.max[0]) < 0.01
                    ):
                        if qubit_data.operator[vertex_paths[0][0].id][0] == self.nul_op:
                            return True
                elif abs(v0_pos[1] - v1_pos[1]) < 0.01:  # horizontal line
                    if (
                        abs(v0_pos[1] - self.cutter_ex.bounds.min[1]) < 0.01
                        or abs(v0_pos[1] - self.cutter_ex.bounds.max[1]) < 0.01
                    ):
                        if qubit_data.operator[vertex_paths[0][0].id][0] == self.ul_op:
                            return True
            return False

        self.exclude = exclude

    def build(self) -> StabSubSystemCode:
        """Builds a rotated surface code code"""
        # Create a code factory
        rotated_surface_code_factory = TileCodeFactory()

        # Configure the code factory
        rotated_surface_code_factory.set_parameters(
            manifold=Plane(),
            tile=CheckerBoardTile,
            tile_optype=self.optype,
            lattice=Lattice(u_vec=CheckerBoardTile.u_vec, v_vec=CheckerBoardTile.v_vec),
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
        rotated_surface_code_factory.update_is_configure()

        # Create the base triangular color code
        return rotated_surface_code_factory.make_code()
