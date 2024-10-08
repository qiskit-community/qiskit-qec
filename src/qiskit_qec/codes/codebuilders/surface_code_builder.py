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
"""Surface code builder example"""

from typing import Optional
from math import sqrt
import numpy as np

from qiskit import QiskitError
from qiskit_qec.codes.codebuilders.builder import Builder
from qiskit_qec.codes.codefactory.tilecodefactory import TileCodeFactory
from qiskit_qec.geometry.shape import Shape
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.geometry.tiles.checkerboardtile import CheckerBoardTile
from qiskit_qec.geometry.lattice import Lattice
from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode
from qiskit_qec.operators.pauli import Pauli


class SurfaceCodeBuilder(Builder):
    """Surface Code Builder Class"""

    # pylint: disable=anomalous-backslash-in-string
    def __init__(
        self,
        d: Optional[int] = None,
        *,
        dx: Optional[int] = None,
        dz: Optional[int] = None,
        ul_op: Optional[Pauli] = Pauli("Z"),
    ) -> None:
        r"""Initializes a surface code builder

        If d is specified then dx and dz are ignored.

        Args:
            d: distance of code
            dx (optional): X distance of code. Default is d.
            dz (optional): Z distance of code. Default is d.
            ul_op (optional): operator type for upper left side corner weight
                three stabilizer operator.

        Examples:
            >>> code = SurfaceCodeBuilder(dx=5, dz=7).build()

            >>> code = SurfaceCodeBuilder(d=9).build()

            >>> code = SurfaceCodeBuilder(dx=5, dz=7, ul_op=Pauli('Z')).build()

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
            origin=(0, -1),
            direction=(1, 1),
            scale1=dx - 1,
            scale2=dz - 1,
            manifold=Plane(),
            delta=0.25,
            dtype=float,
        )

    def build(self) -> StabSubSystemCode:
        """Builds a surface code"""
        # Create a code factory
        surface_code_factory = TileCodeFactory()

        # Configure the code factory
        surface_code_factory.set_parameters(
            manifold=Plane(),
            tile=CheckerBoardTile,
            tile_optype=self.optype,
            lattice=Lattice(u_vec=CheckerBoardTile.u_vec, v_vec=CheckerBoardTile.v_vec),
            cutter=self.cutter,
            expand_value=np.array([1, 1]),
            on_boundary=False,
            boundary_strategy="combine",
            levels=[3, 4],
            integer_snap=True,
            rotate=45,
            scale=sqrt(2),
            exclude=None,
            lattice_view=False,
            precut_tiling_view=False,
            show_qubit_indices=False,
        )

        # Update the factory is_configure check. This is used since we
        # directly updated the TileCodeFactory configuration instead of
        # using the individual TileCodeFactory configuration methods.
        surface_code_factory.update_is_configure()

        # Create the base triangular color code
        return surface_code_factory.make_code()
