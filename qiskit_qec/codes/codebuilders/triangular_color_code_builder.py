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
"""Triangular color code builder example"""

from qiskit import QiskitError
from qiskit_qec.codes.codebuilders.builder import Builder
from qiskit_qec.codes.codefactory.tilecodefactory import TileCodeFactory
from qiskit_qec.geometry.shape import Shape
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.geometry.tiles.hexagontile import HexagonTile
from qiskit_qec.geometry.lattice import Lattice
from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode


class TriangularColorCodeBuilder(Builder):
    """Triangular Color Code Builder Class"""

    # pylint: disable=anomalous-backslash-in-string
    def __init__(self, d: int) -> None:
        """Initializes a triangular color code builder

        Example:d=3
                        o
                       / \
                      /   \
                     /     o
                    /     / \
                   /     /   \
                  o-----o     \
                 /       \     \
                /         \     \
               o-----------o-----o

        Args:
            d: distance of code

        Examples:
            >>> code = TriangularColorCodeBuilder(d=7).build()

        """
        # Create cutter
        if not bool(d % 2) or d < 3:
            raise QiskitError(f"Distance d={d} must be an odd positive integer ≥ 3")
        scale = 3 * (d - 1) / 2
        delta = 0.02
        points = [
            [-delta, -delta],
            [scale * HexagonTile.h, scale * HexagonTile.r + delta],
            [scale + delta, -delta],
        ]
        self.cutter = Shape(points=points)

    def build(self) -> StabSubSystemCode:
        """Builds a triangular color code"""
        # Create a code factory
        triangular_code_factory = TileCodeFactory()

        # Configure the code factory
        triangular_code_factory.set_parameters(
            manifold=Plane(),
            tile=HexagonTile,
            lattice=Lattice(u_vec=HexagonTile.u_vec, v_vec=HexagonTile.v_vec),
            cutter=self.cutter,
            on_boundary=False,
            boundary_strategy="combine",
            levels=[4, 6],
            tile_optype="dXZ",
        )

        # Update the factory is_configure check. This is used since we
        # directly updated the TileCodeFactory configuration instead of
        # using the individual TileCodeFactory configuration methods.
        triangular_code_factory.update_is_configure()

        # Create the base triangular color code
        return triangular_code_factory.make_code()
