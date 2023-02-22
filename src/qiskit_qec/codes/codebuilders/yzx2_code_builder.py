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
"""YZX2 code builder example"""

from qiskit import QiskitError
from qiskit_qec.codes.codebuilders.builder import Builder
from qiskit_qec.codes.codefactory.tilecodefactory import TileCodeFactory
from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode

# Load the appropriate tile: Hexagon Tile
from qiskit_qec.geometry.tiles.hexagontile import HexagonTile

# Load the Shape class to create a cutter
from qiskit_qec.geometry.shape import Shape

# Load the manifold to be tiled
from qiskit_qec.geometry.plane import Plane

# Load the lattice class to tile against
from qiskit_qec.geometry.lattice import Lattice

from qiskit_qec.geometry.model.qubit_data import QubitData


class YZX2CodeBuilder(Builder):
    """YZX2 Code Builder Class"""

    # pylint: disable=anomalous-backslash-in-string
    def __init__(self, d: int) -> None:
        """Initializes a YZX2 code builder

        See https://doi.org/10.22331/q-2022-04-27-698

        Args:
            d: distance of code, d must be an odd positive integer ≥ 3

        Examples:
            >>> code = YZX2CodeBuilder(d=7).build()

        """
        # Set the d parameter
        if not bool(d % 2) or d < 3:
            raise QiskitError(f"Distance d={d} must be an odd positive integer ≥ 3")

        # Define the set of points defining the diamond
        self.scale = d - 1

        # pylint: disable=invalid-name
        h = HexagonTile.h
        r = HexagonTile.r
        self.points = [
            [h, -r / 3],
            [2 * h + self.scale * 3 * h, self.scale * r],
            [h, r / 3 + 2 * r * self.scale],
            [self.scale * (-3 * h), self.scale * r],
        ]

        # Create the diamond cutter
        self.cutter = Shape(points=self.points)

        # Exclude method used to select the required boundary
        # pylint: disable=unused-variable, unused-argument
        def exclude(vertex_paths, qubit_data: QubitData) -> bool:
            """exlude method for Shell class to use for building boundary operators"""

            def _weight_len(path) -> int:
                """Find the weight of the operator from the vertex path listing"""
                length = len(path)
                if path[0] == path[-1] and length > 1:
                    return length - 1
                return length

            weights = [_weight_len(path) for path in vertex_paths]
            weight = sum(weights)

            # Create the equations for the lines bounding the diamond
            # Stored as [A,B,C] for equation Ax+By+C=0
            lines = [0] * 4
            for i in range(4):
                m_grad = (self.points[(i + 1) % 4][1] - self.points[i][1]) / (
                    self.points[(i + 1) % 4][0] - self.points[i][0]
                )
                lines[i] = [1, -m_grad, m_grad * self.points[i][0] - self.points[i][1]]

            def near_line(line, point, epsilon=0.01):
                # Return |Ax+By+C| < epsilon
                return abs(line[0] * point[1] + line[1] * point[0] + line[2]) < epsilon

            def find_indicator(lines_, point):
                return [near_line(line, point) for line in lines_]

            # Ignore any operator that is not of weight 3
            if weight != 3:
                return False
            else:
                triangle_pos = [vertex.pos for vertex in vertex_paths[0]]
                triangle_pos = triangle_pos[:-1]
                indicators = [find_indicator(lines, point) for point in triangle_pos]
                # Find which line the triangle belongs to (if any)
                ind = [i + j + k for i, j, k in zip(*indicators)]
                # Remove the all zero case (not on any lines)
                if sum(ind) == 0:
                    return False
                # Remove the single zero indicator
                if sum(indicators[0]) == 0:
                    on_line_indices = [1, 2]
                elif sum(indicators[1]) == 0:
                    on_line_indices = [0, 2]
                else:
                    on_line_indices = [0, 1]
                # Find which line the triangle is on
                line_indicator = [
                    u & v
                    for u, v in zip(indicators[on_line_indices[0]], indicators[on_line_indices[1]])
                ]
                line_index = [
                    index for index, val in enumerate(line_indicator) if bool(val) is True
                ][0]
                a_y = min(triangle_pos[on_line_indices[0]][1], triangle_pos[on_line_indices[1]][1])
                # TODO: Replace this with a single calculation to avoid branches
                if line_index == 0:
                    k = (a_y) / (2 * r)
                elif line_index == 1:
                    k = (a_y - self.scale * r) / (2 * r)
                elif line_index == 2:
                    k = ((a_y - self.scale * r) / r - 1) / 2
                elif line_index == 3:
                    k = ((a_y / r) - 1) / 2
                if abs(k - round(k)) < 0.01:
                    return False
                return True

        self.exclude = exclude

    def build(self) -> StabSubSystemCode:
        """Builds a YZX2 code"""
        # Create a code factory
        yzx2_code_factory = TileCodeFactory()

        # Configure the code factory
        yzx2_code_factory.set_parameters(
            manifold=Plane(),
            tile=HexagonTile,
            lattice=Lattice(u_vec=HexagonTile.u_vec, v_vec=HexagonTile.v_vec),
            cutter=self.cutter,
            on_boundary=True,
            boundary_strategy="combine",
            inside_levels=[2, 6],
            boundary_levels=[3],
            tile_optype="cYZX2-hXX",
            exclude=self.exclude,
            lattice_view=False,
            precut_tiling_view=False,
            rotate=90,
        )

        # Update the factory is_configure check. This is used since we
        # directly updated the TileCodeFactory configuration instead of
        # using the individual TileCodeFactory configuration methods.
        yzx2_code_factory.update_is_configure()

        # Create the base triangular color code
        return yzx2_code_factory.make_code()
