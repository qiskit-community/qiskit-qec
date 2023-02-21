# This code is part of Qiskit.
#
# (C) Copyright IBM 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
=====================================
Geometry (:mod:`qiskit_qec.geometry`)
=====================================

.. currentmodule:: qiskit_qec.geometry

Geometry module classes and functions
=====================================

.. autosummary::
    :toctree: ../stubs/

    GeometryBounds
    Lattice
    Manifold
    Plane
    Shape
    TwoManifold
    Edge
    Face
    QubitCount
    QubitData
    ShapeObject
    Shell
    Vertex
    WireFrame
    Tile
    CheckerBoardTile
    DiagonalBarTile
    DiagonalHourGlassTile
    HexagonTile
    OctaSquareTile
    SquareDiamondTile
    Tiling
    TileFactory
"""

from .bounds import GeometryBounds
from .lattice import Lattice
from .manifold import Manifold
from .plane import Plane
from .shape import Shape
from .two_manifold import TwoManifold

from .model.edge import Edge
from .model.face import Face
from .model.qubit_count import QubitCount
from .model.qubit_data import QubitData
from .model.shape_object import ShapeObject
from .model.shell import Shell
from .model.vertex import Vertex
from .model.wireframe import WireFrame

from .tiles.tile import Tile
from .tiles.checkerboardtile import CheckerBoardTile
from .tiles.diagonalbartile import DiagonalBarTile
from .tiles.diagonalhourglasstile import DiagonalHourGlassTile
from .tiles.hexagontile import HexagonTile
from .tiles.octasquaretile import OctaSquareTile
from .tiles.squarediamondtile import SquareDiamondTile
from .tiles.tiling import Tiling
from .tiles.tilefactory import TileFactory
