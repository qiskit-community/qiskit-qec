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
# Part of the QEC framework

"""Tile Code Config Class"""

from typing import Optional, List, Union

import numpy as np
from qiskit import QiskitError

from qiskit_qec.geometry.shape import Shape
from qiskit_qec.geometry.manifold import Manifold
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.geometry.tiles.tile import Tile
from qiskit_qec.geometry.lattice import Lattice

class TileCodeConfig:
    """Tile Code Config Class for regalar tiling codes"""
    def __init__(self, name:Optional[str]=None, manifold:Optional[Manifold]=Plane()):
        """Initializes the TileCodeConfig class

        Args:
            name optional): Name for TileCode. Defaults to None.
            manifold (optional): Manifold that tile code will live on. Defaults to Plane().
        """
        self.name = name                    # Name 
        self.manifold = manifold
        self.cutter = None
        self.tile = None
        self.boundary_strategy = None
        self.on_boundary = None
        self.levels = None
        self.rotate = None
        self.scale = None
        self.integer_snap = None
        self.lattice = None
        self.expand_value = None
        self.lattice_view = None
        self.precut_tiling_view = None
        self.final_view = None
        self.show_qubit_indices = None
        self.tiling_view = None
        self.cutter_color = None


    def set_cutter(self,
                   points:Optional[List]=None,
                   cutter:Optional[Shape]=None,
                   )->None:
        """Sets the cutter for the tiling

        Args:
            points (optional): Points that describe the cutter shape. Defaults to None.
            cutter (optional): Shape that descibes the cutter shape. Defaults to None.
        """
        if points is not None:
            self.cutter = Shape(points)
        if cutter is not None:
            self.cutter = cutter


    def set_tile(self, tile:Tile)->None:
        """Set the Tile for the tile code

        Args:
            tile: Tile for tile code
        """
        self.tile = tile

    def set_boundary_strategy(self,
                              boundary_strategy:Optional[str]="combine",
                              on_boundary:Optional[bool]=True,
                              levels:Optional[List]=None)->None:
        """Sets how the boundary for the cutter behaves on stabilizers

        Args:
            boundary_strategy (optional): Stategy for how to deal with an operator
            cut by cutter. Defaults to "combine".
            on_boundary (optional): Whether the boundary is in or out of the inside
            of the cutter. Defaults to True.
            levels (optional): Which weights of curt operators to include. Defaults to [2,3,4]
        """
        self.boundary_strategy = boundary_strategy
        self.on_boundary = on_boundary
        if levels is None:
            levels = [2,3,4]
        self.levels = levels

    def set_cutter_process(self,
                           expand_value:Union[List,np.ndarray]=None)->None:
        """Sets how the cutter is 

        Args:
            expand_value (pptional): Sets how the cutter is expanded. Defaults to [1,1]
        """
        if expand_value is None:
            expand_value = [1,1]
        self.expand_value = np.array(expand_value)
    
    def set_post_process(self,
                         rotate:Optional[int]=0,
                         scale:Optional[float]=1,
                         integer_snap:Optional[bool]=False)->None:
        """Sets how to post process the code

        Args:
            rotate (optional): How much to rotate the code geometry. Defaults to 0.
            scale (optional): How much to scale the code geometry. Defaults to 1.
            integer_snap (optional): Whether to snap the vertex coordinates to
            integer values. Defaults to False.
        """
        self.rotate = rotate
        self.scale = scale
        self.integer_snap = integer_snap

    def set_lattice(self, *, lattice:Optional[Lattice]=None, basis:Optional[List[List]]=None)->None:
        """Set the lattice for the tile code

        Args:
            lattice: Lattice for the tile code
        """
        if lattice is None:
            if basis is None:
                raise QiskitError(f"At least one of lattice or basis must be provided")
            u_vec = np.array(basis[0])
            v_vec = np.array(basis[1])
            lattice = Lattice(u_vec, v_vec)
        self.lattice = lattice

    def set_view_options(self,
                         lattice_view:Optional[bool]=False,
                         precut_tiling_view:Optional[bool]=False,
                         final_view:Optional[bool]=False,
                         show_qubit_indices:Optional[bool]=False,
                         tiling_view:Optional[bool]=False
                         )->None:
        """Sets the viewing options for debugging and viewing final code graphically

        Args:
            lattice_view (optional): View lattice before cutting (debug feature).
            Defaults to False.
            lattice_view (optional): View reult before post processing
            (debug feature). Defaults to False.
            final_view (optional): View the final code automatically. Defaults
            to False.
            show_qubit_indices (optional): Show qubit indices. Defaults to False.
            tiling_view (optional): View tiling before cutter is applied
        """
        self.lattice_view = lattice_view
        self.precut_tiling_view = precut_tiling_view
        self.final_view = final_view
        self.show_qubit_indices = show_qubit_indices
        self.tiling_view = tiling_view

    def set_colors(self,
                   cutter_color:Optional[str]=None)->None:
        self.cutter_color = cutter_color

