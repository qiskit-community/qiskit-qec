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

"""Tile Code Factory Class"""

from typing import Optional, Tuple, Dict

import matplotlib.pyplot as plt

from qiskit_qec.codes.codefactory.tilecodeconfig import TileCodeConfig
from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode
from qiskit_qec.structures.gauge import GaugeGroup
from qiskit_qec.operators.pauli import Pauli

from qiskit_qec.geometry.model.qubit_count import QubitCount
from qiskit_qec.geometry.model.qubit_data import QubitData
from qiskit_qec.geometry.shape import Shape
from qiskit_qec.geometry.lattice import Lattice
from qiskit_qec.geometry.tiles.tiling import Tiling
from qiskit_qec.geometry.model.shell import Shell

class TileCodeFactory:
    """ Tile Code Factory Class"""
    def __new__(cls, tile_config:TileCodeConfig) -> None:

        return cls.make_code(tile_config)

    @classmethod
    def make_code(cls, tile_config:TileCodeConfig)->StabSubSystemCode:
        """_summary_

        Args:
            tile_config (TileCodeConfig): _description_

        Returns:
            StabSubSystemCode: _description_
        """
        # Set up the qubit counter and aux data structures
        qubit_count = QubitCount()
        qubit_data = QubitData()

        # Resrict the lattice for tiling
        lattice = tile_config.lattice.restrict_for_tiling(tile_config.cutter,
                                              tile=tile_config.tile,
                                              expand_value=tile_config.expand_value)

        # Display lattice_view if required
        if tile_config.lattice_view:
            if tile_config.cutter_color is None:
                cutter_color="red"
            else:
                cutter_color=tile_config.cutter_color
            cls.plot(lattice = lattice, cutter=tile_config.cutter, cutter_color=cutter_color)

        # Tile the restriced lattice with SquareDiamondTile
        tiling = Tiling(tile_type=tile_config.tile,
                        lattice=lattice,
                        qubit_count=qubit_count,
                        qubit_data=qubit_data)

        # Find the Tiles inside the cutter
        is_inside = tile_config.cutter.inside(tiling, on_boundary=tile_config.on_boundary)
        
        # Display the precut tiling view if required
        if tile_config.precut_tiling_view:
            cls.plot(shell=tiling,
                     qubit_data=qubit_data,
                     is_inside=is_inside,
                     show_index=False,
                     face_colors=True,
                     show_inside=True,
                     cutter=tile_config.cutter,
                     cutter_color=tile_config.cutter_color)

        # pylint: disable=no-member
        new_shell, new_qubit_data, qubit_count = tiling.extract(is_inside,
                                                               qubit_data,
                                                               qubit_count,
                                                               tile_config.levels,
                                                               boundary_strategy=tile_config.boundary_strategy)


        if tile_config.rotate:
            new_shell.rotate2d(tile_config.rotate)

        if tile_config.scale:
            new_shell.scale(tile_config.scale)

        if tile_config.integer_snap:
            new_shell.integer_snap()

        generators, qubit_data = Shell.shell2symplectic(new_shell,
                                                        new_qubit_data,
                                                        qubit_count)

        gauge_group = GaugeGroup(generators)

        return StabSubSystemCode(gauge_group,
                                 shell=new_shell,
                                 qubit_data=new_qubit_data,
                                 qubit_count=qubit_count)

    @classmethod
    def plot(cls,
             lattice:Optional[Lattice]=None,
             cutter:Optional[Shape]=None,
             shell:Optional[Shell]=None,
             qubit_data:Optional[QubitData]=None,
             is_inside:Dict=None,
             show_index:bool=False,
             face_colors:bool=False,
             show_inside:bool=False,
             figsize:Optional[Tuple[float, float]]=None,
             xcolor:str="red",
             zcolor:str="green",
             ycolor:str="blue",
             cutter_color:str="magenta")->None:
        """Plots the intermediate and final stages of tile code generation

        Args:
            lattice (Optional[Lattice], optional): _description_. Defaults to None.
            cutter (Optional[Shape], optional): _description_. Defaults to None.
            shell (Optional[Shell], optional): _description_. Defaults to None.
            qubit_data (Optional[QubitData], optional): _description_. Defaults to None.
            is_inside (Dict, optional): _description_. Defaults to None.
            show_index (bool, optional): _description_. Defaults to False.
            show_inside (bool, optional): _description_. Defaults to False.
            figsize (Optional[Tuple[float, float]], optional): _description_. Defaults to None.
            xcolor (str, optional): _description_. Defaults to "red".
            zcolor (str, optional): _description_. Defaults to "green".
            ycolor (str, optional): _description_. Defaults to "blue".

        Returns:
            _type_: _description_
        """
        if lattice is None and cutter is None and shell is None:
            return None

        if figsize is None:
            figsize=(10,10)

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)

        if lattice is not None:
            data_in = [list(item) for item in lattice.points]
            in_x = [v[0] for v in data_in]
            in_y = [v[1] for v in data_in]
            
            plt.scatter(in_x, in_y, s=100, label="in")
        
        if cutter is not None:

            border = list(cutter.points)
            border.append(border[0])
            xs, ys = zip(*border)

            plt.plot(xs,ys, color=cutter_color, linestyle="-.")

        if shell is not None:
            stab = []
            for face in shell.faces:
                verts = []
                for vertex in face.vertices:
                    verts.append(list(vertex.pos))
                pauli = qubit_data.operator[face.vertices[0].id]
                if face_colors and qubit_data.face_colors != {}:
                    stab.append([verts, qubit_data.face_colors[face.id]])
                else:
                    if pauli == Pauli('X'):
                        stab.append([verts, xcolor])
                    elif pauli == Pauli('Z'):
                        stab.append([verts, zcolor])
                    else:
                        stab.append([verts, ycolor])

            for stabilizer in stab:
                ax.add_patch(plt.Polygon(stabilizer[0], color=stabilizer[1]))

            if show_inside:
                data_in = [list(vertex.pos) for vertex in shell.vertices if is_inside[vertex]]
                in_x = [v[0] for v in data_in]
                in_y = [v[1] for v in data_in]

                data_out = [list(vertex.pos) for vertex in shell.vertices if not is_inside[vertex]]
                out_x = [v[0] for v in data_out]
                out_y = [v[1] for v in data_out]

                plt.scatter(in_x, in_y, s=100, label="in", color='black')
                plt.scatter(out_x, out_y, s=100, label="out")

            if show_index:
                for vertex in shell.vertices:
                    plt.text(vertex.pos[0],vertex.pos[1],qubit_data.index[vertex.id])

        plt.axis('equal')
