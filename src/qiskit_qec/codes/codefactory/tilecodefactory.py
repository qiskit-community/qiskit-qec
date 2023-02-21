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

import numpy as np

import matplotlib.pyplot as plt
from qiskit import QiskitError

from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode
from qiskit_qec.operators.pauli import Pauli

from qiskit_qec.geometry.model.qubit_count import QubitCount
from qiskit_qec.geometry.model.qubit_data import QubitData
from qiskit_qec.geometry.shape import Shape
from qiskit_qec.geometry.lattice import Lattice
from qiskit_qec.geometry.tiles.tiling import Tiling
from qiskit_qec.geometry.model.shell import Shell


class TileCodeFactory:
    """ "Tile Code Factory"""

    MANIFOLD_PARAMETERS = ["manifold"]
    TILE_PARAMETERS = ["tile", "tile_optype"]
    LATTICE_PARAMETERS = ["lattice"]
    CUTTER_PARAMETERS = ["cutter", "on_boundary"]
    EXCLUDE_PARAMETERS = ["exclude"]
    BOUNDARY_PARAMETERS = ["boundary_strategy", "levels", "inside_levels", "boundary_levels"]
    TILING_PARAMETERS = ["expand_value"]
    DEBUG_VIEWS_PARAMETERS = ["lattice_view", "precut_tiling_view", "final_view"]
    DRAWING_PARAMETERS = [
        "cutter_color",
        "show_qubit_indices",
        "show_face_indices",
        "show_face_colors",
        "show_inside",
        "show_face_ids",
    ]
    POST_PROCESS_PARAMETERS = ["rotate", "scale", "integer_snap"]

    PARAMETER_DICT = {
        "manifold": MANIFOLD_PARAMETERS,
        "tile": TILE_PARAMETERS,
        "lattice": LATTICE_PARAMETERS,
        "cutter": CUTTER_PARAMETERS,
        "exclude": EXCLUDE_PARAMETERS,
        "boundary": BOUNDARY_PARAMETERS,
        "tiling": TILING_PARAMETERS,
        "debug": DEBUG_VIEWS_PARAMETERS,
        "drawing": DRAWING_PARAMETERS,
        "post": POST_PROCESS_PARAMETERS,
    }

    PARAMETERS = (
        MANIFOLD_PARAMETERS
        + TILE_PARAMETERS
        + LATTICE_PARAMETERS
        + CUTTER_PARAMETERS
        + EXCLUDE_PARAMETERS
        + BOUNDARY_PARAMETERS
        + TILING_PARAMETERS
        + DEBUG_VIEWS_PARAMETERS
        + DRAWING_PARAMETERS
        + POST_PROCESS_PARAMETERS
    )

    """ Tile Code Factory Class"""

    def __init__(self, name: Optional[str] = None) -> None:
        """Initializes a TileCodeFactory class

        Args:
            name (optional): Name to associate with Factory.
            Defaults to 'TileCodeFactory'.
        """
        self.name = name if name is not None else "TileCodeFactory"

        self.configured = {}

        # Manifold setting
        self.configured["manifold"] = False
        self.manifold = None

        # Tile setting
        self.configured["tile"] = False
        self.tile = None
        self.tile_optype = None

        # Lattice setting
        self.configured["lattice"] = False
        self.lattice = None

        # Cutter
        self.configured["cutter"] = False
        self.cutter = None
        self.on_boundary = None

        # Exclude
        self.exclude = None

        # Boundary settings
        self.configured["boundary"] = False
        self.boundary_strategy = None
        self.levels = None
        self.inside_levels = None
        self.boundary_levels = None

        # Tiling settings
        self.expand_value = np.array([1, 1])

        # Debug Views
        self.lattice_view = False
        self.precut_tiling_view = False
        self.final_view = False

        # Drawing settings
        self.cutter_color = "black"
        self.show_qubit_indices = False
        self.show_face_indices = False
        self.show_face_colors = True
        self.show_inside = True
        self.show_qubit_ids = False
        self.show_face_ids = True

        # Post process settings
        self.rotate = False
        self.scale = False
        self.integer_snap = False

    def is_configured(self) -> bool:
        """Check if factory is fully configured"""
        return all(value is True for value in self.configured.values())

    def update_is_configure(self):
        """Update the is_configure method for parameters set directly"""
        for para_set in self.configured:
            self.configured[para_set] = all(
                getattr(self, parameter) is not None
                for parameter in TileCodeFactory.PARAMETER_DICT[para_set]
            )
        if not self.configured["boundary"]:
            if self.boundary_strategy is not None:
                if (
                    self.levels is None
                    and self.inside_levels is not None
                    and self.boundary_levels is not None
                ) or (self.levels is not None):
                    self.configured["boundary"] = True

    def set_parameters(self, **kwargs):
        """Set the given Factory parameters

        Raises:
            QiskitError: Unknown parameter
        """
        for key, value in kwargs.items():
            if key not in TileCodeFactory.PARAMETERS:
                raise QiskitError(f"Unknown parameter: {key}")
            setattr(self, key, value)

    def show_parameters(self, parameter_set="all") -> None:
        """Returns the Factory parameters as a string

        Args:
            parameter_set (optional): Sets which groups of parameters to display.
            Defaults to "all". Possible options can be determined using the
            show_parameter_groups() method

        Raises:
            QiskitError: Unknown parameter set

        Returns:
            str: string showing the appropriate factory parameters
        """
        out_str = ""
        if parameter_set == "all":
            for parameter in TileCodeFactory.PARAMETERS:
                out_str += f"{parameter:20} : { getattr(self, parameter)}\n"

        elif parameter_set in TileCodeFactory.PARAMETER_DICT:
            for parameter in TileCodeFactory.PARAMETER_DICT[parameter_set]:
                out_str += f"{parameter:20} : { getattr(self, parameter)}\n"
        else:
            raise QiskitError(f"Unknown parameter set {parameter_set}")

        print(out_str)

    @classmethod
    def show_parameter_groups(cls) -> None:
        """Returns the parameter group names

        Returns:
            str: Parameter group names
        """
        print(list(TileCodeFactory.PARAMETER_DICT.keys()))

    def set_cutter(self, cutter: Shape) -> None:
        """Sets the cutter for the tiling

        Args:
            cutter: Shape that descibes the cutter shape. Defaults to None.
        """
        self.cutter = cutter

    def make_code(self) -> StabSubSystemCode:
        """Make a code using the configured TileCodeFactory

        Returns:
            StabSubSystemCode: Stababilizer Subsystem Code specified by
            factory configuration
        """

        if not self.is_configured():
            not_configured = {
                key: value for key, value in self.configured.items() if value is False
            }
            raise QiskitError(
                f"Factory is not fully configured: {not_configured}. If"
                + " the configuation was updated directly or using the set_parameters"
                + " method then first run the update_is_configure method to update the"
                + " configuration checker manually."
            )

        # Set up the qubit counter and aux data structures
        qubit_count = QubitCount()
        qubit_data = QubitData()

        # Resrict the lattice for tiling
        lattice = self.lattice.restrict_for_tiling(
            self.cutter, tile=self.tile, expand_value=self.expand_value
        )

        # Display lattice_view if required
        if self.lattice_view:
            self.plot(lattice=lattice, cutter=self.cutter, cutter_color=self.cutter_color)

        # Tile the restriced lattice with configured Tile (Shell)
        tiling = Tiling(
            tile_type=self.tile,
            tile_optype=self.tile_optype,
            lattice=lattice,
            qubit_count=qubit_count,
            qubit_data=qubit_data,
        )

        # Find the Tiles inside the cutter
        is_inside = self.cutter.inside(tiling, on_boundary=self.on_boundary)

        # return tiling, qubit_data, qubit_count

        # Display the precut tiling view if required
        if self.precut_tiling_view:
            self.plot(
                shell=tiling,
                qubit_data=qubit_data,
                is_inside=is_inside,
                show_qubit_indices=self.show_qubit_indices,
                show_qubit_ids=self.show_qubit_ids,
                show_face_ids=self.show_face_ids,
                face_colors=self.show_face_colors,
                show_inside=self.show_inside,
                cutter=self.cutter,
                cutter_color=self.cutter_color,
            )

        # pylint: disable=no-member
        new_shell, new_qubit_data, qubit_count = tiling.extract(
            is_inside,
            qubit_data,
            qubit_count,
            levels=self.levels,
            inside_levels=self.inside_levels,
            boundary_levels=self.boundary_levels,
            exclude=self.exclude,
            boundary_strategy=self.boundary_strategy,
        )

        if self.rotate:
            new_shell.rotate2d(self.rotate)

        if self.scale:
            new_shell.scale(self.scale)

        if self.integer_snap:
            new_shell.integer_snap()

        return StabSubSystemCode(
            shell=new_shell, qubit_data=new_qubit_data, qubit_count=qubit_count
        )

    @classmethod
    def plot(
        cls,
        lattice: Optional[Lattice] = None,
        cutter: Optional[Shape] = None,
        shell: Optional[Shell] = None,
        qubit_data: Optional[QubitData] = None,
        is_inside: Dict = None,
        show_qubit_indices: bool = False,
        face_colors: bool = False,
        show_inside: bool = False,
        figsize: Optional[Tuple[float, float]] = None,
        show_qubit_ids: bool = False,
        show_face_ids: bool = True,
        xcolor: str = "red",
        zcolor: str = "green",
        ycolor: str = "blue",
        cutter_color: str = "magenta",
    ) -> None:
        """Plots the intermediate and final stages of tile code generation

        Args:
            lattice (Optional[Lattice], optional): _description_. Defaults to None.
            cutter (Optional[Shape], optional): _description_. Defaults to None.
            shell (Optional[Shell], optional): _description_. Defaults to None.
            qubit_data (Optional[QubitData], optional): _description_. Defaults to None.
            is_inside (Dict, optional): _description_. Defaults to None.
            show_qubit_indices (bool, optional): _description_. Defaults to False.
            face_colors (bool, optional): _description_. Defaults to False.
            show_inside (bool, optional): _description_. Defaults to False.
            figsize (Optional[Tuple[float, float]], optional): _description_. Defaults to None.
            show_qubit_ids: Show qubit ids
            show_face_ids: Show face ids
            xcolor (str, optional): _description_. Defaults to "red".
            zcolor (str, optional): _description_. Defaults to "green".
            ycolor (str, optional): _description_. Defaults to "blue".
            cutter_color (str, optional): _description_. Defaults to "magenta".

        Raises:
            QiskitError: No shell, cutter or lattive provided

        """
        if lattice is None and cutter is None and shell is None:
            raise QiskitError("No shell, cutter or lattive provided")

        if figsize is None:
            figsize = (10, 10)

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

            plt.plot(xs, ys, color=cutter_color, linestyle="-.")

        if shell is not None:
            stab = []
            for face in shell.faces:
                verts = []
                for vertex in face.vertices:
                    verts.append(list(vertex.pos))
                # This is the cheap persons version of deciding on the color.
                pauli = qubit_data.operator[face.vertices[0].id][0]
                if face_colors and qubit_data.face_colors != {}:
                    stab.append([verts, qubit_data.face_colors[face.id]])
                else:
                    if pauli == Pauli("X"):
                        stab.append([verts, xcolor])
                    elif pauli == Pauli("Z"):
                        stab.append([verts, zcolor])
                    else:
                        stab.append([verts, ycolor])

            for stabilizer in stab:
                if len(stabilizer[0]) == 2:
                    ax.add_patch(plt.Polygon(stabilizer[0], color=stabilizer[1], linewidth=7))
                else:
                    ax.add_patch(plt.Polygon(stabilizer[0], color=stabilizer[1]))

            if show_inside:
                data_in = [list(vertex.pos) for vertex in shell.vertices if is_inside[vertex]]
                in_x = [v[0] for v in data_in]
                in_y = [v[1] for v in data_in]

                data_out = [list(vertex.pos) for vertex in shell.vertices if not is_inside[vertex]]
                out_x = [v[0] for v in data_out]
                out_y = [v[1] for v in data_out]

                plt.scatter(in_x, in_y, s=100, label="in", color="black")
                plt.scatter(out_x, out_y, s=100, label="out")

            if show_qubit_indices and qubit_data.index != {}:
                for vertex in shell.vertices:
                    plt.text(vertex.pos[0], vertex.pos[1], qubit_data.index[vertex.id])
            if show_qubit_ids:
                for vertex in shell.vertices:
                    plt.text(vertex.pos[0], vertex.pos[1], qubit_data.qubit[vertex.id])

            if show_face_ids and qubit_data is not None:

                def get_representative_point(points):
                    # Just uses the centroid for the moment
                    # This needs to be replaced for irregular
                    # polygons
                    xs = [point[0] for point in points]
                    zs = [point[1] for point in points]
                    return sum(xs) / len(xs), sum(zs) / len(zs)

                # Shapely version
                # from shapely.geometry import Polygon
                # def get_representative_point(points):
                #     poly = Polygon(points)
                #     # pylint: disable=no-member
                #     cx = poly.representative_point().x
                #     cy = poly.representative_point().y
                #     return cx, cy

                for face in shell.faces:
                    points = [vertex.pos for vertex in face.vertices]
                    if len(points) > 1:
                        if len(points) == 2:
                            mid_point = (points[0] + points[1]) / 2
                            # points.append(np.array([points[0][0], points[1][1]]))
                            plt.text(mid_point[0], mid_point[1], face.id)
                        else:
                            plt.text(*get_representative_point(points), face.id)

        plt.axis("equal")
