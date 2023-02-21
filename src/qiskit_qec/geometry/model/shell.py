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
"""Module for Shell"""

from typing import Callable, Dict, List, Optional, Tuple, Union

import logging

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit_qec.geometry.model.edge import Edge
from qiskit_qec.geometry.model.face import Face
from qiskit_qec.geometry.model.qubit_count import QubitCount
from qiskit_qec.geometry.model.qubit_data import QubitData
from qiskit_qec.geometry.model.shape_object import ShapeObject
from qiskit_qec.geometry.model.vertex import Vertex
from qiskit_qec.geometry.model.wireframe import WireFrame
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.operators.pauli import Pauli
from qiskit_qec.operators.pauli_list import PauliList


class Shell(ShapeObject):
    """`Shell` inherits from `ShapeObject`"""

    def __init__(self, faces: List[Face]) -> None:
        """Inits Shell.
        Each shell keeps track internally of several lists for
        each type of subcomponent (ex: faces, wireframes, edges, vertices)
        Args:
            faces (List[Face]): Faces that make up the `Shell` instance
        """
        super().__init__()

        self.faces = faces
        self.wireframes = []
        self.edges = []
        self.vertices = []

        for face in self.faces:
            self.wireframes += face.wireframes
            self.edges += face.edges
            self.vertices += face.vertices
            face.add_parent(self)

    def union(self, other_shell: "Shell"):
        # TODO control for overlapping subcomponents of other_shell and self
        """Add all subcomponents of other_shell to self.

        Must be disjoint at moment.

        Args:
            other_shell (Shell): A different shell
        """
        self.faces += other_shell.faces
        self.wireframes += other_shell.wireframes
        self.edges += other_shell.edges
        self.vertices += other_shell.vertices

    def delete_subtree(self, roots: ShapeObject):
        """Delete subcomponents of self"""
        pass

    def extract(
        self,
        is_inside: Dict[Vertex, bool],
        qubit_data: QubitData,
        qubit_count: QubitCount,
        levels: Optional[List[int]] = None,
        inside_levels: Optional[List[int]] = None,
        boundary_levels: Optional[List[int]] = None,
        exclude: Optional[Callable] = None,
        boundary_strategy: str = "combine",
        debug: str = False,
    ) -> Tuple["Shell", QubitData, QubitCount]:
        """Extracts the subshell defined by in_vertices and is_inside

        Args:
            is_inside (Dict[Vertex,bool]): _description_
            qubit_data (QubitData): _description_
            qubit_count (QubitCount): _description_
            levels (Optional[List[int]], optional): _description_. Defaults to [2,3,4]
            inside_levels (optional): Which operator weights to include for operators
                inside the cutter region. Default is `levels`
            boundary_levels (optional): Which operator weights to include for operators
                cut by the cutter region (on the boundary). Default is `levels`
            exclude (optional): method used to determine which boundary operators
                are to be used. Must be of the form

                def exclude(vertex_paths: List[List[Vertex]], qubit_data: QubitData) -> bool:
                    pass

                Will return a True value if the given operator (given by the vertex_paths)
                is to be excluded. Deflault always returns False. That is exlude nothing
                extra.
            boundary_strategy (str, optional): _description_. Defaults to "combine".
            debug (optional) : With debug set to True a range of intermediate details are printed. Used
            for debuging the addition of features to the method.

        Returns:
            Tuple[Shell, QubitData, QubitCount]: _description_
        """

        if debug:
            logger = logging.getLogger(__name__)

        if levels is None:
            levels = [2, 3, 4]

        if inside_levels is None:
            inside_levels = levels

        if boundary_levels is None:
            boundary_levels = levels

        if exclude is None:
            # pylint: disable=unused-argument
            def exclude(vertex_paths: List[List[Vertex]], qubit_data: QubitData) -> bool:
                return False

        def _other_vertex(parent: Edge, try_avoid: Vertex) -> Vertex:
            next_vertex = parent.vertices[0]
            if next_vertex == try_avoid:
                try:
                    return parent.vertices[1]
                except KeyError:
                    pass
            return next_vertex

        # Set up new data store
        new_qubit_data = QubitData()

        if not isinstance(levels, list):
            raise QiskitError(f"levels must be None of a list on integers: {levels}")

        def _find_restricted_path(s_vertex: Vertex) -> List[Vertex]:
            vertex_path = [s_vertex]
            edges_of_s_vertex = s_vertex.parents
            if len(edges_of_s_vertex) == 1:
                if len(edges_of_s_vertex[0].vertices) == 1:
                    return vertex_path
            v_start = s_vertex
            e_sn = edges_of_s_vertex[0]
            n_vertex = _other_vertex(e_sn, try_avoid=s_vertex)

            while True:
                if is_inside[n_vertex]:
                    vertex_path.append(n_vertex)

                if n_vertex == v_start:
                    # Path is a cycle returing to start
                    return vertex_path

                edges_of_n_vertex = n_vertex.parents
                if len(edges_of_n_vertex) == 1:
                    # End of path reached (not a loop)
                    return vertex_path

                e_nr = edges_of_n_vertex[0]
                if e_nr == e_sn:
                    e_nr = edges_of_n_vertex[1]

                r_vertex = _other_vertex(e_nr, try_avoid=n_vertex)

                if r_vertex == n_vertex:
                    if is_inside[r_vertex]:
                        # self loop at end of path (should not really happend but just in case)
                        vertex_path.append(r_vertex)
                    return vertex_path

                s_vertex = n_vertex
                n_vertex = r_vertex
                e_sn = e_nr

        def _weight_len(path: List) -> int:
            length = len(path)
            if path[0] == path[-1] and length > 1:
                return length - 1
            return length

        face_list = []

        for face in self.faces:
            face_inside = False
            edges = []
            vertex_paths = []
            rd_vertices = {item for item in face.vertices if is_inside[item]}

            if len(rd_vertices) == len(face.vertices):
                face_inside = True

            if debug:
                logger.debug(
                    "\n\nrd_vertices start = %s for face: %s : %s",
                    rd_vertices,
                    face.id,
                    face.vertices,
                )

            while len(rd_vertices) > 0:
                s_vertex = rd_vertices.pop()
                path = _find_restricted_path(s_vertex)
                if debug:
                    logger.debug("Found path= %s", path)
                rd_vertices.difference_update(path)
                if debug:
                    logger.debug("Remaining rd_vertices = %s", rd_vertices)
                vertex_paths.append(path)

            if debug:
                logger.debug("Final Vertex paths for face=%s:\n%s", face.id, vertex_paths)

            # Process the set of vertex paths to create edges
            if boundary_strategy == "combine":
                weights = [_weight_len(path) for path in vertex_paths]
                weight = sum(weights)
                if debug:
                    logger.debug("Weight=%s", weight)
                    logger.debug("vertex_paths=%s", vertex_paths)
                    logger.debug("exclude=%s", exclude)

                if (
                    (face_inside and weight in inside_levels)
                    or (not face_inside and weight in boundary_levels)
                ) and not exclude(vertex_paths, qubit_data):
                    # if weight in levels and not exclude(vertex_paths, qubit_data):
                    for path in vertex_paths:
                        if debug:
                            logger.debug("Working on v/e for path %s", path)

                        # Create the new vertices
                        if len(path) > 1:
                            if path[0] == path[-1]:
                                # have a loop
                                if debug:
                                    logger.debug("Have a loop!")
                                new_path = [vertex.shallowcopy() for vertex in path[:-1]]
                                short_path = path[:-1]
                                if debug:
                                    logger.debug("short_path is ", short_path)
                            else:
                                # Have non loop path of edges
                                if debug:
                                    logger.debug("Have a non loop path of edges")
                                new_path = [vertex.shallowcopy() for vertex in path]
                                short_path = path
                                if debug:
                                    logger.debug("short_path is ", short_path)
                        else:
                            new_path = [vertex.shallowcopy() for vertex in path]
                            short_path = path
                            if debug:
                                logger.debug("short_path is ", short_path)

                        # Copy over the qubit data for the vertices
                        for new_vertex, old_vertex in zip(new_path, short_path):
                            # Set the qubit for new vertex to the same qubit as old vertex
                            new_qubit_data.qubit[new_vertex.id] = qubit_data.qubit[old_vertex.id]

                            # Record the use of the qubit
                            qubit_count.increment_qubit(new_qubit_data.qubit[new_vertex.id])

                            # Set the vertex operator
                            new_qubit_data.operator[new_vertex.id] = qubit_data.operator[
                                old_vertex.id
                            ]

                        # Create the new edges
                        if len(path) == 1:
                            if debug:
                                logger.debug("Path is length 1")
                            edges.append(Edge([new_path[0]]))
                            if debug:
                                logger.debug("Edges: %s", edges[-1])
                        else:
                            if debug:
                                logger.debug("Path is length > 1")
                            for index, vertex in enumerate(new_path[:-1]):
                                edges.append(Edge([vertex, new_path[index + 1]]))
                            if len(path) == len(short_path) + 1:
                                if debug:
                                    logger.debug(
                                        "Have a loop with len(path) == len(short_path) + 1: \
                                            - adding extra edge"
                                    )
                                if len(short_path) > 2:
                                    # Loop
                                    edges.append(Edge([new_path[-1], new_path[0]]))

            else:
                raise QiskitError(f"Unknown boundary strategy {boundary_strategy}")

            # Create the wireframe and face
            if len(edges) > 0:
                if debug:
                    logger.debug("Create Wireframe and face with edges:")
                    for edge in edges:
                        print(edge)
                new_wf = WireFrame(edges)
                if debug:
                    logger.debug("New wireframe created with id=%s", new_wf.id)
                new_face = Face([new_wf])
                if debug:
                    logger.debug("New face created with id=%s", new_face.id)

                # Copy over face colors (if there any)
                try:
                    new_qubit_data.face_colors[new_face.id] = qubit_data.face_colors[face.id]
                except KeyError:
                    pass

                face_list.append(new_face)

                # Set face data (Keep orientation even if face is strict subface)
                try:
                    new_qubit_data.orientation[new_face.id] = qubit_data.orientation[face.id]
                except AttributeError:
                    pass

        # Create the new shell
        new_shell = Shell(face_list)

        # Update the qubit_counts
        for vertex in self.vertices:
            qubit_count.decrement_qubit(qubit_data.qubit[vertex.id])

        return new_shell, new_qubit_data, qubit_count

    @staticmethod
    def shell2symplectic(
        shell: "Shell",
        qubit_data: QubitData,
        qubit_count: QubitCount,
        from_index: Optional[Dict[int, int]] = None,
        from_qubit: Optional[Dict[int, int]] = None,
    ) -> Tuple[PauliList, QubitData]:
        """Converts a shell into a symplectic matrix given the qubit data and counts

        Args:
            shell (Shell): _description_
            qubit_data (QubitData): _description_
            qubit_count (QubitCount): _description_
            from_index (Optional[Dict[int,int]], optional): _description_. Defaults to None.
            from_qubit (Optional[Dict[int,int]], optional): _description_. Defaults to None.

        Returns:
            PauliList: _description_
        """
        if from_index is None:
            from_index = [
                qubit_id for qubit_id, count in qubit_count.qubits_count.items() if count != 0
            ]
            from_index = dict(enumerate(from_index))
            # from_index = {index: qubit_id for index, qubit_id in enumerate(from_index)}
        if from_qubit is None:
            from_qubit = {qubit_id: index for index, qubit_id in from_index.items()}

        # TODO: Make this a copy that is changed and then returned instead of
        # doing in place and returning
        qubit_data.qubit_to_index = from_qubit
        qubit_data.index_to_qubit = from_index
        qubit_data.index = {
            vertex.id: from_qubit[qubit_data.qubit[vertex.id]] for vertex in shell.vertices
        }

        pauli_str_list = []
        for face in shell.faces:
            num_stacked_operators = len(qubit_data.operator[face.vertices[0].id])
            for i in range(num_stacked_operators):
                pauli_str = ""
                for vertex in face.vertices:
                    pauli_str += str(qubit_data.operator[vertex.id][i]) + str(
                        qubit_data.index[vertex.id]
                    )
                pauli_str_list.append(pauli_str)
        return PauliList(pauli_str_list), qubit_data

    # pylint: disable=unused-argument
    def draw(
        self,
        qubit_data: Optional[QubitData] = None,
        show_index: bool = False,
        show_face_ids: bool = False,
        show_axis: bool = True,
        face_colors: bool = False,
        show_qubits: bool = True,
        figsize: Optional[Tuple[float, float]] = None,
        point_size: Optional[int] = 50,
        xcolor: str = "tomato",
        zcolor: str = "yellowgreen",
        ycolor: str = "steelblue",
        **kwargs,
    ) -> None:
        """Draws the shell

        Coloring only works at the moment for CSS codes.

        Args:
            qubit_data (Optional[QubitData], optional): _description_. Defaults to None.
            show_index (bool, optional): _description_. Defaults to False.
            show_face_ids (bool, optional): _description_. Defaults to False.
            show_axis (bool, optional): _description_. Defaults to True.
            face_colors (bool, optional): _description_. Defaults to False.
            show_qubits (bool, optional): _description_. Defaults to True.
            figsize (Optional[Tuple[float, float]], optional): _description_. Defaults to None.
            point_size (optional): size of points for qubits or vertices. Default is 50
            xcolor (str, optional): _description_. Defaults to "tomato".
            zcolor (str, optional): _description_. Defaults to "yellowgreen".
            ycolor (str, optional): _description_. Defaults to "steelblue".
            **kwargs: other options

        Returns:
            _type_: _description_
        """
        if figsize is None:
            figsize = (10, 10)
        stab = []
        for face in self.faces:
            verts = []
            for vertex in face.vertices:
                verts.append(list(vertex.pos))
            # The following only works for CSS codes. If there are two
            # operators for a given face then only on of the operators
            # will be used.
            pauli = qubit_data.operator[face.vertices[0].id][0]
            if face_colors:
                stab.append([verts, qubit_data.face_colors[face.id]])
            else:
                if pauli == Pauli("X"):
                    stab.append([verts, xcolor])
                elif pauli == Pauli("Z"):
                    stab.append([verts, zcolor])
                else:
                    stab.append([verts, ycolor])

        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)

        for stabilizer in stab:
            if len(stabilizer[0]) == 2:
                ax.add_patch(plt.Polygon(stabilizer[0], color=stabilizer[1], linewidth=7))
            else:
                ax.add_patch(plt.Polygon(stabilizer[0], facecolor=stabilizer[1], edgecolor="black"))

        if show_index and qubit_data is not None:
            offset = min(figsize) * 0.01
            # Change to run over quibits and not vertices as not to
            # print multple times on the same spot
            for vertex in self.vertices:
                plt.text(
                    vertex.pos[0] + offset, vertex.pos[1] + offset, qubit_data.index[vertex.id]
                )

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

            for face in self.faces:
                points = [vertex.pos for vertex in face.vertices]
                if len(points) > 1:
                    if len(points) == 2:
                        mid_point = (points[0] + points[1]) / 2
                        # points.append(np.array([points[0][0], points[1][1]]))
                        plt.text(mid_point[0], mid_point[1], face.id)
                    else:
                        plt.text(*get_representative_point(points), face.id)

        if not show_axis:
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

        if show_qubits:
            qubits = [list(vertex.pos) for vertex in self.vertices]
            qubits_x = [v[0] for v in qubits]
            qubits_y = [v[1] for v in qubits]

            plt.scatter(qubits_x, qubits_y, s=point_size, label="qubits", color="black")

        plt.axis("equal")

    def rotate2d(
        self, angle: Optional[float] = 90, about_point: Optional[Tuple[float]] = None
    ) -> None:
        """Inplace rotation of shell

        Args:
            angle (Optional[float], optional): _description_. Defaults to 90.
            about_point (Optional[Tuple[float]], optional): _description_. Defaults to None.
        """
        if about_point is None:
            about_point = np.array((0, 0))
        for vertex in self.vertices:
            vertex.pos = Plane.rotate(angle, vertex.pos - about_point) + about_point

    def scale(self, scale: float = 1) -> None:
        """Scale the shell by scale

        Args:
            scale (float, optional): Scaling value. Defaults to 1.
        """
        for vertex in self.vertices:
            vertex.pos = scale * vertex.pos

    def integer_snap(self) -> None:
        """Snap vertex components to integers"""
        for vertex in self.vertices:
            vertex.pos = np.rint(vertex.pos)

    def shift(self, vector: Union[None, np.ndarray, Tuple[int]] = None) -> None:
        """Shift/traslate the shell by vector

        Args:
            vector (optional): Shift/traslate the shell by vector. Defaults to None.
        """
        if vector is None:
            vector = (0, 0)
        vector = np.asarray(vector)
        for vertex in self.vertices:
            vertex.pos = vertex.pos + vector
