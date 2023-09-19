"""rustworkx matching object."""

import logging
from copy import deepcopy
from typing import Dict, List, Set, Tuple

import rustworkx as rx
from qiskit_qec.decoders.base_matcher import BaseMatcher
from qiskit_qec.utils import DecodingGraphEdge


class RustworkxMatcher(BaseMatcher):
    """Matching subroutines using rustworkx.

    The input rustworkx graph is expected to have decoding_graph.Node as the type of the node payload
    and decoding_graph.Edge as the type of the edge payload.

    Additionally the edges are expected to have the following properties:
        - edge.properties["measurement_error"] (bool): Whether or not the error
            corresponds to a measurement error.

    The annotated graph will also have "highlighted" properties on edges and vertices.
    """

    def __init__(self, annotate: bool = False):
        """Create the matcher."""
        self.length = {}
        self.path = {}
        self.annotate = annotate
        self.annotated_graph = None
        super().__init__()

    def preprocess(self, graph: rx.PyGraph):
        """Compute shortest paths between vertex pairs in decoding graph.

        Updates sets self.length and self.path.
        """

        # edge_cost_fn = lambda edge: edge["weight"]
        def edge_cost_fn(edge: DecodingGraphEdge):
            return edge.weight

        length = rx.all_pairs_dijkstra_path_lengths(graph, edge_cost_fn)
        self.length = {s: dict(length[s]) for s in length}
        path = rx.all_pairs_dijkstra_shortest_paths(graph, edge_cost_fn)
        self.path = {s: {t: list(path[s][t]) for t in path[s]} for s in path}

    def find_errors(
        self,
        graph: rx.PyGraph,
        idxmap: Dict[Tuple[int, List[int]], int],
        highlighted: List[Tuple[int, Tuple[int]]],
    ) -> Tuple[Set[int], Set[Tuple[int, Tuple[int]]]]:
        """Process a set of highlighted vertices and return error locations.

        Be sure to have called recompute_paths if needed.
        """
        matching = self._compute_matching(idxmap, highlighted)
        logging.info("process: matching = %s", matching)
        qubit_errors, measurement_errors = self._compute_error_correction(
            graph, idxmap, matching, highlighted
        )
        logging.info("process: qubit_errors = %s", qubit_errors)
        logging.debug("process: measurement_errors = %s", measurement_errors)
        return qubit_errors, measurement_errors

    def _compute_matching(
        self,
        idxmap: Dict[Tuple[int, List[int]], int],
        highlighted: List[Tuple[int, Tuple[int]]],
    ) -> Set[Tuple[int, int]]:
        """Compute a min. weight perfect matching of highlighted vertices.

        highlighted is a list of highlighted vertices given as tuples
        (t, qubit_set).
        Return the matching.
        """
        gm = rx.PyGraph(multigraph=False)  # matching graph
        idx = 0  # vertex index in matching graph
        midxmap = {}  # map from (t, qubit_tuple) to vertex index
        for v in highlighted:
            gm.add_node({"dvertex": v})
            midxmap[v] = idx
            idx += 1
        for i, high_i in enumerate(highlighted):
            for j in range(i + 1, len(highlighted)):
                vi = midxmap[high_i]
                vj = midxmap[highlighted[j]]
                vip = idxmap[high_i]
                vjp = idxmap[highlighted[j]]
                gm.add_edge(vi, vj, {"weight": -self.length[vip][vjp]})

        def weight_fn(edge):
            return int(edge["weight"])

        matching = rx.max_weight_matching(gm, max_cardinality=True, weight_fn=weight_fn)
        return matching

    @staticmethod
    def _error_chain_from_vertex_path(
        graph: rx.PyGraph, vertex_path: List[int]
    ) -> Tuple[Set[int], Set[Tuple[int, Tuple[int]]]]:
        """Return a chain of qubit and measurement errors from a vertex path.

        Examine the edges along the path to extract the error chain.
        Store error chains as sets and merge using symmetric difference.
        The vertex_path is a list of rustworkx node indices.
        """
        qubit_errors = set([])
        measurement_errors = set([])
        logging.debug("_error_chain_from_vertex_path %s", vertex_path)
        for i in range(len(vertex_path) - 1):
            v0 = vertex_path[i]
            v1 = vertex_path[i + 1]
            if graph.get_edge_data(v0, v1).properties["measurement_error"] == 1:
                measurement_errors ^= set(
                    [(graph.nodes()[v0].time, tuple(graph.nodes()[v0].qubits))]
                )
            qubit_errors ^= set(graph.get_edge_data(v0, v1).qubits)
            logging.debug(
                "_error_chain_for_vertex_path q = %s, m = %s",
                qubit_errors,
                measurement_errors,
            )
        return qubit_errors, measurement_errors

    def _compute_error_correction(
        self,
        graph: rx.PyGraph,
        idxmap: Dict[Tuple[int, List[int]], int],
        matching: Set[Tuple[int, int]],
        highlighted: List[Tuple[int, Tuple[int]]],
    ) -> Tuple[Set[int], Set[Tuple[int, Tuple[int]]]]:
        """Compute the qubit and measurement corrections.

        graph : the decoding graph
        idxmap : maps (t, qubit_idx) to vertex index
        matching : perfect matching computed by _compute_matching
        highlighted : list of highlighted vertices

        Returns a tuple of sets, (qubit_errors, measurement_errors) where
        qubit_errors contains the indices of qubits with errors and
        measurement_errors contains tuples (t, qubit_set) indicating
        failed measurements.
        """
        used_paths = []
        qubit_errors = set([])
        measurement_errors = set([])
        for p in matching:
            v0 = idxmap[highlighted[p[0]]]
            v1 = idxmap[highlighted[p[1]]]
            # Use the shortest paths between the matched vertices to
            # identify all of the qubits in the error chains
            path = self.path[v0][v1]
            q, m = self._error_chain_from_vertex_path(graph, path)
            # Add the error chains modulo two to get the total correction
            # (uses set symmetric difference)
            qubit_errors ^= q
            measurement_errors ^= m
            used_paths.append(path)
        if self.annotate:
            self.annotated_graph = self._make_annotated_graph(graph, used_paths)
        return qubit_errors, measurement_errors

    @staticmethod
    def _make_annotated_graph(gin: rx.PyGraph, paths: List[List[int]]) -> rx.PyGraph:
        """Highlight the vertex paths and return annotated graph.

        gin : decoding graph
        paths : list of vertex paths, each given as a list of
        vertex indices in the decoding graph.
        """
        graph = deepcopy(gin)
        for path in paths:
            # Highlight the endpoints of the path
            for i in [0, -1]:
                graph.nodes()[path[i]].properties["highlighted"] = True
            # Highlight the edges along the path
            for i in range(len(path) - 1):
                try:
                    idx = list(graph.edge_list()).index((path[i], path[i + 1]))
                except ValueError:
                    idx = list(graph.edge_list()).index((path[i + 1], path[i]))
                edge = graph.edges()[idx]
                edge.properties["highlighted"] = True
        return graph
