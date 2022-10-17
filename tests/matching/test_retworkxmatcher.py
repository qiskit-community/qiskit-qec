"""Tests for the rustworkx matcher subroutines."""
import unittest

from typing import Dict, Tuple
import rustworkx as rx
from qiskit_qec.decoders.rustworkx_matcher import RustworkxMatcher


class TestRustworkxMatcher(unittest.TestCase):
    """Tests for the rustworkx matcher subroutines."""

    def make_test_graph(self) -> Tuple[rx.PyGraph, Dict[Tuple[int, Tuple[int]], int]]:
        """Make a basic decoding graph.

        4 -- 0 -- 1 -- 2 -- 3 -- (4)
        """
        graph = rx.PyGraph(multigraph=False)
        idxmap = {}
        for i, q in enumerate([[0, 1], [1, 2], [2, 3], [3, 4]]):
            node = {"time": 0, "qubits": q, "highlighted": False}
            graph.add_node(node)
            idxmap[(0, tuple(q))] = i
        node = {"time": 0, "qubits": [], "highlighted": False}
        graph.add_node(node)
        idxmap[(0, tuple([]))] = 4
        for dat in [[[0], 0, 4], [[1], 0, 1], [[2], 1, 2], [[3], 2, 3], [[4], 3, 4]]:
            edge = {"qubits": dat[0], "measurement_error": False, "weight": 1, "highlighted": False}
            graph.add_edge(dat[1], dat[2], edge)
        return graph, idxmap

    def setUp(self) -> None:
        self.rxm = RustworkxMatcher(annotate=True)

    def test_preprocess(self):
        """Test preprocessing example."""
        graph, _ = self.make_test_graph()
        self.rxm.preprocess(graph)
        self.assertEqual(self.rxm.length[0][1], 1)
        self.assertEqual(self.rxm.length[1][3], 2)
        self.assertEqual(self.rxm.length[2][4], 2)
        self.assertEqual(self.rxm.length[0][3], 2)
        self.assertEqual(self.rxm.path[0][2], [0, 1, 2])
        self.assertEqual(self.rxm.path[3][0], [3, 4, 0])

    def test_match(self):
        """Test matching example."""
        graph, idxmap = self.make_test_graph()
        self.rxm.preprocess(graph)
        highlighted = [(0, (0, 1)), (0, (1, 2)), (0, (3, 4)), (0, ())]  # must be even
        qubit_errors, measurement_errors = self.rxm.find_errors(graph, idxmap, highlighted)
        self.assertEqual(qubit_errors, set([1, 4]))
        self.assertEqual(measurement_errors, set())

    def test_annotate(self):
        """Test the annotated graph."""
        graph, idxmap = self.make_test_graph()
        self.rxm.preprocess(graph)
        highlighted = [(0, (0, 1)), (0, (1, 2)), (0, (3, 4)), (0, ())]  # must be even
        self.rxm.find_errors(graph, idxmap, highlighted)
        self.assertEqual(self.rxm.annotated_graph[0]["highlighted"], True)
        self.assertEqual(self.rxm.annotated_graph[1]["highlighted"], True)
        self.assertEqual(self.rxm.annotated_graph[2]["highlighted"], False)
        self.assertEqual(self.rxm.annotated_graph[3]["highlighted"], True)
        self.assertEqual(self.rxm.annotated_graph[4]["highlighted"], True)
        eim = self.rxm.annotated_graph.edge_index_map()
        self.assertEqual(eim[0][2]["highlighted"], False)
        self.assertEqual(eim[1][2]["highlighted"], True)
        self.assertEqual(eim[2][2]["highlighted"], False)
        self.assertEqual(eim[3][2]["highlighted"], False)
        self.assertEqual(eim[4][2]["highlighted"], True)


if __name__ == "__main__":
    unittest.main()
