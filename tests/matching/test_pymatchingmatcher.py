"""Tests for the rustworkx matcher subroutines."""
import unittest

from typing import Dict, Tuple
import rustworkx as rx
from qiskit_qec.decoders.pymatching_matcher import PyMatchingMatcher


class TestPyMatchingMatcher(unittest.TestCase):
    """Tests for the pymatching matcher subroutines."""

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
        node = {"time": 0, "qubits": [], "highlighted": False, "is_boundary": True}
        graph.add_node(node)
        idxmap[(0, tuple([]))] = 4
        for dat in [[[0], 0, 4], [[1], 0, 1], [[2], 1, 2], [[3], 2, 3], [[4], 3, 4]]:
            edge = {"qubits": dat[0], "measurement_error": False, "weight": 1, "highlighted": False}
            graph.add_edge(dat[1], dat[2], edge)
        return graph, idxmap

    def setUp(self) -> None:
        self.m = PyMatchingMatcher()

    def test_match(self):
        """Test matching example."""
        graph, idxmap = self.make_test_graph()
        self.m.preprocess(graph)
        highlighted = [(0, (0, 1)), (0, (1, 2)), (0, (3, 4)), (0, ())]  # must be even
        qubit_errors, measurement_errors = self.m.find_errors(graph, idxmap, highlighted)
        self.assertEqual(qubit_errors, set([1, 4]))
        self.assertEqual(measurement_errors, set())


if __name__ == "__main__":
    unittest.main()
