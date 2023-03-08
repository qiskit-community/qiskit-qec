"""Tests for the rustworkx matcher subroutines."""
import unittest

from typing import Dict, Tuple
import rustworkx as rx
from qiskit_qec.decoders.rustworkx_matcher import RustworkxMatcher
from qiskit_qec.utils import DecodingGraphNode, DecodingGraphEdge


class TestRustworkxMatcher(unittest.TestCase):
    """Tests for the rustworkx matcher subroutines."""

    @staticmethod
    def make_test_graph() -> Tuple[rx.PyGraph, Dict[Tuple[int, Tuple[int]], int]]:
        """Make a basic decoding graph.

        4 -- 0 -- 1 -- 2 -- 3 -- (4)
        """
        graph = rx.PyGraph(multigraph=False)
        idxmap = {}
        basic_config = [[0, 1], [1, 2], [2, 3], [3, 4]]
        for i, q in enumerate(basic_config):
            node = DecodingGraphNode(time=0, qubits=q, index=i)
            node.properties["highlighted"] = False
            graph.add_node(node)
            idxmap[(0, tuple(q))] = i
        node = DecodingGraphNode(time=0, qubits=[], index=len(basic_config) + 1)
        node.properties["highlighted"] = False
        graph.add_node(node)
        idxmap[(0, tuple([]))] = 4
        for dat in [[[0], 0, 4], [[1], 0, 1], [[2], 1, 2], [[3], 2, 3], [[4], 3, 4]]:
            edge = DecodingGraphEdge(qubits=dat[0], weight=1)
            edge.properties["measurement_error"] = False
            edge.properties["highlighted"] = False
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
        self.assertEqual(self.rxm.annotated_graph[0].properties["highlighted"], True)
        self.assertEqual(self.rxm.annotated_graph[1].properties["highlighted"], True)
        self.assertEqual(self.rxm.annotated_graph[2].properties["highlighted"], False)
        self.assertEqual(self.rxm.annotated_graph[3].properties["highlighted"], True)
        self.assertEqual(self.rxm.annotated_graph[4].properties["highlighted"], True)
        eim = self.rxm.annotated_graph.edge_index_map()
        self.assertEqual(eim[0][2].properties["highlighted"], False)
        self.assertEqual(eim[1][2].properties["highlighted"], True)
        self.assertEqual(eim[2][2].properties["highlighted"], False)
        self.assertEqual(eim[3][2].properties["highlighted"], False)
        self.assertEqual(eim[4][2].properties["highlighted"], True)


if __name__ == "__main__":
    unittest.main()
