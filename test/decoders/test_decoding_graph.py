# This code is part of Qiskit.
#
# (C) Copyright IBM 2019-2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Tests for DecodingGraph."""

import unittest

import rustworkx as rx

from qiskit_qec.circuits.repetition_code import RepetitionCodeCircuit as RepetitionCode
from qiskit_qec.decoders.decoding_graph import DecodingGraph
from qiskit_qec.utils import DecodingGraphEdge, DecodingGraphNode


def _rep_code(d=3, T=1):
    """Return a repetition code circuit for testing."""
    return RepetitionCode(d, T)


class TestDecodingGraphConstruction(unittest.TestCase):
    """Test DecodingGraph construction from a RepetitionCode circuit."""

    def setUp(self):
        self.code = _rep_code(d=3, T=1)
        self.dg = DecodingGraph(self.code)

    def test_returns_decoding_graph_instance(self):
        """DecodingGraph builds without error."""
        self.assertIsInstance(self.dg, DecodingGraph)

    def test_graph_is_pygraph(self):
        """The internal graph is a rustworkx PyGraph."""
        self.assertIsInstance(self.dg.graph, rx.PyGraph)

    def test_graph_has_nodes(self):
        """A non-trivial code produces at least one node in the graph."""
        self.assertGreater(len(self.dg.graph.nodes()), 0)

    def test_graph_has_edges(self):
        """A non-trivial code produces at least one edge in the graph."""
        self.assertGreater(len(self.dg.graph.edges()), 0)

    def test_logical_nodes_list_populated(self):
        """logical_nodes attribute is a list (may be empty for rep-code helper graph)."""
        self.assertIsInstance(self.dg.logical_nodes, list)

    def test_hyperedges_is_list(self):
        """hyperedges attribute is always a list."""
        self.assertIsInstance(self.dg.hyperedges, list)

    def test_brute_mode_builds(self):
        """brute=True also builds a valid DecodingGraph."""
        dg_brute = DecodingGraph(self.code, brute=True)
        self.assertIsInstance(dg_brute, DecodingGraph)
        self.assertIsInstance(dg_brute.graph, rx.PyGraph)

    def test_d5_graph_has_more_nodes_than_d3(self):
        """A d=5 code graph is strictly larger than a d=3 code graph."""
        code5 = _rep_code(d=5, T=1)
        dg5 = DecodingGraph(code5)
        self.assertGreater(
            len(dg5.graph.nodes()),
            len(self.dg.graph.nodes()),
        )


class TestDecodingGraphNodeIndex(unittest.TestCase):
    """Test node_index and edge_in_graph helper methods."""

    def setUp(self):
        self.code = _rep_code(d=3, T=1)
        self.dg = DecodingGraph(self.code)

    def test_node_index_returns_int_for_every_node(self):
        """node_index returns an int for every node already in the graph."""
        for node in self.dg.graph.nodes():
            idx = self.dg.node_index(node)
            self.assertIsInstance(idx, int)

    def test_edge_in_graph_true_for_existing_edges(self):
        """edge_in_graph is True for every edge in the graph's edge list."""
        for edge in self.dg.graph.edge_list():
            self.assertTrue(self.dg.edge_in_graph(edge))

    def test_edge_in_graph_false_for_non_existing_edge(self):
        """edge_in_graph is False for an edge that was never added."""
        self.assertFalse(self.dg.edge_in_graph((9999, 9998)))

    def test_update_attributes_does_not_raise(self):
        """Calling update_attributes explicitly does not raise."""
        self.dg.update_attributes()  # should be idempotent


class TestDecodingGraphNode(unittest.TestCase):
    """Test DecodingGraphNode construction and equality."""

    def test_boundary_node(self):
        """A boundary node is created without a time value."""
        node = DecodingGraphNode(index=0, is_boundary=True)
        self.assertTrue(node.is_boundary)
        self.assertIsNone(node.time)

    def test_logical_node(self):
        """A logical node is created without a time value."""
        node = DecodingGraphNode(index=0, is_logical=True)
        self.assertTrue(node.is_logical)
        self.assertIsNone(node.time)

    def test_bulk_node_requires_time(self):
        """A non-boundary, non-logical node must carry a time value."""
        from qiskit_qec.exceptions import QiskitQECError

        with self.assertRaises(QiskitQECError):
            DecodingGraphNode(index=0)

    def test_bulk_node_with_time(self):
        """A bulk node with time set is constructed correctly."""
        node = DecodingGraphNode(index=1, qubits=[0, 1], time=0)
        self.assertEqual(node.time, 0)
        self.assertEqual(node.qubits, [0, 1])

    def test_default_qubits_empty(self):
        """qubits defaults to an empty list when not specified."""
        node = DecodingGraphNode(index=0, is_boundary=True)
        self.assertEqual(node.qubits, [])

    def test_properties_dict_empty_by_default(self):
        """properties is an empty dict on a fresh node."""
        node = DecodingGraphNode(index=0, is_boundary=True)
        self.assertEqual(node.properties, {})


class TestDecodingGraphEdge(unittest.TestCase):
    """Test DecodingGraphEdge construction."""

    def test_edge_stores_qubits_and_weight(self):
        """DecodingGraphEdge stores qubits and weight correctly."""
        edge = DecodingGraphEdge(qubits=[0, 1], weight=1)
        self.assertEqual(edge.qubits, [0, 1])
        self.assertEqual(edge.weight, 1)

    def test_edge_default_fault_ids_empty(self):
        """fault_ids defaults to an empty set."""
        edge = DecodingGraphEdge(qubits=[0], weight=1)
        self.assertEqual(edge.fault_ids, set())


class TestDecodingGraphEdgeGraph(unittest.TestCase):
    """Test get_edge_graph and get_node_graph conversion methods."""

    def setUp(self):
        self.code = _rep_code(d=3, T=1)
        self.dg = DecodingGraph(self.code)

    def test_get_edge_graph_returns_pygraph(self):
        """get_edge_graph returns a rustworkx PyGraph."""
        eg = self.dg.get_edge_graph()
        self.assertIsInstance(eg, rx.PyGraph)

    def test_get_node_graph_returns_pygraph(self):
        """get_node_graph returns a rustworkx PyGraph."""
        ng = self.dg.get_node_graph()
        self.assertIsInstance(ng, rx.PyGraph)

    def test_edge_graph_has_no_logical_nodes(self):
        """get_edge_graph converts logical nodes to boundary nodes."""
        eg = self.dg.get_edge_graph()
        for node in eg.nodes():
            self.assertFalse(node.is_logical)

    def test_graph_unchanged_after_get_edge_graph(self):
        """Calling get_edge_graph does not mutate the original graph."""
        node_count_before = len(self.dg.graph.nodes())
        self.dg.get_edge_graph()
        self.assertEqual(len(self.dg.graph.nodes()), node_count_before)


if __name__ == "__main__":
    unittest.main()
