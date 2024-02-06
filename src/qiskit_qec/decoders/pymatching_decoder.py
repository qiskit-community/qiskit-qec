# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# pylint: disable=invalid-name, disable=no-name-in-module, disable=no-member

"""PyMatching"""
from typing import List, Union
from pymatching import Matching
from qiskit_qec.decoders.decoding_graph import DecodingGraph
from qiskit_qec.utils import DecodingGraphNode, DecodingGraphEdge


class PyMatchingDecoder:
    """
    Matching decoder using PyMatching.
    """

    def __init__(
        self,
        code_circuit,
        decoding_graph: DecodingGraph = None,
    ):
        """Setting up the matching object"""
        self.code = code_circuit
        if decoding_graph:
            self.decoding_graph = decoding_graph
        else:
            self.decoding_graph = DecodingGraph(self.code)
        if self.decoding_graph.logical_nodes:
            self.graph = self.decoding_graph.get_edge_graph()
        else:
            self.graph = self.decoding_graph.graph
        self.matcher = self._matching()
        self.indexer = None
        super().__init__()

    def _matching(self) -> Matching:
        return Matching(self.graph)

    def logical_flips(self, syndrome: Union[List[DecodingGraphNode], List[int]]) -> List[int]:
        """
        Args:
            syndrome: Either a list of DecodingGraphNode objects returnes by string2nodes,
            or a list of binaries indicating which node is highlighted, e.g.,
            the output of a stim detector sampler
        Returns: list of binaries indicating which logical is flipped
        """
        syndrome_as_nodes = True
        for elem in syndrome:
            syndrome_as_nodes = syndrome_as_nodes and isinstance(elem, DecodingGraphNode)
        if syndrome_as_nodes:
            syndrome = self.nodes_to_detections(syndrome)
        return self.matcher.decode(syndrome)

    def process(self, string: str) -> List[int]:
        """
        Converts qiskit counts string into a list of flipped logicals
        Args: counts string
        Returns: list of corrected logicals (0 or 1)
        """
        nodes = self.code.string2nodes(string)
        raw_logicals = self.code.string2raw_logicals(string)

        logical_flips = self.logical_flips(nodes)

        corrected_logicals = [
            (int(raw) + flip) % 2 for raw, flip in zip(raw_logicals, logical_flips)
        ]

        return corrected_logicals

    def matched_edges(
        self, syndrome: Union[List[DecodingGraphNode], List[int]]
    ) -> List[DecodingGraphEdge]:
        """
        Args:
            syndrome: Either a list of DecodingGraphNode objects returnes by string2nodes,
            or a list of binaries indicating which node is highlighted.
        Returns: list of DecodingGraphEdge-s included in the matching
        """
        if isinstance(syndrome[0], DecodingGraphNode):
            syndrome = self.nodes_to_detections(syndrome)
        edge_dets = list(self.graph.edge_list())
        edges = self.graph.edges()
        matched_det_pairs = self.matcher.decode_to_edges_array(syndrome)
        det_pairs = []
        for pair in matched_det_pairs:
            if pair[1] == -1:
                pair[-1] = pair[-1] + len(self.graph.nodes())
            pair.sort()
            det_pairs.append(tuple(pair))
        mached_edges = [edges[edge_dets.index(det_pair)] for det_pair in det_pairs]
        return mached_edges

    def nodes_to_detections(self, syndrome_nodes: List[DecodingGraphNode]) -> List[int]:
        """Converts nodes to detector indices to be used by pymatching.Matching.decode"""
        graph_nodes = self.graph.nodes()
        detections = [0] * len(graph_nodes)
        for i, node in enumerate(graph_nodes):
            if node in syndrome_nodes:
                detections[i] = 1
        return detections
