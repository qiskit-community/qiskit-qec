"""PyMatching matching object."""

from typing import List, Tuple, Dict, Set
import logging

import rustworkx as rx
from pymatching import Matching

from qiskit_qec.exceptions import QiskitQECError
from qiskit_qec.utils.indexer import Indexer
from qiskit_qec.decoders.base_matcher import BaseMatcher
from qiskit_qec.decoders.temp_graph_util import ret2net


class PyMatchingMatcher(BaseMatcher):
    """Matching subroutines using PyMatching.

    The input rustworkx graph is expected to have the following properties:
    edge["weight"] : real edge weight
    edge["qubits"] : list of qubit ids associated to edge
    vertex["is_boundary"] : bool, true if boundary node
    """

    def __init__(self):
        """Create the matcher."""
        self.pymatching = None
        self.indexer = None
        super().__init__()

    def preprocess(self, graph: rx.PyGraph):
        """Create the pymatching object.
        Add qubit_id properties to the graph.
        """
        self.indexer = Indexer()
        nxgraph = ret2net(graph)
        for edge in nxgraph.edges(data=True):
            if edge[2]["qubits"]:
                qset = set()
                for q in edge[2]["qubits"]:
                    qset.add(self.indexer[q])
                edge[2]["qubit_id"] = tuple(qset)[0] if len(qset) == 1 else tuple(qset)
            else:
                edge[2]["qubit_id"] = -1
        self.pymatching = Matching(nxgraph)

    def find_errors(
        self,
        graph: rx.PyGraph,
        idxmap: Dict[Tuple[int, List[int]], int],
        highlighted: List[Tuple[int, Tuple[int]]],
    ) -> Tuple[Set[int], Set[Tuple[int, Tuple[int]]]]:
        """Process a set of highlighted vertices and return error locations."""
        syndrome = [0] * len(idxmap)
        for vertex in highlighted:
            syndrome[idxmap[vertex]] = 1
        try:
            correction = self.pymatching.decode(syndrome)
        except AttributeError as attrib_error:
            raise QiskitQECError("Did you call preprocess?") from attrib_error
        qubit_errors = []
        for i, corr in enumerate(correction):
            if corr == 1:
                qubit_errors.append(self.indexer.rlookup(i))
        logging.info("qubit_errors = %s", qubit_errors)
        return set(qubit_errors), set()
