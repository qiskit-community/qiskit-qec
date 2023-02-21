"""Base matching object."""
from abc import ABC, abstractmethod
from typing import List, Tuple, Dict, Set
import rustworkx as rx


class BaseMatcher(ABC):
    """Matching subroutine base class."""

    def __init__(self):
        """Create the base matcher."""
        pass

    @abstractmethod
    def preprocess(self, graph: rx.PyGraph):
        """Do any preprocessing using the graph data."""
        raise NotImplementedError("Not implemented.")

    @abstractmethod
    def find_errors(
        self,
        graph: rx.PyGraph,
        idxmap: Dict[Tuple[int, List[int]], int],
        highlighted: List[Tuple[int, Tuple[int]]],
    ) -> Tuple[Set[int], Set[Tuple[int, Tuple[int]]]]:
        """Process a set of highlighted vertices and return error locations."""
        raise NotImplementedError("Not implemented.")
