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

# pylint: disable=invalid-name

"""
Graph used as the basis of decoders.
"""
from dataclasses import dataclass, field
import itertools
from typing import Any, Dict, List, Tuple, Optional


class DecodingGraphNode:
    def __init__(self, qubits: List[int], index: int, is_boundary=False, time=None) -> None:
        if not is_boundary and time == None:
            raise QiskitQECError(
                "DecodingGraph node must either have a time or be a boundary node."
            )

        self.is_boundary: bool = is_boundary
        self.time: Optional[int] = time if not is_boundary else None
        self.qubits: List[int] = qubits
        self.index: int = index
        # TODO: Should code/decoder specific properties be accounted for when comparing nodes
        self.properties: Dict[str, Any] = dict()

    def __eq__(self, rhs):
        if not isinstance(rhs, DecodingGraphNode):
            return NotImplemented

        result = (
            self.index == rhs.index
            and set(self.qubits) == set(rhs.qubits)
            and self.is_boundary == rhs.is_boundary
        )
        if not self.is_boundary:
            result = result and self.time == rhs.time
        return result

    def __hash__(self) -> int:
        return hash(repr(self))

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value


@dataclass
class DecodingGraphEdge:
    qubits: List[int]
    weight: float
    # TODO: Should code/decoder specific properties be accounted for when comparing edges
    properties: Dict[str, Any] = field(default_factory=dict)

    def __eq__(self, rhs) -> bool:
        if not isinstance(rhs, DecodingGraphNode):
            return NotImplemented

        return set(self.qubits) == set(rhs.qubits) and self.weight == rhs.weight

    def __hash__(self) -> int:
        return hash(repr(self))

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value
