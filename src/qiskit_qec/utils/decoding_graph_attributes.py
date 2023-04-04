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
from typing import Any, Dict, List, Optional

from qiskit_qec.exceptions import QiskitQECError


class DecodingGraphNode:
    """
    Class to describe DecodingGraph nodes.

    Attributes:
     - is_boundary (bool): whether or not the node is a boundary node.
     - time (int): what syndrome node the node corrsponds to. Doesn't
        need to be set if it's a boundary node.
     - qubits (List[int]): List of indices which are stabilized by
        this ancilla.
     - index (int): Unique index in measurement round.
     - properties (Dict[str, Any]): Decoder/code specific attributes.
        Are not considered when comparing nodes.
    """

    def __init__(self, index: int, qubits: List[int] = None, is_boundary=False, time=None) -> None:
        if not is_boundary and time is None:
            raise QiskitQECError(
                "DecodingGraph node must either have a time or be a boundary node."
            )

        self.is_boundary: bool = is_boundary
        self.time: Optional[int] = time if not is_boundary else None
        self.qubits: List[int] = qubits if qubits else []
        self.index: int = index
        self.properties: Dict[str, Any] = {}

    def __getitem__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        elif key in self.properties:
            return self.properties[key]
        else:
            raise QiskitQECError("'" + str(key) + "'" + " is not an an attribute or property of the node.")

    def get(self, key, default):
        return self.__getitem__(key)

    def __setitem__(self, key, value):
        if key in self.__dict__:
            self.__dict__[key] = value
        else:
            self.properties[key] = value
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

    def __str__ (self):
        return str(dict(self))

@dataclass
class DecodingGraphEdge:
    """
    Class to describe DecodingGraph edges.

    Attributes:
     - qubits (List[int]): List of indices of code qubits that correspond to this edge.
     - weight (float): Weight of the edge.
     - properties (Dict[str, Any]): Decoder/code specific attributes.
        Are not considered when comparing edges.
    """

    qubits: List[int]
    weight: float
    # TODO: Should code/decoder specific properties be accounted for when comparing edges
    properties: Dict[str, Any] = field(default_factory=dict)

    def __getitem__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        elif key in self.properties:
            return self.properties[key]
        else:
            raise QiskitQECError("'" + str(key) + "'" + " is not an an attribute or property of the edge.")

    def get(self, key, default):
        return self.__getitem__(key)

    def __setitem__(self, key, value):
        if key in self.__dict__:
            self.__dict__[key] = value
        else:
            self.properties[key] = value
            
    def __eq__(self, rhs) -> bool:
        if not isinstance(rhs, DecodingGraphNode):
            return NotImplemented

        return set(self.qubits) == set(rhs.qubits) and self.weight == rhs.weight

    def __hash__(self) -> int:
        return hash(repr(self))

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    def __str__ (self):
        return str(dict(self))
