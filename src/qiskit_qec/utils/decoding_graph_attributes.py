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
from typing import Union, Any, Dict, List, Set, Optional

from qiskit_qec.exceptions import QiskitQECError


class DecodingGraphNode:
    """
    Class to describe DecodingGraph nodes.

    Attributes:
     - is_boundary (bool): whether or not the node is a boundary node.
     - is_logical (bool): whether or not the node is a logical node.
     - time (int): what syndrome node the node corrsponds to. Doesn't
        need to be set if it's a boundary node.
     - qubits (List[int]): List of indices which are stabilized by
        this ancilla.
     - index (int): Unique index in measurement round.
     - properties (Dict[str, Any]): Decoder/code specific attributes.
        Are not considered when comparing nodes.
    """

    def __init__(
        self,
        index: int,
        qubits: List[int] = None,
        is_boundary=False,
        is_logical=False,
        time=None,
    ) -> None:
        if not is_boundary and not is_logical and time is None:
            raise QiskitQECError(
                "DecodingGraph node must either have a time or be a boundary or logical node."
            )

        self.is_boundary: bool = is_boundary
        self.is_logical: bool = is_logical
        self.time: Optional[int] = None if (is_boundary or is_logical) else time
        self.qubits: List[int] = qubits if qubits else []
        self.index: int = index
        self.properties: Dict[str, Any] = {}

    def __getitem__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        elif key in self.properties:
            return self.properties[key]
        else:
            return QiskitQECError(
                "'" + str(key) + "'" + " is not an an attribute or property of the node."
            )

    def get(self, key, default=None):
        """Return value for given key."""
        # pylint: disable=unnecessary-dunder-call
        output = self.__getitem__(key)
        if isinstance(output, QiskitQECError):
            if default:
                return default
            else:
                raise output
        else:
            return output

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
            and self.is_logical == rhs.is_logical
        )
        if not (self.is_boundary or self.is_logical):
            result = result and self.time == rhs.time
        return result

    def __hash__(self) -> int:
        return hash(repr(self))

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    def __repr__(self):
        return str(dict(self))


@dataclass
class DecodingGraphEdge:
    """
    Class to describe DecodingGraph edges.

    Attributes:
     - qubits (List[int]): List of indices of code qubits that correspond to this edge.
     - weight (float): Weight of the edge.
     - fault_ids fault_ids: Union[Set[int],List[int]]: In the style of pymatching.
     - properties (Dict[str, Any]): Decoder/code specific attributes.
        Are not considered when comparing edges.
    """

    qubits: List[int]
    weight: float
    fault_ids: Union[Set[int], List[int]] = field(default_factory=set)
    properties: Dict[str, Any] = field(default_factory=dict)

    def __getitem__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        elif key in self.properties:
            return self.properties[key]
        else:
            return QiskitQECError(
                "'" + str(key) + "'" + " is not an an attribute or property of the edge."
            )

    def get(self, key, default=None):
        """Return value for given key."""
        # pylint: disable=unnecessary-dunder-call
        value = self.__getitem__(key)
        if isinstance(value, QiskitQECError):
            if default is not None:
                return default
            else:
                raise value
        else:
            return value

    def __setitem__(self, key, value):
        if key in self.__dict__:
            self.__dict__[key] = value
        else:
            self.properties[key] = value

    def __eq__(self, rhs) -> bool:
        if not isinstance(rhs, DecodingGraphEdge):
            return NotImplemented

        return set(self.qubits) == set(rhs.qubits) and self.weight == rhs.weight

    def __hash__(self) -> int:
        return hash(repr(self))

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    def __repr__(self):
        return str(dict(self))


def _nodes2cpp(nodes):
    """
    Convert a list of nodes to the form required by C++ functions.
    """
    # nodes are a tuple with (q0, q1,t, extra)
    # extra is 0 if neither logical nor boundary
    # 1 for logical, 2 for boundary, 3 for both
    # if there is no q1 or t, -1 is used
    cnodes = []
    for node in nodes:
        cnode = []
        cnode += node.qubits
        cnode += [-1] * (2 - len(node.qubits))
        if node.time is None:
            cnode.append(-1)
        else:
            cnode.append(node.time)
        cnode.append(1 * node.is_logical + 2 * node.is_boundary)
        cnodes.append(tuple(cnode))
    return cnodes
