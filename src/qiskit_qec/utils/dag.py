# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

""" Utililities for Qiskit DAGs"""

from qiskit.dagcircuit.dagnode import DAGNode


def node_name_label(node: DAGNode) -> str:
    """Form an identifier string for a node's operation.
    Use node.op._label if it exists. Otherwise use node.name.
    Return a string.
    """
    if "_label" in node.op.__dict__ and node.op._label is not None:
        name_label = node.op._label
    else:
        name_label = node.name
    return name_label
