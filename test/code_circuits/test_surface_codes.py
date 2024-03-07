# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2019.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# pylint: disable=invalid-name

"""Run codes and decoders."""

import unittest

from qiskit_qec.circuits.surface_code import SurfaceCodeCircuit
from qiskit_qec.utils import DecodingGraphNode


class TestSurfaceCodes(unittest.TestCase):
    """Test the surface code circuits."""

    def test_string2nodes(self):
        """
        Tests that the correct input emerges from string2nodes for a selection of inputs.
        """

        test_string = [
            "000000000 0000 0000",
            "000010000 0000 0110",
            "000010000 0110 0000",
            "000000001 0000 0000",
            "100000000 0000 0000",
        ]

        test_nodes = {}
        test_nodes["x"] = [
            [],
            [
                DecodingGraphNode(
                    time=1,
                    qubits=[0, 1, 3, 4],
                    index=1,
                ),
                DecodingGraphNode(
                    time=1,
                    qubits=[4, 5, 7, 8],
                    index=2,
                ),
            ],
            [
                DecodingGraphNode(
                    time=0,
                    qubits=[0, 1, 3, 4],
                    index=1,
                ),
                DecodingGraphNode(
                    time=0,
                    qubits=[4, 5, 7, 8],
                    index=2,
                ),
            ],
            [
                DecodingGraphNode(
                    is_logical=True,
                    is_boundary=True,
                    qubits=[0, 3, 6],
                    index=0,
                ),
                DecodingGraphNode(
                    time=1,
                    qubits=[0, 1, 3, 4],
                    index=1,
                ),
            ],
            [
                DecodingGraphNode(
                    is_logical=True,
                    is_boundary=True,
                    qubits=[2, 5, 8],
                    index=1,
                ),
                DecodingGraphNode(
                    time=1,
                    qubits=[4, 5, 7, 8],
                    index=2,
                ),
            ],
        ]
        test_nodes["z"] = [
            [],
            [
                DecodingGraphNode(
                    time=0,
                    qubits=[1, 4, 2, 5],
                    index=1,
                ),
                DecodingGraphNode(
                    time=0,
                    qubits=[3, 6, 4, 7],
                    index=2,
                ),
            ],
            [
                DecodingGraphNode(
                    time=1,
                    qubits=[1, 4, 2, 5],
                    index=1,
                ),
                DecodingGraphNode(
                    time=1,
                    qubits=[3, 6, 4, 7],
                    index=2,
                ),
            ],
            [
                DecodingGraphNode(is_logical=True, is_boundary=True, qubits=[0, 1, 2], index=0),
                DecodingGraphNode(time=1, qubits=[0, 3], index=0),
            ],
            [
                DecodingGraphNode(is_logical=True, is_boundary=True, qubits=[8, 7, 6], index=1),
                DecodingGraphNode(time=1, qubits=[5, 8], index=3),
            ],
        ]

        for basis in ["x", "z"]:
            code = SurfaceCodeCircuit(3, 1, basis=basis)
            for t, string in enumerate(test_string):
                nodes = test_nodes[basis][t]
                generated_nodes = code.string2nodes(string)
                self.assertTrue(
                    generated_nodes == nodes,
                    "Nodes for basis = "
                    + basis
                    + " and string = "
                    + string
                    + " are\
                    \n"
                    + str(generated_nodes)
                    + " not\n"
                    + str(nodes),
                )

    def test_check_nodes(self):
        """
        Tests for correct interpretation of a set of nodes.
        """
        valid = True

        # basis = 'z'
        code = SurfaceCodeCircuit(3, 3, basis="z")
        # trivial
        nodes = []
        valid = valid and code.check_nodes(nodes) == (True, [], 0)
        # on one side
        nodes = [
            DecodingGraphNode(qubits=[0, 1, 2], is_logical=True, is_boundary=True, index=0),
            DecodingGraphNode(time=3, qubits=[0, 3], index=0),
        ]
        valid = valid and code.check_nodes(nodes) == (True, [], 1.0)
        nodes = [DecodingGraphNode(time=3, qubits=[0, 3], index=0)]
        valid = valid and code.check_nodes(nodes) == (
            True,
            [DecodingGraphNode(qubits=[0, 1, 2], is_logical=True, is_boundary=True, index=0)],
            1.0,
        )
        # and the other
        nodes = [
            DecodingGraphNode(qubits=[8, 7, 6], is_logical=True, is_boundary=True, index=1),
            DecodingGraphNode(time=3, qubits=[5, 8], index=3),
        ]
        valid = valid and code.check_nodes(nodes) == (True, [], 1.0)
        nodes = [DecodingGraphNode(time=3, qubits=[5, 8], index=3)]
        valid = valid and code.check_nodes(nodes) == (
            True,
            [DecodingGraphNode(qubits=[8, 7, 6], is_logical=True, is_boundary=True, index=1)],
            1.0,
        )
        # and in the middle
        nodes = [
            DecodingGraphNode(time=3, qubits=[1, 4, 2, 5], index=1),
            DecodingGraphNode(time=3, qubits=[3, 6, 4, 7], index=2),
        ]
        valid = valid and code.check_nodes(nodes) == (True, [], 1)
        nodes = [DecodingGraphNode(time=3, qubits=[3, 6, 4, 7], index=2)]
        valid = valid and code.check_nodes(nodes) == (
            True,
            [DecodingGraphNode(qubits=[8, 7, 6], is_logical=True, is_boundary=True, index=1)],
            1.0,
        )

        # basis = 'x'
        code = SurfaceCodeCircuit(3, 3, basis="x")
        nodes = [
            DecodingGraphNode(time=3, qubits=[0, 1, 3, 4], index=1),
            DecodingGraphNode(time=3, qubits=[4, 5, 7, 8], index=2),
        ]
        valid = valid and code.check_nodes(nodes) == (True, [], 1)
        nodes = [DecodingGraphNode(time=3, qubits=[4, 5, 7, 8], index=2)]
        valid = valid and code.check_nodes(nodes) == (
            True,
            [DecodingGraphNode(qubits=[2, 5, 8], is_logical=True, is_boundary=True, index=1)],
            1.0,
        )

        # large d
        code = SurfaceCodeCircuit(5, 3, basis="z")
        nodes = [
            DecodingGraphNode(time=3, qubits=[7, 12, 8, 13], index=4),
            DecodingGraphNode(time=3, qubits=[11, 16, 12, 17], index=7),
        ]
        valid = valid and code.check_nodes(nodes) == (True, [], 1)
        nodes = [DecodingGraphNode(time=3, qubits=[11, 16, 12, 17], index=7)]
        valid = valid and code.check_nodes(nodes) == (
            True,
            [
                DecodingGraphNode(
                    qubits=[24, 23, 22, 21, 20], is_logical=True, is_boundary=True, index=1
                )
            ],
            2.0,
        )

        # wrong logical
        nodes = [
            DecodingGraphNode(time=3, qubits=[7, 12, 8, 13], index=4),
            DecodingGraphNode(
                qubits=[24, 23, 22, 21, 20], is_logical=True, is_boundary=True, index=1
            ),
        ]
        valid = valid and code.check_nodes(nodes) == (
            False,
            [DecodingGraphNode(qubits=[0, 1, 2, 3, 4], is_logical=True, is_boundary=True, index=0)],
            2,
        )
        # extra logical
        nodes = [
            DecodingGraphNode(time=3, qubits=[7, 12, 8, 13], index=4),
            DecodingGraphNode(time=3, qubits=[11, 16, 12, 17], index=7),
            DecodingGraphNode(
                qubits=[24, 23, 22, 21, 20], is_logical=True, is_boundary=True, index=1
            ),
        ]
        valid = valid and code.check_nodes(nodes) == (False, [], 0)
        # ignoring extra
        nodes = [
            DecodingGraphNode(time=3, qubits=[7, 12, 8, 13], index=4),
            DecodingGraphNode(time=3, qubits=[11, 16, 12, 17], index=7),
            DecodingGraphNode(
                qubits=[24, 23, 22, 21, 20], is_logical=True, is_boundary=True, index=1
            ),
        ]
        valid = valid and code.check_nodes(nodes, ignore_extras=True) == (True, [], 1)

        self.assertTrue(valid, "A set of nodes did not give the expected outcome for check_nodes.")
