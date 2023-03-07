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

from qiskit_qec.utils.decodoku import Decodoku


class TestDecodoku(unittest.TestCase):
    """Test the Decodoku class."""

    def test_graph_construction(self):
        """Test initialization and graph construction."""

        d = 9
        test_errors = [
            (0, 0, 1, 0),
            ((d - 1) / 2, d, (d - 1) / 2 + 1, d),
        ]
        correct_nodes = [
            {"y": 0, "x": 0, "is_boundary": False, "highlighted": True, "value": 1},
            {"y": 9, "x": 3, "is_boundary": False, "highlighted": True, "value": 1},
            {"y": 9, "x": 4, "is_boundary": False, "highlighted": True, "value": 1},
        ]

        game = Decodoku(d=d, k=2, errors=test_errors)
        highlighted_nodes = []
        for node in game.decoding_graph.graph.nodes():
            if node["highlighted"]:
                highlighted_nodes.append(node)

        self.assertTrue(len(highlighted_nodes) == 3, "Wrong number of nodes for example errors.")
        for node in correct_nodes:
            self.assertTrue(
                node in highlighted_nodes, "Expected node not present for example errors."
            )
