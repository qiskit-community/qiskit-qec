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

"""Run codes and decoders."""

import unittest

from qiskit_qec.circuits import ArcCircuit
from qiskit_qec.decoders import DecodingGraph


class TestCleaner(unittest.TestCase):
    """Test predecoders"""

    def test_measurement_cleaner(self):
        """Test measurement cleaner"""
        links = [(0, 1, 2), (2, 3, 4), (0, 5, 8), (2, 6, 10), (4, 7, 12), (8, 9, 10), (10, 11, 12)]
        code = ArcCircuit(links, 3)
        decoding_graph = DecodingGraph(code)

        # test that isolated measurement errors are removed
        string = "010000 0000000 0000000 0100100"
        cleaned_nodes = decoding_graph.clean_measurements(code.string2nodes(string))
        self.assertTrue(
            len(cleaned_nodes) == 3,
            "Wrong number of cleaned nodes for isolated measurement errors.",
        )

        # test that neighbouring ones aren't
        string = "000000 0000100 0000000 0000100"
        cleaned_nodes = decoding_graph.clean_measurements(code.string2nodes(string))
        self.assertTrue(
            len(cleaned_nodes) == 4,
            "Wrong number of cleaned nodes for neighbouring measurement errors.",
        )

    def test_code_cleaner(self):
        """Test code cleaner"""
        links = [(0, 1, 2), (2, 3, 4), (0, 5, 8), (2, 6, 10), (4, 7, 12), (8, 9, 10), (10, 11, 12)]
        code = ArcCircuit(links, 2)

        error_string = "Single error handled incorrectly by code cleaner"
        # test that single errors are corrected
        for j in range(6):
            string = "1" * j + "0" + "1" * (6 - j - 1) + " 0000000 0000000"
            self.assertTrue(code.clean_code(string) == "111111 0000000 0000000", error_string)
        test_strings = [
            "000010 0001011 0001011",
            "000010 0001011 0000000",
            "000001 0000101 0000101",
        ]
        for string in test_strings:
            self.assertTrue(code.clean_code(string) == "000000 0000000 0000000", error_string)
        self.assertTrue(
            code.clean_code("000010 0001011 0001000") == "000000 0000000 0001000", error_string
        )
        error_string = "code cleaer acts non-trivially on ambiguous syndrome"
        # test syndromes that it shouldn't tackle
        self.assertTrue(
            code.clean_code("000000 0000000 1000001") == "000000 0000000 1000001", error_string
        )
        self.assertTrue(
            code.clean_code("000000 0000000 1000001") == "000000 0000000 1000001", error_string
        )
        self.assertTrue(
            code.clean_code("000000 0000011 0000011") == "000000 0000011 0000011", error_string
        )
        # test a case where it causes a logical error
        self.assertTrue(
            code.clean_code("101010 0001011 0000000") == "000000 0000000 0000000",
            "Code cleaner acts incorrectly on complex syndrome",
        )
