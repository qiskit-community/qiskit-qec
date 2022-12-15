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


class TestRepCodes(unittest.TestCase):
    """Test the surface code circuits."""

    def test_string2nodes(self):
        """
        Tests that the correct input emerges from string2nodes for a selection of inputs.
        """

        test_string = [
            '000000000 0000 0000',
            '000010000 0000 0110',
            '000010000 0110 0000',
            '000000001 0000 0000',
            '100000000 0000 0000',
        ]

        test_nodes = {}
        test_nodes['x'] = [
            [],
            [   {'time': 1, 'qubits': [0, 1, 3, 4], 'is_boundary': False, 'element': 1},
                {'time': 1, 'qubits': [4, 5, 7, 8], 'is_boundary': False, 'element': 2}],
            [
                {'time': 0, 'qubits': [0, 1, 3, 4], 'is_boundary': False, 'element': 1},
                {'time': 0, 'qubits': [4, 5, 7, 8], 'is_boundary': False, 'element': 2}
            ],
            [   {'time': 0, 'qubits': [0, 3, 6], 'is_boundary': True, 'element': 1},
                {'time': 1, 'qubits': [0, 1, 3, 4], 'is_boundary': False, 'element': 1}
            ],
            [
                {'time': 0, 'qubits': [2, 5, 8], 'is_boundary': True, 'element': 0},
                {'time': 1, 'qubits': [4, 5, 7, 8], 'is_boundary': False, 'element': 2}
            ]
        ]
        test_nodes['z'] = [
            [],
            [
                {'time': 0, 'qubits': [1, 4, 2, 5], 'is_boundary': False, 'element': 1},
                {'time': 0, 'qubits': [3, 6, 4, 7], 'is_boundary': False, 'element': 2}
            ],
            [
                {'time': 1, 'qubits': [1, 4, 2, 5], 'is_boundary': False, 'element': 1},
                {'time': 1, 'qubits': [3, 6, 4, 7], 'is_boundary': False, 'element': 2}
            ],
            [
                {'time': 0, 'qubits': [0, 1, 2], 'is_boundary': True, 'element': 1},
                {'time': 1, 'qubits': [0, 3], 'is_boundary': False, 'element': 0}
            ],
            [
                {'time': 0, 'qubits': [8, 7, 6], 'is_boundary': True, 'element': 0},
                {'time': 1, 'qubits': [5, 8], 'is_boundary': False, 'element': 3}
            ]
        ]

        for basis in ['x', 'z']:
            code = SurfaceCodeCircuit(3,1,basis=basis)
            for t, string in enumerate(test_string):
                nodes = test_nodes[basis][t]
                self.assertTrue(
                    code.string2nodes(string) == nodes,
                    "Incorrect nodes for basis = " + basis
                    + "for string = " + string + "."
                )