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

"""Tests for template."""

import unittest
from unittest import TestCase
from qiskit_qec.decoders import UnionFindDecoder
from qiskit_qec.circuits import SurfaceCodeCircuit


class UnionFindDecoderTest(TestCase):
    """Tests for the UnionFind decoder not covered elsewhere"""

    def test_surface_code_d3(self):
        """
        Test the union find decoder on a surface code with d=3 and T=3
        with faults inserted into the final readout.
        """
        for logical in ["0", "1"]:
            code = SurfaceCodeCircuit(d=3, T=1)
            decoder = UnionFindDecoder(code)
            for j in range(code.n):
                string = logical * (j) + str((1 + int(logical)) % 2) + logical * (code.n - j - 1)
                string += " 0000 0000"
                corrected_outcome = decoder.process(string)
                self.assertTrue(
                    corrected_outcome[0] == int(logical), "Correction for surface code failed."
                )


if __name__ == "__main__":
    unittest.main()
