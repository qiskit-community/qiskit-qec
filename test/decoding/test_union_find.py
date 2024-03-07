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
from qiskit_qec.circuits.repetition_code import ArcCircuit


class UnionFindDecoderTest(TestCase):
    """Tests for the UnionFind decoder not covered elsewhere"""

    def test_surface_code_d3(self):
        """
        Test the union find decoder on a surface code with d=3 and T=3
        with faults inserted into the final readout.
        """
        for logical in ["0", "1"]:
            code = SurfaceCodeCircuit(d=3, T=1)
            for use_is_cluster_neutral in [True, False]:
                decoder = UnionFindDecoder(code, use_is_cluster_neutral=use_is_cluster_neutral)
                for j in range(code.n):
                    string = (
                        logical * (j) + str((1 + int(logical)) % 2) + logical * (code.n - j - 1)
                    )
                    string += " 0000 0000"
                    corrected_outcome = decoder.process(string)
                    self.assertTrue(
                        corrected_outcome[0] == int(logical), "Correction for surface code failed."
                    )

    def test_hourglass_ARC(self):
        """
        Tests that clustering is done correctly on an ARC with (possibly misleading) loops
        """
        links = [
            (0, 1, 2),
            (0, 3, 4),
            (2, 5, 6),
            (4, 7, 8),
            (6, 9, 8),
            (8, 11, 10),
            (10, 13, 12),
            (12, 15, 14),
            (12, 17, 16),
            (14, 19, 18),
            (16, 21, 20),
            (18, 22, 20),
        ]
        code = ArcCircuit(links, 0)
        decoder = UnionFindDecoder(code)
        cluster = decoder.cluster(code.string2nodes("00001110000"))
        self.assertTrue(
            len(set(cluster.values())) == 1, "Clustering doesn't create single cluster."
        )


if __name__ == "__main__":
    unittest.main()
