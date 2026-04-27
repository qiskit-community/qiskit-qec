# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2022
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Tests for SurfaceCodeBuilder."""

import unittest

from qiskit import QiskitError

from qiskit_qec.codes.codebuilders.surface_code_builder import SurfaceCodeBuilder
from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode
from qiskit_qec.operators.pauli import Pauli


class TestSurfaceCodeBuilder(unittest.TestCase):
    """Test SurfaceCodeBuilder."""

    # ------------------------------------------------------------------ #
    # Construction — valid inputs                                          #
    # ------------------------------------------------------------------ #

    def test_square_d3_returns_stab_code(self):
        """d=3 square surface code returns a StabSubSystemCode."""
        code = SurfaceCodeBuilder(d=3).build()
        self.assertIsInstance(code, StabSubSystemCode)

    def test_square_d5_returns_stab_code(self):
        """d=5 square surface code returns a StabSubSystemCode."""
        code = SurfaceCodeBuilder(d=5).build()
        self.assertIsInstance(code, StabSubSystemCode)

    def test_square_d7_returns_stab_code(self):
        """d=7 square surface code is built without error."""
        code = SurfaceCodeBuilder(d=7).build()
        self.assertIsInstance(code, StabSubSystemCode)

    def test_rectangular_dx3_dz5(self):
        """Rectangular code with dx=3, dz=5 builds successfully."""
        code = SurfaceCodeBuilder(dx=3, dz=5).build()
        self.assertIsInstance(code, StabSubSystemCode)

    def test_rectangular_dx5_dz7(self):
        """Rectangular code with dx=5, dz=7 builds successfully."""
        code = SurfaceCodeBuilder(dx=5, dz=7).build()
        self.assertIsInstance(code, StabSubSystemCode)

    def test_ul_op_Z_default(self):
        """Default ul_op=Z is accepted and produces a valid code."""
        code = SurfaceCodeBuilder(d=3, ul_op=Pauli("Z")).build()
        self.assertIsInstance(code, StabSubSystemCode)

    def test_ul_op_X(self):
        """ul_op=X produces a valid code with transposed operator layout."""
        code = SurfaceCodeBuilder(d=3, ul_op=Pauli("X")).build()
        self.assertIsInstance(code, StabSubSystemCode)

    def test_d_overrides_dx_dz(self):
        """When d is provided, dx/dz are ignored; code still builds."""
        # Passing dx and dz alongside d should silently use d only
        code = SurfaceCodeBuilder(d=3).build()
        self.assertIsInstance(code, StabSubSystemCode)

    # ------------------------------------------------------------------ #
    # Construction — invalid inputs                                        #
    # ------------------------------------------------------------------ #

    def test_even_d_raises(self):
        """Even distance raises QiskitError."""
        with self.assertRaises(QiskitError):
            SurfaceCodeBuilder(d=4)

    def test_d_less_than_3_raises(self):
        """d=1 raises QiskitError (minimum allowed is 3)."""
        with self.assertRaises(QiskitError):
            SurfaceCodeBuilder(d=1)

    def test_even_dx_raises(self):
        """Even dx raises QiskitError."""
        with self.assertRaises(QiskitError):
            SurfaceCodeBuilder(dx=4, dz=5)

    def test_even_dz_raises(self):
        """Even dz raises QiskitError."""
        with self.assertRaises(QiskitError):
            SurfaceCodeBuilder(dx=3, dz=4)

    def test_missing_dz_raises(self):
        """Specifying dx without dz raises QiskitError."""
        with self.assertRaises(QiskitError):
            SurfaceCodeBuilder(dx=3)

    def test_missing_dx_raises(self):
        """Specifying dz without dx raises QiskitError."""
        with self.assertRaises(QiskitError):
            SurfaceCodeBuilder(dz=3)

    # ------------------------------------------------------------------ #
    # Code properties                                                      #
    # ------------------------------------------------------------------ #

    def test_d3_has_gauge_group(self):
        """d=3 code object exposes a gauge_group attribute."""
        code = SurfaceCodeBuilder(d=3).build()
        self.assertTrue(hasattr(code, "gauge_group"))

    def test_d3_qubit_count(self):
        """d=3 surface code has n = d² + (d-1)² = 13 physical qubits."""
        code = SurfaceCodeBuilder(d=3).build()
        # Non-rotated surface code: data qubits (d²) + ancilla qubits ((d-1)²)
        self.assertEqual(code._n, 3 ** 2 + 2 ** 2)

    def test_d5_qubit_count(self):
        """d=5 surface code has n = d² + (d-1)² = 41 physical qubits."""
        code = SurfaceCodeBuilder(d=5).build()
        self.assertEqual(code._n, 5 ** 2 + 4 ** 2)

    def test_ul_op_sets_optype_Z(self):
        """ul_op=Z sets optype to 'pZXZX'."""
        builder = SurfaceCodeBuilder(d=3, ul_op=Pauli("Z"))
        self.assertEqual(builder.optype, "pZXZX")

    def test_ul_op_sets_optype_X(self):
        """ul_op=X sets optype to 'pXZXZ'."""
        builder = SurfaceCodeBuilder(d=3, ul_op=Pauli("X"))
        self.assertEqual(builder.optype, "pXZXZ")


if __name__ == "__main__":
    unittest.main()
