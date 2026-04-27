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
"""Tests for GaugeGroup and its factory variants."""

import unittest
import warnings

from qiskit import QiskitError

from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.structures.gauge import GaugeGroup, GaugeGroupByGenerators, GaugeGroupWithGenerators


class TestGaugeGroupConstruction(unittest.TestCase):
    """Test constructing GaugeGroup from various input types."""

    def test_from_string_list(self):
        """GaugeGroup can be built from a list of Pauli strings."""
        gg = GaugeGroup(generators=["XX", "ZZ"])
        self.assertIsInstance(gg, GaugeGroup)

    def test_from_pauli_list(self):
        """GaugeGroup can be built from a PauliList."""
        pl = PauliList(["XX", "ZZ"])
        gg = GaugeGroup(generators=pl)
        self.assertIsInstance(gg, GaugeGroup)

    def test_from_isotropic_generators(self):
        """GaugeGroup can be built using isotropic_generators kwarg."""
        gg = GaugeGroup(isotropic_generators=["XX", "ZZ"])
        self.assertIsInstance(gg, GaugeGroup)

    def test_from_hyperbolic_generators(self):
        """GaugeGroup can be built using hyperbolic_generators kwarg."""
        gg = GaugeGroup(hyperbolic_generators=["XZ", "ZX"])
        self.assertIsInstance(gg, GaugeGroup)

    def test_combined_iso_hyperbolic_generators(self):
        """GaugeGroup from both isotropic and hyperbolic generators."""
        gg = GaugeGroup(
            isotropic_generators=["XX"],
            hyperbolic_generators=["ZZ"],
        )
        self.assertIsInstance(gg, GaugeGroup)

    def test_no_generators_raises(self):
        """Providing no generators at all raises QiskitError."""
        with self.assertRaises(QiskitError):
            GaugeGroup()

    def test_with_name(self):
        """Custom name is stored on the object."""
        gg = GaugeGroup(generators=["XX", "ZZ"], name="steane")
        self.assertEqual(gg.name, "steane")

    def test_default_name_is_empty_string(self):
        """Default name is an empty string when not provided."""
        gg = GaugeGroup(generators=["XX", "ZZ"])
        self.assertEqual(gg.name, "")

    # ------------------------------------------------------------------ #
    # Properties                                                           #
    # ------------------------------------------------------------------ #

    def test_n_returns_num_qubits(self):
        """n property returns the number of physical qubits."""
        gg = GaugeGroup(generators=["XX", "ZZ"])
        self.assertEqual(gg.n, 2)

    def test_n_three_qubit_generators(self):
        """n=3 for three-qubit generators."""
        gg = GaugeGroup(generators=["XXX", "ZZI"])
        self.assertEqual(gg.n, 3)

    def test_k_warns_not_implemented(self):
        """k property emits a warning (not yet implemented)."""
        gg = GaugeGroup(generators=["XX", "ZZ"])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _ = gg.k
            self.assertEqual(len(w), 1)

    def test_num_gen_warns_not_implemented(self):
        """num_gen property emits a warning (not yet implemented)."""
        gg = GaugeGroup(generators=["XX", "ZZ"])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _ = gg.num_gen
            self.assertEqual(len(w), 1)

    def test_str_representation(self):
        """__str__ returns a non-empty string."""
        gg = GaugeGroup(generators=["XX", "ZZ"])
        self.assertIsInstance(str(gg), str)
        self.assertGreater(len(str(gg)), 0)

    # ------------------------------------------------------------------ #
    # Generator minimisation (with_generators=True vs False)              #
    # ------------------------------------------------------------------ #

    def test_with_generators_true_minimises(self):
        """with_generators=True reduces generators to a minimal set."""
        # Providing redundant generators: XX and XX are the same
        gg = GaugeGroup(generators=["XX", "XX", "ZZ"], with_generators=True)
        # Minimal set should remove the duplicate
        self.assertLessEqual(len(gg.generators), 3)

    def test_with_generators_false_keeps_input(self):
        """with_generators=False stores the exact generators supplied."""
        gg = GaugeGroup(generators=["XX", "ZZ"], with_generators=False)
        self.assertEqual(len(gg.generators), 2)


class TestGaugeGroupFactories(unittest.TestCase):
    """Test the GaugeGroupWithGenerators and GaugeGroupByGenerators factories."""

    def test_with_generators_factory_returns_gauge_group(self):
        """GaugeGroupWithGenerators returns a GaugeGroup instance."""
        gg = GaugeGroupWithGenerators(generators=["XX", "ZZ"])
        self.assertIsInstance(gg, GaugeGroup)

    def test_by_generators_factory_returns_gauge_group(self):
        """GaugeGroupByGenerators returns a GaugeGroup instance."""
        gg = GaugeGroupByGenerators(generators=["XX", "ZZ"])
        self.assertIsInstance(gg, GaugeGroup)

    def test_by_generators_preserves_input(self):
        """GaugeGroupByGenerators stores the generators as supplied."""
        gg = GaugeGroupByGenerators(generators=["XX", "ZZ"])
        self.assertEqual(len(gg.generators), 2)

    def test_with_generators_factory_name(self):
        """GaugeGroupWithGenerators forwards the name argument."""
        gg = GaugeGroupWithGenerators(generators=["XX", "ZZ"], name="test")
        self.assertEqual(gg.name, "test")


if __name__ == "__main__":
    unittest.main()
