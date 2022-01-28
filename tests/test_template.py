# This code is part of Qiskit.
#
# (C) Copyright IBM 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""Tests for template."""

from unittest import TestCase


class TestWillBeHere(TestCase):
    """Tests will be here."""

    def test_template_class(self):
        """Tests will start here."""
        mdrnlv = "MDRN LV"
        secret_message = "long drive"
        test_cases = {
            "MDRN LV": """Late night. Red wine. And you!
                        Late night, long drive, you can't kill our vibe."""
        }

        assert secret_message in test_cases[mdrnlv]
