# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2019-2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Visualization Module"""
from unittest import TestCase

from qiskit_qec.utils.visualizations import QiskitGameEngine, Screen


class TestVisualization(TestCase):
    """Test pixel, screen and game engine."""

    def test_screen(self):
        """Tests basic screen and pixel functionality."""
        screen = Screen((400, 400), False, pixel_width=10)
        screen.pixel[9, 0].set_color("red")
        screen.pixel[0, 9].set_text("hello")
        self.assertTrue(
            screen.pixel[0, 0].button.button_style == "",
            "Default color not correct for pixel on screen.",
        )
        self.assertTrue(
            screen.pixel[9, 0].button.button_style == "danger",
            "Pixel on screen did not turn red when required.",
        )
        self.assertTrue(
            screen.pixel[0, 9].button.description == "hello",
            "Pixel on screen not displaying correct text",
        )

    def test_engine(self):
        """Test basic game engine functionality."""

        def start(engine):
            engine.screen.pixel[engine.size - 1, 0].set_color("red")
            engine.screen.pixel[0, engine.size - 1].set_text("hello")

        def next_frame(_):
            pass

        engine = QiskitGameEngine(start, next_frame, size=8)
        self.assertTrue(
            engine.screen.pixel[0, 0].button.button_style == "",
            "Default color not correct for pixel in game engine.",
        )
        self.assertTrue(
            engine.screen.pixel[engine.size - 1, 0].button.button_style == "danger",
            "Pixel in game engine did not turn red when required.",
        )
        self.assertTrue(
            engine.screen.pixel[0, engine.size - 1].button.description == "hello",
            "Pixel in game engine not displaying correct text",
        )
