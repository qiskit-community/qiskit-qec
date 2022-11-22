"""Test pauli rep."""

from unittest import TestCase

from qiskit_qec.utils.visualizations import Screen
from qiskit_qec.utils.visualizations import QiskitGameEngine


class TestVisualization(TestCase):
    """Test pixel, screen and game engine."""

    def test_screen(self):
        """Tests basic screen and pixel functionality."""
        screen = Screen((400, 400), False, pixel_width=10)
        screen.pixel[9, 0].set_color("red")
        screen.pixel[0, 9].set_text("hello")
        self.assertEqual(
            screen.pixel[0, 0].button.button_style == "",
            "Default color not correct for pixel on screen.",
        )
        self.assertEqual(
            screen.pixel[9, 0].button.button_style == "danger",
            "Pixel on screen did not turn red when required.",
        )
        self.assertEqual(
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
        self.assertEqual(
            engine.screen.pixel[0, 0].button.button_style == "",
            "Default color not correct for pixel in game engine.",
        )
        self.assertEqual(
            engine.screen.pixel[engine.size - 1, 0].button.button_style == "danger",
            "Pixel in game engine did not turn red when required.",
        )
        self.assertEqual(
            engine.screen.pixel[0, engine.size - 1].button.description == "hello",
            "Pixel in game engine not displaying correct text",
        )
