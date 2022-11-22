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

"""Tools for visualization and interactivity."""
from ipywidgets import widgets
from ipywidgets import Layout
from IPython.display import display


class Pixel:
    """Pixel made from Jupyter widget buttons."""

    def __init__(self, layout, active=False):
        """
        Args:
            layout (ipywidgets.Layout): Layout in which which the pixels are embedded.
            active (bool): Whether the pixels should be clickable.
        """
        self.disabled = not active
        self.button = widgets.ToggleButton(
            description="", button_style="", layout=layout, disabled=self.disabled
        )

    def set_color(self, color):
        """Set color of the pixel.

        Args:
            color (str): Must be "grey"/"gray", "green", "blue", "orange" or "red"
        """
        if color in ["grey", "gray"]:
            self.button.button_style = ""
        elif color == "green":
            self.button.button_style = "success"
        elif color == "blue":
            self.button.button_style = "info"
        elif color == "orange":
            self.button.button_style = "warning"
        elif color == "red":
            self.button.button_style = "danger"

    def set_brightness(self, bright):
        """
        Set brightness of the pixel.

        Args:
            bright (bool): Whether the pixel should be bright.
        """
        if self.disabled:
            self.button.value = not bright

    def set_text(self, text):
        """
        Set text for the pixel.

        Args:
            text (str): Text to display.
        """
        self.button.description = text


class Screen:
    """Screen made from a grid of 'Pixel' objects."""

    def __init__(self, size, active, pixel_width=8):
        """
        Args:
            size (tuple): Width and height of the screen.
            active (bool): Whether the pixels should be clickable.
            pixel_width (int): Width and height of the screen in pixels.
        """

        width = int(size[0] / pixel_width)
        wide = str(7 * width + 24) + "px"
        wider = str(pixel_width * width + (pixel_width - 1) * 4) + "px"
        width = str(width) + "px"
        height = str(int(size[1] / pixel_width)) + "px"
        width = str(int(50 * 8 / pixel_width)) + "px"

        self._layout = Layout(width=width, height=height)
        self._wide_layout = Layout(width=wide, height=height)
        self._wider_layout = Layout(width=wider, height=height)

        self.pixel = {}
        for x in range(pixel_width):
            for y in range(pixel_width):
                self.pixel[x, y] = Pixel(self._layout, active)
        self.pixel["text"] = Pixel(self._wider_layout)


class QiskitGameEngine:
    """Class to facilitate interactivity on a 'Screen' object."""

    def __init__(self, start, next_frame, size=8, active_screen=False):
        """
        Args:
            start (Callable): Function to initialize the screen.
            next_frame (Callable): Function to call when a button is pressed.
            size (int): Width and height of the screen in pixels.
            active_screen (bool): Whether the pixels are pressable.
        """

        self.start = start
        self.next_frame = next_frame
        self.size = size
        self.active_screen = active_screen
        self.pressed_pixels = []

        self.screen = Screen((400, 400), active_screen, pixel_width=size)
        layout = self.screen._layout

        controller = {}
        controller["blank"] = widgets.ToggleButton(description="", button_style="", layout=layout)
        controller["up"] = widgets.ToggleButton(description="▲", button_style="", layout=layout)
        controller["down"] = widgets.ToggleButton(description="▼", button_style="", layout=layout)
        controller["left"] = widgets.ToggleButton(description="◀︎", button_style="", layout=layout)
        controller["right"] = widgets.ToggleButton(description="►", button_style="", layout=layout)
        controller["A"] = widgets.ToggleButton(description="A", button_style="", layout=layout)
        controller["B"] = widgets.ToggleButton(description="B", button_style="", layout=layout)
        controller["X"] = widgets.ToggleButton(description="X", button_style="", layout=layout)
        controller["Y"] = widgets.ToggleButton(description="Y", button_style="", layout=layout)
        controller["next"] = widgets.ToggleButton(
            description="Next", button_style="", layout=self.screen._wide_layout
        )

        [blank, up, down, left, right, a_button, b_button, x_button, y_button, next_button] = [
            controller["blank"],
            controller["up"],
            controller["down"],
            controller["left"],
            controller["right"],
            controller["A"],
            controller["B"],
            controller["X"],
            controller["Y"],
            controller["next"],
        ]

        interface = []
        interface.append(
            widgets.HBox(
                [self.screen.pixel[x, 0].button for x in range(size)]
                + [blank, up, blank, blank, blank, x_button, blank]
            )
        )
        interface.append(
            widgets.HBox(
                [self.screen.pixel[x, 1].button for x in range(size)]
                + [left, blank, right, blank, y_button, blank, a_button]
            )
        )
        interface.append(
            widgets.HBox(
                [self.screen.pixel[x, 2].button for x in range(size)]
                + [blank, down, blank, blank, blank, b_button, blank]
            )
        )
        interface.append(
            widgets.HBox([self.screen.pixel[x, 3].button for x in range(size)] + [next_button])
        )
        for y in range(4, size):
            interface.append(widgets.HBox([self.screen.pixel[x, y].button for x in range(size)]))
        interface.append(self.screen.pixel["text"].button)

        self.controller = controller

        # run user-supplied setup function
        start(self)

        display(widgets.VBox(interface))

        blank.observe(self._given_blank)

        for button, button_obj in self.controller.items():
            if button != "blank":
                button_obj.observe(self._given_button)

        if active_screen:
            for pixel in self.screen.pixel.values():
                pixel.button.observe(self._given_screen)

    def _given_blank(self, _):
        if self.controller["blank"].value:
            self.controller["blank"].value = False

    def _given_button(self, _):
        for button in self.controller.values():
            if button.value:
                self.next_frame(self)
            button.value = False

    def _given_screen(self, _):
        if self.active_screen:
            for pos, pixel in self.screen.pixel.items():
                if pixel.button.value:
                    self.pressed_pixels.append(pos)
                    self.next_frame(self)
                pixel.button.value = False
