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
"""Module for ShapeObject"""


class ShapeObject:  # pylint: disable=function-redefined)
    """ShapeObject is the base class for all geometry objects"""

    last_id = 0

    def __init__(self) -> None:
        """Init ShapeObject

        Args:
            children (Union, optional): Defaults to None.
        """
        self.parents = []
        self.id = self.create_id()  # pylint: disable=invalid-name

    def add_parent(self, parent: "ShapeObject") -> None:
        """Adds parent

        Args:
            parent (ShapeObject): Adds parent
        """
        self.parents.append(parent)

    @staticmethod
    def create_id() -> int:
        """Creates int as unique ID for ShapeObject instance.
        ShapeObject class tracks which IDs have previously been used and are thus not available.

        Returns:
            int: Unique ID
        """
        ShapeObject.last_id += 1
        return ShapeObject.last_id
