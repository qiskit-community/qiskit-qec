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


class ShapeObject:
    last_id = 0
    obj = dict()
    def __init__(self, stype=None, child=None) -> None:
        self.stype = stype
        self.id = self.create_id()
        self.parents = []
        if child is not None:
            ShapeObject.obj[self.id]=child
    
    def add_parent(self, parent):
        self.parents.append(parent)

    def create_id(self):
        ShapeObject.last_id += 1
        return ShapeObject.last_id

    def is_null(self):
        if self.otype is None:
            return True
        else:
            return False