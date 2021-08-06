#############################################################################
##
## Copyright (C) 2013 Riverbank Computing Limited.
## Copyright (C) 2021 The Qt Company Ltd.
## Contact: http://www.qt.io/licensing/
##
## This file is part of the Qt for Python examples of the Qt Toolkit.
##
## $QT_BEGIN_LICENSE:BSD$
## You may use this file under the terms of the BSD license as follows:
##
## "Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##   * Neither the name of The Qt Company Ltd nor the names of its
##     contributors may be used to endorse or promote products derived
##     from this software without specific prior written permission.
##
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
##
## $QT_END_LICENSE$
##
#############################################################################

from enum import Enum
from typing import Tuple

from PySide6.QtWidgets import QComboBox, QDialog, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QSpinBox

import qiskit_qec.diagramscene_rc as diagramscene_rc

print(diagramscene_rc)


# Holds all the GraphicsItems put on display
# Holds all the Enums for dealing with ui

class PauliType(Enum):
    X = "X"
    Y = "Y"
    Z = "Z"
    EMPTY = "EMPTY"
    
    def __str__(self):
        return self.value


class GroupType(Enum):
    UNSET = "UNSET"
    INVALID = "INVALID"
    GAUGEGROUP = "Gauge Group"
    STABILIZER = "Stabilizer"
    
    def __str__(self):
        return self.value


class SelectGroupSectionTypeBox(QDialog):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setWindowTitle("Choose Type")
        self._pauli_type = PauliType.EMPTY
        
        self._type_layout = QVBoxLayout()
        self.setLayout(self._type_layout)
        
        self.label = QLabel()
        self._type_layout.addWidget(self.label)
        self.label.setText("Choose color for selected group(s)")
        
        self._pauli_type_box = QComboBox()
        self._type_layout.addWidget(self._pauli_type_box)
        
        for gt in PauliType:
            # Adds an item to the combobox with the given text, and containing the specified userData (stored in the Qt::UserRole). The item is appended to the list of existing items.
            self._pauli_type_box.addItem(str(gt))
        self._pauli_type_box.currentTextChanged.connect(
            self.set_pauli_type
        )
        self._pauli_type_box.setCurrentText(str(PauliType.EMPTY))
        
        self._choice_layout = QHBoxLayout()
        self._type_layout.addLayout(self._choice_layout)
        
        self._okay_button = QPushButton("Okay")
        self._choice_layout.addWidget(self._okay_button)
        self._okay_button.clicked.connect(self.accept)
    
    def set_label_point(self, key_point: Tuple[ float, float ]):
        if key_point is not None:
            self.label.setText(f"Choosing color for vertex: {key_point}")
    
    @property
    def pauli_type(self):
        return self._pauli_type
    
    def set_pauli_type(self):
        self._pauli_type = PauliType(self._pauli_type_box.currentText())
        print(f"just set type only gt to {self._pauli_type} ")
        self.update()


class ChooseSurfaceCode(QDialog):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setWindowTitle("Surface Code Maker")
        
        self._type_layout = QVBoxLayout()
        self.setLayout(self._type_layout)

        self.height_label = QLabel("Height:")
        self.height_box = QSpinBox()
        self._type_layout.addWidget(self.height_label)
        self._type_layout.addWidget(self.height_box)

        self.wid_label = QLabel("Width:")
        self.width_box = QSpinBox()
        self._type_layout.addWidget(self.wid_label)
        self._type_layout.addWidget(self.width_box)
        

        
        self.okay = QPushButton("Okay")
        self.okay.clicked.connect(self.close)
        self._type_layout.addWidget(self.okay)
        

