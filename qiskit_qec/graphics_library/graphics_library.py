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

import math
import array
import sys
from typing import Optional, Tuple

from PySide6.QtWidgets import QDialog, QMainWindow
from PySide6.QtCore import (QLineF, QPointF, QRect, QRectF, QSize, QSizeF, Qt,
                            Signal)
from PySide6.QtGui import (QAction, QColor, QFont, QIcon, QIntValidator,
                           QPainter, QPainterPath, QPen, QPixmap, QPolygonF,
                           QBrush, QKeyEvent, QPolygon)
from PySide6.QtWidgets import (QStyleOptionGraphicsItem,
    QApplication, QButtonGroup, QComboBox, QFontComboBox, QGraphicsAnchorLayout,
    QGraphicsItem, QGraphicsLineItem, QGraphicsPolygonItem, QGraphicsTextItem,
    QGraphicsEllipseItem, QGraphicsScene, QGraphicsView, QGridLayout,
    QHBoxLayout, QLabel, QMainWindow, QMenu, QMessageBox, QSizePolicy, QToolBox,
    QToolButton, QWidget, QPushButton, QVBoxLayout, QGraphicsSceneMouseEvent)
from enum import Enum
from qiskit_qec.grid_tessellations.tessellation import Square
from typing import  Sequence, Union, Dict
import qiskit_qec.diagramscene_rc as diagramscene_rc
import uuid
import random
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


class DiagramTextItem(QGraphicsTextItem):
    lost_focus = Signal(QGraphicsTextItem)
    
    selected_change = Signal(QGraphicsItem)
    
    def __init__(self, parent=None, scene=None):
        super().__init__(parent, scene)
        
        self.setFlag(QGraphicsItem.ItemIsMovable)
        self.setFlag(QGraphicsItem.ItemIsSelectable)
    
    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemSelectedChange:
            self.selected_change.emit(self)
        return value
    
    def focusOutEvent(self, event):
        self.setTextInteractionFlags(Qt.NoTextInteraction)
        self.lost_focus.emit(self)
        super(DiagramTextItem, self).focusOutEvent(event)
    
    def mouseDoubleClickEvent(self, event):
        if self.textInteractionFlags() == Qt.NoTextInteraction:
            self.setTextInteractionFlags(Qt.TextEditorInteraction)
        super(DiagramTextItem, self).mouseDoubleClickEvent(event)


class GaugeGroup(QGraphicsPolygonItem):
    name = "Gauge Group"
    
    def __init__(self, points, pauli_map: Dict[ QPointF, PauliType ] = None, pauli_def: PauliType = PauliType.EMPTY):
        
        super().__init__(points)
        
        self.setFillRule(Qt.OddEvenFill)
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        
        if pauli_map is None:
            self.broken = False
            self._error_groups = dict()  # vertices -- Paulis
            
            if pauli_def:
                for p in points:
                    self.update_error_group(p, pauli_def)
            else:
                for p in points:
                    self.update_error_group(p, PauliType.EMPTY)
        
        else:
            self.broken = True
            self._error_groups = pauli_map
    
    def generate_polygon(self):
        # force paint event
        self.update(self.boundingRect())
    
    def update_error_group(self, point: QPointF, value: PauliType):
        key_point = (point.x(), point.y())
        self.update_error_group_using_key_point(key_point, value)
        
    def update_error_group_using_key_point(self,key_point:Tuple[float, float], value: PauliType):
        self.broken = True
        self._error_groups[key_point] = value
    
    def get_from_error_group(self, point: QPointF):
        key_point = (point.x(), point.y())
        return self._error_groups[ key_point ]
    
    def key_point_from_error_group(self):
        key_point_list = self._error_groups.keys()
        return key_point_list
    
    
    def is_point_in_error_group(self, point: QPointF):
        key_point = (point.x(), point.y())
        return key_point in self._error_groups
    
    def set_random_paulis(self):
        self.broken = True
        for q in self._error_groups.keys():
            options = list(PauliType)
            options.remove(PauliType.EMPTY)
            c = random.choice(options)
            self._error_groups[ q ] = c
        self.generate_polygon()
    
    def setPolygon(self, polygon: Union[ QPolygonF, Sequence[ QPointF ], QPolygon, QRectF ]) -> None:
        
        # update Paulis for rotation
        cur_poly = self.polygon()
        num_vert = len(cur_poly)
        if cur_poly != 0:
            new_values = [ ]
            for indx in range(num_vert):
                if self.is_point_in_error_group(cur_poly[ indx ]):
                    new_values.append((polygon[ indx ], self.get_from_error_group(cur_poly[ indx ])))
                else:
                    new_values.append((polygon[ indx ], PauliType.EMPTY))
            self._error_groups = dict()
            for tup in new_values:
                self.update_error_group(tup[ 0 ], tup[ 1 ])
        
        super().setPolygon(polygon)
    
    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[ QWidget ] = ...) -> None:
        
        pen = self.pen()
        painter.setPen(pen)
        
        # TODO figure out where to decide polygon shape
        
        if not self.broken:
            painter.setBrush(
                QColor(self.scene().get_pauli_type_color(random.choice(list(self._error_groups.values())))))
            painter.drawPolygon(self.polygon())
        
        else:
            x_sum = 0
            y_sum = 0
            count = 0
            for point in self.polygon():
                count += 1
                x_sum += point.x()
                y_sum += point.y()
            centroid = QPointF(x_sum / count, y_sum / count)
            
            qcolor = QColor('black')
            pen.setColor(qcolor)
            painter.setPen(pen)
            painter.drawPolygon(self.polygon())
            
            for indx in range(len(self.polygon())):
                poly = self.polygon()
                vertex = poly[ indx ]
                
                qcolor = QColor('salmon')
                qcolor.setAlpha(0)
                pen.setColor(qcolor)
                pauli = self.get_from_error_group(vertex)
                painter.setBrush(self.scene().get_pauli_type_color(pauli))
                painter.setPen(pen)
                
                p1_p = poly[ (indx + 1) % (len(poly)) ]
                hp1_p = QPointF((p1_p.x() + vertex.x()) / 2, (p1_p.y() + vertex.y()) / 2)
                p2_p = poly[ (indx - 1) % (len(poly)) ]
                hp2_p = QPointF((p2_p.x() + vertex.x()) / 2, (p2_p.y() + vertex.y()) / 2)
                
                poly = QPolygonF([ hp1_p, vertex, hp2_p, centroid ])
                painter.drawPolygon(poly)
                qcolor.setAlpha(255)
                pen.setColor(qcolor)
                painter.setPen(pen)
        
        if self.isSelected():
            pen = QPen(self.scene().HIGHLIGHT_COLOR)
            painter.setPen(pen)
            painter.pen().setWidth(5)
            trans = QColor("white")
            trans.setAlpha(0)
            painter.setBrush(trans)
            painter.drawPolygon(self.polygon())
    
    def set_entire_group_pauli(self, pauli: PauliType = PauliType.X):
        self.broken = False
        for key in self._error_groups.keys():
            self.update_error_group_using_key_point(key, pauli)
        self.generate_polygon()
    
    
    def mouseReleaseEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        current_poly = self.polygon()
        new_point = self.scene()._tiling.find_closest_vertex(self.scenePos())
        self.setPos(new_point)
        super().mouseReleaseEvent(event)


class Stabilizer(GaugeGroup):
    name = "Stabilizer"

    def __init__(self, qubits: list):
        super().__init__(qubits)


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

    def set_label_point(self, key_point: Tuple[float, float]):
        if key_point is not None:
            self.label.setText(f"Choosing color for vertex: {key_point}")
 
    @property
    def pauli_type(self):
        return self._pauli_type
    
    def set_pauli_type(self):
        self._pauli_type = PauliType(self._pauli_type_box.currentText())
        print(f"just set type only gt to {self._pauli_type} ")
        self.update()