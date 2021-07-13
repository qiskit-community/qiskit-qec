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
from typing import Optional

from PySide6.QtWidgets import QDialog, QMainWindow
from PySide6.QtCore import (QLineF, QPointF, QRect, QRectF, QSize, QSizeF, Qt,
                            Signal)
from PySide6.QtGui import (QAction, QColor, QFont, QIcon, QIntValidator,
                           QPainter, QPainterPath, QPen, QPixmap, QPolygonF,
                           QBrush, QKeyEvent)
from PySide6.QtWidgets import (QStyleOptionGraphicsItem,
    QApplication, QButtonGroup, QComboBox, QFontComboBox, QGraphicsAnchorLayout,
    QGraphicsItem, QGraphicsLineItem, QGraphicsPolygonItem, QGraphicsTextItem,
    QGraphicsEllipseItem, QGraphicsScene, QGraphicsView, QGridLayout,
    QHBoxLayout, QLabel, QMainWindow, QMenu, QMessageBox, QSizePolicy, QToolBox,
    QToolButton, QWidget, QPushButton, QVBoxLayout)
from enum import Enum
import diagramscene_rc
import uuid

print(diagramscene_rc)


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


class Qubit(QGraphicsEllipseItem):
    size = 50
    name = "QUBIT"

    def __init__(self, contextMenu, parent=None):
        super().__init__(0, 0, self.size, self.size, parent)

        self.groups = []
        
        self.uuid = uuid.uuid4()

        self._my_context_menu = contextMenu

        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)

    def add_group(self, stab):
        self.groups.append(stab)

    def remove_group(self, group):
        if group in self.groups:
            self.groups.remove(group)  # must come first to avoid recursion
            group.remove_qubit(self)  # TODO switch to signals

    # TODO deal with deleting stabs

    def image(self):
        pixmap = QPixmap(250, 250)
        pixmap.fill(Qt.transparent)
        painter = QPainter(pixmap)
        painter.setPen(QPen(Qt.black, 8))
        painter.translate(125, 125)

        path = QPainterPath(QPointF(-20, -20))
        path.moveTo(100, 0)
        path.arcTo(-100, -100, 200, 200, 0, 800)

        painter.drawPolyline(path.toFillPolygon())
        return pixmap

    def contextMenuEvent(self, event):
        self.scene().clearSelection()
        self.setSelected(True)
        self._my_context_menu.exec(event.screenPos())

    def itemChange(self, change, value):
        for gr in self.groups:
            gr.generate_polygon()
        return super().itemChange(change, value)

    def get_center(self):
        upper_left = self.scenePos()
        upper_left.setX(upper_left.x() + self.size / 2)
        upper_left.setY(upper_left.y() + self.size / 2)
        return upper_left


class GaugeGroup(QGraphicsPolygonItem):

    name = "Gauge Group"
    
    class PauliType(Enum):
        X = "X"
        Y = "Y"
        Z = "Z"
        EMPTY = "EMPTY"
        

    def __init__(self, qubits: list[Qubit]):

        super().__init__()
        self.setFillRule(Qt.OddEvenFill)
        self._error_groups = {}
        self._group_polies = set()

        self._qubits = []
        self.add_qubits(qubits)
        self.cur_points = []

        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        
        self.generate_polygon()

        
    
    @property
    def centroid(self):
        qlen = len(self._qubits)
        x_sum = 0
        y_sum = 0
        for q in self._qubits:
            x_sum += q.x()
            y_sum += q.y()
        centroid = QPointF(x_sum / qlen, y_sum / qlen)
        return centroid

    def less(self, qba: Qubit, qbb: Qubit):
        point_a = qba.get_center()
        point_b = qbb.get_center()
        if (point_a.x() - self.centroid.x() >=
                0) and (point_b.x() - self.centroid.x() < 0):
            return True
        if (point_a.x() - self.centroid.x() <
                0) and (point_b.x() - self.centroid.x() >= 0):
            return False

        if (point_a.x() - self.centroid.x()
                == 0) and (point_b.x() - self.centroid.x() == 0):
            if (point_a.y() - self.centroid.y() >=
                    0) or (point_b.y() - self.centroid.y() >= 0):
                return point_a.y() > point_b.y()
            return point_b.y() > point_a.y()

        # compute the cross product of vectors(center -> a) x(center -> b) int
        det = (point_a.x() -
               self.centroid.x()) * (point_b.y() - self.centroid.y()) - (
                   point_b.x() - self.centroid.x()) * (point_a.y() -
                                                       self.centroid.y())
        if (det < 0):
            return True
        if (det > 0):
            return False

        # points a and b are on the same line from the center

        # check which point is closer to the center
        d1 = (point_a.x() - self.centroid.x()) * (point_a.x() - self.centroid.x(
        )) + (point_a.y() - self.centroid.y()) * (point_a.y() -
                                                  self.centroid.y())
        d2 = (point_b.x() - self.centroid.x()) * (point_b.x() - self.centroid.x(
        )) + (point_b.y() - self.centroid.y()) * (point_b.y() -
                                                  self.centroid.y())
        return d1 > d2

    def bubble_sort_qubits(self):
        #bubble sort
        n = len(self._qubits)

        # Traverse through all array elements
        for i in range(n):

            # Last i elements are already in place
            for j in range(0, n - i - 1):

                # traverse the array from 0 to n-i-1
                # Swap if the element found is greater
                # than the next element
                if self.less(self._qubits[j], self._qubits[j + 1]):
                    self._qubits[j], self._qubits[j + 1] = self._qubits[
                        j + 1], self._qubits[j]

    def generate_polygon(self) -> [QPointF]:
        self.cur_points = []
        for q in self._qubits:
            self.cur_points.append(q.get_center())
        
        # todo how to properly call paint
        poly = QPolygonF(self.cur_points)
        self.setPolygon(poly)

    def remove_qubit(self, qubit):
        if qubit in self._qubits:
            self._qubits.remove(qubit)  # must come first to avoid recursion

            qubit.remove_group(self)  # TODO switch to signals
            self.bubble_sort_qubits()
            self.generate_polygon()
        if qubit in self._error_groups:
            self._error_groups.pop(qubit.uuid)

    def add_qubits(self, qubits):
        for qb in qubits:
            self._qubits.append(qb)
        self.bubble_sort_qubits()
        
        for qb_ind in range(len(self._qubits)):
            # just choose the type of the qubit next to it
            if self._qubits[qb_ind].uuid not in self._error_groups:
                if qb_ind == 0:
                    self._error_groups[self._qubits[qb_ind].uuid] = self.PauliType.EMPTY
                else:
                    self._error_groups[self._qubits[qb_ind].uuid] = self._error_groups[self._qubits[qb_ind - 1].uuid]

        self.generate_polygon()
        
        
    def set_entire_group_pauli(self, pauli: PauliType=PauliType.X):
        for key in self._error_groups.keys():
            self._error_groups[key] = pauli
        self.generate_polygon()
        
    
    def set_qubit_pauli(self, qubits: Qubit, pauli: PauliType=PauliType.X):
        for qu in qubits:
            self._error_groups[qu.uuid] = pauli
        self.generate_polygon()
        
        
    def paint(self, painter:QPainter, option: QStyleOptionGraphicsItem, widget: Optional[QWidget] = ...) -> None:
        pen = self.pen()
        
        painter.setPen(pen)
        qpoints = []
        for q in self._qubits:
            qpoints.append(q.get_center())
        poly = QPolygonF(qpoints)
        painter.drawPolygon(poly)
        print(f"ERROR GRPU:S {self._error_groups}")
        for qu in range(len(self._qubits)):
            pen.setColor('black')
            pauli = self._error_groups[self._qubits[qu].uuid]
            if pauli == self.PauliType.X:
                col = QColor('red')
            elif pauli == self.PauliType.Y:
                col = QColor('yellow')
            elif pauli == self.PauliType.Z:
                col = QColor('blue')
            elif pauli == self.PauliType.EMPTY:
                col = QColor('white')
                pen.setColor('white')
            else:
                continue

            pen.setColor(QColor('black'))
            painter.setBrush(col)

            qu_p = self._qubits[qu].get_center()
            centr_p = self.centroid
            p1_p = self._qubits[(qu + 1)%(len(self._qubits))].get_center()
            hp1_p = QPointF((p1_p.x() + qu_p.x())/2, (p1_p.y() + qu_p.y())/2)
            p2_p = self._qubits[(qu - 1)%(len(self._qubits))].get_center()
            hp2_p = QPointF((p2_p.x() + qu_p.x())/2, (p2_p.y() + qu_p.y())/2)
            
            poly = QPolygonF([hp1_p, qu_p, hp2_p, centr_p])
            painter.drawPolygon(poly)
            self._group_polies.add(poly)
            


class Stabilizer(GaugeGroup):
    name = "Stabilizer"

    def __init__(self, qubits: list):
        super().__init__(qubits)


class DiagramScene(QGraphicsScene):
    RETURN = 16777220
    DELETE = 16777223
    BACKSPACE = 16777219
    C_KEY = 67
    InsertItem, InsertLine, InsertText, MoveItem = range(4)

    item_inserted = Signal(Qubit)

    text_inserted = Signal(QGraphicsTextItem)

    item_selected = Signal(QGraphicsItem)

    def __init__(self, itemMenu, parent=None):
        super().__init__(parent)

        self._my_item_menu = itemMenu
        self._my_mode = self.MoveItem
        self.line = None
        self._text_item = None
        self._my_item_color = Qt.white
        self._my_text_color = Qt.black
        self._my_line_color = Qt.black
        self._my_font = QFont()

    def set_line_color(self, color):
        self._my_line_color = color

    def set_text_color(self, color):
        self._my_text_color = color
        if self.is_item_change(DiagramTextItem):
            item = self.selectedItems()[0]
            item.setDefaultTextColor(self._my_text_color)

    def set_item_color(self, color):
        self._my_item_color = color
        if self.is_item_change(Qubit):
            item = self.selectedItems()[0]
            item.setBrush(self._my_item_color)

    def set_font(self, font):
        self._my_font = font
        if self.is_item_change(DiagramTextItem):
            item = self.selectedItems()[0]
            item.setFont(self._my_font)

    def set_mode(self, mode):
        self._my_mode = mode

    def editor_lost_focus(self, item):
        cursor = item.textCursor()
        cursor.clearSelection()
        item.setTextCursor(cursor)

        if not item.toPlainText():
            self.removeItem(item)
            item.deleteLater()

    def mousePressEvent(self, mouseEvent):
        if (mouseEvent.button() != Qt.LeftButton):
            return

        if self._my_mode == self.InsertItem:
            item = Qubit(self._my_item_menu)
            item.setBrush(self._my_item_color)
            self.addItem(item)
            item.setPos(mouseEvent.scenePos())
            self.item_inserted.emit(item)
        elif self._my_mode == self.InsertText:
            text_item = DiagramTextItem()
            text_item.setFont(self._my_font)
            text_item.setTextInteractionFlags(Qt.TextEditorInteraction)
            text_item.setZValue(1000.0)
            text_item.lost_focus.connect(self.editor_lost_focus)
            text_item.selected_change.connect(self.item_selected)
            self.addItem(text_item)
            text_item.setDefaultTextColor(self._my_text_color)
            text_item.setPos(mouseEvent.scenePos())
            self.text_inserted.emit(text_item)

        super(DiagramScene, self).mousePressEvent(mouseEvent)

    def mouseMoveEvent(self, mouseEvent):
        if self._my_mode == self.InsertLine and self.line:
            new_line = QLineF(self.line.line().p1(), mouseEvent.scenePos())
            self.line.setLine(new_line)
        elif self._my_mode == self.MoveItem:
            super(DiagramScene, self).mouseMoveEvent(mouseEvent)

    def mouseReleaseEvent(self, mouseEvent):
        self.line = None
        super(DiagramScene, self).mouseReleaseEvent(mouseEvent)

    def is_item_change(self, type):
        for item in self.selectedItems():
            if isinstance(item, type):
                return True
        return False

    def keyPressEvent(self, event: QKeyEvent) -> None:
        if event.key() == self.C_KEY:
            for myi in self.items():
                if isinstance(myi, GaugeGroup):
                    for q in myi._error_groups.keys():
                        myi._error_groups[q] = GaugeGroup.PauliType.Y
                    myi.generate_polygon()
        if event.key() == self.DELETE or event.key() == self.BACKSPACE:
            for i in self.selectedItems():
                if isinstance(i, GaugeGroup):
                    for q in i._qubits:
                        q.remove_group(i)
                if isinstance(i, Qubit):
                    for g in i.groups:
                        g.remove_qubit(i)

                self.removeItem(i)

        if event.key() == self.RETURN:
            self.set_mode(self.InsertItem)
            # gb = self.SelectGroupTypeBox()
            # gb.exec()
            # if gb.group_class == gb.groupType.GAUGEGROUP:
            #     group = GaugeGroup
            # elif gb.group_class == gb.groupType.STABILIZER:
            #     group = Stabilizer
            # else:
            #     # TODO add logging
            #     return
            group = GaugeGroup
            points = []
            qubits = []
            for i in self.selectedItems():
                if isinstance(i, Qubit):
                    qubits.append(i)
                    points.append(i.get_center())
            group_item = group(qubits)
            self.addItem(group_item)
            for q in qubits:
                q.add_group(group_item)

        self.set_mode(self.MoveItem)
        for item in self.selectedItems():
            item.setSelected(False)

    class SelectGroupTypeBox(QDialog):

        class groupType(Enum):
            UNSET = "UNSET"
            INVALID = "INVALID"
            GAUGEGROUP = "Gauge Group"
            STABILIZER = "Stabilizer"

            def __str__(self):
                return self.value

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.setWindowTitle("Choose which group you wish to create")
            self._group_class = self.groupType.UNSET

            self.cur_layout = QVBoxLayout()
            self.setLayout(self.cur_layout)

            self._label = QLabel(
                f"Current group choice is: {str(self.group_class)}")
            self.cur_layout.addWidget(self._label)

            self._gauge_button = QPushButton("Gauge Operator", parent=self)
            self.cur_layout.addWidget(self._gauge_button)
            self._gauge_button.clicked.connect(
                lambda: self._set_group(self.groupType.GAUGEGROUP))

            self._stabilizer_button = QPushButton("Stabilizer Operator",
                                                  parent=self)
            self.cur_layout.addWidget(self._stabilizer_button)
            self._stabilizer_button.clicked.connect(
                lambda: self._set_group(self.groupType.STABILIZER))

            self._choice_layout = QHBoxLayout()
            self.cur_layout.addLayout(self._choice_layout)

            self._okay_button = QPushButton("Okay")
            self._choice_layout.addWidget(self._okay_button)
            self._okay_button.clicked.connect(self.accept)

            self._cancel_button = QPushButton("Cancel")
            self._choice_layout.addWidget(self._cancel_button)
            self._cancel_button.clicked.connect(self.reject)

        def update(self):
            self._label.setText(
                f"Current group choice is: {str(self.group_class)}")

        @property
        def group_class(self):
            return self._group_class

        def _set_group(self, value):
            if not isinstance(value, self.groupType):
                value = self.groupType.INVALID
            self._group_class = value
            self.update()
