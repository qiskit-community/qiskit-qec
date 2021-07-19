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
from qiskit_qec.grid_tessellations.tessellation import Square

import qiskit_qec.diagramscene_rc as diagramscene_rc
import uuid

print(diagramscene_rc)



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


class Qubit(QGraphicsEllipseItem):
    size = 50
    name = "QUBIT"

    def __init__(self, contextMenu, parent=None, x_pos=0, y_pos=0):
        super().__init__(x_pos, y_pos, self.size, self.size, parent)

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


    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[QWidget] = ...) -> None:
        super().paint(painter, option, widget)
        if self.isSelected():
            pen = QPen(self.scene().HIGHLIGHT_COLOR)

            painter.setPen(pen)
            painter.pen().setWidth(2)
            trans = QColor("white")
            trans.setAlpha(0)
            painter.setBrush(trans)
            painter.drawPolygon(self.boundingRect())


class GaugeGroup(QGraphicsPolygonItem):
    name = "Gauge Group"
    
    def __init__(self, qubits: list[ Qubit ], group_type: PauliType = None):
        
        super().__init__()
        self.setFillRule(Qt.OddEvenFill)
        self._error_groups = {}
        self.setup_qubits_and_paulis(qubits, group_type)
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        
        self.generate_polygon()
    
    
    def setup_qubits_and_paulis(self, qubits, group_type: PauliType):
        self.broken = False
        self.qubits = [ ]
        self.cur_points = [ ]
        self.add_qubits(qubits)
        self._group_paulies = set()
        if group_type is not None:
            self.set_entire_group_pauli(group_type)
    
    
    @property
    def centroid(self):
        qlen = len(self.qubits)
        x_sum = 0
        y_sum = 0
        for q in self.qubits:
            x_sum += q.get_center().x()
            y_sum += q.get_center().y()
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
        # bubble sort
        n = len(self.qubits)
        
        # Traverse through all array elements
        for i in range(n):
            
            # Last i elements are already in place
            for j in range(0, n - i - 1):
                
                # traverse the array from 0 to n-i-1
                # Swap if the element found is greater
                # than the next element
                if self.less(self.qubits[ j ], self.qubits[ j + 1 ]):
                    self.qubits[ j ], self.qubits[ j + 1 ] = self.qubits[
                                                                 j + 1 ], self.qubits[ j ]
    
    def generate_polygon(self) -> [ QPointF ]:
        self.cur_points = [ ]
        for q in self.qubits:
            self.cur_points.append(q.get_center())
        
        # todo how to properly call paint
        poly = QPolygonF(self.cur_points)
        self.setPolygon(poly)
    
    def remove_qubit(self, qubit):
        if qubit in self.qubits:
            self.qubits.remove(qubit)  # must come first to avoid recursion
            
            qubit.remove_group(self)  # TODO switch to signals
            self.bubble_sort_qubits()
            self.generate_polygon()
        if qubit in self._error_groups:
            self._error_groups.pop(qubit.uuid)
    
    def add_qubits(self, qubits):
        for qb in qubits:
            self.qubits.append(qb)
        self.bubble_sort_qubits()
        
        for qb_ind in range(len(self.qubits)):
            # just choose the type of the qubit next to it
            if self.qubits[ qb_ind ].uuid not in self._error_groups:
                if qb_ind == 0:
                    self._error_groups[ self.qubits[ qb_ind ].uuid ] = PauliType.EMPTY
                else:
                    self._error_groups[ self.qubits[ qb_ind ].uuid ] = self._error_groups[
                        self.qubits[ qb_ind - 1 ].uuid ]
        
        self.generate_polygon()
    
    
    def set_entire_group_pauli(self, pauli: PauliType = PauliType.X):
        self.broken = False
        for key in self._error_groups.keys():
            self._error_groups[ key ] = pauli
        self.generate_polygon()
    
    
    def set_qubit_pauli(self, qubits: Qubit, pauli: PauliType = PauliType.X):
        self.broken = True
        for qu in qubits:
            self._error_groups[ qu.uuid ] = pauli
        self.generate_polygon()
    
    def set_random_paulis(self):
        for q in self._error_groups.keys():
            options = list(PauliType)
            options.remove(PauliType.EMPTY)
            c = random.choice(options)
        self._error_groups[ q ] = c
        self.generate_polygon()
    
    
    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[ QWidget ] = ...) -> None:
        
        pen = self.pen()
        painter.setPen(pen)
        
        if not self.broken:
            painter.setBrush(QColor(self.scene().get_group_type_color(self._error_groups[ self.qubits[ 0 ].uuid ])))
            qpoints = [ ]
            for q in self.qubits:
                qpoints.append(q.get_center())
            poly = QPolygonF(qpoints)
            painter.drawPolygon(poly)
        
        else:
            for qu in range(len(self.qubits)):
                pen.setColor('black')
                pauli = self._error_groups[ self.qubits[ qu ].uuid ]
                if pauli == PauliType.EMPTY:
                    pen.setColor('white')
                painter.setBrush(self.scene().get_group_type_color(pauli))
                
                qu_p = self.qubits[ qu ].get_center()
                centr_p = self.centroid
                p1_p = self.qubits[ (qu + 1) % (len(self.qubits)) ].get_center()
                hp1_p = QPointF((p1_p.x() + qu_p.x()) / 2, (p1_p.y() + qu_p.y()) / 2)
                p2_p = self.qubits[ (qu - 1) % (len(self.qubits)) ].get_center()
                hp2_p = QPointF((p2_p.x() + qu_p.x()) / 2, (p2_p.y() + qu_p.y()) / 2)
                
                poly = QPolygonF([ hp1_p, qu_p, hp2_p, centr_p ])
                painter.drawPolygon(poly)
                self._group_paulies.add(poly)
        
        if self.isSelected():
            pen = QPen(self.scene().HIGHLIGHT_COLOR)
            painter.setPen(pen)
            painter.pen().setWidth(5)
            trans = QColor("white")
            trans.setAlpha(0)
            painter.setBrush(trans)
            qpoints = [ ]
            for q in self.qubits:
                qpoints.append(q.get_center())
            
            painter.drawPolygon(QPolygonF(qpoints))


class Stabilizer(GaugeGroup):
    name = "Stabilizer"

    def __init__(self, qubits: list):
        super().__init__(qubits)


class DiagramScene(QGraphicsScene):
    RETURN = 16777220
    DELETE = 16777223
    BACKSPACE = 16777219
    KEY_C = 67
    KEY_R = 81
    KEY_1 = 49
    KEY_2 = 50
    KEY_3 = 51
    KEY_4 = 52
    C_KEY = 67

    InsertItem, InsertLine, InsertText, MoveItem = range(4)
    HIGHLIGHT_COLOR = QColor('blue')
    item_inserted = Signal(Qubit)

    text_inserted = Signal(QGraphicsTextItem)

    item_selected = Signal(QGraphicsItem)

    def __init__(self, itemMenu, parent=None):
        super().__init__(parent)
        self.setup_group_type_color_map()
        self._my_item_menu = itemMenu
        self._my_mode = self.MoveItem
        self.line = None
        self._text_item = None
        self._my_item_color = Qt.white
        self._my_text_color = Qt.black
        self._my_line_color = Qt.black
        self._my_font = QFont()
        
        
    def setup_group_type_color_map(self):
        self.group_type_color_map = {}
        # todo make random if becomes too many
        self.group_type_color_map[PauliType.X] = QColor('beige')
        self.group_type_color_map[PauliType.Y] = QColor('darkseagreen')
        self.group_type_color_map[ PauliType.Z ] = QColor('darksalmon')
        self.group_type_color_map[PauliType.EMPTY] = QColor('white')
        
    def update_group_type_color_map(self, group_type:PauliType, color_str:str):
        color = QColor(color_str)
        self.group_type_color_map[group_type] = color
        for i in self.items():
            if isinstance(i, GaugeGroup):
                i.generate_polygon()
        
    def get_group_type_color(self, group_type: PauliType):
        print(f"group pault ti get is; {group_type}")
        print(f"colormaps is; {self.group_type_color_map}")
        return self.group_type_color_map[group_type]

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
        if event.key() == self.KEY_R:
            for i in self.items():
                if isinstance(i, GaugeGroup):
                    i.set_random_paulis()
                    
        if event.key() == self.KEY_C:
            gb = self.SelectGroupSectionTypeBox()
            gb.exec()
            group_type = gb.group_type
            print(f"coloring w group tupe: {group_type}")
            for item in self.selectedItems():
                if isinstance(item, GaugeGroup):
                    selected_qubits = []
                    for qubit in item.qubits:
                        if qubit.isSelected():
                            selected_qubits.append(qubit)
                    item.set_qubit_pauli(selected_qubits, group_type)
            self.clearSelection()
            


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
            gb = self.SelectGroupInfoBox()
            gb.exec()
            if gb.group_class == GroupType.GAUGEGROUP:
                group = GaugeGroup
            elif gb.group_class == GroupType.STABILIZER:
                group = Stabilizer
            else:
                # TODO add logging
                return
            gtype = gb.group_type
            qubits = []
            for i in self.selectedItems():
                if isinstance(i, Qubit):
                    qubits.append(i)
            group_item = group(qubits, gtype)
            self.addItem(group_item)
            for q in qubits:
                q.add_group(group_item)

        self.set_mode(self.MoveItem)
        for item in self.selectedItems():
            item.setSelected(False)

    class SelectGroupSectionTypeBox(QDialog):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.setWindowTitle("Choose Type")
            self._group_type = PauliType.EMPTY
        
            self._type_layout = QVBoxLayout()
            self.setLayout(self._type_layout)
        
            self._group_type_box = QComboBox()
            self._type_layout.addWidget(self._group_type_box)
        
            for gt in PauliType:
                print(f"fish sdlkf sdf GT is {gt}")
                # Adds an item to the combobox with the given text, and containing the specified userData (stored in the Qt::UserRole). The item is appended to the list of existing items.
                self._group_type_box.addItem(str(gt))
            self._group_type_box.currentTextChanged.connect(
                self.set_group_type
            )
            self._group_type_box.setCurrentText(str(PauliType.EMPTY))
        
            self._choice_layout = QHBoxLayout()
            self._type_layout.addLayout(self._choice_layout)
        
            self._okay_button = QPushButton("Okay")
            self._choice_layout.addWidget(self._okay_button)
            self._okay_button.clicked.connect(self.accept)

        @property
        def group_type(self):
            return self._group_type

        def set_group_type(self):
            self._group_type = PauliType(self._group_type_box.currentText())
            print(f"just set type only gt to {self._group_type} ")
            self.update()

    class SelectGroupInfoBox(QDialog):
    
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.setWindowTitle("Choose which group you wish to create")
            self._group_class = GroupType.UNSET
            self._group_type = PauliType.EMPTY
            self._group_info_layout = QHBoxLayout()
            self.setLayout(self._group_info_layout)
        
            self.group_type_layout = QVBoxLayout()
            self._group_info_layout.addLayout(self.group_type_layout)
        
            self._label = QLabel(
                f"Current group choice is: {str(self.group_class)}")
            self.group_type_layout.addWidget(self._label)
        
            self._gauge_button = QPushButton("Gauge Operator", parent=self)
            self._gauge_button.setChecked(False)
            self.group_type_layout.addWidget(self._gauge_button)
            self._gauge_button.clicked.connect(
                lambda: self.set_group(GroupType.GAUGEGROUP))
            self._gauge_button.setFocusPolicy(Qt.NoFocus)
        
            self._stabilizer_button = QPushButton("Stabilizer Operator",
                parent=self)
            self.group_type_layout.addWidget(self._stabilizer_button)
            self._stabilizer_button.clicked.connect(
                lambda: self.set_group(GroupType.STABILIZER))
            self._stabilizer_button.setFocusPolicy(Qt.NoFocus)
        
            self._group_type_box = QComboBox()
            self._group_info_layout.addWidget(self._group_type_box)
        
            for gt in PauliType:
                self._group_type_box.addItem(str(gt))
            self._group_type_box.currentTextChanged.connect(
                self.set_group_type
            )
            self._group_type_box.setCurrentText(str(PauliType.EMPTY))
        
            self._choice_layout = QHBoxLayout()
            self.group_type_layout.addLayout(self._choice_layout)
        
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
    
        @property
        def group_type(self):
            print(f"returning group tyupe: {self._group_type}")
            return self._group_type
    
        def set_group(self, value: GroupType):
            if not isinstance(value, GroupType):
                value = GroupType.INVALID
            self._group_class = value
            self.update()
    
        def set_group_type(self):
            self._group_type = PauliType(self._group_type_box.currentText())
            print(f"just set gt to {self._group_type} ")
            self.update()
            
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
