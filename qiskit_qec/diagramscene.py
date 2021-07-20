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
    QToolButton, QWidget, QPushButton, QVBoxLayout, QGraphicsSceneMouseEvent)
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



class GaugeGroup(QGraphicsPolygonItem):
    name = "Gauge Group"
    
    def __init__(self, points):

        super().__init__(points)
        self.setBrush(QBrush('blue'))
        self.setFillRule(Qt.OddEvenFill)
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        
    def mouseReleaseEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        current_poly = self.polygon()
        print(f"current poly: {current_poly} and pos: {self.pos()} and scenepos: {self.scenePos()}")
        new_point = self.scene()._tiling.find_closest_vertex(self.scenePos())
        self.setPos(new_point)
        super().mouseReleaseEvent(event)
        
        

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
    R_KEY = 82

    InsertItem, InsertLine, InsertText, MoveItem = range(4)
    HIGHLIGHT_COLOR = QColor('blue')

    text_inserted = Signal(QGraphicsTextItem)

    item_selected = Signal(QGraphicsItem)

    def __init__(self, itemMenu, parent=None, tiling=Square):
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
        self._tiling = tiling()
        self.addItem(self._tiling)
        
        
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

        if self._my_mode == self.InsertText:
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
            self._tiling.update(self.sceneRect())
            super(DiagramScene, self).mouseMoveEvent(mouseEvent)

    def mouseReleaseEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        if self.mouseGrabberItem() is not None and self.mouseGrabberItem() != 0:
            self._tiling.update(self.sceneRect())
        super(DiagramScene, self).mouseReleaseEvent(event)
            
        

    def is_item_change(self, type):
        for item in self.selectedItems():
            if isinstance(item, type):
                return True
        return False

    def keyPressEvent(self, event: QKeyEvent) -> None:

        if event.key() == self.DELETE or event.key() == self.BACKSPACE:
            for i in self.selectedItems():
                self.removeItem(i)
                
        if event.key() == self.KEY_2:
            self.add_group(2)
        elif event.key() == self.KEY_3:
            self.add_group(3)
        elif event.key() == self.KEY_4:
            self.add_group(4)
        
        if event.key() == self.R_KEY:
            for item in self.selectedItems():
                if isinstance(item, GaugeGroup):
                    cur_polygon = item.polygon()
                    item.setPolygon(QPolygonF(self._tiling.rotate_tile_around_origin(cur_polygon)))
            
    def add_group(self, vertex_num: int):
        tile = GaugeGroup(self._tiling.make_tile(vertex_num))
        self.addItem(tile)