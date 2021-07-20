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
    
    def __init__(self, points, pauli_map:Dict[QPointF,PauliType]=None, pauli_def:PauliType=PauliType.EMPTY):

        super().__init__(points)
        
        self.setBrush(QBrush('blue'))
        self.setFillRule(Qt.OddEvenFill)
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        
        if pauli_map is None:
            self.broken = False
            self._error_groups = dict() # vertices -- Paulis
            
            if pauli_def:
                for p in points:
                    self.update_error_groups(p, pauli_def)
            else:
                for p in points:
                    self.update_error_groups(p,PauliType.EMPTY)
                
        else:
            self.broken = True
            self._error_groups = pauli_map

    def generate_polygon(self):
        # force paint event
        self.update(self.boundingRect())
        
    def update_error_groups(self, point:QPointF, value:PauliType):
        key_point = (point.x(), point.y())
        self._error_groups[key_point] = value
        
    def get_from_error_group(self, point: QPointF):
        key_point = (point.x(), point.y())
        return self._error_groups[key_point]

    def is_point_in_error_group(self, point: QPointF):
        key_point = (point.x(), point.y())
        return key_point in self._error_groups
        
    def set_random_paulis(self):
        self.broken = True
        for q in self._error_groups.keys():
            options = list(PauliType)
            options.remove(PauliType.EMPTY)
            c = random.choice(options)
            self._error_groups[q] = c
        self.generate_polygon()
        
    def setPolygon(self, polygon: Union[QPolygonF, Sequence[QPointF], QPolygon, QRectF]) -> None:
        
        # update Paulis for rotation
        cur_poly = self.polygon()
        num_vert = len(cur_poly)
        if cur_poly != 0:
            new_values = []
            for indx in range(num_vert):
                if self.is_point_in_error_group(cur_poly[indx]):
                    new_values.append((polygon[indx], self.get_from_error_group(cur_poly[indx])))
                else:
                    new_values.append((polygon[indx], PauliType.EMPTY))
            self._error_groups = dict()
            for tup in new_values:
                self.update_error_groups(tup[0], tup[1])
        
        super().setPolygon(polygon)
        
    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[ QWidget ] = ...) -> None:
    
        pen = self.pen()
        painter.setPen(pen)
    
        if not self.broken:
            
            painter.setBrush(QColor(self.scene().get_group_type_color(random.choice(list(self._error_groups.values())))))
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
                vertex = poly[indx]

                qcolor = QColor('salmon')
                qcolor.setAlpha(0)
                pen.setColor(qcolor)
                pauli = self.get_from_error_group(vertex)
                painter.setBrush(self.scene().get_group_type_color(pauli))
                painter.setPen(pen)
                
             
                p1_p = poly[ (indx + 1) % (len(poly)) ]
                hp1_p = QPointF((p1_p.x() + vertex.x()) / 2, (p1_p.y() + vertex.y()) / 2)
                p2_p = poly[ (indx - 1) % (len(poly)) ]
                hp2_p = QPointF((p2_p.x() + vertex.x()) / 2, (p2_p.y() + vertex.y()) / 2)
            
                poly = QPolygonF([ hp1_p, vertex, hp2_p, centroid ])
                painter.drawPolygon(poly)
    
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
            self.update_error_groups(key, pauli)
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


class DiagramScene(QGraphicsScene):
    RETURN = 16777220
    DELETE = 16777223
    BACKSPACE = 16777219
    KEY_1 = 49
    KEY_2 = 50
    KEY_3 = 51
    KEY_4 = 52
    KEY_C = 67
    KEY_R = 82
    KEY_J = 74
    KEY_A = 65

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
        
        if event.key() == self.KEY_R:
            for item in self.selectedItems():
                if isinstance(item, GaugeGroup):
                    cur_polygon = item.polygon()
                    item.setPolygon(QPolygonF(self._tiling.rotate_tile_around_origin(cur_polygon)))
            self._tiling.update(self._tiling.boundingRect())
        
        if event.key() == self.KEY_C:
            scene_polys = []
            scene_error_groups = dict()
            for item in self.selectedItems():
                if isinstance(item, GaugeGroup):
                    scene_polys.append((item.scenePos(), item.polygon()))
                    
                    for k in item._error_groups.keys():
                        ipos = item.pos()
                        polypos = k
                        scene_vert_pos = (polypos[0] + ipos.x(), polypos[1] + ipos.y())
                        scene_error_groups[scene_vert_pos] = item._error_groups[k]
                        
                    self.removeItem(item)
                    self._tiling.update(self._tiling.boundingRect())

             
            new_poly, new_pos = self._tiling.combine_tiles(scene_polys)
            poly_error_groups = dict()
            for k in scene_error_groups.keys():
                scene_pos = k
                polypos = (scene_pos[0] - new_pos.x(), scene_pos[1] - new_pos.y())
                poly_error_groups[polypos] = scene_error_groups[k]
            gg = GaugeGroup(QPolygonF(new_poly), poly_error_groups)
            gg.setPos(new_pos)
            self.addItem(gg)
        
        if event.key() == self.KEY_J:
            for item in self.selectedItems():
                if isinstance(item, GaugeGroup):
                    item.set_random_paulis()
                    
        if event.key() == self.KEY_A:
            options = list(PauliType)
            options.remove(PauliType.EMPTY)
            size = self._tiling.square_size
            for i in range(10):
                for j in range(3):
                    tile = self._tiling.make_tile(4)
                    gg = GaugeGroup(QPolygonF(tile), pauli_def=random.choice(options))
                    gg.setPos(size * i, size * j)
                    self.addItem(gg)


    def add_group(self, vertex_num: int):
        tile = GaugeGroup(self._tiling.make_tile(vertex_num))
        self.addItem(tile)
        
        