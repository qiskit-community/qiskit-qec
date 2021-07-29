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

from qiskit_qec.graphics_library.graphics_library import GaugeGroup, GroupType,PauliType, DiagramTextItem, SelectGroupSectionTypeBox



class DiagramScene(QGraphicsScene):

    InsertItem, InsertLine, InsertText, MoveItem = range(4)
    HIGHLIGHT_COLOR = QColor('blue')

    text_inserted = Signal(QGraphicsTextItem)

    item_selected = Signal(QGraphicsItem)

    def __init__(self, itemMenu, parent=None, tiling=Square):
        super().__init__(parent)
        self.setup_pauli_type_color_map()
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
        
        
    def setup_pauli_type_color_map(self):
        self.group_type_color_map = {}
        # todo make random if becomes too many
        self.group_type_color_map[PauliType.X] = QColor('beige')
        self.group_type_color_map[PauliType.Y] = QColor('darkseagreen')
        self.group_type_color_map[ PauliType.Z ] = QColor('darksalmon')
        self.group_type_color_map[PauliType.EMPTY] = QColor('white')
        
    def update_pauli_type_color_map(self, group_type:PauliType, color_str:str):
        color = QColor(color_str)
        self.group_type_color_map[group_type] = color
        for i in self.items():
            if isinstance(i, GaugeGroup):
                i.generate_polygon()
        
    def get_pauli_type_color(self, group_type: PauliType):
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

        if event.key() == Qt.Key_Delete or event.key() == Qt.Key_Backspace:
            for i in self.selectedItems():
                self.removeItem(i)
                
        if event.key() == Qt.Key_2:
            self.add_group(2)
        elif event.key() == Qt.Key_3:
            self.add_group(3)
        elif event.key() == Qt.Key_4:
            self.add_group(4)
        
        if event.key() == Qt.Key_R:
            for item in self.selectedItems():
                if isinstance(item, GaugeGroup):
                    cur_polygon = item.polygon()
                    item.setPolygon(QPolygonF(self._tiling.rotate_tile_around_origin(cur_polygon)))
            self._tiling.update(self._tiling.boundingRect())
        
        if event.key() == Qt.Key_C:
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
        
        if event.key() == Qt.Key_J:
            for item in self.selectedItems():
                if isinstance(item, GaugeGroup):
                    item.set_random_paulis()
                    
        if event.key() == Qt.Key_A:
            options = list(PauliType)
            options.remove(PauliType.EMPTY)
            size = self._tiling.square_size
            for i in range(2, 9):
                for j in range(2, 6):
                    tile = self._tiling.make_tile(4)
                    gg = GaugeGroup(QPolygonF(tile), pauli_def=random.choice(options))
                    gg.setPos(size * i, size * j)
                    self.addItem(gg)
        
        if event.key() == Qt.Key_P:
            sgstb = SelectGroupSectionTypeBox()
            sgstb.exec()
            pauli = sgstb.pauli_type
            for item in self.selectedItems():
                if isinstance(item, GaugeGroup):
                    item.set_entire_group_pauli(pauli)
            
        if event.key() == Qt.Key_S:
            if len(self.selectedItems()) >1:
                QMessageBox(None, None, "Please only select 1 group").exec()
            else:
                item = self.selectedItems()[0]
                if isinstance(item, GaugeGroup):
                    for key_point in item.key_point_from_error_group():
                        sgstb = SelectGroupSectionTypeBox()
                        sgstb.set_label_point(key_point)
                        sgstb.exec()
                        pauli = sgstb.pauli_type
                        item.update_error_group_using_key_point(key_point, pauli)
                    
                    item.generate_polygon()
            

    def add_group(self, vertex_num: int):
        tile = GaugeGroup(self._tiling.make_tile(vertex_num))
        tile.setPos(200,200)
        self.addItem(tile)


