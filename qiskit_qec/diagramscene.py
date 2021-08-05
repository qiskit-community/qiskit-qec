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

import random
from typing import Set

from PySide6.QtCore import (QLineF, Qt,
                            Signal)
from PySide6.QtGui import (QColor, QFont, QKeyEvent, QPolygonF)
from PySide6.QtWidgets import (QGraphicsItem, QGraphicsScene, QGraphicsSceneMouseEvent, QGraphicsTextItem, QMessageBox)

import qiskit_qec.diagramscene_rc as diagramscene_rc
from qiskit_qec.grid_tessellations.tessellation import Square

print(diagramscene_rc)

from qiskit_qec.graphics_library.graphics_ui_library import GaugeGroup, GaugeGroupFace, DiagramTextItem
from qiskit_qec.graphics_library.graphics_library import PauliType, SelectGroupSectionTypeBox



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
        
        self.all_gauge_groups = [] # this is bc some Gauges have multiple faces
        
        
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
                i.update_on_bounding_rect()
        
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
            # update error group --> map
            # update sharding
            for item in self.selectedItems():
                if isinstance(item, GaugeGroupFace):
                    # TODO check if belongs to multi-faced group and handle accordingly
                    cur_path = item.path_tile
                    item.setPath(self._tiling.rotate_tile_around_origin(cur_path))
            
 
        if event.key() == Qt.Key_C:
            selected_items = set()
            for item in self.selectedItems():
                if isinstance(item, GaugeGroupFace):
                    selected_items.add(item)
            if len(selected_items) > 1:
                self.combine_group_faces(selected_items)
             
        
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
                    tile = self._tiling.create_tile(4)
                    gg = GaugeGroup(QPolygonF(tile), pauli_def=random.choice(options))
                    gg.setPos(size * i, size * j)
                    self.addItem(gg)
        
        if event.key() == Qt.Key_P:
            sgstb = SelectGroupSectionTypeBox()
            sgstb.exec()
            pauli = sgstb.pauli_type
            for item in self.selectedItems():
                if isinstance(item, GaugeGroup):
                    item.set_entire_group_face_pauli(pauli)
            
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
                    
                    item.update_on_bounding_rect()
            
        self._tiling.update(self._tiling.boundingRect()) # TODO cleaner than calling this everywhere
        super().keyPressEvent(event)
        print("hi")
    
    def add_group(self, vertex_num: int):
        gauge_group_gi = GaugeGroup()
        self.all_gauge_groups.append(gauge_group_gi)
        gauge_group_path = self._tiling.create_tile(vertex_num)
        gauge_group_face = GaugeGroupFace(gauge_group_path)
        gauge_group_gi.set_faces([gauge_group_face])
        gauge_group_face.setPos(200,200)
        self.addItem(gauge_group_face)
        
    def combine_group_faces(self, items_to_combine:Set[GaugeGroupFace]):
        scene_path_tiles = []
        used_gauge_groups = []
        for item in items_to_combine:
            scene_path_tiles.append((item.scenePos(), item.path_tile)) # to give to tessellation to combine tiles
            used_gauge_groups.append(item.gauge_group) # will need to combine all these GaugeGroups into 1
            item.gauge_group.remove_faces([item]) # remove item from GaugeGroup
        
        # take all other polys from all gauge groups and combine into 1
        other_faces = []
        for gg in used_gauge_groups:
            other_faces += gg.get_faces() # TODO edge case: 2 different multi stabs (combined by other stab bits) share a face shape, how to avoid doubling?
            self.all_gauge_groups.remove(gg)
        
        new_gauge_group = GaugeGroup()
        self.all_gauge_groups.append(new_gauge_group)
        new_gauge_group.set_faces(other_faces)
        
        new_path_tile, scene_pos = self._tiling.combine_paths_scene(scene_path_tiles) # combine tiles into 1 tile
        
        for item in items_to_combine: # remove original tiles from scene
            self.removeItem(item)
            
        gauge_group_face = GaugeGroupFace(new_path_tile)
        gauge_group_face.setPos(scene_pos)
        
        new_gauge_group.add_face(gauge_group_face)
        
        self.addItem(gauge_group_face)
        
        





