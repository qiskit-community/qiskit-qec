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
from typing import Dict, Optional

from PySide6.QtCore import (QPointF, Qt,
                            Signal)
from PySide6.QtGui import (QColor, QPainter, QPen)
from PySide6.QtWidgets import (QGraphicsItem, QGraphicsPathItem, QGraphicsSceneMouseEvent, QGraphicsTextItem,
                               QStyleOptionGraphicsItem, QWidget)

import qiskit_qec.diagramscene_rc as diagramscene_rc
from qiskit_qec.graphics_library.graphics_library import PauliType
from qiskit_qec.grid_tessellations.tessellation import GroupTile, Square

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


class GaugeGroupFace(QGraphicsPathItem):
    name = "Gauge Group"
    
    def __init__(self, group_tile: GroupTile, gauge_group=None, pauli_map: Dict[ QPointF, PauliType ] = None, pauli_def: PauliType = PauliType.EMPTY):
        
        super().__init__()
        self.gauge_group = gauge_group
        self.path_tile = group_tile
        # TODO make sure error group and group tile don't go out of sync
        
        self.setPath(group_tile)
        
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self.setFlag(QGraphicsItem.ItemSendsGeometryChanges, True)
    
    def set_gauge_group(self, gg):
        self.gauge_group = gg
    
    def setPath(self, path: GroupTile) -> None: # grouptile
        self.path_tile = path
        
        super().setPath(path)
        # getting shards of new path and adding any old Paulis to the new error groups
        if self.gauge_group is not None:
            self.gauge_group.update_qubits()

   
        
    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[ QWidget ] = ...) -> None:
        
        pen = self.pen()
        painter.setPen(pen)
        
        # TODO figure out where to decide polygon shape
        
        path_tile = self.path_tile
        if not path_tile.is_sharded:
            painter.setBrush(
                QColor(self.scene().get_pauli_type_color(path_tile.get_some_pauli_in_use()))
            )
            painter.drawPath(self.path())
            
        else:
            outline_color = QColor('black')
            painter.setPen(outline_color)

            fill_in_color = QColor('white')
            fill_in_color.setAlpha(255)
            painter.setBrush(fill_in_color)
           
            painter.drawPath(self.path())
            
            for shard_vertex in self.path_tile.sharding.keys():
                painter.setBrush(self.scene().get_pauli_type_color(self.path_tile.get_pauli_for_qubit( (shard_vertex))))
                painter.drawPath(self.path_tile.sharding[shard_vertex])
            
        if self.isSelected():
            #TODO figureo out why only 3/4 of square gets highlighted
            pen = QPen(self.scene().HIGHLIGHT_COLOR)
            painter.setPen(pen)
            painter.pen().setWidth(5)
            trans = QColor("white")
            trans.setAlpha(0)
            painter.setBrush(trans)
            painter.drawPath(self.path())
 
  
    def set_entire_group_face_pauli(self, pauli: PauliType = PauliType.X):
        self.path().set_entire_path_pauli(pauli)


    def mouseReleaseEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        current_path = self.path()
        
        # TODO rip _tiling being private. How to best fix?
        new_point = self.scene()._tiling.find_closest_valid_scene_position(self.scenePos())
        self.setPos(new_point)
        super().mouseReleaseEvent(event)

        
    def update_on_bounding_rect(self):
        # force paint event
        self.update(self.boundingRect())
    
    # TODO in order for multiple faces to act as 1 stabilizer, GaugeGroupFace
    # needs to take input and pass it upstream to GaugeGroup
    # up push selections to GaugeGroup
    # up push dragging to Gauge Group
    # on rotation how to push up to gauge group?


class GaugeGroup():
    
    def __init__(self, faces: [ GaugeGroupFace ] = None):
        if faces is None:
            self.group_qubits = [ ]
            self.group_faces = None
            return
            
        self.group_qubits = []
        for f in faces:
            self.group_qubits = self.group_qubits + f.ordered_qubits()
        self.group_faces = faces

    def add_face(self, face: GaugeGroupFace):
        if self.group_faces is None:
            self.group_faces = [face]
        else:
            self.group_faces.append(face)
        self.update_qubits()
    def set_faces(self, faces: [GaugeGroupFace]):
        self.group_qubits = []
        self.group_faces = faces
        for f in faces:
            f.set_gauge_group(self)
        self.update_qubits()

    def update_qubits(self):
        self.group_qubits = [ ]
        for f in self.group_faces:
            self.group_qubits = self.group_qubits + f.path_tile.ordered_qubits


    def set_entire_group_pauli(self, pauli: PauliType = PauliType.X):
        for f in self.group_faces:
            f.set_entire_group_face_pauli(pauli)
            


class Stabilizer(GaugeGroup):
    name = "Stabilizer"

    def __init__(self, qubits: list):
        super().__init__(qubits)

