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
import random
import qiskit_qec.diagramscene_rc as diagramscene_rc
import uuid

print(diagramscene_rc)


class Qubit(QGraphicsEllipseItem):
    size = 50
    name = "QUBIT"
    
    def __init__(self, contextMenu, parent=None, x_pos=0, y_pos=0):
        super().__init__(x_pos, y_pos, self.size, self.size, parent)
        
        self.groups = [ ]
        
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
            gr.update_on_bounding_rect()
        return super().itemChange(change, value)
    
    def get_center(self):
        upper_left = self.scenePos()
        upper_left.setX(upper_left.x() + self.size / 2)
        upper_left.setY(upper_left.y() + self.size / 2)
        return upper_left
    
    
    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[ QWidget ] = ...) -> None:
        super().paint(painter, option, widget)
        if self.isSelected():
            pen = QPen(self.scene().HIGHLIGHT_COLOR)
            
            painter.setPen(pen)
            painter.pen().setWidth(2)
            trans = QColor("white")
            trans.setAlpha(0)
            painter.setBrush(trans)
            painter.drawPolygon(self.boundingRect())


class PauliType():
    pass

class GaugeBoi():
    def __init__(self, qubits: list[ Qubit ], group_type: PauliType = None):
        self._error_groups = {}
        self.setup_qubits_and_paulis(qubits, group_type)


    def setup_qubits_and_paulis(self, qubits, group_type: PauliType):
        self.broken = False
        self.qubits = [ ]
        self.cur_points = [ ]
        self.add_qubits(qubits)
        self._group_paulies = set()
        if group_type is not None:
            self.set_entire_group_face_pauli(group_type)


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

    def remove_qubit(self, qubit):
        if qubit in self.qubits:
            self.qubits.remove(qubit)  # must come first to avoid recursion
        
            qubit.remove_group(self)  # TODO switch to signals
            self.bubble_sort_qubits()
            self.update_on_bounding_rect()
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
    
        self.update_on_bounding_rect()

        self.update_on_bounding_rect()


    def setup_qubits_and_paulis(self, qubits, group_type: PauliType):
        self.broken = False
        self.qubits = [ ]
        self.cur_points = [ ]
        self.add_qubits(qubits)
        self._group_paulies = set()
        if group_type is not None:
            self.set_entire_group_face_pauli(group_type)


    def update_on_bounding_rect(self) -> [ QPointF ]:
        self.cur_points = [ ]
        for q in self.qubits:
            self.cur_points.append(q.get_center())
    
        # todo how to properly call paint
        poly = QPolygonF(self.cur_points)
        self.setPolygon(poly)


    def set_entire_group_face_pauli(self, pauli: PauliType = PauliType.X):
        self.broken = False
        for key in self._error_groups.keys():
            self._error_groups[ key ] = pauli
        self.update_on_bounding_rect()


    def set_qubit_pauli(self, qubits: Qubit, pauli: PauliType = PauliType.X):
        self.broken = True
        for qu in qubits:
            self._error_groups[ qu.uuid ] = pauli
        self.update_on_bounding_rect()

    def set_random_paulis(self):
        for q in self._error_groups.keys():
            options = list(PauliType)
            options.remove(PauliType.EMPTY)
            c = random.choice(options)
        self._error_groups[ q ] = c
        self.update_on_bounding_rect()

    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[ QWidget ] = ...) -> None:
    
        pen = self.pen()
        painter.setPen(pen)
    
        if not self.broken:
            painter.setBrush(QColor(self.scene().get_pauli_type_color(self._error_groups[ self.qubits[ 0 ].uuid ])))
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
                painter.setBrush(self.scene().get_pauli_type_color(pauli))
            
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
            
            
class DiagramSceneBoi():
    item_inserted = Signal(Qubit)

    def __init__(self):
        pass
    
    def set_item_color(self, color):
        self._my_item_color = color
        if self.is_item_change(Qubit):
            item = self.selectedItems()[0]
            item.setBrush(self._my_item_color)
            
            
    def mousePressEvent(self, mouseEvent):
        if (mouseEvent.button() != Qt.LeftButton):
            return

        if self._my_mode == self.InsertItem:
            item = Qubit(self._my_item_menu)
            item.setBrush(self._my_item_color)
            self.addItem(item)
            item.setPos(mouseEvent.scenePos())
            self.item_inserted.emit(item)


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
                        selected_qubits = [ ]
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
