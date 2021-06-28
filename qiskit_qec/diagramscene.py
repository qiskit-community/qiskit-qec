
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
import sys

from PySide6.QtCore import (QLineF, QPointF, QRect, QRectF, QSize, QSizeF, Qt,
                            Signal)
from PySide6.QtGui import (QAction, QColor, QFont, QIcon, QIntValidator,
                           QPainter, QPainterPath, QPen, QPixmap, QPolygonF, QBrush, QKeyEvent)
from PySide6.QtWidgets import (QApplication, QButtonGroup, QComboBox,
                               QFontComboBox, QGraphicsAnchorLayout,
                               QGraphicsItem, QGraphicsLineItem,
                               QGraphicsPolygonItem, QGraphicsTextItem,
                               QGraphicsScene, QGraphicsView, QGridLayout,
                               QHBoxLayout, QLabel, QMainWindow, QMenu,
                               QMessageBox, QSizePolicy, QToolBox, QToolButton,
                               QWidget)

import diagramscene_rc


class Arrow(QGraphicsLineItem):
    def __init__(self, startItem, endItem, parent=None, scene=None):
        super().__init__(parent, scene)

        self._arrow_head = QPolygonF()

        self._my_start_item = startItem
        self._my_end_item = endItem
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self._my_color = Qt.black
        self.setPen(QPen(self._my_color, 2, Qt.SolidLine,
                Qt.RoundCap, Qt.RoundJoin))

    def set_color(self, color):
        self._my_color = color

    def start_item(self):
        return self._my_start_item

    def end_item(self):
        return self._my_end_item

    def boundingRect(self):
        extra = (self.pen().width() + 20) / 2.0
        p1 = self.line().p1()
        p2 = self.line().p2()
        rect = QRectF(p1, QSizeF(p2.x() - p1.x(), p2.y() - p1.y()))
        return rect.normalized().adjusted(-extra, -extra, extra, extra)

    def shape(self):
        path = super(Arrow, self).shape()
        path.addPolygon(self._arrow_head)
        return path

    def update_position(self):
        start = self.mapFromItem(self._my_start_item, 0, 0)
        end = self.mapFromItem(self._my_end_item, 0, 0)
        self.setLine(QLineF(start, end))

    def paint(self, painter, option, widget=None):
        if (self._my_start_item.collidesWithItem(self._my_end_item)):
            return

        my_start_item = self._my_start_item
        my_end_item = self._my_end_item
        my_color = self._my_color
        my_pen = self.pen()
        my_pen.setColor(self._my_color)
        arrow_size = 20.0
        painter.setPen(my_pen)
        painter.setBrush(self._my_color)

        center_line = QLineF(my_start_item.pos(), my_end_item.pos())
        end_polygon = my_end_item.polygon()
        p1 = end_polygon.at(0) + my_end_item.pos()

        intersect_point = QPointF()
        for i in end_polygon:
            p2 = i + my_end_item.pos()
            poly_line = QLineF(p1, p2)
            intersectType, intersect_point = poly_line.intersects(center_line)
            if intersectType == QLineF.BoundedIntersection:
                break
            p1 = p2

        self.setLine(QLineF(intersect_point, my_start_item.pos()))
        line = self.line()

        angle = math.acos(line.dx() / line.length())
        if line.dy() >= 0:
            angle = (math.pi * 2.0) - angle

        arrow_head1 = QPointF(math.sin(angle + math.pi / 3.0) * arrow_size,
                              math.cos(angle + math.pi / 3) * arrow_size)
        arrow_p1 = line.p1() + arrow_head1
        arrow_head2 = QPointF(math.sin(angle + math.pi - math.pi / 3.0) * arrow_size,
                              math.cos(angle + math.pi - math.pi / 3.0) * arrow_size)
        arrow_p2 = line.p1() + arrow_head2

        self._arrow_head.clear()
        for point in [line.p1(), arrow_p1, arrow_p2]:
            self._arrow_head.append(point)

        painter.drawLine(line)
        painter.drawPolygon(self._arrow_head)
        if self.isSelected():
            painter.setPen(QPen(my_color, 1, Qt.DashLine))
            my_line = QLineF(line)
            my_line.translate(0, 4.0)
            painter.drawLine(my_line)
            my_line.translate(0, -8.0)
            painter.drawLine(my_line)


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

class DiagramItem(QGraphicsPolygonItem):
    Step, Conditional, StartEnd, Io = range(4)

    def __init__(self, diagram_type, contextMenu, parent=None, scene=None):
        super().__init__(parent, scene)

        self.arrows = []
        self.stabs = []

        self.diagram_type = diagram_type
        self._my_context_menu = contextMenu

        path = QPainterPath()
        if self.diagram_type == self.StartEnd:
            path.moveTo(200, 50)
            path.arcTo(150, 0, 50, 50, 0, 90)
            path.arcTo(50, 0, 50, 50, 90, 90)
            path.arcTo(50, 50, 50, 50, 180, 90)
            path.arcTo(150, 50, 50, 50, 270, 90)
            path.lineTo(200, 25)
            self._my_polygon = path.toFillPolygon()
        elif self.diagram_type == self.Conditional:
            self._my_polygon = QPolygonF([
                    QPointF(-100, 0), QPointF(0, 100),
                    QPointF(100, 0), QPointF(0, -100),
                    QPointF(-100, 0)])
        elif self.diagram_type == self.Step:
            self._my_polygon = QPolygonF([
                    QPointF(-100, -100), QPointF(100, -100),
                    QPointF(100, 100), QPointF(-100, 100),
                    QPointF(-100, -100)])
        else:
            self._my_polygon = QPolygonF([
                    QPointF(-120, -80), QPointF(-70, 80),
                    QPointF(120, 80), QPointF(70, -80),
                    QPointF(-120, -80)])

        self.setPolygon(self._my_polygon)
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)

    def remove_arrow(self, arrow):
        try:
            self.arrows.remove(arrow)
        except ValueError:
            pass

    def remove_arrows(self):
        for arrow in self.arrows[:]:
            arrow.start_item().remove_arrow(arrow)
            arrow.end_item().remove_arrow(arrow)
            self.scene().removeItem(arrow)

    def add_arrow(self, arrow):
        self.arrows.append(arrow)

    def add_stab(self, stab):
        self.stabs.append(stab)

    # TODO deal with deleting stabs


    def image(self):
        pixmap = QPixmap(250, 250)
        pixmap.fill(Qt.transparent)
        painter = QPainter(pixmap)
        painter.setPen(QPen(Qt.black, 8))
        painter.translate(125, 125)
        painter.drawPolyline(self._my_polygon)
        return pixmap

    def contextMenuEvent(self, event):
        self.scene().clearSelection()
        self.setSelected(True)
        self._my_context_menu.exec(event.screenPos())

    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemPositionChange:
            for arrow in self.arrows:
                arrow.updatePosition()

        return value
    def get_center(self):
        return self.scenePos()



class Stabilizer(QGraphicsPolygonItem):
    def __init__(self, point_list:list[QPointF], qubits:list):
        super().__init__(point_list)
        self.point_list = point_list
        self.qubits = qubits
        self.qbrush = QBrush(QColor('red'), Qt.SolidPattern)
        self.setBrush(self.qbrush)
        self.pen = QPen()
        self.pen.setColor(QColor('red'))
        self.setPen(self.pen)
        # TODO delete stabs



class DiagramScene(QGraphicsScene):
    RETURN = 16777220
    InsertItem, InsertLine, InsertText, MoveItem = range(4)

    item_inserted = Signal(DiagramItem)

    text_inserted = Signal(QGraphicsTextItem)

    item_selected = Signal(QGraphicsItem)

    def __init__(self, itemMenu, parent=None):
        super().__init__(parent)

        self._my_item_menu = itemMenu
        self._my_mode = self.MoveItem
        self._my_item_type = DiagramItem.Step
        self.line = None
        self._text_item = None
        self._my_item_color = Qt.white
        self._my_text_color = Qt.black
        self._my_line_color = Qt.black
        self._my_font = QFont()

    def set_line_color(self, color):
        self._my_line_color = color
        if self.is_item_change(Arrow):
            item = self.selectedItems()[0]
            item.set_color(self._my_line_color)
            self.update()

    def set_text_color(self, color):
        self._my_text_color = color
        if self.is_item_change(DiagramTextItem):
            item = self.selectedItems()[0]
            item.setDefaultTextColor(self._my_text_color)

    def set_item_color(self, color):
        self._my_item_color = color
        if self.is_item_change(DiagramItem):
            item = self.selectedItems()[0]
            item.setBrush(self._my_item_color)

    def set_font(self, font):
        self._my_font = font
        if self.is_item_change(DiagramTextItem):
            item = self.selectedItems()[0]
            item.setFont(self._my_font)

    def set_mode(self, mode):
        self._my_mode = mode

    def set_item_type(self, type):
        self._my_item_type = type

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
            item = DiagramItem(self._my_item_type, self._my_item_menu)
            item.setBrush(self._my_item_color)
            self.addItem(item)
            item.setPos(mouseEvent.scenePos())
            self.item_inserted.emit(item)
        elif self._my_mode == self.InsertLine:
            self.line = QGraphicsLineItem(QLineF(mouseEvent.scenePos(),
                                        mouseEvent.scenePos()))
            self.line.setPen(QPen(self._my_line_color, 2))
            self.addItem(self.line)
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
        if self.line and self._my_mode == self.InsertLine:
            start_items = self.items(self.line.line().p1())
            if len(start_items) and start_items[0] == self.line:
                start_items.pop(0)
            end_items = self.items(self.line.line().p2())
            if len(end_items) and end_items[0] == self.line:
                end_items.pop(0)

            self.removeItem(self.line)
            self.line = None

            if (len(start_items) and len(end_items) and
                    isinstance(start_items[0], DiagramItem) and
                    isinstance(end_items[0], DiagramItem) and
                    start_items[0] != end_items[0]):
                start_item = start_items[0]
                end_item = end_items[0]
                arrow = Arrow(start_item, end_item)
                arrow.set_color(self._my_line_color)
                start_item.add_arrow(arrow)
                end_item.add_arrow(arrow)
                arrow.setZValue(-1000.0)
                self.addItem(arrow)
                arrow.update_position()

        self.line = None
        super(DiagramScene, self).mouseReleaseEvent(mouseEvent)

    def is_item_change(self, type):
        for item in self.selectedItems():
            if isinstance(item, type):
                return True
        return False

    def keyPressEvent(self, event: QKeyEvent) -> None:

        if event.key() == self.RETURN:
            print("return pressed")
            points = []
            qubits = []
            for i in self.selectedItems():
                qubits.append(i)
                points.append(i.get_center())
                print(f"points are {points}")
                stab = Stabilizer(points, qubits)
                self.addItem(stab)
            for q in qubits:
                q.add_stab(stab)


