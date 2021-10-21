
import math
import array
import sys
from typing import List, Optional

import random
from PySide6.QtWidgets import QDialog, QMainWindow, QStyle
from PySide6.QtCore import (QLineF, QPointF, QRect, QRectF, QSize, QSizeF, Qt,
                            Signal, QObject)
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
import uuid



class Tessellation(QGraphicsPolygonItem):
    def __init__(self):
        super().__init__()

    def paint(self):
        raise NotImplementedError
    
    def find_closest_vertex(self, point:QPointF):
        raise NotImplementedError
    
    def make_tile(self, num_vertices:int) -> QGraphicsItem:
        raise NotImplementedError


class Square(Tessellation):
    # grid -- num of x
    # num -- pixels
    def __init__(self, color='grey', square_size=50, grid_height=10, grid_width=10):
        super(Square, self).__init__()
        self.color = color
        self.square_size = square_size
        self.grid_height = grid_height
        self.grid_width = grid_width
        
    def itemChange(self, change, value):
        self.update(self.boundingRect())
        return super().itemChange(change, value)
        
  
    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[QWidget] = ...) -> None:
        pen = self.pen()
        pen.setColor(self.color)
        painter.setPen(pen)
        
        scene_rect = self.scene().sceneRect()
        scene_height = scene_rect.height()
        scene_width = scene_rect.width()
       
        num_height = self.square_size * self.grid_height
        num_width = self.square_size * self.grid_width

        
        for i in range(self.grid_height):
            yax = i * self.square_size
            painter.drawLine(QLineF(0,yax,num_width, yax))
        
        for j in range(self.grid_width):
            xax = j * self.square_size
            painter.drawLine(QLineF(xax, 0, xax,num_height))
            
    def find_closest_vertex(self, point:QPointF) -> QPointF:
        x = point.x()
        y = point.y()
        
        xrem = x % self.square_size
        if xrem/self.square_size > 0.5:
            closest_x = x-xrem+self.square_size
        else:
            closest_x = x-xrem
            
        yrem = y % self.square_size
        if yrem/self.square_size > 0.5:
            closest_y = y-yrem+self.square_size
        else:
            closest_y = y-yrem
            
        return QPointF(closest_x, closest_y)
    
    
    def make_tile(self, num_vertices:int) -> List[QPointF]:
        # TODO handle < 2
        # TODO automate num_vert to valid polygon?
        if num_vertices == 2:
            return [QPointF(0, self.square_size), QPointF(self.square_size, self.square_size)]
        
        elif num_vertices == 3:
            return[
                QPointF(0,0),
                QPointF(0, self.square_size),
                QPointF(self.square_size, 0),
            ]
        else:
            return [
                QPointF(0, 0),
                QPointF(0, self.square_size),
                QPointF(self.square_size, self.square_size),
                QPointF(self.square_size, 0),
            ]
    
    def perimeter_size(self):
        # TODO turn this into boundRect overwrite
        return (self.grid_width * self.square_size, self.grid_height * self.square_size)
    
        
            
        
            
        
        

        
        
    
    
    
    