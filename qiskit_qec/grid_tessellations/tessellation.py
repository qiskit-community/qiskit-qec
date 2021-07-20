
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
import numpy as np

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
        pen.setWidth(10)
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
        
    def rotate_tile_around_origin(self, tile_points: QPolygonF, ninety_mult: int=1):
        if ninety_mult == 0:
            return tile_points
        if ninety_mult == 1:
            new_points = []
            x_mat = []
            x_sum = 0
            
            y_mat = []
            y_sum = 0
            
            count = 0
            for p in tile_points:
                count += 1
                x_sum += p.x()
                x_mat.append(p.x())
                
                
                y_mat.append(p.y())
                y_sum += p.y()
             
            P = np.array([x_mat, y_mat])
            centr_x = self.square_size/2
            centr_y = self.square_size/2
            
            C_top = np.full((1, count), centr_x)[0]
            C_bot = np.full((1, count), centr_y)[0]
            C = np.array([C_top, C_bot])
            
            R = [
                [0, 1],
                [-1, 0]
            ]
            
            P_new = np.matmul(R, (P - C)) + C
            
            for i in range(len(P_new[0])):
                new_points.append(QPointF(P_new[0][i], P_new[1][i]))
                
                print(f"""oldpoint: {P[0][i], P[1][i]}
                          tilepoint: {tile_points[i]}
                           newpoint: {new_points[i]}""")
            return new_points
        
    
    def perimeter_size(self):
        # TODO turn this into boundRect overwrite
        return (self.grid_width * self.square_size, self.grid_height * self.square_size)
    
        
            
        
            
        
        

        
        
    
    
    
    