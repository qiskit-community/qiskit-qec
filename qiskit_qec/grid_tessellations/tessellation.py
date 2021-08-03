
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
from PySide6.QtWidgets import (QStyleOptionGraphicsItem, QGraphicsPathItem,
    QApplication, QButtonGroup, QComboBox, QFontComboBox, QGraphicsAnchorLayout,
    QGraphicsItem, QGraphicsLineItem, QGraphicsPolygonItem, QGraphicsTextItem,
    QGraphicsEllipseItem, QGraphicsScene, QGraphicsView, QGridLayout,
    QHBoxLayout, QLabel, QMainWindow, QMenu, QMessageBox, QSizePolicy, QToolBox,
    QToolButton, QWidget, QPushButton, QVBoxLayout)
from enum import Enum
import uuid
import numpy as np



class GroupTile(QPainterPath):
    # literally just the outline of the path
    # has info needed to self-break because ....
    # have whatever info tessellation needs to rotate/edit object
    # TODO exists only for orientation/ordered_qubits for Tessellation since Tessellation
    # can't just give qubits in order for rendering since tessellation determines tile shape
    
    def __init__(self, qubits: [QPointF], orientation=0, *args, **kwargs):
        super(GroupTile, self).__init__(*args, **kwargs)
        self.ordered_qubits = qubits
        self.num_qubits = len(self.ordered_qubits)
        self.orientation = orientation
        # TODO self outline
        # TODO self broken shards
    # TODO is this too much? will i have to make 1 for each rotate?
    # not this just keep track of order of qubits
    # be the exact same Q..Item except w 1 array of order of qubits
    # vertices (qubits)
    # rotation
    
    def rotate_me(self, degree:int) -> ([QPointF], int):
        # new Qbits and orientation
        raise NotImplementedError
    

class SquareGridGroupTile(GroupTile):
    def __init__(self, qubits: [QPointF], orientation=0, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def break_me(self):
        pass


class StandardSquareGridGroupTile(SquareGridGroupTile):
    #  QPainterPath::addPolygon
    # closeSubpath()
    def __init__(self, qubits: [QPointF], two_qubit_orientation=0, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.qubits = qubits
        self.two_qubit_orientation = two_qubit_orientation
        
        # TODO generate self
    
    def break_me(self):
        if len(self.qubits) == 2:
            pass
            # disgusting
            # sharding depends on orientation and size of arc which is based off square size
        else:
            pass
            
    def rotate_me(self, rotation:int) -> ([QPointF], int):
        rotation = rotation%360
        if len(self.qubits) == 2:
            #disgusting
            if rotation == 180:
                if self.two_qubit_orientation == 0:
                    return self.qubits, 1
                else:
                    return self.qubits, 0
            if rotation == 90:
                pass
                
        # normal polygon rotation

class Tessellation(QGraphicsPolygonItem):
    #TODO should this be PainterPath?
    # tile == polygon
    # group == qubits
    """
    Designates "qubits" by creating a grid whose intersecting lines are qubits
    """
    def __init__(self):
        super().__init__()

    def paint(self):
        raise NotImplementedError
    
    def find_closest_qubit_location(self, point:QPointF):
        raise NotImplementedError
    
    def generate_tile(self, qubits: List[QPointF]) -> GroupTile:
        raise NotImplementedError

# FUNC
# q = qubit func
# t = tile func
# return q, t

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
        #TODO remove me?
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
            
    def find_closest_qubit_location(self, point:QPointF) -> QPointF:
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
   
    def generate_tile(self, qubits: List[QPointF]) -> GroupTile:
        pass
   
    def create_new_group_material(self, num_vertices) -> (List[QPointF], GroupTile):
        # returns qubits, painterpath
    
        # TODO handle < 2
        # TODO GroupTile instead of poly
        # TODO return List, QPainter
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
     
    def rotate_tile_around_origin(self, qubits: [QPointF], ninety_mult: int=1) -> GroupTile:
        # JUST rotates qubits and then chooses GroupTile Type
        # rotation always leads to new qubits, even when it doesn't
        #vertices of a Path are always qubits
        # standardized shapes for N qubits w vertices
        
        # rotations that change qubits
        # rotations that change orientation of path
        
        # rotate the orientation of the path
        # check
        qubits = self.rotate_qubit_group(qubits, ninety_mult)
        painter_path = self.generate_tile(qubits)
        return  painter_path
        
        
    
    
    def rotate_qubit_group(self, qubits: [QPointF], ninety_mult: int=1) -> [QPointF]:
        # TODO return qubits
        if ninety_mult == 0:
            return tile_points
        if ninety_mult == 1:
        
            new_points = [ ]
            x_mat = [ ]
            x_sum = 0
        
            y_mat = [ ]
            y_sum = 0
        
            count = 0
            for p in tile_points:
                count += 1
                x_sum += p.x()
                x_mat.append(p.x())
            
                y_mat.append(p.y())
                y_sum += p.y()
        
            P = np.array([ x_mat, y_mat ])
            poly_centr = self.calculate_centroid(tile_points)
            centr_x = poly_centr.x()
            centr_y = poly_centr.y()
        
            C_top = np.full((1, count), centr_x)[ 0 ]
            C_bot = np.full((1, count), centr_y)[ 0 ]
            C = np.array([ C_top, C_bot ])
        
            R = [
                [ 0, 1 ],
                [ -1, 0 ]
            ]
        
            P_new = np.matmul(R, (P - C)) + C
        
            for i in range(len(P_new[ 0 ])):
                valid_point = self.find_closest_qubit_location(QPointF(P_new[ 0 ][ i ], P_new[ 1 ][ i ]))
                new_points.append(valid_point)
            
                print(f"""oldpoint: {P[ 0 ][ i ], P[ 1 ][ i ]}
                          tilepoint: {tile_points[ i ]}
                           newpoint: {new_points[ i ]}""")
            return new_points


    def perimeter_size(self):
        # TODO turn this into boundRect overwrite
        return (self.grid_width * self.square_size, self.grid_height * self.square_size)
    
    
    
    def combine_tiles(self, polygon_list: [(QPointF, QPolygonF)]) -> ([QPointF], QPointF):
        # SHOULD NOT CHANGE SCENE POS OF VERTICES BY THE END (WILL MESS UP ERROR GROUPS)
        # covert polys to scene then combine to megapoly then convert back to poly
        min_x = sys.maxsize
        min_y = sys.maxsize
        scene_vertices = []
        scene_vert_set = set()
        for tup in polygon_list:
            poly = tup[1]
            pos = tup[0]
            if pos.x() < min_x:
                min_x = pos.x()
            if pos.y() < min_y:
                min_y = pos.y()
            for point in poly:
                pp = (pos.x() + point.x(), pos.y() + point.y())
                if pp not in scene_vert_set:
                    scene_vertices.append(QPointF(pos.x() + point.x(), pos.y() + point.y()))
                    scene_vert_set.add(pp)
            
            
        centroid = self.calculate_centroid(scene_vertices)
        # vertices will also be bubble sorted
        new_vertices = self.bubble_sort_vertices(scene_vertices, centroid)
        print(f"new points: {new_vertices}")
        
        new_poly_vertices = []
        for point in new_vertices:
            new_poly_vertices.append(QPointF(point.x() - min_x, point.y() - min_y))
        print(f"new poly points: {new_poly_vertices}")

        # subtract by min X and min Y?

        return new_poly_vertices, QPointF(min_x, min_y)
        
    def polygon_centroid(self, polygon:QPolygonF):
        x_sum = 0
        y_sum = 0
        count = 0
        for point in polygon:
            count += 1
            x_sum += point.x()
            y_sum += point.y()
        centroid = QPointF(x_sum / count, y_sum / count)
        return centroid
    
    def calculate_centroid(self, vertices:[QPointF]):
        x_sum = 0
        y_sum = 0
        count = 0
        for point in vertices:
                count += 1
                x_sum += point.x()
                y_sum += point.y()
        centroid = QPointF(x_sum/count, y_sum/count)
        return centroid
        

    def less(self, point_a :QPointF, point_b: QPointF, centroid: QPointF):
        if (point_a.x() - centroid.x() >=
            0) and (point_b.x() - centroid.x() < 0):
            return True
        if (point_a.x() - centroid.x() <
            0) and (point_b.x() - centroid.x() >= 0):
            return False
    
        if (point_a.x() - centroid.x()
            == 0) and (point_b.x() - centroid.x() == 0):
            if (point_a.y() - centroid.y() >=
                0) or (point_b.y() - centroid.y() >= 0):
                return point_a.y() > point_b.y()
            return point_b.y() > point_a.y()
    
        # compute the cross product of vectors(center -> a) x(center -> b) int
        det = (point_a.x() -
               centroid.x()) * (point_b.y() - centroid.y()) - (
                  point_b.x() - centroid.x()) * (point_a.y() -
                                                      centroid.y())
        if (det < 0):
            return True
        if (det > 0):
            return False
    
        # points a and b are on the same line from the center
    
        # check which point is closer to the center
        d1 = (point_a.x() - centroid.x()) * (point_a.x() - centroid.x(
        )) + (point_a.y() - centroid.y()) * (point_a.y() -
                                                  centroid.y())
        d2 = (point_b.x() - centroid.x()) * (point_b.x() - centroid.x(
        )) + (point_b.y() - centroid.y()) * (point_b.y() -
                                                  centroid.y())
        return d1 > d2


    def bubble_sort_vertices(self, vertices:[QPointF], centroid: QPointF):
        # bubble sort
        n = len(vertices)
    
        # Traverse through all array elements
        for i in range(n):
        
            # Last i elements are already in place
            for j in range(0, n - i - 1):
            
                # traverse the array from 0 to n-i-1
                # Swap if the element found is greater
                # than the next element
                if self.less(vertices[ j ], vertices[ j + 1 ], centroid):
                    vertices[ j ], vertices[ j + 1 ] = vertices[j + 1 ], vertices[ j ]
                    
        return vertices
        
        
    
    
    
    