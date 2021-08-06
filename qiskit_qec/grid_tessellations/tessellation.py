import sys
from typing import List, Optional

import random
from PySide6.QtCore import (QLineF, QPointF)
from PySide6.QtGui import (QPainter, QPainterPath, QPolygonF)
from PySide6.QtWidgets import (QStyleOptionGraphicsItem, QGraphicsPolygonItem, QWidget)
import numpy as np
from qiskit_qec.graphics_library.graphics_library import PauliType
import random
import sys
from typing import List, Optional

import numpy as np
from PySide6.QtCore import (QLineF, QPointF)
from PySide6.QtGui import (QPainter, QPainterPath, QPolygonF)
from PySide6.QtWidgets import (QGraphicsPolygonItem, QStyleOptionGraphicsItem, QWidget)

from qiskit_qec.graphics_library.graphics_library import PauliType


class GroupTile(QPainterPath):
    # literally just the outline of the path
    # has info needed to self-break because ....
    # have whatever info tessellation needs to rotate/edit object
    # TODO exists only for orientation/ordered_qubits for Tessellation since Tessellation
    # can't just give qubits in order for rendering since tessellation determines tile shape
    
    def __init__(self, qubits: [QPointF], is_sharded:bool=True):
        print("Group Tile thing ")
        super().__init__()
        self.ordered_qubits = qubits
        print("ordering quibts super: ", self.ordered_qubits)
        self.is_sharded = is_sharded


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
       
    @property
    def sharding(self):
        raise NotImplementedError

    def set_entire_tile_pauli(self, pauli: PauliType):
        raise NotImplementedError


    def update_pauli_map(self, point: QPointF, value: PauliType, is_sharded:bool):
        raise NotImplementedError

    def get_pauli_for_qubit(self, point):
        raise NotImplementedError

    def key_point_from_pauli_map(self):
        raise NotImplementedError


    def is_qubit_in_pauli_map(self, point):
        raise NotImplementedError

    def set_random_paulis(self):
        raise NotImplementedError

    def get_some_pauli_in_use(self):
        raise NotImplementedError
    
    def set_entire_tile_pauli(self):
        raise NotImplementedError
    def get_all_qubits_in_pauli_map(self):
        raise NotImplementedError


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
    

# FUNC
# q = qubit func
# t = tile func
# return q, t

class Square(Tessellation):
    # TODO qubits are always ordered based on less algorithm shown here. The ordering of the qubits is what allows for
    # rotation communication between Group and Tessellation
    # grid -- num of x
    # num -- pixels

    class SquareGridGroupTile(GroupTile):
        #  QPainterPath::addPolygon
        # closeSubpath()
        SHARD_PATH = 0
        SHARD_COLOR = 1
        POLYGON = 0
        ARC = 1


        def __init__(self, qubits: [ QPointF ], two_qubit_orientation=None, tile_type=None, pauli_map=None, is_sharded=False):
            super().__init__(qubits, is_sharded)
            # qubits are always sorted counter_clockwise unless other denoted
            print("creating squre grid tile")
            print("sqaregridgroup rodering qubits: ", self.ordered_qubits)

            if tile_type is None:
                self.tile_type = self.POLYGON
            else:
                self.tile_type = tile_type
                
            if self.tile_type == self.POLYGON:
                self.addPolygon(QPolygonF(self.ordered_qubits))
            else:
                if self.tile_type == self.ARC:
                    # closeSubpath?
                    self.two_qubit_orientation = two_qubit_orientation
                    qubit0 = qubits[0]
                    qubit1 = qubits[1]

                    arc = QPainterPath(qubit0)
                    w = qubit1.x() - qubit0.x()
                    h = qubit1.y() - qubit0.y()
                    

                    
                    
                    sz = h if h > w else w
                    
                    start_angle = 0
                    angle = 180
                    if self.two_qubit_orientation == 1:
                        angle = angle * -1
                        
                    start_point = QPointF(qubit0.x(), qubit0.y())
                    if w == 0:
                        start_point.setX(start_point.x() - sz/2)
                        start_angle = 90
                        
                    if h == 0:
                        start_point.setY(start_point.y() - sz/2)
                        start_angle = 0
                        
                    if w != 0 and h != 0:
                        start_point.setY(start_point.y() - sz/4)
                        start_point.setX(start_point.x() - sz/4)
                        start_angle = 90 + 45
                        
                    arc.arcTo(start_point.x(), start_point.y(), sz, sz, start_angle, angle)
                    
                    self.arc_rect_h = h
                    self.arc_rect_w = w
                    self.start_point = start_point
                    
                    self.addPath(arc)
                    
                    
            self.is_sharded = is_sharded
           
           
            # never shard until called
            self.sharding = None
            if pauli_map is not None:
                self._pauli_map = pauli_map
            else:
                self._pauli_map = dict()
                for qu in self.ordered_qubits:
                    self.update_pauli_map(qu, PauliType.EMPTY)
        #TODO _pauli_map into a proper property
                    
     
        @property
        def sharding(self):
            self.is_sharded = True
            if self._sharding is not None:
                return self._sharding
            self._sharding = dict()
            if len(self.ordered_qubits) == 2:
                pass
                # disgusting
                # sharding depends on orientation and size of arc which is based off square size
            elif self.tile_type == self.POLYGON:
                x_sum = 0
                y_sum = 0
                count = 0
                for point in self.ordered_qubits:
                    count += 1
                    x_sum += point.x()
                    y_sum += point.y()
                centroid = QPointF(x_sum / count, y_sum / count)
        
                for indx in range(len(self.ordered_qubits)):
                    vertex = self.ordered_qubits[ indx ]
                    p1_p = self.ordered_qubits[ (indx + 1) % (len(self.ordered_qubits)) ]
                    hp1_p = QPointF((p1_p.x() + vertex.x()) / 2, (p1_p.y() + vertex.y()) / 2)
                    p2_p = self.ordered_qubits[ (indx - 1) % (len(self.ordered_qubits)) ]
                    hp2_p = QPointF((p2_p.x() + vertex.x()) / 2, (p2_p.y() + vertex.y()) / 2)
                    
                    vertex_point = (vertex.x(), vertex.y())
                    tmp_path = QPainterPath()
                    tmp_path.addPolygon([ hp1_p, vertex, hp2_p, centroid ])
                    self.sharding[ vertex_point ] = tmp_path
                return self._sharding
            
            else:
                raise NotImplementedError
            
        @sharding.setter
        def sharding(self, value):
            self._sharding = value


        def set_entire_tile_pauli(self, pauli: PauliType):
            self.sharding = None
            self.is_sharded = False
            for qu in self._pauli_map:
                self.update_pauli_map(qu, pauli)
        

        def update_pauli_map(self, point: QPointF, value: PauliType, is_sharded=None):
            if is_sharded is not None:
                self.is_sharded = is_sharded
            if isinstance(point, QPointF):
                key_point = (point.x(), point.y())
            else:
                key_point = point
            self._pauli_map[ key_point ] = value

        def get_pauli_for_qubit(self, point):
            if isinstance(point, QPointF):
                key_point = (point.x(), point.y())
            else: 
                key_point = point
            if key_point in self._pauli_map:
                return self._pauli_map[ key_point ]
            return None

        def key_point_from_pauli_map(self):
            key_point_list = self._pauli_map.keys()
            return key_point_list


        def is_qubit_in_pauli_map(self, point):
            if isinstance(point, QPointF):
                key_point = (point.x(), point.y())
            else:
                key_point = point
            return key_point in self._error_groups

        def set_random_paulis(self):
            self.is_sharded = True
            for q in self._pauli_map.keys():
                options = list(PauliType)
                options.remove(PauliType.EMPTY)
                c = random.choice(options)
                self.update_pauli_map(q, c)
            
                
        def get_some_pauli_in_use(self):
            return random.choice(list(self._pauli_map.values()))

        def get_all_qubits_in_pauli_map(self):
            return self.ordered_qubits
            
            
    def __init__(self, color='grey', square_size=50, grid_height=10, grid_width=10):
        
        super(Square, self).__init__()
        self.square_size = square_size
        self.color = color
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
   
    def create_tile(self, num_vertices:int=4, two_orientation=None, two_magnitude=0):
        vertices = self.create_new_group_material(num_vertices, two_magnitude=two_magnitude)
        tile = self.convert_qubits_to_tile(vertices, two_qubit_orientation=two_orientation)
        return tile
       
    def convert_qubits_to_tile(self, qubits: [QPointF], two_qubit_orientation=None, tile_type=None, pauli_map=None, is_sharded=False):
        if len(qubits) > 2:
            if tile_type is None:
                tile_type = self.SquareGridGroupTile.POLYGON
            tile = self.SquareGridGroupTile(qubits, tile_type=tile_type, pauli_map=pauli_map, is_sharded=is_sharded)
            return tile
        else:
            tile = self.SquareGridGroupTile(qubits, tile_type=self.SquareGridGroupTile.ARC ,two_qubit_orientation=two_qubit_orientation)
            return tile


    
    def create_new_group_material(self, num_vertices: int, two_magnitude=0) -> [QPointF]:
        # returns qubits, painterpath
    
        # TODO handle < 2
        # TODO GroupTile instead of poly
        # TODO return List, QPainter
        if num_vertices == 2:
            if two_magnitude == 0:
                return [QPointF(0, self.square_size), QPointF(self.square_size, self.square_size)]
            elif two_magnitude == 1:
                return [QPointF(0, 0), QPointF(0, self.square_size)]
            else:
                raise NotImplementedError
                
        
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
     
    
    
    def rotate_tile_around_origin(self, group_tile: SquareGridGroupTile,  ninety_mult: int=1) -> GroupTile:
        # TODO make sure error group stays okay order qubits and reorder
        if ninety_mult == 0:
            return group_tile
        
        if len(group_tile.ordered_qubits) == 2:
            if ninety_mult == 1:
                # TODO
                pass
            
            if ninety_mult == 2:
                if group_tile.two_qubit_orientation == 0:
                    new_group_tile = self.SquareGridGroupTile(group_tile.ordered_qubits, tile_type=group_tile.tile_type,
                        pauli_map=group_tile._pauli_map, is_sharded=group_tile.is_sharded, two_qubit_orientation=1)
                    return new_group_tile
    
                else:
                    new_group_tile = self.SquareGridGroupTile(group_tile.ordered_qubits, tile_type=group_tile.tile_type,
                        pauli_map=group_tile._pauli_map, is_sharded=group_tile.is_sharded, two_qubit_orientation=0)
                    return new_group_tile

            else:
                raise NotImplementedError
        if ninety_mult == 1:
            if group_tile.tile_type == self.SquareGridGroupTile.POLYGON:
                new_points = []
                x_mat = []
                x_sum = 0
        
                y_mat = []
                y_sum = 0
        
                count = 0
                for p in group_tile.ordered_qubits:
                    count += 1
                    x_sum += p.x()
                    x_mat.append(p.x())
                
                    y_mat.append(p.y())
                    y_sum += p.y()
        
                P = np.array([x_mat, y_mat])
                poly_centr = self.calculate_centroid(group_tile.ordered_qubits)
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
                    
                # update colour mapping
                # so long as we're using this matrix, the color of group_tile.ordered_qubits[x] will be the same as new_points[x]
                # TODO make this prettier
                old_qubits = group_tile.ordered_qubits
                new_pauli_map = {}
                for ind in range(len(old_qubits)):
                    old_qu = old_qubits[ind]
                    pauli = group_tile.get_pauli_for_qubit(old_qu)
                    
                    new_qu = new_points[ind]
                    new_pauli_map[(new_qu.x(), new_qu.y())] = pauli
                    
                    

                new_group_tile = self.SquareGridGroupTile(new_points,None,self.SquareGridGroupTile.POLYGON, pauli_map=new_pauli_map, is_sharded=True)
                return new_group_tile


    def perimeter_size(self):
        # TODO turn this into boundRect overwrite
        return (self.grid_width * self.square_size, self.grid_height * self.square_size)
    
   
    
    def combine_paths_scene(self, scenepos_tile_list: [(QPointF, SquareGridGroupTile)]) -> (SquareGridGroupTile, QPointF):
        # # fix scene issues
        # for tup in scenepos_tile_list:
        #     scene_pos = tup[0]
        #     cur_path = tup[1]
        #     if cur_path.tile_type == self.SquareGridGroupTile.ARC:
        #         scene_pos.setX(scene_pos.x() + cur_path.start_point.x())
        #         scene_pos.setX(scene_pos.y() + cur_path.start_point.y())

        
        # SHOULD NOT CHANGE SCENE POS OF VERTICES BY THE END (WILL MESS UP ERROR GROUPS)
        # covert polys to scene then combine to megapoly then convert back to poly
        min_x = sys.maxsize
        min_y = sys.maxsize
        scene_vertices = []
        scene_vert_set = set()
        for tup in scenepos_tile_list:
            poly = tup[1].ordered_qubits
            # TODO handle case when it's not a polygon
            scene_pos = tup[0]
            for point in poly:
                pp = (scene_pos.x() + point.x(), scene_pos.y() + point.y())
                if pp not in scene_vert_set:
                    scene_vertices.append(QPointF(scene_pos.x() + point.x(), scene_pos.y() + point.y()))
                    scene_vert_set.add(pp)
        for vert in scene_vertices:
            if vert.x() < min_x:
                min_x = vert.x()
            if vert.y() < min_y:
                min_y = vert.y()
            
        centroid = self.calculate_centroid(scene_vertices)
        # vertices will also be bubble sorted
        new_vertices = self.bubble_sort_vertices(scene_vertices, centroid)
        print(f"new points: {new_vertices}")
        
        new_poly_vertices = []
        for point in new_vertices:
            new_poly_vertices.append(QPointF(point.x() - min_x, point.y() - min_y))
        print(f"new poly points: {new_poly_vertices}")


        # combine paulis
        scene_qubit_pos = dict()
        new_pauli_map = dict()
        for tup in scenepos_tile_list:
            scene_pos = tup[0]
            qubits = tup[1].ordered_qubits
            cur_path = tup[1]
            for qu in qubits:
                scene_qubit_pos[(qu.x() + scene_pos.x(), qu.y() + scene_pos.y())] = cur_path.get_pauli_for_qubit(qu) # might overwrite and that's fine
        for qu in scene_qubit_pos.keys():
            new_pauli_map[(qu[0] - min_x, qu[1] - min_y)] = scene_qubit_pos[qu] # convert from scene to local coordinates

        tile = self.convert_qubits_to_tile(new_poly_vertices, tile_type=self.SquareGridGroupTile.POLYGON, pauli_map=new_pauli_map, is_sharded=True)
        # subtract by min X and min Y?
        
        # QPainterpath, upperright most point

        return tile, QPointF(min_x, min_y)
        
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
        
        
    
    def find_closest_valid_scene_position(self, scenePos: QPointF):
        # square is easy,

        new_pos = self.find_closest_qubit_location(scenePos)
        return new_pos
        
        
        
    
    