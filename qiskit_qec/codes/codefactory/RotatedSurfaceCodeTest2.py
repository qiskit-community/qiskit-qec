from math import e
#import vtkmodules.vtk as vtk
import numpy
from numpy.core.numeric import extend_all
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


def list_to_vtkPoints(points_list):
    points = vtk.vtkPoints()
    points.Allocate(len(points_list))
    for point in points_list:
        if len(point) == 2:
            points.InsertNextPoint(*point,0)
        else:
            points.InsertNextPoint(*point)
    return points

def list_to_vtkStringArray(string_list):
    string_array = vtk.vtkStringArray()
    string_array.Allocate(len(string_list))
    for string in string_list:
        string_array.InsertNextValue(string)
    return string_array

def vtkPoints_polgon_to_vtkPolyData(points):
    # Create the polygon
    num_points = points.GetNumberOfPoints()
    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(num_points)

    for k in range(num_points):
         polygon.GetPointIds().SetId(k, k)

    # Add the polygon to a list of polygons
    polygons = vtk.vtkCellArray()
    polygons.InsertNextCell(polygon)

    # Create a PolyData
    polygon_PolyData = vtk.vtkPolyData()
    polygon_PolyData.SetPoints(points)
    polygon_PolyData.SetPolys(polygons)

    return polygon_PolyData

def vtkIdList_to_set(vtk_id_list):
    id_set = set()
    for k in range(vtk_id_list.GetNumberOfIds()):
        id_set.add(vtk_id_list.GetId(k))

    return id_set

def vtkIdList_to_list(vtk_id_list):
    id_list = []
    for k in range(vtk_id_list.GetNumberOfIds()):
        id_list.append(vtk_id_list.GetId(k))

    return id_list

def make_int_entries(data):
    def make_int_entry(input):
        return tuple([round(item) for item in input])
    try:
        n = len(data[0])
        output = []
        for item in data:
            output.append(make_int_entry(item))
    except TypeError:
        output = make_int_entry(data)
    return output

def vtkPoints_to_int_list(vtk_points):
    # Return the points  with integer components instead of floating point components
    points_int = []
    for k in range(vtk_points.GetNumberOfPoints()):
        points_int.append(make_int_entries(vtk_points.GetPoint(k)))
    return points_int

def get_center_point(points):
    # Assumes that all points have the same z value
    # This should be done better later
    if len(points) !=4:
        raise("Error: Incorrect number of points. Need 4")
    else:
        x = (sum([item[0] for item in points])-2)/4 + 1/2.0
        y = (sum([item[1] for item in points])+2)/4 - 1/2.0
    return (x,y,points[0][2])
    
    
def id_to_points(object, cell_point_ids):
    return [object.GetPoint(cp_id) for cp_id in cell_point_ids]

def check_cell_arrays(grid):
    print("Checking cell arrays")
    numArrays = grid.GetCellData().GetNumberOfArrays()
    print(f"Number of Arrays is : {numArrays}")
    numCells = grid.GetNumberOfCells()
    print(f"Number of cells is {numCells}")
    for idx  in range(numArrays):

        array = grid.GetCellData().GetAbstractArray(idx)
        #print(array)
        numTuples = array.GetNumberOfTuples()
        name = array.GetName()
        if (name == None):
            name = ""
      
        if (numTuples < numCells):
            print("numTuples < numCells")
            print(f"Array {name} failed")
            print(f"Cell array {name} with indx ={idx} {array.GetNumberOfComponents()} components, only has {numTuples} tuples but there are {numCells} cells")
            for k in range(array.GetNumberOfTuples()):
                print(f"array: {array.GetTuple(k)}")
            return 1
      
        elif (numTuples > numCells):
            print("numTuples > numCells")
            print(f"Array {name} failed")
            print(f"Cell array {name} with {array.GetNumberOfComponents()} components, has {numTuples} tuples but there are only {numCells}  cells")

        
        else:
            print(f"Array {name} passed")


def print_cell_arrays(grid):
    print("Printing cell arrays")
    numArrays = grid.GetCellData().GetNumberOfArrays()
    print(f"Number of Arrays is : {numArrays}")
    numCells = grid.GetNumberOfCells()
    print(f"Number of cells is {numCells}")
    for idx  in range(numArrays):

        array = grid.GetCellData().GetAbstractArray(idx)
        numTuples = array.GetNumberOfTuples()
        name = array.GetName()
        if (name == None):
            name = ""
        print(f"Array: {name} at k={idx}")
        for k in range(array.GetNumberOfTuples()):
            print(f"array: {array.GetTuple(k)}")
        
      


def check_point_arrays(grid):
    print("Checking point arrays")
    numArrays = grid.GetPointData().GetNumberOfArrays()
    print(f"Number of Arrays is : {numArrays}")
    numPoints = grid.GetNumberOfPoints()
    print(f"Number of cells is {numPoints}")
    for idx  in range(numArrays):

        array = grid.GetPointData().GetAbstractArray(idx)
        #print(array)
        numTuples = array.GetNumberOfTuples()
        name = array.GetName()
        if (name == None):
            name = ""
      
        if (numTuples < numPoints):
            print("numTuples < numPoints")
            print(f"Array {name} failed")
            print(f"Point array {name} with indx ={idx} {array.GetNumberOfComponents()} components, only has {numTuples} tuples but there are {numPoints} points")
            for k in range(array.GetNumberOfTuples()):
                print(f"array: {array.GetTuple(k)}")
            return 1
      
        elif (numTuples > numPoints):
            print("numTuples > numPoints")
            print(f"Array {name} failed")
            print(f"Point array {name} with {array.GetNumberOfComponents()} components, has {numTuples} tuples but there are only {numPoints}  points")

        
        else:
            print(f"Array {name} passed")

class Set():
    def __init__(self):
        pass
class Group(Set):
    def __init__(self):
        super().__init__()
class PauliGroup(Group):
    def __init__(self):
        super().__init__()
class Code():
    def __init__(self):
        pass
class SubsystemCode(Code):
    def __init__(self):
        super().__init__()
class RotatedSurfaceCode(SubsystemCode):
    def __init__(self, model):
        pass
class GeometricCodeModelModeler:
    def __init__(self):
        pass

class RotatedSurfaceCodeGeometricModeler(GeometricCodeModelModeler):
    ANTICLOCKWISE  = 0
    CLOCKWISE      = 1

    X              = 0
    Z              = 1

    NOT_SET        = 0
    SET            = 1

    TOP            = 0
    BOTTOM         = 1

    KEEP           = 0
    DELETE         = 1

    NON_QUBIT      = -1
    INACTIVE_QUBIT = 0
    DATA_QUBIT     = 1
    ANCILLA_QUBIT  = 2
    FLAG_QUBIT     = 3

    # Boundary Transition markers

    TRANS_XZ       = 0
    TRANS_XZZX     = 1
    TRANS_ZXXZ     = 2
    TRANS_ZX       = 3

    STRAIGHT       = 0
    INSIDE_CORNER  = 1
    CORNER         = 2
    OUTSIDE_CORNER = 3

    def __init__(self) -> None:
        # Set Default colors
        self.colorX = 'Gray'
        self.colorZ = 'White'

        # Unstructured grid to store stabilizers

        self.grid = vtk.vtkUnstructuredGrid()

        # Note that these Arrays do not need to be class variables directly
        # Point Data

        self.points = vtk.vtkPoints()
        self.points.Allocate(1000)
        
        self.grid.SetPoints(self.points)

        # What type of point is it? One of 
        # 
        # NON_QUBIT      = -1
        # INACTIVE_QUBIT = 0
        # DATA_QUBIT     = 1
        # ANCILLA_QUBIT  = 2
        # FLAG_QUBIT     = 3

        self.point_type_name = "Point type"
        self.point_type = vtk.vtkTypeInt32Array()
        self.point_type.SetNumberOfComponents(1)
        self.point_type.SetName(self.point_type_name)
        #self.point_type.Allocate(500)

        self.grid.GetPointData().AddArray(self.point_type)

        # Cell Data: Is the face a stabilizer or not? One of 
        # X              = 0
        # Z              = 1

        self.cell_stabilizer_type_name = "Stabilizer type"
        self.cell_stabilizer_type = vtk.vtkTypeInt32Array()
        self.cell_stabilizer_type.SetNumberOfComponents(1)
        self.cell_stabilizer_type.SetName(self.cell_stabilizer_type_name)
        #self.cell_stabilizer_type.Allocate(500)

        self.grid.GetCellData().AddArray(self.cell_stabilizer_type)
 

        # Cell/Stabilizer weight (i.e. Number of qubits defingin stabilizer)
        self.cell_weight_name = "Cell weight"
        self.cell_weight = vtk.vtkTypeInt32Array()
        self.cell_weight.SetNumberOfComponents(1)
        self.cell_weight.SetName(self.cell_weight_name)
        #self.grid.GetCellData().Allocate(500)

        self.grid.GetCellData().AddArray(self.cell_weight)


        # Color of cell integers 0,1,2,3 ... using self.lut to translated to RGB color data
        self.cell_colors_name = "Cell colors"
        self.cell_colors = vtk.vtkFloatArray()
        self.cell_colors.SetNumberOfComponents(1)
        self.cell_colors.SetName(self.cell_colors_name)
        self.grid.GetCellData().AddArray(self.cell_colors)

        self.shape = None
        self.shape_flag = self.NOT_SET # No shape has been specified

        self.boundary = None
        self.boundary_flag = self.NOT_SET # No boundary has been specified
        
        self.boundary_types = None
        self.boundary_types_flag = self.NOT_SET # No boundary types has been specified

        self.correct_boundary_flag = False

        self.cellParity = 0

        self.code_model_flag = self.NOT_SET # Code model has not be generated

        # Create the LUT for scalars to colors
        self.lut = vtk.vtkLookupTable()
        self.lut.SetNumberOfTableValues(7)
        self.lut.SetTableRange(0,6)
        self.lut.Build()
        self.lut.SetTableValue(0, 0, 0, 0, 1) # Black
        self.lut.SetTableValue(1, 1, 0, 0, 1) # Red
        self.lut.SetTableValue(2, 0, 1, 0, 1) # Green
        self.lut.SetTableValue(3, 0, 0, 1, 1) # Blue
        self.lut.SetTableValue(4, 1, 1, 1, 1) # White 
        #self.lut.SetTableValue(5,0.788,0.6745,0.4196,1) # Dark Beige
        #self.lut.SetTableValue(6,248/255.0,233/255.0,174/255.0,1) # Light Beige
        self.lut.SetTableValue(5,201/255.0,172/255.0,107/255.0,1) # Dark Beige
        self.lut.SetTableValue(6,248/255.0,233/255.0,174/255.0,1) # Light Beige

        super().__init__()

    # get/set Shape
    def setShape(self, shape_data, override=False):
        if self.shape_flag == self.SET:
            if override == False:
                raise("Shape of code has already been processed. "
                    "Resetting the shape will reset the code to this new shape. "
                    "Use the override=True flag to force a reshape.")
            else:
                raise("Reshaping is not yet implemented")
        # Set the shape of the code. Convert to vtkPolyData if required
        if isinstance(shape_data, list):
            shape_points = list_to_vtkPoints(shape_data)
            self.shape = vtkPoints_polgon_to_vtkPolyData(shape_points)

        elif isinstance(shape_data, vtk.vtkPoints):
            self.shape = vtkPoints_polgon_to_vtkPolyData(shape_data)

        elif isinstance(shape_data, vtk.vtkPolyData):
            self.shape = shape_data

    def getShape(self):
        # At the moment self.shape is returned. Once the shape has been used
        if self.code_model_flag == self.NOT_SET:
            return self.shape
        else:
            return self.getGeneratedShape()

    def getGeneratedShape(self):
        # Generate the shape of the code from the code Model
        raise("Generating the shape has not yet been implemented")

    def getShapePoints(self):
        return self.getShape().GetPoints()

    def getShapeIntegerPoints(self):
        # Return the points defining the shape with integer components instead of floating point components
        shape_points_int = []
        shape_points = self.getShapePoints()
        for k in range(shape_points.GetNumberOfPoints()):
            shape_points_int.append(tuple([round(item) for item in shape_points.GetPoint(k)]))
        return shape_points_int

    # get/set Colors
    def setColors(self, colorX=1, colorZ=4):
        self.colorX = colorX
        self.colorZ = colorZ

    def getColors(self):
        return ['colorX', self.colorX, 'colorZ', self.colorZ]

    # get/set cellParity
    def setCellParity(self, cellParity):
        if cellParity not in [self.X, self.Z]:
            raise("Unknown cell parity: Use one of <class RotatedSurfaceCodeModel>.X (=0) or <class RotatedSurfaceCodeModel>.Z (=1)")
        self.cellParity = cellParity
    
    def getCellParity(self):
        return self.cellParity

    def encodeXZ(self,input):
        def encode(item):
            if item == 'X':
                return self.X
            if item == 'Z':
                return self.Z
        return [encode(item) for item in input]

    # get/set boundary
    def setBoundary(self, boundary, boundary_types):
        if isinstance(boundary, list):
            self.boundary = list_to_vtkPoints(boundary)
        elif isinstance(boundary, vtk.vtkPoints):
            self.boundary = boundary
        else:
            raise("Boundary must be a list or a vtkPoints class")

        self.boundary_flag = self.SET

        self.boundary_types = self.encodeXZ(self.expandBoundaryTypes(boundary_types))

        #self.boundary_types = list_to_vtkStringArray(self.expandBoundaryTypes(boundary_types))

    def expandBoundaryTypes(self, boundary_types):
        int_boundary = self.getBoundaryInteger() 
        
        lengths = []
        n = len(int_boundary)
        for k in range(n):
            start = int_boundary[k]
            finish = int_boundary[(k+1)%n]
            lengths.append(abs(start[0]-finish[0])+abs(start[1]-finish[1]))
        
        list_pairs = zip(lengths, boundary_types)
        boundary_array = []
        for section_length, type in list_pairs:
            if section_length == 0:
                section_length =1
            boundary_array = boundary_array + [type]*section_length

        print(boundary_array)

        return boundary_array

    def getBoundary(self):
        if self.boundary_flag == self.SET:
            return self.boundary
        else:
            raise("Boundary not yet set")

    def getBoundaryInteger(self):
        return vtkPoints_to_int_list(self.getGenerateBoundary())

    def getGenerateBoundary(self):
        if self.code_model_flag == self.SET:
            raise("Generate Boundary method not yet fully implemented")
        elif self.boundary_flag == self.SET:
            return self.boundary
        else:
            raise("Boundary not yet set")

    # Move and make utility more general
    def getBoundaryLength(self):
        point = [0,0]
        length = 0
        for next_point in self.getBoundary()[1:]:
            length = length + abs(point[0]-next_point[0]) + abs(point[1]-next_point[1]) + 1
            point = next_point
        return length

    # get/set boundary types
    def setBoundaryTypes(self, boundary_types):                  
        self.boundary_types = boundary_types
        self.boundary_types_flag = self.SET

    def getBoundaryTypes(self):
        return self.boundary_types

    # get/set boundary direction
    def setBoundaryDirection(self, direction):
        self.boundary_direction = direction

    def getBoundaryDirection(self):
        return self.boundary_direction

    def GetNextBoundaryPoint(self, point):
        pass

    # get/set cellUnit
    def setCellUnit(self, cellUnit):
        self.cellUnit = cellUnit

    def getCellunit(self):
        return self.cellUnit

    def generateShapeMask(self):
        # Put the origin into a PointArray
        points = vtk.vtkPoints()
        points.InsertNextPoint(0,0,0)

        # Make orgin point into a polyData object
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)

        # Set up the method to determine if a point is inside the shape
        # We only store the origin in the points so create a small mask array
        # and then use the backdoor IsInSideSurface routine. We should compare
        # this to putting all the points from the bounding box into the polyData
        # points ojbect. not sure which is faster.

        self.shape_mask  = vtk.vtkSelectEnclosedPoints()
        self.shape_mask.SetInputData(polyData)
        self.shape_mask.SetSurfaceData(self.shape)
        self.shape_mask.Update()

        # We can now use self.shape_mask.IsInsideSurface(point) to determine is 
        # a point is inside the shape

    def getShapeMask(self):
        return self.shape_mask

    def generateModel(self):
        self.generateShapeMask()
        self.generateBulkStabilizers()
        self.grid.BuildLinks()
        self.generateBoundaryStabilizers()

    def extractShapePointsArray(int_shape_boundary):
        import math

        shape_points=[]
        n = len(int_shape_boundary)
        for k in range(n):
            start = int_shape_boundary[k]
            next = int_shape_boundary[(k+1)%n]
            vector = (next[0]-start[0], next[1]-start[1], 0)
            size = vector[0]+vector[1]
            sign = math.copysign(1,size)
            if vector[0] == 0:
                for t in range(abs(size)):
                    shape_points.append((start[0], start[1]+sign*t,0.0))
            else:
                for t in range(abs(size)):
                    shape_points.append((start[0]+sign*t, start[1],0.0))

        return shape_points


    def generateBoundaryPointIds(self):
        '''
        This approach is simple but good enough for small numbers of points. If the 
        number of points on get very large then this approach should be replaced with
        one that restricts itself to the boundary points. For example one could first
        filter out the edge ids using the vtkFeatureEdges() method as follows:

        geometry = vtk.vtkGeometryFilter()
        geometry.SetInputData(self.grid)
        geometry.Update()

        id_filter = vtk.vtkIdFilter()
        id_filter.SetInputData(geometry.GetOutput())
        id_filter.SetPointIds(True)
        id_filter.SetCellIds(False)
        id_filter.SetPointIdsArrayName("Boundary Point Ids")
        id_filter.Update()

        edge_filter = vtk.vtkFeatureEdges()
        edge_filter.SetInputConnection(id_filter.GetOutputPort())
        edge_filter.BoundaryEdgesOn()
        edge_filter.FeatureEdgesOff()
        edge_filter.NonManifoldEdgesOff()
        edge_filter.ManifoldEdgesOff()
        edge_filter.ColoringOn()
        edge_filter.Update()

        self.boundary_point_ids = edge_filter.GetOutput().GetPointData().GetArray("Boundary Point Ids") # Lambda function

        num_of_shape_ids = self.boundary_point_ids.GetNumberOfValues()

        for k in range(num_of_shape_ids):
            print(f"k={k} with self.boundary_point_ids={self.boundary_point_ids.GetValue(k)}")

        These ids are not likely to be in the proper order. They would then need to be put into order.
        '''

        # Get a the shape boundary in integer coordinates
        int_shape_boundary = self.getShapeIntegerPoints()
        int_shape_boundary_points = make_int_entries(RotatedSurfaceCodeGeometricModeler.extractShapePointsArray(int_shape_boundary))

        # Invert the id to Point map from self.grid

        point_to_id = dict()
        for k in range(self.grid.GetNumberOfPoints()):
            point_to_id[make_int_entries(self.grid.GetPoint(k))] = k

        # Convert the boundary points to ids
        n = len(int_shape_boundary_points)
        boundary_point_ids = vtk.vtkIdList()
        boundary_point_ids.Allocate(n)
        for k in range(n):
            boundary_point_ids.InsertNextId(point_to_id[int_shape_boundary_points[k]])
        return  boundary_point_ids

    def generateBoundaryCellIds(self, boundary_point_ids, direction=None):
        boundary_cell_ids = vtk.vtkIdList()
        n = boundary_point_ids.GetNumberOfIds()
        start_cell_ids = vtk.vtkIdList()
        finish_cell_ids = vtk.vtkIdList()
        for k in range(n):
            start_id = boundary_point_ids.GetId(k)
            finish_id = boundary_point_ids.GetId((k+1)%n)
            self.grid.GetPointCells(start_id, start_cell_ids)
            self.grid.GetPointCells(finish_id, finish_cell_ids)
            # Find intersection (will be one id)
            # TODO: Please fix the very lazy code below. It should terminate the loops once
            # the cell as been found. This will do for now.
            for u in range(start_cell_ids.GetNumberOfIds()):
                u_id = start_cell_ids.GetId(u)
                for v in range(finish_cell_ids.GetNumberOfIds()):
                    v_id = finish_cell_ids.GetId(v)
                    if v_id == u_id:
                        edge_cell_id = u_id
        
            boundary_cell_ids.InsertNextId(edge_cell_id)
        return boundary_cell_ids
   
    def generateBulkStabilizers(self):
        if self.getShape() is None:
            raise " No shape has been set"

        boundingBox =  self.getShape().GetBounds() # (xmin,xmax, ymin,ymax, zmin,zmax)
        xmin = int(boundingBox[0])
        xmax = int(boundingBox[1])
        ymin = int(boundingBox[2])
        ymax = int(boundingBox[3])

        width = int(xmax - xmin + 1)
        height = int(ymax - ymin + 1)
        area = (width+1)*(height+1)
        #The following is temp: to be removed
        self.area = area

        #self.grid.GetPointData().GetArray(self.point_type_name).Allocate(area)
        self.grid.GetCellData().GetArray(self.cell_stabilizer_type_name).Allocate(area)
        self.grid.GetCellData().GetArray(self.cell_weight_name).Allocate(area)
        self.grid.GetCellData().GetArray(self.cell_colors_name).Allocate(area)

        #self.cell_stabilizer_type.Allocate(area)
        #self.cell_colors.Allocate(area)
        #self.cell_weight.Allocate(area)

        if not self.grid.AllocateEstimate(area, 4):
            raise("Insufficient memory to allocate for vtkUnstructuredGrid")

        for x in range(xmin, xmax, self.cellUnit):
            for y in range(ymax, ymin, -self.cellUnit):
                point  = (x,y,0)
                center = (x+0.5, y-0.5, 0)
                if self.getShapeMask().IsInsideSurface(center):
                    # Add a cell for the associate weight four stabilizer
                    stab_points = [point, (x,y-1,0), (x+1,y-1,0), (x+1,y,0)]
                    #ids = [self.grid.GetPointData().GetArray(self.point_type_name).InsertNextPoint(stab_point) for stab_point in stab_points]

                    ids = []
                    for stab_point in stab_points:
                        stab_id = self.grid.GetPoints().InsertNextPoint(stab_point)
                        ids.append(stab_id)
                        self.grid.GetPointData().GetArray(self.point_type_name).InsertValue(stab_id, self.DATA_QUBIT)

                    ids = [self.grid.GetPoints().InsertNextPoint(stab_point) for stab_point in stab_points]

                    cell_id = self.grid.InsertNextCell(vtk.VTK_QUAD, 4, ids)
                    # Record stabilizer type and associated color
                    if (x+y+self.cellParity)%2 == 0:
                        # Stabilizer is of type X
                        self.grid.GetCellData().GetArray(self.cell_stabilizer_type_name).InsertValue(cell_id, self.X)
                        #print(type(self.grid.GetCellData().GetArray(self.cell_colors_name)))
                        self.grid.GetCellData().GetArray(self.cell_colors_name).InsertValue(cell_id, self.colorX)
                        #self.cell_stabilizer_type.InsertValue(cell_id, self.X)
                        #self.cell_colors.InsertValue(cell_id, self.colorX)
                    else:
                        # Stabilizer is of type Z
                        self.grid.GetCellData().GetArray(self.cell_stabilizer_type_name).InsertValue(cell_id, self.Z)
                        self.grid.GetCellData().GetArray(self.cell_colors_name).InsertValue(cell_id, self.colorZ)
                        #self.cell_stabilizer_type.InsertValue(cell_id, self.Z)
                        #self.cell_colors.InsertValue(cell_id, self.colorZ)
                    #self.grid.GetCellData().GetArray(self.cell_weight_name).InsertValue(cell_id, 4)

                    self.grid.GetCellData().GetArray(self.cell_weight_name).InsertValue(cell_id, 4)
        
        # This above process does not check for existing points and so duplicates them. We now clean
        # up these duplicates using the CleanPolyData filter
        
        geometry = vtk.vtkGeometryFilter()
        geometry.SetInputData(self.grid)
        geometry.Update()
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(geometry.GetOutputPort())
        clean.PointMergingOn()
        clean.Update()

        filter = vtk.vtkAppendFilter()
        filter.AddInputData(clean.GetOutput())
        filter.Update()

        temp_grid = vtk.vtkUnstructuredGrid()
        temp_grid.ShallowCopy(filter.GetOutput())

        self.grid = temp_grid



    def generateBoundaryOptions(self, boundary_point_ids, boundary_cell_ids):

        boundary_options = []

        n = boundary_point_ids.GetNumberOfIds()
        for k in range(n):
            #previous_point_id = boundary_point_ids.GetId((k-1)%n)
            current_point_id = boundary_point_ids.GetId(k)
            next_point_id = boundary_point_ids.GetId((k+1)%n)
            back_cell_id = boundary_cell_ids.GetId((k-1)%n)
            front_cell_id = boundary_cell_ids.GetId(k)

            if back_cell_id == front_cell_id: # Outside Corner (going around anticlockwise)
                type = self.cell_stabilizer_type.GetValue(front_cell_id)
                boundary_options.append((current_point_id, current_point_id, int(not type), self.CORNER))
                boundary_options.append((current_point_id, next_point_id, type, self.OUTSIDE_CORNER))
   
            elif self.cell_stabilizer_type.GetValue(back_cell_id) == self.cell_stabilizer_type.GetValue(front_cell_id): # Inside corner
                type = self.cell_stabilizer_type.GetValue(front_cell_id)
                boundary_options.append((current_point_id, next_point_id, type, self.INSIDE_CORNER))

            else: # straight run
                type = self.cell_stabilizer_type.GetValue(front_cell_id)
                boundary_options.append((current_point_id, next_point_id, type, self.STRAIGHT))

        return boundary_options

    def generateBoundaryStabilizers(self):
        if self.shape is None:
            raise(" No shape has been set")
        if self.boundary is None:
            raise("No Boundary is set")

        #check_point_arrays(self.grid)
        #check_cell_arrays(self.grid)

        print("\n---------Generate Boundary Options, Ids, Points and Types---------\n")
        
        boundary_point_ids = self.generateBoundaryPointIds() 

        boundary_cell_ids = self.generateBoundaryCellIds(boundary_point_ids)

        boundary_options = self.generateBoundaryOptions(boundary_point_ids, boundary_cell_ids)

        print(f"Boundary Options (length = {len(boundary_options)}) are\n {boundary_options}")
        print(f"Boundary Types (length {len(self.boundary_types)}) are\n {self.boundary_types}")

        # First identify which qubits (if any) need to be removed,

        #check_point_arrays(self.grid)
        #check_cell_arrays(self.grid)

        print("\n---------First identify which qubits (if any) need to be removed---------\n")

        delete_qubits_indices = set()
        delete_qubits_points_to_indices  = dict()
        n = len(boundary_options)
        for k in range(n):
            current_section = boundary_options[k]
            if current_section[3] == self.CORNER: # Checking corners points for changes 
                print(f"Checking corner with k={k}, current_section = {current_section}")
                if current_section[2] != self.boundary_types[k]: # A change in corner type is requested
                    print("...A Change corner is requested")
                    delete_qubits_indices.add(k)
                    delete_qubits_points_to_indices[boundary_options[k][0]]=k

        print("Delete qubits indices")
        print(delete_qubits_indices)
        print("To indices")
        print(delete_qubits_points_to_indices)

        # Determine if any of the qubits belong to the same cell and group together

        #check_point_arrays(self.grid)
        #check_cell_arrays(self.grid)

        print("\n---------Determine if any of the qubits belong to the same cell and group together---------\n")

        partitioned_deleted_qubits_indices = []
        partitioned_deleted_qubits_cell_ids = []
        
        while len(delete_qubits_indices)>0:
            point_id_index = next(iter(delete_qubits_indices))
            print(f"Working with point {point_id_index}")
            point_id = boundary_options[point_id_index][0]
            print(f"Working with point id {point_id}")

            # Find cell_id for the cell that point_id is attached to
            cell_ids = vtk.vtkIdList()
            self.grid.GetPointCells(point_id, cell_ids)
            cell_id = cell_ids.GetId(0)
            print(f"Cell id is {cell_id}")

            # Get the point ids for this cell
            cell_point_ids = vtk.vtkIdList()
            self.grid.GetCellPoints(cell_id, cell_point_ids)

            print(cell_point_ids.GetNumberOfIds())
            cell_point_ids_list = vtkIdList_to_list(cell_point_ids)

            print(cell_point_ids_list)

            cell_point_indices_set = set()
            # Convert points to indices of boundary_options
            for item in cell_point_ids_list:
                try:
                    index = delete_qubits_points_to_indices[item]
                    cell_point_indices_set.add(index)
                except KeyError:
                    continue
            
            print("cell point indices set")
            print(cell_point_indices_set)

            # Grab the other cell points that are marked to be changed
            cell_point_indices_set = cell_point_indices_set.intersection(delete_qubits_indices)

            delete_qubits_indices = delete_qubits_indices - cell_point_indices_set

            partitioned_deleted_qubits_indices.append(list(cell_point_indices_set))
            partitioned_deleted_qubits_cell_ids.append([cell_id,cell_point_ids_list])

        deleted_qubits_ids = []
        for group in partitioned_deleted_qubits_indices:
            for item in group:
                deleted_qubits_ids.append(boundary_options[item][0])

        deleted_cell_ids = []
        for group in partitioned_deleted_qubits_cell_ids:
            deleted_cell_ids.append(group[0])

        print("Deleted qubits ids")
        print(deleted_qubits_ids)       
        print("Deleted cell ids")
        print(deleted_cell_ids) 
        print("Partitioned deleted qubits indices")
        print(partitioned_deleted_qubits_indices)
        print("Partitioned deleted qubits cell ids")
        print(partitioned_deleted_qubits_cell_ids)

        # Before deleteing qubits add in new cells relative to defined boundary

        #check_point_arrays(self.grid)
        #check_cell_arrays(self.grid)


        print("---------Checking number of tuples for cell_colot_name = NumberOfCells---------")
        print(f"Number of cells = {self.grid.GetNumberOfCells()}")
        print(f"Number of color tuples is = {self.grid.GetCellData().GetArray(self.cell_colors_name).GetNumberOfTuples()}")
        if self.grid.GetNumberOfCells() != self.grid.GetCellData().GetArray(self.cell_colors_name).GetNumberOfTuples():
            print("Error: Number of cells and color tuples do not agree")

        # First add in new corner cells
        print("\n---------Added in new corner cells---------\n")

        for k in range(len(partitioned_deleted_qubits_cell_ids)):
            cell_group = partitioned_deleted_qubits_cell_ids[k]
            cell_id = cell_group[0]
            cell_point_ids = cell_group[1]
            cell_point_ind_del = partitioned_deleted_qubits_indices[k]
            cell_point_ids_del = [boundary_options[k][0] for k in cell_point_ind_del]
            new_cell_point_ids = [pid for pid in cell_point_ids if pid not in cell_point_ids_del]

            #print(f"cell_group = {cell_group}")
            #print(f"cell_id = {cell_id}")
            #print(f"cell_point_ids = {cell_point_ids}")
            #print(f"cell_point_ind_del = {cell_point_ind_del}")
            #print(f"cell_point_ids_del = {cell_point_ids_del}")
            #print(f"new_cell_point_ids = {new_cell_point_ids}\n")
            length = len(new_cell_point_ids)
            check_cell_arrays(self.grid)
            if length == 2:
                cp = get_center_point(id_to_points(self.grid.GetPoints(), cell_point_ids))
                cp_id = self.grid.GetPoints().InsertNextPoint(cp)
                self.grid.GetPointData().GetArray(self.point_type_name).InsertValue(cp_id, self.NON_QUBIT)

                cid = self.grid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [new_cell_point_ids[0], new_cell_point_ids[1], cp_id])

                cell_type = self.grid.GetCellData().GetArray(self.cell_stabilizer_type_name).GetValue(cell_id)
                self.grid.GetCellData().GetArray(self.cell_stabilizer_type_name).InsertValue(cid, cell_type)

                cell_color = self.grid.GetCellData().GetArray(self.cell_colors_name).GetValue(cell_id)
                self.grid.GetCellData().GetArray(self.cell_colors_name).InsertValue(cid, cell_color)

                self.grid.GetCellData().GetArray(self.cell_weight_name).InsertValue(cid, 2)

            elif length == 3:

                cid = self.grid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [new_cell_point_ids[0], new_cell_point_ids[1], new_cell_point_ids[2]])

                cell_type = self.grid.GetCellData().GetArray(self.cell_stabilizer_type_name).GetValue(cell_id)
                self.grid.GetCellData().GetArray(self.cell_stabilizer_type_name).InsertValue(cid, cell_type)

                cell_color = self.grid.GetCellData().GetArray(self.cell_colors_name).GetValue(cell_id)
                self.grid.GetCellData().GetArray(self.cell_colors_name).InsertValue(cid, cell_color)

                self.grid.GetCellData().GetArray(self.cell_weight_name).InsertValue(cid, 3)

        #check_point_arrays(self.grid)
        #check_cell_arrays(self.grid)

        # Update boundary_options and self.boundary_types for deleted cells
        
        print("\n-------------Update boundary_options and self.boundary_types for deleted cells------------\n")
        print(boundary_options)
        print(self.boundary_types)
        n = len(boundary_options)
        boundary_options_marker = [0]*n
        print(f"Number of boundary options is {n}")
        print(f"Number of boundary_types is {len(self.boundary_types)}")
        
        
        # Something is wrong here!! Deleting all options. Should only delete the
        # corners that changed. Not all corners
        # Working with detailed boundary being given
        for k in range(n):
            previous_opt = boundary_options[(k-1)%n]
            current_opt = boundary_options[k]
            next_opt = boundary_options[(k+1)%n]

            if current_opt[3] == self.CORNER and current_opt[2] != self.boundary_types[k]:
                boundary_options_marker[(k-1)%n] = self.DELETE
                boundary_options_marker[k] = self.DELETE
                boundary_options_marker[(k+1)%n] = self.DELETE
        
        new_boundary_options = []
        new_boundary_types = []
        print(len(boundary_options))
        print(len(self.boundary_types))
        for k in range(n):
            if boundary_options_marker[k] == self.KEEP:
                new_boundary_options.append(boundary_options[k])
                new_boundary_types.append(self.boundary_types[k])

        print("Boundary markers: boundary_options_marker")
        print(boundary_options_marker)
        
        boundary_options = new_boundary_options
        self.boundary_types = new_boundary_types

        print("New boundaries: boundary_options, self.boundary_types")
        print(boundary_options)
        print(self.boundary_types)



        # Add boundary stabilizers as defined by new boundary conditions

        n = len(boundary_options)
        for k in range(n):
            print(f"k={k} out of {n}")
            back_section = boundary_options[k]
            front_section = boundary_options[(k+1)%n]
            back_type = self.boundary_types[k]
            front_type = self.boundary_types[(k+1)%n]

            previous_point_id = back_section[0]
            current_point_id = front_section[0]
            next_point_id = front_section[1]
  
            if front_section[3] == self.INSIDE_CORNER:
                print("Inside Corner")
                continue # only temp 
            else: 
                # Check for illegal boundary changes XY<->YX
                if (back_type == front_section[2] and front_type == back_section[2]):
                    print(f"Illegal boundary change at {previous_point_id} -> {current_point_id} -> {next_point_id}")
                    print(f"{current_point_id} occurs at coordinates {self.grid.GetPoint(current_point_id)}")
                    if self.correct_boundary_flag == True:
                        print("Correcting boundary condition")
                        pass
                    else:
                        raise("Illegal boundary condition")
                # Need to add mark for stabilizer etc


                points = self.grid.GetPoints()
                point_type = self.grid.GetPointData().GetArray(self.point_type_name)

                cell_colors = self.grid.GetCellData().GetArray(self.cell_colors_name)
                cell_stabilizer_type = self.grid.GetCellData().GetArray(self.cell_stabilizer_type_name)
                cell_weights = self.grid.GetCellData().GetArray(self.cell_weight_name)
                

                if boundary_options[(k+1)%n][2] != front_type:
                    print(f"{boundary_options[(k+1)%n][2]} compared to {front_type}")
                    cp = self.grid.GetPoint(current_point_id)
                    np = self.grid.GetPoint(next_point_id)
                    dx = int(np[0]-cp[0])
                    dy = int(np[1]-cp[1])
                    if dx == 0:
                        if dy > 0: # Delta = (0,1) -> (u+0.5,v+0.5)
                            peek = (cp[0]+0.5, cp[1]+0.5,0)
                        else:  # Delta = (0,-1) -> (u-0.5,v-0.5)
                            peek = (cp[0]-0.5, cp[1]-0.5,0)
                    else:
                        if dx > 0: # Delta = (1,0) -> (u+0.5, v-0.5)
                            peek = (cp[0]+0.5, cp[1]-0.5,0)
                        else: #Delta = (-1,0) -> (u-0.5, v + 0.5)
                            peek = (cp[0]-0.5, cp[1]+0.5,0)
                        
                    print(f"Adding point {peek}")
                    print(f"Number of points is {points.GetNumberOfPoints()}")
                    pt_id = points.InsertNextPoint(peek)
                    print(f"Number of points is {points.GetNumberOfPoints()}")
                    print(f"with id = {pt_id}")
                    cid = self.grid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [current_point_id, next_point_id, pt_id])
                    print(f"Added triange {[current_point_id, next_point_id, pt_id]}")
                    if front_type == self.X:
                        color = self.colorX
                    else:
                        color = self.colorZ
                    cell_colors.InsertValue(cid, color)
                    print(f"color = {color}")
                    cell_weights.InsertValue(cid, 3)
                    cell_stabilizer_type.InsertValue(cid, color)
                    point_type.InsertValue(pt_id, self.NON_QUBIT)

        print(f"Number of points is {points.GetNumberOfPoints()}")
        print(f"{cell_colors.GetNumberOfValues()}")
        print(f"{cell_weights.GetNumberOfValues()}")
        print(f"{cell_stabilizer_type.GetNumberOfValues()}")
        

        #check_point_arrays(self.grid)
        #check_cell_arrays(self.grid)


        # Or Just mark cells and points as keep or delete then
        # run filter to extract data only associated with stuff to
        # Keep. Then run a clean Filter to remove all other points (if possible)

        # Delete marked cells and marked points
        print("\n--------------Delete marked cells and marked points ---------------\n")

        #deleted_qubits_ids
        #deleted_cell_ids

        # first delete cells
        print("Delete cells")

        geometry = vtk.vtkGeometryFilter()
        geometry.SetInputData(self.grid)
        geometry.Update()

        polydata = geometry.GetOutput()
        polydata.EditableOn()
        polydata.BuildCells()
        polydata.BuildLinks()

        for cell_id in deleted_cell_ids:
            polydata.DeleteCell(cell_id)
            
        for pt_id in deleted_qubits_ids:
            polydata.DeletePoint(pt_id)

        polydata.RemoveDeletedCells()
        polydata.RemoveDeletedPoints()

        filter = vtk.vtkAppendFilter()
        filter.AddInputData(polydata)
        filter.Update()
        
      

        #


        '''
        idd = vtk.vtkIdList()
        self.grid.GetPointCells(178, idd)
        sam = []
        for k in range(idd.GetNumberOfIds()):
            sam.append(idd.GetId(k))


        print(self.grid.GetNumberOfPoints())

        geometry = vtk.vtkGeometryFilter()
        geometry.SetInputData(self.grid)
        geometry.BuildCells()
        geometry.Update()

        
        polydata = geometry.GetOutput()

        print(polydata.GetNumberOfPoints())
     
        polydata.EditableOn()
        polydata.BuildLinks()

        new_color = vtk.vtkFloatArray()
        temp = []
        new_color.Allocate(self.area)

        for k in range(self.cell_colors.GetNumberOfValues()):
            if k not in [138]:
                temp.append(self.cell_colors.GetValue(k))

        new_color = numpy_to_vtk(numpy.array(temp))

        self.cell_colors = new_color



        my_ids = vtk.vtkIdList()
        polydata.GetCellPoints(157, my_ids)

        

        sample = []
        for k in range(my_ids.GetNumberOfIds()):
            sample.append(my_ids.GetId(k))

        print(sample)

         # Now remove the necessary qubits and reshape the resulting stabilizers
        polydata.DeleteCell(138)
        polydata.RemoveDeletedCells()
        print(polydata.GetNumberOfPoints())

        filter = vtk.vtkAppendFilter()
        filter.AddInputData(polydata)
        filter.Update()


        new_color = vtk.vtkFloatArray()
        temp = []
        new_color.Allocate(self.area)

        for k in range(self.cell_colors.GetNumberOfValues()):
            if k not in [138]:
                temp.append(self.cell_colors.GetValue(k))

        new_color = numpy_to_vtk(numpy.array(temp))

        self.cell_colors = new_color

        temp_grid = vtk.vtkUnstructuredGrid()
        temp_grid.ShallowCopy(filter.GetOutput())

        self.grid = temp_grid 

        temp_grid.GetCellData().SetScalars(self.cell_colors)

        
        return
        '''

        n = len(partitioned_deleted_qubits_indices)
        for k in range(n):
            cell_id = partitioned_deleted_qubits_cell_ids[k][0]
            qubit_indices = partitioned_deleted_qubits_indices[k]
            cell_type = self.grid.GetCellType(cell_id)
            cell_point_ids = partitioned_deleted_qubits_cell_ids[k][1]
            del_point_id = boundary_options[qubit_indices[0]][0]
            print(f"Working with cell_id={cell_id} and indices {qubit_indices} of cell type {cell_type}")
            cell_weight  = self.cell_weight.GetValue(cell_id)
            print(f"Cell weight is {cell_weight} and cell ids are {cell_point_ids}")
            print(f"Del point id is {del_point_id}")
            
            
            #geometry = vtk.vtkGeometryFilter()
            #geometry.SetInputData(self.grid)
            #geometry.Update()
            #poly_data = geometry.GetOutput().GetCellData()
            #poly_data.DeleteCell(cell_id)
            #poly_data.RemoveDeletedCells()


        

        

        # Secondly, add any inside corner weight three stabilizers



        # Thirdly, change remaining boarders (checking for allowable changes)




#self.grid.GetCellData().AddArray()


def main():
    colors = vtk.vtkNamedColors()

    # Create a shape and boundary for the rotated surface code:
    # Short cuts should be available for standard shapes and boundary: e.g. square, rentangular
    # Need a method to discribe the boundary (probably XXYXXYX and 
    # [P1,P2,P3,...p_k] where p_i are the coordinates of the transitions)

    # Will will start with a square d x d rotated surface code

    # Code Representation (Symplectic Matrix, ...)

    # TranslationMap (Matrix to Qubit index)
    # TranslationMap (Matrix to Local Coordinate)

    code = 3

    # Create shape of rotated surface code
    # Example 1: Square Code 9 qubit with uniform boundary
    if code==1:
        shape = [[0,0], [0,-2], [2,-2], [2,0]]

        full_boundary = \
            [[0, 0], 'Z', 
            [0, 0], 'X',
            [0, -2], 'X',
            [0, -2], 'Z',
            [2,-2], 'Z',
            [2, -2], 'X',
            [2, 0], 'X',
            [2, 0],'Z'
        ]

        #full_boundary = \
        #    [[0,0], 'Z']
    
    # Example 2: Random Polygon shape code

    if code == 2:
        # Shape of code
        shape = [[0,0],[0,-3],[14,-3],[14,8],[16,8],[16,10],[1,10],[1,8],[7,8],[7,3],[2,3],[2,0]]


        # Exmaple Boundary with illigal changes that will result in non-abelian generators

        # Boundary: there needs to be a variety of ways to indicate the boundary conditions
        # Boundary of code. Boundary change points listed in anticlockwise direction starting with origin [0,0]

        full_boundary = \
            [[0, 0], 'Z', 
            [0, 0],  'X', 
            [0, -3], 'Z', 
            [0,-3], 'X',
            [11, -3], 'Z',
            [14,-3], 'Z',
            [14,-3], 'Z',
            [14, -2], 'X',
            [14, 2], 'Z', # Illegal changing boundaries [14,1]-[14,2] Z -> X combined with [14,2]-[14,3] X -> Z
            [14, 4], 'X',
            [15, 8], 'Z',
            [16, 8], 'Z',
            [16, 8], 'Z',
            [16, 10], 'X',
            [16,10], 'X',
            [11, 10], 'Z',
            [6, 10], 'X',
            [1,10], 'X',
            [1,10], 'X',
            [1, 9], 'Z',
            [1,8], 'Z',
            [1,8], 'Z',
            [7, 6], 'X',
            [7, 3], 'Z',
            [4, 3], 'X',
            [2,3], 'X',
            [2,3], 'X',
            [2, 0], 'Z']

        full_boundary = \
            [[0, 0], 'Z', 
            [0, 0],  'X', 
            [0, -3], 'Z', 
            [0,-3], 'X',
            [11, -3], 'Z',
            [14,-3], 'Z',
            [14,-3], 'Z',
            [14, -2], 'X',
            [14, 2], 'X', # Illegal changing boundaries [14,1]-[14,2] Z -> X combined with [14,2]-[14,3] X -> Z
            [14, 4], 'X',
            [15, 8], 'Z',
            [16, 8], 'Z',
            [16, 8], 'Z',
            [16, 10], 'X',
            [16,10], 'X',
            [11, 10], 'Z', # Illegal changing boundaries
            [6, 10], 'X',
            [1,10], 'X',
            [1,10], 'X',
            [1, 9], 'Z',
            [1,8], 'Z',
            [1,8], 'Z',
            [7, 6], 'X',
            [7, 3], 'Z',
            [4, 3], 'X',
            [2,3], 'X',
            [2,3], 'X',
            [2, 0], 'Z']

        full_boundary = \
            [[0, 0], 'Z', 
            [0, 0],  'X', 
            [0, -3], 'Z', 
            [0,-3], 'X',
            [11, -3], 'Z',
            [14,-3], 'Z',
            [14,-3], 'Z',
            [14, 0], 'X',
            [14, 2], 'X', 
            [14, 5], 'Z',
            [15, 8], 'Z',
            [16, 8], 'Z',
            [16, 8], 'Z',
            [16, 10], 'X',
            [16,10], 'X',
            [11, 10], 'X', 
            [6, 10], 'X',
            [1,10], 'X',
            [1,10], 'Z',
            [1,8], 'Z',
            [1,8], 'Z',
            [7,8], 'Z',
            [7, 5], 'X', 
            [7, 4], 'Z',
            [4, 3], 'Z',
            [2, 3], 'Z',
            [2, 3], 'Z',
            [2, 0], 'Z']

    if code == 3:
        shape = \
            [[0,0],
            [0,-1],
            [2,-1],
            [2,1],
            [1,1],
            [1,0]]

        full_boundary = \
             [[0,0], 'X',
            [0,0], 'X',
            [0,-1], 'X',
            [0,-1], 'X',
            [2,-1], 'X',
            [2,-1], 'X',
            [2,1], 'X',
            [2,1], 'X',
            [1,1], 'X',
            [1,1], 'X',
            [1,0], 'X']           

    boundary = full_boundary[::2]

    # Boundary types for code. List in anticlockwise order. It is assume that the last coordinate 
    # links back to the origin point. In the case below the last region aligns with the first
    # region.

    boundary_types = full_boundary[1::2]





    model = RotatedSurfaceCodeGeometricModeler()

    model.setShape(shape)
    model.setBoundary(boundary, boundary_types)
    model.setBoundaryDirection(model.ANTICLOCKWISE)

    model.setColors(colorX=5, colorZ=6)
    model.setCellParity(model.X) 
    model.setCellUnit(1)
    model.generateModel()

    label_mapper = vtk.vtkLabeledDataMapper()
    label_mapper.SetInputData(model.grid)
    label_mapper.GetLabelTextProperty().SetColor(colors.GetColor3d("Black"))
    label_mapper.GetLabelTextProperty().SetFontSize(24)
    label_actor = vtk.vtkActor2D()
    label_actor.SetMapper(label_mapper)

    
    cell_centers = vtk.vtkCellCenters()
    cell_centers.SetInputData(model.grid)


    label_mapper2 = vtk.vtkLabeledDataMapper()
    label_mapper2.SetInputConnection(cell_centers.GetOutputPort())
    label_mapper2.GetLabelTextProperty().SetColor(colors.GetColor3d("Orange"))
    label_mapper2.GetLabelTextProperty().SetFontSize(24)
    label_actor2 = vtk.vtkActor2D()
    label_actor2.SetMapper(label_mapper2)
    

    #code = SubsystemCode()
    #code.generate(model)

    #return
  
    model.grid.GetCellData().SetScalars(model.grid.GetCellData().GetArray(model.cell_colors_name))
    model.grid.GetCellData().SetActiveScalars(model.cell_colors_name)

    mapper = vtk.vtkDataSetMapper()
    mapper.SetScalarModeToUseCellData()
    mapper.UseLookupTableScalarRangeOn()
    
    mapper.SetLookupTable(model.lut)
    mapper.SetInputData(model.grid)

    #print(f"{model.merge.GetNumberOfPoints()}")

    #for i in range(model.merge.GetNumberOfPoints()):
    #    print(f"i={i}, pt={model.merge.GetPoint(i)}")


    #writer = vtk.vtkXMLUnstructuredGridWriter()
    #writer.SetFileName("RSC-3x3.vtu")
    #writer.SetInputData(model.grid)
    #writer.Write()

    
    #geometry = vtk.vtkGeometryFilter()
    #geometry.SetInputData(model.grid)
    #geometry.Update()

    
    '''
    edge_filter = vtk.vtkFeatureEdges()
    edge_filter.SetInputData(geometry.GetOutput())
    edge_filter.BoundaryEdgesOn()
    edge_filter.FeatureEdgesOn()
    edge_filter.NonManifoldEdgesOn()
    edge_filter.ManifoldEdgesOn()
    edge_filter.ColoringOn()
    edge_filter.Update()

    emapper = vtk.vtkPolyDataMapper()
    emapper.SetInputConnection(edge_filter.GetOutputPort())

    eactor = vtk.vtkActor()
    eactor.SetMapper(emapper)

    eactor.GetProperty().SetEdgeVisibility(1)
    eactor.GetProperty().SetEdgeColor(0.9,0.9,0.4)
    eactor.GetProperty().SetLineWidth(10)
    eactor.GetProperty().SetPointSize(12)
    eactor.GetProperty().SetRenderLinesAsTubes(1)
    eactor.GetProperty().SetRenderPointsAsSpheres(1)
    eactor.GetProperty().SetVertexVisibility(1)
    eactor.GetProperty().SetVertexColor(0.5,1.0,0.8)
    #eactor.GetProperty().EdgeVisibilityOn()

    '''
    

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    #actor.GetProperty().SetEdgeVisibility(0)
    #actor.GetProperty().SetEdgeColor(0,0,0)shape ossdf
    #actor.GetProperty().SetLineWidth(0)
    #actor.GetProperty().SetPointSize(12)
    #actor.GetProperty().SetRenderLinesAsTubes(0)
    #actor.GetProperty().SetRenderPointsAsSpheres(1)
    #actor.GetProperty().SetVertexVisibility(0)
    #actor.GetProperty().SetVertexColor(0,0,0)
    #actor.GetProperty().EdgeVisibilityOn()

    renderer = vtk.vtkRenderer()
    renderer.UseFXAAOn()

    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName('RotatedSurfaceCodeExample')
    renderWindow.AddRenderer(renderer)

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(style)

    cam_orient_manipulator = vtk.vtkCameraOrientationWidget()
    cam_orient_manipulator.SetParentRenderer(renderer)
    # Enable the widget.
    cam_orient_manipulator.On()

    renderer.AddActor(actor)
    renderer.AddViewProp(label_actor)
    renderer.AddViewProp(label_actor2)
    renderer.SetBackground(colors.GetColor3d('White'))

    renderWindow.SetSize(640,640)
    renderWindow.Render()
    renderWindowInteractor.Start()


if __name__ == "__main__":
    main()





    



