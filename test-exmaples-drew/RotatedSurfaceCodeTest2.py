from math import e
import vtk

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
class GeometricCodeModelCreator:
    def __init__(self):
        pass

class RotatedSurfaceCodeGeometricModelCreator(GeometricCodeModelCreator):
    ANTICLOCKWISE = 0
    CLOCKWISE     = 1

    X             = 0
    Z             = 1

    NOT_SET       = 0
    SET           = 1

    TOP           = 0
    BOTTOM        = 1

    def __init__(self) -> None:
        # Set Default colors
        self.colorX = 'Gray'
        self.colorZ = 'White'

        # Unstructured grid to store stabilizers
        self.points = vtk.vtkPoints()
        self.grid = vtk.vtkUnstructuredGrid()
        self.grid.SetPoints(self.points)

        # Field Data

        self.cell_stabilizer_type = vtk.vtkTypeInt32Array()
        self.cell_stabilizer_type.SetNumberOfComponents(1)
        self.cell_stabilizer_type.SetName("StabilizerType")

        self.cell_colors = vtk.vtkFloatArray()

        self.shape = None
        self.shape_flag = self.NOT_SET

        self.boundary = None
        self.boundary_flag = self.NOT_SET

        self.cellParity = 0

        self.codeModel_flag = self.NOT_SET

        # Create the LUT for scalars to colors
        self.lut = vtk.vtkLookupTable()
        self.lut.SetNumberOfTableValues(5)
        self.lut.SetTableRange(0,4)
        self.lut.Build()
        self.lut.SetTableValue(0, 0, 0, 0, 1) # Black
        self.lut.SetTableValue(1, 1, 0, 0, 1) # Red
        self.lut.SetTableValue(2, 0, 1, 0, 1) # Green
        self.lut.SetTableValue(3, 0, 0, 1, 1) # Blue
        self.lut.SetTableValue(4, 1, 1, 1, 1) # White 

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
        if self.codeModel_flag == self.NOT_SET:
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

    # get/set boundary
    def setBoundary(self, boundary):
        if isinstance(boundary, list):
            self.boundary = list_to_vtkPoints(boundary)
        elif isinstance(boundary, vtk.vtkPoints):
            self.boundary = boundary
        else:
            raise("Boundary must be a list or a vtkPoints class")

    def getBoundary(self):
        return self.boundary

    # get/set boundary
    def setBoundary_list(self, boundary_list):
        self.boundary_list = boundary_list

    def getBoundary_list(self):
        return self.boundary_list

    def getBoundaryLength(self):
        point = [0,0]
        length = 0
        for next_point in self.getBoundary()[1:]:
            length = length + abs(point[0]-next_point[0]) + abs(point[1]-next_point[1]) + 1
            point = next_point
        return length

    # get/set boundary types
    def setBoundaryTypes(self, boundaryTypes):                  
        self.boundaryTypes = boundaryTypes

    def getBoundaryTypes(self):
        return self.boundaryTypes

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
        int_shape_boundary_points = make_int_entries(RotatedSurfaceCodeGeometricModelCreator.extractShapePointsArray(int_shape_boundary))

        # Invert the id to Point map from self.grid

        point_to_id = dict()
        for k in range(self.grid.GetNumberOfPoints()):
            point_to_id[make_int_entries(self.grid.GetPoint(k))] = k

        # Convert the boundary points to ids
        boundary_ids = []
        for k in range(len(int_shape_boundary_points)):
            boundary_ids.append(point_to_id[int_shape_boundary_points[k]])

        return  boundary_ids

    def generateBoundaryCellIds(self, boundary_point_ids, direction):
        cell_ids = []
        for k in range(len(boundary_point_ids)):
            start = boundary_point_ids
        
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

        self.cell_stabilizer_type.Allocate(area)
        self.cell_colors.Allocate(area)

        if not self.grid.AllocateEstimate(area, 4):
            raise("Insufficient memory to allocate for vtkUnstructuredGrid")

        for x in range(xmin, xmax, self.cellUnit):
            for y in range(ymax, ymin, -self.cellUnit):
                point  = (x,y,0)
                if self.getShapeMask().IsInsideSurface(point):
                    # Add a cell for the associate weight four stabilizer
                    stab_points = [point, (x,y-1,0), (x+1,y-1,0), (x+1,y,0)]
                    ids = [self.points.InsertNextPoint(stab_point) for stab_point in stab_points]
                    cell_id = self.grid.InsertNextCell(vtk.VTK_QUAD, 4, ids)
                    # Record stabilizer type and associated color
                    if (x+y+self.cellParity)%2 == 0:
                        # Stabilizer is of type X
                        self.cell_stabilizer_type.InsertValue(cell_id, 0)
                        self.cell_colors.InsertValue(cell_id, self.colorX)
                    else:
                        # Stabilizer is of type Z
                        self.cell_stabilizer_type.InsertValue(cell_id, 1)
                        self.cell_colors.InsertValue(cell_id, self.colorZ)
        
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

    def generateBoundaryStabilizers(self):
        if self.shape is None:
            raise(" No shape has been set")
        if self.boundary is None:
            raise("No Boundary is set")
        
        boundary_point_ids = self.generateBoundaryPointIds() 
        top_side_types = self.generateBoundaryCellIds(boundary_point_ids, self.TOP)
        bottom_side_types = self.generateBoundaryCellIds(boundary_point_ids, self.BOTTOM)





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

    # Create shape of rotated surface code
    # Example 1: Square Code 9 qubit with uniform boundary
    # shape = [[0,0], [2,0], [2,-2], [0,-2]]
    # boundary = [[0,0], [2,0], [2,-2], [0,-2]]
    # boundary_types = ['X', 'Z', 'X', 'Z']

    # Example 2: Random Polygon shape code
    # Shape of code
    shape = [[0,0],[0,-3],[14,-3],[14,8],[16,8],[16,10],[1,10],[1,8],[7,8],[7,3],[2,3],[2,0]]
    # Boundary of code. Boundary change points listed in anticlockwise direction starting with origin [0,0]
    boundary = \
        [[0, 0], [0, -3], [11, -3], 
        [14, -2], [14, 2], [14, 4], 
        [15, 8], [16, 1], [6, 15], 
        [6, 10], [1, 9], [7, 6], 
        [7, 3], [4, 3], [2, 0]]
    # Boundary types for code. List in anticlockwise order. Order must be the same as what the boundary is 
    boundary_types = \
        ['I', 'X', 'I', 
        'Z', 'I', 'X', 
        'Z', 'X', 'I', 
        'X', 'Z', 'X', 
        'Z', 'X']

    model = RotatedSurfaceCodeGeometricModelCreator()

    model.setShape(shape)

    model.setBoundary(boundary)
    model.setBoundaryTypes(boundary_types)
    model.setBoundaryDirection(model.ANTICLOCKWISE)

    model.setColors(colorX=2, colorZ=3)
    model.setCellParity(model.X) 
    model.setCellUnit(1)
   
    model.generateModel()


    #code = SubsystemCode()
    #code.generate(model)

    model.grid.GetCellData().SetScalars(model.cell_colors)

    mapper = vtk.vtkDataSetMapper()
    mapper.SetScalarModeToUseCellData()
    mapper.UseLookupTableScalarRangeOn()
    mapper.SetLookupTable(model.lut)
    mapper.SetInputData(model.grid)


    #print(f"{model.merge.GetNumberOfPoints()}")

    #for i in range(model.merge.GetNumberOfPoints()):
    #    print(f"i={i}, pt={model.merge.GetPoint(i)}")

    
    geometry = vtk.vtkGeometryFilter()
    geometry.SetInputData(model.grid)
    geometry.Update()

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
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetEdgeVisibility(1)
    actor.GetProperty().SetEdgeColor(0.9,0.9,0.4)
    actor.GetProperty().SetLineWidth(10)
    actor.GetProperty().SetPointSize(12)
    actor.GetProperty().SetRenderLinesAsTubes(1)
    actor.GetProperty().SetRenderPointsAsSpheres(1)
    actor.GetProperty().SetVertexVisibility(1)
    actor.GetProperty().SetVertexColor(0.5,1.0,0.8)
    actor.GetProperty().EdgeVisibilityOn()

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
    renderer.SetBackground(colors.GetColor3d('Black'))

    renderWindow.Render()
    renderWindowInteractor.Start()


    #code = RotatedSurfaceCode()
    #code.setModel(model)
    #code.Modified()


if __name__ == "__main__":
    main()





    



