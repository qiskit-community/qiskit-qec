import vtk

def polygon_points_to_vtkPolyData(points_list):
    points = vtk.vtkPoints()
    for point in points_list:
        if len(point) == 2:
            points.InsertNextPoint(*point,0)
        else:
            points.InsertNextPoint(*point)
    
    # Create the polygon
    num_points = len(points_list)
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

def isInsidePolygon(point, polygon):
    pass

class Code():
    def __init__(self):
        pass

class SubsystemCode(Code):
    def __init__(self):
        super().__init__()

class RotatedSurfaceCode(SubsystemCode):
    def __init__(self, model, cell_type, x_color, z_color, unit=1):
        pass

class CodeModel:
    def __init__(self):
        pass

class RotatedSurfaceCodeModel(CodeModel):
    def __init__(self) -> None:

        # Set Default colors
        self.colorX = 'Gray'
        self.colorZ = 'White'

        self.points = vtk.vtkPoints()
        self.grid = vtk.vtkUnstructuredGrid()
        self.grid.SetPoints(self.points)

        self.shape = None
        self.boundaries = None

        self.cellType = None

        super().__init__()

    # get/set Shape
    def setShape(self, shape='rectangular', height=2, width=2):
        self.shape = shape
        self.shape_height = height
        self.shape_width = width

    def getShape(self):
        return self.shape

    # get/set Colors
    def setColors(self, colorX='Gray', colorZ='White'):
        self.colorX = colorX
        self.colorZ = colorZ

    def getColors(self):
        return ['colorX', self.colorX, 'colorZ', self.colorZ]

    # get/set cellType
    def setCellType(self, cellType):
        self.cellType = cellType
    
    def getCellType(self):
        return self.cellType

    # get/set boundaries
    def setBoundaries(self, boundaries):
        self.boundaries = boundaries

    def getBoundaries(self):
        return self.boundaries

    # get/set boundary types
    def setBoundaryTypes(self, boundaryTypes):
        self.boundaryTypes = boundaryTypes

    def getBoundaryTypes(self):
        return self.boundaryTypes

    #get/set cellUnit
    def setCellUnit(self, cellUnit):
        self.cellUnit = cellUnit

    def getCellunit(self):
        return self.cellUnit

    def createStabilizers(self):
        self.createWeightFourStabilizer()
        #self.createWeightTwoStabilizer()

    def createWeightFourStabilizer(self):
        if self.shape is None:
            raise " No shape has been set"

        # Put the origin into a PointArray
        points = vtk.vtkPoints()
        points.InsertNextPoint(0,0,0)

        # Make points into a polyData object
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)

        # Set up the method to determine if a point is inside the shape
        # We only store the origin in the points so create a small mask array
        # and then use the backdoor IsInSideSurface routine. We should compare
        # this to putting all the points from the bounding box into the polyData
        # points ojbect. not sure which is faster.

        shape_surface  = vtk.vtkSelectEnclosedPoints()
        shape_surface.SetInputData(polyData)
        shape_surface.SetSurfaceData(self.shape)
        shape_surface.Update()

        # We can now use shape_surface.IsInsideSurface(point) to determine is 
        # a point is inside the shape
        
        boundingBox =  self.shape.GetBounds() # (xmin,xmax, ymin,ymax, zmin,zmax)
        xmin = int(boundingBox[0])
        xmax = int(boundingBox[1])
        ymin = int(boundingBox[2])
        ymax = int(boundingBox[3])

        width = int(xmax - xmin + 1)
        height = int(ymax - ymin + 1)

        if not self.grid.AllocateEstimate((width+1)*(height+1), 4):
            raise "Insufficient memory to allocate for vtkUnstructuredGrid"
        square = [[0,0,0], [1,0,0], [1,-1,0], [0,-1,0]]

        for x in range(xmin, xmax, self.cellUnit):
            for y in range(ymax, ymin, -self.cellUnit):
                point  = [x,y,0]
                if shape_surface.IsInsideSurface(point):
                    stab_points = [point, [x,y-1,0], [x+1,y-1,0], [x+1,y,0]]
                    ids = [self.points.InsertNextPoint(stab_point) for stab_point in stab_points]
                    self.grid.InsertNextCell(vtk.VTK_QUAD, 4, ids)
        

    def createWeightTwoStabilizer(self):
        if self.shape is None:
            raise " No shape has been set"
        if self.boundaries is None:
            raise "No Boundaries are set"

def main():
    colors = vtk.vtkNamedColors()

    # Create a shape and boundary for the rotated surface code:
    # Short cuts should be available for standard shapes and boundaries: e.g. square, rentangular
    # Need a method to discribe the boundaries (probably XXYXXYX and [P1,P2,P3,...p_k] where p_i are the coordinates of the transitions)

    # Will will start with a square d x d rotated surface code

    # Code Representation (Symplectic Matrix, ...)

    # TranslationMap (Matrix to Qubit index)
    # TranslationMap (Matrix to Local Coordinate)

    # Create shape of rotated surface code
    # Squae Code
    #shape_points = [[0,0], [2,0], [2,-2], [0,-2]]
    # Random Polygon shape code
    shape_points = [[0,0],[0,-3],[14,-3],[14,8],[16,8],[16,10],[1,10],[1,8],[7,8],[7,3],[2,3],[2,0]]
    # Convert list of boundary cornders to a vtkPolyData object
    shape = polygon_points_to_vtkPolyData(shape_points)

    boundaries = [[0,0],[2,0], [2,0],[2,-2], [2,-2],[0,-2], [0,-2],[0,0]]
    boundary_types = ['X', 'Z', 'X', 'Z']
    origin_operator = 'X'

   

    model = RotatedSurfaceCodeModel()
    model.setShape(shape)
    model.setBoundaries(boundaries)
    model.setBoundaryTypes(boundary_types)
    model.setColors(colorX='Gray', colorZ='Green')
    model.setCellType('X')
    model.setCellUnit(1)
    model.createStabilizers()

    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(model.grid)
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
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





    



