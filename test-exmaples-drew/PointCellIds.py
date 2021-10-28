import vtk

if vtk.VTK_VERSION_NUMBER >= 89000000000:
    VTK890=True

def main():
    colors = vtk.vtkNamedColors()

    # Coordinates for a square stabilizer
    point_list = [(0,0),(1,0),(1,1),(0,1)]
    labels = [0,1,2,3,1] # Last elements a salar representing color of cell

    # Convert coordinates in 3D vtkPoints
    points = vtk.vtkPoints()
    for point in point_list:
        points.InsertNextPoint(*point,0)

    
    # Create an UnstructuredGrid
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.Allocate(100)

    ugrid.InsertNextCell(vtk.VTK_QUAD,4,labels[:4])

    ugrid.SetPoints(points)


    sphereSource = vtk.vtkSphereSource()
    sphereSource.Update()

    print(f"There are {sphereSource.GetOutput().GetNumberOfPoints()} points")
    print(f"There are {sphereSource.GetOutput().GetNumberOfCells()} cells")

    idFilter = vtk.vtkIdFilter()
    idFilter.SetInputConnection(sphereSource.GetOutputPort())
    if VTK890:
        idFilter.SetPointIdsArrayName("ids")
        idFilter.SetCellIdsArrayName("ids")
    else:
        idFilter.SetIdsArrayName("ids")
    idFilter.Update()

    print(f"point arrays:")

    for i in range(idFilter.GetOutput().GetPointData().GetNumberOfArrays()):
        print(f"{idFilter.GetOutput().GetPointData().GetArrayName(i)}")

    print("cell arrays: ")
  
    for i in range(idFilter.GetOutput().GetCellData().GetNumberOfArrays()):
        print(f"{idFilter.GetOutput().GetCellData().GetArrayName(i)}")

    pointIds = vtk.vtkIdTypeArray()
    pointIds = pointIds.SafeDownCast(idFilter.GetOutput().GetPointData().GetArray("ids"))

    print(f"There are {pointIds.GetNumberOfTuples()} point ids")

    cellIds = vtk.vtkIdTypeArray()
    cellIds = cellIds.SafeDownCast(idFilter.GetOutput().GetCellData().GetArray("ids"))

    print(f"There are {cellIds.GetNumberOfTuples()} cell ids")

    color_array = vtk.vtkFloatArray()
    i = color_array.InsertNextTuple((labels[4],))

    ugrid.GetCellData().SetScalars(color_array)

    # Create the color LUT

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(5)
    lut.SetTableRange(0,4)
    lut.Build()
    lut.SetTableValue(0, 0, 0, 0, 1) # Black
    lut.SetTableValue(1, 1, 0, 0, 1) # Red
    lut.SetTableValue(2, 0, 1, 0, 1) # Green
    lut.SetTableValue(3, 0, 0, 1, 1) # Blue
    lut.SetTableValue(4, 1, 1, 1, 1) # White 

    # Create a mapper to map the color

    ugridMapper = vtk.vtkDataSetMapper()
    ugridMapper.SetScalarModeToUseCellData()
    ugridMapper.UseLookupTableScalarRangeOn()
    ugridMapper.SetLookupTable(lut)
    ugridMapper.SetInputData(ugrid)

    # Get the edges

    # Extract the edges of the triangles just found.
    extractEdges = vtk.vtkExtractEdges()
    extractEdges.SetInputData(ugrid)   

    # Create the tube filter

    tubeFilter = vtk.vtkTubeFilter()
    tubeFilter.SetInputConnection(extractEdges.GetOutputPort())
    tubeFilter.SetRadius(0.025)
    tubeFilter.SetNumberOfSides(50)
    tubeFilter.Update()

    # Create the edge mapper

    edgeMapper = vtk.vtkPolyDataMapper()
    edgeMapper.SetInputConnection(tubeFilter.GetOutputPort())
    edgeMapper.SetScalarRange(0, 26)

    # Create edge actor

    edgeActor = vtk.vtkActor()
    edgeActor.SetMapper(edgeMapper)
    edgeActor.GetProperty().SetSpecular(0.6)
    edgeActor.GetProperty().SetSpecularPower(30)

    # Create a renderer, render window, and interactor
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(edgeActor)
    renderer.SetBackground(colors.GetColor3d('SlateGray'))

    aCamera = vtk.vtkCamera()
    aCamera.Azimuth(-40.0)
    aCamera.Elevation(50.0)

    renderer.SetActiveCamera(aCamera)
    renderer.ResetCamera()

    renderWindow.SetSize(640, 640)
    renderWindow.SetWindowName('Stabilizer')
    renderWindow.Render()

    renderWindowInteractor.Start()


if __name__ == "__main__":
    main()