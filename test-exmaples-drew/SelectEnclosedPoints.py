import vtk

def polygonPointsToVtkPolyData(points_list):
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

def main():
    colors = vtk.vtkNamedColors()

    shape_points = [[0,0],[0,-3],[14,-3],[14,8],[16,8],[16,10],[1,10],[1,8],[7,8],[7,3],[2,3],[2,0]]
    shape_points = [[0,0], [2,0], [2,-2], [0,-2]]
    shape_points = [[0,0],[0,-3],[14,-3],[14,8.0],[16,8.0],[16,10.0],[0,10]]
    shape_points = [[0,0],[0,-3],[14,-3],[14,5],[16,5],[16,8],[0,8]]
    shape = polygonPointsToVtkPolyData(shape_points)

    points = vtk.vtkPoints()
    points.InsertNextPoint(-1.0,0.0,0)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)

    selectEnclosed  = vtk.vtkSelectEnclosedPoints()
    selectEnclosed.SetInputData(polyData)
    selectEnclosed.SetSurfaceData(shape)
    selectEnclosed.Update()

    print(f'point included: {selectEnclosed.IsInside(0)}')

    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(shape)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Visualize
    renderer = vtk.vtkRenderer()
    renderer.UseFXAAOn() # Enable anti-aliasing 

    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('DimGray'))

    renderWindow.SetWindowName('SelectEnclosePoints')
    renderWindow.Render()
    renderWindowInteractor.Start()
    

if __name__ == "__main__":
    main()
