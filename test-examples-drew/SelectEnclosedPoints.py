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
    #shape_points = [[0,0], [2,0], [2,-2], [0,-2]]
    shape = polygonPointsToVtkPolyData(shape_points)

    points = vtk.vtkPoints()
    points.InsertNextPoint(-1.0,0.0,0)
    points.InsertNextPoint(0,0,0)
    points.InsertNextPoint(13,7,0)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)

    selectEnclosed  = vtk.vtkSelectEnclosedPoints()
    selectEnclosed.SetInputData(polyData)
    selectEnclosed.SetSurfaceData(shape)
    selectEnclosed.Update()

    #for i in range(3):
    #    print(f'point included: {selectEnclosed.IsInside(i)}')

    test = [0,0,0]

    print(f'point included: {selectEnclosed.IsInsideSurface(test)}')

    # in order to render a irregular polygon it needs to be decomposed into triangles
    filter = vtk.vtkTriangleFilter()
    filter.SetInputData(shape)
    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(filter.GetOutputPort())

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
