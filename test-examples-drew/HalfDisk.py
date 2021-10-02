
import vtk
import math

# Create a simple 1/2 disk cell

def main():
    colors = vtk.vtkNamedColors()

    # Setup four points
    points = vtk.vtkPoints()

    N = 100
    R=2
    angle = math.pi/N

    for k in range(N+1):
        nangle = k*angle
        points.InsertNextPoint(R*math.cos(nangle),R*math.sin(nangle),0)

    # Create the polygon
    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(N+1)

    for k in range(N+1):
         polygon.GetPointIds().SetId(k, k)

    # Add the polygon to a list of polygons
    polygons = vtk.vtkCellArray()
    polygons.InsertNextCell(polygon)

    # Create a PolyData
    polygonPolyData = vtk.vtkPolyData()
    polygonPolyData.SetPoints(points)
    polygonPolyData.SetPolys(polygons)

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polygonPolyData)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d('Blue'))

    # Visualize
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName('Polygon')
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(style)

    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('DimGray'))
    renderWindow.Render()
    renderWindowInteractor.Start()


if __name__ == '__main__':
    main()
