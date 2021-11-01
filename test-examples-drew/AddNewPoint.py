import vtk

def main():
    colors = vtk.vtkNamedColors()  

    points = vtk.vtkPoints()
    points.Allocate(100)
    id1 = points.InsertNextPoint((0,0,0))
    id2 = points.InsertNextPoint((0,1,0))
    id3 = points.InsertNextPoint((1,0,0))
    id4 = points.InsertNextPoint((1,1,0))

    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)

    grid.InsertNextCell(vtk.VTK_QUAD, 4, [id1,id2,id4,id3])

    geometry = vtk.vtkGeometryFilter()
    geometry.SetInputData(grid)
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

    grid = temp_grid   

    idn = grid.GetPoints().InsertNextPoint((2,0,0))
    grid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [id4,id3,idn])


    # now Delete

    



    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(grid)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    renderer = vtk.vtkRenderer()
    
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(style)

    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('Black'))

    renderWindow.Render()
    renderWindowInteractor.Start()


if __name__=="__main__":
    main()