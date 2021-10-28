import vtk

def main():
    points = vtk.vtkPoints()

    points.InsertNextPoint(0, 0, 0)
    points.InsertNextPoint(1, 0, 0)
    points.InsertNextPoint(1, 1, 0)
    points.InsertNextPoint(0, 1, 0)


    square  = vtk.vtkQuad()
    square.GetPointIds().SetId(0, 0)
    square.GetPointIds().SetId(1, 1)
    square.GetPointIds().SetId(2, 2)
    square.GetPointIds().SetId(3, 3)

    line0 = vtk.vtkLine()
    line0.GetPointIds().SetId(0, 0)
    line0.GetPointIds().SetId(1, 1)

    line1 = vtk.vtkLine()
    line1.GetPointIds().SetId(0, 1)
    line1.GetPointIds().SetId(1, 2)

    line2 = vtk.vtkLine()
    line2.GetPointIds().SetId(0, 2)
    line2.GetPointIds().SetId(1, 3)

    #lines = vtk.vtkCellArray()
    #lines.InsertNextCell(line0)
    #lines.InsertNextCell(line1)
    #lines.InsertNextCell(line2)

    #squares = vtk.vtkCellArray()
    #squares.InsertNextCell(square)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    #polydata.SetPolys(squares)
    #polydata.SetLines(lines)

    polydata.EditableOn()
    polydata.Allocate(100)
    id = polydata.InsertNextCell(vtk.VTK_QUAD, 4, [0,1,2,3])

    ids = vtk.vtkIdList()
    polydata.GetCellPoints(id, ids)
    ll = []
    for k in range(ids.GetNumberOfIds()):
        ll.append(ids.GetId(k))
    print(ll)

    # Tell the polydata to build 'upward' links from points to cells.
    polydata.BuildLinks()
    # Mark a cell as deleted.
    polydata.DeleteCell(0)
    # Remove the marked cell.
    polydata.RemoveDeletedCells()

    polydata.InsertNextCell(vtk.VTK_TRIANGLE, 3, [0,1,2])


    # Visualize
    colors = vtk.vtkNamedColors()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputDataObject(polydata)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d("Peacock"))
    actor.GetProperty().SetLineWidth(4)

    renderer= vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetWindowName("DeleteCells")

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d("Silver"))

    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__=="__main__":
    main()