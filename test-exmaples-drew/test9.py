#!/usr/bin/env python

'''
This example shows how to create an unstructured grid.
'''

import vtk


def main():
    color_names = vtk.vtkNamedColors()

    x = [[0,0,0],[0,1,0],[0,2,0],[1,0,0],[1,1,0],[1,2,0],[2,0,0],[2,1,0],[2,2,0]]

    pts = [[0,3,4,1,2],[4,7,8,5,0],[1,4,5,2,1],[3,6,7,4,1]]

    color_array = vtk.vtkFloatArray()
    for item in pts:
        i = color_array.InsertNextTuple((item[4],))


    # Create two points, P0 and P1
    p0 = [0.0, 0.0, 0.0]
    p1 = [2.0, 2.0, 0.0]

    lineSource = vtk.vtkLineSource()
    lineSource.SetPoint1(p0)
    lineSource.SetPoint2(p1)

    line_mapper = vtk.vtkPolyDataMapper()
    line_mapper.SetInputConnection(lineSource.GetOutputPort())
    line_actor = vtk.vtkActor()
    line_actor.SetMapper(line_mapper)
    line_actor.GetProperty().SetLineWidth(4)
    line_actor.GetProperty().SetColor(color_names.GetColor3d("Peacock"))


    #x = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0], [0, 0, 1], [1, 0, 1], [2, 0, 1], [0, 1, 1],
    #    [1, 1, 1], [2, 1, 1], [0, 1, 2], [1, 1, 2], [2, 1, 2], [0, 1, 3], [1, 1, 3], [2, 1, 3], [0, 1, 4], [1, 1, 4],
    #     [2, 1, 4], [0, 1, 5], [1, 1, 5], [2, 1, 5], [0, 1, 6], [1, 1, 6], [2, 1, 6]]
    # Here we have kept consistency with the Cxx example of the same name.
    # This means we will use slicing in ugrid.InsertNextCell to ensure that the correct
    #  number of points are used.
    #pts = [[0, 1, 4, 3, 6, 7, 10, 9], [1, 2, 5, 4, 7, 8, 11, 10], [6, 10, 9, 12, 0, 0, 0, 0],
    #       [8, 11, 10, 14, 0, 0, 0, 0], [16, 17, 14, 13, 12, 15, 0, 0], [18, 15, 19, 16, 20, 17, 0, 0],
    #       [22, 23, 20, 19, 0, 0, 0, 0], [21, 22, 18, 0, 0, 0, 0, 0], [22, 19, 18, 0, 0, 0, 0, 0],
    #       [23, 26, 0, 0, 0, 0, 0, 0], [21, 24, 0, 0, 0, 0, 0, 0], [25, 0, 0, 0, 0, 0, 0, 0]]
    print(len(x), len(pts))

    renderer = vtk.vtkRenderer()

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    points = vtk.vtkPoints()
    for i in range(0, len(x)):
        points.InsertPoint(i, x[i])


    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.Allocate(100)

    ugrid.InsertNextCell(vtk.VTK_QUAD, 4, pts[0][:4])
    ugrid.InsertNextCell(vtk.VTK_QUAD, 4, pts[1][:4])
    ugrid.InsertNextCell(vtk.VTK_QUAD, 4, pts[2][:4])
    ugrid.InsertNextCell(vtk.VTK_QUAD, 4, pts[3][:4])

    '''
    ugrid.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, pts[0])
    ugrid.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, pts[1])
    ugrid.InsertNextCell(vtk.VTK_TETRA, 4, pts[2][:4])
    ugrid.InsertNextCell(vtk.VTK_TETRA, 4, pts[3][:4])
    ugrid.InsertNextCell(vtk.VTK_POLYGON, 6, pts[4][:6])
    ugrid.InsertNextCell(vtk.VTK_TRIANGLE_STRIP, 6, pts[5][:6])
    ugrid.InsertNextCell(vtk.VTK_QUAD, 4, pts[6][:4])
    ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, pts[7][:3])
    ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, pts[8][:3])
    ugrid.InsertNextCell(vtk.VTK_LINE, 2, pts[9][:2])
    ugrid.InsertNextCell(vtk.VTK_LINE, 2, pts[10][:2])
    ugrid.InsertNextCell(vtk.VTK_VERTEX, 1, pts[11][:1])
    '''

    ugrid.SetPoints(points)

    ugrid.GetCellData().SetScalars(color_array)

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(5)
    lut.SetTableRange(0,4)
    lut.Build()
    lut.SetTableValue(0, 0, 0, 0, 1) # Black
    lut.SetTableValue(1, 1, 0, 0, 1) # Red
    lut.SetTableValue(2, 0, 1, 0, 1) # Green
    lut.SetTableValue(3, 0, 0, 1, 1) # Blue
    lut.SetTableValue(4, 1, 1, 1, 1) # White 


    ugridMapper = vtk.vtkDataSetMapper()

    ugridMapper.SetScalarModeToUseCellData()
    ugridMapper.UseLookupTableScalarRangeOn()
    ugridMapper.SetLookupTable(lut)

    ugridMapper.SetInputData(ugrid)



    ugridActor = vtk.vtkActor()
    ugridActor.SetMapper(ugridMapper)
    ugridActor.GetProperty().SetColor(color_names.GetColor3d('Peacock'))
    ugridActor.GetProperty().EdgeVisibilityOn()

    renderer.AddActor(ugridActor)
    renderer.AddActor(line_actor)
    renderer.SetBackground(color_names.GetColor3d('Beige'))

    camera = vtk.vtkCamera()
    camera.SetPosition(1,1,50)
    camera.SetFocalPoint(1,1,0)
    #camera.SetViewUp(0,0,1)
    camera.ComputeViewPlaneNormal()

    renderer.SetActiveCamera(camera)

    renderer.ResetCamera()

    renWin.SetSize(640, 480)
    renWin.SetWindowName('UGrid')

    # Interact with the data.
    renWin.Render()

    iren.Start()


if __name__ == '__main__':
    main()