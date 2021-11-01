#!/usr/bin/env python

"""
This is (almost) a direct C++ to Python transliteration of
 <VTK-root>/Examples/DataManipulation/Cxx/Cube.cxx from the VTK
 source distribution, which "shows how to manually create vtkPolyData"

A convenience function, mkVtkIdList(), has been added and one if/else
 so the example also works in version 6 or later.
If your VTK version is 5.x then remove the line: colors = vtk.vtkNamedColors()
 and replace the set background parameters with (1.0, 0.9688, 0.8594)

"""

import vtk


def mkVtkIdList(it):
    """
    Makes a vtkIdList from a Python iterable. I'm kinda surprised that
     this is necessary, since I assumed that this kind of thing would
     have been built into the wrapper and happen transparently, but it
     seems not.

    :param it: A python iterable.
    :return: A vtkIdList
    """
    vil = vtk.vtkIdList()
    for i in it:
        vil.InsertNextId(int(i))
    return vil


def main():
    colors = vtk.vtkNamedColors()

    # x = array of 8 3-tuples of float representing the vertices of a cube:
    x = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
         (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
         (3, 0, 1), (3, 1, 1), (3,1,0), (0.6, 0.6, 0.6), (0.5, 0.5, 0.5)]

    # pts = array of 6 4-tuples of vtkIdType (int) representing the faces
    #     of the cube in terms of the above vertices
    pts = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
           (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7),
           ]

    # lines = [(11,12, 8), (11,12,9), (11,12,10)]
    
    
    cube_center = (0.5, 0.5, 0.5)
    qubit_points = [(3, 0, 1), (3, 1, 1), (3,1,0)]
    
    
    

    # We'll create the building blocks of polydata including data attributes.
    cube = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    polys = vtk.vtkCellArray()
    scalars = vtk.vtkFloatArray()

    # Load the point, cell, and data attributes.
    for i, xi in enumerate(x):
        points.InsertPoint(i, xi)
    for pt in pts:
        polys.InsertNextCell(mkVtkIdList(pt))
        
    ptrav = polys.NewIterator()
    ptrav.GoToFirstCell()
    print("all my potential")
    while(not ptrav.IsDoneWithTraversal()):
        cell = ptrav.GetCurrentCell()
        print(f"cur cell is: {cell}")
        ptrav.GoToNextCell()
        
    for i, _ in enumerate(x):
        scalars.InsertTuple1(i, 5)
    

    # We now assign the pieces to the vtkPolyData.
    cube.SetPoints(points)
    cube.SetPolys(polys)
    cube.GetPointData().SetScalars(scalars)

    # Now we'll look at it.
    cubeMapper = vtk.vtkPolyDataMapper()
    cubeMapper.SetInputData(cube)
    cubeMapper.SetScalarRange(cube.GetScalarRange())
    cubeActor = vtk.vtkActor()
    cubeActor.SetMapper(cubeMapper)
   
   

    line_actors = []
    for qu in qubit_points:
        lineSource = vtk.vtkLineSource()
        lineSource.SetPoint1(cube_center)
        lineSource.SetPoint2(qu)
        # Setup actor and mapper
        lineMapper = vtk.vtkPolyDataMapper()
        lineMapper.SetInputConnection(lineSource.GetOutputPort())
        lineActor = vtk.vtkActor()
        lineActor.SetMapper(lineMapper)
        lineActor.GetProperty().SetColor(colors.GetColor3d('Red'))
        lineActor.GetProperty().SetLineWidth(4.55)
        line_actors.append(lineActor)
        print(f"line width is:{lineActor.GetProperty().GetLineWidth()} ")
        

    # The usual rendering stuff.
    camera = vtk.vtkCamera()
    camera.SetPosition(1, 1, 1)
    camera.SetFocalPoint(0, 0, 0)

    renderer = vtk.vtkRenderer()

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    
    renderer.AddActor(cubeActor)
    for la in line_actors:
        renderer.AddActor(la)
    renderer.SetActiveCamera(camera)
    renderer.ResetCamera()
    renderer.SetBackground(colors.GetColor3d("Cornsilk"))

    renWin.SetSize(600, 600)
    renWin.SetWindowName("Cube")

    # interact with data
    renWin.Render()
    iren.Start()


if __name__ == "__main__":
    main()