#!/usr/bin/env python

import vtk

def main():
    colors = vtk.vtkNamedColors()
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName("/Users/dsvandet/Software/Private/vtk/Testing/Data/cow.vtp")
    reader.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(reader.GetOutput())

    actor = vtk.vtkActor()
    actor.GetProperty().SetColor(colors.GetColor3d("Green"))
    actor.SetMapper(mapper)

    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetWindowName("CellPicking")

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    renderer.AddActor(actor)

    renderer.ResetCamera()

    renderer.SetBackground(colors.GetColor3d("Blue"))

    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__ == '__main__':
    main()

