#!/usr/bin/env python

import vtk


def main():
    colors = vtk.vtkNamedColors()

    ren1 = vtk.vtkRenderer()

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren1)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    quadric = vtk.vtkQuadric()
    quadric.SetCoefficients(0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3)

    sample = vtk.vtkSampleFunction()
    sample.SetSampleDimensions(50, 50, 50)
    sample.SetImplicitFunction(quadric)
    sample.ComputeNormalsOff()

    trans = vtk.vtkTransform()
    trans.Scale(0.5, 0.5, 0.333)

    sphere = vtk.vtkSphere()
    sphere.SetRadius(0.25)
    sphere.SetTransform(trans)

    trans2 = vtk.vtkTransform()
    trans2.Scale(0.25, 0.5, 1.0)

    sphere2 = vtk.vtkSphere()
    sphere2.SetRadius(0.25)
    sphere2.SetTransform(trans2)

    booleanUnion = vtk.vtkImplicitBoolean()
    booleanUnion.AddFunction(sphere)
    booleanUnion.AddFunction(sphere2)
    booleanUnion.SetOperationType(0)  # boolean Union

    extract = vtk.vtkExtractGeometry()
    extract.SetInputConnection(sample.GetOutputPort())
    extract.SetImplicitFunction(booleanUnion)

    shrink = vtk.vtkShrinkFilter()
    shrink.SetInputConnection(extract.GetOutputPort())
    shrink.SetShrinkFactor(0.5)

    dataMapper = vtk.vtkDataSetMapper()
    dataMapper.SetInputConnection(shrink.GetOutputPort())
    dataActor = vtk.vtkActor()
    dataActor.SetMapper(dataMapper)

    # outline
    outline = vtk.vtkOutlineFilter()
    outline.SetInputConnection(sample.GetOutputPort())

    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInputConnection(outline.GetOutputPort())

    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper(outlineMapper)
    outlineActor.GetProperty().SetColor(0, 0, 0)

    # Add the actors to the renderer, set the background and size
    #
    ren1.AddActor(outlineActor)
    ren1.AddActor(dataActor)
    ren1.SetBackground(colors.GetColor3d("SlateGray"))

    renWin.SetSize(640, 480)
    renWin.SetWindowName('ExtractData')

    renWin.Render()
    ren1.GetActiveCamera().Azimuth(30)
    ren1.GetActiveCamera().Elevation(30)

    renWin.Render()
    iren.Start()


if __name__ == '__main__':
    main()