#!/usr/bin/env python

"""
=========================================================================

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

"""
#!/usr/bin/env python

import vtk


def main():
    colors = vtk.vtkNamedColors()

    coneSource = vtk.vtkConeSource()
    # coneSource.SetResolution(60)
    # coneSource.SetCenter(-2,0,0)

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(coneSource.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetDiffuseColor(colors.GetColor3d('bisque'))

    # Visualize
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('Salmon'))

    renderWindow.SetSize(640, 480)
    renderWindow.SetWindowName('Cone')
    renderWindow.Render()
    renderWindowInteractor.Start()


if __name__ == '__main__':
    main()