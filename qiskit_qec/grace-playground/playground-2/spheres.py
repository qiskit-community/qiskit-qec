import vtk


def main():
    colors = vtk.vtkNamedColors()
    spheres = []
    # Setup a renderer, render window, and interactor
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d("MistyRose"))
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName("Animation")
    renderWindow.AddRenderer(renderer)

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    for x in range(5):
        for y in range(5):
            sphereSource = vtk.vtkSphereSource()
            sphereSource.SetCenter(0.0, x*10, y*10)
            sphereSource.SetRadius(2)
            sphereSource.SetPhiResolution(30)
            sphereSource.SetThetaResolution(30)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(sphereSource.GetOutputPort())
            actor = vtk.vtkActor()
            actor.GetProperty().SetColor(colors.GetColor3d("Peacock"))
            actor.GetProperty().SetSpecular(0.6)
            actor.GetProperty().SetSpecularPower(30)
            actor.SetMapper(mapper)

            # Add the actor to the scene
            renderer.AddActor(actor)





    # Render and interact
    renderWindow.Render()
    #renderer.GetActiveCamera().Zoom(0.8)
    renderWindow.Render()

    # Initialize must be called prior to creating timer events.
    renderWindowInteractor.Initialize()


    # start the interaction and timer
    renderWindow.Render()
    renderWindowInteractor.Start()


if __name__ == '__main__':
    main()