import vtk


# Call back function


def sphereCallback(obj, event):
    print(f'obj is \n{obj}\n event is {event}')


def main():
    colors = vtk.vtkNamedColors()

    # colors.SetColor('bkg', [0.1, 0.2, 0.4, 1.0])

    # A renderer and render window
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d('MidnightBlue'))

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(renderer)
    renwin.SetWindowName("SphereWidget")

    # An interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renwin)

    # A Sphere widget
    sphereWidget = vtk.vtkSphereWidget()
    sphereWidget.SetInteractor(interactor)
    sphereWidget.SetRepresentationToSurface()
    sphereWidget.GetSphereProperty().SetColor(colors.GetColor3d("BurlyWood"))
    #sphereWidget.OnLeftButtonDown(sphereCallback)

    # Connect the event to a function
    sphereWidget.AddObserver("MouseButtonPress", sphereCallback)
    
    
    sphereWidget.AddObserver("TapEvent", sphereCallback)
    sphereWidget.AddObserver("MiddleButtonPressEvent", sphereCallback)
    sphereWidget.AddObserver("RightButtonPressEvent", sphereCallback)
    sphereWidget.AddObserver("LeftButtonPressEvent", sphereCallback)
    sphereWidget.AddObserver("MouseWheelForwardEvent", sphereCallback)
    sphereWidget.AddObserver("MouseWheelBackwardEvent", sphereCallback)


    sphereWidget.AddObserver("MouseMoveEvent", sphereCallback)
    sphereWidget.AddObserver("CursorChangedEvent", sphereCallback)

    sphereWidget.AddObserver("SelectionChangedEvent", sphereCallback)

    sphereWidget.AddObserver("SwipeEvent", sphereCallback)
    
    sphereWidget.AddObserver("LeftButtonDoubleClickEvent", sphereCallback)
    sphereWidget.AddObserver("MiddleButtonDoubleClickEvent", sphereCallback)
    sphereWidget.AddObserver("RightButtonDoubleClickEvent", sphereCallback)
    sphereWidget.AddObserver("MouseWheelLeftEvent", sphereCallback)
    sphereWidget.AddObserver("MouseWheelRightEvent", sphereCallback)

    call_back_command = vtk.vtkCallbackCommand()

    call_back_command.SetCallback(sphereCallback)

    # Start
    interactor.Initialize()
    renwin.Render()
    sphereWidget.On()
    interactor.Start()


if __name__ == '__main__':
    main()
