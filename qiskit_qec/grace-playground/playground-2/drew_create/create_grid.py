import vtkmodules.all as vtk

colors = vtk.vtkNamedColors()
BACKGROUND_COLOR = colors.GetColor3d('SteelBlue')

def runme():

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName("RSC-3x3.vtu")
    reader.Update()
    grid = reader.GetOutput()
    
    print(f"they show their truth 1 single time: {grid}")

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(BACKGROUND_COLOR)

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(renderer)
    renwin.SetSize(640, 480)
    renwin.SetWindowName('HighlightPickedActor')

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renwin)

    # add the custom style
    style = vtk.vtkInteractorStyleTrackballCamera()
    style.SetDefaultRenderer(renderer)
    interactor.SetInteractorStyle(style)

    ugridMapper = vtk.vtkDataSetMapper()
    ugridMapper.SetInputData( grid)

    ugridActor = vtk.vtkActor()
    ugridActor.SetMapper(ugridMapper)
    ugridActor.GetProperty().SetColor(colors.GetColor3d('Peacock'))
    ugridActor.GetProperty().EdgeVisibilityOn()

    renderer.AddActor(ugridActor)

    interactor.Initialize()
    renwin.Render()
    interactor.Start()


if __name__ == '__main__':
    runme()