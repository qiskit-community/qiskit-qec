import vtk
# Define interaction style
class InteractorStyleMoveVertex(vtk.vtkInteractorStyleTrackballActor):
    def __init__(self):
        self.color = vtk.vtkNamedColors()
        self.Data = vtk.vtkUnstructuredGrid()
        self.GlyphData = vtk.vtkPolyData()
        self.MoveMapper = vtk.vtkPolyDataMapper()
        self.MoveActor = vtk.vtkActor()
        self.MovePolyData = vtk.vtkPolyData()
        self.MoveGlyphFilter = vtk.vtkVertexGlyphFilter()
        self.PointPicker = vtk.vtkPointPicker()
        self.Move = False
        self.SelectedPoint = 0
        self.AddObserver("MouseMoveEvent", self.MyOnMouseMove)
        self.AddObserver("LeftButtonPressEvent", self.MyOnMiddleButtonUp)
        self.AddObserver("LeftButtonReleaseEvent", self.MyOnMiddleButtonDown)
    def InteractorStyleMoveVertex(self):
        self.Move = False
        # Setup ghost glyph
        points = vtk.vtkPoints()
        points.InsertNextPoint(0,0,0)
        self.MovePolyData = vtk.vtkPolyData()
        self.MovePolyData.SetPoints(points)
        self.MoveGlyphFilter = vtk.vtkVertexGlyphFilter()
        self.MoveGlyphFilter.SetInputData(self.MovePolyData)
        self.MoveGlyphFilter.Update()
        self.MoveMapper = vtk.vtkPolyDataMapper()
        self.MoveMapper.SetInputConnection(
            self.MoveGlyphFilter.GetOutputPort())
        self.MoveActor = vtk.vtkActor()
        self.MoveActor.SetMapper(self.MoveMapper)
        self.MoveActor.VisibilityOff()
        self.MoveActor.GetProperty().SetPointSize(10)
        self.MoveActor.GetProperty().SetColor(self.color.GetColor3d("Pink"))
    def MyOnMouseMove(self, obj, event):
        if not self.Move:
            return
        # Forward events
        self.OnMouseMove()
    def MyOnMiddleButtonUp(self, obj, event):
        self.EndPan()
        self.Move = False
        self.MoveActor.VisibilityOff()
        self.Data.GetPoints().SetPoint(self.SelectedPoint,self.MoveActor.GetPosition())
        self.Data.Modified()
        self.GetCurrentRenderer().Render()
        self.GetCurrentRenderer().GetRenderWindow().Render()
    def MyOnMiddleButtonDown(self, obj, event):
        # Get the selected point
        x = self.GetInteractor().GetEventPosition()[0]
        y = self.GetInteractor().GetEventPosition()[1]
        self.FindPokedRenderer(x, y)
        self.PointPicker.Pick(self.GetInteractor().GetEventPosition()[0],
                                self.GetInteractor().GetEventPosition()[1],
                                0, # always zero.
                                self.GetInteractor().GetRenderWindow()
                                    .GetRenderers()
                                    .GetFirstRenderer())
        if (self.PointPicker.GetPointId() >= 0):
            self.StartPan()
            self.MoveActor.VisibilityOn()
            self.Move = True
            self.SelectedPoint = self.PointPicker.GetPointId()
            print(f"Dragging point {self.SelectedPoint}")
            p=[0,0,0]
            self.Data.GetPoint(self.SelectedPoint, p)
            print(f"p: {p[0]} {p[1]} {p[2]}")
            self.MoveActor.SetPosition(p)
            self.GetCurrentRenderer().AddActor(self.MoveActor)
            self.InteractionProp = self.MoveActor
def main():
    color = vtk.vtkNamedColors()
    points = vtk.vtkPoints()
    points.InsertNextPoint(0, 0, 0)
    points.InsertNextPoint(1, 0, 0)
    points.InsertNextPoint(2, 0, 0)
    input = vtk.vtkUnstructuredGrid()
    input.SetPoints(points)
    glyphFilter = vtk.vtkVertexGlyphFilter()
    glyphFilter.SetInputData(input)
    glyphFilter.Update()
    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyphFilter.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(10)
    actor.GetProperty().SetColor(color.GetColor3d("Tomato"))
    # Visualize
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetWindowName("MoveAVertexUnstructuredGrid")
    renderWindowInteractor =vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderer.AddActor(actor)
    renderer.SetBackground(color.GetColor3d("Gray"))
    renderWindow.Render()
    style = InteractorStyleMoveVertex()
    renderWindowInteractor.SetInteractorStyle(style)
    style.Data = input
    style.GlyphData = glyphFilter.GetOutput()
    renderer.GetActiveCamera().Zoom(0.9)
    renderWindowInteractor.Start()
if __name__ == '__main__':
    main()