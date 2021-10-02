import vtk

class main:
  colors = vtk.vtkNamedColors()

  # Create three points. We will join (Origin and P0) with a red line and
  # (Origin and P1) with a green line
  origin = [0.0, 0.0, 0.0]
  p0 = [1.0, 0.0, 0.0]
  p1 = [0.0, 1.0, 0.0]
  p2 = [0.0, 1.0, 2.0]
  p3 = [1.0, 2.0, 3.0]

  q0 = [0,0,0]
  qm = [0.5,0.25,0]
  q1 = [1,0,0]

  # Create a vtkPoints object and store the points in it
  points = vtk.vtkPoints()
  #points.InsertNextPoint(origin)
  #points.InsertNextPoint(p0)
  #points.InsertNextPoint(p1)
  #points.InsertNextPoint(p2)
  #points.InsertNextPoint(p3)

  points.InsertNextPoint(q0)
  points.InsertNextPoint(qm)
  points.InsertNextPoint(q1)

  spline = vtk.vtkParametricSpline()
  spline.SetPoints(points)

  functionSource = vtk.vtkParametricFunctionSource() 
  functionSource.SetParametricFunction(spline)
  functionSource.Update()

  sphere = vtk.vtkSphereSource()
  sphere.SetPhiResolution(21)
  sphere.SetThetaResolution(21)
  sphere.SetRadius(0.025)

  # Setup actor and mapper
  mapper= vtk.vtkPolyDataMapper()
  mapper.SetInputConnection(functionSource.GetOutputPort())

  actor= vtk.vtkActor()
  actor.SetMapper(mapper)
  actor.GetProperty().SetColor(colors.GetColor3d("DarkSlateGrey"))
  actor.GetProperty().SetLineWidth(3.0)

  actor.GetProperty().SetEdgeVisibility(1)
  actor.GetProperty().SetEdgeColor(0.9,0.9,0.4)
  actor.GetProperty().SetLineWidth(6)
  actor.GetProperty().SetPointSize(12)
  actor.GetProperty().SetRenderLinesAsTubes(1)
  actor.GetProperty().SetRenderPointsAsSpheres(1)
  actor.GetProperty().SetVertexVisibility(1)
  actor.GetProperty().SetVertexColor(0.5,1.0,0.8)

  # Create a polydata to store everything in
  polyData = vtk.vtkPolyData()
  polyData.SetPoints(points)

  pointMapper = vtk.vtkGlyph3DMapper()
  pointMapper.SetInputData(polyData)
  pointMapper.SetSourceConnection(sphere.GetOutputPort())

  pointActor = vtk.vtkActor()
  pointActor.SetMapper(pointMapper)
  pointActor.GetProperty().SetColor(colors.GetColor3d("Peacock"))

  # Setup render window, renderer, and interactor
  renderer = vtk.vtkRenderer()
  renderWindow = vtk.vtkRenderWindow()
  renderWindow.AddRenderer(renderer)
  renderWindow.SetWindowName("ParametricSpline")

  renderWindowInteractor = vtk.vtkRenderWindowInteractor()
  renderWindowInteractor.SetRenderWindow(renderWindow)
  renderer.AddActor(actor)
  renderer.AddActor(pointActor)
  renderer.SetBackground(colors.GetColor3d("Silver"))

  renderWindow.Render()
  renderWindowInteractor.Start()

if __name__ == '__main__':
    main()