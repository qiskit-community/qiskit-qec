import math 
import vtk
 
# Let’s make a sixth-order tetrahedron 
order = 6
# The number of points for a sixth-order tetrahedron is 
nPoints = int((order + 1) * (order + 2) * (order + 3) / 6)
 
# Create a tetrahedron and set its number of points. Internally, Lagrange cells
 # compute their order according to the number of points they hold. 
tet = vtk.vtkLagrangeTetra() 
tet.GetPointIds().SetNumberOfIds(nPoints) 
tet.GetPoints().SetNumberOfPoints(nPoints)
tet.Initialize()
 
point = [0.,0.,0.]
barycentricIndex = [0, 0, 0, 0]
 
# For each point in the tetrahedron...
for i in range(nPoints):
  # ...we set its id to be equal to its index in the internal point array. 
  tet.GetPointIds().SetId(i, i)
 
  # We compute the barycentric index of the point... 
  tet.ToBarycentricIndex(i, barycentricIndex)
 
  # ...and scale it to unity.
  for j in range(3):
      point[j] = float(barycentricIndex[j]) / order
 
  # A tetrahedron comprised of the above-defined points has straight
  # edges.
  tet.GetPoints().SetPoint(i, point[0], point[1], point[2])
 
# Add the tetrahedron to a cell array 
tets = vtk.vtkCellArray() 
tets.InsertNextCell(tet)
 
# Add the points and tetrahedron to an unstructured grid 
uGrid =vtk.vtkUnstructuredGrid() 
uGrid.SetPoints(tet.GetPoints())
uGrid.InsertNextCell(tet.GetCellType(), tet.GetPointIds())
 
# Visualize
mapper = vtk.vtkDataSetMapper() 
mapper.SetInputData(uGrid)
 
actor = vtk.vtkActor() 
actor.SetMapper(mapper)
 
renderer = vtk.vtkRenderer() 
renderWindow = vtk.vtkRenderWindow() 
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
style = vtk.vtkInteractorStyleTrackballCamera()
renderWindowInteractor.SetInteractorStyle(style)
renderWindowInteractor.SetRenderWindow(renderWindow)
 
renderer.AddActor(actor) 
renderer.SetBackground(.2, .3, .4)
 
renderWindow.Render() 
renderWindowInteractor.Start()

