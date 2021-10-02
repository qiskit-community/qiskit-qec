#!/usr/bin/env python

# This example shows how to visualize polygons, convex or not.

import vtk

# Define a set of points - these are the ordered polygon vertices
polygonPoints = vtk.vtkPoints()
polygonPoints.SetNumberOfPoints(6)
polygonPoints.InsertPoint(0, 0, 0, 0)
polygonPoints.InsertPoint(1,.4,.4, 0)
polygonPoints.InsertPoint(2, 1, 0, 0)
polygonPoints.InsertPoint(3, 1, 1, 0)
polygonPoints.InsertPoint(4,.1,.7, 0)
polygonPoints.InsertPoint(5, 0, 1, 0)

# Make a cell with these points
aPolygon = vtk.vtkPolygon()
aPolygon.GetPointIds().SetNumberOfIds(6)
aPolygon.GetPointIds().SetId(0, 0)
aPolygon.GetPointIds().SetId(1, 1)
aPolygon.GetPointIds().SetId(2, 2)
aPolygon.GetPointIds().SetId(3, 3)
aPolygon.GetPointIds().SetId(4, 4)
aPolygon.GetPointIds().SetId(5, 5)

# The cell is put into a mesh (containing only one cell)
aPolygonGrid = vtk.vtkUnstructuredGrid()
aPolygonGrid.Allocate(1, 1)
aPolygonGrid.InsertNextCell(aPolygon.GetCellType(), aPolygon.GetPointIds())
aPolygonGrid.SetPoints(polygonPoints)

# This part is needed for non-convex polygon rendering
#aPolygonGeomFilter = vtk.vtkGeometryFilter()
#aPolygonGeomFilter.SetInput(aPolygonGrid)
aPolygonTriangleFilter = vtk.vtkTriangleFilter()
aPolygonTriangleFilter.SetInput(aPolygonGrid.GetOutput())
#
# This one is only to check the triangulation (when factor < 1)
aPolygonShrinkFilter = vtk.vtkShrinkFilter()
aPolygonShrinkFilter.SetShrinkFactor( 0.9 )
#aPolygonShrinkFilter.SetShrinkFactor( 1.0 )
aPolygonShrinkFilter.SetInput( aPolygonGrid)

# Make ready for rendering
aPolygonMapper = vtk.vtkDataSetMapper()
aPolygonMapper.SetInput(aPolygonShrinkFilter.GetOutput())
aPolygonActor = vtk.vtkActor()
aPolygonActor.SetMapper(aPolygonMapper)
aPolygonActor.GetProperty().SetDiffuseColor(1, .4, .5)

# Create the usual rendering stuff.
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(300, 150)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

ren.SetBackground(.1, .2, .4)
ren.AddActor(aPolygonActor)
ren.ResetCamera()

# Render the scene and start interaction.
iren.Initialize()
renWin.Render()
iren.Start()