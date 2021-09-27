#!/usr/bin/env python

# This example works to create a unstructured grid of cells. Allow the user to pick a cell
# and then change the color of that cell. Colors are encoded in the CellData scalars.

# Issue. Changing the scalar is not currenlty changing the cell color once we are in the
# interactor loop

import vtk

# Catch mouse events.
class MouseInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, data=None, actor=None):
        self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)
        self.data = data
        #self.actor = actor
        

    def leftButtonPressEvent(self, obj, event):
        colors = vtk.vtkNamedColors()

        # Get the location of the click (in window coordinates)
        pos = self.GetInteractor().GetEventPosition()

        picker = vtk.vtkCellPicker()
        picker.SetTolerance(0.0005)

        # Pick from this location.
        picker.Pick(pos[0], pos[1], 0, self.GetDefaultRenderer())

        worldPosition = picker.GetPickPosition()
        picked_cell_id = picker.GetCellId()
        print(f"Cell id is: {picked_cell_id}")

        if picker.GetCellId() != -1:
            print(f"Pick position is: {worldPosition[0]}, {worldPosition[1]}, {worldPosition[2]}")
            # Change the color of the cell to white (=4 from LUT)
            self.data.InsertTuple(picked_cell_id, (4,))
        
        # Update
        self.data.Modified()

        # Forward events
        self.OnLeftButtonDown()

def main():
    colors = vtk.vtkNamedColors()

    # Create a simple grid of four quad cells

    # Data for colored grid
    x = [[0,0,0],[0,1,0],[0,2,0],[1,0,0],[1,1,0],[1,2,0],[2,0,0],[2,1,0],[2,2,0]]
    pts = [[0,3,4,1],[4,7,8,5],[1,4,5,2],[3,6,7,4]]
    scalars = [2,2,1,1]

    # Color array for Cell scalar data
    color_array = vtk.vtkFloatArray()
    for item in scalars:
        i = color_array.InsertNextTuple((item,))

    points = vtk.vtkPoints()
    for i in range(0, len(x)):
        points.InsertPoint(i, x[i])

    ugrid = vtk.vtkUnstructuredGrid()

    for i in range(4):
        ugrid.InsertNextCell(vtk.VTK_QUAD, 4, pts[i])

    ugrid.SetPoints(points)

    ugrid.GetCellData().SetScalars(color_array)

    # Create the LUT for scalars to colors
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(5)
    lut.SetTableRange(0,4)
    lut.Build()
    lut.SetTableValue(0, 0, 0, 0, 1) # Black
    lut.SetTableValue(1, 1, 0, 0, 1) # Red
    lut.SetTableValue(2, 0, 1, 0, 1) # Green
    lut.SetTableValue(3, 0, 0, 1, 1) # Blue
    lut.SetTableValue(4, 1, 1, 1, 1) # White 

    ugrid_mapper = vtk.vtkDataSetMapper()

    ugrid_mapper.SetScalarModeToUseCellData()
    ugrid_mapper.UseLookupTableScalarRangeOn()
    ugrid_mapper.SetLookupTable(lut)

    ugrid_mapper.SetInputData(ugrid)

    ugrid_actor = vtk.vtkActor()
    ugrid_actor.SetMapper(ugrid_mapper)
    ugrid_actor.GetProperty().EdgeVisibilityOn()

    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetWindowName("CellPicking")

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    # Set the custom style to use for interaction.
    style = MouseInteractorStyle()
    style.data = ugrid.GetCellData().GetScalars()
    style.SetDefaultRenderer(renderer)

    renderWindowInteractor.SetInteractorStyle(style)

    renderer.AddActor(ugrid_actor)
    renderer.ResetCamera()

    renderer.SetBackground(colors.GetColor3d("Blue"))

    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__ == '__main__':
    main()