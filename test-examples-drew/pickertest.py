#!/usr/bin/env python

import vtk
import math

# Catch mouse events.
class MouseInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, data):
        self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)
        self.Data = data
        self.selectedMapper = vtk.vtkDataSetMapper()
        self.selectedActor = vtk.vtkActor()
        super().__init__()

    def leftButtonPressEvent(self, obj, event):
        colors = vtk.vtkNamedColors()

        # Get the location of the click (in window coordinates)
        pos = self.GetInteractor().GetEventPosition()

        picker = vtk.vtkCellPicker()
        picker.SetTolerance(0.0005)

        # Pick from this location.
        picker.Pick(pos[0], pos[1], 0, self.GetDefaultRenderer())

        worldPosition = picker.GetPickPosition()
        print(f"Cell id is: {picker.GetCellId()}")

        if picker.GetCellId() != -1:
            print(f"Pick position is: {worldPosition[0]}, {worldPosition[1]}, {worldPosition[2]}")

            ids = vtk.vtkIdTypeArray()
            ids.SetNumberOfComponents(1)
            ids.InsertNextValue(picker.GetCellId())

            selectionNode = vtk.vtkSelectionNode()
            selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL)
            selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
            selectionNode.SetSelectionList(ids)

            selection = vtk.vtkSelection()
            selection.AddNode(selectionNode)

            extractSelection = vtk.vtkExtractSelection()
            extractSelection.SetInputData(0, self.Data)
            extractSelection.SetInputData(1, selection)
            extractSelection.Update()

            # In selection
            # selected = vtk.vtkPolyData()
            # selected.ShallowCopy(extractSelection.GetOutput())

            print(f"There are  {extractSelection.GetOutput().GetNumberOfPoints()}  points in the selection.")
            print(f"There are  {extractSelection.GetOutput().GetNumberOfCells()}  cells in the selection.")

            self.selectedMapper.SetInputData(extractSelection.GetOutput())
            self.selectedActor.SetMapper(self.selectedMapper)
            self.selectedActor.GetProperty().EdgeVisibilityOn()
            self.selectedActor.GetProperty().SetColor(colors.GetColor3d("Red"))

            self.selectedActor.GetProperty().SetLineWidth(3)

            self.GetInteractor().GetRenderWindow().GetRenderers().GetFirstRenderer().AddActor(self.selectedActor)

        # Forward events
        self.OnLeftButtonDown()

def main():
    colors = vtk.vtkNamedColors()

    planeSource = vtk.vtkPlaneSource()
    planeSource.Update()

    triangleFilter = vtk.vtkTriangleFilter()
    triangleFilter.SetInputConnection(planeSource.GetOutputPort())
    triangleFilter.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(triangleFilter.GetOutputPort())

    actor = vtk.vtkActor()
    actor.GetProperty().SetColor(colors.GetColor3d("Green"))
    actor.SetMapper(mapper)

    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetWindowName("CellPicking")

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    # Set the custom stype to use for interaction.
    #style = MouseInteractorStyle(triangleFilter.GetOutput())
    style = MouseInteractorStyle()
    style.SetDefaultRenderer(renderer)
    style.Data = triangleFilter.GetOutput()

    renderWindowInteractor.SetInteractorStyle(style)

    renderer.AddActor(actor)

    renderer.ResetCamera()

    renderer.SetBackground(colors.GetColor3d("Blue"))

    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__ == '__main__':
    main()