import vtk

class HighlightInteractorStyle(vtk.vtkInteractorStyleRubberBandPick):
    VTKISRBP_ORIENT = 0
    VTKISRBP_SELECT = 1

    def __init__(self):
        self.PolyData = vtk.vtkPolyData()
        self.SelectedActor = vtk.vtkActor()
        self.SelectedMapper = vtk.vtkDataSetMapper()
        self.SelectedActor.SetMapper(self.SelectedMapper)
        self.AddObserver("OnLeftButtonMouseUpEvent", self.OnLeftButtonMouseUpEvent)

    def OnLeftButtonMouseUpEvent(self, obj, event):
        # Forward events
        self.OnLeftButtonUp()

        if (self.CurrentMode == self.VTKISRBP_SELECT):
            colors = vtk.vtkNamedColors

            frustum = self.GetInteractor().GetPicker().GetFrustum()

            extractPolyDataGeometry = vtk.vtkExtractPolyDataGeometry()
            extractPolyDataGeometry.SetInputData(self.PolyData)
            extractPolyDataGeometry.SetImplicitFunction(frustum)
            extractPolyDataGeometry.Update()

            print(f"Extracted {extractPolyDataGeometry.GetOutput().GetNumberOfCells()} cells")
            
            self.SelectedMapper.SetInputData(extractPolyDataGeometry.GetOutput())

            self.SelectedActor.GetProperty().SetColor(colors.GetColor3d("Tomato"))
            self.SelectedActor.GetProperty().SetPointSize(5)

            self.GetInteractor().GetRenderWindow().GetRenderers().GetFirstRenderer().AddActor(self.SelectedActor)

            self.GetInteractor().GetRenderWindow().Render()
            self.HighlightProp(None)

    def setPolyData(self, data):
        self.PolyData=data


def main(argc, argv):

    
    polydata  ReadPolyData()



