import vtkmodules.all as vtk

colors = vtk.vtkNamedColors()
BACKGROUND_COLOR = colors.GetColor3d('SteelBlue')

class MyCode():
        def __init__(self) -> None:
                # Set Default colors
                self.colorX = 'Gray'
                self.colorZ = 'White'
                # Unstructured grid to store stabilizers
                self.grid = vtk.vtkUnstructuredGrid()
                # Note that these Arrays do not need to be class variables directly
                # Point Data
                self.points = vtk.vtkPoints()
                self.points.Allocate(1000)
        
                self.grid.SetPoints(self.points)
                # What type of point is it? One of
                #
                # NON_QUBIT      = -1
                # INACTIVE_QUBIT = 0
                # DATA_QUBIT     = 1
                # ANCILLA_QUBIT  = 2
                # FLAG_QUBIT     = 3
                self.point_type_name = "Point type"
                self.point_type = vtk.vtkTypeInt32Array()
                self.point_type.SetNumberOfComponents(1)
                self.point_type.SetName(self.point_type_name)
                #self.point_type.Allocate(500)
                self.grid.GetPointData().AddArray(self.point_type)
                # Cell Data: Is the face a stabilizer or not? One of
                # X              = 0
                # Z              = 1
                self.cell_stabilizer_type_name = "Stabilizer type"
                self.cell_stabilizer_type = vtk.vtkTypeInt32Array()
                self.cell_stabilizer_type.SetNumberOfComponents(1)
                self.cell_stabilizer_type.SetName(self.cell_stabilizer_type_name)
                #self.cell_stabilizer_type.Allocate(500)
                self.grid.GetCellData().AddArray(self.cell_stabilizer_type)
        
                # Cell/Stabilizer weight (i.e. Number of qubits defingin stabilizer)
                self.cell_weight_name = "Cell weight"
                self.cell_weight = vtk.vtkTypeInt32Array()
                self.cell_weight.SetNumberOfComponents(1)
                self.cell_weight.SetName(self.cell_weight_name)
                #self.grid.GetCellData().Allocate(500)
                self.grid.GetCellData().AddArray(self.cell_weight)
                # Color of cell integers 0,1,2,3 ... using self.lut to translated to RGB color data
                self.cell_colors_name = "Cell colors"
                self.cell_colors = vtk.vtkFloatArray()
                self.cell_colors.SetNumberOfComponents(1)
                self.cell_colors.SetName(self.cell_colors_name)
                self.grid.GetCellData().AddArray(self.cell_colors)
                

if __name__ == '__main__':
        print("mess that you made me")
        myC = MyCode();

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
        ugridMapper.SetInputData(myC.grid)

        ugridActor = vtk.vtkActor()
        ugridActor.SetMapper(ugridMapper)
        ugridActor.GetProperty().SetColor(colors.GetColor3d('Peacock'))
        ugridActor.GetProperty().EdgeVisibilityOn()

        renderer.AddActor(ugridActor)

        interactor.Initialize()
        renwin.Render()
        interactor.Start()


        


