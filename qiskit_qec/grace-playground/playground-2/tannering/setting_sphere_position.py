import vtk

colors = vtk.vtkNamedColors()
NUMBER_OF_SPHERES = 10

# have qubits inherit from vtkActor?
# have them contain their own source?
# have a separate map: Data_qubit: BLACK
# have each code-actor contain it's type?
# source controls location oddly enough
# mapper contains color
# where to store metadata?
import random

selected_qubits = [ ]
# class Code():
#     def __init__(self):
#         qubits = []
#         groups = []
#
#
# class CodePiece(vtk.vtkActor):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#
# class Qubit(CodePiece):
#     def __init__(self, qid, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.id = qid
#
#     def createGraphicSource(self, pos_x, pos_y, pos_z):
#         self.graphic_source = vtk.vtkSphereSource()
#         self.graphic_source.SetCenter(pos_x, pos_y, pos_z)
#         return self.graphic_source
#
#
#

# class Group(CodePiece):
#     def __init__(self, gid, qubits=None, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.id = gid
#         self.qubits = qubits
#
#
#     def createGraphicSource(self):
#         if self.qubits is None:
#             return None
#
#         self.graphic_source = vtk.vtkCubeSource()
#         return self.graphic_source
#

DATA = 0


class Qubit(vtk.vtkActor):
    def __init__(self, qtype, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.qtype = qtype


class MouseInteractorHighLightActor(vtk.vtkInteractorStyleTrackballCamera):
    
    def __init__(self, parent=None):
        self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)
        self.AddObserver("KeyPressEvent", self.keyPressEvent)
        
        self.LastPickedActor = None
        self.LastPickedProperty = vtk.vtkProperty()
    
    def leftButtonPressEvent(self, obj, event):
        print("my bad ....")
        clickPos = self.GetInteractor().GetEventPosition()
        
        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[ 0 ], clickPos[ 1 ], 0, self.GetDefaultRenderer())
        
        # get the new
        self.NewPickedActor = picker.GetActor()
        
        # If something was selected
        if self.NewPickedActor:
            selected_qubits.append(self.NewPickedActor)
            x = random.randrange(-3, 3, 1)
            y = random.randrange(-3, 3, 1)
            z = random.randrange(-3, 3, 1)

            print(f"setting new pos to {x, y, z}")
            self.NewPickedActor.SetPosition(x,y,z)
            # # If we picked something before, reset its property
            if self.LastPickedActor:
                self.LastPickedActor.GetProperty().DeepCopy(self.LastPickedProperty)
            # Save the property of the picked actor so that we can
            # restore it next time
            self.LastPickedProperty.DeepCopy(self.NewPickedActor.GetProperty())
            # Highlight the picked actor by changing its properties
            self.NewPickedActor.GetProperty().SetColor(colors.GetColor3d('black'))
            
            # save the last picked actor
            self.LastPickedActor = self.NewPickedActor
        
        self.OnLeftButtonDown()
        return
    
    def keyPressEvent(self, obj, event):
        key = self.GetInteractor().GetKeyCode()
        print(f"key is: {key}")
        if key == 'g':
            print(key)
        print("We took the long way down")
        for q in selected_qubits:
            sour = q.GetProperty()
            print(f"my bad habits lead to {sour}\n burn till the fire runs out: ")
            # GetProperty().SetColor(colors.GetColor3d('black'))


def main():
    # A renderer and render window
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d('SteelBlue'))
    
    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(renderer)
    renwin.SetSize(640, 480)
    renwin.SetWindowName('HighlightPickedActor')
    
    # An interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renwin)
    
    # add the custom style
    style = MouseInteractorHighLightActor()
    style.SetDefaultRenderer(renderer)
    interactor.SetInteractorStyle(style)
    
    randomSequence = vtk.vtkMinimalStandardRandomSequence()
    # randomSequence.SetSeed(1043618065)
    # randomSequence.SetSeed(5170)
    randomSequence.SetSeed(8775070)
    # Add spheres to play with
    for i in range(NUMBER_OF_SPHERES):
        source = vtk.vtkSphereSource()
        
        # random radius

        radius = randomSequence.GetRangeValue(0.5, 1.0)
        randomSequence.Next()
        
        source.SetRadius(radius)
        source.SetCenter(0,0,0)
        source.SetPhiResolution(11)
        source.SetThetaResolution(21)
        
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        
        r = randomSequence.GetRangeValue(0.4, 1.0)
        randomSequence.Next()
        g = randomSequence.GetRangeValue(0.4, 1.0)
        randomSequence.Next()
        b = randomSequence.GetRangeValue(0.4, 1.0)
        randomSequence.Next()
        
        actor.GetProperty().SetDiffuseColor(r, g, b)
        actor.GetProperty().SetDiffuse(.8)
        actor.GetProperty().SetSpecular(.5)
        actor.GetProperty().SetSpecularColor(colors.GetColor3d('White'))
        actor.GetProperty().SetSpecularPower(30.0)
        
        renderer.AddActor(actor)
    
    # Start
    interactor.Initialize()
    renwin.Render()
    interactor.Start()


if __name__ == '__main__':
    main()
