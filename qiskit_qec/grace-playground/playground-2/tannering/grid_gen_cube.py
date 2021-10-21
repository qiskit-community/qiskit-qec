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
qcolormap = {DATA:colors.GetColor3d('Black')}

ZGROUP = 0
XGROUP  = 1
gcolormap = {ZGROUP: colors.GetColor3d('Red'), XGROUP: colors.GetColor3d('Yellow')}
class Group(vtk.vtkActor):
    def __init__(self, gtype,qubits, *args, **kwargs):
        super().__init__()
        self.gtype = gtype
        self.qubits = qubits
        

class Qubit(vtk.vtkActor):
    def __init__(self, qutype, *args, **kwargs):
        super().__init__()
        self.qutype = qutype


class MouseInteractorHighLightActor(vtk.vtkInteractorStyleTrackballCamera):
    
    def __init__(self, parent=None):
        self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)
        self.AddObserver("KeyPressEvent", self.keyPressEvent)
        
    def leftButtonPressEvent(self, obj, event):
        print("my bad ....")
        clickPos = self.GetInteractor().GetEventPosition()
        
        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[ 0 ], clickPos[ 1 ], 0, self.GetDefaultRenderer())
        
        # get the new
        self.NewPickedActor = picker.GetActor()
        
        # If something was selected
        if isinstance(self.NewPickedActor, Qubit):
            selected_qubits.append(self.NewPickedActor)
            self.NewPickedActor.GetProperty().SetColor(colors.GetColor3d('White'))
        
        self.OnLeftButtonDown()
        return
    
    def keyPressEvent(self, obj, event):
        key = self.GetInteractor().GetKeyCode()
        print(f"key is: {key}")
        if key == 'z' or key == 'x':
            xsum = 0
            ysum = 0
            for q in selected_qubits:
                pos = q.GetPosition()
                print(f"Cur pos is : {pos}")
                xsum += pos[0]
                ysum += pos[1]
            xpos = xsum/len(selected_qubits)
            ypos = ysum/len(selected_qubits)
            zpos = 1
            source = vtk.vtkCubeSource()
            source.SetXLength(0.25)
            source.SetYLength(0.25)
            source.SetZLength(0.25)
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())
            
            gqs = []
            for qubits in selected_qubits: # can't deep copy bc copies qubits, and can't shallow bc qubits will get cleared rip
                gqs.append(qubits)
            if key == 'z':
                group = Group(ZGROUP, gqs)
                group.GetProperty().SetColor(colors.GetColor3d('Yellow'))
            elif key == 'x':
                group = Group(XGROUP, gqs)
                group.GetProperty().SetColor(colors.GetColor3d('Red'))

                
            group.SetMapper(mapper)
            group.GetProperty().SetDiffuseColor(colors.GetColor3d('White'))
            group.GetProperty().SetDiffuse(.8)
            group.GetProperty().SetSpecular(.5)
            group.GetProperty().SetSpecularColor(colors.GetColor3d('White'))
            group.GetProperty().SetSpecularPower(30.0)
            group.GetProperty().SetColor(gcolormap[group.gtype])
            

            group.SetPosition(xpos, ypos, zpos)
            self.GetDefaultRenderer().AddActor(group)
            
            for qubit in selected_qubits:
                qubit.GetProperty().SetColor(qcolormap[qubit.qutype])
                
            selected_qubits.clear()
        
            
            self.GetDefaultRenderer().GetRenderWindow().Render() # make sure screen updates w/ cube: https://stackoverflow.com/questions/48573811/vtk-view-doesnt-update-until-after-user-interaction
        self.OnKeyPress()






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
    
    # Add spheres to play with
    for i in range(4):
        for j in range(4):
            source = vtk.vtkSphereSource()
            source.SetRadius(0.25)
            source.SetCenter(0, 0, 0)
            source.SetPhiResolution(11)
            source.SetThetaResolution(21)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())
            
            actor = Qubit(DATA)
            actor.SetMapper(mapper)
            
            
            actor.GetProperty().SetDiffuseColor(colors.GetColor3d('Black'))
            actor.GetProperty().SetDiffuse(.8)
            actor.GetProperty().SetSpecular(.5)
            actor.GetProperty().SetSpecularColor(colors.GetColor3d('Black'))
            actor.GetProperty().SetColor(qcolormap[actor.qutype])
            actor.GetProperty().SetSpecularPower(30.0)
            
            actor.SetPosition(i, j, 0)
            renderer.AddActor(actor)
    
    # Start
    interactor.Initialize()
    renwin.Render()
    interactor.Start()


if __name__ == '__main__':
    main()
