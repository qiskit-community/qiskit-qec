import vtkmodules.all as vtk

from qiskit_qec.test_examples_grace.tannering.bad_habits import BADHABITS

import random
colors = vtk.vtkNamedColors()
NUM_GRID = 4

selected_qubits = [ ]
selected_stabilizers = []
old_qubit_actors = set()
old_stabilizer_actors = []

pauli_types = {0:colors.GetColor3d('Yellow'), 1: colors.GetColor3d('Green') }

CUBELEN = 0.25
BACKGROUND_COLOR =  colors.GetColor3d('SteelBlue')
randomSequence = vtk.vtkMinimalStandardRandomSequence()
# randomSequence.SetSeed(1043618065)
# randomSequence.SetSeed(5170)
randomSequence.SetSeed(8775070)

class Qubit():
    radius = 0.1

    def __init__(self, id, qubit_graphics):
        self.id = id
        self.qubit_graphics = qubit_graphics
        self.stab_ids = set()
    
    def add_actors(self, actors):
        if self.qubit_graphics is None:
            self.qubit_graphics = actors
        self.qubit_graphics += actors
        
    def get_id(self):
        return self.id
    
    def add_stab_ids(self, stab_id_list):
        for sid in stab_id_list:
            self.stab_ids.add(sid)
        
class QubitActor(vtk.vtkActor):
    def __init__(self, parent_id):
        super().__init__()
        self.parent_id = parent_id
        
    def get_parent_id(self):
        return self.parent_id
    
class Stabilizer():  # should this be update-able?
    def __init__(self, qubit_ids, id, stabilizer_graphics):
        self.qubits = qubit_ids
        self.id = id
        self.stabilizer_graphics = stabilizer_graphics
        self.pauli_type = 0
        self.tagline = BADHABITS[random.randint(0, len(BADHABITS) - 1)]
    
    def add_actors(self, actors):
        if self.stabilizer_graphics is None:
            self.stabilizer_graphics = actors
        self.stabilizer_graphics += actors
        
    def add_qubits(self, more_qubits:[int]):
        self.qubits += more_qubits
        
    def get_id(self):
        return self.id


class StabilizerActor(vtk.vtkActor):
    def __init__(self, parent_id):
        super().__init__()
        self.parent_id = parent_id
        
        
    def get_parent_id(self):
        return self.parent_id

        
class TannerCode():
    # Stabilizer Id = index in stabilizer array
    def __init__(self):
        self.qubits = []
        self.valid_qubits = {}
        self.next_valid_qubit_id = 0

        
        self.stabilizers = []
        self.valid_stabilizers = {} # stabilizerId: valid/Invalid
        self.next_valid_stabilizer_id = 0
        
    def detect_stabilizer_collision(self, cube_actor):
        cur_pos = list(cube_actor.GetPosition())
        for cube in tanner.get_stabilizers():
            if cube.get_id() != cube_actor.get_parent_id() and tanner.is_stabilizer_valid(cube.get_id()):
                print(f"insatbs: cube is {cube}")
                if cube is not None:
                    for graphic in cube.stabilizer_graphics:
                        this_pos = graphic.GetPosition()
                        lensum = 0
                        for x, y in zip(cur_pos, this_pos): lensum += (x - y) ** 2
                        lensum = lensum ** (1 / 2)
                        if lensum < (CUBELEN * 2 * (2 ** 0.5)):
                            cur_pos[ 2 ] = cur_pos[ 2 ] + 1
        
        return cur_pos


    def get_stabilizers(self):
        return self.stabilizers
    
    def get_qubits(self):
        return self.qubits
    
    def is_stabilizer_valid(self, stab_id):
        return self.valid_stabilizers[stab_id]

    
    def create_stabilizer(self, qubits, stabilizer_actors: [StabilizerActor]=None):
        stab_id = self.get_next_stabilizer_id()
        new_stabilizer = Stabilizer(qubits, stab_id, stabilizer_actors)
        self.stabilizers.append(new_stabilizer)
        self.valid_stabilizers[stab_id] = True
        return new_stabilizer
    
    def create_qubit(self, graphics_items=None):
        qid = self.get_next_qubit_id(),
        qubit = Qubit(qid, graphics_items)
        self.qubits.append(qubit)
        self.valid_qubits[qid] = True

        return qubit
    

    def get_next_stabilizer_id(self):
        nid = self.next_valid_stabilizer_id
        self.next_valid_stabilizer_id += 1
        return nid
        
    def get_next_qubit_id(self):
        nid = self.next_valid_qubit_id
        self.next_valid_qubit_id += 1
        return nid

tanner = TannerCode()


class MouseInteractorHighLightActor(vtk.vtkInteractorStyleTrackballCamera):
    
    def __init__(self, parent=None):
        self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)
        self.AddObserver("KeyPressEvent", self.keyPressEvent)
        self.AddObserver("RightButtonPressEvent", self.rightButtonPressEvent)

        
        self.LastPickedActor = None
        self.LastPickedProperty = vtk.vtkProperty()
    
    def leftButtonPressEvent(self, obj, event):
        

        clickPos = self.GetInteractor().GetEventPosition()
        
        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[ 0 ], clickPos[ 1 ], 0, self.GetDefaultRenderer())
        
        # get the new
        self.NewPickedActor = picker.GetActor()
        
        print(f"NewPickedActor: {self.NewPickedActor} and is of type: {type(self.NewPickedActor)} and isinstance qubitactor: {isinstance(self.NewPickedActor, QubitActor)}")
        
        if isinstance(self.NewPickedActor, QubitActor):
            # If something was selected
            if self.NewPickedActor: #if is #qubit
                selected_qubits.append(self.NewPickedActor)
                # # If we picked something before, reset its property
                old_qubit_actor_prop = vtk.vtkProperty()
                old_qubit_actor_prop.DeepCopy(self.NewPickedActor.GetProperty())
                old_qubit_actors.add(self.NewPickedActor)
                # Highlight the picked actor by changing its properties
                self.NewPickedActor.GetProperty().SetColor(colors.GetColor3d('Red'))
                self.NewPickedActor.GetProperty().SetDiffuse(1.0)
                self.NewPickedActor.GetProperty().SetSpecular(0.0)
                self.NewPickedActor.GetProperty().EdgeVisibilityOn()

                # save the last picked actor
                self.LastPickedActor = self.NewPickedActor
        
        elif isinstance(self.NewPickedActor, StabilizerActor):
            # If something was selected
            if self.NewPickedActor: #if is #stabilizer
                selected_stabilizers.append(self.NewPickedActor)
                # # If we picked something before, reset its property
                old_stabilizer_actor_prop = vtk.vtkProperty()
                old_stabilizer_actor_prop.DeepCopy(self.NewPickedActor.GetProperty())
                old_stabilizer_actors.append((self.NewPickedActor, old_stabilizer_actor_prop))
                # Highlight the picked actor by changing its properties
                self.NewPickedActor.GetProperty().SetColor(colors.GetColor3d('Red'))
                self.NewPickedActor.GetProperty().SetDiffuse(1.0)
                self.NewPickedActor.GetProperty().SetSpecular(0.0)
                self.NewPickedActor.GetProperty().EdgeVisibilityOn()

                # save the last picked actor
                self.LastPickedActor = self.NewPickedActor
                
                
        return self.OnLeftButtonDown()


    def rightButtonPressEvent(self, obj, event):
        for qu in tanner.qubits:
            quid = qu.get_id()
            if isinstance(quid, tuple):
                quid = quid[0]
            tanner.valid_qubits[quid] = True
            
            for actor in tanner.qubits[quid].qubit_graphics:
                actor.GetProperty().SetOpacity(1)

                actor.GetProperty().SetDiffuseColor(colors.GetColor3d('White'))
                actor.GetProperty().SetDiffuse(.8)
                actor.GetProperty().SetSpecular(.5)
                actor.GetProperty().SetSpecularColor(colors.GetColor3d('White'))
                actor.GetProperty().SetSpecularPower(30.0)
                actor.GetProperty().SetRepresentationToSurface()
                actor.GetProperty().SetAmbientColor(colors.GetColor3d('White'))
                actor.GetProperty().SetEdgeVisibility(False)

        self.GetDefaultRenderer().Render()
        


    def keyPressEvent(self, obj, event):
        key = self.GetInteractor().GetKeyCode()
        print(f"key is: {key}")
        if key == 'g':
            
            q_ids  = []
            for q in selected_qubits:
                q_ids.append(q.get_parent_id())
                
            stabilizer = tanner.create_stabilizer(q_ids)
            
            for q in selected_qubits:
                quid = q.get_parent_id()
                if isinstance(quid, tuple):
                    quid = quid[0]
                tanner.qubits[quid].add_stab_ids([stabilizer.get_id()])

            stabilizer_actors = []
            
            xsum = 0
            ysum = 0

            print("We took the long way down")
            for q in selected_qubits:
                pos = q.GetCenter()
                
                xsum += pos[ 0 ]
                ysum += pos[ 1 ]
                
            xpos = xsum / len(selected_qubits)
            ypos = ysum / len(selected_qubits)
            zpos = 1
            
            cubesource = vtk.vtkCubeSource()
            cubesource.SetXLength(CUBELEN)
            cubesource.SetYLength(CUBELEN)
            cubesource.SetZLength(CUBELEN)
            cubemapper = vtk.vtkPolyDataMapper()
            cubemapper.SetInputConnection(cubesource.GetOutputPort())
            
            
            cubeActor = StabilizerActor(stabilizer.get_id())
            cubeActor.SetMapper(cubemapper)
            cubeActor.SetPosition(xpos, ypos, zpos)
            
            new_pos = tanner.detect_stabilizer_collision(cubeActor)
            xpos = new_pos[0]
            ypos = new_pos[1]
            zpos = new_pos[2]

            cubeActor.SetPosition(new_pos)
            
            # TODO: make this not O(n^2) :...(
            collide = vtk.vtkCollisionDetectionFilter()
            collide.SetInputConnection(0, cubesource.GetOutputPort())
           
            

            stabilizer_actors.append(cubeActor)
            

            self.GetDefaultRenderer().AddActor(cubeActor)
            
            
            for q in selected_qubits:
                pos = q.GetCenter()
                lineSource = vtk.vtkLineSource()
                lineSource.SetPoint1(xpos, ypos, zpos)
                lineSource.SetPoint2(pos)
                # Setup actor and mapperg
                lineMapper = vtk.vtkPolyDataMapper()
                lineMapper.SetInputConnection(lineSource.GetOutputPort())
                lineActor = StabilizerActor(stabilizer.get_id())
                lineActor.SetMapper(lineMapper)
                lineActor.GetProperty().SetColor(colors.GetColor3d('Red'))
                lineActor.GetProperty().SetLineWidth(4.55)
                stabilizer_actors.append(lineActor)
                self.GetDefaultRenderer().AddActor(lineActor)

            stabilizer.add_actors(stabilizer_actors)

            selected_qubits.clear()
            
            for actor in old_qubit_actors:
                actor.GetProperty().SetDiffuseColor(colors.GetColor3d('White'))
                actor.GetProperty().SetDiffuse(.8)
                actor.GetProperty().SetSpecular(.5)
                actor.GetProperty().SetSpecularColor(colors.GetColor3d('White'))
                actor.GetProperty().SetAmbientColor(colors.GetColor3d('White'))
                actor.GetProperty().EdgeVisibilityOff()
            old_qubit_actors.clear()

        if key == 'd':
            # delete qubits
            # check if qubits are attached to stabilizers
            # remake stabilizer?

            
            for qubit in selected_qubits:
                quid = qubit.get_parent_id()
                if isinstance(quid, tuple):
                    quid = quid[0]
                    for stab_id in tanner.qubits[quid].stab_ids:
                        selected_stabilizers.append(tanner.stabilizers[stab_id].stabilizer_graphics[0])
                if isinstance(quid, tuple):
                    quid = quid[0]
                for actor in tanner.get_qubits()[quid].qubit_graphics:
                    print("seeing actors: in qubitsklasfd")
                    # actor.GetProperty().SetDiffuseColor(BACKGROUND_COLOR)
                    # actor.GetProperty().SetColor(BACKGROUND_COLOR)
                    actor.GetProperty().SetOpacity(0.0)
                    
                    print(f"here are the actors: {actor}")
                    tanner.valid_qubits[quid] = False
                    
                    if quid in old_qubit_actors:
                        print(f"should be removing: actor: {actor}")
                        old_qubit_actors.pop(quid)
                        
                        
            for stab in selected_stabilizers:
                stid = stab.get_parent_id()
                
                for actor in tanner.get_stabilizers()[stid].stabilizer_graphics:
                    self.GetDefaultRenderer().RemoveActor(actor)
            selected_stabilizers.clear()
            
                    
            selected_qubits.clear()
            
        if key == 'y':
            for stab in selected_stabilizers:
                tanner.stabilizers[stab.get_parent_id()].pauli_type = 0
                stab.GetProperty().SetColor(pauli_types[0])
            selected_stabilizers.clear()
        if key == 'x':
            for stab in selected_stabilizers:
                tanner.stabilizers[stab.get_parent_id()].pauli_type = 1
                stab.GetProperty().SetColor(pauli_types[1])
            selected_stabilizers.clear()
            
            
        if key == 'a':
            if len(selected_stabilizers) == 1 and len(selected_qubits) > 0:
                stab_g = selected_stabilizers[0]
                stab = tanner.stabilizers[stab_g.get_parent_id()]
                for qu in selected_qubits:
                    quid = qu.get_parent_id()
                    if isinstance(quid, tuple):
                        quid = quid[0]
                    stab.add_qubits([quid])

                    pos = qu.GetCenter()
                    lineSource = vtk.vtkLineSource()
                    lineSource.SetPoint1(stab_g.GetCenter())
                    lineSource.SetPoint2(pos)
                    # Setup actor and mapperg
                    lineMapper = vtk.vtkPolyDataMapper()
                    lineMapper.SetInputConnection(lineSource.GetOutputPort())
                    lineActor = StabilizerActor(stab_g.get_parent_id())
                    lineActor.SetMapper(lineMapper)
                    lineActor.GetProperty().SetColor(colors.GetColor3d('Red'))
                    lineActor.GetProperty().SetLineWidth(4.55)
                    self.GetDefaultRenderer().AddActor(lineActor)

                    stab.add_actors([lineActor])

                    for actor in old_qubit_actors:
                        actor.GetProperty().SetDiffuseColor(colors.GetColor3d('White'))
                        actor.GetProperty().SetDiffuse(.8)
                        actor.GetProperty().SetSpecular(.5)
                        actor.GetProperty().SetSpecularColor(colors.GetColor3d('White'))
                        actor.GetProperty().SetAmbientColor(colors.GetColor3d('White'))
                        actor.GetProperty().EdgeVisibilityOff()
                    old_qubit_actors.clear()
                selected_qubits.clear()
                selected_stabilizers.clear()
                
            
        
        

            
        return self.OnKeyPress()
            


def main():
    # A renderer and render window
    
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(BACKGROUND_COLOR)
    
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
    for i in range(NUM_GRID):
        for j in range(NUM_GRID):
    
            qubit = tanner.create_qubit()
            
            source = vtk.vtkSphereSource()
            
            # random position and radius
            
            source.SetRadius(Qubit.radius)
            source.SetCenter(i, j, 0)
            source.SetPhiResolution(11)
            source.SetThetaResolution(21)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())
            actor = QubitActor(qubit.get_id())
            actor.SetMapper(mapper)
            
            
            actor.GetProperty().SetDiffuseColor(colors.GetColor3d('White'))
            actor.GetProperty().SetDiffuse(.8)
            actor.GetProperty().SetSpecular(.5)
            actor.GetProperty().SetSpecularColor(colors.GetColor3d('White'))
            actor.GetProperty().SetAmbientColor(colors.GetColor3d('White'))

            actor.GetProperty().SetSpecularPower(30.0)
            
            
            qubit.add_actors([actor])
            
            renderer.AddActor(actor)
    
    # Start
    interactor.Initialize()
    renwin.Render()
    interactor.Start()


if __name__ == '__main__':
    main()