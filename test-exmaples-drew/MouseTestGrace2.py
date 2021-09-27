#!/usr/bin/env python
import vtk


def main():
    colors = vtk.vtkNamedColors()

    source = vtk.vtkSphereSource()
    source.SetCenter(0, 0, 0)
    source.SetRadius(1)
    source.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(source.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d('MistyRose'))

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d('SlateGray'))
    renderer.AddActor(actor)

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(renderer)
    renwin.SetWindowName('MouseEventsObserver')

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
    interactor.SetRenderWindow(renwin)

    def DummyFunc1(obj, ev):
        print('Before Event')

    def DummyFunc2(obj, ev):
        print('After Event')
        print(f'obj is \n{obj}\n event is {ev}')


    # Print interator gives you a list of registered observers of the current
    # interactor style
    # print(interactor)

    # adding priorities allow to control the order of observer execution
    # (highest value first! if equal the first added observer is called first)
    interactor.AddObserver('LeftButtonPressEvent', DummyFunc1, 1.0)
    interactor.AddObserver('LeftButtonPressEvent', DummyFunc2, -1.0)
    interactor.Initialize()
    renwin.Render()
    interactor.Start()


if __name__ == '__main__':
    main()
