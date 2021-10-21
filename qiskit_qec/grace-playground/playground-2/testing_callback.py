import vtk


def main():
    colors = vtk. vtkNamedColors ()
    cone = vtk. vtkConeSource ()
    cone. SetHeight ( 3.0 )
    cone. SetRadius ( 1.0 )
    cone. SetResolution ( 10 )
    coneMapper = vtk. vtkPolyDataMapper ()
    coneMapper . SetInputConnection (cone. GetOutputPort ())
    coneActor = vtk.vtkActor ()
    coneActor . SetMapper ( coneMapper )
    ren1 = vtk. vtkRenderer ()
    ren1.AddActor( coneActor )
    ren1. SetBackground (colors. GetColor3d ('MidnightBlue'))
    renWin = vtk. vtkRenderWindow ()
    renWin. AddRenderer (ren1)
    renWin.SetSize (300 , 300)
    iren = vtk. vtkRenderWindowInteractor ()
    iren. SetRenderWindow (renWin)
    style = vtk. vtkInteractorStyleTrackballCamera ()
    iren. SetInteractorStyle (style)
    iren. Initialize ()
    iren.Start ()
    
if __name__ == '__main__':
    main()
    