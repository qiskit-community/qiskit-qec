# selectedQubits = [] #colour
#


from __future__ import print_function

import vtk


def runme():
    sphere = vtk.vtkSphereSource()
    sphereMapper = vtk.vtkPolyDataMapper()
    sphereMapper.SetInputConnection(sphere.GetOutputPort())
    sphereActor = vtk.vtkLODActor()
    sphereActor.SetMapper(sphereMapper)
    cone = vtk.vtkConeSource()
    glyph = vtk.vtkGlyph3D()
    glyph.SetInputConnection(sphere.GetOutputPort())
    glyph.SetSourceConnection(cone.GetOutputPort())
    
    glyph.SetVectorModeToUseNormal()
    glyph.SetScaleModeToScaleByVector()
    glyph.SetScaleFactor(0.25)
    spikeMapper = vtk.vtkPolyDataMapper()
    spikeMapper.SetInputConnection(glyph.GetOutputPort())
    spikeActor = vtk.vtkLODActor()
    spikeActor.SetMapper(spikeMapper)
    
    textMapper = vtk.vtkTextMapper()
    tprop = textMapper.GetTextProperty()
    tprop.SetFontFamilyToArial()
    tprop.SetFontSize(10)
    tprop.BoldOn()
    tprop.ShadowOn()
    tprop.SetColor(1, 0, 0)
    textActor = vtk.vtkActor2D()
    textActor.VisibilityOff()
    textActor.SetMapper(textMapper)
    picker = vtk.vtkCellPicker()
    
    
    def annotatePick(obj, event, picker, curtextActor, curtextMapper):
        print(f"picking my {obj} \n\nw event {event}")
        if picker.GetCellId() < 0:
            curtextActor.VisibilityOff()
        else:
            selPt = picker.GetSelectionPoint()
            pickPos = picker.GetPickPosition()
            curtextMapper.SetInput("(%.6f, %.6f, %.6f)" % pickPos)
            curtextActor.SetPosition(selPt[ :2 ])
            curtextActor.VisibilityOn()
    
    
    picker.AddObserver("EndPickEvent",
        lambda obj, event, p=picker, ta=textActor, tm=textMapper: annotatePick(obj, event, p, ta, tm))
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    iren.SetPicker(picker)
    ren.AddActor2D(textActor)
    ren.AddActor(sphereActor)
    ren.AddActor(spikeActor)
    ren.SetBackground(1, 1, 1)
    renWin.SetSize(300, 300)
    ren.ResetCamera()
    cam1 = ren.GetActiveCamera()
    cam1.Zoom(1.4)
    iren.Initialize()
    picker.Pick(85, 126, 0, ren)
    renWin.Render()
    iren.Start()


if __name__ == "__main__":
    runme()
