#!/usr/bin/env python

"""
This code is based on the VTK file: Examples/Modelling/Tcl/DelMesh.py.

This example demonstrates how to use 2D Delaunay triangulation.
We create a fancy image of a 2D Delaunay triangulation. Points are
 randomly generated.
"""

import vtk
import math
import numpy


def generate_lattice( #Help me speed up this function, please!
    image_shape, lattice_vectors, center_pix='image', edge_buffer=2):

    ##Preprocessing. Not much of a bottleneck:
    if center_pix == 'image':
        center_pix = numpy.array(image_shape) // 2
    else: ##Express the center pixel in terms of the lattice vectors
        center_pix = numpy.array(center_pix) - (numpy.array(image_shape) // 2)
        lattice_components = numpy.linalg.solve(
            numpy.vstack(lattice_vectors[:2]).T,
            center_pix)
        lattice_components -= lattice_components // 1
        center_pix = (lattice_vectors[0] * lattice_components[0] +
                      lattice_vectors[1] * lattice_components[1] +
                      numpy.array(image_shape)//2)
    num_vectors = int( ##Estimate how many lattice points we need
        max(image_shape) / numpy.sqrt(lattice_vectors[0]**2).sum())
    lattice_points = []
    lower_bounds = numpy.array((edge_buffer, edge_buffer))
    upper_bounds = numpy.array(image_shape) - edge_buffer

    ##SLOW LOOP HERE. 'num_vectors' is often quite large.
    for i in range(-num_vectors, num_vectors):
        for j in range(-num_vectors, num_vectors):
            lp = i * lattice_vectors[0] + j * lattice_vectors[1] + center_pix
            if all(lower_bounds < lp) and all(lp < upper_bounds):
                lattice_points.append(lp)
    return lattice_points


def main():
    colors = vtk.vtkNamedColors()

    points = vtk.vtkPoints()

    lattice_vectors = [
    numpy.array([0.0, 50.0]),
    numpy.array([ 50.0, 0.0])]
    image_shape = (400, 400)
    spots = generate_lattice(image_shape, lattice_vectors)
    spots = [list(item.astype(int)) for item in spots]
    for lattice_point in spots:
        points.InsertNextPoint(*lattice_point,0)

    # Create a polydata with the points we just created.
    profile = vtk.vtkPolyData()
    profile.SetPoints(points)

    

    '''
    # Perform a 2D Delaunay triangulation on them.
    delny = vtk.vtkDelaunay2D()
    delny.SetInputData(profile)
    delny.SetTolerance(0.001)
    mapMesh = vtk.vtkPolyDataMapper()
    mapMesh.SetInputConnection(delny.GetOutputPort())
    meshActor = vtk.vtkActor()
    meshActor.SetMapper(mapMesh)
    meshActor.GetProperty().SetColor(colors.GetColor3d('MidnightBlue'))
    '''

    # We will now create a nice looking mesh by wrapping the edges in tubes,
    # and putting fat spheres at the points.
    extract = vtk.vtkExtractEdges()
    extract.SetInputConnection(delny.GetOutputPort())
    tubes = vtk.vtkTubeFilter()
    tubes.SetInputConnection(extract.GetOutputPort())
    tubes.SetRadius(0.01)
    tubes.SetNumberOfSides(20)
    mapEdges = vtk.vtkPolyDataMapper()
    mapEdges.SetInputConnection(tubes.GetOutputPort())
    edgeActor = vtk.vtkActor()
    edgeActor.SetMapper(mapEdges)
    edgeActor.GetProperty().SetColor(colors.GetColor3d('peacock'))
    edgeActor.GetProperty().SetSpecularColor(1, 1, 1)
    edgeActor.GetProperty().SetSpecular(0.3)
    edgeActor.GetProperty().SetSpecularPower(20)
    edgeActor.GetProperty().SetAmbient(0.2)
    edgeActor.GetProperty().SetDiffuse(0.8)

    ball = vtk.vtkSphereSource()
    ball.SetRadius(0.025)
    ball.SetThetaResolution(12)
    ball.SetPhiResolution(12)
    balls = vtk.vtkGlyph3D()
    balls.SetInputConnection(delny.GetOutputPort())
    balls.SetSourceConnection(ball.GetOutputPort())
    mapBalls = vtk.vtkPolyDataMapper()
    mapBalls.SetInputConnection(balls.GetOutputPort())
    ballActor = vtk.vtkActor()
    ballActor.SetMapper(mapBalls)
    ballActor.GetProperty().SetColor(colors.GetColor3d('hot_pink'))
    ballActor.GetProperty().SetSpecularColor(1, 1, 1)
    ballActor.GetProperty().SetSpecular(0.3)
    ballActor.GetProperty().SetSpecularPower(20)
    ballActor.GetProperty().SetAmbient(0.2)
    ballActor.GetProperty().SetDiffuse(0.8)

    # Create the rendering window, renderer, and interactive renderer.
    ren = vtk.vtkRenderer()
    ren.UseFXAAOn()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # Add the actors to the renderer, set the background and size.
    ren.AddActor(ballActor)
    ren.AddActor(edgeActor)
    ren.SetBackground(colors.GetColor3d('AliceBlue'))
    renWin.SetSize(512, 512)
    renWin.SetWindowName('DelaunayMesh')

    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1.3)

    # Interact with the data.
    iren.Initialize()
    renWin.Render()
    iren.Start()


if __name__ == '__main__':
    main()