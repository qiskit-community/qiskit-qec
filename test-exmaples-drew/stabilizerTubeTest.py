#!/usr/bin/env python

import vtk
import numpy
import math


def generate_lattice(mask, lattice_vectors, center_pix='image', edge_buffer=2):
    '''
    Lattice generator from https://stackoverflow.com/questions/6141955/efficiently-generate-a-lattice-of-points-in-python
    by kiyo
    TODO: This should be rewritten to avoid problems with highly skewed lattices
    '''
    ##Preprocessing. Not much of a bottleneck:
    if center_pix == 'image':
        center_pix = numpy.array(mask) // 2
    else: ##Express the center pixel in terms of the lattice vectors
        center_pix = numpy.array(center_pix) - (numpy.array(mask) // 2)
        lattice_components = numpy.linalg.solve(
            numpy.vstack(lattice_vectors[:2]).T,
            center_pix)
        lattice_components -= lattice_components // 1
        center_pix = (lattice_vectors[0] * lattice_components[0] +
                      lattice_vectors[1] * lattice_components[1] +
                      numpy.array(mask)//2)
    num_vectors = int( ##Estimate how many lattice points we need
        max(mask) / numpy.sqrt(lattice_vectors[0]**2).sum())
    lattice_points = []
    lower_bounds = numpy.array((edge_buffer, edge_buffer))
    upper_bounds = numpy.array(mask) - edge_buffer

    ##SLOW LOOP HERE. 'num_vectors' is often quite large.
    for i in range(-num_vectors, num_vectors):
        for j in range(-num_vectors, num_vectors):
            lp = i * lattice_vectors[0] + j * lattice_vectors[1] + center_pix
            if all(lower_bounds < lp) and all(lp < upper_bounds):
                lattice_points.append(lp)
    return lattice_points

def main():
    colors = vtk.vtkNamedColors()


    # --------- Define lattice ----------
    lat_vec1 = numpy.array([0.0, 1.0])
    lat_vec2 = numpy.array([ 1.0, 0.0])
    lattice_generators = [lat_vec1, lat_vec2]

    mask = (20, 20)
    raw_lattice = generate_lattice(mask, lattice_generators)

    # Convert to integer lists
    # Note: By casting to int the lattice is not be regular for many settings
    # Need to fix this.
    raw_lattice = [list(item.astype(int)) for item in raw_lattice]

    # 
    lattice = vtk.vtkPoints()
    for lat_point in raw_lattice:
        lattice.InsertNextPoint(*lat_point,0)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(lattice)

    # ---------- Create a sphere source ------------
    sphereSource = vtk.vtkSphereSource()
    sphereSource.SetRadius(0.1)
    sphereSource.SetPhiResolution(20)
    sphereSource.SetThetaResolution(20)

    # Display Spheres at the vertices of the structured grid


    # Coordinates for a square stabilizer
    point_list = [(0,0),(1,0),(1,1),(0,1)]
    labels = [0,1,2,3,1] # Last elements a salar representing color of cell

    # Convert coordinates in 3D vtkPoints
    points = vtk.vtkPoints()
    for point in point_list:
        points.InsertNextPoint(*point,0)

    
    # Create an UnstructuredGrid
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.Allocate(100)

    ugrid.InsertNextCell(vtk.VTK_QUAD,4,labels[:4])

    ugrid.SetPoints(points)

    # Setup color lut and use scalar data to color cell

    color_array = vtk.vtkFloatArray()
    i = color_array.InsertNextTuple((labels[4],))

    ugrid.GetCellData().SetScalars(color_array)

    # Create the color LUT

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(5)
    lut.SetTableRange(0,4)
    lut.Build()
    lut.SetTableValue(0, 0, 0, 0, 1) # Black
    lut.SetTableValue(1, 1, 0, 0, 1) # Red
    lut.SetTableValue(2, 0, 1, 0, 1) # Green
    lut.SetTableValue(3, 0, 0, 1, 1) # Blue
    lut.SetTableValue(4, 1, 1, 1, 1) # White 

    # Create a mapper to map the color

    ugridMapper = vtk.vtkDataSetMapper()
    ugridMapper.SetScalarModeToUseCellData()
    ugridMapper.UseLookupTableScalarRangeOn()
    ugridMapper.SetLookupTable(lut)
    ugridMapper.SetInputData(ugrid)

    # Get the edges

    # Extract the edges of the triangles just found.
    extractEdges = vtk.vtkExtractEdges()
    extractEdges.SetInputData(ugrid)   

    # Create the tube filter

    tubeFilter = vtk.vtkTubeFilter()
    tubeFilter.SetInputConnection(extractEdges.GetOutputPort())
    tubeFilter.SetRadius(0.025)
    tubeFilter.SetNumberOfSides(50)
    tubeFilter.Update()

    # Create the edge mapper

    edgeMapper = vtk.vtkPolyDataMapper()
    edgeMapper.SetInputConnection(tubeFilter.GetOutputPort())
    edgeMapper.SetScalarRange(0, 26)

    # Create edge actor

    edgeActor = vtk.vtkActor()
    edgeActor.SetMapper(edgeMapper)
    edgeActor.GetProperty().SetSpecular(0.6)
    edgeActor.GetProperty().SetSpecularPower(30)

    # Create a renderer, render window, and interactor
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(edgeActor)
    renderer.SetBackground(colors.GetColor3d('SlateGray'))

    aCamera = vtk.vtkCamera()
    aCamera.Azimuth(-40.0)
    aCamera.Elevation(50.0)

    renderer.SetActiveCamera(aCamera)
    renderer.ResetCamera()

    renderWindow.SetSize(640, 640)
    renderWindow.SetWindowName('Stabilizer')
    renderWindow.Render()

    renderWindowInteractor.Start()



if __name__ == '__main__':
    main()
