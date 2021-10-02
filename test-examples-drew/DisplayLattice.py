import vtk
import numpy


# Example: Create a lattice of points and display
# set3D=True  : 3D vtkPolyDataMapper and vtkActor : lattice centers in window
# set3D=False : 2D vtkPolyDataMapper2D and vtkActor2D : lattice fixed in position
# set2Dstyle = True : Use 2D Image Style Interactor (no rotation in 3D)

set2Dstyle = False

set3D = True

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
    numpy.array([ 50, -70.0])]
    image_shape = (500, 500)
    spots = generate_lattice(image_shape, lattice_vectors)
    spots = [list(item.astype(int)) for item in spots]
    for lattice_point in spots:
        points.InsertNextPoint(*lattice_point,0)
    # Note: Casting to int will miss up the lattice. Need to fix this issue
    # Scale to corrcet

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)

    # Create anything you want here, we will use a polygon for the demo.
    polygonSource = vtk.vtkRegularPolygonSource()  # default is 6 sides
    polygonSource.SetNumberOfSides(50)
    polygonSource.SetRadius(4)

    glyph2D = vtk.vtkGlyph2D()
    glyph2D.SetSourceConnection(polygonSource.GetOutputPort())
    glyph2D.SetInputData(polydata)
    glyph2D.Update()

    if set3D:
        mapper = vtk.vtkPolyDataMapper()
    else:
        mapper = vtk.vtkPolyDataMapper2D()
    mapper.SetInputConnection(glyph2D.GetOutputPort())
    mapper.Update()

    if set3D:
        actor = vtk.vtkActor()
    else:
        actor = vtk.vtkActor2D()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d('Orange'))

    # Visualize
    renderer = vtk.vtkRenderer()
    renderer.UseFXAAOn() # Enable anti-aliasing 

    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('DimGray'))


    if set2Dstyle:
        style = vtk.vtkInteractorStyleImage()
        renderWindowInteractor.SetInteractorStyle(style)

    renderWindow.SetWindowName('Glyph2D');
    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__ == '__main__':
    main()
