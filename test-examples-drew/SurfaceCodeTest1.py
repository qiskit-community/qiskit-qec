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
    color_names = vtk.vtkNamedColors()

    points = vtk.vtkPoints()
    
    lattice_vectors = [
    numpy.array([0.0, 1.0]),
    numpy.array([ 1.0, 0.0])]
    image_shape = (20, 20)
    spots = generate_lattice(image_shape, lattice_vectors)
    spots = [list(item.astype(int)) for item in spots]
    for lattice_point in spots:
        points.InsertNextPoint(*lattice_point,0)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)

    '''
    # Create anything you want here, we will use a polygon for the demo.
    polygonSource = vtk.vtkRegularPolygonSource()  # default is 6 sides
    polygonSource.SetRadius(0.1)
    polygonSource.SetNumberOfSides(50)

    glyph2D = vtk.vtkGlyph2D()
    glyph2D.SetSourceConnection(polygonSource.GetOutputPort())
    glyph2D.SetInputData(polydata)
    glyph2D.Update()

    gmapper = vtk.vtkPolyDataMapper()
    gmapper.SetInputConnection(glyph2D.GetOutputPort())
    gmapper.Update()

    gactor = vtk.vtkActor()
    gactor.SetMapper(gmapper)
    gactor.GetProperty().SetColor(color_names.GetColor3d('Gold'))
    '''


    sphereSource = vtk.vtkSphereSource()
    sphereSource.SetRadius(0.1)
    sphereSource.SetPhiResolution(11)
    sphereSource.SetThetaResolution(21)
    sphereSource.Update()

    glyph3d = vtk.vtkGlyph3D()
    glyph3d.GeneratePointIdsOn()
    glyph3d.SetSourceConnection(sphereSource.GetOutputPort())
    glyph3d.SetInputData(polydata)
    glyph3d.SetScaleModeToDataScalingOff()
    glyph3d.Update()
    
    
    gmapper = vtk.vtkPolyDataMapper()
    gmapper.SetInputConnection(glyph3d.GetOutputPort())

    gactor = vtk.vtkActor()
    gactor.SetMapper(gmapper)
    gactor.GetProperty().SetColor(color_names.GetColor3d('Gold'))


    # Setup four points
    rpoints = vtk.vtkPoints()

    N = 100
    R=2
    angle = math.pi/N

    for k in range(N+1):
        nangle = k*angle
        rpoints.InsertNextPoint(R*math.cos(nangle),R*math.sin(nangle),0)

    # Create the polygon
    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(N+1)

    for k in range(N+1):
         polygon.GetPointIds().SetId(k, k)

    # Add the polygon to a list of polygons
    polygons = vtk.vtkCellArray()
    polygons.InsertNextCell(polygon)

    # Create a PolyData
    polygonPolyData = vtk.vtkPolyData()
    polygonPolyData.SetPoints(rpoints)
    polygonPolyData.SetPolys(polygons)

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polygonPolyData)

    actor_sd = vtk.vtkActor()
    actor_sd.SetMapper(mapper)
    actor_sd.GetProperty().SetColor(color_names.GetColor3d('Green'))


    x = [[0,0,0],[0,1,0],[0,2,0],[1,0,0],[1,1,0],[1,2,0],[2,0,0],[2,1,0],[2,2,0]]

    pts = [[0,3,4,1,2],[4,7,8,5,2],[1,4,5,2,1],[3,6,7,4,1]]

    color_array = vtk.vtkFloatArray()
    for item in pts:
        i = color_array.InsertNextTuple((item[4],))
    #color_array.InsertNextTuple((2,))

    points = vtk.vtkPoints()
    for i in range(0, len(x)):
        points.InsertPoint(i, x[i])


    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.Allocate(100)

    ugrid.InsertNextCell(vtk.VTK_QUAD, 4, pts[0][:4])
    #ugrid.InsertNextCell(vtk.VTK_QUAD, 4, pts[1][:4])
    ugrid.InsertNextCell(vtk.VTK_QUAD, 4, pts[2][:4])
    ugrid.InsertNextCell(vtk.VTK_QUAD, 4, pts[3][:4])

    ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [0,2,5])

    ugrid.SetPoints(points)

    ugrid.GetCellData().SetScalars(color_array)

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(5)
    lut.SetTableRange(0,4)
    lut.Build()
    lut.SetTableValue(0, 0, 0, 0, 1) # Black
    lut.SetTableValue(1, 1, 0, 0, 1) # Red
    lut.SetTableValue(2, 0, 1, 0, 1) # Green
    lut.SetTableValue(3, 0, 0, 1, 1) # Blue
    lut.SetTableValue(4, 1, 1, 1, 1) # White 


    ugridMapper = vtk.vtkDataSetMapper()

    ugridMapper.SetScalarModeToUseCellData()
    ugridMapper.UseLookupTableScalarRangeOn()
    ugridMapper.SetLookupTable(lut)

    ugridMapper.SetInputData(ugrid)



    ugridActor = vtk.vtkActor()
    ugridActor.SetMapper(ugridMapper)
    ugridActor.GetProperty().SetColor(color_names.GetColor3d('Black'))
    ugridActor.GetProperty().EdgeVisibilityOn()


    # Visualize
    renderer = vtk.vtkRenderer()
    renderer.UseFXAAOn()

    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName('Polygon')
    renderWindow.AddRenderer(renderer)

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(style)

    cam_orient_manipulator = vtk.vtkCameraOrientationWidget()
    cam_orient_manipulator.SetParentRenderer(renderer)
    # Enable the widget.
    cam_orient_manipulator.On()

    renderer.AddActor(actor_sd)
    renderer.AddActor(ugridActor)
    renderer.AddActor(gactor)

    trans0 = vtk.vtkTransform()
    trans1 = vtk.vtkTransform()

    actor_sd.SetUserTransform(trans0)
    ugridActor.SetUserTransform(trans1)

    trans0.Translate(7.5, 14.0, 0.0)
    trans0.Scale(0.25,0.25,0.0)
    trans1.Translate(7.0,12.0,0.0)

    renderer.SetBackground(color_names.GetColor3d('Black'))

    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__ == '__main__':
    main()