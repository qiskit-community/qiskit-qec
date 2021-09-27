import vtk
import numpy

# Example: Create a lattice of 3D spheres then allow mouse selection and de-selection. 
# Change color of sphere that is selected. Selected spheres (actors) are stored in the 
# gloabl list selected

# Cheap global list to store references of actors selected
# This would obviously be done via some class etc if done for real. 
# This works for the demo

# press key 'w' for wireframe mode
# press key 's' for solid mode

selected = list()

# Create a lattice of spheres in 3D interactor style with ability to select a sphere and change its color

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


colors = vtk.vtkNamedColors()
class MouseInteractorHighLightActor(vtk.vtkInteractorStyleTrackballCamera):

    picked = dict()

    def __init__(self, parent=None):
        self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)

    def leftButtonPressEvent(self, obj, event):
        clickPos = self.GetInteractor().GetEventPosition()

        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())

        # get the new
        self.NewPickedActor = picker.GetActor()

        # If something was selected
        if self.NewPickedActor:

            if self.NewPickedActor in selected:
                self.NewPickedActor.GetProperty().SetColor(colors.GetColor3d('Orange'))
                self.NewPickedActor.GetProperty().SetDiffuse(.8)
                self.NewPickedActor.GetProperty().SetSpecular(.5)
                selected.remove(self.NewPickedActor)
            else:
                # Highlight the picked actor by changing its properties
                self.NewPickedActor.GetProperty().SetColor(colors.GetColor3d('Red'))
                self.NewPickedActor.GetProperty().SetDiffuse(1.0)
                self.NewPickedActor.GetProperty().SetSpecular(0.0)
                #self.NewPickedActor.GetProperty().EdgeVisibilityOn()

                selected.append(self.NewPickedActor)
            

        self.OnLeftButtonDown()
        return


def main():
    # A renderer and render window
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d('DimGray'))

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

    # --------- Define lattice ----------
    #lat_vec1 = numpy.array([0.0, 1.0])
    #lat_vec2 = numpy.array([ 1.0, 0.0])
    lat_vec1 = numpy.array([-40., -10.0])
    lat_vec2 = numpy.array([ 18.0, -37.0])
    lattice_generators = [lat_vec1, lat_vec2]

    mask = (400, 400)
    raw_lattice = generate_lattice(mask, lattice_generators)

    # Convert to integer lists
    # Note: By casting to int the lattice is not be regular for many settings
    # Need to fix this.
    raw_lattice = [list(item.astype(int)) for item in raw_lattice]

    # 
    lattice = vtk.vtkPoints()
    for lat_point in raw_lattice:
        lattice.InsertNextPoint(*lat_point,0)


    actor = points_to_glyph(lattice, 0.02)

    # assign actor to the renderer
    renderer.AddActor(actor)

    # Start
    interactor.Initialize()
    renwin.Render()
    interactor.Start()


def points_to_glyph(points, scale):
    """
    Convert points to glyphs.
    :param points: The points to glyph.
    :param scale: The scale, used to determine the size of the
                  glyph representing the point, expressed as a
                  fraction of the largest side of the bounding
                  box surrounding the points. e.g. 0.05
    :return: The actor.
    """

    bounds = points.GetBounds()
    max_len = 0.0
    for i in range(0, 3):
        max_len = max(bounds[i + 1] - bounds[i], max_len)

    sphere_source = vtk.vtkSphereSource()
    sphere_source.SetRadius(scale * max_len)
    sphere_source.SetPhiResolution(11)
    sphere_source.SetThetaResolution(21)

    pd = vtk.vtkPolyData()
    pd.SetPoints(points)

    mapper = vtk.vtkGlyph3DMapper()
    mapper.SetInputData(pd)
    mapper.SetSourceConnection(sphere_source.GetOutputPort())
    mapper.ScalarVisibilityOff()
    mapper.ScalingOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetDiffuseColor(0.5, 0.1, 0.9)
    actor.GetProperty().SetDiffuse(.8)
    actor.GetProperty().SetSpecular(.5)
    actor.GetProperty().SetSpecularColor(colors.GetColor3d('White'))
    actor.GetProperty().SetSpecularPower(30.0)
    actor.GetProperty().SetColor(colors.GetColor3d('Red'))

    return actor


if __name__ == '__main__':
    main()
