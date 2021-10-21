# because vtk is for visualizing...*cry*

#
# Let me start with a small clarification. Actors in VTK are simple container objects that combine geometric data and visualization properties for rendering. VTK does not offer much functionality on the level of actors. If you want to detect colliding objects, you need to solve this for the geometries (vtkPolyData).
#
# There is no generic collision-detection- or body-intersection-engine in VTK. What you can do in VTK is to apply boolean operations. However, this is not easy to achieve robustly. There are two main approaches:
#
# Method A: Boolean operations on the mesh
# Use [vtkBooleanOperationPolyDataFilter][1] to apply boolean operations directly on the mesh. See here for an example. Unfortunately, this will fail in your case because of the mesh properties of your surfaces (hit key W when looking at the surface in the RenderWindow to inspect the wireframe of your mesh). vtkBooleanOperationPolyDataFilter will perform best if the mesh's triangles are small and have decent condition numbers, that is, if the triangles are not too spiky. (See the manual of the Verdict toolbox for some triangle metrics.) The mesh of your extruded disks, however, consists of very long, spiky triangles. What you would need to do is to first remesh the surface. VTK does not offer remeshing facilities out of the box, but related toolboxes such as VMTK do.
#
# Getting approach (A) right for generic geometries is tricky, because you may end up with non-manifold or leaky surfaces. Also, vtkBooleanOperationPolyDataFilter is known to have some bugs (see here or here). Let's hope that these problems will be fixed some day.
#
# Method B: Boolean operations on implicit functions
# A second approach is to work with implicit functions. For instance, you can represent your tubes as implicit cylinders, and intersect these using [vtkImplicitBoolean][7]. See here for an example. The problem with this approach is that you need to convert the implicit representation of the objects to yield the resulting surface mesh. The marching cubes algorithm in VTK is rather slow, so you need to wait a long time for high resolutions. And it is impossible to preserve sharp edges. On the other hand, it is more robust and easier to deal with.
#
# Sample code
# The code below demonstrates both cases. The remeshing feature I cannot share with you here, therefore only the implicit boolean is functional. The screenshot shows, how the result looks like. (Yellow: the input surfaces, Red: the result)
#
# enter image description here
#
# For more details on terminology and alternative problem formulations, refer to "Collision detection between geometric models: a survey" by Ming and Gottschalk, 1999.
#
# Boolean operation on implicit functions
# This code has been written by normanius under the CC BY-SA 4.0 license.
# License: https://creativecommons.org/licenses/by-sa/4.0/
# Author:  https://stackoverflow.com/users/3388962/normanius
# Date:    July 2018

import vtk
import numpy as np

def compute_transform(start, end):
    # Better compute the matrix in numpy!

    normalized_x = [0]*3
    normalized_y = [0]*3
    normalized_z = [0]*3

    # The X axis is a vector from start to end
    vtk.vtkMath.Subtract(end, start, normalized_x)
    # length = vtk.vtkMath.Norm(normalized_x)
    vtk.vtkMath.Normalize(normalized_x)
    
    # The Z axis is an arbitrary vector cross X
    rng = vtk.vtkMinimalStandardRandomSequence()
    rng.SetSeed(8775070)  # For testing.
    arbitrary = [0]*3
    for i in range(0, 3):
        rng.Next()
        arbitrary[i] = rng.GetRangeValue(-10, 10)
    vtk.vtkMath.Cross(normalized_x, arbitrary, normalized_z)
    vtk.vtkMath.Normalize(normalized_z)
    
    # The Y axis is Z cross X
    vtk.vtkMath.Cross(normalized_z, normalized_x, normalized_y)
    matrix = vtk.vtkMatrix4x4()
    # Create the direction cosine matrix
    
    matrix.Identity()
    # print(f"Creating the idenity: {matrix}")
    for i in range(3):
        matrix.SetElement(i, 0, normalized_x[i])
        matrix.SetElement(i, 1, normalized_y[i] * -1)
        matrix.SetElement(i, 2, normalized_z[i] * -1)

   # print(f"x: {normalized_x}, and y: {normalized_y}, and zzz: {normalized_z}\n and amtrix: {matrix}")
    transform = vtk.vtkTransform()
    transform.Translate(start)          # translate to starting point
    transform.Concatenate(matrix)       # apply direction cosines
    transform.RotateY(90.0)             # align cylinder
    # Don't scale! This changes mesh properties (e.g. aspect ratio)
    #transform.Scale(1.0, 1.0, length)   # scale along the height vector
    return transform

def transform_item(item, transform):
    transformed = vtk.vtkTransformPolyDataFilter()
    transformed.SetInputConnection(item.GetOutputPort())
    transformed.SetTransform(transform)
    transformed.Update()
    return transformed

def create_pipe(radius, thickness, height):
    # This type of pipe is not suited for remeshing, because remeshing does not
    # preserve (feature-) edges. See create_pipe2
    assert(radius>thickness)
    disk = vtk.vtkDiskSource()
    disk.SetCircumferentialResolution(128)
    disk.SetRadialResolution(1)
    disk.SetOuterRadius(radius)
    disk.SetInnerRadius(radius - thickness)
    pipe = vtk.vtkLinearExtrusionFilter()
    pipe.SetInputConnection(disk.GetOutputPort())
    pipe.SetExtrusionTypeToNormalExtrusion()
    pipe.SetVector(0, 0, 1)
    pipe.SetScaleFactor(height)
    pipe.Update()
    return pipe

def create_pipe_implicit(radius, thickness, height):
    center = np.array([0,0,0])
    axis = np.array([0,0,1])
    centerTop = center + height*axis
    centerBottom = center

    # Outer cylinder.
    outer = vtk.vtkCylinder()
    outer.SetCenter(center)
    outer.SetAxis(axis)
    outer.SetRadius(radius)
    # Inner cylinder.
    inner = vtk.vtkCylinder()
    inner.SetCenter(center)
    inner.SetAxis(axis)
    inner.SetRadius(radius-thickness)
    # Top face.
    plane1 = vtk.vtkPlane()
    plane1.SetOrigin(centerTop)
    plane1.SetNormal(np.array(outer.GetAxis()))
    # Bottom face.
    plane2 = vtk.vtkPlane()
    plane2.SetOrigin(centerBottom)
    plane2.SetNormal(-np.array(outer.GetAxis()))
    # Put things together.
    difference = vtk.vtkImplicitBoolean()
    difference.AddFunction(outer)
    difference.AddFunction(inner)
    difference.SetOperationTypeToDifference()
    intersection = vtk.vtkImplicitBoolean()
    intersection.AddFunction(difference)
    intersection.AddFunction(plane1)
    intersection.AddFunction(plane2)
    intersection.SetOperationTypeToIntersection()
    pipe = intersection
    # Also return inner and outer cylinder.
    intersection = vtk.vtkImplicitBoolean()
    intersection.AddFunction(inner)
    intersection.AddFunction(plane1)
    intersection.AddFunction(plane2)
    intersection.SetOperationTypeToIntersection()
    inner = intersection
    intersection = vtk.vtkImplicitBoolean()
    intersection.AddFunction(outer)
    intersection.AddFunction(plane1)
    intersection.AddFunction(plane2)
    intersection.SetOperationTypeToIntersection()
    outer = intersection
    return pipe, inner, outer

def add_to_renderer(renderer, item, color, opacity=1., translate=None):
    colors = vtk.vtkNamedColors()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetScalarVisibility(False)
    mapper.SetInputConnection(item.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d(color))
    actor.GetProperty().SetOpacity(opacity)
    if translate:
        trafo = vtk.vtkTransform()
        trafo.Translate(translate)
        actor.SetUserTransform(trafo)
    renderer.AddActor(actor)
    return mapper, actor

def evaluate_implicit(implicit_function, resolution, bounds):
    sampled = vtk.vtkSampleFunction()
    sampled.SetSampleDimensions(resolution, resolution, resolution)
    sampled.SetModelBounds(bounds)
    sampled.SetImplicitFunction(implicit_function)
    iso = vtk.vtkMarchingCubes()
    iso.SetValue(0,0.)
    iso.SetInputConnection(sampled.GetOutputPort())
    iso.Update()
    return iso

def main():
    colors = vtk.vtkNamedColors()

    # Params.
    radius = 2.
    thickness = 0.5
    start_point = np.array([0] * 3)
    end_point = np.array([0, 0, 10])
    length = np.linalg.norm(start_point-end_point)

    radius2 = 2.
    thickness2 = 0.5
    start_point2 = np.array([-10, 0, 0])
    end_point2 = np.array([0, 0, 7])
    length2 = np.linalg.norm(start_point2-end_point2)

    # Compute transforms.
    transform = compute_transform(start_point, end_point)
    transform2 = compute_transform(start_point2, end_point2)

    ############################################################################
    # BOOLEAN OPERATIONS ON MESHES
    ############################################################################
    if False:
        pipe, inner, outer = create_pipe2(radius=radius, thickness=thickness, height=length)
        pipe2, inner2, outer2 = create_pipe2(radius=radius2, thickness=thickness2, height=length2)
        # Apply the transforms.
        pipe = transform_item(pipe, transform)
        inner = transform_item(inner, transform)
        outer = transform_item(outer, transform)
        pipe2 = transform_item(pipe2, transform2)
        inner2 = transform_item(inner2, transform2)
        outer2 = transform_item(outer2, transform2)

        #pipe_2m1 = boolean_combine(pipe2, pipe, 'difference')
        pipe_2m1 = boolean_combine(pipe2, pipe, 'union') # Ugly! There is a bug in vtk!
        result_bool = pipe_2m1
        #result_bool = boolean_combine(pipe, pipe_2m1, 'union')
        #result_bool = remeshSurface(result_bool, targetArea=.1, iterations=10)

        # Add items to renderer.
        renderer = vtk.vtkRenderer()
        opacity=1.0
        #add_to_renderer(renderer=renderer, item=pipe, color='yellow', opacity=opacity)
        #add_to_renderer(renderer=renderer, item=pipe2, color='yellow', opacity=opacity)
        add_to_renderer(renderer=renderer, item=result_bool, color='red')

    ############################################################################
    # IMPLICIT BOOLEAN
    ############################################################################
    else:
        # We need to know the domain where the implicit function will be
        # evaulated. There is certainly other ways to achieve this. Here,
        # we simply get the bounds from the meshes. Also, we add a margin
        # to avoid artifacts close to the domain boundary.
        pipe = create_pipe(radius=radius, thickness=thickness, height=length)
        pipe2 = create_pipe(radius=radius2, thickness=thickness2, height=length2)
        pipe = transform_item(pipe, transform)
        pipe2 = transform_item(pipe2, transform2)
        bounds = pipe.GetOutput().GetBounds()
        bounds2 = pipe2.GetOutput().GetBounds()

        def applyMargin(bounds, margin):
            extent = [ bounds[1]-bounds[0],
                       bounds[3]-bounds[2],
                       bounds[5]-bounds[4] ]
            bounds = [ bounds[0]-extent[0]*margin, bounds[1]+extent[0]*margin,
                       bounds[2]-extent[1]*margin, bounds[3]+extent[1]*margin,
                       bounds[4]-extent[2]*margin, bounds[5]+extent[2]*margin ]
            return bounds
        bounds = applyMargin(bounds, margin=0.1)
        bounds2 = applyMargin(bounds2, margin=0.1)

        # The bounds of the combined object pipe+pipe2
        boundsCombo = [min(bounds[0], bounds2[0]),
                       max(bounds[1], bounds2[1]),
                       min(bounds[2], bounds2[2]),
                       max(bounds[3], bounds2[3]),
                       min(bounds[4], bounds2[4]),
                       max(bounds[5], bounds2[5])]

        # Let's create implicit functions for the pipes.
        pipeImp, innerImp, outerImp = create_pipe_implicit(radius=radius, thickness=thickness, height=length)
        pipeImp2, innerImp2, outerImp2 = create_pipe_implicit(radius=radius2, thickness=thickness2, height=length2)
        pipeImp.SetTransform(transform.GetInverse())
        pipeImp2.SetTransform(transform2.GetInverse())
        innerImp.SetTransform(transform.GetInverse())
        innerImp2.SetTransform(transform2.GetInverse())
        outerImp.SetTransform(transform.GetInverse())
        outerImp2.SetTransform(transform2.GetInverse())

        # Apply the intersection.
        difference = vtk.vtkImplicitBoolean()
        difference.AddFunction(pipeImp2)
        difference.AddFunction(outerImp)
        difference.SetOperationTypeToDifference()
        union = vtk.vtkImplicitBoolean()
        union.AddFunction(difference)
        union.AddFunction(pipeImp)
        union.SetOperationTypeToUnion()
        # This last operation is required to "cut through" the first pipe.
        difference = vtk.vtkImplicitBoolean()
        difference.AddFunction(union)
        difference.AddFunction(innerImp2)
        difference.SetOperationTypeToDifference()

        # Convert the implicit functions into surfaces.
        pipe = evaluate_implicit(implicit_function=pipeImp,
                                 resolution=100,
                                 bounds=bounds)
        pipe2 = evaluate_implicit(implicit_function=pipeImp2,
                                  resolution=100,
                                  bounds=bounds2)
        result = evaluate_implicit(implicit_function=difference,
                                  resolution=100,
                                  bounds=boundsCombo)

        # Add items to renderer.
        renderer = vtk.vtkRenderer()
        opacity=1.
        add_to_renderer(renderer=renderer, item=pipe, color='yellow', opacity=opacity, translate=[0,5,0])
        add_to_renderer(renderer=renderer, item=pipe2, color='yellow', opacity=opacity, translate=[0,5,0])
        add_to_renderer(renderer=renderer, item=result, color='red')

    # Create a renderer, render window, and interactor.
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetWindowName("Overlapping cylinders example")
    render_window.SetSize(1000,1000)
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    style = vtk.vtkInteractorStyleTrackballCamera()
    render_window_interactor.SetInteractorStyle(style)
    render_window_interactor.SetRenderWindow(render_window)
    # Add the actors to the scene.
    renderer.SetBackground(colors.GetColor3d("Gray"))
    # Render and interact.
    render_window.Render()
    render_window_interactor.Start()

if __name__ == '__main__':
    main()
