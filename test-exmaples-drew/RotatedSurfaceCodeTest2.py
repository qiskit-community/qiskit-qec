import vtk

def polygon_points_to_vtkPolyData(points):
    points = vtk.vtkPoints()
    for point in points:
        if len(point) == 2:
            points.InsertNextPoint(*point,0)
        else:
            points.InsertNextPoint(*point)
    
    # Create the polygon
    num_points = len(points)
    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(num_points)

    for k in range(num_points):
         polygon.GetPointIds().SetId(k, k)

    # Add the polygon to a list of polygons
    polygons = vtk.vtkCellArray()
    polygons.InsertNextCell(polygon)

    # Create a PolyData
    polygon_PolyData = vtk.vtkPolyData()
    polygon_PolyData.SetPoints(points)
    polygon_PolyData.SetPolys(polygons)

    return polygon_PolyData

class Code():
    def __init__(self):
        pass

class SubsystemCode(Code):
    def __init__(self):
        super().__init__()

class RotatedSurfaceCode(SubsystemCode):
    def __init__(self, shape, shape_points, cell_type, x_color, z_color, unit=1):
        


def main():
    pass
    # Create a shape and boundary for the rotated surface code:
    # Short cuts should be available for standard shapes and boundaries: e.g. square, rentangular
    # Need a method to discribe the boundaries (probably XXYXXYX and [P1,P2,P3,...p_k] where p_i are the coordinates of the transitions)

    # Will will start with a square d x d rotated surface code

    # Code Representation (Symplectic Matrix, ...)

    # TranslationMap (Matrix to Qubit index)
    # TranslationMap (Matrix to Local Coordinate)

    # Create shape of rotated surface code
    shape_points = [[0,0], [2,0], [2,-2], [0,-2]]
    shape = polygon_points_to_vtkPolyData(shape_points)
    
    bounds = shape.GetBounds()

    # Define the orgin cell type ('X' for and X stabilizer and Z for Z Stabilizer)
    cell_type = 'X'

    # Define the stabilizer colors
    x_color = 'Green'
    z_color = 'Red'

    code = RotatedSurfaceCode(shape, shape_points, cell_type, x_color, z_color)





    



