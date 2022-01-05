# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import numpy as np

from qiskit.exceptions import QiskitError

from qiskit_qec.geometry.shape import Shape
from qiskit_qec.geometry.plane import Plane
from qiskit_qec.factory.subsystem.css.css_subsystem_factory import CSSSubSystemFactory

from qiskit_qec.exceptions import QiskitQECError

class RSSCFactory(CSSSubSystemFactory):
    """Rotated subsystem surface code (RSSC) code factory.

    The rotated subsystem surace code (RSSC) family consists of codes
    constructed by taking a finite subset weight 3 gauge operators from
    the weight thee operators defined by the infinite tiling given by 
    tiling the plane with one of the tiles:

    (Diamond or Rotated orientation)

             o               o    
            /X\             /Z\  
           o---o           o---o  
          /|   |\         /|   |\     
         oZ|   |Zo   or  oX|   |Xo   
          \|   |/         \|   |/  
           o---o           o---o      
            \X/             \Z/    
             o               o   

    where the gauge operators are weight three X and Z type defined
    around the square with X opposite X and Z opposite Z.

    Or equivelently, it may be constructured by one of the rotated tiles:

    (Square or non-Rotated orientation)

        o--o--o          o--o--o
        |X/ \Z|          |X/ \Z|
        o/   \o     or   o/   \o
        |\   /|          |\   /|
        |Z\ /X|          |Z\ /X|
        o--o--o          o--o--o



    Qubits are located on the lattice Z^2 and each qubit is assigned the
    natural coordinate associated with the Z^2 lattice.

    The tile orientation passed to the code factory as a specification.

    factory = RSSCCodeFactory(...)
    factory.set_spec(orientation=...)

    where orientation is one of (diamond or rotated) or (square or non-rotated). The default is
    rotated. 

    Internal computations are done in the non-rotated tiling as the tiles are located 

    Selection of the subtiling

    The shape of the subtiling:

    One approach to selecting a subtiling is to provide a bounding shape from which
    the subtiling is cut from the full tiling (we assume that the subtiling is finite and
    so the inside of the shape is choosen for the subtiling where the inside of the shape is
    defined as the region with finite area/volume)

    The shape is encoded in a Shape class. We assume that the shape and the full latiice live 
    on the same manifold.

    The selection of the subtiling defining the code can be made in various ways:

        1) Setting the distance d. This will assume that the shape of the subtiling is 
        "square", select type is "Intersect" with k=2. A boundary may alos be provided.

        factory.set_spec(d=d, shape="square", k=2,...)
        
        2) Setting a shape. Gauge operators and their quibts within the bounds of the (closed) 
        shape are included. Those gauge opertors and their qubits that are not strictly bounded
        by the provided shape are included (perhaps partially) according to the select style
        choosen:

            "Cut" : Only gauge opertors and there qubits that are fully bounded by the shape
                    are included
            "Intersect" : Include gauge operators and their qubits that have at least k of those 
                     qubits bounded by the shape. Some gauge operators may 
                     have reduced weight as some of their qubits may not meet the selection 
                     critera
            "Other" : Any other selection method that is possible to choose the operators.

        factory.set_spec(shape = shape, select=<select_type>, k= ...)

    The X and Z gauge operator lists are given as lists of supports.
    There is a consistent qubit ordering of these lists so that
    other modules can construct gate schedules for circuits.


    factory = RSSCCodeFactory(...)
    factory.set_spec(orientation=...)

    """
    
    def __init__(self, stype="numpy"):
        super().__init__(stype=stype)

    def set_specs(self, shape=None, d=None,  select=None, k=3, boundary=None, orientation = "rotated"):
        # set parameters for codes to be manufactured

        if shape is None:
            if d is None:
                raise QiskitError("A value of d must be provided if not providing an shape")

            assert d % 2 == 1, "The distance must be odd"
            assert d > 2, "The distance must be at least 3"

            # Create a square shape 
            # Note: Will work in the non-rotated tiling and switch to the rotated orientation after 
            # any computations if required
            # The origin is chosen to be (0,-1) as the tiles are centered on the 2Z+1 2Z+1 integer
            # lattice and this origin positions itself on the middle lower qubit of a tile.
            manifold = Plane()
            shape = Shape.square(origin=(0,-1), direction=(1,1), length=d-1, manifold=manifold, dtype=int)
            ebounds = shape.bounds.copy()
            amount = np.array([2,2])/2
            ebounds.expand(amount)

            

            num_quibits = d*d + (((d-1)*(d-1)) >>1)

    def get(self):
        # get a manufactured code
        pass


    def order(specs):





    
    def __init__(self, distance):
        """Initialize data for RSSC."""
        if distance % 2 == 0:
            raise Exception("require odd distance")
        if distance < 3:
            raise Exception("require positive distance > 2")
        d = distance
        k = 1
        n = int(d ** 2 + (d - 1) ** 2 / 2)
        x_gauges, self.x_orientations = self._x_gauge_generators(d)
        z_gauges, self.z_orientations = self._z_gauge_generators(d)
        x_stabilizers = self._x_stabilizer_generators(d)
        z_stabilizers = self._z_stabilizer_generators(d)
        logical_z = self._logical_z(d)
        logical_x = self._logical_x(d)
        self.x_boundary = self._x_boundary_qubits(d)
        self.z_boundary = self._z_boundary_qubits(d)
        super().__init__(
            x_gauges,
            z_gauges,
            n,
            k,
            d,
            x_stabilizers,
            z_stabilizers,
            logical_z,
            logical_x,
        )

    def to_index(self, row, col, d):
        """Map a coordinate (row, col) to a qubit index.

        Qubits are indexed from left to right across each row,
        beginning in the top row. The qubits on the faces of
        the lattice are indexed starting from d*d in the same way.
        """
        if row < 0 or row > 2 * (d - 1):
            raise Exception("row out of bounds")
        if col < 0 or col > 2 * (d - 1):
            raise Exception("col out of bounds")
        if row % 2 != col % 2:
            raise Exception("(row,col) invalid")
        if row % 2 == 0 and col % 2 == 0:
            return int(row / 2 * d + col / 2)
        else:
            if (row + col) % 4 != 0:
                raise Exception("(row,col) invalid")
            offset = ((row - 1) / 2 + 1) % 2
            return int(d * d + (row - 1) / 2 * (d - 1) / 2 + (col - 1 - 2 * offset) / 4)

    def _x_gauge_generators(self, d):
        """Return a list of X gauge generators and orientations.

        Each X gauge generator is a list of qubit indices in clockwise
        order around the face starting with the left-most qubit.
        The index -1 is substituted when a qubit does not exist because
        the gauge operator is on the boundary of the lattice.
        """
        x_gauges = []
        x_orientations = []  # 0 = up, 1 = down
        for i in range(d):
            for j in range(d - 1):
                if i % 2 == j % 2:  # triangle points up
                    try:
                        tip = self.to_index(2 * i - 1, 2 * j + 1, d)
                    except Exception:
                        tip = -1
                    x_gauges.append(
                        [
                            self.to_index(2 * i, 2 * j, d),
                            tip,
                            self.to_index(2 * i, 2 * j + 2, d),
                        ]
                    )
                    x_orientations.append(0)
                else:  # triangle points down
                    try:
                        tip = self.to_index(2 * i + 1, 2 * j + 1, d)
                    except Exception:
                        tip = -1
                    x_gauges.append(
                        [
                            self.to_index(2 * i, 2 * j, d),
                            self.to_index(2 * i, 2 * j + 2, d),
                            tip,
                        ]
                    )
                    x_orientations.append(1)
        return x_gauges, x_orientations

    def _z_gauge_generators(self, d):
        """Return a list of Z gauge generators and orientations.

        Each Z gauge generator is a list of qubit indices in clockwise
        order around the face starting with the top-most qubit.
        The index -1 is substituted when a qubit does not exist because
        the gauge operator is on the boundary of the lattice.
        """
        z_gauges = []
        z_orientations = []  # 0 = left, 1 = right
        for i in range(d - 1):
            for j in range(d):
                if i % 2 == j % 2:  # triangle points left
                    try:
                        tip = self.to_index(2 * i + 1, 2 * j - 1, d)
                    except Exception:
                        tip = -1
                    z_gauges.append(
                        [
                            self.to_index(2 * i, 2 * j, d),
                            self.to_index(2 * i + 2, 2 * j, d),
                            tip,
                        ]
                    )
                    z_orientations.append(0)
                else:  # triangle points right
                    try:
                        tip = self.to_index(2 * i + 1, 2 * j + 1, d)
                    except Exception:
                        tip = -1
                    z_gauges.append(
                        [
                            self.to_index(2 * i, 2 * j, d),
                            tip,
                            self.to_index(2 * i + 2, 2 * j, d),
                        ]
                    )
                    z_orientations.append(1)
        return z_gauges, z_orientations

    def _x_stabilizer_generators(self, d):
        """Return a list of X stabilizer generators.

        Each X stabilizer generator is a list of qubit indices. The index
        -1 is substituted when a qubit does not exist because the stabilizer
        is on the boundary of the lattice.
        """
        x_gauges, x_orientations = self._x_gauge_generators(d)
        x_stabilizers = []
        for i in range(d):
            for j in range(d - 1):
                if i == 0 and j % 2 == 1:
                    x_stabilizers.append(x_gauges[(d - 1) * i + j])
                elif i == d - 1 and j % 2 == 0:
                    x_stabilizers.append(x_gauges[(d - 1) * i + j])
                elif i % 2 == j % 2:
                    x_stabilizers.append(
                        x_gauges[(d - 1) * i + j] + x_gauges[(d - 1) * (i + 1) + j]
                    )
        return x_stabilizers

    def _z_stabilizer_generators(self, d):
        """Return a list of Z stabilizer generators.

        Each Z stabilizer generator is a list of qubit indices.
        The index -1 is substituted when a qubit does not exist because
        the stabilizer is on the boundary of the lattice.
        """
        z_gauges, z_orientations = self._z_gauge_generators(d)
        z_stabilizers = []
        for i in range(d - 1):
            for j in range(d):
                if i % 2 == 1 and j == 0:
                    z_stabilizers.append(z_gauges[d * i + j])
                elif i % 2 == 0 and j == d - 1:
                    z_stabilizers.append(z_gauges[d * i + j])
                elif i % 2 == j % 2:
                    z_stabilizers.append(z_gauges[d * i + j] + z_gauges[d * i + j + 1])
        return z_stabilizers

    def _logical_z(self, d):
        """Return the support of the logical Z operators."""
        return [[i for i in range(d)]]

    def _logical_x(self, d):
        """Return the support of the logical X operators."""
        return [[d * i for i in range(d)]]

    def _z_boundary_qubits(self, d):
        """Return a list of singletons each containing a Z boundary qubit.

        A boundary qubit is a qubit adjacent to the Z-type boundary on which
        a logical X chain terminates.
        """
        z_boundary = []
        for i in range(d):
            z_boundary.append([i])
        for i in range(d):
            z_boundary.append([d * (d - 1) + i])
        return z_boundary

    def _x_boundary_qubits(self, d):
        """Return a list of singletons each containing an X boundary qubit.

        A boundary qubit is a qubit adjacent to the X-type boundary on which
        a logical Z chain terminates.
        """
        x_boundary = []
        for i in range(d):
            x_boundary.append([d * i])
        for i in range(d):
            x_boundary.append([d * i + d - 1])
        return x_boundary
