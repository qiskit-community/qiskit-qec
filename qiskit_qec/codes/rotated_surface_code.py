"""Module to define the rotated subsystem surface code."""


from typing import List

from qiskit_qec.exceptions import QiskitQECError


class RSSC:
    """Create a new rotated subsystem surface code (RSSC).

    The X and Z gauge operator lists are given as lists of supports.
    There is a consistent qubit ordering of these lists so that
    other modules can construct gate schedules for circuits.
    """

    def __init__(self, distance: int):
        """Initialize data for RSSC."""
        if distance % 2 == 0:
            raise Exception("require odd distance")
        if distance < 3:
            raise Exception("require positive distance > 2")
        d = distance
        k = 1
        n = int(d**2 + (d - 1) ** 2 / 2)
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

    def __str__(self) -> str:
        """Formatted string."""
        return f"[[{self.n}, {self.k}, {self.d}]] rotated surface code"

    def __repr__(self) -> str:
        """String representation."""
        val = str(self)
        val += f"\nx_gauges = {self.x_gauges}"
        val += f"\nz_gauges = {self.z_gauges}"
        val += f"\nx_stabilizers = {self.x_stabilizers}"
        val += f"\nz_stabilizers = {self.z_stabilizers}"
        val += f"\nlogical_x = {self.logical_x}"
        val += f"\nlogical_z = {self.logical_z}"
        val += f"\nx_boundary = {self.x_boundary}"
        val += f"\nz_boundary = {self.z_boundary}"
        return val

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
