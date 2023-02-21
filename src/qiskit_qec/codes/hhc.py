"""Define the heavy hexagon compass code."""

from typing import List

from qiskit_qec.exceptions import QiskitQECError


class HHC:
    """Heavy hexagon compass code (HHC) data.

    The X and Z gauge operator lists are given as lists of supports.
    There is a consistent qubit ordering of these lists so that
    we can construct gate schedules for circuits.
    """

    def __init__(self, distance: int):
        """Initialize data for HHC."""
        if distance % 2 == 0:
            raise QiskitQECError("even distance not implemented")
        if distance < 3:
            raise QiskitQECError("require distance > 2")
        self.d = distance
        self.k = 1
        self.n = int(distance**2)
        self.x_gauges = self._x_gauge_generators(self.d)
        self.z_gauges = self._z_gauge_generators(self.d)
        self.x_stabilizers = self._x_stabilizer_generators(self.d)
        self.z_stabilizers = self._z_stabilizer_generators(self.d)
        self.logical_z = self._logical_z(self.n)
        self.logical_x = self._logical_x(self.n)
        self.x_boundary = self._x_boundary_qubits(self.d)
        self.z_boundary = self._z_boundary_qubits(self.d)

    def __str__(self) -> str:
        """Formatted string."""
        return f"[[{self.n}, {self.k}, {self.d}]] heavy-hexagon compass code"

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

    @staticmethod
    def to_index(row: int, col: int, d: int) -> int:
        """Map a coordinate (row, col) to a qubit index.

        Qubits are indexed from left to right across each row,
        beginning in the top row.
        """
        if row < 0 or row > d - 1:
            raise Exception("row out of bounds")
        if col < 0 or col > d - 1:
            raise Exception("col out of bounds")
        return int(row * d + col)

    def _x_gauge_generators(self, d: int) -> List[List[int]]:
        """Return a list of X gauge generators.

        We choose to define the code such that the X gauge generators
        have weight two and are oriented vertically on the lattice.
        Each X gauge generator is a list of qubit indices ordered from
        top to bottom.
        """
        x_gauges = []
        for i in range(d - 1):
            for j in range(d):
                x_gauges.append([self.to_index(i, j, d), self.to_index(i + 1, j, d)])
        return x_gauges

    def _z_gauge_generators(self, d: int) -> List[List[int]]:
        """Return a list of Z gauge generators.

        We choose to define the code such that the Z gauge generators
        have weight four. Each Z gauge generator is a list of qubit indices
        in row order from the upper left to lower right.
        """
        z_gauges = []
        for i in range(d + 1):
            for j in range(d - 1):
                if i % 2 == j % 2:
                    if i == 0:
                        z_gauges.append([self.to_index(i, j, d), self.to_index(i, j + 1, d)])
                    elif i == d:
                        z_gauges.append(
                            [self.to_index(i - 1, j, d), self.to_index(i - 1, j + 1, d)]
                        )
                    else:
                        z_gauges.append(
                            [
                                self.to_index(i - 1, j, d),
                                self.to_index(i - 1, j + 1, d),
                                self.to_index(i, j, d),
                                self.to_index(i, j + 1, d),
                            ]
                        )
        return z_gauges

    def _x_stabilizer_generators(self, d: int) -> List[List[int]]:
        """Return a list of X stabilizer generators.

        Each X stabilizer generator is a list of qubit indices.
        Each is a product of at most two X gauge generators.
        """
        x_stabilizers = []
        for i in range(d - 1):
            for j in range(d):
                if i % 2 == j % 2:
                    if j == d - 1:
                        x_stabilizers.append([self.to_index(i, j, d), self.to_index(i + 1, j, d)])
                    else:
                        x_stabilizers.append(
                            [
                                self.to_index(i, j, d),
                                self.to_index(i, j + 1, d),
                                self.to_index(i + 1, j, d),
                                self.to_index(i + 1, j + 1, d),
                            ]
                        )
                if i % 2 != j % 2:
                    if j == 0:
                        x_stabilizers.append([self.to_index(i, j, d), self.to_index(i + 1, j, d)])
        return x_stabilizers

    def _z_stabilizer_generators(self, d: int) -> List[List[int]]:
        """Return a list of Z stabilizer generators.

        Each Z stabilizer generator is a list of qubit indices.
        Each is a product of (d+2)/2 Z gauge generators.
        """
        z_stabilizers = []
        for j in range(d - 1):
            new_stabilizer = []
            for i in range(d):
                new_stabilizer.append(self.to_index(i, j, d))
                new_stabilizer.append(self.to_index(i, j + 1, d))
            z_stabilizers.append(new_stabilizer)
        return z_stabilizers

    @staticmethod
    def _logical_z(n: int) -> List[List[int]]:
        """Return the support of the logical Z operators."""
        return [list(range(n))]

    @staticmethod
    def _logical_x(n: int) -> List[List[int]]:
        """Return the support of the logical X operators."""
        return [list(range(n))]

    @staticmethod
    def _z_boundary_qubits(d: int) -> List[List[int]]:
        """Return a list of singletons each containing a Z boundary qubit.

        A boundary qubit is a qubit adjacent to the Z-type boundary on which
        a logical X chain terminates.
        """
        z_boundary = []
        for i in range(d):
            z_boundary.append([d * i])
        for i in range(d):
            z_boundary.append([d * i + d - 1])
        return z_boundary

    @staticmethod
    def _x_boundary_qubits(d: int) -> List[List[int]]:
        """Return a list of singletons each containing an X boundary qubit.

        A boundary qubit is a qubit adjacent to the X-type boundary on which
        a logical Z chain terminates.
        """
        x_boundary = []
        for i in range(d):
            x_boundary.append([i])
        for i in range(d):
            x_boundary.append([d * (d - 1) + i])
        return x_boundary
