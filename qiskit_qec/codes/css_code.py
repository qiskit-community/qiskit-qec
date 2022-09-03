"""Module to define a class for CSS codes."""


class CSSCode:
    """Class defining a CSS code."""

    def __init__(
        self, x_gauges, z_gauges, n, k, d, x_stabilizers, z_stabilizers, logical_z, logical_x
    ):
        """Initialize data for CSS code."""
        self.d = d
        self.k = k
        self.n = n
        self.x_gauges = x_gauges
        self.z_gauges = z_gauges
        self.x_stabilizers = x_stabilizers
        self.z_stabilizers = z_stabilizers
        self.logical_z = logical_z
        self.logical_x = logical_x
        self.x_gauge_matrix = self._to_matrix(self.x_gauges)
        self.z_gauge_matrix = self._to_matrix(self.z_gauges)
        self.x_stabilizer_matrix = self._to_matrix(self.x_stabilizers)
        self.z_stabilizer_matrix = self._to_matrix(self.z_stabilizers)

    def _to_matrix(self, operators):
        """Compute the binary parity check matrix for a list of operators.

        The operators are lists of supports on self.n qubits.
        """
        mat = []
        for g in operators:
            row = [0] * self.n
            for i in g:
                if i != -1:
                    row[i] = 1
            mat.append(row)
        return mat

    def _syndrome(self, bitstring, operators):
        """Compute the syndrome of a length self.n bit string.

        operators is a list of stabilizers.
        """
        if len(bitstring) != self.n:
            raise Exception("bitstring has wrong length")
        syndrome = []
        for s in operators:
            value = 0
            for i in s:
                if i != -1:
                    value += bitstring[i]
            syndrome.append(value % 2)
        return syndrome

    def x_syndrome(self, bitstring):
        """Compute the X syndrome of a length self.n bit string."""
        return self._syndrome(bitstring, self.x_stabilizers)

    def z_syndrome(self, bitstring):
        """Compute the Z syndrome of a length self.n bit string."""
        return self._syndrome(bitstring, self.z_stabilizers)

    def logical_x_error(self, bitstring):
        """Test for a logical X error."""
        return self._syndrome(bitstring, self.logical_z)

    def logical_z_error(self, bitstring):
        """Test for a logical Z error."""
        return self._syndrome(bitstring, self.logical_x)
