import copy
import sys

import numpy as np
import qiskit_qec.utils.pauli_rep as pauli_rep
from qiskit.exceptions import QiskitError
from qiskit.utils.deprecation import deprecate_function
from qiskit_qec.structures.symplectic_matrix_old import SymplecticMatrixBase


class SymplecticMatrixDenseNumPy(SymplecticMatrixBase):
    XZ_ORDER = 0
    ZX_ORDER = 1

    def __init__(self, xz_matrix) -> None:
        """Symplectic Matrix class: None compressed bool matrices using numpy arrays

        Assumption: The assumption is that the X part of the parity check matrix is on
        the left and the Z part is on the right

        Args:
            matrix [SymplecticMatrixBase, ...] : [description]
        """

        # No dimensional check is being done here. Should decided on this
        self._matrix = self.assign(xz_matrix)

        super().__init__()

    def __array__(self):
        return self._matrix

    def __symplectic__(self):
        return self._matrix

    def joinxz(self, x_matrix=None, z_matrix=None):
        """Join matrices into a single symplectic matrix

        Args:
            x_matrix ([type], optional): [description]. Defaults to None.
            z_matrix ([type], optional): [description]. Defaults to None.

            Converts x_matrix and Z_matrix matrices to np.ndarray's and then
            joins matrices together. [x_matrix|z_matrix]
        """
        x_matrix = np.asarray(x_matrix, dtype=bool)
        z_matrix = np.asarray(z_matrix, dtype=bool)
        if x_matrix.shape != z_matrix.shape:
            print("Error: X and Z matrices must have the same shape")
            exit

        return np.hstack(x_matrix, z_matrix)

    def assign(self, matrix):
        """Assign a input symplectic matrix to a Dense NumPy Symplectic matrix
        Includes any necesary checks

        Args:
            matrix ([type], optional): [description].
        """
        # asarray is used so that we do not have to copy the matrix if possible
        # This is needed to allow PauliList to work correctly
        as_matrix = np.asarray(matrix, bool)
        # Do a simple check for a proper conversion.
        if not isinstance(as_matrix[0][0], bool):
            print("Error: unable to convert matruix input")
            exit

        return as_matrix

    def copy(self):
        """Make a deep copy of the current sympletic matrix"""
        ret = copy.copy(self)
        ret._matrix = self._matrix
        ret._properties = self._properties.copy()
        return ret

    @property
    def num_qubits(self):
        return self._matrix.shape[1] >> 1

    @property
    def num_paulis(self):
        return self._matrix.shape[0]

    @property
    def x(self):
        """X part of symplectic matrix

        Returns:
            numpy.ndarray: X part of symplectic matrix
        """
        return SymplecticMatrixDenseNumPy(
            xz_matrix=self._matrix[0 : self.num_paulis, 0 : self.num_qubits]
        )

    @property
    def z(self):
        """Z part of symplectic matrix

        Returns:
            numpy.ndarray: Z part of symplectic matrix
        """
        return SymplecticMatrixDenseNumPy(
            xz_matrix=self._matrix[0 : self.num_paulis, self.num_qubits : 2 * self.num_qubits]
        )

    @property
    def num_y(self):
        """Count the number of Y Paulis in each Pauli and return an array of result"""
        return np.sum(np.logical_and(self.x, self.z), axis=1)

    @staticmethod
    def stypes():
        return SymplecticMatrixBase.STYPES

    @staticmethod
    def classnames(key=None):
        if key is not None:
            return SymplecticMatrixBase.CLASSNAMES[key]
        return SymplecticMatrixBase.CLASSNAMES

    @staticmethod
    def logical_and(x1, x2, **kwargs):
        return SymplecticMatrixDenseNumPy(
            xz_matrix=np.logical_and(x1._matrix, x2._matrix, **kwargs)
        )

    @staticmethod
    def logical_xor(x1, x2, **kwargs):
        return SymplecticMatrixDenseNumPy(
            xz_matrix=np.logical_xor(x1._matrix, x2._matrix, **kwargs)
        )

    @staticmethod
    def sum(a, **kwargs):
        return np.sum(a._matrix, **kwargs)

    def __and__(self, other):
        return SymplecticMatrixDenseNumPy(xz_matrix=(self._matrix & other._matrix))

    def __xor__(self, other):
        return SymplecticMatrixDenseNumPy(xz_matrix=(self._matrix ^ other._matrix))

    def __eq__(self, other):
        self._matrix = other._matrix
        self._properties = other._properties

    @staticmethod
    def hstack(tup):
        """Stack arrays in sequence horizontally (column wise)

        This function mimicks numpys hstack for this specific symplectic matrix representation

        Args:
            Tuple of SymplecticMatrixDenseNumPy matrices

        Returns:
            SymplecticMatrixDenseNumPy: The matrix formed by stacking the given matrices
        """

        return SymplecticMatrixDenseNumPy(xz_matrix=np.hstack((item._matrix for item in tup)))

    @staticmethod
    def istack(matrix, size, interleave=False):
        """Vertically stack array of vectors.

        Args:
            SymplecticMatrixDenseNumPy matrice
            size (int)
            interleave (boolean): defulat is False


        matrix = [v_1
                v_2
                ...
                v_k]

        istack(matrix, r, interleave=False) gives r vertically stacked copies of array with no iterleaving

        output = [v_1
                  v_2
                  ...
                  v_k

                  ... r times

                  v_1
                  v_2
                  ...
                  v_k]

        istack(matrix, r, interleave=True) gives r vertically stacked copies of array with with iterleaving

        output = [v_1
                  v_1
                  ... r copies
                  v_1

                  ...

                  v_k
                  v_k
                  ... r copies
                  v_k]

        """
        if size == 1:
            return matrix

        if interleave:
            return np.hstack(size * [matrix._matrix]).reshape(
                (size * len(matrix._matrix),) + matrix._matrix.shape[1:]
            )

        return np.vstack(size * [matrix._matrix]).reshape(
            (size * len(matrix._matrix),) + matrix._matrix.shape[1:]
        )

    @staticmethod
    def count_y(matrix=None, x_matrix=None, z_matrix=None):
        """Count the number of Y Paulis"""
        if matrix is not None:
            return np.sum(np.logical_and(matrix.x, matrix.z), axis=1)
        else:
            return np.sum(np.logical_and(x_matrix, z_matrix), axis=1)

    # This methiod should be deprecated and replaced with _from_matrix
    @staticmethod
    def _from_array(z, x, phase=0, input_format=pauli_rep.DEFAULT_EXTERNAL_PAULI_REP_FORMAT):
        """Convert array data to BasePauli data."""
        # Moved from BasePauli class

        array_z = np.atleast_2d(np.asarray(z, dtype=bool))
        array_x = np.atleast_2d(np.asarray(x, dtype=bool))

        if array_z.ndim != 2:
            raise QiskitError("Invalid Pauli z vector shape.")
        if array_x.ndim != 2:
            raise QiskitError("Invalid Pauli x vector shape.")
        if array_z.shape != array_x.shape:
            raise QiskitError("z and x vectors are different size.")

        # Convert group phase convention to internal ZX-phase conversion.
        # Internal Pauli representation is '-iZX' frmt
        # External Pauli representation is '-iYZX' frmt and is assumed as input
        # when referencing is not possible for the phase.
        base_phase = pauli_rep.change_rep(
            phase,
            y_count=SymplecticMatrixDenseNumPy.count_y(x_matrix=array_x, z_matrix=array_z),
            input_format=input_format,
            output_format=pauli_rep.INTERNAL_PAULI_REP_FORMAT,
        )

        return array_z, array_x, base_phase

    @staticmethod
    def _from_matrix(matrix, phase=0, input_format=pauli_rep.DEFAULT_EXTERNAL_PAULI_REP_FORMAT):
        """Convert array data to BasePauli data."""

        matrix = np.atleast_2d(np.asarray(matrix))

        if matrix.shape[1] % 2 != 0:
            raise QiskitError("Invalid matrix shape")

        # Convert group phase convention to internal ZX-phase conversion.
        # Internal Pauli representation is '-iZX' frmt
        # External Pauli representation is '-iYZX' frmt and is assumed as input
        # when referencing is not possible for the phase.
        base_phase = pauli_rep.change_rep(
            phase,
            y_count=SymplecticMatrixDenseNumPy.count_y(matrix),
            input_format=input_format,
            output_format=pauli_rep.INTERNAL_PAULI_REP_FORMAT,
        )

        return matrix, base_phase

    @staticmethod
    def _to_matrix(matrix, phase=0, group_phase=False, sparse=False):
        """Return the matrix matrix from symplectic representation.
        The Pauli is defined as :math:`P = (-i)^{phase + z.x} * Z^z.x^x`
        where ``array = [x, z]``.
        Args:
            matrix (SymplecticMatrixDenseNumPy): Symplectic Matrix representation of Pauli
            phase (int): Pauli phase.
            group_phase (bool): Optional. If True use group-phase convention
                                instead of BasePauli ZX-phase convention.
                                (default: False).
            sparse (bool): Optional. Of True return a sparse CSR matrix,
                        otherwise return a dense Numpy array
                        (default: False).
        Returns:
            array: if sparse=False.
            csr_matrix: if sparse=True.
        """
        # Moved from BasePauli
        num_qubits = matrix.z.size
        # Convert to zx_phase

        if group_phase:
            phase += np.sum(x & z)
            phase %= 4
        dim = 2 ** num_qubits
        twos_array = 1 << np.arange(num_qubits)
        x_indices = np.asarray(x).dot(twos_array)
        z_indices = np.asarray(z).dot(twos_array)
        indptr = np.arange(dim + 1, dtype=np.uint)
        indices = indptr ^ x_indices
        if phase:
            coeff = (-1j) ** phase
        else:
            coeff = 1
        data = np.array([coeff * (-1) ** (bin(i).count("1") % 2) for i in z_indices & indptr])
        if sparse:
            # Return sparse matrix
            from scipy.sparse import csr_matrix

            return csr_matrix((data, indices, indptr), shape=(dim, dim), dtype=complex)
        # Build dense matrix using csr frmt
        mat = np.zeros((dim, dim), dtype=complex)
        for i in range(dim):
            mat[i][indices[indptr[i] : indptr[i + 1]]] = data[indptr[i] : indptr[i + 1]]
        return mat

    @staticmethod
    def _to_label(z, x, phase, group_phase=False, full_group=True, return_phase=False):
        """Return the label string for a Pauli.
        Args:
            z (array): The symplectic representation z vector.
            x (array): The symplectic representation x vector.
            phase (int): Pauli phase.
            group_phase (bool): Optional. If True use group-phase convention
                                instead of BasePauli ZX-phase convention.
                                (default: False).
            full_group (bool): If True return the Pauli label from the full Pauli group
                including complex coefficient from [1, -1, 1j, -1j]. If
                False return the unsigned Pauli label with coefficient 1
                (default: True).
            return_phase (bool): If True return the adjusted phase for the coefficient
                of the returned Pauli label. This can be used even if
                ``full_group=False``.
        Returns:
            str: the Pauli label from the full Pauli group (if ``full_group=True``) or
                from the unsigned Pauli group (if ``full_group=False``).
            Tuple[str, int]: if ``return_phase=True`` returns a tuple of the Pauli
                            label (from either the full or unsigned Pauli group) and
                            the phase ``q`` for the coefficient :math:`(-i)^(q + x.z)`
                            for the label from the full Pauli group.
        """
        num_qubits = z.size
        coeff_labels = {0: "", 1: "-i", 2: "-", 3: "i"}
        label = ""
        for i in range(num_qubits):
            if not z[num_qubits - 1 - i]:
                if not x[num_qubits - 1 - i]:
                    label += "I"
                else:
                    label += "X"
            elif not x[num_qubits - 1 - i]:
                label += "Z"
            else:
                label += "Y"
                if not group_phase:
                    phase -= 1
        phase %= 4
        if phase and full_group:
            label = coeff_labels[phase] + label
        if return_phase:
            return label, phase
        return label
