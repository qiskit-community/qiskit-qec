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

"""Matrix ops."""

from typing import List, Tuple, Union
import numpy as np
from qiskit import QiskitError
from qiskit_qec.arithmetic.modn import gcd_ext, quo, ann, unit


def create_lambda_matrix(n: int) -> np.ndarray:
    r"""Creates the GF(2) Lambda matrix

    [0 I_n]
    [I_n 0]

    Args:
        n : Size of identity matrix to use

    Returns:
        GF(2) Lambda matrix used in the GF(2) symplectic product

    Raises:
        QiskitError: n must be a positive integer
        QiskitError: n must be an integer

    Examples:
        >>> x = create_lambda_matrix(2)
        >>> x.astype(int)
        array([[0, 0, 1, 0],
               [0, 0, 0, 1],
               [1, 0, 0, 0],
               [0, 1, 0, 0]])

    See Also: _create_lambda_matrix
    """
    if not n > 0:
        raise QiskitError("n must be a positive integer")
    if not isinstance(n, (int, np.integer)):
        raise QiskitError("n must be an integer")

    return _create_lambda_matrix(n)


def _create_lambda_matrix(n: int):
    r"""Creates the GF(2) Lambda matrix

    [0 I_n]
    [I_n 0]

    Args:
        n : Size of identity matrix to use

    Returns:
        GF(2) Lambda matrix used in the GF(2) symplectic product

    Examples:
        >>> x = _create_lambda_matrix(2)
        >>> x.astype(int)
        array([[0, 0, 1, 0],
               [0, 0, 0, 1],
               [1, 0, 0, 0],
               [0, 1, 0, 0]])

    See Also:
    create_lambda_matrix
    """
    lambda_matrix = np.eye(n << 1, k=n, dtype=bool)
    lambda_matrix[n:, 0:n] = np.identity(n, dtype=bool)

    return lambda_matrix


def augment_mat(matrix: np.ndarray, pos: str = "right") -> np.ndarray:
    """Augments a given matrix with the identify matrix. The position of the
    identity matrix is given by the optional position argument
    left    ->  [M | I]

    right   ->  [I | M]

    top     ->  [I]
                [M]

    bottom  ->  [M]
                [I]

    Args:
        matrix : Matrix to agument
        pos : Position of augmentation. One of "top", "bottom"
            "left" or "right"

    Raises:
        QiskitError: Input array must have two dimensions
        QiskitError: Unknown position

    Returns:
        Matrix augmented with identity

    Examples:
        >>> matrix = numpy.array([[0,1,0],[0,0,1],[1,1,0]], dtype=numpy.bool_)
        >>> x = augment_mat(matrix)
        >>> x
        array([[False,  True, False,  True, False, False],
            [False, False,  True, False,  True, False],
            [ True,  True, False, False, False,  True]])
        >>> matrix = numpy.array([[0,1,0],[0,0,1],[1,1,0]], dtype=numpy.bool_)
        >>> x = augment_mat(matrix, pos="bottom")
        >>> x
        array([[False,  True, False],
            [False, False,  True],
            [ True,  True, False],
            [ True, False, False],
            [False,  True, False],
            [False, False,  True]])

    See Also:
    _augment_mat
    """
    matrix = np.asarray(matrix)
    if not matrix.ndim == 2:
        raise QiskitError("Input array must have two dimensions")
    if pos not in ["top", "bottom", "right", "left"]:
        raise QiskitError("Unknown position")
    return _augment_mat(matrix, pos)


def _augment_mat(matrix: np.array, pos: str) -> np.ndarray:
    """Augments a given matrix with the identify matrix. The position of the
    identity matrix is given by the optional position argument

    left    ->  [M | I]

    right   ->  [I | M]

    top     ->  [I]
                [M]

    bottom  ->  [M]
                [I]

    Args:
        matrix : Matrix to agument
        pos : Position of augmentation. One of "top", "bottom"
            "left" or "right"

    Returns:
        Matrix augmented with identity

    Examples:
        >>> matrix = np.array([[0,1,0],[0,0,1],[1,1,0]], dtype=numpy.bool_)
        >>> x = _augment_mat(matrix)
        >>> x
        array([[False,  True, False,  True, False, False],
            [False, False,  True, False,  True, False],
            [ True,  True, False, False, False,  True]])

        >>> matrix = np.array([[0,1,0],[0,0,1],[1,1,0]], dtype=numpy.bool_)
        >>> x = _augment_mat(matrix, pos="bottom")
        >>> x
        array([[False,  True, False],
            [False, False,  True],
            [ True,  True, False],
            [ True, False, False],
            [False,  True, False],
            [False, False,  True]])

    See Also:
    augment_mat
    """

    result = None
    if pos == "top":
        _id = np.identity(matrix.shape[1], dtype=matrix.dtype)
        result = np.vstack((_id, matrix))
    elif pos == "bottom":
        _id = np.identity(matrix.shape[1], dtype=matrix.dtype)
        result = np.vstack((matrix, _id))
    elif pos == "left":
        _id = np.identity(matrix.shape[0], dtype=matrix.dtype)
        result = np.hstack((_id, matrix))
    elif pos == "right":
        _id = np.identity(matrix.shape[0], dtype=matrix.dtype)
        result = np.hstack((matrix, _id))
    return result


def rref(matrix: np.ndarray) -> np.ndarray:
    """Computes the Row Reduced Echelon Form for a GF(2) matrix.

    Args:
        matrix : Input GF(2) matrix

    Returns:
        Row Reduced Echelon Form of input matrix

    Examples:
        >>> matrix = numpy.array([[1,0,0,1,0,0,1,0],
                                  [0,1,1,1,0,0,0,1],
                                  [1,1,1,0,1,0,0,0],
                                  [1,0,0,1,0,1,0,1]], dtype=np.bool_)
        >>> r = rref(matrix)
        >>> r.astype(int)
        array([[1, 0, 0, 1, 0, 0, 1, 0],
            [0, 1, 1, 1, 0, 0, 0, 1],
            [0, 0, 0, 0, 1, 0, 1, 1],
            [0, 0, 0, 0, 0, 1, 1, 1]])

    See Also:
    _rref, rref_complete, _rref_complete
    """
    _, rref_mat, _, _ = rref_complete(matrix)
    return rref_mat


def _rref(matrix: np.array) -> np.ndarray:
    """Computes the Row Reduced Echelon Form for a GF(2) matrix.

    Args:
        matrix : Input GF(2) matrix

    Returns:
        Row Reduced Echelon Form of input matrix

    Examples:
        >>> matrix = numpy.array([[1,0,0,1,0,0,1,0],
                                  [0,1,1,1,0,0,0,1],
                                  [1,1,1,0,1,0,0,0],
                                  [1,0,0,1,0,1,0,1]], dtype=np.bool_)
        >>> r = _rref(matrix)
        >>> r.astype(int)
        array([[1, 0, 0, 1, 0, 0, 1, 0],
            [0, 1, 1, 1, 0, 0, 0, 1],
            [0, 0, 0, 0, 1, 0, 1, 1],
            [0, 0, 0, 0, 0, 1, 1, 1]])

    See Also:
    rref, _rref_complete, rref_comnplete
    """
    _, rref_mat, _, _ = _rref_complete(matrix)
    return rref_mat


def rank(matrix: np.ndarray) -> int:
    """Computes the rank of the GF(2) matrix

    Args:
        matrix: Input GF(2) matrix

    Returns:
        Rank of input matrix

    Examples:
        >>> matrix = numpy.array([[1,0,0,1,0,0,1,0],
                                  [0,1,1,1,0,0,0,1],
                                  [1,1,1,0,1,0,0,0],
                                  [1,0,0,1,0,1,0,1]], dtype=np.bool_)
        >>> r = rank(matrix)
        >>>
        4

    See Also:
    _rank
    """
    _, _, _, rank_ = rref_complete(matrix)
    return rank_


def _rank(matrix: np.ndarray) -> int:
    """Computes the rank of the GF(2) matrix

    Args:
        matrix: Input GF(2) matrix

    Returns:
        Rank of input matrix

    Examples:
        >>> matrix = numpy.array([[1,0,0,1,0,0,1,0],
                                  [0,1,1,1,0,0,0,1],
                                  [1,1,1,0,1,0,0,0],
                                  [1,0,0,1,0,1,0,1]], dtype=np.bool_)
        >>> r = _rank(matrix)
        >>>
        4
    See Also:
    rank
    """
    _, _, _, rank_ = _rref_complete(matrix)
    return rank_


def rref_complete(matrix: np.ndarray) -> Tuple[List[int], np.ndarray, np.ndarray, int]:
    """Computes the Row Reduced Echelon Form for a GF(2) matrix as well
    as pivots, transformation matrix and rank.

    Args:
        matrix : Input GF(2) matrix.

    Returns:
        heads
            A 0-1 list with value of position k being 1 if the kth column is a pivot
        rref_mat
            Row Reduced Echelon Form of input matrix
        transform_mat
            Transform used to tranform input matrix into RREF form
        rank
            rank of input matrix

    Raises:
        QiskitError: Not a suitable matrix input")
        QiskitError: Not a two dimensional matrix")

    Examples:
        >>> matrix = numpy.array([[1,0,0,1,0,0,1,0],
                                  [0,1,1,1,0,0,0,1],
                                  [1,1,1,0,1,0,0,0],
                                  [1,0,0,1,0,1,0,1]], dtype=np.bool_)
        >>> heads, rref_mat, transform_mat, rank_ = rref_complete(matrix)
        >>> heads
        [1, 1, 0, 0, 1, 1, 0, 0]
        >>> rref_mat.astype(int)
        array([[1, 0, 0, 1, 0, 0, 1, 0],
               [0, 1, 1, 1, 0, 0, 0, 1],
               [0, 0, 0, 0, 1, 0, 1, 1],
               [0, 0, 0, 0, 0, 1, 1, 1]])
        >>> transform_mat.astype(int)
        array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [1, 1, 1, 0],
               [1, 0, 0, 1]])
        >>> rank_
        4

    See Also:
    _rref_complete, rref, _rref
    """
    matrix = np.array(matrix, dtype=bool)
    if not matrix.shape:  # matrix.shape == ():
        raise QiskitError("Not a suitable matrix input")
    matrix = np.atleast_2d(matrix)
    if not matrix.ndim == 2:
        raise QiskitError("Not a 2 dimensional matrix")

    return _rref_complete(matrix)


def _rref_complete(matrix) -> Tuple[List[int], np.ndarray, np.ndarray, int]:
    """Computes the Row Reduced Echelon Form for a GF(2) matrix as well
    as pivots, transformation matrix and rank.

    Args:
        matrix : Input GF(2) matrix.

    Returns:
        heads
            A 0-1 list with value of position k being 1 if the kth column is a pivot
        rref_mat
            Row Reduced Echelon Form of input matrix
        transform_mat
            Transform used to tranform input matrix into RREF form
        rank
            rank of input matrix

    Examples:
        >>> matrix = numpy.array([[1,0,0,1,0,0,1,0],
                                  [0,1,1,1,0,0,0,1],
                                  [1,1,1,0,1,0,0,0],
                                  [1,0,0,1,0,1,0,1]], dtype=np.bool_)
        >>> heads, rref_mat, transform_mat, rank_ = _rref_complete(matrix)
        >>> heads
        [1, 1, 0, 0, 1, 1, 0, 0]
        >>> rref_mat.astype(int)
        array([[1, 0, 0, 1, 0, 0, 1, 0],
               [0, 1, 1, 1, 0, 0, 0, 1],
               [0, 0, 0, 0, 1, 0, 1, 1],
               [0, 0, 0, 0, 0, 1, 1, 1]])
        >>> transform_mat.astype(int)
        array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [1, 1, 1, 0],
               [1, 0, 0, 1]])
        >>> rank_
        4

    See Also:
    rref_complete, rref, _rref

    Notes:
    The rref implemented here is not parallel and is for dense matrices. The current
    implementation will work for relatively small matrices but as matrices grow, as
    in LDPC code check matrices, a much better implementation should be used. In
    general, a series of methods should be inplemented to handle different needs.
    As a start, the heirachy of algorithms found in arXiv:1111.6549 should be
    implemented.
    """

    matrix = augment_mat(matrix, "right")

    nrows = matrix.shape[0]
    ncols = matrix.shape[1]
    hncols = ncols - nrows

    co_row = []
    heads = hncols * [-1]

    # Clear entries below pivots
    for i in range(nrows):
        row = matrix[i]

        for j in range(hncols):
            if heads[j] != -1:
                if row[j]:
                    row = row ^ co_row[heads[j]]

        # Find next pivot
        j = np.nonzero(row)[0][0]
        if j < hncols:
            co_row.append(row)
            heads[j] = len(co_row) - 1

    # Clear entries above pivots
    for j in range(hncols - 1, -1, -1):
        if heads[j] != -1:
            diff = [item for item in range(heads[j]) if item not in heads[j + 1 :]]
            for i in diff:
                row = co_row[i]
                if row[j]:
                    co_row[i] = co_row[i] ^ co_row[heads[j]]

    pivots = [i for i in heads if i != -1]
    co_row = [co_row[i] for i in pivots]
    co_row = np.array(co_row, dtype=bool)

    rref_mat = co_row[:, :hncols]
    transform_mat = co_row[:, hncols:]
    heads = [int(bool(x + 1)) for x in heads]
    rank_ = sum(heads)

    return heads, rref_mat, transform_mat, rank_


def istack(mat: np.ndarray, size: int, interleave: bool = False) -> np.ndarray:
    """Vertically stacks array of copies of vectors with or with interleaving.

    Args:
        mat: array of vectors to stack or interleave stack
        size: Number of copies to stack or interleave
        interleave (Oprional): Interleave copies if True, not if False. Default is False

    mat = [v_1
            v_2
            ...
            v_k]

    istack(mat, r, interleave=False) gives r vertically stacked copies of array with no iterleaving

    output = [v_1
              v_2
              ...
              v_k
              ... size times
              v_1
              v_2
              ...
              v_k]

    istack(mat, r, interleave=True) gives r vertically stacked copies of array with with iterleaving

    output = [v_1
              v_1
              ... size copies
              v_1
              ...
              v_k
              v_k
              ... size copies
              v_k]

    Returns:
        mat: vertically stacked array of size copies of vectors from input
    """
    if size == 1:
        return mat
    if interleave:
        return np.hstack(size * [mat]).reshape((size * len(mat),) + mat.shape[1:])
    return np.vstack(size * [mat]).reshape((size * len(mat),) + mat.shape[1:])


# ---------------------------------------------------------------
# Span-preserving row operations and Howell matrix form
# This code is adapted from XPFpackage: https://github.com/m-webster/XPFpackage, originally developed by Mark Webster. The original code is licensed under the GNU General Public License v3.0 and Mark Webster has given permission to use the code under the Apache License v2.0.


def do_row_op(mat: np.ndarray, row_op: Tuple[str, List[int], List[int]], N: int) -> np.ndarray:
    """Performs span-preserving row operations on a matrix in the ring Z/NZ. These include:

    1. Swap two rows mat[i] and mat[j] (_swap_rows);
    2. Multiply a row mat[i] by a scalar c (valid only when c is a unit), i.e., mat[i] = c * mat[i] (_multiply_unit);
    3. Add a multiple of one row to another, i.e., mat[i] = mat[i] + c * mat[j] (_add_rows);
    4. Append the product of a row by a scalar c to the end of matrix (used when c is a zero divisor), i.e., mat = mat.append(c * mat[i]) (_append_row);
    5. Update two rows by multiplying by a full-rank 2x2 matrix, i.e., mat[i] = a * mat[i] + b * mat[j], mat[j] = c * mat[i] + d * mat[j] (_update_rows).
    
    Args:
        mat: input matrix
        row_op: tuple (op, rows, coeff), where op is the operation to be performed, rows is a list of row indices, and coeff is a list of coefficients for the operation (empty list if op is 'swap')
        N: modulus

    Returns:
        matrix after performing the row operation

    Raises:
        QiskitError: Modulus N must be a positive integer.
        QiskitError: Input matrix must be a 2D array.
        QiskitError: Row operation must be a valid operation ("swap", "unit", "add", "append", or "update").
        QiskitError: Row operation must involve valid row indices for the input matrix.
        QiskitError: Swap operation must involve two rows.
        QiskitError: Unit operation must involve one row and one coefficient.
        QiskitError: Add operation must involve two rows and one coefficient.
        QiskitError: Append operation must involve one row and one coefficient.
        QiskitError: Update operation must involve two rows and one 2x2 matrix (four coefficients).

    Examples:
        >>> mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> do_row_op(mat, ('swap', [0, 2], []), 4)
        array([[7, 8, 9],
               [4, 5, 6],
               [1, 2, 3]])

        >>> do_row_op(mat, ('unit', [0], [2]), 4)
        array([[2, 0, 2],
               [4, 5, 6],
               [7, 8, 9]])

    """
    if not isinstance(N, (int, np.integer)) or N <= 0:
        raise QiskitError("Modulus N must be a positive integer.")
    mat = np.array(mat, dtype=int)
    if not mat.ndim == 2:
        raise QiskitError("Input matrix must be a 2D array")
    
    # get number of rows
    nrows = mat.shape[0]
    op, rows, coeff = row_op

    if op not in ["swap", "unit", "add", "append", "update"]:
        raise QiskitError("Row operation must be a valid operation (\"swap\", \"unit\", \"add\", \"append\", or \"update\").")
    if any([row >= nrows for row in rows]):
        raise QiskitError("Row operation must involve valid row indices for the input matrix.")
    if op == 'swap':
        if len(rows) != 2:
            raise QiskitError("Swap operation must involve two rows.")
        mat = _swap_rows(mat, rows)
    elif op == 'unit':
        if len(rows) != 1 or len(coeff) != 1:
            raise QiskitError("Unit operation must involve one row and one coefficient.")
        mat = _multiply_unit(mat, rows[0], coeff[0], N)
    elif op == 'add':
        if len(rows) != 2 or len(coeff) != 1:
            raise QiskitError("Add operation must involve two rows and one coefficient.")
        mat = _add_rows(mat, rows, coeff[0], N)
    elif op == 'append':
        if len(rows) != 1 or len(coeff) != 1:
            raise QiskitError("Append operation must involve one row and one coefficient.")
        mat = _append_row(mat, rows[0], coeff[0], N)
    elif op == 'update':
        if len(rows) != 2 or len(coeff) != 4:
            raise QiskitError("Update operation must involve two rows and one 2x2 matrix (four coefficients).")
        mat = _update_rows(mat, rows, coeff, N)

    return mat


def _swap_rows(mat: np.ndarray, rows: Union[list, np.ndarray]) -> np.ndarray:
    """Swaps two rows of a matrix.

    Args:
        mat: input matrix
        rows: list of indices of the two rows to swap

    Returns:
        matrix with rows swapped

    Examples:
        >>> mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> _swap_rows(mat, [0, 2])
        array([[7, 8, 9],
               [4, 5, 6],
               [1, 2, 3]])
    """
    mat[[rows[0], rows[1]]] = mat[[rows[1], rows[0]]]
    return mat


def _multiply_unit(mat: np.ndarray, row: int, c: int, N: int) -> np.ndarray:
    """Multiplies a row of a matrix by a scalar c (valid only when c is a unit) in the ring Z/NZ.

    Args:
        mat: input matrix
        row: index of row to multiply
        c: scalar to multiply by
        N: modulus

    Returns:
        matrix with row multiplied by scalar

    Examples:
        >>> mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> _multiply_unit(mat, 0, 2, 4)
        array([[2, 0, 2],
               [4, 5, 6],
               [7, 8, 9]])
    """
    mat[row] = np.mod(c * mat[row], N)
    return mat


def _add_rows(mat: np.ndarray, rows: Union[list, np.ndarray], c: int, N: int) -> np.ndarray:
    """Adds a multiple of one row to another of a matrix in the ring Z/NZ.

    Args:
        mat: input matrix
        rows: list of indices of the two rows in action
        c: scalar to multiply by
        N: modulus

    Returns:
        matrix with rows added

    Examples:
        >>> mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> _add_rows(mat, [0, 2], 2, 4)
        array([[3, 2, 1],
               [4, 5, 6],
               [7, 8, 9]])
    """
    mat[rows[0]] = np.mod(mat[rows[0]] + c * mat[rows[1]], N)
    return mat


def _append_row(mat: np.ndarray, row: int, c: int, N: int) -> np.ndarray:
    """Appends the product of a row by a scalar c to the end of matrix (used when c is a zero divisor) in the ring Z/NZ.

    Args:
        mat: input matrix
        row: index of row to multiply
        c: scalar to multiply by
        N: modulus

    Returns:
        matrix with row multiplied by scalar and appended

    Examples:
        >>> mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> _append_row(mat, 0, 2, 4)
        array([[1, 2, 3],
               [4, 5, 6],
               [7, 8, 9],
               [2, 0, 2]])
    """
    mat = np.vstack((mat, np.mod(c * mat[row], N)))
    return mat


def _update_rows(mat: np.ndarray, rows: Union[list, np.ndarray], c: Union[list, np.ndarray], N: int) -> np.ndarray:
    """Updates two rows by multiplying by a full-rank 2x2 matrix in the ring Z/NZ.

    Args:
        mat: input matrix
        rows: list of indices of the two rows in action
        c: components of the 2x2 matrix to multiply by
        N: modulus

    Returns:
        matrix with rows updated

    Examples:
        >>> mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> _update_rows(mat, [0, 1], [2, 1, 3, 4], 4)
        array([[2, 1, 0],
               [3, 2, 1],
               [7, 8, 9]])
    """
    r1 = np.mod(c[0] * mat[rows[0]] + c[1] * mat[rows[1]], N)
    r2 = np.mod(c[2] * mat[rows[0]] + c[3] * mat[rows[1]], N)
    mat[rows[0]], mat[rows[1]] = r1, r2
    return mat


def howell(mat: np.ndarray, N: int) -> np.ndarray:
    """Computes the Howell form of a matrix in the ring Z/NZ.

    Args:
        mat: input matrix
        N: modulus

    Returns:
        Howell form of input matrix

    Examples:
        >>> mat = numpy.array([[8, 5, 5],
                               [0, 9, 8],
                               [0, 0, 10]])
        >>> N = 12
        >>> H = howell(mat, N)
        >>> H
        array([[4, 1, 0],
               [0, 3, 0],
               [0, 0, 1]])

    See Also:
    _howell, howell_complete, _howell_complete
    """
    H, _, _ = howell_complete(mat, N)
    return H


def _howell(mat: np.ndarray, N: int) -> np.ndarray:
    """Computes the Howell form of a matrix in the ring Z/NZ.

    Args:
        mat: input matrix
        N: modulus

    Returns:
        Howell form of input matrix

    Examples:
        >>> mat = numpy.array([[8, 5, 5],
                               [0, 9, 8],
                               [0, 0, 10]])
        >>> N = 12
        >>> H = _howell(mat, N)
        >>> H
        array([[4, 1, 0],
               [0, 3, 0],
               [0, 0, 1]])

    See Also:
    howell, howell_complete, _howell_complete
    """
    H, _, _ = _howell_complete(mat, N)
    return H


def howell_complete(mat: np.ndarray, N: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Computes the Howell form of a matrix in the ring Z/NZ, the corresponding transformation matrix, and the kernel.

    Args:
        mat: input matrix
        N: modulus

    Returns:
        H: Howell form of mat
        U: transformation matrix (U @ mat = H)
        K: kernel of mat (mat @ K = 0)

    Raises:
        QiskitError: Modulus N must be a positive integer
        QiskitError: Input matrix must be a 2D array

    Examples:
        >>> mat = numpy.array([[8, 5, 5],
                               [0, 9, 8],
                               [0, 0, 10]])
        >>> N = 12
        >>> H, U, K = howell(mat, N)
        >>> H
        array([[4, 1, 0],
               [0, 3, 0],
               [0, 0, 1]])
        >>> U
        array([[8, 1, 0],
               [0, 7, 4],
               [9, 3, 4]])
        >>> K
        array([[6, 6, 6],
               [0, 4, 4]])

    See Also: 
    _howell_complete, howell, _howell
    """
    if not N > 0 or not isinstance(N, (int, np.integer)):
        raise QiskitError("Modulus N must be a positive integer")
    mat = np.array(mat, dtype=int)
    if not mat.ndim == 2:
        raise QiskitError("Input matrix must be a 2D array")
    
    return _howell_complete(mat, N)


def _howell_complete(mat: np.ndarray, N: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Computes the Howell form of a matrix in the ring Z/NZ, the corresponding transformation matrix, and the kernel.

    Args:
        mat: input matrix
        N: modulus

    Returns:
        H: Howell form of mat
        U: transformation matrix (U @ mat = H)
        K: kernel of mat (K @ mat = 0)

    Examples:
        >>> mat = numpy.array([[8, 5, 5],
                               [0, 9, 8],
                               [0, 0, 10]])
        >>> N = 12
        >>> H, U, K = howell(mat, N)
        >>> H
        array([[4, 1, 0],
               [0, 3, 0],
               [0, 0, 1]])
        >>> U
        array([[8, 1, 0],
               [0, 7, 4],
               [9, 3, 4]])
        >>> K
        array([[6, 6, 3],
               [0, 4, 4]])

    See Also:
    howell_complete, howell, _howell
    """
    H = mat.copy()
    U = np.eye(H.shape[0], dtype=int)
    m, n = H.shape
    
    # set row index to 0
    r = 0
    # going through each column
    for c in range(n):
        # find j such that H[j, c] > 0
        j = r
        while j < m and H[j, c] == 0:
            j += 1
        if j < m:
            # found j: if j > r, swap rows r and j
            if j > r:
                H = _swap_rows(H, r, j)
                U = _swap_rows(U, r, j)

            # multiply row r by a unit to ensure that H[r, c] is a minimal representative
            x = unit(H[r, c], N)
            if x > 1:
                H = _multiply_unit(H, r, x, N)
                U = _multiply_unit(U, r, x, N)
            
            # eliminate entries in column c below row r
            for i in range(r + 1, m):
                if H[i, c] % N > 0:
                    _, s, t, u, v = gcd_ext(H[r, c], H[i, c], N)
                    H = _update_rows(H, [r, i], [s, t, u, v], N)
                    U = _update_rows(U, [r, i], [s, t, u, v], N)
            
            # ensure entries in column c above row r are less than H[r, c]
            b = H[r, c]
            for i in range(r):
                if H[i, c] >= b:
                    x = quo(H[i, c], b, N)
                    H = _add_rows(H, [i, r], -x, N)
                    U = _add_rows(U, [i, r], -x, N)
            
            # if b = H[r, c] is a zero divisor, find the annihilator x that eliminates H[r, c] and append a new row x * H[r]
            x = ann(b, N)
            if x > 0:
                H = _append_row(H, r, x, N)
                U = _append_row(U, r, x, N)
                m = len(H)
            r += 1
        
    # remove rows of zeros
    H = H[H.any(axis=1)]

    # compute the transformation matrix and kernel
    k = len(H)
    K = U[k:, :]
    K = K[K.any(axis=1)]
    U = U[:k, :]

    return H, U, K