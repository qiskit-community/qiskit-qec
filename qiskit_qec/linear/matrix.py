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

from typing import List, Tuple
import numpy as np
from qiskit import QiskitError


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
