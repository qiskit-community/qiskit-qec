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

from time import perf_counter
from typing import List, Tuple
import itertools
import numpy as np
from qiskit import QiskitError
from numba import njit, jit
from qiskit_qec.analysis.extensions import _c_solve


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

def ker2(a):
    """ calculates the kernel of the binary matrix 'a' over the field GF(2). Adapted from code from S. Bravyi.
    Returns a basis for the ker2(a) as rows of a 2d numpy.ndarray. """
    m,n = a.shape
    ker = np.identity(n,dtype=int)

    for i in range(m):
        y = np.dot(a[i], ker) % 2 # multiplication of current row with all columns of ker
        good = ker[:,y==0] # columns of ker that are in the kernel of a[i,:] (and thus in the kernel of a[:i+1,:])
        bad = ker[:, y==1] # colums of ker that are in kernel of a[:i,:] but not in kernel of a[i,:]
        if bad.shape[1]>0: # in case there are enough columns not in the kernel
            new_good = (bad[:,:-1] + bad[:,1:]) % 2 # by construction all of these will be in kernel of a[i,:], independent and complete
            ker = np.concatenate((good, new_good), axis=1) # new basis for kernel of a[:i+1,:]
    # now columns of ker span the binary null-space of a
    return np.transpose(ker)

def is_binary(arr: np.array):
    """ Check if a numpy.array is binary"""
    return np.all((arr.astype(int) == 1) + (arr.astype(int) == 0))

class LinAlgError(ValueError):
    pass

def solve2_numba(A: np.array, b: np.array):
    if not is_binary(A) or not is_binary(b):
        raise ValueError('A and b must be binary arrays')
    
    Ap, bp = jitable(A, b)
    if np.array_equal(Ap, np.empty((0, A.shape[1]), dtype=np.bool_)):
        raise LinAlgError('System has no solution')
    
    t0 = perf_counter()
    error, t_part, nullity = _back_substitution_weight_opt(Ap,bp) #, A.astype(int), b.astype(int)
    t_opt = perf_counter() - t0 - t_part

    return error, (t_part, nullity, t_opt)
    

@jit(nopython=True)
def check_all_zero_rows(A):
    """Manually checks each row for all zeros, compatible with Numba."""
    m = A.shape[0]
    all_zero_rows = np.zeros(m, dtype=np.bool_)
    for i in range(m):
        all_zero_rows[i] = np.all(A[i] == 0)
    return all_zero_rows

@jit(nopython=True)
def jitable(A, b):
    A = np.copy(A).astype(np.bool_)
    b = np.copy(b).astype(np.bool_)
    m, n = A.shape
    if len(b) != m:
        return np.empty((0, A.shape[1]), dtype=np.bool_), np.empty(0, dtype=np.bool_)
    # We can check this in the beginning. Is vital if we never actually do a gaussian elimination, because then we never check
    all_zero_rows = check_all_zero_rows(A) # Find all all-zero rows of A (no degree of freedom)
    safe = np.all(b[all_zero_rows] == 0) # If we have the trivial equation 0=0 for all of them we are safe
    if not safe: # otherwise 0=1 for at least one row, abort as no solution exists
        return np.empty((0, A.shape[1]), dtype=np.bool_), np.empty(0, dtype=np.bool_)
    
    g = 0 # how many times gaussian elimination was actually done
    for i in range(n): # go through all columns of A with increment variable i
        idxs = np.where(A[:, i])[0] # at current column find all rows that have a 1
        potentials = idxs[idxs >= g]
        if len(potentials) == 0:
            continue
        idx = idxs[idxs >= g][0]

        for target in idxs: # Perform Gaussian elimination with the row found above targeting all other rows that have a 1 at column i
            if target == idx:
                continue
            A[target] ^= A[idx]
            b[target] ^= b[idx]
        
        all_zero_rows = check_all_zero_rows(A) # Find all all-zero rows of A (no degree of freedom)
        safe = np.all(b[all_zero_rows] == 0) # If we have the trivial equation 0=0 for all of them we are safe
        if not safe: # otherwise 0=1 for at least one row, abort as no solution exists
            return np.empty((0, A.shape[1]), dtype=np.bool_), np.empty(0, dtype=np.bool_)

        # A[[idx,g]] = A[[g,idx]] # Swap the row that was used for elimination with the one at that has index at the current step i
        # b[[idx,g]] = b[[g,idx]] # Swap the row that was used for elimination with the one at that has index at the current step i
        # Manual row swap
        # Direct assignment for swapping
        if g != idx:
            A[g, :], A[idx, :] = A[idx, :].copy(), A[g, :].copy()
            b[g], b[idx] = b[idx], b[g]

        g += 1 # increment g
    return A,b

def solve2_python(A: np.array, b: np.array):
    """
    Solves the system of equations Ax = b mod 2 for binary matrices and vectors b.
    Will raise an exception if no solution exists.
    Return a solution x and the reduced system of equations (A,b),
    i.e. A will be in reduced row echelon form (but with zero rows 
    remaining to preserve shape) and b corresponding.
    NOTE: this could also be achieved with rref_complete and intelligently using transform_mat
    to transform b and to see if there are solutions. Currently rref_complete is faster than this for smaller matrices.
    This implementation is faster for matrices starting with dimensions of several hundred.
    """

    if not is_binary(A) or not is_binary(b):
        raise ValueError('A and b must be binary arrays')

    A = np.copy(A).astype(bool)
    b = np.copy(b).astype(bool)
    m, n = A.shape
    if len(b) != m:
        raise ValueError('Non compatible shapes')
    
    # We can check this in the beginning. Is vital if we never actually do a gaussian elimination, because then we never check
    all_zero_rows = np.all(A == 0, axis=1) # Find all all-zero rows of A (no degree of freedom)
    safe = np.all(b[all_zero_rows] == 0) # If we have the trivial equation 0=0 for all of them we are safe
    if not safe: # otherwise 0=1 for at least one row, abort as no solution exists
        raise LinAlgError('System has no solution')
    
    g = 0 # how many times gaussian elimination was actually done
    for i in range(n): # go through all columns of A with increment variable i
        idxs = np.where(A[:, i])[0] # at current column find all rows that have a 1
        try:
            idx = idxs[idxs >= g][0] # find the first row that is at g or higher. This will not have any other 1's before
        except IndexError:
            continue
        for target in idxs: # Perform Gaussian elimination with the row found above targeting all other rows that have a 1 at column i
            if target == idx:
                continue
            A[target] ^= A[idx]
            b[target] ^= b[idx]
        
        all_zero_rows = np.all(A == 0, axis=1) # Find all all-zero rows of A (no degree of freedom)
        safe = np.all(b[all_zero_rows] == 0) # If we have the trivial equation 0=0 for all of them we are safe
        if not safe: # otherwise 0=1 for at least one row, abort as no solution exists
            raise LinAlgError('System has no solution')

        A[[idx,g]] = A[[g,idx]] # Swap the row that was used for elimination with the one at that has index at the current step i
        b[[idx,g]] = b[[g,idx]] # Swap the row that was used for elimination with the one at that has index at the current step i

        g += 1 # increment g

    t0 = perf_counter()
    error, t_part, nullity = _back_substitution_weight_opt(A,b) #, A.astype(int), b.astype(int)
    t_opt = perf_counter() - t0 - t_part

    return error, (t_part, nullity, t_opt)

def solve2_cpp(A: np.array, b: np.array):
    """
    Solves the system of equations Ax = b mod 2 for binary matrices and vectors b.
    Will raise an exception if no solution exists.
    Return a solution x and the reduced system of equations (A,b),
    i.e. A will be in reduced row echelon form (but with zero rows 
    remaining to preserve shape) and b corresponding.
    NOTE: this could also be achieved with rref_complete and intelligently using transform_mat
    to transform b and to see if there are solutions. Currently rref_complete is faster than this for smaller matrices.
    This implementation is faster for matrices starting with dimensions of several hundred.
    """

    if not is_binary(A) or not is_binary(b):
        raise ValueError('A and b must be binary arrays')
    if A.shape[1] == 0 and np.any(b):
        raise LinAlgError('System has no solution')
    solvable, A, b, _ = _c_solve(A,b)
    A = np.array(A)
    b = np.array(b)
    if not solvable:
        raise LinAlgError('System has no solution')

    t0 = perf_counter()
    error, t_part, nullity = _back_substitution_weight_opt(A,b) #, A.astype(int), b.astype(int)
    t_opt = perf_counter() - t0 - t_part

    return error, (t_part, nullity, t_opt)

def _back_substitution(A, b):
    """
    Backsubstitution step for solving of linear system of equations mod 2.
    Does NOT give minimum weight solution, but just an arbitrary solution (mostly for testing).
    Input: A and b after gaussian elimination step (A is upper triangular and no 0=1 rows exist)
    """
    A = np.copy(A)
    b = np.copy(b)
    m, n = A.shape
    x = np.nan*np.zeros(n)
    r = min(m,n)
    for i in range(r):
        ones = np.where(A[r-1-i])[0] # find all 1 entries in this row
        x[ones[:-1]] = 0 # in x set all these rows to 0, except the last
        if len(ones) > 0:
            x[ones[-1]] = b[r-1-i] # the last of these rows we set to the value of b at this row. Now we have one solution for this row, not dependant on nans
            # if ones is emtpy that means all zero row, since we assume there is a solution is this step, we can just continue
            
        aint_nans = ~np.isnan(x)
        b = (b + A[:,aint_nans] @ x[aint_nans]) % 2
        A[:, aint_nans] = 0

    x[np.isnan(x)] = 0 # set remaining nans (degrees of freedom) to zero
    return x.astype(int)

def _back_substitution_weight_opt(A, b):
    """
    Backsubstitution step for solving of linear system of equations mod 2.
    Finds ALL solutions of the LSE and returns the minimum weight one.
    Input: A and b after gaussian elimination step (A is upper triangular and no 0=1 rows exist)
    """
    t0 = perf_counter()
    xs = _back_substitution(A, b) # an arbitrary inhomogeneous solution
    t_part = perf_counter() - t0
    ker = ker2(A) # basis of nullspace of A 
    # all solutions of Ax = b can be written as xs + ker

    best = xs
    min_weight = xs.sum()
    if ker.shape[0] > 7:
        return best, t_part, ker.shape[0]
    for sel in itertools.product([False,True], repeat=ker.shape[0]):
        x = (xs + ker[list(sel)].sum(axis=0)) % 2
        if x.sum() < min_weight:
            best = x
            min_weight = x.sum()
    
    return best, t_part, ker.shape[0]

