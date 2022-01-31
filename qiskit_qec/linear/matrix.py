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
from typing import List
from qiskit import QiskitError

def create_lambda_matrix(n: int):
    """Create the boolean Lambda Matrix
    
    [0 I_(n/2)]
    [I_(n/2) 0]

    Args:
        n (int): dimension of I_n, n must be even and positive integer

    Returns:
        numpy.ndarray: boolean Lambda matrix for the symplectic product
    """
    assert n>0, QiskitError("n must be positive")
    assert isinstance(n,(int,np.integer)), QiskitError("n must be an integer")
    assert not n%2, QiskitError("n must be an even integer")

    return _create_lambda_matrix(n)


def _create_lambda_matrix(n: int):
    """Create the boolean Lambda Matrix
    
    [0 I_(n/2)]
    [I_(n/2) 0]

    Args:
        n (int): dimension of I_n

    Returns:
        numpy.ndarray: boolean Lambda matrix for the symplectic product
    """
    half_n = n>>1
    lambda_matrix = np.eye(n, k=half_n, dtype=bool)
    lambda_matrix[half_n:n, 0:half_n] = np.identity(half_n, dtype=bool)

    return lambda_matrix

def augment_mat(matrix, pos="left"):
    """Augment a given matrix with the identify matrix. The position of the 
    identity matrix is given by the optional argument <pos>

    left    ->  [M | I]

    right   ->  [I | M]

    top     ->  [I]
                [M]

    bottom  ->  [M]
                [I]

    Args:
        matrix (array): Matrix to agument
        pos (str, optional): Position of augmentation. Defaults to "left".

    Returns:
        [np.ndarray]: Augmented matrix with Identity
    """
    matrix = np.asarray(matrix)
    assert matrix.ndim == 2, QiskitError("Input array must have two dimensions")
    assert pos in ["top", "bottom", "right", "left"], QiskitError("Unknown position")
    return _augment_mat(matrix, pos)

def _augment_mat(matrix, pos):
    """Augment a given matrix with the identify matrix. The position of the 
    identity matrix is given by the argument <pos>

    left    ->  [M | I]

    right   ->  [I | M]

    top     ->  [I]
                [M]

    bottom  ->  [M]
                [I]

    Args:
        matrix (np.array): Matrix to agument
        pos (str): Position of augmentation.

    Returns:
        [np.ndarray]: Augmented matrix with Identity
    """
    
    if pos == "top":
        _id = np.identity(matrix.shape[1], dtype = matrix.dtype)
        return np.vstack((_id, matrix))
    elif pos == "bottom":
        _id = np.identity(matrix.shape[1], dtype = matrix.dtype)
        return np.vstack((matrix, _id))
    elif pos == "left":
        _id = np.identity(matrix.shape[0], dtype = matrix.dtype)
        return np.hstack((_id, matrix))
    elif pos == "right":
        _id = np.identity(matrix.shape[0], dtype = matrix.dtype)
        return np.hstack((matrix, _id))    

def rref(matrix):
    """Compute the Row Reduced Echelon Form for a boolean matrix with checks and matrix conversions.

    Args:
        matrix (2d matrix): Matrix of some form. Will be converted to an numpy.ndarray
        of bools

    Returns:
        rref_mat (numpy.ndarray) : Row Reduced Echelon Form of input matrix
    """
    heads, rref_mat, transform_mat, rank = rref_complete(matrix)
    return rref_mat

def _rref(matrix):
    """Compute the Row Reduced Echelon Form for a boolean matrix

    Args:
        matrix (2d matrix): numpy.ndarry of a boolean matrix

    Returns:
        rref_mat (numpy.ndarray) : Row Reduced Echelon Form of input matrix
    """
    heads, rref_mat, transform_mat, rank = _rref_complete(matrix)
    return rref_mat

def rref_complete(matrix):
    """Compute the Row Reduced Echelon Form for a boolean matrix and associated 
    results with checks and matrix conversions

    Args:
        matrix (2d matrix): Matrix of some form. Will be converted to an numpy.ndarray
        of bools

    Returns:
        heads (list): A 0-1 list with value of position k being 1 if the kth column is a pivot
        rref_mat (numpy.ndarray) : Row Reduced Echelon Form of input matrix
        transform_mat : Transform used to tranform input matrix into RREF form
        rank :  rank of input matrix
    """
    matrix = np.array(matrix, dtype=bool)
    assert matrix.shape != (), QiskitError("Not a suitable matrix input")
    matrix = np.atleast_2d(matrix)
    assert matrix.ndim == 2, QiskitError("Not a 2 dimensional matrix")

    return _rref_complete(matrix)

def _rref_complete(matrix):
    """Compute the Row Reduced Echelon Form for a boolean matrix and associated 
    results.

    Args:
        matrix (2d matrix): numpy.ndarry of a boolean matrix

    Returns:
        heads (list): A 0-1 list with value of position k being 1 if the kth column is a pivot
        rref_mat (numpy.ndarray) : Row Reduced Echelon Form of input matrix
        transform_mat : Transform used to tranform input matrix into RREF form
        rank :  rank of input matrix
    """
    # The rref implemented here is not parallel and is not sophisticated. Will work just
    # fine for relatively small matrices but as matrices grow, as in LDPC code check matrices,
    # we should have a much better implementation. Specifically we should implement
    # a general rref method that calls the appropriate methods depending on the 
    # characteristics of the matrix. I would start by implementing the algorithms in the
    # following paper: arXiv:1111.6549 from 2011 for matrices over GF(2). For matrices over
    # GF(q) and alike the approach below can be easily modified as a start.

    matrix = augment_mat(matrix, "right")

    nrows = matrix.shape[0]
    ncols = matrix.shape[1]
    hncols = ncols-nrows

    co_row = []
    heads = hncols*[-1]

    # Clear entries below pivots
    for i in range(nrows):
        row = matrix[i]

        for j in range(hncols):
            if heads[j] != -1:
                if row[j] == True:
                    row = row ^ co_row[heads[j]]

        # Find next pivot
        j = np.nonzero(row)[0][0]
        if j < hncols:
            co_row.append(row)
            heads[j] = len(co_row)-1

    # Clear entries above pivots
    for j in range(hncols-1,-1,-1):
        if heads[j] != -1:
            diff = [item for item in range(heads[j]) if item not in heads[j+1:]]
            for i in diff:
                row = co_row[i]
                if row[j] == True:
                    co_row[i] = co_row[i] ^ co_row[heads[j]]
                    
    pivots = [i for i in heads if i != -1]
    co_row = [co_row[i] for i in pivots]
    co_row = np.array(co_row, dtype=bool)
    rank = sum(pivots)
    
    rref_mat = co_row[:,:hncols]
    transform_mat = co_row[:,hncols:]
    heads = [int(bool(x+1)) for x in heads]
   
    
    return heads, rref_mat, transform_mat, rank