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

from collections import deque
from typing import List, Any, Tuple
from typing import Union, Optional

import numpy as np
from numpy.typing import ArrayLike

from qiskit import QiskitError
from qiskit_qec.linear import matrix as mt


def all_commute(matrix: ArrayLike) -> bool:
    r"""Determines if each possible pair of different rows of the
    GF(2) symplectic matrix have zero symplectic product. If the rows represent
    Pauli operators then the this method deterimes if the operators
    defined by the matrix generate an abelian subgroup.

    Args:
        matrix : Input GF(2) symplectic matrix

    Returns:
        out: True if operators mutually commute - have zero symplectic product

    Examples:
        >>> matrix = numpy.array([[1,0,0,1,0,0,1,0],
                                  [0,1,1,1,0,0,0,1],
                                  [1,1,1,0,1,0,0,0],
                                  [1,0,0,1,0,1,0,1]], dtype=numpy.bool_)
        >>> all_commute(matrix)
        False

        >>> matrix = np.array([[1,0,0,0,0,0,0,0],
                               [0,1,0,0,0,0,0,0],
                               [0,0,1,0,0,0,0,0],
                               [1,0,0,0,0,0,0,1]], dtype=numpy.bool_)
        >>> all_commute(matrix)
        True
    """
    test_mat = symplectic_product(matrix, matrix)
    return not test_mat.any()


# ---------------------------------------------------------------


def symplectic_product(mat1: ArrayLike, mat2: ArrayLike) -> int:
    r"""Returns the symplectic product of two GF(2) symplectic matrices.

    Let math:'A', math:'B' be two GF(2) symplectic matrices of width math:'2m',
    then the symplectic product is defined as

    .. math::
        (A,B) = A \cdot \Lambda_n \cdot B^T.

    Args:
        mat1, mat2: Input GF(2) symplectic matrixes

    Returns:
        int: Symplectic product of mat1 and mat2

    Raises:
        QiskitError: Input matrices/vectors must be GF(2) symplectic matrices/vectors
        QiskitError: Input matrices must have the same number of dimensions
        QiskitError: Input matrices must be 1 or 2 dimensional

    Examples:
        >>> mat1 = numpy.array([[1,0,1,1],[0,1,1,0]], dtype=np.bool_)
        >>> mat2 = numpy.array([[0,1,1,1],[0,0,1,0]], dtype=np.bool_)
        >>> p = symplectic_product(mat1, mat2)
        >>> p.astype(int)
        array([[0, 1], [1, 0]]

    See Also:
    _symplectic_product_vv, _symplectic_product_dense
    """
    mat1_np_array = np.array(mat1, dtype=np.int8)
    mat2_np_array = np.array(mat2, dtype=np.int8)

    if not is_symplectic_form(mat1) or not is_symplectic_form(mat2):
        raise QiskitError(f"Input matrices/vectors must be GF(2) symplectic matrices/vectors")

    if not mat1_np_array.ndim == mat2_np_array.ndim:
        QiskitError(
            f"Input matrices must have the \
        same dimensions: {mat1_np_array.ndim} is \
        not equal to {mat2_np_array.ndim}"
        )

    if mat1_np_array.ndim == 1:
        if not mat1_np_array.shape[0] == mat2_np_array.shape[0]:
            QiskitError(
                f"Input vectors must have the same \
            dimensions: {mat1_np_array.shape[0]} not equal \
            to {mat2_np_array.shape[0]}"
            )

        return _symplectic_product_vv(mat1, mat2, mat1_np_array.shape[0] >> 1)

    elif mat1_np_array.ndim == 2:
        return _symplectic_product_dense(mat1_np_array, mat2_np_array)

    else:
        raise QiskitError(
            f"Input matrices must be 1 or 2 dimensional:\
            {mat1_np_array.ndim}, {mat2_np_array.ndim}"
        )


def _symplectic_product_vv(vec1: ArrayLike, vec2: ArrayLike, n: int) -> int:
    r"""Finds the sympletic product or two GF(2) symplectic vectors of
    length 2n: vec1 . Lambda . vec2^T where

    lambda = [0 I]
             [I 0]

    Args:
        vec1, vec2: Input GF(2) symplectic vectors
        n: Input size, half of the symplectic vector length

    Returns:
        out: Symplectic product of vec1 and vec2

    Examples:
        >>> a = np.array([1,0,0,0,1,1,1,0,1,0])
        >>> b = np.array([1,1,1,1,0,0,1,0,1,1])
        >>> _symplectic_product_vv(a, b, 5)
        0

    See Also:
    symplectic_product, _symplectic_product_dense

    Notes:
    If changing the method please make sure that the new method
    is faster. Note that this method is faster if the vectors are
    numpy arrays with dtype=int8
    """
    r = 0
    for i in range(n):
        r += vec1[i] * vec2[n + i] + vec1[n + i] * vec2[i]
    return r % 2


def _symplectic_product_dense(mat1: np.ndarray, mat2: np.ndarray) -> Union[int, np.ndarray]:
    r"""Returns the symplectic product of two GF(2) symplectic matrices.

    Let math:'A', math:'B' be two GF(2) symplectic matrices of width math:'2m',
    then the symplectic product is defined as

    .. math::
        (A,B) = A \cdot \Lambda_n \cdot B^T.

    Args:
        mat1, mat2: Input GF(2) symplectic matrixes

    Returns:
        int: Symplectic product of mat1 and mat2

    Examples:
        >>> mat1 = numpy.array([[1,0,1,1],[0,1,1,0]], dtype=np.bool_)
        >>> mat2 = numpy.array([[0,1,1,1],[0,0,1,0]], dtype=np.bool_)
        >>> p = _symplectic_product_dense(mat1, mat2)
        >>> p.astype(int)
        array([[0, 1], [1, 0]]

    See Also:
    _symplectic_product_vv, _symplectic_product_numpy
    """
    m1, m2 = np.hsplit(mat1, 2)
    result = np.hstack((m2, m1)).dot(mat2.transpose()) % 2
    if result.size == 1:
        return int(result.item())
    else:
        return result


# ---------------------------------------------------------------


def make_commute_hyper(
    a: ArrayLike,
    x: ArrayLike,
    z: ArrayLike,
    arange: Optional[ArrayLike] = None,
    xrange: Optional[ArrayLike] = None,
    zrange: Optional[ArrayLike] = None,
) -> np.ndarray:
    r"""Makes an element(s) commute with hyperbolic pair(s)

    Let a = [a_0,...,a_(k-1)] where a_i are GF(2) symplectic vectors. Let
    x = [x_0,...,x_(l-1)] and z =[z_0,...,z_(l-1)] where x_i and z_i are
    GF(2) symplectic vectors such that (x_i,z_i) are hyerbolic pairs from
    the hyperbolic a basis <x_0,...,x_(l-1),z_0,...,x_(l-1)>. It is assumed
    that {a_0,...,a_(k-1),x_0,...,x_(l-1),z_0,...,x_(l-1)} is an independent
    set of symplectic vectors.

    This method returns a set of vectors b = [b_0,...,b_(k-1)] such that

    1) b_0, ..., b_(k-1) each have zero symplectic product with each of the
    hyperbolic vectors x_0,...,x_(l-1),z_0,...,x_(l-1)
    2) span(b_i, x_j,z_j)  = span(a_i, x_j,z_j)
    for i=0,..,k-1 and j=0,...,l-1

    If the symplectic vectors are considered as Pauli operators then the
    method returns a set of operators [op_(b_0),...,op_(b_(k-1))] such that

    1) op_(b_0), ..., op_(b_(k-1)) each commute with the each of the
    hyperbolic operators op_(x_0),...,op_(x_(l-1)),op_(z_0),...,op_(x_(l-1))
    2) <op_(b_i), op_(x_j),op_(z_j)> = <op_(a_i), op_(x_j),op_(z_j)>
    for each i=0,...,k-1 and j=0,...,l-1

    Args:
        a: Input GF(2) symplectic vectors
        x,z: GF(2) hyperbolic pair vectors
        arange (optional): range of indices from a to make commute. Defaults to None.
        xrange (optional): range of indices from x to use. Defaults to None.
        zrange (optional): range of indices from z to use. Defaults to None.

    Raises:
        QiskitError: Input matrices/vectors must bf GF(2) symplectic matrices/vectors
        QiskitError: Input range is not iterable")
        QiskitError: Input matrices/vectors must have the same number of columns/length

    Returns:
        out: GF(2) symplectic vectors that commute with the given hyperbolic pairs

    Examples:
        >>> a = np.array([1,1,1,0,0,0],dtype=np.bool_)
        >>> x = np.array([0,0,1,0,0,0],dtype=np.bool_)
        >>> z = np.array([0,0,0,0,0,1],dtype=np.bool_)
        >>> a = make_commute_hyper(a, x, z)
        >>> a.astype(int)
        array([1, 1, 0, 0, 0, 0])

        >>> a = np.array([1,1,1,0,0,0,0,0], dtype=np.bool_)
        >>> x = np.array([[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]], dtype=np.bool_)
        >>> z = np.array([[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0]], dtype=np.bool_)
        >>> xrange = [0,1]
        >>> zrange = [0,1]
        >>> a = make_commute_hyper(a, x, z, xrange = xrange, zrange=zrange)
        >>> a.astype(int)
        array([1, 0, 0, 0, 0, 0, 0, 0])

        >>> a = np.array([[1,1,1,0,0,0,0,0],[0,1,1,0,0,0,0,0]], dtype=np.bool_) # X1X2X3, X2X3
        >>> x = np.array([0,1,0,0,0,0,0,0], dtype=np.bool_) # X2
        >>> z = np.array([0,0,0,0,0,1,0,0], dtype=np.bool_) # Z2
        >>> arange = [0,1]
        >>> a = make_commute_hyper(a, x, z, arange)
        >>> a.astype(int)
        array([[1, 0, 1, 0, 0, 0, 0, 0],
               [0, 0, 1, 0, 0, 0, 0, 0]])

        >>> a = np.array([[1,1,1,0,0,0,0,0],[0,1,1,1,0,0,0,0]], dtype=np.bool_)
        >>> x = np.array([[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]], dtype=np.bool_)
        >>> z = np.array([[0,0,0,0,0,1,0,0], [0,0,0,0,0,0,1,0]], dtype=np.bool_)
        >>> arange = [0,1]
        >>> a = make_commute_hyper(a, x, z, arange)
        >>> a.astype(int)
        array([[1, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0]])

    See Also:
    _make_commute_hyper

    """
    if not (is_symplectic_form(a) and is_symplectic_form(x) and is_symplectic_form(z)):
        raise QiskitError(f"Input matrices/vectors must be GF(2) symplectic matrices/vectors")

    def make_list(srange):
        if srange is not None:
            try:
                srange = list(srange)
            except TypeError:
                raise QiskitError(f"Input range {srange} is not iterable")

        return srange

    arange = make_list(arange)
    xrange = make_list(xrange)
    zrange = make_list(zrange)

    a = np.array(a)
    # if 'a' is a vector then return result as a vector otherwise return as a matrix
    squeeze = not bool(a.ndim - 1)

    a = np.atleast_2d(a)
    x = np.atleast_2d(np.array(x))
    z = np.atleast_2d(np.array(z))

    if not (a.shape[1] == x.shape[1] == z.shape[1]):
        raise QiskitError(f"Input matrices/vectors must have the same number of columns/length")

    return _make_commute_hyper(a, x, z, arange, xrange, zrange, squeeze)


# ---------------------------------------------------------------


def _make_commute_hyper(
    a: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    arange: Optional[ArrayLike] = None,
    xrange: Optional[ArrayLike] = None,
    zrange: Optional[ArrayLike] = None,
    squeeze: bool = True,
) -> np.ndarray:
    r"""Makes an element(s) commute with hyperbolic pair(s)

    Let a = [a_0,...,a_(k-1)] where a_i are GF(2) symplectic vectors. Let
    x = [x_0,...,x_(l-1)] and z =[z_0,...,z_(l-1)] where x_i and z_i are
    GF(2) symplectic vectors such that (x_i,z_i) are hyerbolic pairs from
    the hyperbolic a basis <x_0,...,x_(l-1),z_0,...,x_(l-1)>. It is assumed
    that {a_0,...,a_(k-1),x_0,...,x_(l-1),z_0,...,x_(l-1)} is an independent
    set of symplectic vectors.

    This method returns a set of vectors b = [b_0,...,b_(k-1)] such that

    1) b_0, ..., b_(k-1) each have zero symplectic product with each of the
    hyperbolic vectors x_0,...,x_(l-1),z_0,...,x_(l-1)
    2) span(b_i, x_j,z_j)  = span(a_i, x_j,z_j)
    for i=0,..,k-1 and j=0,...,l-1

    If the symplectic vectors are considered as Pauli operators then the
    method returns a set of operators [op_(b_0),...,op_(b_(k-1))] such that

    1) op_(b_0), ..., op_(b_(k-1)) each commute with the each of the
    hyperbolic operators op_(x_0),...,op_(x_(l-1)),op_(z_0),...,op_(x_(l-1))
    2) <op_(b_i), op_(x_j),op_(z_j)> = <op_(a_i), op_(x_j),op_(z_j)>
    for each i=0,...,k-1 and j=0,...,l-1

    Args:
        a: Input GF(2) symplectic vectors
        x,z: GF(2) hyperbolic pair vectors
        arange (optional): range of indices from a to make commute. Defaults to None.
        xrange (optional): range of indices from x to use. Defaults to None.
        zrange (optional): range of indices from z to use. Defaults to None.
        squeeze (optional): squeeze = True will return a vector if a vector results

    Returns:
        out: GF(2) symplectic vectors that commute with the given hyperbolic pairs

    Examples:
        >>> a = np.array([[1,1,1,0,0,0]],dtype=np.bool_)
        >>> x = np.array([[0,0,1,0,0,0]],dtype=np.bool_)
        >>> z = np.array([[0,0,0,0,0,1]],dtype=np.bool_)
        >>> a = _make_commute_hyper(a, x, z, squeeze=True)
        >>> a.astype(int)
        array([1, 1, 0, 0, 0, 0])

        >>> a = np.array([[1,1,1,0,0,0,0,0]], dtype=np.bool_)
        >>> x = np.array([[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]], dtype=np.bool_)
        >>> z = np.array([[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0]], dtype=np.bool_)
        >>> xrange = [0,1]
        >>> zrange = [0,1]
        >>> a = _make_commute_hyper(a, x, z, xrange = xrange, zrange=zrange, squeeze=False)
        >>> a.astype(int)
        array([[1, 0, 0, 0, 0, 0, 0, 0]])

        >>> a = np.array([[1,1,1,0,0,0,0,0],[0,1,1,0,0,0,0,0]], dtype=np.bool_)
        >>> x = np.array([[0,1,0,0,0,0,0,0]], dtype=np.bool_)
        >>> z = np.array([[0,0,0,0,0,1,0,0]], dtype=np.bool_)
        >>> arange = [0,1]
        >>> a = _make_commute_hyper(a, x, z, arange)
        >>> a.astype(int)
        array([[1, 0, 1, 0, 0, 0, 0, 0],
               [0, 0, 1, 0, 0, 0, 0, 0]])

        >>> a = np.array([[1,1,1,0,0,0,0,0],[0,1,1,1,0,0,0,0]], dtype=np.bool_)
        >>> x = np.array([[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]], dtype=np.bool_)
        >>> z = np.array([[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0]], dtype=np.bool_)
        >>> arange = [0,1]
        >>> a = _make_commute_hyper(a, x, z, arange)
        >>> a.astype(int)
        array([[1, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0]])

    See Also:
    make_commute_hyper

    """
    # Assuming that a, x, y are all 2d numpy arrays
    a = a.copy()

    if arange is None:
        arange = range(a.shape[0])
    if xrange is None:
        xrange = range(x.shape[0])
    if zrange is None:
        zrange = range(z.shape[0])

    for i, j in zip(xrange, zrange):
        for k in arange:
            # a[k] = _commute_with_pair(a[k], x[i], z[j])
            num_qubits = x[i].shape[0] >> 1
            if _symplectic_product_vv(a[k], x[i], num_qubits):
                a[k] = a[k] ^ z[j]
            if _symplectic_product_vv(a[k], z[j], num_qubits):
                a[k] = a[k] ^ x[i]
    if squeeze:
        a = np.squeeze(a)
    return a


# ---------------------------------------------------------------


def locate_hyper_partner(
    matrix: np.ndarray, vector: np.ndarray
) -> Union[None, Tuple[np.ndarray, int]]:
    """Locate a hyperbolic/noncommutative parter

    Let [m_0,...,m_(k-1)] be the input search matrix and v be the input vector.
    This method search through the vectors m_i to locate a (hyperbolic) vector that has
    a symplectic product with v of 1. If a hyperbolic partner 'av' is found then it
    and its index in the search matrix is turned as a tuple (av, index). If no such
    vector exists then None value is returned.

    Args:
        matrix: GF(2) symplectic matrix (for search)
        vector: GF(2) symplectic vector to find a hyperbolic pair for

    Raises:
        QiskitError: Input matrix must be a GF(2) symplectic matrix
            and input vector must be a GF(2) symplectic vector
        QiskitError: Input {matrix} must be a 2 dimensional array
        QiskitError: Input {vector} must be a 1 dimensional array
        QiskitError: Input matrix and vector must have the same number
            of columns/length

    Returns:
        (av, index) or None: Tuple of the found hyperbolic partner and its index in the
            search matrix. 'None' if no hyperbolic partner is found.

    Examples:
        >>> matrix = np.array([[1,0,1,0,0,0,0,0],[0,1,1,0,0,0,0,0]], dtype=np.bool_)
        >>> vector = np.array([0,0,0,0,0,1,0,0], dtype=np.bool_)
        >>> av, index = locate_hyper_partner(matrix, vector)
        >>> av.astype(int)
        array([0, 1, 1, 0, 0, 0, 0, 0])
        >>> index
        1

    See Also:
    _locate_hyper_partner, build_hyper_partner, _build_hyper_partner
    """
    if not (is_symplectic_form(matrix) and is_symplectic_vector_form(vector)):
        raise QiskitError(
            f"Input {matrix} must be a GF(2) symplectic matrix\
            and input {vector} must be a GF(2) symplectic vector"
        )

    if not matrix.ndim == 2:
        raise QiskitError(f"Input {matrix} must be a 2 dimensional array")

    if not vector.ndim == 1:
        raise QiskitError(f"Input {vector} must be a 1 dimensional array")

    if not matrix.shape[1] == vector.shape[0]:
        raise QiskitError(
            f"Input matrix and vector must have the same number \
            of columns/length {matrix.shape[1]}!={vector.shape[0]}"
        )

    return _locate_hyper_partner(matrix, vector)


def _locate_hyper_partner(matrix: np.ndarray, vector: np.ndarray) -> Union[None, np.ndarray]:
    """Locate a hyperbolic/noncommutative parter

    Let [m_0,...,m_(k-1)] be the input search matrix and v be the input vector.
    This method search through the vectors m_i to locate a (hyperbolic) vector that has
    a symplectic product with v of 1. If a hyperbolic partner 'av' is found then it
    and its index in the search matrix is turned as a tuple (av, index). If no such
    vector exists then None value is returned.

    Args:
        matrix: GF(2) symplectic matrix (for search)
        vector: GF(2) symplectic vector to find a hyperbolic pair for

    Returns:
        (av, index) or None: Tuple of the found hyperbolic partner and its index in the
            search matrix. 'None' if no hyperbolic partner is found.

    Examples:
        >>> matrix = np.array([[1,0,1,0,0,0,0,0],[0,1,1,0,0,0,0,0]], dtype=np.bool_)
        >>> vector = np.array([0,0,0,0,0,1,0,0], dtype=np.bool_)
        >>> av, index = locate_hyper_partner(matrix, vector)
        >>> av.astype(int)
        array([0, 1, 1, 0, 0, 0, 0, 0])
        >>> index
        1

    See Also:
    locate_hyper_partner, build_hyper_partner, _build_hyper_partner
    """
    n = matrix.shape[1] >> 1
    for index, item in enumerate(matrix):
        if _symplectic_product_vv(item, vector, n) == 1:
            return (item.copy(), index)
    return None


# ---------------------------------------------------------------


def build_hyper_partner(matrix, index: int) -> np.ndarray:
    """Builds an independent hyperbolic partner for the input vector indexed

    Let the input matrix be [m_0,...,m_(k-1)] and v = m_index. It is assumed
    that the vectors m_0,...,m_(k-1) have a zero pairwise
    symplectic product (i.e. represent a set of Pauli operators that pairwise
    commmute). This method will find a GF(2) symplectic vector that
    has zero symplectic product with each m_i != v and a symplectic product of 1
    with the vector v. If the vectors represent Pauli operators then this
    method will find a Pauli operator that commute with the operators
    represented by the vectors m_i != v and that anticommutes with the Pauli
    operator that is represented by v.

    Args:
        matrix: GF(2) symplectic matrix representing a set of independent
            commuting generators
        index: index of generator to build a hyperbolic partner for

    Raises:
        QiskitError: Input matrix must be a GF(2) symplectic matrix
        QiskitError: Input matrix must represent a set of commuting operators
        QiskitError: Input matrix does not represent a set of independent
            operators, it does not have have full rank
        QiskitError: Input index out or range

    Returns:
        out: a hyperbolic partner for the given vector wrt the set of commuting
            generators

    Examples:
        >>> matrix = np.array(
            [[1,0,0,0,0,0,0,0],
             [0,1,0,0,0,0,0,0],
             [0,0,1,0,0,0,0,0],
             [0,0,0,1,0,0,0,0]], dtype=np.bool_)
        >>> av = build_hyper_partner(matrix, 0)
        >>> av.astype(int)
        array([0, 0, 0, 0, 1, 0, 0, 0])

    See Also:
    _build_hyper_partner, locate_hyper_partner, _locate_hyper_partner

    Notes: This method is the implementation of Proposition 10.4 from Nielsen
    and Chuang's Quantum Computation and Quantum Information

    """
    matrix = np.atleast_2d(np.array(matrix, dtype=bool))

    if not is_symplectic_form(matrix):
        raise QiskitError(f"Input {matrix} must be a GF(2) symplectic matrix")

    # matrix -> all associated operators must commute
    if not all_commute(matrix) == True:
        raise QiskitError("Input matrix must represent a set of commuting operators")

    rank_ = mt.rank(matrix)
    if rank_ != matrix.shape[0]:
        raise QiskitError(
            f"Input matrix does not represent a set of independent \
            operators, it does not have have full rank: {rank_}"
        )

    if index not in range(matrix.shape[1] >> 1):
        raise QiskitError(f"Input index out or range: {index}>={matrix.shape[1]>>1}")

    return _build_hyper_partner(matrix, index)


def _build_hyper_partner(matrix, index: int) -> np.ndarray:
    """Builds an independent hyperbolic partner for the input vector indexed

    Let the input matrix be [m_0,...,m_(k-1)] and v = m_index. It is assumed
    that the vectors m_0,...,m_(k-1) have a zero pairwise
    symplectic product (i.e. represent a set of Pauli operators that pairwise
    commmute). This method will find a GF(2) symplectic vector that
    has zero symplectic product with each m_i != v and a symplectic product of 1
    with the vector v. If the vectors represent Pauli operators then this
    method will find a Pauli operator that commute with the operators
    represented by the vectors m_i != v and that anticommutes with the Pauli
    operator that is represented by v.

    Args:
        matrix: GF(2) symplectic matrix representing a set of independent
            commuting generators
        index: index of generator to build a hyperbolic partner for

    Returns:
        out: a hyperbolic partner for the given vector wrt the set of commuting
            generators

    Examples:
        >>> matrix = np.array(
            [[1,0,0,0,0,0,0,0],
             [0,1,0,0,0,0,0,0],
             [0,0,1,0,0,0,0,0],
             [0,0,0,1,0,0,0,0]], dtype=np.bool_)
        >>> av = _build_hyper_partner(matrix, 0)
        >>> av.astype(int)
        array([0, 0, 0, 0, 1, 0, 0, 0])

    See Also:
    build_hyper_partner, locate_hyper_partner, _locate_hyper_partner

    Notes: This method is the implementation of Proposition 10.4 from Nielsen
    and Chuang's Quantum Computation and Quantum Information

    """

    nrows = matrix.shape[0]
    ncols = matrix.shape[1]

    _lambda = mt.create_lambda_matrix(ncols >> 1)
    slambda = np.matmul(matrix, _lambda)

    heads, rref_mat, transform_mat, rank = mt._rref_complete(slambda)

    e_index = np.zeros(nrows, dtype=bool)
    e_index[index] = True

    trans_e_index = np.matmul(transform_mat, e_index)

    pivot = 0
    result = np.zeros(ncols, dtype=bool)
    for i in range(ncols):
        if heads[i] == 1:
            result[i] = trans_e_index[pivot]
            pivot = +1

    return result


# ---------------------------------------------------------------


def symplectic_gram_schmidt(
    a: ArrayLike, x: Optional[np.ndarray] = None, z: Optional[np.ndarray] = None
):
    """Applies the sympletic Gram-Schmidt process to the input matrix

    Apply the symplectic GramSchmidt process to the input symplectic matrix. Resulting
    hyperbolic pairs are added to x and z arrays. Elements of the center will we added to the
    center array.

    Args:
        a: Symplectic matrix
        x, z (optional): GF(2) Symplectic matrices representing hyperbolic pairs to
        build upon. Default is None.

    Raises:
        QiskitError Input matric not a GF(2) symplectic matrix
        QiskitError: Input hyperbolic array x is not a GF(2) sympletic matrix
        QiskitError: Input hyperbolic array z is not a GF(2) sympletic matrix
        QiskitError: Input hyperbolic arrays have different dimensions
        QiskitError: Input hyperbolic matrices do not represent a hyperbolic basis

    Returns:
        center, x, z: Center array and hyperbolic pairs split accross x and z

    Examples:
        >>> a = np.array([[0,1,0,0,1,0,1,0],
                          [0,0,0,0,1,1,0,1],
                          [1,1,1,0,0,1,0,0],
                          [1,1,0,1,0,0,0,0]], dtype=np.bool_)
        >>> center_, x, z = symplectic_gram_schmidt(a)
        >>> center_.astype(int)
        array([[1, 1, 1, 0, 1, 0, 0, 1],
               [1, 0, 0, 1, 0, 1, 1, 1]])
       >>> x.astype(int)
       array([[0, 1, 0, 0, 1, 0, 1, 0]])
       >>> z.astype(int)
       array([[0, 0, 0, 0, 1, 1, 0, 1]])

    Also See:
    _symplectic_gram_schmid

    TODO: Add an example that shows using the optional x and z arrays
    """

    a = np.atleast_2d(np.array(a))
    if not is_symplectic_matrix_form(a):
        raise QiskitError(f"Input matrix not a GF(2) symplectic matrix")
    try:
        x = [item for item in x]
        if not is_symplectic_vector_form(x[0]):
            raise QiskitError(f"Input hyperbolic array x is not a GF(2) sympletic matrix")
    except TypeError:
        x = []

    try:
        z = [item for item in z]
        if not is_symplectic_vector_form(z[0]):
            raise QiskitError(f"Input hyperbolic array z is not a GF(2) sympletic matrix")
    except TypeError:
        z = []

    if not len(x) == len(z):
        raise QiskitError(f"Input hyperbolic arrays have different dimensions")

    if len(x) > 0 and x[0].shape[0] != z[0].shape[0]:
        raise QiskitError(f"Input hyperbolic arrays have different dimensions")

    if not is_hyper_form(x, z):
        raise QiskitError(f"Input hyperbolic matrices do not represent a hyperbolic basis")

    return _symplectic_gram_schmidt(a, x, z)


def _symplectic_gram_schmidt(a: np.ndarray, x: List[np.ndarray], z: List[np.ndarray]):
    """Applies the sympletic Gram-Schmidt process to the input matrix

    Apply the symplectic GramSchmidt process to the input symplectic matrix. Resulting
    hyperbolic pairs are added to x and z arrays. Elements of the center will we added to the
    center array.

    Args:
        a: Symplectic matrix
        x, z (optional): GF(2) Symplectic matrices representing hyperbolic pairs to
        build upon. Default is None.

    Raises:
        QiskitError Input matric not a GF(2) symplectic matrix
        QiskitError: Input hyperbolic array x is not a GF(2) sympletic matrix
        QiskitError: Input hyperbolic array z is not a GF(2) sympletic matrix
        QiskitError: Input hyperbolic arrays have different dimensions

    Returns:
        center, x, z: Center array and hyperbolic pairs split accross x and z

    Examples:
        >>> a = np.array([[0,1,0,0,1,0,1,0],
                          [0,0,0,0,1,1,0,1],
                          [1,1,1,0,0,1,0,0],
                          [1,1,0,1,0,0,0,0]], dtype=np.bool_)
        >>> center_, x, z = _symplectic_gram_schmidt(a, [], [])
        >>> center_.astype(int)
        array([[1, 1, 1, 0, 1, 0, 0, 1],
               [1, 0, 0, 1, 0, 1, 1, 1]])
        >>> x.astype(int)
        array([[0, 1, 0, 0, 1, 0, 1, 0]])
        >>> z.astype(int)
        array([[0, 0, 0, 0, 1, 1, 0, 1]])

    Also See:
    _symplectic_gram_schmid

    TODO: Add an example that shows using the optional x and z arrays
    """

    a = a.copy()
    a_view = a
    center_ = []

    while a_view.shape[0] > 0:
        elem = a_view[0]
        # Remove elem from a_view
        a_view = a_view[1:]
        try:
            elem_p, index = _locate_hyper_partner(a_view, elem)
            x.append(elem)
            z.append(elem_p)

            # Revove elem_p from a_view
            temp_view = a_view[:-1]
            temp_view[index:] = a_view[index + 1 :]
            a_view = temp_view

            a_view = make_commute_hyper(a_view, elem, elem_p)
            # TODO: Change the above line to use the faster version
            # a_view = _make_elements_commute_with_hyper_pair(
            #    a_view, range(a_view.shape[0]), elem, elem_p
            # )
        except TypeError:
            center_.append(elem)

    x = np.asarray(x)
    z = np.asarray(z)
    if len(center_) == 0:
        center_ = np.zeros(shape=(0, x.shape[1]), dtype=np.bool_)
    else:
        center_ = np.asarray(center_)

    return center_, x, z


# ---------------------------------------------------------------


def is_symplectic_matrix_form(
    matrix: ArrayLike, dtype: Optional[Union[bool, np.bool_, int, np.integer]] = None
) -> bool:
    """Is the input matrix GF(2) symplectic

    Checks if the given array like matrix is in the form of a symplectic matrix:
    two dimensional, even number of columns, 0/1 or boolean entries. The optional
    argument dtype can be given to check if entries are a specific dtype

    Args:
        matrix: Input matrix to be checked
        dtype: Optional. Check if given matrix is of type dtype. Default: None

    Returns:
        out: True if the input matrix is GF(2) symplectic. False otherwise.

    Examples:
        >>> matrix = numpy.array([[1,0,0,1,0,0,1,0],[0,1,1,1,0,0,0,1]], dtype=numpy.bool_)
        >>> is_symplectic_matrix_form(matrix)
        True

        >>> matrix = numpy.array([[1,0,0,1,0,0,1],[0,1,1,1,0,0,0]], dtype=numpy.bool_)
        >>> is_symplectic_matrix_form(matrix)
        False

        >>> matrix = numpy.array([[1,0,0,1],[0,1,1,1]], dtype=numpy.bool_)
        >>> is_symplectic_matrix_form(matrix, dtype=int)
        False

    See Also:
    is_symplectic_vector_form, is_symplectic_form
    """
    matrix = np.array(matrix)

    if matrix.ndim != 2:
        return False
    if matrix.shape[1] % 2:
        return False
    if not np.array_equal(matrix, matrix % 2):
        return False
    if dtype is None:
        return True
    if not isinstance(matrix[0][0], dtype):
        return False
    return True


def is_symplectic_vector_form(
    vector: ArrayLike, dtype: Optional[Union[bool, np.bool_, int, np.integer]] = None
) -> bool:
    """Is the input vector GF(2) symplectic

    Checks if the given array like vector is in the form of a symplectic vector:
    two dimensional, even number of columns, 0/1 or boolean entries. The optional
    argument dtype can be given to check if entries are a specific dtype

    Args:
        vector: Input vector to be checked
        dtype: Optional. Check if given vector is of type dtype. Default: None

    Returns:
        out: True if the input vector is GF(2) symplectic. False otherwise.

    Examples:
        >>> vector = numpy.array([[1,0,0,1,0,0,1,0],[0,1,1,1,0,0,0,1]], dtype=numpy.bool_)
        >>> is_symplectic_vector_form(vector)
        False

        >>> vector = numpy.array([[1,0,0,1,0,0,1,0]], dtype=numpy.bool_)
        >>> is_symplectic_vector_form(vector)
        False

        >>> vector = numpy.array([[1,0,0,1,0,0,1,0]], dtype=numpy.bool_)
        >>> is_symplectic_vector_form(vector)
        True

        >>> vector = numpy.array([1,0,0,1], dtype=numpy.int8)
        >>> is_symplectic_vector_form(vector, dtype=numpy.int8)
        False

    See Also:
    is_symplectic_matrix_form, is_symplectic_form
    """
    vector = np.array(vector)
    if not vector.ndim == 1:
        return False
    if vector.shape[0] % 2:
        return False
    if not np.array_equal(vector, vector % 2):
        return False
    if dtype is None:
        return True
    if not isinstance(vector[0], dtype):
        return False
    return True


def is_symplectic_form(
    a: Any, dtype: Optional[Union[bool, np.bool_, int, np.integer]] = None
) -> bool:
    """Is the input a GF(2) symplectic matrix or vector

    Args:
        a: Input to be checked
        dtype (Optional): Check if given matrix/vector is of type dtype. Defaults to None.

    Returns:
        out: True if the input matrix/vector is GF(2) symplectic. False otherwise.

    Examples:
        >>> vector = numpy.array([[1,0,0,1,0,0,1,0]], dtype=numpy.bool_)
        >>> is_symplectic_form(vector)
        True

        >>> matrix = numpy.array([[1,0,0,1,0,0,1],[0,1,1,1,0,0,0]], dtype=numpy.bool_)
        >>> is_symplectic_form(matrix)
        False

    See Also:
    is_symplectic_vector_form, is_symplectic_matrix_form
    """

    a = np.array(a)
    if a.ndim == 1:
        return is_symplectic_vector_form(a, dtype=dtype)
    elif a.ndim == 2:
        return is_symplectic_matrix_form(a, dtype=dtype)
    else:
        return False


# ---------------------------------------------------------------


def is_center(center_: ArrayLike, matrix: ArrayLike) -> bool:
    """Does the input center matrix represent the center of the supplied matrix?

    Let op(center_) = [op(c_0),op(c_1),...,op(c_(k-1))] be the Pauli operators
    represented by the input center_matrix. Let op(matrix_) = [op(m_0),op(m_1),...,op(m_(t-1))]
    be the Pauli operators represented by the input matrix. This method
    returns True if

    <op(center_)> = Z(<op(matrix)>)

    and False otherwise.

    Args:
        center_ (ArrayLike): Generators of center to be checked
        matrix (ArrayLike): Generators of full group

    Raises:
        QiskitError: Not all inputs are not GF(2) symplectic matrices/vectors

    Returns:
        out: True if <op(center_)> = Z(<op(matrix)>), False otherwise

    Examples:
        >>> matrix = np.array(
            [[0,1,0,0,1,0,1,0],
            [0,0,0,0,1,1,0,1],
            [1,1,1,0,0,1,0,0],
            [1,1,0,1,0,0,0,0]], dtype=np.bool_)
        >>> center_ = np.array([[1, 1, 1, 0, 1, 0, 0, 1],
                                [1, 0, 0, 1, 0, 1, 1, 1]], dtype=np.bool_)
        >>> is_center(center_, matrix)
        True

    See Also:
    center
    """
    matrix = np.atleast_2d(np.array(matrix))
    center_ = np.atleast_2d(np.array(center_))
    if not (is_symplectic_matrix_form(center_) and is_symplectic_matrix_form(matrix)):
        QiskitError(f"Not all inputs are not GF(2) symplectic matrices")
    cal_center = center(matrix)

    return is_same_span(center_, cal_center)


# ---------------------------------------------------------------


def is_same_span(matrix1: ArrayLike, matrix2: ArrayLike) -> bool:
    """Does span(rows of matrix1) = span(rows of matrix2)?

    Args:
        matrix1: First set of vectors
        matrix2: Second set of vectors

    Returns:
        out: True if span(rows of matrix1) = span(rows of matrix2). False otherwise

    Raises:
        QiskitError: Inpiut matrices must by GF(2) symplectic matrices

    Examples:
        >>> matrix1 = numpy.array([[1,1,0,0],[0,0,1,0]], dtype=numpy.bool_)
        >>> matrix2 = numpy.array([[1,1,1,0],[0,0,1,0]], dtype=numpy.bool_)
        >>> is_same_span(matrix1, matrix2)
        True

        >>> matrix1 = numpy.array([[1,1,0,0],[0,0,1,0]], dtype=numpy.bool_)
        >>> matrix2 = numpy.array([[1,1,1,0],[0,0,0,1]], dtype=numpy.bool_)
        >>> is_same_span(matrix1, matrix2)
        False
    """
    matrix1 = np.atleast_2d(np.array(matrix1))
    matrix2 = np.atleast_2d(np.array(matrix2))

    if not (is_symplectic_form(matrix1) and is_symplectic_form(matrix2)):
        raise QiskitError(f"Inpiut matrices must by GF(2) symplectic matrices")

    if matrix1.shape[1] != matrix2.shape[1]:
        return False

    _, rref_matrix1, _, rank_matrix1 = mt.rref_complete(matrix1)
    _, rref_matrix2, _, rank_matrix2 = mt.rref_complete(matrix2)

    if rank_matrix1 != rank_matrix2:
        return False

    rref_matrix1 = rref_matrix1[:rank_matrix1]
    rref_matrix2 = rref_matrix2[:rank_matrix2]
    return np.array_equal(rref_matrix1, rref_matrix2)


def is_hyper_form(x: ArrayLike, z: ArrayLike) -> bool:
    """Do the input matrices form a hyperbolic/symnplectic basis?

    Args:
        x,z: Pairs to test if they form a hyperbolic/symnplectic basis

    Returns:
        out: True if input matrices form a hyperbolic/symplectic basis

    Examples:
        >>> x = numpy.array([[1,0,0,0,],[0,1,0,0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0,0,1,0,],[0,0,0,1]], dtype=numpy.bool_)
        >>> is_hyper_form(x,z)
        True
    """

    matrix = np.vstack((x, z))
    test = symplectic_product(matrix, matrix)
    return np.array_equal(test, mt._create_lambda_matrix(matrix.shape[0] >> 1))


def is_stabilizer_group(matrix: ArrayLike) -> bool:
    """Do the rows of the input matrix represent the generators of an
    abelian Pauli subgroup?

    Args:
        matrix: GF(2) symplectic matrix

    Raises:
        QiskitError: Input matrix not a GF(2) symplectic matrix

    Returns:
        out: True is matrix represents a Stabilizer group

    Examples:
        >>> matrix = numpy.array([[1,0,0,0,0,0],
                               [0,1,0,0,0,0],
                               [1,0,0,0,0,1]], dtype=np.bool_)
        >>> is_stabilizer_group(matrix)
        True
    """
    matrix = np.atleast_2d(np.array(matrix))
    if not is_symplectic_matrix_form(matrix):
        raise QiskitError(f"Input matrix not a GF(2) symplectic matrix")
    return all_commute(matrix)


# ---------------------------------------------------------------


def center(matrix: ArrayLike, preserve: bool = False) -> np.ndarray:
    """Find the center of the group with generators given by the symplectic matrix

    Args:
        matrix: GF(2) symplectic matrix
        preserve: If True then an attempt will be made to preserve then generators form

    Raises:
        QiskitError: Input matrix is not a symplectic matrix

    Returns:
        out : Generators for the center, represented as a symplectic matrix,
            of the group with generators given by the input symplectic matrix

    Examples:
        >>> matrix = numpy.array([[1,1,0,0,1,0],
                               [0,0,0,1,1,0],
                               [0,1,0,1,0,1]], dtype=numpy.bool_)
       >>> center_ = center(matrix)
       >>> center_.astype(int)
       array([[1, 1, 0, 0, 1, 0]])

    See Also:
    _center, _center_preserve

    TODO: Add in example with preserve=True is useful
    """
    matrix = np.atleast_2d(np.array(matrix))
    if not is_symplectic_matrix_form(matrix):
        raise QiskitError("Input matrix is not a symplectic matrix")
    if preserve:
        return _center_preserve(matrix)
    else:
        return _center(matrix)


def _center(matrix: np.ndarray) -> bool:
    """Find the center of the group with generators given by the symplectic matrix

    Args:
        matrix: GF(2) symplectic matrix

    Returns:
        out : Generators for the center, represented as a symplectic matrix,
            of the group with generators given by the input symplectic matrix

    Examples:
        >>> matrix = numpy.array([[1,1,0,0,1,0],
                               [0,0,0,1,1,0],
                               [0,1,0,1,0,1]], dtype=numpy.bool_)
       >>> center_ = _center(matrix)
       >>> center_.astype(int)
       array([[1, 1, 0, 0, 1, 0]])

    Notes:
    This method may nbot preserve any input vectors

    See Also:
    center

    TODO: Add in example with preserve=True is useful
    """
    center, hyper1, hyper2 = _symplectic_gram_schmidt(matrix, [], [])
    return center


def _center_preserve(matrix: np.ndarray):
    """Find the center of the group with generators given by the symplectic matrix

    Args:
        matrix: GF(2) symplectic matrix

    Returns:
        out : Generators for the center, represented as a symplectic matrix,
            of the group with generators given by the input symplectic matrix

    Examples:
        >>> matrix = numpy.array([[1,1,0,0,1,0],
                               [0,0,0,1,1,0],
                               [0,1,0,1,0,1]], dtype=numpy.bool_)
       >>> center_ = _center_preserve(matrix)
       >>> center_.astype(int)
       array([[1, 1, 0, 0, 1, 0]])

    Notes:
    This method attempots to preserve any input vectors

    See Also:
    _center, center

    TODO: Add in example with preserve is actually needed
    """
    # First move any generator that is in the center to the front of the list
    rematrix = deque()
    num = matrix.shape[1] >> 1
    for i, opi in enumerate(reversed(matrix)):
        break_flag = False
        for j, opj in enumerate(matrix):
            if _symplectic_product_vv(opi, opj, num) == 1:
                break_flag = True
                break
        if break_flag:
            rematrix.append(opi)
        else:
            rematrix.appendleft(opi)
    # Now calculate the center as normal
    rematrix = np.array(rematrix)
    return _center(rematrix)


# ---------------------------------------------------------------


def extend_basis_to_pauli_group(matrix: ArrayLike) -> np.ndarray:
    """Given a set of generators (not necessarily independent) find
    a full basis using as many of the provided generators.

    Args:
        matrix (ArrayLike): Set of generators (symplectic vectors)

    Returns:
        np.ndarray: A maximal independant set
    """
    matrix = np.atleast_2d(np.asarray(matrix))
    assert is_symplectic_matrix_form(matrix), QiskitError(f"Input matrix not in a symplectic form")
    return _extend_basis_to_pauli_group(matrix)


def _extend_basis_to_pauli_group(matrix: np.ndarray) -> np.ndarray:
    """Given a set of generators (not necessarily independent) find
    a full basis using as many of the provided generators.

    Args:
        matrix (np.ndarray): Set of genera*tors (symplectic vectors)

    Returns:
        np.ndarray: A maximal independant set

    """
    aug_matrix = mt.augment_mat(matrix, "bottom")
    heads, rref_mat, transform_mat, rank = mt._rref_complete(aug_matrix.T)
    shape = (rank, aug_matrix.shape[1])
    ext_matrix = np.zeros(shape, dtype=np.bool_)
    posns = np.flatnonzero(heads)
    for k, index in enumerate(posns):
        ext_matrix[k] = aug_matrix[index]
    return ext_matrix


def extend_center_into_hyper_form(in_hyper1, in_hyper2, in_center):
    """Given a basis for a group in the form of center, hyper1 and hyper2 find hyperbolic partners for
    each element of the center (acenter) such that hyper1 = <hyper1, center> and hyper2 = <hyper2, acenter>
    are hyperbolic pairs. Thus <hyper1, hyper2> will be isomorphic some P_m

    Args:
        in_hyper1 ([type]): [description]
        in_hyper2 ([type]): [description]
        in_center ([type]): [description]
    """
    if in_center is None:
        in_center = []

    in_center = np.atleast_2d(np.array(in_center))
    hyper1 = np.atleast_2d(in_hyper1)
    hyper2 = np.atleast_2d(in_hyper2)

    if in_center.shape[1] == 0:
        in_center = np.zeros((0, in_hyper1.shape[1]))

    assert is_center(center, np.vstack((hyper1, hyper2, center))), QiskitError(
        f"Input center is not center"
    )
    assert is_hyper_form(hyper1, hyper2), QiskitError(
        f"Inpur hyper1,hyper2 are not in hyperbolic pairs"
    )

    return _extend_center_into_hyper_form(in_hyper1, in_hyper2, in_center)


def _extend_center_into_hyper_form(in_hyper1, in_hyper2, in_center):
    """Given a basis for a group in the form of center, hyper1 and hyper2 find hyperbolic partners for
    each element of the center (acenter) such that hyper1 = <hyper1, center> and hyper2 = <hyper2, acenter>
    are hyperbolic pairs. Thus <hyper1, hyper2> will be isomorphic some P_m

    # Change the name to something like _extend_iso_hyper_form_to_hyper_form

    Args:
        in_hyper1 ([type]): [description]
        in_hyper2 ([type]): [description]
        in_center ([type]): [description]

    Returns:
        hyper1 (numpy.ndarray), hyper2 (numpy.ndarray)
    """

    center_size = in_center.shape[0]
    hyper_size = in_hyper1.shape[0]
    shape = (center_size, in_hyper1.shape[1])
    spacer = np.zeros(shape, dtype=np.bool_)
    hyper1 = np.vstack((in_hyper1, spacer))
    hyper2 = np.vstack((in_hyper2, spacer))
    center_ = in_center.copy()

    while center_size > 0:
        hop = _build_hyper_partner(center_[:center_size], center_size - 1)
        hop = _make_element_commute_with_hyper_pairs(
            hop, hyper1, hyper2, range(hyper_size), range(hyper_size)
        )

        center_size -= 1
        hyper1[hyper_size] = center_[center_size]
        hyper2[hyper_size] = hop
        hyper_size += 1

    return hyper1, hyper2


def isotropic_hyperbolic_form(matrix):
    return _symplectic_gram_schmidt(matrix, [], [])


def isotropic_hyperbolic_basis(
    matrix: Union[None, ArrayLike], hyper1: Union[None, ArrayLike], hyper2: Union[None, ArrayLike]
):

    if (hyper1 is None) ^ (hyper2 is None):
        raise QiskitError(f"hyper1 and hyper2 must be both be None or both be array like")

    if hyper1 is not None:
        hyper1 = np.array(hyper1)
        hyper2 = np.array(hyper2)
        assert is_symplectic_matrix_form(hyper1), QiskitError(f"{hyper1} not a symplectic matrix")
        assert is_symplectic_matrix_form(hyper2), QiskitError(f"{hyper2} not a symplectic matrix")

        if matrix is not None:
            matrix = np.array(matrix)
            assert is_symplectic_matrix_form(matrix), QiskitError(
                f"{matrix} not a symplectic matrix"
            )
            assert is_center(matrix, np.vstack(matrix, hyper1, hyper2))
        else:
            matrix = None
    else:
        assert matrix is not None, QiskitError(
            "At one of matrix, hyper1 and hyper2 cannot all be None"
        )
        matrix = np.array(matrix)
        assert is_symplectic_matrix_form(matrix), QiskitError(f"{matrix} not a symplectic matrix")
        hyper1 = []
        hyper2 = []

    return _isotropic_hyperbolic_basis(matrix, hyper1, hyper2)


def _isotropic_hyperbolic_basis(
    matrix: Union[None, np.ndarray], hyper1: List[np.ndarray], hyper2: List[np.ndarray]
):
    """[summary]

    Args:
        matrix (np.ndarray): matrix if hyper1 and hyper2 are None, center otherwise
        hyper1 (List[np.ndarray]): [description]
        hyper2 (List[np.ndarray]): [description]

    Returns:
    """

    if len(hyper1) == 0:
        center_, hyper1, hyper2 = _symplectic_gram_schmidt(matrix, [], [])
    else:
        center_ = matrix

    if center_ is not None:
        hyper1, hyper2 = _extend_center_into_hyper_form(hyper1, hyper2, center_)

    basis = _extend_basis_to_pauli_group(np.vstack((hyper1, hyper2)))

    added = hyper1.shape[1] - 2 * hyper1.shape[0]

    basis_com = _make_elements_commute_with_hyper_pairs(
        basis[0:added],
        range(2 * hyper1.shape[0], basis.shape[0]),
        hyper1,
        range(hyper1.shape[0]),
        hyper2,
        range(hyper1.shape[0]),
    )

    _, hyper1_ans, hyper2_ans = symplectic_gram_schmidt(basis_com, hyper1, hyper2)

    return hyper1_ans, hyper2_ans


def remove_hyper_elements_from_hyper_form(
    hyper1: ArrayLike, hyper2: ArrayLike, center_: ArrayLike, indices: ArrayLike
):
    """Transfer those elements/vectors from hyper1 with an index in indices into center_
    and delete the corresponding hyperbolic partner from hyper2.

    Args:
        hyper1 (ArrayLike): First array of hyperbolic pairs and source
        hyper2 (ArrayLike): Second array of hyperbolic pairs
        center_ (ArrayLike): center array and sink
        indices (ArrayLike): indices indicating which rows to transfer from source to sink

    Returns:
        Returns: hyper1, hyper2, center_
    """
    hyper1 = np.array(hyper1)
    hyper2 = np.array(hyper2)
    indices = list(indices)

    assert hyper1.ndim > 1, QiskitError(f"hyper1, hyper2 must have the same number of dimensions")
    assert hyper1.shape == hyper2.shape, QiskitError(
        f"hyper1 (shape={hyper1.shape})and hyper2 (shape={hyper2.shape}) \
            must have the same shape"
    )
    if center_ is not None:
        center_ = np.atleast_2d(np.array(center_))
        assert hyper1.ndim == center_.ndim, QiskitError(
            f"center must have the same number of dimensions as hyper1 and hyper2"
        )
        assert hyper1.shape[1] == center_.shape[1], QiskitError(
            f"hyper1 and hyper2 must have the same size in the second dimension as the center_"
        )

    return _remove_hyper_elements_from_hyper_form(hyper1, hyper2, center_, indices)


def _remove_hyper_elements_from_hyper_form(
    hyper1: np.ndarray, hyper2: np.ndarray, center_: np.ndarray, indices: List[int]
):
    """Transfer those elements/vectors from hyper1 with an index in indices into center_
    and delete the corresponding hyperbolic partner from hyper2.

    Args:
        hyper1 (np.ndarray): First array of hyperbolic pairs and source
        hyper2 (np.ndarray): Second array of hyperbolic pairs
        center_ (np.ndarray): center array and sink
        indices (List[int]): indices indicating which rows to transfer from source to sink

    Returns: hyper1, hyper2, center_
    """

    rm_size = len(indices)
    size = hyper1.shape[0] - rm_size
    shape = (size, hyper1.shape[1])
    new_hyper1 = np.zeros(shape, dtype=np.bool_)
    new_hyper2 = np.zeros(shape, dtype=np.bool_)
    pos = 0
    part_center = deque()
    for i in range(hyper1.shape[0]):
        if i in indices:
            part_center.appendleft(hyper1[i].copy())
        else:
            new_hyper1[pos] = hyper1[i].copy()
            new_hyper2[pos] = hyper2[i].copy()
            pos += 1

    part_center = np.array(part_center)
    if center_ is None:
        new_center = part_center
    else:
        new_center = np.vstack((center_, part_center))

    return new_center, new_hyper1, new_hyper2


def normalizer(
    matrix: Union[None, ArrayLike] = None,
    hyper1: Union[None, ArrayLike] = None,
    hyper2: Union[None, ArrayLike] = None,
):
    """Return the normalizer of the group generated by the generators represented in the
    symplectic matrix(s):

    Args:
        matrix (Union[None, ArrayLike], optional): [description]. Defaults to None.
        hyper1 (Union[None, ArrayLike], optional): [description]. Defaults to None.
        hyper2 (Union[None, ArrayLike], optional): [description]. Defaults to None.

    Raises:
        QiskitError: [description]

    Returns:
        [type]: [description]
    """

    assert not (matrix is None and hyper1 is None and hyper2 is None), QiskitError(
        "All inputs should not be None"
    )
    if matrix is not None:
        matrix = np.array(matrix)
        assert is_symplectic_matrix_form(matrix), QiskitError(
            f"{matrix} must be a symplectic matrix"
        )

    if (hyper1 is None) ^ (hyper2 is None):
        raise QiskitError(f"hyper1 and hyper2 must be both be None or both be array like")

    if matrix is not None and hyper1 is None:
        if is_stabilizer_group(matrix):
            return _normalizer_stabilizer_group(matrix)
        else:
            return _normalizer_gauge_group(matrix)

    if matrix is None:
        matrix = np.zeros(shape=(0, hyper1.shape[1]), dtype=np.bool_)

    hyper1 = np.array(hyper1)
    hyper2 = np.array(hyper2)
    assert is_symplectic_matrix_form(hyper1) and is_symplectic_matrix_form(hyper2)
    assert hyper1.shape == hyper2.shape, QiskitError(f"hyper1 and hyper2 must have the same shape")
    assert matrix.shape[1] == hyper1.shape[1], QiskitError(
        f"All inputs must have the same number of columns/length"
    )
    return _normalizer_gauge_group_preserve(matrix, hyper1, hyper2)


def _normalizer_stabilizer_group(matrix):
    """[summary]

    Args:
        matrix ([type]): [description]

    Returns:
        [type]: [description]
    """

    dist_center = mt.rank(matrix)

    matrix_ext = _extend_basis_to_pauli_group(matrix)
    center_, hyper1, hyper2 = _symplectic_gram_schmidt(matrix_ext, [], [])
    center_, hyper1, hyper2 = _remove_hyper_elements_from_hyper_form(
        hyper1, hyper2, center_, list(range(dist_center))
    )
    return center_, hyper1, hyper2


def _normalizer_gauge_group_preserve(center_: np.ndarray, hyper1: np.ndarray, hyper2: np.ndarray):
    """[summary]

    Args:
        center_ (np.ndarray): [description]
        hyper1 (np.ndarray): [description]
        hyper2 (np.ndarray): [description]

    Returns:
        [type]: [description]
    """
    center_size = center_.shape[0]
    gauge_degree = hyper1.shape[0]

    hyper1, hyper2 = _extend_center_into_hyper_form(hyper1, hyper2, center_)
    matrix_ext = _extend_basis_to_pauli_group(np.vstack((hyper1, hyper2)))
    matrix_ext = _make_elements_commute_with_hyper_pairs(
        matrix_ext,
        range(hyper1.shape[0] << 1, matrix_ext.shape[0]),
        hyper1,
        range(hyper1.shape[0]),
        hyper2,
        range(hyper2.shape[0]),
    )
    matrix = matrix_ext[hyper1.shape[0] << 1 :]
    lhyper1 = [item.copy() for item in matrix_ext[: hyper1.shape[0]]]
    lhyper2 = [item.copy() for item in matrix_ext[hyper1.shape[0] : hyper1.shape[0] << 1]]
    center_, hyper1, hyper2 = _symplectic_gram_schmidt(matrix, lhyper1, lhyper2)
    indices = list(range(gauge_degree, gauge_degree + center_size))
    return _remove_hyper_elements_from_hyper_form(hyper1, hyper2, center_, indices)


def _normalizer_gauge_group(matrix):
    """[summary]

    Args:
        matrix ([type]): [description]

    Returns:
        [type]: [description]
    """
    center_, hyper1, hyper2 = _symplectic_gram_schmidt(matrix, [], [])

    return _normalizer_gauge_group_preserve(center_, hyper1, hyper2)


# ---------------------------------------------------------------
# ---------------------------------------------------------------
#
# Methods below have been replaced with
#
# make_commute_hyper and _make_commute_hyper
#
# TODO: Change code to reflect these changes
# ---------------------------------------------------------------
# ---------------------------------------------------------------


def make_element_commute_with_hyper_pair(
    vector: ArrayLike, hyper1: ArrayLike, hyper2: ArrayLike
) -> np.ndarray:
    """Makes an element commute with a hyperbolic pair

    A new vector new_vector is formed from the input vector such that
    this new vector has a zero symplectic inner product with the two
    input hyperbolic vectors such that

    span(vector hyper1, hyper2) = span(new_vector, hyper1, hyper2)

    In terms of Pauli operators, a new Pauli operator is created such
    that it commutes with the input hyperbolic operator pair such that

    <op_vector, op_hyper1, op_hyper2> = <op_new_vector, op_hyper1, op_hyper2>

    Args:
        vector: Symplectic vector encoding a Pauli
        hyper1, hyper2: Symplectic vectors encoding a hyperbolic pair of vectors

    Returns:
        out : A vector that has zero symplectic product with the given hyperbolic pair

    Raises:
        QiskitError: All inputs must be one dimensional
        QiskitError: Input vectors must be GF(2) symplectic vectors

    Examples:
        >>> m = np.array([1,1,1,0,0,1], dtype=np.bool_)
        >>> h1 = np.array([0,0,1,0,0,0], dtype=np.bool_)
        >>> h2 = np.array([0,0,0,0,0,1], dtype=np.bool_)
        >>> x = make_element_commute_with_hyper_pair(test_matrix[1], test_matrix[2], test_matrix[3])
        >>> x.astype(int)
        array([1,1,0,0,0,0])

    See Also:
    _make_element_commute_with_hyper_pair, make_element_commute_with_hyper_pairs,
    _make_element_commute_with_hyper_pairs, make_elements_commute_with_hyper_pair,
    _make_elements_commute_with_hyper_pair, make_elements_commute_with_hyper_pairs,
    _make_elements_commute_with_hyper_pairs

    Notes:
    The this and the associated other functions should be made into a single functon
    list make_commute_hyper(matrix, hyper1, hyper2, mrange=None, h1range=None, h2range=None)
    or something like this
    """
    vector = np.array(vector)
    hyper1 = np.array(hyper1)
    hyper2 = np.array(hyper2)
    if not (
        is_symplectic_vector_form(vector)
        and is_symplectic_vector_form(hyper1)
        and is_symplectic_vector_form(hyper2)
    ):
        raise QiskitError(f"Input vectors must be GF(2) symplectic vectors")

    return _make_element_commute_with_hyper_pair(vector, hyper1, hyper2)


def _make_element_commute_with_hyper_pair(
    vector: np.ndarray, hyper1: np.ndarray, hyper2: np.ndarray
) -> np.ndarray:
    """Modify the operator given by vector so that it commutes
    with the hyperbolic pair (hyper1,hyper2) such that this new vector
    is in <op(vector), op(hyper1), op(hyper2)>

    Args:
        vector (numpy.ndarray dtype=bool): Symplectic vector encoding a Pauli operator
        hyper1 (numpy.ndarray dtype=bool): Symplectic vector encoding one of
            the operators in the hyperbolic elements
        hyper2 (numpy.ndarray dtype=bool): Symplectic vector encoding other one of
            the operators in the hyperbolic elements

    Returns:
        numpy.ndarray dtype=bool : Symplectic vector encoding a Pauli operator that commutes
            with the hyperbolic pair and is in <op(vector), op(hyper1), op(hyper2)>
    """
    new_vector = vector.copy()
    num_qubits = hyper1.shape[0] >> 1
    if _symplectic_product_vv(new_vector, hyper1, num_qubits):
        new_vector = new_vector ^ hyper2
    if _symplectic_product_vv(new_vector, hyper2, num_qubits):
        new_vector = new_vector ^ hyper1
    return new_vector


def make_element_commute_with_hyper_pairs(
    vector: ArrayLike, hyper1: ArrayLike, hyper2: ArrayLike, range1: ArrayLike, range2: ArrayLike
) -> np.ndarray:
    """Modify the operator elem so that it commutes with all the hyperbolic pairs given.

    Args:
        vector (numpy.ndarray dtype=bool): Symplectic vector encoding a Pauli operator
        hyper1 (numpy.ndarray dtype=bool): Array of symplectic vectors encoding one element of
            the operators in the hyperbolic pairs. Vectors to be used will be in range1
        hyper2 (numpy.ndarray dtype=bool): Array of symplectic vectors encoding the other element of
            the operators in the hyperbolic pairs. Vectors to be used will be in range2
        range1 (list,...): Indices for hyper1 to use
        range2 (list,...): Indices for hyper2 to use

    Returns:
        numpy.ndarray dtype=bool: Operator that commutes with provided hyperbolic pairs that
        is derived from input operator/vector
    """
    # Note: This may better be called to_isotropic(op, hyper1, hyper2)
    vector = np.array(vector)
    hyper1 = np.array(hyper1)
    hyper2 = np.array(hyper2)
    assert vector.ndim == 1 and hyper1.ndim > 1 and hyper2.ndim > 1, QiskitError(
        f"All inputs must be one dimensional: {vector.ndim},{hyper1.ndim}, {hyper2.ndim}"
    )
    assert vector.shape[0] == hyper1.shape[1] == hyper2.shape[1], QiskitError(
        "Inputs matrices/vectors must have the same number of columns/length"
    )
    assert not (vector.shape[0] % 2 or hyper1.shape[1] % 2 or hyper2.shape[1] % 2), QiskitError(
        f"All vectors must have an even length: \
        {vector.shape[0]},{hyper1.shape[0]},{hyper2.shape[0]}"
    )
    try:
        range1 = list(range1)
    except TypeError:
        QiskitError("Input range1 not iterable")
    try:
        range2 = list(range2)
    except TypeError:
        QiskitError("Input range2 not iterable")

    assert set(range1).issubset(range(hyper1.shape[0])) and set(range2).issubset(
        range(hyper2.shape[0])
    ), QiskitError(f"Input ranges not valid per hyper1 and hyper2 inputs ranges")

    return _make_element_commute_with_hyper_pairs(vector, hyper1, hyper2, range1, range2)


def _make_element_commute_with_hyper_pairs(
    vector: np.ndarray, hyper1: np.ndarray, hyper2: np.ndarray, range1: List, range2: List
) -> np.ndarray:
    """Modify the operator elem so that it commutes with all the hyperbolic pairs given.

    Args:
        vector (numpy.ndarray dtype=bool): Symplectic vector encoding a Pauli operator
        hyper1 (numpy.ndarray dtype=bool): Array of symplectic vectors encoding one element of
            the operators in the hyperbolic pairs. Vectors to be used will be in range1
        hyper2 (numpy.ndarray dtype=bool): Array of symplectic vectors encoding the other element of
            the operators in the hyperbolic pairs. Vectors to be used will be in range2
        range1 (list): Indices for hyper1 to use
        range2 (list): Indices for hyper2 to use

    Returns:
        numpy.ndarray dtype=bool: Operator that commutes with provided hyperbolic pairs that
        is derived from input operator/vector
    """
    new_vector = vector.copy()

    for index_1, index_2 in zip(range1, range2):
        new_vector = _make_element_commute_with_hyper_pair(
            new_vector, hyper1[index_1], hyper2[index_2]
        )

    return new_vector


# ---------------------------------------------------------------


def make_elements_commute_with_hyper_pair(
    matrix: ArrayLike,
    mrange: Union[range, List[int], np.ndarray],
    hyper1: ArrayLike,
    hyper2: ArrayLike,
):
    """Make the operators given by the range in the input matrix commute
    with the hyperbolic pair (hyper1, hyper2).

    Args:
        matrix (ArrayLike): Symplectic matrix (check). Source of operators to make commute
        mrange (Union[range, List[int], np.ndarray]): Range inicating which operators (rows) of the
            input matrix to make commute with hyperbolic pair
        hyper1 (ArrayLike): First part of hyperbolic pair
        hyper2 (ArrayLike): Second part of hyperbolic pair

    Returns:
        matrix: Modifed inoput matrix with operatiors (rows) given in range now computing with provided
        hyperbolic pair.
    """

    matrix = np.array(matrix, dtype=np.bool_)
    try:
        mrange = list(mrange)
    except TypeError:
        QiskitError("Input range1 not iterable")

    hyper1 = np.array(hyper1)
    hyper2 = np.array(hyper2)

    assert matrix.ndim == 2 and hyper1.ndim == 1 and hyper2.ndim == 1, QiskitError(
        f"input matrix must be 2 dimensional and hyper1, \
        hyper2 must be one dimernsional: {matrix.ndim}(2), {hyper1.ndim}(1), {hyper2.ndim}(1)"
    )
    assert not (matrix.shape[1] % 2 or hyper1.shape[0] % 2 or hyper2.shape[0] % 2), QiskitError(
        f"Inputs must have an even number of columns/length:\
        {matrix.shape[1]},{hyper1.shape[0]},{hyper2.shape[0]}"
    )
    assert matrix.shape[1] == hyper1.shape[0] == hyper2.shape[0], QiskitError(
        f"Input matices/vectors must have the same number of columns/length"
    )

    assert set(mrange).issubset(range(matrix.shape[1])), QiskitError(
        f"Input range not a valid range for input matrix: {mrange}"
    )

    return _make_elements_commute_with_hyper_pair(matrix, mrange, hyper1, hyper2)


def _make_elements_commute_with_hyper_pair(
    matrix: np.ndarray, mrange: List, hyper1: np.ndarray, hyper2: np.ndarray
):
    """Make the operators given by the range in the input matrix commute
    with the hyperbolic pair (hyper1, hyper2).

    Args:
        matrix (np.ndarray): Symplectic matrix (check). Source of operators to make commute
        range (List): Range inicating which operators (rows) of the
            input matrix to make commute with hyperbolic pair
        hyper1 (np.ndarray): First part of hyperbolic pair
        hyper2 (np.ndarray): Second part of hyperbolic pair

    Returns:
        matrix: Modifed inoput matrix with operatiors (rows) given in range now computing with provided
        hyperbolic pair.
    """

    new_matrix = matrix.copy()

    for index in mrange:
        new_matrix[index] = _make_element_commute_with_hyper_pair(new_matrix[index], hyper1, hyper2)

    return new_matrix


def make_elements_commute_with_hyper_pairs(in_matrix, mrange, hyper1, h1range, hyper2, h2range):
    """[summary]

    Args:
        matrix ([type]): [description]
        mrange ([type]): [description]
        hyper1 ([type]): [description]
        h1range ([type]): [description]
        hyper2 ([type]): [description]
        h2range ([type]): [description]

    Returns:
        [type]: [description]

    """

    return _make_elements_commute_with_hyper_pairs(
        in_matrix, mrange, hyper1, h1range, hyper2, h2range
    )


def _make_elements_commute_with_hyper_pairs(in_matrix, mrange, hyper1, h1range, hyper2, h2range):
    matrix = in_matrix.copy()

    for i, j in zip(h1range, h2range):
        matrix = _make_elements_commute_with_hyper_pair(matrix, mrange, hyper1[i], hyper2[j])

    return matrix
