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

"""Symplectic functions."""

from collections import deque
from typing import Any, List, Optional, Tuple, Union

import numpy as np
from qiskit import QiskitError

from qiskit_qec.linear import matrix as mt


# pylint: disable=invalid-name
def all_commute(matrix: np.ndarray) -> bool:
    r"""Determines if each possible pair of different rows of the
    GF(2) symplectic matrix have zero symplectic product. If the rows represent
    Pauli operators then the this method determines if the operators
    defined by the matrix generate an abelian subgroup.

    Args:
        matrix: Input GF(2) symplectic matrix

    Returns:
        True if operators mutually commute - have zero symplectic product

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
    if matrix.shape[0] == 1:
        return True
    test_mat = np.asarray(symplectic_product(matrix, matrix))
    return not test_mat.any()


# ---------------------------------------------------------------


def symplectic_product(mat1: np.ndarray, mat2: np.ndarray) -> Union[int, np.ndarray]:
    r"""Returns the symplectic product of two GF(2) symplectic matrices.

    Let math:'A', math:'B' be two GF(2) symplectic matrices of width math:'2m',
    then the symplectic product is defined as

    .. math::
        (A,B) = A \cdot \Lambda_n \cdot B^T.

    Args:
        mat1, mat2: Input GF(2) symplectic matrixes

    Returns:
        Symplectic product of mat1 and mat2

    Raises:
        QiskitError: Input matrices/vectors must be GF(2) symplectic matrices/vectors
        QiskitError: Input matrices must have the same number of dimensions
        QiskitError: Input matrices must be 1 or 2 dimensional

    Examples:
        >>> mat1 = numpy.array([[1,0,1,1],[0,1,1,0]], dtype=numpy.bool_)
        >>> mat2 = numpy.array([[0,1,1,1],[0,0,1,0]], dtype=numpy.bool_)
        >>> p = symplectic_product(mat1, mat2)
        >>> p.astype(int)
        array([[0, 1], [1, 0]]

    See Also:
    _symplectic_product_vv, _symplectic_product_dense
    """
    mat1_np_array = np.array(mat1, dtype=np.int8)
    mat2_np_array = np.array(mat2, dtype=np.int8)

    if not is_symplectic_form(mat1) or not is_symplectic_form(mat2):
        raise QiskitError("Input matrices/vectors must be GF(2) symplectic matrices/vectors")

    if not mat1_np_array.ndim == mat2_np_array.ndim:
        raise QiskitError(
            f"Input matrices must have the \
            same dimensions: {mat1_np_array.ndim} is \
            not equal to {mat2_np_array.ndim}"
        )

    if mat1_np_array.ndim == 1:
        if not mat1_np_array.shape[0] == mat2_np_array.shape[0]:
            raise QiskitError(
                f"Input vectors must have the same \
            dimensions: {mat1_np_array.shape[0]} not equal \
            to {mat2_np_array.shape[0]}"
            )
        if mat1_np_array.dtype == bool:
            return _symplectic_product_vv_boolean(
                mat1_np_array, mat2_np_array, mat1_np_array.shape[0] >> 1
            )
        return _symplectic_product_vv(mat1_np_array, mat2_np_array, mat1_np_array.shape[0] >> 1)

    elif mat1_np_array.ndim == 2:
        return _symplectic_product_dense(mat1_np_array, mat2_np_array)

    else:
        raise QiskitError(
            f"Input matrices must be 1 or 2 dimensional:\
            {mat1_np_array.ndim}, {mat2_np_array.ndim}"
        )


def _symplectic_product_vv(vec1: np.ndarray, vec2: np.ndarray, n: int) -> int:
    r"""Finds the sympletic product or two GF(2) symplectic vectors of
    length 2n: vec1 . Lambda . vec2^T where

    lambda = [0 I]
             [I 0]

    Warning: This method requires integer components. Python bool and nump bool_
    type will give incorrect results. Use _symplectic_product_vv_boolean  or the
    more general symplectic_product

    Args:
        vec1, vec2: Input GF(2) symplectic vectors
        n: Input size, half of the symplectic vector length

    Returns:
        Symplectic product of vec1 and vec2

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
    assert vec1.dtype != bool
    r = 0
    for i in range(n):
        r += vec1[i] * vec2[n + i] + vec1[n + i] * vec2[i]
    return r % 2


def _symplectic_product_vv_boolean(vec1: np.ndarray, vec2: np.ndarray, n: int) -> int:
    r"""Finds the sympletic product or two GF(2) (boolean types) symplectic vectors of
    length 2n: vec1 . Lambda . vec2^T where

    lambda = [0 I]
             [I 0]

    Warning: This method requires boolean components.

    Args:
        vec1, vec2: Input GF(2) symplectic vectors
        n: Input size, half of the symplectic vector length

    Returns:
        out: Symplectic product of vec1 and vec2

    Examples:
        >>> a = np.array([1,0,0,0,1,1,1,0,1,0], dtype=numpy.bool_)
        >>> b = numpy.array([1,1,1,1,0,0,1,0,1,1], dtype=numpy.bool_)
        >>> _symplectic_product_vv_boolean(a, b, 5)
        0

    See Also:
    symplectic_product, _symplectic_product_dense, _symplectic_product_vv

    Notes:
    If changing the method please make sure that the new method
    is faster. Note that this method is faster if the vectors are
    numpy arrays with dtype=int8
    """
    assert vec1.dtype == bool
    r = False
    for i in range(n):
        r += vec1[i] & vec2[n + i] ^ vec1[n + i] & vec2[i]
    return int(r)


def _symplectic_product_dense(mat1: np.ndarray, mat2: np.ndarray) -> Union[int, np.ndarray]:
    r"""Returns the symplectic product of two GF(2) symplectic matrices.

    Let math:'A', math:'B' be two GF(2) (as integer not boolean) symplectic matrices of width math:'2m',
    then the symplectic product is defined as

    .. math::
        (A,B) = A \cdot \Lambda_n \cdot B^T.

    Args:
        mat1, mat2: Input GF(2) symplectic matrixes

    Returns:
        Symplectic product of mat1 and mat2

    Examples:
        >>> mat1 = numpy.array([[1,0,1,1],[0,1,1,0]], dtype=numpy.bool_)
        >>> mat2 = numpy.array([[0,1,1,1],[0,0,1,0]], dtype=numpy.bool_)
        >>> p = _symplectic_product_dense(mat1, mat2)
        >>> p.astype(int)
        array([[0, 1], [1, 0]]

    See Also:
    _symplectic_product_vv, symplectic_product
    """
    assert mat1.dtype != bool
    m1, m2 = np.hsplit(mat1, 2)  # pylint: disable=unbalanced-tuple-unpacking
    result = np.hstack((m2, m1)).dot(mat2.transpose()) % 2
    if result.size == 1:
        return int(result.item())
    else:
        return result


# ---------------------------------------------------------------


def make_commute_hyper(
    a: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    arange: Optional[np.ndarray] = None,
    xrange: Optional[np.ndarray] = None,
    zrange: Optional[np.ndarray] = None,
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
        x: GF(2) hyperbolic pair vector
        z: GF(2) hyperbolic pair vector
        arange (optional): range of indices from a to make commute. Defaults to None.
        xrange (optional): range of indices from x to use. Defaults to None.
        zrange (optional): range of indices from z to use. Defaults to None.

    Raises:
        QiskitError: Input matrices/vectors must bf GF(2) symplectic matrices/vectors
        QiskitError: Input range is not iterable")
        QiskitError: Input matrices/vectors must have the same number of columns/length

    Returns:
        GF(2) symplectic vectors that commute with the given hyperbolic pairs

    Examples:
        >>> a = numpy.array([1,1,1,0,0,0],dtype=numpy.bool_)
        >>> x = numpy.array([0,0,1,0,0,0],dtype=numpy.bool_)
        >>> z = numpy.array([0,0,0,0,0,1],dtype=numpy.bool_)
        >>> a = make_commute_hyper(a, x, z)
        >>> a.astype(int)
        array([1, 1, 0, 0, 0, 0])

        >>> a = numpy.array([1,1,1,0,0,0,0,0], dtype=numpy.bool_)
        >>> x = numpy.array([[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0]], dtype=numpy.bool_)
        >>> xrange = [0,1]
        >>> zrange = [0,1]
        >>> a = make_commute_hyper(a, x, z, xrange = xrange, zrange=zrange)
        >>> a.astype(int)
        array([1, 0, 0, 0, 0, 0, 0, 0])

        >>> a = numpy.array([[1,1,1,0,0,0,0,0],[0,1,1,0,0,0,0,0]], dtype=numpy.bool_) # X1X2X3, X2X3
        >>> x = numpy.array([0,1,0,0,0,0,0,0], dtype=numpy.bool_) # X2
        >>> z = numpy.array([0,0,0,0,0,1,0,0], dtype=numpy.bool_) # Z2
        >>> arange = [0,1]
        >>> a = make_commute_hyper(a, x, z, arange)
        >>> a.astype(int)
        array([[1, 0, 1, 0, 0, 0, 0, 0],
               [0, 0, 1, 0, 0, 0, 0, 0]])

        >>> a = numpy.array([[1,1,1,0,0,0,0,0],[0,1,1,1,0,0,0,0]], dtype=numpy.bool_)
        >>> x = numpy.array([[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0,0,0,0,0,1,0,0], [0,0,0,0,0,0,1,0]], dtype=numpy.bool_)
        >>> arange = [0,1]
        >>> a = make_commute_hyper(a, x, z, arange)
        >>> a.astype(int)
        array([[1, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0]])

    See Also:
    _make_commute_hyper
    """
    if not (is_symplectic_form(a) and is_symplectic_form(x) and is_symplectic_form(z)):
        raise QiskitError("Input matrices/vectors must be GF(2) symplectic matrices/vectors")

    def make_list(srange):
        if srange is not None:
            try:
                srange = list(srange)
            except TypeError as terror:
                raise QiskitError(f"Input range {srange} is not iterable") from terror

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

    if not a.shape[1] == x.shape[1] == z.shape[1]:
        raise QiskitError("Input matrices/vectors must have the same number of columns/length")

    return _make_commute_hyper(a, x, z, arange, xrange, zrange, squeeze)


# ---------------------------------------------------------------


def _make_commute_hyper(
    a: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    arange: Optional[np.ndarray] = None,
    xrange: Optional[np.ndarray] = None,
    zrange: Optional[np.ndarray] = None,
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
        GF(2) symplectic vectors that commute with the given hyperbolic pairs

    Examples:
        >>> a = numpy.array([[1,1,1,0,0,0]],dtype=numpy.bool_)
        >>> x = numpy.array([[0,0,1,0,0,0]],dtype=numpy.bool_)
        >>> z = numpy.array([[0,0,0,0,0,1]],dtype=numpy.bool_)
        >>> a = _make_commute_hyper(a, x, z, squeeze=True)
        >>> a.astype(int)
        array([1, 1, 0, 0, 0, 0])

        >>> a = numpy.array([[1,1,1,0,0,0,0,0]], dtype=numpy.bool_)
        >>> x = numpy.array([[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0]], dtype=numpy.bool_)
        >>> xrange = [0,1]
        >>> zrange = [0,1]
        >>> a = _make_commute_hyper(a, x, z, xrange = xrange, zrange=zrange, squeeze=False)
        >>> a.astype(int)
        array([[1, 0, 0, 0, 0, 0, 0, 0]])

        >>> a = numpy.array([[1,1,1,0,0,0,0,0],[0,1,1,0,0,0,0,0]], dtype=numpy.bool_)
        >>> x = numpy.array([[0,1,0,0,0,0,0,0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0,0,0,0,0,1,0,0]], dtype=numpy.bool_)
        >>> arange = [0,1]
        >>> a = _make_commute_hyper(a, x, z, arange)
        >>> a.astype(int)
        array([[1, 0, 1, 0, 0, 0, 0, 0],
               [0, 0, 1, 0, 0, 0, 0, 0]])

        >>> a = numpy.array([[1,1,1,0,0,0,0,0],[0,1,1,1,0,0,0,0]], dtype=numpy.bool_)
        >>> x = numpy.array([[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0]], dtype=numpy.bool_)
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
            if _symplectic_product_vv(a[k].astype(int), x[i].astype(int), num_qubits):
                a[k] = a[k] ^ z[j]
            if _symplectic_product_vv(a[k].astype(int), z[j].astype(int), num_qubits):
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
        Tuple of the found hyperbolic partner (av) and its index in the
            search matrix. 'None' if no hyperbolic partner is found.

    Examples:
        >>> matrix = numpy.array([[1,0,1,0,0,0,0,0],[0,1,1,0,0,0,0,0]], dtype=numpy.bool_)
        >>> vector = numpy.array([0,0,0,0,0,1,0,0], dtype=numpy.bool_)
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
        Tuple of the found hyperbolic partner (av) and its index in the
            search matrix. 'None' if no hyperbolic partner is found.

    Examples:
        >>> matrix = numpy.array([[1,0,1,0,0,0,0,0],[0,1,1,0,0,0,0,0]], dtype=numpy.bool_)
        >>> vector = numpy.array([0,0,0,0,0,1,0,0], dtype=numpy.bool_)
        >>> av, index = _locate_hyper_partner(matrix, vector)
        >>> av.astype(int)
        array([0, 1, 1, 0, 0, 0, 0, 0])
        >>> index
        1

    See Also:
    locate_hyper_partner, build_hyper_partner, _build_hyper_partner
    """
    n = matrix.shape[1] >> 1
    for index, item in enumerate(matrix):
        if _symplectic_product_vv(item.astype(int), vector.astype(int), n) == 1:
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
        matrix (np.array, SimplecticMatrix): GF(2) symplectic matrix representing a set of independent
            commuting generators
        index: index of generator to build a hyperbolic partner for

    Raises:
        QiskitError: Input matrix must be a GF(2) symplectic matrix
        QiskitError: Input matrix must represent a set of commuting operators
        QiskitError: Input matrix does not represent a set of independent
            operators, it does not have have full rank
        QiskitError: Input index out or range

    Returns:
        a hyperbolic partner for the given vector wrt the set of commuting
            generators

    Examples:
        >>> matrix = numpy.array(
            [[1,0,0,0,0,0,0,0],
             [0,1,0,0,0,0,0,0],
             [0,0,1,0,0,0,0,0],
             [0,0,0,1,0,0,0,0]], dtype=numpy.bool_)
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
    if not all_commute(matrix):
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
        a hyperbolic partner for the given vector wrt the set of commuting
            generators

    Examples:
        >>> matrix = numpy.array(
            [[1,0,0,0,0,0,0,0],
             [0,1,0,0,0,0,0,0],
             [0,0,1,0,0,0,0,0],
             [0,0,0,1,0,0,0,0]], dtype=numpy.bool_)
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

    heads, _, transform_mat, _ = mt._rref_complete(slambda)

    e_index = np.zeros(nrows, dtype=bool)
    e_index[index] = True

    trans_e_index = np.matmul(transform_mat, e_index)

    pivot = 0
    result = np.zeros(ncols, dtype=bool)
    for i in range(ncols):
        if heads[i] == 1:
            result[i] = trans_e_index[pivot]
            pivot += 1

    return result


# ---------------------------------------------------------------


def symplectic_gram_schmidt(
    a: np.ndarray, x: Optional[np.ndarray] = None, z: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Applies the sympletic Gram-Schmidt process to the input matrix

    Apply the symplectic GramSchmidt process to the input symplectic matrix. Resulting
    hyperbolic pairs are added to x and z arrays. Elements of the center will be added to the
    center array.

    Args:
        a: Symplectic matrix
        x (optional): GF(2) Symplectic matrices representing hyperbolic pairs to
        z (optional): GF(2) Symplectic matrices representing hyperbolic pairs to
        build upon. Default is None.

    Raises:
        QiskitError Input matric not a GF(2) symplectic matrix
        QiskitError: Input hyperbolic array x is not a GF(2) sympletic matrix
        QiskitError: Input hyperbolic array z is not a GF(2) sympletic matrix
        QiskitError: Input hyperbolic arrays have different dimensions
        QiskitError: Input hyperbolic matrices do not represent a hyperbolic basis

    Returns:
        Center array and hyperbolic pairs split accross x and z

    Examples:
        >>> a = numpy.array([[0,1,0,0,1,0,1,0],
                          [0,0,0,0,1,1,0,1],
                          [1,1,1,0,0,1,0,0],
                          [1,1,0,1,0,0,0,0]], dtype=numpy.bool_)
        >>> center_, x, z = symplectic_gram_schmidt(a)
        >>> center_.astype(int)
        array([[1, 1, 1, 0, 1, 0, 0, 1],
               [1, 0, 0, 1, 0, 1, 1, 1]])
       >>> x.astype(int)
       array([[0, 1, 0, 0, 1, 0, 1, 0]])
       >>> z.astype(int)
       array([[0, 0, 0, 0, 1, 1, 0, 1]])

    Also See:
    _symplectic_gram_schmidt

    TODO: Add an example that shows using the optional x and z arrays
    """

    a = np.atleast_2d(np.array(a))
    if not is_symplectic_matrix_form(a):
        raise QiskitError("Input matrix not a GF(2) symplectic matrix")
    if x is None:
        x = []
    else:
        x = np.atleast_2d(x)
        x = list(x)
        if not is_symplectic_vector_form(x[0]):
            raise QiskitError("Input hyperbolic array x is not a GF(2) sympletic matrix")

    if z is None:
        z = []
    else:
        z = np.atleast_2d(z)
        z = list(z)
        if not is_symplectic_vector_form(z[0]):
            raise QiskitError("Input hyperbolic array z is not a GF(2) sympletic matrix")

    if not len(x) == len(z):
        raise QiskitError("Input hyperbolic arrays have different dimensions")

    if len(x) > 0 and x[0].shape[0] != z[0].shape[0]:
        raise QiskitError("Input hyperbolic arrays have different dimensions")

    if x != []:
        if not is_hyper_form(x, z):
            raise QiskitError("Input hyperbolic matrices do not represent a hyperbolic basis")

    return _symplectic_gram_schmidt(a, x, z)


def _symplectic_gram_schmidt(
    a: np.ndarray, x: List[np.ndarray], z: List[np.ndarray]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
        Center array and hyperbolic pairs split accross x and z

    Examples:
        >>> a = numpy.array([[0,1,0,0,1,0,1,0],
                          [0,0,0,0,1,1,0,1],
                          [1,1,1,0,0,1,0,0],
                          [1,1,0,1,0,0,0,0]], dtype=numpy.bool_)
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
            if elem.any():
                center_.append(elem)

    x = np.asarray(x)
    z = np.asarray(z)
    if len(center_) == 0:
        center_ = np.zeros(shape=(1, x.shape[1]), dtype=np.bool_)
    else:
        center_ = np.asarray(center_)

    return center_, x, z


# ---------------------------------------------------------------


def count_num_y(matrix: np.ndarray, scalar: bool = True) -> Union[np.ndarray, int]:
    """Returns the number of positions with 1's in k and n+k positions
    for matrices/vectors of width 2n for all k

    Args:
        matrix: Input GF(2) symplectic matrix/vector
        scalar: If scalar is True and a vector is input then a scalar will
            be output else a vector will be output. Default is True

    Raises:
        QiskitError: Input matrix/vector not a GF(2) symplectic matrix

    Returns:
        result: number of positions with 1's in k and n+k positions
    for matrices/vectors of width 2n for all k.

    Examples:
        >>> a = np.array([1,0,1,1], dtype=np.bool_)
        >>> count_num_y(a)
        1

        >>> a = np.array([1,0,1,1], dtype=np.bool_)
        >>> count_num_y(a, scalar=False)
        array([1])

        >>> b = np.array([[1,0,1,1],[0,0,1,1]], dtype=np.bool_)
        >>> count_num_y(b)
        array([1,0])

    """
    matrix = np.atleast_2d(matrix)
    if not is_symplectic_form(matrix):
        raise QiskitError("Input matrix/vector not a GF(2) symplectic matrix")
    num_qubits = matrix.shape[1] >> 1
    result = _count_num_y(matrix, num_qubits)
    if scalar and matrix.shape[0] == 1:
        return result[0]
    else:
        return result


def _count_num_y(matrix: np.ndarray, n: int) -> np.ndarray:
    """Returns the number of positions with 1's in k and n+k positions
    for matrices/vectors of width 2n for all k

    Args:
        matrix (np.ndarray): Input GF(2) symplectic matrix/vector
        n (int): half of the number of columns of input matrix

    Returns:
        np.ndarray: number of positions with 1's in k and n+k positions
    for matrices/vectors of width 2n for all k.

    Examples:
        >>> a = np.array([1,0,1,1], dtype=np.bool_)
        >>> _count_num_y(a)
        array([1])

        >>> b = np.array([[1,0,1,1],[0,0,1,1]], dtype=np.bool_)
        >>> count_num_y(a, scalar=False)
        array([1,0])
    """
    return np.sum(np.logical_and(matrix[:, :n], matrix[:, n:]), axis=1, dtype=int)


# ---------------------------------------------------------------


def is_symplectic_matrix_form(
    matrix: np.ndarray, dtype: Optional[Union[bool, np.bool_, int, np.integer]] = None
) -> bool:
    """Is the input matrix GF(2) symplectic

    Checks if the given array-like matrix is in the form of a symplectic matrix:
    two dimensional, even number of columns, 0/1 or boolean entries. The optional
    argument dtype can be given to check if entries are a specific dtype

    Args:
        matrix: Input matrix to be checked
        dtype: Optional. Check if given matrix is of type dtype. Default: None

    Returns:
        True if the input matrix is GF(2) symplectic. False otherwise.

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
    vector: np.ndarray, dtype: Optional[Union[bool, np.bool_, int, np.integer]] = None
) -> bool:
    """Is the input vector GF(2) symplectic

    Checks if the given array like vector is in the form of a symplectic vector:
    two dimensional, even number of columns, 0/1 or boolean entries. The optional
    argument dtype can be given to check if entries are a specific dtype

    Args:
        vector: Input vector to be checked
        dtype: Optional. Check if given vector is of type dtype. Default: None

    Returns:
        True if the input vector is GF(2) symplectic. False otherwise.

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
        True if the input matrix/vector is GF(2) symplectic. False otherwise.

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


def is_center(cntr: np.ndarray, matrix: np.ndarray) -> bool:
    """Does the input center matrix represent the center of the supplied matrix?

    Let op(cntr) = [op(c_0),op(c_1),...,op(c_(k-1))] be the Pauli operators
    represented by the input center_matrix. Let op(matrix) = [op(m_0),op(m_1),...,op(m_(t-1))]
    be the Pauli operators represented by the input matrix. This method
    returns True if <op(cntr)> = Z(<op(matrix)>) and False otherwise.

    Args:
        cntr (np.ndarray): Generators of center to be checked
        matrix (np.ndarray): Generators of full group

    Raises:
        QiskitError: Not all inputs are not GF(2) symplectic matrices/vectors

    Returns:
        True if <op(cntr)> = Z(<op(matrix)>), False otherwise

    Examples:
        >>> matrix = numpy.array(
            [[0,1,0,0,1,0,1,0],
            [0,0,0,0,1,1,0,1],
            [1,1,1,0,0,1,0,0],
            [1,1,0,1,0,0,0,0]], dtype=numpy.bool_)
        >>> cntr = numpy.array([[1, 1, 1, 0, 1, 0, 0, 1],
                                [1, 0, 0, 1, 0, 1, 1, 1]], dtype=numpy.bool_)
        >>> is_center(cntr, matrix)
        True
    """
    matrix = np.atleast_2d(np.array(matrix))
    cntr = np.atleast_2d(np.array(cntr))
    if not (is_symplectic_matrix_form(cntr) and is_symplectic_matrix_form(matrix)):
        raise QiskitError("Not all inputs are not GF(2) symplectic matrices")
    cal_center = center(matrix)

    return is_same_span(cntr, cal_center)


# ---------------------------------------------------------------


def is_same_span(matrix1: np.ndarray, matrix2: np.ndarray) -> bool:
    """Does span(rows of matrix1) = span(rows of matrix2)?

    Args:
        matrix1: First set of vectors
        matrix2: Second set of vectors

    Returns:
        True if span(rows of matrix1) = span(rows of matrix2). False otherwise

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
        raise QiskitError("Inpiut matrices must by GF(2) symplectic matrices")

    if matrix1.shape[1] != matrix2.shape[1]:
        return False

    _, rref_matrix1, _, rank_matrix1 = mt.rref_complete(matrix1)
    _, rref_matrix2, _, rank_matrix2 = mt.rref_complete(matrix2)

    if rank_matrix1 != rank_matrix2:
        return False

    rref_matrix1 = rref_matrix1[:rank_matrix1]
    rref_matrix2 = rref_matrix2[:rank_matrix2]
    return np.array_equal(rref_matrix1, rref_matrix2)


def is_hyper_form(x: Union[list, np.ndarray], z: Union[list, np.ndarray]) -> bool:
    """Do the input matrices form a hyperbolic/symplectic basis?

    Args:
        x,z: Pairs to test if they form a hyperbolic/symnplectic basis

    Returns:
        True if input matrices form a hyperbolic/symplectic basis

    Examples:
        >>> x = numpy.array([[1,0,0,0,],[0,1,0,0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0,0,1,0,],[0,0,0,1]], dtype=numpy.bool_)
        >>> is_hyper_form(x,z)
        True
    """
    if isinstance(x, list):
        x = np.array(x)
    if isinstance(z, list):
        z = np.array(z)

    # Check for empty x and z: Null comdition -> True
    if x.shape[0] == 0 and z.shape[0] == 0:
        return True

    matrix = np.vstack((x, z))
    test = symplectic_product(matrix, matrix)
    return np.array_equal(test, mt._create_lambda_matrix(matrix.shape[0] >> 1))


def is_stabilizer_group(matrix: np.ndarray) -> bool:
    """Do the rows of the input matrix represent the generators of an
    abelian Pauli subgroup?

    Args:
        matrix: GF(2) symplectic matrix

    Raises:
        QiskitError: Input matrix not a GF(2) symplectic matrix

    Returns:
        True is matrix represents a Stabilizer group

    Examples:
        >>> matrix = numpy.array([[1,0,0,0,0,0],
                               [0,1,0,0,0,0],
                               [1,0,0,0,0,1]], dtype=numpy.bool_)
        >>> is_stabilizer_group(matrix)
        True
    """
    matrix = np.atleast_2d(np.array(matrix))
    if not is_symplectic_matrix_form(matrix):
        raise QiskitError("Input matrix not a GF(2) symplectic matrix")
    return all_commute(matrix)


# ---------------------------------------------------------------


def center(matrix: np.ndarray, preserve: bool = False) -> np.ndarray:
    """Find the center of the group with generators given by the symplectic matrix

    Args:
        matrix: GF(2) symplectic matrix
        preserve: If True then an attempt will be made to preserve then generators form

    Raises:
        QiskitError: Input matrix is not a symplectic matrix

    Returns:
        Generators for the center, represented as a symplectic matrix,
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
        Generators for the center, represented as a symplectic matrix,
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
    return _symplectic_gram_schmidt(matrix, [], [])[0]


def _center_preserve(matrix: np.ndarray) -> np.ndarray:
    """Find the center of the group with generators given by the symplectic matrix

    Args:
        matrix: GF(2) symplectic matrix

    Returns:
        Generators for the center, represented as a symplectic matrix,
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
    for opi in reversed(matrix):
        break_flag = False
        for opj in matrix:
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


def basis_for_pauli_group(matrix: np.ndarray) -> np.ndarray:
    """Given a set of generators (not necessarily independent) find
    a full basis using as many of the provided generators as possible.

    Args:
        matrix: Set of generators (in GF(2) symplectic form)

    Raises:
        QiskitError: Input matrix not in a GF(2) symplectic matrix

    Returns:
        A maximal independant set

    Examples:
        >>> matrix = numpy.array([[1,1,0,0,1,0],[0,0,0,1,1,0],[0,1,0,1,0,1]], dtype=numpy.bool_)
        >>> basis = basis_for_pauli_group(matrix)
        >>> basis.astype(int)
        array([[1, 1, 0, 0, 1, 0],
                [0, 0, 0, 1, 1, 0],
                [0, 1, 0, 1, 0, 1],
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0]])
    """
    matrix = np.atleast_2d(np.asarray(matrix))
    if not is_symplectic_matrix_form(matrix):
        raise QiskitError("Input matrix not in a GF(2) symplectic matrix")

    return _basis_for_pauli_group(matrix)


def _basis_for_pauli_group(matrix: np.ndarray) -> np.ndarray:
    """Given a set of generators (not necessarily independent) find
    a full basis using as many of the provided generators as possible.

    Args:
        matrix: Set of generators (in GF(2) symplectic form)

    Returns:
        A maximal independant set

    Examples:
        >>> matrix = numpy.array([[1,1,0,0,1,0],[0,0,0,1,1,0],[0,1,0,1,0,1]], dtype=numpy.bool_)
        >>> basis = _basis_for_pauli_group(matrix)
        >>> basis.astype(int)
        array([[1, 1, 0, 0, 1, 0],
                [0, 0, 0, 1, 1, 0],
                [0, 1, 0, 1, 0, 1],
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0]])
    """
    aug_matrix = mt.augment_mat(matrix, "bottom")
    heads, _, _, rank = mt._rref_complete(aug_matrix.T)
    shape = (rank, aug_matrix.shape[1])
    ext_matrix = np.zeros(shape, dtype=np.bool_)
    posns = np.flatnonzero(heads)
    for k, index in enumerate(posns):
        ext_matrix[k] = aug_matrix[index]
    return ext_matrix


def make_hyperbolic(
    center_: np.ndarray, x: np.ndarray, z: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Makes a isotropic hyperbolic basis into an hyperbolic basis

    Let center, x and z (where x, z and center are represented by GF(2)
    symplectic matrices/vectors) be a isotropic-hyperbolic basis with
    2m+k total generators. This method makes/extends this basis into
    a hyperbolic basis for the (m+k)-th Pauli group by constructing hyperbolic
    partners (acenter) for each element in the center such that

    P_(n+k) = <iI, x_new, z_new>

    where x_new = <x, center>  and z_new = <z, acenter>

    Args:
        center_: isotropic part of basis (center)
        x: hyperbolic part of basis
        z: hyperbolic part of basis

    Raises:
        QiskitError: Input center is not center of group generated by vectors from center, x and z
        QiskitError: Input matrices x, z are not in hyperbolic pairs

    Examples:
        >>> center_ = numpy.array([[1, 1, 1, 0, 1, 0, 0, 1],
                                [1, 0, 0, 1, 0, 1, 1, 1]], dtype=numpy.bool_)
        >>> x = numpy.array([[0, 1, 0, 0, 1, 0, 1, 0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0, 0, 0, 0, 1, 1, 0, 1]], dtype=numpy.bool_)
        >>> center_, x, z = symplectic_gram_schmidt(a, x, z)
        >>> x, z = make_hyperbolic(center_, x, z)
        >>> x
        array([[0, 1, 0, 0, 1, 0, 1, 0],
               [1, 0, 0, 1, 0, 1, 1, 1],
               [1, 1, 1, 0, 1, 0, 0, 1]])
        >>> z
        array([[0, 0, 0, 0, 1, 1, 0, 1],
               [0, 0, 0, 0, 1, 0, 1, 0],
               [0, 1, 0, 1, 0, 0, 0, 0]])

    See Also:
    _make_hyperbolic

    """
    if center_ is None:
        return x, z

    center_ = np.atleast_2d(np.array(center_))
    if center_.shape[0] == 1 and not center_.any():
        return x, z

    center_ = min_generating(center_)

    x = np.atleast_2d(np.array(x))
    z = np.atleast_2d(np.array(z))

    if not is_center(center_, np.vstack((x, z, center_))):
        raise QiskitError(
            "Input center is not center of group generated by vectors from center, x and z"
        )

    if not is_hyper_form(x, z):
        raise QiskitError("Input matrices x, z are not in hyperbolic pairs")

    return _make_hyperbolic(center_, x, z)


def _make_hyperbolic(
    center_: np.ndarray, x: np.ndarray, z: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Makes a isotropic hyperbolic basis into an hyperbolic basis

    Let center, x and z (where x, z and center are represented by GF(2)
    symplectic matrices/vectors) be a isotropic-hyperbolic basis with
    2m+k total generators. This method makes/extends this basis into
    a hyperbolic basis for the (m+k)-th Pauli group by constructing hyperbolic
    partners (acenter) for each element in the center such that

    P_(n+k) = <iI, x_new, z_new>

    where x_new = <x, center>  and z_new = <z, acenter>

    Args:
        center_: isotropic part of basis (center)
        x, z: hyperbolic part of basis

    Examples:
        >>> center_ = numpy.array([[1, 1, 1, 0, 1, 0, 0, 1],
                                [1, 0, 0, 1, 0, 1, 1, 1]], dtype=numpy.bool_)
        >>> x = numpy.array([[0, 1, 0, 0, 1, 0, 1, 0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0, 0, 0, 0, 1, 1, 0, 1]], dtype=numpy.bool_)
        >>> center_, x, z = symplectic_gram_schmidt(a, x, z)
        >>> x, z = _make_hyperbolic(center_, x, z)
        >>> x
        array([[0, 1, 0, 0, 1, 0, 1, 0],
               [1, 0, 0, 1, 0, 1, 1, 1],
               [1, 1, 1, 0, 1, 0, 0, 1]])
        >>> z
        array([[0, 0, 0, 0, 1, 1, 0, 1],
               [0, 0, 0, 0, 1, 0, 1, 0],
               [0, 1, 0, 1, 0, 0, 0, 0]])

    See Also:
    make_hyperbolic

    """

    center_size = center_.shape[0]
    hyper_size = x.shape[0]
    shape = (center_size, x.shape[1])
    spacer = np.zeros(shape, dtype=np.bool_)
    x = np.vstack((x, spacer))
    z = np.vstack((z, spacer))
    center_ = center_.copy()

    while center_size > 0:
        hop = _build_hyper_partner(center_[:center_size], center_size - 1)
        # TODO: Change the use of make_commute_hyper to _make_commute_hyper
        hop = make_commute_hyper(hop, x, z, xrange=range(hyper_size), zrange=range(hyper_size))
        # hop = _make_element_commute_with_hyper_pairs(
        #    hop, x, z, range(hyper_size), range(hyper_size)
        # )

        center_size -= 1
        x[hyper_size] = center_[center_size]
        z[hyper_size] = hop
        hyper_size += 1

    return x, z


def make_isotropic_hyperbolic_form(
    matrix: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Creates a isotrophic hyperbolic basis from a set of generators

    Args:
        matrix (np.ndarray): GF(2) symplectic matrix

    Returns:
        An isotrophic hyperbolic basis (center, x, z)

    Examples:
        >>> matrix = np.array([[0,1,0,0,1,0,1,0],
                               [0,0,0,0,1,1,0,1],
                               [1,1,1,0,0,1,0,0],
                               [1,1,0,1,0,0,0,0]], dtype=np.bool_)
        >>> center_, x, z = make_isotropic_hyperbolic_form(matrix)
        >>> center_.astype(int)
        array([[1, 1, 1, 0, 1, 0, 0, 1],
               [1, 0, 0, 1, 0, 1, 1, 1]])
        >>> x.astype(int)
        array([[0, 1, 0, 0, 1, 0, 1, 0]])
        >>> z.astype(int)
        array([[0, 0, 0, 0, 1, 1, 0, 1]])

    See Also:
    symplectic_gram_schmidt, _symplectic_gram_schmidt
    """
    return _symplectic_gram_schmidt(matrix, [], [])


def hyperbolic_basis_for_pauli_group(
    matrix: Optional[np.ndarray] = None,
    x: Optional[np.ndarray] = None,
    z: Optional[np.ndarray] = None,
    n: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Builds a hyperbolic basis for the associated Pauli group

    Args:
        matrix (Optional): Input GF(2) symplectic matrix
        x: (Optional): Input GF(2) hyperbolic pairs
        z (Optional): Input GF(2) hyperbolic pairs
        n: Size of standard Pauli group basis to create (if matrix,x,z are None)

    Raises:
        QiskitError: x and z must be both be None or both be array like
        QiskitError: x not a GF(2) symplectic matrix
        QiskitError: z not a symplectic matrix
        QiskitError: matrix not a GF(2) symplectic matrix
        QiskitError: When providing an input hyperbolic set the input matrix
            must be the center of the full set of generators
        QiskitError: If matrix, x and z are None then n must be provided

    Returns:
        Hyperbolic basis for the associated Pauli group (new_x, new_z)

    Examples:
        >>> matrix = numpy.array([[0,1,0,0,1,0,1,0],
                               [0,0,0,0,1,1,0,1],
                               [1,1,1,0,0,1,0,0],
                               [1,1,0,1,0,0,0,0]], dtype=numpy.bool_)
        >>> center_, x, z = make_isotropic_hyperbolic_form(matrix)
        >>> nx, nz = sysp.hyperbolic_basis_for_pauli_group(center_, x, z)
        >>> nx.astype(int)
        array([[0, 1, 0, 0, 1, 0, 1, 0],
               [1, 0, 0, 1, 0, 1, 1, 1],
               [1, 1, 1, 0, 1, 0, 0, 1],
               [1, 0, 1, 1, 0, 0, 0, 0]])
        >>> nz.astype(int)
        array([[0, 0, 0, 0, 1, 1, 0, 1],
               [0, 0, 0, 0, 1, 0, 1, 0],
               [0, 1, 0, 1, 0, 0, 0, 0],
               [0, 1, 0, 1, 0, 0, 1, 0]])

        >>> matrix = numpy.array([[0,1,0,0,1,0,1,0],
                               [0,0,0,0,1,1,0,1],
                               [1,1,1,0,0,1,0,0],
                               [1,1,0,1,0,0,0,0]], dtype=numpy.bool_)
        >>> nx, nz = hyperbolic_basis_for_pauli_group(matrix)
        array([[0, 1, 0, 0, 1, 0, 1, 0],
               [1, 0, 0, 1, 0, 1, 1, 1],
               [1, 1, 1, 0, 1, 0, 0, 1],
               [1, 0, 1, 1, 0, 0, 0, 0]])
        >>> nz.astype(int)
        array([[0, 0, 0, 0, 1, 1, 0, 1],
               [0, 0, 0, 0, 1, 0, 1, 0],
               [0, 1, 0, 1, 0, 0, 0, 0],
               [0, 1, 0, 1, 0, 0, 1, 0]])

        >>> x = numpy.array([[0, 1, 0, 0, 1, 0, 1, 0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0, 0, 0, 0, 1, 1, 0, 1]], dtype=numpy.bool_)
        >>> nx, nz = hyperbolic_basis_for_pauli_group(x=x, z=z)
        >>> nx.astype(int)
        array([[0, 1, 0, 0, 1, 0, 1, 0],
               [1, 1, 0, 0, 0, 1, 1, 1],
               [1, 1, 1, 0, 0, 0, 0, 0],
               [0, 1, 0, 1, 0, 0, 0, 0]])
        >>> nz.astype(int)
        array([[0, 0, 0, 0, 1, 1, 0, 1],
               [0, 0, 0, 0, 1, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 0, 1]])

        >>> x, z = hyperbolic_basis_for_pauli_group(n=5)
        >>> x.astype(int)
        array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 1, 0, 0, 0, 0, 0]])
        >>> z.astype(int)
        array([[0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]])

    See Also:
    _hyperbolic_basis_for_pauli_group
    """

    if (x is None) ^ (z is None):
        raise QiskitError("x and z must be both be None or both be array like")

    if x is not None:
        x = np.atleast_2d(np.array(x))
        z = np.atleast_2d(np.array(z))
        if not is_symplectic_matrix_form(x):
            raise QiskitError(f"{x} not a GF(2) symplectic matrix")
        if not is_symplectic_matrix_form(z):
            raise QiskitError(f"{z} not a symplectic matrix")

        if matrix is not None:
            matrix = np.atleast_2d(np.array(matrix))
            if not is_symplectic_matrix_form(matrix):
                raise QiskitError(f"{matrix} not a GF(2) symplectic matrix")
            if not is_center(matrix, np.vstack((matrix, x, z))):
                raise QiskitError(
                    "When providing an input hyperbolic set \
                    the input matrix must be the center of the full set of \
                        generators"
                )
        else:
            matrix = None
    else:
        if matrix is None:
            if n is None:
                raise QiskitError("If matrix, x and z are None then n must be provided")
            zero = np.zeros(shape=(n, n), dtype=np.bool_)
            x = mt.augment_mat(zero, "left")
            z = mt.augment_mat(zero, "right")
            return x, z
        matrix = np.atleast_2d(np.array(matrix))
        if not is_symplectic_matrix_form(matrix):
            raise QiskitError(f"{matrix} not a GF(2) symplectic matrix")
        x = []
        z = []

    return _hyperbolic_basis_for_pauli_group(matrix, x, z)


def _hyperbolic_basis_for_pauli_group(
    matrix: Optional[np.ndarray], x: Optional[np.ndarray], z: Optional[np.ndarray]
) -> Tuple[np.ndarray, np.ndarray]:
    """Builds a hyperbolic basis for the associated Pauli group

    Args:
        matrix (Optional): Input GF(2) symplectic matrix
        x, z (Optional): Input GF(2) hyperbolic pairs

    Returns:
        Hyperbolic basis for the associated Pauli group (new_x, new_z)

    Examples:
        >>> matrix = numpy.array([[0,1,0,0,1,0,1,0],
                               [0,0,0,0,1,1,0,1],
                               [1,1,1,0,0,1,0,0],
                               [1,1,0,1,0,0,0,0]], dtype=numpy.bool_)
        >>> center_, x, z = make_isotropic_hyperbolic_form(matrix)
        >>> nx, nz = _hyperbolic_basis_for_pauli_group(center_, x, z)
        >>> nx.astype(int)
        array([[0, 1, 0, 0, 1, 0, 1, 0],
               [1, 0, 0, 1, 0, 1, 1, 1],
               [1, 1, 1, 0, 1, 0, 0, 1],
               [1, 0, 1, 1, 0, 0, 0, 0]])
        >>> nz.astype(int)
        array([[0, 0, 0, 0, 1, 1, 0, 1],
               [0, 0, 0, 0, 1, 0, 1, 0],
               [0, 1, 0, 1, 0, 0, 0, 0],
               [0, 1, 0, 1, 0, 0, 1, 0]])

        >>> matrix = numpy.array([[0,1,0,0,1,0,1,0],
                               [0,0,0,0,1,1,0,1],
                               [1,1,1,0,0,1,0,0],
                               [1,1,0,1,0,0,0,0]], dtype=numpy.bool_)
        >>> nx, nz = _hyperbolic_basis_for_pauli_group(matrix)
        array([[0, 1, 0, 0, 1, 0, 1, 0],
               [1, 0, 0, 1, 0, 1, 1, 1],
               [1, 1, 1, 0, 1, 0, 0, 1],
               [1, 0, 1, 1, 0, 0, 0, 0]])
        >>> nz.astype(int)
        array([[0, 0, 0, 0, 1, 1, 0, 1],
               [0, 0, 0, 0, 1, 0, 1, 0],
               [0, 1, 0, 1, 0, 0, 0, 0],
               [0, 1, 0, 1, 0, 0, 1, 0]])

        >>> x = numpy.array([[0, 1, 0, 0, 1, 0, 1, 0]], dtype=numpy.bool_)
        >>> z = numpy.array([[0, 0, 0, 0, 1, 1, 0, 1]], dtype=numpy.bool_)
        >>> nx, nz = _hyperbolic_basis_for_pauli_group(x=x, z=z)
        >>> nx.astype(int)
        array([[0, 1, 0, 0, 1, 0, 1, 0],
               [1, 1, 0, 0, 0, 1, 1, 1],
               [1, 1, 1, 0, 0, 0, 0, 0],
               [0, 1, 0, 1, 0, 0, 0, 0]])
        >>> nz.astype(int)
        array([[0, 0, 0, 0, 1, 1, 0, 1],
               [0, 0, 0, 0, 1, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 0, 1]])

    See Also:
    hyperbolic_basis_for_pauli_group
    """

    if len(x) == 0:
        center_, x, z = _symplectic_gram_schmidt(matrix, [], [])
    else:
        center_ = matrix

    if center_ is not None:
        x, z = _make_hyperbolic(center_, x, z)

    basis = _basis_for_pauli_group(np.vstack((x, z)))
    added = x.shape[1] - 2 * x.shape[0]
    basis_com = _make_commute_hyper(basis[-added:], x, z)

    _, x_new, z_new = symplectic_gram_schmidt(basis_com, x, z)

    return x_new, z_new


def remove_hyper_elements_from_hyper_form(
    center_: Optional[np.ndarray], x: np.ndarray, z: np.ndarray, indices: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transfers those elements/vectors from x with an index in indices into center_
    and deletes the corresponding hyperbolic partner from z.

    Args:
        center_: center array
        x: Hyperbolic pairs, x being the source
        z: Hyperbolic pairs, x being the source
        indices: indices indicating which rows to transfer from source to sink

    Raises:
        QiskitError: x and z must be GF(2) symplectic matrices/vectors
        QiskitError: x and z must have the same shape
        QiskitError: x and z must have the same size in the second
                dimension as the center

    Returns:
        isotropic hyperbolic basis (center, x, z)

    Examples:
        >>> x, z = sysp.hyperbolic_basis_for_pauli_group(n=5)
        >>> center_, x_new, z_new = remove_hyper_elements_from_hyper_form(None,x,z,[0,1])
        >>> center_.astype(int)
        array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 0, 0, 0, 0, 0, 0, 0, 0]])
        >>> x_new.astype(int)
        array([[0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
              [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]])
        >>> z_new.astype(int)
        array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
              [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 0, 0, 1, 0, 0]])
    """
    x = np.atleast_2d(np.array(x))
    z = np.atleast_2d(np.array(z))
    indices = list(indices)

    if not (is_symplectic_matrix_form(x) and is_symplectic_matrix_form(z)):
        raise QiskitError("x and z must be GF(2) symplectic matrices/vectors")

    if not x.shape == z.shape:
        raise QiskitError(
            f"x (shape={x.shape})and z (shape={z.shape}) \
            must have the same shape"
        )

    if center_ is not None:
        center_ = np.atleast_2d(np.array(center_))
        if not is_symplectic_matrix_form(center_):
            raise QiskitError("Input center is not a GF(2) symplectiv matrix/vector")
        if not x.shape[1] == center_.shape[1]:
            raise QiskitError(
                "x and z must have the same size in the second \
                dimension as the center"
            )

    return _remove_hyper_elements_from_hyper_form(center_, x, z, indices)


def _remove_hyper_elements_from_hyper_form(
    center_: np.ndarray, x: np.ndarray, z: np.ndarray, indices: List[int]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transfers those elements/vectors from x with an index in indices into center_
    and delete the corresponding hyperbolic partner from z.

    Args:
        center_: center array
        x, z: Hyperbolic pairs, x being the source
        indices: indices indicating which rows to transfer from source to sink

    Returns:
        isotropic hyperbolic basis (center, x, z)

    Examples:
        >>> x, z = sysp.hyperbolic_basis_for_pauli_group(n=5)
        >>> center_, x_new, z_new = _remove_hyper_elements_from_hyper_form(None,x,z,[0,1])
        >>> center_.astype(int)
        array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 0, 0, 0, 0, 0, 0, 0, 0]])
        >>> x_new.astype(int)
        array([[0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
              [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]])
        >>> z_new.astype(int)
        array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
              [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 0, 0, 1, 0, 0]])
    """

    rm_size = len(indices)
    if rm_size == 0:
        return np.zeros(shape=(0, x.shape[1]), dtype=np.bool_), x, z
    size = x.shape[0] - rm_size
    shape = (size, x.shape[1])
    new_x = np.zeros(shape, dtype=np.bool_)
    new_z = np.zeros(shape, dtype=np.bool_)
    pos = 0
    part_center = deque()
    for i in reversed(range(x.shape[0])):
        if i in indices:
            part_center.appendleft(x[i].copy())
        else:
            new_x[pos] = x[i].copy()
            new_z[pos] = z[i].copy()
            pos += 1

    part_center = np.array(part_center)
    if center_ is None:
        new_center = part_center
    else:
        new_center = np.vstack((center_, part_center))

    if new_center.shape[0] > 1:
        new_center = new_center[np.where(new_center.any(axis=1))[0]]

    return new_center, new_x, new_z


def min_generating(
    matrix: Optional[np.ndarray] = None,
    x: Optional[np.ndarray] = None,
    z: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Returns a minimal generating/linearily independent set of rows.

    If only a matrix is provided, the method returns a submatrix with maximally
    independent set of rows. If a matrix and a set of hyperbolic pairs are provided then
    the method returns a submatrix such that the rows of the submatrix and the
    hyperbolic pairs are maximally linearily independent.

    Note: This method returns matrix as a 2d matrix (relative to numpy)

    Args:
        matrix (Optional[np.ndarray]): Input GF(2) symplectic matrix
        x (Optional[np.ndarray]): Input hyperbolic set - pair of GF(2) symplectic matrices
        z (Optional[np.ndarray]): Input hyperbolic set - pair of GF(2) symplectic matrices

    Raises:
        QiskitError: An input matrix is required
        QiskitError: Inputs x and z matrices do not have the same shape
        QiskitError: Input hyper pair matrices do not represent a set of Hyperbolic pairs
        QiskitError: Input matrix and x and z components must have the same number of columns

    Returns:
        matrix: minimal generating set or rows

    Examples:
        >>> mat = numpy.array([[1,0,1,0],[1,1,0,1],[0,1,1,1]])
        >>> res = min_generating(mat)
        >>> res.astype(int)
        array([[1,0,1,0],
               [1,1,0,1]])

        >>> mat = numpy.array([[1,0,1,0],[1,1,0,1],[0,1,1,1]])
        >>> x = numpy.array([[0,1,0,0]])
        >>> z = numpy.array([[0,0,0,1]])
        >>> res = symp.min_generating(mat,x,z)
        >>> res.astype(int)
        array([[1, 0, 1, 0],
               [1, 1, 0, 1]])

    See Also:
        _min_generating_matrix, _min_generating_matrix_xz
    """

    # all inputs are None

    if matrix is None and x is None and z is None:
        raise QiskitError("An input matrix is required")

    # Only matrix is provided
    if x is None and z is None:
        matrix = np.atleast_2d(matrix)
        return _min_generating_matrix(matrix)

    # Only x and z are provided
    if matrix is None:
        raise QiskitError("An input matrix is required")

    # All matrix, x and z are provided
    if x.shape != z.shape:
        raise QiskitError("Inputs x and z matrices do not have the same shape")
    matrix = np.atleast_2d(matrix)
    x = np.atleast_2d(x)
    z = np.atleast_2d(z)
    if matrix.shape[1] != x.shape[1]:
        raise QiskitError(
            "Input matrix and x and z components must have the same number of columns"
        )
    return _min_generating_matrix_xz(matrix, x, z)


def _min_generating_matrix(matrix: np.ndarray) -> np.ndarray:
    """Returns a matrix consisting of rows from the input matrix such that the resulting
    matrix and the input matrix have the same rank.

    Note: This method returns everything as 2d matrices

    Args:
        matrix: Input GF(2) symplectic matrix

    Examples:
        >>> mat = numpy.array([[1,0,1,0],[1,1,0,1],[0,1,1,1],[1,1,1,1],[0,0,0,0]])
        >>> _min_generating_matrix(mat)
        array([[1, 0, 1, 0],
               [1, 1, 0, 1],
               [0, 1, 1, 1]])

    Returns:
        matrix: submatrix of input matrix with full rank
    """
    heads, _, _, rank_ = mt.rref_complete(matrix.T)
    if rank_ == matrix.shape[0]:
        return matrix

    posns = np.flatnonzero(heads)
    ext_matrix = np.zeros(shape=(posns.shape[0], matrix.shape[1]), dtype=np.bool_)
    for k, index in enumerate(posns):
        ext_matrix[k] = matrix[index]
    return ext_matrix


def _min_generating_matrix_xz(matrix, x, z) -> np.ndarray:
    """Returns a submatrix of the input matrix such that combined rows of the submatrix and the
    hyperbolic pairs are linearily independent of GF(2)

    Args:
        matrix: GF(2) symplectic matrix
        x, y: hyperbolic basis pairs

    Returns:
        matrix: Submatrix of input matrix
    """
    # Priority is given to preserving the hyperbolic set (x,z) and so matrix is stack on the bottom.
    tmp_matrix = np.vstack((x, z, matrix))
    heads, _, _, rank_ = mt.rref_complete(tmp_matrix)
    if rank_ < tmp_matrix.shape[0]:
        # Since (x,z) has already been reduced any removal will appear in the matrix part
        posns = np.flatnonzero(heads[2 * x.shape[0] :])
        ext_matrix = np.zeros(shape=(posns.shape[0], matrix.shape[1]), dtype=np.bool_)
        for k, index in enumerate(posns):
            ext_matrix[k] = matrix[index]
        matrix = ext_matrix
    return matrix


def normalizer(
    matrix: Optional[np.ndarray] = None,
    x: Optional[np.ndarray] = None,
    z: Optional[np.ndarray] = None,
    min_gen: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns the normalizer of the group generated by the generators represented in the
    symplectic matrix(s):

    Args:
        matrix (optional): GF(2) symplectic matrix. Defaults to None.
        x (optional): Hyperbolic pairs. Defaults to None.
        z (optional): Hyperbolic pairs. Defaults to None.
        min_gen (optional): If True then the matrix will be reduced to be full
            rank (i.e. the rows will be a minimal generating set). Default is False

    Raises:
        QiskitError: All inputs should not be None
        QiskitError: matrix must be a GF(2) symplectic matrix
        QiskitError: x and z must be both be None or both be array like
        QiskitError: x and z must be GF(2) symplectic matrices/vectors
        QiskitError: x and z must have the same shape
        QiskitError: All inputs must have the same number of columns/length

    Returns:
        Isotropic hyperbolic form/basis of normalizer (center, x_new, y_new)

    Examples:
        >>> a = np.array([[1, 1, 1, 0, 1, 0, 0, 1]], dtype=np.bool_)
        >>> x = np.array([[0, 1, 0, 0, 1, 0, 1, 0]], dtype=np.bool_)
        >>> z = np.array([[0, 0, 0, 0, 1, 1, 0, 1]], dtype=np.bool_)
        >>> center_, x_new, z_new = sysp.normalizer(a, x, z)
        >>> center_.astype(int)
        array([[1, 1, 1, 0, 1, 0, 0, 1]])
        >>> x_new.astype(int)
        array([[1, 0, 1, 1, 0, 0, 0, 0],
               [1, 1, 1, 0, 0, 0, 1, 1],
               [0, 1, 0, 0, 1, 0, 1, 0]])
        >>> z_new.astype(int)
        array([[1, 1, 1, 0, 1, 0, 0, 0],
               [1, 1, 0, 0, 0, 1, 0, 0],
               [0, 0, 0, 0, 1, 1, 0, 1]])
    """

    if matrix is None and x is None and z is None:
        raise QiskitError("All inputs should not be None")
    if matrix is not None:
        matrix = np.atleast_2d(np.array(matrix))
        if not is_symplectic_matrix_form(matrix):
            raise QiskitError(f"{matrix} must be a GF(2) symplectic matrix")

    if (x is None) ^ (z is None):
        raise QiskitError("x and z must be both be None or both be array like")

    if matrix is not None and x is None:
        if min_gen:
            matrix = min_generating(matrix)
        if is_stabilizer_group(matrix):
            return _normalizer_abelian_group(matrix)
        else:
            return _normailzer_group(matrix)

    x = np.atleast_2d(np.array(x))
    z = np.atleast_2d(np.array(z))
    if not (is_symplectic_matrix_form(x) and is_symplectic_matrix_form(z)):
        raise QiskitError("x and z must be GF(2) symplectic matrices/vectors")

    zero_mat = False
    if matrix is None:
        matrix = np.zeros(shape=(0, x.shape[1]), dtype=np.bool_)
        zero_mat = True

    if not x.shape == z.shape:
        raise QiskitError("x and z must have the same shape")

    if not matrix.shape[1] == x.shape[1]:
        raise QiskitError("All inputs must have the same number of columns/length")

    if not is_center(matrix, np.vstack((matrix, x, z))):
        raise QiskitError(
            "When a matrix and hyperbolic parts are "
            "provided the matrix must represent the center of the union"
        )
    if not zero_mat:
        matrix = min_generating(matrix, x, z)
    return _normalizer_group_preserve(matrix, x, z)


def _normalizer_abelian_group(
    matrix: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns the normalizer of the abelian group generated by the generators represented in the
    symplectic matrix(s):

    Args:
        matrix (optional): GF(2) symplectic matrix. Defaults to None.
        x, z  (optional): Hyperbolic pairs. Defaults to None.

    Returns:
        Isotropic hyperbolic form/basis of normalizer (center, x_new, y_new)

    TODO: Add examples
    """
    dist_center = mt.rank(matrix)

    matrix_ext = _basis_for_pauli_group(matrix)
    center_, x, z = _symplectic_gram_schmidt(matrix_ext, [], [])
    center_, x, z = _remove_hyper_elements_from_hyper_form(center_, x, z, list(range(dist_center)))

    if center_.shape[0] > 1:
        center_ = center_[np.where(center_.any(axis=1))[0]]

    return center_, x, z


def _normalizer_group_preserve(
    center_: np.ndarray, x: np.ndarray, z: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns the normalizer of the group generated by the generators represented in the
    symplectic matrix(s) and trys to preserve elements:

    This method assumes that the matrix has full rank. That is that the
    set of generators is independent.

    Args:
        matrix (optional): GF(2) symplectic matrix. Defaults to None.
        x, z  (optional): Hyperbolic pairs. Defaults to None.

    Returns:
        Isotropic hyperbolic form/basis of normalizer (center, x_new, y_new)

    TODO: Add examples
    """
    center_size = center_.shape[0]
    if center_size == 1:
        if not center_[0].any():
            center_size = 0

    gauge_degree = x.shape[0]

    if center_size > 0:
        x, z = _make_hyperbolic(center_, x, z)

    matrix_ext = _basis_for_pauli_group(np.vstack((x, z)))
    matrix_ext = make_commute_hyper(matrix_ext, x, z, range(x.shape[0] << 1, matrix_ext.shape[0]))
    # matrix_ext = _make_elements_commute_with_hyper_pairs(
    #    matrix_ext,
    #    range(x.shape[0] << 1, matrix_ext.shape[0]),
    #    x,
    #    range(x.shape[0]),
    #    z,
    #    range(z.shape[0]),
    # )
    matrix = matrix_ext[x.shape[0] << 1 :]
    lx = [item.copy() for item in matrix_ext[: x.shape[0]]]
    lz = [item.copy() for item in matrix_ext[x.shape[0] : x.shape[0] << 1]]
    center_, x, z = _symplectic_gram_schmidt(matrix, lx, lz)
    indices = list(range(gauge_degree, gauge_degree + center_size))
    return _remove_hyper_elements_from_hyper_form(center_, x, z, indices)


def _normailzer_group(matrix):
    """Returns the normalizer of the group generated by the generators represented in the
    symplectic matrix(s):

    This method assumes that the matrix has full rank. That is that the
    set of generators is independent.

    Args:
        matrix: GF(2) Sympletic matrix

    Returns:
        Isotropic hyperbolic form/basis of normalizer (center, x_new, y_new)

    TODO: Add examples
    """
    center_, hyper1, hyper2 = _symplectic_gram_schmidt(matrix, [], [])

    return _normalizer_group_preserve(center_, hyper1, hyper2)
