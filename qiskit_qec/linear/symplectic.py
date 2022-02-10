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

from matplotlib.contour import QuadContourSet
import numpy as np

from collections import deque

from numpy.typing import ArrayLike
from typing import List, Any
from typing import Union

from qiskit import QiskitError
from sympy import Q
from qiskit_qec.linear import matrix as mt

def all_commute(check_matrix):
    """Determine is the operators of a check_matrix generate an abelian group"""
    test_mat = symplectic_product(check_matrix, check_matrix)
    return not test_mat.any()

# ---------------------------------------------------------------

def symplectic_product(
        mat1: ArrayLike, 
        mat2: ArrayLike):
    """General Symplectic Product of two binary matrices

    Args:
        mat1 (ArrayLike): binary matrix
        mat2 (ArrayLike): binary matrix

    Returns:
        int (0,1): Symplectic product of mat1 and mat2
    """
    mat1_np_array = np.array(mat1, dtype=np.int8)
    mat2_np_array = np.array(mat2, dtype=np.int8)

    assert np.array_equal(mat1_np_array, mat1_np_array % 2), \
        QiskitError(f"Input matrices must be binary entries: {mat1}")
    assert np.array_equal(mat2_np_array, mat2_np_array % 2), \
        QiskitError(f"Input matrices must be binary entries: {mat2}")

    assert mat1_np_array.ndim == mat2_np_array.ndim, \
        QiskitError(f"Input matrices must be of the \
        same dimensions: {mat1_np_array.ndim} is \
        not equal to {mat2_np_array.ndim}")

    if mat1_np_array.ndim == 1:
        assert hasattr(mat1, "__getitem__") \
            and hasattr(mat2, "__getitem__"), \
            QiskitError(f"Input matrices must be \
            scriptable objects")
        assert mat1_np_array.shape[0] == mat2_np_array.shape[0], \
            QiskitError(f"Input vectors must have the same \
            dimensions: {mat1_np_array.shape[0]} not equal \
            to {mat2_np_array.shape[0]}")
        
        assert mat1_np_array.shape[0]%2 == 0, \
            QiskitError(f"Input vectors must be of an even length: {mat1_np_array.shape[0]}")
        assert mat2_np_array.shape[0]%2 == 0, \
            QiskitError(f"Input vectors must be of an even length: {mat2_np_array.shape[0]}")

        return _symplectic_product_vv(mat1, mat2, mat1_np_array.shape[0]>>1)

    elif mat1_np_array.ndim == 2:
        assert mat1_np_array.shape[1]%2 == 0, \
            QiskitError(f"Input matrices must have an even number of columns: {mat1_np_array.shape[1]}")
        assert mat2_np_array.shape[1]%2 == 0, \
            QiskitError(f"Input vectors must have an even number of columns: {mat2_np_array.shape[1]}")
        
        return _symplectic_product_numpy(mat1_np_array, mat2_np_array)

    else:
        raise QiskitError(f"Input matrices must be 1 or 2 dimensional:\
            {mat1_np_array.ndim}, {mat2_np_array.ndim}")

def _symplectic_product_vv(
        vec1: ArrayLike, 
        vec2: ArrayLike, 
        n: int):
    """Find the sympletic product or two binary vectors of 
    length 2n: vec1 . Lambda . vec2^T where

    lambda = [0 I]
             [I 0]

    Args:
        vec1 (ArrayLike): [description]
        vec2 (ArrayLike): [description]
        n (int): [description]

    Returns:
        int (0,1): Symplectic product of vec1 and vec2

    Optimization:
    If changing the method please make sure that the new method
    is faster than the method currently used. Note that this is 
    faster if the vectors are numpy arrays with dtype=int8
    """
    r=0
    for i in range(n):
        r += vec1[i]*vec2[n+i] + vec1[n+i]*vec2[i]
    return r%2
        
def _symplectic_product_numpy(
        mat1: np.ndarray, 
        mat2: np.ndarray):
    """Find the sympletic product or two binary matrices. Matrices 
    must have numpy.int8 values. numpy.bool_ and other bool entries must be 
    converted to numpy.int8 before using - as numpy.hstack((m2, m1)).dot(mat2.transpose())
    will proced a numpy.bool_ matrix which will distroy the needed p[arity
    measurement

    lambda = [0 I]
             [I 0]

    Args:
        mat1 (numpy.ndarray): Binary symplectic matrix
        mat2 (numpy.ndarray): Binary symplectic matrix

    Returns:
        [numpy.ndarray (dtype=numpy.int8)]: Symplectic products of the rows of of mat1 and mat2
    """
    m1, m2 = np.hsplit(mat1, 2)
    return np.hstack((m2, m1)).dot(mat2.transpose()) % 2

# ---------------------------------------------------------------

def make_element_commute_with_hyper_pair(
        vector: np.ndarray, 
        hyper1: np.ndarray, 
        hyper2: np.ndarray):
    """Modify the operator given by vector so that it commutes
    with the hyperbolic pair (hyper1,hyper2) such that this new vector
    is in <vector, hyper1, hyper2>

    Args:
        vector (numpy.ndarray dtype=bool): Symplectic vector encoding a Pauli
        hyper1 (numpy.ndarray dtype=bool): Symplectic vector encoding one of 
            the operators in the hyperbolic elements 
        hyper2 (numpy.ndarray dtype=bool): Symplectic vector encoding other one of 
            the operators in the hyperbolic elements

    Returns:
        numpy.ndarray dtype=bool : Symplectic vector encoding a modified vector that commutes
            with the hyperbolic pair and is in <vector, hyper1, hyper2>
    """
    assert vector.ndim == 1 and hyper1.ndim == 1 and hyper2.ndim == 1, \
        QiskitError(f"All inputs must be one dimensional: {vector.ndim},{hyper1.ndim}, {hyper2.ndim}")
    assert not (vector.shape[0]%2 or hyper1.shape[0]%2 or hyper2.shape[0]%2), \
        QiskitError(f"All vectors must have an even length: \
        {vector.shape[0]},{hyper1.shape[0]},{hyper2.shape[0]}")
    return _make_element_commute_with_hyper_pair(vector, hyper1, hyper2)


def _make_element_commute_with_hyper_pair(
        vector: np.ndarray, 
        hyper1: np.ndarray, 
        hyper2: np.ndarray):
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
    num_qubits = hyper1.shape[0]>>1
    if _symplectic_product_vv(new_vector, hyper1, num_qubits):
        new_vector = new_vector ^ hyper2
    if _symplectic_product_vv(new_vector, hyper2, num_qubits):
        new_vector = new_vector ^ hyper1
    return new_vector

# ---------------------------------------------------------------

def make_element_commute_with_hyper_pairs(
        vector: np.ndarray, 
        hyper1: np.ndarray, 
        hyper2: np.ndarray, 
        range1: Union[List,np.ndarray], 
        range2: Union[List,np.ndarray]):
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
    assert vector.ndim == 1 and hyper1.ndim > 1 and hyper2.ndim > 1, \
        QiskitError(f"All inputs must be one dimensional: {vector.ndim},{hyper1.ndim}, {hyper2.ndim}")
    assert vector.shape[0] == hyper1.shape[1] == hyper2.shape[1], \
        QiskitError("Inputs matrices/vectors must have the same number of columns/length")
    assert not (vector.shape[0]%2 or hyper1.shape[1]%2 or hyper2.shape[1]%2), \
        QiskitError(f"All vectors must have an even length: \
        {vector.shape[0]},{hyper1.shape[0]},{hyper2.shape[0]}")
    try:    
        range1 = list(range1)
    except TypeError:
        QiskitError("Input range1 not iterable")
    try:
        range2 = list(range2)
    except TypeError:
        QiskitError("Input range2 not iterable")
  
    assert set(range1).issubset(range(hyper1.shape[0])) and set(range2).issubset(range(hyper2.shape[0])), \
        QiskitError(f"Input ranges not valid per hyper1 and hyper2 inputs ranges")

    return _make_element_commute_with_hyper_pairs(vector, hyper1, hyper2, range1, range2)

def _make_element_commute_with_hyper_pairs(
        vector: np.ndarray, 
        hyper1: np.ndarray, 
        hyper2: np.ndarray, 
        range1: List, 
        range2: List):
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
        new_vector = _make_element_commute_with_hyper_pair(new_vector, hyper1[index_1], hyper2[index_2])

    return new_vector

# ---------------------------------------------------------------

def make_elements_commute_with_hyper_pair(
        matrix: ArrayLike, 
        mrange: Union[range, List[int], np.ndarray], 
        hyper1: ArrayLike, 
        hyper2: ArrayLike):
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

    assert matrix.ndim == 2 and hyper1.ndim == 1 and hyper2.ndim == 1, \
        QiskitError(f"input matrix must be 2 dimensional and hyper1, \
        hyper2 must be one dimernsional: {matrix.ndim}(2), {hyper1.ndim}(1), {hyper2.ndim}(1)")
    assert not (matrix.shape[1]%2 or hyper1.shape[0]%2 or hyper2.shape[0]%2), \
        QiskitError(f"Inputs must have an even number of columns/length:\
        {matrix.shape[1]},{hyper1.shape[0]},{hyper2.shape[0]}")
    assert matrix.shape[1]==hyper1.shape[0]==hyper2.shape[0], \
        QiskitError(f"Input matices/vectors must have the same number of columns/length")

    assert set(mrange).issubset(range(matrix.shape[1])), \
        QiskitError(f"Input range not a valid range for input matrix: {mrange}")

    return _make_elements_commute_with_hyper_pair(matrix, mrange, hyper1, hyper2)


def _make_elements_commute_with_hyper_pair(
        matrix: np.ndarray, 
        mrange: List, 
        hyper1: np.ndarray, 
        hyper2: np.ndarray):
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

    # ---------------------------------------------------------------

def find_noncommutative_partner(matrix: np.ndarray, vector: np.ndarray):
    """Find a operator/vector in the input matrix that anticommutes
    with the input vector (writh respect to the sympletic product)

    Args:
        matrix (numpy.ndarray): Check matrix
        vector (numpy.ndarray): Vector/operator to find a noncommuting pair

    Returns:
        numpy.ndarray: a operator/vector in the checlk matrix that anticommutes
        with the input vectir re;lative to the symplectic product. If no such
        vector is found then None is returned.
    """
    assert matrix.ndim == 2, QiskitError(f"{matrix} must be a 2 dimensional array")
    assert vector.ndim == 1, QiskitError(f"{vector} must be a 1 dimensional array")
    assert matrix.shape[1] == vector.shape[0], \
        QiskitError(f"Input matrix and vector must have the same number \
            of columns/length {matrix.shape[1]}!={vector.shape[0]}")

    return _find_noncommutative_partner(matrix, vector)

def _find_noncommutative_partner(matrix: np.ndarray, vector: np.ndarray):
    """Find a operator/vector in the input matrix that anticommutes
    with the input vector (writh respect to the sympletic product)

    Args:
        matrix (numpy.ndarray): Check matrix
        vector (numpy.ndarray): Vector/operator to find a noncommuting pair

    Returns:
        numpy.ndarray: a operator/vector in the checlk matrix that anticommutes
        with the input vectir re;lative to the symplectic product. If no such
        vector is found then None is returned.

    Note: operator returned is a copy. Using index to maintain a reference to
    original input matrix
    """
    n = matrix.shape[1]>>1
    for index, item in enumerate(matrix):
        if _symplectic_product_vv(item, vector, n) == 1: 
            return (item.copy(), index)
    return None

 # ---------------------------------------------------------------

def find_hyperbolic_partner(check_matrix, index: int):
    check_matrix = np.atleast_2d(np.array(check_matrix,  dtype=bool))

    # check_matrix -> all associated operators must commute
    assert all_commute(check_matrix) == True, \
        QiskitError("Check matrix must represent a set of commuting operators")

    assert check_matrix.shape[1] % 2 == 0, \
        QiskitError(f"Input check matrix must have an even number of columns: {check_matrix.shape[1]}")

    assert index in range(check_matrix.shape[1]>>1), \
        QiskitError(f"index out or range: {index}>={check_matrix.shape[1]>>1}")

    return _find_hyperbolic_partner(check_matrix, index)


def _find_hyperbolic_partner(check_matrix, index: int):
    """Find Hyperbolic partner to a given generator of a set of commuting generators.

    Given an set of generators S=<g_1,g_2,...,g_k> in symbolic matrix form,
    find a hyperbolic element g in P_n to g_i
    So g commutes with every generator g_j for i different from i but anticommutes with g_i
    Reference Prop 10.4. Nielsen and Chuang
    Assumption is that S contains order two elements all of which commute

    Args:
        generators (numpy.ndarry): Symplectic representation of generators (number of rows = number of generators)
        index (int): index of generator to find partner for

    Returns:
        numpy.ndarray: Symplectic representation of hyperbolic partner for g_index
    """

    nrows = check_matrix.shape[0]
    ncols = check_matrix.shape[1]

    _lambda = mt.create_lambda_matrix(ncols)
    slambda = np.matmul(check_matrix, _lambda)

    heads, rref_mat, transform_mat, rank = mt._rref_complete(slambda)

    e_index = np.zeros(nrows , dtype=bool)
    e_index[index] = True

    trans_e_index = np.matmul(transform_mat, e_index)

    pivot = 0
    result = np.zeros(ncols, dtype=bool)
    for i in range(ncols):
        if heads[i] == 1:
            result[i] = trans_e_index[pivot]
            pivot =+ 1
    
    return result

# ---------------------------------------------------------------

def symplectic_gram_schmidt(
        in_matrix: np.ndarray, 
        hyper1: np.ndarray=None, 
        hyper2: np.ndarray=None):
    """Apply the symplectic GramSchmidt process to the input symplectic matrix. Resulting
    hyperbolic pairs are added to hyper1 and hyper2. Resulting center operators/vecrtors
    are added to center.

    Args:
        in_matrix (np.ndarray): Symplectic matrix
        hyper1 (np.ndarray, optional): Symplectic matrix representing one element of
            hyperbolic pairs. Defaults to None.
        hyper2 (np.ndarray, optional): Symplectic matrix representing other element of\
            hyperbolic pairs. Defaults to None.

    Returns:
        center (np.ndarray), hyper1 (np.ndarray), hyper2 (np.ndarray): Center array, 
        hyperbolic pairs with one element in hyper1 and other element in hyper2
    """
    try:
        hyper1 = [item for item in hyper1]
    except TypeError:
        hyper1 = []

    try:
        hyper2 = [item for item in hyper2]
    except TypeError:
        hyper2 = []

    return _symplectic_gram_schmidt(in_matrix, hyper1, hyper2)

def _symplectic_gram_schmidt(
        in_matrix: np.ndarray, 
        hyper1: List[np.ndarray], 
        hyper2: List[np.ndarray]):
    """Apply the symplectic GramSchmidt process to the input symplectic matrix. Resulting
    hyperbolic pairs are added to hyper1 and hyper2. Resulting center operators/vectors
    are added to center.

    Note: Center elements are added to hyper1

    Args:
        matrix (nd.ndarray bool): Symplectic matrix
        hyper1 (List[nd.ndarray bool]): Symplectic matrix representing one element of hyperbolic pairs
        hyper2 (List[nd.ndarray bool]): Symplectic matrix representing other element of hyperbolic pairs

    Returns:
        center (np.ndarray), hyper1, hyper2: 
    """

    matrix = in_matrix.copy()
    matrix_view = matrix
    center_ = []

    while matrix_view.shape[0] > 0:
        elem = matrix_view[0]  
        # Remove elem from matrix_view
        matrix_view = matrix_view[1:]
        try:
            elem_p, index = _find_noncommutative_partner(matrix_view, elem)
            hyper1.append(elem)
            hyper2.append(elem_p)

            # Revove elem_p from matrix_view
            temp_view = matrix_view[:-1]
            temp_view[index:] = matrix_view[index+1:]
            matrix_view = temp_view

            matrix_view = _make_elements_commute_with_hyper_pair(
                            matrix_view, 
                            range(matrix_view.shape[0]), 
                            elem, 
                            elem_p)
        except TypeError:
            center_.append(elem)

    hyper1 = np.asarray(hyper1)
    hyper2 = np.asarray(hyper2)
    if len(center_) == 0:
        center_ = np.zeros(shape=(0,hyper1.shape[1]), dtype=np.bool_)
    else:
        center_ = np.asarray(center_)
            
    return center_, hyper1, hyper2

# ---------------------------------------------------------------

def is_symplectic_matrix_form(
        matrix: np.ndarray, 
        dtype: Union[None, bool, np.bool_, int, np.integer]=None) -> bool:
    """Check if the given numpy matrix is in the form of a symplectic matrix: two dimensional, 
    even number of columns, 0/1 or boolean entries. Optional argument dtype can be given to check if
    entries are a specific dtype

    Args:
        matrix (np.ndarray): Input matrix to be checked
        dtype (Union[bool, np.bool_, int, np.integer]): Optional. Check if given matrix 
        is of type dtypeDefault: None

    Returns:
        bool: True if matrix is of symplectic form
    """
    if matrix.ndim != 2:
        return False
    if matrix.shape[1]%2:
        return False
    if not np.array_equal(matrix, matrix%2):
        return False
    if dtype is None:
        return True
    if not isinstance(matrix[0][0], dtype):
        return False
    return True

# ---------------------------------------------------------------

def is_center(center_: ArrayLike, matrix: ArrayLike):
    """Does <center_> = Z(<matrix>)?

    Args:
        center_ (ArrayLike): Generators of center to be checked
        matrix (ArrayLike): Generators of full group

    Returns:
        boolen: True if <center_> = Z(<matrix>), False otherwise

    """
    assert is_symplectic_matrix_form(center_) and is_symplectic_matrix_form(matrix), \
        QiskitError(f"Not all inputs are not in symplectic form")
    cal_center = center(matrix)

    return is_same_span(center_, cal_center)

# ---------------------------------------------------------------  

def is_same_span(matrix1: ArrayLike, matrix2: ArrayLike):
    """Does <matrix1> = <matrix2>

    Args:
        matrix1 (ArrayLike): First set of generators
        matrix2 (ArrayLike): Second set of generators

    Returns:
        boolean: True if <matrix1> = <matrix2>, False otherwise
    """
    _, rref_matrix1, _, rank_matrix1 = mt.rref_complete(matrix1)
    _, rref_matrix2, _, rank_matrix2 = mt.rref_complete(matrix2)

    if rank_matrix1 != rank_matrix2:
        return False

    rref_matrix1 = rref_matrix1[:rank_matrix1]
    rref_matrix2 = rref_matrix2[:rank_matrix2]
    return np.array_equal(rref_matrix1, rref_matrix2)

def is_hyper_form(hyper1, hyper2):
    """Do hyper1 and hyper2 contains hyperbolic/symnplectic pairs?

    Args:
        hyper1 ([type]): [description]
        hyper2 ([type]): [description]

    Returns:
        [type]: [description]
    """
    matrix = np.vstack((hyper1, hyper2))
    test = symplectic_product(matrix, matrix)
    return np.array_equal(test, mt._create_lambda_matrix(matrix.shape[0]))

def is_stabilizer_group(matrix: ArrayLike):
    matrix = np.array(matrix)
    assert is_symplectic_matrix_form(matrix), QiskitError(f"matrix not symplectic matrix")
    return all_commute(matrix)


# ---------------------------------------------------------------

def center(matrix: ArrayLike, preserve: bool=False):
    """Find the center of the group with generators given by the symplectic matrix

    Args:
        matrix (ArrayLike): Symplectic matrix
        preserve (bool): If True then an attempt will be made to preserve generators

    Returns:
        numpy.ndarray : Generators for the center, represented as a symplectic matrix,
        of the group with generators given by the input  symplectic matrix


    GetCenter := function(G)
    return GramSchmidt(G).center;
    end;

    """
    matrix = np.atleast_2d(np.array(matrix))
    assert is_symplectic_matrix_form(matrix), QiskitError("Input matrix is not a symplectic matrix")
    if preserve:
        return _center_preserve(matrix)
    else:
        return _center(matrix)

def _center(matrix: np.ndarray):
    """Find the center of the group with generators given by the symplectic matrix

    Args:
        matrix (numpy.ndarray): Symplectic matrix

    Returns:
        numpy.ndarray : Generators for the center, represented as a symplectic matrix,
        of the group with generators given by the input  symplectic matrix

    Method:
        Uses the symplectic Gram Schmidt process to find the center.
    """
    center, hyper1, hyper2 = _symplectic_gram_schmidt(matrix, [], [])
    return center

def _center_preserve(matrix: np.ndarray ):
    """Find the center but try to maintain any generator that is in the center as a 
    generator of the center

    Args:
        matrix (np.ndarray): Symplectic matrix of generators

    Returns:
        nd.ndarray: center
    """
    # First move any generator that is in the center to the front of the list
    rematrix = deque()
    num = matrix.shape[1]>>1
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

    return _make_elements_commute_with_hyper_pairs(in_matrix, mrange, hyper1, h1range, hyper2, h2range)
    

def _make_elements_commute_with_hyper_pairs(in_matrix, mrange, hyper1, h1range, hyper2, h2range):
    matrix = in_matrix.copy()

    for i,j in zip(h1range, h2range):
        matrix = _make_elements_commute_with_hyper_pair(matrix, mrange, hyper1[i], hyper2[j])

    return matrix

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
        in_center = np.zeros((0,in_hyper1.shape[1]))

    assert is_center(center, np.vstack((hyper1, hyper2, center))), QiskitError(f"Input center is not center")
    assert is_hyper_form(hyper1, hyper2), QiskitError(f"Inpur hyper1,hyper2 are not in hyperbolic pairs")

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
        hop = _find_hyperbolic_partner(center_[:center_size], center_size-1)
        hop = _make_element_commute_with_hyper_pairs(
            hop, 
            hyper1,  
            hyper2,
            range(hyper_size),
            range(hyper_size))
    
        center_size -= 1
        hyper1[hyper_size] = center_[center_size]
        hyper2[hyper_size] = hop
        hyper_size += 1

    return hyper1, hyper2

def isotropic_hyperbolic_form(matrix):
    return _symplectic_gram_schmidt(matrix, [], [])

def isotropic_hyperbolic_basis(matrix: Union[None,ArrayLike], hyper1: Union[None, ArrayLike], hyper2: Union[None, ArrayLike]):  

    if (hyper1 is None) ^ (hyper2 is None):
        raise QiskitError(f"hyper1 and hyper2 must be both be None or both be array like")

    if hyper1 is not None:
        hyper1 = np.array(hyper1)
        hyper2 = np.array(hyper2)
        assert is_symplectic_matrix_form(hyper1), QiskitError(f"{hyper1} not a symplectic matrix")
        assert is_symplectic_matrix_form(hyper2), QiskitError(f"{hyper2} not a symplectic matrix")

        if matrix is not None:
            matrix = np.array(matrix)
            assert is_symplectic_matrix_form(matrix), QiskitError(f"{matrix} not a symplectic matrix")
            assert is_center(matrix, np.vstack(matrix, hyper1, hyper2))
        else:
            matrix = None
    else:
        assert matrix is not None, QiskitError("At one of matrix, hyper1 and hyper2 cannot all be None")
        matrix = np.array(matrix)
        assert is_symplectic_matrix_form(matrix), QiskitError(f"{matrix} not a symplectic matrix")
        hyper1 = []
        hyper2 = []

    return _isotropic_hyperbolic_basis(matrix, hyper1, hyper2)

def _isotropic_hyperbolic_basis(matrix: Union[None,np.ndarray], hyper1: List[np.ndarray], hyper2: List[np.ndarray]):
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
    
    added = hyper1.shape[1]-2*hyper1.shape[0]

    basis_com = _make_elements_commute_with_hyper_pairs(
        basis[0:added], 
        range(2*hyper1.shape[0],basis.shape[0]),
        hyper1,
        range(hyper1.shape[0]),
        hyper2,
        range(hyper1.shape[0]))
    
    _, hyper1_ans, hyper2_ans = symplectic_gram_schmidt(basis_com, hyper1, hyper2)

    return hyper1_ans, hyper2_ans

def remove_hyper_elements_from_hyper_form(
        hyper1: ArrayLike, 
        hyper2: ArrayLike,
        center_: ArrayLike,
        indices: ArrayLike):
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
    assert hyper1.shape == hyper2.shape, \
        QiskitError(f"hyper1 (shape={hyper1.shape})and hyper2 (shape={hyper2.shape}) \
            must have the same shape")
    if center_ is not None:
        center_ = np.atleast_2d(np.array(center_))
        assert hyper1.ndim == center_.ndim, \
            QiskitError(f"center must have the same number of dimensions as hyper1 and hyper2")
        assert hyper1.shape[1] == center_.shape[1], \
            QiskitError(f"hyper1 and hyper2 must have the same size in the second dimension as the center_")
        
    return _remove_hyper_elements_from_hyper_form(hyper1, hyper2, center_, indices)

def _remove_hyper_elements_from_hyper_form(
        hyper1: np.ndarray, 
        hyper2: np.ndarray,
        center_: np.ndarray,
        indices: List[int]):
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
    shape = (size,hyper1.shape[1])
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
    matrix: Union[None, ArrayLike]=None, 
    hyper1: Union[None, ArrayLike]=None,
    hyper2: Union[None, ArrayLike]=None):
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
    
    assert not (matrix is None and hyper1 is None and hyper2 is None), \
        QiskitError("All inputs should not be None")
    if matrix is not None:    
        matrix = np.array(matrix)
        assert is_symplectic_matrix_form(matrix), QiskitError(f"{matrix} must be a symplectic matrix")

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
    assert matrix.shape[1] == hyper1.shape[1], \
        QiskitError(f"All inputs must have the same number of columns/length")
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
                                hyper1, 
                                hyper2, 
                                center_, 
                                list(range(dist_center)))
    return center_, hyper1, hyper2

def _normalizer_gauge_group_preserve(
        center_: np.ndarray, 
        hyper1: np.ndarray, 
        hyper2: np.ndarray):
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
                    range(hyper1.shape[0]<<1, matrix_ext.shape[0]),
                    hyper1,
                    range(hyper1.shape[0]),
                    hyper2,
                    range(hyper2.shape[0]))
    matrix = matrix_ext[hyper1.shape[0]<<1:]
    lhyper1 = [item.copy() for item in matrix_ext[:hyper1.shape[0]]]
    lhyper2 = [item.copy() for item in matrix_ext[hyper1.shape[0]: hyper1.shape[0]<<1]]
    center_, hyper1, hyper2 = _symplectic_gram_schmidt(matrix, lhyper1, lhyper2)
    indices = list(range(gauge_degree,gauge_degree+center_size))
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

