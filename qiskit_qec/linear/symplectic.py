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

def all_commute(check_matrix):
    """Determine is the operators of a check_matrix generate an abelian group"""
    test_mat = symplectic_product(check_matrix, check_matrix)
    return not test_mat.any()

def symplectic_product(mat1, mat2):
    """General Symplectic Product of two binary matrices

    Args:
        mat1 (any): binary matrix
        mat2 (any): binary matrix

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

def _symplectic_product_vv(vec1, vec2, n):
    """Find the sympletic product or two binary vectors of 
    length 2n: vec1 . Lambda . vec2^T where

    lambda = [0 I]
             [I 0]

    Args:
        vec1 (array): [description]
        vec2 (array): [description]
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
        
def _symplectic_product_numpy(mat1, mat2):
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

def make_element_commute_with_hyper_pair(vector: np.ndarray, hyper1: np.ndarray, hyper2: np.ndarray):
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


def _make_element_commute_with_hyper_pair(vector, hyper1, hyper2):
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
    new_vector = vector.copy()
    num_qubits = hyper1.shape[0]>>1
    if _symplectic_product_vv(new_vector, hyper1, num_qubits):
        new_vector = new_vector ^ hyper2
    if _symplectic_product_vv(new_vector, hyper2, num_qubits):
        new_vector = new_vector ^ hyper1
    return new_vector


def _make_element_commute_with_hyper_pairs(
        vector: np.ndarray, 
        hyper1: np.ndarray, 
        hyper2: np.ndarray, 
        range1: np.ndarray, 
        range2: np.ndarray):
    """[summary]

    # Make the operator elem commute with the hyperbolic pairs Xe(xrange) and Ze(zrange)
    MakeElementCommuteGeneral := function(elem, Xe, xrange, Ze, zrange)
    local  i, elem_mod;

    elem_mod := StructuralCopy(elem);

    for i in [1..Size(xrange)] do
        elem_mod := MakeElementCommute(elem_mod, Xe[xrange[i]], Ze[zrange[i]]);
    od;

    return elem_mod;
    end;

    Args:
        vector ([type]): [description]
        hyper1 ([type]): [description]
        hyper2 ([type]): [description]
        range1 ([type]): [description]
        range2 ([type]): [description]

    Returns:
        [type]: [description]
    """
    assert vector.ndim == 1 and hyper1.ndim > 1 and hyper2.ndim > 1, \
        QiskitError(f"All inputs must be one dimensional: {vector.ndim},{hyper1.ndim}, {hyper2.ndim}")
    assert vector.shape[0] == hyper1.shape[1] == hyper2.shape[1], \
        QiskitError("Inputs matrices/vectors must have the same number of columns/length")
    assert not (vector.shape[0]%2 or hyper1.shape[0]%2 or hyper2.shape[0]%2), \
        QiskitError(f"All vectors must have an even length: \
        {vector.shape[0]},{hyper1.shape[0]},{hyper2.shape[0]}")
    assert set(range1).issubset(range(hyper1.shape[0])) and set(range2).issubset(range(hyper2.shape[0])) 



def _make_element_commute_with_hyper_pairs(vector, hyper1, hyper2, range1, range2):
    """[summary]

    # Make the operator elem commute with the hyperbolic pairs Xe(xrange) and Ze(zrange)
    MakeElementCommuteGeneral := function(elem, Xe, xrange, Ze, zrange)
    local  i, elem_mod;

    elem_mod := StructuralCopy(elem);

    for i in [1..Size(xrange)] do
        elem_mod := MakeElementCommute(elem_mod, Xe[xrange[i]], Ze[zrange[i]]);
    od;

    return elem_mod;
    end;

    Args:
        vector ([type]): [description]
        hyper1 ([type]): [description]
        hyper2 ([type]): [description]
        range1 ([type]): [description]
        range2 ([type]): [description]

    Returns:
        [type]: [description]
    """
    new_vector = vector.copy()

    for index_1, index_2 in zip(range1, range2):
        new_vector = _make_element_commute_with_hyper_pair(new_vector, hyper1, hyper2)

    return new_vector
    

def _make_elements_commute_with_hyper_pair(matrix, range, hyper1, hyper2):
    """[summary]

    # Make the elements in S(range) commute with the hyperbolic pair Xe, Ze
    MakeElementsCommuteGeneral := function(S, range, Xe, Ze)
    local k;

    for k in range do
        S[k] := MakeElementCommute(S[k], Xe, Ze);
    od;

    return S;
    end;

    Args:
        matrix ([type]): [description]
        range ([type]): [description]
        hyper1 ([type]): [description]
        hyper2 ([type]): [description]

    Returns:
        [type]: [description]
    """

    new_matrix = matrix.copy()

    for index in range:
        new_matrix[index] = _make_element_commute_with_hyper_pair(new_matrix[index], hyper1, hyper2)

    return new_matrix

def find_noncommutative_partner(matrix, vector):
    """[summary]

    Args:
        matrix ([type]): [description]
        vector ([type]): [description]

    Returns:
        [type]: [description]

    # Find generator in S that does not commute with operator el
    FindNonCommunativePartner := function(S, el)
    local i;

    for i in [1..Size(S)] do
        if SymplecticProduct(S[i],el)=1 then
        return [S[i],i];
        fi;
    od;
    return [];
    end;
    """

    for index, item in enumerate(matrix):
        if _symplectic_product_vv(item, vector) == 1: 
            return (item, index)
    return None

    
