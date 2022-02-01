
import numpy as np
from typing import List
from qiskit import QiskitError
from qiskit_qec.linear import matrix
from qiskit_qec.linear import symplectic

def find_hyperbolic_partner(check_matrix, index: int):
    check_matrix = np.atleast_2d(np.array(check_matrix,  dtype=bool))

    # check_matrix -> all associated operators must commute
    assert symplectic.all_commute(check_matrix) == True, \
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

    _lambda = matrix.create_lambda_matrix(ncols)
    slambda = np.matmul(check_matrix, _lambda)

    heads, rref_mat, transform_mat, rank = matrix._rref_complete(slambda)

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

