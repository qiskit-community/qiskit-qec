from time import perf_counter
from typing import List, Tuple
import itertools
import numpy as np
from qiskit import QiskitError
from numba import njit, jit
from qiskit_qec.analysis.extensions import _c_solve
from qiskit_qec.analysis.extensions import _c_solve_sparse, _c_solve_optimal_sparse

# Solvers

def solve(A: np.array, b: np.array, minimize_weight: bool = False, max_dof: int = 9, stats: bool = False, lse_solver: str = None):
    if not is_binary(A) or not is_binary(b):
        raise ValueError('A and b must be binary arrays')
    if lse_solver is None:
        return solve_python(A=A, b=b, minimize_weight=minimize_weight, max_dof=max_dof, stats=stats)
    elif lse_solver == 'python':
        return solve_python(A=A, b=b, minimize_weight=minimize_weight, max_dof=max_dof, stats=stats)
    elif lse_solver == 'numba':
        return solve_numba(A=A, b=b, minimize_weight=minimize_weight, max_dof=max_dof, stats=stats)
    elif lse_solver == 'c++':
        return solve_cpp(A=A, b=b, minimize_weight=minimize_weight, max_dof=max_dof, stats=stats)
    else:
        raise ValueError(f'lse_solve has to be one of [None, "python", "numba", "c++"], not "{lse_solver}"')

def solve_cpp(A: np.array, b: np.array, minimize_weight: bool = False, max_dof: int = 9, stats: bool = False):
    t0 = perf_counter()
    A_sparse = [np.where(row)[0] for row in A]
    n = A.shape[1]
    t1 = perf_counter()
    if minimize_weight:
        solvable, x_sparse, dof, minimized = _c_solve_optimal_sparse(A_sparse, n, b, max_dof)
    else:
        solvable, x_sparse = _c_solve_sparse(A_sparse, n, b)
        dof = np.nan
        minimized = False

    t2 = perf_counter()
    if solvable:
        x = np.zeros(n, dtype=bool)
        x[x_sparse] = True
        t3 = perf_counter()

        t_to_sparse = t1 - t0
        t_solve_sparse = t2 - t1
        t_from_sparse = t3 - t2

    else:
        x = None

        t_to_sparse = t1 - t0
        t_solve_sparse = t2 - t1
        t_from_sparse = 0

    if stats:
        stats = {'t_to_sparse': t_to_sparse,
                 't_solve_sparse': t_solve_sparse,
                 't_from_sparse': t_from_sparse,
                 'dof': dof,
                 'minimized': minimized}
        return solvable, x, stats
    return solvable, x

def solve_python(A: np.array, b: np.array, minimize_weight: bool = False, max_dof: int = 9, stats: bool = False):
    t0 = perf_counter()
    try:
        Ap, bp = gaussian_elimination_python(A, b)
        solvable = True
    except LinAlgError:
        solvable = False
    t1 = perf_counter()
    if not solvable:
        if stats:
            stats = {'t_gauss_el': t1-t0,
                    't_backsub': 0,
                    't_weight_minimization': 0,
                    'dof': np.nan,
                    'minimized': False}
            return False, None, stats
        return False, None
    x_part = _back_substitution_python(Ap, bp)
    t2 = perf_counter()

    if minimize_weight:
        x, dof, minimized = _minimize_weight_python(x_part, Ap, max_dof)
        t_weight_minimization = perf_counter() - t2

    else:
        x = x_part
        dof = np.nan
        minimized = False
        t_weight_minimization = 0

    t_gauss_el = t1 - t0
    t_backsub = t2 - t1

    if stats:
        stats = {'t_gauss_el': t_gauss_el,
                 't_backsub': t_backsub,
                 't_weight_minimization': t_weight_minimization,
                 'dof': dof,
                 'minimized': minimized}
        return solvable, x, stats
    return solvable, x

def solve_numba(A: np.array, b: np.array, minimize_weight: bool = False, max_dof: int = 9, stats: bool = False):
    t0 = perf_counter()
    Ap, bp = gaussian_elimination_numba(A, b)
    if np.array_equal(Ap, np.empty((0, A.shape[1]), dtype=np.bool_)): # not solvable
        if stats:
            stats = {'t_gauss_el': perf_counter()-t0,
                    't_backsub': 0,
                    't_weight_minimization': 0,
                    'dof': np.nan,
                    'minimized': False}
            return False, None, stats
        return False, None
    
    t1 = perf_counter()
    x_part = _back_substitution_numba(Ap, bp)
    t2 = perf_counter()

    if minimize_weight:
        x, dof, minimized = _minimize_weight_numba(x_part, Ap, max_dof)
        t_weight_minimization = perf_counter() - t2

    else:
        x = x_part
        dof = np.nan
        minimized = False
        t_weight_minimization = 0

    t_gauss_el = t1 - t0
    t_backsub = t2 - t1

    if stats:
        stats = {'t_gauss_el': t_gauss_el,
                 't_backsub': t_backsub,
                 't_weight_minimization': t_weight_minimization,
                 'dof': dof,
                 'minimized': minimized}
        return True, x, stats
    return True, x

    # Ap, bp = gaussian_elimination_numba(A, b)
    # if np.array_equal(Ap, np.empty((0, A.shape[1]), dtype=np.bool_)):
    #     raise LinAlgError('System has no solution')
    
    # t0 = perf_counter()
    # error, t_part, nullity = _back_substitution_weight_opt(Ap,bp) #, A.astype(int), b.astype(int)
    # t_opt = perf_counter() - t0 - t_part

    # return error, (t_part, nullity, t_opt)

# Helper functions

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
    
    Ap, bp = gaussian_elimination_numba(A, b)
    if np.array_equal(Ap, np.empty((0, A.shape[1]), dtype=np.bool_)):
        raise LinAlgError('System has no solution')
    
    t0 = perf_counter()
    error, t_part, nullity = _back_substitution_weight_opt(Ap,bp) #, A.astype(int), b.astype(int)
    t_opt = perf_counter() - t0 - t_part

    return error, (t_part, nullity, t_opt)

# gaussian elimination section

@jit(nopython=True)
def check_all_zero_rows(A):
    """Manually checks each row for all zeros, compatible with Numba."""
    m = A.shape[0]
    all_zero_rows = np.zeros(m, dtype=np.bool_)
    for i in range(m):
        all_zero_rows[i] = np.all(A[i] == 0)
    return all_zero_rows

@jit(nopython=True)
def gaussian_elimination_numba(A, b):
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

def gaussian_elimination_python(A, b):
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

    return A, b

# deprecated
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

    Ap, bp = gaussian_elimination_python(A,b) 

    t0 = perf_counter()
    error, t_part, nullity = _back_substitution_weight_opt(Ap,bp) #, A.astype(int), b.astype(int)
    t_opt = perf_counter() - t0 - t_part

    return error, (t_part, nullity, t_opt)

# deprecated
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

# deprecated
def solve2(A: np.array, b: np.array, lse_solver: str = None):
    if lse_solver is None:
        return solve2_python(A=A, b=b)
    elif lse_solver == 'python':
        return solve2_python(A=A, b=b)
    elif lse_solver == 'numba':
        return solve2_numba(A=A, b=b)
    elif lse_solver == 'c++':
        return solve2_cpp(A=A, b=b)
    else:
        raise ValueError(f'lse_solve has to be one of [None, "python", "numba", "c++"], not "{lse_solver}"')

# back substitution section

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

def _back_substitution_python(A, b):
    return _back_substitution(A, b)

def _back_substitution_numba(A, b):
    return _back_substitution(A, b)

# deprecated
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
    return xs, 0, 0

    best = xs
    min_weight = xs.sum()
    if ker.shape[0] > 9:
        return best, t_part, ker.shape[0]
    for sel in itertools.product([False,True], repeat=ker.shape[0]):
        x = (xs + ker[list(sel)].sum(axis=0)) % 2
        if x.sum() < min_weight:
            best = x
            min_weight = x.sum()
    
    return best, t_part, ker.shape[0]

# minimizing weight section

def _minimize_weight(x_part, A, max_dof):
    ker = ker2(A) # basis of nullspace of A 
    # all solutions of Ax = b can be written as xs + ker

    best = x_part
    min_weight = x_part.sum()
    if ker.shape[0] > max_dof:
        return best, ker.shape[0], False
    for sel in itertools.product([False,True], repeat=ker.shape[0]):
        x = (x_part + ker[list(sel)].sum(axis=0)) % 2
        if x.sum() < min_weight:
            best = x
            min_weight = x.sum()
    
    return best, ker.shape[0], True

def _minimize_weight_python(x_part, A, max_dof):
    return _minimize_weight(x_part, A, max_dof)

def _minimize_weight_numba(x_part, A, max_dof):
    return _minimize_weight(x_part, A, max_dof)