"""Compute code distance."""
from typing import List
from itertools import combinations, product
import functools
import logging
import numpy as np

from qiskit_qec.analysis.extensions import C_MIN_DISTANCE, C_MIN_DISTANCE_BY_TESTS

from qiskit_qec.exceptions import QiskitQECError
from qiskit_qec.linear.matrix import rank
from qiskit_qec.linear import symplectic

if C_MIN_DISTANCE:
    from qiskit_qec.analysis.extensions import _c_minimum_distance, _c_minimum_distance_by_tests

logger = logging.getLogger(__name__)


# pylint: disable=too-many-locals
def _distance_test(stab: np.ndarray, logic_op: np.ndarray, weight: int) -> bool:
    """Tests if a low weight logical Pauli operator exists.

    Tests whether there exists a Pauli operator with Hamming weight <= wt
    that commutes with every row of the stabilizer table (stab) and
    anti-commutes with a given logical operator (logicOp).

    Pauli operators on n qubits are represented as integer numpy arrays
    of length 2n in the order [X-part, Z-part].

    Returns True or False.
    """
    # number of qubits
    n = int(stab.shape[1] / 2)
    m = stab.shape[0]
    pow2 = np.array([2**i for i in range(m + 1)], dtype=int)
    if (weight % 2) == 0:
        w1 = int(weight / 2)
        w2 = int(weight / 2)
    else:
        w1 = int((weight + 1) / 2)
        w2 = int((weight - 1) / 2)
    assert (w1 + w2) == weight

    # Compute syndromes of all single-qubit errors
    single_qubit_synd = []
    for q in range(n):
        # single-qubit X error
        syndx = np.append(stab[:, n + q], logic_op[n + q])
        single_qubit_synd.append(np.dot(pow2, syndx))
        # single-qubit Z error
        syndz = np.append(stab[:, q], logic_op[q])
        single_qubit_synd.append(np.dot(pow2, syndz))
        # single-qubit Y error
        single_qubit_synd.append(np.dot(pow2, syndx) ^ np.dot(pow2, syndz))

    # examine all errors with the weight w1
    mask1 = 2**m
    t1c = []
    t1a = []
    if w1 > 0:
        for err in combinations(single_qubit_synd, w1):
            synd = functools.reduce(lambda i, j: int(i) ^ int(j), err)
            if synd & mask1:
                t1a.append(synd ^ mask1)
            else:
                t1c.append(synd)
    t1c = set(t1c)
    t1a = set(t1a)

    if w1 != w2:
        # examine all errors with the weight w2
        t2c = []
        t2a = []
        if w2 > 0:
            for err in combinations(single_qubit_synd, w2):
                synd = functools.reduce(lambda i, j: int(i) ^ int(j), err)
                if synd & mask1:
                    t2a.append(synd ^ mask1)
                else:
                    t2c.append(synd)
    else:
        t2c = t1c
        t2a = t1a

    return any(elem in t1c for elem in t2a) or any(elem in t1a for elem in t2c)


def _minimum_distance_2_python(stabilizer: np.ndarray, gauge: np.ndarray, max_weight) -> int:
    """Minimum distance of (subsystem) stabilizer code.

    stabilizer is a symplectic matrix generating the stabilizer group.
    gauge is a symplectic matrix generating the gauge group, if any.

    Second method based on parititioning errors.

    Returns the minimum distance of the code, or 0 if greater than max_weight.
    """
    if symplectic.is_stabilizer_group(gauge):
        _, xl, zl = symplectic.normalizer(stabilizer.astype(bool))
    else:
        center, x, z = symplectic.symplectic_gram_schmidt(gauge)
        _, xp, zp = symplectic.normalizer(center, x, z)
        x_rows = x.view([("", x.dtype)] * x.shape[1])
        xp_rows = xp.view([("", xp.dtype)] * xp.shape[1])
        xl = np.setdiff1d(xp_rows, x_rows).view(xp.dtype).reshape(-1, xp.shape[1])
        z_rows = z.view([("", z.dtype)] * z.shape[1])
        zp_rows = zp.view([("", zp.dtype)] * zp.shape[1])
        zl = np.setdiff1d(zp_rows, z_rows).view(zp.dtype).reshape(-1, zp.shape[1])
    if xl.shape[0] == 0:  # k = 0, fall back to first method
        return _minimum_distance_1_python(stabilizer, gauge, max_weight)
    weight = max_weight + 1
    for row in range(xl.shape[0]):
        for w in range(1, max_weight + 1):
            if _distance_test(stabilizer.astype(int), xl[row].astype(int), w):
                weight = min(weight, w)
                break
    for row in range(zl.shape[0]):
        for w in range(1, max_weight + 1):
            if _distance_test(stabilizer.astype(int), zl[row].astype(int), w):
                weight = min(weight, w)
                break
    if weight < max_weight + 1:
        return weight
    else:
        return 0


def _minimum_distance_1_python_core(
    n: int,
    n_minus_k_plus_r: int,
    k: int,
    pauli: List[str],
    max_weight: int,
    stabilizer: np.ndarray,
    gauge: np.ndarray,
) -> int:
    """Core of minimum distance algorithm that enumerates errors."""
    for weight in range(1, min(n, max_weight) + 1):
        iterable = [pauli] * weight
        for combination in combinations(list(range(n)), weight):
            for paulis in product(*iterable):
                error = np.zeros((2 * n,))
                for i in range(weight):
                    if paulis[i] in ("x", "y"):
                        error[combination[i]] = 1
                    if paulis[i] in ("z", "y"):
                        error[combination[i] + n] = 1
                test_matrix = np.vstack([stabilizer, error])
                test_matrix_2 = np.vstack([gauge, error])
                commutes = symplectic.all_commute(test_matrix)
                in_gauge = rank(test_matrix_2) == n_minus_k_plus_r
                if commutes:
                    if (k > 0 and not in_gauge) or (k == 0 and in_gauge):
                        distance = weight
                        return distance
    return 0


def _minimum_distance_1_python(stabilizer: np.ndarray, gauge: np.ndarray, max_weight) -> int:
    """Minimum distance of (subsystem) stabilizer code.

    stabilizer is a symplectic matrix generating the stabilizer group.
    gauge is a symplectic matrix generating the gauge group.

    Method enumerates all errors up to weight d.

    Returns the minimum distance of the code, or 0 if greater than max_weight.
    """
    n = int(stabilizer.shape[1] / 2)
    n_minus_k_minus_r = rank(stabilizer)
    n_minus_k_plus_r = rank(gauge)
    n_minus_k = (n_minus_k_minus_r + n_minus_k_plus_r) / 2
    k = n - n_minus_k
    pauli = ["x", "y", "z"]
    distance = _minimum_distance_1_python_core(
        n, n_minus_k_plus_r, k, pauli, max_weight, stabilizer, gauge
    )
    return distance


def _minimum_distance_2_compiled(stabilizer: np.ndarray, gauge: np.ndarray, max_weight) -> int:
    """Minimum distance of (subsystem) stabilizer code.

    stabilizer is a symplectic matrix generating the stabilizer group.
    gauge is a symplectic matrix generating the gauge group, if any.

    Second method based on parititioning errors, using compiled implementation.

    Returns the minimum distance of the code, or 0 if greater than max_weight.
    """

    if symplectic.is_stabilizer_group(gauge):
        _, xl, zl = symplectic.normalizer(stabilizer.astype(bool))
    else:
        center, x, z = symplectic.symplectic_gram_schmidt(gauge)
        _, xp, zp = symplectic.normalizer(center, x, z)
        x_rows = x.view([("", x.dtype)] * x.shape[1])
        xp_rows = xp.view([("", xp.dtype)] * xp.shape[1])
        xl = np.setdiff1d(xp_rows, x_rows).view(xp.dtype).reshape(-1, xp.shape[1])
        z_rows = z.view([("", z.dtype)] * z.shape[1])
        zp_rows = zp.view([("", zp.dtype)] * zp.shape[1])
        zl = np.setdiff1d(zp_rows, z_rows).view(zp.dtype).reshape(-1, zp.shape[1])

    inputform1 = stabilizer.astype(np.int32).tolist()
    inputform1p = gauge.astype(np.int32).tolist()
    inputform2 = xl.astype(np.int32).tolist()
    inputform3 = zl.astype(np.int32).tolist()
    # pylint: disable=c-extension-no-member
    if xl.shape[0] == 0:  # k = 0, fall back to first method
        # pylint: disable=c-extension-no-member
        return _c_minimum_distance(inputform1, inputform1p, max_weight)
    else:
        if C_MIN_DISTANCE_BY_TESTS:
            return _c_minimum_distance_by_tests(  # pylint: disable=c-extension-no-member
                inputform1, inputform2, inputform3, max_weight
            )
        else:
            logger.exception("from _c_minimum_distance_by_tests did not load - check import logs.")
            raise ImportError("from _c_minimum_distance_by_tests did not load - check import logs.")


def minimum_distance(
    stabilizer_or_gauge: np.ndarray,
    max_weight: int = 10,
    method: str = "enumerate",
    try_compiled: bool = True,
) -> int:
    """Minimum distance of (subsystem) stabilizer code.

    stabilizer_or_gauge is a symplectic matrix generating
    either the stabilizer or gauge group.

    method and try_compiled select an algorithm and implementation.

    Returns the minimum distance of the code, or 0 if greater than max_weight.
    """
    method_enumerate: str = "enumerate"
    method_partition: str = "partition"
    available_methods = {method_enumerate, method_partition}
    if method not in available_methods:
        raise QiskitQECError(f"method {method} is not supported.")
    if symplectic.is_stabilizer_group(stabilizer_or_gauge):
        stabilizer = stabilizer_or_gauge
        gauge = stabilizer_or_gauge
    else:
        stabilizer = symplectic.center(stabilizer_or_gauge)
        gauge = stabilizer_or_gauge

    if C_MIN_DISTANCE and try_compiled is True:
        inputform1 = stabilizer.astype(np.int32).tolist()
        inputform2 = gauge.astype(np.int32).tolist()
        if method == method_enumerate:
            distance = _c_minimum_distance(inputform1, inputform2, max_weight)
        elif method == method_partition:
            distance = _minimum_distance_2_compiled(stabilizer, gauge, max_weight)
    else:
        logger.exception("from compiled extension was not loaded: switching to no compiled")
        if method == method_enumerate:
            distance = _minimum_distance_1_python(stabilizer, gauge, max_weight)
        elif method == method_partition:
            distance = _minimum_distance_2_python(stabilizer, gauge, max_weight)
    return distance
