"""Compute code properties."""
from typing import List
from itertools import combinations, product
import functools
import logging
import numpy as np

from qiskit_qec.linear.matrix import rank
from qiskit_qec.linear.symplectic import all_commute, is_stabilizer_group

logger = logging.getLogger(__name__)


def attempt_import(module_name: str):
    """Try to load a module."""
    try:
        __import__(module_name)
    except ImportError as e:  # pylint: disable=invalid-name
        logger.exception(  # pylint: disable=logging-fstring-interpolation
            f"__import__({module_name}) failed, raising {e}"
        )
        return False
    else:
        return True


# TODO: wrap minimum distance test into a min. dist. method
# TODO: from symplectic import normalizer
# TODO: implement distance_test in C


def distance_test(stab: np.ndarray, logicOp: np.ndarray, wt: int) -> bool:
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
    if (wt % 2) == 0:
        w1 = int(wt / 2)
        w2 = int(wt / 2)
    else:
        w1 = int((wt + 1) / 2)
        w2 = int((wt - 1) / 2)
    assert (w1 + w2) == wt

    # Compute syndromes of all single-qubit errors
    single_qubit_synd = []
    for q in range(n):
        # single-qubit X error
        syndX = np.append(stab[:, n + q], logicOp[n + q])
        single_qubit_synd.append(np.dot(pow2, syndX))
        # single-qubit Z error
        syndZ = np.append(stab[:, q], logicOp[q])
        single_qubit_synd.append(np.dot(pow2, syndZ))
        # single-qubit Y error
        single_qubit_synd.append(np.dot(pow2, syndX) ^ np.dot(pow2, syndZ))

    # examine all errors with the weight w1
    mask1 = 2**m
    T1c = []
    T1a = []
    if w1 > 0:
        for err in combinations(single_qubit_synd, w1):
            synd = functools.reduce(lambda i, j: int(i) ^ int(j), err)
            if synd & mask1:
                T1a.append(synd ^ mask1)
            else:
                T1c.append(synd)
    T1c = set(T1c)
    T1a = set(T1a)

    if not (w1 == w2):
        # examine all errors with the weight w2
        T2c = []
        T2a = []
        if w2 > 0:
            for err in combinations(single_qubit_synd, w2):
                synd = functools.reduce(lambda i, j: int(i) ^ int(j), err)
                if synd & mask1:
                    T2a.append(synd ^ mask1)
                else:
                    T2c.append(synd)
    else:
        T2c = T1c
        T2a = T1a

    return any(elem in T1c for elem in T2a) or any(elem in T1a for elem in T2c)


def min_distance_core(
    n: int,
    n_minus_k_plus_r: int,
    k: int,
    pauli: List[str],
    max_weight: int,
    stabilizer: np.ndarray,
    gauge: np.ndarray,
):
    """Core of minimum distance algorithm."""
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
                commutes = all_commute(test_matrix)
                in_gauge = rank(test_matrix_2) == n_minus_k_plus_r
                if commutes:
                    if (k > 0 and not in_gauge) or (k == 0 and in_gauge):
                        distance = weight
                        return distance
    return 0


def minimum_distance_python(
    stabilizer: np.ndarray, gauge: np.ndarray = None, max_weight: int = 10
) -> int:
    """Minimum distance of (subsystem) stabilizer code.

    stabilizer is a symplectic matrix generating the stabilizer group.
    gauge is a symplectic matrix generating the gauge group, if any.

    Returns the minimum distance of the code, or 0 if greater than max_weight.
    """
    assert is_stabilizer_group(stabilizer)
    if gauge is None:
        gauge = stabilizer
    n = int(stabilizer.shape[1] / 2)
    n_minus_k_minus_r = rank(stabilizer)
    n_minus_k_plus_r = rank(gauge)
    n_minus_k = (n_minus_k_minus_r + n_minus_k_plus_r) / 2
    k = n - n_minus_k
    pauli = ["x", "y", "z"]
    distance = min_distance_core(n, n_minus_k_plus_r, k, pauli, max_weight, stabilizer, gauge)
    return distance


def minimum_distance(stabilizer: np.ndarray, gauge: np.ndarray = None, max_weight: int = 10) -> int:
    """Minimum distance of (subsystem) stabilizer code.

    stabilizer is a symplectic matrix generating the stabilizer group.
    gauge is a symplectic matrix generating the gauge group, if any.

    Returns the minimum distance of the code, or 0 if greater than max_weight.
    """
    if gauge is None:
        gauge = stabilizer
    if attempt_import("qiskit_qec.extensions.compiledextension"):
        from qiskit_qec.extensions import compiledextension

        inputform1 = stabilizer.astype(np.int32).tolist()
        inputform2 = gauge.astype(np.int32).tolist()
        distance = compiledextension.minimum_distance(inputform1, inputform2, max_weight)
    else:
        distance = minimum_distance_python(stabilizer, gauge, max_weight)
    return distance
