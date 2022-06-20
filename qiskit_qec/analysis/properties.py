"""Compute code properties."""
from typing import List
from itertools import combinations, product
import numpy as np

from qiskit_qec.linear.matrix import rank
from qiskit_qec.linear.symplectic import all_commute, is_stabilizer_group


def min_distance_helper(
    n: int,
    n_minus_k_plus_r: int,
    k: int,
    pauli: List[str],
    max_weight: int,
    stabilizer: np.ndarray,
    gauge: np.ndarray,
):
    """Helper for min distance."""
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


def minimum_distance(stabilizer: np.ndarray, gauge: np.ndarray = None, max_weight: int = 10) -> int:
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
    distance = min_distance_helper(n, n_minus_k_plus_r, k, pauli, max_weight, stabilizer, gauge)
    return distance
