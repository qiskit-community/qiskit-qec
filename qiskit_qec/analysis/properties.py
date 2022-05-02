"""Compute code properties."""
from itertools import combinations, product
import numpy as np

from qiskit_qec.linear.matrix import rank
from qiskit_qec.linear.symplectic import all_commute, is_stabilizer_group
import logging


def minimum_distance(stabilizer: np.ndarray, max_weight: int = 10) -> int:
    logger = logging.getLogger(__name__)
    """Minimum distance of stabilizer code.

    stabilizer is a symplectic matrix generating the stabilizer group.

    Returns the minimum distance of the code, or 0 if greater than max_weight.
    """
    assert is_stabilizer_group(stabilizer)
    n = int(stabilizer.shape[1] / 2)
    n_minus_k = rank(stabilizer)
    k = n - n_minus_k
    pauli = ["x", "y", "z"]
    distance = 0
    try:
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
                    commutes = all_commute(test_matrix)
                    in_stabilizer = rank(test_matrix) == n_minus_k
                    if commutes:
                        if (k > 0 and not in_stabilizer) or (k == 0 and in_stabilizer):
                            distance = weight
                            raise StopIteration
    except StopIteration as e:
        logger.exception(f"Raised {e}")
    return distance
