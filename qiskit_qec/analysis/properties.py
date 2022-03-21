"""Compute code properties."""
from itertools import combinations, product
import numpy as np

from qiskit_qec.linear.matrix import rank
from qiskit_qec.linear.symplectic import all_commute, is_stabilizer_group


def minimum_distance(stabilizer: np.ndarray) -> int:
    """Minimum distance of stabilizer code.

    stabilizer is symplectic matrix with full rank generating set.

    Returns the minimum distance of the code, or -1 if greater than max_weight.
    """
    max_weight = 10
    assert is_stabilizer_group(stabilizer)
    n = int(stabilizer.shape[1] / 2)
    n_minus_k = int(stabilizer.shape[0])
    assert n_minus_k == rank(stabilizer)
    k = n - n_minus_k
    pauli = ["x", "y", "z"]
    quit = False
    distance = -1
    try:
        for weight in range(1, n + 1):
            iterable = [pauli] * weight
            for combination in combinations(list(range(n)), weight):
                for paulis in product(*iterable):
                    error = np.zeros((2 * n,))
                    for i in range(weight):
                        if paulis[i] == "x" or paulis[i] == "y":
                            error[combination[i]] = 1
                        if paulis[i] == "z" or paulis[i] == "y":
                            error[combination[i] + n] = 1
                    test_matrix = np.vstack([stabilizer, error])
                    commutes = all_commute(test_matrix)
                    in_stabilizer = rank(test_matrix) == n_minus_k
                    if commutes:
                        if (k > 0 and not in_stabilizer) or (k == 0 and in_stabilizer):
                            distance = weight
                            raise StopIteration
    except StopIteration:
        pass
    return distance
