import numpy as np
import itertools
import functools

def distance_test(stab, logicOp, wt):
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
    pow2 = np.array([2 ** i for i in range(m + 1)], dtype=int)
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
    mask1 = 2 ** m
    T1c = []
    T1a = []
    if w1 > 0:
        for err in itertools.combinations(single_qubit_synd, w1):
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
            for err in itertools.combinations(single_qubit_synd, w2):
                synd = functools.reduce(lambda i, j: int(i) ^ int(j), err)
                if synd & mask1:
                    T2a.append(synd ^ mask1)
                else:
                    T2c.append(synd)
    else:
        T2c = T1c
        T2a = T1a

    return any(elem in T1c for elem in T2a) or any(elem in T1a for elem in T2c)
