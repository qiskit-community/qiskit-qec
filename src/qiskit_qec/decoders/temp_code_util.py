"""Temporary module with methods for codes."""
from functools import partial
from typing import List


def temp_gauge_products(stabilizers: List[List[int]], gauges: List[List[int]]) -> List[List[int]]:
    """Compute product of gauge operators for each stabilizer.

    stabilizers = list of stabilizer operator supports
    gauges = list of gauge operator supports

    For each stabilizer operator, record the indices of the
    gauge operators in 'gauges' whose product is that stabilizer operator.

    Return the list of indices.
    """
    gauge_products = []
    for _, stab in enumerate(stabilizers):
        # is_contained = lambda x, j: set(gauges[j]).intersection(set(x)) == set(gauges[j])
        def is_contained(x, j):
            return set(gauges[j]).intersection(set(x)) == set(gauges[j])

        products = filter(partial(is_contained, stab), range(len(gauges)))
        gauge_products.append(list(products))
    return gauge_products


def temp_syndrome(bitstring: List[int], operators: List[List[int]]) -> List[int]:
    """Compute the syndrome of a bit string.

    The operators are given as lists of supports. Negative values ignored.
    """
    syndrome = []
    for supp in operators:
        value = 0
        for i in supp:
            if i >= 0:
                value += bitstring[i]
        syndrome.append(value % 2)
    return syndrome
