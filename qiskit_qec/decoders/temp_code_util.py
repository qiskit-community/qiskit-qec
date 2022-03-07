"""Temporary module with methods for codes."""
from typing import List


def temp_syndrome(bitstring: List[int], operators: List[List[int]]) -> List[int]:
    """Compute the syndrome of a bit string.

    The operators are given as lists of supports. Negative values ignored.
    """
    syndrome = []
    for s in operators:
        value = 0
        for i in s:
            if i >= 0:
                value += bitstring[i]
        syndrome.append(value % 2)
    return syndrome
