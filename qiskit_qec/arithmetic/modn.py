# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
#
# This code is adapted from XPFpackage: https://github.com/m-webster/XPFpackage, originally developed by Mark Webster. The original code is licensed under the GNU General Public License v3.0 and Mark Webster has given permission to use the code under the Apache License v2.0.

"""Modular arithmetic Z/nZ."""

from typing import Tuple
import math
import numpy as np
from qiskit import QiskitError


def gcd_ext(a: int, b: int, n: int) -> Tuple[int, int, int, int, int]:
    """Implements the extended Euclidean algorithm in the ring Z/nZ: for any two integers a & b, find g, s, t, u, v that satisfy

    1. g = gcd_n(a, b) = s * a + t * b, where gcd_n stands for the greatest common divisor in the ring Z/nZ, s & t are called the Bezout coefficients for a & b;
    2. (u * a + v * b) mod n = 0;
    3. (s * v - t * u) mod n = 1.

    Args:
        a, b: input integers
        n: modulus

    Returns:
        g, s, t, u, v: g is the greatest common divisor of (a mod n) and (b mod n); s, t, u, and v are integer coefficients satisfying (s*a + t*b) mod n = g, (u*a + v*b) mod n = 0, (s*v - t*u) mod n = 1

    Raises:
        QiskitError: Input must be integers.
        QiskitError: n must be a positive integer.

    Examples:
        >>> gcd_ext(15, 6, 4)
        (1, 1, -1, -2, 3)

    See Also: 
    _gcd_ext
    """
    if not isinstance(a, (int, np.integer)) or not isinstance(b, (int, np.integer)):
        raise QiskitError("Input must be integers.")
    if not n > 0 or not isinstance(n, (int, np.integer)):
        raise QiskitError("n must be a positive integer")

    return _gcd_ext(a, b, n)


def _gcd_ext(a: int, b: int, n: int) -> Tuple[int, int, int, int, int]:
    """Implements the extended Euclidean algorithm in the ring Z/nZ: for any two integers a & b, find g, s, t, u, v that satisfy

    1. g = gcd_n(a, b) = s * a + t * b, where gcd_n stands for the greatest common divisor in the ring Z/nZ, s & t are called the Bezout coefficients for a & b;
    2. (u * a + v * b) mod n = 0;
    3. (s * v - t * u) mod n = 1.

    Args:
        a, b: input integers
        n: modulus

    Returns:
        g, s, t, u, v: g is the greatest common divisor of (a mod n) and (b mod n); s, t, u, and v are integer coefficients satisfying (s*a + t*b) mod n = g, (u*a + v*b) mod n = 0, (s*v - t*u) mod n = 1

    Examples:
        >>> _gcd_ext(15, 6, 4)
        (1, 1, -1, -2, 3)

    See Also:
    gcd_ext
    """
    old_s, s  = 1, 0
    old_t, t  = 0, 1
    old_r, r  = a % n, b % n

    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t

    p = int(math.copysign(1, t * old_s - s * old_t))
    u, v = p * s, p * t
    g, s, t = old_r, old_s, old_t

    return (g, s, t, u, v)


def quo(a: int, b: int, n: int) -> int:
    """Computes the quotient of a/b in the ring Z/nZ, i.e. returns integer q such that a = (b * q) mod n. Returns None if b mod n = 0.

    Args:
        a: numerator
        b: denominator
        n: modulus

    Returns:
        quotient of a/b in the ring Z/nZ

    Raises:
        QiskitError: Input must be integers.
        QiskitError: n must be a positive integer.

    Examples:
        >>> quo(25, 5, 4)
        1
        
        >>> quo(25, 8, 4)
        None

    See Also:
    _quo
    """
    if not isinstance(a, (int, np.integer)) or not isinstance(b, (int, np.integer)):
        raise QiskitError("Input must be integers.")
    if not n > 0 or not isinstance(n, (int, np.integer)):
        raise QiskitError("n must be a positive integer")

    return _quo(a, b, n)


def _quo(a: int, b: int, n: int) -> int:
    """Computes the quotient of a/b in the ring Z/nZ, i.e. returns integer q such that a = (b * q) mod n. Returns None if b mod n = 0.

    Args:
        a: numerator
        b: denominator
        n: modulus

    Returns:
        quotient of a/b in the ring Z/nZ

    Examples:
        >>> _quo(25, 5, 4)
        1

        >>> _quo(25, 8, 4)
        None

    See Also:
    quo
    """
    a, b = a % n, b % n
    if b == 0:
        return None
    return (a // b) % n


def div(a: int, b: int, n: int) -> int:
    """Computes the divisor of a/b in the ring Z/nZ, i.e., returns integer d such that b * d = a mod n. Returns None if no such d exists.

    Args:
        a: numerator
        b: denominator
        n: modulus

    Returns:
        divisor of a/b in the ring Z/nZ

    Raises:
        QiskitError: Input must be integers.
        QiskitError: n must be a positive integer.

    Examples:
        >>> div(24, 8, 5)
        3

        >>> div(24, 8, 4)
        None

        >>> div(23, 10, 4)
        None

    See Also:
    _div
    """
    if not isinstance(a, (int, np.integer)) or not isinstance(b, (int, np.integer)):
        raise QiskitError("Input must be integers.")
    if not n > 0 or not isinstance(n, (int, np.integer)):
        raise QiskitError("n must be a positive integer")

    return _div(a, b, n)


def _div(a: int, b: int, n: int) -> int:
    """Computes the divisor of a/b in the ring Z/nZ, i.e., returns integer d such that b * d = a mod n. Returns None if no such d exists.

    Args:
        a: numerator
        b: denominator
        n: modulus

    Returns:
        divisor of a/b in the ring Z/nZ

    Examples:
        >>> _div(24, 8, 5)
        3

        >>> _div(24, 8, 4)
        None

        >>> _div(23, 10, 4)
        None

    See Also:
    div
    """
    a, b = a % n, b % n
    if b == 0:
        return None
    g = math.gcd(b, n)
    if a % g != 0:
        return None
    else:
        r = a % b
        while r > 0:
            a += n
            r = a % b
        return a // b % n


def ann(a: int, n: int) -> int:
    """Computes the annihilator of a in the ring Z/nZ, i.e., returns integer b such that (b * a) mod N = 0.

    Args:
        a: input integer
        n: modulus

    Returns:
        annihilator of a in the ring Z/nZ

    Raises:
        QiskitError: Input must be integers.
        QiskitError: n must be a positive integer.

    Examples:
        >>> ann(10, 5)
        1

        >>> ann(3, 5)
        0

    See Also:
    _ann
    """
    if not isinstance(a, (int, np.integer)):
        raise QiskitError("Input must be integers.")
    if not n > 0 or not isinstance(n, (int, np.integer)):
        raise QiskitError("n must be a positive integer")

    return _ann(a, n)


def _ann(a: int, n: int) -> int:
    """Computes the annihilator of a in the ring Z/nZ, i.e., returns integer b such that (b * a) mod n = 0.

    Args:
        a: input integer
        n: modulus

    Returns:
        annihilator of a in the ring Z/nZ

    Examples:
        >>> _ann(10, 5)
        1

        >>> _ann(3, 5)
        0

    See Also:
    ann
    """
    a = a % n
    if a == 0:
        return 1
    u = n // math.gcd(a, n)
    return u % n


def stab(a: int, b: int, n: int) -> int:
    """Returns a ring element c such that gcd(a + b * c, n) = gcd(a, b, n) in the ring Z/nZ.

    Args:
        a, b: input integers
        n: modulus

    Returns:
        ring element c such that gcd(a+b*c, N) = gcd(a, b, N)

    Raises:
        QiskitError: Input must be integers.
        QiskitError: n must be a positive integer.

    Examples:
        >>> stab(25, 8, 6)
        0

        >>> stab(24, 8, 6)
        1

    See Also:
    _stab
    """
    if not isinstance(a, (int, np.integer)) or not isinstance(b, (int, np.integer)):
        raise QiskitError("Input must be integers.")
    if not n > 0 or not isinstance(n, (int, np.integer)):
        raise QiskitError("n must be a positive integer")

    return _stab(a, b, n)


def _stab(a: int, b: int, n: int) -> int:
    """Returns a ring element c such that gcd(a + b * c, n) = gcd(a, b, n) in the ring Z/nZ.

    Args:
        a, b: input integers
        n: modulus

    Returns:
        ring element c such that gcd(a+b*c, N) = gcd(a, b, N)

    Examples:
        >>> _stab(25, 8, 6)
        0

        >>> _stab(24, 8, 6)
        1

    See Also:
    stab
    """
    a, b = a % n, b % n
    g = math.gcd(math.gcd(a, b), n)
    n_old = n
    a, n = a // g, n // g
    if n == 0:
        c = 0
    else:
        a = a % n
        if a == 0:
            c = 1
        else:
            r = int(math.ceil(math.log2(math.log2(n)))) if n > 1 else 1
            for _ in range(r):
                a = a * a % n
            c = n // math.gcd(a, n)
    return c % n_old


def unit(a: int, n: int) -> int:
    """Computes a unit u such that for element a in the ring Z/nZ, i.e., (a * u) mod n = gcd(a, n).

    Args:
        a: input integer
        n: modulus

    Returns:
        unit of a in the ring Z/nZ

    Raises:
        QiskitError: Input must be integers.
        QiskitError: n must be a positive integer.

    Examples:
        >>> unit(10, 5)
        1

        >>> unit(3, 5)
        2

    See Also:
    _unit
    """
    if not isinstance(a, (int, np.integer)):
        raise QiskitError("Input must be integers.")
    if not n > 0 or not isinstance(n, (int, np.integer)):
        raise QiskitError("n must be a positive integer")

    return _unit(a, n)


def _unit(a: int, n: int) -> int:
    """Computes a unit u such that for element a in the ring Z/nZ, i.e., (a * u) mod n = gcd(a, n).

    Args:
        a: input integer
        n: modulus

    Returns:
        unit of a in the ring Z/nZ

    Examples:
        >>> _unit(10, 5)
        1

        >>> _unit(3, 5)
        2

    See Also:
    unit
    """
    a = a % n
    if a == 0:
        return 1
    g = math.gcd(a, n)
    s = div(g, a, n)
    if g == 1:
        return s
    d = stab(s, n // g, n)
    return (s + d * n // g) % n