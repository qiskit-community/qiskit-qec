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

# Part of the QEC framework

"""Sm numpy."""

import numpy as np


def asarray(*args, **kargs):
    """To array."""
    return np.asarray(*args, **kargs)


def atleast_2d(*args, **kargs):
    """View inputs as arrays with at least two dimensions."""
    return np.atleast_2d(*args, **kargs)


def concatenate(*args, **kargs):
    """Concatenate."""
    return np.concatenate(*args, **kargs)


def sum(*args, **kargs):  # pylint: disable=redefined-builtin
    """Sum."""
    return np.sum(*args, **kargs)


def zeros(*args, **kargs):
    """Return a new array of given shape and type, filled with zeros."""
    return np.zeros(*args, **kargs)
