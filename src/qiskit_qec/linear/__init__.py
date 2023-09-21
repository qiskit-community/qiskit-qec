# This code is part of Qiskit.
#
# (C) Copyright IBM 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
=================================
Linear (:mod:`qiskit_qec.linear`)
=================================

.. currentmodule:: qiskit_qec.linear


Linear module classes and functions
===================================

.. autosummary::
    :toctree: ../stubs/

    Bit
    SMatrix

Linear matrix functions
=======================

.. autosummary::
    :toctree: ../stubs/

    create_lambda_matrix
    augment_mat
    rref
    rank
    rref_complete

Linear symplectic functions
===========================

.. autosummary::
    :toctree: ../stubs/

    all_commute
    symplectic_product
    make_commute_hyper
    locate_hyper_partner
    build_hyper_partner
    symplectic_gram_schmidt
    is_symplectic_matrix_form
    is_symplectic_vector_form
    is_symplectic_form
    is_hyper_form
    is_center
    is_same_span
    is_stabilizer_group
"""

from .matrix import create_lambda_matrix, augment_mat, rref, rank, rref_complete
from .symplectic import (
    all_commute,
    symplectic_product,
    make_commute_hyper,
    locate_hyper_partner,
    build_hyper_partner,
    symplectic_gram_schmidt,
    is_symplectic_matrix_form,
    is_symplectic_vector_form,
    is_symplectic_form,
    is_hyper_form,
    is_center,
    is_same_span,
    is_stabilizer_group,
)
