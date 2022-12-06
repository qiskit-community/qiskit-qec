"""Test xp pauli rep."""

from unittest import TestCase
import numpy as np

from qiskit_qec.utils.xp_pauli_rep import (
    xp_symplectic2str,
    INDEX_SYNTAX,
    PRODUCT_SYNTAX,
    LATEX_SYNTAX,
    XP_SYMPLECTIC_SYNTAX
)

class TestXPPauliRep(TestCase):
    """Test xp pauli rep."""

    def test_xp_symplectic2str(self):
        """Tests xp_symplectic2str function."""

        precision = 8
        matrix1 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0], dtype=np.int64)
        phase_exp1 = 12
        matrix2 = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0], dtype=np.int64)
        phase_exp2 = 2
        matrix = np.array([matrix1, matrix2])
        phase_exp = np.array([phase_exp1, phase_exp2])

        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision), np.array(['XP8((w,12)(XP4)2(X)1(X)0)', 'XP8((w,2)(P3)3(XP2)2(X)1(X)0)']))
        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision, qubit_order="left-to-right"), np.array(['XP8((w,12)(X)0(X)1(XP4)2)', 'XP8((w,2)(X)0(X)1(XP2)2(P3)3)']))
        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision, syntax=XP_SYMPLECTIC_SYNTAX), np.array(['XP8(12|1 1 1 0 0 0 0|0 0 4 0 0 0 0)', 'XP8(2|1 1 1 0 0 0 0|0 0 2 3 0 0 0)']))
        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision, no_phase=True), np.array(['XP8((XP4)2(X)1(X)0)', 'XP8((P3)3(XP2)2(X)1(X)0)']))
        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision, syntax=PRODUCT_SYNTAX), np.array(['XP8((w,12)(I)(I)(I)(I)(XP4)(X)(X))', 'XP8((w,2)(I)(I)(I)(P3)(XP2)(X)(X))']))
        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision, syntax=LATEX_SYNTAX), np.array(['XP_{8}((w,12)(XP^{4})_{2}(X)_{1}(X)_{0})', 'XP_{8}((w,2)(P^{3})_{3}(XP^{2})_{2}(X)_{1}(X)_{0})']))

        # Tests conversion to unique vector format for different formats
        precision = 4
        matrix1 = np.array([1, 2, 3, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0], dtype=np.int64)
        matrix2 = np.array([1, 2, 3, 0, 0, 0, 0, 0, 0, 2, 6, 0, 0, 0], dtype=np.int64)
        matrix = np.array([matrix1, matrix2])

        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision), np.array(['XP4((w,4)(XP)2(X)0)', 'XP4((w,2)(P2)3(XP2)2(X)0)']))
        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision, qubit_order="left-to-right"), np.array(['XP4((w,4)(X)0(XP)2)', 'XP4((w,2)(X)0(XP2)2(P2)3)']))
        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision, syntax=XP_SYMPLECTIC_SYNTAX), np.array(['XP4(4|1 0 1 0 0 0 0|0 0 1 0 0 0 0)', 'XP4(2|1 0 1 0 0 0 0|0 0 2 2 0 0 0)']))
        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision, no_phase=True), np.array(['XP4((XP)2(X)0)', 'XP4((P2)3(XP2)2(X)0)']))
        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision, syntax=PRODUCT_SYNTAX), np.array(['XP4((w,4)(I)(I)(I)(I)(XP)(I)(X))', 'XP4((w,2)(I)(I)(I)(P2)(XP2)(I)(X))']))
        np.testing.assert_equal(xp_symplectic2str(matrix, phase_exp, precision, syntax=LATEX_SYNTAX), np.array(['XP_{4}((w,4)(XP)_{2}(X)_{0})', 'XP_{4}((w,2)(P^{2})_{3}(XP^{2})_{2}(X)_{0})']))

