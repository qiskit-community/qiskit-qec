"""Test pauli rep."""

from unittest import TestCase
import numpy as np

from qiskit_qec.utils.pauli_rep import (
    split_pauli_enc,
    change_pauli_encoding,
    stand_phase_str,
    is_scalar,
    squeeze,
    is_exp_type,
    cpxstr2cpx,
    cpx2cpxstr,
    exp2cpx,
    cpx2exp,
    expstr2exp,
    exp2expstr,
    exp2exp,
    cpxstr2expstr,
    expstr2cpxstr,
    cpxstr2exp,
    expstr2cpx,
    cpx2expstr,
    str2exp,
    split_pauli,
    encode_of_phase_str,
    encode_of_tensor_str,
    str2symplectic,
    symplectic2str,
    str2str,
    INDEX_SYNTAX,
    PRODUCT_SYNTAX,
    LATEX_SYNTAX,
)


class TestPauliRep(TestCase):
    """Test pauli rep."""

    def test_split_pauli_enc(self):
        """Tests split_pauli_enc function."""
        encoding = "iXZ"
        result = split_pauli_enc(encoding)
        expecting = ("i", "XZ")

        self.assertEqual(result, expecting)

    def test_change_pauli_encoding(self):
        """Tests change_pauli_encoding function."""
        result1 = change_pauli_encoding(
            2, y_count=1, input_pauli_encoding="-iXZY", output_pauli_encoding="-iXZ"
        )
        self.assertEqual(result1, 1)

        result2 = change_pauli_encoding(
            np.array([2, 1]),
            y_count=np.array([1, 3]),
            input_pauli_encoding="-iXZY",
            output_pauli_encoding="-isXZ",
        )
        self.assertTrue(np.array_equal(result2, np.array([[1, 0], [0, 1]])))

        result3 = change_pauli_encoding(
            np.array([1, 1]),
            y_count=1,
            input_pauli_encoding="isXZY",
            output_pauli_encoding="iXZ",
        )
        self.assertEqual(result3, 0)

    def test_stand_phase_str(self):
        """Tests stand_phase_str function."""
        self.assertEqual(stand_phase_str("1j"), "i")
        self.assertEqual(stand_phase_str("(-1j,0)(-1,1)"), "(-i,0)(-1,1)")

    def test_is_scalar(self):
        """Tests is_scalar function."""
        self.assertFalse(is_scalar(np.array([1, 2, 3])))
        self.assertTrue(is_scalar(0 - 1j))
        self.assertTrue(is_scalar("xxxyxyii"))
        self.assertTrue(is_scalar([1, 0], scalar_pair=True))

    def test_squeeze(self):
        """Tests squeeze function."""
        self.assertTrue(np.array_equal(squeeze(np.array([1, 2, 3])), np.array([1, 2, 3])))

        self.assertEqual(squeeze(np.array(1), scalar=True), 1)

        self.assertTrue(np.array_equal(squeeze(np.array(1)), np.array(1)))

        self.assertTrue(np.array_equal(squeeze(np.array([[1, 2, 3]])), np.array([1, 2, 3])))

        self.assertTrue(np.array_equal(squeeze(np.array([[[[1, 2]]]])), np.array([1, 2])))

    def test_is_exp_type(self):
        """Tests is_exp_type function."""
        self.assertFalse(is_exp_type(np.array([1]), "is"))
        self.assertTrue(is_exp_type(np.array([1]), "i"))

    def test_cpxstr2cpx(self):
        """Tests cpxstr2cpx function."""
        self.assertEqual(cpxstr2cpx("+1j"), 1j)
        self.assertEqual(cpxstr2cpx("-1i"), (-0 - 1j))

    def test_cpx2cpxstr(self):
        """Tests cpx2cpxstr function."""
        self.assertEqual(cpx2cpxstr(1j), "i")
        self.assertTrue(
            np.array_equal(
                cpx2cpxstr([1j, +1j, -1, 1j, (0 - 1j)]),
                np.array(["i", "i", "-", "i", "-i"], dtype="<U2"),
            )
        )
        self.assertTrue(
            np.array_equal(
                cpx2cpxstr([1j, +1j, -1, 1j, (0 - 1j)], ones=True),
                np.array(["i", "i", "-1", "i", "-i"], dtype="<U2"),
            )
        )
        self.assertTrue(
            np.array_equal(cpx2cpxstr(1j, same_type=False), np.array(["i"], dtype="<U2"))
        )

    def test_exp2cpx(self):
        """Tests exp2cpx function."""
        self.assertTrue(np.array_equal(exp2cpx([[0, 1]], "-is"), np.array([-1.0 + 0.0j])))

        self.assertTrue(np.array_equal(exp2cpx([0, 1], "-i"), np.array([1.0 + 0.0j, -0.0 - 1.0j])))

        self.assertTrue(
            np.array_equal(
                exp2cpx(np.array([[0, 0], [1, 0]]), "is"),
                np.array([1.0 + 0.0j, 0.0 + 1.0j]),
            )
        )

    def test_cpx2exp(self):
        """Tests cpx2exp function."""
        self.assertEqual(cpx2exp(1j, "-i"), 3)
        self.assertTrue(np.array_equal(cpx2exp(1j, "-is"), np.array([1, 1])))

    def test_expstr2exp(self):
        """Tests expstr2exp function."""
        self.assertTrue(np.array_equal(expstr2exp("(-i,1)(-1,1)"), np.array([1, 1])))

        self.assertTrue(
            np.array_equal(
                expstr2exp(np.array(["(-i,1)", "(-i,2)", "(-i,0)"])),
                np.array([1, 2, 0]),
            )
        )

        self.assertEqual(expstr2exp("(i,3)"), 3)

    def test_exp2expstr(self):
        """Tests exp2expstr function."""
        self.assertEqual(exp2expstr([0, 1], "-is"), "(-i,0)(-1,1)")
        self.assertTrue(
            np.array_equal(
                exp2expstr([0, 1, 2, 3], "-i"),
                np.array(["(-i,0)", "(-i,1)", "(-i,2)", "(-i,3)"], dtype="<U6"),
            )
        )

    def test_exp2exp(self):
        """Tests exp2exp function."""
        self.assertTrue(
            np.array_equal(
                exp2exp(np.array([1, 2, 1, 3, 2]), "i", "is"),
                np.array([[1, 0], [0, 1], [1, 0], [1, 1], [0, 1]]),
            )
        )
        self.assertTrue(
            np.array_equal(
                exp2exp(np.array([(0, 1), (1, 1), (1, 0), (1, 1)]), "is", "-i"),
                np.array([2, 1, 3, 1]),
            )
        )

    def test_cpxstr2expstr(self):
        """Tests cpxstr2expstr function."""
        self.assertEqual(cpxstr2expstr("1", "i"), "(i,0)")
        self.assertEqual(cpxstr2expstr("-i", "is"), "(i,1)(-1,1)")
        self.assertTrue(
            np.array_equal(
                cpxstr2expstr(["1", "-1", "i", "-i"], "is"),
                np.array(
                    ["(i,0)(-1,0)", "(i,0)(-1,1)", "(i,1)(-1,0)", "(i,1)(-1,1)"],
                    dtype="<U11",
                ),
            )
        )

    def test_expstr2cpxstr(self):
        """Tests expstr2cpxstr function."""
        self.assertEqual(expstr2cpxstr("(-i,1)(-1,1)", "-is"), "i")
        self.assertEqual(expstr2cpxstr("(-i,2)", "-i", ones=True), "-1")

    def test_cpxstr2exp(self):
        """Tests cpxstr2exp function."""
        self.assertEqual(cpxstr2exp("-i", "i"), 3)
        self.assertTrue(np.array_equal(cpxstr2exp("-i", "-is"), np.array([1, 0])))

    def test_exp2cpxstr(self):
        """Tests exp2cpxstr function."""
        self.assertEqual(expstr2cpx("(i,1)", "i"), 1j)
        self.assertEqual(expstr2cpx("(i,1)(-1,1)", "is"), (-0 - 1j))
        self.assertTrue(
            np.array_equal(
                expstr2cpx(["(i,0)", "(i,1)", "(i,1)", "(i,3)"], "i"),
                np.array([1.0 + 0.0j, 0.0 + 1.0j, 0.0 + 1.0j, -0.0 - 1.0j]),
            )
        )

    def test_cpx2expstr(self):
        """Tests cpx2expstr function."""
        self.assertEqual(cpx2expstr(1j, "-i"), "(-i,3)")
        self.assertTrue(
            np.array_equal(
                cpx2expstr([1j, 1, -1j, -1], "-is"),
                np.array(
                    ["(-i,1)(-1,1)", "(-i,0)(-1,0)", "(-i,1)(-1,0)", "(-i,0)(-1,1)"],
                    dtype="<U12",
                ),
            )
        )

    def test_str2exp(self):
        """Tests str2exp function."""
        self.assertEqual(str2exp("i", "i"), 1)
        self.assertEqual(str2exp("(-i,3)", "i"), 1)
        self.assertEqual(str2exp("", "i"), 0)
        self.assertTrue(
            np.array_equal(
                str2exp(["i", "-i", "+1", "1"], "-is"),
                np.array([[1, 1], [1, 0], [0, 0], [0, 0]]),
            )
        )

    def test_split_pauli(self):
        """Tests split_pauli function."""
        self.assertEqual(("(-i,3)", "XIIXYZZX"), split_pauli("(-i,3)XIIXYZZX"))
        self.assertEqual(("(-i,1)(-1,0)", "XXIXZ"), split_pauli("(-i,1)(-1,0)XXIXZ"))
        self.assertEqual(("i", "X1Z3Y7"), split_pauli("iX1Z3Y7"))

    def test_encode_of_phase_str(self):
        """Tests encode_of_phase_str function."""
        self.assertEqual(encode_of_phase_str("(i,0)(-1,1)"), "is")
        self.assertEqual(encode_of_phase_str("+1j"), "complex")
        self.assertEqual(encode_of_phase_str("(i,0)"), "i")

    def test_encode_of_tensor_str(self):
        """Tests encode_of_tensor_str function."""
        rep, syn = encode_of_tensor_str("X1Y2Z5Y4", encoded=False)
        self.assertEqual(rep, ["XZY", "YZX"])
        self.assertEqual(syn, "index")

        rep, syn = encode_of_tensor_str("XYXZZIX", encoded=False)
        self.assertEqual(rep, ["XZY", "YZX"])

        rep, syn = encode_of_tensor_str("(X1)(Z2)(Z3X3)(Z10)", encoded=False)
        self.assertEqual(rep, ["ZX"])
        self.assertEqual(syn, "index")

    def test_str2symplectic(self):
        """Tests str2symplectic function."""
        matrix, phase_exp = str2symplectic(
            "iXXIZY", qubit_order="left-to-right", output_encoding="-isXZ"
        )
        self.assertTrue(np.array_equal(np.array([0, 1]), phase_exp))
        self.assertTrue(np.array_equal(np.array([[1, 1, 0, 0, 1, 0, 0, 0, 1, 1]]), matrix))

        matrix, phase_exp = str2symplectic(
            ["iYII", "-iX0Z2", "X1Z2"],
            qubit_order="left-to-right",
            output_encoding="-iXZY",
        )
        self.assertTrue(np.array_equal(np.array([3, 1, 0]), phase_exp))
        self.assertTrue(
            np.array_equal(
                np.array([[1, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 1]]),
                matrix,
            )
        )

        matrix, phase_exp = str2symplectic(
            "iX1X3Y4Z9",
            qubit_order="left-to-right",
            output_encoding="-isXZ",
            index_start=1,
        )
        self.assertTrue(np.array_equal(np.array([0, 1]), phase_exp))
        self.assertTrue(
            np.array_equal(
                np.array([[1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]]),
                matrix,
            )
        )

        matrix, phase_exp = str2symplectic(
            "iX_{1}",
            qubit_order="left-to-right",
            output_encoding="-isXZ",
            index_start=0,
        )
        self.assertTrue(np.array_equal(np.array([1, 1]), phase_exp))
        self.assertTrue(
            np.array_equal(
                np.array([[0, 1, 0, 0]]),
                matrix,
            )
        )

        matrix, phase_exp = str2symplectic(
            "iX_{1}X_{3}Y_{4}Z_{9}",
            qubit_order="left-to-right",
            output_encoding="-isXZ",
            index_start=1,
        )

        self.assertTrue(np.array_equal(np.array([0, 1]), phase_exp))
        self.assertTrue(
            np.array_equal(
                np.array([[1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]]),
                matrix,
            )
        )

        matrix, phase_exp = str2symplectic("iXXXIZYZ")
        self.assertEqual(phase_exp, 0)
        self.assertTrue(
            np.array_equal(np.array([[0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]]), matrix)
        )

    def test_symplectic2str(self):
        """Tests symplectic2str function."""
        matrix = np.array(
            [
                [0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1],
            ],
            dtype=np.bool_,
        )
        phase_exp = np.array([3, 1], dtype=np.int8)

        self.assertTrue(
            np.array_equal(
                symplectic2str(matrix, phase_exp),
                np.array(["-iY9Y4X3X1", "iZ10X8Y7Y6X3"], dtype="<U12"),
            )
        )
        self.assertTrue(
            np.array_equal(
                symplectic2str(matrix, phase_exp, qubit_order="left-to-right"),
                np.array(["-iX1X3Y4Y9", "iX3Y6Y7X8Z10"], dtype="<U12"),
            )
        )

        self.assertTrue(
            np.array_equal(
                symplectic2str(matrix, phase_exp, qubit_order="left-to-right", syntax=LATEX_SYNTAX),
                np.array(["-iX_{1}X_{3}Y_{4}Y_{9}", "iX_{3}Y_{6}Y_{7}X_{8}Z_{10}"], dtype="<U30"),
            )
        )

        self.assertTrue(
            np.array_equal(
                symplectic2str(matrix),
                np.array(["-Y9Y4X3X1", "-Z10X8Y7Y6X3"], dtype="<U12"),
            )
        )
        self.assertEqual(symplectic2str(matrix[0]), "-Y9Y4X3X1")
        self.assertEqual(symplectic2str(matrix[0], no_phase=True), "Y9Y4X3X1")
        self.assertEqual(
            symplectic2str(matrix[0], no_phase=True, syntax=LATEX_SYNTAX), "Y_{9}Y_{4}X_{3}X_{1}"
        )
        self.assertTrue(
            np.array_equal(
                symplectic2str(matrix),
                np.array(["-Y9Y4X3X1", "-Z10X8Y7Y6X3"], dtype="<U12"),
            )
        )

        self.assertTrue(
            np.array_equal(
                symplectic2str(
                    matrix,
                    phase_exp,
                    input_encoding="iXZ",
                    output_phase_encoding=None,
                    output_tensor_encoding="XZY",
                    syntax=INDEX_SYNTAX,
                ),
                np.array(["iY9Y4X3X1", "-iZ10X8Y7Y6X3"], dtype="<U13"),
            )
        )

        self.assertTrue(
            np.array_equal(
                symplectic2str(
                    matrix,
                    phase_exp,
                    input_encoding="iXZ",
                    output_phase_encoding=None,
                    output_tensor_encoding="XZY",
                    syntax=LATEX_SYNTAX,
                ),
                np.array(["iY_{9}Y_{4}X_{3}X_{1}", "-iZ_{10}X_{8}Y_{7}Y_{6}X_{3}"], dtype="<U40"),
            )
        )

        self.assertTrue(
            np.array_equal(
                symplectic2str(
                    matrix,
                    phase_exp,
                    input_encoding="iXZ",
                    output_phase_encoding=None,
                    output_tensor_encoding="XZ",
                    syntax=INDEX_SYNTAX,
                ),
                np.array(
                    ["-i(X9Z9)(X4Z4)(X3)(X1)", "i(Z10)(X8)(X7Z7)(X6Z6)(X3)"],
                    dtype="<U26",
                ),
            )
        )

        self.assertTrue(
            np.array_equal(
                symplectic2str(
                    matrix,
                    phase_exp,
                    input_encoding="iXZ",
                    output_phase_encoding=None,
                    output_tensor_encoding="XZ",
                    syntax=LATEX_SYNTAX,
                ),
                np.array(
                    [
                        "-i(X_{9}Z_{9})(X_{4}Z_{4})(X_{3})(X_{1})",
                        "i(Z_{10})(X_{8})(X_{7}Z_{7})(X_{6}Z_{6})(X_{3})",
                    ],
                    dtype="<U50",
                ),
            )
        )

        self.assertTrue(
            np.array_equal(
                symplectic2str(
                    matrix,
                    phase_exp,
                    input_encoding="iXZ",
                    output_phase_encoding="-is",
                    output_tensor_encoding="XZ",
                    syntax=PRODUCT_SYNTAX,
                ),
                np.array(
                    [
                        "(-i,1)(-1,0)(I)(XZ)(I)(I)(I)(I)(XZ)(X)(I)(X)(I)",
                        "(-i,1)(-1,1)(Z)(I)(X)(XZ)(XZ)(I)(I)(X)(I)(I)(I)",
                    ],
                    dtype="<U47",
                ),
            )
        )
        self.assertTrue(
            np.array_equal(
                symplectic2str(
                    matrix,
                    phase_exp,
                    input_encoding="iXZ",
                    output_phase_encoding="-is",
                    output_tensor_encoding="XZY",
                    syntax=INDEX_SYNTAX,
                    index_start=2,
                ),
                np.array(["(-i,1)(-1,1)Y11Y6X5X3", "(-i,1)(-1,0)Z12X10Y9Y8X5"], dtype="<U24"),
            )
        )
        self.assertTrue(
            np.array_equal(
                symplectic2str(
                    matrix,
                    phase_exp,
                    input_encoding="iXZ",
                    output_phase_encoding="-is",
                    qubit_order="left-to-right",
                    output_tensor_encoding="XZY",
                    syntax=INDEX_SYNTAX,
                    index_start=1,
                    index_str="_",
                ),
                np.array(
                    ["(-i,1)(-1,1)X_2X_4Y_5Y_10", "(-i,1)(-1,0)X_4Y_7Y_8X_9Z_11"],
                    dtype="<U28",
                ),
            )
        )

    def test_str2str(self):
        """Tests str2str function"""
        self.assertEqual(
            str2str("iX1X3Y4Z9", LATEX_SYNTAX, qubit_order_output="left-to-right"),
            "iX_{1}X_{3}Y_{4}Z_{9}",
        )
        self.assertEqual(
            str2str("iX_{1}X_{3}Y_{4}Z_{9}", INDEX_SYNTAX, qubit_order_output="left-to-right"),
            "iX1X3Y4Z9",
        )
        self.assertEqual(str2str("IIX", INDEX_SYNTAX), "X0")
        self.assertEqual(str2str("Y_{3}Y_{4}", INDEX_SYNTAX), "Y4Y3")
        self.assertEqual(
            str2str("Y_{3}Y_{4}", PRODUCT_SYNTAX, qubit_order_output="left-to-right"), "IIIYY"
        )
        self.assertEqual(
            str2str(
                "Y_{3}Y_{4}",
                PRODUCT_SYNTAX,
                index_start_input=3,
                index_start_output=0,
                qubit_order_output="left-to-right",
            ),
            "YY",
        )
        self.assertEqual(
            str2str(
                "iIYY",
                LATEX_SYNTAX,
                index_start_input=0,
                index_start_output=0,
                qubit_order_input="left-to-right",
                qubit_order_output="left-to-right",
            ),
            "iY_{1}Y_{2}",
        )
        self.assertEqual(
            str2str(
                "iIYY",
                LATEX_SYNTAX,
                index_start_input=0,
                index_start_output=0,
                qubit_order_input="right-to-left",
                qubit_order_output="left-to-right",
            ),
            "iY_{0}Y_{1}",
        )
        self.assertEqual(
            str2str(
                "Y3Y4",
                PRODUCT_SYNTAX,
                index_start_input=0,
                index_start_output=0,
                qubit_order_output="left-to-right",
            ),
            "IIIYY",
        )

        self.assertEqual(
            str2str(
                "iX_{1}",
                phase_encoding_output_string="-is",
                tensor_encoding_output_string="XZ",
                syntax_output=INDEX_SYNTAX,
            ),
            "(-i,1)(-1,1)(X1)",
        )

        self.assertEqual(
            str2str(
                "iX_{1}Y_{2}",
                phase_encoding_output_string="is",
                tensor_encoding_output_string="XZY",
                syntax_output=LATEX_SYNTAX,
            ),
            "(i,1)(-1,1)Y_{2}X_{1}",
        )

        self.assertEqual(
            str2str(
                "ZY",
                phase_encoding_output_string="-i",
                tensor_encoding_output_string="ZX",
                syntax_output=LATEX_SYNTAX,
            ),
            "(-i,1)(Z_{1})(Z_{0}X_{0})",
        )
