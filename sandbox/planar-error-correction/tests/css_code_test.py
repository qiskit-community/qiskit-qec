import unittest
from plerco.pec_python.qec_code.css_code import CSSCode


class TestCSSCode(unittest.TestCase):
    def test_steane_code(self):
        # Same X/Z stabilizers and S = G
        stabilizers = [[3, 4, 5, 6], [1, 2, 5, 6], [0, 2, 4, 6]]
        logical = [[0, 1, 2]]
        steane = CSSCode(
            stabilizers,
            stabilizers,
            7,
            1,
            3,
            stabilizers,
            stabilizers,
            logical,
            logical,
        )
        # All single-bit errors
        errors = [0, 1, 2, 4, 8, 16, 32, 64]
        synds = [0, 1, 2, 3, 4, 5, 6, 7]
        for i in range(8):
            s = format(errors[i], "07b")
            b = [int(x) for x in s]
            b.reverse()
            syn = steane.x_syndrome(b)
            self.assertEqual(syn, steane.z_syndrome(b))
            isyn = int("".join(map(str, syn)), 2)
            self.assertEqual(isyn, synds[i])
            # print(errors[i], s, b, syn, isyn)
            if isyn > 0:
                b[isyn - 1] ^= 1  # correct error
            self.assertEqual(steane.x_syndrome(b), [0, 0, 0])
            self.assertEqual(steane.logical_z_error(b), [0])
            self.assertEqual(steane.logical_x_error(b), [0])
            for i in logical[0]:
                b[i] ^= 1
            self.assertEqual(steane.logical_z_error(b), [1])
            self.assertEqual(steane.logical_x_error(b), [1])


if __name__ == "__main__":
    unittest.main()
