"""Test python error propagator selector."""
import unittest
import copy
import os
from qiskit_qec.codes import SubSystemCode
from qiskit_qec.qec_codebase.qec_codebase_interactor import QECCodeBase
from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.structures.gauge import GaugeGroup
from qiskit_qec.exceptions import QiskitQECError


class TestQECCodebaseInteractor(unittest.TestCase):
    """Test QECCodebase interactor"""

    SOLN = {
        (2, 0): {
            "b5dfa20f-8d48-4fd1-aa6f-0556f3021558": {
                "stabilizer": ["xi", "ix"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 8,
                "n": 2,
                "k": 0,
                "is_degenerate": 0,
                "weight_enumerator": [1, 2, 1],
                "d": 1,
                "uuid": "b5dfa20f-8d48-4fd1-aa6f-0556f3021558",
                "is_css": 1,
                "is_decomposable": 1,
                "is_gf4linear": 0,
            },
            "8ce0c4f9-135f-41b1-8609-ffcd9bc49362": {
                "stabilizer": ["xx", "zz"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 12,
                "n": 2,
                "k": 0,
                "is_degenerate": 0,
                "is_triorthogonal": 0,
                "weight_enumerator": [1, 0, 3],
                "d": 2,
                "uuid": "8ce0c4f9-135f-41b1-8609-ffcd9bc49362",
                "is_css": 1,
                "is_decomposable": 0,
                "is_gf4linear": 1,
            },
        },
        (2, 1): {
            "9f95c997-4048-4001-b2fd-aa919be5866f": {
                "stabilizer": ["xx"],
                "logical_ops": ["ix", "zz"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 8,
                "n": 2,
                "k": 1,
                "is_degenerate": 0,
                "weight_enumerator": [1, 0, 1],
                "d": 1,
                "uuid": "9f95c997-4048-4001-b2fd-aa919be5866f",
                "is_css": 1,
                "is_decomposable": 0,
                "is_gf4linear": 0,
            },
            "6b47540a-46a3-4e6f-bc97-886fe08caf44": {
                "stabilizer": ["xi"],
                "logical_ops": ["ix", "iz"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 12,
                "n": 2,
                "k": 1,
                "is_degenerate": 0,
                "weight_enumerator": [1, 1, 0],
                "d": 1,
                "uuid": "6b47540a-46a3-4e6f-bc97-886fe08caf44",
                "is_css": 1,
                "is_decomposable": 0,
                "is_gf4linear": 0,
            },
        },
        (3, 0): {
            "a02d2b58-79c0-4c9e-8122-90cd57a57e6a": {
                "stabilizer": ["xii", "ixi", "iix"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 48,
                "n": 3,
                "k": 0,
                "is_degenerate": 0,
                "weight_enumerator": [1, 3, 3, 1],
                "d": 1,
                "uuid": "a02d2b58-79c0-4c9e-8122-90cd57a57e6a",
                "is_css": 1,
                "is_decomposable": 1,
                "is_gf4linear": 0,
            },
            "a814c308-296f-4ec3-91d6-9509408e56da": {
                "stabilizer": ["xix", "ixx", "zzz"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 24,
                "n": 3,
                "k": 0,
                "is_degenerate": 0,
                "is_triorthogonal": 0,
                "weight_enumerator": [1, 0, 3, 4],
                "d": 2,
                "uuid": "a814c308-296f-4ec3-91d6-9509408e56da",
                "is_css": 1,
                "is_decomposable": 0,
                "is_gf4linear": 0,
            },
            "12d90740-cd73-44c1-8f33-79c6cad4596a": {
                "stabilizer": ["ixi", "ziz", "yiy"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 24,
                "n": 3,
                "k": 0,
                "is_degenerate": 0,
                "weight_enumerator": [1, 1, 3, 3],
                "d": 1,
                "uuid": "12d90740-cd73-44c1-8f33-79c6cad4596a",
                "is_css": 1,
                "is_decomposable": 1,
                "is_gf4linear": 0,
            },
        },
        (3, 1): {
            "a7753bb2-01f8-461b-aedd-24d42a7322da": {
                "stabilizer": ["xii", "ixi"],
                "logical_ops": ["iix", "iiz"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 48,
                "n": 3,
                "k": 1,
                "is_degenerate": 0,
                "weight_enumerator": [1, 2, 1, 0],
                "d": 1,
                "uuid": "a7753bb2-01f8-461b-aedd-24d42a7322da",
                "is_css": 1,
                "is_decomposable": 1,
                "is_gf4linear": 0,
            },
            "d340c69c-0310-4ce9-af44-1fd83c3ce67b": {
                "stabilizer": ["ixi", "xxx"],
                "logical_ops": ["iix", "ziz"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 16,
                "n": 3,
                "k": 1,
                "is_degenerate": 0,
                "weight_enumerator": [1, 1, 1, 1],
                "d": 1,
                "uuid": "d340c69c-0310-4ce9-af44-1fd83c3ce67b",
                "is_css": 1,
                "is_decomposable": 1,
                "is_gf4linear": 0,
            },
            "103f7c5b-8601-4175-844d-695eee247b71": {
                "stabilizer": ["xix", "ixx"],
                "logical_ops": ["iix", "zzz"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 48,
                "n": 3,
                "k": 1,
                "is_degenerate": 0,
                "weight_enumerator": [1, 0, 3, 0],
                "d": 1,
                "uuid": "103f7c5b-8601-4175-844d-695eee247b71",
                "is_css": 1,
                "is_decomposable": 0,
                "is_gf4linear": 0,
            },
            "d6ade957-fdac-45d5-814c-d386c3f4604a": {
                "stabilizer": ["xxi", "zzi"],
                "logical_ops": ["iix", "iiz"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 72,
                "n": 3,
                "k": 1,
                "is_degenerate": 0,
                "weight_enumerator": [1, 0, 3, 0],
                "d": 1,
                "uuid": "d6ade957-fdac-45d5-814c-d386c3f4604a",
                "is_css": 1,
                "is_decomposable": 0,
                "is_gf4linear": 1,
            },
            "0a80b836-3c62-47ec-9efa-a417b8b01b12": {
                "stabilizer": ["xix", "zzz"],
                "logical_ops": ["ixx", "ziz"],
                "is_subsystem": 1,
                "belongs_to_standard_codebase": 1,
                "aut_group_size": 8,
                "n": 3,
                "k": 1,
                "is_degenerate": 0,
                "weight_enumerator": [1, 0, 1, 2],
                "d": 1,
                "uuid": "0a80b836-3c62-47ec-9efa-a417b8b01b12",
                "is_css": 1,
                "is_decomposable": 0,
                "is_gf4linear": 0,
            },
        },
    }

    @classmethod
    def setUpClass(cls) -> None:
        """Set up interactor for class"""
        codebase_name = "test_codebase"
        codebase_dir_path = os.path.join(
            *([os.sep] + os.path.abspath(__file__).split(os.sep)[:-1] + [codebase_name])
        )
        cls.db = QECCodeBase(codebase_dir_path)

    def tearDown(self) -> None:
        """Delete everything in the playground"""
        self.db.delete_playground_codebase(True)

    def test_get_basic_code(self):
        """Test getting a basic code"""
        basic_code_soln = self.SOLN[(2, 1)]
        basic_codes = self.db.get_subsystem_code(2, 1)
        self.assertEqual(len(basic_codes), 2)
        for code in basic_codes:
            self.assertEqual(basic_code_soln[code.parameters[self.db.QEC_UUID]], code.parameters)

    def test_multi_codes(self):
        """Test getting many codes with different n,k"""
        multi_codes = self.db.get_subsystem_code([3, 2], [0, 1])
        multi_code_soln = {}
        for i in [2, 3]:
            for j in [0, 1]:
                multi_code_soln = {**multi_code_soln, **self.SOLN[i, j]}
        for code in multi_codes:
            self.assertEqual(multi_code_soln[code.parameters[self.db.QEC_UUID]], code.parameters)

    def test_empty_query(self):
        """Test getting no codes"""
        self.assertEqual(len(self.db.get_subsystem_code(5, 4)), 0)

    def test_query_with_singular_param(self):
        """Test w/ param"""
        test_codes = [[(3, 2), {"logical_ops": "ixi"}, 2], [(2, 1), {"aut_group_size": 12}, 1]]
        for test in test_codes:
            returned_codes = self.db.get_subsystem_code(*test[0], **test[1])
            self.assertEqual(len(returned_codes), test[2])
            for code in returned_codes:
                code_val = code.parameters[list(test[1].keys())[0]]
                test_val = list(test[1].values())[0]
                self.assertTrue(code_val == test_val or test_val in code_val)

    def test_query_multiple_param(self):
        """Query with multiple parameters"""
        qargs = [
            (3, 1),
            {self.db.WEIGHT_ENUMERATOR: 3, self.db.IS_GF4LINEAR: 1},
            1,
            "d6ade957-fdac-45d5-814c-d386c3f4604a",
        ]
        retrieved_codes = self.db.get_subsystem_code(*qargs[0], **qargs[1])
        self.assertEqual(len(retrieved_codes), qargs[2])
        self.assertEqual(retrieved_codes[0].parameters[self.db.QEC_UUID], qargs[3])

        # none
        self.assertEqual(
            len(
                self.db.get_subsystem_code(
                    *(3, 1), **{self.db.AUT_GROUP_SIZE: 8, self.db.IS_DEGENERATE: 1}
                )
            ),
            0,
        )

    def test_create_query(self):
        """Test making new code"""
        qargs = (2, 1)
        prev_codes = self.db.get_subsystem_code(*qargs)
        original_len = len(prev_codes)
        one_off_id = list(self.SOLN[qargs].keys())[0]
        code_info = copy.deepcopy(self.SOLN[qargs][one_off_id])
        my_code = SubSystemCode(GaugeGroup(PauliList(code_info[self.db.STABILIZER][0].upper())))
        my_code.parameters = code_info
        my_code.parameters["diffle_wumpus"] = 123
        with self.assertRaises(QiskitQECError) as cur_context:
            self.db.store_new_subsystem_code(my_code)
            self.assertTrue(
                "is a new field and you have set allow_new_fields to False" in cur_context.exception
            )

        self.db.store_new_subsystem_code(my_code, allow_new_fields=True)
        new_codes = self.db.get_subsystem_code(*qargs)
        self.assertEqual(len(new_codes), original_len + 1)

        uuids = set()
        for code in new_codes:
            cur_uuid = code.parameters[self.db.QEC_UUID]
            self.assertTrue(cur_uuid not in uuids)
            uuids.add(cur_uuid)

    def test_create_query_nonunique_code(self):
        """Test that db yells at you if you try to save a code that already exists without force"""
        qargs = (2, 1)
        one_off_id = list(self.SOLN[qargs].keys())[0]
        code_info = copy.deepcopy(self.SOLN[qargs][one_off_id])
        code_info.pop(self.db.QEC_UUID)
        my_code = SubSystemCode(GaugeGroup(PauliList(code_info[self.db.STABILIZER][0].upper())))
        my_code.parameters = code_info
        with self.assertRaises(QiskitQECError) as cur_context:
            self.db.store_new_subsystem_code(my_code)
            self.assertTrue(
                "is not completely unique from other code in the database:" in cur_context.exception
            )

    def test_create_query_force_nonunique_code(self):
        """Test making new code that is not unique but saving it anyway"""
        qargs = (2, 1)
        one_off_id = list(self.SOLN[qargs].keys())[0]
        code_info = copy.deepcopy(self.SOLN[qargs][one_off_id])
        my_code = SubSystemCode(GaugeGroup(PauliList(code_info[self.db.STABILIZER][0].upper())))
        my_code.parameters = code_info
        code_info.pop(self.db.QEC_UUID)
        code_info.pop(self.db.BELONGS_TO_STANDARD_CODEBASE)
        self.db.store_new_subsystem_code(my_code, force=True)
        retrieved_codes = self.db.get_subsystem_code(*qargs, **code_info)
        self.assertEqual(len(retrieved_codes), 2)
        self.assertNotEqual(
            retrieved_codes[0].parameters[self.db.QEC_UUID],
            retrieved_codes[1].parameters[self.db.QEC_UUID],
        )

    def test_cache(self):
        """Test Cache"""
        my_code = SubSystemCode(GaugeGroup(PauliList("XX")))
        my_code.parameters = {self.db.IS_GF4LINEAR: 1, self.db.N: 4, self.db.K: 2}

        for _ in range(300):
            self.db.store_new_subsystem_code(my_code, flush_cache=False, force=True)
        self.db.flush_cache()
        self.assertEqual(len(self.db.playground_cache), 0)
        retrieved_codes = self.db.get_subsystem_code(4, 2)
        self.assertEqual(len(retrieved_codes), 300)


if __name__ == "__main__":
    unittest.main()
