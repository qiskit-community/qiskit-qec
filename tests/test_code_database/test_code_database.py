"""Indexer object test."""
import unittest
from qiskit_qec.qec_codebase.qec_codebase_interactor import QECCodeBase


class TestCodeBase(unittest.TestCase):
    """Test Indexer class."""

    def test_retrieve_small_query(self):
        qec_db = QECCodeBase()

        small_code_sets = [
            ((1, 0), set(["b8ad169a-2d48-454e-80bc-957f7d45b807"])),
            (
                (6, 0),
                set(
                    [
                        "dd40352d-7e6b-465c-bcec-4a28ef8b46b3",
                        "debd6be2-c623-4bfd-80a3-f77a4a337456",
                        "439ed374-3d68-43ec-8fc4-d629952b4648",
                        "438114b1-44d0-4105-90b8-e38458cc9654",
                        "67eca164-5071-4cc6-abb3-4470f4a1afae",
                        "391ab6bf-6bac-4881-9320-30428ea225eb",
                        "16527a8a-d824-4b7d-a933-4175b4d87919",
                        "8f04b1c2-03ab-44a2-9f43-bafeca57de13",
                        "c4a9ef5f-331c-4f67-afe7-7a23dfa26b7f",
                        "64f371c7-0fe9-48b6-ba61-fcb435668fe5",
                        "707bd7b0-c034-4633-b43d-9db20fa9aafa",
                        "6ee5a0ef-5b61-48c6-94e5-fd4147f317c7",
                        "2fa2df4f-9b1d-41dc-9229-dc500866b30d",
                        "f85849f0-52fa-4ebe-88a7-da7bd72c945a",
                        "bb3c6e57-0913-4fa6-bac8-baa67d7fb8f0",
                        "1a4fe6ef-23c2-4f19-a344-ed5d940413ea",
                        "f1877244-e620-4906-8b60-eb01c72a8e86",
                        "8d9abcf9-82e2-4f36-aae7-842426eab115",
                        "3bba9f88-527a-4c61-aac3-bcfa4df64466",
                        "63c53914-031d-4040-a6a6-182dac147074",
                        "6a8d7438-1683-498d-a17e-2022a053d3a6",
                        "e8d45e4c-4a85-4ff9-b441-5a4f0f016545",
                        "02568544-a737-4084-b9ed-2942deb92b9b",
                        "28a10c13-6ee2-4aff-860b-87f58d633734",
                        "7f752ebc-cedc-4af4-92a7-05c9992a3764",
                        "70838bd3-dda1-498d-9138-cd5626e6eef9",
                    ]
                ),
            ),
        ]

        for code_test in small_code_sets:
            search_input = code_test[0]
            soln = code_test[1]

            queried_codes_json = qec_db._load_code_jsons_from_n_k(*search_input)

            assert len(soln) == len(queried_codes_json)
            for q_code_uuid in queried_codes_json.keys():
                assert q_code_uuid in soln, f"quiid: {q_code_uuid} not in soln {soln}"

    def test_retrieve_small_codes(self):
        qec_db = QECCodeBase()

        small_code_sets = [
            ((1, 0), set(["b8ad169a-2d48-454e-80bc-957f7d45b807"])),
            (
                (6, 0),
                set(
                    [
                        "dd40352d-7e6b-465c-bcec-4a28ef8b46b3",
                        "debd6be2-c623-4bfd-80a3-f77a4a337456",
                        "439ed374-3d68-43ec-8fc4-d629952b4648",
                        "438114b1-44d0-4105-90b8-e38458cc9654",
                        "67eca164-5071-4cc6-abb3-4470f4a1afae",
                        "391ab6bf-6bac-4881-9320-30428ea225eb",
                        "16527a8a-d824-4b7d-a933-4175b4d87919",
                        "8f04b1c2-03ab-44a2-9f43-bafeca57de13",
                        "c4a9ef5f-331c-4f67-afe7-7a23dfa26b7f",
                        "64f371c7-0fe9-48b6-ba61-fcb435668fe5",
                        "707bd7b0-c034-4633-b43d-9db20fa9aafa",
                        "6ee5a0ef-5b61-48c6-94e5-fd4147f317c7",
                        "2fa2df4f-9b1d-41dc-9229-dc500866b30d",
                        "f85849f0-52fa-4ebe-88a7-da7bd72c945a",
                        "bb3c6e57-0913-4fa6-bac8-baa67d7fb8f0",
                        "1a4fe6ef-23c2-4f19-a344-ed5d940413ea",
                        "f1877244-e620-4906-8b60-eb01c72a8e86",
                        "8d9abcf9-82e2-4f36-aae7-842426eab115",
                        "3bba9f88-527a-4c61-aac3-bcfa4df64466",
                        "63c53914-031d-4040-a6a6-182dac147074",
                        "6a8d7438-1683-498d-a17e-2022a053d3a6",
                        "e8d45e4c-4a85-4ff9-b441-5a4f0f016545",
                        "02568544-a737-4084-b9ed-2942deb92b9b",
                        "28a10c13-6ee2-4aff-860b-87f58d633734",
                        "7f752ebc-cedc-4af4-92a7-05c9992a3764",
                        "70838bd3-dda1-498d-9138-cd5626e6eef9",
                    ]
                ),
            ),
        ]

        for code_test in small_code_sets:
            search_input = code_test[0]
            soln = code_test[1]

            queried_codes = qec_db.get_subsystem_code(*search_input)
            assert len(soln) == len(queried_codes)
            for q_code in queried_codes:
                assert (
                    q_code.parameters["uuid"] in soln
                ), f"quiid: {q_code.parameters['uuid']} not in soln {soln}"

    def test_unrecorded_code(self):
        qec_db = QECCodeBase()
        invalid_code = (10, 5)
        with self.assertRaises(FileNotFoundError):
            qec_db.get_subsystem_code(*invalid_code)

    def test_additional_params(self):
        qec_db = QECCodeBase()

        addtl_param_tests = [
            (((1, 0), {"gottesman_form": "y"}, {}), []),
            (
                ((3, 2), {"css_logicals": "ixi"}, {}),
                ["8adffac4-a07a-4f07-ad14-a936d73b0d8e", "d9fbf163-c283-49d2-a326-caa304c3de21"],
            ),
            (
                (
                    (2, 1),
                    {"definitely_not_a_real_param": 7},
                    {"definitely_not_a_real_param"},
                ),
                ["9f95c997-4048-4001-b2fd-aa919be5866f", "6b47540a-46a3-4e6f-bc97-886fe08caf44"],
            ),
            (
                (
                    (2, 1),
                    {"definitely_not_a_real_param": 7},
                    {},
                ),
                [],
            ),
        ]

        for code_test in addtl_param_tests:
            search_input = code_test[0]
            n_k_input = search_input[0]
            addtl_params = search_input[1]
            addtl_params_that_neednt_exist = search_input[2]
            soln = code_test[1]
            queried_codes = qec_db.get_subsystem_code(
                *n_k_input, addtl_params_that_neednt_exist, **addtl_params
            )
            assert len(soln) == len(queried_codes)
            for q_code in queried_codes:
                assert (
                    q_code.parameters["uuid"] in soln
                ), f" q_code_uuid: {q_code.parameters[ 'uuid' ]} not in soln {soln}"
                for key, value in addtl_params.items():
                    if key not in addtl_params_that_neednt_exist:
                        assert (
                            value in q_code.parameters[key] or q_code.parameters[key] == value
                        ), f"q_code.parameters[{key}]: {q_code.parameters[key]} fails for value: {value}"


if __name__ == "__main__":
    unittest.main()
