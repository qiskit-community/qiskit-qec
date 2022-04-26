from qiskit_qec.qec_codebase.qec_codebase_interactor import QECCodeBase
from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.structures.gauge import GaugeGroup
from qiskit_qec.codes import SubSystemCode
from qiskit_qec.exceptions import QiskitQECError

# download codebase from:
# then give QECCodeBase path to codebase: https://github.ibm.com/Grace-Harper/ibm-codebase-repo

if __name__ == "__main__":
    codebase_dir_path = "/Users/graceharperibm/correcting/QISKITQEC/ibm-codebase-repo/codebase"  # put your path here
    db = QECCodeBase(codebase_dir_path)

    # BASIC QUERY
    derulo_codes = db.get_subsystem_code(2, 1)  # all 2,1 codes
    print("derulo", len(derulo_codes))
    print(derulo_codes[0].gauge_group.generators)
    print(derulo_codes[0].parameters)

    # MULTI-QUERY
    swift_codes = db.get_subsystem_code([3, 2], [0, 1])  # all codes: (2,0), (2,1), (3,0), (3,1)
    print("swift", len(swift_codes))

    # QUERY sans results
    nope_codes = db.get_subsystem_code(3, 4)
    print("nope", len(nope_codes))

    # QUERY EXAMPLES w/ params
    the1975_additional_params = {"logical_ops": "ixi"}
    the1975_codes = db.get_subsystem_code(
        3, 2, **the1975_additional_params
    )  # all (3,2) codes that have css_logicals === "ixi"
    print("1975", len(the1975_codes))

    glass_animals_params = {"weight_enumerator": 1}
    glass_animals_codes = db.get_subsystem_code(
        2, 1, acceptable_missing_params={"weight_enumerator"}, **glass_animals_params
    )  # all codes (2,1) such that they either have "1" in their weight_enumerator list OR don't have weight_enumerator as a variable
    print("glass", len(glass_animals_codes))

    # Deriving your own code
    new_code = derulo_codes[0]
    new_code.parameters["my_new_param"] = 123
    print(new_code.parameters)

    db.store_new_subsystem_code(new_code, allow_new_fields=True)
    new_derulo_codes = db.get_subsystem_code(2, 1)  # all 2,1 codes
    print("new_derulo_codes", len(new_derulo_codes))

    # You can limit your query results to only those in the standard codebase (Andrew's codebase)
    #  by telling the database not to query the playground codebase

    only_standard_codes_derulo = db.get_subsystem_code(2, 1, allow_playground=False)
    print("only standard derulo codes", len(only_standard_codes_derulo))

    # If you're going to add a lot of codes at once, set cache=True and flush after:
    my_code = SubSystemCode(GaugeGroup(PauliList("XX")))
    my_code.parameters = {db.IS_GF4LINEAR: 1, db.N: 4, db.K: 2}

    for i in range(300):
        db.store_new_subsystem_code(my_code, flush_cache=False, force=True)

    db.flush_cache()  # actually create/write to the json file

    retrieved_codes = db.get_subsystem_code(4, 2)

    db.delete_playground_codebase(True)  # kill all the codes you've made
