from qiskit_qec.qec_codebase.qec_codebase_interactor import QECCodeBase
if __name__ == "__main__":
    db = QECCodeBase()
    derulo_codes = db.get_subsystem_code(2,1) # all 2,1 codes
    swift_codes = db.get_subsystem_code([3,2], [0,1]) # all codes: (2,0), (2,1), (3,0), (3,1)
    
    try:
        nope_codes = db.get_subsystem_code(3,4)
    except FileNotFoundError as e:
        print(f"n=3, k=4 is doesn't exist in our db: {e}")
    
    
    the1975_additional_params = {"css_logicals": "ixi"}
    the1975_codes = db.get_subsystem_code(3,2,**the1975_additional_params) # all (3,2) codes that have css_logicals === "ixi"
    
    glass_animals_params = {"weight_enumerator":1}
    glass_animals_codes = db.get_subsystem_code(2,1,addtl_params_that_neednt_exist={"weight_enumerator"} ,**glass_animals_params)  # all codes (2,1) such that they either have "1" in their weight_enumerator list OR don't have weight_enumerator as a variable

    
    print(derulo_codes)
    print(swift_codes)
    print(the1975_codes)
    print(glass_animals_codes)
    
    print(derulo_codes[0].gauge_group.generators)
    print(derulo_codes[0].parameters)