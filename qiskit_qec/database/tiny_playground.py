import json
import os
import uuid
from timeit import default_timer as timer
from tinydb import Query, TinyDB
# from qiskit_qec.codes.subsystemcodes import SubSystemCode
# from qiskit_qec.operators.pauli_list import PauliList
# from qiskit_qec.structures.gauge import GaugeGroup


# # Custom test with parameters:
# class TempDB:
#     """Not for extended use"""

#     count = 0

#     def __init__(self, database_file):
#         if self.count > 0:
#             print("WARNING. You're opening > 1 instance of TempDB.")
#         self._db = TinyDB(database_file)

#     def get_codes_with_attributes(self, **frag_dict):
#         return self._db.search(Query().fragment(frag_dict))

#     def get_codes_with_dict(self, frag_dict: dict):
#         return self._db.search(Query().fragment(frag_dict))

#     def get_a_code_with_attributes(self, **frag_dict):
#         return self._codify(self._db.search(Query().fragment(frag_dict))[0])

#     def get_a_code_with_dict(self, frag_dict: dict):
#         return self._codify(self._db.search(Query().fragment(frag_dict))[0])

#     def _codify(self, code_info):  # TAKES code from list and inserts into new object
#         sscode = SubSystemCode(GaugeGroup(PauliList(code_info["gottesman_form"][0].upper())))
#         sscode.parameters = code_info

#         return sscode

#     #
#     # def in_range(val, m, n):
#     #     return m <= val <= n

#     # def insert_code(self):
#     #     raise Exception("PLEASE DO NOT INSERT NEW CODE ATM")


# if __name__ == "__main__":
#     db = TempDB("tiniest_db.json")  # put in db json
#     my_code = db.get_a_code_with_attributes(n=10, k=0, group_size_1=96.0)

#     print(my_code.gauge_group.generators)
#     print(my_code.parameters)




# class UUIDEncoder(json.JSONEncoder):
#     def default(self, obj):
#         if isinstance(obj, uuid.UUID):
#             return obj.hex
#         return json.JSONEncoder.default(self, obj)
DATADIR = "/Users/graceharperibm/correcting/stabilizer-codes/database"
test_fx = lambda f : f.endswith(".json")
NEWDB = "code_db.json"
used_uuids = set()

def unique_id_static():
    cur_uuid = uuid.uuid4().hex
    while cur_uuid in used_uuids:
        cur_uuid = uuid.uuid4().hex
    used_uuids.add(cur_uuid)
    return cur_uuid

def ensure_no_duplicate(db,code_info):
    print("...ensuring no duplicates...")
    return  len(db.search(Query().fragment(code_info))) == 0

def reading():
    db = TinyDB(NEWDB)
    
    json_files = [os.path.join(DATADIR, f) for f in os.listdir(DATADIR) if test_fx(f)]

    for jsonf in json_files:
        print(f"running... {jsonf}")
        with open(jsonf) as f:
            data = json.load(f)

        for k, code_info in data.items():
            # if ensure_no_duplicate(db, code_info):
            #     print("no duplicate found")
            cur_uuid = unique_id_static()
            code_info[ "UUID" ] = cur_uuid
            db.insert(code_info)
            # else:
            #     print(f"PANIC {code_info} has duplicated")
            #     return
        print(f"completed: {jsonf}")

# from typing import Callable


# def time_function(my_function:Callable):
    
#     def inner(*args, **kwargs):
#         start = timer()
#         func_ret = my_function(*args, **kwargs)
#         end = timer()
#         return (end - start, func_ret) 
    
#     return inner 


# @time_function
# def get_all_records(field:str, value):
#     return db.search(where(field) == value)
# TODO verify
# go through every json object and make sure it exists in the database once
def verify():
    # TODO Manually verify # of files
    # verify each json obj exists in our db
    pass
    
if __name__ == "__main__":
    start = timer()
    reading()
    end = timer()
    print(f"processing took {(end -start)/(60 * 60)} hours")
    
    # tenzero = "/Users/graceharperibm/correcting/stabilizer-codes/database/database_10_0.json"

    
    # # reading()
    
    # # # TODO run length off keys 
    # time_per_3990_field = 140.646
    
    # tot = 0
    # json_files = [os.path.join(DATADIR, f) for f in os.listdir(DATADIR) if test_fx(f)]
    # print(f"num files: {len(json_files)}")
    # for jsonf in json_files:
    #     print(f"jsonf is {jsonf}")
    #     with open(jsonf) as f:
    #         data = json.load(f)
    #     tot += len(data.values())
    # print(f" total docs: {tot}, for total time of: {(((tot/3990) * time_per_3990_field)/60)/60} hours")
    
