from tinydb import TinyDB, Query
from qiskit_qec.codes.subsystemcodes import SubSystemCode
from qiskit_qec.structures.gauge import GaugeGroup
from qiskit_qec.operators.pauli_list import PauliList

from typing import List

# Custom test with parameters:
class TempDB:
    """Not for extended use"""

    count = 0

    def __init__(self, database_file):
        if self.count > 0:
            print("WARNING. You're opening > 1 instance of TempDB.")
        self._db = TinyDB(database_file)

    def get_codes_with_attributes(self, **frag_dict):
        return self._db.search(Query().fragment(frag_dict))

    def get_codes_with_dict(self, frag_dict: dict):
        return self._db.search(Query().fragment(frag_dict))

    def get_a_code_with_attributes(self, **frag_dict):
        return self._codify(self._db.search(Query().fragment(frag_dict))[0])

    def get_a_code_with_dict(self, frag_dict: dict):
        return self._codify(self._db.search(Query().fragment(frag_dict))[0])

    def _codify(self, code_info):  # TAKES code from list and inserts into new object
        sscode = SubSystemCode(GaugeGroup(PauliList(code_info["gottesman_form"][0].upper())))
        sscode.parameters = code_info

        return sscode

    #
    # def in_range(val, m, n):
    #     return m <= val <= n

    # def insert_code(self):
    #     raise Exception("PLEASE DO NOT INSERT NEW CODE ATM")


if __name__ == "__main__":
    db = TempDB("tiniest_db.json")  # put in db json
    my_code = db.get_a_code_with_attributes(n=10, k=0, group_size_1=96.0)

    print(my_code.gauge_group.generators)
    print(my_code.parameters)


#
# class UUIDEncoder(json.JSONEncoder):
#     def default(self, obj):
#         if isinstance(obj, uuid.UUID):
#             return obj.hex
#         return json.JSONEncoder.default(self, obj)
# DATADIR = "/Users/graceharperibm/correcting/qiskit-qec/raw_data"
# test_fx = lambda f : f.endswith(".json")
# NEWDB = "code_db.json"
# db = TinyDB(NEWDB)
#
# used_uuids = set()
# def unique_id_static():
#     cur_uuid = uuid.uuid4().hex
#     while cur_uuid in used_uuids:
#
#         cur_uuid = uuid.uuid4().hex
#
#     used_uuids.add(cur_uuid)
#     return cur_uuid
#
# def ensure_no_duplicate(db,code_info):
#     return True
#     return  len(db.search(Query().fragment(code_info))) == 0
#
# def reading():
#     json_files = [os.path.join(DATADIR, f) for f in os.listdir(DATADIR) if test_fx(f)]
#     print(json_files)
#
#     for jsonf in json_files:
#         with open(jsonf) as f:
#             data = json.load(f)
#         for code_info in data.values():
#             if ensure_no_duplicate(db, code_info):
#                 cur_uuid = unique_id_static()
#                 code_info[ "UUID" ] = cur_uuid
#                 db.insert(code_info)
#             else:
#                 print(f"\n\n\n\n\n\n\n\n\n\n\nERROR PREEXISTED: {code_info}")
#         print(f"completed: {jsonf}")
#     for item in db:
#         print(item)
#     print(used_uuids)
#
