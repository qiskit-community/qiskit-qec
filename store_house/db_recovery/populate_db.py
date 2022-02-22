# import shutil
# import json
# import os, uuid
#
#
# # create dirs and copy original db into db dir
# ORIG_ABS = "<path to original database>"
# orig_file_name_format_d_k = (
#     ORIG_ABS + "database_{}_{}.json"
# )
# new_file_name = "/qec_codebase/n_length_{}_codes/n_length_{}_k_dim_{}_codes/"  # n,n,k
#
# for i in range(1, 11):
#     os.mkdir("n_length_" + str(i) + "_codes")
#     for j in range(i):
#         os.mkdir(
#             "n_length_" + str(i) + "_codes/" + "n_length_" + str(i) + "_k_dim_" + str(j) + "_codes"
#         )
#
# for n_length in range(1, 11):
#     for k_dim in range(n_length):
#         cur_file_name = orig_file_name_format_d_k.format(n_length, k_dim)
#         if os.path.exists(cur_file_name):
#             cur_new_file = new_file_name.format(n_length, n_length, k_dim)
#             shutil.copy(cur_file_name, cur_new_file)
#
#
#
#
# # convert codes into db format w/ str(uuid)
# ORIG_ABS = "qiskit-qec/qiskit_qec/qec_codebase"
# db_file_name = (
#     "/qec_codebase/n_length_{}_codes/n_length_{}_k_dim_{}_codes/database_{}_{}.json"  # n,n,k
# )
# github_file_name = (
#     "/qec_codebase/n_length_{}_codes/n_length_{}_k_dim_{}_codes/codes_{}_{}.json"  # n,n,k
# )
#
#
#
# uuid_dict = {}
#
# for n_length in range(1, 11):
#     for k_dim in range(n_length):
#         cur_file_name = ORIG_ABS + db_file_name.format(
#             n_length, n_length, k_dim, n_length, k_dim
#         )
#         new_file_name = ORIG_ABS + github_file_name.format(
#             n_length, n_length, k_dim, n_length, k_dim
#         )
#         if os.path.exists(cur_file_name):
#             with open(cur_file_name, "r") as f:
#                 data = json.load(f)
#
#             new_data = {}
#
#             for k, v in data.items():
#                 cur_uuid = str(uuid.uuid4())  # not sure what problems str-ing will cause
#                 if cur_uuid in uuid_dict:  # in the super small chance there is a uuid collision
#                     raise Exception("panic. probability is a lie")
#
#                 v["uuid"] = cur_uuid
#                 new_data[cur_uuid] = v
#                 uuid_dict[cur_uuid] = 1
#
#             # new data should now == old data w. uuid
#             # delete old file and rewrite with new data
#             with open(new_file_name, "w") as f:
#                 json.dump(new_data, f, indent=4)
#
#             with open(new_file_name, "rb") as f:
#                 redata = json.load(f)
#
#             for k in new_data:
#                 if k not in redata:
#                     raise Exception("UUID error!")
