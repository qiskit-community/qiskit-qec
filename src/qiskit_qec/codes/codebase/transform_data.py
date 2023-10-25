#
# (C) Copyright IBM 2017, 2020
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
""" This is the code to transform the stabilizer code database """

import csv
import json
import os

import qiskit_qec.utils.pauli_rep as rep
from qiskit_qec.linear.symplectic import symplectic_product

# standard arguments
N = "n"  # pylint: disable=invalid-name
K = "k"  # pylint: disable=invalid-name
INDEX = "index"
D = "d"  # pylint: disable=invalid-name
AUT_GROUP_SIZE = "aut_group_size"
IS_CSS = "is_css"
IS_DECOMPOSABLE = "is_decomposable"
IS_DEGENERATE = "is_degenerate"
IS_GF4LINEAR = "is_gf4linear"
IS_TRIORTHOGONAL = "is_triorthogonal"
IS_SUBSYSTEM = "is_subsystem"
LOGICAL_OPS = "logical_ops"
ISOTROPIC_GEN = "isotropic_generators"
WEIGHT_ENUMERATOR = "weight_enumerator"
GAUGE_GROUP = ("gauge_group",)
CITATION = "citation"
NAME = "name"
UUID = "uuid"
BELONGS_TO_STANDARD_CODEBASE = "belongs_to_standard_codebase"
CODE_TYPE = "code_type"

REQUIRED_CODE_INFO = {N, K, ISOTROPIC_GEN}
REQUIRED_CODE_STORAGE_INFO = {N, K, UUID, ISOTROPIC_GEN}
VALID_CODE_INFO_FIELDS = {
    N,
    K,
    INDEX,
    D,
    AUT_GROUP_SIZE,
    IS_CSS,
    IS_DECOMPOSABLE,
    IS_DEGENERATE,
    IS_GF4LINEAR,
    IS_TRIORTHOGONAL,
    IS_SUBSYSTEM,
    LOGICAL_OPS,
    ISOTROPIC_GEN,
    WEIGHT_ENUMERATOR,
    GAUGE_GROUP,
    CITATION,
    NAME,
    CODE_TYPE,
    UUID,
}

PARENT = os.path.dirname(os.path.realpath(__file__))
file_name_template = os.path.join(
    PARENT, "data", "base", "base_data", "n_{}", "codes_n_{}_k_{}.json"
)  # n,n,k

OLD_GROUPSIZE_2 = "group_size_2"
OLD_GROUPSIZE_1 = "group_size_1"
OLD_CSS_FORM = "css_form"
OLD_GF4_LINEAR_FORM = "gf4linear_form"
OLD_LOW_WEIGHT_FORM = "low_weight_form"
OLD_GOTTESMAN_FORM = "gottesman_form"
OLD_CSS_LOGICALS = "css_logicals"
OLD_GF4_LINEAR_LOGICALS = "gf4linear_logicals"
OLD_LOGICALS = "logicals"
OLD_STABILIZER = "stabilizer"

COPY_OVER_DIRECTLY = {
    N,
    K,
    D,
    IS_CSS,
    IS_DECOMPOSABLE,
    IS_DEGENERATE,
    IS_GF4LINEAR,
    IS_TRIORTHOGONAL,
    WEIGHT_ENUMERATOR,
    UUID,
}


def convert_jsons(data, i, csvfile="fails.csv"):
    """Convert the data"""
    new_data = {}

    for index, code_info in data.items():
        n = code_info[N]
        k = code_info[K]
        code_str = f"[[{n},{k},{index}]]"
        new_code_info = {}
        # create any new fields that are missing/require old fields

        if OLD_LOW_WEIGHT_FORM not in code_info:
            print(f"LOW WEIGHT FORM MISSING FOR: {code_str}")
            problem = "LOW WEIGHT FORM MISSING"
            with open(csvfile, "a", newline="", encoding="utf-8") as cur_csvfile:
                csvwriter = csv.writer(
                    cur_csvfile, delimiter=" ", quotechar="|", quoting=csv.QUOTE_MINIMAL
                )

                csvwriter.writerow([code_info[N], code_info[K], index, problem])
            continue

        new_code_info[ISOTROPIC_GEN] = code_info[OLD_LOW_WEIGHT_FORM]
        if code_info.get(IS_CSS, False) == 1 and OLD_CSS_LOGICALS in code_info:
            new_code_info[LOGICAL_OPS] = code_info[OLD_CSS_LOGICALS]
        elif code_info.get(IS_GF4LINEAR, False) == 1 and OLD_GF4_LINEAR_LOGICALS in code_info:
            new_code_info[LOGICAL_OPS] = code_info[OLD_GF4_LINEAR_LOGICALS]
        elif OLD_LOGICALS in code_info:
            new_code_info[LOGICAL_OPS] = code_info[OLD_LOGICALS]

        # Convert strings to upper case and then into index syntax with XZY
        # (also YZX as phases are ignored) tensor form
        new_code_info[ISOTROPIC_GEN] = [item.upper() for item in new_code_info[ISOTROPIC_GEN]]
        matrix_iso, phase_exp = rep.str2symplectic(
            new_code_info[ISOTROPIC_GEN], qubit_order="left-to-right"
        )
        temp = rep.symplectic2str(
            matrix_iso, phase_exp, syntax=1, qubit_order="left-to-right", same_type=False
        )
        new_code_info[ISOTROPIC_GEN] = list(temp)

        if LOGICAL_OPS in new_code_info:
            # Convert strings to upper case and then into index syntax with XZY
            # (also YZX as phases are ignored) tensor form
            new_code_info[LOGICAL_OPS] = [item.upper() for item in new_code_info[LOGICAL_OPS]]
            matrix_logical, _ = rep.str2symplectic(
                new_code_info[LOGICAL_OPS], qubit_order="left-to-right"
            )
            new_code_info[LOGICAL_OPS] = list(
                rep.symplectic2str(matrix_logical, syntax=1, qubit_order="left-to-right")
            )

            for s_index, stab in enumerate(matrix_iso):
                for l_index, log in enumerate(matrix_logical):
                    if symplectic_product(stab, log):
                        problem = f"CODE {code_str}: Logical operators not commuting with stabilizers \
                            - stabilizer: {new_code_info[ISOTROPIC_GEN][s_index]} and \
                            logical: {new_code_info[LOGICAL_OPS][l_index]} do not commute"
                        print(problem)
                        with open(csvfile, "a", newline="", encoding="utf-8") as cur_csvfile:
                            csvwriter = csv.writer(
                                cur_csvfile,
                                delimiter=" ",
                                quotechar="|",
                                quoting=csv.QUOTE_MINIMAL,
                            )

                            csvwriter.writerow(
                                [
                                    code_info[N],
                                    code_info[K],
                                    index,
                                    problem,
                                ]
                            )
                        continue

        if i % 10000 == 0:
            print(f"i is {i}")
        new_code_info[IS_SUBSYSTEM] = 1
        new_code_info[INDEX] = index
        new_code_info[CODE_TYPE] = "StabSubSystemCode"

        if code_info.get(OLD_GROUPSIZE_2, 129834675987349857) == 0:
            # auto_group size
            new_code_info[AUT_GROUP_SIZE] = int(
                code_info[OLD_GROUPSIZE_1]
            )  # if float isn't whole number already then group_size_2 shouldn't == 0 then there was already corruption of the data. #pylint: disable=line-too-long

        for key_field in COPY_OVER_DIRECTLY:
            if key_field in code_info:
                new_code_info[key_field] = code_info[key_field]

        new_data[index] = new_code_info
        i += 1

    return new_data, i


def remove_new_dirs():
    """Remove directories to start again"""
    for n_len in range(11):  # test
        for k_d in range(n_len):
            curnt_file_name = file_name_template.format(n_len, n_len, k_d)
            new_json_file_name = curnt_file_name.split(".")[0] + "-NEW.json"
            try:
                os.remove(new_json_file_name)
            except FileNotFoundError:
                pass


def process_files(csvfile="fails.csv"):
    """Process Files"""
    with open(csvfile, "w", encoding="utf-8"):
        pass

    i = 0
    for n_length in range(6):
        for k_dim in range(n_length):
            cur_file_name = file_name_template.format(n_length, n_length, k_dim)
            if os.path.exists(cur_file_name):
                with open(cur_file_name, "r", encoding="utf-8") as curfile:
                    data = json.load(curfile)

                updated_data, i = convert_jsons(data, i, csvfile)
                new_file_name = cur_file_name.split(".")[0] + "-NEW.json"
                with open(new_file_name, "w", encoding="utf-8") as file_to_overwrite:
                    json.dump(updated_data, file_to_overwrite, indent=1)


if __name__ == "__main__":
    remove_new_dirs()
    process_files()
