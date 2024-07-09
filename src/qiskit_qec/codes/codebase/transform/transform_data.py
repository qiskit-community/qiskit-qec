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
import re


# This program will:
# Change the filenames : data_N_K.json -> codes_n_N_k_K.json
# Change lists in data_N_K.json to dictionary with index from 0 to # of codes - 1
# Shift indices -1 : so index from 1 to index from 0
#   for aut_group_generators (perms and qubits)
#   isotropic_generators
#   logical_ops

# standard arguments
N = "n"  # pylint: disable=invalid-name
K = "k"  # pylint: disable=invalid-name
INDEX = "index"
D = "d"  # pylint: disable=invalid-name
AUT_GROUP_SIZE = "aut_group_size"
AUT_GROUP_GENERATORS = "aut_group_generators"
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

VALID_CODE_INFO_FIELDS = {
    N,
    K,
    INDEX,
    D,
    AUT_GROUP_SIZE,
    AUT_GROUP_GENERATORS,
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
file_name_template = os.path.join(PARENT, "input_data", "n_{}", "codes_n_{}_k_{}.json")  # n, n, k

COPY_OVER_DIRECTLY = {N, K, D, IS_SUBSYSTEM, WEIGHT_ENUMERATOR, UUID, AUT_GROUP_SIZE}

REDUCE_INDEX = {AUT_GROUP_GENERATORS, ISOTROPIC_GEN, LOGICAL_OPS}

BOOLEAN_KEYS = {IS_CSS, IS_DECOMPOSABLE, IS_DEGENERATE, IS_GF4LINEAR, IS_SUBSYSTEM}


def indices_minus_one(in_string: str) -> str:

    # Use regular expression to find numeric and non-numeric sequences
    segments = re.findall(r"(\d+|[^0-9]+)", in_string)

    # Decrease numbers by 1 and throw an exception for negatives
    for i in range(len(segments)):
        if segments[i].isdigit():
            segments[i] = int(segments[i]) - 1
            if segments[i] < 0:
                raise ValueError("Error: Negative number encountered after decreasing.")

    # Convert numbers back to strings
    segments = [str(segment) for segment in segments]

    # Join all entries into a single string
    new_string = "".join(segments)

    return new_string


def convert_jsons(data, i, csvfile="fails.csv"):
    """Convert the data"""

    # Convert list to dictionary with index
    keys = range(len(data))
    data = dict(zip(keys, data))

    for index, code_info in data.items():
        n = code_info[N]
        k = code_info[K]
        code_str = f"[[{n},{k}]]:{index}"

        # Update aut_group_generators
        # shift index by -1

        for element in REDUCE_INDEX:
            try:
                code_info[element] = [indices_minus_one(item) for item in code_info[element]]
            except ValueError as err:
                print(f"Error on index reduction for code {code_str}")
                problem = "NEGATIVE INDEX"
                with open(csvfile, "a", newline="", encoding="utf-8") as cur_csvfile:
                    csvwriter = csv.writer(
                        cur_csvfile, delimiter=" ", quotechar="|", quoting=csv.QUOTE_MINIMAL
                    )
                    csvwriter.writerow([code_info[N], code_info[K], index, problem])
                continue

        # Update Booleans to integer Boolean

        for element in BOOLEAN_KEYS:
            code_info[element] = int(code_info.get(element))

        if i % 10000 == 0:
            print(f"Processed {i} records")
        i += 1

        data[index] = code_info

    return data, i


def remove_new_files(test=False):
    """Remove directories to start again"""
    print("Remove directories to start again")
    for n_len in range(10):
        for k_d in range(n_len):
            curnt_file_name = file_name_template.format(n_len, n_len, k_d)
            new_json_file_name = (
                curnt_file_name[:-5] + "-NEW.json"
            )  # Assumes file name ends with .json
            try:
                print(f"curnt_file_name = {curnt_file_name}")
                print(f"Attempting to remove file : {new_json_file_name}")
                os.remove(new_json_file_name)
            except FileNotFoundError:
                pass


def process_files(csvfile="fails.csv", test=False):
    """Process Files"""
    print("Processing files")
    with open(csvfile, "w", encoding="utf-8"):
        pass

    i = 0  # Only used for showing the status of the process, counts total records processed
    for n_length in range(10):
        for k_dim in range(n_length):
            cur_file_name = file_name_template.format(n_length, n_length, k_dim)
            # print(f"Checking file : {cur_file_name}")
            if os.path.exists(cur_file_name):
                # print(f"Working with found file:")
                with open(cur_file_name, "r", encoding="utf-8") as curfile:
                    data = json.load(curfile)
                updated_data, i = convert_jsons(data, i, csvfile)
                new_file_name = (
                    cur_file_name[:-5] + "-NEW.json"
                )  # Assumes file name ends with .json
                with open(new_file_name, "w", encoding="utf-8") as file_to_overwrite:
                    json.dump(updated_data, file_to_overwrite, indent=1)


if __name__ == "__main__":
    test_prog = False
    print(f"Transform database: test mode is: {test_prog}")
    if test_prog:
        print(f"PARENT = {PARENT}")
    remove_new_files(test=test_prog)
    process_files(test=test_prog)
