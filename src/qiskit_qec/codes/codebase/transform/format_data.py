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

# Pull the data files in the data_input directory in the src/qiskit_qec/codes/codebase directory
# run python format_data.py
# Will not modify .gz files - manually do or first unzip


import os
import shutil

# Assuming the current directory contains files like 'data_A_B.json'
current_directory = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(current_directory, "input_data")
files = os.listdir(data_dir)

# Loop through files in the current directory
for filename in files:
    if filename.endswith(".json"):
        # Extract integers n and k from the filename
        try:
            parts = filename.split("_")
            N = int(parts[1])
            K = int(parts[2][:-5])

            # New filename
            new_filename = f"codes_n_{N}_k_{K}.json"

            # Rename the file
            os.rename(os.path.join(data_dir, filename), os.path.join(data_dir, new_filename))

            # Create subdirectory if it doesn't exist
            subdir_path = os.path.join(data_dir, f"n_{N}")
            os.makedirs(subdir_path, exist_ok=True)

            # Move the renamed file to the subdirectory
            shutil.move(
                os.path.join(data_dir, new_filename), os.path.join(subdir_path, new_filename)
            )

        except Exception as e:
            print(f"Error processing file {filename}: {e}")

print("Files renamed and moved successfully.")


if __name__ == "__main__":
    test_prog = True
