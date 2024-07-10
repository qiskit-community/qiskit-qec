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

import os
import shutil
import re

debug = False

current_directory = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(current_directory, "input_data")
dest_dir = os.path.join(current_directory, "data", "base", "base_data")

files_and_dirs = os.listdir(data_dir)

# Filter out files
dirs = [d for d in files_and_dirs if os.path.isdir(os.path.join(data_dir, d))]


def is_valid_dir(s):
    # Define a regular expression pattern to match "n_xxx" where xxx is an integer
    pattern = r"^n_\d+$"

    # Use re.match to check if the string matches the pattern
    if re.match(pattern, s):
        return True
    else:
        return False


def move_files():
    # Loop through directories in the current directory
    for dir in dirs:
        if is_valid_dir(dir):

            # Create new directory

            new_dir = os.path.join(dest_dir, str(dir))
            if debug is False:
                os.makedirs(new_dir, exist_ok=True)

            # move files and rename
            source_dir = os.path.join(data_dir, dir)
            files = os.listdir(source_dir)

            for file in files:
                print(f"Testing {file} : {file[-9:]}")
                if file[-9:] == "-NEW.json":
                    try:
                        old_filename = os.path.join(source_dir, file)

                        # New filename
                        filename = file[:-9] + ".json"

                        # New filename + path
                        new_filename = os.path.join(new_dir, filename)

                        # Rename the file
                        if debug is True:
                            print(f"mv {old_filename} {new_filename}")
                        else:
                            os.rename(old_filename, new_filename)

                    except Exception as e:
                        print(f"Error processing file {file}: {e}")

    print("Files renamed and moved successfully.")


if __name__ == "__main__":
    move_files()
