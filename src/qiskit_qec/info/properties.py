# This code is part of Qiskit.
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

"""Code Properties Class."""

from typing import Dict


class Properties:
    """Code Properties Class."""

    N = "n"  # pylint: disable=invalid-name
    K = "k"  # pylint: disable=invalid-name
    INDEX = "index"
    D = "d"  # pylint: disable=invalid-name

    # The following fields should eventually be moved to an external schema
    # that is used by the web server database as well as this class. Here
    # for now to get things moving. Also these fields are likely to change
    # as the final taxonomy is developed. Therefore code should not rely
    # on the exact naming etc used below. So for example a used should
    # not use CodeLibrarian.IS_CSS but instead use the string "is_css".

    # Boolean Fields

    IS_CSS = "is_css"
    IS_DECOMPOSABLE = "is_decomposable"
    IS_DEGENERATE = "is_degenerate"
    IS_GF4LINEAR = "is_gf4linear"
    IS_TRIORTHOGONAL = "is_triorthogonal"

    IS_CSS_KEY = 1 << 0
    IS_DECOMPOSABLE_KEY = 1 << 1
    IS_DEGENERATE_KEY = 1 << 2
    IS_GF4LINEAR_KEY = 1 << 3
    IS_TRIORTHOGONAL_KEY = 1 << 4

    _name_to_key = {
        IS_CSS: IS_CSS_KEY,
        IS_DECOMPOSABLE: IS_DECOMPOSABLE_KEY,
        IS_DEGENERATE: IS_DEGENERATE_KEY,
        IS_GF4LINEAR: IS_GF4LINEAR_KEY,
        IS_TRIORTHOGONAL: IS_TRIORTHOGONAL_KEY,
    }

    _key_to_name = {
        IS_CSS_KEY: IS_CSS,
        IS_DECOMPOSABLE_KEY: IS_DECOMPOSABLE,
        IS_DEGENERATE_KEY: IS_DEGENERATE,
        IS_GF4LINEAR_KEY: IS_GF4LINEAR,
        IS_TRIORTHOGONAL_KEY: IS_TRIORTHOGONAL,
    }

    # Non-Boolean Fields

    LOGICAL_OPS = "logical_ops"
    STABILIZER = "stabilizer"
    ISOTROPIC_GEN = "isotropic_generators"
    HYPERBOLIC_GEN = "hyperbolic_generators"
    WEIGHT_ENUMERATOR = "weight_enumerator"
    GAUGE_GROUP = "gauge_group"
    AUT_GROUP_SIZE = "aut_group_size"
    CITATION = "citation"
    NAME = "name"
    UUID = "uuid"

    # Code Types : Used to reconstruct the code
    TYPE = "code_type"
    CODE_TYPES = {"StabSubSystemCode": "A symplectic matrix is necessary to create a Code"}

    def __init__(self, **prop_dict: Dict):
        self.properties = prop_dict

    def __str__(self) -> str:
        info_str = ""
        for key, value in self.properties.items():
            info_str = info_str + f"{key:20} : {value}" + "\n"
        return info_str

    def __repr__(self):
        info_str = "{"
        for key, value in self.properties.items():
            info_str = info_str + f"{key:20} : {value}," + "\n"
        return info_str + "}"

    def __getitem__(self, value):
        return self.properties[value]

    def __setitem__(self, key, value):
        self.properties[key] = value

    @property
    def info(self):
        """Print info"""
        print(
            f"[[{self.properties[Properties.N]},"
            + f"{self.properties[Properties.K]}]]-{self.properties[Properties.INDEX]}"
            + f" of type {self.properties[Properties.TYPE]}"
        )
        line = "-" * 79
        print(line)
        print(str(self))
