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
"""Codes Librarian and Library module"""

import os
import json
import warnings
from typing import Union, List, Optional, Tuple, Dict, Callable
import configparser
from itertools import product

from qiskit.exceptions import QiskitError

from qiskit_qec.codes.code import Code
from qiskit_qec.codes.stabsubsystemcodes import StabSubSystemCode
from qiskit_qec.structures.gauge import GaugeGroup
from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.info.properties import Properties


class CodeLibrary:
    """CodeLibrary class"""

    def __init__(self, name: str, path: str, config_filename: str) -> None:
        """Code library class.

        Args:
            name (str): name of codebase
            path (str): path to codebase
            config_file (str): configuration file nme
        """
        # Set the name of the library
        self.name = name

        # Set the absolute path to the library, config_filename, ...
        if not os.path.isabs(path):
            raise QiskitError(
                f"Path provided must be an absolute path to library directory: {path}"
            )
        self.path = path

        self.config_file = os.path.join(self.path, config_filename)
        self.config_filename = config_filename
        if not os.path.isfile(self.config_file):
            raise QiskitError(f"Input config file: {self.config_file} is not a file")

        # Load the config file and parser
        _library_config = configparser.ConfigParser()
        _library_config.read(self.config_file)

        # Set the location fort the data files (data_path)
        try:
            data_dir = _library_config["Data Directory"]["data_dir"]
        except KeyError:
            data_dir = "."

        if data_dir == ".":
            self.data_path = self.path
        else:
            self.data_path = os.path.join(self.path, data_dir)

        # Set the fetch method
        try:
            self.fetch_method = _library_config["Fetch Method"]["fetch_method"]
        except KeyError:
            self.fetch_methods = CodeLibrarian._DEFAULT_FETCH

        # Set the min/max ranges for n, k and index
        try:
            self.n_min = int(_library_config["Ranges"]["n_min"])
            self.n_max = int(_library_config["Ranges"]["n_max"])
        except KeyError as n_minmax_error:
            raise QiskitError(
                "Max and min n not provided in \
                library config.ini"
            ) from n_minmax_error

        try:
            self.k_min = int(_library_config["Ranges"]["k_min"])
            self.k_max = int(_library_config["Ranges"]["k_max"])
        except KeyError as k_minmax_error:
            raise QiskitError(
                "Max and min k not provided \
                in library config.ini"
            ) from k_minmax_error

        try:
            self.index_min = int(_library_config["Ranges"]["index_min"])
            self.index_max = int(_library_config["Ranges"]["index_max"])
        except KeyError as index_minmax_error:
            raise QiskitError(
                "Max and min index not provided \
                in library config.ini"
            ) from index_minmax_error

        # Set the n data director name template and the n,k data file name templates
        try:
            self.n_dir_format = _library_config["File Formats"]["n_dir_format"]
        except KeyError:
            self.n_dir_format = CodeLibrarian.n_dir_format

        try:
            self.n_k_codes_file_format = _library_config["File Formats"]["n_k_codes_file_format"]
        except KeyError:
            self.n_k_codes_file_format = CodeLibrarian.n_k_codes_file_format

        # Load data into dictionary object if required
        if self.fetch_method == "memory":
            self.data = {}
            for n in range(self.n_min, self.n_max + 1, 1):
                self.data[n] = {}
                for k in range(self.k_min, self.k_max + 1, 1):
                    file_to_load = os.path.join(
                        self.data_path,
                        self.n_dir_format.format(n),
                        self.n_k_codes_file_format.format(n, k),
                    )
                    # check that we have this file
                    if os.path.exists(file_to_load):
                        with open(file_to_load, mode="r", encoding="utf-8") as codes_json_file:
                            self.data[n][k] = {
                                int(index): data
                                for index, data in json.load(codes_json_file).items()
                            }

    @staticmethod
    def data2code(**record) -> Code:
        """_summary_

        Returns:
            Code: _description_
        """
        if record[Properties.TYPE] == "StabSubSystemCode":
            gauge_group = GaugeGroup(
                isotropic_generators=record.get(Properties.ISOTROPIC_GEN, None),
                hyperbolic_generators=record.get(Properties.HYPERBOLIC_GEN, None),
            )
            return StabSubSystemCode(gauge_group)
        else:
            raise QiskitError(f"Unsupported code type: {record[Properties.TYPE]}")

    def search(
        self,
        n: Union[int, List[int], None] = None,
        k: Union[int, List[int], None] = None,
        index: Union[int, List[int], None] = None,
        info_only=False,
        **kwargs,
    ) -> Union[List[Code], List[Tuple[int, int, int]]]:
        """Search method for library"""

        if isinstance(n, int) and isinstance(k, int) and isinstance(index, int):
            single_search = True
            n_range = n
            k_range = k
            index_range = index
        else:
            single_search = False

            if n is None:
                n_range = range(self.n_min, self.n_max + 1, 1)
            elif isinstance(n, int):
                n_range = [n]
            else:
                n_range = n

            if k is None:
                k_range = range(self.k_min, self.k_max + 1, 1)
            elif isinstance(k, int):
                k_range = [k]
            else:
                k_range = k

            if isinstance(index, int):
                index_range = [index]
            else:
                index_range = index

        return self._search(
            n_range,
            k_range,
            index_range,
            info_only=info_only,
            single_search=single_search,
            **kwargs,
        )

    def _search(
        self,
        n_range: List[int],
        k_range: List[int],
        index_range: List[int],
        info_only: bool,
        single_search: bool,
        **kwargs,
    ) -> Union[List[Code], List[Tuple[int, int, int]]]:
        """Search method for library"""

        codes = []

        if self.fetch_method == "memory":
            # If we have on (n,k,index) pair then simply return the requested code data
            if single_search:
                try:
                    code_data = self.data[n_range][k_range][index_range]
                    if info_only:
                        codes.append(Properties(**code_data))
                    else:
                        codes.append(self.data2code(**code_data))
                except KeyError:
                    pass
                return codes

            for n, k in product(n_range, k_range):
                try:
                    data = self.data[n][k]
                    if index_range is not None:
                        data = {ind: data[ind] for ind in index_range if ind in data.keys()}

                    code_data = [
                        entry for entry in data.values() if kwargs.items() <= entry.items()
                    ]

                    if info_only:
                        codes = codes + [Properties(**entry) for entry in code_data]
                    else:
                        codes = codes + [self.data2code(**entry) for entry in code_data]
                except KeyError:
                    continue
            return codes
        else:
            raise QiskitError("Only mempory method currently implemented")

    def in_range(self, n, k, index) -> bool:
        """_summary_

        Args:
            n (_type_): _description_
            k (_type_): _description_
            index (_type_): _description_

        Returns:
            bool: _description_
        """
        if n is None:
            n_value = True
        else:
            if isinstance(n, List):
                n_min = min(n)
                n_max = max(n)
            else:
                n_min = n
                n_max = n

            n_value = self.n_min <= n_min and n_max <= self.n_max

        if k is None:
            k_value = True
        else:
            if isinstance(k, List):
                k_min = min(k)
                k_max = max(k)
            else:
                k_min = k
                k_max = k

            k_value = self.k_min <= k_min and k_max <= self.k_max

        if index is None:
            index_value = True
        else:
            if isinstance(index, List):
                index_min = min(index)
                index_max = max(index)
            else:
                index_min = index
                index_max = index

            index_value = self.index_min <= index_min and index_max <= self.index_max

        return n_value and k_value and index_value


class CodeLibrarian:
    """CodeLibrarian class."""

    _codebase_path = ""
    _data_path = ""
    _library_dirs = []

    _fetch_methods = {}

    n_dir_format = ""
    n_k_codes_file_format = ""

    _CONFIG_FILE = "config.ini"
    _DEFAULT_FETCH = "http"
    _DEFAULT_ORDER = "hfm"
    _ORDER_DICT = {"h": "http", "f": "file", "m": "memory"}

    def __init__(self, config: Optional[str] = None) -> None:
        """CodeLibrarian class.

        Args:
            config: configuration
        """

        # Get absolute path to codebase directory
        CodeLibrarian._codebase_path = os.path.abspath(os.path.dirname(__file__))
        default_config = os.path.join(CodeLibrarian._codebase_path, CodeLibrarian._CONFIG_FILE)

        # Read the config file to determine which method
        # (memory, file or http) to use for each
        # code library. Default is order is http, file, memory
        _config = configparser.ConfigParser()
        if config is None:
            _config.read(default_config)
        else:
            _config.read(config)

        # Set the path to the codebase directory
        try:
            _dir_name = _config["paths"]["codebase_path"]
            CodeLibrarian._codebase_path = _dir_name
        except KeyError:
            pass

        try:
            CodeLibrarian._data_path = os.path.join(
                CodeLibrarian._codebase_path, _config["paths"]["data_path"]
            )
        except KeyError as data_path_error:
            raise QiskitError(
                "A data path (data_path) must be set in \
                the config.ini file in the [paths] section"
            ) from data_path_error

        # Set the order in which to search for code libraries
        try:
            _order = _config["order"]["order"]
        except KeyError:
            _order = CodeLibrarian._DEFAULT_ORDER

        CodeLibrarian._order = [CodeLibrarian._ORDER_DICT[letter] for letter in _order]

        # Get list of library code directories
        joined = [
            os.path.join(CodeLibrarian._data_path, item)
            for item in os.listdir(CodeLibrarian._data_path)
        ]
        CodeLibrarian._library_dirs = [path for path in joined if os.path.isdir(path)]

        try:
            CodeLibrarian.n_dir_format = _config["Default File Formats"]["n_dir_format"]
        except KeyError as n_dir_format_error:
            raise QiskitError(
                "Default n_dir_format not specified in config.ini file. \
                See [File Formats] section."
            ) from n_dir_format_error

        # Set the file format templates
        try:
            CodeLibrarian.n_k_codes_file_format = _config["Default File Formats"][
                "n_k_codes_file_format"
            ]
        except KeyError as n_k_codes_file_format_error:
            raise QiskitError(
                "Default n_k_codes_file_format not specified in config.ini \
                file. See [File Formats] section."
            ) from n_k_codes_file_format_error

        # Create the CodeLibrarys
        self.libraries = {}
        for library_dir in CodeLibrarian._library_dirs:
            library_path = os.path.join(CodeLibrarian._data_path, library_dir)
            self.libraries[library_dir] = CodeLibrary(
                library_dir, library_path, CodeLibrarian._CONFIG_FILE
            )

    def get(
        self,
        n: Union[int, List[int], None] = None,
        k: Union[int, List[int], None] = None,
        index: Union[int, List[int], None] = None,
        info_only=False,
        **kwargs,
    ) -> Union[List[Code], List[Tuple[int, int, int]]]:
        """_summary_

        Args:
            n (Optional[int,List[int]], optional): _description_. Defaults to None.
            k (Optional[int,List[int]], optional): _description_. Defaults to None.
            index (Optional[int,List[int]], optional): _description_. Defaults to None.
            info_only (bool, optional): _description_. Defaults to False.
            **kwargs: _de_

        Returns:
            Union[List[Code],List[Tuple[int,int,int]]]: _description_
        """
        codes = []
        for library in self.libraries.values():
            if library.in_range(n, k, index):
                codes = codes + library.search(n, k, index, info_only=info_only, **kwargs)

        return codes


_codebase = CodeLibrarian()


def small_code(n: int, k: int, index: int, info_only: bool = False) -> Union[Code, Properties]:
    """Returns the [[n,k]]-index code in the library.

    Args:
        n: number of physical qubits
        k: number of encoded qubits
        index: Unique index for the [[n,k]] code
        info_only: If True then only a Properties object will be returned containing
        information about the code requested. If False then a Code will be returned in
        a suitable Code subclass. Defaults to False.

    Raises:
        QiskitError: n, number of physical qubits, must be specified.
        QiskitError: k, number of encoded qubits, must be specified.
        QiskitError: index, index assigned to the requested [[n,k]] code, must be specified.

    Returns:
        Union[Code, Properties]: Returns the [[n,k]]-index Code if info_only is False and a Properties
        object if the info_only is True.

    Example:
        >>> code_info = small_code(5,2,0, info_only=True)
        >>> code_info.info
            [[5,2]]-0 StabSubSystemCode
            -------------------------------------------------------------------------------
            isotropic_generators : ['X0X2', 'Z1Z4', 'Z0Z2']
            logical_ops          : ['X3', 'X1X4', 'Z3', 'Z4']
            is_subsystem         : 1
            index                : 0
            code_type            : StabSubSystemCode
            aut_group_size       : 576
            is_decomposable      : 1
            weight_enumerator    : [1, 0, 4, 0, 3, 0]
            is_css               : 1
            uuid                 : 47e52342-8ffc-48e5-be51-c60aab28e9b8
            is_degenerate        : 0
            d                    : 1
            is_gf4linear         : 0
            k                    : 2
            n                    : 5

        >>> code = small_code(5,2,0)
        >>> code.gaugegroup.generators
            PauliList(['X0X2', 'Z1Z4', 'Z0Z2'])

    Note:
        This is the equivalent to the GAP method SmallGroup but for codes
    """
    if not isinstance(n, int):
        raise QiskitError("n, number of physical qubits, must be specified.")

    if not isinstance(k, int):
        raise QiskitError("k, number of encoded qubits, must be specified.")

    if not isinstance(index, int):
        raise QiskitError("index, index assigned to the requested [[n,k]] code, must be specified.")

    return _small_code(n, k, index, info_only=info_only)


def _small_code(n: int, k: int, index: int, info_only: bool) -> Union[Code, Properties]:
    """Returns the [[n,k]]-index code in the library.

    Args:
        n: number of physical qubits
        k: number of encoded qubits
        index: Unique index for the [[n,k]] code
        info_only: If True then only a Properties object will be returned containing
        information about the code requested. If False then a Code will be returned in
        a suitable Code subclass.

    Returns:
        Union[Code, Properties]: Returns the [[n,k]]-index Code if info_only is False and a Properties
        object if the info_only is True.

    """

    return _all_small_codes(n=n, k=k, index=index, info_only=info_only, list_only=False)


def all_small_codes(
    n: Union[int, List[int], None] = None,
    k: Union[int, List[int], None] = None,
    index: Union[int, List[int], None] = None,
    info_only: bool = False,
    list_only: bool = False,
    **kwargs,
) -> Union[Code, List[Code], Properties, List[Properties]]:
    """Returns all small codes.

    Args:
        n (Union[int,List[int],None], optional): _description_. Defaults to None.
        k (Union[int,List[int],None], optional): _description_. Defaults to None.
        index (Union[int,List[int],None], optional): _description_. Defaults to None.
        info_only (bool, optional): _description_. Defaults to False.
        list_only (optional): If True then a list will always be returned. Defaults to False
        **kwargs: desc

    Returns:
        Union[Code, List[Code], Properties, List[Properties]]: codes or properties matches query

    Example:
        >>> all_small_codes(4)
            [<qiskit_qec.codes.stabsubsystemcodes.StabSubSystemCode at 0x13e304f40>,
             <qiskit_qec.codes.stabsubsystemcodes.StabSubSystemCode at 0x13e30a490>,
             <qiskit_qec.codes.stabsubsystemcodes.StabSubSystemCode at 0x13e30a880>,
             <qiskit_qec.codes.stabsubsystemcodes.StabSubSystemCode at 0x13e30aa60>]

        >>> code = all_small_codes(4,is_css=False)
        >>> code.gauge_group.generators
            PauliList(['Z0X1Z2', 'Y0Y1X3', 'Z1X2Z3'])

        >>> code_info = all_small_codes(4,is_css=False, info_only=True)
        >>> code_info.info
            [[4,1]]-8 StabSubSystemCode
            -------------------------------------------------------------------------------
            isotropic_generators : ['Z0X1Z2', 'Y0Y1X3', 'Z1X2Z3']
            logical_ops          : ['Z2X3', 'Z0Z3']
            is_subsystem         : 1
            index                : 8
            code_type            : StabSubSystemCode
            aut_group_size       : 24
            is_decomposable      : 0
            weight_enumerator    : [1, 0, 0, 4, 3]
            is_css               : 0
            uuid                 : e6e4edd6-2ec6-467f-9187-0cc64bc51f1a
            is_degenerate        : 0
            d                    : 2
            is_gf4linear         : 0
            k                    : 1
            n                    : 4

    Note:
        This is the equivalent to the GAP method AllSmallGroups but for codes
    """
    if n is not None:
        if not isinstance(n, (int, list)):
            raise QiskitError(f"n must be an integer or a list of integers: {n}")

    if k is not None:
        if not isinstance(k, (int, list)):
            raise QiskitError(f"k must be an integer or a list of integers: {k}")

    if index is not None:
        if not isinstance(index, (int, list)):
            raise QiskitError(f"index must be an integer or a list of integers: {index}")

    info_only = bool(info_only)
    list_only = bool(list_only)

    return _all_small_codes(
        n=n, k=k, index=index, info_only=info_only, list_only=list_only, **kwargs
    )


def _all_small_codes(
    n: Union[int, List[int], None],
    k: Union[int, List[int], None],
    index: Union[int, List[int], None],
    info_only,
    list_only,
    **kwargs,
) -> Union[Code, List[Code], Properties, List[Properties]]:
    """_summary_

    Args:
        n (Union[int,List[int],None], optional): _description_. Defaults to None.
        k (Union[int,List[int],None], optional): _description_. Defaults to None.
        index (Union[int,List[int],None], optional): _description_. Defaults to None.
        info_only (bool, optional): _description_. Defaults to False.
        list_only (optional): If True then a list will always be returned. Defaults to False
        **kwargs: desc

    Returns:
        Union[Code, List[Code], Properties, List[Properties]]: _description_
    """
    try:
        result = _codebase.get(n, k, index, info_only=info_only, **kwargs)
    except NameError as err:
        raise QiskitError("Codebase not initilized") from err

    if list_only or len(result) > 1 or len(result) == 0:
        return result
    else:
        return result[0]
