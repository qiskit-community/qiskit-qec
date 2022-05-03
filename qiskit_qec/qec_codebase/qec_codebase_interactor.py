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
# Part of the QEC framework
"""Module for interacting with codebase"""

import copy
import json
import logging
import os
import shutil
import uuid
from typing import Any, Dict, List, Set, Union

from qiskit_qec.codes.subsystemcodes import SubSystemCode
from qiskit_qec.operators.pauli_list import PauliList
from qiskit_qec.py_utils.exceptions import QiskitQECError
from qiskit_qec.structures.gauge import GaugeGroup


class QECCodeBase:
    """QEC Codebase"""

    INSTANCE_NUMBER = [0]

    # standard arguments
    N = "n"  # pylint: disable=invalid-name
    K = "k"  # pylint: disable=invalid-name
    D = "d"  # pylint: disable=invalid-name
    AUT_GROUP_SIZE = "aut_group_size"
    IS_CSS = "is_css"
    IS_DECOMPOSABLE = "is_decomposable"
    IS_DEGENERATE = "is_degenerate"
    IS_GF4LINEAR = "is_gf4linear"
    IS_TRIORTHOGONAL = "is_triorthogonal"
    IS_SUBSYSTEM = "is_subsystem"
    LOGICAL_OPS = "logical_ops"
    STABILIZER = "stabilizer"
    WEIGHT_ENUMERATOR = "weight_enumerator"
    GAUGE_GROUP = "gauge_group"
    CITATION = "citation"
    NAME = "name"
    QEC_UUID = "uuid"
    BELONGS_TO_STANDARD_CODEBASE = "belongs_to_standard_codebase"

    REQUIRED_CODE_INFO = {N, K, STABILIZER}
    REQUIRED_CODE_STORAGE_INFO = {N, K, QEC_UUID, STABILIZER}
    VALID_CODE_INFO_FIELDS = {
        N,
        K,
        D,
        AUT_GROUP_SIZE,
        IS_CSS,
        IS_DECOMPOSABLE,
        IS_DEGENERATE,
        IS_GF4LINEAR,
        IS_TRIORTHOGONAL,
        IS_SUBSYSTEM,
        LOGICAL_OPS,
        STABILIZER,
        WEIGHT_ENUMERATOR,
        GAUGE_GROUP,
        CITATION,
        NAME,
        QEC_UUID,
        BELONGS_TO_STANDARD_CODEBASE,
    }

    CODE_INFO_TYPING = {
        N: (int,),
        K: (int,),
        D: (int,),
        AUT_GROUP_SIZE: (int,),
        IS_CSS: (int,),
        IS_DECOMPOSABLE: (int,),
        IS_DEGENERATE: (int,),
        IS_GF4LINEAR: (int,),
        IS_TRIORTHOGONAL: (int,),
        IS_SUBSYSTEM: (int,),
        LOGICAL_OPS: (
            list,
            str,
        ),
        STABILIZER: (
            list,
            str,
        ),
        WEIGHT_ENUMERATOR: (
            list,
            int,
        ),
        CITATION: (str,),
        NAME: (str,),
        QEC_UUID: (str,),
        BELONGS_TO_STANDARD_CODEBASE: (int,),
    }

    def __init__(self, codebase_parent_path: str):
        """Interface for QEC Codebase.
        Currently the actual directories making up the
        codebase reside in https://github.ibm.com/Grace-Harper/ibm-codebase-repo
        It must be downloaded and it's local path given
        as the input parameter, codebase_parent_path
        in order for QECCodeBase to work.

        Args:
            codebase_parent_path (str): Path to where codebase is being stored locally

        Raises:
            QiskitError: Something went wrong.
        """
        self.logger = logging.getLogger(__name__)
        self.INSTANCE_NUMBER[0] += 1
        if self.INSTANCE_NUMBER[0] > 1:
            self.logger.warning(  # pylint: disable=logging-fstring-interpolation
                f"YOU SHOULD PROBABLY *NOT* CREATE MORE THAN 1 QECCodeBase. "
                f"You are trying to create {self.INSTANCE_NUMBER} instances. "
                f"THERE IS NO CONCURRENCY PROTECTION!"
            )
        self.logger.warning(
            "Database is a work in progress. Please read the DATABASE_NOTES in the top folder of the"
            " qiskit-qec software"
        )
        # parent_name = pathlib.Path(__file__).parent.resolve()
        self.parent_name = codebase_parent_path
        self.qec_standard_subsystem_codebase_name = (  # pylint: disable=invalid-name
            "qec_standard_subsystem_codebase"
        )
        self.codebase_path = os.path.join(
            self.parent_name, self.qec_standard_subsystem_codebase_name
        )

        self.qec_playground_codebase_name = "qec_playground_codebase"
        self.playground_path = os.path.join(self.parent_name, self.qec_playground_codebase_name)

        self.code_path_folder = "n_length_{}_codes"
        self.code_path_format = os.path.join(
            self.code_path_folder, "codes_n_{}_k_{}.json"
        )  # n,n,k"

        self.playground_cache = {}

    def get_subsystem_code(
        self,
        n_len: Union[List[int], int],
        k_dim: Union[List[int], int],
        acceptable_missing_params: Set[str] = None,
        allow_standard_db=True,
        allow_playground=True,
        **additional_params,  # pylint: disable=missing-param-doc
    ) -> (List[SubSystemCode]):
        """
        Args:
            n_len (Union[List[int]],int]): length of the code
            k_dim (Union[List[int]],int]): dimension of the code
            acceptable_missing_params (Set[str]): set of parameter
                    keys (example: "is_decomposable")
                    that the user doesn't require exists in the current code.
                    Example: If you have is_decomposable=True
                    in additional_params and an empty acceptable_missing_params
                    then a code without the "is_decomposable" will be discarded.
                    However, if "is_decomposable" is in the acceptable_missing_params set then
                    while codes that have  is_decomposable=False will be discarded,
                    codes that are missing the
                    "is_decomposable" field will not be discarded based solely
                    on their lacking of this field
            allow_standard_db (bool): Include results from standard codebase
            allow_playground (bool): Include results from playground codebase
            **additional_params (**kwargs): additional parameters to make the search more strict

        Returns: List of Subsystem Codes matching the input criteria

        """
        # housekeeping for a cleaner interface
        if isinstance(n_len, int):
            n_len = [n_len]
        if isinstance(k_dim, int):
            k_dim = [k_dim]

        # if additional_params[self.IS_SUBSYSTEM] != 1:
        #     raise QiskitQECError("CANNOT get subsystem code that isn't a subsystem")

        if self.BELONGS_TO_STANDARD_CODEBASE in additional_params:
            self.logger.warning(
                "Are you sure you mean to query whether code belongs to standard codebase?"
            )
        if self.QEC_UUID in additional_params:
            self.logger.warning("Are you sure you mean to query on uuid?")

        qec_code_jsons = {}
        for n in n_len:
            for k in k_dim:
                qec_code_jsons = {
                    **qec_code_jsons,
                    **self._load_code_storage_format_from_n_k(
                        n,
                        k,
                        allow_standard_db=allow_standard_db,
                        allow_playground=allow_playground,
                    ),
                }

        if len(additional_params) > 0:
            self.logger.warning(  # pylint: disable=logging-fstring-interpolation
                f"You have {len(additional_params)} additional params. Additional params may"
                " cause really slow results"
            )

            if acceptable_missing_params is None:
                acceptable_missing_params = set()
            qec_code_jsons = self._additional_param_query(
                qec_code_jsons,
                acceptable_missing_params,
                **additional_params,
            )

        qec_code_list = self._convert_codes_storage_format_to_subsystem_codes(qec_code_jsons)
        return qec_code_list

    def store_new_subsystem_code(
        self,
        code: SubSystemCode,
        force=False,
        keep_uuid=False,
        allow_new_fields=False,
        verify_code_info_for_storage=True,
        validate_code_uniqueness=True,
        flush_cache=True,
    ) -> str:
        """Store new subsystem code in playground codebase

        Args:
            code (SubSystemCode): code to be saved
            force (bool, optional): Override all checks and force code
            to be saved to playground. Defaults to False.
            keep_uuid (bool, optional): Use the same the uuid already stored
            in the code.parameters for storing in the codebase. Defaults to False.
            allow_new_fields (bool, optional): Allow code.parameters to
            contain new fields not currently in codebase schema. Defaults to False.
            verify_code_info_for_storage (bool, optional): Make sure all standard code.parameters
            information follows codebase schema. Defaults to True.
            validate_code_uniqueness (bool, optional): Make sure there does not exist
            another code in the database that is equal or a superset of code. Defaults to True.
            flush_cache (bool): Flush cache by writing codes to playground and emptying cache

        Returns:
            str: str(uuid) of saved code
        """
        storage_formatted_code = self._convert_code_subsystem_to_storage_format(code, keep_uuid)
        retval = self._protected_store_new_code_storage_format(
            storage_formatted_code,
            force_input=force,
            verify_code_info_for_storage=verify_code_info_for_storage,
            validate_code_uniqueness=validate_code_uniqueness,
            allow_new_fields=allow_new_fields,
        )
        if flush_cache:
            self.flush_cache()
        return retval

    def delete_subsystem_code_from_playground(self, code: SubSystemCode):
        """Delete subsystem code from playground codebase

        Args:
            code (SubSystemCode): code to be deleted

        Raises:
            QiskitQECError: if deletion is unsuccessful

        Returns:
            str: uuid of recently deleted code
        """
        storage_code = self._convert_code_subsystem_to_storage_format(
            code
        )  # will create a new uuid
        code_info = next(iter(storage_code.values()))
        # take care of uuid
        if self.QEC_UUID in code.parameters:
            code_info[self.QEC_UUID] = code.parameters[self.QEC_UUID]
        else:
            code_info.pop(self.QEC_UUID)
        n_len = code_info[self.N]
        k_dim = code_info[self.K]
        new_playground_nk_codes = self._load_code_storage_format_from_n_k(
            n_len, k_dim, allow_standard_db=False
        )
        params_optional_existence = set()
        matching_code = self._additional_param_query(
            new_playground_nk_codes,
            params_optional_existence,
            **code_info,
        )

        if len(matching_code) == 0:
            return ""
        if len(matching_code) > 1:
            raise QiskitQECError(
                f"{code} with parameters: {code.parameters} matched more than 1 code in the codebase:"
                f" {matching_code}"
            )

        dead_uuid = next(iter(matching_code.keys()))
        new_playground_nk_codes.pop(dead_uuid)  # remove it from dict

        self._overwrite_codes_to_playground_file(n_len, k_dim, new_playground_nk_codes)
        return dead_uuid

        # return UUID if deleted, else return ""

    def _overwrite_codes_to_playground_file(self, n_len, k_dim, new_playground_nk_codes):
        playground_path = self._evaluate_playground_path(n_len, k_dim)

        if len(new_playground_nk_codes) > 0:
            with open(playground_path, "w", encoding="utf-8") as playground_file:
                json.dump(new_playground_nk_codes, playground_file, indent=4)
        else:
            # delete
            if os.path.exists(playground_path):
                os.remove(playground_path)
            enclosing_folder = self._evaluate_playground_folder_path(n_len)
            if os.path.exists(enclosing_folder):
                if len(os.listdir(enclosing_folder)) == 0:
                    os.rmdir(enclosing_folder)

    def _convert_code_subsystem_to_storage_format(
        self,
        code: SubSystemCode,
        keep_uuid=False,
    ):
        self.logger.warning(
            "NOTICE. All uppercase letters will be forced into lowercase for storage"
        )
        assert code.parameters is not None, "code.parameters should not be empty"
        param_dict_deep_copy = copy.deepcopy(code.parameters)
        stabilizer_info = code.gauge_group.generators.stabilizer_info
        if isinstance(stabilizer_info, list):
            stabilizer_str = stabilizer_info[0].lower()
        elif isinstance(stabilizer_info, str):
            stabilizer_str = stabilizer_info.lower()

        if not isinstance(stabilizer_str, str):
            raise NotImplementedError(
                "Currently can only handle subsystem codes that store pauli strings as pdata, ex:"
                " 'IIIXXX'"
            )
        if self.STABILIZER in param_dict_deep_copy:
            if stabilizer_str not in param_dict_deep_copy[self.STABILIZER]:
                param_dict_deep_copy[self.STABILIZER].append(stabilizer_str)
        else:
            param_dict_deep_copy[self.STABILIZER] = [stabilizer_str]
        if keep_uuid:
            if self.QEC_UUID not in param_dict_deep_copy:
                raise QiskitQECError(f" {self.QEC_UUID} not in subsystem code's parameters")
        else:
            new_uuid = str(uuid.uuid4())
            param_dict_deep_copy[self.QEC_UUID] = new_uuid
        code_uuid = param_dict_deep_copy[self.QEC_UUID]
        code_storage = {code_uuid: param_dict_deep_copy}
        return code_storage

    def _evaluate_codepath(
        self,
        n_len: int,
        k_dim: int,
        allow_standard_db=True,
        allow_playground=True,
    ) -> List[str]:
        """

        Args:
            n_len (int): length of the code
            k_dim (int): dimension of the code

        Returns: List of codebase files that contain n,k codes

        """
        paths = []
        if allow_standard_db:
            paths.append(self._evaluate_standard_codebase_path(n_len, k_dim))

        if allow_playground:
            paths.append(self._evaluate_playground_path(n_len, k_dim))
        return paths

    def _evaluate_playground_path(self, n_len: int, k_dim: int) -> str:
        """"""
        return os.path.join(self.playground_path, self.code_path_format).format(n_len, n_len, k_dim)

    def _evaluate_playground_folder_path(self, n_len: int) -> str:
        return os.path.join(self.playground_path, self.code_path_folder).format(n_len)

    def _evaluate_standard_codebase_path(self, n_len: int, k_dim: int) -> str:

        return os.path.join(self.codebase_path, self.code_path_format).format(n_len, n_len, k_dim)

    def _load_code_storage_format_from_n_k(
        self,
        n_len: int,
        k_dim: int,
        allow_standard_db=True,
        allow_playground=True,
        allow_cache=True,
    ) -> Dict:
        """
        Args:
            n_len (int): length of the code
            k_dim (int): dimension of the code

        Returns: Dictionary of code dictionaries from the codebase
        """
        qec_code_jsons = {}
        code_paths = self._evaluate_codepath(n_len, k_dim, allow_standard_db, allow_playground)
        for path in code_paths:
            try:
                with open(path, "r", encoding="utf-8") as code_file:  # should except if invalid
                    qec_code_jsons = {**json.load(code_file), **qec_code_jsons}
            except FileNotFoundError:
                self.logger.exception(  # pylint: disable=logging-fstring-interpolation
                    f"{path} does not exist"
                )  # pylint: disable=logging-fstring-interpolation
        if allow_cache:
            if (n_len, k_dim) in self.playground_cache:
                qec_code_jsons = {**qec_code_jsons, **self.playground_cache[(n_len, k_dim)]}
        return qec_code_jsons

    def _convert_codes_storage_format_to_subsystem_codes(
        self, qec_code_jsons: Dict
    ) -> List[SubSystemCode]:
        """

        Args:
            qec_code_jsons (Dictionary): Dictionary of code dictionaries from database

        Returns: List of Subsystem codes from dictionaries

        """
        qec_code_list = []
        for cur_code_json in qec_code_jsons.values():
            a_stabilizer = cur_code_json[self.STABILIZER][0].upper()
            cur_code = SubSystemCode(GaugeGroup(PauliList(a_stabilizer)), parameters=cur_code_json)

            qec_code_list.append(cur_code)
        return qec_code_list

    def _additional_param_query(
        self,
        json_codes: Dict[str, Dict],
        acceptable_missing_params: Set[str],
        **additional_params,
    ) -> Dict[str, Dict]:
        """

        Args:
            json_codes (Dict): Dictionary of code dictionaries loaded from the database
            acceptable_missing_params (Set[str]): Set of additional parameters
            such, when trying to find valid matches in the codebase,
            if a codebase code is missing that parameter, it won't be immediately
            invalidated.
            **additional_params (**kwargs): Additional parameters to use for whittling down code results

        Returns: Dictionary of code dictionaries matching the additional parameter specifications

        """
        results = {}
        for k, cur_code_storage_format in json_codes.items():
            if self._is_code_valid_match_capitalization_agnostic(
                cur_code_storage_format,
                acceptable_missing_params,
                additional_params,
            ):
                results[k] = json_codes[cur_code_storage_format["uuid"]]

        return results

    def _is_code_valid_match_capitalization_agnostic(
        self,
        code_info_storage_format: Dict[str, Any],
        acceptable_missing_params: Set,
        additional_params: Dict[str, Any],
    ) -> bool:
        """Checks that every value in (field, value) in the addtl_params
        is either equal to or an element of the code_info_storage_format's field
        Args:
            code_info_storage_format (Dict[str,Any]): The code information dictionary
            acceptable_missing_params (Set): Parameters that can be
            ignored if not in the code information dictionary
            additional_params (Dict[str,Any]): The values that the code information
            needs to contain to be a match

        Returns:
            bool: True is code is valid match. Otherwise, False

        """
        for param, value in additional_params.items():
            if param in code_info_storage_format:
                retrieved_code_value = code_info_storage_format[param]

                if isinstance(retrieved_code_value, type(value)):
                    try:
                        if not retrieved_code_value.lower() == value.lower():
                            return False
                    except AttributeError:
                        if not retrieved_code_value == value:
                            return False
                else:
                    try:
                        if not (
                            value in retrieved_code_value
                            or value.lower() in retrieved_code_value
                            or value.upper() in retrieved_code_value
                        ):
                            return False
                    except (
                        TypeError,
                        AttributeError,
                    ):  # retrieved_code_value probably wasn't iterable
                        return False
            else:
                if (
                    param not in acceptable_missing_params
                ):  # if it's not in the dictionary but needs to be
                    return False
        return True

    def _validate_code_storage_fields(self, code_info_storage_format, allow_new_fields=False):
        """Performs some basic checks
        to make your life a little easier by alerting you if something is
        glaringly off. Side Note: if you have a new field that you feel
        would be a useful addition to the codebase schema, include this
        in the PR."""
        assert (
            len(code_info_storage_format) == 1
        ), "function takes in 1 code in storage format at a time"
        code_uuid = next(
            iter(code_info_storage_format.keys())
        )  # get a random key from dict, but dict len is 1
        assert isinstance(code_uuid, str), "uuid must be strings to be json compatible"
        code_field_dict = code_info_storage_format[code_uuid]
        for field in self.REQUIRED_CODE_STORAGE_INFO:
            if field not in code_field_dict:
                raise QiskitQECError(f"{field} is a required field for all codes")

        if not allow_new_fields:
            for field in code_field_dict.keys():
                if field not in self.VALID_CODE_INFO_FIELDS:
                    raise QiskitQECError(
                        f"{field} is a new field and you have set allow_new_fields to {allow_new_fields}"
                    )

        # verify value types
        for field, value in code_field_dict.items():
            if field in self.VALID_CODE_INFO_FIELDS:
                valid_type = self.CODE_INFO_TYPING[field]
                ctype = valid_type[0]
                if ctype is list:
                    elm_type = valid_type[1]
                    try:
                        assert all(isinstance(elm, elm_type) for elm in value)
                    except AssertionError as assert_error:
                        raise QiskitQECError(
                            f"{field} must be" + "list of {elm_type}"
                        ) from assert_error

                    if isinstance(value[0], str):  # force all strs to be lowercase
                        try:
                            assert all(elm == elm.lower() for elm in value)
                        except AssertionError as assert_error:
                            raise QiskitQECError(
                                f" {field} contains strings that are uppercase"
                            ) from assert_error
                else:
                    if not isinstance(value, ctype):
                        raise QiskitQECError(f"{field} must be of type {ctype}")
                    if (
                        isinstance(value, str)
                        and not field == self.CITATION
                        and not field == self.NAME
                    ):
                        if not value.lower() == value:
                            raise QiskitQECError(f"{field} contains strings that are uppercase")

    def _validate_code_storage_format_uniqueness(self, code_storage_format) -> None:
        """

        Args:
            code_storage_format: code_info in the same format it's used to be stored on disk

        Returns: if code already exists via n,k, uuid

        """
        # check if code already exists via n,k, uuid
        # if already exists, return "already exists" & uuid -> uuid, True
        # check if code already exists via n,k, properties
        # if already exists, return "already exists" & uuid -> uuid, True
        # if existing db code is strict subset of properties (excluding uuid/name):
        # return existing db code -> uuid Subset
        # if existing db code is strict subset of properties w/ "force_new"
        #            -------------------> ask Andrew

        def _is_subset(search_code_info, prex_code):

            for key in code_info.keys():
                if key is not self.QEC_UUID:
                    if key not in prex_code or search_code_info[key] != prex_code[key]:
                        return False
            return True

        self.logger.warning(  # pylint: disable=logging-not-lazy
            "Checking that code is unique' is not sophisticated!"
            + "It checks only if a code in the codebase/playground has an exact match"
            + "for each field (except uuid) in the code being added"
        )
        code_info = next(
            iter(code_storage_format.values())
        )  # get a random key from dict, but dict len is 1
        n_len = code_info[self.N]
        k_dim = code_info[self.K]

        pre_existing_codes = self._load_code_storage_format_from_n_k(n_len, k_dim)
        for _, prex_code_info in pre_existing_codes.items():
            match = _is_subset(code_info, prex_code_info)
            if match:
                raise QiskitQECError(
                    f"Code: {code_info} is not completely unique from other code in the database:"
                    f" {prex_code_info}\n"
                    + "To still save this code rerun with"
                    + "force_input=True which will insert without checking"
                )

    def _store_code_storage_format_in_playground(self, cached_codes_to_be_stored):
        """CHECKS NOTHING. Literally shoves it into database."""

        for code_nk, codes in cached_codes_to_be_stored.items():
            n_len = code_nk[0]
            k_dim = code_nk[1]

            for code_info in iter(codes.values()):
                code_info[self.BELONGS_TO_STANDARD_CODEBASE] = 0

            playground_path = self._evaluate_playground_path(n_len, k_dim)

            new_playground_nk_codes = self._load_code_storage_format_from_n_k(
                n_len, k_dim, allow_standard_db=False, allow_cache=False
            )
            if len(new_playground_nk_codes) == 0:
                new_playground_nk_codes = {}
                enclosing_folder = playground_path[: playground_path.rfind(os.path.sep)]
                if not os.path.exists(enclosing_folder):
                    os.makedirs(enclosing_folder)
            new_playground_nk_codes = {
                **new_playground_nk_codes,
                **codes,
            }
            self._overwrite_codes_to_playground_file(n_len, k_dim, new_playground_nk_codes)

    def _protected_store_new_code_storage_format(
        self,
        code_storage_format: dict,
        verify_code_info_for_storage=False,
        validate_code_uniqueness=True,
        allow_new_fields=False,
        force_input=False,
    ) -> str:
        """Will store to PLAYGROUND only"""
        if not force_input:
            if verify_code_info_for_storage:
                self._validate_code_storage_fields(
                    code_storage_format, allow_new_fields=allow_new_fields
                )
            if validate_code_uniqueness:
                self._validate_code_storage_format_uniqueness(code_storage_format)

        code_info = next(iter(code_storage_format.values()))
        code_nk_key = (code_info[self.N], code_info[self.K])
        if code_nk_key in self.playground_cache:
            self.playground_cache[code_nk_key] = {
                **self.playground_cache[code_nk_key],
                **code_storage_format,
            }
        else:
            self.playground_cache[code_nk_key] = code_storage_format

        # TODO rewriting everything will be very slow -- figure if there's a way to append json
        # hopefully not a problem given playground might be small

        # return uuid of newly stored code
        return next(iter(code_storage_format.keys()))

    def flush_cache(self):
        """Write everything in cache to disk and clear cache"""
        try:
            self._store_code_storage_format_in_playground(self.playground_cache)
            self.playground_cache = {}
        except Exception as cur_except:  # pylint: disable=broad-except
            self.logger.exception(  # pylint: disable=logging-fstring-interpolation
                f"{cur_except}. Please try to reflush-cache. "
                f"Playground might have data partially "
                f"written to it."
            )

    def delete_playground_codebase(self, force=False):
        """Delete all codes in playground codebase

        Args:
            force (bool, optional): If True, actually delete playground codebase. Defaults to False.
        """
        if not force:
            self.logger.warning(
                "Consider committing to git before deleting playground codebase;\n To"
                " COMPLETELY delete playground codebase, rerun with force=True."
            )
            return
        self.playground_cache = {}
        for file in os.listdir(self.playground_path):
            shutil.rmtree(os.path.join(self.playground_path, file))
