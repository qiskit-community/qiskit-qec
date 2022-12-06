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
#
# This code is based on the paper: "The XP Stabiliser Formalism: a
# Generalisation of the Pauli Stabiliser Formalism with Arbitrary Phases", Mark
# A. Webster, Benjamin J. Brown, and Stephen D. Bartlett. Quantum 6, 815
# (2022).
"""
N-qubit XPPauli Representation Encodings and Conversion Module
"""

# pylint: disable=invalid-name,anomalous-backslash-in-string
# pylint: disable=bad-docstring-quotes  # for deprecate_function decorator

import numbers
import re
from typing import Any, Iterable, List, Optional, Tuple, Union

import numpy as np
from qiskit.circuit import Gate
from qiskit.quantum_info.operators.scalar_op import ScalarOp
from qiskit.exceptions import QiskitError
from scipy.sparse import csr_matrix


# -------------------------------------------------------------------------------
# Module Variables/States
# -------------------------------------------------------------------------------

# Set the internal encodings
# The internal encoding cannot be changed by changing this constant
# These constants are for reference only and do not change the behavior of
# the XPPauli methods. See [ref] for details on the different encodings
# TODO: Include ref for above.

INTERNAL_TENSOR_ENCODING = "XP"
INTERNAL_PHASE_ENCODING = "w"
INTERNAL_XP_PAULI_ENCODING = INTERNAL_PHASE_ENCODING + INTERNAL_TENSOR_ENCODING

DEFAULT_EXTERNAL_TENSOR_ENCODING = "XP"
DEFAULT_EXTERNAL_PHASE_ENCODING = "w"
DEFAULT_EXTERNAL_XP_PAULI_ENCODING = (
    DEFAULT_EXTERNAL_PHASE_ENCODING + DEFAULT_EXTERNAL_TENSOR_ENCODING
)

# Set the external encodings
# The external encodings may be changed via the phase_rep_format,
# symp_rep_format, or pauli_rep_format methods. Any method changing
# these values must make sure to update the others
# Tensor encodings are: 'XZ', 'XZY', 'ZX', 'YZX'
# Phase encodings are: 'i', '-i', 'is', '-is'
# See [ref] for details on the different encodings
# TODO: Include ref for above.
# TODO update this comment after formats have been finalized.

# w is exp(pi*i/N), as defined in XP Formalism paper. If the phase exponent is
# p, then the operator's phase is w**p, where p is an integer between and
# including 0 and 2N-1. P is diag(1,w**2).
PHASE_ENCODINGS = ["w"]
PHASE_ENCODINGS_IMI = []
PHASE_ENCODINGS_ISMIS = []
TENSOR_ENCODINGS = ["XP"]
Y_TENSOR_ENCODINGS = []
XP_PAULI_ENCODINGS = [
    "wXP",
]
XP_PAULI_ENCODINGS_SPLIT = {
    "wXP": ("w", "XP"),
}

# Different string syntax formats are available. The they are "product" syntax and
# "index" syntax. "product" syntax represents a XPPauli operator of the form
# :math: $p * T_1 \otimes T_2 \otimes ... \otimes T_n$ as :math" $pT1T2T3...Tn$. See the
# following exmaples:
#
# -iX \otimes Y \otimes Z -> -iXYZ
# X \otimes Y \otimes Z \otimes I \otimes I \otimes I  -> XYZII
#
# The index syntax only represents the non identity XPPaulis. Index syntax specifically
# indiciates the index that the XPPauli's are acting on. Following Qiskit's current internal
# indexing:
#
# -iX \otimes Y \otimes Z -> -iZ0Y1X2
# X \otimes Y \otimes Z \otimes I \otimes I \otimes I  -> Z3Y4X5
# TODO update this comment/following REGEX after formats have been finalized.

# -i, -1j, +1, ...
COMPLEX_REGEX = r"[\-+]?1?[ij]?"
# (i,2), (1j, 0), ( +i , 2 ), ...
ENC_I_REGEX = r"\(\s*[+]?[ij]\s*\,\s*[0123]\s*\)"
# (-i,2), (-1j, 0), ( -i , 2 ), ...
ENC_MI_REGEX = r"\(\s*\-1?[ij]\s*\,\s*[0123]\s*\)"
# (i,0)(-1,1), (1j,1)(-1 , 0), ...
ENC_IS_REGEX = r"\(\s*1?[ij]\s*\,\s*[01]\s*\)\s*\(\s*\-1\s*\,\s*[01]\s*\)"
# (-i,0)(-1,1), (-1j,1)(-1, 0), ...
ENC_MIS_REGEX = r"\(\s*\-1?[ij]\s*\,\s*[01]\s*\)\s*\(\s*\-1\s*\,\s*[01]\s*\)"

CAP_STAND_ENC_I_REGEX = r"\(i\,([0123])\)"
CAP_STAND_ENC_MI_REGEX = r"\(-i\,([0123])\)"
CAP_STAND_ENC_IS_REGEX = r"\(i\,([01])\)\(-1\,([01])\)"
CAP_STAND_ENC_MIS_REGEX = r"\(-i\,([01])\)\(-1\,([01])\)"

CAP_STAND_I_PATTERN = re.compile(f"^{CAP_STAND_ENC_I_REGEX}$")
CAP_STAND_MI_PATTERN = re.compile(f"^{CAP_STAND_ENC_MI_REGEX}$")
CAP_STAND_IS_PATTERN = re.compile(f"^{CAP_STAND_ENC_IS_REGEX}$")
CAP_STAND_MIS_PATTERN = re.compile(f"^{CAP_STAND_ENC_MIS_REGEX}$")

CAP_PHASE_PATTERN = {
    "i": CAP_STAND_I_PATTERN,
    "-i": CAP_STAND_MI_PATTERN,
    "is": CAP_STAND_IS_PATTERN,
    "-is": CAP_STAND_MIS_PATTERN,
}

COMPLEX_PATTERN = re.compile(f"^({COMPLEX_REGEX})$")

I_PATTERN = re.compile(f"^({ENC_I_REGEX})$")
MI_PATTERN = re.compile(f"^({ENC_MI_REGEX})$")
IS_PATTERN = re.compile(f"^({ENC_IS_REGEX})$")
MIS_PATTERN = re.compile(f"^({ENC_MIS_REGEX})$")

PHASE_PATTERN = {
    "complex": COMPLEX_PATTERN,
    "i": I_PATTERN,
    "-i": MI_PATTERN,
    "is": IS_PATTERN,
    "-is": MIS_PATTERN,
}

XP_PAULI_START_REGEX = r"\(?[IXZY].*"
SPLIT_PATTERN = re.compile(f"^(.*?)({XP_PAULI_START_REGEX})")

PHASE_REGEX = r"[\-+]?1?[ij]?"
XP_PAULI_REGEX = r"[IXZY]"

INDEX_REGEX = r".*?[0-9].*"
INDEX_INDICATE_PATTERN = re.compile(f"^{INDEX_REGEX}$")

LATEX_REGEX = r".*?_\{[0-9].*"
LATEX_INDICATE_PATTERN = re.compile(f"^{LATEX_REGEX}$")


ENC_INDEX_XZ_REGEX = r"(\((X[0-9]+|Z[0-9]+|X([0-9]+)Z\3)\))+"
ENC_INDEX_XZY_REGEX = r"([XZY][0-9]+)+"
ENC_INDEX_ZX_REGEX = r"(\((X[0-9]+|Z[0-9]+|Z([0-9]+)X\3)\))+"
ENC_INDEX_YZX_REGEX = r"([XZY][0-9]+)+"


INDEX_XZ_PATTERN_CAP = re.compile(r"\((X[0-9]+|Z[0-9]+|X(?:[0-9]+)Z[0-9]+)\)")
INDEX_ZX_PATTERN_CAP = re.compile(r"\((Z[0-9]+|X[0-9]+|Z(?:[0-9]+)X[0-9]+)\)")

INDEX_XZ_PATTERN = re.compile(f"^{ENC_INDEX_XZ_REGEX}$")
INDEX_XZY_PATTERN = re.compile(f"^{ENC_INDEX_XZY_REGEX}$")
INDEX_ZX_PATTERN = re.compile(f"^{ENC_INDEX_ZX_REGEX}$")
INDEX_YZX_PATTERN = re.compile(f"^{ENC_INDEX_YZX_REGEX}$")

TENSOR_INDEX_PATTERN = {
    "XZ": INDEX_XZ_PATTERN,
    "XZY": INDEX_XZY_PATTERN,
    "ZX": INDEX_ZX_PATTERN,
    "YZX": INDEX_YZX_PATTERN,
}

ENC_LATEX_XZ_REGEX = r"(\((X_\{[0-9]+\}|Z_\{[0-9]+\}|X_\{([0-9]+\})Z\3)\))+"
ENC_LATEX_XZY_REGEX = r"([XZY]_\{[0-9]+\})+"
ENC_LATEX_ZX_REGEX = r"(\((X_\{[0-9]+\}|Z_\{[0-9]+\}|Z_\{([0-9]+\})X\3)\))+"
ENC_LATEX_YZX_REGEX = r"([XZY]_\{[0-9]+\})+"


LATEX_XZ_PATTERN_CAP = re.compile(r"\((X_\{[0-9]+\}|Z_\{[0-9]+\}|X(?:_\{[0-9]+\})Z_\{[0-9]+\})\)")
LATEX_ZX_PATTERN_CAP = re.compile(r"\((Z_\{[0-9]+\}|X_\{[0-9]+\}|Z(?:_\{[0-9]+\})X_\{[0-9]+\})\)")

LATEX_XZ_PATTERN = re.compile(f"^{ENC_LATEX_XZ_REGEX}$")
LATEX_XZY_PATTERN = re.compile(f"^{ENC_LATEX_XZY_REGEX}$")
LATEX_ZX_PATTERN = re.compile(f"^{ENC_LATEX_ZX_REGEX}$")
LATEX_YZX_PATTERN = re.compile(f"^{ENC_LATEX_YZX_REGEX}$")

TENSOR_LATEX_PATTERN = {
    "XZ": LATEX_XZ_PATTERN,
    "XZY": LATEX_XZY_PATTERN,
    "ZX": LATEX_ZX_PATTERN,
    "YZX": LATEX_YZX_PATTERN,
}


ENC_PRODUCT_XZ_CAP = r"\(([XZ]|XZ|I)\)"
ENC_PRODUCT_ZX_CAP = r"\(([ZX]|ZX|I)\)"

PRODUCT_XZ_PATTERN_CAP = re.compile(f"{ENC_PRODUCT_XZ_CAP}")
PRODUCT_ZX_PATTERN_CAP = re.compile(f"{ENC_PRODUCT_ZX_CAP}")

ENC_PRODUCT_XZ_REGEX = r"(\(([XZ]|XZ|I)\))+"
ENC_PRODUCT_XZY_REGEX = r"[XZYI]+"
ENC_PRODUCT_ZX_REGEX = r"(\(([ZX]|ZX|I)\))+"
ENC_PRODUCT_YZX_REGEX = r"[XZYI]+"

PRODUCT_XZ_PATTERN = re.compile(f"^{ENC_PRODUCT_XZ_REGEX}$")
PRODUCT_XZY_PATTERN = re.compile(f"^{ENC_PRODUCT_XZY_REGEX}$")
PRODUCT_ZX_PATTERN = re.compile(f"^{ENC_PRODUCT_ZX_REGEX}$")
PRODUCT_YZX_PATTERN = re.compile(f"^{ENC_PRODUCT_YZX_REGEX}$")

TENSOR_PRODUCT_PATTERN = {
    "XZ": PRODUCT_XZ_PATTERN,
    "XZY": PRODUCT_XZY_PATTERN,
    "ZX": PRODUCT_ZX_PATTERN,
    "YZX": PRODUCT_YZX_PATTERN,
}

PRODUCT_SYNTAX = 0
INDEX_SYNTAX = 1
LATEX_SYNTAX = 2
XP_SYMPLECTIC_SYNTAX = 3
DEFAULT_SYNTAX = 0
SYNTAX_TO_TEXT = ["product", "index", "latex", "XP symplectic"]

DEFAULT_QUBIT_ORDER = "right-to-left"
QUBIT_ORDERS = ["right-to-left", "left-to-right"]


def _is_pattern(string, pattern):
    """_summary_"""
    return bool(pattern.search(string))


# -------------------------------------------------------------------------------
# Encoding lists
# -------------------------------------------------------------------------------


def get_phase_encodings() -> List[str]:
    """Returns the availble phase encodings

    Returns:
        encoding: List of available phase encodings

    Examples:
        >>> get_phase_encodings()
        ['w']

    See Also
        get_tensor_encodings, get_pauli_encodings
    """
    return PHASE_ENCODINGS


def get_tensor_encodings() -> List[str]:
    """Returns the available tensor encodings

    Returns:
        encoding: List of available tensor encodings

    Examples:
        >>> get_tensor_encodings()
        ['XP']

    See Also:
        get_phase_encodings, get_pauli_encodings
    """
    return TENSOR_ENCODINGS


def get_xp_pauli_encodings() -> List[str]:
    """Returns the available XPPauli encodings

    Returns:
        encodings : List of available XPPauli encodings

    Example:
        >>> get_xp_pauli_encodings()
        ['wXP']

    See Also:
        get_phase_encodings, get_tensor_encodings

    """
    return XP_PAULI_ENCODINGS


# -------------------------------------------------------------------------------
# Encoding Methods and Conversions
# -------------------------------------------------------------------------------

def split_xp_pauli_enc(encoding: str) -> Tuple[str, str]:
    """Splits the XPPauli encoding into the phase and tensor encodings

    Args:
        encoding: XPPauli encoding

    Raises:
        QiskitError: Encoding not valid

    Returns:
        phase_enc, tensor_enc: phase encoding and tensor encoding

    Exampes:
        >>> encoding = "wXP"
        >>> split_xp_pauli_encoding(encoding)
        ('w', 'XP')

    See Also:
        _split_xp_pauli_encoding
    """
    if encoding not in XP_PAULI_ENCODINGS:
        raise QiskitError(f"Encoding not valid: {encoding}")
    return _split_xp_pauli_enc(encoding)


def _split_xp_pauli_enc(encoding: str) -> Tuple[str, str]:
    """Splits the XPPauli encoding into the phase and tensor encodings

    Args:
        encoding: XPPauli encoding string
    """
    return XP_PAULI_ENCODINGS_SPLIT[encoding]


def get_phase_enc(encoding: str) -> str:
    """Returns the phase encoding part of the XPPauli encoding string

    Args:
        encoding: XPPauli encoding string

    Returns:
        phase_enc: phase encoding
    """
    phase_part, _ = split_xp_pauli_enc(encoding)
    return phase_part


def get_tensor_enc(encoding: str) -> str:
    """Returns the tensor encoding part of the XPPauli encoding string

    Args:
        encoding: XPPauli encoding string

    Returns:
        phase_enc: tensor encoding
    """
    _, tensor_part = split_xp_pauli_enc(encoding)
    return tensor_part


# TODO depending on what encoding formats are decided, these need to be
# implemented or removed.
# pylint: disable=unused-argument
def change_xp_pauli_encoding(
    phase_exp: Any,
    y_count: Union[np.array, int] = 0,
    *,
    input_xp_pauli_encoding: str = INTERNAL_XP_PAULI_ENCODING,
    output_xp_pauli_encoding: str = DEFAULT_EXTERNAL_XP_PAULI_ENCODING,
    same_type=True,
) -> Any:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _change_xp_pauli_encoding(
    phase_exponent: np.ndarray,
    y_count: np.ndarray,
    input_xp_pauli_encoding: str,
    output_xp_pauli_encoding: str,
) -> Any:
    """_summary_"""
    pass


# TODO depending on what encoding formats are decided, these need to be
# implemented or removed.
# pylint: disable=unused-argument
def stand_phase_str(
    phase_str: str, same_type: bool = True, imaginary: str = "i"
) -> Union[np.ndarray, str]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _stand_phase_str(phase_string: np.ndarray, imaginary: str) -> np.ndarray:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def is_scalar(obj: Any, scalar_pair: bool = False):
    """_summary_"""
    pass


# pylint: disable=unused-argument
def squeeze(array_: Any, scalar: bool = False) -> bool:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def is_exp_type(phase_exp: Any, input_phase_encoding: str) -> bool:
    """_summary_"""
    pass


# Conversion function between XPPauli phases as

# complex string: "0-1j", ...
# complex type: 0-1j, ...
# encoded string: "(-i, 1)", ...
# encoded type: (-i,0)(-1,1), ...

# Methods convering between type do not change the encoding scheme if present
# Encoding schemes can be change via the exp2exp method

# This method are all vector enabled and have two forms: method and _method.
# The method versions have all the checks and conversions in them. The raw
# methods assume correct input formats and generally do very little checking
# or conversions.
#
# The same_type parameter is used to allow scalar values to be output
# when scalar values are input

# ----------------------------------------------------------------------
# pylint: disable=unused-argument
def cpxstr2cpx(
    cpx_str: Union[str, np.ndarray, List[str]],
    same_type: bool = True,
) -> Union[np.ndarray, numbers.Complex]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _cpxstr2cpx(cpx_string: np.ndarray) -> np.ndarray:
    """_summary_"""
    pass


# ----------------------------------------------------------------------
# pylint: disable=unused-argument
def cpx2cpxstr(
    cpx: Union[numbers.Complex, np.ndarray], same_type: bool = True, ones: bool = False
) -> Union[str, np.ndarray]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _cpx2cpxstr(cpx: np.ndarray, one_str: str) -> Union[str, np.ndarray]:
    """_summary_"""
    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def exp2cpx(
    phase_exp: Any, input_encoding: str, same_type: bool = True
) -> Union[np.ndarray, numbers.Complex]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _exp2cpx(phase_exp: np.ndarray, input_encoding: str) -> np.ndarray:
    """_summary_"""
    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def cpx2exp(
    cpx: numbers.Complex, output_encoding: str, same_type: bool = True, roundit: bool = True
) -> Union[np.ndarray, Tuple[numbers.Integral, numbers.Integral], numbers.Integral]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _cpx2exp(cpx: numbers.Complex, encoding: str) -> np.ndarray:
    """_summary_"""
    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def expstr2exp(exp_str: Union[np.ndarray, str], same_type: bool = True) -> Any:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _expstr2exp(exp_string, encoding: str) -> np.ndarray:
    """_summary_"""
    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def exp2expstr(
    phase_exp: Any,
    input_encoding: str = DEFAULT_EXTERNAL_XP_PAULI_ENCODING,
    same_type: bool = True,
) -> Union[np.ndarray, str]:
    """Converts encoded phases (exponents) to their string representations

    Note: This method does more than apply str method as the string representations of the
    different encodings have a specific syntaxs.

    Args:
        phase_exp: Phase encosings to convert to string representations
        input_encoding: Encoding of the input phase exponents. Defaults to 
        DEFAULT_EXTERNAL_XP_PAULI_ENCODING.
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Raises:
        QiskitError: Invalid phase exponent encoding

    Returns:
        exp_str: string representations of given phase exponents

    Examples:
        TODO

    See Also:
        _exp2expstr
    """
    if input_encoding not in PHASE_ENCODINGS:
        raise QiskitError(f"Invalid phase exponent encoding: {input_encoding}")

    phase_exp = np.atleast_1d(phase_exp)

    return _exp2expstr(phase_exp, input_encoding)


def _exp2expstr(phase_exp: np.ndarray, encoding: str) -> np.ndarray:
    """Converts encoded phases (exponents) to their string representations

    Note: This method does more than apply str as the string representations of the
    different encodings have a specific syntax.

    Args:
        phase_exp: Phase encosings to convert to string representations
        encoding: Encoding of the input phase exponents

    Raises:
        QiskitError: Invalid phase exponent encoding

    Returns:
        exp_str: The encoding is not supported

    Examples:
        TODO

    See Also:
        exp2expstr
    """
    if encoding == "w":
        return np.array(["(w," + str(item) + ")" for item in phase_exp])
    else:
        raise QiskitError(f"The encoding {encoding} is not supported.")


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def exp2exp(
    phase_exp: Union[np.ndarray, Any],
    input_encoding=INTERNAL_PHASE_ENCODING,
    output_encoding=DEFAULT_EXTERNAL_PHASE_ENCODING,
    same_type: bool = True,
):
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _exp2exp(phase_exp, input_encoding, output_encoding):
    """_summary_"""
    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def cpxstr2expstr(
    cpx_str: Union[np.ndarray, str], encoding: str, same_type: bool = True
) -> Union[np.ndarray, str]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def expstr2cpxstr(
    exp_str: Union[np.ndarray, str], encoding: str, same_type: bool = True, ones: bool = False
) -> Union[str, np.ndarray]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def cpxstr2exp(cpx_str: Union[np.ndarray, str], encoding: str, same_type: bool = True) -> Any:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def exp2cpxstr(
    phase_exp: Any, encoding: str, same_type: bool = True, ones: bool = False
) -> Union[str, np.ndarray]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def expstr2cpx(
    phase_str: Union[np.ndarray, str], encoding: str, same_type: bool = True
) -> Union[np.ndarray, numbers.Complex]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def cpx2expstr(
    cpx: Union[np.ndarray, numbers.Complex], encoding: str, same_type: bool = True
) -> Union[np.ndarray, str]:
    """_summary_"""
    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def str2exp(
    phase_str: Union[np.ndarray, str], encoding: str, same_type: bool = True
) -> Union[np.ndarray, str]:
    """_summary_"""
    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def _str2exp(phase_str: np.ndarray, encoding: str) -> np.ndarray:
    """_summary_"""
    pass


# ---------------------------------------------------------------------
# Label parsing helper functions
# ---------------------------------------------------------------------


# pylint: disable=unused-argument
def split_pauli(pauli_str: str, same_type: bool = True) -> Tuple[str, str]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _split_pauli(pauli_str: str):
    """_summary_"""

    def _split(p_str):
        pass

    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def encode_of_phase_str(phase_str: str, same_type: bool = True) -> Union[np.ndarray, str]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _encode_of_phase_str(phase_str: str) -> np.ndarray:
    """_summary_"""

    # pylint: disable=unused-variable
    def _find_encode(p_str):
        pass

    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def encode_of_tensor_str(
    tensor_str: str, encoded: bool = True, same_type: bool = True
) -> List[Tuple[List, Union[str, int]]]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _encode_of_tensor_str(
    tensor_str: np.ndarray, encoded: bool
) -> List[Tuple[List, Union[str, int]]]:
    """_summary_"""

    # pylint: disable=unused-variable
    def _find_encode(t_str: str, ecode: bool):
        pass

    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def str2symplectic(
    pauli_str: Union[np.ndarray, str],
    qubit_order: str = "right-to-left",
    output_encoding: Optional[str] = INTERNAL_XP_PAULI_ENCODING,
    index_start: int = 0,
    same_type: bool = True,
) -> Tuple[np.ndarray, Union[np.array, Any]]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _str2symplectic(
    pauli_str: np.ndarray,
    qubit_order: str,
    output_encoding: str,
    index_start: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """_summary_"""
    pass


# ----------------------------------------------------------------------
def xp_symplectic2str(
    matrix: np.ndarray,
    phase_exp: Any = None,
    precision:int = None,
    input_encoding: str = INTERNAL_XP_PAULI_ENCODING,
    output_phase_encoding: str = None,
    no_phase=False,
    output_tensor_encoding: str = DEFAULT_EXTERNAL_TENSOR_ENCODING,
    syntax: str = INDEX_SYNTAX,
    qubit_order: str = "right-to-left",
    index_start: int = 0,
    same_type=True,
    index_str="",
) -> Union[np.ndarray, str]:
    """Converts a symplectic matrix and phase to string representations

    Args:
        matrix: Generalized symplectic matrix for XP operator
        phase_exp (optional): Phase exponent(s) for matrix. A value of
            None will lead to unity phases. Defaults to None.
        precision: Precision of XP operator.
        input_encoding (optional): XPPauli encoding of phase relative to
            matrix. Defaults to INTERNAL_XP_PAULI_ENCODING.
        output_phase_encoding (optional): Encoding used to represent phases.
            A value of None will result in complex phases notation. Defaults
            to None.
        no_phase (optional): When set to True, no phase will appear no matter
            what encoding is selected. So the symplectic matrix [1, 1] will produce
            the operator Y in 'XZY' encoding but also (XZ) in the 'XZ' encoding which
            are different operators if phases are considered.
        output_tensor_encoding (optional): Encoding of XPPauli tensor
            (without phase). Defaults to DEFAULT_EXTERNAL_TENSOR_ENCODING.
        syntax (optional): Syntax of pauli tensor. Values are
            PRODUCT_SYNTAX = 0, INDEX_SYNTAX=1, LATEX_SYNTAX=2 and XP_SYMPLECTIC_SYNTAX=3.
            Defaults to INDEX_SYNTAX.
        qubit_order (optional): Order in which qubits are read. options are
            "right-to-left" and "left-to-right". Defaults to "right-to-left".
        index_start (optional): Lowest value for index in index syntax tensors.
            Defaults to 0
        same_type (optional): Scalar/Vector return flag. Defaults to True.
        index_str (optional): String that get inserted between operator and numbers in
            index format. Default is "".

    Raises:
        QiskitError: Unsupport syntax

    Returns:
        xp_pauli_str: XPPauli strings

    Examples:
        TODO
    """
    matrix = np.atleast_2d(matrix)
    num_qubits = matrix.shape[1] >> 1
    matrix[:,0:num_qubits] = np.mod(matrix[:,0:num_qubits], 2)
    matrix[:,num_qubits:] = np.mod(matrix[:,num_qubits:], precision)
    if no_phase:
        phase_str = np.full((matrix.shape[0],), "")
    else:
        if phase_exp is None:
            phase_exp = np.zeros(shape=(matrix.shape[0],), dtype=np.int64)
        else:
            phase_exp = np.atleast_1d(phase_exp)
            phase_exp = np.mod(phase_exp, 2*precision)

        # If multiple phase/tensor encodings are implemented, the conversion
        # needs to go here.

        if output_phase_encoding is None:
            if syntax != XP_SYMPLECTIC_SYNTAX:
                phase_str = exp2expstr(phase_exp, "w")

    tensor_str = []

    _XPENC = ["(I)", "(X)", "(P{zexp})", "(XP{zexp})"]
    _ENC = {"XP": _XPENC}

    if syntax == PRODUCT_SYNTAX:
        for xppauli in matrix:
            tmp_tensor_str = ""
            for index in range(num_qubits):
                tmp_enc = ""
                rep = ""
                if xppauli[index+num_qubits] > 1:
                    rep = str(xppauli[index+num_qubits])
                if xppauli[index+num_qubits] == 0:
                   tmp_enc = _ENC[output_tensor_encoding][xppauli[index]]
                else:
                   tmp_enc = _ENC[output_tensor_encoding][2+xppauli[index]].replace("{zexp}", rep)
                if tmp_enc:
                    if qubit_order == "left-to-right":
                        tmp_tensor_str += tmp_enc
                    else:
                        tmp_tensor_str = tmp_enc + tmp_tensor_str

            tensor_str.append(tmp_tensor_str)
    elif syntax == XP_SYMPLECTIC_SYNTAX:
        for i, xppauli in enumerate(matrix):
            tmp_tensor_str = "XP"+str(precision)
            tmp_tensor_str += "("+str(phase_exp[i])+"|"+str(xppauli[:num_qubits])[1:-1]+"|"+str(xppauli[num_qubits:])[1:-1]+")"

            tensor_str.append(tmp_tensor_str)
    elif syntax in (INDEX_SYNTAX, LATEX_SYNTAX):
        if syntax == LATEX_SYNTAX:
            ind_str_repr = _ind_to_latex_repr
            sup_str_repr = _sup_to_latex_repr
        else:
            ind_str_repr = str
            sup_str_repr = str

        for xppauli in matrix:
            tmp_tensor_str = ""
            for index in range(num_qubits):
                if output_tensor_encoding == "XP":
                    tmp_enc = ""
                    rep = ""
                    if xppauli[index+num_qubits] > 1:
                        rep = sup_str_repr(xppauli[index+num_qubits])
                    if xppauli[index] == 1 and xppauli[index+num_qubits] == 0:
                        tmp_enc = _ENC[output_tensor_encoding][1]
                    elif xppauli[index+num_qubits] > 0:
                        tmp_enc = _ENC[output_tensor_encoding][2+xppauli[index]].replace("{zexp}", rep)
                    
                    if tmp_enc:
                        if qubit_order == "left-to-right":
                            tmp_tensor_str += (tmp_enc + index_str + ind_str_repr(index + index_start))
                        else:
                            tmp_tensor_str = (tmp_enc + index_str + ind_str_repr(index + index_start) + tmp_tensor_str)

            tensor_str.append(tmp_tensor_str)

    else:
        raise QiskitError(f"Unsupported syntax: {syntax}")

    if syntax != XP_SYMPLECTIC_SYNTAX:
        if syntax in (PRODUCT_SYNTAX, INDEX_SYNTAX):
            result = ["XP"+str(precision)+"("+ p_str + t_str+")" for p_str, t_str in zip(phase_str, tensor_str)]
        else:
            result = ["XP"+_ind_to_latex_repr(precision)+"("+ p_str + t_str+")" for p_str, t_str in zip(phase_str, tensor_str)]
    else:
        result = tensor_str
    if matrix.shape[0] == 1 and same_type:
        return result[0]
    else:
        return np.array(result)


# ----------------------------------------------------------------------
# pylint: disable=unused-argument
def str2str(
    pauli_str: Union[np.ndarray, str],
    syntax_output: str,
    phase_encoding_output_string: str = None,
    tensor_encoding_output_string: str = DEFAULT_EXTERNAL_TENSOR_ENCODING,
    qubit_order_input: str = "right-to-left",
    qubit_order_output: str = "right-to-left",
    index_start_input: int = 0,
    index_start_output: int = 0,
) -> Union[np.ndarray, str]:
    """_summary_"""
    pass


# ----------------------------------------------------------------------
# Array(s) to Symplectic
# ----------------------------------------------------------------------


def from_array(
    matrix: Union[List, Tuple, np.ndarray],
    phase_exp: Union[int, List, Tuple, np.ndarray] = None,
    precision: Union[int, List, Tuple, np.ndarray] = None,
    input_pauli_encoding: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert array and phase_exp to the matrix, phase_exp and precision in BaseXPPauli's internal
    XPPauli encoding (xp_pauli_rep.INTERNAL_XP_PAULI_ENCODING)
    Args:
        matrix (_type_): _description_
        phase_exp (_type_): _description_
        precision (_type_): Precision of XP operator
        input_pauli_encoding: input XPPauli encoding
    Returns:
        _type_: _description_
    """
    if input_pauli_encoding is None:
        input_pauli_encoding = DEFAULT_EXTERNAL_XP_PAULI_ENCODING
    if isinstance(matrix, np.ndarray) and matrix.dtype == np.int64:
        matrix_data = matrix
    else:
        matrix_data = np.asarray(matrix, dtype=np.int64)
    matrix_data = np.atleast_2d(matrix_data)
    # TODO
    # if not is_symplectic_matrix_form(matrix_data):
    #     raise QiskitError("Input matrix not a symplectic matrix or symplectic vector")
    if phase_exp is None or (isinstance(phase_exp, numbers.Integral) and phase_exp == 0):
        phase_exp = np.zeros(shape=(matrix_data.shape[0],))
    # TODO may need to implement change_pauli_encoding
    return matrix_data, phase_exp, precision


# pylint: disable=unused-argument
def from_split_array(
    x: Union[List, Tuple, np.ndarray],
    z: Union[List, Tuple, np.ndarray],
    phase_exp: Union[int, List, Tuple, np.ndarray],
    input_pauli_encoding: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """_summary_"""
    pass


# ----------------------------------------------------------------------
# Symplectic to Complex Matrix conversion
# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def to_cpx_matrix(
    matrix: np.ndarray,
    phase_exp: Optional[np.ndarray],
    pauli_encoding: str = INTERNAL_XP_PAULI_ENCODING,
    sparse: bool = False,
) -> Union[np.ndarray, csr_matrix]:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def _to_cpx_matrix(
    matrix: np.ndarray, phase_exp: int, num_qubits: int, sparse: bool = False
) -> Union[np.ndarray, csr_matrix]:
    """_summary_"""
    pass


# ----------------------------------------------------------------------
# Utility functions
# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def indices_to_boolean(indices: Iterable[int], dim: int) -> np.ndarray:
    """_summary_"""
    pass


# ----------------------------------------------------------------------
def _ind_to_latex_repr(index: int) -> str:
    """Adds curly braces and an underscore to an index.

    Args:
        index: An integer

    Returns:
        str: string in LaTeX syntax
    """
    return f"_{{{index}}}"


# ----------------------------------------------------------------------
def _sup_to_latex_repr(superscript: int) -> str:
    """Adds curly braces and a caret to a superscript.

    Args:
        superscript: An integer

    Returns:
        str: string in LaTeX syntax
    """
    return f"^{{{superscript}}}"


# pylint: disable=unused-argument
def boolean_to_indices(booleans: Iterable[bool]) -> np.ndarray:
    """_summary_"""
    pass


# pylint: disable=unused-argument
def scalar_op2symplectic(
    op: ScalarOp, output_encoding: str = DEFAULT_EXTERNAL_PHASE_ENCODING
) -> Tuple[np.ndarray, Union[np.array, Any]]:
    """_summary_"""
    pass


# ----------------------------------------------------------------------


# pylint: disable=unused-argument
def gate2symplectic(
    gate: Gate, encoding: str = INTERNAL_XP_PAULI_ENCODING
) -> Tuple[np.ndarray, Union[np.array, Any]]:
    """_summary_"""
    pass


# ----------------------------------------------------------------------
