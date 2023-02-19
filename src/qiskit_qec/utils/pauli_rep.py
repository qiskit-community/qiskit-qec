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
"""
N-qubit Pauli Representation Encodings and Conversion Module
"""

# pylint: disable=invalid-name,anomalous-backslash-in-string

import numbers
import re
from typing import Any, Iterable, List, Optional, Tuple, Union

import numpy as np
from qiskit.circuit import Gate
from qiskit.circuit.library.generalized_gates import PauliGate
from qiskit.circuit.library.standard_gates import IGate, XGate, YGate, ZGate
from qiskit.exceptions import QiskitError
from qiskit.quantum_info.operators.scalar_op import ScalarOp
from scipy.sparse import csr_matrix
from qiskit_qec.linear.symplectic import count_num_y, is_symplectic_matrix_form


# -------------------------------------------------------------------------------
# Module Variables/States
# -------------------------------------------------------------------------------

# Set the interal encodings
# The internal encoding cannot be changed by changing this constant
# These constants are for reference only and do not change the behavior of
# the Pauli methods. See [ref] for details on the different encodings
# TODO: Include ref for above.

INTERNAL_TENSOR_ENCODING = "ZX"
INTERNAL_PHASE_ENCODING = "-i"
INTERNAL_PAULI_ENCODING = INTERNAL_PHASE_ENCODING + INTERNAL_TENSOR_ENCODING

DEFAULT_EXTERNAL_TENSOR_ENCODING = "YZX"
DEFAULT_EXTERNAL_PHASE_ENCODING = "-i"
DEFAULT_EXTERNAL_PAULI_ENCODING = DEFAULT_EXTERNAL_PHASE_ENCODING + DEFAULT_EXTERNAL_TENSOR_ENCODING

# Set the external encodings
# The external encodings may be changed via the phase_rep_format,
# symp_rep_format, or pauli_rep_format methods. Any method changing
# these values must make sure to update the others
# Tensor encodings are: 'XZ', 'XZY', 'ZX', 'YZX'
# Phase encodings are: 'i', '-i', 'is', '-is'
# See [ref] for details on the different encodings
# TODO: Include ref for above.

PHASE_ENCODINGS = ["i", "-i", "is", "-is"]
PHASE_ENCODINGS_IMI = ["i", "-i"]
PHASE_ENCODINGS_ISMIS = ["is", "-is"]
TENSOR_ENCODINGS = ["XZ", "XZY", "ZX", "YZX"]
Y_TENSOR_ENCODINGS = ["XZY", "YZX"]
PAULI_ENCODINGS = [
    "iXZ",
    "iXZY",
    "iZX",
    "iYZX",
    "-iXZ",
    "-iXZY",
    "-iZX",
    "-iYZX",
    "isXZ",
    "isXZY",
    "isZX",
    "isYZX",
    "-isXZ",
    "-isXZY",
    "-isZX",
    "-isYZX",
]
PAULI_ENCODINGS_SPLIT = {
    "iXZ": ("i", "XZ"),
    "iXZY": ("i", "XZY"),
    "iZX": ("i", "ZX"),
    "iYZX": ("i", "YZX"),
    "-iXZ": ("-i", "XZ"),
    "-iXZY": ("-i", "XZY"),
    "-iZX": ("-i", "ZX"),
    "-iYZX": ("-i", "YZX"),
    "isXZ": ("is", "XZ"),
    "isXZY": ("is", "XZY"),
    "isZX": ("is", "ZX"),
    "isYZX": ("is", "YZX"),
    "-isXZ": ("-is", "XZ"),
    "-isXZY": ("-is", "XZY"),
    "-isZX": ("-is", "ZX"),
    "-isYZX": ("-is", "YZX"),
}

# Different string syntax formats are available. The they are "product" syntax and
# "index" syntax. "product" syntax represents a Pauli operator of the form
# :math: $p * T_1 \otimes T_2 \otimes ... \otimes T_n$ as :math" $pT1T2T3...Tn$. See the
# following exmaples:
#
# -iX \otimes Y \otimes Z -> -iXYZ
# X \otimes Y \otimes Z \otimes I \otimes I \otimes I  -> XYZII
#
# The index syntax only represents the non identity Paulis. Index syntax specifically
# indiciates the index that the Pauli's are acting on. Following Qiskit's current internal
# indexing:
#
# -iX \otimes Y \otimes Z -> -iZ0Y1X2
# X \otimes Y \otimes Z \otimes I \otimes I \otimes I  -> Z3Y4X5

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

PAULI_START_REGEX = r"\(?[IXZY].*"
SPLIT_PATTERN = re.compile(f"^(.*?)({PAULI_START_REGEX})")

PHASE_REGEX = r"[\-+]?1?[ij]?"
PAULI_REGEX = r"[IXZY]"

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
DEFAULT_SYNTAX = 0
SYNTAX_TO_TEXT = ["product", "index"]

DEFAULT_QUBIT_ORDER = "right-to-left"
QUBIT_ORDERS = ["right-to-left", "left-to-right"]


def _is_pattern(string, pattern):
    return bool(pattern.search(string))


# -------------------------------------------------------------------------------
# Encoding lists
# -------------------------------------------------------------------------------


def get_phase_encodings() -> List[str]:
    """Returns the availble phase encodings

    Returns:
        encoding: List of available phase encodings

    Examples:
        >>> phase_encodings()
        ['i', '-i', 'is', '-is']

    See Also
        get_tensor_encodings, get_pauli_encldings
    """
    return PHASE_ENCODINGS


def get_tensor_encodings() -> List[str]:
    """Returns the available tensor encodings

    Returns:
        encoding: List of available tensor encodings

    Examples:
        >>> tensor_encodings()
        ['XZ', 'XZY', 'ZX', 'YZX']

    See Also:
        get_phase_encodings, get_pauli_encodings
    """
    return TENSOR_ENCODINGS


def get_pauli_encodings() -> List[str]:
    """Returns the available Pauli encodings

    Returns:
        encodings : List of available Pauli encodings

    Example:
        >>> pauli_encodings()
        ['iXZ',
         'iXZY',
         'iZX',
         'iYZX',
         '-iXZ',
         '-iXZY',
         '-iZX',
         '-iYZX',
         'isXZ',
         'isXZY',
         'isZX',
         'isYZX',
         '-isXZ',
         '-isXZY',
         '-isZX',
         '-isYZX']

    See Also:
        get_phase_encodings, get_tensor_encodings

    """
    return PAULI_ENCODINGS


# -------------------------------------------------------------------------------
# Encoding Methods and Conversions
# -------------------------------------------------------------------------------


def split_pauli_enc(encoding: str) -> Tuple[str, str]:
    """Splits the Pauli encoding into the phase and tensor encodings

    Args:
        encoding: Pauli encpoding

    Raises:
        QiskitError: Encoding not valid

    Returns:
        phase_enc, tensor_enc: phase encoding and tensor encoding

    Exampes:
        >>> encoding = "iXZ'
        >>> split_pauli_encoding(encoding)
        ('i', 'XZ')

    See Also:
        _split_pauli_encoding
    """
    if encoding not in PAULI_ENCODINGS:
        raise QiskitError(f"Encoding not valid: {encoding}")
    return _split_pauli_enc(encoding)


def _split_pauli_enc(encoding: str) -> Tuple[str, str]:
    """Splits the Pauli encoding into the phase and tensor encodings

    Args:
        encoding: Pauli encoding string
    """
    return PAULI_ENCODINGS_SPLIT[encoding]


def get_phase_enc(encoding: str) -> str:
    """Returns the phase encodeing part of the Pauli encoding string

    Args:
        encoding: Pauli encoding string

    Returns:
        phase_enc: phase encoding
    """
    phase_part, _ = split_pauli_enc(encoding)
    return phase_part


def get_tensor_enc(encoding: str) -> str:
    """Returns the tensor encodeing part of the Pauli encoding string

    Args:
        encoding: Pauli encoding string

    Returns:
        phase_enc: tensor encoding
    """
    _, tensor_part = split_pauli_enc(encoding)
    return tensor_part


def change_pauli_encoding(
    phase_exp: Any,
    y_count: Union[np.array, int] = 0,
    *,
    input_pauli_encoding: str = INTERNAL_PAULI_ENCODING,
    output_pauli_encoding: str = DEFAULT_EXTERNAL_PAULI_ENCODING,
    same_type=True,
) -> Any:
    """Converts the input phase exponent of a Pauli operator encoded in the input_pauli_encoding
    into the phase exponent of the same Pauli operator encoded in the output_pauli_encoding.

    This method should not be confused with the exp2exp method that changes the encoding of a
    phase exponent independently of a Pauli operator - i.e. assumes that y_count = 0.

    So for example a the Pauli operator "(-i,2)XYIX" ( = -XYIX ) when encoded using the "-iXZY"
    encoding is equal to the Pauli operator "(-i,1)(-1,0)(X)(XZ)(I)(X)" when encoded using the
    "-isXZ" encoding since:

    -XYIX = (-i)^2 XYIX
          = -i (X)(XZ)(I)(X)
          = (-i)^1(-1)^0 (X)(XZ)(I)(X)

    Args:
        phase_exp : phase exponent(s) to convert relative to a Pauli(s) with y_count Y's
        y_count : number of Ys (XZ,ZX) factors in associated Pauli's
        input_pauli_encoding: Pauli encoding of input
        output_pauli_encoding: Pauli encoding of output
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Raises:
        QiskitError: Unsupported Pauli encoding
        QiskitError: The phase exponent must be an element in Z_4 or a pair of elements
            in GF(2) x GF(2)

    Examples:
        >>> change_pauli_encoding(2,
                                  y_count=1,
                                  input_pauli_encoding="-iXZY",
                                  output_pauli_encoding="-iXZ")
        1

        >>> change_pauli_encoding(np.array([2,1]),
                                  y_count=np.array([1,3]),
                                  input_pauli_encoding="-iXZY",
                                  output_pauli_encoding="-isXZ")
        array([[1, 0],
               [0, 1]])

        >>> change_pauli_encoding(numpy.array([1,1]),
                                  y_count=1,
                                  input_pauli_encoding="isXZY",
                                  output_pauli_encoding="iXZ")
        0

        See Also:
            _change_pauli_encoding, exp2exp, _exp2exp
    """

    if input_pauli_encoding not in PAULI_ENCODINGS:
        raise QiskitError(f"Unsupported Pauli encoding {input_pauli_encoding}")
    if output_pauli_encoding not in PAULI_ENCODINGS:
        raise QiskitError(f"Unsupported Pauli endoding {output_pauli_encoding}")
    if not isinstance(phase_exp, (int, np.integer)) and not isinstance(
        phase_exp, (tuple, np.ndarray, list)
    ):
        raise QiskitError(
            "The phase exponent must be an element in Z_4 or a pair of elements in GF(2) x GF(2)"
        )

    input_phase_encoding, _ = _split_pauli_enc(input_pauli_encoding)

    scalar_pair = input_phase_encoding in PHASE_ENCODINGS_ISMIS
    scalar = is_scalar(phase_exp, scalar_pair=scalar_pair) and same_type

    if scalar_pair:
        phase_exp = np.atleast_2d(phase_exp)
    else:
        phase_exp = np.atleast_1d(phase_exp)

    y_count = np.atleast_1d(y_count)

    if scalar:
        return squeeze(
            _change_pauli_encoding(phase_exp, y_count, input_pauli_encoding, output_pauli_encoding),
            scalar=scalar,
        )

    return _change_pauli_encoding(phase_exp, y_count, input_pauli_encoding, output_pauli_encoding)


def _change_pauli_encoding(
    phase_exponent: np.ndarray,
    y_count: np.ndarray,
    input_pauli_encoding: str,
    output_pauli_encoding: str,
) -> Any:
    """Converts the input phase exponent of a Pauli operator encoded in the input_pauli_encoding
    into the phase exponent of the same Pauli operator encoded in the output_pauli_encoding.

    This method should not be confused with the exp2exp method that changes the encoding of a
    phase exponent independently of a Pauli operator - i.e. assumes that y_count = 0.

    So for example a the Pauli operator "(-i,2)XYIX" ( = -XYIX ) when encoded using the "-iXZY"
    encoding is equal to the Pauli operator "(-i,1)(-1,0)(X)(XZ)(I)(X)" when encoded using the
    "-isXZ" encoding since:

    -XYIX = (-i)^2 XYIX
          = -i (X)(XZ)(I)(X)
          = (-i)^1(-1)^0 (X)(XZ)(I)(X)

    Args:
        phase_exponent: phase exponent(s) to convert
        y_count: number of Y (XZ,ZX) factors in Pauli's
        input_pauli_encoding: Pauli encoding of input
        output_pauli_encoding: Pauli encoding of output

    Examples:
        >>> _change_pauli_encoding(np.array([2]),
                                   y_count=np.array([1]),
                                   input_pauli_encoding="-iXZY",
                                   output_pauli_encoding="-iXZ")
        array([1])

        >>> _change_pauli_encoding(np.array([[0,1],[1,1]]),
                                   y_count=np.array([1,3]),
                                   input_pauli_encoding="-isXZY",
                                   output_pauli_encoding="-iXZ")
        array([1, 0])

    See Also:
        change_pauli_encoding, exp2exp. _exp2exp

    """
    # phases change with changing symplectic formats via a multiple of i. This
    # multiple is given by the converter table: S_1->S_2 has multiplier i^converter[S_1][S_2]

    converter = {
        "XZ": {"ZX": 2, "XZ": 0, "XZY": 1, "YZX": 1},
        "ZX": {"ZX": 0, "XZ": 2, "XZY": 3, "YZX": 3},
        "XZY": {"XZY": 0, "YZX": 0, "XZ": 3, "ZX": 1},
        "YZX": {"YZX": 0, "XZY": 0, "XZ": 3, "ZX": 1},
    }

    input_phase_encoding, input_tensor_encoding = _split_pauli_enc(input_pauli_encoding)
    output_phase_encoding, output_tensor_encoding = _split_pauli_enc(output_pauli_encoding)

    multiplier = converter[input_tensor_encoding][output_tensor_encoding]
    phase_exponent = exp2exp(phase_exponent, input_phase_encoding, output_phase_encoding)

    # modify phase exponents to align with the Y count of the Paulis

    def _cal_phase(exp, marker):
        if marker < 2:
            return (marker, exp[1])
        else:
            return (marker % 2, (exp[1] + 1) % 2)

    if output_phase_encoding == "i":
        phase_exponent = np.mod(phase_exponent + 3 * multiplier * y_count, 4)
    elif output_phase_encoding == "-i":
        phase_exponent = np.mod(phase_exponent + multiplier * y_count, 4)
    elif output_phase_encoding == "is":
        res = np.mod(phase_exponent.T[0] + multiplier * y_count, 4)
        phase_exponent = np.asarray(
            [_cal_phase(exp, marker) for exp, marker in zip(phase_exponent, res)]
        )
    else:
        # multiplier = (4 - multiplier) % 4
        res = np.mod(phase_exponent.T[0] + multiplier * y_count, 4)
        phase_exponent = np.asarray(
            [_cal_phase(exp, marker) for exp, marker in zip(phase_exponent, res)]
        )
    return phase_exponent


def stand_phase_str(
    phase_str: str, same_type: bool = True, imaginary: str = "i"
) -> Union[np.ndarray, str]:
    """Standardize the phase string

    The string representing a phase may use different
    representations for the complex number i such
    as "i', "1i", "1j", "j", etc. Note that the string may be
    encoded or not encoded. Likewise other aspects of the
    string have different variations. These method
    converts the input string representing a phase into a
    new string that represents the same phase but uses a
    standard set of symbols and conventions that make other
    utility functions simpler.

    Args:
        phase_str: string representing a phase of a Pauli operator
        same_type (optional): is same type. Default is True
        imaginary: determines what symbol is used to represent the complex number sqrt(-1).
            Defaults to 'i'

    Examples:
        >>> phase_str = "1j"
        >>> stand_phase_str(phase_str)
        "i"

        >>> phase_str = "(-1j,0)(-1,1)"
        >>> stand_phase_str(phase_str)
        "(-i,0)(-1,1)

    Returns:
        out: standardized phase string(s)

    See Also:
        _stand_phase_str
    """
    scalar = is_scalar(phase_str) and same_type
    phase_str = np.atleast_1d(phase_str)

    if scalar:
        return squeeze(_stand_phase_str(phase_str, imaginary), scalar=scalar)

    return _stand_phase_str(phase_str, imaginary)


def _stand_phase_str(phase_string: np.ndarray, imaginary: str) -> np.ndarray:
    """Standardize the phase string

    The string representing a phase may use different
    representations for the complex number i such
    as "i', "1i", "1j", "j", etc. Note that the sting may be
    encoded or not encoded. Likewise other aspects of the
    string have different variations. These method
    converts the input string representing a phase into a
    new string that represents the same phase but uses a
    standard set of symbols and conventions that make other
    utility functions simpler.

    This method does not make any checks.

    Args:
        phase_str: numpy array of string(s) representing a phase of a Pauli operator(s)
        imaginary: determines what symbol is used to represent the complex number sqrt(-1).

    Examples:
        >>> phase_str = np.array(["1j"])
        >>> _stand_phase_str(phase_str)
        array(["i"])

        >>> phase_str = np.array(["(-1j,0)(-1,1)"])
        >>> stand_phase_str(phase_str)
        array(["(-i,0)(-1,1)"])

    Returns:
        out:Standardized phase string(s)

    See Also:
        stand_phase_str
    """
    res = []
    for string in phase_string:
        ans = re.sub("\+?1?[ij]", imaginary, string)
        ans = re.sub("\+?1", "1", ans)
        ans = re.sub("0-", "-", ans)
        res.append(ans)
    return np.array(res)


def is_scalar(obj: Any, scalar_pair: bool = False):
    """Deteremines if an obj is a scalar or not (i.e. not an array of some sort).

    This method is not equivalent to checking if the object is not iterable.
    For this method a string is a scalar even though a string is iterable.

    Args:
        obj: object to test if it is a scalar (i.e. not iterable)
        scalar_pair: If scalar elements are pairs then set scalar_pair to True. Default is False

    Returns:
        out: True if obj is a scalar, False if it is iterable

    Examples:
        >>> a=numpy.array([1,2,3])
        >>> is_scalar(a)
        False

        >>> a=0-1j
        >>> is_scalar(a)
        True

        >>> a="xxxyxyii"
        >>> is_scalar(a)
        True

        >>> a=[1,0]
        >>> is_scalar(a, scalar_pair=True)
        True
    """
    if isinstance(obj, str):
        return True
    else:
        if isinstance(obj, Iterable):
            if scalar_pair:
                return not isinstance(obj[0], Iterable)
            else:
                return False
        else:
            return True


def squeeze(array_: Any, scalar: bool = False) -> bool:
    """Squeeze an numpy array with the option to convert a resulting
    0d array to a scalar (scalar=True)

    Args:
        array_ : Object to see if it is a scalar
        scalar (optional): Will return the scalar if True. Defaults to False.

    Returns:
        bool: a squeezed object

    Examples:
        >>> squeeze(numpy.array([1,2,3]))
        array([1,2,3])

        >>> squeeze(numpy.array(1), scalar=True)
        1

        >>> squeeze(numpy.array(1))
        array(1)

        >>> squeeze(numpy.array([[1,2,3]]))
        array([1,2,3])

        >>> squeeze(numpy.array([[[[1,2,3,4,5]]]]))
        array([1,2,3,4,5])
    """
    array_ = np.squeeze(np.array(array_, dtype=object))
    if array_.shape == () and scalar is True:
        return array_.item()
    else:
        return array_


def is_exp_type(phase_exp: Any, input_phase_encoding: str) -> bool:
    """Determines if a given input is an encoded phase exponent of a Pauli operator

    Care must be taken when phase_exponents are pairs in order to distinguish
    between [0,1] being a single exponent or two exponents. For pairs, phase_exponent
    should be [[0,1]]

    Args:
        phase_exp: input to check is an encoded phase exponent
        input_phase_encoding: phase encoding to check against

    Raises:
        QiskitError: Invalid phase exponent encoding

    Returns:
        out: True if input is a encoded phase

    Examples:
        >>> a = np.array([1])
        >>> is_exp_type(a,"is")
        False

        >>> is_exp_type(a,"i")
        True

    """
    if input_phase_encoding not in PHASE_ENCODINGS:
        raise QiskitError(f"Invalid phase exponent encoding {input_phase_encoding}")
    if isinstance(phase_exp, np.ndarray):
        phase_exp = phase_exp.tolist()
    if not isinstance(phase_exp, list):
        phase_exp = [phase_exp]
    if input_phase_encoding in ["i", "-i"]:
        return all((x in [0, 1, 2, 3]) for x in phase_exp)
    if input_phase_encoding in ["is", "-is"]:
        return all((item in [[0, 0], [0, 1], [1, 0], [1, 1]]) for item in phase_exp)

    return False


# Conversion function between Pauli phases as

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
def cpxstr2cpx(
    cpx_str: Union[str, np.ndarray, List[str]],
    same_type: bool = True,
) -> Union[np.ndarray, numbers.Complex]:
    """Converts the complex number defined by the input string(s) to a complex number

    Args:
        cpx_str: Complex number string(s)
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Returns:
        cpx: complex number from input string of complex number

    Examples:
        >>> cpx_str = "+1j"
        >>> cpxstr2cpx(cpx_str)
        1j

        >>> cpx_str = "-1i"
        >>> cpxstr2cpx(cpx_str)
        (-0-1j)

    See Also:
        Range of conversions between encoded phases... [ref]
    """
    scalar = is_scalar(cpx_str) and same_type
    cpx_str = np.atleast_1d(cpx_str)

    cpx_str = stand_phase_str(cpx_str, same_type=False)

    if scalar:
        return squeeze(_cpxstr2cpx(cpx_str), scalar=scalar)
    else:
        return _cpxstr2cpx(cpx_str)


def _cpxstr2cpx(cpx_string: np.ndarray) -> np.ndarray:
    """Converts the complex number defined by the input string(s) to a complex number

    Complex strings must be in standard representation. See stand_phase_str()

    Args:
        cpx_string : Array of complex number string(s)

    Returns:
        cpx: complex number from input string of complex number

    Raises:
        QiskitError: Not all strings where of strings of complex phases.

    Examples:
        >>> cpx_str = numpy.array(["i"])
        >>> _cpxstr2cpx(cpx_str)
        array([1j])

        >>> cpx_str = numpy.array(["-i"])
        >>> _cpxstr2cpx(cpx_str)
        array([-0.-1.j])

    See Also:
        Range of conversions between encoded phases... [ref]
    """
    CONV_ = {
        "": 1,
        "1": 1,
        "-": -1,
        "i": 0 + 1j,
        "-i": 0 - 1j,
        "-1": -1,
        "()": 1,
        "(1)": 1,
        "(i)": 1j,
        "(-i)": 0 - 1j,
        "(-1)": -1,
    }
    try:
        return np.array([CONV_[item] for item in cpx_string])
    except KeyError as not_all_complex:
        raise QiskitError(
            "Not all strings where of strings of complex phases. \
            Phases must be in standard representation."
        ) from not_all_complex


# ----------------------------------------------------------------------
def cpx2cpxstr(
    cpx: Union[numbers.Complex, np.ndarray], same_type: bool = True, ones: bool = False
) -> Union[str, np.ndarray]:
    """Converts the complex number(s) to strings

    Args:
        cpx: Complex number(s)
        same_type (optional): Scalar/Vector return flag. Defaults to True.
        ones (optional): Defines whether to present ones as "1" (True) or
            the empyty string (False). Defaults to False.

    Returns:
        cpx_str: complex number(s) as string(s)

    Examples:
        >>> cpx = 1j
        >>> cpx2cpxstr(cpx)
        'i'

        >>> cpx = [1j,+1j, -1, 1j, (0-1j)]
        >>> cpx2cpxstr(cpx)
        array(['i', 'i', '-', 'i', '-i'], dtype='<U2')

        >>> cpx = [1j,+1j, -1, 1j, (0-1j)]
        >>> cpx2cpxstr(cpx, ones=True)
        array(['i', 'i', '-1', 'i', '-i'], dtype='<U2')

        >>> cpx = 1j
        >>> cpx2cpxstr(cpx, same_type=False)
        array(['i'], dtype='<U1')

    See Also:
        _cpx2cpxstr
    """
    scalar = is_scalar(cpx) and same_type
    cpx = np.atleast_1d(cpx)

    if ones:
        one_str = "1"
    else:
        one_str = ""

    if scalar:
        return squeeze(_cpx2cpxstr(cpx.round(), one_str=one_str), scalar=scalar)
    else:
        return _cpx2cpxstr(cpx.round(), one_str=one_str)


def _cpx2cpxstr(cpx: np.ndarray, one_str: str) -> Union[str, np.ndarray]:
    """Converts the complex number(s) to strings for phases of the Pauli group.

    Note: This is not a general complex number to string method which
    can be done with the use of str. This method is used to create string
    representations of Pauli operators.

    Args:
        cpx: Array of complex number(s)
        one_str: String to use for 1s

    Returns:
        cpx_str: Array of complex strings

    Examples:
        >>> cpx = numpy.array([1])
        >>> _cpx2cpxstr(cpx, one_str='+')
        array(['+'], dtype='<U1')

        >>> cpx = np.array([1, -1j])
        >>> _cpx2cpxstr(cpx, one_str='')
        array(['', '-i'], dtype='<U2')

    See Also:
        cpx2cpxstr
    """

    cpx_str = [str(item) for item in cpx]
    cpx_str = stand_phase_str(cpx_str)

    _ENC = {1: one_str, 1j: "i", -1: "-" + one_str, 0 - 1j: "-i"}
    return np.asarray([_ENC[item] for item in cpx])


# ----------------------------------------------------------------------


def exp2cpx(
    phase_exp: Any, input_encoding: str, same_type: bool = True
) -> Union[np.ndarray, numbers.Complex]:
    """Converts a encoded phase(s) to a complex number(s)

    Args:
        phase_exp: Encoded phase
        input_encoding: Encoding used to encode the phase
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Raises:
        QiskitError: Invalid phase exponent encoding
        QiskitError: Phase exponent not encoded using given encoding

    Returns:
        cpx: encoded complex phase(s) as complex number(s)

    Examples:
        >>> phase_exp = [[0,1]]
        >>> exp2cpx(phase_exp, '-is')
        array([-1.+0.j])

        >>> phase_exp=[0,1]
        >>> exp2cpx(phase_exp, '-i')
        array([ 1.+0.j, -0.-1.j])

        >>> phase_exp=numpy.array([[0,0],[1,0]])
        >>> exp2cpx(phase_exp, 'is')
        array([1.+0.j, 0.+1.j])

    See Also:
        _exp2cpx, ...
    """
    if input_encoding not in PHASE_ENCODINGS:
        raise QiskitError(f"Invalid phase exponent encoding {input_encoding}")

    scalar_pair = input_encoding in PHASE_ENCODINGS_ISMIS
    scalar = is_scalar(phase_exp, scalar_pair=scalar_pair) and same_type

    if scalar_pair:
        phase_exp = np.atleast_2d(phase_exp)
    else:
        phase_exp = np.atleast_1d(phase_exp)

    if not is_exp_type(phase_exp, input_encoding):
        raise QiskitError(f"{phase_exp} not encoded using '{input_encoding}' encoding")

    if scalar:
        return squeeze(_exp2cpx(phase_exp, input_encoding), scalar=scalar)
    else:
        return _exp2cpx(phase_exp, input_encoding)


def _exp2cpx(phase_exp: np.ndarray, input_encoding: str) -> np.ndarray:
    """Converts a encoded phase(s) to a complex number(s)

    Args:
        phase_exp: Encoded phase
        encoding: Encoding used to encode the phase

    Raises:
        QiskitError: Encoding is not supported

    Returns:
        cpx: array of encoded complex phase(s) as complex number(s)

    Examples:
        >>> phase_exp = numpy.array([[0,1]])
        >>> _exp2cpx(phase_exp, '-is')
        array([-1.+0.j])

        >>> phase_exp=numpy.array([0,1])
        >>> _exp2cpx(phase_exp, '-i')
        array([ 1.+0.j, -0.-1.j])

        >>> phase_exp=numpy.array([[0,0],[1,0]])
        >>> _exp2cpx(phase_exp, 'is')
        array([1.+0.j, 0.+1.j])

    See Also:
        exp2cpx, ...
    """
    if input_encoding == "i":
        return (1j) ** phase_exp
    if input_encoding == "-i":
        return (-1j) ** phase_exp
    if input_encoding == "is":
        trans = phase_exp.T
        return np.multiply(1j ** trans[0], (-1) ** trans[1])
    if input_encoding == "-is":
        trans = phase_exp.T
        return ((-1j) ** trans[0]) * (-1) ** trans[1]
    else:
        raise QiskitError(f"{input_encoding} encoding is not supported.")


# ----------------------------------------------------------------------


def cpx2exp(
    cpx: numbers.Complex, output_encoding: str, same_type: bool = True, roundit: bool = True
) -> Union[np.ndarray, Tuple[numbers.Integral, numbers.Integral], numbers.Integral]:
    """Converts complex phases to encoded exponents for Pauli group phases

    Args:
        cpx: complex pauli phases to encode
        output_encoding: phase encoding to use
        same_type (optional): Scalar/Vector return flag. Defaults to True.
        roundit (optional): round the complex numbers if desired befor encoding. Defaults to True.

    Raises:
        QiskitError: _description_

    Returns:
        encoded_phase : encoded phases

    Examples:
        >>> cpx = 1j
        >>> encoding = '-i'
        >>> cpx2exp(cpx, encoding)
        3

        >>> cpx = 1j
        >>> encoding = '-is'
        >>> cpx2exp(cpx, encoding)
        array([1, 1])

    See Also:
        _cpx2exp
    """
    if output_encoding not in PHASE_ENCODINGS:
        raise QiskitError(f"Invalid phase exponent encoding: {output_encoding}")

    scalar = is_scalar(cpx) and same_type
    cpx = np.atleast_1d(cpx)

    if roundit:
        cpx = cpx.round()

    if scalar:
        return squeeze(_cpx2exp(cpx, output_encoding), scalar=scalar)
    else:
        return _cpx2exp(cpx, output_encoding)


def _cpx2exp(cpx: numbers.Complex, encoding: str) -> np.ndarray:
    """Converts complex phases to encoded exponents for Pauli group phases

    Args:
        cpx: array pf complex pauli phases to encode
        encoding: phase encoding to use

    Raises:
        QiskitError: Complex type phases must be complex numbers in ``[1, -1j, -1, 1j]

    Returns:
        encoded_phase : array of encoded phases

    Examples:
        >>> cpx = numpy.array([1j])
        >>> encoding = '-i'
        >>> _cpx2exp(cpx, encoding)
        array([3])

        >>> cpx = 1j
        >>> encoding = '-is'
        >>> c_px2exp(cpx, encoding)
        array([1, 1])

    """
    _ENC = {
        "i": {1: 0, 1j: 1, -1: 2, 0 - 1j: 3},
        "-i": {1: 0, 0 - 1j: 1, -1: 2, 1j: 3},
        "is": {1: (0, 0), -1: (0, 1), 1j: (1, 0), 0 - 1j: (1, 1)},
        "-is": {1: (0, 0), -1: (0, 1), 0 - 1j: (1, 0), 1j: (1, 1)},
    }
    try:
        return np.array([_ENC[encoding][i] for i in cpx])
    except Exception as exception:
        raise QiskitError(
            "Complex type phases must be complex numbers in ``[1, -1j, -1, 1j]"
        ) from exception


# ----------------------------------------------------------------------


def expstr2exp(exp_str: Union[np.ndarray, str], same_type: bool = True) -> Any:
    """Converts string representations of phase encoded phases to encoded phases

    Args:
        exp_str: string representations of phase exponents. Each string representation
            must be encoded with the same encoding.
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Raises:
        QiskitError: String encoding is not implemented.

    Returns:
        phase_exp: encoded phases

    Examples:
        >>> exp_str = '(-i,1)(-1,1)'
        >>> expstr2exp(exp_str)
        array([1, 1], dtype=int8)

        >>> exp_str = '(i,3)'
        >>> expstr2exp(exp_str)
        3

        >>> exp_str = np.array(['(-i,1)', '(-i,2)', '(-i,0)'])
        >>> expstr2exp(exp_str)
        array([1, 2, 0], dtype=int8)

    See Also:
        _expstr2exp
    """
    scalar = is_scalar(exp_str) and same_type
    exp_str = np.atleast_1d(exp_str)
    exp_str = stand_phase_str(exp_str, same_type=False)
    # It is assumed that all phase strings are of the same encoding/representation
    encoding = encode_of_phase_str(exp_str[0])

    if encoding in PHASE_ENCODINGS:
        if scalar:
            return squeeze(_expstr2exp(exp_str, encoding), scalar=scalar)
        return _expstr2exp(exp_str, encoding)

    raise QiskitError(f"String encoding {encoding} is not implemented.")


def _expstr2exp(exp_string, encoding: str) -> np.ndarray:
    """Converts string representations of phase encoded phases to encoded phases

    Args:
        exp_str: string representations of phase exponents. Each string representation
            must be encoded with the same encoding.
        encoding: encoding used to encode phases

    Returns:
        phase_exp: encoded phases

    Examples:
        >>> exp_str = numpy.array(['(-i,1)(-1,1)'])
        >>> _expstr2exp(exp_str, '-is')
        array([1, 1], dtype=int8)

        >>> exp_str = numpy.array(['(i,3)'])
        >>> _expstr2exp(exp_str, 'i')
        3

        >>> exp_str = np.array(['(-i,1)', '(-i,2)', '(-i,0)'])
        >>> _expstr2exp(exp_str, '-i')
        array([1, 2, 0], dtype=int8)

    See Also:
        expstr2exp
    """
    if encoding in PHASE_ENCODINGS_IMI:
        result = np.zeros(shape=(exp_string.shape[0],), dtype=np.int8)
        for i in range(exp_string.shape[0]):
            val = re.findall(CAP_PHASE_PATTERN[encoding], exp_string[i])
            result[i] = val[0]
    elif encoding in PHASE_ENCODINGS_ISMIS:
        result = np.zeros(shape=(exp_string.shape[0], 2), dtype=np.int8)
        for i in range(exp_string.shape[0]):
            val = re.findall(CAP_PHASE_PATTERN[encoding], exp_string[i])
            result[i][0] = val[0][0]
            result[i][1] = val[0][1]
    return result


# ----------------------------------------------------------------------


def exp2expstr(
    phase_exp: Any,
    input_encoding: str,
    same_type: bool = True,
) -> Union[np.ndarray, str]:
    """Converts encoded phases (exponents) to their string representations

    Note: This method does more than apply str method as the string representations of the
    different encodings have a specific syntaxs.

    Args:
        phase_exp: Phase encosings to convert to string representations
        input_encoding: Encoding of the input phase exponents
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Raises:
        QiskitError: Invalid phase exponent encoding

    Returns:
        exp_str: string representations of given phase exponents

    Examples:
        >>> phase_exp = [0,1]
        >>> exp2expstr(phase_exp, '-is')
        '(-i,0)(-1,1)'

        >>> phase_exp = [0,1,2,3]
        >>> exp2expstr(phase_exp, '-i')
        array(['(-i,0)', '(-i,1)', '(-i,2)', '(-i,3)'], dtype='<U6')

    See Also:
        _exp2expstr, ...
    """
    if input_encoding not in PHASE_ENCODINGS:
        raise QiskitError(f"Invalid phase exponent encoding: {input_encoding}")

    scalar_pair = input_encoding in PHASE_ENCODINGS_ISMIS
    scalar = is_scalar(phase_exp, scalar_pair=scalar_pair) and same_type

    if input_encoding in PHASE_ENCODINGS_ISMIS:
        phase_exp = np.atleast_2d(phase_exp)
    else:
        phase_exp = np.atleast_1d(phase_exp)

    if scalar:
        return squeeze(_exp2expstr(phase_exp, input_encoding), scalar=scalar)
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
        >>> phase_exp = numpy.array([[0,1]])
        >>> _exp2expstr(phase_exp, '-is')
        array(['(-i,0)(-1,1)'], dtype='<U12')

        >>> phase_exp = numpy.array([0,1,2,3])
        >>> _exp2expstr(phase_exp, '-i')
        array(['(-i,0)', '(-i,1)', '(-i,2)', '(-i,3)'], dtype='<U6')

    See Also;
        exp2expstr
    """
    if encoding == "i":
        return np.array(["(i," + str(item) + ")" for item in phase_exp])
    if encoding == "-i":
        return np.array(["(-i," + str(item) + ")" for item in phase_exp])
    if encoding == "is":
        return np.array(["(i," + str(item[0]) + ")(-1," + str(item[1]) + ")" for item in phase_exp])
    if encoding == "-is":
        return np.array(
            ["(-i," + str(item[0]) + ")(-1," + str(item[1]) + ")" for item in phase_exp]
        )
    else:
        raise QiskitError(f"The encoding {encoding} is not supported.")


# ----------------------------------------------------------------------


def exp2exp(
    phase_exp: Union[np.ndarray, Any],
    input_encoding=INTERNAL_PHASE_ENCODING,
    output_encoding=DEFAULT_EXTERNAL_PHASE_ENCODING,
    same_type: bool = True,
):
    """Convert between the different phase exponents of encoded phase

    The possible phase encodings are:
        a) ['i'] :math:`i^r` where :math:`r \in \mathbb{Z}_4`
        b) ['-i'] :math:`(-i)^r` where :math:`r \in \mathbb{Z}_4`
        c) ['is'] :math:`i^r (-1)^s` where :math:`r,s \in \mathbb{Z}_2`
        d) ['-is'] :math:`(-i)^r (-1)^s` where :math:`r,s, \in \mathbb{Z}_2`

    Args:
        phase_exp: phase exponent(s)
        input_encoding (optional): Encoding for inputs phase exponents
            Defaults to INTERNAL_PHASE_ENCODING
        output_encoding (optional): Encoding for output phase exponents
            Defaults to DEFAULT_EXTERNAL_ENCODING
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Raises:
        QiskitError: Phase encoding is not supported or invalid

    Returns:
        output: phase exponent(s) encoded into new encoding

    Examples:
        >>> phase_exp = np.array([1,2,1,3,2])
        >>> exp2exp(phase_exp,'i','is')
        array([[1, 0],
               [0, 1],
               [1, 0],
               [1, 1],
               [0, 1]])

        >>> phase_exp = np.array([(0,1),(1,1),(1,0),(1,1)])
        >>> exp2exp(phase_exp,'is','-i')
        array([2, 1, 3, 1])

    """

    if input_encoding == output_encoding:
        return phase_exp
    if (input_encoding not in PHASE_ENCODINGS) or (output_encoding not in PHASE_ENCODINGS):
        raise QiskitError("Phase encoding is not supported or invalid")

    scalar_pair = input_encoding in PHASE_ENCODINGS_ISMIS
    scalar = is_scalar(phase_exp, scalar_pair=scalar_pair) and same_type
    if scalar_pair:
        phase_exp = np.atleast_2d(phase_exp)
    else:
        phase_exp = np.atleast_1d(phase_exp)

    if scalar:
        return squeeze(_exp2exp(phase_exp, input_encoding, output_encoding), scalar=scalar)
    return _exp2exp(phase_exp, input_encoding, output_encoding)


def _exp2exp(phase_exp, input_encoding, output_encoding):
    """Convert between the different phase exponents of encoded phase

    The possible phase encodings are:
        a) ['i'] :math:`i^r` where :math:`r \in \mathbb{Z}_4`
        b) ['-i'] :math:`(-i)^r` where :math:`r \in \mathbb{Z}_4`
        c) ['is'] :math:`i^r (-1)^s` where :math:`r,s \in \mathbb{Z}_2`
        d) ['-is'] :math:`(-i)^r (-1)^s` where :math:`r,s, \in \mathbb{Z}_2`

    Args:
        phase_exp: array of phase exponents
        input_encoding (optional): Encoding for inputs phase exponents
            Defaults to INTERNAL_PHASE_ENCODING
        output_encoding (optional): Encoding for output phase exponents
            Defaults to DEFAULT_EXTERNAL_ENCODING

    Raises:
        QiskitError: Phase encoding is not supported or invalid

    Returns:
        output: array of phase exponent(s) encoded into new encoding

    Examples:
        >>> phase_exp = np.array([1,2,1,3,2])
        >>> _exp2exp(phase_exp,'i','is')
        array([[1, 0],
               [0, 1],
               [1, 0],
               [1, 1],
               [0, 1]])

        >>> phase_exp = np.array([(0,1),(1,1),(1,0),(1,1)])
        >>> _exp2exp(phase_exp,'is','-i')
        array([2, 1, 3, 1])

    """
    #  Method: Time memory tradeoff. Basically a compressed lookup table
    #  The conversions are done via precomputed matrices that describe how the various
    #  indices change when converting between different formats. Tuple formats are
    #  linearized to use single integer indices for lookup and are then expanded
    #  at the end if neeed.

    #  Converstion example:
    #  :math:`i^r` to :math:`(-i)^s`  is encoded as 0 to 1 via _ENC matrix
    #  Find which transformation to use via _TRANS[0][1] = 1
    #  _CN[1]= [0,3,2,1] and so :math:`i^{[0,1,2,3]} = (-i)^{[0,3,2,1]}` etc
    #  This method works on vectors of encoded phase exponents

    # Binary expansion of an index
    _BI = [[0, 0], [0, 1], [1, 0], [1, 1]]
    # Format encoding
    _ENC = {"i": 0, "-i": 1, "is": 2, "-is": 3}
    # Conversion matrix split and compressed into two index matrices
    # Transformation indices
    _TI = [[0, 1, 2, 3], [0, 3, 2, 1], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2], [0, 1, 3, 2]]
    # Index to transformation matrices via (input_encoding, output_encoding) pairs
    _TRANS = [[0, 1, 2, 4], [1, 0, 4, 2], [2, 3, 0, 5], [3, 2, 5, 0]]
    # Conversion is done via precalculated tables that are stored in _CN, _DD and _TRANS
    # with the frmt encoded into _ENC
    input_encoding = _ENC[input_encoding]
    output_encoding = _ENC[output_encoding]
    # Linearize: Convert pairs to an index if needed

    if input_encoding > 1:
        # input frmt is in ['is', '-is']
        linear_phase_exp = 2 * phase_exp.T[0] + phase_exp.T[1]
    else:
        # input frmt is in ['i', '-i']
        linear_phase_exp = phase_exp
    # Calculate and then apply the transformation matrix
    trans = _TRANS[input_encoding][output_encoding]
    partial_encoded_phase_exp = np.asarray([_TI[trans][t] for t in linear_phase_exp])
    # Delinearize: Convert indexes back to pairs if needed
    if output_encoding > 1:
        encoded_phase_exp = np.asarray([_BI[t] for t in partial_encoded_phase_exp])
    else:
        encoded_phase_exp = partial_encoded_phase_exp
    return encoded_phase_exp


# ----------------------------------------------------------------------


def cpxstr2expstr(
    cpx_str: Union[np.ndarray, str], encoding: str, same_type: bool = True
) -> Union[np.ndarray, str]:
    """Converts a complex string represnetation into a exponent string representation.

    Args:
        cpx_str: String representating a complex pauli group phase
        encoding: Encoding to be used to encode the phase
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Returns:
        exp_str: String representing the encoded phase

    Examples:
        >>> cpxstr2expstr('1','i')
        '(i,0)'

        >>> cpxstr2expstr('-i','is')
        '(i,1)(-1,1)'

        >>> cpx_str = ['1','-1','i','-i']
        >>> cpxstr2expstr(cpx_str,'is')
        array(['(i,0)(-1,0)', '(i,0)(-1,1)',
              '(i,1)(-1,0)', '(i,1)(-1,1)'], dtype='<U11')

    See Also:
        TODO
    """
    cpx = cpxstr2cpx(cpx_str, same_type=same_type)
    phase_exp = cpx2exp(cpx, encoding, same_type=same_type)
    return exp2expstr(phase_exp, encoding, same_type=same_type)


def expstr2cpxstr(
    exp_str: Union[np.ndarray, str], encoding: str, same_type: bool = True, ones: bool = False
) -> Union[str, np.ndarray]:
    """Converts a string(s) representing a phase exponent into a string representing
    the equivalent complex number.

    Args:
        exp_str: String(s) representing the encoded Pauli group phases
        encoding (str): How the phase exponents are encoded
        same_type (optional): Scalar/Vector return flag. Defaults to True.
        ones (optional): ones (optional): Defines whether to present ones as "1" (True) or
            the empyty string (False). Defaults to False.

    Returns:
        cpx_str: String(s) representing the input phase exponent(s)

    Examples:
        >>> expstr2cpxstr('(-i,1)(-1,1)','-is')
        'i'

        >>> expstr2cpxstr('(-i,2)','-i',ones=True)
        '-1'

    See Also:
        TODO
    """
    phase_exp = expstr2exp(exp_str, same_type=same_type)
    cpx = exp2cpx(phase_exp, encoding, same_type=same_type)
    return cpx2cpxstr(cpx, same_type=same_type, ones=ones)


def cpxstr2exp(cpx_str: Union[np.ndarray, str], encoding: str, same_type: bool = True) -> Any:
    """Converts strings representing a complex pauli grioup phase to an
    encoded phase.

    Args:
        cpx_str: String(s) representing complex number Pauli groups phases as strings
        encoding (str): How the phase exponents are encoded
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Returns:
        phase_exp: Returns a encoded phase

    Examples:
        >>> cpxstr2exp('-i','i')
        3

        >>> cpxstr2exp('-i','-is')
        array([1, 0])

    See Also:
        TODO
    """
    cpx = cpxstr2cpx(cpx_str, same_type=same_type)
    return cpx2exp(cpx, encoding, same_type=same_type)


def exp2cpxstr(
    phase_exp: Any, encoding: str, same_type: bool = True, ones: bool = False
) -> Union[str, np.ndarray]:
    """Converts encoded Pauli groups phases into strings that represent
    the associated complex numbers

    Args:
        phase_exp: Phase exponents to convfert to strings representing complex numbers
        encoding: Encoding used to encoded phase exponents
        same_type (optional): Scalar/Vector return flag. Defaults to True.
        ones (optional): ones (optional): Defines whether to present ones as "1" (True) or
            the empyty string (False). Defaults to False.

    Returns:
        cpx_str: String representations of complex pauli group phases

    Examples:
        >>> exp2cpxstr(1,'i')
        'i'

        >>> rep.exp2cpxstr([1,1],'is')
        '-i'

        >>> phase_exp=np.array([[0,1],[1,1]])
        >>> exp2cpxstr(phase_exp,'is', ones=True)
        array(['-1', '-i'], dtype='<U2')

    See Also:
        TODO
    """
    cpx = exp2cpx(phase_exp, encoding, same_type=same_type)
    return cpx2cpxstr(cpx, same_type=same_type, ones=ones)


def expstr2cpx(
    phase_str: Union[np.ndarray, str], encoding: str, same_type: bool = True
) -> Union[np.ndarray, numbers.Complex]:
    """Converts strings representing phase exponents to complex numbers

    Args:
        phase_str: Strings representing phase exponents
        encoding: Encoding used to encoded phase exponents
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Returns:
        cpx : Complex number(s) representing the phase exponent string(s)

    Examples:
        >>> expstr2cpx('(i,1)','i')
        1j

        >>> expstr2cpx('(i,1)(-1,1)','is')
        (-0-1j)

        >>> exp_str = ['(i,0)','(i,1)','(i,1)','(i,3)']
        >>> expstr2cpx(exp_str,'i')
        array([ 1.+0.j,  0.+1.j,  0.+1.j, -0.-1.j])

    See Also:
        TODO
    """
    phase_exp = expstr2exp(phase_str, same_type=same_type)
    return exp2cpx(phase_exp, encoding, same_type=same_type)


def cpx2expstr(
    cpx: Union[np.ndarray, numbers.Complex], encoding: str, same_type: bool = True
) -> Union[np.ndarray, str]:
    """Converts complex Pauli group phases into strings representing phase exponents

    Args:
        cpx: Complex Pauli group phases
        encoding: Encoding to use to encode phases
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Returns:
        exp_str: Strings representing the encoded phases

    Examples:
        >>> cpx = [1j,1,-1j,-1]
        >>> cpx2expstr(cpx, '-is')
        array(['(-i,1)(-1,1)', '(-i,0)(-1,0)',
              '(-i,1)(-1,0)', '(-i,0)(-1,1)'], dtype='<U12')

        >>> cpx2expstr(1j, '-i')
        '(-i,3)'
    """
    phase_exp = cpx2exp(cpx, encoding, same_type=same_type)
    return exp2expstr(phase_exp, encoding, same_type=same_type)


# ----------------------------------------------------------------------


def str2exp(
    phase_str: Union[np.ndarray, str], encoding: str, same_type: bool = True
) -> Union[np.ndarray, str]:
    """Converts string representations of a Pauli group phases to and encoded phase exponents

    Args:
        phase_str: string(s) representing Pauli group phases
        encoding: Encoding to use to encode the phases
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Returns:
        phase_exp: Return phase exponent(s) encoded with the provided encoding

    Examples:
        >>> phase_str = "i"
        >>> str2exp(phase_str, 'i')
        1

        >>> phase_str = "(-i,3)"
        >>> str2exp(phase_str, 'i')
        1

        >>> phase_str = ""
        >>> str2exp(phase_str, 'i')
        0

        >>> phase_str = ["i", "-i", "+1", "1"]
        >>> str2exp(phase_str, '-is')
        array([[1, 1],
               [1, 0],
               [0, 0],
               [0, 0]])

        See Also:
            TODO
    """
    scalar = is_scalar(phase_str) and same_type
    phase_str = np.atleast_1d(phase_str)

    if scalar:
        return squeeze(_str2exp(phase_str, encoding), scalar=scalar)
    else:
        return _str2exp(phase_str, encoding)


# ----------------------------------------------------------------------


def _str2exp(phase_str: np.ndarray, encoding: str) -> np.ndarray:
    """Converts string representations of a Pauli group phases to and encoded phase exponents

    Args:
        phase_str: string(s) representing Pauli group phases
        encoding: Encoding to use to encode the phases

    Returns:
        phase_exp: Return phase exponent(s) encoded with the provided encoding

    Examples:
        >>> phase_str = numpy.array(["i"])
        >>> _str2exp(phase_str, 'i')
        array([1])

        >>> phase_str = numpy.array(["(-i,3)"])
        >>> _str2exp(phase_str, 'i')
        array([1])

        >>> phase_str = numpy.array([""])
        >>> _str2exp(phase_str, 'i')
        array([0])

        >>> phase_str = numpy.array(["i", "-i", "+1", "1"])
        >>> _str2exp(phase_str, '-is')
        array([[1, 1],
               [1, 0],
               [0, 0],
               [0, 0]])

        See Also:
            TODO
    """
    phase_enc = encode_of_phase_str(phase_str)

    def _trans(p_str, input_encoding, n, output_encoding):
        if input_encoding in PHASE_ENCODINGS:
            p_exp = expstr2exp(p_str)
            p_exp = exp2exp(p_exp, input_encoding, output_encoding)
        elif input_encoding == "complex":
            p_exp = cpxstr2exp(p_str, output_encoding)

        else:
            p_exp = np.zeros(shape=(1, n), dtype=np.int8)
        return p_exp

    n = phase_str.shape[0]
    return np.array(
        [_trans(item, input_enc, n, encoding) for item, input_enc in zip(phase_str, phase_enc)]
    )


# ---------------------------------------------------------------------
# Label parsing helper functions
# ---------------------------------------------------------------------


def split_pauli(pauli_str: str, same_type: bool = True) -> Tuple[str, str]:
    """Splits up a string reprenting a Pauli operator into the phase part and the Pauli part.

    Given a string representing a Pauli operator (with or without a phase) two strings
    will be returned. One for what is believed to be the phase part and one that is
    believed to the the phase free Pauli part. For will formed string representations
    the resulting string will be as expected.

    Phases can be in any of the following forms (as described in TODO: ref)

    'complex': i, -i, j, -j, 1i, 1j, -1j, -1i, ...
          'i': (i,0), (i,1), (i,2), (i,3), (j,0), ...
         '-i': (-i,0), (-i,1), (-i,2), (-i,3), (-j,0), ...
         'is': (i,0)(-1,1), (j,1)(-1,0), ...
        '-is': (-i,1)(-1,1), (-1j,0)(-1,1), ...

    The non-phase part of the Pauli string can be in one of the following forms

        Product Form Representations:

            'XZ' : '(X)(XZ)(Z)(I)(XZ)'
            'XZY': 'XYZIIIIIZYYX'
            'ZX' : '(Z)(I)(I)(ZX)(Z)'
            'YZX': 'XYYYXZIIIIXZ'

        Index Form Representations:

            'XZ' : '(X1)(Z3)(X4Z4)(Z7)
            'XZY': 'Z3Z6Y7Y7X4'
            'ZX' : '(Z3)(X4)(Z5X5)(Z10)
            'YZX': 'X5Z20Y21'

    Args:
        pauli_str: string describing the Pauli operator
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Returns:
        phase_str, pauli_str: Pauli string split into phase and Pauli part

    Examples:
        >>> string = "(-i,3)XIIXYZZX"
        >>> phase_str, pauli_str = split_pauli(string)
        >>> phase_str
        '(-i,3)'
        >>> pauli_str
        'XIIXYZZX'

        >>> string = "iX1Z3Y7"
        >>> phase_str, pauli_str = split_pauli(string)
        >>> phase_str
        'i'
        >>> pauli_str
        'X1Z3Y7'

        >>> string = "(-i,1)(-1,0)XXIXZ"
        >>> phase_str, pauli_str = split_pauli(string)
        >>> phase_str
        '(-i,1)(-1,0)'
        >>> pauli_str
        'XXIXZ'
    """

    scalar = is_scalar(pauli_str) and same_type

    pauli_str = np.atleast_1d(pauli_str)

    if scalar:
        phase_part, tensor_part = _split_pauli(pauli_str)
        return squeeze(phase_part, scalar=scalar), squeeze(tensor_part, scalar=scalar)
    return _split_pauli(pauli_str)


def _split_pauli(pauli_str: str):
    """Splits up a strings reprentings a Pauli operators into the phase and tensor parts

    See split_pauli for full details.

    Args:
        pauli_str: array of strings representing Pauli group elements

    Returns:
        out: array of split Pauli strings

    Examples:
        >>> pauli_str = ["(-i,1)X1Y2Z5Y4", "iX1Z5"]
        >>> phase_str, tensor_str = split_pauli(pauli_str, same_type=True)
        >>> phase_str
        array(['(-i,1)', 'i'], dtype='<U8')
        >>> tensor_str
        array(['X1Y2Z5Y4', 'X1Z5'], dtype='<U8')

    See Also:
        split_pauli
    """

    def _split(p_str):
        result = re.findall(SPLIT_PATTERN, p_str)
        if not result:
            raise QiskitError(f"Pauli tensor {pauli_str} string is not valid")
        return result[0][0], result[0][1]

    result = np.array([_split(item) for item in pauli_str])
    return result[:, 0], result[:, 1]


# ----------------------------------------------------------------------


def encode_of_phase_str(phase_str: str, same_type: bool = True) -> Union[np.ndarray, str]:
    """Returns how a string(s) representing a phase(s) are encoded

    Phases can be in any of the following forms (as described in TODO: ref)

    'complex': i, -i, j, -j, 1i, 1j, -1j, -1i, ...
          'i': (i,0), (i,1), (i,2), (i,3), (j,0), ...
         '-i': (-i,0), (-i,1), (-i,2), (-i,3), (-j,0), ...
         'is': (i,0)(-1,1), (j,1)(-1,0), ...
        '-is': (-i,1)(-1,1), (-1j,0)(-1,1), ...

    Args:
        phase_str: string representing a coefficient of a Pauli operator
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Raises:
        QiskitError: Unknown phase encoding

    Returns:
        endoding: How the input string(s) are encoded

    Examples:
        >>> string = "(i,0)(-1,1)"
        >>> enc = encode_of_phase_str(string)
        >>> enc
        'is'

        >>> string = "+1j"
        >>> enc = encode_of_phase_str(string)
        >>> enc
        'complex'

        >>> string = "(i,0)"
        >>> enc = encode_of_phase_str(string)
        >>> enc
        'i'

        >>> phase_str = array(['(-i,1)', 'i'])
        >>> rep.encode_of_phase_str(phase_str)
        array(['-i', 'complex'], dtype='<U7')
    """

    scalar = is_scalar(phase_str) and same_type

    phase_str = np.atleast_1d(phase_str)

    if scalar:
        return squeeze(_encode_of_phase_str(phase_str), scalar=scalar)
    return _encode_of_phase_str(phase_str).tolist()


def _encode_of_phase_str(phase_str: str) -> np.ndarray:
    """Returns how a string representing a phase is encoded

    Args:
        phase_str (str): _description_

    Returns:
        np.ndarray: _description_
    """

    def _find_encode(p_str):
        # Returns how a single string representing a phase is encoded
        if _is_pattern(p_str, pattern=PHASE_PATTERN["complex"]):
            return "complex"
        for encoding in PHASE_ENCODINGS:
            if _is_pattern(p_str, pattern=PHASE_PATTERN[encoding]):
                return encoding
        raise QiskitError(f"Unknown phase encoding: {encoding}")

    return np.array([_find_encode(item) for item in phase_str])


# ----------------------------------------------------------------------


def encode_of_tensor_str(
    tensor_str: str, encoded: bool = True, same_type: bool = True
) -> List[Tuple[List, Union[str, int]]]:
    """Returns how a string representing a Pauli tensor without a phase component
    is encoded.

    The non-phase part of the Pauli string can be in one of the following forms

        Product (syntax) Form Representations and examples:

            'XZ' : '(X)(XZ)(Z)(I)(XZ)'
            'XZY': 'XYZIIIIIZYYX'
            'ZX' : '(Z)(I)(I)(ZX)(Z)'
            'YZX': 'XYYYXZIIIIXZ'

        Index (syntax) Form Representations and exanples:

            'XZ' : '(X1)(Z3)(X4Z4)(Z7)
            'XZY': 'Z3Z6Y7Y7X4'
            'ZX' : '(Z3)(X4)(Z5X5)(Z10)
            'YZX': 'X5Z20Y21'

    If multiple representions fit the input string then all possible representations
    are returned.

    Args:
        tensor_str: string representation a pauli operator with no phase
        encoded (optional): The syntax type of the input string is either
            returned as a name (str) or an encoded integer (pauli_rep.PRODUCT
            or pauli_rep.INDEX). If encoded is set to True then the encoded
            integer is returned. If encoded is set to False then the name (str)
            is returned. Defaults to True.
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Raises:
        QiskitError: Unknown Pauli syntax type

    Returns:
        encoding_rep, syntax_type

    Examples:
        >>> string = "X1Y2Z5Y4"
        >>> encoding_rep, syn = encode_of_tensor_str(string, encoded=False)
        >>> encoding_rep
        ['XZY', 'YZX']
        >>> syn
        'index'

        >>> string = "XYXZZIX"
        >>> encoding_rep, syn = encode_of_tensor_str(string, encoded=False)
        >>> encoding_rep
        ['XZY', 'YZX']
        >>> syn = pauli_rep.PRODUCT
        True

        >>> string = (X1)(Z2)(Z3X3)(Z10)
        >>> encoding_rep, syn = encode_of_tensor_str(string, encoded=False)
        >>> encoding_rep
        ['ZX']
        >>> syn
        'index'
    """
    scalar = is_scalar(tensor_str) and same_type

    tensor_str = np.atleast_1d(tensor_str)

    if scalar:
        return squeeze(_encode_of_tensor_str(tensor_str, encoded), scalar=scalar)
    return _encode_of_tensor_str(tensor_str, encoded)


def _encode_of_tensor_str(
    tensor_str: np.ndarray, encoded: bool
) -> List[Tuple[List, Union[str, int]]]:
    """_summary_

    Args:
        tensor_str (np.ndarray): _description_
        encoded (bool): _description_

    Returns:
        Tuple[List, Union[str, int]]: _description_
    """

    def _find_encode(t_str: str, ecode: bool):
        encoding_rep = []
        syntax_type = None
        if _is_pattern(t_str, pattern=LATEX_INDICATE_PATTERN):
            syntax_type = LATEX_SYNTAX
            for encoding in TENSOR_ENCODINGS:
                if _is_pattern(t_str, pattern=TENSOR_LATEX_PATTERN[encoding]):
                    encoding_rep.append(encoding)

        elif _is_pattern(t_str, pattern=INDEX_INDICATE_PATTERN):
            # String using index syntax
            syntax_type = INDEX_SYNTAX
            for encoding in TENSOR_ENCODINGS:
                if _is_pattern(t_str, pattern=TENSOR_INDEX_PATTERN[encoding]):
                    encoding_rep.append(encoding)

        else:
            # string using product syntax
            syntax_type = PRODUCT_SYNTAX
            for encoding in TENSOR_ENCODINGS:
                if _is_pattern(t_str, pattern=TENSOR_PRODUCT_PATTERN[encoding]):
                    encoding_rep.append(encoding)

        try:
            _ = encoding_rep[0]
            if ecode:
                return encoding_rep, syntax_type
            else:
                return encoding_rep, SYNTAX_TO_TEXT[syntax_type]
        except IndexError as unknown_syntax:
            raise QiskitError(f"Unknown Pauli syntax type: {t_str}") from unknown_syntax

    return [_find_encode(item, encoded) for item in tensor_str]


# ----------------------------------------------------------------------


def str2symplectic(
    pauli_str: Union[np.ndarray, str],
    qubit_order: str = "right-to-left",
    output_encoding: Optional[str] = INTERNAL_PAULI_ENCODING,
    index_start: int = 0,
    same_type: bool = True,
) -> Tuple[np.ndarray, Union[np.array, Any]]:
    """Converts strings of Pauli group elements into phase exponents and their
    symplectic matrix representation.

    Args:
        pauli_str: Strings representing Pauli group elements
        qubit_order (optional): order in which to read product representation Paulis.
            Defaults to "right-to-left". Alternative is "left-to-right".
        output_encoding (optional): Pauli encoding for phase_exponent matrix pair.
            Defaults to INTERNAL_PAULI_ENCODING.
        index_start (optional): Lowest value for index in index syntax tensors.
            Defaults to 0
        same_type (optional): Scalar/Vector return flag. Defaults to True.

    Raises:
        QiskitError: Negative index_start values are not supported
        QiakitError: Pauli encoding is not supported or invalid

    Returns:
        symplectic_matrix, phase_exp : Returns the phase exponent(s) and s
        ymplectic matrix of the inputs

    Examples:
        >>> pauli_str = "iXXIZY"
        >>> matrix, phase_exp = str2symplectic(pauli_str,
                                               qubit_order="left-to-right",
                                               output_encoding='-isXZ')
        >>> phase_exp
        array([0, 1], dtype=int8)
        >>> matrix.astype(int)
        array([[1, 1, 0, 0, 1, 0, 0, 0, 1, 1]])

        >>> pauli_str = ["iYII", "-iX0Z2", "X1Z2"]
        >>> matrix, phase_exp = str2symplectic(pauli_str,
                                               qubit_order="left-to-right",
                                               output_encoding='-iXZY')
        >>> phase_exp
        array([3, 1, 0], dtype=int8)
        >>> matrix.astype(int)
        array([[1, 0, 0, 1, 0, 0],
               [1, 0, 0, 0, 0, 1],
               [0, 1, 0, 0, 0, 1]])

        >>> pauli_str = "iX1X3Y4Z9"
        >>> matrix, phase_exp = str2symplectic(pauli_str,
                                               qubit_order="left-to-right",
                                               output_encoding='-isXZ',
                                               index_start=1)
        >>> phase_exp
        array([0, 1], dtype=int8)
        >>> matrix.astype(int)
        array([[1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]])

        >>> matrix, phase_exp = rep.str2symplectic("iXXXIZYZ")
        >>> phase
        0
        >>> matrix.astype(int)
        array([[0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]])

    See Also:
        TODO
    """

    if index_start < 0:
        raise QiskitError("Negative index_start values are not supported")

    if output_encoding not in PAULI_ENCODINGS:
        raise QiskitError("Pauli encoding is not supported or invalid")

    scalar = is_scalar(pauli_str) and same_type

    pauli_str = np.atleast_1d(pauli_str)

    if scalar:
        matrix, phase_exp = _str2symplectic(pauli_str, qubit_order, output_encoding, index_start)
        return matrix, squeeze(phase_exp, scalar=scalar)
    return _str2symplectic(pauli_str, qubit_order, output_encoding, index_start)


def _str2symplectic(
    pauli_str: np.ndarray,
    qubit_order: str,
    output_encoding: str,
    index_start: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Converts strings of Pauli group elements into phase exponents and their
    symplectic matrix representation.

    Args:
        pauli_str: Strings representing Pauli group elements
        qubit_order: order in which to read product representation Paulis.
            Defaults to "right-to-left". Alternative is "left-to-right".
        output_encoding: Pauli encoding for phase_exponent matrix pair.
        index_start: Lowest value for index in index syntax tensors. Defaults to 0

    Returns:
        symplectic_matrix, phase_exp : Returns the phase exponent(s) and
        symplectic matrix of the inputs
    """

    phase_str, tensor_str = split_pauli(pauli_str, same_type=False)
    phase_exp = str2exp(phase_str, INTERNAL_PHASE_ENCODING)
    tensor_enc = encode_of_tensor_str(tensor_str)
    index_store = []
    tensor_store = []
    num_qubits_store = []

    _TRANS = {"I": "I", "X": "X", "Z": "Z", "XZ": "Y", "ZX": "Y"}

    for tensor, (tensor_encodings, syntax) in zip(tensor_str, tensor_enc):
        if syntax == PRODUCT_SYNTAX:
            if "XZ" in tensor_encodings:
                # Make this better as below
                new_tensor = re.findall(ENC_PRODUCT_XZ_CAP, tensor)
                new_tensor = [_TRANS[item] for item in new_tensor]
            elif "ZX" in tensor_encodings:
                # Make this better as below
                new_tensor = re.findall(ENC_PRODUCT_ZX_CAP, tensor)
                new_tensor = [_TRANS[item] for item in new_tensor]
            else:
                new_tensor = tensor

            num_qubits = len(new_tensor)

            if qubit_order == "right-to-left":
                indices = list(reversed(range(num_qubits)))
            else:
                indices = list(range(num_qubits))

        elif syntax == INDEX_SYNTAX:
            if "XZ" in tensor_encodings:
                new_tensor = re.sub(r"X([0-9]+)Z[0-9]+", r"Y\1", tensor)
                new_tensor = re.sub(r"[\(\)]", "", new_tensor)
            elif "ZX" in tensor_encodings:
                new_tensor = re.sub(r"Z([0-9]+)X[0-9]+", r"Y\1", tensor)
                new_tensor = re.sub(r"[\(\)]", "", new_tensor)
            else:
                new_tensor = tensor

            indices = list(map(int, re.findall("\d+", tensor)))
            num_qubits = max(indices) - index_start + 1
            new_tensor = re.findall(f"{PAULI_REGEX}+", tensor)
        elif syntax == LATEX_SYNTAX:
            if "XZ" in tensor_encodings:
                new_tensor = re.sub(r"X_\{([0-9]+\{)Z_\{[0-9]+\}", r"Y\1", tensor)
                new_tensor = re.sub(r"[\(\)]", "", new_tensor)
            elif "ZX" in tensor_encodings:
                new_tensor = re.sub(r"Z_\{([0-9]+\{)X_\{[0-9]+\}", r"Y\1", tensor)
                new_tensor = re.sub(r"[\(\)]", "", new_tensor)
            else:
                new_tensor = tensor

            indices = list(map(int, re.findall("\d+", tensor)))
            num_qubits = max(indices) - index_start + 1
            new_tensor = re.findall(f"{PAULI_REGEX}+", tensor)
        else:
            raise QiskitError(f"Unknown input syntax: {syntax}")
        index_store.append(indices)
        tensor_store.append(new_tensor)
        num_qubits_store.append(num_qubits)

    num_qubits = max(num_qubits_store)
    matrix = np.zeros(shape=(len(tensor_store), 2 * num_qubits), dtype=np.bool_)
    for row_index, (indices, tensor, (tensor_encodings, syntax)) in enumerate(
        zip(index_store, tensor_store, tensor_enc)
    ):
        # Fill in the symplectic representation
        for i, char in enumerate(tensor):
            if syntax == PRODUCT_SYNTAX:
                index = indices[i]
            else:
                index = indices[i] - index_start
            if char == "X":
                matrix[row_index, index] = True
            elif char == "Z":
                matrix[row_index, index + num_qubits] = True
            elif char == "Y":
                matrix[row_index, index] = True
                matrix[row_index, index + num_qubits] = True

    out_phase_enc, _ = split_pauli_enc(output_encoding)
    if out_phase_enc in PHASE_ENCODINGS_ISMIS:
        new_phase_exp = np.zeros(shape=(phase_exp.shape[0], 2), dtype=np.int8)
    else:
        new_phase_exp = np.zeros(shape=phase_exp.shape, dtype=np.int8)

    y_count = count_num_y(matrix, scalar=False)

    for index, (p_exp, count, (tensor_encodings, syntax)) in enumerate(
        zip(phase_exp, y_count, tensor_enc)
    ):
        n_exp = change_pauli_encoding(
            p_exp,
            count,
            input_pauli_encoding=INTERNAL_PHASE_ENCODING + tensor_encodings[0],
            output_pauli_encoding=output_encoding,
        )
        new_phase_exp[index] = n_exp

    return matrix, new_phase_exp


# ----------------------------------------------------------------------
def symplectic2str(
    matrix: np.ndarray,
    phase_exp: Any = None,
    input_encoding: str = INTERNAL_PAULI_ENCODING,
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
        matrix: GF(2) symplectic matrix
        phase_exp (optional): Phase exponent(s) for matrix. A value of
            None will lead to unity phases. Defaults to None.
        input_encoding (optional): Pauli encoding of phase relative to
            matrix. Defaults to INTERNAL_PAULI_ENCODING.
        output_phase_encoding (optional): Encoding used to represent phases.
            A value of None will result in complex phases notation. Defaults
            to None.
        no_phase (optional): When set to True, no phase will appear no matter
            what encoding is selected. So the symplectic matrix [1, 1] will produce
            the operator Y in 'XZY' encoding but also (XZ) in the 'XZ' encoding which
            are different operators if phases are considered.
        output_tensor_encoding (optional): Encoding of pauli tensor
            (without phase). Defaults to DEFAULT_EXTERNAL_TENSOR_ENCODING.
        syntax (optional): Syntax of pauli tensor. Values are
            PRODUCT_SYNTAX = 0, INDEX_SYNTAX=1 and LATEX_SYNTAX=2. Defaults to INDEX_SYNTAX.
        qubit_order (optional): Order in which qubits are read. options aree
            "right-to-left" and "left-to-right". Defaults to "right-to-left".
        index_start (optional): Lowest value for index in index syntax tensors.
            Defaults to 0
        same_type (optional): Scalar/Vector return flag. Defaults to True.
        index_str (optional): String that get inserted between operator and numbers in
            index format. Default is "".

    Raises:
        QiskitError: Unsupport syntax

    Returns:
        pauli_str: Pauli strings

    Examples:
        >>> matrix = numpy.array([[0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, \
                                   0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
                                  [0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, \
                                   0, 0, 0, 0, 0, 1, 1, 0, 0, 1]], dtype=numpy.bool_)
        >>> phase_exp = numpy.array([3, 1], dtype=numpy.int8)

        >>> symplectic2str(matrix, phase_exp)
        array(['-iY9Y4X3X1', 'iZ10X8Y7Y6X3'], dtype='<U12')

        >>> symplectic2str(matrix, phase_exp, qubit_order="left-to-right")
        array(['-iX1X3Y4Y9', 'iX3Y6Y7X8Z10'], dtype='<U12')

        >>> symplectic2str(matrix)
        array(['-Y9Y4X3X1', '-Z10X8Y7Y6X3'], dtype='<U12')

        >>> symplectic2str(matrix[0])
       '-Y9Y4X3X1'

        >>> symplectic2str(matrix[0], no_phase=True)
        'Y9Y4X3X1'

        >>> symplectic2str(matrix,
                           phase_exp,
                           input_encoding='iXZ',
                           output_phase_encoding=None,
                           output_tensor_encoding='XZY',
                           syntax=rep.INDEX_SYNTAX)
        array(['iY9Y4X3X1', '-iZ10X8Y7Y6X3'], dtype='<U13')

        >>> symplectic2str(matrix,
                           phase_exp,
                           input_encoding='iXZ',
                           output_phase_encoding=None,
                           output_tensor_encoding='XZ',
                           syntax=rep.INDEX_SYNTAX)
        array(['-i(X9Z9)(X4Z4)(X3)(X1)', 'i(Z10)(X8)(X7Z7)(X6Z6)(X3)'], dtype='<U26')

        >>> symplectic2str(matrix,
                           phase_exp,
                           input_encoding='iXZ',
                           output_phase_encoding='-is',
                           output_tensor_encoding='XZ',
                           syntax=rep.PRODUCT_SYNTAX)
        array(['(-i,1)(-1,0)(I)(XZ)(I)(I)(I)(I)(XZ)(X)(I)(X)(I)',
               '(-i,1)(-1,1)(Z)(I)(X)(XZ)(XZ)(I)(I)(X)(I)(I)(I)'], dtype='<U47')

        >>> symplectic2str(matrix,
                           phase_exp,
                           input_encoding='iXZ',
                           output_phase_encoding='-is',
                           output_tensor_encoding='XZY',
                           syntax=rep.INDEX_SYNTAX,
                           index_start=2)
        >>> array(['(-i,1)(-1,1)Y11Y6X5X3', '(-i,1)(-1,0)Z12X10Y9Y8X5'], dtype='<U24')

        >>> symplectic2str(matrix,
                           phase_exp,
                           input_encoding='iXZ',
                           output_phase_encoding='-is',
                           qubit_order="left-to-right",
                           output_tensor_encoding='XZY',
                           syntax=rep.INDEX_SYNTAX,
                           index_start=1,
                           index_str="_")
        array(['(-i,1)(-1,1)X_2X_4Y_5Y_10',
              '(-i,1)(-1,0)X_4Y_7Y_8X_9Z_11'], dtype='<U28')
    """
    matrix = np.atleast_2d(matrix)
    if no_phase:
        phase_str = np.full((matrix.shape[0],), "")
    else:
        if phase_exp is None:
            phase_exp = np.zeros(shape=(matrix.shape[0],), dtype=np.int8)
        else:
            phase_exp = np.atleast_1d(phase_exp)

        y_count = count_num_y(matrix, scalar=False)

        if output_phase_encoding is None:
            temp_enc = "-i" + output_tensor_encoding

            phase_exp = change_pauli_encoding(
                phase_exp,
                y_count,
                input_pauli_encoding=input_encoding,
                output_pauli_encoding=temp_enc,
            )

            phase_str = exp2cpxstr(phase_exp, get_phase_enc(temp_enc))
        else:
            phase_exp = change_pauli_encoding(
                phase_exp,
                y_count,
                input_pauli_encoding=input_encoding,
                output_pauli_encoding=output_phase_encoding + output_tensor_encoding,
            )
            phase_str = exp2expstr(phase_exp, output_phase_encoding)

    tensor_str = []

    _YENC = ["I", "X", "Z", "Y"]
    _XZENC = ["(I)", "(X)", "(Z)", "(XZ)"]
    _ZXENC = ["(I)", "(X)", "(Z)", "(ZX)"]
    _ENC = {"XZY": _YENC, "YZX": _YENC, "XZ": _XZENC, "ZX": _ZXENC}

    num_qubits = matrix.shape[1] >> 1
    if syntax == PRODUCT_SYNTAX:
        for pauli in matrix:
            tmp_tensor_str = ""
            marker = pauli[:num_qubits].astype(np.int8) + 2 * pauli[num_qubits:].astype(np.int8)
            for mark in marker:
                if qubit_order == "left-to-right":
                    tmp_tensor_str += _ENC[output_tensor_encoding][mark]
                else:
                    tmp_tensor_str = _ENC[output_tensor_encoding][mark] + tmp_tensor_str
            tensor_str.append(tmp_tensor_str)
    elif syntax in (INDEX_SYNTAX, LATEX_SYNTAX):
        if syntax == LATEX_SYNTAX:
            str_repr = _ind_to_latex_repr
        else:
            str_repr = str

        for pauli in matrix:
            tmp_tensor_str = ""
            marker = pauli[:num_qubits].astype(np.int8) + 2 * pauli[num_qubits:].astype(np.int8)
            if output_tensor_encoding in Y_TENSOR_ENCODINGS:
                for index, mark in enumerate(marker):
                    if mark != 0:
                        if qubit_order == "left-to-right":
                            tmp_tensor_str += (
                                _ENC[output_tensor_encoding][mark]
                                + index_str
                                + str_repr(index + index_start)
                            )
                        else:
                            tmp_tensor_str = (
                                _ENC[output_tensor_encoding][mark]
                                + index_str
                                + str_repr(index + index_start)
                                + tmp_tensor_str
                            )
            else:
                for index, mark in enumerate(marker):
                    if mark != 0:
                        if mark == 3:
                            if output_tensor_encoding == "XZ":
                                if qubit_order == "left-to-right":
                                    tmp_tensor_str += (
                                        "(X"
                                        + index_str
                                        + str_repr(index + index_start)
                                        + "Z"
                                        + index_str
                                        + str_repr(index + index_start)
                                        + ")"
                                    )
                                else:
                                    tmp_tensor_str = (
                                        "(X"
                                        + index_str
                                        + str_repr(index + index_start)
                                        + "Z"
                                        + index_str
                                        + str_repr(index + index_start)
                                        + ")"
                                        + tmp_tensor_str
                                    )
                            else:
                                if qubit_order == "left-to-right":
                                    tmp_tensor_str += (
                                        "(Z"
                                        + index_str
                                        + str_repr(index + index_start)
                                        + "X"
                                        + index_str
                                        + str_repr(index + index_start)
                                        + ")"
                                    )
                                else:
                                    tmp_tensor_str = (
                                        "(Z"
                                        + index_str
                                        + str_repr(index + index_start)
                                        + "X"
                                        + index_str
                                        + str_repr(index + index_start)
                                        + ")"
                                        + tmp_tensor_str
                                    )
                        else:
                            if qubit_order == "left-to-right":
                                tmp_tensor_str += (
                                    "("
                                    + _YENC[mark]
                                    + index_str
                                    + str_repr(index + index_start)
                                    + ")"
                                )
                            else:
                                tmp_tensor_str = (
                                    "("
                                    + _YENC[mark]
                                    + index_str
                                    + str_repr(index + index_start)
                                    + ")"
                                    + tmp_tensor_str
                                )

            tensor_str.append(tmp_tensor_str)
    else:
        raise QiskitError(f"Unsupport syntax: {syntax}")

    result = [p_str + t_str for p_str, t_str in zip(phase_str, tensor_str)]
    if matrix.shape[0] == 1 and same_type:
        return result[0]
    else:
        return np.array(result)


# ----------------------------------------------------------------------
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
    """Converts between different string representations of Pauli operators

    Args:
        pauli_str: Strings representing Pauli group elements
        syntax_output (optional): Syntax of output pauli tensor. Values are
            PRODUCT_SYNTAX = 0, INDEX_SYNTAX=1 and LATEX_SYNTAX=2. Defaults to INDEX_SYNTAX.
        phase_encoding_output_string (optional): Encoding used to represent phases of the output.
            A value of None will result in complex phases notation. Defaults
            to None.
        tensor_encoding_output_string (optional): Encoding of output pauli tensor
            (without phase). Defaults to DEFAULT_EXTERNAL_TENSOR_ENCODING.
        qubit_order_input (optional): order in which to read product representation Paulis of the input.
            Defaults to "right-to-left". Alternative is "left-to-right". Only relevatent if the syntax of
            the input is in PRODUCT_SYNTAX.
        qubit_order_output (optional): order in which to read product representation Paulis of the
            output. Defaults to "right-to-left". Alternative is "left-to-right".
        index_start_input (optional): Lowest value for index in index syntax tensors for the input.
            Defaults to 0
        index_start_output (optional): Lowest value for index in index syntax tensors for the output.
            Defaults to 0


    Returns:
        pauli_str: Pauli strings
    """

    matrix, new_phase_exp = str2symplectic(
        pauli_str, index_start=index_start_input, qubit_order=qubit_order_input
    )
    return symplectic2str(
        matrix,
        new_phase_exp,
        qubit_order=qubit_order_output,
        syntax=syntax_output,
        index_start=index_start_output,
        output_phase_encoding=phase_encoding_output_string,
        output_tensor_encoding=tensor_encoding_output_string,
    )


# ----------------------------------------------------------------------
# Array(s) to Symplectic
# ----------------------------------------------------------------------


def from_array(
    matrix: Union[List, Tuple, np.ndarray],
    phase_exp: Union[int, List, Tuple, np.ndarray] = None,
    input_pauli_encoding: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert array and phase_exp to the matrix and phase_exp in BasePauli's internal
    Pauli encoding (pauli_rep.INTERNAL_PAULI_ENCODING)
    Args:
        matrix (_type_): _description_
        phase_exp (_type_): _description_
        input_pauli_encoding: input Pauli encoding
    Returns:
        _type_: _description_
    """
    if input_pauli_encoding is None:
        input_pauli_encoding = DEFAULT_EXTERNAL_PAULI_ENCODING
    if isinstance(matrix, np.ndarray) and matrix.dtype == bool:
        matrix_data = matrix
    else:
        matrix_data = np.asarray(matrix, dtype=bool)
    matrix_data = np.atleast_2d(matrix_data)
    if not is_symplectic_matrix_form(matrix_data):
        raise QiskitError("Input matrix not a symplectic matrix or symplectic vector")
    if phase_exp is None or (isinstance(phase_exp, numbers.Integral) and phase_exp == 0):
        phase_exp = np.zeros(shape=(matrix_data.shape[0],))
    y_count = count_num_y(matrix_data)
    in_phase_exp = change_pauli_encoding(
        phase_exp,
        y_count,
        input_pauli_encoding=input_pauli_encoding,
        output_pauli_encoding=INTERNAL_PAULI_ENCODING,
        same_type=False,
    )
    return matrix_data, in_phase_exp


def from_split_array(
    x: Union[List, Tuple, np.ndarray],
    z: Union[List, Tuple, np.ndarray],
    phase_exp: Union[int, List, Tuple, np.ndarray],
    input_pauli_encoding: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert split array (separate X and Z arrays) and phase_exp to the matrix
    and phase_exp in BasePauli's internal
    Pauli encoding (pauli_rep.INTERNAL_PAULI_ENCODING)

    Args:
        x (Union[List, Tuple, np.ndarray]): _description_
        z (Union[List, Tuple, np.ndarray]): _description_
        phase_exp (Union[int, List, Tuple, np.ndarray]): _description_
        input_pauli_encoding (Optional[str], optional): _description_. Defaults to None.

    Raises:
        QiskitError: _description_

    Returns:
        Tuple[np.ndarray, np.ndarray]: _description_
    """

    if input_pauli_encoding is None:
        input_pauli_encoding = DEFAULT_EXTERNAL_PAULI_ENCODING
    if isinstance(x, np.ndarray) and x.dtype == bool:
        x_data = x
    else:
        x_data = np.asarray(x, dtype=bool)
    x_data = np.atleast_2d(x_data)
    if isinstance(z, np.ndarray) and z.dtype == bool:
        z_data = z
    else:
        z_data = np.asarray(z, dtype=bool)
    z_data = np.atleast_2d(z_data)
    matrix = np.hstack((x_data, z_data))
    if not is_symplectic_matrix_form(matrix):
        raise QiskitError("Input matrix not a symplectic matrix or symplectic vector")
    y_count = count_num_y(matrix)
    in_phase_exp = change_pauli_encoding(
        phase_exp,
        y_count,
        input_pauli_encoding=input_pauli_encoding,
        output_pauli_encoding=INTERNAL_PAULI_ENCODING,
        same_type=False,
    )
    return matrix, in_phase_exp


# ----------------------------------------------------------------------
# Symplectic to Complex Matrix conversion
# ----------------------------------------------------------------------


def to_cpx_matrix(
    matrix: np.ndarray,
    phase_exp: Optional[np.ndarray],
    pauli_encoding: str = INTERNAL_PAULI_ENCODING,
    sparse: bool = False,
) -> Union[np.ndarray, csr_matrix]:
    """Return the complex matrix representation from its symplectic representation with phase.
    Args:
        matrix (np.ndarray): _description_
        phase_exp (Optional[np.ndarray]): _description_
        pauli_encoding (str, optional): _description_. Defaults to pauli_rep.INTERNAL_PAULI_ENCODING.
        sparse (bool, optional): _description_. Defaults to False.
    Raises:
        QiskitError: Input matrix is not a symplectic GF(2) matrix
        QiskitError: Can only convert a single Pauli operator to a complex matrix representation
    Returns:
        matrix: complex matrix representation of Pauli operator.
    """
    if not is_symplectic_matrix_form(matrix):
        raise QiskitError("Input matrix is not a symplectic GF(2) matrix")

    matrix = np.atleast_2d(matrix)
    if matrix.shape[0] != 1:
        raise QiskitError(
            "Can only convert a single Pauli operator to a complex matrix representation"
        )

    num_qubits = matrix.shape[1] >> 1

    if phase_exp is None:
        phase_exp = 0
    else:
        phase_exp = np.atleast_1d(phase_exp)

    # Convert to internal Pauli encoding if needed
    if pauli_encoding != INTERNAL_PAULI_ENCODING:
        num_y = count_num_y(matrix)
        phase_exp = change_pauli_encoding(
            phase_exp,
            num_y,
            input_pauli_encoding=pauli_encoding,
            output_pauli_encoding=INTERNAL_PAULI_ENCODING,
        )
    phase_exp = phase_exp[0]
    matrix = np.squeeze(matrix)
    return _to_cpx_matrix(matrix, phase_exp, num_qubits, sparse=sparse)


def _to_cpx_matrix(
    matrix: np.ndarray, phase_exp: int, num_qubits: int, sparse: bool = False
) -> Union[np.ndarray, csr_matrix]:
    """Return the complex matrix representation from its symplectic representation with phase.
    Args:
        z (array): The symplectic representation z vector.
        x (array): The symplectic representation x vector.
        phase_exp (int): Pauli phase.
        sparse (bool): Optional. Of True return a sparse CSR matrix,
                       otherwise return a dense Numpy array
                       (default: False).
    Returns:
        array: if sparse=False.
        csr_matrix: if sparse=True.
    """
    dim = 2**num_qubits
    twos_array = 1 << np.arange(num_qubits)
    x_indices = np.asarray(matrix[:num_qubits]).dot(twos_array)
    z_indices = np.asarray(matrix[num_qubits:]).dot(twos_array)
    indptr = np.arange(dim + 1, dtype=np.uint)
    indices = indptr ^ x_indices
    if phase_exp:
        coeff = (-1j) ** phase_exp
    else:
        coeff = 1
    data = np.array([coeff * (-1) ** (bin(i).count("1") % 2) for i in z_indices & indptr])
    if sparse:
        # Return sparse matrix
        return csr_matrix((data, indices, indptr), shape=(dim, dim), dtype=complex)
    # Build dense matrix using csr format
    mat = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        mat[i][indices[indptr[i] : indptr[i + 1]]] = data[indptr[i] : indptr[i + 1]]
    return mat


# ----------------------------------------------------------------------
# Utility functions
# ----------------------------------------------------------------------


def indices_to_boolean(indices: Iterable[int], dim: int) -> np.ndarray:
    """Converts an array/list of indices to a numpy boolean vector

    Given an input iterator of indices [x_0,x_2,...,X_k] and a vector dimension
    n create the binary/bool vector b such that b[x_j]=True for all j with b[h]=False
    for h not x_j where the length of b is dim.

    Args:
        indices (list): Non negative integers (indices)
        dim (integer): Dimension/length of vector to be returned

    Returns:
        np.array: Boolean vector of dimension dim
    """
    row = np.zeros(dim, dtype=bool)
    for k in indices:
        row[k] = True
    return row


# ----------------------------------------------------------------------
def _ind_to_latex_repr(index: int) -> str:
    """Adds curly braces and an underscore to an index.

    Args:
        index: An integer

    Returns:
        str: string in LaTeX syntax
    """
    return f"_{{{index}}}"


def boolean_to_indices(booleans: Iterable[bool]) -> np.ndarray:
    """Converts a array of boolean values into a index numpy array of\
        non-zero-elements

    Args:
        booleans (boolean array): An array of boolean values

    Returns:
        numpy.ndarray: Array of indices of non-zero values
    """
    return np.nonzero(np.asarray(booleans))


def scalar_op2symplectic(
    op: ScalarOp, output_encoding: str = DEFAULT_EXTERNAL_PHASE_ENCODING
) -> Tuple[np.ndarray, Union[np.array, Any]]:
    """Convert a ScalarOp to symplectic representation with phase.

    TODO: Allow this to work on arrays of ScalarOps

    Args:
        op: Input scalarOp
        output_encoding: Phase encoding to use to encode phase from ScalarOp.
            Default is INTERNAL_PHASE_ENCODING='-i'

    Raises:
        QiskitError: Operator is not an N-qubit identity

    Returns:
        matrix, phase_exponent: GF(2) symplectic matrix and phase_exponent
        representing ScalarOp
    """
    if op.num_qubits is None:
        raise QiskitError(f"{op} is not an N-qubit identity")
    matrix = np.zeros(shape=(1, 2 * op.num_qubits), dtype=np.bool_)
    phase_exp = cpx2exp(op.coeff, output_encoding=output_encoding)
    return matrix, phase_exp


# ----------------------------------------------------------------------


def gate2symplectic(
    gate: Gate, encoding: str = INTERNAL_PAULI_ENCODING
) -> Tuple[np.ndarray, Union[np.array, Any]]:
    """Converts a Pauli gate to a symplectic matrix with phase

    Args:
        gate: Gate
        encoding (optional): Pauli encoding to encode symplectic matrix with phase.
            Defaults to DEFAULT_EXTERNAL_PAULI_ENCODING='-iYZX';

    Raises:
        QiskitError: Invalid Pauli instruction

    Returns:
        matrix, phase_exp: phase exponent and symplectic matrix
    """
    if isinstance(gate, PauliGate):
        return str2symplectic(gate.params[0], output_encoding=encoding)
    if isinstance(gate, IGate):
        return str2symplectic("I", output_encoding=encoding)
    if isinstance(gate, XGate):
        return str2symplectic("X", output_encoding=encoding)
    if isinstance(gate, YGate):
        return str2symplectic("Y", output_encoding=encoding)
    if isinstance(gate, ZGate):
        return str2symplectic("Z", output_encoding=encoding)
    raise QiskitError("Invalid Pauli instruction.")


# ----------------------------------------------------------------------
