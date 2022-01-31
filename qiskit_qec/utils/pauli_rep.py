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

# TODO finish documenting

# pylint: disable=invalid-name,anomalous-backslash-in-string
# pylint: disable=bad-docstring-quotes  # for deprecate_function decorator

import numbers
import re

import numpy as np
from qiskit.exceptions import QiskitError

# -------------------------------------------------------------------------------
# Module Variables/States
# -------------------------------------------------------------------------------

# Set the interal representations
# The internal representation cannot be changed by changing this constant
# This constants are for reference only and do not change the behavio(u_vec)r of
# the PauliListBase methods. See [ref] for details on the different representation
# formats
# TODO: Include ref for above.

INTERNAL_SYMP_REP_FORMAT = "ZX"
INTERNAL_PHASE_REP_FORMAT = "-i"
INTERNAL_PAULI_REP_FORMAT = INTERNAL_PHASE_REP_FORMAT + INTERNAL_SYMP_REP_FORMAT

DEFAULT_EXTERNAL_SYMP_FORMAT = "YZX"
DEFAULT_EXTERNAL_PHASE_REP_FORMAT = "-i"
DEFAULT_EXTERNAL_PAULI_REP_FORMAT = DEFAULT_EXTERNAL_PHASE_REP_FORMAT + DEFAULT_EXTERNAL_SYMP_FORMAT

# Set the external representations
# The external representation may be changed via the phase_rep_format,
# symp_rep_format, or pauli_rep_format methods properties. Any method changing
# these values must make sure to update the others
# Symplectic representation formats are: 'XZ', 'XZY', 'ZX', 'YZX'
# Phase representation formats are: 'i', '-i', 'is', '-is'
# See [ref] for details on the different representation formats
# TODO: Include ref for above.

PHASE_REP_FORMATS = ["i", "-i", "is", "-is"]
PHASE_REP_FORMATS_IMI = ["i", "-i"]
PHASE_REP_FORMATS_ISMIS = ["is", "-is"]
SYMP_REP_FORMATS = ["XZ", "XZY", "ZX", "YZX"]
PAULI_REP_FORMATS = [
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
PAULI_REP_FORMATS_SPLIT = {
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
# indiciates the index tha the Pauli's are acting on. Following Qiskit's current internal
# indexing:
#
# -iX \otimes Y \otimes Z -> -iZ0Y1X2
# X \otimes Y \otimes Z \otimes I \otimes I \otimes I  -> Z3Y4X5

PHASE_REGEX = "[\-+]?1?[ij]?"
PAULI_REGEX = "[IXZY]"

# Product Syntax
#
PHASE_REGEX = "[\-+]?1?[ij]?"
PAULI_REGEX = "[IXZY]"

PRODUCT_PATTERN_SYNTAX = re.compile(f"^{PHASE_REGEX}({PAULI_REGEX})+$")
PRODUCT_CAPTURE_BASE_PATTERN = re.compile(f"({PAULI_REGEX})")
PRODUCT_SPLIT_PATTERN = re.compile(f"^({PHASE_REGEX})(({PAULI_REGEX})+$)")

INDEX_PATTERN_SYNTAX = re.compile(f"^{PHASE_REGEX}({PAULI_REGEX}[0-9]+)+$")
INDEX_CAPTURE_BASE_PATTERN = PRODUCT_CAPTURE_BASE_PATTERN
INDEX_SPLIT_PATTERN = re.compile(f"^({PHASE_REGEX})(({PAULI_REGEX}[0-9]+)+$)")

PRODUCT_SYNTAX = 0
INDEX_SYNTAX = 1
DEFAULT_PATTERN_SYNTAX = PRODUCT_PATTERN_SYNTAX
DEFAULT_SYNTAX = 0

PATTERN_LIST = [PRODUCT_PATTERN_SYNTAX, INDEX_PATTERN_SYNTAX]
CAPTURE_LIST = [PRODUCT_CAPTURE_BASE_PATTERN, INDEX_CAPTURE_BASE_PATTERN]
PATTERN = 0
CAPTURE = 1
FORM_LIST = [PATTERN_LIST, CAPTURE_LIST]
SPLIT_LIST = [PRODUCT_SPLIT_PATTERN, INDEX_SPLIT_PATTERN]


def _load_pattern(identifier=0):
    return _base_load_pattern(identifier, PATTERN)


def _load_capture(identifier=0):
    return _base_load_pattern(identifier, CAPTURE)


def _base_load_pattern(identifier, scope):
    return FORM_LIST[scope][identifier]


def _is_pattern(string, syntax=0):
    pattern = _load_pattern(identifier=syntax)
    return bool(pattern.search(string))


def _syntax_type(string):
    """Determines which syntax type the given Pauli string is in:
    producttype: e.g. -iXYXXIIIIIIXXIX with return PRODUCT_SYNTAX
    index type: e.g. -iX2Y3Z45X11 will return INDEX_SYNTAX
    Args:
        string (str): A Pauli string
    Returns:
        integer: PRODUCT_SYNTAX for product type and INDEX_SYNTAX for index type
    """
    if _is_pattern(string, PRODUCT_SYNTAX):
        return PRODUCT_SYNTAX
    elif _is_pattern(string, INDEX_SYNTAX):
        return INDEX_SYNTAX
    else:
        raise QiskitError(f"Unknown Pauli syntax type: {string}")


# Formats lists
def phase_formats():
    """[summary]
    Returns:
        [type]: [description]
    """
    return PHASE_REP_FORMATS


def symp_formats():
    """[summary]
    Returns:
        [type]: [description]
    """
    return SYMP_REP_FORMATS


def pauli_formats():
    """[summary]
    Returns:
        [type]: [description]
    """
    return PAULI_REP_FORMATS


# -------------------------------------------------------------------------------
# Representation Methods and Conversions
# -------------------------------------------------------------------------------


def _split_rep(rep):
    """Split the Pauli representation format into the phase and symplectic representation
    formats
    Args:
        rep (rep): Pauli representation format
    """
    return PAULI_REP_FORMATS_SPLIT[rep]


def convert_phase_exp(
    phase, input_format=INTERNAL_PHASE_REP_FORMAT, output_format=DEFAULT_EXTERNAL_PHASE_REP_FORMAT
):
    """Convert between the different phase exponents of encoded phase
    Phase Representation/Encodings:
        a) ['i' format] :math:`i^r` where :math:`r \in \mathbb{Z}_4`
        b) ['-i' format] :math:`(-i)^r` where :math:`r \in \mathbb{Z}_4`
        c) ['is' format] :math:`i^r (-1)^s` where :math:`r,s \in \mathbb{Z}_2`
        d) ['-is' format] :math:`(-i)^r (-1)^s` where :math:`r,s, \in \mathbb{Z}_2`
    Args:
        phase (numpy.ndarray or list): array of phase exponents
        input_format (str, optional): Format that input phase is encoded in.
            Defaults to INTERNAL_PHASE_REP_FORMAT
        output_format (str, optional): Format that output phase will be encoded in.
            Defaults to DEFAULT_EXTERNAL_PHASE_REP_FORMAT
    Returns:
        np.array: phase exponent encoded into new format
    Raises:
        QiskitError: Phase representation format not supported or invalid
    Examples:
        phase_exp = np.array([1,2,1,3,2])
        convert_phase_exponent(phase_exp,'i','is')
        array([[1, 0],
               [0, 1],
               [1, 0],
               [1, 1],
               [0, 1]])
        phase_exp = np.array([(0,1),(1,1),(1,0),(1,1)])
        convert_phase_exponent(phase_exp,'is','-i')
        array([2, 1, 3, 1])
    Method:
        The conversions are done via precomputed matrices that describe how the various
        indices change when converting between different formats. Tuple formats are
        linearized to use single integer indices for lookup and are then expanded
        at the end if neeed.
        Converstion example:
        :math:`i^r` to :math:`(-i)^s`  is encoded as 0 to 1 via _ENC matrix
        Find which transformation to use via _TRANS[0][1] = 1
        _CN[1]= [0,3,2,1] and so :math:`i^{[0,1,2,3]} = (-i)^{[0,3,2,1]}` etc
        This method works on vectors of encoded phase exponents
    Raises:
        QiskitError
    """
    if input_format == output_format:
        return phase
    if (input_format not in phase_formats()) or (output_format not in phase_formats()):
        raise QiskitError("Phase representation format not supported or invalid")
    if not isinstance(phase, np.ndarray):
        phase = np.asarray(phase)
    # Binary expansion of an index
    _BI = [[0, 0], [0, 1], [1, 0], [1, 1]]
    # Format encoding
    _ENC = {"i": 0, "-i": 1, "is": 2, "-is": 3}
    # Conversion matrix split and compressed into two index matrices
    # Transformation indices
    _TI = [[0, 1, 2, 3], [0, 3, 2, 1], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2], [0, 1, 3, 2]]
    # Index to transformation matrices via (input_format, output_format) pairs
    _TRANS = [[0, 1, 2, 4], [1, 0, 4, 2], [2, 3, 0, 5], [3, 2, 5, 0]]
    # Conversion is done via precalculated tables that are stored in _CN, _DD and _TRANS
    # with the format encoded into _ENC
    input_format = _ENC[input_format]
    output_format = _ENC[output_format]
    # Linearize: Convert pairs to an index if needed

    if input_format > 1:
        # input format is in ['is', '-is']
        linear_phase = 2 * phase.T[0] + phase.T[1]
    else:
        # input format is in ['i', '-i']
        linear_phase = phase
    # Calculate and then apply the transformation matrix
    trans = _TRANS[input_format][output_format]
    partial_encoded_phase = np.asarray([_TI[trans][t] for t in linear_phase])
    # Delinearize: Convert indexes back to pairs if needed
    if output_format > 1:
        encoded_phase = np.asarray([_BI[t] for t in partial_encoded_phase])
    else:
        encoded_phase = partial_encoded_phase
    return encoded_phase


def change_rep(
    phase_exponent,
    y_count=0,
    *,
    input_format=INTERNAL_PAULI_REP_FORMAT,
    output_format=DEFAULT_EXTERNAL_PAULI_REP_FORMAT,
):
    """Convert a phase exponent from input_format representation to output_format
    representation with respect to a Pauli wth y_count number of Y's (XZ's, ZX's).

    Args:
        y_count(numpy.ndarray of int): number of Y (XZ,ZX) factors in Pauli's
        phase_exponent(str): phase exponent to convert
        input_format(str): Pauli format representation of phase input
        output_format(str): Pauli format representation of phase output
    Raises:
        QiskitError
    """

    if input_format not in PAULI_REP_FORMATS:
        raise QiskitError(f"Unsupported Pauli represenations {input_format}")
    if output_format not in PAULI_REP_FORMATS:
        raise QiskitError(f"Unsupported Pauli represenations {output_format}")
    if not isinstance(phase_exponent, (int, np.integer)) and not isinstance(
        phase_exponent, (tuple, np.ndarray, list)
    ):
        raise QiskitError(
            "The phase exponent must be an element in Z_4" " or a pair of elements in GF(2) x GF(2)"
        )
    return _change_rep(phase_exponent, y_count, input_format, output_format)


def _change_rep(phase_exponent, y_count, input_format, output_format):
    """Convert a phase exponent from input_format representation to output_format
    representation. This method is vectorized.

    Args:
        phase_exponent(numpy.ndarray of int or int tuples): phase exponent to convert
        y_count(numpy.ndarray of int): number of Y (XZ,ZX) factors in Pauli's
        input_format(str): Pauli format representation of input
        output_format(str): Pauli format representation of output
    Method:
    """
    # phases change with changing symplectic formats via a multiple of i. This
    # multiple is given by the converter table: S_1->S_2 has multiplier i^converter[S_1][S_2]
    converter = {
        "XZ": {"ZX": 0, "XZ": 2, "XZY": 3, "YZX": 3},
        "ZX": {"ZX": 0, "XZ": 2, "XZY": 3, "YZX": 3},
        "XZY": {"XZY": 0, "YZX": 0, "XZ": 1, "ZX": 1},
        "YZX": {"YZX": 0, "XZY": 0, "XZ": 1, "ZX": 1},
    }

    input_phase_format, input_symp_format = _split_rep(input_format)
    output_phase_format, output_symp_format = _split_rep(output_format)
    multiplier = converter[input_symp_format][output_symp_format]
    phase_exponent = convert_phase_exp(phase_exponent, input_phase_format, output_phase_format)

    def _cal_phase(exp, marker):
        if marker < 2:
            return (marker, exp[1])
        else:
            return (marker % 2, (exp[1] + 1) % 2)

    # lookup table converter between powers of i to powers of (-i)
    _I_TO_NI = [0, 3, 2, 1]

    if output_phase_format == "i":
        phase_exponent = np.mod(phase_exponent + multiplier * y_count, 4)
    elif output_phase_format == "-i":
        multiplier = _I_TO_NI[multiplier]
        phase_exponent = np.mod(phase_exponent + multiplier * y_count, 4)
    elif output_phase_format == "is":
        res = np.mod(phase_exponent[0] + multiplier * y_count, 4)
        phase_exponent = np.asarray(
            [_cal_phase(exp, marker) for exp, marker in zip(phase_exponent, res)]
        )
    else:
        # output_phase_format == '-is'
        multiplier = _I_TO_NI[multiplier]
        res = np.mod(phase_exponent[0] + multiplier * y_count, 4)
        phase_exponent = np.asarray(
            [_cal_phase(exp, marker) for exp, marker in zip(phase_exponent, res)]
        )
    return phase_exponent


# Conversion function between Pauli phases (str, complex(phase), exp)

# phase2exp
# exp2phase
# phase2str
# str2phase
# exp2str
# str2exp

# This method are all vector enabled and have two forms: method and r_method.
# The method versions have all the checks and conversions in them. The raw
# methods (r_) assume correct input formats and generally do very little checking
# or conversions. The same_type parameter is used to allow scalar values to be output
# when scalar values are input


def norm_phase_str(phase_string, same_type=True):
    """Norm phase

    Args:
        phase_string (str): phase tring
        same_type (bool): is same type

    Returns:

    """
    scalar = is_scalar(phase_string) and same_type
    phase_string = np.atleast_1d(phase_string)

    if scalar:
        return squeeze(r_norm_phase_str(phase_string), scalar=scalar)

    return r_norm_phase_str(phase_string)


def r_norm_phase_str(phase_string):
    """TODO: docstring"""
    return np.array(
        [
            string.replace("+", "", 1).replace("1", "", 1).replace("j", "i", 1)
            for string in phase_string
        ]
    )


def is_scalar(obj):
    """Deteremine if an obj is a scalar. This is not fully robust in that it
    does not take into account all possible ways to encode a scalar but it
    does work for lists, tuples, and np.ndarrays
    """
    return (  # pylint: disable=consider-using-ternary
        isinstance(obj, np.ndarray) and obj.shape == ()
    ) or (not isinstance(obj, (list, tuple)))


def squeeze(array_, scalar=False):
    """Squeeze an numpy array with the option to convert a resulting
    0d array to a scalar (scalar=True)

    Examples:

    np.array([1,2,3]) -> np.array([1,2,3])
    np.array(1) -> 1 (scalar = True)
    np.array(1) -> np.array(1) (scalar = False)
    np.array([[1,2,3]]) -> np.array([1,2,3])
    np.array([[[[1,2,3,4,5]]]]) -> np.array([1,2,3,4,5])
    """
    array_ = np.squeeze(array_)
    if array_.shape == () and scalar is True:
        return array_.item()
    else:
        return array_


def is_exp_type(phase_exponent, input_phase_format):
    """Care must be taken when phase_exponents are pairs in order to distinguish
    between [0,1] being a single exponent or two exponents. For pairs phase_exponent
    should be [[0,1]]
    """
    if input_phase_format not in PHASE_REP_FORMATS:
        raise QiskitError(f"Invalid phase exponent format {input_phase_format}")
    if isinstance(phase_exponent, np.ndarray):
        phase_exponent = phase_exponent.tolist()
    if not isinstance(phase_exponent, list):
        phase_exponent = [phase_exponent]
    if input_phase_format in ["i", "-i"]:
        return all((x in [0, 1, 2, 3]) for x in phase_exponent)
    if input_phase_format in ["is", "-is"]:
        return all((item in [[0, 0], [0, 1], [1, 0], [1, 1]]) for item in phase_exponent)

    return False


def phase2exp(
    phase, output_phase_format=DEFAULT_EXTERNAL_PHASE_REP_FORMAT, roundit=True, same_type=True
):
    """Convert an array of phases to an array of phase exponents"""
    if output_phase_format not in PHASE_REP_FORMATS:
        raise QiskitError(f"Invalid phase exponent format {output_phase_format}")

    scalar = is_scalar(phase) and same_type
    phase = np.atleast_1d(phase)

    if roundit:
        phase = phase.round()

    if scalar:
        return squeeze(r_phase2exp(phase, output_phase_format), scalar=scalar)
    else:
        return r_phase2exp(phase, output_phase_format)


def r_phase2exp(phase, output_phase_format):
    """Convert an array of phases to an array of phase exponents"""
    _ENC = {
        "i": {1: 0, 1j: 1, -1: 2, 0 - 1j: 3},
        "-i": {1: 0, 0 - 1j: 1, -1: 2, 1j: 3},
        "is": {1: (0, 0), -1: (0, 1), 1j: (1, 0), 0 - 1j: (1, 1)},
        "-is": {1: (0, 0), -1: (0, 1), 0 - 1j: (1, 0), 1j: (1, 1)},
    }
    try:
        return np.array([_ENC[output_phase_format][i] for i in phase])
    except Exception as exception:
        raise QiskitError("phases must be complex numbers in ``[1, -1j, -1, 1j]") from exception


def exp2phase(phase_exponent, input_phase_format, same_type=True):
    """Convert an array of phase exponents to an array of phases

    Args:
        phase_exponent (numpy.ndarray, phase exponet): phase exponent(s), array or scalar
        input_phase_format (str): phase format of exponents
        same_type (bool, optional): Controls behaviour when phase_exponent is a single 
        phase expoent that is not in a list or array. Defaults to True.

    Raises:
        QiskitError: Invalid phase format provided
        QiskitError: Input phase_exponent is not in provided phase_format

    Returns:
        complex: complex global phase of input_exponent
    """
    if input_phase_format not in PHASE_REP_FORMATS:
        raise QiskitError(f"Invalid phase exponent format {input_phase_format}")

    scalar = is_scalar(phase_exponent) and same_type
    phase_exponent = np.atleast_1d(phase_exponent)

    if not is_exp_type(phase_exponent, input_phase_format):
        raise QiskitError(f"{phase_exponent} not in {input_phase_format} format")

    if scalar:
        return squeeze(r_exp2phase(phase_exponent, input_phase_format), scalar=scalar)
    else:
        return r_exp2phase(phase_exponent, input_phase_format)


def r_exp2phase(phase_exponent, input_phase_format):
    """Convert an array of phase exponents to an array of phases"""
    if input_phase_format == "i":
        return 1j ** phase_exponent
    if input_phase_format == "-i":
        return (-1j) ** phase_exponent
    if input_phase_format == "is":
        trans = phase_exponent.T
        return np.multiply(1j ** trans[0], (-1) ** trans[1])
    if input_phase_format == "-is":
        trans = phase_exponent.T
        return ((-1j) ** trans[0]) * (-1) ** trans[1]
    else:
        raise QiskitError(f"{input_phase_format} format is not supported yet.")


def phase2str(phase, same_type=True):
    """Comvert an array of phases to an array of strings

    Args:
        phase ([type]): [description]

    Raises:
        QiskitError: [description]

    Returns:
        [type]: [description]
    """
    scalar = is_scalar(phase) and same_type
    phase = np.atleast_1d(phase)

    # We only bother to check the first entry is a complex number (inc real numbers)
    if not isinstance(phase[0], numbers.Complex):
        raise QiskitError(f"{phase}: is not of type complex")
    if scalar:
        return squeeze(r_phase2str(phase.round()), scalar=scalar)
    else:
        return r_phase2str(phase.round())


def r_phase2str(phase):
    """Comvert an array of phases to an array of strings

    Args:
        phase ([type]): [description]

    Raises:
        QiskitError: [description]

    Returns:
        [type]: [description]
    """
    _ENC = {1: "", 1j: "i", -1: "-", 0 - 1j: "-i"}
    return np.asarray([_ENC[item] for item in phase])


def str2phase(phase_string, same_type=True):
    """Convert the string representing the phase of a Pauli to a complex type phase
    Assumes that the string representing the phase is of the form [\-+]?1?[ij]?
    Args:
        same_type (bool): same type flag
        phase_string (str): string representating Pauli phase

    Output:
        phase (complex): phase of Pauli
    """
    scalar = is_scalar(phase_string) and same_type

    phase_string = np.atleast_1d(phase_string)
    phase_string = norm_phase_str(phase_string, same_type=False)

    if scalar:
        return squeeze(r_str2phase(phase_string), scalar=scalar)
    else:
        return r_str2phase(phase_string)


def r_str2phase(phase_string):
    """Convert the string representing the phase of a Pauli to a complex type phase
    Assumes that the string representing the phase is normalized.
    Args:
        phase_string (str): string representating Pauli phase

    Output:
        phase (complex): phase of Pauli

    See Also:
    ---------
    phase_from_phase_string: Convert the string representing the phase of a Pauli to a complex type phase

    """
    CONV_ = {"": 1, "i": 1j, "-i": -1j, "-": -1}
    return np.array([CONV_[item] for item in phase_string])


def str2exp(phase_string, output_format=DEFAULT_EXTERNAL_PHASE_REP_FORMAT, same_type=True):
    """Return the phase from a str representaing a Pauli phase"""

    scalar = is_scalar(phase_string) and same_type

    phase_string = np.atleast_1d(phase_string)
    phase_string = norm_phase_str(phase_string, same_type=False)

    if output_format not in PHASE_REP_FORMATS:
        raise QiskitError(f"str2exp {output_format} ")
    if scalar:
        return squeeze(r_str2exp(phase_string, output_format=output_format), scalar=scalar)
    else:
        return r_str2exp(phase_string, output_format=output_format)


def r_str2exp(phase_string, output_format=DEFAULT_EXTERNAL_PHASE_REP_FORMAT):
    """TODO: description

    Args:
        phase_string (str): phase string
        output_format (str): format
    """
    CONV_ = {
        "i": {"": 0, "-i": 3, "-": 2, "i": 1},
        "-i": {"": 0, "-i": 1, "-": 2, "i": 3},
        "is": {"": (0, 0), "-i": (1, 1), "-": (0, 1), "i": (1, 0)},
        "-is": {"": (0, 0), "-i": (1, 0), "-": (0, 1), "i": (1, 1)},
    }

    return np.array([CONV_[output_format][item] for item in phase_string])


# TODO(Drew Vandeth): This method should be deprecated. Use phase_to_exponent instead
# which are more general and handle vectors
# _phase_from_complex(coeff) -> phase2exp(coeff, output_phase_format="-i")


def _phase_from_complex(coeff):
    """Return the phase exponent ('-i' format) from a phase - method not array capable"""
    if np.isclose(coeff, 1):
        return 0
    if np.isclose(coeff, -1j):
        return 1
    if np.isclose(coeff, -1):
        return 2
    if np.isclose(coeff, 1j):
        return 3
    raise QiskitError("Pauli can only be multiplied by 1, -1j, -1, 1j.")


# ---------------------------------------------------------------------
# Label parsing helper functions
# ---------------------------------------------------------------------


def split_pauli(label, syntax=None):
    """Split Pauli label into unsigned group label and coefficient label.
    Uses _syntax_type to check or determine
    which syntax is being used if not provided

    Assumes that phase is expanded and one of 1,i,-i,-1 in some suitable representation
    (e.g. as in  X (->1) jX (->i) +X (->1) etc.)

    Allows for product and index syntaxs with XZY and YZX representations.
    """

    if syntax is None:
        try:
            syntax = _syntax_type(label)
        except QiskitError:
            print("Unknown string format for Pauli")
    return _split_pauli(label, syntax)


def _split_pauli(label, syntax):
    """Split Pauli label into unsigned group label and coefficient label.

    See split_pauli_label for basis information

    """
    result = re.findall(SPLIT_LIST[syntax], label)
    phase = result[0][0]
    pauli = result[0][1]
    return phase, pauli


#  Moved from Pauli
def from_label(label, qubit_order="right-to-left"):
    """Return the symplectic representation of Pauli string.
    Args:
        label (str): the Pauli string label.
        qubit_order (str): the order in which qubits are assigned to Paulis in the string.
            Has not effect on INDEX_SYNTAX as the order is defined in the string.
            Default: "right-to-left"
    Returns:
        BasePauli: the BasePauli corresponding to the label.
    Raises:
        QiskitError: if Pauli string is not valid.
    """
    # TODO: This method works on building two seperate arrays, one
    # for x and one for z. Then it combines them into one. This can be
    # done without splitting them.
    # Determine which syntax is being used and if valid
    syntax = _syntax_type(label)
    # Split string into coefficient and non signed Pauli Operators
    coeff, pauli = _split_pauli(label, syntax=syntax)
    # Convert coefficient to internal phase exponent
    phase = 0 if not coeff else phase2exp(str2phase(coeff))
    if syntax == PRODUCT_SYNTAX:
        num_qubits = len(pauli)
        if qubit_order == "right-to-left":
            indices = list(reversed(range(num_qubits)))
        else:
            indices = list(range(num_qubits))
    elif syntax == INDEX_SYNTAX:
        indices = list(map(int, re.findall("\d+", pauli)))
        num_qubits = max(indices) + 1
        pauli = re.findall(f"{PAULI_REGEX}+", pauli)
    else:
        QiskitError(f"Conversion to label not implemented for syntax type {syntax}")
    array_z = np.zeros((1, num_qubits), dtype=bool)
    array_x = np.zeros((1, num_qubits), dtype=bool)
    array_phase = np.array([phase], dtype=int)
    # Creae a symplectic representation
    for i, char in enumerate(pauli):
        index = indices[i]
        if char == "X":
            array_x[0, index] = True
        elif char == "Z":
            array_z[0, index] = True
        elif char == "Y":
            array_x[0, index] = True
            array_z[0, index] = True
            array_phase += 1

    array = np.hstack((array_x, array_z))
    return array, array_phase % 4


# Utility functions


def indices_to_boolean(indices, dim):
    """Convert an array/list of indices to a boolean vector

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


def boolean_to_indices(booleans):
    """Convert a array of boolean values into a index list of non-zero-elements

    Args:
        booleans (boolean array): An array of boolean values

    Returns:
        numpy.ndarray: Array of indices of non-zero values
    """
    return np.nonzero(np.asarray(booleans))
