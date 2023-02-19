# This code is part of Qiskit.
#
# (C) Copyright IBM 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Analysis extensions"""

import logging  # for logging!

logger = logging.getLogger(__name__)

# Load extensions if available and set appriate indicator flags

try:
    from qiskit_qec.analysis._c_analysis import _CErrorPropagator  # pylint: disable=unused-import

    C_ERROR_PROPAGATOR = True
except ImportError as import_error:
    logger.exception(  # pylint: disable=logging-fstring-interpolation
        f"from qiskit_qec.analysis._c_analysis import _CErrorPropagator \
            failed, raising {import_error}"
    )
    C_ERROR_PROPAGATOR = False

try:
    from qiskit_qec.analysis._c_analysis import _CFaultEnumerator  # pylint: disable=unused-import

    C_FAULT_ENUMERATOR = True
except ImportError as import_error:
    logger.exception(  # pylint: disable=logging-fstring-interpolation
        f"from qiskit_qec.analysis._c_analysis import _CFaultEnumerator \
            failed, raising {import_error}"
    )
    C_FAULT_ENUMERATOR = False

try:
    from qiskit_qec.analysis._c_analysis import _CFaultSampler  # pylint: disable=unused-import

    C_FAULT_SAMPLER = True
except ImportError as import_error:
    logger.exception(  # pylint: disable=logging-fstring-interpolation
        f"from qiskit_qec.analysis._c_analysis import _CFaultSampler failed, \
            raising {import_error}"
    )
    C_FAULT_SAMPLER = False

try:
    from qiskit_qec.analysis._c_analysis import _c_minimum_distance  # pylint: disable=unused-import

    C_MIN_DISTANCE = True
except ImportError as import_error:
    logger.exception(  # pylint: disable=logging-fstring-interpolation
        f"from qiskit_qec.analysis._c_analysis import _c_minimum_distance \
            failed, raising {import_error}"
    )
    C_MIN_DISTANCE = False

try:
    # pylint: disable=unused-import
    from qiskit_qec.analysis._c_analysis import _c_minimum_distance_by_tests

    C_MIN_DISTANCE_BY_TESTS = True
except ImportError as import_error:
    logger.exception(  # pylint: disable=logging-fstring-interpolation
        f"from qiskit_qec.analysis._c_analysis import _c_minimum_distance_by_tests \
            failed, raising {import_error}"
    )
    C_MIN_DISTANCE_BY_TESTS = False

try:
    from qiskit_qec.analysis._c_analysis import _c_rank  # pylint: disable=unused-import

    C_RANK = True
except ImportError as import_error:
    logger.exception(  # pylint: disable=logging-fstring-interpolation
        f"from qiskit_qec.analysis._c_analysis import _c_rank failed, \
            raising {import_error}"
    )
    C_RANK = False

try:
    from qiskit_qec.analysis._c_analysis import _c_isotropic  # pylint: disable=unused-import

    C_ISOTROPIC = True
except ImportError as import_error:
    logger.exception(  # pylint: disable=logging-fstring-interpolation
        f"from qiskit_qec.analysis._c_analysis import _c_isotropic failed, \
            raising {import_error}"
    )
    C_ISOTROPIC = False
