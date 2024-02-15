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
# pylint: disable=unused-import

"""Code circuit extensions"""

import logging  # for logging!

logger = logging.getLogger(__name__)

# Load extensions if available and set appriate indicator flags

try:
    from qiskit_qec.analysis._c_circuits import _c_check_nodes

    C_CHECK_NODES = True
except ImportError as import_error:
    logger.exception(  # pylint: disable=logging-fstring-interpolation
        f"from qiskit_qec.analysis._c_circuits import _c_check_nodes \
            failed, raising {import_error}"
    )
    C_CHECK_NODES = False

try:
    from qiskit_qec.analysis._c_circuits import _c_is_cluster_neutral

    C_IS_CLUSTER_NEUTRAL = True
except ImportError as import_error:
    logger.exception(  # pylint: disable=logging-fstring-interpolation
        f"from qiskit_qec.analysis._c_circuits import _c_is_cluster_neutral \
            failed, raising {import_error}"
    )
    C_IS_CLUSTER_NEUTRAL = False
