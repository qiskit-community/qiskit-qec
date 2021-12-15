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

import sys

import numpy as np
from scipy import sparse

import qiskit_qec.linear.smatrix_api.sm_numpy as sm_numpy
import qiskit_qec.linear.smatrix_api.sm_sparse as sm_sparse

class SMatrix:
    def __new__(cls, matrix, stype='numpy'):
        """Return a matrix of a specific type depending on the value of stype

        Args:
            matrix ([type]): input matrix to be converted to an stype matrix
            stype (str, optional): type of matrix to return. Defaults to 'numpy'.

        Raises:
            TypeError: An error is raised if an unknown stype is provided.

        Returns:
            stype matrix: a matrix of type stype initialized by the input matrix
        """
        if stype == 'numpy':
            return np.asarray(matrix)
        elif stype == "sparse":
            # Note: Need to check that sparse is referencing than copying
            return sparse.csr_matrix(matrix)
        else:
            raise TypeError("Unknown stype")
            
    def get_methods(stype='numpy'):
        """Return the module that contains the api methods for a given stype matrix

        Args:
            stype (str, optional): type of matrix module api to retrun. Defaults to 'numpy'.

        Raises:
            KeyError: A KeyError will be raised if a needed module that cannot be loaded 
            is not loaded.

        Returns:
            module : A module that includes the api for the stype matrices
        """
        # Assumes that the needed modules are of a specific name sm.stype 
        # 
        # numpy -> sm_numpy
        # sparse -> sm_sparse
        #
        # This needs to be done better
        module_name = f"qiskit_qec.linear.smatrix_api.sm_{stype}"
        try:
            module = sys.modules[module_name]
        except KeyError:
            # TODO: If not loaded -> load via __import__
            raise KeyError(f"Module {module_name} not loaded.")
        return module