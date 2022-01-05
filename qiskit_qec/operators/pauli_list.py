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

import numpy as np
from scipy import sparse

from qiskit.exceptions import QiskitError
from qiskit.quantum_info.operators.symplectic.pauli_table import PauliTable
from qiskit.quantum_info.operators.symplectic.stabilizer_table import StabilizerTable

import qiskit_qec.utils.pauli_rep as rep

from qiskit_qec.operators.base_pauli import BasePauli
from qiskit_qec.operators.pauli import Pauli

class PauliList(BasePauli):

    # Set the max number of qubits * paulis before string truncation
    __truncate__ = 2000

    def __init__(self, pdata, stype='numpy', phase_exponent=0, input_qubit_order="right-to-left"):

        if isinstance(pdata, BasePauli):
            matrix = pdata.matrix
            p_ext  = pdata.phase_exponent
            stype  = pdata.stype
        elif isinstance(pdata, StabilizerTable):
            # Conversion from legacy StabilizerTable
            raise QiskitError("Input method not yet implemented")
        elif isinstance(pdata, PauliTable):
            # Conversion from legacy PauliTable
            raise QiskitError("Input method not yet implemented")
        elif isinstance(pdata, np.ndarray):
            if isinstance(pdata[0], (bool, np.bool_, int, np.integer)):
                matrix = np.atleast_2d(pdata)
                p_ext = np.asarray(phase_exponent)
            elif isinstance(pdata[0], (list, tuple, np.ndarray)): 
                # TODO: (conditional part and approach below) This should be done better (lazy at the moment)
                # TODO: Check for stype setting
                matrix = np.asarray(pdata)
                if matrix.shape[0] == 1:
                    p_ext = phase_exponent
                elif phase_exponent==0:
                    p_ext = np.zeros(matrix.shape[0], dtype=int)
                else:
                    p_ext = np.asarray(phase_exponent)
            else:
                # Conversion as iterable of Paulis
                matrix, p_ext, ex_stype = self._from_paulis(pdata, input_qubit_order)
                if ex_stype is not None:
                    stype = ex_stype
        else:
            # Conversion as iterable of Paulis
            matrix, p_ext, ex_stype = self._from_paulis(pdata, input_qubit_order)
        
            if ex_stype is not None:
                stype = ex_stype

        super().__init__(matrix, p_ext, stype)     
        self.paulis = self.init_paulis()

    def init_paulis(self):
        """Create the initial list of Pauli objects that 
        reference the primary symplectic matrix. This is used
        as a time memory tradeoff to significantly speed of access
        of individual Pauli operators.

        Returns:
            List[Pauli]: List of Paulis that reference the primary symplectic matrix
            representation
        """
        #for i in range(self.matrix.shape[0]):
        #    print(f"self.matrix[i]={self.matrix[i]}")
        paulis = [Pauli(self.sm.atleast_2d(self.matrix[i]),phase_exponent=self.phase_exponent[i]) for i in range(self.matrix.shape[0])]
        return paulis

    def __getitem__(self, slc):
        """General getitem method. Probably not that fast. Use othe getitem methods
        for fast results. getitem, ...

        TODO: This needs to be expanded with checks and other selection capabilities

        Args:
            slc (int): [description]

        Returns:
            [type]: [description]
        """
        # if slc is a single index use this form
        return self.paulis[slc]
        # if slc is a slice of something else use other forms
        # To be added.

    def getitem(self, i):
        """Specific getitem to get a single Pauli using the precomputed Paulis.

        Args:
            i (int): Index of required Pauli

        Returns:
            Pauli: Pauli at index i
        """
        return PauliList(self.matrix[i])

    #def __repr__(self):
    #    return np.array2string(self.matrix)

    def __repr__(self):
        """Display representation."""
        return self._truncated_str(True)

    def __str__(self):
        """Print representation."""
        return self._truncated_str(False)

    def _truncated_str(self, show_class):
        stop = self.num_paulis
        if self.__truncate__:
            max_paulis = self.__truncate__ // self.num_qubits
            if self.num_paulis > max_paulis:
                stop = max_paulis
        labels = [str(self[i]) for i in range(stop)]
        prefix = "PauliList(" if show_class else ""
        tail = ")" if show_class else ""
        if stop != self.num_paulis:
            suffix = ", ...]" + tail
        else:
            suffix = "]" + tail
        list_str = np.array2string(
            np.array(labels), threshold=stop + 1, separator=", ", prefix=prefix, suffix=suffix
        )
        return prefix + list_str[:-1] + suffix


    def _from_paulis(self, data, stype=None, input_qubit_order="right-to-left"):
        """Construct a PauliList from a list of Pauli data.

        Args:
            data (iterable): list of Pauli data.
            stype (): Override which stype is used. Default is None

        Returns:
            PauliList: the constructed PauliList.

        Raises:
            QiskitError: If the input list is empty or contains invalid
            Pauli strings.
        """
        if not isinstance(data, (list, tuple, set, np.ndarray)):
            data = [data]
        num_paulis = len(data)
        if num_paulis == 0:
            raise QiskitError("Input Pauli list is empty.")
        paulis = []
        for i in data:
            if not isinstance(i, Pauli):
                #TODO: Should check if matrix and then use appropriate tools
                # Will not work at the moment with non-numpy arrays
                if isinstance(i, (np.ndarray,list)):
                    paulis.append(Pauli(np.atleast_2d(i), input_qubit_order=input_qubit_order))
                else:
                    paulis.append(Pauli(i, input_qubit_order=input_qubit_order))
            else:
                paulis.append(i)

        # Determine the num_qubits (max number of the individual Paulis)
        qubit_count = [pauli.num_qubits for pauli in paulis]
        num_qubits = max(qubit_count)

        #TODO: This will not work when using other type of matrices. Needs
        # to be generalized
        matrix = np.zeros((num_paulis, 2*num_qubits), dtype=bool)
        phase_exponent = np.zeros(num_paulis, dtype=int)
    
        # Note that the next step implicitly broadcasts an input array from two different shapes
        # This is take extra time but is probably okay for most applications
        stype = paulis[0].stype
        for i, pauli in enumerate(paulis):
            if pauli.stype != stype:
                # TODO: Handle situations when a list of Paulis is passed 
                # using more than one internal symplectic matrix representation
                raise QiskitError("Different stype lists are not yet implemented")

            matrix[i][:2*qubit_count[i]] = pauli.matrix
            #print(f"pauli.phase_exponent={pauli.phase_exponent}")
            phase_exponent[i] = pauli.phase_exponent
        return matrix, phase_exponent, stype