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

from qiskit.exceptions import QiskitError

from qiskit.circuit import Instruction, QuantumCircuit
from qiskit.quantum_info.operators.scalar_op import ScalarOp

import qiskit_qec.utils.pauli_rep as rep
from qiskit_qec.operators.base_pauli import BasePauli

class Pauli(BasePauli):
    # Pauli (x,z) lookup table to string
    pltb_str = {(0,0):"I", (1,0):"X", (0,1):"Z", (1,1):"Y"}
    # Pauli (x,z) lookup table to int
    pltb_int = {(0,0):0, (1,0):1, (0,1):2, (1,1):3}

    def __init__(
        self, 
        data, 
        *,
        x=None,
        z=None,
        phase_exponent=None,
        stype='numpy', 
        label=None, 
        input_qubit_order="right-to-left"):

        #print("Pauli")
        #print(f"data={data}")
        #print(f"phase_exponent={phase_exponent}")
        # TODO: Update to include other possible smatrix types
        # Should call a SMatrix function rather than list them
        # directly here
        if isinstance(data, np.ndarray):
            #print("Choice 1")
            matrix         = np.atleast_2d(data)
            phase_exponent = phase_exponent
        elif isinstance(data, BasePauli):
            #print("Choice 2")
            matrix         = data.matrix
            phase_exponent = data.phase_exponent
        elif isinstance(data, tuple):
            #print("Choice 3")
            if len(data) not in [1, 2]:
                if len(data) == 3: # DEPRECATED
                    # Convert (z,x),p) to (matrix, p)
                    # Note old specific order
                    raise QiskitError("Not yet implemented")
                raise QiskitError(
                    "Invalid input tuple for Pauli, input tuple must be"
                    " `(matrix, phase)` or `(matrix,)`"
                )
            raise QiskitError("Not yet implemented")
            matrix, phase_exponent = self._from_array(*data)
        elif isinstance(data, str):
            matrix, phase_exponent = rep.from_label(data, qubit_order = input_qubit_order)
        elif isinstance(data, ScalarOp):
            raise QiskitError("Not yet implemented")
            matrix, phase_exponent = self._from_scalar_op(data)
        elif isinstance(data, (QuantumCircuit, Instruction)):
            raise QiskitError("Not yet implemented")
            matrix, phase_exponent = self._from_circuit(data)
        elif x is not None:  
            #print("Choice 4")
            if z is None: # DEPRECATED
                # Using old Pauli initialization with positional args instead of kwargs
                z = data
                raise QiskitError("Not yet implemented")
                matrix, phase_exponent = self._from_array_deprecated(z, x)
            else:
                raise QiskitError("Not yet implemented")
                matrix, phase_exponent = self._from_array(x=x, z=z, phase_exponent=phase_exponent)

        elif label is not None:  # DEPRECATED
            raise QiskitError("Not yet implemented")
            matrix, phase_exponent = self._from_label_deprecated(label)
        else:
            raise QiskitError("Invalid input data for Pauli.")
        #print("End Choice")
        #print("to BasePauli")
        #print(f"matrix={matrix}")
        #print(f"shape={matrix.shape}")
        #print(f"phase_exponent={phase_exponent}")
        # Initialize BasePauli
        if matrix.shape[0] != 1:
            raise QiskitError("Input is not a single Pauli")      

        
        super().__init__(matrix, phase_exponent, stype=stype)
        self.vlist = self.matrix[0].tolist()

    def __getitem__(self, i):
        """Get Pauli for qubit i

        Args:
            i (int): index of qubit

        Returns:
            Pauli: Pauli acting on qubit i
        """
        # This is expensive ~6us compared to ~100ns for a Pauli get from a PauliList. Better
        # than the ~25mu that the existing Pauli takes in quantum_info
        # Does allow referencing but limited to single indexing. 
        # 
        # There are pros and cons associated with how we represent the symplectic
        # matrix. If you represent this matrix as a X and Z matrix then more slicing options 
        # are possible here that would allow referening but only to a limited extent.
        # 
        # Linear slicing with steps could be down with full vector representations. Linear slicing
        # with steps is chosing a contigous block of qubits to select and steping in that block.
        # for expample say [l:k:s]. To do this you can reshape the full vector from
        # (1, 2*num_qubits) to (2, num_qubits) and then slicing that matrix with [l:k:s]. This 
        # gives some matrix of size (2, selected_num_qubits) with the 0th row representing the 
        # X part and the 1st row representing the Z part. Yoiu cannot reshape this as it will distroy
        # the refering.
        # 
        # This could be used by a second form of the vector (a matrix) would need to be stored and
        # interpreted appropriately.
        # 
        # If refercing is not necessary then it is easy to extract the necessary value in what
        # ever form is desired.
        # 
        # Suggest that __getitem__ be as general as possible to allow all sorts of possible selctions
        # some of which will and some will not provide referencing. Then have some specific fast functions
        # that check less and can do less.
        # 
        # With relative speeds
        #
        # __getitem__  6us (will get slower as more checks and capability added)
        # getitem 6us
        # fast_getitem_int 1.5us
        # fast_getitem_str 1.5us
        # etc 
        # If you just reture the tuple (x,z) value then you get
        # getxz 
        #
        # One can get faster results if for Paulis we also store a list version
        # 
        # vlist_getitem_raw 530ns
        # vlist_getitem_int 580ns
        # 
        # These results are fore small dimensions (and using nump arrays). For different situtations 
        # different methods may should be used.
    

        return self.getitem(i)

    def getitem(self, i):
        """Get Pauli for qubit i

        Args:
            i (int): index of qubit

        Returns:
            Pauli: Pauli acting on qubit i
        """
        return Pauli(self.matrix[0][i:self.matrix.shape[1]:self.num_qubits])

    def fast_getitem_str(self, i):
        """Get Pauli for qubit i

        Args:
            i (int): index of qubit

        Returns:
            str: Streing representing the Pauli acting on qubit i, (0,0):"I", (1,0):"X", (0,1):"Z", (1,1):"Y"
        """
        return Pauli.pltb_str[(self.matrix[0][i], self.matrix[0][i+self.num_qubits])]

    def fast_getitem_int(self, i):
        """Get Pauli for qubit i

        Args:
            i (int): index of qubit

        Returns:
            int: Integer representing the Pauli acting on qubit i, (I:0, X:1, Z:2, Y:4)
        """
        
        return Pauli.pltb_int[(self.matrix[0][i], self.matrix[0][i+self.num_qubits])]   

    def getitem_raw(self,i):
        """Get Pauli symplectic (x,z)-tuple for qubit i

        Args:
            i (int): index of qubit

        Returns:
            [type]: Symplectic (x,z)-tuple for the Pauli on qubit i
        """
        return (self.matrix[0][i], self.matrix[0][i+self.num_qubits])


    def vlist_getitem_raw(self, i):
        """Get Pauli symplectic (x,z)-tuple for qubit i. Requires extra storage 
        for the list representation of the Pauli

        Args:
            i (int): index of qubit

        Returns:
            [type]: Symplectic (x,z)-tuple for the Pauli on qubit i
        """
        return (self.vlist[i], self.vlist[i+self.num_qubits])

    def vlist_getitem_int(self, i):
        """Get Pauli for qubit i using a stored list of the Pauli.

        Args:
            i (int): index of qubit

        Returns:
            int: Integer representing the Pauli acting on qubit i, (I:0, X:1, Z:2, Y:4)
        """
        return Pauli.pltb_int[(self.vlist[i], self.vlist[i+self.num_qubits])]
        
    def __repr__(self):
        return np.array2string(self.matrix)