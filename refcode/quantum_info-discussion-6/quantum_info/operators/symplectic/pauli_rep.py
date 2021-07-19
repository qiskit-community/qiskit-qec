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
N-qubit Pauli Representation Format Class
"""
# pylint: disable=invalid-name
# pylint: disable=bad-docstring-quotes  # for deprecate_function decorator

import numpy as np

from qiskit.utils.deprecation import deprecate_function
from qiskit.exceptions import QiskitError

class PauliRep():
    r"""N-Qubit Pauli Representation Format Class
    
    """

    #-------------------------------------------------------------------------------
    # Class Variables/States
    #-------------------------------------------------------------------------------

    # Set the interal representations
    # The internal representation cannot be changed by changing this constant
    # This constants are for reference only and do not change the behavio(u)r of
    # the PauliList methods. See [ref] for details on the different representation
    # formats
    # TODO: Include ref for above.

    __INTERNAL_SYMP_REP_FORMAT__ = 'ZX'
    __INTERNAL_PHASE_REP_FORMAT__ = '-i'
    __INTERNAL_PAULI_REP_FORMAT__ = __INTERNAL_PHASE_REP_FORMAT__ \
        + __INTERNAL_SYMP_REP_FORMAT__
    
    __DEFAULT_EXTERNAL_SYMP_FORMAT__ = 'YZX'
    __DEFAULT_EXTERNAL_PHASE_REP_FORMAT__ = '-i'
    

    # Set the external representations
    # The external representation may be changed via the phase_rep_format,
    # symp_rep_format, or pauli_rep_format methods properties. Any method changing
    # these values must make sure to update the others
    # Symplectic representation formats are: 'XZ', 'XZY', 'ZX', 'YZX'
    # Phase representation formats are: 'i', '-i', 'is', '-is'
    # See [ref] for details on the different representation formats
    # TODO: Include ref for above.

    __external_symp_rep_format__ = __DEFAULT_EXTERNAL_SYMP_FORMAT__
    __external_phase_rep_format__ = __DEFAULT_EXTERNAL_PHASE_REP_FORMAT__
    __external_pauli_rep_format__ = __external_phase_rep_format__ \
        + __external_symp_rep_format__

    __PHASE_REP_FORMATS__ = ['i','-i','is','-is']
    __SYMP_REP_FORMATS__ = ['XZ','XZY', 'ZX', 'YZX']
    __PAULI_REP_FORMATS__ =['iXZ','iXZY','iZX','iYZX',
                               '-iXZ','-iXZY','-iZX','-iYZX',
                               'isXZ','isXZY','isZX','isYZX',
                               '-isXZ','-isXZY','-isZX','-isYZX']
    __PAULI_REP_FORMATS_SPLIT__ = {
        'iXZ': ('i','XZ'), 'iXZY':('i','XZY'), 'iZX':('i','ZX'), 'iYZX':('i','YZX'),
        '-iXZ': ('-i','XZ'), '-iXZY':('-i','XZY'), '-iZX':('-i','ZX'),
        '-iYZX':('-i','YZX'), 'isXZ': ('is','XZ'), 'isXZY':('is','XZY'),
        'isZX':('is','ZX'), 'isYZX':('is','YZX'), '-isXZ': ('-is','XZ'),
        '-isXZY':('-is','XZY'), '-isZX':('-is','ZX'), '-isYZX':('-is','YZX')
        }
    
    @classmethod
    def _set_formats(cls, *, phase_format=None, symp_format=None):
        """Update the external phase and symplectic represenations
        
        Calling with no parameters will reset to the default external
        representations

        Args:
            phase_format (str) : Phase format string from __PHASE_REP_FORMATS__. Default = '-i'
            symp_format (str): Symplectic format string from __SYMP_REP_FORMATS__. Default = 'YZX'

        Raises:
            QiskitError: If formats are not implemented
        """
            
        if phase_format is not None:
            if phase_format in PauliRep.__PHASE_REP_FORMATS__:
                PauliRep.__external_phase_rep_format__ = phase_format
            else:
                raise QiskitError("Invalid phase format")
                
        if symp_format is not None:
            if symp_format  in PauliRep.__SYMP_REP_FORMATS__:
                PauliRep.__external_symp_rep_format__ = symp_format
            else:
                raise QiskitError("Invalid symplectic format")
                    
        if phase_format is None:
            phase_format = PauliRep.__DEFAULT_EXTERNAL_PHASE_REP_FORMAT__
        if symp_format is None:
            symp_format = PauliRep.__DEFAULT_EXTERNAL_SYMP_FORMAT__
        
        PauliRep.__external_phase_rep_format__ = phase_format
        PauliRep.__external_symp_rep_format__ = symp_format
        PauliRep.__external_pauli_rep_format__ = phase_format + symp_format

    # ---------------------------------------------------------------------
    # Property methods
    # ---------------------------------------------------------------------

    @classmethod
    def set_formats(cls,*, phase_format=None, symp_format=None, pauli_format=None):
        """Update the external phase and symplectic represenations
        
        Calling with no parameters will reset to the default external
        representations

        You can either set the formats individually via phase_format and
        symp_format or in one combined string via pauli_format

        Args:
            phase_format (str) : Phase format string from __PHASE_REP_FORMATS__. Default = '-i'
            symp_format (str): Symplectic format string from __SYMP_REP_FORMATS__. Default = 'YZX'
            pauli_format (str): Pauli format string from __PAULI_REP_FORMATS__. Default = None

        Raises:
            QiskitError: If formats are not implemented
        """

        if pauli_format is not None:
            if pauli_format in PauliRep.__PAULI_REP_FORMATS__:
                phase_format, symp_format = cls._split_rep(pauli_format)
        PauliRep._set_formats(phase_format=phase_format, symp_format=symp_format)

    # External formats

    @classmethod
    def external_pauli_format(cls, in_format=None):
        """Display or set the external Pauli operator format

        Calling with no parameter will reset to the default external
        Pauli representations

        Args:
            format (str or None): External Pauli format to use. Default is None

        Returns:
            str: External Pauli operator format

        Raises:
            QiskitError: If unsupported format is requested
        """
        if in_format is None:
            return PauliRep.__external_pauli_rep_format__

        if in_format in PauliRep.__PAULI_REP_FORMATS__:
            phase_format, symp_format = cls._split_rep(in_format)
            PauliRep.set_formats(phase_format=phase_format, symp_format=symp_format)
            return PauliRep.__external_pauli_rep_format__
        else:
            raise QiskitError("Invalid Pauli representation format or unsupported format")

    @classmethod
    def external_phase_format(cls, phase_format=None):
        """Display the external phase representation format for Pauli's operators

        Returns:
            str: Phase representation format for Pauli operator

        """
        if phase_format is None:
            return PauliRep.__external_phase_rep_format__

        PauliRep.set_formats(phase_format=phase_format)

    @classmethod
    def external_symp_format(cls, symp_format=None):
        """Display the external symplectic representation format for Pauli's

        Returns:
            str: Symplectic representation format for Pauli operator
        """
        if symp_format is None:
            return PauliRep.__external_symp_rep_format__

        PauliRep.set_formats(symp_format=symp_format)
        return PauliRep.__external_symp_rep_format__

    # Internal formats

    @classmethod
    def internal_phase_format(cls):
        """Display the internal phase format for Pauli Operators
        
        Returns:
            str: Internal phase format
        """
        return PauliRep.__INTERNAL_PHASE_REP_FORMAT__

    @classmethod
    def internal_symp_format(cls):
        """Display the internal symplectic format for Pauli Operators

        Returns:
            str: internal symplectic format
        """
        return PauliRep.__INTERNAL_SYMP_REP_FORMAT__

    @classmethod
    def internal_pauli_format(cls):
        """Display the internal Pauli operator format

        Returns:
            str: Internal Pauli format
        """
        return PauliRep.__INTERNAL_PAULI_REP_FORMAT__

    # Formats lists

    @classmethod
    def phase_formats(cls):
        """Display the available phase formats for Pauli operators

        Returns:
            list(str): available phase formats
        """
        return PauliRep.__PHASE_REP_FORMATS__

    @classmethod
    def symp_formats(cls):
        """Display the available symplectic formats for Pauli operators

        Returns:
            list(str): aailable symplectic formats
        """
        return PauliRep.__SYMP_REP_FORMATS__
    
    @classmethod
    def pauli_formats(cls):
        """Display the available Pauli formats

        Returns:
            list(str): available Pauli formats
        """
        return PauliRep.__PAULI_REP_FORMATS__
    
    #-------------------------------------------------------------------------------
    # Representation Methods and Conversions
    #------------------------------------------------------------------------------- 
    
    @staticmethod
    def _count_y(array_1, array_2):
        """Count the number of I Pauli's"""
        return np.sum(np.logical_and(array_1, array_2), axis=1)
    
    @staticmethod
    def _split_rep(rep):
        """Split the Pauli representation format into the phase and symplectic representation
        formats

        Args:
            rep (rep): Pauli representation format
        """
        return PauliRep.__PAULI_REP_FORMATS_SPLIT__[rep]

    @staticmethod
    def convert_phase_exponent(
            phase,
            input_format=None, 
            output_format=None):
        r"""Convert between the different phase representations/encodings of phase exponents

        Phase Representation/Encodings:

            a) ['i' format] :math:`i^r` where :math:`r \in \mathbb{Z}_4`
            b) ['-i' format] :math:`(-i)^r` where :math:`r \in \mathbb{Z}_4`
            c) ['is' format] :math:`i^r (-1)^s` where :math:`r,s \in \mathbb{Z}_2`
            d) ['-is' format] :math:`(-i)^r (-1)^s` where :math:`r,s, \in \mathbb{Z}_2`

        Args:
            phase (numpy.ndarray or list): array of phase exponents
            input_format (str, optional): Format that input phase is encoded in.
                Defaults to PauliRep.__INTERNAL_PHASE_REP_FORMAT__
            output_format (str, optional): Format that output phase will be encoded in.
                Defaults to PauliRep.__external_phase_rep_format__

        Returns:
            np.array: phase exponent encoded into new format

        Raises:
            QiskitError: Phase representation format not supported or invalid

        Examples:

            phase_exp = np.array([1,2,1,3,2])
            convert_encoded_phase(phase_exp,'i','is')

            array([[1, 0],
                   [0, 1],
                   [1, 0],
                   [1, 1],
                   [0, 1]])

            phase_exp = np.array([(0,1),(1,1),(1,0),(1,1)])
            convert_encoded_phase(phase_exp,'is','-i')

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
        """
        if input_format is None:
            input_format=PauliRep.internal_phase_format()
        if output_format is None:
            output_format=PauliRep.external_phase_format()
        
        if input_format == output_format:
            return phase

        if (input_format not in PauliRep.phase_formats()) \
             or (output_format not in PauliRep.phase_formats()):
            raise QiskitError("Phase representation format not supported or invalid")

        if not isinstance(phase, np.ndarray):
            phase = np.asarray([phase])

        # Binary expansion of an index
        _BI = [[0,0],
               [0,1],
               [1,0],
               [1,1]]

        # Format encoding
        _ENC = {"i":0, "-i":1, "is":2, "-is":3}

        # Conversion matrix split and compressed into two index matrices
        # Transformation indices
        _TI = [[0,1,2,3],
               [0,3,2,1],
               [0,2,1,3],
               [0,2,3,1],
               [0,3,1,2],
               [0,1,3,2]]

        # Index to transformation matrices via (input_format, output_format) pairs
        _TRANS = [[0,1,2,4],
                  [1,0,4,2],
                  [2,3,0,5],
                  [3,2,5,0]]

        # Conversion is done via precalculated tables that are stored in _CN, _DD and _TRANS
        # with the format encoded into _ENC

        input_format = _ENC[input_format]
        output_format = _ENC[output_format]

        # Linearize: Convert pairs to an index if needed
    
        if input_format > 1:
            # input format is in ['is', '-is']
            linear_phase = 2*phase.T[0]+phase.T[1]
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

    @staticmethod
    def change_representation(
            phase,
            y_count=0,
            *,
            input_format=None,
            output_format=None,
            direction='out'):
        """Convert a phase exponent from input_format representation to output_format
        representation with respect to a Pauli wth y_count number of Y's (XZ's, ZX's).

        Setting one of input_format or output_format will override any direction given
        via the direction input.
        
        Args:
            phase(str): phase exponent to convert
            input_format(str): Pauli format representation of phase input. Default = External Rep Format
            output_format(str): Pauli format representation of phase output. Default = Internal Rep Format
            direction(str): 'in' for external to internal and 'out' for internal to external representations. Default 'out'
        """
        if input_format is None and output_format is None:
            if direction == 'in':
                input_format=PauliRep.external_pauli_format()
                output_format=PauliRep.internal_pauli_format()
            elif direction == 'out':
                input_format=PauliRep.internal_pauli_format()
                output_format=PauliRep.external_pauli_format()
            else:
                raise QiskitError("Invalid direction for conversion")

        if input_format is None:
            input_format=PauliRep.external_pauli_format()
        if output_format is None:
            output_format=PauliRep.internal_pauli_format()
        
        if input_format not in PauliRep.__PAULI_REP_FORMATS__:
            raise QiskitError(f"Unsupported Pauli represenations {input_format}")
        if output_format not in PauliRep.__PAULI_REP_FORMATS__:
            raise QiskitError(f"Unsupported Pauli represenations {output_format}")
        if not isinstance(phase, (int, np.integer)) and not isinstance(phase, (tuple, np.ndarray, list)):
            raise QiskitError("The phase exponent must be an element in Z_4"
                " or a pair of elements in F_2 x F_2"
                )

        if not isinstance(y_count, np.ndarray):
            if isinstance(y_count, list):
                y_count = np.asarray(y_count)
            elif isinstance(y_count, (int, np.integer)):
                y_count = np.asarray([y_count])
            else:
                raise QiskitError("y_count not an integer or numpy.array or list of integers")

        return PauliRep._change_representation(phase, y_count, input_format, output_format)

    @staticmethod
    def _change_representation(
            phase,
            y_count,
            input_format, 
            output_format):
        """Convert a phase exponent from input_format representation to output_format
        representation.
        
        Args:
            phase(numpy.ndarray of int or int tuples): phase exponent to convert
            y_count(numpy.ndarray of int): number of Y (XZ,ZX) factors in Pauli's
            input_format(str): Pauli format representation of input
            output_format(str): Pauli format representation of output

        Method:

        """
        input_phase_format, input_symp_format = PauliRep._split_rep(input_format)
        output_phase_format, output_symp_format = PauliRep._split_rep(output_format)

        # phases change with changing symplectic formats via a multiple of i. This
        # multiple is given by the converter table: S_1->S_2 has multiplier i^converter[S_1][S_2]

        # The converter table is used to convert to powers of i. This is flipped to -i if needed.
        
        converter={'XZ':{'XZ':0, 'ZX':2, 'XZY':3, 'YZX':3},
                   'ZX':{'ZX':0, 'XZ':2, 'XZY':1, 'YZX':1},
                  'XZY':{'XZY':0, 'YZX':0, 'XZ':1, 'ZX':1},
                  'YZX':{'YZX':0, 'XZY':0, 'XZ':3, 'ZX':3}}
        #print(input_symp_format, output_symp_format)
        multiplier = converter[input_symp_format][output_symp_format]
        #print(f"multiplier={multiplier}")
        if output_phase_format in ['-i','-is'] and multiplier%2 != 0:
            multiplier = np.mod(multiplier + 2, 4)

        #print(f"multiplier={multiplier}")

        phase = PauliRep.convert_phase_exponent(phase, input_phase_format, output_phase_format)

        #print(phase)

        def _cal_phase(exp, marker):
            if marker < 2:
                return (marker, exp[1])
            else:
                return (marker%2, (exp[1]+1)%2)

        if output_phase_format in ['i','-i']:
            phase = np.mod(phase + multiplier * y_count, 4)
        else:
            res = np.mod(phase.T[0] + multiplier * y_count, 4)
            phase  = np.asarray([_cal_phase(exp, marker) for exp, marker in zip(phase,res)])

        #print(phase)
        return phase

    @staticmethod
    def exponent_to_coeff(phase, input_phase_format):
        """Convert an array of phase exponents to an array of phase coefficients """
        if input_phase_format not in PauliRep.__PHASE_REP_FORMATS__:
            raise QiskitError(f"Invalid phase exponent format {input_phase_format}")
        return PauliRep._exponent_to_coeff(phase, input_phase_format)

    @staticmethod
    def _exponent_to_coeff(phase, input_phase_format):
        """Convert an array of phase exponents to an array of phase coefficients"""
        if isinstance(phase, (int, np.integer)):
            phase = np.asarray([phase])
        if input_phase_format == 'i':
            return 1j**phase
        if input_phase_format == '-i':
            return (0-1j)**phase
        if not isinstance(phase, (list,np.ndarray)):
            raise QiskitError("phase needs to be a list or nump.ndarray for 'is' and '-is' phase formats")
        if isinstance(phase, list):
            phase = np.asarray(phase)
        if phase.ndim !=2:
            raise QiskitError("phase has wrong dimsions for a 'is', '-is' phase vector")
        if input_phase_format == 'is':
            trans = phase.T
            return np.multiply(1j**trans[0],(-1)**trans[1])
        if input_phase_format == '-is':
            trans = phase.T
            return ((0-1j)**trans[0])*(-1)**trans[1]

    @staticmethod
    def coeff_to_exponent(coeff, output_phase_format, roundit=True):
        """Convert an array of phase coefficients to an array of phase exponents"""
        if output_phase_format not in PauliRep.__PHASE_REP_FORMATS__:
            raise QiskitError(f"Invalid phase exponent format {output_phase_format}")
        if not isinstance(coeff, np.ndarray):
            if isinstance(coeff, list):
                coeff = np.asarray(coeff)
            elif isinstance(coeff, (int, float, complex)):
                coeff = np.asarray([coeff])
            else:
                raise QiskitError("Coefficient provided is not a number")
        if roundit:
            return PauliRep._coeff_to_exponent(coeff.round(), output_phase_format)
        else:
            return PauliRep._coeff_to_exponent(coeff, output_phase_format)

    @staticmethod
    def _coeff_to_exponent(coeff, output_phase_format):
        """Convert an array of phase coefficients to an array of phase exponents"""
        _ENC = {'i':{1:0, 1j:1, -1:2, 0-1j:3}, '-i':{1:0, 0-1j:1, -1:2, 1j:3}, 
        'is':{1:(0,0), -1:(0,1), 1j:(1,0), 0-1j:(1,1)}, '-is':{1:(0,0), -1:(0,1), 0-1j:(1,0), 1j:(1,1)}}
        try:
            return np.asarray([_ENC[output_phase_format][i] for i in coeff ])
        except:
            raise QiskitError("coefficients must be complex numbers in ``[1, -1j, -1, 1j]")


    # TODO: This method should be deprecated. Use coeff_to_exponent or _coeff_to_exponent instead
    # which are more general and handle vectors
    @staticmethod
    def _phase_from_complex(coeff):
        """Return the phase from a label - not array capable"""
        if np.isclose(coeff, 1):
            return 0
        if np.isclose(coeff, -1j):
            return 1
        if np.isclose(coeff, -1):
            return 2
        if np.isclose(coeff, 1j):
            return 3
        raise QiskitError("Pauli can only be multiplied by 1, -1j, -1, 1j.")

    @staticmethod
    def _from_array(z, x, phase=0, input_format=None):
        """Convert array data to BasePauli data."""
        
        if input_format is None:
            input_format = PauliRep.external_pauli_format()
        
        if isinstance(z, np.ndarray) and z.dtype == bool:
            array_z = z
        else:
            array_z = np.asarray(z, dtype=bool)
        if array_z.ndim == 1:
            array_z = array_z.reshape((1, array_z.size))
        elif array_z.ndim != 2:
            raise QiskitError("Invalid Pauli z vector shape.")

        if isinstance(x, np.ndarray) and x.dtype == bool:
            array_x = x
        else:
            array_x = np.asarray(x, dtype=bool)
        if array_x.ndim == 1:
            array_x = array_x.reshape((1, array_x.size))
        elif array_x.ndim != 2:
            raise QiskitError("Invalid Pauli x vector shape.")

        if array_z.shape != array_x.shape:
            raise QiskitError("z and x vectors are different size.")

        # Convert group phase convention to internal ZX-phase conversion.
        # Internal Pauli representation is '-iZX' format
        # External Pauli representation is '-iYZX' format and is assumed as input
        # when referencing is not possible for the phase.

        base_phase = PauliRep.change_representation(
                phase,
                y_count=PauliRep._count_y(array_z, array_x),
                input_format=input_format,
                output_format=PauliRep.internal_pauli_format())

        return array_z, array_x, base_phase

    @staticmethod
    def _to_matrix(z, x, phase=0, group_phase=False, sparse=False):
        """Return the matrix matrix from symplectic representation.

        The Pauli is defined as :math:`P = (-i)^{phase + z.x} * Z^z.x^x`
        where ``array = [x, z]``.

        Args:
            z (array): The symplectic representation z vector.
            x (array): The symplectic representation x vector.
            phase (int): Pauli phase.
            group_phase (bool): Optional. If True use group-phase convention
                                instead of BasePauli ZX-phase convention.
                                (default: False).
            sparse (bool): Optional. Of True return a sparse CSR matrix,
                           otherwise return a dense Numpy array
                           (default: False).

        Returns:
            array: if sparse=False.
            csr_matrix: if sparse=True.
        """
        num_qubits = z.size

        # Convert to zx_phase
        if group_phase:
            phase += np.sum(x & z)
            phase %= 4

        dim = 2**num_qubits
        twos_array = 1 << np.arange(num_qubits)
        x_indices = np.asarray(x).dot(twos_array)
        z_indices = np.asarray(z).dot(twos_array)

        indptr = np.arange(dim + 1, dtype=np.uint)
        indices = indptr ^ x_indices
        if phase:
            coeff = (-1j)**phase
        else:
            coeff = 1
        data = np.array([coeff * (-1) ** (bin(i).count('1') % 2)
                         for i in z_indices & indptr])
        if sparse:
            # Return sparse matrix
            from scipy.sparse import csr_matrix
            return csr_matrix((data, indices, indptr), shape=(dim, dim),
                              dtype=complex)

        # Build dense matrix using csr format
        mat = np.zeros((dim, dim), dtype=complex)
        for i in range(dim):
            mat[i][indices[indptr[i]:indptr[i + 1]]] = data[indptr[i]:indptr[i + 1]]
        return mat

    @staticmethod
    def _to_label(z, x, phase, group_phase=False,
                  full_group=True, return_phase=False):
        """Return the label string for a Pauli.

        Args:
            z (array): The symplectic representation z vector.
            x (array): The symplectic representation x vector.
            phase (int): Pauli phase.
            group_phase (bool): Optional. If True use group-phase convention
                                instead of BasePauli ZX-phase convention.
                                (default: False).
            full_group (bool): If True return the Pauli label from the full Pauli group
                including complex coefficient from [1, -1, 1j, -1j]. If
                False return the unsigned Pauli label with coefficient 1
                (default: True).
            return_phase (bool): If True return the adjusted phase for the coefficient
                of the returned Pauli label. This can be used even if
                ``full_group=False``.

        Returns:
            str: the Pauli label from the full Pauli group (if ``full_group=True``) or
                from the unsigned Pauli group (if ``full_group=False``).
            Tuple[str, int]: if ``return_phase=True`` returns a tuple of the Pauli
                            label (from either the full or unsigned Pauli group) and
                            the phase ``q`` for the coefficient :math:`(-i)^(q + x.z)`
                            for the label from the full Pauli group.
        """
        num_qubits = z.size
        coeff_labels = {0: '', 1: '-i', 2: '-', 3: 'i'}
        label = ''
        for i in range(num_qubits):
            if not z[num_qubits - 1 - i]:
                if not x[num_qubits - 1 - i]:
                    label += 'I'
                else:
                    label += 'X'
            elif not x[num_qubits - 1 - i]:
                label += 'Z'
            else:
                label += 'Y'
                if not group_phase:
                    phase -= 1
        phase %= 4
        if phase and full_group:
            label = coeff_labels[phase] + label
        if return_phase:
            return label, phase
        return label

