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
"""Module fo Pauli List"""
import numbers
from collections import defaultdict
from typing import Iterable, List, Tuple, Union

import numpy as np
import rustworkx as rx
from qiskit.exceptions import QiskitError
from qiskit.quantum_info.operators.custom_iterator import CustomIterator
from qiskit.quantum_info.operators.mixins import GroupMixin, LinearMixin

from qiskit_qec.operators.base_pauli import BasePauli
from qiskit_qec.operators.pauli import Pauli
from qiskit_qec.utils import pauli_rep


class PauliList(BasePauli, LinearMixin, GroupMixin):
    """`PauliList` inherits from `BasePauli`"""

    # Set the max number of qubits * paulis before string truncation
    _truncate__ = 2000

    def __init__(
        self,
        data: Union[BasePauli, np.ndarray, Tuple[np.ndarray], Iterable, None] = None,
        phase_exp: Union[int, np.ndarray, None] = None,
        *,
        input_pauli_encoding: str = BasePauli.EXTERNAL_PAULI_ENCODING,
        input_qubit_order: str = "right-to-left",
        tuple_order: str = "zx",
    ) -> None:
        """Inits a PauliList

        Args:
            data (str): List of Pauli Operators. Ex: 'IIXXZ'
            phase_exp (int, optional): i**phase_exp. Defaults to 0.
            input_qubit_order (str, optional): Order to read pdata. Defaults to "right-to-left".

        Raises:
            QiskitError: Something went wrong.
        """
        if data is None:
            matrix = np.empty(shape=(0, 0), dtype=np.bool_)
            phase_exp = np.empty(shape=(0,), dtype=np.int8)
        elif isinstance(data, BasePauli):
            matrix = data.matrix
            phase_exp = data._phase_exp

        elif isinstance(data, np.ndarray):
            if data.size == 0:
                matrix = np.empty(shape=(0, 0), dtype=np.bool_)
                phase_exp = np.empty(shape=(0,), dtype=np.int8)
            elif isinstance(data[0], str):
                matrix, phase_exp = self._from_paulis(data, input_qubit_order)
            else:
                if phase_exp is None:
                    phase_exp = 0
                matrix, phase_exp = pauli_rep.from_array(
                    data, phase_exp, input_pauli_encoding=input_pauli_encoding
                )
        elif isinstance(data, tuple):
            if len(data) not in [2, 3]:
                raise QiskitError(
                    "Invalid input tuple for PauliList, input tuple must be `(z, x, phase)` or `(z, x)`"
                )
            if tuple_order not in ["zx", "xz"]:
                raise QiskitError(f"`tuple_order` {tuple_order} not valid")
            if len(data) == 3:
                if tuple_order == "zx":
                    matrix, phase_exp = pauli_rep.from_split_array(
                        data[1], data[0], data[2], input_pauli_encoding=input_pauli_encoding
                    )
                else:
                    matrix, phase_exp = pauli_rep.from_split_array(
                        *data, input_pauli_encoding=input_pauli_encoding
                    )
            elif tuple_order == "zx":
                matrix, phase_exp = pauli_rep.from_split_array(
                    data[1], data[0], 0, input_pauli_encoding=input_pauli_encoding
                )
            else:
                matrix, phase_exp = pauli_rep.from_array(
                    *data, 0, input_pauli_encoding=input_pauli_encoding
                )
        else:
            # Conversion as iterable of Paulis
            if len(data) == 0:
                matrix = np.empty(shape=(0, 0), dtype=np.bool_)
                phase_exp = np.empty(shape=(0,), dtype=np.int8)
            else:
                matrix, phase_exp = self._from_paulis(data, input_qubit_order)

        super().__init__(matrix, phase_exp)

        self.paulis = [
            Pauli(self.matrix[i], phase_exp=self._phase_exp[i]) for i in range(self.matrix.shape[0])
        ]

    # ---------------------------------------------------------------------
    # Init Methods
    # ---------------------------------------------------------------------

    @staticmethod
    def _from_paulis(data: Iterable, input_qubit_order: str = "right-to-left"):
        """Construct a PauliList from a list of Pauli data.

        Args:
            data (iterable): list of Pauli data.

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
        for pauli_data in data:
            if not isinstance(pauli_data, Pauli):
                paulis.append(Pauli(pauli_data, input_qubit_order=input_qubit_order))
            else:
                paulis.append(pauli_data)

        # Determine the num_qubits (max number of the individual Paulis)
        qubit_count = [pauli.num_qubits for pauli in paulis]
        num_qubits = max(qubit_count)

        matrix = np.zeros((num_paulis, num_qubits << 1), dtype=bool)
        phase_exp = np.zeros(num_paulis, dtype=int)

        for index, pauli in enumerate(paulis):
            matrix[index][: pauli.num_qubits] = pauli.matrix[0][: pauli.num_qubits]
            matrix[index][num_qubits : num_qubits + pauli.num_qubits] = pauli.matrix[0][
                pauli.num_qubits :
            ]
            phase_exp[index] = pauli._phase_exp[0]
        return matrix, phase_exp

    # ---------------------------------------------------------------------
    # Property Methods
    # ---------------------------------------------------------------------

    @property
    def phase(self):
        """Return the phase vector of the PauliList.

        Note: This is different from the quantum_info phase property which
        instead returns the phase_exp
        """
        # Convert internal exponent frmt to complex number representation
        return pauli_rep.exp2cpx(
            pauli_rep.change_pauli_encoding(
                self._phase_exp, self.num_y, output_pauli_encoding=BasePauli.EXTERNAL_PAULI_ENCODING
            ),
            input_encoding=BasePauli.EXTERNAL_PHASE_ENCODING,
        )

    @phase.setter
    def phase(self, phase):
        """Set the phase vector of the PauliList

        Args:
            phase (numpy.ndarray or complex numbers): Array of phases,
                phases must be one of [1,-1, 1j, -1j]
        """
        phase_exp = pauli_rep.cpx2exp(phase, output_encoding=pauli_rep.INTERNAL_PHASE_ENCODING)
        self._phase_exp[:] = phase_exp

    @property
    def shape(self):
        """The full shape of the :meth:`array`"""
        return self._num_paulis, self.num_qubits

    @property
    def size(self):
        """The number of Pauli rows in the table."""
        return self._num_paulis

    @property
    def num_paulis(self):
        """Returns the number of Pauli's in List"""
        return self._num_paulis

    @property
    def phase_exp(self):
        """Return the phase exponent vector of the PauliList"""
        return pauli_rep._change_pauli_encoding(
            self._phase_exp,
            self.num_y,
            pauli_rep.INTERNAL_PAULI_ENCODING,
            BasePauli.EXTERNAL_PAULI_ENCODING,
        )

    @phase_exp.setter
    def phase_exp(self, phase_exp, input_phase_encoding=BasePauli.EXTERNAL_PHASE_ENCODING):
        """Set the phase exponent vector of the PauliList. Note that this method
        converts the phase exponents directly and does not take into account the
        number of Y paulis in the representation.

        Args:
            phase_exp (_type_): _description_
            input_phase_encoding (_type_, optional): _description_. Defaults to
                BasePauli.EXTERNAL_PHASE_ENCODING.
        """
        self._phase_exp[:] = pauli_rep.exp2exp(
            phase_exp,
            input_encoding=input_phase_encoding,
            output_encoding=pauli_rep.INTERNAL_PHASE_ENCODING,
        )

    @property
    def settings(self):
        """Return settings."""
        return {"data": self.to_labels()}

    # ---------------------------------------------------------------------
    # Magic Methods and related methods
    # ---------------------------------------------------------------------

    def __getitem__(self, index):
        """Return a view of the PauliList."""
        # Returns a view of specified rows of the PauliList
        # This supports all slicing operations the underlying array supports.
        if isinstance(index, tuple):
            if len(index) == 1:
                index = index[0]
            elif len(index) > 2:
                raise IndexError(f"Invalid PauliList index {index}")

        # Row-only indexing
        if isinstance(index, (int, np.integer)):
            # Single Pauli
            return self.paulis[index]
        elif isinstance(index, (slice, list, np.ndarray)):
            # Sub-Table view
            return PauliList(BasePauli(self.matrix[index], self._phase_exp[index]))

        # Row and Qubit indexing
        return PauliList(self.matrix[index])

    def getaslist(self, slc: Union[numbers.Integral, slice]) -> List["Pauli"]:
        """_summary_

        Returns:
            _type_: _description_
        """
        return self.paulis[slc]

    def __setitem__(self, index, value):
        """Update PauliList."""
        if isinstance(index, tuple):
            if len(index) == 1:
                index = index[0]
            elif len(index) > 2:
                raise IndexError(f"Invalid PauliList index {index}")

        # Modify specified rows of the PauliList
        if not isinstance(value, PauliList):
            value = PauliList(value)

        self.matrix[index] = value.matrix

        if not isinstance(index, tuple):
            # Row-only indexing
            self._phase_exp[index] = value._phase_exp[0]
        else:
            # Row and Qubit indexing
            self._phase_exp[index[0]] += value._phase_exp
            self._phase_exp %= 4

    def __repr__(self):
        """Display representation."""
        return self._truncated_str(True)

    def __str__(self):
        """Print representation."""
        return self._truncated_str(False)

    def _truncated_str(self, show_class):
        if self.num_paulis == 0:
            if show_class:
                return "PauliList([])"
            else:
                return "[]"
        stop = self._num_paulis
        if self._truncate__:
            max_paulis = self._truncate__ // self.num_qubits
            if self._num_paulis > max_paulis:
                stop = max_paulis
        labels = [str(self[i]) for i in range(stop)]
        prefix = "PauliList(" if show_class else ""
        tail = ")" if show_class else ""
        if stop != self._num_paulis:
            suffix = ", ...]" + tail
        else:
            suffix = "]" + tail
        list_str = np.array2string(
            np.array(labels), threshold=stop + 1, separator=", ", prefix=prefix, suffix=suffix
        )
        return prefix + list_str[:-1] + suffix

    def __array__(self, dtype=None):
        """Convert to numpy array"""
        # pylint: disable=unused-argument
        shape = (len(self),) + 2 * (2**self.num_qubits,)
        ret = np.zeros(shape, dtype=complex)
        for i, mat in enumerate(self.matrix_iter()):
            ret[i] = mat
        return ret

    def __eq__(self, other):
        """Entrywise comparison of Pauli equality."""
        if not isinstance(other, PauliList):
            other = PauliList(other)
        if not isinstance(other, BasePauli):
            return False
        return self._eq(other)

    def __len__(self):
        """Return the number of Pauli rows in the table."""
        return self._num_paulis

    # ----
    #
    # ----

    def delete(self, ind, qubit=False):
        """Return a copy with Pauli rows deleted from table.

        When deleting qubits the qubit index is the same as the
        column index of the underlying :attr:`X` and :attr:`Z` arrays.

        Args:
            ind (int or list): index(es) to delete.
            qubit (bool): if True delete qubit columns, otherwise delete
                          Pauli rows (Default: False).

        Returns:
            PauliList: the resulting table with the entries removed.

        Raises:
            QiskitError: if ind is out of bounds for the array size or
                         number of qubits.

        Note: Update this method to work with other encodings (assumes Y type encoding in phase_exp)
        """
        if isinstance(ind, int):
            ind = [ind]

        # Row deletion
        if not qubit:
            if max(ind) >= len(self):
                raise QiskitError(
                    "fIndices {ind} are not all less than the size \
                     of the PauliList ({len(self)})"
                )
            matrix = np.delete(self.matrix, ind, axis=0)
            phase_exp = np.delete(self.phase_exp, ind)

            return PauliList(matrix, phase_exp=phase_exp)

        # Column (qubit) deletion
        if max(ind) >= self.num_qubits:
            raise QiskitError(
                "Indices {ind} are not all less than the number of\
                 qubits in the PauliList ({self.num_qubits})"
            )
        ind = ind + [item + self.num_qubits for item in ind]
        matrix = np.delete(self.matrix, ind, axis=1)
        # Use self.phase, not self._phase as deleting qubits can change the
        # ZX phase convention
        return PauliList(
            matrix, phase_exp=self.phase_exp, input_pauli_encoding=BasePauli.EXTERNAL_PAULI_ENCODING
        )

    def insert(self, ind, value, qubit=False):
        """Insert Pauli's into the table.

        When inserting qubits the qubit index is the same as the
        column index of the underlying :attr:`X` and :attr:`Z` arrays.

        Args:
            ind (int): index to insert at.
            value (PauliList): values to insert.
            qubit (bool): if True delete qubit columns, otherwise delete
                          Pauli rows (Default: False).

        Returns:
            PauliList: the resulting table with the entries inserted.

        Raises:
            QiskitError: if the insertion index is invalid.

        Note: Update this method to work with other encodings (assumes Y type encoding in phase_exp)
        """
        if not isinstance(ind, int):
            raise QiskitError("Insert index must be an integer.")

        if not isinstance(value, PauliList):
            value = PauliList(value)
        # Row insertion
        size = self._num_paulis
        if not qubit:
            if ind > size:
                raise QiskitError(
                    f"Index {ind} is larger than the number of rows in the" " PauliList ({size})."
                )
            matrix = np.insert(self.matrix, ind, value.matrix, axis=0)

            base_phase_exp = np.insert(self._phase_exp, ind, value._phase_exp)

            return PauliList(BasePauli(matrix, phase_exp=base_phase_exp))

        # Column insertion
        if ind > self.num_qubits:
            raise QiskitError(
                f"Index {ind} is greater than number of qubits"
                " in the PauliList ({self.num_qubits})"
            )
        if len(value) == 1:
            # Pad blocks to correct size
            value_x = np.vstack(size * [value._x])
            value_z = np.vstack(size * [value._z])
            value_phase_exp = np.repeat(value.phase_exp, size)
        elif len(value) == size:
            #  Blocks are already correct size
            value_x = value._x
            value_z = value._z
            value_phase_exp = value.phase_exp
        else:
            # Blocks are incorrect size
            raise QiskitError(
                f"Input PauliList must have a single row, or \
                 the same number of rows as the Pauli Table \
                 ({size})."
            )

        # Build new array by blocks
        z = np.hstack([self._z[:, :ind], value_z, self._z[:, ind:]])
        x = np.hstack([self._x[:, :ind], value_x, self._x[:, ind:]])
        phase_exp = self.phase_exp + value_phase_exp

        return PauliList((z, x, phase_exp), input_pauli_encoding=BasePauli.EXTERNAL_PAULI_ENCODING)

    def argsort(self, weight=False, phase=False):
        """Return indices for sorting the rows of the table.

        The default sort method is lexicographic sorting by qubit number.
        By using the `weight` kwarg the output can additionally be sorted
        by the number of non-identity terms in the Pauli, where the set of
        all Pauli's of a given weight are still ordered lexicographically.

        Args:
            weight (bool): Optionally sort by weight if True (Default: False).
            phase (bool): Optionally sort by phase before weight or order
                          (Default: False).

        Returns:
            array: the indices for sorting the table.
        """
        # Get order of each Pauli using
        # I => 0, X => 1, Y => 2, Z => 3
        x = self._x
        z = self._z
        order = 1 * (x & ~z) + 2 * (x & z) + 3 * (~x & z)
        phases_exp = self.phase_exp
        # Optionally get the weight of Pauli
        # This is the number of non identity terms
        if weight:
            weights = np.sum(x | z, axis=1)

        # To preserve ordering between successive sorts we
        # are use the 'stable' sort method
        indices = np.arange(self._num_paulis)

        # Initial sort by phases
        sort_inds = phases_exp.argsort(kind="stable")
        indices = indices[sort_inds]
        order = order[sort_inds]
        if phase:
            phases_exp = phases_exp[sort_inds]
        if weight:
            weights = weights[sort_inds]

        # Sort by order
        for i in range(self.num_qubits):
            sort_inds = order[:, i].argsort(kind="stable")
            order = order[sort_inds]
            indices = indices[sort_inds]
            if weight:
                weights = weights[sort_inds]
            if phase:
                phases_exp = phases_exp[sort_inds]

        # If using weights we implement a sort by total number
        # of non-identity Paulis
        if weight:
            sort_inds = weights.argsort(kind="stable")
            indices = indices[sort_inds]
            phases_exp = phases_exp[sort_inds]

        # If sorting by phase we perform a final sort by the phase value
        # of each pauli
        if phase:
            indices = indices[phases_exp.argsort(kind="stable")]
        return indices

    def sort(self, weight=False, phase=False):
        """Sort the rows of the table.

        The default sort method is lexicographic sorting by qubit number.
        By using the `weight` kwarg the output can additionally be sorted
        by the number of non-identity terms in the Pauli, where the set of
        all Pauli's of a given weight are still ordered lexicographically.

        **Example**

        Consider sorting all a random ordering of all 2-qubit Paulis

        .. jupyter-execute::

            from numpy.random import shuffle
            from qiskit.quantum_info.operators import PauliList

            # 2-qubit labels
            labels = ['II', 'IX', 'IY', 'IZ', 'XI', 'XX', 'XY', 'XZ',
                      'YI', 'YX', 'YY', 'YZ', 'ZI', 'ZX', 'ZY', 'ZZ']
            # Shuffle Labels
            shuffle(labels)
            pt = PauliList(labels)
            print('Initial Ordering')
            print(pt)

            # Lexicographic Ordering
            srt = pt.sort()
            print('Lexicographically sorted')
            print(srt)

            # Weight Ordering
            srt = pt.sort(weight=True)
            print('Weight sorted')
            print(srt)

        Args:
            weight (bool): optionally sort by weight if True (Default: False).
            phase (bool): Optionally sort by phase before weight or order
                          (Default: False).

        Returns:
            PauliList: a sorted copy of the original table.
        """
        return self[self.argsort(weight=weight, phase=phase)]

    def unique(self, return_index=False, return_counts=False):
        """Return unique Paulis from the table.

        **Example**

        .. jupyter-execute::

            from qiskit.quantum_info.operators import PauliList

            pt = PauliList(['X', 'Y', '-X', 'I', 'I', 'Z', 'X', 'iZ'])
            unique = pt.unique()
            print(unique)

        Args:
            return_index (bool): If True, also return the indices that
                                 result in the unique array.
                                 (Default: False)
            return_counts (bool): If True, also return the number of times
                                  each unique item appears in the table.

        Returns:
            PauliList: unique
                the table of the unique rows.

            unique_indices: np.ndarray, optional
                The indices of the first occurrences of the unique values in
                the original array. Only provided if ``return_index`` is True.

            unique_counts: np.array, optional
                The number of times each of the unique values comes up in the
                original array. Only provided if ``return_counts`` is True.

        # Check is all phases used are correct (_phase_exp versus phase_exp)
        """
        # Check if we need to stack the phase array
        if np.any(self._phase_exp != self._phase_exp[0]):
            # Create a single array of Pauli's and phases for calling np.unique on
            # so that we treat different phased Pauli's as unique
            array = np.hstack(
                [self._z, self._x, self.phase_exp.reshape((self.phase_exp.shape[0], 1))]
            )
        else:
            # All Pauli's have the same phase so we only need to sort the array
            array = np.hstack([self._z, self._x])

        # Get indexes of unique entries
        if return_counts:
            _, index, counts = np.unique(array, return_index=True, return_counts=True, axis=0)
        else:
            _, index = np.unique(array, return_index=True, axis=0)

        # Sort the index so we return unique rows in the original array order
        sort_inds = index.argsort()
        index = index[sort_inds]
        unique = PauliList((self._z[index], self._x[index], self.phase_exp[index]))

        # Concatinate return tuples
        ret = (unique,)
        if return_index:
            ret += (index,)
        if return_counts:
            ret += (counts[sort_inds],)
        if len(ret) == 1:
            return ret[0]
        return ret

    # ---------------------------------------------------------------------
    # BaseOperator methods
    # ---------------------------------------------------------------------

    def tensor(self, other):
        """Return the tensor product with each Pauli in the list.

        Args:
            other (PauliList): another PauliList.

        Returns:
            PauliList: the list of tensor product Paulis.

        Raises:
            QiskitError: if other cannot be converted to a PauliList, does
                         not have either 1 or the same number of Paulis as
                         the current list.
        """
        if not isinstance(other, PauliList):
            other = PauliList(other)
        return PauliList(super().tensor(other))

    def expand(self, other):
        """Return the expand product of each Pauli in the list.

        Args:
            other (PauliList): another PauliList.

        Returns:
            PauliList: the list of tensor product Paulis.

        Raises:
            QiskitError: if other cannot be converted to a PauliList, does
                         not have either 1 or the same number of Paulis as
                         the current list.
        """
        if not isinstance(other, PauliList):
            other = PauliList(other)
        if len(other) not in [1, len(self)]:
            raise QiskitError(
                "Incompatible PauliLists. Other list must "
                "have either 1 or the same number of Paulis."
            )
        return PauliList(super().expand(other))

    def compose(self, other, qargs=None, front=False, inplace=False):
        """Return the composition self∘other for each Pauli in the list.

        Args:
            other (PauliList): another PauliList.
            qargs (None or list): qubits to apply dot product on (Default: None).
            front (bool): If True use `dot` composition method [default: False].
            inplace (bool): If True update in-place (default: False).

        Returns:
            PauliList: the list of composed Paulis.

        Raises:
            QiskitError: if other cannot be converted to a PauliList, does
                         not have either 1 or the same number of Paulis as
                         the current list, or has the wrong number of qubits
                         for the specified qargs.
        """
        if qargs is None:
            qargs = getattr(other, "qargs", None)
        if not isinstance(other, PauliList):
            other = PauliList(other)
        if len(other) not in [1, len(self)]:
            raise QiskitError(
                "Incompatible PauliLists. Other list must "
                "have either 1 or the same number of Paulis."
            )
        return PauliList(super().compose(other, qargs=qargs, front=front, inplace=inplace))

    # pylint: disable=arguments-differ
    def dot(self, other, qargs=None, inplace=False):
        """Return the composition other∘self for each Pauli in the list.

        Args:
            other (PauliList): another PauliList.
            qargs (None or list): qubits to apply dot product on (Default: None).
            inplace (bool): If True update in-place (default: False).

        Returns:
            PauliList: the list of composed Paulis.

        Raises:
            QiskitError: if other cannot be converted to a PauliList, does
                         not have either 1 or the same number of Paulis as
                         the current list, or has the wrong number of qubits
                         for the specified qargs.
        """
        return self.compose(other, qargs=qargs, front=True, inplace=inplace)

    def _add(self, other, qargs=None):
        """Append two PauliLists.

        If ``qargs`` are specified the other operator will be added
        assuming it is identity on all other subsystems.

        Args:
            other (PauliList): another table.
            qargs (None or list): optional subsystems to add on
                                  (Default: None)

        Returns:
            PauliList: the concatenated list self + other.
        """

        if self.num_qubits == 0:
            return other.copy()
        if other.num_qubits == 0:
            return self.copy()

        if qargs is None:
            qargs = getattr(other, "qargs", None)

        if not isinstance(other, PauliList):
            other = PauliList(other)

        self._op_shape._validate_add(other._op_shape, qargs)

        base_phase_exp = np.hstack((self._phase_exp, other._phase_exp))

        if qargs is None or (sorted(qargs) == qargs and len(qargs) == self.num_qubits):
            matrix = np.vstack((self.matrix, other.matrix))

        else:
            # Pad other with identity and then add
            padded = BasePauli(
                np.zeros((other.size, self.num_qubits << 1), dtype=bool),
                phase_exp=np.zeros(other.size, dtype=int),
            )
            padded = padded.compose(other, qargs=qargs, inplace=True)
            matrix = np.vstack([self.matrix, padded.matrix])

        return PauliList(BasePauli(matrix, phase_exp=base_phase_exp))

    def _multiply(self, phase):
        """Multiply each Pauli in the list by a phase.

        Args:
            other (complex or array): a complex number in [1, -1j, -1, 1j]

        Returns:
            PauliList: the list of Paulis other * self.

        Raises:
            QiskitError: if the phase is not in the set [1, -1j, -1, 1j].
        """
        return PauliList(super()._multiply(phase))

    def conjugate(self):
        """Return the conjugate of each Pauli in the list."""
        return PauliList(super().conjugate())

    def transpose(self):
        """Return the transpose of each Pauli in the list."""
        return PauliList(super().transpose())

    def adjoint(self):
        """Return the adjoint of each Pauli in the list."""
        return PauliList(super().adjoint())

    def inverse(self):
        """Return the inverse of each Pauli in the list."""
        return PauliList(super().adjoint())

    # ---------------------------------------------------------------------
    # Utility methods
    # ---------------------------------------------------------------------

    def commutes(self, other, qargs=None):
        """Return True for each Pauli that commutes with other.

        Args:
            other (PauliList): another PauliList operator.
            qargs (list): qubits to apply dot product on (default: None).

        Returns:
            bool: True if Pauli's commute, False if they anti-commute.
        """
        if qargs is None:
            qargs = getattr(other, "qargs", None)
        if not isinstance(other, BasePauli):
            other = PauliList(other)
        return super().commutes(other, qargs=qargs)

    def anticommutes(self, other, qargs=None):
        """Return True if other Pauli that anticommutes with other.

        Args:
            other (PauliList): another PauliList operator.
            qargs (list): qubits to apply dot product on (default: None).

        Returns:
            bool: True if Pauli's anticommute, False if they commute.
        """
        return np.logical_not(self.commutes(other, qargs=qargs))

    def commutes_with_all(self, other):
        """Return indexes of rows that commute other.

        If other is a multi-row Pauli list the returned vector indexes rows
        of the current PauliList that commute with *all* Pauli's in other.
        If no rows satisfy the condition the returned array will be empty.

        Args:
            other (PauliList): a single Pauli or multi-row PauliList.

        Returns:
            array: index array of the commuting rows.
        """
        return self._commutes_with_all(other)

    def anticommutes_with_all(self, other):
        """Return indexes of rows that commute other.

        If other is a multi-row Pauli list the returned vector indexes rows
        of the current PauliList that anti-commute with *all* Pauli's in other.
        If no rows satisfy the condition the returned array will be empty.

        Args:
            other (PauliList): a single Pauli or multi-row PauliList.

        Returns:
            array: index array of the anti-commuting rows.
        """
        return self._commutes_with_all(other, anti=True)

    def _commutes_with_all(self, other, anti=False):
        """Return row indexes that commute with all rows in another PauliList.

        Args:
            other (PauliList): a PauliList.
            anti (bool): if True return rows that anti-commute, otherwise
                         return rows that commute (Default: False).

        Returns:
            array: index array of commuting or anti-commuting row.
        """
        # Update/compare speed with new symplectic methods
        if not isinstance(other, PauliList):
            other = PauliList(other)
        comms = self.commutes(other[0])
        (inds,) = np.where(comms == int(not anti))
        for pauli in other[1:]:
            comms = self[inds].commutes(pauli)
            (new_inds,) = np.where(comms == int(not anti))
            if new_inds.size == 0:
                # No commuting rows
                return new_inds
            inds = inds[new_inds]
        return inds

    def evolve(self, other, qargs=None, frame="h"):
        r"""Evolve the Pauli by a Clifford.

        This returns the Pauli :math:`P^\prime = C.P.C^\dagger`.

        By choosing the parameter frame='s', this function returns the Schrödinger evolution of the Pauli
        :math:`P^\prime = C.P.C^\dagger`. This option yields a faster calculation.

        Args:
            other (Pauli or Clifford or QuantumCircuit): The Clifford operator to evolve by.
            qargs (list): a list of qubits to apply the Clifford to.
            frame (string): 'h' for Heisenberg or 's' for Schrödinger framework.

        Returns:
            Pauli: the Pauli :math:`C.P.C^\dagger`.

        Raises:
            QiskitError: if the Clifford number of qubits and qargs don't match.
        """
        from qiskit.circuit import Instruction, QuantumCircuit

        if qargs is None:
            qargs = getattr(other, "qargs", None)

        if not isinstance(other, (BasePauli, Instruction, QuantumCircuit)):
            # Convert to a PauliList
            other = PauliList(other)

        return PauliList(super().evolve(other, qargs=qargs, frame=frame))

    def to_labels(self, array=False):
        r"""Convert a PauliList to a list Pauli string labels.

        For large PauliLists converting using the ``array=True``
        kwarg will be more efficient since it allocates memory for
        the full Numpy array of labels in advance.

        .. list-table:: Pauli Representations
            :header-rows: 1

            * - Label
              - Symplectic
              - Matrix
            * - ``"I"``
              - :math:`[0, 0]`
              - :math:`\begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}`
            * - ``"X"``
              - :math:`[1, 0]`
              - :math:`\begin{bmatrix} 0 & 1 \\ 1 & 0  \end{bmatrix}`
            * - ``"Y"``
              - :math:`[1, 1]`
              - :math:`\begin{bmatrix} 0 & -i \\ i & 0  \end{bmatrix}`
            * - ``"Z"``
              - :math:`[0, 1]`
              - :math:`\begin{bmatrix} 1 & 0 \\ 0 & -1  \end{bmatrix}`

        Args:
            array (bool): return a Numpy array if True, otherwise
                          return a list (Default: False).

        Returns:
            list or array: The rows of the PauliList in label form.
        """
        if array:
            return self.to_labels()
        else:
            return self.to_label().tolist()

    def to_matrix(self, sparse=False, array=False):
        r"""Convert to a list or array of Pauli matrices.

        For large PauliLists converting using the ``array=True``
        kwarg will be more efficient since it allocates memory a full
        rank-3 Numpy array of matrices in advance.

        .. list-table:: Pauli Representations
            :header-rows: 1

            * - Label
              - Symplectic
              - Matrix
            * - ``"I"``
              - :math:`[0, 0]`
              - :math:`\begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}`
            * - ``"X"``
              - :math:`[1, 0]`
              - :math:`\begin{bmatrix} 0 & 1 \\ 1 & 0  \end{bmatrix}`
            * - ``"Y"``
              - :math:`[1, 1]`
              - :math:`\begin{bmatrix} 0 & -i \\ i & 0  \end{bmatrix}`
            * - ``"Z"``
              - :math:`[0, 1]`
              - :math:`\begin{bmatrix} 1 & 0 \\ 0 & -1  \end{bmatrix}`

        Args:
            sparse (bool): if True return sparse CSR matrices, otherwise
                           return dense Numpy arrays (Default: False).
            array (bool): return as rank-3 numpy array if True, otherwise
                          return a list of Numpy arrays (Default: False).

        Returns:
            list: A list of dense Pauli matrices if `array=False` and `sparse=False`.
            list: A list of sparse Pauli matrices if `array=False` and `sparse=True`.
            array: A dense rank-3 array of Pauli matrices if `array=True`.
        """
        if not array:
            # We return a list of Numpy array matrices
            return list(self.matrix_iter(sparse=sparse))
        # For efficiency we also allow returning a single rank-3
        # array where first index is the Pauli row, and second two
        # indices are the matrix indices
        dim = 2**self.num_qubits
        ret = np.zeros((self.size, dim, dim), dtype=complex)
        iterator = self.matrix_iter(sparse=sparse)
        for i in range(self.size):
            ret[i] = next(iterator)
        return ret

    # ---------------------------------------------------------------------
    # Custom Iterators
    # ---------------------------------------------------------------------

    def label_iter(self):
        """Return a label representation iterator.

        This is a lazy iterator that converts each row into the string
        label only as it is used. To convert the entire table to labels use
        the :meth:`to_labels` method.

        Returns:
            LabelIterator: label iterator object for the PauliList.
        """

        class LabelIterator(CustomIterator):
            """Label representation iteration and item access."""

            def __repr__(self):
                return f"<PauliList_label_iterator at {hex(id(self))}>"

            def __getitem__(self, key):
                return self.obj[key].to_label()

        return LabelIterator(self)

    def matrix_iter(self, sparse=False):
        """Return a matrix representation iterator.

        This is a lazy iterator that converts each row into the Pauli matrix
        representation only as it is used. To convert the entire table to
        matrices use the :meth:`to_matrix` method.

        Args:
            sparse (bool): optionally return sparse CSR matrices if True,
                           otherwise return Numpy array matrices
                           (Default: False)

        Returns:
            MatrixIterator: matrix iterator object for the PauliList.
        """

        class MatrixIterator(CustomIterator):
            """Matrix representation iteration and item access."""

            def __repr__(self):
                return f"<PauliList_matrix_iterator at {hex(id(self))}>"

            def __getitem__(self, key):
                return self.obj[key].to_matrix(sparse=sparse)

        return MatrixIterator(self)

    # ---------------------------------------------------------------------
    # Class methods
    # ---------------------------------------------------------------------

    @classmethod
    def from_symplectic(cls, z, x, phase_exp=0):
        """Construct a PauliList from a symplectic data.

        Args:
            z (np.ndarray): 2D boolean Numpy array.
            x (np.ndarray): 2D boolean Numpy array.
            phase_exp (np.ndarray or None): Optional, 1D integer array from Z_4.

        Returns:
            PauliList: the constructed PauliList.

        Note: Initialization this way will copy matrices and not reference them.

        TODO: Fix this method to be more general and not in old form only
            (i.e. include matrix inputs ...)
        """

        matrix, phase_exp = pauli_rep.from_split_array(
            x, z, phase_exp, input_pauli_encoding=BasePauli.EXTERNAL_PAULI_ENCODING
        )
        return PauliList(BasePauli(matrix, phase_exp=phase_exp))

    def _noncommutation_graph(self):
        """Create an edge list representing the qubit-wise non-commutation graph.

        An edge (i, j) is present if i and j are not commutable.

        Returns:
            List[Tuple(int,int)]: A list of pairs of indices of the PauliList that
                are not commutable.
        """
        # convert a Pauli operator into int vector where {I: 0, X: 2, Y: 3, Z: 1}
        mat1 = np.array(
            [op.z + 2 * op.x for op in self],
            dtype=np.int8,
        )
        mat2 = mat1[:, None]
        # mat3[i, j] is True if i and j are qubit-wise commutable
        mat3 = (((mat1 * mat2) * (mat1 - mat2)) == 0).all(axis=2)
        # convert into list where tuple elements are qubit-wise non-commuting operators
        return list(zip(*np.where(np.triu(np.logical_not(mat3), k=1))))

    def group_qubit_wise_commuting(self):
        """Partition a PauliList into sets of mutually qubit-wise commuting Pauli strings.

        Returns:
            List[PauliList]: List of PauliLists where each PauliList contains commutable Pauli operators.
        """
        nodes = range(self._num_paulis)
        edges = self._noncommutation_graph()
        graph = rx.PyGraph()
        graph.add_nodes_from(nodes)
        graph.add_edges_from_no_data(edges)
        # Keys in coloring_dict are nodes, values are colors
        coloring_dict = rx.graph_greedy_color(graph)
        groups = defaultdict(list)
        for idx, color in coloring_dict.items():
            groups[color].append(idx)
        return [PauliList([self[i] for i in x]) for x in groups.values()]
