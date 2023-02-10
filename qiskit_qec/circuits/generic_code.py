# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# pylint: disable=invalid-name

"""Creates a CodeCircuit object for a pre-existing circuit."""

from qiskit_aer import AerSimulator


class GenericCodeCircuit:
    """GenericCodeCircuit Class"""

    def __init__(self, circuits, ese_list, ese_info=None):
        """
        circuits (list, tuple or QuantumCircuit): `QuantumCircuit` or collection thereof.
        The output strings of these circuits should have the same format.
        ese_list: List of tuples, with each tuple containing a set of indices for the
        output strings (numbered right to left), such that the parity of the corresponding
        output bits should have a fixed value.
        """

        # store circuit in a dictionary with the index as the key
        if isinstance(circuits, (list, tuple, set)):
            self.circuit = {}
            for j, circuit in enumerate(circuits):
                self.circuit[j] = circuit
        else:
            self.circuit = {0: circuits}
        # choose the first to be the base circuit
        self.base = 0

        self.ese_list = ese_list
        self.ese_info = ese_info
        self._get_parities()

    def _get_ese_string(self, string, ese):
        """
        For the given string, determine the parities of all error sensitive events.
        """
        ese_string = ""
        for index in ese:
            ese_string += string[-1 - index]
        return ese_string

    def _get_parities(self):
        """
        Simulates circuit to determine what values the parities should ideally be
        for each circuit.
        """
        qasm_sim = AerSimulator()
        parities = {}
        for basis, qc in self.circuit.items():
            parities[basis] = []
            counts = qasm_sim.run(qc).result().get_counts()
            string = list(counts.keys())[0]
            for ese in self.ese_list:
                ese_string = self._get_ese_string(string, ese)
                parity = ese_string.count("1") % 2
                parities[basis].append(parity)
        self.parities = parities

    def string2nodes(self, string, **kwargs):
        """
        Args:
            string (string): Results string to convert.
            kwargs (dict): Additional keyword arguments.
                index (int): Index of circuit for which nodes are
                calculated.
        Returns:
            dict: List of nodes corresponding to to the non-trivial
            elements in the string.
        """
        index = kwargs.get("index")
        if index is None:
            index = self.base
        nodes = []
        for e, ese in enumerate(self.ese_list):
            ese_string = self._get_ese_string(string, ese)
            if ese_string.count("1") % 2 != self.parities[index][e]:
                if self.ese_info:
                    nodes.append(
                        {
                            "element": e,
                            "is_boundary": self.ese_info["is_boundary"],
                            "qubits": self.ese_info["qubits"],
                            "time": self.ese_info["time"],
                        }
                    )
                else:
                    nodes.append({"element": e, "is_boundary": False, "qubits": [], "time": 0})
        return nodes
