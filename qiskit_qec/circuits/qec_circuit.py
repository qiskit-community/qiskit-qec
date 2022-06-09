# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2019.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# pylint: disable=invalid-name

"""Classes that create and manage circuits for codes."""
from typing import List, Tuple

from qiskit import QuantumCircuit


class CSSCircuit:
    """
    Create and manage circuits for generic CSS codes.
    """

    def __init__(
        self,
        n: int,
        css_x_gauge_ops: List[Tuple[int]],
        css_x_stabilizer_ops: List[Tuple[int]],
        css_x_boundary: List[int],
        css_z_gauge_ops: List[Tuple[int]],
        css_z_stabilizer_ops: List[Tuple[int]],
        css_z_boundary: List[int],
        basis: str,
        round_schedule: str,
        blocks: int,
        resets: bool,
        delay: float,
    ):
        """

        Args:
            n : number of code qubits
            css_x_gauge_ops : list of supports of X gauge operators
            css_x_stabilizer_ops : list of supports of X stabilizers
            css_x_boundary : list of qubits along the X-type boundary
            css_z_gauge_ops : list of supports of Z gauge operators
            css_z_stabilizer_ops : list of supports of Z stabilizers
            css_x_boundary : list of qubits along the Z-type boundary
            basis : initializaton and measurement basis ("x" or "z")
            round_schedule : gauge measurements in each block
            blocks : number of measurement blocks
            resets : Whether to include a reset gate after mid-circuit measurements.
            delay: Time (in dt) to delay after mid-circuit measurements (and reset).
        """
        self.n = n
        self.css_x_gauge_ops = css_x_gauge_ops
        self.css_x_stabilizer_ops = css_x_stabilizer_ops
        self.css_x_boundary = css_x_boundary
        self.css_z_gauge_ops = css_z_gauge_ops
        self.css_z_stabilizer_ops = css_z_stabilizer_ops
        self.css_z_boundary = css_z_boundary
        self.basis = basis
        self.round_schedule = round_schedule
        self.blocks = blocks
        self.resets = resets
        self.delay = delay

        self.circuit = QuantumCircuit()

    def x(self, barrier=False):
        """
        Applies a logical x to the circuit.

        Args:
            barrier (bool): Boolean denoting whether to include a barrier at
                the end.
        """
        pass

    def z(self, barrier=False):
        """
        Applies a logical x to the circuit.

        Args:
            barrier (bool): Boolean denoting whether to include a barrier at
                the end.
        """
        pass

    def syndrome_measurement(self, final: bool = False, barrier: bool = False, delay: int = 0):
        """
        Application of a syndrome measurement round.

        Args:
            final (bool): Whether to disregard the reset (if applicable) due to this
            being the final syndrome measurement round.
            barrier (bool): Whether to include a barrier at the end.
            delay (float): Time (in dt) to delay after mid-circuit measurements (and reset).
        """
        pass

    def readout(self):
        """
        Readout of all code qubits, which corresponds to a logical measurement
        as well as allowing for a measurement of the syndrome to be inferred.
        """
        pass

    def string2nodes(self, string, all_logicals=False):
        """
        Convert output string from running the circuit into a set of nodes.
        Args:
            string (string): Results string to convert.
            all_logicals (bool): Whether to include logical nodes
            irrespective of value.
        Returns:
            dict: List of nodes corresponding to to the non-trivial
            elements in the string.
        """
        pass

    def string2raw_logicals(self, string):
        """
        Extracts raw logical measurement outcomes from output string.
        Args:
            string (string): Results string from which to extract logicals
        Returns:
            list: Raw values for logical operators that correspond to nodes.
        """
        pass
