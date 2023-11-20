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

"""Class that manage circuits for codes."""
from abc import ABC, abstractmethod


class CodeCircuit(ABC):
    """Abstract class to manage circuits for codes, as well
    as other fault-tolerant circuits.

    A CodeCircuit requires the methods `string2nodes`,
    `check_nodes` and `is_cluster_neutral` in order to
    interface with its `DecodingGraph` and decoders.
    """

    def __init__(self):
        """Initialization of classes that inherent from CodeCircuit can
        be done in various ways, depending on the code or code family
        to be initialized. In all cases, the initialization must define
        the following attributes:

        Args:
            circuit (dict): A dictionary with circuits as values, and
                labels (typically strings) as keys.
            base (string) The label for the above regarded as the base case,
                used in decoding graph generation.
            d (int): Code distance.
            T (int): number of syndrome measurement rounds.
        """
        pass

    @abstractmethod
    def string2nodes(self, string, **kwargs):
        """
        Convert output string from circuits into a set of nodes for
        `DecodingGraph`.
        Args:
            string (string): Results string to convert.
            kwargs (dict): Any additional keyword arguments.
        """
        pass

    @abstractmethod
    def measured_logicals(self):
        """
        Returns a list of logical operators, each expressed as a list of qubits for which
        the parity of the final readouts corresponds to the raw logical readout.
        """
        pass

    @abstractmethod
    def check_nodes(self, nodes, ignore_extra_boundary=False, minimal=False):
        """
        Determines whether a given set of nodes are neutral. If so, also
        determines any additional logical readout qubits that would be
        flipped by the errors creating such a cluster and how many errors
        would be required to make the cluster.
        Args:
            nodes (list): List of nodes, of the type produced by `string2nodes`.
            ignore_extra_boundary (bool): If `True`, undeeded boundary nodes are
            ignored.
            minimal (bool): Whether output should only reflect the minimal error
            case.
        Returns:
            neutral (bool): Whether the nodes independently correspond to a valid
            set of errors.
            flipped_logical_nodes (list): List of qubits nodes for logical
            operators that are flipped by the errors, that were not included
            in the original nodes.
            num_errors (int): Minimum number of errors required to create nodes.
        """
        pass

    @abstractmethod
    def is_cluster_neutral(self, atypical_nodes):
        """
        Determines whether or not the cluster is neutral, meaning that one or more
        errors could have caused the set of atypical nodes (syndrome changes) passed
        to the method.

        Default version here assumes that it is as simple as an an even/odd assessment
        (as for repetition codes, surface codes, etc). This should be overwritten for
        more complex codes. It also should be used with care, by only supplying sets
        of nodes for which the even/odd assessment is valid.
        Args:
            atypical_nodes (dictionary in the form of the return value of string2nodes)
        """
        return not bool(len(atypical_nodes) % 2)
