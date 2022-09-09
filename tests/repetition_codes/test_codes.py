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

"""Run codes and decoders."""

import sys
import unittest

from qiskit import Aer, QuantumCircuit, execute
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors import depolarizing_error
from qiskit_qec.circuits.repetition_code import RepetitionCodeCircuit as RepetitionCode
from qiskit_qec.decoders.decoding_graph import DecodingGraph

sys.path.append("../../")


def get_syndrome(code, noise_model, shots=1024):
    """Runs a code to get required results."""
    circuits = [code.circuit[log] for log in ["0", "1"]]

    job = execute(
        circuits,
        Aer.get_backend("aer_simulator_stabilizer"),
        noise_model=noise_model,
        shots=shots,
    )
    results = {}
    for log in ["0", "1"]:
        results[log] = job.result().get_counts(log)

    return results


def get_noise(p_meas, p_gate):
    """Define a noise model."""
    error_gate1 = depolarizing_error(p_gate, 1)
    error_gate2 = error_gate1.tensor(error_gate1)

    noise_model = NoiseModel()
    noise_model.add_all_qubit_readout_error([[1 - p_meas, p_meas], [p_meas, 1 - p_meas]])
    noise_model.add_all_qubit_quantum_error(error_gate1, ["h"])
    noise_model.add_all_qubit_quantum_error(error_gate2, ["cx"])

    return noise_model


class TestCodes(unittest.TestCase):
    """Test the topological codes module."""

    def single_error_test(
        self, code
    ):  # NOT run directly by unittest; called by test_graph_constructions
        """
        Insert all possible single qubit errors into the given code,
        and check that each creates a pair of syndrome nodes.
        """

        for logical in ["0", "1"]:
            qc = code.circuit[logical]
            blank_qc = QuantumCircuit()
            for qreg in qc.qregs:
                blank_qc.add_register(qreg)
            for creg in qc.cregs:
                blank_qc.add_register(creg)
            error_circuit = {}
            circuit_name = {}
            depth = len(qc)
            for j in range(depth):
                qubits = qc.data[j][1]
                for qubit in qubits:
                    for error in ["x", "y", "z"]:
                        temp_qc = blank_qc.copy()
                        temp_qc.name = str((j, qubit, error))
                        temp_qc.data = qc.data[0:j]
                        getattr(temp_qc, error)(qubit)
                        temp_qc.data += qc.data[j : depth + 1]
                        circuit_name[(j, qubit, error)] = temp_qc.name
                        error_circuit[temp_qc.name] = temp_qc

            simulator = Aer.get_backend("aer_simulator")
            job = simulator.run(list(error_circuit.values()), shots=1)

            for j in range(depth):
                qubits = qc.data[j][1]
                for qubit in qubits:
                    for error in ["x", "y", "z"]:
                        results = job.result().get_counts(str((j, qubit, error)))
                        for string in results:
                            nodes = code.string2nodes(string, logical=logical)
                            self.assertIn(
                                len(nodes),
                                [0, 2],
                                "Error of type "
                                + error
                                + " on qubit "
                                + str(qubit)
                                + " at depth "
                                + str(j)
                                + " creates "
                                + str(len(nodes))
                                + " nodes in syndrome graph, instead of 2.",
                            )

    def test_string2nodes_1(self):
        """Test string2nodes with different logical values."""
        code = RepetitionCode(3, 2)
        s0 = "0 0  01 00 01"
        s1 = "1 1  01 00 01"
        self.assertTrue(
            code.string2nodes(s0, logical="0") == code.string2nodes(s1, logical="1"),
            "Error: Incorrect nodes from results string",
        )

    def test_string2nodes_2(self):
        """Assert string2nodes produces correct output"""
        test_info = [
            [
                (5, 0),
                "00001",
                [
                    {"time": 0, "qubits": [0], "is_boundary": True, "element": 0},
                    {"time": 0, "qubits": [0, 1], "is_boundary": False, "element": 0},
                ],
            ]
        ]
        for cur_test in test_info:
            d, T = cur_test[0]
            outstring = cur_test[1]
            results = cur_test[2]
            code = RepetitionCode(d, T, resets=True)
            code.readout()
            self.assertTrue(
                results == code.string2nodes(outstring),
                f"string2nodes on {outstring}"
                + f"produced bad output: { code.string2nodes(outstring)} instead of {results}",
            )

    def test_graph_construction(self):
        """Check that single errors create a pair of nodes for all types of code."""
        codes = {}
        for d in [2, 3]:
            for T in [1, 2]:
                for xbasis in [False, True]:
                    for resets in [False, True]:
                        for delay in [0, 16]:
                            codes[d, T, xbasis, resets, delay] = RepetitionCode(
                                d, T, xbasis=xbasis, resets=resets, delay=delay
                            )
        for params, code in codes.items():
            d, T = params[0:2]
            self.single_error_test(code)
            if len(params) == 5:
                d, T, xbasis, resets, delay = params
                if delay > 0 and T > 1:
                    num_delays = code.circuit["0"].count_ops()["delay"]
                    self.assertTrue(
                        num_delays == (d - 1) * (T - 1),
                        "Error: wrong number of delay gates.",
                    )

    def test_weight(self):
        """Error weighting code test."""
        error = (
            "Error: Calculated error probability not correct for "
            + "test result '0 0  11 00' in d=3, T=1 repetition code."
        )
        code = RepetitionCode(3, 1)
        dec = DecodingGraph(code)
        test_results = {"000 00": 1024, "010 11": 512}
        p = dec.get_error_probs(test_results)
        n0 = dec.graph.nodes().index(
            {"time": 0, "is_boundary": False, "qubits": [0, 1], "element": 0}
        )
        n1 = dec.graph.nodes().index(
            {"time": 0, "is_boundary": False, "qubits": [1, 2], "element": 1}
        )
        # edges in graph aren't directed and could be in any order
        if (n0, n1) in p:
            self.assertTrue(round(p[n0, n1], 2) == 0.33, error)
        else:
            self.assertTrue(round(p[n1, n0], 2) == 0.33, error)


if __name__ == "__main__":
    unittest.main()
