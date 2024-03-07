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

import itertools
import unittest
from random import choices

from qiskit import QuantumCircuit, execute
from qiskit.providers.fake_provider import FakeJakarta
from qiskit_aer import Aer, AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit_aer.noise.errors import depolarizing_error

from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.circuits.repetition_code import ArcCircuit
from qiskit_qec.circuits.repetition_code import RepetitionCodeCircuit as RepetitionCode
from qiskit_qec.decoders.decoding_graph import DecodingGraph
from qiskit_qec.decoders.hdrg_decoders import BravyiHaahDecoder, UnionFindDecoder
from qiskit_qec.utils import DecodingGraphNode


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


class TestRepCodes(unittest.TestCase):
    """Test the repetition code circuits."""

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
                    DecodingGraphNode(is_logical=True, is_boundary=True, qubits=[0], index=0),
                    DecodingGraphNode(
                        time=0,
                        qubits=[0, 1],
                        index=0,
                    ),
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
            for T in [1, 2, 3]:
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
        code = RepetitionCode(3, 1)
        dec = DecodingGraph(code)
        test_results = {"000 00": 1024, "010 11": 512}
        for method in ["spitz", "naive"]:
            error = (
                "Error: Calculated error probability not correct in"
                + " d=3, T=1 repetition code using calculation method '"
                + method
                + "'."
            )
            p = dec.get_error_probs(test_results, method=method)
            n0 = dec.graph.nodes().index(DecodingGraphNode(time=0, qubits=[0, 1], index=0))
            n1 = dec.graph.nodes().index(DecodingGraphNode(time=0, qubits=[1, 2], index=1))
            # edges in graph aren't directed and could be in any order
            if (n0, n1) in p:
                self.assertTrue(round(p[n0, n1], 2) == 0.33, error)
            else:
                self.assertTrue(round(p[n1, n0], 2) == 0.33, error)


class TestARCCodes(unittest.TestCase):
    """Test the ARC code circuits."""

    def single_error_test(
        self, code
    ):  # NOT run directly by unittest; called by test_graph_constructions
        """
        Insert all possible single qubit errors into the given code,
        and check that each creates a pair of syndrome nodes.
        """

        # determine the neighbourhood of each code qubit
        incident_links = {None: set()}
        link_graph = code._get_link_graph()
        for n, node in enumerate(link_graph.nodes()):
            edges = link_graph.incident_edges(n)
            incident_links[node] = set()
            for edge in edges:
                incident_links[node].add(link_graph.edges()[edge]["link qubit"])
            if node in code.z_logicals:
                incident_links[node].add(None)

        minimal = True
        for basis in [code.basis, code.basis[::-1]]:
            # use the fault enumerator to test possible faults
            qc = code.circuit[basis]
            fe = FaultEnumerator(qc, method="stabilizer")
            blocks = list(fe.generate_blocks())
            fault_paths = list(itertools.chain(*blocks))

            # check that the syndrome for single faults is as it should be
            for _, _, _, output in fault_paths:
                string = "".join([str(c) for c in output[::-1]])
                nodes = code.string2nodes(string)
                # check that it doesn't extend over more than two rounds
                ts = [node.time for node in nodes if not node.is_logical]
                if ts:
                    minimal = minimal and (max(ts) - min(ts)) <= 1
                # check that it corresponds to more than one node (or none)
                self.assertTrue(
                    len(nodes) != 1,
                    "Error: Single error creates only one node",
                )
                # check that it doesn't extend beyond the neigbourhood of a code qubit
                flat_nodes = code.flatten_nodes(nodes)
                link_qubits = set(node.properties["link qubit"] for node in flat_nodes)
                minimal = minimal and link_qubits in incident_links.values()
                self.assertTrue(
                    minimal,
                    "Error: Single error creates too many nodes",
                )
                neutral, flipped_logicals, num = code.check_nodes(nodes)
                # check that the nodes are neutral
                self.assertTrue(
                    neutral and flipped_logicals == [],
                    "Error: Single error nodes are not neutral for string "
                    + string
                    + " which yields "
                    + str(code.check_nodes(nodes))
                    + " for nodes "
                    + str(nodes),
                )
                # and caused by at most a single error
                self.assertTrue(
                    num <= 1,
                    "Error: Nodes seem to be caused by more than one error for "
                    + string
                    + " which yields "
                    + str(code.check_nodes(nodes)),
                )

    def test_graph_construction(self):
        """Test single errors for a range of layouts"""
        triangle = [(0, 1, 2), (2, 3, 4), (4, 5, 0)]
        tadpole = [(0, 1, 2), (2, 3, 4), (4, 5, 0), (2, 7, 6)]
        t_pose = [(0, 1, 2), (2, 3, 4), (2, 5, 6), (6, 7, 8)]
        for links in [triangle, tadpole, t_pose]:
            for resets in [True, False]:
                conditional_resets = [False]
                if resets:
                    conditional_resets.append(True)
                for conditional_reset in conditional_resets:
                    for logical in ["0", "1"]:
                        code = ArcCircuit(
                            links,
                            T=2,
                            barriers=True,
                            delay=1,
                            basis="xy",
                            run_202=False,
                            resets=resets,
                            conditional_reset=conditional_reset,
                            logical=logical,
                        )
                        self.assertTrue(
                            code.resets == resets,
                            "Code has resets="
                            + str(code.resets)
                            + " when it should be "
                            + str(resets)
                            + ".",
                        )
                        self.single_error_test(code)

    def test_202s(self):
        """Test that [[2,0,2]] codes appear when needed and act as required."""
        links = [(0, 1, 2), (2, 3, 4), (4, 5, 6), (6, 7, 0)]
        T = 15
        # first, do they appear when needed
        for run_202 in [True, False]:
            code = ArcCircuit(links, T=T, run_202=run_202, rounds_per_202=5)
            running_202 = False
            for t in range(T):
                tau, _, _ = code._get_202(t)
                if tau:
                    running_202 = True
            self.assertTrue(
                running_202 == run_202,
                "Error: [[2,0,2]] codes not present when required." * run_202
                + "Error: [[2,0,2]] codes present when not required." * (not run_202),
            )
        # second, do they yield non-trivial outputs yet trivial nodes
        code = ArcCircuit(links, T=T, run_202=True, logical="1", rounds_per_202=5)
        backend = Aer.get_backend("aer_simulator")
        counts = backend.run(code.circuit[code.basis]).result().get_counts()
        self.assertTrue(len(counts) > 1, "No randomness in the results for [[2,0,2]] circuits.")
        nodeless = True
        for string in counts:
            nodeless = nodeless and not code.string2nodes(string)
        self.assertTrue(nodeless, "Non-trivial nodes found for noiseless [[2,0,2]] circuits.")

    def test_single_error_202s(self):
        """Test a range of single errors for a code with [[2,0,2]] codes."""
        links = [(0, 1, 2), (2, 3, 4), (4, 5, 0), (2, 7, 6)]
        for T in [21, 25]:
            code = ArcCircuit(links, T, run_202=True, barriers=True, logical="1", rounds_per_202=5)
            assert code.run_202
            # insert errors on a selection of qubits during a selection of rounds
            qc = code.circuit[code.base]
            for q in [0, 2, 5, 6, 1, 7]:
                for t in [2 + 5 * l for l in range(len(links))]:
                    error_qc = QuantumCircuit()
                    for qreg in qc.qregs:
                        error_qc.add_register(qreg)
                    for creg in qc.cregs:
                        error_qc.add_register(creg)
                    barrier_num = 0
                    for gate in qc.data:
                        if gate[0].name == "barrier":
                            barrier_num += 1
                            if barrier_num == 2 * t + 1:
                                if q % 2 == 0:
                                    error_qc.z(code.code_qubit[code.code_index[q]])
                                else:
                                    error_qc.x(code.link_qubit[code.link_index[q]])
                        error_qc.append(gate)
                    counts = AerSimulator().run(error_qc).result().get_counts()
                    for string in counts:
                        # look at only bulk non-conjugate nodes
                        nodes = [
                            node
                            for node in code.string2nodes(string)
                            if "conjugate" not in node.properties and not node.is_logical
                        ]
                        # require at most two (or three for the trivalent vertex or neighbouring aux)
                        self.assertTrue(
                            len(nodes) <= (2 + int(q in (2, 7))),
                            "Too many nodes for a single error in a [[2,0,2]] code for T="
                            + str(T)
                            + ".",
                        )

    def test_feedforward(self):
        """Test that the correct behaviour is seen with feedforward for [[2,0,2]] codes."""
        links = [(0, 1, 2), (2, 3, 4), (4, 5, 6)]
        T = 10
        # try codes with and without feedforward correction
        code = ArcCircuit(
            links, T, barriers=True, basis="xy", color={0: 0, 2: 1, 4: 0, 6: 1}, logical="1"
        )
        correct = True
        # insert an initial bitflip on qubit 2
        test_qcs = []
        for basis in [code.basis, code.basis[::-1]]:
            qc = code.circuit[basis]

            test_qc = QuantumCircuit()
            for qreg in qc.qregs:
                test_qc.add_register(qreg)
            for creg in qc.cregs:
                test_qc.add_register(creg)
            test_qc.x(code.code_qubit[2])
            for gate in qc:
                test_qc.append(gate)
            test_qcs.append(test_qc)
        result = AerSimulator().run(test_qcs).result()
        # check result strings are correct
        for j in range(2):
            counts = result.get_counts(j)
            for string in counts:
                # final result should be same as initial
                correct_final = code.logical + str((int(code.logical) + 1) % 2) + code.logical * 2
                correct = correct and string[0:4] == correct_final
        self.assertTrue(correct, "Result string not as required")

    def test_bases(self):
        """Test that correct rotations are used for basis changes."""
        links = [(0, 1, 2), (2, 3, 4)]
        rightops = True
        for ba in ["x", "y", "z"]:
            for sis in ["x", "y", "z"]:
                basis = ba + sis
                code = ArcCircuit(links, T=3, basis=basis)
                ops = code.circuit[code.base].count_ops()
                if "x" in basis or "y" in basis:
                    rightops = rightops and "h" in ops
                else:
                    rightops = rightops and "h" not in ops
                if "y" in basis:
                    for op in ["s", "sdg"]:
                        rightops = rightops and op in ops
                else:
                    for op in ["s", "sdg"]:
                        rightops = rightops and op not in ops
        self.assertTrue(rightops, "Error: Required rotations for basis changes not present.")

    def test_anisotropy(self):
        """Test that code qubits have neighbors with the opposite color."""
        link_num = 10
        links = [(2 * j, 2 * j + 1, 2 * (j + 1)) for j in range(link_num)]
        code = ArcCircuit(links, T=2)
        color = code.color
        for j in range(1, link_num - 1):
            self.assertTrue(
                color[2 * j] != color[2 * (j - 1)] or color[2 * j] != color[2 * (j + 1)],
                "Error: Code qubit does not have neighbor of oppposite color.",
            )

    def test_transpilation(self):
        """Test correct transpilation to a backend."""
        backend = FakeJakarta()
        links = [(0, 1, 3), (3, 5, 6)]
        schedule = [[(0, 1), (3, 5)], [(3, 1), (6, 5)]]
        code = ArcCircuit(links, schedule=schedule, T=2, delay=1000, logical="0")
        circuit = code.transpile(backend)
        self.assertTrue(code.schedule == schedule, "Error: Given schedule not used.")
        circuit = code.transpile(backend, echo_num=(0, 2))
        self.assertTrue(
            circuit[code.base].count_ops()["x"] == 4, "Error: Wrong echo sequence for link qubits."
        )
        circuit = code.transpile(backend, echo_num=(2, 0))
        self.assertTrue(
            circuit[code.base].count_ops()["x"] == 26, "Error: Wrong echo sequence for code qubits."
        )
        self.assertTrue(
            circuit[code.base].count_ops()["cx"] == 8,
            "Error: Wrong number of cx gates after transpilation.",
        )

    def test_weight(self):
        """Error weighting code test."""
        code = ArcCircuit(
            [
                [0, 1, 2],
                [2, 3, 4],
            ],
            1,
            schedule=[[[0, 1]], [[2, 3]], [[2, 1]], [[4, 3]]],
        )
        dec = DecodingGraph(code)
        test_results = {"000 00": 1024, "010 11": 512}
        for method in ["spitz", "naive"]:
            error = (
                "Error: Calculated error probability not correct for in"
                + " d=3, T=1 repetition code using calculation method '"
                + method
                + "'."
            )
            p = dec.get_error_probs(test_results, method=method)
            node = DecodingGraphNode(time=0, qubits=[0, 2], index=1)
            node.properties["link qubits"] = 1
            n0 = dec.graph.nodes().index(node)
            node = DecodingGraphNode(time=0, qubits=[2, 4], index=0)
            node.properties["link qubits"] = 3
            n1 = dec.graph.nodes().index(node)
            # edges in graph aren't directed and could be in any order
            if (n0, n1) in p:
                self.assertTrue(round(p[n0, n1], 2) == 0.33, error)
                ns = (n0, n1)
            else:
                self.assertTrue(round(p[n1, n0], 2) == 0.33, error)
                ns = (n1, n0)
            pc = code.get_error_coords(test_results, dec)
            self.assertTrue(round(pc[2, 0, 0.2][ns], 2) == 0.33, error)


class TestDecoding(unittest.TestCase):
    """Test decoders for repetition codes"""

    @staticmethod
    def test_empty_decoding_graph():
        """Test initializtion of decoding graphs with None"""
        DecodingGraph(None)

    def clustering_decoder_test(self, Decoder):  # NOT run directly by unittest
        """Test decoding of ARCs and RCCs with clustering decoders"""

        # parameters for test
        d = 8
        p = 0.1
        N = 1000

        # first an RCC
        codes = []
        codes.append(RepetitionCode(d, 1))
        # then a linear ARC
        links = [(2 * j, 2 * j + 1, 2 * (j + 1)) for j in range(d - 1)]
        codes.append(ArcCircuit(links, 0, logical="1"))
        self.assertTrue(codes[-1]._linear, "Linear ARC not recognised as such")
        # then make a bunch of non-linear ARCs
        links_cross = [(2 * j, 2 * j + 1, 2 * (j + 1)) for j in range(d - 2)]
        links_cross.append((2 * (d - 2), 2 * (d - 2) + 1, 2 * (int(d / 2))))
        links_cross.append(((2 * (int(d / 2))), 2 * (d - 1), 2 * (d - 1) + 1))
        codes.append(ArcCircuit(links_cross, 0))
        self.assertTrue(not codes[-1]._linear, "Non-inear ARC not recognised as such")
        # ladder (works for even d)
        half_d = int(d / 2)
        links_ladder = []
        for row in [0, 1]:
            for j in range(half_d - 1):
                delta = row * (2 * half_d - 1)
                links_ladder.append((delta + 2 * j, delta + 2 * j + 1, delta + 2 * (j + 1)))
        q = links_ladder[-1][2] + 1
        for j in range(half_d):
            delta = 2 * half_d - 1
            links_ladder.append((2 * j, q, delta + 2 * j))
            q += 1
        codes.append(ArcCircuit(links_ladder, 0))
        # now run them all and check it works
        for c, code in enumerate(codes):
            decoding_graph = DecodingGraph(code)
            if c >= 0 and Decoder is UnionFindDecoder:
                decoder = Decoder(code, decoding_graph=decoding_graph, use_peeling=False)
            else:
                decoder = Decoder(code, decoding_graph=decoding_graph)
            min_error_num = code.d
            min_error_string = ""
            for _ in range(N):
                # generate random string
                string = "".join([choices(["1", "0"], [1 - p, p])[0] for _ in range(d)])
                for _ in range(code.T):
                    string = string + " " + "0" * (d - 1)
                # get and check corrected_z_logicals
                corrected_z_logicals = decoder.process(string)
                for node in decoder.decoding_graph.logical_nodes:
                    if node.index < len(corrected_z_logicals):
                        error = corrected_z_logicals[node.index] != 1
                        if error:
                            error_num = string.split(" ", maxsplit=1)[0].count("0")
                            if error_num < min_error_num:
                                min_error_num = error_num
                                min_error_string = string
            # check that min num errors to cause logical errors >d/3
            self.assertTrue(
                min_error_num > d / 3,
                str(min_error_num)
                + " errors cause logical error despite d="
                + str(code.d)
                + " for code "
                + str(c)
                + " with "
                + min_error_string
                + "."
                + " Corresponding clusters are "
                + str(decoder.cluster(code.string2nodes(string, all_logicals=True)))
                + ".",
            )

    def heavy_hex_test(self, Decoder):  # NOT run directly by unittest
        """Test decoding of heavy hex ARC"""
        links = [
            (0, 1, 2),
            (2, 3, 4),
            (4, 5, 6),
            (6, 7, 8),
            (8, 9, 10),
            (10, 11, 12),
            (0, 14, 18),
            (4, 15, 22),
            (8, 16, 26),
            (12, 17, 30),
            (18, 19, 20),
            (20, 21, 22),
            (22, 23, 24),
            (24, 25, 26),
            (26, 27, 28),
            (28, 29, 30),
            (30, 31, 32),
            (20, 33, 39),
            (24, 34, 43),
            (28, 35, 47),
            (32, 36, 51),
            (37, 38, 39),
            (39, 40, 41),
            (41, 42, 43),
            (43, 44, 45),
            (45, 46, 47),
            (47, 48, 49),
            (49, 50, 51),
            (37, 52, 56),
            (41, 53, 60),
            (45, 54, 64),
            (49, 55, 68),
            (56, 57, 58),
            (58, 59, 60),
            (60, 61, 62),
            (62, 63, 64),
            (64, 65, 66),
            (66, 67, 68),
            (68, 69, 70),
            (58, 71, 77),
            (62, 72, 81),
            (66, 73, 85),
            (70, 74, 89),
            (75, 76, 77),
            (77, 78, 79),
            (79, 80, 81),
            (81, 82, 83),
            (83, 84, 85),
            (85, 86, 87),
            (87, 88, 89),
            (75, 90, 94),
            (79, 91, 98),
            (83, 92, 102),
            (87, 93, 106),
            (94, 95, 96),
            (96, 97, 98),
            (98, 99, 100),
            (100, 101, 102),
            (102, 103, 104),
            (104, 105, 106),
            (106, 107, 108),
            (96, 109, 114),
            (100, 110, 118),
            (104, 111, 122),
            (108, 112, 126),
            (114, 115, 116),
            (116, 117, 118),
            (118, 119, 120),
            (120, 121, 122),
            (122, 123, 124),
            (124, 125, 126),
        ]
        schedule = [
            [
                (0, 14),
                (2, 3),
                (4, 15),
                (6, 7),
                (8, 16),
                (10, 11),
                (12, 17),
                (18, 19),
                (22, 23),
                (26, 27),
                (30, 31),
                (20, 33),
                (24, 34),
                (28, 35),
                (32, 36),
                (39, 40),
                (43, 44),
                (47, 48),
                (37, 52),
                (41, 53),
                (45, 54),
                (49, 55),
                (56, 57),
                (60, 61),
                (64, 65),
                (68, 69),
                (58, 71),
                (62, 72),
                (66, 73),
                (70, 74),
                (77, 78),
                (81, 82),
                (85, 86),
                (75, 90),
                (79, 91),
                (83, 92),
                (87, 93),
                (94, 95),
                (98, 99),
                (102, 103),
                (106, 107),
                (96, 109),
                (100, 110),
                (104, 111),
                (108, 112),
                (114, 115),
                (118, 119),
                (122, 123),
            ],
            [
                (0, 1),
                (4, 5),
                (8, 9),
                (18, 14),
                (22, 15),
                (26, 16),
                (30, 17),
                (20, 21),
                (24, 25),
                (28, 29),
                (39, 33),
                (43, 34),
                (47, 35),
                (51, 36),
                (37, 38),
                (41, 42),
                (45, 46),
                (49, 50),
                (56, 52),
                (60, 53),
                (64, 54),
                (68, 55),
                (58, 59),
                (62, 63),
                (66, 67),
                (77, 71),
                (81, 72),
                (85, 73),
                (89, 74),
                (75, 76),
                (79, 80),
                (83, 84),
                (87, 88),
                (94, 90),
                (98, 91),
                (102, 92),
                (106, 93),
                (96, 97),
                (100, 101),
                (104, 105),
                (114, 109),
                (118, 110),
                (122, 111),
                (126, 112),
                (116, 117),
                (120, 121),
                (124, 125),
            ],
            [
                (2, 1),
                (4, 3),
                (6, 5),
                (8, 7),
                (10, 9),
                (12, 11),
                (22, 21),
                (26, 25),
                (30, 29),
                (20, 19),
                (24, 23),
                (28, 27),
                (32, 31),
                (39, 38),
                (43, 42),
                (47, 46),
                (51, 50),
                (41, 40),
                (45, 44),
                (49, 48),
                (60, 59),
                (64, 63),
                (68, 67),
                (58, 57),
                (62, 61),
                (66, 65),
                (70, 69),
                (77, 76),
                (81, 80),
                (85, 84),
                (89, 88),
                (79, 78),
                (83, 82),
                (87, 86),
                (98, 97),
                (102, 101),
                (106, 105),
                (96, 95),
                (100, 99),
                (104, 103),
                (108, 107),
                (118, 117),
                (122, 121),
                (126, 125),
                (116, 115),
                (120, 119),
                (124, 123),
            ],
        ]
        code = ArcCircuit(
            links, 10, schedule=schedule, run_202=False, basis="zx", logical="0", resets=True
        )
        if Decoder is UnionFindDecoder:
            decoder = Decoder(code, use_peeling=False)
        else:
            decoder = Decoder(code)
        string = (
            "110100001100010110010011110110100011111100000101100101 "
            + "01111010000111111110111110111111000010001000111111101110111111101010001 "
            + "11110100001000110001011110110111100000111111011011011100011001000110111 "
            + "11110100010000010110010100110110000011010110101000010101011100001000111 "
            + "11110101010001001010001110001111000011011100111001011100001010001000111 "
            + "01010000110001100011010110001110010011000000111010011000000100011010011 "
            + "11001000110001110011010011101101010101000110101010010000000000011111011 "
            + "11000000100011111001010101101011010101000100111110010000001010001101011 "
            + "11000000001000100000010000001011001101000100100111010000110001101111100 "
            + "11000001000000001001010001111100010000100111011110000000011000010000101 "
            + "01000010000000000000110000001100011000000100000010000010000100010000000"
        )
        self.assertTrue(
            decoder.process(string)[0] == 0,
            "Incorrect decoding for example string with heavy-hex ARC.",
        )

    def test_bravyi_haah(self):
        """Test decoding of ARCs and RCCs with Bravyi Haah"""
        self.clustering_decoder_test(BravyiHaahDecoder)
        self.heavy_hex_test(BravyiHaahDecoder)

    def test_union_find(self):
        """Test decoding of ARCs and RCCs with Union Find"""
        self.clustering_decoder_test(UnionFindDecoder)
        self.heavy_hex_test(UnionFindDecoder)


if __name__ == "__main__":
    unittest.main()
