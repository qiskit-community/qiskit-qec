"""Abstract object to construct CSS matching decoders for circuit noise."""

from copy import deepcopy, copy
from math import log
from abc import ABC, abstractmethod
from typing import List, Tuple, Dict, Set
import logging

import json

from sympy import Poly, symbols, Symbol
from qiskit import QuantumCircuit
import retworkx as rx
import networkx as nx
from pymatching import Matching

from qiskit_qec.exceptions import QiskitQECError
from qiskit_qec.utils.indexer import Indexer
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel
from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.decoders.decoding_graph import CSSDecodingGraph
from qiskit_qec.decoders.temp_code_util import temp_syndrome, temp_gauge_products


class CircuitModelMatchingDecoder(ABC):
    """Matching decoder for circuit noise."""

    METHOD_RETWORKX: str = "retworkx"
    METHOD_PYMATCHING: str = "pymatching"
    AVAILABLE_METHODS = {METHOD_RETWORKX, METHOD_PYMATCHING}

    def __init__(
        self,
        css_x_gauge_ops: List[Tuple[int]],
        css_x_stabilizer_ops: List[Tuple[int]],
        css_x_boundary: List[int],
        css_z_gauge_ops: List[Tuple[int]],
        css_z_stabilizer_ops: List[Tuple[int]],
        css_z_boundary: List[int],
        circuit: QuantumCircuit,
        model: PauliNoiseModel,
        basis: str,
        round_schedule: str,
        blocks: int,
        method: str,
        uniform: bool,
        graph: rx.PyGraph = None,
    ):
        """Create a matching decoder.

        Specialized to (subsystem) CSS codes encoding one logical qubit,
        and quantum circuits that prepare and measure in the Z or X basis
        and do repeated Pauli measurements.

        css_x_gauge_ops : list of supports of X gauge operators
        css_x_stabilizer_ops : list of supports of X stabilizers
        css_x_boundary : list of qubits along the X-type boundary
        css_z_gauge_ops : list of supports of Z gauge operators
        css_z_stabilizer_ops : list of supports of Z stabilizers
        css_x_boundary : list of qubits along the Z-type boundary
        circuit : entire quantum circuit to build the decoding graph
        model : noise for operations in the circuit
        basis : initializaton and measurement basis ("x" or "z")
        round_schedule : gauge measurements in each block
        blocks : number of measurement blocks
        method : matching implementation
        uniform : use same edge weight everywhere?

        """
        self.css_x_gauge_ops = css_x_gauge_ops
        self.css_x_stabilizer_ops = css_x_stabilizer_ops
        self.css_x_boundary = css_x_boundary
        self.css_z_gauge_ops = css_z_gauge_ops
        self.css_z_stabilizer_ops = css_z_stabilizer_ops
        self.css_z_boundary = css_z_boundary
        self.model = model
        self.blocks = blocks

        if self.blocks < 1:
            raise QiskitQECError("expected positive integer for blocks")
        self.round_schedule = round_schedule
        if set(self.round_schedule) > set("xyz"):
            raise QiskitQECError("expected round schedule of 'x', 'y', 'z' chars")
        self.basis = basis
        if not self.basis in ("x", "z"):
            raise QiskitQECError("expected basis to be 'x' or 'z'")

        self.uniform = uniform

        if method not in self.AVAILABLE_METHODS:
            raise QiskitQECError("fmethod {method} is not supported.")
        self.method = method
        if self.method == "METHOD_PYMATCHING":
            self.usepymatching = True
        else:
            self.usepymatching = False

        self.z_gauge_products = temp_gauge_products(self.css_z_stabilizer_ops, self.css_z_gauge_ops)
        self.x_gauge_products = temp_gauge_products(self.css_x_stabilizer_ops, self.css_x_gauge_ops)

        if graph:
            (
                self.idxmap,
                self.node_layers,
                self.graph,
                self.pymatching_indexer,
                self.layer_types,
            ) = self._process_graph(graph, blocks, round_schedule, basis)
        else:
            dg = CSSDecodingGraph(
                css_x_gauge_ops,
                css_x_stabilizer_ops,
                css_x_boundary,
                css_z_gauge_ops,
                css_z_stabilizer_ops,
                css_z_boundary,
                blocks,
                round_schedule,
                basis,
            )

            (
                self.idxmap,
                self.node_layers,
                self.graph,
                self.pymatching_indexer,
                self.layer_types,
            ) = (dg.idxmap, dg.node_layers, dg.graph, dg.pymatching_indexer, dg.layer_types)

        logging.debug("layer_types = %s", self.layer_types)

        self.ridxmap = {v: k for k, v in self.idxmap.items()}

        self.circuit = circuit
        self.event_map = {}
        self.parameters = model.get_operations()
        self.edge_weight_polynomials = {}
        self.symbols = None
        if not self.uniform:
            fe = FaultEnumerator(circuit, order=1, method="stabilizer", model=self.model)
            self.event_map = self._enumerate_events(
                self.css_x_gauge_ops,
                self.css_x_stabilizer_ops,
                self.css_x_boundary,
                self.x_gauge_products,
                self.css_z_gauge_ops,
                self.css_z_stabilizer_ops,
                self.css_z_boundary,
                self.z_gauge_products,
                self.blocks,
                self.round_schedule,
                self.basis,
                self.layer_types,
                fe,
            )
            logging.debug("event_map = %s", self.event_map)
            self.symbols, self.edge_weight_polynomials = self._edge_weight_polynomials(
                self.model, self.event_map
            )
            logging.debug("symbols = %s", self.symbols)
            logging.debug("edge_weight_polynomials = %s", self.edge_weight_polynomials)
            self.graph = self._revise_decoding_graph(
                self.idxmap, self.graph, self.edge_weight_polynomials
            )

        self.length = {}  # recomputed in update_edge_weights
        self.path = {}  # recomputed in update_edge_weights
        self.pymatching = None  # constructed in update_edge_weights

    def _process_graph(self, graph, blocks, round_schedule, basis):

        pymatching_indexer = Indexer()

        # symmetrize hook errors
        for j, edge in enumerate(graph.edges()):
            n0, n1 = graph.edge_list()[j]
            source = graph.nodes()[n0]
            target = graph.nodes()[n1]
            if source["time"] != target["time"]:
                if source["is_boundary"] == target["is_boundary"] == False:
                    new_source = source.copy()
                    new_source["time"] = target["time"]
                    nn0 = graph.nodes().index(new_source)
                    new_target = target.copy()
                    new_target["time"] = source["time"]
                    nn1 = graph.nodes().index(new_target)
                    graph.add_edge(nn0, nn1, edge)

        edges_to_remove = []
        for j, edge in enumerate(graph.edges()):

            n0, n1 = graph.edge_list()[j]
            source = graph.nodes()[n0]
            target = graph.nodes()[n1]

            # add the required attributes
            # highlighted', 'measurement_error','qubit_id' and 'error_probability'
            edge["highlighted"] = False
            edge["measurement_error"] = int(source["time"] != target["time"])
            if edge["qubits"]:
                edge["qubit_id"] = pymatching_indexer[edge["qubits"][0]]
            else:
                edge["qubit_id"] = -1
            edge["error_probability"] = 0.01 * edge["weight"]

            # make it so times of boundary/boundary nodes agree
            if source["is_boundary"] and not target["is_boundary"]:
                if source["time"] != target["time"]:
                    new_source = source.copy()
                    new_source["time"] = target["time"]
                    n = graph.add_node(new_source)
                    edge["measurement_error"] = 0
                    edges_to_remove.append((n0, n1))
                    graph.add_edge(n, n1, edge)

        # remove old boundary/boundary nodes
        for n0, n1 in edges_to_remove:
            graph.remove_edge(n0, n1)

        for n0, source in enumerate(graph.nodes()):
            for n1, target in enumerate(graph.nodes()):

                # add weightless nodes connecting different boundaries
                if source["time"] == target["time"]:
                    if source["is_boundary"] and target["is_boundary"]:
                        if source["qubits"] != target["qubits"]:
                            edge = {
                                "measurement_error": 0,
                                "weight": 0,
                                "highlighted": False,
                                "error_probability": 0.0,
                            }
                            edge["qubits"] = list(
                                set(source["qubits"]).intersection((set(target["qubits"])))
                            )
                            if edge["qubits"]:
                                edge["qubit_id"] = pymatching_indexer[edge["qubit_id"][0]]
                            else:
                                edge["qubit_id"] = -1
                            if (n0, n1) not in graph.edge_list():
                                graph.add_edge(n0, n1, edge)

                # connect one of the boundaries at different times
                if target["time"] == source["time"] + 1:
                    if source["qubits"] == target["qubits"] == [0]:
                        edge = {
                            "qubits": [],
                            "qubit_id": -1,
                            "measurement_error": 0,
                            "weight": 0,
                            "highlighted": False,
                            "error_probability": 0.0,
                        }
                        if (n0, n1) not in graph.edge_list():
                            graph.add_edge(n0, n1, edge)

        # symmetrize edges
        for j, edge in enumerate(graph.edges()):
            n0, n1 = graph.edge_list()[j]
            if (n1, n0) not in graph.edge_list():
                graph.add_edge(n1, n0, edge)

        idxmap = {}
        for n, node in enumerate(graph.nodes()):
            idxmap[node["time"], tuple(node["qubits"])] = n

        node_layers = []
        for node in graph.nodes():
            time = node["time"]
            if len(node_layers) < time + 1:
                node_layers += [[]] * (time + 1 - len(node_layers))
            node_layers[time].append(node["qubits"])

        # create a list of decoding graph layer types
        # the entries are 'g' for gauge and 's' for stabilizer
        layer_types = []
        last_step = basis
        for _ in range(blocks):
            for step in round_schedule:
                if basis == "z" and step == "z" and last_step == "z":
                    layer_types.append("g")
                elif basis == "z" and step == "z" and last_step == "x":
                    layer_types.append("s")
                elif basis == "x" and step == "x" and last_step == "x":
                    layer_types.append("g")
                elif basis == "x" and step == "x" and last_step == "z":
                    layer_types.append("s")
                last_step = step
        if last_step == basis:
            layer_types.append("g")
        else:
            layer_types.append("s")

        return idxmap, node_layers, graph, pymatching_indexer, layer_types

    def _revise_decoding_graph(
        self,
        idxmap: Dict[Tuple[int, List[int]], int],
        graph: rx.PyGraph,
        edge_weight_polynomials: Dict[Tuple[int, Tuple[int]], Dict[Tuple[int, Tuple[int]], Poly]],
    ) -> rx.PyGraph:
        """Add edge weight polynomials to the decoding graph g and prune it.

        Update attribute "weight_poly" on decoding graph edges contained in
        edge_weight_polynomials. Remove all other edges from the decoding graph
        that have non-zero weight.
        """
        for s1, sub in edge_weight_polynomials.items():
            for s2, wpoly in sub.items():
                if s1 not in idxmap:
                    raise QiskitQECError(f"vertex {s1} not in decoding graph")
                if s2 not in idxmap:
                    raise QiskitQECError(f"vertex {s2} not in decoding graph")
                if not graph.has_edge(idxmap[s1], idxmap[s2]):
                    raise QiskitQECError("edge {s1} - {s2} not in decoding graph")
                    # TODO: new edges may be needed for hooks, but raise exception if we're surprised

                data = graph.get_edge_data(idxmap[s1], idxmap[s2])
                data["weight_poly"] = wpoly
        remove_list = []
        for source, target in graph.edge_list():
            edge_data = graph.get_edge_data(source, target)
            if "weight_poly" not in edge_data and edge_data["weight"] != 0:
                if (
                    "is_boundary" in graph.nodes()[source] and graph.nodes()[source]["is_boundary"]
                ) or (
                    "is_boundary" in graph.nodes()[target] and graph.nodes()[target]["is_boundary"]
                ):
                    # pymatching requires all qubit_ids from 0 to
                    # n-1 to appear as edge properties. If we
                    # remove too many edges, not all qubit_ids will
                    # be present. We solve this by keeping unweighted
                    # that are edges connected to the boundary. Instead
                    # we increase their weight to a large value.
                    edge_data["weight"] = 1000000
                    logging.info("increase edge weight (%d, %d)", source, target)
                else:
                    # Otherwise, remove the edge
                    remove_list.append((source, target))
                    logging.info("remove edge (%d, %d)", source, target)
        graph.remove_edges_from(remove_list)
        return graph

    def update_edge_weights(self, model: PauliNoiseModel):
        """Evaluate the numerical edge weights and update graph data.

        For each edge in the decoding graph that has a "weight_poly"
        property, evaluate the polynomial at the given model parameters
        and set the corresponding "weight" property. Once this is done,
        recompute the shortest paths between pairs of vertices
        in the decoding graph.

        model is a PauliNoiseModel whose error probabilities have been
        previously assigned. The probabilities are then assigned to
        the variables in self.symbols.

        Updates properties of self.graph.
        If not using pymatching, updates sets self.length and self.path.
        If using pymatching, constructs a pymatching object self.pymatching.

        Note that if self.uniform is True, it is still necessary to call
        this function to construct the matching decoder object (pymatching)
        or compute the shortest paths between vertices in the decoding
        graph (retworkx).

        Args:
            model: moise model
        """
        parameter_values = [model.get_error_probability(name) for name in self.parameters]
        if not self.uniform:
            if len(parameter_values) != len(self.parameters):
                raise QiskitQECError("wrong number of error rate parameters")
            symbol_list = [self.symbols[s] for s in self.parameters]
            assignment = dict(zip(symbol_list, parameter_values))
            logging.info("update_edge_weights %s", str(assignment))
            # P(chain) = \prod_i (1-p_i)^{1-l(i)}*p_i^{l(i)}
            #          \propto \prod_i ((1-p_i)/p_i)^{l(i)}
            # -log P(chain) \propto \sum_i -log[((1-p_i)/p_i)^{l(i)}]
            # p_i is the probability that edge i carries an error
            # l(i) is 1 if the link belongs to the chain and 0 otherwise
            for source, target in self.graph.edge_list():
                edge_data = self.graph.get_edge_data(source, target)
                if "weight_poly" in edge_data:
                    logging.info(
                        "update_edge_weights (%d, %d) %s",
                        source,
                        target,
                        str(edge_data["weight_poly"]),
                    )
                    restriction = {x: assignment[x] for x in edge_data["weight_poly"].gens}
                    p = edge_data["weight_poly"].eval(restriction).evalf()
                    assert p < 0.5, "edge flip probability too large"
                    edge_data["weight"] = log((1 - p) / p)
                    edge_data["error_probability"] = p
        if not self.usepymatching:
            # Recompute the shortest paths between pairs of vertices
            # in the decoding graph
            edge_cost_fn = lambda edge: edge["weight"]
            length = rx.all_pairs_dijkstra_path_lengths(self.graph, edge_cost_fn)
            self.length = {s: dict(length[s]) for s in length}
            path = rx.all_pairs_dijkstra_shortest_paths(self.graph, edge_cost_fn)
            self.path = {s: {t: list(path[s][t]) for t in path[s]} for s in path}
        else:
            nxgraph = ret2net(self.graph)
            self.pymatching = Matching(nxgraph)

    def _enumerate_events(
        self,
        css_x_gauge_ops: List[Tuple[int]],
        css_x_stabilizer_ops: List[Tuple[int]],
        css_x_boundary: List[int],
        x_gauge_products: List[int],
        css_z_gauge_ops: List[Tuple[int]],
        css_z_stabilizer_ops: List[Tuple[int]],
        css_z_boundary: List[int],
        z_gauge_products: List[int],
        blocks: int,
        round_schedule: str,
        basis: str,
        layer_types: List[str],
        fault_enumerator: FaultEnumerator,
    ) -> Dict[Tuple[int, Tuple[int]], Dict[Tuple[int, Tuple[int]], Dict[List[str], int]]]:
        """Enumerate fault events in the input circuit.

        Use the code definition to identify highlighted edges
        in a decoding graph and return a dict containing the events
        that highlight each edge.

        The basis input value 'x' or 'z' informs whether to
        look at the Z error or X error syndrome, respectively.

        fault_enumerator is a FaultEnumerator object.

        Return a dict containing the total number of events of each
        type: event_map[v0][v1][name][pauli] contains the number of
        events where a gate "name" fails with error "pauli" and
        the edge between v0 and v1 is highlighted.

        Args:
            css_x_gauge_ops: x gauge ops
            css_x_stabilizer_ops: x stabilizer ops
            css_x_boundary: x boundary
            x_gauge_products: x gauge products
            css_z_gauge_ops: z gauge ops
            css_z_stabilizer_ops: z stabilizer ops
            css_z_boundary: z boundary
            z_gauge_products: z gauge products
            blocks: blocks
            round_schedule:
            basis: basis
            layer_types: layer types
            fault_enumerator: fault enumerator

        Returns:
            Events map
        """
        event_map = {}
        for event in fault_enumerator.generate():
            # Unpack the event data
            # Select the first element since order = 1
            ctr = event[0]  # event counter
            comb = event[1][0]  # combination of faulty operations
            pauli = event[2][0]  # Pauli error string
            outcome = event[3]  # result of simulation
            logging.debug("event %d %s %s %s", ctr, comb, pauli, outcome)

            (
                x_gauge_outcomes,
                z_gauge_outcomes,
                final_outcomes,
            ) = self._partition_outcomes(blocks, round_schedule, outcome)

            # Compute the highlighted vertices
            # Note that this only depends on the stabilizers at each
            # time and does not require an explicit decoding graph
            gauge_outcomes, highlighted = self._highlighted_vertices(
                css_x_gauge_ops,
                css_x_stabilizer_ops,
                css_x_boundary,
                x_gauge_products,
                css_z_gauge_ops,
                css_z_stabilizer_ops,
                css_z_boundary,
                z_gauge_products,
                basis,
                layer_types,
                x_gauge_outcomes,
                z_gauge_outcomes,
                final_outcomes,
            )
            logging.debug("gauge_outcomes %s", gauge_outcomes)
            logging.debug("highlighted %s", highlighted)
            # Examine the highlighted vertices to find the edge of the
            # decoding graph that corresponds with this fault event
            if len(highlighted) > 2:
                raise QiskitQECError("too many highlighted vertices for a " + "single fault event")
            if len(highlighted) == 1:  # _highlighted_vertices highlights the boundary
                raise QiskitQECError("only one highlighted vertex for a " + "single fault event")
            if len(highlighted) == 2:
                v0 = highlighted[0]
                v1 = highlighted[1]
                if basis == "z":
                    boundary = css_z_boundary
                elif basis == "x":
                    boundary = css_x_boundary
                # Is the special boundary vertex highlighted?
                if v1 == (0, tuple(boundary[0])):
                    # Replace it with an adjacent vertex
                    for b in boundary:
                        assert len(b) == 1  # Assume each b has one element
                        isect = set(b).intersection(set(v0[1]))
                        if len(isect) > 0:
                            v1 = (v0[0], tuple(b))
                            break
                submap1 = event_map.setdefault(v0, {})
                submap2 = submap1.setdefault(v1, {})
                submap3 = submap2.setdefault(comb, {})
                eventcount = submap3.setdefault(pauli, 0)
                submap3[pauli] = eventcount + 1
        return event_map

    @abstractmethod
    def _partition_outcomes(
        self, blocks: int, round_schedule: str, outcome: List[int]
    ) -> Tuple[List[List[int]], List[List[int]], List[int]]:
        """Process the raw outcome and return results.

        blocks = number of repetition of round_schedule
        round_schedule = string of z and x characters
        outcome = list of 0, 1 outcomes

        Return lists x_gauge_outcomes, z_gauge_outcomes, final_outcomes.
        """
        raise NotImplementedError("Not implemented.")

    def _edge_weight_polynomials(
        self,
        model: PauliNoiseModel,
        event_map: Dict[Tuple[int, Tuple[int]], Dict[Tuple[int, Tuple[int]], Dict[List[str], int]]],
    ) -> Tuple[
        Dict[str, Symbol],
        Dict[Tuple[int, Tuple[int]], Dict[Tuple[int, Tuple[int]], Poly]],
    ]:
        """Compute edge weight polynomials given the error events.

        event_map is the output of _enumerate_events
        """
        symbs = {n: symbols(n) for n in model.get_operations()}
        edge_weight_expressions = {}
        for n1, submap1 in event_map.items():
            for n2 in submap1.keys():
                # check the names in the event map and warn if symbs
                # does not contain one of the names
                for name in event_map[n1][n2].keys():
                    if name not in symbs:
                        logging.warning("%s in event_map but not in model", name)
                # construct a linear approximation to the edge probability
                # using the weights from the noise model
                expr = 0
                for name in self.model.get_operations():
                    if name in event_map[n1][n2]:
                        for pauli, count in event_map[n1][n2][name].items():
                            expr += count * model.get_pauli_weight(name, pauli) * symbs[name]
                map1 = edge_weight_expressions.setdefault(n1, {})
                map1[n2] = Poly(expr)
        return symbs, edge_weight_expressions

    def process(self, outcomes, export=False, filename="decoding.json"):
        """Process a set of outcomes and return corrected final outcomes.

        Be sure to have called update_edge_weights for the
        noise parameters so that the edge weights are updated.

        The result is a list of code.n integers that are 0 or 1.
        These are the corrected values of the final transversal
        measurement in the basis given by self.basis.
        """
        logging.debug("process: outcomes = %s", outcomes)

        x_gauge_outcomes, z_gauge_outcomes, final_outcomes = self._partition_outcomes(
            self.blocks, self.round_schedule, outcomes
        )

        gauge_outcomes, highlighted = self._highlighted_vertices(
            self.css_x_gauge_ops,
            self.css_x_stabilizer_ops,
            self.css_x_boundary,
            self.x_gauge_products,
            self.css_z_gauge_ops,
            self.css_z_stabilizer_ops,
            self.css_z_boundary,
            self.z_gauge_products,
            self.basis,
            self.layer_types,
            x_gauge_outcomes,
            z_gauge_outcomes,
            final_outcomes,
        )
        logging.info("process: gauge_outcomes = %s", gauge_outcomes)
        logging.info("process: final_outcomes = %s", final_outcomes)
        logging.info("process: highlighted = %s", highlighted)

        if not self.usepymatching:
            matching = self._compute_matching(self.idxmap, self.length, highlighted)
            logging.info("process: matching = %s", matching)
            qubit_errors, measurement_errors = self._compute_error_correction(
                self.graph,
                self.idxmap,
                self.path,
                matching,
                highlighted,
                export,
                filename,
            )
            logging.info("process: qubit_errors = %s", qubit_errors)
            logging.debug("process: measurement_errors = %s", measurement_errors)
        else:
            # Input: highlighted (list of highlighted vertices)
            # Output: qubit_errors (list of qubits to flip)
            syndrome = [0] * len(self.idxmap)
            for vertex in highlighted:
                syndrome[self.idxmap[vertex]] = 1
            try:
                correction = self.pymatching.decode(syndrome)
            except AttributeError as attrib_error:
                raise QiskitQECError("Did you call update_edge_weights?") from attrib_error
            qubit_errors = []
            for i, corr in enumerate(correction):
                if corr == 1:
                    qubit_errors.append(self.pymatching_indexer.rlookup(i))
            logging.info("process: qubit_errors = %s", qubit_errors)

        corrected_outcomes = copy(final_outcomes)
        for i in qubit_errors:
            if i != -1:
                corrected_outcomes[i] = (corrected_outcomes[i] + 1) % 2
        logging.info("process: corrected_outcomes = %s", corrected_outcomes)
        if self.basis == "z":
            test = temp_syndrome(corrected_outcomes, self.css_z_stabilizer_ops)
        elif self.basis == "x":
            test = temp_syndrome(corrected_outcomes, self.css_x_stabilizer_ops)
        logging.debug("process: test syndrome = %s", test)
        if sum(test) != 0:
            raise QiskitQECError("decoder failure: syndrome should be trivial!")
        return corrected_outcomes

    def _highlighted_vertices(
        self,
        css_x_gauge_ops: List[Tuple[int]],
        css_x_stabilizer_ops: List[Tuple[int]],
        css_x_boundary: List[int],
        x_gauge_products: List[int],
        css_z_gauge_ops: List[Tuple[int]],
        css_z_stabilizer_ops: List[Tuple[int]],
        css_z_boundary: List[int],
        z_gauge_products: List[int],
        basis: str,
        layer_types: List[str],
        x_gauge_outcomes: List[List[int]],
        z_gauge_outcomes: List[List[int]],
        final_outcomes: List[int],
    ) -> Tuple[List[List[int]], List[Tuple[int, Tuple[int]]]]:
        """Identify highlighted vertices in the decoding graph for an outcome.

        Gauge operator measurement outcomes are lists of integers 0, 1.
        """
        if basis == "z":
            gauge_outcomes = z_gauge_outcomes
            gauges = css_z_gauge_ops
            stabilizers = css_z_stabilizer_ops
            boundary = css_z_boundary
            gauge_products = z_gauge_products
        elif basis == "x":
            gauge_outcomes = x_gauge_outcomes
            gauges = css_x_gauge_ops
            stabilizers = css_x_stabilizer_ops
            boundary = css_x_boundary
            gauge_products = x_gauge_products
        final_gauges = []
        for supp in gauges:
            parity = 0
            for i in supp:
                if i != -1:  # supp can contain -1 if no qubit at that site
                    parity += final_outcomes[i]
            final_gauges.append(parity % 2)
        gauge_outcomes.append(final_gauges)

        highlighted = []
        # Now iterate over the layers and look at appropriate
        # syndrome differences
        for i, ltype in enumerate(layer_types):
            if ltype == "g":
                # compare current and past gauge measurements
                # if a bit differs, the vertex (i, g) is highlighted
                for j, gauge_op in enumerate(gauges):
                    if (i == 0 and gauge_outcomes[i][j] == 1) or (
                        i > 0 and gauge_outcomes[i][j] != gauge_outcomes[i - 1][j]
                    ):
                        highlighted.append((i, tuple(gauge_op)))
            elif ltype == "s":
                # compare current and past stabilizer measurements
                # if a bit differs, the vertex (i, s) is highlighted
                for j, stab_op in enumerate(stabilizers):
                    outcome = 0
                    prior_outcome = 0
                    for k in gauge_products[j]:
                        outcome += gauge_outcomes[i][k]
                        if i > 0:
                            prior_outcome += gauge_outcomes[i - 1][k]
                    outcome %= 2
                    prior_outcome %= 2
                    if outcome != prior_outcome:
                        highlighted.append((i, tuple(stab_op)))
        logging.debug("|highlighted| = %d", len(highlighted))
        # If the total number of highlighted vertices is odd,
        # add a single special highlighted vertex at the boundary
        if len(highlighted) % 2 == 1:
            highlighted.append((0, tuple(boundary[0])))
        logging.debug("highlighted = %s", highlighted)
        return gauge_outcomes, highlighted

    def _compute_matching(
        self,
        idxmap: Dict[Tuple[int, List[int]], int],
        length: Dict[int, Dict[int, int]],
        highlighted: List[Tuple[int, Tuple[int]]],
    ) -> Set[Tuple[int, int]]:
        """Compute a min. weight perfect matching of highlighted vertices.

        highlighted is a list of highlighted vertices given as tuples
        (t, qubit_set).
        Return the matching.
        """
        gm = rx.PyGraph(multigraph=False)  # matching graph
        idx = 0  # vertex index in matching graph
        midxmap = {}  # map from (t, qubit_tuple) to vertex index
        for v in highlighted:
            gm.add_node({"dvertex": v})
            midxmap[v] = idx
            idx += 1
        for i, high_i in enumerate(highlighted):
            for j in range(i + 1, len(highlighted)):
                vi = midxmap[high_i]
                vj = midxmap[highlighted[j]]
                vip = idxmap[high_i]
                vjp = idxmap[highlighted[j]]
                gm.add_edge(vi, vj, {"weight": -length[vip][vjp]})

        def weight_fn(edge):
            return int(edge["weight"])

        matching = rx.max_weight_matching(gm, max_cardinality=True, weight_fn=weight_fn)
        return matching

    def _error_chain_for_vertex_path(
        self, graph: rx.PyGraph, vertex_path: List[int]
    ) -> Tuple[Set[int], Set[int]]:
        """Return a chain of qubit and measurement errors for a vertex path.

        Examine the edges along the path to extract the error chain.
        Store error chains as sets and merge using symmetric difference.
        The vertex_path is a list of retworkx node indices.
        """
        qubit_errors = set([])
        measurement_errors = set([])
        logging.debug("_error_chain_for_vertex_path %s", vertex_path)
        for i in range(len(vertex_path) - 1):
            v0 = vertex_path[i]
            v1 = vertex_path[i + 1]
            if graph.get_edge_data(v0, v1)["measurement_error"] == 1:
                measurement_errors ^= set(
                    [(graph.nodes()[v0]["time"], tuple(graph.nodes()[v0]["qubits"]))]
                )
            qubit_errors ^= set(graph.get_edge_data(v0, v1)["qubits"])
            logging.debug(
                "_error_chain_for_vertex_path q = %s, m = %s",
                qubit_errors,
                measurement_errors,
            )
        return qubit_errors, measurement_errors

    def _compute_error_correction(
        self,
        gin: rx.PyGraph,
        idxmap: Dict[Tuple[int, List[int]], int],
        paths: Dict[int, Dict[int, List[int]]],
        matching,
        highlighted,
        export: bool = False,
        filename: str = "graphFile.json",
    ) -> Tuple[Set[int], Set[int]]:
        """Compute the qubit and measurement corrections.

        gin is the decoding graph.
        idxmap maps (t, qubit_idx) to vertex index.
        paths is all pairs shortest paths between vertices.
        matching is the perfect matching computed by _compute_matching.
        highlighted is the list of highlighted vertices computed by
        _highlighted_vertices.
        if export is True, write the decoding graph to filename with
        the corrective paths highlighted.
        Returns a tuple of sets, (qubit_errors, measurement_errors) where
        qubit_errors contains the indices of qubits with errors and
        measurement_errors contains tuples (t, qubit_set) indicating the
        failed measurement.
        """
        used_paths = []
        qubit_errors = set([])
        measurement_errors = set([])
        for p in matching:
            v0 = idxmap[highlighted[p[0]]]
            v1 = idxmap[highlighted[p[1]]]
            # Use the shortest paths between the matched vertices to
            # identify all of the qubits in the error chains
            path = paths[v0][v1]
            q, m = self._error_chain_for_vertex_path(gin, path)
            # Add the error chains modulo two to get the total correction
            # (uses set symmetric difference)
            qubit_errors ^= q
            measurement_errors ^= m
            used_paths.append(path)
        if export:
            self._highlight_vertex_paths_export_json(gin, used_paths, filename)
        return qubit_errors, measurement_errors

    def _highlight_vertex_paths_export_json(
        self, gin: rx.PyGraph, paths: List[List[int]], filename: str
    ):
        """Highlight the vertex paths and export JSON to file.

        paths is a list of vertex paths, each given as a list of
        vertex indices in the decoding graph.
        """
        graph = deepcopy(gin)
        for path in paths:
            # Highlight the endpoints of the path
            for i in [0, -1]:
                graph.nodes()[i]["highlighted"] = True
            # Highlight the edges along the path
            for i in range(len(path) - 1):
                idx = list(graph.graph.edge_list()).index((i, i + 1))
                graph.edges[idx]["highlighted"] = True
        self._export_json(graph, filename)

    def _export_json(self, graph: nx.Graph, filename: str):
        """Export a JSON formatted file with the graph data."""
        with open(filename, "w", encoding="utf-8") as fp:
            # The sympy fields are not serializable, so we
            # include the default function that casts to a string
            json.dump(nx.node_link_data(graph), fp, indent=4, default=str)
        fp.close()


def ret2net(graph):
    """Convert retworkx graph to equivalent networx graph."""
    nx_graph = nx.Graph()
    for j, node in enumerate(graph.nodes()):
        nx_graph.add_node(j)
        for k, v in node.items():
            nx.set_node_attributes(nx_graph, {j: v}, k)
    for j, (n0, n1) in enumerate(graph.edge_list()):
        nx_graph.add_edge(n0, n1)
        for k, v in graph.edges()[j].items():
            nx.set_edge_attributes(nx_graph, {(n0, n1): v}, k)
    return nx_graph
