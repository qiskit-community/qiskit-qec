"""Abstract object for matching decoders for CSS codes and circuit noise."""

import logging
from abc import ABC, abstractmethod
from copy import copy
from math import log
from typing import Dict, List, Tuple
from sympy import Poly, Symbol, symbols

import rustworkx as rx
from qiskit import QuantumCircuit
from qiskit_qec.analysis.faultenumerator import FaultEnumerator
from qiskit_qec.decoders.decoding_graph import CSSDecodingGraph, DecodingGraph
from qiskit_qec.utils import DecodingGraphEdge
from qiskit_qec.decoders.pymatching_matcher import PyMatchingMatcher
from qiskit_qec.decoders.rustworkx_matcher import RustworkxMatcher
from qiskit_qec.decoders.temp_code_util import temp_gauge_products, temp_syndrome
from qiskit_qec.exceptions import QiskitQECError
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel


class CircuitModelMatchingDecoder(ABC):
    """Matching decoder for circuit noise."""

    METHOD_RETWORKX: str = "rustworkx"
    METHOD_PYMATCHING: str = "pymatching"
    AVAILABLE_METHODS = {METHOD_RETWORKX, METHOD_PYMATCHING}

    def __init__(
        self,
        n: int,
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
        decoding_graph: DecodingGraph = None,
        annotate: bool = False,
    ):
        """Create a matching decoder.

        Specialized to (subsystem) CSS codes encoding one logical qubit,
        and quantum circuits that prepare and measure in the Z or X basis
        and do repeated Pauli measurements.

        n : block size of quantum code
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
        annotate : for rustworkx method, compute self.matcher.annotated_graph
        """
        self.n = n
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
        self.matcher = None
        if self.method == self.METHOD_PYMATCHING:
            self.matcher = PyMatchingMatcher()
        else:
            self.matcher = RustworkxMatcher(annotate)

        self.z_gauge_products = temp_gauge_products(self.css_z_stabilizer_ops, self.css_z_gauge_ops)
        self.x_gauge_products = temp_gauge_products(self.css_x_stabilizer_ops, self.css_x_gauge_ops)

        if decoding_graph:
            (
                self.idxmap,
                self.node_layers,
                self.graph,
                self.layer_types,
            ) = self._process_graph(decoding_graph.graph, blocks, round_schedule, basis)
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
                self.layer_types,
            ) = (dg.idxmap, dg.node_layers, dg.graph, dg.layer_types)

        logging.debug("layer_types = %s", self.layer_types)

        self.ridxmap = {v: k for k, v in self.idxmap.items()}

        self.circuit = circuit
        self.event_map = {}
        self.parameters = model.get_operations()
        self.edge_weight_polynomials = {}
        self.symbols = None
        if not self.uniform:
            fe = FaultEnumerator(circuit, order=1, method="propagator", model=self.model)
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

    @staticmethod
    def _process_graph(
        graph: DecodingGraph, blocks: int, round_schedule: str, basis: str
    ) -> Tuple[Dict[Tuple[int, List[int]], int], List[List[int]], DecodingGraph, List[str]]:
        """Process a decoding graph to add required attributes."""

        # symmetrize hook errors
        for j, edge in enumerate(graph.edges()):
            n0, n1 = graph.edge_list()[j]
            source = graph.nodes()[n0]
            target = graph.nodes()[n1]
            if source.time != target.time:
                if source.is_boundary == target.is_boundary == False:
                    new_source = copy(source)
                    new_source.time = target.time
                    nn0 = graph.nodes().index(new_source)
                    new_target = copy(target)
                    new_target.time = source.time
                    nn1 = graph.nodes().index(new_target)
                    graph.add_edge(nn0, nn1, edge)

        edges_to_remove = []
        for j, edge in enumerate(graph.edges()):
            n0, n1 = graph.edge_list()[j]
            source = graph.nodes()[n0]
            target = graph.nodes()[n1]

            # add the required attributes
            # highlighted', 'measurement_error','qubit_id' and 'error_probability'
            edge.properties["highlighted"] = False
            edge.properties["measurement_error"] = int(source.time != target.time)

            # make it so times of boundary/boundary nodes agree
            if source.is_boundary and not target.is_boundary:
                if source.time != target.time:
                    new_source = copy(source)
                    new_source.time = target.time
                    n = graph.add_node(new_source)
                    edge.properties["measurement_error"] = 0
                    edges_to_remove.append((n0, n1))
                    graph.add_edge(n, n1, edge)

        # remove old boundary/boundary nodes
        for n0, n1 in edges_to_remove:
            graph.remove_edge(n0, n1)

        for n0, source in enumerate(graph.nodes()):
            for n1, target in enumerate(graph.nodes()):
                # add weightless nodes connecting different boundaries
                if source.time == target.time:
                    if source.is_boundary and target.is_boundary:
                        if source.qubits != target.qubits:
                            edge = DecodingGraphEdge(
                                weight=0,
                                qubits=list(set(source.qubits).intersection((set(target.qubits)))),
                            )
                            edge.properties["highlighted"] = False
                            edge.properties["measurement_error"] = 0
                            if (n0, n1) not in graph.edge_list():
                                graph.add_edge(n0, n1, edge)

                # connect one of the boundaries at different times
                if target.time == (source.time or 0) + 1:
                    if source.qubits == target.qubits == [0]:
                        edge = DecodingGraphEdge(weight=0, qubits=[])
                        edge.properties["highlighted"] = False
                        edge.properties["measurement_error"] = 0
                        if (n0, n1) not in graph.edge_list():
                            graph.add_edge(n0, n1, edge)

        # symmetrize edges
        for j, edge in enumerate(graph.edges()):
            n0, n1 = graph.edge_list()[j]
            if (n1, n0) not in graph.edge_list():
                graph.add_edge(n1, n0, edge)

        idxmap = {}
        for n, node in enumerate(graph.nodes()):
            idxmap[node.time, tuple(node.qubits)] = n

        node_layers = []
        for node in graph.nodes():
            time = node.time or 0
            if len(node_layers) < time + 1:
                node_layers += [[]] * (time + 1 - len(node_layers))
            node_layers[time].append(node.qubits)

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

        return idxmap, node_layers, graph, layer_types

    @staticmethod
    def _revise_decoding_graph(
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
                    # TODO: new edges may be needed for hooks, but raise exception for now
                    raise QiskitQECError("edge {s1} - {s2} not in decoding graph")
                data = graph.get_edge_data(idxmap[s1], idxmap[s2])
                data.properties["weight_poly"] = wpoly
        remove_list = []
        for source, target in graph.edge_list():
            edge_data = graph.get_edge_data(source, target)
            if "weight_poly" not in edge_data.properties and edge_data.weight != 0:
                # Remove the edge
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

        Updates properties of matcher.

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
                edge_data = self.graph.get_edge_data(source, target).properties
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
        self.matcher.preprocess(self.graph)

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

    def process(self, outcomes: List[int]) -> List[int]:
        """Process a set of outcomes and return corrected final outcomes.

        Be sure to have called update_edge_weights for the
        noise parameters.

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

        qubit_errors, _ = self.matcher.find_errors(self.graph, self.idxmap, highlighted)

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

    @staticmethod
    def _highlighted_vertices(
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
