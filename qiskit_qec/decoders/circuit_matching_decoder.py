"""Abstract object to construct CSS matching decoders for circuit noise."""

from copy import deepcopy, copy
from math import log
from abc import ABC, abstractmethod
from typing import List, Tuple, Dict, Set
import logging

import json
from sympy import Poly, symbols, Symbol
from qiskit import QuantumCircuit
import networkx as nx
from pymatching import Matching

from qiskit_qec.utils.indexer import Indexer
from qiskit_qec.noise.paulinoisemodel import PauliNoiseModel
from qiskit_qec.analysis.faultenumerator import FaultEnumerator

from qiskit_qec.decoders.temp_code_util import temp_syndrome


class CircuitModelMatchingDecoder(ABC):
    """Matching decoder for circuit noise."""

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
    ):  # TODO defaults?
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
            raise Exception("expected positive integer for blocks")
        self.round_schedule = round_schedule
        if set(self.round_schedule) > set("xyz"):
            raise Exception("expected round schedule of 'x', 'y', 'z' chars")
        self.basis = basis
        if not self.basis in ("x", "z"):
            raise Exception("expected basis to be 'x' or 'z'")

        self.uniform = uniform
        self.method = method
        if not self.method in ("networkx", "pymatching"):
            raise Exception("unsupported implementation")
        if self.method == "pymatching":
            self.usepymatching = True
        else:
            self.usepymatching = False

        self.z_gauge_products = self._gauge_products(
            self.css_z_stabilizer_ops, self.css_z_gauge_ops
        )
        self.x_gauge_products = self._gauge_products(
            self.css_x_stabilizer_ops, self.css_x_gauge_ops
        )
        self.layer_types = self._layer_types(self.blocks, self.round_schedule, self.basis)
        logging.debug("layer_types = %s", self.layer_types)

        self.idxmap, self.node_layers, self.graph, self.pymatching_indexer = self._decoding_graph(
            self.css_x_gauge_ops,
            self.css_x_stabilizer_ops,
            self.css_x_boundary,
            self.css_z_gauge_ops,
            self.css_z_stabilizer_ops,
            self.css_z_boundary,
            self.layer_types,
            self.basis,
        )
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

        self.length = {}  # recomputed in update_edge_weights
        self.path = {}  # recomputed in update_edge_weights
        self.pymatching = None  # constructed in update_edge_weights

    def _gauge_products(self, stabilizers: List[List[int]], gauges: List[List[int]]):
        """Record gauge products for each stabilizer.

        stabilizers = list of stabilizer operator supports
        gauges = list of gauge operator supports

        For each stabilizer operator, record the indices of the
        gauge operators in 'gauges' whose product is that stabilizer operator.

        Return the list of indices.
        """
        gauge_products = []
        for stab in stabilizers:
            is_contained = lambda x: set(x).intersection(set(stab)) == set(x)
            products = filter(is_contained, gauges)
            gauge_products.append(list(products))
        return gauge_products

    def _revise_decoding_graph(
        self,
        idxmap: Dict[Tuple[int, List[int]], int],
        graph: nx.Graph,
        edge_weight_polynomials: Dict[Tuple[int, Tuple[int]], Dict[Tuple[int, Tuple[int]], Poly]],
    ) -> nx.Graph:
        """Add edge weight polynomials to the decoding graph g and prune it.

        Update attribute "weight_poly" on decoding graph edges contained in
        edge_weight_polynomials. Remove all other edges from the decoding graph
        that have non-zero weight.
        """
        for s1, sub in edge_weight_polynomials.items():
            for s2, wpoly in sub.items():
                if s1 not in idxmap:
                    raise Exception(f"vertex {s1} not in decoding graph")
                if s2 not in idxmap:
                    raise Exception(f"vertex {s2} not in decoding graph")
                if not graph.has_edge(idxmap[s1], idxmap[s2]):
                    raise Exception("edge {s1} - {s2} not in decoding graph")
                    # TODO: new edges may be needed for hooks, but raise exception if we're surprised
                graph.edges[idxmap[s1], idxmap[s2]]["weight_poly"] = wpoly
        remove_list = []
        for edge in graph.edges(data=True):
            if "weight_poly" not in edge[2] and edge[2]["weight"] != 0:
                if (
                    "is_boundary" in graph.nodes(data=True)[edge[0]]
                    and graph.nodes(data=True)[edge[0]]["is_boundary"]
                ) or (
                    "is_boundary" in graph.nodes(data=True)[edge[1]]
                    and graph.nodes(data=True)[edge[1]]["is_boundary"]
                ):
                    # pymatching requires all qubit_ids from 0 to
                    # n-1 to appear as edge properties. If we
                    # remove too many edges, not all qubit_ids will
                    # be present. We solve this by keeping unweighted
                    # that are edges connected to the boundary. Instead
                    # we increase their weight to a large value.
                    edge[2]["weight"] = 1000000
                    logging.info("increase edge weight (%d, %d)", edge[0], edge[1])
                else:
                    # Otherwise, remove the edge
                    remove_list.append((edge[0], edge[1]))
                    logging.info("remove edge (%d, %d)", edge[0], edge[1])
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
        graph (networkx).
        """
        parameter_values = [model.get_error_probability(name) for name in self.parameters]
        if not self.uniform:
            if len(parameter_values) != len(self.parameters):
                raise Exception("wrong number of error rate parameters")
            symbol_list = [self.symbols[s] for s in self.parameters]
            assignment = dict(zip(symbol_list, parameter_values))
            logging.info("update_edge_weights %s", str(assignment))
            # P(chain) = \prod_i (1-p_i)^{1-l(i)}*p_i^{l(i)}
            #          \propto \prod_i ((1-p_i)/p_i)^{l(i)}
            # -log P(chain) \propto \sum_i -log[((1-p_i)/p_i)^{l(i)}]
            # p_i is the probability that edge i carries an error
            # l(i) is 1 if the link belongs to the chain and 0 otherwise
            for edge in self.graph.edges(data=True):
                if "weight_poly" in edge[2]:
                    logging.info(
                        "update_edge_weights (%d, %d) %s",
                        edge[0],
                        edge[1],
                        str(edge[2]["weight_poly"]),
                    )
                    restriction = {x: assignment[x] for x in edge[2]["weight_poly"].gens}
                    p = edge[2]["weight_poly"].eval(restriction).evalf()
                    assert p < 0.5, "edge flip probability too large"
                    edge[2]["weight"] = log((1 - p) / p)
                    edge[2]["error_probability"] = p
        if not self.usepymatching:
            # Recompute the shortest paths between pairs of vertices
            # in the decoding graph.
            self.length = dict(nx.all_pairs_dijkstra_path_length(self.graph))
            self.path = dict(nx.all_pairs_dijkstra_path(self.graph))  # note "cutoff" argument
        else:
            self.pymatching = Matching(self.graph)

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

            x_gauge_outcomes, z_gauge_outcomes, final_outcomes = self._partition_outcomes(
                blocks, round_schedule, outcome
            )

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
                raise Exception("too many highlighted vertices for a " + "single fault event")
            if len(highlighted) == 1:  # _highlighted_vertices highlights the boundary
                raise Exception("only one highlighted vertex for a " + "single fault event")
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
        pass

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
                self.graph, self.idxmap, self.path, matching, highlighted, export, filename
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
                raise Exception("Did you call update_edge_weights?") from attrib_error
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
            raise Exception("decoder failure: syndrome should be trivial!")
        return corrected_outcomes

    def export_decoding_graph_json(self, filename: str):
        """Write a file containing the decoding graph."""
        self._export_json(self.graph, filename)

    def _layer_types(self, blocks: int, round_schedule: str, basis: str) -> List[str]:
        """Return a list of decoding graph layer types.

        The entries are 'g' for gauge and 's' for stabilizer.
        """
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
        return layer_types

    def _decoding_graph(
        self,
        css_x_gauge_ops: List[Tuple[int]],
        css_x_stabilizer_ops: List[Tuple[int]],
        css_x_boundary: List[Tuple[int]],
        css_z_gauge_ops: List[Tuple[int]],
        css_z_stabilizer_ops: List[Tuple[int]],
        css_z_boundary: List[Tuple[int]],
        layer_types: List[str],
        basis: str,
    ) -> Tuple[Dict[Tuple[int, List[int]], int], List[List[int]], nx.Graph, Indexer]:
        """Construct the decoding graph for the given basis.

        This method sets edge weights all to 1 and is based on
        computing intersections of operator supports.

        Returns a tuple (idxmap, node_layers, G, pymatching_indexer)
        where idxmap is a dict
        mapping tuples (t, qubit_set) to integer vertex indices in the
        decoding graph G. The list node_layers contains lists of nodes
        for each time step.
        """
        graph = nx.Graph()
        gauges = []
        stabilizers = []
        boundary = []
        if basis == "z":
            gauges = css_z_gauge_ops
            stabilizers = css_z_stabilizer_ops
            boundary = css_z_boundary
        elif basis == "x":
            gauges = css_x_gauge_ops
            stabilizers = css_x_stabilizer_ops
            boundary = css_x_boundary

        pymatching_indexer = Indexer()

        # Construct the decoding graph
        idx = 0  # vertex index counter
        idxmap = {}  # map from vertex data (t, qubits) to vertex index
        node_layers = []
        for t, layer in enumerate(layer_types):
            # Add vertices at time t
            node_layer = []
            if layer == "g":
                all_z = gauges
            elif layer == "s":
                all_z = stabilizers
            for g in all_z:
                graph.add_node(idx, time=t, qubits=g, highlighted=False)
                logging.debug("node %d t=%d %s", idx, t, g)
                idxmap[(t, tuple(g))] = idx
                node_layer.append(idx)
                idx += 1
            for g in boundary:
                # Add optional is_boundary property for pymatching
                graph.add_node(idx, time=t, qubits=g, highlighted=False, is_boundary=True)
                logging.debug("boundary %d t=%d %s", idx, t, g)
                idxmap[(t, tuple(g))] = idx
                node_layer.append(idx)
                idx += 1
            node_layers.append(node_layer)
            if layer == "g":
                all_z = gauges + boundary
            elif layer == "s":
                all_z = stabilizers + boundary
            # Add space-like edges at time t
            # The qubit sets of any pair of vertices at time
            # t can intersect on multiple qubits.
            # If they intersect, we add an edge and label it by
            # one of the common qubits. This makes an assumption
            # that the intersection operator is equivalent to a single
            # qubit operator modulo the gauge group.
            # Space-like edges do not correspond to syndrome errors, so the
            # syndrome property is an empty list.
            for i, g in enumerate(all_z):
                for j in range(i + 1, len(all_z)):
                    h = all_z[j]
                    u = list(set(g).intersection(set(h)))
                    if -1 in u:
                        u.remove(-1)
                    if len(u) > 0:
                        # Include properties for use with pymatching:
                        # qubit_id is an integer or set of integers
                        # weight is a floating point number
                        # error_probability is a floating point number
                        graph.add_edge(
                            idxmap[(t, tuple(g))],
                            idxmap[(t, tuple(h))],
                            qubits=[u[0]],
                            measurement_error=0,
                            weight=1,
                            highlighted=False,
                            qubit_id=pymatching_indexer[u[0]],
                            error_probability=0.01,
                        )
                        logging.debug("spacelike t=%d (%s, %s)", t, g, h)
                        logging.debug(" qubits %s qubit_id %s", [u[0]], pymatching_indexer[u[0]])

            # Add boundary space-like edges
            for i in range(len(boundary) - 1):
                g = boundary[i]
                h = boundary[i + 1]
                # Include properties for use with pymatching:
                # qubit_id is an integer or set of integers
                # weight is a floating point number
                # error_probability is a floating point number
                graph.add_edge(
                    idxmap[(t, tuple(g))],
                    idxmap[(t, tuple(h))],
                    qubits=[],
                    measurement_error=0,
                    weight=0,
                    highlighted=False,
                    qubit_id=-1,
                    error_probability=0.0,
                )
                logging.debug("spacelike boundary t=%d (%s, %s)", t, g, h)

            # Add (space)time-like edges from t to t-1
            # By construction, the qubit sets of pairs of vertices at S and T
            # at times t-1 and t respectively
            # either (a) contain each other (S subset T or T subset S) and
            # |S|,|T|>1,
            # (b) intersect on one or more qubits, or (c) are disjoint.
            # In case (a), we add an edge that corresponds to a syndrome bit
            # error at time t-1.
            # In case (b), we add an edge that corresponds to a spacetime hook
            # error, i.e. a syndrome bit error at time t-1
            # together with an error on one of the common qubits. Again
            # this makes an assumption that all such errors are equivalent.
            # In case (c), we do not add an edge.
            # Important: some space-like hooks are not accounted for.
            # They can have longer paths between non-intersecting operators.
            # We will account for these in _revise_decoding_graph if needed.
            if t > 0:
                current_sets = gauges
                prior_sets = gauges
                if layer_types[t] == "s":
                    current_sets = stabilizers
                if layer_types[t - 1] == "s":
                    prior_sets = stabilizers
                for g in current_sets:
                    for h in prior_sets:
                        u = list(set(g).intersection(set(h)))
                        if -1 in u:
                            u.remove(-1)
                        if len(u) > 0:  # not Case (c)
                            # Include properties for use with pymatching:
                            # qubit_id is an integer or set of integers
                            # weight is a floating point number
                            # error_probability is a floating point number
                            # Case (a)
                            if set(u) == set(h) or set(u) == set(g):
                                graph.add_edge(
                                    idxmap[(t - 1, tuple(h))],
                                    idxmap[(t, tuple(g))],
                                    qubits=[],
                                    measurement_error=1,
                                    weight=1,
                                    highlighted=False,
                                    qubit_id=-1,
                                    error_probability=0.01,
                                )
                                logging.debug("timelike t=%d (%s, %s)", t, g, h)
                            else:  # Case (b)
                                q_idx = pymatching_indexer[u[0]]
                                graph.add_edge(
                                    idxmap[(t - 1, tuple(h))],
                                    idxmap[(t, tuple(g))],
                                    qubits=[u[0]],
                                    measurement_error=1,
                                    weight=1,
                                    highlighted=False,
                                    qubit_id=q_idx,
                                    error_probability=0.01,
                                )
                                logging.debug("spacetime hook t=%d (%s, %s)", t, g, h)
                                logging.debug(" qubits %s qubit_id %s", [u[0]], q_idx)
                # Add a single time-like edge between boundary vertices at
                # time t-1 and t
                graph.add_edge(
                    idxmap[(t - 1, tuple(boundary[0]))],
                    idxmap[(t, tuple(boundary[0]))],
                    qubits=[],
                    measurement_error=0,
                    weight=0,
                    highlighted=False,
                    qubit_id=-1,
                    error_probability=0.0,
                )
                logging.debug("boundarylink t=%d (%s, %s)", t, g, h)
        logging.debug("indexer %s", pymatching_indexer)
        return idxmap, node_layers, graph, pymatching_indexer

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
        for g in gauges:
            parity = 0
            for i in g:
                if i != -1:  # g can contain -1 if no qubit at that site
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
        gm = nx.Graph()  # matching graph
        idx = 0  # vertex index in matching graph
        midxmap = {}  # map from (t, qubit_tuple) to vertex index
        for v in highlighted:
            gm.add_node(idx, dvertex=v)
            midxmap[v] = idx
            idx += 1
        for i, high_i in enumerate(highlighted):
            for j in range(i + 1, len(highlighted)):
                vi = midxmap[high_i]
                vj = midxmap[highlighted[j]]
                vip = idxmap[high_i]
                vjp = idxmap[highlighted[j]]
                gm.add_edge(vi, vj, weight=-length[vip][vjp])
        matching = nx.max_weight_matching(gm, maxcardinality=True, weight="weight")
        return matching

    def _error_chain_for_vertex_path(
        self, graph: nx.Graph, vertex_path: List[int]
    ) -> Tuple[Set[int], Set[int]]:
        """Return a chain of qubit and measurement errors for a vertex path.

        Examine the edges along the path to extract the error chain.
        Store error chains as sets and merge using symmetric difference.
        The vertex_path is a list of networkx node indices.
        """
        qubit_errors = set([])
        measurement_errors = set([])
        logging.debug("_error_chain_for_vertex_path %s", vertex_path)
        for i in range(len(vertex_path) - 1):
            v0 = vertex_path[i]
            v1 = vertex_path[i + 1]
            if graph[v0][v1]["measurement_error"] == 1:
                measurement_errors ^= set(
                    [(graph.nodes[v0]["time"], tuple(graph.nodes[v0]["qubits"]))]
                )
            qubit_errors ^= set(graph[v0][v1]["qubits"])
            logging.debug(
                "_error_chain_for_vertex_path q = %s, m = %s", qubit_errors, measurement_errors
            )
        return qubit_errors, measurement_errors

    def _compute_error_correction(
        self,
        gin: nx.Graph,
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
        self, gin: nx.Graph, paths: List[List[int]], filename: str
    ):
        """Highlight the vertex paths and export JSON to file.

        paths is a list of vertex paths, each given as a list of
        vertex indices in the decoding graph.
        """
        g = deepcopy(gin)
        for path in paths:
            # Highlight the endpoints of the path
            nx.set_node_attributes(
                g, {path[0]: {"highlighted": True}, path[-1]: {"highlighted": True}}
            )
            # Highlight the edges along the path
            for i in range(len(path) - 1):
                nx.set_edge_attributes(g, {(path[i], path[i + 1]): {"highlighted": True}})
        self._export_json(g, filename)

    def _export_json(self, g: nx.Graph, filename: str):
        """Export a JSON formatted file with the graph data."""
        with open(filename, "w", encoding="utf-8") as fp:
            # The sympy fields are not serializable, so we
            # include the default function that casts to a string
            json.dump(nx.node_link_data(g), fp, indent=4, default=str)
        fp.close()
