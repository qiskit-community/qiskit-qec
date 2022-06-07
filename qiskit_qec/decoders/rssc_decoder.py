"""Object to construct decoders and decode the RSSC."""

import logging
import networkx as nx

from qiskit_qec.decoder.circuit_matching_decoder  import CircuitModelMatchingDecoder
from qiskit_qec.utils.indexer import Indexer


class RSSCDecoder(CircuitModelMatchingDecoder):
    """Decoder for the rotated subsystem surface code."""

    def __init__(self, code, circuit, model, round_schedule):
        """Create a decoder object."""
        # Sum the total number of bits per round
        self.bits_per_round = 0
        self.round_schedule = round_schedule
        if not isinstance(self.round_schedule, str):
            raise Exception("expected round_schedule to be a string")
        if set(self.round_schedule) > set("xz"):
            raise Exception("expected round schedule of 'x', 'z' chars")
        for step in self.round_schedule:
            if step == 'z':
                self.bits_per_round += len(code.z_gauges)
            elif step == 'x':
                self.bits_per_round += len(code.x_gauges)
        super().__init__(code, circuit, model, config)

    def _partition_outcomes(self, rounds, round_schedule, outcomes):
        """Process the raw outcome and return results.

        rounds = number of rounds
        round_schedule = string of z and x characters
        outcomes = list of 0, 1 outcomes

        Return lists x_gauge_outcomes, z_gauge_outcomes, final_outcomes.
        """
        x_gauge_outcomes = []
        z_gauge_outcomes = []
        for r in range(rounds):
            bits_into_round = 0
            for s in range(len(round_schedule)):
                if round_schedule[s] == 'z':
                    z_gauge_outcomes.append(
                            outcomes[r*self.bits_per_round + bits_into_round:
                                     r*self.bits_per_round + bits_into_round
                                     + len(self.code.z_gauges)])
                    bits_into_round += len(self.code.z_gauges)
                elif round_schedule[s] == 'x':
                    x_gauge_outcomes.append(
                            outcomes[r*self.bits_per_round + bits_into_round:
                                     r*self.bits_per_round + bits_into_round
                                     + len(self.code.x_gauges)])
                    bits_into_round += len(self.code.x_gauges)
        final_outcomes = outcomes[-self.code.n:]
        return x_gauge_outcomes, z_gauge_outcomes, final_outcomes

    def _decoding_graph(self, code, rounds, round_schedule, basis):
        """Construct the decoding graph for the given basis.

        This method assumes a model where edge weights are uniformly 1.

        Returns a tuple (idxmap, node_layers, G) where idxmap is a dict
        mapping tuples (t, qubit_set) to integer vertex indices in the
        decoding graph G. The list node_layers contains lists of nodes
        for each time step.
        """
        G = nx.Graph()
        gauges = []
        stabilizers = []
        boundary = []
        if basis == 'z':
            gauges = self.code.z_gauges
            stabilizers = self.code.z_stabilizers
            boundary = self.code.z_boundary
        elif basis == 'x':
            gauges = self.code.x_gauges
            stabilizers = self.code.x_stabilizers
            boundary = self.code.x_boundary
        layer_types = self.layer_types

        logging.debug('_decoding_graph gauges %s' % gauges)
        logging.debug('_decoding_graph stabilizers %s' % stabilizers)
        logging.debug('_decoding_graph boundary %s' % boundary)
        logging.debug('_decoding_graph layer_types %s' % layer_types)

        self.pymatching_indexer = Indexer()

        # Construct the decoding graph
        idx = 0
        idxmap = {}
        node_layers = []
        for t in range(len(layer_types)):
            layer = layer_types[t]
            # Add vertices at time t
            node_layer = []
            if layer == 'g':
                all_z = gauges
            elif layer == 's':
                all_z = stabilizers
            for g in all_z:
                G.add_node(idx, time=t, qubits=g, highlighted=False)
                logging.debug('node %d t=%d %s' % (idx, t, g))
                idxmap[(t, tuple(g))] = idx
                node_layer.append(idx)
                idx += 1
            for g in boundary:
                # Add optional is_boundary property for pymatching
                G.add_node(idx, time=t, qubits=g, highlighted=False,
                           is_boundary=True)
                logging.debug('boundary %d t=%d %s' % (idx, t, g))
                idxmap[(t, tuple(g))] = idx
                node_layer.append(idx)
                idx += 1
            node_layers.append(node_layer)
            if layer == 'g':
                all_z = gauges + boundary
            elif layer == 's':
                all_z = stabilizers + boundary
            # Add space-like edges at time t
            # By construction, the qubit sets of any pair of vertices at time
            # t intersect on at most one qubit.
            # (this assumption is specific to the RSSC)
            # If they intersect, we add an edge and label it by the common
            # qubit.
            # Space-like edges do not correspond to syndrome errors, so the
            # syndrome property is an empty list.
            for i in range(len(all_z)):
                for j in range(i+1, len(all_z)):
                    g = all_z[i]
                    h = all_z[j]
                    u = list(set(g).intersection(set(h)))
                    if -1 in u:
                        u.remove(-1)
                    if len(u) > 0:
                        # Include properties for use with pymatching:
                        # qubit_id is an integer or set of integers
                        # weight is a floating point number
                        # error_probability is a floating point number
                        G.add_edge(idxmap[(t, tuple(g))],
                                   idxmap[(t, tuple(h))],
                                   qubits=u, measurement_error=0,
                                   weight=1, highlighted=False,
                                   qubit_id=self.pymatching_indexer[u[0]],
                                   error_probability=0.01)
                        logging.debug('spacelike t=%d (%s, %s)'
                                      % (t, g, h))

            # Add boundary space-like edges
            for i in range(len(boundary)-1):
                g = boundary[i]
                h = boundary[i+1]
                # Include properties for use with pymatching:
                # qubit_id is an integer or set of integers
                # weight is a floating point number
                # error_probability is a floating point number
                G.add_edge(idxmap[(t, tuple(g))],
                           idxmap[(t, tuple(h))],
                           qubits=[], measurement_error=0,
                           weight=0, highlighted=False,
                           qubit_id=-1,
                           error_probability=0.0)
                logging.debug(
                    'spacelike boundary t=%d (%s, %s)'
                    % (t, g, h))

            # Add (space)time-like edges from t to t-1
            # By construction, the qubit sets of pairs of vertices at S and T
            # at times t-1 and t respectively
            # either (a) contain each other (S subset T or T subset S) and
            # |S|,|T|>1,
            # (b) intersect on a qubit, or (c) are disjoint.
            # In case (a), we add an edge that corresponds to a syndrome bit
            # error at time t-1.
            # In case (b), we add an edge that corresponds to a hook error,
            # i.e. a syndrome bit error at time t-1
            # together with an error on the common qubit.
            # In case (c), we do not add an edge.
            if t > 0:
                current_sets = gauges
                prior_sets = gauges
                if layer_types[t] == 's':
                    current_sets = stabilizers
                if layer_types[t-1] == 's':
                    prior_sets = stabilizers
                for g in current_sets:
                    for h in prior_sets:
                        u = list(set(g).intersection(set(h)))
                        if -1 in u:
                            u.remove(-1)
                        if len(u) > 0:
                            # Include properties for use with pymatching:
                            # qubit_id is an integer or set of integers
                            # weight is a floating point number
                            # error_probability is a floating point number
                            if len(u) == 1:
                                q_idx = self.pymatching_indexer[u[0]]
                                G.add_edge(idxmap[(t-1, tuple(h))],
                                           idxmap[(t, tuple(g))],
                                           qubits=u, measurement_error=1,
                                           weight=1, highlighted=False,
                                           qubit_id=q_idx,
                                           error_probability=0.01)
                                logging.debug('hook t=%d (%s, %s)'
                                              % (t, g, h))
                            else:
                                G.add_edge(idxmap[(t-1, tuple(h))],
                                           idxmap[(t, tuple(g))],
                                           qubits=[], measurement_error=1,
                                           weight=1, highlighted=False,
                                           qubit_id=-1,
                                           error_probability=0.01)
                                logging.debug('timelike t=%d (%s, %s)'
                                              % (t, g, h))
                # Add a single time-like edge between boundary vertices at
                # time t-1 and t
                G.add_edge(idxmap[(t-1, tuple(boundary[0]))],
                           idxmap[(t, tuple(boundary[0]))],
                           qubits=[], measurement_error=0,
                           weight=0, highlighted=False,
                           qubit_id=-1, error_probability=0.0)
                logging.debug('boundarylink t=%d (%s, %s)'
                              % (t, g, h))
        logging.debug('indexer %s' % self.pymatching_indexer)
        return idxmap, node_layers, G
