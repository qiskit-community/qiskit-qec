Module circuit_matching_decoder
===============================
Abstract object to construct matching decoders for circuit noise.

Classes
-------

`CircuitModelMatchingDecoder(code, circuit, model, config)`
:   Matching decoder for circuit noise.
    
    Create a decoder object.
    
    Decoder for one logical qubit with repeated syndrome measurements.
    
    code is a CSSCode containing the error-correcting code definition.
    
    circuit is a QuantumCircuit to build the decoder's graphical model
    for circuit noise.
    
    model is a PauliNoiseModel whose operations coincide with
    the operations in the circuit.
    
    config is a Config object containing parameter values:
    
        config['circuit']['basis']:  state is initialized and
            measured in the X or Z basis, 'x' or 'z'.
    
        config['circuit']['round_schedule']: Round schedule for
            X or Z gauge rounds. For example, 'zx' means a round
            is 1 Z gauge measurement followed by 1 X gauge
            measurement.
    
        config['circuit']['rounds']: Rounds is the number of
            repetitions of the round schedule between initialization
            and measurement. For example, 2 rounds with round
            schedule 'zx' means a total of 2 Z measurement cycles
            and 2 X measurement cycles.
    
        config['decoder']['method']: Method selects the decoder
            implementation. It can be either 'matching_pymatching'
            or 'matching_networkx'. The first selects the pymatching
            library and the second selects a pure python implementation.
    
        config['decoder']['uniform']: Boolean that, if True, turns
            off fault propagation through the circuit and instead
            uses equal edge weights for all edges in the decoding
            graph. Even if uniform is True, it is still necessary
            to call update_edge_weights at least once before
            calling process.

    ### Ancestors (in MRO)

    * abc.ABC

    ### Methods

    `export_decoding_graph_json(self, filename)`
    :   Write a file containing the decoding graph.

    `process(self, outcomes, export=False, filename='decoding.json')`
    :   Process a set of outcomes and return corrected final outcomes.
        
        Be sure to have called update_edge_weights for the
        noise parameters so that the edge weights are updated.
        
        The result is a list of code.n integers that are 0 or 1.
        These are the corrected values of the final transversal
        measurement in the basis given by self.basis.

    `update_edge_weights(self, model)`
    :   Evaluate the numerical edge weights and update graph data.
        
        For each edge in the decoding graph that has a "weight_poly"
        property, evaluate the polynomial at the given model parameters
        and set the corresponding "weight" property. Once this is done,
        recompute the shortest paths between pairs of vertices
        in the decoding graph.
        
        model is a PauliNoiseModel whose error probabilities have been
        previously assigned. The probabilities are then assigned to
        the variables in self.symbols.
        
        Updates properties of self.g.
        If not using pymatching, updates sets self.length and self.path.
        If using pymatching, constructs a pymatching object self.pymatching.
        
        Note that if self.uniform is True, it is still necessary to call
        this function to construct the matching decoder object (pymatching)
        or compute the shortest paths between vertices in the decoding
        graph (networkx).