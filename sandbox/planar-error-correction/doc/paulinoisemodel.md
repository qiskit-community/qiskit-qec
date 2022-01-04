Module paulinoisemodel
======================
Pauli circuit-level noise model.

Classes
-------

`PauliNoiseModel(fromdict=None)`
:   Encapsulate a circuit-level Pauli noise model.
    
    Create a new Pauli noise model.
    
    Optionally provide a dictionary for the model.

    ### Methods

    `add_operation(self, name, paulichanneldict)`
    :   Add a new faulty operation.
        
        name = string label for operation
        paulichanneldict = dictionary {paulistring: weight, ...}
        
        Each paulistring contains "i", "x", "y", and "z".
        The weights do not need to be normalized.

    `as_aer_noise_model(self)`
    :   Return corresponding Aer noise model.

    `get_error_probability(self, name)`
    :   Get the error probability of an operation.

    `get_operations(self)`
    :   Return the list of defined names.

    `get_pauli_error_probability(self, name, paulistring)`
    :   Get the error probability for a particular term.
        
        name = string label for operation
        paulistring = string containing only "i", "x", "y", and "z".

    `get_pauli_error_types(self)`
    :   Return a dict of error types for each operation.

    `get_pauli_weight(self, name, paulistring)`
    :   Get the weight for a particular term.
        
        name = string label for operation
        paulistring = string containing only "i", "x", "y", and "z".

    `set_error_probability(self, name, p)`
    :   Assign an error probability to an operation.

    `set_scale_factor(self, name, factor)`
    :   Assign a scaling factor to an operation.

    `set_scaled_error_probabilities(self, p)`
    :   Scale and assign error probabilities to operations.
        
        Only operations with a scale factor property will be changed.

    `to_dict(self)`
    :   Return corresponding dictionary.