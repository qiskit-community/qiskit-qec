Module faultenumerator
======================
Quantum circuit fault path enumerator.

Classes
-------

`FaultEnumerator(circ, order=1, method='stabilizer', model=None, sim_seed=0)`
:   Enumerates faults in a circuit according to a noise model.
    
    Construct a fault enumerator object.
    
    circ = QuantumCircuit object
    order = number of faults to insert
    method = simulation method
    model = PauliNoiseModel object
    
    Faulty operations are identified by a string
    containing the Qiskit gate name or label
    if the latter exists.

    ### Methods

    `generate(self)`
    :   Generator to iterate over faults in a quantum circuit.
        
        Note: The bit order of outcome is reversed from Qiskit.
        
        Yields (index, labels, error, outcome).

    `generate_blocks(self, blocksize=10000)`
    :   Generator to iterate over sequences of faults in a quantum circuit.
        
        blocksize = number of faults to process per call
        
        Note: The bit order of outcome is reversed from Qiskit.
        
        Yields [(index, labels, error, outcome), ...] with approximately
        blocksize elements. May be slightly more per call, or less on the
        final call.