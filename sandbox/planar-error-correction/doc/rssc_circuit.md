Module rssc_circuit
===================
Object to construct quantum circuits for the RSSC.

Classes
-------

`RSSCCircuit(rssc, config)`
:   Create quantum circuits for syndrome measurements.
    
    Specialized to rotated subsystem surface code.
    
    Create an object associated to a RSSC.
    
    rssc = RSSC code object
    config = Config object

    ### Methods

    `syndrome_measurement(self, rounds=None, round_schedule=None, basis=None, initial_state=None, logical_paulis=None)`
    :   Construct repeated syndrome measurement circuit.
        
        Method parameters will override the configuration that
        was provided to the constructor.
        
        rounds = number of repetitions of round_schedule (int)
        round_schedule = schedule of x/z gauge rounds (str)
        basis = initialization and measurement basis, x or z (str)
        initial_state = eigenvalue + or - of basis (str)
        logical_paulis = what logical Pauli to apply after each
          X or Z gauge round in round_schedule (str)
        
        Additional options from self.config:
        
        self.config['circuit']['schedule']: (str) If this equals
        "higgott-breuckmann", we will the schedule in
        Phys. Rev. X 11, 031039 (2021). If this equals
        "heavy-hex", we will use our own circuits.
        
        self.config['circuit']['barriers']: (bool) If this is True,
        insert barrier commands between steps of the syndrome circuit.
        
        self.config['circuit']['idles']: (bool) If this is True,
        insert identity gates at time steps where qubits are idle.
        Use a timing model where two-qubit gates and measurements
        have the same duration and single-qubit gates have negligible
        duration.
        
        self.config['circuit']['distinct_measurement_idle']: (bool) If
        this is True, insert an 'idm' labeled identity gate on the data
        qubits while the ancillas are measured, so the idle errors during
        measurement can be changed independently from idle errors elsewhere
        in the circuit.