Module hhc_circuit
==================
Object to construct quantum circuits for the HHC.

Classes
-------

`HHCCircuit(hhc, config)`
:   Create quantum circuits for syndrome measurements.
    
    Specialized to heavy hexagon subsystem code.
    
    Create an object associated to a HHC.
    
    hhc = HHC object with code definition
    config = Config object with circuit parameters

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
        "heavy-hex-yoder", we will use Ted Yoder's d=3 schedule.
        
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
        
        self.config['circuit']['group_meas']: (bool) If this is True,
        group right and left flag measurements, and group final syndrome
        and data qubit measurements.
        
        self.config['circuit']['num_initialize']: (None or int) If None,
        initialization happens in the beginning of gauge operation.
        If int, initialize all qubits (datas, syndromes, and flags) this
        number of times, and reset qubit after measurement rather than
        initialize before each gauge operations, except for final measurement.