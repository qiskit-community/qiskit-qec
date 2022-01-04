Module rssc
===========
Module to define the rotated subsystem surface code.

Classes
-------

`RSSC(distance)`
:   Create a new rotated subsystem surface code (RSSC).
    
    The X and Z gauge operator lists are given as lists of supports.
    There is a consistent qubit ordering of these lists so that
    other modules can construct gate schedules for circuits.
    
    Initialize data for RSSC.

    ### Ancestors (in MRO)

    * css_code.CSSCode

    ### Methods

    `to_index(self, row, col, d)`
    :   Map a coordinate (row, col) to a qubit index.
        
        Qubits are indexed from left to right across each row,
        beginning in the top row. The qubits on the faces of
        the lattice are indexed starting from d*d in the same way.