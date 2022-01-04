Module hhc
==========
Module to define the heavy hexagon compass code.

Classes
-------

`HHC(distance)`
:   Create a new heavy hexagon compass code (HHC).
    
    The X and Z gauge operator lists are given as lists of supports.
    There is a consistent qubit ordering of these lists so that
    other modules can construct gate schedules for circuits.
    
    Initialize data for HHC.

    ### Ancestors (in MRO)

    * css_code.CSSCode

    ### Methods

    `to_index(self, row, col, d)`
    :   Map a coordinate (row, col) to a qubit index.
        
        Qubits are indexed from left to right across each row,
        beginning in the top row.