Module css_code
===============
Module to define a class for CSS codes.

Classes
-------

`CSSCode(x_gauges, z_gauges, n, k, d, x_stabilizers, z_stabilizers, logical_z, logical_x)`
:   Class defining a CSS code.
    
    Initialize data for CSS code.

    ### Methods

    `logical_x_error(self, bitstring)`
    :   Test for a logical X error.

    `logical_z_error(self, bitstring)`
    :   Test for a logical Z error.

    `x_syndrome(self, bitstring)`
    :   Compute the X syndrome of a length self.n bit string.

    `z_syndrome(self, bitstring)`
    :   Compute the Z syndrome of a length self.n bit string.