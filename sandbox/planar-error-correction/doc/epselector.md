Module epselector
=================
Factory for error propagators.

Classes
-------

`EPSelector()`
:   Select and return an ErrorPropagator.
    
    Create a new factory.

    ### Methods

    `get_error_propagator(self, eptype='auto')`
    :   Return an error propagator.
        
        auto = try to automatically detect extensions
        c = force the compiled version
        py = force the python version