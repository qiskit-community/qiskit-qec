import numpy as np
import copy
#from properties import Properties

class SymplecticMatrixBase:
    DENSELIST = 0
    DENSENUMPY  = 1
    SPARSESCIPY = 2
    DENSEBIT = 3
    SPARSEBIT = 4

    STYPES = [
        DENSELIST,
        DENSENUMPY,
        SPARSESCIPY,
        DENSEBIT,
        SPARSEBIT]

    CLASSNAMES = {
        DENSELIST:"SymplecticMatrixDenseList",
        DENSENUMPY: "SymplecticMatrixDenseNumPy",
        SPARSESCIPY:"SymplecticMatrixSparseSciPy",
        DENSEBIT:"SymplecticMatrixDenseBit",
        SPARSEBIT:"SymplecticMatrixSparseBit"
    }

    def __init__(self) -> None:
        #self.properties = Properties(__class__.__name__)
        pass


