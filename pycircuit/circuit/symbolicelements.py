# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import circuit
from circuit import *
from constants_sympy import *
from sympy import exp

class Diode(circuit.Diode, Circuit):
    def i(self, x, epar=circuit.defaultepar):
        VD = x[0]-x[1]
        VT = kboltzmann*epar.T / qelectron
        I = self.mpar.IS*(exp(VD/VT)-1.0)
        return self.toolkit.array([I, -I])

if __name__ == "__main__":
    import doctest
    doctest.testmod()
