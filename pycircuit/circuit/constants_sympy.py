# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from sympy import Symbol
"""
Physical constants
"""

## Symbolics for constants
kboltzmann=sympy.Symbol('k', real=True, positive=True)         # Boltzmann's constant
eps0 = sympy.Symbol('eps0', real=True, positive=True)          # Vacuum permittivity
epsRSi = sympy.Symbol('epsRSi', real=True, positive=True)      # Relative permittivity of Si
epsRSiO2 = sympy.Symbol('epsRSiO2', real=True, positive=True)  # Relative permittivity of SiO2 
qelectron = sympy.Symbol('qelectron', positive=True, real=True)# Elementary charge


