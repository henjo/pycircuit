# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""
Physical constants
"""

## ⚠ SI 2019 exact value.  This was 1.38e-23 until 2026-09-05 -- 4.7e-4 low,
## which the Spectre comparison suite had to carry as a PARAMETER on both
## sides of every noise test to keep the tools from disagreeing over a
## constant.  `qelectron` below is the same class of number (1.602e-19
## against the exact 1.602176634e-19, 1.1e-4) and is NOT changed here.
kboltzmann=1.380649e-23   # Boltzmann's constant, exact (SI 2019)
eps0 = 8.8542e-12         # Vacuum permittivity
epsRSi = 11.7             # Relative permittivity of Si
epsRSiO2 = 3.9            # Relative permittivity of SiO2 
qelectron=1.602e-19       # Elementary charge
