# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.circuit import SubCircuit, VCCS, G, C, IS, VS, Parameter, gnd, \
    R, symbolic, TwoPortAnalysis, AC, Noise
from pycircuit.circuit.mos import MOS
from sympy import Symbol, simplify, ratsimp, sympify, factor, limit
from numpy import array, zeros
from copy import copy
from pycircuit.circuit.symbolicapprox import *
    
## Create regulated cascode circuit

c = SubCircuit()

nin, nout, n1 = c.add_nodes('in', 'out', 'n1')

gm1, gm2, gds1, gds2, Cgs1, Cgs2= [Symbol(symname, real=True) for symname in 'gm1,gm2,gds1,gds2,Cgs1,Cgs2'.split(',')]

c['M2'] = MOS(nin, n1, gnd, gnd, gm = gm2, gds = gds2, Cgs=0*Cgs1)
c['M1'] = MOS(n1, nout, nin, nin, gm = gm1, gds = gds1, Cgs=0*Cgs2)
#c['r'] = R(nin, gnd, r = Symbol('Rs', real=True))

## Perform twoport analysis with noise

twoportana = TwoPortAnalysis(c, nin, gnd, nout, gnd, noise=True, noise_outquantity='i', toolkit=symbolic)

res2port = twoportana.solve(Symbol('s'), complexfreq=True)

y11 = res2port['twoport'].Y[0,0]

print('Input impedance:', 1/y11)
#print('Approx. input impedance', approx(1/y11, ['gds'], n = 1))
print('Input referred current noise PSD, Sin:', ratsimp(res2port['Sin']))
print('Approx. input referred current noise PSD, Sin:', approx(res2port['Sin'], ['gds'], n=1))
print('Input referred voltage noise PSD, Svn:', ratsimp(res2port['Svn']))
