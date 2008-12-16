# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.circuit.circuit import SubCircuit, VCCS, G, C, IS, VS, Parameter, gnd, R, Nullor, VCVS
from pycircuit.circuit.mos import MOS
from pycircuit.circuit.symbolicanalysis import SymbolicTwoPortAnalysis, SymbolicAC, SymbolicNoise
from sympy import Symbol, simplify, ratsimp, sympify, factor, limit, solve
from numpy import array, zeros
from copy import copy
from pycircuit.circuit.symbolicapprox import *

## Multi FeedBack (MFB) Filter

circuitMFB = SubCircuit()

nin, nout, n1, n2 = circuitMFB.add_nodes('in', 'out', 'n1', 'n2')

R1, R2, R3, C1, C2, Gain = [Symbol(symname, real=True) for symname in 'R1,R2,R3,C1,C2,Gain'.split(',')]

circuitMFB['R3'] = R(nin, n1, r = R3)
circuitMFB['R2'] = R(n1, n2, r = R2)
circuitMFB['R1'] = R(n1, nout, r = R1)
circuitMFB['C1'] = C(n1, gnd, c = C1)
circuitMFB['C2'] = C(n2, nout, c = C2)
circuitMFB['Nullor'] = Nullor(n2, gnd, nout, gnd)
#circuitMFB['VCVS'] = VCVS(n2, gnd, nout, gnd, g = Gain)


circuitMFB['VSource'] = VS(nin, gnd, vac=1)

res = SymbolicAC(circuitMFB).solve(Symbol('s'), complexfreq=True)

## DC Gain
print simplify(res['out'])

## AC Transfer function
print simplify(res['out']).limit('s',0)

## Poles
#print solve(1/simplify(res['out']),'s')

# # Remove soure to able to do an two port analysis
del circuitMFB['VSource']
# ## Perform twoport analysis with noise
twoportana = SymbolicTwoPortAnalysis(circuitMFB, nin, gnd, nout, gnd, noise=True, noise_outquantity='i')

res2port = twoportana.solve(Symbol('s'), complexfreq=True)

y11 = res2port['twoport'].Y[0,0]

# print 'Input impedance:', 1/y11
# #print 'Approx. input impedance', approx(1/y11, ['gds'], n = 1)
# print 'Input referred current noise PSD, Sin:', ratsimp(res2port['Sin'])
# print 'Approx. input referred current noise PSD, Sin:', approx(res2port['Sin'], ['C1'], n=1)
# print 'Input referred voltage noise PSD, Svn:', ratsimp(res2port['Svn'])
