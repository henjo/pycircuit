# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.circuit.circuit import SubCircuit, VCCS, G, C, IS, VS, Parameter, gnd, R, Nullor, VCVS
from pycircuit.circuit.mos import MOS
from pycircuit.circuit.symbolicanalysis import SymbolicTwoPortAnalysis, SymbolicAC, SymbolicNoise
from sympy import Symbol, simplify, ratsimp, sympify, factor, limit, solve, pprint, fraction, collect
from numpy import array, zeros
from copy import copy
from pycircuit.circuit.symbolicapprox import *

## Multi FeedBack (MFB) Filter
print("Multi FeedBack (MFB) Filter Example")

circuit_MFB = SubCircuit()

nin, nout, n1, n2, n3 = circuit_MFB.add_nodes('in', 'out', 'n1', 'n2','n3')

R1, R2, R3, R4, C1, C2, Gain = [Symbol(symname, real=True) for symname in 'R1,R2,R3,R4,C1,C2,Gain'.split(',')]

circuit_MFB['R3'] = R(nin, n1, r = R3)
circuit_MFB['R2'] = R(n1, n2, r = R2)
circuit_MFB['R1'] = R(n1, n3, r = R1)
circuit_MFB['C1'] = C(n1, gnd, c = C1)
circuit_MFB['C2'] = C(n2, n3, c = C2)
circuit_MFB['R4'] = R(n3, nout, r = R4) # Added to output to test error in twoport analysis
circuit_MFB['Nullor'] = Nullor(n2, gnd, n3, gnd)

# Voltage source for AC analysis
circuit_MFB['VSource'] = VS(nin, gnd, vac=1)

res = SymbolicAC(circuit_MFB).solve(Symbol('s'), complexfreq=True) # What assumptions are added to symbol 's'?

res_out = simplify(res['out'])

s = Symbol('s', complex = True)

res_simp = simplify(res_out.subs('s',s))

## DC Gain
dc_gain = simplify(res_simp).limit('s',0)
print("")
print("DC Gain:")
pprint(dc_gain)
print("")

## AC Transfer function
tf = collect(res_simp,s)
print("AC Transfer function:")
pprint(tf)
print("")

## Denominator of transfer function
tf_denom = fraction(tf)[1]
tf_denom = tf_denom.expand()
print("Denominator of transfer function times R[1]:")
pprint(tf_denom)
print("")

## Poles
tf_poles = solve(tf_denom,s) 
tf_poles = simplify(tf_poles)
print("Poles of transfer function:")
print(tf_poles)
p1 = tf_poles # Need to exctract one pole
p2 = tf_poles # Need to extract one pole
pprint(p1)
pprint(p2)
print("")

## Remove soure to able to do an two port analysis
del circuit_MFB['VSource']
## Perform twoport analysis with noise
twoportana = SymbolicTwoPortAnalysis(circuit_MFB, nin, gnd, nout, gnd, noise=True, noise_outquantity='i')

res2port = twoportana.solve(Symbol('s'), complexfreq=True)

y11 = res2port['twoport'].Y[0,0]

# print('Input impedance:', 1/y11)
# #print('Approx. input impedance', approx(1/y11, ['gds'], n = 1))
# print('Input referred current noise PSD, Sin:', ratsimp(res2port['Sin']))
# print('Approx. input referred current noise PSD, Sin:', approx(res2port['Sin'], ['C1'], n=1))
# print('Input referred voltage noise PSD, Svn:', ratsimp(res2port['Svn']))
