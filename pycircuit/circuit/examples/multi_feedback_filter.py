# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from sympy import *
from pycircuit.circuit import *
from pycircuit.circuit.toolkit import symbolic

## Multi FeedBack (MFB) Filter
print("Multi FeedBack (MFB) Filter Example")

circuit_MFB = SubCircuit(toolkit=symbolic)

nin, nout, n1, n2, n3 = circuit_MFB.add_nodes('in', 'out', 'n1', 'n2','n3')

R1, R2, R3, R4, C1, C2 = [Symbol(symname, real=True, positive=True)
                          for symname in 'R1,R2,R3,R4,C1,C2'.split(',')]

circuit_MFB['R3'] = R(nin, n1, r = R3)
circuit_MFB['R2'] = R(n1, n2, r = R2)
circuit_MFB['R1'] = R(n1, n3, r = R1)
circuit_MFB['C1'] = C(n1, gnd, c = C1)
circuit_MFB['C2'] = C(n2, n3, c = C2)
circuit_MFB['R4'] = R(n3, nout, r = R4) # Added to output to test error in twoport analysis
circuit_MFB['Nullor'] = Nullor(n2, gnd, n3, gnd)

# Voltage source for AC analysis
circuit_MFB['VSource'] = VS(nin, gnd, vac=1)

s = Symbol('s', complex = True)
res = AC(circuit_MFB, toolkit=symbolic).solve(s, complexfreq=True)
## res.v(node, refnode) is the modern accessor -- the 2008-era dict-style
## res['out'] never worked against AC's result object (empty by design; see
## doc/architecture.md P14).
res_out = simplify(res.v(nout, gnd))

## DC Gain
dc_gain = simplify(res_out).limit(s, 0)
print("")
print("DC Gain:")
pprint(dc_gain)
print("")

## AC Transfer function
tf = collect(res_out, s)
print("AC Transfer function:")
pprint(tf)
print("")

## Denominator of transfer function
tf_denom = fraction(tf)[1]
tf_denom = tf_denom.expand()
print("Denominator of transfer function:")
pprint(tf_denom)
print("")

## Poles
tf_poles = solve(tf_denom, s)
print("Poles of transfer function:")
pprint(tf_poles)
print("")

## Remove source to be able to do a twoport analysis
del circuit_MFB['VSource']
## Perform twoport analysis with noise
twoportana = TwoPortAnalysis(circuit_MFB, nin, gnd, nout, gnd, noise=True,
                             noise_outquantity='i', toolkit=symbolic)

res2port = twoportana.solve(Symbol('s'), complexfreq=True)

y11 = res2port['twoport'].Y[0,0]

## The denominator here is the same characteristic polynomial as the
## transfer function above -- input impedance and the transfer function of
## the same circuit share poles, a useful cross-check that this is actually
## right, not merely non-crashing.
print('Input impedance:')
pprint(simplify(1/y11))
print("")

## Sin/Svn (input-referred noise PSDs) are also available from res2port but
## are omitted here: ratsimp on them does not return quickly for this
## circuit's noise expressions, which would make this example slow to run.
