from pycircuit.circuit.circuit import SubCircuit, VCCS, G, C, IS, VS, Parameter, gnd, R
from pycircuit.circuit.mos import MOS
from pycircuit.circuit.symbolicanalysis import SymbolicTwoPortAnalysis, SymbolicAC, SymbolicNoise
from sympy import Symbol, simplify, ratsimp, sympify, factor
from numpy import array, zeros
from copy import copy

## Create regulated cascode circuit

c = SubCircuit()

nin, nout, n1 = c.addNodes('in', 'out', 'n1')

gm1, gm2, gds1, gds2 = [Symbol(symname, real=True) for symname in 'gm1,gm2,gds1,gds2'.split(',')]

c['M2'] = MOS(nin, n1, gnd, gnd, gm = gm2, gds = gds2)
c['M1'] = MOS(n1, nout, nin, nin, gm = gm1, gds = gds1)
#c['r'] = R(nin, gnd, r = Symbol('Rs', real=True))

## Perform twoport analysis with noise

twoportana = SymbolicTwoPortAnalysis(c, nin, gnd, nout, gnd, noise=True, noise_outquantity='i')

res2port = twoportana.run(Symbol('s'), complexfreq=True)

y11 = res2port['twoport'].Y[0,0]

## Approximate by replacing gdsX by gdsX*t and do a taylor series around t=0
t = Symbol('t')
zint = 1/y11.subs({'gds1': gds1*t, 'gds2': gds2*t})
zinapprox = zint.series('t', n=2).subs({'t': 1})

print 'Input impedance:', 1/y11
print 'Approx. input impedance', zinapprox
print 'Input referred current noise PSD, Sin:', ratsimp(res2port['Sin'])
print 'Input referred voltage noise PSD, Svn:', ratsimp(res2port['Svn'])

