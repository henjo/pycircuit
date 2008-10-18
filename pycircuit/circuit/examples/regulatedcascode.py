from pycircuit.circuit.circuit import SubCircuit, VCCS, G, C, IS, VS, Parameter, gnd, R
from pycircuit.circuit.mos import MOSSmallSignal
from pycircuit.circuit.symbolicanalysis import SymbolicTwoPortAnalysis, SymbolicAC, SymbolicNoise
from sympy import Symbol, simplify, ratsimp, sympify
from numpy import array, zeros
from copy import copy

c = SubCircuit()

nin, nout, n1 = [c.addNode(name) for name in ('in', 'out', 'n1')]

gm1, gm2, gds1, gds2 = [Symbol(s, real=True) for s in 'gm1,gm2,gds1,gds2'.split(',')]

c['M2'] = MOSSmallSignal(nin, n1, gnd, gnd, gm = gm2, gds = gds2)
c['M1'] = MOSSmallSignal(n1, nout, nin, nin, gm = gm1, gds = gds1)
c['r'] = R(nin, gnd, r = Symbol('Rs', real=True))

res2port = SymbolicTwoPortAnalysis(c, nin, gnd, nout, gnd).run(freqs = array([Symbol('s')]), complexfreq=True)

y11 = res2port['twoport'].Y[0,0]

print 'Input impedance:', (1/y11).y[0]

## Calculate noise

## current noise
c['vl'] = VS(nout, gnd, vac = 0)
c2 = copy(c)
c['is'] = IS(nin, gnd, iac = Symbol('iin'))

resnoise = SymbolicNoise(c, inputsrc=c['is'], outputsrc=c['vl']).run()
print 'Gain: ', resnoise['gain']
print 'Input referred current noise: ', resnoise['in2in']

## voltage noise
c = c2
c['vs'] = VS(nin, gnd, vac = Symbol('vin'))

resnoise = SymbolicNoise(c, inputsrc=c['vs'], outputsrc=c['vl']).run()

print 'Input referred voltage noise: ', resnoise['vn2in']
