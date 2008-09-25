"""Circuit element tests
"""

from numpy import *
from scipy import *
from pycircuit.circuit.circuit import VS, R, Nullor, SubCircuit, gnd
from pycircuit.circuit.symbolic import SymbolicAC
from sympy import Symbol, Matrix, symbols, simplify, together, factor, cancel


def test_nullor_vva():
    """Test nullor element by building a V-V amplifier"""
    
    c = SubCircuit()

    Vin = Symbol('Vin')
    R1 =Symbol('R1')
    R2 = Symbol('R2')
    
    nin = c.addNode('in')
    n1 = c.addNode('n1')
    n2 = c.addNode('n2')
    nout = c.addNode('out')
    
    c['vin'] = VS(nin, gnd, v=Vin)
    c['R1'] = R(n1, gnd, r=R1)
    c['R2'] = R(n2, n1, r=R2)
    c['nullor'] = Nullor(nin, n1, nout, n2)
    
    result = SymbolicAC(c).solve()

    vout = result.getSignal('out')

    assert vout - Vin * R1 / (R1 + R2) == 0

if __name__ == '__main__':
    test_nullor_vva()
