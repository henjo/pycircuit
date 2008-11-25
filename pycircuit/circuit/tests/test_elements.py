"""Circuit element tests
"""

from pycircuit.circuit.circuit import VS, R, Nullor, SubCircuit, gnd
from pycircuit.circuit.symbolic import SymbolicAC
from sympy import Symbol, Matrix, symbols, simplify, together, factor, cancel

def test_nullor_vva():
    """Test nullor element by building a V-V amplifier"""
    
    c = SubCircuit()

    Vin = Symbol('Vin')
    R1 =Symbol('R1')
    R2 = Symbol('R2')
    
    nin = c.add_node('in')
    n1 = c.add_node('n1')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, v=Vin)
    c['R1'] = R(n1, gnd, r=R1)
    c['R2'] = R(nout, n1, r=R2)
    c['nullor'] = Nullor(n1, nin, gnd, nout)
    
    result = SymbolicAC(c).run()
    
    vout = result['out']
    
    assert simplify(vout - Vin * (R1 + R2) / R1) == 0, \
        'Did not get the expected result, %s != 0'% \
        str(simplify(vout - Vin * (R1 + R2) / R1))

if __name__ == '__main__':
    test_nullor_vva()
