# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from nose.tools import *
import pycircuit.circuit.circuit 
from pycircuit.circuit import *
from pycircuit.circuit import symbolic
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from copy import copy
from test_circuit import create_current_divider
from sympy import var, simplify, integrate, oo, limit, gruntz, pi

def test_symbolic_ac():
    pycircuit.circuit.circuit.default_toolkit = symbolic
    cir=SubCircuit()

    var('v0 R1 C1 s')

    cir['R1']=R(1, 2, r=R1)
    cir['R2']=C(2, gnd, c=C1)
    cir['VS']=VS(1, gnd, vac=v0)

    res = AC(cir, toolkit = symbolic).solve(freqs = s, complexfreq=True)
    assert_equal(simplify(res.v(2,gnd)-v0/(1+s*R1*C1)), 0)

def test_symbolic_noise_vin_vout():
    pycircuit.circuit.circuit.default_toolkit = symbolic
    c = SubCircuit()

    var('R1 R2 k T V', real=True, positive=True)

    c['vs'] = VS(1, gnd, vac=V)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)

    noise = Noise(c, inputsrc='vs', outputnodes=('2', gnd), 
                  toolkit=symbolic)
    res = noise.solve(s, complexfreq=True)

    assert_equal(simplify(res['Svnout']), simplify(4*R1*R2*k*T/(R1 + R2)))
    assert_equal(simplify(res['Svninp']), (4*R1*R2*k*T + 4*k*T*R1**2)/R2)
    assert_equal(simplify(res['gain'] - R2 / (R1 + R2)), 0)

def test_symbolic_noise_vin_iout():
    pycircuit.circuit.circuit.default_toolkit = symbolic
    c = SubCircuit()
    
    var('R1 R2 R3 k T V', real=True, positive=True)

    c['vs'] = VS(1, gnd, vac=V)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)
    c['vl'] = VS(2, gnd)

    noise = Noise(c, inputsrc='vs', outputsrc='vl', 
                  toolkit=symbolic)
    res = noise.solve(s, complexfreq=True)
    
    assert_equal(simplify(res['Sinout']), simplify(4*k*T*(R1+R2)/(R1*R2)))
    assert_equal(simplify(res['Svninp']), (4*R1*R2*k*T + 4*k*T*R1**2)/R2)
    assert_equal(simplify(res['gain']), 1/R1)

def test_symbolic_noise_iin_vout():
    pycircuit.circuit.circuit.default_toolkit = symbolic
    c = SubCircuit()
    
    var('R1 R2', real=True)
    var('k T I')

    c['is'] = IS(1, gnd, iac=I)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)

    noise = Noise(c, inputsrc='is', outputnodes=('2', gnd), 
                  toolkit=symbolic)
    res = noise.solve(s, complexfreq=True)

    assert_equal(simplify(res['Svnout']), 4*R2*k*T)
    assert_equal(simplify(res['Sininp']), 4*k*T/R2)
    assert_equal(simplify(res['gain']), R2)


def test_symbolic_noise_iin_iout():
    pycircuit.circuit.circuit.default_toolkit = symbolic
    c = SubCircuit()
    
    var('R1 R2 R3', real=True)
    var('k T I s')

    c['is'] = IS(1, gnd, iac=I)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)
    c['vl'] = VS(2, gnd)

    noise = Noise(c, inputsrc='is', outputsrc='vl', 
                  toolkit=symbolic)
    res = noise.solve(s, complexfreq=True)

    assert_equal(simplify(res['Sinout']), 4*k*T/R2)
    assert_equal(simplify(res['Sininp']), 4*k*T/R2)
    assert_equal(simplify(res['gain']), 1)

def test_symbolic_noise_kt_over_C():
    pycircuit.circuit.circuit.default_toolkit = symbolic
    cir = SubCircuit(toolkit = symbolic)

    var('r c w w1 k T V', real=True, positive=True)

    s = 1j * w

    cir['vs'] = VS(1, gnd, vac=V)
    cir['R'] = R(1, 2, r=r)
    cir['C'] = C(2, gnd, c=c)

    noise = Noise(cir, inputsrc='vs', outputnodes=('2', gnd), 
                  toolkit = symbolic)
    res = noise.solve(s, complexfreq=True)


    noise_voltage_power = limit(integrate(simplify(res['Svnout']), 
                                          (w, 0, w1)), w1, oo)

    assert_equal(noise_voltage_power, 2*pi*k*T/c)

    assert_equal(simplify(res['gain'] - 1/(1 + s*r*c)), 0)

def test_integer_component_values():
    """Test dc analysis with integer component values
    
       As python per default uses integer arithmetics integer
       component values can lead to problems. (Python < 3.0)
    """
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit(toolkit=numeric)

    n1,n2 = c.add_nodes('net1', 'net2')

    c['vs'] = VS(n1, gnd, v = 9)
    c['R1'] = R( n1,  n2, r = 50)
    c['R2'] = R( n2, gnd, r = 50)

    dc = DC(c)
    res = dc.solve()
    
    assert_equal(res.v('net1'), 9.0)
