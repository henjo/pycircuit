# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from nose.tools import *
from pycircuit.circuit import *
from pycircuit.circuit import symbolic
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from copy import copy
from test_circuit import create_current_divider
from sympy import var, simplify, integrate, oo, limit, gruntz, pi

def test_symbolic_ac():
    cir=SubCircuit()

    var('v0 R1 C1 s')

    cir['R1']=R(1, 2, r=R1)
    cir['R2']=C(2, gnd, c=C1)
    cir['VS']=VS(1, gnd, vac=v0)

    res = AC(cir, toolkit = symbolic).solve(freqs = s, complexfreq=True)
    assert_equal(simplify(res.v(2,gnd)-v0/(1+s*R1*C1)), 0)

def test_symbolic_noise_vin_vout():
    c = SubCircuit()

    var('R1 R2 kT V', real=True, positive=True)

    c['vs'] = VS(1, gnd, vac=V)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)

    noise = Noise(c, inputsrc='vs', outputnodes=('2', gnd), 
                  toolkit=symbolic)
    res = noise.solve(s, complexfreq=True)

    assert_equal(simplify(res['Svnout']), simplify(4*R1*R2*kT/(R1 + R2)))
    assert_equal(simplify(res['Svninp']), (4*R1*R2*kT + 4*kT*R1**2)/R2)
    assert_equal(simplify(res['gain'] - R2 / (R1 + R2)), 0)

def test_symbolic_noise_vin_iout():
    c = SubCircuit()
    
    var('R1 R2 R3 kT V', real=True, positive=True)

    c['vs'] = VS(1, gnd, vac=V)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)
    c['vl'] = VS(2, gnd)

    noise = Noise(c, inputsrc='vs', outputsrc='vl', 
                  toolkit=symbolic)
    res = noise.solve(s, complexfreq=True)
    
    assert_equal(simplify(res['Sinout']), simplify(4*kT*(R1+R2)/(R1*R2)))
    assert_equal(simplify(res['Svninp']), (4*R1*R2*kT + 4*kT*R1**2)/R2)
    assert_equal(simplify(res['gain']), 1/R1)

def test_symbolic_noise_iin_vout():
    c = SubCircuit()
    
    var('R1 R2', real=True)
    var('kT I')

    c['is'] = IS(1, gnd, iac=I)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)

    noise = Noise(c, inputsrc='is', outputnodes=('2', gnd), 
                  toolkit=symbolic)
    res = noise.solve(s, complexfreq=True)

    assert_equal(simplify(res['Svnout']), 4*R2*kT)
    assert_equal(simplify(res['Sininp']), 4*kT/R2)
    assert_equal(simplify(res['gain']), R2)


def test_symbolic_noise_iin_iout():
    c = SubCircuit()
    
    var('R1 R2 R3', real=True)
    var('kT I')

    c['is'] = IS(1, gnd, iac=I)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)
    c['vl'] = VS(2, gnd)

    noise = Noise(c, inputsrc='is', outputsrc='vl', 
                  toolkit=symbolic)
    res = noise.solve(s, complexfreq=True)

    assert_equal(simplify(res['Sinout']), 4*kT/R2)
    assert_equal(simplify(res['Sininp']), 4*kT/R2)
    assert_equal(simplify(res['gain']), 1)

def test_symbolic_noise_kt_over_C():
    cir = SubCircuit()

    var('r c w w1 kT V', real=True, positive=True)

    s = 1j * w

    cir['vs'] = VS(1, gnd, vac=V)
    cir['R'] = R(1, 2, r=r)
    cir['C'] = C(2, gnd, c=c)

    noise = Noise(cir, inputsrc='vs', outputnodes=('2', gnd), 
                  toolkit=symbolic)
    res = noise.solve(s, complexfreq=True)


    noise_voltage_power = gruntz(integrate(simplify(res['Svnout']), 
                                           (w, 0, w1)), w1, oo)

    assert_equal(noise_voltage_power, 2*pi*kT/c)

    assert_equal(simplify(res['gain'] - 1/(1 + s*r*c)), 0)
