# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.circuit.circuit 
from pycircuit.circuit import *
from pycircuit.circuit import symbolic
from pycircuit.circuit import symbolic_poly
import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from copy import copy
from pycircuit.circuit.tests.test_circuit import create_current_divider
from .test_circuit import create_current_divider
from sympy import var, simplify, integrate, oo, limit, gruntz, pi, I
import unittest
from numpy.testing import assert_equal

## Run every symbolic test on both the stock and the polynomial-domain toolkit
TOOLKITS = [symbolic, symbolic_poly]

@pytest.fixture(autouse=True)
def reset_default_toolkit():
    """Reset the default toolkit after each test to prevent leaking symbolic."""
    from pycircuit.circuit.toolkit import numeric
    yield
    pycircuit.circuit.circuit.default_toolkit = numeric

def check_symbolic_ac(toolkit):
    pycircuit.circuit.circuit.default_toolkit = toolkit
    cir = SubCircuit()

    var('v0 R1 C1 s')

    cir['R1'] = R(1, 2, r=R1)
    cir['R2'] = C(2, gnd, c=C1)
    cir['VS'] = VS(1, gnd, vac=v0)

    res = AC(cir, toolkit = toolkit).solve(freqs = s, complexfreq=True)
    assert simplify(res.v(2,gnd)-v0/(1+s*R1*C1)) == 0

def check_symbolic_noise_vin_vout(toolkit):
    pycircuit.circuit.circuit.default_toolkit = toolkit
    c = SubCircuit()

    var('R1 R2 V', real=True, positive=True)
    var('s')

    c['vs'] = VS(1, gnd, vac=V)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)

    noise = Noise(c, inputsrc='vs', outputnodes=('2', gnd),
                  toolkit= toolkit)
    res = noise.solve(s, complexfreq=True)

    assert simplify(res['Svnout']) == simplify(4*R1*R2*noise.toolkit.kboltzmann*noise.par.epar.T/(R1 + R2))
    assert simplify(res['Svninp']) == simplify(4*noise.toolkit.kboltzmann*noise.par.epar.T*R1*(R2 + R1)/R2)
    assert simplify(res['gain'] - R2 / (R1 + R2)) == 0

def check_symbolic_noise_vin_iout(toolkit):
    pycircuit.circuit.circuit.default_toolkit = toolkit
    c = SubCircuit()
    
    var('R1 R2 R3 V', real=True, positive=True)
    var('s')

    c['vs'] = VS(1, gnd, vac=V)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)
    c['vl'] = VS(2, gnd)

    noise = Noise(c, inputsrc='vs', outputsrc='vl', 
                  toolkit= toolkit)
    res = noise.solve(s, complexfreq=True)
    
    assert simplify(res['Sinout']) == simplify(4*noise.toolkit.kboltzmann*noise.par.epar.T*(R1+R2)/(R1*R2))
    assert simplify(res['Svninp']) == simplify(4*noise.toolkit.kboltzmann*noise.par.epar.T*R1*(R2+R1)/R2)
    assert simplify(res['gain']) == 1/R1

def check_symbolic_noise_iin_vout(toolkit):
    pycircuit.circuit.circuit.default_toolkit = toolkit
    c = SubCircuit()
    
    var('R1 R2', real=True)
    var('Iin s')

    c['is'] = IS(1, gnd, iac=Iin)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)

    noise = Noise(c, inputsrc='is', outputnodes=('2', gnd), 
                  toolkit= toolkit)
    res = noise.solve(s, complexfreq=True)

    assert simplify(res['Svnout']) == 4*R2*noise.toolkit.kboltzmann*noise.par.epar.T
    assert simplify(res['Sininp']) == 4*noise.toolkit.kboltzmann*noise.par.epar.T/R2
    assert simplify(res['gain']) == R2


def check_symbolic_noise_iin_iout(toolkit):
    pycircuit.circuit.circuit.default_toolkit = toolkit
    c = SubCircuit()
    
    var('R1 R2 R3', real=True)
    var('Iin s')
    k = toolkit.kboltzmann

    c['is'] = IS(1, gnd, iac=Iin)
    c['R1'] = R(1, 2, r=R1)
    c['R2'] = R(2, gnd, r=R2)
    c['vl'] = VS(2, gnd)

    noise = Noise(c, inputsrc='is', outputsrc='vl', 
                  toolkit= toolkit)
    res = noise.solve(s, complexfreq=True)

    T = noise.par.epar.T

    assert simplify(res['Sinout']) == 4*k*T/R2
    assert simplify(res['Sininp']) == 4*k*T/R2
    assert simplify(res['gain']) == 1

def check_symbolic_noise_kt_over_C(toolkit):
    pycircuit.circuit.circuit.default_toolkit = toolkit
    cir = SubCircuit(toolkit = toolkit)

    k = toolkit.kboltzmann
    var('r c w w1 V', real=True, positive=True)

    s = I * w

    cir['vs'] = VS(1, gnd, vac=V)
    cir['R'] = R(1, 2, r=r)
    cir['C'] = C(2, gnd, c=c)

    noise = Noise(cir, inputsrc='vs', outputnodes=('2', gnd), 
                  toolkit = toolkit)
    res = noise.solve(s, complexfreq=True)

    T = noise.par.epar.T

    svnout = simplify(res['Svnout'])

    noise_voltage_power = simplify(integrate(svnout, (w, 0, oo)))

    assert noise_voltage_power == 2*pi*k*T/c

    assert simplify(res['gain'] - 1/(1 + s*r*c)) == 0

import pytest

@pytest.mark.parametrize('t', TOOLKITS)
def test_symbolic_ac(t):
    check_symbolic_ac(t)

@pytest.mark.parametrize('t', TOOLKITS)
def test_symbolic_noise_vin_vout(t):
    check_symbolic_noise_vin_vout(t)

@pytest.mark.parametrize('t', TOOLKITS)
def test_symbolic_noise_vin_iout(t):
    check_symbolic_noise_vin_iout(t)

@pytest.mark.parametrize('t', TOOLKITS)
def test_symbolic_noise_iin_vout(t):
    check_symbolic_noise_iin_vout(t)

@pytest.mark.parametrize('t', TOOLKITS)
def test_symbolic_noise_iin_iout(t):
    check_symbolic_noise_iin_iout(t)

@pytest.mark.parametrize('t', TOOLKITS)
def test_symbolic_noise_kt_over_C(t):
    check_symbolic_noise_kt_over_C(t)
