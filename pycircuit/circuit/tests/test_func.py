# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from nose.tools import *

import pycircuit.circuit.func as func
from pycircuit.circuit import symbolic, numeric
import sympy
import numpy as np

from numpy.testing import assert_array_equal

def test_timefunction():
    for toolkit in symbolic, numeric:
        f = func.TimeFunction(toolkit=toolkit)
        assert f.f(0) == 0
        assert f.next_event(0) == toolkit.inf
    
def test_sin():
    sympy.var('vo va freq td theta t')
    phase = sympy.Symbol('phase')

    sin = func.Sin(toolkit = symbolic,
                   offset=vo, amplitude=va, freq=freq, td=td, 
                   theta=theta, phase=phase)

    v = vo + va*sympy.exp(-theta*(t - td)) * \
        sympy.sin(2*sympy.pi*freq*(t-td)+phase*sympy.pi/180)
    assert_equal(sin.f(t), v)

    ## Test next event, phase = 0
    sin = func.Sin(toolkit = symbolic,
                   offset=vo, amplitude=va, freq=freq, td=td, 
                   theta=theta, phase=0)
    period = 1/freq
    assert_equal(sin.next_event(period+td), period+td + period/4)
    assert_equal(sin.next_event(period+td + period / 8), period+td + period/4)
    assert_equal(sin.next_event(period+td - period / 16), period+td)

    ## Test next event, phase = phase
    phase = 1
    sin = func.Sin(toolkit = symbolic,
                   offset=vo, amplitude=va, freq=freq, td=td, 
                   theta=theta, phase=phase)
    period = 1/freq
    t_nextevent = sin.next_event(period+td + period / 8 - phase*period/360)
    assert_equal(t_nextevent.expand(),
                 (period+td + period/4 - phase*period/360).expand())
    
def test_pulse():
    t = sympy.Symbol('t')

    v1 = 1.1
    v2 = -0.9

    td = 0.4
    tr = 0.1
    tf = 0.1
    pw = 0.5
    per = 2.0
    
    eps = 1e-6
    
    pulse = func.Pulse(toolkit = symbolic,
                       v1=v1, v2=v2, td=td, tr=tr, tf=tf, pw=pw, per=per)
    
    tpoints = np.array((0,td,td+tr,td+tr+pw,td+tr+pw+tf,10))
    vpoints = np.array((v1,v1,v2,v2,v1,v1))
    
    tref = np.arange(0,per, 0.005)
    
    for tstart in 0,per:
        for t in tref:
            vref = np.interp(t,tpoints,vpoints)
            assert_almost_equal(pulse.f(t + tstart), vref)

    assert_almost_equal(pulse.next_event(0), 0)
    assert_almost_equal(pulse.next_event(td/2), td)
    assert_almost_equal(pulse.next_event(td), td+tr)
    assert_almost_equal(pulse.next_event(td+tr/2), td+tr)
    assert_almost_equal(pulse.next_event(td+tr+pw), td+tr+pw+tf)
    assert_almost_equal(pulse.next_event(td+tr+pw-eps), td+tr+pw)
    assert_almost_equal(pulse.next_event(td+tr+pw+tf), per)
    assert_almost_equal(pulse.next_event(td+tr+pw+tf-eps), td+tr+pw+tf)
    assert_almost_equal(pulse.next_event(per+td/2), per+td)
