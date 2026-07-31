# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.


import pycircuit.circuit.func as func
from pycircuit.circuit import symbolic, numeric
import sympy
import numpy as np

from numpy.testing import assert_array_equal, assert_almost_equal

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
    assert sin.f(t) == v

    ## NEXT EVENT.  A sine has no discontinuity after `td`, so it contributes no
    ## breakpoints there.  This test used to require an event every quarter period
    ## -- at the peaks and zero crossings -- i.e. it encoded the defect stage 4g(a)
    ## removes.  A breakpoint makes the transient truncate the step to land on it
    ## and then re-arm an order drop across it; for the trapezoidal rule each such
    ## restart re-seeds the alternating (-1)^n mode in the companion current.
    ##
    ## Measured effect of removing them, on a VSin-driven RC over five periods:
    ## evaluations forced to order 1 fell from 120 of 1236 to 3 of 1596, and the
    ## three that remain are evaluations 1, 2 and 3 -- the genuine opening of the
    ## run, which has no history and therefore cannot be higher order.
    ## See doc/transient_work_plan.md stage 4g.
    numeric_sin = func.Sin(offset=0.0, amplitude=1.0, freq=1e6, td=1e-6,
                           theta=0.0, phase=0)

    ## `td` IS a real event: the waveform is only smooth *after* the source starts,
    ## and the derivative is generally discontinuous there.  So it is reported once.
    assert numeric_sin.next_event(0.0) == 1e-6
    assert numeric_sin.next_event(0.5e-6) == 1e-6

    ## After the source has started, nothing -- at a peak, at a zero crossing, and
    ## in between alike.
    assert numeric_sin.next_event(1e-6) == float('inf')
    assert numeric_sin.next_event(1.25e-6) == float('inf')   # peak
    assert numeric_sin.next_event(1.5e-6) == float('inf')    # zero crossing
    assert numeric_sin.next_event(5e-6) == float('inf')

    ## A source with no delay has no events at all.
    undelayed = func.Sin(offset=0.0, amplitude=1.0, freq=1e6, td=0.0,
                         theta=0.0, phase=0)
    assert undelayed.next_event(0.0) == float('inf')

    ## Under a symbolic toolkit `t < td` is a relational rather than a bool, so
    ## asking for its truth value raises.  There is no meaningful breakpoint
    ## schedule for a symbolic time; it must degrade to "no events", not blow up.
    assert sin.next_event(t) == symbolic.inf

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

    ## `next_event(0)` IS `td`, NOT 0.  This line asserted 0 and so pinned the
    ## defect in place: a caller enumerating events by feeding the result back in
    ## never advanced, which hung `jaxtransient.solve_batched` and left
    ## `transient.py:762` re-calling at a nudged `t` to work around it.  Every
    ## other assertion in this block already reads "strictly the next edge after
    ## t" -- this one was the odd one out.  Stage 9(d).
    assert_almost_equal(pulse.next_event(0), td)
    assert_almost_equal(pulse.next_event(td/2), td)
    assert_almost_equal(pulse.next_event(td), td+tr)
    assert_almost_equal(pulse.next_event(td+tr/2), td+tr)
    assert_almost_equal(pulse.next_event(td+tr+pw), td+tr+pw+tf)
    assert_almost_equal(pulse.next_event(td+tr+pw-eps), td+tr+pw)
    assert_almost_equal(pulse.next_event(td+tr+pw+tf), per)
    assert_almost_equal(pulse.next_event(td+tr+pw+tf-eps), td+tr+pw+tf)
    assert_almost_equal(pulse.next_event(per+td/2), per+td)

    ## The invariant itself, not just the table above: every TimeFunction must
    ## return strictly more than its argument, or an enumerating caller stalls.
    ## Walked rather than sampled, because the failure showed up only after the
    ## accumulated value drifted off the edge grid -- see Pulse.next_event.
    t_walk = 0.0
    for _ in range(40):
        nxt = float(pulse.next_event(t_walk))
        assert nxt > t_walk, \
            'next_event(%r) returned %r, which does not advance' % (t_walk, nxt)
        t_walk = nxt
