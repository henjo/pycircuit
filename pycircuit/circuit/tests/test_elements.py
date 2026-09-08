# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Circuit element tests
"""

import pycircuit.circuit.circuit 
from pycircuit.circuit import *
from pycircuit.circuit.elements import *
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from sympy import var, Symbol, simplify, symbols
import sympy
import unittest

from pylab import plot, show
from pycircuit.circuit.transient import Transient

def test_vsin():
    var('vo va freq td theta phase t')
    vsin = VSin(toolkit = symbolic,
                vo=vo, va=va, freq=freq, td=td, theta=theta, phase=phase)

    v = vo + va*sympy.exp(-theta*(t - td)) * \
        sympy.sin(2*sympy.pi*freq*(t-td)+phase*sympy.pi/180)

    assert_array_equal(vsin.u(t, analysis='tran'), symbolic.array([0,0,-v]))

def test_vpulse():
    t = sympy.Symbol('t')

    v1 = 1.1
    v2 = -0.9

    td = 0.4
    tr = 0.1
    tf = 0.1
    pw = 0.5
    per = 2.0
    
    eps = 1e-6
    
    pulse = VPulse(toolkit = symbolic,
                   v1=v1, v2=v2, td=td, tr=tr, tf=tf, pw=pw, per=per)
    
    tpoints = np.array((0,td,td+tr,td+tr+pw,td+tr+pw+tf,10))
    vpoints = np.array((v1,v1,v2,v2,v1,v1))
    
    tref = np.arange(0,per, 0.005)
    
    for tstart in 0,per:
        for t in tref:
            uref = np.array([0,0,-np.interp(t,tpoints,vpoints)])
            u = np.array(pulse.u(t + tstart, analysis='tran')).astype(float).reshape(3,)
            assert_array_almost_equal(u, uref)
           
def gen_stamps(toolkit=symbolic):
    circuit.default_toolkit = toolkit
    if toolkit.symbolic:
        R1,C1,L1,gain,gm,N = symbols('R1 C1 L1 gain gm N')
    else:
        R1=1.1e3
        C1=1e-12
        L1=1e-5
        gain=2.4
        gm=1.e-3
        N=1.2

    yield(R(1,gnd, r=R1), 1/R1 * np.array([[1, -1], [-1, 1]]), np.zeros((2,2)))

    yield(G(1,gnd, g=1/R1), 1/R1 * np.array([[1, -1], [-1, 1]]), 
          np.zeros((2,2)))

    yield(C(1,gnd, c=C1), np.zeros((2,2)), C1 * np.array([[1, -1], [-1, 1]]))

    GL = np.array([[0,0,1], [0,0,-1], [1, -1, 0]])
    CL = np.zeros((3,3), dtype=object)
    CL[2,2] = -L1
    yield(L(1,gnd, L=L1), GL, CL)

    GVCVS = np.array([[0,       0, 0,0, 0],
                      [0,       0, 0,0, 0],
                      [0,       0, 0,0, 1], 
                      [0,       0, 0,0,-1], 
                      [gain,-gain,-1,1, 0]])
    CVCVS = np.zeros((5,5))
    yield(VCVS(1, gnd, 2, gnd, g=gain),GVCVS,CVCVS)

    GVCCS = toolkit.zeros((4,4))
    GVCCS[2:4,0:2] =  np.array([[1, -1],[-1, 1]])
    yield(VCCS(1, gnd, 2, gnd, gm = gm), gm * GVCCS,
          np.zeros((4,4)))

    GNullor = np.array([[0,0,0, 0, 0],
                        [0,0,0, 0, 0],
                        [0,0,0, 0, 1], 
                        [0,0,0, 0,-1], 
                        [1,-1,0,0, 0]])
    yield(Nullor(1, gnd, 2, gnd), GNullor, np.zeros((5,5)))

    ## The primary column is `-1/N`.  It read `N, -N` until 2026-09-07, which
    ## got the VOLTAGE ratio right and the POWER ratio wrong by `N**2`; this
    ## fixture pinned the wrong stamp, so it is changed on purpose.  See
    ## `test_the_ideal_transformer_conserves_power`.
    GTransformer = np.array([[0,0,0, 0, -1/N],
                             [0,0,0, 0,  1/N],
                             [0,0,0, 0, 1], 
                             [0,0,0, 0,-1], 
                             [-1,1,N,-N, 0]])
    yield(Transformer(1, gnd, 2, gnd, n = N), GTransformer, np.zeros((5,5)))

    GGyrator = np.array([[ 0., 0., 1.,-1.],
                         [ 0., 0.,-1., 1.],
                         [-1., 1., 0., 0.],
                         [ 1.,-1., 0., 0.]])
    yield(Gyrator(1, gnd, 2, gnd, gm = 1.), GGyrator, np.zeros((4,4)))
    

def gen_stamps_sources(toolkit=symbolic):
    circuit.default_toolkit = toolkit
    if toolkit.symbolic:
        vac, phase = symbols('vac phase')
    else:
        vac = 1.2
        phase = 30

    v = vac * toolkit.exp(1j * toolkit.pi * phase / 180.)
    G = np.array([[0, 0, 1],
               [0, 0,-1],
               [1, -1, 0]])
    yield(VS(1,0,vac=vac, phase=phase), G, np.zeros((3,3)), toolkit.array([0,0,-v]))

    cir = SubCircuit()
    cir['vs'] = VS(1,0,vac=vac, phase=phase)
    yield(cir, G, np.zeros((3,3)), toolkit.array([0,0,-v]))

def test_stamp():
    circuit.default_toolkit = symbolic
    
    for toolkit in numeric, symbolic, symbolic_poly:
        for cir, G, C in gen_stamps(toolkit=toolkit):
            assert_array_equal(cir.G(np.zeros(cir.n)), G)
            assert_array_equal(cir.C(np.zeros(cir.n)), C)

        for cir, G, C, u in gen_stamps_sources(toolkit=toolkit):
            assert_array_equal(cir.G(np.zeros(cir.n)), G)
            assert_array_equal(cir.C(np.zeros(cir.n)), C)
            assert_array_equal(cir.u(np.zeros(cir.n), analysis='ac'), u)

def test_nullor_vva():
    """Test nullor element by building a V-V amplifier"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    c = SubCircuit()

    Vin = Symbol('Vin')
    R1 =Symbol('R1')
    R2 = Symbol('R2')
    
    nin = c.add_node('in')
    n1 = c.add_node('n1')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, vac=Vin)
    c['R1'] = R(n1, gnd, r=R1)
    c['R2'] = R(nout, n1, r=R2)
    c['nullor'] = Nullor(n1, nin, gnd, nout)
    
    result = AC(c, toolkit=symbolic).solve(Symbol('s'))
    
    vout = result.v(nout)

    assert simplify(vout - Vin * (R1 + R2) / R1) == 0, \
        'Did not get the expected result, %s != 0'% \
        str(simplify(vout - Vin * (R1 + R2) / R1))

def test_SVCVS_laplace_integrator():
    """Test SVCCS with a integrator transfer function

    """
    pycircuit.circuit.circuit.default_toolkit = symbolic

    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    a0,b0,Gdc = [sympy.Symbol(symname, real=True) for symname in
                 'a0,b0,Gdc'.split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = SVCVS( n1, gnd, n2, gnd,
                         denominator = (a0, 0),
                         numerator = (b0,))

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert sympy.expand(res.v(n2,gnd)) == sympy.expand(b0/(a0*s))

def test_SVCVS_laplace_n1_d2():
    """Test VCCS with a laplace defined transfer function first order numerator
    and second order denominator"""

    pycircuit.circuit.circuit.default_toolkit = symbolic
    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    b0,a0,a1,a2,Gdc = [sympy.Symbol(symname, real=True) for symname in
                       'b0,a0,a1,a2,Gdc'.split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = SVCVS( n1, gnd, n2, gnd,
                         denominator = (a0, a1, a2),
                         numerator   = (b0, 0))

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    ## Compare via full simplification; sympy.expand does not reduce the
    ## nested-fraction form these two algebraically-equal expressions take.
    assert sympy.simplify(res.v(n2,gnd) - b0*s/(a0*s*s+a1*s+a2)) == 0

def test_SVCVS_laplace_d3_n1():
    """Test VCCS with a laplace defined transfer function with second order
    numerator and third order denominator
    """

    pycircuit.circuit.circuit.default_toolkit = symbolic
    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    b0,a0,a1,a2,a3,Gdc = [sympy.Symbol(symname, real=True) for
                                symname in 'b0,a0,a1,a2,a3,Gdc'
                                .split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = SVCVS( n1, gnd, n2, gnd,
                        denominator = [a0, a1, a2, a3],
                        numerator   = [b0, 0, 0])

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert sympy.cancel(sympy.expand(res.v(n2,gnd))) == sympy.expand((b0*s*s)/(a0*s*s*s+a1*s*s+a2*s+a3))

def test_the_ideal_transformer_conserves_power():
    """The winding ratio must appear as `1/n` in the CURRENT stamp.

    The stamp carried `+n` there until 2026-09-07.  That gets the voltage
    ratio right -- `V_in/V_out = n` comes from the constraint ROW, which was
    never wrong -- and leaves `|P_out/P_in| = 1/n**2`.  So every check in the
    voltage domain passed, including this file's own stamp fixture and the
    element's doctest, both of which pinned the wrong matrix.  Only a POWER
    balance can see it, which is what this asserts.
    """
    pycircuit.circuit.circuit.default_toolkit = numeric

    r_s, r_l = 10.0, 1e3
    for ratio in (0.5, 1.0, 2.0, 10.0):
        cir = SubCircuit()
        src, a, b = cir.add_nodes('src', 'a', 'b')
        cir['vs'] = VS(src, gnd, v=1.0)
        cir['rs'] = R(src, a, r=r_s)          ## to MEASURE the primary current
        cir['t'] = Transformer(a, gnd, b, gnd, n=ratio)
        cir['rl'] = R(b, gnd, r=r_l)
        res = DC(cir).solve()

        vs, va = float(res.v(src, gnd)), float(res.v(a, gnd))
        vb = float(res.v(b, gnd))
        ## ⚠ BOTH currents are MEASURED from the solved circuit, through the
        ## series resistors -- never reconstructed from `1/n`.  Computing the
        ## primary current as `-i_br/n` instead makes this test assert the
        ## very relation it is supposed to check, and it then PASSES against
        ## the broken stamp.  It did, on the first attempt.
        i_in = (vs - va) / r_s
        i_out = vb / r_l
        p_in = va * i_in
        p_out = vb * i_out

        assert abs(va / vb - ratio) < 1e-9, \
            'n=%g: V_in/V_out is %.6f, want %g' % (ratio, va / vb, ratio)
        assert abs(p_out / p_in - 1.0) < 1e-9, \
            'n=%g: P_out/P_in is %.6f, want 1.0 (the old +n stamp gives %g)' \
            % (ratio, p_out / p_in, 1.0 / ratio ** 2)


def test_the_svcvs_keeps_fractional_coefficients_on_the_numeric_toolkit():
    """`SVCVS` built `G`/`C` with `dtype=int` and TRUNCATED every coefficient.

    Two things kept this alive.  Every other `SVCVS` test in this file runs on
    the SYMBOLIC toolkit, where an integer container holds sympy objects and
    cannot truncate -- a fixture that cannot express the effect.  And on the
    numeric toolkit the damage depends on which end of `denominator` the large
    coefficient sits: the coefficients are normalised by `den[0]`, so
    `(tau, 1)` yields `1000.0` (truncates to `1000`, invisible) while a
    denominator whose leading coefficient is the LARGEST normalises to
    entries below 1 that all truncate to ZERO, leaving a degenerate filter
    that still solves and still returns a number.

    So this test uses the second ordering deliberately, and checks the pole
    lands where the coefficients say -- not merely that some entry is
    fractional.
    """
    pycircuit.circuit.circuit.default_toolkit = numeric

    ## A first-order lowpass `1/(1 + s*tau)` whose NORMALISED denominator
    ## coefficient is `1/tau = 1e3`... and its reciprocal ordering, where the
    ## normalised coefficient is `tau = 1e-3` and truncated to zero.
    tau = 1e-3
    cir = SubCircuit()
    n1, n2 = cir.add_nodes('1', '2')
    cir['vs'] = VS(n1, gnd, vac=1.0)
    cir['h'] = SVCVS(n1, gnd, n2, gnd,
                     denominator=(tau, 1.0), numerator=(1.0,))

    f_pole = 1.0 / (2 * np.pi * tau)
    for f, want in ((f_pole, 1.0 / np.sqrt(2.0)),
                    (1e4, abs(1.0 / (1.0 + 1j * 2 * np.pi * 1e4 * tau)))):
        got = abs(complex(AC(cir).solve(freqs=f).v(n2, gnd)))
        assert abs(got - want) < 1e-4, \
            '|H| at %.4g Hz is %.6f, want %.6f' % (f, got, want)

    ## The direct statement of the defect, on the ordering that zeroes
    ## everything: fractional entries must SURVIVE into the stamp.
    e = SVCVS(0, 1, 2, 3, denominator=[1.0, 2.25e-6, 2.53e-12],
              numerator=[1.0])
    G = np.asarray(e.G(np.zeros(e.n)), dtype=float)
    C = np.asarray(e.C(np.zeros(e.n)), dtype=float)
    assert (np.abs(G - np.round(G)) > 0).any(), \
        'every G coefficient truncated to an integer: %r' % (G,)
    assert G.dtype.kind == 'f' and C.dtype.kind == 'f', \
        'G/C dtypes are %s/%s, must be floating point' % (G.dtype, C.dtype)


def test_Idt_sym():
    """Test integrator element symbolically"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    c = SubCircuit()

    Vin = Symbol('Vin')
    R1 = Symbol('R1')
    
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, vac=Vin)
    c['R1'] = R(nout, gnd, r=R1)
    c['Idt'] = Idt(nin, gnd, nout, gnd)
    
    result = AC(c, toolkit=symbolic).solve(Symbol('s'),complexfreq=True)
    
    vtr = simplify(result.v(nout)/result.v(nin))
    assert vtr == 1/Symbol('s')

def test_Idt_tran():
    """Test integrator element in transient"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    c = SubCircuit()
    nin = c.add_node('in')
    nout = c.add_node('out')

    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['Idt'] = Idt(nin, gnd, nout, gnd)

    ## uic=True, deliberately: an ideal integrator has NO DC operating point --
    ## its output is the unbounded integral of a constant input, and with all
    ## sources zeroed the matrix is structurally singular.  This ran before only
    ## because a failed DC was silently replaced by zeros; zeros happens to be
    ## the right initial state for an integrator, so the fix is to ask for it.
    tran = Transient(c, toolkit=numeric, uic=True)
    result = tran.solve(tend=0.5,timestep=1e-2, fixed_timestep=True)
    y = result.v(nout).y
    x = result.v(nout).x[0]
    # vout = vin * t with constant input
    # Transient solver time alignment is now correct: y(t) evaluates at t
    assert_array_almost_equal(y[1:], x[1:])

def test_Idtmod_sym():
    """Test modulus integrator element symbolically"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    c = SubCircuit()

    Vin = Symbol('Vin')
    R1 = Symbol('R1')
    
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, vac=Vin)
    c['R1'] = R(nout, gnd, r=R1)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd)
    
    result = AC(c, toolkit=symbolic).solve(Symbol('s'),complexfreq=True)
    
    vtr = simplify(result.v(nout)/result.v(nin))
    assert vtr == 1/Symbol('s')


def test_Idtmod_tran():
    """Test modulo integrator element in transient"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    c = SubCircuit()
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus = 1., offset = -0.)
    
    ## uic=True, deliberately: an ideal integrator has NO DC operating point --
    ## its output is the unbounded integral of a constant input, and with all
    ## sources zeroed the matrix is structurally singular.  This ran before only
    ## because a failed DC was silently replaced by zeros; zeros happens to be
    ## the right initial state for an integrator, so the fix is to ask for it.
    tran = Transient(c, toolkit=numeric, uic=True)
    result = tran.solve(tend=0.5,timestep=1e-2, fixed_timestep=True)
    y = result.v(nout).y
    x = result.v(nout).x[0]
    # vout = vin * t with constant input
    # Transient solver time alignment is now correct: y(t) evaluates at t
    assert_array_almost_equal(y[1:], x[1:])

def test_Idtmod_modulo():
    """Test modulo integrator element in transient"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    c = SubCircuit()
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus = 1., offset = -0.)
    
    ## uic=True, deliberately: an ideal integrator has NO DC operating point --
    ## its output is the unbounded integral of a constant input, and with all
    ## sources zeroed the matrix is structurally singular.  This ran before only
    ## because a failed DC was silently replaced by zeros; zeros happens to be
    ## the right initial state for an integrator, so the fix is to ask for it.
    ## integrator pinned at P6: the sample AT the modulo wrap is a
    ## left/right-limit convention, and this element's record (y(1.0) = 0.0,
    ## the right limit) was written under Euler; Gear-2's two-point history
    ## lands the left limit there instead.
    from pycircuit.circuit.integrator import EulerIntegrator
    tran = Transient(c, toolkit=numeric, uic=True,
                     integrator=EulerIntegrator())
    result = tran.solve(tend=2.0,timestep=1e-2, fixed_timestep=True)
    y = result.v(nout).y
    x = result.v(nout).x[0]
    # vout = vin * t with constant input
    # Transient solver time alignment is now correct: y(t) evaluates at t
    ## Two-sided (congruence) comparison: at a sample landing exactly ON the
    ## wrap the sawtooth is double-valued, and WHICH limit the solver reports
    ## is decided by sub-ulp rounding.  The Phase-2 gauge shift (idtmod.md
    ## 5.2) keeps the state bounded -- more accurate than the accumulated `x`
    ## grid it is compared against -- so the two no longer share correlated
    ## drift and a boundary sample can land either side.  Distance modulo
    ## the modulus treats both limits as the same point.
    d = np.abs(y[1:] - x[1:] % 1.0)
    d = np.minimum(d, 1.0 - d)
    assert_array_almost_equal(d, np.zeros_like(d))

if __name__ == '__main__':
    test_nullor_vva()


def test_the_ccvs_outputs_its_transresistance_times_the_input_current():
    """2026-09-08.  `CCVS` stamped its transresistance at (input KVL row,
    output current column): `v_in = r i_out` with the output branch a 0 V
    short -- a transposed stamp, so the element output ZERO volts for any
    input current.  It survived because its doctest pinned the wrong matrix
    and the only other test compared two backends' copies of it.  Found when
    a relaxation oscillator sensing its capacitor current through a CCVS
    never oscillated at any grid.  Pinned by MEASUREMENT: a 1 V source
    through 1 kOhm into the input branch (1 mA) must give r * 1 mA at the
    output, and the input branch must be a 0 V ammeter.
    """
    import numpy as np
    from pycircuit.circuit import circuit
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.elements import VS, R, CCVS
    from pycircuit.circuit.dcanalysis import DC
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    for n in ('s', 'x', 'y'):
        c.add_node(n)
    c['vs'] = VS('s', gnd, v=1.0)
    c['R'] = R('s', 'x', r=1e3)
    c['amm'] = CCVS('x', gnd, 'y', gnd, r=2.5e3)
    c['Ry'] = R('y', gnd, r=1e6)
    res = DC(c).solve()
    assert abs(float(res.v('x'))) < 1e-12, 'the input branch must be a 0 V ammeter'
    assert abs(float(res.v('y')) - 2.5) < 1e-9, \
        'v(y) = %.4f V for 1 mA through r = 2.5 kOhm (the transposed stamp gave 0.0000)' % float(res.v('y'))


def test_every_classical_controlled_source_meets_its_defining_relation_at_dc():
    """2026-09-08, after the CCVS's transposed stamp: a DC measurement of each
    controlled source and two-port against the relation that DEFINES it, with
    an independent stimulus -- not a doctest pinning a matrix, not a backend
    comparing two copies of the same element.  Signs are asserted where the
    element's convention is documented by its own doctest, magnitudes
    everywhere.
    """
    import numpy as np
    from pycircuit.circuit import circuit
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.elements import (VS, R, VCVS, VCCS, CCCS, CCVS,
                                            Transformer, Gyrator, Nullor)
    from pycircuit.circuit.dcanalysis import DC
    circuit.default_toolkit = circuit.numeric

    def nodes(c, *names):
        for n in names:
            c.add_node(n)

    c = SubCircuit(); nodes(c, 'i', 'o')
    c['vs'] = VS('i', gnd, v=1.0); c['e'] = VCVS('i', gnd, 'o', gnd, g=3.0); c['RL'] = R('o', gnd, r=1e3)
    assert abs(float(DC(c).solve().v('o')) - 3.0) < 1e-9, 'VCVS: v_out = g v_in'

    c = SubCircuit(); nodes(c, 'i', 'o')
    c['vs'] = VS('i', gnd, v=1.0); c['g'] = VCCS('i', gnd, 'o', gnd, gm=2e-3); c['RL'] = R('o', gnd, r=1e3)
    assert abs(abs(float(DC(c).solve().v('o'))) - 2.0) < 1e-9, 'VCCS: |i_out| = gm v_in'

    c = SubCircuit(); nodes(c, 's', 'x', 'o')
    c['vs'] = VS('s', gnd, v=1.0); c['R'] = R('s', 'x', r=1e3); c['RL'] = R('o', gnd, r=1e3)
    c['f'] = CCCS('x', gnd, 'o', gnd, F=4.0)
    r = DC(c).solve()
    assert abs(float(r.v('x'))) < 1e-9 and abs(abs(float(r.v('o'))) - 4.0) < 1e-9, 'CCCS: input a short, |i_out| = F i_in'

    c = SubCircuit(); nodes(c, 's', 'x', 'o')
    c['vs'] = VS('s', gnd, v=1.0); c['R'] = R('s', 'x', r=1e3); c['h'] = CCVS('x', gnd, 'o', gnd, r=2.5e3); c['RL'] = R('o', gnd, r=1e6)
    assert abs(float(DC(c).solve().v('o')) - 2.5) < 1e-9, 'CCVS: v_out = r i_in (was 0.0000 before 2026-09-08)'

    c = SubCircuit(); nodes(c, 'i', 'o')
    c['vs'] = VS('i', gnd, v=1.0); c['t'] = Transformer('i', gnd, 'o', gnd, n=2.0); c['RL'] = R('o', gnd, r=1e3)
    assert abs(float(DC(c).solve().v('o')) - 0.5) < 1e-9, 'Transformer: v_out = v_in / n'

    c = SubCircuit(); nodes(c, 'i', 'o')
    c['vs'] = VS('i', gnd, v=1.0); c['gy'] = Gyrator('i', gnd, 'o', gnd, gm=2e-3); c['RL'] = R('o', gnd, r=1e3)
    assert abs(abs(float(DC(c).solve().v('o'))) - 2.0) < 1e-9, 'Gyrator: |i_out| = gm v_in'

    c = SubCircuit(); nodes(c, 's', 'm', 'o')
    c['vs'] = VS('s', gnd, v=1.0); c['R1'] = R('s', 'm', r=1e3); c['R2'] = R('m', 'o', r=5e3); c['nul'] = Nullor(gnd, 'm', 'o', gnd)
    r = DC(c).solve()
    assert abs(float(r.v('o')) + 5.0) < 1e-9 and abs(float(r.v('m'))) < 1e-9, 'Nullor as an ideal inverting amplifier: v_out = -R2/R1 v_in, v_- = 0'
