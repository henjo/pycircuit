# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.circuit.circuit 
from pycircuit.circuit import *
from pycircuit.circuit.tests._warnpolicy import quiet
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from .test_circuit import create_current_divider
import unittest

def test_integer_component_values():
    """Test dc analysis with integer component values
    
       As python per default uses integer arithmetics integer
       component values can lead to problems. (Python < 3.0)
    """
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit(toolkit=numeric)

    c['vs'] = VS('net1', gnd, v = 9)
    c['R1'] = R( 'net1',  'net2', r = 50)
    c['R2'] = R( 'net2', gnd, r = 50)

    dc = DC(c)
    res = dc.solve()
    
    assert res.v('net2') == 4.5

    assert res.i('R2.plus') == 0.09

def TODOtest_noise_dc_steady_state():
    """Test that dc-steady state is accounted for in noise simulations
    """
    pass

def test_noise_with_frequency_vector():
    """Test that noise analysis support an array as input argument for frequency

    """
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit(toolkit=numeric)

    n1,n2 = c.add_nodes('net1', 'net2')

    c['vs'] = VS(n1, gnd, v = 9.)
    c['R1'] = R( n1,  n2, r = 50.)
    c['R2'] = R( n2, gnd, r = 50.)
    
    noise = Noise(c, inputsrc='vs', outputnodes=(n2, gnd))
    should = np.array([noise.solve(0)['Svnout'],noise.solve(1)['Svnout']])
    res = noise.solve(np.array([0,1]))
    assert_array_equal(res['Svnout'], should)


def test_a_frequency_dependent_CY_sweeps():
    """⚠ FOUND BY THE PRIVATE COMPARISON SUITE (2026-09-05): a compact
    model's `CY` came out RAGGED under an array frequency -- the flicker
    entry array-valued, the thermal entries scalar, at `kf = 0` too since
    the term is emitted unconditionally -- and the small-signal analysis
    handed `CY` the whole sweep, so `Noise(...).solve(freqs=<array>)`
    failed on every compact model while a scalar worked.  Two things were
    wrong and both are fixed: the generated `CY` broadcasts its entries,
    and the analysis evaluates `CY` PER FREQUENCY (the assembly takes
    scalar entries, so a handwritten coloured source was blocked by the
    second one on its own).  Gates: a coloured `IS` and a level-1 MOS
    sweep to exactly the per-frequency values, nonzero, and the three
    compact models return `(n, n, nf)` under an array.
    """
    import warnings
    import numpy as np
    from pycircuit.circuit import SubCircuit, R, C, gnd, AC, Noise
    from pycircuit.circuit.simwarnings import SimulationWarning
    from pycircuit.circuit.elements import VS, IS
    from pycircuit.circuit.elements_hdl import (MosLevel1Hdl, GummelPoonNpnHdl,
                                                EkvNmosHdl)
    circuit.default_toolkit = circuit.numeric
    warnings.simplefilter('ignore')     ## (scoped to this test by pytest)
    warnings.simplefilter('error', SimulationWarning)
    fr = np.array([1e3, 1e5, 1e7])

    def sweep(cir, src, out):
        ra = Noise(cir, inputsrc=src, outputnodes=(out, gnd)).solve(freqs=fr)
        sa = np.asarray(ra['Svnout'], dtype=float).ravel()
        ss = np.array([float(np.asarray(Noise(cir, inputsrc=src,
                                              outputnodes=(out, gnd))
                                        .solve(freqs=f)['Svnout']).ravel()[0])
                       for f in fr])
        assert np.all(ss > 0.0), 'a degenerate fixture proves nothing'
        assert np.max(np.abs(sa / ss - 1.0)) < 1e-12, (sa, ss)
        return sa
    cir = SubCircuit()
    cir.add_node('in')
    cir.add_node('a')
    cir['V'] = VS('in', gnd, v=1.0, vac=1.0)
    cir['R'] = R('in', 'a', r=1e3)
    cir['C'] = C('a', gnd, c=1e-9)
    cir['n'] = IS('a', gnd, i=0.0, noisePSD=1e-20, noiseTau=1e-6)
    sa = sweep(cir, 'V', 'a')
    assert sa[0] > 10.0 * sa[2], 'the coloured source must roll off'
    mos = SubCircuit()
    mos.add_node('d')
    mos.add_node('g')
    mos.add_node('vdd')
    mos['Vg'] = VS('g', gnd, v=2.0, vac=1.0)
    mos['Vd'] = VS('vdd', gnd, v=3.0)
    mos['RL'] = R('vdd', 'd', r=1e3)
    mos['M'] = MosLevel1Hdl('d', 'g', gnd, gnd)
    AC(mos).solve(freqs=fr)
    sweep(mos, 'Vg', 'd')
    for cls, args in ((MosLevel1Hdl, ('d', 'g', 's', 'b')),
                      (GummelPoonNpnHdl, ('c', 'b', 'e')),
                      (EkvNmosHdl, ('d', 'g', 's', 'b'))):
        el = cls(*args)
        n = len(el.terminals)
        cy = np.asarray(el.CY(np.linspace(0.1, 1.0, n), np.array([1.0, 10.0])))
        assert cy.shape == (n, n, 2), (cls.__name__, cy.shape)


def test_ac_and_noise_at_a_given_operating_point_read_the_devices_there():
    """AC and noise at an operating point handed in (`dcx`) read a stateful
    limiter's device AT that point (`dc_steady_state`), as they do at the one
    their own DC solve finds.  ⚠ A `Diode` reads `G` as the tangent at its
    stored `_vlim`, which a DC sweep leaves at its LAST point: AC at an
    earlier bias point of the sweep read 0.21x the gain and 0.045x the
    output noise (2026-09-28).  The transimpedance analysis evaluates at
    its own fixed point, and reads the devices there too."""
    from pycircuit.circuit.dcanalysis import DC, DCSweep
    from pycircuit.circuit.analysis_ss import AC, Noise, TransimpedanceAnalysis
    circuit.default_toolkit = numeric

    def build(v=0.0):
        c = SubCircuit()
        c['vs'] = VS(1, gnd, v=v, vac=1.0)
        c['R'] = R(1, 2, r=1e3)
        c['D'] = Diode(2, gnd)
        return c

    def small_signal(c, dcx=None):
        kw = {} if dcx is None else {'dcx': dcx}
        g = complex(np.asarray(AC(c, **kw).solve(np.array([1e3])).v(2, gnd)).ravel()[0])
        n = float(np.real(Noise(c, inputsrc='vs', outputnodes=(2, gnd),
                                **kw).solve(1e3)['Svnout']))
        return g, n
    vals = np.linspace(0.0, 5.0, 11)
    c = build()
    with quiet():
        x = np.asarray(DCSweep(c).solve('vs', 'v', vals).x, dtype=float)
    ## the reference: the same point on a circuit no analysis has touched
    g, n = small_signal(c, dcx=x[:, 3])
    g_ref, n_ref = small_signal(build(vals[3]), dcx=x[:, 3])
    assert abs(g / g_ref - 1.0) < 1e-12 and abs(n / n_ref - 1.0) < 1e-12, \
        (abs(g) / abs(g_ref), n / n_ref)

    def zin(c):
        res = TransimpedanceAnalysis(c).solve(1e3, [Branch(2, gnd)])
        return np.asarray(res[0], dtype=complex) if isinstance(res, (list, tuple)) \
            else np.asarray(res, dtype=complex)
    c = build(5.0)
    DC(c).solve()                 # the diode's state now at a forward bias
    assert np.allclose(zin(c), zin(build(5.0)), rtol=1e-12, atol=0.0)
