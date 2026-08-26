# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""The fifth Phase-6 library batch: the VCO with real phase noise, the PLL's
small behavioural blocks, the optical devices, the MESFET/HEMT family and
SPICE MOS level 3 -- the first models written with ``params_as``.

====================  ====================================================
model                 what it is, and the reference it is held to
====================  ====================================================
`VcoHdl`              noise injected into the FREQUENCY; the phase-noise
                      PSD is asserted against ``sf/(2*pi*f)^2`` and
                      ``kff/f/(2*pi*f)^2`` -- Leeson's slopes from the
                      integral, to 1e-9 -- and the free-running frequency
                      and phase accumulation against ``f0 + kvco*vc`` and
                      ``idtmod``'s modulus off a transient
`DividerHdl`          phase-domain ``/N``; output frequency by counting
                      crossings, and `Cross` measured against the same
                      model without it (the second production user)
`MixerHdl`            the trigonometric identity, by projection
`ChargePumpHdl`       UP/DN currents, both off, both on
`PhotodiodeHdl`       the single-diode solar cell closed form: ``Isc``,
                      ``Voc = n*Vt*ln(1 + Iph/IS)``, and the implicit
                      I-V with ``rs`` and ``rsh`` solved independently
`LedHdl`              optical output linear in current above threshold;
                      an optocoupler's CTR as ``resp*eta``
`MesfetCurticeHdl`    Curtice 1980 eq. (1), transcribed
`MesfetStatzHdl`      Statz 1987 via `mesload.c`, both modes, gate diodes
`HemtAngelovHdl`      Angelov 1992 eqs. (1)-(3), and the peak-gm identity
`MosLevel3Hdl`        level 1 in the limit where the two coincide
`MosLevel3PmosHdl`    (``gamma = 0``, every level-3 knob off), and an
                      independent numpy transcription of `mos3load.c`
                      with every knob on, both modes, forward bulk
====================  ====================================================

Every Jacobian claim goes through `hdl.check_jacobians`; every limiter
is measured against its absence on the same circuit; the mutation log is
in the batch report.
"""
import math
import warnings

import numpy as np
import pytest
import sympy
from numpy.testing import assert_allclose

import pycircuit.circuit.circuit
from pycircuit.circuit.circuit import defaultepar, SubCircuit
from pycircuit.circuit import numeric, gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.analysis_ss import Noise
from pycircuit.circuit.elements import VS, R, IS as ISRC
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution, Node,
                                   Cross, TEMP, check_jacobians, explain,
                                   x_layout, var, idtmod, KBOLTZMANN,
                                   QELECTRON)
from pycircuit.utilities.param import Parameter

_KB = float(KBOLTZMANN)
_QE = float(QELECTRON)
_T0 = float(defaultepar.T)
_UT = _KB * _T0 / _QE
_EPSOX = 3.9 * 8.854187817e-12
_EPSSI = 11.7 * 8.854187817e-12
_NI = 1.45e10


@pytest.fixture(autouse=True)
def _numeric_toolkit():
    old = pycircuit.circuit.circuit.default_toolkit
    pycircuit.circuit.circuit.default_toolkit = numeric
    yield
    pycircuit.circuit.circuit.default_toolkit = old


def _mk(cls, *nodes, **kw):
    el = cls(*[Node(nm) for nm in nodes], **kw)
    el.update_iparv()
    return el


def _dc(c):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return DC(c).solve()


def _crossings(t, v):
    """Linearly interpolated zero-crossing times of a sampled waveform."""
    k = np.nonzero(np.diff(np.sign(v)) != 0)[0]
    return t[k] - v[k] * (t[k + 1] - t[k]) / (v[k + 1] - v[k])


def _freq_of(t, v):
    z = _crossings(t, v)
    return (len(z) - 1) / (2.0 * (z[-1] - z[0]))


## ======================================================================
## 1.  The VCO
## ======================================================================

def _vco_transient(vc, f0=1e4, kvco=1e3, tend=2e-3, **kw):
    c = SubCircuit()
    c['vc'] = VS('vc', gnd, v=vc)
    c['X'] = eh.VcoHdl('vc', gnd, 'out', gnd, 'ph', f0=f0, kvco=kvco, **kw)
    c['R'] = R('out', gnd, r=1e3)
    tr = Transient(c, toolkit=numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(tend=tend, timestep=1e-7)
    return (tr.statistics, np.asarray(res.v('out').x[0], float),
            np.asarray(res.v('out').y, float), np.asarray(res.v('ph').y, float))


def test_vco_free_running_frequency_follows_kvco():
    """``f = f0 + kvco*vc``, measured by counting output crossings at two
    control voltages, and ``kvco`` recovered as their difference.  The
    controller's steps are coarse against the period (about 15 per
    cycle), so the crossing-count frequency is good to 0.1%; asserted
    at 0.5%, and ``kvco`` from the difference at 3%."""
    _s, t1, v1, _p = _vco_transient(0.5)
    _s, t2, v2, _p = _vco_transient(-0.5)
    f1, f2 = _freq_of(t1, v1), _freq_of(t2, v2)
    assert_allclose(f1, 10500.0, rtol=5e-3)
    assert_allclose(f2, 9500.0, rtol=5e-3)
    assert_allclose(f1 - f2, 1e3, rtol=3e-2)


def test_vco_phase_accumulates_and_wraps_at_the_modulus():
    """The phase terminal is ``(f*t) mod modulus`` in cycles: a sawtooth
    whose value at ``t`` is the accumulated cycle count.  With
    ``modulus = 4`` at 10.5 kHz, ``t = 1 ms`` is 10.5 cycles, i.e. 2.5
    after the wrap; ``modulus = 1`` gives 0.5.  Both asserted, so a
    modulus that was ignored (always 1) fails the first."""
    for modulus, expect in ((4.0, 2.5), (1.0, 0.5)):
        _s, t, _v, ph = _vco_transient(0.5, modulus=modulus)
        assert ph.min() >= 0.0 and ph.max() < modulus
        assert_allclose(np.interp(1e-3, t, ph), expect, atol=2e-3)
        ## the sawtooth's slope is the frequency, between wraps
        k = (t > 0.2e-3) & (t < 0.3e-3)
        slope = np.polyfit(t[k], ph[k], 1)[0]
        assert_allclose(slope, 10500.0, rtol=1e-3)


def test_vco_noise_is_on_the_phase_integral():
    """THE point of the model.  A noise analysis at the output gives a
    PSD of ``(2*pi*va)^2 * sf/(2*pi*f)^2`` for white frequency noise
    and ``(2*pi*va)^2 * kff/f/(2*pi*f)^2`` for flicker -- the -20 and
    -30 dB/decade slopes of Leeson's model, produced by the integral
    and not by a shaped source.  Asserted to 1e-9 at three decades,
    from the same adjoint noise analysis every other element uses.

    The control: the noise's ONLY entry in ``CY`` is the internal
    frequency node's row -- nothing on the output, nothing on the
    phase terminal -- so a source that had been added to the output
    (flat) or to the phase (flat in phase) would fail."""
    va = 1.0
    for sf, kff in ((1e-3, 0.0), (0.0, 1.0), (2e-3, 0.5)):
        c = SubCircuit()
        c['vc'] = VS('vc', gnd, v=0.5, vac=1.0)
        c['X'] = eh.VcoHdl('vc', gnd, 'out', gnd, 'ph', f0=1e4, kvco=1e3,
                           va=va, sf=sf, kff=kff)
        c['R'] = R('out', gnd, r=1e3)
        for f in (1e2, 1e3, 1e4):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                r = Noise(c, inputsrc='vc', outputnodes=('out', gnd)).solve(f)
            w = 2 * np.pi * f
            expect = (2 * np.pi * va) ** 2 * (sf + kff / f) / w ** 2
            assert_allclose(float(np.real(r['Svnout'])), expect, rtol=1e-9)
    el = _mk(eh.VcoHdl, 'cp', 'cn', 'outp', 'outn', 'ph', sf=1e-3, kff=1.0)
    lay = {nm: k for k, nm, _kind in x_layout(el)}
    CY = np.asarray(el.CY(np.zeros(len(lay)), 2 * np.pi * 10.0), float)
    nz = {(i, j) for i, j in zip(*np.nonzero(CY))}
    assert nz == {(lay['fn'], lay['fn']), (lay['fn'], lay['outn']),
                  (lay['outn'], lay['fn']), (lay['outn'], lay['outn'])}, nz
    assert_allclose(CY[lay['fn'], lay['fn']], 1e-3 + 1.0 / 10.0, rtol=1e-12)


class _ChainedIdtmod(Behavioural):
    """A `var()` model with a pinned `idtmod`: the DSL bug the VCO found
    (the chained path ignored the DC pin) in its smallest form."""
    instparams = [Parameter(name='k', desc='', unit='', default=2.0)]

    @staticmethod
    def analog(ip, im, op, om):
        bi, bo = Branch(ip, im), Branch(op, om)
        f = var(k * bi.V, 'f')                                         # noqa
        return Contribution(bo.V, idtmod(f, 0.25, 1.0, 0.0))


def test_a_chained_idtmod_is_pinned_at_dc():
    """Regression for the bug `VcoHdl` found: on the chained path
    ``i``/``G`` returned the transient stamp at DC, so the state row was
    ``ds/dt = f`` instead of ``s - ic`` and the operating point was
    singular.  `IdtmodHdl` is flat and never saw it.  The DC solve must
    converge and land the state on its (wrapped) ``ic``."""
    c = SubCircuit()
    c['vs'] = VS('a', gnd, v=3.0)
    c['X'] = _ChainedIdtmod('a', gnd, 'o', gnd)
    c['R'] = R('o', gnd, r=1e3)
    r = _dc(c)
    assert_allclose(float(r.v('o')), 0.25, atol=1e-9)


## ======================================================================
## 2.  Divider, mixer, charge pump
## ======================================================================

class _DividerNoCross(Behavioural):
    """`DividerHdl` without the `Cross`: the control for what `Cross`
    is worth."""
    params_as = 'p'
    instparams = [Parameter(name='n', desc='', unit='', default=4.0),
                  Parameter(name='gain', desc='', unit='', default=50.0),
                  Parameter(name='voh', desc='', unit='V', default=1.0),
                  Parameter(name='vol', desc='', unit='V', default=-1.0)]

    @staticmethod
    def analog(p, inp, inn, outp, outn):
        bi, bo = Branch(inp, inn), Branch(outp, outn)
        s = var(sympy.sin(2 * sympy.pi * bi.V / p.n), 's')
        return Contribution(bo.V, p.vol + (p.voh - p.vol) / 2
                            * (1 + sympy.tanh(p.gain * s)))


def _divider_transient(cls, n=4.0, f=10500.0, tend=2e-3):
    c = SubCircuit()
    c['vc'] = VS('vc', gnd, v=0.5)
    c['X'] = eh.VcoHdl('vc', gnd, 'out', gnd, 'ph', f0=f - 500.0, kvco=1e3,
                       modulus=n)
    c['R'] = R('out', gnd, r=1e3)
    c['D'] = cls('ph', gnd, 'div', gnd, n=n)
    c['R2'] = R('div', gnd, r=1e3)
    tr = Transient(c, toolkit=numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(tend=tend, timestep=1e-7)
    return (tr.statistics, np.asarray(res.v('div').x[0], float),
            np.asarray(res.v('div').y, float))


def test_divider_output_frequency_is_the_input_over_n():
    """Crossings counted on the divided output: ``f/n`` at two ratios.
    The rising and falling edges both count, so ``2*f/n`` per second."""
    for n in (4.0, 7.0):
        _s, t, v = _divider_transient(eh.DividerHdl, n=n)
        assert_allclose(_freq_of(t, v), 10500.0 / n, rtol=1e-2)


def test_divider_cross_is_measured_against_its_absence():
    """`Cross`'s second production user, measured the way the first was
    (`test_cross_lands_the_comparator_edges`): the distance from each
    analytic edge time ``k*n/(2f)`` to the nearest accepted timepoint,
    with and without the declaration, and the step cost.

    Recorded at introduction (n = 4, f = 10.5 kHz, 2 ms), distance to
    the nearest accepted point over the first eight edges::

        without Cross   6.0e-6  1.6e-5  1.5e-5  5.4e-6  4.1e-6  1.4e-5  1.7e-5  7.3e-6
        with Cross      6.0e-6  8.2e-6  6.7e-6  7.7e-6  7.0e-6  7.5e-6  7.1e-6  7.4e-6

    Nine breakpoints for eleven edges, 75 accepted steps against 73.
    The worst edge improves 2.3x and the scatter closes to 6..8 us --
    ONE first-order prediction over a 2.7e-5 s step on a 95 us period,
    never refined.  Batch 2's comparator had the same predictor land
    3e-16 s on three edges of four: there the controller polled a
    second time after the breakpoint, here it never does, so the
    ten-decade accuracy scatter batch 2 recorded is, on this model, a
    uniform 7 us -- and the best no-Cross edge (4.1e-6) beats the best
    Cross edge.  Cost: 1.03x steps, not 3x; the comparator's 3x came
    from a latching output restarting the controller small at each of
    its four edges, and the divider's tanh edge does not.  So `Cross`
    on its second user: brackets to first order, no convergence, cheap.
    Assertions are bounds on what was measured, so a regression in
    either direction shows.
    """
    n, f = 4.0, 10500.0
    s_on, t_on, _v = _divider_transient(eh.DividerHdl, n=n, f=f)
    s_off, t_off, _v = _divider_transient(_DividerNoCross, n=n, f=f)
    edges = np.arange(1, 9) * n / (2 * f)
    near_on = [float(np.min(np.abs(t_on - x))) for x in edges]
    near_off = [float(np.min(np.abs(t_off - x))) for x in edges]
    assert s_off.breakpoints_hit == 0
    assert s_on.breakpoints_hit >= 6, s_on.breakpoints_hit
    assert max(near_on) < max(near_off) / 1.5, (near_on, near_off)
    assert max(near_on) < 1.5e-5, near_on
    assert s_on.accepted_steps < 1.5 * s_off.accepted_steps, (
        s_on.accepted_steps, s_off.accepted_steps)


def test_mixer_products_are_the_trigonometric_identity():
    """``A*sin(w1 t) * B*sin(w2 t) = AB/2 [cos((w1-w2)t) - cos((w1+w2)t)]``:
    the transient output projected on the four candidate frequencies.
    Both products at ``k*A*B/2`` to 1% on fixed 2 us steps, nothing at
    the inputs' own frequencies (the leakage a wrong model would show)
    below 1% of it.
    """
    from pycircuit.circuit.elements import VSin
    f1, f2, A, B, k = 5e3, 3e3, 1.0, 0.5, 2.0
    c = SubCircuit()
    c['v1'] = VSin('rf', gnd, va=A, freq=f1)
    c['v2'] = VSin('lo', gnd, va=B, freq=f2)
    c['X'] = eh.MixerHdl('rf', gnd, 'lo', gnd, 'if', gnd, k=k)
    c['R'] = R('if', gnd, r=1e3)
    tr = Transient(c, toolkit=numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ## fixed 2 us steps: 62 per period of the 8 kHz sum, so the
        ## projection measures the model and not the step controller
        res = tr.solve(tend=2e-3, timestep=2e-6, fixed_timestep=True)
    t = np.asarray(res.v('if').x[0], float)
    v = np.asarray(res.v('if').y, float)
    ## one period of the difference frequency (2 kHz), uniform resample
    tu = np.linspace(0.5e-3, 1.5e-3, 4001)[:-1]
    vu = np.interp(tu, t, v)

    def amp(fr):
        return 2 * abs(np.mean(vu * np.exp(-2j * np.pi * fr * tu)))
    expect = k * A * B / 2
    assert_allclose(amp(f1 - f2), expect, rtol=1e-2)
    assert_allclose(amp(f1 + f2), expect, rtol=1e-2)
    assert amp(f1) < 0.01 * expect and amp(f2) < 0.01 * expect
    assert abs(np.mean(vu)) < 0.01 * expect


def test_charge_pump_up_and_down_currents():
    """UP sources ``iup`` into the output, DN sinks ``idn``, both on is
    the difference, both off is nothing: measured as the voltage across
    a 1 kOhm load.  ``iup != idn`` deliberately, so that a pump with
    the two swapped is caught."""
    iup, idn = 100e-6, 80e-6
    got = {}
    for vu, vd in ((0, 0), (1, 0), (0, 1), (1, 1)):
        c = SubCircuit()
        c['vu'] = VS('u', gnd, v=vu)
        c['vd'] = VS('d', gnd, v=vd)
        c['X'] = eh.ChargePumpHdl('u', 'd', 'o', gnd, iup=iup, idn=idn)
        c['R'] = R('o', gnd, r=1e3)
        got[(vu, vd)] = float(_dc(c).v('o')) / 1e3
    assert abs(got[(0, 0)]) < 1e-12
    assert_allclose(got[(1, 0)], iup, rtol=1e-6)
    assert_allclose(got[(0, 1)], -idn, rtol=1e-6)
    assert_allclose(got[(1, 1)], iup - idn, rtol=1e-6)


## ======================================================================
## 3.  Photodiode / solar cell / LED
## ======================================================================

#: A silicon cell's card: nA saturation current, n = 1.3, a small
#: series resistance and a finite shunt.  ``resp`` scaled so that
#: ``V(opt) = 1`` (one sun) gives ``Iph = 0.4 A``.  ``tnom`` is the
#: ambient, so the temperature path is the identity here.
CELL = dict(IS=2e-9, n=1.3, rs=0.02, rsh=50.0, resp=0.4, tnom=_T0)


def _cell_iv(v, sun, card=CELL):
    """The single-diode cell, solved independently by bisection on ``I``:
    ``I = Iph - IS*(exp((V + I*rs)/(n*Vt)) - 1) - (V + I*rs)/rsh``,
    ``I`` positive OUT of the anode."""
    iph = card['resp'] * sun
    n_, is_, rs, rsh = card['n'], card['IS'], card['rs'], card['rsh']

    def f(i):
        vj = v + i * rs
        return iph - is_ * (math.exp(vj / (n_ * _UT)) - 1) - vj / rsh - i
    lo, hi = -10.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _cell_current(v, sun, card=CELL):
    """The model's anode current, positive out of the anode, at bias v."""
    c = SubCircuit()
    c['vs'] = VS('a', gnd, v=v)
    c['opt'] = VS('opt', gnd, v=sun)
    c['X'] = eh.PhotodiodeHdl('a', gnd, 'opt', gnd, **card)
    r = _dc(c)
    ## the source's branch current flows into its plus terminal; a lit
    ## cell pushes current out of the anode INTO the source
    return float(r.i('vs.plus'))


def test_photodiode_dark_is_the_spice_diode():
    """With the optical port at zero the model must be `DiodeSpiceHdl`
    with a 1e30 ohm shunt, i.e. the same diode: currents equal to 1e-12
    at four biases, both signs."""
    for v in (-2.0, 0.0, 0.4, 0.6):
        c = SubCircuit()
        c['vs'] = VS('a', gnd, v=v)
        c['D'] = eh.DiodeSpiceHdl('a', gnd, IS=2e-9, n=1.3, rs=0.02,
                                  tnom=_T0)
        ref = float(_dc(c).i('vs.plus'))
        got = _cell_current(v, 0.0, dict(CELL, rsh=1e30))
        assert_allclose(got, ref, rtol=1e-12, atol=1e-30)


def test_solar_cell_isc_and_voc_are_the_closed_forms():
    """The ideal cell (``rs = 0``, ``rsh = inf``): ``Isc = Iph`` exactly
    and ``Voc = n*Vt*ln(1 + Iph/IS)``, found by bisecting the MODEL's
    current for zero.  Then the real card's ``Voc``, against the
    independent implicit solve."""
    ideal = dict(CELL, rs=0.0, rsh=1e30)
    for sun in (0.2, 1.0):
        iph = CELL['resp'] * sun
        assert_allclose(_cell_current(0.0, sun, ideal), iph, rtol=1e-9)
        voc = CELL['n'] * _UT * math.log(1 + iph / CELL['IS'])
        lo, hi = 0.0, 1.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if _cell_current(mid, sun, ideal) > 0:
                lo = mid
            else:
                hi = mid
        assert_allclose(0.5 * (lo + hi), voc, rtol=1e-6)
    ## the real card, at three biases across the knee and beyond Voc
    for v in (0.0, 0.45, 0.55, 0.62):
        assert_allclose(_cell_current(v, 1.0), _cell_iv(v, 1.0),
                        rtol=1e-8, atol=1e-12)


def test_solar_cell_fill_factor_and_maximum_power():
    """Maximum-power point of the model's I-V against the same point on
    the independent solve, and the fill factor inside the band a
    silicon cell with ``n = 1.3`` and these parasitics has (0.6..0.85
    -- Green's approximation gives 0.80 for the ideal cell here).  The
    band is a plausibility check; the point-to-point comparison is
    the evidence."""
    vs = np.linspace(0.0, 0.72, 73)
    im = np.array([_cell_current(v, 1.0) for v in vs])
    ir = np.array([_cell_iv(v, 1.0) for v in vs])
    assert_allclose(im, ir, rtol=1e-7, atol=1e-10)
    pm = vs * im
    k = int(np.argmax(pm))
    assert 0 < k < len(vs) - 1
    isc = im[0]
    voc = vs[np.nonzero(im <= 0)[0][0]]
    ff = pm[k] / (isc * voc)
    assert 0.6 < ff < 0.85, ff


def test_led_output_is_linear_in_current_above_threshold():
    """``V(opt) = eta*(I - ith)`` above ``ith``, zero below, driven by a
    current source: five currents, the slope by regression equal to
    ``eta`` to 1e-9, and nothing below threshold."""
    eta, ith = 0.1, 5e-3
    got = []
    for i in (1e-3, 4e-3, 6e-3, 10e-3, 20e-3):
        c = SubCircuit()
        c['is'] = ISRC(gnd, 'a', i=i)
        c['X'] = eh.LedHdl('a', gnd, 'opt', gnd, IS=1e-18, n=2.0, eta=eta,
                           ith=ith, tnom=_T0)
        c['R'] = R('opt', gnd, r=1e6)
        got.append(float(_dc(c).v('opt')))
    assert got[0] == 0.0 and got[1] == 0.0
    assert_allclose(got[2:], [eta * (i - ith) for i in (6e-3, 10e-3, 20e-3)],
                    rtol=1e-9)


def test_optocoupler_current_transfer_ratio_is_resp_times_eta():
    """An LED driving a photodiode's optical port: the photodiode's
    short-circuit current is ``resp*eta*I_led``.  The optical branch
    carries no current -- the LED's V-contributed optical branch
    current is zero to abstol, because the photodiode only READS the
    port -- which is the "free in MNA" claim of the roadmap, checked."""
    eta, resp, iled = 0.05, 0.6, 10e-3
    c = SubCircuit()
    c['is'] = ISRC(gnd, 'a', i=iled)
    c['L'] = eh.LedHdl('a', gnd, 'opt', gnd, IS=1e-18, n=2.0, eta=eta,
                       tnom=_T0)
    c['P'] = eh.PhotodiodeHdl('pa', gnd, 'opt', gnd, IS=1e-12, n=1.0,
                              resp=resp, tnom=_T0)
    c['vs'] = VS('pa', gnd, v=0.0)          # short circuit
    r = _dc(c)
    assert_allclose(float(r.i('vs.plus')), resp * eta * iled, rtol=1e-9)
    ## the optical port's own branch current: the LED's V-contributed
    ## optical branch carries an unknown, and nothing draws from it
    k = [str(n) for n in c.nodes].index('opt')
    assert abs(float(r.x[k])) > 0                   # the port IS driven
    assert_allclose(float(r.x[k]), eta * iled, rtol=1e-9)


## ======================================================================
## 4.  MESFET / HEMT
## ======================================================================

def _ids(cls, vd, vg, vs=0.0, **kw):
    c = SubCircuit()
    c['vd'] = VS('d', gnd, v=vd)
    c['vg'] = VS('g', gnd, v=vg)
    c['vs'] = VS('s', gnd, v=vs)
    c['X'] = cls('d', 'g', 's', **kw)
    return -float(_dc(c).i('vd.plus'))


def _curtice(vd, vg, vto=-2.0, beta=2.5e-3, alpha=2.0, lam=0.0):
    return beta * max(vg - vto, 0.0) ** 2 * (1 + lam * vd) * math.tanh(
        alpha * vd)


def _statz(vd, vg, vs, vto=-2.0, beta=2.5e-3, alpha=2.0, lam=0.0, b=0.3):
    """`mesload.c`, forward and reverse arms."""
    vgs, vgd, vds = vg - vs, vg - vd, vd - vs
    if vds >= 0:
        vgst = vgs - vto
        if vgst <= 0:
            return 0.0
        betap = beta * (1 + lam * vds)
        core = betap * vgst ** 2 / (1 + b * vgst)
        if vds >= 3 / alpha:
            return core
        return core * (1 - (1 - alpha * vds / 3) ** 3)
    vgdt = vgd - vto
    if vgdt <= 0:
        return 0.0
    betap = beta * (1 - lam * vds)
    core = betap * vgdt ** 2 / (1 + b * vgdt)
    if -vds >= 3 / alpha:
        return -core
    return -core * (1 - (1 + alpha * vds / 3) ** 3)


def _angelov(vd, vg, ipk=50e-3, vpk=-0.2, p1=1.5, p2=0.0, p3=0.0,
             alphar=1.0, alphas=0.5, lam=0.02):
    dv = vg - vpk
    psi = p1 * dv + p2 * dv ** 2 + p3 * dv ** 3
    tp = 1 + math.tanh(psi)
    return ipk * tp * (1 + lam * vd) * math.tanh((alphar + alphas * tp) * vd)


def test_curtice_is_the_published_equation():
    card = dict(vto=-1.5, beta=3e-3, alpha=2.5)
    card['lambda'] = 0.05
    for vd, vg in ((3.0, -1.0), (0.3, -0.5), (1.0, 0.0), (3.0, -2.0),
                   (0.05, -0.2)):
        got = _ids(eh.MesfetCurticeHdl, vd, vg, **card)
        assert_allclose(got, _curtice(vd, vg, -1.5, 3e-3, 2.5, 0.05),
                        rtol=1e-9, atol=1e-20)


def test_statz_is_mesload_in_both_modes():
    """Nine biases: linear, cubic-saturation, above ``3/alpha``, cutoff,
    and four with the drain below the source so the reverse arm and the
    ``lambda`` sign flip are exercised.  Then the gate diodes: the gate
    current at 0.7 V forward is ``2*IS*(exp(v/Vt) - 1)`` when both
    junctions see it."""
    card = dict(vto=-2.0, beta=2.5e-3, alpha=2.0, b=0.3)
    card['lambda'] = 0.04
    for vd, vg, vs in ((3.0, -1.0, 0.0), (1.0, -1.0, 0.0), (0.2, -0.5, 0.0),
                       (3.0, -2.5, 0.0), (0.0, -1.0, 3.0), (0.0, -1.0, 1.0),
                       (0.0, -0.5, 0.2), (2.0, 0.2, 3.0), (1.5, -1.0, 1.5)):
        got = _ids(eh.MesfetStatzHdl, vd, vg, vs, **card)
        ## the drain terminal also carries the gate-drain Schottky
        ## junction's current (IS = 1e-14 by default), reverse here
        igd = 1e-14 * (math.exp((vg - vd) / _UT) - 1)
        assert_allclose(got, _statz(vd, vg, vs, -2.0, 2.5e-3, 2.0, 0.04, 0.3)
                        - igd, rtol=1e-9, atol=1e-30)
    c = SubCircuit()
    c['vg'] = VS('g', gnd, v=0.7)
    c['X'] = eh.MesfetStatzHdl('d', 'g', 's', IS=1e-14)
    c['vd'] = VS('d', gnd, v=0.0)
    c['vs'] = VS('s', gnd, v=0.0)
    ig = -float(_dc(c).i('vg.plus'))
    assert_allclose(ig, 2 * 1e-14 * (math.exp(0.7 / _UT) - 1), rtol=1e-9)


def test_angelov_is_the_published_equation_and_gm_peaks_at_vpk():
    """Angelov 1992 eqs. (1)-(3) at six biases with every coefficient
    non-zero; then the two identities of the parameterisation with
    ``p2 = p3 = alphas = 0``: the saturated current at ``vgs = vpk`` is ``ipk``
    (``tanh(0) = 0``) and ``gm`` -- measured by central difference on
    DC solves -- is maximal there, at ``ipk*p1*tanh(alpha*vds)``."""
    card = dict(ipk=60e-3, vpk=-0.3, p1=1.2, p2=0.3, p3=-0.1, alphar=0.8,
                alphas=0.6)
    card['lambda'] = 0.03
    for vd, vg in ((3.0, -0.3), (0.2, 0.0), (1.0, -0.8), (3.0, 0.3),
                   (0.05, -0.3), (2.0, -1.5)):
        got = _ids(eh.HemtAngelovHdl, vd, vg, **card)
        assert_allclose(got, _angelov(vd, vg, 60e-3, -0.3, 1.2, 0.3, -0.1,
                                      0.8, 0.6, 0.03), rtol=1e-9)
    ## `alphas = 0` too: with the gate-dependent saturation parameter on,
    ## `tanh(alpha*vds)` moves with `vgs` and the peak is displaced by
    ## 7e-4 (measured) -- the identity is of the p1-only, fixed-alpha
    ## parameterisation
    card = dict(ipk=50e-3, vpk=-0.2, p1=1.5, p2=0.0, p3=0.0, alphar=1.5,
                alphas=0.0)
    card['lambda'] = 0.0
    vd = 3.0
    assert_allclose(_ids(eh.HemtAngelovHdl, vd, -0.2, **card),
                    50e-3 * math.tanh(1.5 * vd), rtol=1e-9)
    ## gm by central differences on the ELEMENT's current (no solver
    ## tolerance in the way; a DC solve's 1e-4 reltol showed up as 7e-4
    ## on the difference when this was first written through `_ids`)
    el = _mk(eh.HemtAngelovHdl, 'd', 'g', 's', **card)
    h = 1e-4

    def ids_el(vg):
        return float(el.i(np.array([vd, vg, 0.0]), defaultepar)[0])
    gm = {vg: (ids_el(vg + h) - ids_el(vg - h)) / (2 * h)
          for vg in (-0.6, -0.4, -0.2, 0.0, 0.2)}
    assert max(gm, key=gm.get) == -0.2, gm
    assert_allclose(gm[-0.2], 50e-3 * 1.5 * math.tanh(1.5 * vd), rtol=1e-7)


def test_gm_over_id_shapes():
    """``gm/Id`` as a function of gate drive, from central differences on
    DC solves: Curtice's square law gives ``2/(vgs - vto)`` (falling as
    the drive grows), Angelov's gives ``p1*(1 - tanh(psi))`` -- equal
    to ``p1`` at ``vpk`` and falling to zero as the channel saturates
    -- both to 1e-4."""
    h = 1e-4
    for vg in (-1.5, -1.0, -0.5):
        i0 = _ids(eh.MesfetCurticeHdl, 3.0, vg)
        gm = (_ids(eh.MesfetCurticeHdl, 3.0, vg + h)
              - _ids(eh.MesfetCurticeHdl, 3.0, vg - h)) / (2 * h)
        assert_allclose(gm / i0, 2.0 / (vg + 2.0), rtol=1e-4)
    for vg in (-0.6, -0.2, 0.4):
        i0 = _ids(eh.HemtAngelovHdl, 3.0, vg, alphas=0.0)
        gm = (_ids(eh.HemtAngelovHdl, 3.0, vg + h, alphas=0.0)
              - _ids(eh.HemtAngelovHdl, 3.0, vg - h, alphas=0.0)) / (2 * h)
        psi = 1.5 * (vg + 0.2)
        assert_allclose(gm / i0, 1.5 * (1 - math.tanh(psi)), rtol=1e-4)


@pytest.mark.parametrize('cls,card,x', [
    (eh.MesfetCurticeHdl, {'lambda': 0.05}, [3.0, -1.0, 0.0]),
    (eh.MesfetCurticeHdl, {}, [0.2, -0.5, 0.0]),
    (eh.MesfetCurticeHdl, {}, [3.0, -2.0, 0.0]),        # at cutoff
    (eh.HemtAngelovHdl, dict(p2=0.3, p3=-0.1), [3.0, -0.3, 0.0]),
    (eh.HemtAngelovHdl, {}, [0.1, 0.2, 0.0]),
    (eh.MesfetStatzHdl, dict(cgs=1e-12, cgd=0.5e-12, rd=5.0, rs=3.0),
     [3.0, -1.0, 0.0, 2.9, 0.1]),
    (eh.MesfetStatzHdl, dict(cgs=1e-12, cgd=0.5e-12), [0.0, -1.0, 3.0]),
    (eh.MesfetStatzHdl, dict(cgs=1e-12, cgd=0.5e-12), [1.0, 0.6, 0.0]),
])
def test_mesfet_jacobians_by_finite_differences(cls, card, x):
    el = _mk(cls, *('d', 'g', 's'), **card)
    res = check_jacobians(el, np.array(x, float), rtol=3e-5)
    assert res.ok, res
    for e in res.unresolved:
        assert e.reason in ('roundoff', 'truncation'), e


class _CurticeNoLimit(Behavioural):
    """`MesfetCurticeHdl` without its limiters, the control."""
    params_as = 'p'
    instparams = [Parameter(name='vto', desc='', unit='V', default=-2.0),
                  Parameter(name='beta', desc='', unit='', default=2.5e-3),
                  Parameter(name='alpha', desc='', unit='', default=2.0)]

    @staticmethod
    def analog(p, d, g, s):
        bgs, bds = Branch(g, s), Branch(d, s)
        vgst = var(eh._maxc(bgs.V - p.vto, 0.0), 'vgst')
        return Contribution(bds.I, p.beta * vgst * vgst
                            * sympy.tanh(p.alpha * bds.V))


class _StatzNoLimit(Behavioural):
    """`MesfetStatzHdl`'s intrinsic channel and gate diodes with NO
    limiter -- the control for the group."""
    params_as = 'p'
    instparams = [Parameter(name='vto', desc='', unit='V', default=-2.0),
                  Parameter(name='beta', desc='', unit='', default=2.5e-3),
                  Parameter(name='alpha', desc='', unit='', default=2.0),
                  Parameter(name='b', desc='', unit='', default=0.3),
                  Parameter(name='IS', desc='', unit='A', default=1e-14)]

    @staticmethod
    def analog(p, d, g, s):
        bgs, bgd, bds = Branch(g, s), Branch(g, d), Branch(d, s)
        vgs, vgd, vds = var(bgs.V, 'vgs'), var(bgd.V, 'vgd'), var(bds.V, 'vds')
        vtT = var(eh._vt(TEMP), 'vtT')

        def channel(vg, vd, tag):
            vgst = var(eh._maxc(vg - p.vto, 0.0), 'vgst' + tag)
            core = var(p.beta * vgst * vgst / (1 + p.b * vgst), 'core' + tag)
            afact = var(1 - p.alpha * eh._minc(vd, 3.0 / p.alpha) / 3,
                        'afact' + tag)
            return var(core * (1 - afact ** 3), 'ids' + tag)
        ids = var(sympy.Piecewise((channel(vgs, vds, 'f'), vds >= 0.0),
                                  (-channel(vgd, -vds, 'r'), True)), 'ids')
        return (Contribution(bds.I, ids),
                Contribution(bgs.I, p.IS * (eh._expl(vgs / vtT) - 1)),
                Contribution(bgd.I, p.IS * (eh._expl(vgd / vtT) - 1)))


def _cs_stage(cls, vdd=5.0, vg=-1.0, rl=500.0, **kw):
    """A common-source stage: ``rl`` from ``vdd`` to the drain, the gate
    driven through 1 kOhm, so a wild iterate can forward-bias the gate."""
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=vdd)
    c['vg'] = VS('vgg', gnd, v=vg)
    c['rg'] = R('vgg', 'g', r=1e3)
    c['rl'] = R('vdd', 'd', r=rl)
    c['X'] = cls('d', 'g', gnd, **kw)
    return c


def _newton_from(c, x0, maxiter=100):
    from pycircuit.circuit.tests.test_limit_identity import _plain_newton_from
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            x, its = _plain_newton_from(c, x0, maxiter)
        return its, x
    except (NoConvergenceError, SingularMatrix, FloatingPointError,
            np.linalg.LinAlgError, OverflowError, ValueError):
        return None, None


def test_mesfet_limiters_measured_against_their_absence():
    """`newton-limiting`: the intervention against its absence, on a real
    solve, from wild starts, counting Jacobian evaluations, with the
    answer checked against ``DC()``.

    Recorded at introduction, common-source stage, starts uniform at
    0 / +20 / -20 V, ``[limited, unlimited]``::

        Curtice   0 V: [5, 5]     +20 V: [4, 3]     -20 V: [6, 3]
        Statz     0 V: [5, 4]     +20 V: [14, F]    -20 V: [6, 3]

    The tanh channel cannot overflow, so on Curtice the
    ``fetlim``/``limvds`` group buys nothing and COSTS one to three
    iterations from a wild start (as `limit_together`'s first user
    found on level 1: no exponential, nothing to compress -- and here
    the clamp holds a bounded device back).  Statz has two Schottky
    junctions: from +20 V the unlimited device overflows them and
    fails, and the ``pnjlim`` group brings it home in 14; from -20 V
    (junctions reverse) the unlimited device is faster.  Both halves
    asserted: the limited device converges from every start to
    ``DC()``'s drain voltage, and the unlimited Statz does NOT from
    +20 V -- so a limiter that stopped biting would be visible.
    """
    ref = {}
    for cls, twin, card in ((eh.MesfetCurticeHdl, _CurticeNoLimit, {}),
                            (eh.MesfetStatzHdl, _StatzNoLimit, {})):
        c = _cs_stage(cls, **card)
        ref[cls] = float(_dc(c).v('d'))
        for start in (0.0, 20.0, -20.0):
            c = _cs_stage(cls, **card)
            n = len(c.nodes) + len(c.branches)
            x0 = np.full(n, start)
            x0[c.get_node_index(gnd)] = 0.0
            its, x = _newton_from(c, x0)
            assert its is not None and its <= 30, (cls.__name__, start, its)
            assert_allclose(x[c.get_node_index('d')], ref[cls], rtol=1e-5)
            c2 = _cs_stage(twin, **card)
            its2, x2 = _newton_from(c2, x0)
            if its2 is not None:
                assert_allclose(x2[c2.get_node_index('d')], ref[cls],
                                rtol=1e-4)
            if cls is eh.MesfetStatzHdl and start == 20.0:
                assert its2 is None, its2


def test_mesfet_limiter_declarations_are_what_the_docstrings_say():
    txt = explain(eh.MesfetStatzHdl, source=False, symbolic=False)
    assert ('3 $limit probes (pnjlim on (g,si) [group 1], pnjlim on (g,di) '
            '[group 1], limvds on (di,si) [group 1])') in txt, txt
    txt = explain(eh.HemtAngelovHdl, source=False, symbolic=False)
    assert ('2 $limit probes (fetlim on (g,s) [group 1], limvds on (d,s) '
            '[group 1])') in txt, txt
    assert 'PCNR: vector route, 2 probes over 3 rows' in txt, txt
    ## and `lambda` is the parameter's real name on every one of them
    for cls in (eh.MesfetStatzHdl, eh.MesfetCurticeHdl, eh.HemtAngelovHdl):
        assert 'lambda' in [q.name for q in cls.instparams]


## ======================================================================
## 5.  MOS level 3
## ======================================================================

#: The coincidence card: no body effect, every level-3 knob off.  Here
#: level 3's bulk-charge triode term has ``fbody = 0`` and its ``vdsat``
#: is ``vgs - vth``, which is Shichman-Hodges exactly.
COINCIDE = dict(vto=0.75, kp=8e-5, gamma=0.0, phi=0.70, w=20e-6, l=2e-6,
                tnom=_T0)

#: Every knob on, at values a 1 um process card would carry.  ``as`` is
#: SPICE's own spelling.
FULL = dict(vto=0.7, u0=500.0, tox=2e-8, nsub=1e16, xj=3e-7, nfs=1e11,
            eta=0.05, delta=0.5, theta=0.1, vmax=1e5, kappa=0.3,
            w=10e-6, l=1e-6, ld=0.05e-6, tnom=_T0,
            IS=1e-14, cj=3e-4, cjsw=3e-10, ad=1.6e-11, pd=1.2e-5,
            ps=1.2e-5, cgso=3.5e-10, cgdo=2.5e-10, cgbo=2.0e-10)
FULL['as'] = 1.6e-11


def _m3_ref(vd, vg, vs, vb, card):
    """SPICE level 3, transcribed from `mos3load.c` (via `mos3.va`),
    written as the C is: ``if`` on the solution, the mode swap by
    comparing the drain and source potentials, ``goto`` labels as
    early returns.  At ``T = tnom`` the temperature path is the
    identity, so ``vbi = vto - gamma*sqrt(phi)``.  Returns the channel
    current in the n-channel convention plus the two bulk currents."""
    p = dict(card)
    leff = p['l'] - 2 * p.get('ld', 0.0)
    weff = p['w'] - 2 * p.get('wd', 0.0)
    cox = _EPSOX / p['tox'] if 'tox' in p else _EPSOX / 1e-7
    nsub = p.get('nsub', 0.0) * 1e6
    if 'gamma' in p:
        gamma = p['gamma']
    else:
        gamma = math.sqrt(2 * _QE * _EPSSI * nsub) / cox
    if 'phi' in p:
        phi = p['phi']
    else:
        phi = max(2 * _UT * math.log(nsub / (_NI * 1e6)), 0.1)
    kp = p.get('kp', 0.0) or p.get('u0', 600.0) * cox * 1e-4
    u0 = p.get('u0', 600.0)
    vbi = p['vto'] - gamma * math.sqrt(phi)
    alpha = 2 * _EPSSI / (_QE * nsub) if nsub > 0 else 0.0
    xd = math.sqrt(alpha)
    narrow = p.get('delta', 0.0) * 0.5 * math.pi * _EPSSI / cox
    etal = p.get('eta', 0.0) * 8.15e-22 / (cox * leff ** 3)
    xj, ld, nfs = p.get('xj', 0.0), p.get('ld', 0.0), p.get('nfs', 0.0)
    theta, vmax, kappa = p.get('theta', 0.0), p.get('vmax', 0.0), \
        p.get('kappa', 0.2)
    if vd >= vs:
        mode, vgs, vds, vbs = 1.0, vg - vs, vd - vs, vb - vs
    else:
        mode, vgs, vds, vbs = -1.0, vg - vd, vs - vd, vb - vd

    def channel():
        if vbs <= 0:
            sqphbs = math.sqrt(phi - vbs)
            phibs = phi - vbs
        else:
            sqphbs = math.sqrt(phi) / (1 + vbs / (2 * phi))
            phibs = sqphbs ** 2
        if xj != 0 and xd != 0:
            wps = xd * sqphbs
            wponxj = wps / xj
            wconxj = 0.0631353 + 0.8013292 * wponxj - 0.01110777 * wponxj ** 2
            arga = wconxj + ld / xj
            argc = wponxj / (1 + wponxj)
            argb = math.sqrt(1 - argc * argc)
            fshort = 1 - (xj / leff) * (arga * argb - ld / xj)
        else:
            fshort = 1.0
        gammas = gamma * fshort
        fbody = 0.5 * gammas / (2 * sqphbs) + narrow / weff
        qbonco = gammas * sqphbs + narrow * phibs / weff
        vth = vbi - etal * vds + qbonco
        von = vth
        xn = 0.0
        if nfs != 0:
            xn = 1 + _QE * nfs * 1e4 / cox + qbonco / (2 * phibs)
            von = vth + _UT * xn
        elif vgs <= von:
            return 0.0
        vgsx = max(vgs, von)
        onfg = 1 + theta * (vgsx - vth)
        fgate = 1 / onfg
        us = u0 * 1e-4 * fgate
        vdsat = (vgsx - vth) / (1 + fbody)
        vdsc = None
        if vmax > 0:
            vdsc = leff * vmax / us
            vdsat = vdsat + vdsc - math.sqrt(vdsat ** 2 + vdsc ** 2)
        vdsx = min(vds, vdsat)
        if vdsx == 0:
            return 0.0
        cdo = vgsx - vth - 0.5 * (1 + fbody) * vdsx
        beta = kp * weff / leff * fgate
        cdrain = beta * cdo * vdsx
        fdrain = 1.0
        if vmax > 0:
            fdrain = 1 / (1 + vdsx / vdsc)
            cdrain *= fdrain
        delxl = None
        if vds <= vdsat:
            if not (vmax > 0 or alpha == 0):
                delxl = math.sqrt(kappa * alpha * vdsat / 8) * (vds / vdsat) ** 4
        elif vmax <= 0:
            delxl = math.sqrt(kappa * alpha * (vds - vdsat + vdsat / 8))
        elif alpha != 0:
            cdsat = cdrain
            gdsat = max(cdsat * (1 - fdrain) / vdsc, 1e-12)
            emax = kappa * cdsat / (leff * gdsat)
            arga = 0.5 * emax * alpha
            argc = kappa * alpha
            argb = math.sqrt(arga * arga + argc * (vds - vdsat))
            delxl = argb - arga
        if delxl is not None:
            if delxl > 0.5 * leff:
                delxl = leff - leff * leff / (4 * delxl)
            cdrain = cdrain / (1 - delxl / leff)
        if vgs < von:
            cdrain *= math.exp((vgs - von) / (xn * _UT))
        return cdrain
    ich = mode * channel()
    is_ = p.get('IS', 1e-14)
    ibs = is_ * (math.exp(min((vb - vs) / _UT, 300.0)) - 1)
    ibd = is_ * (math.exp(min((vb - vd) / _UT, 300.0)) - 1)
    return ich, ibs, ibd


def _m3_terminal_currents(vd, vg, vs, vb, card):
    ich, ibs, ibd = _m3_ref(vd, vg, vs, vb, card)
    return np.array([ich - ibd, 0.0, -ich - ibs, ibs + ibd])


M3_BIAS = [(3.0, 2.0, 0.0, 0.0), (0.3, 2.0, 0.0, 0.0), (1.2, 2.0, 0.0, 0.0),
           (3.0, 0.5, 0.0, 0.0), (3.0, 2.0, 0.0, -2.0), (0.3, 2.0, 0.0, -2.0),
           (0.0, 2.0, 3.0, -2.0), (0.0, 2.0, 0.3, 0.0), (5.0, 5.0, 0.0, -1.0),
           (0.5, 2.0, 0.0, 0.55), (2.0, 2.5, 0.0, 0.65), (3.0, 0.9, 0.0, 0.0),
           (0.05, 1.5, 0.0, 0.0), (0.0, 0.9, 3.0, -0.3)]


def test_level3_coincides_with_level1_where_it_should():
    """At ``gamma = 0`` with every level-3 knob off, the two models are
    the same equation.  Nine biases, both signs of ``vds``, to 1e-9 on
    the channel current -- and the CONTROL: with ``gamma = 0.5`` they
    differ by the bulk-charge factor, which is asserted to be at least
    5% in saturation so the coincidence is not the trivial one."""
    m1 = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **COINCIDE)
    m3 = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **dict(COINCIDE, kappa=0.0))
    for vd, vg, vs, vb in M3_BIAS[:9]:
        x = np.array([vd, vg, vs, vb])
        i1 = m1.i(x, defaultepar)
        i3 = m3.i(x, defaultepar)
        assert_allclose(i3, i1, rtol=1e-9, atol=1e-13)
    m1 = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **dict(COINCIDE, gamma=0.5))
    m3 = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b',
             **dict(COINCIDE, gamma=0.5, kappa=0.0))
    x = np.array([3.0, 2.0, 0.0, 0.0])
    i1, i3 = m1.i(x, defaultepar)[0], m3.i(x, defaultepar)[0]
    assert abs(i3 - i1) > 0.05 * i1, (i1, i3)


def test_level3_matches_the_mos3load_transcription_with_every_knob_on():
    """Fourteen biases on the full card against `_m3_ref`: linear,
    saturation, subthreshold (``nfs`` on), cutoff, both signs of
    ``vds``, forward bulk (the continuation arm), and ``vds`` right at
    the ``vdsat`` seam neighbourhood.  1e-9 relative on all four
    terminal currents."""
    m3 = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **FULL)
    for vd, vg, vs, vb in M3_BIAS:
        x = np.array([vd, vg, vs, vb])
        got = m3.i(x, defaultepar)
        ref = _m3_terminal_currents(vd, vg, vs, vb, FULL)
        assert_allclose(got, ref, rtol=1e-9, atol=1e-18,
                        err_msg=str((vd, vg, vs, vb)))


@pytest.mark.parametrize('knob', ['xj', 'nfs', 'eta', 'delta', 'theta',
                                  'vmax', 'kappa'])
def test_each_level3_knob_does_something_and_the_transcription_agrees(knob):
    """Each knob switched off in turn: the model still agrees with the
    transcription (so the ``Piecewise`` on the parameter selects the
    right arm), and the current CHANGES by at least 0.5% at a bias the
    knob reaches -- a knob wired to nothing would pass the first half
    alone."""
    off = dict(FULL)
    off[knob] = 0.0
    bias = {'nfs': (3.0, 0.3, 0.0, 0.0), 'kappa': (3.0, 2.0, 0.0, 0.0),
            'eta': (5.0, 1.2, 0.0, 0.0)}.get(knob, (3.0, 2.0, 0.0, 0.0))
    x = np.array(bias)
    m_on = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **FULL)
    m_off = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **off)
    i_on = m_on.i(x, defaultepar)[0]
    i_off = m_off.i(x, defaultepar)[0]
    assert_allclose(i_off, _m3_terminal_currents(*bias, off)[0], rtol=1e-9)
    assert abs(i_on - i_off) > 5e-3 * abs(i_off), (knob, i_on, i_off)


def test_level3_subthreshold_slope_is_xn_times_vt():
    """With ``nfs`` on, the current below ``von`` is exponential with
    slope ``ln(10)*xn*Vt`` per decade; measured on the model between
    two gate voltages 0.2 V apart in weak inversion and compared with
    the ``xn`` the transcription computes.  ``eta`` off so ``von``
    does not move with ``vds``, and ``IS`` off: the drain terminal
    carries the reverse bulk-drain junction's 1e-14 A too, which at
    1e-12 A of channel current bends the log slope by 3e-4."""
    card = dict(FULL, eta=0.0, IS=1e-30)
    m3 = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **card)
    i1 = m3.i(np.array([1.0, 0.3, 0.0, 0.0]), defaultepar)[0]
    i2 = m3.i(np.array([1.0, 0.5, 0.0, 0.0]), defaultepar)[0]
    assert 0 < i1 < i2 < 1e-6
    slope = 0.2 / math.log10(i2 / i1)          # V/decade
    ## xn from the transcription's ingredients at vbs = 0
    cox = _EPSOX / card['tox']
    nsub = card['nsub'] * 1e6
    gamma = math.sqrt(2 * _QE * _EPSSI * nsub) / cox
    phi = 2 * _UT * math.log(nsub / (_NI * 1e6))
    leff = card['l'] - 2 * card['ld']
    xd = math.sqrt(2 * _EPSSI / (_QE * nsub))
    wponxj = xd * math.sqrt(phi) / card['xj']
    wconxj = 0.0631353 + 0.8013292 * wponxj - 0.01110777 * wponxj ** 2
    argc = wponxj / (1 + wponxj)
    fshort = 1 - (card['xj'] / leff) * (
        (wconxj + card['ld'] / card['xj']) * math.sqrt(1 - argc ** 2)
        - card['ld'] / card['xj'])
    narrow = card['delta'] * 0.5 * math.pi * _EPSSI / cox
    qbonco = gamma * fshort * math.sqrt(phi) + narrow * phi / card['w']
    xn = 1 + _QE * card['nfs'] * 1e4 / cox + qbonco / (2 * phi)
    assert_allclose(slope, math.log(10) * xn * _UT, rtol=1e-6)


def test_level3_pmos_is_the_reversed_nmos():
    """``MosLevel3PmosHdl`` at ``(vd, vg, vs, vb)`` equals the n-channel
    at ``(-vd, -vg, -vs, -vb)`` with every current negated."""
    n = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **FULL)
    p = _mk(eh.MosLevel3PmosHdl, 'd', 'g', 's', 'b', **FULL)
    for vd, vg, vs, vb in M3_BIAS:
        x = np.array([vd, vg, vs, vb])
        assert_allclose(p.i(-x, defaultepar), -n.i(x, defaultepar),
                        rtol=1e-12, atol=1e-25)


@pytest.mark.parametrize('name,card,x', [
    ('sat', FULL, [3.0, 2.0, 0.0, -0.5]),
    ('lin', FULL, [0.3, 2.0, 0.0, 0.0]),
    ('sub', FULL, [3.0, 0.5, 0.0, 0.0]),
    ('rev', FULL, [0.0, 2.0, 3.0, -1.0]),
    ('fwd-bulk', FULL, [1.5, 2.0, 0.0, 0.3]),
    ('no-vmax', dict(FULL, vmax=0.0), [3.0, 2.0, 0.0, -0.5]),
    ('no-vmax-lin', dict(FULL, vmax=0.0), [0.3, 2.0, 0.0, 0.0]),
    ('no-nsub', dict(FULL, nsub=0.0, xj=0.0), [3.0, 2.0, 0.0, -0.5]),
    ('no-nfs-cut', dict(FULL, nfs=0.0), [3.0, 0.5, 0.0, 0.0]),
    ('rsh', dict(FULL, rsh=20.0, nrd=2.0, nrs=2.0),
     [3.0, 2.0, 0.0, -0.5, 2.9, 0.05]),
])
def test_level3_jacobians_by_finite_differences(name, card, x):
    el = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **card)
    res = check_jacobians(el, np.array(x, float), rtol=3e-5)
    assert res.ok, '%s\n%s' % (name, res)
    assert not res.nonfinite, (name, res.nonfinite)
    for e in res.unresolved:
        assert e.reason in ('roundoff', 'truncation'), (name, e)
        if e.reason == 'roundoff':
            assert abs(e.ana) < e.floor, (name, e)


def test_level3_jacobian_check_fails_on_a_corrupted_derivative():
    el = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **FULL)
    x = np.array([3.0, 2.0, 0.0, -0.5])
    assert check_jacobians(el, x, rtol=3e-5).ok
    good = el.G
    try:
        el.G = lambda *a, **kw: 3.0 * good(*a, **kw)
        bad = check_jacobians(el, x, rtol=3e-5)
    finally:
        del el.G
    assert not bad.ok and bad.verdict == 'FAILED', bad


def test_level3_as_is_declared_under_spices_own_name():
    """``as`` is a Python keyword and the first parameter in the library
    declared under it.  It is accepted on the instance, read as
    ``p['as']``, and it does something: with ``js`` given the source
    junction's saturation current is ``js*as``."""
    assert 'as' in [q.name for q in eh.MosLevel3Hdl.instparams]
    card = dict(FULL, js=1e-4)
    a = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **dict(card, **{'as': 1e-10}))
    b = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **dict(card, **{'as': 2e-10}))
    x = np.array([0.0, 0.0, 0.0, 0.5])          # source junction forward
    ia = a.i(x, defaultepar)[2]
    ib = b.i(x, defaultepar)[2]
    assert_allclose(ib / ia, 2.0, rtol=1e-9)


def test_level3_ohmic_resistance_collapses_or_makes_internal_nodes():
    bare = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **FULL)
    withr = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **dict(FULL, rd=10.0,
                                                          rs=5.0))
    assert len(x_layout(bare)) == 4
    assert len(x_layout(withr)) == 6
    ## and the resistance is where it says: a 1 mA drain current through
    ## rd = 10 ohm drops 10 mV between d and di
    x = np.array([3.0, 2.0, 0.0, 0.0, 2.99, 0.0])
    i = withr.i(x, defaultepar)
    assert_allclose(i[0], 0.01 / 10.0, rtol=1e-12)


def test_level3_inverter_switches_where_a_dc_sweep_says():
    """A CMOS inverter from the two level-3 devices: output high at
    ``vin = 0``, low at ``vin = vdd``, and the switching point found by
    bisection lies between the two thresholds -- a solve through the
    whole DC chain on a circuit with both models, both modes and
    the bulk junctions live."""
    nm = dict(FULL)
    pm = dict(FULL, w=20e-6)

    def vout(vin):
        c = SubCircuit()
        c['vdd'] = VS('vdd', gnd, v=3.0)
        c['vin'] = VS('in', gnd, v=vin)
        c['MN'] = eh.MosLevel3Hdl('out', 'in', gnd, gnd, **nm)
        c['MP'] = eh.MosLevel3PmosHdl('out', 'in', 'vdd', 'vdd', **pm)
        return float(_dc(c).v('out'))
    assert vout(0.0) > 2.99 and vout(3.0) < 0.01
    lo, hi = 0.0, 3.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if vout(mid) > 1.5:
            lo = mid
        else:
            hi = mid
    assert 1.0 < 0.5 * (lo + hi) < 2.0, (lo, hi)


## ======================================================================
## 6.  The batch's own claims: params_as everywhere, and the cache
## ======================================================================

NEW = [eh.VcoHdl, eh.DividerHdl, eh.MixerHdl, eh.ChargePumpHdl,
       eh.PhotodiodeHdl, eh.LedHdl, eh.MesfetStatzHdl, eh.MesfetCurticeHdl,
       eh.HemtAngelovHdl, eh.MosLevel3Hdl, eh.MosLevel3PmosHdl]


def test_every_model_in_this_batch_uses_the_parameter_namespace():
    for cls in NEW:
        assert cls.params_as == 'p', cls.__name__
        assert cls._hdl_cache_status in ('hit', 'miss', 'disabled'), (
            cls.__name__, cls._hdl_cache_status)


def test_level3_no_vmax_card_matches_the_transcription_in_triode():
    """Mutation-driven: the below-``vdsat`` channel-length-modulation
    ramp ``sqrt(kappa*alpha*vdsat/8)*(vds/vdsat)^4`` exists only on a
    ``vmax = 0`` card, and the full-card comparison never selects it.
    Changing its exponent from 4 to 3 survived every test above --
    measured -- so here is the card and the biases that reach it, plus
    the plain-root arm above ``vdsat`` on the same card."""
    card = dict(FULL, vmax=0.0)
    m3 = _mk(eh.MosLevel3Hdl, 'd', 'g', 's', 'b', **card)
    for vd, vg, vs, vb in ((0.3, 2.0, 0.0, 0.0), (0.8, 2.0, 0.0, 0.0),
                           (1.2, 2.0, 0.0, -1.0), (3.0, 2.0, 0.0, 0.0),
                           (0.0, 2.0, 0.6, 0.0)):
        x = np.array([vd, vg, vs, vb])
        assert_allclose(m3.i(x, defaultepar),
                        _m3_terminal_currents(vd, vg, vs, vb, card),
                        rtol=1e-9, atol=1e-18, err_msg=str((vd, vg, vs, vb)))
