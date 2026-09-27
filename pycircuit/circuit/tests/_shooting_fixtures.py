"""Circuits, noise elements and helpers shared by the shooting test
files (`test_shooting_*.py`).  Split out of test_analysis_shooting.py
on 2026-09-27; a helper used by one file lives in that file.
"""
from pycircuit.circuit import *
from pycircuit.circuit.shooting import (PAC, algebraic_conditioning,
                                        topological_index)
import warnings
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter as _HdlParameter, white_noise)
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest
import functools as _functools
from pycircuit.circuit.tests._shooting_elements import (_NuModNoise,
    _NuMult,
    _SgnAmpFlicker,
    _SgnPsdFlicker,
    _SwitchHdl)


## ---------------------------------------------------------------------------
## STAGE 11 -- PSS: `method` now selects something, and the inverse is a solve.
## ---------------------------------------------------------------------------

def _series_rlc(Lv=1e-3, Cv=1e-9, Rv=50.0, va=1.0):
    """Series RLC driven AT resonance, where |v(C)| = Q * va analytically."""
    import numpy as _np
    circuit.default_toolkit = circuit.numeric
    f0 = 1.0 / (2 * _np.pi * _np.sqrt(Lv * Cv))
    c = SubCircuit()
    c['vs'] = VSin(1, gnd, va=va, freq=f0)
    c['R'] = R(1, 2, r=Rv)
    c['L'] = L(2, 3, L=Lv)
    c['C'] = C(3, gnd, c=Cv)
    return c, f0, (1.0 / Rv) * _np.sqrt(Lv / Cv)


def _pss_peak(method, steps=20):
    import warnings
    cir, f0, Q = _series_rlc()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = PSS(cir, method=method).solve(period=1.0 / f0,
                                            timestep=1.0 / (f0 * steps))
    v = np.asarray(res['tpss'].v(3, gnd), dtype=float)
    return 0.5 * (v.max() - v.min()), Q


# ---------------------------------------------------------------------------
# Phase 1: the shooting Newton was not a Newton
# ---------------------------------------------------------------------------

def _q20_rlc(f0=1e3, Q=20.0):
    """A resonator whose per-period decay is exp(-pi/Q) = 0.8546.

    That number is the whole diagnostic: successive substitution converges at
    exactly the circuit's own decay rate, so observing 0.855 per iteration is
    how the missing Jacobian was found.
    """
    L_, C_ = 1e-3, 1.0 / ((2 * np.pi * f0) ** 2 * 1e-3)
    c = SubCircuit()
    c.add_node('a'); c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, freq=f0)
    c['R'] = R('a', 'b', r=(1.0 / Q) * np.sqrt(L_ / C_))
    c['L'] = L('b', 'c', L=L_)
    c['C'] = C('c', gnd, c=C_)
    return c


def _shooting_trace(method, reltol=1e-4, maxiterations=30):
    """Run PSS, returning (the shooting Newton's trace as `(max |F|,
    max |I - J|)` per residual evaluation, non-converged?, result)."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_q20_rlc(), method=method, reltol=reltol)
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter('always')
        res = pss.solve(period=1e-3, timestep=1e-5,
                        maxiterations=maxiterations, trace=True)
    nonconv = any('did not converge' in str(c.message) for c in caught)
    trace = [(float(np.max(np.abs(F))),
              float(np.max(np.abs(np.eye(len(z)) - J))))
             for z, F, J in pss.shooting_trace]
    return trace, nonconv, res


def _pss_lte(method, timestep=1e-5, reltol=1e-3, **kw):
    """Run the Q=20 resonator and return (peak amplitude, the pss object)."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_q20_rlc(), method=method, reltol=reltol, **kw)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        res = pss.solve(period=1e-3, timestep=timestep, maxiterations=40)
    pss._caught = [str(x.message) for x in caught
                   if issubclass(x.category, RuntimeWarning)]
    assert pss.converged, '%s did not converge' % method
    peak = float(np.max(np.abs(
        np.asarray(res['tpss'].v('c'), dtype=float).ravel())))
    return peak, pss


def _force_plain_map(pss):
    """Run `pss` on the PLAIN map where it would solve for gear's pair --
    the OLD formulation, one entering unknown with `x_{-1}` manufactured by
    an order-dropped step -- and leave every other kind alone.  Overrides
    `PSS._map_kind`, the one place the map's kind is decided, so everything
    that asks (`_solves_history`, the Newton, `factored_period`) sees the
    same answer.  A test-level override on purpose, not a Parameter: a user
    has no reason to ask for the formulation that measured 1.266e-01 V of
    avoidable error."""
    kind = pss._map_kind()
    pss._map_kind = lambda: 'plain' if kind == 'pair' else kind
    return pss


def _rc_ladder(sections):
    """A driven RC ladder whose `m` grows linearly with `sections`."""
    c = SubCircuit()
    c.add_node('n0')
    c['vs'] = VSin('n0', gnd, va=1.0, freq=1e3)
    for k in range(sections):
        a, b = 'n%d' % k, 'n%d' % (k + 1)
        c.add_node(b)
        c['R%d' % k] = R(a, b, r=1e3)
        c['C%d' % k] = C(b, gnd, c=1e-9)
    return c


def _varying_c_ladder(sections=6):
    """The ladder with a STATE-DEPENDENT capacitance on its last node.

    ⚠ WITHOUT THIS THE MATVEC TESTS CANNOT FAIL, and that was measured, not
    supposed.  `_step_sensitivity` carries a RING of capacitances -- `Cs[0]`
    and `Cs[1]`, the two steps a Gear-2 companion reaches back -- and on a
    linear RC ladder every `C` along the period is the SAME matrix, so
    corrupting the ring (`Cs = [C_new, C_new]` instead of
    `[C_new, Cs[0]]`) changes nothing and the tests stayed green through
    the mutation.  A `q_func` makes `C` a function of the solution, the
    ring entries genuinely differ (measured: 2.99e-09 against a 1e-09
    linear part), and the same mutation is caught.
    """
    c = _rc_ladder(sections)
    last = 'n%d' % sections
    c['Q'] = BSource(last, gnd, last, gnd,
                     q_func=lambda v: 2e-9 * np.tanh(v / 0.5))
    return c


def _scaled_vdp(mu=1.0, gain=1e4):
    """van der Pol with a VCVS copy of `v` scaled by `gain`.

    The copy is perfectly slaved, so the ORBIT is unchanged and only the
    arithmetic inside the phase-pin's `argmax` moves.  That makes it the
    instrument for asking whether comparing coordinates in different units
    is a defect.
    """
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    c.add_node('big')
    c['E'] = VCVS('v', gnd, 'big', gnd, g=gain)
    c['Rb'] = R('big', gnd, r=1e9)
    return c


def _pac_circuit(f0=1e3, Q=20.0):
    """The Q=20 resonator with an AC amplitude on its source.

    `vac` and `va` are different knobs: `va` drives the LARGE signal that
    the periodic operating point is a response to, `vac` is the small
    signal PAC linearises for. A circuit with only `va` set has nothing
    for PAC to analyse, which `PAC.solve` refuses rather than returning
    zeros.
    """
    L_, C_ = 1e-3, 1.0 / ((2 * np.pi * f0) ** 2 * 1e-3)
    c = SubCircuit()
    c.add_node('a'); c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, vac=1.0, freq=f0)
    c['R'] = R('a', 'b', r=(1.0 / Q) * np.sqrt(L_ / C_))
    c['L'] = L('b', 'c', L=L_)
    c['C'] = C('c', gnd, c=C_)
    return c


def _adjoint_ladder(sections=8):
    c = SubCircuit()
    c.add_node('n0')
    c['vs'] = VSin('n0', gnd, va=1.0, vac=1.0, freq=1e3)
    for k in range(sections):
        a, b = 'n%d' % k, 'n%d' % (k + 1)
        c.add_node(b)
        c['R%d' % k] = R(a, b, r=1e3)
        c['C%d' % k] = C(b, gnd, c=1e-7)
    return c


def _diode_mixer():
    c = SubCircuit()
    c['vs'] = VSin(1, gnd, vac=1.0, va=2.0, freq=1e6, phase=20)
    c['R'] = R(1, 2, r=1e4)
    c['D'] = Diode(2, gnd)
    c['C'] = C(2, gnd, c=1e-12)
    return c


def _vdp_with_slow_node(tau_over_T=None, T=6.6634, mu=1.0):
    """van der Pol, optionally with one weakly coupled slow RC node.

    The coupling resistor is large, so the node barely loads the orbit —
    what it adds is a Floquet multiplier at `exp(-T/tau)`, which is the
    only thing under test.
    """
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    if tau_over_T is not None:
        c.add_node('w')
        rbig = 1e6
        c['Rs'] = R('v', 'w', r=rbig)
        c['Cs'] = C('w', gnd, c=tau_over_T * T / rbig)
    return c


def _solve_slow(tau_over_T):
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _vdp_with_slow_node(tau_over_T)
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-11)
    x0 = np.zeros(m)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=6.6634, timestep=6.6634 / 300, x0=x0,
                  maxiterations=60)
    assert pss.converged, 'tau/T=%r did not converge' % (tau_over_T,)
    return cir, pss


def _vdp_with_noise(psd=1e-6, mu=1.0):
    """van der Pol with a stationary white current source at its core.

    `i=0` so the orbit is untouched; only `CY` changes. This is the exact
    configuration the Monte Carlo gate measured, which is what lets the
    closed form be checked against a physical number.
    """
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    c['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
    return c


def _solve_vdp_noise(npts=240, psd=1e-6):
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _vdp_with_noise(psd)
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-12)
    x0 = np.zeros(m)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=6.6634, timestep=6.6634 / npts, x0=x0,
                  maxiterations=60)
    assert pss.converged
    return cir, pss, PAC(cir, toolkit=circuit.numeric)


def _rc_noisy(Cval=1e-7, Rval=1e3, per=1e-3):
    """An RC lowpass with a noisy resistor, driven so a PSS exists.

    Linear, so its linearisation is time-invariant and the covariance is
    constant — and its capacitor-voltage variance is the exact `kT/C`,
    famously independent of `R`.
    """
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, vac=1.0, freq=1.0 / per)
    c['R'] = R('a', 'b', r=Rval)
    c['C'] = C('b', gnd, c=Cval)
    return c


class _Flicker(IS):
    """A current source whose PSD is `noisePSD * (fref / f)` — coloured.

    ⚠ EVERY SOURCE IN THE DISCRETE LIBRARY IS WHITE, so nothing in the
    tree can exercise the coloured path. This exists for that, and for
    nothing else. `CY(x, w)` already takes `w`; no element used it.
    """

    instparams = IS.instparams + [
        Parameter(name='fref', desc='Corner of the 1/f law', unit='Hz',
                  default=1.0)]

    def CY(self, x, w, epar=None):
        f = abs(float(w)) / (2.0 * np.pi)
        scale = self.iparv.fref / max(f, 1e-300)
        p = self.iparv.noisePSD * scale
        return self.toolkit.array([[p, -p], [-p, p]])


def _lc_osc(a=0.0, rs=0.0, psd=1e-6, npts=240, flicker=False, fref=1.0,
            white=0.0):
    """An LC oscillator with optional even nonlinearity and tank loss.

    `a` breaks the waveform's half-wave symmetry; `rs` breaks the LOSSLESS
    tank's structural identity. Both are needed — see the test below.
    `white` adds a white source beside the flicker one.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: 1.0 * (u - u ** 3 / 3.0)
                       + a * (u ** 2 - 2.0))
    if rs > 0.0:
        cir.add_node('x')
        cir['L'] = L('v', 'x', L=1.0)
        cir['Rs'] = R('x', gnd, r=rs)
    else:
        cir['L'] = L('v', gnd, L=1.0)
    if flicker:
        cir['n'] = _Flicker('v', gnd, i=0.0, noisePSD=psd, fref=fref)
    else:
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
    if white > 0.0:
        cir['nw'] = IS('v', gnd, i=0.0, noisePSD=white)
    pss = PSS(cir, method='gear', reltol=1e-12)
    m = cir.n - 1
    x0 = np.zeros(m)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=6.66, timestep=6.66 / npts, x0=x0, maxiterations=80)
    assert pss.converged, 'a=%r rs=%r did not converge' % (a, rs)
    return cir, pss, PAC(cir, toolkit=circuit.numeric)


class _SquareMixer(Circuit):
    """`i_out = g·(v_lo)²·(v_in)` — transconductance `g·cos²(ω₀t)`.

    A four-quadrant-style element built for one purpose: give the
    noise-to-output transfer an exactly known periodic shape. Driving `lo`
    with a unit cosine makes the small-signal transconductance `g·cos²`,
    whose Fourier coefficients are `H₀ = ½`, `H_±2 = ¼` and nothing else.
    """

    terminals = ('outp', 'outn', 'inp', 'inn', 'lop', 'lon')
    instparams = [Parameter(name='g', desc='Transconductance coefficient',
                            unit='A/V^3', default=1.0)]

    def i(self, x, epar=None):
        g = self.iparv.g
        f = g * (x[4] - x[5]) ** 2 * (x[2] - x[3])
        return self.toolkit.array([f, -f, 0.0, 0.0, 0.0, 0.0])

    def G(self, x, epar=None):
        g = self.iparv.g
        vi = x[2] - x[3]
        vlo = x[4] - x[5]
        J = np.zeros((6, 6))
        dvi = g * vlo * vlo
        dlo = 2.0 * g * vlo * vi
        for r, sgn in ((0, 1.0), (1, -1.0)):
            J[r, 2], J[r, 3] = sgn * dvi, -sgn * dvi
            J[r, 4], J[r, 5] = sgn * dlo, -sgn * dlo
        return self.toolkit.array(J)

    def C(self, x, epar=None):
        return self.toolkit.zeros((6, 6))

    def CY(self, x, w, epar=None):
        return self.toolkit.zeros((6, 6))


def _cos2_mixer(f0=1e3, psd=1e-18, rin=1e3, rout=1e3, g=1e-3, npts=400):
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    for nn in ('lo', 'vin', 'out'):
        cir.add_node(nn)
    cir['vlo'] = VSin('lo', gnd, va=1.0, freq=f0)
    cir['nsrc'] = IS('vin', gnd, i=0.0, noisePSD=psd)
    cir['rin'] = R('vin', gnd, r=rin)
    cir['M'] = _SquareMixer('out', gnd, 'vin', gnd, 'lo', gnd, g=g)
    cir['rout'] = R('out', gnd, r=rout)
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1.0 / f0, timestep=1.0 / (f0 * npts), refnode=gnd,
                  maxiterations=60)
    assert pss.converged
    irn = pss.irefnode
    k = cir.get_node_index('out')
    k = k - 1 if k > irn else k
    d = np.zeros(cir.n - 1)
    d[k] = 1.0
    return cir, pss, PAC(cir, toolkit=circuit.numeric), d, (rout * g * rin) ** 2



def _loss_osc(kind, Q=8.0, npts=480, a=0.0, idc_node=None, idc=0.0,
              period=None):
    """A van der Pol tank whose ONLY noisy element is its loss resistor.

    ⚠ THE TWO FORMS ARE THE SAME PHYSICS AND DIFFERENT MNA ROWS, which is
    the whole point of the pair.  `series` puts the loss in the inductor
    branch, so node `x` has no capacitance and its KCL row is PURELY
    ALGEBRAIC -- and that is where the resistor's noise current lands.
    `parallel` puts the equivalent loss `Rp = L/(C*Rs)` across the
    capacitor, a DIFFERENTIAL row.  For a high-`Q` tank the two agree to
    `O(1/Q^2)`, and they are matched here to 3e-6 in amplitude.

    A series loss resistor is not an exotic topology -- it is where a real
    inductor's loss physically sits.

    `a` adds the even term. A LINEAR functional of the PPV -- which is what
    a DC injection measures -- vanishes by half-wave symmetry without it,
    while the QUADRATIC one does not care; so the sensitivity tests need it
    and the `c` comparison does not.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2 * np.pi * Q)
    rs = 0.2 * mu
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0)
                       + a * mu * (u ** 2 - 2.0))
    if kind == 'series':
        cir.add_node('x')
        cir['L'] = L('v', 'x', L=1.0)
        cir['Rs'] = R('x', gnd, r=rs)
    else:
        cir['L'] = L('v', gnd, L=1.0)
        cir['Rp'] = R('v', gnd, r=1.0 / rs)
    ## ⚠ a DC current keeps the circuit AUTONOMOUS, so the period stays an
    ## unknown -- a time-varying source would make it driven and there
    ## would be no `dT` to measure at all (roadmap section 0c)
    if idc_node is not None:
        cir['Idc'] = IS(idc_node, gnd, i=idc)
    T = 2 * np.pi if period is None else float(period)
    pss = PSS(cir, method='gear', reltol=1e-12)
    x0 = np.zeros(cir.n - 1)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=250)
    assert pss.converged, '%s did not converge' % kind
    return cir, pss, PAC(cir, toolkit=circuit.numeric), rs



def _raw_pair_integrals(pss, cir):
    """Orbit integrals of the RAW pair block's equation-row PPV.

    `samples` is the pair-consistent contraction (see `ppv`), second order
    everywhere with an ADDITIVE floor of ~1e-6 |v| on its mean (the
    `h G^T z` term's own mean).  The raw first block's DC content is the
    consistent one's times `1.5 s`, `s` the pair-consistency scale
    (exactly 2/3 for an isochronous pair): a MULTIPLICATIVE error of
    `1.5 s - 1`, which is +8.1% on the bias-sensitive core (`s = 0.7207`),
    -0.16% on the series-loss tank, +0.014% on the divider -- and on a row
    whose true mean is 4e-6 |v|, as the divider's node v, 0.014% of it is
    4e-11 absolute, which is why the raw block read as "exact to 3e-11"
    there while the consistent object's additive floor exceeded the
    signal.  So on a TINY-mean row the raw block is the better DC
    estimator (its error scales with the signal), and the gates that live
    on the same-grid DC identity read it through this; on a row with a
    real mean the consistent object is (second order, no scale error).
    Measured 2026-09-05; see `test_the_raw_pair_dc_is_the_consistent_dc_times_1p5_s`.
    """
    m = cir.n - 1
    _v, info = pss.ppv()
    h = np.diff(np.asarray(info['times'], dtype=float))
    n = len(h)
    Xf = np.asarray(pss.waveform[1], dtype=float)
    ar, ac = pss._algebraic_adjoint_pattern(Xf[:, 0])
    S = np.array([pss._equation_row_ppv(np.asarray(st)[:m],
                                        Xf[:, j if j < Xf.shape[1] else -1],
                                        ar, ac)
                  for j, st in enumerate(info['samples_pair'])])
    return [float((S[:n, j] * h).sum()) for j in range(m)], info


def _vdp_ppv_method(method, npts, Q=8.0):
    """The same van der Pol under a named integrator, converged."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * Q)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    pss = PSS(cir, method=method, reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=np.array([2.0, 0.0]),
                  maxiterations=200)
    assert pss.converged, '%s/%d did not converge' % (method, npts)
    pss.monodromy = 'native'   # its callers measure the native plain path
    return cir, pss


def _coloured_vdp(kind, Q=8.0, npts=400):
    """The same physical noise two ways: `IS(noiseTau)` on the tank node,
    or a white `IS` through an RC into a linear `BSource` on that node.
    `P = g^2 Pw Rf^2`, `tau = Rf Cf`, chosen at `tau ~ 0.3 T` so that the
    source-side and output-side frequencies differ by a visible factor.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * Q)
    T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    Rf, Cf, g, Pw = 1.0, 0.3 * T, 1e-2, 1e-4
    P = g * g * Pw * Rf * Rf
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    if kind == 'coloured':
        c['n'] = IS('v', gnd, i=0.0, noisePSD=P, noiseTau=Rf * Cf)
    elif kind == 'filtered':
        c.add_node('f')
        c['nw'] = IS('f', gnd, i=0.0, noisePSD=Pw)
        c['rf'] = R('f', gnd, r=Rf)
        c['cf'] = C('f', gnd, c=Cf)
        c['gm'] = BSource('f', gnd, gnd, 'v', i_func=lambda u, _g=g: _g * u)
    else:
        c['n'] = IS('v', gnd, i=0.0, noisePSD=P)
    pss = PSS(c, method='gear', reltol=1e-12)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=300)
    assert pss.converged
    ov = [str(n) for n in c.nodes].index('v')
    return c, pss, PAC(c, toolkit=circuit.numeric), ov


def _orbit_modulated_vdp(kind, method='gear', a=0.0, kk=0.05):
    """van der Pol (Q = 8, 400 points) with a noise source whose level
    follows the orbit, and its REALISATION: the same physics as a
    STATIONARY source into an algebraic node `n`, times the modulating
    voltage through `_NuMult` (the multiplier adds no state, and `n` is 0
    on the orbit, so the modes are the same).

      'white'        `white_noise((k V_v)^2)` on the tank (`_NuModNoise`)
      'flicker'      `k V_v flicker_noise(1)`, SIGNED (`_SgnAmpFlicker`)
      'flicker_psd'  the same PSD without its sign (`_SgnPsdFlicker`)
      'lorentz'      a Lorentzian at `(k V(v, b))^2`, `b` a 3 V DC source
                     (the level keeps its sign), `_ModLorentzCtl`
      'lorentz_moving'  and its corner moving with V (`shape`)
      '<kind>_ref'   the realisation of 'white' / 'flicker' / 'lorentz'
      'white_across' a white source between `n` and the tank (1e-4), so
                     a sign error in the adjoint's algebraic entry shows
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0) + a * u * u)
    ctl = ('v', gnd)
    if kind.startswith('lorentz'):
        c.add_node('b')
        c['vb'] = VS('b', gnd, v=3.0)
        ctl = ('v', 'b')
    tau = 0.3 * 2.0 * np.pi
    if kind.endswith('_ref'):
        c.add_node('n')
        c['xi'] = {'white_ref': lambda: IS('n', gnd, i=0.0, noisePSD=1.0),
                   'flicker_ref': lambda: _Flicker('n', gnd, i=0.0, noisePSD=1.0),
                   'lorentz_ref': lambda: IS('n', gnd, i=0.0, noisePSD=1.0,
                                             noiseTau=tau)}[kind]()
        c['Rn'] = R('n', gnd, r=1.0)
        c['M'] = _NuMult('v', gnd, 'n', gnd, *ctl, k=kk)
    elif kind == 'white_across':
        c.add_node('n')
        c['xi'] = IS('n', 'v', i=0.0, noisePSD=1e-4)
        c['Rn'] = R('n', gnd, r=1.0)
        c['M'] = _NuMult('v', gnd, 'n', gnd, *ctl, k=kk)
    elif kind == 'white':
        c['src'] = _NuModNoise('v', gnd, *ctl, k=kk)
    elif kind == 'flicker':
        c['src'] = _SgnAmpFlicker('v', gnd, *ctl, k=kk)
    elif kind == 'flicker_psd':
        c['src'] = _SgnPsdFlicker('v', gnd, *ctl, k=kk)
    else:
        c['src'] = _ModLorentzCtl('v', gnd, *ctl, noisePSD=1.0, tau=tau, k=kk,
                                  shape=0.02 if kind == 'lorentz_moving' else 0.0)
    T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    pss = PSS(c, method=method, reltol=1e-12)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 400, x0=x0, maxiterations=300)
    assert pss.converged
    ov = [str(n) for n in c.nodes].index('v')
    return c, pss, PAC(c, toolkit=circuit.numeric), ov


class _ModLorentzCtl(Circuit):
    """A Lorentzian current p -> n whose LEVEL follows V(cp, cn):
    ``P (k V)^2 / (1 + (w tau(V))^2)``, ``tau(V) = tau (1 + shape V^2)``.
    `shape = 0`: a level under a fixed spectral shape (separable); `shape >
    0`: the corner moves with V too."""
    terminals = ('p', 'n', 'cp', 'cn')
    instparams = [Parameter(name='noisePSD', desc='', unit='', default=0.0),
                  Parameter(name='tau', desc='', unit='s', default=1e-7),
                  Parameter(name='k', desc='', unit='', default=1.0),
                  Parameter(name='shape', desc='', unit='', default=0.0)]

    def CY(self, x, w, epar=None):
        v = float(x[2] - x[3])
        tau = self.iparv.tau * (1.0 + self.iparv.shape * v * v)
        p = (self.iparv.noisePSD * (self.iparv.k * v) ** 2
             / (1.0 + (float(w) * tau) ** 2))
        out = np.zeros((4, 4))
        out[0, 0] = out[1, 1] = p
        out[0, 1] = out[1, 0] = -p
        return self.toolkit.array(out)


class _ModLorentzSigned(_ModLorentzCtl):
    """`_ModLorentzCtl` stating its SIGNED amplitude, ``k V sqrt(P / (1 +
    (w tau)^2))``, `shape = 0` (`Element.noise_amplitudes`)."""

    def noise_amplitudes(self, x, w=0, epar=None):
        a = self.iparv.k * float(x[2] - x[3]) * np.sqrt(
            self.iparv.noisePSD / (1.0 + (float(w) * self.iparv.tau) ** 2))
        return np.array([[a], [-a], [0.0], [0.0]])


class _TwoModLorentz(Circuit):
    """Two INDEPENDENT Lorentzian currents p -> n in one element, under two
    modulations: ``V(ap, an) sqrt(P1 L1(w))`` and ``V(bp, bn) sqrt(P2
    L2(w))``, ``L(w) = 1 / (1 + (w tau)^2)``.  Its `CY` is their sum, whose
    spectral SHAPE moves along an orbit where the two voltages' ratio does:
    a per-band source on the moving-shape path.  `signed`: whether it states
    its two columns (`Element.noise_amplitudes`)."""
    terminals = ('p', 'n', 'ap', 'an', 'bp', 'bn')
    instparams = [Parameter(name='P1', desc='', unit='', default=0.0),
                  Parameter(name='tau1', desc='', unit='s', default=1e-7),
                  Parameter(name='P2', desc='', unit='', default=0.0),
                  Parameter(name='tau2', desc='', unit='s', default=1e-8)]
    signed = False

    def _columns(self, x, w):
        p = self.iparv
        c1 = float(x[2] - x[3]) * np.sqrt(p.P1 / (1.0 + (float(w) * p.tau1) ** 2))
        c2 = float(x[4] - x[5]) * np.sqrt(p.P2 / (1.0 + (float(w) * p.tau2) ** 2))
        W = np.zeros((6, 2))
        W[0] = (c1, c2)
        W[1] = (-c1, -c2)
        return W

    def CY(self, x, w, epar=None):
        W = self._columns(x, w)
        return self.toolkit.array(W @ W.T)

    def noise_amplitudes(self, x, w=0, epar=None):
        return self._columns(x, w) if self.signed else None


def _a9_vdp(Q=8.0, psd=1e-6, cval=1.0, lval=1.0, a=0.0):
    """van der Pol with the reactances AND the half-wave symmetry as KNOBS.

    The unit-reactance, symmetric default is exactly what hid two defects on
    2026-09-07: the `q^T p` normalisation (needs `C != 1` to show) and the
    replayed adjoint being `C^T q` (needs an ASYMMETRIC orbit to show).  `a`
    adds `a*u^2` to the nonlinearity and breaks the symmetry; `cval`/`lval`
    move the reactances at fixed `w0` when `lval = 1/cval`.
    See `_hostile_oscillator` for the configuration that has BOTH.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * Q)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=cval)
    cir['L'] = L('v', gnd, L=lval)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0) + a * u * u)
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
    w0 = 1.0 / np.sqrt(cval * lval)
    T = 2.0 * np.pi / w0 / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    pss = PSS(cir, method='gear', reltol=1e-12)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 400, x0=np.array([2.0, 0.0]),
                  maxiterations=300)
    assert pss.converged
    return cir, pss


def _injection_lock_edge(cir_fn, ratio=0.2, steps=None, npts=400,
                         maxiterations=40, timing=None):
    """Continuation along the LOCKED branch of a driven oscillator; returns
    the last multiple of the named averaging prediction at which the branch
    is still stable, plus the trace.  A6's instrument.

    ⚠ NOT convergence: a forced circuit always has a periodic solution at the
    drive period, so the shooting converges on both sides of the lock edge --
    outside it merely lands on the UNSTABLE suppressed branch (amplitude
    ~0.4 V, |lam| ~1.015 on the hostile fixture).  Lock is a Floquet
    stability boundary: the forced orbit's dominant multiplier reaching 1.
    ⚠ NOT a fixed seed either: it hops branches (locked / non-convergent /
    suppressed).  Seed each step from the previous LOCKED orbit.
    """
    import warnings as _w
    from pycircuit.circuit.elements import ISin
    circuit.default_toolkit = circuit.numeric
    ## 0.2x steps: 0.1x cost 11 min per gate; the edge is reported as the
    ## last LOCKED multiple, so the coarser grid only rounds it down by at
    ## most one step, inside every bound below.
    ## ⚠ maxiterations=40, NOT 200 (measured 2026-09-08 after a py-spy dump
    ## found the hostile gate holding the full suite for 30 CPU-minutes):
    ## every LOCKED step is seeded from the previous locked orbit and
    ## converges in 1.6-3.3 s; every step PAST the edge is a solve that
    ## fails only when the cap is reached, at ~1 s per iteration (a 400-step
    ## transient plus a 400-step sensitivity traversal each), and the loop
    ## runs one or two of them by design.  The cap IS the cost of the failed
    ## steps.  At 40 the hostile gate takes 113 s and returns the same edge
    ## (2.0x) and the same |lam| trace to five digits; `timing` collects
    ## (r, seconds, solved) per step for the next time this is asked.
    steps = np.arange(0.0, 2.01, 0.2) if steps is None else steps

    def amp(cir, p):
        iv = cir.get_node_index('v')
        W = np.asarray(p.waveform[1], dtype=float)
        return (W[iv].max() - W[iv].min()) / 2

    def dom(p):
        fp = p.factored_period()
        M = np.column_stack([np.asarray(fp.matvec(e), float).ravel()
                             for e in np.eye(fp.width)])
        return float(np.max(np.abs(np.linalg.eigvals(M))))

    with _w.catch_warnings():
        _w.simplefilter('ignore')
        cir0, cval, mu = cir_fn(0.0, None)
        T0 = 2.0 * np.pi * np.sqrt(cval * cir0['L'].ipar.L)
        p = PSS(cir0, method='gear', reltol=1e-11)
        p.solve(period=T0, timestep=T0 / npts, x0=np.array([2.0, 0.0]),
                maxiterations=400)
        f0 = 1.0 / float(p.period)
        A = amp(cir0, p)
    iinj = ratio * mu * A
    ## THE PREDICTION, NAMED: averaging the vdP equation gives
    ## dphi/dt = dw - (I/(2 C A)) sin(phi), so the lock half-range is
    ## I/(2 C A) rad/s.  No Q in it -- Adler's Q and I_osc both collapse.
    pred_hz = iinj / (2.0 * cval * A) / (2.0 * np.pi)

    seed, last, trace = None, None, []
    for r in steps:
        finj = f0 + r * pred_hz
        cir, _, _ = cir_fn(iinj, finj)
        pp = PSS(cir, method='gear', reltol=1e-9)
        x0 = np.array([2.0, 0.0]) if seed is None else seed
        import time as _time
        _t0 = _time.perf_counter()
        try:
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                pp.solve(period=1.0 / finj, timestep=1.0 / finj / npts, x0=x0,
                         maxiterations=maxiterations, x0_unknown=False)
            if timing is not None:
                timing.append((float(r), _time.perf_counter() - _t0, True))
            lam, am = dom(pp), amp(cir, pp)
            locked = lam < 1.0 and am > 0.6 * A
            trace.append((float(r), am, lam, locked))
            if locked:
                seed = np.delete(np.asarray(pp.waveform[1], float)[:, 0],
                                 pp.irefnode)
                last = float(r)
            elif last is not None:
                break
        except Exception:
            if timing is not None:
                timing.append((float(r), _time.perf_counter() - _t0, False))
            trace.append((float(r), np.nan, np.nan, False))
            if last is not None and r > last + 0.35:
                break
    return last, pred_hz, f0, trace


def _vdp_injected(cval=1.0, lval=1.0, a=0.0, Q=8.0):
    from pycircuit.circuit.elements import ISin
    mu = 1.0 / (2.0 * np.pi * Q)

    def build(iinj, finj):
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=cval)
        cir['L'] = L('v', gnd, L=lval)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0) + a * u * u)
        if iinj > 0:
            cir['inj'] = ISin('v', gnd, ia=iinj, freq=finj)
        return cir, cval, mu
    return build


def _a10_vdp(Q=1e4):
    mu = 1.0 / (2.0 * np.pi * Q)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    return cir, mu


def _e5_circuit(g22, cscale=1.0, r1=1e3, r2=1e3):
    """A VCCS cancelling node 2's self conductance, so `G_22` passes through
    ZERO with `rank C` FIXED -- an index change with no rank change and no
    topology change.  `theta_0 = 1` exactly when `G_22 = 0`.  (Construction
    relayed from the docs session, 2026-09-11, as the `e.5` case.)"""
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('1')
    c.add_node('2')
    c['i'] = IS(gnd, '1', i=1e-3)
    c['c'] = C('1', gnd, c=1e-9 * cscale)
    c['r1'] = R('1', '2', r=r1)
    c['r2'] = R('2', gnd, r=r2)
    c['gm'] = VCCS('2', gnd, '2', gnd, gm=g22 - (1.0 / r1 + 1.0 / r2))
    return c


def _diode_fixture():
    circuit.default_toolkit = circuit.numeric
    per = 1e-3
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=0.8, freq=1.0 / per, vo=0.6)
    c['rs'] = R('a', 'b', r=50.0)
    c['d'] = Diode('b', gnd)
    c['cl'] = C('b', gnd, c=1e-9)
    c['rl'] = R('b', gnd, r=1e4)
    return c


def _sampler_fixture(elements, npts=400):
    """The switched-capacitor sampler of
    `test_a_switched_capacitor_holds_kTC_with_per_step_CY` with the noisy
    elements supplied: `elements(cir)` adds them between 'in'/'out'/'ck'."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    fclk = 100e3
    T = 1.0 / fclk
    cir = SubCircuit()
    for nd in ('in', 'out', 'ck'):
        cir.add_node(nd)
    cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=fclk, phase=0.0)
    cir['Vck'] = VSin('ck', gnd, vo=0.0, va=1.0, freq=fclk, phase=90.0)
    cir['C0'] = C('out', gnd, c=100e-12)
    elements(cir)
    pss = PSS(cir, method='gear', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=np.zeros(cir.n - 1),
                  maxiterations=100)
    assert pss.converged
    io = [str(nd) for nd in cir.nodes if str(nd) != 'gnd!'].index('out')
    return cir, pss, io, PAC(cir, toolkit=circuit.numeric), T


_KB, _TEMP = 1.38e-23, 300.0


def _sw(kb=_KB, **kw):
    kw.setdefault('gon', 1e-3)
    kw.setdefault('goff', 1e-9)
    return _SwitchHdl('in', 'out', 'ck', gnd, vth=0.0, vs=50e-3, temp=_TEMP,
                      kb=kb, **kw)


def _sampler_fixture_method(elements, method, npts=400):
    """`_sampler_fixture` under another integrator."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    fclk = 100e3
    T = 1.0 / fclk
    cir = SubCircuit()
    for nd in ('in', 'out', 'ck'):
        cir.add_node(nd)
    cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=fclk, phase=0.0)
    cir['Vck'] = VSin('ck', gnd, vo=0.0, va=1.0, freq=fclk, phase=90.0)
    cir['C0'] = C('out', gnd, c=100e-12)
    elements(cir)
    pss = PSS(cir, method=method, reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=np.zeros(cir.n - 1),
                  maxiterations=100)
    assert pss.converged
    io = [str(nd) for nd in cir.nodes if str(nd) != 'gnd!'].index('out')
    return cir, pss, io, PAC(cir, toolkit=circuit.numeric), T


def _pulse_clocked_sampler(T, cval=100e-12, kb=1.38e-23, temp=300.0):
    """The switched-capacitor sampler with a PULSE clock: the switch turns
    on at k T (rising ramp `tr` = T/200 at phase 0) and off at T/2."""
    cir = SubCircuit()
    cir.add_node('in')
    cir.add_node('out')
    cir.add_node('ck')
    cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=1.0 / T, phase=0.0)
    cir['Vck'] = VPulse('ck', gnd, v1=-1.0, v2=1.0, td=0.0, tr=T / 200,
                        tf=T / 200, pw=T / 2 - T / 200, per=T)
    cir['S0'] = _SwitchHdl('in', 'out', 'ck', gnd, gon=1e-3, goff=1e-9,
                           vth=0.0, vs=50e-3, temp=temp, kb=kb)
    cir['C0'] = C('out', gnd, c=cval)
    return cir


def _pwm_loop(T):
    """A voltage-mode PWM loop with a STATE-dependent switching instant: a
    sawtooth against the output filtered by an RC of 10 T, a `VSwitch`
    (20 mV `tanh` window) feeding an LC filter, a 20 ohm path for the
    inductor current in the off phase."""
    from pycircuit.circuit.elements import VSwitch
    cir = SubCircuit()
    for n_ in ('vin', 'ramp', 'sw', 'out', 'fb'):
        cir.add_node(n_)
    cir['Vin'] = VS('vin', gnd, v=5.0)
    cir['Vramp'] = VPulse('ramp', gnd, v1=0.0, v2=4.0, td=0.0, tr=0.98 * T,
                          pw=0.01 * T, tf=0.01 * T, per=T)
    cir['S'] = VSwitch('vin', 'sw', 'ramp', 'fb', Ron=0.1, Roff=1e6,
                       Von=0.01, Voff=-0.01)
    cir['Rf'] = R('out', 'fb', r=1e4)
    cir['Cf'] = C('fb', gnd, c=1e-8)
    cir['L'] = L('sw', 'out', L=2e-6)
    cir['C'] = C('out', gnd, c=1e-6)
    cir['Rl'] = R('out', gnd, r=5.0)
    cir['Rd'] = R('sw', gnd, r=20.0)
    return cir


def _comparator_relaxation_oscillator():
    """An autonomous circuit with a STATE event: RC from a rail, a sharp
    comparator (`VSwitch`, 0.2 mV window) across C driven by the capacitor
    voltage through two RC lags against a reference.  One lag settles into
    a sliding equilibrium at the threshold; two oscillate at 1.39 us."""
    from pycircuit.circuit.elements import VSwitch
    cir = SubCircuit()
    for n_ in ('vdd', 'c', 'fb0', 'fb1', 'ref'):
        cir.add_node(n_)
    cir['Vdd'] = VS('vdd', gnd, v=5.0, vac=0.0)
    cir['Vref'] = VS('ref', gnd, v=2.5, vac=0.0)
    cir['R1'] = R('vdd', 'c', r=1e3)
    cir['C1'] = C('c', gnd, c=1e-9)
    cir['R2'] = R('c', 'fb0', r=1e3)
    cir['C2'] = C('fb0', gnd, c=3e-10)
    cir['R3'] = R('fb0', 'fb1', r=1e3)
    cir['C3'] = C('fb1', gnd, c=3e-10)
    cir['S'] = VSwitch('c', gnd, 'fb1', 'ref', Ron=10.0, Roff=1e7, Von=1e-4, Voff=-1e-4)
    return cir


def _jitter_sampler(T=1e-6, V1=5.0, V2=10.0):
    """A comparator-jitter sampler -- the analytic gate for the bordered
    Lyapunov closure.  A sawtooth (slope `s1 = V1 / 0.9T`) crosses a NOISY
    threshold node `n` (`Vref` through `R_n` into `C_n`: `kT/C_n`); the
    crossing opens a VSwitch that disconnects a hold capacitor from a
    SECOND ramp (slope `s2 = V2 / 0.9T`).  The held value is `ramp2(t_c)`
    and the crossing jitters by `-dn / s1`, so ``Var(hold) = (s2/s1)^2
    kT/C_n`` -- reachable only through the crossing's motion: the switch
    has no noise model and the sources none."""
    from pycircuit.circuit import VSwitch
    cir = SubCircuit()
    for nm in ('saw', 'r2', 'ref', 'n', 'hold'):
        cir.add_node(nm)
    cir['Vsaw'] = VPulse('saw', gnd, v1=0.0, v2=V1, td=0.0, tr=0.9 * T,
                         tf=0.05 * T, pw=0.0, per=T, vac=0.0)
    cir['Vr2'] = VPulse('r2', gnd, v1=0.0, v2=V2, td=0.0, tr=0.9 * T,
                        tf=0.05 * T, pw=0.0, per=T, vac=0.0)
    cir['Vref'] = VS('ref', gnd, v=2.5, vac=0.0)
    cir['Rn'] = R('ref', 'n', r=1e4)
    cir['Cn'] = C('n', gnd, c=1e-12)
    ## ON while saw < n; the 0.2 mV window is 36 ps on the slope
    cir['S'] = VSwitch('hold', 'r2', 'n', 'saw', Ron=10.0, Roff=1e9,
                       Von=1e-4, Voff=-1e-4)
    cir['Ch'] = C('hold', gnd, c=1e-11)
    return cir





@_functools.lru_cache(maxsize=1)
def _exact_relaxation_oscillator_model():
    """The EXACT piecewise-linear model of `_comparator_relaxation_oscillator`
    with an IDEAL comparator -- the reference the staged-events tests are
    gated against.  The flow is linear in each switch state (c shorted
    through `Ron`, or `Roff`), the states (c, fb0, fb1); the orbit is found
    by shooting on the two crossings of ``h . x = fb1 = 2.5``; the period
    map is the two matrix exponentials joined by saltation matrices ``I +
    (f+ - f-) h^T / (h . f-)``.  Returned as a namespace: `A_on, b_on,
    A_off, b_off` (the flows `x' = A x + b`), `Bv` (a unit tone on the
    rail), `flow(A, b, x, t) -> (x(t), e^{At})`, `xa`/`xb` (the ON->OFF and
    OFF->ON switching states), `t_off, t_on, T`, `E_off, E_on, S0, S1` and
    `M` (the period map from just before `xa`'s switching), `orbit_at(t)`
    and `left_at(v)` -- the closure sampling a LEFT vector `v` of `M`
    along the orbit, ``(map from t to T)^T v``.  Cached: the three users
    and every frequency of the forced one share one shooting solve.
    (2026-09-23, refactor E9 item 4 -- it was copied three times.)"""
    from types import SimpleNamespace
    from scipy.linalg import expm
    from scipy.optimize import fsolve
    R1 = R2 = R3 = 1e3
    C1, C2, C3 = 1e-9, 3e-10, 3e-10
    VDD, VREF, RON, ROFF = 5.0, 2.5, 10.0, 1e7

    def sysm(g):
        A = np.array([[-(1 / R1 + 1 / R2 + g) / C1, 1 / (R2 * C1), 0.0],
                      [1 / (R2 * C2), -(1 / R2 + 1 / R3) / C2, 1 / (R3 * C2)],
                      [0.0, 1 / (R3 * C3), -1 / (R3 * C3)]])
        return A, np.array([VDD / (R1 * C1), 0.0, 0.0])

    A_on, b_on = sysm(1 / RON)
    A_off, b_off = sysm(1 / ROFF)
    Bv = np.array([1 / (R1 * C1), 0.0, 0.0])

    def flow(A, b, x, t):
        E = expm(A * t)
        return E @ x + np.linalg.solve(A, (E - np.eye(3)) @ b), E

    h = np.array([0.0, 0.0, 1.0])

    def resid(z):
        xa, t_off, t_on = z[:3], z[3], z[4]
        xb, _ = flow(A_off, b_off, xa, t_off)
        xc, _ = flow(A_on, b_on, xb, t_on)
        return np.concatenate((xc - xa, [xa[2] - VREF, xb[2] - VREF]))

    z = fsolve(resid, np.array([1.0, 2.4, 2.5, 1.1e-6, 0.3e-6]), xtol=1e-13)
    assert np.linalg.norm(resid(z)) < 1e-9
    xa, t_off, t_on = z[:3], float(z[3]), float(z[4])
    T = t_off + t_on
    xb, E_off = flow(A_off, b_off, xa, t_off)
    _xc, E_on = flow(A_on, b_on, xb, t_on)

    def salt(A_pre, b_pre, A_post, b_post, x):
        f_pre, f_post = A_pre @ x + b_pre, A_post @ x + b_post
        return np.eye(3) + np.outer(f_post - f_pre, h) / float(h @ f_pre)

    S0 = salt(A_on, b_on, A_off, b_off, xa)
    S1 = salt(A_off, b_off, A_on, b_on, xb)
    M = E_on @ S1 @ E_off @ S0

    def orbit_at(t):
        if t < t_off:
            return flow(A_off, b_off, xa, t)[0]
        return flow(A_on, b_on, xb, t - t_off)[0]

    def left_at(v):
        def at(t):
            if t < t_off:
                xt, _ = flow(A_off, b_off, xa, t)
                _, Et = flow(A_off, b_off, xt, t_off - t)
                return xt, (E_on @ S1 @ Et).T @ v
            xt, _ = flow(A_on, b_on, xb, t - t_off)
            _, Et = flow(A_on, b_on, xt, T - t)
            return xt, Et.T @ v
        return at

    return SimpleNamespace(A_on=A_on, b_on=b_on, A_off=A_off, b_off=b_off, Bv=Bv,
                           flow=flow, h=h, xa=xa, xb=xb, t_off=t_off, t_on=t_on, T=T,
                           E_off=E_off, E_on=E_on, S0=S0, S1=S1, M=M,
                           orbit_at=orbit_at, left_at=left_at)


def _relaxation_oscillator_seed(cir, phase=0.8):
    """A seed and a period hint for a PSS of `_comparator_relaxation_oscillator`
    (or a copy of it with an AC source), taken from the EXACT model's orbit at
    `phase` of its period -- where the OFF branch has settled, away from both
    crossings, so the staged solve's events land at ~0.14 T and ~0.20 T.
    Returns ``(x0_reduced, T)``; the source nodes are set, the branch
    currents left for the first Newton step.

    ⚠ This replaced six identical `lte_grid` seed solves (2026-09-23), 72.6 s
    each, which nothing else in those tests used: from either seed the staged
    solve converges to the same period (5e-9 apart at 200 points), and the
    unstaged solves stay unconverged (100 vs 200 points 8.6e-3 apart), which
    is what the autonomous-stage test pins."""
    mdl = _exact_relaxation_oscillator_model()
    names = [str(n_) for n_ in cir.nodes]
    x = np.zeros(cir.n)
    xs = mdl.orbit_at(phase * mdl.T)
    for nm, val in (('c', xs[0]), ('fb0', xs[1]), ('fb1', xs[2]),
                    ('vdd', 5.0), ('ref', 2.5)):
        x[names.index(nm)] = val
    return np.delete(x, cir.get_node_index(gnd)), mdl.T


def _exact_relaxation_oscillator_ppv():
    """The EXACT PPV of `_comparator_relaxation_oscillator` with an ideal
    comparator (`_exact_relaxation_oscillator_model`): the left null vector
    of the saltation-joined period map, normalised ``v . f = 1`` on the
    flow just before the ON->OFF switching, sampled along the orbit.
    Returns ``(T, sample)`` with ``sample(x)`` the PPV at the orbit point
    nearest the state `x = (c, fb0, fb1)`."""
    from scipy.linalg import null_space
    mdl = _exact_relaxation_oscillator_model()
    v = null_space(mdl.M.T - np.eye(3), rcond=1e-9)
    assert v.shape[1] == 1
    v = v[:, 0]
    v = v / float(v @ (mdl.A_on @ mdl.xa + mdl.b_on))
    at = mdl.left_at(v)
    tgrid = np.linspace(0.0, mdl.T, 4001)[:-1]
    orbit = np.array([at(t)[0] for t in tgrid])

    def sample(x):
        k = int(np.argmin(np.linalg.norm(orbit - np.asarray(x, dtype=float), axis=1)))
        return at(tgrid[k])[1]

    return mdl.T, sample


def _pow_flicker(ef):
    """A 1/f^ef current, `flicker_noise(k, ef)` (one-sided `k / f^ef`)."""
    from pycircuit.circuit.hdl import flicker_noise as _fn

    class _PowFlicker(Behavioural):
        params_as = 'p'
        instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

        @staticmethod
        def analog(p, outp, outn):
            return Contribution(Branch(outp, outn).I, _fn(p.k, ef))
    return _PowFlicker


def _mixed_exponent_rc(kind):
    """Two RC branches `o1`, `o2` off one driven node, and a 1/f^0.8 plus a
    1/f^2 source on them: 'mixed' one element on DISJOINT branches (a
    component whose exponent differs between entries), 'split' two
    elements, 'corr' one element whose two sources are CORRELATED (cross
    entry `sqrt(p1 p2)`, of a third slope: no split makes them
    independent), 'xcorr' the same correlated noise as the two elements
    of 'split' and a `CY` override adding the cross entry."""
    import warnings
    from pycircuit.circuit.hdl import flicker_noise as _fn
    circuit.default_toolkit = circuit.numeric
    T, Rv, Cv, k = 1e-6, 1e3, 1e-9, 1e-20

    class _Mixed(Behavioural):
        params_as = 'p'
        instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

        @staticmethod
        def analog(p, a, an, b, bn):
            return (Contribution(Branch(a, an).I, _fn(p.k, 0.8)),
                    Contribution(Branch(b, bn).I, _fn(p.k, 2.0)))

    class _Corr(Circuit):
        terminals = ('a', 'b')
        instparams = [Parameter(name='k', desc='', unit='', default=1.0)]

        def CY(self, x, w, epar=None):
            f = abs(float(w)) / (2.0 * np.pi)
            p1, p2 = self.iparv.k / f ** 0.8, self.iparv.k / f ** 2.0
            xx = np.sqrt(p1 * p2)
            return self.toolkit.array(np.array([[p1, xx], [xx, p2]]))

    class _XCorr(SubCircuit):
        def CY(self, x, w, epar=circuit.defaultepar):
            out = np.array(SubCircuit.CY(self, x, w, epar))
            f = abs(float(w)) / (2.0 * np.pi)
            xx = np.sqrt((k / f ** 0.8) * (k / f ** 2.0))
            i, j = self.get_node_index('o1'), self.get_node_index('o2')
            out[i, j] += xx
            out[j, i] += xx
            return self.toolkit.array(out)
    c = _XCorr() if kind == 'xcorr' else SubCircuit()
    for nd in ('in', 'o1', 'o2'):
        c.add_node(nd)
    c['V'] = VSin('in', gnd, va=0.1, vo=0.0, freq=1.0 / T)
    for b in ('1', '2'):
        c['R' + b] = R('in', 'o' + b, r=Rv, noisy=False)
        c['C' + b] = C('o' + b, gnd, c=Cv)
    if kind == 'mixed':
        c['n'] = _Mixed('o1', gnd, 'o2', gnd, k=k)
    elif kind == 'corr':
        c['n'] = _Corr('o1', 'o2', k=k)
    else:
        c['n1'] = _pow_flicker(0.8)('o1', gnd, k=k)
        c['n2'] = _pow_flicker(2.0)('o2', gnd, k=k)
    pss = PSS(c, method='radau', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 100, maxiterations=40)
    names = [str(x) for x in c.nodes if str(x) != 'gnd!']
    return pss, PAC(c, toolkit=circuit.numeric), names, (T, Rv, Cv, k)



def _vdp_asym():
    """Van der Pol in LC form with an asymmetric term (period `T_REF` =
    6.330195895892 s, radau at 800 points)."""
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=4.0)
    cir['L'] = L('v', gnd, L=0.25)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: (u - u ** 3 / 3.0) + 0.3 * u * u)
    return cir
