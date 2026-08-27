# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""`limit_identity`: a ``$limit`` probe with no law (roadmap sec. 15, the
proposal at the end of Stage 2).

Vector PCNR admits a device only when every resistive current reaches the
solution through its declared probes.  A MOSFET's current reads three
branch voltages and SPICE's limiter set covers two of them, so
`EkvNmosHdl` was refused with ``var(vsb) reads b, g, which no $limit probe
limits``.  An identity probe gives the third branch a PCNR unknown whose
law is "leave it alone".

Three claims, each with its own control:

1. **On the ordinary path it is exactly nothing.**  A model whose only
   probe is an identity probe returns its input bit for bit from
   ``limit()``; beside a real probe it neither moves nor is moved by the
   write-back beyond what the real probe does alone (section 1).
2. **Under vector PCNR it qualifies the device**, and the gmin pair view
   -- the pnj-only list every ladder reads -- does not see it (section 2).
3. **Measured on the diff pair**, EKV, plain Newton against PCNR, both
   instance orders, from a zero start and from +20 V (section 3).
   Order-independence and convergence are asserted separately; every
   converged tail equals ``DC()``'s.

The controls are mutation-checked: the identity kind's `apply_limit`
arm replaced by ``vold`` breaks 1; the EKV declaration reverted to a raw
``bsb.V`` breaks 2 and 3; the pair view filter widened to all kinds breaks
the gmin test.
"""
import warnings

import numpy as np
import pytest
import sympy
from numpy.testing import assert_allclose

import pycircuit.circuit.circuit
from pycircuit.circuit import numeric, gnd
from pycircuit.circuit.circuit import SubCircuit, defaultepar
from pycircuit.circuit.elements import VS, IS as ISRC, R
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import pcnr as P
from pycircuit.circuit._limiting import apply_limit
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   explain, limit_identity, limit_fet,
                                   limit_pnj, limit_together, softplus,
                                   var, vt)
from pycircuit.utilities.param import Parameter


@pytest.fixture(autouse=True)
def _numeric_toolkit():
    old = pycircuit.circuit.circuit.default_toolkit
    pycircuit.circuit.circuit.default_toolkit = numeric
    yield
    pycircuit.circuit.circuit.default_toolkit = old


## ======================================================================
## 1.  The ordinary path: nothing, exactly.
## ======================================================================

class _IdOnly(Behavioural):
    """A three-terminal device whose current reads two branches and
    declares an identity probe on one of them."""
    instparams = [Parameter(name='g', desc='', unit='S', default=1e-3)]

    @staticmethod
    def analog(a, b, c):
        bab, bbc = Branch(a, b), Branch(b, c)
        v1 = var(limit_identity(bab.V), 'v1')
        return Contribution(bab.I, g * v1 * (1 + bbc.V ** 2))           # noqa


class _IdBesideFet(Behavioural):
    """An exponential 'FET': `fetlim` on (g,s), identity on (s,b)."""
    instparams = [Parameter(name='vto', desc='', unit='V', default=0.5)]

    @staticmethod
    def analog(d, g, s, b):
        bgs, bds, bsb = Branch(g, s), Branch(d, s), Branch(s, b)
        vgs = var(limit_fet(bgs.V, vto), 'vgs')                        # noqa
        vsb = var(limit_identity(bsb.V), 'vsb')
        return Contribution(bds.I, 1e-6 * softplus((vgs - vto        # noqa
                                                    + 0.3 * vsb) / vt()) ** 2
                            * sympy.tanh(bds.V))


class _FetOnly(Behavioural):
    """`_IdBesideFet` without the identity probe: the CONTROL for
    "beside a real probe the identity changes nothing"."""
    instparams = [Parameter(name='vto', desc='', unit='V', default=0.5)]

    @staticmethod
    def analog(d, g, s, b):
        bgs, bds, bsb = Branch(g, s), Branch(d, s), Branch(s, b)
        vgs = var(limit_fet(bgs.V, vto), 'vgs')                        # noqa
        vsb = var(bsb.V, 'vsb')
        return Contribution(bds.I, 1e-6 * softplus((vgs - vto        # noqa
                                                    + 0.3 * vsb) / vt()) ** 2
                            * sympy.tanh(bds.V))


def test_apply_limit_identity_returns_the_proposal_itself():
    """The kind's whole law.  ``vnew`` back, the object itself -- not
    ``vold + (vnew - vold)``, which is not ``vnew`` in floating point and
    would make the write-back's ``vlim == vnew`` test fail at random."""
    rng = np.random.default_rng(3)
    for _ in range(200):
        vn, vo = float(rng.normal(scale=1e3)), float(rng.normal())
        out = apply_limit('id', vn, vo, (), numeric)
        assert out == vn and out is vn


def test_an_identity_only_device_limits_nothing_bit_for_bit():
    """`limit(x, x0)` is the identity on the whole sub-vector, for wild
    steps from any last-accepted point.  Bitwise, because the state-free
    convention lets "did limiting fire?" be a convergence signal and a
    single ulp of drift would fire it forever."""
    el = _IdOnly(*[eh.Node(nm) for nm in 'abc'])
    el.update_iparv()
    txt = explain(el, source=False, symbolic=False)
    assert '1 $limit probe (identity (no law) on (a,b))' in txt, txt
    rng = np.random.default_rng(11)
    for _ in range(100):
        x0 = rng.normal(size=3)
        x = x0 + rng.normal(scale=10 ** rng.uniform(-6, 8), size=3)
        out = el.limit(x, x0, defaultepar)
        assert np.array_equal(out, x), (x, out)


def test_beside_a_real_probe_the_identity_probe_changes_no_write_back():
    """The same device with and without the identity probe: `limit()`
    returns the same vector to the bit wherever the real probe bites,
    and the identity probe's own branch is never the one written.

    The control is the `_FetOnly` twin, compiled from the same text with
    the declaration removed."""
    a = _IdBesideFet(*[eh.Node(nm) for nm in 'dgsb'])
    b = _FetOnly(*[eh.Node(nm) for nm in 'dgsb'])
    for el in (a, b):
        el.update_iparv()
    assert 'identity (no law) on (s,b)' in explain(a, source=False)
    assert 'identity (no law)' not in explain(b, source=False)
    rng = np.random.default_rng(5)
    bit = 0
    for _ in range(300):
        x0 = rng.normal(size=4)
        x = x0 + rng.normal(scale=10 ** rng.uniform(-2, 3), size=4)
        oa = a.limit(x, x0, defaultepar)
        ob = b.limit(x, x0, defaultepar)
        assert np.array_equal(oa, ob), (x, oa, ob)
        bit += int(not np.array_equal(oa, x))
    assert bit > 100, bit          # the real probe did fire, often


def test_an_identity_probe_cannot_be_grouped():
    """A group holds a probe that did not bite as a constraint; an
    identity probe has nothing to hold."""
    b1, b2 = Branch(eh.Node('a'), eh.Node('b')), Branch(eh.Node('b'),
                                                        eh.Node('c'))
    with pytest.raises(ValueError, match='no law'):
        limit_together(limit_fet(b1.V, 0.5), limit_identity(b2.V))
    with pytest.raises(ValueError, match='branch potential'):
        limit_identity(b1.I)


## ======================================================================
## 2.  Under vector PCNR: the device qualifies; the ladders never see it.
## ======================================================================

def _ekv_explain():
    return explain(eh.EkvNmosHdl, source=False, symbolic=False)


def test_ekv_qualifies_for_vector_pcnr_through_its_bulk_probe():
    """The motivating case.  Without a probe on the third branch:
    ``PCNR: does not qualify -- var(vsb) reads b, g, which no $limit
    probe limits``.  With one, the three probes are the device's block.

    ⚠ **The bulk probe was `limit_identity` until 2026-08-26 and is now
    `limit_delta` (roadmap sec. 34).**  What qualifies the device is
    having a probe on that branch at all -- which is what this test is
    for, and it is unchanged.  What the identity probe could NOT do is
    limit the seed, and that was the last circuit on the PCNR fallback.
    The identity mechanism itself is still supported and still tested,
    in section 1, on a synthetic device; it simply has no library user
    any more.
    """
    txt = _ekv_explain()
    assert ('3 $limit probes (deltalim on (s,b), fetlim on (g,s) '
            '[params at last iterate], limvds on (d,s))') in txt, txt
    assert ('PCNR: vector route, 3 probes over 4 rows (deltalim on (s,b), '
            'fetlim on (g,s), limvds on (d,s))') in txt, txt
    assert 'does not qualify' not in txt
    el = eh.EkvNmosHdl(*[eh.Node(nm) for nm in 'dgsb'])
    assert tuple(k for _p, k in el.pcnr_probes) == ('delta', 'fet', 'vds')
    ## and the PMOS, whose branches are declared reversed
    txt = explain(eh.EkvPmosHdl, source=False, symbolic=False)
    assert 'PCNR: vector route, 3 probes over 4 rows' in txt, txt
    assert 'deltalim on (b,s)' in txt, txt


def test_the_refusal_still_fires_without_a_bulk_probe():
    """The control for the test above: the SAME EKV body with the raw
    ``bsb.V`` -- no probe of any kind -- is refused by name.  A qualification test that could not
    fail would be decoration."""
    from pycircuit.circuit.hdl import TEMP
    import pycircuit.circuit.elements_hdl as m
    src_ok = m._ekv_analog(TEMP, +1)
    cls = type('EkvRaw', (Behavioural,),
               dict(instparams=m._ekv_params(),
                    analog=staticmethod(_ekv_without_identity())))
    txt = explain(cls, source=False, symbolic=False)
    assert 'does not qualify' in txt and 'which no $limit probe limits' \
        in txt, txt
    assert 'var(vsb) reads b, g' in txt, txt
    del src_ok


def _ekv_without_identity():
    """`_ekv_analog(TEMP, +1)` with the bulk probe made a no-op, built
    by shadowing the declaring name in the body's globals."""
    import types
    import pycircuit.circuit.elements_hdl as m
    from pycircuit.circuit.hdl import TEMP
    g = dict(m._ekv_analog.__globals__)
    g['limit_delta'] = lambda probe, vmax: probe
    f = types.FunctionType(m._ekv_analog.__code__, g, '_ekv_analog',
                           m._ekv_analog.__defaults__,
                           m._ekv_analog.__closure__)
    return f(TEMP, +1)


def test_the_gmin_pair_view_stays_pnj_only():
    """`pcnr_junction_pairs` is the GMIN-target list the ladders read on
    both paths.  A MOSFET pair contributes nothing to it -- an identity
    probe on ``(s,b)`` there would put a conductance across the bulk
    junction of every EKV device in the circuit.  The device IS a PCNR
    participant (`pcnr_devices`) all the same."""
    c = _diffpair(0.3, False)
    assert P.pcnr_junction_pairs(c) == []
    devs = P.pcnr_devices(c)
    assert [d.instance for d in devs] == ['M1', 'M2']
    assert all(d.m == 3 and d.kinds == ('delta', 'fet', 'vds')
               for d in devs)
    ## and a real junction in the same circuit is still listed
    c['D'] = eh.DiodeHdl('d1', gnd)
    pairs = P.pcnr_junction_pairs(c)
    assert [p[0] for p in pairs] == ['D'], pairs


## ======================================================================
## 3.  The measurement: the diff pair, EKV, both orders, two starts.
## ======================================================================

#: 20 um wide EKV: a device that sits in moderate inversion at 100 uA per
#: side, which is where the softplus channel has an exponential to
#: compress.  The 1 um default card puts 2.4 V of overdrive on the
#: tail and every solver takes seven iterations from anywhere.
CARD = dict(w=20e-6)


def _dev(*nodes):
    return eh.EkvNmosHdl(*nodes, gnd, **CARD)


def _diffpair(vin, m2_first, vdd=5.0):
    """`test_device_limiter._diffpair` with the bulk on ground."""
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=vdd)
    c['vp'] = VS('inp', gnd, v=2.5 + vin)
    c['vn'] = VS('inn', gnd, v=2.5 - vin)
    c['rl1'] = R('vdd', 'd1', r=5e3)
    c['rl2'] = R('vdd', 'd2', r=5e3)
    c['itail'] = ISRC('tail', gnd, i=200e-6)
    devs = [('M1', ('d1', 'inp', 'tail')), ('M2', ('d2', 'inn', 'tail'))]
    for name, nodes in (devs[::-1] if m2_first else devs):
        c[name] = _dev(*nodes)
    return c


def _x_start(c, level):
    n = len(c.nodes) + len(c.branches)
    x0 = np.full(n, float(level))
    x0[c.get_node_index(gnd)] = 0.0
    return x0


def _plain_newton_from(c, x0, maxiter=100, reltol=1e-4, vabstol=1e-6,
                       iabstol=1e-12):
    """`test_limit_fet._plain_newton` with a starting point: `StandardNewton`
    with the element limiters and nothing else -- no ladder, no anchor."""
    from pycircuit.circuit.nrsolver import StandardNewton
    iref = c.get_node_index(gnd)
    epar = defaultepar

    def rfunc(xr):
        x = np.insert(xr, iref, 0.0)
        F = c.i(x, epar) + c.u(0, analysis='dc', epar=epar)
        J = c.G(x, epar)
        return (np.delete(F, iref),
                np.delete(np.delete(J, iref, 0), iref, 1))

    def limiter(xr, x0r):
        x = c.limit(np.insert(xr, iref, 0.0), np.insert(x0r, iref, 0.0),
                    epar)
        return np.delete(x, iref)

    nn, nb = len(c.nodes), len(c.branches)
    abstol = np.delete(np.concatenate((iabstol * np.ones(nn),
                                       vabstol * np.ones(nb))), iref)
    xtol = np.delete(np.concatenate((vabstol * np.ones(nn),
                                     iabstol * np.ones(nb))), iref)
    x, its = StandardNewton().solve_system(
        np.delete(x0, iref), rfunc, numeric, reltol, abstol, xtol, maxiter,
        limiter=limiter)
    return np.insert(x, iref, 0.0), its


def _row(vin, start, maxiter=100):
    """`[M1 first, M2 first]` of ``(plain, pcnr)`` -- each ``(its, tail)``,
    ``(None, None)`` on failure."""
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix
    row = []
    for m2 in (False, True):
        c = _diffpair(vin, m2)
        it = c.get_node_index('tail')
        try:
            x, its = _plain_newton_from(c, _x_start(c, start), maxiter)
            plain = (its, float(x[it]))
        except (NoConvergenceError, SingularMatrix, FloatingPointError,
                np.linalg.LinAlgError, OverflowError, ValueError):
            plain = (None, None)
        c = _diffpair(vin, m2)
        try:
            x, _v, its = P.solve_dc(c, gnd, x0=_x_start(c, start),
                                    reltol=1e-4, abstol=1e-6,
                                    maxiter=maxiter)
            pc = (its, float(x[it]))
        except (RuntimeError, np.linalg.LinAlgError, FloatingPointError,
                OverflowError, ValueError):
            pc = (None, None)
        row.append((plain, pc))
    return row


def _dc_tail(vin):
    c = _diffpair(vin, False)
    r = DC(c).solve()
    return float(r.x[c.get_node_index('tail')])


VIN = (-1.0, -0.3, 0.0, 0.3, 1.0)


@pytest.mark.filterwarnings('ignore:overflow encountered')
@pytest.mark.filterwarnings('ignore:invalid value')
def test_the_ekv_diff_pair_from_zero_plain_newton_against_pcnr():
    """Roadmap sec. 15 Stage 0's first signature was ``diff pair, EKV,
    start 0: 5 vs 7, order-dependent count``.  Measured here on the
    shipped EKV with the identity probe, `[M1 first, M2 first]`,
    Jacobian evaluations (recorded at introduction)::

        vin           -1.0      -0.3       0       0.3      1.0
        plain        [7, 5]    [7, 5]   [7, 7]   [5, 7]   [5, 7]
        PCNR         [7, 7]    [7, 7]   [7, 7]   [7, 7]   [7, 7]

    Two claims, separately:

    1. **Order-independence**: PCNR's counts are equal in both orders at
       every ``vin`` and the table is symmetric in ``vin``.  Plain
       Newton's are NOT -- the Stage 0 signature reproduces -- and that
       asymmetry is asserted too, so that this test cannot pass on a
       circuit where there was nothing to fix.
    2. **Convergence**: both converge everywhere from a zero start, to
       ``DC()``'s tail.  There is no failure to rescue here; the price
       of order-independence is the lucky order's two iterations.
    """
    rows = {vin: _row(vin, 0.0) for vin in VIN}
    for vin, row in rows.items():
        ref = _dc_tail(vin)
        (pl1, pc1), (pl2, pc2) = row
        ## claim 1
        assert pc1[0] == pc2[0], (vin, row)
        ## claim 2
        for its, tail in (pl1, pl2, pc1, pc2):
            assert its is not None and its <= 20, (vin, row)
            assert_allclose(tail, ref, rtol=1e-5)
    assert rows[1.0][0][1][0] == rows[-1.0][0][1][0], rows
    assert rows[0.3][0][1][0] == rows[-0.3][0][1][0], rows
    ## plain Newton is order-DEPENDENT on this circuit, as Stage 0 found
    plain_asym = [vin for vin in VIN
                  if rows[vin][0][0][0] != rows[vin][1][0][0]]
    assert plain_asym, rows


@pytest.mark.filterwarnings('ignore:overflow encountered')
@pytest.mark.filterwarnings('ignore:invalid value')
def test_the_ekv_diff_pair_from_twenty_volts():
    """The same pair from a uniform +20 V start, Stage 1's wild-start
    case.  **PCNR now converges here, and the test that used to record
    its failure asserts the convergence instead.**

    ::

        vin           -1.0      -0.3       0       0.3      1.0
        plain        [6, 6]    [6, 6]   [6, 6]   [6, 6]   [6, 6]
        PCNR         [8, 8]    [7, 7]   [7, 7]   [7, 7]   [8, 8]
        PCNR before  [F, F]    [F, F]   [F, F]   [F, F]   [F, F]

    The failure it used to record, kept because the mechanism is the
    interesting part: PCNR failed at its FIRST `predict`,
    ``LinAlgError: Singular matrix``, in both orders and from any
    uniform start of 2 V or more.  A uniform start has both gates 17 V
    below the tail, both channels cut off with ``iff = irr = 0`` to the
    last bit, and the tail node held by nothing but the ideal current
    source -- its whole conductance is the two devices' 2e-19 S.  The
    plain path's LU happens to pivot on that 2e-19 and the limiter
    repairs the step; PCNR's Schur complement rounded the same row to a
    zero pivot.

    What fixed it was **not** a damped `predict`.  Sec. 27 found the
    wild-start failures were a SEED problem -- `v_lim_init` handed the
    first Jacobian raw branch voltages -- and the fix passes the seed
    through each device's own law.  This pair was the one circuit that
    fix could not reach, because its bulk probe was `limit_identity`,
    whose law is "leave it alone".  Sec. 34's `limit_delta` gave that
    probe a real law and the failure went with it.

    The old docstring said "a future damped `predict` that converges
    here should make this docstring stale, not the test red" -- and the
    test was written so that it did exactly that.  It was right about the
    outcome and wrong about the cause.

    Asserted now: both orders agree, **both converge**, and every tail
    equals ``DC()``'s.
    """
    rows = {vin: _row(vin, 20.0, maxiter=200) for vin in VIN}
    for vin, row in rows.items():
        ref = _dc_tail(vin)
        (pl1, pc1), (pl2, pc2) = row
        assert pc1[0] == pc2[0], (vin, row)
        assert pl1[0] is not None and pl2[0] is not None, (vin, row)
        ## PCNR converges here since sec. 34.  Asserted, not merely
        ## recorded -- otherwise a regression back to [F, F] would pass.
        assert pc1[0] is not None and pc2[0] is not None, (vin, row)
        assert pc1[0] < 20, (vin, row)
        for its, tail in (pl1, pl2, pc1, pc2):
            if its is not None:
                assert_allclose(tail, ref, rtol=1e-4)


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_dc_pcnr_takes_the_ekv_pair_through_pcnr_and_agrees_with_dc():
    """`DC(pcnr=True)` gates on `pcnr_devices()` (roadmap sec. 15,
    Stage 2's second defect), so the pair -- which has no junction in the
    pnj-only view -- must go through `pcnr.solve_dc` and not fall
    through.  Counted by wrapping `solve_dc`."""
    calls = []
    real = P.solve_dc

    def counted(*a, **kw):
        calls.append(1)
        return real(*a, **kw)

    c = _diffpair(0.3, False)
    ref = _dc_tail(0.3)
    P.solve_dc = counted
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            dc = DC(c, pcnr=True)
            r = dc.solve()
    finally:
        P.solve_dc = real
    assert calls == [1], calls
    assert dc.pcnr_fell_back is False
    assert_allclose(float(r.x[c.get_node_index('tail')]), ref, rtol=1e-6)
