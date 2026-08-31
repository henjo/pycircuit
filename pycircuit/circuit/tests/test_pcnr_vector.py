# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Vector PCNR, Stage 1: m > 1 limited unknowns per device (roadmap sec. 15).

The unit is a DEVICE'S PROBE VECTOR: a Gummel-Poon BJT owns ``v = (vbe,
vbc)`` and answers ``pcnr_i(v) -> R^3``, ``pcnr_didv(v) -> R^{3x2}``.
Scalar junctions are the ``m = 1`` case behind `pcnr.PcnrDevice`, and
section 1 below runs the ORIGINAL formulas (copied from the pre-Stage-1
`pcnr.py`) beside the new ones on every scalar participant in the tree.

**What was measured while building it, and what the design got wrong.**
Two premises of the design of record did not survive contact:

* *"assemble everything, then subtract the participant's own stamp"* is
  exact in arithmetic only.  A participant's current at the NODE
  voltages is unbounded (PCNR limits ``v_lim``, not the nodes), and on a
  BJT mirror from a zero start the iterate visits vbe = 5.2 V on the
  nodes with `expl` keeping the current finite at ~1e72: the subtraction
  left ~1e56 of noise in the base row, and WHICH noise depended on the
  assembly's summation order -- 11 iterations in one instance order, 179
  in the other, from identical states.  The layer now EXCLUDES the
  participant from the assembly (section 1 pins that the participant's
  ``i`` is never called at the node voltage).  The same defect had been
  making the charge-diode DC test take 14 iterations instead of 8: at
  every other iterate the node sat at 5 V, ``(1e-3 + 4.6e72) - 4.6e72``
  lost the resistor's conductance, and the step was nonsense.
* *"under PCNR the Stage 0 mirror converges in both orders from +20 V"*
  -- no.  Section 5's table: PCNR is order-INDEPENDENT at every start
  (the structural claim, and it holds), converges where plain Newton
  fails in BOTH orders (+5 V), and fails identically in both orders
  where plain Newton was `[FAIL, 7]` (+20 V).  The mechanism is traced
  in the docstring there; it is not the layer's.
"""
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pycircuit.circuit.circuit
from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit, defaultepar
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import R, VS, Diode
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution, explain,
                                   limexp, limit_pnj, var, vt)
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import pcnr as P
from pycircuit.circuit._limiting import _pnjlim
from pycircuit.circuit.nrsolver import NoConvergenceError, StandardNewton
from pycircuit.circuit.analysis import SingularMatrix
from pycircuit.circuit.tests.test_elements_hdl_library3 import NPN, NPN_IDEAL
from pycircuit.circuit.tests.test_pcnr import _fig1
from pycircuit.circuit.tests.test_pcnr_charge import (
    _two_junction_circuit, _circuit as _charge_diode)
from pycircuit.utilities.param import Parameter


@pytest.fixture(autouse=True)
def _numeric_toolkit():
    old = pycircuit.circuit.circuit.default_toolkit
    pycircuit.circuit.circuit.default_toolkit = numeric
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        yield
    pycircuit.circuit.circuit.default_toolkit = old


## ======================================================================
## 1.  The adapter is invisible: the ORIGINAL formulas, copied, beside
##     the new ones on every scalar participant.
## ======================================================================

def _old_augmented_system(circuit, x, v_lim, junctions, epar=defaultepar,
                          u_extra=None):
    """`pcnr.augmented_system` as it stood before Stage 1 (dense blocks),
    including its assemble-then-subtract of each participant."""
    n, k = circuit.n, len(junctions)
    g_mna = np.array(circuit.i(x, epar), dtype=float)
    J_mm = np.array(circuit.G(x, epar), dtype=float)
    if u_extra is not None:
        g_mna = g_mna + np.asarray(u_extra, dtype=float)
    nodemap = circuit.elementnodemap
    seen = set()
    for instance, element, _ra, _rb in junctions:
        if instance in seen:
            continue
        seen.add(instance)
        rows = nodemap[instance]
        sub = x[rows]
        g_mna[rows] -= np.asarray(element.i(sub, epar), dtype=float)
        J_mm[np.ix_(rows, rows)] -= np.asarray(element.G(sub, epar),
                                               dtype=float)
    J_ml, J_lm, g_lim, didv = np.zeros((n, k)), np.zeros((k, n)), \
        np.zeros(k), []
    _seen_of = {}
    for idx, (instance, element, ra, rb) in enumerate(junctions):
        jn = _seen_of.get(instance, 0)
        _seen_of[instance] = jn + 1
        v = float(v_lim[idx])
        params = {p.name: getattr(element.iparv, p.name)
                  for p in element.instparams}
        i_terms = element.pcnr_i(v, params, epar, element.toolkit, jn=jn)
        di_terms = element.pcnr_didv(v, params, epar, element.toolkit, jn=jn)
        g_mna[ra] += float(i_terms[0])
        g_mna[rb] += float(i_terms[1])
        J_ml[ra, idx] += float(di_terms[0])
        J_ml[rb, idx] += float(di_terms[1])
        didv.append((float(di_terms[0]), float(di_terms[1])))
        g_lim[idx] = v - (float(x[ra]) - float(x[rb]))
        J_lm[idx, ra] = -1.0
        J_lm[idx, rb] = +1.0
    return g_mna, g_lim, J_mm, J_ml, J_lm, didv


def _old_schur(g_mna, g_lim, J_mm, junctions, didv):
    schur = np.array(J_mm, copy=True)
    rhs_corr = np.zeros(len(g_mna))
    for idx, (_i, _e, ra, rb) in enumerate(junctions):
        dia, dib = didv[idx]
        schur[ra, ra] += dia
        schur[ra, rb] -= dia
        schur[rb, ra] += dib
        schur[rb, rb] -= dib
        rhs_corr[ra] += dia * g_lim[idx]
        rhs_corr[rb] += dib * g_lim[idx]
    return g_mna - rhs_corr, schur


def _old_refine(junctions, v_old, v_new, epar=defaultepar):
    out = np.array(v_new, dtype=float, copy=True)
    _seen_of = {}
    for idx, (inst, element, _ra, _rb) in enumerate(junctions):
        jn = _seen_of.get(inst, 0)
        _seen_of[inst] = jn + 1
        limiter = getattr(element, 'pcnr_limit', None)
        if limiter is not None:
            params = {q.name: getattr(element.iparv, q.name)
                      for q in element.instparams}
            out[idx] = float(limiter(float(v_new[idx]), float(v_old[idx]),
                                     params, epar, element.toolkit, jn=jn))
            continue
        VT = element.toolkit.kboltzmann * epar.T / element.toolkit.qelectron
        out[idx] = _pnjlim(float(v_new[idx]), float(v_old[idx]), VT,
                           getattr(element.iparv, 'IS', 0.0), element.toolkit)
    return out


def _old_solve_dc(circuit, refnode, reltol=1e-6, abstol=1e-12, maxiter=200):
    junctions = P.pcnr_junctions(circuit)
    n, iref = circuit.n, circuit.get_node_index(refnode)
    x = np.zeros(n)
    v_lim = np.array([float(x[ra] - x[rb]) for _i, _e, ra, rb in junctions])
    u = np.asarray(circuit.u(0.0, defaultepar, analysis='dc'), dtype=float)
    for it in range(maxiter):
        g_mna, g_lim, J_mm, _Jml, _Jlm, didv = _old_augmented_system(
            circuit, x, v_lim, junctions, u_extra=u)
        f_eff, schur = _old_schur(g_mna, g_lim, J_mm, junctions, didv)
        keep = [i for i in range(n) if i != iref]
        dx = np.zeros(n)
        dx[keep] = np.linalg.solve(schur[np.ix_(keep, keep)], -f_eff[keep])
        dx_lim = np.array([-(g_lim[k] + (-dx[ra] + dx[rb]))
                           for k, (_i, _e, ra, rb) in enumerate(junctions)])
        x_new = x + dx
        x_new[iref] = 0.0
        v_new = _old_refine(junctions, v_lim, v_lim + dx_lim)
        conv = (np.max(np.abs(dx)) < reltol * np.max(np.abs(x_new)) + abstol
                and np.max(np.abs(g_lim))
                < reltol * max(np.max(np.abs(v_new)), 1.0) + abstol)
        x, v_lim = x_new, v_new
        if conv:
            return x, v_lim, it + 1
    raise RuntimeError('old solve_dc did not converge')


SCALAR_CIRCUITS = [
    ('fig1 two diodes', lambda: _fig1()),
    ('fig1 IS 1e-15/1e-12', lambda: _fig1(is1=1e-15, is2=1e-12)),
    ('two-junction', _two_junction_circuit),
    ('charge diode 1 V', lambda: _charge_diode(vsrc=1.0)),
]


@pytest.mark.parametrize('name, mk', SCALAR_CIRCUITS,
                         ids=[n for n, _m in SCALAR_CIRCUITS])
def test_the_adapter_reproduces_the_original_blocks_on_scalar_participants(
        name, mk):
    """Every block the old code produced, entry for entry, at a sane
    iterate.  `g_lim`, `J_ml`, `J_lm` and the per-probe derivatives are
    asserted BIT-IDENTICAL; `g_mna` and `J_mm` too, at these points,
    although the claim there is only "to rounding": the old code
    assembled the participant and subtracted it again, the new one never
    assembles it, and ``(a + d) - d`` is not ``a`` in floating point
    once ``d`` dominates ``a``.  Section 1's docstring says what that
    cost the old code."""
    c = mk()
    x = np.linspace(0.1, 0.6, c.n)
    pairs = P.pcnr_junctions(c)
    v = np.full(len(pairs), 0.55)
    u = np.asarray(c.u(0.0, defaultepar, analysis='dc'), dtype=float)
    go = _old_augmented_system(c, x, v, pairs, u_extra=u)
    gn = P.augmented_system(c, x, v, pairs, defaultepar, u_extra=u,
                            dense_blocks=True)
    assert np.array_equal(go[1], gn[1]), 'g_lim'
    assert np.array_equal(go[3], gn[3]), 'J_ml'
    assert np.array_equal(go[4], gn[4]), 'J_lm'
    assert np.array_equal(go[0], gn[0]), 'g_mna'
    assert np.array_equal(go[2], gn[2]), 'J_mm'
    ## The per-probe derivatives, placed into the new (rows, probes, block)
    ## shape: block[local a, j] == dia, block[local b, j] == dib, and
    ## nothing else in the column.
    flat = 0
    for rows, probes, block in gn[5]:
        for j, (ra, rb) in enumerate(probes):
            dia, dib = go[5][flat]
            flat += 1
            la, lb = rows.index(ra), rows.index(rb)
            col = np.zeros(len(rows))
            col[la] += dia
            col[lb] += dib
            assert np.array_equal(block[:, j], col), (name, j)
    assert flat == len(pairs)
    ## And the reduced system the two hand to the linear solver.
    fo, Jo = _old_schur(go[0], go[1], go[2], pairs, go[5])
    fn, Jn = P.schur_reduce(gn[0], gn[1], gn[2], junctions=pairs, didv=gn[5])
    assert np.array_equal(fo, fn) and np.array_equal(Jo, Jn)
    ## `refine` per device == the old per-junction loop.
    v_new = v + np.linspace(0.0, 0.3, len(pairs))
    assert np.array_equal(_old_refine(pairs, v, v_new),
                          P.refine(pairs, v, v_new, defaultepar, x_old=x))


@pytest.mark.parametrize('name, mk, its', [
    ('fig1 two diodes', lambda: _fig1(), 8),
    ('fig1 IS 1e-15/1e-12', lambda: _fig1(is1=1e-15, is2=1e-12), 14),
    ('two-junction', _two_junction_circuit, 15),
], ids=['fig1', 'fig1-1e-12', 'two-junction'])
def test_scalar_participants_take_the_same_iterations_as_before(name, mk, its):
    """Same count, same answer, old solver against new, on the circuits
    whose participant never dwarfs its row.  The counts are pinned as
    measured against the pre-Stage-1 module (8, 14, 15)."""
    xo, vo, io = _old_solve_dc(mk(), gnd)
    xn, vn, i_n = P.solve_dc(mk(), gnd)
    assert io == i_n == its, (name, io, i_n)
    assert_allclose(xn, xo, rtol=0, atol=1e-12)
    assert_allclose(vn, vo, rtol=0, atol=1e-12)


def test_the_charge_diode_at_5_v_took_14_iterations_because_of_the_subtraction():
    """THE NUMBER THAT MOVED, and why it was wrong before.

    `test_pcnr.py::test_a_charge_storing_pcnr_device_is_accepted` records
    "the same 14 DC iterations".  Traced: at every other iterate the
    node sat at the 5 V source, the diode's `exp` there is 1e71 A / 4.6e72
    S, and the old ``J_mm[rows, rows] -= G(sub)`` computed
    ``(1e-3 + 4.6e72) - 4.6e72 = 0`` -- the 1 kOhm source resistor's
    conductance was gone from the row, the step was nonsense, and the
    node bounced between 5 V and ``v_lim`` for seven round trips.  With
    the participant excluded from the assembly the resistor survives and
    the same solve takes 9 monotone iterations to the same point.

    (9, not the 8 first recorded: Stage 2 gave `solve_dc` a RESIDUAL
    test beside the step test, and at iteration 8 the residual is
    5.3e-7 A, 4.7x over `StandardNewton`'s tolerance for that row --
    the old criterion never looked at it.  Iteration 9 brings it to
    3e-11.  The old solver's 14 is unchanged because its criterion is
    the one it always had.)
    """
    xo, _vo, io = _old_solve_dc(_charge_diode(vsrc=5.0), gnd)
    xn, _vn, i_n = P.solve_dc(_charge_diode(vsrc=5.0), gnd)
    assert io == 14 and i_n == 9, (io, i_n)
    assert_allclose(xn, xo, rtol=0, atol=1e-8)


def test_a_participant_is_excluded_from_the_assembly_not_subtracted():
    """The mechanism that replaced the subtraction, on a scalar diode:
    during `augmented_system` the participant's own `i` and `G` are not
    called, and afterwards the element is exactly as it was (no instance
    attribute left behind)."""
    c = _fig1()
    d1 = c.elements['D1']
    calls = []
    real_i, real_G = Diode.i, Diode.G
    Diode.i = lambda self, *a, **k: calls.append('i') or real_i(self, *a, **k)
    Diode.G = lambda self, *a, **k: calls.append('G') or real_G(self, *a, **k)
    try:
        pairs = P.pcnr_junctions(c)
        P.augmented_system(c, np.zeros(c.n), np.zeros(len(pairs)), pairs,
                           dense_blocks=False)
    finally:
        Diode.i, Diode.G = real_i, real_G
    assert calls == [], calls
    assert 'i' not in d1.__dict__ and 'G' not in d1.__dict__
    ## and the ordinary path is untouched: the class methods answer again.
    assert np.asarray(c.i(np.zeros(c.n)), float).shape == (c.n,)


## ======================================================================
## 2.  The vector protocol on the Gummel-Poon BJT.
## ======================================================================

def _bjt(**card):
    el = eh.GummelPoonNpnHdl('c', 'b', 'e', **card)
    el.update_iparv()
    return el


def _params(el):
    return {p.name: getattr(el.iparv, p.name) for p in el.instparams}


def test_the_default_card_bjt_declares_two_pnj_probes_over_three_rows():
    el = _bjt()
    assert el.pcnr_probes == (((1, 2), 'pnj'), ((1, 0), 'pnj'))
    assert not hasattr(el, 'pcnr_junctions')      # NOT the scalar protocol
    assert 'PCNR: vector route, 2 probes over 3 rows (pnjlim on (b,e), ' \
           'pnjlim on (b,c))' in explain(el, source=False, symbolic=False)
    ## The p-n-p is the same model with its branches reversed.
    pnp = eh.GummelPoonPnpHdl('c', 'b', 'e')
    pnp.update_iparv()
    assert pnp.pcnr_probes == (((2, 1), 'pnj'), ((0, 1), 'pnj'))


@pytest.mark.parametrize('card', [dict(), dict(NPN)], ids=['default', 'NPN'])
def test_pcnr_i_is_the_resistive_current_and_didv_its_derivative(card):
    """At ``v`` equal to the branch voltages, ``pcnr_i(v)`` IS the
    element's ``i(x)`` (there is no charge at DC), to the bit; and the
    ``3 x 2`` block is its derivative by finite differences."""
    el = _bjt(**card)
    pr = _params(el)
    for x in (np.array([2.0, 0.7, 0.0]), np.array([0.3, 0.75, 0.0]),
              np.array([5.0, -0.2, 0.0])):
        v = np.array([x[1] - x[2], x[1] - x[0]])
        i_el = np.asarray(el.i(x), float)
        i_pc = np.asarray(el.pcnr_i(v, pr, defaultepar, numeric), float)
        assert i_pc.shape == (3,)
        assert np.array_equal(i_el, i_pc), (x, i_el, i_pc)
        blk = np.asarray(el.pcnr_didv(v, pr, defaultepar, numeric), float)
        assert blk.shape == (3, 2)
        eps = 1e-7
        fd = np.zeros_like(blk)
        for j in range(2):
            dv = np.zeros(2)
            dv[j] = eps
            fd[:, j] = (np.asarray(el.pcnr_i(v + dv, pr, defaultepar,
                                             numeric), float)
                        - np.asarray(el.pcnr_i(v - dv, pr, defaultepar,
                                               numeric), float)) / (2 * eps)
        scale = max(np.max(np.abs(blk)), 1e-30)
        assert_allclose(blk, fd, rtol=1e-5, atol=1e-7 * scale)
        ## Kirchhoff: the three terminal currents sum to zero.
        assert abs(np.sum(i_pc)) <= 1e-12 * max(np.max(np.abs(i_pc)), 1e-30)


def test_a_dropped_didv_column_is_caught_by_the_finite_difference_check():
    """Mutation control for the test above: a block missing its second
    column disagrees with FD by the whole ``d/dvbc`` column, which at
    this bias is the transport current's slope -- not a rounding-level
    disagreement."""
    el = _bjt(**NPN)
    pr = _params(el)
    v = np.array([0.75, 0.75 - 0.3])            # saturated: both slopes big
    blk = np.asarray(el.pcnr_didv(v, pr, defaultepar, numeric), float)
    bad = blk.copy()
    bad[:, 1] = 0.0
    eps = 1e-7
    dv = np.array([0.0, eps])
    fd1 = (np.asarray(el.pcnr_i(v + dv, pr, defaultepar, numeric), float)
           - np.asarray(el.pcnr_i(v - dv, pr, defaultepar, numeric),
                        float)) / (2 * eps)
    assert not np.allclose(bad[:, 1], fd1, rtol=1e-5,
                           atol=1e-7 * np.max(np.abs(blk)))
    assert_allclose(blk[:, 1], fd1, rtol=1e-5, atol=1e-7 * np.max(np.abs(blk)))


def test_pcnr_limit_is_the_declared_pnjlim_per_probe_and_nothing_is_read_from_x():
    """On the default card neither limiter parameter reads the solution
    (`isT` and `nf*vtT` are parameter-and-temperature only), so the
    block is limited probe by probe with SPICE's `pnjlim` on the
    device's own ``IS``/``VT`` -- and `x_old_sub` is accepted and
    ignored.  The paper's independence: each probe from its own
    ``vold``."""
    el = _bjt()
    pr = _params(el)
    ((_r0, k0, _m0, p0), (_r1, k1, _m1, p1)) = el._hdl_info['limit_spec']
    assert k0 == k1 == 'pnj'
    assert not any(getattr(f, '_wants_x', False) for f in p0 + p1)
    args = [pr[q] for q in el._hdl_paramnames] + [defaultepar.T]
    IS0, VT0 = [float(f(*args)) for f in p0]
    IS1, VT1 = [float(f(*args)) for f in p1]
    v_new, v_old = np.array([5.0, 0.3]), np.array([0.7, 0.2])
    want = np.array([_pnjlim(5.0, 0.7, VT0, IS0, numeric),
                     _pnjlim(0.3, 0.2, VT1, IS1, numeric)])
    got = el.pcnr_limit(v_new, v_old, pr, defaultepar, numeric,
                        np.array([9.0, 9.0, 9.0]))
    assert np.array_equal(got, want), (got, want)
    assert got[0] < 1.0 and got[1] == 0.3           # one bit, one passed


def test_a_limiter_parameter_that_reads_the_solution_sees_x_old_sub():
    """`_wants_x` through the vector route: the parameter is evaluated at
    the LAST ACCEPTED sub-vector handed to `pcnr_limit`, SPICE's `von`
    semantics -- so two different `x_old` give two different limits."""
    class _Dx(Behavioural):
        instparams = [Parameter(name='IS', desc='Is', unit='A',
                                default=1e-14)]

        @staticmethod
        def analog(p, n):
            b = Branch(p, n)
            vtx = var(vt() * (1.0 + 0.5 * b.V), 'vtx')
            v = var(limit_pnj(b.V, IS, vtx), 'vlim')             # noqa
            return Contribution(b.I, IS * (limexp(v / vt()) - 1))  # noqa

    el = _Dx('p', 'n')
    el.update_iparv()
    assert '[params at last iterate]' in explain(el, source=False,
                                                 symbolic=False)
    assert el.pcnr_probes == (((0, 1), 'pnj'),)
    pr = {'IS': 1e-14}
    a = el.pcnr_limit(np.array([5.0]), np.array([0.7]), pr, defaultepar,
                      numeric, np.array([0.0, 0.0]))
    b = el.pcnr_limit(np.array([5.0]), np.array([0.7]), pr, defaultepar,
                      numeric, np.array([1.0, 0.0]))
    assert a[0] != b[0] and 0.5 < a[0] < 1.0 and 0.5 < b[0] < 1.0
    ## and `refine` routes the circuit's `x_old` to it.
    c = SubCircuit()
    c['vs'] = VS('a', gnd, v=1.0)
    c['R'] = R('a', 'b', r=1e3)
    c['D'] = _Dx('b', gnd)
    c.update_iparv()
    devs = P.pcnr_devices(c)
    x0 = np.zeros(c.n)
    x1 = np.zeros(c.n)
    x1[c.get_node_index('b')] = 1.0
    r0 = P.refine(devs, np.array([0.7]), np.array([5.0]), defaultepar, x_old=x0)
    r1 = P.refine(devs, np.array([0.7]), np.array([5.0]), defaultepar, x_old=x1)
    assert r0[0] == a[0] and r1[0] == b[0]


## ======================================================================
## 3.  Refusals, and the gmin pair view.
## ======================================================================

def test_the_refusals_name_their_rule():
    from pycircuit.circuit.tests.test_limit_fet import _fet
    ## ⚠ THIS ASSERTION EXPIRED, 2026-08-27, and is inverted rather than
    ## deleted.  It used to read:
    ##
    ##   assert 'PCNR: does not qualify -- I(b) reads b, bi, which no
    ##           $limit probe limits' in txt
    ##
    ## and it was the whole reason PCNR could not take a BJT with a base
    ## resistance -- which is to say, any real BJT card.  The base
    ## resistance branch is now declared as a `limit_identity` probe
    ## (roadmap sec. 42), so its current reaches the solution through a
    ## declared unknown and the rule is satisfied rather than waived.
    el = _bjt(rb=10.0)
    txt = explain(el, source=False, symbolic=False)
    assert 'PCNR: does not qualify' not in txt, txt
    assert [k for _p, k in el.pcnr_probes].count('id') == 1, el.pcnr_probes
    ## ⚠ EXPIRED 2026-08-27 (roadmap sec. 44), inverted not deleted.  It
    ## read:
    ##
    ##   assert 'PCNR: does not qualify -- a branch-current unknown'
    ##          in explain(th, ...)
    ##
    ## A V-contributed branch no longer refuses the device.  What still
    ## does is a GENERATED STATE -- device memory -- and that is what
    ## this now pins, because it is the half of the old rule that
    ## survived.
    th = eh.GummelPoonNpnThermalHdl('c', 'b', 'e', 'th', 'tha')
    th.update_iparv()
    assert 'does not qualify' not in explain(th, source=False,
                                             symbolic=False)
    idt = eh.IdtHdl(*['p%d' % i for i in range(len(eh.IdtHdl.terminals))])
    idt.update_iparv()
    assert 'PCNR: does not qualify -- a generated state' in explain(
        idt, source=False, symbolic=False)
    ## ⚠ EXPIRED 2026-08-27, inverted rather than deleted.  It read:
    ##
    ##   assert 'var(dT) reads th, tha, which no $limit probe limits'
    ##          in explain(th, ...)
    ##
    ## `SelfHeating` now declares the thermal branch as a
    ## `limit_identity` probe (roadmap sec. 43), so `dT` is a limited
    ## unknown and every temperature-dependent current -- which is all
    ## of them -- reaches the solution through declared probes.
    th = eh.GummelPoonNpnThermalHdl('c', 'b', 'e', 'th', 'tha', rth=100.0)
    th.update_iparv()
    assert 'does not qualify' not in explain(th, source=False,
                                             symbolic=False)
    assert [k for _p, k in th.pcnr_probes].count('id') >= 1, th.pcnr_probes
    ## Stage 2: a `fet`/`vds` probe is a PCNR unknown like any other.
    f = _fet('both')('d', 'g', 's')
    f.update_iparv()
    assert 'PCNR: vector route, 2 probes over 3 rows (fetlim on (g,s), ' \
           'limvds on (d,s))' in explain(f, source=False, symbolic=False)
    assert hasattr(f, 'pcnr_i')
    ## EKV used to read its bulk-source potential RAW and was refused by
    ## the unlimited-branch rule ("var(vsb) reads b, g, which no $limit
    ## probe limits").  It has declared a probe on that branch since
    ## 2026-08-26 and qualifies with three; the refusal is still pinned,
    ## on the raw twin, in `test_limit_identity.py`.  That probe was
    ## `limit_identity` for a day and is `limit_delta` now (sec. 34) --
    ## what qualifies the device is having a probe there at all, and
    ## what the identity one could not do was limit the seed.
    ekv = eh.EkvNmosHdl('d', 'g', 's', 'b')
    ekv.update_iparv()
    assert 'PCNR: vector route, 3 probes over 4 rows (deltalim ' \
           'on (s,b), fetlim on (g,s), limvds on (d,s))' \
           in explain(ekv, source=False, symbolic=False)
    assert hasattr(ekv, 'pcnr_i')


def test_the_gmin_pair_view_lists_pnj_probes_only_in_v_lim_order():
    """`pcnr_junctions` is what `dcanalysis._jrows` and the JAX ladder
    put a gmin across, on the ordinary path too.  A `fet`/`vds` probe
    must never reach it (a gmin across `vgs` is a gate leak); the pairs
    it does list are in flat ``v_lim`` order so that a pair-shaped
    consumer seeds ``v_lim`` correctly."""
    c = SubCircuit()
    c['vcc'] = VS('vcc', gnd, v=5.0)
    c['rin'] = R('vcc', 'nb', r=4.3e3)
    c['rl'] = R('vcc', 'no', r=2e3)
    c['q1'] = eh.GummelPoonNpnHdl('nb', 'nb', gnd)
    c['q2'] = eh.GummelPoonNpnHdl('no', 'nb', gnd)
    c['d'] = Diode('no', gnd, IS=1e-15)
    devs = P.pcnr_devices(c)
    assert [d.instance for d in devs] == ['q1', 'q2', 'd']
    assert [d.off for d in devs] == [0, 2, 4]
    assert [d.scalar for d in devs] == [False, False, True]
    pairs = P.pcnr_junctions(c)
    assert [(i, ra, rb) for i, _e, ra, rb in pairs] == \
        [(d.instance, ra, rb) for d in devs for ra, rb in d.probes]
    assert P.flat_probes(devs) == [(ra, rb) for _i, _e, ra, rb in pairs]
    assert P.pcnr_junction_pairs is P.pcnr_junctions

    ## A synthetic vector device with a non-pnj probe: kept by
    ## `pcnr_devices`, EXCLUDED from the pair view.
    class _Fake(object):
        pcnr_probes = (((1, 2), 'fet'), ((0, 2), 'vds'), ((3, 2), 'pnj'))

    dev = P._device_of('m', _Fake(), [10, 11, 12, 13])
    assert dev.kinds == ('fet', 'vds', 'pnj') and dev.m == 3
    only = [(ra, rb) for (ra, rb), k in zip(dev.probes, dev.kinds)
            if k == 'pnj']
    assert only == [(13, 12)]
    ## Mutation control: the view built WITHOUT the tag filter would list
    ## the gate probe (11, 12), which is precisely what must not happen.
    assert (11, 12) in dev.probes and (11, 12) not in only


def test_the_ordinary_dc_path_still_solves_bjt_circuits_with_the_new_gmin_rows():
    """`dcanalysis._jrows` now receives the BJT's two junctions as gmin
    targets on the ordinary path.  The library's DC tests cover the
    answer; this pins that a diode-connected device (two rows equal) is
    handled and that the ladder is not tripped on a plain solve."""
    c = SubCircuit()
    c['vin'] = VS('nvi', gnd, v=5.0)
    c['rin'] = R('nvi', 'nb', r=4.25e3)
    c['vo'] = VS('no', gnd, v=3.0)
    c['q1'] = eh.GummelPoonNpnHdl('nb', 'nb', gnd, **NPN_IDEAL)
    c['q2'] = eh.GummelPoonNpnHdl('no', 'nb', gnd, **NPN_IDEAL)
    res = DC(c, toolkit=numeric).solve()
    assert 0.6 < float(res.v('nb')) < 0.85


## ======================================================================
## 4.  Same answer, with and without.
## ======================================================================

def _mirror(q2_first=False, card=None, load=2e3, vout=None):
    """Diode-connected q1 driving q2 from one base node: 5 V through 4.25 k
    into the base, and either a resistive load to the supply or a
    voltage source holding q2's collector."""
    card = dict(NPN_IDEAL, bf=80.0, vaf=60.0) if card is None else card
    c = SubCircuit()
    c['vin'] = VS('nvi', gnd, v=5.0)
    c['rin'] = R('nvi', 'nb', r=4.25e3)
    if vout is None:
        c['rl'] = R('nvi', 'no', r=load)
    else:
        c['vo'] = VS('no', gnd, v=vout)
    devs = [('q1', ('nb', 'nb', gnd)), ('q2', ('no', 'nb', gnd))]
    for nm, nodes in (devs[::-1] if q2_first else devs):
        c[nm] = eh.GummelPoonNpnHdl(*nodes, **card)
    return c


def _ce(card=None, vin=1.0):
    """The common-emitter stage of `test_elements_hdl_library3`'s
    operating-point test, at ``rb = re = rc = 0`` so that it qualifies."""
    c = SubCircuit()
    c['vcc'] = VS('nvcc', gnd, v=5.0)
    c['vin'] = VS('nvi', gnd, v=vin)
    c['rb'] = R('nvi', 'nb', r=47e3)
    c['rc'] = R('nvcc', 'nc', r=4.7e3)
    c['q1'] = eh.GummelPoonNpnHdl('nc', 'nb', gnd, **(card or {}))
    return c


SAME_ANSWER = [
    ('mirror default card, VS output', lambda: _mirror(card={}, vout=3.0)),
    ('mirror ideal card, R load', lambda: _mirror()),
    ('mirror NPN card, R load', lambda: _mirror(card=dict(NPN))),
    ('CE default card', lambda: _ce()),
    ('CE NPN card', lambda: _ce(card=dict(NPN))),
]


@pytest.mark.parametrize('name, mk', SAME_ANSWER, ids=[n for n, _ in SAME_ANSWER])
def test_pcnr_and_the_ordinary_path_find_the_same_operating_point(name, mk):
    """To 1e-9 -- with BOTH solvers at ``reltol = 1e-9``.  At the default
    1e-6 the two stop 2.1e-9 apart on the first circuit, which is inside
    their own stopping band and says nothing about the equations."""
    off = np.asarray(DC(mk(), toolkit=numeric, reltol=1e-9).solve().x, float)
    on = np.asarray(DC(mk(), toolkit=numeric, pcnr=True,
                       reltol=1e-9).solve().x, float)
    x, v_lim, _its = P.solve_dc(mk(), gnd, reltol=1e-9)
    assert_allclose(on, off, rtol=0, atol=1e-9)
    assert_allclose(x, off, rtol=0, atol=1e-9)
    ## The limited unknowns agree with the branch voltages they stand for.
    c = mk()
    for dev, k in zip(P.pcnr_devices(c), range(10)):
        for j, (ra, rb) in enumerate(dev.probes):
            assert abs(v_lim[dev.off + j] - (x[ra] - x[rb])) < 1e-9
    ## and it is a real bias point, not cutoff.
    assert 0.6 < float(off[c.get_node_index('nb')]) < 0.85


def test_the_o_k_schur_trick_survives_vector_unknowns():
    """The corrected premise: with ``t = 3, m = 2`` the sparse rank-one
    form equals the dense ``J_ml @ J_lm`` -- only the column widened."""
    c = _mirror(card=dict(NPN))
    devs = P.pcnr_devices(c)
    x = np.zeros(c.n)
    x[c.get_node_index('nvi')] = 5.0
    x[c.get_node_index('nb')] = 0.7
    x[c.get_node_index('no')] = 2.0
    v = P.v_lim_init(devs, x) + np.array([0.01, -0.02, 0.03, -0.04])
    u = np.asarray(c.u(0.0, defaultepar, analysis='dc'), dtype=float)
    g_mna, g_lim, J_mm, J_ml, J_lm, didv = P.augmented_system(
        c, x, v, devs, defaultepar, u_extra=u, dense_blocks=True)
    assert J_ml.shape == (c.n, 4) and J_lm.shape == (4, c.n)
    ## a BJT's column has THREE nonzeros (its rows), a row still two.  q2's
    ## column: q1 is diode-connected, so two of its rows are ONE row.
    assert np.count_nonzero(J_ml[:, 2]) == 3
    assert np.count_nonzero(J_ml[:, 0]) == 2
    assert all(np.count_nonzero(J_lm[k]) in (0, 2) for k in range(4))
    f_s, J_s = P.schur_reduce(g_mna, g_lim, J_mm, junctions=devs, didv=didv)
    f_d, J_d = P.schur_reduce(g_mna, g_lim, J_mm, J_ml, J_lm)
    assert_allclose(J_s, J_d, rtol=1e-12, atol=1e-12 * np.max(np.abs(J_d)))
    assert_allclose(f_s, f_d, rtol=1e-12, atol=1e-12 * np.max(np.abs(f_d)))
    dx = np.linspace(-0.1, 0.1, c.n)
    assert_allclose(P.dx_lim_of(devs, g_lim, dx), -(g_lim + J_lm @ dx),
                    rtol=0, atol=1e-15)


## ======================================================================
## 5.  The intervention against its absence: plain Newton, gmin = 0,
##     both instance orders, several starts.
## ======================================================================

def _plain(c, x0, reltol=1e-4, vabstol=1e-6, iabstol=1e-12, maxiter=100):
    """`StandardNewton` with the ordinary device limiting and NO rescue
    ladder, from ``x0``; ``(iterations, x)`` or ``(None, None)``."""
    iref = c.get_node_index(gnd)
    epar = defaultepar

    def rfunc(xr):
        x = np.insert(xr, iref, 0.0)
        F = c.i(x, epar) + c.u(0, analysis='dc', epar=epar)
        J = c.G(x, epar)
        return np.delete(F, iref), np.delete(np.delete(J, iref, 0), iref, 1)

    def limiter(xr, x0r):
        return np.delete(c.limit(np.insert(xr, iref, 0.0),
                                 np.insert(x0r, iref, 0.0), epar), iref)

    nn, nb = len(c.nodes), len(c.branches)
    abstol = np.delete(np.concatenate((iabstol * np.ones(nn),
                                       vabstol * np.ones(nb))), iref)
    xtol = np.delete(np.concatenate((vabstol * np.ones(nn),
                                     iabstol * np.ones(nb))), iref)
    try:
        x, its = StandardNewton().solve_system(
            np.delete(np.asarray(x0, float), iref), rfunc, numeric, reltol,
            abstol, xtol, maxiter, limiter=limiter)
    except (NoConvergenceError, SingularMatrix, FloatingPointError,
            OverflowError, ValueError):
        return None, None
    return its, np.insert(x, iref, 0.0)


def _pcnr(c, x0):
    try:
        x, _v, its = P.solve_dc(c, gnd, x0=x0, reltol=1e-4, abstol=1e-6)
    except (RuntimeError, np.linalg.LinAlgError, FloatingPointError,
            OverflowError, ValueError):
        return None, None
    return its, x


def _rows(mk, start):
    """`([plain q1-first, plain q2-first], [pcnr q1-first, pcnr q2-first])`
    as ``(iterations, max|x - reference|)``; a failure is ``(None, None)``."""
    ref = np.asarray(DC(mk(), toolkit=numeric).solve().x, float)
    plain, pcnr = [], []
    for q2_first in (False, True):
        for solver, out in ((_plain, plain), (_pcnr, pcnr)):
            c = mk(q2_first)
            x0 = np.full(c.n, start)
            x0[c.get_node_index(gnd)] = 0.0
            its, x = solver(c, x0)
            out.append((its, None if x is None
                        else float(np.max(np.abs(x - ref)))))
    return plain, pcnr


def test_the_stage_0_bjt_mirror_under_pcnr_is_order_independent_at_every_start():
    """THE MEASUREMENT, full NPN card at ``rb = 0``, resistive load, plain
    Newton at ``gmin = 0`` with the ordinary limiter against `pcnr.solve_dc`
    (which has no ladder and no anchor either), `[q1 first, q2 first]`,
    Jacobian evaluations::

        start   plain          pcnr          pcnr before the
                                             limited seed (2026-08-26)
         0 V    [9, 9]         [9, 9]        [9, 9]
        -5 V    [9, 9]         [9, 9]        [9, 9]
        +5 V    [FAIL, FAIL]   [8, 8]        [164, 164]
       +10 V    [7, 7]         [8, 8]        [7, 7]
       +20 V    [FAIL, 7]      [8, 8]        [FAIL, FAIL]

    **Claim 3 has expired, as the last paragraph below anticipated.**
    `v_lim_init` used to seed each unknown with the RAW branch voltage, so
    a uniform 20 V start seeded ``vbe = 20 V`` and the first Jacobian was
    built at ``exp(20/vt)``: ``cond = 4.6e94``, first step 4.5e45 V.  The
    seed is now passed through the device's own limiter (`solve_dc` only
    -- the transient seeds from a CONVERGED point and still uses the raw
    value).  +20 V converges in 8, +5 V costs 164 -> 8, and +10 V pays one
    iteration for it.  The diagnosis recorded below was correct and is
    kept: it is why the seed, not the step, was the thing to fix.

    The +10 V PCNR entry was 6 until Stage 2 gave `solve_dc` a residual
    test: at iteration 6 the base row's residual is 1.3e-5 A, 4.2x over
    the tolerance `StandardNewton` applies to the same row, and the old
    step-only criterion accepted it.  Every other entry, +5 V's 164
    included, is unchanged by the stricter test.

    Three claims, each asserted below:

    1. **PCNR is order-independent at every start** -- equal counts, or
       the same failure, in both orders.  That is the structural claim
       (each device limits its OWN unknowns; nothing is shared) and it
       holds without exception.
    2. **Plain Newton is order-dependent at +20 V**: `[FAIL, 7]`.  Traced,
       the "win" is an accident.  q2's base-collector probe maps its 9e45
       V proposal to 2.83 V (SPICE's `pnjlim`, ``vold <= 0`` branch,
       ``VT*log(vnew/VT)``); visited FIRST, q2 writes that 2.83 V relative
       to the still-UNLIMITED base (16.5 V) and then moves the base to
       0.85 V, leaving the collector at 13.7 V -- reverse-biased by a
       stale anchor.  Visited second, the same 2.83 V lands forward of
       the already-limited base: a 1e33 S junction, singular matrix.
    3. **PCNR fails identically in both orders at +20 V** -- and the
       design of record expected it to converge.  It reproduces the
       declared law honestly: ``vbc = 2.83 V`` in both orders, then the
       singular matrix.  Capping that branch at ``vcrit`` (tried as a
       diagnostic, not shipped) exposes the second, structural obstacle:
       PCNR leaves the NODES unlimited, so the accepted iterate carries
       ``no = -9.08e45``, and the +9.08e45 correction the next step
       computes has an ulp of ~1e30 -- the node cannot be brought back
       to 15 V by addition.  The ordinary path never carries such a node
       because its limiter writes the bounded value INTO it.  From +10
       and +40 V the capped law converged order-independently; +20 V is
       where the rounding lands badly.  A start of +20 V on every node
       is not a bias any circuit reaches; +5 V is, and there PCNR is the
       only one of the two that converges.

    The +20 V row is asserted as EQUALITY between the orders, so a
    later fix that makes both converge passes unchanged.
    """
    mk = lambda q=False: _mirror(q, card=dict(NPN))          # noqa: E731
    table = {s: _rows(mk, s) for s in (0.0, -5.0, 5.0, 10.0, 20.0)}
    for start, (plain, pcnr) in table.items():
        assert pcnr[0][0] == pcnr[1][0], (start, pcnr)      # claim 1
        for its, err in pcnr:
            if its is not None:
                assert err < 1e-5, (start, pcnr)
    assert table[0.0][1][0][0] == 9 and table[-5.0][1][0][0] == 9, table
    assert table[10.0][1][0][0] == 8 and table[10.0][0][0][0] == 7, table
    ## claim 2 -- the control, as measured (the winning order at 7).
    assert table[20.0][0][0][0] is None and table[20.0][0][1][0] == 7, table
    ## +5 V: plain fails both orders, PCNR converges both.
    assert table[5.0][0] == [(None, None), (None, None)], table[5.0]
    assert all(its is not None and its < 20 for its, _e in table[5.0][1]), \
        table[5.0]
    ## claim 3, INVERTED 2026-08-26 (the docstring anticipated this).  It
    ## used to read "PCNR fails identically in both orders at +20 V" and
    ## was asserted as equality precisely so a fix would pass it unchanged.
    ## Equality still holds and is still the structural claim; what is new
    ## is that both orders now CONVERGE, so assert that too -- otherwise
    ## this row would silently accept a regression back to [FAIL, FAIL].
    assert table[20.0][1][0] == table[20.0][1][1], table[20.0]
    assert all(its is not None and its < 20 for its, _e in table[20.0][1]), \
        table[20.0]


def test_the_papers_shape_two_bjts_with_is_four_decades_apart():
    """The paper's Fig. 1 with BJTs: two devices in parallel on one base
    branch and one collector node, ``IS`` 1e-16 and 1e-12.  Plain Newton
    shows no clash here (the shared-gate case composes, as the earlier
    search found); PCNR is order-independent and never slower."""
    def mk(q2_first=False):
        c = SubCircuit()
        c['vcc'] = VS('vcc', gnd, v=5.0)
        c['rb'] = R('vcc', 'nb', r=47e3)
        c['rl'] = R('vcc', 'nc', r=2e3)
        devs = [('q1', dict(IS=1e-16)), ('q2', dict(IS=1e-12))]
        for nm, card in (devs[::-1] if q2_first else devs):
            c[nm] = eh.GummelPoonNpnHdl('nc', 'nb', gnd,
                                        **dict(NPN_IDEAL, **card))
        return c
    for start in (0.0, -5.0, 20.0):
        plain, pcnr = _rows(mk, start)
        assert pcnr[0][0] == pcnr[1][0] and pcnr[0][0] is not None, \
            (start, pcnr)
        assert plain[0][0] == plain[1][0] and plain[0][0] is not None, \
            (start, plain)
        assert pcnr[0][0] <= plain[0][0], (start, plain, pcnr)
        for its, err in plain + pcnr:
            assert err < 1e-5


## ======================================================================
## 6.  Per-component convergence.
## ======================================================================

def test_convergence_is_judged_per_component():
    """`|g_lim,j| < reltol*max(|v_j|, 1) + abstol` for EVERY j.  The old
    criterion scaled every component by ``max|v_new|``, so a 40 V
    component (a `vds`, Stage 2) would have loosened a 0.7 V junction's
    residual forty-fold.  Mutation control: the old criterion accepts
    this vector; the new one refuses it."""
    g_lim = np.array([2e-5, 0.0])
    v_new = np.array([0.7, 40.0])
    reltol, abstol = 1e-6, 1e-12
    old = np.max(np.abs(g_lim)) < reltol * max(np.max(np.abs(v_new)), 1.0) \
        + abstol
    assert old is True or old == True                     # noqa: E712
    assert P.lim_converged(g_lim, v_new, reltol, abstol) is False
    ## and where every component is one scalar junction the two agree.
    assert P.lim_converged(np.array([5e-7]), np.array([0.7]), reltol, abstol)
    assert not P.lim_converged(np.array([2e-6]), np.array([0.7]), reltol,
                               abstol)
    assert P.lim_converged(np.array([2e-6, 0.0]), np.array([40.0, 0.7]),
                           reltol, abstol)


def test_solve_dc_uses_the_per_component_criterion():
    """`solve_dc` stops on the shared function, not a private formula: a
    spy that refuses once forces one more iteration on the BJT stage."""
    calls = []
    real = P.lim_converged

    def spy(g, v, r, a):
        calls.append(np.array(g, copy=True))
        return real(g, v, r, a)

    P.lim_converged = spy
    try:
        _x, _v, its = P.solve_dc(_ce(), gnd)
    finally:
        P.lim_converged = real
    ## Behind a short-circuit `and`: consulted only once the MNA step has
    ## converged, and the last consultation is the one that stopped it.
    assert 1 <= len(calls) <= its and its > 3
    assert np.all(np.abs(calls[-1]) < 1e-6 + 1e-12)


## ======================================================================
## 7.  Stage 2: FET probes (`fetlim`/`limvds`) and a redundant probe.
## ======================================================================

def _aug_fd(c, x, v, devs, eps=1e-6):
    """Finite-difference ``d[g_mna; g_lim]/d[x; v]`` against the assembled
    ``[[J_mm, J_ml], [J_lm, I]]``.  Returns (assembled, fd)."""
    u = np.asarray(c.u(0.0, defaultepar, analysis='dc'), dtype=float)

    def g(xx, vv):
        gm, gl, *_r = P.augmented_system(c, xx, vv, devs, defaultepar,
                                         u_extra=u, dense_blocks=False)
        return np.concatenate((gm, gl))

    gm, gl, J_mm, J_ml, J_lm, _d = P.augmented_system(
        c, x, v, devs, defaultepar, u_extra=u, dense_blocks=True)
    k = len(v)
    J = np.block([[J_mm, J_ml], [J_lm, np.eye(k)]])
    z = np.concatenate((x, v))
    fd = np.zeros_like(J)
    for j in range(len(z)):
        dz = np.zeros(len(z))
        dz[j] = eps
        fd[:, j] = (g(*np.split(z + dz, [c.n])) - g(*np.split(z - dz, [c.n]))) \
            / (2 * eps)
    return J, fd


def _mos4_cascode():
    from pycircuit.circuit.tests.test_device_limiter import _mos4, _cascode4
    return _cascode4(_mos4('group'), 20.0, 2.0, 0.8)


def _l1_cascode(vb=-0.5, **card):
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=20.0)
    c['vg2'] = VS('g2', gnd, v=2.0)
    c['vg1'] = VS('g1', gnd, v=0.8)
    c['vb'] = VS('bulk', gnd, v=vb)
    c['M2'] = eh.MosLevel1Hdl('vdd', 'g2', 'mid', 'bulk', **card)
    c['M1'] = eh.MosLevel1Hdl('mid', 'g1', gnd, 'bulk', **card)
    return c


@pytest.mark.parametrize('mk', [_mos4_cascode, _l1_cascode],
                         ids=['mos4-group', 'level1'])
def test_the_augmented_jacobian_is_the_derivative_with_fet_probes(mk):
    """``J == d[g_mna; g_lim]/d[x; v]`` by central differences on a
    cascode of four-probe MOSFETs -- three unknowns each (`fetlim` on
    `vgs`, `limvds` on `vds`, `pnjlim` on `vbs`), the redundant `(b,d)`
    read off the tree.  Level 1 carries its charge; at DC it is inert
    and the block is the resistive derivative."""
    c = mk()
    devs = P.pcnr_devices(c)
    assert [d.kinds for d in devs] == [('fet', 'vds', 'pnj')] * 2
    x = np.zeros(c.n)
    for nm, val in (('vdd', 20.0), ('g2', 2.0), ('g1', 0.8), ('bulk', -0.5),
                    ('mid', 1.1)):
        x[c.get_node_index(nm)] = val
    v = P.v_lim_init(devs, x) + np.array([0.02, -0.1, 0.01, -0.03, 0.2, 0.0])
    J, fd = _aug_fd(c, x, v, devs)
    scale = np.max(np.abs(J))
    assert_allclose(J, fd, rtol=1e-5, atol=2e-6 * scale)
    ## Mutation control: a `fet` probe's `didv` column dropped from J_ml
    ## is caught -- that column is the drain current's slope in `vgs`.
    bad = J.copy()
    bad[:c.n, c.n] = 0.0
    assert not np.allclose(bad, fd, rtol=1e-5, atol=2e-6 * scale)


def test_the_cycle_is_carried_as_a_tree_and_the_redundant_probe_is_named():
    """`_mos4` declares `(g,s)`, `(d,s)`, `(b,s)`, `(b,d)`; the last closes
    a triangle.  The device's unknowns are the first three, `(b,d)` is
    recorded with its combination `+(b,s) -(d,s)` and its `pnjlim`, and
    `explain()` says so.  At convergence every unknown equals its branch
    voltage and the implied `vbd` equals the branch it stands for."""
    from pycircuit.circuit.tests.test_device_limiter import _mos4
    el = _mos4('group')('d', 'g', 's', 'b')
    el.update_iparv()
    assert el.pcnr_probes == (((1, 2), 'fet'), ((0, 2), 'vds'),
                              ((3, 2), 'pnj'))
    assert el.pcnr_redundant == (((3, 0), 'pnj', (0, -1, 1)),)
    txt = explain(el, source=False, symbolic=False)
    assert 'PCNR: vector route, 3 probes over 4 rows (fetlim on (g,s), ' \
           'limvds on (d,s), pnjlim on (b,s))' in txt
    assert 'PCNR redundant probe: pnjlim on (b,d) = -(d,s) +(b,s)' in txt
    c = _mos4_cascode()
    x, v_lim, _its = P.solve_dc(c, gnd, reltol=1e-9, abstol=1e-12)
    for dev in P.pcnr_devices(c):
        for j, (ra, rb) in enumerate(dev.probes):
            assert abs(v_lim[dev.off + j] - (x[ra] - x[rb])) < 1e-12
        for (la, lb), _k, coeffs in dev.redundant:
            implied = float(np.dot(coeffs, v_lim[dev.off:dev.off + dev.m]))
            assert abs(implied - (x[dev.rows[la]] - x[dev.rows[lb]])) < 1e-12


def test_limit_block_is_the_per_probe_loop_without_a_redundant_probe():
    """Bit for bit: a tree of probes is limited probe by probe, and a
    block no law touched is returned as it came."""
    probes = [(1, 2), (0, 2)]
    ident = [[1, 0], [0, 1]]
    laws = [lambda vn, vo: min(vn, vo + 1.0), lambda vn, vo: vn]
    out = P.limit_block(probes, ident, [5.0, 3.3], [0.5, 3.0], laws)
    assert list(out) == [1.5, 3.3]
    same = P.limit_block(probes, ident, [0.7, 3.3], [0.5, 3.0], laws)
    assert list(same) == [0.7, 3.3]


def test_limit_block_resolves_a_cycle_by_dropping_the_smallest_correction():
    """Three laws on two unknowns (`vbd = vbs - vds`): the kept pair is
    the two asking for MORE, the third follows from them.  Which probe
    was declared redundant changes nothing about the resolved branch
    voltages -- the rule is a function of the data (`limit_together`'s
    write-back rule, for a tree of branch voltages)."""
    clamp = lambda hi: (lambda vn, vo: min(vn, hi))          # noqa: E731
    floor = lambda lo: (lambda vn, vo: max(vn, lo))          # noqa: E731
    keep = lambda vn, vo: vn                                 # noqa: E731
    ## unknowns (vds, vbs); redundant vbd = vbs - vds.  Raw: vds = 10,
    ## vbs = 0.9, vbd = -9.1.  limvds-like clamp on vds to 4 (corr 6),
    ## pnj-like clamp on vbs to 0.8 (corr 0.1), vbd's law: no bite.
    a = P.limit_block([(0, 2), (3, 2), (3, 0)], [[1, 0], [0, 1], [-1, 1]],
                      [10.0, 0.9], [2.0, 0.7], [clamp(4.0), clamp(0.8), keep])
    ## vds = 4 and vbs = 0.8 are both honoured; vbd = -3.2 follows.
    assert list(a) == [4.0, 0.8]
    ## Now let vbd's law bite HARDER than vbs's: it asks for -5 (corr 4.1)
    ## against vbs's 0.1.  vbs is the one dropped: vbs = vbd + vds = -1.
    b = P.limit_block([(0, 2), (3, 2), (3, 0)], [[1, 0], [0, 1], [-1, 1]],
                      [10.0, 0.9], [2.0, 0.7],
                      [clamp(4.0), clamp(0.8), floor(-5.0)])
    assert_allclose(b, [4.0, -1.0], rtol=0, atol=1e-15)
    ## Same laws, `vbd` carried as the unknown and `vbs` redundant:
    ## unknowns (vds, vbd), vbs = vbd + vds.  The four branch voltages
    ## resolve identically.
    b2 = P.limit_block([(0, 2), (3, 0), (3, 2)], [[1, 0], [0, 1], [1, 1]],
                       [10.0, -9.1], [2.0, -1.3],
                       [clamp(4.0), floor(-5.0), clamp(0.8)])
    assert_allclose([b2[0], b2[1] + b2[0]], b, rtol=0, atol=1e-15)


def test_von_reads_x_old_sub_under_pcnr_on_level_1():
    """SPICE's body-biased `von` through the vector route: `fetlim`'s
    threshold is evaluated at the device's LAST ACCEPTED sub-vector, so
    two `x_old` bulk biases clamp the same `vgs` proposal to two
    different values -- and the clamp moves the way the body effect
    says (a more negative `vbs` raises `von`)."""
    el = eh.MosLevel1Hdl('d', 'g', 's', 'b', gamma=0.9)
    el.update_iparv()
    pr = {p.name: getattr(el.iparv, p.name) for p in el.instparams}
    pr['$given:gamma'] = 1.0
    v_new = np.array([6.0, 1.0, -0.5])     # vgs proposal far above von
    v_old = np.array([0.0, 1.0, -0.5])
    x_a = np.array([1.0, 0.0, 0.0, 0.0])   # d, g, s, b: vbs = 0
    x_b = np.array([1.0, 0.0, 0.0, -3.0])  # vbs = -3
    a = el.pcnr_limit(v_new, v_old, pr, defaultepar, numeric, x_a)
    b = el.pcnr_limit(v_new, v_old, pr, defaultepar, numeric, x_b)
    ## `fetlim` off-and-turning-on lands at `von + 0.5`.
    assert a[0] != b[0] and b[0] > a[0] + 0.5, (a, b)
    assert a[1] == b[1] == 1.0 and a[2] == b[2] == -0.5
    ## and `refine` routes the circuit's `x_old` to it.  (`gamma` must be
    ## GIVEN: the default card's is 0, no body effect, and the two
    ## `x_old` then clamp identically -- which is not a defect, it is
    ## `von` with nothing to follow.)
    c = _l1_cascode(gamma=0.9)
    for dev in P.pcnr_devices(c):
        if dev.instance == 'M1':
            break
    ib = c.get_node_index('bulk')
    v0 = np.zeros(6)
    vn = np.zeros(6)
    vn[dev.off] = 6.0
    x0 = np.zeros(c.n)
    x1 = np.zeros(c.n)
    x1[ib] = -3.0
    r0 = P.refine(P.pcnr_devices(c), v0, vn, defaultepar, x_old=x0)
    r1 = P.refine(P.pcnr_devices(c), v0, vn, defaultepar, x_old=x1)
    assert r1[dev.off] > r0[dev.off] + 0.4, (r0, r1)


def test_the_gmin_pair_view_lists_only_the_junctions_of_a_mosfet():
    """With `_mos4` in the circuit the pair view holds the two bulk
    junctions per device -- `(b,s)` from the tree and the redundant
    `(b,d)` -- and never `(g,s)` or `(d,s)`: a gmin across `vgs` is a
    gate leak."""
    c = _mos4_cascode()
    pairs = P.pcnr_junctions(c)
    devs = {d.instance: d for d in P.pcnr_devices(c)}
    assert len(pairs) == 4
    for inst, _el, ra, rb in pairs:
        d = devs[inst]
        assert ra == d.rows[3] and rb in (d.rows[2], d.rows[0])
        assert (ra, rb) not in [(d.rows[1], d.rows[2])]     # never (g,s)
    ## Mutation control: the tree's `fet`/`vds` pairs exist and are not
    ## listed.
    all_tree = [(ra, rb) for d in devs.values() for ra, rb in d.probes]
    listed = [(ra, rb) for _i, _e, ra, rb in pairs]
    assert all(p not in listed for p, k in zip(all_tree, [k for d in
               devs.values() for k in d.kinds]) if k != 'pnj')


@pytest.mark.parametrize('which', ['mos4', 'fet'])
def test_pcnr_and_the_ordinary_path_agree_on_the_cascode(which):
    """To 1e-9 on the node voltages at `reltol = 1e-9` on both paths; the
    sources' branch currents to the same relative accuracy."""
    from pycircuit.circuit.tests.test_device_limiter import _mos4, _cascode4
    from pycircuit.circuit.tests.test_limit_fet import _fet, _cascode
    if which == 'mos4':
        mk = lambda: _cascode4(_mos4('group'), 20.0, 2.0, 0.8)   # noqa
    else:
        mk = lambda: _cascode(_fet('both'), 20.0, 2.0, 0.8)      # noqa
    c = mk()
    off = np.asarray(DC(c, toolkit=numeric, reltol=1e-9).solve().x, float)
    x, _v, _its = P.solve_dc(mk(), gnd, reltol=1e-9, abstol=1e-12)
    nn = len(c.nodes)
    assert_allclose(x[:nn], off[:nn], rtol=0, atol=1e-9)
    assert_allclose(x[nn:], off[nn:], rtol=1e-9, atol=1e-15)
    assert 0.5 < float(x[c.get_node_index('mid')]) < 2.0


def test_solve_dc_judges_the_residual_and_not_only_the_step():
    """THE FAKE CONVERGENCE, pinned.  `_fet('both')` cascode at
    `(40, 1, 1.2)`: the first iterate puts the supply's branch current at
    7.5e27 A.  The old criterion, ``max|dx| < reltol*max|x_new|`` over
    the WHOLE vector, then tolerated 7e23 V on every node and declared
    convergence in 5 iterations at ``mid = -2.12`` with a KCL residual
    of 7.5e27 -- 688 V (in the branch current, 2 V at the node) from
    the answer.  With `StandardNewton`'s per-component step AND
    residual tests the solve lands where `DC()` does."""
    from pycircuit.circuit.tests.test_limit_fet import _fet, _cascode
    c = _cascode(_fet('both'), 40.0, 1.0, 1.2)
    ref = np.asarray(DC(_cascode(_fet('none'), 40.0, 1.0, 1.2),
                        toolkit=numeric).solve().x, float)
    x, _v, its = P.solve_dc(c, gnd, reltol=1e-4, abstol=1e-6, maxiter=100)
    im = c.get_node_index('mid')
    assert abs(x[im] - ref[im]) < 1e-3, (x[im], ref[im])
    F = np.asarray(c.i(x) + c.u(0, analysis='dc'), float)
    F[c.get_node_index(gnd)] = 0.0
    assert np.max(np.abs(F)) < 1e-6, np.max(np.abs(F))
    assert its > 5


def test_the_redundant_law_is_load_bearing_on_level_1():
    """THE CYCLE DECISION'S MEASUREMENT.  The alternative -- carry the
    tree and DROP the redundant `(b,d)` law -- was priced: with it a
    no-op, Level 1's cascode from a uniform 20 V start fails at all 48
    grid points (`vbd = vbs - vds` is then unbounded forward and the
    bulk-drain diode overflows); with the law applied over the tree 38
    converge.  `_mos4` from the same start: 1676 Jacobians and 4
    failures against 2838 and 5.  One point of that, as the pin."""
    real = P.limit_block

    def noop_redundant(probes, coeffs, v_new, v_old, laws):
        m = len(v_new)
        return real(probes, coeffs, v_new, v_old,
                    list(laws[:m]) + [lambda a, b: a] * (len(laws) - m))

    def run():
        c = _l1_cascode()
        x0 = np.full(c.n, 20.0)
        x0[c.get_node_index(gnd)] = 0.0
        try:
            _x, _v, its = P.solve_dc(c, gnd, x0=x0, reltol=1e-4,
                                     abstol=1e-6, maxiter=100)
        except (RuntimeError, np.linalg.LinAlgError):
            return None
        return its

    applied = run()
    P.limit_block = noop_redundant
    try:
        dropped = run()
    finally:
        P.limit_block = real
    assert applied is not None and applied <= 20, applied
    assert dropped is None, dropped


## ---------------------------------------------------------------------------
## The limited seed, and the damper that is off by default (roadmap sec. 15).
## ---------------------------------------------------------------------------

def test_a_sane_branch_voltage_is_seeded_unchanged():
    """The clamp must be INERT on any start a circuit actually reaches.

    This is what makes limiting the seed a free change rather than a
    trade: every branch voltage a Gummel-Poon junction sees in normal
    operation passes through untouched, and only values no junction can
    hold are pulled in.
    """
    from pycircuit.circuit import pcnr
    c = _mirror(card=dict(NPN))
    devs = pcnr.pcnr_devices(c)
    ib = c.get_node_index('nb')
    for vbe in (0.0, 0.3, 0.7, 0.8):
        x = np.zeros(c.n)
        x[ib] = vbe
        raw = np.array([float(x[a] - x[b])
                        for a, b in pcnr.flat_probes(devs)])
        got = pcnr.v_lim_init(devs, x, limit=True)
        assert_allclose(got, raw, rtol=0, atol=0), vbe
    ## ...and a value no junction can hold IS pulled in.
    x = np.zeros(c.n)
    x[ib] = 20.0
    pulled = pcnr.v_lim_init(devs, x, limit=True)
    assert float(np.max(np.abs(pulled))) < 1.0, pulled


def test_the_seed_is_limited_only_where_the_start_is_arbitrary():
    """`solve_dc` clamps its seed; the transient must not.

    The transient seeds from the previous time point's ACCEPTED solution,
    which is a real operating point and the best information available --
    clamping it would discard that.  `limit` defaults to off for exactly
    that reason, and the two transient call sites rely on the default.
    """
    from pycircuit.circuit import pcnr
    c = _mirror(card=dict(NPN))
    devs = pcnr.pcnr_devices(c)
    x = np.zeros(c.n)
    x[c.get_node_index('nb')] = 20.0
    raw = pcnr.v_lim_init(devs, x)                    # default: no clamp
    assert float(np.max(np.abs(raw))) == 20.0, raw
    assert float(np.max(np.abs(
        pcnr.v_lim_init(devs, x, limit=True)))) < 1.0


def test_the_gmin_damper_is_off_by_default_and_has_lost_its_case():
    """`solve_dc(gmin=...)` regularises the Schur complement.

    It was built in sec. 31 for one circuit: the EKV differential pair
    from a wild start, which raised `LinAlgError` undamped and converged
    with it.  It shipped OFF anyway, because it traded that circuit for
    cascode grid point (2, 2, 2) at every value from 1e-12 to 1e-6.

    ⚠ **That justification is gone.** `limit_delta` on EKV's bulk (sec.
    34) makes the pair converge undamped, and a 72-configuration probe
    afterwards found nothing that `solve_dc` fails on. So the damper now
    rescues **no circuit anyone can name**, while still costing the grid
    point when switched on.

    It is kept rather than removed -- "nothing I can find fails" is not
    "nothing fails", and it is the only tool for a genuinely rank-
    deficient start -- but its default must stay off, and what this test
    can still assert is the property that made it safe to ship at all:
    **the anchor is on the JACOBIAN only, never the residual, so it
    cannot move the answer.**
    """
    import inspect
    from pycircuit.circuit import pcnr
    from pycircuit.circuit.tests.test_limit_identity import _diffpair
    assert inspect.signature(pcnr.solve_dc).parameters['gmin'].default == 0.0

    c = _diffpair(0.0, False)
    x0 = np.full(c.n, 20.0)
    x0[c.get_node_index(gnd)] = 0.0

    ## undamped: converges now (it did not when the damper was built)
    x_off, _, its_off = pcnr.solve_dc(_diffpair(0.0, False), gnd, x0=x0)
    assert its_off < 20

    ## damped: same answer, because gmin never touches the residual
    x_on, _, its_on = pcnr.solve_dc(c, gnd, x0=x0, gmin=1e-9)
    assert its_on < 20
    i_tail = c.get_node_index('tail')
    assert abs(float(x_on[i_tail]) - float(x_off[i_tail])) < 1e-6

    ## ...and against the ordinary path too
    ref = DC(_diffpair(0.0, False)).solve(x0=x0)
    assert abs(float(x_on[i_tail]) - float(ref.v('tail', gnd))) < 1e-6


## ------------------------------------------------------------------------
## Vector PCNR Stage 3 (roadmap sec. 49.2): the traced limiter twin.


def _vec_dev(name, nargs, toolkit):
    import pycircuit.circuit.circuit as _cm
    from pycircuit.circuit import elements_hdl as _eh
    _cm.default_toolkit = toolkit
    cls = getattr(_eh, name)
    el = cls(*([None] * nargs), toolkit=toolkit)
    params = {q.name: getattr(el.iparv, q.name) for q in el.instparams}
    return cls, params


def test_pcnr_limit_branchless_is_the_same_law_as_the_cpu_form():
    """The traced twin must be the device's OWN law, not a second one.

    PCNR's modularity claim is that the device supplies the limiter and the
    architecture supplies only the one-unknown-per-quantity structure. A twin
    that drifted from `pcnr_limit` would make the JAX backend's operating
    point a different device's answer.
    """
    from pycircuit.circuit.circuit import defaultepar
    for name, nargs in (('EkvNmosHdl', 4), ('GummelPoonNpnHdl', 3)):
        cls, params = _vec_dev(name, nargs, numeric)
        m = len(cls.pcnr_probes)
        v_new = np.linspace(0.6, 1.1, m)
        v_old = np.linspace(0.5, 0.9, m)
        x_old = np.zeros(8)
        a = np.asarray(cls.pcnr_limit(v_new, v_old, params, defaultepar,
                                      numeric, x_old), dtype=float)
        b = np.asarray(cls.pcnr_limit_branchless(v_new, v_old, params,
                                                 defaultepar, numeric,
                                                 x_old), dtype=float)
        assert np.array_equal(a, b), (name, a, b)


def test_a_redundant_probe_device_is_refused_by_name():
    """`MosLevel1Hdl` has a redundant probe; the traced twin says so.

    `pcnr.limit_block`'s rule for an over-determined set is a data-dependent
    sort -- decreasing correction, drop what closes a cycle -- and that is the
    one thing in this path that does not become a select. Refusing by name is
    the alternative to silently applying a DIFFERENT law, which is the failure
    mode sec. 47 spent its length on.
    """
    from pycircuit.circuit.circuit import defaultepar
    cls, params = _vec_dev('MosLevel1Hdl', 4, numeric)
    m = len(cls.pcnr_probes)
    with pytest.raises(NotImplementedError, match='REDUNDANT probe'):
        cls.pcnr_limit_branchless(np.linspace(0.6, 1.1, m),
                                  np.linspace(0.5, 0.9, m), params,
                                  defaultepar, numeric, np.zeros(8))


def test_the_traced_twin_traces_and_differentiates():
    """`jit` for the Newton loop, a finite Jacobian for PCNR's own solve."""
    pytest.importorskip('jax')
    import jax
    import jax.numpy as jnp
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.circuit import defaultepar
    cls, params = _vec_dev('GummelPoonNpnHdl', 3, jaxtoolkit)
    m = len(cls.pcnr_probes)
    v_new = jnp.linspace(0.6, 1.1, m)
    v_old = jnp.linspace(0.5, 0.9, m)

    def f(a, b):
        return cls.pcnr_limit_branchless(a, b, params, defaultepar,
                                         jaxtoolkit, jnp.zeros(8))

    eager = np.asarray(f(v_new, v_old), dtype=float)
    jitted = np.asarray(jax.jit(f)(v_new, v_old), dtype=float)
    assert np.array_equal(eager, jitted), (eager, jitted)
    jac = np.asarray(jax.jacobian(lambda a: f(a, v_old))(v_new))
    assert np.all(np.isfinite(jac)), jac


def test_a_solution_reading_limiter_parameter_traces_from_cache_too():
    """The positive form of a test that expired twice in one day.

    `EkvNmosHdl`'s `fet` law takes a parameter that reads the solution
    (SPICE's `von`, `_wants_x`), so under tracing it is handed a TRACER --
    and its chain was compiled against the numpy kernel. `_limit_par_for`
    now recompiles it from sympy ingredients that `_limit_par_fn` attaches
    and, since `CACHE_FORMAT` 4, the frozen function record carries through
    the compile cache.

    ⚠ `x_old_sub` is built INSIDE the jitted function on purpose. Closing
    over a `jnp.zeros` created outside makes it a compile-time constant that
    a numpy chain can consume, and a probe written that way passed while the
    feature did not work. The tracer has to be genuine for this to mean
    anything.
    """
    pytest.importorskip('jax')
    import jax
    import jax.numpy as jnp
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.circuit import defaultepar

    cls, params = _vec_dev('EkvNmosHdl', 4, jaxtoolkit)
    laws = cls.pcnr_limit_branchless.__defaults__[0]
    assert any(getattr(pf, '_wants_x', False) for _k, pfs in laws
               for pf in pfs), 'EKV should read the solution somewhere'

    def f(a, b):
        ## Traced, not closed over -- see the docstring.
        return cls.pcnr_limit_branchless(a, b, params, defaultepar,
                                         jaxtoolkit, jnp.zeros(8))

    m = len(cls.pcnr_probes)
    v_new, v_old = jnp.linspace(0.6, 1.1, m), jnp.linspace(0.5, 0.9, m)
    eager = np.asarray(f(v_new, v_old), dtype=float)
    jitted = np.asarray(jax.jit(f)(v_new, v_old), dtype=float)
    assert np.array_equal(eager, jitted), (eager, jitted)
    jac = np.asarray(jax.jacobian(lambda a: f(a, v_old))(v_new))
    assert np.all(np.isfinite(jac)), jac

    ## Still the CPU's own answer, bit for bit.
    cls_n, params_n = _vec_dev('EkvNmosHdl', 4, numeric)
    cpu = np.asarray(cls_n.pcnr_limit(np.asarray(v_new), np.asarray(v_old),
                                      params_n, defaultepar, numeric,
                                      np.zeros(8)), dtype=float)
    assert np.array_equal(eager, cpu), (eager, cpu)


def test_the_cache_carries_what_a_traced_recompile_needs():
    """A cache hit must not silently cost a capability.

    The failure this guards is specific and was measured: with the
    ingredients absent from the frozen record, the same class reported them
    present on a first run into an empty cache dir and absent on the second,
    so the device kept working and quietly could not be traced. Bumping the
    format alone did NOT fix it -- that invalidates stale entries while the
    newly written ones drop the ingredients just the same, which is the
    difference between invalidating a cache and fixing what it stores.
    """
    from pycircuit.circuit import _hdl_cache
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.hdl import _limit_par_for

    cls, _params = _vec_dev('EkvNmosHdl', 4, jaxtoolkit)
    laws = cls.pcnr_limit_branchless.__defaults__[0]
    wants_x = [pf for _k, pfs in laws for pf in pfs
               if getattr(pf, '_wants_x', False)]
    assert wants_x

    for pf in wants_x:
        assert hasattr(pf, '_hdl_limit_par'), (
            'a solution-reading limiter parameter reached the class without '
            'its recompile ingredients -- a cache hit that costs tracing')
        assert _limit_par_for(pf, jaxtoolkit) is not pf

    ## And the record format is the one that carries them.
    assert _hdl_cache.CACHE_FORMAT >= 4


def test_the_freeze_thaw_round_trip_keeps_the_recompile_ingredients():
    """The round trip DIRECTLY, so this bites on a cold cache too.

    ⚠ The test above only fails when the cache is warm -- on a first run the
    class is built rather than restored, so the ingredients are present
    however the record is written. That is a hole: CI with a fresh cache
    would pass over exactly the regression this work fixed. Freezing and
    thawing here removes the dependence on cache state, which is the only
    version of this assertion that holds in both conditions.
    """
    from pycircuit.circuit import _hdl_cache
    from pycircuit.circuit.toolkit import jaxtoolkit

    cls, _params = _vec_dev('EkvNmosHdl', 4, jaxtoolkit)
    laws = cls.pcnr_limit_branchless.__defaults__[0]
    pf = next(pf for _k, pfs in laws for pf in pfs
              if getattr(pf, '_wants_x', False))
    if not hasattr(pf, '_hdl_limit_par'):
        pytest.skip('class was restored from cache without ingredients; the '
                    'test above is the one that catches that')

    thawed = _hdl_cache.thaw(_hdl_cache.freeze({'f': pf}))['f']
    assert getattr(thawed, '_wants_x', None) is True
    assert hasattr(thawed, '_hdl_limit_par'), (
        'freeze/thaw dropped the recompile ingredients: a function restored '
        'from cache would work and silently not be traceable')
    assert thawed._hdl_limit_par is not None


def test_gummel_poon_traces_regardless_of_the_cache():
    """The device whose limiter parameters are CONSTANTS needs no twin.

    Kept beside the EKV test deliberately: if a change ever breaks only one
    of the two, which one says whether the fault is in the parameter chain or
    in the law.
    """
    pytest.importorskip('jax')
    import jax
    import jax.numpy as jnp
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.circuit import defaultepar

    cls, params = _vec_dev('GummelPoonNpnHdl', 3, jaxtoolkit)
    m = len(cls.pcnr_probes)
    v_new, v_old = jnp.linspace(0.6, 1.1, m), jnp.linspace(0.5, 0.9, m)

    def f(a, b):
        return cls.pcnr_limit_branchless(a, b, params, defaultepar,
                                         jaxtoolkit, jnp.zeros(8))

    eager = np.asarray(f(v_new, v_old), dtype=float)
    assert np.array_equal(eager, np.asarray(jax.jit(f)(v_new, v_old),
                                            dtype=float))
    cls_n, params_n = _vec_dev('GummelPoonNpnHdl', 3, numeric)
    cpu = np.asarray(cls_n.pcnr_limit(np.asarray(v_new), np.asarray(v_old),
                                      params_n, defaultepar, numeric,
                                      np.zeros(8)), dtype=float)
    assert np.array_equal(eager, cpu), (eager, cpu)


## ------------------------------------------------------------------------
## Stage 3 (roadmap sec. 49): the traced DEVICE-shaped path.


def _jax_bits():
    pytest.importorskip('jax')
    import jax.numpy as jnp
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import (_device_arrays,
                                                pcnr_vector_blocks,
                                                pcnr_vector_inner_loop)
    return jnp, jaxtoolkit, _device_arrays, pcnr_vector_blocks, \
        pcnr_vector_inner_loop


def _mos_probe_circuit(toolkit):
    import pycircuit.circuit.circuit as _cm
    from pycircuit.circuit import elements_hdl as _eh
    _cm.default_toolkit = toolkit
    c = SubCircuit(toolkit=toolkit)
    nd, ng = c.add_node('d'), c.add_node('g')
    c['vg'] = VS(ng, gnd, v=1.2, toolkit=toolkit)
    c['vd'] = VS(nd, gnd, v=1.8, toolkit=toolkit)
    c['Rd'] = R(nd, gnd, r=1e4, toolkit=toolkit)
    c['M'] = _eh.EkvNmosHdl(nd, ng, gnd, gnd, toolkit=toolkit)
    return c


def test_the_traced_blocks_equal_the_cpus_augmented_system():
    """The assembly, before any Newton loop is trusted.

    If these agree, the rest of Stage 3 is the loop rather than the maths --
    which is why this is asserted separately from the solve below, at a
    RANDOM (x, v_lim) rather than at a solution, where a wrong block could
    still look right.
    """
    jnp, jaxtoolkit, _device_arrays, blocks, _loop = _jax_bits()
    from pycircuit.circuit.circuit import defaultepar
    from pycircuit.circuit.pcnr import pcnr_devices, augmented_system

    c = _mos_probe_circuit(numeric)
    devs = pcnr_devices(c)
    n, k = c.n, sum(d.m for d in devs)
    rng = np.random.default_rng(0)
    x = rng.uniform(-0.3, 1.2, n)
    v = rng.uniform(0.2, 0.9, k)
    ref = augmented_system(c, x, v, devs, defaultepar)[:5]

    cj = _mos_probe_circuit(jaxtoolkit)
    got = blocks(cj, _device_arrays(cj, defaultepar), jnp.asarray(x),
                 jnp.asarray(v), defaultepar, cj.n)
    for name, a, b in zip(('g_mna', 'g_lim', 'J_mm', 'J_ml', 'J_lm'),
                          ref, got):
        a, b = np.asarray(a, float), np.asarray(b, float)
        assert a.shape == b.shape, (name, a.shape, b.shape)
        assert np.array_equal(a, b), (name, np.max(np.abs(a - b)))


def test_vector_pcnr_finds_the_cpus_operating_point_on_the_diff_pair():
    """G1 and G2 together, on the circuit Stage 2 was built for.

    The EKV differential pair is what took Stage 2 from `[14, FAIL]` to
    `[22, 22]` on the CPU -- two devices, three probes each, which is the
    inter-device clash vector PCNR was commissioned to remove.

    ⚠ Gated on AGREEMENT, not on the solve succeeding. Sec. 40 records a
    PCNR answer that converged happily with an internal node out by a whole
    gate voltage and the terminal current still matching to seven digits;
    convergence saw nothing.
    """
    jnp, jaxtoolkit, _device_arrays, _blocks, loop = _jax_bits()
    import pycircuit.circuit.circuit as _cm
    from pycircuit.circuit.circuit import defaultepar
    from pycircuit.circuit import pcnr as _pcnr
    from pycircuit.circuit.tests.test_limit_identity import _diffpair

    _cm.default_toolkit = numeric
    x_ref = np.asarray(_pcnr.solve_dc(_diffpair(0.0, False), gnd,
                                      x0=np.zeros(_diffpair(0.0, False).n))[0],
                       dtype=float)

    _cm.default_toolkit = jaxtoolkit
    cj = _diffpair(0.0, False)
    devs = _device_arrays(cj, defaultepar)
    assert len(devs) == 2 and [d.m for d in devs] == [3, 3], \
        [(d.instance, d.m) for d in devs]
    u = jnp.asarray(np.asarray(cj.u(0.0, defaultepar, analysis='dc'),
                               dtype=float))
    st = loop(cj, devs, jnp.zeros(cj.n), defaultepar, cj.n,
              cj.get_node_index(gnd), reltol=1e-6, abstol=1e-12,
              xtol=1e-12, maxiter=200, u_extra=u)
    _cm.default_toolkit = numeric

    assert bool(st.converged), 'traced vector PCNR did not converge'
    assert np.max(np.abs(np.asarray(st.x, float) - x_ref)) < 1e-12, \
        np.max(np.abs(np.asarray(st.x, float) - x_ref))


def test_vector_pcnr_converges_from_several_starts():
    """A limiter's job is to make the start not matter."""
    jnp, jaxtoolkit, _device_arrays, _blocks, loop = _jax_bits()
    import pycircuit.circuit.circuit as _cm
    from pycircuit.circuit.circuit import defaultepar
    from pycircuit.circuit.dcanalysis import DC

    c = _mos_probe_circuit(numeric)
    x_ref = np.asarray(DC(c, toolkit=numeric).solve().x, dtype=float)

    cj = _mos_probe_circuit(jaxtoolkit)
    devs = _device_arrays(cj, defaultepar)
    n = cj.n
    iref = cj.get_node_index(gnd)
    u = jnp.asarray(np.asarray(cj.u(0.0, defaultepar, analysis='dc'),
                               dtype=float))
    for level in (0.0, 1.0, 20.0):
        start = np.full(n, level)
        ## The reference row is pinned, not solved: a start that puts a
        ## voltage on it keeps it for ever, which is the solver being right
        ## and the caller being wrong.
        start[iref] = 0.0
        st = loop(cj, devs, jnp.asarray(start), defaultepar, n, iref,
                  reltol=1e-6, abstol=1e-12, xtol=1e-12, maxiter=200,
                  u_extra=u)
        _cm.default_toolkit = numeric
        assert bool(st.converged), level
        assert np.max(np.abs(np.asarray(st.x, float) - x_ref)) < 1e-9, level


def test_a_redundant_probe_device_is_refused_at_setup():
    """Named, before the compile, not somewhere inside a jitted loop."""
    jnp, jaxtoolkit, _device_arrays, _blocks, _loop = _jax_bits()
    import pycircuit.circuit.circuit as _cm
    from pycircuit.circuit import elements_hdl as _eh
    from pycircuit.circuit.circuit import defaultepar

    _cm.default_toolkit = jaxtoolkit
    c = SubCircuit(toolkit=jaxtoolkit)
    nd, ng = c.add_node('d'), c.add_node('g')
    c['vd'] = VS(nd, gnd, v=1.0, toolkit=jaxtoolkit)
    c['M'] = _eh.MosLevel1Hdl(nd, ng, gnd, gnd, toolkit=jaxtoolkit)
    try:
        with pytest.raises(NotImplementedError, match='REDUNDANT probe'):
            _device_arrays(c, defaultepar)
    finally:
        _cm.default_toolkit = numeric


def test_a_scalar_only_circuit_falls_through_to_the_junction_path():
    """G4: nothing that worked before is routed somewhere new.

    `_device_arrays` returns None when no device declares `pcnr_probes`, so a
    circuit of scalar-protocol junctions keeps P19's ported path exactly.
    """
    jnp, jaxtoolkit, _device_arrays, _blocks, _loop = _jax_bits()
    import pycircuit.circuit.circuit as _cm
    from pycircuit.circuit.elements import Diode
    from pycircuit.circuit.circuit import defaultepar

    _cm.default_toolkit = jaxtoolkit
    try:
        c = SubCircuit(toolkit=jaxtoolkit)
        n1 = c.add_node('a')
        c['vs'] = VS(n1, gnd, v=0.7, toolkit=jaxtoolkit)
        c['d'] = Diode(n1, gnd, toolkit=jaxtoolkit)
        assert _device_arrays(c, defaultepar) is None
    finally:
        _cm.default_toolkit = numeric
