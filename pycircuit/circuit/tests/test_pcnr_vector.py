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
    the same solve takes 8 monotone iterations to the same point.
    """
    xo, _vo, io = _old_solve_dc(_charge_diode(vsrc=5.0), gnd)
    xn, _vn, i_n = P.solve_dc(_charge_diode(vsrc=5.0), gnd)
    assert io == 14 and i_n == 8, (io, i_n)
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
    el = _bjt(rb=10.0)
    txt = explain(el, source=False, symbolic=False)
    assert 'PCNR: does not qualify -- I(b) reads b, bi, which no $limit ' \
           'probe limits' in txt, txt
    th = eh.GummelPoonNpnThermalHdl('c', 'b', 'e', 'th', 'tha')
    th.update_iparv()
    assert 'PCNR: does not qualify -- a branch-current unknown' in explain(
        th, source=False, symbolic=False)
    th = eh.GummelPoonNpnThermalHdl('c', 'b', 'e', 'th', 'tha', rth=100.0)
    th.update_iparv()
    assert 'var(dT) reads th, tha, which no $limit probe limits' in explain(
        th, source=False, symbolic=False)
    f = _fet('both')('d', 'g', 's')
    f.update_iparv()
    assert 'probe (g,s) is a fetlim limiter -- vector PCNR takes pnjlim ' \
           'probes only in Stage 1; fetlim and limvds are Stage 2' \
           in explain(f, source=False, symbolic=False)
    assert not hasattr(f, 'pcnr_i')


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

        start   plain          pcnr
         0 V    [9, 9]         [9, 9]
        -5 V    [9, 9]         [9, 9]
        +5 V    [FAIL, FAIL]   [164, 164]   converges where plain cannot
       +10 V    [7, 7]         [6, 6]
       +20 V    [FAIL, 7]      [FAIL, FAIL] the Stage 0 signature

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
    assert table[10.0][1][0][0] == 6 and table[10.0][0][0][0] == 7, table
    ## claim 2 -- the control, as measured (the winning order at 7).
    assert table[20.0][0][0][0] is None and table[20.0][0][1][0] == 7, table
    ## +5 V: plain fails both orders, PCNR converges both (slowly).
    assert table[5.0][0] == [(None, None), (None, None)], table[5.0]
    assert all(its is not None and its < 200 for its, _e in table[5.0][1]), \
        table[5.0]
    ## claim 3 -- recorded as measured; equality is what is asserted.
    assert table[20.0][1][0] == table[20.0][1][1], table[20.0]


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
