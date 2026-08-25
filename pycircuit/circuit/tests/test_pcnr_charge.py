# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""PCNR with a charge-storing participant (hdl.md sec. 9, phase A1).

A junction device with capacitance is the normal case, not an exotic
one, and until now `pcnr.augmented_system` refused it outright.  It is
allowed because the resulting Newton system is exactly consistent:

* `g_MNA`'s charge part is a function of ``x_MNA`` and its derivative
  ``Geq`` stays in ``J_MNA/MNA``;
* the device's resistive current is a function of ``v_lim`` alone and
  its derivative goes to ``J_MNA/lim``;
* ``J_lim/lim`` is the identity, ``J_lim/MNA`` the incidence row.

This module proves that by finite-differencing the whole augmented
residual, and then checks the thing a user cares about: PCNR on and off
produce the same transient.
"""

import warnings

import numpy as np
import pytest
import sympy
from numpy.testing import assert_allclose

import pycircuit.circuit.circuit
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.elements import SubCircuit, VS, VSin, R
from pycircuit.circuit import pcnr as _pcnr
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution, ddt,
                                   limexp, vt)
from pycircuit.utilities.param import Parameter


class DiodeWithCj(Behavioural):
    """A diode that also stores charge -- the case A1 unblocks."""
    instparams = [Parameter(name='IS', desc='Sat current', unit='A',
                            default=1e-13),
                  Parameter(name='cj', desc='Junction capacitance', unit='F',
                            default=1e-12)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I,
                            IS * (sympy.exp(b.V / vt()) - 1)   # noqa: F821
                            + ddt(cj * b.V))                   # noqa: F821


def _circuit(vsrc=1.0, sinusoid=False):
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    na, nb = c.add_node('a'), c.add_node('b')
    if sinusoid:
        c['vs'] = VSin(na, gnd, va=vsrc, freq=1e6)
    else:
        c['vs'] = VS(na, gnd, v=vsrc)
    c['R'] = R(na, nb, r=1e3)
    c['D'] = DiodeWithCj(nb, gnd, IS=1e-13, cj=1e-11)
    c.update_iparv()
    return c


def test_charge_storing_element_participates():
    """It qualifies now -- and it really does store charge, so this is
    the case the layer used to refuse."""
    assert DiodeWithCj.pcnr_junctions == ((0, 1),)
    c = _circuit()
    js = _pcnr.pcnr_junctions(c)
    assert len(js) == 1 and js[0][0] == 'D'
    Cm = np.asarray(c['D'].C(np.zeros(2)), float)
    assert np.any(np.abs(Cm) > 0.0), Cm


def _augmented(c, x, v_lim, h, q_last):
    """Assemble (g, J) of the augmented system the way the transient
    does: backward Euler companion terms recomputed at THIS x."""
    from pycircuit.circuit.circuit import defaultepar
    epar = defaultepar
    q = np.asarray(c.q(x, epar), float)
    Cm = np.asarray(c.C(x, epar), float)
    iq = (q - q_last) / h
    Geq = Cm / h
    u = np.asarray(c.u(0.0, epar, analysis='tran'), float)
    js = _pcnr.pcnr_junctions(c)
    g_mna, g_lim, J_mm, J_ml, J_lm, _didv = _pcnr.augmented_system(
        c, x, v_lim, js, epar, u_extra=iq + u, dense_blocks=True,
        J_extra=Geq)
    n, k = len(x), len(js)
    g = np.concatenate([np.asarray(g_mna, float), np.asarray(g_lim, float)])
    J = np.zeros((n + k, n + k))
    J[:n, :n] = J_mm
    J[:n, n:] = J_ml
    J[n:, :n] = J_lm
    J[n:, n:] = np.eye(k)
    return g, J


def test_augmented_jacobian_is_the_derivative_with_charge():
    """The claim in pcnr.augmented_system's comment, finite-differenced:
    with a charge-storing participant, J == dg/dx over BOTH blocks --
    the MNA unknowns and the limited unknown."""
    c = _circuit()
    n = c.n
    x = np.array([1.0, 0.6, -4e-4][:n]) if n == 3 else np.linspace(0.2, 0.6, n)
    v_lim = np.array([0.55])
    h = 1e-9
    q_last = np.zeros(n)

    g0, J = _augmented(c, x, v_lim, h, q_last)
    k = len(v_lim)
    eps = 1e-7
    Jfd = np.zeros_like(J)
    for j in range(n):                      # d/dx_MNA
        dx = np.zeros(n); dx[j] = eps
        gp, _ = _augmented(c, x + dx, v_lim, h, q_last)
        gm, _ = _augmented(c, x - dx, v_lim, h, q_last)
        Jfd[:, j] = (gp - gm) / (2 * eps)
    for j in range(k):                      # d/dv_lim
        dv = np.zeros(k); dv[j] = eps
        gp, _ = _augmented(c, x, v_lim + dv, h, q_last)
        gm, _ = _augmented(c, x, v_lim - dv, h, q_last)
        Jfd[:, n + j] = (gp - gm) / (2 * eps)

    scale = max(1.0, float(np.max(np.abs(J))))
    assert_allclose(J, Jfd, rtol=2e-4, atol=1e-6 * scale)


def test_charge_participant_dc_matches_non_pcnr():
    c = _circuit(vsrc=5.0)
    v_off = float(DC(c, toolkit=numeric, pcnr=False).solve().v('b'))
    v_on = float(DC(c, toolkit=numeric, pcnr=True).solve().v('b'))
    ## PCNR's own convergence tolerance, not a modelling difference: the
    ## two solvers stop at slightly different points on the same curve
    ## (the hand-written Diode shows the same ~2e-6 spread).
    assert_allclose(v_on, v_off, rtol=1e-5)
    assert 0.4 < v_off < 0.8, v_off        # a real forward-biased junction


def test_charge_participant_transient_matches_non_pcnr():
    """The user-visible claim: with capacitance present, a PCNR transient
    reproduces the ordinary one.  This raised NotImplementedError before
    phase A1."""
    def run(pcnr):
        c = _circuit(vsrc=2.0, sinusoid=True)
        tran = Transient(c, toolkit=numeric, pcnr=pcnr)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(tend=2e-6, timestep=2e-8, fixed_timestep=True)
        return np.asarray(res.v('b').y, float)

    y_off = run(False)
    y_on = run(True)
    assert y_on.shape == y_off.shape
    assert_allclose(y_on, y_off, rtol=1e-5, atol=1e-9)
    assert np.ptp(y_off) > 0.1             # the waveform actually swings


def test_charge_participant_transient_used_to_be_refused():
    """Guard the refusal's removal: the layer must no longer raise for a
    charge-storing participant (and the reason it may not is that the
    Jacobian test above passes)."""
    c = _circuit()
    from pycircuit.circuit.circuit import defaultepar
    js = _pcnr.pcnr_junctions(c)
    n = c.n
    x = np.zeros(n)
    Cm = np.asarray(c.C(x, defaultepar), float)
    _pcnr.augmented_system(c, x, np.zeros(len(js)), js, defaultepar,
                           u_extra=np.zeros(n), dense_blocks=True,
                           J_extra=Cm / 1e-9)      # J_extra set == transient


def test_a_pcnr_participant_is_never_evaluated_at_the_node_voltages():
    """A correction to a correction (2026-08-25, vector PCNR Stage 1).

    This test used to be ``test_limexp_is_what_makes_a_pcnr_participant_
    robust`` and asserted that a raw-``exp`` diode DIES under PCNR while
    the ``limexp`` one survives, because ``pcnr.augmented_system``
    assembled ``cir.i(x)`` -- participant included, at the NODE voltage
    -- and subtracted the device's own ``i(sub)`` again: ``inf - inf =
    nan``.  That was true, and it was the visible edge of a defect that
    did not need ``inf`` to bite: with ``expl`` keeping the current
    finite, a BJT mirror's iterate at vbe = 5.2 V on the nodes left
    cancellation NOISE of 1e56 in the base row, and which noise depended
    on the assembly's summation order -- 11 iterations in one instance
    order, 179 in the other, from identical states.

    The layer now EXCLUDES a participant from the assembly (the paper's
    loading: the device is evaluated at its own ``v`` and nowhere else),
    so there is nothing to cancel and the node voltage at which the
    device would have overflowed is never visited by its own ``i``.
    Both diodes therefore converge under PCNR, to the same answer the
    ordinary path finds.  ``limexp`` still earns its place on the
    ``pcnr=False`` path and for the CHARGE, which stays in the MNA block
    at the node voltage.
    """
    from pycircuit.circuit.hdl import limexp

    def diode(expfn):
        class _D(Behavioural):
            instparams = [Parameter(name='IS', desc='Is', unit='A',
                                    default=1e-13)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I,
                                    IS * (expfn(b.V / vt()) - 1))  # noqa
        return _D

    pycircuit.circuit.circuit.default_toolkit = numeric
    raw, lim = diode(sympy.exp)('p', 'n'), diode(limexp)('p', 'n')
    raw.update_iparv(); lim.update_iparv()
    x20 = np.array([20.0, 0.0])
    assert not np.isfinite(np.asarray(raw.i(x20), float)[0])
    assert np.all(np.isfinite(np.asarray(lim.i(x20), float)))

    def solve(cls, pcnr):
        c = SubCircuit()
        na, nb = c.add_node('a'), c.add_node('b')
        c['vs'] = VS(na, gnd, v=20.0)      # hard forward drive
        c['R'] = R(na, nb, r=1.0)
        c['D'] = cls(nb, gnd, IS=1e-13)
        c.update_iparv()
        return float(DC(c, toolkit=numeric, pcnr=pcnr).solve().v('b'))

    ## Both models solve fine WITHOUT pcnr -- the ordinary path's
    ## continuation ladders handle them.
    v_ref = solve(diode(limexp), False)
    assert 0.7 < v_ref < 1.0, v_ref

    ## With pcnr, BOTH survive now: the participant's own `i` is never
    ## called at the node voltage.  The raw-exp model's `i` at 20 V is
    ## still `inf` (asserted above) -- it is simply not asked.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert_allclose(solve(diode(sympy.exp), True), v_ref, rtol=1e-5)
        assert_allclose(solve(diode(limexp), True), v_ref, rtol=1e-5)

    ## The mechanism, pinned: during `augmented_system` the participant's
    ## `i` and `G` are not consulted at all.
    c = SubCircuit()
    na, nb = c.add_node('a'), c.add_node('b')
    c['vs'] = VS(na, gnd, v=20.0)
    c['R'] = R(na, nb, r=1.0)
    c['D'] = diode(sympy.exp)(nb, gnd, IS=1e-13)
    c.update_iparv()
    calls = []
    real_i = type(c['D']).i

    def spy(self, *a, **k):
        calls.append(1)
        return real_i(self, *a, **k)

    type(c['D']).i = spy
    try:
        js = _pcnr.pcnr_junctions(c)
        x = np.zeros(c.n)
        x[c.get_node_index(na)] = x[c.get_node_index(nb)] = 20.0
        _pcnr.augmented_system(c, x, np.array([0.7]), js, dense_blocks=False)
    finally:
        type(c['D']).i = real_i
    assert calls == [], 'the participant was evaluated at the node voltage'


## ------------------------------------------------------------------------
## Phase A1, multi-junction half: a device that owns SEVERAL limited
## quantities -- which is the case PCNR's clash argument is actually about.


class TwoJunction(Behavioural):
    """Base-emitter and base-collector junctions sharing the base.

    The BJT shape, and the case that could not be expressed at all before
    the junction index reached the device protocol: the layer called
    ``pcnr_i(v, ...)`` with no way to say WHICH junction it was asking
    about.  Both junctions share the base node, so this is precisely the
    situation PCNR exists for -- two limited quantities over overlapping
    terminals, which under classical limiting fight over the same node.
    """
    instparams = [Parameter(name='ISE', desc='Ise', unit='A', default=1e-14),
                  Parameter(name='ISC', desc='Isc', unit='A', default=1e-15)]

    @staticmethod
    def analog(base, emitter, collector):
        be, bc = Branch(base, emitter), Branch(base, collector)
        return (Contribution(be.I, ISE * (limexp(be.V / vt()) - 1)),  # noqa
                Contribution(bc.I, ISC * (limexp(bc.V / vt()) - 1)))  # noqa


def _two_junction_circuit():
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    nb, ne, nc = c.add_node('b'), c.add_node('e'), c.add_node('c')
    c['vb'] = VS(nb, gnd, v=0.75)
    c['Re'] = R(ne, gnd, r=100.0)
    c['Rc'] = R(nc, gnd, r=100.0)
    c['Q'] = TwoJunction(nb, ne, nc)
    c.update_iparv()
    return c


def test_two_junction_device_declares_both():
    """Both junctions are declared, sharing terminal 0 (the base)."""
    assert TwoJunction.pcnr_junctions == ((0, 1), (0, 2))
    c = _two_junction_circuit()
    js = _pcnr.pcnr_junctions(c)
    assert len(js) == 2
    assert all(j[0] == 'Q' for j in js)
    ## They share the base row and differ in the other.
    assert js[0][2] == js[1][2] and js[0][3] != js[1][3]


def test_two_junction_device_answers_per_junction():
    """``jn`` selects which junction: the two saturation currents differ
    by 10x, so the returned currents must too."""
    from pycircuit.circuit.circuit import defaultepar
    params = {'ISE': 1e-14, 'ISC': 1e-15}
    i0 = np.asarray(TwoJunction.pcnr_i(0.6, params, defaultepar, numeric,
                                       jn=0), float)
    i1 = np.asarray(TwoJunction.pcnr_i(0.6, params, defaultepar, numeric,
                                       jn=1), float)
    assert_allclose(i0[0] / i1[0], 10.0, rtol=1e-9)
    assert_allclose(i0[1], -i0[0], rtol=1e-12)
    g0 = np.asarray(TwoJunction.pcnr_didv(0.6, params, defaultepar, numeric,
                                          jn=0), float)
    g1 = np.asarray(TwoJunction.pcnr_didv(0.6, params, defaultepar, numeric,
                                          jn=1), float)
    assert_allclose(g0[0] / g1[0], 10.0, rtol=1e-9)
    ## And the limiter is asked per junction too (different IS -> different
    ## vc, so the limited value differs).
    l0 = float(TwoJunction.pcnr_limit(5.0, 0.7, params, defaultepar,
                                      numeric, jn=0))
    l1 = float(TwoJunction.pcnr_limit(5.0, 0.7, params, defaultepar,
                                      numeric, jn=1))
    assert np.isfinite(l0) and np.isfinite(l1)


def test_two_junction_device_solves_with_and_without_pcnr():
    c = _two_junction_circuit()
    off = DC(c, toolkit=numeric, pcnr=False).solve()
    on = DC(c, toolkit=numeric, pcnr=True).solve()
    for node in ('e', 'c'):
        assert_allclose(float(on.v(node)), float(off.v(node)), rtol=1e-5)
    ## Both junctions are actually conducting, so this exercises the
    ## shared-base case rather than a degenerate one.
    assert float(off.v('e')) > 0.05
    assert float(off.v('c')) > 0.02


def test_jax_pcnr_calls_the_device_instead_of_rebuilding_it():
    """The traced backend now ASKS each device for its junction.

    It used to rebuild every junction itself as ``IS*(exp(v/VT)-1)`` with
    one global VT, reading ``IS`` by name -- so a device whose saturation
    current is ``ISE``/``ISC`` got ``IS = 0``, a junction carrying no
    current, and a confident wrong answer; and multi-junction and
    charge-storing participants were impossible here while the CPU
    accepted them.  Calling the device removes all three limits at once.
    """
    pytest.importorskip('jax')
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import _junction_arrays, JAXTransient
    from pycircuit.circuit.transient import Transient as CpuTransient

    saved = pycircuit.circuit.circuit.default_toolkit
    try:
        pycircuit.circuit.circuit.default_toolkit = jaxtoolkit

        def build(tk):
            pycircuit.circuit.circuit.default_toolkit = tk
            c = SubCircuit(toolkit=tk)
            nb, ne, nc = (c.add_node('b'), c.add_node('e'), c.add_node('c'))
            c['vb'] = VS(nb, gnd, v=0.75, toolkit=tk)
            c['Re'] = R(ne, gnd, r=100.0, toolkit=tk)
            c['Rc'] = R(nc, gnd, r=100.0, toolkit=tk)
            c['Q'] = TwoJunction(nb, ne, nc, toolkit=tk)
            c.update_iparv()
            return c

        cj = build(jaxtoolkit)
        meta = _junction_arrays(cj)
        assert meta is not None
        assert len(meta[0]) == 2            # BOTH junctions, not one
        assert meta[4] is not None          # device-supplied evaluation

        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = JAXTransient(cj, pcnr=True, reltol=1e-5).solve(
                gnd, tend=1e-6, timestep=2e-8, uic=True,
                fixed_timestep=True)
            y_jax = float(np.asarray(res.v('e'), float).reshape(-1)[-1])

            cc = build(numeric)
            y_cpu = float(np.asarray(CpuTransient(
                cc, toolkit=numeric, pcnr=True).solve(
                    tend=1e-6, timestep=2e-8,
                    fixed_timestep=True).v('e').y, float)[-1])
        ## The device is evaluated by the same code on both backends now,
        ## so they must agree closely rather than merely both converge.
        assert_allclose(y_jax, y_cpu, rtol=1e-5)
    finally:
        pycircuit.circuit.circuit.default_toolkit = saved


def test_jax_fallback_charge_check_asks_the_C_matrix():
    """The FALLBACK path -- taken only when a device supplies no
    ``pcnr_i`` for the traced backend to call -- must test whether the
    device ACTUALLY stores charge, not whether it happens to have an
    ``eval_q_pure`` method: every generated element has one, charge or
    not, so the method test refused charge-free devices and blamed a
    charge they did not have."""
    pytest.importorskip('jax')
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import _junction_arrays

    saved = pycircuit.circuit.circuit.default_toolkit
    try:
        pycircuit.circuit.circuit.default_toolkit = jaxtoolkit
        c = SubCircuit(toolkit=jaxtoolkit)
        na, nb = c.add_node('a'), c.add_node('b')
        c['vs'] = VS(na, gnd, v=1.0, toolkit=jaxtoolkit)
        c['R'] = R(na, nb, r=1e3, toolkit=jaxtoolkit)
        ## A generated, charge-FREE diode: has eval_q_pure, stores nothing.
        from pycircuit.circuit.hdl import Behavioural as _B
        class _D(_B):
            instparams = [Parameter(name='IS', desc='Is', unit='A',
                                    default=1e-13)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I,
                                    IS * (sympy.exp(b.V / vt()) - 1))  # noqa
        el = _D(nb, gnd, IS=1e-13, toolkit=jaxtoolkit)
        c['D'] = el
        c.update_iparv()
        assert hasattr(el, 'eval_q_pure')          # the misleading signal
        assert not np.any(np.abs(np.asarray(el.C(np.zeros(el.n)),
                                            float)) > 0)
        assert _junction_arrays(c) is not None     # accepted, correctly
    finally:
        pycircuit.circuit.circuit.default_toolkit = saved


## ------------------------------------------------------------------------
## Vector PCNR, Stage 1 (roadmap sec. 15): the same finite-difference
## claim on a participant with TWO limited unknowns and charge.


def _bjt_circuit():
    """A common-emitter stage on the full NPN card at ``rb = re = rc = 0``
    -- every charge block on (cje, cjc, tf/xtf/vtf, tr, xcjc), two
    limited junctions sharing the base."""
    from pycircuit.circuit.tests.test_elements_hdl_library3 import NPN
    from pycircuit.circuit import elements_hdl as eh
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    nb, nc, nvi, nvcc = (c.add_node('nb'), c.add_node('nc'),
                         c.add_node('nvi'), c.add_node('nvcc'))
    c['vcc'] = VS(nvcc, gnd, v=5.0)
    c['vin'] = VS(nvi, gnd, v=1.0)
    c['rb'] = R(nvi, nb, r=47e3)
    c['rc'] = R(nvcc, nc, r=4.7e3)
    c['q1'] = eh.GummelPoonNpnHdl(nc, nb, gnd, **NPN)
    c.update_iparv()
    return c


def _bjt_fd_point(c):
    x = np.zeros(c.n)
    x[c.get_node_index('nvcc')] = 5.0
    x[c.get_node_index('nvi')] = 1.0
    x[c.get_node_index('nb')] = 0.76
    x[c.get_node_index('nc')] = 2.4
    ## Deliberately OFF the branch voltages, so g_lim != 0 and the
    ## cross-coupling of v into g_mna is exercised, not just the diagonal.
    v_lim = np.array([0.76 + 0.01, 0.76 - 2.4 - 0.02])
    return x, v_lim


def test_augmented_jacobian_is_the_derivative_on_a_vector_participant_with_charge():
    """`J_eff == df_eff/dx` and the full augmented ``J == dg/d[x; v]``,
    finite-differenced, on the BJT: the ``3 x 2`` block in ``J_MNA/lim``,
    the incidence rows in ``J_lim/MNA``, the identity in ``J_lim/lim``,
    and the companion ``Geq`` from every charge block in ``J_MNA/MNA``."""
    c = _bjt_circuit()
    assert c['q1'].pcnr_probes == (((1, 2), 'pnj'), ((1, 0), 'pnj'))
    assert np.any(np.abs(np.asarray(c['q1'].C(np.array([2.4, 0.76, 0.0])),
                                    float)) > 0.0)
    n = c.n
    from pycircuit.circuit.circuit import defaultepar
    x, v_lim = _bjt_fd_point(c)
    h = 1e-9
    q_last = np.asarray(c.q(x * 0.9, defaultepar), float)

    g0, J = _augmented(c, x, v_lim, h, q_last)
    k = len(v_lim)
    Jfd = np.zeros_like(J)
    for j in range(n):
        eps = 1e-7 * max(1.0, abs(x[j]))
        dx = np.zeros(n)
        dx[j] = eps
        gp, _ = _augmented(c, x + dx, v_lim, h, q_last)
        gm, _ = _augmented(c, x - dx, v_lim, h, q_last)
        Jfd[:, j] = (gp - gm) / (2 * eps)
    for j in range(k):
        eps = 1e-7
        dv = np.zeros(k)
        dv[j] = eps
        gp, _ = _augmented(c, x, v_lim + dv, h, q_last)
        gm, _ = _augmented(c, x, v_lim - dv, h, q_last)
        Jfd[:, n + j] = (gp - gm) / (2 * eps)

    scale = max(1.0, float(np.max(np.abs(J))))
    assert_allclose(J, Jfd, rtol=2e-4, atol=1e-6 * scale)
    ## The block is where the device's conductance lives: the MNA/MNA
    ## block carries NO resistive conductance of the BJT, only Geq.
    assert np.abs(J[n:, n:] - np.eye(k)).max() == 0.0
    nb, nc = c.get_node_index('nb'), c.get_node_index('nc')
    assert np.count_nonzero(J[:n, n:]) == 6           # 3 rows x 2 probes
    Cm = np.asarray(c.C(x, defaultepar), float)
    assert_allclose(J[nb, nc], Cm[nb, nc] / h + 0.0, rtol=1e-9)
    ## and `schur_reduce` is the Schur complement of that J.
    from pycircuit.circuit import pcnr as _p
    js = _p.pcnr_devices(c)
    q = np.asarray(c.q(x, defaultepar), float)
    u = np.asarray(c.u(0.0, defaultepar, analysis='tran'), float)
    g_mna, g_lim, J_mm, _a, _b, didv = _p.augmented_system(
        c, x, v_lim, js, defaultepar,
        u_extra=(q - q_last) / h + u, dense_blocks=False,
        J_extra=Cm / h)
    f_eff, J_eff = _p.schur_reduce(g_mna, g_lim, J_mm, junctions=js, didv=didv)
    J_ml, J_lm = J[:n, n:], J[n:, :n]
    assert_allclose(J_eff, J[:n, :n] - J_ml @ J_lm, rtol=1e-12,
                    atol=1e-12 * scale)
    assert_allclose(f_eff, g0[:n] - J_ml @ g0[n:], rtol=1e-12,
                    atol=1e-12 * max(1.0, np.max(np.abs(g0))))


def test_a_dropped_didv_column_fails_the_finite_difference_check():
    """Mutation control for the test above, on the same point: with the
    BJT's ``d/dvbc`` column zeroed the augmented Jacobian disagrees with
    FD by the transport current's slope, far outside the tolerance."""
    c = _bjt_circuit()
    cls = type(c['q1'])
    real = cls.pcnr_didv

    def dropped(v, params, epar, toolkit):
        blk = np.array(real(v, params, epar, toolkit), dtype=float)
        blk[:, 1] = 0.0
        return blk

    n = c.n
    x, v_lim = _bjt_fd_point(c)
    ## SATURATED, so the base-collector column carries a real slope: at
    ## the forward-active point above `d/dvbc` is the Early term alone
    ## (~1e-5 S) and a dropped column would be caught by only a decade.
    x[c.get_node_index('nc')] = 0.05
    v_lim[1] = 0.76 - 0.05
    h, q_last = 1e-9, np.zeros(n)
    cls.pcnr_didv = staticmethod(dropped)
    try:
        _g, J_bad = _augmented(c, x, v_lim, h, q_last)
    finally:
        cls.pcnr_didv = staticmethod(real)
    _g, J_ok = _augmented(c, x, v_lim, h, q_last)
    eps = 1e-7
    dv = np.zeros(2)
    dv[1] = eps
    gp, _ = _augmented(c, x, v_lim + dv, h, q_last)
    gm, _ = _augmented(c, x, v_lim - dv, h, q_last)
    col = (gp - gm) / (2 * eps)
    scale = max(1.0, float(np.max(np.abs(J_ok))))
    assert_allclose(J_ok[:, n + 1], col, rtol=2e-4, atol=1e-6 * scale)
    assert not np.allclose(J_bad[:, n + 1], col, rtol=2e-4,
                           atol=1e-6 * scale)
    assert np.max(np.abs(J_bad[:, n + 1] - col)) > 1e-4 * scale
