# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Roadmap item 12.3 -- the simulator-level `gmin` anchor.

WHAT THE ITEM IS.  EKV's channel conductance underflows to *exactly*
zero in deep cutoff (``softplus(-800)`` is ``0.0``, not a denormal), so
a stacked pair with no other DC path to the intermediate node hands the
solver a Jacobian with an empty column and gets
``SingularMatrix: 'nm' appears in no equation``.
``compact.PspMosLongChannel`` carries a private ``GLEAK = 1e-12`` for
exactly this; ``elements_hdl.EkvNmosHdl`` deliberately does not, because
1 pS at a volt is 1 pA and that is the size of the very currents the
model exists to measure -- ``test_weak_inversion_...`` below puts a
number on that refusal: **15%** at ``vgs = 0``.

WHAT AN ANCHOR HAS TO GET RIGHT, and what these tests are organised
around: **gmin must rescue a NUMERICALLY empty row and must not
disguise a STRUCTURALLY empty one.**  Those two look identical in the
matrix -- a column of exact zeros -- and differ only in whether the
column is still empty *at the answer*.

So every test here comes in the shape `validation-design`'s
"test the INTERVENTION against its absence" section asks for: the case
that only works with the anchor, the same case proving it fails
without, and -- the assertion that actually finds bugs -- the answer
being *unchanged*.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit.elements import VS, R, C, IS
from pycircuit.circuit.analysis import NoConvergenceError, SingularMatrix
from pycircuit.circuit.dcanalysis import DC
import pycircuit.circuit.elements_hdl as eh


## The card `test_elements_hdl_library3` uses for its limiter work, copied
## rather than imported so this file does not depend on a test module that
## another campaign owns.
LIM_CARD = dict(vto=0.5, gamma=0.7, phi=0.7, kp=1.5e-4, cox=6.9e-3,
                w=20e-6, l=2e-6)


def _stacked_pair(vdd=40.0, vg=0.2, rleak=None):
    """The circuit item 12.3 was written against.

    Two EKV NMOS in series from a 40 V rail with a common gate at 0.2 V.
    Nothing else touches ``nm``, so at the origin -- where every solve
    starts -- both devices are in deep cutoff and the ``nm`` column is
    bit-exactly zero.
    """
    c = SubCircuit()
    for n in ('nd', 'nm', 'ng', 'nvdd'):
        c.add_node(n)
    c['vdd'] = VS(c.get_node('nvdd'), gnd, v=vdd)
    c['vg'] = VS(c.get_node('ng'), gnd, v=vg)
    c['rd'] = R(c.get_node('nvdd'), c.get_node('nd'), r=2e3)
    c['m1'] = eh.EkvNmosHdl(c.get_node('nd'), c.get_node('ng'),
                            c.get_node('nm'), gnd, **LIM_CARD)
    c['m2'] = eh.EkvNmosHdl(c.get_node('nm'), c.get_node('ng'),
                            gnd, gnd, **LIM_CARD)
    if rleak is not None:
        c['rl'] = R(c.get_node('nm'), gnd, r=rleak)
    return c


def _solve_counting_jacobians(cir, **kvargs):
    """Solve, and count Jacobian assemblies -- the quantity the rescue
    costs.  Wall clock is not counted, per `validation-design`."""
    dc = DC(cir, **kvargs)
    n = [0]
    G0 = cir.G

    def G(*a, **kw):
        n[0] += 1
        return G0(*a, **kw)

    cir.G = G
    try:
        res = dc.solve()
    except Exception as exc:                       # returned, not raised
        res = exc
    finally:
        del cir.G
    return n[0], res, dc


# ---------------------------------------------------------------------------
# 1.  The intervention, and its absence
# ---------------------------------------------------------------------------

def test_the_stacked_pair_converges_with_the_anchor_and_fails_without_it():
    """Both halves, on ONE circuit, differing only in ``gmin``.

    Without the second half this test measures nothing: the pair might
    have started converging for some unrelated reason.
    """
    n_on, res_on, dc_on = _solve_counting_jacobians(_stacked_pair())
    assert not isinstance(res_on, Exception), res_on
    assert np.isfinite(float(res_on.v('nm', gnd)))

    n_off, res_off, _ = _solve_counting_jacobians(_stacked_pair(), gmin=0.0)
    assert isinstance(res_off, SingularMatrix), res_off
    assert "'nm'" in str(res_off), str(res_off)

    ## The cost, recorded rather than asserted tightly -- it is a
    ## continuation ladder and its rung count is allowed to move.
    ## Measured at the landing of item 12.3: 24 with, 5 before the
    ## unanchored solve gives up.
    assert n_off < n_on <= 4 * n_off + 40, (n_off, n_on)


def test_gmin_zero_is_the_pre_item_behaviour_and_nothing_else_changed():
    """The off switch is a switch, not a mode.

    ``gmin=0`` must reproduce the exact exception text item 12.3
    recorded, so a caller who wants the old diagnosis keeps it.
    """
    with pytest.raises(SingularMatrix, match="appears in no equation"):
        DC(_stacked_pair(), gmin=0.0).solve()


# ---------------------------------------------------------------------------
# 2.  The answer must not move.  This is the assertion that finds bugs.
# ---------------------------------------------------------------------------

def test_the_anchored_answer_is_the_leak_free_limit():
    """An anchor that moves the answer is a bug, so measure the move.

    The reference is built by *taking the anchor away in the limit*: an
    explicit leak resistor from ``nm`` to ground, swept up until its
    conductance is far below anything in the circuit.  If the anchored
    answer is the leak-free answer, the sequence converges ON it.

    Measured at landing:

        rleak      v(nm)               relative to the anchored answer
        1e9        0.017588037157      2.57e-2
        1e12       0.018051385528      2.66e-5
        1e14       0.018051861024      2.69e-7
        1e16       0.018051865779      5.60e-9

    Three decades of leak, three decades of agreement: the residual is
    the *leak*, not the anchor.  Note what that says about the 1 GOhm
    testbench leak the library tests use -- it is 2.6% away from the
    answer, i.e. **the anchor is the more faithful of the two**.
    """
    anchored = float(DC(_stacked_pair()).solve().v('nm', gnd))

    prev = None
    for rleak, tol in ((1e9, 3e-2), (1e12, 3e-5), (1e14, 3e-7), (1e16, 1e-8)):
        v = float(DC(_stacked_pair(rleak=rleak)).solve().v('nm', gnd))
        assert_allclose(v, anchored, rtol=tol,
                        err_msg='leak %g disagrees by more than %g' % (rleak, tol))
        if prev is not None:
            ## and it is a monotone approach, not a coincidence at one value
            assert abs(v - anchored) < abs(prev - anchored), (rleak, v, prev)
        prev = v

    ## The task's stated comparison, with the number it actually earns:
    ## a 1 GOhm leak is a 2.6% perturbation of this stage.
    v1g = float(DC(_stacked_pair(rleak=1e9)).solve().v('nm', gnd))
    assert 1e-2 < abs(v1g / anchored - 1) < 5e-2, v1g


def test_the_answer_does_not_depend_on_the_value_of_gmin():
    """Four decades of ``gmin``, and the answer moves by 4e-14 V.

    Because the returned point comes from a final ``gmin = 0`` solve --
    the anchor is a path, not a term -- this is not "small", it is
    "absent".  If someone ever makes the anchor permanent by default,
    this is the test that fails first.
    """
    ref = np.asarray(DC(_stacked_pair()).solve().x, dtype=float)
    for g in (1e-13, 1e-11, 1e-10):
        dc = DC(_stacked_pair(), gmin=g)
        x = np.asarray(dc.solve().x, dtype=float)
        assert dc.gmin_anchor_retained is False
        assert float(np.max(abs(x - ref))) < 1e-12, (g, np.max(abs(x - ref)))


def test_a_bias_that_needs_no_anchor_is_bit_for_bit_unchanged():
    """The no-op case, asserted at the strongest tolerance there is.

    The anchor engages only on a singularity, so on a solve that
    succeeds it is not in the matrix at any iterate -- and "unchanged to
    1e-9" would not have distinguished that from a 1 pS shunt that
    happens to be small here.  ``array_equal`` does.
    """
    a = np.asarray(DC(_stacked_pair(vdd=5.0, vg=2.5)).solve().x, dtype=float)
    b = np.asarray(DC(_stacked_pair(vdd=5.0, vg=2.5),
                      gmin=0.0).solve().x, dtype=float)
    assert np.array_equal(a, b)


# ---------------------------------------------------------------------------
# 3.  The measurement the model refused an anchor FOR
# ---------------------------------------------------------------------------

def _subthreshold(vg, vd=1.0, rleak=None):
    c = SubCircuit()
    for n in ('nd', 'ng'):
        c.add_node(n)
    c['vd'] = VS(c.get_node('nd'), gnd, v=vd)
    c['vg'] = VS(c.get_node('ng'), gnd, v=vg)
    c['m1'] = eh.EkvNmosHdl(c.get_node('nd'), c.get_node('ng'), gnd, gnd,
                            **LIM_CARD)
    if rleak is not None:
        c['rl'] = R(c.get_node('nd'), gnd, r=rleak)
    return c


def test_weak_inversion_current_is_unperturbed_where_a_permanent_anchor_would_bite():
    """The reason ``EkvNmosHdl`` refused a private ``GLEAK``, in numbers.

    Two things have to be true for this test to be evidence:

    1. the measurement is SENSITIVE at this bias.  A permanent 1 pS to
       ground -- exactly what a default ``gmin`` term would be -- moves
       the drain current by 15% at ``vgs = 0`` and 4.4% at 0.05 V.  That
       is the counterfactual, and it is asserted, so a future change
       that makes the anchor permanent cannot pass by claiming the
       effect is negligible;
    2. the real answer does not move at all -- ``array_equal``, not a
       tolerance.

    Measured at landing (``vd = 1 V``):

        vgs    Id            with a permanent 1 pS shunt   ratio
        0.00   6.7076e-12    7.7076e-12                    1.149
        0.05   2.2755e-11    2.3755e-11                    1.044
        0.10   7.8572e-11    7.9572e-11                    1.013
    """
    ratios = []
    for vg in (0.0, 0.05, 0.10):
        anchored = DC(_subthreshold(vg)).solve()
        plain = DC(_subthreshold(vg), gmin=0.0).solve()
        assert np.array_equal(np.asarray(anchored.x, dtype=float),
                              np.asarray(plain.x, dtype=float)), vg

        ## the counterfactual: a PERMANENT 1 pS, which is what a naive
        ## always-on gmin would be
        shunted = DC(_subthreshold(vg, rleak=1e12)).solve()
        i_ref = float(anchored.i('vd.plus'))
        i_shunt = float(shunted.i('vd.plus'))
        ratios.append(abs(i_shunt / i_ref))
        ## and the current really is in the regime that makes this hurt
        assert abs(i_ref) < 1e-10, (vg, i_ref)

    assert ratios[0] > 1.10, ratios     # 15% at vgs = 0
    assert ratios[1] > 1.03, ratios
    ## and the sensitivity falls away as the current rises, which is the
    ## shape a fixed 1 pA error has and no other error does
    assert ratios[0] > ratios[1] > ratios[2] > 1.0, ratios


# ---------------------------------------------------------------------------
# 4.  Gate 1 -- an unknown still in no equation AT THE ANSWER
# ---------------------------------------------------------------------------

def test_a_node_reachable_only_through_a_capacitor_still_fails_and_is_named():
    """The classic floating node: nothing but a capacitor touches it, so
    it has no DC path at any bias and never will.

    This is `test_transient_repairs.test_gate_6_1_a_floating_node_is_named`'s
    circuit, kept working by the anchor's gate 1 rather than by the
    anchor never running -- the anchor DOES run here, converges, and the
    answer is then thrown away because the ``floaty`` column is empty in
    the ``gmin = 0`` Jacobian at that answer too.
    """
    c = SubCircuit()
    c['VS'] = VS('n1', gnd, v=1.0)
    c['R1'] = R('n1', 'n2', r=1e3)
    c['Cfl'] = C('floaty', gnd, c=1e-9)

    with pytest.raises(SingularMatrix) as exc:
        DC(c).solve()
    msg = str(exc.value)
    assert 'floaty' in msg, msg
    assert 'REJECTED' in msg, msg
    assert 'manufactured answer' in msg, msg


def test_a_node_with_only_a_current_source_still_fails_and_is_named():
    """The case that shows why terminal incidence is NOT the test.

    ``dead`` has an element attached, so any "is a device connected to
    it?" rule calls it anchorable -- and anchoring it would return
    ``v = -I/gmin = 1e9 V`` with a straight face.  A current source
    contributes to the right-hand side and never to the Jacobian, so the
    column is empty at the answer, and gate 1 catches it there.
    """
    c = SubCircuit()
    c['VS'] = VS('n1', gnd, v=1.0)
    c['R1'] = R('n1', gnd, r=1e3)
    c['I1'] = IS('dead', gnd, i=1e-3)

    with pytest.raises(SingularMatrix) as exc:
        DC(c).solve()
    assert 'dead' in str(exc.value), str(exc.value)


# ---------------------------------------------------------------------------
# 5.  Gate 2 -- the answer must not depend on the anchor
# ---------------------------------------------------------------------------

def _fed_island():
    """An island {a, b} joined by a resistor and fed 1 mA from outside
    with no return path.  Every column is occupied, so gate 1 sees
    nothing wrong; the only place for the milliamp to go is the anchor,
    which puts the island at ``I / (2 * gmin) = 5e8 V``."""
    c = SubCircuit()
    c['VS'] = VS('n1', gnd, v=1.0)
    c['Rg'] = R('n1', gnd, r=1e3)
    c['Rab'] = R('a', 'b', r=1e3)
    c['I1'] = IS(gnd, 'a', i=1e-3)
    return c


def test_an_anchor_that_holds_the_answer_up_is_rejected():
    """Gate 2, on the case gate 1 cannot see."""
    with pytest.raises(SingularMatrix) as exc:
        DC(_fed_island()).solve()
    msg = str(exc.value)
    assert 'DETERMINES the answer' in msg, msg
    assert ("'a'" in msg or "'b'" in msg), msg


def test_the_residual_test_gate_2_replaced_would_have_passed_this_circuit():
    """Why gate 2 is a SENSITIVITY test and not a residual test.

    The obvious gate is "does the anchored answer satisfy the unanchored
    KCL to ``reltol * I_scale + abstol``?" -- the solver's own ``conv_f``
    on the ``gmin = 0`` system.  It reads well and it is fooled here,
    because ``I_scale`` is ``|J|.|x|`` and the anchor inflated ``|x|`` to
    5e8 V.  The tolerance is computed from the corrupted answer, so it
    grows exactly as fast as the corruption.

    This reconstructs the anchor as two explicit 1 TOhm resistors --
    numerically the same matrix -- and measures the miss against the
    tolerance that gate would have used.  It passes by four orders of
    magnitude, which is the whole reason it is not the gate.
    """
    gmin = 1e-12
    anchored = _fed_island()
    anchored['ga'] = R('a', gnd, r=1.0 / gmin)
    anchored['gb'] = R('b', gnd, r=1.0 / gmin)
    dc = DC(anchored)
    x = np.asarray(dc.solve().x, dtype=float)
    assert abs(x[anchored.get_node_index(anchored.get_node('a'))]) > 1e7, x

    pure = _fed_island()                       # the SAME system, gmin = 0
    F = np.asarray(pure.i(x, dc.epar)
                   + pure.u(0, analysis='dc', epar=dc.epar), dtype=float)
    J = abs(np.asarray(pure.G(x, dc.epar), dtype=float))
    ## The gate runs on the REDUCED system -- the reference row is removed
    ## before the solve, and it is the row that carries the imbalance.
    iref = pure.get_node_index(gnd)
    F = np.delete(F, iref)
    J = np.delete(np.delete(J, iref, axis=0), iref, axis=1)
    x = np.delete(x, iref)
    I_scale = J.dot(abs(x)) + abs(F)
    tol = 1e-4 * I_scale + 1e-12
    ## The would-be gate: every row inside tolerance -> "accept".
    assert np.all(abs(F) <= tol), (abs(F) - tol)
    ## and by a margin that makes the point rather than by a whisker
    assert float(np.max(abs(F) / tol)) < 1e-3, float(np.max(abs(F) / tol))


def test_an_undetermined_common_mode_is_anchored_and_says_so():
    """The one shape where the anchor legitimately STAYS.

    An island with no injection has infinitely many solutions -- every
    common mode satisfies KCL exactly -- and no amount of continuation
    picks one.  gmin picks the minimum-norm member, which is what every
    SPICE does, and both gates pass honestly: no column is empty, and
    moving gmin a decade does not move the answer, because the answer
    does not depend on it.

    ``gmin_anchor_retained`` is how a caller finds out, and it is the
    only circuit found in this campaign that sets it.
    """
    c = SubCircuit()
    c['VS'] = VS('n1', gnd, v=1.0)
    c['Rg'] = R('n1', gnd, r=1e3)
    c['Rab'] = R('a', 'b', r=1e3)

    dc = DC(c)
    res = dc.solve()
    assert dc.gmin_anchor_retained is True
    assert_allclose(float(res.v('n1', gnd)), 1.0, rtol=1e-9)
    assert_allclose(float(res.v('a', gnd)), 0.0, atol=1e-9)
    assert_allclose(float(res.v('b', gnd)), 0.0, atol=1e-9)

    ## and it is genuinely a solution: KCL holds on the UNANCHORED system
    x = np.asarray(res.x, dtype=float)
    F = np.asarray(c.i(x, dc.epar) + c.u(0, analysis='dc', epar=dc.epar),
                   dtype=float)
    assert float(np.max(abs(F))) < 1e-12, F


# ---------------------------------------------------------------------------
# 6.  The anchor goes on node rows only
# ---------------------------------------------------------------------------

def test_a_branch_row_is_never_anchored():
    """A conductance in a voltage source's KVL equation is not a leaky
    element, it is a wrong equation -- so the anchor is masked to node
    rows.  ``GminSteppingNewton``'s gshunt ladder does NOT mask, which is
    harmless there only because it always finishes at ``g = 0``.

    Driven directly, because a circuit-level test cannot see this: the
    DC path finishes with a pure solve, so the anchored matrices never
    reach the answer.
    """
    from pycircuit.circuit.nrsolver import (GminAnchorNewton, NonLinearSolver,
                                            StandardNewton)

    ## Rows 0, 1 are node rows; row 2 is a branch row.  Rows 0 and 1 are
    ## linearly dependent, so the pure matrix is singular with no empty
    ## row or column -- the shape gate 1 cannot see.
    J0 = np.array([[0.0, 0.0, 1.0],
                   [0.0, 0.0, -1.0],
                   [1.0, -1.0, 0.0]])
    rhs = np.array([0.0, 0.0, 1.0])

    def eval_FJ(x):
        return J0.dot(x) - rhs, J0.copy()

    class _Recorder(NonLinearSolver):
        def __init__(self, inner):
            self.inner = inner
            self.diags = []

        def solve_system(self, x0, ev, toolkit, *a, **kw):
            self.diags.append(np.diag(np.asarray(ev(np.zeros(3))[1]) - J0))
            return self.inner.solve_system(x0, ev, toolkit, *a, **kw)

    rec = _Recorder(StandardNewton())

    class _AlwaysSingular(NonLinearSolver):
        def solve_system(self, *a, **kw):
            raise SingularMatrix('by construction')

    solver = GminAnchorNewton(_AlwaysSingular(), node_rows=[0, 1],
                              gmin=1e-12, rung_solver=rec)
    tol = np.full(3, 1e-12)
    x, _its = solver.solve_system(np.zeros(3), eval_FJ, numeric, 1e-4,
                                  tol, tol, 100)

    assert solver.anchor_retained is True
    ## x0 - x1 = 1 and the anchor splits it symmetrically; x2 = -g/2, so
    ## it is the branch unknown that shows an unmasked anchor most
    ## clearly -- anchoring row 2 as well changes it outright.
    assert_allclose(x[0], 0.5, rtol=1e-6)
    assert_allclose(x[1], -0.5, rtol=1e-6)

    ## Every matrix the rung solver ever saw -- the pure ones among them,
    ## which is why the branch check is over all of them and the anchored
    ## check is over the ones that carry an anchor at all.
    assert rec.diags, 'the rung solver was never called'
    anchored = [d for d in rec.diags if d.any()]
    assert len(anchored) >= 5, len(rec.diags)     # a ladder ran, not one solve
    for d in rec.diags:
        assert d[2] == 0.0, d          # the branch row, NEVER anchored
    for d in anchored:
        assert d[0] > 0.0 and d[0] == d[1], d

