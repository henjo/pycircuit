"""Shooting tests: shooting dae.  Split out of test_analysis_shooting.py on
2026-09-27 (in its original order); shared helpers are in
`_shooting_fixtures.py`, the HDL elements in `_shooting_elements.py`.
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
from pycircuit.circuit.tests._shooting_fixtures import (_diode_fixture,
    _e5_circuit,
    _loss_osc,
    _rc_noisy)


def test_an_index_2_solve_that_converges_is_correct():
    """Index-2 circuits: what actually happens, measured on four topologies.

    An external review reported that index-2 circuits break the DEFAULT
    method, with "89% error on the v-source branch current", and proposed
    detecting index > 1 and refusing.  Reproduced and then measured
    further, and the conclusion does not survive:

        circuit      index>1   trap   euler   gear
        CV-loop        yes       F      F      T
        CV-loop2       yes       X      F      F
        LI-cutset      yes       T      T      T
        LI-cutset2     yes       X      X      F

    ⚠ INDEX > 1 IS NOT PREDICTIVE OF FAILURE.  Every one of those is index
    2 by the projector test (`N^T G N` singular, `N` a null basis of `C`),
    and on `LI-cutset` all three methods converge.  A refusal keyed on index
    would reject circuits that work.

    ⚠ AND GEAR IS NOT THE WORKAROUND.  It succeeds on `CV-loop` -- the one
    circuit the review tested -- and FAILS on two of the other three.  Had
    the refusal named `method='gear'` as the remedy, as proposed, it would
    have sent users to a method that does not generalise.

    ⚠ WHAT THIS PINS is the part that would matter if it broke: when an
    index-2 solve reports CONVERGED, its answer is right.  Measured against
    a settled transient over every state, algebraic ones included, the worst
    relative error was 6.6e-05 (CV-loop/gear) and 1.3e-05 / 1.7e-04 /
    7.4e-05 (LI-cutset trap/euler/gear).  The failures are LOUD --
    `converged = False` or an exception -- not silent wrong answers, which
    is why this is documented rather than refused.  A regression to a
    silently wrong algebraic variable is exactly what this test would catch.

    ⚠ WHAT IS PINNED IS PINNED FOR THE LINEAR CASE.  All four circuits are
    LINEAR, so `rank(C)` and the constraint structure are constant along the
    orbit.  On a NONLINEAR index-2 circuit both can vary with the operating
    point, so the index itself can change around the period -- and a solve
    could then be consistent at the points it checks and inconsistent
    between them, which is the one shape a converged-and-wrong answer could
    still take.  Offered as the gap rather than as a finding: no such
    circuit has been built here, and neither reviewer nor this test has
    evidence one exists in this tree.
    """
    import warnings
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    def cv_loop():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
        c['c1'] = C('a', 'b', c=1e-9)
        c['c2'] = C('b', gnd, c=1e-9)
        c['r'] = R('b', gnd, r=1e9)
        return c

    ## it really is index 2: N^T G N is singular for N a null basis of C
    cir = cv_loop()
    n = cir.n
    Cm = np.asarray(cir.C(np.zeros(n)), dtype=float)
    Gm = np.asarray(cir.G(np.zeros(n)), dtype=float)
    iref = cir.get_node_index(gnd)
    keep = [i for i in range(n) if i != iref]
    Cm, Gm = Cm[np.ix_(keep, keep)], Gm[np.ix_(keep, keep)]
    _U, sv, Vt = np.linalg.svd(Cm)
    d = int(np.sum(sv > len(keep) * sv[0] * np.finfo(float).eps))
    N = Vt[d:].T
    s2 = np.linalg.svd(N.T @ Gm @ N, compute_uv=False)
    ## ⚠ the projection can be identically ZERO, which is singular in the
    ## strongest sense -- so the ratio needs a guarded denominator rather
    ## than `s2[-1]/s2[0]`, which is 0/0 there.  Measured: this circuit
    ## gives exactly that.
    ratio = float(s2[-1] / max(s2[0], 1e-300))
    assert ratio < 1e-10, \
        'this circuit is no longer index 2 (sigma_min/sigma_max = %.3e), ' \
        'so the test has lost its subject' % ratio

    ## reference: a settled transient
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ## ⚠ a SHORT reference on purpose.  R*C here is ~0.5 s against a
        ## 1 ms period, so no affordable transient settles the slow DC
        ## drift -- but the quantity compared is the per-period AMPLITUDE,
        ## which the capacitive divider sets within a couple of periods.
        ## Sixty periods cost 52 s and bought no accuracy over eight.
        rt = Transient(cv_loop(), reltol=1e-10).solve(
            refnode=gnd, tend=per * 8, timestep=per / 200)
    t = np.asarray(rt.sweep_values, dtype=float).ravel()
    X = np.asarray(rt.x, dtype=float)
    W = X[:, t > t[-1] - per]
    ref = np.array([0.5 * (W[i].max() - W[i].min()) for i in range(W.shape[0])])

    pss = PSS(cv_loop(), method='gear', reltol=1e-8)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = pss.solve(period=per, timestep=per / 200, maxiterations=40)
    assert pss.converged, 'gear no longer converges on the CV-loop'
    Xp = np.asarray(res['tpss'].x, dtype=float)
    got = np.array([0.5 * (Xp[i].max() - Xp[i].min())
                    for i in range(Xp.shape[0])])
    k = min(len(got), len(ref))
    worst = np.max(np.abs(got[:k] - ref[:k])
                   / np.maximum(np.abs(ref[:k]), 1e-12))
    assert worst < 1e-3, \
        'a CONVERGED index-2 solve disagrees with a settled transient by ' \
        '%.3e over its states -- that is the silent wrong answer this ' \
        'exists to catch' % worst


def test_the_fill_is_skipped_entirely_without_an_algebraic_row():
    """The algebraic pattern gates the FILL; the `C^-T` conversion is not
    gated by it and runs regardless.

    ⚠ THAT DISTINCTION IS THE C^2 FIX. An earlier version skipped all
    per-sample work when there were no algebraic rows, which was right for
    the fill and WRONG for the conversion: a plain ODE circuit still needs
    `C^-T` whenever its capacitance is not 1 F. The pattern is asserted
    here because it still decides whether the fill runs.
    """
    cir, pss, _pac, _rs = _loss_osc('parallel', Q=8.0, npts=240)
    x0f = np.asarray(pss.waveform[1], dtype=float)[:, 0]
    rows, cols = pss._algebraic_adjoint_pattern(x0f)
    assert rows == [] and cols == [], \
        'the parallel-loss tank has no algebraic row; got rows=%r cols=%r' \
        % (rows, cols)

    cs, ps, _p2, _r2 = _loss_osc('series', Q=8.0, npts=240)
    x0s = np.asarray(ps.waveform[1], dtype=float)[:, 0]
    rows_s, cols_s = ps._algebraic_adjoint_pattern(x0s)
    assert len(rows_s) == 1 and len(cols_s) == 1, \
        'the series-loss tank has exactly one algebraic row and one ' \
        'algebraic state; got %r / %r' % (rows_s, cols_s)


def _measured_index(cir):
    """The index by direct computation: `N^T G N` singular for `N` a null

    ⚠⚠ THIS TEST IS COMPUTING A CANONICAL CHARACTERISTIC UNDER A LOCAL NAME.
    Hanke & März define `theta = dim(N ∩ S)` with `N = ker E` and
    `S = {z : Fz in im E}`.  For MNA with `E = C`, `F = G` and symmetric `C`,
    `Gz in im C` iff `Gz ⊥ ker C` iff `N^T G z = 0` -- so "`N^T G N` singular
    for `N` a null basis of `C`" IS `theta_0 > 0`.  (Relayed from a source
    reading, 2026-09-11.)

    ⚠ And `theta_0` is what makes the pole count come out.  Common Ground's
    Cor. c.pencil: for a constant pair, `theta_i` is the NUMBER OF JORDAN
    BLOCKS of order >= 2+i in the Weierstrass-Kronecker nilpotent, so
    `theta_0` COUNTS the index-2 constraints and `d = r - sum theta_i`.  That
    is the coefficient behind this tree's "rank(C) overcounts the finite poles
    by one per index-2 constraint": the one IS `theta_0`, and a circuit with
    two independent index-2 constraints would overcount by two.

    ⚠⚠ IT READS AN INDEX AND DOES NOT CERTIFY REGULARITY -- see
    `topological_index`'s docstring for the Brenan/Campbell/Petzold
    counterexample where every local pencil is regular and the solution family
    is infinite-dimensional.  No frozen-`t` test can see that.

    basis of `C` means index 2. The reference the topological criterion is
    gated against, and the one that overruled a relayed quote."""
    from pycircuit.circuit.analysis import remove_row_col
    n = cir.n
    Cm = np.asarray(cir.C(np.zeros(n)), dtype=float)
    Gm = np.asarray(cir.G(np.zeros(n)), dtype=float)
    irn = cir.get_node_index(gnd)
    Cr, Gr = remove_row_col((Cm, Gm), irn, circuit.numeric)
    Cr = np.asarray(Cr, dtype=float)
    Gr = np.asarray(Gr, dtype=float)
    _u, sv, vt = np.linalg.svd(Cr)
    tol = max(Cr.shape) * np.finfo(float).eps * (sv[0] if len(sv) else 1.0)
    if not np.any(sv <= tol):
        return 1
    N = vt[sv <= tol].T
    s2 = np.linalg.svd(N.T @ Gr @ N, compute_uv=False)
    if not len(s2):
        return 1
    ## ⚠ `N^T G N` IDENTICALLY ZERO IS THE MOST SINGULAR CASE, NOT THE LEAST.
    ## A first version guarded with `s2.max() > 0` and so returned 1 for the
    ## cv_loop fixture, whose block IS zero -- the reference disagreeing with
    ## the criterion because the REFERENCE was wrong.
    ##
    ## ⚠⚠ AND THE SECOND VERSION WAS BLIND WHENEVER `dim(ker C) == 1`, which
    ## is the commonest case.  It read `s2.min()/s2.max() < 1e-10`, and for a
    ## ONE-dimensional null space that ratio is IDENTICALLY 1 -- so the test
    ## could only ever fire on an EXACTLY zero block.  On a circuit whose
    ## algebraic block is 1e-20 it returned 1.  Blind, not trigger-happy, and
    ## no fixture here has a 1-D near-singular block so nothing failed.
    ##
    ## The fix needs an EXTERNAL scale, because a single singular value has no
    ## internal one to compare against.  `N^T G N` and `G` are both
    ## conductances and `N` is orthonormal, so `s2.min() <= tol * ||G||` is
    ## dimensionless AND free of the capacitance unit -- which is the property
    ## that matters, since an absolute tolerance here would give a false
    ## index-2 window widening as `1/||C||`.  Same reasoning as
    ## `algebraic_conditioning`, arrived at independently and kept inline: this
    ## is the reference the topological criterion is gated against, and a
    ## reference that calls shipped code is not independent of it.
    gscale = float(np.linalg.norm(Gr, 2))
    if gscale == 0.0:
        return 1
    return 2 if float(s2.min()) <= 1e-10 * gscale else 1


def _index_fixtures():
    circuit.default_toolkit = circuit.numeric
    per = 1e-3
    out = {}

    def cv():
        c = SubCircuit(); c.add_node('a'); c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
        c['c1'] = C('a', 'b', c=1e-9); c['c2'] = C('b', gnd, c=1e-9)
        c['r'] = R('b', gnd, r=1e9)
        return c
    out['cv_loop'] = (cv, 2, 'loop')

    def vac():
        c = SubCircuit(); c.add_node('a')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
        c['c1'] = C('a', gnd, c=1e-9); c['r'] = R('a', gnd, r=1e9)
        return c
    out['v_across_c'] = (vac, 2, 'loop')

    def lic():
        c = SubCircuit(); c.add_node('a'); c.add_node('b')
        c['is'] = IS(gnd, 'a', i=1e-3); c['l1'] = L('a', 'b', L=1e-3)
        c['c1'] = C('b', gnd, c=1e-9); c['r'] = R('b', gnd, r=1e3)
        return c
    out['li_cutset'] = (lic, 2, 'cutset')

    def conly():
        """⚠ a C-ONLY loop: relayed theory says index 2, MEASUREMENT says 1."""
        c = SubCircuit()
        for nn in ('a', 'b', 'cc'):
            c.add_node(nn)
        c['c1'] = C('a', 'b', c=1e-9); c['c2'] = C('b', 'cc', c=1e-9)
        c['c3'] = C('cc', 'a', c=1e-9); c['r'] = R('a', gnd, r=1e6)
        return c
    out['c_only_loop'] = (conly, 1, None)

    def both():
        """⚠ a C-only loop AND a C-V loop, which catches a union-find that
        reports only the first closing edge."""
        c = SubCircuit()
        for nn in ('a', 'b', 'cc', 'd'):
            c.add_node(nn)
        c['c1'] = C('a', 'b', c=1e-9); c['c2'] = C('b', 'cc', c=1e-9)
        c['c3'] = C('cc', 'a', c=1e-9)
        c['vs'] = VSin('d', gnd, va=1.0, freq=1.0 / per)
        c['c4'] = C('d', gnd, c=1e-9)
        c['r'] = R('a', gnd, r=1e6); c['r2'] = R('d', gnd, r=1e6)
        return c
    out['both_loops'] = (both, 2, 'loop')

    def i1rc():
        c = SubCircuit(); c.add_node('a'); c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
        c['r'] = R('a', 'b', r=1e3); c['c1'] = C('b', gnd, c=1e-9)
        return c
    out['index1_rc'] = (i1rc, 1, None)

    def i1rlc():
        c = SubCircuit(); c.add_node('a'); c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
        c['r1'] = R('a', 'b', r=1e3); c['l1'] = L('a', 'b', L=1e-3)
        c['c1'] = C('b', gnd, c=1e-9)
        return c
    out['index1_rlc'] = (i1rlc, 1, None)
    return out


def test_the_topological_index_agrees_with_the_measured_one():
    """⚠ Estevez Schwarz & Tischendorf's criterion, gated against a direct
    computation on the MNA matrices.

    Index 2 IFF a C-V loop or an L-I cutset. ⚠⚠ IT IS A DIAGNOSTIC, NOT A
    REFUSAL -- C4 stays closed, because `index > 1` is not predictive of
    convergence. What this buys is that a failure can NAME the offending
    elements, which is the criterion's own design goal.

    ⚠⚠⚠ AND C-ONLY LOOPS ARE EXCLUDED HERE AGAINST THE RELAYED QUOTE. A
    first version counted them and disagreed with the measurement on three
    separate C-only topologies. The arithmetic is checkable: a grounded
    capacitor ring has `det C = c1 c2 + c1 c3 + c2 c3 != 0`, not even a DAE;
    a floating triangle has `C` singular but `N^T G N = (1/R)/3 != 0`, so the
    constraint is uniquely solvable and the index is 1. A C-only loop makes
    `C` singular WITHOUT making the index 2.
    """
    from pycircuit.circuit.shooting import topological_index
    for name, (build, expect, where) in _index_fixtures().items():
        cir = build()
        idx, info = topological_index(cir)
        assert idx == _measured_index(cir), \
            '%s: topological index %d against a measured %d' \
            % (name, idx, _measured_index(cir))
        assert idx == expect, '%s: expected index %d, got %d' % (name, expect, idx)
        if where == 'loop':
            assert info['loop'], '%s: index 2 with no loop named' % name
            assert any(info['kinds'][nm] == 'V' for nm in info['loop']), \
                '%s: a C-V loop must contain a voltage source; got %r' \
                % (name, info['loop'])
        elif where == 'cutset':
            assert info['cutset'], '%s: index 2 with no cutset named' % name
        else:
            assert not info['loop'] and not info['cutset'], \
                '%s: index 1 but something was named: %r / %r' \
                % (name, info['loop'], info['cutset'])


def test_the_index_criterion_localises_and_flags_what_it_cannot_classify():
    """The criterion's whole point is LOCALISATION, and its honesty is the
    `provisional` flag.

    ⚠ The theorem excludes CONTROLLED SOURCES, so a netlist carrying one
    gets an answer WITH a caveat rather than silence -- a named assumption
    beats a withheld verdict.
    """
    from pycircuit.circuit.shooting import topological_index
    circuit.default_toolkit = circuit.numeric
    fx = _index_fixtures()
    cir = fx['both_loops'][0]()
    _idx, info = topological_index(cir)
    ## the C-V loop is on node d -- vs with c4 -- NOT the a/b/cc triangle
    assert set(info['loop']) == {'vs', 'c4'}, \
        'the loop must be localised to the offending elements; got %r' \
        % (info['loop'],)
    assert not info['provisional'] and not info['unclassified']

    ## and a controlled source makes the verdict provisional
    mu = 1.0 / (2 * np.pi * 8.0)
    c2 = SubCircuit()
    c2.add_node('v')
    c2['C'] = C('v', gnd, c=1.0)
    c2['L'] = L('v', gnd, L=1.0)
    c2['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: mu * (u - u ** 3 / 3.0))
    _i2, info2 = topological_index(c2)
    assert info2['provisional'], 'a controlled source must be flagged'
    assert any('BSource' in u for u in info2['unclassified'])


def test_noise_in_the_constraints_is_detected():
    """⚠ Winkler's index-1 SDAE precondition, `im B subset im C`.

    A resistor's noise injected at a node with no capacitance puts noise
    into a CONSTRAINT row, which makes the circuit an SDAE WITH DIRECT NOISE
    -- outside the class the theory covers, and the reason that node has no
    finite variance for `K_orb` to report (roadmap section 0j).
    """
    from pycircuit.circuit.shooting import noise_enters_constraints
    from pycircuit.circuit.analysis import remove_row_col
    import warnings
    circuit.default_toolkit = circuit.numeric

    def check(kind):
        cir, pss, pac, _rs = _loss_osc(kind, Q=8.0, npts=240)
        cy = np.real(pac._cy_reduced(pss, 2 * np.pi / float(pss.period)))
        X = np.asarray(pss.waveform[1], dtype=float)
        Cr, = remove_row_col((np.asarray(cir.C(X[:, 0]), dtype=float),),
                             pss.irefnode, pss.toolkit)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return noise_enters_constraints(np.asarray(Cr, dtype=float), cy)

    bad_p, res_p = check('parallel')
    bad_s, res_s = check('series')
    assert not bad_p, \
        'the parallel-loss tank puts its noise on a differential row; ' \
        'residual %.3e' % res_p
    assert bad_s, \
        'the series-loss tank puts its noise on an ALGEBRAIC row and must ' \
        'be flagged; residual %.3e' % res_s
    assert res_s > 0.5, \
        'the violation is total, not marginal -- the whole CY column lies ' \
        'outside im C; got %.3e' % res_s


def test_v_loops_and_i_cutsets_are_reported_as_ILL_POSED_not_as_index_2():
    """⚠ A DIFFERENT CATEGORY, and the distinction is the point.

    A loop of voltage sources over-determines KVL; a cutset of current
    sources over-determines KCL. Either way the MNA system is STRUCTURALLY
    SINGULAR and has no solution at all, barring an exact cancellation --
    it does not have a higher index. The index criterion presumes a
    well-posed network.

    Reporting these as "index 2" would send a reader hunting a solver
    problem when the netlist is the error, so they come back separately and
    `index` keeps its own meaning.
    """
    from pycircuit.circuit.shooting import topological_index
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    ## two voltage sources in parallel -- a V-only loop
    c = SubCircuit()
    c.add_node('a')
    c['vs1'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
    c['vs2'] = VSin('a', gnd, va=2.0, freq=1.0 / per)
    c['r'] = R('a', gnd, r=1e3)
    idx, info = topological_index(c)
    assert info['ill_posed'], 'two parallel voltage sources form a V loop'
    assert info['v_loop'], 'the offending source must be named'
    ## ⚠ and there is NO index to report -- the DAE index presumes a solvable
    ## system, and returning 2 here would point at the solver
    assert idx is None, 'an ill-posed netlist must not be given an index'
    assert not info['loop'] and not info['cutset'], \
        'a pure-V loop must not be reported as a C-V loop: %r' % (info['loop'],)

    ## a node reachable only through current sources -- an I-only cutset
    c2 = SubCircuit()
    c2.add_node('a'); c2.add_node('b')
    c2['is1'] = IS(gnd, 'a', i=1e-3)
    c2['is2'] = IS('a', 'b', i=1e-3)
    c2['r'] = R('b', gnd, r=1e3)
    i2, info2 = topological_index(c2)
    assert info2['ill_posed'], 'node a is isolated by current sources'
    assert i2 is None
    assert set(info2['i_cutset']) >= {'is1', 'is2'}, \
        'both current sources bound the cutset; got %r' % (info2['i_cutset'],)

    ## and a well-posed netlist must NOT be flagged, or this shows nothing
    c3 = SubCircuit()
    c3.add_node('a'); c3.add_node('b')
    c3['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
    c3['r'] = R('a', 'b', r=1e3)
    c3['c1'] = C('b', gnd, c=1e-9)
    _i3, info3 = topological_index(c3)
    assert not info3['ill_posed'], \
        'a plain RC ladder must not be flagged ill-posed: %r / %r' \
        % (info3['v_loop'], info3['i_cutset'])


@pytest.mark.slow
def test_radau_keeps_order_on_differential_and_loses_two_on_algebraic_index2():
    """On an INDEX-2 MNA, Radau IIA(3) holds classical order in the
    DIFFERENTIAL components and drops to 3 in the ALGEBRAIC ones.

    Hairer, Lubich & Roche (LNM 1409, 1989) sec. 5, Thm 5.9: with `det A != 0`
    and stiff accuracy -- exactly the pair Radau IIA(3) has -- "there is NO
    ORDER REDUCTION IN THE Y-COMPONENT", while the z-component estimate "is in
    general optimal", i.e. NOT improved.  Measured here, capacitor-loop
    fixture, analytic reference::

        npts   err(differential)   ord     err(algebraic)   ord
           5   2.4738e-04          --      1.1349e-07       --
          10   4.4801e-06         5.79     1.3163e-08      3.11
          20   1.0875e-07         5.36     1.4153e-09      3.22
          40   3.0134e-09         5.17     1.6408e-10      3.11
          80   8.8831e-11         5.08     1.9753e-11      3.05
         160   2.6968e-12         5.04     2.4232e-12      3.03

    ⚠⚠ WHY THIS TEST EXISTS: **outcome (b) looks exactly like a tableau bug
    and is not one.**  A reviewer measuring only the algebraic component would
    see Radau IIA(3) -- order 5 -- converging at 3 and open a defect against
    the tableau.  The reduction is the DAE INDEX.  Splitting the error by
    subspace is what makes the two distinguishable.

    ⚠ THREE WAYS THIS MEASUREMENT CAN LIE, all closed here rather than argued:

      * **the circuit might not be index 2**, in which case both components
        sit at classical order and the test says nothing.  Asserted first, via
        `sigma_min/sigma_max` of `N^T G N` on a null basis of `C`;
      * **the reference might be wrong**, which shows up as a CONSTANT error
        and reads as "order 0".  The first attempt did exactly that -- 1.5747
        at every npts -- because `analysis='ac'` gave a phasor in the COSINE
        convention while `VSin` drives a SINE, and `vac` is a separate
        parameter defaulting to 1.  ⚠ The residual check that passed it,
        `|C jwX + G X + U| = 6e-22`, only verified the LINEAR SOLVE.  The
        phasor is now derived from the time function the transient actually
        integrates, and validated against `C xdot + G x + u(t)`;
      * **the sweep might be floor-limited**, which also reads as order 0.
        The first fixture had `tau = 1 s` against a 1 ms period -- nearly
        quasi-static -- so the differential error started at 1.9e-12, already
        at the floor, and flattened.  `tau ~ PER/10` puts the coarsest point
        at 2.5e-04, eight decades above it.

    ⚠ The fixture is LINEAR and NOISELESS by choice.  A noisy sweep has three
    regimes -- deterministic-dominated ~2, noise-dominated 1/2, floor-limited 0
    -- so refining `h` makes the observed order WORSE and no order statement
    is meaningful without measuring the floor first.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    per = 1e-3
    w = 2 * np.pi / per

    def cv_loop():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
        c['c1'] = C('a', 'b', c=1e-9)
        c['c2'] = C('b', gnd, c=1e-9)
        c['r'] = R('b', gnd, r=1e5)
        return c

    cir = cv_loop()
    n = cir.n
    iref = cir.get_node_index(gnd)
    keep = [i for i in range(n) if i != iref]
    Cm = np.asarray(cir.C(np.zeros(n)), dtype=float)[np.ix_(keep, keep)]
    Gm = np.asarray(cir.G(np.zeros(n)), dtype=float)[np.ix_(keep, keep)]
    m = len(keep)

    ## (1) it must BE index 2, or both components sit at classical order and
    ## the measurement is vacuous
    _U, sv, Vt = np.linalg.svd(Cm)
    d = int(np.sum(sv > m * sv[0] * np.finfo(float).eps))
    N = Vt[d:].T
    s2 = np.linalg.svd(N.T @ Gm @ N, compute_uv=False)
    assert float(s2[-1] / max(s2[0], 1e-300)) < 1e-10, \
        'the fixture is not index 2, so an order split would prove nothing'
    Vdiff, Valg = Vt[:d].T, Vt[d:].T
    assert Vdiff.shape[1] >= 1 and Valg.shape[1] >= 1

    ## (2) the reference, from the SOURCE FUNCTION the transient integrates
    NS = 2048
    ts = np.arange(NS) * per / NS
    Us = np.array([np.asarray(cir.u(t, analysis='tran'),
                              dtype=float)[keep] for t in ts])
    Uph = (2.0 / NS) * np.sum(Us * np.exp(-1j * w * ts)[:, None], axis=0)
    X = np.linalg.solve(Gm + 1j * w * Cm, -Uph)

    def x_exact(t):
        return np.real(X * np.exp(1j * w * t))

    ## ⚠ validated in the TIME DOMAIN, not against the system it was built from
    chk = 0.0
    for t in (0.0, per * 0.137, per * 0.41, per * 0.76):
        xd = np.real(1j * w * X * np.exp(1j * w * t))
        ut = np.asarray(cir.u(t, analysis='tran'), dtype=float)[keep]
        chk = max(chk, float(np.max(np.abs(Cm @ xd + Gm @ x_exact(t) + ut))))
    assert chk < 1e-12, 'the analytic reference does not satisfy the ODE: %.3e' % chk

    ## (3) the sweep
    errs = []
    for npts in (5, 10, 20, 40, 80):
        p = PSS(cv_loop(), method='radau', reltol=1e-13)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=per, timestep=per / npts, maxiterations=40)
        assert p.converged, 'radau did not converge at npts=%d' % npts
        t = np.asarray(p.waveform[0], dtype=float)
        Xr = np.asarray(p.waveform[1], dtype=float)[keep, :]
        E = np.array([Xr[:, k] - x_exact(t[k]) for k in range(len(t))])
        errs.append((float(np.max(np.abs(E @ Vdiff))),
                     float(np.max(np.abs(E @ Valg)))))

    ## ⚠ the coarsest point must be FAR above the floor, or the orders below
    ## are measuring roundoff.  This is the check the first fixture failed.
    assert errs[0][0] > 1e-6, \
        'the differential error starts at %.3e -- too close to the floor for ' \
        'the sweep to resolve an order' % errs[0][0]

    od = [np.log2(errs[i][0] / errs[i + 1][0]) for i in range(len(errs) - 1)]
    oa = [np.log2(errs[i][1] / errs[i + 1][1]) for i in range(len(errs) - 1)]
    assert od[-1] > 4.5, \
        'the DIFFERENTIAL components should keep classical order (~5), got ' \
        '%.2f (all: %s)' % (od[-1], ['%.2f' % z for z in od])
    assert 2.6 < oa[-1] < 3.5, \
        'the ALGEBRAIC components should drop to ~3, got %.2f (all: %s)' \
        % (oa[-1], ['%.2f' % z for z in oa])
    assert od[-1] - oa[-1] > 1.5, \
        'the whole point is the SPLIT: differential %.2f against algebraic ' \
        '%.2f' % (od[-1], oa[-1])


def test_the_topological_index_agrees_with_an_incidence_RANK_criterion():
    """A SECOND, ALGEBRAICALLY INDEPENDENT route to the same verdict.

    `topological_index` implements Estevez Schwarz & Tischendorf by UNION-FIND
    on the netlist graph.  Lamour, Marz & Tischendorf's Lemma 3.45 (p. 243,
    relayed by the docs session, 2026-09-09) states the same two criteria as
    RANK CONDITIONS on incidence matrices:

        [A_C A_R A_V] has full ROW rank   iff no cutset of inductances and
                                              current sources only
        Q_C^T A_V   has full COLUMN rank  iff no loop of capacitances and
                                              voltage sources with at least
                                              one voltage source

    with `Q_C` a basis of `ker(A_C^T)`.  ⚠ THE POINT IS THE INDEPENDENCE: a
    second graph algorithm could share a traversal bug with the first, and
    linear algebra over the incidence matrices cannot.  Both routes are run
    on four topologies and each rank condition is made to FAIL on the one it
    is meant to catch -- an all-pass comparison would prove nothing.

    ⚠ The rank route is a CHECK, not a replacement: `topological_index` stays
    authoritative (it names the offending elements, which a rank cannot), and
    both rest on the same hypothesis -- Theorem 3.47's "let all current and
    voltage sources be INDEPENDENT", which is this tree's controlled-source
    caveat reached from a second source.
    """
    from pycircuit.circuit.shooting import diagnostics as _diag
    circuit.default_toolkit = circuit.numeric

    def incidence(cir):
        nodes = list(cir.nodes)
        nn = len(nodes)
        iref = cir.get_node_index(gnd)
        keep = [i for i in range(nn) if i != iref]
        pos = {n: k for k, n in enumerate(keep)}
        nmap = cir.elementnodemap
        cols = {'C': [], 'L': [], 'V': [], 'I': [], 'R': []}
        for name in cir.elements:
            cls = type(cir[name]).__name__
            kind = ('C' if cls in _diag._TI_CAPACITIVE else
                    'V' if cls in _diag._TI_VOLTAGE else
                    'L' if cls in _diag._TI_INDUCTIVE else
                    'I' if cls in _diag._TI_CURRENT else
                    'R' if cls in _diag._TI_RESISTIVE else '?')
            if kind == '?':
                continue
            idx = sorted(set(int(i) for i in np.asarray(nmap[name]).ravel()
                             if int(i) < nn))
            ends = [i for i in idx if i != iref]
            col = np.zeros(len(keep))
            if len(ends) == 2:
                col[pos[ends[0]]] = 1.0
                col[pos[ends[1]]] = -1.0
            elif len(ends) == 1:
                col[pos[ends[0]]] = 1.0
            else:
                continue
            cols[kind].append(col)
        return ({k: (np.array(v).T if v else np.zeros((len(keep), 0)))
                 for k, v in cols.items()}, len(keep))

    def rank_index(cir):
        A, m = incidence(cir)
        CRV = np.hstack([A['C'], A['R'], A['V']])
        no_cutset = int(np.linalg.matrix_rank(CRV)) == m if CRV.size else m == 0
        if A['C'].size:
            _u, sv, vt = np.linalg.svd(A['C'].T, full_matrices=True)
            rk = int(np.sum(sv > max(A['C'].shape) * sv[0] * np.finfo(float).eps))
            QC = vt[rk:].T
        else:
            QC = np.eye(m)
        ## ⚠ AN EMPTY `Q_C` MAKES THE CONDITION FAIL, NOT PASS.  When `A_C`
        ## has full row rank there is no kernel, `Q_C` is `m x 0`, and
        ## `Q_C^T A_V` is `0 x nV`: a matrix with no rows cannot have full
        ## COLUMN rank unless it has no columns either.  Guarding with
        ## `if QC.size: ... else: no_loop = True` inverts exactly the C-V
        ## loop this is meant to catch -- measured, it returned index 1 on
        ## the C-V fixture.
        M = QC.T @ A['V']
        no_loop = (M.shape[1] == 0) or (int(np.linalg.matrix_rank(M)) == M.shape[1])
        return (1 if (no_cutset and no_loop) else 2), no_cutset, no_loop

    def tank():
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        return c

    def tank_r_only():
        ## ⚠⚠ THE CASE THIS TEST CLAIMED TO COVER AND NEVER DID.  `tank` and
        ## `tank+RC` both carry a capacitor at EVERY node, so both are index
        ## 0; with only those and the two index-2 fixtures, the traversal side
        ## of this cross-check spanned {0, 0, 2, 2} and the `1` in
        ## `max(traversal, 1) == rank` was ONLY ever reached through the floor.
        ## A genuine index-1 circuit had never been compared against the rank
        ## criterion here.  Node `w` carries no capacitor: rank C = 2 of 3.
        c = tank()
        c.add_node('w')
        c['Rs'] = R('v', 'w', r=100.0)
        c['Rl'] = R('w', gnd, r=1e3)
        return c

    def tank_rc():
        c = tank()
        c.add_node('w')
        c['Rs'] = R('v', 'w', r=100.0)
        c['Cs'] = C('w', gnd, c=1e-3)
        return c

    def li_cutset():
        c = SubCircuit()
        c.add_node('v')
        c.add_node('w')
        c['C'] = C('v', gnd, c=1.0)
        c['L1'] = L('v', 'w', L=0.5)
        c['L2'] = L('w', gnd, L=0.5)
        return c

    def cv_loop():
        c = SubCircuit()
        c.add_node('v')
        c.add_node('b')
        c['C'] = C('v', gnd, c=0.5)
        c['C1'] = C('v', 'b', c=0.5)
        c['Vb'] = VS('b', gnd, v=1.0)
        c['L'] = L('v', gnd, L=1.0)
        return c

    seen = {}
    spanned = set()
    for label, build in (('tank', tank), ('tank+RC', tank_rc),
                         ('tank+R-only', tank_r_only),
                         ('L-I cutset', li_cutset), ('C-V loop', cv_loop)):
        cir = build()
        traversal, _info = topological_index(cir)
        rank, no_cutset, no_loop = rank_index(cir)
        ## ⚠⚠ COMPARED ON THE 1-vs-2 AXIS, because LEMMA 3.45 IS ALSO FLOORED
        ## AT 1.  Its two conditions -- `[A_C A_R A_V]` full row rank iff no
        ## L-I cutset, `Q_C^T A_V` full column rank iff no C-V loop --
        ## distinguish 1 from 2 and cannot express 0, the same limitation
        ## Estevez Schwarz & Tischendorf has.  When the index-0 rung was added
        ## to `topological_index` (2026-09-11) this cross-check FAILED on the
        ## `tank` fixture, reading (0, 1), and it was RIGHT TO: the tank is
        ## C + L with no voltage source, its reduced `C` is NONSINGULAR (2x2,
        ## rank 2, sigma_min 1.0), and it is genuinely index 0.  The traversal
        ## is correct and the rank criterion cannot see it.
        ## ⚠ The cross-check for THAT rung is a different instrument,
        ## `sigma_min(C)`, in `test_topological_index_reports_the_index_0_rung
        ## _and_agrees_with_the_rank_test`.  Folding it in here would compare
        ## two floored criteria and prove nothing.
        assert max(traversal, 1) == rank, (label, traversal, rank)
        seen[label] = (rank, no_cutset, no_loop)
        spanned.add(traversal)
    ## ⚠⚠ THE FIXTURE SET MUST SPAN ALL THREE INDICES, asserted rather than
    ## assumed.  A FLOORED INSTRUMENT COLLAPSES TWO VALUES INTO ONE NAME, so
    ## a case can go missing without any assertion failing -- every assertion
    ## here was written against the floored criterion and was therefore
    ## correct while index 1 was untested.  No test catches that; only
    ## reading does.  This line is what makes the next rung fail loudly
    ## instead of a case vanishing into the floor.
    ## (Consequence-to-look-for named by docs-46, 2026-09-11, which found the
    ## same hole in its own pole-count check: not just a wrong label, a
    ## MISSING CASE.)
    assert spanned == {0, 1, 2}, spanned
    ## ⚠ A STALE LABEL FROM BEFORE THE INDEX-0 RUNG: this line used to read
    ## "the two index-1 topologies", and `tank` is index 0 -- the name was
    ## minted while `topological_index` was floored at 1.  `seen` holds RANK
    ## values, and the rank criterion is ALSO floored, so the 1 here is the
    ## floor and not the tank's index.  The assertion is right; the label was
    ## describing the instrument's ceiling as if it were the circuit.
    ## (Class flagged by docs-46, 2026-09-11, which hit the same stale name on
    ## its own 'index-1 tank' fixture.)
    ## the three topologies the RANK criterion reads as 1 pass BOTH conditions
    ## -- two of them index 0 by the floor, and `tank+R-only` a GENUINE index 1
    assert seen['tank'] == (1, True, True) and seen['tank+RC'] == (1, True, True), seen
    assert seen['tank+R-only'] == (1, True, True), seen['tank+R-only']
    ## and each index-2 one fails exactly the condition named for it -- without
    ## this the agreement above could be four passes of a criterion that never
    ## fires
    assert seen['L-I cutset'] == (2, False, True), seen['L-I cutset']
    assert seen['C-V loop'] == (2, True, False), seen['C-V loop']


def test_topological_index_reports_the_index_0_rung_and_agrees_with_the_rank_test():
    """⚠⚠ THE FUNCTION USED TO BE FLOORED AT 1 AND ANSWERED 1 FOR AN IMPLICIT
    ODE, SILENTLY.

    Estevez Schwarz & Tischendorf's criterion is "index 2 if and only if the
    network contains a C-V loop or an L-I cutset, otherwise 1" — it cannot
    return 0. Theorem 3.47's index-0 case (a capacitive path from every node to
    datum AND no voltage sources) was recorded in the docstring as something
    the theory "adds" and was never implemented, so a van der Pol read 1.

    How it was found is the part worth keeping: a NUMERICAL probe read index 0
    for that circuit and was CLAMPED to agree with this function — agreeing
    for the wrong reason, against a reference that could not represent the
    answer. The clamp is gone and the rung is implemented.

    ⚠ THE INDEPENDENT CROSS-CHECK is the rank condition, and it is what makes
    this more than a restatement: **index 0 IFF the reduced `C` is
    NONSINGULAR** — that is exactly what an implicit ODE is. Asserted in both
    directions on all three fixtures.

    ⚠ AN INDUCTOR DOES NOT SPOIL INDEX 0, which a reading of the theorem's
    wording alone might get wrong: the flux term makes that branch row
    differential, not algebraic. The van der Pol carries one and is the case
    that shows it.
    """
    import os
    import sys
    import warnings
    import numpy as np
    from pycircuit.circuit.circuit import gnd as _gnd, defaultepar
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.analysis import remove_row_col
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.shooting import topological_index
    from pycircuit.circuit.integrator import Gear2Integrator
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    '..', '..', '..', 'benchmarks'))
    from noise_floor_sources_off import vdp
    from pycircuit.circuit.tests.test_stage_predictor import (_expg_fixture,
                                                              PER)
    from pycircuit.circuit.tests.test_glm import _cv_loop

    def c_is_nonsingular(cir):
        tr = Transient(cir, integrator=Gear2Integrator(), reltol=1e-12)
        iref = cir.get_node_index(_gnd)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            x = np.asarray(DC(cir, refnode=_gnd).solve().x,
                           dtype=float).ravel()
        C = np.asarray(cir.C(x, defaultepar), dtype=float)
        (Cr,) = remove_row_col((C,), iref, tr.toolkit)
        Cr = np.asarray(Cr, dtype=float)
        sv = np.linalg.svd(Cr, compute_uv=False)
        return int((sv > 1e-12 * max(sv.max(), 1e-300)).sum()) == Cr.shape[0]

    cases = [(lambda: vdp()[0], 'van der Pol', 0),
             (lambda: _expg_fixture(PER), 'ExpG', 1),
             (lambda: _cv_loop(1e-3), 'C-V loop', 2)]
    for build, name, want in cases:
        cir = build()
        idx, info = topological_index(cir)
        assert idx == want, (name, idx, want)
        ## and the two routes must agree, in BOTH directions
        assert (idx == 0) == c_is_nonsingular(build()), (name, idx)

    ## the index-0 verdict reports what it tested
    _i, info = topological_index(vdp()[0])
    assert info['cap_path_to_datum'] is True, info

    ## ⚠ A PINNED LIMITATION, not an aspiration: `_TI_CAPACITIVE` matches the
    ## built-in `C` class only, so a BEHAVIOURAL charge is INVISIBLE to the
    ## topological route -- including to the C-V loop test, which predates
    ## this.  Measured on the branch fixture, whose `CubicCap` carries a real
    ## `q`: `cap_path_to_datum` is False because the capacitor cannot be seen.
    ## Its agreement with the rank test there is COINCIDENTAL -- topological
    ## says 1 because it is blind, rank says 1 because `C` vanishes at that
    ## operating point -- which is why the agreement above is asserted only on
    ## fixtures both routes can see.
    from branch_selection import build as branch_build
    cir = branch_build(-1.0)
    assert hasattr(cir['cq'], 'q')
    _i, info = topological_index(cir)
    assert info['cap_path_to_datum'] is False, info


def _ac_fixtures():
    """RC (index 1), C-V loop (index 2), and a circuit with NO algebraic
    block -- one per verdict `algebraic_conditioning` can return."""
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    def rc(k=1.0):
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VSin('a', gnd, va=0.8, freq=1.0 / per)
        c['rs'] = R('a', 'b', r=50.0)
        c['cl'] = C('b', gnd, c=1e-9 * k)
        c['rl'] = R('b', gnd, r=1e4)
        return c

    def cv():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
        c['c1'] = C('a', 'b', c=1e-9)
        c['c2'] = C('b', gnd, c=2e-9)
        c['rl'] = R('b', gnd, r=1e4)
        return c

    def ode():
        ## Current-driven, capacitor to datum on the only node: `C` is
        ## nonsingular, so there is no algebraic block at all.
        c = SubCircuit()
        c.add_node('a')
        c['i'] = IS(gnd, 'a', i=1e-3)
        c['c'] = C('a', gnd, c=1e-9)
        c['r'] = R('a', gnd, r=1e3)
        return c

    return rc, cv, ode


def _explicit_block_sigma(cir):
    """`sigma_min(Z^T G N)` the long way -- `N = ker C`, `Z = ker C^T`, both
    by SVD.  The reference `algebraic_conditioning` must reproduce WITHOUT
    forming either basis."""
    from pycircuit.circuit.analysis import remove_row_col
    import pycircuit.circuit.analysis as _an
    n = cir.n
    Cm = np.asarray(cir.C(np.zeros(n), defaultepar), dtype=float)
    Gm = np.asarray(cir.G(np.zeros(n), defaultepar), dtype=float)
    Cm, Gm = [np.asarray(m, dtype=float) for m in
              remove_row_col((Cm, Gm), cir.get_node_index(gnd), _an.numeric)]
    u, sv, vt = np.linalg.svd(Cm)
    tol = max(Cm.shape) * np.finfo(float).eps * sv[0]
    if not np.any(sv <= tol):
        return None
    N = vt[sv <= tol].T
    Z = u[:, sv <= tol]
    sb = np.linalg.svd(Z.T @ Gm @ N, compute_uv=False)
    return float(sb.min()) if len(sb) else None


def test_algebraic_conditioning_reproduces_the_block_without_forming_it():
    """`sigma_min(C + hG)/h -> sigma_min(d g_2/d y)`, and the verdict agrees
    with the topological index on all three fixtures.

    The point of the routine is that it needs NO index-1 `(x, y)` splitting
    and no null basis: two SVDs at different `h` from `C` and `G` alone.  So
    the gate is that it reproduces the explicitly formed `Z^T G N` block.
    """
    rc, cv, ode = _ac_fixtures()

    sigma, info = algebraic_conditioning(rc())
    assert info['verdict'] == 'well-conditioned', info['verdict']
    ## the documented contract: `spread - 1` bounds the relative error
    assert abs(sigma - _explicit_block_sigma(rc())) <= \
        (info['spread'] - 1.0) * _explicit_block_sigma(rc())
    ref = _explicit_block_sigma(rc())
    ## ⚠ THE BOUND IS SET BY THE LADDER, NOT BY WHAT PASSES.  Three readouts
    ## were tried and the numbers are the argument: the MEDIAN of the ladder
    ## gives 1.1e-09, its LAST point 2.0e-13 here but 1.0e-04 on the swallowed
    ## -capacitor fixture below, and the FLATTEST 3-POINT WINDOW 2.0e-12 here
    ## and 2.0e-10 there.  1e-10 passes only the last of those.
    assert abs(sigma - ref) <= 1e-10 * ref, (sigma, ref)
    assert topological_index(rc())[0] == 1

    ## An index-2 circuit is exactly `theta_0 > 0`, i.e. a SINGULAR block.
    sigma, info = algebraic_conditioning(cv())
    assert info['verdict'] == 'singular', info['verdict']
    assert sigma == 0.0
    assert _explicit_block_sigma(cv()) == 0.0
    assert topological_index(cv())[0] == 2

    ## `C` nonsingular: nothing to condition, and the ratio would grow like
    ## `1/h` forever.  Reporting a number here would be reporting noise.
    sigma, info = algebraic_conditioning(ode())
    assert info['verdict'] == 'no-algebraic-block', info['verdict']
    assert sigma is None
    assert _explicit_block_sigma(ode()) is None


def test_algebraic_conditioning_is_invariant_to_the_CAPACITANCE_UNIT():
    """⚠⚠ THE HAZARD THIS EXISTS TO EXCLUDE.  A rank test on these blocks with
    an ABSOLUTE tolerance smears the index-2 crossing into a false window,
    because the relevant singular value goes as `|G_22| * ||C|| / ||G||` --
    so the window widens as `1/||C||` and at picofarads it is enormous.  A
    circuit does not change its index when its capacitors are restated in
    different units, so the verdict must not move either.  (Hazard relayed
    from the docs session, 2026-09-11, which hit it in its own first version.)

    Two things make it not apply here and BOTH are load-bearing: the null
    space of `C` is taken by a RELATIVE tolerance, and the verdict is the
    FLATNESS of the ratio, which is scale-free, rather than a magnitude
    compared against a fixed number.
    """
    rc, _cv, _ode = _ac_fixtures()
    seen = []
    for k in (1e6, 1e3, 1.0, 1e-3, 1e-6):
        sigma, info = algebraic_conditioning(rc(k))
        assert info['verdict'] == 'well-conditioned', (k, info['verdict'])
        seen.append(sigma)
    assert max(seen) / min(seen) <= 1 + 1e-9, seen


def test_algebraic_conditioning_says_no_window_rather_than_a_plausible_number():
    """⚠ THE WINDOW CAN BE EMPTY, AND THAT MUST NOT BE ROUNDED TO AN ANSWER.

    The ratio has to be read above the roundoff floor (`eps*||C||`) and below
    the turn (`sigma_r(C)/sigma`), and a badly conditioned `C` leaves nothing
    in between.  A routine that fits a plateau anyway would return a confident
    number built out of roundoff.

    ⚠ The response is NOT monotone in the conditioning of `C`, and that reads
    like a bug until you see why.  Widening the capacitance spread to 1e10 and
    1e14 closes the window, but at 1e16 the small capacitor drops BELOW the
    relative rank tolerance, is correctly reclassified as part of the
    algebraic block, and the answer comes back -- as `1/r2`, which is the
    right answer for a node whose capacitance is numerically absent.
    """
    per = 1e-3

    def ladder(spread):
        c = SubCircuit()
        for nm in ('a', 'b', 'd'):
            c.add_node(nm)
        c['vs'] = VSin('a', gnd, va=0.8, freq=1.0 / per)
        c['rs'] = R('a', 'b', r=50.0)
        c['cl'] = C('b', gnd, c=1e-9)
        c['rl'] = R('b', gnd, r=1e4)
        c['r2'] = R('b', 'd', r=1e3)
        c['c2'] = C('d', gnd, c=1e-9 / spread)
        return c

    sigma, info = algebraic_conditioning(ladder(1e6))
    assert info['verdict'] == 'well-conditioned', info['verdict']

    for spread in (1e10, 1e12):
        sigma, info = algebraic_conditioning(ladder(spread))
        assert info['verdict'] == 'no-window', (spread, info['verdict'])
        assert sigma is None

    ## and the non-monotone tail: the capacitor vanishes into the null space
    sigma, info = algebraic_conditioning(ladder(1e16))
    assert info['verdict'] == 'well-conditioned', info['verdict']
    ## ⚠ AND IT IS NOT EXACT HERE, WHICH IS THE POINT OF REPORTING `spread`.
    ## A first version of this gate asserted `sigma == 1/r2` to 1e-9 and
    ## failed at 1.01e-04 -- the plateau is FLAT to 1e-4 and WRONG by 1e-4,
    ## because `C` still carries the 1e-25 F singular value that the rank
    ## tolerance swallowed.  Flatness certifies coherence, not convergence.
    ref = _explicit_block_sigma(ladder(1e16))
    assert abs(sigma - ref) <= (info['spread'] - 1.0) * ref, (sigma, ref)
    ## ⚠ AND THE LADDER IS U-SHAPED HERE, WHICH IS WHY THE READOUT IS A
    ## WINDOW AND NOT AN END.  The ratio descends THROUGH the true value and
    ## climbs again as the 1e-25 F singular value starts to tell, so the last
    ## point is the WORST: reading it gives 1.0e-04, the flattest window gives
    ## 2.0e-10.  Pin the gap so the readout cannot quietly move to an end.
    good = [r for r, u in zip(info['ratio'], info['usable']) if u]
    assert abs(good[-1] - ref) > 1e-5 * ref, (
        'the last ladder point is supposed to be BAD on this fixture; if it '
        'went good the gate no longer distinguishes the readouts', good[-1])


def test_algebraic_conditioning_decides_on_FLATNESS_not_on_descent():
    """⚠ A DESCENT RULE IS THE WRONG RULE AND THE ROUNDOFF FLOOR IS WHY.

    On a singular block the ratio falls like `h` while the leading term is
    still resolved -- on the C-V loop, four clean decades of exactly that --
    and then TURNS UPWARD below the floor, because what is left is roundoff
    divided by `h`.  Measured on random systems: 1.3e-07, 5.2e-09, 7.7e-07,
    3.0e-05 across four decades.  So "is it decreasing" reports a singular
    block as healthy at small enough `h`, and the floor guard is what keeps
    the turn out of the usable ladder.

    This gate asserts both halves: the usable ladder on a singular block
    descends by decades, and the points the floor guard REJECTED are the ones
    that would have broken a descent rule.
    """
    _rc, cv, _ode = _ac_fixtures()
    sigma, info = algebraic_conditioning(cv())
    usable = [r for r, u in zip(info['ratio'], info['usable']) if u]
    assert len(usable) >= 3, info
    ## a decade of `h` is a decade of ratio -- the signature of a zero block
    for a, b in zip(usable, usable[1:]):
        assert 8.0 <= a / b <= 12.0, usable
    ## and the spread is what condemns it, not the direction.  `spread` is
    ## now the FLATTEST 3-POINT WINDOW, so on a block falling a decade per
    ## probe its best window still spans two decades: the arithmetic says 100
    ## and it measures 99.0.  Bound at 50 -- naming the predicted value first,
    ## because a gate written as `> 100` fails on 99.0 and invites widening.
    assert info['spread'] > 50.0, info['spread']
    ## the guard actually rejected points -- otherwise this proves nothing
    assert not all(info['usable']), info['usable']


def _p1_vcvs_in_a_cv_loop(g, c1=1e-9, c2=2e-9):
    """The documented P1 fixture: a VCVS of gain `g` inside a C-V loop.
    Index 2 off `g* = 1 + c2/c1` and index 3 on it, and it DC-solves cleanly
    and silently at every gain.  `topological_index` cannot classify the VCVS."""
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('v')
    c.add_node('b')
    c['i'] = IS(gnd, 'v', i=1e-3)
    c['r'] = R('v', gnd, r=10.0)
    c['c1'] = C('v', 'b', c=c1)
    c['c2'] = C('v', gnd, c=c2)
    c['e'] = VCVS('v', gnd, 'b', gnd, g=g)
    return c


def test_the_index_0_rung_refuses_to_fire_on_an_unclassified_element():
    """⚠⚠ A NEGATIVE CLAIM CANNOT BE MADE PROVISIONALLY.

    Theorem 3.47's index-0 test asserts an ABSENCE -- a capacitive path from
    every node to datum AND NO VOLTAGE SOURCES.  An element outside the
    classifier's covered class is `'?'`, so `kinds[nm] == 'V'` is FALSE for it
    and the absence test passes VACUOUSLY.  A VCVS is exactly that: a voltage
    source the classifier does not recognise.

    The rung shipped WITH this defect on 2026-09-11 and returned INDEX 0 for
    the P1 fixture at every gain -- a circuit that is index 2 or 3.  The
    docstring of `_resolve_x0_unknown` records that every provisional fixture
    returned 1, which is what this restores.

    The other criteria need no such guard because they assert PRESENCE: "I
    found a C-V loop" stands whatever else is in the netlist, and `provisional`
    then says only that there may be more.  An absence cannot be established
    from a partial reading at all.
    """
    for g in (0.0, 2.999, 3.0, 4.0):
        idx, info = topological_index(_p1_vcvs_in_a_cv_loop(g))
        assert info['provisional'], g
        assert idx == 1, ('index-0 rung fired on an unclassified voltage '
                          'source', g, idx)

    ## ⚠ AND THE FIX MUST NOT BE "NEVER RETURN 0" -- the rung still fires
    ## where it should, on a netlist the classifier reads completely.
    tank = SubCircuit()
    tank.add_node('v')
    tank['C'] = C('v', gnd, c=1.0)
    tank['L'] = L('v', gnd, L=1.0)
    idx, info = topological_index(tank)
    assert idx == 0 and not info['unclassified'], (idx, info['unclassified'])


def test_the_reference_index_sees_a_ONE_DIMENSIONAL_near_singular_block():
    """⚠⚠ `_measured_index` was BLIND whenever `dim(ker C) == 1`.

    It tested `s2.min()/s2.max() < 1e-10`, and for a ONE-dimensional null
    space that ratio is IDENTICALLY 1 -- so it could only ever fire on an
    EXACTLY zero block, and returned 1 for a block of 1e-20.  No fixture here
    had a 1-D near-singular block, so nothing failed and nothing caught it.

    The fix needs an EXTERNAL scale: `s2.min() <= tol * ||G||`, dimensionless
    and free of the capacitance unit, since an absolute tolerance would give a
    false index-2 window widening as `1/||C||`.
    """
    ## the null space really is one-dimensional -- that is what made it blind
    cir = _e5_circuit(2e-14)
    Cm = np.asarray(cir.C(np.zeros(cir.n), defaultepar), dtype=float)
    from pycircuit.circuit.analysis import remove_row_col
    import pycircuit.circuit.analysis as _an
    Cm, = remove_row_col((Cm,), cir.get_node_index(gnd), _an.numeric)
    Cm = np.asarray(Cm, dtype=float)
    sv = np.linalg.svd(Cm, compute_uv=False)
    tol = max(Cm.shape) * np.finfo(float).eps * sv[0]
    assert int(np.sum(sv <= tol)) == 1, sv

    assert _measured_index(_e5_circuit(2e-14)) == 2
    assert _measured_index(_e5_circuit(0.0)) == 2
    ## and a healthy block still reads 1 -- the fix is not "always 2"
    assert _measured_index(_e5_circuit(2e-3)) == 1


def test_pss_warns_when_the_numeric_block_disagrees_with_the_topological_index():
    """The consumer: say so when the two instruments disagree and the netlist
    would otherwise solve cleanly and report nothing.

    ⚠⚠ THE CONTROL IS THE WHOLE TEST.  A warning that fires whenever the
    reading is PROVISIONAL would be noise -- most provisional netlists are
    perfectly ordinary.  It must fire on a provisional netlist whose block is
    SINGULAR and stay silent on a provisional netlist whose block is not.

    ⚠ A first version of that control put a VCVS across a capacitor and
    expected silence.  That is a C-V LOOP: the fixture was genuinely index 2
    and the warning was right.  A VCCS is the correct control -- it adds no
    branch equation, so it makes the reading provisional without changing the
    index.
    """
    def _warned(cir):
        p = PSS(cir)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            resolved = p._resolve_x0_unknown(None)
            hits = [str(w.message) for w in caught
                    if 'ALGEBRAIC BLOCK' in str(w.message)]
        ## never a behaviour change -- the remedy is not known to apply at
        ## index 3, so this reports and does not act
        assert resolved is False, resolved
        return hits

    for g in (2.999, 3.0):
        hits = _warned(_p1_vcvs_in_a_cv_loop(g))
        assert len(hits) == 1, (g, hits)
        assert 'index >= 2' in hits[0] and 'x0_unknown=True' in hits[0]

    per = 1e-3

    def base():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VSin('a', gnd, va=0.8, freq=1.0 / per)
        c['rs'] = R('a', 'b', r=50.0)
        c['cl'] = C('b', gnd, c=1e-9)
        c['rl'] = R('b', gnd, r=1e4)
        return c

    ## a plain, fully classified, well-conditioned netlist
    assert _warned(base()) == []

    ## ⚠ THE CONTROL THAT MATTERS: PROVISIONAL, and still silent
    c = base()
    c.add_node('d')
    c['g'] = VCCS('b', gnd, 'd', gnd, gm=1e-3)
    c['rd'] = R('d', gnd, r=1e3)
    c['cd'] = C('d', gnd, c=1e-9)
    assert topological_index(c)[1]['provisional'], 'control is not provisional'
    assert _warned(c) == []


def test_algebraic_conditioning_is_history_free_above_the_critical_voltage():
    """The re-sync `limit(x, x)` clamps against the STORED state, so from a
    stale one it lands short of `x` wherever the junction is above its
    critical voltage (about 0.73 V here) -- and the gate above reads at
    `x = 0`, where nothing limits, so it cannot see this.  Measured before
    the fix, a diode on an ALGEBRAIC node at 0.8 V: sigma 0.990 on a fresh
    circuit, 0.728 from a stored 5 V, 0.0210 from a stored 0 V.  On
    `_diode_fixture` the diode's node is differential and sigma does not
    move, which is why this needs its own circuit.
    """
    circuit.default_toolkit = circuit.numeric

    def build(stale=None):
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c.add_node('o')
        c['vs'] = VSin('a', gnd, va=0.8, freq=1e3, vo=0.6)
        c['rs'] = R('a', 'b', r=50.0)
        c['d'] = Diode('b', gnd)          # no capacitor on `b`: algebraic
        c['rl'] = R('b', 'o', r=1e3)
        c['cl'] = C('o', gnd, c=1e-9)
        if stale is not None:
            c['d'].__dict__['_vlim'] = stale
        return c

    c0 = build()
    x = np.zeros(c0.n)
    x[c0.get_node_index('a')] = 0.81
    x[c0.get_node_index('b')] = 0.8
    fresh, info = algebraic_conditioning(c0, x=x)
    assert info['verdict'] == 'well-conditioned' and abs(fresh - 0.990) < 1e-2, fresh
    for stale in (0.0, 5.0, -5.0):
        got, _info = algebraic_conditioning(build(stale), x=x)
        assert abs(got / fresh - 1.0) < 1e-12, (stale, got, fresh)
    ## and the state it found is the state it leaves
    c = build(0.0)
    algebraic_conditioning(c, x=x)
    assert c['d'].__dict__['_vlim'] == 0.0


def test_algebraic_conditioning_leaves_the_limiting_state_as_it_found_it():
    """⚠⚠ "A DIAGNOSTIC THAT CHANGES THE SIMULATION IS A DEFECT, AND THIS ONE
    DID" -- recorded on `transient.py`'s branch check (`branch_check`),
    which left `_vlim` at a speculative solve's value and moved the NEXT step's
    Jacobian.  Making this routine's reading PURE required re-syncing the
    limiting state with `limit(x, x)`, so it now has the same obligation.
    """
    cir = _diode_fixture()
    cir['d'].__dict__['_vlim'] = 5.0
    algebraic_conditioning(cir)
    assert cir['d'].__dict__.get('_vlim') == 5.0, \
        cir['d'].__dict__.get('_vlim')

    ## and it must not CREATE state on a circuit that had none
    fresh = _diode_fixture()
    assert '_vlim' not in fresh['d'].__dict__
    algebraic_conditioning(fresh)
    assert '_vlim' not in fresh['d'].__dict__, fresh['d'].__dict__.get('_vlim')


def test_the_orbit_rate_is_the_daes_own_and_its_stencil_fallback_is_live_and_second_order():
    """`PAC._orbit_rate` -- the rate that turns a node's motion in time into a
    state change for every fixed-time consumer -- and its fallback, given a
    test that fires it (refactor E9 item 5, 2026-09-23).

    On the driven noisy RC the source node's rate is ``w cos(w t)`` exactly,
    and the DAE form reads it to the last digit at every node.  ⚠ Including
    NODE 0, which it read as 0 until today: `VSin` clamps `t - td` at zero
    (SPICE's rule for a transient), so its derivative at exactly t = td = 0
    is the LEFT one, where the periodic steady state (t = 0 == t = T) has
    `w`.  Node 0 is now evaluated as node N.  It never moved a consumer --
    the fixed-time correction multiplies the rate by the node's own motion,
    zero at node 0 -- which is why it survived.

    The stencil fallback (the warned path for a singular differential /
    algebraic split, an index above one) agrees with the DAE rate to its own
    second order: 6.7e-4 / 1.66e-4 / 4.1e-5 at 100 / 200 / 400 points,
    ratios 4.04 and 4.02.  Forcing the split singular fires it, warns, and
    returns exactly the stencil."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    w = 2.0 * np.pi * 1e3
    errs = []
    for N in (100, 200, 400):
        cir = _rc_noisy()
        p = PSS(cir, method='radau', reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=1e-3, timestep=1e-3 / N, maxiterations=40)
        pac = PAC(cir, toolkit=circuit.numeric)
        rd = pac._orbit_rate(p, [])
        rs = pac._orbit_rate_stencil(p, [])
        ts = np.asarray(p.waveform[0], dtype=float)
        ia = [str(n_) for i, n_ in enumerate(cir.nodes) if i != p.irefnode].index('a')
        assert np.max(np.abs(rd[:, ia] - w * np.cos(w * ts))) < 1e-9 * w
        assert abs(rd[0, ia] - w) < 1e-9 * w, rd[0, ia]           # node 0 is node N
        errs.append(float(np.max(np.abs(rs - rd)) / np.max(np.abs(rd))))
    assert errs[0] < 1e-3 and 3.5 < errs[0] / errs[1] < 4.5 and 3.5 < errs[1] / errs[2] < 4.5, errs
    ## the fallback fires on a singular split, warns, and IS the stencil
    zero = lambda x: np.zeros((len(x), len(x)))
    p._C_at, p._G_at = zero, zero
    with pytest.warns(RuntimeWarning, match='three-node stencil'):
        rf = pac._orbit_rate(p, [])
    assert np.array_equal(rf, pac._orbit_rate_stencil(p, []))
