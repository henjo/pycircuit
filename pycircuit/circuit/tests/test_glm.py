"""The Nordsieck general linear method: tableau checks and the index-2 no-split gate.

Voigtmann (PhD, HU Berlin 2006) Theorem 9.5: a stiffly accurate GLM with stage order q = p, V power
bounded and M_inf nilpotent converges at order p on an index-2 DAE at constant stepsize -- in every
component.  Radau IIA(3) reads 5 / 3 and ESDIRK43 4 / 2 on the C-V loop fixture here
(test_radau_keeps_order_on_differential_and_loses_two_on_algebraic_index2); a q = p method must
read p / p.  Hypothesis (c) of the theorem -- the input vector exact to O(h^p) -- is the starting
machinery this file gates: the method is run from the EXACT Nordsieck vector and from the computed
one, and the orders must agree.
"""
import warnings
import numpy as np
import pytest

from pycircuit.circuit import circuit
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import C, L, R, VSin
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import (GLM2Integrator, GLM3Integrator,
                                          GLM4Integrator, RadauIIA3Integrator)


@pytest.mark.parametrize('cls', [GLM2Integrator, GLM3Integrator, GLM4Integrator])
def test_nordsieck_glm_tableau_is_order_p_stiffly_accurate_and_stable(cls):
    """Order and stage order are STRUCTURAL (U, V from the order conditions):
    one step on y = t^k is exact for k <= p in output and stages, and not at
    k = p + 1 (or the check saw nothing).  Stability is what a tableau has
    to be checked for: eig(V) = {1, ~0, ~0}, rho(M_inf) ~ 0 (nilpotent to
    the numerical construction's residual), rho(M(z)) <= 1 on the left
    half-plane; and stiff accuracy (B[0] == A[-1], c_s == 1, V[0] == U[-1])."""
    g = cls()
    v = g.verify()
    p = g.order
    for k, e_out, e_stage in v['exactness']:
        if k <= p:
            assert e_out < 1e-12 and e_stage < 1e-12, (k, e_out, e_stage)
        else:
            assert e_out > 1e-3 and e_stage > 1e-3, (k, e_out, e_stage)
    assert v['stiffly_accurate']
    assert abs(v['eig_V'][0] - 1.0) < 1e-6 and np.all(v['eig_V'][1:] < 5e-3), v['eig_V']
    ## the pair/triple at eps^(1/k) is a defective block's computed spectrum,
    ## not a residual to tighten -- the assertion is `Vdot^p = 0` above
    assert v['rho_M_inf'] < 5e-3, v['rho_M_inf']
    ## ⚠ RELATIVE to the tableau's own scale: GLM4's B carries entries of
    ## order 2500, so an absolute floor here would be meaningless (the docs
    ## session's point -- a residual that is too small is as uninformative as
    ## one that is too large).
    scale = max(float(np.max(np.abs(np.asarray(g.B, dtype=float)))), 1.0)
    assert v['nilpotency_residual'] < 1e-8 * scale, (v['nilpotency_residual'], scale)
    assert v['rho_lhp'] <= 1.0 + 1e-9, v['rho_lhp']


def _cv_loop(per):
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
    c['c1'] = C('a', 'b', c=1e-9)
    c['c2'] = C('b', gnd, c=1e-9)
    c['r'] = R('b', gnd, r=1e5)
    return c


def _reference(cir, per):
    """The analytic periodic solution and the subspaces, as the radau
    index-2 test builds them (validated in the time domain)."""
    n = cir.n
    iref = cir.get_node_index(gnd)
    keep = [i for i in range(n) if i != iref]
    Cfull = np.asarray(cir.C(np.zeros(n)), dtype=float)
    Gfull = np.asarray(cir.G(np.zeros(n)), dtype=float)
    Cm = Cfull[np.ix_(keep, keep)]
    Gm = Gfull[np.ix_(keep, keep)]
    m = len(keep)
    w = 2 * np.pi / per
    _U, sv, Vt = np.linalg.svd(Cm)
    d = int(np.sum(sv > m * sv[0] * np.finfo(float).eps))
    N = Vt[d:].T
    s2 = np.linalg.svd(N.T @ Gm @ N, compute_uv=False)
    assert float(s2[-1] / max(s2[0], 1e-300)) < 1e-10, 'not index 2'
    Vdiff, Valg = Vt[:d].T, Vt[d:].T
    NS = 2048
    ts = np.arange(NS) * per / NS
    Us = np.array([np.asarray(cir.u(t, analysis='tran'), dtype=float)[keep] for t in ts])
    Uph = (2.0 / NS) * np.sum(Us * np.exp(-1j * w * ts)[:, None], axis=0)
    X = np.linalg.solve(Gm + 1j * w * Cm, -Uph)
    chk = 0.0
    for t in (0.0, per * 0.137, per * 0.41, per * 0.76):
        xd = np.real(1j * w * X * np.exp(1j * w * t))
        ut = np.asarray(cir.u(t, analysis='tran'), dtype=float)[keep]
        chk = max(chk, float(np.max(np.abs(Cm @ xd + Gm @ np.real(X * np.exp(1j * w * t)) + ut))))
    assert chk < 1e-12, chk
    Xfull = np.zeros(n, dtype=complex)
    Xfull[keep] = X

    def x_exact(t):
        return np.real(Xfull * np.exp(1j * w * t))

    def q_nordsieck_exact(tn, h, p):
        ## Q_k = h^k d^k q/dt^k, q = C x (linear capacitors), full-size
        return np.array([h ** k * (Cfull @ np.real((1j * w) ** k * Xfull * np.exp(1j * w * tn)))
                         for k in range(p + 1)])
    return keep, Vdiff, Valg, x_exact, q_nordsieck_exact


def _orders(cir, per, integ, keep, Vdiff, Valg, x_exact, override=None, npts_list=(10, 20, 40, 80)):
    errs = []
    for npts in npts_list:
        tran = Transient(cir, toolkit=circuit.numeric, integrator=integ, reltol=1e-12,
                         iabstol=1e-16, vabstol=1e-14)
        if override is not None:
            tran._glm_startup_override = override
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(refnode=gnd, tend=per, x0=x_exact(0.0), timestep=per / npts,
                             fixed_timestep=True)
        t = np.asarray(res.sweep_values, dtype=float)
        X = np.asarray(res.x, dtype=float)
        X = X if X.shape[0] == cir.n else X.T
        E = np.array([X[keep, k] - x_exact(t[k])[keep] for k in range(len(t))])
        errs.append((float(np.max(np.abs(E @ Vdiff))), float(np.max(np.abs(E @ Valg)))))
    od = [np.log2(errs[i][0] / errs[i + 1][0]) for i in range(len(errs) - 1)]
    oa = [np.log2(errs[i][1] / errs[i + 1][1]) for i in range(len(errs) - 1)]
    return errs, od, oa


@pytest.mark.parametrize('cls', [GLM2Integrator, GLM3Integrator, GLM4Integrator])
def test_nordsieck_glm_has_no_index2_order_split_and_the_starting_vector_does_not_limit_it(cls):
    """On the index-2 C-V loop, from the exact periodic state, one period at
    constant step: the GLM's DIFFERENTIAL and ALGEBRAIC errors both fall at
    order p (no split), first with the EXACT Nordsieck starting vector fed
    through `_glm_startup_override` (the method alone), then with the
    computed one (`Transient._glm_startup`: p Radau IIA(3) substeps and a
    degree-p interpolant) -- the same orders, so the starting machinery
    delivers hypothesis (c).  Radau on the same harness is printed as the
    comparison (its algebraic component is the reduced one)."""
    per = 1e-3
    cir = _cv_loop(per)
    keep, Vdiff, Valg, x_exact, q_exact = _reference(cir, per)
    g = cls()
    p = g.order
    exact = lambda tn, x0, h: q_exact(tn, h, p)
    e_ex, od_ex, oa_ex = _orders(cir, per, g, keep, Vdiff, Valg, x_exact, override=exact)
    e_cp, od_cp, oa_cp = _orders(cir, per, cls(), keep, Vdiff, Valg, x_exact)
    assert e_ex[0][0] > 1e-9 and e_ex[0][1] > 1e-12, ('coarsest point near the floor', e_ex)
    for od, oa, label in ((od_ex, oa_ex, 'exact start'), (od_cp, oa_cp, 'computed start')):
        assert od[-1] > p - 0.4, (label, od)
        assert oa[-1] > p - 0.4, (label, oa)
        ## ⚠ THE SPLIT THAT MATTERS IS ONE-SIDED: the ALGEBRAIC component
        ## falling below p is the index-2 order reduction (radau: 5 / 3).
        ## The differential component coming out ABOVE p is superconvergence
        ## -- GLM4's c-bounded tableau reads 4.85 / 4.00 -- and a two-sided
        ## |od - oa| < 0.5 would fail on it, which is the wrong verdict.
        assert oa[-1] > od[-1] - 1.0, ('the ALGEBRAIC component must not lag',
                                       label, od, oa)
    ## the computed start must not cost order or accuracy against the exact one
    assert abs(od_cp[-1] - od_ex[-1]) < 0.5 and abs(oa_cp[-1] - oa_ex[-1]) < 0.5, (od_cp, od_ex, oa_cp, oa_ex)
    assert e_cp[-1][0] < 3.0 * e_ex[-1][0] and e_cp[-1][1] < 3.0 * e_ex[-1][1], (e_cp[-1], e_ex[-1])


def _pss_errors(method, per, keep, Vdiff, Valg, x_exact, npts_list=(20, 40, 80)):
    from pycircuit.circuit.shooting import PSS
    errs = []
    for npts in npts_list:
        p = PSS(_cv_loop(per), method=method, reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p.solve(period=per, timestep=per / npts, maxiterations=40)
        assert p.converged, (method, npts)
        t = np.asarray(p.waveform[0], dtype=float)
        X = np.asarray(p.waveform[1], dtype=float)[keep, :]
        E = np.array([X[:, k] - x_exact(t[k])[keep] for k in range(len(t))])
        errs.append((float(np.max(np.abs(E @ Vdiff))), float(np.max(np.abs(E @ Valg)))))
    return errs


def test_shooting_with_a_glm_keeps_the_algebraic_order_at_one_factorisation_per_step():
    """`PSS(method='glm3')` on the index-2 C-V loop: the shooting solve
    converges and the ALGEBRAIC component holds order 3 -- Radau IIA(3)'s
    order on this fixture -- with one real factorisation per step, where
    esdirk43 (also one per step, stage order 2) drops to 2.

    ⚠ The Newton's Jacobian here is APPROXIMATE by construction
    (`_traverse_glm`: the startup's derivative is dropped) and the residual
    is not, so the fixed point is the method's own; this test is what says
    the approximation does not cost the solve.  Measured at npts = 80:
    glm3 8.2e-11 against esdirk43's 5.7e-10 and radau's 2.0e-11.
    """
    per = 1e-3
    cir = _cv_loop(per)
    keep, Vdiff, Valg, x_exact, _q = _reference(cir, per)
    e_glm = _pss_errors('glm3', per, keep, Vdiff, Valg, x_exact)
    e_esd = _pss_errors('esdirk43', per, keep, Vdiff, Valg, x_exact)
    oa_glm = np.log2(e_glm[-2][1] / e_glm[-1][1])
    oa_esd = np.log2(e_esd[-2][1] / e_esd[-1][1])
    assert oa_glm > 2.6, (oa_glm, e_glm)
    assert oa_esd < 2.4, (oa_esd, e_esd)
    assert e_glm[-1][1] < 0.5 * e_esd[-1][1], (e_glm[-1], e_esd[-1])


def test_the_glm_period_map_is_on_the_nordsieck_state_and_carries_the_circuits_multiplier():
    """A multivalue method's period map acts on the whole Nordsieck vector,
    so `factored_period()` comes back with width `r*m = (p+1)*m`, not `m`.
    Its spectrum must then be the circuit's `m` multipliers plus `(r-1)*m`
    of the METHOD's own, which sit at zero because `V`'s lower block is
    nilpotent -- the same shape as a DAE's structural zeros.

    Measured on the index-2 C-V loop at 40 points: the dominant multiplier
    is 6.737685e-03 against radau's 6.737685e-03 (six digits, two entirely
    different maps), and the other eleven are below 1e-18.
    """
    from pycircuit.circuit.shooting import PSS
    per = 1e-3
    m = _cv_loop(per).n - 1
    mults = {}
    for method in ('glm3', 'radau'):
        p = PSS(_cv_loop(per), method=method, reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p.solve(period=per, timestep=per / 40, maxiterations=40)
        fp = p.factored_period()
        M = np.column_stack([fp.matvec(e) for e in np.eye(fp.width)])
        mults[method] = (fp, np.sort(np.abs(np.linalg.eigvals(M)))[::-1])
    fp_glm, lam_glm = mults['glm3']
    _fp_r, lam_rad = mults['radau']
    assert fp_glm.kind == 'glm'
    assert fp_glm.width == 4 * m, (fp_glm.width, m)      # r = p + 1 = 4
    assert abs(lam_glm[0] / lam_rad[0] - 1.0) < 1e-4, (lam_glm[0], lam_rad[0])
    assert np.all(lam_glm[1:] < 1e-12), lam_glm
    assert np.all(lam_rad[1:] < 1e-12), lam_rad


def _vdp(Q=15.9):
    from pycircuit.circuit.elements import BSource
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2 * np.pi * Q)
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: mu * (u - u ** 3 / 3.0))
    return c, 2 * np.pi / np.sqrt(1 - mu ** 2 / 4)


def test_a_glm_finds_the_free_period_and_the_orbits_multipliers():
    """The autonomous path: `(x_0, T)` unknown, with the period column from
    `_traverse_glm(want_dT=True)`.

    ⚠ A MULTIVALUE METHOD HAS TWO EXPLICIT `T` DEPENDENCES, not one: the
    grid's `h = frac T` inside every step, AND the starting vector's own
    `Q_k = h^k q^(k)` scaling, `dQ_k/dT = (k/T) Q_k`.  Both are carried; what
    is dropped is the Radau substeps inside the startup, the same term the
    driven Jacobian drops.

    Measured on van der Pol at Q = 15.9 against radau's period
    (6.283224654, converged to the digit at both grids):

        method    60 points        120 points
        glm3      6.283226014      6.283224730
        trbdf2    6.286107414      6.283933050

    i.e. RELATIVE period errors 2.2e-07 / 1.2e-08 for glm3 (ratio 17.7, so
    order ~4.1 in the period) against TR-BDF2's 4.6e-04 / 1.1e-04 -- four
    orders below the same-cost method at the finer grid -- and its Floquet
    pair comes back at 1.000000 / 0.939050 against radau's 1.000000 /
    0.939043.  ⚠ The assertions are on RELATIVE errors; the table above is
    absolute periods.
    """
    from pycircuit.circuit.shooting import PSS
    _c, T0 = _vdp()
    per = {}
    lam_glm = lam_rad = None
    for method in ('radau', 'glm3', 'trbdf2'):
        for npts in (60, 120):
            c, _T = _vdp()
            p = PSS(c, method=method, reltol=1e-12)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                p.solve(period=T0, timestep=T0 / npts, x0=np.array([2.0, 0.0]),
                        maxiterations=200)
            assert p.converged, (method, npts)
            per[(method, npts)] = float(p.period)
            if npts == 120 and method in ('glm3', 'radau'):
                fp = p.factored_period()
                M = np.column_stack([fp.matvec(e) for e in np.eye(fp.width)])
                lam = np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
                if method == 'glm3':
                    lam_glm = lam
                else:
                    lam_rad = lam
    ref = per[('radau', 120)]
    e_glm = [abs(per[('glm3', n)] / ref - 1.0) for n in (60, 120)]
    e_trb = [abs(per[('trbdf2', n)] / ref - 1.0) for n in (60, 120)]
    assert e_glm[1] < 1e-7, e_glm
    assert e_glm[0] / e_glm[1] > 8.0, ('the period should converge, not sit', e_glm)
    assert e_glm[1] < 1e-3 * e_trb[1], (e_glm, e_trb)
    assert abs(lam_glm[0] - 1.0) < 1e-5, lam_glm[:3]
    assert abs(lam_glm[1] / lam_rad[1] - 1.0) < 1e-3, (lam_glm[1], lam_rad[1])
    assert np.all(lam_glm[2:] < 1e-6), lam_glm


def test_the_glm_adjoint_is_the_transpose_of_its_period_map():
    """`matvec_transposed` for `kind='glm'` -- the reverse-mode adjoint of the
    multivalue recursion, which is what a phase-sensitivity or small-signal
    surface would replay over a GLM operating point.

    Asserted against the dense transpose of the forward map (built by a
    different recursion, not by transposing this one), the adjoint identity
    on random vectors, the complex path, and the `collect` contract (one
    adjoint state per step, at the map's own `r*m` width).
    """
    from pycircuit.circuit.shooting import PSS
    per = 1e-3
    p = PSS(_cv_loop(per), method='glm3', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p.solve(period=per, timestep=per / 40, maxiterations=40)
    fp = p.factored_period()
    w = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(w)])
    MT = np.column_stack([fp.matvec_transposed(e) for e in np.eye(w)])
    scale = max(float(np.max(np.abs(M))), 1e-300)
    assert np.max(np.abs(MT - M.T)) < 1e-12 * scale, np.max(np.abs(MT - M.T))
    rng = np.random.default_rng(0)
    a, b = rng.standard_normal(w), rng.standard_normal(w)
    assert abs(float(a @ (MT @ b)) - float((M @ a) @ b)) < 1e-12 * scale
    z = a + 1j * b
    assert np.max(np.abs(fp.matvec_transposed(z) - MT @ z)) < 1e-12 * scale
    _end, ts, states = fp.matvec_transposed(np.eye(w)[0], collect=True)
    assert len(states) == len(fp.steps) and states[0].shape == (w,)
    assert len(ts[0]) == 4                       # one reverse solve per stage


def test_the_order_four_glm_beats_radau_on_the_algebraic_component_at_one_lu_per_step():
    """⚠⚠ THE POINT OF THE WHOLE GLM ARC, measured.

    On an index-2 DAE Radau IIA(3) splits: classical order 5 in the
    differential components, 3 in the ALGEBRAIC ones (`min(p, q)` with
    `q = 3` the stage order).  A `q = p` method has no split -- so an
    order-4 one should carry order 4 in BOTH, and beat Radau exactly where
    the reduction bites, at one real factorisation per step against Radau's
    coupled `3n` solve.

    Driven PSS on the C-V loop, algebraic component, npts 20 / 40 / 80::

        glm4   3.9e-11  2.2e-12  1.3e-13   order 4.07
        radau  1.4e-09  1.6e-10  2.0e-11   order 3.05

    150x at the finest grid.  ⚠ Radau still wins the DIFFERENTIAL component
    (order 5 against 4) and its tableau is far better scaled; this test
    asserts the algebraic claim only, which is the one the theorem makes.
    """
    per = 1e-3
    cir = _cv_loop(per)
    keep, Vdiff, Valg, x_exact, _q = _reference(cir, per)
    e_glm = _pss_errors('glm4', per, keep, Vdiff, Valg, x_exact)
    e_rad = _pss_errors('radau', per, keep, Vdiff, Valg, x_exact)
    oa_glm = np.log2(e_glm[-2][1] / e_glm[-1][1])
    oa_rad = np.log2(e_rad[-2][1] / e_rad[-1][1])
    assert oa_glm > 3.6, (oa_glm, e_glm)
    assert oa_rad < 3.4, (oa_rad, e_rad)
    assert e_glm[-1][1] < 0.2 * e_rad[-1][1], (e_glm[-1][1], e_rad[-1][1])


def test_small_signal_surfaces_work_over_a_glm_operating_point_through_the_twin():
    """⚠⚠ `ppv` USED TO RETURN A WRONG ANSWER SILENTLY over a GLM orbit, and
    that is what this gate exists to stop.

    A multivalue method's own period map acts on the NORDSIECK state (width
    `r*m`), while `ppv`, `oscillator_covariance` and the PAC family all want
    a map on `x`.  The two are NOT related by taking the first block: a
    state kick `dx` perturbs the higher Nordsieck components too
    (`dQ_k = h^k d^k(C dx)/dt^k`).  Measured on van der Pol at Q = 15.9,
    both naive extractions -- `w[:m]` and `C^T w_0` -- are wrong by 2 % in
    norm and by 13x in the small component, and `ppv` returned the first
    of them with no complaint because it renormalised to `v . xdot = 1`,
    which any scaling satisfies.

    The fix is the machinery `trap` and `euler` already use for the same
    shape of reason: `carries_own_monodromy()` is False for a multivalue
    method, so the state-space surfaces take a TR-BDF2 TWIN re-solved on
    this orbit's own grid.  The orbit stays the GLM's; only the map is
    borrowed.  Measured here against radau on the same fixture: PPV 4e-4,
    the phase diffusion constant 1.3e-3, the oscillator spectrum 5e-4.
    """
    from pycircuit.circuit.shooting import PSS, PAC
    from pycircuit.circuit.elements import IS
    got = {}
    for method in ('radau', 'glm3'):
        cir, T0 = _vdp()
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        p = PSS(cir, method=method, reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p.solve(period=T0, timestep=T0 / 60, x0=np.array([2.0, 0.0]),
                    maxiterations=200)
        assert p.converged
        m = cir.n - 1
        pac = PAC(cir, toolkit=circuit.numeric)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            v, _info = p.ppv()
            _K, d, _i2 = pac.oscillator_covariance(p)
            Sv, _ = pac.oscillator_spectrum(p, np.array([0.05 / float(p.period)]), 0)
        got[method] = (np.asarray(v)[:m], d / float(p.period), float(Sv[0]),
                       p.monodromy_twin())
    (v_r, c_r, s_r, tw_r), (v_g, c_g, s_g, tw_g) = got['radau'], got['glm3']
    ## radau's own map is on `x`, so it is its own twin; the GLM borrows one
    assert tw_r is not None and tw_g is not None
    assert getattr(tw_g.par, 'method', None) == 'trbdf2', tw_g
    assert np.max(np.abs(v_g - v_r)) < 3e-3 * np.max(np.abs(v_r)), (v_g, v_r)
    assert abs(c_g / c_r - 1.0) < 5e-3, (c_g, c_r)
    assert abs(s_g / s_r - 1.0) < 5e-3, (s_g, s_r)


def test_floquet_modes_off_the_nordsieck_map_drops_the_methods_own_multipliers():
    """`floquet_modes` on a multivalue map: its null filter
    (`FLOQUET_NULL_TOL`) drops the `(r-1)*m` multipliers that are the
    METHOD's own -- they sit at zero because `V`'s lower block is nilpotent
    -- so what comes back is the CIRCUIT's pair, matching radau's to 1e-4
    on van der Pol.  ⚠ The MULTIPLIERS transfer; the EIGENVECTORS do not,
    they live in the `r*m` Nordsieck space, which is why the state-space
    surfaces take a twin instead (see the test above).
    """
    from pycircuit.circuit.shooting import PSS
    lam = {}
    for method in ('radau', 'glm3'):
        cir, T0 = _vdp()
        p = PSS(cir, method=method, reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p.solve(period=T0, timestep=T0 / 60, x0=np.array([2.0, 0.0]),
                    maxiterations=200)
        fm = p.floquet_modes(fp=p.factored_period())
        lam[method] = sorted((abs(x['lam']) for x in fm), reverse=True)
    assert len(lam['glm3']) == len(lam['radau']) == 2, lam
    assert abs(lam['glm3'][0] - 1.0) < 1e-5, lam['glm3']
    assert abs(lam['glm3'][1] / lam['radau'][1] - 1.0) < 1e-3, lam


@pytest.mark.parametrize('method', ['trap', 'gear', 'radau', 'glm3'])
def test_floquet_modes_works_when_called_with_no_arguments(method):
    """⚠ `PSS.floquet_modes()` -- its own documented default call -- used to
    raise `AttributeError: 'NoneType' object has no attribute
    'factored_period'` on EVERY method, because the body dereferenced the
    first parameter, which is named `pss_unused` and defaults to `None`.
    Every call site inside the package passes `self`, so nothing was
    failing in the suite and the defect only showed when a caller believed
    the signature.  Reading `self` instead changes no existing answer:
    asserted here by running both forms and requiring bit-identical
    multipliers, across a twin-taking method (trap, glm3), a method that is
    its own twin (gear, radau), and both map kinds.
    """
    from pycircuit.circuit.shooting import PSS
    cir, T0 = _vdp()
    p = PSS(cir, method=method, reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p.solve(period=T0, timestep=T0 / 60, x0=np.array([2.0, 0.0]),
                maxiterations=200)
    assert p.converged
    with_arg = p.floquet_modes(p)
    no_arg = p.floquet_modes()
    assert len(no_arg) == len(with_arg) and len(no_arg) >= 2
    for a, b in zip(with_arg, no_arg):
        assert abs(a['lam'] - b['lam']) < 1e-14, (method, a['lam'], b['lam'])


## ------------------------------------------------------------------------
## The stage predictor

def _expg_fixture(per, va=0.8):
    """A driven RC with a STATE-FREE exponential conductance.

    ⚠ `Diode` cannot measure a stage predictor at all: its `G` linearises
    around a STORED `_vlim`, so its Newton is seed-blind by construction --
    measured, an ALL-ZEROS seed gives the same iteration histogram, to the
    count, as the exact one.  This element's `i` and `G` are both functions
    of the passed `x`, with no state and no limiting, so the iteration count
    responds to the seed.
    """
    import sympy
    import pycircuit.circuit.circuit as _cc
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Parameter)

    class ExpG(Behavioural):
        instparams = [Parameter(name='IS', desc='sat', unit='A',
                                default=1e-12),
                      Parameter(name='VT', desc='thermal', unit='V',
                                default=0.026)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return (Contribution(b.I, IS * (sympy.exp(b.V / VT) - 1)),)  # noqa: F821

    _cc.default_toolkit = numeric
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=va, freq=1.0 / per)
    c['rs'] = R('a', 'b', r=50.0)
    c['nl'] = ExpG('b', gnd, IS=1e-12, VT=0.026)
    c['cl'] = C('b', gnd, c=1e-9)
    c['rl'] = R('b', gnd, r=1e4)
    return c


def _run_counted(cls, mode, per, npts, va=0.8):
    """One fixed-step transient; returns (device i evaluations, worst seed
    error, worst stage iteration count, waveform)."""
    from pycircuit.circuit import nrsolver as NR
    from pycircuit.circuit import transient as TR
    log = []
    solve_orig = NR.StandardNewton.solve_system

    def wrapped(self, x0, eval_FJ, *a, **k):
        x, it = solve_orig(self, x0, eval_FJ, *a, **k)
        log.append((float(np.max(np.abs(np.asarray(x0) - np.asarray(x)))), it))
        return x, it
    NR.StandardNewton.solve_system = wrapped
    TR.Transient.glm_predictor = mode
    n = {'i': 0}
    cir = _expg_fixture(per, va)
    fi = cir.i
    cir.i = lambda *a, **k: (n.__setitem__('i', n['i'] + 1), fi(*a, **k))[1]
    try:
        tr = Transient(cir, integrator=cls(), reltol=1e-9)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tr.solve(refnode=gnd, tend=per, timestep=per / npts,
                           fixed_timestep=True)
    finally:
        NR.StandardNewton.solve_system = solve_orig
        TR.Transient.glm_predictor = 'local'
    return (n['i'], float(np.max([l[0] for l in log])),
            int(np.max([l[1] for l in log])),
            np.asarray(res.v('b'), dtype=float))


@pytest.mark.parametrize('cls', [GLM3Integrator, GLM4Integrator])
def test_glm_stage_predictor_cuts_the_newton_work_without_moving_the_answer(cls):
    """The stage predictor is a COST change and must be nothing else.

    Three properties, each with the pre-predictor behaviour
    (`glm_predictor='none'`) as its own control in the same run:

    1. FEWER device evaluations.  Measured 13-29% on this fixture; gated at
       8% so a modest regression still fails.
    2. The SAME answer.  A seed change must not select a different root --
       the risk is real on a multi-rooted stage system (the coupled path's
       line-search note), so this is asserted to rounding, not to a tolerance.
    3. NO WORSE TAIL.  ⚠ This is the property that rejected the two simpler
       predictors: "previous stage shifted" and "the previous step's stages
       continued a whole step" both improved the MEAN and made the worst
       stage cost 13 and 33 Newton iterations against the fallback's 6.  A
       gate on the mean alone would have passed both.
    """
    per = 1e-3
    base_i, base_seed, base_it, base_x = _run_counted(cls, 'none', per, 200)
    pred_i, pred_seed, pred_it, pred_x = _run_counted(cls, 'local', per, 200)
    assert pred_i < 0.92 * base_i, (pred_i, base_i)
    assert np.max(np.abs(pred_x - base_x)) < 1e-11 * max(
        float(np.max(np.abs(base_x))), 1e-3), float(np.max(np.abs(pred_x - base_x)))
    assert pred_seed <= base_seed * 1.05, (pred_seed, base_seed)
    assert pred_it <= base_it, (pred_it, base_it)


@pytest.mark.parametrize('cls', [GLM2Integrator, GLM3Integrator, GLM4Integrator])
def test_glm_stage_predictor_is_exact_on_a_polynomial(cls):
    """STRUCTURAL, no circuit: the predictor interpolates the nodes it is
    given, so on a trajectory that IS a polynomial of its degree it returns
    the stage value exactly.  That is what says the weights are a Lagrange
    evaluation at `c_i` and not an approximation with a fitted constant."""
    g = cls()
    A, U, B, V, c, p = g.tableau()
    s = A.shape[0]
    deg = g.order
    W = Transient._glm_local_weights(c, deg)
    for k in range(deg + 1):
        ## y(tau) = tau^k in the step's local coordinate; the previous step's
        ## stages sit at c_j - 1 and this step's at c_j
        Yp = [np.array([(c[j] - 1.0) ** k]) for j in range(s)]
        Y = [np.array([c[j] ** k]) for j in range(s)]
        for i in range(s):
            ip, ic, w = W[i]
            got = sum(w[a] * Yp[j] for a, j in enumerate(ip)) \
                + sum(w[len(ip) + a] * Y[j] for a, j in enumerate(ic))
            assert abs(float(np.asarray(got).ravel()[0]) - c[i] ** k) < 1e-9, \
                (cls.__name__, k, i)


def test_glm_stage_predictor_declines_a_changed_step_and_a_period_seam():
    """The predictor reads the PREVIOUS step's stages, which are only at the
    right places if that step used the same `h` and ended where this one
    starts.  Both guards are asserted directly, because a stale record is a
    silently wrong guess, not an error."""
    tr = Transient(_expg_fixture(1e-3), integrator=GLM3Integrator())
    A, U, B, V, c, p = GLM3Integrator().tableau()
    s = A.shape[0]
    xn = np.zeros(4)
    tr.base_integrator = GLM3Integrator()
    ## no predecessor at all -- the state at a period seam
    tr._glm_prev = None
    assert tr._glm_stage_predictor(xn, 1e-5, 0.0, c) is None
    Yp = [np.zeros(4) for _ in range(s)]
    ## a good record predicts
    tr._glm_prev = (Yp, np.zeros(4), 1e-5, 0.0)
    assert tr._glm_stage_predictor(xn, 1e-5, 0.0, c) is not None
    ## a CHANGED STEP declines
    tr._glm_prev = (Yp, np.zeros(4), 2e-5, 0.0)
    assert tr._glm_stage_predictor(xn, 1e-5, 0.0, c) is None
    ## a record that does not END where this step STARTS declines
    tr._glm_prev = (Yp, np.zeros(4), 1e-5, 5e-5)
    assert tr._glm_stage_predictor(xn, 1e-5, 0.0, c) is None
