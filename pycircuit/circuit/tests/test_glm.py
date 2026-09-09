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
from pycircuit.circuit.elements import C, R, VSin
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import GLM2Integrator, RadauIIA3Integrator


def test_nordsieck_glm_tableau_is_order_p_stiffly_accurate_and_stable():
    """Order and stage order are STRUCTURAL (U, V from the order conditions):
    one step on y = t^k is exact for k <= p in output and stages, and not at
    k = p + 1 (or the check saw nothing).  Stability is what a tableau has
    to be checked for: eig(V) = {1, ~0, ~0}, rho(M_inf) ~ 0 (nilpotent to
    the numerical construction's residual), rho(M(z)) <= 1 on the left
    half-plane; and stiff accuracy (B[0] == A[-1], c_s == 1, V[0] == U[-1])."""
    g = GLM2Integrator()
    v = g.verify()
    p = g.order
    for k, e_out, e_stage in v['exactness']:
        if k <= p:
            assert e_out < 1e-12 and e_stage < 1e-12, (k, e_out, e_stage)
        else:
            assert e_out > 1e-3 and e_stage > 1e-3, (k, e_out, e_stage)
    assert v['stiffly_accurate']
    assert abs(v['eig_V'][0] - 1.0) < 1e-6 and np.all(v['eig_V'][1:] < 5e-3), v['eig_V']
    assert v['rho_M_inf'] < 5e-3, v['rho_M_inf']
    assert v['nilpotency_residual'] < 1e-4, v['nilpotency_residual']
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


def test_nordsieck_glm_has_no_index2_order_split_and_the_starting_vector_does_not_limit_it():
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
    g = GLM2Integrator()
    p = g.order
    exact = lambda tn, x0, h: q_exact(tn, h, p)
    e_ex, od_ex, oa_ex = _orders(cir, per, g, keep, Vdiff, Valg, x_exact, override=exact)
    e_cp, od_cp, oa_cp = _orders(cir, per, GLM2Integrator(), keep, Vdiff, Valg, x_exact)
    assert e_ex[0][0] > 1e-9 and e_ex[0][1] > 1e-11, ('coarsest point near the floor', e_ex)
    for od, oa, label in ((od_ex, oa_ex, 'exact start'), (od_cp, oa_cp, 'computed start')):
        assert od[-1] > p - 0.4, (label, od)
        assert oa[-1] > p - 0.4, (label, oa)
        assert abs(od[-1] - oa[-1]) < 0.5, ('the point is NO split', label, od, oa)
    ## the computed start must not cost order or accuracy against the exact one
    assert abs(od_cp[-1] - od_ex[-1]) < 0.3 and abs(oa_cp[-1] - oa_ex[-1]) < 0.3, (od_cp, od_ex, oa_cp, oa_ex)
    assert e_cp[-1][0] < 3.0 * e_ex[-1][0] and e_cp[-1][1] < 3.0 * e_ex[-1][1], (e_cp[-1], e_ex[-1])
