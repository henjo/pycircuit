"""STAGE 12C -- Fang's coupled (x, h) stepping around a STAGE method's (or a
Nordsieck GLM's) step: a PROTOTYPE, measured against the standard adaptive
stage run, which decided that `coupled_lte=True` stays refused for them
(2026-09-25; the numbers are at the refusal in `Transient._solve` and in
`doc/pss_log_260902.md`).

Fang's loop per time point (DAC 2013 Fig. 4, sec. 3.4 'approx'):
  solve the step at h (the method's own stage solve, converged);
  eq (6): the solution-space LTE, extrapolation of degree = the method's
  order (capped by the accepted history), normalised by the same etol;
  in the band [gamma_min, gamma_max]: accept; else h <- step_for_error_ratio
  (clamped to 1 -/+ eta) and re-solve (warm start: the step re-solved).
The accepted h carries forward as the next guess (no separate predictor).
Metrics: accepted steps, stage solves (each re-solve counts), wall clock,
error against the closed form (max / median, past the opening two points as
`stage12b_coupled` reads it).  ⚠ The `it` column is Newton iterations where
the transient counts them: the prototype's inner transient does not, and
radau's coupled step does not either -- a 0 there is not a measurement.
Run:
    PYTHONPATH=. python benchmarks/transient_review/stage12c_fang_stage_methods.py [methods...]
"""
import sys, time, warnings as _w
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pycircuit.circuit import circuit
circuit.default_toolkit = circuit.numeric
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.shooting import PSS
from pycircuit.circuit._lte_kernels import solution_lte, step_for_error_ratio
import stage12b_coupled as B

RELTOLS = (1e-4, 1e-5, 1e-6)
methods = sys.argv[1:] or ['radau', 'trbdf2']


def standard(builder, analytic, method, reltol):
    cir, kw, node = builder()
    integ = PSS(cir, method=method)._integrator_for(method)
    tr = Transient(cir, integrator=integ, reltol=reltol, toolkit=circuit.numeric)
    t0 = time.time()
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        res = tr.solve(**kw)
    wall = time.time() - t0
    w = res.v(node)
    t = np.asarray(w.x, dtype=float).ravel()
    v = np.asarray(w.y, dtype=float).ravel()
    st = tr.statistics
    err = np.abs(v - analytic(t))[2:]
    return dict(steps=st.accepted_steps, solves=st.accepted_steps + st.rejected_steps,
                rejected=st.rejected_steps, iters=st.newton_iterations, wall=wall,
                emax=float(err.max()), emed=float(np.median(err)))


def fang(builder, analytic, method, reltol, gamma=(0.7, 3.0), eta=0.15):
    cir, kw, node = builder()
    tend, h = float(kw['tend']), float(kw['timestep'])
    p = PSS(cir, method=method, reltol=reltol)
    tr = p._new_transient(p._integrator_for(method))
    p._tran = tr
    m = cir.n - 1
    iref = p.irefnode
    x0 = kw.get('x0')
    x = np.zeros(m) if x0 is None else np.delete(np.asarray(x0, dtype=float), iref)
    p._begin_period(x)
    integ = tr.base_integrator
    order = int(getattr(integ, 'ORDER', 2))
    ## the same tolerance the coupled LMM path uses (LTERATIO, reltol ref + lte_abstol)
    from pycircuit.circuit.stepcontroller import SolutionLTEController
    ctrl = SolutionLTEController().set_relref(tr.par.relref)
    t = 0.0
    xs_hist = [x.copy()]          # accepted states, newest first
    h_hist = []                   # accepted steps, newest first
    ts, vs = [0.0], [x.copy()]
    solves = 0
    it0 = tr.statistics.newton_iterations if hasattr(tr, 'statistics') else 0
    t0 = time.time()
    lo, hi = 1.0 - eta, 1.0 + eta
    target = (gamma[0] * gamma[1]) ** 0.5
    inode = [str(n_) for n_ in cir.nodes].index(node)
    io = inode if inode < iref else inode - 1
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        while t < tend * (1 - 1e-12):
            h = min(h, tend - t)
            h_entry = h
            for _it in range(20):
                xn = np.asarray(p.solve_timestep(x, t + h, h), dtype=float)
                solves += 1
                if len(xs_hist) < 2 or not h_hist:
                    break
                deg = min(order, len(xs_hist) - 1, len(h_hist))
                lte = solution_lte(xn, xs_hist[:deg + 1], h_hist[:deg], h)
                xf = np.insert(xn, iref, 0.0)
                xl = np.insert(xs_hist[0], iref, 0.0)
                ref = ctrl._reference(xf, xl, False, len(cir.nodes), tr.toolkit)
                etol = tr.LTERATIO * (tr.par.reltol * ref + tr._lte_abstol_vector())
                etol = np.delete(np.asarray(etol, dtype=float), iref)
                err = float(np.max(np.abs(lte) / etol))
                if gamma[0] <= err <= gamma[1]:
                    break
                if h >= tend - t and err < gamma[0]:
                    break           # landing on tend: cannot grow
                ratio = target / max(err, 1e-300)
                h_new = step_for_error_ratio(h, h_hist[:deg], ratio, lo, hi)
                h_new = min(max(h_new, h_entry * 0.1), h_entry * 10.0, tend - t)
                if abs(h_new - h) <= 1e-12 * h:
                    break
                h = h_new
            t = t + h
            x = xn
            xs_hist.insert(0, x.copy())
            h_hist.insert(0, h)
            del xs_hist[order + 2:]
            del h_hist[order + 1:]
            ts.append(t)
            vs.append(x.copy())
    wall = time.time() - t0
    ts = np.asarray(ts)
    v = np.asarray([vv[io] for vv in vs])
    err = np.abs(v - analytic(ts))[2:]
    it1 = tr.statistics.newton_iterations if hasattr(tr, 'statistics') else 0
    return dict(steps=len(ts) - 1, solves=solves, rejected=solves - (len(ts) - 1),
                iters=it1 - it0, wall=wall, emax=float(err.max()), emed=float(np.median(err)))


cases = [('rc-vsin', B.rc_vsin, B.rc_vsin_analytic), ('stiff-rlc', B.stiff_rlc, B.stiff_rlc_analytic)]
for method in methods:
    for name, bld, ana in cases:
        for rt in RELTOLS:
            s = standard(bld, ana, method, rt)
            try:
                f = fang(bld, ana, method, rt)
            except Exception as e:
                f = None
                print('%-6s %-9s rt %.0e  FANG EXC %s' % (method, name, rt, str(e)[:120]), flush=True)
            fmt = lambda d: 'steps %5d solves %5d rej %4d it %6d wall %6.2f emax %.2e emed %.2e' % (
                d['steps'], d['solves'], d['rejected'], d['iters'], d['wall'], d['emax'], d['emed'])
            print('%-6s %-9s rt %.0e  STD  %s' % (method, name, rt, fmt(s)), flush=True)
            if f is not None:
                print('%-6s %-9s rt %.0e  FANG %s' % (method, name, rt, fmt(f)), flush=True)
