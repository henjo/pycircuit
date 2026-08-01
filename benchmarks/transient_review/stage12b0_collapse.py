"""GATE 12B-0 -- reproduce the recorded step-size collapse before rebuilding on it.

`transient.py::_solve_coupled` carries this note:

    This is the paper's robust "approximate Newton" (sec. 3.4); it replaces the
    exact (N+1) Schur update, which is very sensitive to step changes and
    collapses the step size.

and commit 9e8fb74 says "Trying to enable the exact (N+1) Schur update
destabilised it".  Stage 12B is about to build on the exact update, so the rule
is: a fix reused from a previous failure needs its failure mode checked first --
two failures that look alike in a summary need not be the same failure.

READING THE CODE THAT PRODUCED THE NOTE (git 9e8fb74^) turns up two defects that
have nothing to do with the Schur update:

  1. THE SOLVED STEP WAS THROWN AWAY.  The loop does `h = h_next` after a
     converged coupled solve, then four lines later `t += h_curr; h = h_curr`.
     So the step size the (N+1) system solved for was overwritten by the one it
     started from, every single step.
  2. THE LTE WAS IDENTICALLY ZERO FOR THE DEFAULT METHOD.  `calc_E` branches on
     'trapezoidal'/'trap' and 'gear2' with `else: lte = zeros`; the default was
     'euler', so E was the constant -TRTOL and there was nothing to solve.

Both are recorded in 9e8fb74's own message.  On top of that the LTE formulas
there predate stage 4 -- the gear2 residual is `h^2 (h+h_last)/3 * dd2`, which
carries the 3/4 optimism stage 4i removed -- and the tolerance is charge
flavoured (`reltol*max(|q|,|q_last|) + abstol`) rather than solution flavoured.

So the note attributes to the exact update a collapse observed on a loop that
discarded its own answer, using an estimator later proved wrong.  This script
separates those: the old loop is reconstructed with each confound as a switch,
so the collapse can be attributed to something rather than assumed.

Run:
    MPLBACKEND=Agg PYTHONPATH=.:benchmarks python benchmarks/transient_review/stage12b0_collapse.py
"""
import warnings
from copy import copy

import numpy as np

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.nrsolver import SchurCoupledNewton, NoConvergenceError
from pycircuit.circuit.analysis import CircuitResult

import transient_stage4 as S

TRTOL = 7.0
## The band the original used, verbatim from git 9e8fb74^ -- and the paper's own
## ltemin/ltemax.  Inside it the step size is left alone: the method solves for a
## step CLOSE to the largest one meeting the LTE, not for an exact root.
BAND_MIN, BAND_MAX = 0.7, 3.0
MAX_STEPS = 20000
## Consecutive steps sitting on the floor before the run is called stalled.
STALL_STEPS = 50


class OldCoupled(Transient):
    """The pre-9e8fb74 coupled loop, with its confounds exposed as flags.

    Transcribed from `git show 9e8fb74^:pycircuit/circuit/transient.py` rather
    than reimplemented, so that what is measured is the code the note was written
    about.  The flags below are the only departures from it.
    """

    def run_old(self, tend, timestep, x0=None, refnode=gnd,
                keep_solved_h=True, analytical_eh=True, new_estimator=False,
                use_band=True):
        self.irefnode = self.cir.get_node_index(refnode)
        n = self.cir.n
        if x0 is None:
            from pycircuit.circuit.dcanalysis import DC
            x0 = DC(self.cir).solve().x

        x = x0
        self.base_integrator = self._get_integrator()
        hist_len = max(2, self.base_integrator.get_required_history())
        self._qlast = self.toolkit.array([self.cir.q(x) for _ in range(hist_len)])
        self._iqlast = self.toolkit.zeros((hist_len, n))

        X = [copy(x)]
        if hasattr(self.cir, 'accept_step'):
            self.cir.accept_step(0.0, X[-1], self.epar)
        timelist = []

        self._is_first_step = True
        t = 0.0
        h = timestep
        reltol = self.par.reltol
        xtol = self.par.vabstol
        abstol = self.toolkit.concatenate(
            (self.par.iabstol * self.toolkit.ones(len(self.cir.nodes)),
             self.par.vabstol * self.toolkit.ones(len(self.cir.branches))))

        h_trace = []
        dh_trace = []
        n_nonconv = 0
        n_at_floor = 0
        hmin = self.par.minstep

        while t < tend and len(h_trace) < MAX_STEPS:
            if t + h > tend:
                h = tend - t

            x_curr = copy(X[-1])
            h_curr = h

            def eval_FJ(xr, h_curr, _t=t):
                x_full = self.toolkit.zeros(n)
                x_full[:self.irefnode] = xr[:self.irefnode]
                x_full[self.irefnode + 1:] = xr[self.irefnode:]

                self._dt = h_curr
                C = self.cir.C(x_full)
                q = self.cir.q(x_full)
                iq, Geq = self.get_diff(q, C)
                F = self.cir.i(x_full) + iq + \
                    self.cir.u(_t + h_curr, analysis=self.par.analysis)
                J = self.cir.G(x_full) + Geq

                F_r = self.toolkit.delete(F, self.irefnode)
                J_r = self.toolkit.delete(J, self.irefnode, axis=0)
                J_r = self.toolkit.delete(J_r, self.irefnode, axis=1)

                ## p = df_ckt/dh, by central difference on the step size.
                eps = max(1e-8 * h_curr, 1e-15)
                self._dt = h_curr + eps
                iq_p, _ = self.get_diff(q, C)
                F_p = self.cir.i(x_full) + iq_p + \
                    self.cir.u(_t + h_curr + eps, analysis=self.par.analysis)
                self._dt = h_curr - eps
                iq_m, _ = self.get_diff(q, C)
                F_m = self.cir.i(x_full) + iq_m + \
                    self.cir.u(_t + h_curr - eps, analysis=self.par.analysis)
                J_h_r = self.toolkit.delete((F_p - F_m) / (2 * eps), self.irefnode)
                self._dt = h_curr

                E = self._calc_E(x_full, h_curr, n, reltol, new_estimator)

                ## THE BAND SHORT-CIRCUIT, and the point of the whole method: the
                ## step size is NOT solved exactly.  Once the LTE is anywhere
                ## inside [gamma_min, gamma_max]*TRTOL the step is good enough, so
                ## the LTE equation is replaced by `0 = 1 * dh`, which pins dh at
                ## zero and lets the Newton iteration finish on the circuit
                ## equations alone.  Without it the iteration keeps chasing an
                ## exact solve of f_lte and walks h down until it underflows --
                ## which is what the first version of this script measured, having
                ## dropped this block in transcription.
                if (use_band and analytical_eh and not self._is_first_step
                        and (BAND_MIN - 1.0) * TRTOL <= E <= (BAND_MAX - 1.0) * TRTOL):
                    return F_r, J_r, J_h_r, 0.0, self.toolkit.zeros(len(F_r)), 1.0

                if self._is_first_step:
                    return F_r, J_r, J_h_r, 0.0, self.toolkit.zeros(len(F_r)), 1.0

                ## d = df_lte/dh.  The LTE goes like h^p, so this is analytic.
                p_ord = 3.0 if self._effective_method in ('trapezoidal', 'trap') else 2.0
                E_h = p_ord * (E + TRTOL) / h_curr
                E_x_r = self.toolkit.zeros(len(F_r))
                if not analytical_eh:
                    ## q^T = df_lte/dv by finite difference -- O(n) extra residual
                    ## evaluations, which is why the default dropped it.
                    for i in range(len(F_r)):
                        idx = i if i < self.irefnode else i + 1
                        xp = copy(x_full); xp[idx] += eps
                        xm = copy(x_full); xm[idx] -= eps
                        E_x_r[i] = (self._calc_E(xp, h_curr, n, reltol, new_estimator)
                                    - self._calc_E(xm, h_curr, n, reltol, new_estimator)) \
                            / (2 * eps)
                return F_r, J_r, J_h_r, E, E_x_r, E_h

            def limiter_func(xr_next, xr_curr):
                a = self.toolkit.insert(xr_next, self.irefnode, 0.0)
                b = self.toolkit.insert(xr_curr, self.irefnode, 0.0)
                a = self.cir.limit(a, b, self.epar)
                return self.toolkit.concatenate(
                    (a[:self.irefnode], a[self.irefnode + 1:]))

            x_curr_r = self.toolkit.concatenate(
                (X[-1][:self.irefnode], X[-1][self.irefnode + 1:]))
            abstol_r = self.toolkit.concatenate(
                (abstol[:self.irefnode], abstol[self.irefnode + 1:]))

            try:
                (x_next_r, h_next), _ = SchurCoupledNewton().solve_system(
                    (x_curr_r, h), eval_FJ, self.toolkit, reltol, abstol_r,
                    xtol, self.par.maxiter, limiter=limiter_func,
                    hmin=hmin)
                converged = True
            except NoConvergenceError:
                converged = False

            if not converged:
                n_nonconv += 1
                h *= 0.25
                if h < 1e-15:
                    return dict(h=np.array(h_trace), dh=np.array(dh_trace),
                                t=t, tend=tend,
                                nonconv=n_nonconv, died='h < 1e-15')
                continue

            x_next = self.toolkit.zeros(n)
            x_next[:self.irefnode] = x_next_r[:self.irefnode]
            x_next[self.irefnode + 1:] = x_next_r[self.irefnode:]
            x_curr = x_next

            t += h_curr
            h_trace.append(h_curr)
            dh_trace.append(h_next - h_curr)
            ## Flooring the solved step stops it underflowing, but a step pinned
            ## AT the floor is not progress -- it is the same collapse, bounded.
            ## Reported as its own outcome so the two are not confused.
            n_at_floor = n_at_floor + 1 if h_curr <= hmin * (1.0 + 1e-9) else 0
            if n_at_floor >= STALL_STEPS:
                return dict(h=np.array(h_trace), dh=np.array(dh_trace),
                            t=t, tend=tend, nonconv=n_nonconv,
                            died='stalled at minstep')
            ## THE CONFOUND.  The original wrote `h = h_curr` here, discarding
            ## the step the coupled system had just solved for.
            h = h_next if keep_solved_h else h_curr
            if h <= 0.0:
                return dict(h=np.array(h_trace), dh=np.array(dh_trace),
                            t=t, tend=tend,
                            nonconv=n_nonconv, died='h <= 0 (solved h was %g)' % h_next)

            timelist.append(t)
            X.append(copy(x_curr))
            if hasattr(self.cir, 'accept_step'):
                self.cir.accept_step(t, X[-1], self.epar)

            self._dt = h_curr
            self._dt_last = h_curr
            self._is_first_step = False
            self._iqlast = self.toolkit.concatenate(
                (self.toolkit.array([self._iq]), self._iqlast))[:-1]
            self._qlast = self.toolkit.concatenate(
                (self.toolkit.array([self.cir.q(x_curr)]), self._qlast))[:-1]

        return dict(h=np.array(h_trace), dh=np.array(dh_trace), t=t, tend=tend,
                    nonconv=n_nonconv,
                    died=None if t >= tend else 'hit MAX_STEPS')

    def _calc_E(self, x_val, h_val, n, reltol, new_estimator):
        """Normalised LTE minus TRTOL, in the old code's convention.

        `new_estimator=True` swaps in the post-stage-4 charge-based third divided
        difference and its corrected Gear-2 constant, keeping everything else the
        same -- so the two runs differ in the estimator alone.
        """
        if self._is_first_step:
            return 0.0
        q_val = self.cir.q(x_val)
        abstol = self.par.iabstol

        if new_estimator:
            from pycircuit.circuit import _lte_kernels as K
            h_last = getattr(self, '_dt_last', h_val)
            ## `_dt_last2` EXISTS and is None before three points accumulate,
            ## so `getattr(..., default)` returns None rather than the default.
            h_last2 = getattr(self, '_dt_last2', None) or h_last
            dd3 = K.third_divided_difference(
                q_val, self._qlast, h_val, h_last, h_last2)
            if self._effective_method in ('trapezoidal', 'trap'):
                lte = K.trapezoidal_lte(h_val, dd3)
            else:
                lte = K.gear2_lte(h_val, h_last, dd3)
        elif self._effective_method in ('trapezoidal', 'trap'):
            dd2 = ((q_val - self._qlast[0]) / h_val - self._iqlast[0]) * 2.0 / h_val
            lte = 1.0 / 12.0 * (h_val ** 3) * dd2
        elif self._effective_method == 'gear2':
            dd1_n = (q_val - self._qlast[0]) / h_val
            dd1_nm1 = (self._qlast[0] - self._qlast[1]) / getattr(self, '_dt_last', h_val)
            dd2_n = (dd1_n - dd1_nm1) / (h_val + getattr(self, '_dt_last', h_val))
            lte = (h_val ** 2) * (h_val + getattr(self, '_dt_last', h_val)) / 3.0 * dd2_n
        else:
            ## THE OTHER CONFOUND: the default method fell here, so the LTE was
            ## identically zero and E was the constant -TRTOL.
            lte = self.toolkit.zeros(n)

        etol = reltol * self.toolkit.maximum(np.abs(q_val), np.abs(self._qlast[0])) + abstol
        return float(np.max(np.abs(lte) / etol)) - TRTOL


def _report(label, r):
    h = r['h']
    if len(h) == 0:
        print('%-52s  no steps taken (%s)' % (label, r['died']))
        return
    frac = 100.0 * r['t'] / r['tend']
    ## The FINAL step is `tend - t`, i.e. whatever residue is left, and on a run
    ## that finishes it is routinely 15 orders of magnitude below the rest.
    ## Reporting it as "h last" made a completed run look like a collapse in the
    ## first version of this script -- so the interior steps are what is summarised.
    interior = h[:-1] if len(h) > 1 else h
    ## Does the coupled solve move the step at all?  If `dh` is identically zero
    ## the (N+1) system is not adapting and nothing about collapse can be read
    ## from the run either way.
    dh = r['dh']
    moved = int(np.count_nonzero(np.abs(dh) > 0.0))
    print('%-52s %7d %8.1f%% %10.3e %10.3e %6d  %s'
          % (label, len(h), frac, np.min(interior), np.median(interior),
             moved, r['died'] or 'completed'))


def main():
    circuits = dict(S._fa_circuits())
    print(__doc__.strip().splitlines()[0])
    print()
    print('interior steps only (the final step is the tend residue, not a choice)')
    print('%-52s %7s %8s %10s %10s %6s  %s'
          % ('configuration', 'steps', 'reached', 'h min', 'h median',
             'dh!=0', 'outcome'))

    for cname in ('rc-vsin', 'stiff-rlc'):
        print('--- %s' % cname)
        from pycircuit.circuit.integrator import (TrapezoidalIntegrator,
                                                  Gear2Integrator)
        for method, make in (('trap', TrapezoidalIntegrator),
                             ('gear2', Gear2Integrator)):
            ## `analytical_eh=True` DROPS q^T (the old default set E_x to zero),
            ## so it is not the exact (N+1) system at all -- it is Fang sec. 3.4
            ## with the LTE gradient omitted.  Both are run, or "the exact update
            ## collapses" would be a claim about a system missing a row.
            ## `band=no` is the configuration the first run of this script
            ## measured by accident; it is kept as the control.
            for keep, aeh, new_est, band in (
                    (True, True, True, False),
                    (True, True, True, True),
                    (True, True, False, True),
                    (True, False, True, True),
                    (False, True, True, True)):
                    cir, kw = circuits[cname]()
                    tran = OldCoupled(cir, toolkit=numeric, reltol=1e-5,
                                      integrator=make())
                    label = '  %-5s band=%-4s keep_h=%-5s q^T=%-4s est=%s' % (
                        method, 'yes' if band else 'NO', keep,
                        'no' if aeh else 'FD', 'stage4' if new_est else 'old')
                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore')
                            r = tran.run_old(kw['tend'], kw['timestep'],
                                             x0=kw.get('x0'),
                                             keep_solved_h=keep,
                                             analytical_eh=aeh,
                                             new_estimator=new_est,
                                             use_band=band)
                        _report(label, r)
                    except Exception as e:
                        print('%-52s  raised %s: %s'
                              % (label, type(e).__name__, str(e)[:52]))
        print()

    print('GATE 12B-0 asks whether "the exact (N+1) Schur update collapses the')
    print('step size" is a property of the update or of the loop around it.')
    print('Compare keep_h=False (the shipped behaviour, which discarded the')
    print('solved step) against keep_h=True (actually using it), and the old')
    print('estimator against the post-stage-4 one.')


if __name__ == '__main__':
    main()
