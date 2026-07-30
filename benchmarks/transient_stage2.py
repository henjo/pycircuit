"""STAGE 2 gates: performance, without touching numerics.

`doc/transient_work_plan.md` stage 2.  Every sub-stage is independently gated so a
regression is attributable, and the whole stage is defined as **behaviour
preserving**: the declared stop condition is that waveform drift must be exactly
`0.00e+00`, not "small".

Usage::

    # record a reference waveform before changing anything
    python -u benchmarks/transient_stage2.py --save-ref

    # after a change: timing + drift against that reference
    python -u benchmarks/transient_stage2.py

    # per-step assembly counts (gate 2c)
    python -u benchmarks/transient_stage2.py --counts

    # exact-equality check of i/q/u/G/C against a reference implementation (gate 2b)
    python -u benchmarks/transient_stage2.py --stamps

TIMING METHODOLOGY, and it is not optional on this machine.  Stage 1 recorded up to
+-57% run-to-run variation on individual tests, and a single-sample comparison there
read as a 35% regression that turned out to be 3%.  Every timing here is therefore
**min of N repetitions** (min, not mean: it is the sample least polluted by
scheduler noise), and any claimed speedup should be sanity-checked against a
non-timing invariant -- step count, assembly count, or profile call count -- before
it is believed.
"""

import argparse
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
for _p in (_REPO, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import pycircuit.circuit.circuit as circuit_module
from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit.benchmark_circuits import build_leapfrog_network
from pycircuit.circuit.elements import BSource
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import TrapezoidalIntegrator

F1 = 8e3
G3 = 1e-3
STAGES = 5

## Long enough to be dominated by the steady per-step cost rather than by setup,
## short enough to repeat several times.  The IM3 harness's own tolerances, so the
## step-size behaviour matches the circuit this stage is measured on elsewhere.
TEND = 3.0e-4
TIMESTEP = 1.0 / (F1 * 128)
TRAN_OPTS = dict(vabstol=1e-9, reltol=1e-6)

REF_PATH = os.path.join(
    os.environ.get('STAGE2_REF_DIR', '/tmp'), 'stage2_reference.npz')


def build():
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = numeric
    try:
        cir = SubCircuit(toolkit=numeric)
        node, amp, _names = build_leapfrog_network(
            cir, stages=STAGES, tones=((1e-3, F1),))
        cir['nl'] = BSource(amp[0]['e1'], gnd, amp[0]['e1'], gnd,
                            i_func=lambda v: G3 * v ** 3)
        return cir, amp
    finally:
        circuit_module.default_toolkit = saved


def run_once():
    cir, amp = build()
    out = amp[STAGES - 1]['out']
    opts = dict(TRAN_OPTS)
    opts['integrator'] = TrapezoidalIntegrator()
    t0 = time.perf_counter()
    res = Transient(cir, toolkit=numeric, **opts).solve(
        refnode=gnd, tend=TEND, timestep=TIMESTEP)
    secs = time.perf_counter() - t0
    wave = res.v(out)
    return secs, np.asarray(wave.x[0]), np.asarray(wave.y)


def timed(reps):
    best = None
    times = []
    for _ in range(reps):
        secs, t, v = run_once()
        times.append(secs)
        if best is None:
            best = (t, v)
    return min(times), times, best[0], best[1]


def cmd_save_ref(args):
    secs, times, t, v = timed(args.reps)
    np.savez(REF_PATH, t=t, v=v, secs=secs)
    print('reference saved to %s' % REF_PATH)
    print('  steps      %d' % len(t))
    print('  runtime    %.3f s  (min of %d: %s)'
          % (secs, args.reps, ' '.join('%.2f' % x for x in times)))
    print('  v(out)[-1] %.17g' % v[-1])
    return 0


def cmd_compare(args):
    if not os.path.exists(REF_PATH):
        print('no reference at %s -- run with --save-ref first' % REF_PATH)
        return 1
    ref = np.load(REF_PATH)
    secs, times, t, v = timed(args.reps)

    print('STAGE 2 -- timing and drift against the recorded reference')
    print('  reference  %.3f s, %d steps' % (float(ref['secs']), len(ref['t'])))
    print('  now        %.3f s, %d steps  (min of %d: %s)'
          % (secs, len(t), args.reps, ' '.join('%.2f' % x for x in times)))
    print('  SPEEDUP    %.2fx' % (float(ref['secs']) / secs))
    print()

    ## Drift is only meaningful if the solver took the same path.  A different
    ## step count means the numerics changed, which this stage forbids outright --
    ## report it rather than interpolating the difference away.
    if len(t) != len(ref['t']):
        print('  STEP COUNT CHANGED: %d -> %d' % (len(ref['t']), len(t)))
        print('  Stage 2 is defined as behaviour-preserving; a changed step count')
        print('  means something other than performance changed.  STOP.')
        return 1

    dt = np.max(np.abs(t - ref['t']))
    dv = np.max(np.abs(v - ref['v']))
    print('  max |t - t_ref|  %.2e s' % dt)
    print('  max |v - v_ref|  %.2e V' % dv)
    exact = (dt == 0.0) and (dv == 0.0)
    print('  DRIFT EXACTLY ZERO: %s' % ('YES' if exact else 'NO'))
    if not exact:
        print('  The stage gate says STOP if drift is not exactly 0.00e+00.')
        return 1
    return 0


def cmd_counts(args):
    """Gate 2c: assemblies per accepted step."""
    cir, amp = build()
    counts = {}

    def wrap(name):
        real = getattr(cir, name)

        def counting(*a, **k):
            counts[name] = counts.get(name, 0) + 1
            return real(*a, **k)
        setattr(cir, name, counting)

    for name in ('i', 'q', 'u', 'G', 'C'):
        wrap(name)

    opts = dict(TRAN_OPTS)
    opts['integrator'] = TrapezoidalIntegrator()
    res = Transient(cir, toolkit=numeric, **opts).solve(
        refnode=gnd, tend=TEND, timestep=TIMESTEP)
    steps = len(np.asarray(res.v(amp[STAGES - 1]['out']).x[0]))

    print('STAGE 2c -- assemblies per accepted step  (%d steps)' % steps)
    for name in ('G', 'C', 'i', 'u', 'q'):
        print('  %-2s  %7d calls   %6.2f per step'
              % (name, counts.get(name, 0), counts.get(name, 0) / steps))
    return 0


## ---------------------------------------------------------------------------
## The pre-stage-2 assembly, verbatim from `circuit.py` at 4156c0a.  It lives here
## rather than in the library because it is scaffolding for one gate: the point of
## gate 2b is that the rewritten assembly produces *bit-identical* stamps, and that
## can only be shown against the code it replaced.
## ---------------------------------------------------------------------------

def reference_subvectors(self, methodname, x, args, dtype=None, params_tree=None):
    n = self.n
    lhs = self.toolkit.zeros(n, dtype=dtype)

    batched = self.toolkit.batched_contributions(
        self, getattr(self, '_eval_groups', None), methodname, x, args,
        params_tree, ndim=1)
    if batched is not None:
        values, indices = batched
        lhs = self.toolkit.add_at(lhs, indices, values)

    for instance, element in self.elements.items():
        if getattr(self, '_eval_groups', None) and element.__class__ in self._eval_groups and x is not None:
            continue
        if x is not None:
            subx = x[self.elementnodemap[instance]]
            rhs = getattr(element, methodname)(subx, *args)
        else:
            rhs = getattr(element, methodname)(*args)

        if instance in self._map_indices_1d:
            indices = self._map_indices_1d[instance]
            if hasattr(self.toolkit, 'add_at'):
                rhs_flat = self.toolkit.reshape(rhs, (-1,)).flatten()
                lhs = self.toolkit.add_at(lhs, indices, rhs_flat)
            else:
                rhs_flat = np.asarray(rhs).flatten()
                try:
                    np.add.at(lhs, indices, rhs_flat)
                except TypeError:
                    lhs = lhs.astype(object)
                    np.add.at(lhs, indices, rhs_flat)
    return lhs


def reference_submatrices(self, methodname, x, args, params_tree=None):
    n = self.n
    build_sparse = hasattr(self.toolkit, 'build_sparse')
    if build_sparse:
        all_data, all_rows, all_cols = [], [], []
    else:
        lhs = self.toolkit.zeros((n, n))

    batched = self.toolkit.batched_contributions(
        self, getattr(self, '_eval_groups', None), methodname, x, args,
        params_tree, ndim=2)
    if batched is not None:
        values, (rows, cols) = batched
        if build_sparse:
            all_data.append(values); all_rows.append(rows); all_cols.append(cols)
        else:
            lhs = self.toolkit.add_at(lhs, (rows, cols), values)

    for instance, element in self.elements.items():
        if getattr(self, '_eval_groups', None) and element.__class__ in self._eval_groups:
            continue
        nodemap = self.elementnodemap[instance]
        if x is not None:
            rhs = getattr(element, methodname)(x[nodemap], *args)
        else:
            rhs = getattr(element, methodname)(*((None,) + tuple(args)))

        if instance in self._map_indices_2d:
            rows, cols = self._map_indices_2d[instance]
            if build_sparse:
                all_data.append(np.asarray(rhs).flatten())
                all_rows.append(rows); all_cols.append(cols)
            elif hasattr(self.toolkit, 'add_at'):
                rhs_flat = self.toolkit.reshape(rhs, (-1,)).flatten()
                lhs = self.toolkit.add_at(lhs, (rows, cols), rhs_flat)
            else:
                rhs_flat = np.asarray(rhs).flatten()
                try:
                    np.add.at(lhs, (rows, cols), rhs_flat)
                except TypeError:
                    lhs = lhs.astype(object)
                    np.add.at(lhs, (rows, cols), rhs_flat)

    if build_sparse:
        if not all_data:
            return self.toolkit.build_sparse([], [], [], shape=(n, n))
        return self.toolkit.build_sparse(
            np.concatenate(all_data), np.concatenate(all_rows),
            np.concatenate(all_cols), shape=(n, n))
    return lhs


def cmd_stamps(args):
    """Gate 2b: every stamp must be bit-identical to the code it replaced."""
    from pycircuit.circuit.circuit import defaultepar

    cir, _amp = build()
    rng = np.random.default_rng(0)
    worst = {}
    nonzero = {}

    ## Random states rather than one: a single state can agree by accident (all
    ## contributions landing on distinct indices, so no accumulation order exists
    ## to get wrong). Several states with a nonlinear element in the circuit
    ## exercise the duplicate-index path that the scatter rewrite actually changes.
    for _trial in range(8):
        x = rng.normal(scale=0.3, size=cir.n)
        cases = [
            ('i', 1, x, (defaultepar,)),
            ('q', 1, x, (defaultepar,)),
            ## Sampled at a quarter period, NOT at t=0: a sine source is zero at
            ## the origin, so `u` there is an all-zero vector and compares equal
            ## no matter what the assembly does.  The first version of this gate
            ## did exactly that and reported "EXACT" against 0 nonzeros.
            ('u', 1, None, (0.25 / F1, defaultepar, 'tran')),
            ('G', 2, x, (defaultepar,)),
            ('C', 2, x, (defaultepar,)),
        ]
        for name, ndim, state, extra in cases:
            if ndim == 1:
                fast = cir._add_element_subvectors(name, state, extra)
                slow = reference_subvectors(cir, name, state, extra)
            else:
                fast = cir._add_element_submatrices(name, state, extra)
                slow = reference_submatrices(cir, name, state, extra)
            fast = np.asarray(fast, dtype=float)
            slow = np.asarray(slow, dtype=float)
            d = float(np.max(np.abs(fast - slow)))
            worst[name] = max(worst.get(name, 0.0), d)
            nonzero[name] = max(nonzero.get(name, 0), int(np.count_nonzero(slow)))

    print('STAGE 2b -- rewritten vs pre-stage-2 assembly, 8 random states')
    print('%-4s %14s %12s  %s' % ('', 'max abs diff', 'nonzeros', 'verdict'))
    ok = True
    for name in ('i', 'q', 'u', 'G', 'C'):
        print('%-4s %14.3e %12d  %s'
              % (name, worst[name], nonzero[name],
                 'EXACT' if worst[name] == 0.0 else 'DIFFERS'))
        ok = ok and worst[name] == 0.0
    ## A stamp of all zeros would compare equal trivially; report the nonzero
    ## count so "exact" cannot be read off an empty matrix.
    print('  ALL EXACTLY 0.0: %s' % ('YES' if ok else 'NO'))
    return 0 if ok else 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--save-ref', action='store_true')
    p.add_argument('--counts', action='store_true')
    p.add_argument('--stamps', action='store_true')
    p.add_argument('--reps', type=int, default=3)
    args = p.parse_args()
    if args.save_ref:
        return cmd_save_ref(args)
    if args.counts:
        return cmd_counts(args)
    if args.stamps:
        return cmd_stamps(args)
    return cmd_compare(args)


if __name__ == '__main__':
    sys.exit(main())
