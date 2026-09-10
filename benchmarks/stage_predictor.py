"""The stage predictor's cost measurement, every integrator family.

⚠⚠ THE INSTRUMENT IS THE HARD PART, and two obvious fixtures cannot measure a
stage predictor AT ALL:

* `Diode` linearises its `G` around a STORED `_vlim` (its own docstring, stage
  13-2), so its Newton is seed-blind -- measured, an ALL-ZEROS seed gives the
  iteration histogram [436, 244, 108], to the count, that the exact seed gives;
* a LINEAR circuit forces a one-step Newton, so the index-2 C-V loop the GLM
  work used elsewhere reads no change either.

`ExpG` below is a state-free exponential: `i` and `G` are both functions of the
passed `x`, with no stored state and no limiting.

Run it: `python benchmarks/stage_predictor.py`.
"""
import warnings

import numpy as np
import sympy

import pycircuit.circuit.circuit
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import C, R, VSin
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, Parameter
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import transient as TR
from pycircuit.circuit import nrsolver as NR
from pycircuit.circuit.integrator import (RadauIIA3Integrator,
                                          ESDIRK43Integrator,
                                          TRBDF2Integrator, Gear2Integrator,
                                          TrapezoidalIntegrator,
                                          GLM3Integrator, GLM4Integrator)
from pycircuit.circuit.transient import Transient  # noqa: F811

warnings.simplefilter('ignore')
pycircuit.circuit.circuit.default_toolkit = numeric
PER = 1e-3


class ExpG(Behavioural):
    """i = IS (exp(V/VT) - 1), state-free, no limiting."""
    instparams = [Parameter(name='IS', desc='sat', unit='A', default=1e-12),
                  Parameter(name='VT', desc='thermal', unit='V', default=0.026)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return (Contribution(b.I, IS * (sympy.exp(b.V / VT) - 1)),)  # noqa: F821


def fixture(va=0.8):
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=va, freq=1.0 / PER)
    c['rs'] = R('a', 'b', r=50.0)
    c['nl'] = ExpG('b', gnd, IS=1e-12, VT=0.026)
    c['cl'] = C('b', gnd, c=1e-9)
    c['rl'] = R('b', gnd, r=1e4)
    return c


def measure(cls, mode, N, va=0.8, fixed=True):
    log = []
    solve_orig = NR.StandardNewton.solve_system

    def wrapped(self, x0, eval_FJ, *a, **k):
        x, it = solve_orig(self, x0, eval_FJ, *a, **k)
        log.append(it)
        return x, it
    NR.StandardNewton.solve_system = wrapped
    TR.Transient.stage_predictor = mode
    n = {'i': 0}
    cir = fixture(va)
    fi = cir.i
    cir.i = lambda *a, **k: (n.__setitem__('i', n['i'] + 1), fi(*a, **k))[1]
    try:
        tr = Transient(cir, integrator=cls(), reltol=1e-9)
        res = tr.solve(refnode=gnd, tend=PER, timestep=PER / N,
                       fixed_timestep=fixed)
    finally:
        NR.StandardNewton.solve_system = solve_orig
        TR.Transient.stage_predictor = 'on'
    wv = res.v('b')
    return (n['i'], int(max(log)) if log else 0,
            tr.statistics.accepted_steps,
            np.asarray(wv.x[0], dtype=float), np.asarray(wv.y, dtype=float))


FAMILIES = [(RadauIIA3Integrator, 'radau'), (ESDIRK43Integrator, 'esdirk43'),
            (TRBDF2Integrator, 'trbdf2'), (Gear2Integrator, 'gear2'),
            (TrapezoidalIntegrator, 'trap'), (GLM3Integrator, 'glm3'),
            (GLM4Integrator, 'glm4')]


def table():
    for fixed in (False, True):
        print('=== %s ===' % ('ADAPTIVE step (the default path)' if not fixed
                              else 'FIXED step'))
        print('%-9s %5s %4s %8s %8s %8s %9s %10s'
              % ('method', 'va', 'N', 'i off', 'i on', 'change', 'maxit',
                 'drift'))
        for cls, name in FAMILIES:
            if not fixed and name.startswith('glm'):
                continue          # the Nordsieck state needs a constant step
            for va, N in ((0.8, 200), (2.0, 40)):
                i0, m0, s0, t0, y0 = measure(cls, 'off', N, va, fixed)
                i1, m1, s1, t1, y1 = measure(cls, 'on', N, va, fixed)
                if t0.shape == t1.shape and np.allclose(t0, t1):
                    d = float(np.max(np.abs(y1 - y0)))
                else:
                    d = float(np.max(np.abs(np.interp(t0, t1, y1) - y0)))
                print('%-9s %5.1f %4d %8d %8d %7.1f%% %4d/%-4d %10.2e%s'
                      % (name, va, N, i0, i1, 100.0 * (i1 - i0) / i0, m0, m1,
                         d, '' if s0 == s1 else '  steps %d->%d' % (s0, s1)))
        print()


def clamp_sweep():
    """Why `PRED_CLAMP` is 1.5.  The prediction may move at most this many
    linear-extrapolation displacements from the newest node."""
    print('=== PRED_CLAMP sweep (change in device evaluations / worst iters) ===')
    cases = [(Gear2Integrator, 'gear2', 2.0, 40, True),
             (Gear2Integrator, 'gear2', 0.8, 200, True),
             (Gear2Integrator, 'gear2', 0.8, 200, False),
             (RadauIIA3Integrator, 'radau', 0.8, 200, False),
             (ESDIRK43Integrator, 'esdirk43', 2.0, 40, True),
             (GLM4Integrator, 'glm4', 0.8, 200, True)]
    Ks = (1.0, 1.5, 2.0, 3.0, 1e9)
    print('%-9s %4s %4s %-4s %9s | %s'
          % ('method', 'va', 'N', 'stp', 'off',
             '  '.join('K=%-9s' % ('inf' if k > 1e8 else k) for k in Ks)))
    for cls, name, va, N, fixed in cases:
        i0, m0 = measure(cls, 'off', N, va, fixed)[:2]
        row = []
        for K in Ks:
            TR.Transient.PRED_CLAMP = K
            try:
                i1, m1 = measure(cls, 'on', N, va, fixed)[:2]
            finally:
                TR.Transient.PRED_CLAMP = 1.5
            row.append('%+6.1f%%/%-3d' % (100.0 * (i1 - i0) / i0, m1))
        print('%-9s %4.1f %4d %-4s %5d(%2d) | %s'
              % (name, va, N, 'fix' if fixed else 'adp', i0, m0,
                 '  '.join(row)))
    print()


if __name__ == '__main__':
    seed_headroom()
    table()
    clamp_sweep()


def seed_headroom():
    """How far each family's seed sits from the value it converges to, as a
    FRACTION OF ONE STEP'S OWN STATE MOTION.

    This is the measurement that says a stage predictor was worth building at
    all, and it is what the roadmap's "every other integrator seeds just as
    crudely" line is made of.  A fraction near 1 means the Newton starts a
    whole step behind; near 0.5 means half a step.

    ⚠ The number has to be RELATIVE to the step's motion, or it just measures
    how fast the circuit is moving.  Measured with the predictor OFF, since
    the point is what the old seed did.
    """
    print('=== seed distance / one step of state motion (predictor OFF) ===')
    TR.Transient.stage_predictor = 'off'
    print('%-9s %8s %8s %-34s %s'
          % ('method', 'i-calls', 'i/step', 'seed', 'median   max'))
    try:
        for cls, name, coupled, seed in (
                (RadauIIA3Integrator, 'radau', True, 'x_n for ALL stages'),
                (ESDIRK43Integrator, 'esdirk43', False, 'previous stage'),
                (TRBDF2Integrator, 'trbdf2', False, 'previous stage')):
            rec = []
            meth = '_rk_step_coupled' if coupled else '_rk_step_dirk'
            fn = getattr(Transient, meth)

            def wrapped(self, x0, t, pf=None, _fn=fn, _c=coupled, _r=rec):
                out = _fn(self, x0, t, pf)
                Y = getattr(self, '_rk_Y', None)
                if Y:
                    xn = np.asarray(x0, dtype=float)
                    seeds = ([xn] * len(Y) if _c else
                             [xn] + [np.asarray(Y[j], dtype=float)
                                     for j in range(len(Y) - 1)])
                    sd = max(float(np.max(np.abs(np.asarray(Y[i], dtype=float)
                                                 - seeds[i])))
                             for i in range(len(Y)))
                    mv = float(np.max(np.abs(np.asarray(out[0], dtype=float)
                                             - xn)))
                    if mv > 0:
                        _r.append(sd / mv)
                return out
            setattr(Transient, meth, wrapped)
            n = {'i': 0}
            cir = fixture(0.8)
            fi = cir.i
            cir.i = lambda *a, **k: (n.__setitem__('i', n['i'] + 1),
                                     fi(*a, **k))[1]
            tr = Transient(cir, integrator=cls(), reltol=1e-9)
            try:
                tr.solve(refnode=gnd, tend=PER, timestep=PER / 200,
                         fixed_timestep=True)
            finally:
                setattr(Transient, meth, fn)
            r = np.array(rec)
            print('%-9s %8d %8.2f %-34s %.2f     %.2f'
                  % (name, n['i'], n['i'] / tr.statistics.accepted_steps, seed,
                     np.median(r), r.max()))
        ## the multistep family has no stages: its seed IS x_n, one step behind
        for cls, name in ((Gear2Integrator, 'gear2'),
                          (TrapezoidalIntegrator, 'trap')):
            n = {'i': 0}
            cir = fixture(0.8)
            fi = cir.i
            cir.i = lambda *a, **k: (n.__setitem__('i', n['i'] + 1),
                                     fi(*a, **k))[1]
            tr = Transient(cir, integrator=cls(), reltol=1e-9)
            tr.solve(refnode=gnd, tend=PER, timestep=PER / 200,
                     fixed_timestep=True)
            print('%-9s %8d %8.2f %-34s 1.00 by definition'
                  % (name, n['i'], n['i'] / tr.statistics.accepted_steps,
                     'x_n, one solve per step'))
    finally:
        TR.Transient.stage_predictor = 'on'
    print()
