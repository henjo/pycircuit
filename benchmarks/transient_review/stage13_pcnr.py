"""GATES 13-4 and 13-5 -- what PCNR costs, and whether it converges as often.

13-4 declared: the `Ax=b` actually solved stays the size of the original MNA
system, and the per-iteration cost stays within 15% of the classic path.
13-5 declared: everything that converges with limiting must still converge.

**Cost is measured per NEWTON ITERATION, not per solve.** That unit correction
was forced on gate 12-3, where dividing wall clock by time points mixed
cost-per-unit-work with work-per-point and could not separate them.

Run:
    MPLBACKEND=Agg PYTHONPATH=. python benchmarks/transient_review/stage13_pcnr.py
"""
import time
import warnings

import numpy as np

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, VS, IS, Diode
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit import pcnr as P
import pycircuit.circuit.nrsolver as nrs


def fig1(is1=1e-15, is2=1e-9):
    """The paper's own example: two parallel diodes, six decades apart."""
    c = SubCircuit()
    c['vs'] = VS('e1', gnd, v=1.0)
    c['R'] = R('e1', 'e2', r=1e3)
    c['D1'] = Diode('e2', gnd, IS=is1)
    c['D2'] = Diode('e2', gnd, IS=is2)
    return c


def current_driven():
    """1 A into a diode -- `test_dc_pcnr`'s circuit, which needs limiting."""
    c = SubCircuit()
    c['is'] = IS(gnd, 1, i=1.0)
    c['D1'] = Diode(1, gnd)
    return c


def hard_drive(v=20.0):
    """A large forward drive: without limiting the first step overshoots hugely."""
    c = SubCircuit()
    c['vs'] = VS('a', gnd, v=v)
    c['R'] = R('a', 'b', r=10.0)
    c['D'] = Diode('b', gnd, IS=1e-15)
    return c


def chain(n=4):
    """`n` diodes in series -- each one limits, and they interact."""
    c = SubCircuit()
    c['vs'] = VS('n0', gnd, v=5.0)
    c['R'] = R('n0', 'n1', r=1e3)
    for k in range(n):
        c['D%d' % k] = Diode('n%d' % (k + 1), 'n%d' % (k + 2) if k < n - 1 else gnd,
                             IS=1e-15)
    return c


def parallel_many(m=6):
    """`m` parallel diodes with spread saturation currents -- the clash case,
    scaled up: every one of them limits the same branch voltage."""
    c = SubCircuit()
    c['vs'] = VS('e1', gnd, v=1.0)
    c['R'] = R('e1', 'e2', r=100.0)
    for k in range(m):
        c['D%d' % k] = Diode('e2', gnd, IS=10.0 ** (-16 + k))
    return c


def parallel_big(m=60):
    """Enough devices that the linear algebra is not lost in Python overhead.

    The five small circuits above answer gate 13-5 well and gate 13-4 badly: at a
    handful of unknowns the per-iteration time is dominated by interpreter
    overhead, so it says nothing about what the Schur complement costs.
    """
    c = SubCircuit()
    c['vs'] = VS('e1', gnd, v=1.0)
    for k in range(m):
        c['R%d' % k] = R('e1', 'n%d' % k, r=100.0)
        c['D%d' % k] = Diode('n%d' % k, gnd, IS=10.0 ** (-16 + (k % 8)))
    return c


CASES = [('fig1 (2 parallel)', fig1),
         ('current driven', current_driven),
         ('hard drive 20V', hard_drive),
         ('series chain x4', chain),
         ('6 parallel', parallel_many),
         ('60 diodes', parallel_big)]


REPEATS = 5


def run_classic(build, repeats=REPEATS):
    """Returns (ok, iterations, wall, value).  Best-of-N: the minimum is the run
    least disturbed by whatever else the machine was doing, and scheduling noise
    only ever adds time."""
    best = None
    for _ in range(repeats):
        r = _run_classic_once(build)
        if best is None or (r[0] and r[2] < best[2]):
            best = r
    return best


def _run_classic_once(build):
    count = {'n': 0}
    real = nrs.StandardNewton.solve_system

    def patched(self, x0, func, *a, **k):
        def spy(x, *aa, **kk):
            count['n'] += 1
            return func(x, *aa, **kk)
        return real(self, x0, spy, *a, **k)

    nrs.StandardNewton.solve_system = patched
    try:
        c = build()
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = DC(c, toolkit=numeric).solve()
        wall = time.perf_counter() - t0
        node = [n for n in c.nodes if n.name != 'gnd!'][-1]
        return True, count['n'], wall, float(res.v(node))
    except Exception:
        return False, count['n'], float('nan'), float('nan')
    finally:
        nrs.StandardNewton.solve_system = real


def run_pcnr(build, repeats=REPEATS):
    best = None
    for _ in range(repeats):
        r = _run_pcnr_once(build)
        if best is None or (r[0] and r[2] < best[2]):
            best = r
    return best


def _run_pcnr_once(build):
    try:
        c = build()
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            x, _v, its = P.solve_dc(c, gnd, maxiter=500)
        wall = time.perf_counter() - t0
        node = [n for n in c.nodes if n.name != 'gnd!'][-1]
        return True, its, wall, float(x[c.get_node_index(node)])
    except Exception:
        return False, -1, float('nan'), float('nan')


def main():
    print('GATE 13-5 (does it converge where limiting does?) and 13-4 (what it costs)')
    print()
    print('%-20s %-8s %6s %8s %10s %14s' %
          ('circuit', 'method', 'ok', 'NRiter', 'us/iter', 'answer'))
    agree = disagree = 0
    for name, build in CASES:
        vals = {}
        for label, fn in (('limiting', run_classic), ('pcnr', run_pcnr)):
            ok, its, wall, val = fn(build)
            vals[label] = (ok, val)
            print('%-20s %-8s %6s %8s %10s %14s'
                  % (name, label, 'yes' if ok else 'NO',
                     its if its >= 0 else '-',
                     '%.1f' % (1e6 * wall / max(its, 1)) if ok else '-',
                     '%.9f' % val if ok else '-'))
        (ok_a, va), (ok_b, vb) = vals['limiting'], vals['pcnr']
        if ok_a and ok_b:
            if abs(va - vb) <= 1e-6 * max(abs(va), 1.0):
                agree += 1
            else:
                disagree += 1
                print('%-20s %s' % ('', '^^ ANSWERS DIFFER by %.3e' % abs(va - vb)))
        print()
    print('answers agreeing: %d   disagreeing: %d' % (agree, disagree))
    print()
    print('Read `us/iter` only on the largest case. At a handful of unknowns it is')
    print('interpreter overhead, not linear algebra, and says nothing about what')
    print('the Schur complement costs. Timings are best-of-%d.' % REPEATS)


if __name__ == '__main__':
    main()
