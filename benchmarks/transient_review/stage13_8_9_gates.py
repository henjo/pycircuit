"""GATES 13-8 and 13-9 -- the skip set, then the Diode fix it makes safe.

13-8: `augmented_system` stops building the PCNR devices into the ordinary
assembly and subtracting them out again; `SubCircuit.i`/`G` take a `skip` set.
13-9: outside a classic-limiting scope `Diode.G`/`Diode.i` evaluate at the node
voltage instead of a stale `_vlim`, so `cir.G(x)` is an honest Jacobian.

See doc/pcnr_native_design.md section 5 (correction of 2026-11-16) for why the
order is forced, and doc/transient_work_plan.md gates 13-8/13-9 for what each
measurement is gating.

Three modes:

    snapshot OUT.npz        run every case, classic and PCNR, and save it
    compare  A.npz B.npz    bit-identity / max rel diff, case by case
    probes                  call counts, assembly cost, Jacobian honesty

Run from the repo root:
    PYTHONPATH=. MPLBACKEND=Agg python3 -u benchmarks/transient_review/stage13_8_9_gates.py probes
"""
import sys
import time
import warnings

import numpy as np

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VS, VSin, IS, Diode
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import pcnr as P

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import stage13_pcnr as S13  # noqa: E402  -- the six DC circuits of gates 13-4/13-5


def rect():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=5.0, freq=1e3)
    c['D'] = Diode('a', 'b', IS=1e-15)
    c['R'] = R('b', gnd, r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c


def mains():
    c = SubCircuit()
    c['VS'] = VSin(1, gnd, va=10.0, freq=50.0)
    c['D1'] = Diode(1, 2)
    c['C1'] = C(2, gnd, c=100e-6)
    c['R1'] = R(2, gnd, r=1e3)
    return c


## ---------------------------------------------------------------- counters

COUNT = {'G': 0, 'i': 0, 'limit': 0}
_real = {'G': Diode.G, 'i': Diode.i, 'limit': Diode.limit}


def _counting(name):
    real = _real[name]

    def wrapper(self, *a, **k):
        COUNT[name] += 1
        return real(self, *a, **k)
    return wrapper


def counting(on):
    for name in COUNT:
        setattr(Diode, name, _counting(name) if on else _real[name])
        COUNT[name] = 0


## ---------------------------------------------------------------- cases

def dc_classic(build):
    c = build()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = DC(c, toolkit=numeric).solve()
    return np.asarray(res.x, float).ravel()


def dc_pcnr(build):
    c = build()
    x, v_lim, its = P.solve_dc(c, gnd, maxiter=500)
    return np.concatenate([np.asarray(x, float), np.asarray(v_lim, float), [its]])


def tran(build, node, tend, timestep, pcnr, coupled=False):
    t = Transient(build(), toolkit=numeric, reltol=1e-5, pcnr=pcnr)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = t.solve(tend=tend, timestep=timestep, coupled_lte=coupled)
    w = res.v(node)
    return np.concatenate([np.asarray(w.x, float).ravel(),
                           np.asarray(w.y, float).ravel(),
                           [t.statistics.accepted_steps]])


def cases():
    out = []
    for name, build in S13.CASES:
        key = name.replace(' ', '_').replace('(', '').replace(')', '')
        out.append(('dc_classic_' + key, lambda b=build: dc_classic(b)))
        out.append(('dc_pcnr_' + key, lambda b=build: dc_pcnr(b)))
    for pc in (False, True):
        tag = 'pcnr' if pc else 'classic'
        out.append(('tran_%s_rect' % tag, lambda p=pc: tran(rect, 'b', 2e-3, 1e-5, p)))
        out.append(('tran_%s_mains' % tag, lambda p=pc: tran(mains, 2, 40e-3, 1e-4, p)))
        out.append(('coupled_%s_rect' % tag,
                    lambda p=pc: tran(rect, 'b', 2e-3, 1e-5, p, coupled=True)))
    return out


def snapshot(path):
    data = {}
    counting(True)
    for name, fn in cases():
        for k in COUNT:
            COUNT[k] = 0
        t0 = time.perf_counter()
        data[name] = fn()
        data[name + '__counts'] = np.array([COUNT['G'], COUNT['i'], COUNT['limit']])
        print('%-34s %8.2f s  len %6d  G/i/limit calls %s'
              % (name, time.perf_counter() - t0, data[name].size,
                 data[name + '__counts'].tolist()), flush=True)
    counting(False)
    np.savez(path, **data)


def compare(pa, pb):
    a, b = np.load(pa), np.load(pb)
    for name in sorted(k for k in a.files if not k.endswith('__counts')):
        x, y = a[name], b[name]
        ca, cb = a[name + '__counts'].tolist(), b[name + '__counts'].tolist()
        if x.shape != y.shape:
            verdict = 'SHAPE %s -> %s' % (x.shape, y.shape)
        elif np.array_equal(x, y):
            verdict = 'bit-identical'
        else:
            verdict = 'max rel diff %.3e' % (np.max(np.abs(x - y))
                                             / max(np.max(np.abs(x)), 1e-300))
        print('%-34s %-22s G/i/limit %s -> %s' % (name, verdict, ca, cb))


## ---------------------------------------------------------------- probes

def probe_counts():
    """13-8b: during a PCNR DC solve, no Diode.i/G call at all should remain --
    the device enters only through `pcnr_i`/`pcnr_didv` at its own unknown."""
    print('13-8b  Diode.G / Diode.i calls during P.solve_dc')
    counting(True)
    for name, build in S13.CASES:
        c = build()
        for k in COUNT:
            COUNT[k] = 0
        P.solve_dc(c, gnd, maxiter=500)
        print('   %-20s G %5d   i %5d' % (name, COUNT['G'], COUNT['i']))
    counting(False)


def probe_cost(repeats=200):
    """13-8c: `augmented_system` per call on 60 diodes, against classic assembly."""
    c = S13.parallel_big()
    junctions = P.pcnr_junctions(c)
    from pycircuit.circuit.circuit import defaultepar
    epar = defaultepar
    rng = np.random.default_rng(0)
    x = rng.uniform(0.0, 0.7, c.n)
    x[c.get_node_index(gnd)] = 0.0
    v_lim = np.array([x[ra] - x[rb] for _i, _e, ra, rb in junctions])

    def best(fn):
        fn()
        b = np.inf
        for _ in range(repeats):
            t0 = time.perf_counter()
            fn()
            b = min(b, time.perf_counter() - t0)
        return 1e6 * b

    classic = best(lambda: (c.i(x, epar), c.G(x, epar)))
    aug = best(lambda: P.augmented_system(c, x, v_lim, junctions, epar,
                                          dense_blocks=False))
    print('13-8c  60 diodes, best of %d: classic i+G %.0f us, augmented_system %.0f us, '
          'ratio %.2f' % (repeats, classic, aug, aug / classic))


def probe_honest():
    """13-9a: after a classic DC op has set `_vlim`, is `cir.G(x)` the true
    device Jacobian at a DIFFERENT x?  Before the fix it is g(_vlim)."""
    from pycircuit.circuit.circuit import defaultepar
    c = mains()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        DC(c, toolkit=numeric).solve()
    d = c['D1']
    print('13-9a  mains rectifier after a classic DC op; _vlim = %r'
          % getattr(d, '_vlim', None))
    VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron
    rows = c.elementnodemap['D1']
    worst = 0.0
    for v in (-18.0, -1.0, 0.3, 0.6, 0.75):
        x = np.zeros(c.n)
        x[rows[0]] = v
        x[rows[1]] = 0.0
        G = np.asarray(c.G(x, defaultepar), float)
        g_asm = G[rows[0], rows[0]]
        g_true = d.iparv.IS * np.exp(v / VT) / VT
        ## Node 1's diagonal holds only the diode: VS stamps node-branch
        ## off-diagonals, and C1/R1 sit on node 2.
        rel = abs(g_asm - g_true) / g_true
        worst = max(worst, rel)
        print('   v = %6.2f V   assembled g %.6e   true g %.6e   rel err %.3e'
              % (v, g_asm, g_true, rel))
    print('   worst rel err %.3e' % worst)


def probe_honest_after_pcnr():
    """13-9b: the same question at the END of a PCNR transient, which is the
    case gate 13-6 found: `_vlim` frozen while the junction swung."""
    from pycircuit.circuit.circuit import defaultepar
    c = mains()
    t = Transient(c, toolkit=numeric, reltol=1e-5, pcnr=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t.solve(tend=40e-3, timestep=1e-4)
    d = c['D1']
    VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron
    rows = c.elementnodemap['D1']
    x = np.zeros(c.n)
    x[rows[0]] = 0.7
    G = np.asarray(c.G(x, defaultepar), float)
    g_true = d.iparv.IS * np.exp(0.7 / VT) / VT
    print('13-9b  after a PCNR transient: _vlim = %r; at v = 0.7 V assembled g %.6e, '
          'true g %.6e, rel err %.3e'
          % (getattr(d, '_vlim', None), G[rows[0], rows[0]], g_true,
             abs(G[rows[0], rows[0]] - g_true) / g_true))


def probe_pss():
    """13-9f: PSS runs its own Newton with NO limiter, so `Diode.limit` is never
    called there and the diode is evaluated at whatever `_vlim` its first `G`/`i`
    call wrote -- frozen for the whole solve.  Compared against an independent
    reference: a classic transient run to steady state, last period."""
    from pycircuit.circuit.shooting import PSS
    seen = set()
    calls = {'n': 0, 'with_vlim': 0}
    real_G = Diode.G

    def spy_G(self, x, *a, **k):
        r = real_G(self, x, *a, **k)
        ## After the call, so the first call's own write is counted -- that
        ## write is the defect.  `nan` is never used as a sentinel: nan != nan,
        ## so a set of them counts every one as distinct.
        calls['n'] += 1
        if '_vlim' in self.__dict__:
            calls['with_vlim'] += 1
            seen.add(round(float(self._vlim), 9))
        return r

    for label, build, node, period, steps in (
            ('5 V half-wave rectifier', rect, 'b', 1e-3, 100),):
        Diode.G = spy_G
        seen.clear()
        calls.update(n=0, with_vlim=0)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = PSS(build()).solve(period=period, timestep=period / steps)
            v = np.asarray(res['tpss'].v(node, gnd), dtype=float).ravel()
            pss = (float(v.max()), float(v.min()))
            err = None
        except Exception as e:  # noqa: BLE001
            pss, err = None, '%s: %s' % (type(e).__name__, str(e)[:100])
        finally:
            Diode.G = real_G
        t = Transient(build(), toolkit=numeric, reltol=1e-6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r = t.solve(tend=20 * period, timestep=period / steps)
        w = r.v(node)
        tt = np.asarray(w.x, float).ravel()
        yy = np.asarray(w.y, float).ravel()
        last = yy[tt >= 19 * period]
        print('13-9f  %s: PSS %s; Diode.G calls %d, of which %d had a _vlim '
              '(distinct values: %s)'
              % (label, err or 'max %.6f min %.6f' % pss, calls['n'],
                 calls['with_vlim'], sorted(seen)[:4]))
        print('       reference (classic transient, period 20): max %.6f min %.6f'
              % (float(last.max()), float(last.min())))


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'snapshot':
        snapshot(sys.argv[2])
    elif mode == 'compare':
        compare(sys.argv[2], sys.argv[3])
    elif mode == 'probes':
        probe_counts()
        probe_cost()
        probe_honest()
        probe_honest_after_pcnr()
        probe_pss()
    elif mode == 'pss':
        probe_pss()
    else:
        raise SystemExit(__doc__)
