"""GATES 13-10 and 13-11 -- `_pnjlim`'s missing condition, and skip on the batched path.

13-10: `_pnjlim` compresses every step above `vc`.  SPICE limits only when the
step also exceeds 2 VT (ngspice-47 `devsup.c:54`:
`if((vnew > vcrit) && (fabs(vnew - vold) > (vt + vt)))`).  Gate 13-9c traced a
3.2e-14 waveform difference to exactly this.
13-11: a `skip` of a class the JAX toolkit stamps in one `vmap` is refused.
That refusal also made PCNR unusable on a JAX-toolkit circuit.

See doc/transient_work_plan.md gates 13-10 and 13-11 for what each measurement
gates.

Modes:

    table                   `_pnjlim` on a grid of (vnew, vold)
    snapshot OUT.npz        DC and transient cases, classic and PCNR
    compare  A.npz B.npz    case by case
    mismatch                converged classic solves with `_vlim != v`
    jax                     PCNR DC on JAX-toolkit circuits vs numeric

Run from the repo root (`jax` needs the repo's .venv):
    PYTHONPATH=. MPLBACKEND=Agg python3 -u benchmarks/transient_review/stage13_10_11_gates.py table
"""
import sys
import warnings

import numpy as np

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit, defaultepar
from pycircuit.circuit.elements import R, C, VS, VSin, VPulse, Diode
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import pcnr as P
from pycircuit.circuit._limiting import _pnjlim

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import stage13_pcnr as S13         # noqa: E402  the six DC circuits of 13-4/13-5
import stage13_8_9_gates as G89    # noqa: E402  rect, mains

VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron


## ------------------------------------------------------------------ table

def table():
    IS = 1e-14
    vc = VT * np.log(VT / (IS * 1.414213562))
    print('vc = %.6f V, 2 VT = %.6f V' % (vc, 2 * VT))
    print('%10s %10s %12s %14s %10s' % ('vold', 'vnew', 'step/VT', '_pnjlim', 'limited'))
    for vold in (0.0, 0.5, vc - 0.01, vc + 0.01):
        for dv in (0.001, 0.02, 0.05, 0.06, 0.2, 2.0):
            vnew = vold + dv
            out = _pnjlim(vnew, vold, VT, IS, numeric)
            print('%10.4f %10.4f %12.3f %14.9f %10s'
                  % (vold, vnew, dv / VT, out, 'yes' if out != vnew else ''))


## ------------------------------------------------------------------ cases

def ce_stage(rb):
    from pycircuit.circuit.semiconductors import BJT
    cir = SubCircuit(toolkit=numeric)
    cir['VCC'] = VS('vcc', gnd, v=5.0)
    cir['VIN'] = VS('vin', gnd, v=0.8)
    cir['RB'] = R('vin', 'b', r=rb)
    cir['RC'] = R('vcc', 'c', r=1e3)
    cir['Q1'] = BJT('c', 'b', gnd)
    return cir


def switching_stage(per=1e-6):
    from pycircuit.circuit.semiconductors import BJT
    cir = SubCircuit(toolkit=numeric)
    cir['VCC'] = VS('vcc', gnd, v=5.0)
    cir['VB'] = VPulse('vb', gnd, v1=0.0, v2=0.9, td=per / 10, tr=per / 1000,
                       tf=per / 1000, pw=per / 2, per=per)
    cir['RB'] = R('vb', 'b', r=1e3)
    cir['RC'] = R('vcc', 'c', r=1e3)
    cir['Q1'] = BJT('c', 'b', gnd)
    return cir


def bridge():
    c = SubCircuit()
    c['VS'] = VSin(1, 2, va=10, freq=50)
    c['R_src'] = R(2, gnd, r=1e-3)
    c['D1'] = Diode(1, 3)
    c['D2'] = Diode(gnd, 1)
    c['D3'] = Diode(2, 3)
    c['D4'] = Diode(gnd, 2)
    c['C1'] = C(3, gnd, c=10e-6)
    c['R1'] = R(3, gnd, r=100)
    return c


def pump():
    c = SubCircuit()
    c['VP'] = VPulse(1, gnd, v1=-5, v2=5, tr=1e-6, tf=1e-6, pw=5e-6, per=10e-6)
    c['C1'] = C(1, 2, c=1e-6)
    c['D1'] = Diode(gnd, 2)
    c['D2'] = Diode(2, 3)
    c['C2'] = C(3, gnd, c=1e-6)
    return c


def clipper():
    c = SubCircuit()
    c['VS'] = VSin(1, gnd, va=100, freq=1e3)
    c['R1'] = R(1, 2, r=10)
    c['D1'] = Diode(2, gnd)
    c['D2'] = Diode(gnd, 2)
    return c


class _NewtonCount:
    """Counts Newton iterations of every StandardNewton solve, as 13-4 did."""

    def __init__(self):
        import pycircuit.circuit.nrsolver as nrs
        self.nrs, self.n = nrs, 0
        self.real = nrs.StandardNewton.solve_system

    def __enter__(self):
        real, me = self.real, self

        def patched(solver, x0, func, *a, **k):
            def spy(x, *aa, **kk):
                me.n += 1
                return func(x, *aa, **kk)
            return real(solver, x0, spy, *a, **k)
        self.nrs.StandardNewton.solve_system = patched
        return self

    def __exit__(self, *exc):
        self.nrs.StandardNewton.solve_system = self.real
        return False


def dc_classic(build):
    c = build()
    with _NewtonCount() as nc, warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = DC(c, toolkit=numeric).solve()
    x = np.asarray(res.x, float).ravel()
    ## the unlimited residual at the answer, as gate 5's check does
    f = (np.asarray(c.i(x, defaultepar), float)
         + np.asarray(c.u(0, analysis='dc', epar=defaultepar), float))
    f = np.delete(f, c.get_node_index(gnd))
    return x, [nc.n, float(np.max(np.abs(f)))]


def dc_pcnr(build):
    x, v, its = P.solve_dc(build(), gnd, maxiter=500)
    return np.concatenate([np.asarray(x, float), np.asarray(v, float)]), [its, 0.0]


def tran(build, node, tend, ts, pcnr=False, coupled=False, uic=False):
    t = Transient(build(), toolkit=numeric, reltol=1e-5, pcnr=pcnr, uic=uic)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = t.solve(tend=tend, timestep=ts, coupled_lte=coupled)
    w = res.v(node)
    s = t.statistics
    return (np.concatenate([np.asarray(w.x, float).ravel(), np.asarray(w.y, float).ravel()]),
            [s.accepted_steps, s.rejected_steps, s.newton_iterations])


def cases():
    out = []
    for name, build in S13.CASES:
        key = name.replace(' ', '_').replace('(', '').replace(')', '')
        out.append(('dc_classic_' + key, lambda b=build: dc_classic(b)))
        out.append(('dc_pcnr_' + key, lambda b=build: dc_pcnr(b)))
    for rb in (10.0, 100.0, 1000.0):
        out.append(('dc_classic_ce_rb%g' % rb, lambda r=rb: dc_classic(lambda: ce_stage(r))))
    out += [
        ('tran_classic_rect', lambda: tran(G89.rect, 'b', 2e-3, 1e-5)),
        ('tran_classic_mains', lambda: tran(G89.mains, 2, 40e-3, 1e-4)),
        ('tran_classic_bridge', lambda: tran(bridge, 3, 0.04, 1e-4)),
        ('tran_classic_pump', lambda: tran(pump, 3, 50e-6, 1e-7, uic=True)),
        ('tran_classic_clipper', lambda: tran(clipper, 2, 2e-3, 1e-5)),
        ('tran_classic_bjt_switch', lambda: tran(switching_stage, 'c', 1.6e-6, 2.5e-9)),
        ('coupled_classic_rect', lambda: tran(G89.rect, 'b', 2e-3, 1e-5, coupled=True)),
        ('tran_pcnr_rect', lambda: tran(G89.rect, 'b', 2e-3, 1e-5, pcnr=True)),
        ('tran_pcnr_mains', lambda: tran(G89.mains, 2, 40e-3, 1e-4, pcnr=True)),
        ('tran_pcnr_clipper', lambda: tran(clipper, 2, 2e-3, 1e-5, pcnr=True)),
    ]
    return out


def snapshot(path):
    data = {}
    for name, fn in cases():
        try:
            vec, counts = fn()
        except Exception as e:  # noqa: BLE001  -- a failure IS a result here
            print('%-30s FAILED %s: %s' % (name, type(e).__name__, str(e)[:90]), flush=True)
            continue
        data[name] = vec
        data[name + '__counts'] = np.array(counts, dtype=float)
        print('%-30s counts %s' % (name, counts), flush=True)
    np.savez(path, **data)


def compare(pa, pb):
    a, b = np.load(pa), np.load(pb)
    names = sorted(set(k for k in a.files + b.files if not k.endswith('__counts')))
    for name in names:
        if name not in a.files or name not in b.files:
            print('%-30s present only in %s' % (name, 'A' if name in a.files else 'B'))
            continue
        x, y = a[name], b[name]
        ca, cb = a[name + '__counts'].tolist(), b[name + '__counts'].tolist()
        if x.shape != y.shape:
            verdict = 'shape %s -> %s' % (x.shape, y.shape)
        elif np.array_equal(x, y):
            verdict = 'bit-identical'
        else:
            verdict = 'max rel diff %.3e' % (np.max(np.abs(x - y))
                                             / max(np.max(np.abs(x)), 1e-300))
        print('%-30s %-24s %s -> %s' % (name, verdict,
                                        ['%g' % v for v in ca], ['%g' % v for v in cb]))


## ------------------------------------------------------------------ 13-9c recheck

def mismatch():
    """Gate 13-9c traced its 3.2e-14 to converged solves leaving `_vlim != v`."""
    for label, build, name in (('mains', G89.mains, 'D1'), ('rect', G89.rect, 'D')):
        c = build()
        d = c[name]
        rows = c.elementnodemap[name]
        t = Transient(c, toolkit=numeric, reltol=1e-5)
        real = Transient._newton
        st = {'n': 0, 'bad': 0, 'maxd': 0.0}

        def spy(self, func, x0):
            x = real(self, func, x0)
            v = float(x[rows[0]] - x[rows[1]])
            st['n'] += 1
            if float(d._vlim) != v:
                st['bad'] += 1
                st['maxd'] = max(st['maxd'], abs(float(d._vlim) - v))
            return x

        Transient._newton = spy
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                t.solve(tend=40e-3 if label == 'mains' else 2e-3,
                        timestep=1e-4 if label == 'mains' else 1e-5)
        finally:
            Transient._newton = real
        print('%-5s converged solves %d, with _vlim != v: %d, max |_vlim - v| %.3e V'
              % (label, st['n'], st['bad'], st['maxd']))


## ------------------------------------------------------------------ 13-11

def _jax_build(builder):
    from pycircuit.circuit import circuit as cm
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = cm.default_toolkit
    cm.default_toolkit = jaxtoolkit
    try:
        c = builder()
        c.toolkit = jaxtoolkit
        c.update_iparv()
        return c
    finally:
        cm.default_toolkit = saved


def jax_probe():
    """PCNR DC on a JAX-toolkit circuit, against the same circuit on numeric."""
    from pycircuit.circuit.toolkit import jaxtoolkit
    sizes = []
    real_gen = jaxtoolkit.generate_batched_eval

    def spy_gen(cls, method='i'):
        f = real_gen(cls, method)

        def g(X, p, epar):
            sizes.append((cls.__name__, int(np.shape(X)[0])))
            return f(X, p, epar)
        return g

    jaxtoolkit.generate_batched_eval = spy_gen
    try:
        for name, build in S13.CASES:
            sizes.clear()
            try:
                cj = _jax_build(build)
                grouped = sorted(k.__name__ for k in (cj._eval_groups or {}))
                xj, vj, itj = P.solve_dc(cj, gnd, maxiter=500)
                xn, vn, itn = P.solve_dc(build(), gnd, maxiter=500)
                d = float(np.max(np.abs(np.asarray(xj, float) - xn))
                          / max(np.max(np.abs(xn)), 1e-300))
                diode_batches = sorted(set(n for c_, n in sizes if c_ == 'Diode'))
                print('%-20s groups %s | iters jax %d numeric %d | max rel diff %.3e '
                      '| Diode batch sizes seen %s'
                      % (name, grouped, itj, itn, d, diode_batches or 'none'))
            except Exception as e:  # noqa: BLE001
                print('%-20s %s: %s' % (name, type(e).__name__, str(e)[:110]))
    finally:
        jaxtoolkit.generate_batched_eval = real_gen


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'table':
        table()
    elif mode == 'snapshot':
        snapshot(sys.argv[2])
    elif mode == 'compare':
        compare(sys.argv[2], sys.argv[3])
    elif mode == 'mismatch':
        mismatch()
    elif mode == 'jax':
        jax_probe()
    else:
        raise SystemExit(__doc__)
