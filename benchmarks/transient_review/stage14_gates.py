"""STAGE 14 -- PCNR for a BJT: multi-junction devices, and charge on limited unknowns.

See doc/transient_work_plan.md, stage 14, for what each mode gates.

    matrix     14-2a  BJT.pcnr_didv against an analytic Ebers-Moll derivative
    dc         14-2b/c/d  the CE stage: classic vs PCNR, call counts
    jax        14-2e  the CE stage under PCNR on JAX vs numeric (repo .venv)
    tran       14-3a/c  BJT switch: classic vs PCNR steps and waveform, call counts
    jacobian   14-3b  controller Jacobian vs a finite difference of the residual
    coupled    14-3d  coupled path against a tight standard-path reference

Run from the repo root:
    PYTHONPATH=. MPLBACKEND=Agg python3 -u benchmarks/transient_review/stage14_gates.py matrix
"""
import sys
import warnings

import numpy as np

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import pcnr as P
from pycircuit.circuit.semiconductors import BJT

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import stage13_10_11_gates as Gx   # noqa: E402  ce_stage, switching_stage

VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron
PARAMS = {'IS': 1e-14, 'BF': 100.0, 'BR': 1.0, 'VA': 100.0}


def ebers_moll_didv(vbe, vbc, p=PARAMS):
    """Analytic d(ic, ib, ie)/d(vbe, vbc) for BJT.eval_i_pure, Early factor included."""
    IS, BF, BR, VA = p['IS'], p['BF'], p['BR'], p['VA']
    ebe, ebc = np.exp(vbe / VT), np.exp(vbc / VT)
    i_be, i_bc = IS * (ebe - 1.0), IS * (ebc - 1.0)
    g_be, g_bc = IS * ebe / VT, IS * ebc / VT
    early = 1.0 + (vbe - vbc) / VA
    dicc_be = g_be * early + (i_be - i_bc) / VA
    dicc_bc = -g_bc * early - (i_be - i_bc) / VA
    return np.array([[dicc_be, dicc_bc - g_bc / BR],
                     [g_be / BF, g_bc / BR],
                     [-dicc_be - g_be / BF, -dicc_bc]])


def matrix():
    """Three instruments, because the first one alone is wrong in cut-off.

    1. Against the analytic derivative, relative to its largest entry.  In cut-off
       every conductance is ~1e-29 S and the numeric toolkit's central difference
       returns exactly 0 -- `IS*(exp(-38.7) - 1)` cannot resolve the exponential
       next to 1 -- so this reads 1.0 there while being 6e-30 S absolute.
    2. `D == jacobian(eval_i_pure) @ M` exactly: the adapter adds nothing to the
       toolkit's own differentiation.
    3. With EXACT differentiation (JAX autodiff, x64) against the analytic form:
       the adapter's chain rule, with no difference quotient in the way.
    """
    try:
        import jax
        jax.config.update('jax_enable_x64', True)
        from pycircuit.circuit.toolkit import jaxtoolkit
    except Exception:  # noqa: BLE001
        jaxtoolkit = None
    M = BJT._pcnr_map()
    for label, vbe, vbc in (('forward active', 0.70, -4.30), ('saturated', 0.75, 0.60),
                            ('cut off', -1.0, -5.0), ('reverse active', -3.0, 0.70)):
        v = np.array([vbe, vbc])
        A = ebers_moll_didv(vbe, vbc)
        D = np.asarray(BJT.pcnr_didv(v, PARAMS, defaultepar, numeric), dtype=float)
        GM = np.asarray(numeric.jacobian(BJT.eval_i_pure, M @ v, PARAMS, defaultepar),
                        dtype=float) @ M
        rel = float(np.max(np.abs(D - A)) / np.max(np.abs(A)))
        xb = 1.234                                       # arbitrary common offset
        xt = np.array([xb - vbc, xb, xb - vbe])          # c, b, e
        i_ref = np.asarray(BJT.eval_i_pure(xt, PARAMS, defaultepar, numeric), dtype=float)
        i_got = np.asarray(BJT.pcnr_i(v, PARAMS, defaultepar, numeric), dtype=float)
        irel = float(np.max(np.abs(i_got - i_ref)) / max(np.max(np.abs(i_ref)), 1e-300))
        line = ('%-15s | numeric vs analytic: rel %.2e, abs %.2e S (max|A| %.2e) '
                '| D == G@M %s | pcnr_i vs eval_i_pure %.1e'
                % (label, rel, float(np.max(np.abs(D - A))), float(np.max(np.abs(A))),
                   bool(np.array_equal(D, GM)), irel))
        if jaxtoolkit is not None:
            DJ = np.asarray(BJT.pcnr_didv(v, PARAMS, defaultepar, jaxtoolkit), dtype=float)
            line += ' | autodiff vs analytic rel %.2e' % (
                float(np.max(np.abs(DJ - A)) / np.max(np.abs(A))))
        print(line)


class Count:
    def __init__(self, cls, names):
        self.cls, self.names = cls, names
        self.n = dict.fromkeys(names, 0)
        self.real = {k: getattr(cls, k) for k in names}

    def __enter__(self):
        for k in self.names:
            real = self.real[k]

            def spy(*a, _k=k, _real=real, **kw):
                self.n[_k] += 1
                return _real(*a, **kw)
            setattr(self.cls, k, spy)
        return self

    def __exit__(self, *exc):
        for k in self.names:
            setattr(self.cls, k, self.real[k])
        return False


def dc():
    for rb in (10.0, 100.0, 1000.0):
        xc, its_c = Gx.dc_classic(lambda r=rb: Gx.ce_stage(r))
        c = Gx.ce_stage(rb)
        with Count(BJT, ('limit', 'G', 'i')) as cnt:
            xp, vp, its_p = P.solve_dc(c, gnd, maxiter=500)
        xp = np.asarray(xp, float)
        f = (np.asarray(c.i(xp, defaultepar), float)
             + np.asarray(c.u(0, analysis='dc', epar=defaultepar), float))
        f = np.delete(f, c.get_node_index(gnd))
        rel = float(np.max(np.abs(xp - xc)) / np.max(np.abs(xc)))
        vb, vcol = float(xp[c.get_node_index('b')]), float(xp[c.get_node_index('c')])
        print('CE rb=%-6g classic %2d its | PCNR %2d its | max rel diff %.3e | residual %.3e '
              '| v(b) %.6f v(c) %.6f | BJT calls during PCNR %s'
              % (rb, its_c[0], its_p, rel, float(np.max(np.abs(f))), vb, vcol, cnt.n))
    hits = {'aug': 0}
    real = P.augmented_system

    def spy(*a, **k):
        hits['aug'] += 1
        return real(*a, **k)
    P.augmented_system = spy
    try:
        with Count(BJT, ('limit',)) as cnt, warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r = DC(Gx.ce_stage(100.0), toolkit=numeric, pcnr=True).solve()
    finally:
        P.augmented_system = real
    print('DC(pcnr=True) on CE rb=100: augmented_system %d, BJT.limit %d, v(b) %.6f'
          % (hits['aug'], cnt.n['limit'], float(r.v('b'))))


def jax():
    for rb in (10.0, 100.0, 1000.0):
        cj = Gx_jax(lambda r=rb: Gx.ce_stage(r))
        xj, _vj, itj = P.solve_dc(cj, gnd, maxiter=500)
        xn, _vn, itn = P.solve_dc(Gx.ce_stage(rb), gnd, maxiter=500)
        rel = float(np.max(np.abs(np.asarray(xj, float) - xn)) / np.max(np.abs(xn)))
        print('CE rb=%-6g JAX %2d its, numeric %2d its, max rel diff %.3e' % (rb, itj, itn, rel))


def Gx_jax(builder):
    from pycircuit.circuit import circuit as cm
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = cm.default_toolkit
    cm.default_toolkit = jaxtoolkit
    try:
        c = builder()
        c.toolkit = jaxtoolkit
        for e in c.elements.values():
            e.toolkit = jaxtoolkit
        c.update_iparv()
        return c
    finally:
        cm.default_toolkit = saved


def _tran(pcnr, coupled=False, reltol=1e-5, build=None):
    t = Transient((build or Gx.switching_stage)(), toolkit=numeric, reltol=reltol, pcnr=pcnr)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = t.solve(tend=1.6e-6, timestep=2.5e-9, coupled_lte=coupled)
    w = res.v('d' if build is jfet_stage else 'c')
    s = t.statistics
    return (np.asarray(w.x, float).ravel(), np.asarray(w.y, float).ravel(),
            (s.accepted_steps, s.rejected_steps, s.newton_iterations))


def tran():
    """14-3a (13-6's standard) and 14-3c."""
    tc, yc, nc = _tran(False)
    names = ('limit', 'G', 'i', 'C', 'q')
    calls = dict.fromkeys(names, 0)
    real_m = {k: getattr(BJT, k) for k in names}
    real_step = Transient._solve_timestep_pcnr
    armed = {'on': False}

    def step(self, *a, **k):
        armed['on'] = True
        return real_step(self, *a, **k)

    def counted(name):
        def spy(*a, **k):
            if armed['on']:
                calls[name] += 1
            return real_m[name](*a, **k)
        return spy
    Transient._solve_timestep_pcnr = step
    for k in names:
        setattr(BJT, k, counted(k))
    try:
        tp, yp, np_ = _tran(True)
    finally:
        Transient._solve_timestep_pcnr = real_step
        for k in names:
            setattr(BJT, k, real_m[k])
    print('BJT switch, standard path: classic accepted/rejected/newton %s | PCNR %s' % (nc, np_))
    if tc.shape == tp.shape and np.array_equal(tc, tp):
        d = float(np.max(np.abs(yc - yp)) / np.max(np.abs(yc)))
        print('  identical time grids; waveform max rel diff %.3e' % d)
    else:
        d = float(np.max(np.abs(np.interp(tc, tp, yp) - yc)) / np.max(np.abs(yc)))
        print('  grids DIFFER (%d vs %d points); interpolated max rel diff %.3e'
              % (tc.size, tp.size, d))
    print('  BJT calls from the first PCNR timestep on: %s' % calls)


def jacobian(samples=40):
    """14-3b: the controller Jacobian against a central difference of the full
    transient residual F(x) = i(x) + iq(x) + u(t), with every device evaluated
    at the node voltages -- which at convergence is where v_lim is."""
    real = Transient._solve_timestep_pcnr
    worst = {'rel': 0.0, 'n': 0, 'seen': 0}

    def spy(self, x0, t, *a, **k):
        x, feval, J, f = real(self, x0, t, *a, **k)
        worst['seen'] += 1
        if worst['seen'] % 7 or worst['n'] >= samples:
            return x, feval, J, f
        cir, ep = self.cir, self.epar
        h_last = getattr(self, '_dt_last', self._dt)

        def F(xx):
            q = np.asarray(cir.q(xx, ep), float)
            iq, _g = self.active_integrator.compute_derivatives(
                q_curr=q, C_curr=1.0, h_curr=self._dt, q_last=self._qlast,
                iq_last=self._iqlast, h_last=h_last,
                is_first_step=getattr(self, '_is_first_step', False),
                toolkit=self.toolkit)
            return (np.asarray(cir.i(xx, ep), float) + np.asarray(iq, float)
                    + np.asarray(cir.u(t, ep, analysis=self.par.analysis), float))
        xx = np.asarray(x, float)
        n = len(xx)
        FD = np.zeros((n, n))
        for j in range(n):
            dl = 1e-7 * max(1.0, abs(xx[j]))
            xp, xm = xx.copy(), xx.copy()
            xp[j] += dl
            xm[j] -= dl
            FD[:, j] = (F(xp) - F(xm)) / (2 * dl)
        keep = [i for i in range(n) if i != self.irefnode]
        Jr = np.asarray(J, float)[np.ix_(keep, keep)]
        FDr = FD[np.ix_(keep, keep)]
        rel = float(np.max(np.abs(Jr - FDr)) / np.max(np.abs(FDr)))
        worst['rel'] = max(worst['rel'], rel)
        worst['n'] += 1
        return x, feval, J, f
    Transient._solve_timestep_pcnr = spy
    try:
        _tran(True)
    finally:
        Transient._solve_timestep_pcnr = real
    print('controller Jacobian vs central difference of the residual: worst rel %.3e over '
          '%d converged time points' % (worst['rel'], worst['n']))


def coupled():
    """14-3d: both coupled runs against a tight standard-path reference."""
    tr, yr, nr = _tran(False, reltol=1e-8)
    print('reference (standard, reltol 1e-8): %d accepted steps' % nr[0])
    for label, pc in (('classic', False), ('PCNR', True)):
        try:
            t, y, n = _tran(pc, coupled=True)
        except Exception as e:  # noqa: BLE001
            print('coupled %-7s FAILED %s: %s' % (label, type(e).__name__, str(e)[:100]))
            continue
        e = np.abs(y - np.interp(t, tr, yr))
        print('coupled %-7s steps %s | max err %.4e V | median err %.3e V | last t - tend %+.1e'
              % (label, n, float(e.max()), float(np.median(e)), t[-1] - 1.6e-6))


## ------------------------------------------------------------ 14-5: JFET

JP = {'VTO': -2.0, 'BETA': 1e-3, 'LAMBDA': 0.01, 'IS': 1e-14}


def jfet_sat_didv(vgs, vgd, p=JP):
    """Analytic d(id, ig, is)/d(vgs, vgd) in forward saturation (vds > 0)."""
    VTO, B, L, IS = p['VTO'], p['BETA'], p['LAMBDA'], p['IS']
    vds = vgs - vgd
    ov = vgs - VTO
    did_gs = 2 * B * ov * (1 + L * vds) + B * ov ** 2 * L
    did_gd = -B * ov ** 2 * L
    ggs, ggd = IS * np.exp(vgs / VT) / VT, IS * np.exp(vgd / VT) / VT
    return np.array([[did_gs, did_gd - ggd],
                     [ggs, ggd],
                     [-did_gs - ggs, -did_gd]])


def jfet_cs(vg, vdd=10.0, rd=2e3, rg=1e6):
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit.elements import R, VS
    from pycircuit.circuit.semiconductors import JFET
    c = SubCircuit(toolkit=numeric)
    c['VDD'] = VS('vdd', gnd, v=vdd)
    c['VG'] = VS('vg', gnd, v=vg)
    c['RD'] = R('vdd', 'd', r=rd)
    c['RG'] = R('vg', 'g', r=rg)
    c['J1'] = JFET('d', 'g', gnd, **JP)
    return c


JFET_CASES = [('saturation', lambda: jfet_cs(-1.0)),
              ('triode', lambda: jfet_cs(0.0, vdd=1.0, rd=100.0)),
              ('gate forward', lambda: jfet_cs(2.0, rg=1e3)),
              ('cut off', lambda: jfet_cs(-3.0))]


def jfet_stage():
    """A JFET stage with a capacitive load and a pulsed gate -- the device has no
    charge model, so the LTE comes from the load."""
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit.elements import R, C, VS, VPulse
    from pycircuit.circuit.semiconductors import JFET
    per = 1e-6
    c = SubCircuit(toolkit=numeric)
    c['VDD'] = VS('vdd', gnd, v=10.0)
    c['VG'] = VPulse('vg', gnd, v1=-3.0, v2=1.0, td=per / 10, tr=per / 100,
                     tf=per / 100, pw=per / 2, per=per)
    c['RG'] = R('vg', 'g', r=1e3)
    c['RD'] = R('vdd', 'd', r=2e3)
    c['CL'] = C('d', gnd, c=1e-11)
    c['J1'] = JFET('d', 'g', gnd, **JP)
    return c


def jfet():
    from pycircuit.circuit.semiconductors import JFET
    try:
        import jax  # noqa: F401
        from pycircuit.circuit.toolkit import jaxtoolkit
    except Exception:  # noqa: BLE001
        jaxtoolkit = None
    M = JFET._pcnr_map()
    print('14-5a')
    for label, vgs, vgd in (('saturation', -1.0, -8.5), ('triode', 0.0, -0.05),
                            ('cut off', -3.0, -13.0), ('gate forward', 0.7, -7.0)):
        v = np.array([vgs, vgd])
        D = np.asarray(JFET.pcnr_didv(v, JP, defaultepar, numeric), dtype=float)
        GM = np.asarray(numeric.jacobian(JFET.eval_i_pure, M @ v, JP, defaultepar),
                        dtype=float) @ M
        line = '  %-13s D == G@M %s' % (label, bool(np.array_equal(D, GM)))
        if jaxtoolkit is not None:
            DJ = np.asarray(JFET.pcnr_didv(v, JP, defaultepar, jaxtoolkit), dtype=float)
            line += ' | autodiff vs numeric rel %.2e' % (
                float(np.max(np.abs(DJ - D)) / np.max(np.abs(DJ))))
            if label == 'saturation':
                A = jfet_sat_didv(vgs, vgd)
                line += ' | autodiff vs analytic rel %.2e' % (
                    float(np.max(np.abs(DJ - A)) / np.max(np.abs(A))))
        print(line)
    print('14-5b')
    for label, b in JFET_CASES:
        xc, its_c = Gx.dc_classic(b)
        c = b()
        with Count(JFET, ('limit', 'G', 'i')) as cnt:
            xp, _vp, its_p = P.solve_dc(c, gnd, maxiter=500)
        rel = float(np.max(np.abs(np.asarray(xp, float) - xc)) / np.max(np.abs(xc)))
        print('  %-13s classic %2d its | PCNR %2d its | max rel diff %.3e | JFET calls %s'
              % (label, its_c[0], its_p, rel, cnt.n))


def jfet_tran():
    """14-5c."""
    from pycircuit.circuit.semiconductors import JFET
    out = {}
    for pc in (False, True):
        names = ('limit', 'G', 'i')
        calls = dict.fromkeys(names, 0)
        real_m = {k: getattr(JFET, k) for k in names}
        real_step = Transient._solve_timestep_pcnr
        armed = {'on': False}

        def step(self, *a, **k):
            armed['on'] = True
            return real_step(self, *a, **k)

        def counted(name):
            def spy(*a, **k):
                if armed['on']:
                    calls[name] += 1
                return real_m[name](*a, **k)
            return spy
        Transient._solve_timestep_pcnr = step
        for k in names:
            setattr(JFET, k, counted(k))
        try:
            t, y, n = _tran(pc, build=jfet_stage)
        finally:
            Transient._solve_timestep_pcnr = real_step
            for k in names:
                setattr(JFET, k, real_m[k])
        out[pc] = (t, y, n, dict(calls), armed['on'])
    (tc, yc, nc, _c, _a), (tp, yp, np_, cp, ap) = out[False], out[True]
    same = tc.shape == tp.shape
    d = float(np.max(np.abs(yc - yp)) / np.max(np.abs(yc))) if same else float('nan')
    dt = float(np.max(np.abs(tc - tp))) if same else float('nan')
    print('JFET stage: classic %s | PCNR %s (PCNR path ran: %s) | waveform max rel diff %.3e '
          '| max |dt| %.2e | JFET calls after the op %s' % (nc, np_, ap, d, dt, cp))


## ------------------------------------------------------------ 14-6

TRAN_CASES = [('rect', Gx.G89.rect, 'b', 2e-3, 1e-5), ('mains', Gx.G89.mains, 2, 40e-3, 1e-4),
              ('clipper', Gx.clipper, 2, 2e-3, 1e-5),
              ('bjt_switch', Gx.switching_stage, 'c', 1.6e-6, 2.5e-9),
              ('jfet_stage', jfet_stage, 'd', 1.6e-6, 2.5e-9)]


def _run_counted(build, node, tend, ts, pcnr):
    """Waveform and counts, plus classic limiter calls over the WHOLE run and PCNR
    assemblies before the first timestep (i.e. inside the operating point)."""
    from pycircuit.circuit.elements import Diode
    from pycircuit.circuit.semiconductors import JFET
    lim = {'n': 0}
    pre = {'aug': 0, 'stepping': False}
    reals = [(cls, cls.limit) for cls in (Diode, BJT, JFET)]
    real_aug, real_step = P.augmented_system, Transient._solve_timestep_pcnr
    real_ts = Transient.solve_timestep

    def mk(real):
        def spy(*a, **k):
            lim['n'] += 1
            return real(*a, **k)
        return spy

    def aug(*a, **k):
        if not pre['stepping']:
            pre['aug'] += 1
        return real_aug(*a, **k)

    def ts_spy(self, *a, **k):
        pre['stepping'] = True
        return real_ts(self, *a, **k)
    for cls, real in reals:
        cls.limit = mk(real)
    P.augmented_system, Transient.solve_timestep = aug, ts_spy
    try:
        t = Transient(build(), toolkit=numeric, reltol=1e-5, pcnr=pcnr)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = t.solve(tend=tend, timestep=ts)
    finally:
        for cls, real in reals:
            cls.limit = real
        P.augmented_system, Transient.solve_timestep = real_aug, real_ts
    w = res.v(node)
    st = t.statistics
    return (np.concatenate([np.asarray(w.x, float).ravel(), np.asarray(w.y, float).ravel()]),
            [st.accepted_steps, st.rejected_steps, st.newton_iterations, lim['n'], pre['aug']])


def snapshot(path):
    data = {}
    for name, b, node, tend, ts in TRAN_CASES:
        for pc in (False, True):
            key = '%s_%s' % (name, 'pcnr' if pc else 'classic')
            vec, counts = _run_counted(b, node, tend, ts, pc)
            data[key], data[key + '__counts'] = vec, np.array(counts, dtype=float)
            print('%-20s accepted/rejected/newton/limiter calls/op assemblies %s'
                  % (key, counts), flush=True)
    np.savez(path, **data)


def compare(pa, pb):
    a, b = np.load(pa), np.load(pb)
    for name in sorted(k for k in a.files if not k.endswith('__counts')):
        x, y = a[name], b[name]
        ca, cb = a[name + '__counts'].astype(int).tolist(), b[name + '__counts'].astype(int).tolist()
        if x.shape != y.shape:
            v = 'shape %s -> %s' % (x.shape, y.shape)
        elif np.array_equal(x, y):
            v = 'bit-identical'
        else:
            v = 'max rel diff %.3e' % (np.max(np.abs(x - y)) / np.max(np.abs(x)))
        print('%-20s %-24s %s -> %s' % (name, v, ca, cb))


if __name__ == '__main__':
    {'matrix': matrix, 'dc': dc, 'jax': jax, 'tran': tran, 'jacobian': jacobian,
     'coupled': coupled, 'jfet': jfet, 'jfet_tran': jfet_tran,
     'snapshot': lambda: snapshot(sys.argv[2]),
     'compare': lambda: compare(sys.argv[2], sys.argv[3])}[sys.argv[1]]()
