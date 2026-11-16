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


if __name__ == '__main__':
    {'matrix': matrix, 'dc': dc, 'jax': jax}[sys.argv[1]]()
