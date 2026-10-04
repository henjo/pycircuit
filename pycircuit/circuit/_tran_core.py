"""THE EVALUATE CORE: a transient's residual, Jacobian and companion
evaluated in one C call (2026-10-03; speed round 7, commit 1).

A Newton iterate of the numeric transient spends most of its time in the
Python around C that already exists: the four circuit passes (`G`, `C`,
`i`, `q`), each a plan with its batches, templates and bincount, then the
companion term and the residual in numpy.  Measured warm on a 20-MosLevel1
chain: the four passes 45 us, the residual 48.5, of which the C work is
12.7.  Here ONE C call does the passes the caller asked for -- the batches
through the kernels' own pointers, the constant vector groups through
numpy's OWN `cblas_dgemv` (the address of the routine numpy loaded, so the
same kernel runs: measured 1000 of 1000 products equal to `np.dot` and to
the plan's `matmul`), the templates, the bincounts from +0.0 in array order
-- then the companion in the integrator's own operation order, `f = (i +
iq) + u` and `J = G + Geq`.  Bit for bit the Python path: measured on every
trial of the stage-0 prototype (pss_log 2026-10-03), 20-MosLevel1 357 ->
311 us a step in the run.  The reference-row reduction stays Python's.

THREE SITES, ONE ENTRY (`evaluate`, `want` in 'fj' | 'j' | 'f'):
`_residual_and_jacobian` (G, C, i, q; `(f, J)`); `jacobian_only` at the
converged point (G, C, q; `(None, J)`); the chord's `residual_only` (i, q;
the conductance from the seed's held `C`; `f`).  Each leaves the state the
Python path leaves: `_C_cache`, `_q_cache`, the memo record `_C_at_state`
makes, and `get_diff`'s `_iq`, `_Geq`, `_Cmat`, `_companion_coeffs`,
`_effective_method`, `active_integrator`.  Python still decides the active
integrator (`check_order_drop`) and computes its scalars (`bdf2_alphas`,
`theta_at`), and still asks `_C_lookup` first: on a hit the C pass is
skipped, as `_C_at_state` skips it.

THE FORMULA ENUM, each as `_lte_kernels` writes it: gear2 `(a0*q + a1*q1)
+ a2*q2`, `a0*C`; euler `(q - q1)/h`, `C/h`; trap `2.0*(q - q1)/h - iq1`,
`2.0*C/h`; theta `th = theta*h`, `(q - q1)/th - ((1.0 - theta)/theta)*iq1`,
`C/th`.  The stage methods, Radau, PCNR, DC and every toolkit but the
numeric one are not served.

A DECLINE BEFORE ANY STATE IS TOUCHED (None; the Python path runs): the
switch off, a toolkit that is not numeric, a state that is not one finite
float64 vector of the circuit's size, a circuit with an element neither
constant nor C-bound (a hand-written limiter, a Python hdl class, a nested
circuit: the core does not exist for it until its plan rebuilds), a batch
not ready (a detached class, a patched method, an instance shadow of any
of the four methods, a pack that failed, a temperature that is not one
number), a history row of the wrong shape, an integrator outside the enum.
The core's tables are keyed by the stamp plan (rebuilt when it is), the
pack mirrors are the core's own (`_hdl_batch`'s rule: mirror the tuple the
element holds, never repack one), the readiness check runs once per call
over the unique elements rather than once per pass.

`CORE` (env `PYCIRCUIT_TRAN_CORE=0`) is read on every call; `STATUS` says
why the core is off where it is.  Pinned against its Python reference in
`tests/test_pinned_pairs.py`.  History: `doc/transient_history.md`,
`_tran_core.py`; `doc/hdl_roadmap_260824.md` sec. 65.
"""
import ctypes
import glob
import os

import numpy as np

from pycircuit.circuit import _paths

## the counters (`_paths`): each decline by its reason, each served call
_PC = _paths.COUNTS
_no = _paths.no
_K = {w: {r: f'core.{w}:{r}' for r in (
    'off', 'toolkit', 'shadow_cir', 'shadow_tr', 'unservable', 'x', 'formula',
    'history', 'cmat', 'clookup', 'ready', 'u', 'served')} for w in ('fj', 'j', 'f', 'c')}

CORE = os.environ.get('PYCIRCUIT_TRAN_CORE', '1') != '0'
#: where nothing reads a step's Jacobian (a fixed multistep run:
#: `transient._SteppingLoop.j_unread`), the converged point is evaluated
#: without G -- `want` 'c' here, the C Newton's fold likewise (speed round
#: 9, C0); env `PYCIRCUIT_SKIP_UNREAD_J=0` turns it off
SKIP_UNREAD_J = os.environ.get('PYCIRCUIT_SKIP_UNREAD_J', '1') != '0'
STATUS = 'not loaded'

CORE_C = r"""
#include <stdint.h>
#include <string.h>
typedef void (*hdl_fn_t)(const double *, const double *, double *);
typedef void (*dgemv_t)(int, int, int64_t, int64_t, double, const double *, int64_t,
                        const double *, int64_t, double, double *, int64_t);
typedef struct {
    long n;
    const int64_t *p_L; const int64_t *p_toff; const double *tmplall;
    const int64_t *flatall; const int64_t *p_nbins;
    long nb; const int64_t *b_pass; void **b_fn; const int64_t *b_ne;
    const int64_t *b_k; const int64_t *b_so; const int64_t *b_ti;
    const int64_t *NMall; const int64_t *b_nmoff; double **PRall;
    const int64_t *b_proff; const int64_t *bdstall; const int64_t *b_dstoff;
    long ng; const int64_t *g_pass; const int64_t *g_E; const int64_t *g_k;
    const double *Sall; const int64_t *g_soff; const int64_t *Iall;
    const int64_t *g_ioff; const int64_t *g_nsel; const int64_t *selall;
    const int64_t *gdstall; const int64_t *g_seloff;
    void *dgemv;
    double *work; double *ywork; double *binG; double *bini;
} hdl_core_t;
/* passes 0 G, 1 C, 2 i, 3 q; want bits 1 G, 2 C, 4 i (q always) */
long hdl_fn(const hdl_core_t *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J)
{
    double X[64], o[4096], y[64], xg[64];
    double *bins[4] = {c->binG, C, c->bini, q};
    const long bits[4] = {1, 2, 4, 0};
    dgemv_t gemv = (dgemv_t) c->dgemv;
    long n = c->n, m, b, g, e, j, s;
    const double *Cs;
    for (m = 0; m < 4; m++) {
        double *buf = c->work;
        long L = c->p_L[m];
        if (bits[m] && !(want & bits[m])) continue;
        if (m < 2) memcpy(buf, c->tmplall + c->p_toff[m], (size_t) L * sizeof(double));
        else memset(buf, 0, (size_t) L * sizeof(double));
        for (b = 0; b < c->nb; b++) {
            if (c->b_pass[b] != m) continue;
            hdl_fn_t f = (hdl_fn_t) c->b_fn[b];
            long ne = c->b_ne[b], k = c->b_k[b], so = c->b_so[b], ti = c->b_ti[b];
            const int64_t *NM = c->NMall + c->b_nmoff[b];
            const int64_t *dst = c->bdstall + c->b_dstoff[b];
            double **PR = c->PRall + c->b_proff[b];
            for (e = 0; e < ne; e++) {
                for (j = 0; j < k; j++) X[j] = x[NM[e * k + j]];
                if (ti >= 0) PR[e][ti] = T;
                f(X, PR[e], o);
                for (j = 0; j < so; j++) buf[dst[e * so + j]] = o[j];
            }
        }
        for (g = 0; g < c->ng; g++) {
            if (c->g_pass[g] != m) continue;
            long E = c->g_E[g], k = c->g_k[g], nsel = c->g_nsel[g];
            const double *S = c->Sall + c->g_soff[g];
            const int64_t *I = c->Iall + c->g_ioff[g];
            const int64_t *sel = c->selall + c->g_seloff[g];
            const int64_t *gdst = c->gdstall + c->g_seloff[g];
            for (e = 0; e < E; e++) {
                for (j = 0; j < k; j++) xg[j] = x[I[e * k + j]];
                gemv(101, 111, k, k, 1.0, S + e * k * k, k, xg, 1, 0.0, y, 1);
                for (j = 0; j < k; j++) c->ywork[e * k + j] = y[j];
            }
            for (s = 0; s < nsel; s++) buf[gdst[s]] = c->ywork[sel[s]];
        }
        {
            const int64_t *flat = c->flatall + c->p_toff[m];
            double *bin = bins[m];
            long nbins = c->p_nbins[m];
            memset(bin, 0, (size_t) nbins * sizeof(double));
            for (s = 0; s < L; s++) bin[flat[s]] += buf[s];
        }
    }
    Cs = (want & 2) ? C : Cin;
    if (formula == 0) {
        for (j = 0; j < n; j++) iq[j] = (a0 * q[j] + a1 * q1[j]) + a2 * q2[j];
        for (s = 0; s < n * n; s++) Geq[s] = a0 * Cs[s];
    } else if (formula == 1) {
        for (j = 0; j < n; j++) iq[j] = (q[j] - q1[j]) / h;
        for (s = 0; s < n * n; s++) Geq[s] = Cs[s] / h;
    } else if (formula == 2) {
        for (j = 0; j < n; j++) iq[j] = (2.0 * (q[j] - q1[j])) / h - iq1[j];
        for (s = 0; s < n * n; s++) Geq[s] = (2.0 * Cs[s]) / h;
    } else {
        double th = theta * h, w = (1.0 - theta) / theta;
        for (j = 0; j < n; j++) iq[j] = (q[j] - q1[j]) / th - w * iq1[j];
        for (s = 0; s < n * n; s++) Geq[s] = Cs[s] / th;
    }
    if (want & 4) for (j = 0; j < n; j++) F[j] = (c->bini[j] + iq[j]) + u[j];
    if (want & 1) for (s = 0; s < n * n; s++) J[s] = c->binG[s] + Geq[s];
    return 0;
}
"""
CORE_CDEF = """
typedef struct {
    long n;
    const int64_t *p_L; const int64_t *p_toff; const double *tmplall;
    const int64_t *flatall; const int64_t *p_nbins;
    long nb; const int64_t *b_pass; void **b_fn; const int64_t *b_ne;
    const int64_t *b_k; const int64_t *b_so; const int64_t *b_ti;
    const int64_t *NMall; const int64_t *b_nmoff; double **PRall;
    const int64_t *b_proff; const int64_t *bdstall; const int64_t *b_dstoff;
    long ng; const int64_t *g_pass; const int64_t *g_E; const int64_t *g_k;
    const double *Sall; const int64_t *g_soff; const int64_t *Iall;
    const int64_t *g_ioff; const int64_t *g_nsel; const int64_t *selall;
    const int64_t *gdstall; const int64_t *g_seloff;
    void *dgemv;
    double *work; double *ywork; double *binG; double *bini;
} hdl_core_t;
long hdl_fn(const hdl_core_t *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J);
"""

PASSES = ('G', 'C', 'i', 'q')
#: the integrators the core's companion enum covers, by class name
FORMULA = {'Gear2Integrator': 0, 'EulerIntegrator': 1,
           'TrapezoidalIntegrator': 2, 'ThetaIntegrator': 3}
WANT = {'fj': 7, 'j': 3, 'f': 4, 'c': 2}

_driver = None


def _numpy_dgemv():
    """The address of `cblas_dgemv` in the OpenBLAS numpy itself loaded
    (`numpy.libs/libscipy_openblas64_*.so`, the instance already mapped:
    `CDLL` returns it), or None where the wheel has no such library."""
    try:
        d = os.path.join(os.path.dirname(np.__file__), '..', 'numpy.libs')
        paths = glob.glob(os.path.join(d, 'libscipy_openblas64_*.so'))
        if not paths:
            return None
        lib = ctypes.CDLL(os.path.realpath(paths[0]))
        return ctypes.cast(lib.scipy_cblas_dgemv64_, ctypes.c_void_p).value
    except (OSError, AttributeError):
        return None


def driver():
    """`(ffi, cfn, dgemv address)` once loaded, or None (`STATUS` says why)."""
    global _driver, STATUS
    if _driver is None:
        from pycircuit.circuit import _hdl_cbackend as cb
        dgemv = _numpy_dgemv()
        if dgemv is None:
            _driver = False
            STATUS = "off (numpy's own cblas_dgemv not found)"
        else:
            try:
                ffi, cfn, _key, _cold, _secs = cb.load_kernel(CORE_C, CORE_CDEF)
            except (cb.CompileError, OSError) as e:
                _driver = False
                STATUS = f'off ({e})'
            else:
                _driver = (ffi, cfn, dgemv)
                STATUS = 'c'
    return _driver or None


class Unservable(Exception):
    """The circuit has an element the core cannot evaluate."""


class _Core:
    """One circuit's four plans flattened into the C struct, the batches'
    pointer tables and the core's own pack mirrors."""

    def __init__(self, cir, plan, ffi, dgemv):
        from pycircuit.circuit import _hdl_batch, _stamp_plan
        from pycircuit.circuit import _hdl_cbackend as cb
        self.plan = plan
        self.n = n = cir.n
        for m in PASSES:
            if plan.methods.get(m) is None:
                plan.methods[m] = (_stamp_plan._MatrixPlan(cir, m) if m in ('G', 'C')
                                   else _stamp_plan._VectorPlan(cir, m))
        p_L, p_toff, tmpl, flat, p_nbins = [], [], [], [], []
        batches, b_pass, NM, b_nmoff, dst, b_dstoff, b_proff = [], [], [], [], [], [], []
        g_pass, S, g_soff, I, g_ioff, sel, gdst, g_seloff = [], [], [], [], [], [], [], []
        off = nmoff = dstoff = pr_len = soff = ioff = seloff = 0
        maxy = 1
        for mi, m in enumerate(PASSES):
            mp = plan.methods[m]
            mat = m in ('G', 'C')
            extra = []
            for c_ in mp.calls:
                ## a lone C-bound element is not batched by `split` (two or
                ## more): the core serves it as a batch of one
                _inst, el, nm, a, b = c_
                cls = type(el)
                info = getattr(cls, '_hdl_info', None)
                fn = (info or {}).get('funcs', {}).get(m)
                kern = fn.__dict__.get('_hdl_c') if fn is not None else None
                k = len(nm)
                if not (a >= 0 and info is not None and info.get('chained')
                        and not info['state_meta']['dc_pins'] and el.toolkit is cir.toolkit
                        and m not in el.__dict__ and _hdl_batch.is_generated(cls, m)
                        and isinstance(kern, cb.CKernel) and kern.nx == k
                        and kern.shape == ((k, k) if mat else (k,))
                        and b - a == (k * k if mat else k)):
                    raise Unservable(f'{_inst} in the {m} pass')
                extra.append(_hdl_batch.Batch(cls, m, [c_], kern))
            if mat:
                L = mp.template.size
                tmpl.append(mp.template)
                flat.append(mp.flat)
                p_nbins.append(n * n)
            else:
                L = mp.length
                tmpl.append(np.zeros(L))
                flat.append(mp.idx)
                p_nbins.append(n)
                for Sg, Ig, selg, dstg, ok in mp.groups:
                    if not ok:
                        raise Unservable(f'a constant group of the {m} pass is not matmul-exact')
                    E, k = Ig.shape
                    g_pass.append(mi)
                    S.append(Sg.reshape(-1))
                    g_soff.append((soff, E, k))
                    soff += Sg.size
                    I.append(Ig.reshape(-1))
                    g_ioff.append(ioff)
                    ioff += Ig.size
                    sel.append(selg)
                    gdst.append(dstg)
                    g_seloff.append((seloff, selg.size))
                    seloff += selg.size
                    maxy = max(maxy, E * k)
            p_L.append(L)
            p_toff.append(off)
            off += L
            for bt in list(mp.batches) + extra:
                batches.append(bt)
                b_pass.append(mi)
                NM.append(bt.NM.reshape(-1))
                b_nmoff.append(nmoff)
                nmoff += bt.NM.size
                dst.append(bt.dst)
                b_dstoff.append(dstoff)
                dstoff += bt.dst.size
                b_proff.append(pr_len)
                pr_len += bt.n
        def i64(a):
            return np.ascontiguousarray(np.asarray(a, dtype=np.int64))

        def cat64(parts):
            return (np.ascontiguousarray(np.concatenate(parts).astype(np.int64))
                    if parts else np.zeros(1, np.int64))
        self.batches = batches
        self.nb, self.ng = len(batches), len(g_pass)
        self.b_fn = np.zeros(max(self.nb, 1), dtype=np.uintp)
        self.PRall = np.zeros(max(pr_len, 1), dtype=np.uintp)
        self.b_proff = b_proff
        arrays = {
            'p_L': i64(p_L), 'p_toff': i64(p_toff), 'p_nbins': i64(p_nbins),
            'tmplall': np.ascontiguousarray(np.concatenate(tmpl)), 'flatall': cat64(flat),
            'b_pass': i64(b_pass), 'b_ne': i64([bt.n for bt in batches]),
            'b_k': i64([bt.k for bt in batches]), 'b_so': i64([bt.so for bt in batches]),
            'b_ti': i64([bt.ti for bt in batches]), 'NMall': cat64(NM), 'b_nmoff': i64(b_nmoff),
            'b_proff': i64(b_proff), 'bdstall': cat64(dst), 'b_dstoff': i64(b_dstoff),
            'g_pass': i64(g_pass), 'g_E': i64([e for _o, e, _k in g_soff]),
            'g_k': i64([k for _o, _e, k in g_soff]),
            'Sall': np.ascontiguousarray(np.concatenate(S)) if S else np.zeros(1),
            'g_soff': i64([o for o, _e, _k in g_soff]), 'Iall': cat64(I), 'g_ioff': i64(g_ioff),
            'g_nsel': i64([sz for _o, sz in g_seloff]), 'selall': cat64(sel), 'gdstall': cat64(gdst),
            'g_seloff': i64([o for o, _sz in g_seloff]),
            'work': np.empty(max(int(max(p_L)), 1)), 'ywork': np.empty(maxy),
            'binG': np.empty(n * n), 'bini': np.empty(n),
        }
        if any(bt.k > 64 or bt.so > 4096 for bt in batches) or any(k > 64 for _o, _e, k in g_soff):
            raise Unservable('an element wider than the core holds')
        self.arrays = arrays
        cs = ffi.new('hdl_core_t *')
        cs.n = n
        cs.nb, cs.ng = self.nb, self.ng
        self.handles = []
        for name, arr in arrays.items():
            h = ffi.from_buffer('double *' if arr.dtype == np.float64 else 'int64_t *', arr)
            self.handles.append(h)
            setattr(cs, name, h)
        self.handles.append(ffi.from_buffer('void **', self.b_fn))
        cs.b_fn = self.handles[-1]
        self.handles.append(ffi.from_buffer('double **', self.PRall))
        cs.PRall = self.handles[-1]
        cs.dgemv = ffi.cast('void *', dgemv)
        self.cs = cs
        self.ffi = ffi
        ## the unique elements across the batches, each with its positions in
        ## the pointer table: one shadow-and-pack check per call
        pos = {}
        for b, bt in enumerate(batches):
            for e, el in enumerate(bt.els):
                pos.setdefault(id(el), (el, []))[1].append(b_proff[b] + e)
        self.uniq = [(el, tuple(ps)) for el, ps in pos.values()]
        self.mirror = [None] * len(self.uniq)
        self.dptr = ffi.typeof('double *')

    def ready(self, epar):
        """`Batch.run`'s checks, once: per batch the class (bound, the
        kernel the one cast, the method generated), per unique element the
        shadows of all four methods and the pack mirror.  The temperature,
        or None (the Python path), counted."""
        T = self.probe(epar)
        if T.__class__ is str:
            return _no('core.ready:' + T)
        return T

    def probe(self, epar):
        """`ready` without the count: the temperature, or why not (a str).
        For a caller that falls back to a path which asks `ready` itself
        (`_tran_newton_c`): the decline is counted once, there."""
        T = getattr(epar, 'T', 300.0)
        if type(T) is not float:
            if type(T) is not int and np.ndim(T) != 0:
                return 'T'
            T = float(T)
        from pycircuit.circuit import _hdl_batch
        ffi, b_fn = self.ffi, self.b_fn
        kern0 = None
        for b, bt in enumerate(self.batches):
            info = bt.info
            if not info.get('_c_bound'):
                return 'unbound'
            kern = info['funcs'][bt.m].__dict__.get('_hdl_c')
            if kern is not bt.kern:
                if not bt._take(kern):
                    return 'kernel'
                b_fn[b] = int(ffi.cast('uintptr_t', kern.cfn))
            elif not b_fn[b]:
                b_fn[b] = int(ffi.cast('uintptr_t', kern.cfn))
            if not _hdl_batch.is_generated(bt.cls, bt.m):
                return 'generated'
            if kern0 is None:
                kern0 = kern
        mirror, PR = self.mirror, self.PRall
        for ui, (el, ps) in enumerate(self.uniq):
            d = el.__dict__
            if 'G' in d or 'C' in d or 'i' in d or 'q' in d:
                return 'shadow'
            cp = d.get('_hdl_cp')
            if cp is None:
                try:
                    cp = kern0.pack(el)
                except (TypeError, ValueError):
                    cp = False
                d['_hdl_cp'] = cp
            if cp is False:
                return 'pack'
            if cp is not mirror[ui]:
                mirror[ui] = cp
                a = cp[0].ctypes.data
                for p_ in ps:
                    PR[p_] = a
        return T


def core_for(tr):
    """The transient's core for its circuit's current plan, or None where
    the circuit cannot be served (remembered until the plan rebuilds)."""
    from pycircuit.circuit import _stamp_plan
    cir = tr.cir
    plan = _stamp_plan._plan_for(cir)
    rec = tr.__dict__.get('_tran_core')
    if rec is not None and rec[0] is plan:
        return rec[1]
    drv = driver()
    core = None
    if drv is not None:
        ffi, _cfn, dgemv = drv
        try:
            core = _Core(cir, plan, ffi, dgemv)
        except Unservable:
            core = None
    _PC['once:core.build:' + ('nodriver' if drv is None else
                              'unservable' if core is None else 'built')] += 1
    tr.__dict__['_tran_core'] = (plan, core)
    return core


def _row(a, n):
    """A history row as the C reads it, or None."""
    if type(a) is np.ndarray and a.dtype == np.float64 and a.ndim == 1 and a.shape[0] == n:
        return a if a.flags.c_contiguous else np.ascontiguousarray(a)
    return None


def evaluate(tr, x, t, provided_function, want):
    """`want` 'fj': `(f, J)` as `_residual_and_jacobian`; 'j': `(None, J)`
    as `jacobian_only`; 'c': `(None, None)`, `jacobian_only`'s state where
    nothing reads its `J` (C and q; no G); 'f': `f` as `residual_only` -- or
    None where the core does not serve the call (nothing touched; the
    Python path runs)."""
    K = _K[want]
    if not CORE:
        return _no(K['off'])
    from pycircuit.circuit.toolkit import NumericToolkit
    if type(tr.toolkit) is not NumericToolkit:
        return _no(K['toolkit'])
    ## THE METHODS THE CORE STANDS IN FOR MUST BE THE ONES IT MIRRORS: an
    ## instance shadow of the circuit's passes (a test counting `cir.G`, a
    ## harness timing it) or of the transient's companion machinery (a spy
    ## on `get_diff`) is a caller that expects those calls -- the same rule
    ## as a batch's on an element's shadow (`_hdl_batch`)
    cd = tr.cir.__dict__
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd:
        return _no(K['shadow_cir'])
    td = tr.__dict__
    if ('get_diff' in td or '_companion_at' in td or '_C_at_state' in td
            or '_C_lookup' in td or '_source_at' in td):
        return _no(K['shadow_tr'])
    core = core_for(tr)
    if core is None:
        return _no(K['unservable'])
    n = core.n
    if not (type(x) is np.ndarray and x.dtype == np.float64 and x.ndim == 1
            and x.flags.c_contiguous and x.shape[0] == n and np.isfinite(x).all()):
        return _no(K['x'])
    h = tr._dt
    h_last = tr._dt_last if tr._dt_last is not None else h
    active = tr.base_integrator.check_order_drop(h, h_last, tr._is_first_step)
    formula = FORMULA.get(type(active).__name__)
    if formula is None:
        return _no(K['formula'])
    q1 = _row(tr._qlast[0], n)
    q2 = _row(tr._qlast[1], n) if formula == 0 else q1
    iq1 = _row(tr._iqlast[0], n) if formula in (2, 3) else q1
    if q1 is None or q2 is None or iq1 is None:
        return _no(K['history'])
    a0 = a1 = a2 = theta = 0.0
    if formula == 0:
        from pycircuit.circuit._lte_kernels import bdf2_alphas
        a0, a1, a2 = bdf2_alphas(h, h_last)
    elif formula == 3:
        theta = active.theta_at(h)
    bits = WANT[want]
    if want == 'f':
        Cin = tr._Cmat
        if not (type(Cin) is np.ndarray and Cin.dtype == np.float64 and Cin.shape == (n, n)
                and Cin.flags.c_contiguous):
            return _no(K['cmat'])
        C = Cin
    else:
        Cin = tr._C_lookup(x)
        if Cin is not None:
            if not (type(Cin) is np.ndarray and Cin.dtype == np.float64
                    and Cin.shape == (n, n) and Cin.flags.c_contiguous):
                return _no(K['clookup'])
            bits &= ~2
            C = Cin
        else:
            C = np.empty((n, n))
            Cin = C
    T = core.ready(tr.epar)
    if T is None:
        return _no(K['ready'])
    u = None
    if bits & 4:
        u = tr._source_at(t, provided_function)
        if not (type(u) is np.ndarray and u.dtype == np.float64 and u.ndim == 1
                and u.shape[0] == n and u.flags.c_contiguous):
            return _no(K['u'])
    _PC[K['served']] += 1
    ffi, cfn, _dgemv = driver()
    fb, dptr = ffi.from_buffer, core.dptr
    q, iq, Geq = np.empty(n), np.empty(n), np.empty((n, n))
    F = np.empty(n) if bits & 4 else q
    J = np.empty((n, n)) if bits & 1 else Geq
    cfn(core.cs, fb(dptr, x), T, bits, formula, a0, a1, a2, float(h), theta,
        fb(dptr, q1), fb(dptr, q2), fb(dptr, iq1), fb(dptr, u if u is not None else q1),
        fb(dptr, Cin), fb(dptr, C), fb(dptr, q), fb(dptr, iq), fb(dptr, Geq),
        fb(dptr, F), fb(dptr, J))
    ## the state the Python path leaves behind
    if want != 'f':
        if bits & 2 and tr._memo_ok():
            tr._memo_put(x, {'C': C})
        tr._q_cache = (x, q)
        tr._C_cache = (x, C)
    tr.active_integrator = active
    tr._iq = iq
    tr._Geq = Geq
    tr._Cmat = C
    tr._companion_coeffs = active.companion_coefficients(h, h_last)
    tr._effective_method = type(active).__name__
    if want == 'fj':
        return F, J
    if want == 'j':
        return None, J
    if want == 'c':
        return None, None
    return F
