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

THE FUSED KERNELS (speed round 9, stage 3).  A class's four passes print
four kernels that each recompute their whole chain; its fused kernel
(`_hdl_cbackend.fuse_csrc`: the printed statements unioned by name, each
run under the bits of the passes that print it) is called ONCE a call for
every pass the call wants, before the passes, into a staging buffer --
where the class's batches cover two or more passes over the same elements
on the same nodes (`_fuse_groups`) -- and each batch then scatters its
outputs from there in its own order: the same values in the same places,
the passes' bins as before.  The readiness check takes the class's fused
kernel only while it is bound and was printed from the passes' bound
kernels (and the switch, `_hdl_cbackend.FUSE`, is on); else the batches
call their own.

READY AS STAMPED (speed round 9, stage 4).  The readiness check (`probe`:
every batch's kernel the class's bound one, its pass the generated one,
every element without a shadow, every pack mirrored) and the lookup
(`core_for`: the circuit's plan) read the same dicts every call, and their
answer does not change between steps: a check that passed watches what it
read (`_watch`: one counter any change to a watched dict moves) and passes
again unchecked while the counter stands -- what no dict watcher sees
compared on every call: the elements' `__dict__` objects, the passes'
code, `ParameterDict`'s epoch, the circuit's size.

`CORE` (env `PYCIRCUIT_TRAN_CORE=0`) is read on every call; `STATUS` says
why the core is off where it is.  Pinned against its Python reference in
`tests/test_pinned_pairs.py`.  History: `doc/transient_history.md`,
`_tran_core.py`; `doc/hdl_roadmap_260824.md` sec. 65.
"""
import ctypes
import glob
import operator
import os

import numpy as np

from pycircuit.circuit import _lte_kernels, _paths, _watch
from pycircuit.circuit import toolkit as _toolkit

## the counters (`_paths`): each decline by its reason, each served call
_PC = _paths.COUNTS
_GETDICT = operator.attrgetter('__dict__')
_no = _paths.no
_K = {w: {r: f'core.{w}:{r}' for r in (
    'off', 'toolkit', 'shadow_cir', 'shadow_tr', 'unservable', 'x', 'formula',
    'history', 'cmat', 'clookup', 'ready', 'u', 'served')} for w in ('fj', 'j', 'f', 'c', 'p')}

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
typedef void (*hdl_fz_t)(const double *, const double *, long, double *, double *, double *,
                         double *);
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
    long nz; void **z_fn; const int64_t *ztab; double *zst;
    const int64_t *bftab; const int64_t *fidxall;
} hdl_core_t;
/* passes 0 G, 1 C, 2 i, 3 q; want bits 1 G, 2 C, 4 i (q always).  A class
   with a fused kernel (z_fn) is evaluated once a call, before the passes,
   for every pass the call wants, into its staging (zst); its batches then
   scatter from there -- the passes' writes and bins as before.  ztab, a
   group's row: ne, k, ti, mask, nmoff, proff, then per pass G C i q its
   output length and its staging offset; bftab, a batch's row: its group
   (-1: none) and its offset in fidxall (each element's place in it; -1:
   the group's own order, the usual case) */
#define ZW 14
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
    long n = c->n, m, b, g, e, j, s, z;
    const double *Cs;
    for (z = 0; z < c->nz; z++) {
        const int64_t *zt = c->ztab + ZW * z;
        hdl_fz_t ff = (hdl_fz_t) c->z_fn[z];
        long zw = ((want & 7) | 8) & zt[3];
        if (!ff || !zw) continue;
        {
            long ne = zt[0], k = zt[1], ti = zt[2];
            const int64_t *NM = c->NMall + zt[4];
            double **PR = c->PRall + zt[5];
            double *sG = c->zst + zt[7], *sC = c->zst + zt[9], *si = c->zst + zt[11];
            double *sq = c->zst + zt[13];
            for (e = 0; e < ne; e++) {
                for (j = 0; j < k; j++) X[j] = x[NM[e * k + j]];
                if (ti >= 0) PR[e][ti] = T;
                ff(X, PR[e], zw, sG + e * zt[6], sC + e * zt[8], si + e * zt[10], sq + e * zt[12]);
            }
        }
    }
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
            long fz = c->nz ? c->bftab[2 * b] : -1;
            const int64_t *NM = c->NMall + c->b_nmoff[b];
            const int64_t *dst = c->bdstall + c->b_dstoff[b];
            double **PR = c->PRall + c->b_proff[b];
            if (fz >= 0 && c->z_fn[fz]) {
                const double *sz = c->zst + c->ztab[ZW * fz + 7 + 2 * m];
                long fo = c->bftab[2 * b + 1];
                const int64_t *fi = fo >= 0 ? c->fidxall + fo : 0;
                for (e = 0; e < ne; e++) {
                    const double *oz = sz + (fi ? fi[e] : e) * so;
                    for (j = 0; j < so; j++) buf[dst[e * so + j]] = oz[j];
                }
                continue;
            }
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
    /* formula 4: the passes alone (`passes`: G in binG, i in bini, C and q
       in the caller's buffers) -- no companion, no residual, no Jacobian */
    if (formula == 4) return 0;
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
    long nz; void **z_fn; const int64_t *ztab; double *zst;
    const int64_t *bftab; const int64_t *fidxall;
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


def _fuse_groups(batches, b_pass, cb):
    """The classes the core evaluates with their fused kernel:
    `(groups, b_fz, b_fidx)` -- each group `(info, kernel, elements, pass
    bits, node maps)`, per batch its group (-1: none) and each of its
    elements' place in the group (None: the group's own order).  A class is fused where its fused kernel
    is bound and its batches cover two or more of its passes, the same
    elements in each, once each, on the same nodes, with the kernel's
    widths; each refusal counted once a build (`once:core.fuse:<why>`).
    (Every pass lists a class's elements in one order, the plan's: the
    node maps are then compared a batch at a time -- an element at a time
    cost a 20-element chain's core build 1.7 M instructions.)"""
    by_cls = {}
    for b, bt in enumerate(batches):
        by_cls.setdefault(bt.cls, []).append(b)
    groups, b_fz, b_fidx = [], [-1] * len(batches), [None] * len(batches)
    for bs in by_cls.values():
        first = batches[bs[0]]
        info = first.info
        fk = info.get('_c_fused') if cb.FUSE else None
        if fk is None and cb.FUSE and info.get('_c_bound') and info.get('_c_fused_status') == 'off':
            ## (bound while the switch was off: built now)
            fk = cb.bind_fused(info, cb.fuse_source(info))
        passes = [PASSES[b_pass[b]] for b in bs]
        why = None
        if fk is None:
            why = 'nokernel'
        elif len(set(passes)) != len(passes):
            why = 'passes'          # (a class split into two batches in one pass)
        elif len(passes) < 2 or not set(passes) <= set(fk.parts):
            why = 'passes'
        else:
            ti = -1 if fk.layout[1] is None else fk.layout[1]
            for b, m in zip(bs, passes):
                bt = batches[b]
                if bt.k != fk.nx or bt.so != fk.sizes[m] or bt.ti != ti:
                    why = 'widths'
                    break
        els, places = first.els, {}
        if why is None:
            nm0 = first.NM.tobytes()
            index = None
            for b in bs:
                bt = batches[b]
                if bt.els is els or (len(bt.els) == len(els)
                                     and all(map(operator.is_, bt.els, els))):
                    ## (the plan's one order: the node maps compared whole,
                    ## both int64 and C-ordered as `Batch` makes them)
                    if bt is not first and not (bt.NM.shape == first.NM.shape
                                                and bt.NM.tobytes() == nm0):
                        why = 'nodes'
                        break
                    places[b] = None
                    continue
                if index is None:
                    index = {id(el): i for i, el in enumerate(els)}
                rows = [index.get(id(el), -1) for el in bt.els]
                if (len(index) != len(els) or len(rows) != len(els) or -1 in rows
                        or len(set(rows)) != len(rows)):
                    why = 'elements'
                    break
                rows = np.asarray(rows, dtype=np.int64)
                if not np.array_equal(bt.NM, first.NM[rows]):
                    why = 'nodes'
                    break
                places[b] = rows
        if why is not None:
            _PC['once:core.fuse:' + why] += 1
            continue
        z = len(groups)
        mask = sum(1 << PASSES.index(m) for m in passes)
        groups.append((info, fk, els, mask, first.NM))
        for b in bs:
            b_fz[b] = z
            b_fidx[b] = places[b]
        _PC['once:core.fuse:fused'] += 1
    return groups, b_fz, b_fidx


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
        ## THE FUSED KERNELS (speed round 9, stage 3): a class whose fused
        ## kernel is bound (`_hdl_cbackend.bind_fused`) and whose batches
        ## cover two or more passes over the same elements -- each once a
        ## pass, on the same nodes -- is evaluated once a call into a
        ## staging buffer; its batches scatter from there
        zgroups, b_fz, b_fidx = _fuse_groups(batches, b_pass, cb)
        self._cb = cb
        z_proff, ztab, fidx, bftab = [], [], [], []
        st_len = fo = 0
        for _info, fk, els, mask, znm in zgroups:
            ne = len(els)
            ztab += (ne, fk.nx, -1 if fk.layout[1] is None else fk.layout[1], mask,
                     nmoff, pr_len)
            z_proff.append(pr_len)
            pr_len += ne
            NM.append(znm.reshape(-1))
            nmoff += znm.size
            for m in PASSES:
                sz = fk.sizes.get(m, 0)
                ztab += (sz, st_len)
                st_len += sz * ne
        if zgroups:
            for b in range(len(batches)):
                rows = b_fidx[b]
                if rows is None:
                    ## (the group's own order: no table)
                    bftab += (b_fz[b], -1)
                else:
                    bftab += (b_fz[b], fo)
                    fidx.append(rows)
                    fo += len(rows)
        self.zgroups = zgroups
        self.nz = len(zgroups)
        self.b_fz = b_fz
        ## (`probe`'s stamp: none yet)
        self.stamp, self.uniq_els, self.uniq_dicts, self.codes, self.arms = -1, (), (), (), 0
        ## (the checks the stamp passed since it was armed: `_watch.counted`)
        self.held = 0
        self.z_fn = np.zeros(max(self.nz, 1), dtype=np.uintp)
        self.zk = [None] * self.nz
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
        if zgroups:
            ## (the fused groups' tables: none at all where nothing fuses)
            arrays.update(ztab=np.array(ztab, dtype=np.int64), zst=np.empty(max(st_len, 1)),
                          bftab=np.array(bftab, dtype=np.int64), fidxall=cat64(fidx))
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
        cs.nz = self.nz
        self.handles.append(ffi.from_buffer('void **', self.z_fn))
        cs.z_fn = self.handles[-1]
        cs.dgemv = ffi.cast('void *', dgemv)
        self.cs = cs
        self.ffi = ffi
        ## the unique elements across the batches, each with its positions in
        ## the pointer table (a fused group's too): one shadow-and-pack check
        ## per call
        pos = {}
        for b, bt in enumerate(batches):
            for e, el in enumerate(bt.els):
                pos.setdefault(id(el), (el, []))[1].append(b_proff[b] + e)
        for z, g in enumerate(zgroups):
            for e, el in enumerate(g[2]):
                pos[id(el)][1].append(z_proff[z] + e)
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
        ## AS STAMPED (speed round 9, stage 4): nothing the check below read
        ## has changed since it passed -- its dicts and types are watched
        ## (`_watch`); what no watcher sees, compared here: the elements'
        ## `__dict__` objects and the passes' code
        ep = _watch.EPOCH
        if (ep is not None and self.stamp == ep.value
                and tuple(map(_GETDICT, self.uniq_els)) == self.uniq_dicts
                and all(map(_same_code, self.codes))):
            self.held += 1
            return T
        before = _watch.now()
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
                if self.b_fz[b] >= 0:
                    ## (its group's fused kernel is checked again below)
                    self.zk[self.b_fz[b]] = None
            elif not b_fn[b]:
                b_fn[b] = int(ffi.cast('uintptr_t', kern.cfn))
            if not _hdl_batch.is_generated(bt.cls, bt.m):
                return 'generated'
            if kern0 is None:
                kern0 = kern
        if self.nz:
            ## A FUSED KERNEL SERVES while it is the class's bound one and
            ## was printed from the passes' bound kernels; else the batches
            ## call their own (the switch `_hdl_cbackend.FUSE` too)
            zk, z_fn, on = self.zk, self.z_fn, self._cb.FUSE
            for z, g in enumerate(self.zgroups):
                info = g[0]
                fk = info.get('_c_fused') if on else None
                if fk is not zk[z] or (fk is None and z_fn[z]):
                    zk[z] = fk
                    funcs, fk0 = info['funcs'], self.zgroups[z][1]
                    ok = (fk is not None and fk.sizes == fk0.sizes and fk.nx == fk0.nx
                          and fk.layout == fk0.layout and all(
                              k is funcs[m].__dict__.get('_hdl_c')
                              for m, k in fk.parts.items()))
                    z_fn[z] = int(ffi.cast('uintptr_t', fk.cfn)) if ok else 0
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
        self._stamp(before, _hdl_batch)
        return T

    def _stamp(self, before, _hdl_batch):
        """Watch what the check read and stamp it (`_watch.arm`): the
        classes' info dicts and pass lists, the passes' function dicts (their
        bound kernels), the classes, the elements' dicts, the backend's and
        the batches' modules (their switches and drivers)."""
        infos = {id(bt.info): bt.info for bt in self.batches}
        dicts = list(infos.values())
        dicts += [info['funcs'] for info in infos.values()]
        dicts += [info['funcs'][bt.m].__dict__ for bt in self.batches
                  for info in (bt.info,)]
        self.uniq_els = tuple(el for el, _ps in self.uniq)
        self.uniq_dicts = tuple(map(_GETDICT, self.uniq_els))
        dicts += self.uniq_dicts
        dicts += [_hdl_batch.__dict__, self._cb.__dict__]
        self.codes = tuple({(bt.cls, bt.m): (bt.cls, bt.m, getattr(bt.cls, bt.m).__code__)
                            for bt in self.batches}.values())
        self.arms, self.held = _watch.counted(self.arms, self.held), 0
        self.stamp = _watch.arm(dicts, before) if self.arms <= _watch.MAX_ARMS else -1


def _same_code(c):
    """A class's pass still the code it was at the stamp (a function's
    `__code__` is no dict)."""
    return getattr(c[0], c[1]).__code__ is c[2]


def core_for(tr):
    """The transient's core for its circuit's current plan, or None where
    the circuit cannot be served (remembered until the plan rebuilds).
    As stamped (`_watch`): the circuit's dict, its elements and its node map
    unchanged since, the circuit's plan the one recorded, `ParameterDict`'s
    epoch the plan's and the circuit's size the plan's -- the record stands
    without asking `_plan_for`."""
    cir = tr.cir
    rec = tr.__dict__.get('_tran_core')
    ep = _watch.EPOCH
    if (rec is not None and ep is not None and rec[2] == ep.value
            and cir.__dict__.get('_stamp_plan') is rec[0] and rec[3]._epoch is rec[0].epoch
            and rec[0].n == len(cir.nodes) + len(cir.branches)):
        ## (the record `[plan, core, stamp, ParameterDict, arms, held]`: its
        ## checks passed since armed counted, `_watch.counted`)
        rec[5] += 1
        return rec[1]
    before = _watch.now()
    from pycircuit.circuit import _stamp_plan
    from pycircuit.utilities.param import ParameterDict
    plan = _stamp_plan._plan_for(cir)
    if rec is not None and rec[0] is plan:
        arms = _watch.counted(rec[4], rec[5])
        tr.__dict__['_tran_core'] = [plan, rec[1], _plan_stamp(cir, before, arms),
                                     ParameterDict, arms, 0]
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
    tr.__dict__['_tran_core'] = [plan, core, _plan_stamp(cir, before, 1), ParameterDict, 1, 0]
    return core


def _plan_stamp(cir, before, arms):
    """`_plan_for`'s dict reads, watched: the circuit's dict (its plan, its
    element and node-map attributes), the elements and the node map
    themselves.  (`ParameterDict`'s epoch, a class attribute, is compared
    on every call.)"""
    if arms > _watch.MAX_ARMS:
        return -1
    nodemap = cir.__dict__.get('elementnodemap')
    dicts = [cir.__dict__, cir.elements]
    if isinstance(nodemap, dict):
        dicts.append(nodemap)
    return _watch.arm(dicts, before)


def _row(a, n):
    """A history row as the C reads it, or None."""
    if type(a) is np.ndarray and a.dtype == np.float64 and a.ndim == 1 and a.shape[0] == n:
        return a if a.flags.c_contiguous else np.ascontiguousarray(a)
    return None


_NTK = []


def _ntk():
    """`NumericToolkit`, imported at the first call (a function-level import
    runs importlib's `_handle_fromlist` on every call)."""
    from pycircuit.circuit.toolkit import NumericToolkit
    _NTK.append(NumericToolkit)
    return NumericToolkit


#: the circuit's passes for the stage paths from the core (`passes`, speed
#: round 9, stage 8); env `PYCIRCUIT_CORE_PASSES=0` turns it off
CORE_PASSES = os.environ.get('PYCIRCUIT_CORE_PASSES', '1') != '0'
_PASS_BITS = {'G': 1, 'C': 2, 'i': 4, 'q': 0}


def passes(tr, x, which):
    """`{m: cir.m(x, epar)}` for each pass `m` of `which` (letters of
    'GCiq') -- fresh float64 arrays, bit for bit the circuit's passes --
    from ONE core call (formula 4: the passes alone, through the classes'
    fused kernels), or None where the core does not serve the call (nothing
    touched; the caller makes its own calls): the switches, a toolkit that
    is not numeric, an instance shadow of a pass on the circuit, a circuit
    the core cannot serve, a state that is not one finite float64 vector of
    its size, a batch not ready -- each counted (`core.p:<why>`) as
    `evaluate` counts its own.  For the stage methods (Radau, the DIRKs),
    whose step evaluates the passes at several points a step (speed round
    9, stage 8)."""
    K = _K['p']
    if not (CORE and CORE_PASSES):
        return _no(K['off'])
    if type(tr.toolkit) is not (_NTK[0] if _NTK else _ntk()):
        return _no(K['toolkit'])
    cd = tr.cir.__dict__
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd:
        return _no(K['shadow_cir'])
    ## A CIRCUIT THE CORE CANNOT SERVE declines at once under the plan its
    ## dict holds (as `core_for` remembers it per plan): `core_for`'s full
    ## check ran on every call of such a circuit -- a diode ladder by
    ## Radau, whose stamps never hold, paid 9 k instructions a call, +0.95 %
    ## of its run.  A plan gone stale is replaced at the circuit's next
    ## evaluation, which follows; until then the decline is the circuit's
    ## own passes.
    td = tr.__dict__
    kept = td.get('_passes_no')
    if kept is not None and cd.get('_stamp_plan') is kept:
        return _no(K['unservable'])
    core = core_for(tr)
    if core is None:
        plan = cd.get('_stamp_plan')
        if plan is not None:
            td['_passes_no'] = plan
        return _no(K['unservable'])
    n = core.n
    if not (type(x) is np.ndarray and x.dtype == np.float64 and x.ndim == 1
            and x.flags.c_contiguous and x.shape[0] == n and np.isfinite(x).all()):
        return _no(K['x'])
    T = core.ready(tr.epar)
    if T is None:
        return _no(K['ready'])
    bits = 0
    for m in which:
        bits |= _PASS_BITS[m]
    ffi, cfn, _dgemv = driver()
    fb, dptr = ffi.from_buffer, core.dptr
    dummy = core.__dict__.get('_pass_dummy')
    if dummy is None:
        dummy = core._pass_dummy = (np.zeros(max(n * n, 1)), )
        dummy = core._pass_dummy = (dummy[0], fb(dptr, dummy[0]))
    q = np.empty(n)
    C = np.empty((n, n)) if bits & 2 else None
    dC = fb(dptr, C) if C is not None else dummy[1]
    _PC[K['served']] += 1
    cfn(core.cs, fb(dptr, x), T, bits, 4, 0.0, 0.0, 0.0, 1.0, 0.0,
        dummy[1], dummy[1], dummy[1], dummy[1], dC, dC, fb(dptr, q),
        dummy[1], dummy[1], dummy[1], dummy[1])
    out = {}
    for m in which:
        if m == 'G':
            out['G'] = core.arrays['binG'].reshape(n, n).copy()
        elif m == 'C':
            out['C'] = C
        elif m == 'i':
            out['i'] = core.arrays['bini'].copy()
        else:
            out['q'] = q
    return out


def evaluate(tr, x, t, provided_function, want):
    """`want` 'fj': `(f, J)` as `_residual_and_jacobian`; 'j': `(None, J)`
    as `jacobian_only`; 'c': `(None, None)`, `jacobian_only`'s state where
    nothing reads its `J` (C and q; no G); 'f': `f` as `residual_only` -- or
    None where the core does not serve the call (nothing touched; the
    Python path runs)."""
    K = _K[want]
    if not CORE:
        return _no(K['off'])
    if type(tr.toolkit) is not _toolkit.NumericToolkit:
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
        a0, a1, a2 = _lte_kernels.bdf2_alphas(h, h_last)
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
