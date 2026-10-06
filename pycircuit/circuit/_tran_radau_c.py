"""RADAU'S DENSE STAGE NEWTON IN C (speed round 9, B2; 2026-10-05).

A Radau IIA(3) step's coupled ``3m`` Newton (`_coupled_stage_solver`'s
`_stage_newton`, undamped and without a shunt) spent ~5 M instructions a
step on a 20-MosLevel1 chain in Python and numpy around ~1.5 M of C: per
iteration the block system (`_coupled_stage_system`: nine blocks, each
through `remove_row_col`, ~0.7 M), the three stage evaluations' wrappers,
the limiter's wrapper per stage, numpy's solve wrapper and the convergence
test.  Here ONE C call does the Newton: per assembly the three stages'
passes through the evaluate core's own C (`_tran_core.CORE_C`, formula 4,
as `_tran_core.passes`), `K_j = -(i_j + u_j)`, ``F_i = (q_i - q_n) - h
(((0 + A_i0 K_0) + A_i1 K_1) + A_i2 K_2)`` and the blocks ``C_i + (h
A_ij) G_j`` (``(h A_ij) G_j`` off the diagonal) over the full width, then
reduced and written column-major as numpy's solve copies them; numpy's own
`scipy_dgesv_64_` on ``-R``; per stage ``Y_prev + insert(dY_i, iref,
0.0)`` (the reference row's ``+ 0.0``), the limiter walk's own driver
(`_hdl_climit._WALK_C`) against ``Y_prev`` and the limited step's largest
reduced entry; the convergence test (`_stages_converged`: the step within
``reltol`` of the largest reduced stage plus ``abstol``) after the
assembly at the new stages, as the Python assembles before it tests.
Bit for bit the Python Newton (pinned in `tests/test_pinned_pairs.py`).

TRANSACTIONAL.  The C writes nothing outside its buffers.  Where it
converges, what the Python Newton leaves is made here: every assembly's
three stages recorded in the device memo (`_memo_put`, in the Python's
order: `_finish_stage_step`, the next step's ``q_n``, the transform path
and PSS read them), the step's source memo holding the three stage times
(their first assembly's calls are made here, before the C) and its hit
count for every later assembly's.  Anything unusual -- a non-finite value,
a floating-point exception in its own arithmetic (numpy would warn), a
residual too large for its ``sum |R|`` (the undamped Newton's only use of
it), a LAPACK info, the walk stopping, more assemblies than its buffers
hold, maxiter -- hands the solve back: the source memo's entries and counts
rolled back, and the Python Newton runs from the same seed (the same
answer, exception, warning).

DECLINES, before anything is touched (`_paths.no`): the switch (and the
evaluate core's two), no build, a toolkit that is not numeric, a subclass
of the transient or a patched piece of the Newton (`_chain`), an instance
shadow of the circuit's passes or limiter or of the transient's pieces, a
`provided_function`, a stateful limiter or circuit-level limiting, a bypass,
no source memo, a circuit the core cannot serve, a limiting element that
is not its C kernel, a batch not ready, inputs of another type or size.

`ENABLED` (env `PYCIRCUIT_RADAU_C=0`) is read on every call; `STATUS` says
why the path is off where it is.  History: `doc/transient_history.md`,
`_tran_radau_c.py`.
"""
import os

import numpy as np

from pycircuit.circuit import _paths

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_RADAU_C', '1') != '0'
STATUS = 'not loaded'
#: the assemblies one call holds (a Newton of more iterations is handed back)
MAXA = 8

RADAU_C = r"""
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
typedef long (*core_fn_t)(const void *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J);
typedef long (*walk_fn_t)(void **F, double **PR, const int64_t *TI, double T,
            const int64_t *NM, const int64_t *OFF, const int64_t *K,
            double *xg, const double *x0g, long start, long n);
typedef void (*dgesv_t)(const int64_t *, const int64_t *, double *, const int64_t *,
                        int64_t *, double *, const int64_t *, int64_t *);
typedef struct {
    const void *core; void *core_fn; double T;
    const double *binG, *bini; double *dummy;
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    void *dgesv;
    long n, iref, maxiter, maxa;
    double h, reltol, abstol;
    double A[9];
    const double *qn, *u, *seed;
    double *Ya, *qa, *ia, *Ca, *Ga;
    double *K, *F, *R, *Jf, *b, *blk, *stp, *tolv;
    int64_t *ipiv;
    long assemblies, iters, status;
} radau_t;

#define RC_FLAGS (FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO)

static long rc_bail(radau_t *s, long why)
{
    feclearexcept(FE_ALL_EXCEPT);
    s->status = why;
    return why;
}

/* one assembly at the stages Ya[a]: each stage's passes through the core
   (formula 4: G in binG, i in bini, C and q in place), then the residual and
   the blocks in `_coupled_stage_system`'s operations over the full width
   (each result stored, so each operation precedes the flag test) */
static long assemble(radau_t *s, long a)
{
    long n = s->n, m = n - 1, r = s->iref, M = 3 * m, nn = n * n;
    long i, j, k, rr, cc, ri, ci;
    core_fn_t core = (core_fn_t) s->core_fn;
    double *Y = s->Ya + a * 3 * n, *q = s->qa + a * 3 * n, *iv = s->ia + a * 3 * n;
    double *C = s->Ca + a * 3 * nn, *G = s->Ga + a * 3 * nn;
    double h = s->h;
    for (j = 0; j < 3; j++) {
        /* (a non-finite stage: `passes` declines and the Python evaluates) */
        for (k = 0; k < n; k++) if (!isfinite(Y[j * n + k])) return 2;
        core(s->core, Y + j * n, s->T, 7, 4, 0.0, 0.0, 0.0, 1.0, 0.0,
             s->dummy, s->dummy, s->dummy, s->dummy, s->dummy, C + j * nn, q + j * n,
             s->dummy, s->dummy, s->dummy, s->dummy);
        memcpy(G + j * nn, s->binG, (size_t) nn * sizeof(double));
        memcpy(iv + j * n, s->bini, (size_t) n * sizeof(double));
    }
    feclearexcept(FE_ALL_EXCEPT);
    for (j = 0; j < 3; j++)
        for (k = 0; k < n; k++) s->K[j * n + k] = -(iv[j * n + k] + s->u[j * n + k]);
    for (i = 0; i < 3; i++) {
        const double a0 = s->A[3 * i], a1 = s->A[3 * i + 1], a2 = s->A[3 * i + 2];
        double *f = s->F + i * n;
        for (k = 0; k < n; k++) {
            double sm = 0.0 + a0 * s->K[k];
            sm = sm + a1 * s->K[n + k];
            sm = sm + a2 * s->K[2 * n + k];
            f[k] = (q[i * n + k] - s->qn[k]) - h * sm;
        }
        for (k = 0, rr = 0; k < n; k++) {
            if (k == r) continue;
            s->R[i * m + rr] = f[k];
            rr++;
        }
    }
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            const double hA = h * s->A[3 * i + j];
            const double *Gj = G + j * nn, *Ci = C + i * nn;
            double *blk = s->blk;
            if (i == j) for (k = 0; k < nn; k++) blk[k] = Ci[k] + hA * Gj[k];
            else for (k = 0; k < nn; k++) blk[k] = hA * Gj[k];
            for (ri = 0, rr = 0; ri < n; ri++) {
                if (ri == r) continue;
                for (ci = 0, cc = 0; ci < n; ci++) {
                    if (ci == r) continue;
                    s->Jf[(j * m + cc) * M + (i * m + rr)] = blk[ri * n + ci];
                    cc++;
                }
                rr++;
            }
        }
    }
    if (fetestexcept(RC_FLAGS)) return 3;
    return 0;
}

long hdl_fn(radau_t *s)
{
    long n = s->n, m = n - 1, r = s->iref, M = 3 * m, a = 0, it, i, k, e;
    int64_t M64 = M, one64 = 1, info64 = 0;
    walk_fn_t walk = (walk_fn_t) s->walk_fn;
    s->assemblies = 0;
    s->iters = 0;
    s->status = 0;
    /* every input finite: anything else is the Python Newton's */
    if (!(isfinite(s->h) && isfinite(s->reltol) && isfinite(s->abstol))) return rc_bail(s, 2);
    for (k = 0; k < 9; k++) if (!isfinite(s->A[k])) return rc_bail(s, 2);
    for (k = 0; k < n; k++) if (!isfinite(s->qn[k])) return rc_bail(s, 2);
    for (k = 0; k < 3 * n; k++) if (!isfinite(s->u[k])) return rc_bail(s, 2);
    memcpy(s->Ya, s->seed, (size_t) (3 * n) * sizeof(double));
    e = assemble(s, 0);
    if (e) return rc_bail(s, e);
    for (it = 0; it < s->maxiter; it++) {
        double scale = 0.0, ynorm = 0.0;
        double *Yp = s->Ya + a * 3 * n, *Yn;
        if (a + 1 >= s->maxa) return rc_bail(s, 6);
        Yn = Yp + 3 * n;
        /* np.linalg.solve(Jbig, -R), then sum|R| (the undamped Newton reads
           it for its overflow alone: a bound, numpy's pairwise sum not
           mirrored) */
        for (k = 0; k < M; k++) {
            if (!(fabs(s->R[k]) < 0x1p1000)) return rc_bail(s, 4);
            s->b[k] = -s->R[k];
        }
        for (k = 0; k < M * M; k++) if (!isfinite(s->Jf[k])) return rc_bail(s, 4);
        ((dgesv_t) s->dgesv)(&M64, &one64, s->Jf, &M64, s->ipiv, s->b, &M64, &info64);
        /* (numpy's solve clears what LAPACK raised; an info is its LinAlgError) */
        feclearexcept(FE_ALL_EXCEPT);
        if (info64 != 0) return rc_bail(s, 5);
        for (k = 0; k < M; k++) if (!isfinite(s->b[k])) return rc_bail(s, 9);
        for (i = 0; i < 3; i++) {
            const double *yp = Yp + i * n, *d = s->b + i * m;
            double *yn = Yn + i * n, mx = 0.0;
            for (k = 0; k < r; k++) yn[k] = yp[k] + d[k];
            yn[r] = yp[r] + 0.0;
            for (k = r + 1; k < n; k++) yn[k] = yp[k] + d[k - 1];
            if (fetestexcept(RC_FLAGS)) return rc_bail(s, 3);
            if (walk) {
                e = walk(s->wF, s->wPR, s->wTI, s->T, s->wNM, s->wOFF, s->wK, yn, yp, 0, s->wn);
                if (e < s->wn) return rc_bail(s, 8);
                feclearexcept(FE_ALL_EXCEPT);
            }
            for (k = 0; k < n; k++) s->stp[k] = yn[k] - yp[k];
            if (fetestexcept(RC_FLAGS)) return rc_bail(s, 3);
            for (k = 0; k < n; k++) {
                double v;
                if (k == r) continue;
                v = fabs(s->stp[k]);
                if (!isfinite(v)) return rc_bail(s, 9);
                if (v > mx) mx = v;
            }
            if (mx > scale) scale = mx;
        }
        a++;
        e = assemble(s, a);
        if (e) return rc_bail(s, e);
        for (i = 0; i < 3; i++)
            for (k = 0; k < n; k++) {
                double v;
                if (k == r) continue;
                v = fabs(Yn[i * n + k]);
                if (v > ynorm) ynorm = v;
            }
        feclearexcept(FE_ALL_EXCEPT);
        s->tolv[0] = s->reltol * ynorm;
        s->tolv[1] = s->tolv[0] + s->abstol;
        if (fetestexcept(RC_FLAGS)) return rc_bail(s, 3);
        s->iters = it + 1;
        if (scale <= s->tolv[1]) {
            s->assemblies = a + 1;
            s->status = 1;
            return 1;
        }
    }
    return rc_bail(s, 11);
}
"""
RADAU_CDEF = """
typedef struct {
    const void *core; void *core_fn; double T;
    const double *binG, *bini; double *dummy;
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    void *dgesv;
    long n, iref, maxiter, maxa;
    double h, reltol, abstol;
    double A[9];
    const double *qn, *u, *seed;
    double *Ya, *qa, *ia, *Ca, *Ga;
    double *K, *F, *R, *Jf, *b, *blk, *stp, *tolv;
    int64_t *ipiv;
    long assemblies, iters, status;
} radau_t;
long hdl_fn(radau_t *s);
"""

#: the C's status codes past 1 (converged): why it handed the solve back
BAIL = {2: 'nonfinite', 3: 'flags', 4: 'rbound', 5: 'dgesv', 6: 'assemblies', 8: 'walkstop',
        9: 'step', 11: 'maxiter'}

_driver = None
_MOD = {}
#: numpy's solve as imported: a caller's stand-in (a spy, an injected
#: failure) takes the Python Newton
_SOLVE = np.linalg.solve


def driver():
    """`(ffi, cfn, dgesv address)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        from pycircuit.circuit import _hdl_cbackend as cb
        from pycircuit.circuit import _tran_newton_c
        nd = _tran_newton_c.driver()
        if nd is None:
            STATUS = f'off (the C Newton: {_tran_newton_c.STATUS})'
            return None
        try:
            ffi, cfn, _key, _cold, _secs = cb.load_kernel(RADAU_C, RADAU_CDEF)
        except (cb.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, nd[2]['dgesv'])
        STATUS = 'c'
    return _driver or None


def _chain():
    """The machinery `solve` stands in for, as defined (`_paths.genuine`),
    read once every piece is its module's own; None: not yet.  A call then
    compares what the Python Newton would call now (`_now`)."""
    if 'chain' not in _MOD:
        from pycircuit.circuit import _hdl_climit, _tran_core, _tran_radau
        from pycircuit.circuit.circuit import SubCircuit
        from pycircuit.circuit.transient import Transient as T
        R, CM = 'pycircuit.circuit._tran_radau', 'pycircuit.circuit._tran_companion'
        own = ((T, 'Transient', 'pycircuit.circuit.transient'),
               (T._coupled_stage_context, '_RadauStages._coupled_stage_context', R),
               (T._coupled_stage_system, '_RadauStages._coupled_stage_system', R),
               (T._stage_limiter_states, '_RadauStages._stage_limiter_states', R),
               (T._stages_converged, '_RadauStages._stages_converged', R),
               (T._limiters_at_rest, '_RadauStages._limiters_at_rest', R),
               (T._memo_put, '_CompanionModel._memo_put', CM),
               (T._source_at, '_CompanionModel._source_at', CM),
               (T._stage_source, '_SequentialStages._stage_source',
                'pycircuit.circuit._tran_stages'),
               (T._newton_limiter, '_StepNewton._newton_limiter',
                'pycircuit.circuit._tran_newton'),
               (_tran_radau.remove_row_col, 'remove_row_col', 'pycircuit.circuit.analysis'),
               (_tran_core.passes, 'passes', 'pycircuit.circuit._tran_core'),
               (SubCircuit.limit, 'SubCircuit.limit', 'pycircuit.circuit.circuit'),
               (_hdl_climit.limit_walk, 'limit_walk', 'pycircuit.circuit._hdl_climit'))
        if not all(_paths.genuine(*o) for o in own):
            return None
        _MOD['chain'] = tuple(o for o, _q, _m in own)
        _MOD['mods'] = (_tran_radau, _tran_core, _hdl_climit)
    return _MOD['chain']


def _now(T, mods, cir):
    """What the Python Newton calls now, in `_chain`'s order."""
    tr_, tc, hc = mods
    return (T._coupled_stage_context, T._coupled_stage_system, T._stage_limiter_states,
            T._stages_converged, T._limiters_at_rest, T._memo_put, T._source_at,
            T._stage_source, T._newton_limiter, tr_.remove_row_col, tc.passes,
            getattr(type(cir), 'limit', None), hc.limit_walk)


class _Ctx:
    """One transient's struct and buffers for one core, size and reference
    row (the walk's fields named as `_tran_newton_c._Ctx`'s, whose
    `_walk_ready` sets them)."""

    __slots__ = (
        'buf',
        'core',
        'ffi',
        'iref',
        'keep',
        'n',
        's',
        'walk',
        'walk_arms',
        'walk_dicts',
        'walk_els',
        'walk_held',
        'walk_hw_els',
        'walk_hwk',
        'walk_kern',
        'walk_ok',
        'walk_packs',
        'walk_stamp',
        'walk_watch',
    )

    def __init__(self, core, n, iref, ffi, cfn_core, cffi_core, dgesv):
        m = n - 1
        M = 3 * m
        A = MAXA
        self.core, self.n, self.iref, self.ffi = core, n, iref, ffi
        self.s = s = ffi.new('radau_t *')
        self.buf = buf = {
            'Ya': np.empty((A, 3, n)), 'qa': np.empty((A, 3, n)), 'ia': np.empty((A, 3, n)),
            'Ca': np.empty((A, 3, n, n)), 'Ga': np.empty((A, 3, n, n)),
            'qn': np.empty(n), 'u': np.empty((3, n)), 'seed': np.empty((3, n)),
            'K': np.empty((3, n)), 'F': np.empty((3, n)), 'R': np.empty(M),
            'Jf': np.empty(M * M), 'b': np.empty(M), 'blk': np.empty(n * n),
            'stp': np.empty(n), 'tolv': np.empty(2), 'dummy': np.zeros(max(n * n, 1)),
        }
        self.keep = []
        fb = ffi.from_buffer
        for k, arr in buf.items():
            h = fb('double *', arr)
            self.keep.append(h)
            setattr(s, k, h)
        ipiv = np.empty(max(M, 1), dtype=np.int64)
        buf['ipiv'] = ipiv
        self.keep.append(fb('int64_t *', ipiv))
        s.ipiv = self.keep[-1]
        s.dgesv = ffi.cast('void *', dgesv)
        s.core = ffi.cast('void *', int(cffi_core.cast('uintptr_t', core.cs)))
        s.core_fn = ffi.cast('void *', int(cffi_core.cast('uintptr_t', cfn_core)))
        ## the core's own bins: G and i of the pass just run
        for k in ('binG', 'bini'):
            a = core.arrays[k]
            self.keep.append(a)
            setattr(s, k, ffi.cast('double *', a.ctypes.data))
        s.n, s.iref, s.maxa = n, iref, A
        s.walk_fn = ffi.NULL
        s.wn = 0
        self.walk = self.walk_ok = None
        self.walk_els = self.walk_dicts = self.walk_packs = self.walk_kern = ()
        self.walk_hw_els = self.walk_hwk = ()
        self.walk_stamp, self.walk_arms, self.walk_held = -1, 0, 0
        self.walk_watch = ()


def _undo_source(memo, keys, had, counts):
    """The attempt's one trace, removed: the source memo's entries its
    calls made and the counts they added -- the Python Newton makes them
    again."""
    for key, h in zip(keys, had):
        if not h:
            memo.pop(key, None)
    for k, v in counts.items():
        _PC[k] = v


def solve(tr, ctx, seed, src, provided_function, lims, nobypass, reltol, abstol, maxit):
    """`_stage_newton(seed)` (undamped, no shunt) in C: the three converged
    stages as fresh full-width arrays, the device memo and the source memo
    as the Python Newton leaves them -- or None: declined or handed back,
    nothing left behind, the Python Newton runs."""
    if not ENABLED:
        return _no('radau_c:off')
    if provided_function is not None:
        ## (the Python calls it at every assembly's stage times: a caller may
        ## count it)
        return _no('radau_c:pf')
    if lims:
        return _no('radau_c:stateful')
    if not nobypass:
        return _no('radau_c:bypass')
    drv = _driver or driver()
    if not drv:
        return _no('radau_c:driver')
    ffi, cfn, dgesv = drv
    chain = _MOD.get('chain') or _chain()
    if chain is None:
        return _no('radau_c:patched')
    T = chain[0]
    if type(tr) is not T:
        return _no('radau_c:class')
    cir = tr.cir
    if _now(T, _MOD['mods'], cir) != chain[1:] or np.linalg.solve is not _SOLVE:
        return _no('radau_c:patched')
    from pycircuit.circuit import _tran_newton_c as NC
    M = NC._MOD or NC._mods()
    if type(tr.toolkit) is not M['NumericToolkit']:
        return _no('radau_c:toolkit')
    cd, td = cir.__dict__, tr.__dict__
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd or 'limit' in cd:
        return _no('radau_c:shadow_cir')
    if ('_memo_put' in td or '_source_at' in td or '_coupled_stage_system' in td
            or '_stages_converged' in td or '_stage_limiter_states' in td
            or '_limiters_at_rest' in td or '_newton_limiter' in td):
        return _no('radau_c:shadow_tr')
    if M['limiting'].CIRCUIT_LEVEL:
        return _no('radau_c:stateful')
    memo = td.get('_u_memo')
    if memo is None or not M['companion'].U_MEMO:
        return _no('radau_c:umemo')
    tc = M['core']
    if not (tc.CORE and tc.CORE_PASSES):
        return _no('radau_c:core')
    core = tc.core_for(tr)
    if core is None:
        return _no('radau_c:unservable')
    n, iref = core.n, ctx.iref
    if not (0 <= iref < n) or n < 2 or ctx.m != n - 1:
        return _no('radau_c:iref')
    h = ctx.h
    if type(h) is not float and type(h) is not np.float64:
        return _no('radau_c:h')
    if (type(reltol) is not float and type(reltol) is not int
            and type(reltol) is not np.float64) or type(reltol) is bool:
        return _no('radau_c:tol')
    qn, Amat = ctx.qn, ctx.Amat
    if not (type(qn) is np.ndarray and qn.dtype == np.float64 and qn.shape == (n,)
            and type(Amat) is np.ndarray and Amat.dtype == np.float64 and Amat.shape == (3, 3)):
        return _no('radau_c:ctx')
    Y0 = [np.array(y, dtype=float) for y in seed]
    if len(Y0) != 3 or any(y.shape != (n,) for y in Y0):
        return _no('radau_c:seed')
    has_limiter = tr._newton_limiter() is not None
    if has_limiter:
        Tw = getattr(tr.epar, 'T', 300.0)
        if type(Tw) is not float and type(Tw) is not int and type(Tw) is not np.float64:
            ## (`limit_walk` declines such a temperature: the Python loop limits)
            return _no('radau_c:T')
    rec = td.get('_radau_c')
    if rec is None or rec.core is not core or rec.n != n or rec.iref != iref:
        cffi_core, cfn_core, _dg = tc.driver()
        rec = _Ctx(core, n, iref, ffi, cfn_core, cffi_core, dgesv)
        td['_radau_c'] = rec
    s, buf = rec.s, rec.buf
    if has_limiter:
        if not NC._walk_ready(tr, rec):
            return _no('radau_c:limit')
    elif s.wn:
        s.walk_fn = ffi.NULL
        s.wn = 0
        rec.walk = rec.walk_ok = None
    Tc = core.probe(tr.epar)
    if Tc.__class__ is str:
        return _no('radau_c:ready')
    ## ---- from here the attempt leaves one trace: the step's source memo ----
    analysis = tr.par.analysis
    keys = [(tt, analysis) for tt in ctx.tstage]
    try:
        had = [key in memo for key in keys]
    except TypeError:
        return _no('radau_c:u')
    counts = {k: _PC.get(k, 0) for k in M['u_keys']}
    try:
        us = [src(tt) for tt in ctx.tstage]
    except Exception:                                          # noqa: BLE001
        ## (the call raises again in the Python Newton, where it raised before)
        _undo_source(memo, keys, had, counts)
        return _no('radau_c:u')
    if not all(type(u) is np.ndarray and u.dtype == np.float64 and u.shape == (n,)
               for u in us):
        _undo_source(memo, keys, had, counts)
        return _no('radau_c:u')
    for j in range(3):
        buf['u'][j] = us[j]
        buf['seed'][j] = Y0[j]
    np.copyto(buf['qn'], qn)
    for k, v in enumerate(Amat.ravel().tolist()):
        s.A[k] = v
    s.T = Tc
    s.h, s.reltol, s.abstol, s.maxiter = float(h), float(reltol), float(abstol), int(maxit)
    status = cfn(s)
    if status != 1:
        _undo_source(memo, keys, had, counts)
        if status == 8:
            ## (the walk stopped: its tables are taken again next time)
            rec.walk_ok = None
        return _no('radau_c:bail:' + BAIL.get(status, str(status)))
    ## ---- converged: the memo as the Python Newton leaves it ----
    A = int(s.assemblies)
    Ya, qa, ia, Ca, Ga = buf['Ya'], buf['qa'], buf['ia'], buf['Ca'], buf['Ga']
    put = tr._memo_put
    for a in range(A):
        for j in range(3):
            put(Ya[a, j], {'q': qa[a, j].copy(), 'i': ia[a, j].copy(),
                           'C': Ca[a, j].copy(), 'G': Ga[a, j].copy()})
    _PC['umemo:hit'] += 3 * (A - 1)
    _PC['radau_c:served'] += 1
    _PC['radau_c:assemblies'] += A
    return [Ya[A - 1, j].copy() for j in range(3)]
