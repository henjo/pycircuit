"""THE SHOOTING'S COUPLED STAGE STEP MAP IN C (speed round 11, item 4;
2026-10-06).

A fully implicit stage method's dense period walk (`_walk_stage`, Radau
IIA) propagates the monodromy through every step: the step's stage system
``J Z = B``, ``J[i][j] = delta_ij C(Y_i) + h A_ij G(Y_j)``, assembled
(`_stage_block`), LU-factored by the caller's dense solver (SciPy's
`lu_factor`) and solved for ``B = [C_n P]_i`` (`_StageStep.solve`, SciPy's
`getrs`) -- the last stage block of the solution the new `P`.  On the PSP
stage's radau PSS that was ~530 k instructions a step around ~50 k of
LAPACK: the readers' memo lookups and reductions, a dozen small numpy
calls assembling the block, `lu_factor`'s and `lu_solve`'s checks.  Here
ONE C call does the step map: the block assembled column-major -- as
`lu_factor` copies it for LAPACK -- from the stages' FULL `C` and `G` (the
step's device memo, the records the readers would read; the reference row
and column skipped as the readers reduce them) in `_stage_block`'s IEEE
operations (``h A_ij`` first, then its product with `G_j`, `C_i` added on
the diagonal blocks), SciPy's own `scipy_dgetrf_`, the right-hand side
stacked from numpy's ``C_n @ P`` (made in Python as `solve` makes it:
numpy's own BLAS) and SciPy's own `scipy_dgetrs_`; the new `P` is the last
stage block of the solution laid out as `getrs` returns it (Fortran
order), so the next step's product is numpy's same call.  Bit for bit the
Python step (the stage-0 prototype: 120 of 120 steps' `P`, LU and pivots;
pinned in `tests/test_pinned_pairs.py`).

ONLY THE PLAIN WALK: dense, no step kept, no period or event column -- the
shooting Newton's monodromy (`_walk_stage` asks only then).

DECLINES, counted, before anything is kept or counted otherwise
(`_paths.no`): the switch; no build or no SciPy OpenBLAS; a piece the C
stands in for not its module's own (the stage step, the block, the
readers, `_factorise`, the solvers' `factor` and choice, `_Factored`,
`lu_solve`, SciPy's `lu_factor`, `_StageStep.solve` and the solves it
calls) or shadowed on the walk or the transient; a toolkit that is not
numeric; a linear solver that is not `AutoSolver`'s dense choice or a
`DenseSolver`; a tableau that is not float64 ``s x s``, ``s <= MAXS``; a stage system
past `MAXN` unknowns;
what the readers would read otherwise -- a junction (PCNR's `G`),
limiting state kept or the trims off (the sync would run), a stage whose
reference entry is not +0.0 or that the memo lacks; and the C's own
hand-backs: a non-finite input (`_stage_block` declines and `lu_factor`
raises), a floating-point flag in the assembly (`_stage_block`'s
`errstate`), a non-finite ``C_n P`` (`lu_solve` raises), a singular block
(`lu_factor` warns), a LAPACK `info`.  The Python step then runs, as
before: the same answer, exception, warning and counts.

History: `doc/shooting_history.md`, `_sens_c.py`.
"""
import ctypes
import os

import numpy as np

from pycircuit.circuit import _limiting, _paths

from . import _pss_inner

#: `PYCIRCUIT_SENS_C=0` makes the Python step every time
SENS_C = os.environ.get('PYCIRCUIT_SENS_C', '1') != '0'
#: the most stages, and the most unknowns of the stage system, it takes:
#: its block's buffer stays with the transient between steps (the Python
#: step's is freed) -- past 600 unknowns 2.9 MB, and LAPACK's cost, not
#: the Python around it, is the step's
MAXS = 8
MAXN = 600
STATUS = 'not loaded'
_driver = None
_MOD = {}
_PC = _paths.COUNTS
_no = _paths.no
_F64 = np.dtype(np.float64)
#: +0.0's bytes: a full stage whose reference entry holds them is the
#: state the readers rebuild from the reduced one (`_stage_reads`)
_PLUS_ZERO = np.float64(0.0).tobytes()

SENS_C_SRC = r"""
#include <fenv.h>
#include <stddef.h>
#include <stdint.h>
typedef void (*dgetrf_t)(const int *, const int *, double *, const int *, int *, int *);
typedef void (*dgetrs_t)(const char *, const int *, const int *, const double *, const int *,
                         const int *, double *, const int *, int *, size_t);
typedef struct {
    void *dgetrf, *dgetrs;
    long s, m, k, nf, iref;
    double h;
    const double *A, *Cf, *Gf, *base;
    double *J; int *ipiv; double *Z;
    long status, info;
} sens_t;
#define SENS_FLAGS (FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO)
/* the full row (column) of reduced index r: the reference one skipped */
#define FULL(r) ((r) < t->iref ? (r) : (r) + 1)
long hdl_fn(sens_t *t)
{
    long s = t->s, m = t->m, k = t->k, nf = t->nf, n = s * m, nn = nf * nf;
    long i, j, r, c;
    double h = t->h;
    /* every reduced entry finite: `_stage_block` declines otherwise, and the
       loop's block makes `lu_factor` raise */
    if (!isfinite(h)) return t->status = 2;
    for (i = 0; i < s * s; i++) if (!isfinite(t->A[i])) return t->status = 2;
    for (i = 0; i < s; i++)
        for (r = 0; r < m; r++)
            for (c = 0; c < m; c++) {
                long f = i * nn + FULL(r) * nf + FULL(c);
                if (!isfinite(t->Cf[f]) || !isfinite(t->Gf[f])) return t->status = 2;
            }
    /* C_n P finite: `lu_solve` refuses it otherwise */
    for (i = 0; i < m * k; i++) if (!isfinite(t->base[i])) return t->status = 4;
    /* the block, column-major: J[(j m + c) n + i m + r] = [i == j] C_i[r, c]
       + (h A_ij) G_j[r, c] -- `_stage_block`'s operations, under its errstate */
    feclearexcept(FE_ALL_EXCEPT);
    for (i = 0; i < s; i++)
        for (j = 0; j < s; j++) {
            double hA = h * t->A[i * s + j];
            const double *G = t->Gf + j * nn, *C = t->Cf + i * nn;
            for (c = 0; c < m; c++) {
                double *col = t->J + (j * m + c) * n + i * m;
                long fc = FULL(c);
                for (r = 0; r < m; r++) {
                    long f = FULL(r) * nf + fc;
                    double v = hA * G[f];
                    if (i == j) v = C[f] + v;
                    col[r] = v;
                }
            }
        }
    if (fetestexcept(SENS_FLAGS)) return t->status = 3;
    int n32 = (int) n, k32 = (int) k, info = 0;
    ((dgetrf_t) t->dgetrf)(&n32, &n32, t->J, &n32, t->ipiv, &info);
    t->info = info;
    if (info != 0) return t->status = 5;
    /* the right-hand side, column-major: every stage's block C_n P */
    for (c = 0; c < k; c++)
        for (i = 0; i < s; i++)
            for (r = 0; r < m; r++)
                t->Z[c * n + i * m + r] = t->base[r * k + c];
    char trans = 'N';
    ((dgetrs_t) t->dgetrs)(&trans, &n32, &k32, t->J, &n32, t->ipiv, t->Z, &n32, &info, 1);
    t->info = info;
    if (info != 0) return t->status = 6;
    return t->status = 1;
}
"""

SENS_CDEF = """
typedef struct {
    void *dgetrf, *dgetrs;
    long s, m, k, nf, iref;
    double h;
    const double *A, *Cf, *Gf, *base;
    double *J; int *ipiv; double *Z;
    long status, info;
} sens_t;
long hdl_fn(sens_t *t);
"""

#: the C's status codes past 1 (served): why it handed the step back
BAIL = {2: 'nonfinite', 3: 'flags', 4: 'rhs', 5: 'singular', 6: 'lapack'}


def driver():
    """`(ffi, cfn, dgetrf, dgetrs)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            import scipy
            import scipy.linalg  # (its OpenBLAS mapped: `lu_factor`'s)

            from pycircuit.circuit._tran_newton_c import _openblas
            spl = _openblas(scipy, 'libscipy_openblas-*.so')
            if spl is None:
                STATUS = "off (SciPy's own OpenBLAS not found)"
                return None
            getrf = ctypes.cast(spl.scipy_dgetrf_, ctypes.c_void_p).value
            getrs = ctypes.cast(spl.scipy_dgetrs_, ctypes.c_void_p).value
        except (ImportError, OSError, AttributeError) as e:
            STATUS = f'off ({e})'
            return None
        from pycircuit.circuit import _hdl_cbackend as cb
        try:
            ffi, cfn, _key, _cold, _secs = cb.load_kernel(SENS_C_SRC, SENS_CDEF)
        except (cb.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, getrf, getrs)
        STATUS = 'c'
    return _driver or None


def _scipy_own(f, name, mod):
    """`f` is SciPy's `name` as `mod` defines it -- or SciPy's own batch
    wrapper around it (`scipy._lib._util._apply_over_batch`, SciPy 1.16+),
    which calls it unchanged on an array of its core rank."""
    if _paths.genuine(f, name, mod):
        return True
    w = getattr(f, '__wrapped__', None)
    code = getattr(f, '__code__', None)
    return (w is not None and code is not None and _paths.genuine(w, name, mod)
            and code.co_name == 'wrapper'
            and code.co_filename.replace(os.sep, '/').endswith('scipy/_lib/_util.py'))


def _pieces():
    """What the Python step calls, as their modules define them -- read once
    every one is (`_paths.genuine`); None: one is a stand-in, read again at
    the next call."""
    g = _MOD.get('g')
    if g is None:
        import scipy.linalg

        from pycircuit.circuit import linearsolver as LS
        from pycircuit.circuit.toolkit import NumericToolkit

        from . import _numerics, _pss_walks, _steps
        IT, PW = _pss_inner._InnerTransient, _pss_walks._PeriodWalks
        M_IN, M_W = 'pycircuit.circuit.shooting._pss_inner', 'pycircuit.circuit.shooting._pss_walks'
        M_ST, M_NU = 'pycircuit.circuit.shooting._steps', 'pycircuit.circuit.shooting._numerics'
        M_LS = 'pycircuit.circuit.linearsolver'
        pieces = ((PW._stage_step, '_PeriodWalks._stage_step', M_W),
                  (_pss_walks._stage_block, '_stage_block', M_W),
                  (_pss_walks._stage_reads, '_stage_reads', M_W),
                  (IT._factorise, '_InnerTransient._factorise', M_IN),
                  (_steps._StageStep.solve, '_StageStep.solve', M_ST),
                  (_numerics._lu_solve_split, '_lu_solve_split', M_NU),
                  (_numerics._complex_solve, '_complex_solve', M_NU),
                  (LS.AutoSolver.factor, 'AutoSolver.factor', M_LS),
                  (LS.AutoSolver._select, 'AutoSolver._select', M_LS),
                  (LS.DenseSolver.factor, 'DenseSolver.factor', M_LS),
                  (LS._Factored.solve, '_Factored.solve', M_LS),
                  (LS.lu_solve, 'lu_solve', M_LS))
        if (not all(_paths.genuine(*p) for p in pieces)
                or not _scipy_own(scipy.linalg.lu_factor, 'lu_factor', 'scipy.linalg._decomp_lu')):
            return None
        pieces += ((scipy.linalg.lu_factor, 'lu_factor', 'scipy.linalg._decomp_lu'),)
        g = tuple(p[0] for p in pieces)
        _MOD.update(g=g, mods=(_pss_walks, _steps, _numerics, LS, scipy.linalg),
                    types=(LS.AutoSolver, LS.DenseSolver, NumericToolkit))
    return g


def _now():
    """What the Python step calls now, in `_pieces`' order."""
    W, ST, NU, LS, SL = _MOD['mods']
    return (W._PeriodWalks._stage_step, W._stage_block, W._stage_reads,
            _pss_inner._InnerTransient._factorise, ST._StageStep.solve,
            ST._lu_solve_split, NU._complex_solve, LS.AutoSolver.factor,
            LS.AutoSolver._select, LS.DenseSolver.factor, LS._Factored.solve,
            LS.lu_solve, SL.lu_factor)


#: the walk's and the transient's names an instance shadow of would send the
#: Python step elsewhere
_WALK_NAMES = ('_stage_step', '_factorise', '_C_at', '_G_at', '_sync_limit_at',
               '_insert_refnode', '_get_linearsolver', '_pcnr_junctions', '_transient')


class _Ctx:
    """One inner transient's struct and buffers for one size: the tableau, the stages'
    full `C` and `G`, ``C_n P`` and the block (its LU) and pivots -- the
    solution's array is fresh every step (it becomes `P`)."""

    __slots__ = ('A', 'A_bytes', 'Cf', 'Gf', 'J', 'base', 'ffi', 'ipiv', 'keep', 'key', 't')

    def __init__(self, ffi, getrf, getrs, s, nf, k, iref):
        m = nf - 1
        n = s * m
        self.key = (s, nf, k, iref)
        self.ffi = ffi
        self.A = np.empty((s, s))
        self.A_bytes = None
        self.Cf = np.empty((s, nf, nf))
        self.Gf = np.empty((s, nf, nf))
        self.base = np.empty((m, k))
        self.J = np.empty(n * n)
        self.ipiv = np.empty(n, dtype=np.int32)
        t = self.t = ffi.new('sens_t *')
        t.dgetrf = ffi.cast('void *', getrf)
        t.dgetrs = ffi.cast('void *', getrs)
        t.s, t.m, t.k, t.nf, t.iref = s, m, k, nf, iref
        self.keep = (ffi.from_buffer('double *', self.A), ffi.from_buffer('double *', self.Cf),
                     ffi.from_buffer('double *', self.Gf), ffi.from_buffer('double *', self.base),
                     ffi.from_buffer('double *', self.J), ffi.from_buffer('int *', self.ipiv))
        t.A, t.Cf, t.Gf, t.base, t.J, t.ipiv = self.keep


def step(walks, xn, h, tab, P):
    """The plain dense walk's step map from `xn` through the step just
    taken -- the new monodromy `P`, as `_stage_step` and `_StageStep.solve`
    make it -- or None (counted): the Python step runs."""
    if not SENS_C:
        return _no('pss.sens:off')
    drv = _driver or driver()
    if not drv:
        return _no('pss.sens:driver')
    g = _MOD.get('g') or _pieces()
    if g is None:
        return _no('pss.sens:patched')
    W = _MOD['mods'][0]
    rd = W._READERS.get('g') or W._readers()
    T, d = type(walks), walks.__dict__
    if (rd is None or _now() != g
            or (T._stage_step, T._factorise) != (g[0], g[3])
            or (T._C_at, T._G_at, T._sync_limit_at, T._insert_refnode) != rd[:4]
            or any(name in d for name in _WALK_NAMES)):
        return _no('pss.sens:patched')
    AS, DS, NTK = _MOD['types']
    tk = walks.toolkit
    if type(tk) is not NTK:
        return _no('pss.sens:toolkit')
    ## (the readers' own conditions: `_stage_reads`)
    if not _pss_inner.SHOOT_TRIM or _limiting.CIRCUIT_LEVEL:
        return _no('pss.sens:limits')
    if walks._pcnr_junctions():
        return _no('pss.sens:junctions')
    tr = walks._transient()
    td = tr.__dict__
    if type(tr)._memo_get is not rd[4] or '_memo_get' in td:
        return _no('pss.sens:patched')
    lims = td.get('_stateful_lims')
    if lims is None:
        lims = _limiting.stateful_limiters(tr.cir)
    if lims:
        return _no('pss.sens:limits')
    ## the caller's solver: `AutoSolver`'s dense choice, or a `DenseSolver`
    ## (`_factorise`: its `factor` is SciPy's `lu_factor`)
    ls = walks._get_linearsolver()
    if type(ls) is AS:
        dense = ls._choice
        if (dense is None or type(dense) is not DS or 'factor' in ls.__dict__
                or '_select' in ls.__dict__ or 'factor' in dense.__dict__):
            return _no('pss.sens:solver')
    elif type(ls) is not DS or 'factor' in ls.__dict__:
        return _no('pss.sens:solver')
    A = tab[0]
    if not (type(A) is np.ndarray and (A.dtype is _F64 or A.dtype == np.float64)
            and A.ndim == 2 and A.shape[0] == A.shape[1] and 0 < A.shape[0] <= MAXS):
        return _no('pss.sens:tableau')
    s = A.shape[0]
    if type(h) is not float and type(h) is not np.float64:
        return _no('pss.sens:tableau')
    ## the stages and `x_n` as the readers read them: one memo record each
    iref = walks.irefnode
    Y = tr.last_step.Y
    if len(Y) != s:
        return _no('pss.sens:state')
    recs = []
    for yf in Y:
        if (type(yf) is not np.ndarray or yf.dtype is not _F64 or yf.ndim != 1
                or not 0 <= iref < yf.shape[0] or yf[iref:iref + 1].tobytes() != _PLUS_ZERO):
            return _no('pss.sens:state')
        rec = tr._memo_get(yf)
        if rec is None or 'C' not in rec or 'G' not in rec:
            return _no('pss.sens:memo')
        recs.append(rec)
    nf = Y[0].shape[0]
    if s * (nf - 1) > MAXN:
        return _no('pss.sens:size')
    recn = tr._memo_get(walks._insert_refnode(xn))
    if recn is None or 'C' not in recn:
        return _no('pss.sens:memo')
    for M in (recn['C'], *(r['C'] for r in recs), *(r['G'] for r in recs)):
        if type(M) is not np.ndarray or M.dtype is not _F64 or M.shape != (nf, nf):
            return _no('pss.sens:kind')
    ## `C_n P` as `_StageStep.solve` makes it, from `_C_at(x_n)`'s reduction
    try:
        Cn = np.asarray(_pss_inner.remove_row_col((recn['C'],), iref, tk)[0])
        base = Cn @ P
    except Exception:                                          # noqa: BLE001
        ## (the Python step raises it again, where it raised before)
        return _no('pss.sens:kind')
    m = nf - 1
    if not (type(base) is np.ndarray and base.dtype is _F64 and base.ndim == 2
            and base.shape[0] == m):
        return _no('pss.sens:kind')
    out = _map(td, drv, A, [r['C'] for r in recs], [r['G'] for r in recs], iref, h, base)
    if type(out) is int:
        return _no('pss.sens:' + BAIL.get(out, str(out)))
    ## the readers' sync, skipped as `_G_at` skips it, counted as it counts
    _PC['pss.sync:skipped'] += s
    _PC['pss.sens:served'] += 1
    return out


def _map(td, drv, A, Cfs, Gfs, iref, h, base):
    """The C step map on the stages' FULL `C` and `G` (reference row and
    column `iref`): the new `P`, the last stage block of a fresh solution
    in `getrs`' Fortran order -- or the C's status (an int: `BAIL`).  The
    struct and buffers are the inner transient's (`td`), kept per size."""
    s, nf, k = A.shape[0], Cfs[0].shape[0], base.shape[1]
    m = nf - 1
    ctx = td.get('_sens_c')
    if ctx is None or ctx.key != (s, nf, k, iref):
        ffi, _cfn, getrf, getrs = drv
        ctx = _Ctx(ffi, getrf, getrs, s, nf, k, iref)
        td['_sens_c'] = ctx
    ab = A.tobytes()
    if ab != ctx.A_bytes:
        np.copyto(ctx.A, A)
        ctx.A_bytes = ab
    for i in range(s):
        np.copyto(ctx.Cf[i], Cfs[i])
        np.copyto(ctx.Gf[i], Gfs[i])
    np.copyto(ctx.base, base)
    ## (the solution in a fresh array, `getrs`' Fortran order: its last
    ## stage block, a view of it, is the new `P`)
    Zc = np.empty((k, s * m))
    t = ctx.t
    t.Z = ctx.ffi.from_buffer('double *', Zc)
    t.h = float(h)
    status = drv[1](t)
    t.Z = ctx.ffi.NULL
    if status != 1:
        return int(status)
    return Zc.T[(s - 1) * m:s * m]
