"""THE DENSE MULTISTEP WALK'S MONODROMY STEP IN C: `_steps._LiveMonodromy.step`
on its live columns -- `_lmm_recursion` with each capacitance read through
its nonzero columns (`_ColumnsOf`) and the step's solve through the
caller's dense solver, numpy's `linalg.solve` -- in ONE C call, bit for bit
the Python step (pinned in `tests/test_pinned_pairs.py`).

OWNS.  On the walk object (the PSS), `_lmmsens_c` (`_Ctx`: the
`lmmsens_t` struct and its work buffers, kept per size and widest
monodromy).  `_paths`: `pss.lmm:` -- `pss.lmm:<why>` a decline or a
hand-back (`BAIL`), `pss.lmm:served`.  Switch: `PYCIRCUIT_LMMSENS_C=0`
(`LMMSENS_C`, read on every call); `STATUS` says why the path is off where
it is.

THE STEP.  ``S = sum_k a_k C_{n-k} P_{n-k} (+ b Pq)``, ``P_n = -Jf^-1 S``,
``Pq_n = a_0 C_n P_n + S`` on the live columns.  Each product is numpy's
matmul of the gathered arrays -- a capacitance's nonzero columns and the
monodromy's rows they read, both C-contiguous -- which is numpy's own
`cblas_dgemm` (row-major, neither transposed) from two columns, its own
loop for one, zeros for none; the solve is numpy's own `dgesv` on the
column-major copies `linalg.solve` makes; every elementwise operation is
the Python step's, in its order.  Pn and Pqn are fresh arrays; the walk's
`_LiveMonodromy.advance` keeps them as the Python step keeps its own.

ONLY THE DENSE WALK WITH NO STEP KEPT: a kept step solves through its own
factorisation (`_factorise`), which is not numpy's solve.

DECLINES (`pss.lmm:<why>`), with nothing written that the Python step
reads: the switch; no build or no numpy OpenBLAS; a piece the C stands in
for not its module's own (the step, the recursion, `_ColumnsOf`, the
solvers' `solve` and choice, the toolkit's `linearsolver`, numpy's
`linalg.solve`); a toolkit that is not numeric; a linear solver that is not
`AutoSolver`'s dense choice or a `DenseSolver`; a walk carrying every
column (`LIVE_COLUMNS` off, or after a non-finite step); coefficients that
are not two or three real numbers; a size outside ``2 <= m <= MAXM``, fewer
than two live columns; an array not float64 C-contiguous of its size; and
the C's own hand-backs -- a singular `Jf` (numpy's solve raises), a result
not finite (the Python step's full-width fallback).  The Python step then
runs, as before: the same answer, exception and counts.

History: `doc/shooting_history.md`, `_lmmsens_c.py`.
"""
import os

import numpy as np

from pycircuit.circuit import _cbuild, _paths

#: `PYCIRCUIT_LMMSENS_C=0` makes the Python step every time
LMMSENS_C = os.environ.get('PYCIRCUIT_LMMSENS_C', '1') != '0'
#: the most unknowns it takes: its dense work buffers stay with the walk
MAXM = 2000
STATUS = 'not loaded'
_driver = None
_MOD = {}
_PC = _paths.COUNTS
_no = _paths.no
_F64 = np.dtype(np.float64)
_REAL = (float, int, np.floating)

LMMSENS_C_SRC = _cbuild.csrc('lmmsens.c')

LMMSENS_CDEF = _cbuild.csrc('lmmsens.h')

#: the C's results, named in lmmsens.h (`_CODES`); `BAIL`: {code: why} of
#: the codes past 1 (served) -- why it handed the step back
_CODES = _cbuild.constants(LMMSENS_CDEF)
BAIL = _cbuild.bails(_CODES, 'LMMSENS_BAIL_')


def driver():
    """`(ffi, cfn, dgesv, dgemm)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            npl = _cbuild.openblas(np)
            if npl is None:
                STATUS = "off (numpy's own OpenBLAS not found)"
                return None
            gesv = _cbuild.address(npl.scipy_dgesv_64_)
            gemm = _cbuild.address(npl.scipy_cblas_dgemm64_)
        except (OSError, AttributeError) as e:
            STATUS = f'off ({e})'
            return None
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(LMMSENS_C_SRC, LMMSENS_CDEF,
                                                            'pcx_lmmsens')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, gesv, gemm)
        STATUS = 'c'
    return _driver or None


def _pieces():
    """What the Python step calls, as their modules define them -- read once
    every one is (`_paths.genuine`); None: one is a stand-in, read again at
    the next call."""
    g = _MOD.get('g')
    if g is None:
        from pycircuit.circuit import _numeric
        from pycircuit.circuit import linearsolver as LS
        from pycircuit.circuit.toolkit import NumericToolkit

        from . import _steps
        M_ST, M_LS = 'pycircuit.circuit.shooting._steps', 'pycircuit.circuit.linearsolver'
        pieces = ((_steps._LiveMonodromy.step, '_LiveMonodromy.step', M_ST),
                  (_steps._LiveMonodromy._read, '_LiveMonodromy._read', M_ST),
                  (_steps._ColumnsOf.__init__, '_ColumnsOf.__init__', M_ST),
                  (_steps._ColumnsOf.__matmul__, '_ColumnsOf.__matmul__', M_ST),
                  (_steps._lmm_recursion, '_lmm_recursion', M_ST),
                  (LS.AutoSolver.solve, 'AutoSolver.solve', M_LS),
                  (LS.AutoSolver._select, 'AutoSolver._select', M_LS),
                  (LS.DenseSolver.solve, 'DenseSolver.solve', M_LS),
                  (_numeric.linearsolver, 'linearsolver', 'pycircuit.circuit._numeric'))
        if not (all(_paths.genuine(*p) for p in pieces) and _numpy_solve(np.linalg.solve)):
            return None
        g = tuple(p[0] for p in pieces) + (np.linalg.solve,)
        _MOD.update(g=g, mods=(_steps, LS, _numeric),
                    types=(LS.AutoSolver, LS.DenseSolver, NumericToolkit))
    return g


def _numpy_solve(f):
    """`f` is numpy's own `linalg.solve`: its array-function dispatcher (a
    type of numpy's, wrapping the function) -- not a caller's function, nor
    a wrapper made like one (`functools.wraps` copies the names, not the
    type)."""
    return (type(f).__module__ == 'numpy'
            and getattr(f, '__qualname__', None) == 'solve'
            and str(getattr(f, '__module__', '')).startswith('numpy.linalg'))


def _now():
    """What the Python step calls now, in `_pieces`' order."""
    ST, LS, NU = _MOD['mods']
    return (ST._LiveMonodromy.step, ST._LiveMonodromy._read, ST._ColumnsOf.__init__,
            ST._ColumnsOf.__matmul__, ST._lmm_recursion, LS.AutoSolver.solve,
            LS.AutoSolver._select, LS.DenseSolver.solve, NU.linearsolver, np.linalg.solve)


class _Ctx:
    """One walk's struct and work buffers for one size and widest monodromy;
    `last`, the last call's arrays by role with their pointers -- this
    call's inputs are that call's arrays moved on one place (the ring), so
    their pointers are kept, and so are the capacitances' column lists
    (`reuse`)."""

    __slots__ = ('A', 'Ck', 'Pk', 'S', 'T', 'X', 'ffi', 'ipiv', 'kcs', 'keep', 'last', 'm', 'ok',
                 't', 'wcap')

    def __init__(self, ffi, gesv, gemm, m, wcap):
        self.ffi, self.m, self.wcap = ffi, m, wcap
        self.Ck = np.empty(m * m)
        self.A = np.empty(m * m)
        self.Pk = np.empty(m * wcap)
        self.T = np.empty(m * wcap)
        self.S = np.empty(m * wcap)
        self.X = np.empty(m * wcap)
        self.ipiv = np.empty(m, dtype=np.int64)
        self.kcs = np.empty(3 * m, dtype=np.int64)
        self.last = {}
        self.ok = None
        t = self.t = ffi.new('lmmsens_t *')
        t.dgesv = ffi.cast('void *', gesv)
        t.dgemm = ffi.cast('void *', gemm)
        t.m, t.wcap = m, wcap
        fb = ffi.from_buffer
        self.keep = (fb('double *', self.Ck), fb('double *', self.Pk), fb('double *', self.T),
                     fb('double *', self.S), fb('double *', self.A), fb('double *', self.X),
                     fb('int64_t *', self.ipiv), fb('int64_t *', self.kcs))
        t.Ck, t.Pk, t.T, t.S, t.A, t.X, t.ipiv, t.kcs = self.keep


def _array(a, shape):
    return (type(a) is np.ndarray and a.dtype is _F64 and a.shape == shape
            and a.flags.c_contiguous)


def _seen(last, role, a, shape):
    """`a` is the last call's array in `role` (made or checked there), or
    passes `_array`."""
    got = last.get(role)
    return (got is not None and got[0] is a) or _array(a, shape)


def _which(tk, ls, g):
    """Why the toolkit or the solver is not the one the C stands in for
    (a `pss.lmm:` key), or ''."""
    AS, DS, NTK = _MOD['types']
    if type(tk) is not NTK or tk.linearsolver is not g[8]:
        return 'pss.lmm:toolkit'
    if type(ls) is AS:
        dense = ls._choice
        if (dense is None or type(dense) is not DS or 'solve' in ls.__dict__
                or '_select' in ls.__dict__ or 'solve' in dense.__dict__):
            return 'pss.lmm:solver'
    elif type(ls) is not DS or 'solve' in ls.__dict__:
        return 'pss.lmm:solver'
    return ''


def step(walks, mono, Cs, C_new, alphas, b, Jf):
    """`mono.step(Cs, C_new, alphas, b, solve)` for a dense walk keeping no
    step, `solve` the caller's solver on `Jf` -- its ``(P_new, Pq_new)``, or
    None (counted): the Python step runs."""
    if not LMMSENS_C:
        return _no('pss.lmm:off')
    drv = _driver or driver()
    if not drv:
        return _no('pss.lmm:driver')
    g = _MOD.get('g') or _pieces()
    if g is None or _now() != g:
        return _no('pss.lmm:patched')
    wd = walks.__dict__
    tk = walks.toolkit
    ls = walks._get_linearsolver()
    ctx = wd.get('_lmmsens_c')
    ## (the toolkit and solver as the last call found them: the same
    ## objects, the solver's choice unchanged -- else checked again)
    if ctx is None or ctx.ok != (tk, ls, getattr(ls, '_choice', None), g):
        why = _which(tk, ls, g)
        if why:
            return _no(why)
        if ctx is not None:
            ctx.ok = (tk, ls, getattr(ls, '_choice', None), g)
    if mono.off:
        return _no('pss.lmm:every_column')
    L = len(alphas)
    if not (2 <= L <= 3 and all(isinstance(a, _REAL) for a in alphas)
            and isinstance(b, _REAL)):
        return _no('pss.lmm:coeffs')
    P0, P1 = mono.P
    Pq = mono.Pq
    m, w = np.shape(P0) if np.ndim(P0) == 2 else (0, 0)
    if not (2 <= m <= MAXM and w >= 2):
        return _no('pss.lmm:size')
    ## (an input the last call made or read is that call's array, moved on
    ## one place: its kind and layout are known)
    last = ctx.last if ctx is not None and ctx.m == m else {}
    if not (_seen(last, 'Pn', P0, (m, w)) and _seen(last, 'P0', P1, (m, w))
            and _seen(last, 'Pqn', Pq, (m, w)) and _seen(last, 'Cn', Cs[0], (m, m))
            and _seen(last, 'C1', Cs[1], (m, m)) and _array(C_new, (m, m))
            and type(Jf) is np.ndarray and Jf.dtype is _F64 and Jf.shape == (m, m)):
        return _no('pss.lmm:layout')
    out = _call(wd, drv, m, w, alphas, b, P0, P1, Pq, Cs, C_new,
                Jf if Jf.flags.c_contiguous else np.ascontiguousarray(Jf))
    if type(out) is int:
        return _no('pss.lmm:' + BAIL.get(out, str(out)))
    _PC['pss.lmm:served'] += 1
    return out


def _call(wd, drv, m, w, alphas, b, P0, P1, Pq, Cs, C_new, Jf):
    """The C step on the walk's struct (`wd`'s, kept per size): fresh
    ``(Pn, Pqn)``, or the C's status (an int: `BAIL`)."""
    ffi, cfn, gesv, gemm = drv
    ctx = wd.get('_lmmsens_c')
    if ctx is None or ctx.m != m or ctx.wcap < w:
        ctx = _Ctx(ffi, gesv, gemm, m, max(w, ctx.wcap if ctx is not None and ctx.m == m else 0))
        wd['_lmmsens_c'] = ctx
    t = ctx.t
    t.w, t.L = w, len(alphas)
    for k, a in enumerate(alphas):
        t.alpha[k] = float(a)
    t.b = float(b)
    Pn = np.empty((m, w))
    Pqn = np.empty((m, w))
    fb = ffi.from_buffer
    last = ctx.last

    def ptr(a, role):
        got = last.get(role)
        return got[1] if got is not None and got[0] is a else fb('double *', a)
    C1, C2 = Cs[0], Cs[1]
    p = (ptr(P0, 'Pn'), ptr(P1, 'P0'), ptr(Pq, 'Pqn'), ptr(C1, 'Cn'), ptr(C2, 'C1'),
         fb('double *', C_new), fb('double *', Pn), fb('double *', Pqn))
    ## (a capacitance's column list from the last call where it is the same
    ## array, moved on one place in the ring: the walk never writes one)
    old1, oldn = last.get('C1'), last.get('Cn')
    t.reuse = ((1 if oldn is not None and oldn[0] is C1 else 0)
               | (2 if old1 is not None and old1[0] is C2 else 0))
    t.P0, t.P1, t.Pq, t.C1, t.C2, t.Cn, t.Pn, t.Pqn = p
    t.Jf = fb('double *', Jf)
    status = cfn(t)
    null = ffi.NULL
    t.P0 = t.P1 = t.Pq = t.C1 = t.C2 = t.Cn = t.Jf = t.Pn = t.Pqn = null
    ## (the C made this call's column lists whatever it returned)
    ctx.last = {'Pn': (Pn, p[6]), 'P0': (P0, p[0]), 'Pqn': (Pqn, p[7]), 'Cn': (C_new, p[5]),
                'C1': (C1, p[3])}
    if status != _CODES['LMMSENS_DONE']:
        return int(status)
    return Pn, Pqn
