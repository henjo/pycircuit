"""THE SHOOTING'S COUPLED STAGE STEP MAP ON THE WALK'S KLU IN C: a plain
walk's step on its StructuralKLU (`_PeriodWalks._stage_step` with the
walk's numeric, `StructuralKLU.walk`, and `_StageStep.solve`, as
`_walk_stage` calls them for a fully implicit stage method) in ONE C call
(`csrc/sensklu.c`), bit for bit the Python step (pinned in
`tests/test_pinned_pairs.py`).

THE STEP.  The block is the C stage Newton's last assembly at the stages it
returned (`_tran_radau_c.final_block`: the very values, in the KLU's
order); the C makes `_KLUWalk.factor_values` of it -- the walk's numeric
refactored and checked by its probe's residual, a fresh factor where the
refactor or the check fails, the first block factored -- then the stage
system's right-hand side, ``C_n P`` once a stage (``C_n P`` numpy's
product, made here as `_StageStep.solve` makes it), copied column-major
as `_KLUNumeric._run` copies it, and klu_solve of its columns at once.
The new `P` is the last stage block of that solution, laid out as the
Python's.  The walk's numeric, generation and counts (and the solver's
`numerics`) are left as `factor_values` leaves them; a fresh factor that
fails raises its LinAlgError, a solve that fails `_run`'s.

OWNS.  On each StructuralKLU it serves, `_sensklu_c` (a `_Ctx`: the struct,
the probe's buffers, the binding).  `_paths`: `pss.sensklu:served`,
`pss.sensklu:<why>` a decline or a hand-back (`BAIL`).  Switch:
`PYCIRCUIT_SENSKLU_C=0` (`SENSKLU_C`, read on every call); `STATUS` says
why the path is off where it is.

DECLINES (nothing written that the Python step reads): the switch; no
build; a piece the C stands in for not its module's own or shadowed on the
walk; a monodromy that is not a float64 `m x k` array (k >= 1); a KLU that
is not the walk's for this block or not analysed; no C last assembly for
these stages (the Python step factors `_stage_block`'s); a flag in the
check's arithmetic that numpy would warn of (`BAIL`, handed back after a
refactor in place, which the Python's repeat makes again).

History: `doc/shooting_history.md`, `_sensklu_c.py`.
"""
import os

import numpy as np

from pycircuit.circuit import _cbuild, _paths

#: `PYCIRCUIT_SENSKLU_C=0` makes the Python step every time
SENSKLU_C = os.environ.get('PYCIRCUIT_SENSKLU_C', '1') != '0'
STATUS = 'not loaded'
_driver = None
_MOD = {}
_PC = _paths.COUNTS
_no = _paths.no
_F64 = np.dtype(np.float64)

SENSKLU_SRC = _cbuild.csrc('sensklu.c')

SENSKLU_CDEF = _cbuild.csrc('sensklu.h')

#: the C's results, named in sensklu.h (`_CODES`); `BAIL`: {code: why} of
#: the hand-backs
_CODES = _cbuild.constants(SENSKLU_CDEF)
BAIL = _cbuild.bails(_CODES, 'SENSKLU_BAIL_')


def driver():
    """`(ffi, cfn)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(SENSKLU_SRC, SENSKLU_CDEF,
                                                             'pcx_sensklu')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn)
        STATUS = 'c'
    return _driver or None


def _pieces():
    """What the Python step calls, as their modules define them -- read once
    every one is (`_paths.genuine`); None: one is a stand-in, read again at
    the next call."""
    g = _MOD.get('g')
    if g is None:
        from pycircuit.circuit import _tran_radau_c
        from pycircuit.circuit import linearsolver as LS

        from . import _numerics, _pss_inner, _pss_walks, _steps
        PW, IT = _pss_walks._PeriodWalks, _pss_inner._InnerTransient
        M_W, M_IN = 'pycircuit.circuit.shooting._pss_walks', 'pycircuit.circuit.shooting._pss_inner'
        M_ST, M_NU = 'pycircuit.circuit.shooting._steps', 'pycircuit.circuit.shooting._numerics'
        M_LS = 'pycircuit.circuit.linearsolver'
        pieces = ((PW._stage_step, '_PeriodWalks._stage_step', M_W),
                  (_pss_walks._final_lu, '_final_lu', M_W),
                  (_pss_walks._own_factorise, '_own_factorise', M_W),
                  (IT._C_at, '_InnerTransient._C_at', M_IN),
                  (_tran_radau_c.final_block, 'final_block', 'pycircuit.circuit._tran_radau_c'),
                  (LS._KLUWalk.factor_values, '_KLUWalk.factor_values', M_LS),
                  (LS._KLUWalk._holds, '_KLUWalk._holds', M_LS),
                  (LS._KLUBorrowed._run, '_KLUBorrowed._run', M_LS),
                  (LS._KLUNumeric._run, '_KLUNumeric._run', M_LS),
                  (LS._Factored.solve, '_Factored.solve', M_LS),
                  (LS.StructuralKLU._factor, 'StructuralKLU._factor', M_LS),
                  (LS.StructuralKLU._free, 'StructuralKLU._free', M_LS),
                  (_steps._StageStep.solve, '_StageStep.solve', M_ST),
                  (_numerics._lu_solve_split, '_lu_solve_split', M_NU),
                  (_numerics._complex_solve, '_complex_solve', M_NU))
        if not all(_paths.genuine(*p) for p in pieces):
            return None
        g = tuple(p[0] for p in pieces)
        _MOD.update(g=g, mods=(_pss_walks, _pss_inner, _tran_radau_c, LS, _steps, _numerics))
    return g


def _now():
    """What the Python step calls now, in `_pieces`' order."""
    W, IN, RC, LS, ST, NU = _MOD['mods']
    return (W._PeriodWalks._stage_step, W._final_lu, W._own_factorise,
            IN._InnerTransient._C_at, RC.final_block, LS._KLUWalk.factor_values,
            LS._KLUWalk._holds, LS._KLUBorrowed._run, LS._KLUNumeric._run, LS._Factored.solve,
            LS.StructuralKLU._factor, LS.StructuralKLU._free, ST._StageStep.solve,
            NU._lu_solve_split, NU._complex_solve)


#: the walk's names an instance shadow of would send the Python step elsewhere
_WALK_NAMES = ('_stage_step', '_factorise', '_C_at', '_transient', '_step_klu')


class _Ctx:
    """One StructuralKLU's struct, bound to its pattern, analysis, common
    and KLU's entries, and the probe's buffers."""

    __slots__ = ('cols', 'ffi', 'keep', 'pb', 'pr', 'px', 'sym', 't')

    def __init__(self, ffi, sk):
        t = self.t = ffi.new('sensklu_t *')
        self.ffi = ffi
        M = sk.m
        self.cols = np.repeat(np.arange(M), np.diff(sk.Ap)).astype(np.int64)
        self.pb, self.px, self.pr = np.empty(M), np.empty(M), np.empty(M)
        fb = ffi.from_buffer
        self.keep = (fb('int *', sk.Ap), fb('int *', sk.Ai), fb('int64_t *', self.cols),
                     fb('double *', self.pb), fb('double *', self.px), fb('double *', self.pr))
        t.kAp, t.kAi, t.cols, t.pb, t.px, t.pr = self.keep
        lib = sk._lib
        t.klu_factor, t.klu_refactor, t.klu_solve, t.klu_free = (
            ffi.cast('void *', _cbuild.address(f))
            for f in (lib.klu_factor, lib.klu_refactor, lib.klu_solve, lib.klu_free_numeric))
        t.klu_common = ffi.cast('void *', sk._ct.addressof(sk._common))
        self.sym = sk._symbolic
        t.klu_sym = ffi.cast('void *', self.sym)
        t.resid_tol = float(sk.REFACTOR_RESIDUAL_TOL)
        t.M, t.nnz = M, sk.nnz


def block_step(walk, Ax, base, s):
    """`walk.factor_values(Ax)` and the stage system's solve of `base`
    (``C_n P``, `m x k`, C-ordered float64) once a stage, for `s` stages,
    in one C call: the solution (`s m x k`, Fortran order: its last stage
    block the new `P`) with the walk's numeric, generation and counts as
    `factor_values` leaves them -- or None (counted): handed back, the
    Python step runs.  Raises the Python's LinAlgError where it raises one
    (a fresh factor that fails, a solve that fails)."""
    drv = _driver or driver()
    if not drv:
        return _no('pss.sensklu:driver')
    sk = walk.solver
    if not sk._symbolic:
        return _no('pss.sensklu:unanalysed')
    ffi, cfn = drv
    ctx = sk.__dict__.get('_sensklu_c')
    if ctx is None or ctx.sym is not sk._symbolic:
        ctx = sk._sensklu_c = _Ctx(ffi, sk)
    m, k = base.shape
    M = s * m
    if M != sk.m or Ax.shape != (sk.nnz,):
        return _no('pss.sensklu:size')
    ## (column-major `M x k`: the transpose of a C-ordered `k x M`, the
    ## strides of `_run`'s Fortran copy -- cffi takes a C-contiguous buffer)
    Xt = np.empty((k, M))
    X = Xt.T
    t = ctx.t
    keep = (ffi.from_buffer('double *', Ax), ffi.from_buffer('double *', base),
            ffi.from_buffer('double *', Xt))
    t.Ax, t.base, t.X = keep
    t.num = ffi.cast('void *', walk._num or 0)
    t.m, t.k, t.s = m, k, s
    status = cfn(t)
    t.Ax = t.base = t.X = ffi.NULL
    if status == _CODES['SENSKLU_BAIL_FLAGS']:
        return _no('pss.sensklu:bail:' + BAIL.get(status, str(status)))
    ## (`factor_values`'s course, made: its generation, numeric and counts)
    walk.gen += 1
    walk._num = int(ffi.cast('uintptr_t', t.num)) or None
    walk.refactors += t.refactors
    walk.fallbacks += t.fallbacks
    walk.factors += t.factors
    sk.numerics += t.factors
    if status == _CODES['SENSKLU_SINGULAR']:
        raise np.linalg.LinAlgError('Singular matrix')
    if status == _CODES['SENSKLU_SOLVE_FAILED']:
        raise np.linalg.LinAlgError('klu solve failed')
    return X


def step(walks, xn, h, tab, P, walk_lu):
    """The plain walk's step map on the walk's KLU from `xn` through the
    step just taken: `(the new P, the step on a full-width P)` -- the second
    for `_LiveColumns.advance`'s non-finite case, on this step's
    factorisation -- or None (counted): the Python step runs."""
    if not SENSKLU_C:
        return _no('pss.sensklu:off')
    g = _MOD.get('g') or _pieces()
    if g is None or _now() != g:
        return _no('pss.sensklu:patched')
    W, _IN, RC, LS, ST, _NU = _MOD['mods']
    d = walks.__dict__
    if any(name in d for name in _WALK_NAMES) or type(walk_lu) is not LS._KLUWalk:
        return _no('pss.sensklu:patched')
    A = tab[0]
    s = A.shape[0]
    m = walks.cir.n - 1
    if not (type(P) is np.ndarray and P.dtype == _F64 and P.ndim == 2 and P.shape[0] == m
            and P.shape[1] >= 1):
        return _no('pss.sensklu:p')
    sk = walk_lu.solver
    if not W._own_factorise(walks) or sk is not walks._step_klu(s * m):
        return _no('pss.sensklu:klu')
    tr = walks._transient()
    got = RC.final_block(tr, tr.last_step.Y, h, A)
    if got is None or got[1] is not sk:
        return _no('pss.sensklu:final')
    Ax = got[0]
    if not (type(Ax) is np.ndarray and Ax.dtype == _F64 and Ax.flags.c_contiguous):
        return _no('pss.sensklu:final')
    Cn = np.asarray(walks._C_at(xn))
    base = np.ascontiguousarray(Cn @ P)
    X = block_step(walk_lu, Ax, base, s)
    if X is None:
        return None    # paths: counted by block_step
    _PC['pss.sensklu:served'] += 1
    _A, b, c = tab

    def full(Pf):
        ## (this step's factorisation, as `st.solve` holds it)
        lu = LS._Factored(LS._KLUBorrowed(walk_lu, walk_lu._num), lambda f, v: f.solve(v),
                          lambda f, v: f.solve_transposed(v), sparse=True)
        return ST._StageStep(Cn, None, h, A, b, c, lu=lu).solve(Pf)
    return X[(s - 1) * m:s * m], full
