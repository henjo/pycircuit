"""THE BRANCH CHECK'S CHEAP SCREEN IN C (the C stepping loop plan's stage 3,
step d; 2026-10-11).  After every converged step the branch check screens
`C` at the converged point (`_BranchCheck._branch_screen_C`): |C|'s largest
entry, the run's running maximum of it, the diagonal's largest and its
smallest positive entry -- six numpy reductions over n x n, ~6 us a step at
15 unknowns and ~60 us at 269 -- and only where that fires the SVD decides.
Here one C call makes the cheap level (`csrc/branch.c`), bit for bit the
same comparisons on the same values; the SVD stays the Python's
(`_BranchCheck._branch_svd`).

THE PATTERN.  Where `C` is the one the evaluation core wrote on its
pattern -- the C Newton's converged point or `_tran_core.evaluate`'s, which
record it with the pattern (`_C_pattern`) -- only the pattern's entries can
be nonzero, and only they are read (the largest of the others is +0.0).
Any other `C` is read whole, up to `_cbuild.INPLACE_MIN_ENTRIES` entries;
a larger one is numpy's (vectorised reductions beat a plain C loop there:
the error test's lesson, roadmap 113).

EXACT OR NOT AT ALL: an entry that is not finite hands the screen back
(numpy's reductions and the SVD then say what they said); the running
maximum is the caller's to assign.  DECLINES (`_paths.no`): the switch, no
build, a patched or shadowed `_branch_screen_C`, a `C` that is not a
C-contiguous float64 square array, a large `C` without the pattern.
`ENABLED` (env `PYCIRCUIT_BRANCH_C=0`) is read on every call.
"""
import os

import numpy as np

from pycircuit.circuit import _cbuild, _paths

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_BRANCH_C', '1') != '0'
STATUS = 'not loaded'

BRANCH_C = _cbuild.csrc('branch.c')
BRANCH_CDEF = _cbuild.csrc('branch.h')
_CODES = _cbuild.constants(BRANCH_CDEF)
BAIL = _cbuild.bails(_CODES, 'BRANCH_BAIL_')
QUIET, COLLAPSED, SVD = (_CODES['BRANCH_' + k] for k in ('QUIET', 'COLLAPSED', 'SVD'))

_F64 = np.dtype(np.float64)
_driver = None
_OWN = None


def driver():
    """`(ffi, cfn, the struct, the pattern it points at)` once loaded, or
    None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(BRANCH_C, BRANCH_CDEF)
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = [ffi, cfn, ffi.new('branch_t *'), None, None]
        STATUS = 'c'
    return _driver or None


def _own():
    """`_BranchCheck._branch_screen_C` as defined (read once it is its
    module's own), or None."""
    global _OWN
    if _OWN is None:
        from pycircuit.circuit._tran_branch import _BranchCheck
        f = _BranchCheck._branch_screen_C
        if not _paths.genuine(f, '_BranchCheck._branch_screen_C',
                              'pycircuit.circuit._tran_branch'):
            return None
        _OWN = f
    return _OWN


def screen(tr, C):
    """`(verdict, nrm, ref)` -- `_branch_screen_C`'s cheap level for `tr` on
    `C`: QUIET, COLLAPSED or SVD, |C|'s largest entry, the new running
    maximum -- or None: declined or handed back, nothing written."""
    if not ENABLED:
        return _no('branch_c:off')
    drv = _driver or driver()
    if not drv:
        return _no('branch_c:driver')
    own = _OWN or _own()
    if own is None or type(tr)._branch_screen_C is not own or '_branch_screen_C' in tr.__dict__:
        return _no('branch_c:patched')
    if not (type(C) is np.ndarray and C.dtype is _F64 and C.ndim == 2
            and C.shape[0] == C.shape[1] and C.shape[0] > 0 and C.flags.c_contiguous):
        return _no('branch_c:C')
    n = C.shape[0]
    td = tr.__dict__
    mark = td.get('_C_pattern')
    pat = mark[1] if mark is not None and mark[0] is C else None
    if pat is None and n * n > _cbuild.INPLACE_MIN_ENTRIES:
        return _no('branch_c:large')
    cmax = td.get('_branch_cmax', 0.0)
    tol = tr.BRANCH_SCREEN_TOL
    if type(cmax) is not float or type(tol) is not float:
        return _no('branch_c:x')
    ffi, cfn, s = drv[0], drv[1], drv[2]
    if pat is None:
        s.npat = 0
    else:
        if drv[3] is not pat:
            ## (the pattern pointed at once, while it stands)
            drv[4] = ffi.from_buffer('int64_t *', pat)
            drv[3] = pat
            s.pat = drv[4]
        s.npat = len(pat)
    s.n = n
    s.C = ffi.from_buffer('double *', C)
    s.cmax, s.tol_rel = cmax, tol
    code = cfn(s)
    if code > SVD:
        return _no('branch_c:bail:' + BAIL.get(code, str(code)))
    _PC['branch_c:served'] += 1
    if pat is not None:
        _PC['branch_c:pattern'] += 1
    return code, s.nrm, s.ref
