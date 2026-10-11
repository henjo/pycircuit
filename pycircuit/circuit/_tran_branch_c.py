"""THE BRANCH CHECK'S CHEAP SCREEN IN C: the cheap level of
`_BranchCheck._branch_screen_C` -- after every converged step, `C` at the
converged point: |C|'s largest entry, the run's running maximum of it, the
diagonal's largest and its smallest positive entry -- as one call of
`csrc/branch.c`, bit for bit the same comparisons on the same values.  Only
where that fires does the SVD decide, and the SVD stays the Python's
(`_BranchCheck._branch_svd`).

OWNS nothing on others' objects: `screen` reads the transient's
`_C_pattern`, `_branch_cmax` and `BRANCH_SCREEN_TOL` and returns the new
running maximum for its caller to assign (`_tran_branch._branch_screen`);
its struct and the pattern's pointer are the module's (`driver`).  `_paths`
prefix `branch_c:`: `branch_c:served`, `branch_c:pattern` (served on the
pattern), the declines below and the hand-backs (`branch_c:bail:<why>`);
the C loop re-creates `served` and `pattern` for the screens it makes.
Switch: `ENABLED` (env `PYCIRCUIT_BRANCH_C=0`), read by `screen` on every
call and by the C loop on every entry (`_tran_loop_c._screen`).

THE PATTERN.  Where `C` is the one the evaluation core wrote on its
pattern -- the C Newton's converged point or `_tran_core.evaluate`'s, which
record it with the pattern (`_C_pattern`) -- only the pattern's entries can
be nonzero, and only they are read (the largest of the others is +0.0).
Any other `C` is read whole, up to `_cbuild.INPLACE_MIN_ENTRIES` entries;
a larger one declines, and numpy's reductions screen it.

EXACT OR NOT AT ALL: an entry that is not finite hands the screen back
(numpy's reductions and the SVD then say what they said).  DECLINES
(`_paths.no`): the switch, no build, a patched or shadowed
`_branch_screen_C`, a `C` that is not a C-contiguous float64 square array,
a large `C` without the pattern, a running maximum or tolerance that is not
a float.  History: `doc/transient_history.md`.
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
    """`[ffi, cfn, the struct, the pattern it points at, its pointer]` once
    loaded, or None (`STATUS`)."""
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
