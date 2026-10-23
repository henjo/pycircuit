"""THE BRANCH CHECK'S CHEAP SCREEN IN C: the cheap level of
`_BranchCheck._branch_screen_C` -- after every converged step, `C` at the
converged point: |C|'s largest entry, the run's running maximum of it, the
diagonal's largest and its smallest positive entry -- and, where that
fires, THE CERTIFICATE (`_BranchCheck._branch_certify`: rank C >= r0
proved without the SVD), as one call of `csrc/branch.c`, bit for bit
the same comparisons on the same values.  Only where neither decides does
the SVD, and the SVD stays the Python's (`_BranchCheck._branch_svd`).

OWNS nothing on others' objects: `screen` reads the transient's
`_C_pattern`, `_branch_cmax`, `BRANCH_SCREEN_TOL`, `BRANCH_CERT_GROWTH` and
the run's certificate data (`_branch_cert_kept`) and returns the new
running maximum for its caller to assign (`_tran_branch._branch_screen`);
at the run's first screen the cheap level fires on it has the data made
(`_branch_cert_data`, as the Python makes them there) and calls again.  Its
struct, the pattern's pointer and the certificate's binding are the
module's (`driver`).  `_paths` prefix `branch_c:`: `branch_c:served`,
`branch_c:pattern` (served on the pattern), `branch_c:cert` (the screens
the certificate proved), the declines below and the hand-backs
(`branch_c:bail:<why>`); the C loop re-creates `served`, `pattern` and
`cert` for the screens it makes.  Switch: `ENABLED` (env
`PYCIRCUIT_BRANCH_C=0`), read by `screen` on every call and by the C loop on
every entry (`_tran_loop_c._screen`).

THE PATTERN.  Where `C` is the one the evaluation core wrote on its
pattern -- the C Newton's converged point or `_tran_core.evaluate`'s, which
record it with the pattern (`_C_pattern`) -- only the pattern's entries can
be nonzero, and only they are read (the largest of the others is +0.0).
Any other `C` is read whole, up to `_cbuild.INPLACE_MIN_ENTRIES` entries;
a larger one declines, and numpy's reductions screen it.

THE CERTIFICATE'S BLOCK (`cert_bind`): the run's data (`_tran_branch._Cert`:
B's rows and columns, its blocks, the work for its largest) and numpy's own
`scipy_dgesv_64_` -- numpy's inv is that dgesv on the identity -- bound on
the module's struct when the data or the growth change, and on the C
loop's own once a run (`_tran_loop_c._bind_screen`).

EXACT OR NOT AT ALL: an entry that is not finite hands the screen back
(numpy's reductions and the SVD then say what they said).  DECLINES
(`_paths.no`): the switch, no build, a patched or shadowed
`_branch_screen_C` or piece of the certificate, a `C` that is not a
C-contiguous float64 square array, a large `C` without the pattern, a
running maximum, tolerance or growth that is not a float, certificate data
for another size, or none of numpy's dgesv to prove them with.  History:
`doc/transient_history.md`.
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
    """`[ffi, cfn, the struct, the pattern it points at, its pointer, the
    certificate bound, its pointers, numpy's dgesv or None]` once loaded, or
    None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(BRANCH_C, BRANCH_CDEF, 'pcx_branch')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        ## (numpy's own dgesv, the certificate's: None where it is not
        ## found -- a run with certificate data then declines)
        try:
            npl = _cbuild.openblas(np)
            dgesv = None if npl is None else _cbuild.address(npl.scipy_dgesv_64_)
        except (OSError, AttributeError):
            dgesv = None
        ## [ffi, cfn, the struct, its pattern, the pattern's pointer, the
        ## certificate bound (data, growth), its pointers, dgesv]
        _driver = [ffi, cfn, ffi.new('branch_t *'), None, None, None, (), dgesv]
        STATUS = 'c'
    return _driver or None


#: the pieces `screen` stands in for, read from the transient's class
#: (`_own`); none may be shadowed on the instance
_METHODS = ('_branch_screen_C', '_branch_certify', '_branch_cert_data', '_branch_cert_kept')


def _own():
    """`(the methods _METHODS names, the module's _certify and _cert_data)`
    as defined (read once each is its module's own), or None."""
    global _OWN
    if _OWN is None:
        from pycircuit.circuit import _tran_branch as tb
        mod = 'pycircuit.circuit._tran_branch'
        meths = tuple(getattr(tb._BranchCheck, m) for m in _METHODS)
        if not (all(_paths.genuine(f, '_BranchCheck.' + m, mod) for f, m in zip(meths, _METHODS))
                and _paths.genuine(tb._certify, '_certify', mod)
                and _paths.genuine(tb._cert_data, '_cert_data', mod)):
            return None
        _OWN = (meths, tb._certify, tb._cert_data)
    return _OWN


def _patched(tr, own):
    """A piece of the screen or the certificate not as defined: on the
    transient's class, its instance, or the module."""
    from pycircuit.circuit import _tran_branch as tb
    T, td = type(tr), tr.__dict__
    return (tuple(getattr(T, m) for m in _METHODS) != own[0]
            or any(m in td for m in _METHODS)
            or tb._certify is not own[1] or tb._cert_data is not own[2])


def cert_bind(ffi, s, cert, growth, dgesv):
    """THE CERTIFICATE's block on the struct `s` (made by `ffi`): the run's
    data `cert` (`_tran_branch._Cert`), the growth, numpy's dgesv -- the
    pointers returned for the caller to keep; none (`r0` 0) where `cert` is
    None."""
    if cert is None:
        s.r0 = 0
        return ()
    fb = ffi.from_buffer
    keep = (fb('int64_t *', cert.ci), fb('int64_t *', cert.cj), fb('int64_t *', cert.bstart),
            fb('int64_t *', cert.bid), fb('double *', cert.work), fb('int64_t *', cert.ipiv))
    s.ci, s.cj, s.bstart, s.bid, s.work, s.ipiv = keep
    s.nb, s.growth = cert.nb, growth
    s.dgesv = ffi.cast('void *', dgesv)
    s.r0 = cert.r0
    return keep + (cert,)


def _cert(drv, tr, n, data):
    """The module's struct bound to the run's certificate data `data` (or
    none), re-bound only where the data or the growth moved: False where it
    cannot be bound (growth not a float, data for another size, no
    dgesv)."""
    growth = tr.BRANCH_CERT_GROWTH
    if type(growth) is not float:
        return _no('branch_c:x')
    if data is not None and (data.n != n or drv[7] is None):
        return _no('branch_c:cert')
    key = drv[5]
    if key is None or key[0] is not data or key[1] != growth:
        drv[6] = cert_bind(drv[0], drv[2], data, growth, drv[7])
        drv[5] = (data, growth)
    return True


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
    if own is None or _patched(tr, own):
        return _no('branch_c:patched')
    if not (type(C) is np.ndarray and C.dtype is _F64 and C.ndim == 2
            and C.shape[0] == C.shape[1] and C.shape[0] > 0 and C.flags.c_contiguous):
        return _no('branch_c:C')
    n = C.shape[0]
    td = tr.__dict__
    ## (the core's pattern where `C` is the C the core wrote on it -- or one
    ## of a stage step's C's, `_tran_radau_c.solve`'s mark)
    mark = td.get('_C_pattern')
    pat = None
    if mark is not None:
        m0 = mark[0]
        if m0 is C or (m0.__class__ is tuple and any(o is C for o in m0)):
            pat = mark[1]
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
    known, data = tr._branch_cert_kept()
    if not _cert(drv, tr, n, data):
        return None    # paths: counted by _cert
    s.ncert = 0
    code = cfn(s)
    if code == SVD and not known:
        ## (the run's first screen the cheap level fires on: the
        ## certificate's data made, as `_branch_certify` makes them, and the
        ## screen made again with them)
        data = tr._branch_cert_data()
        if data is not None:
            if not _cert(drv, tr, n, data):
                return None    # paths: counted by _cert
            code = cfn(s)
    if code > SVD:
        return _no('branch_c:bail:' + BAIL.get(code, str(code)))
    _PC['branch_c:served'] += 1
    if pat is not None:
        _PC['branch_c:pattern'] += 1
    if s.ncert:
        _PC['branch_c:cert'] += 1
    return code, s.nrm, s.ref
