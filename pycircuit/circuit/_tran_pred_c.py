"""THE STAGE PREDICTOR'S FAST PATH IN C: `_StagePredictor._predict_fast` --
the seed every multistep step starts its Newton from, extrapolated from the
accepted history -- as `csrc/pred.c`, which the stepping loop in C calls
(`csrc/loop.c`), the shooting's stage walk (`shooting/_stagewalk_c`) and,
for a coupled stage step's three seeds, `_StagePredictor._pred_stages`
(`stages`).

OWNS on the transient `_pred_c` (`stages`' struct and history buffers);
the C loop keeps its `pred_t` (made by `new_state` in `_tran_loop_c._bind`)
and the history in its own buffers; `predict` reads `tr._pred_hist` and
writes nothing.  `_paths` prefix `pred_c:`, counted by `predict` and
`stages`: `pred_c:served` (`predict`'s), the declines (`pred_c:off`,
`driver`, `hist`, `deg`, `periodic`), the general path (`pred_c:general`)
and the hand-backs (`pred_c:bail:<why>`); a seed `stages` makes counts as
the Python's `pred:fast` (its caller's), and inside the C loop a prediction
counts so too, re-created by `_tran_loop_c._write_back`.  Switch: `ENABLED`
(env `PYCIRCUIT_PRED_C=0`), read by `predict` and `stages` on every call;
the C loop does not read it.

EXACT OR NOT AT ALL.  Bit for bit the Python's: the nodes newest first with
the 1e-13 dedupe and the tie test, the fit's times in Python's float order,
`np.vander`'s running products, numpy's own `dgesv` (the weights) and
`dgemv` (`w @ X`: `cblas_dgemv(RowMajor, Trans, ...)` is that product), the
clamp in numpy's `clip` order.  A floating-point flag numpy would report,
and every case the fast path sends to the general one, hand the prediction
back -- an underflow only where numpy's error state reports it, for
`stages` and the stage walk (`pred_t.under_ok`); the C loop and `predict`
hand it back always.

`predict(tr, ttarget, deg)` packs `tr._pred_hist` and answers as
`_predict_fast` does -- the prediction, None, or `_tran_predictor._GENERAL`
-- for the twin's tests.  History: `doc/transient_history.md`.
"""
import os

import numpy as np

from pycircuit.circuit import _cbuild, _paths, func

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_PRED_C', '1') != '0'
STATUS = 'not loaded'

PRED_C = _cbuild.csrc('pred.c')
PRED_CDEF = _cbuild.csrc('pred.h')
_CODES = _cbuild.constants(PRED_CDEF)
BAIL = _cbuild.bails(_CODES, 'PRED_BAIL_')
DONE, NONE, GENERAL = (_CODES['PRED_' + k] for k in ('DONE', 'NONE', 'GENERAL'))
#: the capacities (pred.h): the nodes a history holds, the fit's nodes -- so
#: the largest degree the struct holds
_CAP = _cbuild.constants(PRED_CDEF)
HIST_MAX, FIT_MAX = _CAP['PRED_HIST_MAX'], _CAP['PRED_FIT_MAX']
MAX_DEG = FIT_MAX - 1

_driver = None


def driver():
    """`(ffi, cfn, numpy's dgesv and dgemv addresses)` once loaded, or None
    (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        from pycircuit.circuit import _tran_newton_c
        nd = _tran_newton_c.driver()
        if nd is None:
            STATUS = f'off (numpy\'s LAPACK: {_tran_newton_c.STATUS})'
            return None
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(PRED_C, PRED_CDEF, 'pcx_pred')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, nd[2]['dgesv'], nd[2]['dgemv'])
        STATUS = 'c'
    return _driver or None


def new_state(ffi, dgesv, dgemv, n):
    """A predictor struct for `n` unknowns and its buffers (kept with it)."""
    s = ffi.new('pred_t *')
    buf = {'out': np.empty(n), 'xt': np.empty(FIT_MAX * n), 'motion': np.empty(n)}
    keep = [s, buf]
    for k, arr in buf.items():
        h = ffi.from_buffer('double *', arr)
        keep.append(h)
        setattr(s, k, h)
    s.n = n
    s.dgesv = ffi.cast('void *', dgesv)
    s.dgemv = ffi.cast('void *', dgemv)
    s.wk = 0
    s.under_ok = 0
    return s, buf, keep


def predict(tr, ttarget, deg=None):
    """`tr._predict_fast(ttarget, deg)` in C: the prediction (a fresh
    array), None, or `_tran_predictor._GENERAL` -- also for a flag or a
    history the struct does not hold."""
    from pycircuit.circuit._tran_predictor import _GENERAL
    if not ENABLED:
        _no('pred_c:off')
        return _GENERAL
    drv = _driver or driver()
    if not drv:
        _no('pred_c:driver')
        return _GENERAL
    hist = getattr(tr, '_pred_hist', None) or ()
    if type(hist) is not list and len(hist):
        _no('pred_c:hist')
        return _GENERAL
    if deg is None:
        deg = tr._pred_degree()
    if not (1 <= int(deg) <= MAX_DEG) or len(hist) > HIST_MAX:
        _no('pred_c:deg')
        return _GENERAL
    if getattr(tr, '_periodic_rows', None):
        _no('pred_c:periodic')
        return _GENERAL
    hk = len(hist)
    if hk < 2:
        return None
    n = len(hist[-1][1])
    ht = np.array([float(e[0]) for e in hist])
    hx = np.empty((hk, n))
    for i, e in enumerate(hist):
        x = e[1]
        if not (type(x) is np.ndarray and x.dtype == np.float64 and x.shape == (n,)):
            _no('pred_c:hist')
            return _GENERAL
        hx[i] = x
    ffi, cfn, dgesv, dgemv = drv
    s, buf, _keep = new_state(ffi, dgesv, dgemv, n)
    s.hk, s.deg, s.ttarget, s.clamp = hk, int(deg), float(ttarget), float(tr.PRED_CLAMP)
    s.ht = ffi.from_buffer('double *', ht)
    s.hx = ffi.from_buffer('double *', hx)
    status = cfn(s)
    if status == DONE:
        _PC['pred_c:served'] += 1
        return buf['out'].copy()
    if status == NONE:
        return None
    _no('pred_c:general' if status == GENERAL else 'pred_c:bail:' + BAIL.get(status, str(status)))
    return _GENERAL


def stages(tr, hist, x, ttargets, deg):
    """`_pred_stages`' seeds for `ttargets` in C where its fast path serves
    them: pcx_pred a target (`_predict_fast`'s twin) on a state kept on the
    transient (`_pred_c`), the history's nodes copied in once a call --
    each seed a fresh array, a target with no prediction the point `x` --
    or None (counted): the Python's loop makes them."""
    if not ENABLED:
        return _no('pred_c:off')
    drv = _driver or driver()
    if not drv:
        return _no('pred_c:driver')
    hk = len(hist)
    if not (2 <= hk <= HIST_MAX and 1 <= int(deg) <= MAX_DEG):
        return _no('pred_c:deg')
    if getattr(tr, '_periodic_rows', None):
        return _no('pred_c:periodic')
    n = len(hist[-1][1])
    td = tr.__dict__
    st = td.get('_pred_c')
    if st is None or st[0] != n:
        ffi, _cfn, dgesv, dgemv = drv
        s, buf, keep = new_state(ffi, dgesv, dgemv, n)
        ht, hx = np.empty(HIST_MAX), np.empty((HIST_MAX, n))
        keep += [ht, hx]
        s.ht = ffi.from_buffer('double *', ht)
        s.hx = ffi.from_buffer('double *', hx)
        st = td['_pred_c'] = (n, s, buf, keep, ht, hx)
    _n, s, buf, _keep, ht, hx = st
    for i, e in enumerate(hist):
        xi = e[1]
        if not (type(xi) is np.ndarray and xi.dtype == np.float64 and xi.shape == (n,)):
            return _no('pred_c:hist')
        ht[i] = e[0]
        hx[i] = xi
    s.hk, s.deg, s.clamp = hk, int(deg), float(tr.PRED_CLAMP)
    ## (an underflow numpy ignores is no hand-back)
    s.under_ok = 1 if func._underflow_ignored() else 0
    cfn = drv[1]
    out = []
    for tt in ttargets:
        s.ttarget = float(tt)
        status = cfn(s)
        if status == DONE:
            out.append(buf['out'].copy())
        elif status == NONE:
            out.append(np.array(x, dtype=float))
        else:
            return _no('pred_c:general' if status == GENERAL
                       else 'pred_c:bail:' + BAIL.get(status, str(status)))
    return out
