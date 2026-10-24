"""RADAU'S TRANSFORM NEWTON IN C: `_RadauStages._transform_loop`, the
simplified Newton of the cost transform (`_RadauStages._rk_step_transformed`:
PSS's default on PSP-class circuits under 'auto'), in ONE C call
(`solve`); and `_RadauStages._radau_frozen`'s two factors in another
(`fold`).  Bit for bit the Python (pinned in `tests/test_pinned_pairs.py`).

OWNS.  On the transient, `tr._radau_tc` (`_Ctx`: the `radau_tf_t` struct,
its buffers, what was last set in it -- the tableau, `P` and `V`, the kept
LU's addresses, the pattern's arrays, the KLU handles -- and the limiter
walk's fields, which `_tran_newton_c._walk_ready` sets) and `tr._tc_end`
(the step's end passes, handed to `_RadauStages._stage_end_passes`, which
pops it).  On a call it writes what the Python loop leaves: the step's
source memo (`tr._u_memo`) and its count `umemo:hit`, and the frozen
factors' state -- the kept LU's `NumpyLU._info` (or, the real factor on
the analysis solver's StructuralKLU, its `_numeric` and counts), the
complex solver's `ComplexKLUSolver.refactors` and `_fresh`.  The module keeps `fold`'s
driver, buffers and last pattern (`_FOLD`).  `_paths`: `radau_tc:` --
`radau_tc:<why>` a decline, `radau_tc:bail:<why>` a hand-back (`BAIL`),
`radau_tc:served` (the loop's outcome made), `radau_tc:noconv` (of them,
maxiter reached), `radau_tc:probe` (the early abort's first iteration alone,
not converged: `solve(..., probe=True)`), `radau_tc:iterations`,
`radau_tc:end` -- and `fold`'s
`radau.frozen:fold` (made), `radau.frozen:fold.same` (its packed values
too) and `radau.frozen:fold.<why>` (not made).  Switches, each read on
every call: `PYCIRCUIT_RADAU_TC=0` (`ENABLED`: the loop),
`PYCIRCUIT_RADAU_TC_END=0` (`END`: the end passes), `PYCIRCUIT_RADAU_FOLD=0`
(`FOLD`: the factors); `STATUS` says why the loop's path is off where it
is, `_FOLD['status']` the factors'.

THE LOOP.  `C` and `G` are frozen at ``x_n`` and the coupled solve goes
through ``A^{-1}``'s eigenbasis, its factors made once a step
(`_radau_frozen`).  Per iteration the Python loop makes the three stages'
passes through `_tran_core.passes`, the residual (`_coupled_stage_system`),
the frozen solve (`_FrozenTransform.solve`: the right-hand sides, numpy's
kept LU, the complex factor's KLU solve and its residual check, the
update), and per stage the step, the limiter and the convergence test;
False after `maxiter` iterations, and the step raises (the dense coupled
Newton takes it).
The C makes, per iteration, the stages' ``q`` and ``i`` through the
evaluate core's own C (formula 4, as `passes`), ``K_j = -(i_j + u_j)``,
``F_i = (q_i - q_n) - h (((0 + A_i0 K_0) + A_i1 K_1) + A_i2 K_2)`` reduced
(`_tran_radau_c`'s assembly, without the blocks); the right-hand sides
``-((P_k0 F0 + P_k1 F1) + P_k2 F2)`` -- a complex scalar times a real
vector, as numpy multiplies it: the vector cast to complex, numpy's
complex product; numpy's kept LU (`linearsolver.NumpyLU`'s protocol: its
OpenBLAS's `dgesv` at the first solve, `dgetrs` on the kept factors
after) -- or, where `_radau_frozen` keeps none (the analysis solver an
AutoSolver whose choice is the circuit's StructuralKLU, `_real_klu`),
`StructuralKLU.solve`'s sequence (pcx_kluseq.h: the refactor on the
numeric it holds, a fresh factor at a zero pivot or a failed residual --
the sparse verdict, numpy's dgemv where it does not decide) on the
factor's values at the pattern, taken once a call with `_ours`'s count (a
nonzero off the pattern declines: `klu_foreign`); the complex factor's KLU on the solver's own handles (the refactor
`solve_prepared` would make, at most once a record, then `klu_z_solve`);
``dY_i = V_i0 w0 + 2 Re(V_i1 w1)``; per stage ``Y_prev + insert(dY_i,
iref, 0.0)``, the limiter walk's own driver and the limited step's largest
reduced entry; the convergence test after the update, as the Python tests.
Its two outcomes are served: converged, and `maxiter` iterations without
(the Python loop's own iterates, so its False).  THE EARLY ABORT
(`Transient.radau_early_abort`) calls it for one iteration (`probe`),
reads the iterate in the buffer (`iterate`) and, not abandoned, goes on
from it (`resume`: no readiness checks, no source call -- nothing the C
reads moved -- the iterates one call makes; handed back, the Python loop
runs from the iterate, `seed`).

⚠ NUMPY'S COMPLEX PRODUCT IS READ, NOT ASSUMED.  With numpy's X86_V3
loops (AVX2 and FMA3) it FUSES -- ``re = fma(ar, br, -(ai*bi))``, ``im =
fma(ar, bi, ai*br)``, at every position of every length -- and a build or
CPU without FMA multiplies separately.  `cmul_mode` asks numpy itself once,
on inputs where the two forms round apart (and on signed zeros and
underflow, where only the fused form's sign rules hold), and the C takes
the form numpy showed; neither: the path is off.

⚠ THE RESIDUAL CHECK IS DECIDED WITH A MARGIN.  `solve_prepared` falls
back to a fresh factor where ``max|A x - b| / max|b|`` exceeds its
tolerance (1e-8), read with numpy's complex `abs` (its own algorithm, not
`hypot`) of SciPy's product.  The C computes the same quantities in its own
rounding and decides only outside the band their difference can reach
(``32 (m + 4) eps`` of ``max_k sum_j |A_kj| |x_j|`` over ``max|b|``, plus
``16 eps`` of the residual): a residual nearer the tolerance than that, or
past it, hands the iteration back and the Python decides, as it does.

TRANSACTIONAL.  The C writes nothing outside its buffers but the two
factorisations' state, which the Python would have left the same: numpy's
kept LU after its first `dgesv` (`NumpyLU._info`) -- or the real factor's
StructuralKLU, its sequence's numeric and counts where the outcome is
served, the one it found where handed back (any other freed:
`pcx_klu_finish`) -- and the KLU numeric holding a refactor of the step's
record (`ComplexKLUSolver._fresh`, `refactors`).  A refactor of the same values with the same pivots repeats
its bits and `dgetrs` on `dgesv`'s factors is `dgesv`'s solve, so a Python
loop run after a hand-back makes the C's iterates again.  The step's
source memo holds the three stage times (their calls made here, before
the C) and its hit count for every later iteration's.  Anything unusual --
a non-finite value, a floating-point exception in its own arithmetic
(numpy would warn), a non-finite real right-hand side (the analysis
solver's), a singular real factor (numpy's error), the KLU refactor or
solve failing, the residual check not decided or failing (the fallback's
fresh factor), the walk stopping -- hands the loop back: the source memo's
entries and counts rolled back, and the Python loop runs from the same
seed (the same answer, exception, warning).

THE STEP'S END PASSES (`END`).  Once converged the C makes, in the same
call, the passes `_stage_end_passes` makes next: all four at the last
stage, and `C` and `G` at the first two where the shooting reads them
(`_stage_G_read`) -- each `_tran_core.passes`'s core call with its bits,
outputs and dummies, in its order (the core's own arrays are its scratch,
each call writing what it reads: a call the Python then does not need
leaves nothing behind) -- and hands them over with the very stage list it
returns (`_tc_end`, the watch counter with it).  `_stage_end_passes` takes
them for those stages only, the counter unmoved, all four passes wanted at
the last stage and `passes`'s readiness holding as stamped
(`_tran_core.passes_held`), making the counters each of its calls makes
(`passes_count`); else it makes the calls (`radau.fuse:unhanded`).

THE FROZEN FACTORS (`FOLD`).  `_radau_frozen` makes its two factors --
``(lam0/h) Cr + Gr`` real, ``(lam1/h) Cr + Gr`` complex -- and the complex
one's packed values by `fold`, one small C call of its own: numpy's
operations in its order (two roundings for the real one; for the complex
one ``Cr`` cast to complex, numpy's product in `cmul_mode`'s form, then
``Gr`` cast to complex), the four flags its `errstate(all='raise')` raises
on, and the nonzeros where they fall in the last record's pattern
(`_csc_of_dense`'s order: column by column, rows ascending).  Every
decision stays in `_radau_frozen`, in its order: a flag declines as
numpy's error does, the kept LU is made as before, the record is
`prepare`'s own (`ComplexKLUSolver.prepare_values`) where the values fell
in the pattern and `prepare`'s otherwise.  A non-finite input is numpy's
(a NaN's payload follows its loops).  Its own call, not inside the
transform's: inside, a factor that would raise surfaces only after the
call's readiness checks have run, and the step would be made again around
them.  Its buffers are kept per size, their handles made once.

DECLINES, before anything is touched (`_paths.no`): the switch (and the
evaluate core's two), no frozen factors or no kept LU or marshalled record
(`_radau_frozen`'s own), no build, numpy's complex product of neither form,
numpy's error state not ignoring underflow (the C does not test it: numpy's
default ignores it), a toolkit that is not numeric, a subclass of the
transient or a patched piece of the loop, an instance shadow of the
circuit's passes or limiter or of the transient's pieces, a
`provided_function`, a stateful limiter or circuit-level limiting, a
bypass, no source memo, a circuit the core cannot serve, a limiting element
that is not its C kernel, a batch not ready, a KLU solver whose numeric is
not of the record's pattern (the Python analyses and factors), a singular
kept LU, inputs of another type or size.

History: `doc/transient_history.md`, `_tran_radau_tc.py`.
"""
import ctypes
import itertools
import math
import os

import numpy as np

from pycircuit.circuit import _cbuild, _paths, _watch

_PC = _paths.COUNTS
#: (float64's dtype: a check by identity first, `==` where it is not the
#: object)
_F64 = np.dtype(np.float64)
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_RADAU_TC', '1') != '0'
#: THE STEP'S END PASSES in the same call: the passes `_stage_end_passes`
#: makes next, made by the C after it converges and handed over
#: (`_tc_end`); env `PYCIRCUIT_RADAU_TC_END=0` turns it off
END = os.environ.get('PYCIRCUIT_RADAU_TC_END', '1') != '0'
#: THE FROZEN FACTORS IN C: `_radau_frozen`'s two factors and the complex
#: one's packed values made by one small C call (`fold`); env
#: `PYCIRCUIT_RADAU_FOLD=0` turns it off
FOLD = os.environ.get('PYCIRCUIT_RADAU_FOLD', '1') != '0'
STATUS = 'not loaded'

RADAU_TC = _cbuild.csrc('radau_tc.c')
RADAU_TC_CDEF = _cbuild.csrc('pcx_klu.h') + _cbuild.csrc('radau_tc.h')

FOLD_C = _cbuild.csrc('radau_fold.c')
FOLD_CDEF = _cbuild.csrc('radau_fold.h')

#: the C's results, named in radau_tc.h and radau_fold.h (`_CODES`,
#: `_FOLD_CODES`); `BAIL`: {code: why} of the solve's codes past 1
#: (converged) -- why it handed the loop back
_CODES = _cbuild.constants(RADAU_TC_CDEF)
_FOLD_CODES = _cbuild.constants(FOLD_CDEF)
BAIL = _cbuild.bails(_CODES, 'RADAU_TC_BAIL_')

_driver = None
_MOD = {}
#: `fold`'s loaded kernel and its status
_FOLD = {}


def _cmul_samples():
    """Pairs ``(a, b)`` of complex operands on which numpy's two possible
    complex products round apart -- drawn where the separate products
    differ from the fused ones -- plus signed zeros and underflow (the
    right-hand sides' ``P F``: a real vector cast to complex), as
    deterministic arrays."""
    rng = np.random.default_rng(20261005)
    ar, ai, br, bi = (rng.standard_normal(4000) * 10.0 ** rng.integers(-20, 20, size=4000)
                      for _ in range(4))
    keep = np.array([math.fma(a, c, -(d * e)) != a * c - d * e
                     for a, c, d, e in zip(ar, br, ai, bi)])
    a = (ar + 1j * ai)[keep][:64]
    b = (br + 1j * bi)[keep][:64]
    ## the real right-hand sides' corner: a tiny exact product rounding to
    ## a signed zero, against a zero imaginary product of either sign
    tiny = np.array([3e-200, -3e-200, 0.0, -0.0, 1.5, -2.5])
    za = np.array([complex(x, y) for x in (1e-200, -1e-200, 2.0, -0.0)
                   for y in (1.0, -1.0, 0.0, -0.0, 1e-200, -1e-200)])
    return a, b, za, tiny


def cmul_mode():
    """1 where numpy's complex product is the fused form, 0 where it is the
    separate one, None where it is neither -- read from numpy itself on
    `_cmul_samples`: a complex scalar times complex vectors (the update)
    and times real vectors (the right-hand sides), every position of
    lengths 1 to 17."""
    if 'cmul' in _MOD:
        return _MOD['cmul']
    a, b, za, tiny = _cmul_samples()

    def forms(x, y):
        return ((math.fma(x.real, y.real, -(x.imag * y.imag)),
                 math.fma(x.real, y.imag, x.imag * y.real)),
                (x.real * y.real - x.imag * y.imag, x.real * y.imag + x.imag * y.real))
    seen = set()
    for L in range(1, 18):
        for j in range(0, len(a), 4):
            x = np.complex128(a[j])
            ys = np.resize(b, L)
            got = x * ys
            for t in range(L):
                f, sp = forms(complex(x), complex(ys[t]))
                g = (float(got[t].real), float(got[t].imag))
                seen.add('f' if _same(g, f) else ('s' if _same(g, sp) else 'x'))
    for x in za:
        xc = np.complex128(x)
        got = xc * tiny
        for t in range(len(tiny)):
            f, sp = forms(complex(x), complex(float(tiny[t]), 0.0))
            g = (float(got[t].real), float(got[t].imag))
            ok_f, ok_s = _same(g, f), _same(g, sp)
            seen.add('f' if ok_f and not ok_s else ('s' if ok_s and not ok_f else
                                                       ('b' if ok_f else 'x')))
    seen.discard('b')
    _MOD['cmul'] = 1 if seen == {'f'} else (0 if seen == {'s'} else None)
    return _MOD['cmul']


def _same(g, f):
    """Two (re, im) pairs the same bits (the sign of a zero, a NaN's
    payload)."""
    return (np.float64(g[0]).tobytes() == np.float64(f[0]).tobytes()
            and np.float64(g[1]).tobytes() == np.float64(f[1]).tobytes())


def driver():
    """`(ffi, cfn, dgesv, dgetrs)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        from pycircuit.circuit import _tran_radau_c, linearsolver
        if _tran_radau_c.driver() is None:
            STATUS = f'off (the radau C: {_tran_radau_c.STATUS})'
            return None
        lap = linearsolver._numpy_lapack()
        if lap is None:
            STATUS = "off (numpy's own OpenBLAS not found)"
            return None
        mode = cmul_mode()
        if mode is None:
            STATUS = "off (numpy's complex product is neither form)"
            return None
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(RADAU_TC, RADAU_TC_CDEF, 'pcx_radau_tc')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _ct, gesv, getrs = lap
        from pycircuit.circuit import _tran_newton_c, _tran_radau_c, func
        _MOD['NC'], _MOD['RC'], _MOD['func'] = _tran_newton_c, _tran_radau_c, func
        _driver = (ffi, cfn, _cbuild.address(gesv), _cbuild.address(getrs))
        STATUS = 'c'
    return _driver or None


def _fold_driver():
    """`(ffi, cfn, struct)` for `fold` once loaded, or None (the reason in
    `_FOLD['status']`)."""
    d = _FOLD.get('drv')
    if d is None:
        d = _FOLD['drv'] = False
        if cmul_mode() is None:
            _FOLD['status'] = "off (numpy's complex product is neither form)"
            return None
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(FOLD_C, FOLD_CDEF, 'pcx_radau_fold')
        except (_cbuild.CompileError, OSError) as e:
            _FOLD['status'] = f'off ({e})'
            return None
        d = _FOLD['drv'] = (ffi, cfn, ffi.new('radau_fold_t *'))
        d[2].fused = _MOD['cmul']
        _FOLD['status'] = 'c'
    return d or None


def _fold_bufs(m):
    """`fold`'s buffers for size `m` (the inputs copied in, the outputs
    copied out), kept in `_FOLD['m']`: each buffer's handle made once and
    set in the struct."""
    ffi, _cfn, f = _FOLD['drv']
    bufs = (np.empty((m, m)), np.empty((m, m)), np.empty((m, m)),
            np.empty((m, m), dtype=np.complex128), np.empty(m * m, dtype=np.complex128))
    hs = [ffi.from_buffer('double *', a) for a in bufs]
    f.Cr, f.Gr, f.rf, f.cf, f.ax = hs
    f.m = m
    b = _FOLD['m'] = (m, *bufs, hs)
    return b


def _fold_nz(nz, m):
    """The last record's pattern as `fold`'s C reads it -- `(nz, at, handle)`,
    `at` 1 where it is laid out as `A.T`, 0 as `A` -- or None (another
    kind); its handle made once for each pattern."""
    if not (type(nz) is np.ndarray and nz.dtype == np.bool_ and nz.shape == (m, m)):
        return None
    if nz.flags.c_contiguous:
        at = 1
    elif nz.flags.f_contiguous:
        ## (`A.T != 0` is laid out as `A` where `A` is C-ordered: its
        ## transpose is the C-ordered view of that memory)
        at = 0
    else:
        return None
    ffi = _FOLD['drv'][0]
    c = _FOLD['nz'] = (nz, at, ffi.from_buffer('uint8_t *', nz if at else nz.T))
    return c


def fold(Cr, Gr, s0, c1, last):
    """`_RadauStages._radau_frozen`'s two factors from ONE C call, or None
    where the C does not make them -- the switch, no build, inputs other
    than two float64 square arrays of one size, a non-finite entry (each
    counted, `radau.frozen:fold.<why>`): the caller makes them in numpy, as
    before.  Else `(flagged, rf, cf, data)`: `flagged`, an operation raised
    a flag numpy raises on under `errstate(all='raise')` (the caller
    declines as it does; nothing else returned); `rf` and `cf`, ``s0 * Cr +
    Gr`` and ``c1 * Cr + Gr`` -- fresh arrays, numpy's bytes, for the
    scalars as the caller's expression makes them, finite; `data`, the
    complex factor's nonzeros where they fall in `last`'s pattern
    (`ComplexKLUSolver._csc_last`), in its order
    (`ComplexKLUSolver.prepare_values`) -- None where they fall elsewhere,
    there is no pattern or it has no nonzeros.  It stands in for the numpy
    calls, their errstate, the mask and the gather."""
    if not FOLD:
        return _no('radau.frozen:fold.off')
    d = _FOLD.get('drv')
    if d is None:
        d = _fold_driver()
    if not d:
        return _no('radau.frozen:fold.driver')
    m = Cr.shape[0] if type(Cr) is np.ndarray and Cr.ndim == 2 else 0
    if not (m and type(Gr) is np.ndarray and Cr.shape == (m, m) and Gr.shape == (m, m)
            and (Cr.dtype is _F64 or Cr.dtype == np.float64)
            and (Gr.dtype is _F64 or Gr.dtype == np.float64)):
        return _no('radau.frozen:fold.input')
    f = d[2]
    b = _FOLD.get('m')
    if b is None or b[0] != m:
        b = _fold_bufs(m)
    np.copyto(b[1], Cr)
    np.copyto(b[2], Gr)
    nnz = 0
    if last is not None:
        c = _FOLD.get('nz')
        if c is None or c[0] is not last[0]:
            c = _fold_nz(last[0], m)
        nnz = last[2].shape[0] if c is not None else 0
    if nnz:
        f.nz, f.nz_at, f.nnz = c[2], c[1], nnz
    else:
        f.nz = d[0].NULL
    f.s0, f.c1r, f.c1i = s0, c1.real, c1.imag
    d[1](f)
    if f.status == _FOLD_CODES['FOLD_STATUS_NONFINITE']:
        return _no('radau.frozen:fold.nonfinite')
    _PC['radau.frozen:fold'] += 1
    if f.status:
        return (True, None, None, None)
    rf, cf = b[3].copy(), b[4].copy()
    if nnz and f.same:
        _PC['radau.frozen:fold.same'] += 1
        return (False, rf, cf, b[5][:nnz].copy())
    return (False, rf, cf, None)


def _chain():
    """The loop's pieces the C stands in for, as defined (`_paths.genuine`),
    read once every piece is its module's own; None: not yet.  A call
    compares what the Python loop would call now (`_now`)."""
    if 'chain' not in _MOD:
        from pycircuit.circuit import _tran_radau, _tran_radau_c
        from pycircuit.circuit.linearsolver import ComplexKLUSolver, NumpyLU
        base = _tran_radau_c._MOD.get('chain') or _tran_radau_c._chain()
        if base is None:
            return None
        R, L = 'pycircuit.circuit._tran_radau', 'pycircuit.circuit.linearsolver'
        own = ((_tran_radau._RadauStages._transform_loop, '_RadauStages._transform_loop', R),
               (_tran_radau._FrozenTransform.solve, '_FrozenTransform.solve', R),
               (_tran_radau._transform_rhs, '_transform_rhs', R),
               (_tran_radau._transform_back, '_transform_back', R),
               (NumpyLU.solve, 'NumpyLU.solve', L),
               (ComplexKLUSolver.solve_prepared, 'ComplexKLUSolver.solve_prepared', L))
        if not all(_paths.genuine(*o) for o in own):
            return None
        _MOD['chain'] = base + tuple(o for o, _q, _m in own)
        _MOD['mods'] = (_tran_radau, _tran_radau_c, NumpyLU, ComplexKLUSolver)
    return _MOD['chain']


def _now(T, cir):
    """What the Python loop calls now, in `_chain`'s order."""
    tr_, rc, LU, ZS = _MOD['mods']
    return (rc._now(T, rc._MOD['mods'], cir)
            + (T._transform_loop, tr_._FrozenTransform.solve, tr_._transform_rhs,
               tr_._transform_back, LU.solve, ZS.solve_prepared))


class _Ctx:
    """One transient's struct and buffers for one core, size and reference
    row (the walk's fields named as `_tran_newton_c._Ctx`'s, whose
    `_walk_ready` sets them)."""

    __slots__ = (
        'A_bytes',
        'Ap_obj',
        'buf',
        'core',
        'ffi',
        'iref',
        'keep',
        'keep_ap',
        'klu_kk',
        'lu_addr',
        'n',
        'p_of',
        'rk_keep',
        's',
        'served',
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

    def __init__(self, core, n, iref, ffi, cfn_core, cffi_core, drv):
        m = n - 1
        self.core, self.n, self.iref, self.ffi = core, n, iref, ffi
        self.s = s = ffi.new('radau_tf_t *')
        self.buf = buf = {
            'Y': np.empty((3, n)), 'Yp': np.empty((3, n)), 'q': np.empty((3, n)),
            'iv': np.empty((3, n)), 'K': np.empty((3, n)), 'F': np.empty((3, n)),
            'R': np.empty(3 * m), 'b0': np.empty(m), 'x0': np.empty(m),
            'r1': np.empty(2 * m), 'x1': np.empty(2 * m), 'y1': np.empty(2 * m),
            'ab': np.empty(m), 'dY': np.empty(3 * m), 'stp': np.empty(n), 'tolv': np.empty(2),
            'qn': np.empty(n), 'u': np.empty((3, n)), 'seed': np.empty((3, n)),
            'dummy': np.zeros(max(n * n, 1)),
            ## (the step's end passes' outputs: `END`)
            'e_q': np.empty(n), 'e_i': np.empty(n), 'e_C': np.empty(n * n),
            'e_G': np.empty(n * n), 'e_C01': np.empty(2 * n * n), 'e_G01': np.empty(2 * n * n),
            'e_q01': np.empty(n),
            'rk_y': np.empty(max(m, 1)), 'rk_t': np.empty(max(2 * m, 2)),
        }
        self.keep = []
        fb = ffi.from_buffer
        for k, arr in buf.items():
            h = fb('double *', arr)
            self.keep.append(h)
            setattr(s, k, h)
        _f, _c, dgesv, dgetrs = drv
        s.dgesv = ffi.cast('void *', dgesv)
        s.dgetrs = ffi.cast('void *', dgetrs)
        s.core = ffi.cast('void *', int(cffi_core.cast('uintptr_t', core.cs)))
        s.core_fn = ffi.cast('void *', int(cffi_core.cast('uintptr_t', cfn_core)))
        a = core.arrays['bini']
        self.keep.append(a)
        s.bini = ffi.cast('double *', a.ctypes.data)
        g = core.arrays['binG']
        self.keep.append(g)
        s.binG = ffi.cast('double *', g.ctypes.data)
        s.n, s.iref = n, iref
        s.rk_on = 0
        s.rk_dgemv = ffi.cast('void *', _MOD['RC'].driver()[3])
        self.rk_keep = None
        s.walk_fn = ffi.NULL
        s.wn = 0
        ## (what does not change from step to step is set once: the product
        ## form, the passes' bits; the rest by `solve`, when it changes)
        s.fused = _MOD['cmul']
        s.bits = _MOD['bits']
        self.A_bytes = self.p_of = self.Ap_obj = self.keep_ap = self.klu_kk = None
        self.lu_addr = self.served = None
        self.walk = self.walk_ok = None
        self.walk_els = self.walk_dicts = self.walk_packs = self.walk_kern = ()
        self.walk_hw_els = self.walk_hwk = ()
        self.walk_stamp, self.walk_arms, self.walk_held = -1, 0, 0
        self.walk_watch = ()


def solve(tr, ctx, fz, seed, src, provided_function, lims, nobypass, reltol, abstol, maxit,
          probe=False):
    """The transform's Newton loop from the stages `seed` in C: the three
    converged stages as fresh full-width arrays, or False where `maxit`
    iterations did not converge (the Python loop's False), the source memo
    and the two factorisations' state as the Python loop leaves them -- or
    None: declined (`radau_tc:<why>`) or handed back
    (`radau_tc:bail:<why>`), nothing left behind that the Python loop would
    not make again, and the Python loop runs.  `fz` is the step's frozen factors (`_radau_frozen`'s
    `_FrozenTransform`); where `END` asks, the step's end passes are left
    in `tr._tc_end` with the very stage list returned.  `probe`: the early
    abort's first iteration (`_RadauStages._rk_step_transformed`) -- its
    False counted `radau_tc:probe`, not as the loop's outcome, and the
    iterate it reached read by `iterate`."""
    if not ENABLED:
        return _no('radau_tc:off')
    if fz is None or fz.prep is None:
        return _no('radau_tc:frozen')
    if provided_function is not None:
        return _no('radau_tc:pf')
    if lims:
        return _no('radau_tc:stateful')
    if not nobypass:
        return _no('radau_tc:bypass')
    drv = _driver or driver()
    if not drv:
        return _no('radau_tc:driver')
    if not _MOD['func']._underflow_ignored():
        return _no('radau_tc:underflow')
    ffi, cfn = drv[0], drv[1]
    chain = _MOD.get('chain') or _chain()
    if chain is None:
        return _no('radau_tc:patched')
    T = chain[0]
    if type(tr) is not T:
        return _no('radau_tc:class')
    cir = tr.cir
    if _now(T, cir) != chain[1:]:
        return _no('radau_tc:patched')
    NC = _MOD['NC']
    M = NC._MOD or NC._mods()
    if type(tr.toolkit) is not M['NumericToolkit']:
        return _no('radau_tc:toolkit')
    cd, td = cir.__dict__, tr.__dict__
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd or 'limit' in cd:
        return _no('radau_tc:shadow_cir')
    if ('_memo_put' in td or '_source_at' in td or '_coupled_stage_system' in td
            or '_stages_converged' in td or '_stage_limiter_states' in td
            or '_limiters_at_rest' in td or '_newton_limiter' in td
            or '_transform_loop' in td):
        return _no('radau_tc:shadow_tr')
    if M['limiting'].CIRCUIT_LEVEL:
        return _no('radau_tc:stateful')
    memo = td.get('_u_memo')
    if memo is None or not M['companion'].U_MEMO:
        return _no('radau_tc:umemo')
    tc = M['core']
    if not (tc.CORE and tc.CORE_PASSES):
        return _no('radau_tc:core')
    core = tc.core_for(tr)
    if core is None:
        return _no('radau_tc:unservable')
    n, iref = core.n, ctx.iref
    m = n - 1
    if not (0 <= iref < n) or n < 2 or ctx.m != m or fz.m != m:
        return _no('radau_tc:iref')
    h = ctx.h
    if type(h) is not float and type(h) is not np.float64:
        return _no('radau_tc:h')
    if (type(reltol) is not float and type(reltol) is not int
            and type(reltol) is not np.float64) or type(reltol) is bool:
        return _no('radau_tc:tol')
    qn, Amat = ctx.qn, ctx.Amat
    if not (type(qn) is np.ndarray and (qn.dtype is _F64 or qn.dtype == np.float64)
            and qn.shape == (n,) and type(Amat) is np.ndarray
            and (Amat.dtype is _F64 or Amat.dtype == np.float64) and Amat.shape == (3, 3)):
        return _no('radau_tc:ctx')
    if len(seed) != 3:
        return _no('radau_tc:seed')
    y0, y1, y2 = seed
    if not (type(y0) is np.ndarray and type(y1) is np.ndarray and type(y2) is np.ndarray
            and (y0.dtype is _F64 or y0.dtype == np.float64)
            and (y1.dtype is _F64 or y1.dtype == np.float64)
            and (y2.dtype is _F64 or y2.dtype == np.float64)):
        seed = [np.array(y, dtype=float) for y in seed]
        y0, y1, y2 = seed
    if y0.shape != (n,) or y1.shape != (n,) or y2.shape != (n,):
        return _no('radau_tc:seed')
    ## the two factorisations as `_FrozenTransform.solve` would use them
    lu, zs, prep = fz.lu, fz.zs, fz.prep
    _ZS, _LU = _MOD['mods'][3], _MOD['mods'][2]
    sk = None
    if lu is None:
        ## (no dense LU kept: the analysis solver's StructuralKLU, `_real_klu`)
        sk = _real_klu(fz.ls, m)
        if sk.__class__ is str:
            return _no('radau_tc:frozen')
    elif type(lu) is not _LU or lu.n != m:
        return _no('radau_tc:frozen')
    if (type(zs) is not _ZS or prep[1] != m or 'solve' in zs.__dict__
            or 'solve_prepared' in zs.__dict__):
        return _no('radau_tc:frozen')
    if lu is not None and lu._info is not None and lu._info != 0:
        ## (singular: the Python raises numpy's error at its first solve)
        return _no('radau_tc:lu')
    if not (zs._symbolic and zs._numeric and prep[5] == zs._pattern):
        ## (the record's pattern not analysed and factored: the Python does it)
        return _no('radau_tc:klu')
    has_limiter = tr._newton_limiter() is not None
    if has_limiter:
        Tw = getattr(tr.epar, 'T', 300.0)
        if type(Tw) is not float and type(Tw) is not int and type(Tw) is not np.float64:
            return _no('radau_tc:T')
    rec = td.get('_radau_tc')
    if rec is None or rec.core is not core or rec.n != n or rec.iref != iref:
        _MOD['bits'] = tc._PASS_BITS['q'] | tc._PASS_BITS['i']
        cffi_core, cfn_core, _dg = tc.driver()
        rec = _Ctx(core, n, iref, ffi, cfn_core, cffi_core, drv)
        td['_radau_tc'] = rec
    s, buf = rec.s, rec.buf
    rec.served = None
    if has_limiter:
        if not NC._walk_ready(tr, rec):
            return _no('radau_tc:limit')
    elif s.wn:
        s.walk_fn = ffi.NULL
        s.wn = 0
        rec.walk = rec.walk_ok = None
    Tc = core.probe(tr.epar)
    if Tc.__class__ is str:
        return _no('radau_tc:ready')
    zk = _MOD.get('klu')
    if zk is None or zk[0] is not zs._lib:
        zk = _MOD['klu'] = (zs._lib, _cbuild.address(zs._lib.klu_z_refactor),
                            _cbuild.address(zs._lib.klu_z_solve))
    ## ---- from here the attempt leaves one trace: the step's source memo ----
    analysis = tr.par.analysis
    keys = [(tt, analysis) for tt in ctx.tstage]
    try:
        had = [key in memo for key in keys]
    except TypeError:
        return _no('radau_tc:u')
    u_keys = M['u_keys']
    cvals = tuple(map(_PC.get, u_keys, itertools.repeat(0)))
    try:
        us = [src(tt) for tt in ctx.tstage]
    except Exception:                                          # noqa: BLE001
        ## (the call raises again in the Python loop, where it raised before)
        NC._undo_source(memo, keys, had, dict(zip(u_keys, cvals)))
        return _no('radau_tc:u')
    u0, u1, u2 = us
    if not (type(u0) is np.ndarray and type(u1) is np.ndarray and type(u2) is np.ndarray
            and (u0.dtype is _F64 or u0.dtype == np.float64)
            and (u1.dtype is _F64 or u1.dtype == np.float64)
            and (u2.dtype is _F64 or u2.dtype == np.float64)
            and u0.shape == (n,) and u1.shape == (n,) and u2.shape == (n,)):
        NC._undo_source(memo, keys, had, dict(zip(u_keys, cvals)))
        return _no('radau_tc:u')
    buf['u'][...] = us
    buf['seed'][...] = seed
    np.copyto(buf['qn'], qn)
    ## (the tableau once its values change; `P`'s and `V`'s entries once the
    ## step's cached tuple changes -- `_radau_frozen` keeps it per step size)
    ab = Amat.tobytes()
    if ab != rec.A_bytes:
        for k, v in enumerate(Amat.ravel().tolist()):
            s.A[k] = v
        rec.A_bytes = ab
    if fz.p is not rec.p_of:
        for k, z in enumerate(fz.p):
            s.p[2 * k] = float(z.real)
            s.p[2 * k + 1] = float(z.imag)
        for i, (a, b) in enumerate(fz.v):
            s.v[3 * i] = float(a)
            s.v[3 * i + 1] = float(b.real)
            s.v[3 * i + 2] = float(b.imag)
        rec.p_of = fz.p
    ## (the kept LU by its addresses -- a refilled one keeps them; the
    ## pattern's arrays by their identity -- a kept pattern keeps them; or
    ## the real factor on KLU, `_rk_bind`)
    if sk is not None:
        if not _rk_bind(rec, sk, fz):
            NC._undo_source(memo, keys, had, dict(zip(u_keys, cvals)))
            return _no('radau_tc:klu_foreign')
    else:
        s.rk_on = 0
        addr = lu.addr
        if addr != rec.lu_addr:
            s.lu_a = ffi.cast('double *', addr[0])
            s.lu_ipiv = ffi.cast('int64_t *', addr[1])
            rec.lu_addr = addr
        s.lu_info = -1 if lu._info is None else 0
    Ap, Ai = prep[2], prep[3]
    if Ap is not rec.Ap_obj or Ai is not rec.keep_ap[1]:
        rec.keep_ap = (Ap, Ai, ffi.from_buffer('int *', Ap), ffi.from_buffer('int *', Ai))
        s.Ap, s.Ai = rec.keep_ap[2], rec.keep_ap[3]
        rec.Ap_obj = Ap
    ax = ffi.from_buffer('double *', prep[4])
    s.Ax = ax
    kk = (zk[1], zk[2], zs._symbolic, zs._numeric, ctypes.addressof(zs._common))
    if kk != rec.klu_kk:
        s.klu_refactor = ffi.cast('void *', kk[0])
        s.klu_solve = ffi.cast('void *', kk[1])
        s.Symbolic = ffi.cast('void *', kk[2])
        s.Numeric = ffi.cast('void *', kk[3])
        s.Common = ffi.cast('void *', kk[4])
        rec.klu_kk = kk
    fresh = zs._fresh is prep
    s.klu_fresh = 1 if fresh else 0
    s.T = Tc
    s.h, s.reltol, s.abstol = float(h), float(reltol), float(abstol)
    s.rtol, s.maxiter = float(zs.REFACTOR_RESIDUAL_TOL), int(maxit)
    ## THE STEP'S END PASSES (`END`): asked where `_stage_end_passes`
    ## would make them -- its switch on and the step's memo rolling (it
    ## reads and records only there); `C` and `G` at the first two stages
    ## where the shooting reads them (`_stage_G_read`)
    end = END and _MOD['mods'][0].STAGE_FUSE and bool(td.get('_memo_rolling'))
    s.end_bits = 7 if end else 0
    s.end_cg = 3 if end and td.get('_stage_G_read') else 0
    status = cfn(s)
    return _outcome(tr, rec, fz, status, probe, 1,
                    lambda: NC._undo_source(memo, keys, had, dict(zip(u_keys, cvals))))


def _real_klu(ls, m):
    """The real factor's solver where `_radau_frozen` keeps no dense LU: the
    analysis solver `ls` an AutoSolver whose choice is the circuit's
    `StructuralKLU` -- each piece `_FrozenTransform.real_solve` calls as
    defined and unshadowed, the solver analysed and of this size -- or why
    not (`klu_real`)."""
    g = _MOD.get('klu_real')
    if g is None:
        from pycircuit.circuit import _tran_radau
        from pycircuit.circuit.linearsolver import AutoSolver, StructuralKLU
        L = 'pycircuit.circuit.linearsolver'
        own = ((_tran_radau._FrozenTransform.real_solve, '_FrozenTransform.real_solve',
                'pycircuit.circuit._tran_radau'),
               (AutoSolver.solve, 'AutoSolver.solve', L),
               (AutoSolver._select, 'AutoSolver._select', L),
               (StructuralKLU.solve, 'StructuralKLU.solve', L),
               (StructuralKLU._ours, 'StructuralKLU._ours', L),
               (StructuralKLU._factor, 'StructuralKLU._factor', L),
               (StructuralKLU._analysed, 'StructuralKLU._analysed', L),
               (StructuralKLU._free, 'StructuralKLU._free', L))
        if not all(_paths.genuine(*o) for o in own):
            return 'klu_real'
        g = _MOD['klu_real'] = (_tran_radau, AutoSolver, StructuralKLU,
                                tuple(o for o, _q, _m in own))
    tr_, AS, SK, pieces = g
    if (tr_._FrozenTransform.real_solve, AS.solve, AS._select, SK.solve, SK._ours, SK._factor,
            SK._analysed, SK._free) != pieces:
        return 'klu_real'
    if type(ls) is not AS or 'solve' in ls.__dict__ or '_select' in ls.__dict__:
        return 'klu_real'
    sk = ls._choice
    if (type(sk) is not SK or sk.m != m or not sk._symbolic
            or any(k in sk.__dict__ for k in ('solve', '_ours', '_factor', '_analysed', '_free',
                                              '_other'))):
        return 'klu_real'
    return sk


def _rk_bind(rec, sk, fz):
    """The real factor on KLU bound for this call (radau_tc.h's `rk`): the
    factor's values at the solver's pattern (`_ours`'s take, once: it is
    frozen) -- False where a nonzero lies off it (`_ours` hands that to the
    fallback solver) -- the residual's dense work, the solver's pattern,
    analysis, common, library entries, and its numeric and tolerance now."""
    rf = fz.rf
    if not (type(rf) is np.ndarray and rf.dtype == np.float64 and rf.flags.c_contiguous):
        return False
    ax = rf.ravel().take(sk.take)
    if np.count_nonzero(rf) != np.count_nonzero(ax):
        return False
    ffi, s, m = rec.ffi, rec.s, rec.n - 1
    k = s.rk
    keep = rec.rk_keep
    if keep is None or keep[0] is not sk:
        lib = sk._lib
        work = np.empty(max(m * m, 1))
        keep = rec.rk_keep = [sk, ffi.from_buffer('int *', sk.Ap), ffi.from_buffer('int *', sk.Ai),
                              None, work, ffi.from_buffer('double *', work)]
        k.kAp, k.kAi = keep[1], keep[2]
        s.rk_A = keep[5]
        k.klu_factor, k.klu_refactor, k.klu_solve, k.klu_free = (
            ffi.cast('void *', _cbuild.address(f))
            for f in (lib.klu_factor, lib.klu_refactor, lib.klu_solve, lib.klu_free_numeric))
        k.klu_common = ffi.cast('void *', sk._ct.addressof(sk._common))
    keep[3] = ax
    k.resid_tol = float(sk.REFACTOR_RESIDUAL_TOL)
    k.kAx = ffi.from_buffer('double *', ax)
    k.klu_sym = ffi.cast('void *', sk._symbolic)
    num = sk._numeric
    k.klu_num = ffi.cast('void *', num) if num else ffi.NULL
    s.rk_on = 1
    return True


def _rk_back(rec, sk):
    """The analysis solver as the C's real solves leave it, as its own
    `solve` would: the numeric the sequence ends on, the counts."""
    k = rec.s.rk
    sk._numeric = int(rec.ffi.cast('uintptr_t', k.klu_out)) or None
    sk.refactors += int(k.refactors)
    sk.factors += int(k.factors)
    sk.pivot_fallbacks += int(k.pivot_fallbacks)
    sk.residual_fallbacks += int(k.residual_fallbacks)


def _outcome(tr, rec, fz, status, probe, first, undo):
    """What a call of the C (`solve`, `resume`) made of the loop, as `solve`
    returns it: the factorisations' state written back, the counts, the
    converged stages and the end passes.  `first`: the call's iterations
    whose source values were taken before it (1 for `solve`'s first
    iteration; 0 for `resume`, the step's later iterations all); `undo`:
    the source memo's rollback on a hand-back, or None."""
    s, buf, n = rec.s, rec.buf, rec.n
    lu, zs, prep = fz.lu, fz.zs, fz.prep
    ## the factorisations' state, as the Python would have left it: the
    ## kept LU after its first `dgesv` -- or the real factor on KLU, its
    ## sequence's numeric and counts where it stands (a bail took it back) --
    ## the numeric holding the record's refactor (a failed refactor is the
    ## Python's to make and count again)
    if s.rk_on:
        if status in (_CODES['RADAU_TC_CONVERGED'], _CODES['RADAU_TC_MAXITER']):
            _rk_back(rec, fz.ls._choice)
    elif s.lu_info >= 0:
        lu._info = int(s.lu_info)
    if s.refactored:
        zs.refactors += int(s.refactored)
        zs._fresh = prep
    if status == _CODES['RADAU_TC_MAXITER']:
        ## (the Python loop's False: its iterations made, the source calls
        ## after the first its memo's hits)
        it = int(s.iters)
        _PC['umemo:hit'] += 3 * (it - first)
        _PC['radau_tc:iterations'] += it
        if probe:
            _PC['radau_tc:probe'] += 1
            return False
        _PC['radau_tc:served'] += 1
        _PC['radau_tc:noconv'] += 1
        ## (its checks made at this watch counter all the same: the C
        ## stepping loop may enter after the dense fallback, `_tran_tcloop_c`)
        rec.served = _watch.now()
        return False
    if status != _CODES['RADAU_TC_CONVERGED']:
        if undo is not None:
            undo()
        if status == _CODES['RADAU_TC_BAIL_WALKSTOP']:
            ## (the walk stopped: its tables are taken again next time)
            rec.walk_ok = None
        return _no('radau_tc:bail:' + BAIL.get(status, str(status)))
    it = int(s.iters)
    _PC['umemo:hit'] += 3 * (it - first)
    _PC['radau_tc:served'] += 1
    _PC['radau_tc:iterations'] += it
    ## (the stamp the C stepping loop enters on: its checks made at this
    ## watch counter, `_tran_tcloop_c`)
    rec.served = _watch.now()
    td = tr.__dict__
    Y = buf['Y']
    out = [Y[j].copy() for j in range(3)]
    if s.end_done:
        ## (fresh arrays, as `passes` returns them, keyed in its order; the
        ## stage list itself the token: `_stage_end_passes` takes them for
        ## these stages only, with the watch counter unmoved)
        nn = (n, n)
        last = {'q': buf['e_q'].copy(), 'i': buf['e_i'].copy(),
                'C': buf['e_C'].reshape(nn).copy(), 'G': buf['e_G'].reshape(nn).copy()}
        cg = None
        if s.end_cg:
            C01, G01 = buf['e_C01'], buf['e_G01']
            cg = [{'C': C01[j * n * n:(j + 1) * n * n].reshape(nn).copy(),
                   'G': G01[j * n * n:(j + 1) * n * n].reshape(nn).copy()} for j in (0, 1)]
        td['_tc_end'] = (out, _watch.now(), last, cg, rec.core)
        _PC['radau_tc:end'] += 1
    return out


def iterate(tr):
    """The iterate a `solve(..., probe=True)` on `tr` that returned False
    reached -- the loop's after its iteration, bit for bit the Python
    loop's -- as the three rows of the C's buffer: read them before the
    next call."""
    Y = tr.__dict__['_radau_tc'].buf['Y']
    return [Y[0], Y[1], Y[2]]


def resume(tr, fz, maxit):
    """The loop after a `solve(..., probe=True)` on this step returned False,
    resumed from the iterate it reached for `maxit` more iterations: what
    one `solve` call of ``maxit + 1`` would have returned -- the converged
    stages, or False -- or None, handed back, and the Python loop runs from
    the iterate (`seed`).  Nothing the C reads moved between the two calls:
    the early abort's check between them (`_RadauStages._radau_early_check`)
    reads the iterate, its source call at the step's start is another key
    of the step's memo, and a device pass at ``x_n`` writes the core's own
    scratch, which the C's calls write before they read; the stages' source
    values are in the buffer since the probe, so no source call is made
    and none is rolled back."""
    rec = tr.__dict__['_radau_tc']
    s, buf = rec.s, rec.buf
    np.copyto(buf['seed'], buf['Y'])
    if s.rk_on:
        num = fz.ls._choice._numeric
        s.rk.klu_num = rec.ffi.cast('void *', num) if num else rec.ffi.NULL
    else:
        s.lu_info = -1 if fz.lu._info is None else 0
    s.klu_fresh = 1 if fz.zs._fresh is fz.prep else 0
    s.maxiter = int(maxit)
    return _outcome(tr, rec, fz, _driver[1](s), False, 0, None)


def seed(tr):
    """The stages a `resume` that handed back started from (the probe's
    iterate), as fresh full-width arrays."""
    Y = tr.__dict__['_radau_tc'].buf['seed']
    return [Y[j].copy() for j in range(3)]
