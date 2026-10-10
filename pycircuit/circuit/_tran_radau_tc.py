"""RADAU'S TRANSFORM NEWTON IN C (speed round 10, B3.3; 2026-10-05).

The cost transform (`_RadauStages._rk_step_transformed`: PSS's default on
PSP-class circuits under 'auto') is a SIMPLIFIED Newton -- `C` and `G`
frozen at ``x_n``, the coupled solve through ``A^{-1}``'s eigenbasis, its
factors made once a step (`_radau_frozen`, B3.1).  Its loop spent ~3.5 M
instructions a step on the PSP stage's radau PSS in Python and numpy
around ~0.4 M of device kernels: per iteration the three stages' passes
through `_tran_core.passes`, the residual (`_coupled_stage_system`), the
frozen solve (`_FrozenTransform.solve`: the right-hand sides, numpy's kept
LU, the complex factor's KLU solve and its residual check, the update),
and per stage the step, the limiter and the convergence test.  Here ONE C
call runs the loop: per iteration the stages' ``q`` and ``i`` through the
evaluate core's own C (formula 4, as `passes`), ``K_j = -(i_j + u_j)``,
``F_i = (q_i - q_n) - h (((0 + A_i0 K_0) + A_i1 K_1) + A_i2 K_2)`` reduced
(`_tran_radau_c`'s assembly, without the blocks); the right-hand sides
``-((P_k0 F0 + P_k1 F1) + P_k2 F2)`` -- a complex scalar times a real
vector, as numpy multiplies it: the vector cast to complex, numpy's
complex product; numpy's kept LU (`linearsolver.NumpyLU`'s protocol: its
OpenBLAS's `dgesv` at the first solve, `dgetrs` on the kept factors
after); the complex factor's KLU on the solver's own handles (the refactor
`solve_prepared` would make, at most once a record, then `klu_z_solve`);
``dY_i = V_i0 w0 + 2 Re(V_i1 w1)``; per stage ``Y_prev + insert(dY_i,
iref, 0.0)``, the limiter walk's own driver and the limited step's largest
reduced entry; the convergence test after the update, as the Python tests.
Bit for bit the Python loop (pinned in `tests/test_pinned_pairs.py`).

⚠ NUMPY'S COMPLEX PRODUCT IS READ, NOT ASSUMED.  On this box (numpy's
X86_V3 loops: AVX2 and FMA3) it FUSES -- ``re = fma(ar, br, -(ai*bi))``,
``im = fma(ar, bi, ai*br)``, at every position of every length -- and a
build or CPU without FMA multiplies separately.  `cmul_mode` asks numpy
itself once, on inputs where the two forms round apart (and on signed
zeros and underflow, where only the fused form's sign rules hold), and
the C takes the form numpy showed; neither: the path is off.

⚠ THE RESIDUAL CHECK IS DECIDED WITH A MARGIN.  `solve_prepared` falls
back to a fresh factor where ``max|A x - b| / max|b|`` exceeds its
tolerance (1e-8), read with numpy's complex `abs` (its own algorithm, not
`hypot`) of SciPy's product.  The C computes the same quantities in its own
rounding and decides only outside the band their difference can reach
(``32 (m + 4) eps`` of ``max_k sum_j |A_kj| |x_j|`` over ``max|b|``, plus
``16 eps`` of the residual): a residual nearer the tolerance than that, or
past it, hands the iteration back and the Python decides, as it does --
in practice the residual sits ~1e-15 below a 1e-8 tolerance.

TRANSACTIONAL.  The C writes nothing outside its buffers but the two
factorisations' state, which the Python would have left the same: numpy's
kept LU after its first `dgesv` (`NumpyLU._info`), and the KLU numeric
holding a refactor of the step's record (`ComplexKLUSolver._fresh`,
`refactors`).  A refactor of the same values with the same pivots repeats
its bits and `dgetrs` on `dgesv`'s factors is `dgesv`'s solve, so a Python
loop run after a hand-back makes the C's iterates again.  The step's
source memo holds the three stage times (their calls made here, before
the C) and its hit count for every later iteration's.  Anything unusual --
a non-finite value, a floating-point exception in its own arithmetic
(numpy would warn), a non-finite real right-hand side (the analysis
solver's), a singular real factor (numpy's error), the KLU refactor or
solve failing, the residual check not decided or failing (the fallback's
fresh factor), the walk stopping, maxiter -- hands the loop back: the
source memo's entries and counts rolled back, and the Python loop runs
from the same seed (the same answer, exception, warning).

THE STEP'S END PASSES (speed round 12, stage 5a; `END`, env
`PYCIRCUIT_RADAU_TC_END`).  Once converged the C makes, in the same call,
the passes `_stage_end_passes` makes next: all four at the last stage, and
`C` and `G` at the first two where the shooting reads them (`_stage_G_read`)
-- each `_tran_core.passes`'s core call with its bits, outputs and dummies,
in its order (the core's own arrays are its scratch, each call writing what
it reads: a call the Python then does not need leaves nothing behind) -- and
hands them over with the very stage list it returns (`_tc_end`, the watch
counter with it).  `_stage_end_passes` takes them for those stages only,
the counter unmoved, all four passes wanted at the last stage and
`passes`'s readiness holding as stamped (`_tran_core.passes_held`), making
the counters each of its calls makes (`passes_count`); else it makes the
calls (`radau.fuse:unhanded`).  Each hand-off's own `passes` call cost
~57 k instructions of readiness and marshalling around its kernels.

THE FROZEN FACTORS (speed round 12, stage 5b; `FOLD`, env
`PYCIRCUIT_RADAU_FOLD`).  `_radau_frozen` makes its two factors -- ``(lam0/h)
Cr + Gr`` real, ``(lam1/h) Cr + Gr`` complex -- and the complex one's packed
values by `fold`, one small C call of its own: numpy's operations in its
order (two roundings for the real one; for the complex one ``Cr`` cast to
complex, numpy's product in `cmul_mode`'s form, then ``Gr`` cast to
complex), the four flags its `errstate(all='raise')` raises on, and the
nonzeros where they fall in the last record's pattern (`_csc_of_dense`'s
order: column by column, rows ascending).  Every decision stays in
`_radau_frozen`, in its order: a flag declines as numpy's error does, the
kept LU is made as before, the record is `prepare`'s own
(`ComplexKLUSolver.prepare_values`) where the values fell in the pattern and
`prepare`'s otherwise.  A non-finite input is numpy's (a NaN's payload
follows its loops).  Its own call, not inside the transform's: inside, a
factor that would raise surfaces only after the call's readiness checks
have run, and the step would be made again around them.  The C's work is
nothing on the PSP stage's 6x6 factors; the call's marshalling is the cost
(buffers kept per size: `ffi.from_buffer` costs ~3.5 k a call).

DECLINES, before anything is touched (`_paths.no`): the switch (and the
evaluate core's two), no frozen factors or no kept LU or marshalled record
(`_radau_frozen`'s own), no build, numpy's complex product of neither
form, a toolkit that is not numeric, a subclass of the transient or a
patched piece of the loop, an instance shadow of the circuit's passes or
limiter or of the transient's pieces, a `provided_function`, a stateful
limiter or circuit-level limiting, a bypass, no source memo, a circuit the
core cannot serve, a limiting element that is not its C kernel, a batch
not ready, a KLU solver whose numeric is not of the record's pattern (the
Python analyses and factors), a singular kept LU, inputs of another type
or size.

`ENABLED` (env `PYCIRCUIT_RADAU_TC=0`) is read on every call; `STATUS` says
why the path is off where it is.  History: `doc/transient_history.md`,
`_tran_radau_tc.py`.
"""
import ctypes
import itertools
import math
import os

import numpy as np

from pycircuit.circuit import _cbuild, _paths, _watch

_PC = _paths.COUNTS
#: (float64's dtype: a check by identity first, `==` where it is not the
#: object -- speed round 11)
_F64 = np.dtype(np.float64)
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_RADAU_TC', '1') != '0'
#: THE STEP'S END PASSES in the same call (speed round 12, stage 5a): the
#: passes `_stage_end_passes` makes next, made by the C after it converges
#: and handed over (`_tc_end`); env `PYCIRCUIT_RADAU_TC_END=0` turns it off
END = os.environ.get('PYCIRCUIT_RADAU_TC_END', '1') != '0'
#: THE FROZEN FACTORS IN C (speed round 12, stage 5b): `_radau_frozen`'s two
#: factors and the complex one's packed values made by one small C call
#: (`fold`); env `PYCIRCUIT_RADAU_FOLD=0` turns it off
FOLD = os.environ.get('PYCIRCUIT_RADAU_FOLD', '1') != '0'
STATUS = 'not loaded'

RADAU_TC = _cbuild.csrc('radau_tc.c')
RADAU_TC_CDEF = _cbuild.csrc('radau_tc.h')

FOLD_C = _cbuild.csrc('radau_fold.c')
FOLD_CDEF = _cbuild.csrc('radau_fold.h')

#: the C's status codes past 1 (converged): why it handed the loop back
BAIL = {2: 'nonfinite', 3: 'flags', 4: 'rhs', 5: 'lu', 6: 'refactor', 7: 'klusolve',
        8: 'residual', 9: 'step', 10: 'walkstop', 11: 'maxiter'}

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
        from pycircuit.circuit import _hdl_cbackend as cb
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
            ffi, cfn, _key, _cold, _secs = cb.load_kernel(RADAU_TC, RADAU_TC_CDEF)
        except (cb.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _ct, gesv, getrs = lap
        from pycircuit.circuit import _tran_newton_c, _tran_radau_c
        _MOD['NC'], _MOD['RC'] = _tran_newton_c, _tran_radau_c
        _driver = (ffi, cfn, ctypes.cast(gesv, ctypes.c_void_p).value,
                   ctypes.cast(getrs, ctypes.c_void_p).value)
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
        from pycircuit.circuit import _hdl_cbackend as cb
        try:
            ffi, cfn, _key, _cold, _secs = cb.load_kernel(FOLD_C, FOLD_CDEF)
        except (cb.CompileError, OSError) as e:
            _FOLD['status'] = f'off ({e})'
            return None
        d = _FOLD['drv'] = (ffi, cfn, ffi.new('radau_fold_t *'))
        d[2].fused = _MOD['cmul']
        _FOLD['status'] = 'c'
    return d or None


def _fold_bufs(m):
    """`fold`'s buffers for size `m` -- the inputs copied in, the outputs
    copied out: a buffer's handle is made once (`ffi.from_buffer` costs
    ~3.5 k instructions a call, more than the C's work on a small factor)
    -- set in the struct."""
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
    """`_RadauStages._radau_frozen`'s two factors from ONE C call (speed
    round 12, stage 5b), or None where the C does not make them -- the
    switch, no build, inputs other than two float64 square arrays of one
    size, a non-finite entry (each counted, `radau.frozen:fold.<why>`): the caller
    makes them in numpy, as before.  Else `(flagged, rf, cf, data)`:
    `flagged`, an operation raised a flag numpy raises on under
    `errstate(all='raise')` (the caller declines as it does; nothing else
    returned); `rf` and `cf`, ``s0 * Cr + Gr`` and ``c1 * Cr + Gr`` -- fresh
    arrays, numpy's bytes, for the scalars as the caller's expression makes
    them, finite; `data`, the complex factor's nonzeros where they fall in
    `last`'s pattern (`ComplexKLUSolver._csc_last`), in its order
    (`ComplexKLUSolver.prepare_values`) -- None where they fall elsewhere,
    there is no pattern or it has no nonzeros.  It stands in for the numpy calls, their errstate,
    the mask and the gather: ~85 k instructions a step on the PSP stage's
    radau PSS."""
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
    if f.status == 2:
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
            ## (the step's end passes' outputs, stage 5a)
            'e_q': np.empty(n), 'e_i': np.empty(n), 'e_C': np.empty(n * n),
            'e_G': np.empty(n * n), 'e_C01': np.empty(2 * n * n), 'e_G01': np.empty(2 * n * n),
            'e_q01': np.empty(n),
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
        s.walk_fn = ffi.NULL
        s.wn = 0
        ## (what does not change from step to step is set once: the product
        ## form, the passes' bits; the rest when it changes -- speed round 10,
        ## B3.5: a cast from `arr.ctypes.data` cost 13.5 k instructions, 21
        ## of `P`'s and `V`'s entries ~30 k a step)
        s.fused = _MOD['cmul']
        s.bits = _MOD['bits']
        self.A_bytes = self.p_of = self.Ap_obj = self.keep_ap = self.klu_kk = None
        self.lu_addr = None
        self.walk = self.walk_ok = None
        self.walk_els = self.walk_dicts = self.walk_packs = self.walk_kern = ()
        self.walk_hw_els = self.walk_hwk = ()
        self.walk_stamp, self.walk_arms, self.walk_held = -1, 0, 0
        self.walk_watch = ()


def solve(tr, ctx, fz, seed, src, provided_function, lims, nobypass, reltol, abstol, maxit):
    """The transform's Newton loop from the stages `seed` in C: the three
    converged stages as fresh full-width arrays, the source memo and the
    two factorisations' state as the Python loop leaves them -- or None:
    declined or handed back, nothing left behind that the Python loop
    would not make again, and the Python loop runs."""
    if not ENABLED:
        return _no('radau_tc:off')
    if fz is None or fz.lu is None or fz.prep is None:
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
    if (type(lu) is not _LU or type(zs) is not _ZS or lu.n != m or prep[1] != m
            or 'solve' in zs.__dict__ or 'solve_prepared' in zs.__dict__):
        return _no('radau_tc:frozen')
    if lu._info is not None and lu._info != 0:
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
        zk = _MOD['klu'] = (zs._lib, ctypes.cast(zs._lib.klu_z_refactor, ctypes.c_void_p).value,
                            ctypes.cast(zs._lib.klu_z_solve, ctypes.c_void_p).value)
    ## ---- from here the attempt leaves one trace: the step's source memo ----
    RC = _MOD['RC']
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
        RC._undo_source(memo, keys, had, dict(zip(u_keys, cvals)))
        return _no('radau_tc:u')
    u0, u1, u2 = us
    if not (type(u0) is np.ndarray and type(u1) is np.ndarray and type(u2) is np.ndarray
            and (u0.dtype is _F64 or u0.dtype == np.float64)
            and (u1.dtype is _F64 or u1.dtype == np.float64)
            and (u2.dtype is _F64 or u2.dtype == np.float64)
            and u0.shape == (n,) and u1.shape == (n,) and u2.shape == (n,)):
        RC._undo_source(memo, keys, had, dict(zip(u_keys, cvals)))
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
    ## pattern's arrays by their identity -- a kept pattern keeps them)
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
    ## THE STEP'S END PASSES (stage 5a): asked where `_stage_end_passes`
    ## would make them -- its switch on and the step's memo rolling (it
    ## reads and records only there); `C` and `G` at the first two stages
    ## where the shooting reads them (`_stage_G_read`)
    end = END and _MOD['mods'][0].STAGE_FUSE and bool(td.get('_memo_rolling'))
    s.end_bits = 7 if end else 0
    s.end_cg = 3 if end and td.get('_stage_G_read') else 0
    status = cfn(s)
    ## the factorisations' state, as the Python would have left it: the
    ## kept LU after its first `dgesv`, the numeric holding the record's
    ## refactor (a failed refactor is the Python's to make and count again)
    if s.lu_info >= 0:
        lu._info = int(s.lu_info)
    if s.refactored:
        zs.refactors += int(s.refactored)
        zs._fresh = prep
    if status != 1:
        RC._undo_source(memo, keys, had, dict(zip(u_keys, cvals)))
        if status == 10:
            ## (the walk stopped: its tables are taken again next time)
            rec.walk_ok = None
        return _no('radau_tc:bail:' + BAIL.get(status, str(status)))
    it = int(s.iters)
    _PC['umemo:hit'] += 3 * (it - 1)
    _PC['radau_tc:served'] += 1
    _PC['radau_tc:iterations'] += it
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
        td['_tc_end'] = (out, _watch.now(), last, cg, core)
        _PC['radau_tc:end'] += 1
    return out
