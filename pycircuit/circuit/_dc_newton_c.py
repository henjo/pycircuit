"""THE DC OPERATING POINT'S BASE NEWTON SOLVE IN C: one
`StandardNewton.solve_system` of the DC's continuation chain -- on the
evaluation its ladders compose -- as one call of `csrc/dcnewton.c`.

The operating point is found by Newton solves on `dcanalysis`'s own
residual (`func`: `cut(i + u)`; source stepping's `source_callback`:
`cut(i + lam * u)`; `cut(G)`), each under the ladders' deformations: gshunt
(`nrsolver.GminSteppingNewton`: `F + x g`, `J + 0.0`, g on the diagonal)
and junction gmin over it (`JunctionGminSteppingNewton`: the per-junction
terms).  Here ONE C call makes the whole solve: the C core's G and i
(`_tran_core`), the residual and the layers, the solve -- numpy's own
dgesv, or the run's StructuralKLU in its own sequence -- the step, the
limiter walk (`_hdl_climit`'s, set up by `_tran_newton_c._walk_ready`) and
nrsolver's convergence test, its rounding floor's too.  Bit for bit `StandardNewton.solve_system` on
that evaluation (pinned in `tests/test_pinned_pairs.py`).

THE COMPOSITION.  The wrappers say what they are: each layer's closure its
`dc_layer` -- `('gshunt', g, inner)`, `('jgmin', g, inner, rows)`,
`('source', lam, inner)` -- the base (`dcanalysis.refnode_removed`'s) its
`dc_base` `(func, irefnode, reduced_G)`, the DC's own functions their
`dc_owner` `(the DC analysis, 'func' or 'source')` and its limiter its
`dc_owner`, the analysis.  Served: junction gmin over gshunt over source
stepping over the base, any layer absent, the source layer exactly where
the base is the source callback.  Anything else -- the pseudo-transient
ladder, the anchor, another analysis's evaluation -- is not a DC
composition and runs as it did, uncounted.

OWNS, on the DC analysis: `_dc_newton_c` (a `_Rec`: the struct, its
buffers, the core, the KLU binding, the walk).  A serve writes, on the
run's StructuralKLU, its `_numeric` and its `refactors`, `factors`,
`pivot_fallbacks` and `residual_fallbacks` as its own solve would have.
`_paths` prefix `dcnewton_c:`: `dcnewton_c:served`, `dcnewton_c:iterations`,
`dcnewton_c:maxiter`, the declines and the hand-backs
(`dcnewton_c:bail:<why>`).  Switch: `ENABLED` (env `PYCIRCUIT_DCNEWTON_C=0`),
read on every call.

THE PACKED LAYOUT (`_packed`).  Inside the scope of its circuit the record
holds no dense matrix: its own core (`_tran_core._Core(packed=True)`: no
n x n bin) makes G in its packed bin, J is made at KLU's entries alone, the
layers at theirs, and the convergence scale and the refactor's residual by
the reference's CSC kernel on KLU's arrays -- the Python's on the pattern,
bit for bit; KLU on the structure whatever the size (`AutoSolver.
choose_packed`, `_select`'s own choice for a PatternMatrix).

EXACT OR NOT AT ALL.  A flag numpy would report, a non-finite value, a
failed factor or solve, or a limiter the walk stops at hands the solve back
-- the numeric StructuralKLU held at entry kept, nothing of the analysis's
written -- and the Python's solve runs from x0: its warnings, its
exceptions, its refactors on the same pivot order.  The iteration limit is
served: the last iteration's arrays, and the Python's own exception
(`nrsolver._standard_failure`).

DECLINES (`_paths.no`): before this module is imported (`nrsolver._dc_solve`)
a solve below `nrsolver.DC_C_MIN_N` unknowns and the run's first (its
solver not chosen yet) below its structural floor; here a first solve the
structure does not choose for (`first_choice`), the switch, no build, the rounding floor's
helpers patched, a scaler that is not NoneScaler, a limiter that is not the DC's, a linear solver that is not
the DC's AutoSolver over its dense solve or an analysed StructuralKLU whose
pattern holds the core's, a toolkit that is not numeric, numpy not ignoring
underflow, the circuit-level limiter or a limiting element the walk cannot
take, a circuit the core cannot serve or its readiness (`_Core.probe`), x0
and the tolerances not float64 vectors of the reduced size; on the pattern
a solve no structure chooses KLU for, a junction entry off KLU's pattern;
inside another circuit's packed scope (`dcnewton_c:packed`).

`STATUS` says why the path is off where it is."""
import os

import numpy as np

from pycircuit.circuit import _cbuild, _packed, _paths

ENABLED = os.environ.get('PYCIRCUIT_DCNEWTON_C', '1') != '0'
STATUS = 'not loaded'

_SOURCES = None


def _sources():
    """`(the unit's source, its declaration, its codes, its bails' names,
    pcx_abi.h's codes)`, read at the first need: importing the module reads
    no file (it is imported with nrsolver, inside a run's first operating
    point, where a file read is 0.35 ms)."""
    global _SOURCES
    if _SOURCES is None:
        c = _cbuild.csrc('dcnewton.c')
        h = _cbuild.csrc('pcx_klu.h') + _cbuild.csrc('dcnewton.h')
        codes = _cbuild.constants(h)
        _SOURCES = (c, h, codes, _cbuild.bails(codes, 'DCNEWTON_BAIL_'),
                    _cbuild.constants(_cbuild.csrc('pcx_abi.h')))
    return _SOURCES


def __getattr__(name):
    """`DCNEWTON_C` and `DCNEWTON_CDEF`, the unit's source and declaration
    (`_sources`)."""
    if name == 'DCNEWTON_C':
        return _sources()[0]
    if name == 'DCNEWTON_CDEF':
        return _sources()[1]
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

_no = _paths.no
_PC = _paths.COUNTS
_F64 = np.dtype(np.float64)
_driver = None


def driver():
    """`(ffi, cfn, dgesv, dgemv, ddot)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            npl = _cbuild.openblas(np)
            if npl is None:
                STATUS = "off (numpy's own OpenBLAS not found)"
                return None
            dgesv = _cbuild.address(npl.scipy_dgesv_64_)
            dgemv = _cbuild.address(npl.scipy_cblas_dgemv64_)
            ## (the rounding floor's: np.dot of two vectors)
            ddot = _cbuild.address(npl.scipy_cblas_ddot64_)
        except (OSError, AttributeError) as e:
            STATUS = f'off ({e})'
            return None
        try:
            src, cdef = _sources()[:2]
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(src, cdef,
                                                             'pcx_dcnewton')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, dgesv, dgemv, ddot)
        STATUS = 'c'
    return _driver or None


def composition(eval_FJ):
    """The DC composition `eval_FJ` is -- `(the DC analysis, lam or None,
    gshunt's g or None, (junction gmin's g, its rows) or None)` -- or None
    where it is not one (THE COMPOSITION, above)."""
    layers, order, f = {}, [], eval_FJ
    while True:
        lay = getattr(f, 'dc_layer', None)
        if lay is None:
            break
        if lay[0] in layers:
            return None
        layers[lay[0]] = lay
        order.append(lay[0])
        f = lay[2]
    base = getattr(f, 'dc_base', None)
    owner = getattr(base[0], 'dc_owner', None) if base is not None else None
    if owner is None:
        return None
    dc, kind = owner
    if (order != [k for k in ('jgmin', 'gshunt', 'source') if k in layers]
            or ('source' in layers) != (kind == 'source') or base[1] != dc.irefnode):
        return None
    lam = layers['source'][1] if 'source' in layers else None
    gs = layers['gshunt'][1] if 'gshunt' in layers else None
    jg = (layers['jgmin'][1], layers['jgmin'][3]) if 'jgmin' in layers else None
    return dc, lam, gs, jg


class _Rec:
    """One DC analysis's struct (`s`, a `dcnewton_t`) and buffers for one
    core, size and reference row, set once: a solve copies its inputs in.
    Its walk fields are `_tran_newton_c._walk_ready`'s (`s`, `ffi`, `walk*`);
    `klu` the StructuralKLU it is bound to (`_klu_bind`), None: the dense
    solve; `limits` whether any element of the circuit limits."""

    __slots__ = ('act', 'buf', 'core', 'ffi', 'iref', 'jrows', 'keep', 'klu', 'klu_keep',
                 'limits', 'n', 'packed', 'pk', 'plan', 's', 'walk', 'walk_arms', 'walk_dicts',
                 'walk_els', 'walk_held', 'walk_kern', 'walk_ok', 'walk_packs', 'walk_stamp',
                 'walk_watch')

    def __init__(self, core, plan, n, iref, drv, core_drv, packed=False):
        from pycircuit.circuit import _tran_core
        from pycircuit.circuit.nrsolver import FLOOR_ROWS
        ffi, _cfn, dgesv, dgemv, ddot = drv
        cffi_core, cfn_core, _dg = core_drv
        m = n - 1
        self.core, self.plan, self.n, self.iref, self.ffi = core, plan, n, iref, ffi
        self.packed, self.pk = packed, ()
        self.s = s = ffi.new('dcnewton_t *')
        E = np.empty
        ## (R and |J| start at zero: only their active entries are written --
        ## on pattern pages, resident only where written: `_tran_core.pattern_zeros`;
        ## THE PACKED LAYOUT reads none of the dense ones)
        Z = _tran_core.pattern_zeros
        if packed:
            Z = D = lambda _size: E(1)
        else:
            D = E
        self.buf = buf = {
            'x': E(m), 'Fr': E(m), 'xn': E(m), 'xd': E(m), 'xf': E(n), 'Ff': E(n), 'xnf': E(n),
            'x0f': E(n), 'Cs': Z(n * n), 'q': E(n), 'R': Z(m * m), 'absJ': Z(m * m),
            'A': D(m * m), 'b': E(m), 'kx': E(m), 'ax': E(m), 'Is': E(m), 'ky': E(m),
            'u': E(n), 'abstol': E(m), 'xtol': E(m),
            ## (the rounding floor's, `nrsolver.FLOOR_ROWS` as the record is made)
            'Tf': D(m * m), 'Ef': E(m * FLOOR_ROWS), 'af': E(m), 'Isf': E(m),
            ## (the sparse verdicts' scratch -- pcx_verdict.h)
            'vt': E(2 * m),
        }
        fb = ffi.from_buffer
        self.keep = [fb('double *', a) for a in buf.values()]
        for k, h in zip(buf, self.keep):
            setattr(s, k, h)
        ipiv = np.empty(max(m, 1), dtype=np.int64)
        self.keep += [ipiv, fb('int64_t *', ipiv)]
        s.ipiv64 = self.keep[-1]
        floor = (np.empty(max(m, 1), dtype=np.int64), np.empty(max(m, 1), dtype=np.int64))
        self.keep += [*floor, fb('int64_t *', floor[0]), fb('int64_t *', floor[1])]
        s.ipivf, s.ffail = self.keep[-2:]
        s.floor_rows = FLOOR_ROWS
        s.dgesv, s.dgemv = ffi.cast('void *', dgesv), ffi.cast('void *', dgemv)
        s.ddot = ffi.cast('void *', ddot)
        s.klu_tsolve = ffi.NULL
        s.core = ffi.cast('void *', int(cffi_core.cast('uintptr_t', core.cs)))
        s.core_fn = ffi.cast('void *', int(cffi_core.cast('uintptr_t', cfn_core)))
        abi = _sources()[4]
        s.pattern = abi['PCX_WANT_PATTERN'] if core.npat else 0
        ## (THE PACKED LAYOUT: the core's packed sets, G in its packed bin)
        s.packed = int(packed)
        if packed:
            s.pattern |= abi['PCX_WANT_SETS']
        s.gmap = s.kdiag = s.jk = ffi.NULL
        s.n, s.iref = n, iref
        s.walk_fn, s.wn = ffi.NULL, 0
        self.walk = self.walk_ok = None
        self.walk_els = self.walk_dicts = self.walk_packs = self.walk_kern = ()
        self.walk_watch = ()
        self.walk_stamp, self.walk_arms, self.walk_held = -1, 0, 0
        self.limits = None
        self.jrows = None
        self.klu = None
        self.klu_keep = ()
        ## (the active entries as the solve binds them -- KLU's pattern, or
        ## every entry for the dense solve -- in `_solver`: m^2 of them made
        ## only where the dense solve is chosen)
        self.act = None
        s.nact = 0

    def _dense_active(self):
        """Every entry active (the dense solve): the reduced row-major
        index and its place in the full n x n G."""
        n, iref, m = self.n, self.iref, self.n - 1
        r, c = np.divmod(np.arange(m * m, dtype=np.int64), m)
        self._set_active(np.arange(m * m, dtype=np.int64),
                         (r + (r >= iref)) * n + (c + (c >= iref)))

    def _set_active(self, act_r, act_f):
        ffi, s = self.ffi, self.s
        act_r = np.ascontiguousarray(act_r, dtype=np.int64)
        act_f = np.ascontiguousarray(act_f, dtype=np.int64)
        kAx = np.empty(max(len(act_r), 1))
        self.act = (act_r, act_f, kAx, ffi.from_buffer('int64_t *', act_r),
                    ffi.from_buffer('int64_t *', act_f), ffi.from_buffer('double *', kAx))
        s.nact = len(act_r)
        s.act_r, s.act_f, s.klu.kAx = self.act[3:]


def _klu_bind(rec, sk):
    """Point the record at `sk`, the DC's StructuralKLU: its pattern the
    active entries (in its CSC order), its arrays, analysis, common and the
    library's entries.  False -- declined -- where it is not analysed yet
    (the Python's first solve does that), the core writes every entry, or
    the core's pattern does not lie within its own (a nonzero off it would
    go to its fallback in the Python)."""
    n, iref, m = rec.n, rec.iref, rec.n - 1
    core = rec.core
    if sk.m != m or not sk._symbolic or not core.npat:
        return False
    r, c = np.divmod(np.asarray(sk.take, dtype=np.int64), m)
    act_f = (r + (r >= iref)) * n + (c + (c >= iref))
    pat = core.arrays['pat']
    pr, pc = np.divmod(pat, n)
    if not np.isin(pat[(pr != iref) & (pc != iref)], act_f).all():
        return False
    if rec.packed and not _packed_maps(rec, sk, act_f):
        return False
    ## (packed, the C reads `gmap` for the core's places: no act_f kept)
    rec._set_active(sk.take, np.zeros(1, np.int64) if rec.packed else act_f)
    ffi, s, lib = rec.ffi, rec.s, sk._lib
    ## (and its held values, pcx_held.h: the same buffers the Python's
    ## `solve` reads and writes)
    keep = (sk, ffi.from_buffer('int *', sk.Ap), ffi.from_buffer('int *', sk.Ai),
            ffi.from_buffer('double *', sk._held_v), ffi.from_buffer('int64_t *', sk._held_st))
    k = s.klu
    k.kAp, k.kAi = keep[1:3]
    s.held_v, s.held_st = keep[3:]
    k.klu_factor, k.klu_refactor, k.klu_solve, k.klu_free = (
        ffi.cast('void *', _cbuild.address(f))
        for f in (lib.klu_factor, lib.klu_refactor, lib.klu_solve, lib.klu_free_numeric))
    s.klu_tsolve = ffi.cast('void *', _cbuild.address(lib.klu_tsolve))
    k.klu_common = ffi.cast('void *', sk._ct.addressof(sk._common))
    k.klu_sym = ffi.cast('void *', sk._symbolic)
    k.resid_tol = float(sk.REFACTOR_RESIDUAL_TOL)
    rec.klu, rec.klu_keep = sk, keep
    return True


def _packed_maps(rec, sk, act_f):
    """THE PACKED LAYOUT's maps (dcnewton.h) for `sk`: each KLU entry's index
    in the core's pattern (-1 off it), each row's diagonal entry in KLU's
    arrays -- False where a row has none (the core's pattern, within KLU's,
    holds the whole diagonal: never so)."""
    m, ffi = rec.n - 1, rec.ffi
    pat = rec.core.arrays['pat']
    g = np.searchsorted(pat, act_f)
    on = g < len(pat)
    on[on] = pat[g[on]] == act_f[on]
    gmap = np.ascontiguousarray(np.where(on, g, -1), dtype=np.int64)
    r, c = np.divmod(np.asarray(sk.take, dtype=np.int64), m)
    kdiag = np.full(m, -1, dtype=np.int64)
    at = np.flatnonzero(r == c)
    kdiag[r[at]] = at
    if (kdiag < 0).any():
        return False
    rec.pk = (gmap, kdiag, ffi.from_buffer('int64_t *', gmap),
              ffi.from_buffer('int64_t *', kdiag), sk)
    rec.s.gmap, rec.s.kdiag = rec.pk[2], rec.pk[3]
    return True


def _junction_entries(rec, ja, jb):
    """Each junction's entries (a, a), (a, b), (b, a), (b, b) in KLU's
    arrays, in the loop's order (THE PACKED LAYOUT) -- or None where one is
    off its pattern (the Python's on the pattern raises there)."""
    m = rec.n - 1
    key = np.asarray(rec.pk[4].take, dtype=np.int64)
    order = np.argsort(key, kind='stable')
    skey = key[order]
    want = np.stack((ja * m + ja, ja * m + jb, jb * m + ja, jb * m + jb), axis=1).reshape(-1)
    i = np.minimum(np.searchsorted(skey, want), max(skey.size - 1, 0))
    if want.size and not (skey[i] == want).all():
        return None
    return np.ascontiguousarray(order[i], dtype=np.int64)


#: the rounding floor's helpers as `nrsolver` defines them, read once
#: genuine (`_settings`)
_FLOOR = []


def _settings(dc, toolkit, limiter, scaler):
    """Why the solve's settings take the Python (a reason), or None: the
    toolkit, the scaler, the limiter, numpy's error state, the rounding
    floor's helpers patched (`nrsolver._step_size`, `_within_rounding`)."""
    from pycircuit.circuit import _limiting, nrsolver
    from pycircuit.circuit import toolkit as _toolkit
    from pycircuit.circuit.scaler import NoneScaler
    floor = (nrsolver._step_size, nrsolver._within_rounding)
    if not _FLOOR:
        if not all(_paths.genuine(f, f.__name__, 'pycircuit.circuit.nrsolver') for f in floor):
            return 'patched'
        _FLOOR.extend(floor)
    if floor != tuple(_FLOOR):
        return 'patched'
    if type(toolkit) is not _toolkit.NumericToolkit:
        return 'toolkit'
    if scaler is not None and type(scaler) is not NoneScaler:
        return 'scaler'
    if getattr(limiter, 'dc_owner', None) is not dc:
        return 'limiter'
    if _limiting.CIRCUIT_LEVEL:
        return 'circuit_limit'
    if np.geterr()['under'] != 'ignore':
        return 'underflow'
    return None


def _record(dc, drv, linsolver, packed=False):
    """The DC analysis's record for its circuit as it stands -- its core
    (built again where the stamp plan was), the solver bound, the walk
    ready -- or why not (a reason).  `packed`: on THE PACKED LAYOUT (a record
    and core of its own kind)."""
    from pycircuit.circuit import _stamp_plan, _tran_core
    cir = dc.cir
    n, iref = cir.n, dc.irefnode
    if not _tran_core.CORE:
        return 'core'
    plan = _stamp_plan._plan_for(cir)
    rec = dc.__dict__.get('_dc_newton_c')
    if (rec is None or rec.plan is not plan or rec.n != n or rec.iref != iref
            or rec.packed != packed):
        core_drv = _tran_core.driver()
        if core_drv is None:
            return 'core'
        try:
            core = _tran_core._Core(cir, plan, core_drv[0], core_drv[2], packed=packed)
        except _tran_core.Unservable:
            return 'unservable'
        rec = dc.__dict__['_dc_newton_c'] = _Rec(core, plan, n, iref, drv, core_drv, packed)
    why = _solver(rec, linsolver)
    if why is None:
        why = _walk(rec, cir)
    return why or rec


def _solver(rec, ls):
    """The record bound to the solve `ls` makes: its dense solve, or its
    analysed StructuralKLU -- or why not."""
    from pycircuit.circuit.linearsolver import AutoSolver, StructuralKLU
    ch = None
    if type(ls) is AutoSolver and rec.packed:
        ch = _packed_choice(ls, rec.n - 1)
        if type(ch) is not StructuralKLU:
            return 'linsolver'
    elif type(ls) is AutoSolver:
        ch = ls._choice or first_choice(ls, rec.n - 1)
    if ch is None:
        return 'linsolver'
    if ch is ls._dense:
        if rec.klu is not None:
            rec.klu, rec.klu_keep = None, ()
            rec.s.klu.klu_sym = rec.ffi.NULL
            rec.act = None
        if rec.act is None:
            rec._dense_active()
        return None
    if type(ch) is not StructuralKLU:
        return 'linsolver'
    if rec.klu is ch and rec.s.klu.klu_sym == rec.ffi.cast('void *', ch._symbolic or 0):
        return None
    return None if _klu_bind(rec, ch) else 'klu'


def _packed_choice(ls, m):
    """THE PACKED LAYOUT's solver: KLU on the structure whatever the size
    (`AutoSolver.choose_packed`, `_select`'s own choice for a PatternMatrix),
    analysed as its first factor analyses it -- or None."""
    ch = ls.choose_packed(m)
    if ch is None or ch is ls._dense:
        return None
    try:
        ch._analysed()
    except np.linalg.LinAlgError:
        return None
    return ch


def structure_chooses(ls, m):
    """Whether the AutoSolver `ls`, not chosen yet, is at its structural
    floor for an m x m matrix or above, where the circuit's structure may
    choose before the first matrix does (`first_choice`)."""
    return (getattr(ls, 'structure', None) is not None
            and m >= getattr(ls, 'structural_min_n', m + 1))


def first_choice(ls, m):
    """The run's first solve, its AutoSolver not chosen yet: the
    StructuralKLU its first matrix would choose where the structure decides
    it, from the structural floor -- analysed as its first factor analyses
    it (`_tran_newton_c._first_choice`) -- or None: below the floor (the size
    alone chooses there, on the first matrix) or where a value could decide.
    pcx_kluseq.h then factors afresh with no numeric yet, as
    `StructuralKLU._solve`'s first solve does."""
    if not structure_chooses(ls, m):
        return None
    from pycircuit.circuit import _tran_newton_c
    return _tran_newton_c._first_choice(ls, m)


def _walk(rec, cir):
    """The limiter as `cir.limit` makes it: the walk, every limiting
    element its C kernel (`_tran_newton_c._walk_ready`); or no walk, where
    nothing limits -- or 'limit'."""
    from pycircuit.circuit import _hdl_climit, _tran_newton_c
    if rec.limits is None:
        rec.limits = any(hasattr(el, 'limit') for el in cir.elements.values())
    if not rec.limits:
        rec.s.walk_fn, rec.s.wn = rec.ffi.NULL, 0
        return None
    if not _tran_newton_c._MOD:
        _tran_newton_c._mods()
    if not (_hdl_climit.WALK and _tran_newton_c._walk_ready(_Of(cir), rec)):
        return 'limit'
    return None


class _Of:
    """What `_tran_newton_c._walk_ready` reads of a transient: its circuit."""

    __slots__ = ('cir',)

    def __init__(self, cir):
        self.cir = cir


def _vector(a, m):
    return type(a) is np.ndarray and a.dtype is _F64 and a.shape == (m,)


def _load(rec, dc, x0, lam, gs, jg, maxiter, reltol, abstol, xtol):
    """The call's inputs into the record -- the temperature (the core's
    readiness), u, the layers, the tolerances, x0, the numeric -- or why
    not (a reason)."""
    s, buf, m = rec.s, rec.buf, rec.n - 1
    T = rec.core.probe(dc.epar)
    if T.__class__ is str:
        return 'ready'
    u = dc.cir.u(0, analysis='dc', epar=dc.epar)
    if not (_vector(u, rec.n) and _vector(x0, m) and _vector(abstol, m) and _vector(xtol, m)):
        return 'x'
    np.copyto(buf['u'], u)
    np.copyto(buf['x'], x0)
    np.copyto(buf['abstol'], abstol)
    np.copyto(buf['xtol'], xtol)
    s.T, s.maxiter, s.reltol = T, int(maxiter), float(reltol)
    s.has_lam, s.lam = (0, 0.0) if lam is None else (1, float(lam))
    s.has_gshunt, s.gshunt = (0, 0.0) if gs is None else (1, float(gs))
    if jg is None:
        s.has_jgmin, s.gjgmin, s.nj = 0, 0.0, 0
    else:
        g, rows = jg
        if rec.jrows is None or rec.jrows[0] is not rows or rec.jrows[5] is not rec.pk:
            ja = np.ascontiguousarray([a for a, _b in rows], dtype=np.int64)
            jb = np.ascontiguousarray([b for _a, b in rows], dtype=np.int64)
            if len(ja) and not ((ja >= 0) & (ja < m) & (jb >= 0) & (jb < m)).all():
                return 'junctions'
            ffi = rec.ffi
            jk = _junction_entries(rec, ja, jb) if rec.packed else np.zeros(1, np.int64)
            if jk is None:
                return 'junctions'
            rec.jrows = (rows, ja, jb, ffi.from_buffer('int64_t *', ja),
                         ffi.from_buffer('int64_t *', jb), rec.pk, jk,
                         ffi.from_buffer('int64_t *', jk))
        s.has_jgmin, s.gjgmin, s.nj = 1, float(g), len(rows)
        s.ja, s.jb, s.jk = rec.jrows[3], rec.jrows[4], rec.jrows[7]
    if rec.klu is not None:
        num = rec.klu._numeric
        s.klu.klu_num = rec.ffi.cast('void *', num) if num else rec.ffi.NULL
    return None


def _write_back(rec):
    """StructuralKLU as its own solve leaves it: the numeric the sequence
    ends with, its counters."""
    sk, s = rec.klu, rec.s
    if sk is None:
        return
    out = int(rec.ffi.cast('uintptr_t', s.klu.klu_out))
    sk._numeric = out or None
    k = s.klu
    sk.refactors += k.refactors
    sk.factors += k.factors
    sk.pivot_fallbacks += k.pivot_fallbacks
    sk.residual_fallbacks += k.residual_fallbacks


def solve(x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names,
          linsolver):
    """`StandardNewton.solve_system` in C for a DC composition: `(x_next,
    iterations)`, or the Python's NoConvergenceError raised (the iteration
    limit) -- or None: not a DC composition (uncounted), declined or handed
    back, nothing of the analysis's written; the Python's solve runs."""
    comp = composition(eval_FJ)
    if comp is None:
        return None    # paths: not a decline (not a DC composition)
    if not ENABLED:
        return _no('dcnewton_c:off')
    dc, lam, gs, jg = comp
    ## (THE PACKED LAYOUT, `_packed`: its circuit's scope the record on the
    ## pattern; another circuit's scope declines)
    packed = _packed.active(dc.cir)
    if not packed and _packed.any_active():
        return _no('dcnewton_c:packed')
    ## (the cheap checks first: a solve they decline never loads the unit)
    why = _settings(dc, toolkit, limiter, scaler)
    if why is not None:
        return _no('dcnewton_c:' + why)
    drv = _driver or driver()
    if not drv:
        return _no('dcnewton_c:driver')
    rec = _record(dc, drv, linsolver, packed)
    if type(rec) is str:
        return _no('dcnewton_c:' + rec)
    why = _load(rec, dc, x0, lam, gs, jg, maxiter, reltol, abstol, xtol)
    if why is not None:
        return _no('dcnewton_c:' + why)
    status = drv[1](rec.s)
    codes, bails = _sources()[2:4]
    if rec.s.dense_conv:
        _PC['dcnewton_c:dense_conv'] += int(rec.s.dense_conv)
    if rec.s.dense_resid:
        _PC['dcnewton_c:dense_resid'] += int(rec.s.dense_resid)
    if status == codes['DCNEWTON_CONVERGED'] or status == codes['DCNEWTON_MAXITER']:
        _write_back(rec)
        _PC['dcnewton_c:served'] += 1
        _PC['dcnewton_c:iterations'] += int(rec.s.iters)
        buf = rec.buf
        if status == codes['DCNEWTON_CONVERGED']:
            return buf['x'].copy(), int(rec.s.iters)
        _PC['dcnewton_c:maxiter'] += 1
        from pycircuit.circuit.nrsolver import _standard_failure
        x_next = buf['xn'].copy()
        ## (the residual's scale `|J| |x_next| + |F|`: the C's last test's `Is`
        ## and F)
        raise _standard_failure(buf['Fr'].copy(), buf['xd'].copy(), x_next, x_next, maxiter,
                                reltol, abstol, xtol, row_names, toolkit,
                                I_scale=buf['Is'] + np.abs(buf['Fr']))
    return _no('dcnewton_c:bail:' + bails.get(status, str(status)))
