"""THE EVALUATE CORE: a transient's residual, Jacobian and companion --
`_CompanionModel._residual_and_jacobian` and `Transient._solve_timestep`'s
`jacobian_only` and `residual_only`, with the four circuit passes they call
(`G`, `C`, `i`, `q`) -- evaluated in ONE C call (`csrc/core.c`), bit for bit; and the circuit's
passes alone for the stage methods (`passes`).  Pinned against its Python
reference in `tests/test_pinned_pairs.py`.

OWNS.  On the transient: `_tran_core` (the record `[plan, core, stamp,
ParameterDict, arms, held]`: the `_Core` for the circuit's stamp plan, None
where the circuit cannot be served, and `core_for`'s stamp) and
`_passes_no` (`passes`'s kept decline: the stamp plan of a circuit the core
cannot serve).  On each HDL element: `_hdl_cp` (its packed parameters, made
by `_Core.probe` where it has none).  On the stamp plan: a pass's
plan where `plan.methods` has none.  `_paths`: `core.<want>:<why>` (`want`
'fj', 'j', 'c', 'f' for `evaluate`, 'p' for `passes`; `<why>` a decline or
`served`), `core.ready:<why>` (`_Core.ready`), `once:core.build:<built |
unservable | nodriver>`, `once:core.fuse:<fused | why not>`.  Switches:
`PYCIRCUIT_TRAN_CORE=0` (`CORE`, read on every call; the C Newton's too),
`PYCIRCUIT_CORE_PASSES=0` (`CORE_PASSES`: `passes` alone, read on every
call), `PYCIRCUIT_CORE_PATTERN=0` (`PATTERN`, read when a core is built),
`PYCIRCUIT_SKIP_UNREAD_J=0` (`SKIP_UNREAD_J`, read by the transient's
stepping loop); `STATUS` says why the core is off where it is.

ONE CALL does the passes the caller asked for -- the batches through the
kernels' own pointers, the constant vector groups through numpy's OWN
`cblas_dgemv` (the address of the routine numpy loaded, so the same kernel
runs as `np.dot` and the plan's `matmul`), the templates, the bincounts
from +0.0 in array order -- then the companion in the integrator's own
operation order, `f = (i + iq) + u` and `J = G + Geq`.  The reference-row
reduction stays Python's.

FOUR SITES, ONE ENTRY (`evaluate`, `want` in 'fj' | 'j' | 'c' | 'f'):
`_residual_and_jacobian` (G, C, i, q; `(f, J)`); `jacobian_only` at the
converged point (G, C, q; `(None, J)`), or its state alone where nothing
reads its `J` (C, q; `(None, None)`); the chord's `residual_only` (i, q;
the conductance from the seed's held `C`; `f`).  Each leaves the state the
Python path leaves: `_C_cache`, `_q_cache`, the memo record `_C_at_state`
makes, and `get_diff`'s `_iq`, `_Geq`, `_Cmat`, `_companion_coeffs`,
`_effective_method`, `active_integrator` -- and `_C_pattern` where it wrote
`C` on the pattern.  Python still decides the active integrator
(`check_order_drop`) and computes its scalars (`bdf2_alphas`, `theta_at`),
and still asks `_C_lookup` first: on a hit the C pass is skipped, as
`_C_at_state` skips it.

THE FORMULA ENUM, each as `_lte_kernels` writes it: gear2 `(a0*q + a1*q1)
+ a2*q2`, `a0*C`; euler `(q - q1)/h`, `C/h`; trap `2.0*(q - q1)/h - iq1`,
`2.0*C/h`; theta `th = theta*h`, `(q - q1)/th - ((1.0 - theta)/theta)*iq1`,
`C/th`; and `PCX_FORMULA_PASSES` (4), the passes alone (`passes`).
`evaluate` does not serve the stage methods (Radau, the DIRKs: their passes
come from `passes`), PCNR, DC or any toolkit but the numeric one.

A DECLINE BEFORE ANY STATE IS TOUCHED (None; the Python path runs): the
switch off, a toolkit that is not numeric, a state that is not one finite
float64 vector of the circuit's size, a circuit with an element neither
constant nor C-bound (a hand-written limiter, a Python hdl class, a nested
circuit: the core does not exist for it until its plan rebuilds), a batch
not ready (a detached class, a patched method, an instance shadow of any
of the four methods, a pack that failed, a temperature that is not one
number), a history row of the wrong shape, an integrator outside the enum.
The core's tables are keyed by the stamp plan (rebuilt when it is), the
pack mirrors are the core's own (`_hdl_batch`'s rule: mirror the tuple the
element holds, never repack one), the readiness check runs once per call
over the unique elements rather than once per pass.

THE PATTERN.  The core's pattern is every entry a matrix pass can write and
the diagonal; where it is at most `PATTERN_MAX_FILL` of n x n (and
`PATTERN` is on) a caller whose matrices start at zero asks the core to
write only its entries (`PCX_WANT_PATTERN`; `_Core.npat` nonzero), and
every entry otherwise.

THE THREADS.  A call's fused-group elements may be made on a pool of
threads (core.c; `THREADS`, `PYCIRCUIT_CORE_THREADS`): each element's
kernel is a function of its gathered inputs and its parameter row alone,
writing its own rows of the staging, and the bins are summed after by the
caller in their order -- every byte, and every floating-point flag, as one
thread leaves them.

THE FUSED KERNELS.  A class's four passes print four kernels that each
recompute their whole chain; its fused kernel (`_hdl_cbackend.fuse_csrc`:
the printed statements unioned by name, each run under the bits of the
passes that print it) is called ONCE a call for every pass the call wants,
before the passes, into a staging buffer -- where the class's batches cover
two or more passes over the same elements on the same nodes
(`_fuse_groups`) -- and each batch then scatters its outputs from there in
its own order: the same values in the same places, the passes' bins as
without it.  The readiness check takes the class's fused kernel only while
it is bound and was printed from the passes' bound kernels (and the switch,
`_hdl_cbackend.FUSE`, is on); else the batches call their own.

READY AS STAMPED.  The readiness check (`probe`: every batch's kernel the
class's bound one, its pass the generated one, every element without a
shadow, every pack mirrored) and the lookup (`core_for`: the circuit's
plan) read the same dicts every call, and their answer does not change
between steps: a check that passed watches what it read (`_watch`: one
counter any change to a watched dict moves) and passes again unchecked
while the counter stands -- what no dict watcher sees compared on every
call: the elements' `__dict__` objects, the passes' code, `ParameterDict`'s
epoch, the circuit's size.

History: `doc/transient_history.md`, `_tran_core.py`.
"""
import operator
import os

import numpy as np

from pycircuit.circuit import _cbuild, _lte_kernels, _paths, _watch
from pycircuit.circuit import toolkit as _toolkit

_PC = _paths.COUNTS
_GETDICT = operator.attrgetter('__dict__')
#: float64 (an exact array's dtype is this object: `is`, then `==`)
_F64 = np.dtype(np.float64)
_no = _paths.no
#: the counters' keys (`_paths`), by `want` and reason: each decline by its
#: reason, each served call (`core.<want>:<why>`)
_K = {w: {r: f'core.{w}:{r}' for r in (
    'off', 'toolkit', 'shadow_cir', 'shadow_tr', 'unservable', 'x', 'formula',
    'history', 'cmat', 'clookup', 'ready', 'u', 'served')} for w in ('fj', 'j', 'f', 'c', 'p')}

CORE = os.environ.get('PYCIRCUIT_TRAN_CORE', '1') != '0'
#: where nothing reads a step's Jacobian (a fixed multistep run:
#: `transient._SteppingLoop.j_unread`), the converged point is evaluated
#: without G -- `want` 'c' here, the C Newton's fold likewise; env
#: `PYCIRCUIT_SKIP_UNREAD_J=0` turns it off
SKIP_UNREAD_J = os.environ.get('PYCIRCUIT_SKIP_UNREAD_J', '1') != '0'
STATUS = 'not loaded'

CORE_C = _cbuild.csrc('core.c')
CORE_CDEF = _cbuild.csrc('core.h')
#: what an element may be for the core's stack buffers (core.h)
_CAP = _cbuild.constants(CORE_CDEF)
_MAX_IN, _MAX_OUT = _CAP['CORE_MAX_INPUTS'], _CAP['CORE_MAX_OUTPUTS']
#: a fused group's row in `ztab` (core.h's `core_ztab`)
_ZT_WIDTH = _CAP['ZT_WIDTH']

#: the circuit passes, in the core's order (its tables are laid out pass
#: by pass; a fused group's pass mask is `1 << PASSES.index(m)`)
PASSES = ('G', 'C', 'i', 'q')
#: the core's call codes, named in pcx_abi.h
_ABI = _cbuild.constants(_cbuild.csrc('pcx_abi.h'))
#: THE CIRCUIT'S PATTERN: a caller whose matrices are zero off it asks the
#: core to write only the pattern's entries (`PCX_WANT_PATTERN`) -- for a
#: core whose pattern is at most `PATTERN_MAX_FILL` of n x n (else every
#: entry); env `PYCIRCUIT_CORE_PATTERN=0` turns it off (read when a core is
#: built)
PATTERN = os.environ.get('PYCIRCUIT_CORE_PATTERN', '1') != '0'
PATTERN_MAX_FILL = 0.25
def _default_threads():
    """Four, or half the CPUs this process may run on where that is fewer:
    jge's Radau transient (348 MOSFETs a call) 1.44x on two threads, 1.88x on
    four, 2.22x on eight; mux8 (64) 1.58x on four, 1.78x on eight -- the
    workers spinning between calls, half the box held for most of the gain."""
    try:
        ncpu = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        ncpu = os.cpu_count() or 1
    return max(1, min(4, ncpu // 2))


#: THE DEVICE KERNELS ON THREADS (core.c's pool): a call's fused groups on
#: this many threads, the caller one of them, where its elements reach the
#: pool's minimum -- every byte as one thread makes them
#: (`PYCIRCUIT_CORE_THREADS`; the test suite's 1, its workers filling the box)
THREADS = max(1, int(os.environ.get('PYCIRCUIT_CORE_THREADS', '') or _default_threads()))
_WANT_PATTERN = _ABI['PCX_WANT_PATTERN']
#: the integrators the core's companion enum covers, by class name
FORMULA = {'Gear2Integrator': _ABI['PCX_FORMULA_GEAR2'],
           'EulerIntegrator': _ABI['PCX_FORMULA_EULER'],
           'TrapezoidalIntegrator': _ABI['PCX_FORMULA_TRAPEZOIDAL'],
           'ThetaIntegrator': _ABI['PCX_FORMULA_THETA']}
_G, _C, _I = _ABI['PCX_PASS_G'], _ABI['PCX_PASS_C'], _ABI['PCX_PASS_I']
#: the companions by name: Gear-2's (two charges back), theta's (its theta), and
#: the two whose companion carries the current (`iq`)
_GEAR2, _THETA = _ABI['PCX_FORMULA_GEAR2'], _ABI['PCX_FORMULA_THETA']
_IQ_FORMS = (_ABI['PCX_FORMULA_TRAPEZOIDAL'], _THETA)
#: each `want`'s pass bits (q is always made)
WANT = {'fj': _G | _C | _I, 'j': _G | _C, 'f': _I, 'c': _C}

_driver = None


def _numpy_dgemv():
    """The address of numpy's own `cblas_dgemv` (`_cbuild.openblas`), or
    None where the wheel has no such library."""
    try:
        lib = _cbuild.openblas(np)
        return None if lib is None else _cbuild.address(lib.scipy_cblas_dgemv64_)
    except (OSError, AttributeError):
        return None


def driver():
    """`(ffi, cfn, dgemv address)` once loaded, or None (`STATUS` says why)."""
    global _driver, STATUS
    if _driver is None:
        from pycircuit.circuit import _hdl_cbackend as cb
        dgemv = _numpy_dgemv()
        if dgemv is None:
            _driver = False
            STATUS = "off (numpy's own cblas_dgemv not found)"
        else:
            try:
                ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(CORE_C, CORE_CDEF, 'pcx_core')
            except (cb.CompileError, OSError) as e:
                _driver = False
                STATUS = f'off ({e})'
            else:
                _driver = (ffi, cfn, dgemv)
                STATUS = 'c'
    return _driver or None


class Unservable(Exception):
    """The circuit has an element the core cannot evaluate."""


def _fuse_groups(batches, b_pass, cb):
    """The classes the core evaluates with their fused kernel:
    `(groups, b_fz, b_fidx)` -- each group `(info, kernel, elements, pass
    bits, node maps)`, per batch its group (-1: none) and each of its
    elements' place in the group (None: the group's own order).  A class
    is fused where its fused kernel is bound and its batches cover two or
    more of its passes, the same elements in each, once each, on the same
    nodes, with the kernel's widths; each refusal counted once a build
    (`once:core.fuse:<why>`).  (Every pass lists a class's elements in one
    order, the plan's: the node maps are then compared whole; a batch in
    another order is mapped through the group's element index first.)"""
    by_cls = {}
    for b, bt in enumerate(batches):
        by_cls.setdefault(bt.cls, []).append(b)
    groups, b_fz, b_fidx = [], [-1] * len(batches), [None] * len(batches)
    for bs in by_cls.values():
        first = batches[bs[0]]
        info = first.info
        fk = info.get('_c_fused') if cb.FUSE else None
        if fk is None and cb.FUSE and info.get('_c_bound') and info.get('_c_fused_status') == 'off':
            ## (bound while the switch was off: built now)
            fk = cb.bind_fused(info, cb.fuse_source(info))
        passes = [PASSES[b_pass[b]] for b in bs]
        why = None
        if fk is None:
            why = 'nokernel'
        elif len(set(passes)) != len(passes):
            why = 'passes'          # (a class split into two batches in one pass)
        elif len(passes) < 2 or not set(passes) <= set(fk.parts):
            why = 'passes'
        else:
            ti = -1 if fk.layout[1] is None else fk.layout[1]
            for b, m in zip(bs, passes):
                bt = batches[b]
                if bt.k != fk.nx or bt.so != fk.sizes[m] or bt.ti != ti:
                    why = 'widths'
                    break
        els, places = first.els, {}
        if why is None:
            nm0 = first.NM.tobytes()
            index = None
            for b in bs:
                bt = batches[b]
                if bt.els is els or (len(bt.els) == len(els)
                                     and all(map(operator.is_, bt.els, els))):
                    ## (the plan's one order: the node maps compared whole,
                    ## both int64 and C-ordered as `Batch` makes them)
                    if bt is not first and not (bt.NM.shape == first.NM.shape
                                                and bt.NM.tobytes() == nm0):
                        why = 'nodes'
                        break
                    places[b] = None
                    continue
                if index is None:
                    index = {id(el): i for i, el in enumerate(els)}
                rows = [index.get(id(el), -1) for el in bt.els]
                if (len(index) != len(els) or len(rows) != len(els) or -1 in rows
                        or len(set(rows)) != len(rows)):
                    why = 'elements'
                    break
                rows = np.asarray(rows, dtype=np.int64)
                if not np.array_equal(bt.NM, first.NM[rows]):
                    why = 'nodes'
                    break
                places[b] = rows
        if why is not None:
            _PC['once:core.fuse:' + why] += 1
            continue
        z = len(groups)
        mask = sum(1 << PASSES.index(m) for m in passes)
        groups.append((info, fk, els, mask, first.NM))
        for b in bs:
            b_fz[b] = z
            b_fidx[b] = places[b]
        _PC['once:core.fuse:fused'] += 1
    return groups, b_fz, b_fidx


def _i64(a):
    return np.ascontiguousarray(np.asarray(a, dtype=np.int64))


def _cat64(parts):
    return (np.ascontiguousarray(np.concatenate(parts).astype(np.int64))
            if parts else np.zeros(1, np.int64))


def _lone_batches(cir, m, mp, cb):
    """Pass `m`'s lone C-bound elements (`mp.calls`), each a batch of one --
    `split` batches two or more, the core serves a lone one too -- or
    `Unservable` where one is not a chained, C-bound element the core can
    take whole.  A class with DC pins is taken: it evaluates
    its pinned functions only in a DC analysis, which the core declines
    (`_Core.probe`); `split` keeps such a class out of its batches, which
    serve the DC analyses too."""
    from pycircuit.circuit import _hdl_batch
    mat = m in ('G', 'C')
    extra = []
    for c_ in mp.calls:
        _inst, el, nm, a, b = c_
        cls = type(el)
        info = getattr(cls, '_hdl_info', None)
        fn = (info or {}).get('funcs', {}).get(m)
        kern = fn.__dict__.get('_hdl_c') if fn is not None else None
        k = len(nm)
        if not (a >= 0 and info is not None and info.get('chained')
                and el.toolkit is cir.toolkit
                and m not in el.__dict__ and _hdl_batch.is_generated(cls, m)
                and isinstance(kern, cb.CKernel) and kern.nx == k
                and kern.shape == ((k, k) if mat else (k,))
                and b - a == (k * k if mat else k)):
            raise Unservable(f'{_inst} in the {m} pass')
        extra.append(_hdl_batch.Batch(cls, m, [c_], kern))
    return extra


class _Flat:
    """The four plans flattened, as `_Core` lays them out for the struct:
    per pass its template, its slots (`flat`) and its bins; per batch its
    pass, node map, destinations and first place in the pointer table; per
    constant group (a vector pass's) its pass, matrix, inputs, selection
    and destinations; per fused group (`add_fused`) its `core_ztab` row.
    Each list's running offset beside it."""

    def __init__(self):
        self.p_L, self.p_toff, self.tmpl, self.flat, self.p_nbins = [], [], [], [], []
        self.batches, self.b_pass, self.NM, self.b_nmoff = [], [], [], []
        self.dst, self.b_dstoff, self.b_proff = [], [], []
        self.g_pass, self.S, self.g_soff, self.I, self.g_ioff = [], [], [], [], []
        self.sel, self.gdst, self.g_seloff = [], [], []
        self.off = self.nmoff = self.dstoff = self.pr_len = 0
        self.soff = self.ioff = self.seloff = 0
        self.maxy = 1

    def add_pass(self, mi, m, mp, n, extra):
        """Pass `m` (`PASSES[mi]`), its plan `mp`: its template and slots,
        its constant groups (a vector pass's), its batches and then `extra`
        (its lone elements, `_lone_batches`)."""
        if m in ('G', 'C'):
            L = mp.template.size
            self.tmpl.append(mp.template)
            self.flat.append(mp.flat)
            self.p_nbins.append(n * n)
        else:
            L = mp.length
            self.tmpl.append(np.zeros(L))
            self.flat.append(mp.idx)
            self.p_nbins.append(n)
            for Sg, Ig, selg, dstg, ok in mp.groups:
                if not ok:
                    raise Unservable(f'a constant group of the {m} pass is not matmul-exact')
                self._add_group(mi, Sg, Ig, selg, dstg)
        self.p_L.append(L)
        self.p_toff.append(self.off)
        self.off += L
        for bt in list(mp.batches) + extra:
            self._add_batch(mi, bt)

    def _add_group(self, mi, Sg, Ig, selg, dstg):
        E, k = Ig.shape
        self.g_pass.append(mi)
        self.S.append(Sg.reshape(-1))
        self.g_soff.append((self.soff, E, k))
        self.soff += Sg.size
        self.I.append(Ig.reshape(-1))
        self.g_ioff.append(self.ioff)
        self.ioff += Ig.size
        self.sel.append(selg)
        self.gdst.append(dstg)
        self.g_seloff.append((self.seloff, selg.size))
        self.seloff += selg.size
        self.maxy = max(self.maxy, E * k)

    def _add_batch(self, mi, bt):
        self.batches.append(bt)
        self.b_pass.append(mi)
        self.NM.append(bt.NM.reshape(-1))
        self.b_nmoff.append(self.nmoff)
        self.nmoff += bt.NM.size
        self.dst.append(bt.dst)
        self.b_dstoff.append(self.dstoff)
        self.dstoff += bt.dst.size
        self.b_proff.append(self.pr_len)
        self.pr_len += bt.n

    def add_fused(self, zgroups, b_fz, b_fidx):
        """The fused groups (`_fuse_groups`'s): each one's `core_ztab` row,
        its node map and its elements' places in the pointer table (after
        every batch's), its staging; per batch its group and its elements'
        order in it (`bftab`, `fidx`).  Returns each group's first place in
        the pointer table."""
        z_proff, self.ztab, self.fidx, self.bftab = [], [], [], []
        self.st_len = fo = 0
        for _info, fk, els, mask, znm in zgroups:
            ## (a row in core.h's `core_ztab` order: ZT_NE, ZT_K, ZT_TI,
            ## ZT_PASSES, ZT_NM, ZT_PR, then each pass's _LEN and _AT)
            ne = len(els)
            row = (ne, fk.nx, -1 if fk.layout[1] is None else fk.layout[1], mask, self.nmoff,
                   self.pr_len)
            z_proff.append(self.pr_len)
            self.pr_len += ne
            self.NM.append(znm.reshape(-1))
            self.nmoff += znm.size
            for m in PASSES:
                sz = fk.sizes.get(m, 0)
                row += (sz, self.st_len)
                self.st_len += sz * ne
            assert len(row) == _ZT_WIDTH
            self.ztab += row
        if zgroups:
            for b in range(len(self.batches)):
                rows = b_fidx[b]
                if rows is None:
                    ## (the group's own order: no table)
                    self.bftab += (b_fz[b], -1)
                else:
                    self.bftab += (b_fz[b], fo)
                    self.fidx.append(rows)
                    fo += len(rows)
        return z_proff

    def arrays(self, n, fused):
        """The struct's arrays by their `hdl_core_t` names, the pattern among
        them (and, `fused`, the fused groups' tables): `(arrays, npat)`."""
        g_soff, batches = self.g_soff, self.batches
        arrays = {
            'p_L': _i64(self.p_L), 'p_toff': _i64(self.p_toff), 'p_nbins': _i64(self.p_nbins),
            'tmplall': np.ascontiguousarray(np.concatenate(self.tmpl)), 'flatall': _cat64(self.flat),
            'b_pass': _i64(self.b_pass), 'b_ne': _i64([bt.n for bt in batches]),
            'b_k': _i64([bt.k for bt in batches]), 'b_so': _i64([bt.so for bt in batches]),
            'b_ti': _i64([bt.ti for bt in batches]), 'NMall': _cat64(self.NM),
            'b_nmoff': _i64(self.b_nmoff), 'b_proff': _i64(self.b_proff),
            'bdstall': _cat64(self.dst), 'b_dstoff': _i64(self.b_dstoff),
            'g_pass': _i64(self.g_pass), 'g_E': _i64([e for _o, e, _k in g_soff]),
            'g_k': _i64([k for _o, _e, k in g_soff]),
            'Sall': np.ascontiguousarray(np.concatenate(self.S)) if self.S else np.zeros(1),
            'g_soff': _i64([o for o, _e, _k in g_soff]), 'Iall': _cat64(self.I),
            'g_ioff': _i64(self.g_ioff), 'g_nsel': _i64([sz for _o, sz in self.g_seloff]),
            'selall': _cat64(self.sel), 'gdstall': _cat64(self.gdst),
            'g_seloff': _i64([o for o, _sz in self.g_seloff]),
            'work': np.empty(max(int(max(self.p_L)), 1)), 'ywork': np.empty(self.maxy),
            'binG': np.zeros(n * n), 'bini': np.empty(n),
        }
        ## THE PATTERN: every entry a matrix pass can write (G's and C's
        ## slots) and the diagonal -- `binG` above starts at zero, and a
        ## caller's matrices are zero off it where it asks for it
        mats = [f for f, nb in zip(self.flat, self.p_nbins) if nb == n * n]
        pat = np.unique(np.concatenate(mats + [np.arange(n, dtype=np.int64) * (n + 1)]))
        npat = len(pat) if PATTERN and len(pat) <= PATTERN_MAX_FILL * n * n else 0
        arrays['pat'] = np.ascontiguousarray(pat, dtype=np.int64)
        if fused:
            ## (the fused groups' tables: none at all where nothing fuses)
            arrays.update(ztab=np.array(self.ztab, dtype=np.int64),
                          zst=np.empty(max(self.st_len, 1)),
                          bftab=np.array(self.bftab, dtype=np.int64), fidxall=_cat64(self.fidx))
        return arrays, npat


class _Core:
    """One circuit's four plans flattened into the C struct (`cs`), the
    batches' pointer tables and the core's own pack mirrors.  Building one
    fills a pass's plan the stamp plan lacks (`plan.methods`) and raises
    `Unservable` where an element is neither constant nor C-bound, a
    constant group is not matmul-exact, or an element is wider than the
    core holds.  `npat` is the pattern's size, 0 where the core writes
    every entry; `arrays['pat']` the pattern."""

    def __init__(self, cir, plan, ffi, dgemv):
        from pycircuit.circuit import _hdl_cbackend as cb
        from pycircuit.circuit import _stamp_plan
        self.plan = plan
        self.n = n = cir.n
        for m in PASSES:
            if plan.methods.get(m) is None:
                plan.methods[m] = (_stamp_plan._MatrixPlan(cir, m) if m in ('G', 'C')
                                   else _stamp_plan._VectorPlan(cir, m))
        fl = _Flat()
        for mi, m in enumerate(PASSES):
            mp = plan.methods[m]
            fl.add_pass(mi, m, mp, n, _lone_batches(cir, m, mp, cb))
        batches = self.batches = fl.batches
        self.nb, self.ng = len(batches), len(fl.g_pass)
        ## a class with DC pins among them (`_lone_batches`): `probe` declines
        ## a DC analysis's epar
        self.pins = any(b.info['state_meta']['dc_pins'] for b in batches)
        ## THE FUSED KERNELS: a class whose fused kernel is bound
        ## (`_hdl_cbackend.bind_fused`) and whose batches cover two or more
        ## passes over the same elements -- each once a pass, on the same
        ## nodes -- is evaluated once a call into a staging buffer; its
        ## batches scatter from there
        zgroups, b_fz, b_fidx = _fuse_groups(batches, fl.b_pass, cb)
        self._cb = cb
        z_proff = fl.add_fused(zgroups, b_fz, b_fidx)
        self.zgroups = zgroups
        self.nz = len(zgroups)
        self.b_fz = b_fz
        ## (`probe`'s stamp: none yet)
        self.stamp, self.uniq_els, self.uniq_dicts, self.codes, self.arms = -1, (), (), (), 0
        ## (the checks the stamp passed since it was armed: `_watch.counted`)
        self.held = 0
        self.z_fn = np.zeros(max(self.nz, 1), dtype=np.uintp)
        self.zk = [None] * self.nz
        self.b_fn = np.zeros(max(self.nb, 1), dtype=np.uintp)
        self.PRall = np.zeros(max(fl.pr_len, 1), dtype=np.uintp)
        self.b_proff = fl.b_proff
        arrays, self.npat = fl.arrays(n, bool(zgroups))
        if any(bt.k > _MAX_IN or bt.so > _MAX_OUT for bt in batches) \
                or any(k > _MAX_IN for _o, _e, k in fl.g_soff):
            raise Unservable('an element wider than the core holds')
        self.arrays = arrays
        self._struct(ffi, dgemv)
        self._index_elements(z_proff)
        self.dptr = ffi.typeof('double *')

    def _struct(self, ffi, dgemv):
        """The C struct (`cs`): the sizes, every array by its field's name,
        the pointer tables, numpy's dgemv -- each buffer's handle kept."""
        cs = ffi.new('hdl_core_t *')
        cs.n = self.n
        cs.nb, cs.ng = self.nb, self.ng
        cs.npat = self.npat
        self.handles = []
        for name, arr in self.arrays.items():
            h = ffi.from_buffer('double *' if arr.dtype == np.float64 else 'int64_t *', arr)
            self.handles.append(h)
            setattr(cs, name, h)
        self.handles.append(ffi.from_buffer('void **', self.b_fn))
        cs.b_fn = self.handles[-1]
        self.handles.append(ffi.from_buffer('double **', self.PRall))
        cs.PRall = self.handles[-1]
        cs.nz = self.nz
        cs.nthreads = THREADS
        self.handles.append(ffi.from_buffer('void **', self.z_fn))
        cs.z_fn = self.handles[-1]
        cs.dgemv = ffi.cast('void *', dgemv)
        self.cs = cs
        self.ffi = ffi

    def _index_elements(self, z_proff):
        """The unique elements across the batches, each with its places in
        the pointer table (a fused group's too: `z_proff`, its first): one
        shadow-and-pack check per call (`uniq`, `mirror`)."""
        pos = {}
        for b, bt in enumerate(self.batches):
            for e, el in enumerate(bt.els):
                pos.setdefault(id(el), (el, []))[1].append(self.b_proff[b] + e)
        for z, g in enumerate(self.zgroups):
            for e, el in enumerate(g[2]):
                pos[id(el)][1].append(z_proff[z] + e)
        self.uniq = [(el, tuple(ps)) for el, ps in pos.values()]
        self.mirror = [None] * len(self.uniq)

    def ready(self, epar):
        """`Batch.run`'s checks, once: per batch the class (bound, the
        kernel the one cast, the method generated), per unique element the
        shadows of all four methods and the pack mirror.  The temperature,
        or None (the Python path), counted."""
        T = self.probe(epar)
        if T.__class__ is str:
            return _no('core.ready:' + T)
        return T

    def probe(self, epar):
        """`ready` without the count: the temperature, or why not (a str).
        For a caller that falls back to a path which asks `ready` itself
        (`_tran_newton_c`): the decline is counted once, there.  As stamped,
        the stamp's checks alone; else the batches' kernels
        (`_check_batches`), the fused kernels (`_check_fused`), the elements
        (`_check_elements`), and a new stamp.  Declined first: a DC
        analysis's epar where the core holds a class with DC pins."""
        T = getattr(epar, 'T', 300.0)
        if type(T) is not float:
            if type(T) is not int and np.ndim(T) != 0:
                return 'T'
            T = float(T)
        ## (the pinned functions are a DC analysis's, and the core evaluates
        ## the transient's: `_lone_batches`)
        if self.pins and getattr(epar, 'analysis_kind', None) == 'dc':
            return 'pins_dc'
        ## AS STAMPED: nothing the check below read has changed since it
        ## passed -- its dicts are watched (`_watch`); what no watcher sees,
        ## compared here: the elements' `__dict__` objects and the passes'
        ## code
        ep = _watch.EPOCH
        if (ep is not None and self.stamp == ep.value
                and tuple(map(_GETDICT, self.uniq_els)) == self.uniq_dicts
                and all(map(_same_code, self.codes))):
            self.held += 1
            return T
        before = _watch.now()
        from pycircuit.circuit import _hdl_batch
        kern_of = self._check_batches(_hdl_batch)
        if kern_of.__class__ is str:
            return kern_of
        if self.nz:
            self._check_fused()
        why = self._check_elements(kern_of)
        if why is not None:
            return why
        self._stamp(before, _hdl_batch)
        return T

    def _check_batches(self, _hdl_batch):
        """Each batch's class bound, its kernel the one taken (`b_fn`), its
        method generated: per class one of its kernels (`_check_elements`
        packs by it), or why not."""
        ffi, b_fn = self.ffi, self.b_fn
        kern_of = {}
        for b, bt in enumerate(self.batches):
            info = bt.info
            if not info.get('_c_bound'):
                return 'unbound'
            kern = info['funcs'][bt.m].__dict__.get('_hdl_c')
            if kern is not bt.kern:
                if not bt._take(kern):
                    return 'kernel'
                b_fn[b] = int(ffi.cast('uintptr_t', kern.cfn))
                if self.b_fz[b] >= 0:
                    ## (its group's fused kernel is checked again below)
                    self.zk[self.b_fz[b]] = None
            elif not b_fn[b]:
                b_fn[b] = int(ffi.cast('uintptr_t', kern.cfn))
            if not _hdl_batch.is_generated(bt.cls, bt.m):
                return 'generated'
            ## (a class's kernels share its layout: any of them packs it)
            kern_of.setdefault(bt.cls, kern)
        return kern_of

    def _check_fused(self):
        """A FUSED KERNEL SERVES while it is the class's bound one and was
        printed from the passes' bound kernels; else the batches call their
        own (the switch `_hdl_cbackend.FUSE` too): `z_fn` 0."""
        ffi, zk, z_fn, on = self.ffi, self.zk, self.z_fn, self._cb.FUSE
        for z, g in enumerate(self.zgroups):
            info = g[0]
            fk = info.get('_c_fused') if on else None
            if fk is not zk[z] or (fk is None and z_fn[z]):
                zk[z] = fk
                funcs, fk0 = info['funcs'], self.zgroups[z][1]
                ok = (fk is not None and fk.sizes == fk0.sizes and fk.nx == fk0.nx
                      and fk.layout == fk0.layout and all(
                          k is funcs[m].__dict__.get('_hdl_c')
                          for m, k in fk.parts.items()))
                z_fn[z] = int(ffi.cast('uintptr_t', fk.cfn)) if ok else 0

    def _check_elements(self, kern_of):
        """Each unique element unshadowed and packed (`_hdl_cp`, packed by
        its class's kernel where it is not yet), its pack's address in the
        pointer table wherever the pack moved: None, or why not."""
        mirror, PR = self.mirror, self.PRall
        for ui, (el, ps) in enumerate(self.uniq):
            d = el.__dict__
            if 'G' in d or 'C' in d or 'i' in d or 'q' in d:
                return 'shadow'
            cp = d.get('_hdl_cp')
            if cp is None:
                ## (packed by ITS class's kernel, as its batch packs it: the
                ## array's size and the temperature's slot are the class's)
                k = kern_of.get(el.__class__)
                if k is None:
                    return 'pack'
                try:
                    cp = k.pack(el)
                except (TypeError, ValueError):
                    cp = False
                d['_hdl_cp'] = cp
            if cp is False:
                return 'pack'
            if cp is not mirror[ui]:
                mirror[ui] = cp
                a = cp[0].ctypes.data
                for p_ in ps:
                    PR[p_] = a
        return None

    def _stamp(self, before, _hdl_batch):
        """Watch what the check read and stamp it (`_watch.arm`): the
        classes' info dicts and pass lists, the passes' function dicts (their
        bound kernels), the classes, the elements' dicts, the backend's and
        the batches' modules (their switches and drivers)."""
        infos = {id(bt.info): bt.info for bt in self.batches}
        dicts = list(infos.values())
        dicts += [info['funcs'] for info in infos.values()]
        dicts += [info['funcs'][bt.m].__dict__ for bt in self.batches
                  for info in (bt.info,)]
        self.uniq_els = tuple(el for el, _ps in self.uniq)
        self.uniq_dicts = tuple(map(_GETDICT, self.uniq_els))
        dicts += self.uniq_dicts
        dicts += [_hdl_batch.__dict__, self._cb.__dict__]
        self.codes = tuple({(bt.cls, bt.m): (bt.cls, bt.m, getattr(bt.cls, bt.m).__code__)
                            for bt in self.batches}.values())
        self.arms, self.held = _watch.counted(self.arms, self.held), 0
        self.stamp = _watch.arm(dicts, before) if self.arms <= _watch.MAX_ARMS else -1


def _same_code(c):
    """A class's pass still the code it was at the stamp (a function's
    `__code__` is no dict)."""
    return getattr(c[0], c[1]).__code__ is c[2]


def core_for(tr):
    """The transient's core for its circuit's current plan, or None where
    the circuit cannot be served (remembered until the plan rebuilds).
    As stamped (`_watch`): the circuit's dict, its elements and its node map
    unchanged since, the circuit's plan the one recorded, `ParameterDict`'s
    epoch the plan's and the circuit's size the plan's -- the record stands
    without asking `_plan_for`."""
    cir = tr.cir
    rec = tr.__dict__.get('_tran_core')
    ep = _watch.EPOCH
    if (rec is not None and ep is not None and rec[2] == ep.value
            and cir.__dict__.get('_stamp_plan') is rec[0] and rec[3]._epoch is rec[0].epoch
            and rec[0].n == len(cir.nodes) + len(cir.branches)):
        ## (the record `[plan, core, stamp, ParameterDict, arms, held]`: its
        ## checks passed since armed counted, `_watch.counted`)
        rec[5] += 1
        return rec[1]
    before = _watch.now()
    from pycircuit.circuit import _stamp_plan
    from pycircuit.utilities.param import ParameterDict
    plan = _stamp_plan._plan_for(cir)
    if rec is not None and rec[0] is plan:
        arms = _watch.counted(rec[4], rec[5])
        tr.__dict__['_tran_core'] = [plan, rec[1], _plan_stamp(cir, before, arms),
                                     ParameterDict, arms, 0]
        return rec[1]
    drv = driver()
    core = None
    if drv is not None:
        ffi, _cfn, dgemv = drv
        try:
            core = _Core(cir, plan, ffi, dgemv)
        except Unservable:
            core = None
    _PC['once:core.build:' + ('nodriver' if drv is None else
                              'unservable' if core is None else 'built')] += 1
    tr.__dict__['_tran_core'] = [plan, core, _plan_stamp(cir, before, 1), ParameterDict, 1, 0]
    return core


def _plan_stamp(cir, before, arms):
    """`_plan_for`'s dict reads, watched: the circuit's dict (its plan, its
    element and node-map attributes), the elements and the node map
    themselves.  (`ParameterDict`'s epoch, a class attribute, is compared
    on every call.)"""
    if arms > _watch.MAX_ARMS:
        return -1
    nodemap = cir.__dict__.get('elementnodemap')
    dicts = [cir.__dict__, cir.elements]
    if isinstance(nodemap, dict):
        dicts.append(nodemap)
    return _watch.arm(dicts, before)


def _row(a, n):
    """A history row as the C reads it, or None."""
    if type(a) is np.ndarray and (a.dtype is _F64 or a.dtype == np.float64) and a.ndim == 1 and a.shape[0] == n:
        return a if a.flags.c_contiguous else np.ascontiguousarray(a)
    return None


_NTK = []


def _ntk():
    """`NumericToolkit`, imported at the first call (a function-level import
    runs importlib's `_handle_fromlist` on every call)."""
    from pycircuit.circuit.toolkit import NumericToolkit
    _NTK.append(NumericToolkit)
    return NumericToolkit


#: the circuit's passes for the stage methods from the core (`passes`);
#: env `PYCIRCUIT_CORE_PASSES=0` turns it off
CORE_PASSES = os.environ.get('PYCIRCUIT_CORE_PASSES', '1') != '0'
#: each pass's want bit (q is always made)
_PASS_BITS = {'G': _G, 'C': _C, 'i': _I, 'q': 0}


def passes(tr, x, which):
    """`{m: cir.m(x, epar)}` for each pass `m` of `which` (letters of
    'GCiq') -- fresh float64 arrays, bit for bit the circuit's passes --
    from ONE core call (formula 4: the passes alone, through the classes'
    fused kernels), or None where the core does not serve the call (nothing
    touched; the caller makes its own calls): the switches, a toolkit that
    is not numeric, an instance shadow of a pass on the circuit, a circuit
    the core cannot serve, a state that is not one finite float64 vector of
    its size, a batch not ready -- each counted (`core.p:<why>`) as
    `evaluate` counts its own.  For the stage methods (Radau, the DIRKs),
    whose step evaluates the passes at several points a step."""
    K = _K['p']
    if not (CORE and CORE_PASSES):
        return _no(K['off'])
    if type(tr.toolkit) is not (_NTK[0] if _NTK else _ntk()):
        return _no(K['toolkit'])
    cd = tr.cir.__dict__
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd:
        return _no(K['shadow_cir'])
    ## A CIRCUIT THE CORE CANNOT SERVE declines at once under the plan its
    ## dict holds (`_passes_no`, as `core_for` remembers it per plan),
    ## without `core_for`'s full check, which such a circuit's stamps need
    ## not spare.  A plan gone stale is replaced at the circuit's next
    ## evaluation, which follows; until then the decline is the circuit's
    ## own passes.
    td = tr.__dict__
    kept = td.get('_passes_no')
    if kept is not None and cd.get('_stamp_plan') is kept:
        return _no(K['unservable'])
    core = core_for(tr)
    if core is None:
        plan = cd.get('_stamp_plan')
        if plan is not None:
            td['_passes_no'] = plan
        return _no(K['unservable'])
    n = core.n
    if not (type(x) is np.ndarray and (x.dtype is _F64 or x.dtype == np.float64) and x.ndim == 1
            and x.flags.c_contiguous and x.shape[0] == n and np.isfinite(x).all()):
        return _no(K['x'])
    T = core.ready(tr.epar)
    if T is None:
        return _no(K['ready'])
    bits = 0
    for m in which:
        bits |= _PASS_BITS[m]
    ffi, cfn, _dgemv = driver()
    fb, dptr = ffi.from_buffer, core.dptr
    dummy = core.__dict__.get('_pass_dummy')
    if dummy is None:
        dummy = core._pass_dummy = (np.zeros(max(n * n, 1)), )
        dummy = core._pass_dummy = (dummy[0], fb(dptr, dummy[0]))
    q = np.empty(n)
    C = ((np.zeros if core.npat else np.empty)((n, n))) if bits & _C else None
    dC = fb(dptr, C) if C is not None else dummy[1]
    _PC[K['served']] += 1
    cfn(core.cs, fb(dptr, x), T, bits | (_WANT_PATTERN if core.npat else 0), 4,
        0.0, 0.0, 0.0, 1.0, 0.0,
        dummy[1], dummy[1], dummy[1], dummy[1], dC, dC, fb(dptr, q),
        dummy[1], dummy[1], dummy[1], dummy[1])
    out = {}
    for m in which:
        if m == 'G':
            out['G'] = core.arrays['binG'].reshape(n, n).copy()
        elif m == 'C':
            out['C'] = C
        elif m == 'i':
            out['i'] = core.arrays['bini'].copy()
        else:
            out['q'] = q
    return out


def passes_held(tr, core):
    """Whether a call of `passes` on this transient would find its readiness
    as stamped -- `core_for`'s record and `core.ready`'s probe both holding,
    the switches, the toolkit and the circuit's passes as `passes` checks
    them (for passes another C made in the same step: `_tran_radau_tc`).
    Nothing counted: `passes_count` makes the
    counters of each such call -- a check that holds holds for every call
    the caller makes before anything it watches can change."""
    if not (CORE and CORE_PASSES) or type(tr.toolkit) is not (_NTK[0] if _NTK else _ntk()):
        return False
    cir = tr.cir
    cd = cir.__dict__
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd:
        return False
    rec = tr.__dict__.get('_tran_core')
    ep = _watch.EPOCH
    if not (rec is not None and rec[1] is core and ep is not None and rec[2] == ep.value
            and cd.get('_stamp_plan') is rec[0] and rec[3]._epoch is rec[0].epoch
            and rec[0].n == len(cir.nodes) + len(cir.branches)):
        return False
    ## (`probe`'s temperature test: a float, an int, or one number)
    T = getattr(tr.epar, 'T', 300.0)
    if type(T) is not float and type(T) is not int and np.ndim(T) != 0:
        return False
    ## (`probe`'s DC-pin test)
    if core.pins and getattr(tr.epar, 'analysis_kind', None) == 'dc':
        return False
    return (core.stamp == ep.value
            and tuple(map(_GETDICT, core.uniq_els)) == core.uniq_dicts
            and all(map(_same_code, core.codes)))


def passes_count(tr, core):
    """The counters one call of `passes` makes on its stamped path (after
    `passes_held`): a held `core_for`, a held probe, a served pass."""
    tr.__dict__['_tran_core'][5] += 1
    core.held += 1
    _PC[_K['p']['served']] += 1


def _companion(tr, n, K):
    """THE STEP'S COMPANION SET UP as `get_diff` sets it up -- the same lines
    as `_tran_newton_c._companion`'s, kept twice on purpose (why:
    `doc/transient_history.md`, `_tran_core.evaluate`): a change here is a
    change there.  `(h, h_last, the active integrator, its formula, (q1,
    q2, iq1), (a0, a1, a2, theta))`, or None: declined (counted, `K`)."""
    h = tr._dt
    h_last = tr._dt_last if tr._dt_last is not None else h
    active = tr.base_integrator.check_order_drop(h, h_last, tr._is_first_step)
    formula = FORMULA.get(type(active).__name__)
    if formula is None:
        return _no(K['formula'])
    q1 = _row(tr._qlast[0], n)
    q2 = _row(tr._qlast[1], n) if formula == _GEAR2 else q1
    iq1 = _row(tr._iqlast[0], n) if formula in _IQ_FORMS else q1
    if q1 is None or q2 is None or iq1 is None:
        return _no(K['history'])
    a0 = a1 = a2 = theta = 0.0
    if formula == _GEAR2:
        a0, a1, a2 = _lte_kernels.bdf2_alphas(h, h_last)
    elif formula == _THETA:
        theta = active.theta_at(h)
    return h, h_last, active, formula, (q1, q2, iq1), (a0, a1, a2, theta)


def _c_inputs(tr, x, want, n, core, K):
    """The call's passes (`WANT`'s bits) and its C: 'f' the held `_Cmat`;
    else `_C_lookup`'s at x where it has one (no C pass), or a fresh buffer
    the C pass writes.  `(bits, the companion's C in, C out)`, or None:
    declined (counted, `K`)."""
    bits = WANT[want]
    if want == 'f':
        Cin = tr._Cmat
        if not (type(Cin) is np.ndarray and (Cin.dtype is _F64 or Cin.dtype == np.float64) and Cin.shape == (n, n)
                and Cin.flags.c_contiguous):
            return _no(K['cmat'])
        C = Cin
    else:
        Cin = tr._C_lookup(x)
        if Cin is not None:
            if not (type(Cin) is np.ndarray and (Cin.dtype is _F64 or Cin.dtype == np.float64)
                    and Cin.shape == (n, n) and Cin.flags.c_contiguous):
                return _no(K['clookup'])
            bits &= ~_C
            C = Cin
        else:
            C = (np.zeros if core.npat else np.empty)((n, n))
            Cin = C
    return bits, Cin, C


def _leave(tr, x, want, bits, core, C, q, iq, Geq, active, h, h_last):
    """The state the Python path leaves behind: the memo's C, the q and C
    caches, the C's pattern (all but 'f'), and the step's companion as
    `get_diff` leaves it."""
    ## the state the Python path leaves behind
    if want != 'f':
        if bits & _C and tr._memo_ok():
            tr._memo_put(x, {'C': C})
        tr._q_cache = (x, q)
        tr._C_cache = (x, C)
        if bits & _C and core.npat:
            ## (this C the core wrote on its pattern: the branch screen reads
            ## only the pattern's entries of it -- `_tran_branch_c`)
            tr._C_pattern = (C, core.arrays['pat'])
    ## (the step's companion left as `get_diff` leaves it -- the same writes
    ## as `_tran_newton_c._serve`'s; `transient._RUN_STATE`'s 'step' block)
    tr.active_integrator = active
    tr._iq = iq
    tr._Geq = Geq
    tr._Cmat = C
    tr._companion_coeffs = active.companion_coefficients(h, h_last)
    tr._effective_method = type(active).__name__


def evaluate(tr, x, t, provided_function, want):
    """`want` 'fj': `(f, J)` as `_residual_and_jacobian`; 'j': `(None, J)`
    as `jacobian_only`; 'c': `(None, None)`, `jacobian_only`'s state where
    nothing reads its `J` (C and q; no G); 'f': `f` as `residual_only` -- or
    None where the core does not serve the call (nothing touched; the
    Python path runs).  A served call leaves the state the Python path
    leaves (module docstring: FOUR SITES)."""
    K = _K[want]
    if not CORE:
        return _no(K['off'])
    if type(tr.toolkit) is not _toolkit.NumericToolkit:
        return _no(K['toolkit'])
    ## THE METHODS THE CORE STANDS IN FOR MUST BE THE ONES IT MIRRORS: an
    ## instance shadow of the circuit's passes (a test counting `cir.G`, a
    ## harness timing it) or of the transient's companion machinery (a spy
    ## on `get_diff`) is a caller that expects those calls -- the same rule
    ## as a batch's on an element's shadow (`_hdl_batch`)
    cd = tr.cir.__dict__
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd:
        return _no(K['shadow_cir'])
    td = tr.__dict__
    if ('get_diff' in td or '_companion_at' in td or '_C_at_state' in td
            or '_C_lookup' in td or '_source_at' in td):
        return _no(K['shadow_tr'])
    core = core_for(tr)
    if core is None:
        return _no(K['unservable'])
    n = core.n
    if not (type(x) is np.ndarray and (x.dtype is _F64 or x.dtype == np.float64) and x.ndim == 1
            and x.flags.c_contiguous and x.shape[0] == n and np.isfinite(x).all()):
        return _no(K['x'])
    comp = _companion(tr, n, K)
    if comp is None:
        return None    # paths: counted by _companion
    h, h_last, active, formula, (q1, q2, iq1), (a0, a1, a2, theta) = comp
    cin = _c_inputs(tr, x, want, n, core, K)
    if cin is None:
        return None    # paths: counted by _c_inputs
    bits, Cin, C = cin
    T = core.ready(tr.epar)
    if T is None:
        return _no(K['ready'])
    u = None
    if bits & _I:
        u = tr._source_at(t, provided_function)
        if not (type(u) is np.ndarray and (u.dtype is _F64 or u.dtype == np.float64) and u.ndim == 1
                and u.shape[0] == n and u.flags.c_contiguous):
            return _no(K['u'])
    _PC[K['served']] += 1
    ffi, cfn, _dgemv = driver()
    fb, dptr = ffi.from_buffer, core.dptr
    ## (the matrices zero off the pattern where the core writes only that)
    mat = np.zeros if core.npat else np.empty
    q, iq, Geq = np.empty(n), np.empty(n), mat((n, n))
    F = np.empty(n) if bits & _I else q
    J = mat((n, n)) if bits & _G else Geq
    cfn(core.cs, fb(dptr, x), T, bits | (_WANT_PATTERN if core.npat else 0), formula, a0, a1, a2,
        float(h), theta,
        fb(dptr, q1), fb(dptr, q2), fb(dptr, iq1), fb(dptr, u if u is not None else q1),
        fb(dptr, Cin), fb(dptr, C), fb(dptr, q), fb(dptr, iq), fb(dptr, Geq),
        fb(dptr, F), fb(dptr, J))
    _leave(tr, x, want, bits, core, C, q, iq, Geq, active, h, h_last)
    if want == 'fj':
        return F, J
    if want == 'j':
        return None, J
    if want == 'c':
        return None, None
    return F
