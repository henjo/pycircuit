"""PCNR'S DC OPERATING POINT IN C: one `pcnr.solve_dc` as one call of
`csrc/pcnrdc.c` -- this module is its glue, and the printer of the one
piece the unit cannot share with the rest of the C runtime, each class's
REFINE (`refine_kernel`).

THE SOLVE.  `solve` is called by `pcnr.solve_dc` once its start is set
(x, the limited unknowns, the sources scaled) and makes the whole loop in
C: the ordinary part on the circuit's VIEW (`pcnr.pcnr_view`) by the C
core built on it (`_tran_core`), each participant's `i` and `G` -- the
very kernels its `el.i` / `el.G` call (`_hdl_cbackend.CKernel`) -- at the
local point its probe tree builds (`PcnrDevice._kernel_map`), the Schur
step, numpy's own dgesv, the refine and the three tests: bit for bit the
Python's loop (pinned in `tests/test_pinned_pairs.py`).  OWNS, on the
circuit: `_pcnr_c` (a `_Rec`, rebuilt where the view or its stamp plan
is).

THE REFINE.  `PcnrDevice.limit` runs a class's generated `pcnr_limit`:
each law probe's parameters (SPICE's `von` among them: a parameter that
reads the solution, at the device's last accepted sub-vector), its law
(`_limiting.apply_limit`) and `pcnr.limit_block`'s rule over the device's
unknowns.  Here that is ONE C kernel per class, bit for bit the Python's:
the parameter chains printed by `_hdl_climit._parameter_chains` (the same
printer, the same symbol map), the laws of `csrc/limit.c`, and
`limit_block` -- each probe's value as a loop over the unknowns in order,
the per-probe laws, the Kruskal and potentials of `csrc/pcnr_refine.c`
where a redundant probe over-determines them.  `int hdl_fn(const double
*x, const double *p, double *out)`, the chain functions' entry: `x` the
device's last accepted sub-vector, `p` its packed parameters (the
temperature at its slot), `out` the new unknowns (limited in place) and
then the old ones; 1 declines -- a NaN among the keys to sort, whose order
the Python's `sorted` decides -- and the Python's refine answers.

EXACT OR NOT AT ALL.  A flag numpy would report, a non-finite value, a
failed solve or a declined refine hands the solve back -- nothing of the
caller's written -- and the Python's loop runs from the start: its
warnings, its exceptions.  The iteration limit is served: the Python's
own RuntimeError, from the last iteration's arrays.

DECLINES (`_paths.no`, prefix `pcnr_c:`): the switch, no build, gmin (the
damper), a solver that is not numpy's dense solve (none given, or the
DC's AutoSolver choosing its dense solve), a toolkit that is not numeric,
numpy not ignoring underflow, a participant that is scalar, whose class
the kernels do not serve (`_class_ok`), or whose instance holds its own
`i`, `G` or `pcnr_limit` (`_OWN`: the Python's `el.i(pot)` reads it, the C
its class's kernels and law), an evaluation session over a
class not C-bound, a circuit the core cannot serve or its readiness, and
inputs that are not float64 vectors of their sizes.  A singular solve is
served too: numpy's LinAlgError, its own message (`_singular_text`).  Switch:
`PYCIRCUIT_PCNR_C=0` (`ENABLED`), read on every call.
"""
import os
import types

import numpy as np

from pycircuit.circuit import _cbuild, _paths
from pycircuit.circuit import _hdl_cbackend as _cb

ENABLED = os.environ.get('PYCIRCUIT_PCNR_C', '1') != '0'
STATUS = 'not loaded'

#: what the C stands in for on a participant: its `i` and `G` (the class's
#: kernels) and its limiter (the class's refine law) -- an instance holding
#: its own is the Python's, which reads it (`PcnrDevice._kernel_stamp`'s
#: `el.i`, `el.G`; `PcnrDevice.limit`'s `pcnr_limit`)
_OWN = ('i', 'G', 'pcnr_limit')

_no = _paths.no
_PC = _paths.COUNTS
_F64 = np.dtype(np.float64)
_SOURCES = None
_driver = None


def _sources():
    """`(the unit's source, its declaration, its codes, its bails' names,
    pcx_abi.h's codes)`, read at the first need."""
    global _SOURCES
    if _SOURCES is None:
        c = _cbuild.csrc('pcnrdc.c')
        h = _cbuild.csrc('pcx_klu.h') + _cbuild.csrc('pcnrdc.h')
        codes = _cbuild.constants(h)
        _SOURCES = (c, h, codes, _cbuild.bails(codes, 'PCNRDC_BAIL_'),
                    _cbuild.constants(_cbuild.csrc('pcx_abi.h')))
    return _SOURCES


def __getattr__(name):
    """`PCNRDC_C` and `PCNRDC_CDEF`, the unit's source and declaration
    (`_sources`), and `REFINE_C`, the refine kernels' prelude."""
    if name == 'PCNRDC_C':
        return _sources()[0]
    if name == 'PCNRDC_CDEF':
        return _sources()[1]
    if name == 'REFINE_C':
        return _refine_prelude()
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')


def driver():
    """`(ffi, cfn, dgesv, dgemv)` once loaded, or None (`STATUS`)."""
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
        except (OSError, AttributeError) as e:
            STATUS = f'off ({e})'
            return None
        try:
            src, cdef = _sources()[:2]
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(src, cdef, 'pcx_pcnrdc')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, dgesv, dgemv)
        STATUS = 'c'
    return _driver or None


## ---------------------------------------------------------------------------
## The refine kernel.
## ---------------------------------------------------------------------------

#: texts read once: 'refine' the refine kernels' prelude, 'singular'
#: numpy's message (a dict, not module scalars: the leak detector reports a
#: module scalar that moves)
_TEXT = {}


def _refine_prelude():
    """`csrc/limit.c` (the laws, the ranking order, union-find) and
    `csrc/pcnr_refine.c` (the Kruskal and potentials), read once."""
    if 'refine' not in _TEXT:
        from pycircuit.circuit import _hdl_climit
        _TEXT['refine'] = _hdl_climit._LIMIT_C + _cbuild.csrc('pcnr_refine.c')
    return _TEXT['refine']


class Refused(Exception):
    """A class whose refine the kernel does not take (the reason)."""


def _laws_of(cls):
    """`(laws, probes, coeffs)` of the class's generated `pcnr_limit`: per
    law probe `(kind, parameter functions)`, its local row pair, its
    coefficients over the unknowns."""
    fn = getattr(cls, 'pcnr_limit', None)
    if fn is None or not fn.__defaults__:
        raise Refused('the class has no generated pcnr_limit')
    code = fn.__code__
    names = code.co_varnames[code.co_argcount - len(fn.__defaults__):code.co_argcount]
    d = dict(zip(names, fn.__defaults__))
    if not {'_laws', '_probes', '_coeffs'} <= set(d):
        raise Refused('pcnr_limit carries no laws')
    return d['_laws'], d['_probes'], d['_coeffs']


def _block_c(probes, coeffs, m):
    """Each probe's value over the unknowns -- `raw[p] = 0.0; raw[p] +=
    c * vn[j]` for every nonzero coefficient in the unknowns' order, `old`
    likewise (`pcnr.limit_block`'s loop: the sum's order and a zero's sign
    the Python's)."""
    lines = []
    for p in range(len(probes)):
        for name, v in (('raw', 'vn'), ('old', 'vo')):
            lines.append(f'  {name}[{p}] = 0.0;')
            for j in range(m):
                c = float(coeffs[p][j])
                if c != 0.0:
                    lines.append(f'  {name}[{p}] += ({c!r}) * {v}[{j}];')
    return lines


def _chains(cls, spec, n_p, t_index):
    """`_hdl_climit._parameter_chains`'s text -- its SymPy printing, most of
    a small deck's first PCNR solve -- kept beside the class's compiled
    info (`_hdl_cache.side_text`), keyed by the printers' own sources."""
    from pycircuit.circuit import _hdl_cache, _hdl_climit
    tag = (f'pcnr refine chains {_hdl_cache._file_hash(_hdl_climit.__file__)} '
           f'{_hdl_cache._file_hash(__file__)} {(n_p, t_index)!r}')
    return _hdl_cache.side_text(
        cls, tag, lambda: _hdl_climit._parameter_chains(spec, n_p, t_index))


def render_refine(cls):
    """`(csrc, nx, layout)` of the class's refine kernel, or `Refused`."""
    from pycircuit.circuit import _hdl_climit, hdl
    info = getattr(cls, '_hdl_info', None) or {}
    fi = (info.get('funcs') or {}).get('i')
    if fi is None or getattr(fi, '_csrc', None) is None:
        raise Refused('the class carries no C source')
    laws, probes, coeffs = _laws_of(cls)
    m, npr = len(coeffs[0]) if coeffs else 0, len(probes)
    nx = int(fi._cshape[0])
    n_p, t_index = fi._clayout
    if not 0 < m <= npr <= 32 or nx > 32:
        raise Refused('more rows or probes than the kernel holds')
    spec = [(tuple(probes[j]), kind, 0, pfs) for j, (kind, pfs) in enumerate(laws)]
    try:
        text = _chains(cls, spec, n_p, t_index)
    except _hdl_climit.Refused as e:
        raise Refused(str(e))
    npf = [len(pfs) for _k, pfs in laws]
    out = ['/* PCNR\'s refine for one device of this class: x its last accepted',
           '   sub-vector, p its packed parameters, out its new unknowns (limited',
           '   in place) then its old ones; 1 declines (a NaN key) */',
           f'int {hdl._C_ENTRY}(const double *x, const double *p, double *out) {{',
           f'  enum {{ M = {m}, NP = {npr}, NR = {nx} }};',
           text,
           f'  static const int A[NP] = {_hdl_climit._c_list(r[0] for r in probes)};',
           f'  static const int B[NP] = {_hdl_climit._c_list(r[1] for r in probes)};',
           f'  static const int KIND[NP] = {_hdl_climit._c_list(_hdl_climit.KINDS[k] for k, _p in laws)};',
           '  const double P0[NP] = {' + ', '.join(
               f'pv_{j}_0' if npf[j] > 0 else '0.0' for j in range(npr)) + '};',
           '  const double P1[NP] = {' + ', '.join(
               f'pv_{j}_1' if npf[j] > 1 else '0.0' for j in range(npr)) + '};',
           '  double vn[M], vo[M], raw[NP], old[NP], lim[NP];',
           '  int bit = 0;',
           '  for (int j = 0; j < M; j++) { vn[j] = out[j]; vo[j] = out[M + j]; }']
    out += _block_c(probes, coeffs, m)
    out += ['  for (int q = 0; q < NP; q++) {',
            '    lim[q] = _lim(KIND[q], raw[q], old[q], P0[q], P1[q]);',
            '    if (lim[q] != raw[q]) bit = 1;',
            '  }',
            '  if (!bit) return 0;']
    if npr == m:
        out += ['  for (int j = 0; j < M; j++) out[j] = lim[j];', '  return 0;']
    else:
        out += ['  return _pcnr_tree(M, NP, A, B, raw, lim, NR, out);']
    out.append('}')
    return _refine_prelude() + '\n'.join(out) + '\n', nx, (n_p, t_index)


def refine_kernel(cls):
    """The class's refine kernel -- `(ffi, cfn, nx, layout)` -- bound once on
    its info (`_c_pcnr_refine`, `_c_pcnr_refine_status`), or None where it
    is refused or does not build."""
    from pycircuit.circuit import _hdl_climit
    info = cls._hdl_info
    if '_c_pcnr_refine' in info:
        return info['_c_pcnr_refine']
    try:
        csrc, nx, layout = render_refine(cls)
        ffi, cfn, _key, _cold, _secs = _cb.load_kernel(csrc, _hdl_climit.limit_cdef())
    except Refused as e:
        info['_c_pcnr_refine'], info['_c_pcnr_refine_status'] = None, f'refused ({e})'
        return None
    except (_cb.CompileError, OSError) as e:
        info['_c_pcnr_refine'], info['_c_pcnr_refine_status'] = None, f'build failed ({e})'
        return None
    kern = info['_c_pcnr_refine'] = types.SimpleNamespace(ffi=ffi, cfn=cfn, nx=nx,
                                                          layout=layout)
    info['_c_pcnr_refine_status'] = 'c'
    return kern


def _packed(el):
    """The element's packed parameters as its kernels read them (`CKernel`'s
    `_hdl_cp`, made where absent), or None."""
    packed = el.__dict__.get('_hdl_cp')
    if packed is None:
        try:
            packed = type(el)._hdl_info['funcs']['i']._hdl_c.pack(el)
        except (TypeError, ValueError, AttributeError):
            ## (AttributeError: a class whose `i` has no C binding -- one
            ## the C solve declines, asked for its refine alone)
            packed = False
        el.__dict__['_hdl_cp'] = packed
    return packed or None


def refine(dev, v_new, v_old, epar, x_old_sub):
    """One device's refine through its class's kernel: the limited unknowns
    (a fresh array) -- `PcnrDevice.limit`'s answer bit for bit -- or None
    (the kernel declined or is absent).  The test hook of the kernel; the
    C unit calls it directly."""
    kern = refine_kernel(type(dev.element))
    packed = _packed(dev.element)
    if kern is None or packed is None:
        return None
    p, pcast = packed
    t_index = kern.layout[1]
    if t_index is not None:
        p[t_index] = getattr(epar, 'T', 300.0)
    m = dev.m
    io = np.empty(2 * m)
    io[:m], io[m:] = v_new, v_old
    xs = np.ascontiguousarray(x_old_sub, dtype=float)
    ffi = kern.ffi
    if kern.cfn(ffi.from_buffer('double *', xs), pcast, ffi.from_buffer('double *', io)):
        return _no('pcnr_c:refine_nan_key')
    return io[:m].copy()


## ---------------------------------------------------------------------------
## The solve.
## ---------------------------------------------------------------------------

def _class_ok(cls, t):
    """The class's kernels -- `(i, G, refine)`, each `(ffi, cfn)` -- where
    `el.i` / `el.G` reach them on every DC call of a device with `t` rows
    (the chained path, no DC pin, C-bound, the hdl methods themselves, one
    packed layout), else a reason (a str).  Cached on the class's info
    per row count, for the `i` and `G` kernel objects bound now: a detach
    or a re-attach (`_hdl_cbackend`) binds others, and the answer is
    made again."""
    from pycircuit.circuit import hdl
    info = getattr(cls, '_hdl_info', None)
    if info is None:
        return 'class'
    funcs = info.get('funcs') or {}
    fi, fG = funcs.get('i'), funcs.get('G')
    ki = fi.__dict__.get('_hdl_c') if fi is not None else None
    kG = fG.__dict__.get('_hdl_c') if fG is not None else None
    cache = info.setdefault('_c_pcnr_kernels', {})
    hit = cache.get(t)
    if hit is not None and hit[0] is ki and hit[1] is kG:
        return hit[2]
    got = 'class'
    meth = [cls.__dict__.get(nm) for nm in ('i', 'G')]
    if not (info.get('chained') and info.get('_c_bound') and ki is not None and kG is not None
            and not info['state_meta']['dc_pins']
            and all(getattr(f, '__module__', None) == hdl.__name__ for f in meth)):
        pass
    elif not (ki.shape == (t,) and tuple(kG.shape) in ((t, t), (t * t,)) and ki.nx <= t
              and kG.nx <= t and (ki.n_p, ki.t_index) == (kG.n_p, kG.t_index)):
        got = 'shape'
    else:
        rk = refine_kernel(cls)
        if rk is None or rk.layout != (ki.n_p, ki.t_index):
            got = 'refine'
        else:
            got = ((ki.ffi, ki.cfn), (kG.ffi, kG.cfn), (rk.ffi, rk.cfn), ki.t_index)
    cache[t] = (ki, kG, got)
    return got


def _addr(ffi, kffi, cfn):
    return ffi.cast('void *', int(kffi.cast('uintptr_t', cfn)))


def _tables(devices, kerns, epar):
    """The participants' tables, each device's slice at its offset
    (pcnrdc.h's THE PARTICIPANTS): `(the int64 tables, the float ones --
    T and the remainders --, the blocks' length)`."""
    I64 = np.int64
    tabs = {nm: [] for nm in ('d_t', 'd_m', 'd_off', 'd_rows', 'rows', 'prb', 'd_walk',
                              'n_walk', 'walk', 'd_mir', 'n_mir', 'mir', 'd_tm', 'd_lin',
                              'ti')}
    Tm, lin = [], []
    ntm = nwalk = nmir = nrows = nlin = 0
    for d, kern in zip(devices, kerns):
        G_lin, walk, mirror, T = d._kernel_map(d.params(), epar)
        tabs['d_t'].append(d.t)
        tabs['d_m'].append(d.m)
        tabs['d_off'].append(d.off)
        tabs['d_rows'].append(nrows)
        tabs['rows'] += d.rows
        nrows += d.t
        for ra, rb in d.probes:
            tabs['prb'] += [ra, rb]
        tabs['d_walk'].append(nwalk)
        tabs['n_walk'].append(len(walk))
        for row, parent, j, sign in walk:
            tabs['walk'] += [row, parent, j, int(sign)]
        nwalk += len(walk)
        tabs['d_mir'].append(nmir)
        tabs['n_mir'].append(len(mirror))
        for row, src in mirror:
            tabs['mir'] += [row, src]
        nmir += len(mirror)
        tabs['d_tm'].append(ntm)
        Tm.append(np.asarray(T, dtype=float).ravel())
        ntm += d.t * d.m
        if G_lin is None:
            tabs['d_lin'].append(-1)
        else:
            tabs['d_lin'].append(nlin)
            lin.append(np.asarray(G_lin, dtype=float).ravel())
            nlin += d.t * d.t
        tabs['ti'].append(-1 if kern[3] is None else kern[3])
    arrs = {nm: np.ascontiguousarray(v if v else [0], dtype=I64) for nm, v in tabs.items()}
    arrs_f = {'Tm': np.concatenate(Tm) if Tm else np.zeros(1),
              'lin': np.concatenate(lin) if lin else np.zeros(1)}
    return arrs, arrs_f, ntm


class _Rec:
    """One circuit's struct (`s`, a `pcnrdc_t`), buffers and device tables
    for one view, stamp plan, participant set, reference row and set of
    class kernels (`kerns`, `_class_ok`'s), set once: a solve copies its
    inputs in and points `prow` at the packs."""

    __slots__ = (
        'P',
        'buf',
        'core',
        'els',
        'ffi',
        'insts',
        'iref',
        'k',
        'keep',
        'kerns',
        'klu',
        'klu_keep',
        'kmaps',
        'n',
        'packs',
        'plan',
        'prow',
        'red',
        's',
        'view',
    )

    def __init__(self, cir, view, plan, core, devices, kerns, iref, drv, core_drv, epar):
        ffi, _cfn, dgesv, dgemv = drv
        cffi_core, cfn_core, _dg = core_drv
        n = cir.n
        k = sum(d.m for d in devices)
        nd = len(devices)
        self.view, self.plan, self.core, self.iref, self.n, self.k = view, plan, core, iref, n, k
        self.insts = tuple(d.instance for d in devices)
        self.kerns = tuple(kerns)
        self.ffi = ffi
        self.s = s = ffi.new('pcnrdc_t *')
        I64 = np.int64
        arrs, arrs_f, ntm = _tables(devices, kerns, epar)
        tmax = max(d.t for d in devices)
        mmax = max(d.m for d in devices)
        m = n - 1
        E, Z = np.empty, np.zeros
        self.buf = buf = {
            'u': E(n), 'x': E(n), 'v': E(max(k, 1)), 'fe': E(n), 'dx': E(n), 'gl': E(max(k, 1)),
            ## (J zero: off the pattern nothing writes it -- pcnrdc.h)
            'F': E(n), 'J': Z(n * n), 'Jabs': E(n * n), 'rc': E(n), 'A': E(max(m * m, 1)),
            'b': E(max(m, 1)), 'xn': E(n), 'ax': E(n), 'Is': E(n), 'vn': E(max(k, 1)),
            'R': E(max(m * m, 1)), 'kx': E(max(m, 1)), 'ky': E(max(m, 1)), 'kt': E(max(2 * m, 1)),
            'blk': E(max(ntm, 1)), 'pot': E(tmax), 'cur': E(tmax), 'Gl': E(tmax * tmax),
            'xs': E(tmax), 'rio': E(2 * mmax), 'Cs': Z(n * n), 'q': E(n),
        }
        fb = ffi.from_buffer
        keep = []
        for nm, a in buf.items():
            h = fb('double *', a)
            keep += [a, h]
            setattr(s, nm, h)
        for nm, a in arrs.items():
            h = fb('int64_t *', a)
            keep += [a, h]
            setattr(s, nm, h)
        for nm, a in arrs_f.items():
            h = fb('double *', a)
            keep += [a, h]
            setattr(s, nm, h)
        ipiv = np.empty(max(m, 1), dtype=I64)
        keep += [ipiv, fb('int64_t *', ipiv)]
        s.ipiv64 = keep[-1]
        fns = {nm: ffi.new('void *[]', max(nd, 1)) for nm in ('ifn', 'Gfn', 'rfn')}
        for j, kern in enumerate(kerns):
            for nm, (kffi, cfn) in zip(('ifn', 'Gfn', 'rfn'), kern[:3]):
                fns[nm][j] = _addr(ffi, kffi, cfn)
        for nm, h in fns.items():
            setattr(s, nm, h)
        self.prow = ffi.new('double *[]', max(nd, 1))
        s.prow = self.prow
        keep += list(fns.values())
        self.keep = keep
        self.els = tuple(d.element for d in devices)
        self.packs = [None] * nd
        s.dgesv, s.dgemv = ffi.cast('void *', dgesv), ffi.cast('void *', dgemv)
        s.core = ffi.cast('void *', int(cffi_core.cast('uintptr_t', core.cs)))
        s.core_fn = ffi.cast('void *', int(cffi_core.cast('uintptr_t', cfn_core)))
        s.pattern = _sources()[4]['PCX_WANT_PATTERN'] if core.npat else 0
        s.n, s.iref, s.nn, s.nd, s.k = n, iref, len(cir.nodes), nd, k
        ## each row's reduced index (`solve` sets it per call: the reference
        ## and the held rows out)
        self.red = np.empty(n, dtype=I64)
        keep += [self.red, fb('int64_t *', self.red)]
        s.red = keep[-1]
        ## the junction pairs `jgmin` is put across (`pcnr_junction_pairs`)
        from pycircuit.circuit import pcnr as _pcnr
        pairs = _pcnr.pcnr_junction_pairs(cir)
        for nm, col in (('ja', 2), ('jb', 3)):
            a = np.array([p[col] for p in pairs], dtype=I64).reshape(-1)
            h = fb('int64_t *', a)
            keep += [a, h]
            setattr(s, nm, h)
        s.nj, s.gjgmin = len(pairs), 0.0
        ## THE PATTERN (pcnrdc.h): the core's, each Schur position -- a
        ## device row and a probe's ra or rb -- and the junction pairs'
        ## entries; none where the core writes every entry
        self.P, self.kmaps = None, {}
        s.np, s.P = 0, ffi.NULL
        if core.npat:
            parts = [np.asarray(core.arrays['pat'], dtype=I64),
                     np.arange(n, dtype=I64) * (n + 1)]
            for d in devices:
                rows = np.asarray(d.rows, dtype=I64) * n
                for ra, rb in d.probes:
                    parts += [rows + ra, rows + rb]
            for p in pairs:
                a, b = int(p[2]), int(p[3])
                parts.append(np.array([a * n + a, a * n + b, b * n + a, b * n + b], dtype=I64))
            self.P = np.ascontiguousarray(np.unique(np.concatenate(parts)), dtype=I64)
            keep += [self.P, fb('int64_t *', self.P)]
            s.np, s.P = len(self.P), keep[-1]
        self.klu, self.klu_keep = None, ()
        s.klu.klu_sym = ffi.NULL
        s.kfull, s.nfor, s.ffull = ffi.NULL, 0, ffi.NULL

    def bind(self, ch):
        """The solve: None the dense one, else `ch` -- an analysed
        StructuralKLU of this size: its pattern (the entries `take` of the
        reduced matrix), arrays, analysis, common and the library's entries
        (D4's `_klu_bind`)."""
        ffi, s = self.ffi, self.s
        s.kfull, s.nfor, s.ffull = ffi.NULL, 0, ffi.NULL
        if ch is None:
            self.klu, self.klu_keep = None, ()
            s.klu.klu_sym = ffi.NULL
            return
        if self.klu is ch and s.klu.klu_sym == ffi.cast('void *', ch._symbolic):
            return
        take = np.ascontiguousarray(ch.take, dtype=np.int64)
        kAx = np.empty(max(len(take), 1))
        keep = (ch, take, kAx, ffi.from_buffer('int64_t *', take),
                ffi.from_buffer('double *', kAx), ffi.from_buffer('int *', ch.Ap),
                ffi.from_buffer('int *', ch.Ai))
        s.nact, s.act_r, s.klu.kAx = len(take), keep[3], keep[4]
        k, lib = s.klu, ch._lib
        k.kAp, k.kAi = keep[5:]
        k.klu_factor, k.klu_refactor, k.klu_solve, k.klu_free = (
            ffi.cast('void *', _cbuild.address(f))
            for f in (lib.klu_factor, lib.klu_refactor, lib.klu_solve, lib.klu_free_numeric))
        k.klu_common = ffi.cast('void *', ch._ct.addressof(ch._common))
        k.klu_sym = ffi.cast('void *', ch._symbolic)
        k.resid_tol = float(ch.REFACTOR_RESIDUAL_TOL)
        self.klu, self.klu_keep = ch, keep

    def on_pattern(self, held):
        """KLU on the pattern (pcnrdc.h's KLU ON THE PATTERN) for the bound
        solve and this call's reduction: each KLU entry's index in J, and the
        pattern's entries in the reduced matrix off KLU's -- made once per
        solver and held set.  Nothing where there is no pattern or no KLU."""
        s, ch = self.s, self.klu
        if ch is None or self.P is None:
            return
        key = (id(ch), tuple(held))
        got = self.kmaps.get(key)
        if got is None or got[0] is not ch:
            n, m, P, red = self.n, int(s.m), self.P, self.red
            take = self.klu_keep[1]
            unred = np.flatnonzero(red >= 0).astype(np.int64)
            kfull = np.ascontiguousarray(unred[take // m] * n + unred[take % m], dtype=np.int64)
            r, c = P // n, P % n
            act = (red[r] >= 0) & (red[c] >= 0)
            pos = red[r[act]] * m + red[c[act]]
            ffull = np.ascontiguousarray(P[act][~np.isin(pos, take)], dtype=np.int64)
            fb = self.ffi.from_buffer
            got = self.kmaps[key] = (ch, kfull, ffull, fb('int64_t *', kfull),
                                     fb('int64_t *', ffull) if len(ffull) else self.ffi.NULL)
        s.kfull, s.nfor, s.ffull = got[3], len(got[2]), got[4]

    def write_back(self):
        """StructuralKLU as its own solve leaves it after a served outcome:
        the numeric the sequence ends with (None after a failed factor), its
        counters."""
        sk, k = self.klu, self.s.klu
        if sk is None:
            return
        sk._numeric = int(self.ffi.cast('uintptr_t', k.klu_out)) or None
        sk.refactors += k.refactors
        sk.factors += k.factors
        sk.pivot_fallbacks += k.pivot_fallbacks
        sk.residual_fallbacks += k.residual_fallbacks

    def point_packs(self):
        """`prow` at each element's packed parameters as they stand (a
        parameter change drops the pack: made again here), or False."""
        packs, prow, ffi = self.packs, self.prow, self.ffi
        for j, el in enumerate(self.els):
            packed = el.__dict__.get('_hdl_cp')
            if packed is None or packed is False:
                packed = _packed(el)
                if packed is None:
                    return False
            if packed is not packs[j]:
                packs[j] = packed
                prow[j] = ffi.cast('double *', packed[0].ctypes.data)
        return True


def _settings(circuit, gmin):
    """Why the solve's settings take the Python (a reason), or None."""
    from pycircuit.circuit import toolkit as _toolkit
    if gmin:
        return 'gmin'
    if type(circuit.toolkit) is not _toolkit.NumericToolkit:
        return 'toolkit'
    if np.geterr()['under'] != 'ignore':
        return 'underflow'
    return None


def _solver(ls, m):
    """The solve `ls` makes on an m x m matrix: None (numpy's dense solve),
    an analysed StructuralKLU, 'choose' (an AutoSolver that has not chosen:
    its first matrix decides) -- or why not (a reason)."""
    from pycircuit.circuit.linearsolver import AutoSolver, StructuralKLU
    if ls is None:
        return None
    if type(ls) is not AutoSolver:
        return 'linsolver'
    ch = ls._choice
    floor = min(ls.min_n, ls.structural_min_n) if ls.structure is not None else ls.min_n
    if ch is None and m < floor:
        ## (the AutoSolver's own first choice below its floor -- `_select`'s
        ## first test, which reads no value: dense, made as its first solve
        ## would make it, the count included.  Above it the choice reads the
        ## matrix (its fill), so the first Schur matrix makes it: 'choose' --
        ## between a structural floor and `min_n` too, where a zero matrix
        ## would choose a StructuralKLU the real one might not)
        ch = ls._select(np.zeros((m, m)))
    if ch is None:
        return 'choose'
    if ch is ls._dense:
        return None
    if type(ch) is not StructuralKLU or ch.m != m or not ch._symbolic:
        return 'linsolver'
    return ch


def _choose(rec, ls):
    """After PCNRDC_CHOOSE: the AutoSolver's choice on the first reduced
    Schur matrix, as the Python's first solve makes it -- `_select`, and
    for a StructuralKLU its `_ours` test and its analysis (`_factor` calls
    `_analysed` first) -- then the record bound to it.  None, or why not:
    the Python's solve runs ('foreign': `_ours` declines the matrix)."""
    from pycircuit.circuit.linearsolver import StructuralKLU
    n = rec.n
    keep = rec.red >= 0
    R = rec.buf['J'].reshape(n, n)[np.ix_(keep, keep)]
    ch = ls._select(R)
    if ch is ls._dense:
        rec.bind(None)
        return None
    if type(ch) is not StructuralKLU or ch.m != int(rec.s.m):
        return 'linsolver'
    if ch._ours(R) is None:
        return 'foreign'
    ch._analysed()
    rec.bind(ch)
    return None


def _record(circuit, devices, iref, epar, drv):
    """The circuit's record for its view as it stands -- or why not."""
    from pycircuit.circuit import _evalhint, _stamp_plan, _tran_core
    from pycircuit.circuit import pcnr as _pcnr
    if not _tran_core.CORE:
        return 'core'
    kerns = []
    for d in devices:
        if d.scalar or d.rows is None:
            return 'scalar'
        ed = d.element.__dict__
        if any(nm in ed for nm in _OWN):
            return 'own'
        if type(d.element.toolkit) is not type(circuit.toolkit):
            return 'toolkit'
        kern = _class_ok(type(d.element), d.t)
        if type(kern) is str:
            return kern
        kerns.append(kern)
    if _evalhint.current() is not None and not all(
            type(d.element)._hdl_info.get('_c_bound') for d in devices):
        return 'evalhint'
    view = _pcnr.pcnr_view(circuit, devices, epar)
    plan = _stamp_plan._plan_for(view)
    rec = circuit.__dict__.get('_pcnr_c')
    insts = tuple(d.instance for d in devices)
    if (rec is None or rec.view is not view or rec.plan is not plan or rec.iref != iref
            or rec.insts != insts or rec.n != circuit.n
            or any(a is not b for a, b in zip(rec.kerns, kerns))):
        core_drv = _tran_core.driver()
        if core_drv is None:
            return 'core'
        try:
            core = _tran_core._Core(view, plan, core_drv[0], core_drv[2])
        except _tran_core.Unservable:
            return 'unservable'
        rec = circuit.__dict__['_pcnr_c'] = _Rec(circuit, view, plan, core, devices, kerns,
                                                  iref, drv, core_drv, epar)
    return rec


def _vector(a, m):
    return type(a) is np.ndarray and a.dtype is _F64 and a.shape == (m,)


def solve(circuit, irefnode, devices, x, v_lim, u, epar, maxiter, reltol, abstol, iabstol,
          gmin, gshunt, linsolver, jgmin=0.0, held=()):
    """`pcnr.solve_dc`'s loop in C from its start: `(x, v_lim, iterations)`,
    or its exception raised (the iteration limit's RuntimeError, a non-finite
    iterate's, a singular solve's LinAlgError) -- or None: declined or
    handed back, nothing written; the Python's loop runs."""
    if not ENABLED:
        return _no('pcnr_c:off')
    if maxiter < 1:
        return _no('pcnr_c:maxiter')
    why = _settings(circuit, gmin)
    if why is None:
        ch = _solver(linsolver, circuit.n - 1 - len(held))
        why = ch if type(ch) is str and ch != 'choose' else None
    if why is not None:
        return _no('pcnr_c:' + why)
    drv = _driver or driver()
    if not drv:
        return _no('pcnr_c:driver')
    rec = _record(circuit, devices, irefnode, epar, drv)
    if type(rec) is str:
        return _no('pcnr_c:' + rec)
    T = rec.core.probe(epar)
    if T.__class__ is str:
        return _no('pcnr_c:ready')
    if not (_vector(x, rec.n) and _vector(v_lim, rec.k) and _vector(u, rec.n)):
        return _no('pcnr_c:x')
    if not rec.point_packs():
        return _no('pcnr_c:pack')
    s, buf = rec.s, rec.buf
    ## the reduction: the reference row and each held row out (`predict`'s)
    red = rec.red
    red.fill(-1)
    kept = np.ones(rec.n, dtype=bool)
    kept[[rec.iref, *held]] = False
    red[kept] = np.arange(int(kept.sum()), dtype=np.int64)
    s.m = int(kept.sum())
    np.copyto(buf['x'], x)
    np.copyto(buf['v'][:rec.k], v_lim)
    np.copyto(buf['u'], u)
    s.T, s.maxiter = T, int(maxiter)
    s.reltol, s.abstol, s.iabstol, s.gshunt = (float(reltol), float(abstol), float(iabstol),
                                               float(gshunt))
    s.gjgmin = float(jgmin)
    codes = _sources()[2]
    s.choose = ch == 'choose'
    if s.choose:
        rec.bind(None)
        status = drv[1](s)
        if status != codes['PCNRDC_CHOOSE']:
            return _outcome(rec, status, maxiter)
        why = _choose(rec, linsolver)
        if why is not None:
            return _no('pcnr_c:' + why)
        s.choose = 0
    else:
        rec.bind(ch)
    if rec.klu is not None:
        num = rec.klu._numeric
        s.klu.klu_num = rec.ffi.cast('void *', num) if num else rec.ffi.NULL
    rec.on_pattern(held)
    return _outcome(rec, drv[1](s), maxiter)


def _outcome(rec, status, maxiter):
    """The call's outcome as the Python's loop gives it: the solution, its
    exception raised -- StructuralKLU written back as its solve leaves it --
    or None (handed back)."""
    codes, bails = _sources()[2:4]
    s, buf = rec.s, rec.buf
    ## (the verdicts the dense tests made: pcnrdc.h's OUT)
    if s.dense_conv:
        _PC['pcnr_c:dense_conv'] += int(s.dense_conv)
    if s.dense_resid:
        _PC['pcnr_c:dense_resid'] += int(s.dense_resid)
    if status == codes['PCNRDC_CONVERGED']:
        rec.write_back()
        _PC['pcnr_c:served'] += 1
        _PC['pcnr_c:iterations'] += int(s.iters)
        return buf['x'].copy(), buf['v'][:rec.k].copy(), int(s.iters)
    if status == codes['PCNRDC_MAXITER']:
        rec.write_back()
        _PC['pcnr_c:served'] += 1
        _PC['pcnr_c:maxiter'] += 1
        _PC['pcnr_c:iterations'] += int(maxiter)
        g, dx, f = (float(np.max(np.abs(buf['gl'][:rec.k]))), float(np.max(np.abs(buf['dx']))),
                    float(np.max(np.abs(buf['fe']))))
        raise RuntimeError(f'PCNR did not converge in {maxiter:d} iterations: max|g_lim| = '
                           f'{g:g}, max|dx| = {dx:g}, max|f| = {f:g}')
    if status == codes['PCNRDC_SINGULAR']:
        rec.write_back()
        _PC['pcnr_c:served'] += 1
        _PC['pcnr_c:singular'] += 1
        _PC['pcnr_c:iterations'] += int(s.iters) + 1
        ## (StructuralKLU's `_factor` says it in its own words)
        raise np.linalg.LinAlgError('Singular matrix' if rec.klu is not None
                                    else _singular_text())
    if status == codes['PCNRDC_NONFINITE']:
        from pycircuit.circuit.pcnr import NONFINITE
        rec.write_back()
        _PC['pcnr_c:served'] += 1
        _PC['pcnr_c:nonfinite'] += 1
        _PC['pcnr_c:iterations'] += int(s.iters) + 1
        raise RuntimeError(NONFINITE % (int(s.iters) + 1))
    return _no('pcnr_c:bail:' + bails.get(status, str(status)))


def _singular_text():
    """numpy.linalg.solve's own message for a singular matrix, read once
    from numpy itself."""
    if 'singular' not in _TEXT:
        try:
            np.linalg.solve(np.zeros((1, 1)), np.zeros(1))
        except np.linalg.LinAlgError as e:
            _TEXT['singular'] = str(e)
    return _TEXT['singular']
