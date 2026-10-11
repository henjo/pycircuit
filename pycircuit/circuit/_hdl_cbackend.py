# -*- coding: utf-8 -*-
# Copyright (c) 2008-2026 Pycircuit Development Team
# See LICENSE for details.
"""The C backend for chained `Behavioural` elements: each chained
function's numpy source (`fn._src`, the reference) compiled from its C
rendering (`fn._csrc`) and called in its place, bit for bit but for the
named exceptions under Fidelity.

`generate_code` prints every x-taking let-chain twice: the numpy source
(`fn._src`, the reference -- what `explain()` shows, what the compile
cache stores, what every test compares against) and a C rendering
(`fn._csrc`, from `hdl._CChainPrinter`).  This module owns everything
after that print: finding a compiler, building the shared object,
caching it, loading it, and calling it.

OWNS
    On the class, `_hdl_backend_status` (what happened, below) and
    `_hdl_backend_ulp` (the listed functions its C calls, `libm_check`);
    on the class info, `_c_bound`, `_c_fused`, `_c_fused_status`,
    `_backend_pending` and `_backend_seen` (and, through `_hdl_climit`,
    `_c_limit` and `_c_limit_status`); on each chain function, `_hdl_c`
    (its `CKernel`); on each element, `_hdl_cp` (the packed parameters,
    `CKernel.pack`'s, dropped by the element's `update` when a parameter
    moves).  The module keeps `_LIBM` (the C library's check).  `_paths`:
    `once:backend:` (`c`, `numpy`, `compile_failed`: a class's resolve).
    Switches: `PYCIRCUIT_HDL_BACKEND` (`c`, `numpy` or `auto`, read by
    `hdl._backend_requested`) and `PYCIRCUIT_HDL_CFUSE=0` (`FUSE`); the
    build's are `_cbuild`'s and `_hdl_cache`'s (`PYCIRCUIT_HDL_CC`,
    `PYCIRCUIT_C_SANITIZE`, `PYCIRCUIT_HDL_CACHE`).

Selection
    The default is `'auto'`: every chained class runs C where it can be
    served -- the kernels printed (`_csrc`), cffi installed, the compile
    cache on (without a store every process would rebuild every kernel), a
    compiler or the stored objects -- and numpy otherwise, QUIETLY, saying
    why in its status.  It is resolved at the class's FIRST INSTANCE
    (`ensure`, from `Behavioural.__init__`), not at class creation:
    importing the library compiles nothing, a base class that only spawns
    collapse variants never builds, and an instance runs one backend from
    its first evaluation.  An instance on the JAX or symbolic toolkit
    leaves the class unresolved (its own path is the jax twin).
    `PYCIRCUIT_HDL_BACKEND=c` / `numpy` (or `hdl.set_backend`, or a class
    attribute `hdl_backend`) pins it, applied at class creation.
    `cls._hdl_backend_status` always says what actually happened: `'c'`,
    `'numpy'`, `'numpy (eager path: ...)'`, `'numpy (compile failed:
    ...)'` (a missing compiler too, when C is asked for), `'numpy (cffi
    not installed)'`, `'numpy (auto: <why>)'`, `'auto (resolved at the
    first instance)'`.  When C was REQUESTED and cannot be delivered, a
    warning is issued as well -- requesting C and silently getting numpy
    is the failure this refuses; under 'auto' only a build that FAILS with
    a compiler present warns (a defect, either way).

The limiter
    A class that declares `$limit` gets its generated `limit()` as a
    kernel of its own (`_hdl_climit`): rendered from the spec and the
    parameter functions' kept statements, built and stored beside the
    chain functions' objects, bound by `_resolve` after them and unbound
    by `detach`.  It refuses itself -- an unprintable node, `tanh`, a
    parameter without ingredients -- and never the class;
    `info['_c_limit_status']` says.  On a C-bound class the closure's
    parameter functions are not called: a test that counts them switches
    the limiter kernel off (`_hdl_climit.ENABLED`).  In a circuit's `limit`
    pass every run of such elements is walked in ONE C call on the live
    state (`_hdl_climit.limit_walk`), the loop's order and write-backs
    kept.

Fidelity
    The `.so` is compiled with `-O2` and no fast-math, no FMA
    contraction, and every transcendental kept a real libm call
    (`-fno-builtin-*`), because gcc otherwise folds `pow(x, 2.0)` to
    `x*x` -- more accurate than numpy's scalar `**2`, which is glibc
    `pow`, and therefore WRONG for a backend whose contract is
    bit-identity with the numpy path.  The kernel prelude declares
    those functions `__attribute__((const))` (`hdl._KERNEL_C`), so gcc
    calls an identical one once instead of at every repeat -- still
    glibc's call, still its value.  Two named exceptions
    (`test_hdl_cbackend.py` pins both):

    * `tanh`: numpy ships its own vectorised tanh, which differs from
      libm's by an ulp on some arguments.  A chain that uses `tanh`
      agrees with the numpy path to that ulp (amplified where the
      model cancels), and agrees BITWISE with the same numpy source
      run with libm's tanh.
    * the sign of exact zeros: the numpy path computes integer-typed
      subchains (`where(c, 1, 0)` is int64, and integer zeros carry no
      sign) where C is all doubles.  Values compare equal; bytes can
      differ on zero entries.

    And one that depends on the CPU: numpy picks its loops by the CPU it
    runs on, and where it brings its own float64 `exp`, `log`, ...
    (AVX-512 builds, for one) they are not the C library's.  `libm_check`
    compares each function of the kernel prelude with numpy's, in the
    forms the numpy path calls it, once a process on a fixed probe set.
    A class whose kernels or limiter call a function that differs KEEPS C
    (the speed, with the difference named) and names it:
    `cls._hdl_backend_ulp`, `explain`'s backend line, and one
    `PlatformWarning` a process for the functions beyond tanh -- the class
    then agrees with the numpy path to that function's ulp (amplified
    where the model cancels), and bitwise with its numpy source run with
    the C library's functions, as a tanh class does.

Cache
    The `.so` lives beside the pickled compile-cache entries
    (`_hdl_cache.cache_dir()`), keyed by SHA-256 of the complete C
    source, the compiler's version line and the flags.  Cold: one `cc`
    run per distinct function, launched in parallel across the class's
    functions.  Warm: `dlopen` only.  A corrupt or unloadable `.so` is
    rebuilt in place.  Writes are `mkstemp` + `os.replace`, so two
    processes on one key cannot tear the file.  With the compile cache
    disabled (`PYCIRCUIT_HDL_CACHE=0`) objects are built in a
    per-process temporary directory instead and discarded at exit.

Calling convention
    `void hdl_fn(const double *x, const double *p, double *out)` --
    `x` the solution vector, `p` the packed trailing arguments
    (parameter values, temperature, givenness flags) in the numpy
    signature's order, `out` the preallocated result.  `p` is packed
    once per parameter change (`update_iparv` drops it via the
    element's observer), not per call; the temperature slot is written
    per call because `epar` carries it.  A fresh `out` is allocated per
    call so the returned array is independently owned, exactly as
    `np.asarray(list)` made it on the numpy path.  cffi releases the
    GIL for the call and the function has no globals, so it is
    reentrant.

The batches
    In a circuit's stamp plan a C-bound class's elements run as ONE call
    per class per pass (`_hdl_batch`): a driver of its own (`csrc/pass.c`,
    a runtime unit loaded through `_cbuild.load_unit`) loops over the
    elements calling the kernel's function pointer on each element's
    `x[nm]` and its own pack, the temperature written into the pack's
    slot exactly as `CKernel.__call__` writes it.  What a batch cannot
    serve -- an instance shadow, a patched method, a detached class, a
    temperature that is not one number -- takes `CKernel.__call__` per
    element as before.

The fused kernel
    For the transient's evaluate core a class also gets ONE kernel for
    its four passes (`fuse_csrc`): the printed statements unioned by
    name, each run under the bits of the passes that print it, each
    pass's outputs to its own array -- the same statements, so the same
    values (a NaN's sign aside).  Built with the others, bound as
    `info['_c_fused']`, refused on its own (`info['_c_fused_status']`);
    `PYCIRCUIT_HDL_CFUSE=0` leaves it unbuilt and unused.
"""

import contextlib
import importlib.util
import os
import re
import subprocess
import tempfile
import types

import numpy as np

from pycircuit.circuit import _cbuild
from pycircuit.circuit import simwarnings as _sw

## (the build and load machinery lives in `_cbuild`, re-exported for this
## module's callers -- PATCH THEM ON `_cbuild`: `CFLAGS`, `_compiler` and
## `_loaded` are read there)
from pycircuit.circuit._cbuild import (  # noqa: F401
    _FALSE,
    CFLAGS,
    SANITIZE,
    SANITIZE_FLAGS,
    SO_FORMAT,
    CompileError,
    _asan_preloaded,
    _build,
    _cc_tag,
    _compiler_env,
    _dlopen,
    _fn_cdef,
    _key_lock,
    _so_dir,
    find_compiler,
    find_so,
    load_kernel,
    so_path,
    source_key,
)

## ----------------------------------------------------------------------
## Kernels.

#: numpy's float64 dtype object, one per process: `x.dtype is _F64` is a
#: pointer test where `x.dtype == np.float64` is a comparison.
_F64 = np.dtype('float64')


class CKernel(object):
    """One compiled chain function, callable as the element methods
    need it: `kernel(element, x, epar)` -> a fresh float64 array.

    THE CALL IS THE COST for a mid-sized model, so it is kept small:
    `ffi.from_buffer` on the arrays (the buffer protocol; the cdata keeps
    the array alive for the call), a pre-resolved pointer type (`dptr`), a
    type test for the usual float `T`.  A fresh `out` every call (callers
    keep and even write into what they get); None where the numpy function
    must answer; `ValueError` for a state shorter than the kernel reads."""

    __slots__ = ('ffi', 'cfn', 'shape', 'nx', 'n_p', 't_index', 'key',
                 'built_s', 'dptr')

    def __init__(self, ffi, cfn, shape, nx, layout, key, built_s):
        self.ffi, self.cfn, self.shape = ffi, cfn, shape
        self.dptr = ffi.typeof('double *')
        self.nx = nx
        self.n_p, self.t_index = layout
        self.key, self.built_s = key, built_s

    def pack(self, element):
        """The packed trailing-argument vector for `element`, built
        once per parameter state (the iparv observer drops it)."""
        from pycircuit.circuit import hdl
        vals = hdl._params_of(element)
        n_par = len(element._hdl_paramnames)
        p = np.empty(self.n_p)
        p[:n_par] = vals[:n_par]
        if self.t_index is not None:
            p[self.t_index] = 0.0
        p[n_par + (1 if self.t_index is not None else 0):] = vals[n_par:]
        cast = self.ffi.cast('double *', p.ctypes.data)
        return (p, cast)

    def __call__(self, element, x, epar):
        """The kernel's output at `x`, or None where it cannot serve the
        call -- a parameter or temperature that is not one number, a state
        that is not one vector of numbers -- and the numpy function then
        does (`hdl._chained_eval`), as it broadcasts there."""
        d = element.__dict__
        packed = d.get('_hdl_cp')
        if packed is None:
            try:
                packed = self.pack(element)
            except (TypeError, ValueError):
                packed = False
            d['_hdl_cp'] = packed
        if packed is False:
            return None
        p, pcast = packed
        t_index = self.t_index
        if t_index is not None:
            T = getattr(epar, 'T', 300.0)
            ## (a float or an int is one number; anything else -- a numpy
            ## scalar, a 0-d array -- is asked, an array steps aside)
            if type(T) is not float and type(T) is not int \
                    and np.ndim(T) != 0:
                return None
            p[t_index] = T
        if (type(x) is np.ndarray and x.dtype is _F64
                and x.flags.c_contiguous):
            xa = x
        else:
            try:
                xa = np.ascontiguousarray(x, dtype=float)
            except (TypeError, ValueError):
                return None
        if xa.ndim != 1:
            return None
        if xa.shape[0] < self.nx:
            raise ValueError('x has shape %r; this kernel reads x[0..%d]'
                             % (np.shape(x), self.nx - 1))
        out = np.empty(self.shape)
        ffi = self.ffi
        self.cfn(ffi.from_buffer(self.dptr, xa), pcast,
                 ffi.from_buffer(self.dptr, out))
        return out


def kernel_for(fn, nx, rebuild_corrupt=True):
    """The `CKernel` for a chain function carrying `_csrc`: `load_kernel`
    on its source, wrapped.  Returns `(kernel, cold)`."""
    ffi, cfn, key, cold, secs = _cbuild.load_kernel(fn._csrc, None, rebuild_corrupt)
    kern = CKernel(ffi, cfn, tuple(fn._cshape), nx, tuple(fn._clayout),
                   key, secs)
    return kern, cold


## ----------------------------------------------------------------------
## The fused kernel.

#: ONE MASKED KERNEL A CLASS for the evaluate core (`_tran_core`): a
#: class's four passes print four functions that each recompute their
#: whole chain -- `i` is inside `G`, `q` inside `C`, and `G` and `C` share
#: about half their statements -- where one function evaluating what a
#: call wants computes each statement once.  Printed from the passes'
#: C (`fuse_csrc`), built beside them and bound as `info['_c_fused']`;
#: refused on its own, never the class.  env `PYCIRCUIT_HDL_CFUSE=0`: not
#: built, and the core calls the passes' own kernels (read at the bind and
#: at every readiness check).
FUSE = os.environ.get('PYCIRCUIT_HDL_CFUSE', '1').strip().lower() not in _FALSE

#: the passes a fused kernel computes: bit `1 << k` is `FUSE_PASSES[k]`
#: (the evaluate core's order)
FUSE_PASSES = ('G', 'C', 'i', 'q')

_FUSE_STMT = re.compile(r'  const double L_(\w+) = ([^;]*);\Z')
_FUSE_OUT = re.compile(r'  out\[(\d+)\] = ([^;]*);\Z')
_FUSE_REF = re.compile(r'\bL_(\w+)')


class FuseRefused(Exception):
    """The passes' printed C cannot be unioned into one masked kernel."""


def fused_cdef():
    """The fused kernel's entry: `void hdl_fn(x, p, want, oG, oC, oi, oq)`."""
    from pycircuit.circuit import hdl
    return (f'void {hdl._C_ENTRY}(const double *x, const double *p, long want, '
            'double *oG, double *oC, double *oi, double *oq);')


def fuse_csrc(srcs):
    """One masked C function from a class's printed pass kernels (`srcs`:
    a pass of `FUSE_PASSES` -> its `_csrc`): the union of their
    statements by name, each function's order kept (every statement
    follows what it reads), each statement carrying the bits of the passes
    that print it and run under `if (want & bits)`; each pass's outputs
    written to its own array (`oG`, `oC`, `oi`, `oq`) under its bit.  The
    statements are the printed ones, unchanged: the same IEEE operations,
    only grouped (a NaN's sign may differ where gcc orders a commutative
    operation's operands otherwise -- the caveat the backend carries
    against numpy).  `FuseRefused` for a line the printer does not make
    (`hdl._render_c` is the only printer), one name standing for two
    texts, or a statement or output reading a name outside its passes.

    Returns `(text, stats)`."""
    from pycircuit.circuit import hdl
    head = f'void {hdl._C_ENTRY}(const double *x, const double *p, double *out) {{'
    order, text, bits, outs = [], {}, {}, {}
    for b, m in enumerate(FUSE_PASSES):
        src = srcs.get(m)
        if src is None:
            continue
        lines = src.split('\n')
        if len(lines) < 3 or lines[0] != head or lines[-2:] != ['}', '']:
            raise FuseRefused(f'{m}: not the printed frame')
        o = []
        for ln in lines[1:-2]:
            mt = _FUSE_STMT.match(ln)
            if mt is not None and not o:
                name, expr = mt.groups()
                prev = text.get(name)
                if prev is None:
                    order.append(name)
                    text[name], bits[name] = expr, 1 << b
                elif prev != expr:
                    raise FuseRefused(f'L_{name} stands for two texts')
                else:
                    bits[name] |= 1 << b
                continue
            mt = _FUSE_OUT.match(ln)
            if mt is None or int(mt.group(1)) != len(o):
                raise FuseRefused(f'{m}: a line the printer does not make')
            o.append(mt.group(2))
        outs[m] = o
    if len(outs) < 2:
        raise FuseRefused('fewer than two passes')
    ## every name read: printed before its reader and in each of its
    ## reader's passes (so a run its guard skips is read by no run it lets
    ## through)
    pos = {n: k for k, n in enumerate(order)}
    for n in order:
        for d in _FUSE_REF.findall(text[n]):
            if d not in pos or pos[d] >= pos[n] or bits[n] & ~bits[d]:
                raise FuseRefused(f'L_{n} reads L_{d} outside its passes')
    for b, m in enumerate(FUSE_PASSES):
        for e in outs.get(m, ()):
            for d in _FUSE_REF.findall(e):
                if d not in pos or not bits[d] & (1 << b):
                    raise FuseRefused(f'an output of {m} reads L_{d} outside it')
    lines = [(f'void {hdl._C_ENTRY}(const double *x, const double *p, long want, '
              'double *oG, double *oC, double *oi, double *oq) {')]
    for k in range(0, len(order), 32):
        lines.append('  double {};'.format(', '.join('L_' + n for n in order[k:k + 32])))
    run, runs = None, 0
    for n in order:
        if bits[n] != run:
            if run is not None:
                lines.append('  }')
            run, runs = bits[n], runs + 1
            lines.append(f'  if (want & {run}) {{')
        lines.append(f'    L_{n} = {text[n]};')
    if run is not None:
        lines.append('  }')
    for b, m in enumerate(FUSE_PASSES):
        if m in outs:
            lines.append(f'  if (want & {1 << b}) {{')
            lines.extend(f'    o{m}[{k}] = {e};' for k, e in enumerate(outs[m]))
            lines.append('  }')
    lines.append('}')
    stats = {'union': len(order), 'runs': runs,
             'printed': sum(1 for m in outs for ln in srcs[m].split('\n')
                            if _FUSE_STMT.match(ln))}
    return '\n'.join(lines) + '\n', stats


class FusedKernel:
    """A class's fused kernel: the loaded function, the pass kernels it was
    printed from (`parts`: the core calls it only while they are the
    class's bound kernels), each pass's output length (`sizes`), the bits
    of its passes (`mask`), the x length it reads and the pack layout."""

    __slots__ = ('built_s', 'cfn', 'ffi', 'key', 'layout', 'mask', 'nx', 'parts', 'sizes',
                 'stats')

    def __init__(self, ffi, cfn, parts, sizes, nx, layout, key, built_s, stats):
        self.ffi, self.cfn, self.parts, self.sizes = ffi, cfn, parts, sizes
        self.mask = sum(1 << FUSE_PASSES.index(m) for m in parts)
        self.nx, self.layout, self.key, self.built_s, self.stats = nx, layout, key, built_s, stats


def fuse_source(info):
    """The class's fused source as `_build_missing_parallel` takes it, or
    None with `info['_c_fused_status']` saying why: the passes' kernels
    must agree on the x they read and the pack they take."""
    if not FUSE:
        info['_c_fused_status'] = 'off'
        return None
    srcs, shapes, nxs, layouts = {}, {}, set(), set()
    for m in FUSE_PASSES:
        fn = info['funcs'].get(m)
        if fn is not None and getattr(fn, '_csrc', None) is not None:
            srcs[m], shapes[m] = fn._csrc, tuple(fn._cshape)
            nxs.add(fn._cshape[0])
            layouts.add(tuple(fn._clayout))
    if len(nxs) > 1 or len(layouts) > 1:
        info['_c_fused_status'] = 'refused (the passes read different x or packs)'
        return None
    try:
        text, stats = fuse_csrc(srcs)
    except FuseRefused as e:
        info['_c_fused_status'] = f'refused ({e})'
        return None
    return types.SimpleNamespace(_csrc=text, shapes=shapes, nx=nxs.pop(),
                                 layout=layouts.pop(), stats=stats)


def bind_fused(info, src):
    """Bind the class's fused kernel from `fuse_source`'s object (None:
    nothing to bind) once its passes' kernels are bound: `info['_c_fused']`
    the kernel or None, `info['_c_fused_status']` 'c' or why not."""
    info['_c_fused'] = None
    if src is None:
        return None
    parts = {m: info['funcs'][m].__dict__.get('_hdl_c') for m in src.shapes}
    if not all(isinstance(k, CKernel) for k in parts.values()):
        info['_c_fused_status'] = 'refused (a pass without its kernel)'
        return None
    try:
        ffi, cfn, key, _cold, secs = _cbuild.load_kernel(src._csrc, fused_cdef())
    except (_cbuild.CompileError, OSError) as e:
        info['_c_fused_status'] = f'compile failed ({e})'
        return None
    sizes = {m: int(np.prod(s)) for m, s in src.shapes.items()}
    fk = FusedKernel(ffi, cfn, parts, sizes, src.nx, src.layout, key, secs, src.stats)
    info['_c_fused'] = fk
    info['_c_fused_status'] = 'c'
    return fk


def _build_missing_parallel(fns):
    """Compile, in parallel, every distinct `_csrc` in `fns` whose
    object is not on disk yet.  Failures fall through to `kernel_for`,
    which reports them; this only warms the store."""
    from pycircuit.circuit import hdl
    cc, _version = _cbuild.find_compiler()
    if cc is None:
        return
    jobs = []
    seen = set()
    with contextlib.ExitStack() as held:
        for fn in fns:
            key = _cbuild.source_key(fn._csrc)
            if key in seen or key in _cbuild._loaded:
                continue
            seen.add(key)
            if _cbuild.find_so(key) is not None:
                continue
            ## a key another process is building is left to it:
            ## `kernel_for` waits on the lock and then finds its object
            if not held.enter_context(_cbuild._key_lock(key, blocking=False)):
                continue
            if _cbuild.find_so(key) is not None:
                continue
            dest = _cbuild.so_path(key)
            d = os.path.dirname(dest)
            fd, cpath = tempfile.mkstemp(suffix='.c', dir=d)
            with os.fdopen(fd, 'w') as fh:
                fh.write(hdl._KERNEL_C)
                fh.write(fn._csrc)
            sopath = cpath[:-2] + '.so.tmp'
            proc = subprocess.Popen([cc] + list(_cbuild.CFLAGS) +
                                    ['-o', sopath, cpath, '-lm'],
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL, env=_cbuild._compiler_env())
            jobs.append((proc, cpath, sopath, dest))
        for proc, cpath, sopath, dest in jobs:
            rc = proc.wait()
            if rc == 0:
                try:
                    os.replace(sopath, dest)
                except OSError:
                    pass
            for p in (cpath, sopath):
                try:
                    os.unlink(p)
                except OSError:
                    pass


## ----------------------------------------------------------------------
## Attaching to a class.

#: The functions worth a C body: the per-iteration ones.  `CY` (AC),
#: `u`/`dudt`/`uac` (x-free, cheap) stay numpy.  (The limiter is a kernel
#: of its own, `_hdl_climit`, bound beside these.)
C_FUNCS = ('i', 'G', 'q', 'C', 'i_dc', 'G_dc')


## -- THE C LIBRARY AGAINST NUMPY -----------------------------------------------

#: The transcendentals the kernel prelude declares (`hdl._KERNEL_C`), each
#: with the numpy function the numpy path calls for it (`hdl._C_NAMES`).
#: `pow` is checked as the `**` of a 0-d array (`numpy.where`'s result:
#: numpy's own `power` loop) and of a numpy scalar (the C library's `pow`
#: on every CPU, checked anyway).  `sqrt`, `fabs`, `floor` and `ceil` are
#: exact in every library: not here.
LIBM_PAIRS = (('exp', 'exp'), ('log', 'log'), ('log1p', 'log1p'), ('expm1', 'expm1'),
              ('sin', 'sin'), ('cos', 'cos'), ('tan', 'tan'), ('asin', 'arcsin'),
              ('acos', 'arccos'), ('atan', 'arctan'), ('sinh', 'sinh'), ('cosh', 'cosh'),
              ('tanh', 'tanh'), ('atan2', 'arctan2'), ('hypot', 'hypot'), ('pow', 'power'))
_LIBM_BINARY = frozenset(('atan2', 'hypot', 'pow'))
_LIBM_CALL = re.compile(r'\b(' + '|'.join(c for c, _n in LIBM_PAIRS) + r')\s*\(')
#: the probe's outcome: `differ` None until `libm_check` ran, then the
#: functions that differ; `detail` per function (differing, probed);
#: `status` what was found, or why nothing could be checked (`libm_status`).
#: A dict, not module scalars: the leak detector reports a module scalar
#: that moves, and this fills at a process's first C resolve.
_LIBM = {'differ': None, 'detail': {}, 'status': 'not checked yet'}


def _libm_probe_values():
    """The fixed probe set: the special values (signed zeros, subnormals,
    the smallest normal, the overflow and underflow edges of `exp`, +-1e300,
    the largest double, infinities, NaN), magnitudes from 1e-310 to 1e308
    of both signs, and dense bands where the functions curve (+-40, +-1)."""
    import math
    specials = [0.0, -0.0, 5e-324, -5e-324, 2.2250738585072014e-308,
                -2.2250738585072014e-308, 1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 3.0, 10.0,
                -10.0, 100.0, 709.78, 709.79, 710.0, -745.1, -745.2, -746.0, 1e300, -1e300,
                1.7976931348623157e308, -1.7976931348623157e308, math.inf, -math.inf,
                math.nan, math.pi, math.pi / 2, 1e-8, -1e-8, 0.999999, -0.999999, 1e22,
                2.0 ** 52, 2.0 ** 60]
    mags = list(10.0 ** np.linspace(-310.0, 308.0, 220))
    bands = list(np.linspace(-40.0, 40.0, 161)) + list(np.linspace(-1.0, 1.0, 81))
    return np.array(specials + mags + [-m for m in mags] + bands, dtype=float)


def _libm_differ(a, b):
    """Entries whose bits differ, NaN against NaN not counted (a NaN's
    payload and sign are not the contract; the tests' `_compare`)."""
    a = np.ascontiguousarray(a, dtype=float)
    b = np.ascontiguousarray(b, dtype=float)
    return (~(np.isnan(a) & np.isnan(b))) & (a.view(np.int64) != b.view(np.int64))


def _libm_probe(numpy_funcs=None):
    """`(the C names whose numpy function differs, {name: (differing,
    probed)})` -- or None where the C library cannot be reached from here.
    Each function along an array of the probe set and on every third value
    as a numpy scalar (a CPU's numpy may run the two through different
    loops); `numpy_funcs` (tests: a planted difference) replaces numpy's
    functions by their numpy name."""
    import ctypes
    import ctypes.util
    try:
        lib = ctypes.CDLL(ctypes.util.find_library('m') or None)
        cfun = {}
        for c, _n in LIBM_PAIRS:
            f = getattr(lib, c)
            f.restype = ctypes.c_double
            f.argtypes = [ctypes.c_double] * (2 if c in _LIBM_BINARY else 1)
            cfun[c] = f
    except (OSError, AttributeError) as e:
        _LIBM['status'] = f'not checked: the C library is not reachable ({e})'
        return None
    over = numpy_funcs or {}
    vals = _libm_probe_values()
    third = vals[::3]
    a2 = vals[::12]
    b2 = np.concatenate([vals[::29], [2.0, 3.0, 0.5, -1.0, -0.5, 1.5, -2.0]])
    A, B = (m.ravel() for m in np.meshgrid(a2, b2))
    differ, detail = set(), {}
    with np.errstate(all='ignore'):
        for c, n in LIBM_PAIRS:
            uf = over.get(n) or getattr(np, n)
            f = cfun[c]
            if c in _LIBM_BINARY:
                ref = np.array([f(float(a), float(b)) for a, b in zip(A, B)])
                got = [_libm_differ(uf(A, B), ref)]
                if c == 'pow':
                    k = slice(None, None, 4)
                    got.append(_libm_differ([uf(np.asarray(a), b) for a, b in zip(A[k], B[k])],
                                            ref[k]))
                    got.append(_libm_differ([np.float64(a) ** np.float64(b)
                                             for a, b in zip(A[k], B[k])], ref[k]))
            else:
                ref = np.array([f(float(v)) for v in vals])
                got = [_libm_differ(uf(vals), ref),
                       _libm_differ([uf(np.float64(v)) for v in third], ref[::3])]
            d, m = sum(int(g.sum()) for g in got), sum(g.size for g in got)
            detail[c] = (d, m)
            if d:
                differ.add(c)
    _LIBM['status'] = ('checked: ' + (', '.join(f'{c} differs in {detail[c][0]} of {detail[c][1]}'
                                                for c in sorted(differ))
                                      or 'every function agrees'))
    return frozenset(differ), detail


def libm_check():
    """The functions of the kernel prelude whose numpy twin gives other
    bits than the C library's on this CPU, as a frozenset of C names;
    measured once a process (`_libm_probe`).  Empty where nothing
    could be checked (`libm_status` says so).  Where it finds one beyond
    tanh, the process gets ONE `PlatformWarning` naming them."""
    if _LIBM['differ'] is None:
        out = _libm_probe()
        _LIBM['differ'], _LIBM['detail'] = (frozenset(), {}) if out is None else out
        beyond = sorted(_LIBM['differ'] - {'tanh'})
        if beyond:
            _sw.warn(f"hdl C backend: on this CPU numpy's {', '.join(beyond)} differ from the "
                     "C library's: a C-bound class calling them agrees with the numpy path to "
                     "an ulp, not bitwise, as tanh does on every CPU -- `cls._hdl_backend_ulp` "
                     f"names a class's ({_LIBM['status']})", _sw.PlatformWarning)
    return _LIBM['differ']


def libm_status():
    """What `libm_check` found ('checked: ...'), or why it found nothing."""
    return _LIBM['status']


def _limit_helper_calls():
    """For each helper of the limiter's own prelude (`_hdl_climit._LIMIT_C`):
    the transcendentals it calls, through the helpers it calls."""
    from pycircuit.circuit import _hdl_climit
    got = _LIBM.get('helpers')
    if got is None:
        pre = _hdl_climit._LIMIT_C
        heads = list(re.finditer(r'^static\s+\w+\s+(\w+)\s*\(', pre, re.MULTILINE))
        bodies = {m.group(1): pre[m.end():(heads[k + 1].start() if k + 1 < len(heads)
                                            else len(pre))]
                  for k, m in enumerate(heads)}
        direct = {h: set(_LIBM_CALL.findall(b)) for h, b in bodies.items()}
        got = {}
        for h in bodies:
            seen, todo, calls = {h}, [h], set()
            while todo:
                g = todo.pop()
                calls |= direct[g]
                for other in bodies:
                    if other not in seen and re.search(r'\b' + other + r'\s*\(', bodies[g]):
                        seen.add(other)
                        todo.append(other)
            got[h] = frozenset(calls)
        _LIBM['helpers'] = got
    return got


def libm_calls(csrc):
    """The functions of `LIBM_PAIRS` a kernel's C source calls -- a
    limiter's helper prelude counted only through the helpers its own code
    calls (every limit kernel carries the whole prelude)."""
    from pycircuit.circuit import _hdl_climit
    pre = _hdl_climit._LIMIT_C
    if pre not in csrc:
        return frozenset(_LIBM_CALL.findall(csrc))
    body = csrc.replace(pre, '')
    out = set(_LIBM_CALL.findall(body))
    for h, calls in _limit_helper_calls().items():
        if re.search(r'\b' + h + r'\s*\(', body):
            out |= calls
    return frozenset(out)


def _name_ulp(cls, info, fns, lsrc):
    """`cls._hdl_backend_ulp`: the listed functions (`libm_check`) its bound
    C calls -- the chain functions and, where it bound, its limiter (`lsrc`,
    `_hdl_climit.source_for`'s; the fused kernel's statements are the chain
    functions')."""
    srcs = [fn._csrc for fn in fns]
    if lsrc is not None and info.get('_c_limit') is not None:
        srcs.append(lsrc._csrc)
    listed = libm_check()
    called = frozenset().union(*(libm_calls(s) for s in srcs)) if srcs else frozenset()
    cls._hdl_backend_ulp = tuple(sorted(called & listed))


def _note(cls, status, warn=None):
    cls._hdl_backend_status = status
    ## (named again where the class binds C: `_name_ulp`)
    cls._hdl_backend_ulp = ()
    if warn:
        _sw.warn(f'hdl C backend: {cls.__name__}: {warn}', _sw.CostWarning)


def detach(cls, info):
    """Remove any C bindings so the class runs numpy again."""
    for name in C_FUNCS:
        fn = info['funcs'].get(name)
        if fn is not None and hasattr(fn, '_hdl_c'):
            del fn._hdl_c
    info['_c_bound'] = False
    info['_c_fused'] = None
    info['_c_fused_status'] = 'detached'
    if info.get('limit_spec'):
        from pycircuit.circuit import _hdl_climit
        _hdl_climit.unbind(info)


def attach(cls, info):
    """Bind C kernels to `cls`'s chained functions, per the requested
    backend -- or, for 'auto' (the default), mark the class to be
    resolved at its first instance (`ensure`).  Always sets
    `cls._hdl_backend_status`; never raises for a build problem (it warns
    and falls back to numpy), only for a misconfigured backend name."""
    from pycircuit.circuit import hdl
    which = hdl._backend_requested(cls)
    if which == 'numpy':
        detach(cls, info)
        info['_backend_pending'] = False
        _note(cls, 'numpy')
        return
    if not info['chained']:
        detach(cls, info)
        info['_backend_pending'] = False
        _note(cls, 'numpy (eager path: the C backend serves let-chain '
                   'models only)')
        return
    if which == 'auto':
        detach(cls, info)
        if info.get('_backend_seen'):
            ## instances exist (a pin was lifted): resolve now, so they run
            ## what they ran before it
            info['_backend_pending'] = False
            _resolve(cls, info, explicit=False)
        else:
            info['_backend_pending'] = True
            _note(cls, 'auto (resolved at the first instance)')
        return
    info['_backend_pending'] = False
    _resolve(cls, info, explicit=True)


def ensure(cls, info, toolkit=None):
    """At an instance's construction (`Behavioural.__init__`): resolve an
    'auto' class's backend, once.  An instance on the JAX or symbolic
    toolkit leaves it unresolved -- it runs its own path (the jax twin)."""
    if getattr(toolkit, 'jax', False) or getattr(toolkit, 'symbolic', False):
        return
    info['_backend_seen'] = True
    if info.get('_backend_pending'):
        info['_backend_pending'] = False
        _resolve(cls, info, explicit=False)


def _resolve(cls, info, explicit):
    """Bind the C kernels for an explicit request or an 'auto' class.
    Where C cannot be served: numpy, with the reason in the status -- and a
    warning when C was asked for; under 'auto' only a build that fails
    with a compiler present warns."""
    from pycircuit.circuit import _paths

    def refuse(status, why):
        _paths.COUNTS['once:backend:numpy'] += 1
        detach(cls, info)
        if explicit:
            _note(cls, f'numpy ({status})',
                  warn=f'requested backend "c" but {why}; running numpy')
        else:
            _note(cls, f'numpy (auto: {status})')

    funcs = info['funcs']
    todo = {}
    for name in C_FUNCS:
        fn = funcs.get(name)
        if fn is None or id(fn) in todo:
            continue
        if getattr(fn, '_csrc', None) is not None:
            todo[id(fn)] = fn
        elif hasattr(fn, '_creason'):
            return refuse(f'C rendering unavailable: {fn._creason}',
                          'the chain could not be printed to C '
                          f'({fn._creason})')
    if not todo:
        return refuse('no C source on this class; recompile with '
                      'hdl.EMIT_C_SOURCE on',
                      'the compiled functions carry no C source')
    ## The kernels load through cffi, which numpy does not need: without
    ## it the class is refused here (numpy), never an ImportError out of
    ## class creation.
    if importlib.util.find_spec('cffi') is None:
        return refuse('cffi not installed', 'cffi is not installed')
    if not explicit:
        from pycircuit.circuit import _hdl_cache
        if not _hdl_cache.enabled():
            ## no store: every process would build every kernel again
            return refuse('the compile cache is off', '')
        try:
            _cbuild._so_dir()
        except OSError:
            return refuse('the compile cache directory is not usable', '')
        if _cbuild.find_compiler()[0] is None and any(
                _cbuild.find_so(_cbuild.source_key(fn._csrc)) is None
                for fn in todo.values()):
            return refuse('no C compiler', '')
    ## (An explicit request has no compiler check: a warm `.so` store
    ## serves a machine with no compiler at all -- `kernel_for` only builds
    ## on a miss, and raises `CompileError` naming the missing compiler
    ## when it must build.  'auto' looked above, and stays quiet.)
    cold = False
    ## (the class's limiter, printed from its spec and built beside the
    ## chain functions; refused on its own, never the class: `_hdl_climit`)
    from pycircuit.circuit import _hdl_climit
    lsrc = _hdl_climit.source_for(info) if info.get('limit_spec') else None
    ## (and the evaluate core's fused kernel, printed from the passes' C:
    ## refused on its own too, `fuse_source`)
    fsrc = fuse_source(info)
    try:
        _build_missing_parallel(list(todo.values())
                                + [s for s in (lsrc, fsrc) if s is not None])
        for fn in todo.values():
            ## The x length the kernel reads equals its output row
            ## count: `i`, `G`, `q`, `C` are all n-per-side.  Taken from
            ## the function itself because `cls._hdl_info` is not
            ## assigned yet when the metaclass calls attach.
            kern, was_cold = kernel_for(fn, fn._cshape[0])
            cold = cold or was_cold
            fn._hdl_c = kern
    except (_cbuild.CompileError, OSError) as e:
        ## (an OSError too: the store's directory -- a file in its place, no
        ## permission -- refuses the class here, never out of class creation)
        _paths.COUNTS['once:backend:compile_failed'] += 1
        detach(cls, info)
        _note(cls, 'numpy (compile failed: %s)' % e,
              warn=('requested backend "c" but the build failed (%s); '
                    'running numpy' if explicit else
                    'the C build failed (%s); running numpy') % e)
        return
    ## (read by `_hdl_cse.take`: a C-bound class's calls never fuse)
    info['_c_bound'] = True
    _paths.COUNTS['once:backend:c'] += 1
    _note(cls, 'c')
    if lsrc is not None:
        _hdl_climit.bind(cls, info, lsrc)
    bind_fused(info, fsrc)
    _name_ulp(cls, info, todo.values(), lsrc)
