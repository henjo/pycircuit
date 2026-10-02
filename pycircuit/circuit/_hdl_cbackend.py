# -*- coding: utf-8 -*-
# Copyright (c) 2008-2026 Pycircuit Development Team
# See LICENSE for details.
"""The C backend for chained `Behavioural` elements.

`generate_code` prints every x-taking let-chain twice: the numpy source
(`fn._src`, the reference -- what `explain()` shows, what the compile
cache stores, what every test compares against) and a C rendering
(`fn._csrc`, from `hdl._CChainPrinter`).  This module owns everything
after that print: finding a compiler, building the shared object,
caching it, loading it, and calling it.

Why it exists: the chained numpy function is straight-line scalar
Python -- for the PSP MOSFET's Jacobian, ~21 400 interpreted calls per
evaluation.  The same statements compiled by the system C compiler run
219-229x faster at bitwise-identical results (`doc/
backend_spike_260826.md`); this module is that measurement built out.

Selection
    `PYCIRCUIT_HDL_BACKEND=c` (or `hdl.set_backend('c')`, or a class
    attribute `hdl_backend = 'c'`), read once per class at compile /
    attach time.  The default is numpy.  The choice is per CLASS, and
    `cls._hdl_backend_status` always says what actually happened:
    `'c'`, `'numpy'`, `'numpy (eager path: ...)'`, `'numpy (no C
    compiler ...)'`, `'numpy (compile failed: ...)'`.  When C was
    REQUESTED and cannot be delivered, a warning is issued as well --
    requesting a 200x backend and silently getting the 1x one is the
    failure mode this refuses to have.

Fidelity
    The `.so` is compiled with `-O2` and no fast-math, no FMA
    contraction, and every transcendental kept a real libm call
    (`-fno-builtin-*`), because gcc otherwise folds `pow(x, 2.0)` to
    `x*x` -- more accurate than numpy's scalar `**2`, which is glibc
    `pow`, and therefore WRONG for a backend whose contract is
    bit-identity with the numpy path.  The kernel prelude declares
    those functions `__attribute__((const))` (`hdl._KERNEL_C`), so gcc
    calls an identical one once instead of at every repeat -- still
    glibc's call, still its value; PSP's `G` 50 -> 22 us (2026-10-02).
    Two named, measured exceptions (`test_hdl_cbackend.py` pins both):

    * `tanh`: numpy ships its own vectorised tanh, which differs from
      libm's in ~30% of arguments by an ulp.  A chain that uses `tanh`
      agrees with the numpy path to that ulp (amplified where the
      model cancels), and agrees BITWISE with the same numpy source
      run with libm's tanh.
    * the sign of exact zeros: the numpy path computes integer-typed
      subchains (`where(c, 1, 0)` is int64, and integer zeros carry no
      sign) where C is all doubles.  Values compare equal; bytes can
      differ on zero entries.

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
"""

import contextlib
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from pycircuit.circuit import simwarnings as _sw

#: Bump when the key must invalidate every built object (a flags change
#: already does; this is for changes to the convention itself).
SO_FORMAT = 1

#: Compilation flags.  Strict IEEE doubles: no fast-math, no FMA
#: contraction, no `-march`.  Every `-fno-builtin-<f>` keeps `<f>` a
#: real libm call so the compiler cannot substitute its own (usually
#: MORE accurate) folding for glibc's runtime result -- `pow(x, 2.0)`
#: is the measured case: glibc is within 0.52 ulp, `x*x` is exact, and
#: the numpy path computes glibc's answer.  Constant folding of pure
#: literals is fine (gcc folds via MPFR, correctly rounded, same as
#: Python's own parse-time arithmetic); folding of CALLS is not.
CFLAGS = ('-O2', '-shared', '-fPIC', '-fno-fast-math', '-ffp-contract=off',
          '-fno-builtin-pow', '-fno-builtin-exp', '-fno-builtin-log',
          '-fno-builtin-log1p', '-fno-builtin-expm1', '-fno-builtin-sin',
          '-fno-builtin-cos', '-fno-builtin-tan', '-fno-builtin-asin',
          '-fno-builtin-acos', '-fno-builtin-atan', '-fno-builtin-atan2',
          '-fno-builtin-sinh', '-fno-builtin-cosh', '-fno-builtin-tanh',
          '-fno-builtin-hypot')

_FALSE = ('0', 'off', 'false', 'no', '')


def enabled_env():
    """The backend the environment asks for, `'numpy'` unless set."""
    return os.environ.get('PYCIRCUIT_HDL_BACKEND', '').strip() or 'numpy'


## ----------------------------------------------------------------------
## Compiler discovery.

_compiler = None


def find_compiler():
    """`(path, version-line)` of the C compiler, or `(None, why)`.

    `$PYCIRCUIT_HDL_CC` wins, then `$CC`, then `cc`/`gcc`/`clang` on
    the PATH.  Memoised per process; the version line is part of every
    cache key, so a compiler upgrade rebuilds.
    """
    global _compiler
    if _compiler is not None:
        return _compiler
    names = []
    for env in ('PYCIRCUIT_HDL_CC', 'CC'):
        v = os.environ.get(env, '').strip()
        if v:
            names.append(v)
    names += ['cc', 'gcc', 'clang']
    for name in names:
        path = shutil.which(name)
        if path is None:
            continue
        try:
            out = subprocess.run([path, '--version'], capture_output=True,
                                 text=True, timeout=60)
        except (OSError, subprocess.SubprocessError) as e:
            _compiler = (None, '%s --version failed: %s' % (path, e))
            return _compiler
        if out.returncode == 0 and out.stdout:
            _compiler = (path, out.stdout.splitlines()[0].strip())
            return _compiler
    _compiler = (None, 'no C compiler on PATH (tried %s)'
                 % ', '.join(names))
    return _compiler


## ----------------------------------------------------------------------
## The .so store.

_tmp_dir = None


def _so_dir():
    """Where built objects live: the compile-cache directory, or a
    per-process temporary directory when the cache is disabled."""
    from pycircuit.circuit import _hdl_cache
    if _hdl_cache.enabled():
        d = _hdl_cache.cache_dir()
        os.makedirs(d, exist_ok=True)
        return d
    global _tmp_dir
    if _tmp_dir is None:
        _tmp_dir = tempfile.mkdtemp(prefix='pycircuit-hdl-c-')
        import atexit
        atexit.register(shutil.rmtree, _tmp_dir, ignore_errors=True)
    return _tmp_dir


def source_key(csrc):
    """SHA-256 of the SOURCE side of the key: format, flags, the C
    kernel prelude and the function itself.  The compiler's identity is
    part of the FILENAME instead (`so_path`), so that a machine with no
    compiler at all can still find and serve an object built from this
    exact source by some compiler (`find_so`)."""
    from pycircuit.circuit import hdl
    h = hashlib.sha256()
    for part in ('format=%d' % SO_FORMAT, ' '.join(CFLAGS),
                 hdl._KERNEL_C, csrc):
        h.update(part.encode('utf-8'))
        h.update(b'\0')
    return h.hexdigest()


def _cc_tag():
    _path, version = find_compiler()
    if version is None:
        return None
    return hashlib.sha256(version.encode('utf-8')).hexdigest()[:12]


def so_path(key):
    """The exact object file for `key` under the CURRENT compiler --
    the build target.  Requires a compiler to name."""
    tag = _cc_tag()
    if tag is None:
        raise CompileError(find_compiler()[1])
    return os.path.join(_so_dir(), 'c-%s-%s.so' % (key, tag))


def find_so(key):
    """An existing object for `key`: the current compiler's if it is
    there, else any compiler's -- same source, same flags, same kernel,
    so the semantics are pinned even when the codegen is another
    compiler's.  None when nothing is stored."""
    d = _so_dir()
    tag = _cc_tag()
    if tag is not None:
        exact = os.path.join(d, 'c-%s-%s.so' % (key, tag))
        if os.path.exists(exact):
            return exact
    try:
        names = sorted(n for n in os.listdir(d)
                       if n.startswith('c-%s-' % key) and
                       n.endswith('.so'))
    except OSError:
        return None
    return os.path.join(d, names[0]) if names else None


class CompileError(Exception):
    pass


@contextlib.contextmanager
def _key_lock(key, blocking=True):
    """Exclusive, across processes, the right to build `key`'s object
    (`fcntl.flock` on `locks/<key>.lock` in the store): yields whether it
    was acquired -- with `blocking=False` it may not be, and the caller
    leaves the build to the holder.  So eight test workers asking for one
    model's kernels compile it once, not eight times (each a `cc` run of
    up to a minute for PSP).  Where `flock` is unavailable it yields True
    and builds proceed unlocked, as before (each still atomic)."""
    try:
        import fcntl
    except ImportError:                            # pragma: no cover
        yield True
        return
    d = os.path.join(_so_dir(), 'locks')
    try:
        os.makedirs(d, exist_ok=True)
        fd = os.open(os.path.join(d, f'c-{key}.lock'),
                     os.O_RDWR | os.O_CREAT, 0o644)
    except OSError:
        yield True
        return
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | (0 if blocking
                                             else fcntl.LOCK_NB))
        except OSError:
            yield False
            return
        try:
            yield True
        finally:
            fcntl.flock(fd, fcntl.LOCK_UN)
    finally:
        os.close(fd)


def _build(csrc, dest):
    """Compile `csrc` (prelude added here) to `dest`, atomically."""
    from pycircuit.circuit import hdl
    cc, _version = find_compiler()
    if cc is None:
        raise CompileError(_version)
    d = os.path.dirname(dest)
    fd, cpath = tempfile.mkstemp(suffix='.c', dir=d)
    sopath = cpath[:-2] + '.so.tmp'
    try:
        with os.fdopen(fd, 'w') as fh:
            fh.write(hdl._KERNEL_C)
            fh.write(csrc)
        proc = subprocess.run([cc] + list(CFLAGS) + ['-o', sopath, cpath,
                                                     '-lm'],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            err = (proc.stderr or proc.stdout or '').strip()
            raise CompileError('%s exited %d: %s'
                               % (cc, proc.returncode,
                                  err.splitlines()[0] if err else '?'))
        os.replace(sopath, dest)
    finally:
        for p in (cpath, sopath):
            try:
                os.unlink(p)
            except OSError:
                pass


def _dlopen(path):
    import cffi
    ffi = cffi.FFI()
    from pycircuit.circuit import hdl
    ffi.cdef('void %s(const double *x, const double *p, double *out);'
             % hdl._C_ENTRY)
    lib = ffi.dlopen(path)
    return ffi, getattr(lib, hdl._C_ENTRY)


## ----------------------------------------------------------------------
## Kernels.

class CKernel(object):
    """One compiled chain function, callable as the element methods
    need it: `kernel(element, x, epar)` -> a fresh float64 array."""

    __slots__ = ('ffi', 'cfn', 'shape', 'nx', 'n_p', 't_index', 'key',
                 'built_s')

    def __init__(self, ffi, cfn, shape, nx, layout, key, built_s):
        self.ffi, self.cfn, self.shape = ffi, cfn, shape
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
        if self.t_index is not None:
            from pycircuit.circuit import hdl
            T = hdl._epar_T(epar)
            if np.ndim(T) != 0:
                return None
            p[self.t_index] = T
        if (type(x) is np.ndarray and x.dtype == np.float64
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
        self.cfn(ffi.cast('double *', xa.ctypes.data), pcast,
                 ffi.cast('double *', out.ctypes.data))
        return out


#: In-process memo: key -> (ffi, cfn), so one function shared between
#: entries (`i_dc` IS `i`) or between classes is loaded once.
_loaded = {}


def kernel_for(fn, nx, rebuild_corrupt=True):
    """The `CKernel` for a chain function carrying `_csrc`.

    Raises `CompileError` when there is no compiler or the build fails;
    the caller turns that into a status, never into a broken class.
    Returns `(kernel, cold)` where `cold` says whether a compile ran.
    """
    csrc = fn._csrc
    key = source_key(csrc)
    cold = False
    import time
    t0 = time.perf_counter()
    if key in _loaded:
        ffi, cfn = _loaded[key]
    else:
        path = find_so(key)
        if path is None:
            ## (under the key's lock, and looked for again inside it: a
            ## process that held it may have built the object meanwhile)
            with _key_lock(key):
                path = find_so(key)
                if path is None:
                    path = so_path(key)     # raises without a compiler
                    _build(csrc, path)
                    cold = True
        try:
            ffi, cfn = _dlopen(path)
        except OSError as e:
            if not rebuild_corrupt:
                raise CompileError('unloadable object: %s' % e)
            ## A truncated or foreign file at our key: rebuild once.
            try:
                os.unlink(path)
            except OSError:
                pass
            path = so_path(key)
            _build(csrc, path)
            cold = True
            try:
                ffi, cfn = _dlopen(path)
            except OSError as e2:
                raise CompileError('rebuilt object still unloadable: %s'
                                   % e2)
        _loaded[key] = (ffi, cfn)
    kern = CKernel(ffi, cfn, tuple(fn._cshape), nx, tuple(fn._clayout),
                   key, time.perf_counter() - t0)
    kern.built_s = kern.built_s if cold else 0.0
    return kern, cold


def _build_missing_parallel(fns):
    """Compile, in parallel, every distinct `_csrc` in `fns` whose
    object is not on disk yet.  Failures fall through to `kernel_for`,
    which reports them; this only warms the store."""
    from pycircuit.circuit import hdl
    cc, _version = find_compiler()
    if cc is None:
        return
    jobs = []
    seen = set()
    with contextlib.ExitStack() as held:
        for fn in fns:
            key = source_key(fn._csrc)
            if key in seen or key in _loaded:
                continue
            seen.add(key)
            if find_so(key) is not None:
                continue
            ## a key another process is building is left to it:
            ## `kernel_for` waits on the lock and then finds its object
            if not held.enter_context(_key_lock(key, blocking=False)):
                continue
            if find_so(key) is not None:
                continue
            dest = so_path(key)
            d = os.path.dirname(dest)
            fd, cpath = tempfile.mkstemp(suffix='.c', dir=d)
            with os.fdopen(fd, 'w') as fh:
                fh.write(hdl._KERNEL_C)
                fh.write(fn._csrc)
            sopath = cpath[:-2] + '.so.tmp'
            proc = subprocess.Popen([cc] + list(CFLAGS) +
                                    ['-o', sopath, cpath, '-lm'],
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL)
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
#: `u`/`dudt`/`uac` (x-free, cheap) stay numpy.
C_FUNCS = ('i', 'G', 'q', 'C', 'i_dc', 'G_dc')


def _note(cls, status, warn=None):
    cls._hdl_backend_status = status
    if warn:
        _sw.warn(f'hdl C backend: {cls.__name__}: {warn}', _sw.CostWarning)


def detach(cls, info):
    """Remove any C bindings so the class runs numpy again."""
    for name in C_FUNCS:
        fn = info['funcs'].get(name)
        if fn is not None and hasattr(fn, '_hdl_c'):
            del fn._hdl_c
    info['_c_bound'] = False


def attach(cls, info):
    """Bind C kernels to `cls`'s chained functions, per the requested
    backend.  Always sets `cls._hdl_backend_status`; never raises for a
    build problem (it warns and falls back to numpy), only for a
    misconfigured backend name."""
    from pycircuit.circuit import hdl
    which = hdl._backend_requested(cls)
    if which != 'c':
        detach(cls, info)
        _note(cls, 'numpy')
        return
    if not info['chained']:
        _note(cls, 'numpy (eager path: the C backend serves let-chain '
                   'models only)')
        return
    funcs = info['funcs']
    todo = {}
    for name in C_FUNCS:
        fn = funcs.get(name)
        if fn is None or id(fn) in todo:
            continue
        if getattr(fn, '_csrc', None) is not None:
            todo[id(fn)] = fn
        elif hasattr(fn, '_creason'):
            _note(cls, 'numpy (C rendering unavailable: %s)' % fn._creason,
                  warn='requested backend "c" but the chain could not be '
                       'printed to C (%s); running numpy' % fn._creason)
            return
    if not todo:
        _note(cls, 'numpy (no C source on this class; recompile with '
                   'hdl.EMIT_C_SOURCE on)',
              warn='requested backend "c" but the compiled functions carry '
                   'no C source; running numpy')
        return
    ## The kernels load through cffi, which numpy does not need: without
    ## it a request for C ran into an ImportError out of class creation
    ## (until 2026-10-02).
    if importlib.util.find_spec('cffi') is None:
        detach(cls, info)
        _note(cls, 'numpy (cffi not installed)',
              warn='requested backend "c" but cffi is not installed; '
                   'running numpy')
        return
    ## No compiler check here: a warm `.so` store serves a machine
    ## with no compiler at all (`kernel_for` only builds on a miss, and
    ## raises `CompileError` naming the missing compiler when it must
    ## build).
    _build_missing_parallel(list(todo.values()))
    cold = False
    try:
        for fn in todo.values():
            ## The x length the kernel reads equals its output row
            ## count: `i`, `G`, `q`, `C` are all n-per-side.  Taken from
            ## the function itself because `cls._hdl_info` is not
            ## assigned yet when the metaclass calls attach.
            kern, was_cold = kernel_for(fn, fn._cshape[0])
            cold = cold or was_cold
            fn._hdl_c = kern
    except CompileError as e:
        detach(cls, info)
        _note(cls, 'numpy (compile failed: %s)' % e,
              warn='requested backend "c" but the build failed (%s); '
                   'running numpy' % e)
        return
    ## (read by `_hdl_cse.take`: a C-bound class's calls never fuse)
    info['_c_bound'] = True
    _note(cls, 'c')
