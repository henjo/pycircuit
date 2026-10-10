# Copyright (c) 2008-2026 Pycircuit Development Team
# See LICENSE for details.
"""Building and loading pycircuit's C objects -- the HDL kernels and the
transient's runtime units alike.

CONTRACT
    `load_kernel(csrc, cdef)` returns `(ffi, cfn, key, cold, seconds)`
    for the C source `csrc`: the object stored under `source_key(csrc)`
    (format, `CFLAGS`, the kernel prelude `hdl._KERNEL_C` and the
    source), built under that key's lock when the store has none, loaded
    through cffi's ABI mode with `cdef` declaring the entry `hdl_fn` (the
    chain kernels' `_fn_cdef()` when None), memoised per process in
    `_loaded`.  `CompileError` when no compiler is found or a build
    fails; a caller turns it into a status.
    `CFLAGS` are strict IEEE doubles: no fast-math, no FMA contraction,
    no `-march`, every transcendental a real libm call -- the bit
    identity with the numpy path rests on them.  With
    `PYCIRCUIT_C_SANITIZE=1` the sanitizer flags are added (their own
    keys) and the import refuses without the ASan runtime preloaded.
    The store is the compile cache's directory, or a per-process temporary
    one with the cache off; an object's file name carries the compiler's
    version tag, and another compiler's object for the same source is
    served where this one is absent (`find_so`).

    Patch these names HERE, not on `_hdl_cbackend`, which re-exports them
    for its callers: `CFLAGS`, `_compiler`, `_loaded` are read from this
    module.

Moved out of `_hdl_cbackend` (the C stepping loop plan's stage R1b,
2026-10-10) byte for byte: every key, and so every stored object, is the
one it was.  History of the pieces: `_hdl_cbackend`'s module docstring.
"""

import contextlib
import hashlib
import os
import re
import shutil
import subprocess
import tempfile

#: The runtime units' C sources, as files (the C loop plan's stage R1c):
#: `csrc(name)` is the text of `csrc/<name>`, read without newline
#: translation, its `#include "x.h"` lines replaced by those headers' text
#: (R1d) -- so what is keyed, compiled and digested is every line compiled.
CSRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc')

#: a line including a header of `csrc/` (a header includes none of them)
_LOCAL_INCLUDE = re.compile(r'^#include "([\w.]+)"\n', re.MULTILINE)


def _read(name):
    with open(os.path.join(CSRC_DIR, name), encoding='utf-8', newline='') as fh:
        return fh.read()


def csrc(name):
    """The text of the C source file `csrc/<name>`, its local headers
    expanded in place."""
    return _LOCAL_INCLUDE.sub(lambda m: _read(m.group(1)), _read(name))


def constants(text):
    """`{name: value}` of the enumerators of a C text (a header's `enum
    {...}` bodies): how Python reads the codes the C says by name.  Each
    enumerator spells its value -- an integer, or earlier names joined by
    `|`."""
    text = re.sub(r'/\*.*?\*/', ' ', text, flags=re.DOTALL)
    out = {}
    for body in re.findall(r'\benum\b[^{;]*\{(.*?)\}', text, flags=re.DOTALL):
        for item in body.split(','):
            if not item.strip():
                continue
            name, eq, expr = (part.strip() for part in item.partition('='))
            if not eq:
                raise ValueError(f'enumerator {name} spells no value')
            value = 0
            for part in expr.split('|'):
                part = part.strip()
                value |= int(part) if part.isdigit() else out[part]
            out[name] = value
    return out


def bails(codes, prefix):
    """`{code: why}` of the `codes` named `prefix<WHY>`: a unit's bails as
    its Python side reports them (`why` lowered)."""
    return {v: k[len(prefix):].lower() for k, v in codes.items() if k.startswith(prefix)}


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

#: THE SANITIZER BUILD (2026-10-04, testing for development, stage 4).
#: With `PYCIRCUIT_C_SANITIZE=1` every C object -- the chain kernels, the
#: limiters, the pass driver, the walk, the evaluate core: all build here --
#: is compiled with AddressSanitizer and UndefinedBehaviorSanitizer, so an
#: out-of-bounds access or undefined behaviour stops the process with the
#: file and line instead of corrupting a neighbour silently.  No float
#: checks: the kernels make inf and NaN on purpose, where Python raises.
#: The flags are part of `source_key`, so sanitized objects sit in the same
#: store under their own keys and a normal process never loads one.  The
#: interpreter is not instrumented: the ASan runtime must be PRELOADED
#: (`scripts/sanitize_suite.sh` sets everything).  Without it every object
#: would fail to load, `load_kernel` would rebuild and fail again, every
#: class would fall back to numpy and every test would pass having tested
#: nothing -- so the import refuses instead.
SANITIZE = os.environ.get('PYCIRCUIT_C_SANITIZE', '').strip().lower() not in _FALSE
SANITIZE_FLAGS = ('-g', '-fno-omit-frame-pointer', '-fsanitize=address,undefined',
                  '-fno-sanitize-recover=all')


def _asan_preloaded():
    import ctypes
    try:
        return hasattr(ctypes.CDLL(None), '__asan_init')
    except OSError:
        return False


if SANITIZE:
    if not _asan_preloaded():
        raise RuntimeError(
            'PYCIRCUIT_C_SANITIZE=1 needs the AddressSanitizer runtime preloaded '
            '(LD_PRELOAD="$(gcc -print-file-name=libasan.so) $(gcc '
            '-print-file-name=libstdc++.so)"): without it no '
            'sanitized object loads and every class would run numpy. '
            'Use scripts/sanitize_suite.sh.')
    CFLAGS = CFLAGS + SANITIZE_FLAGS


def _compiler_env():
    """The environment a compiler runs in: this process's, without the
    sanitizer runtime preload (gcc is not built for it)."""
    if not SANITIZE:
        return None
    env = dict(os.environ)
    env.pop('LD_PRELOAD', None)
    return env


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
            out = subprocess.run([path, '--version'], capture_output=True, env=_compiler_env(),
                                 text=True, timeout=60, check=False)
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
                              capture_output=True, text=True, env=_compiler_env())
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


def _fn_cdef():
    """The chain functions' entry: `void hdl_fn(x, p, out)`."""
    from pycircuit.circuit import hdl
    return f'void {hdl._C_ENTRY}(const double *x, const double *p, double *out);'


def _dlopen(path, cdef=None):
    """`(ffi, cfn)` for the object at `path`: the chain functions' entry
    by default, or the entry `cdef` declares (the limiter's `int hdl_fn`,
    `_hdl_climit`) -- the one name, another signature; the object is
    keyed by its source, so the two never meet in one file."""
    import cffi
    ffi = cffi.FFI()
    from pycircuit.circuit import hdl
    ffi.cdef(cdef or _fn_cdef())
    lib = ffi.dlopen(path)
    return ffi, getattr(lib, hdl._C_ENTRY)


#: In-process memo: key -> (ffi, cfn), so one function shared between
#: entries (`i_dc` IS `i`) or between classes is loaded once.
_loaded = {}


def load_kernel(csrc, cdef=None, rebuild_corrupt=True):
    """The loaded object for the C source `csrc` (the prelude is added):
    `(ffi, cfn, key, cold, seconds)` -- built under the key's lock when
    the store has no object for it, loaded through `cdef`'s entry (the
    chain functions' by default; `_hdl_climit`'s otherwise), memoised per
    process by the source's key.

    Raises `CompileError` when there is no compiler or the build fails;
    the caller turns that into a status, never into a broken class.
    `cold` says whether a compile ran; `seconds` is its time, 0 otherwise.
    """
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
            ffi, cfn = _dlopen(path, cdef)
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
                ffi, cfn = _dlopen(path, cdef)
            except OSError as e2:
                raise CompileError('rebuilt object still unloadable: %s'
                                   % e2)
        _loaded[key] = (ffi, cfn)
    return ffi, cfn, key, cold, (time.perf_counter() - t0 if cold else 0.0)
