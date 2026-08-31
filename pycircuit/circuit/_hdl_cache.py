# -*- coding: utf-8 -*-
# Copyright (c) 2008-2026 Pycircuit Development Team
# See LICENSE for details.
"""On-disk cache of compiled `Behavioural` elements.

Compiling a `Behavioural` class -- running `analog()`, assembling the
residual vectors symbolically, differentiating them and printing the
result to Python -- costs from milliseconds for a resistor to about a
second for a Gummel-Poon BJT and a minute for a surface-potential
MOSFET, and every interpreter that imports the library pays it again.
Everything that work produces is text and metadata: generated Python
source for each compiled function, and a dict of plain data and sympy
expressions.  This module stores that output, keyed on everything it
depends on, and hands it back on the next start so that no symbolic
work runs at all.

What is cached
    The whole `_hdl_info` dict `generate_code` returns, with every
    compiled function replaced by a record of its source text and the
    names it needs from its namespace.  On a hit the source is executed
    again in a namespace rebuilt the same way, so the functions are the
    same code bound to the same objects and their results are
    bit-identical to a cold compile's.  The symbolic vectors (`sym_spec`
    and `pure_spec`, the eager path's exact forms) are pickled as sympy
    expressions and come back equal.  The per-mask classes a collapsing
    model builds through `_collapse_variant` are separate classes with
    their own key, so they are cached too.

What is not cached
    The JAX twins (`eval_i_pure`, the PCNR `_jax` entries) are compiled
    lazily at first traced use, never at import, and are not stored.  A
    class whose `analog()` source `inspect` cannot recover (defined in a
    REPL or an `exec` string) compiles cold every time and says so in
    `cls._hdl_cache_status`.  A compile that raises is not cached, so
    the error is raised again on the next start.  Warnings the compile
    emits are not replayed on a hit; `hdl.py` itself emits none.

Key
    SHA-256 over: the cache format number, the source of `hdl.py` and of
    this module, sympy/numpy/Python versions, the class's module and
    qualified name, `inspect.getsource(cls.analog)`, every `instparams`
    declaration (name, default, range, unit, description), the collapse
    mask, a class-body `terminals` declaration, and the source of every
    module the analog body reaches by name (the class's own module and
    any pycircuit module a global it references was defined in).  Any
    of these changing misses.  What the key does NOT see: a helper two
    hops away in a module the analog body does not name, and code
    outside the `pycircuit` package.

Location
    `$PYCIRCUIT_HDL_CACHE_DIR`, else `$XDG_CACHE_HOME/pycircuit/hdl`,
    else `~/.cache/pycircuit/hdl`.  Never the source tree.  One pickle
    per key, written to a temporary file in the same directory and
    moved into place with `os.replace`, so two processes compiling the
    same class cannot leave a torn file and a reader never sees a
    partial one.  An unreadable, truncated or foreign file is a miss.

Switches
    `PYCIRCUIT_HDL_CACHE=0` (or `off`, `false`, `no`) disables the
    cache for the process; so does setting `ENABLED = False` here.
    `PYCIRCUIT_HDL_CACHE_DEBUG=1` warns whenever a class could not be
    cached, with the reason.

The pickle is loaded only from the user's own cache directory.  It is
not a format for exchanging compiled models between machines or users.
"""

import builtins
import dis
import hashlib
import inspect
import linecache
import opcode
import os
import pickle
import sys
import tempfile
import types
import warnings
import weakref

from pycircuit.utilities import param

#: Bump to invalidate every entry written by an older layout of the
#: records below.  (The hdl.py source hash already invalidates on any
#: change to the compiler; this is for changes to THIS module's format.)
#: 2: chain records carry the C rendering (`_csrc`/`_cshape`/`_clayout`
#: or `_creason`) beside the numpy source, so a cache hit can feed the
#: C backend without re-running sympy.
## 3 (2026-08-31): the generated `pcnr_limit_branchless` and the
## `_hdl_limit_par` ingredients that let a limiter parameter be recompiled
## for a traced toolkit (roadmap sec. 49.3).  Both are EMISSION changes, and
## a format bump is how emission changes invalidate: the `_dependency_hashes`
## do cover `hdl.py`'s source, but an entry restored from cache does not run
## `_limit_par_fn` again, so the attribute it now attaches was simply absent
## on every cache hit -- the generated device kept working and silently could
## not be traced.  Same shape as the PCNR_LIFT_AFFINE trap recorded below,
## reached from the other direction.
## 4 (2026-08-31): the frozen function record carries `_hdl_limit_par` as
## well as `_wants_x`, so a limiter parameter restored from cache can still
## be recompiled for a traced toolkit.  Bumping 2 -> 3 invalidated the stale
## entries but the newly written ones dropped the ingredients just the same,
## which is the difference between invalidating a cache and fixing what it
## stores.
CACHE_FORMAT = 4

#: Module-level switch; the environment variable is consulted too.
ENABLED = True

_FALSE = ('0', 'off', 'false', 'no', '')


def enabled():
    if not ENABLED:
        return False
    return os.environ.get('PYCIRCUIT_HDL_CACHE', '1').strip().lower() \
        not in _FALSE


def _debug():
    return os.environ.get('PYCIRCUIT_HDL_CACHE_DEBUG', '').strip().lower() \
        not in _FALSE


def cache_dir():
    d = os.environ.get('PYCIRCUIT_HDL_CACHE_DIR')
    if d:
        return os.path.expanduser(d)
    xdg = os.environ.get('XDG_CACHE_HOME')
    base = xdg if xdg else os.path.join(os.path.expanduser('~'), '.cache')
    return os.path.join(base, 'pycircuit', 'hdl')


class Uncacheable(Exception):
    """This compile output cannot be represented in the cache."""


## ----------------------------------------------------------------------
## The key.

_file_hash_cache = {}


def _file_hash(path):
    """SHA-256 of a source file, computed once per process per path."""
    try:
        st = os.stat(path)
        tag = (path, st.st_mtime_ns, st.st_size)
    except OSError:
        return None
    if tag in _file_hash_cache:
        return _file_hash_cache[tag]
    try:
        with open(path, 'rb') as fh:
            h = hashlib.sha256(fh.read()).hexdigest()
    except OSError:
        h = None
    _file_hash_cache[tag] = h
    return h


def _compiler_hash():
    """Source of hdl.py and of this module."""
    from pycircuit.circuit import hdl
    parts = []
    for mod in (hdl, sys.modules[__name__]):
        h = _file_hash(mod.__file__)
        if h is None:
            return None
        parts.append(h)
    return '/'.join(parts)


def _dependency_hashes(func):
    """Hashes of the pycircuit source files the analog body reaches.

    The body's own module always; and the module of every pycircuit
    function or class it refers to by a global name.  A helper the body
    calls lives in one of these, so editing it misses.
    """
    files = set()
    try:
        files.add(inspect.getsourcefile(func))
    except TypeError:
        return None
    g = func.__globals__
    objs = [g.get(name) for name in func.__code__.co_names]
    for cell in func.__closure__ or ():
        try:
            objs.append(cell.cell_contents)
        except ValueError:
            pass
    for obj in objs:
        if obj is None:
            continue
        modname = getattr(obj, '__module__', None)
        if not (isinstance(modname, str) and
                modname.startswith('pycircuit')):
            continue
        mod = sys.modules.get(modname)
        path = getattr(mod, '__file__', None)
        if path:
            files.add(path)
    out = []
    for path in sorted(p for p in files if p):
        h = _file_hash(path)
        if h is None:
            return None
        out.append('%s=%s' % (os.path.basename(path), h))
    return out


def _param_decl(p):
    return repr((p.name, p.default, p.minval, p.maxval, p.unit, p.desc))


_ATOMS = (int, float, complex, str, bytes, bool, type(None))
_CONTAINERS = (dict, list, set)


def _value_sig(name, val, depth, seen, in_closure):
    """A stable text for one value the analog body can see.

    The source of `analog()` is not its whole identity: a body that
    reads `which` from an enclosing function, or `SCALE` from its
    module, compiles differently for each value with the same text --
    and the test suite builds classes exactly that way.  Anything
    whose value has no stable text is refused, naming the variable.
    """
    if isinstance(val, types.ModuleType):
        return 'module:%s' % val.__name__
    if isinstance(val, type):
        return 'class:%s.%s' % (val.__module__, val.__qualname__)
    if isinstance(val, types.FunctionType):
        if id(val) in seen or depth > 3 or _is_library(val):
            ## Library code is keyed by its module's source (the
            ## dependency hashes); walking its globals would drag in
            ## counters and registries that are not part of any model.
            return 'func:%s.%s' % (val.__module__, val.__qualname__)
        seen.add(id(val))
        try:
            src = inspect.getsource(val)
        except (OSError, TypeError):
            raise Uncacheable('analog() calls %s(), whose source is not '
                              'recoverable' % name)
        return 'func:%s.%s:%s:{%s}' % (
            val.__module__, val.__qualname__,
            hashlib.sha256(src.encode('utf-8')).hexdigest(),
            ';'.join(_bindings(val, depth + 1, seen)))
    if isinstance(val, (types.BuiltinFunctionType, types.MethodType)):
        return 'callable:%s' % getattr(val, '__qualname__', repr(val))
    if isinstance(val, _ATOMS):
        return repr(val)
    if isinstance(val, (tuple, frozenset)):
        return '(%s)' % ','.join(_value_sig(name, v, depth, seen, in_closure)
                                 for v in val)
    if isinstance(val, param.Parameter):
        return _param_decl(val)
    try:
        import sympy
        if isinstance(val, sympy.Basic):
            return 'sympy:' + sympy.srepr(val)
    except ImportError:
        pass
    try:
        import numpy
        if isinstance(val, (numpy.ndarray, numpy.generic)):
            arr = numpy.asarray(val)
            return 'ndarray:%s%r:%s' % (arr.dtype, arr.shape,
                                        arr.tobytes().hex())
    except ImportError:
        pass
    if isinstance(val, _CONTAINERS):
        ## Keyed by content.  A body that WRITES into one is refused
        ## separately (`_writes_outside`), because no key can replay a
        ## side effect.
        text = repr(val)
        if ' at 0x' in text:
            raise Uncacheable('analog() reads %s, whose repr has no stable '
                              'text' % name)
        return '%s:%s' % (type(val).__name__, text)
    raise Uncacheable('analog() reads %s, a %s, whose value has no stable '
                      'text for the cache key' % (name, type(val).__name__))


def _is_library(func):
    mod = getattr(func, '__module__', None) or ''
    return mod.startswith('pycircuit') and '.tests' not in mod


_STORES = {'STORE_SUBSCR': 3, 'STORE_ATTR': 2, 'DELETE_SUBSCR': 3,
           'DELETE_ATTR': 2}
_OUTER_LOADS = ('LOAD_DEREF', 'LOAD_GLOBAL', 'LOAD_NAME', 'LOAD_CLASSDEREF',
                'LOAD_FROM_DICT_OR_DEREF')


def _writes_outside(code):
    """The name of something outside the body's own locals that the
    body stores into, or None.

    `seen['g'] = ...` compiles to a load of `seen` two instructions
    before the `STORE_SUBSCR`; an attribute store, one before.  A body
    that does this is instrumenting itself -- a diagnostics test does,
    to see the parameter injection -- and a cache hit, which never runs
    the body, cannot replay it.  Stores into the body's own locals are
    its business.
    """
    ins = list(dis.get_instructions(code))
    for k, i in enumerate(ins):
        back = _STORES.get(i.opname)
        if back is None:
            continue
        for j in range(max(0, k - back), k):
            if ins[j].opname in _OUTER_LOADS:
                return ins[j].argval
    for c in code.co_consts:
        if isinstance(c, types.CodeType):
            found = _writes_outside(c)
            if found is not None:
                return found
    return None


def _bindings(func, depth=0, seen=None):
    """Every name `func` resolves outside its own locals, with its
    value's signature: closure cells, defaults, and the globals its code
    names."""
    if seen is None:
        seen = set()
    out = []
    if func.__closure__:
        for name, cell in zip(func.__code__.co_freevars, func.__closure__):
            try:
                val = cell.cell_contents
            except ValueError:
                out.append('cell:%s=<empty>' % name)
                continue
            out.append('cell:%s=%s' % (name, _value_sig(name, val, depth,
                                                        seen, True)))
    for name, val in zip(reversed(func.__code__.co_varnames),
                         reversed(func.__defaults__ or ())):
        out.append('default:%s=%s' % (name, _value_sig(name, val, depth,
                                                       seen, True)))
    for name, val in sorted((func.__kwdefaults__ or {}).items()):
        out.append('kwdefault:%s=%s' % (name, _value_sig(name, val, depth,
                                                         seen, True)))
    g = func.__globals__
    named = set()

    def globals_of(code):
        for name in code.co_names:
            if name in g and name not in named:
                named.add(name)
                out.append('global:%s=%s' % (
                    name, _value_sig(name, g[name], depth, seen, False)))
        for c in code.co_consts:
            if isinstance(c, types.CodeType):
                ## A nested function or comprehension inside the body
                ## reads the same globals.
                globals_of(c)
    globals_of(func.__code__)
    return out


def key_for(cls):
    """The cache key of `cls`.  Raises `Uncacheable` when the class
    cannot be keyed, saying why."""
    func = cls.__dict__.get('analog', getattr(cls, 'analog', None))
    if isinstance(func, (staticmethod, classmethod)):
        func = func.__func__
    if not isinstance(func, types.FunctionType):
        raise Uncacheable('analog() is not a plain function')
    try:
        src = inspect.getsource(func)
    except (OSError, TypeError):
        raise Uncacheable('analog() source not recoverable')
    written = _writes_outside(func.__code__)
    if written is not None:
        raise Uncacheable('analog() writes into %s, which is not its own '
                          'local -- a cache hit would never run the body, '
                          'so that write would not happen' % written)
    comp = _compiler_hash()
    deps = _dependency_hashes(func)
    if comp is None or deps is None:
        raise Uncacheable('a source file could not be read')
    bindings = _bindings(func)
    import numpy
    import sympy
    from pycircuit.circuit import hdl as _hdl_mod
    parts = [
        'format=%d' % CACHE_FORMAT,
        'compiler=%s' % comp,
        'sympy=%s' % sympy.__version__,
        'numpy=%s' % numpy.__version__,
        'python=%s' % sys.version,
        'class=%s.%s' % (cls.__module__, cls.__qualname__),
        'terminals=%r' % (cls.__dict__.get('terminals'),),
        'params_as=%r' % (getattr(cls, 'params_as', None),),
        'mask=%r' % (getattr(cls, '_hdl_collapse_mask', None),),
        ## A RUNTIME flag that changes what is compiled must be in the
        ## key.  The dependency hashes cover `hdl.py`'s SOURCE, so an
        ## edit invalidates -- but flipping `PCNR_LIFT_AFFINE` in a
        ## running process changes the emitted model while leaving every
        ## file byte-identical.  Without this line the two builds share
        ## a key, and the second silently gets the first one's code.
        ## Found the hard way: a probe run with the lift on wrote 924
        ## entries that then served the lift-off suite, and the symptom
        ## was a model that stopped being refused.
        'pcnr_lift=%r' % (getattr(_hdl_mod, 'PCNR_LIFT_AFFINE', False),),
        ## Card-constant folding (roadmap sec. 30): two cards fold to two
        ## DIFFERENT compiled models from identical `analog` source and
        ## identical `instparams`, so without this they would share a key
        ## and the second card would silently run the first card's code.
        ## `repr` of a float is exact and round-trips, so the key is
        ## stable across runs.
        'fold=%r' % (tuple(sorted(
            (getattr(cls, '_hdl_fold_values', None) or {}).items())),),
        'params=' + '|'.join(_param_decl(p)
                             for p in getattr(cls, 'instparams', ())),
        'deps=' + '|'.join(deps),
        'bindings=' + '|'.join(bindings),
        'analog=' + src,
    ]
    return hashlib.sha256('\n'.join(parts).encode('utf-8')).hexdigest()


## ----------------------------------------------------------------------
## Freezing compiled functions to records, and thawing them back.

_LOAD_GLOBAL = opcode.opmap['LOAD_GLOBAL']
_LOAD_NAME = opcode.opmap.get('LOAD_NAME')
_EXTENDED_ARG = opcode.opmap['EXTENDED_ARG']


def _global_loads(code):
    """Every name a code object (and its nested ones) loads as a global.

    A raw scan of the bytecode rather than `dis.get_instructions`,
    which builds an Instruction per opcode and cost 1.4 s over the
    library's 500 generated functions -- a quarter of the compile it
    was meant to save.  `co_names` alone will not do: it also lists
    attribute names, and `numpy.where` must not be mistaken for a
    global `where`.  `test_hdl_cache` checks this scan against `dis`.
    """
    out = set()
    names = code.co_names
    raw = code.co_code
    ext = 0
    for k in range(0, len(raw), 2):
        op = raw[k]
        arg = raw[k + 1] | ext
        if op == _EXTENDED_ARG:
            ext = arg << 8
            continue
        ext = 0
        if op == _LOAD_GLOBAL:
            ## Since 3.11 the low bit asks for a NULL push; the name
            ## index is the rest.
            out.add(names[arg >> 1])
        elif op == _LOAD_NAME:
            out.add(names[arg])
    for c in code.co_consts:
        if isinstance(c, types.CodeType):
            out |= _global_loads(c)
    return out


def _global_loads_dis(code):
    """The same set via `dis`; the reference `_global_loads` is tested
    against."""
    out = set()
    for ins in dis.get_instructions(code):
        if ins.opname in ('LOAD_GLOBAL', 'LOAD_NAME'):
            out.add(ins.argval)
    for c in code.co_consts:
        if isinstance(c, types.CodeType):
            out |= _global_loads_dis(c)
    return out


def _check_namespace(fn, ns, what):
    """Refuse a record whose rebuilt namespace would bind any global the
    function loads to a different object than the compile did."""
    g = fn.__globals__
    b = vars(builtins)
    for name in _global_loads(fn.__code__):
        if name in ns:
            if name in g and g[name] is not ns[name]:
                raise Uncacheable('%s binds %r to a different object than '
                                  'the rebuilt namespace would' % (what, name))
            if name not in g and name not in b:
                raise Uncacheable('%s loads %r, unbound at compile'
                                  % (what, name))
        elif name in g:
            raise Uncacheable('%s needs %r, which the rebuilt namespace '
                              'does not provide' % (what, name))
        elif name not in b:
            raise Uncacheable('%s loads unbound name %r' % (what, name))


def _imp_registry():
    """Name -> `_imp_` of every DSL function class that has one."""
    from pycircuit.circuit import hdl
    out = {}
    for name, obj in vars(hdl).items():
        imp = getattr(obj, '_imp_', None)
        if isinstance(obj, type) and imp is not None:
            out[name] = imp
    return out


def _lambdify_namespace(imports, extras):
    """The namespace `sympy.lambdify(modules=NUMPY_MODULES)` builds.

    Same order as lambdify: numpy's namespace, then the DSL's kernel map
    on top, then the implemented functions the expression carried
    (`_imp_`), then the printer's module imports, then `builtins` and
    `range`.
    """
    from sympy.utilities.lambdify import _import, MODULES
    from pycircuit.circuit import hdl
    import numpy as np
    _import('numpy')
    ns = {}
    ns.update(MODULES['numpy'][0])
    ns.update(dict(hdl._KERNEL_NUMPY, _wrapfloor=np.floor))
    if extras:
        reg = _imp_registry()
        for name in extras:
            if name not in reg:
                raise Uncacheable('no implemented function named %r' % name)
            ns[name] = reg[name]
    for ln in imports:
        exec(ln, {}, ns)
    ns.update({'builtins': builtins, 'range': range})
    return ns


def _lambdify_imports(fn):
    doc = fn.__doc__ or ''
    marker = 'Imported modules:'
    k = doc.find(marker)
    if k < 0:
        return []
    return [ln.strip() for ln in doc[k + len(marker):].splitlines()
            if ln.strip()]


_counter = [0]


def _freeze_function(fn):
    inner = getattr(fn, '_hdl_inner', None)
    if inner is not None:
        return ('first', _freeze_function(inner))
    src = getattr(fn, '_src', None)
    if src is not None:
        from pycircuit.circuit import hdl
        _check_namespace(fn, hdl._chain_namespace(None), 'chain function')
        ## The C rendering rides along as plain data: text and two small
        ## tuples.  Without it a hit could never serve the C backend --
        ## the sympy lists it is printed from no longer exist.
        cmeta = None
        if getattr(fn, '_csrc', None) is not None:
            cmeta = (fn._csrc, tuple(fn._cshape), tuple(fn._clayout))
        elif getattr(fn, '_creason', None) is not None:
            cmeta = ('', getattr(fn, '_creason'), None)
        return ('chain', src, fn.__name__, cmeta)
    rec = getattr(fn, '_hdl_lambdify', None)
    if rec is not None:             # rebuilt from a record: record again
        return rec
    if (fn.__doc__ or '').startswith('Created with lambdify'):
        try:
            src = inspect.getsource(fn)
        except (OSError, TypeError):
            raise Uncacheable('lambdify source not recoverable')
        imports = _lambdify_imports(fn)
        base = _lambdify_namespace(imports, ())
        g = fn.__globals__
        ## The implemented functions (`_imp_`) lambdify bound from the
        ## expression: every global the function loads that the base
        ## namespace lacks or binds to something else.  They override
        ## the base, as they did in lambdify's own namespace order.
        extras = sorted(n for n in _global_loads(fn.__code__)
                        if n in g and (n not in base or base[n] is not g[n]))
        ns = _lambdify_namespace(imports, extras)
        _check_namespace(fn, ns, 'lambdify function')
        return ('lambdify', src, imports, extras)
    raise Uncacheable('function %r is neither chain-compiled nor '
                      'lambdified' % (fn.__name__,))


def _register_source(src, filename):
    linecache.cache[filename] = (len(src), None, src.splitlines(True),
                                 filename)


def _thaw_function(rec):
    from pycircuit.circuit import hdl
    kind = rec[0]
    if kind == 'first':
        return hdl._first_of(_thaw_function(rec[1]))
    if kind == 'chain':
        _, src, name, cmeta = rec
        ns = hdl._chain_namespace(None)
        exec(compile(src, '<hdl-chain>', 'exec'), ns)
        fn = ns[name]
        fn._src = src
        if cmeta is not None:
            if cmeta[0]:
                fn._csrc, fn._cshape, fn._clayout = cmeta
            else:
                fn._creason = cmeta[1]
        return fn
    if kind == 'lambdify':
        _, src, imports, extras = rec
        ns = _lambdify_namespace(imports, extras)
        _counter[0] += 1
        filename = '<hdl-cache-lambdify-%d>' % _counter[0]
        _register_source(src, filename)
        loc = {}
        exec(compile(src, filename, 'exec'), ns, loc)
        fn = loc['_lambdifygenerated']
        fn._hdl_lambdify = rec
        weakref.finalize(fn, linecache.cache.pop, filename, None)
        return fn
    raise Uncacheable('unknown record kind %r' % (kind,))


_FN = '__hdl_fn__'


def freeze(obj, memo=None):
    """`info` with every function replaced by a picklable record.

    A function that appears more than once (`i_dc` IS `i` when nothing
    is pinned) is recorded once and referred to by number after that,
    so the thawed dict shares it the same way.
    """
    if memo is None:
        memo = {}
    if isinstance(obj, types.FunctionType):
        token = memo.get(id(obj))
        if token is not None:
            return (_FN, token, None, None, None)
        token = len(memo)
        memo[id(obj)] = token
        ## TWO attributes the compiler hangs on a function that its callers
        ## read.  `_wants_x` is read by `limit()`.  `_hdl_limit_par` is the
        ## sympy ingredients `_limit_par_for` recompiles a limiter parameter
        ## from for a traced toolkit -- and it MUST travel, because the thaw
        ## path below rebuilds the function without running `_limit_par_fn`,
        ## so a warm process would otherwise hold working code that cannot be
        ## traced.  Measured before this was carried: the same class reported
        ## the ingredients present on a first run into an empty cache dir and
        ## absent on the second (roadmap sec. 49.3).
        ##
        ## The COMPILED jax twin is deliberately not carried -- it is derived,
        ## and `_jax` is skipped below for the same reason.
        return (_FN, token, _freeze_function(obj),
                getattr(obj, '_wants_x', None),
                getattr(obj, '_hdl_limit_par', None))
    if isinstance(obj, dict):
        return {k: freeze(v, memo) for k, v in obj.items() if k != '_jax'}
    if isinstance(obj, list):
        return [freeze(v, memo) for v in obj]
    if isinstance(obj, tuple):
        return tuple(freeze(v, memo) for v in obj)
    if isinstance(obj, (types.BuiltinFunctionType, types.MethodType)):
        raise Uncacheable('callable %r cannot be recorded' % (obj,))
    return obj


def thaw(obj, memo=None):
    if memo is None:
        memo = {}
    if isinstance(obj, tuple) and len(obj) == 5 and obj[0] == _FN:
        _, token, rec, wants_x, limit_par = obj
        if rec is None:
            return memo[token]
        fn = _thaw_function(rec)
        if wants_x is not None:
            fn._wants_x = wants_x
        if limit_par is not None:
            try:
                fn._hdl_limit_par = limit_par
            except AttributeError:              # pragma: no cover
                pass
        memo[token] = fn
        return fn
    if isinstance(obj, dict):
        return {k: thaw(v, memo) for k, v in obj.items()}
    if isinstance(obj, list):
        return [thaw(v, memo) for v in obj]
    if isinstance(obj, tuple):
        return tuple(thaw(v, memo) for v in obj)
    return obj


## ----------------------------------------------------------------------
## Files.

def _path(key):
    return os.path.join(cache_dir(), key + '.pkl')


def _read(key):
    """`(payload, reason)`: the frozen payload for `key`, or None and
    why not.  Every failure is a miss; the reason is for the status."""
    path = _path(key)
    if not os.path.exists(path):
        return None, 'no entry'
    try:
        with open(path, 'rb') as fh:
            entry = pickle.load(fh)
    except Exception as e:
        return None, 'unreadable entry: %s' % _brief(e)
    if not isinstance(entry, dict) or entry.get('key') != key or \
            entry.get('format') != CACHE_FORMAT:
        return None, 'foreign entry'
    return entry['payload'], 'hit'


def _brief(e):
    return '%s: %s' % (type(e).__name__, e)


def _write(key, payload):
    d = cache_dir()
    os.makedirs(d, exist_ok=True)
    data = pickle.dumps(dict(format=CACHE_FORMAT, key=key, payload=payload),
                        protocol=pickle.HIGHEST_PROTOCOL)
    fd, tmp = tempfile.mkstemp(prefix='.' + key[:16] + '-', suffix='.tmp',
                               dir=d)
    try:
        with os.fdopen(fd, 'wb') as fh:
            fh.write(data)
        os.replace(tmp, _path(key))
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def clear():
    """Delete every entry in the cache directory."""
    d = cache_dir()
    if not os.path.isdir(d):
        return 0
    n = 0
    for name in os.listdir(d):
        if name.endswith(('.pkl', '.so')):
            try:
                os.unlink(os.path.join(d, name))
                n += 1
            except OSError:
                pass
    return n


## ----------------------------------------------------------------------
## The entry point the metaclass calls.

def _note(cls, status):
    cls._hdl_cache_status = status
    if _debug() and status not in ('hit', 'miss'):
        warnings.warn('hdl cache: %s: %s' % (cls.__name__, status))


def compiled_info(cls, compile_fn):
    """`compile_fn(cls)`'s result, from the cache when it is there.

    Sets `cls._hdl_cache_status` to one of ``'hit'``, ``'miss'``,
    ``'miss (<reason>)'`` (an entry was there and could not be used --
    unreadable, foreign, or one that would not rebuild -- and has been
    replaced), ``'disabled'``, ``'uncacheable: <reason>'`` (the class
    was compiled but its output could not be recorded) or
    ``'unsaved: <reason>'`` (recorded, but the file could not be
    written).
    """
    if not enabled():
        info = compile_fn(cls)
        _note(cls, 'disabled')
        return info
    try:
        key = key_for(cls)
    except Uncacheable as e:
        info = compile_fn(cls)
        _note(cls, 'uncacheable: %s' % e)
        return info
    except Exception as e:      # a key that cannot be computed at all
        info = compile_fn(cls)
        _note(cls, 'uncacheable: %s' % _brief(e))
        return info
    payload, why = _read(key)
    if payload is not None:
        try:
            info = thaw(payload)
        except Exception as e:
            ## An entry that unpickled but will not rebuild -- the
            ## reason goes in the status, because a silent miss here
            ## looks exactly like a key that never matches.
            why = 'entry would not rebuild: %s' % _brief(e)
        else:
            _note(cls, 'hit')
            return info
    info = compile_fn(cls)
    try:
        payload = freeze(info)
    except Uncacheable as e:
        _note(cls, 'uncacheable: %s' % e)
        return info
    try:
        _write(key, payload)
    except Exception as e:
        _note(cls, 'unsaved: %s' % _brief(e))
        return info
    _note(cls, 'miss' if why == 'no entry' else 'miss (%s)' % why)
    return info
