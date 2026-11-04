"""Names that become sympy objects on their first symbolic use.

A numeric run of a cached model builds no sympy expression, and importing
sympy is ~0.15 s of a short deck's whole process.  So the modules a numeric
run imports bind sympy -- the package itself (`module`), names taken from it,
and the HDL language's own sympy classes and symbols (`_hdl_sym`) -- as
STAND-INS (`standin`): objects that import their target on the first use
that needs it (a call, an attribute, arithmetic, a comparison, a hash,
`_sympy_`) and forward to it.

The moment sympy is imported -- by anyone: a stand-in's use, a cold compile,
a test, a user's own ``import sympy`` -- every stand-in is resolved and
replaced by its target in every module's namespace (`settle`).  Code that
runs on sympy objects therefore meets the real names: no sympy object can
exist while a module-level stand-in is left in place.  A stand-in captured
before that -- in a closure, a default, a container -- stays, and forwards;
`resolve_function` replaces one in a function's closure and defaults.

``isinstance(x, standin)`` is False while sympy is not imported (no object is
an instance of a sympy class then) and the real test after.

The compile cache keys a stand-in as its target (`describe`): a class keyed
by a process that has imported sympy and by one that has not is one key.
"""

import importlib
import importlib.abc
import operator
import sys
import threading
import types

_LOCK = threading.RLock()
#: The stand-ins not yet resolved.
_PENDING = []
_UNSET = object()
#: Set once sympy's own module has run (`_Then`): it is still marked as
#: importing until its import returns, and the stand-ins settle before that.
_EXECUTED = [False]


def _import(name):
    return importlib.import_module(name)


def _ready(after):
    """Whether the module `after` names (or nothing) has finished importing."""
    if after is None:
        return True
    mod = sys.modules.get(after)
    spec = getattr(mod, '__spec__', None)
    return mod is not None and not getattr(spec, '_initializing', False)


def _value_of(target):
    """The object a stand-in's `target` names: ``(module, attr, call)``."""
    mod, attr, call = target
    obj = _import(mod)
    if attr is not None:
        obj = getattr(obj, attr)
    if call:
        obj = obj()
    return obj


class Standin:
    """The object ``getattr(import_module(module), attr)`` (called, with
    `call`) stands for, imported on first use.  `sig` and `dep` are what the
    compile cache keys it as: `_hdl_cache._value_sig`'s text for the target
    and the module its dependency hashes name (or None)."""
    __slots__ = ('__weakref__', '_after', '_dep', '_obj', '_sig', '_target')

    def __init__(self, target, sig, dep, after):
        self._target = target
        self._obj = _UNSET
        self._sig = sig
        self._dep = dep
        self._after = after

    def _get(self):
        obj = self._obj
        if obj is _UNSET:
            import sympy  # noqa: F401 -- its import settles every stand-in
            settle()
            obj = self._obj
            if obj is _UNSET:       # its module still importing
                obj = _value_of(self._target)
        return obj

    def __getattr__(self, name):
        if name.startswith('__') and name.endswith('__'):
            raise AttributeError(name)
        return getattr(self._get(), name)

    def __reduce__(self):
        return (_value_of, (self._target,))

    def __repr__(self):
        return repr(self._get())

    def __str__(self):
        return str(self._get())

    def __call__(self, *args, **kwargs):
        return self._get()(*args, **kwargs)

    def _sympy_(self):
        return self._get()

    def __instancecheck__(self, obj):
        if self._obj is _UNSET and 'sympy' not in sys.modules:
            return False
        return isinstance(obj, self._get())

    def __subclasscheck__(self, cls):
        if self._obj is _UNSET and 'sympy' not in sys.modules:
            return False
        return issubclass(cls, self._get())

    def __hash__(self):
        return hash(self._get())

    def __eq__(self, other):
        return self._get() == other

    def __ne__(self, other):
        return self._get() != other

    def __lt__(self, other):
        return self._get() < other

    def __le__(self, other):
        return self._get() <= other

    def __gt__(self, other):
        return self._get() > other

    def __ge__(self, other):
        return self._get() >= other

    def __neg__(self):
        return -self._get()

    def __pos__(self):
        return +self._get()

    def __abs__(self):
        return abs(self._get())

    def __float__(self):
        return float(self._get())

    def __complex__(self):
        return complex(self._get())


def _binary(name, reflected):
    fn = getattr(operator, name)
    if reflected:
        def op(self, other):
            return fn(other, self._get())
    else:
        def op(self, other):
            return fn(self._get(), other)
    return op


for _name in ('add', 'sub', 'mul', 'truediv', 'floordiv', 'mod', 'pow'):
    setattr(Standin, f'__{_name}__', _binary(_name, False))
    setattr(Standin, f'__r{_name}__', _binary(_name, True))
del _name


class LazyModule(types.ModuleType):
    """The module `name`, imported on the first attribute it is asked for."""

    def __init__(self, name):
        super().__init__(name)
        self.__dict__['_lazy_obj'] = _UNSET

    def _get(self):
        obj = self.__dict__['_lazy_obj']
        if obj is _UNSET:
            import sympy  # noqa: F401 -- its import settles every stand-in
            settle()
            obj = self.__dict__['_lazy_obj']
            if obj is _UNSET:
                obj = _import(self.__name__)
        return obj

    def __getattr__(self, name):
        if name == '_lazy_obj':
            raise AttributeError(name)
        return getattr(self._get(), name)

    def __repr__(self):
        return f'<lazy module {self.__name__!r}>'


def _loaded():
    return 'sympy' in sys.modules and (_EXECUTED[0] or _ready('sympy'))


def module(name):
    """The module `name` (a sympy module), or a `LazyModule` for it while
    sympy is not imported."""
    if _loaded():
        return _import(name)
    with _LOCK:
        lm = LazyModule(name)
        _PENDING.append(lm)
    _watch()
    return lm


def standin(mod, attr, sig=None, call=False, after=None, cls_module=None):
    """`attr` of the module `mod` (its value called, with `call`), or a
    `Standin` for it while sympy is not imported -- or while `after`, a
    module the target's own import needs, is still importing.  `sig`: its
    cache-key text where it is not a class; a class -- defined in
    `cls_module`, where that is not `mod` -- is keyed as `_hdl_cache` keys a
    class, by module and name, and a pycircuit one by its module's source."""
    target = (mod, attr, call)
    if _loaded() and _ready(after):
        return _value_of(target)
    dep = None
    if sig is None:
        home = cls_module or mod
        sig = f'class:{home}.{attr}'
        if home.startswith('pycircuit'):
            dep = home
    with _LOCK:
        s = Standin(target, sig, dep, after)
        _PENDING.append(s)
    _watch()
    return s


def describe(val):
    """``(sig, dep)`` -- the cache key's text and dependency module -- for a
    stand-in, else None."""
    t = type(val)
    if t is Standin:
        return val._sig, val._dep
    if t is LazyModule:
        return f'module:{val.__name__}', None
    return None


def is_standin(val):
    return type(val) is Standin or type(val) is LazyModule


def resolved(val):
    """`val`, or the object it stands for."""
    if type(val) is Standin or type(val) is LazyModule:
        return val._get()
    return val


def settle(ready=()):
    """Resolve every pending stand-in whose module can be imported now and
    put its target in its place in every module's namespace -- where sympy
    is imported (else nothing).  `ready`: modules taken as imported (a
    module settling at its own end, its names all defined)."""
    if not _loaded():
        return
    with _LOCK:
        done = {}
        keep = []
        for s in _PENDING:
            if type(s) is LazyModule:
                s.__dict__['_lazy_obj'] = obj = _import(s.__name__)
            elif s._after in ready or _ready(s._after):
                s._obj = obj = _value_of(s._target)
            else:
                keep.append(s)
                continue
            done[id(s)] = (s, obj)
        _PENDING[:] = keep
        if not done:
            return
        for m in list(sys.modules.values()):
            d = getattr(m, '__dict__', None)
            if type(d) is not dict:
                continue
            swap = [(k, done[id(v)][1]) for k, v in list(d.items())
                    if (type(v) is Standin or type(v) is LazyModule) and id(v) in done]
            for k, obj in swap:
                d[k] = obj


def resolve_function(func, depth=0):
    """Replace the stand-ins in `func`'s closure cells and defaults -- and in
    those of the functions its cells hold -- by their targets."""
    if not isinstance(func, types.FunctionType) or depth > 3:
        return
    for cell in func.__closure__ or ():
        try:
            val = cell.cell_contents
        except ValueError:
            continue
        if is_standin(val):
            cell.cell_contents = resolved(val)
        elif isinstance(val, types.FunctionType):
            resolve_function(val, depth + 1)
    if func.__defaults__ and any(is_standin(v) for v in func.__defaults__):
        func.__defaults__ = tuple(resolved(v) for v in func.__defaults__)
    if func.__kwdefaults__ and any(is_standin(v) for v in func.__kwdefaults__.values()):
        func.__kwdefaults__ = {k: resolved(v) for k, v in func.__kwdefaults__.items()}


## -- sympy's import settles the stand-ins --------------------------------------

class _Then(importlib.abc.Loader):
    """sympy's own loader, then `settle`."""

    def __init__(self, loader):
        self._loader = loader

    def create_module(self, spec):
        return self._loader.create_module(spec)

    def exec_module(self, mod):
        loader = self._loader
        mod.__loader__ = loader
        if mod.__spec__ is not None:
            mod.__spec__.loader = loader
        loader.exec_module(mod)
        _EXECUTED[0] = True
        _unwatch()
        settle()


class _Watch(importlib.abc.MetaPathFinder):
    """On sympy's import, its spec with `_Then` for loader (once)."""

    def find_spec(self, fullname, path, target=None):
        if fullname != 'sympy':
            return None
        for finder in sys.meta_path:
            if finder is self or not hasattr(finder, 'find_spec'):
                continue
            spec = finder.find_spec(fullname, path, target)
            if spec is not None:
                break
        else:
            return None
        if spec.loader is not None and not isinstance(spec.loader, _Then):
            spec.loader = _Then(spec.loader)
        return spec


_WATCH = [None]


def _watch():
    with _LOCK:
        if _WATCH[0] is None and 'sympy' not in sys.modules:
            _WATCH[0] = w = _Watch()
            sys.meta_path.insert(0, w)


def _unwatch():
    with _LOCK:
        w = _WATCH[0]
        if w is not None:
            _WATCH[0] = None
            try:
                sys.meta_path.remove(w)
            except ValueError:
                pass


def load():
    """Import sympy, every stand-in settled: the compiler's entry, where the
    language's names must be its objects before any expression is built."""
    import sympy  # noqa: F401 -- its import settles the stand-ins
    settle()
