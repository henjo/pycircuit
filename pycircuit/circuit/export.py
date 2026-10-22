"""A flat circuit written out as standalone Python: the source of a module
whose function rebuilds it, element for element, as the solvers see it.

    >>> from pycircuit.circuit import SubCircuit, gnd
    >>> from pycircuit.circuit.elements import R, VS
    >>> from pycircuit.circuit.export import to_python, rebuild, selfcheck
    >>> c = SubCircuit()
    >>> c['V1'] = VS('in', gnd, v=1.5)
    >>> c['R1'] = R('in', gnd, r=2e3)
    >>> src = to_python(c, 'divider')
    >>> c2 = rebuild(src, 'divider')
    >>> [str(n) for n in c2.nodes] == [str(n) for n in c.nodes]
    True
    >>> selfcheck(c, src, name='divider') is None
    True

WHAT A REBUILD MUST KEEP, and how (each is a way a hand-written copy of a
circuit took other solver paths than the original):

* THE NODE ORDER.  The unknowns' order steers the solver -- pivots, PCNR's
  routes, the paths from zeros -- and a start vector is in that order: the
  nodes are created in `cir.nodes` order first, gnd in its slot, before any
  element adds them in first-use order.
* GIVENNESS.  A parameter left at its default is not a parameter given at
  its default value: a model reads `is_given` (an HDL class packs it with
  its parameters -- `hdl._params_of` -- and SPICE's level 3 MOS computes
  PHI from NSUB only where PHI is not given).  Only the given parameters
  are written, raw from `ipar` (a string expression stays one).
* THE CLASS.  A collapse variant (`X_collapseNN`, `hdl._collapse_variant`)
  is written as its base: the parameters collapse it again.  A
  `BSourceHdl` class (one a traced expression) is written as
  `BSourceHdl(i_func=..., q_func=...)` from `bsource_exprs` by
  `sympy.srepr` -- exact, and the class is keyed by the expression.  A
  class outside the library is written from its source
  (`inspect.getsource`); one local to a function, or with no source, is
  refused.
* THE NAMES.  A parameter that is a Python keyword (MOS's `as`) goes
  through `**{...}`.
* THE VALUES.  Python floats by `repr` (an exact round trip), ints, bools,
  strings, None, tuples and lists of those, numpy arrays with their dtype;
  anything else -- a callable, an object -- is refused.

`to_python` then CHECKS what it wrote (`check=True`): each element rebuilt
beside the original, i/q/G/C at random x and u (dc, tran and ac),
`next_event` and the noise matrix CY compared bit for bit, the stateful
limiters' analysis state (elements.Diode's `_vlim`) cleared on both sides
for the comparison and put back after -- so state the parameters do not
hold (a source whose time function was replaced after construction, a
value written to `iparv` alone) raises `ExportError` naming the element,
never a silently different circuit.  `selfcheck` compares the whole
rebuilt circuit the same way.

Not exported (an `ExportError` names the element and why): hierarchy (a
subcircuit instance, a SubCircuit subclass, a circuit with parameters of
its own), an element on another toolkit, a callable or object parameter
(`elements.BSource`'s functions: `BSourceHdl` is the exportable form).
"""
import inspect
import keyword
import math
import textwrap

import numpy as np
import sympy

#: the module the exported source is executed as (`rebuild`)
MODULE_NAME = 'pycircuit_export'

__all__ = ['ExportError', 'rebuild', 'selfcheck', 'to_python', 'x0_to_python']


class ExportError(TypeError):
    """What the exporter cannot write, or wrote differently: the element,
    its class and why."""


def to_python(cir, name='circuit', check=True, allow_source=True):
    """The source of a module defining ``def <name>():``, which returns the
    rebuilt flat circuit `cir` (the module's notes say what is kept and
    how).  `check`: each element rebuilt beside the original and compared
    (`ExportError` where one differs).  `allow_source`: a class outside the
    library is written from its source (else refused)."""
    src = _to_python(cir, name, allow_source)
    if check:
        _check_elements(cir, src, name)
    return src


def rebuild(src, name='circuit'):
    """The circuit the exported source `src` builds: its module executed,
    ``<name>()`` called."""
    ns = {'__name__': MODULE_NAME}
    exec(compile(src, f'<export {name}>', 'exec'), ns)  # noqa: S102  (the exported module)
    return ns[name]()


def selfcheck(cir, src=None, n=5, seed=0, epar=None, name='circuit', scale=1.0):
    """None where the rebuilt circuit is `cir` as the solvers see it, else
    the first difference (a string): the node, branch and element order,
    each element's class, then i/q/G/C at `n` random x (uniform in
    [-scale, scale]; CY at the first), u (dc, tran, ac) and `next_event` at
    `n` times,
    bit for bit (a NaN equal to a NaN, a signed zero its sign).  `src`
    defaults to `to_python(cir, name)`."""
    from pycircuit.circuit import _limiting
    from pycircuit.circuit.circuit import defaultepar
    epar = defaultepar if epar is None else epar
    if src is None:
        src = to_python(cir, name)
    new = rebuild(src, name)
    if _node_key(new.nodes) != _node_key(cir.nodes):
        return 'node order'
    if _branch_key(new.branches) != _branch_key(cir.branches):
        return 'branch order'
    if list(new.elements) != list(cir.elements):
        return 'element order'
    for k, e in cir.elements.items():
        if not _same_class(type(e), type(new.elements[k])):
            return f'element {k!r}: class {type(new.elements[k]).__name__}, not {type(e).__name__}'
    snap = _cleared(cir, new)
    try:
        return _first_difference(cir, new, cir.n, n, seed, epar, scale)
    finally:
        _limiting.state_restore(snap)


def x0_to_python(x, name='x0'):
    """The source of ``<name> = np.array([...])`` for the vector `x`
    (float64: every entry by `repr`, exact; a NaN, an infinity, a signed
    zero kept)."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ExportError(f'a start vector of shape {x.shape}')
    return f'{name} = np.array({_value(x.tolist(), name)})\n'


def _to_python(cir, name, allow_source, given_only=True, node_order=True):
    """`to_python`'s text, unchecked.  `given_only`, `node_order` False:
    every parameter written / the nodes left to first use (what a hand copy
    does; the tests show the check catches both)."""
    from pycircuit.circuit import circuit as _circuit
    from pycircuit.circuit.circuit import SubCircuit
    if type(cir) is not SubCircuit:
        raise ExportError(f'the circuit is a {type(cir).__name__}, not a SubCircuit: hierarchy '
                          'is not exported')
    if cir.instparams:
        raise ExportError('the circuit has instance parameters of its own')
    if cir.toolkit is not _circuit.default_toolkit:
        raise ExportError(f'the circuit\'s toolkit is {cir.toolkit!r}, not the default its '
                          'rebuild would take')
    user, lines, bs_used = {}, [f'def {name}():', '    c = SubCircuit()'], False
    if node_order:
        nodes = [(str(nd.name), bool(nd.isglobal)) for nd in cir.nodes]
        lines.append(f'    for _name, _isglobal in {nodes!r}:')
        lines.append('        c.append_node(N(_name, _isglobal))')
    for iname, e in cir.elements.items():
        where = f'element {iname!r} ({type(e).__name__})'
        if isinstance(e, SubCircuit):
            raise ExportError(f'{where}: a subcircuit instance -- hierarchy is not exported')
        if e.toolkit is not cir.toolkit:
            raise ExportError(f'{where}: a toolkit of its own')
        tmap = cir.term_node_map[iname]
        args = [f'N({str(tmap[t].name)!r}, {bool(tmap[t].isglobal)!r})' for t in e.terminals]
        t = type(e)
        exprs = getattr(t, 'bsource_exprs', None)
        if exprs is not None:
            bs_used = True
            args += _bsource_args(exprs, where)
            ctor = 'BSourceHdl'
        else:
            t = getattr(t, '_hdl_collapse_base', t)
            ctor = _library_path(t)
            if ctor is None:
                if t.__qualname__ not in user:
                    user[t.__qualname__] = _class_source(t, where, allow_source)
                ctor = t.__qualname__
            args += _param_args(e, where, given_only)
        lines.append(f'    c[{iname!r}] = {ctor}({", ".join(args)})')
    lines.append('    return c')
    return '\n'.join(_head(bs_used, user) + lines) + '\n'


def _head(bs_used, user):
    """The module's imports and helpers: `N` (a node, gnd itself), `_bs` (a
    BSourceHdl function from its traced expression), the user classes."""
    head = ['"""A circuit written by pycircuit.circuit.export.to_python."""',
            'import numpy as np  # noqa: F401  (array parameters)',
            'import sympy  # noqa: F401  (BSourceHdl expressions)',
            '',
            'from pycircuit.circuit import elements, elements_hdl  # noqa: F401',
            'from pycircuit.circuit.circuit import Node, SubCircuit, gnd',
            'from pycircuit.circuit.elements_hdl import BSourceHdl  # noqa: F401',
            '', '',
            'def N(name, isglobal):',
            '    """A node by name -- gnd itself where it is gnd."""',
            ('    return gnd if (name, isglobal) == (gnd.name, True) else Node(name, '
             'isglobal=isglobal)'),
            '']
    if bs_used:
        head += ['',
                 'def _bs(expr, sym):',
                 '    """A BSourceHdl function from its traced expression (sympy srepr)."""',
                 '    e = eval(expr, vars(sympy))',
                 '    s = eval(sym, vars(sympy)) if sym is not None else None',
                 '    return (lambda v: e.subs(s, v)) if s is not None else (lambda v: e)',
                 '']
    if user:
        head += ['', 'from pycircuit.circuit.hdl import *  # noqa: F401,F403  (the classes below)',
                 'from pycircuit.utilities.param import Parameter  # noqa: F401', '']
        for src in user.values():
            head += ['', src.rstrip(), '']
    return head + ['', '']


def _library_path(t):
    """`elements.X` / `elements_hdl.X` where `t` is that library class."""
    from pycircuit.circuit import elements, elements_hdl
    for modname, mod in (('elements', elements), ('elements_hdl', elements_hdl)):
        if getattr(mod, t.__name__, None) is t:
            return f'{modname}.{t.__name__}'
    return None


def _class_source(t, where, allow_source):
    if not allow_source:
        raise ExportError(f'{where}: class {t.__module__}.{t.__qualname__} is not pycircuit\'s')
    if '<locals>' in t.__qualname__:
        raise ExportError(f'{where}: class {t.__qualname__} is local to a function')
    try:
        src = inspect.getsource(t)
    except (OSError, TypeError) as e:
        raise ExportError(f'{where}: class {t.__qualname__} has no source ({e})') from e
    return textwrap.dedent(src)


def _param_args(e, where, given_only):
    """The constructor's parameter arguments: the given ones (all where not
    `given_only`), raw from `ipar`; a keyword or non-identifier through
    `**{...}`."""
    names = list(dict.fromkeys(p.name for p in e.instparams))
    if given_only:
        names = [p for p in names if e.ipar.is_given(p)]
    plain, kws = [], []
    for p in names:
        v = _value(e.ipar.get(p), f'{where}, parameter {p!r}')
        (kws if keyword.iskeyword(p) or not p.isidentifier() else plain).append((p, v))
    out = [f'{p}={v}' for p, v in plain]
    if kws:
        out.append('**{' + ', '.join(f'{p!r}: {v}' for p, v in kws) + '}')
    return out


def _bsource_args(exprs, where):
    out = []
    for kw, expr in zip(('i_func', 'q_func'), exprs):
        if expr == 0:
            continue
        syms = sorted(expr.free_symbols, key=str)
        if len(syms) > 1:
            raise ExportError(f'{where}: a BSourceHdl expression in {len(syms)} symbols')
        sym = sympy.srepr(syms[0]) if syms else None
        out.append(f'{kw}=_bs({sympy.srepr(expr)!r}, {sym!r})')
    return out


def _value(v, where):
    """The exact source of a parameter value, or ExportError."""
    if v is None or isinstance(v, (bool, str)):
        return repr(v)
    if isinstance(v, (int, np.integer)):
        return repr(int(v))
    if isinstance(v, (float, np.floating)):
        f = float(v)
        if math.isnan(f):
            return "float('nan')"
        if math.isinf(f):
            return "float('inf')" if f > 0 else "float('-inf')"
        return repr(f)
    if isinstance(v, np.ndarray):
        if v.dtype.kind not in 'biuf':
            raise ExportError(f'{where}: an array of dtype {v.dtype}')
        return f'np.array({_value(v.tolist(), where)}, dtype={v.dtype.str!r})'
    if isinstance(v, (list, tuple)):
        body = ', '.join(_value(x, where) for x in v)
        if isinstance(v, tuple):
            return '(' + body + (',' if len(v) == 1 else '') + ')'
        return '[' + body + ']'
    raise ExportError(f'{where}: a value of type {type(v).__name__} ({v!r:.60}) cannot be '
                      'written')


def _check_elements(cir, src, name, n=3, seed=1):
    """Each element rebuilt beside the original and compared (`_differs`):
    ExportError naming the first that differs."""
    from pycircuit.circuit import _limiting
    new = rebuild(src, name)
    snap = _cleared(cir, new)
    try:
        rng = np.random.default_rng(seed)
        for k, e in cir.elements.items():
            f = new.elements[k]
            where = f'element {k!r} ({type(e).__name__})'
            if not _same_class(type(e), type(f)):
                raise ExportError(f'{where}: rebuilt as {type(f).__name__}')
            why = _first_difference(e, f, e.n, n, rng, None, 1.0)
            if why is not None:
                raise ExportError(f'{where}: its {why} differs once rebuilt (state its '
                                  'parameters do not hold)')
    finally:
        _limiting.state_restore(snap)


def _first_difference(a, b, size, n, seed, epar, scale):
    """The first of i/q/G/C at `n` random x (CY at the first), u (dc, tran,
    ac) and next_event at `n` times on which the circuits (or elements) `a`
    and `b` differ -- a value or the exception raised -- else None."""
    from pycircuit.circuit.circuit import defaultepar
    epar = defaultepar if epar is None else epar
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    for j in range(n):
        x = rng.uniform(-1.0, 1.0, size) * scale
        for f in ('i', 'q', 'G', 'C'):
            if not _same(_outcome(getattr(a, f), x, epar), _outcome(getattr(b, f), x, epar)):
                return f'{f} at x #{j}'
        ## (the noise at the first x alone: an HDL class's CY is Python's, and
        ## one point shows a parameter it reads)
        if j == 0 and not _same(_outcome(a.CY, x, 2e3 * np.pi, epar),
                                _outcome(b.CY, x, 2e3 * np.pi, epar)):
            return 'CY at x #0'
    for j in range(n):
        t = float(rng.uniform(0.0, 1e-3)) if j else 0.0
        for an in ('dc', 'tran', 'ac'):
            if not _same(_outcome(a.u, t, epar, an), _outcome(b.u, t, epar, an)):
                return f'u({an}) at t #{j} = {t!r}'
        if not _same(_outcome(a.next_event, t), _outcome(b.next_event, t)):
            return f'next_event at t #{j} = {t!r}'
    return None


def _node_key(nodes):
    return [(str(nd.name), bool(nd.isglobal)) for nd in nodes]


def _branch_key(branches):
    return [(str(b.plus.name), str(b.minus.name)) for b in branches]


def _same_class(a, b):
    """The same class -- or, for one written from its source, the same
    qualified name and collapse mask (the rebuilt one is a new class
    object; the numbers decide the rest)."""
    if a is b:
        return True
    if _library_path(getattr(a, '_hdl_collapse_base', a)) is not None:
        return False
    return a.__qualname__ == b.__qualname__ and (getattr(a, '_hdl_collapse_mask', None)
                                                 == getattr(b, '_hdl_collapse_mask', None))


def _outcome(f, *args):
    """A call's value, or the exception it raised (its type and message):
    both sides must agree on either."""
    try:
        return f(*args)
    except Exception as e:                                      # noqa: BLE001
        return ('raised', type(e).__name__, str(e))


def _cleared(*circuits):
    """The stateful limiters' analysis state (elements.Diode's `_vlim`: the
    last solve's, not the circuit's) snapshotted, then cleared -- the
    comparison starts from no state on both sides; `state_restore` puts it
    back."""
    from pycircuit.circuit import _limiting
    lims = []
    for c in circuits:
        if hasattr(c, 'elements'):
            lims += _limiting.stateful_limiters(c)
    snap = _limiting.limiter_snapshot(lims)
    for e in lims:
        e.reset_limit_state()
    return snap


def _same(a, b):
    if (isinstance(a, tuple) and a[:1] == ('raised',)) or (isinstance(b, tuple)
                                                         and b[:1] == ('raised',)):
        return a == b
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or a.dtype != b.dtype:
        return False
    if a.dtype.kind in 'fc':
        return bool(np.array_equal(a, b, equal_nan=True)) and bool(
            np.array_equal(np.signbit(a.real), np.signbit(b.real)))
    return bool(np.array_equal(a, b))
