"""Read SPICE model cards -- enough of the format to ingest a real PDK.

A foundry model card is not a list of numbers.  The IHP PSP103 card is
359 parameters, most of them quoted expressions referencing ``.param``
multipliers that a *corner section* defines, `.include`d from inside a
``.subckt`` so the card can also see instance parameters like ``w``,
``ng`` and ``pre_layout``.  None of it resolves without following that
whole chain.

So this is not a netlist parser and does not try to be one.  It reads
the declarative part -- ``.lib`` sections, ``.include``, ``.param``,
``.subckt``, ``.model`` -- and answers one question:

    what are the numeric parameters of model X, in corner Y, for an
    instance with geometry Z?

Everything else (device lines, analyses, control blocks) is recognised
well enough to be skipped.

Usage::

    deck = spicecard.read('cornerMOSlv.lib', section='mos_tt')
    p = deck.model_params('sg13g2_lv_nmos_psp', w=1e-6, l=0.13e-6, ng=1)

**Statistical functions return their nominal value.**  ``agauss``,
``gauss``, ``aunif`` and ``unif`` describe a distribution; without a
random draw the only defensible answer is the centre, which is also what
a nominal-corner simulation uses.  Monte-Carlo would need a generator
threaded through, and is not implemented rather than faked.
"""

import math
import os
import re


#: SPICE engineering suffixes.  ``meg`` must be tried before ``m``.
_SUFFIX = [('meg', 1e6), ('mil', 25.4e-6), ('t', 1e12), ('g', 1e9),
           ('k', 1e3), ('m', 1e-3), ('u', 1e-6), ('n', 1e-9),
           ('p', 1e-12), ('f', 1e-15), ('a', 1e-18)]

_NUM = re.compile(
    r'\b(\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?)'
    r'(meg|mil|t|g|k|m|u|n|p|f|a)?([a-zA-Z_]*)\b')

#: What an expression may call.  Deliberately closed: an expression is
#: data from a vendor file, and `eval` over an open namespace would make
#: reading a model card equivalent to running it.
_FUNCS = {
    'sqrt': math.sqrt, 'exp': math.exp, 'ln': math.log,
    'log': math.log10, 'log10': math.log10, 'abs': abs,
    'pow': lambda a, b: a ** b, 'max': max, 'min': min,
    'int': int, 'floor': math.floor, 'ceil': math.ceil,
    'sin': math.sin, 'cos': math.cos, 'tan': math.tan,
    'atan': math.atan, 'sinh': math.sinh, 'cosh': math.cosh,
    'tanh': math.tanh, 'sgn': lambda x: (x > 0) - (x < 0),
    ## Statistical: nominal value, see the module docstring.
    'agauss': lambda nom, *a: nom,
    'gauss': lambda nom, *a: nom,
    'aunif': lambda nom, *a: nom,
    'unif': lambda nom, *a: nom,
    'limit': lambda x, lo, hi: max(lo, min(hi, x)),
    'if': lambda c, a, b: a if c else b,
}


class SpiceCardError(Exception):
    """Anything wrong with a card: syntax, a missing name, a cycle."""


def _strip_comments(line):
    """Drop `*` full-line comments and `;` / `$` trailing ones.

    Quote-aware, because a `;` can legitimately appear inside a quoted
    expression.
    """
    s = line.rstrip('\n')
    if s.lstrip().startswith('*'):
        return ''
    out, quote = [], None
    for ch in s:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in "'\"":
            quote = ch
            out.append(ch)
            continue
        if ch in ';$':
            break
        out.append(ch)
    return ''.join(out)


def _logical_lines(path, seen=None):
    """Yield (text, dirname) with continuations joined and comments gone.

    `.include` is followed here rather than later, so the caller sees one
    flat stream; `dirname` travels with each line because an include path
    is relative to the file that names it.
    """
    seen = seen if seen is not None else set()
    real = os.path.realpath(path)
    if real in seen:
        raise SpiceCardError('circular .include of %s' % path)
    if not os.path.exists(path):
        raise SpiceCardError('no such file: %s' % path)

    here = os.path.dirname(os.path.abspath(path))
    with open(path, errors='replace') as fh:
        raw = fh.readlines()

    pending = None
    for line in raw + ['\n']:
        text = _strip_comments(line)
        stripped = text.strip()
        if not stripped:
            ## A comment or a blank line does NOT end a continued
            ## statement -- real cards put commentary between `+` lines,
            ## and flushing here silently truncated the card.
            continue
        if stripped.startswith('+'):
            if pending is None:
                raise SpiceCardError('continuation with nothing to continue '
                                     'in %s: %r' % (path, line))
            pending += ' ' + stripped[1:]
            continue
        if pending is not None and pending.strip():
            yield pending.strip(), here
        pending = text
    if pending is not None and pending.strip():
        yield pending.strip(), here


_ASSIGN = re.compile(r"([A-Za-z_][\w.\[\]]*)\s*=\s*"
                     r"('[^']*'|\"[^\"]*\"|\{[^}]*\}|[^\s=]+)")


def _assignments(text):
    """`a=1 b='x*2' c={y}` -> [(name, raw_value), ...], order preserved."""
    return [(m.group(1).lower(), m.group(2)) for m in _ASSIGN.finditer(text)]


class _Scope(object):
    """A parameter namespace: global, or one `.subckt`."""

    def __init__(self, name, parent):
        self.name = name
        self.parent = parent
        #: name -> list of (condition_or_None, raw_expression), in file
        #: order.  A list because `.if` can define one name several ways.
        self.params = {}

    def define(self, name, raw, cond):
        self.params.setdefault(name, []).append((cond, raw))

    def chain(self):
        s, out = self, []
        while s is not None:
            out.append(s)
            s = s.parent
        return out


class Model(object):
    """One `.model` card: its type, its raw parameters, and its scope."""

    def __init__(self, name, mtype, scope):
        self.name = name
        self.type = mtype
        self.scope = scope
        self.raw = {}

    def __repr__(self):
        return 'Model(%r, %r, %d params)' % (self.name, self.type,
                                             len(self.raw))


class Deck(object):
    """A parsed deck: parameters, models and subcircuit definitions."""

    def __init__(self):
        self.global_scope = _Scope(None, None)
        self.models = {}
        self.subckt_ports = {}

    ## ---------------------------------------------------------------
    ## Expression evaluation
    ## ---------------------------------------------------------------

    @staticmethod
    def _pythonise(expr):
        """Rewrite a SPICE expression as Python.

        Suffix handling is the fiddly part: `1u` is 1e-6, but `1e-6` must
        not have its `e` eaten, and a bare identifier like `ns1` is a
        name rather than a number with a suffix.
        """
        e = expr.strip()
        if e[:1] in "'\"" and e[-1:] == e[:1]:
            e = e[1:-1]
        elif e.startswith('{') and e.endswith('}'):
            e = e[1:-1]
        e = e.replace('^', '**')

        def num(m):
            base, suf, tail = m.group(1), m.group(2), m.group(3)
            if tail:
                ## `1nsomething` -- not a suffixed number, leave alone.
                return m.group(0)
            if not suf:
                return base
            for name, mult in _SUFFIX:
                if suf.lower() == name:
                    return '(%s*%g)' % (base, mult)
            return m.group(0)

        return _NUM.sub(num, e)

    def _evaluate(self, raw, resolver):
        code = self._pythonise(raw)
        try:
            return float(eval(code, {'__builtins__': {}}, resolver))
        except SpiceCardError:
            raise
        except ZeroDivisionError:
            raise SpiceCardError('division by zero evaluating %r' % raw)
        except Exception as exc:
            raise SpiceCardError('cannot evaluate %r (as %r): %s'
                                 % (raw, code, exc))

    ## ---------------------------------------------------------------
    ## Resolution
    ## ---------------------------------------------------------------

    def _resolver(self, scope, overrides):
        """A mapping that resolves names lazily, with cycle detection.

        Lookup order: caller overrides, then the scope chain innermost
        first, then the function table.  Overrides win because that is
        what an instance parameter IS -- the subckt's `.param w=0.5u` is
        a default the instantiation replaces.
        """
        deck = self
        cache = dict(overrides)
        active = set()

        class _Res(dict):
            def __missing__(self, key):
                k = key.lower()
                if k in cache:
                    return cache[k]
                if k in _FUNCS:
                    return _FUNCS[k]
                if k in active:
                    raise SpiceCardError(
                        'circular parameter definition through %r' % k)
                for sc in scope.chain():
                    if k not in sc.params:
                        continue
                    active.add(k)
                    try:
                        for cond, raw in sc.params[k]:
                            if cond is not None and not deck._evaluate(
                                    cond, self):
                                continue
                            val = deck._evaluate(raw, self)
                            cache[k] = val
                            return val
                    finally:
                        active.discard(k)
                raise SpiceCardError('undefined parameter %r' % key)

        return _Res()

    def model_params(self, name, **overrides):
        """Resolved numeric parameters of one `.model`.

        `overrides` supply whatever the card reads from outside itself --
        instance geometry, subcircuit flags -- by name, case-insensitive.
        """
        key = name.lower()
        if key not in self.models:
            raise SpiceCardError(
                'no model %r in this deck (have: %s)'
                % (name, ', '.join(sorted(self.models)) or 'none'))
        model = self.models[key]
        res = self._resolver(model.scope,
                             {k.lower(): v for k, v in overrides.items()})
        out = {}
        for pname, raw in model.raw.items():
            out[pname] = self._evaluate(raw, res)
        return out

    def param(self, name, **overrides):
        """Resolve one global `.param`."""
        res = self._resolver(self.global_scope,
                             {k.lower(): v for k, v in overrides.items()})
        return res[name.lower()]


def read(path, section=None):
    """Parse `path`, entering `.LIB section` if one is named.

    A `.LIB`/`.ENDL` block is opt-in: without `section` every one of them
    is skipped, which is what makes a corner file safe to read (its
    sections define the same names differently, and concatenating them
    would give whichever came last).
    """
    deck = Deck()
    _parse_file(deck, path, section, deck.global_scope, set())
    return deck


def _parse_file(deck, path, section, scope, seen):
    if not os.path.exists(path):
        raise SpiceCardError('no such file: %s' % path)
    real = os.path.realpath(path)
    if real in seen:
        raise SpiceCardError('circular .include of %s' % path)
    seen = set(seen) | {real}
    ## `.LIB` blocks we are not collecting; nonzero means skipping.
    lib_skip = 0
    in_wanted_lib = False
    ## `.if` nesting: a stack of conditions, or None once a branch of the
    ## chain has been taken.
    cond_stack = []
    scopes = [scope]
    model = None

    for text, here in _logical_lines(path, set()):
        low = text.lower()
        head = low.split(None, 1)[0] if low.split() else ''

        if head in ('.lib', '.endl'):
            parts = text.split()
            if head == '.endl':
                if in_wanted_lib and lib_skip == 0:
                    in_wanted_lib = False
                elif lib_skip:
                    lib_skip -= 1
                continue
            ## `.lib <file> <section>` is an include; `.lib <section>`
            ## opens a block.
            if len(parts) >= 3:
                if lib_skip == 0:
                    _parse_file(deck, os.path.join(here, parts[1]),
                                parts[2], scopes[-1], seen)
                continue
            want = parts[1].lower() if len(parts) > 1 else None
            if lib_skip == 0 and section is not None and want == section.lower():
                in_wanted_lib = True
            else:
                lib_skip += 1
            continue
        if lib_skip:
            continue
        if section is not None and not in_wanted_lib and head not in (
                '.include', '.inc'):
            ## Outside the requested section: only follow includes that
            ## sit at file level, nothing else counts.
            if head not in ('.subckt', '.ends', '.model', '.param', '.if',
                            '.else', '.elseif', '.endif'):
                continue

        if head in ('.include', '.inc'):
            parts = text.split()
            if len(parts) >= 2:
                _parse_file(deck, os.path.join(here, parts[1].strip('"\'')),
                            None, scopes[-1], seen)
            continue

        if head == '.if':
            cond_stack.append(text[text.find('(') + 1:text.rfind(')')])
            continue
        if head == '.elseif':
            if cond_stack:
                cond_stack[-1] = text[text.find('(') + 1:text.rfind(')')]
            continue
        if head == '.else':
            if cond_stack:
                cond_stack[-1] = 'not (%s)' % cond_stack[-1]
            continue
        if head == '.endif':
            if cond_stack:
                cond_stack.pop()
            continue

        cond = ' and '.join('(%s)' % c for c in cond_stack) or None

        if head in ('.param', '.params'):
            for pname, raw in _assignments(text[len(head):]):
                scopes[-1].define(pname, raw, cond)
            continue

        if head == '.subckt':
            parts = text.split()
            if len(parts) < 2:
                raise SpiceCardError('malformed .subckt: %r' % text)
            sub = _Scope(parts[1].lower(), scopes[-1])
            ports = [p for p in parts[2:] if '=' not in p]
            deck.subckt_ports[parts[1].lower()] = ports
            for pname, raw in _assignments(text):
                sub.define(pname, raw, None)
            scopes.append(sub)
            continue

        if head == '.ends':
            if len(scopes) > 1:
                scopes.pop()
            continue

        if head == '.model':
            parts = text.split()
            if len(parts) < 3:
                raise SpiceCardError('malformed .model: %r' % text[:80])
            model = Model(parts[1].lower(), parts[2].lower(), scopes[-1])
            deck.models[model.name] = model
            ## Parameters may follow on the same logical line.
            after = text.split(None, 3)
            if len(after) > 3:
                for pname, raw in _assignments(after[3]):
                    model.raw[pname] = raw
            continue

        ## Anything else -- device lines, analyses, control blocks -- is
        ## not this reader's business.
        model = None
