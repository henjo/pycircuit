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

**The lexical rules are SPICE's** (2026-10-07), for netlists as for
cards: a scale factor in any case (``1MEG``, ``1.5P``; ``meg`` and
``mil`` before ``m``, which is milli), and the letters after it a unit,
ignored (``10pF``, ``1.8mA``, ``5V`` -- so ``1F`` is a femto and
``1Mohm`` a milliohm, as in SPICE); a ``.model`` card's parameters may sit
in parentheses; a ``.subckt`` header's ports end at ``params:`` or its
first assignment; a ``$`` starts a comment only where it starts a token
(``v$d5`` is a name).  `number` reads one value.

**Expressions are parsed, not rewritten**: SPICE's tokens, Python's
operators with Python's precedence -- what the regex translation they
replaced handed to Python, so an expression it read means what it meant
(pinned on the IHP cards) -- and SPICE's ``c ? a : b``, ``&&``, ``||``,
``!`` and ``<>`` with C's.  Only numbers, names, calls and operators
parse: no attribute, subscript, string or lambda reaches `eval`.
"""

import functools
import keyword
import math
import os
import re

#: SPICE scale factors, matched in any case at the start of the letters
#: after a number.  ``meg`` and ``mil`` must be tried before ``m``.
_SUFFIX = [('meg', 1e6), ('mil', 25.4e-6), ('t', 1e12), ('g', 1e9),
           ('k', 1e3), ('m', 1e-3), ('u', 1e-6), ('n', 1e-9),
           ('p', 1e-12), ('f', 1e-15), ('a', 1e-18)]

#: A number: its digits (an exponent included), then its letters -- a
#: scale factor, a unit, or both.
_DIGITS = r'(\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?)([A-Za-z_]*)'
_NUMBER = re.compile(r'\s*([-+]?)' + _DIGITS + r'\s*\Z')

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


def _scale(letters):
    """The multiplier a number's `letters` begin with, or None (a unit)."""
    low = letters.lower()
    for name, mult in _SUFFIX:
        if low.startswith(name):
            return mult
    return None


def _number_text(digits, letters):
    """A number as Python text: its digits (an integer's leading zeros
    dropped -- Python refuses them), times its scale factor written as the
    regex translation wrote it; the unit gone."""
    if digits.isdigit():
        digits = str(int(digits))
    mult = _scale(letters)
    return digits if mult is None else f'({digits}*{mult:g})'


def _literal(text):
    return int(text) if text.isdigit() else float(text)


def number(text):
    """One SPICE value as a float: ``1MEG``, ``10pF``, ``-2.5e-3``, ``5V``.
    The value its text has in an expression (the digits times the scale
    factor, in Python's arithmetic)."""
    m = _NUMBER.match(text)
    if m is None:
        raise SpiceCardError(f'not a number: {text!r}')
    sign, digits, letters = m.groups()
    if digits.isdigit():
        digits = str(int(digits))
    v = _literal(digits)
    mult = _scale(letters)
    if mult is not None:
        v = v * _literal(f'{mult:g}')
    return float(-v if sign == '-' else v)


## -------------------------------------------------------------------
## Expressions: tokens, a Pratt parser, Python text
## -------------------------------------------------------------------

_TOKEN = re.compile(r"""
    (?P<num>\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?)(?P<letters>[A-Za-z_]*)
  | (?P<name>[A-Za-z_]\w*)
  | (?P<op>\*\*|==|!=|<>|<=|>=|&&|\|\||[-+*/%^()<>!?:,])
""", re.VERBOSE)

#: Python's words an expression may use: its boolean operators (an `.if`
#: condition is built with them) and its constants.
_WORD_OPS = ('and', 'or', 'not')
_WORD_OK = _WORD_OPS + ('True', 'False', 'None')
#: A name that is another of Python's words -- `as`, `is`, `lambda` (a
#: source area, a saturation current, a channel-length modulation), or
#: SPICE's `if(c, a, b)` -- reaches Python under this prefix, which the
#: resolver strips.
_KW = '_kw_'

#: Infix operators: left binding power (higher binds tighter).  Python's
#: precedence for Python's operators; SPICE's `?:`, `&&`, `||` with C's.
_LBP = {'?': 10, 'or': 20, '||': 20, 'and': 30, '&&': 30,
        '<': 50, '>': 50, '<=': 50, '>=': 50, '==': 50, '!=': 50, '<>': 50,
        '+': 60, '-': 60, '*': 70, '/': 70, '%': 70, '**': 90, '^': 90}
_CMP = ('<', '>', '<=', '>=', '==', '!=', '<>')
#: Prefix operators: the binding power of their operand.  Python's `not`
#: is looser than a comparison; C's `!` is a unary operator.  A unary
#: minus is looser than `**` on its right (`-2**2` is -4, as in Python).
_PREFIX = {'-': 80, '+': 80, '!': 80, 'not': 40}


def _tokens(text):
    out, pos, n = [], 0, len(text)
    while True:
        while pos < n and text[pos].isspace():
            pos += 1
        if pos == n:
            return out
        m = _TOKEN.match(text, pos)
        if m is None:
            raise SpiceCardError(f'cannot parse {text!r}: unexpected {text[pos:pos + 10]!r}')
        if m.group('num') is not None:
            out.append(('num', _number_text(m.group('num'), m.group('letters'))))
        elif m.group('name') is not None:
            name = m.group('name')
            if name in _WORD_OPS:
                out.append(('op', name))
            elif keyword.iskeyword(name) and name not in _WORD_OK:
                out.append(('name', _KW + name))
            else:
                out.append(('name', name))
        else:
            out.append(('op', m.group('op')))
        pos = m.end()


class _Parser:
    """`text` to a tree: ('num' | 'name', text), ('call', name, args),
    ('group', node), ('unary', op, node), ('bin', op, a, b), ('cmp',
    operands, ops) (a chain, as Python's), ('bool', 'and' | 'or',
    operands) (flattened as Python flattens an unparenthesised chain),
    ('ifexp', c, a, b)."""

    def __init__(self, text):
        self.text = text
        self.toks = _tokens(text)
        self.i = 0

    def fail(self, what):
        raise SpiceCardError(f'cannot parse {self.text!r}: {what}')

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else (None, None)

    def take(self):
        tok = self.peek()
        self.i += 1
        return tok

    def expect(self, op):
        if self.take() != ('op', op):
            self.fail(f'expected {op!r}')

    def parse(self):
        if not self.toks:
            self.fail('empty')
        node = self.expr(0)
        if self.i != len(self.toks):
            self.fail(f'unexpected {self.peek()[1]!r}')
        return node

    def expr(self, rbp):
        left = self.prefix()
        while True:
            kind, op = self.peek()
            if kind != 'op' or _LBP.get(op, 0) <= rbp:
                return left
            self.take()
            left = self.infix(left, op)

    def prefix(self):
        kind, v = self.take()
        if kind == 'num':
            return ('num', v)
        if kind == 'name':
            if self.peek() != ('op', '('):
                return ('name', v)
            self.take()
            args = []
            if self.peek() != ('op', ')'):
                args.append(self.expr(0))
                while self.peek() == ('op', ','):
                    self.take()
                    args.append(self.expr(0))
            self.expect(')')
            return ('call', v, args)
        if kind == 'op' and v == '(':
            inner = self.expr(0)
            self.expect(')')
            return ('group', inner)
        if kind == 'op' and v in _PREFIX:
            return ('unary', v, self.expr(_PREFIX[v]))
        self.fail(f'unexpected {v!r}' if kind else 'unexpected end')

    def infix(self, left, op):
        if op == '?':
            a = self.expr(0)
            self.expect(':')
            return ('ifexp', left, a, self.expr(_LBP['?'] - 1))
        if op in _CMP:
            items, ops = [left, self.expr(_LBP[op])], [op]
            while self.peek()[0] == 'op' and self.peek()[1] in _CMP:
                ops.append(self.take()[1])
                items.append(self.expr(_LBP[op]))
            return ('cmp', items, ops)
        if op in ('and', '&&', 'or', '||'):
            py = 'and' if op in ('and', '&&') else 'or'
            right = self.expr(_LBP[op])
            if left[0] == 'bool' and left[1] == py:
                return ('bool', py, left[2] + [right])
            return ('bool', py, [left, right])
        ## (`**` is right-associative)
        return ('bin', op, left, self.expr(_LBP[op] - (op in ('**', '^'))))


def _emit(node):
    kind = node[0]
    if kind in ('num', 'name'):
        return node[1]
    if kind == 'group':
        return f'({_emit(node[1])})'
    if kind == 'call':
        return f"{node[1]}({', '.join(_emit(a) for a in node[2])})"
    if kind == 'unary':
        op = 'not ' if node[1] in ('!', 'not') else node[1]
        return f'({op}{_emit(node[2])})'
    if kind == 'bin':
        op = '**' if node[1] == '^' else node[1]
        return f'({_emit(node[2])} {op} {_emit(node[3])})'
    if kind == 'cmp':
        out = [_emit(node[1][0])]
        for op, item in zip(node[2], node[1][1:], strict=True):
            out += ['!=' if op == '<>' else op, _emit(item)]
        return f"({' '.join(out)})"
    if kind == 'bool':
        return f"({f' {node[1]} '.join(_emit(a) for a in node[2])})"
    return f'({_emit(node[2])} if {_emit(node[1])} else {_emit(node[3])})'


@functools.lru_cache(maxsize=65536)
def _python_text(expr):
    """A SPICE expression -- quoted, braced or bare -- as Python text."""
    e = expr.strip()
    if e[:1] in "'\"" and e[-1:] == e[:1]:
        e = e[1:-1]
    elif e.startswith('{') and e.endswith('}'):
        e = e[1:-1]
    return _emit(_Parser(e).parse())


def _strip_comments(line):
    """Drop `*` full-line comments and `;` / `$` trailing ones -- a `$`
    only where it starts a token (the line's start, or after a blank or a
    comma, as SPICE reads it: `v$d5` is a name).

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
        if ch == ';' or (ch == '$' and (out[-1] if out else ' ') in ' \t,'):
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


_MODEL_TYPE = re.compile(r'([A-Za-z_]\w*)\s*(.*)\Z', re.DOTALL)


def _model_card(text):
    """A `.model` line as (name, type, its parameters' text): the
    parameters may sit in parentheses, glued to the type
    (`NPN(BF=100 ...)`) or not, the closing one on a line of its own."""
    parts = text.split(None, 2)
    m = _MODEL_TYPE.match(parts[2]) if len(parts) == 3 else None
    if m is None:
        raise SpiceCardError('malformed .model: %r' % text[:80])
    rest = m.group(2).strip()
    if rest.startswith('('):
        rest = rest[1:-1] if rest.endswith(')') else rest[1:]
    return parts[1].lower(), m.group(1).lower(), rest


_PARAMS_KW = re.compile(r'(?<!\S)params:', re.IGNORECASE)


def _subckt_ports(text):
    """A `.subckt` line's ports: the words after its name, up to `params:`
    or its first assignment (`w=1u`, `w = 1u`)."""
    cut = len(text)
    for m in (_PARAMS_KW.search(text), _ASSIGN.search(text)):
        if m is not None:
            cut = min(cut, m.start())
    return text[:cut].split()[2:]


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
        """A SPICE expression as Python (see the module note): `1u` is
        1e-6, `1e-6` keeps its `e`, `ns1` is a name, `10pF` is 1e-11."""
        return _python_text(expr)

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
                if k.startswith(_KW):
                    key = k = k[len(_KW):]
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
            deck.subckt_ports[parts[1].lower()] = _subckt_ports(text)
            for pname, raw in _assignments(text):
                sub.define(pname, raw, None)
            scopes.append(sub)
            continue

        if head == '.ends':
            if len(scopes) > 1:
                scopes.pop()
            continue

        if head == '.model':
            name, mtype, params = _model_card(text)
            model = Model(name, mtype, scopes[-1])
            deck.models[model.name] = model
            for pname, raw in _assignments(params):
                model.raw[pname] = raw
            continue

        ## Anything else -- device lines, analyses, control blocks -- is
        ## not this reader's business.
        model = None
