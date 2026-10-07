"""Read a SPICE netlist -- the circuit as SPICE wrote it.

A `Netlist` holds a deck's element lines (at the top level and per
subcircuit, each with its file and line), its subcircuits, its models and
parameters (scoped: a `.model` or `.param` inside a `.subckt` is that
subcircuit's), its analyses, `.ic` and `.nodeset`, options, prints and
temperatures, and the regression suite's `*COMP` directives.  Nothing is
evaluated that need not be and nothing simulated:
`pycircuit.circuit.spice_import` maps a `Netlist` onto pycircuit.  This
module imports nothing from the simulator (the utilities' layering).

The lexical rules and the expressions are `spicecard`'s: SPICE's scale
factors and units, continuation lines, `*`, `;` and token-start `$`
comments.  Names -- of elements, nodes, models, subcircuits -- are kept in
lower case, as SPICE compares them.  The first line is the title.
`.include` reads another file in place; `.lib file section` reads that
section of a model library, as `spicecard` reads a corner; a `.control`
block (ngspice's scripts) is passed over whole.  A subcircuit or a model
may be used before it is defined.

What the reader does not read it says: `unsupported` lists every dot
command it passed over, with its file and line.
"""
import os
import re
from dataclasses import dataclass, field

from pycircuit.utilities import spicecard
from pycircuit.utilities.spicecard import SpiceCardError


class SpiceNetlistError(SpiceCardError):
    """A netlist that cannot be read as written."""


@dataclass(frozen=True)
class Where:
    """A logical line's origin: its file and its first physical line."""
    file: str
    line: int

    def __str__(self):
        return f'{os.path.basename(self.file)}:{self.line}'


@dataclass
class Card:
    """One element line: its name (lower case), its positional words and
    its `name=value` parameters, raw -- parentheses and commas separate
    words, as SPICE reads an element line; a `{...}` or quoted expression
    is one word."""
    name: str
    words: list
    params: dict
    where: Where

    @property
    def kind(self):
        return self.name[0]


@dataclass
class Directive:
    """A dot command (its name without the dot, lower case): its words and
    its `name=value` parameters, raw."""
    name: str
    words: list
    params: dict
    where: Where


@dataclass
class Subckt:
    """A `.subckt`: its ports, the scope of its parameters (their defaults)
    and models, its element lines, the subcircuit it is defined in."""
    name: str
    ports: list
    scope: object
    parent: object = None
    cards: list = field(default_factory=list)
    models: dict = field(default_factory=dict)
    where: Where = None


#: An element line's or a dot command's items.
_ITEM = re.compile(r"""
    (?P<expr>\{[^}]*\}|'[^']*'|"[^"]*")
  | (?P<eq>=)
  | (?P<word>[^\s=(),{}'"]+)
  | (?P<sep>[\s(),]+)
""", re.VERBOSE)

#: A `.print` line's items: `v(a,b)` and `{...}` whole.
_PRINT_ITEM = re.compile(r"""
    (?P<expr>\{[^}]*\})
  | (?P<assign>[A-Za-z_]\w*\s*=\s*[^\s={}]+)
  | (?P<call>[A-Za-z_]\w*\s*\([^)]*\))
  | (?P<word>[^\s{}]+)
  | (?P<sep>\s+)
""", re.VERBOSE)

#: `.ic` / `.nodeset`: `v(node)=value`.
_NODE_VALUE = re.compile(r'v\s*\(\s*([^\s(),]+)\s*\)\s*=\s*([^\s=()]+)', re.IGNORECASE)

ANALYSES = ('tran', 'op', 'dc', 'ac', 'hb', 'noise', 'tf', 'sens', 'pz', 'four', 'disto')
_OPTIONS = ('options', 'option', 'opt', 'opts')


def _items(text, where):
    """`text`'s words and its `name=value` parameters (a name lower case,
    a value raw)."""
    toks, pos = [], 0
    while pos < len(text):
        m = _ITEM.match(text, pos)
        if m is None:
            raise SpiceNetlistError(f'{where}: cannot read {text[pos:pos + 20]!r}')
        pos = m.end()
        if m.lastgroup != 'sep':
            toks.append((m.lastgroup, m.group()))
    words, params, i = [], {}, 0
    while i < len(toks):
        if i + 1 < len(toks) and toks[i + 1][0] == 'eq':
            if toks[i][0] != 'word' or i + 2 >= len(toks) or toks[i + 2][0] == 'eq':
                raise SpiceNetlistError(f'{where}: a malformed assignment in {text!r}')
            params[toks[i][1].lower()] = toks[i + 2][1]
            i += 3
        elif toks[i][0] == 'eq':
            raise SpiceNetlistError(f'{where}: a malformed assignment in {text!r}')
        else:
            words.append(toks[i][1])
            i += 1
    return words, params


def _lower(words):
    """Words in lower case, as SPICE compares them -- an expression as
    written."""
    return [w if w[0] in "{'\"" else w.lower() for w in words]


def _directive(name, body, where):
    words, params = _items(body, where)
    return Directive(name, _lower(words), params, where)


def _print_items(text, where):
    words, params, pos = [], {}, 0
    while pos < len(text):
        m = _PRINT_ITEM.match(text, pos)
        if m is None:
            raise SpiceNetlistError(f'{where}: cannot read {text[pos:pos + 20]!r}')
        pos = m.end()
        if m.lastgroup == 'assign':
            name, value = m.group().split('=', 1)
            params[name.strip().lower()] = value.strip()
        elif m.lastgroup != 'sep':
            words.append(m.group())
    return words, params


class Netlist:
    """A read netlist (see the module note).  `cards` are the top level's
    element lines; `subckts` maps a name to its `Subckt`; `models` the top
    level's models (a subcircuit's are its own); `deck` the `spicecard`
    deck that resolves parameters and evaluates expressions; `analyses`,
    `options`, `prints` and `temps` are `Directive`s in file order; `ic`
    and `nodeset` map a node to its raw value (`ic_where`,
    `nodeset_where`: the first such line); `comps` are the `*COMP`
    lines' texts; `unsupported` lists (what, where) the reader passed
    over."""

    def __init__(self):
        self.title = ''
        self.deck = spicecard.Deck()
        self.cards = []
        self.subckts = {}
        self.models = {}
        self.analyses = []
        self.options = []
        self.prints = []
        self.temps = []
        self.ic = {}
        self.nodeset = {}
        self.ic_where = self.nodeset_where = None
        self.comps = []
        self.unsupported = []
        self.files = []

    def model(self, name, subckt=None):
        """The `.model` `name` as seen from inside `subckt` (a `Subckt`, or
        None for the top level): its own, then its definers', then the top
        level's.  None where there is none."""
        key, s = name.lower(), subckt
        while s is not None:
            if key in s.models:
                return s.models[key]
            s = s.parent
        return self.models.get(key)


def read(path):
    """Read the netlist `path` (see the module note)."""
    net = Netlist()
    with open(path, errors='replace') as fh:
        net.title = fh.readline().strip()
    _read_file(net, path, title=True, seen=frozenset(), stack=[])
    return net


def _read_file(net, path, title, seen, stack):
    real = os.path.realpath(path)
    if real in seen:
        raise SpiceNetlistError(f'circular .include of {path}')
    seen = seen | {real}
    net.files.append(path)
    with open(path, errors='replace') as fh:
        for n, line in enumerate(fh, 1):
            m = re.match(r'\s*\*\s*comp\b(.*)', line, re.IGNORECASE)
            if m:
                net.comps.append((m.group(1).strip(), Where(path, n)))
    skipping = None                     # a `.control` block's start
    for text, here, lineno in spicecard._numbered_lines(path, title=title):
        where = Where(path, lineno)
        head = text.split(None, 1)[0].lower()
        if skipping is not None:
            if head == '.endc':
                skipping = None
            continue
        if head == '.end':
            break
        if not head.startswith('.'):
            words, params = _items(text, where)
            card = Card(words[0].lower(), _lower(words[1:]), params, where)
            (stack[-1].cards if stack else net.cards).append(card)
            continue
        name = head[1:]
        body = text[len(head):]
        if name == 'control':
            skipping = where
        elif name == 'subckt':
            _open_subckt(net, text, where, stack)
        elif name == 'ends':
            if not stack:
                raise SpiceNetlistError(f'{where}: .ends without .subckt')
            stack.pop()
        elif name == 'model':
            mname, mtype, mparams = spicecard._model_card(text)
            scope = stack[-1].scope if stack else net.deck.global_scope
            model = spicecard.Model(mname, mtype, scope)
            for pname, raw in spicecard._assignments(mparams):
                model.raw[pname] = raw
            (stack[-1].models if stack else net.models)[mname] = model
        elif name in ('param', 'params'):
            scope = stack[-1].scope if stack else net.deck.global_scope
            for pname, raw in spicecard._assignments(body):
                scope.define(pname, raw, None)
        elif name in ('include', 'inc'):
            words = body.split()
            if not words:
                raise SpiceNetlistError(f'{where}: .include names no file')
            _read_file(net, os.path.join(here, words[0].strip('"\'')), False, seen, stack)
        elif name == 'lib' and len(body.split()) >= 2:
            _read_library(net, os.path.join(here, body.split()[0].strip('"\'')),
                          body.split()[1], stack, seen, where)
        elif name in ANALYSES:
            net.analyses.append(_directive(name, body, where))
        elif name in _OPTIONS:
            net.options.append(_directive('options', body, where))
        elif name == 'temp':
            net.temps.append(_directive(name, body, where))
        elif name == 'print':
            net.prints.append(Directive(name, *_print_items(body, where), where))
        elif name in ('ic', 'nodeset'):
            pairs = _NODE_VALUE.findall(body)
            if not pairs or _NODE_VALUE.sub('', body).strip():
                raise SpiceNetlistError(f'{where}: cannot read .{name} {body.strip()!r}')
            target = net.ic if name == 'ic' else net.nodeset
            if getattr(net, name + '_where') is None:
                setattr(net, name + '_where', where)
            for node, value in pairs:
                target[node.lower()] = value
        else:
            net.unsupported.append((text.split(None, 1)[0], where))
    if skipping is not None:
        raise SpiceNetlistError(f'{skipping}: .control without .endc')
    if title and stack:
        raise SpiceNetlistError(f'{stack[-1].where}: .subckt {stack[-1].name} without .ends')


def _open_subckt(net, text, where, stack):
    parts = text.split()
    if len(parts) < 2:
        raise SpiceNetlistError(f'{where}: malformed .subckt')
    sname = parts[1].lower()
    parent = stack[-1] if stack else None
    scope = spicecard._Scope(sname, parent.scope if parent else net.deck.global_scope)
    for pname, raw in spicecard._assignments(text):
        scope.define(pname, raw, None)
    if sname in net.subckts:
        raise SpiceNetlistError(f'{where}: .subckt {sname} defined again '
                                f'(first at {net.subckts[sname].where})')
    sub = Subckt(sname, [p.lower() for p in spicecard._subckt_ports(text)], scope,
                 parent=parent, where=where)
    net.subckts[sname] = sub
    stack.append(sub)


def _read_library(net, path, section, stack, seen, where):
    """`.lib file section`: the section's models and parameters, read as
    `spicecard` reads a corner, into the current scope."""
    before = set(net.deck.models)
    scope = stack[-1].scope if stack else net.deck.global_scope
    try:
        spicecard._parse_file(net.deck, path, section, scope, set(seen))
    except SpiceCardError as e:
        raise SpiceNetlistError(f'{where}: {e}') from e
    target = stack[-1].models if stack else net.models
    for mname in set(net.deck.models) - before:
        target[mname] = net.deck.models[mname]
