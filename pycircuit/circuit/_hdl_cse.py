"""BIT-IDENTICAL COMMON-SUBEXPRESSION ELIMINATION of the generated chain
functions (2026-10-02; the fused-evaluation plan's F1a).

A chained model's `G` (and `C`) is printed by `hdl._chain_compile` as the
value chain followed by one forward-mode derivative statement per definition
and unknown, each re-printing the local partials inline -- so a partial is
evaluated once PER UNKNOWN, and the values `i` (`q`) already computed are
computed again.  PSP's `G` makes ~24.7k calls of which ~3.75k are distinct.

This pass removes the repetition from the GENERATED PYTHON, after the
compile cache, and it is exact by construction rather than by luck:

* Python evaluates every syntax-tree subtree as a unit, in the same order,
  whether it is written inline or computed into a variable first -- so a
  subtree hoisted into a temporary computes the same bits.  (Doing it in
  sympy would not be: replacing a subexpression by a symbol changes the
  canonical ordering and flattening of `Add`/`Mul` -- the association and
  the summation order -- and the last bits move.)
* the generated function is SSA (every name assigned once) and its
  primitives are pure (`numpy.*`, `minc/maxc/_step/_rdiv/_recip2`), so the
  same subtree means the same value wherever it appears after its inputs;
* nothing is evaluated lazily -- no `x if c else y`, `and`/`or` or
  comprehensions (`numpy.where` evaluates both arms) -- and the pass REFUSES
  a function that has one, leaving it as it was.

So it changes HOW MANY TIMES an operation runs, never what it computes.
Each result is verified before use: every variable and the return value,
fully expanded into its inputs, must be the same expression as the
original's (`_verify`).  Constants are keyed by type and `repr`, so `0`,
`0.0` and `-0.0` stay distinct (an integer sub-chain is what flips the sign
of a zero).

The text is produced by splicing the original source (each statement is one
line; a subtree's span is its exact text), never by `ast.unparse`, whose
recursion a PSP statement of ~300 nested operations would exhaust.  The
reference function stays reachable (`fn._hdl_ref`) and `fn._src` stays the
REFERENCE text -- what `explain()` shows, what the JAX twin re-executes and
what the C backend prints its kernels from -- with the optimised text in
`fn._src_cse`.  `fn._hdl_codelen` keeps the reference's bytecode size for
the 'auto' Newton options (`_tran_newton.compiled_jacobian_size`), whose
calibration is in that size.

Measured (F0, 2026-10-02): PSP's `G` 13.9 -> 3.2 ms, `C` 7.2 -> 1.6 ms, `i`
1.19 -> 0.90 ms; across the 22 chained library classes `G` 1.2-5.2x (median
3.0x); every function of every class byte-identical over its sweep, the
raise-mode behaviour unchanged; the PSP stage's PSS 3.5x end to end, its
waveform byte-identical.  `PYCIRCUIT_HDL_CSE=0` turns the pass off.
"""
import ast
import hashlib
import os
import sys
import tempfile

import numpy as np

#: the derived store's format; bump when the emitted text changes
FORMAT = 1

#: the pass runs at class creation unless `PYCIRCUIT_HDL_CSE=0`
ENABLED = os.environ.get('PYCIRCUIT_HDL_CSE', '1').strip() != '0'

#: the chain functions it applies to (those taking the state `x`)
FUNCS = ('i', 'q', 'G', 'C', 'i_dc', 'G_dc')

_REFUSE = (ast.IfExp, ast.BoolOp, ast.Lambda, ast.ListComp, ast.SetComp,
           ast.DictComp, ast.GeneratorExp, ast.NamedExpr, ast.Starred,
           ast.Await, ast.Yield, ast.YieldFrom, ast.JoinedStr)
_HOIST = (ast.BinOp, ast.UnaryOp, ast.Call, ast.Compare)
#: the label length of each key kind (the rest of a key is child ids)
_LABEL = {'B': 2, 'U': 2, 'C': 3, 'Cmp': 2, 'A': 2, 'S': 1, 'Tuple': 1,
          'List': 1, 'N': 2, 'K': 3}


class Refused(Exception):
    """The function has a form the pass does not handle; it is left as is."""


class _Interner:
    """Hash-consing of subtree keys: equal keys, equal ints."""

    def __init__(self):
        self.ids = {}
        self.keys = []
        self.kind = []
        self.const_only = []

    def intern(self, key, kind, const_only):
        k = self.ids.get(key)
        if k is None:
            k = self.ids[key] = len(self.keys)
            self.keys.append(key)
            self.kind.append(kind)
            self.const_only.append(const_only)
        return k


def _children(n):
    """`(children, label)`: the child expressions in a fixed order and the
    node's own part of its key."""
    if isinstance(n, ast.BinOp):
        return [n.left, n.right], ('B', type(n.op).__name__)
    if isinstance(n, ast.UnaryOp):
        return [n.operand], ('U', type(n.op).__name__)
    if isinstance(n, ast.Call):
        return ([n.func] + list(n.args) + [kw.value for kw in n.keywords],
                ('C', len(n.args), tuple(kw.arg for kw in n.keywords)))
    if isinstance(n, ast.Compare):
        return ([n.left] + list(n.comparators),
                ('Cmp', tuple(type(o).__name__ for o in n.ops)))
    if isinstance(n, ast.Attribute):
        return [n.value], ('A', n.attr)
    if isinstance(n, ast.Subscript):
        return [n.value, n.slice], ('S',)
    if isinstance(n, (ast.Tuple, ast.List)):
        return list(n.elts), (type(n).__name__,)
    if isinstance(n, ast.Name):
        return [], ('N', n.id)
    if isinstance(n, ast.Constant):
        return [], ('K', type(n.value).__name__, repr(n.value))
    raise Refused(f'unsupported node {type(n).__name__}')


def _key_tree(root, it, kid, env=None):
    """Post-order, iterative: `kid[id(node)]` for every node under `root`.
    With `env` (name -> key), a Name bound there takes that key -- the
    expansion `_verify` compares under."""
    stack = [(root, False)]
    while stack:
        n, done = stack.pop()
        if isinstance(n, _REFUSE):
            raise Refused(type(n).__name__)
        if env is not None and isinstance(n, ast.Name) and n.id in env:
            kid[id(n)] = env[n.id]
            continue
        ch, label = _children(n)
        if not done:
            stack.append((n, True))
            for c in reversed(ch):
                stack.append((c, False))
            continue
        cks = tuple(kid[id(c)] for c in ch)
        const_only = (isinstance(n, ast.Constant)
                      or (bool(ch) and not isinstance(n, ast.Name)
                          and all(it.const_only[k] for k in cks)))
        kid[id(n)] = it.intern(label + cks, type(n), const_only)
    return kid[id(root)]


def _parse(src):
    if not src.isascii():
        raise Refused('non-ASCII source (the splice works in byte offsets)')
    mod = ast.parse(src)
    if len(mod.body) != 1 or not isinstance(mod.body[0], ast.FunctionDef):
        raise Refused('not one function')
    fn = mod.body[0]
    for st in fn.body:
        if isinstance(st, ast.Assign):
            if (len(st.targets) != 1
                    or not isinstance(st.targets[0], ast.Name)):
                raise Refused('assignment form')
        elif not isinstance(st, ast.Return):
            raise Refused(f'statement {type(st).__name__}')
        if st.lineno != st.end_lineno:
            raise Refused('a statement over several lines')
    targets = [st.targets[0].id for st in fn.body
               if isinstance(st, ast.Assign)]
    if len(set(targets)) != len(targets):
        raise Refused('not SSA')
    return fn


def cse_source(src):
    """`(text, stats)`: `src` (one generated `def`) with every repeated
    non-trivial subtree computed once, verified equal (`_verify`).  Raises
    `Refused` for a function it does not handle."""
    fn = _parse(src)
    lines = src.split('\n')
    it = _Interner()
    kid = {}
    for st in fn.body:
        _key_tree(st.value, it, kid)
    ## DAG reference counts: each distinct subtree's children (with
    ## multiplicity), plus the statement roots
    refs = [0] * len(it.keys)
    for key in it.keys:
        for c in key[_LABEL[key[0]]:]:
            refs[c] += 1
    for st in fn.body:
        refs[kid[id(st.value)]] += 1

    def hoistable(k):
        return (refs[k] >= 2 and issubclass(it.kind[k], _HOIST)
                and not it.const_only[k])

    ## temporaries' names cannot collide with any name in the function
    used = {n.id for n in ast.walk(fn) if isinstance(n, ast.Name)}
    used |= {a.arg for a in fn.args.args}
    prefix = '_cse'
    while any(u.startswith(prefix) for u in used):
        prefix += '_'

    var_of = {}       # key -> the name holding it, once defined
    ntemp = [0]

    def span(n):
        return lines[n.lineno - 1][n.col_offset:n.end_col_offset]

    def splice(n, subs):
        base = span(n)
        off = n.col_offset
        out, pos = [], 0
        for c, txt in sorted(subs, key=lambda s: s[0].col_offset):
            a, b = c.col_offset - off, c.end_col_offset - off
            out.append(base[pos:a])
            out.append(txt)
            pos = b
        out.append(base[pos:])
        return ''.join(out)

    def render(root, pre):
        """The text of `root` with every repeated subtree replaced by its
        variable, appending the new temporaries' definitions to `pre` in
        evaluation order (children before parents)."""
        out = {}
        stack = [(root, False)]
        while stack:
            n, done = stack.pop()
            if not done:
                if n is not root:
                    nm = var_of.get(kid[id(n)])
                    if nm is not None:
                        out[id(n)] = nm
                        continue
                stack.append((n, True))
                for c in reversed(_children(n)[0]):
                    stack.append((c, False))
                continue
            ch = _children(n)[0]
            subs = [(c, out[id(c)]) for c in ch if out.get(id(c)) is not None]
            txt = splice(n, subs) if subs else None
            k = kid[id(n)]
            if n is not root and hoistable(k):
                t = f'{prefix}{ntemp[0]}'
                ntemp[0] += 1
                pre.append((t, span(n) if txt is None else txt))
                var_of[k] = t
                out[id(n)] = t
            else:
                out[id(n)] = txt
        txt = out[id(root)]
        return span(root) if txt is None else txt

    body = []
    for st in fn.body:
        pre = []
        k = kid[id(st.value)]
        if isinstance(st, ast.Assign):
            tgt = st.targets[0].id
            if k in var_of:
                rhs = var_of[k]
            else:
                rhs = render(st.value, pre)
                if (issubclass(it.kind[k], _HOIST)
                        and not it.const_only[k]):
                    var_of[k] = tgt
            body.extend(f'    {t} = {txt}' for t, txt in pre)
            body.append(f'    {tgt} = {rhs}')
        else:
            rhs = render(st.value, pre)
            body.extend(f'    {t} = {txt}' for t, txt in pre)
            body.append(f'    return {rhs}')
    head = lines[fn.lineno - 1:fn.body[0].lineno - 1]
    text = '\n'.join(head + body)
    _verify(fn, text)
    return text, {'statements': len(fn.body), 'temporaries': ntemp[0]}


def _canonical(fn, it):
    """Each assigned name's and the return value's expression with every
    name the function assigns expanded into its definition -- the function
    as expressions of its inputs (hash-consed, so linear)."""
    env, kid = {}, {}
    out = {}
    for st in fn.body:
        k = _key_tree(st.value, it, kid, env)
        if isinstance(st, ast.Assign):
            env[st.targets[0].id] = k
            out[st.targets[0].id] = k
        else:
            out[None] = k
    return out


def _verify(fn, text):
    """Every variable of the original and its return value are the SAME
    expressions of the inputs in the optimised text (the temporaries
    expanded), or this raises -- the guarantee the pass rests on."""
    it = _Interner()
    a, b = _canonical(fn, it), _canonical(_parse(text), it)
    for name, k in a.items():
        if b.get(name) != k:
            raise AssertionError(f'CSE changed {name!r}')


## -- the derived store ------------------------------------------------------

def _own_hash():
    with open(__file__, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


_OWN = None


def _store_path(src):
    global _OWN
    from pycircuit.circuit import _hdl_cache
    if not _hdl_cache.enabled():
        return None
    if _OWN is None:
        _OWN = _own_hash()
    h = hashlib.sha256()
    for part in (f'format={FORMAT}', _OWN, sys.version, src):
        h.update(part.encode('utf-8'))
        h.update(b'\0')
    return os.path.join(_hdl_cache.cache_dir(), 'cse', h.hexdigest() + '.py')


def _optimised_text(src):
    """The optimised text of `src`, from the store when it is there; None
    when the pass refuses the function."""
    path = _store_path(src)
    if path is not None:
        try:
            with open(path, 'r', encoding='utf-8') as fh:
                text = fh.read()
            return text or None
        except OSError:
            pass
    try:
        text, _stats = cse_source(src)
    except (Refused, AssertionError, SyntaxError, RecursionError):
        ## (an AssertionError is `_verify` refusing its own output: the
        ## reference function is kept, which is always correct)
        text = ''
    if path is not None:
        try:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path),
                                       suffix='.tmp')
            with os.fdopen(fd, 'w', encoding='utf-8') as fh:
                fh.write(text)
            os.replace(tmp, path)
        except OSError:
            pass
    return text or None


def optimise(info):
    """Replace the chained functions of `info['funcs']` (`FUNCS`) by their
    optimised twins, keeping shared identities (`i_dc is i`) and every
    attribute the rest of the package reads (`_src` the reference text,
    `_csrc`/`_cshape`/`_clayout`/`_creason` for the C backend), plus
    `_hdl_ref` (the reference function), `_src_cse` and `_hdl_codelen`.
    A function the pass refuses stays as it is.  No-op unless chained and
    `ENABLED`."""
    if not ENABLED or not info.get('chained'):
        return
    funcs = info['funcs']
    done = {}
    for name in FUNCS:
        f = funcs.get(name)
        if f is None or not hasattr(f, '_src'):
            continue
        if '_hdl_ref' in f.__dict__:
            continue                       # already optimised
        g = done.get(id(f))
        if g is None:
            text = _optimised_text(f._src)
            if text is None:
                done[id(f)] = f
                continue
            ns = f.__globals__
            loc = {}
            exec(compile(text, '<hdl-chain-cse>', 'exec'), ns, loc)  # noqa: S102 -- the generated chain, as `_chain_compile` runs it
            g = loc['_f']
            g.__dict__.update(f.__dict__)
            g._hdl_ref = f
            g._src_cse = text
            g._hdl_codelen = len(f.__code__.co_code)
            done[id(f)] = g
        funcs[name] = g


## -- fused variants (the plan's F2) ------------------------------------------
##
## Several of a chained class's functions at one state, computed by ONE
## function: their statements merged by name -- the same name is the same
## text in every function `_chain_compile` printed for the class (asserted;
## the same sympy objects through the same printer), so each output is
## computed exactly as its own function computes it -- then the union
## optimised as above.  PSP's i, q, G and C at one state: 6.3 ms as four
## optimised calls, 3.4 ms fused (F0, 2026-10-02); G and C share 1222
## derivative statements, i and q 341 values.  Built on first use for the set
## an evaluation session asks for (`take`), stored like the twins.

#: the order a fused function returns its outputs in
ORDER = ('i', 'q', 'G', 'C')


def fuse_source(srcs):
    """One generated function from several of the same signature: the union
    of their statements by target name, each function's order kept (every
    statement follows what it reads), the outputs returned as a tuple in
    the order given.  `Refused` if the signatures differ or one name stands
    for two texts."""
    head, seen, body, rets = None, {}, [], []
    for src in srcs:
        lines = src.split('\n')
        if head is None:
            head = lines[0]
        elif lines[0] != head:
            raise Refused('different signatures')
        if not lines[-1].startswith('    return '):
            raise Refused('no single-line return')
        for ln in lines[1:-1]:
            name, sep, _rhs = ln.strip().partition(' = ')
            if not sep or not ln.startswith('    '):
                raise Refused('unexpected line')
            prev = seen.get(name)
            if prev is None:
                seen[name] = ln
                body.append(ln)
            elif prev != ln:
                raise Refused(f'{name} stands for two texts')
        rets.append(lines[-1][len('    return '):])
    return '\n'.join([head] + body
                     + [f"    return ({', '.join(rets)},)"])


def fused(info, which):
    """The function computing exactly `which` (names of `ORDER`) of the
    chained class whose `info` this is, returning them in `ORDER` -- built
    on first use and kept in `info['_fused']`; None where it is not
    available (fewer than two names, a function missing, a refusal)."""
    names = tuple(k for k in ORDER if k in which)
    cache = info.setdefault('_fused', {})
    if names in cache:
        return cache[names]
    fn = None
    funcs = info.get('funcs', {})
    refs = [funcs.get(k) for k in names]
    if len(names) >= 2 and all(r is not None for r in refs):
        refs = [r.__dict__.get('_hdl_ref', r) for r in refs]
        if all(hasattr(r, '_src') for r in refs):
            try:
                src = fuse_source([r._src for r in refs])
            except Refused:
                src = None
            if src is not None:
                text = _optimised_text(src) or src
                loc = {}
                exec(compile(text, '<hdl-chain-fused>', 'exec'),  # noqa: S102 -- the generated chain
                     refs[0].__globals__, loc)
                fn = loc['_f']
                fn._hdl_names = names
    cache[names] = fn
    return fn


## -- evaluation sessions (the plan's F3) -------------------------------------

#: evaluation sessions fuse unless `PYCIRCUIT_HDL_FUSE=0` (the optimised
#: separate functions then answer every call, as outside a session)
FUSE_ENABLED = os.environ.get('PYCIRCUIT_HDL_FUSE', '1').strip() != '0'

#: a chained class's sessions engage only when its reference `G` has at least
#: this much bytecode: below it the memo's bookkeeping is the evaluation
FUSE_MIN_CODE = 5000

#: the states one element keeps per session: one site interleaves several
#: (Radau's stages, the PSS walk's points)
MEMO_ENTRIES = 8


def take(el, method, x, epar, info, args_of):
    """`el`'s `method` at `x`, from one fused pass over the session in
    force's methods (`_evalhint`), the others kept for the calls that
    follow at the same state; None for the separate call -- no session,
    `method` not in it, or a bypass: a JAX or symbolic toolkit, an instance
    attribute shadowing one of the methods (PCNR's participants), a
    function bound to the C backend, a class below `FUSE_MIN_CODE`, a state
    that is not a float vector, a fused pass that raises (the separate call
    then raises, or not, exactly as before).  Each output is handed out
    once, so no array is shared; the memo is keyed on the session, the
    state's bytes and the parameter list `args_of` returns (a new list
    after `update()` or a temperature change)."""
    from pycircuit.circuit import _evalhint
    s = _evalhint.current()
    if s is None or method not in s.which or len(s.which) < 2:
        return None
    ok = info.get('_fuse_ok')
    if ok is None:
        g = info['funcs'].get('G')
        ref = None if g is None else g.__dict__.get('_hdl_ref', g)
        ok = info['_fuse_ok'] = bool(
            ENABLED and FUSE_ENABLED and info.get('chained') and ref is not None
            and len(ref.__code__.co_code) >= FUSE_MIN_CODE)
    if not ok:
        return None
    d = el.__dict__
    if 'i' in d or 'q' in d or 'G' in d or 'C' in d:
        return None
    tk = el.toolkit
    if getattr(tk, 'jax', False) or getattr(tk, 'symbolic', False):
        return None
    funcs = info['funcs']
    for k in s.which:
        f = funcs.get(k)
        if f is None or f.__dict__.get('_hdl_c') is not None:
            return None
    if (type(x) is not np.ndarray or x.dtype != np.float64
            or x.ndim != 1):
        return None
    fm = d.get('_hdl_fm')
    if fm is None or fm[0] is not s:
        fm = d['_hdl_fm'] = (s, {})
    entries = fm[1]
    key = x.tobytes()
    args = args_of(el, epar)
    rec = entries.get(key)
    if rec is not None and rec[0] is args:
        out = rec[1].pop(method, None)
        if out is not None:
            return out
    fn = fused(info, s.which)
    if fn is None:
        return None
    try:
        vals = fn(x, *args)
    except Exception:  # noqa: BLE001 -- the separate call decides, as before
        return None
    outs = {k: np.asarray(v, dtype=float) for k, v in zip(fn._hdl_names, vals)}
    out = outs.pop(method)
    entries.pop(key, None)
    entries[key] = [args, outs]
    while len(entries) > MEMO_ENTRIES:
        entries.pop(next(iter(entries)))
    return out
