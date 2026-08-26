"""The backend spike of doc/hdl_roadmap_260824.md section 14.5.

A MEASUREMENT, not a build.  Nothing here ships as a backend; the point
is to let the numbers choose between numba, C-via-cffi and llvmlite for
the chained code generator's output, on the real flagship model, before
anyone commits weeks to one of them.

Input: `PspMosLongChannel._hdl_info['funcs']['G']._src` (and `i`) --
the straight-line scalar source `_chain_compile` emits.  It is taken
from the LIVE class, so this harness measures whatever the printer
emits today, not a copy that can drift.

Candidates, cheapest first:

  numba   the emitted source, `@njit`, with ONE mechanical rewrite: the
          heterogeneous nested-list `return` becomes stores into a
          preallocated array (numba cannot type `[[float, 0], ...]`).
          `numpy.where`, the ufunc comparisons, `abs`, `**` and `sign`
          all compile as emitted; the kernel helpers (`maxc`, `minc`,
          `_step`, `_rdiv`, `_recip2`) are re-declared as `@njit`
          twins of their numpy runtime forms.
  C       the same AST transpiled to C: every `numpy.where(c, a, b)`
          becomes `_sel(c, a, b)` -- a function CALL, so both arms are
          evaluated as arguments exactly as Python evaluates them
          before `where` picks; `numpy.maximum`/`minimum` become
          NaN-propagating twins of numpy's own definition; `**2` is
          `x*x` and any other power is libm `pow`.  Compiled with the
          system `cc` at -O2, NO fast-math, called through cffi.

Both are given the parameter vector two ways: the 155-argument
signature as emitted, and packed into one array (`p[k]`), because the
argument-marshalling cost is part of what a backend would pay.

Measured for each (section 14.5): per-call `G` and `i` at the section
13 bias, compile time, agreement with the numpy path over a bias sweep
that INCLUDES the extremes (+-1e30, +-100) the safety primitives exist
for, and finiteness wherever the numpy path is finite.

Findings (2026-08-26, doc/backend_spike_260826.md has the table):

  - C, bitwise-faithful (`-fno-builtin-pow`): G 75 us (229x), i 7.7 us,
    0 ulp against numpy at all 4985 sweep points, nothing lost finite.
    With `x*x` for `**2` (`--pow2-mul`): G 14 us (1235x) and a worst
    2.4e-10 relative -- ONE ulp in one square, amplified by a
    cancellation in the derivative chain.  numpy's scalar `**2` is
    glibc `pow`, which is not correctly rounded; `x*x` is.
  - numba: fine on `i` (9 s JIT, 6.6 us); on `G` the default JIT did
    not return in 600 s, and `NUMBA_OPT=0` gave 68 s / 622 us (28x).
    Don't.  llvmlite skipped for the same reason.
  - the numpy path itself is non-finite at 96 of the 4096 extreme
    grid points (Vb = +-1e30 with Vs or Vd at +-1e30); no candidate
    changed that in either direction.

Run:  python benchmarks/backend_spike.py                 (everything)
      python benchmarks/backend_spike.py --no-numba      (C only)
      python benchmarks/backend_spike.py --no-numba --pow2-mul
      NUMBA_OPT=0 python benchmarks/backend_spike.py --no-c --only G \
          --numba-pack-only
      python benchmarks/backend_spike.py --trace G -0.05 -0.525 0 0
                                (first differing intermediate, C vs numpy)
      python benchmarks/backend_spike.py --json out.json --sweep-n 400

Per-call timings are only comparable on an idle machine; the numba
default run on G is expected not to finish.
"""

import argparse
import ast
import itertools
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hdl_model_cost import PDK, T27          # noqa: E402

EXTREMES = (-1e30, -100.0, -1.0, 0.0, 0.7, 1.2, 100.0, 1e30)


def _bench(fn, n, warmup=3):
    for _ in range(warmup):
        fn()
    t0 = time.perf_counter()
    for _ in range(n):
        fn()
    return (time.perf_counter() - t0) / n


def _rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


## ----------------------------------------------------------------------
## The model, and its emitted source.

def build_psp():
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import psp_scaling
    from pycircuit.utilities import spicecard
    cm.default_toolkit = numeric
    deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                          section='mos_tt')
    w, l = 10e-6, 1e-6
    kw = psp_scaling.to_long_channel(
        deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1, m=1,
                          pre_layout=1), w=w, l=l, T=T27)
    t0 = time.perf_counter()
    from pycircuit.circuit.compact import PspMosLongChannel
    e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                          cm.Node('b'), **kw)
    e.update_iparv()
    return e, time.perf_counter() - t0


def emitted_source(e, which):
    f = type(e)._hdl_info['funcs'][which]
    src = getattr(f, '_src', None) or getattr(f, '_hdl_inner')._src
    return f, src


def arg_vector(e):
    from pycircuit.circuit import hdl
    from pycircuit.circuit.circuit import defaultepar
    return [float(v) for v in hdl._args_of(e, defaultepar)]


## ----------------------------------------------------------------------
## The mechanical rewrites, on the AST.

class _FoldReduce(ast.NodeTransformer):
    """`numpy.logical_or.reduce((a, b, c))` -> nested binary calls.

    Sympy's `NumPyPrinter` prints `Or`/`And` this way (section 8 noted
    the same escape for `Min`/`Max`).  The tuple is built -- every
    element evaluated -- before `reduce` runs, and the nested form
    evaluates exactly the same elements in the same order, so this
    changes nothing but the call shape.  numba cannot type
    `ufunc.reduce` on a tuple; C has no such thing.
    """

    def visit_Call(self, node):
        self.generic_visit(node)
        f = node.func
        if isinstance(f, ast.Attribute) and f.attr == 'reduce' and \
                isinstance(f.value, ast.Attribute) and \
                isinstance(f.value.value, ast.Name) and \
                f.value.value.id == 'numpy' and len(node.args) == 1 and \
                isinstance(node.args[0], (ast.Tuple, ast.List)):
            elts = node.args[0].elts
            out = elts[0]
            for nxt in elts[1:]:
                out = ast.Call(func=f.value, args=[out, nxt], keywords=[])
            return ast.copy_location(out, node)
        return node


def parse_fn(src):
    tree = ast.parse(src)
    fdef = _FoldReduce().visit(tree.body[0])
    ast.fix_missing_locations(fdef)
    assert isinstance(fdef, ast.FunctionDef) and \
        isinstance(fdef.body[-1], ast.Return)
    return fdef


def return_shape(fdef):
    """The `return` as a list of (index-tuple, expr-node)."""
    val = fdef.body[-1].value
    assert isinstance(val, ast.List)
    if val.elts and all(isinstance(r, ast.List) for r in val.elts):
        cells = [((r, c), x) for r, row in enumerate(val.elts)
                 for c, x in enumerate(row.elts)]
        shape = (len(val.elts), len(val.elts[0].elts))
    else:
        cells = [((c,), x) for c, x in enumerate(val.elts)]
        shape = (len(val.elts),)
    return shape, cells


class _WhereCall(ast.NodeTransformer):
    """`numpy.where(c, a, b)` -> `_where(c, a, b)`, for numba only.

    numba types `numpy.where` on scalars as numpy does -- a 0-d ARRAY
    -- and then has no `abs()` for that type, so the emitted source
    fails to type as it stands.  `_where` is an njit scalar select.
    Both arms are still evaluated: they are ARGUMENTS, computed before
    the call, exactly as they are for `numpy.where`.  The value is the
    same by definition (`where` on scalars returns one of its two
    arguments), and numpy's truth test (nonzero or NaN) is what `if c`
    does in numba.
    """

    def visit_Call(self, node):
        self.generic_visit(node)
        f = node.func
        if isinstance(f, ast.Attribute) and f.attr == 'where' and \
                isinstance(f.value, ast.Name) and f.value.id == 'numpy':
            node.func = ast.copy_location(ast.Name(id='_where',
                                                   ctx=ast.Load()), f)
        return node


class _PackParams(ast.NodeTransformer):
    """`name` -> `p[k]` for every parameter name.  Nothing else moves."""

    def __init__(self, names):
        self.idx = {n: k for k, n in enumerate(names)}

    def visit_Name(self, node):
        k = self.idx.get(node.id)
        if k is None or not isinstance(node.ctx, ast.Load):
            return node
        return ast.copy_location(
            ast.Subscript(value=ast.Name(id='p', ctx=ast.Load()),
                          slice=ast.Constant(value=k), ctx=ast.Load()),
            node)


def python_array_return(src, pack=False, where_call=False):
    """The emitted source with its `return [[...]]` turned into stores.

    This is the ONE change numba needs, and it changes no arithmetic:
    every cell expression is unparsed verbatim (`ast.unparse` prints
    floats by `repr`, which round-trips) into `_out[r, c] = <expr>`.
    """
    fdef = parse_fn(src)
    argnames = [a.arg for a in fdef.args.args]
    params = argnames[1:]
    shape, cells = return_shape(fdef)
    body = fdef.body[:-1]
    if where_call:
        tw = _WhereCall()
        body = [tw.visit(s) for s in body]
        cells = [(ix, tw.visit(x)) for ix, x in cells]
    if pack:
        tr = _PackParams(params)
        body = [tr.visit(s) for s in body]
        cells = [(ix, tr.visit(x)) for ix, x in cells]
        sig = 'x, p'
    else:
        sig = ', '.join(argnames)
    lines = ['def _f(%s):' % sig]
    lines += ['    ' + ast.unparse(s) for s in body]
    lines.append('    _out = numpy.empty(%r)' % (shape,))
    for ix, x in cells:
        lines.append('    _out[%s] = %s' % (', '.join(map(str, ix)),
                                            ast.unparse(x)))
    lines.append('    return _out')
    return '\n'.join(lines), shape, params


## ----------------------------------------------------------------------
## Candidate 1: numba.

def numba_kernel_namespace():
    """The chain namespace, with the kernel helpers as njit twins.

    Each is the numpy runtime form from hdl.py with the sympy/jax
    dispatch guards removed -- numba cannot compile a try/except around
    a jax import, and inside an njit function there is nothing to
    dispatch on.  The ARITHMETIC is identical, which is the property
    the agreement sweep then checks.
    """
    import numba
    ns = dict(numpy=np)
    ns['maxc'] = numba.njit(lambda a, b: np.maximum(a, b))
    ns['minc'] = numba.njit(lambda a, b: np.minimum(a, b))
    ns['_step'] = numba.njit(lambda a, b: 1.0 * (a >= b))
    ns['_rdiv'] = numba.njit(lambda b, e: b / (b * b + e * e))
    ns['_recip2'] = numba.njit(lambda b, e: 1.0 / (b * b + e * e))
    ns['_wrapfloor'] = numba.njit(lambda a: np.floor(a))

    @numba.njit
    def _where(c, a, b):
        return a if c else b
    ns['_where'] = _where
    return ns


def numba_compile(src, x, args, pack):
    import numba
    psrc, shape, params = python_array_return(src, pack=pack,
                                              where_call=True)
    ns = numba_kernel_namespace()
    exec(compile(psrc, '<spike-numba>', 'exec'), ns)
    jf = numba.njit(ns['_f'])
    call = (lambda x_: jf(x_, np.asarray(args))) if pack else \
        (lambda x_: jf(x_, *args))
    rss0 = _rss_mb()
    t0 = time.perf_counter()
    call(x)                                   # JIT happens here
    t_compile = time.perf_counter() - t0
    return call, dict(compile_s=t_compile, rss_delta_mb=_rss_mb() - rss0,
                      src_lines=psrc.count('\n') + 1)


## ----------------------------------------------------------------------
## Candidate 2: C via cffi.

C_PRELUDE = r"""
#include <math.h>
/* Every helper is a FUNCTION: C evaluates all arguments before the
   call, so `_sel(c, a, b)` computes both arms exactly as Python does
   for `numpy.where(c, a, b)`, and only then picks.  `c != 0.0` is
   numpy's truth test: nonzero, and NaN, are true. */
static inline double _sel(double c, double a, double b) { return (c != 0.0) ? a : b; }
/* numpy's maximum/minimum, verbatim: NaN in either argument wins. */
static inline double _npmax(double a, double b) { return (a >= b || a != a) ? a : b; }
static inline double _npmin(double a, double b) { return (a <= b || a != a) ? a : b; }
static inline double _step(double a, double b) { return 1.0 * (a >= b); }
static inline double _rdiv(double b, double e) { return b / (b * b + e * e); }
static inline double _recip2(double b, double e) { return 1.0 / (b * b + e * e); }
static inline double _sq(double a) { return a * a; }
static inline double _sign(double a) { return (a != a) ? a : (a > 0.0) ? 1.0 : (a < 0.0) ? -1.0 : 0.0; }
static inline double _land(double a, double b) { return (a != 0.0) && (b != 0.0); }
static inline double _lor(double a, double b) { return (a != 0.0) || (b != 0.0); }
"""

_C_CALLS = {
    'numpy.where': '_sel', 'numpy.maximum': '_npmax',
    'numpy.minimum': '_npmin', 'maxc': '_npmax', 'minc': '_npmin',
    '_step': '_step', '_rdiv': '_rdiv', '_recip2': '_recip2',
    'numpy.sqrt': 'sqrt', 'numpy.exp': 'exp', 'numpy.log': 'log',
    'numpy.sign': '_sign', 'abs': 'fabs', 'numpy.abs': 'fabs',
    'numpy.floor': 'floor', '_wrapfloor': 'floor',
    'numpy.logical_and': '_land', 'numpy.logical_or': '_lor',
    'numpy.tanh': 'tanh', 'numpy.sinh': 'sinh', 'numpy.cosh': 'cosh',
}
_C_CMP = {'numpy.greater': '>', 'numpy.greater_equal': '>=',
          'numpy.less': '<', 'numpy.less_equal': '<=',
          'numpy.equal': '==', 'numpy.not_equal': '!='}
_C_CONST = {'numpy.nan': 'NAN', 'numpy.inf': 'INFINITY',
            'numpy.pi': 'M_PI'}
_C_BINOP = {ast.Add: '+', ast.Sub: '-', ast.Mult: '*', ast.Div: '/'}

## numpy's scalar `x ** 2` is libm `pow(x, 2.0)` -- measured: it differs
## from `x * x` in 13 of 20000 random doubles, because glibc's pow is
## within 0.52 ulp, not correctly rounded, where `x * x` is.  So the
## bitwise-faithful emission is `pow`; `x * x` is the (very slightly
## more accurate) alternative, kept as an option to price it.
POW2_AS_MUL = False
_C_CMPOP = {ast.Gt: '>', ast.GtE: '>=', ast.Lt: '<', ast.LtE: '<=',
            ast.Eq: '==', ast.NotEq: '!='}


def _c_const(v):
    if isinstance(v, bool):
        return '1.0' if v else '0.0'
    if isinstance(v, int):
        return '%d.0' % v
    r = repr(float(v))
    if r == 'inf':
        return 'INFINITY'
    if r == '-inf':
        return '(-INFINITY)'
    if r == 'nan':
        return 'NAN'
    if '.' not in r and 'e' not in r:
        r += '.0'
    return r


def _callee(node):
    f = node.func
    if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name):
        return '%s.%s' % (f.value.id, f.attr)
    if isinstance(f, ast.Name):
        return f.id
    raise NotImplementedError(ast.dump(f))


def cexpr(n):
    """Python expression AST -> C, in double arithmetic."""
    if isinstance(n, ast.Constant):
        return _c_const(n.value)
    if isinstance(n, ast.Name):
        return n.id
    if isinstance(n, ast.Attribute):
        key = '%s.%s' % (n.value.id, n.attr)
        return _C_CONST[key]
    if isinstance(n, ast.Subscript):
        assert isinstance(n.slice, ast.Constant)
        return '%s[%d]' % (cexpr(n.value), n.slice.value)
    if isinstance(n, ast.UnaryOp):
        if isinstance(n.op, ast.USub):
            return '(-%s)' % cexpr(n.operand)
        if isinstance(n.op, ast.UAdd):
            return cexpr(n.operand)
        raise NotImplementedError(ast.dump(n))
    if isinstance(n, ast.BinOp):
        if isinstance(n.op, ast.Pow):
            r = n.right
            if POW2_AS_MUL and isinstance(r, ast.Constant) and \
                    not isinstance(r.value, bool) and float(r.value) == 2.0:
                return '_sq(%s)' % cexpr(n.left)
            return 'pow(%s, %s)' % (cexpr(n.left), cexpr(r))
        return '(%s %s %s)' % (cexpr(n.left), _C_BINOP[type(n.op)],
                               cexpr(n.right))
    if isinstance(n, ast.Compare):
        assert len(n.ops) == 1
        return '(%s %s %s)' % (cexpr(n.left), _C_CMPOP[type(n.ops[0])],
                               cexpr(n.comparators[0]))
    if isinstance(n, ast.Call):
        name = _callee(n)
        a = [cexpr(x) for x in n.args]
        if name in _C_CMP:
            assert len(a) == 2
            return '(%s %s %s)' % (a[0], _C_CMP[name], a[1])
        if name in _C_CALLS:
            return '%s(%s)' % (_C_CALLS[name], ', '.join(a))
        raise NotImplementedError('call %s' % name)
    raise NotImplementedError(ast.dump(n))


def c_source(src, fname):
    """The emitted Python, as one C function `fname(x, p, out)`."""
    fdef = parse_fn(src)
    argnames = [a.arg for a in fdef.args.args]
    params = argnames[1:]
    shape, cells = return_shape(fdef)
    tr = _PackParams(params)
    out = [C_PRELUDE,
           'void %s(const double *x, const double *p, double *out) {'
           % fname]
    for s in fdef.body[:-1]:
        assert isinstance(s, ast.Assign) and len(s.targets) == 1 and \
            isinstance(s.targets[0], ast.Name)
        out.append('  const double %s = %s;'
                   % (s.targets[0].id, cexpr(tr.visit(s.value))))
    ncol = shape[1] if len(shape) == 2 else 1
    for ix, x in cells:
        flat = ix[0] * ncol + (ix[1] if len(ix) == 2 else 0)
        out.append('  out[%d] = %s;' % (flat, cexpr(tr.visit(x))))
    out.append('}')
    return '\n'.join(out), shape, params


def c_compile(src, fname, workdir, cc='cc', opt='-O2'):
    csrc, shape, params = c_source(src, fname)
    cpath = os.path.join(workdir, fname + '.c')
    sopath = os.path.join(workdir, 'lib%s.so' % fname)
    with open(cpath, 'w') as fh:
        fh.write(csrc)
    ## No -ffast-math, no -march=native: strict IEEE double, no FMA
    ## contraction, so the float path is the one the numpy code walks.
    ## `-fno-builtin-pow` because gcc folds `pow(x, 2.0)` to `x*x` at
    ## -O1 and above even without fast-math -- and numpy's scalar
    ## `x**2` is glibc's pow, which is NOT correctly rounded, so the
    ## folded form is one ulp off numpy in ~0.07% of cases.  Measured
    ## (dissect at bias (-0.05, -0.525, 0, 0)): that single ulp became
    ## a 2.4e-10 relative difference in G[0, 2] through a cancellation.
    cmd = [cc, opt, '-shared', '-fPIC', '-fno-fast-math',
           '-ffp-contract=off', '-fno-builtin-pow', '-o', sopath, cpath,
           '-lm']
    t0 = time.perf_counter()
    subprocess.run(cmd, check=True)
    t_compile = time.perf_counter() - t0
    import cffi
    ffi = cffi.FFI()
    ffi.cdef('void %s(const double *x, const double *p, double *out);'
             % fname)
    lib = ffi.dlopen(sopath)
    cfn = getattr(lib, fname)
    return cfn, ffi, dict(compile_s=t_compile, c_lines=csrc.count('\n') + 1,
                          shape=shape, cmd=' '.join(cmd),
                          so_bytes=os.path.getsize(sopath))


def c_caller(cfn, ffi, args, shape):
    p = np.asarray(args, float)
    pp = ffi.cast('double *', p.ctypes.data)
    out = np.empty(shape)
    op = ffi.cast('double *', out.ctypes.data)
    ## The bias changes every call, so its pointer is taken per call --
    ## that cost belongs to the candidate.  `out` is REUSED: the caller
    ## must copy if it keeps it, as it must with any C stamp.
    keep = (p, out)

    def call(x_):
        cfn(ffi.cast('double *', x_.ctypes.data), pp, op)
        return out
    call._keep = keep
    return call


## ----------------------------------------------------------------------
## Agreement.

def compare(ref, cand):
    """Per-entry agreement of `cand` against the numpy result `ref`.

    Returns (worst_rel, worst_abs, worst_ulp, n_lost_finite, n_gained,
    n_both_nonfinite, argmax entry) over entries where both are finite;
    an entry finite in numpy and not in the candidate is a LOSS and the
    section 14.6 failure.
    """
    ref = np.asarray(ref, float).ravel()
    cand = np.asarray(cand, float).ravel()
    fr, fc = np.isfinite(ref), np.isfinite(cand)
    both = fr & fc
    lost = int(np.sum(fr & ~fc))
    gained = int(np.sum(~fr & fc))
    nonf = int(np.sum(~fr & ~fc))
    if not both.any():
        return dict(rel=0.0, abs=0.0, ulp=0, lost=lost, gained=gained,
                    both_nonfinite=nonf, at=None)
    a, b = ref[both], cand[both]
    d = np.abs(a - b)
    den = np.maximum(np.abs(a), np.abs(b))
    rel = np.where(den > 0, d / np.where(den > 0, den, 1.0), 0.0)
    ia = a.view(np.int64)
    ib = b.view(np.int64)
    same_sign = (np.signbit(a) == np.signbit(b))
    ulp = np.where(same_sign, np.abs(ia - ib), np.iinfo(np.int64).max)
    ulp = np.where(d == 0, 0, ulp)
    k = int(np.argmax(rel))
    return dict(rel=float(rel[k]), abs=float(d.max()), ulp=int(ulp.max()),
                lost=lost, gained=gained, both_nonfinite=nonf,
                at=int(np.flatnonzero(both)[k]),
                ref_at=float(a[k]), cand_at=float(b[k]))


def sweep_points(e, n_extremes):
    pts = []
    combos = list(itertools.product(EXTREMES, repeat=4))
    if n_extremes and n_extremes < len(combos):
        rng = np.random.default_rng(0)
        idx = rng.choice(len(combos), n_extremes, replace=False)
        combos = [combos[i] for i in sorted(idx)]
    for vd, vg, vs, vb in combos:
        pts.append(('extreme', (vd, vg, vs, vb)))
    try:
        from psp_reference import SWEEPS
        for s in SWEEPS:
            vals = np.arange(s['start'], s['stop'] + s['step'] / 2,
                             s['step'])
            for v in vals:
                b = dict(s['bias'])
                b[s['sweep']] = float(v)
                pts.append((s['name'], (b.get('Vd', 0.0), b.get('Vg', 0.0),
                                        b.get('Vs', 0.0), b.get('Vb', 0.0))))
    except Exception as exc:                          # pragma: no cover
        print('  (reference sweeps not enumerated: %s)' % exc)
    return pts


def run_sweep(e, pts, ref_fn, cands):
    """`cands`: {name: callable(x) -> array}.  Returns per-candidate worst."""
    worst = {n: dict(rel=0.0, abs=0.0, ulp=0, lost=0, gained=0,
                     both_nonfinite=0, at=None, where=None,
                     lost_where=[]) for n in cands}
    ref_raised = 0
    ref_nonfinite_pts = 0
    ref_finite_entries = 0
    t0 = time.perf_counter()
    for k, (label, (vd, vg, vs, vb)) in enumerate(pts):
        x = np.asarray(e.bias(vd, vg, vs, vb), float)
        with np.errstate(all='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                ref = np.asarray(ref_fn(x), float)
            except Exception as exc:
                ref_raised += 1
                ref = None
            outs = {n: np.array(f(x), float, copy=True)
                    for n, f in cands.items()}
        if ref is None:
            continue
        if not np.isfinite(ref).all():
            ref_nonfinite_pts += 1
        ref_finite_entries += int(np.isfinite(ref).sum())
        for n, o in outs.items():
            c = compare(ref, o)
            w = worst[n]
            if c['lost']:
                w['lost'] += c['lost']
                if len(w['lost_where']) < 5:
                    w['lost_where'].append(
                        dict(label=label, bias=(vd, vg, vs, vb),
                             entries=[int(i) for i in np.flatnonzero(
                                 np.isfinite(ref.ravel())
                                 & ~np.isfinite(o.ravel()))][:6]))
            w['gained'] += c['gained']
            w['both_nonfinite'] += c['both_nonfinite']
            w['abs'] = max(w['abs'], c['abs'])
            w['ulp'] = max(w['ulp'], c['ulp'])
            if c['rel'] > w['rel']:
                w['rel'] = c['rel']
                w['at'] = c['at']
                w['where'] = dict(label=label, bias=(vd, vg, vs, vb),
                                  entry=c['at'], ref=c['ref_at'],
                                  cand=c['cand_at'])
        if (k + 1) % 500 == 0:
            print('    %d/%d points, %.0f s' % (k + 1, len(pts),
                                                time.perf_counter() - t0),
                  flush=True)
    return dict(points=len(pts), ref_raised=ref_raised,
                ref_nonfinite_points=ref_nonfinite_pts,
                ref_finite_entries=ref_finite_entries,
                elapsed_s=time.perf_counter() - t0, worst=worst)


## ----------------------------------------------------------------------
## Trace: where does C first leave numpy?

def trace(e, which, bias, workdir, cc='cc', opt='-O2'):
    """Every intermediate of `which` at one bias, numpy against C.

    The worst OUTPUT entry is far from the cause: at the bias this was
    written for, G[0, 2] differs by 2.4e-10 and the cause is one ulp
    in one `**2` two thousand lines earlier.  So compare the chain,
    not the result, and report the first line that differs.
    """
    import cffi
    from pycircuit.circuit import hdl
    f, src = emitted_source(e, which)
    fdef = parse_fn(src)
    names = [s.targets[0].id for s in fdef.body[:-1]]
    psrc, shape, params = python_array_return(src, pack=True)
    lines = [l for l in psrc.split('\n')
             if not l.startswith('    _out') and not l.startswith('    return')]
    lines.append('    return [%s]' % ', '.join(names))
    ns = hdl._chain_namespace(dict(hdl._KERNEL_NUMPY, _wrapfloor=np.floor))
    exec('\n'.join(lines), ns)
    args = arg_vector(e)
    x = np.asarray(e.bias(*bias), float)
    with np.errstate(all='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore')
        py = np.array([float(v) for v in ns['_f'](x, np.asarray(args))])
    csrc, _, _ = c_source(src, 'tr')
    clines = [l for l in csrc.split('\n') if not l.startswith('  out[')]
    clines = clines[:-1] + ['  out[%d] = %s;' % (k, n)
                            for k, n in enumerate(names)] + ['}']
    cpath = os.path.join(workdir, 'tr_%s.c' % which)
    sopath = os.path.join(workdir, 'libtr_%s.so' % which)
    with open(cpath, 'w') as fh:
        fh.write('\n'.join(clines))
    subprocess.run([cc, opt, '-shared', '-fPIC', '-fno-fast-math',
                    '-ffp-contract=off', '-fno-builtin-pow', '-o', sopath,
                    cpath, '-lm'], check=True)
    ffi = cffi.FFI()
    ffi.cdef('void tr(const double*, const double*, double*);')
    lib = ffi.dlopen(sopath)
    out = np.empty(len(names))
    p = np.asarray(args, float)
    lib.tr(ffi.cast('double*', x.ctypes.data), ffi.cast('double*', p.ctypes.data),
           ffi.cast('double*', out.ctypes.data))
    d = np.abs(py - out)
    den = np.maximum(np.abs(py), np.abs(out))
    rel = np.where(den > 0, d / np.where(den > 0, den, 1), 0)
    rel = np.where((np.isfinite(py) != np.isfinite(out)), np.inf, rel)
    print('%s at bias %s: %d intermediates' % (which, bias, len(names)))
    first = np.flatnonzero(rel > 0)
    if not len(first):
        print('  bitwise identical, every intermediate')
        return
    src_of = {s.targets[0].id: ast.unparse(s.value) for s in fdef.body[:-1]}
    print('  first differing:')
    for k in first[:8]:
        print('    %-28s py=%.17g c=%.17g rel=%.2e'
              % (names[k], py[k], out[k], rel[k]))
    print('  worst:')
    for k in np.argsort(-rel)[:8]:
        print('    %-28s py=%.17g c=%.17g rel=%.2e'
              % (names[k], py[k], out[k], rel[k]))
    k = first[0]
    print('  %s = %s' % (names[k], src_of[names[k]][:500]))


def verdict(speedup, w):
    if w is None:
        return 'not measured'
    if w['lost']:
        return ('FAILS 14.6: %d entries finite in numpy and not here'
                % w['lost'])
    if speedup >= 50:
        return 'worth building out (>= 50x, full-sweep agreement)'
    if speedup < 10:
        return 'not worth the surface area (< 10x)'
    return 'between 10x and 50x: report and decide'


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--no-numba', action='store_true')
    ap.add_argument('--no-c', action='store_true')
    ap.add_argument('--no-sweep', action='store_true')
    ap.add_argument('--sweep-n', type=int, default=0,
                    help='extreme-grid points (0 = all 4096)')
    ap.add_argument('--reps', type=int, default=200)
    ap.add_argument('--cc', default='cc')
    ap.add_argument('--opt', default='-O2')
    ap.add_argument('--json', metavar='PATH')
    ap.add_argument('--workdir', default=None)
    ap.add_argument('--only', choices=('G', 'i'), default=None,
                    help='one function only')
    ap.add_argument('--numba-pack-only', action='store_true',
                    help='skip the 155-argument numba variant')
    ap.add_argument('--pow2-mul', action='store_true',
                    help='emit x*x for **2 in C instead of pow(x, 2.0)')
    ap.add_argument('--trace', nargs=5, metavar=('FN', 'VD', 'VG', 'VS', 'VB'),
                    help='intermediates of FN at one bias, numpy vs C')
    args = ap.parse_args(argv)
    global POW2_AS_MUL
    POW2_AS_MUL = args.pow2_mul
    WHICH = (args.only,) if args.only else ('G', 'i')

    if not os.path.isdir(PDK):
        print('IHP Open PDK not at %s' % PDK)
        return 1
    workdir = args.workdir or tempfile.mkdtemp(prefix='backend-spike-')
    print('workdir', workdir)

    e, t_inst = build_psp()
    print('PspMosLongChannel instantiated in %.2f s (%s)'
          % (t_inst, type(e)._hdl_cache_status))
    if args.trace:
        trace(e, args.trace[0], tuple(float(v) for v in args.trace[1:]),
              workdir, cc=args.cc, opt=args.opt)
        return 0
    x0 = np.asarray(e.bias(1.2, 0.8, 0.0, 0.0), float)
    argv_ = arg_vector(e)
    print('argument vector: %d values' % len(argv_))

    result = dict(instantiate_s=t_inst, python=sys.version.split()[0],
                  numpy=np.__version__, rows={})
    fns = {}
    for which in WHICH:
        f, src = emitted_source(e, which)
        fns[which] = dict(fn=f, src=src)
        with open(os.path.join(workdir, '%s_src.py' % which), 'w') as fh:
            fh.write(src)
        print('%s: %d emitted lines' % (which, src.count('\n') + 1))

    ## --- baseline -----------------------------------------------------
    rows = {}
    for which in WHICH:
        f = fns[which]['fn']
        raw = lambda x_, _f=f: _f(x_, *argv_)
        elem = getattr(e, which)
        t_raw = _bench(lambda: raw(x0), args.reps)
        t_elem = _bench(lambda: elem(x0), args.reps)
        rows[which] = dict(numpy_raw_us=t_raw * 1e6, numpy_elem_us=t_elem * 1e6)
        fns[which]['raw'] = raw
        print('numpy %s: %.1f us raw, %.1f us through the element'
              % (which, t_raw * 1e6, t_elem * 1e6))
    ## The section-13 numbers are the element path; that is the baseline.

    cands = {'G': {}, 'i': {}}
    meta = {}

    ## --- numba ----------------------------------------------------------
    if not args.no_numba:
        import numba
        result['numba'] = numba.__version__
        for which in WHICH:
            for pack in ((True,) if args.numba_pack_only
                         else (True, False)):
                name = 'numba' + ('' if pack else '-155args')
                print('compiling %s %s ...' % (name, which), flush=True)
                try:
                    call, m = numba_compile(fns[which]['src'], x0, argv_,
                                            pack)
                except Exception as exc:
                    print('  FAILED: %s: %s' % (type(exc).__name__,
                                                str(exc)[:400]))
                    meta[(name, which)] = dict(error=str(exc)[:400])
                    continue
                t = _bench(lambda: call(x0), args.reps)
                m['us'] = t * 1e6
                meta[(name, which)] = m
                cands[which][name] = call
                print('  compile %.1f s, %.2f us per call, +%.0f MB RSS'
                      % (m['compile_s'], m['us'], m['rss_delta_mb']))

    ## --- C ---------------------------------------------------------------
    if not args.no_c:
        for which in WHICH:
            name = 'C-cffi'
            print('compiling %s %s (%s %s) ...' % (name, which, args.cc,
                                                   args.opt), flush=True)
            try:
                cfn, ffi, m = c_compile(fns[which]['src'], 'psp_' + which,
                                        workdir, cc=args.cc, opt=args.opt)
            except Exception as exc:
                print('  FAILED: %s: %s' % (type(exc).__name__,
                                            str(exc)[:400]))
                meta[(name, which)] = dict(error=str(exc)[:400])
                continue
            call = c_caller(cfn, ffi, argv_, m['shape'])
            t = _bench(lambda: call(x0), args.reps)
            m['us'] = t * 1e6
            meta[(name, which)] = m
            cands[which][name] = call
            print('  compile %.1f s, %.2f us per call, %d C lines, %d B .so'
                  % (m['compile_s'], m['us'], m['c_lines'], m['so_bytes']))

    ## --- sweep -------------------------------------------------------------
    sweeps = {}
    if not args.no_sweep:
        pts = sweep_points(e, args.sweep_n)
        print('agreement sweep: %d bias points (%d extreme-grid)'
              % (len(pts), sum(1 for p in pts if p[0] == 'extreme')),
              flush=True)
        for which in WHICH:
            if not cands[which]:
                continue
            print('  %s ...' % which, flush=True)
            sweeps[which] = run_sweep(e, pts, fns[which]['raw'],
                                      cands[which])

    ## --- report ---------------------------------------------------------------
    print()
    print('%-16s %-3s %10s %8s %10s %12s %8s %8s  %s'
          % ('candidate', 'fn', 'us/call', 'speedup', 'compile s',
             'worst rel', 'ulp', 'lost', 'verdict'))
    print('-' * 110)
    for which in WHICH:
        base = rows[which]['numpy_elem_us']
        print('%-16s %-3s %10.1f %8s %10s' % ('numpy (element)', which,
                                             base, '1x', '-'))
        for (name, w), m in meta.items():
            if w != which:
                continue
            if 'error' in m:
                print('%-16s %-3s %10s  %s' % (name, which, 'FAILED',
                                               m['error'][:60]))
                continue
            sp = base / m['us']
            ws = sweeps.get(which, {}).get('worst', {}).get(name)
            vd = verdict(sp, ws) if which == 'G' else ''
            print('%-16s %-3s %10.2f %7.0fx %10.1f %12s %8s %8s  %s'
                  % (name, which, m['us'], sp, m['compile_s'],
                     ('%.2e' % ws['rel']) if ws else '-',
                     ws['ulp'] if ws else '-', ws['lost'] if ws else '-',
                     vd))
            rows[which][name] = dict(m, speedup=sp, worst=ws)
    print()
    for which, s in sweeps.items():
        print('%s sweep: %d points, numpy raised %d, numpy non-finite at %d '
              'points, %d finite reference entries, %.0f s'
              % (which, s['points'], s['ref_raised'],
                 s['ref_nonfinite_points'], s['ref_finite_entries'],
                 s['elapsed_s']))
        for name, w in s['worst'].items():
            print('  %-14s worst rel %.3e (abs %.3e, %d ulp) at %s'
                  % (name, w['rel'], w['abs'], w['ulp'], w['where']))
            print('  %-14s lost %d  gained %d  both-nonfinite %d'
                  % ('', w['lost'], w['gained'], w['both_nonfinite']))
            for lw in w['lost_where']:
                print('      LOST at %s' % lw)
    result['rows'] = rows
    result['sweeps'] = sweeps
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(result, fh, indent=1, default=str)
        print('wrote %s' % args.json)
    return 0


if __name__ == '__main__':
    sys.exit(main())
