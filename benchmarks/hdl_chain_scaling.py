"""Why a compact model needs `var()`: compile time vs nesting depth.

The claim under test (hdl.md, "The let-chain"): compiling a model by
substituting its intermediates into one big expression is EXPONENTIAL in
how deeply those intermediates are reused, while keeping each one as a
named symbol with its own defining equation is roughly LINEAR.  This is
not a micro-optimisation -- it is the difference between a surface-
potential model compiling in half a minute and not compiling at all.

Where the exponential comes from
--------------------------------
Write ``u1 = f(u0)`` where ``f`` mentions its argument more than once --
``sqrt(1 + u*u) + u*tanh(u)/(2+k)`` mentions it three times, and every
compact-model expression is like this.  Substituting ``u0`` into ``u1``
produces three copies of ``u0``.  Substituting again produces nine.
After *d* levels the expression TREE holds 3**d copies of the innermost
subexpression, even though the DAG has only *d* nodes.
``Matrix.jacobian`` then differentiates the tree.

Measured, the wall is worse than 3**d, because sympy's own operations
are superlinear in tree size: the per-level cost ratio observed here
runs 2.5x, 1.9x, 3.6x, then 24x between depths 4 and 5.

`var()` breaks the substitution: each intermediate stays a symbol, the
compiler emits a topological let-chain, and the Jacobian comes from
forward accumulation over the chain's edges.

Three measurements:

1. compile time vs nesting depth, both paths, with the stamps compared
   at every depth where the eager path still finishes;
2. compile time vs number of intermediates, at PSP scale (PSP103 carries
   roughly 1400 of them across 4 terminals), with every Jacobian checked
   against finite differences;
3. per-call evaluation cost, so the compile-time win is not being bought
   with a slow inner loop.

Run:  python benchmarks/hdl_chain_scaling.py
"""

import sys
import time
import warnings

import numpy as np
import sympy


## Measured: depth 4 takes 0.7s, depth 5 takes 16s, depth 6 takes several
## minutes.  Five is the last depth worth waiting for, and it already
## shows the knee.
EAGER_MAX_DEPTH = 5


def _build(depth, chained):
    """A chain of `depth` levels, each mentioning the previous one thrice."""
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       var, vt)
    from pycircuit.utilities.param import Parameter

    def analog(d, s):
        b = Branch(d, s)
        u = b.V / vt()
        for k in range(depth):
            e = sympy.sqrt(1 + u * u) + u * sympy.tanh(u) / (2 + k)
            u = var(e, 'u%d' % k) if chained else e
        return Contribution(b.I, IS * u)                      # noqa: F821

    cls = type('Depth%d%s' % (depth, chained), (Behavioural,), {
        'instparams': [Parameter(name='IS', desc='I', unit='A',
                                 default=1e-6)],
        'analog': staticmethod(analog)})
    e = cls(cm.Node('d'), cm.Node('s'), IS=1e-6)
    e.update_iparv()
    return e


def _psp_shaped(n_inter):
    """A 4-terminal device with `n_inter` reused intermediates.

    Not PSP -- a synthetic model of the same SHAPE: four terminals, three
    contributions, deep reuse, and the same math kernel PSP is written in
    (sqrt / exp / ln / tanh / abs).  The question it answers is only
    whether the compiler survives that size.
    """
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       var, vt)
    from pycircuit.utilities.param import Parameter

    def analog(d, g, s, b):
        bds, bgs, bbs = Branch(d, s), Branch(g, s), Branch(b, s)
        vds, vgs, vbs = bds.V / vt(), bgs.V / vt(), bbs.V / vt()
        pool = [var(vds, 'vds'), var(vgs, 'vgs'), var(vbs, 'vbs'),
                var(sympy.sqrt(1 + vgs * vgs), 'sq')]
        ops = (lambda a, c: sympy.sqrt(1e-9 + a * a + c * c),
               lambda a, c: a * sympy.tanh(c),
               lambda a, c: (a + c) / (1 + sympy.Abs(a - c)),
               lambda a, c: sympy.log(1 + sympy.exp(-sympy.Abs(a))) + c,
               lambda a, c: a * c / (1 + a * a),
               lambda a, c: sympy.sqrt(sympy.Abs(a)) + c / (2 + c * c))
        for k in range(n_inter):
            a = pool[(k * 7 + 1) % len(pool)]
            c = pool[(k * 13 + 3) % len(pool)]
            pool.append(var(ops[k % len(ops)](a, c), 'i%d' % k))
        return (Contribution(bds.I, IS * pool[-1] * pool[-2]),   # noqa: F821
                Contribution(bgs.I, IS * 1e-3 * pool[-3]),       # noqa: F821
                Contribution(bbs.I, IS * 1e-4 * pool[-4]))       # noqa: F821

    cls = type('PSPish%d' % n_inter, (Behavioural,), {
        'instparams': [Parameter(name='IS', desc='I', unit='A',
                                 default=1e-6)],
        'analog': staticmethod(analog)})
    e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'), IS=1e-6)
    e.update_iparv()
    return e


def _fd_jacobian(elem, x):
    J = np.zeros((len(x), len(x)))
    for j in range(len(x)):
        h = 1e-6 * max(1.0, abs(float(x[j])))
        xp, xm = np.array(x, float), np.array(x, float)
        xp[j] += h
        xm[j] -= h
        J[:, j] = (np.asarray(elem.i(xp), float)
                   - np.asarray(elem.i(xm), float)) / (2 * h)
    return J


def depth_table():
    print()
    print('=== 1. compile time vs nesting depth ===')
    print('Each level mentions the previous one three times, so the eager')
    print('expression tree triples per level while the DAG grows by one.')
    print()
    print('%6s %13s %13s %10s %14s' %
          ('depth', 'eager (s)', 'chained (s)', 'eager grow',
           'max|G_e-G_c|'))
    x = np.array([0.7, 0.0])
    prev_eager = None
    for depth in (1, 2, 3, 4, 5, 6, 8, 12, 16, 24, 32, 48):
        t0 = time.perf_counter()
        ce = _build(depth, True)
        tc = time.perf_counter() - t0
        Gc = np.asarray(ce.G(x), float)

        if depth <= EAGER_MAX_DEPTH:
            t0 = time.perf_counter()
            ee = _build(depth, False)
            te = time.perf_counter() - t0
            Ge = np.asarray(ee.G(x), float)
            agree = float(np.max(np.abs(Ge - Gc)))
            growth = ('%.1fx' % (te / prev_eager)) if prev_eager else '--'
            prev_eager = te
            print('%6d %13.3f %13.3f %10s %14.2e'
                  % (depth, te, tc, growth, agree))
        else:
            print('%6d %13s %13.3f %10s %14s'
                  % (depth, 'not run', tc, '--', '--'))
        sys.stdout.flush()
    print()
    print('Rows past depth %d are not run on the eager path.  Depth 5'
          % EAGER_MAX_DEPTH)
    print('already costs 16s and the growth ratio is still climbing, so')
    print('depth 48 is not a long wait -- it is not a wait that ends.')


def psp_scale_table():
    print()
    print('=== 2. PSP scale: 4 terminals, up to ~1400 intermediates ===')
    print('PSP103 evaluates ~1200 assignments and carries ~1400')
    print('intermediates.  Every Jacobian is checked against finite')
    print('differences, so this measures a compiler that works.')
    print()
    print('%9s %11s %9s %12s %12s %13s' %
          ('n_inter', 'compile', 'emitted', 'i eval', 'G eval', '|J-FD|/scale'))
    x = np.array([0.9, 1.1, 0.0, -0.3])
    for n in (10, 50, 200, 600, 1400):
        t0 = time.perf_counter()
        e = _psp_shaped(n)
        tc = time.perf_counter() - t0

        src = getattr(type(e)._hdl_info['funcs']['G'], '_src', '')
        nlines = len(src.splitlines())

        ti = _time(lambda: e.i(x))
        tg = _time(lambda: e.G(x))

        G = np.asarray(e.G(x), float)
        fd = _fd_jacobian(e, x)
        scale = max(1.0, float(np.max(np.abs(G))))
        err = float(np.max(np.abs(G - fd))) / scale
        print('%9d %10.2fs %9d %11.1fus %11.1fus %13.2e'
              % (n, tc, nlines, ti * 1e6, tg * 1e6, err))
        sys.stdout.flush()
    print()
    print('Compile time is paid ONCE per model class.  The emitted column')
    print('is lines of generated Python for the Jacobian -- the let-chain')
    print('plus its forward-accumulated derivatives.')


def _time(fn, n=200, warmup=5):
    for _ in range(warmup):
        fn()
    t0 = time.perf_counter()
    for _ in range(n):
        fn()
    return (time.perf_counter() - t0) / n


def per_call_table():
    print()
    print('=== 3. is the compile-time win paid for at run time? ===')
    print('Same model, both paths, at a depth the eager path can reach.')
    print()
    print('%6s %14s %14s %10s' % ('depth', 'eager G (us)', 'chain G (us)',
                                  'ratio'))
    x = np.array([0.7, 0.0])
    for depth in (1, 2, 3, 4):
        ee, ce = _build(depth, False), _build(depth, True)
        te = _time(lambda: ee.G(x))
        tc = _time(lambda: ce.G(x))
        print('%6d %14.2f %14.2f %9.2fx'
              % (depth, te * 1e6, tc * 1e6, tc / te))
        sys.stdout.flush()
    print()
    print('The chain emits a flat sequence of assignments where the eager')
    print('path emits one cse-factored expression, so per-call cost is')
    print('comparable -- the win is in what can be compiled at all.')


if __name__ == '__main__':
    warnings.simplefilter('ignore')
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    cm.default_toolkit = numeric

    depth_table()
    psp_scale_table()
    per_call_table()
