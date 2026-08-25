"""What does a model cost to compile, and to evaluate, as it gets bigger?

`hdl_overhead.py` answers a different question -- how a generated element
compares to a HAND-WRITTEN one on the same device.  That comparison stops
being available exactly where it starts to matter: there is no
hand-written PSP103 to compare against.

This is the ladder instead.  The same two measurements (compile once,
evaluate per Newton iteration) taken across the whole range of model
complexity the DSL is asked to cover, from a resistor to a
359-parameter surface-potential MOSFET.  The point is the SHAPE of the
growth, not any single row:

  - `i` and `G` grow at different rates.  For small elements the
    Jacobian is roughly the cost of the current; for the MOSFET it is
    ~19x, because `G` is 4x4 partial derivatives of an expression that
    already took a millisecond once.  `G/i` is therefore the column to
    read, and it is the one a code-generator change has to move.

  - `chained` records which of the two code generators ran.  A model
    that calls `var()` takes the chained path (`hdl.py:1701`), which
    does not produce a `pure_spec` -- so `solve_batched` cannot take it.
    Every production model is chained.  See doc/hdl_roadmap_260824.md
    section 2.

This exists to be a BASELINE.  The roadmap gates its Phase 4 (a pure
form for the chained path) on a measurement rather than on an argument,
and this is that measurement -- run it before the change and after, and
compare the same rows.

Costs are dominated by the MOSFET's ~90 s compile, so that row is built
ONCE with no repetition, and is skipped entirely without the IHP PDK.

Run:  python benchmarks/hdl_model_cost.py
      python benchmarks/hdl_model_cost.py --json out.json
      python benchmarks/hdl_model_cost.py --no-psp     (seconds, not minutes)
"""

import argparse
import json
import os
import sys
import time
import warnings

import numpy as np


PDK = os.path.expanduser(
    '~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models')

T27 = 273.15 + 27.0


def _bench(fn, n, warmup=3):
    for _ in range(warmup):
        fn()
    t0 = time.perf_counter()
    for _ in range(n):
        fn()
    return (time.perf_counter() - t0) / n


def _time_once(fn):
    """One shot, no warm-up -- for things that cost a minute to build."""
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def _shape(cls):
    """What the compiler decided about this class."""
    info = cls._hdl_info
    return dict(chained=bool(info['chained']),
                pure=info['pure_spec'] is not None,
                const_G=bool(info['const_G']),
                nodes=len(info['terminalnames']) + len(info['internalnames']),
                params=len(info['paramnames']))


def _measure(name, build, bias, reps, compile_reps=3):
    """Compile it, then time `i` and `G` at one bias."""
    if compile_reps > 1:
        cls, _ = _time_once(build)               # warm any import cost
        t_compile = _bench(build, compile_reps, warmup=0)
    else:
        cls, t_compile = _time_once(build)

    e = cls
    ## Without this the instance re-resolves its parameter vector on every
    ## call and every row is ~10x too slow -- which is how the first run of
    ## this file disagreed with hdl_overhead.py by an order of magnitude.
    e.update_iparv()
    row = dict(name=name, compile_s=t_compile)
    row.update(_shape(type(e)))

    x = np.asarray(bias, float)
    t_i = _bench(lambda: e.i(x), reps)
    t_G = _bench(lambda: e.G(x), reps)
    row.update(i_us=t_i * 1e6, G_us=t_G * 1e6,
               ratio=(t_G / t_i if t_i else float('nan')))
    return row


def ladder(with_psp=True):
    import sympy
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       var, expl)
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric
    rows = []

    ## --- rung 1: linear, eager path, constant Jacobian ------------------
    def build_r():
        class _R(Behavioural):
            instparams = [Parameter(name='rr', desc='R', unit='ohm',
                                    default=1e3)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, b.V / rr)          # noqa: F821
        return _R('p', 'n', rr=1e3)

    ## --- rung 2: nonlinear with a BARE exp ------------------------------
    ## The reference point, and the one hdl_overhead.py measures.  No
    ## production model may use this -- a bare exp overflows -- but it is
    ## the floor the safe primitive is charged against.
    def build_d_bare():
        class _Db(Behavioural):
            instparams = [Parameter(name='iss', desc='IS', unit='A',
                                    default=1e-13),
                          Parameter(name='vtt', desc='VT', unit='V',
                                    default=0.02585)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(
                    b.I, iss * (sympy.exp(b.V / vtt) - 1))   # noqa: F821
        return _Db('p', 'n')

    ## --- rung 3: the same diode through `expl` --------------------------
    ## Differs from rung 2 by exactly one thing: the exponential is the
    ## range-safe one.  So this pair prices the SAFETY, with the physics,
    ## the node count and the code generator all held fixed.
    def build_d():
        class _D(Behavioural):
            instparams = [Parameter(name='iss', desc='IS', unit='A',
                                    default=1e-13),
                          Parameter(name='vtt', desc='VT', unit='V',
                                    default=0.02585)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, iss * (expl(b.V / vtt) - 1))  # noqa: F821
        return _D('p', 'n')

    ## --- rung 4: the SAME diode written with var() ----------------------
    ## The only difference is the let-chain, so this row isolates what
    ## taking the chained code generator costs, with the physics held
    ## fixed.  It is the cleanest read on the section-2 finding.
    def build_d_chained():
        class _Dc(Behavioural):
            instparams = [Parameter(name='iss', desc='IS', unit='A',
                                    default=1e-13),
                          Parameter(name='vtt', desc='VT', unit='V',
                                    default=0.02585)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                z = var(b.V / vtt, 'z')                      # noqa: F821
                return Contribution(b.I, iss * (expl(z) - 1))  # noqa: F821
        return _Dc('p', 'n')

    ## --- rung 4: a multi-node nonlinear element -------------------------
    def build_bjt():
        class _Q(Behavioural):
            instparams = [Parameter(name='iss', desc='IS', unit='A',
                                    default=1e-16),
                          Parameter(name='bf', desc='beta', unit='',
                                    default=100.0),
                          Parameter(name='vtt', desc='VT', unit='V',
                                    default=0.02585)]

            @staticmethod
            def analog(c, b, e):
                bbe, bbc = Branch(b, e), Branch(b, c)
                ibe = var(iss * (expl(bbe.V / vtt) - 1), 'ibe')   # noqa: F821
                ibc = var(iss * (expl(bbc.V / vtt) - 1), 'ibc')   # noqa: F821
                ict = var(ibe - ibc, 'ict')
                return (Contribution(Branch(c, e).I, ict),
                        Contribution(bbe.I, ibe / bf),            # noqa: F821
                        Contribution(bbc.I, ibc / bf))            # noqa: F821
        return _Q('c', 'b', 'e')

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        rows.append(_measure('R (linear, eager)', build_r, [0.3, 0.0], 20000))
        rows.append(_measure('diode, bare exp (eager)', build_d_bare,
                             [0.4, 0.0], 20000))
        rows.append(_measure('diode, expl (eager)', build_d,
                             [0.4, 0.0], 20000))
        rows.append(_measure('diode, expl + var (chained)', build_d_chained,
                             [0.4, 0.0], 20000))
        rows.append(_measure('BJT (3 terminals, chained)', build_bjt,
                             [0.0, 0.7, 0.0], 5000))

        if with_psp:
            psp = psp_row()
            if psp is not None:
                rows.append(psp)
    return rows


def psp_row():
    """The flagship: 359 parameters, ~90 s to compile.  Built once."""
    if not os.path.isdir(PDK):
        print('  (PSP row skipped: IHP Open PDK not at %s)\n' % PDK)
        return None

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

    print('  building the PSP MOSFET (this takes ~90 s) ...', flush=True)

    def build():
        ## Import inside, so the class-definition cost is counted.
        from pycircuit.circuit.compact import PspMosLongChannel
        e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        e.update_iparv()
        return e

    e, t_compile = _time_once(build)
    row = dict(name='PSP103 nmos (chained)', compile_s=t_compile)
    row.update(_shape(type(e)))

    ## Saturation, well above the leakage floor -- a bias the model is
    ## actually used at, not a corner.
    x = np.asarray(e.bias(1.2, 0.8, 0.0, 0.0), float)
    t_i = _bench(lambda: e.i(x), 200, warmup=2)
    t_G = _bench(lambda: e.G(x), 50, warmup=2)
    row.update(i_us=t_i * 1e6, G_us=t_G * 1e6, ratio=t_G / t_i)
    return row


def dispatch_costs():
    """Where the evaluation time actually goes: printer dispatch.

    The Newton inner loop calls an element at a SCALAR bias, and on
    scalars numpy's array machinery is almost all overhead.  These are
    the primitives sympy's numpy printer can emit, priced per call:

    `numpy.select` is the expensive one, and it is what sympy lowers
    `Piecewise`, `Min`, `Max` AND `Heaviside` to.  `numpy.minimum` has
    the same semantics as a two-argument `Min` and costs ~9x less, so
    every Min/Max in a model pays a dispatch penalty that has nothing
    to do with its arithmetic.

    This is why a 5-operation hand-written `Piecewise` measures the
    same as a 130-operation one: you are timing the dispatch, not the
    expression.  Op counts are the wrong proxy for cost here, which is
    worth knowing before optimising any of this symbolically.
    """
    prims = [
        ('numpy.exp', lambda: np.exp(0.4)),
        ('numpy.minimum', lambda: np.minimum(1.0, 2.0)),
        ('numpy.heaviside', lambda: np.heaviside(0.4, 0.5)),
        ('numpy.where', lambda: np.where(True, 1.0, 2.0)),
        ('numpy.amin', lambda: np.amin([1.0, 2.0])),
        ('numpy.select', lambda: np.select([True, False], [1.0, 2.0])),
    ]
    print('Per-call dispatch cost on SCALAR inputs:')
    for name, fn in prims:
        print('  %-18s %6.2f us' % (name, _bench(fn, 20000, warmup=500) * 1e6))
    print()


def emitted(expr, wrt):
    """Which numpy primitives does this expression's derivative emit?

    The instrument for the above: a primitive is only as cheap as what
    the printer turns it into, and that is invisible from the sympy
    side.  Assert on this, not on the value alone.
    """
    import inspect
    import re
    import sympy
    f = sympy.lambdify([wrt], sympy.diff(expr, wrt), modules=['numpy'],
                       cse=True)
    src = inspect.getsource(f)
    return {n: len(re.findall(r'\b' + n + r'\(', src))
            for n in ('select', 'amin', 'amax', 'minimum', 'maximum',
                      'heaviside', 'where', 'exp')
            if re.search(r'\b' + n + r'\(', src)}


def primitive_table():
    """What each kernel primitive costs, and what it emits."""
    import sympy
    from pycircuit.circuit.hdl import expl, hypsmooth, safe_sqrt, safe_ln

    v = sympy.Symbol('v', real=True)
    cases = [('bare exp', sympy.exp(v / 0.02585)),
             ('expl', expl(v / 0.02585)),
             ('hypsmooth', hypsmooth(v, 1e-6)),
             ('safe_sqrt', safe_sqrt(v)),
             ('safe_ln', safe_ln(v))]

    print('%-14s %9s %9s   %s' % ('primitive', 'val us', 'deriv us',
                                  'emitted by d/dv'))
    for name, e in cases:
        try:
            fv = sympy.lambdify([v], e, modules=['numpy'], cse=True)
            fd = sympy.lambdify([v], sympy.diff(e, v), modules=['numpy'],
                                cse=True)
            tv = _bench(lambda: fv(0.4), 5000, warmup=200) * 1e6
            td = _bench(lambda: fd(0.4), 5000, warmup=200) * 1e6
            em = ', '.join('%s x%d' % (k, n)
                           for k, n in sorted(emitted(e, v).items()))
            print('%-14s %9.2f %9.2f   %s' % (name, tv, td, em))
        except Exception as exc:                      # pragma: no cover
            print('%-14s   (skipped: %s)' % (name, exc))
    print()


def report(rows):
    print()
    print('%-28s %10s %10s %10s %7s  %-7s %s'
          % ('model', 'compile', 'i (us)', 'G (us)', 'G/i', 'path', 'pure?'))
    print('-' * 90)
    for r in rows:
        comp = ('%8.1f ms' % (r['compile_s'] * 1e3) if r['compile_s'] < 10
                else '%8.1f s ' % r['compile_s'])
        print('%-28s %10s %10.2f %10.2f %6.1fx  %-7s %s'
              % (r['name'], comp, r['i_us'], r['G_us'], r['ratio'],
                 'chained' if r['chained'] else 'eager',
                 'yes' if r['pure'] else 'NO'))
    print()

    ## The number that decides whether the DSL is usable at circuit scale.
    print('At circuit scale -- one Newton iteration, 100 devices:')
    for r in rows:
        print('  %-28s %10.2f ms' % (r['name'], r['G_us'] * 100 / 1e3))
    print()

    nopure = [r['name'] for r in rows if r['chained'] and not r['pure']]
    if nopure:
        print('No pure_spec (so solve_batched cannot take these):')
        for n in nopure:
            print('  - %s' % n)
        print()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--no-psp', action='store_true',
                    help='skip the ~90 s MOSFET compile')
    ap.add_argument('--dispatch-only', action='store_true',
                    help='just the printer-dispatch and primitive tables')
    ap.add_argument('--json', metavar='PATH',
                    help='also write the rows as JSON, for later comparison')
    args = ap.parse_args()

    dispatch_costs()
    primitive_table()
    if args.dispatch_only:
        return 0
    rows = ladder(with_psp=not args.no_psp)
    report(rows)

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(dict(rows=rows), fh, indent=2)
        print('wrote %s' % args.json)
    return 0


if __name__ == '__main__':
    sys.exit(main())
