"""What does the Behavioural HDL cost against hand-written elements?

The claim under test (hdl.md, "Cost"): a compiled contribution equation is
not a slow interpreter -- ``sympy.lambdify(..., cse=True)`` emits a plain
numpy function, so a generated element should land within a small factor
of a hand-written stamp, and the gap should come from *shape* (one call
per element returning a dense vector, vs a stored matrix returned by
reference) rather than from expression walking.

Three measurements, all on the same circuits:

1. per-call element evaluation (i/G/q/C) -- the inner loop of Newton;
2. compile time -- what you pay ONCE, at class definition;
3. an end-to-end RC transient, HDL elements vs elements.py elements,
   with the waveforms compared so the timing is of equal work.

Bounded on purpose, like phase0_baseline: one circuit per row, single
timed repetition after a warm-up -- a yardstick, not a benchmark suite.

Run:  python benchmarks/hdl_overhead.py
"""

import time
import warnings

import numpy as np


def _bench(fn, n, warmup=3):
    for _ in range(warmup):
        fn()
    t0 = time.perf_counter()
    for _ in range(n):
        fn()
    return (time.perf_counter() - t0) / n


def per_call_table():
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.elements import R, C, Diode
    from pycircuit.circuit import elements_hdl as eh

    cm.default_toolkit = numeric
    rows = []
    cases = [
        ('R  i', R('p', 'n', r=1e3), eh.RHdl('p', 'n', r=1e3), 'i',
         np.array([0.3, 0.0])),
        ('R  G', R('p', 'n', r=1e3), eh.RHdl('p', 'n', r=1e3), 'G',
         np.array([0.3, 0.0])),
        ('C  q', C('p', 'n', c=1e-9), eh.CHdl('p', 'n', c=1e-9), 'q',
         np.array([0.3, 0.0])),
        ('C  C', C('p', 'n', c=1e-9), eh.CHdl('p', 'n', c=1e-9), 'C',
         np.array([0.3, 0.0])),
        ('D  i', Diode('p', 'n', IS=1e-13),
         eh.DiodeHdl('p', 'n', IS=1e-13), 'i', np.array([0.4, 0.0])),
        ('D  G', Diode('p', 'n', IS=1e-13),
         eh.DiodeHdl('p', 'n', IS=1e-13), 'G', np.array([0.4, 0.0])),
    ]
    for label, ref, hdl_el, meth, x in cases:
        ref.update_iparv(); hdl_el.update_iparv()
        t_ref = _bench(lambda: getattr(ref, meth)(x), 20000)
        t_hdl = _bench(lambda: getattr(hdl_el, meth)(x), 20000)
        rows.append((label, t_ref * 1e6, t_hdl * 1e6, t_hdl / t_ref))

    print('%-6s %12s %12s %8s' % ('call', 'hand (us)', 'hdl (us)', 'ratio'))
    for label, a, b, r in rows:
        print('%-6s %12.3f %12.3f %7.2fx' % (label, a, b, r))
    print()


def compile_time():
    """The once-per-class cost: symbolic assembly + Jacobians + lambdify."""
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, ddt

    def build_resistor():
        class _R(Behavioural):
            instparams = [Parameter(name='rr', desc='R', unit='ohm',
                                    default=1e3)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, b.V / rr)          # noqa: F821
        return _R

    def build_diode():
        class _D(Behavioural):
            instparams = [Parameter(name='iss', desc='IS', unit='A',
                                    default=1e-13),
                          Parameter(name='vtt', desc='VT', unit='V',
                                    default=0.02585)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(
                    b.I, iss * (sympy.exp(b.V / vtt) - 1))  # noqa: F821
        return _D

    for name, fn in (('2-terminal linear', build_resistor),
                     ('2-terminal nonlinear', build_diode)):
        t = _bench(fn, 3, warmup=1)
        print('compile %-22s %8.1f ms (once, at class definition)'
              % (name, t * 1e3))
    print()


def transient_end_to_end():
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import gnd
    from pycircuit.circuit.elements import SubCircuit, VSin, R, C
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.transient import Transient

    def build(use_hdl, stages=8):
        cm.default_toolkit = numeric
        c = SubCircuit()
        prev = c.add_node('n0')
        c['vs'] = VSin(prev, gnd, va=1.0, freq=2e3)
        Rc, Cc = (eh.RHdl, eh.CHdl) if use_hdl else (R, C)
        for k in range(stages):
            nxt = c.add_node('n%d' % (k + 1))
            c['R%d' % k] = Rc(prev, nxt, r=1e3)
            c['C%d' % k] = Cc(nxt, gnd, c=1e-8)
            prev = nxt
        return c, 'n%d' % stages

    out = {}
    for use_hdl in (False, True):
        c, node = build(use_hdl)
        tran = Transient(c, toolkit=numeric, uic=True)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            t0 = time.perf_counter()
            res = tran.solve(tend=1e-3, timestep=1e-6, fixed_timestep=True)
            wall = time.perf_counter() - t0
        out[use_hdl] = (wall, np.asarray(res.v(node).y, float),
                        res.statistics.accepted_steps)

    w_ref, y_ref, n_ref = out[False]
    w_hdl, y_hdl, n_hdl = out[True]
    err = float(np.max(np.abs(y_hdl - y_ref)))
    print('8-stage RC ladder transient (%d steps both):' % n_ref)
    print('   hand-written elements : %7.3f s' % w_ref)
    print('   HDL-generated elements: %7.3f s   (%.2fx)'
          % (w_hdl, w_hdl / w_ref))
    print('   max |v_hdl - v_hand|  : %.2e  (equal work compared)' % err)
    assert n_ref == n_hdl, (n_ref, n_hdl)
    assert err < 1e-12, err
    print()


def main():
    print(__doc__.splitlines()[0])
    print()
    per_call_table()
    compile_time()
    transient_end_to_end()
    print('The generated elements are numpy functions, not an interpreter: '
          'the per-call gap is\nconstant-factor call overhead, and it buys '
          'exact symbolic Jacobians for free.')


if __name__ == '__main__':
    main()
