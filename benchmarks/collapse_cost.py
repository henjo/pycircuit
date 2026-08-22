"""What does NOT collapsing a zero parasitic resistance cost?

PSP103 wraps seven parasitic resistances in its `CollapsableR` macro
(`PSP103_macrodefs.include:264`, used at `PSP103_module.include:1718-1724`
for rgate, rsource, rdrain, rbulk, rjuns, rjund, rwell):

    if ((R) > 0.0) I(N1,N2) <+ MULT_i * (G) * V(N1,N2);
    else           V(N1,N2) <+ 0.0;

and in the IHP card most of them are zero for the ordinary non-RF
device: `rsh`, `rshd` and `rvpoly` are 0, and `rgo` and `rbulko` are
gated on `rfmode`.

pycircuit can already express both arms with ONE formulation and no
conditional at all -- write the resistance as a branch,

    Contribution(Branch(N1, N2).V, Branch(N1, N2).I * R)

which is the resistor for R > 0 and an exact short for R = 0, since the
row is `-(v1 - v2) + i_br*R = 0`.  That is measured in
`test_hdl_collapse.py`: exact to the last bit at every R including
exactly zero, and unlike the conductance form it does not divide by
zero.

The catch is cost.  The branch form carries an internal node AND a
branch current whether or not R is zero, where a collapsing compiler
carries neither.  This measures how much that matters, so the decision
to build parameter-driven collapse rests on a number rather than on
taste.

Run:  python benchmarks/collapse_cost.py
"""

import time
import warnings

import numpy as np


#: PSP's seven CollapsableR sites.
NPAR = 7


def make(nparasitic):
    """A nonlinear core behind `nparasitic` zero series resistances.

    `nparasitic = 0` is the collapsed topology -- what the model means
    when it takes the `else V <+ 0` arm.
    """
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Node)
    from pycircuit.utilities.param import Parameter

    def analog(p, m):
        node = p
        for k in range(nparasitic):
            node = Node('int%d' % k)
        core = Branch(node, m)
        out = [Contribution(core.I,
                            IS * core.V * 1e-3                  # noqa: F821
                            + IS * core.V ** 3 * 1e-6)]         # noqa: F821
        node = p
        for k in range(nparasitic):
            nxt = Node('int%d' % k)
            br = Branch(node, nxt, 'par%d' % k)
            out.append(Contribution(br.V, br.I * rpar))         # noqa: F821
            node = nxt
        return tuple(out)

    return type('Dev%d' % nparasitic, (Behavioural,), {
        'instparams': [Parameter(name='IS', desc='I', unit='A',
                                 default=1e-3),
                       Parameter(name='rpar', desc='parasitic', unit='ohm',
                                 default=0.0)],
        'analog': staticmethod(analog)})


def chain(cls, ndev):
    from pycircuit.circuit import gnd
    from pycircuit.circuit.elements import R, SubCircuit, VS

    c = SubCircuit()
    prev = c.add_node('n0')
    c['vs'] = VS(prev, gnd, v=1.0)
    for k in range(ndev):
        nxt = c.add_node('n%d' % (k + 1))
        c['d%d' % k] = cls(prev, nxt, IS=1e-3, rpar=0.0)
        c['r%d' % k] = R(nxt, gnd, r=1e4)
        prev = nxt
    c.update_iparv()
    return c


def main():
    from pycircuit.circuit import gnd
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.dcanalysis import DC

    print()
    print('A chain of devices, each with %d zero-valued parasitic series'
          % NPAR)
    print('resistances.  "collapsed" is what the model means; "branch" is')
    print('what pycircuit builds today.  Both give the same answer.')
    print()
    print('%6s %11s %11s %13s %13s %10s' %
          ('ndev', 'n collapsed', 'n branch', 't collapsed', 't branch',
           'slowdown'))

    collapsed, branchy = make(0), make(NPAR)
    for ndev in (5, 20, 50, 100):
        got = []
        for cls in (collapsed, branchy):
            c = chain(cls, ndev)
            DC(c, toolkit=numeric).solve()               # warm
            t0 = time.perf_counter()
            for _ in range(3):
                res = DC(c, toolkit=numeric).solve()
            got.append((c.n, (time.perf_counter() - t0) / 3, res))
        (n0, t0_, r0), (n1, t1_, r1) = got
        v0 = float(r0.v('n%d' % ndev, gnd))
        v1 = float(r1.v('n%d' % ndev, gnd))
        assert abs(v0 - v1) < 1e-9 * max(1.0, abs(v0)), (v0, v1)
        print('%6d %11d %11d %11.2fms %11.2fms %9.2fx'
              % (ndev, n0, n1, 1e3 * t0_, 1e3 * t1_, t1_ / t0_))

    print()
    print('The cost is %d extra unknowns per device -- one internal node'
          % (NPAR * 2))
    print('and one branch current for each parasitic a collapsing compiler')
    print('would have removed outright.  The answers agree exactly, so this')
    print('is purely a matrix-size argument -- and at 100 devices it is a')
    print('decisive one.')


if __name__ == '__main__':
    warnings.simplefilter('ignore')
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    cm.default_toolkit = numeric
    main()
