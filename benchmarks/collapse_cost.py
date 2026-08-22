"""Parameter-driven node collapse: what it costs not to have it.

PSP103 wraps seven parasitic resistances in its `CollapsableR` macro
(`PSP103_macrodefs.include:264`, used at `PSP103_module.include:1718-1724`
for rgate, rsource, rdrain, rbulk, rjuns, rjund, rwell):

    if ((R) > 0.0) I(N1,N2) <+ MULT_i * (G) * V(N1,N2);
    else           V(N1,N2) <+ 0.0;

and in the IHP card most of them are zero for the ordinary non-RF
device: `rsh`, `rshd` and `rvpoly` are 0, and `rgo` and `rbulko` are
gated on `rfmode`.

Three ways to build that, compared on the same circuit:

* **hand** -- the collapsed topology written out directly.  What the
  model MEANS when `R` is zero, and the target to match.
* **branch** -- both arms as one relation,
  `Contribution(Branch(N1,N2).V, Branch(N1,N2).I * R)`.  Exact for every
  R including zero (the row degenerates to `v1 = v2`), and it needs no
  conditional at all -- but it carries an internal node AND a branch
  current whether or not the resistance exists.
* **Collapse** -- the model written as PSP writes it, with
  `Collapse(br, R <= 0)` alongside the resistor contribution, letting the
  compiler take the collapsed arm.

This is what justified building `Collapse`: the branch form is correct
and, at 100 devices, 14x slower.

Run:  python benchmarks/collapse_cost.py
"""

import time
import warnings

import numpy as np


#: PSP's seven CollapsableR sites.
NPAR = 7


def make(nparasitic, declare_collapse=False):
    """A nonlinear core behind `nparasitic` zero series resistances.

    `nparasitic = 0` is the collapsed topology -- what the model means
    when it takes the `else V <+ 0` arm, hand-written.
    `declare_collapse` writes the model the way PSP does, with a
    `Collapse` per parasitic, and lets the compiler reach the same
    topology from `rpar = 0`.
    """
    from pycircuit.circuit.hdl import (Behavioural, Branch, Collapse,
                                       Contribution, Node)
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
            if declare_collapse:
                out.append(Contribution(br.I, br.V / rpar))     # noqa: F821
                out.append(Collapse(br, rpar <= 0))             # noqa: F821
            else:
                out.append(Contribution(br.V, br.I * rpar))     # noqa: F821
            node = nxt
        return tuple(out)

    return type('Dev%d%s' % (nparasitic, declare_collapse), (Behavioural,), {
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
    print('resistances, built three ways.  All give the same answer.')
    print()
    print('%6s %8s %8s %8s %11s %11s %11s' %
          ('ndev', 'n hand', 'n branch', 'n Coll.', 't hand', 't branch',
           't Collapse'))

    hand, branchy, collapsing = make(0), make(NPAR), make(NPAR, True)
    for ndev in (5, 20, 50, 100):
        got = []
        for cls in (hand, branchy, collapsing):
            c = chain(cls, ndev)
            DC(c, toolkit=numeric).solve()               # warm
            t0 = time.perf_counter()
            for _ in range(3):
                res = DC(c, toolkit=numeric).solve()
            got.append((c.n, (time.perf_counter() - t0) / 3, res))
        vs = [float(r.v('n%d' % ndev, gnd)) for _, _, r in got]
        for v in vs[1:]:
            assert abs(v - vs[0]) < 1e-9 * max(1.0, abs(vs[0])), vs
        print('%6d %8d %8d %8d %9.2fms %9.2fms %9.2fms'
              % (ndev, got[0][0], got[1][0], got[2][0],
                 1e3 * got[0][1], 1e3 * got[1][1], 1e3 * got[2][1]))

    print()
    print('All three give the same answer (asserted).  "hand" is the')
    print('collapsed topology written out by hand; "branch" carries every')
    print('parasitic as V = I*R, costing %d extra unknowns per device;'
          % (NPAR * 2))
    print('"Collapse" is the model written as PSP writes it, with the')
    print('compiler taking the collapsed arm from rpar = 0.  Collapse')
    print('should track "hand" -- that is the point of the feature.')


if __name__ == '__main__':
    warnings.simplefilter('ignore')
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    cm.default_toolkit = numeric
    main()
