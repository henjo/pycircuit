"""Stage-predictor benchmark: a state-free exponential nonlinearity.

The `Diode` element is useless as an instrument here -- it linearises `G`
around a STORED `_vlim`, so its Newton is seed-independent by construction
(measured: an all-zeros seed gives the same iteration histogram as an exact
one).  This element is a plain behavioural exponential: `i` and `G` are both
functions of the passed `x`, no state, no limiting.
"""
import warnings, numpy as np, sympy
warnings.simplefilter('ignore')
import pycircuit.circuit.circuit
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import C, R, VSin
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, Parameter
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import transient as TR
from pycircuit.circuit import nrsolver as NR
from pycircuit.circuit.integrator import (GLM2Integrator, GLM3Integrator,
                                          GLM4Integrator, RadauIIA3Integrator,
                                          ESDIRK43Integrator, TRBDF2Integrator,
                                          Gear2Integrator, TrapezoidalIntegrator)

pycircuit.circuit.circuit.default_toolkit = numeric


class ExpG(Behavioural):
    """i = IS (exp(V/VT) - 1), state-free, no limiting."""
    instparams = [Parameter(name='IS', desc='sat', unit='A', default=1e-12),
                  Parameter(name='VT', desc='thermal', unit='V', default=0.026)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return (Contribution(b.I, IS * (sympy.exp(b.V / VT) - 1)),)  # noqa: F821


PER = 1e-3


def fixture(va=0.8):
    c = SubCircuit(); c.add_node('a'); c.add_node('b')
    c['vs'] = VSin('a', gnd, va=va, freq=1.0 / PER)
    c['rs'] = R('a', 'b', r=50.0)
    c['nl'] = ExpG('b', gnd, IS=1e-12, VT=0.026)
    c['cl'] = C('b', gnd, c=1e-9)
    c['rl'] = R('b', gnd, r=1e4)
    return c


def measure(integ, mode, N, va=0.8):
    log = []
    solve_orig = NR.StandardNewton.solve_system

    def wrapped(self, x0, eval_FJ, *a, **k):
        x, it = solve_orig(self, x0, eval_FJ, *a, **k)
        log.append((float(np.max(np.abs(np.asarray(x0) - np.asarray(x)))), it))
        return x, it
    NR.StandardNewton.solve_system = wrapped
    TR.Transient.glm_predictor = mode
    n = {'i': 0}
    cir = fixture(va)
    fi = cir.i
    cir.i = lambda *a, **k: (n.__setitem__('i', n['i'] + 1), fi(*a, **k))[1]
    tr = Transient(cir, integrator=integ, reltol=1e-9)
    try:
        r = tr.solve(refnode=gnd, tend=PER, timestep=PER / N, fixed_timestep=True)
    finally:
        NR.StandardNewton.solve_system = solve_orig
    it = np.array([l[1] for l in log], dtype=int)
    d = np.array([l[0] for l in log])
    ns = tr.statistics.accepted_steps
    return dict(iters=int(it.sum()), mean=float(it.mean()), steps=ns,
                icalls=n['i'], seed_err=float(np.median(d)),
                x=np.asarray(r.v('b'), dtype=float))


def table():
    """The measurement the predictor's docstring quotes."""
    from pycircuit.circuit import nrsolver as NR
    print('%-5s %5s %4s %-8s %11s %8s %7s %8s %8s'
          % ('meth', 'va', 'N', 'pred', 'worst seed', 'mean it',
             'max it', 'i-calls', 'vs none'))
    for cls, name in ((GLM3Integrator, 'glm3'), (GLM4Integrator, 'glm4')):
        for va, N in ((2.0, 40), (0.8, 40), (0.8, 200), (2.0, 200)):
            base = None
            for mode in ('none', 'local'):
                log = []
                so = NR.StandardNewton.solve_system

                def w(self, x0, ev, *a, **k):
                    x, it = so(self, x0, ev, *a, **k)
                    log.append((float(np.max(np.abs(np.asarray(x0)
                                                    - np.asarray(x)))), it))
                    return x, it
                NR.StandardNewton.solve_system = w
                TR.Transient.glm_predictor = mode
                n = {'i': 0}
                cir = fixture(va)
                fi = cir.i
                cir.i = lambda *a, **k: (n.__setitem__('i', n['i'] + 1),
                                         fi(*a, **k))[1]
                tr = Transient(cir, integrator=cls(), reltol=1e-9)
                try:
                    tr.solve(refnode=gnd, tend=PER, timestep=PER / N,
                             fixed_timestep=True)
                finally:
                    NR.StandardNewton.solve_system = so
                d = np.array([l[0] for l in log])
                it = np.array([l[1] for l in log])
                if base is None:
                    base = n['i']
                print('%-5s %5.1f %4d %-8s %11.3e %8.3f %7d %8d %7.1f%%'
                      % (name, va, N, mode, d.max(), it.mean(), it.max(),
                         n['i'], 100.0 * (n['i'] - base) / base))
            print()
    TR.Transient.glm_predictor = 'local'


if __name__ == '__main__':
    table()
