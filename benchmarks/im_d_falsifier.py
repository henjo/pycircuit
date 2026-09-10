"""Does violating the `im D(t)` hypothesis actually cost anything HERE?

THE PREMISE.  Lamour, März & Tischendorf (2013) put one hypothesis under
IRK(DAE) convergence (Thm 5.7), GLM convergence at stage order (Thm 5.9 --
what `NordsieckGLMIntegrator` rests on) and contractivity transfer (Thm 6.9):
`im D(t)` time-invariant.  It is one term in one equation -- the IERODE's field
is `u' = R'(t)u + D(t)omega(u,t)` and the hypothesis exists to kill `R'(t)u`.
For charge-oriented MNA it means `im C(x)` constant along the orbit.

⚠⚠ WHY THIS FILE EXISTS.  Measured 2026-09-10, every fixture behind this
tree's GLM order, stage-predictor and noise-floor results has `C(x)` LITERALLY
CONSTANT -- `max|C(x1) - C(x2)| == 0` over random `x` -- because every
reactance in them is linear.  So those results sit inside Thm 5.9's scope
VACUOUSLY.  And the bias is not particular to them: a smoothly varying
nonlinear `C(v) > 0` is rank-1 throughout and does not exercise the condition
either, so ANY ordinary circuit fixture satisfies the hypothesis trivially and
teaches nothing about it.  `CubicCap` below is built to violate it.

⚠⚠ AND THE INSTRUMENT IS NOT AN ORDER SWEEP, which is what this file was
first going to be.  Example 3.34 and Theorem 3.53: the zeros of `det G1` split
the domain into open regularity regions, index-1 regular on each, with exactly
one solution through each consistent point INSIDE a region -- but at a
critical point "there are TWO solutions passing through". **The failure mode
at a regularity boundary is UNIQUENESS, not accuracy.**  A convergence or
order metric asks a question that has no yes there and would return a clean
null meaning nothing.

So the test is a BRANCH test.  Integrate through the crossing several times,
each run individually legitimate and differing only in where the crossing
lands inside a step (a grid offset), then refine.  A unique solution makes
those runs agree to discretisation error, shrinking as `h` shrinks.  A
non-unique one makes them separate by `O(1)` and STAY separated under
refinement.  The constant-rank control must show the shrinking, or the test
cannot tell the two apart.

THE ANSWER (2026-09-10).  **The hypothesis bites, and it costs the high-order
methods their entire advantage.**  Spread over grid placements, at 800 points
per period:

    method    LINEAR C     NONLINEAR C>0    rank C DROPS
    gear      3.00e-06     2.93e-06         1.91e-05     rate 3.99x -> 2.83x
    trbdf2    5.10e-09     7.55e-09         8.63e-06     rate 7.97x -> 4.86x
    glm3      2.26e-12     3.65e-12         1.46e-05     rate 16.3x -> 3.59x
    radau     4.44e-16     6.66e-16         9.75e-06     rate ~exact -> 2.64x

On both controls the spread shrinks at the method's own LOCAL order (h^(p+1):
gear 4x, TR-BDF2 8x, GLM3 16x per doubling) and the four methods separate by
ten orders of magnitude.  Where rank C drops they all collapse onto the SAME
magnitude and the same slow rate -- radau loses EIGHT ORDERS and ends up no
better than gear.

⚠ CONTROL B IS WHAT MAKES THAT MEAN ANYTHING.  A nonlinear `C(v) > 0` that
varies by 2x across the orbit behaves IDENTICALLY to a linear one -- 3.97x /
7.94x / 16.97x against 3.97x / 7.94x / 16.77x -- so the loss is the RANK
CHANGE and not the nonlinearity.  Without this fixture the comparison would
have been confounded and the conclusion unsupported.

⚠ AND `q = c0 V^3/3` IS A POLYNOMIAL, infinitely differentiable.  The order
loss is therefore not a smoothness failure of the model.

⚠⚠ WHAT THIS DOES NOT SHOW.  The spread still SHRINKS on the violated
fixture, so this is an ORDER COLLAPSE and NOT the uniqueness failure the
theory points at.  This orbit crosses the boundary TRANSVERSALLY at isolated
points; Thm 3.53's non-uniqueness concerns solutions AT or ALONG the border.
That case is not built here.
"""
import warnings

import numpy as np
import sympy

import pycircuit.circuit.circuit as _cc
from pycircuit.circuit.circuit import SubCircuit, gnd, defaultepar
from pycircuit.circuit.elements import R, VSin, C as Cap
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter, ddt)
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.integrator import (Gear2Integrator, TRBDF2Integrator,
                                          RadauIIA3Integrator,
                                          GLM3Integrator)

## `zero_order_sweep` needs these; `CubicCap` above already uses the rest
import sympy  # noqa: F811,E402

warnings.simplefilter('ignore')
_cc.default_toolkit = numeric
PER = 1e-3


class CubicCap(Behavioural):
    """`q = c0 V^3 / 3`, so `C = dq/dV = c0 V^2` VANISHES at `V = 0`.

    rank C is 1 away from zero and 0 at it, so an orbit crossing `V = 0`
    crosses a regularity boundary -- which is what the hypothesis forbids.
    """
    instparams = [Parameter(name='c0', desc='coefficient', unit='F/V^2',
                            default=1e-6)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return (Contribution(b.I, ddt(c0 * b.V ** 3 / 3)),)   # noqa: F821


def rank_changing(va=1.0):
    c = SubCircuit()
    c.add_node('s')
    c.add_node('a')
    c['vs'] = VSin('s', gnd, va=va, freq=1.0 / PER)
    c['rs'] = R('s', 'a', r=1e3)
    c['cq'] = CubicCap('a', gnd, c0=1e-6)
    c['rl'] = R('a', gnd, r=1e4)
    return c


def constant_rank(va=1.0):
    """The CONTROL: same topology, a LINEAR capacitor, so `C` is constant and
    `im D` is time-invariant.  Without this the branch test cannot tell
    non-uniqueness from ordinary discretisation error."""
    c = SubCircuit()
    c.add_node('s')
    c.add_node('a')
    c['vs'] = VSin('s', gnd, va=va, freq=1.0 / PER)
    c['rs'] = R('s', 'a', r=1e3)
    c['cq'] = Cap('a', gnd, c=1e-7)
    c['rl'] = R('a', gnd, r=1e4)
    return c


class PositiveNonlinearCap(Behavioural):
    """`q = c0 (V + V^3/3)`, so `C = c0 (1 + V^2)` -- STRICTLY POSITIVE, so
    rank C is 1 everywhere and `im D` is time-invariant DESPITE C varying by
    a factor of 2 across the orbit.

    ⚠⚠ THIS IS THE CONTROL THAT DECIDES WHAT THE RESULT MEANS.  Without it the
    comparison is confounded: the rank-changing fixture has a NONLINEAR `C`
    and the linear-capacitor control does not, so an order loss could be
    nonlinearity rather than the rank change.  This fixture is nonlinear and
    rank-constant, so it separates the two.  (It is also, per the source
    reading, exactly the shape that does NOT exercise the hypothesis -- a
    smoothly varying `C(v) > 0` is rank 1 throughout.)
    """
    instparams = [Parameter(name='c0', desc='coefficient', unit='F',
                            default=1e-7)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return (Contribution(b.I, ddt(c0 * (b.V + b.V ** 3 / 3))),)  # noqa: F821


def nonlinear_constant_rank(va=1.0):
    c = SubCircuit()
    c.add_node('s')
    c.add_node('a')
    c['vs'] = VSin('s', gnd, va=va, freq=1.0 / PER)
    c['rs'] = R('s', 'a', r=1e3)
    c['cq'] = PositiveNonlinearCap('a', gnd, c0=1e-7)
    c['rl'] = R('a', gnd, r=1e4)
    return c


def rank_probe(build):
    """The instrument check: does the fixture ACTUALLY change rank?  A fixture
    that does not exercise the condition is the whole trap this file is about,
    so this runs before anything is concluded."""
    cir = build()
    ia = cir.get_node_index('a')
    ranks = []
    for v in (-1.0, -0.1, -0.01, 0.0, 0.01, 0.1, 1.0):
        x = np.zeros(cir.n)
        x[ia] = v
        C = np.asarray(cir.C(x, defaultepar), dtype=float)
        tol = 1e-12 * max(abs(C).max(), 1.0)
        ranks.append((v, float(abs(C).max()),
                      int(np.linalg.matrix_rank(C, tol=tol))))
    return ranks


def endpoint(cls, build, npts, offset=0.0):
    """One transient from 0 to PER over a grid whose INTERIOR nodes are
    shifted by `offset` steps, so the zero crossing lands somewhere different
    inside a step -- while the run still ENDS AT EXACTLY `PER`.

    ⚠⚠ THE OBVIOUS VERSION OF THIS IS AN ARTEFACT, and it was this file's
    sixth instrument failure, caught by the same tell as the noise floor's
    four: shifting `tend` by `offset*h` makes the runs end at DIFFERENT TIMES,
    so the endpoint differs by `|dV/dt| * h * doffset` -- an `O(h)` spread that
    is METHOD-INDEPENDENT.  Measured, it read 1.7129e-02 / 1.7130e-02 /
    1.7130e-02 / 1.7130e-02 for gear / TR-BDF2 / radau / GLM3 and shrank at
    exactly 2.00x per doubling for all of them, order 2 and order 5 alike.
    Two different integrators agreeing to five significant figures on a
    quantity that is supposed to be their own error.
    """
    cir = build()
    h = PER / npts
    ## endpoints pinned at 0 and PER; interior nodes shifted by offset*h
    inner = [(k - offset) * h for k in range(1, npts + 1)]
    times = [0.0] + [t for t in inner if 1e-15 < t < PER * (1 - 1e-15)] + [PER]
    tr = Transient(cir, integrator=cls(), reltol=1e-11)
    ia = cir.get_node_index('a')
    x = np.asarray(DC(cir, refnode=gnd).solve().x, dtype=float).ravel()
    tr.irefnode = cir.get_node_index(gnd)
    tr.epar.t = 0.0
    tr._begin_run(x, cir.n)
    for j in range(1, len(times)):
        tr._dt_last = tr._dt if j > 1 else None
        tr._dt = times[j] - times[j - 1]
        tr.epar.t = times[j]
        x, _f, _J, _ = tr.solve_timestep(x, times[j])
        tr._push_history(x)
    return float(np.asarray(x, dtype=float)[ia])


def branch_test(build, label):
    """Do individually legitimate runs that differ only in grid PLACEMENT
    agree, and does their spread SHRINK under refinement?"""
    print('=== %s ===' % label)
    print('%-8s %6s %14s %16s %10s'
          % ('method', 'npts', 'spread over', 'endpoint', 'shrinks by'))
    print('%-8s %6s %14s %16s' % ('', '', 'grid offsets', 'mean'))
    for cls, name in ((Gear2Integrator, 'gear'), (TRBDF2Integrator, 'trbdf2'),
                      (RadauIIA3Integrator, 'radau'),
                      (GLM3Integrator, 'glm3')):
        prev = None
        for npts in (200, 400, 800):
            vals = [endpoint(cls, build, npts, off)
                    for off in (0.0, 0.17, 0.31, 0.53, 0.79)]
            v = np.array(vals)
            spread = float(v.max() - v.min())
            tag = ''
            if prev is not None:
                tag = ('%.2fx' % (prev / spread)) if spread > 0 else 'exact'
            print('%-8s %6d %14.4e %16.9e %10s'
                  % (name, npts, spread, float(v.mean()), tag))
            prev = spread
        print()


if __name__ == '__main__':
    for build, label in ((rank_changing, 'rank-CHANGING (CubicCap)'),
                         (nonlinear_constant_rank, 'NONLINEAR, rank-constant'),
                         (constant_rank, 'LINEAR control')):
        print('rank C along V(a) -- %s' % label)
        for v, mx, r in rank_probe(build):
            print('   V=%+6.2f  max|C| %.3e  rank %d' % (v, mx, r))
        print()
    branch_test(constant_rank,
                'CONTROL A: LINEAR C, im D time-invariant')
    branch_test(nonlinear_constant_rank,
                'CONTROL B: NONLINEAR C > 0, im D STILL time-invariant '
                '-- separates the rank change from the nonlinearity')
    branch_test(rank_changing, 'im D VIOLATED: rank C drops at V = 0')
    zero_order_sweep()


## ------------------------------------------------------------------------
## Does the collapse carry DAE content, or only the order of C's zero?

def _cap_of_order(k):
    """A capacitance with a zero of order `k` at `V = 0`: `C = c0 |V|^k`.
    `k = 0` is the linear control, `k = 2` is `CubicCap`."""
    class Cap(Behavioural):
        instparams = [Parameter(name='c0', desc='c', unit='F', default=1e-6)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            if k == 1:
                q = c0 * b.V * sympy.Abs(b.V) / 2                # noqa: F821
            elif k == 2:
                q = c0 * b.V ** 3 / 3                            # noqa: F821
            else:
                q = c0 * b.V * sympy.Abs(b.V) ** k / (k + 1)     # noqa: F821
            return (Contribution(b.I, ddt(q)),)

    def build(va=1.0):
        c = SubCircuit()
        c.add_node('s')
        c.add_node('a')
        c['vs'] = VSin('s', gnd, va=va, freq=1.0 / PER)
        c['rs'] = R('s', 'a', r=1e3)
        c['cq'] = Cap('a', gnd, c0=1e-6)
        c['rl'] = R('a', gnd, r=1e4)
        return c
    return build


def zero_order_sweep():
    """⚠⚠ THE MEASUREMENT THAT NARROWS THIS FILE'S OWN CONCLUSION.

    docs-46 reduced the collapse to a ONE-NODE SCALAR model with no DAE
    structure at all and found the exponent set by the ORDER OF THE ZERO of
    `C`, not by the method's order -- predicting `spread ~ h^p` with `p = 3/2`
    at `k = 1` and `4/3` at `k = 2`.  If that holds here, the collapse is a
    property of a vanishing capacitance and carries no index content, and
    `im D(t)` is the right DESCRIPTION of when it happens without being the
    MECHANISM.

    MEASURED (N = 200/400/800/1600): TR-BDF2 1.496 -> 1.413 and GLM3 1.527 ->
    1.411 as `k` goes 1 -> 2, both landing on the prediction at `k = 1` and
    moving the right way.  ⚠ radau reads 1.611 -> 1.709, noisier and moving
    the WRONG way with `k`, which is not explained and is recorded rather than
    fitted.
    """
    print('=== the collapse exponent against the order of C\'s zero ===')
    print('    predicted p: 1.500 at k=1, 1.333 at k=2, method-INDEPENDENT')
    print('%-8s %-4s %s   fitted p'
          % ('method', 'k', '  '.join('N=%-9d' % n
                                      for n in (200, 400, 800, 1600))))
    for cls, name in ((TRBDF2Integrator, 'trbdf2'),
                      (RadauIIA3Integrator, 'radau'),
                      (GLM3Integrator, 'glm3')):
        for k in (1, 2):
            build = _cap_of_order(k)
            sp = []
            for N in (200, 400, 800, 1600):
                v = [endpoint(cls, build, N, o)
                     for o in (0.0, 0.17, 0.31, 0.53, 0.79)]
                sp.append(max(v) - min(v))
            ns = np.array([200.0, 400.0, 800.0, 1600.0])
            p = -np.polyfit(np.log(ns), np.log(np.array(sp)), 1)[0]
            print('%-8s %-4d %s   %.3f'
                  % (name, k, '  '.join('%.3e' % x for x in sp), p))
    print()
