"""Is the `1/h` amplification of a defect LOCAL, or does it PROPAGATE?

THE QUESTION.  2026-09-10: the filtered local-error estimate reads one order
BELOW its declared `EMBEDDED_ORDER + 1` on a DAE, for every method, because
`J = C + a h G` is singular in `C` so `J^-1` behaves like `1/h` on the
algebraic subspace.  Whether that matters depends entirely on whether the
amplified part is LOCAL.

THE PUBLISHED FORM (relayed from docs-46, CITED not re-derived here).
Proposition 8.10 for BDF on a regular index-2 DAE with properly stated leading
term carries exactly that `1/h` term, and Prop 8.6 generalises it to
`h^-(mu-1)` at index `mu`.  Section 8.4 note (6): "In the case of higher index
constant coefficient DAEs the errors (1/h^i)delta_l are LOCAL; they are not
propagated. This situation is also given in the case of linear variable
coefficient index-2 DAEs with properly stated leading term. Unfortunately,
already in the case of linear index-3 DAEs those bad error terms can be
propagated."

⚠⚠ AN EARLIER ATTEMPT AT THIS WAS DELETED RATHER THAN COMMITTED.  It measured
global error against INTERVAL LENGTH at fixed `h` and read `alpha = 0.00` for
every method on both an index-1 and an index-2 fixture -- not because there is
no propagation, but because these are CONTRACTING driven systems whose error
saturates regardless, so "flat in interval length" was guaranteed before it
started.  A gate that cannot fail is not a gate.  The interval sweep was the
right instinct aimed at the wrong quantity: note (6) does not claim the error
is SMALL, it claims the AMPLIFIED PART IS LOCAL -- one defect's fate over the
next few steps, not accumulated error over an interval.

⚠ AND THE OBVIOUS FIX IS ALSO WRONG.  Asserting `err(n0)/err(n0+2)` hard-codes
an offset of 2, but a defect does not die at a fixed offset -- it dies after
MU steps.  An index-3 defect that vanishes at `n0+3` would be scored as
"propagated" by an offset-2 rule.  The discriminator is "is it EXACTLY ZERO
after `mu` steps", with `mu` supplied by the fixture.

THE DESIGN.  Inject a defect of known size at ONE step and difference against
the undisturbed run of the SAME discretisation, so truncation cancels exactly
rather than swamping the defect at the `h` where the question gets interesting.
Then:

  * INSTRUMENT ALIVE: `|e|` at the injection step must scale as `h^-(mu-1)` --
    it must GROW as the grid refines.  Without that the amplification is not
    being measured at all and anything downstream is meaningless.
  * LOCAL: `|e|` is exactly zero from step `n0 + mu` on.
  * PROPAGATED: it survives `mu` steps and still scales as `h^-(mu-1)`.

⚠⚠⚠ THIS MEASUREMENT DID NOT COME OFF, AND THE FILE IS KEPT FOR THAT REASON.
ITS OWN INSTRUMENT-ALIVE CHECK FAILS: `|e|` at the injection step reads
EXACTLY the injected `delta` at every grid, amplification exponent +0.00 where
the theorem wants -1 at index 2.  The `h^-(mu-1)` amplification is not
reproduced here, so the "exactly zero after mu steps" readings below prove
NOTHING and must not be quoted as an answer.  What is recorded is which
injections are already known to be wrong:

  * INTO THE STATE (add `delta` to `x` after the step converges).  Then `|e|`
    at that step is trivially `|delta|`, at every `h`.  It measures the kick.
  * INTO THE SOURCE VECTOR via `provided_function`, which is added to `u`.
    Closer -- it does perturb the equation -- but still reads `|e| = delta`
    exactly, and on the index-2 C-V loop the tail is ~100x SMALLER than
    `delta` and decays slowly rather than being amplified and then vanishing.

So the next attempt should start by asking WHERE the theorem's `q_ni` enters
this assembly, rather than by choosing a plausible-looking hook.  A defect
that enters as a current on an algebraic row is not obviously the same object
as a perturbation of the stage equations in the theorem's coordinates, and
that identification is the step that was skipped here.

The QUESTION remains open: the arm of note (6) that covers a
variable-coefficient nonlinear MNA is the index-2 one, it is stated to hold,
and the margin to the propagating index-3 case is one index level.

Run: `python benchmarks/defect_locality.py`
"""
import warnings

import numpy as np

from pycircuit.circuit.circuit import gnd
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.integrator import (Gear2Integrator, TRBDF2Integrator,
                                          RadauIIA3Integrator, GLM3Integrator)

warnings.simplefilter('ignore')
PER = 1e-3


def index2():
    """The C-V loop: index 2 (a capacitor loop closed by a voltage source)."""
    from pycircuit.circuit.tests.test_glm import _cv_loop
    return _cv_loop(PER)


def index1():
    """The state-free exponential: index 1."""
    from pycircuit.circuit.tests.test_stage_predictor import _expg_fixture
    return _expg_fixture(PER)


def march(cls, build, h, nsteps, kick_at=None, kick=None):
    """`nsteps` steps of `h`.  A defect, if given, perturbs the STEP EQUATION
    at step `kick_at` -- not the state.

    ⚠⚠ INJECTING INTO THE STATE MEASURES NOTHING.  Adding `delta` to `x` after
    a step converges makes `|e|` at that step exactly `|delta|`, at every `h` --
    measured, amplification exponent +0.00 where the theorem wants -1.  The
    quantity Prop 8.10 amplifies is `q_ni`, a perturbation of the EQUATION, and
    its effect `delta/h` appears in the SOLUTION.  `provided_function` is
    exactly that hook: it is added to the source vector `u`, so a constant
    vector there is a defect in the residual.
    """
    cir = build()
    tr = Transient(cls and cir, integrator=cls(), reltol=1e-13)
    tr.irefnode = cir.get_node_index(gnd)
    x = np.asarray(DC(cir, refnode=gnd).solve().x, dtype=float).ravel()
    tr.epar.t = 0.0
    tr._begin_run(x, cir.n)
    xs = [x.copy()]
    zero = np.zeros(cir.n)
    for j in range(1, nsteps + 1):
        tr._dt_last = tr._dt if j > 1 else None
        tr._dt = h
        tr.epar.t = j * h
        if kick_at is not None and j == kick_at:
            pf = (lambda _t, _k=np.asarray(kick, dtype=float): _k)
        else:
            pf = (lambda _t, _z=zero: _z)
        x, _f, _J, _ = tr.solve_timestep(x, j * h, provided_function=pf)
        tr._push_history(x)
        xs.append(np.asarray(x, dtype=float).copy())
    return np.array(xs)


def subspaces(build, x):
    """`(differential, algebraic)` unit directions from an SVD of `C`.

    `im C` is the differential subspace -- the one whose information is
    carried from step to step -- and its orthogonal complement is the
    algebraic one, where `J^-1`'s `1/h` lives.
    """
    cir = build()
    from pycircuit.circuit.circuit import defaultepar
    C = np.asarray(cir.C(x, defaultepar), dtype=float)
    U, sv, _Vt = np.linalg.svd(C)
    tol = max(sv) * 1e-10 if len(sv) else 0.0
    r = int((sv > tol).sum())
    return U[:, 0], (U[:, -1] if r < C.shape[0] else None), r


def locality(cls, name, build, label, mu, npts=(200, 400, 800), delta=1e-6):
    """Inject on the ALGEBRAIC direction and follow the difference."""
    print('  %-8s %-9s mu=%d' % (name, label, mu))
    print('    %6s %12s %s' % ('N', '|e| at n0', 'then, step by step'))
    amps = []
    for N in npts:
        h = PER / N
        n = min(N, 120)
        n0 = n // 2
        base = march(cls, build, h, n)
        d_dir, a_dir, _r = subspaces(build, base[n0])
        direction = a_dir if a_dir is not None else d_dir
        pert = march(cls, build, h, n, kick_at=n0, kick=delta * direction)
        e = np.max(np.abs(pert - base), axis=1)
        amps.append(float(e[n0]))
        tail = '  '.join('%.3e' % v for v in e[n0:n0 + 5])
        print('    %6d %12.4e  %s' % (N, e[n0], tail))
    sl = [np.log(amps[i] / amps[i - 1]) / np.log(2.0)
          for i in range(1, len(amps))]
    print('    amplification exponent in h: %s   (want %+d for mu=%d)'
          % (' '.join('%+.2f' % x for x in sl), -(mu - 1), mu))
    print()


def main():
    print('A defect injected at ONE step, differenced against the undisturbed')
    print('run of the SAME discretisation.  LOCAL means exactly zero from')
    print('step n0+mu on; the amplification at n0 must scale as h^-(mu-1),')
    print('or the instrument is not measuring the amplification at all.')
    print()
    for cls, name in ((Gear2Integrator, 'gear'), (TRBDF2Integrator, 'trbdf2'),
                      (RadauIIA3Integrator, 'radau'),
                      (GLM3Integrator, 'glm3')):
        locality(cls, name, index1, 'index-1', 1)
        locality(cls, name, index2, 'index-2', 2)


if __name__ == '__main__':
    main()


def jinv_column_exponents(build, label, known_index, hs=(1e-6, 1e-7, 1e-8)):
    """docs-46's route: the amplifying directions are the COLUMNS OF `J^-1`
    that grow as `h` shrinks, and `max exponent + 1` is the tractability index.

    At the injection step the previous error is zero, so `e = J^-1(-delta)`
    with `J = A D/h + B` -- the step matrix already assembled and factored.
    Probing unit vectors therefore needs no `G_2`, no `Q_1` and no decoupling.
    Verified by them to every digit on the book's index-1/2/3 fixtures, with
    the caveat that they could not run it on a real MNA matrix.

    ⚠⚠ MEASURED HERE, AND IT DOES NOT TRANSFER.  On charge-oriented MNA:

        index-1 (ExpG)      a_i = [+0.00, -0.38, +0.00]   max+1 = 1   MATCH
        index-2 (C-V loop)  a_i = [+0.00, -1.00, +0.00]   max+1 = 1   NO

    `topological_index` independently calls these 1 and 2, and names the loop
    (`vs`, `c1`, `c2`).  No column of `J^-1` grows on the index-2 circuit --
    one SHRINKS exactly as `h` (the differential direction) and the rest are
    flat -- so the recipe reads index 1 where the circuit is index 2.

    ⚠ AND `J` STAYS WELL-CONDITIONED, which is the sharper form of the same
    fact: cond(J) reads 4.98e+02 / 5.02e+01 / 5.34e+00 / 3.73e+00 at
    h = 1e-6 / 1e-7 / 1e-8 / 1e-9.  It DECREASES as the grid refines, where an
    index-2 DAE is expected to give `1/h^2`.  So the index-2 character of this
    circuit is not visible in `C/h + G` at all, and no probing of that matrix
    can recover it.

    The identification gap is therefore still open and is the same one: `C` in
    this assembly is the product `A D`, and the theorem's per-row structure
    lives in coordinates where `A` and `D` are separate.  "Needs nothing from
    the theorem" holds in the book's coordinates and not in these.
    """
    from pycircuit.circuit.circuit import defaultepar
    from pycircuit.circuit.analysis import remove_row_col
    from pycircuit.circuit.dcanalysis import DC as _DC
    cir = build()
    tr = Transient(cir, integrator=Gear2Integrator(), reltol=1e-12)
    iref = cir.get_node_index(gnd)
    x = np.asarray(_DC(cir, refnode=gnd).solve().x, dtype=float).ravel()
    norms = []
    for h in hs:
        C = np.asarray(cir.C(x, defaultepar), dtype=float)
        G = np.asarray(cir.G(x, defaultepar), dtype=float)
        (Jr,) = remove_row_col((C / h + G,), iref, tr.toolkit)
        Ji = np.linalg.inv(np.asarray(Jr, dtype=float))
        norms.append(np.abs(Ji).max(axis=0))
        print('    h=%.0e  |J^-1| cols %s   cond %.2e'
              % (h, '  '.join('%.4e' % v for v in norms[-1]),
                 np.linalg.cond(np.asarray(Jr, dtype=float))))
    a = [np.log(norms[-1][k] / norms[0][k]) / np.log(hs[0] / hs[-1])
         for k in range(len(norms[0]))]
    print('  %-22s a_i = [%s]   max+1 = %.0f   index = %d   %s'
          % (label, ', '.join('%+.2f' % v for v in a), max(a) + 1,
             known_index,
             'MATCH' if abs(max(a) + 1 - known_index) < 0.25 else 'NO'))
    return a
