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

    ⚠⚠⚠ THE CONCLUSION FIRST RECORDED HERE WAS WRONG, AND IT WAS WRONG IN THE
    WAY THIS SESSION KEEPS BEING WRONG: measured in the PRE-ASYMPTOTIC WINDOW.
    It read

        index-1 (ExpG)      a_i = [+0.00, -0.38, +0.00]   max+1 = 1   MATCH
        index-2 (C-V loop)  a_i = [+0.00, -1.00, +0.00]   max+1 = 1   NO

    and concluded the route "does not transfer to MNA".  It does.  See
    :func:`sigma_min_index`.  What was wrong was the STATISTIC and the RANGE:

    * `cond(J)` mixes `sigma_max ~ C/h`, which grows trivially for any circuit
      containing a capacitor and has nothing to do with the index, with
      `sigma_min`, where the index actually lives.  Use `sigma_min` or
      equivalently `||J^-1||_2` (docs-46's correction, confirmed here).
    * the per-COLUMN probe is BASIS-DEPENDENT and the which-row information
      does not survive a change of coordinates -- withdrawn by its author
      after testing it under random invertible transforms.
    * and the sweep 1e-6..1e-9 sat AT THE TURNING POINT, where the reading is
      neither value.

    ⚠ `cond(J)` falling in that window is real but is not evidence of
    anything: its author could not reproduce the FALL with an honest row
    scaling and withdrew the near-match as coming from a broken fixture.
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


def sigma_min_index(build, label, known_index, decades=22, per_decade=1):
    """The index from `sigma_min(J)`, `J = C/h + G`, with the window LOCATED
    rather than assumed.

    `||J^-1||_2 ~ h^-(mu-1)`, so the exponent plus one is the tractability
    index.  Two factorisations and an SVD; no projector sequence, no symbolic
    structure.  ⚠ NOT per row -- that part of the original proposal is
    basis-dependent and was withdrawn by its author after testing it under
    random invertible transforms.

    ⚠⚠ THE WINDOW IS THE WHOLE DIFFICULTY, AND A FIXED RANGE IS THE TRAP.
    Above the window the reactive term `C/h` is negligible against `G`, `J` is
    effectively resistive, and the index character is ABSENT -- not weak,
    absent -- so the probe correctly reads 1 for everything.  MEASURED on the
    index-2 C-V loop, exponent of `||J^-1||` per decade:

        1e-3 .. 1e-9    -0.62  -0.99  -0.78     reads index 1   WRONG
        1e-9 .. 1e-15   +0.78  +1.00  +1.00     reads index 2   right

    An earlier version of this function hard-coded `1e-9 .. 1e-17`, which is
    the same mistake one layer up: the turn tracks `||C||/||G||` FROM THE
    ASSEMBLED MATRICES (docs-46, measured over 15 decades with the two scales
    varied independently -- no appeal to `RC` or to any physical time
    constant), but THE CONSTANT IS FIXTURE-DEPENDENT AND LANDS ON EITHER SIDE:
    theirs sits 7.6x ABOVE `||C||/||G||`, this tree's C-V loop about 0.05x
    BELOW it.  So the window is SWEPT and the turn LOCATED, never computed and
    probed at.

    ⚠ A single decade cannot distinguish "flat because index 1" from "flat
    because pre-asymptotic".  That is exactly the trap both sessions fell into,
    and it is why the sweep here is wide by default.

    ⚠ There is a roundoff floor at `eps * ||C||/||G||`, below which `C/h`
    swamps `G` in the sum and `J` goes numerically singular though it is
    mathematically fine.  For these fixtures that is ~4e-25, so the usable
    window is about 16 decades wide and CENTRED on the crossover, not
    unbounded.  Rows below the floor are dropped.

    ⚠⚠⚠ DO NOT CLAMP THE EXPONENT AT ZERO.  A previous version did, to make a
    van der Pol's `a = -1.000` agree with `topological_index`'s 1 -- and that
    was AGREEING FOR THE WRONG REASON.  `topological_index` implements
    Estevez Schwarz & Tischendorf, which is FLOORED AT 1 BY CONSTRUCTION
    ("index 2 iff the network contains a C-V loop or an L-I cutset, otherwise
    1"); its own docstring records that Theorem 3.47's INDEX-0 case is not
    implemented.  So the probe was being clamped to match a reference that
    cannot represent the answer the probe was giving.

    MEASURED directly, which settles it without the topological criterion at
    all -- INDEX 0 MEANS `C` IS NONSINGULAR, an implicit ODE:

        fixture       rank C (reduced)   sigma_min(C)   a        index
        van der Pol        2 of 2          1.0000e+00   -1.000     0
        ExpG               1 of 3          0.0000e+00   +0.000     1
        C-V loop           2 of 3          0.0000e+00   +1.000     2

    The van der Pol has C, L and a BSource and NO voltage source, so it also
    meets Thm 3.47's stated condition.  The rule is UNCLAMPED,
    `index = max_i a_i + 1`, and it extends DOWN as well as up: -1 -> 0,
    0 -> 1, +1 -> 2, +2 -> 3.  The probe is strictly more capable than the
    topological criterion here, and clamping threw that away.
    """
    from pycircuit.circuit.circuit import defaultepar
    from pycircuit.circuit.analysis import remove_row_col
    from pycircuit.circuit.dcanalysis import DC as _DC
    cir = build()
    tr = Transient(cir, integrator=Gear2Integrator(), reltol=1e-12)
    iref = cir.get_node_index(gnd)
    x = np.asarray(_DC(cir, refnode=gnd).solve().x, dtype=float).ravel()
    C = np.asarray(cir.C(x, defaultepar), dtype=float)
    G = np.asarray(cir.G(x, defaultepar), dtype=float)
    ratio = np.abs(C).max() / max(np.abs(G).max(), 1e-300)
    floor = np.finfo(float).eps * ratio
    ## centred on the crossover, wide on both sides, stopping at the floor
    top = np.log10(ratio) + 3.0
    hs = [10.0 ** e for e in np.arange(top, top - decades, -1.0 / per_decade)
          if 10.0 ** e > 1e2 * floor]
    smin = []
    for h in hs:
        (Jr,) = remove_row_col((C / h + G,), iref, tr.toolkit)
        smin.append(float(np.linalg.svd(np.asarray(Jr, dtype=float),
                                        compute_uv=False).min()))
    a = [np.log(smin[i - 1] / smin[i]) / np.log(hs[i - 1] / hs[i])
         for i in range(1, len(hs))]
    asym = a[-3:]
    ## ⚠ UNCLAMPED: -1 means index 0 (C nonsingular, an implicit ODE), and
    ## clamping it to 1 only agreed with a reference that is floored at 1
    idx_meas = round(max(asym)) + 1
    turn = None
    for i in range(len(a) - 1, 0, -1):
        if a[i - 1] < 0.5 <= a[i]:
            turn = hs[i]
            break
    print('  %-22s ||C||/||G|| = %.2e  floor = %.1e  turn = %s'
          % (label, ratio, floor,
             ('%.1e (%.3gx the ratio)' % (turn, turn / ratio)) if turn
             else 'none found (flat throughout)'))
    print('  %-22s asymptotic %s   max(a)+1 = %.0f   index = %d   %s'
          % ('', ' '.join('%+.3f' % v for v in asym), idx_meas, known_index,
             'MATCH' if abs(idx_meas - known_index) < 0.5 else 'NO'))
    return asym
