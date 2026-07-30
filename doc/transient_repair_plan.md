# Repairing the transient analysis — the plan

**Status: written 2026-07-30; stages 1-5 executed 2026-07-30. Every OUTCOME line
below was filled in from a real run. Stage 6 remains blocked and unattempted.**

**Same-box baseline for every runtime and tally figure below**, re-measured at
`b2ab5bb` immediately before stage 1 rather than taken from the header:
**715 passed, 6 skipped, 0 failed in 578.49 s**, `test_stress_stiff_rlc_pulse`
91.35 s. This box is slower than the one the 492 s figure came from (+17.6%), so
runtime gates are judged against 578.49 s; both numbers are quoted where it
matters.

Reasoning, scope and rejections: `transient_repair_reasoning.md`. That document is
where an argument about *whether* to do this belongs. This one is the order of work and
the gates.

Baselines to beat, both measured at `b2ab5bb`:

- full suite `-m "" --timeout=400`: **715 passed, 6 skipped, 0 failed**, ~492 s clean
  (617 s when other jobs share the box)
- doc build: **build succeeded, 2 warnings, 0 ERROR lines**

Running the suite: `pytest.ini` deselects 17 `slow` tests by default. **`-m ""` is
mandatory before reporting any result here** — the deselected ones are the only tests
comparing an analysis against a time-domain reference, which is exactly the evidence
class this work is about.

---

## Stage 1 — fix `Gear2Integrator.compute_lte`

Replace the `'classic'` branch (`integrator.py:192-206`) so it estimates `q'''` scaled
by `h**2`: second divided difference of `g = dq/dt` taken from the `iq` companion-current
history, with the VSS coefficients
`alpha0 = (2h1+h2)/(h1(h1+h2))`, `alpha1 = -(h1+h2)/(h1 h2)`, `alpha2 = h1/(h2(h1+h2))`
and `lte = -(1/3) h1 (h1+h2) * dd2_g`. Keep the returned `p = 3.0` (it is order+1 and is
already correct).

Do **not** touch the default `lte_formula`, and do **not** touch `stepcontroller.py`, in
this stage. One change, one measurement.

**Gate 1-1 (the estimate is asymptotically exact).** Unit-level, no circuit: on an
analytic `q(t)` with uniform `h` and full supplied history, `estimate/true` must
approach 1 as `h` halves, and the observed order in `h` must be 2.0 +/- 0.1 (currently
2.99). Declared success: ratio within 2% of 1 at the smallest `h`, monotonically
approaching it.
OUTCOME: **PASSED.** `q(t) = [sin(2pi 1e6 t), 0.5 cos(0.7*2pi 1e6 t)]`, uniform
`h` swept 4 ns -> 0.25 ns, history built by running the integrator's own
companion recursion forward from an exact start (which is what a run stores).
Reference = `iq - q'(t_n)`, the residual the controller maps through `J^-1`.

- before: ratio 4.7e-16 -> 1.2e-16, i.e. flat zero; order **3.001**
- after: ratio 0.9795 -> 0.9901 -> 0.9951 -> 0.9976 -> **0.9988** (monotone,
  0.12% off 1); order **1.993**

Independent cross-check, against the *textbook one-step* LTE instead (exact past
derivatives supplied): the patched branch gives **0.66545** where the algebra
predicts exactly 2/3 -- and 0.00000 before. Euler gives 0.50008 (predicted 1/2),
Gear2-YWR 0.49909 (predicted 1/2), Trapezoidal 0.83091 (predicted 5/6). Four
independently derived constants all hit, so the derivation is not being fitted.

**Recorded negative result, found while building this gate:** the
companion-consistent reference is *ill-posed for Trapezoidal*. Trap is the only
method whose companion current depends on its own past value (`iq = 2 dq/h -
iq_last`), a recursion with eigenvalue -1, so its companion error carries an
undamped alternating component. Starting the warm-up from an exact derivative
puts that component at a node, and the reference collapses to O(h^3) noise
(ratio 22.0 at even warm-up lengths, 0.75 at odd ones). This is why the stage-4
unit test pins ratios against the one-step reference, which is well-conditioned
for all three methods, and not against this one.

**Gate 1-2 (the controller is alive).** On a small stiff circuit, Gear2 `'classic'` must
reject at least one step, and its step count must change by >20% when `reltol` moves
1e-4 -> 1e-6. Declared success: both true. Currently: zero rejections ever, bit-identical
across tolerances.
OUTCOME: **PASSED.** Reference case: series RLC loop `C(1,gnd)-R(1,2)-L(2,gnd)`,
R=1k L=1uH C=1uF, poles -1.000e+03 and -1.000e+09 (stiffness 1e6), released from
v1 = 1 V, `tend` 5e-3 (5 slow tau), `max_step` 2e-4.

| Gear2 `'classic'` | steps @1e-4 | rejections @1e-4 | steps @1e-6 |
|---|---|---|---|
| before | 26 | **0** | 26 (+0.0%) |
| after | 40 | **4** | 171 (**+327.5%**) |

Both conditions met. For scale, `'ywr'` on the same case is 37 steps / 4
rejections / +321.6% -- the patched `'classic'` now behaves like the formula that
was already working, which is the outcome the fix predicts.

**Gate 1-3 (accuracy, against an independent reference).** On a linear stiff case with
an analytic matrix-exponential solution, patched Gear2's max error must be within 2x of
Trapezoidal's at the same tolerance. Declared success: within 2x. This is stated against
a *small* circuit deliberately — the leapfrog is not a valid reference here, see the
reasoning document.
OUTCOME: **PASSED as declared, with a negative sub-result that must be read with
it.** Same reference case; reference solution `y(t) = expm(A t) y0` from
`scipy.linalg.expm`, not another integration.

Declared metric -- max error over the run, Gear2 vs Trapezoidal at reltol 1e-4:

| `max_step` | gear2-classic | trap-classic | ratio |
|---|---|---|---|
| 2e-4 | 1.4603e-02 | 1.4603e-02 | **1.000** |
| 1e-5 | 4.9544e-05 | 4.9176e-05 | **1.007** |

Within 2x, so the gate passes. **But the metric is degenerate, and saying why is
the more useful result:** at *both* `max_step` values the max error equals
`0.5 h^2 q''` for the genuine first step (5e-5 at h=1e-5 to two digits) and is
identical across methods to three digits, because every method drops to Backward
Euler for a step that is then accepted with no error check at all. Shrinking
`max_step` does not remove it -- it is structural, and it is defect (D) at the one
place stage 3 deliberately leaves it, the genuine first step.

Against a start-up-free metric (exact solution re-propagated from the simulated
state at the third time point) the ratio is **3.94** at `max_step` 1e-5 and
**3.86** at 2e-4 -- i.e. a stricter reading of this gate *fails*. Three facts say
that shortfall is not attributable to stage 1:

- `'ywr'`, untouched by stage 1, scores **3.940** against classic's 3.936 -- so it
  is a BDF-2-versus-Trapezoidal property, not a property of the estimate;
- stage 1 *improved* the number, 6.64 -> 3.86 at `max_step` 2e-4;
- Gear-2's error constant is 2/9 against Trapezoidal's 1/12, a 2.7x handicap
  before any step control enters.

Recorded rather than resolved: the plan's 2x threshold was written against the
degenerate metric and is not the right threshold for the sharp one. Also seen and
not investigated: `trap` + `'ywr'` at `max_step` 1e-5 took 5348 rejections for
1339 steps -- a rejection storm in a path nothing here touches.

**Gate 1-4 (blast radius, the one that decides whether this ships).** Full suite `-m ""`.
Declared success: 715 passed with **exactly one** expected failure,
`test_lte_formula_ywr`, which stage 2 then rewrites. Any *other* failure is a stop —
report it, do not fix it by loosening the test that caught it.
OUTCOME: **PASSED.** `1 failed, 715 passed, 6 skipped` — the one failure is
`test_analysis_transient.py::test_lte_formula_ywr`, exactly as predicted, and it
is the *only* one. Nothing else in the suite moved.

**Confound, recorded because it changes how every tally below must be read:** a
concurrent agent working the `leapfrog_5th_order` instability (the investigation
stage 6 is blocked on) committed `9cfa357`, `3fe5468` and `ff5c6e6` onto this
branch *while these stages were running*, moving HEAD off `b2ab5bb` and adding
one test (`test_benchmark_circuits.py::test_leapfrog_has_no_right_half_plane_poles`).
That is why 715 passed + 1 failed is the correct pass count here rather than 714:
the collected total went 721 -> 722. None of the files this repair touches
(`integrator.py`, `stepcontroller.py`, `transient.py`, `jaxtransient.py`,
`lte_dae.rst`, `test_analysis_transient.py`) were modified by those commits —
verified with `git diff b2ab5bb HEAD --` on that path list — so the attribution of
every measurement above and below is intact. From stage 2 on, the "all green"
tally is **716 passed, 6 skipped, 0 failed**.

**Gate 1-5 (suite runtime).** Record total suite time and
`test_stress_stiff_rlc_pulse` alone. Declared success: total within 20% of the 492 s
baseline. **Maintainer's decision, 2026-07-30: if it regresses beyond 20%, mark the
worst offender `@pytest.mark.slow`** rather than tuning tolerances — a fast test that no
longer stresses the controller is worse than a slow one that does. Record which test was
marked and both runtimes (default selection and `-m ""`), since this changes suite
behaviour for everyone.
OUTCOME: **PASSED, and it went the other way — the suite got faster.**

| | baseline (`b2ab5bb`) | after stage 1 | change |
|---|---|---|---|
| full suite `-m ""` | 578.49 s | **541.45 s** | **-6.4%** |
| `test_stress_stiff_rlc_pulse` | 91.35 s | 96.02 s | +5.1% |

Well inside the 20% band, so **no test was marked `@pytest.mark.slow`** and the
maintainer's fallback decision was not exercised. Worth recording anyway, because
it would not have worked as written: the heaviest transient tests, including
`test_stress_stiff_rlc_pulse`, live in `test_analysis_transient_stress.py`, which
already carries a module-level `pytestmark = pytest.mark.slow`. Marking "the worst
offender slow" would have been a no-op for exactly the tests it was aimed at, and
would have had no effect on an `-m ""` run at all. If runtime ever does become the
binding constraint, the lever has to be something else.

Why a correct estimate is *cheaper* here rather than dearer: the repair mostly
changes where the steps go, not how many. `test_stress_stiff_rlc_pulse` runs with
`timestep=1e-9`, so its ~25000 steps were already `max_step`-bound and adding
rejections costs it only ~5%; elsewhere the controller now grows the step where
the old estimate had it thrashing. The reasoning document's worry that this test
would blow up is a **refuted premise**, and P15 in `architecture.md` is also now
stale: it records this test at ~266 s of a ~492 s total, where it currently sits at
96 s of 541 s (the `e37ddad` tolerance-flavour fix, not this one, did that).

---

## Stage 2 — rewrite the test that asserts the bug

`test_analysis_transient.py::test_lte_formula_ywr` asserts `n_gy > n_gc` and
`e_gy < e_gc`, which hold only because `'classic'` is blind. Replace with assertions
true of a correct implementation: both formulas reject steps, both respond to `reltol`,
and both agree with an accurate reference to within a stated factor. Keep a test that
distinguishes the two formulas — `'ywr'` has a constant 0.75 bias, a correct `'classic'`
does not — so the file still has a reason to exist.

**Gate 2-1.** Full suite `-m ""` returns to **715 passed, 6 skipped, 0 failed**.
OUTCOME: **PASSED — 716 passed, 6 skipped, 0 failed in 481.08 s** (716, not 715,
for the concurrent-agent reason recorded under gate 1-4). Runtime now **-16.8%**
against the 578.49 s baseline.

Recorded because it nearly wasted an hour: the *first* attempt at this gate came
back `3 failed, 713 passed`, all three in `test_benchmark_circuits.py`, all with
`TypeError: build_leapfrog_network() got multiple values for argument 'stages'`
from `benchmark_circuits.py:857` — a file this repair does not touch, caught
mid-write by the concurrent agent (its mtime moved twice inside the collection
window). Re-running that module alone with nothing of ours changed gave
`27 passed`. The lesson is the one already in the working method: read the
failure, do not read the count. A tally that includes somebody else's
half-saved file looks exactly like a regression.

**Gate 2-2.** The rewritten test must FAIL if stage 1's change is reverted. A test that
passes against both implementations is not testing the thing that was broken; verify by
actually reverting, running it, and restoring.
OUTCOME: **PASSED.** Verified by actually reverting `integrator.py` with
`git checkout`, running the single test, and restoring from a saved copy:

```
E   AssertionError: gear2-classic rejected no step at all on a 10 us RC charge
E   assert 0 >= 1
pycircuit/circuit/tests/test_analysis_transient.py:468: AssertionError
1 failed in 0.46s
```

So the first assertion to bite is the one that names the actual defect. The test
also carries three assertions that pin the *constants* (`2/3`, `1/2`, and the
`4/3` between them, against the one-step LTE); those too are unreachable with the
old branch, whose ratio is ~1e-15.

---

## Stage 2b — switch the default `lte_formula` to `'ywr'`

**Maintainer's decision, 2026-07-30**: belt-and-braces. Stage 1 makes `'classic'`
mathematically correct; this ships `'ywr'` as the default anyway, on the grounds that it
has the longer track record. Accepted cost: `'ywr'` carries a constant 0.75 underestimate
of the true LTE where a corrected `'classic'` has none, so the default is slightly
optimistic about its own error. `TRTOL = 7.0` already absorbs more than that factor.

A **separate stage from 1 on purpose.** Stage 1 changes what the estimate computes; this
changes which estimate runs by default. Bundled, a step-count or accuracy change could
not be attributed to either. `Gear2Integrator.__init__`, `integrator.py:135`.

**Gate 2b-1 (the switch is what moved).** Record step count and accuracy on stage 1-3's
reference case for all four combinations: `'classic'` and `'ywr'`, before and after this
change. Declared success: post-change default behaviour is *identical* to explicitly
passing `lte_formula='ywr'`, and the corrected `'classic'` remains reachable and correct.
OUTCOME:

**Gate 2b-2.** Full suite `-m ""` at 715 passed, 6 skipped, 0 failed. Anything that moves
here is a test that was implicitly depending on the default, which is worth naming.
OUTCOME:

---

## Stage 3 — the unconditional post-breakpoint accept

`stepcontroller.py:26-28` exempts `is_first_step` from any error check;
`transient.py:314` re-arms that flag at every breakpoint; `Sin.next_event` fires every
quarter period. Distinguish "no history exists" (genuine first step of the run —
accepting is correct) from "history was just reset at a breakpoint" (a step that can and
should be bounded). Also remove the dead `err = 0.5`.

Separate from stages 1-2 on purpose: the 40% leapfrog waveform discrepancy is probably
Gear2 *and* this, and doing them together makes the attribution unrecoverable.

**Gate 3-1 (the hole is closed).** On a `VSin`-driven stiff circuit, count steps
accepted without an LTE evaluation. Declared success: exactly 1 for the whole run
(the genuine first step), against one per quarter period now.
OUTCOME:

**Gate 3-2 (discontinuities still work).** A `VPulse` edge must still be integrated
without a rejection loop — the breakpoint is there because the discontinuity is real,
and a small step is genuinely wanted. Declared success: `test_breakpoints.py` passes and
the pulse edge takes no more than 3 rejections.
OUTCOME:

**Gate 3-3.** Full suite `-m ""` at 715/6/0, runtime within 20% of baseline.
OUTCOME:

---

## Stage 4 — the tests that would have caught all of this

Three classes, none of which exist (see reasoning (E)):

1. `compute_lte` scales with the correct power of `h`, for **every** integrator. Unit
   test, analytic `q(t)`, no circuit. This is the cheap one and it is the one that would
   have caught the original defect.
2. A step is actually rejected on a stiff case — for every integrator. The current
   suite's `adapt_steps < fixed_steps` assertion is satisfied *most* emphatically by a
   blind controller, so it cannot serve.
3. Step count is monotone non-decreasing, and error monotone non-increasing, as `reltol`
   tightens — for every integrator.

**Gate 4-1.** Each new test must fail against the pre-stage-1 code and pass after.
Demonstrate both directions; a test only known to pass proves nothing about its power.
OUTCOME:

**Gate 4-2.** Full suite `-m ""`, new tally recorded (it will exceed 715).
OUTCOME:

---

## Stage 5 — the JAX copy, and the docs

`jaxtransient.py:247-254` duplicates the defective expression for `method='gear'`; lines
257-260 have the same shape for Trapezoidal. Either fix both to match stage 1 or delete
the dead branches — but do not leave two mutually-confirming wrong implementations,
which is how this survived in two places.

Correct `doc/src/circuit/lte_dae.rst:154-165`: it says the classic estimate costs
"roughly an order of magnitude" of accuracy, when in fact it disables the controller
outright (zero rejections, pinned at `max_step`, needing ~3e10x tighter tolerances to
react). The document was right that something was wrong and wrong about how much, which
is the more dangerous of the two.

**Gate 5-1.** `test_jaxtransient.py` passes; if jax is importable, the JAX and numeric
LTE estimates agree to a stated tolerance on the same input.
OUTCOME:

**Gate 5-2.** Doc build: **build succeeded, 2 warnings, 0 ERROR lines**. Check
`grep -cE 'ERROR'` separately from the warning grep — an `ERROR:` line cannot be matched
by `grep -i warn`, which is how two of them stayed invisible before.
OUTCOME:

---

## Stage 6 — only now, re-run the leapfrog harness

`benchmarks/nonlinear_leapfrog_sweep.py` (two-tone IM3, `b2ab5bb`). This stage is
**blocked** on the separate `leapfrog_5th_order` instability investigation: with two
right-half-plane poles there is no steady state and no transient can succeed, so running
it earlier only produces a number that is growth, not distortion.

**Gate 6-1 (the structural gate that failed before).** Two-tone spectrum must show
discrete lines at bins 9/10/11 of a 100 us window, with `|vout|` at bin 10 within 10x of
the perturbation's 6.7e-05 V — not the 6.3e+03 V and smooth pedestal seen at `b2ab5bb`.
OUTCOME:

**Gate 6-2 (the artifact floor).** IM3 read at the *input* node must be far below the
IM3 read at the output. The input is two ideal series sources and cannot intermodulate,
so whatever appears there is artifact. Declared success: input-node IM3/fund at least
10x below the output's. At default tolerances it was 1.58e-04 against a 2.34e-04 signal
— i.e. useless — so this gate is doing real work.
OUTCOME:

**Gate 6-3 (transient convergence).** IM3 must agree between two tolerance settings
before any value is quoted. A single-tolerance transient number is not evidence.
OUTCOME:

**Gate 6-4 (the actual deliverable).** Transient IM3 vs perturbation IM3 at U^13.
Declared success: agreement to a few percent. **Agreement to 1e-6 would be a red flag,
not a triumph** — it would mean the two paths are not independent and should be
suspected of sharing code.
OUTCOME:

---

## Order, and what ships with what

1 -> 2 must be adjacent (1 knowingly breaks a test that 2 repairs; do not leave the tree
red between commits). 2b follows 2. 3, 4, 5 are independent of each other and of 1-2b,
except that stage 4's tests are most convincing when written against the pre-1 code. 6 is
blocked externally.

**Approved scope, 2026-07-30: stages 1 through 5.** Stage 6 stays blocked on the
`leapfrog_5th_order` instability investigation.

Each stage commits with its explanation and its measured gate outcomes in the message —
not a documentation pass afterwards. Negative results get recorded in the same voice as
positive ones; a refuted premise is a result, and an unrecorded one gets re-attempted.

**Do not `git push`.** 58 commits are already unpushed on `cna-jax-vectorization`;
pushing is the maintainer's call.
