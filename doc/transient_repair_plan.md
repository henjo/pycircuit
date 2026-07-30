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
OUTCOME: **PASSED.** Stage 1-3's reference case, reltol 1e-4. `err` is the
start-up-free error against `expm`.

| | steps | rejections | err |
|---|---|---|---|
| explicit `'classic'`, before 2b | 40 | 4 | 1.6115e-03 |
| explicit `'classic'`, after 2b | 40 | 4 | 1.6115e-03 |
| explicit `'ywr'`, before 2b | 37 | 4 | 2.0018e-03 |
| explicit `'ywr'`, after 2b | 37 | 4 | 2.0018e-03 |
| `Gear2Integrator()`, before 2b | 40 | 4 | 1.6115e-03 |
| `Gear2Integrator()`, after 2b | **37** | 4 | **2.0018e-03** |

Both explicit forms are bit-identical across the change, so the *only* thing that
moved is which estimate the default selects — the attribution stage 2b was split
out to protect. The default is now bit-identical to explicit `'ywr'`, and
corrected `'classic'` remains reachable and, on this case, slightly *more*
accurate for slightly more steps (1.61e-03 in 40 against 2.00e-03 in 37) — which
is the accepted cost of the decision, visible in the measurement.

**Gate 2b-2.** Full suite `-m ""` at 715 passed, 6 skipped, 0 failed. Anything that moves
here is a test that was implicitly depending on the default, which is worth naming.
OUTCOME: **PASSED — 716 passed, 6 skipped, 0 failed in 450.78 s** (-22.1% against
the 578.49 s baseline). **Nothing moved**, so there is nothing to name: no test in
the suite was depending on `Gear2Integrator()`'s default LTE formula. Which is
itself worth recording, since it means the suite would not have caught the switch
either way.

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
OUTCOME: **PASSED.** `VSin(va=1, freq=1e4)` -> `R 1k` -> `C 1uF` with `L 1uH` and
`R 1` to ground; `tend` 200 us, `max_step` 5 us, so `Sin.next_event` fires 8 times
in the run.

| unchecked accepts | before | after |
|---|---|---|
| gear2 | 11 | **1** |
| trap | 8 | **1** |
| euler | 9 | **1** |

Exactly one, for every integrator — the genuine first step. Note the before-counts
exceed the 8 breakpoints: rejected steps re-enter the loop with the flag still
armed, so the hole was slightly wider than "one per quarter period".

**Gate 3-2 (discontinuities still work).** A `VPulse` edge must still be integrated
without a rejection loop — the breakpoint is there because the discontinuity is real,
and a small step is genuinely wanted. Declared success: `test_breakpoints.py` passes and
the pulse edge takes no more than 3 rejections.
OUTCOME: **PASSED.** `test_breakpoints.py` passes (and so do
`test_analysis_transient.py`, `test_minstep.py`, `test_uic.py`: 15 passed). On its
`VPulse` circuit (1 us edges, `max_step` 10 us):

| | steps | total rejections | max rejections at a pulse edge |
|---|---|---|---|
| gear2 before | 347 | 6 | 3 |
| gear2 after | 366 (+5.5%) | 16 | **3** |
| euler before | 1217 | 7 | 3 |
| euler after | 1347 (+10.7%) | 13 | **3** |

No rejection loop: 3 is the `MAX_REJECT` cap and it is reached at an edge both
before and after, so the discontinuity is being handled the same way — the edge
step is simply now *checked* before being taken. The 5-11% step-count cost is the
price of closing the hole, and it is paid where a discontinuity actually is.

**Recorded because it is counter-intuitive:** on the `VSin` case rejections went
*down*, 98 -> 47 for gear2. Checking the post-breakpoint step means the controller
enters the next step with a sane `h` instead of a full `max_step` one it then has
to fight back down.

**Gate 3-3.** Full suite `-m ""` at 715/6/0, runtime within 20% of baseline.
OUTCOME: **PASSED — 716 passed, 6 skipped, 0 failed in 453.01 s** (-21.7% against
the 578.49 s baseline; `test_stress_stiff_rlc_pulse` 70.78 s against 91.35 s).
Closing the hole costs nothing at suite level: no test anywhere depended on a
post-breakpoint step being accepted unevaluated.

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
OUTCOME: **PASSED, both directions demonstrated by actually reverting.** 18
parametrized cases across the three new test functions (4 integrator configurations:
`euler`, `trap`, `gear2-classic`, `gear2-ywr`).

| `pycircuit/circuit/integrator.py` | result |
|---|---|
| at HEAD | **18 passed** in 12.46 s |
| restored to `3f2627e^` (pre-stage-1), file otherwise untouched | **3 failed, 15 passed** in 12.41 s |

The three failures are exactly:

```
test_compute_lte_order_and_scale[gear2-classic-Gear2Integrator-classic-2.0-0.666…]
test_step_is_actually_rejected_on_stiff_case[gear2-classic]
test_step_count_and_error_respond_to_reltol[gear2-classic]
```

**One failure per new test function, and every one of them on the `gear2-classic`
parametrization alone.** That is the outcome worth having rather than a bare "they
fail": `euler`, `trap` and `gear2-ywr` pass in *both* states, so the tests discriminate
the specific defect instead of being broadly sensitive to the file changing. A test
that had failed on all four configurations would have proved much less.

The reverted file was restored to HEAD and the working tree confirmed clean afterwards.

**Gate 4-2.** Full suite `-m ""`, new tally recorded (it will exceed 715).
OUTCOME: **PASSED. 734 passed, 6 skipped, 0 failed in 511.93 s.** Up 18 from the
716 after stage 2 -- the 18 new parametrized cases (3 test functions x 4 integrator
configurations, minus the combinations the `_LTE_CASES` table does not enumerate).
Runtime is -11.5% against the same-box 578.49 s baseline, so nothing needed marking
`slow` and the maintainer's runtime decision never had to be invoked.

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
OUTCOME: **PASSED.** `jaxtransient.estimate_lte` carried the same defect for BOTH
2nd-order methods -- a second divided difference of the charge (which is `q''`) scaled
by `h^3`, where `q'''` scaled by `h^2` is required.  Its comment even asserted "h^3
q''' via a second divided difference of the charge", which is self-contradictory.
Fixed to match stage 1: `q'''` as twice the second divided difference of `g` from
`iq_history`.  The `euler` branch was correct and is untouched; it sets the
charge-domain convention the other two now follow.

Both implementations fed the SAME q/iq/h history, `q(t)=sin(2*pi*1e6*t)`, `h` swept
4 ns -> 0.25 ns:

| method | JAX/CPU at h=4e-9 | at h=2.5e-10 | `(JAX/CPU)/(2h)` |
|---|---|---|---|
| euler | 8.025e-09 | 5.001e-10 | 1.00309 -> **1.00021** |
| trap | 2.000e-09 | 1.250e-10 | **0.25000 at every h** |
| gear | 2.667e-09 | 1.667e-10 | **0.33333 at every h** |

`p+1` agrees for all three (2.0, 3.0, 3.0).

**The ratio must carry a factor of h, and an earlier version of this check got that
wrong and reported FAIL.** The CPU returns `Eg` in CURRENT units (the controller
applies `J^-1` afterwards) while JAX is charge-domain by design; charge = current x
time.  Testing the raw ratio for constancy therefore tests the wrong thing.  Divided
by `2h` it is constant to **0.0000%** for trap and gear over a 16x range in `h`, and
converges to 1 for euler (the residual drift is the divided difference's own O(h)
error).  Fixed per-method constants with zero drift is far stronger evidence than a
single-point match would have been.

**Out-of-scope defect found here and NOT fixed -- reported instead.**
`outer_time_loop` receives `reltol`/`abstol` and threads them into
`newton_inner_loop` (jaxtransient.py:360), but **neither LTE call site passes them**:
`ywr_error_ratio(...)` at :372 and `lte_error_ratio(lte, q_curr)` at :378 both fall
back to hard-coded `trtol=7.0, lte_rel=1e-3, lte_abs=1e-6`.  So in the JAX transient
the user's tolerances reach the Newton convergence test but **never reach the
step-size controller**, on BOTH paths including the default `'ywr'` one.  This is the
same class of defect as the Gear2 bug -- a controller that cannot be tuned -- and
unlike the copied LTE expression it is LIVE rather than dormant.  Outside the
approved stage-5 scope ("fix the JAX copy"), so it is recorded here as a decision for
the maintainer rather than actioned.

Also worth noting for whoever takes that on: `lte_formula` does not mean the same
thing in the two implementations.  On the CPU it selects between two algebraic
formulas inside one `compute_lte`, both returning `Eg`.  In JAX it selects between two
*different functions returning different kinds of quantity* -- `ywr_error_ratio`
returns an already-normalised ratio and needs an extra `G` evaluation plus a solve per
step, while `estimate_lte` returns an un-normalised charge-domain LTE.  They share no
contract, which is exactly why they were free to drift.

**Gate 5-2.** Doc build: **build succeeded, 2 warnings, 0 ERROR lines**. Check
`grep -cE 'ERROR'` separately from the warning grep — an `ERROR:` line cannot be matched
by `grep -i warn`, which is how two of them stayed invisible before.
OUTCOME: **PASSED on its declared terms -- and it exposed a stale-number trap that
matters more than the gate did.**

Declared part: `build succeeded, 2 warnings`, `grep -cE 'ERROR'` = **0**. A
`grep -ciE warning` returns **5**, which disagrees with the summary; resolving that
disagreement rather than trusting either number: 2 are venv `RuntimeWarning`s about
`sys.prefix` (not sphinx diagnostics), 2 are the known `pycircuit.post.cds` /
`cds.skill` autodoc import failures that form the clean baseline, and 1 is the summary
line itself. The summary is authoritative, as the standing rule says.

**THE TRAP.** `doc/src/circuit/distortion_ddd.rst:500-560` holds a live `exec-rst`
block calling `bc.leapfrog_5th_order()`, and `leapfrog_redo_plan.md` listed it under
"Documents that DO self-update" -- the "never type a measured number into prose" rule
apparently paying for itself. **It does not self-update.** The block did execute (the
HTML contains a rendered table, and `grep -c 'exec-rst|bc.leapfrog_5th_order()'` on the
built page is 0, so it did not fall back to echoing its source). But sphinx reuses
cached doctrees for source files that have not themselves changed, and
`distortion_ddd.rst` did not change when the fixture did. So the page still shows the
values computed from the **unstable** fixture.

Verified by recomputing the same quantities against the fixed fixture rather than by
assuming:

| drive | node volts (built page) | recomputed | HD3 (built page) | recomputed |
|---|---|---|---|---|
| 0.01 | 6.329e-05 | **7.1404e-05** (+12.8%) | 5.919e-11 | **6.2704e-11** (+5.9%) |
| 3 | 1.847e-02 | **2.0665e-02** (+11.9%) | 4.665e-06 | **4.8100e-06** (+3.1%) |

So a live block is only self-updating with respect to **its own source**, never with
respect to the code it calls. A fixture change elsewhere leaves it silently stale, and
"build succeeded, 2 warnings" is emitted either way. `leapfrog_redo_plan.md` has been
corrected: stage T3 must force a clean rebuild (`sphinx -E`, or delete
`doc/build/doctrees`) and then diff the numbers, not merely rebuild.

Incidental, and reassuring for the redo: the shifts are 3-13%, not orders of magnitude,
so the qualitative conclusions are likely to survive T2-2 even though every number in
them moves.

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
