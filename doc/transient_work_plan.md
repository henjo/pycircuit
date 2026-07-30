# Transient subsystem — the work plan

**Status: written 2026-07-30. NOTHING HAS RUN. Every OUTCOME line is deliberately blank
and must be filled in from a real run, never predicted.**

Findings and evidence: `doc/transient_review.md`. That document is where an argument about
*whether* something is broken belongs, and it carries the measurements. This one is the
order of work and the gates.

Evidence scripts: `benchmarks/transient_review/` (probes, not gates — see its README).

Prior work this builds on, all complete: `doc/transient_repair_plan.md` (the Gear2 LTE
repair, stages 1-5, gates recorded).

---

## How this plan is executed

These are not optional and they are what the plan is shaped around.

1. **Gates are declared before the stage runs.** Each has a success criterion written in
   advance and a blank `OUTCOME:`. Fill it from a real run. **A refuted premise is a
   result** — record it in the same voice as a positive one. Never pre-fill, never
   predict, never write a number that was not measured.
2. **The explanation ships in the same commit as the code.** Not a documentation pass
   afterwards. Writing the explanation is what forces you to have one, and the prose then
   becomes an anchor to check the code against — defects surface when the two disagree.
3. **Verify the output, not the exit code.** "Tests ran", "build succeeded", "no errors"
   are not verification. Read the tally line. For a doc build, require the positive
   `build succeeded` line and grep the built HTML for a *computed* value.
4. **Never type a measured number into prose** where the build can generate it. A pasted
   number is correct exactly once.
5. **Commit messages explain *why* and name the evidence.**
6. **Never `git push`.** That is the maintainer's call.

**Baselines to beat.** Re-measure both on the executing machine before stage 1 rather
than trusting these; the review recorded that this box's suite is ~17% slower than the
one the 492 s figure came from.

- full suite `-m "" --timeout=400`: **734 passed, 6 skipped, 0 failed** (at `e9dd894`)
- doc build: **build succeeded, 2 warnings, 0 ERROR lines**
- `pytest.ini` deselects 17 `slow` tests. **`-m "" ` is mandatory** before reporting any
  suite result — the deselected ones are the only tests comparing an analysis against a
  time-domain reference.

**A standing risk, stated once.** Several stages change step counts and therefore
runtimes. `test_stress_stiff_rlc_pulse` is already ~54% of suite runtime (`architecture.md`
P15). Every stage's suite gate records runtime, and the maintainer's standing decision
applies: if runtime regresses past 20%, mark the worst offender `@pytest.mark.slow` rather
than loosening a test — a fast test that no longer stresses the controller is worse than a
slow one that does.

---

# STAGE 0 — close the open questions before implementing

The review left specific things unreviewed, and made decisions that are the maintainer's
rather than the implementer's. Doing these first is cheaper than discovering them in the
middle of stage 7.

## 0.1 Reviews still to do

Each as a **named lens**, like the four that produced `transient_review.md`. One specific
lens finds more than three generic passes.

**0.1a — `jaxtransient.py`, as a JAX/compiler engineer.** The four-lens review read its
Newton loop, LTE and step control but explicitly *not* its chunking/batching machinery:
`solve_batched`, `run_chunk`, the TLine ring buffer under `vmap`, the `10000` literal
repeated six times including in modulo arithmetic. Stage 9 proposes consolidating this
file; that cannot be planned without knowing what is in it. Ask specifically: what
genuinely requires static shapes and traceable control flow, and what is duplicated only
because nobody factored it?
OUTCOME:

**0.1b — `semiconductors.py` and `mos.py`, as a device-modelling expert.** Stage 5 adds
junction limiting to `BJT`/`JFET`/`ZenerDiode` on the strength of a proposed diff that was
never tested. Also open: no charge model on any semiconductor (`BJT.C(x) == 0`), the
`Varactor` clamp that makes C fall to zero above the knee instead of extrapolating,
`MOS_ACM.__init__` calling a class not in its MRO, and the absence of any large-signal
MOSFET. Ask: what is the minimum device set that makes CMOS and bipolar transients
expressible, and is the `Semiconductor` base the right seam?
OUTCOME:

**0.1c — `shooting.py`, as an RF/steady-state expert.** PSS/PAC are advertised in the
analysis inventory and the review found them unusable beyond toys: `method` Parameter
never read, no limiter, a dense `np.linalg.inv` per timestep inside the shooting Newton,
and a PAC matrix that is ~150 TB at N=137/M=1000. Ask: repair, rewrite against the
`Integrator`/`LinearSolver` seams stage 7 creates, or withdraw and document as
unimplemented? A thin advertised feature is worse than an absent one.
OUTCOME:

**0.1d — `_solve_coupled`, as a numerical analyst.** The third transcription of the time
loop (`transient.py:429-551`). Read for consistency but its Fang-2013 co-determination
loop and `MAX_LTE_ITERS = 10` were never analysed. It ignores breakpoints, `uic`,
`fixed_timestep` and any injected step controller. Ask: is this worth keeping at all, or
is it a research prototype that should be deleted or moved behind an explicit flag?
OUTCOME:

## 0.2 Measurements that gate later stages

**0.2a — Does the trapezoidal `(-1)^n` contamination actually affect the IM3 harness?**
`benchmarks/nonlinear_leapfrog_sweep.py` sets `TrapezoidalIntegrator()` and drives two
`VSin` tones, i.e. it hits the exact seeding mechanism (order drop at a quarter-period
breakpoint). The integrator was chosen on a step-count comparison a contaminated estimator
could have produced.
**Gate:** compare the harness's step count and output against the same run with
`Gear2Integrator()` and with a suppressed-breakpoint `VSin`. Declared success: the step
counts agree within 20% and the waveform within the IM3 measurement's own tolerance. If
they do not, **the recorded 10x integrator speedup is void** and must be re-measured after
stage 4.
OUTCOME:

**0.2b — What does the JAX charge-domain LTE path cost, and should it exist?**
`estimate_lte` + `lte_error_ratio` avoid a per-step `G` evaluation and linear solve, which
is their one genuine advantage over `ywr_error_ratio`. But `lte_abs = 1e-6` is a *voltage*
tolerance applied to a **charge**, so that path currently never rejects a step.
**Gate:** measure the per-step cost of the extra `G` + solve as a fraction of the whole
step, on a circuit of realistic size. Declared decision rule: **if under 10%, delete
`estimate_lte` and make `lte_formula` select an algebra only** (as it does on the CPU);
if over, fix its tolerance to a charge-referenced one and keep both.
OUTCOME:

**0.2c — Re-measure the suite baseline on the executing machine.** Both tallies and both
runtimes, `-m ""`, before any change.
OUTCOME:

## 0.3 Decisions for the maintainer

**0.3a — `vabstol` serves two roles.** It is Newton's x-tolerance on node rows
(`transient.py:142`) *and* the LTE tolerance on node rows (`:294`). The 1e-12 -> 1e-6 fix
was reasoned about only as a step-control knob, so it also loosened Newton's node
convergence by 10^6, unmeasured; and `DC.vabstol` is still 1e-12, so the operating point
seeding every transient is solved 10^6 times tighter than any step after it. All four
lenses flagged this.
**Options:** (i) two Parameters — `vabstol` for Newton, `lte_vabstol` for the controller,
with DC and Transient sharing the former; (ii) one Parameter with the controller's ratio
made explicit rather than coincidental; (iii) adopt SPICE's split properly, adding
`chgtol` and making the LTE criterion charge-flavoured, which is the largest change and
the most standard.
**Recommendation: (i) now, (iii) as a later stage if `chgtol` is wanted.**
DECISION:

**0.3b — Does `Gear2` keep `'ywr'` as its default?** It was set to `'ywr'` belt-and-braces
when `'classic'` was repaired. The evidence now runs the other way: `'classic'` measures
asymptotically exact (1.000282 against 2/9) while `'ywr'` is 3/4 biased *and* step-ratio
dependent (swinging 16x across ratios against classic's 4x), and end-to-end `'ywr'` needed
57 rejections and 4 force-accepts where `'classic'` needed 29 and 0.
**Recommendation: flip the default to `'classic'` in stage 4**, and delete the "5/6"
claim in the `integrator.py` comment, which is measurably wrong.
DECISION:

**0.3c — Scope of stage 10.** The missing-analyses list (DC sweep, `.ic`/`.nodeset`,
netlist import, large-signal MOSFET) is weeks of work and is a product decision, not an
engineering one.
DECISION:

**Stage 0 exit criterion:** every OUTCOME above filled, every DECISION answered. Stages
1-3 may start before 0.1a-d land (they touch none of that code); **stage 4 is blocked on
0.2a, stage 5 on 0.1b, stage 7 on nothing, stage 9 on 0.1a, and stage 11 on 0.1c.**

---

# STAGE 1 — stop the silent failures

The two defects that return confidently wrong answers. Nothing else in this plan matters
if a run can quietly report a circuit that was never biased.

**Work.** (a) `transient.py:248-252` and `:443-444` — do not substitute zeros on a failed
DC; raise, with a message naming the escape (`uic=True` or an explicit `x0=`). (b) Narrow
both `except Exception` clauses (`dcanalysis.py:130-133`, `transient.py:169-172`) so a
`ZeroDivisionError` from a source model is not reported as a convergence failure. (c) Pass
`self.epar` to `cir.C/q/i/G/u` in the transient (`transient.py:211-215` and the four
equivalents in `_solve_coupled`); add `bypasstol` to `defaultepar`. (d) Construct the inner
`DC` with the transient's own toolkit, epar, tolerances, maxiter, nrsolver, scaler and
refnode.

**Docs in the same commit:** a section in `doc/transient.rst` (or the transient module
docstring) stating what happens when the operating point fails and how to ask for a
zero start deliberately.

**Gate 1-1 (the failure is visible).** A circuit with no DC solution must raise, with a
message that names the circuit condition. Declared success: it raises; the message
contains the word `uic`; and no waveform is returned.
OUTCOME:

**Gate 1-2 (epar actually arrives).** Instrument a device's `i()` during a transient with
`epar.T = 400`. Declared success: the device sees T = 400 in **both** the inner DC and
every transient step. Currently it sees 300 in both.
OUTCOME:

**Gate 1-3 (bypass is connected).** `bypass=True, bypasstol=1e-2` and `bypass=False` must
produce *different* step/evaluation counts. Currently they are identical because
`bypasstol` is missing from `defaultepar` and every device takes its `except
AttributeError` branch.
OUTCOME:

**Gate 1-4 (blast radius).** Full suite `-m ""`. Declared success: 734/6/0, **or** a list
of failures each individually explained as a test that depended on the silent fallback —
which is information worth having, not a reason to restore it.
OUTCOME:

---

# STAGE 2 — performance, without touching numerics

10.5x measured and prototyped, waveform drift exactly `0.00e+00`. Do it before the linear
solver: at n=137 the solve is 2.1% of runtime and assembly is 96%, so replacing the solver
first would optimise the wrong thing.

**Work, in this order** (each independently gated, so a regression is attributable):

**2a. Single-thread the BLAS around the solve.** `threadpoolctl.threadpool_limits(1)`, or
document loudly. Threaded LAPACK on a 136x136 matrix measured **0.238 ms single-threaded
against 4.462 ms with 4 threads** — thread-spawn overhead against 1.7 MFLOP of work.
**Gate 2a:** the same transient, same machine, timed both ways. Declared success: >= 2x,
waveform drift `0.00e+00`. Record whether `threadpoolctl` is an acceptable dependency; if
not, the fallback is documentation plus a warning when >1 thread is detected.
OUTCOME:

**2b. Assembly: hoist the per-element probes, batch the scatter.** `circuit.py:1290-1402`.
`hasattr(self.toolkit, 'add_at')` fires once per element per stamp; `NumericToolkit` lacks
it, so `Toolkit.__getattr__` builds a formatted error string and **raises 255 000 times**
— 34% of total runtime is attribute-lookup machinery. Hoist the probes out of the loop and
scatter once via `np.bincount` (keep `np.add.at` for object/complex dtypes).
**Gate 2b:** `i`, `q`, `u`, `G`, `C` each compared against the current implementation on a
non-trivial circuit. Declared success: **max abs difference exactly 0.0**, and >= 1.8x on
assembly.
OUTCOME:

**2c. Stop re-assembling after Newton converges.** `transient.py:218-219` re-evaluates
`func(x)` at the point Newton just converged to, and `StandardNewton` discards the `(F, J)`
it already has. `:364` and `:409` each recompute `cir.q(x)` at that same `x`. Measured
3.17 assemblies per accepted step where ~2.17 are needed.
**Gate 2c:** instrument the per-step assembly count. Declared success: `G`/`C`/`i`/`u` at
~2.17 per step, `q` correspondingly reduced, waveform drift `0.00e+00`. Callers using
`provided_function` keep an opt-in exact-at-x recompute.
OUTCOME:

**Gate 2-final (compound).** Full suite `-m ""` at 734/6/0, and the end-to-end speedup
recorded against the stage-0 baseline. Declared success: >= 5x on the review's benchmark
circuit with drift `0.00e+00`. **If drift is not exactly zero, stop** — this stage is
defined as behaviour-preserving, and a nonzero drift means something else changed.
OUTCOME:

**Docs in the same commit:** a short "performance notes" section recording *why* assembly
dominates and what the three changes do, so the next person does not re-derive it.

---

# STAGE 3 — the first step

`reltol` currently does not control transient accuracy: the run opens at `max_step` and
that step is accepted unevaluated, and its error is the maximum of the whole run —
identical to four digits across every integrator and 1e-3..1e-6 reltol.

**Work.** Add a `firststep` Parameter (default `None`). When unset, open at
`timestep * 1e-3`. Apply to `_solve_coupled` identically. The principled version is a
Hairer-style estimate from `q'`/`q''` at the operating point; **implement the ramp first
and record whether the Hairer estimate is worth the complexity** — that is a measurement,
not an assumption.

**Gate 3-1 (reltol becomes a knob).** On a circuit with an analytic solution, global error
must fall monotonically with reltol over 1e-3..1e-6. Declared success: at least 20x total
reduction. Currently: 0x (bit-identical).
OUTCOME:

**Gate 3-2 (cost).** Declared success: step count rises no more than 20%. Measured on the
review's circuit: +5% for Trapezoidal.
OUTCOME:

**Gate 3-3 (the order effect appears).** With the first step ramped, a 2nd-order method
must beat backward Euler on the same circuit and tolerance. Currently they are identical,
which is the clearest single symptom of the defect.
OUTCOME:

**Gate 3-4.** Full suite `-m ""`, runtime recorded.
OUTCOME:

---

# STAGE 4 — step-control correctness

Blocked on 0.2a. Four defects in the controller and the estimators, grouped because they
interact and separating them would make the interactions unattributable — but **each gets
its own measurement before the next is applied.**

**4a. `PIController` gains.** `stepcontroller.py:155` uses Gustafsson's numerators
undivided; they must be `0.3/k` and `0.4/k`. Spectral radius **1.12 (k=2), 1.78 (k=3)** —
unstable, converted by the clamp into a permanent period-2 limit cycle. The computed
`exponent = 1.0/p` at `:139` is dead code. Also update `last_err` on the rejection path,
or fall back to pure-I for the step after a rejection.
**Gate 4a:** the step-size recursion must converge rather than cycle. Declared success: on
a smooth problem, h settles to within 5% of a fixed point instead of alternating 2:1.
OUTCOME:

**4b. `MAX_REJECT` force-accept.** `transient.py:390` grows `dt` **10x** after accepting an
over-tolerance step. Variable-step BDF-2 is zero-stable only for ratio < 1+sqrt(2) =
2.414214; at 10x the parasitic root is 4.76. `Gear2(ywr)`, the shipped default, is the
only configuration measured reaching this path. Replace growth with an order drop, and
**warn** — an unbounded accepted truncation error must not be invisible.
**Gate 4b:** on the review's circuit, no step ratio exceeds 2.414, and any force-accept
emits a `RuntimeWarning` naming `t` and `h`.
OUTCOME:

**4c. Backward Euler's variable-step bias.** `integrator.py:83` — `gn - gn_1` approximates
`((h1+h2)/2) q''`, not `h1 q''`. Rescale by `2*h_curr/(h_curr+h_last)`.
**Gate 4c:** est/true flat across step ratio 0.25..4. Declared success: within 5% of 1
across the range. Currently 0.62..2.47.
OUTCOME:

**4d. `Trapezoidal(lte_formula='ywr')` is order-inconsistent on non-uniform grids.**
`integrator.py:123` — the plain second difference leaves an O(h) term where the truncation
error is O(h^2); measured 112x too large at ratio 0.25, -436x at 4.0, scaling as 1/h. Its
correct divided-difference generalisation is algebraically identical to the `classic`
branch, so **delete the branch and keep `'ywr'` as an alias** rather than creating a
duplicate.
**Gate 4d:** est/true within 5% of the `classic` column across ratio 0.25..4.
OUTCOME:

**4e. `Gear2.check_order_drop` guards the wrong direction.** `integrator.py:158` fires on
*shrink*, which is unconditionally zero-stable; there is no guard on *growth*, which is
the unstable one — which is how 4b slipped through. Replace with an upper-ratio guard.
**Gate 4e:** a growth ratio above 2.414 triggers the order drop; a shrink does not.
OUTCOME:

**4f. Default `lte_formula`.** Per decision 0.3b.
**Gate 4f:** whichever default is chosen, record rejections and force-accepts for all four
integrator/formula combinations on the same circuit, so the choice is evidenced.
OUTCOME:

**4g. The trapezoidal `(-1)^n` contamination.** The largest item in this stage and the one
that blocks the IM3 harness. Three tiers, and **tier (a) alone may suffice — measure
before doing (b)**:
  - (a) `Sin.next_event` (`func.py:29-35`) plants a breakpoint every quarter period on a
    C-infinity waveform, and each one re-arms the order drop that seeds the parasitic
    mode. Return `td` then `inf`. Breakpoints are for discontinuities.
  - (b) If contamination persists, difference a mode-free quantity: `d_n = (q_n -
    q_{n-1})/h` is `(g_n + g_{n-1})/2` and carries no alternating component. Needs one
    more charge in history and an `h_last2` interface change — that interface change is
    the real cost.
  - (c) A `TRBDF2Integrator` (gamma = 2-sqrt(2)): L-stable, and its embedded estimator is
    genuinely `q'''`-based so it dispenses with the problem entirely. The `Integrator` ABC
    docstring already names it.
**Gate 4g-1:** with (a) alone, count Newton evaluations forced to order 1 on a plain
`VSin`-driven RC. Declared success: **0**, against the measured 120 of 1236.
OUTCOME:
**Gate 4g-2:** est/true for Trapezoidal must not scale as 1/h. Declared success: ratio
bounded within 2x of 1 as h falls 1e-2 -> 1e-4. Currently -10 -> -2976.
OUTCOME:
**Gate 4g-3:** re-run 0.2a. Declared success: the harness's integrator choice is either
confirmed or corrected, and the recorded 10x is either upheld or replaced.
OUTCOME:

**4h. `fixed_timestep=True` does not fix the timestep.** `transient.py:415-416` restores
`dt` only when *not* fixed-step, so breakpoint truncation is permanent. Measured: expected
~20 steps, got **19 002**, dt collapsing to 3.276e-22 s.
**Gate 4h:** a `VPulse`-driven fixed-step run takes `tend/timestep` +- 1 steps.
OUTCOME:

**Gate 4-final.** Full suite `-m ""`, runtime recorded. **Expect test churn here** — this
stage changes step counts by design, and any test asserting a step count is exposed.
OUTCOME:

**Docs in the same commit:** `doc/src/circuit/lte_dae.rst` gains the variable-step story
(error constants, the step-ratio bound, why `'ywr'` trapezoidal was withdrawn). This is
the document that already understated the Gear2 defect once; it should not understate
these.

---

# STAGE 5 — convergence: limiting and the aid ladder

Blocked on 0.1b.

**Work.** (a) A `junctions` class attribute plus a shared `limit()` on `Semiconductor`, so
`BJT`, `JFET` and `ZenerDiode` get `pnjlim`. (b) An `_expl()` clamp (SPICE's `EXPMAX`
treatment) so `exp()` overflow cannot reach the central-difference Jacobian and become
`nan`. (c) **Fix `SubCircuit.limit`** (`circuit.py:1194-1200`): `subx = x[nodemap]` is
fancy indexing, therefore a copy, and the return value is discarded — a limiter that
writes `x` has no effect. `Diode` does not notice because it keeps `_vlim` and linearises
around it, which is exactly why the bug survived. (d) Consider making `limit` return the
limited `x` (state-free) rather than mutating device state — **survey this separately**, it
has a wide blast radius into `elements.py` and is the prerequisite for the JAX path ever
having limiting.

**Gate 5-1 (the failing case converges).** A common-emitter BJT stage with a
voltage-driven base must reach a DC operating point at base resistances 10, 100 and 1000
ohm. Currently all three fail.
OUTCOME:

**Gate 5-2 (no `nan` reaches the Jacobian).** Instrument for non-finite entries during the
continuation. Declared success: none, at any gmin or source-stepping rung.
OUTCOME:

**Gate 5-3 (the limiter is actually applied).** A test that fails if `SubCircuit.limit`'s
return value is discarded — i.e. a device whose limiter clamps `x` rather than keeping
private state.
OUTCOME:

**Gate 5-4 (no regression on circuits that already converge).** Iteration counts on the
existing nonlinear tests must not increase by more than 20%. **This is the gate the
review explicitly did not measure**, so it is the one most likely to bite.
OUTCOME:

---

# STAGE 6 — diagnostics and statistics

The stage that makes every future failure self-explaining. It is deliberately after 1-5,
because those generate the failures worth reporting well.

**Work.** (a) Classify a structurally singular matrix **before** continuation is attempted
— an LU zero pivot names the row, and the row maps to a node via `cir.get_node_name`, so
"no DC path to ground" is reachable. Currently three layers of re-wrapping turn it into
`Source Stepping failed at lambda=0.0`. (b) Name the offending node on a convergence
failure: `xdiff`, `F`, `conv_x`, `conv_f` are all in scope one `argmax` away, and
`Circuit.name_state_vector` already exists. (c) A statistics object on the result:
accepted/rejected steps, Newton iterations, LTE rejections, force-accepts, order drops,
min/max step, breakpoints hit, and time in load vs solve. `solve_system` already returns
the iteration count and `transient.py:158` discards it.

**Gate 6-1.** A floating node produces a message naming the node and the condition.
OUTCOME:

**Gate 6-2.** A non-convergent circuit produces a message naming the worst node, its
residual, and the tolerance it missed.
OUTCOME:

**Gate 6-3.** Statistics are populated and non-zero on a run that rejects steps, and the
force-accept counter from 4b appears there.
OUTCOME:

**Docs in the same commit:** a "diagnosing a failed run" section — this is the
documentation with the highest ratio of user value to effort in the whole plan.

---

# STAGE 7 — the linear solver

Independent of stages 1-6 and may run in parallel. The case is **not** the 2.1% at
n=137 — it is the wall: dense MNA is 200 MB at n=5000 and 3.2 GB at n=20 000, so
pycircuit cannot run an ordinary mixed-signal block at all.

**7a. Take the reference node out of the matrix.** `analysis.py:63-69` `remove_row_col`
does an O(n^2) `np.delete` copy of J on **every Newton iteration** purely to drop the
ground row and column. This is what forces a dense copy and destroys any cached pattern —
so it comes first, before any solver work.
**Gate 7a:** identical results on the transient tests; `np.delete` absent from the
profile.
OUTCOME:

**7b. A `LinearSolver` strategy object**, in the shape the codebase already uses for
`Integrator`/`StepController`/`Scaler`/`NonLinearSolver` — `analyze` / `factor` / `solve`
split so the factorisation survives the Newton iteration. `DenseLU` (scipy
`lu_factor`/`lu_solve`) below ~200 unknowns, `SuperLU` above; measure the crossover rather
than guessing it (dense wins at n=29, sparse wins 5x at n=137, 15x at n=542).
**Gate 7b-1:** results identical to the dense path on every existing transient test.
**Gate 7b-2:** the step controller's `J^-1 Eg` solve (`stepcontroller.py:60`) reuses the
factors — a third of all linear solves are currently redundant re-factorisations.
**Gate 7b-3:** a circuit at n ~ 2000 that currently exhausts memory now runs.
OUTCOME:

**7c. KLU with `klu_refactor`, behind an optional import.** scipy's `splu` recomputes
COLAMD on every call — **94% of its cost**. Reusing the symbolic phase needs KLU's
analyze/factor/refactor split; that is what Xyce does and it is the single biggest
transient lever in production tools. Discover it the way `_ginac_ext` and `symengine`
already are.
**Gate 7c:** measured against 7b on the same circuits; declared success >= 3x on the
factor+solve at n >= 500. **Rejected alternative:** a native Markowitz in pure Python —
the ordering machinery already exists in `ddd.py:1023-1196`, but a Python LU over 3000
nonzeros will lose to compiled SuperLU. *Reconsider if* a Cython/numba build step becomes
acceptable anyway for stage 2's assembly work.
OUTCOME:

**7d. Delete `pybsmatrix.py`** (340 unreferenced lines, no pivoting, and a `fbsub` whose
division sits inside the wrong loop so it cannot be correct), and **fix
`test_sparse_toolkit.py` to construct the circuit with the toolkit under test** — it
currently passes while never exercising the sparse path. Fix the test first, then decide
whether `_sparse_numeric` is worth keeping given it is 4x slower than dense.
OUTCOME:

**Docs in the same commit:** the sparsity and ordering story, including the measured fill
table, and an explicit note that `ddd.py` has had Markowitz all along.

---

# STAGE 8 — source models and `TLine`

**Work.** (a) `VPulse`/`IPulse` crash on their own class defaults — `per=0` (SPICE's one
shot) and `tr=0`/`tf=0` (SPICE substitutes TSTEP) both divide by zero, the latter because
`Pulse.f` evaluates all branches eagerly before `where()` selects. (b) `VSin`/`ISin`
ignore `td` and produce a *growing* waveform before the source starts — measured **2835 V
from a 1 V source**. (c) `Sin.next_event` — already done in 4g(a). (d) `TLine`: the DC
stamp is chosen by `len(self.history) == 0`, so the DC answer depends on whether a
transient ran first (`v(b) = 1` before, `0` after); history is never reset between runs;
there is no `next_event` and no `max_step` cap at TD, so the delay is measured as **4x too
long at dt=2e-9 and 10x at dt=5e-9**, silently.

**Gate 8-1:** `VPulse()` and `VSin()` with class defaults run to completion.
OUTCOME:
**Gate 8-2:** `VSin(td=...)` holds at the offset for `t < td`.
OUTCOME:
**Gate 8-3:** a `TLine` DC solve gives the same answer before and after a transient, and
`max_step <= TD/2` is enforced per line.
OUTCOME:

---

# STAGE 9 — consolidate `jaxtransient`

Blocked on 0.1a. The Gear2 LTE defect had to be found and fixed **twice** because this
file is a second transcription rather than a backend; the open JAX tolerance defect is the
same shape.

**Work.** (a) A `_lte_kernels.py` of backend-agnostic scalar algebra — the LTE formulas
are pure elementwise arithmetic on `(g_n, g_nm1, g_nm2, h1, h2)` with no control flow, so
they are *already* traceable. This also removes **five copies of the same three VSS alpha
coefficients**. (b) Give `JAXTransient` a `parameters` list: it currently declares none, so
`JAXTransient(cir, reltol=1e-6)` **raises `KeyError`** and there is no supported way to set
a tolerance; `eval_method='gear'` is hard-coded at both call sites; `nrsolver` and `scaler`
are accepted and silently ignored. (c) Thread tolerances to the LTE call sites — but as
**separate** `reltol`/`iabstol`/`vabstol`, or this re-creates the residual-vs-solution
flavour bug already fixed on the CPU. (d) Fix `for elem in self.cir.elements` (a dict, so
it yields string keys and `hasattr` is always False) — **0 breakpoints are found, always**
— together with the enumeration loop that is an *infinite loop* for pulse sources, since
fixing either alone converts a silent wrong answer into a hang. (e) Make the JAX Newton
report non-convergence instead of committing the unconverged point and computing its LTE
from it.

**Gate 9-1:** the three CPU gates from `transient_repair_plan.md` stage 4 ported to JAX —
LTE scales with the right power of h, a step is actually rejected, step count and error
respond to `reltol`. **None of these is currently expressible**, which is the asymmetry
that let the copied LTE survive.
OUTCOME:
**Gate 9-2:** a CPU/JAX agreement test in the suite. The stage-5 cross-check was run by
hand and written into prose; it is not in the suite, so the next divergence is invisible.
OUTCOME:
**Gate 9-3:** a `VPulse` transient under JAX hits the pulse edges.
OUTCOME:

---

# STAGE 10 — missing analyses

Scope per decision 0.3c. Ranked by the review's assessment of value:

1. DC sweep (`.dc`) — the most conspicuous absence; no `DCSweep` class exists.
2. `.tran` output control — `timestep` is currently *both* the initial step and `max_step`,
   and output points are the solver's own non-uniform steps, so every FFT needs hand
   resampling. `nonlinear_leapfrog_sweep.py` does exactly this, and interpolating
   non-uniform data before an FFT is a real accuracy hazard the user is forced to own.
3. `.ic` / `.nodeset` — `uic=True` currently means "start from zeros", not SPICE's
   per-element initial conditions; oscillators and latches are unstartable.
4. A SPICE-subset netlist reader — everything else in interop is downstream of getting a
   circuit *in*.
5. Large-signal MOSFET — no CMOS transient is expressible today.
6. Waveform export (raw/PSF/CSV).

**Also, one line:** `Transient` is not exported from `pycircuit/circuit/__init__.py`.

---

# STAGE 11 — `shooting.py`

Blocked on 0.1c. Repair, rewrite against the seams stage 7 creates, or withdraw.

---

## Order and dependencies

```
0.1a ────────────────────────────────────► 9
0.1b ──────────────────► 5
0.1c ──────────────────────────────────────────► 11
0.2a ──────────► 4
0.2c ─► (all suite gates)
0.3a ──► 1, 4
0.3b ──► 4

1 ─► 2 ─► 3 ─► 4 ─► 5 ─► 6          (7 in parallel throughout)
                              8, 10 after 6
```

**Stages 1-3 are the ones to do if only one thread is available**: they stop the silent
failures, give 10.5x, and make `reltol` mean something. Roughly a week including gates and
documentation.

Each stage commits with its explanation and its measured gate outcomes in the message.
Negative results are recorded in the same voice as positive ones.

**Do not `git push`.**
