# The transient analyses: the history kept out of the code

On 2026-09-27 (Andreas: "Move") the comments and docstrings of
`pycircuit/circuit/transient.py`, `jaxtransient.py` and `integrator.py`
were split the way `doc/shooting_history.md` split the shooting package's on
2026-09-24.  The code keeps what it does and what someone calling or
changing it must know now, in the present tense.  This document keeps the
rest, VERBATIM: the dated narratives, the measurements behind each
decision told as a story, the alternatives that were tried and rejected,
the defects and how they were found, and the original wording of every
paragraph that was condensed.

Each section is one module (and class), each subsection one member.  A
comment or docstring whose text moved says so with a line
`History: \`doc/transient_history.md\`, \`<Class.member>\``; search this file
for the member's name.  The text is the code's own, line for line, with the
comment markers and indentation removed -- read the "before" of a member
as its comment or docstring at commit `a2753239` or earlier.  Where a
sentence in the code had become WRONG about the current code it was
rewritten there, and its original is here.  The chronological record of
the work is `doc/pss_log_260902.md`.


## `transient.py` -- module level

### `MAX_GROWTH_RATIO`

2026-09-27 (moved from the code):

The comment on the `stepcontroller` import, before the move:

The clamp the step controller applies to every accepted step.  The force-accept
path in `solve()` is the one place that used to bypass it, and 4b's whole point
is that it must not: one bound, named once.  `stepcontroller` imports nothing
from this package, so this is import-safe at module level.

### `_threadpool_limits`

2026-09-27 (moved from the code):

The comment above the `threadpoolctl` import, before the move:

STAGE 2a -- BLAS thread control, discovered rather than required.

Circuit matrices are small (n ~ 10^2), so a threaded LAPACK spends more time
spawning and synchronising threads than doing the 1.7 MFLOP of work.  Measured
on the 139-unknown leapfrog: the whole transient runs **1.72x faster** with BLAS
limited to one thread, on a 4-core box.

`threadpoolctl` IS NOW A DEPENDENCY (2026-08-31, maintainer's decision).

It was optional, on the strength of a 1.72x measured on a 4-CORE box, and
this comment used to say that making it a hard dependency for 1.72x was the
maintainer's call.  The overhead scales with the core count the thread pool
spans, so that number was a property of that machine: re-measured on 24
cores the same comparison gives **14-20x** (minima 8.2 s against 145.1 s
over three interleaved pairs on the leapfrog, identical step counts; 14.2x
with a nonlinearity live).  The call was re-taken against 14-20x and went
the other way.

The import stays guarded anyway, and `blas_single_thread_available()` stays
in the API, so a stripped environment degrades to the status quo ante rather
than failing to import a circuit simulator.  What changed is that the
speedup no longer depends on someone remembering `OMP_NUM_THREADS=1` at the
shell -- which is how it was silently forfeited on this machine for weeks.

⚠ THIS LIMIT BELONGS TO THE TRANSIENT AND SHOULD NOT BE COPIED TO `DC` OR
`AC`.  Measured 2026-08-31 (work plan 2a-bis): the penalty is per-call
thread-pool overhead, crossing over around 50-100 solves, NOT a property of
problem size.  A transient runs thousands of small assemblies and solves in
a Python loop and so is destroyed by it; DC performs a handful of large
operations and PREFERS threads at every size measured (0.57x at n=28 down to
**0.43x at n=2503** -- limiting it would cost 2.3x), and AC never flips even
over a 500-point sweep (0.80x-0.98x).  Extending this context manager to
them is a measured regression, not an oversight.

### `resample_uniform`

2026-09-27 (moved from the code):

The docstring below its summary line, before the move:

STAGE 10.2.  A transient returns the solver's own adaptive points, so their
spacing is whatever the step controller chose -- measured at a **1000x**
spread on an ordinary RC driven by a sine.  Anything that needs a uniform
grid, above all an FFT, therefore has to resample, and every caller has been
left to do that themselves: `benchmarks/nonlinear_leapfrog_sweep.py` reads

    spec = np.fft.rfft(np.interp(grid, t, v)) / npt

-- `np.interp` is LINEAR, applied to a solution the integrator computed to
second order.  That throws away an order of accuracy before the transform,
which is a real hazard the caller was forced to own rather than a stylistic
preference.

QUADRATIC, to match the integrator.  Three-point Lagrange through the
interval's neighbours, the same reasoning `TLine._interpolate_history` gives
for its own history lookup: a first-order interpolant feeding a second-order
method injects error the method never made.  Measured on the solver's real
(non-uniform) grid, interpolating a known signal so only the interpolation
error is present:

    signal                  linear      quadratic    ratio
    fundamental           1.12e-03      8.11e-05     13.8x
    5th harmonic          1.55e-02      1.20e-02      1.3x
    decaying exponential  2.00e-04      4.90e-06     40.8x

**The 5th-harmonic row is the honest one**: where the adaptive grid is barely
resolving the signal, no interpolant recovers what was not sampled, and the
right fix there is a smaller `max_step`, not a cleverer resample.

Give exactly one of ``npoints``, ``step`` or ``grid``.

``grid`` takes an explicit set of output times, for a caller whose window
this function cannot otherwise express.  It was added 2026-08-31 for the
caller this docstring already named: an FFT that needs a TRAILING window
with ``endpoint=False``, because dropping the duplicated period boundary is
what keeps the tones exactly on bins.  Neither ``npoints`` (which spans
``t[0]..t[-1]`` inclusive) nor ``step`` (which starts at ``t[0]``) can say
that, which is the likely reason the named caller kept its own `np.interp`
for months after this function existed to replace it.


## `transient.py` -- `TransientStatistics`

### `class docstring`

2026-09-27 (moved from the code):

The docstring below its summary line, before the move:

STAGE 6(c).  Every number here was already being computed and thrown away --
`solve_system` returns its iteration count and the call site bound it to `_`;
the step controller knows what it rejected; the force-accept path from 4b
counts nothing.  A run that takes 40x more steps than expected is currently
indistinguishable, from the outside, from one that does not.

On the COUPLED path a persistently failing point raises BY DESIGN (F13):
its steps are solved, not rejected for error, so its only force-accept
is a step the excursion veto (`max_dv_step`) still refused after
`_CoupledSteps.max_reject` retries.  `rejected_steps` counts failed Newton
attempts too, on every path (one stepping loop since 2026-09-23; the LMM
loop alone used to drop them).

The force-accept counter is the one to read first.  It counts steps accepted
with an unbounded truncation error, and after 4d it should be zero on every
circuit measured -- so a non-zero value is the run telling you that part of
its own result is not error-controlled.

### `__init__`

2026-09-27 (moved from the code):

The comment on `branch_screens`, before the move:

`branch_check`: steps whose `rank C` fell below its structural
value (screened), and those a re-solve found a second root of.
⚠ Slots since 2026-09-27: without them `_branch_count` raised on
the first screen of every `solve`, the check caught its own error
and switched itself off -- it never reported on a full solve.


## `transient.py` -- `TransientStepError`

### `class docstring`

2026-09-27 (moved from the code):

The docstring, before the move:

A time point that could not be solved even at `minstep`, after the
continuation rescue.  Both a `NoConvergenceError` (what the Runge-Kutta
and coupled loops raised) and a `RuntimeError` (what the LMM loop raised)
-- one stepping loop since 2026-09-23, so a caller written against either
keeps catching it.


## `transient.py` -- module level

### `the step families`

2026-09-27 (moved from the code):

The comment above `_LMMSteps`, before the move:

THE STEP FAMILIES -- what differs between integrators in the one stepping
loop of `Transient._solve` (2026-09-23; it was three loops: the LMM loop,
`_run_rk_adaptive` and `_solve_coupled`, which had drifted apart -- the
excursion bound in two of them, the continuation rescue in two, three
Newton-failure ladders, three exception types, two setups).  The loop owns
everything else: breakpoints, `tend`, the failure ladder and rescue, the
excursion veto, the rejection budget and force-accept, state events, the
bookkeeping of an accepted step.  A family supplies:

  attempt(X, t, h, hold, provided_function) -> (x_new, h_taken, J)
      one step from X[-1] (`hold`: its size is imposed), or raises
      NoConvergenceError
  judge(X, x_new, h, J, clamped)  ->  (ok, h_next)  its error test
  h_after_force(h, h_next)  ->  the step after a force-accept
  next_breakpoint(t), after_accept(landing), finish()

and its retry budget at one time point: `max_retries` attempts in all, of
which `max_reject` may be rejections for error -- then the step is
force-accepted.


## `transient.py` -- `_LMMSteps`

### `__init__`

2026-09-27 (moved from the code):

The comment on `_step_controller_is_auto`, before the move:

Marked so the coupled family can tell this apart from a
controller the CALLER injected.  Without the distinction, any
object that ran an LMM run first presents this auto-created
controller to a coupled run, which then refuses a controller
nobody asked for -- 11 tests failed exactly that way.


## `transient.py` -- `_StageSteps`

### `class docstring`

2026-09-28 (the stage family honours `relref`; Andreas: "Do as you
recommend", on item 3 of the exact-diode audit's open items).

The class docstring before the change:

Runge-Kutta methods (radau, trbdf2, esdirk43) and Nordsieck GLMs:
self-starting, so no divided-difference LTE -- a step is judged by the
method's own embedded estimate, the filtered `_rk_est` the step leaves
(a GLM delivers `_glm_error_estimate` through the same slot).  Accept at
``err <= 1``; the next step is ``h * clamp(0.9 err^(-1/(p+1)), 0.5, 2)``
with ``p = EMBEDDED_ORDER`` (TR-BDF2 2(3) -> 1/3, Radau 5(3) -> 1/4), so
one anomalous estimate cannot swing the step wildly.  No history to
freeze at a corner, so more rejections are meaningful than for an LMM.

The weight was ``reltol |x_new| + abstol`` per unknown: pointwise, and
blind to `relref`, which the multistep family had honoured since D3
(default 'sigglobal', as in a commercial simulator).  How it was found:
TR-BDF2 looked wasteful on the hard-driven diode (20 V / 1 ohm, 1 kHz;
1530 attempts where esdirk43 took 509 and radau 364 for the same error).
Its rejection RATE was not the anomaly -- radau rejected 29 % there
against TR-BDF2's 27 % -- its step count was, and the rows that set the
step were the source's branch current: 406 of 414 rejections and 1004 of
1116 accepted steps, at the current's zero crossings, where the pointwise
tolerance collapses to `lte_iabstol`.  Adaptive GLM4 on the same circuit
crawled at h ~ 1e-20 s near t = 0 for the same reason.  Three better
controllers (no growth after a rejection; halving on a second; Gustafsson's
predictive) cut rejections 2-3x and saved no work: not the cause.

Measured on the real code (reltol 1e-6, attempts and max error against a
radau reltol 1e-12 fixed-step reference; 'pointlocal' / 'alllocal' /
'sigglobal'), the stateful Diode -- its state-free twin takes the same
steps, but for 635 in place of 633 in one cell:

    hard      trbdf2    1537 3e-5 |  368 3e-5 |  324 3e-5
    hard      radau      294 5e-5 |  179 5e-5 |  190 3e-5
    hard      esdirk43   488 3e-5 |  132 3e-5 |  130 6e-5
    hard      glm3      1105 3e-5 |  411 3e-5 |  388 3e-5
    hard      glm4       649 3e-5 |  383 3e-5 |  364 3e-5
    rectifier trbdf2     896 6e-5 |  733 6e-5 |  437 7e-5
    rectifier radau      324 2e-6 |  273 2e-6 |  198 2e-6
    rectifier esdirk43   207 2e-6 |  176 2e-5 |  143 4e-5
    rectifier glm3       986 2e-6 |  815 2e-6 |  544 2e-6
    rectifier glm4       633 2e-6 |  538 2e-6 |  399 2e-6

and at reltol 1e-4 the looser reference costs more where the pointwise
one had over-delivered: radau on the rectifier 1e-5 -> 4e-5, glm4 on the
rectifier 2e-6 -> 2e-4, glm4 on a smooth diode 1e-5 -> 2e-4.  Andreas
chose 'sigglobal' for every method knowing that price: `reltol` now means
the same thing whatever the integrator.

'pointlocal' is ``max(|x_new|, |x_n|)``, the multistep family's and
radau5's, not the old ``|x_new|``: within 1-5 % of the old step counts on
TR-BDF2, radau and esdirk43 where the old counts were measured (1537
against 1530 on the hard diode), but
GLM4 does not crawl under it (649 attempts) -- the crawl needed the
candidate's own zero.  The running maximum takes accepted points only,
where the multistep controller folds every judged candidate into it; the
counts moved by at most ~5 % against a candidate-folding prototype
(TR-BDF2 hard 325 -> 324, radau rectifier 197 -> 198, esdirk43 rectifier
136 -> 143, GLM4 hard 382 -> 364).
Tests: `test_the_stage_reference_follows_relref_over_accepted_points`,
`test_the_stage_methods_measure_their_error_against_relref`,
`test_a_stage_method_refuses_an_unknown_relref`,
`test_adaptive_glm4_does_not_crawl_at_a_source_current_zero` (each fails
on the old code; the last hangs there).


## `transient.py` -- `Transient`

### `_get_integrator`

2026-09-27 (moved from the code):

The comment on the Gear-2 default, before the move:

P6 OWNER DECISION (2026-08-21): Gear-2 is the shipped default,
matching the JAX backend at last -- identical scripts used to
get different methods by backend.  Chosen on the phase-0
measurement (half the steps and half the wall-clock of Euler
on the same circuit at the same tolerance) and because the
estimator work of stages 4g/4i was built for it; the
conformance harness pins the pair.

P22 RETIRED THE COUPLED EULER CARVE-OUT (same day): the
Gear-2 coupled livelock traced to eq (6) being measured on
ALGEBRAIC rows -- the rectifier's source-current row carries
the diode's dq/dt through KCL, its accepted value holds the
OLD grid's derivative convention, and the deviation floor
(2.5e-6 A against etol 3.6e-7) was h-independent.  The
state-row mask (_state_row_mask) retires that whole class;
with it, coupled+Gear-2 completes the rectifier in 259
points against Euler's 769 at the same accuracy (9.7e-3 vs
9.9e-3 against a fine reference).  The `coupled` parameter
the carve-out introduced is GONE -- the dead-knob scan
flagged it the moment it stopped being read.

### `parameters`

2026-09-27 (moved from the code):

The comment on `vabstol`, before the move:

NEWTON's x-tolerance on node rows, and nothing else.  Shared with `DC`,
which uses the same value, so the operating point and the steps after it
are solved to the same accuracy.

This used to be 1e-6 and used for BOTH roles.  The 1e-12 -> 1e-6 change
was reasoned about purely as a step-control knob (see `lte_vabstol`
below) and its effect on Newton was never measured -- it loosened node
convergence by 10^6 as a side effect, while `DC.vabstol` stayed at 1e-12,
so every transient was seeded by an operating point solved a million
times tighter than any step that followed it.  Decision 0.3a/0.3d in
`doc/transient_work_plan.md` split the two roles; this is the Newton one.
⚠ 1e-6 SINCE 2026-09-19 (Andreas), in DC, Transient, JAXTransient and
PSS TOGETHER -- the four share one meaning and one default.  It was
1e-12, below what double precision can deliver on a node once ANY
unknown in the circuit is large: a PLL shooting solve reached its
solution in two iterations and then failed the STEP test for 170 more
on a node at its zero crossing (`reltol |x|` gone, 1e-12 left)
against rounding of 1e-09 -- 142 s and "not converged" for an answer
that 5.7 s delivers at 1e-9 and at 1e-6 alike, digit for digit.
`lte_vabstol` is NOT this quantity and stays where it is.

The comment on `lte_vabstol`, before the move:

THE STEP CONTROLLER's tolerances, which are a different quantity: they
apply to `lte = J^-1 Eg`, not to Newton's residual or its x-update.

BACK TO 1e-12, matching `vabstol`, at gate D3-e.  The history is worth
keeping because it is a clean example of a workaround outliving its
defect.  This was 1e-12; it was raised to 1e-6 (a commercial simulator's `vabstol`,
SPICE's VNTOL) because on the 127-unknown leapfrog the timestep
collapsed to 5 ns against a 39 ns cap -- the controller accepts on
max(|lte|/etol) over ALL unknowns, most of that circuit's nodes carry no
signal, and under `pointlocal` the relative reference on a quiet node
tends to zero, so `etol` degenerated to TRTOL*abstol on numerical noise.
Raising the floor cut the step count 5.4x.

**That was treating the symptom.**  The cause was `pointlocal`, and
`sigglobal` -- the default since decision D3 -- references every unknown
to the largest signal in the circuit, so the reference cannot degenerate
and the floor is never reached.  Measured at gate D3-e, `lte_vabstol` at
1e-6, 1e-9 and 1e-12 give **bit-identical** results under `sigglobal`:
403 steps on a pulsed RC and 601 on a circuit with a quiet node, at every
one of the three values.  Under `pointlocal` the same change costs
+8.5% and +9.2% -- which is what "load-bearing" looks like, and why the
workaround was needed then and is not now.

So the principled value returns at measured zero cost.  **If you select
`relref='pointlocal'`, this floor becomes load-bearing again** and 1e-6
may be the better choice for your circuit.

The comment on `relref`, before the move:

What the RELATIVE part of the LTE tolerance is measured against --
A commercial simulator's parameter of the same name, and a commercial simulator's default.

`pointlocal` is what pycircuit did for its whole history: each unknown
referenced to itself, at this instant.  On a node carrying no signal
that reference tends to zero, so the tolerance collapses to the absolute
floor and the controller chases numerical noise on an idle node -- which
is the defect `lte_vabstol` was raised a millionfold to work around.
Measured on the leapfrog: at `lte_vabstol = 1e-12`, `pointlocal` needs
3.53x the steps of the shipped configuration where `sigglobal` needs
1.49x, i.e. it removes 81% of the excess.

DECISION D3, SECOND ATTEMPT -- `sigglobal` IS NOW THE DEFAULT, as in
A commercial simulator.  It was adopted, sent back by gate D3-a, and re-run.

What sent it back: referencing the tolerance to the largest signal lets
steps grow, and on an estimator carrying the trapezoidal `(-1)^n` mode
that was enough to break the controller's response to `reltol` --
accuracy stopped falling monotonically as the tolerance tightened.
4g(b) and 4i removed that contamination, and on re-run **all six
integrator/formula combinations are monotone in both step count and
error**.

What it buys, measured at MATCHED ACCURACY rather than at matched
`reltol` -- the two are not the same thing, and comparing step counts at
equal `reltol` overstates the win because `sigglobal`'s error at a given
`reltol` is ~1.5x larger:

  euler 1.48-2.06x fewer steps, gear2 1.44-1.60x, trapezoidal 1.31-1.47x

The figure previously recorded here was "1.7-2.5x fewer steps", taken at
equal `reltol`.  That was a relabelling of the tolerance as much as a
speedup; the honest worst case is 1.31x.

The comment on `lte_gamma_min`, before the move:

STAGE 12A -- Fang's acceptance band (DAC 2013 eq 15) and step-change
damper (eq 16).  The defaults are the historical one-sided test, so
nothing changes until a caller asks; see `StepController.set_lte_band`
for why the paper's own 0.7/3.0 are not adopted as defaults.
THE DEFAULT IS THE STRING 'auto', NOT A NUMBER, AND NOT None --
F5's lesson (doc/transient_review_260820.md).  Every documented
value is meaningful (0.0 disables the lower bound, 1.0 is the
historical threshold, None disables the damper), so a numeric or
None default cannot be told apart from an explicit request for it:
the coupled path used `par.lte_gamma_min or 0.7`, which silently
replaced an explicit 0.0 with 0.7 and made the documented settings
unreachable there.  'auto' resolves per path: the standard
controller maps it to (0.0, 1.0, None) inside set_lte_band, the
coupled path maps it to Fang's (0.7, 3.0, 0.15) in _coupled_band.
Any explicit value is honoured verbatim on both paths.

The comment on `coupled_method`, before the move:

STAGE 12B -- how the coupled path corrects the step size.

'approx'   Fang sec. 3.4: the new step comes from the error RATIO
           (eq 17) and the solution is corrected by eq (18).  The
           default, and the one with the measured record.
'bordered' Fang eq (12)/(14), RETIRED 2026-09-27: once its
           double-counted `q^T dv0` term was removed it took the
           same steps as 'approx' to every printed digit; asking
           for it raises, naming this.

The comment on `radau_transform`, before the move:

THE RADAU COST TRANSFORM, promoted 2026-09-27 (Andreas: "Promote
but measure first") from a private switch only tests set.  radau5's
eig(A^-1) split: one real and one complex m x m solve per
SIMPLIFIED-Newton iteration instead of the dense 3m system, and the
dense full Newton when it stalls.  Measured on a diode-loaded RC
ladder: fixed step 1.74x at m = 102 and 3.05x at m = 402, the answer
the same to the Newton tolerance (1.2e-5 at the defaults; the
dense-vs-transform test pins 1e-9 at tight ones); adaptive 1.5x
(m = 52) to 3.1x (m = 402), at most one fallback per run.  Off by
default: simplified Newton can stall on a strongly nonlinear step,
and the dense solve is the correctness reference.

The comment on `ic`, before the move:

STAGE 10.3 -- SPICE's `.ic`, for `uic=True`.

`uic=True` used to mean "start from a vector of zeros", which is not
what SPICE means by it and leaves a whole class of circuit
unstartable: an LC tank at zero is AT an equilibrium and stays there,
and a latch at zero sits on its metastable point.  Neither can be
simulated at all without a way to say where it starts.

Node voltages only, and that is a scope decision rather than an
oversight -- see `_initial_state` for what is deferred and why.

The comment on the step cap, before the move:

OWNER DECISION (2026-08-21): the cap is DECOUPLED from
`timestep` and renamed.  `timestep` used to double as the largest
accepted step, which silently made the step count on gentle
circuits a property of the requested output density rather than
of the error control -- the phase-0 "matched order" step
agreement (211 vs 210) was partly both backends sitting on the
same cap.  `timestep` now only sets the opening-step scale and
the fixed_timestep grid.

### `__init__`

2026-09-27 (moved from the code):

The comment above `__init__`, before the move:

`irefnode` was accepted here and never read -- the same shape as the dropped-
`toolkit` defect recorded below, found by the dead-argument scan
(doc/transient_review_260820.md, F18).  Passing it now fails loudly.

The comment on the `super().__init__` call, before the move:

(The class attribute already includes Analysis.parameters; the old
re-concatenation here double-included the base list -- hygiene.)
`toolkit` was accepted and then DROPPED -- it was never forwarded, so
`Transient(cir, toolkit=X)` silently ran on `cir.toolkit` instead.  It
went unnoticed because callers pass the toolkit the circuit already has,
which makes the two agree by coincidence rather than by construction.

The comment on the PCNR bookkeeping, before the move:

⚠ PCNR BOOKKEEPING IS INITIALISED HERE, NOT ONLY IN `_solve`.
`_solve` resets these per analysis, which is right -- but SHOOTING
never calls `_solve`: it drives `solve_timestep` directly on its own
grid, so a transient built for it reached the PCNR paths with the
attributes absent.  The stage-method fallback tripped over it
immediately (`AttributeError: 'Transient' object has no attribute
'pcnr_solves'`); the LMM path carried the same latent bug and had
simply never been reached that way.  Defining them at construction
makes every entry point safe, and `_solve`'s reset still gives each
analysis a clean count.

The comment on `irefnode`, before the move:

THE REFERENCE NODE IS `self.irefnode`, SET HERE AND NOWHERE ELSE
SPONTANEOUSLY -- `solve()` overwrites it from its `refnode`
argument.  It used to be set only inside `_solve`, which is the
hidden-ordering hazard this file itself documents at
`_apply_voltage_ics`; and `solve_timestep`'s own `refnode=gnd`
default meant the PCNR path pinned gnd while step control stripped
the caller's row -- two reference nodes in one solve
(doc/transient_review_260820.md, F7; measured: pcnr=True with
refnode='b' held gnd at 0 and let 'b' swing 4 V).  One fact, one
home.

### `_honours_continuation_rescue`

2026-09-27 (moved from the code):

The last paragraph of the docstring, before the move:

So this is `True` today for everything.  `_rescue_step` consults it
(since 2026-09-27; until then nothing did), because the honest-
diagnostic bug it was written for is easy to reintroduce: `_solve` used to report that a continuation "could not
rescue the point" on a path where none had been attempted.  A new step
path that reaches neither a ladder nor a fallback should return `False`
here rather than inherit a message that claims a rescue it never tried.

### `BRANCH_SCREEN_TOL`

2026-09-27 (moved from the code):

The paragraph on the negative-eigenvalue screen, in the comment above `BRANCH_SCREEN_TOL`, before the move:

⚠⚠ THE OBVIOUS SCREEN DOES NOT WORK HERE, and it was measured before it
was rejected.  A peer session proposed firing on the STEP MATRIX
ACQUIRING A NEGATIVE EIGENVALUE -- `C/h + G < 0` needs `C` small AND `G`
negative, so one number carries both ingredients -- and verified it 6/6
against root counts on scalar and 2x2 systems.  On real MNA it fires on
EVERYTHING: measured, min Re eig(J) is negative on the index-2 C-V loop
(-9.99e-07), on the exponential fixture (-9.90e-07) and on a van der Pol
(-1.00e+00).  Two independent reasons, both structural: MNA WITH A
VOLTAGE SOURCE IS A SADDLE-POINT SYSTEM and is indefinite by
construction, and an OSCILLATOR'S `G` HAS A NEGATIVE EIGENVALUE BY
DESIGN with no rank drop anywhere.  4 false fires out of 4 ordinary
circuits, so it cannot gate anything.

### `_branch_screen`

2026-09-27 (moved from the code):

The comment on the reference scale, before the move:

⚠ THE SCALE IS THE LARGEST `C` THIS RUN HAS SEEN, not the random
states' (`scale0`, until 2026-09-27).  Those sit anywhere in +-1 V,
where a forward junction's diffusion capacitance is e^38 above any
state the circuit visits: measured 0.217 F on a rectifier, so a
real 1.7e-10 F junction read as ZERO against 2e-6 F of load, the
screen fired on 2086 of 2500 steps and the confirmation re-solved
every one (3.8-4.9x the run; hidden while the check disabled
itself, see `TransientStatistics`).  The RANK stays structural --
a running maximum of the rank reads quiet where `C` is zero for the
whole run -- and for a linear `C` the two scales are one number.

### `_branch_confirm`

2026-09-27 (moved from the code):

The comment on the root test, before the move:

⚠⚠ AND IT MUST ACTUALLY BE A ROOT.  A solver that hands the
SEED BACK -- converged at iteration zero, or bailed -- looks
exactly like a second solution to a pure distance test, and
that is a FALSE ALARM on a default-on diagnostic.  Measured
on the coupled path before this check existed: the
"alternative" was 0.7071067811865475 in every component,
which is precisely the perturbed seed.

### `_branch_after_solve`

2026-09-27 (moved from the code):

The comment on the limiting-state snapshot, before the move:

⚠⚠ A DIAGNOSTIC THAT CHANGES THE SIMULATION IS A DEFECT.  Every
speculative solve writes the devices' limiting state (`Diode`'s
`_vlim`, which its `i` and `G` read), so the NEXT step would
start from the alternative's linearisation -- measured once:
`_vlim` 0.0 with the check off, 0.1017 with it on.  It goes
back EXACTLY, from a snapshot; the `limit(x, x)` re-sync this
used to do clamps against the stored state and lands short
above a junction's critical voltage (`state_snapshot`).

The comment on the exception handler, before the move:

⚠⚠ A DIAGNOSTIC MUST NOT BE ABLE TO FAIL A SOLVE -- it runs
after the answer is in hand and only reports.  BUT A BARE
`pass` HERE HID ITS OWN FIRST BUG: `self.statistics` does not
exist on a hand-driven march, the AttributeError was swallowed,
and the whole check silently did nothing while looking healthy.
So the failure is recorded and announced ONCE.

### `_branch_after_coupled`

2026-09-27 (moved from the code):

The comment on "fired" and the direction, before the move:

⚠ "FIRED" AND "GAVE A DIRECTION" ARE DIFFERENT ANSWERS, and
conflating them is why the first version of this read ZERO
screens on the very fixture it was written for: when `C`
collapses ENTIRELY the screen returns `(True, None)` -- there is
no null direction to name because every direction is one -- and
a `fired = direction` idiom then reads it as "did not fire".

The comment on the gauge shift, before the move:

⚠⚠ THE PERTURBATION MUST NOT BE A GAUGE SHIFT.  These are
FULL-WIDTH stage vectors, so a direction of `ones` moves the
REFERENCE NODE too -- a common-mode shift the circuit cannot
see.  The solve leaves the pinned row alone and hands it back
unchanged, the "alternative" differs from the base only in that
row, and its residual is EXACTLY ZERO because it is the same
physical solution.  Measured: `alt` came back as
[0.7071, 0.7071] with `r_alt = 0.0`, and that was 5 false alarms
out of 5 on the attracting control -- which the residual check
could not catch, the residual being genuinely zero.

### `_newton`

2026-09-27 (moved from the code):

The comment on the narrow `except`, before the move:

NARROW, deliberately.  This used to be `except Exception`, which turned
every failure inside a device model into "the circuit did not converge" --
a `ZeroDivisionError` from a source with `tr=0`, a `TypeError` from a
mis-specified parameter, an `AttributeError` from a typo in a subclass.
All of them were reported as a convergence failure, which sends the reader
to look at the bias point of a circuit whose real problem is a bug three
frames down.  The solvers already classify what they mean, so only their
own exceptions and genuine linear-algebra failures are translated here;
everything else propagates with its original type and traceback.

The comment on promoting to `SingularMatrix`, before the move:

The solvers wrap a singular factorisation as NoConvergenceError
("Singular Jacobian: ..."); promote that to SingularMatrix so callers
can tell "no solution here" from "could not get there".  Matching on
the message is weak and stage 6 replaces it with a real classification
off the zero pivot -- but it is what the previous code did, and
changing the taxonomy is not this stage's job.

The comment on the iteration count, before the move:

Stage 6(c): this count was bound to `_` and discarded.

### `get_diff`

2026-09-27 (moved from the code):

The comment above `get_diff`, before the move:

`method` was accepted and never read (doc/transient_review_260820.md, F18);
the integrator is selected by the `integrator` Parameter, not per call.

### `_opening_step`

2026-09-27 (moved from the code):

The docstring below its summary line, before the move:

STAGE 3.  A run used to open at `timestep`, which is also `max_step` -- the
largest step the controller is ever allowed to take.  The step controller
accepts the first step unevaluated, because with no history there is nothing
to difference and no truncation error can be estimated, so that opening step
was both the **largest** and the **only unchecked** step in the run.

Its error then dominated everything after it.  Measured on an RC step
response against the analytic solution, trapezoidal, before this method
existed: the global error was **1.3212e-01 at reltol 1e-3, 1e-4, 1e-5 AND
1e-6** -- identical to five digits -- while the step count went from 24 to
195.  Eight times the work for the same answer.  Backward Euler and the
two second-order methods also agreed to five digits, which is why the
integrator choice appeared not to matter.

Opening at `timestep * 1e-3` costs one cheap step and leaves the controller
to grow the step from there, which it does geometrically, so the ramp is
paid off within a handful of steps.

The principled alternative is a Hairer-style estimate from `q'`/`q''` at the
operating point.  The plan asks for the ramp first and a *measurement* of
whether the estimate is worth the complexity, rather than an assumption --
see the outcome recorded under gate 3-1.

### `_begin_run`

2026-09-27 (moved from the code):

The first paragraph of the docstring, before the move:

The standard and coupled loops (one loop since 2026-09-23) each
carried a byte-identical copy of this, and `PSS` needs it too -- it re-integrates one period from a
fresh state on every shooting iteration, so "begin a run" happens
many times per analysis there.

### `_pred_promote`

2026-09-27 (moved from the code):

The docstring paragraph on the accept site, before the move:

⚠ Called from the ACCEPT site (through :meth:`_push_history`) and
nowhere else, so a REJECTED step's stages never become nodes -- they
are samples of a trajectory the run then threw away.  (The stage
methods' adaptive driver used to be a second accept site that skipped
:meth:`_push_history`; one stepping loop since 2026-09-23, and every
family accepts through it -- a stage method reads no charge rings, but
its idtmod rows need the periodic shifts.)

### `_pred_note`

2026-09-27 (moved from the code):

The comment on copying the nodes, before the move:

⚠⚠ COPY, DO NOT VIEW.  `np.asarray` on an array that is already
float64 returns THE SAME OBJECT, so these entries would alias the
live state and stage vectors -- and the periodic gauge shift
subtracts `n*modulus` from every live history it knows about, so an
aliased entry takes the shift TWICE and the accepted state is
corrupted.  Measured as a real failure of
`test_a_state_fold_breaks_the_period_map_at_the_ENDPOINT_not_on_the_grid`,
and it failed with the predictor switched OFF -- the bookkeeping
runs either way, so 'off' is a control for the SEED, not for this.

### `_predict_state`

2026-09-27 (moved from the code):

The docstring paragraph on the clamp, before the move:

⚠⚠ THE CLAMP IS NOT OPTIONAL, and it is what the two rejected
predictors of the GLM measurement lacked.  A polynomial continued past
its last node can leave the region the circuit actually visits, and on
an exponential device a 3x overshoot is ``exp(3 dV / VT)``: measured,
the worst stage then cost 33 Newton iterations against the old seed's
6, while the MEAN still improved -- which is how such a heuristic
passes its own gate and fails in use.  The prediction is therefore
confined COMPONENTWISE to the range its own nodes span, widened by the
motion between the two newest.  The scale comes from the data.

The comment on the self-validation gate, before the move:

⚠ A SELF-VALIDATION GATE WAS BUILT HERE AND REMOVED, because it
changed nothing measurable: refit at the same degree from the nodes
one older, retrodict the newest node, and decline when the miss
exceeded a fraction of the step's motion.  Swept over thresholds
from 0.1 to infinity it moved ONE reading by 0.8% and every other
by nothing -- because it looks BACKWARD, and the case it was built
for is a knee that has not happened yet.  What actually bounds the
damage is the clamp below.

The comment on wrapping states, before the move:

⚠⚠ A WRAPPING STATE IS NOT A TRAJECTORY THIS CAN FIT.  A periodic
row folds by its modulus, and that is a DISCONTINUITY in exactly the
curve a polynomial is being put through -- the gauge shift keeps the
recorded nodes in one gauge, but the fold can also fall between the
newest node and the target, and then the fit runs straight across
it.  MEASURED as a real failure: on `Idtmod` with the wrap landing
exactly ON a grid point, where the period map is genuinely
discontinuous and the PSS is supposed to converge anyway
(`test_a_state_reset_needs_no_saltation_but_grid_alignment_is_a_cliff`),
predicting these rows stopped the shooting Newton converging at all.
Those rows keep the old seed -- the newest node's value -- and every
other row still gets the prediction.

### `_push_history`

2026-09-27 (moved from the code):

The first paragraph of the docstring, before the move:

Both accept paths in this file carried their own copy of these two
lines, and a third caller is arriving (`PSS`, which imposes its own
grid and so cannot use either loop).  Three transcriptions of a ring
push is how the trailing-window gauge shift below gets applied to two
of them and forgotten in the third.

### `_initial_state`

2026-09-27 (moved from the code):

The docstring paragraphs on `uic` and on element initial conditions, before the move:

STAGE 10.3.  `uic=True` previously meant a vector of zeros, which is not
SPICE's meaning and makes a class of circuit unsimulable rather than
merely inconvenient: **an LC tank at zero is at an equilibrium** and will
sit there forever, and a latch at zero is on its metastable point. There
was no way to start either.

**Node voltages only.** Element-level initial conditions -- SPICE's
``C ... IC=v`` and ``L ... IC=i`` -- are NOT implemented here, and the two
are deferred for different reasons:

* ``L``'s is a branch current, and its unknown exists in the MNA vector,
  so it needs only a reliable element-to-branch-index mapping. That is
  mechanical but not free: `SubCircuit.branches` is a flattened list and
  the element that owns each entry is not recorded.
* ``C``'s is a branch *voltage*, which constrains a DIFFERENCE of two
  node unknowns rather than either of them. A set of such constraints is
  a spanning-tree problem, not an assignment, and a floating capacitor
  chain has no unique node-voltage solution without one.

**Reconsider if** a circuit needs a floating capacitor's initial voltage
or a nonzero starting inductor current -- both are real requirements that
this does not cover, and neither is expressible by naming node voltages.

### `_apply_element_ics`

2026-09-27 (moved from the code):

The comment on the `state_ic` gate, before the move:

GATED ON `state_ic`, NOT ON A PARAMETER SPELLED `ic`.  This
test used to read `iparv.ic` first and then never use the
value on this path -- a NAME check.  A generated model whose
state seed is called anything else (`x0`, `phi0`) fell through
it, and `uic=True` started that state at zero while the DC pin
was perfectly correct: no error, no warning, a wrong waveform.
`IC_KIND == 'state'` and `state_ic` are installed under the
same condition (hdl.py, `state_meta['dc_pins']`), so this is
the same question asked of the thing that answers it.

### `_descendant_has_ic`

2026-09-27 (moved from the code):

The comment in the loop, before the move:

Same correction as `_apply_element_ics`: a state element
declares its seed through `state_ic`, not through a parameter
named `ic`.  Asking the old question here made the guard
UNDER-detect, which is the direction that silently drops an
initial condition instead of refusing it.

### `_solve_operating_point`

2026-09-27 (moved from the code):

The docstring below its summary line, before the move:

A failure here **raises**.  It used to substitute a vector of zeros, so a
circuit that had no operating point at all -- or one whose solve hit a bug in
a device model -- returned a complete, plausible-looking waveform computed
from a bias point that was never found.  Nothing in the result distinguished
that from a successful run, which makes it the most expensive class of defect
this module had: it does not fail, it lies.

The inner `DC` is constructed from *this* analysis's configuration rather
than from `DC`'s defaults.  Before, `DC(self.cir)` inherited none of the
transient's toolkit, environment parameters, tolerances, solver or scaler, so
the operating point could be solved at a different temperature, to a
different accuracy, and with a different Newton strategy than every step that
followed it -- and the mismatch was invisible.

### `LTERATIO`

2026-09-27 (moved from the code):

The comment above `LTERATIO`, before the move:

The LTE tolerance multiplier.  `TRTOL` in this module, `lteratio` in
A commercial simulator: the LTE estimate is deliberately conservative, so the allowed
truncation error is this many times the Newton-solve tolerance.
P2 (doc/backend_parity_260821.md): settable at last -- the asymmetry
was REVERSED, JAX declaring `TRTOL` as a Parameter while this side
hardcoded a class constant, so a user tuning a commercial simulator's `lteratio`
could do it on one backend only.  A property rather than the old
class attribute, so every existing `self.LTERATIO` read follows the
Parameter and the two cannot drift.

### `_newton_xtol_vector`

2026-09-27 (moved from the code):

The last paragraph of the docstring, before the move:

An increment on a node is a voltage and on a branch is a current, so the
two vectors are transposed with respect to each other.  Getting this
backwards is the same class of error stage 0.3d separated for the LTE
tolerances: the numbers are dimensionally different quantities that
happened to share a default (both 1e-12 until 2026-09-19, when
`vabstol` became 1e-6 -- which is what makes a swap visible).

### `_companion_at`

2026-09-27 (moved from the code):

The docstring, before the move:

``(iq, Geq)``: the step's companion current and conductance at `x`
(the current ``self._dt``), with the charge cached against the state
it belongs to.  One assembly for every step's residual and Jacobian
(2026-09-24: it was six copies).

`self.epar`, not the module-level `defaultepar`.  Omitting it meant
every device in a transient was evaluated at defaultepar's T = 300 K
whatever the caller asked for, and -- because `Analysis.__init__`
attaches `bypasstol` to the analysis's own epar and nowhere else --
every device took its `except AttributeError` branch and the `bypass`
parameter did nothing at all.

STAGE 2c.  The charge vector is stashed alongside the state it belongs
to.  `solve()` needs `q` at the converged point twice more -- once for
the step controller and once for the history roll -- and was
recomputing the whole assembly both times at an x it had already
evaluated.  Measured 5.08 `q` assemblies per accepted step against
3.06 for every other stamp; the difference is exactly those two.
Keyed by the state so a stale value can never be served: the check is
identity-then-equality on x, not a bare "did we cache".

### `_source_at`

2026-09-27 (moved from the code):

The docstring paragraph on the contract, before the move:

ONE CONTRACT: `provided_function(t)` is an extra source term, on every
path.  The standard path used to treat it as a post-solve callback
`provided_function(f, J, C)` whose result was unpacked and never read,
while the coupled and PCNR paths added it to `u` -- two contradictory
meanings behind one parameter, flag-selected
(doc/transient_review_260820.md, F4).  The callback contract was born
dead: its introducing commit says "currently is calculated but returns
no value to solve method", and no consumer ever appeared.  The live
semantics wins; callback callers break loudly on arity.

### `_fang_timestep_inner`

2026-09-27 (moved from the code):

The paragraph on refusing an injected controller, in the comment on GATE 12-4, before the move:

So an injected controller is REFUSED unless it is one whose law this
path actually implements.  Accepting it and using it only for
`relref` would make `tran.step_controller = IntegralController()` look
honoured while doing nothing -- the same class of defect as a
documented feature that does not exist, which is what this path was
until now (it silently built its own and ignored the caller's).

The comment on PCNR, before the move:

PCNR ON THE COUPLED PATH.  `pcnr=True` was SILENTLY IGNORED here:
`_solve` dispatches on `coupled_lte` before it ever looks at `pcnr`,
so the run took the classic limiter and the parameter did nothing --
measured as 0 PCNR steps against 4869 `Diode.limit` calls, and results
bit-identical to `pcnr=False`.

The comment on `hold_h`, before the move:

`hold_h` -- the step size is IMPOSED, not free.  A step truncated
onto a breakpoint or onto `tend` has its size decided by where it
must land, so there is nothing for the coupled system to SOLVE.

Without it the truncation was pointless -- `fang_timestep` solved
for its own `h` and walked straight off the edge again: 0 of 10
pulse edges landed on, worst miss 1.24e-7 s, the whole rise time.

BUT "DO NOT SOLVE FOR h" IS NOT "DO NOT CHECK THE ERROR", and
conflating the two was a defect worth the same scrutiny as the one
it replaced.  A held step was accepted blind, so its truncation
error was governed by nothing: on the pulsed RC the maximum error
sat at 1.465e-2 at BOTH reltol 1e-5 and 1e-6 -- identical across a
decade of tolerance, the signature of a quantity no tolerance
controls -- and the mean was 5.9x the standard path's at 1e-6,
getting worse as the tolerance tightened.

A held step whose error is over the band is reported so the
caller can shrink and retry -- UNLESS the grid is locked.

`fixed_timestep` is the caller stating that the output points are
theirs, so shrinking is not an option available to us: the honest
response to an over-tolerance step on a locked grid is to take it
and let the run's accuracy be what the caller asked for, exactly
as the standard path does. Conflating "truncated onto a
breakpoint" with "grid imposed by the caller" broke
`test_fixed_timestep_keeps_the_grid_on_the_coupled_path`: the
retry shrank `h` and the uniform grid disappeared.

The paragraph on eq (12), in the comment on a failed LTE condition, before the move:

SEC. 3.4's APPROXIMATE NEWTON, NOT EQ (12), AND THE REASON IS
MEASURED.  Eq (12) recovers `dh` from eq (14), whose denominator
is `q^T dxh + d`.  Those two terms are the solution's sensitivity
to the step size and the extrapolation's slope, and BOTH are
approximately `dv/dt`: their difference is the truncation error's
derivative, which is tiny by construction.  Measured on a driven
RC at h = 1.6e-7: `q^T dxh = +1.818e9`, `d = -1.820e9`, denominator
-2e6.  Three digits lost, and the SIGN of the denominator decided
by the cancellation -- so `dh` saturated at the eta limit with an
essentially arbitrary sign and the step drifted down four decades
while `err` sat at 0.2, far BELOW the band that should have grown
it.  Eq (12) computes a small quantity as the difference of two
large ones; this is very likely what sec. 3.4 means by "the
coupled nonlinear system sometimes is very sensitive to the change
of step size".

The comment on the retired bordered method, before the move:

(`coupled_method='bordered'`, Fang eq (12)/(14) -- a Newton step
on the LTE equation with an analytic denominator -- stood here
until 2026-09-27; retired, see the check above.  History:
doc/pss_log_260902.md, 2026-09-27, and git.)

The second paragraph of the comment on `h_want`, before the move:

That hole is what let the first step after a pulse edge be
accepted at 2.0e-7 s when it needed 3.55e-9 s, a factor of 56.
The step came out of `fang_timestep` at exactly 0.2x its entry
value -- MIN_SHRINK_RATIO, the within-time-point floor -- after 12
iterations, reporting converged. It produced v = 0.033333 against
an analytic 0.018731, a 78% single-step error, and the resulting
1.465e-2 was IDENTICAL at reltol 1e-5 and 1e-6 because nothing
about it was tolerance-controlled.

The measured paragraph of the comment on saturation, before the move:

Measured on the pulsed RC: a step that needed to shrink tenfold
just after a rising edge declared itself converged after a single
15% shrink, ran at h = 4.0e-8 s where the standard path used
4.4e-9 s, and left a maximum error of 1.465e-2 that was IDENTICAL
at reltol 1e-5 and 1e-6 -- the signature of an error no tolerance
governs.

### `_lte_tolerance`

2026-09-27 (moved from the code):

The comment on the state-row mask, before the move:

P22: eq (6) over the STATE rows only -- an infinite tolerance on
algebraic rows removes them from the band test, the controlling-
node argmax through this one mechanism (and, until its retirement on
2026-09-27, the 'bordered' branch's lte_gradients).  See
_state_row_mask for the derivation and the measured livelock this
retires.

### `_excursion_ratio`

2026-09-27 (moved from the code):

The last paragraph of the docstring, before the move:

⚠ THE RUNGE-KUTTA LOOP NEVER HAD IT (fixed 2026-09-23, before the
three stepping loops became one).  The check was typed out in the LMM
loop and again in the coupled loop: every Runge-Kutta and GLM method
ignored both knobs SILENTLY -- on a pulsed RC with 'auto' (bound
0.098 V) gear and trap went 70 -> 86 steps and held 0.097 V while
radau stayed at 69 steps and 0.0996 V, trbdf2 at 71.

### `_pcnr_attempt`

2026-09-27 (moved from the code):

The docstring from its second sentence, before the move:

`pcnr_status` is 'used' / 'partial' / 'fell-back' from the two
counts.  (One place since 2026-09-27; the LMM, stage and coupled
paths each carried the four lines.)

⚠ ANY EXCEPTION FALLS BACK, on every path (Andreas, 2026-09-27), as
the LMM step and DC always did: a PCNR failure on one point must not
end the run.  The stage and coupled paths used to catch only a
non-convergence, and a singular Jacobian then ENDED the transient --
measured, Radau + PCNR on a FET cascode with a vanishing capacitor,
`LinAlgError: Singular matrix`, where device limiting ran.  Every
path's warning names the exception's type, so a genuine bug on the
PCNR path is still visible in the log.

### `_pcnr_newton`

2026-09-27 (moved from the code):

The last line of the docstring, before the move:

(One loop since 2026-09-27; each family had its own copy.)

The comment on both residuals, before the move:

BOTH residuals, not just the MNA one.

`solve_dc` checks `g_lim` and this path did not, which means it
could return with `v_lim != e_a - e_b` -- the diode evaluated at a
voltage that is not the node voltage, so the returned vector is
not a solution of the circuit at all.  Everything downstream then
inherits it: the charge history is wrong, and the LTE estimate
built from that history reads low, so the step controller takes
large steps believing they are accurate.  Measured on a half-wave
rectifier as a median accepted error of 0.0066 against a target of
0.81, while the actual waveform error was 2.5x worse than the
classic path's.

### `_solve_timestep_pcnr`

2026-09-27 (moved from the code):

The comment on the stage predictor, before the move:

STAGE PREDICTOR -- and it has to be here, not only on the limiting
path, or the two stop agreeing.  ⚠ MEASURED as a real failure of
`test_gate_13_6_pcnr_and_limiting_take_the_same_steps`: with the
predictor on one path only, the two converge to values that differ
in the last digits, that moves the LTE estimate, and the step
sequences part company at 5e-7 by the end of the run.  The gate is
right and the asymmetry was the defect.  It also seeds `v_lim`, and
limiting the seed is what fixed PCNR's one documented failure.

The comment on the Jacobian handed to the step controller, before the move:

THE JACOBIAN HANDED TO THE STEP CONTROLLER MUST BE THE ONE
THIS PATH ACTUALLY SOLVED, and `cir.G(x) + Geq` is not it.

`Diode.G` linearises around `_vlim`, which only `Diode.limit`
updates -- and PCNR never calls it, because limiting is the
thing PCNR replaces.  So `_vlim` stays at whatever it was
first set to and the diode's conductance is frozen there:
measured on a half-wave rectifier as `_vlim` stuck at 0.0 V
across 2283 `G` evaluations while the junction actually swung
-18.47 to 0.75 V.  `cir.G(x)` therefore carries NO diode
conductance at all.

Inside `augmented_system` that cancels -- the same wrong value
is added by `cir.G` and subtracted again -- but the controller
computes `lte = J^-1 Eg`, so handing it that matrix maps the
truncation error through a Jacobian missing the diode.

The right matrix is the one `predict` factorises: the non-PCNR
part plus each probe's `didv` column as a rank-one update.  At
convergence `v_lim == e_a - e_b`, so it is exactly the
Jacobian of the residual with respect to `x` -- and it is
`schur_reduce`'s matrix, taken from there rather than
written out a second time (the copy that used to live here
knew only the two-terminal `(dia, dib)` shape).

2026-09-28: the condensed comment on the controller's Jacobian, as the
2026-09-27 move left it, before it was rewritten -- `augmented_system`
excludes the PCNR devices from the ordinary assembly (the skip set), so
their stamp is no longer "added and subtracted", and the multistep step
now syncs the limiting state after its Newton (`limit_sync`):

`Diode.G` linearises around `_vlim`, which only `Diode.limit`
updates -- and PCNR never calls it, because limiting is the
thing PCNR replaces.  So `_vlim` stays at whatever it was
first set to and the diode's conductance is frozen there:
`cir.G(x)` carries NO diode conductance at all.

Inside `augmented_system` that cancels -- the same wrong value
is added by `cir.G` and subtracted again -- but the controller
computes `lte = J^-1 Eg`, so handing it that matrix maps the
truncation error through a Jacobian missing the diode.

### `_stage_source`

2026-09-27 (moved from the code):

The comment above `_stage_source`, before the move:

-- what every stage step shares (2026-09-23: the DIRK, GLM and the three
Radau steps each carried a copy of the source closure and of the
epilogue; the DIRK and GLM steps of the implicit stage solve too) ------

### `_solve_implicit_stage`

2026-09-27 (moved from the code):

The last sentence of the docstring, before the move:

(Under a GLM this is E1, 2026-09-16: before it, `pcnr=True` did device
limiting there and SAID it had used PCNR -- 3 PCNR solves against
39987 `Diode.limit` calls on a half-wave rectifier.)

The comment on the branch check, before the move:

THE BRANCH CHECK, CONFIRMED on the stage equation `_newton`
would have solved (until 2026-09-24 this path ran NONE:
neither the screen nor the confirmation)

### `_rk_step_dirk`

2026-09-27 (moved from the code):

The last paragraph of the docstring, before the move:

⚠ THIS IS THE GENERIC FORM OF THE OLD `_solve_timestep_trbdf2`.  For
TR-BDF2's ESDIRK tableau the two implicit stages share the diagonal
``d = STAGE_DIAG`` (the one-LU property), and the stage form is
algebraically identical to the old TR + BDF2-companion writing
(verified before the bespoke step was removed).  Any SDIRK/ESDIRK to
come reuses this untouched.

### `_glm_stage_predictor`

2026-09-27 (moved from the code):

The docstring paragraphs on the rejected predictors, before the move:

⚠⚠ TWO SIMPLER PREDICTORS WERE BUILT FIRST AND BOTH LOST, on a
state-free exponential at 40 and 200 points per period
(``benchmarks/stage_predictor.py``; device evaluations against the
old guess, then the worst seed error, then the worst stage's Newton
iterations):

==================  ==============  ============  ===========
predictor           device evals    worst seed    worst iters
==================  ==============  ============  ===========
old guess           --              7.5e-02       6
``Y_i^prev + dx``   +0.6% to -6%    2.0e-01       13
full-step poly      -2% to -18%     7.3e-01       33
this one            -12% to -29%    7.5e-02       5
==================  ==============  ============  ===========

A gate on the MEAN passes all four.  The mechanism is structural: the
old guess is always a value the circuit ACTUALLY ATTAINED, so it can
never sit in a device's overflow region, while a polynomial continued
a whole step can, and on an exponential a 3x overshoot is
``exp(3 dV / VT)``.  ⚠ A ratio test against the step's own motion does
NOT screen the bad case -- measured, the bad prediction's displacement
is 2.98 of that motion and the TRUE stage spread reaches 2.98 too.

### `_solve_timestep_glm`

2026-09-27 (moved from the code):

The comment on the two Nordsieck slots, before the move:

⚠⚠ TWO SLOTS, AND THE SECOND ONE IS WHAT MAKES ADAPTIVE STEPPING
POSSIBLE AT ALL.  `_glm_Q` holds the vector this method last
PRODUCED, valid at the end of that step; `_glm_Q_at_entry` holds the
one it last CONSUMED, valid at the start.  A step that the
controller REJECTS has already overwritten the first, and the retry
-- which starts from the same `x` at the same `tn`, only with a
smaller `h` -- then finds no vector valid at `tn` and runs the full
startup.  MEASURED before this existed: every single step rejected
once and every step paying a startup (2732 accepted, 2732 rejected,
2733 startups on GLM2), 102034 device evaluations against radau's
2419 on the same problem -- 42x, and a vicious cycle rather than a
slow path, because a fresh startup's top component makes the next
estimate spurious too.

### `_rk_step_coupled_pcnr`

2026-09-27 (moved from the code):

The comment on the residual current, before the move:

⚠ THE RESIDUAL CURRENT IS `g_mna`, NOT the Schur RHS
`f_eff`.  `augmented_system` STAMPS the junction current at
`v_lim` into `g_mna` (via `dev.stamp`), so `g_mna` is the
PHYSICAL MNA current with junctions at `v_lim` -- exactly
`cir.i` once `v_lim` == the branch voltage.  `f_eff =
g_mna - J_ml g_lim` folds the junction current into the
Newton-STEP right-hand side instead, where it VANISHES as
`g_lim -> 0`; using it as the residual dropped the junction
current at convergence and converged the coupled Newton to a
neighbouring, wrong root (measured: step-1 node error 8.5e-5,
true stage residual 4e-14 vs 8e-25 for device limiting).  Only
`G_eff = J_eff` (the Schur-reduced Jacobian) is taken here;
the junction is eliminated from the step by the correct phase
(`dx_lim_of` + `refine`) below.

The comment on the missing continuation ladder, before the move:

⚠ THERE IS DELIBERATELY NO CONTINUATION LADDER HERE.  The step
instead FALLS BACK to the device-limiting coupled solve (see the
caller, `_rk_step_coupled`), which carries one.  A ladder was
built for this path three times and removed each time; the
design and the reason are recorded so a fourth attempt starts
from the evidence rather than repeating it.

THE DESIGN, IF IT IS EVER NEEDED.  A CAPACITANCE ACROSS THE
LIMITED JUNCTION, anchored at the last accepted state -- the
two-node incidence stamp of `JunctionGminSteppingNewton` but
carrying `g (v_j - v_j,n)` instead of `g v_j`, i.e. the
backward-Euler form of `C = g h a_ii` in parallel with the
junction.  It rides `pcnr.augmented_system`'s existing
`u_extra`/`J_extra` hooks (`u_extra` the branch current,
`J_extra` the 2x2 pattern), so it needs no new plumbing, and the
Schur reduction carries it by construction.  Two rules that cost
measurements to learn: the anchor must go ACROSS THE JUNCTION and
not on every row -- `g * eye(n)` also anchors a voltage source's
BRANCH-CURRENT row, where `g (i - i_n)` is a conductance on a
current unknown; and its schedule must start ABOVE the circuit's
own conductance (`||G(x_n)||inf`), since a first attempt marched
`g <= 1 S` against a 10 mOhm source and never bit.

IT WORKS AS A MECHANISM: with the anchor present the strong rungs
converge in TWO iterations with `|g_lim| = 0`, and at the rung
where a whole-diagonal anchor let the junction gap snap back to
359 V the two-node stamp held it to 50 V.

⚠⚠ WHY IT IS NOT BUILT.  No circuit is known where PCNR fails at
a normal iteration budget.  PCNR's ONE documented failure -- the
BJT mirror of `test_dc_pcnr.py` from a uniform 20 V start -- was
fixed at its source by LIMITING THE SEED (`pcnr.v_lim_init`,
+20 V: LinAlgError -> 8 iterations), and that docstring already
says a continuation could never have fixed it: *"No ladder around
the solve could help, because every rung began by building the
same Jacobian at the same unlimited seed."*  Re-measured: that
mirror solves at DC with `pcnr_status='used'`, and driven as a
transient (pulsed 0->5 V, rise times to 1 ps, steps to 1 ps,
radau and trbdf2) every combination converges with 0 rungs.
Transient also suppresses the mode structurally -- PCNR fails on
a far-off INITIAL GUESS, and every step here starts from the last
accepted state.

SO THE TRIGGER TO WATCH FOR is a PCNR stage failure whose
`g_lim` stays large while the MNA state diverges (the signature:
`|g_lim|` ~ hundreds of volts, `ynorm` running to 1e17) on a
circuit at a DEFAULT `maxiter`.  ⚠ Do NOT validate a candidate by
starving `maxiter`: a ladder must end with a PURE solve of the
original system, so a starved budget defeats the final rung
whatever the deformation -- two earlier rungs were rejected on
evidence from exactly that broken instrument.

### `_coupled_stage_solver`

2026-09-27 (moved from the code):

The comment on the line search, in the nested `_stage_newton`, before the move:

⚠ BACKTRACKING LINE SEARCH, AS A RETRY (owner decision
2026-09-08, "Do 2").  This Newton had no damping at all, and a
whole class of circuits cannot be integrated undamped at ANY
grid: a nonlinearity fed by a branch current through a
capacitor sees the companion difference quotient, whose stage
sensitivity is tau/h, and the undamped basin is 0.94 h/(k tau)
(measured on a scalar model of one stage, READING-LOG 2.156) --
refining the grid SHRINKS it in proportion.  A comparator on a
current sense through a capacitor failed here at every grid
tried, on every method; with the retry its steps converge.
The rule and floor are `DampedNewton`'s (Armijo 1e-4, alpha >=
0.05).
⚠⚠ WHY A RETRY AND NOT A BLANKET SEARCH: a blanket line search
CHANGED CONVERGED ANSWERS.  At h = 100 tau on a tanh charge
circuit the coupled stage system has more than one root (its
a_23 < 0 makes the coupled equations non-monotone where the
scalar map is monotone); the undamped-plus-limited path found
the physical one, v in [0, 1], and the damped path a spurious
one at v = -0.010, and two shooting tests moved with it.  So
`damped=False` (the first attempt) is the old algorithm bit
for bit -- the full step always, the old convergence test on
it -- and the search is tried only after it has failed, before
the shunt ladder.  Every case that converged before converges
to the same numbers; the suite is the control.

The comment on the convergence test, in the nested `_stage_newton`, before the move:

⚠ THE OLD CONVERGENCE TEST, ON THE FULL STEP, BEFORE
THE LINE SEARCH SEES IT.  Two wrong versions preceded
this one, both caught by the identity control: testing
the size of a DAMPED step let a step at the floor stop
the iteration short of the root (a tanh charge circuit
at h = 100 tau returned v in [-0.010, 0.986] against
the undamped [0, 1]); requiring a full step to pass
Armijo hung at the roundoff floor, where both residual
norms are noise and the test fails by chance, so no
full step was ever accepted and every step burned
`maxit` iterations.  A converged full step is accepted
here exactly as the undamped Newton accepted it, so
every case that converged before converges to the same
numbers; the search engages only when the full step is
neither small nor residual-reducing.

### `_rk_step_coupled`

2026-09-27 (moved from the code):

The first paragraph of the comment on the continuation rescue, before the move:

⚠ THE CONTINUATION RESCUE REACHES THE COUPLED PATH TOO.  It is
armed by `_solve` once the step has shrunk to `minstep`, and that
arming used to do NOTHING here -- `_continuation_rescue` is read
inside `self._newton`, which the coupled `sm` solve does not go
through (measured: with the flag set over a 40-step run TR-BDF2
wrapped the rescue solver 80 times and Radau 0).  So the last
resort before `_solve` gives up was silently absent on every
fully-implicit method, and the failure then claimed a ladder had
been tried when none had.

The last paragraph of the comment on the continuation rescue, before the move:

⚠ THE LINE SEARCH IS THE LAST RESORT, AFTER THE LADDER (owner
decision 2026-09-08, "Do 2").  Tried as the FIRST retry it
pre-empted the ladder and beat it to a SPURIOUS root on a tanh
charge circuit at h = 100 tau (v = -0.010 where the ladder's
homotopy gives the physical [0, 1]); the ladder tracks a
physical branch and the search does not.  Where no ladder is
available -- the PSS path never arms one, it drives
`solve_timestep` on its own grid -- there is nothing to
pre-empt (an undamped failure used to raise straight out), so
the search is the one retry there, opt-in through
`_damped_last_resort`, which PSS sets; a plain adaptive
transient keeps its step sequence unchanged.

The comment on branch detection, before the move:

BRANCH DETECTION on the COUPLED path.  ⚠ This path does NOT go
through `self._newton`, so the check wired there reached the LMM,
DIRK-sequential and GLM stage paths and left the fully-implicit one
-- which is the PSS default -- unscreened.  Same screen, same
confirmation, but the solve being re-run is the coupled `3m`
`_stage_newton` rather than a single-stage residual, so it needs its
own call: a stage of the block can be at a rank drop while the
others are not, and it is the BLOCK that has to be re-solved.

### `solve_timestep`

2026-09-27 (moved from the code):

The comment on Runge-Kutta methods, before the move:

ANY Runge-Kutta method: the one tableau-driven stage step, which
picks the DIRK-sequential or fully-implicit-coupled path from the
tableau's structure.  PCNR (when `par.pcnr`) is the per-stage
limiting on the DIRK-sequential path -- see `_rk_stage_pcnr`; it
flows to shooting too, since the inner transient calls this same
method.  (The FULL coupled path still limits; PCNR-on-coupled is
a follow-up.)

The comment on PCNR participation, before the move:

Gate PARTICIPATION on the device records, not on the
pnj-only pair view: that view exists for the gmin ladders
and is empty for a circuit of pure fetlim/limvds devices,
so `pcnr=True` on a MOSFET differential pair used to fall
through to the ordinary solver SILENTLY (vector PCNR
Stage 2, 2026-08-26).

The comment on a circuit with no participants, before the move:

Asked for, and no device declares a probe.  Falling
through is right -- refusing would be a worse answer --
but it must SAY SO, or `pcnr=True` and `pcnr=False` are
indistinguishable from outside.  Same defect DC carried
until sec. 47.

The second paragraph of the docstring of the nested `jacobian_only`, before the move:

So `cir.i(x)` and `cir.u(t)` were assembled once per accepted step and
discarded.  `C`, `q` and `G` are all still needed -- `C` and `G` build
`J`, and `q` feeds the charge cache and the history roll -- so this is
not a cheaper approximation of the same work, it is the same work minus
two vectors that had no consumer.

### `solve`

2026-09-27 (moved from the code):

The comment above `solve`, before the move:

`analytical_eh` was accepted through all three signatures and read by nothing
(doc/transient_review_260820.md, F8) -- superseded by `coupled_method`,
which carries the measured record for both branches.  Deleted, so passing
it raises TypeError instead of being silently discarded.

### `_finish_result`

2026-09-27 (moved from the code):

The docstring, before the move:

The run's `CircuitResult` (one helper since 2026-09-23, when the
three stepping loops each had a copy -- then one loop).

⚠ The t=0 point IS part of the result -- SPICE convention, and what
the JAX backend already does.  `X[0]` is the operating point (or the
uic vector) the run worked to compute; dropping it made every
index-aligned backend comparison off by one and left
resample_uniform unable to reproduce the initial value
(doc/transient_review_260820.md, F12(a)).

STAGE 10.2 -- resampled onto a uniform grid if one was asked for,
after the run rather than inside it, deliberately: the adaptive grid
is what the error control is defined on, so the solver keeps
choosing its own steps and only the REPORTED points change.  The
statistics are reachable from the result, not only from the
analysis (Stage 6(c), F13), so a caller who kept only the waveform
can still ask what produced it.

### `_solve`

2026-09-27 (moved from the code):

The comment on refusing `coupled_lte` for stage methods, before the move:

⚠ FANG'S COUPLED PATH IS BUILT ON A LINEAR MULTISTEP COMPANION: it
solves the step from eq (6), a solution-space LTE over the step
history.  A stage method or a GLM judges its step by its own
embedded estimate; under `coupled_lte=True` they died inside
`compute_derivatives` with a message that never named the flag.
Refused here, by name, before any work.
⚠ AND KEPT REFUSED ON MEASUREMENT (2026-09-25, plan item 4): a
prototype of the loop around the stage step (eq (6) of the
method's order, eq (17)'s error-ratio step, the step re-solved)
against the standard adaptive stage run, on the stage-12
benchmarks against their closed forms.  It does not keep Fang's
no-rejection property -- its re-solves outnumber the standard
run's rejections -- and on the stiff RLC its error is set by no
tolerance (flat across two decades of reltol: radau 1.21e-3 where
the standard run reads 1.6e-7 .. 6.5e-10, trbdf2 2.33e-3 against
4e-4 .. 1.9e-5, glm2 9.7e-5 against 3.5e-5 .. 1.2e-5), because
eq (6) needs accepted history and so cannot judge the opening
steps an embedded estimate judges from the first.  On the driven
RC it saves steps by integrating less accurately (glm2: 3x fewer
steps, 8.5x the error).  Script:
`benchmarks/transient_review/stage12c_fang_stage_methods.py`.

The comment on resetting element state, before the move:

Position matters and cost a test to learn: placed after the initial
`accept_step(0.0, ...)` this wiped the very history that call had just
seeded, so `TLine.G` saw an empty buffer and stamped the line as a DC
SHORT -- v(p1) came out 1.0 where 0.5 is correct.  Elements that carry
state must be reset before the run seeds them, not after.

The comment on the step cap, before the move:

DECISION D2, 2026-08-01.  The clamp on how large an ACCEPTED step may
grow.  It defaults to `timestep` -- the historical behaviour, and what
`.tran tstep` means to most callers -- but it is now reachable.

Measured before deciding.  On a run of ~5 tau the clamp is mostly
irrelevant: above `timestep ~ 3e-4` the ERROR CONTROLLER becomes binding
and `max dt` stalls at 2.97e-4 however much slack the clamp is given.
But on a run of 100 tau, ~99% of it quiescent, it is clamped at every
setting -- **1027 steps to traverse a dead-flat solution whose error is
4.4e-16**.  That cost is paid by exactly the circuits that idle, which is
most mixed-signal ones.

The DEFAULT is not changed: doing so would move every waveform in the
package for a benefit only idle-heavy runs see.  SPICE's own default for
the equivalent knob (`TMAX`) is `(tstop - tstart)/50`, which a caller can
now ask for directly.
Decoupled from `timestep` by owner decision (see the Parameter):
None means SPICE's TMAX default, tend/50.  The old
max_step-below-timestep refusal is moot -- `timestep` no longer
promises an output density, so a cap below it just clamps the
opening step like any other.

The comment on the delay-element cap, before the move:

`TLine` interpolates its history at `t - TD`; with `dt` comparable to
`TD` there is nothing to interpolate between, and the delay simply comes
out wrong -- measured 2.00x TD at dt = 1e-9, 4.00x at 2e-9 and 8.00x at
5e-9 under `fixed_timestep`, with no warning of any kind.  The adaptive
controller usually rescues it (1.01-1.05x), which is exactly why it went
unnoticed: the defect only bites the configuration nobody checks.

The comment on the LTE `abstol`, before the move:

`abstol`/`xtol`.  This is the `xtol` one; using `_newton`'s `abstol` here
silently applied iabstol (1 pA) as a *voltage* tolerance to every node,
which is what made a larger `vabstol` have no effect on node rows at all.

It reads `lte_vabstol`/`lte_iabstol`, NOT `vabstol`/`iabstol`.  Sharing
the Newton parameters was decision 0.3a's defect: one knob could not be
moved for the controller without silently moving Newton's convergence
criterion with it.

The comment on the step family, before the move:

for every method.  A Nordsieck GLM is not a Runge-Kutta method, but
it is self-starting, keeps no charge ring and delivers its own
estimate through `_rk_est`, so it is a stage family (it once
reached the LMM controller and died in `get_diff`).

The comment on the opening ramp, before the move:

The opening ramp exists to stop the ONE step the controller cannot check
from dominating the run.  Under `fixed_timestep` there is no controller
and the step is never adapted, so ramping would not open small and grow
-- it would run the ENTIRE simulation at `timestep*1e-3`, a thousand
times more steps for a result the caller explicitly asked to be
uniform.  Caught by `test_transient_RLC` and three others.

The comment on growth retries (F14), before the move:

F14 (doc/transient_review_260820.md): a lower-band GROWTH retry --
the controller redoing a too-ACCURATE step larger -- is a voluntary
redo, not a failure, and must not trip the force-accept.  Before
the split, three consecutive growth retries during the opening ramp
of a QUIESCENT circuit reached the force-accept path: 2 spurious
warnings on a settled RC with the band at (0.5, 3.0).  Over-tolerance
rejections strictly shrink in every controller, growth retries return
only behind a strict-growth guard, so `h_next > h` tells them apart;
the family's `max_retries` bounds both against pathological
alternation.

The comment on breakpoints under `fixed_timestep`, before the move:

STAGE 4h -- UNDER `fixed_timestep` THE GRID WINS: a
breakpoint no longer moves it (a truncation used to be
permanent and collapse the step geometrically: 292 steps
for 30), but crossing one still drops the order.  `<=`,
TOLERANCED: `t` accumulates by `+= h`, so an edge exactly
on a grid point is a float knife-edge (measured: the drop
fired at the first edge and missed every later one).

### `_rescue_step`

2026-09-27 (moved from the code):

The last paragraph of the docstring, before the move:

⚠ IT USED TO BE UNREACHABLE FOR EVERY STAGE METHOD on the default
path: the Runge-Kutta loop raised where the LMM loop rescued, so the
ladder `_rk_step_coupled` carries could only fire under
`fixed_timestep=True`.  One loop, one ladder.

The comment on `_honours_continuation_rescue`, before the move:

⚠ WHETHER THIS PATH REACHES A LADDER AT ALL (`_honours_continuation_
rescue`, consulted here since 2026-09-27 -- it was written to guard
this message and had no caller): the error must not claim a rescue
that was never attempted, and a success is a rescue only where one
could run


## `jaxtransient.py` -- module level

### `TLINE_HISTORY_DEPTH`

2026-09-27 (moved from the code):

The comment above `TLINE_HISTORY_DEPTH`, before the move:

Depth of the per-TLine delay-line ring buffer, in accepted steps.  It was the
bare literal `10000` in seven places, three of them inside modulo arithmetic,
which made it impossible to see that the buffer is a RING: past this many
accepted steps `tline_head` wraps and `interp_tlines` starts interpolating
against entries from a previous lap.  `cond_fun` stops at the end of the buffer
and returns whatever is there, so the result is a plausible wrong waveform with
no error and no warning.  `JAXTransient.solve` now checks for the wrap and
raises; sizing the buffer from the run instead of fixing it belongs to stage 9.


## `jaxtransient.py` -- `NewtonState`

### `converged`

2026-09-27 (moved from the code):

The comment above the `converged` field, before the move:

STAGE 9(e).  The loop exits on `F_norm <= conv_tol` OR `iters >= maxiter`
and the caller took `x` either way, so an unconverged iterate was committed
as the step's solution and its LTE computed from it.  Measured on a LINEAR
RC at maxiter=1: 4.97e-2 max error against 4.28e-3, reported as nothing at
all.  A traced loop cannot raise, so the fact travels out as a flag.


## `jaxtransient.py` -- `TransientState`

### `sig_max`

2026-09-27 (moved from the code):

The comment above the `sig_max` field, before the move:

STAGE 9(c) -- the running maximum of |x| over ALL past steps and ALL
unknowns.  This is a commercial simulator's `sigglobal`, which `Transient` already ships as
its default `relref`, and it is what makes an absolute LTE floor of 1e-12
safe: under `pointlocal` -- each unknown against itself, now -- a node
carrying no signal drives `ref -> 0`, the tolerance collapses to the floor,
and the controller chases numerical noise on a quiet node.  The CPU hid that
for a while by raising the floor a millionfold; `sigglobal` is the fix, and
this field is what carries it through a traced loop.

### `n_rejected`

2026-09-27 (moved from the code):

The comment above the `n_rejected` field, before the move:

STAGE 9, gate 9-1(b) -- how many steps the controller REJECTED.  Nothing
reported this, so "a step is actually rejected" -- one of the three CPU
step-control gates -- was not expressible on this backend, which is the
asymmetry that let a copied LTE defect survive being fixed twice.  A
rejected step advances neither `t` nor `step_idx`, so it leaves no trace in
the output buffers; it has to be counted where it happens.

### `n_nonconverged`

2026-09-27 (moved from the code):

The comment above the `n_nonconverged` field, before the move:

Steps rejected because the Newton solve did not converge, and steps
ACCEPTED anyway because `dt` was already at the floor.  The second is the
one to read: it counts places where the returned waveform is not a solution
of the circuit equations, which is exactly what 9(e) found happening
silently.

### `n_forced_lte`

2026-09-27 (moved from the code):

The comment above the `n_forced_lte` field, before the move:

F19(c) (doc/transient_review_260820.md): a step at the dt floor whose
Newton CONVERGED but whose LTE still failed used to be accepted with no
trace at all -- the CPU counts these as force_accepts and warns that
the accepted error is unbounded.  Counted here for the same reason.


## `jaxtransient.py` -- module level

### (the Phase 1 integrator note)

2026-09-27 (moved from the code):

The comment above `backward_euler_step`, before the move:

STAGE 9(a) -- these three were transcriptions of `integrator.py`'s.  They now
call the one definition in `_lte_kernels`, which is plain arithmetic and so
traces under `jax.jit` unchanged.

Euler and Gear-2 were already bit-for-bit the same expression, so those runs do
not move.

THE TRAPEZOIDAL BRANCH WAS DELETED once (review hygiene): no production
path could reach it (at the time, `eval_method` was hardcoded 'gear' at
both call sites; P6 has since exposed the choice as the `integrator`
Parameter), and its LTE formula was the uniform-grid one -- wrong the
moment the step changed, had anyone ever wired it up.  Dead-but-plausible
solver branches are exactly how the 3/4-optimism defect survived twice.

REBUILT 2026-09-01, and the deletion is why the rebuild is trustworthy:
restoring the branch meant restoring only the COMPANION (one shared kernel
call), while the estimator had to be written against the charge, not
transcribed from the g-based form the old branch used.  Had the branch
survived unreachable, the wrong formula would have come back with it.

The record said for weeks that "a variable-step trap estimator has not
been written".  That was wrong: `trapezoidal_lte` and
`third_divided_difference` were already in `_lte_kernels`, already used by
the CPU, already plain arithmetic that traces.  What had not been written
was the WIRING -- a much smaller claim, and the difference between them is
about a day of work someone did not do.

### `effective_first_order`

2026-09-27 (moved from the code):

The docstring body, before the move:

True while the run has fewer than two completed steps recorded
(``h_history[0] == 0`` or ``h_history[1] == 0`` -- RUN-global facts, the
buffers carry across chunk boundaries) and when the previous accepted
step landed on a breakpoint (``force_first_order``, F11).  ONE
definition, consumed by the integration dispatch, the LTE estimator, and
the step-size exponent alike -- F19's point
(doc/transient_review_260820.md): the integration used to fall back to
Euler dynamically while the estimator stayed statically Gear-2, so an
order-1 step was scored by the order-2 formula on a zero-seeded history
with the wrong exponent handed to the step controller.

Deliberately NOT ``step_idx < 2``: step_idx is CHUNK-local and resets at
every chunk boundary, so the old predicate re-dropped the order (and,
once the estimator followed it, re-scored steps) at each boundary --
measured as chunking changing a 59-step run to 61 steps.  The history
buffers are the run-global truth.

### `compute_integration`

2026-09-27 (moved from the code):

The comment in `do_trap`, before the move:

Trapezoidal needs the previous COMPANION CURRENT, which Euler and
Gear-2 do not -- `iq_history` has been maintained all along (rolled
on accept only, which is exactly the semantics the recursion
needs) and until now was read by nobody but the LTE estimator.

### `_tline_emfs`

2026-09-27 (moved from the code):

The second paragraph of the docstring, before the move:

Extracted from newton_inner_loop's closures (TLine-under-coupled/PCNR
port): the delay-line history lives in TransientState and is updated by
the shared accept machinery, so the only thing that kept the other
assemblies TLine-blind was this code being trapped in one function.  The
derivatives are the linear interpolant's segment slope -- free from the
same lookup -- and are what Fang's ``p = df/dh`` needs from a source
whose value depends on the step size through ``t - TD``.

### `tline_stamp_correction`

2026-09-27 (moved from the code):

The first paragraph of the docstring, before the move:

THE STANDARD-PATH DEFECT THIS FIXES (found while porting TLine under the
coupled/PCNR paths, present since the JAX TLine landed): ``TLine.G``
selects its stamp on ``len(self.history) == 0`` -- Python state that only
the CPU's accept_step writes.  On this backend the traced ring buffer
replaced accept_step, the Python history stays empty forever, and every
assembly therefore stamped the DC form (the line as an ideal short:
v1 - v2 = 0, i1 + i2 = 0) while the injected EMFs assumed the transient
form (v - Z0*i = e).  Measured on a matched 1 V line: |v(far end)| =
24.5 V, delay 0 -- and NO e2e JAX TLine test existed to see it.

### `_adaptive_ladder_traced`

2026-09-27 (moved from the code):

The failure-bookkeeping comment in `body`, before the move:

Failure bookkeeping, in three phases: escalate, then DESCEND,
then refine.

⚠ THE DESCENT PHASE WAS MISSING, AND `e_end` WAS THEREFORE A LIE.
With no rung yet landed the driver escalated toward `e_max` and,
once that was spent, refined just below `e_start` -- halving from
`step/2` and giving up at `min_step`.  So the exponents it could
ever try were

    {e_start, e_start+step, ... e_max}  and
    {e_start-1, e_start-0.5, e_start-0.25}

and NOTHING BELOW, however low `e_end` was set.  Measured on a
synthetic rung that converges only for g <= 10^k: the ladder
landed for k >= -1 and failed for every k <= -2, against an
`e_end` of -12 advertising a search down to 1e-12.

That is why Psi-tc looked "DC-calibrated" in transient: its first
rung is at g = 1, and if that one does not converge the ladder
could not reach the g that would.  Re-scaling the exponents (the
recorded guess) would only have moved the same one-decade window
somewhere else -- measured, a half-decade shift rescued 4 of 4
test circuits while the shipped grid rescued 1 of 4, which looks
like scale sensitivity and is actually this.

Descending keeps the step and walks `e_low` down to `e_end`, so
the whole declared range is searched before giving up.  It can
only make the driver reach further: the phase runs exactly where
the old code was already about to fail.

### `newton_inner_loop`

2026-09-27 (moved from the code):

The F6(b) comment above `cond_fun`, before the move:

F6(b) (doc/transient_review_260820.md): THE CPU'S PER-ROW CRITERIA,
replacing `conv_tol = abstol + reltol * F_norm0` -- which was wrong
three ways at once.  Flavour: the scalar floor threaded here was
vabstol, a voltage, against a residual of KCL currents (F6(a) picked
the majority flavour; the per-row vectors finish the job -- iabstol on
node rows, vabstol on branch rows, and the transposed pair for the
update test).  Reference: relative to the INITIAL residual, so a bad
predictor loosened the target -- false convergence -- while a good
one collapsed it toward the absolute floor -- spurious
non-convergence.  Norm: a summed L1, so one badly-failed row diluted
with circuit size.  The body below mirrors StandardNewton: per-row
conv_f against reltol*I_scale + abstol, per-row conv_x against the
step actually taken, both computed on the consistent (F(x), dx) pair
-- which also retires stage 9(e)'s trailing re-evaluation, since the
converged flag now travels in state and refers to the pair that
produced the returned x, exactly as on the CPU.

The F16 comment above the update clamp, before the move:

F16 (doc/transient_review_260820.md): the update clamp is a
PARAMETER, default disabled (jnp.inf -> alpha 1.0).  The hardcoded
0.5 V made any swing beyond ~maxiter*0.5 V non-convergent by
construction with a false "failed to converge" diagnosis -- a 48 V
rail was unreachable in 100 iterations -- and under the per-row
criteria (F6(b)) it even made a LINEAR 1 V step need multiple
clamped iterations.  A flat voltage clamp punishes the linear part
of the circuit for the nonlinear part's sins; per-junction
limiting (pnjlim's shape) is the eventual replacement if the JAX
element set grows junctions.

The comment above `conv_f`, before the move:

THE REDUCED ROWS, not all of them (P11's port found this): the
reference row is not an equation of the solved system -- its
residual is whatever KCL imbalance the source terms carry (an
unbalanced provided_function put the full injection there), and
scoring it can NEVER be satisfied by any x, which livelocked the
run at maxiter on every step, silently force-accepting at the dt
floor.  The CPU's Newton has always tested the reduced system.

### `ywr_error_ratio`

2026-09-27 (moved from the code):

The comment in `_gear_branch`, before the move:

STAGE 9, gate 9-1(a) -- THE 3/4 OPTIMISM, FOUND A THIRD TIME.

This was YWR's Table I GEAR2 residual,
    -(1/8) ((h1+h2)/(h1 h2)) (h2 g_n - (h1+h2) g_{n-1} + h1 g_{n-2}),
which on a uniform grid reduces to -(1/4) h^2 q''' against a true BDF-2
local truncation error of -(1/3) h^2 q'''.  So it reported 3/4 of the
error at every step -- the solver was 25% optimistic about its own
accuracy, on the default eval_method of both entry points.

The CPU found and fixed this in stage 4i and the fix never crossed to
this file, which is precisely the divergence stage 9 exists to close:
the same defect has now been found three times in two transcriptions.
Measured here at 2.5e5 against the CPU's 3.333e5 for q''' = 1e6.

The form below is 4i's: estimate q''' from the second divided difference
of the companion current and multiply by the method's own error
constant, so the coefficient is derived rather than transcribed.

The comment at the head of the trapezoidal branch, before the move:

⚠ THIS IS NOT THE GEAR-2 BRANCH WITH A DIFFERENT CONSTANT, and the
difference is the whole reason the old trapezoidal branch was
deleted rather than repaired.

The Gear-2 branch above differences the COMPANION CURRENT `g`.
Trapezoidal must not: its recursion
    iq_n = 2 (q_n - q_{n-1})/h - iq_{n-1}
carries an UNDAMPED (-1)^n homogeneous mode, so differencing `g`
measures that mode rather than the truncation error -- recorded in
`_lte_kernels` and `integrator.py` as a 1.9x swing produced by step
history alone.  The estimator therefore differences the CHARGE:
`third_divided_difference` returns q'''/6 and `trapezoidal_lte`
consumes exactly that, so the pairing is stated once (the same
discipline that kept the 3/4 optimism from recurring here).

The deleted branch used YWR Table I's TRAP entry,
`Eg = -(1/6)(g_n - 2 g_{n-1} + g_{n-2})` -- a uniform-grid formula,
and `g`-based, so it was wrong twice over.

The P7 comment, before the move:

P7: the CPU's relref modes, selected by a trace-static string.  The
reference is a per-row VECTOR now; 'sigglobal' collapses each UNIT
GROUP (node voltages vs branch currents) to its own maximum -- the
old scalar max over everything mixed volts and amps, which the CPU
comment warns silently disables node error control on circuits with
large branch currents.  `local` includes the current iterate, exactly
as the CPU's `_reference` does.

The comment above `etol`, before the move:

`lte_abstol` is a per-row VECTOR (volts on node rows, amps on branch rows),
built by `JAXTransient._lte_abstol`.  A scalar here applies one physical
kind of tolerance to rows of another -- 0.3a's residual-vs-solution defect,
which 9(c) exists to avoid repeating.  A scalar still broadcasts, so a
direct caller that passes one gets the old behaviour rather than an error.

### `collect_breakpoints`

2026-09-27 (moved from the code):

The second and third paragraphs of the docstring, before the move:

STAGE 9(d).  This replaces two copies that disagreed.  ``solve`` iterated
``for elem in cir.elements`` -- a **dict**, so it yielded string keys and
``hasattr('V1', 'next_event')`` was False, giving **0 breakpoints, always**.
``solve_batched`` iterated ``.items()`` and was correct, and therefore hit the
second bug instead: ``Pulse.next_event`` returned ``t`` itself at ``t = 0``, so
the enumeration never advanced and the call **hung**.  Two bugs cancelling, one
per copy, which is the argument for having one copy.

``Pulse.next_event`` is fixed at the source, but the progress guard below stays:
it is the difference between a wrong breakpoint list and a wall-clock hang, and
only one of those is diagnosable from a stack trace.

The TLINE WAVEFRONT ARRIVALS comment, before the move:

TLINE WAVEFRONT ARRIVALS: a source corner reaching a delay line
re-emerges at the far end TD later -- a from-zero kink in an ALGEBRAIC
variable that no element reports.  Registering {corner + TD, corner +
2*TD} makes those steps breakpoint-landings, which the coupled path's
post-breakpoint grace (see fang_inner_loop) can then absorb -- the two
mechanisms only work TOGETHER: arrivals-as-breakpoints alone was tried
first and falsified, because a registered kink without the grace still
livelocks on the h-independent relative LTE.  Deeper bounce ancestry
is truncated, as every SPICE truncates it.

### `PI_K_I`

2026-09-27 (moved from the code):

The comment above `PI_K_I, PI_K_P`, before the move:

Gustafsson's gains, as `PIController.__init__` defaults them on the CPU.
⚠ PER UNIT ORDER -- `k_I/p` and `k_P/p`, not the bare numerators.  Used
undivided the loop is linearly UNSTABLE (spectral radius 1.12 at p=2, 1.78
at p=3), and it does not look like divergence: the growth clamp converts the
growing oscillation into a permanent period-2 limit cycle, measured on the
CPU as h alternating 0.857/0.429 for the length of the run.  The only test
there asserted `len(steps) > 10`, so it was invisible.  A test pins these
against the CPU object so the two definitions cannot drift.

### `calculate_next_dt`

2026-09-27 (moved from the code):

The F17 comment, before the move:

THE CPU'S LAW, F17 (doc/transient_review_260820.md): aim at
target = safety**p, not at err = 1.0 -- the rejection threshold
itself.  Aiming at the edge meant every successor step that landed a
hair over 1.0 was a rejected re-solve; measured on rc-vsin at
reltol 1e-4/1e-6 as a 44%/40% rejection rate, cut to a few percent by
the safety margin.  `(0.9**p / err)**(1/p) = 0.9 * err**(-1/p)` --
identical arithmetic to a bare safety multiplier, written in the CPU's
vocabulary (stepcontroller._band_target) so the backends read as one
law.  The clamps are the CPU's named constants rather than re-derived
literals: the shrink floor was 0.1 here against MIN_SHRINK_RATIO=0.2.

### (the header of Fang's coupled time step, above `FangState`)

2026-09-27 (moved from the code):

The header comment, before the move:

Fang's coupled (x, h) time step -- the 'approx' branch, traced (P19).

A faithful port of Transient._fang_timestep_inner's DEFAULT path, scoped the
way doc/backend_parity_260821.md P19 records: sec. 3.4's error-ratio step
correction with the eq (18) solution update, hold_h for imposed step sizes,
the within-point excursion clamps and the thwarted-shrink saturation test.
(The 'bordered' eq (12) branch was never ported, and was retired on the
CPU on 2026-09-27.)  The LTE degree follows the effective order (F19): both degrees are
computed and selected, which keeps every shape static under the trace.

⚠ CORRECTED 2026-09-01.  This list used to name three more items, and all
three had since been done -- PCNR-inside-Fang (it is right below, keyed on
`pcnr_meta`), grid_locked "no fixed_timestep on this backend yet"
(`fixed_timestep` is a Parameter and gives bit-equal fixed-grid waveforms;
only the COMBINATION with coupled_lte is refused, in `solve`), and cir.limit
"PCNR's job, next stage" (PCNR shipped, both views).  A scope note written
beside the code it scopes is read as current, so it outranks the docstring
in a reader's mind while being the thing nobody updates.  **When this and
the JAXTransient class docstring disagree, the docstring is the ledger and
wins.**  What is genuinely still CPU-only lives there.

### `fang_inner_loop`

2026-09-27 (moved from the code):

The comment above `no_hist`, before the move:

The band cannot be evaluated before two accepted points exist.  The
post-breakpoint grace on delay-line circuits arrives through this same
test: the accept path EMPTIES h_history on a breakpoint landing when
the circuit has TLines (see do_accept), so the next step reads as
history-free here.  An explicit `or force_first_order` term was tried
and reverted -- it duplicated the reset's effect on TLine circuits and
extended a band-blind step to every breakpoint on circuits that never
needed it, the exact accuracy cost the CPU measured at 9.8e-4 ->
5.4e-3 median on the pulsed RC when its reset ran ungated.

The comment above `_vector`, before the move:

PCNR comes in two views and fang now takes both.  The PAIR view
(`pcnr_junctions`) lists pn-junction probes; the DEVICE view (sec. 49)
lets a device own `m` unknowns of mixed kind, which is what every real
compact model needs.  Until 2026-09-01 only the pair view was wired
here and the device view reached this function as the literal string
'vector' used as an array index.

The comment above `p_e` in `solve_h`, before the move:

THE dh-DERIVATIVE MUST DESCRIBE THE METHOD THE RESIDUAL USED.
This used to be euler-or-gear2 unconditionally, while the
residual above is assembled with `method=eval_method` -- so
`integrator='euler'` with `coupled_lte=True` integrated with
Euler and then, on every non-first-order step, differentiated a
GEAR-2 companion with respect to h.  Fang's step-size Newton was
solving a slightly wrong equation; it still converged on x, which
is why nothing ever failed and the mismatch survived.

The comment above `eps`, before the move:

eps scales with the STEP, not with absolute time: an
absolute-scaled eps (1e-9 at t << 1) straddled the first
pulse edge from inside the opening ramp, so the derivative
saw the future corner and eq (18) corrupted the solve --
measured as the coupled path dying at t ~ 1e-12 on any
driven TLine circuit.

### (the VECTOR PCNR header, above `PcnrVectorDevice`)

2026-09-27 (moved from the code):

The header comment, before the move:

VECTOR PCNR (roadmap sec. 49): the DEVICE-shaped path.

`_junction_arrays` below is the pnj-only PAIR view, and it keeps its own
consumer -- the gmin ladders put a conductance across each junction pair on
the ordinary `pcnr=False` path, where a gmin across a FET's `vgs` would be a
gate leak.  Vector PCNR's unit is the DEVICE: one unknown per limited
quantity, `m` of them, of mixed kind.  Measured on a MosLevel1Hdl,
`pcnr_devices` reports 1 device with m=3 where `pcnr_junction_pairs`
reports 2 pairs -- different objects, different laws -- which is the
distinction Stage 2 turned on when `pcnr=True` on a MOSFET differential
pair fell through to the ordinary solver in silence.

### `_junction_arrays`

2026-09-27 (moved from the code):

The second paragraph of the docstring, before the move:

The traced loop used to rebuild EVERY junction itself as
``IS*(exp(v/VT)-1)`` with a single global VT, reading ``IS`` by name.
That silently gave a device whose saturation current is not called
``IS`` a junction carrying no current, and it made multi-junction and
charge-storing participants impossible here while the CPU path
accepted them.  Asking the device instead removes all three limits at
once: any shape, any number of junctions, and charge simply stays in
the MNA block exactly as it does on the CPU (this assembly subtracts
the junction term at the node voltage and adds it at ``v_lim``; it
never touched the charge).

### `outer_time_loop`

2026-09-27 (moved from the code):

The comment above `t_eps`, before the move:

The same epsilon `calculate_next_dt` uses to decide a breakpoint is "already
reached".  They disagreed: after 500 steps of 1e-5 the accumulated `t` sits
~1e-18 short of `tend`, which the breakpoint filter treats as arrived (so it
drops the clamp) while `t < tend` treats as not arrived (so it takes another
FULL step).  That is the residual overshoot -- exactly one timestep, measured
at t[-1] = 5.010e-3 for a requested 5e-3.

The comment in `time_cond`, before the move:

A forced NON-converged accept is already an unconditional raise
after the chunk (stage 9(e)) -- no run continues past one -- so
the chunk exits the moment it happens instead of marching to
chunk_size at the dt floor first.  P22's mask port measured why
this matters: with the algebraic rows masked, a cold-start
circuit whose Newton fails at every h reached the floor honestly
and then ground out 500 forced steps at ~1 s each on GPU (the
coupled point is ~100 serial kernel launches per attempt) before
the raise -- an effective hang.  Pre-mask, the algebraic-row
error was ACCIDENTALLY load-bearing as a step governor that
rescued the Newton; the mask removed the accident, this exit
replaces it with the deliberate escape.

The STAGE 9(f) comment above the controller Jacobian, before the move:

STAGE 9(f) -- ONE ESTIMATOR.  There used to be a charge-domain branch
here, selected by `lte_formula='classic'`.  It was deleted rather than
repaired: its tolerance applied `lte_abs = 1e-6`, a VOLTAGE floor, to a
CHARGE -- one microcoulomb, against node charges of pico- to
femtocoulombs -- so the normalized error could not reach 1 and no step
was ever rejected.  The controller ran open-loop.  Under `solve()` that
degenerates to a fixed-step run (`dt_max = timestep`), which is why it
looked plausible; under `solve_batched` (`dt_max = tend/10`) it also
costs accuracy.  0.2b measured the `J^-1` mapping this branch avoided at
1-3% of a step, well under its own 10% keep-it threshold.

The F19(c) comment above `forced_lte`, before the move:

F19(c): converged Newton, failing LTE, at the floor -- accepted
with an unbounded truncation error, which the CPU force-accept
warns about and this path used to swallow silently.

The comment at the head of `do_accept`, before the move:

`state.t + state.dt`, NOT `state.t`.  This step has been accepted, so
the next one starts where this one ENDS; passing the old time sized the
next step to reach the breakpoint from the PREVIOUS position and so
overshot it by about one step.  Measured before the fix: t[-1] of
5.0559e-3 against a requested tend of 5e-3, where `Transient` lands on
it exactly.  `tend` is itself in `t_breaks_array`, which is why the
overshoot showed up at the end of every run.

The comment in the `coupled and fixed_timestep` branch, before the move:

Under a locked grid fang solved no h -- the step was held --
so there is no "solved step carries forward" to honour and
the grid branch below is the right one.  Ordering matters:
`coupled` used to win unconditionally, which would have
carried the HELD step forward as if it had been solved and
quietly re-derived the grid from itself.

The comment in the `coupled` branch, before the move:

The SOLVED step carries forward; breakpoint truncation for
the NEXT point happens at the next fang entry.  Solved is
the operative word (P22's mask port measured the hole): a
NON-converged fang's h is not a solved h -- with the band
below tolerance it GROWS within the point, so a forced
floor-accept came back at ~2x dt_min, the at_floor escape
went un-sticky, and the run crawled at ~100 fang
iterations per 1e-18 s -- an effective hang where the
unmasked path reached its post-chunk non-convergence
raise in seconds.  A forced accept keeps the floor.

The ARC 7 comment, before the move:

ARC 7 (idtmod.md sec. 5.3): the in-trace wrap-crossing dt
cap, which that section listed as future work because
`t_breaks_array` is static and a state-dependent event
cannot enter it.  It enters HERE instead: not as a
breakpoint, but as a cap on the next step, computed from
state the trace already holds.  Branchless -- no host
control flow, so it compiles inside `lax.while_loop`.

⚠ WHAT THIS DOES AND DOES NOT BUY, measured before it was
built.  A wrap is a discontinuity in the OUTPUT MAP, not in
the ODE: the gauge shift above keeps the state continuous
and the sawtooth is an exact `floored_wrap` of it, so there
is no kink for the integrator to resolve.  Accordingly it
buys NOTHING in accuracy or cost -- against a tight
reference on a varying integrand, JAX and the CPU agreed to
within 1-5% of each other's error (6.171e-03 vs 6.112e-03 at
reltol 1e-4; 2.925e-04 vs 2.775e-04 at 1e-6) at the same
step counts, and sec. 5.3's predicted step-collapse at wraps
DID NOT REPRODUCE (70 points against the CPU's 74, no
rejection storm).

What it buys is SAMPLE PLACEMENT.  Measured on the exact
ramp, distance from the true corners: CPU (which has
breakpoints) 3.55e-15, JAX 3.29e-02 -- up to three output
timesteps past the corner.  That matters to a consumer that
RESAMPLES or edge-detects the output, because interpolating
a sawtooth across an unmarked corner returns values the
signal never takes (sec. 3.4's consumer discontinuity).

THE PREDICTION is linear, from the step just accepted --
the CPU's own rule in sec. 5.3 ("the linearly predicted
crossing time").  It only has to BRACKET the corner; the
step is capped, not landed exactly, and the next step
re-predicts from there.

The comment above `_eps` in the wrap-crossing cap, before the move:

A point sitting exactly ON the boundary it is leaving
must travel a whole modulus, not zero -- otherwise the
cap collapses the next step to dt_min at every wrap,
turning an accuracy-neutral change into a cost one.

⚠ DEFENSIVE AND UNEXERCISED, and said plainly rather
than left to look tested.  The case needs a DESCENDING
state landing exactly on `offs`: ascending, a landing
on the top is shifted to `offs` and the next boundary
going up is a full modulus away, so no zero arises.
Descending was built and driven (v = -1 into modulus 1,
three wraps) and the branch still never fires, because
the cap lands just SHORT of the boundary in floating
point -- the gauge shift has then already wrapped the
state to the top and the distance is again a modulus.
Removing this branch changed that run not at all
(64 points, corners at 8e-15, both ways).  It is kept
because "unreachable in the arithmetic I could
construct" is not "unreachable", and the failure it
prevents is a step collapse at every wrap.

The comment above the carried counters, before the move:

CARRIED, not defaulted.  `do_accept` builds a fresh
TransientState rather than `_replace`-ing, so every field it
omits silently reverts to the NamedTuple default -- and these are
cumulative counters.  Omitting them reset the rejection count to
zero on each accepted step, which made "16 configurations, zero
rejections" a measurement of the bug rather than of the solver.

The F11 comment on `force_first_order`, before the move:

F11: the CPU's breakpoint discipline, ported.  A step that
LANDS on a breakpoint (calculate_next_dt truncates onto
them, so landing is exact to t_eps) marks the NEXT step
"do not trust a 2nd-order polynomial through this point":
effective_first_order then drops both the integration and
the error estimate to order 1 for exactly one step, instead
of differencing a g-history that straddles the corner --
which cost a rejection burst at every edge (measured on a
VPulse RC before this line: 38 rejections in 183 accepted
steps, edge-synchronous).

The comment in the coupled branch of `do_reject`, before the move:

The CPU coupled path's outer retry: h *= 0.25 whatever the
failure flavour (Newton or held-step LTE) -- applied to
the ENTRY h, not to fang's returned h.  fang GROWS h
within the point when the (masked) band reads the error
as tiny, by up to MAX_GROWTH=5x; shrinking the grown h by
4x nets a 1.25x GROWTH per reject cycle, so a cold-start
circuit whose Newton fails at every h never reached the
dt floor where the forced-accept escape lives -- an
in-trace livelock the CPU cannot express (its retry loop
is Python-bounded at 10 attempts).  Measured behind
P22's mask port; the mask only made it reachable.


## `jaxtransient.py` -- `JAXTransientStatistics`

### (class docstring)

2026-09-27 (moved from the code):

The class docstring, before the move:

What a JAX transient run did, as opposed to what it returned.

The subset of `transient.TransientStatistics` a traced `while_loop` can count.
`rejected_steps` is the load-bearing one: a rejected step advances neither `t`
nor `step_idx`, so it leaves no trace in the output buffers, and without this
counter the CPU gate "a step is actually rejected" could not be stated on this
backend at all.  That is the asymmetry stage 9 exists to remove -- it is why a
copied LTE defect had to be found and fixed twice.


## `jaxtransient.py` -- `JAXTransient`

### (class docstring)

2026-09-27 (moved from the code):

The CPU-only bullet of the class docstring, before the move:

* **CPU-only**: ⚠ EVERY STAGE AND MULTIVALUE METHOD (corrected
  2026-09-16; the line below had said "nothing" since 2026-09-01 and was
  read as such by the roadmap).  This backend takes 'gear', 'euler' and
  'trap' -- the three LMM companions -- and refuses anything else in
  three places; `radau`, `trbdf2`, `esdirk43` and the GLM family appear
  nowhere in this file.  So a per-stage Newton has never run here, and
  putting one of them on this backend is a NEW capability rather than a
  port of a missing case.  (The coupled 'bordered' eq (12) branch, once
  refused here, was retired on the CPU on 2026-09-27: after its
  double-counted term was removed it took the same steps as 'approx'.)
  Its whole purpose is to trade a few more time points for fewer Newton
  iterations; under Gear-2 it loses on both.  Porting it would import a
  defect and call it parity.  See the P19 note for what would reopen it.

  Two items left this list on 2026-09-01.  *Trapezoidal* is now
  `integrator='trap'` -- the record had said a variable-step estimator
  "has not been written" when the kernels existed and only the wiring
  did not; the part that genuinely had to be written is that the
  estimator differences the CHARGE, since differencing the companion
  current measures the trap recursion's undamped (-1)^n mode.  *coupled
  + fixed_timestep* now works: `grid_locked` reduced to one flag, since
  an over-band HELD step is normally reported unconverged so the caller
  shrinks, and under a caller-imposed grid shrinking is not available.

### `parameters`

2026-09-27 (moved from the code):

The comment above `parameters`, before the move:

STAGE 9(b)/(c) -- TOLERANCES ARE SETTABLE, AND THEY ARE THE CPU'S.

This class declared no tolerances at all, so `JAXTransient(cir, reltol=1e-6)`
raised `KeyError: 'parameter reltol not in parameter dictionary'` and there
was no supported way to ask for a tighter run.  The values were hard-coded
at the kernel's own defaults (`reltol=1e-3, abstol=1e-6` for Newton;
`trtol=7.0, lte_rel=1e-3, lte_abs=1e-6` for the LTE, never threaded from
the caller at all).  Two of those disagreed with `Transient`'s shipped
defaults, so the two backends were solving to different accuracies while
presenting as the same analysis.

THE NAMES AND FLAVOURS ARE `Transient`'s DELIBERATELY.  0.3a's residual-vs-
solution defect on the CPU came from applying one scalar absolute tolerance
to rows of different physical kinds; the plan's 9(c) warns in as many words
that threading a scalar here would re-create it.  So the absolute
tolerances are VECTORS built the same way `Transient` builds them --
`lte_vabstol` on node rows, `lte_iabstol` on branch rows, because the LTE
lives in the SOLUTION domain where nodes carry volts and branches carry
amps.  `relref` is not offered yet: the CPU's `sigglobal` needs a reduction
over the whole vector inside the traced loop, and that is 9(c)'s remaining
work rather than something to fake with a scalar.

The comment on the LTE acceptance band, before the move:

The LTE acceptance band, with the CPU's 'auto' sentinel semantics
(F5): 'auto' resolves to Fang's (0.7, 3.0, 0.15) on the coupled
path; explicit values pass through verbatim.  The standard JAX
path does not read the band yet (parity item P8).

The comment on `max_dv`, before the move:

F16: opt-in, OFF by default -- the old hardcoded 0.5 V made any
swing beyond ~maxiter*0.5 V non-convergent by construction.

The comment on `max_dv_step` (its first eight lines are about `timestep_max`), before the move:

STAGE 9(g).  Ported from `Transient`, same name, default and validation.
DECISION D2 -- the same knob as the CPU's, same name and same
default, because 9(d) and 9(g) both exist because these two backends
had drifted apart on exactly this kind of detail.  OWNER DECISION
(2026-08-21): decoupled from `timestep` and renamed -- `timestep`
doubling as the cap made gentle circuits step-cap-limited, where no
tolerance knob could move the run (measured: identical 209-step
rc-vsin runs at reltol 1e-4 and 1e-6).

The P1 comment, before the move:

P1: `uic`/`minstep` were **kwargs reads -- `solve(uicc=True)` ran
silently with defaults, the dead-knob defect class at the call
boundary.  Declared as Parameters (CPU names, CPU defaults); the
solve()/solve_batched() arguments remain as explicit per-call
overrides, None meaning "use the Parameter".

The P5 comment on `analysis`, before the move:

P5's other half: the inherited `analysis` default is not 'tran',
and a threaded self.par.analysis of anything else makes every
source's u() return zeros -- caught by the rc-charging gate as an
all-zero waveform the moment the hardcode was removed.  Same
re-declaration the CPU Transient makes.

The P6 comment on `integrator`, before the move:

P6: the integrator choice, reachable at last -- the traced loop
implemented euler and gear all along, but eval_method was
hardcoded at both call sites.  String-valued on this backend (the
traced kernels select by name; there is no Integrator instance to
hold state), which is also why P17 can refuse strategy objects
without refusing this.

'trap' joined them 2026-09-01.  This comment said until then that
trapezoidal "stays CPU-only until someone ports a VARIABLE-STEP
trap estimator" -- kept here because the sentence was wrong in an
expensive way: the estimator existed, in `_lte_kernels`, already
used by the CPU and already traceable.  Only the wiring was
missing.  ⚠ The trap estimator differences the CHARGE, not the
companion current, because the trap recursion carries an undamped
(-1)^n mode; see `ywr_error_ratio`.

The P7 comment on `relref`, before the move:

P7: the CPU's relref, same values, same default; 'sigglobal' on
this backend now splits unit groups (node voltages vs branch
currents) exactly as the CPU does -- the pre-P7 scalar reference
mixed volts and amps.  The COUPLED path keeps its scalar
reference (its measured record depends on it; see fang's note).

### `__init__`

2026-09-27 (moved from the code):

The P17 comment, before the move:

P17 -- and this refusal is the PERMANENT contract, not a "not
yet": `nrsolver`/`scaler`/`linearsolver` are Python strategy
objects dispatched per iteration, and a traced `jax.lax.while_loop`
cannot call into them -- ever; that is what tracing means.  They
used to be accepted in silence (the "thin advertised feature" 0.1c
warns about, the same defect shape as the `lte_formula` knob
removed in 9(f)); `linearsolver` joined the refusal at Phase C --
the traced loop solves with jnp.linalg.solve and a passed solver
object was still being swallowed.

### `_initial_state`

2026-09-27 (moved from the code):

The P12 comment above the shared initial-state methods, before the move:

P12: the CPU's initial-state machinery, SHARED rather than ported --
`_initial_state` and its helpers are pre-loop Python over names and
indices (they build a numpy vector; the chunk loop converts), so the
bound functions run here unchanged, spanning-tree capacitor solve
included.  The old `node.ic` attribute walk this replaces was dead
code posing as a feature: nothing in the package or tests ever SET
node.ic, and a misspelled node name in it could not even be detected.

### `_timestep_max`

2026-09-27 (moved from the code):

The docstring body, before the move:

`Transient` resolves identically: None means tend/50, SPICE's TMAX
default.  Decoupled from `timestep` by owner decision -- the old
None-means-timestep coupling made gentle circuits step-cap-limited,
with no tolerance knob able to move the run.

### `_element_cap`

2026-09-27 (moved from the code):

The docstring, before the move:

Clamp ``dt_max`` to the tightest element step cap (stage 8(d)
parity).  A TLine cannot be resolved by steps as long as its delay --
the CPU measured the observed delay at 2.00x TD when dt = TD, with no
warning -- and this backend never applied the cap at all: standard-path
runs were correct only because the adaptive controller happened to
keep steps small, and the coupled path solves h freely and would walk
straight past TD/2 on any quiet stretch.

### `_opening_step`

2026-09-27 (moved from the code):

The first two paragraphs of the docstring, before the move:

STAGE 9(g), ported from `Transient._opening_step`, which records the
measurement: a run opening at `timestep` -- which is also `dt_max`, the
largest step the controller may ever take -- makes the first step both the
**largest** and the **only unchecked** one, because with no history there is
nothing to difference and no truncation error can be estimated.  Its error
then dominates everything after it.

On the CPU that showed up as a global error of 1.3212e-01 at reltol 1e-3,
1e-4, 1e-5 AND 1e-6 -- identical to five digits -- while the step count went
from 24 to 195.  **Gate 9-1(c) measured the identical signature here**:
4.2535e-3 across the same four decades, always at index 1, step count 53 to
85.  Eight times the work, or 1.6x, for the same answer.

### `_periodic_state_arrays_batched`

2026-09-27 (moved from the code):

The second paragraph of the docstring, before the move:

A swept ``modulus``/``offset`` used to DROP the element from the
shift (correct Phase-1 fallback, unbounded state); now the
declaration is re-evaluated per lane with the lane's parameter
column applied, so swept lanes KEEP the bounded-state property.

### `_pcnr_setup`

2026-09-27 (moved from the code):

The comment in the vector branch, before the move:

⚠ THE CHARGE IS THE LIMIT HERE, not the current.  Vector
PCNR shadows a participant's `i`/`G` out of the ordinary
assembly and re-stamps it at `v_lim` through `pcnr_i`, which
traces -- so the DEVICE's own `i` never runs, and a device
whose `i` cannot be traced is still fine.  Its `q` is not
shadowed: charge stays in the MNA block at the node voltages,
which is the CPU's documented trade, so the transient DOES
call it.  Measured 2026-08-31: no class declaring
`pcnr_probes` has a traceable `q` today.

Probed once, here, rather than left to fail as a
TracerArrayConversionError several frames inside a compiled
chain -- which is what this refusal replaced, and is worse
than the NotImplementedError it replaced in turn.

### `solve`

2026-09-27 (moved from the code):

The P3 comment (above `_tline_setup`), before the move:

P3: `dt_min=1e-15` was a second, differently-named, three-decades-
looser floor for the same physical quantity `solve` calls `minstep`
-- both entry points now read the one `minstep` Parameter (per-call
override under the same name).  The first argument accepted a node
object and called get_node_index on it, so its old name `irefnode`
lied; renamed `refnode`, defaulting to the same object as everywhere
else (P4).

The comment above the chunk state, before the move:

`sig_max` and `n_rejected` are RUNNING TOTALS and must cross the chunk
boundary.  Rebuilding the state without them reset the `sigglobal`
reference to zero every CHUNK_SIZE steps -- so a long run silently
reverted to a `pointlocal`-like tolerance at each boundary -- and threw
the rejection count away.  Same shape as the `_dt_last2` reset the CPU
side had: a per-run quantity re-seeded by a per-call constructor.

The STAGE 9(e) comment in the chunk loop, before the move:

STAGE 9(e) -- CHECKED PER CHUNK, NOT AT THE END, BECAUSE THE END MAY
NEVER ARRIVE.  Rejecting a non-converged step shrinks `dt`; a circuit
that cannot converge at any step size drives `dt` to `dt_min` and then
advances by `dt_min` forever.  Measured while writing this: a linear RC
at maxiter=1 needs ~5e12 steps to reach tend that way, so warning
"after the run" is a warning that never prints.  That is the trade the
plan flags for 9(d) -- turning a silent wrong answer into a hang -- and
it is avoided by raising here, where the loop is bounded by CHUNK_SIZE.

### `_tline_setup`

2026-09-27 (moved from the code):

The docstring, before the move:

The circuit's TLines for a run starting at `x0s`, one state per
lane: `(params (T, 2) [TD, Z0], indices (T, 6), history (L, T,
depth, 5), tlines)`, `tlines` the `(name, element)` pairs.  Each line's delay history holds its DC state at
every slot -- `t = -1` so nothing collides with `t = 0`, which sits
at head 0 -- so the delayed terms read the operating point until the
run has been going for TD.  `solve` passes one state, `solve_batched`
one per lane (one setup since 2026-09-27; the batched copy left the
history at zero).

### `solve_batched`

2026-09-27 (moved from the code):

The comment above `batchable`, before the move:

An override for a class that is not in the vmap evaluation groups was
SILENTLY IGNORED: `params_tree` is only consumed for classes with
`eval_i_pure`/`eval_q_pure`, so `{'R': {'r': ...}}` produced N
bit-identical lanes presented as N samples -- a parameter sweep that
sweeps nothing, with no symptom (doc/transient_review_260820.md, F2;
measured with r = 1 ohm vs 1 kohm).  Refuse loudly until the class is
made batchable.

The comment on the class-name keying, before the move:

THE TREE IS KEYED BY CLASS NAME -- the vmapped groups are
per-class stacks, and `batched_contributions` consumes
`params_tree[cls.__name__]`.  Every test happened to name its
instance after its class (`c['R'] = R(...)`), which hid this from
users until an instance named 'R1' was refused with a message
about classes (roadmap item 1, idtmod.md sec. 8).  An INSTANCE
key is remapped when unambiguous -- its class has exactly one
instance -- and refused with the correct spelling otherwise,
because with several instances the override's per-instance column
order is the group order, which the caller must address by class.

The STAGE 9(g) comment on `dt_max`, before the move:

STAGE 9(g) -- `timestep`, not `tend/10`.  `solve` and `solve_batched`
are two entry points of one class and disagreed about this by ~50x,
so the same circuit at the same requested timestep was error-
controlled to two different standards depending on which was called.
`tend/50` matches `solve` and matches the CPU's resolution.

The comment above the delay-line setup, before the move:

⚠ each lane's delay history from ITS OWN starting state, as `solve`
fills it -- until 2026-09-27 it was zeros ("Init with zero for
now"), and a line carrying DC read 0 V for its first TD: measured,
a matched 1 V line started from its operating point came out 0.5 V
off at the first point, where `solve` held it exactly

The PER-LANE COLLECTION comment, before the move:

PER-LANE COLLECTION, NO PADDING.  The old code trimmed every lane
to the batch's rectangular [0, max_steps) window and FILLED
SHORTER LANES FORWARD -- duplicating both the state AND the
timestamp, so shorter lanes' results carried repeated abscissae
(which break interpolation downstream), and a lane that had
already reached tend in a previous chunk had b_len == 0, making
the fill source `x_chunk[b, -1:0, :]` an EMPTY slice: numpy
raised "could not broadcast (0,n) into (max_steps,n)" and any
heterogeneous sweep spanning more than one chunk crashed.

Nothing needs the lanes to be rectangular: chunk-to-chunk state
flows through `final_state`, and the method returns a separate
Result per lane, each with its own time base.  So each lane
simply keeps its own valid slice and the padding -- the only
code that could crash or corrupt -- is deleted rather than
guarded (doc/transient_review_260820.md, F1(c); measured with
c = 1e-9 vs 1e-6 F, CHUNK_SIZE=10).

The comment above the ring-buffer check, before the move:

the delay line's ring buffer, per lane, as `solve` checks it
(never checked here before 2026-09-27: a wrapped buffer reads
a previous lap's entries with no other symptom)

The STAGE 9(e) comment in the chunk loop, before the move:

STAGE 9(e) ON THE BATCH, per chunk as `solve` checks it: a lane
whose Newton failed even at dt_min stops itself (`time_cond`'s
`alive`), and until 2026-09-27 came back as a waveform cut
short of `tend` with NO error -- measured, 2 points ending at
t = 1e-18 of 5e-3 on the linear RC where `solve` raises.


## `integrator.py` -- module level

### `ZERO_STABILITY_RATIO`

2026-09-27 (moved from the code):

The comment above `ZERO_STABILITY_RATIO`, before the move:

STAGE 4e -- the zero-stability bound on the step-size ratio for variable-step
BDF-2.  Grigorieff's result is that the homogeneous recursion's parasitic root
stays inside the unit disc only while `h_n/h_{n-1} < 1 + sqrt(2)`:

    ratio    parasitic root   growth over 20 steps
    2.414214       1.000000                      1     <- the bound
    2.5            1.041667                  2.262
    3.0            1.285714                  152.4
    10.0           4.761905               3.59e+13

It bounds *growth* only.  Shrinking is unconditionally zero-stable, which is
the whole content of 4e: the guard this constant now protects used to fire on
the shrink and leave the growth unwatched.

### (the `lte_formula` note)

2026-09-27 (moved from the code):

The module note "WHY THERE IS NO `lte_formula` PARAMETER", before the move:

WHY THERE IS NO `lte_formula` PARAMETER.  Removed in stage 9(f), 2026-07-31.

It chose between the classic divided-difference estimates and the
Yao-Wang-Roychowdhury Table I residuals.  Three changes removed its effect on
this backend before it was removed as API:

  4g(b)  the trapezoidal estimator stopped differencing `g` (the companion
         current), which carries an undamped (-1)^n mode;
  4i     both second-order estimators moved onto a shared third-derivative
         estimate taken from a divided difference of the CHARGE, which reads
         neither formula;
  4d     the one-step fallback -- the last place either branch ran -- now takes
         the divided-difference form unconditionally, because YWR's TRAP entry
         is a uniform-grid formula and its GEAR2 residual is 3/4 of the true
         truncation error.

So `'ywr'` and `'classic'` produced bit-identical runs for every integrator, and
the parameter was kept for a while, accepted and documented as inert.  What
settled its removal was the OTHER backend: `jaxtransient.py` carried its own
`lte_formula`, where `'classic'` selected a charge-domain estimator whose
tolerance applied `lte_abs = 1e-6` -- a VOLTAGE floor -- to a CHARGE.  One
microcoulomb, against node charges of pico- to femtocoulombs, so the normalized
error could never reach 1 and no step was ever rejected: the controller ran
open-loop.  One parameter name meant "selects nothing" here and "selects a
broken estimator" there, which is worse than either alone.

Both are gone.  Each backend now has exactly one estimator, and passing
`lte_formula=` raises TypeError rather than being silently ignored -- a kwarg
accepted and discarded is how the JAX defect stayed invisible.


## `integrator.py` -- `Integrator`

### (capability queries)

2026-09-27 (moved from the code):

The comment above the capability queries, before the move:

--- CAPABILITY QUERIES (polymorphic dispatch) ---
The shooting/transient stacks used to branch on `isinstance(...)` and
method-name strings at ~35 sites; these let a caller ask the METHOD
instead, so a new integrator arrives with the right answers rather than
needing an edit at each site.  See doc/integrator_architecture_260906.md.

### `compute_lte`

2026-09-27 (moved from the code):

The UNITS paragraph of the docstring, before the move:

The controller consumes it by mapping it through ``J^-1`` into the solution
domain and comparing *that* against ``reltol``/``vabstol``/``iabstol``, which is
dimensionally sound.  Comparing the raw return value against a charge tolerance
is not: decision 0.3d's option (D) was designed around exactly that and was
refuted on it, having reproduced the units defect gate 0.2b recorded on the JAX
backend.  If a future estimator needs a charge, it is ``h`` times this.


## `integrator.py` -- `EulerIntegrator`

### `__init__`

2026-09-27 (moved from the code):

The comment in `__init__`, before the move:

No `lte_formula`: see the module note above.  For
Backward Euler the two formulas always coincided, so this class never
had a choice to preserve across an order drop in the first place.

### `compute_lte`

2026-09-27 (moved from the code):

The STAGE 4c comment, before the move:

STAGE 4c -- THE VARIABLE-STEP CORRECTION.

Backward Euler's companion current is `(q_n - q_{n-1})/h`, which is a
centred approximation of `q'` at the MIDPOINT of the step, not at the
node.  So `g_n - g_{n-1}` differences two midpoint derivatives separated
by `(h_curr + h_last)/2`, not by `h_curr`, and it therefore estimates
`((h1+h2)/2) q''` where the truncation error is `(h1/2) q''`.

On a uniform grid the two coincide and the estimator is exact, which is
why this went unnoticed.  Off it, the estimate is wrong by
`(h1+h2)/(2 h1)` -- measured est/true, before this correction:

    ratio  0.25    0.5     1.0     2.0     4.0
    est/true  2.5246  1.5089  1.0040  0.7522  0.6265

which is a 4.03x spread across the sweep, and it is the wrong direction
twice over: on a shrinking step the error is OVERstated (so the
controller shrinks further than it needs to) and on a growing step it is
UNDERstated (so the controller grows past what the tolerance allows).

Rescaling by `2 h1 / (h1 + h2)` converts the midpoint spacing back to the
step.  After: 1.0098 / 1.0059 / 1.0040 / 1.0029 / 1.0024 -- flat to
within 1% across the same sweep.


## `integrator.py` -- `TrapezoidalIntegrator`

### `compute_lte`

2026-09-27 (moved from the code):

The STAGE 4g(b) comment, before the move:

STAGE 4g(b) -- DIFFERENCE A MODE-FREE QUANTITY.

What the estimator must produce, from YWR eq (22) with p=1, k=2 and the
trapezoidal coefficients (alpha = [1/h, -1/h], beta = [1/2, 1/2]), is

    Eg = -(h^2/6) q'''

once the controller's own J^-1 has absorbed the (q_x + 0.5 h f_x)^-1
factor.  Table I approximates q''' by a second difference of the
companion current g; eq (22) does not require that, and the paper's own
wording is "(22) AND FINITE DIFFERENCE APPROXIMATION" -- the choice of
difference is free, and for TRAP the g-based choice is what goes wrong.

g carries an undamped parasitic mode.  The trapezoidal companion
`iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}` has homogeneous solution
`iq_n = -iq_{n-1}`, i.e. `(-1)^n`, and nothing damps it.  Differencing g
therefore differences that mode, and the estimate depends on the step
history that preceded it: measured est/true_local at h=1e-9, ratio 1, on
two different prefixes of the SAME problem, 1.3176 and 0.6780 -- a 1.9x
swing from history alone.

`d` has no such component.  The trapezoidal relation gives
`(iq_n + iq_{n-1})/2 = (q_n - q_{n-1})/h = d_n`, and the mode flips sign
every step, so it CANCELS EXACTLY in that sum.  Expanding about the
interval midpoint `m_n = (t_n + t_{n-1})/2`,

    d_n = q'(m_n) + (h_n^2/24) q'''(m_n) + O(h^4)

so d samples q' at midpoints and a second divided difference of d OVER
THE MIDPOINTS estimates q'''/2.  Measured against the local truncation
error, est/true at ratio 1: 0.9273 / 0.9933 / 0.9993 as h falls
1e-8 -> 1e-10, against the g-based form's 0.8067 / 0.6780 / 0.6678 --
asymptotically exact where the old one holds a 33% underestimate.

The midpoint spacings are the part that must not be got wrong; using
h_curr/h_last (the NODE spacings) is the same class of error as the
backward-Euler defect stage 4c fixed.

The STAGE 4i comment in the `h_last2` branch, before the move:

STAGE 4i.  Eg = -(h^2/6) q''' and the shared helper returns q'''/6,
so the whole formula is -h^2 times it.

4g(b) differenced `d_k = (q_k - q_{k-1})/h_k` over the interval
midpoints instead, which removed the parasitic mode and was
asymptotically exact but kept a +-12% bias at the extremes of the
reachable step-ratio range: `d_k` carries `(h_k^2/24) q'''`, and that
term cancels only on a uniform grid.  The charge carries no method
error at all, so differencing it directly removes the residual --
measured spread over ratio 0.008..2.414 falls from 1.26x to 1.008x.

The comment above the one-step fallback, before the move:

`h_last2 is None` means q_last[2] is not yet a real past point, which is
true for exactly ONE step of a run -- the second, where the ring buffer
still holds the initial charge twice.  Falling back to the g-based form
for that single step is better than returning zeros, which would make it
an unchecked step of the kind stage 3 exists to remove.
THE DIVIDED-DIFFERENCE FORM, AND `lte_formula` DOES NOT SELECT HERE.
YWR's Table I TRAP entry -- `Eg = -(1/6)(g_n - 2 g_{n-1} + g_{n-2})` --
carries a single `h` and an UNWEIGHTED second difference, i.e. it is a
uniform-grid formula, where the same table's GEAR2 entry carries h1 and
h2 explicitly.  Off a uniform grid it is wrong by O(1/h), and the grid is
not uniform here: stage 3's opening ramp is still growing the step at the
one point this fallback runs.  Taking the divided-difference form
unconditionally is stage 4d's stated fix.


## `integrator.py` -- `ThetaIntegrator`

### `needs_consistent_iq0`

2026-09-27 (moved from the code):

The comment in `needs_consistent_iq0`, before the move:

⚠ THE PREREQUISITE THE LINEAR GATE COULD NOT SEE.  B2's gate solved
the periodic state in closed form, so it never had an opening step
and never exercised the seed.  Refusing the Euler opener means this
method reads `iq_{-1}` on its first step, where the run seeds zero:
measured, that alone costs a full order (0.97 against 2.00).

⚠⚠ AND IT IS A JACOBIAN STATEMENT AS WELL AS AN ACCURACY ONE.  The
seed `-(i(x_0) + u(t_0))` is a FUNCTION OF `x_0`, so any Jacobian
taken with respect to `x_0` -- the shooting monodromy above all --
carries `d(iq_0)/d(x_0) = -G(x_0)`.  `shooting.py` seeded that at
ZERO for every method (nothing before this one formed a companion
current at `x_0`), which ANNIHILATED `null(C)` in the monodromy:
the L-stable opener this class exists to remove, reintroduced in
the derivative.  Cost 99 shooting evaluations against trap's 3 on a
LINEAR circuit, where an exact Newton must land in one step; the
answer was unchanged throughout.  See `PSS._pq_seed_at_x0`.  A NEW
METHOD THAT RETURNS TRUE HERE INHERITS THAT FIX AND NEEDS NOTHING.


## `integrator.py` -- `Gear2Integrator`

### `__init__`

2026-09-27 (moved from the code):

The comment in `__init__`, before the move:

No `lte_formula`: removed in 9(f) -- see the module note above.
The history is kept because both entries are recorded results, and both
are about a choice that no longer exists rather than about this class.

THE DEFAULT USED TO BE 'ywr', chosen belt-and-braces when 'classic' was
repaired, on the grounds that it had the longer track record.  The price
was that the YWR GEAR2 residual estimates (1/4) h^2 q''' against a true
(1/3) h^2 q''', so it reported 3/4 of the truncation error at every step
where a corrected 'classic' is asymptotically exact.  Stage 4i moved both
variants onto the divided-difference form, which is the 'classic' one, so
that optimism is gone rather than merely defaulted around.

THE "5/6" AN EARLIER COMMENT CLAIMED WAS AN ARTEFACT, worth keeping
because the number is so clean.  5/6 = 0.8333 is what the trapezoidal
estimator reads when it is handed EXACT derivatives as its `g` history
instead of the companion currents a real run produces.  Measured against
the local truncation error with the real history it reads 1.09 / 1.31 /
1.33 as h falls 1e-8 -> 1e-10 -- it does not converge at all, let alone
to 5/6.  Decision 0.3b called the claim measurably wrong; this is the
measurement, and the mechanism.

### `get_required_history`

2026-09-27 (moved from the code):

The comment in `get_required_history`, before the move:

THREE since stage 4i, for the same reason trapezoidal needs three: the
METHOD looks back two steps -- `compute_derivatives` uses q_last[0] and
q_last[1] -- but the ESTIMATOR takes a third divided difference of the
charge and so needs q_{n-3}.  Until 4g(b) built the `h_last2` plumbing
this was not available, and the comment in `compute_lte` below recorded
it as the reason the g-based form had to be used.

### `check_order_drop`

2026-09-27 (moved from the code):

The comment above the growth guard, before the move:

STAGE 4e -- THE GUARD USED TO WATCH THE WRONG DIRECTION.

The only test here was `if h_curr / h_last < 0.1`, and it was labelled
as protecting the validity of the high-order polynomial.  It does not:
variable-step BDF-2's parasitic root leaves the unit disc only on
*growth*, past `ZERO_STABILITY_RATIO`, and any ratio below 1 is
unconditionally zero-stable.  So the one ratio that can actually
destabilise the recursion was unwatched -- which is how the 10x growth
on `transient.py`'s force-accept path (4b) survived: nothing downstream
would have caught it.  (The shrink test itself is kept, for a different
and measured reason; see below.)  Measured on the stiff RLC, this new
branch fires 0 times, because the controller's own clamp
(`MAX_GROWTH_RATIO` = 2.0) keeps every normal step inside the bound.

**That makes this a backstop, and a backstop that never fires in a
healthy run is the point of it** -- it is what turns "no accepted step
ratio exceeds the bound" from an accident of two clamps agreeing into
something the integrator enforces for itself.  Dropping to Euler is the
right response rather than refusing the step: order 1 has no parasitic
root to amplify, so the ratio becomes harmless instead of forbidden.

The comment above the shrink branch, before the move:

THE SHRINK BRANCH IS KEPT, AND RE-LABELLED.  The plan said replace; the
measurement said add, so it is added and the reason is written down.

Removing it outright took `Gear2('ywr')` and `Gear2('classic')` from 0
force-accepts to 1 each on the stiff RLC at reltol 1e-5, because it is
not idle: it fires 3-6 times a run there.  What it is doing is nothing
to do with zero-stability -- it is a STALLED-ESTIMATE heuristic.  A step
only shrinks 10x below the last accepted one after several consecutive
rejections, and what rejects repeatedly is a 2nd-order estimate built on
a third difference of a solution that is not three times differentiable
-- i.e. a discontinuity.  Dropping to order 1 there is what every
simulator does across a corner, and it is the same medicine 4b now
administers at the rejection cap, one retry later and after having
accepted an over-tolerance step to get there.  Deleting it would have
traded a controlled Euler step for a force-accepted 2nd-order one.

So: the guard above is the stability bound, this one is economics, and
the defect 4e names was never that this branch existed -- it was that
this branch was ALL there was, and it was labelled as protecting a
stability property it has nothing to do with.  **Reconsider if** the
rejection cap ever becomes rejection-count-aware: `h_curr/h_last < 0.1`
is a proxy for "we have rejected three times at this time point", and
the thing it is a proxy for is known exactly one level up in
`transient.py`, where it would not need a threshold at all.

### `compute_lte`

2026-09-27 (moved from the code):

The STAGE 4i comment, before the move:

STAGE 4i -- THE ESTIMATOR USED TO DIFFERENCE THE METHOD'S OWN ERROR.

Gear-2's local truncation error is `-(1/6) h1 (h1+h2) q'''`, so every
companion current in the history carries an error of exactly that shape.
Both branches below take a second divided difference of `g` at the
nodes, which differences those errors along with the signal.  The
damage was computed by hand before it was measured, and the two agree to
0.3% at every step ratio:

    h1/h2      0.008    0.05     0.1    0.25       1       2       4
    predicted  83.34  13.365   6.727   2.800  1.0000   0.778   0.700
    measured   83.06   13.32    6.71    2.79   0.998   0.775   0.695

It vanishes exactly at h1 = h2 = h3, which is why the estimator measured
asymptotically exact (1.000282 against 2/9) on a uniform grid: that
measurement was taken at the one ratio where the defect is zero.  A step
ratio of 0.008 is reached after three consecutive rejections, so the
worst case is not hypothetical -- it is the step where the controller
has just collapsed the step size and is told the error is 83x worse than
it is.

The fix is to estimate q''' from the CHARGES, which carry no method
error.  The obstacle used to be real and was recorded here: Gear-2 kept
only two past charges, so a third divided difference was unavailable.
4g(b) lifted it.

The comment above the one-step fallback, before the move:

`h_last2 is None` for exactly one step of a run -- the second, before the
ring buffer holds four real charges.  The g-based form below serves that
step; returning zeros would make it unchecked, which is the defect stage
3 removed from the first step.

IT IS THE DIVIDED-DIFFERENCE FORM, NOT YWR's, AND `lte_formula` DOES NOT
SELECT HERE.  YWR's Table I GEAR2 residual estimates `(1/4) h^2 q'''`
against a true `(1/3)`, so it reports 3/4 of the truncation error --
measured on this exact fallback as -2.827659e+01 where the correct value
is -3.770212e+01.  After 4i this was the ONLY step of a run where the
choice still had any effect on the CPU path, so taking the accurate one
unconditionally is what finishes 4d: "delete the branch and keep 'ywr' as
an alias".
(The YWR Table I GEAR2 residual that used to be selectable here,
 `Eg = -(1/8)((h1+h2)/(h1 h2))(h2 g_n - (h1+h2) g_{n-1} + h1 g_{n-2})`,
 is derived and compared against this one in doc/src/circuit/lte_dae.rst;
 it is not kept as dead code.)

The CLASSIC GEAR-2 comment, before the move:

--- CLASSIC GEAR-2 LOCAL TRUNCATION ERROR ---
Taylor-expanding the VSS companion current above about t_n (the alpha
coefficients kill the q' and q'' terms by construction) leaves

    iq - q'(t_n) = -(1/6) h1 (h1 + h2) q'''(t_n) + O(h^3)

-- equal steps: -(1/3) h^2 q''', the textbook BDF-2 result.  So what
has to be estimated here is the THIRD derivative of the charge, scaled
by h^2.

A second divided difference of q yields only q'', and Gear-2 keeps just
two past charges (get_required_history() == 2), so a third divided
difference of q is not available at all.  The third derivative is
therefore taken as the second divided difference of g = dq/dt, read off
the companion-current history -- the same information the YWR branch
above uses.  Estimating q'' here and multiplying by h^3 (as this branch
did until 2026-07) is dimensionally not a current: it undershoots the
truncation error by a factor of order h*omega, which on a 1 MHz signal
at nanosecond steps is ~1e-15.  The controller then never rejects a
step, saturates the growth limiter every step, pins h at max_step and
stops responding to reltol/abstol altogether.


## `integrator.py` -- `RadauIIA3Integrator`

### (class docstring)

2026-09-27 (moved from the code):

The stability paragraphs of the class docstring, before the move:

⚠⚠ EVERY STABILITY LINE ABOVE IS A PROPERTY OF THE TABLEAU, AND ON A DAE
THAT IS ONLY HALF THE PREMISE.  Lamour, Marz & Tischendorf (2013)
Example 5.1: applying the IMPLICIT Euler method to a DAE in STANDARD FORM
induces the EXPLICIT Euler method on the inner variable -- "the
A-stability gets lost when the method is applied to a DAE in standard
form", unstable at ``h = 0.0202`` against a limit of ``0.02`` at
``lambda = -100``.  The formulation, not the method, decides whether
A-stability transfers.

This tree is on the right side of that because charge-oriented MNA IS a
properly stated leading term -- but that premise was UNSTATED here until
2026-09-10, and an unstated premise is an unprotected one: nothing would
notice if a formulation change moved it.  The condition it turns on is
``im D(t)`` time-invariant, which for charge-oriented MNA is ``im C(x)``
constant along the orbit; it moves only on a RANK change (a switch, a
device leaving conduction), NOT on a smoothly varying ``C(v) > 0``.

⚠ THE SAME TERM GOVERNS THE GLM CONVERGENCE RESULT.  Prop 4.7 and the
proof of Thm 5.7 make it one term in one equation: the IERODE's field is
``u' = R'(t)u + D(t)omega(u,t)``, and the hypothesis exists to kill
``R'(t)u``.  It is the same premise under IRK(DAE) convergence (5.7), GLM
convergence at stage order (5.9, which is what
:class:`NordsieckGLMIntegrator` rests on) and contractivity transfer
(6.9).  Relayed from a source reading, not verified here.

⚠⚠ AND OUR FIXTURES DO NOT EXERCISE IT.  Measured 2026-09-10: on the
index-2 C-V loop, the state-free exponential and the van der Pol -- the
three fixtures behind the GLM order, stage-predictor and noise-floor
results -- ``C(x)`` is not merely constant-rank but LITERALLY CONSTANT
(``max|C(x1) - C(x2)| == 0`` over random ``x``), because every reactance
in them is linear.  So ``im D`` is time-invariant TRIVIALLY and those
results sit inside the theorem's scope VACUOUSLY.  Whether violating the
hypothesis actually costs order in this implementation is NOT MEASURED and
would need a fixture whose ``rank C(x)`` genuinely changes along the orbit.

What it does NOT have (2026-09-07, measured after the method became the
``PSS`` default): a radius of absolute monotonicity.  Kraaijevanger's
``R(A, b)`` -- computed here from the definition, on the SAME script that
returns ``2`` for Crank-Nicolson, ``inf`` for backward Euler and
``1 + sqrt(2)`` for TR-BDF2 at every ``gamma`` against Bonaventura &
Della Rocca's closed form -- is ``R(A, b) = 0`` for Radau IIA(3), and
for Radau IIA(2) too.  The cause is one tableau entry: ``a_23 =
-2/225 - sqrt(6)/75 < 0``, and ``R > 0`` requires ``A >= 0`` entrywise
(Kraaijevanger 1991), so NO step size, however small, carries a
monotonicity / positivity / TVD guarantee under this method; TR-BDF2's
``~21%`` margin over trapezoidal has no Radau counterpart at all.  This
is the theorem, not an accident of the tableau: Kraaijevanger's order
barrier for unconditional contractivity is ``p <= 1`` (the "Radau IIA"
named there as unconditionally contractive is the ONE-stage member,
i.e. backward Euler), and conditional contractivity at ``p = 5`` needs
``A >= 0``, which the collocation tableau does not give.  Practical
reading: the order win that reaches 1 ppb at 30-80 points per period
is bought with LARGE steps, and large steps are exactly where a method
with ``R = 0`` may overshoot or go negative on a stiff switching
circuit -- nothing in this tree has yet MEASURED such a failure, so
this is a recorded absence of a guarantee, not a recorded failure.
``ESDIRK43`` is in the same position (``R = 0``, its ``a_32 = -1743/31250``
and ``min(A) = -0.59``); computed on the coded ``A``/``B`` of all three
classes, TR-BDF2 is the ONLY stage method in this file with a positive
radius (``2.41421``, the closed form to the digit).
AND THE ZERO IS STRUCTURAL, NOT AN ACCIDENT OF THE TABLEAU (peer
reading of the same paper, same day): Thm 8.5 of Kraaijevanger 1991,
which is BUTCHER's result (headed "J. C. Butcher; private communication
1989" -- checked at the source 2026-09-08, on disk as
`Kraaijevanger-1991-Contractivity of Runge-Kutta methods.pdf`): "Let
(A,b) be an ARBITRARY coefficient scheme with A >= 0. Then the stage
order is at most 2. Further, if it equals 2 then A has a zero row" --
no irreducibility needed, and the proof says WHICH row: c1 = 0, an
EXPLICIT FIRST STAGE, the only shape A >= 0 permits at stage order 2.
``A >= 0`` is NECESSARY for ``R > 0`` (Thm 4.2, Kraaijevanger's own,
verbatim "for irreducible coefficient schemes ... R(A,b) > 0 if and
only if A >= 0, b > 0 and Inc(A^2) <= Inc(A)"), so a positive radius
forces stage order <= 2 AND an explicit first stage for EVERY
Runge-Kutta method; Radau IIA(3) has stage order 3, so its
negative entry is the theorem showing its face and no better
collocation tableau exists to look for.  Chained with Voigtmann's
Theorem 5 (index-2 convergence order = min(p, q)): on an index-2
circuit a method can carry an ABSOLUTE-MONOTONICITY guarantee OR order
above 2, never both -- a Runge-Kutta limitation, not a DIRK one, binding
Radau exactly as hard as ESDIRK.  And the split has its theorem
(Voigtmann's thesis, Thm 9.5, verified 2026-09-08): a GLM with stage
order q = p, nilpotent M_inf, stiffly accurate, keeps full order p on
an index-2 DAE at CONSTANT stepsize -- Radau IIA has q = s, p = 2s-1,
so it can never meet q = p (measured 5 / 3.05 here); TR-BDF2 (q = p =
2) meets it and shows no split (2.04).  The reduction is a property of
methods with q < p, not of index-2 as such.
⚠⚠ BUT "CONTRACTIVITY" IS TWO DIFFERENT GUARANTEES, and this method
HAS the other one unconditionally (peer correction, same night, from
Hairer & Wanner Thm 12.9 p.210 -- on disk all along -- and Lamour Ch.6
§6.1): "The methods Gauss, Radau IA, Radau IIA and Lobatto IIIC are
algebraically stable and therefore also B-stable."  B-stability is
contraction of a one-sided-Lipschitz flow in an inner-product norm,
with NO stepsize restriction ("B-stable Runge-Kutta methods reflect
contractivity devoid of stepsize restrictions", Lamour) -- the natural
notion for a DISSIPATIVE circuit of passive elements.  Absolute
monotonicity (Kraaijevanger, SSP) protects componentwise positivity /
TVD / bounds and needs ``A >= 0``, hence stage order <= 2.  Each method
here has exactly one: Radau IIA(3) is B-stable with ``R(A, b) = 0``;
TR-BDF2 has ``R = 1 + sqrt(2)`` and is NOT algebraically stable (Kennedy
& Carpenter: a DIRK may have algebraic stability OR stage order two,
not both).  So the default is BETTER protected than the paragraph above
reads, on the notion that matters most for circuits; the ``R = 0`` gap
is the narrower componentwise one.  A componentwise-bounds measurement
(an RC ladder under a square wave) tests ONLY absolute monotonicity: a
violation there does not contradict B-stability, and a clean result
does not establish it.
⚠⚠ AND B-STABILITY IS AN ODE PROPERTY THAT DOES NOT CARRY TO A DAE FOR
FREE (peer qualification of the correction above, same night, Lamour
Ch.6 §6.2 verbatim): "in general, we cannot expect that algebraically
stable Runge-Kutta methods, in particular the implicit Euler method,
preserve the decay behavior of the exact DAE solution without strong
stepsize restrictions, not even when we restrict the class of DAEs to
linear ones.  This depends on how the DAE is formulated."  The
condition (their eq 6.16): a properly stated leading term whose
``im D(t)`` -- the IMAGE SPACE of the charge Jacobian, not ``D`` itself
-- is independent of ``t``; then the IRK reaches the inherent ODE
unchanged and algebraic stability + a contractive DAE give contractivity
with no step restriction.  So the default's protection is a THREE-PART
CONDITIONAL: B-stable (yes, by theorem) + contractive DAE (a property of
the circuit) + constant ``im D(t)`` (a property of the formulation and
the circuit -- a smooth nonlinear capacitance of constant rank is fine;
a switch, or a device entering/leaving a region where it contributes a
state, is what breaks it).  Two of the three are unverified for every
fixture in this tree.  ⚠ This tree's ``Diode`` has ``G``/``i`` only, no
charge, so a diode peak detector keeps ``im D`` constant by
construction and is NOT the structure-change fixture; a ``VSwitch`` in
series with a capacitor would be.  TR-BDF2 sits at the corner (stage order 2,
``R = 1 + sqrt(2)``) and its explicit first stage is the zero row the
theorem requires -- checked on the coded ``A``: row 0 is the only zero
row and ``A >= 0`` holds; Radau IIA(3) has no zero row and ``A >= 0``
fails; ESDIRK43 has the explicit stage (clears the NECESSARY condition)
and fails ``A >= 0`` on ``a_32 < 0``, so ``R = 0`` anyway.  The check
measured Butcher's construction before either reader knew what it was
testing: the "zero-row tell" is an explicit-first-stage tell.  The one unexplored exit is the GLM class (Voigtmann: "diagonally
implicit methods with high stage order are possible"), where
Kraaijevanger's RK theorems do not bind -- whether a GLM can carry high
stage order AND a positive radius is OPEN here, not hinted.


## `integrator.py` -- `NordsieckGLMIntegrator`

### (class docstring)

2026-09-27 (moved from the code):

The SCOPE paragraph of the class docstring, before the move:

⚠ SCOPE, and this paragraph has been WRONG TWICE, so read the dates.  The
shooting period map on the `r*m` Nordsieck state IS built (`PSS(method=
'glm3')` and friends), and so is ADAPTIVE STEPPING (2026-09-10): the
estimate is the change in the top Nordsieck component,
`Transient._glm_error_estimate`, delivered through the same `_rk_est` slot
the Runge-Kutta methods use, with `EMBEDDED_ORDER = p`.

### `carries_own_monodromy`

2026-09-27 (moved from the code):

The docstring, before the move:

False -- and NOT because the map is inaccurate.  This method's own
period map acts on the NORDSIECK state (width `r*m`), and every
state-space consumer (`ppv`, `floquet_modes`, `oscillator_covariance`,
the PAC/pnoise family) wants a map on `x`.  ⚠⚠ THE TWO ARE NOT RELATED
BY TAKING THE FIRST BLOCK: a state kick `dx` at `t = 0` perturbs the
HIGHER Nordsieck components too (`dQ_k = h^k d^k(C dx)/dt^k`), so
neither `w[:m]` nor `C^T w_0` is the state PPV -- MEASURED against
radau on van der Pol at Q = 15.9, both are wrong by 2 % in norm and by
13x in the small component, and `ppv` used to return the first of them
SILENTLY.

Since 2026-09-24/25 the map ON THE STATE is built (the startup
linearised, `_GLMStartup`; `_GLMPeriod.state_map`), and `ppv`,
`floquet_modes`, PAC, the adjoint rows, pnoise and `sampled_noise`
read it (`PSS._state_map`).  Still False, because two things keep a
GLM from BEING a twin, both measured: the startup that opens the map
breaks the discrete phase symmetry (the unit multiplier sits
``O(h^p)`` off 1 -- an oscillator's small-signal response near a
harmonic is therefore read with the pole carried analytically,
`PAC._deflated_solve`), and its white-noise injection is one shared
sample per step, first order (`PAC._lyapunov_pieces_glm`).  The
covariance surfaces take a radau twin by default
(`PSS._state_twin`); `monodromy='native'` keeps the GLM's.


## `integrator.py` -- `GLM4Integrator`

### (class docstring)

2026-09-27 (moved from the code):

The last three paragraphs of the class docstring, before the move:

⚠⚠ BUT AT EQUAL WORK RADAU DOMINATES, and the "one factorisation per
step" argument for this class is REFUTED on that fixture: timed, glm4
takes 2.4x-4.7x radau's wall clock at the SAME grid.  `s = r = p + 1`
is five sequential `m x m` stage solves against radau's ONE coupled
`3m` solve, and at `m = 3` a dense 9x9 factorisation is trivial while
five Newtons are not.  ⚠ THAT EXPLANATION IS WRONG AND IS WITHDRAWN:
dense LU is ~n^3/3, so five `m` solves against one `3m` solve is the
same ratio at EVERY `m` -- there is no size threshold, and on flops
alone glm4 should have been cheaper at m = 3.  MEASURED mechanism
(device-evaluation counts, 40 points): glm4 116.2 `i` calls per step
against radau's 39.0, i.e. 3.0x the nonlinear work -- more than its
stage-count ratio of 5/3, because each of its five stages runs its OWN
Newton to convergence while radau solves three stages as one coupled
system.  That ratio does not shrink with circuit size and dominates on
a circuit with expensive compact models.  (The timing was fair: radau's
opt-in cost transform was OFF, so it paid the full dense 3m solve and
glm4 still lost; with it ON radau is ~five real m-solve equivalents,
factorisation parity.)  The lever, unbuilt, is the STAGE INITIAL GUESS:
23 `i` evaluations per stage says the per-stage Newton starts far from
its root, and the Nordsieck vector carries the scaled derivatives a
Taylor predictor would need.  So this method's claim is accuracy per
STEP on the algebraic component, not accuracy per second.

⚠ THE ABSCISSAE ARE INSIDE [0, 1] BY CONSTRUCTION -- c = [0.173, 0.402,
0.641, 0.794, 1] -- which matters more for a circuit than for a general
ODE code (docs session): a stage at `c > 1` evaluates the device models
PAST the end of the step, further from the operating point the Jacobian
was formed at, and samples a time-dependent source in a DIFFERENT PWL or
pulse SEGMENT rather than extrapolating the same one.  An earlier
tableau had c = 1.02 and 1.219; bounding `c` costs nothing theoretically
(it is a free parameter, Wright p. 80, and his canonical choice is
inside) and, MEASURED, nothing in practice either -- the bounded solve
converges in the same ten minutes as the unbounded one.

⚠⚠ STILL BADLY SCALED, AND SHIPPED WITH THAT STATED: `B`'s last rows
carry entries of order 3900 where GLM3's are O(1).  That is a DIFFERENT
problem from the abscissae and a magnitude penalty bolted onto the
feasibility solve does not fix it (measured: the two constraints
together do not converge in 65 minutes where either alone takes ten).
Wright §3.11 -- *"methods where all coefficients have magnitude less
than or equal to one CAN BE FOUND for methods of high order"*, findable
because the construction is linear in the free parameters -- is the
machinery for it, not more restarts.  Measure before preferring this to
GLM3 on a circuit with sharp sources; see the roadmap.

## `stepcontroller.py` -- `StepController`

### `StepController.relref`

2026-10-01 (moved from the code, review C21):

The comment above `relref` before the move:

ITEM 2+.3 -- what the RELATIVE part of the LTE tolerance is measured against.

The tolerance is `lteratio * (reltol*ref + abstol)`.  Until now `ref` was
hard-coded to `max(|x_curr|, |x_last|)` -- each unknown against itself, at
this instant.  That is a commercial simulator's `pointlocal`, and it has a failure mode that
is easy to miss: on a node carrying no signal `ref -> 0`, so the tolerance
collapses to `abstol` and the controller starts chasing numerical noise on a
quiet node.  On the leapfrog that alone cut the step size 5.4x, and the fix
applied at the time was to raise the absolute floor a millionfold
(`lte_vabstol` 1e-12 -> 1e-6), which treats the symptom.

A commercial simulator's answer is `relref`, and its default is `sigglobal`: measure each
signal against the largest signal anywhere in the circuit, over all past
time, so a quiet node inherits a sane reference instead of degenerating.

  pointlocal  each unknown against itself, now.  (pycircuit's historical
              behaviour, and still selectable.)
  alllocal    each unknown against its OWN largest value so far.
  sigglobal   each unknown against the largest value of ANY unknown so far.

THE WORKAROUND ABOVE IS NOW GONE.  With `sigglobal` shipped, `lte_vabstol` is
back to 1e-12: measured at gate D3-e, 1e-6 / 1e-9 / 1e-12 give bit-identical
runs under `sigglobal` (403 steps on a pulsed RC, 601 with a quiet node, at
every value), where under `pointlocal` the same change costs 8.5-9.2%.  That
difference IS the symptom, and it is what the floor was raised to hide.

DEFAULT IS `sigglobal` SINCE DECISION D3's SECOND ATTEMPT, matching a commercial simulator.
It was adopted, sent back by its own gate, and re-run once the reason for the
failure was removed -- see the D3 gates in `doc/transient_work_plan.md`.

### `StepController.lte_gamma_min`

2026-10-01 (moved from the code, review C21):

The paragraph on the band's defaults before the move:

THE DEFAULTS BELOW REPRODUCE THE PREVIOUS BEHAVIOUR EXACTLY: `gamma_min=0`
makes the lower test vacuous and `gamma_max=1` is the historical `err > 1`
rejection.  That is deliberate -- stage 12 is behind a flag until 12D, and
a band that changed the default path would make its own gate unreadable.

## `transient.py` -- `Transient` (review C10-C12, 2026-10-01)

### class docstring

2026-10-01 (moved from the code, review C10/C11):

The CPU-only sentence before it was corrected (the JAX backend has a
variable-step trapezoidal estimator; what it lacks is every stage and
multivalue method):

CPU-only, with cause: trapezoidal integration (a correct VARIABLE-step
trap estimator exists only here) and the
`nrsolver`/`scaler`/`linearsolver` strategy objects -- ...

The backward-Euler sketch and the two examples, before the move (neither
example had ever run: `test_doctests.py` did not collect `transient.py`.
The RC one started from the DC operating point, where the capacitor
already sits at 9.90 V, so it read 9.90 against "6.3"; with `uic=True` it
reads 6.29.  The RLC one compared one instant of a 0.07 V sinusoid with
0.0063 and was deleted):

    i(t) = c*dv/dt
    v(t) = L*di/dt

    The usual companion models are used.
    backward euler:
    i(n+1) = c/dt*(v(n+1) - v(n)) = geq*v(n+1) + Ieq
    v(n+1) = L/dt*(i(n+1) - i(n)) = req*i(n+1) + Veq

    def F(x): return i(x)+Geq(x)*x+u(x)+ueq(x0), G(x)+Geq(x)
    x0=x(n)
    x(n+1) = fsolve(F, x0, fprime=J)

    Linear circuit example:
    >>> tran = Transient(c)
    >>> res = tran.solve(tend=10e-3,timestep=1e-4)
    >>> expected = 6.3
    >>> abs(res.v(n2, gnd)[-1] - expected) < 1e-2*expected #node 2 of last x
    True

    Linear circuit example:
    >>> from pycircuit.circuit.elements import ISin
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> c['Isin'] = ISin(gnd, n1, ia=1e-3, freq=16e3)
    >>> c['R'] = R(n1, gnd, r=200)
    >>> c['C'] = C(n1, gnd, c=1e-6)
    >>> c['L'] = L(n1, gnd, L=1e-4)
    >>> tran = Transient(c)
    >>> res = tran.solve(tend=260e-6,timestep=1e-6)
    >>> expected = 0.063
    >>> abs(res.v(n1,gnd)[-1]) < 1e-1*expected #node 2 of last x
    True

The TODO under the class constants, long done (the step is LTE-adaptive):

TODO:
* Implement automatic timestep adjustment, using difference between
  BE and trapezoidal as a measure of the error.
  Reference: "Time Step Control in Transient Analysis", by SHUBHA VIJAYCHAND

### `Transient._opening_step`

2026-10-01 (review C12): the comment read "An opening step larger than
max_step would be capped on the very next step anyway, and asking for one is
more likely a mistake than an intent." -- the code caps at `timestep`
(`min(firststep, timestep)`), not at the step cap.

## `integrator.py` -- corrected claims (2026-10-01, review C1-C3)

Messages and docstrings the review found WRONG about the current code
(named functions that no longer exist, "not yet built" for what is
built), rewritten on 2026-10-01; the removed text, verbatim:


`integrator.py`:

```
:meth:`companion_coefficients`.  Those three methods therefore raise here.
```

```
'd(iq)/dT shared by the LMMs does not apply; the autonomous '
'shooting dT for a stage method is not yet built.')
```

```
dedicated coupled solve (``_solve_timestep_radau``), not through the
```

```
'loop runs it via _solve_timestep_radau (one coupled 3n Newton '
```

```
'estimate is computed in Transient._solve_timestep_radau and '
'consumed by _run_radau_adaptive.')
```

```
Still not built: no PCNR stage path; not on the JAX backend.
```

## `nrsolver.py` -- `ChordNewton` (2026-10-01)

### (class docstring)

The review's S15 asked for a chord Jacobian on compact models.  Measured
before built: on a compact MOSFET (PSP) common-source stage, gear PSS at 40
points, `G` is 60.7 % and `C` 31.4 % of the solve (one `G` ~19 ms against
`i` ~0.4 ms), and a step evaluates them 3.18 times -- 2.18 Newton
evaluations plus the branch check's `C` and `jacobian_only`'s `G` at the
converged point.  The step's converged-point Jacobian is what the step
controller, the branch check and the shooting's monodromy read, so a chord
keeps that evaluation and holds the iterations' Jacobian at the seed: two
a step.  Built opt-in (`chord_jacobian`), the iterations' residual from
`i`, `q` and the companion current (every multistep companion's current is
a function of the charges alone), the full Newton from the same seed where
the increments stop contracting.  Measured, full Newton against chord:

| case | time | `G` | answer apart | fallbacks |
|---|---|---|---|---|
| PSP CS gear PSS, 40 points | 23.66 -> 15.54 s | 1017 -> 640 | 1.3e-13 (map 4e-16) | 0 |
| PSP CS gear transient, 3 periods | 5.49 -> 3.35 s | 236 -> 138 | 1.7e-13 | 0 |
| van der Pol gear PSS, 400 | 1.32 -> 1.25 s | 8381 -> 4804 | 2.6e-14 | 0 |
| comparator oscillator, staged gear, 200 | 9.84 -> 9.71 s | 43572 -> 28984 | 5.5e-12 | 5 |
| PWM loop, staged gear, 60 | 1.62 -> 1.86 s | 7409 -> 4941 | 1.6e-10 | 29 |
| diode mixer, gear, fixed grid | -- | 652 -> 397 | < 1e-9 | 9 (Newton iterations 492 -> 1052) |

On the compact model the iteration counts did not move: the predictor's
seed is close enough that the held Jacobian contracts as fast as a fresh
one.  Where the Jacobian is cheap or the steps switch it buys nothing or
costs, hence opt-in.

### `AUTO_JACOBIAN_CODE`: where 'auto' turns the options on (2026-10-01)

Andreas: "Do 1 and 2" (2: `chord_jacobian` and `radau_transform` on where
the Jacobian is expensive, without asking).  The criterion is STRUCTURAL --
the bytecode of the circuit's compiled hdl `G` and `C`
(`compiled_jacobian_size`) -- because a timing probe would make a run's
last bits depend on the machine's load.  Per element it tracks the
evaluation time: `DiodeHdl` 0.2 KB (1.6 us a `G`), HEMT 2.2 KB (8 us),
SPICE diode 20 KB (0.13 ms), MosLevel1 26 KB (0.2 ms), MosLevel3 92 KB (1.2
ms), PSP 1.8 MB (14 ms); hand-written elements and `BSource` count 0.

Calibrated on a driven stage's PSS at 40 points, gear chord / radau
transform against the full Newton (predicted: only PSP and MosLevel3
would gain -- WRONG, every compiled device from 20 KB did):

| device (code) | gear chord | radau transform | answer apart | fallbacks |
|---|---|---|---|---|
| PSP (1.8 MB), switching | -47 % | -43 % | 6e-11 / 4e-8 | 0 |
| MosLevel3 (92 KB) | -27 % (-36 % switching) | -47 % (-42 %) | 0..2e-9 | 0 |
| Gummel-Poon (57 KB) | -33 % | -28 % | 2e-11 / 2e-9 | 0 |
| EKV (26 KB) | -19 % | -29 % | 0 / 2e-16 | 0 |
| MosLevel1 (26 KB) | -17 % (-22 %) | -31 % (-20 %) | 0..2e-9 | 0 |
| SPICE diode (20 KB) | -24 % | -30 % | 2e-12 / 6e-9 | 0 |
| HEMT (2.2 KB) | -5 % | -7 % | -- | 0 |
| `DiodeHdl` (0.2 KB) | -8 % | -2 % | -- | 0 |

and where the options LOST (S15, the chord's own table) every circuit was
of hand-written elements: van der Pol -6 / +18 %, the switching PWM loop
+15 / +79 %.  The threshold, 10 KB, sits in the gap between 2.2 and 20 KB.
The default became 'auto' for both; True / False force them, and only an
ASKED-for transform is warned when PCNR takes precedence.

### `C_KERNEL_SHARE`: 'auto' on the hdl C backend (2026-10-01)

Andreas: "Do 1" (the gap left by 'auto': `compiled_jacobian_size` read a
C-bound function by its Python bytecode).  On the C backend
(`hdl.set_backend('c')`) a model's `G` costs 3.9 us (MosLevel1, from 206)
to 49 us (PSP, from 14 ms).  Re-calibrated there, the same stages
(predicted: PSP little gain, the mid-sized models break-even or a loss --
half right):

| device | gear chord | radau transform |
|---|---|---|
| PSP (small / switching) | -19 / -25 % | -64 / -15 % |
| MosLevel3 (small / switching) | -3 / -6 % | -3 / +8 % |
| MosLevel1 (small / switching) | -5 / -5 % | -5 / +11 % |
| Gummel-Poon | -6 % | +10 % |
| EKV | -6 % | -7 % |
| SPICE diode | -7 % | -6 % |

So a C-bound function counts `C_KERNEL_SHARE` (1/100) of its bytecode:
PSP stays on (1.8 MB -> 18 KB), the mid-sized models go off (under 1 KB) --
forgoing the chord's few per cent to avoid the transform's 8-11 % losses.

### 'auto' re-calibrated: the CSE twins, the fused passes, the const libm (2026-10-02)

Andreas: "Do both F1b and F4" (the fused-evaluation plan's last two
items; F4: re-tune 'auto' now that `G` and `C` are 3-5x cheaper on numpy,
and the C kernels' libm calls are declared const).  The same stages, the
minimum of 3 interleaved runs (PREDICTED: smaller gains, the mid-sized
models near break-even -- the gains did shrink; no model broke):

| device (reference / twin size) | numpy chord | numpy transform | C chord | C transform |
|---|---|---|---|---|
| PSP (1.78 MB / 584 KB), switching | -25 / -32 % | -68 / -25 % | -12 / -18 % | -63 / -9 % |
| MosLevel3 (92 / 43 KB), switching | -14 / -19 % | -27 / -16 % | -2 / -4 % | -3 / +9 % |
| Gummel-Poon (57 / 20 KB) | -10 % | 0 % | -3 % | +15 % |
| EKV (26 / 12 KB) | -5 % | -12 % | -2 % | -5 % |
| MosLevel1 (26 / 13 KB), switching | -6 / -9 % | -15 / +3 % | -3 / -4 % | -3 / +16 % |
| SPICE diode (20 / 8 KB) | -6 % | -10 % | -3 % | -5 % |
| HEMT (2.2 / 1.2 KB) | 0 % | -3 % | | |
| `DiodeHdl` (0.2 KB) | -1 % | +2 % | | |

No crossover moved, so `AUTO_JACOBIAN_CODE` (10 KB) and `C_KERNEL_SHARE`
(1/100) stand.  The proxy stays the REFERENCE function's size
(`_hdl_codelen`), not the twin's that runs: the twin's would place every
model the same with the threshold moved into the 1.2-8 KB gap, but it
moves with `PYCIRCUIT_HDL_CSE`, and the switch is meant to be
bit-identical end to end -- a model near the threshold would change its
Newton with it.  The option descriptions' percentages corrected: chord
"17-47 %" -> "5-32 %", transform "20-68 %" -> "up to 68 %".

Later the same day, after the numpy twins' exact scalar fast paths (round 2,
stage C): PSP and MosLevel3 still gain clearly (chord -12 to -35 %,
transform -11 to -70 %); the mid-sized models' chord -3 to -8 %, their
transform near break-even (-11 to +6 %: Gummel-Poon +5 %, switching
MosLevel1 +6 %) -- no crossover moved.  And the C backend became the
DEFAULT for chained models (stage D): those mid-sized models now run C
where a compiler is, at a hundredth of their size, with chord and transform
off -- so on the test suite every moved result but the `tanh` models' came
from this option flipping, the answers apart by the Newton tolerance (the
rectifier 1.4e-4 at the default reltol, the rest 3e-9 and below).

## `_stamp_plan.py` -- the constant-stamp plan (2026-10-01, the speed plan's P2)

### (module docstring)

Andreas (`/plan`, 2026-10-01): "What can we do to speed up transient,
shooting, pac, pnoise and related analyses ... memory allocation or how the
solutions are shared".  The profiles put the circuit assembly at ~75 % of
a transient or PSS on an ordinary circuit: a Python loop over every element
on every pass, five passes per Newton iterate (C, q, u, i, G) -- and on a
20-section RC ladder 41 of the 42 elements are resistors and capacitors
whose stamps never change, re-stamped every time.  P1 trimmed the loop
(the toolkit probes, the scatter, the default sources); this plan drops it
for the constant elements, bit for bit.

The design rests on two measured facts.  A batched `np.matmul(S (E,k,k),
x[IDX][..., None])` reproduced per-element `np.dot(G_e, x_e)` in 0 of
20000 trials (k = 2, 3; both reach OpenBLAS `gemv_t`), where `einsum` and
hand-written sums did not; hence the per-process self-check per stamp
size, which falls back to per-element `np.dot`.  And a `bincount` bin
starts at +0.0, so it is never -0.0 and adding +-0.0 to it changes no bit:
the constant stamps' exact zeros are left out of the matrix passes, and of
the vector passes for a finite `x` only (`0 * inf` is NaN).

THE CLASS TABLE, by `pair_kind` (method identity; G with i, C with q) on
2026-10-01 -- 'cached' a stamp built in `update()`, 'zero' `Circuit`'s
zero method with its default partner, None re-stamped every pass:

| class | G | C | why |
|---|---|---|---|
| R, G, VS (+ VSin, VPulse, VPWL, VExp, VAM, VSFFM), VCVS, CCVS, VCCS, CCCS, Nullor, Transformer, Gyrator, IProbe | cached | zero | declare `_constant_stamps = ('G',)` |
| L, SVCVS, CoupledInductors | cached | cached | declare `('G', 'C')` |
| C | zero | cached | declares `('C',)` |
| IS, ISin, IPulse, IPWL, IExp, IAM, ISFFM | zero | zero | sources only (`u`) |
| Idt, Idtmod, IdtmodCircular, IdtmodQuadrature (`_IdtBase`) | None | cached | `C` declared; `G` / `i` their own |
| Diode, VCVS_limited, VSwitch, ISwitch, TLine | None | zero | `G` computed per state |
| BSource, NonLinearVCCS | None | None | own `i` (and `q`) -- both report `Circuit.linear` True |
| SubCircuit, ProbeWrapper, CircuitProxy | None | None | assembled as a child (a nested SubCircuit has its own plan) |
| every hdl `Behavioural` | None | None | generated `i` / `q`, not `dot(G, x)` |

`Circuit.linear` was not used: BSource, NonLinearVCCS, VSwitch and TLine
report True, and `volterra.py` still finds its nonlinear elements by it (a
latent defect noted in the plan, out of its scope).

The gate's one failure on the first build (G67): a doctest's `c.G(...)` on
a lone capacitor came out `array([[0, 0], [0, 0]])` -- every entry of that
pass is an exact zero, so the plan's `bincount` was EMPTY, and an empty
`bincount` is int64 even with float weights.  An empty pass now returns
the loop's float zeros (`test_a_pass_with_nothing_to_stamp_is_the_loops_
float_zeros`).

### The C-bound classes' batches (`_hdl_batch.py`, 2026-10-03, speed round 6)

Rounds 2 and 5 made every device evaluation of a chained hdl class a C
kernel; what a mid-sized model's pass was made of after them was the call
around each kernel: on a 20-MosLevel1 chain `cir.G` cost 56 us, of which
the kernels' C work was 6.7 -- twenty wrapped calls at 2.0 us through the
plan's loop at 0.8-1.1 us each.  `_stamp_plan.split` now gathers, when
the plan is built, the non-constant elements of each C-bound chained
class (no DC pins, the generated method by code identity, the circuit's
toolkit, a kernel of the element's size and the slots' shape, two or more
of them) into a `Batch`; `assemble_matrix` / `assemble_vector` call the
rest per element first, as before, then each batch once: ONE C driver
call (`PASS_C`, built and loaded through `_hdl_cbackend.load_kernel`
under its own key) looping over the elements with the kernel's own
function pointer on `x[NM]` and each element's own pack, the temperature
written into the pack's slot as `CKernel.__call__` writes it, the outputs
written into the elements' slots -- so the `bincount` sums the same
values in the same order.  The batch mirrors the pack tuple the element
holds at the pass and never repacks an existing one (`state_restore`
puts an old `__dict__` back with no epoch move).  Checked every pass:
`ENABLED` (`PYCIRCUIT_HDL_BATCH=0`), the binding, the kernel's identity,
the method's code; per element an instance shadow (PCNR's) or a failed
pack hands that element back to the loop and the rest of its class stays
one call.  The legacy fallback receives the batches' outputs only when
it runs (`_run_batches`, `_every_call`).

Measured (parent dc2214da, 5 interleaved rounds, bytes and statistics
the same): 20 MosLevel1 802.3 -> 445.9 us a step (-44.4 %), 20
Gummel-Poon -40.2 %, 20 PSP -16.5 %, the PSP stage and its PSS
unchanged.  The gate: >= 25 %.

The round's second commit: every chained library class's `u`, `u_dc` and
`dudt` return literal zeros (22 of 22), and `cir.u` called all of them --
forty-four calls at 1.3 us on the 20-MosLevel1 chain, 19 % of its step
after the batches -- to add exact +0.0 to bins that start at +0.0.
`_add_element_subvectors` now skips such an element as it skips
`Circuit.u`'s default (`_hdl_batch.zero_source`: an `ast` read of the
compiled function's source, once per class, `info['_u_zero']`; the
method's code identity per call), on the same conditions plus no `dtype`
and not the 'ac' analysis.  Measured: 20 MosLevel1 446.9 -> 412.7 us a
step (-7.6 %), 20 Gummel-Poon -7.1 %, bit for bit.

The round's third commit, the limiter: `SubCircuit.limit` is a sequential
loop (the next element reads the earlier write-back; `x0` often IS `x`;
duplicate rows take the last value), so the batch is not its shape --
`_hdl_climit.limit_walk` walks every run of C-kernel elements in ONE C
call on the live state, in the loop's order with its gathers and
scatters, and the driver returns at the first element Python must answer
(a hand-written limiter, a nested circuit, a shadow, a failed pack, a
declined call), which the caller handles with the loop's own statement
before the walk resumes.  A circuit with no kernel-capable element (PSP's
limiter is its own Python) is answered from the walk's cache before any
other check.  Measured: 20 MosLevel1 409.7 -> 352.1 us a step (-14.1 %),
20 Gummel-Poon -18.5 %, the PSP cases unchanged, bit for bit.  The round
closed: since dc2214da the 20-MosLevel1 chain is 802 -> 352 us a step.

## `_tran_core.py` -- the evaluate core (2026-10-03, speed round 7's commit 1)

### (module docstring)

Andreas: "How could we take a big speed step?" -> the Newton iterate in C,
planned with a two-day ceiling first; the ceiling refused the iterate
(1.15x against a 1.4x gate: the plan's premise was the in-run tree's
inclusive timers, which rank pieces and do not size them) and he chose
to build its evaluate half for real.  One C call evaluates the passes a
site asks for and the companion in the integrator's own operation order;
the LMM step's three evaluation sites take it and fall back to the Python
path wherever it declines, which happens before any state is touched.
The constant vector groups' products go through the address of the
`cblas_dgemv` numpy itself loaded, so the same kernel runs (a C loop,
with or without FMA, differs on most products); the bincounts are
sequential from +0.0 as numpy's; the companion divides where the kernel
divides.  The state the Python path leaves (`_q_cache`, `_C_cache`, the
memo, `get_diff`'s six) is written from the C outputs; `_C_lookup` is
asked first and a hit skips the C pass, as `_C_at_state` does.

The pinned pair "the transient evaluation and the core" fired during the
build, when a lint fix changed the twin after the record: re-recorded in
the same commit.  Measured (parent 82819d61): 20 MosLevel1 357.8 ->
277.3 us a step (-22.5 %), 20 Gummel-Poon -18.9 %, the PSP cases -3..-8
%, bit for bit.

## `_tran_newton_c.py` -- the Newton solve in C (2026-10-04, speed round 8's stage 4)

### (module docstring)

Speed round 8 sized the step by instruction counts, standalone and warm,
before planning (round 7's lesson): with the evaluate core in place, a
20-MosLevel1 Newton iteration was ~750 k instructions of Python and numpy
around ~250 k of C, at about one iteration a step.  A scratch prototype
(stage 0) put the whole solve in one C call -- the core's own C through its
struct pointer, the reduced `J` column-major as numpy's solve copies it,
numpy's own `scipy_dgesv_64_` (and, for the chord, SciPy's own
`scipy_dgetrf_`/`scipy_dgetrs_`: the chord factors with SciPy's OpenBLAS,
whose LUs differ from numpy's), the walk's driver, the convergence test in
`nrsolver`'s order -- and measured it bit for bit on every served case.  Its
first cut gained 10 %: it made its buffers, handles and struct fields per
solve and checked readiness in two loops, and with one iteration a step a
solve's setup IS an iteration's cost.  Persistent buffers and the converged
point's evaluation folded into the same call (the largest single gain)
gave 1.39x / 1.59x on the 20-MosLevel1 / 20-GP marginal step: BUILD.

The build adds what the prototype lacked: every decline counted before
anything is touched, the circuit's own reasons (a core that cannot serve
it, a stateful limiter, a limiting element without a C kernel -- PSP's)
kept with its stamp plan so it declines at its first check (the
prototype's late decline cost the PSP stage +3.7 % a step) -- the plan the
circuit's dict holds, compared by identity: `_plan_for`'s own check cost
~1.5 us a step inside a run, and a plan gone stale is replaced at the
circuit's next evaluation, a decline being today's Newton meanwhile; the core's
readiness through `_Core.probe`, uncounted, so a decline is counted once,
by the path that follows; the walk's readiness at C speed (tuple compares
of the elements' dicts and packs); the evaluate core's switch honoured.  A
bail rolls back the attempt's one trace, the step's source memo, and
`_newton` repeats the arithmetic from the same seed.  The branch screen
reads `C` at the converged point from `_C_cache` and skips its own write
when its lookup served it; where it fires, the confirmation's speculative
solves overwrite the step's state and `jacobian_only` runs again after it.
Measured (parent d371a594): the 20-MosLevel1 marginal step -33.4 %, 20-GP
-38.1 %, the adaptive MOS run -29.9 %; declining circuits within 0.5 %.

Later the same day: the solve declined on an INSTANCE shadow of the
machinery it stands in for, but not on a subclass overriding it or a patch
on its class or module -- either would have been bypassed on every served
step.  The C error test (stage 3) had to learn the same thing first: it
compares its chain's pieces by identity against originals read once every
one is its module's own (`_paths.genuine`: the name, the module, no
`__wrapped__`; read lazily without that, a test that ran first under a
patch had the patch taken for the original).  The Newton solve now checks
the transient's class exactly (`newton_c:class`) and `_newton`,
`_newton_limiter`, `_residual_and_jacobian`, `_get_nrsolver`,
`_get_scaler`, the reference-row helpers `_newton` calls, `nrsolver`'s two
loops and the evaluate core's Python entry (`newton_c:patched`): ~4 k
instructions a served step.  No subclass of Transient exists in the
package or its tests, and no run was served past such a patch: the one
test helper that wraps `StandardNewton.solve_system` on the class (to log
iterations, `test_stage_predictor`) also shadows `cir.i`, which the C
declined -- the gate shows the reason move.

## `_tran_lte_c.py` -- the adaptive error test in C (2026-10-04, speed round 8's stage 3)

### (module docstring)

The default transient is adaptive gear, and every step attempt is judged:
on the 20-MosLevel1 chain after stage 4, 475 k instructions an attempt, 20 %
of the adaptive run -- `bind` 56 k (a `fields()` walk per call), the charge
LTE 233 k (`compute_lte` 77 k, numpy's solve 95 k of which LAPACK is 82 k,
the reference row cut and restored), the normalised error 133 k (the
reference 58 k, `normalised_error` 59 k).  A scratch prototype put
`_charge_lte`, `_normalised` and their maximum into one C call behind one
new method both controllers call (`StepController._max_error`), checked
three ways (the repo's code, the prototype off, on) on nine cases -- every
`relref`, both controllers, the trapezoid, a pulse with rejections and
order drops, a stateful diode: every attempt's verdict, next step, error
and running reference identical -- and counted -10.0 % on the adaptive MOS
run, -6.5 % on an adaptive Gummel-Poon chain, in fresh processes.  (In one
process the counts were bimodal, a persistent +4 % after some runs in any
mode: the interpreter's state, which `--count`'s fresh processes avoid.)

Exactness is by refusal: every input finite; the floating-point flags
cleared on entry and tested after each part whose numpy calls would check
them, every result stored before its test; LAPACK's own flags cleared as
numpy's solve clears them; a tolerance not positive, a non-finite solution
or a negative zero in the running reference (where numpy's SIMD reductions
may pair signed zeros differently) handed back -- the Python chain then
makes numpy's own warnings and errors and the solve's fallback, as before.
The integrator's scalars -- its coefficient (`-h (h + h_last)`,
`-(h**2)`: no `pow` in C) and the divided difference's step sums -- are
computed in the glue in Python floats, numpy's step values converted: the
same IEEE operations numpy's scalars make (0 mismatches in 2 M draws),
without their warnings; one that overflows takes the Python chain, which
warns as before.  (The first build took Python floats only, and the gate
showed the decline in 372 tests: a period's grid gives numpy's.)  The
chain's pieces are read once a call finds every one its module's own
(the name, the module, no `__wrapped__`): read at
the first call as built first, a test that ran first under a patch of
`third_divided_difference` had the patch taken for the original.  The
plan's Python trim: `StepLTEInputs.bind` skips the field names when called
with keywords only (the transient's way), -0.83 %.  Measured (parent
dbd4e1c3): the adaptive MOS run -9.9 %, an adaptive diode ladder -7.7 %,
the PSS cases -1.0..-1.2 %.

## `tests/test_pinned_pairs.py` -- the pinned pairs (2026-10-03)

Andreas, after speed round 6: "How would we keep the python code and c
code synced?"  A Python object that DEFINES a behaviour and the C text,
or the second Python path, that REPRODUCES it are recorded together as a
pair of source digests; a change to either side fails the test naming
the moved side and its twin until the record is re-made -- after the
twin was checked -- in the same commit.  Four pairs now: `CKernel.
__call__` + `pack` / the pass driver and `Batch`; `SubCircuit.limit` +
`CLimitKernel.__call__` / the walk; `apply_limit` + `device_writeback` /
the limiter prelude and renderer; the assembly loops / the plan, `split`
and `zero_source`.  It also pins that the metaclass defines i/G/q/C/limit/
u/dudt exactly once, which the code-identity marker rests on.  A digest
proves only that something changed; agreement stays the sweeps', the
recorded gate's and the private comparison's job.  Every new C twin --
the Newton iterate above all -- gets a pair in its first commit.

## `_tran_companion.py` -- the step's source memo (2026-10-02, speed round 4's stage B)

### `_source_at`, `Transient.solve_timestep`, `_stage_source`

Andreas, after speed round 3: "Plan for the per-step solver machinery".
The round's measurement (`benchmarks/step_machinery.py --tree`, the log
of 2026-10-02) put the sources at 2.03 assemblies per gear step on the
PSP stage and 1.89 on a 20-PSP chain -- the same `t` each time, since
the Newton's iterations, the chord's residual-only ones and the branch
confirmation's re-solve all ask `u` at the step's end -- and the stage
methods' `_stage_source` closure asked `cir.u` at every stage of every
coupled iteration (Radau 10.24 a step).

Now `solve_timestep` opens a memo (`_u_memo = {}`) and closes it in a
`finally`; `_source_at` serves a `(t, analysis)` it has seen and
assembles the others; `_stage_source` goes through it.  The same
function on the same inputs is the same bits, and nothing a source
reads moves inside a step: no numeric path writes `epar.t` (grep: reads
only), `analysis_kind` is scoped around a whole solve, an element's
state moves only at `accept_step` / `reset_state`.  Outside a step there
is no memo (`_begin_run`, `_memo_clear` and the `finally` all set None),
so a later caller at a time already seen -- the shooting re-entering a
period, a source whose state an accept moved -- assembles afresh.
`provided_function(t)` is still called at every request and added into a
new vector: a caller may count it (the F4 tests do).  `U_MEMO` switches
the memo off, for the byte-identity tests and as an escape.

Measured against the parent (4af78ddd), each tree in its own
interpreter, five interleaved rounds: the PSP stage -4.0 % a step, its
PSS -4.2 %, 20 PSP -3.0 %, 20 MosLevel1 -1.0 % (it seldom iterates
twice); bytes and statistics identical.  The round's other finding, for
the record: the assembly's remaining per-call cost is numpy's own floor
(a minimal rewrite saves 1-5 us a call, 1.5-2 % of a step), so no
assembly stage was built.

## `_tran_newton.py`, `analysis.insert_row`, `nrsolver.py`, `circuit._hook_elements` -- the Newton's plumbing (2026-10-02, speed round 4's stage C)

### `insert_row`, `_newton_tolerances`, `_reduced_row_names`, `_as_float`, `_ring_push`, `_hook_elements`

The round's in-run tree put the plumbing around a step's Newton at a
third of the PSP stage's step: the reference row in and out through
`np.insert` (4.8 us a call, two an iteration in the limiter alone) and
one-element `concatenate`s, the tolerance vectors built in both flavours
and reduced at every step, the row names for a failure message a step
almost never writes, `toolkit.array` copies of sums nobody else held,
every magnitude taken twice in the Newton's test, the chord's held
Jacobian's magnitudes every iteration, a generator context manager per
evaluation session, both element hooks polled on every element, the
history rings rebuilt through a one-row array and a view.

Each item is now the same arithmetic on the same values or a cache of a
pure function keyed on everything it reads.  `analysis.insert_row` copies
a float64 vector into a fresh one through three slices (the toolkit form
for anything else; the inserted row +0.0 both ways).  `_newton_tolerances`
builds both flavours and both widths once per `(sizes, iabstol, vabstol,
irefnode)` and keeps them READ-ONLY -- no consumer writes into one, and
one that did would raise instead of poisoning every later step; a
symbolic tolerance is built every time.  `_reduced_row_names` is kept per
circuit shape.  `_as_float` returns a fresh float64 sum as it is.
`_companion_at(x, C)` takes the `C` its caller looked up.  `_memo_get`
builds no key on an empty memo.  `abs(F)`, `abs(x_next)` once an
iteration; the chord's `abs(J)` once.  `_evalhint.evaluating` is a
`__slots__` class.  `SubCircuit._hook_elements(name)` lists the elements
whose hook is not `Circuit`'s no-op or is shadowed on the instance, kept
per topology under the stamp plan's contract (an instance patch after the
list was built is seen at the next topology change); `next_event` keeps
`maximum(t, min)` with `inf` where none declares one.  `_ring_push`
copies the ring into a fresh array.

Measured against the parent (bfaf3ecd), five interleaved rounds, bytes
and statistics identical: the PSP stage -13.3 % a step, its PSS -11.4 %,
20 MosLevel1 -5.9 %, 20 PSP -2.8 % -- above the predicted 5-7 %: the
Python removed ran cold between the kernel's calls, so the in-run
inflation the plan noted cut both ways.

## `_tran_predictor.py` -- the multistep fast path, and a tie (2026-10-04, speed round 8's 2b)

### `_predict_fast`, `_pred_fast_nodes`

A prediction with no extra nodes from a history wholly behind the target
-- every multistep step -- in 90 k instructions where the general path
took 168 k: behind the target, nearest-first IS newest-first, so the order
and the 1e-13 dedupe depend on the history alone and are kept per history
list; the fit's times in Python floats; the clamp the same ufunc.

Except where two distances to the target ROUND EQUAL.  The round's deep
hypothesis run (3000 draws a twin) found it: times 0 and 1.6e-45 behind a
target 1.4e-3 away are one distance in floating point; the general path's
stable sort lists the older of the two first, its dedupe keeps that one,
and its fit took another node than the fast path's -- -2.73 against 0 in
one row.  A tie needs two times closer than one ulp of the larger
distance; the transient's own history keeps its times 1e-14 apart and
never holds one, but the fast path's contract is the general path's bytes
for any history behind the target.  It keeps the smallest gap between two
history times with its per-list cache and steps aside (`pred:tie`) where
that gap is within two ulps of the largest distance -- every tie, by
construction.

## `_tran_predictor.py` -- the predictor's bookkeeping, its weights kept (2026-10-02, speed round 4's stage D)

### `_predict_state`, `_fit`

The round's in-run tree put the predictor at 7 % of the PSP stage's
step: the node times deduplicated through a generator per pair, the
newest nodes sorted twice, and the Vandermonde system solved at every
step for weights that depend on the normalised node times alone -- which
the shooting repeats at every period walk and a uniform grid repeats
too.  The dedupe is a plain loop with the same rule and first-wins order,
the newest are sorted once, and `_fit` (a method) keeps its weights on
the instance keyed by `tau`'s bits (`PRED_WEIGHT_MEMO`; bounded,
forgotten at `_memo_clear` so a `solve` starts afresh while a shooting
solve's walks share it).  `np.vander` stays: its `multiply.accumulate`
powers are the bits.  The clamp and the periodic rows are untouched.
Measured against the parent (6f87af70): the PSP stage -5.2 % a step,
its PSS -3.8 %; bytes and statistics identical; the old code
transliterated is the test's oracle over 300 random node histories.

## `_tran_radau_c.py` -- Radau's dense stage Newton in C (2026-10-05, speed round 9's B2)

### (module docstring)

Radau IIA(3) is PSS's default method, and with stage 8 (B1) its stages
took the evaluate core's passes; what remained of a 20-MosLevel1 step
(5.3 M instructions) was the Newton's Python: `_coupled_stage_system`
(1.40 M a step, two calls: nine blocks each through `remove_row_col`), the
passes' wrappers (6 calls), the limiter's wrapper per stage (0.42 M),
numpy's solve wrapper, the convergence test, the closure's glue.  B0, the
same system vectorised in guarded Python, was refused at -6.95 % (its line
-8 %); here the whole undamped, unshunted Newton is one C call -- the
assembly in `_coupled_stage_system`'s operations and order over the full
width (Python's `sum` starts at the int 0: `0 + A_i0 K_0` makes a -0.0 a
+0.0, so the C adds `0.0 +` too), numpy's own `scipy_dgesv_64_` on the
reduced system written column-major as numpy's solve copies it, the
reference row's `Y_prev + 0.0`, the walk against the previous stage, the
convergence test after the next assembly, as the Python tests.

What it does not mirror, it bounds or hands back: `Rnorm = sum |R|` (read
by the undamped Newton only for its overflow warning) by a bound, every
`|R| < 2^1000`, not numpy's pairwise sum; a floating-point exception in its
own arithmetic (overflow, invalid, division, underflow -- numpy would warn
or might) by handing back, after clearing what the core, the walk and
LAPACK raised (numpy clears before each operation); a non-finite value, a
LAPACK info, the walk stopping, maxiter, and more assemblies than its
buffers hold (`MAXA`, 8) likewise.  A hand-back rolls back the source
memo's entries and counts the C's calls made -- the C Newton's keys plus
the source plan's (`_tran_newton_c`'s `u_keys`; stage 5's `src.u:*` were
missing from that list until its switches-off check found a bail's trace)
-- and the Python Newton runs from the same seed.  Where it converges,
every assembly's three stages go into the device memo through
`_memo_put` in the Python's order, and the source memo's hits for every
later assembly are counted: the transform path, `_finish_stage_step`, the
next step's `q_n` and PSS read the same records.

The tests drive one captured step's `_stage_newton` (the closure kept by
a spy on `_coupled_stage_solver`, the memos copied as its Newton found
them) with drawn seeds, step sizes and entering charges, the C on and
off: the stages' bytes or the exception, both memos, the counts and the
warnings the same, and on a hand-back every count but the C's own.
Measured against the parent (ccd7bb33): the radau 20-MosLevel1 marginal
step 5.46 M -> 3.20 M instructions (-41.4 %); `--check`, the same step -43.1 %
in wall time.  The round's lesson from stage 6 (refused: fewer
instructions, slower in situ) does not touch it: one call stands in for
milliseconds of Python and numpy a step, not microseconds.

## `_tran_radau.py` -- the transform's factors once a step (2026-10-05, speed round 10's B3.1)

### `_RadauStages._radau_frozen`, `_FrozenTransform`, `linearsolver.NumpyLU`, `ComplexKLUSolver.prepare` / `solve_prepared`

The cost transform (`_rk_step_transformed`; PSS's default on PSP-class
circuits under 'auto') is a SIMPLIFIED Newton: `C` and `G` frozen at
`x_n` for the step.  Its two factors, `(gamma_r/h) C + G` and
`((alpha + i beta)/h) C + G`, are therefore the same at every iteration,
and so are their factorisations -- yet `_radau_transform_solve` rebuilt
`P`, both factors, numpy's LU of the real one, SciPy's CSC copy of the
complex one and its KLU refactor at every iteration (870 k instructions
an iteration on the PSP stage, 2.03 iterations a step on its PSS).
`_radau_frozen` makes them once a step and `_FrozenTransform.solve` is the
per-iteration solve's bytes; `_radau_transform_solve` stays, as the path
of a declined step.

Three traps, each found by checking the piece against what it stands in
for before it was used:

- NUMPY'S SOLVE IS `dgesv`, AND `dgesv` IS NOT `dgetrf` + `dgetrs` ON A
  THREADED OPENBLAS.  With one right-hand side `dgesv` factors on one
  thread below 10000 unknowns; a standalone `dgetrf` of 100 unknowns or
  more uses several, and rounds differently (210 of 360 solves differed
  from 100 to 1000 unknowns; below 100, 9000 of 9000 agreed -- the first
  check drew sizes below 50 and would have shipped the defect).
  `NumpyLU` makes numpy's own `dgesv` call at its first solve and keeps
  the factors for `dgetrs`.
- A REFACTOR IS NOT A FACTOR.  KLU's refactor of the same values with the
  same pivots repeats its own bits, so `solve_prepared` skips the refactor
  `solve` would repeat for the same record (`_fresh`) -- but after a full
  factor it still refactors once, as `solve` did: skipping that one moved
  the answers of 15 tests.
- NUMPY SHOWS A WARNING ONCE A LINE.  The right-hand sides and the update
  are helpers both solves call (`_transform_rhs`, `_transform_back`); a
  frozen solve forming them at its own line would show, after a declined
  step, a warning the per-iteration line had already shown.

SciPy's `csc_matrix(A)` of the dense complex factor cost 382 k
instructions; `_csc_of_dense` builds the same arrays with numpy (copies:
the nonzeros column by column, rows ascending, a NaN counted, a signed
zero not) and the residual check's product is SciPy's own `csc_matvec`,
as `_matmul_vector` calls it.  Without it the stage read -5.8 % on the
radau PSS (its line -7 %); with it -10.3 %.  Measured against the parent
(3fb7ec26): the radau PSS of the PSP stage 1.303 G -> 1.169 G instructions
(-10.3 %), -13.2 % in wall time; nothing else moved.

## `_tran_radau_tc.py` -- the transform's Newton in C (2026-10-05, speed round 10's B3.3)

### (module docstring), `_RadauStages._transform_loop`

With its factors made once a step (B3.1), the transform's Newton loop was
half of the PSP stage's radau PSS step: 3.5 M instructions of Python and
numpy around ~0.4 M of device kernels.  The loop became a method of its
own (`_transform_loop`, its statements unchanged) so a C call could stand
in for it and be tested against it, and `_tran_radau_tc.solve` runs it in
C: the stages' passes through the evaluate core, B2's residual, the
transform's right-hand sides and update, numpy's kept LU and the complex
factor's KLU on their own handles, the walk, the convergence test.

Two things it had to READ rather than assume, each found before the C was
written:

- NUMPY'S COMPLEX PRODUCT FUSES ON THIS BOX.  numpy 2.5.2's X86_V3 loops
  (AVX2, FMA3) compute `re = fma(ar, br, -(ai*bi))` and `im = fma(ar, bi,
  ai*br)`, at every position of every length (20 000 of 20 000; the
  separate products differed in 958).  A build or CPU without FMA would
  multiply separately, so `cmul_mode` asks numpy itself, once, on inputs
  where the two forms round apart -- and on signed zeros and underflow,
  where a real vector cast to complex (the right-hand sides' `P F`) tells
  them apart only by the sign of a zero -- and the C takes the form numpy
  showed (the backend compiles with `-ffp-contract=off`: `fma` is called,
  never inferred).  Neither form: the path is off.
- NUMPY'S COMPLEX `abs` IS ITS OWN.  It differs from `hypot` in ~10 % of
  values, so the residual check `solve_prepared` makes (fall back to a
  fresh factor past 1e-8) cannot be reproduced bit for bit; it can be
  DECIDED outside the band the two roundings can reach (`32 (m + 4) eps`
  of `max_k sum_j |A_kj||x_j|` over `max|b|`, plus `16 eps` of the
  residual), and anything inside the band or past the tolerance is handed
  back.  The residual sits ~1e-15 below a 1e-8 tolerance in practice.

The hand-back protocol leans on B3.1's two facts: a refactor of the same
values and pivots repeats its bits, and `dgetrs` on `dgesv`'s factors is
`dgesv`'s solve -- so the factorisations' state the C leaves (the kept LU
after its first `dgesv`, the numeric holding the record's refactor) is
the state the Python loop would have left, and a Python loop run after a
hand-back makes the C's iterates again.  A test planted a zeroed record
expecting KLU's refactor to fail; on this block structure the refactor
SUCCEEDS (a 3x3 probe's failed), the solve is not finite and the residual
check hands back -- the expectation was wrong, not the C, and the refactor
and solve failures are planted in the C's function pointers instead.
Measured against the parent
(1e4b6c67): the radau PSS of the PSP stage 1.168 G -> 0.805 G instructions
(-31.1 %), -26.5 % in wall time; nothing else moved.

## `_tran_radau.py`, `_tran_branch.py` -- a Radau step's stage passes fused; the screen's proxy (2026-10-05, speed round 10's B3.2)

### `_RadauStages._stage_end_passes`, `_PeriodWalks._walk_stage`'s hint, `Transient._branch_screen`

A Radau step's stages are read by three readers after its Newton: the
branch screen (`C` at every stage), the step end (`_finish_stage_step`:
`q`, `i`, `C`, `G` at `x_{n+1}`) and, in a PSS walk that factors the step
at once, the shooting (`_C_at`, `_G_at` at every stage).  Each evaluated
on its own, through the circuit's passes: the device kernels ran twice
at the same state.  The memo (`_memo_get`/`_memo_put`, rolled per step)
was already the place where one reader leaves an evaluation for the
next; `_stage_end_passes` fills it first, in one core call a stage --
`qiCG` at `x_{n+1}`, `CG` at the others where the walk says it will read
`G` there (`_stage_G_read`, set around `solve_timestep` alone, so a step
that raises leaves it down) -- and every reader then finds its own.
Bit for bit: the core's passes are the circuit's own, and the memo
records only where a recorded value is the one re-evaluating gives
(`_memo_ok`).  Not on the dense path, whose Newton records all four
passes at every assembly's stages already.

What made the gain FUSION rather than the core: for PSP the core's single
passes cost what the circuit's own do -- the kernel is the cost -- while
the kernels' shared statements make `CG` 301 k against 476 k apart.

The screen's proxy took `np.max(np.abs(C))` three times, the diagonal's
`abs` separately and the positive mask twice; each is now made once, the
same reductions of the same arrays.  Its oracle is the old code, verbatim,
in the test, on every branch enumerated: a drawn test missed a planted
`>=` in the positive mask (the branch wants a zero beside large diagonal
entries and no NaN, rare in draws).  Measured against the
parent (584411d3): the radau PSS of the PSP stage -14.9 % instructions,
-18.7 % in wall time; gear, through the proxy, -4..-6 %.

## `shooting/_pss_inner.py`, `_pss_walks.py` -- a stage walk's per-step trims (2026-10-05, speed round 10's B3.4)

### `_InnerTransient._sync_limit_at`, `_stage_block`, `_skip_step_mats`, `_insert_refnode`

A stage walk factors every step it takes: `_stage_step` reads `C` and `G`
at the step's stages (from the device memo since B3.2) and assembles the
coupled `3m x 3m` system.  Four of its costs were nothing's: the limit
sync before each `G` read (`limit(x, x)` -- it moves only a limiter that
KEEPS state, `Diode`'s, and the PSP circuits keep none), the block loop's
nine slice assignments, the inner transient's reduced `Jf`/`Geq`/`C` after
every step (only `_walk_lmm` reads them) and the generic `concatenate`
inserting the reference node (`analysis.insert_row` had replaced it in the
Newton in speed round 4).

The one pass's lesson is the NaN payload: IEEE addition commutes but for
which NaN survives a NaN + NaN, and numpy decides that by its loop -- the
SIMD body or the scalar tail -- so a rewrite that changes an array's
shape can keep the other NaN even in the same operand order.  A drawn
test with only `np.nan` cannot see it; NaNs of three payloads did.  The
one pass declines non-finite blocks and step; the loop keeps them as it
always did.  Measured against the
parent (4c25985b): the radau PSS -4.15 % instructions, -1.2..-2.0 % in
wall time.

## `_tran_radau.py`, `_tran_radau_tc.py`, `linearsolver.py` -- what a radau step keeps for the next (2026-10-05, speed round 10's B3.5)

### `_radau_frozen`, `NumpyLU.make(reuse=)`, `_csc_of_dense(last=)`, the C's inputs

B3.1 made the transform's two factorisations once a step; each step still
made them from nothing: a new LU (buffers, LAPACK arguments), a new CSC
pattern for the complex factor, `P = diag(lam) Tinv / h` and `V`'s entries,
and the C call (B3.3) re-marshalled every pointer and entry into its
struct.  Now a step keeps what the last one made where the values say it
is the same: the LU refilled in place once the frozen transform that read
it is gone (`NumpyLU.owner`, a weak reference -- never while it lives, so
a nested step cannot overwrite an outer step's factors), the pattern's
`Ap`, `Ai` and key where the nonzeros fall where they fell, `P` and `V`
per step size and type, and each C input set where it changes.  The tests
that pin the keeps plant each one stale and require the C to differ from
the Python loop: a keep whose bookkeeping is right can still serve wrong
values, and only the answer shows it.  Measured against the parent
(be80acf2): the radau PSS -6.0 % instructions, -9.2 % in wall time.

## `_stamp_plan.py`, `func.py` -- the sources evaluated directly; `Sin.f` at a scalar time (2026-10-06, speed round 10's B3.7)

### `source_direct`, `_SourceDirect`, `_added`, `_genuine`; `_sin_scalar`

A `u` pass called each source's own `u` -- an `iparv.v` read through
`ParameterDict.__getattr__`, the time function's `f(t)`, a 3-entry array --
and the loop gathered and bincounted them: ~117 k instructions on a
four-element circuit, ~60 k of it around two additions.  Where every
element a pass calls is a VS- or IS-family source, the pass is now those
additions and the loop's bincount in Python floats -- bins from +0.0,
contributions in element order -- with the add made as `u` makes it
wherever numpy could warn about it (its warnings re-emitted from `u`'s line
and module), and its checks stamped: the sources' dicts are watched, not
the other elements' (a limiter writes its element's dict every
iteration), and what no watcher sees is compared on every call against the
definitions captured at first use.  `Sin.f` at a scalar time is its
expression in Python floats around the toolkit's own `exp` and `sin`,
served only where nothing the expression does could warn.  Measured
against the parent (0a05e7ae): the radau PSS -4.4 % instructions and -4.0 % in wall time, the gear
device chains -3.4 % and -6.3..-6.9 %.

## `_watch.py`, `_tran_newton_c.py`, `_tran_core.py`, `_stamp_plan.py` -- the watch counter between solves (2026-10-06)

### `_ParRec.bp_obj`, `_watch.counted` and `HOLD`, `_walk_arm`

One counter (`_watch.EPOCH`) stands for every watched dict, and every
stamped check compares its stamp with it.  Every transient solve moved it
twice: its operating point is made inside `analysis.analysis_kind`, which
writes `analysis_kind` into `epar`'s values and restores it, and the C
Newton's parameter record watched those values.  So every stamp broke
twice a solve, and with re-arms capped at 16 a transient solved again and
again -- a sweep, a Monte Carlo loop -- had every stamped check unstamped
from about its 8th solve on: +13..+20 % a solve.  The record now watches
`epar`'s own dict and its parameter table and compares the one value it
keeps from them, `bypasstol`, by its object; a re-arm after a stamp that
held 8 checks or more no longer counts toward the cap; and the limiter
walk, stamped only by its full setup before, is stamped again after its
tuple check.  Measured: one transient solved again, -22.6 % instructions
a solve; a single solve unchanged.

## `integrator.py` -- a stage method's tableau classified once (2026-10-06, speed round 11)

### `RungeKuttaIntegrator.is_stiffly_accurate`, `stage_structure`, `_classify`

The stage step asks, on every step, whether its tableau is stiffly
accurate and what structure it has (`_solve_timestep_rk`); each answer
was `np.allclose` over the tableau's rows, ~130 k instructions apiece --
0.29 M of a Radau step and 0.43 M of a TR-BDF2 one, for the same answer
every time.  Each is now kept with the tableau `butcher()` returned for it,
and given again while `butcher()` returns that object.  Measured: the
radau PSS of the PSP stage -7.8 % instructions, a radau transient of a
20-MosLevel1 chain -8.3 %.

## `analysis.py` -- small reduces and inserts as one indexed copy (2026-10-06, speed round 11)

### `_reduce_small`, `_insert_small`, `_TAKE_2D`, `_TAKE_1D`

`_reduce_ndarray` drops the reference row and column with four slice
copies, the layout fastest for the large matrices it was written for; on
a circuit of a few nodes each copy is its own numpy call -- ~28 k
instructions for a 7 x 7 -- and `insert_row`'s slices ~12 k.  Up to a
32 x 32 matrix or a 256-entry vector both are now one indexed copy
through a cached index, the same entries bit for bit in a fresh array:
~7 k and ~3 k.  Past those sizes the slices stay.  Measured: the radau
PSS of the PSP stage -4.9 % instructions and -1.8 % in time, the gear PSS
-5.5 % and -1.0 %, the van der Pol PSS -4.7 % and -2.4 %.

## `shooting/_pss_walks.py`, `_tran_radau_tc.py` -- the stage step reads each stage once; the transform solve's checks (2026-10-06, speed round 11)

### `_stage_reads`, `_readers`; `solve`'s checks

The shooting's coupled stage step reduced each of the step's stages and
then read `C` and `G` there through `_C_at` and `_G_at`, each of which
inserted the reference row back, looked the transient's device memo up
and reduced the matrix again -- seven readers a sensitivity step, for
matrices the step had just evaluated at those very states.  Where the
readers would read exactly that (they are their modules' own, nothing
keeps limiting state, no junction goes through PCNR, each full stage's
reference entry is +0.0, the memo holds both), each stage is now looked
up once at the full state and its two matrices reduced as the readers
reduce them, the limit sync's skip counted as theirs; elsewhere the
readers run.  The transform solve's array checks compare the float64
dtype by identity first, and its snapshot of the source counters is a
tuple.  Measured: the radau PSS of the PSP stage -2.3 % instructions
and -1.3 % in time.


## `_stamp_plan.py`, `circuit.py`, `_tran_core.py`, `toolkit.py` -- the per-call imports once at module level (2026-10-06, speed round 12)

### `_plan_for`, `_ineligible`, `_add_element_submatrices`, `_add_element_subvectors`, `evaluate`, `matrix_from_entries`, `SubCircuit.limit`

A function-level import costs ~1.6 k instructions every time the function
runs (`import numpy as np` ~0.8 k), and seven functions on the assembly's
and the step's paths ran one or two on every call: the stamp plan's two
checks (every assembly), the two assembly passes (numpy again: `np` is the
module's), the evaluate core (the toolkit class; the BDF2 kernel under
gear), the numeric toolkit's matrix builder, the limiting loop.  Each is
now imported once at module level -- as the module, read by attribute
where a name was taken from it, so a monkeypatch of that attribute is
still seen -- except `_hdl_climit`, which imports `hdl` and through it
`circuit`.  `_ineligible` compares the dtype by identity first.  Measured:
the van der Pol PSS -2.8 % instructions and -1.7 % in time, the diode
ladders -2.9..-3.1 % and -1.5..-1.9 %.


## `_numeric.py`, `_tran_companion.py` -- numpy's module reductions at their floor (2026-10-07, speed round 12)

### `alltrue`, `_C_lookup`, `_q_at`, `_memo_put`, `_memo_get`

`numpy.all` costs 11.6 k instructions a call where the array's own `.all()`
costs 4.5 k -- the same `logical_and.reduce(a, None, bool)`, without the
module function's dispatch -- and the numeric toolkit's `alltrue` was
`numpy.all` itself: 16.6 k calls a van der Pol PSS (the C lookup's
comparison and the Python Newton's convergence tests).  It now takes an
exact ndarray's own method where nothing else is asked, numpy's function
otherwise.  The C lookup missed 83 % of its lookups after a ~20 k
comparison: a miss is decided at the first entry where two exact float64
vectors differ (NaN never equal, signed zeros equal, as `==` has them),
there and in the charge's identical comparison (`_q_at`); the full
comparison otherwise.  The device memo's key skips `asarray` for an exact
float64 array (the same bytes).  Measured: the van der Pol PSS -3.95 %
instructions and -3.8 % in time, the gear ladder -1.06 %
instructions.


## `_tran_radau_tc.py`, `_tran_radau.py`, `_tran_core.py` -- the step's end passes in the transform's C call (2026-10-07, speed round 12)

### `RADAU_TC`, `solve`, `_RadauStages._stage_end_passes`, `passes_held`, `passes_count`

After the transform Newton's C call returned the converged stages, the
step made its end passes through `passes` -- 2.7 core calls a step, each
~57 k instructions of readiness and marshalling around its kernels.  The C
now makes those calls itself once converged, in `passes`'s form and order,
and hands them over with the stage list it returns; `_stage_end_passes`
takes them for those stages where every call it replaces would find its
readiness as stamped and nothing watched has changed since, making the
counters those calls make, and makes the calls otherwise.  Measured: the
PSP stage's radau PSS -3.88 % instructions and -4.1 % in time.


## `_tran_radau.py`, `_tran_radau_tc.py`, `linearsolver.py` -- the frozen factors in C (2026-10-07, speed round 12)

### `_RadauStages._radau_frozen`, `fold`, `FOLD_C`, `ComplexKLUSolver.prepare_values`

The transform's frozen factors cost 131 k instructions a step on the PSP
stage's radau PSS -- numpy calls on 6x6 arrays under an errstate, the
complex factor's mask, comparison and gather for its KLU record, the
kept LU's finiteness test.  One small C call now makes both factors as
numpy does (its operations in its order, its complex product's form, the
flags its errstate raises on) and the complex one's packed values where
they fall in the last record's pattern; `_radau_frozen` makes every
decision as before, in its order, and `prepare_values` is `prepare`'s own
record of those values.  Its own call, not inside the transform's: a
factor that would raise must decline before the transform's readiness
checks run, as it did.  The C's work is nothing at this size; the
marshalling is everything -- buffers kept per size.  Measured: the PSP
stage's radau PSS -1.55 % instructions and -3.8 % in time.
