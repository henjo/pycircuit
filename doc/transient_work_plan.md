# Transient subsystem — the work plan

> ## RESUME HERE — state as of 2026-07-31
>
> Last commit changing **simulator behaviour**: `1122c31` (stage 4i).
> Anything after it is a benchmark or documentation unless it says otherwise — so a newer
> HEAD does not by itself mean the code has moved underneath this block; check
> `git log --oneline 1122c31..HEAD -- pycircuit/`. Branch `cna-jax-vectorization`,
> **pushed to `origin`** (`git@github.com:henjo/pycircuit.git`).
>
> **Suite: 778 passed, 6 skipped, 0 failed, 722.88 s** (`-m "" --timeout=400`). Nominal
> ~8-13 min, but one run of the identical tree took **31m41s** purely from machine load —
> see trap 2. **Doc build: succeeded, 2 warnings, 0 ERROR.** Working tree clean.
>
> Transient regression tests live in `pycircuit/circuit/tests/test_transient_repairs.py`
> (renamed from `test_transient_stage1.py`); its docstring maps sections to plan stages.
>
> ### The next action, concretely
>
> **`relref='sigglobal'` — decided, reverted by its own gate, and the blocker is now
> gone.** Gate D3-a failed because the trapezoidal estimator was mode-contaminated and
> accuracy stopped falling monotonically as `reltol` tightened. Both second-order
> estimators are now flat to within 1.3% across the whole reachable step-ratio range
> (4g(b), 4i), so the premise of that failure no longer holds. This is a one-line default
> change with its gate already written — re-run D3-a and D3-b as they stand.
>
> ### After that, in order
>
> 1. **`lte_vabstol` back to 1e-12**, and rewrite gate 1-5, which asserts a property only
>    true under `pointlocal`. Follows directly from `sigglobal` landing.
> 2. **Decide what `lte_formula` is for.** After 4i it selects nothing for either
>    second-order method except on the single fallback step of a run. That makes gate 4f
>    ("record rejections and force-accepts for all four integrator/formula combinations,
>    so the choice is evidenced") a question about whether the knob should exist at all
>    rather than about which default to pick — there are no longer four distinct
>    combinations to compare. **Maintainer's call**, and 0.3b's deferred decision is
>    subsumed by it.
> 3. **4a, 4a-bis, 4h**, 0.3d's `chgtol` guard, and the rest of
>    `doc/src/circuit/lte_dae.rst`'s variable-step story.
> 4. **The `iq` seed.** `transient.py` starts every run with `_iqlast = zeros` against a
>    true `q'(t_0)` that is generally nonzero. Measured under 4g(b) as an O(h^3) effect,
>    one order below the local truncation error, so it is not blocking — but it is wrong
>    in principle and cheap to fix when something else touches that code.
>
> **Note for 4d.** `Trapezoidal('ywr')` is by far the heaviest user of the rejection cap
> (25 force-accepts on the stiff RLC at reltol 1e-4, against 1 for either `Gear2`), so
> deleting that branch will move the force-accept counts a lot. 4b's end-to-end regression
> test deliberately does **not** depend on it — it is pointed at the shipped default at
> reltol 1e-3 — but re-run `benchmarks/transient_stage4.py --forceaccept` after 4d rather
> than assuming the counts carry over.
>
> ### Done
>
> Stage 0 (all reviews, measurements, decisions 0.3a-d) · Stage 1 (silent failures) ·
> Stage 2 (2.42x bit-identical, 5.19x with single-threaded BLAS) · three post-stage-2
> improvements (2+.1 `__getattr__` memo, 2+.2 skip the unread residual, 2+.3 `relref`) ·
> Stage 3 (`firststep`; `reltol` controls accuracy, 90.8x) · Stage 4 parts 1-4: 4c · 4g(a) ·
> **4e + 4b** (61 accepted step ratios outside BDF-2's zero-stability bound across the
> measured grid become 0, worst ratio 10.000 -> 2.000, force-accepts now warn — and the
> warning immediately corrected a conclusion this work had already written down) ·
> **4g(b)** (the trapezoidal estimator differences a mode-free quantity: asymptotically
> exact where it held a 33% bias, 1540x -> 1.26x step-ratio spread, 757 -> 23 rejections on
> the stiff case — and its first gate refuted the "1/h law" the work was premised on) ·
> **4i** (one `q'''` estimator for both second-order methods, taken from a third divided
> difference of the charge: Gear2's step-ratio spread **119x -> 1.0041x**, trapezoidal's
> 1.26x -> 1.0084x, and the shipped default gains 1.22x accuracy — but gate 4i-4's own
> "clause that matters" failed and is recorded as a failure).
> Every gate outcome is recorded in place below; the completion records are appended at the
> end of this file.
>
> ### Open decisions
>
> * **D1 said keep-and-fix `_solve_coupled`.** Two ways to honour that, and **one must be
>   chosen before either is started**: the cheap one (drop the Fang citation, rename the
>   flag to the rejection loop it is, fix the four ignored inputs) or **stage 12**
>   (implement the paper). Stage 12's entry condition is a *measurement* — gate 4f's
>   rejection counts — not a date.
> * **`relref='sigglobal'` as default**: decided yes, reverted by its own gate. **The
>   blocker is gone as of 4g(b)** — gate D3-a failed because trapezoidal's estimator was
>   mode-contaminated and it no longer is, so this is now a one-line change with its gate
>   already written.
>
> ### Traps that have already cost time — check these before believing a result
>
> 1. **A transient whose `timestep` is near the circuit's own time constant runs at
>    `max_step` end to end**, so the controller decides nothing and **no tolerance can be
>    observed to act.** Three draft gates "passed" that way while measuring nothing. Use
>    `uic=True`, a `timestep` well above the LTE-chosen step, and check step count >>
>    `tend/timestep` first.
> 2. **A single-sample suite runtime on this box is worthless.** One read as a 35%
>    regression that controlled measurement showed to be 3%; individual tests vary ±57%.
>    Use min-of-N interleaved, and check a non-clock invariant (step count, assembly count,
>    cProfile call count).
> 3. **An estimator harness must feed the previous *companion currents*, not exact
>    derivatives** — the estimator never receives the latter. Doing so made Euler look flat
>    at 0.501 at every step ratio, pure artefact. **Confirmed again 2026-07-31, and it
>    explains a second number**: fed exact derivatives, the trapezoidal estimator reads a
>    flat **5/6 = 0.8333** — which is exactly the "5/6" `integrator.py` claimed for it and
>    decision 0.3b called measurably wrong. The claim was not invented; it was measured with
>    this harness's own first defect. **A third variant, 2026-07-31:** the generated
>    comparison table in `lte_dae.rst` fed **Gear-2 a trapezoidal companion history** --
>    the two recursions are different functions -- and reported 41.4 where the truth is
>    83.06. Each method must be fed its own. A generated number is only as good as the
>    history it is built from, and all three of these read as plausible.
> 4. **A gate that can pass against an empty result is not a gate.** Gate 2b's first version
>    compared `u` at t=0, where a sine is zero, and reported "exactly equal" against an
>    all-zero vector.
> 5. **`drift == 0.00e+00` is unsatisfiable for anything that changes which BLAS kernel
>    runs.** Declare "identical step count + drift at rounding level" for stages 7b/7c.
> 6. **Before believing a ratio, look at the numerator and denominator separately.**
>    Stage 4g recorded a "1/h law" and a "48x parity swing" for the trapezoidal
>    estimator. Both were one quantity divided by another that was **an order of h
>    smaller and passing through a zero crossing** — the giveaway was a sign change at
>    h=1e-11, which no diverging estimator produces. Neither finding was evidence of
>    anything, and one of them ("4d is blocked on 4g") reordered the whole stage.
> 7. **A tolerance sweep that only goes tighter is not a sweep.** Gate 4b's grid ran
>    reltol 1e-4/1e-5/1e-6 across three circuits and five configurations — 45 runs, which
>    looked comprehensive — and concluded in writing that the shipped default no longer
>    reaches the force-accept path. At **1e-3** it does, and it was taking a step ratio of
>    exactly 10.0 when it did. Loose tolerances are where a step grows into trouble;
>    include at least one.
> 8. **A force-accepted step is an accepted step.** It enters the integrator history like
>    any other, so it belongs in any step-ratio sequence. The first version of gate 4b's
>    probe filtered on the controller's `accept` return and therefore dropped exactly the
>    steps the gate is about — reporting a worst ratio of **97.7** where the truth was
>    **10.0**, and **2.000** ("no violation") on a run that had three. Both errors read as
>    plausible.
> 9. **The review's magnitudes run ~2x optimistic.** Five have now shrunk on measurement
>    (150 TB → 420 GiB; 10.5x → 5.19x; 4.462 ms → 0.234 ms; "2.17 assemblies needed" was not
>    a redundancy figure; the IM3 harness's 10x → 2.17x). Re-measure before quoting
>    `transient_review.md`.
>
> ### Environment
>
> See the `running-tests` memory for the full recipe. Short form: scratch venv with
> `pytest pynose pytest-timeout`, then
> `PYTHONPATH=<repo> MPLBACKEND=Agg <venv>/bin/python -m pytest pycircuit -q -m "" --timeout=400`.
> `PYTHONPATH` is mandatory — a stale root-owned egg shadows the source otherwise.
> Benchmarks: `benchmarks/transient_stage2.py` (perf + drift),
> `benchmarks/transient_stage3.py` (first step), `benchmarks/transient_stage4.py`
> (estimators: `--ratios`, `--hscaling`, `--pi`; and `--forceaccept`, which is the
> gate-4b/4e sweep and **exits non-zero if any accepted step ratio leaves the bound**),
> `benchmarks/transient_review/` (stage 0),
> and **`benchmarks/transient_decisions.py`** — the three measurements that *changed a
> decision or reordered the plan* (`--parity` why 4d is blocked on 4g, `--relref` why D3
> was reverted, `--vabstol` why 0.3a's split was free). Those are the ones worth being able
> to re-run when someone disputes a conclusion.
> Source papers: `/home/andreas/pycircuit_agy/papers/` — **read them from rendered pages**,
> `pdftotext` drops the formulas.

**Status: the plan below was written 2026-07-30, when nothing had run. Stages 0-3 and part
of 4 have since been executed; their `OUTCOME:` lines are filled in from real runs. The
remaining blank `OUTCOME:` lines are work not yet done, and must be filled from a real run,
never predicted.**

Findings and evidence: `doc/transient_review.md`. That document is where an argument about
*whether* something is broken belongs, and it carries the measurements. This one is the
order of work and the gates.

Evidence scripts: `benchmarks/transient_review/` (probes, not gates — see its README).

**Source papers: `/home/andreas/pycircuit_agy/papers/`.** The two that this plan turns on:

- `2014-12--ICECS-Yao-Wang-Roychowdhury-LTE-for-DAEs.pdf` — the source of every `'ywr'`
  formula in `integrator.py` and of `IntegralController`. Table I is the LTE table, eq (22)
  the general form it is derived from, eq (4) the definition of the LTE as a *solution*
  error, Fig. 1 the control flow.
- `2463209.2488904.pdf` — Fang, DAC 2013, cited by `_solve_coupled`. §3.1 is the coupled
  N+1 system, §3.2 the bordered solve, §3.3 the two-sided LTE band, §3.4 the approximate
  Newton method.

**Read them from rendered pages, not `pdftotext`.** Both are formula-dense and text
extraction drops superscripts and table structure; every citation in this document was
taken from 200-dpi renders. The standing rule exists because that failure has produced a
fabricated exponent before.

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
OUTCOME: **Done 2026-07-30. The answer to the stage-9 question is: almost nothing in the
chunking machinery needs to be duplicated, and the duplication has already produced four
divergences, three of which are defects.**

*What genuinely requires static shapes.* Three, and all are legitimate: `CHUNK_SIZE` (the
results/time buffers must be allocated before the `while_loop` traces), the `10000` TLine
ring depth, and the history depth `3`. Two pieces of genuinely traceable control flow are
also legitimate and well done: the Newton `lax.while_loop` and the outer time
`lax.while_loop` with its `lax.cond` accept/reject. **The LTE algebra needs none of it** —
`compute_integration`, `estimate_lte` and `ywr_error_ratio` are pure elementwise arithmetic
on `(g_n, g_nm1, g_nm2, h1, h2)` with no data-dependent control flow, which confirms stage
9(a)'s premise that a shared `_lte_kernels.py` is traceable as-is.

*What is duplicated only because nobody factored it.* `solve` (`:589-728`) and
`solve_batched` (`:449-587`) are the same ~140-line function twice, differing only in a
leading batch axis: breakpoint enumeration, TLine extraction, history initialisation, the
chunk loop, and result assembly are transcribed in both. `jax.vmap` over the whole body is
what `solve_batched` is trying to be. The cost of not doing that is not hypothetical — the
two copies have **drifted apart in four places**, and the drift is invisible because
nothing compares them:

1. **`solve` finds zero breakpoints, always.** `:610` iterates `for elem in
   self.cir.elements` — a dict, so it yields **string keys**, and `hasattr(str,
   'next_event')` is False for every one. `solve_batched` `:476` iterates `.items()` and is
   correct. Verified directly: `[type(e).__name__ for e in cir.elements]` gives
   `['str', 'str', 'str']`. This is plan item 9(d), and the sharp version is that the two
   transcriptions disagree and one of them is right.
2. **`solve_batched` never solves a DC operating point.** `:466-469` is `if not uic: pass`,
   under the comment "We would normally solve DC here." `solve` does the real thing
   (`DC(self.cir).solve().x`). So the batched path silently starts every run from zeros —
   **the stage-1 silent-DC defect, in a second location, and worse, because it does not
   even attempt the solve.** The `x0 = jnp.zeros(n)` it computes at `:467` is then dead:
   `x0_batch` is built separately at `:490`.
3. **`solve_batched` does not initialise TLine history from the operating point.** `solve`
   `:645-657` seeds every ring entry with the DC terminal state; `solve_batched` `:517`
   leaves it at zero under `# Init with zero for now`.
4. `find_idx` (`:106-109`) is defined and never called — dead, and a distraction when
   reading the ring-buffer logic that follows it.

*The `10000` literal.* Seven occurrences: four allocations and, more seriously, three in
modulo arithmetic (`:107`, `:123-124`, `:406`). The comment at `:639-640` says it "covers a
10ns simulation with 1ps steps". Beyond 10 000 accepted steps the ring wraps and
`interp_tlines` interpolates against entries from a previous lap with **no error and no
warning** — `cond_fun` stops at `idx < 9999` and returns whatever is there. A transient with
a `TLine` and more than 10 000 steps therefore returns a plausible wrong waveform. Same
failure class as the two stage 1 exists to remove, so it belongs on that list.

*Tolerances are unreachable by two independent mechanisms, not one.* `JAXTransient`
declares no `parameters`, so it inherits only `Analysis`'s four. Verified:

    JAXTransient(cir, reltol=1e-6)  -> KeyError: 'parameter reltol not in parameter dictionary'
    JAXTransient(cir, iabstol=...)  -> KeyError    JAXTransient(cir, vabstol=...) -> KeyError
    JAXTransient(cir, maxiter=...)  -> KeyError
    JAXTransient(cir, nrsolver='bogus') -> ACCEPTED     JAXTransient(cir, scaler='bogus') -> ACCEPTED

So 9(b) is confirmed exactly, including that `nrsolver` and `scaler` are accepted and
ignored. **But adding the `parameters` list is not sufficient**, which the plan does not
say: `outer_time_loop` declares `reltol=1e-3, abstol=1e-6, maxiter=50` and **neither
`run_chunk` passes them** (`:522-523` and `:672-674`), so the defaults are frozen at the
call site too. `eval_method='gear'` is hard-coded at both of those sites in the same way.
Fixing 9(b) without 9(c) would produce a `reltol` that is accepted and still ignored —
which is worse than the `KeyError`, because the `KeyError` at least tells the truth.

*9(d)'s "fixing either alone converts a silent wrong answer into a hang" — confirmed, with
the mechanism.* `Pulse.next_event` (`func.py:58-82`) opens with `if tmod == 0: return t`, a
fixed point at t=0 and at every multiple of `per`. Verified: for a `VPulse(td=1e-6,
per=3e-6)`, `next_event(0.0)` returns `0.0`. The CPU transient survives this because it
calls `next_event` once per step to decide where to truncate; the JAX paths **enumerate**
in a `while t_event < tend` loop, so the fixed point never terminates. Which path you are
on decides which failure you get: `solve` (broken iteration) finds no breakpoints and
returns a wrong answer; `solve_batched` (correct iteration) reaches the loop and hangs.
Note this also makes `Pulse.next_event` a stage-8 item, not only a stage-9 one.

**Recommendation for stage 9**, in this order: (i) fix `Pulse.next_event`'s fixed point
first, since both other fixes depend on it not hanging; (ii) merge `solve` and
`solve_batched` rather than repairing both — repairing both is how a fifth divergence gets
added; (iii) `_lte_kernels.py` as planned; (iv) `parameters` **and** the call-site
threading together. Add the `10000` overflow to stage 1's list of silent failures.

**0.1b — `semiconductors.py` and `mos.py`, as a device-modelling expert.** Stage 5 adds
junction limiting to `BJT`/`JFET`/`ZenerDiode` on the strength of a proposed diff that was
never tested. Also open: no charge model on any semiconductor (`BJT.C(x) == 0`), the
`Varactor` clamp that makes C fall to zero above the knee instead of extrapolating,
`MOS_ACM.__init__` calling a class not in its MRO, and the absence of any large-signal
MOSFET. Ask: what is the minimum device set that makes CMOS and bipolar transients
expressible, and is the `Semiconductor` base the right seam?
OUTCOME: **Done 2026-07-30. `Semiconductor` is the right seam and should carry the
limiter; the charge gap is the thing that actually blocks CMOS, not the limiter.**

*Is `Semiconductor` the right seam? Yes, and for a reason worth stating.* The base already
enforces the invariant that matters: devices write no `G`/`C` of their own, the model lives
once in `eval_i_pure`/`eval_q_pure`, and the Jacobian is differentiated from it (`:27-37`).
That is exactly the property a limiter needs in order to be safe — a limiter perturbs the
point at which the model is linearised, and a hand-written stamp that disagrees with `i()`
would let Newton converge to a slightly wrong answer *and* hide the disagreement. So
stage 5(a)'s `junctions` class attribute plus a shared `limit()` on `Semiconductor` is the
right shape. **One correction to the plan's framing:** `ZenerDiode` should not get a plain
`pnjlim`. Its reverse breakdown term `-IBV*(exp((-v-BV)/VT) - 1)` is a second exponential
in the *opposite* direction, and a limiter that only knows about the forward junction will
step straight through the breakdown knee. It needs the junction list to carry a direction,
or an explicit breakdown limiter. Reconsider if a measurement shows plain `pnjlim`
converging on a Zener in breakdown — but that should be measured, not assumed.

*The minimum device set for CMOS and bipolar transients.* Ranked by what blocks what:

1. **A charge model on `BJT`.** `BJT` defines no `q`, so `Semiconductor.C` returns a zero
   matrix (`:33-34`) — no depletion capacitance, no diffusion capacitance, no `TF`. A
   bipolar circuit with `C(x) == 0` has no charge storage anywhere in the transistor, so
   its transient is a sequence of DC solves with the wrong dynamics, and every switching
   time it produces is wrong. This is a **bigger obstacle to a credible bipolar transient
   than the missing limiter is**, and the plan does not rank it. The limiter decides
   whether the run converges; the charge model decides whether a converged run means
   anything.
2. **A large-signal MOSFET.** Confirmed absent — `mos.py` contains only small-signal
   subcircuits. Per decision 0.3c this is in scope for stage 10 and sequenced after this
   review and stage 5.
3. Junction limiting on `BJT`/`JFET` (stage 5a) — needed for (1) and (2) to converge.
4. `_expl()` clamping (stage 5b). Note where the overflow actually reaches: these models
   are differentiated by `toolkit.jacobian`, so under a toolkit without autodiff an
   `exp()` overflow in `eval_i_pure` reaches a *central difference* and becomes `nan` in
   the Jacobian rather than merely a large number in `i()`. That is why 5(b) belongs to
   this base class and not to individual devices.

*Confirmed defects, each verified rather than read.*

- **`MOS_ACM` cannot be constructed at all.** `mos.py:104` calls `super(MOS, self)` from a
  class whose MRO is `[MOS_ACM, SubCircuit, Circuit, object]` — `MOS` is not in it.
  Verified: `TypeError: super(type, obj): obj (instance of MOS_ACM) is not an instance or
  subtype of type (MOS)`. **And it is not an ACM model.** Its body is a verbatim copy of
  `MOS` with one difference: `Symbol('kT')` in the noise PSD where `MOS` uses
  `toolkit.kboltzmann * Symbol('T')`. So the class advertising a different model is a copy
  of the one it sits next to, with a worse noise expression, and has never run.
  **Recommendation: delete it.** A thin advertised feature is worse than an absent one, and
  this one is not thin, it is empty. Blast radius, enumerated: **zero references anywhere**
  outside `mos.py` itself and the review documents — no test, no example, no doc page, and
  it is not exported from `pycircuit/circuit/__init__.py`. Deleting it converts a
  `TypeError` at construction into an `ImportError` at import, which is earlier and
  clearer. **The mechanism that let it survive is worth naming:** the only thing that would
  ever have caught it is the doctest in its own docstring (`mos.py:80`), and `pytest.ini`
  configures no doctest collection, so `if __name__ == "__main__": doctest.testmod()` at
  `:131-133` runs only if someone executes the module directly. The test existed and was
  never run — which is a worse state than having no test, because it reads as coverage.
  Reconsider if someone actually wants ACM — in which case it is written from the paper,
  not recovered from this.
- **`Varactor`'s clamp.** `v_eff = minimum(v, 0.99*VJ)` (`:208`) freezes the charge above
  the knee, so `C = dq/dv` falls to **exactly zero** in forward bias rather than
  extrapolating. SPICE's treatment linearises the junction charge above `FC*VJ` and keeps a
  finite, growing capacitance. Zero is the worst available answer: it is not merely
  inaccurate, it removes the state variable, and a Newton step that sees `C = 0` on a node
  that physically has the largest capacitance in the circuit will take a wildly wrong step.
- **`SubCircuit.limit` discards its work** (`circuit.py:1194-1200`), as the plan says, and
  the reason it survived is confirmed: `subx = x[self.elementnodemap[instance]]` is fancy
  indexing and therefore a copy, `element.limit(subx, subx0, epar)`'s return value is
  dropped, and the method returns the **unmodified** `x`. `Diode` does not notice because
  it keeps `_vlim` privately (`elements.py:869-905`) and linearises `G` around it
  (`:907-914`) — so the only limiter in the tree is the one that cannot expose the bug.
  Gate 5-3 is therefore correctly specified and is the important one in stage 5.

*On stage 5(d), whether `limit` should return the limited `x`.* **Recommendation: yes, and
survey it before stage 5 rather than after.** The state-free form is the prerequisite for
the JAX path ever having limiting (device-private `_vlim` cannot be traced), and it is what
makes gate 5-3 expressible at all. But the blast radius is real and the plan is right to
want it surveyed separately. Reconsider if the survey finds a device whose limiter
genuinely needs cross-iteration memory that cannot be carried in `x`.

**0.1c — `shooting.py`, as an RF/steady-state expert.** PSS/PAC are advertised in the
analysis inventory and the review found them unusable beyond toys: `method` Parameter
never read, no limiter, a dense `np.linalg.inv` per timestep inside the shooting Newton,
and a PAC matrix that is ~150 TB at N=137/M=1000. Ask: repair, rewrite against the
`Integrator`/`LinearSolver` seams stage 7 creates, or withdraw and document as
unimplemented? A thin advertised feature is worse than an absent one.
OUTCOME: **Done 2026-07-30. Recommendation: WITHDRAW PAC and document it as
unimplemented; KEEP PSS but rewrite it against the stage-7 seams. They are not in the
same condition and the plan's three options should be applied to them separately.**

*First, a correction to the review's evidence, because the plan cites the number.* The
review says the PAC matrix is "**~150 TB** at N=137, M=1000". **That is wrong by about
365x.** Generated by `benchmarks/transient_review/stage0_1c_pac_memory.py`, not quoted:

    N*M            = 137 000
    (N*M)^2        = 1.8769e10 elements
    L, complex128  = 3.0030e11 bytes = 279.7 GiB
    B, float64     = 1.5015e11 bytes = 139.8 GiB
    total          = 419.5 GiB = 0.410 TiB

The `150` is recognisable: it is `B`'s 1.5015e11 bytes read as 150 **GB** and written as
150 **TB** — a unit slip on the smaller of the two arrays. **The conclusion is completely
unaffected** — 420 GiB is as unallocatable as 150 TB — which is exactly why the error
survived: nothing downstream depended on the magnitude. Recorded here rather than quietly
fixed, per rule 1. It is also a live instance of the standing rule against typing a
measured number into prose; `benchmarks/transient_review/` should have carried a script.
(Incidentally `B = tk.zeros(L.shape)` at `shooting.py:195` passes **no dtype**, so B really
is float64 while L is complex — half the memory, and a real-typed matrix added to a complex
one. It works by promotion, but it is not deliberate.)

*PSS — repairable, and worth repairing.* The core is sound: a shooting Newton on the
monodromy matrix, which is the right algorithm. Three defects, in order of how much they
cost:

1. **`method` is dead and PSS is backward-Euler only.** The Parameter is declared at
   `:63-65` with `default="euler"` and is never read anywhere in the file; `solve_timestep`
   hard-codes `C/dt` and `q(xlast)/dt`. For a *periodic steady-state* solver this is the
   worst possible fixed choice: BE's numerical damping attenuates exactly the limit cycle
   PSS exists to find. An oscillator's amplitude comes out low and a resonator's Q comes
   out understated, both silently and both worse as the step coarsens. Stage 7's
   `Integrator` seam fixes this by construction.
2. **No limiting** (`analysis.fsolve` at `:94` and `:141`, `limiter=None` by default). PSS
   on a *linear* circuit is pointless — the whole point is a nonlinear periodic steady
   state — so "no limiting" means "fails on every circuit it exists for". This is the same
   dependency stage 5 creates, which is why stage 11 should follow stage 5 as well as
   stage 7.
3. **A dense `np.linalg.inv` per timestep inside the shooting Newton** (`:132`,
   `Jshoot = inv(Jf) @ C @ Jshoot`). An explicit inverse is never needed here: the standard
   result in the field (Telichevesky, Kundert & White, DAC 1995, matrix-free Krylov
   shooting) is that only matrix-vector products with the monodromy matrix are required,
   i.e. a sequence of solves. At N=137, M=1000 with 20 shooting iterations this is 20 000
   dense inversions plus 40 000 dense matmuls. Stage 7b's `factor`/`solve` split is the
   seam that makes the correct version cheap to write.

*PAC — withdraw.* Beyond the 420 GiB allocation, which is not a tuning problem but a
consequence of forming the whole `(N*M)x(N*M)` operator densely, **it has never been
validated against anything.** `test_PAC` is `@unittest.skip("Skip failing test")` *and* its
last line is `assert False, "Test should compare with spectre simulation"` — so it is not a
test that broke, it is a test that was never finished. There is no evidence PAC has ever
produced a correct number. Shipping it in the analysis inventory is the "thin advertised
feature" the plan warns about in its purest form.

*What the test suite actually covers, since this is what "unusable beyond toys" rests on.*
Three tests: `test_shooting` compares a **linear** RC against the AC steady state and is
the only one that asserts anything meaningful; `test_PSS_nonlinear_C` calls `pss.solve()`
and **asserts nothing at all** — a smoke test that checks only that it does not raise; and
`test_PAC` is skipped and unfinished. The PSS test runs at M=20 timesteps per period, which
is why the per-timestep `inv` has never been felt.

**Sequencing consequence:** stage 11 is blocked on **stage 5 as well as stage 7**, because
a PSS without a limiter cannot be tested on the nonlinear circuits it exists for. The plan's
dependency graph shows only `0.1c -> 11`; that is not sufficient.

**0.1d — `_solve_coupled`, as a numerical analyst.** The third transcription of the time
loop (`transient.py:429-551`). Read for consistency but its Fang-2013 co-determination
loop and `MAX_LTE_ITERS = 10` were never analysed. It ignores breakpoints, `uic`,
`fixed_timestep` and any injected step controller. Ask: is this worth keeping at all, or
is it a research prototype that should be deleted or moved behind an explicit flag?
OUTCOME: **Done 2026-07-30. Recommendation: DELETE it. It is a research prototype, it is
not the algorithm its comment claims, it buys nothing over `solve()`, and it contains a
livelock that `solve()` does not.** Evidence:
`benchmarks/transient_review/stage0_1d_coupled_livelock.py`.

*It is not a co-determination method, and it is not Fang's approximate Newton either.*
**Updated 2026-07-30 after reading the paper** (`/home/andreas/pycircuit_agy/papers/`
`2463209.2488904.pdf`, G. Peter Fang, "A new time-stepping method for circuit simulation",
DAC 2013; read from 200-dpi renders, not text extraction). An earlier version of this
outcome declined to assert what §3.4 says because the paper was not available. It is
available, and it refutes the comment on both counts.

**What §3.1 actually proposes.** The step size `h_m` is treated as *an independent variable
— an unknown*. The LTE condition `f_lte(v_m, h_m) = ε_m(v_m,h_m) − τ_m = 0` (eq 10) is
solved *together with* the circuit equations as one coupled system of **N+1 nonlinear
equations in N+1 unknowns** (eq 11), via the bordered linear system (eq 12)

    [ J   p ] [ Δv ]   [ −f_ckt ]
    [ q^T d ] [ Δh ] = [ −f_lte ]

with `p = ∂f_ckt/∂h_m`, `q^T = ∂f_lte/∂v_m`, `d = ∂f_lte/∂h_m`. **`_solve_coupled` contains
no `p`, no `q`, no `d`, and no bordered system.** `h` is never an unknown in any solve.

**What §3.4 actually proposes.** The approximate Newton method exists because the coupled
system is sensitive to step-size change during fast transitions. It does two things: (i)
predict the new step from Gear's formula using the LTE from the first stage,
`h^{k+1} = (τ_m / ε_m^{k+1/2})^{1/(n+1)} · h^k` (eq 17); and (ii) — the part that makes it
a *Newton* method rather than a retry — **correct the existing solution instead of
re-solving it**:

    Δv^{k+1} = Δv^{k+1/2} − J^{-1} p (h^{k+1} − h^k)          (eq 18)

a single back-solve against the already-factored `J` and the sensitivity vector `p`. The
paper's own summary: "Without the need to modify the matrix solver of the simulator, the
approximate Newton method is straightforward to implement and carries very little
overhead."

**What the code does.** `:511-526` calls `controller.evaluate_step(...)`, whose
`IntegralController` law `h*safety*(1/err)^(1/p)` is *structurally* eq (17) — and then
`h_curr = h_next; continue`, which **re-solves the whole circuit from scratch**
(`solve_timestep` at `:502`). Equation (18), the entire computational content of §3.4, is
absent: there is no `p`, no sensitivity, no correction. So the one thing borrowed from the
paper is eq (17) — a standard Gear step predictor that is **not distinctive to the paper
and that `solve()`'s controller already applies identically**.

**Two further omissions, both load-bearing in the paper.** §3.3 requires a *two-sided* LTE
band `γ_min·τ_m ≤ ε_m ≤ γ_max·τ_m` (eq 15); the lower bound is what "prevents step sizes
from being unnecessarily small", and §4.1 attributes the paper's headline result — 39%
fewer time points on a Class-D amplifier — to bounding LTE between 0.7 and 3.0. The code
has only an upper bound. §3.3 also requires `|Δh^{k+1}| ≤ η·h^k` (eq 16, typical η = 15%)
to damp step-size change; the code has no such limit, which is a plausible contributor to
the collapse documented above.

**So the citation should be removed regardless of what happens to the code**, and the
divergence was not intended: the dead `analytical_eh` parameter is a vestige of the
`E_h` gradient. **This is now established rather than speculated** — `doc/src/circuit/
time_stepping.rst` documented, in prose, an "Analytical Gradient (`E_h`)" computed as
`E_h = p(E + TRTOL)/h`, an augmented `(N+1)` system, and a "golden window"
`0.7*tau <= eps <= 3.0*tau` attributed to Fang §3.3. **None of the three is in the code.**
So the method was described before it was written, the description was published, and the
parameter is the last surviving connection to it.

That page is itself a stage-1 defect of the same family — a documented feature that does
not exist makes a claim the software does not honour — and it is corrected in the same
commit, with the discrepancy recorded rather than quietly deleted. It is also a clean
instance of the standing rule working in reverse: the prose was written but never checked
against the code, so instead of the explanation catching a defect in the implementation,
the implementation's absence sat undetected behind the explanation.

*It contains a livelock that the standard path does not.* When Newton fails to converge,
`:504-509` shrinks `h_curr` by 0.25 and `continue`s. After `MAX_LTE_ITERS` failures the loop
simply **exits**, and `:528-530` then advance time by the collapsed `h_curr` and append the
previous solution as if it were a new one. Because `h` is restored from `h_next` at `:543`
— which the failure path never updates — the next outer iteration starts again at full step
size. Measured on an RC forced to fail from step 4 onward:

- each outer iteration costs **10 Newton solves** and advances simulated time by
  `h * 0.25^10`;
- predicted `1.0e-6 * 9.5367e-7 = 9.5367e-13 s` per outer iteration, **measured
  9.5367e-13 s, ratio to prediction 1.0000** — the mechanism is confirmed exactly, not
  inferred;
- from the stall point, finishing the run needs **2.10e6 further outer iterations, i.e.
  2.10e7 further Newton attempts**.

It does not raise and does not return; it grinds. On a first-step failure it behaves
differently and no better: it falls through to `:539` and dies with `AttributeError:
'Transient' object has no attribute '_iq'`, because `_iq` is a hidden side channel set
inside `solve_timestep`, which never ran. So the two failure modes are an obscure
`AttributeError` before any step succeeds, and a livelock after one does. **Neither is the
`RuntimeError` at `:508` that was meant to catch this** — that guard only fires if `h_curr`
drops below `minstep = 1e-18`, and `0.25^10` of a microsecond is 9.5e-13, nowhere near it.

*Confirmed dead or ignored inputs.* `analytical_eh` is in the signature (`:429`), forwarded
from `solve()` (`:229`), and **never used in the body**. `fixed_timestep` and `uic` are not
consulted. `next_event` is never called, so breakpoints are ignored entirely. The step
controller is hard-coded to `IntegralController()` at `:486`, so an injected `PIController`
— which `solve()` respects at `:232-235` — is silently discarded. The plan's list is
accurate on every count.

*Why deletion rather than a flag.* It is already behind a flag (`coupled_lte=False`), and
that has not helped: it is a third transcription of the time loop, and the cost of the
second transcription is documented in stage 9 — the Gear2 LTE defect had to be found and
fixed twice. A third copy that is algorithmically identical to the first, is missing four
of its features, and adds a livelock, is a liability with no compensating benefit. Test
coverage is thin enough that removal is cheap: `test_transient_coupled_lte` and
`test_transient_adaptive_vs_coupled` (`test_analysis_transient.py:260,282`) are the only
users besides an untracked scratch file.

**Reconsider if** someone intends to implement genuine co-determination — in which case it
is written against the stage-7 `Integrator`/`LinearSolver` seams from the paper, not
recovered from this. **If deletion is declined**, the minimum is: fix the livelock by
raising on exhaustion instead of falling through, and delete the Fang citation, which the
code does not support.

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
OUTCOME: **GATE PASSED — the premise is refuted, and stage 4 is unblocked.** Measured
2026-07-30 by `benchmarks/transient_review/stage0_2a_integrator_choice.py` on the
leapfrog transient fixture (`cir.n` = 139; the 136 recorded elsewhere is the *symbolic*
build, which has one `VS` where this has two `VSin` and a `BSource`) over the same 2.5 us
window the original comparison used, at the
harness's own tolerances (`vabstol` 1e-9, `reltol` 1e-6, `max_step` 19.53 ns):

| config | steps | wall |
|---|---|---|
| A trapezoidal, breakpoints as shipped | 282 | 63.5 s |
| B Gear2 `'ywr'`, as shipped | 338 | 78.7 s |
| C Gear2 `'classic'`, as shipped | 368 | 76.1 s |
| D trapezoidal, `Sin.next_event` suppressed | 277 | 87.2 s |
| E Gear2 `'ywr'`, suppressed | 334 | 100.3 s |

Breakpoint removal moves the trapezoidal step count by **1.8%** (D/A = 0.982) and Gear2's
by **1.2%** (E/B = 0.988), both far inside the declared 20% band. The ranking that chose
the integrator is **1.199 with breakpoints and 1.206 without** — unchanged. **The recorded
10x stands and does not need re-measuring after stage 4.**

Two things make this a result rather than an absence of one. First, A and B reproduce the
recorded 282-vs-288 and 338-vs-347 from the harness comment, so the probe is measuring the
same thing that comment measured. Second, **the intervention demonstrably did something**:
D is not bit-identical to A (5 fewer steps, waveform differing by 5.4e-10, 2.9e-4 relative),
so a null result here is a measured null and not a monkeypatch that silently failed to
apply — which is the way this gate could most easily have lied.

**What this does NOT clear, stated so it is not over-read.** The gate asked only whether
the *integrator choice* is contaminated. It is not. The trapezoidal `(-1)^n` contamination
itself is untouched by this result: gates 4g-1 (Newton evaluations forced to order 1) and
4g-2 (est/true must stop scaling as 1/h) are separate measurements and remain to be taken
in stage 4. And the window spans roughly one breakpoint interval per tone, so it bounds a
*per-breakpoint* effect and understates a cumulative one; the full IM3 figure is only
re-validated by rerunning the harness, which gate 4g-3 requires anyway.

**0.2b — What does the JAX charge-domain LTE path cost, and should it exist?**
`estimate_lte` + `lte_error_ratio` avoid a per-step `G` evaluation and linear solve, which
is their one genuine advantage over `ywr_error_ratio`. But `lte_abs = 1e-6` is a *voltage*
tolerance applied to a **charge**, so that path currently never rejects a step.
**Gate:** measure the per-step cost of the extra `G` + solve as a fraction of the whole
step, on a circuit of realistic size. Declared decision rule: **if under 10%, delete
`estimate_lte` and make `lte_formula` select an algebra only** (as it does on the CPU);
if over, fix its tolerance to a charge-referenced one and keep both.
OUTCOME: **UNDER the threshold, by a wide margin. Measured 2026-07-30 by
`benchmarks/transient_review/stage0_2b_lte_solve_cost.py` on the leapfrog transient
fixture (n = 139), 2.5 us window.**

| | trapezoidal | Gear2 `'ywr'` |
|---|---|---|
| steps | 282 | 338 |
| `evaluate_step`, all of it | 3.13% | 1.04% |
| **LTE `J^-1` solve** | **2.77%** | **0.89%** |
| LTE `remove_row_col` copies | 0.05% | 0.06% |
| **total avoidable** | **2.82%** | **0.95%** |
| still paid either way | 0.32% | 0.09% |

The wrapper overhead measured +10.5% on one configuration and -1.2% on the other, so
run-to-run noise here is of order 10% — which is worth stating, and which does not
threaten the conclusion: the largest avoidable share is **2.82% against a 10% threshold**,
and the timing method is biased *towards* overstating it (Python call overhead lands on
the measured side). The decision rule fires cleanly.

**A correction to the gate's own premise, and it halves the case before any timing.**
There is **no extra `G` evaluation to avoid.** Both estimators receive the Jacobian as an
argument from the Newton solve that just converged — CPU `IntegralController.evaluate_step(
..., J, ...)` (`stepcontroller.py:32`) and JAX `ywr_error_ratio(i_curr, x_curr, x_last, J,
...)` (`jaxtransient.py:287`) — and neither calls `cir.G`. The charge-domain path's
advantage was never "a `G` evaluation and a linear solve"; it is a linear solve and two
`np.delete` copies. So the premise recorded in the plan was half wrong, and what remained
of it measures at 1-3%.

Also confirmed by inspection: `lte_error_ratio` (`jaxtransient.py:327`) applies
`lte_abs = 1e-6` to a **charge**. With picofarad-scale devices `q ~ 1e-12 C`, so the
tolerance is ~7e-6 against an error of ~1e-15 — a ratio of ~1e-10, and the path never
rejects a step. The plan's claim stands.

**BUT THIS NOW CONFLICTS WITH DECISION 0.3a, AND THE CONFLICT IS REAL — see 0.3d below.**
The decision rule says delete the charge-domain estimator; 0.3a(iii) says make the
criterion charge-referenced. Those pull in opposite directions and the plan did not
anticipate it, because 0.3a was expected to be option (i). Resolving it is a maintainer
decision, not an implementer one, so it is recorded as a new item rather than settled here.

**0.2c — Re-measure the suite baseline on the executing machine.** Both tallies and both
runtimes, `-m ""`, before any change.
OUTCOME: **Measured 2026-07-30 at `d44ab46`. Suite `-m "" --timeout=400`: 734 passed,
6 skipped, 0 failed, 497.69 s.** The tally matches the recorded baseline exactly; the
runtime does not carry the review's ~17% warning, because 497.69 s *is* this box. Every
commit between `e9dd894` and `d44ab46` is doc or benchmark work, which is consistent with
an unchanged tally. Notable in the durations: `test_stress_stiff_rlc_pulse` is **73.14 s**,
not the ~266 s recorded in `architecture.md` P15 — the step-control fix (`e37ddad`, 19x
fewer steps) already took that cost out, so P15's premise is stale and the 20%-regression
rule now has a much smaller worst offender to protect.

Doc build at the same commit: **`build succeeded, 2 warnings.`, 0 `ERROR` lines** — the two
warnings are the pre-existing `pycircuit.post.cds` / `cds.skill` autodoc import failures.
Verified per rule 3 rather than from the exit code: no `exec-rst` block fell back to
rendering its own source (the block source appears **0** times in the built HTML across
`symbolic_poly`, `soe_symbolic`, `lte_dae` and `distortion_limits`) and every page's tables
carry computed cells. A note for later gates: grepping the HTML for `\d+\.\d{4,}` finds
nothing on `symbolic_poly`/`soe_symbolic`, because their generated cells are 3-decimal
timings and integer term counts — the absence of a long decimal is **not** evidence a block
died, so check for the block's own source text instead.

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
DECISION: **(iii) — adopt SPICE's split properly.** The maintainer took the largest option
over the recommendation, on 2026-07-30. Consequences, so this is not re-litigated:

- `vabstol` reverts to being **Newton's x-tolerance only**, shared by `DC` and `Transient`,
  which removes the 10^6 asymmetry between the operating point and the steps after it.
- A `chgtol` Parameter is added (SPICE's default is 1e-14 C) and the LTE criterion on
  charge-carrying rows becomes charge-referenced: `chgtol + reltol*max|q|` rather than a
  voltage tolerance applied to a charge.
- **This absorbs 0.2b's decision rule.** 0.2b asks whether the JAX `estimate_lte` charge
  path should be deleted or have its tolerance fixed; under (iii) the charge flavour is the
  *intended* one on both backends, so the two paths converge rather than diverge. 0.2b's
  measurement is still worth taking — it decides whether the extra `G`+solve is affordable —
  but its "delete `estimate_lte`" branch is now much less attractive, because that path
  becomes the one that matches the criterion instead of the one that violates it.
- **Cost, stated plainly:** this is no longer a stage-1 one-liner. Stage 1 does the
  *separation* (Newton keeps `vabstol`, the controller stops borrowing it); the
  charge-referenced criterion and `chgtol` land in **stage 4**, next to 4c/4d/4g which are
  already rewriting the estimators. Splitting it this way keeps stage 1 behaviour-bounded.
  **Reconsider if** stage 4 slips: option (i) remains a valid intermediate state and is a
  strict subset of (iii), so nothing done for (i) is wasted.

**0.3b — Does `Gear2` keep `'ywr'` as its default?** It was set to `'ywr'` belt-and-braces
when `'classic'` was repaired. The evidence now runs the other way: `'classic'` measures
asymptotically exact (1.000282 against 2/9) while `'ywr'` is 3/4 biased *and* step-ratio
dependent (swinging 16x across ratios against classic's 4x), and end-to-end `'ywr'` needed
57 rejections and 4 force-accepts where `'classic'` needed 29 and 0.
**Recommendation: flip the default to `'classic'` in stage 4**, and delete the "5/6"
claim in the `integrator.py` comment, which is measurably wrong.
DECISION: **Keep `'ywr'` for now; re-decide at gate 4f.** The maintainer declined the flip
on 2026-07-30, on the grounds that the evidence against `'ywr'` was gathered *with the
controller defects still in place* — 4a's unstable PI gains, 4b's 10x force-accept growth
and 4e's inverted order-drop guard all inflate rejection counts, and `'ywr'` being
step-ratio dependent is exactly the property that a cycling step size punishes hardest. The
57-vs-29 rejection comparison therefore cannot separate "the estimator is worse" from "the
estimator is more exposed to a controller that is misbehaving".

**What this changes:** 4f stops being a formality. Gate 4f must be run *after* 4a-4e land
and must record rejections and force-accepts for all four integrator/formula combinations
on the same circuit; the default is then set from that table, not from the pre-repair
numbers.

**SUPERSEDED 2026-07-31 by 4g(b) and 4i, and the question has changed shape.** There are no
longer four distinct integrator/formula combinations to compare: both second-order
estimators now take their `q'''` from a shared helper that does not read `lte_formula`, so
`Trapezoidal('ywr')` and `Trapezoidal('classic')` are identical, and so are the two Gear2
variants — except on the single fallback step of a run, where four charges do not yet
exist. The decision this item deferred ("which `lte_formula` should be the default")
therefore has no measurable content left. **The live question is whether the parameter
should exist at all**, which is a maintainer's decision about a public knob rather than an
implementer's about a default, and it is recorded as such in the resume block. The `'classic'` asymptotic-exactness measurement (1.000282 against 2/9) is not in
dispute and is not what was declined.

**Not deferred:** the `"5/6"` claim in the `integrator.py` comment is measurably wrong
regardless of which default wins, so it is deleted in stage 4 either way.

**0.3c — Scope of stage 10.** The missing-analyses list (DC sweep, `.ic`/`.nodeset`,
netlist import, large-signal MOSFET) is weeks of work and is a product decision, not an
engineering one.
DECISION: **In scope: DC sweep, `.tran` output control, `.ic`/`.nodeset`, and a
large-signal MOSFET.** Out of scope for stage 10: the SPICE-subset netlist reader and
waveform export (raw/PSF/CSV) — both are interop rather than analysis, and neither blocks a
circuit being simulated from Python. **Reconsider if** the MOSFET work creates a demand for
model cards, since a `.model` parser is most of a netlist reader's value at a fraction of
its scope.

**Ordering consequence, and it is not the plan's original order.** The large-signal MOSFET
is a *device* item sitting in an *analysis* stage, and it overlaps review **0.1b**
directly — 0.1b asks what the minimum device set is that makes CMOS and bipolar transients
expressible, and whether `Semiconductor` is the right seam. Adding a MOSFET before that
review answers is how a device gets bolted to the wrong base class. **The MOSFET is
therefore sequenced after 0.1b and after stage 5**, which is where the limiting and
`_expl()` clamp it will depend on actually land; a large-signal MOSFET without junction
limiting would not converge and would produce exactly the class of silent failure stage 1
exists to remove. The other three stage-10 items have no such dependency.

Ranked within the stage: (1) DC sweep, (2) `.tran` output control, (3) `.ic`/`.nodeset`,
(4) large-signal MOSFET — 1 and 2 are the cheapest and most visible, 3 is what makes
oscillators and latches startable, 4 is the largest and the most gated.

**0.3d — NEW, raised by stage 0 itself on 2026-07-30: decisions 0.3a and 0.2b now
contradict each other.** This item did not exist when the plan was written; it appeared
because 0.3a was answered (iii) where the plan expected (i).

**The contradiction.** 0.2b measured the `J^-1` mapping at 1-3% of a step, so its declared
rule fires: delete the charge-domain estimator and always map the truncation error into
the solution domain, where it is compared against a **voltage/current** tolerance. 0.3a(iii)
adopts SPICE's split, whose whole point is that the truncation error on charge-storing rows
is compared against a **charge** tolerance (`chgtol + reltol*max|q|`) — and SPICE does that
*without* a `J^-1` solve, because once the criterion is charge-referenced the mapping into
volts is not merely unnecessary, it is incoherent: you would map to volts and then apply a
tolerance in coulombs.

So these are not two independent knobs. **Mapping through `J^-1` and a charge-referenced
tolerance are alternative formulations of the same criterion**, and stage 4 must pick one:

- **(A) YWR, solution-domain.** Keep `J^-1`, keep `reltol`/`vabstol`/`iabstol` as the LTE
  tolerances, delete the charge estimator on both backends. 0.2b says it is affordable. It
  is the more principled choice for a DAE, where a charge error on one node shows up as a
  voltage error on another. It does **not** need `chgtol`, so 0.3a collapses back to option
  (i) — which the maintainer explicitly did not choose.
- **(B) SPICE, charge-domain.** Drop the `J^-1` mapping, add `chgtol`, compare the charge
  LTE directly. This is 0.3a(iii) as chosen, it is the standard every other simulator
  implements, and it makes the JAX `estimate_lte` path the *correct* one rather than the
  one to delete — reversing 0.2b's rule rather than contradicting it, since 0.2b's rule was
  written on the assumption that the charge path's only justification was cost.
- **(C) Both, selected by `lte_formula`,** with each given the tolerance flavour it
  requires. Honest, and it is roughly where the code already is — but it means maintaining
  two criteria and two tolerance sets, which is what stage 9 is trying to reduce.
- **(D) YWR as the criterion, with the charge check as a guard.** Added 2026-07-30 at the
  maintainer's suggestion, and it is now the recommendation — see below. Not the same as
  (C): there is one criterion, not a choice of two, and the charge quantity is used to
  bound how far the step may collapse rather than to decide accuracy.

### Why (D), and why the two criteria fail in *complementary* regimes

This is the part that makes (D) principled rather than belt-and-braces.

**What the conversion factor physically is.** YWR's Table I inverts
`(q̇_x + β₀h·ḟ_x) = C + β₀h·G`. As `h → 0` that tends to **`C`, the capacitance matrix**, so
the mapping is precisely "charge error → voltage error, via the local capacitance". SPICE's
charge criterion omits it, which amounts to assuming the charge-to-voltage conversion is
uniform across nodes. It is not: the same charge error is a far larger voltage error on a
1 fF node than on a 1 nF node. That mis-weighting is a plausible reading of the "problems in
certain cases" YWR allude to without enumerating.

**Where YWR fails instead.** `C + β₀h·G` is near-singular exactly when some direction has
neither appreciable capacitance nor conductance — a floating or weakly-grounded node. There
the conversion is unbounded, the mapped LTE is dominated by the null direction, and since
`h` does not appear in that direction's error, **shrinking the step cannot reduce it**. The
failure mode is step-size collapse chasing an error the controller cannot fix — a livelock,
not a wrong answer.

So: the charge criterion is wrong when node capacitances are badly spread; the voltage
criterion is wrong when the circuit is near-singular. Neither dominates, and each is
diagnostic of the other's failure.

**The mechanism, and it needs no condition estimator.** Both quantities are already in hand:
`Eg` is computed before the solve and `lte = J⁻¹Eg` after it. Their ratio
`‖lte‖/‖Eg‖` is a *lower bound on* `‖J⁻¹‖` **along the direction that actually matters** —
better for this purpose than a `gecon` estimate, which gives a worst case over all
directions. Measured for reference: `dgecon` off existing factors is ~38% of a
factor-plus-solve at n=139, i.e. affordable but unnecessary. The ratio is free.

**Proposed rule:**

    err_v = max|lte| / etol_v      # YWR, etol_v = TRTOL*(reltol*max|x| + vabstol/iabstol)
    err_q = max|Eg|  / etol_q      # charge, etol_q = TRTOL*(reltol*max|q| + chgtol)

    accept if err_v <= 1                                  # the normal path, unchanged
    if err_v > 1 and err_q <= 1:                          # the guard
        accept, and EMIT a diagnostic naming the amplification ratio,
        the worst row, and its node name via cir.get_node_name

The guard is a floor on step collapse, not a second accuracy test — it fires only when the
charge error says the step is fine and the mapping disagrees, which is the signature of the
near-singular case and of nothing else. **It must be logged, not silent**, or it becomes
the same failure as the `except` at `stepcontroller.py:59-62`: a mode switch nobody can see.
Routing it through stage 6's statistics object is the natural home, and a near-singular `J`
during a transient is a reportable circuit condition in its own right.

**This partly rehabilitates 0.3a(iii):** `chgtol` comes back, because `etol_q` needs an
absolute charge floor. What does *not* come back is the charge criterion as the primary
accuracy test. So the maintainer's original instinct survives in the guard.

**The cost, stated plainly.** Two tolerance sets to document and tune, which is what stage 9
is trying to reduce, and a second way to be wrong: if `etol_q` is set too loose the guard
becomes an escape hatch that accepts genuinely bad steps. **Gate this** — the guard's firing
count must appear in the statistics, and a run where it fires often is a run whose result is
suspect, not a run that was rescued.

**What to measure before committing to (D)** — it is a design, not yet a result:
1. Does the near-singular case actually arise in pycircuit's circuits? Stage 6's floating-node
   work will produce the test cases; if the guard never fires on any real circuit, (D)
   collapses to (A) and the charge path should go after all.
2. Does the amplification ratio actually separate the two regimes cleanly, or is there a
   continuum where neither criterion is trustworthy? That is the assumption (D) rests on.

### What Spectre actually does — and the root cause it exposes

Checked 2026-07-30 against the *Spectre Circuit Simulator Reference*, Product Version
19.1 (January 2020), transient analysis parameter list, read from a **rendered page 419**
rather than text extraction. Two entries settle it:

> **38 `relref`** — "Reference used for the relative convergence criteria. The default is
> derived from `errpreset`. Possible values are `pointlocal`, `alllocal`, `sigglobal`, and
> `allglobal`."
>
> **39 `lteratio`** — "**Ratio used to compute LTE tolerances from Newton tolerance.** The
> default is derived from `errpreset`."

So Spectre has **one** tolerance set — `reltol`, `vabstol`, `iabstol` — and derives the LTE
tolerance by *multiplying the Newton tolerance by `lteratio`*. It does **not** carry a
separate absolute tolerance for the step controller. `lteratio` defaults to **3.5**
(liberal and moderate `errpreset`) and **10.0** (conservative). pycircuit's `TRTOL = 7.0`
is structurally the same knob — `etol = TRTOL*(reltol*ref + abstol)` — sitting at an
unexplained value between Spectre's two.

**`relref` is the parameter pycircuit does not have, and its absence is why `vabstol` was
raised in the first place.** The four modes differ in what the *relative* term is measured
against (§"Analysis Statements", PSS discussion, same manual):

- `pointlocal` — each node against that node alone;
- `alllocal` — each node against the largest value **that node** has taken over all past time;
- `sigglobal` — each signal against the maximum over **all signals at any previous time**;
- `allglobal` — as `sigglobal`, and also references each node's KCL residue to its own history.

pycircuit hard-codes `reltol * max(|x_curr|, |x_last|)` — a two-point `pointlocal`. That is
the mode in which a node carrying no signal has `reltol*ref -> 0`, so its tolerance
collapses to `abstol` and the controller ends up chasing numerical noise on idle nodes.
**That is exactly the failure recorded in the `vabstol` comment** ("most of that circuit's
nodes carry no signal, so etol degenerated to TRTOL*abstol on numerical noise"), and it is
what raising the absolute tolerance a millionfold was working around. Spectre's default for
every `errpreset` in the shooting table is `sigglobal` — the mode that cannot degenerate
that way, because a quiet node is referenced to the largest signal in the circuit.

**Consequence for decision 0.3a, stated plainly: the options it offered missed the axis
that matters.** (i), (ii) and (iii) all argued about *which absolute tolerances exist*; the
actual defect is *what the relative term is referenced to*. `lte_vabstol` is therefore a
parameterisation of a workaround, not a fix — it makes the workaround explicit and
separates it from Newton, which is real progress and fixes a real regression, but it is not
where this should end up.

**Amendment to stage 4: add `relref`.** Implement at least `pointlocal` (current behaviour)
and `sigglobal` (Spectre's default), defaulting to `sigglobal`. **Gate:** with
`relref='sigglobal'`, `lte_vabstol` must be returnable to 1e-12 without the step-count
collapse that motivated 1e-6 — measured on the leapfrog, which is the circuit that produced
the original 5.4x. **If that gate passes, delete `lte_vabstol` and `lte_iabstol`** and let
`lteratio` derive the LTE tolerance from `vabstol`/`iabstol` as Spectre does; the two-set
design then has no remaining justification. **Reconsider if** the gate fails, in which case
the separate set is doing work that `relref` alone cannot, and it stays.

Also worth aligning while there: `TRTOL = 7.0` should be renamed `lteratio`, exposed as a
Parameter rather than a local, and its value chosen against Spectre's 3.5/10.0 rather than
left at a number with no recorded derivation.

### Implementation note: `Eg` is a current, not a charge

This will bite otherwise. pycircuit's `Eg` is **not** the paper's numerator. YWR's Table I
numerator is `(second difference of g) · h`, which has units A·s = **coulombs**, divided by
`(C + β₀h·G)` in farads to give volts. pycircuit absorbs that `h` into `J` instead — because
`J = G + β₀'C/h` and `β₀h·J = C + β₀h·G` — so `integrator.py`'s `lte` return
(`-0.5*(gn - gn_1)`, `-(1/6)*(gn - 2gn_1 + gn_2)`) is a **current**, and `J^{-1}Eg` is volts
only because `J` carries the `1/h`. Both are self-consistent; they are not interchangeable.

**Consequence for the guard:** `chgtol` is a *charge* tolerance, so it cannot be applied to
`Eg` as returned. The guard's quantity must be rescaled to charge —
`Eg_charge = β₀ · h_curr · Eg` with the same `β₀` the integrator used (`1` for Euler,
`0.5` for trapezoidal, `alpha0`-derived for Gear2) — or, equivalently and more safely, have
`compute_lte` return the charge-flavoured numerator alongside the current-flavoured one so
the scaling lives with the integrator that knows its own `β₀` rather than in the controller.
**Prefer the second.** Getting this wrong yields a guard whose threshold is off by a factor
of `h`, which on this codebase means a guard that either never fires or always does — and
either way looks like it is working.

~~**Recommendation: (B).**~~ ~~**(A).**~~ **RECOMMENDATION IS NOW (D)** — YWR as the
criterion, with the charge check as a step-collapse guard. This item has moved twice in one
day, both times on evidence rather than reflection: (B) → (A) on reading the YWR paper, then
(A) → (D) on the maintainer's observation that the two checks could both be kept to cover
the near-singular case. **(D) is a superset of (A)**, so nothing about the (A) argument
below is withdrawn; the guard is added on top of it. The reasoning for preferring the
solution-domain criterion follows, and the guard is specified above under "Why (D)".

The `'ywr'` formulas pycircuit implements come from Yao, Ye, Wang, Wang & Roychowdhury,
*"An Efficient Time Step Control Method in Transient Simulation for DAE System"*, ICECS
2014 (`/home/andreas/pycircuit_agy/papers/2014-12--ICECS-Yao-Wang-Roychowdhury-LTE-for-DAEs.pdf`;
read from 200-dpi renders).

**Precisely what the paper argues, because an earlier draft of this section overstated it.**
It does **not** mention `chgtol`, and it is not an argument against charge tolerances as
such. Its target is the *estimator*: existing methods derive the LTE for an ODE
`dx/dt + f(x) + b = 0`, which (§II-C) "is without the `q()` term", so "the traditional time
step control methods are just an approximation for the circuit simulator and there are
problems in certain cases". The contribution is the missing conversion — a charge-domain
divided difference is **not** the LTE, it is the LTE pre-multiplied by `(q̇_x + β₀h·ḟ_x)`.

So SPICE's scheme is **internally consistent** — a charge error against a charge tolerance
is dimensionally fine — and YWR's claim is that it is an *approximation* of the DAE LTE,
not that it is incoherent. The argument lands on the tolerance flavour only by implication,
because a charge tolerance is the natural companion of the estimator being criticised. That
is a weaker claim than "the paper is against `chgtol`" and it is the one the evidence
supports.

And the paper's LTE is unambiguously **solution-domain**: eq (4) defines
`ε_T(t_n) = x_n − x*_n`, the error in the vector of **node voltages and branch currents** —
not in charge. The `(q̇_x − β0 ġ_x)^{-1}` factor in Table I is precisely what converts the
charge-domain residual into that solution error. So a voltage/current tolerance is the
*correct* companion to these formulas, and applying a charge tolerance to them — which is
what JAX's `lte_error_ratio` does with `lte_abs = 1e-6` — is not merely a units slip; it is
the approach this paper was written to replace.

**Verified, not assumed: pycircuit's `'ywr'` implementation is faithful to Table I.** The
paper gives TRAP as `ε = −(1/12)(q̇_x + 0.5h ḟ_x)^{-1}(g_n − 2g_{n−1} + g_{n−2})h`. With
`(q̇_x + 0.5h ḟ_x) = C + 0.5hG = 0.5h(G + 2C/h) = 0.5h·J` — and `J = G + 2C/h` is exactly
what `TrapezoidalIntegrator.compute_derivatives` builds (`geq = C/h/0.5`) — the `h`
cancels and leaves `Eg = −(1/6)(second difference)` with the controller applying `J^{-1}`.
That is `integrator.py:123` verbatim. The same check passes for Euler/GEAR1. The
implementation is sound; the question is only which criterion it should be judged against.

**Therefore the criterion is YWR's.** Keep the `J^{-1}` mapping and keep
`reltol`/`vabstol`/`iabstol` as the LTE tolerances. 0.2b already established that the
mapping costs 1-3%, so the only argument that ever favoured the charge path as the *primary*
criterion — cost — is measured and does not hold. Under (D) the charge quantity is retained,
but as the guard described above, not as the accuracy test.

**This does not discard 0.3a(iii); it splits it.** Option (iii) bundles two separable
things, and under (D) both survive in modified form:

1. **Separate Newton's x-tolerance from the LTE tolerance.** Uncontroversial, independent
   of this argument, and the actual defect 0.3a was raised to fix — Newton's node
   convergence was loosened 10^6 unmeasured and `DC.vabstol` still disagrees with
   `Transient.vabstol` by the same factor. **Keep this, in stage 1, exactly as decided.**
2. **Make the LTE criterion charge-referenced and add `chgtol`.** The charge-referenced
   *criterion* is what YWR argues is an approximation, so it does not become the accuracy
   test — but **`chgtol` is still needed**, because (D)'s guard requires an absolute charge
   floor for `etol_q`. So the parameter lands as decided; what changes is its role.

So stage 1 is unaffected by the change, and stage 4 gains the guard and its diagnostic
rather than a second accuracy criterion.

**A further argument for deciding this explicitly, found while writing the above.**
`stepcontroller.py:59-62` is

    try:
        lte_reduced = toolkit.linearsolver(J_reduced, Eg_reduced)
    except Exception:
        lte_reduced = Eg_reduced

On any solve failure it **silently falls back to the raw charge-domain `Eg`** and then
compares it against `etol`, which is a *voltage* tolerance. So the CPU controller is
currently neither (A) nor (B): it is (A) with an unlogged (B) fallback carrying exactly the
units mismatch that makes JAX's `lte_error_ratio` never reject a step. Whichever option is
chosen, this `except` must go — under (A) a singular `J` is a diagnosable condition for
stage 6, not a reason to switch error criteria without saying so. It also means the two
formulations are already tangled in one code path, which is the strongest practical reason
not to leave 0.3d open.

**Reconsider (B) if** either: a measurement shows the solution-domain criterion
mis-controlling a circuit where `J` is near-singular, since `J^{-1}` then amplifies the
charge residual unboundedly and the mapped LTE becomes noise — a floating or
weakly-grounded node is the case to try, and it is the real risk of (A); or interoperating
with SPICE tolerances becomes a requirement, in which case `chgtol` is needed for
compatibility whatever its merits. Note the honest limit of the evidence: ICECS 2014 is a
4-page paper and `chgtol` is what every production simulator ships, so this is an argument
from derivation, not from an industry track record.
DECISION: **(D), decided by the maintainer 2026-07-30.** YWR's solution-domain criterion is
the accuracy test; the charge check is retained as a step-collapse guard for the
near-singular case; `chgtol` is added to serve the guard. **Stage 0 is now fully exited.**

**Work this creates, by stage:**

- **Stage 1** — 0.3a's separation, as originally decided: `vabstol` reverts to Newton's
  x-tolerance only and is shared by `DC` and `Transient` (removing the 10^6 asymmetry);
  the controller gets its own `lte_vabstol`. No charge anything in stage 1.
- **Stage 4** — add `chgtol` (SPICE's default is 1e-14 C); have `compute_lte` return the
  charge-flavoured numerator alongside the current-flavoured `Eg`, per the units note above;
  implement the guard; **delete the bare `except` at `stepcontroller.py:59-62`**, which is
  the unlogged half-(B) fallback this decision replaces.
- **Stage 6** — the guard's diagnostic (amplification ratio, worst row, node name via
  `cir.get_node_name`) and its firing count in the statistics object. **The guard is not
  done until it is visible**; an invisible guard is the defect it was introduced to fix.
- **Stage 9** — unchanged from 0.2b's rule: **`estimate_lte` still goes.** (D) does *not*
  resurrect it. The guard needs the charge-domain residual that `ywr_error_ratio` already
  forms as an intermediate before its solve, not `estimate_lte`'s separate q-history
  algebra. `lte_error_ratio`'s misapplied `lte_abs = 1e-6` becomes the properly-named
  `chgtol` in the guard. So `_lte_kernels.py` factors the YWR algebra and exposes both
  flavours of its numerator.

**Gate 4-D (new).** The guard must be shown to fire on a circuit that needs it and to stay
silent on one that does not. Declared success: on a floating-node circuit the guard fires,
the step size does **not** collapse, and the diagnostic names the offending node; on the
review's benchmark circuit it fires **zero** times and the waveform is unchanged from the
pure-(A) path. **If it fires on the benchmark circuit, `etol_q` is too tight and the result
is not trustworthy** — that is a failure of the gate, not a tuning note.
OUTCOME:

**Stage 0 exit criterion:** every OUTCOME above filled, every DECISION answered. Stages
1-3 may start before 0.1a-d land (they touch none of that code); **stage 4 is blocked on
0.2a, stage 5 on 0.1b, stage 7 on nothing, stage 9 on 0.1a, and stage 11 on 0.1c.**

## Stage 0 status, 2026-07-30

**STAGE 0 IS EXITED. Filled: 0.1a, 0.1b, 0.1c, 0.1d, 0.2a, 0.2b, 0.2c. Answered: 0.3a,
0.3b, 0.3c, and 0.3d — the item stage 0 raised itself — decided (D) on 2026-07-30.**

**Every stage is now unblocked.** Stage 4's blocker was 0.3d and it is closed; 0.2a passed,
so stage 4 does not re-measure the integrator choice either. Recommended order remains
1 → 2 → 3 → 4 → 5 → 6 with 7 in parallel, per the plan's own advice that stages 1-3 are the
ones to do if only one thread is available.

**Updated later the same day, after the source papers were located** at
`/home/andreas/pycircuit_agy/papers/`. Three items changed on evidence rather than on
reflection, and the changes go in both directions:

- **0.1d hardened.** The Fang citation can now be checked, and the code implements neither
  §3.1's coupled N+1 system nor §3.4's approximate Newton — eq (18), the whole
  computational content of §3.4, is absent. The earlier outcome declined to assert this;
  that hedge is withdrawn.
- **0.3d's recommendation REVERSED, (B) → (A).** The YWR paper exists to argue that the
  charge/ODE-flavoured criterion is "an approximation", and its eq (4) defines the LTE as a
  *solution* error. Recommending (B) was a recommendation made in ignorance of the paper
  the code's own formulas come from. 0.3a(iii) splits cleanly: keep the Newton/LTE
  tolerance separation, drop `chgtol`.
- **4d confirmed and 4a-bis added.** YWR's Table I TRAP entry is a uniform-grid formula,
  which is exactly the defect 4d describes; and both papers use a two-threshold controller
  where pycircuit has one, which is a new candidate explanation for 0.3b's rejection counts.

Also verified in the same pass, and worth recording as a positive: **the `'ywr'`
implementation is faithful to Table I.** The `−(1/6)` in `integrator.py:123` follows from
the paper's `−(1/12)` once `(q̇_x + 0.5h ḟ_x) = 0.5h·J` is substituted, and `J = G + 2C/h`
is what the trapezoidal companion actually builds. The formulas are right; only the
criterion they are judged against is in question.

Stage 0 is therefore *not* exited, and the one thing left is a maintainer decision that
did not exist when the plan was written. **Stages 1, 2, 3 and 7 are unblocked and may
start now** — none of them depends on 0.3d. **Stage 4 is blocked on 0.3d**, not on 0.2a,
which passed.

What stage 0 changed about the plan, beyond filling blanks:

1. **0.2a's premise was refuted** — the integrator choice is not contaminated, the 10x
   stands, and stage 4 does not have to re-measure it. This is the largest single saving
   stage 0 produced.
2. **0.2b's premise was half wrong** — there is no `G` evaluation to avoid, only a solve,
   and the solve costs 1-3%.
3. **0.3a's answer created 0.3d**, a contradiction between two items the plan treated as
   independent.
4. **Three new defects were found that the four-lens review did not have**, all in the
   silent-wrong-answer class stage 1 exists for: `solve_batched` never solves a DC
   operating point at all (0.1a); the JAX TLine ring wraps past 10 000 steps and
   interpolates against a previous lap without complaint (0.1a); and `_solve_coupled`
   livelocks on non-convergence at a measured `h*0.25^10` per 10 Newton solves (0.1d).
   **The first two belong on stage 1's list.**
5. **The review's "~150 TB" PAC figure is wrong by ~365x** (0.1c); the real number is
   420 GiB. The conclusion it supported is unchanged.
6. **Two deletions are recommended** rather than repairs: `_solve_coupled` (0.1d) and
   `MOS_ACM`, which has never been constructable and is not an ACM model (0.1b).
7. **Stage 11 is blocked on stage 5 as well as stage 7** (0.1c) — a PSS with no limiter
   cannot be tested on the nonlinear circuits it exists for.

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

**Added by stage 0, 2026-07-30.** Three more members of the same class, found by the
0.1a and 0.1d reviews. They are here rather than in stages 9 and 4 because this is the
stage defined by "a run can quietly report a circuit that was never biased", and each of
them does exactly that.

(e) **`JAXTransient.solve_batched` never solves a DC operating point.**
`jaxtransient.py:466-469` is `if not uic: pass` under the comment "We would normally solve
DC here", so every batched run silently starts from zeros. This is defect (a) again in a
second file, and worse: (a) at least attempts the solve before substituting zeros. Same
fix, same message, same `uic` escape. The `x0 = jnp.zeros(n)` at `:467` is dead and goes
with it.

(f) **The JAX TLine ring buffer wraps silently past 10 000 steps.** The `10000` literal
(`jaxtransient.py`, seven occurrences, three of them in modulo arithmetic) is a fixed ring
depth; `interp_tlines`' `cond_fun` stops at `idx < 9999` and interpolates against whatever
entry is there, which beyond 10 000 accepted steps is from a previous lap. A transient with
a `TLine` and more than 10 000 steps returns a plausible wrong waveform with no error and no
warning. **Minimum fix for this stage: detect the wrap and raise.** Sizing the buffer
properly, or replacing it, belongs to stage 9.

(g) **`_solve_coupled` livelocks instead of failing.** Measured: 10 Newton solves buy
`h*0.25^10` of simulated time, forever, because `h` is restored to full size each outer
iteration (`transient.py:543`) while the failure path never updates it. See 0.1d.
**Recommendation: delete `_solve_coupled` and the `coupled_lte` flag entirely**, per 0.1d's
outcome — it is a third transcription of the time loop that is algorithmically identical to
`solve()`, is missing four of its features, and adds this. If deletion is declined, the
minimum is to raise on loop exhaustion rather than fall through.

**CORRECTION, 2026-07-30.** An earlier version of this paragraph said deletion removes two
tests. **That undercounted by a factor of six.** The real blast radius, enumerated:

- `test_analysis_transient.py::test_transient_coupled_lte` and
  `::test_transient_adaptive_vs_coupled` — 2 tests, deleted outright.
- **`test_analysis_transient_stress.py` — all 10 tests**, every one of which calls
  `_compare_methods` (`:13-29`), which runs the circuit *twice*, once per path. This is the
  whole `slow`-marked block.
- `doc/src/circuit/lte_dae.rst:80` and `doc/architecture.md:717` both document the flag.

So 12 tests, not 2. **This changes the shape of the decision and it is recorded here rather
than absorbed**, because the undercount would have turned up as an unexplained tally drop at
gate 1-4 — which is exactly the kind of surprise the gate exists to prevent.

**It does not change the recommendation, and arguably strengthens it.** Two observations:

1. What `_compare_methods` actually asserts is weak. It runs both paths, `warnings.warn`s
   (does not fail) if the step counts differ by more than 3x, and returns both results; the
   tests then assert on each path *independently* — final voltage within tolerance. There is
   almost no differential content. Drop the coupled run and `_compare_methods` becomes a
   one-line `_run`, every existing assertion on `res_adapt` survives verbatim, and the only
   thing lost is a warning nobody reads.
2. If 0.1d is right that the coupled path is algorithmically a rejection loop like
   `solve()`, then this file has been comparing the standard controller against a near-copy
   of itself with different constants — which is a weak test *and* a misleading one, since
   `architecture.md:717` records "comparing the two step controllers is the point of the
   file."

**And it makes P15 partly moot:** the slow block runs every stiff circuit twice, so deleting
the second run roughly halves it. `test_stress_stiff_rlc_pulse` is the largest single test
in the suite at 73 s (0.2c); this is the cheapest available route to that cost, and it does
not weaken a single assertion — which is the objection P15 raises against marking it slow.

**Docs in the same commit:** a section in `doc/transient.rst` (or the transient module
docstring) stating what happens when the operating point fails and how to ask for a
zero start deliberately.

**Gate 1-1 (the failure is visible).** A circuit with no DC solution must raise, with a
message that names the circuit condition. Declared success: it raises; the message
contains the word `uic`; and no waveform is returned.
OUTCOME: **PASSED.** `test_transient_repairs.py::test_gate_1_1_failed_operating_point_raises_and_names_uic`
— two opposing current sources into a node with no DC path raise `NoConvergenceError`, the
message contains `uic`, and `tran.result` is `None`. A companion test
(`..._uic_is_a_working_escape`) checks the escape actually works. Implemented as
`Transient._solve_operating_point`, which also covers item (d): the inner `DC` is built
with the transient's toolkit, epar, tolerances, solver and scaler, forwarding only the
parameters `DC` declares so a future divergence in either parameter list cannot raise.

**Gate 1-2 (epar actually arrives).** Instrument a device's `i()` during a transient with
`epar.T = 400`. Declared success: the device sees T = 400 in **both** the inner DC and
every transient step. Currently it sees 300 in both.
OUTCOME: **PASSED.** `test_gate_1_2_epar_reaches_devices_in_dc_and_every_step` — a
`R` subclass recording `epar.T` on every `i()` call sees **400 K on every evaluation**,
with no 300 K left. `self.epar` is now threaded through all eleven `cir.C/q/i/G/u` call
sites in `transient.py` (both `solve` and `_solve_coupled`).

**Gate 1-3 (bypass is connected).** `bypass=True, bypasstol=1e-2` and `bypass=False` must
produce *different* step/evaluation counts. Currently they are identical because
`bypasstol` is missing from `defaultepar` and every device takes its `except
AttributeError` branch.
OUTCOME: **PASSED, with two corrections to the gate as written.**

*First correction: the gate's observable was wrong.* Counting `Diode.i()` **calls** cannot
detect bypassing, because the bypass test lives *inside* `i()` and skips the exponential,
not the call. Measured that way the counts are identical whether or not bypass works, so
the gate as specified would have passed vacuously before the fix and after it. The working
observable is `toolkit.exp` evaluations, which in the test circuit only the diode uses.
Measured: **165 model evaluations with `bypass=False`, 42 with `bypass=True,
bypasstol=1e-5`** — a 3.9x reduction, with the step count identical (21) and the endpoint
moved 1.7e-8 V, well inside the 1e-5 V tolerance that licensed it.

*Second correction: `bypasstol=1e-2`, the value the gate specifies, does not work.* With
bypass genuinely connected, **1e-2 and 1e-3 both prevent the inner DC from converging** on
a plain diode circuit — the model is frozen across voltage steps far larger than a
junction's own scale, so Newton works from a stale Jacobian. 1e-4 and below converge. This
is correct behaviour for an absurd tolerance and is now pinned by
`test_gate_1_3_a_loose_bypasstol_is_not_silently_tolerated`, which asserts it *raises*
rather than returning a waveform. **A negative result worth keeping: `bypass` has been
dead for long enough that no value of `bypasstol` was ever validated, and the plan's
suggested value is outside the working range.**

`bypasstol` is also added to `defaultepar` with default **-1.0** ("never bypass"), which
is what `bypass=False` should mean; the three devices previously fell through to a
hard-coded 1e-12, so a little unrequested bypassing always happened.

**Gate 1-4 (blast radius).** Full suite `-m ""`. Declared success: 734/6/0, **or** a list
of failures each individually explained as a test that depended on the silent fallback —
which is information worth having, not a reason to restore it.
OUTCOME: **PASSED via the second branch, and the list is the interesting part.** The first
full run after the code change gave **738 passed, 6 failed, 6 skipped** (586 s). With 10
new gate tests added to a 734-test baseline the expected total is 744, and 738 + 6 = 744 —
so every new test passed and exactly **six pre-existing tests** broke.

**All six were consuming the silent fallback, and all six are the same defect.** Five
report `Source Stepping failed at lambda=0.0` — failure at the *first* continuation rung,
with every source zeroed, which means structurally singular rather than merely hard:

| test | circuit | why there is no operating point |
|---|---|---|
| `test_elements.py::test_Idt_tran` | `Idt` | ideal integrator: output is the unbounded integral of a constant input |
| `test_elements.py::test_Idtmod_tran` | `Idtmod` | as above |
| `test_elements.py::test_Idtmod_modulo` | `Idtmod` | as above |
| `test_doctests.py::test_elements_module_doctests` | `Idtmod` docstring, 2 of 92 | as above |
| `test_stress_delayed_stiff_avalanche` | `Idtmod`, node 3 driven only by the integrator | no DC path to ground at all |
| `test_stress_charge_pump` | `C1`/`C2` open at DC, nodes reach ground only through diodes | numerically singular |

**Resolution: `uic=True` at each of the six sites, not a code change.** Zeros is the
physically correct initial state in every case — an integrator starts at zero, a charge
pump starts discharged — so no test's intent changed; what changed is that the intent is
now stated instead of arriving via a failure path. This is the outcome the gate anticipated
and it is recorded rather than smoothed over: **six tests in this suite have been
exercising a fallback rather than the behaviour they appear to test.**

**One further defect fell out of the fix**, which is why the second run was needed:
`_solve_coupled` **ignores `uic`** (one of the four inputs 0.1d found it dropping). That was
invisible while a failed DC silently produced zeros — ignoring `uic` and failing the DC gave
the same answer. With the fallback gone, `uic=True` raised on exactly the circuits that
need it. Fixed by honouring `self.par.uic` there too.

**Final tally: 744 passed, 6 skipped, 0 failed** (673.80 s) — 734 baseline + 10 new gates,
nothing else changed.

### Runtime, and a caution about measuring it on this box

The suite ran **497.69 s at baseline and 673.80 s after**, which reads as **+35%** and would
trip the standing 20% rule. **It is not a regression, and chasing it is the more useful
result.** What the investigation found:

1. Tests that touch **no transient code** moved just as much between runs —
   `test_SVCVS_laplace_d3_n1` went 19.52 s (baseline suite) -> 30.67 s (post-change suite)
   -> **21.59 s** when re-run on a quiet machine. Nothing in stage 1 can affect that test.
2. The one change that could plausibly cost time, `vabstol` 1e-6 -> 1e-12, is **free**:
   same circuit, **4094 steps both ways and an endpoint identical to every digit
   (delta exactly 0.000e+00)**. Newton converges quadratically, so the node x-tolerance was
   never the binding constraint. *This is the measurement decision 0.3a asked for and never
   had — the 10^6 loosening bought nothing, so reverting it costs nothing.*
3. Against a worktree at `8ce58c2`, the stress circuit gives **25095 steps and
   v(3) = 5.048556078 on both sides** — bit-identical.
4. cProfile on both: **4,880,532 calls (after) against 4,880,547 (before)** — fifteen calls
   apart out of 4.88 million, with the same hot spots in the same order. There is no extra
   per-step work.
5. Interleaved timing, three alternating repetitions: min **12.70 s vs 12.30 s (+3.3%)**,
   with a within-configuration spread of 24%.

**Conclusion: roughly 2-4% overhead, well inside the threshold**, and the honest reading of
the +35% is that a single-sample suite comparison is not a usable measurement on this
machine. **Recorded as a caution for every later stage's suite gate:** this box has shown up
to ±57% run-to-run variation on individual tests, so a runtime change of less than about
50% between two single suite runs means nothing. Attribute with interleaved repeats, step
counts, and call counts — not with one number against another.

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
OUTCOME: **PARTLY REFUTED. The win is real but smaller than recorded, and the gate's drift
criterion is unsatisfiable in principle for this change.**

*The measurement.* On the 139-unknown leapfrog, `benchmarks/transient_stage2.py`, min of 3:
**71.367 s -> 41.396 s, i.e. 1.72x**, step count identical (324 both). **Below the declared
>= 2x**, so the gate fails as written.

*The review's supporting figure does not reproduce.* It recorded "0.238 ms single-threaded
against 4.462 ms with 4 threads" for a 136x136 solve. Measured here at n=139, 300 reps:
**0.234 ms with the default thread count against 0.182 ms single-threaded** — a 1.29x
ratio, not 18x. The 4.462 ms figure is not reproducible on this box and should not be
cited again. What *is* real is that the saving is not confined to the solve: the solve is
~2% of runtime, so a solve-only effect could not produce 1.72x end-to-end at all. The cost
is thread-pool overhead spread across the many small array operations in assembly, which
is why the limit is applied around the whole transient rather than around `linearsolver`.

*The drift criterion cannot be met, and that is a defect in the gate, not the work.*
Measured drift single- vs multi-threaded: `max|dt| = 8.65e-17 s`, `max|dv| = 2.16e-15 V`.
**Verified this is the BLAS and not our code**: re-running with the *same* thread setting
reproduces the reference at exactly `0.00e+00` on both axes, so the solver is
deterministic; changing the thread count changes which LAPACK reduction order runs, and
floating-point addition is not associative. `0.00e+00` is the right criterion for 2b and
2c, which only reorganise our own arithmetic — it is **unachievable for any change that
alters which BLAS kernel executes**. The correct criterion here is: identical step count
(324) plus drift at BLAS rounding level, both met.

*Dependency decision.* **`threadpoolctl` is not installed on this machine.** Implemented
as an optional import, discovered the way `symengine` and `_ginac_ext` already are: when
present the transient limits BLAS to one thread itself; when absent nothing changes and
`blas_single_thread_available()` returns False. Documented, with the environment-variable
equivalent (`OMP_NUM_THREADS=1` etc.) which is how the 1.72x above was obtained. **Whether
to make it a hard dependency for 1.72x is a maintainer decision and is left open.** No
warning is emitted when >1 thread is detected — without `threadpoolctl` the thread count
cannot be reliably determined, and a warning on every run that might be wrong is worse
than the documentation.

**2b. Assembly: hoist the per-element probes, batch the scatter.** `circuit.py:1290-1402`.
`hasattr(self.toolkit, 'add_at')` fires once per element per stamp; `NumericToolkit` lacks
it, so `Toolkit.__getattr__` builds a formatted error string and **raises 255 000 times**
— 34% of total runtime is attribute-lookup machinery. Hoist the probes out of the loop and
scatter once via `np.bincount` (keep `np.add.at` for object/complex dtypes).
**Gate 2b:** `i`, `q`, `u`, `G`, `C` each compared against the current implementation on a
non-trivial circuit. Declared success: **max abs difference exactly 0.0**, and >= 1.8x on
assembly.
OUTCOME: **PASSED on both halves.** `benchmarks/transient_stage2.py --stamps` compares the
rewritten assembly against the pre-stage-2 code (retained verbatim in the gate script, not
in the library) over 8 random states:

| stamp | max abs diff | nonzeros |
|---|---|---|
| `i` | **0.000e+00** | 137 |
| `q` | **0.000e+00** | 135 |
| `u` | **0.000e+00** | 1 |
| `G` | **0.000e+00** | 701 |
| `C` | **0.000e+00** | 523 |

Timing: **1.86x end-to-end** with drift **exactly `0.00e+00`** and an identical step count,
against a gate asking 1.8x on assembly alone — assembly was ~96% of runtime, so 1.86x
end-to-end is more than 1.8x on assembly. Profile after: assembly is **~48%** of runtime
and attribute-lookup machinery is **under 5%**, from ~96% and ~34%.

**The nonzero column is not decoration.** The first version of this gate sampled `u` at
`t = 0`, where a sine source is zero — so `u` compared "exactly equal" against an all-zero
vector and proved nothing. It now samples at a quarter period. **A gate that can pass
against an empty result is not a gate**, and this one nearly shipped that way.

*Why bit-identical rather than merely close:* `np.add.at(lhs, idx, v)` and
`bincount(idx, weights=v)` both walk the index array in order and sum duplicates as they
are met, so the floating-point operations occur in the same sequence. The single `lhs +=`
at the end is exact because `lhs` is provably still zero at that point — the bincount path
is only reached when the toolkit has no `add_at`, and such a toolkit also returns `None`
from `batched_contributions`, so nothing has been accumulated yet. That reasoning was
checked against the code before relying on it, not assumed.

**2c. Stop re-assembling after Newton converges.** `transient.py:218-219` re-evaluates
`func(x)` at the point Newton just converged to, and `StandardNewton` discards the `(F, J)`
it already has. `:364` and `:409` each recompute `cir.q(x)` at that same `x`. Measured
3.17 assemblies per accepted step where ~2.17 are needed.
**Gate 2c:** instrument the per-step assembly count. Declared success: `G`/`C`/`i`/`u` at
~2.17 per step, `q` correspondingly reduced, waveform drift `0.00e+00`. Callers using
`provided_function` keep an opt-in exact-at-x recompute.
OUTCOME: **HALF PASSED, HALF REFUTED — and the refuted half is refuted on principle, not on
effort.**

| stamp | before | after | target |
|---|---|---|---|
| `G` | 3.06 | **3.06** | 2.17 |
| `C` | 3.06 | **3.06** | 2.17 |
| `i` | 3.06 | **3.06** | 2.17 |
| `u` | 3.06 | **3.06** | 2.17 |
| `q` | **5.08** | **3.06** | reduced |

*The `q` half passed.* The step controller and the history roll each recomputed `cir.q(x)`
at a state the assembly had just evaluated; those two calls are exactly the 5.08 - 3.06
gap, and they are now served from a value cached against the state that produced it (guard
is identity, then full equality on `x` — a stale charge vector would corrupt the LTE
silently, which is the failure class stage 1 exists to remove). Drift **exactly
`0.00e+00`**.

*The `G`/`C`/`i`/`u` half cannot be done without changing the answer.* The plan's premise
is that "`StandardNewton` discards the `(F, J)` it already has". It does — **but that
`(F, J)` is evaluated at `x`, while the solver returns `x_next`** (`nrsolver.py:32-71`:
`F, J = eval_FJ(x)`, then `x_next = x + xdiff`, then `return x_next`). The `(F, J)` in hand
is at the *previous iterate*, not at the converged point. Reusing it would feed the step
controller a Jacobian from the wrong state, changing the LTE and therefore the step
sequence — which stage 2 forbids outright.

Nor can it be recovered by restructuring: making the solver evaluate at `x_next` before
returning costs exactly the one assembly it would save. With `k` Newton iterations the
count is `k + 1` either way, and the measured 3.06 is `k ~ 2.06` plus the final evaluation.
**Reaching 2.17 would require Newton to converge in ~1.17 iterations — a convergence
improvement, not the removal of redundant work.** The review's "~2.17 are needed" is
therefore not a redundancy figure and should not be carried forward as one.

**Gate 2-final (compound).** Full suite `-m ""` at 734/6/0, and the end-to-end speedup
recorded against the stage-0 baseline. Declared success: >= 5x on the review's benchmark
circuit with drift `0.00e+00`. **If drift is not exactly zero, stop** — this stage is
defined as behaviour-preserving, and a nonzero drift means something else changed.
OUTCOME: **PASSED. Suite 744 passed, 6 skipped, 0 failed** (654 s) — unchanged from stage 1,
no test touched. End-to-end on the leapfrog, min of 5, 324 steps in every configuration:

| configuration | runtime | speedup | drift |
|---|---|---|---|
| stage-0 baseline | 71.367 s | — | — |
| **2b + 2c (as shipped)** | **29.518 s** | **2.42x** | **exactly 0.00e+00** |
| **2a + 2b + 2c** (BLAS held to 1 thread) | **13.738 s** | **5.19x** | 2.16e-15 V (BLAS only) |

**The >= 5x is met (5.19x), and the drift condition is met in the only sense available.**
The stage splits cleanly along the stop condition: everything pycircuit does itself — the
hoisted probes, the batched scatter, the charge cache — is **bit-identical, exactly
`0.00e+00`**, at 2.42x. The remaining 2.1x comes from the BLAS thread count, where exact
reproducibility is not achievable at all because a threaded LAPACK reduces in a different
order (established under gate 2a by reproducing the reference exactly at a *fixed* thread
count). Step count is 324 in every configuration, which is the invariant that actually
shows nothing about the method changed.

**The review's 10.5x does not reproduce.** `compound.py` recorded 8.80 s -> 0.84 s; the same
three changes measure **5.19x** here, and only 2.42x without the thread limit. The
individual components are also smaller than recorded (1.72x against a claimed 2.3-2.5x for
threads; the 4.462 ms solve figure not reproducible at all). **The direction and the
mechanism were right; the magnitude was overstated by about 2x.**

**Timing caution, restated because it is severe on this box.** The min-of-5 samples for the
shipped configuration were 29.52 / 37.45 / 43.95 / 44.87 / 55.01 s — an **86% spread**. The
single-threaded samples are tighter (13.74-19.85, 44%), which is itself evidence that
thread contention is the noise source. Every number above is a minimum of five; a single
sample here is worthless.

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
OUTCOME: **PASSED, 90.8x.** RC step response against the analytic solution, trapezoidal,
`uic=True`, `timestep` well above the controller's own step so `max_step` is not what
binds (233 steps against a `tend/timestep` of 10):

| `reltol` | steps before | err before | steps after | err after |
|---|---|---|---|---|
| 1e-3 | 24 | 1.3212e-01 | 36 | **3.1264e-03** |
| 1e-4 | 50 | 1.3212e-01 | 66 | **7.1688e-04** |
| 1e-5 | 103 | 1.3212e-01 | 127 | **1.5812e-04** |
| 1e-6 | 195 | 1.3212e-01 | 233 | **3.4419e-05** |

**The "before" column is the finding.** The error is identical to five digits at every
tolerance while the step count grows 8x — the solver did eight times the work for exactly
the same answer, because the one step the controller cannot check was also the largest one
in the run. After the ramp the error falls monotonically and by 90.8x.

**The default is now derived rather than proposed.** Sweeping the opening step:

| `firststep` | steps @1e-4 | reltol's effect (err 1e-3 / err 1e-6) |
|---|---|---|
| `timestep` | 50 | **1.0x** — reltol does nothing |
| `timestep*1e-1` | 58 | **1.0x** — still nothing |
| `timestep*1e-2` | 62 | 63.9x |
| **`timestep*1e-3`** | **66** | **90.8x** |
| `timestep*1e-4` | 69 | 90.0x |
| `timestep*1e-6` | 75 | 91.7x |

The effect saturates at 1e-3: smaller openings buy no accuracy and cost steps, and 1e-1 is
useless. **1e-3 is the knee**, which is why it is the default.

**On the Hairer estimate, which the plan asked to be measured rather than assumed:** the
ramp costs a *fixed* ~15 steps (see gate 3-2), spent climbing geometrically from
`timestep*1e-3` back to `max_step`. A Hairer-style estimate from `q'`/`q''` would open near
the right size and save most of those 15 steps — **worth it only for short runs**, since on
the review's benchmark the whole ramp is 3.1% and on a 300-tau run it is 4.2%. Recommend
**not** implementing it until a use case appears where startup dominates; recorded here so
the question is closed rather than open.

**Gate 3-2 (cost).** Declared success: step count rises no more than 20%. Measured on the
review's circuit: +5% for Trapezoidal.
OUTCOME: **PASSED on the circuit the gate names (+3.1%), FAILED on a short RC (+32%), and
the mechanism is understood.** On the review's benchmark circuit — the one the +5%
expectation came from — the leapfrog goes **324 -> 334 steps, +3.1%**, close to the
prediction.

The ramp is a **fixed startup cost of ~15 steps**, not a proportional one, so its relative
cost is set by run length:

| run length | open at max_step | ramped | growth |
|---|---|---|---|
| 10 tau | 50 | 66 | **+32.0%** |
| 30 tau | 85 | 100 | +17.6% |
| 100 tau | 155 | 170 | +9.7% |
| 300 tau | 355 | 370 | +4.2% |

**Recorded as a partial failure rather than rounded off**: a run of only a few time
constants pays over the declared 20%. That is the correct trade — those are exactly the
runs where the unevaluated opening step did the most damage — but it is a real cost and a
short-run user will see it. The 15 steps are the geometric climb from `timestep*1e-3` back
to `max_step` at the controller's 2x-per-step growth cap, which is `log2(1000) ~ 10` plus
the accept/reject overhead.

**Gate 3-3 (the order effect appears).** With the first step ramped, a 2nd-order method
must beat backward Euler on the same circuit and tolerance. Currently they are identical,
which is the clearest single symptom of the defect.
OUTCOME: **PASSED.** RC step response at `reltol=1e-5`:

| integrator | steps before | err before | steps after | err after |
|---|---|---|---|---|
| backward Euler | 789 | 1.3212e-01 | 949 | 1.5059e-03 |
| **trapezoidal** | 103 | 1.3212e-01 | 127 | **1.5812e-04** |
| Gear2 | 131 | 1.3212e-01 | 155 | 4.1607e-04 |

Before: all three agree **to five digits** while trapezoidal uses 7.7x fewer steps than
Euler — the integrator choice was invisible in the answer. After: trapezoidal is **9.5x
more accurate than Euler** *and* uses 7.5x fewer steps, which is the behaviour a
second-order method is supposed to have.

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

**4a-bis. NEW, from the source papers (2026-07-30): the controller has one threshold where
both papers have two.** `IntegralController` rejects whenever `err > 1.0`
(`stepcontroller.py:84`) — a single test that decides *both* "shrink the next step" and
"redo this one". Both papers separate those:

- **YWR Fig. 1** shrinks the next step when `ε_T > ε_spec`, but only **redoes** the current
  step when `ε_T ≥ F_redo · ε_spec`, with a distinct `F_RedoCut`.
- **Fang §4** describes the same structure in the method he compares against: "if the
  normalized local truncation error is greater than a given value (**we used 4.63**), the
  current time point will be recomputed with a smaller time step."

pycircuit folds `TRTOL = 7.0` into `etol`, which scales the single threshold but does not
*separate* the two decisions: the accept test and the predictor deliberately aim at the
same target (`stepcontroller.py:66-72`). So a step whose LTE is 1.01x over tolerance is
re-solved from scratch, where both papers would accept it and take a smaller next step.

**Hypothesis, to be measured and not assumed:** this is a candidate explanation for the
rejection counts in 0.3b — 57 rejections for `Gear2('ywr')` against 29 for `'classic'`. If
the single threshold is what drives rejection count, the estimator comparison in gate 4f is
confounded by it, which is a second reason 0.3b was right to defer the default. **Measure
this before 4f**, by adding an `F_redo` band and re-running the four-way comparison; if
rejection counts converge, the `'ywr'`-vs-`'classic'` gap was mostly the controller.

**4b. `MAX_REJECT` force-accept.** `transient.py:390` grows `dt` **10x** after accepting an
over-tolerance step. Variable-step BDF-2 is zero-stable only for ratio < 1+sqrt(2) =
2.414214; at 10x the parasitic root is 4.76. `Gear2(ywr)`, the shipped default, is the
only configuration measured reaching this path. Replace growth with an order drop, and
**warn** — an unbounded accepted truncation error must not be invisible.
**Gate 4b:** on the review's circuit, no step ratio exceeds 2.414, and any force-accept
emits a `RuntimeWarning` naming `t` and `h`.
OUTCOME: **PASSED on both halves.** Measured 2026-07-31 by
`benchmarks/transient_stage4.py --forceaccept`. The answer to the question the resume block
said to ask first is **yes — the path is still reached, by the shipped default among
others, and it was taking a ratio of exactly 10.0 when it did.**

**A CONCLUSION THIS GATE GOT WRONG FIRST, AND WHAT CORRECTED IT.** The initial sweep
covered reltol 1e-4/1e-5/1e-6 across 3 circuits and 5 integrator/formula combinations. In
all 9 of those runs `Gear2('ywr')` reaches the cap **zero** times and only
`Trapezoidal('ywr')` reaches it, 5-159 times per run — so this section was first written to
say the review's "`Gear2(ywr)` is the only configuration reaching the force-accept path"
had been *refuted*, and that stages 3 and 4c had moved the default off the path entirely.

**That was wrong, and the thing that refuted it was the warning added by 4b itself, on its
first full-suite run.** It fired in five places the sweep never looked, including
`test_step_count_and_error_respond_to_reltol` for **both** `Gear2('classic')` and
`Gear2('ywr')`. Re-measuring at **reltol 1e-3** — a looser tolerance than the sweep
covered — on the stiff RLC:

| config, stiff RLC @ 1e-3 | steps | force-accepts | worst ratio | ratios > 2.414 |
|---|---|---|---|---|
| **gear2-ywr (the default)** | 179 -> 181 | 1 -> 1 | **10.0000 -> 2.0000** | **1 -> 0** |
| gear2-classic | 193 -> 195 | 1 -> 1 | **10.0000 -> 2.0000** | **1 -> 0** |
| trap-ywr | 147 -> 159 | 41 -> 9 | 5.8987 -> 2.000 | 12 -> 0 |

So the review's §4.4 was **not** stale: the shipped default does reach the path and did
take the 10x. What is true is only the weaker claim that it is not the *only* configuration
and no longer the most frequent one — `Trapezoidal('ywr')` reaches it one to two orders of
magnitude more often. **The lesson is about the measurement, not the code:** a sweep of
three tight tolerances looked comprehensive and missed the regime where the default
misbehaves, because a loose tolerance is what lets the step grow into trouble. This is
recorded rather than quietly fixed because the corrected conclusion is the *less*
convenient one, and because the mechanism that caught it is the one 4b was adding anyway.

**The defect is confirmed live and is not an edge case.** On the stiff RLC at reltol 1e-5,
`Trapezoidal('ywr')` force-accepted **78 times in 873 accepted steps**, and the accepted-step
sequence contained **9 ratios above 2.414, the largest exactly 10.0** — all nine sitting
*immediately after* a force-accept. Nothing anywhere in any run produced a ratio above 2.0
by any other route, because the controller's own clamp is 2.0; the escape hatch was the
sole source.

Before and after, obtained by re-running the same script against the stashed source:

| circuit, reltol | steps | force-accepts | worst ratio | ratios > 2.414 |
|---|---|---|---|---|
| rc-vsin 1e-6 | 1419 -> 1423 | 5 -> 4 | 5.308 -> **2.000** | 3 -> **0** |
| stiff-rlc 1e-4 | 344 -> 393 | 36 -> 25 | 6.467 -> **2.000** | 8 -> **0** |
| stiff-rlc 1e-5 | 873 -> 923 | 78 -> 26 | 10.000 -> **2.000** | 9 -> **0** |
| stiff-rlc 1e-6 | 2301 -> 2322 | 20 -> 15 | 5.916 -> **2.000** | 9 -> **0** |
| rc-pulse 1e-4 | 1165 -> 1419 | 159 -> 121 | 3.878 -> **2.000** | 13 -> **0** |
| rc-pulse 1e-5 | 2484 -> 2757 | 113 -> 101 | 2.885 -> **2.000** | 2 -> **0** |
| rc-pulse 1e-6 | 6643 -> 6730 | 35 -> 35 | 2.945 -> **2.000** | 3 -> **0** |

(all `Trapezoidal('ywr')`; the Gear2 rows are in the 1e-3 table above.) Over the full
4-tolerance x 3-circuit x 5-configuration grid, **61 accepted step ratios outside the bound
become 0**, and the worst ratio anywhere in any run is now exactly the controller's own
clamp.

**The force-accept count falls as well, which was not predicted.** 78 -> 26, 159 -> 121,
36 -> 25. The order drop actually escapes the stalled regime; the 10x jump did not — it was
simply rejected again on the next step, so the old code paid four Newton solves to leave a
collapsed step size and then collapsed straight back. The step count rises 2-22% in
exchange, which is the honest cost of no longer taking illegal jumps.

**Blast radius: small, and confined to runs that were reaching the escape hatch.** Of the
45 configurations at the three tight tolerances, **36 are bit-identical in step count** —
every Euler run, every `trap-classic` run and every `Gear2` run. At 1e-3 the two Gear2 runs
that do reach the cap move by 2 steps each (179 -> 181, 193 -> 195). Full suite:
**761 passed, 6 skipped, 0 failed, 670.26 s** (`-m "" --timeout=400`), against a 755-passed
baseline plus the 6 tests added here — no test needed changing.

**The warning.** 26 `RuntimeWarning`s from the 26 force-accepts of the stiff 1e-5 run, each
naming `t` and `h`, and attributed to the caller's `solve()` line (`stacklevel=3`, not 2 —
the loop lives in `_solve`, so 2 blames `solve`'s own body and tells the caller nothing).
Asserted in `test_transient_repairs.py`. **It earned its place immediately**: see the
correction above, and note that it also located the path in `test_stress_charge_pump`,
`test_distortion_vs_transient` and gate 1-5's pulsed RC — three more places nobody had
looked.

**Where the end-to-end regression test points, and why.** At the shipped default
`Gear2('ywr')`, reltol 1e-3, not at `Trapezoidal('ywr')` — even though the latter is a much
louder case. 4d deletes the `'ywr'` trapezoidal branch as an alias of `'classic'`, so a
test built on it would quietly stop testing anything the day 4d lands. The default
configuration will still be there.

**4c. Backward Euler's variable-step bias.** `integrator.py:83` — `gn - gn_1` approximates
`((h1+h2)/2) q''`, not `h1 q''`. Rescale by `2*h_curr/(h_curr+h_last)`.
**Gate 4c:** est/true flat across step ratio 0.25..4. Declared success: within 5% of 1
across the range. Currently 0.62..2.47.
OUTCOME: **PASSED, and the plan's diagnosis was exactly right.** Measured by
`benchmarks/transient_stage4.py --ratios`, which drives the estimator directly with an
analytic charge so that `true = iq - q'(t_n)` is exactly computable:

| ratio `h_curr/h_last` | 0.25 | 0.5 | 1.0 | 2.0 | 4.0 |
|---|---|---|---|---|---|
| before | 2.5246 | 1.5089 | 1.0040 | 0.7522 | 0.6265 |
| **after** | **1.0098** | **1.0059** | **1.0040** | **1.0030** | **1.0025** |

The recorded range was 0.62..2.47 and it measures **0.6265..2.5246** — the plan's number
reproduces. Spread falls from 4.03x to **1.01x**. The rescaling by `2*h1/(h1+h2)` is
exactly the right correction, as proposed.

**The bias runs the wrong way twice**, which is worth stating because it explains the
symptom: on a *shrinking* step the error was overstated, so the controller shrank further
than it needed to; on a *growing* step it was understated, so the controller grew past what
the tolerance allowed. Both push toward the ratio extremes rather than damping them.

**4d. `Trapezoidal(lte_formula='ywr')` is order-inconsistent on non-uniform grids.**
`integrator.py:123` — the plain second difference leaves an O(h) term where the truncation
error is O(h^2); measured 112x too large at ratio 0.25, -436x at 4.0, scaling as 1/h. Its
correct divided-difference generalisation is algebraically identical to the `classic`
branch, so **delete the branch and keep `'ywr'` as an alias** rather than creating a
duplicate.

**CONFIRMED BY THE SOURCE PAPER, 2026-07-30.** YWR Table I gives TRAP as
`ε = −(1/12)(q̇_x + 0.5h ḟ_x)^{-1}(g_n − 2g_{n−1} + g_{n−2})h` — **a single `h` and an
unweighted second difference, i.e. a uniform-grid formula.** Its GEAR2 entry, by contrast,
carries `h1` and `h2` explicitly and weights the differences
(`h2·g_n − (h1+h2)·g_{n−1} + h1·g_{n−2}`). So the paper itself is non-uniform for GEAR2 and
uniform-only for TRAP, and pycircuit copied the TRAP entry verbatim **including its
uniform-grid assumption**. This is not a transcription error — it is a correct
transcription of a formula whose domain of validity was not stated in the table.

This also supplies the fix and a check on the plan's proposal: the paper's **eq (22) is the
general form** from which Table I's entries are derived by finite-difference approximation,
so the non-uniform TRAP formula can be derived rather than guessed. **Do that derivation
before deleting the branch** — the plan asserts the generalisation is algebraically
identical to `classic`, and eq (22) is how to confirm or refute that claim instead of
assuming it.
**Gate 4d:** est/true within 5% of the `classic` column across ratio 0.25..4.
OUTCOME: **BLOCKED ON 4g — and this reorders the stage.** The gate divides by the true
truncation error of the trapezoidal rule. On trapezoidal that quantity is **dominated by
the alternating mode 4g exists to remove**, so the ratio measures 4g and not 4d.

Demonstrated directly: holding the step size and ratio fixed at 1.0 and varying only the
number of steps taken beforehand, the *true* error at the final step is

| prefix length | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|
| euler | 1.4333e+04 | 1.4333e+04 | 1.4333e+04 | 1.4333e+04 | 1.4333e+04 |
| **trap-classic** | **-2.8012e+01** | **-5.7469e-01** | **-2.7817e+01** | **-7.7062e-01** | **-2.7620e+01** |
| gear2-classic | -5.6695e+01 | -5.6695e+01 | -5.6695e+01 | -5.6695e+01 | -5.6695e+01 |

Euler and Gear2 are identical to five digits regardless of history; **trapezoidal swings
48x on the parity of the step count alone.** That is the `(-1)^n` mode, and it is ~48x
larger than the smooth truncation error it is supposed to be measured against. Any
est/true figure for trapezoidal — including the plan's recorded "112x too large at ratio
0.25, -436x at 4.0" — is therefore a measurement of the contamination.

**Consequence: 4g must be done before 4d**, which the plan does not say; it lists 4d first
and treats the two as independent. 4d's gate should be re-run only once 4g-1 and 4g-2 pass.

> **CORRECTED 2026-07-31 by gate 4g-b0.** The parity table above is measured against the
> *propagated* companion-current error. Against the **local** truncation error --
> the quantity a step controller must estimate -- trapezoidal reads **-2.8395e+01 at
> every prefix length from 3 to 7**, i.e. no parity dependence at all. The 48x swing was
> the O(h^3) propagated mode passing near its own cancellation.
>
> **So 4d was never blocked on 4g.** It was blocked on a contaminated reference, and it
> measures cleanly the moment the reference is fixed. 4g(b) then resolved it for
> trapezoidal by construction anyway -- see gate 4g-b3 -- but that is a coincidence of
> sequencing, not a dependency. **The ordering claim in this OUTCOME is withdrawn.**

**Also observed, and NOT in the plan:** `gear2-classic` has its own step-ratio dependence —
2.7751 / 1.5504 / 0.9983 / 0.7766 / 0.6988 across ratios 0.25..4, a 3.97x spread, exact at
ratio 1 and biased off it, i.e. the same *shape* as the Euler defect 4c just fixed. It is
free of the parity contamination (see the table above), so unlike trapezoidal it **can** be
measured now. Recorded as a new item rather than folded silently into 4d; the `2*h1/(h1+h2)`
rescaling does not fix it (it leaves 1.110 / 1.118 at the extremes), so it needs its own
derivation.

**4e. `Gear2.check_order_drop` guards the wrong direction.** `integrator.py:158` fires on
*shrink*, which is unconditionally zero-stable; there is no guard on *growth*, which is
the unstable one — which is how 4b slipped through. Replace with an upper-ratio guard.
**Gate 4e:** a growth ratio above 2.414 triggers the order drop; a shrink does not.
OUTCOME: **First half PASSED. Second half REFUTED by measurement, and deliberately not
implemented.** The growth guard is added; the shrink guard is **kept**, against the plan's
"replace", and re-labelled for what it actually does.

*The growth guard, and the honest thing about it.* `Gear2.check_order_drop` now drops to
Euler above `ZERO_STABILITY_RATIO = 1 + sqrt(2)`. Measured over the same runs, **it fires
zero times** — because the controller's own clamp (`MAX_GROWTH_RATIO = 2.0`) keeps every
normal accepted step inside the bound, and 4b has stopped the one path that bypassed it.
The maximum ratio reaching the guard in any Gear2 run is exactly 2.000. **A backstop that
never fires in a healthy run is the point of it**: it converts "no accepted step ratio
exceeds the bound" from an accident of two clamps happening to agree into something the
integrator enforces for itself. It is asserted as a unit test rather than end to end,
because there is no longer any end-to-end path that reaches it.

*Why the shrink half was refuted.* The old guard is not idle: on the stiff RLC it fires
**3-6 times per run** (0.16-0.31% of Newton evaluations); on the RC + `VSin` it fires 0
times. Deleting it, as the plan specified, took `Gear2('ywr')` and `Gear2('classic')` from
**0 force-accepts to 1 each** on the stiff RLC at reltol 1e-5 — i.e. it converted a
controlled Euler step into a force-accepted over-tolerance one. Keeping it returns both to
0 and to bit-identical step counts (565 and 619, unchanged from before any of this work).

*What it is actually doing, which is why the plan's diagnosis was half right.* It is not a
stability guard and the comment claiming it was is deleted. It is a **stalled-estimate
heuristic**: a step only shrinks 10x below the last accepted one after several consecutive
rejections, and what rejects repeatedly is a 2nd-order estimate built on a third difference
of something that is not three times differentiable. Dropping to order 1 there is the same
medicine 4b now administers at the rejection cap, one retry earlier and without having to
accept an over-tolerance step to get there. **So 4e's defect was never that this branch
existed — it was that this branch was all there was, and was labelled as protecting a
stability property it has nothing to do with.**

**Reconsider if** the rejection cap ever becomes rejection-count-aware: `h_curr/h_last <
0.1` is a proxy for "we have rejected three times at this time point", and the thing it is
a proxy for is known exactly one level up in `transient.py`, where it would need no
threshold at all.

**Also landed here:** `MAX_GROWTH_RATIO` and `MIN_SHRINK_RATIO` in `stepcontroller.py`. The
2.0 and 0.2 were unexplained literals repeated across both controllers, and the force-accept
path's 10x was invisible precisely because there was no named bound for it to be compared
against. One bound, named once, imported by `transient.py`.

**A note for stage 9, checked while here.** `jaxtransient.calculate_next_dt` (`:351`)
clamps to `(0.1, 2.0)` — a **third** copy of the same two numbers, and its lower bound
disagrees with the CPU path's 0.2. The JAX path has no rejection cap and therefore does
**not** carry 4b's defect, so this is a consolidation item rather than a correctness one;
it belongs with 9(a)'s shared kernels.

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
OUTCOME: **PASSED in substance; the literal "0" was unachievable and the gate was
mis-specified.** `VSin`-driven RC, trapezoidal, five periods: **3 order-1 evaluations out
of 1596**, against the recorded 120 of 1236 — a 40x reduction.

The three remaining are evaluations **1, 2 and 3** — the genuine opening of the run, and
**zero** after it. The first step of a run has no history to difference, so it *must* run
at order 1; that is correct behaviour and the same fact stage 1 recorded in the
`_no_history` comment. **Declaring 0 was an error in the gate, not a target that was
missed**: no configuration of any correct integrator can achieve it. The measurable claim
is "no drive-synchronous order drops", and that is met exactly.
**Gate 4g-2:** est/true for Trapezoidal must not scale as 1/h. Declared success: ratio
bounded within 2x of 1 as h falls 1e-2 -> 1e-4. Currently -10 -> -2976.
OUTCOME: **FAILED — tier (a) is necessary but NOT sufficient, which is the measurement the
plan asked for before committing to (b).** est/true at a fixed step ratio of 1.0:

| | h=1e-8 | h=1e-9 | h=1e-10 |
|---|---|---|---|
| euler | 1.0371 | 1.0040 | 1.0004 |
| gear2-classic | 0.9829 | 0.9983 | 0.9998 |
| **trapezoidal** | **4.8652** | **65.1020** | **681.3979** |

Euler and Gear2 converge to 1 as the step shrinks, which is what a consistent estimator
does. Trapezoidal grows by **exactly 10x for every 10x reduction in h** — the 1/h law the
plan recorded as "-10 -> -2976", reproduced.

**Why tier (a) cannot fix it, stated so (b) is not attempted twice:** the harness that
produces the table above contains **no sources, no breakpoints and no order drops at all**
— it is the bare companion recursion driven by an analytic charge. The alternating mode is
*intrinsic* to the trapezoidal companion `iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}`, whose
homogeneous solution is `iq_n = -iq_{n-1}`, i.e. `(-1)^n`, undamped. Breakpoints **re-seed**
it — which is why tier (a) is still worth having, and why gate 4g-1 improved 40x — but
removing every breakpoint leaves the mode running from whatever seeds it at t=0.

**Tier (b) is therefore required, and its design is settled by the same algebra.**

> **CORRECTED 2026-07-31 by gate 4g-b0, before tier (b) was built.** The table above is
> `est / true_propagated`, and `true_propagated` is **O(h^3)** where `est` is O(h^2) --
> so the ratio grows as 1/h whatever the estimator does. It also changes sign at
> h=1e-11, which a diverging estimator cannot. **There was no 1/h divergence; the
> denominator was vanishing.** Against the local truncation error the same estimator
> reads 0.8067 / 0.6780 / 0.6678 -- a bounded ~33% underestimate.
>
> The conclusion "tier (b) is required" happens to survive, but **not for this reason**:
> the real defect is a 4x step-ratio spread and a 1.9x dependence on step history, not a
> divergence. See the 4g(b) section below for what was actually measured, and note that
> the reasoning quoted here -- the algebra of `d` -- was correct even though the
> measurement motivating it was not. The
trapezoidal relation gives `(iq_n + iq_{n-1})/2 = (q_n - q_{n-1})/h = d_n`, and the
alternating component **cancels exactly in that sum** because it flips sign each step. So
`d` is mode-free and the estimator should difference `d`, not `g`:

    LTE ~ -(1/6) (d_n - 2 d_{n-1} + d_{n-2})      [uniform grid; the divided-difference
                                                   generalisation is needed off it]

with `d_n = (q_n - q_{n-1})/h_curr`, `d_{n-1} = (q_{n-1} - q_{n-2})/h_last`,
`d_{n-2} = (q_{n-2} - q_{n-3})/h_last2`.

**The cost is the interface, exactly as the plan warned.** It needs a third past charge
(`get_required_history()` 1 -> 3 for Trapezoidal) and `h_last2` threaded through
`Integrator.compute_lte`, its three implementations, both `StepController` subclasses and
the transient's history bookkeeping. **Not attempted at the end of a long session**: stage
3 has already shown once that a hasty change to this code path breaks the two QUCS
reference tests, and an invasive interface change deserves a fresh start rather than a
tired one. The measurement that decides it is done and recorded; the change itself is next.
**Gate 4g-3:** re-run 0.2a. Declared success: the harness's integrator choice is either
confirmed or corrected, and the recorded 10x is either upheld or replaced.
OUTCOME: **Choice CONFIRMED. The 10x is REPLACED by 2.17x — and the difference was Euler's
broken estimator, not the trapezoidal rule's merit.** Same window, same tolerances:

| integrator | recorded | now | change |
|---|---|---|---|
| euler (the inherited default) | 2896 | **341** | **8.5x fewer steps** |
| trapezoidal | 288 | **157** | 1.8x |
| gear2 `'ywr'` | 347 | **163** | 2.1x |
| gear2 `'classic'` | 376 | **167** | 2.3x |
| **euler / trapezoidal** | **10.1x** | **2.17x** | |

Trapezoidal is still the fastest, so the harness's choice stands. But **the 10x that
justified it was mostly Euler's variable-step LTE bias** — gate 4c — which made Euler
demand roughly eight times more steps than its own accuracy required. Once the estimator is
correct the honest margin is 2.17x. This is a fourth review magnitude that shrank on
measurement, and unlike the others it shrank because a *defect elsewhere* was inflating it.

Re-running 0.2a itself: **D/A = E/B = 1.000 exactly**, because the probe's
"suppress breakpoints" monkeypatch is now a **no-op** — `Sin.next_event` already returns
`inf`. The ranking is 1.038 with and without, so 0.2a's original conclusion (breakpoints do
not decide the integrator choice) is upheld, now trivially. Note the gear2-vs-trapezoidal
gap has narrowed from 1.205x to **1.038x**: the three second-order methods are now within
4% of each other, which makes gate 4f's choice of default a much finer decision than the
pre-repair numbers suggested.

### 4g(b) — the mode-free estimator: derivation, scope and gates

**Written 2026-07-31, before any code.** Gate 4g-2 established that tier (a) is necessary
and not sufficient and that tier (b) is required; this section is the design it needs, and
its gates are declared here so they cannot be fitted afterwards.

#### The derivation, checked against the paper rather than guessed

The plan says to derive the non-uniform form from YWR eq (22) instead of assuming it is
"algebraically identical to `classic`". Read from a 200-dpi render of page 3:

    eps_T = SUM_i [ alpha_i (t_{n-i} - t_{n-p})^{k+1}
                    - beta_i (k+1)(t_{n-i} - t_{n-p})^k ]
            (qdot_x - beta_0 gdot_x)^-1  q*^{(k+1)}(zeta) / (k+1)!     (22)

and Table I's entries are stated to follow from "(22) **and finite difference
approximation**". **That phrase is the whole opening.** Eq (22) fixes the coefficient and
the quantity — `q'''` for a 2nd-order method — but the *choice of finite difference used to
approximate `q'''`* is free. The paper approximates it by a second difference of
`g = dq/dt`; nothing in eq (22) requires that, and for TRAP it is the choice that imports
the parasitic mode.

Evaluating (22) for TRAP (p=1, k=2, alpha = [1/h, -1/h], beta = [1/2, 1/2]) gives
`SUM = -h^2/2`, hence

    eps_T = -(h^2/6) (qdot_x + 0.5 h fdot_x)^-1 q'''(zeta)

and since `(qdot_x + 0.5 h fdot_x) = 0.5 h J`, the controller's own `J^-1` absorbs the
factor and the residual the integrator must return is **`Eg = -(h^2/6) q'''`**. Two
independent checks that this is right: it reproduces the paper's Table I TRAP entry exactly
once `(g_n - 2g_{n-1} + g_{n-2}) -> h^2 q'''` on a uniform grid, and it reproduces the
one-step companion error derived directly from `iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}` with
exact history, which is `-(h^2/6) q'''` as well.

**So the task is to estimate `q'''` without touching `iq`.** Let

    d_n = (q_n - q_{n-1}) / h_n

Two facts make `d` the right quantity. The trapezoidal relation gives
`(iq_n + iq_{n-1})/2 = d_n`, and the parasitic component flips sign every step, so it
**cancels exactly in that sum** — `d` is mode-free by construction. And Taylor-expanding
about the interval **midpoint** `m_n = (t_n + t_{n-1})/2` gives

    d_n = q'(m_n) + (h_n^2/24) q'''(m_n) + O(h^4)

i.e. `d` samples `q'` at midpoints, second-order accurate. A second divided difference of
`d` **over the midpoints** therefore estimates `q'''/2`:

    delta1 = (h_curr + h_last)/2          delta2 = (h_last + h_last2)/2
    DD2    = [ (d_n - d_{n-1})/delta1 - (d_{n-1} - d_{n-2})/delta2 ] / (delta1 + delta2)
    Eg     = -(h_curr^2/3) * DD2

On a uniform grid this collapses to `Eg = -(1/6)(d_n - 2 d_{n-1} + d_{n-2})`, which is
**exactly the form recorded under gate 4g-2** — the derivation was done independently and
lands on the recorded design, which is the check that it is the right one.

**A known, quantified imperfection, stated in advance so it is not discovered as a
surprise.** The `(h_n^2/24) q'''` term in `d_n` depends on `h_n`, which is *not* a smooth
function of `t`. On a uniform grid it is identical at all three points and cancels from
`DD2` exactly. Off it, hand-calculation at a sustained growth ratio of 2 gives a
contribution of `q'''/54` against a target of `q'''/2` — a **3.7% bias**. That is why gate
4g-b2 below declares 10% rather than the 5% used for 4c, and the number to check against is
3.7%, not zero.

#### Scope

**In:** `TrapezoidalIntegrator.compute_lte` rewritten on `d`; `get_required_history()`
1 -> 3; `h_last2` threaded through `Integrator.compute_lte`, its three implementations, both
`StepController.evaluate_step` subclasses, and the transient's history bookkeeping in
**both** transcriptions (`transient.py` `solve()` and `_solve_coupled()`); a fallback for
the one step of a run that has only two genuine `d` values.

**Out, each with the fact that would change it:**

* **Tier (c), `TRBDF2Integrator`.** Still the structurally cleanest answer — L-stable, no
  parasitic mode at all, embedded `q'''` estimator. **Reconsider if** 4g-b2 or 4g-b4 fails,
  or if trapezoidal ringing on stiff modes (review 4.6) turns out to matter for a real
  circuit; tier (b) fixes the *estimator*, not the method's ringing.
* **Seeding `iq` history from the operating point instead of zeros.** `transient.py:561`
  sets `_iqlast = zeros`, so the run starts with `iq_0 = 0` against a true `q'(t_0)` that
  is generally nonzero — a direct excitation of the parasitic mode. This is a real and
  separate defect, and tier (b) makes the *estimator* immune to it but does not remove the
  mode from the *solution*. **Reconsider if** the parity swing in the true error survives
  4g-b4; then the seed is the next thing to fix, and it is cheap.
  **RESOLVED, and the reconsider-if did not fire.** The parity swing was itself an artefact
  of the propagated reference (gate 4g-b0): the *local* error shows no parity dependence at
  all, and the propagated mode measures O(h^3), one order below the local truncation error.
  So the zero seed excites a mode that is real but an order smaller than the error the
  controller is already bounding. It remains wrong in principle and is worth fixing when
  something else touches that code, but it is not blocking anything and 4g(b) does not need
  it.
* Changing `Gear2`/`Euler` estimators. They measure correct already (4c, and 4g-2's own
  table), and touching them would make any change here unattributable.

#### Gates, declared before implementation

**Gate 4g-b0 (the reference is not fitted).** This is the gate that has to come first,
because 4g-b1 needs a reference and the obvious one is contaminated. The harness's `true`
is `iq_computed - q'(t_n)`, which contains error *propagated* from earlier steps as well as
the error committed at this one. A per-step estimator cannot and must not track the
propagated part — the parasitic mode is undamped and step-ratio independent, so shrinking
`h` does not reduce it, and a controller chasing it collapses the step for no accuracy.
The reference must therefore be the **local** truncation error, `-(h^2/6) q'''(t_n)`,
computed analytically.

**Changing the reference to make a number look better is exactly what rule 8 forbids**, so
this gate is what makes the change legitimate rather than convenient. Declared success:
for `euler`, `gear2-classic` and `gear2-ywr` — the estimators already known good — the
`est/true` figures under the analytic local reference agree with those under the existing
propagated reference to **within 5%** at h = 1e-8/1e-9/1e-10. If they agree, the new
reference is not fitted: it reproduces the old wherever the old was trustworthy and differs
only where the old was measuring accumulation. **If they disagree, tier (b) does not
proceed on this reference** and the disagreement is the finding.
OUTCOME: **PASSED, in the strongest available form — the disagreement is not "within 5%",
it is exactly zero.** Measured 2026-07-31:

| config | h=1e-8 | h=1e-9 | h=1e-10 | relative difference |
|---|---|---|---|---|
| euler | 1.03708 | 1.00396 | 1.00040 | **0.00e+00** |
| gear2-classic | 0.982878 | 0.998344 | 0.999840 | **0.00e+00** |
| gear2-ywr | 0.737159 | 0.748758 | 0.749880 | **0.00e+00** |

**And the reason is structural, not lucky.** `EulerIntegrator.compute_derivatives` and
`Gear2Integrator.compute_derivatives` are **pure functions of the charge history** — neither
reads `iq_last`. So for those methods "the error standing in the companion current" and "the
residual when exact history is substituted" are *the same expression*, and no choice of
reference can separate them. **Trapezoidal is the only method in the tree whose companion
current feeds back on itself** (`iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}`), and it is therefore
the only one for which the distinction exists at all. That is the cleanest possible answer
to "is the new reference fitted?": it cannot be, because it is identical to the old one
everywhere except the single method under investigation.

Cross-check of the derivation while here: the numerically-computed local error agrees with
the analytic `-(h^2/6) q'''` to ratio **0.96638 / 0.99669 / 0.999669** as h falls
1e-8 -> 1e-10, converging to 1 as a leading-order term must.

### THE PREMISE OF GATE 4g-2 IS REFUTED: there was never a 1/h divergence

**This is the most important result in 4g(b), and it arrived before a line of estimator
code was written.** Gate 4g-2 recorded trapezoidal `est/true` as
**4.8652 / 65.1020 / 681.3979** at h = 1e-8/1e-9/1e-10 and concluded "trapezoidal grows by
**exactly 10x for every 10x reduction in h** — the 1/h law". Taking the three quantities
separately instead of as a ratio:

| h | est | true_local | true_propagated |
|---|---|---|---|
| 1e-8 | -3.17072e+03 | -2.75313e+03 | -6.51709e+02 |
| 1e-9 | -3.74136e+01 | -2.83949e+01 | -5.74693e-01 |
| 1e-10 | -3.79295e-01 | -2.84797e-01 | -5.56643e-04 |
| 1e-11 | -3.84566e-03 | -2.86982e-03 | **+5.98188e-05** |

Per-decade scaling: `est` **98.6x**, `true_local` **99.7x** — both O(h^2), as they must be.
`true_propagated` scales **1032x**, i.e. **O(h^3)**, and at h=1e-11 it **changes sign**.

**So the estimator never diverged. The denominator vanished.** `est/true_prop` is
O(h^2)/O(h^3) = O(1/h) — the recorded law, produced entirely by dividing by a quantity one
order smaller that was on its way through a zero crossing. The sign flip at 1e-11 is the
tell: a diverging estimator does not change sign.

The same correction applies to the **parity table** recorded under gate 4d, which showed
trapezoidal's true error swinging 48x on the parity of the step count while Euler and Gear2
were flat. Re-measured against the local reference, trapezoidal's local truncation error is
**-2.8395e+01 at every prefix length from 3 to 7** — no swing at all. The 48x was the
propagated O(h^3) mode passing near its own cancellation, not the method's error.

**What this means for the plan, stated plainly.** Two of stage 4's recorded findings were
artefacts of dividing by a small number: the "1/h law" and the "48x parity swing". Neither
is evidence of anything. **4d was never blocked on 4g** — it was blocked on a contaminated
reference, and with the local reference it measures cleanly and immediately (see 4g-b3).

### What IS real, and why tier (b) was still built

Having refuted the stated premise, the honest next question was whether to build tier (b) at
all. The measurement that decided it: `est/true_local` across step ratio, with the estimator
fed the real companion currents, versus the same estimator fed exact derivatives as a
*diagnostic*:

| trap-classic | r=0.25 | r=0.5 | r=1 | r=2 | r=4 | spread |
|---|---|---|---|---|---|---|
| fed real companion currents | 3.2850 | 1.9737 | 1.3176 | 0.9887 | 0.8225 | **4.0x** |
| fed exact derivatives (diagnostic) | 0.9308 | 0.8861 | 0.8300 | 0.7733 | 0.7265 | **1.3x** |

The spread collapses from 4.0x to 1.3x when the mode is removed from what is differenced.
**So the contamination is real — it just shows up as a 4x step-ratio spread, not a 1/h
divergence.** That is the same magnitude, and the same shape, as the backward-Euler defect
gate 4c fixed (4.03x), which was judged worth fixing. Two further symptoms confirm it: the
same estimator on the same problem returns **1.3176 on one prefix and 0.6780 on another** —
a 1.9x swing from step history alone — and it does not converge to 1 as h falls, sitting at
**0.8067 / 0.6780 / 0.6678**, a persistent ~33% underestimate.

(The diagnostic column's *absolute* values are not meaningful — feeding exact derivatives is
the artefact recorded as trap 3, and it reproduces that trap exactly: Euler reads 0.5010.
Only the *spread* is being read from it.)

**Gate 4g-b1 (the 1/h law is gone).** `est/true` for Trapezoidal at fixed step ratio 1.0,
against the 4g-b0 reference. Recorded before: **4.8652 / 65.1020 / 681.3979** at
h = 1e-8/1e-9/1e-10, growing exactly 10x per decade. Declared success: bounded within
**2x of 1** across the same three h, and not monotone in h.
OUTCOME: **PASSED.** est/true at step ratio 1.0, under the 4g-b0 reference:

| | h=1e-8 | h=1e-9 | h=1e-10 | h=1e-11 |
|---|---|---|---|---|
| euler | 1.0371 | 1.0040 | 1.0004 | 1.0000 |
| gear2-classic | 0.9829 | 0.9983 | 0.9998 | 1.0054 |
| **trapezoidal, before** | 0.8067 | 0.6780 | 0.6678 | 0.6904 |
| **trapezoidal, after** | **0.9273** | **0.9933** | **0.9993** | **0.9952** |

The d-based estimator is **asymptotically exact**, converging to 1 like Euler and Gear2 do.
The g-based one it replaces sits at a persistent ~33% underestimate that does *not* improve
with h — which is the real defect, and is not the defect the gate was declared against.
Note the "before" row is itself prefix-dependent (it reads 1.3176 rather than 0.6780 on the
other prefix layout used earlier in this section); the "after" row is not.

**Gate 4g-b2 (it is correct off a uniform grid too).** `est/true` across step ratio
0.25..4, as for 4c and 4d. Declared success: within **10%** of 1 across the range — 10 and
not 5 because the midpoint-frame term above contributes a calculated 3.7% at ratio 2.
OUTCOME: **PASSED over the range the controller can reach; FAILED at r=4 as declared, by
12.9%.** Both halves are reported because the gate was declared 0.25..4:

| trap (d-based) | r=0.008 | r=0.05 | r=0.1 | r=0.25 | r=0.5 | r=1 | r=2 | r=2.414 | r=4 |
|---|---|---|---|---|---|---|---|---|---|
| est/true | 0.8782 | 0.8916 | 0.8986 | 0.9181 | 0.9468 | 0.9933 | 1.0577 | 1.0771 | 1.1292 |

Spread over 0.25..4 is **1.26x**, against **1540x** for the g-based `'ywr'` form and 4.0x
for the g-based `'classic'` form. Within 10% of 1 at r = 0.25, 0.5, 1, 2 and 2.414; outside
it at r=4 (**+12.9%**) and at deep shrinks below ~0.1 (**-12.2%** at r=0.008).

**Why the r=4 failure is reported rather than the gate widened.** `MAX_GROWTH_RATIO` is 2.0
and 4b removed the only path that bypassed it, so **no accepted step ratio above 2.0 can
occur**; r=4 is unreachable by construction and r=2.414 is the outermost value the bound
itself permits. The deep-shrink end *is* reachable — three consecutive rejections give
0.2^3 = 0.008 — and -12.2% there is a real, if mild, over-optimism. Both are the predicted
midpoint-frame term: `d_n` carries `(h_n^2/24) q'''`, which cancels exactly on a uniform
grid and not otherwise. The prediction written before the measurement was 3.7% at r=2; the
measurement is **5.8%**, same term, same order, larger coefficient.

> **RETIRED 2026-07-31 by gate 4i-3.** The midpoint-frame term is gone: 4i differences the
> charge itself, which carries no `h`-dependent term at all, and trapezoidal now measures
> **within 1.3% across the whole reachable range** (spread 1.0084x against the 1.26x this
> gate recorded). The partial failure above is no longer a live limitation, and the
> measurement of it is what motivated 4i to supersede the estimator this gate was written
> for — one commit later.

**A finding about the SHIPPED DEFAULT that this sweep turned up and 4g does not fix.**
Extending the ratio sweep to the deep-shrink end, which no previous gate did:

| | r=0.008 | r=0.05 | r=0.1 | r=0.25 | r=1 | r=4 |
|---|---|---|---|---|---|---|
| euler | 1.0020 | 1.0021 | 1.0022 | 1.0025 | 1.0040 | 1.0097 |
| trapezoidal (after 4g-b) | 0.8782 | 0.8916 | 0.8986 | 0.9181 | 0.9933 | 1.1292 |
| **gear2-classic** | **83.06** | **13.32** | **6.71** | **2.79** | 0.9983 | 0.6952 |

`Gear2` — the default integrator — **overestimates its own truncation error by 83x at a step
ratio of 0.008**, which is reachable after three consecutive rejections. An estimator that
inflates by 83x when the step has just collapsed will demand a further collapse, which is a
plausible mechanism for the rejection cascades 4b's escape hatch exists to break. This is
the "`gear2-classic` has its own step-ratio dependence" item already recorded under 4d, but
recorded there as a **3.97x** spread over 0.25..4; over the *reachable* range it is **119x**.
Not fixed here — touching Gear2 in the same change would make 4g unattributable — but it is
now the largest known estimator defect in the tree.

**Gate 4g-b3 (4d becomes measurable, and is then answered).** With 4g-b1 passing, re-run
gate 4d: `Trapezoidal('ywr')` est/true within 5% of the `classic` column across ratio
0.25..4. This gate exists to be *run*, not to be assumed — 4d's own OUTCOME records that it
was blocked precisely because the ratio measured 4g. Note this may now compare *three*
formulas rather than two.
OUTCOME: **PASSED, trivially and by construction — and 4d is thereby resolved for
trapezoidal rather than merely unblocked.** The d-based branch does not read
`lte_formula` at all, so `Trapezoidal('ywr')` and `Trapezoidal('classic')` now return
**identical** values at every step ratio: 0.9029 / 0.9412 / 0.9933 / 1.0622 / 1.1395. The
difference between them is 0 by construction, not within 5%.

The `'ywr'` and `'classic'` branches survive only on the single fallback step at the start
of a run, where three past charges do not yet exist. That is visible end to end: in gate
4g-b4's sweep `trap-classic` and `trap-ywr` are bit-identical on every circuit except the
stiff RLC, where they differ by 2-3 steps out of 150-460 — the one step of divergence,
propagated.

**What this means for 4d as a plan item.** Its stated fix was "delete the `'ywr'` branch and
keep it as an alias of `'classic'`". For **trapezoidal** that is now nearly true already and
the remaining work is cosmetic: make the fallback use one formula instead of two, and the
distinction disappears entirely. 4d's *other* half — the `gear2-classic` step-ratio
dependence recorded under its OUTCOME — is untouched and, per 4g-b2 above, is much larger
than recorded.

**Gate 4g-b4 (end to end, on the circuits that were suffering).** `Trapezoidal('ywr')` and
`Trapezoidal('classic')` on the gate-4b sweep. Declared success: rejection counts fall by
at least **2x** on the stiff RLC at reltol 1e-5 (recorded: 757 rejections and 26
force-accepts for `trap-ywr`, 18 and 0 for `trap-classic`), **no** accepted step ratio
leaves the zero-stability bound, and no configuration regresses in step count by more than
20%.
OUTCOME: **PASSED on all three clauses, and the rejection fall is 33x, not the 2x
declared.** `Trapezoidal('ywr')`, before -> after:

| circuit, reltol | steps | rejections | force-accepts |
|---|---|---|---|
| rc-vsin 1e-4 | 212 -> 212 | 3 -> **0** | 0 -> 0 |
| rc-vsin 1e-5 | 298 -> 223 | 108 -> **1** | 6 -> **0** |
| stiff-rlc 1e-3 | 159 -> 155 | 73 -> **16** | 9 -> **1** |
| stiff-rlc 1e-4 | 393 -> 262 | 209 -> **19** | 25 -> **1** |
| stiff-rlc 1e-5 | 923 -> 464 | 757 -> **23** | 26 -> **1** |
| rc-pulse 1e-3 | 1419 -> 712 | 875 -> **40** | 121 -> **0** |
| rc-pulse 1e-4 | 1419 -> 1178 | 875 -> **43** | 121 -> **0** |
| rc-pulse 1e-5 | 2757 -> 1942 | 1868 -> **76** | 101 -> **0** |

**33x fewer rejections on the stiff RLC at 1e-5 against a declared 2x, and the step count
HALVES** (923 -> 464) rather than regressing — the run was spending its rejections on an
estimator that could not settle. `worst accepted ratio 2.0000, ratios outside the bound 0`
across the whole sweep, so 4b's invariant survives. `trap-classic` moves by at most +6.9%
(rc-pulse 1e-3, 666 -> 712), well inside the declared 20%.

**Blast radius: exactly zero outside trapezoidal.** Every `euler`, `gear2-classic` and
`gear2-ywr` row in the 45-configuration sweep is **identical** to the gate-4b run, because
neither reads `h_last2`.

**Gate 4g-b5 (the two QUCS reference tests still hold).** Named explicitly because stage 3
broke exactly these with a hasty change to this code path, and because they are the only
tests comparing a transient against an external simulator.
OUTCOME: **PASSED.** `test_transient_RC`, `test_transient_RLC` and
`test_transient_nonlinear_C` all pass (29 passed in the selected subset). Worth noting why
they were never at much risk this time: all three run `fixed_timestep=True` or a step the
controller does not bind, so they exercise `compute_derivatives` — which 4g(b) does not
touch — rather than `compute_lte`. Stage 3's change touched the former; this one does not.

**Gate 4g-b6.** Full suite `-m ""`, runtime recorded; doc build verified by content.
OUTCOME: **PASSED. 770 passed, 6 skipped, 0 failed, 722.89 s** (`-m "" --timeout=400`),
against 761 before 4g(b) plus the 9 tests added here — **no existing test needed
changing**, on a change that moves trapezoidal step counts by up to 2x. Runtime 670 -> 723 s
is +7.9%, inside the standing 20% rule and inside this box's own run-to-run variation
(trap 2).

Doc build: **build succeeded, 2 warnings, 0 ERROR**, verified per rule 3 rather than from
the exit code — no `exec-rst` block fell back to rendering its own source (searched for
each block's own function definitions: 0 occurrences) and every table cell is computed.

**One defect the doc build caught, worth recording because it is trap 3 for the third
time.** The first version of the new page's comparison block seeded the g-history with an
exact derivative rather than recursing it, and produced a flat **0.8354 / 0.8334 /
0.8333** — i.e. 5/6, the artefact — instead of the real 1.0920 / 1.3130 / 1.3314. It was
caught only because 5/6 was recognisable from the `integrator.py` comment. A generated
table is not automatically a trustworthy one; it is only as good as the history it builds.

---

---

# 4i — the Gear2 estimator's step-ratio bias

**NEW ITEM, raised by 4g(b)'s sweep on 2026-07-31 and written before any code.** The plan
carried this as "the new `gear2-classic` ratio dependence", recorded under 4d's OUTCOME as
a **3.97x spread over ratio 0.25..4**. Extending the sweep to the *reachable* shrink end
shows it is far larger than that, and it is on the **shipped default integrator**:

| `gear2-classic` est/true | r=0.008 | r=0.05 | r=0.1 | r=0.25 | r=1 | r=2 | r=4 |
|---|---|---|---|---|---|---|---|
| measured | **83.06** | 13.32 | 6.71 | 2.79 | 0.998 | 0.775 | 0.695 |

Ratio 0.008 is reachable: three consecutive rejections shrink by `0.2^3`. So on the step
where the controller has just collapsed the step size, the estimator tells it the error is
**83x worse than it is**, and it collapses further.

## The mechanism, predicted before it was measured

Gear2's exact local truncation error is obtained from the interpolation-error form: the
quadratic through `t_n, t_{n-1}, t_{n-2}` differentiated at `t_n` leaves

    iq_n - q'(t_n) = -(1/6) h1 (h1 + h2) q'''(t_n)

(equal steps: `-(1/3) h^2 q'''`, the textbook BDF-2 constant). So **every companion current
in the history carries an error of that shape**, `e_k = -(1/6) h_k (h_k + h_{k-1}) q'''`,
and the estimator — which takes a second divided difference of `g` at the nodes —
differences those errors along with the signal. Computing `DD2(e)` by hand against a target
of `q'''/2`:

| h1/h2 | 0.008 | 0.05 | 0.1 | 0.25 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|---|
| **predicted** est/true | **83.34** | 13.365 | 6.727 | 2.800 | 1.0000 | 0.778 | 0.700 |
| **measured** est/true | **83.06** | 13.32 | 6.71 | 2.79 | 0.998 | 0.775 | 0.695 |

**Agreement to 0.3% at every ratio.** The defect is explained, not merely observed, and the
prediction is what makes the fix's target unambiguous. It also explains the shape: `DD2(e)`
is exactly zero when `h1 = h2 = h3`, which is why the estimator measures asymptotically
exact on a uniform grid (1.000282 against 2/9, recorded under 0.3b) and is wrong everywhere
else. **That measurement was never wrong; it was taken at the one ratio where the defect
vanishes.**

This is the same family as Euler's (4c, 4.03x, fixed by a rescale) and trapezoidal's
(4g(b), 4.0x, fixed by differencing something else) — the third instance of "the estimator
differences a quantity that carries the method's own error".

## The fix, derived rather than guessed

The obstacle was stated in `integrator.py` itself: *"Gear-2 keeps just two past charges
(`get_required_history() == 2`), so a third divided difference of q is not available at
all."* **4g(b) lifted exactly that constraint** — it added `h_last2` and a third past charge
for trapezoidal, and the plumbing is already threaded through `compute_lte`, both
`StepController` subclasses and both of the transient's history transcriptions. Gear2 needs
no new interface, only `get_required_history() -> 3` and a different formula.

With four charges, `q'''` is available directly and **exactly**: the third divided
difference of the charge is, by the mean-value form of the interpolation error,

    q[t_n, t_{n-1}, t_{n-2}, t_{n-3}] = q'''(zeta) / 6      for some zeta in the span

with **no coefficient that depends on the step ratios** — which is precisely what both
previous formulations lacked. Substituting:

    Gear2:        Eg = -(1/6) h1 (h1+h2) q'''  =  -h1 (h1+h2) * DD3
    Trapezoidal:  Eg = -(1/6) h1^2       q'''  =  -h1^2       * DD3

so **one shared `q'''` estimator serves both**, each scaled by its own error constant, and
the constants come straight from YWR eq (22) rather than from a difference formula.

## Scope: this also supersedes the estimator 4g(b) landed one commit ago

**Stated plainly because it widens what was asked for.** Prototyped against the exact local
truncation error, over the ratios the controller can actually produce:

| | current | with shared DD3 |
|---|---|---|
| `gear2` spread over 0.008..2.414 | **119x** | **1.004x** |
| `trapezoidal` spread over 0.008..2.414 | 1.26x | **1.008x** |

The DD3 form is better for trapezoidal too, and by enough to matter: it removes the ±12%
residual at the extremes that gate 4g-b2 had to report as a partial failure. Leaving
trapezoidal on the midpoint-`d` form would mean keeping a measurably worse formula *and* a
second `q'''` estimator in the same file, for no reason other than that it shipped first.
So both are moved onto the shared helper, **with separate gates** so the two changes stay
attributable.

The reason DD3 wins is worth stating: the midpoint form differences
`d_k = q'(m_k) + (h_k^2/24) q'''(m_k)`, whose second term depends on `h_k` and therefore
does not cancel off a uniform grid — the same *kind* of contamination as the `g`-based
forms, one order smaller. DD3 differences the charge itself, which carries no method error
at all.

**Out of scope, with the fact that would change it:** the `lte_formula` parameter. Once both
integrators use DD3, `'ywr'` and `'classic'` select nothing except on the single fallback
step, for Gear2 as well as trapezoidal. **Reconsider — and this is now a live question for
gate 4f** — whether the parameter should be removed outright; it is currently a knob that
advertises a choice it no longer makes. Not done here because deleting a public parameter is
a maintainer's call, not an implementer's.

## Gates, declared before implementation

**Gate 4i-1 (Gear2 is flat across the reachable ratio range).** `est/true` against the
exact local truncation error, at ratios 0.008 .. 2.414 — 0.008 because three consecutive
rejections reach it, 2.414 because `ZERO_STABILITY_RATIO` is the largest ratio the
integrator will run at second order. Declared success: **within 5% of 1 at every ratio**,
against the measured 83.06 .. 0.695.
OUTCOME: **PASSED, with 7x more margin than declared.** Measured 2026-07-31 by
`benchmarks/transient_stage4.py --ratios` extended to the reachable range:

| est/true | r=0.008 | r=0.05 | r=0.1 | r=0.25 | r=0.5 | r=1 | r=2 | r=2.414 |
|---|---|---|---|---|---|---|---|---|
| **gear2, before** | **83.06** | 13.32 | 6.71 | 2.79 | 1.55 | 0.998 | 0.775 | 0.745 |
| **gear2, after** | **0.9966** | 0.9966 | 0.9965 | 0.9963 | 0.9958 | 0.9950 | 0.9933 | 0.9925 |

Worst deviation **0.7%** against the declared 5%. Spread over the range falls from
**119x to 1.0041x**. `gear2-classic` and `gear2-ywr` are now identical, exactly as
trapezoidal's two became under 4g-b3 — the shared helper does not read `lte_formula`.

**Gate 4i-2 (it is consistent).** `est/true` for Gear2 at ratio 1 as h falls
1e-8 -> 1e-11. Declared success: converges toward 1 and is within 2% at h = 1e-10.
OUTCOME: **PASSED.** 0.9442 / 0.9950 / **0.9995** / 1.0031 at h = 1e-8/1e-9/1e-10/1e-11 —
within **0.05%** at 1e-10, against the declared 2%.

**Gate 4i-3 (trapezoidal does not regress).** The same sweep for Trapezoidal. Declared
success: spread over 0.008..2.414 **no worse than the 1.26x it has today**, and within 5%
of 1 at every ratio — i.e. strictly better than gate 4g-b2, which had to report +12.9% at
r=4 and -12.2% at r=0.008.
OUTCOME: **PASSED, and it retires gate 4g-b2's partial failure.** Trapezoidal now measures
0.9867 / 0.9949 / 0.9949 / 0.9946 / 0.9942 / 0.9933 / 0.9916 / 0.9909 across
0.008..2.414 — spread **1.0084x** against the 1.26x it had one commit ago, worst deviation
**1.3%** against 4g-b2's -12.2%. The `(h_k^2/24) q'''` term that 4g(b) had to declare as a
known residual is gone, because the charge carries no such term.

**Gate 4i-4 (end to end).** The gate-4b sweep, all three circuits, four tolerances, five
configurations. Declared success: **no accepted step ratio outside the zero-stability
bound**; **no configuration regresses in step count by more than 20%**; and Gear2's
rejection and force-accept counts do not rise. The last clause is the one that matters —
an estimator that no longer inflates 83x on a collapsed step should reject *less*, and if it
rejects more the fix has done something other than what this section claims.
OUTCOME: **TWO CLAUSES PASSED, THE THIRD FAILED — and the clause was mis-specified, which
is a different thing from the fix being wrong.** Recorded in full because the failing clause
was declared as "the one that matters".

*What passed.* `worst accepted step ratio 2.0000, ratios outside the bound 0` across the
whole sweep. Step counts rise **2-9%** for Gear2 (largest: rc-pulse at reltol 1e-3,
808 -> 880), well inside the declared 20%. Euler is bit-identical everywhere; trapezoidal
moves by at most 1 step.

*What failed.* Gear2's rejection counts **rose** on 5 of 9 configurations — rc-pulse 1e-3
36 -> 53, 1e-4 42 -> 60, 1e-5 43 -> 72; stiff-rlc 1e-3 19 -> 20, 1e-4 15 -> 17.

*Why, and why the clause was wrong rather than the code.* The clause assumed the defect was
one-directional — that the estimator only ever *over*-reported. It does not. From the table
under 4i-1, the old estimator over-reports on **shrink** (83x at 0.008) and **under-reports
on growth** (0.775 at r=2, 0.695 at r=4). A controller spends most of a run growing toward
`max_step`, so it was being told the error was ~22% smaller than it is, and **it rarely
rejected because it was under-reporting, not because its steps were good.** An honest
estimator runs closer to the true tolerance boundary and therefore crosses it more often.
More rejections is what correct error control looks like here.

*The compensating measurement, taken because the assertion above is not self-evidently
true.* Global error on the stiff RLC against the analytic solution, before -> after:

| config | reltol | steps | rejections | global error |
|---|---|---|---|---|
| **gear2-ywr** (the default) | 1e-4 | 63 -> 69 | 7 -> 7 | 2.4325e-03 -> **2.0090e-03**, 1.21x better |
| **gear2-ywr** | 1e-5 | 102 -> 112 | 7 -> 7 | 5.4326e-04 -> **4.4443e-04**, 1.22x better |
| **gear2-ywr** | 1e-6 | 185 -> 204 | 7 -> 7 | 1.1917e-04 -> **9.7563e-05**, 1.22x better |
| gear2-classic | 1e-5 | 111 -> 112 | 7 -> 7 | 4.4660e-04 -> 4.4443e-04, 1.00x |
| euler | all | unchanged | unchanged | **bit-identical** |

The shipped default gains **1.22x accuracy for ~10% more steps**. `gear2-classic` is
neutral, and the reason the two differ is that they now *agree*: YWR's GEAR2 residual
estimates `(1/4) h^2 q'''` against a true `(1/3)`, so `'ywr'` was optimistic by a further
3/4 on top of the ratio bias and had more to gain. Both rows now read identically.

**STATED PLAINLY, BECAUSE IT IS THE LESS FLATTERING READING: this is not a speedup.**
Interpolating the `gear2-ywr` rows, reaching the *old* error of 1.1917e-04 now costs about
195 steps against the old 185 — so at matched accuracy the new estimator is marginally more
expensive on this circuit. What it buys is that **`reltol` now means what it says**, which
is the same value proposition as 4c and stage 3, plus robustness in the regime the 83x
lives in: a step that has just collapsed after three rejections is no longer told its error
is 83x worse than it is. The sweep circuits reach ratios below 0.1 only 3-6 times per run,
so the sweep understates that benefit rather than demonstrating it.

**Gate 4i-5 (the fallback still cannot go blind).** With `h_last2 is None` — the one step of
a run that has fewer than four charges — both integrators must return a non-zero, finite
estimate rather than zeros. Zeros would make that step unchecked, which is the defect stage
3 removed from the first step.
OUTCOME: **PASSED.** All four configurations return finite, non-zero estimates on the
fallback path: trap `-2.356777e+01` (both formulas), gear2-classic `-3.770212e+01`,
gear2-ywr `-2.827659e+01`. Note the two Gear2 fallbacks still differ, which is the only
place `lte_formula` now has any effect at all for either second-order method.

**Gate 4i-6.** Full suite `-m ""`, runtime recorded; doc build verified by content, not by
exit code.
OUTCOME: **PASSED. 778 passed, 6 skipped, 0 failed, 722.88 s** (`-m "" --timeout=400`),
against 770 before 4i plus the 8 tests added here. **Exactly one existing test needed
changing** — `get_required_history()` for Gear2, 2 -> 3 — which is the churn 4i was
expected to cause and is the very constraint `integrator.py` used to record as the reason
the g-based form had to be used. Runtime is unchanged from the 4g(b) run to within a
second (722.89 -> 722.88 s), which is coincidence rather than a measurement: this box
varies far more than that between runs (trap 2).

Doc build: **build succeeded, 2 warnings, 0 ERROR**, verified per rule 3 — no `exec-rst`
block fell back to rendering its own source, and the regenerated tables reproduce the
benchmark numbers exactly (`GEAR2, g` at ratio 0.008 reads **83.0625**, matching both the
sweep and the hand calculation).

**A defect in the generated documentation, caught before commit and recorded because it
is trap 3 for the third time.** The new comparison table fed **Gear-2 a trapezoidal
companion history** — the two recursions are different functions — and reported **41.4**
where the truth is 83.06. It was caught only by noticing the number disagreed with the
benchmark. Each method now builds its own companion history: trapezoidal by recursion
over a prefix, Gear-2 evaluated directly at each past node, since its companion is a pure
function of three charges.

**A second doc-build gotcha, new and worth recording.** The first attempt to regenerate
this page produced *stale* tables and still reported "build succeeded": the `.rst` edit
had failed silently (an `index()` anchor matched an earlier prose mention of the same
heading inside a code comment), so sphinx saw an unchanged source and skipped the page.
**A doc build that reports success proves nothing about a page whose source it did not
read.** Forcing a full rebuild with `-E` re-runs every `exec-rst` block in the project and
exceeded ten minutes, so the workable check is to confirm the source actually changed and
then read the regenerated values.

---

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

# STAGE 12 — implement Fang's method for real

Added 2026-07-31, after decision **D1** chose "keep and fix" for `_solve_coupled`. That
decision leaves a module whose comment cites a paper it does not implement; this stage is
the option that resolves it by making the citation true. **The alternative remains open and
is cheaper** — drop the citation, rename the flag to describe the rejection loop it
actually is, and fix the four ignored inputs. Do not start this stage without deciding
which of the two is wanted, because doing the cheap one first is wasted if this lands.

**Source:** G. Peter Fang, *"A new time-stepping method for circuit simulation"*, DAC 2013 —
`/home/andreas/pycircuit_agy/papers/2463209.2488904.pdf`. Read it from rendered pages; the
formulas do not survive text extraction. What 0.1d established, so it is not re-derived:
the code implements **neither** §3.1 nor §3.4, and `analytical_eh` is a vestige of the
`E_h` gradient that `doc/src/circuit/time_stepping.rst` once documented but that was never
written.

**What the method actually is.** The step size `h_m` becomes an **unknown**, solved together
with the circuit equations as `N+1` nonlinear equations (eq 11), through the bordered system

    [ J    p ] [ Δv ]   [ −f_ckt ]
    [ q^T  d ] [ Δh ] = [ −f_lte ]                                    (eq 12)

with `p = ∂f_ckt/∂h_m`, `q^T = ∂f_lte/∂v_m`, `d = ∂f_lte/∂h_m`. §3.2 gives two ways to
solve it without a full `(N+1)` factorisation: partial LU keeping the last row undetermined,
or the reduction to an `N`-system `(J − p q^T / d) Δv = −f_ckt + f_lte p / d` (eq 13) with
`Δh` recovered from eq (14) — a rank-one update, so the existing factors can be reused.
§3.4's approximate Newton avoids re-solving altogether: predict `h` from Gear's formula
(eq 17) and **correct** the solution already computed,
`Δv^{k+1} = Δv^{k+1/2} − J^{-1} p (h^{k+1} − h^k)` (eq 18).

**ENTRY CONDITION — a measurement, not a date. Do not start this stage until it is met.**
Fang's benefit is the elimination of *rejected* steps. Stage 4 already attacks rejections
from four directions (4a's unstable PI gains, 4b's 10x force-accept growth, 4e's inverted
order-drop guard) and **4a-bis is nearly the same benefit for a fraction of the cost**:
both source papers use a *two-threshold* controller that only redoes a step when
`ε >= F_redo * ε_spec` — Fang's own comparison method used **4.63** — where pycircuit
re-solves whenever `err > 1.0`. A step 1.01x over tolerance is currently recomputed from
scratch.

So: **run gate 4f, and read the rejection counts.** If stage 4 leaves rejections near zero,
this stage has almost nothing left to buy and D1's cheap resolution is the honest one. If
rejections are still material — say above 10% of accepted steps on the stress circuits —
then the bordered system is buying something a threshold cannot, and this stage is
justified. **That number does not exist yet; it is what gate 4f produces.**

**Dependencies, and they are real:**

- **Stage 4 must be complete.** The method's entire premise is that the LTE is a quantity
  worth solving *for*. Making `h` an unknown driven by an estimator that is 681x wrong at
  small `h` (4g-2) would be worse than the rejection loop it replaces, not better.
- **Stage 7b** gives the `factor`/`solve` split that eq (13)'s rank-one update and §3.2's
  factor reuse both need. Without it the "very little overhead" the paper claims is not
  available and the method is a pessimisation.
- **Item 2+.2's seam** — the solver returning `(F, J)` at the converged point — was deferred
  to stage 7a for the reduced/full reconciliation. It is a prerequisite here too.

**Work.** (a) `p`, `q^T`, `d`: `p = ∂f_ckt/∂h` is analytic given the integration formula
(`f_ckt` contains `q/h` terms), and `d = ∂f_lte/∂h` follows from the LTE formula's own
`h`-dependence — this is what `analytical_eh` was named for. (b) The reduced `N`-system of
eq (13) rather than a literal `(N+1)` factorisation. (c) §3.3's **two-sided** band
`γ_min τ ≤ ε ≤ γ_max τ` (eq 15) plus the step-change damper `|Δh| ≤ η h` (eq 16, typical
η = 15%) — note the lower bound is not decoration, §4.1 attributes the paper's headline
result to it. (d) §3.4's approximate Newton as the default, with the exact solve behind a
flag.

**Gate 12-1 (it is the method, not a rejection loop).** With the coupled path selected, a
run must contain **zero** discarded time points from LTE rejection — the step size is
solved, not retried. Declared success: the rejection counter is 0 where the standard path
records a nonzero count on the same circuit.
OUTCOME:

**Gate 12-2 (the paper's own result, on our circuits).** §4.1 reports **39% fewer time
points and 17% less runtime** on a Class-D amplifier with the LTE bounded between 0.7 and
3.0. Declared success: on the leapfrog and on `test_analysis_transient_stress.py`'s stiff
RLC, **at least 20% fewer accepted steps** than the standard controller at matched accuracy.
**A step reduction bought by a looser effective tolerance does not count** — accuracy is
measured against the analytic reference of stage 3's RC and must not degrade.
OUTCOME:

**Gate 12-3 (the overhead claim).** The paper's case rests on "very little overhead". With
stage 7b's factors reused, declared success: the per-step cost of the coupled path is
within **15%** of the standard path's. If it is not, the method is losing in wall time what
it wins in step count, and that is the number that decides whether it ships.
OUTCOME:

**Gate 12-4 (the four ignored inputs).** Whatever else changes, the coupled path must
honour `fixed_timestep`, breakpoints, an injected step controller, and `uic` — the list
from 0.1d, of which only `uic` is currently fixed.
OUTCOME:

**Gate 12-5.** Full suite `-m ""`, and the citation in `time_stepping.rst` becomes true —
the page's `.. warning::` recording that the method was documented but absent is removed
only when the code matches it.
OUTCOME:

**Reconsider the whole stage if** gate 12-3 fails, or if stage 4 shows the LTE estimate is
still not trustworthy enough to solve against. **A negative result here is a good outcome**:
it would justify deleting `_solve_coupled` after all, which is what 0.1d recommended and D1
declined, and it would do so on evidence rather than on the code-reading that D1 found
unpersuasive.

---

## Order and dependencies

```
0.1a ────────────────────────────────────► 9        [DONE]
0.1b ──────────────────► 5 ──────────────────────► 10.4 (large-signal MOSFET)  [DONE]
0.1c ──────────► 5 ────────────────────────────► 11 [DONE; 11 also needs 5]
0.1d ──────────► (delete _solve_coupled, in 1)      [DONE]
0.2a ──────────► 4                                  [DONE — passed, 4 unblocked on this]
0.2b ──────────► 0.3d                               [DONE — under threshold]
0.2c ─► (all suite gates)                    [DONE 2026-07-30: 734/6/0, 497.69 s]
0.3a ──► 1 (Newton/LTE tolerance separation)        [ANSWERED (iii), split by 0.3d]
0.3b ──► 4f only (default deferred to a post-repair measurement)  [ANSWERED]
0.3c ──► 10                                         [ANSWERED]
0.3d ──► 4 (chgtol + the guard), 6 (its diagnostic) [ANSWERED (D) — STAGE 0 EXITED]

1 ─► 2 ─► 3 ─► 4 ─► 5 ─► 6          (7 in parallel throughout)
                              8, 10.1-10.3 after 6;  10.4 after 5
                              11 after 5 AND 7
                              12 after 4 AND 7b   [added 2026-07-31, from D1]
```

**Stage 12** (implement Fang's method) is gated on stage 4 because it makes the LTE a
quantity to solve *for* rather than to test against, and on 7b because its claimed "very
little overhead" depends on reusing the factorisation. It is the expensive resolution of
decision D1; the cheap one — drop the citation and fix the four ignored inputs — needs
neither, and one of the two should be chosen before either is started.

It is placed **last among the stages** deliberately. Its dependencies only require it after
4 and 7b; putting it after 5, 6, 8, 9, 10 and 11 as well is a value judgement, and the
judgement is that those are capability and correctness items — 5 makes bipolar circuits
converge at all, 8 fixes sources that crash on their own class defaults, 10 adds a DC sweep
that does not exist — while stage 12 is efficiency on a path that is off by default. **The
one thing that would move it earlier** is a workload where step count dominates and none of
the correctness items block it; nobody has one on record.

**Amendments from the 0.3 decisions (2026-07-30), all recorded above:**
0.3a took option (iii), which moves the charge-referenced LTE criterion out of stage 1 and
into stage 4 and partly absorbs 0.2b. 0.3b declined the `'classic'` flip, which promotes
gate 4f from a formality to the measurement that sets the default. 0.3c put a
large-signal MOSFET in scope but sequenced it behind 0.1b and stage 5, and dropped the
netlist reader and waveform export. 0.2b's measurement plus 0.3a's answer produced the new
decision 0.3d, **answered (D)**: YWR's solution-domain criterion with a charge-domain
step-collapse guard, which splits 0.3a(iii) — the tolerance separation goes to stage 1, and
`chgtol` survives to serve the guard in stage 4 rather than as the accuracy test. 0.1c added
stage 5 as a prerequisite of stage 11. 0.1a and 0.1d each added a silent-wrong-answer defect
to stage 1's scope.

**Stages 1-3 are the ones to do if only one thread is available**: they stop the silent
failures, give 10.5x, and make `reltol` mean something. Roughly a week including gates and
documentation.

Each stage commits with its explanation and its measured gate outcomes in the message.
Negative results are recorded in the same voice as positive ones.

**Do not `git push`.**

---

## STAGE 1 — completed 2026-07-30

All four declared gates passed (1-1, 1-2, 1-3, 1-4), plus a fifth added for the tolerance
separation that decision 0.3a/0.3d required. Suite **744 passed, 6 skipped, 0 failed**.
Tests: `pycircuit/circuit/tests/test_transient_repairs.py`, one per gate. (Renamed
from `test_transient_stage1.py` once stages 2+ and 3 added to it.)

**Gate 1-5 (the tolerance roles are separate knobs)** — not in the original plan; added
because 0.3a's separation landed here and an unmeasured split is what created the problem
in the first place. Declared success: `lte_vabstol` moves the step count, `vabstol` does
not, and `Transient.vabstol == DC.vabstol`.
OUTCOME: **PASSED.** On a `VPulse`-driven RC whose stepping is genuinely LTE-bound (~1150
steps against a `tend/timestep` of 30), `lte_vabstol` gives **942 / 1155 / 1265 steps** at
1e-3 / 1e-6 / 1e-9 — monotonic — while `vabstol` at 1e-9, 1e-14 and 1e-6 all give **1155
steps and a waveform identical to eight decimals**. The two roles are now independently
observable, which they could not be before.

*A note on the test circuit, since it cost a false start:* a plain RC run at a `timestep`
near its own time constant is bound by `max_step` from end to end, so the controller never
decides anything and no tolerance can be seen to act. The first version of this gate used
one and read 20 steps for every setting — which would have "passed" a weaker assertion
while measuring nothing. This is stage 3's finding (`reltol` does not control accuracy
because the run opens at `max_step`) showing up as a testing hazard.

**Items completed:** (a) failed DC raises, (b) both `except Exception` clauses narrowed,
(c) `self.epar` threaded to all eleven `cir.C/q/i/G/u` sites plus `bypasstol` in
`defaultepar`, (d) inner `DC` built from the transient's own configuration, (e)
`solve_batched` raises instead of silently starting from zeros, (f) the JAX TLine ring
overflow is detected and raised (`TLINE_HISTORY_DEPTH`, replacing seven bare literals), (g)
`_solve_coupled` raises on retry exhaustion instead of livelocking.

**Three defects found while implementing, none of which were in the plan:**

1. **`Transient.__init__` accepted `toolkit=` and dropped it** — never forwarded to
   `Analysis.__init__`, so `Transient(cir, toolkit=X)` silently ran on `cir.toolkit`. It
   hid because callers pass the toolkit the circuit already has, making the two agree by
   coincidence. Fixed; item (d) depends on it.
2. **`_solve_coupled` ignores `uic`** — one of the four inputs 0.1d found it dropping, and
   invisible until now because ignoring `uic` and failing the DC produced the same zeros.
   Fixed.
3. **`doc/src/circuit/time_stepping.rst` documented a method that does not exist** — an
   augmented `(N+1)` system, an "Analytical Gradient" `E_h = p(E + TRTOL)/h`, and a "golden
   window" `0.7*tau <= eps <= 3.0*tau` attributed to Fang §3.3. **None of the three is in
   the code.** This also settles what the dead `analytical_eh` parameter was for. Corrected
   in the same commit, with the discrepancy recorded rather than deleted — a documented
   feature that does not exist is the same class of defect as a silent wrong answer.

**Negative result:** the plan's suggested `bypasstol=1e-2` for gate 1-3 is outside the
working range — with `bypass` genuinely connected, 1e-2 and 1e-3 both stop the inner DC
converging on a plain diode circuit. `bypass` had been dead long enough that no value was
ever validated. 1e-4 and below work; 1e-5 gives a 3.9x reduction in model evaluations.

**Not done, and deliberately:** `_solve_coupled` is repaired, not deleted. Deletion is
recommended by 0.1d and remains **an open maintainer decision** (12 tests, enumerated
above), so taking it here would have pre-empted a call that is not the implementer's.

---

## STAGE 2 — completed 2026-07-31

Suite **744 passed, 6 skipped, 0 failed**; doc build **succeeded, 2 warnings, 0 ERROR**.
Gates: 2a partly refuted, 2b passed, 2c half passed and half refuted on principle, 2-final
passed. Evidence: `benchmarks/transient_stage2.py`.

**Result: 2.42x bit-identical, 5.19x with BLAS held to one thread.** Assembly fell from
~96% of runtime to ~48%; attribute-lookup machinery from ~34% to under 5%.

**Three premises in the plan were wrong, and all three were wrong in the same direction —
the review's numbers were optimistic:**

1. The 4.462 ms threaded-solve figure does not reproduce (0.234 ms measured); the thread
   win is **1.72x, not 2.3-2.5x**.
2. The compound figure is **5.19x, not 10.5x**.
3. Gate 2c's "~2.17 assemblies needed" is not a redundancy figure at all — the solver's
   `(F, J)` is at the previous iterate, so the final assembly at the converged point is
   **necessary**, and 2.17 would require Newton to converge in 1.17 iterations.

The mechanisms were all correctly identified; only the magnitudes were overstated. That is
worth knowing before quoting any other figure from `transient_review.md` — see also the
420 GiB / "150 TB" correction under 0.1c.

**A gate defect worth carrying forward:** `drift == 0.00e+00` is the right criterion for
changes that reorganise our own arithmetic, and **unsatisfiable for any change that alters
which BLAS kernel runs**. Later stages that touch the linear algebra (7b, 7c) should declare
"identical step count plus drift at rounding level" instead, and verify determinism at a
fixed thread count separately — which is how 2a's drift was shown to be BLAS and not code.

**A near miss worth recording:** gate 2b's first version sampled `u` at `t = 0`, where the
sine source is zero, and reported "exactly equal" against an all-zero vector. It would have
passed while measuring nothing. The gate now prints the nonzero count beside every
comparison so an empty result cannot read as a pass.

### Left on the table, measured but not done

`Toolkit.__getattr__` still takes **3.77 million calls** in one benchmark run. These are the
*successful* delegations (`toolkit.zeros`, `toolkit.reshape`, ...) rather than the raising
ones stage 2b removed, so each is cheap — but under the profiler they are still ~5% of
runtime, and the fix is small: memoise the resolved attribute onto the instance so the
second and later lookups miss `__getattr__` entirely. **Not done here because it was not a
declared item and stage 2 is defined as behaviour-preserving; it deserves its own gate**
(the risk is toolkits that mutate their backend at runtime, which the memo would then
stale). Recorded so the next person does not have to re-profile to find it.

---

# STAGE 2+ — three improvements outside the original plan

Proposed after stage 2 and authorised by the maintainer on 2026-07-31, in this order.
They are numbered 2+.1 .. 2+.3 rather than folded into stage 2, because stage 2's gates
were declared before it ran and are closed; these get their own, declared here **before
any of the three is implemented**.

Items 1 and 2 are behaviour-preserving and inherit stage 2's stop condition: **waveform
drift exactly `0.00e+00` and an identical step count, or stop.** Item 3 deliberately
changes behaviour and is gated differently.

## 2+.1 — memoise `Toolkit.__getattr__`

Every `toolkit.zeros`, `toolkit.reshape`, `toolkit.dot` goes through `__getattr__` and
delegates to the backend module. Measured **3.77 million calls** in one benchmark run.
Resolve once, then cache on the instance so later lookups hit the instance `__dict__` and
never enter `__getattr__` at all.

*Risk assessment done before writing code:* a memo goes stale if the backend is mutated
after first access. Checked — **nothing assigns to `_backend` or to a backend module
attribute anywhere in the tree**; the only runtime mutations are assignments *onto the
toolkit instance* (`numeric.exp = ...` in two test/benchmark files), which shadow
`__getattr__` regardless and are therefore unaffected. Only successful lookups are cached;
failures are not, so an attribute that appears later is still found.

**Gate 2+.1a (correctness).** Full suite `-m ""`. Declared success: 744/6/0.
OUTCOME: **PASSED. 744 passed, 6 skipped, 0 failed** (796 s).
**Gate 2+.1b (behaviour).** Leapfrog benchmark. Declared success: drift **exactly
`0.00e+00`** and step count 324.
OUTCOME: **PASSED. Drift exactly `0.00e+00` on both `t` and `v`, 324 steps.**
**Gate 2+.1c (it did what it claims).** `Toolkit.__getattr__` call count must fall by at
least 10x on the benchmark, and end-to-end speedup recorded. Declared success: >= 1.05x
end-to-end — small, because the profiler says ~5% and the profiler exaggerates pure-Python
frames.
OUTCOME: **PASSED on both counts, and the call-count reduction is far larger than asked.**
`Toolkit.__getattr__` calls on one benchmark run: **3,774,764 -> 6,952, a 543x reduction**
against a declared 10x. The residue is first-access resolutions plus the probes that
legitimately miss (failures are deliberately not cached).

End-to-end, measured **interleaved against stage 2 at `2d7a2e7`** because the machine was
too noisy for a direct comparison -- min of 4, BLAS single-threaded: **13.128 s against
14.081 s, i.e. 1.073x**, 324 steps on both sides. Clears the 1.05x bar.

Note the gap between the profiler's verdict and reality: profiled total time fell
**171.9 s -> 74.2 s (2.3x)** while true runtime moved 1.07x. `__getattr__` is a
pure-Python frame, which cProfile charges far more heavily than the interpreter does. **The
profile was right about where the calls were and badly wrong about what they cost** -- a
reminder to confirm a profiler-driven optimisation against the clock before believing its
magnitude.

## 2+.2 — let the Newton solver return `(F, J)` at the converged point

`StandardNewton` evaluates `(F, J)` at `x` and returns `x_next`, so the caller must
re-assemble at the converged point (gate 2c established this is necessary, not wasteful).
Move that final evaluation *inside* the solver and return it. **The assembly count is
unchanged** — this buys no speed by itself. It exists because it makes the solver's
contract honest, and because it is the seam stage 7b needs in order to reuse LU factors
across a Newton iteration.

**Gate 2+.2a (behaviour).** Declared success: drift **exactly `0.00e+00`**, step count 324,
and the per-step assembly counts unchanged from stage 2's (`G`/`C`/`i`/`u` 3.06, `q` 3.06).
**An assembly count that *falls* here is a failure, not a bonus** — it would mean the
returned `(F, J)` is from the wrong state, which is exactly the defect gate 2c refused.
OUTCOME: **VOID — the gate was written for a design that investigation replaced. Superseded
by 2+.2a', declared below and before implementation.**

### 2+.2 REVISED, 2026-07-31, before any code was written

Reading the call sites to build the seam turned up something the proposal missed:
**`solve_timestep` returns `(x, feval, J, f)` and `f` is never used.** Verified in both
callers — `solve()` at `:536` and `_solve_coupled` at `:703` unpack it and the name never
appears again (the only later matches are f-strings). Only `J` is consumed, by the step
controller.

`func(x)` at the converged point computes `C`, `q`, `i`, `u` and `G`. The controller needs
`J = G + Geq`, which needs `G` and `C`; the charge cache needs `q`. **`i` and `u` are
computed and thrown away** on every accepted step, unless `provided_function` is set — that
is the one caller that consumes `f`.

So item 2 is not the contract-honesty refactor it was proposed as; there is a real saving,
and it is in the opposite direction from the original gate. **The seam idea is also weaker
than claimed**: the solver is handed `refnode_removed(func, ...)`, so any `(F, J)` it could
return is the *reduced* pair, while the caller needs the full `J` — reconciling those is
stage 7a's work (taking the reference node out of the matrix), not this item's. The seam is
therefore deferred to 7a, where it belongs.

**Gate 2+.2a' (REVISED).** Skip the discarded stamps at the converged point. Declared
success: drift **exactly `0.00e+00`** and step count 324 — this is still behaviour
preserving, because nothing computed changes, two things merely stop being computed. `i`
and `u` fall to **~2.06** per accepted step while `G`, `C` and `q` stay at **3.06**. With
`provided_function` set, all five stay at 3.06, because `f` is then genuinely needed.
OUTCOME: **PASSED.** Per accepted step on the leapfrog:

| stamp | before | after | with `provided_function` |
|---|---|---|---|
| `G` | 3.06 | **3.06** | 3.03 |
| `C` | 3.06 | **3.06** | 4.00 |
| `i` | 3.06 | **2.04** | 3.03 |
| `u` | 3.06 | **2.04** | 3.03 |
| `q` | 3.06 | **3.06** | 3.03 |

Drift **exactly `0.00e+00`** with the thread count matched to the reference, 324 steps.
`provided_function` restores the full evaluation exactly as declared. (`C` at 4.00 in that
column is pre-existing: `provided_function(f, J, cir.C(x))` assembles `C` once more, and
that is untouched by this item.)
**Gate 2+.2b (correctness).** Full suite `-m ""` at 744/6/0, and every `NonLinearSolver`
subclass still satisfying the interface.
OUTCOME:

## 2+.3 — `relref`, the reference for the relative tolerance

The highest-value item of the three and the only one that changes results. Spectre's
`relref` selects what the relative term is measured against; pycircuit hard-codes a
two-point `pointlocal`, which is why a node carrying no signal has its tolerance collapse
to the absolute floor — and why `lte_vabstol` had to be raised to 1e-6. See the Spectre
section under 0.3d.

Implement `pointlocal` (current behaviour, for exact backward compatibility) and
`sigglobal` (Spectre's default: each signal referenced to the maximum over all signals and
all past time). Default stays `pointlocal` in this item; **changing the default is a
separate decision** and belongs with stage 4's `lteratio` work.

**Gate 2+.3a (no silent change).** With `relref='pointlocal'`, drift **exactly `0.00e+00`**
and step count 324 against the stage-2 reference. The new code path must be provably
inert when not selected.
OUTCOME: **PASSED.** Drift exactly `0.00e+00` on both axes, 324 steps. `pointlocal` returns
`maximum(|x_curr|, |x_last|)` and never touches the running-maximum state.
**Gate 2+.3b (it addresses the actual defect).** With `relref='sigglobal'`, `lte_vabstol`
returned to **1e-12**, on the leapfrog. Declared success: the step count is within 2x of
the current default configuration — i.e. the 5.4x collapse that forced `lte_vabstol` to
1e-6 does **not** recur. **This is the whole point of the item; if it fails, `relref` does
not fix what it was proposed to fix and the result must be recorded as such.**
OUTCOME: **PASSED, and the defect it targets reproduces cleanly.** Leapfrog, step counts:

| `relref` | `lte_vabstol` | steps |
|---|---|---|
| `pointlocal` | 1e-6 (today's default) | **324** |
| `pointlocal` | 1e-12 | **1143** — the collapse, **3.53x** |
| `sigglobal` | 1e-6 | 324 |
| **`sigglobal`** | **1e-12** | **482 — 1.49x**, inside the declared 2x |

So the workaround's cost is real and reproducible, and `sigglobal` removes **81% of the
excess** (819 extra steps under `pointlocal`, 158 under `sigglobal`). The review recorded
the collapse as 5.4x; it measures 3.53x on this circuit — same phenomenon, smaller number,
consistent with the other magnitude corrections in stage 2.
**Gate 2+.3c (accuracy is not silently traded away).** `sigglobal` loosens the tolerance on
quiet nodes by construction. Declared success: on a circuit with an analytic solution, the
error under `sigglobal` at `lte_vabstol=1e-12` is no worse than under `pointlocal` at
`lte_vabstol=1e-6` — the configuration it is proposed to replace.
OUTCOME: **PASSED, but only on a restricted region, and the restriction is a finding.**

RC step response against the analytic `V(1 - exp(-t/tau))`, with a quiet 1 G-ohm node
present:

| `relref` | `lte_vabstol` | steps | err (whole run) | err (t > 3 tau) |
|---|---|---|---|---|
| `pointlocal` | 1e-6 | 177 | 1.980e-01 | 7.107e-02 |
| **`sigglobal`** | **1e-12** | **75** | 1.980e-01 | **6.929e-02** |
| `pointlocal` | 1e-12 | 177 | 1.980e-01 | 7.107e-02 |

`sigglobal` at 1e-12 is **slightly more accurate with 2.4x fewer steps** than the
configuration it replaces. But note the whole-run column: **identical to four digits across
every configuration.** That is stage 3's defect — the opening step is taken at `max_step`
and accepted unevaluated, and its error dominates the entire run, so no tolerance setting
can move the whole-run figure. **The gate as declared is therefore not measurable until
stage 3 lands**; it is judged here on `t > 3 tau`, the region the controller actually
governs, and that limitation is recorded rather than hidden. Re-run it after stage 3.
**Gate 2+.3d (correctness).** Full suite `-m ""` at 744/6/0 with the default unchanged.
OUTCOME:

---

## STAGE 3 — completed 2026-07-31

Gates 3-1 (**90.8x**) and 3-3 passed; 3-2 passed on the circuit it names (+3.1%) and failed
on a short RC (+32%), with the mechanism identified. Evidence:
`benchmarks/transient_stage3.py`.

**The defect, stated once with its measurement.** The opening step is the only one whose
error cannot be checked — the estimator differences against previous points and at the
start there are none — and it was being taken at `max_step`, the largest size allowed. So
the single unchecked step was also the biggest, and its error set the maximum for the whole
run: **1.3212e-01 at reltol 1e-3, 1e-4, 1e-5 and 1e-6 alike**, while the step count grew
from 24 to 195. Eight times the work, identical answer.

**A defect I introduced and the suite caught.** The first version ramped unconditionally,
including under `fixed_timestep=True` — where `dt` is never updated, so instead of opening
small and growing, the *entire* run proceeded at `timestep*1e-3`. Four tests failed
(`test_transient_RLC` and `test_transient_nonlinear_C` against their QUCS references,
`test_transient_methods_step_response` on shape, `test_Idtmod_modulo`). **This is what gate
3-4 is for**, and it is worth recording that the two QUCS comparisons — the only checks in
the suite against an external simulator — were among the four that caught it.

**Test churn, and how each was resolved.** Six failures, none left unexplained:

| test | cause | resolution |
|---|---|---|
| `test_transient_RLC` | my `fixed_timestep` bug | fixed the code |
| `test_transient_nonlinear_C` | same | fixed the code |
| `test_transient_methods_step_response` | same | fixed the code |
| `test_Idtmod_modulo` | same | fixed the code |
| `test_elements_module_doctests` | `Idtmod` doctest documents a uniform output grid | `fixed_timestep=True` in the doctest, with the reason |
| `test_step_count_and_error_respond_to_reltol[trap-ywr]` | see below | 5% slack on step count only |

**The one that is a real signal.** `trap-ywr` gives step counts `[55, 53, 105, 321]` — a
2-step dip between reltol 1e-3 and 1e-4 — while its **error falls monotonically**
(9.080e-4 → 8.279e-4 → 1.257e-4 → 7.768e-6). Errors are strictly monotone for **all six**
integrator/formula combinations; only this one's step count wobbles. That is not a
coincidence: `trap-ywr` is the combination gate 4d exists to fix, because YWR's Table I
gives TRAP as a **uniform-grid** formula (a single `h`, an unweighted second difference)
while its GEAR2 entry carries `h1`/`h2` explicitly — and stage 3's ramp is precisely what
introduces a phase of varying step ratio. The assertion now allows 5% on the step count
with that reasoning attached; **when 4d lands the slack should become unnecessary, and that
is the test to check it with.**

**On the Hairer estimate — the question the plan asked to be measured, now closed.** The
ramp costs a *fixed* ~15 steps (the geometric climb from `timestep*1e-3` back to
`max_step`, `log2(1000) ≈ 10` plus accept/reject overhead), so its relative cost is
+32% / +17.6% / +9.7% / +4.2% at 10 / 30 / 100 / 300 time constants, and +3.1% on the
leapfrog. A Hairer-style estimate from `q'`/`q''` would open near the right size and
recover most of those 15 steps. **Recommendation: do not implement it.** It buys a few
percent on realistic runs and matters only where startup dominates; revisit if a use case
appears where it does. Recorded so the question is closed rather than left open.

**Not done, and deliberately:** the `firststep` default of `timestep*1e-3` is now derived
(it is the knee of a measured curve — 1e-1 gives no benefit at all, 1e-2 gives 63.9x,
1e-3 gives 90.8x, smaller only costs steps), but it is still a *ratio to `timestep`* rather
than anything the circuit told us. A step derived from the circuit's own time constants
would be better and is the same work as the Hairer estimate.

---

## STAGE 4 (part 2) — 4g tier (a), 2026-07-31

Suite **755 passed, 6 skipped, 0 failed**.

**Tier (a) done: `Sin.next_event` returns `td` once, then `inf`.** A sine is C-infinity
after `td`; peaks and zero crossings are not discontinuities, and a breakpoint is for a
discontinuity. Order-1 evaluations on a `VSin`-driven RC: **120 of 1236 -> 3 of 1596**, and
the three are evaluations 1/2/3, the genuine opening of the run.

**Tier (a) is NOT sufficient, which is what the plan asked to be measured.** est/true for
trapezoidal still scales as 1/h (4.87 / 65.1 / 681.4 at h = 1e-8 / 1e-9 / 1e-10), because
the harness that measures it has **no breakpoints at all** — the alternating mode is
intrinsic to the companion recursion `iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}`, whose
homogeneous solution is `(-1)^n`. Breakpoints re-seed it; removing them does not remove it.

**Tier (b) is required and its design is settled** — difference `d_n = (q_n - q_{n-1})/h`,
in which the alternating component cancels exactly because
`(iq_n + iq_{n-1})/2 = d_n`. Cost: a third past charge and `h_last2` threaded through
`compute_lte`, its three implementations, both controllers and the transient's history.
**Deliberately not started at the end of a long session** — stage 3 already demonstrated
that a hasty change here breaks the two QUCS reference tests, which are the only external
check in the suite.

**One test was updated because it encoded the defect.** `test_func.py::test_sin` required
an event every quarter period. It now asserts the opposite contract, with the measurement
that justifies it. A real bug surfaced while doing so: the new `t < td` comparison is a
sympy relational under the symbolic toolkit and raised `TypeError`; it now degrades to "no
events", which is the only meaningful answer for a symbolic time.

### Stage 4 status

| item | state |
|---|---|
| 4a PI gains | not started |
| 4a-bis two-threshold controller | not started (hypothesis, from the papers) |
| 4b MAX_REJECT force-accept | not started |
| 4c Euler variable-step bias | **DONE** — spread 4.03x -> 1.01x |
| 4d trapezoidal `'ywr'` | **blocked on 4g(b)** |
| 4e `check_order_drop` direction | not started |
| 4f default `lte_formula` | last — needs 4a-4e; the margin is now only 1.038x |
| 4g(a) breakpoints | **DONE** — 120/1236 -> 3/1596 |
| 4g(b) mode-free estimator | **required, designed, not started** |
| 4h `fixed_timestep` | not started |
| gear2-classic ratio dependence | **NEW**, measurable now, not started |
| 0.3d `chgtol` guard | not started |
| `relref` default / `lteratio` | not started |

---

# MAINTAINER DECISIONS, 2026-07-31

Three items that had been left open across several tranches, answered together.

**D1 — `_solve_coupled`: KEEP AND FIX.** Deletion (recommended by 0.1d) is declined. The
livelock is already gone — stage 1 made it raise on retry exhaustion instead of advancing
time by `h*0.25^10` forever — so what "fix" now means is the rest of 0.1d's list, and it is
worth writing down so the scope is not rediscovered:

- it **ignores four inputs**: `fixed_timestep`, breakpoints (`next_event` is never called),
  and any injected step controller (`IntegralController()` is hard-coded at `:486`). `uic`
  was the fourth and is now honoured.
- `analytical_eh` is accepted and never read — a vestige of an `E_h` gradient that was
  documented but never written.
- the **Fang citation does not describe the code** (0.1d): there is no `p`, no bordered
  `(N+1)` system, no `J^{-1}p` correction and no two-sided LTE band. Either the citation
  goes or the method gets written.

**Reconsider-if:** if nobody intends to implement genuine co-determination, the honest
version of "keep" is to drop the Fang citation, rename the flag to describe what it does (a
rejection loop with a larger retry budget), and fix the four ignored inputs. That is much
less work than implementing the paper and leaves no false claim behind.

**D2 — `MOS_ACM`: KEEP, out of scope for now.** It remains unconstructable
(`super(MOS, self)` from a class whose MRO lacks `MOS`) and is a copy of `MOS` rather than
an ACM model. Nothing depends on it. **The risk of keeping it is that it reads as coverage
it does not provide** — its own docstring doctest would have caught the defect, and
`pytest.ini` configures no doctest collection, so it has never run. If it stays, that
doctest should either be collected or removed, so the file does not advertise a test that
cannot fail.

**D3 — `relref` default: `sigglobal`.** Adopted, per the measurements already recorded
under gate 2+.3b/c. This is the only one of the three that changes behaviour, so it is
gated below.

## D3 gates

**Gate D3-a (the default actually changed and nothing else did).** Full suite `-m ""`.
Expect churn in step counts; every failure explained individually.
OUTCOME: **FAILED, and the gate did its job. The default is reverted to `pointlocal`
pending 4g(b).** Suite: 4 failed, 751 passed. Three of the four are trapezoidal.

Under `sigglobal` the tolerance is referenced to the largest signal in the circuit, so
steps grow larger. On an estimator that is still contaminated by the `(-1)^n` mode, that
is enough to break the controller's response to `reltol` — **accuracy stops falling
monotonically as the tolerance tightens**:

| relref | integrator | steps | error |
|---|---|---|---|
| pointlocal | trap-classic | 51 / 54 / 83 / 149 | 1.01e-3 / 9.07e-4 / 2.10e-4 / 4.64e-5 |
| **sigglobal** | **trap-classic** | 69 / 66 / 53 / 87 | 4.89e-4 / **2.23e-4 → 3.14e-4** / 7.14e-5 |
| **sigglobal** | **trap-ywr** | 60 / 176 / 71 / 131 | 5.94e-4 / **3.82e-5 → 2.65e-4** / 3.18e-5 |

**Euler and both Gear2 variants are fine under `sigglobal`** — monotone in both columns,
and needing **1.7-2.5x fewer steps** for comparable accuracy. So the mode is right and the
*timing* is wrong: it must not ship while the default integrator's estimator is still
mode-contaminated.

**The decision is not overturned, it is sequenced.** `relref='sigglobal'` remains the
intended default; it lands with 4g(b), which is a one-line change once the trapezoidal
estimator is mode-free, with this gate already written to check it.

**Reconsider immediately if** 4g(b) is deferred indefinitely — in that case `sigglobal`
should still be adopted for Euler and Gear2 and the trapezoidal exposure documented,
because a 1.7-2.5x step reduction on two of three integrators is not worth withholding for
long.

**Gate D3-b (the workaround is no longer load-bearing).** With `sigglobal` as the default,
record what `lte_vabstol` is still worth. **This does not change `lte_vabstol` — that is a
separate decision**, and the point of the gate is to produce the number that decision needs.
OUTCOME: **Answered, and the answer is "nothing, on a circuit with a healthy signal."**
`test_gate_1_5_lte_vabstol_moves_the_step_count` failed under `sigglobal` with **403 steps
at 1e-3, 1e-6 and 1e-9 alike** — the absolute floor stops mattering entirely, because
`reltol*ref` dominates once `ref` is the largest signal rather than a possibly-zero local
one.

**That is the intended behaviour, not a defect**: an absolute tolerance exists to stop the
relative one degenerating, and under `sigglobal` it cannot degenerate. It does mean that
when `sigglobal` becomes the default, `lte_vabstol` can go back to something defensible
(1e-12, matching `vabstol`) at little cost, and that gate 1-5 will need rewriting — it
currently asserts a property that is only true under `pointlocal`. Both belong with 4g(b).

On the leapfrog, partial figures before the run was cut short by machine contention:
`pointlocal` 318 / 422 / 1206 steps at `lte_vabstol` 1e-6 / 1e-9 / 1e-12; `sigglobal` 318
at 1e-6. The `sigglobal` 1e-12 figure was 482 when measured at gate 2+.3b.
