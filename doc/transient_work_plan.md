# Transient subsystem — the work plan

> ## RESUME HERE — state as of 2026-08-01
>
> Last commit changing **default simulator behaviour**: `3a56c45` (9(g) — the JAX opening
> step and `dt_max`). Since then `def248c` (7a) and `433ead2` (7b) changed the CPU path
> bit-identically: 7b's sparse solver is **opt-in**, and the default is still
> `numpy.linalg.solve`.
>
> **The suite has 4 pre-existing COLLECTION ERRORS** — `pycircuit/post/cds/test/*` cannot
> import `pexpect` — which `pytest -q`'s tail does not show. Those four modules have never
> run; "N passed" excludes them. See traps 15 and 16.
> A newer HEAD does not by itself mean the code has moved underneath this block; check
> `git log --oneline 3a56c45..HEAD -- pycircuit/`. Branch `cna-jax-vectorization`
> (`git@github.com:henjo/pycircuit.git`) — **check `git status -sb` before assuming it is
> pushed**; several commits have sat unpushed at a time in this work.
>
> **Suite: 844 passed, 6 skipped, 0 failed** (`-m "" --timeout=400`). Runtimes this session
> ranged 676 s to 1941 s on near-identical source, entirely from other jobs on the box —
> see trap 2 before reading anything into one. Nominal
> ~8-13 min, but one run of the identical tree took **31m41s** purely from machine load —
> see trap 2. **Doc build: succeeded, 2 warnings, 0 ERROR, 0 degraded `exec-rst` blocks** —
> re-verified 2026-08-01 after D2, by reading the built page rather than the exit code: both
> generated tables in `lte_dae.html` hold live numbers. The two warning *lines*
> were read, not just counted: both are pre-existing autodoc failures importing
> `pycircuit.post.cds`. It was 3 warnings until 0.3d, when `example6.rst` was found broken
> by an earlier change in this session; see trap 11. Working tree clean.
>
> Transient regression tests live in `pycircuit/circuit/tests/test_transient_repairs.py`
> (renamed from `test_transient_stage1.py`); its docstring maps sections to plan stages.
>
> ### The next action, concretely
>
> **Stage 9 is COMPLETE** — (a) shared `_lte_kernels.py`, (b) tolerances settable, (c)
> per-row LTE tolerance + `sigglobal`, (d) the breakpoint hang, (e) non-converged Newton,
> (f) `lte_formula` removed, (g) the opening step + `dt_max`, and gates 9-1..9-3. Its own
> thesis held up: **every defect it found was a fix that existed on one backend and not the
> other** — 4i's Gear-2 constant, stage 3's opening step, 9(d)'s breakpoint scan — never a
> new bug.
>
> **Stage 7 is part-done.** 7a (reference-node removal, 2.0-4.6x, bit-identical) and 7b (a
> `LinearSolver` strategy: `DenseSolver` default, `SuperLUSolver`, `AutoSolver` selecting on
> **fill, not `n`**) are in. **Remaining: 7c (KLU) and 7d (delete `pybsmatrix.py`, and fix
> `test_sparse_toolkit.py`, which passes today while never exercising the sparse path).**
>
> **Recommended next: 7d before 7c.** 7d is cheap and removes a test that currently gives
> false assurance. 7c's declared ">= 3x on factor+solve" would land as roughly another
> 1.1-1.3x end to end, because **the solve is only ~8% of a transient** — 7b measured
> 1.09x at n=402, 1.20x at n=802, 1.25x at n=1202 against isolated solver wins of 28-74x.
>
> **Both open decisions are now TAKEN — see "DECISIONS TAKEN, 2026-08-01".** D1: factor
> reuse is **deferred**, not rejected — it saves ~2.7%% of runtime at today's sizes and costs
> bit-identity, and becomes worth ~20%% only once 2+.5 makes n=5000 circuits buildable. D2:
> `max_step` is now a parameter on both backends, defaulting to `None` (= `timestep`), so
> the default moves nothing while a quiescent run can ask for SPICE's `tend/50` and get
> **8.9x fewer steps at equal accuracy**.
>
> **Also queued, with gates already declared:** 2+.4 (assemble the reduced system directly)
> and **2+.5 (circuit construction is O(N^2.27) — 24.8 s to build an 800-element ladder)**.
> **Do 2+.5 first**: it is smaller, and 2+.4's gates are unmeasurable until a circuit that
> size can be built at all. 2+.5 is also why 7b's `n=2000` row was killed rather than
> measured.
>
> Still open by design: 5+.4 (large-signal MOSFET, sequenced into stage 10), the JAX Newton
> loop's scalar-norm convergence test (the per-row flavour split is done only for the LTE),
> and the P13/P5 architecture-review residuals. Stages 8, 10, 11, 12 untouched.

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

**These four are now WORK ITEMS WITH GATES, in `STAGE 5+` below.** They spent from
2026-07-30 to 2026-07-31 as prose inside this outcome — findings that nothing scheduled and
nothing checked — and stage 5 changed which of them binds. See there.

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
  maintainer's suggestion. **It was the recommendation from 2026-07-30 22:07 until its own
  entry measurements refuted it on 2026-07-31 — see the OUTCOME below. The standing
  resolution is (A).** Not the same as
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

**OUTCOME, 2026-07-31 — BOTH MEASUREMENTS TAKEN, AND (D) IS REFUTED. 0.3d RESOLVES TO
(A): NO CHARGE CRITERION, NO GUARD, NO `chgtol`.** The maintainer's suggestion is not
overruled by argument; it is refuted by its own entry measurements, which is why they were
declared. `chgtol` was never implemented, so this resolution is a decision not to write
code rather than a change to any.

*Measurement 1 — does the near-singular case arise?* Instrumented probe over five circuits
× two tolerances, 9510 accepted steps, computing both criteria without changing behaviour.
The guard fires **15 times**, all on `rc-pulse`. Taken alone this reads as a pass — and it
was recorded as one for about ten minutes. **It is not**, for the reason trap 6 exists:

- The two circuits built *specifically* to create the near-singular direction — a node
  grounded only through 1 GΩ, and a node held only by a 1 aF capacitor — **never fire it
  once**. The 1 GΩ appendage reproduced the unmodified circuit's step count and worst ratio
  to four significant figures, i.e. it perturbed nothing. A resistive appendage carries no
  charge, so `Eg` is zero on those rows and the guard cannot see them.
- The firings are on an ordinary RC pulse circuit, at the margin (`err_v` 1.49–4.16 against
  `err_q` 0.41–0.83), not in a distinguished regime.

*Measurement 2 — does the amplification separate the regimes cleanly?* **No. It is a
continuum spanning five orders of magnitude with no gap**, which is the assumption (D)
rests on, stated as such when (D) was written:

| circuit | median ‖lte‖/‖Eg‖ | max ‖lte‖/‖Eg‖ |
|---|---|---|
| stiff-rlc | 0.0064 – 0.066 | 1.4 – 2.0 |
| rc-vsin | 4.5 – 43.7 | 62 – 91 |
| weak-node | 5.7 – 53.5 | 91 – 247 |
| rc-pulse | 5.9 – 70.4 | 500 |

*Why it fired at all, which is the finding.* `err_v/err_q` is **not** the amplification.
It factors exactly as `(‖lte‖/‖Eg‖) × (etol_q/etol_v)`, and on **every one of the 15
firing steps `etol_q/etol_v` is exactly 0.01** — both tolerances sitting on their absolute
floors, where the ratio degenerates to `chgtol/vabstol = 1e-14/1e-12`. So the guard's
threshold is not a property of the circuit: it fires when `‖lte‖/‖Eg‖ > vabstol/chgtol`,
a number chosen by whoever last set the two floors. **pycircuit's `lte_vabstol = 1e-12` —
set by this plan — makes that threshold 10⁶ times more permissive than SPICE's `VNTOL =
1e-6` would.** A guard whose sensitivity moves by six orders of magnitude when an unrelated
tolerance is retuned is not a guard.

*And the criterion is dimensionally unsound.* `compute_lte` returns `Eg = -(h²/6)·q'''`.
With `q'''` in C/s³ that is **C/s — a current**, verified by scaling rather than by
derivation: halving `h` holds `Eg/h²` exactly constant (166667 = q'''/6 trapezoidal,
333333 = q'''/3 Gear-2) while `Eg/h³` doubles. (D)'s `err_q = max|Eg| / (TRTOL(reltol·max|q|
+ chgtol))` therefore **divides amperes by coulombs**; `err_q ≤ 1` is a test on a quantity
with units of 1/s. This is the same defect gate 0.2b recorded on the JAX backend — a
tolerance applied to a quantity of different units — arriving through a different door,
and it was about to be written into the CPU backend deliberately.

*The repair does not survive either, and this is what closes the question.* Comparing
`h·Eg` — a genuine charge — against `etol_q` is dimensionally correct and **fires on
28–99.7% of all rejected steps**:

| circuit | rejected steps | (D) fires | repaired guard fires |
|---|---|---|---|
| rc-vsin 1e-4 / 1e-6 | 376 / 3668 | 0 / 0 | 254 (67.6%) / 3069 (83.7%) |
| stiff-rlc 1e-4 / 1e-6 | 125 / 1235 | 0 / 0 | 36 (28.8%) / 351 (28.4%) |
| rc-pulse 1e-4 / 1e-6 | 1020 / 7992 | 5 / 10 | 463 (45.4%) / 7961 (99.6%) |
| weak-node 1e-4 / 1e-6 | 406 / 3853 | 0 / 0 | 190 (46.8%) / 3843 (99.7%) |

A guard that overrides the step controller on 99.6% of its rejections **is** the "never
rejects a step" failure mode 0.2b found on the JAX charge path. So the dimensionally wrong
version is inert-by-accident and the dimensionally right version is catastrophic; there is
no `chgtol` in between, because the quantity being thresholded is a continuum and the
scale separating the regimes was never physical.

**Resolution: (A).** One criterion — YWR, in the solution domain, which is what the paper
argues for. `chgtol` does not come back. The near-singular case is real but is not
detectable this way; it is already handled where it belongs, by stage 6's structural
singularity classification, which names the offending row instead of silently accepting a
step. **Reconsider if** someone exhibits a circuit where a step is rejected forever with a
well-conditioned charge residual *and* stage 6's diagnostic does not identify it — that is
the observation (D) was reaching for, and none of the five circuits here produces it.

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

OUTCOME (2026-08-02): **THE GUARD IS WITHDRAWN — its premise does not reproduce.**
`benchmarks/transient_review/stage0_3d_guard.py`.

The gate asks for "a circuit that needs it". There is not one. The premise was that a
near-singular `J` makes `J^{-1}` amplify the charge residual unboundedly, so the
solution-domain LTE becomes noise and the step collapses. Measured on four circuits:

| circuit | steps | fallback fired | max cond(J) | max abs lte | outcome |
|---|---|---|---|---|---|
| healthy RC | 212 | 0 | 1.00e+03 | 3.91e-02 | completes |
| series caps (no DC path) | — | 1 | — | — | **fails at the operating point** |
| dangling cap | — | 1 | — | — | **fails at the operating point** |
| weak ground (1 Tohm) | 30 | **0** | **1.00e+12** | 6.25e-02 | completes |

Two things fall out. **The truly degenerate cases never reach step control**: they fail at
the DC operating point, with the diagnostic stage 1 added — *"singular Jacobian: 'm' appears
in no equation, so nothing determines it — for a node that means no DC path to ground"*.
That is a better place to catch it and a better message than a step-control guard could give.

**And the genuinely near-singular case is handled.** At `cond(J) = 1.0e12` the run completes,
the LTE stays the same order of magnitude as the healthy circuit's, and the fallback fires
zero times. The step does not collapse.

So the work this decision created for stages 4, 6 and 9 — `chgtol`, the charge-flavoured
numerator from `compute_lte`, the guard, its diagnostic and its statistics counter — is
**not done and should not be done**. Option (A) alone is what shipped, and it is sufficient
on the evidence.

**Reconsider if** a circuit is found where the LTE solve genuinely fails or the step
collapses on a near-singular `J`. That is now detectable rather than silent: the fallback
warns (see below), so such a circuit will announce itself instead of quietly producing a
step size that is not error-controlled.

**One part of the decision IS done, and it was the part that mattered.** 0.3d also said to
"delete the bare `except` at `stepcontroller.py:59-62`, the unlogged half-(B) fallback".
The defect there was never the fallback — it is that `Eg` is a *current* and `J^{-1}Eg` a
*voltage*, so it silently substitutes one flavour for the other and then compares the result
against a voltage tolerance. It is now loud. It is not removed, because removing it would
turn a degenerate Jacobian into an exception raised from inside step control instead of from
the operating point, where the diagnostic is far better.

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
- **0.3d's recommendation REVERSED, (B) → (A).** *(Option (D) was added 31 minutes after
  this paragraph and superseded it for a day, so a reader arriving here between
  2026-07-30 22:07 and 2026-07-31 finds the plan contradicting itself. It no longer does:
  (D) was refuted by its own entry measurements and 0.3d is back at (A), for a stronger
  reason than this paragraph gives.)* The YWR paper exists to argue that the
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

**2b. Assembly: hoist the per-element probes, batch the scatter.** `_add_element_submatrices` / `_add_element_subvectors`.
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
OUTCOME: **PASSED, and superseded many times since.** The runtime this gate asked for was
never written down at the time; the suite has been run to completion after every stage from
4 onwards, most recently **943 passed, 6 skipped, 3 xfailed, 0 failed in 497 s**
(2026-08-02). Recorded now rather than left blank, because a blank gate is
indistinguishable from an unrun one — the point of this sweep.

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

### 4a scope and gates, written 2026-07-31 before the fix

**The plan's spectral radii are confirmed by derivation, not just reproduced.** Taking
`x_n = ln(h_n/h*)` and `err_n = (h_n/h*)^p` — which is what the harness drives — the PI
update `h_{n+1} = h_n · err_n^{-k_I} · (err_{n-1}/err_n)^{k_P}` linearises to

    x_{n+1} = x_n (1 - p(k_I + k_P)) + x_{n-1} (p k_P)
    char:    z^2 - (1 - p(k_I + k_P)) z - p k_P = 0

| gains | p | roots | radius |
|---|---|---|---|
| shipped (0.3, 0.4) | 2 | -1.1165, +0.7165 | **1.1165** |
| shipped (0.3, 0.4) | 3 | -1.7758, +0.6758 | **1.7758** |
| corrected (0.3/p, 0.4/p) | 2 | +0.8000, -0.5000 | **0.8000** |
| corrected (0.3/p, 0.4/p) | 3 | +0.8000, -0.5000 | **0.8000** |

matching the recorded 1.12 and 1.78. **The corrected radius is identical for both p**, which
is the entire reason the gains are written `0.3/k` and `0.4/k` rather than as two smaller
constants: dividing by the order makes the closed loop order-independent. A fix that merely
picked smaller numbers would pass gate 4a at one order and fail at another, so that is
gated separately below.

**A claim in the plan that is REFUTED.** 4a states "The computed `exponent = 1.0/p` at
`:139` is dead code." In the current file it is line 244 and it **is used**, at line 248, on
the rejection branch — and the same is true of `IntegralController`'s copy, used on both
branches. Either the claim was true of an older revision or it was mistaken; it is not true
now, and nothing is deleted on its account.

**The rejection path, and why the fix is one line.** After a rejection the controller
returns without touching `self.last_err`, so the next accepted step differences against an
error two steps stale — and worse, one measured at a *different* step size at the same time
point. Textbook practice (Hairer & Wanner II.4, Gustafsson) is to take the step after a
rejection with the elementary (pure-I) controller and resume PI afterwards. That behaviour
already exists here: `if self.last_err is None: self.last_err = err` makes the P factor
`(err_last/err)^{k_P}` exactly 1. So setting `last_err = None` on the rejection path *is*
the pure-I fallback, using machinery already present rather than adding a mode.

**Gate 4a-1 (the declared one).** On a smooth problem h settles to within 5% of a fixed
point instead of alternating 2:1. Recorded before: a permanent period-2 cycle,
h = 0.8572 / 0.4286, spread exactly 2.0000 — i.e. running against the growth clamp every
other step.
OUTCOME: **PASSED.** h converges to the fixed point: last eight values 0.9993 .. 0.9999, tail spread **1.0006** against the declared 5% and against the previous permanent 2.0000.

**Gate 4a-2 (the fix is order-independent, not just smaller).** The closed-loop spectral
radius must be **identical for p = 2 and p = 3** and below 1 for both. This distinguishes
dividing by the order from picking two smaller constants, which would pass 4a-1 at whichever
order it was tuned on.
OUTCOME: **PASSED, exactly.** Driving the real update law from the same start, the observed convergence rate is **0.8000 at p = 2, 3 and 4** — identical to four figures and matching the analytic root. A fix that had merely chosen smaller constants would have shown three different rates.

**Gate 4a-3 (a rejection does not leave a stale error behind).** After a rejected step, the
next accepted step must use the elementary update rather than differencing against an error
measured at a different step size. Verified by instrumenting the factor actually applied,
not by reading the code.
OUTCOME: **PASSED.** With `last_err = 0.9` and `err = 0.5` at p=3 the factor is **1.159149**; after a rejection it is **1.071773**, which equals `err^(-k_I/p)` to machine precision — the P term is exactly 1, i.e. the elementary update.

**Gate 4a-4 (end to end, and it must not be worse than the default controller).**
`PIController` on the gate-4b sweep. Declared success: no accepted step ratio outside the
zero-stability bound; and against `IntegralController` on the same circuits, **rejections no
more than 1.5x and step counts no more than 1.2x**. A PI controller that is stable but
needs half again as many steps is not obviously an improvement, and if that is what it does
the honest outcome is to say so rather than to ship it as a win.
OUTCOME: **The stability clause PASSED; the rejections clause FAILED at 5.00x against a
declared 1.5x — and the clause was measuring the wrong comparison, which is a criticism of
how it was declared, not a reason to move it.** Both readings are below.

*What passed.* `worst accepted step ratio 2.0000, ratios outside the bound 0` across all
nine configurations, so `PIController` respects 4b's invariant. Step counts are **1.00-1.13x**
those of `IntegralController`, inside the declared 1.2x.

*What failed, with the raw numbers rather than the ratio.* `PIController` vs
`IntegralController`, rejections:

| circuit, reltol | integral | PI |
|---|---|---|
| rc-vsin 1e-3 / 1e-4 | 0 / 0 | 0 / 0 |
| **rc-vsin 1e-5** | **0** | **4** |
| stiff-rlc 1e-3 / 1e-4 / 1e-5 | 4 / 2 / 2 | 4 / 2 / 3 |
| rc-pulse 1e-3 | 9 | 8 |
| **rc-pulse 1e-4** | **36** | **24** |
| **rc-pulse 1e-5** | **48** | **61** |

The 5.00x is the `rc-vsin 1e-5` row under the declared `(n+1)` smoothing: 0 -> 4 rejections
in a **238-step** run, i.e. 1.7% of steps. The smoothing does not rescue a ratio when the
baseline is zero, and the gate should have used an absolute floor. PI is better on one row
(36 -> 24), worse on two, and identical on the rest.

**The substantive finding, which the clause obscured: `PIController` is not an improvement
over the default `IntegralController` on these circuits.** It is comparable, at 0-13% more
steps. So 4a fixes a defect in an opt-in path and **does not** justify promoting PI to the
default — and no such promotion is made.

*The comparison that actually judges 4a*, taken because the declared one answers a
different question — PI **before vs after** the gains fix, same circuits, same integrator:

| circuit, reltol | steps | rejections |
|---|---|---|
| rc-vsin 1e-3 | 214 -> 214 | 0 -> 0 |
| rc-vsin 1e-4 | 215 -> 218 | 0 -> 0 |
| **rc-vsin 1e-5** | 237 -> 238 | **34 -> 4** |
| stiff-rlc 1e-3 | 61 -> 66 | 5 -> 4 |
| **stiff-rlc 1e-4** | 87 -> 94 | **13 -> 2** |
| **stiff-rlc 1e-5** | 142 -> 150 | **29 -> 3** |
| rc-pulse 1e-3 | 493 -> 475 | 11 -> 8 |
| **rc-pulse 1e-4** | 691 -> 650 | **35 -> 24** |
| **rc-pulse 1e-5** | 979 -> 975 | **132 -> 61** |

**Up to 9.7x fewer rejections at essentially unchanged step counts.** That is the limit
cycle being paid for: a controller alternating against its growth clamp every other step
overshoots the tolerance half the time, and each overshoot is a rejected step and a wasted
Newton solve. The fix removes the oscillation and the rejections go with it.

**Gate 4a-5.** Full suite `-m ""`; doc build verified by content.
OUTCOME: **PASSED. 788 passed, 6 skipped, 0 failed, 676.30 s** (`-m "" --timeout=400`),
against 784 before plus the 4 tests added here — no existing test needed changing, which
is expected rather than reassuring: `PIController` is opt-in and the suite exercises it in
one place. Doc build: **build succeeded, 2 warnings, 0 ERROR**, verified by content.

**That is also the honest limit of this gate.** A full suite is weak evidence for a change
to a controller almost nothing uses; the load-bearing measurements for 4a are the harness
(4a-1, 4a-2), the unit checks (4a-3), and the before/after sweep under 4a-4.

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

**THE HYPOTHESIS IS NOW UNTESTABLE AS STATED, and answered another way (2026-07-31).** The
57-vs-29 gap it exists to explain was between `Gear2('ywr')` and `Gear2('classic')`, and
after 4i those two are the *same code* — `lte_formula` selects nothing for either
second-order method except on one fallback step. There is no gap left to attribute. The
four-way comparison gate 4f asked for has no four distinct configurations to run.

What the intervening work says about the underlying question — *is the single accept
threshold driving rejection counts?* — is that it was not the dominant term. Rejections fell
by 33x on the stiff case from fixing the trapezoidal **estimator** (4g(b)), by up to 9.7x on
`PIController` from fixing its **gains** (4a), and Gear2's rejections *rose* slightly when
its estimator was made honest (4i) because the old one under-reported. Three separate causes,
none of them the threshold.

**Reconsider if** a circuit is found where rejections stay high with every estimator now
correct — `F_redo` is still a real feature both papers have and pycircuit does not, and the
argument for it (a step 1.01x over tolerance is re-solved from scratch where both papers
would accept it and shorten the next step) is untouched by any of the above. It is simply no
longer justified by the 57-vs-29 evidence, because that evidence has evaporated.

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
OUTCOME: **CLOSED AS NOT ANSWERABLE AS POSED — there are no longer four combinations.**
4g(b) and 4i moved both second-order estimators onto a shared `q'''` helper that does not
read `lte_formula`, so `Trapezoidal('ywr')` and `Trapezoidal('classic')` return bit-identical
values, and so do the two `Gear2` variants. A four-way comparison has two distinct members.

Measured, CPU path, with the four charges a normal step has:

| integrator | `classic` vs `ywr` |
|---|---|
| Euler | bit-identical (the two formulas coincide, and always did) |
| Trapezoidal | bit-identical |
| Gear2 | bit-identical |

**The question this gate existed to settle — which `lte_formula` should be the default —
therefore has no measurable content left.** 0.3b's deferred decision is subsumed. What
replaced it is a different question, recorded below and staged.

### The `lte_formula` decision, 2026-07-31

**Where it still had an effect before this change.** Exactly one place on the CPU: the
single fallback step of a run, before the ring buffer holds four charges. There Euler is
identical, Trapezoidal differs only in the low bits (same formula, different operation
order), and **Gear2 differs by exactly 4/3** — the recorded YWR optimism, its GEAR2 residual
estimating `(1/4) h^2 q'''` against a true `(1/3)`.

**And a second consumer that none of this work touched.** `jaxtransient.py` carries its own
`lte_formula`, and there it still selects a genuinely different algorithm: `'ywr'` maps a
g-difference through `J^-1`, `'classic'` takes a charge-domain estimate with no `J^-1`.
Gate 0.2b established that the charge path applies `lte_abs = 1e-6` — a **voltage**
tolerance — to a **charge**, so it never rejects a step. **On the JAX backend
`lte_formula='classic'` selects a broken estimator.** Verified 2026-07-31 that no commit in
this session modified that file.

So one parameter name meant a real choice with a known defect on one backend and almost
nothing on the other. Staged rather than settled in one move:

**(D) NOW — remove the last place the choice bites on the CPU, and finish 4d.** Both
second-order fallbacks take the divided-difference form unconditionally, and the dead `ywr`
fallback branches go. For Gear2 because `ywr` is 3/4 optimistic there; for Trapezoidal
because YWR's TRAP entry is a **uniform-grid** formula (established from the paper under 4d)
being used on a grid the opening ramp makes non-uniform. This is exactly 4d's stated fix —
*"delete the branch and keep `'ywr'` as an alias"* — so **4d is closed by it**.

**(B) NOW — keep the parameter, and say what it is.** It stays accepted, so nothing breaks,
and the docstrings state plainly that it selects nothing on the CPU path and why. A knob
that advertises a choice it does not make is the "thin advertised feature" 0.1c warns
about; documenting it is the honest interim, and removing public API for a knob with no
remaining CPU effect is not worth a breaking change on its own.

**(C) WITH STAGE 9 — remove it from both backends.** Stage 9 already owns merging the two
transient paths, 0.2b's decision rule already says `lte_formula` should "select an algebra
only", and the JAX charge path's tolerance defect belongs to the same work. Doing it there
is one change instead of two, and avoids leaving the name meaning different things on the
two backends in the interim. **Recorded as a stage-9 item so it is not lost.**

**Gate 4f-D1 (the parameter is inert on the CPU, end to end).** Not at the estimator —
end to end. Two full transient runs differing only in `lte_formula` must produce
**bit-identical** waveforms and step counts, for every integrator, on a circuit whose
opening ramp makes the early grid non-uniform. Declared at bit-identical rather than "close"
because anything else means a branch is still live somewhere.
OUTCOME: **PASSED — bit-identical, not merely close.** Two full runs differing only in `lte_formula`, on a `VPulse` circuit whose opening ramp makes the early grid non-uniform: euler 403/403 steps, trapezoidal 153/153, gear2 183/183, and `max |dv| = 0.000e+00` for all three.

**Gate 4f-D2 (the fallback got the better formula, not just a consistent one).** Gear2's
fallback estimate must change from the `ywr` value to the `classic` one — a factor 4/3 —
rather than the two merely being made to agree on whichever was there.
OUTCOME: **PASSED.** The Gear2 fallback estimate is now **-3.770212e+01**, the divided-difference value, where `'ywr'` gave **-2.827659e+01** — the 4/3 the YWR GEAR2 residual is short by. The two were made to agree on the *accurate* one.

**Gate 4f-D3 (no regression).** Full suite `-m ""`, and
`benchmarks/transient_stage4.py --forceaccept` unchanged except where the fallback step
moves a run.
OUTCOME: **PASSED. 788 passed, 6 skipped, 0 failed, 1372.99 s** (`-m "" --timeout=400`)
— the same tally as before, with three tests rewritten rather than added. The 1373 s is
machine load, not this change (trap 2 records a 31m41s run of an identical tree); the
previous run of nearly the same source took 676 s. Doc build: **build succeeded, 2
warnings, 0 ERROR**, and the generated `classic`/`ywr` columns are now identical, which is
the page demonstrating gate 4f-D1 rather than asserting it.

**And it did far more than "move a run" — see below.**

### THE FORCE-ACCEPT PATH IS NOW UNREACHABLE

`--forceaccept` after this change records **zero force-accepts across all 45
configurations** — 3 circuits x 3 tolerances x 5 integrator/formula combinations — where
before this session it was reached up to 159 times in a single run. The progression, all
measured on the same sweep:

| after | worst force-accepts in one run | accepted ratios outside the bound |
|---|---|---|
| (before this session) | 159 | **61** |
| 4b + 4e | 121 | 0 |
| 4g(b) | 26 | 0 |
| 4i | 1 | 0 |
| **4d / 4f-D** | **0** | **0** |

**One step per run was seeding the rest.** The fallback estimate ran once, before the ring
buffer held four charges — and under `'ywr'` it reported 3/4 of the truncation error. An
under-report there lets the controller grow the step further than the tolerance allows, and
the run then spends rejections climbing back out. Correcting that single step removed the
last configuration that reached the rejection cap at all.

**This is the more useful reading of 4b.** Its escape hatch was never the disease: it was
the visible symptom of estimators that under-reported, and it fired progressively less as
each of them was corrected — 4c, 4g(b), 4i, and now 4d. The hatch stays, and its warning
stays, because an estimator can always be wrong on a circuit nobody has run yet; but on
everything measured here it is now dead code, which is what a correct error controller
should make it.

**Consequence for the 4b regression test**, whose `n_forced > 0` precondition existed so it
would fail loudly rather than pass against an empty result: it did exactly that, and there
is no configuration left to re-point it at. The precondition is **inverted** — it now
asserts `n_forced == 0` and says what to do if that ever stops being true. The mechanism
itself is still tested deterministically through `_RejectionInjector`, which was built for
this eventuality.

**Test churn, three tests, all pinning the removed distinction:**

1. `test_compute_lte_order_and_scale[gear2-ywr]` pinned the ratio at **1/2**; it is now 2/3,
   the same as `gear2-classic`, because the fallback both selections reach is the same code.
2. `test_lte_formula_ywr` asserted `classic/ywr = 4/3` as *"what still separates the
   formulas"*. Nothing separates them now, and that is the assertion instead.
3. `test_gate_4b_no_accepted_step_ratio_leaves_the_stability_bound` — the precondition
   above.

None was adjusted to make a number fit: each pinned a real property that this change
deliberately removed, and each now pins what replaced it.

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

### 4h scope and gates, written 2026-07-31 before the fix

**The mechanism, and it is one line.** `dt` is a loop-carried variable that breakpoint
truncation *overwrites* (`dt = float(next_t_break - t)`), and the restore at the bottom of
the loop is guarded by `if not fixed_timestep`. So under fixed-step a truncation is
permanent, and because each subsequent breakpoint truncates the already-shrunken `dt`
again, the step collapses geometrically. Reproduced on a `VPulse` at `tend=3e-6`,
`timestep=1e-7`: **292 steps against an expected 30**, final `dt` **1.241e-19 s**.

**A second defect, found while reproducing the first and NOT in the plan.** With no
breakpoints at all — a `VSin` drive, whose `next_event` returns `inf` since 4g(a) — a fixed
run still ends with a degenerate step: `41` steps for an expected 40, final `dt`
**2.033e-20 s**. A uniform grid divides `tend` exactly, so `tend - t` at the end is pure
floating-point residue, and `if t + dt > tend: dt = tend - t` turns that residue into a
step. Measured on the adaptive path for comparison: **0 degenerate steps** on both
circuits, because a controller-chosen `dt` does not land on `tend` to within 1e-20. So this
one is specific to fixed-step and the fix does not touch the adaptive path.

**The design question the declared gate settles.** `tend/timestep +- 1` cannot be met while
still truncating at breakpoints: a `VPulse` at `per=1e-6` over `tend=3e-6` has ~12 edges, so
truncate-and-restore would give `tend/timestep + 12`. **So under `fixed_timestep` the grid
wins and breakpoints do not move it.** That is the right reading of the flag — it exists so
a caller can ask for exactly these output points (an FFT-friendly grid, or comparison
against a reference at known times), and a grid that is silently not uniform is the defect,
not the feature.

**What is kept.** Crossing a breakpoint *inside* a step still drops the integration order
for the following step. That costs nothing, needs no grid change, and is the part of
breakpoint handling that protects against fitting a polynomial across a discontinuity.
**Reconsider if** someone needs edges resolved exactly under a fixed step — the answer then
is the adaptive path with `max_step`, not a fixed grid that quietly is not one.

**Gate 4h-1 (the declared one).** A `VPulse`-driven fixed-step run takes `tend/timestep`
+- 1 steps. Recorded before: 292 against 30.
OUTCOME: **PASSED. 292 -> 30 steps against an expected 30 (`tend=3e-6`, `timestep=1e-7`), i.e. exact rather than within the declared +-1.** The `VSin` control, which has no breakpoints at all since 4g(a), goes 41 -> 40.

**Gate 4h-2 (the grid is actually uniform).** Every step of a fixed run equals `timestep`
to within rounding, except at most the last, which must be in `(0, timestep]`. Declared
because 4h-1 counts steps and a step count can be right while the spacing is not — the
`VSin` case above had a step count inside +-1 *and* a 2e-20 s step.
OUTCOME: **PASSED.** Max deviation from `timestep` across the whole VPulse run is **1.06e-22 s** on a 1e-7 s grid, and the smallest step is a full `timestep` rather than the 1.241e-19 s the run used to end on. The `VSin` control is uniform to 2.75e-21 s and its 2.033e-20 s final step is gone. **This gate earned its separate existence**: the `VSin` case passed 4h-1 before the fix (41 steps against 40, inside +-1) while ending on a step 14 orders of magnitude below the others.

**Gate 4h-3 (the discontinuity is still handled).** A breakpoint crossed within a fixed step
must still force the order drop on the next step. Verified by instrumenting the effective
integrator, not by inspection.
OUTCOME: **PASSED — 29 order drops in 89 Newton evaluations** on the fixed VPulse run, against 0 if breakpoints had stopped being noticed and 89 if the run had gone first-order throughout. Both bounds are asserted, because either would be a way for this to look fine while being wrong.

**Running this gate changed the code, which is the point of running it.** The first version tested `next_t_break < t + dt` and scored 23. A breakpoint landing *exactly* on a grid point is not a rounding curiosity when the grid is uniform: with `td = 1e-7` against a `1e-7` grid, **2 of the 9 edges in a 30-step run land exactly on one**, and a strict `<` gives them no order drop at all — the step ending on the edge does not see it, and the next iteration asks `next_event(t)` from the edge itself, whose fixed-point guard skips past it. `<=` catches both; the count goes 23 -> 29 and the step count is unchanged at 30.

**Gate 4h-4 (a fixed grid that cannot be honoured says so).** `transient.py` shrinks `dt` by
0.25 on a `NoConvergenceError` and retries — which is the only way to make progress, but
under `fixed_timestep` it silently abandons the grid the caller asked for. Declared: that
path emits a `RuntimeWarning` under fixed-step naming `t` and the step it fell back to.
Same failure class as 4b's force-accept, and the same treatment.
OUTCOME: **PASSED.** With a `NoConvergenceError` injected at the fifth step, exactly one `RuntimeWarning` is raised, naming the time and the step fallen back to: *"Newton did not converge at t=4e-07 s with the requested fixed timestep 1e-07 s; falling back to 2.5e-08 s for this step. The output grid is no longer uniform."* The shrink itself is kept -- it is the only way to make progress -- but it is no longer silent.

**Gate 4h-5.** Full suite `-m ""`; doc build verified by content.
OUTCOME: **PASSED. 784 passed, 6 skipped, 0 failed, 805.34 s** (`-m "" --timeout=400`),
against 779 before plus the 5 tests added here — **no existing test needed changing**,
which is worth noting for a change that alters the step sequence of every fixed-step run.
The reason is that the two QUCS reference tests use `fixed_timestep=True` with a `VSin`
drive, whose `next_event` has returned `inf` since 4g(a), so they never had a breakpoint
to truncate on — only the degenerate final step, whose removal changes no sample they
assert against. Runtime 670 -> 805 s is machine load, not this change (trap 2): the
fixed-step runs in the suite got *shorter*.

Doc build: **build succeeded, 2 warnings, 0 ERROR**, verified by content.

**Gate 4-final.** Full suite `-m ""`, runtime recorded. **Expect test churn here** — this
stage changes step counts by design, and any test asserting a step count is exposed.
OUTCOME: **PASSED.** The churn the gate predicted did arrive and is recorded in the stage-4
sections themselves. Same note as gate 3-4: the runtime was not written down at the time and
the suite has been green after every stage since — currently **943 passed, 6 skipped, 3
xfailed, 0 failed in 497 s** (2026-08-02).

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

### Stage 5 entry measurement and sub-gates, written 2026-07-31 before any code

**The entry condition holds exactly as recorded.** A common-emitter NPN stage — `VCC` 5 V
through 1 k to the collector, `VIN` 0.8 V through `RB` to the base, emitter grounded — fails
to reach a DC operating point at **all three** base resistances:

| RB | result |
|---|---|
| 10 | `NoConvergenceError`: Source Stepping failed at lambda=1.0 |
| 100 | `NoConvergenceError`: Source Stepping failed at lambda=1.0 |
| 1000 | `NoConvergenceError`: Source Stepping failed at lambda=1.0 |

**And the mechanism is the one the plan assumes, verified rather than inferred.**
Instrumenting the model over the failing solve at RB = 100:

    model evaluations           405
    non-finite i() returns       12
    non-finite G() returns       12
    evaluations with v_be/VT > 709 (exp overflow)   2 of 405

`exp()` overflows to `inf`, the Ebers-Moll expression then computes `inf - inf` = `nan`, and
once `nan` reaches `x` every subsequent evaluation is `nan` — which is why the recorded
`v_be` range degenerates entirely. **Note the operating voltages themselves never overflow**:
`exp(5.0/VT) = exp(193)` is finite. The excursion comes from Newton steps, not from the
supply, which is why a clamp on the *model* is needed and not merely a sensible bias point.

**Two fixes are indicated and they are measured separately**, because applying both at once
would make neither attributable — the discipline this plan applies everywhere else:

**Gate 5-0a (the clamp alone: what does removing `nan` buy?).** Apply only the `_expl()`
clamp — SPICE's `EXPMAX` treatment — and re-run the three-resistance sweep. Declared: record
convergence and iteration counts. **No pass/fail**: this exists to attribute the two fixes,
and either answer is informative. If the clamp alone converges, the limiter's justification
changes from "makes it work" to "makes it converge faster", and that should be said.
OUTCOME: **ANSWERED, and the answer is worse than either branch the gate anticipated.** With `_expl` alone and no limiting:

| RB | result | v(c) | v(b) |
|---|---|---|---|
| 10 | **CONVERGED** | **-201.0 V** | **-200.1 V** |
| 100 | `NoConvergenceError` | - | - |
| 1000 | `NoConvergenceError` | - | - |

Non-finite `i()` and `G()` returns fall from 12 to **0** on all three, so the clamp does exactly what it claims. But it does not fix convergence, and at RB=10 it produces something worse than a failure: a **converged answer 200 V below the rails**.

**That answer is not a solver bug — it is a genuine solution**, and checking this mattered. The base-collector junction sits forward-biased at 0.91 V passing `IS*exp(0.91/VT) = 19.8 A`, against the `20.09 A` the base resistor delivers; the residual is 4.7e-12. It satisfies the circuit equations. It is simply a spurious operating point that no physical stage would occupy, and **only limiting keeps Newton out of it**. (An intermediate reading here claimed the solver had returned a non-solution; that was wrong — it checked only the base-emitter junction and missed the 20 A flowing in the base-collector one.)

So the limiter's justification is neither of the two the gate offered: it is not that the clamp is sufficient, nor that limiting merely speeds things up. **The clamp removes `nan`; limiting is what keeps the answer physical.** Both are needed and they do different jobs.

**Gate 5-0b (`limit` must return the limited vector).** `SubCircuit.limit` does
`subx = x[self.elementnodemap[instance]]`, which is fancy indexing and therefore **a copy**;
it calls `element.limit(subx, ...)` and discards the result, then returns the unmodified
`x`. So a limiter that writes its argument has no effect whatsoever. `Diode` does not notice
because it keeps `_vlim` privately and linearises `G` around it — the only limiter in the
tree is the one that cannot expose the bug.

Item 5(d) asked whether `limit` should become state-free. **It is not optional: gate 5-3
cannot be satisfied without it**, because "a test that fails if the return value is
discarded" presupposes there is a return value. Declared: `Semiconductor.limit` returns the
limited sub-vector, `SubCircuit.limit` writes it back, and `Diode`'s existing `None`-returning
form keeps working unchanged during the transition. Migrating `Diode` is a separate decision
with its own measurement — it is the only device whose `G` reads `_vlim`, so moving it
changes numbers on circuits that already converge.
OUTCOME: **IMPLEMENTED AS DECLARED.** `Semiconductor.limit` returns a limited copy; `SubCircuit.limit` writes it back with `x[nodemap] = limited`, which reaches `x` where the fancy-indexed *read* did not. `Diode`'s `None`-returning form is accepted unchanged and has its own test, so the transition breaks nothing. Migrating `Diode` off `_vlim` stays a separate decision — it is the only device whose `G` reads that state, so moving it changes numbers on circuits that already converge, and there is no reason to spend that in the same change.

**On `ZenerDiode`, per review 0.1b:** it does **not** get a plain `pnjlim`. Its reverse
breakdown term `-IBV*(exp((-v-BV)/VT) - 1)` is a second exponential in the opposite
direction, and a limiter that only knows about the forward junction steps straight through
the breakdown knee. It gets the `_expl()` clamp, which is direction-agnostic, and its
junction list stays empty until someone measures a breakdown limiter.
OUTCOME: **HONOURED.** `ZenerDiode.junctions` is left empty, so it gets the direction-agnostic `_expl` clamp and no `pnjlim`. `BJT` and `JFET` get junction lists. Recorded because it is a deliberate omission rather than an oversight: a forward-junction limiter would step straight through the breakdown knee, which is a second exponential running the other way.

**Gate 5-1 (the failing case converges).** A common-emitter BJT stage with a
voltage-driven base must reach a DC operating point at base resistances 10, 100 and 1000
ohm. Currently all three fail.
OUTCOME: **PASSED.** All three converge, in **17 model evaluations** against 405 and 421 spent failing:

| RB | v(b) | v(c) |
|---|---|---|
| 10 | 0.7301 | 0.0260 |
| 100 | 0.7038 | 0.0522 |
| 1000 | 0.6961 | 0.1169 |

`v_be` of 0.70-0.73 V is a forward-biased silicon junction and the collector sits in saturation, which is what 0.8 V of drive through these resistances should give. **The test asserts the bias point, not just convergence** — gate 5-0a showed that "it converged" is satisfied by an operating point 200 V below the rails.

**Gate 5-2 (no `nan` reaches the Jacobian).** Instrument for non-finite entries during the
continuation. Declared success: none, at any gmin or source-stepping rung.
OUTCOME: **PASSED. Zero non-finite `i()` or `G()` returns**, across the whole continuation, against 12 of 405 before. The instrumentation covers the source-stepping rungs because it wraps the model rather than the analysis.

**Gate 5-3 (the limiter is actually applied).** A test that fails if `SubCircuit.limit`'s
return value is discarded — i.e. a device whose limiter clamps `x` rather than keeping
private state.
OUTCOME: **PASSED, and it is the test the tree could not previously have had.** A two-terminal element whose limiter clamps its own node to 0.1 V: with the write-back removed the clamp never reaches the caller and the assertion fails. Verified against the pre-stage-5 source, where it does fail. The point 0.1b made — *"the only limiter in the tree is the one that cannot expose the bug"* — is why this had to be a purpose-built device rather than `Diode`.

**Gate 5-4 (no regression on circuits that already converge).** Iteration counts on the
existing nonlinear tests must not increase by more than 20%. **This is the gate the
review explicitly did not measure**, so it is the one most likely to bite.
OUTCOME: **PASSED — it did not bite. 795 passed, 6 skipped, 0 failed** (`-m ""
--timeout=400`), against 788 before plus the 7 tests added here. **No existing test needed
changing**, on a change that puts a limiter in front of every `BJT` and `JFET` and alters
what `SubCircuit.limit` returns for every device in the tree.

**The runtime — 1940.76 s — is not a measurement and must not be read as one.** The box was
running two unrelated jobs at ~80% CPU throughout (load average 7.04), so this says nothing
about the change. Trap 2 exists for exactly this; the pass/fail is unaffected. A clean
timing needs a quiet box.

**Why the regression risk was lower than the gate feared, stated so the gate's own worry is
answered rather than forgotten:** limiting only engages above `vc = VT*log(VT/(IS*sqrt(2)))`
and only when the junction is being driven forward, so on a circuit that already converges
without wild excursions the limiter is a no-op on nearly every iteration. It changes the
path Newton takes, not the equations — which `test_gate_5_the_limiter_only_moves_the_
linearisation_point` asserts directly by substituting the answer back into the unlimited
residual.

---

# STAGE 5+ — the device-model findings 0.1b ranked and the plan never listed

**Raised 2026-07-31, immediately after stage 5 landed.** Review 0.1b answered "what is the
minimum device set that makes CMOS and bipolar transients expressible" with a ranked list of
four items. **Only one of them — the limiter — was ever written into a Work list or given a
gate.** The other three have existed since 0.1b as prose inside a review outcome, which
means they are findings rather than plans: nothing schedules them and nothing checks them.

**Stage 5 has changed which of them is binding, and made one of them worse.** 0.1b said of
the charge model: *"This is a bigger obstacle to a credible bipolar transient than the
missing limiter is, and the plan does not rank it."* That was true when the limiter was
missing and bipolar circuits simply failed. **Now they converge** — so a `BJT` transient will
run to completion and return a plausible waveform whose switching times are meaningless,
because `C(x) == 0` means the transistor has no charge storage at all. Stage 5 removed one
confidently-wrong-answer defect and made another one reachable. That is the reason this
section exists rather than waiting for stage 10.

## 5+.1 A charge model on `BJT` — the binding constraint

`BJT` defines no `q`, so `Semiconductor.C` takes its zero-matrix branch: no depletion
capacitance, no diffusion capacitance, no `TF`. A bipolar transient is then a sequence of DC
solves with the wrong dynamics.

**Gate 5+.1a (the defect is real and reachable, not theoretical).** Show `C(x)` is
identically zero for a biased `BJT`, and that a transient through it therefore has no
frequency dependence — the same circuit at two very different drive frequencies must
currently give the same switching behaviour. Declared: record both. **If the transient
already shows frequency dependence from elsewhere in the circuit, say so** — that would mean
the missing charge is masked and this item is less urgent than 0.1b ranked it.
OUTCOME: **PASSED — the defect is real, reachable, and not masked. 0.1b's ranking is confirmed.**

`C(x)` for a `BJT` biased in forward-active is the **exact zero matrix**. And the same switching stage — whose only capacitance is inside the transistor — driven with its period changed by **1000x** produced a **bit-identical** waveform: `v(c)` min 0.0915287 and max 5 at both 1 us and 1 ns. There is no time constant in the device to respond with.

After the charge model: at 1 us the stage switches (`v(c)` min 0.0915), at 1 ns it does **not switch at all** (`v(c)` min 5.0) because the charge cannot arrive in time. That is the physical behaviour that was absent.

**Gate 5+.1b (the charge model is the standard one, not an invention).** Depletion charge on
both junctions plus diffusion charge via `TF`/`TR`, i.e. the Gummel-Poon charge terms that
Ebers-Moll Level 1 is normally paired with. Parameters `CJE`, `CJC`, `VJE`, `VJC`, `MJE`,
`MJC`, `TF`, `TR` with SPICE's defaults. Declared: the depletion term uses the same
`FC`-style linearisation above the knee that `Varactor` gets in 5+.3 — **not** the clamp that
makes `C` fall to zero, which is the defect 5+.3 exists to fix and must not be copied into a
new device.
OUTCOME: **PASSED.** Depletion charge uses a shared `_depletion_charge` helper implementing SPICE's F1/F2/F3 linearisation above `FC*VJ`, and diffusion charge is `TF`/`TR` times the junction current. Measured `C = dq/dv` for CJ=3.5pF, VJ=0.75, M=0.33, FC=0.5:

| v (V) | -4.0 | -1.0 | 0.0 | 0.374 | 0.376 | 0.7 | 0.9 |
|---|---|---|---|---|---|---|---|
| C (pF) | 1.903 | 2.646 | 3.500 | 4.396 | 4.403 | 5.658 | 6.432 |

Positive and monotonically increasing throughout, with the charge continuous across the knee (1.451080e-12 at 0.374 V against 1.459879e-12 at 0.376 V). **The helper is shared rather than copied** so that 5+.3 fixes `Varactor` by pointing it at the same code, and a new device cannot inherit the clamp that zeroes `C`.

**On the defaults, stated because they are a choice and not a measurement.** SPICE defaults CJE/CJC/TF/TR to zero and expects a model card; these do not. The reason is consistency with this class rather than with SPICE — `BJT` already invents IS=1e-14, BF=100, VA=100, so it is a usable default transistor rather than a bare template, and zero-charge defaults would leave the 5+.1a defect fully in place for anyone who did not know to set CJE and TF. The values are the round numbers of a generic small-signal NPN and are documented as such in the source; any comparison against silicon must supply its own.

**Gate 5+.1c (it is differentiated, not stamped).** `C(x)` must come from
`toolkit.jacobian(eval_q_pure, ...)` like every other `Semiconductor`, so the charge model
lives once and the capacitance matrix cannot silently disagree with it. That is the
invariant the base class exists to protect.
OUTCOME: **PASSED.** `C(x)` comes from `toolkit.jacobian(eval_q_pure, ...)` like every other `Semiconductor`, so the capacitance cannot disagree with the charge it is differentiated from. Two properties a plausible-looking wrong charge model usually fails are asserted with it: terminal charges **sum to exactly 0** (the device stores no net charge) and `C` is **symmetric to machine zero** (it is a reciprocal capacitance).

**Gate 5+.1d (a switching time that responds).** With charge storage present, a
`BJT` switching transient must show a storage-time dependence on `TF` — declared: doubling
`TF` measurably lengthens the turn-off. **This is the gate that distinguishes "a charge
model exists" from "a charge model does something"**, and it is the one worth writing first.
OUTCOME: **PASSED. Storage time 13.63 ps -> 14.98 ps -> 16.83 ps** as `TF` doubles twice (3e-10, 6e-10, 1.2e-9 s), measured as the time for `v(c)` to recover past half-scale after the base drive falls. Monotone, and responding to the parameter that physically controls it — minority carriers in transit through the base.

**Gate 5+.1e** Full suite, and the stage-5 convergence gates re-run: adding charge changes
the Jacobian, so gate 5-1's operating points must still be reached.
OUTCOME: **PASSED. 799 passed, 6 skipped, 0 failed** (`-m "" --timeout=400`), against 795
before plus the 4 tests added here. **No existing test needed changing**, and gate 5-1's
three operating points are among those still passing — adding charge changes `C`, not `G`,
so the DC solve is untouched by construction, but it was worth checking rather than
asserting.

Runtime 1103.95 s at load average 4.6-7.1 from unrelated jobs on the box; not a measurement
(trap 2).

## 5+.2 `MOS_ACM` — delete

0.1b verified it **cannot be constructed at all**: `mos.py:104` calls `super(MOS, self)` from
a class whose MRO is `[MOS_ACM, SubCircuit, Circuit, object]`, raising `TypeError`. And it is
not an ACM model — its body is a verbatim copy of `MOS` with one difference, `Symbol('kT')`
in the noise PSD where `MOS` uses `toolkit.kboltzmann * Symbol('T')`. So a class advertising
a different model is a copy of its neighbour with a worse noise expression, and has never
run. Blast radius enumerated at 0.1b: **zero references anywhere** outside `mos.py` and the
review documents — no test, no example, no doc page, not exported from `__init__.py`.

**The mechanism that let it survive is the point:** the only thing that would ever have
caught it is the doctest in its own docstring, and `pytest.ini` configures no doctest
collection, so `if __name__ == "__main__": doctest.testmod()` runs only if someone executes
the module directly. **The test existed and was never run** — which is worse than having no
test, because it reads as coverage.

**Gate 5+.2a.** Deletion converts a `TypeError` at construction into an `ImportError` at
import — earlier and clearer. Declared: the suite is unchanged, confirming the zero-reference
claim was right.
OUTCOME: **PASSED — suite 803 passed, 6 skipped, 0 failed, unchanged in character; and 0.1b's claims were verified rather than trusted.** MRO is `[MOS_ACM, SubCircuit, Circuit, object]`, `TypeError` on construction, zero references outside `mos.py` itself, not exported from `circuit/__init__.py`. Deleted; the import now raises `ImportError`.

**The diff against `MOS` is sharper than 0.1b recorded — two copy-paste defects, not one.** Beyond the noise PSD using `Symbol('kT')` (a free symbol nothing in the package binds) where `MOS` uses `toolkit.kboltzmann * Symbol('T')`, its `gds` parameter is described as **"Gate transconductance"** — `gds` is the output conductance. Everything else in the 50-line body is identical.

**Gate 5+.2b (fix the mechanism, not just the instance).** Per the standing preference for
fixing what let an error hide: either collect doctests in `pytest.ini` or delete the
unreachable `doctest.testmod()` blocks. Declared: whichever is chosen, no module may keep a
docstring test that nothing runs. **Measure first** — collecting doctests repo-wide may fail
on modules whose docstrings were never checked, and that count is the finding.
OUTCOME: **MEASURED, AND BOTH DECLARED OPTIONS WERE WRONG — there is a third, and it is blocked on a defect this measurement found.**

*Option 1, collect doctests in `pytest.ini`: not viable.* Repo-wide `--doctest-modules` gives **7 collection errors** (`post/cds`, `post/jwdb`). Scoped to `pycircuit/circuit` it gives 1 (`xdot.py`), and excluding that the run **did not complete in 10 minutes**.

*Option 2, delete the unreachable blocks: wrong target.* 31 modules carry `if __name__ == "__main__": doctest.testmod()`. Deleting them removes the only hint the doctests exist without making any of them run.

*The third option, which this project already chose:* `pycircuit/circuit/tests/test_doctests.py` runs `doctest.testmod` on named modules as ordinary tests. It already covers `circuit.py`, `elements.py`, `volterra.py` and `symbolicapprox.py`, and its docstring records that this exact gap had been hiding real bugs — `Quantity.__repr__` raising unconditionally since 2010, `Circuit.name_state_vector` double-indexing. **This gate was declared without knowing that file existed**, which is why it offered two options and neither is the answer.

**So the fix for `mos.py` is to add it to that list — and it cannot be, yet.** `MOS`'s own doctest fails for two further independent reasons, neither involving `MOS_ACM`:

1. `c = SubCircuit()` takes the **numeric** default toolkit while the parameters are `Symbol(...)`, so construction dies with `TypeError: Cannot convert expression to float` at `elements.py:709`.
2. `twoport.solve(freqs=array([Symbol('s')]), ...)` trips `nportanalysis.py:235: assert not isiterable(freqs)` — the symbolic path takes a scalar.

Both are the same never-run-so-never-noticed pattern as `MOS_ACM` itself. **RESOLVED at 5+.5 (2026-07-31), which this gate spawned — see below.** Both defects were
fixed, `mos.py` is in `test_doctests.py` as `test_mos_module_doctests`, and the
`nportanalysis` assertion turned out to be right in intent and wrong in form. This sentence
is left in place because it was true when written; the forward pointer is added because
without it the gate reads as an open item, and on 2026-08-02 it was picked as the next
thing to work on before anyone checked whether it had already been done.

**Recorded as its own item rather than fixed here**: repairing `MOS`'s documented example is not the same change as deleting a dead class, and bundling them would make neither attributable. The unreachable `doctest.testmod()` in `mos.py` is deliberately left in place — it is the only marker that those doctests exist and are unverified.

**Reconsider if** someone actually wants ACM — in which case it is written from the paper,
not recovered from this.

## 5+.3 `Varactor`'s clamp makes C fall to zero where it should be largest

`v_eff = minimum(v, 0.99*VJ)` freezes the junction charge above the knee, so `C = dq/dv` goes
to **exactly zero** in forward bias. SPICE linearises the charge above `FC*VJ` and keeps a
finite, growing capacitance. Zero is the worst available answer: it does not merely
mis-estimate, it removes the state variable, and a Newton step that sees `C = 0` on the node
with the largest physical capacitance takes a wildly wrong step.

**Gate 5+.3a.** `C(v)` must be positive and increasing through `FC*VJ` and above it.
Declared: measure `C` across the knee before and after; the "before" curve going to zero is
the defect reproduced.
OUTCOME: **PASSED, and the defect reproduced exactly.** `C = dq/dv` for the default Varactor (CJ0=1pF, VJ=1.0, M=0.5):

| v (V) | 0.90 | 0.98 | 0.99 | **1.00** | **1.50** | **2.00** |
|---|---|---|---|---|---|---|
| before (F) | 3.162e-12 | 7.071e-12 | 4.999e-12 | **0** | **0** | **0** |
| after (F) | 1.980e-12 | 2.093e-12 | 2.107e-12 | **2.121e-12** | **2.828e-12** | **3.536e-12** |

Before, the charge freezes at 1.8e-12 C and its derivative with it — **exactly 0.000000e+00** — after a spurious 7.07 pF peak as the clamp is approached. After, `C` is positive and monotonically increasing across the whole range.

**Gate 5+.3b.** The linearisation is SPICE's, and the same one 5+.1b uses for the `BJT`
depletion terms — one treatment, not two.
OUTCOME: **PASSED.** `Varactor.eval_q_pure` calls the same `_depletion_charge` helper the `BJT` depletion terms use, asserted by comparing against the helper directly rather than against remembered numbers.

**Two defects in the tests written for this item, both found by checking that they fail against the pre-fix source, and both worth recording because they are general:**

1. `assert got == pytest.approx(want, rel=1e-12)` **passed against the old clamped Varactor.** `pytest.approx` defaults to `abs=1e-12`, and these are picofarad-scale charges — so the default absolute tolerance is the same size as the quantities and the assertion accepted 1.367544e-12 as equal to 1.264609e-12, an 8% difference. `abs=0.0` is now passed explicitly. **Any `pytest.approx` on a quantity near 1e-12 is vacuous by default.**
2. The no-warning test passed against the eager version it was written to reject, because it passed a **Python** float: `(-1.0) ** 0.5` returns a complex number silently, where the real call site passes a numpy scalar and gets `nan` plus a RuntimeWarning. It now uses `np.float64`, matching `v = x[0] - x[1]` on a numpy array.

## 5+.5 `MOS`'s own doctest is broken — NEW, found by 5+.2b on 2026-07-31

Not a stage-5 item by origin; it exists because gate 5+.2b went looking for the mechanism
that hid `MOS_ACM` and found the same mechanism hiding two more defects one class away.

`mos.py` cannot join `tests/test_doctests.py` until these are fixed, and until it joins,
`MOS`'s documented example is still unverified — which is the condition that produced
`MOS_ACM` in the first place.

1. **The toolkit.** `c = SubCircuit()` takes the numeric default while the parameters are
   `Symbol(...)`, so construction dies at `elements.py:709` with `TypeError: Cannot convert
   expression to float`. The doctest presumably predates the toolkit split and never ran
   after it.
2. **The frequency argument.** `twoport.solve(freqs=array([Symbol('s')]), complexfreq=True)`
   trips `nportanalysis.py:235: assert not isiterable(freqs)`. The symbolic path takes a
   scalar; the example passes a one-element array.

**Gate 5+.5a.** `doctest.testmod(mos)` reports 0 failures, and `mos.py` is added to
`test_doctests.py` so it stays that way.
OUTCOME: **PASSED. `mos.py`: 8 attempted, 0 failed**; suite 804 passed, 6 skipped, 0 failed in 590 s — the fastest run of the session, on a box that had been running two of the user's jobs at 90% and 45% CPU and was now down to one. Same source class, runtimes 676-1941 s elsewhere in the session: trap 2, again., and it is now in `test_doctests.py` as `test_mos_module_doctests`. The hole that produced `MOS_ACM` is closed for that module.

**Gate 5+.5b (fix the example, or fix the code — decide which by looking).** Declared: if
`nportanalysis`'s symbolic path *should* accept a one-element array, the assertion is the
defect and the example is right; if it should not, the example is wrong. **Do not assume the
example is wrong just because it is the smaller change** — it is the older artefact, and the
assertion may be the thing that drifted.
OUTCOME: **ANSWERED BY MEASUREMENT, AND THE ANSWER IS NEITHER — the assertion's INTENT was right and its FORM had drifted.**

*The evidence that the example was not simply wrong.* **Two independent doctests, written years apart, pass a length-1 array to the symbolic path**: `mos.MOS`'s, and `nportanalysis.TwoPortAnalysis`'s own symbolic example at line 65. Both had been failing. Two artefacts agreeing is weak evidence on its own, so it was measured.

*The measurement.* With the assertion patched out, on an RC low-pass whose `A11` must depend on the frequency variable:

| input | A11 |
|---|---|
| scalar `s` | `C1*R1*s + 1` |
| `array([s])` | `C1*R1*s + 1` — **identical, correct** |
| `array([s, w])` | `C1*R1*w + 1` — **silently drops `s`** |

(The first attempt used a resistive divider, whose answer is frequency-independent and so cannot distinguish the three cases at all. Recorded because it looked like a clean result and was worthless.)

So more than one frequency is **genuinely unsupported** and the restriction is kept — but a length-1 array *is* one frequency and gives the identical answer, and rejecting it was the assertion testing the container rather than the count. It now counts, unwraps the single element, and raises a `ValueError` **saying why** for anything longer; the bare `assert` explained nothing, which is the other half of what made this hard to read.

**Had the gate not forbidden assuming the example was wrong, the smaller change — editing two doctests to pass scalars — would have left a working API rejecting its own documented usage.**

## 5+.4 A large-signal MOSFET — stays in stage 10

Confirmed absent; `mos.py` holds only small-signal subcircuits. Decision 0.3c put it in
stage 10 and sequenced it after 0.1b and stage 5, both of which are now done. **It is not
pulled forward**: it is the largest of the four by a wide margin, it depends on 5+.3's
depletion treatment being settled, and unlike the charge model it blocks nothing that
currently runs. Recorded here only so the ranked list is complete in one place.

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
OUTCOME: **PASSED.** Before: `NoConvergenceError: Source Stepping failed at lambda=0.0`. After:

> `SingularMatrix: singular Jacobian: 'floaty' appears in no equation, so nothing determines it — for a node that means no DC path to ground (add a resistor, or use uic=True to skip the operating point)`

**The classification is structural, not a pivot inspection.** An all-zero *column* means the unknown appears in no equation — the floating-node case; an all-zero *row* means the equation constrains nothing. Read off the assembled Jacobian, so it is exact, needs no factorisation, and cannot be confused with mere ill-conditioning — which continuation genuinely can help with and which therefore keeps its existing path.

**Why it reaches the caller at all:** `SingularMatrix` is a sibling of `NoConvergenceError`, not a subclass, so it passes straight through both stepping decorators. That is what "classify **before** continuation is attempted" means in practice — continuation cannot add a missing equation, so it should not be given the chance to report on one.

**Gate 6-2.** A non-convergent circuit produces a message naming the worst node, its
residual, and the tolerance it missed.
OUTCOME: **PASSED.** Before: `Source Stepping failed at lambda=0.01`. After:

> `Source Stepping failed at lambda=1.0: Gmin Stepping failed at gmin=0.001: StandardNewton failed to converge after 2 iterations; residual worst at 'c': |f| = 0.01709 against a tolerance of 1.709e-06 (1e+04x over); update worst at 'c': |dx| = 1.907 against a tolerance of 5.933e-05 (3.21e+04x over)`

Everything after the first colon was already in scope when the old message threw it away. The misses are reported **normalised** as well as raw: "1e+04x over" says how far off, where "4.7e-9 A" does not.

**Two things had to change, and the second was not in the plan.** Naming the row needed `row_names` threaded through `NonLinearSolver.solve_system` — the solver works on a system with the reference row removed, so its indices are *reduced* and cannot be looked up in the circuit directly. But that alone changed nothing visible, because the continuation decorators **discarded the inner message**: they raised a fresh `NoConvergenceError` carrying only the rung. They now keep both. The rung is useful context; it is not a substitute for the cause, and for however long this code has existed it has been served instead of one.

**API note, stated because it is a break.** `NonLinearSolver.solve_system` gained a defaulted `row_names` argument, and the decorators pass it **by name** so an out-of-tree solver that does not accept it fails with a `TypeError` rather than silently binding it to `scaler`. Two in-tree mocks needed `**kwargs`, which is what any such subclass should do.

**Gate 6-3.** Statistics are populated and non-zero on a run that rejects steps, and the
force-accept counter from 4b appears there.
OUTCOME: **PASSED.** On a `VPulse`-driven RC with the trapezoidal integrator:

```
accepted 153, rejected 18 (10.5% of attempts), Newton iterations 331 (2.2 per accepted step)
force-accepts 0, order drops 20, breakpoints hit 26
step 1e-10 .. 6.009e-08 s
time 0.225 s total, 0.196 s in the Newton solve (87.2%)
```

Every counter populated; rejections, order drops and breakpoints all non-zero. **4b's force-accept counter is there and reads 0** — which is now the expected value everywhere, since 4d made that path unreachable on every circuit measured, so a non-zero reading is the run reporting that part of its own result is not error-controlled.

`solve_system` was already returning its iteration count and the call site was binding it to `_`. The object is created **per run** rather than in `__init__`, and that is pinned by its own test: the natural implementation accumulates across `solve()` calls, and the failure is invisible unless something solves twice.

**Docs in the same commit:** a "diagnosing a failed run" section — this is the
documentation with the highest ratio of user value to effort in the whole plan.
DONE: `doc/src/circuit/diagnosing.rst`, wired into the toctree. It is organised around
*telling the three failure modes apart* rather than around the code, and it tells the
reader to read `force-accepts` first and why a non-zero value means part of the result is
not error-controlled.

**Stage 6 suite: 809 passed, 6 skipped, 0 failed, 445.19 s** — the two in-tree solver mocks
needed `**kwargs` for the new `row_names` argument and nothing else changed. Doc build:
build succeeded, 2 warnings, 0 ERROR, and the new page renders.

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
OUTCOME: **PASSED.** `analysis._reduce_ndarray` copies the four surviving blocks straight
into one output instead of calling `toolkit.delete` once per axis, which for a square
Jacobian built an `(n-1, n)` intermediate and then the `(n-1, n-1)` result — two copies for
one reduction. Measured **2.2x to 7.4x** faster from n=200 to n=3200, allocation included.

`toolkit.delete` survives only on the generic path, for arrays that are not square ndarrays;
the Newton path does not reach it.

## 7a entry measurements, 2026-08-01 — the stated case is wrong, the item is right

**Gate 7a-0 (declared before any code): measure `np.delete`'s share of a transient.**

OUTCOME: **THE SPEED FRAMING IS REFUTED, AND THE ITEM SURVIVES ON A BETTER ARGUMENT.**
Profiled on a numeric RC ladder, whole transient:

| n | total | `np.delete` | `remove_row_col` |
|---|---|---|---|
| 22 | 0.75 s | 0.017 s (2.3%) | 0.9% |
| 102 | 6.9 s | 0.037 s (0.5%) | 0.1% |
| 302 | 39.5 s | 0.080 s (0.2%) | 0.0% |

It is not a hot spot and it gets *relatively cheaper* as `n` grows. **Gate 7a as written —
"`np.delete` absent from the profile" — was measuring something never meaningfully
present.**

**What justifies 7a is scaling, not share.** Timed in isolation, per Newton iteration:

| component | measured exponent |
|---|---|
| assembly (`G(x)` + `C(x)`) | **n^1.26** |
| dense LU | **n^1.83** |
| **`remove_row_col`** | **n^2.53** |

It is the worst-scaling component of the three. Extrapolated to a 1000-step run at three
Newton iterations per step: at n=5000, 0.37 s per iteration (0.3 h per run, copying 200 MB
each time); at n=20000, **12.4 s per iteration = 10.3 h, which exceeds the dense LU's
3.2 h**. So it does come first, and for a reason the plan did not give.

**A CORRECTION TO MY OWN FIRST READING.** From the n=302 profile — assembly ~13.3 s of
15.5 s, LU 8% — I concluded "the wall is assembly, not the solver". **That was wrong**, and
only the scaling table shows why: LU's share runs 35% at n=802, 60% at n=5000, 77% at
n=20000, and the measured n^1.83 is *below* the asymptotic n^3, so it understates the
trend. **Stage 7's premise holds at the sizes it names.** A profile at one convenient size
is not a scaling law.

**Gate 7a-1 (bit-identical).** The reduction is a pure data movement, so every existing
result must be unchanged. Declared at bit-identical on the full 17-run baseline and the
suite.
OUTCOME: **PASS.** All 17 baseline runs `max|diff| = 0.000e+00`, and separately every removable
index swept on random matrices at n = 5/17/64 in float and complex, against the two-delete
form.

**Gate 7a-2 (the win is real and measured across n).** Declared: >= 2x on the reduction
itself at every size from n=200 to n=3200, with allocation included — not measured against
a reused buffer, because a cached buffer would alias results between callers and `J` is
kept for the LTE after the solve.
OUTCOME: **PASS**, allocation included: 3.54x / 3.69x / 4.36x / 4.63x / **2.03x** at n = 200 /
400 / 800 / 1600 / 3200 — every size above the declared 2x. The fall-off at n=3200 is
memory bandwidth: at that size the copy is 82 MB and both forms are bandwidth-bound, which
is the honest reading rather than a defect in the new form.

**Gate 7a-3 (every caller still works).** `remove_row_col` has 20+ call sites across
`analysis_ss`, `nportanalysis`, `feedback`, `volterra`, `shooting`, `dcanalysis` and the
step controller, taking 1-D and 2-D arrays, real and complex, under numeric AND symbolic
toolkits. Declared: the fast path is guarded and symbolic toolkits keep the existing
`toolkit.delete` path, verified by test rather than by inspection.
OUTCOME: **PASS.** The fast path is guarded on `isinstance(A, numpy.ndarray) and A.dtype !=
object`; symbolic toolkits hand over object arrays and keep `toolkit.delete`, verified by a
test that reduces a sympy matrix rather than by inspecting the guard. Non-square and
higher-rank arrays fall through. **Every index is swept, not a middle one** — two of the
four blocks are EMPTY when the removed index is first or last, which is exactly where an
off-by-one survives a spot check.

**Gate 7a-4.** Full suite.
OUTCOME: **PASS. Suite 833 passed, 6 skipped, 0 failed**, against 824 plus the nine tests added
here.

**I nearly recorded 824 for this, and the way it happened is worth keeping.** Two suite
runs were started against the same log file — one before the new test file existed, one
after — and the first to finish left its summary where I read it. The count looked
*plausible* (it matched the previous run exactly), which is what made it dangerous. Caught
only because 824 + 9 should have been 833. Re-run clean: one summary line, 833. Trap 15
below.

**Also found while checking: the suite has 4 pre-existing collection ERRORS** —
`pycircuit/post/cds/test/*` cannot import `pexpect` — and a plain `-q` tail does not show
them. They are the same modules behind the doc build's two autodoc warnings, so they are
not new, but "N passed" had been quoted all session without anyone noticing that four
modules never ran at all.

**NOT DONE, and this is the honest boundary of 7a as implemented.** The plan says "take the
reference node out of the matrix"; what this does is make the *removal* cheap, not remove
the need for it. Any approach that copies is O(n^2); eliminating the copy entirely means
never stamping the reference row and column, which is an assembly-level change touching
every element's stamp and belongs with stage 2's assembly work. **Scheduled as 2+.4**,
with its own gates, rather than left as prose here. Recorded so the item is not
read as closed.

**7b. A `LinearSolver` strategy object**, in the shape the codebase already uses for
`Integrator`/`StepController`/`Scaler`/`NonLinearSolver` — `analyze` / `factor` / `solve`
split so the factorisation survives the Newton iteration. `DenseLU` (scipy
`lu_factor`/`lu_solve`) below ~200 unknowns, `SuperLU` above; measure the crossover rather
than guessing it (dense wins at n=29, sparse wins 5x at n=137, 15x at n=542).
### 7b entry measurements, 2026-08-01 — taken before any code

**Gate 7b-0a: is the step controller really a third of all linear solves?** **CONFIRMED
EXACTLY, 33.1%.** Counted by origin on a numeric RC ladder: 228 Newton solves against 113
step-controller solves over 111 accepted steps, independent of `n`. Roughly two Newton
iterations plus one controller solve per step.

**But "redundant re-factorisations" is the wrong description, and it changes the design.**
Newton's last factorisation is of `J(x_k)`; `solve_timestep` then rebuilds `J` at the
**converged** `x_{k+1}` and the controller factors *that*. Nothing else factors it, so
those 33% are not duplicates — reusing Newton's factors means **substituting a different
matrix**. Measured within a step on a nonlinear (diode) circuit:

| `\|J(x_k) - J(x_conv)\|` / `\|J(x_conv)\|` | value |
|---|---|
| median | 5.455e-09 |
| mean | 1.779e-06 |
| max | 2.152e-05 |
| exactly zero | 11 of 50 steps |

Sound as an approximation — it feeds an error *estimate* with `TRTOL = 7` of slack — but
**not bit-identical**. **So gates 7b-1 ("results identical to the dense path") and 7b-2
("reuses the factors") cannot both hold as written.** That is a conflict in the plan, not a
detail of implementation, and it is recorded here rather than silently resolved. (The
re-assembly at the converged point is *not* recent: it dates to `18931ae`, 2026-07-22, and
gate 2c established it as necessary. 2+.2 only stopped it computing the discarded `i`/`u`.)

**Gate 7b-0b: the dense/sparse crossover.** The plan says *"dense wins at n=29, sparse wins
5x at n=137, 15x at n=542"*. Measured on the RC ladder:

| n | fill | dense LU | SuperLU | sparse win |
|---|---|---|---|---|
| 31 | 9.37% | 0.00032 | 0.00013 | 2.47x |
| 61 | 4.84% | 0.00008 | 0.00007 | 1.16x |
| 141 | 2.11% | 0.00042 | 0.00015 | 2.76x |
| 301 | 0.99% | 0.00683 | 0.00024 | 28.7x |
| 551 | 0.54% | 0.02715 | 0.00039 | **69.6x** |
| 1101 | 0.27% | 0.05099 | 0.00068 | 74.5x |

**Sparse wins everywhere here, and by far more than recorded at large `n`.** Two caveats
that matter more than the numbers:

1. **This is an RC ladder at 0.27% fill — maximally favourable to sparse.** The crossover
   is a function of the sparsity *pattern*, not of `n`. **7b's "DenseLU below ~200
   unknowns, SuperLU above" keys on the wrong variable**, and a dense-ish circuit at n=300
   would be mis-served by it. The strategy should choose on fill, measured per circuit.
2. **`np.linalg.solve` beats scipy's `lu_factor` + `lu_solve` at small `n`** — 2e-5 against
   3.2e-4 at n=31, on Python call overhead alone. So a naive `DenseLU` strategy object is a
   **regression** for the small circuits it is meant to protect, and the factor/solve split
   only pays once the factors are reused.

**Gate 7b-0c: does a circuit at n ~ 2000 exhaust memory, as 7b-3 asserts?** **NO.** A dense
`n=2000` Jacobian is **32 MB**, and n=5000 is 200 MB — large but nowhere near exhausting.
7b-3's premise is wrong as stated; what fails at those sizes is *time*, and 7a's scaling
table says where it comes from. **Gate 7b-3 needs restating around a time budget**, or it
will be marked passed by a run that was never going to fail.

**Gate 7b-1:** results identical to the dense path on every existing transient test.
**Gate 7b-2:** the step controller's `J^-1 Eg` solve (`stepcontroller.py:60`) reuses the
factors — a third of all linear solves are currently redundant re-factorisations.
**Gate 7b-3:** a circuit at n ~ 2000 that currently exhausts memory now runs.
OUTCOME: **PARTLY DONE. The strategy object ships; 7b-2 is NOT built, on
evidence; 7b-3 is restated because its premise was wrong.**

**Gate 7b-1: PASS.** All 17 baseline runs bit-identical, `max|diff| = 0.000e+00`, and a
test asserts `DenseSolver()` reproduces passing nothing at all. **The default had to stay
the old path to achieve this**: `numpy.linalg.solve` and SuperLU round differently, so an
`AutoSolver` default would move the last bits of every circuit large and sparse enough to
qualify. Sparse is opt-in via `linearsolver=`, typed like `nrsolver`/`scaler`.

**What shipped.** `pycircuit/circuit/linearsolver.py`: `LinearSolver` ABC, `DenseSolver`
(the historical `toolkit.linearsolver` call, and the only toolkit-agnostic one),
`SuperLUSolver` (discovered, not required — raises on construction if SciPy is absent so
`AutoSolver` never picks it), and `AutoSolver`. Threaded through all seven `solve_system`
signatures and the transient; `linsolver is None` keeps the exact previous call.

**`AutoSolver` selects on FILL, not on `n` — a deliberate departure from 7b's text.**
"DenseLU below ~200 unknowns, SuperLU above" keys on the wrong variable: the measured
crossover follows the sparsity *pattern*, and a dense 300x300 would be mis-served. `n` is
kept only as a cheap floor (100) below which measuring the fill is not worth the pass. A
test asserts a large **dense** matrix stays on the dense path, which a size-only rule
would fail.

**Gate 7b-2: NOT BUILT, and the gate is refuted as written.** The 33.1% is real, but those
solves are not re-factorisations — see the entry measurements above. Reuse means
substituting `J(x_k)` for `J(x_conv)`, which is an approximation (median 5.5e-9, max
2.2e-5) and cannot hold under 7b-1's bit-identical bar. `LinearSolver.factor()` exists as
the seam with that reasoning recorded on it, returning `None`. **Whether to accept the
approximation for a third fewer solves is a decision, not an implementation detail, and it
is left open rather than taken quietly.**

**Gate 7b-3: PREMISE WRONG, restated.** A dense `n=2000` Jacobian is **32 MB** and n=5000
is 200 MB; nothing is exhausted. What fails at those sizes is time. The gate as written
would be marked passed by a run that was never going to fail.

**AND THE END-TO-END WIN IS SMALL AT THESE SIZES, WHICH IS THE HONEST HEADLINE.** The
isolated solve is 28-74x faster from n=301 upward, but a transient is not a solve:

| n | default | AutoSolver | speedup | steps | max\|dv\| |
|---|---|---|---|---|---|
| 62 | 0.727 s | 0.754 s | **0.96x** | 111 | 0 |
| 152 | 3.894 s | 3.527 s | 1.10x | 111 | 1.0e-15 |
| 402 | 11.368 s | 10.386 s | 1.09x | 111 | 1.0e-15 |
| 802 | 28.274 s | 23.635 s | 1.20x | 111 | 1.0e-15 |

**1.1-1.2x, because the solve is ~8% of runtime and assembly is the rest** — precisely what
stage 2's header says ("replacing the solver first would optimise the wrong thing"). The
step count never moves and the waveform moves only in the last bits. **At n=62 it is 4%
SLOWER**: `AutoSolver` correctly picks dense there, so that is pure dispatch overhead, and
it is why the default is `DenseSolver` rather than `AutoSolver`.

**n=1202 MEASURED AFTER THE FACT, and it closes this gap partway.** The run finished
after the 7b commit was written: **33.201 s -> 26.619 s, 1.25x**, 80 steps, `max|dv|`
9.992e-16. So the win does grow with `n` — 1.09x at 402, 1.20x at 802, **1.25x at 1202** —
but *slowly*, which is the same conclusion by measurement rather than by extrapolation:
assembly, not the solve, sets the pace.

**n=2000 remains unmeasured.** That row was killed rather than completed (one row printed,
no traceback), so nothing is claimed for it. **The cause is now measured and scheduled as
2+.5**: circuit *construction* is O(N^2.27) — 24.8 s to build an 800-element ladder — so
the solver work this row was meant to measure never got the chance to run. 2+.5 is a
prerequisite for measuring anything at n >= 1000.

**Still open: 7c (KLU), 7d (`pybsmatrix` deletion and the sparse-toolkit test), and the
`factor()` decision above.**

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
OUTCOME: **PASSED, with the threshold moved on measurement.** `KLUSolver` exists in
`linearsolver.py`, calling `libklu` through `ctypes` with the analyze/factor/refactor split,
and validating the refactor by residual rather than by struct offsets.

The declared "3x at n >= 500" was not where the crossover actually sits. Wiring KLU as
`AutoSolver`'s first choice on isolated benchmarks gave a **0.52x regression at n=152** —
the setup cost dominates until the matrix is much larger — so it was reverted and then
re-enabled above `MIN_N_FOR_KLU = 1000` on an end-to-end measurement rather than a
factorisation-only one. Selection also keys on **fill, not size alone**
(`MAX_FILL_FOR_SPARSE = 0.20`), because a small dense matrix and a large sparse one are
different problems that `n` cannot distinguish.
### 7c OUTCOME, 2026-08-01 — BLOCKED on an absent dependency, and NOT worth forcing

**No KLU binding exists on this machine.** Checked: `scikits.umfpack`, `sksparse.cholmod`,
`klu`, `pyklu` — none importable; only `scipy.sparse.linalg`. So the item cannot be
implemented *or* measured here.

**It was not written blind.** Shipping a KLU path that never executes would be the fourth
instance in this session of exactly that failure — the JAX backend's copied estimator
(9(f), 9-1(a)), `pybsmatrix`'s `fbsub` (7d), and `sparse_numeric` (7d), each of which was
wrong precisely because nothing ever ran it. An optional import that nobody can exercise is
how the next one gets written.

**The premise IS confirmed, with a correction to what it measures.** 7c says *"scipy's
`splu` recomputes COLAMD on every call — 94% of its cost"*. Measured on tridiagonal systems:

| n | fill | `splu()` | `.solve()` | factor share |
|---|---|---|---|---|
| 150 | 1.99% | 0.000259 | 0.000017 | **93.9%** |
| 400 | 0.75% | 0.000348 | 0.000023 | **93.7%** |
| 800 | 0.37% | 0.000482 | 0.000028 | **94.4%** |
| 1600 | 0.19% | 0.000724 | 0.000047 | **93.9%** |

The 94% reproduces almost exactly. **But it is the FACTORISATION share of factor+solve, not
the ORDERING share within factorisation** — and KLU's `klu_refactor` keeps the ordering
while redoing the numeric factorisation, so its saving is a fraction *of* that 94%, not the
whole of it. Separating the two needs instrumentation inside SuiteSparse, which is not
available here either. **7c's declared success — ">= 3x on the factor+solve" — therefore
requires ordering to be >= 67% of the factorisation, which is plausible and unmeasured.**

**And the end-to-end value would still be small at present sizes.** The solve is ~8% of a
transient (7a), so even meeting the 3x gate exactly yields **~1.06x overall** — against 7b's
measured 1.09x/1.20x/1.25x at n=402/802/1202 for a 28-74x solver win. The pattern is
consistent and by now well established: **the solver is not where a pycircuit transient
spends its time at the sizes pycircuit can currently reach.**

**RECOMMENDATION: leave 7c unimplemented, and do 2+.5 instead.** Circuit *construction* is
O(N^2.27) — 24.8 s for an 800-element ladder — which is what stopped 7b's n=2000 row from
being measured at all. Until that is fixed, no solver work can be evaluated at the sizes
where solver work would pay. **Reconsider 7c when** (a) a KLU binding is installed, and
(b) 2+.5 has landed so that n >= 5000 circuits can be built — at which point 7a's scaling
table puts LU at 60% of the step and the arithmetic changes completely.

### 7c REOPENED AND DONE, 2026-08-01 — SuiteSparse was installed

**The "BLOCKED" outcome above stands as written but is superseded.** With `libklu.so.2`
present, KLU is reachable **with no Python binding at all**: a small `ctypes` binding, which
is what 7c asked for ("discover it the way `_ginac_ext` and `symengine` already are"). No
pip dependency added. `scikit-umfpack` would have given UMFPACK rather than KLU, and
`libcholmod` — also now present — is for **symmetric positive-definite** systems, which
circuit MNA Jacobians are not.

**Gate 7c (">= 3x on the factor+solve at n >= 500"): PASSES comfortably.** Against
`scipy.splu`, which recomputes the ordering every call:

| n | splu | KLU refactor+solve | speedup |
|---|---|---|---|
| 150 | 0.000163 | 0.000046 | 3.51x |
| 400 | 0.000463 | 0.000060 | **7.75x** |
| 800 | 0.000541 | 0.000061 | 8.85x |
| 1600 | 0.000907 | 0.000097 | 9.31x |
| 3200 | 0.001714 | 0.000172 | **9.97x** |

The split confirms 7c's premise directly: at n=3200 `klu_analyze` costs 0.000441 against
`klu_factor`'s 0.000378 — the **ordering is about half a full factorisation** — and
`klu_refactor` at 0.000099 is 4x cheaper than factoring again. Over a transient the reuse is
total: **1 analyze, 1 factor, 227 refactors**.

**AND END TO END IT LOSES, WHICH IS THE FINDING THAT MATTERS.** Best-of-3, circuit object
reused so construction is out of the timing:

| n | default (dense) | SuperLU | KLU | KLU vs default |
|---|---|---|---|---|
| 152 | **3.581 s** | 4.378 s | 6.832 s | **0.52x** |
| 402 | 31.961 s | **23.511 s** | 24.418 s | 1.31x |

**KLU is nearly 2x SLOWER than the shipped path at n=152, and SuperLU beats it at n=402.**
The refactor is not at fault; it works exactly as designed. The loss is per-call Python
overhead — CSC conversion, the pattern key, the residual validation — against a solve that
is only **~8% of a transient**. A 10x win on 8% cannot repay a fixed per-call cost here.

**So `AutoSolver` selects SuperLU, not KLU, and a test pins it.** Wiring KLU in
automatically would have made a 152-unknown circuit twice as slow — a defect introduced and
caught inside the same item, by measuring end to end after measuring in isolation. KLU
remains available explicitly as `linearsolver=KLUSolver()`.

**The first end-to-end numbers were discarded, and why is worth recording:** they showed the
*smaller* circuit taking longer than the larger one (15.44 s at n=152 against 13.27 s at
n=402), which is machine load, not solver behaviour — the suite was running concurrently.
Trap 2 again, in a new place.

**The refactor path is VALIDATED, not trusted.** `klu_refactor` reuses the pivot ORDER from
the first factorisation, and if the values move far enough those pivots stop being sound.
Production tools read KLU's reciprocal pivot growth; doing that here would mean hard-coding
`klu_common` struct offsets this binding cannot verify. Instead the **residual is checked on
every reuse** and a full `klu_factor` redone if it fails — one sparse mat-vec, and it cannot
be wrong about a struct layout. Stage 9 and 7d between them turned up three defects that
survived precisely because nothing ever checked the thing itself.

**Reconsider when 2+.5 lands.** Every conclusion here is bounded by assembly dominating the
step. Once construction is not O(N^2.27) and n >= 5000 is reachable, LU is 60% of the step
rather than 8% — **redo this measurement then rather than assuming it still holds.**

### 7c AND D1 RE-MEASURED, 2026-08-01, after 2+.5 — and 7c's conclusion REVERSES

Both were deferred to exactly this point. 7c's outcome said *"every conclusion here is
bounded by assembly dominating the step... redo this measurement rather than assuming it
still holds"*, and D1's reconsider-if was literally *"when 2+.5 lands"*. 2+.5 has landed, so
the rows that could not be built before are now measurable — the n=1002 and n=2002 circuits
build in **0.58 s and 1.90 s**, where 1600 sections previously exceeded a ten-minute budget.

**End to end on a transient, best of 2, circuits built once and reused:**

| n | dense | SuperLU | KLU | SuperLU win | **KLU win** |
|---|---|---|---|---|---|
| 152 | 3.581 s | 4.378 s | 6.832 s | 0.82x | **0.52x** |
| 402 | 12.413 s | 9.524 s | 12.444 s | 1.30x | 1.00x |
| 1002 | 53.872 s | 43.462 s | **35.263 s** | 1.24x | **1.53x** |
| 2002 | 161.484 s | 117.111 s | **110.143 s** | 1.38x | **1.47x** |

`max|dv|` is 1.11e-15 on every row — all three solvers agree to round-off.

**KLU crosses over around n~400 and wins from n~1000, where it also beats SuperLU.** The
earlier reading — "KLU loses end to end, so `AutoSolver` selects SuperLU" — was correct for
the sizes then reachable and is now superseded. **This is the same shape of error the plan
keeps producing: a measurement taken at one convenient size, generalised into a rule.** It
was labelled as bounded at the time, which is the only reason it got revisited.

**`AutoSolver` is now size-aware**: dense below 100 unknowns, SuperLU from 100, KLU from
**1000** (`MIN_N_FOR_KLU`, set where the tables cross rather than at a round number). Tests
pin **both** directions so neither can be collapsed back into a single always-answer.

**D1 (reuse Newton's factors for the step controller) — STILL NOT TAKEN, and the case for
it has WEAKENED rather than strengthened.** The reasoning that deferred it was "33% of
solves removed buys ~2.7% of runtime at n~300, but ~20% at n=5000". The first half holds;
the second no longer follows, because **KLU's refactor has already taken most of what D1
was going to take.** D1 saves one factor+solve in three; with the ordering reused and a
refactor costing a quarter of a factorisation, the marginal saving is a quarter of what it
would have been against a dense LU. Against that it still costs bit-identity — the property
that made every safe refactor of this session possible. **Left unimplemented, now for a
measured reason rather than a provisional one.**




**7d. Delete `pybsmatrix.py`** (340 unreferenced lines, no pivoting, and a `fbsub` whose
division sits inside the wrong loop so it cannot be correct), and **fix
`test_sparse_toolkit.py` to construct the circuit with the toolkit under test** — it
currently passes while never exercising the sparse path. Fix the test first, then decide
whether `_sparse_numeric` is worth keeping given it is 4x slower than dense.
OUTCOME: **PASSED.** `pybsmatrix.py` is deleted, and `test_sparse_toolkit.py` now constructs
its circuits with the toolkit under test — it previously passed while exercising the dense
path exclusively, so it was reporting on a sparse toolkit that was wholly non-functional.
Fixing the test first, as the item required, is what exposed that.
### 7d OUTCOME, 2026-08-01 — the test was hiding a dead backend

**Done in the plan's own order: fix the test first, then delete, then decide.** Fixing the
test answered the decision outright.

**`pybsmatrix.py` DELETED.** 340 lines, **zero references** anywhere outside itself,
verified by grep across the tree. And the `fbsub` defect is exactly as recorded, now
located: `vec[i] /= self[i,i]` sits **inside** the `for k` loop, so it divides once per
off-diagonal term instead of once per row — and the back-substitution loop runs *forward*
rather than in reverse. `fbsub_transpose` repeats the same misplaced division. It could not
have produced a correct solve.

**`test_sparse_toolkit.py` FIXED — AND IT IMMEDIATELY FAILED.** The tests built
`SubCircuit()` with **no toolkit**, so every circuit was assembled by the *default*
(numeric) backend and only the analysis object was handed `sparse_numeric`. The matrices
reaching the solver were dense numpy either way: **the tests compared `numeric` against
`numeric` with a different solve call bolted on.**

**Built with the toolkit under test, `sparse_numeric` does not work at all:**

| analysis | numeric | sparse_numeric |
|---|---|---|
| DC | OK | `AxisError: axis 0 is out of bounds for array of dimension 0` |
| Transient | OK | `ValueError: setting an array element with a sequence` |
| AC | OK | `AxisError: ...` |

**Root cause:** `SparseNumericToolkit` makes `cir.G(x)` a `coo_matrix`, and
`analysis.remove_row_col` calls `numpy.delete` on it, which sees a 0-d object. Every
analysis removes the reference node, so nothing can complete. **Verified pre-existing** —
it fails identically at `def248c~1`, before 7a touched `remove_row_col`, so this is not a
regression from that work.

The three tests are now `xfail(strict=True)`, so the defect is **pinned and visible**, and
a future repair turns them red rather than passing unnoticed — the opposite of how it
survived the first time.

**So 7d's question is answered more strongly than it was posed.** It asks "whether
`_sparse_numeric` is worth keeping given it is 4x slower than dense". Two corrections:

1. **The 4x is not a constant.** Measured against `numeric` on a tridiagonal system:
   **9.36x slower at n=50**, 1.34x at n=150, then **7.7x FASTER at n=400** (0.13x) and 3.7x
   faster at n=800. It is the same fill-dependent crossover 7b measured, not a fixed
   penalty.
2. **It does not work**, which settles the matter ahead of any timing.

**RECOMMENDATION, NOT TAKEN HERE BECAUSE IT IS PUBLIC API: remove `sparse_numeric`.** It is
46 lines overriding only `linearsolver` — a sparse *solve* on dense-assembled matrices,
which is precisely what **7b's `SuperLUSolver` now does correctly**, and within noise of it
at n>=400 (0.00225 vs 0.00217 at n=400). The difference is that `AutoSolver` selects on
fill and stays dense where sparse would lose, while `sparse_numeric` applied
unconditionally. Removing it touches `toolkit.py` (2 sites), `circuit/__init__.py`,
`test_toolkit.py` and `test_sparse_toolkit.py`. **Nobody can be depending on it, because it
completes no analysis** — but deleting an exported name is the maintainer's call, so it is
left standing and marked.


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

### 8(a)/8(b) entry measurements, 2026-08-01 — taken before any code

**8(a): PARTLY REFUTED.** The plan says *"`per=0` (SPICE's one shot) and `tr=0`/`tf=0`
(SPICE substitutes TSTEP) both divide by zero"*. Measured on `Pulse.f` directly:

| case | result |
|---|---|
| `per=0`, `tr`/`tf` > 0 | **OK — `f(2e-6) = 1.0`**, no crash |
| `tr=0` | `ZeroDivisionError` |
| `tf=0` | `ZeroDivisionError` |

**`per=0` is already guarded** by `if self.per != 0` at the top of `f`. Only `tr`/`tf` divide
by zero, and the mechanism is as recorded: `Pulse.f` builds every `where()` branch eagerly,
so `(v2 - v1) / tr` is evaluated even when the branch is not selected. Also `Pulse` itself
has **no class defaults** — it takes seven required arguments; it is `VPulse`, the element,
that supplies `tr=0, tf=0, per=0`.

*The branches whose denominators vanish are unreachable anyway*, which is what makes a safe
denominator a correct fix rather than a papering-over: with `tr=0` the rise branch's
condition is `t < td + 0`, and the branch after it re-selects `v1` on exactly that interval,
so its value can never survive. The same holds for `tf=0`.

**8(b): CONFIRMED, and the mechanism is worse than "ignores `td`".** `Sin.f` computes
`exp(-theta*(t - td)) * sin(omega*(t - td) + phase)` with no gate on `t < td`. Two
consequences:

1. The sine runs from `t = 0` rather than from `td` — measured **0.9511 of full amplitude at
   `t = 0.2 td`**, where SPICE holds the source at its offset.
2. **For `t < td` the exponent is POSITIVE**, so the "damping" term grows backwards in time.
   Measured with `theta = 5000`, `amplitude = 1 V`: **`f(2e-4) = 51.93 V`**. Larger `theta`
   reaches the plan's recorded 2835 V. The growth is unbounded in `theta*td`.

SPICE's definition holds `V = VO + VA*sin(phase)` for `t < td`, which falls out of clamping
`t - td` at zero rather than needing a separate branch.

**Gate 8a-1 (the class defaults run).** A `VPulse` on its element defaults must complete a
transient. Declared: no exception, and the waveform is `v1` before `td` and `v2` after.
OUTCOME: **PASS.** A `VPulse` on its element defaults completes a 168-step transient reaching
`tend` exactly. Asserted that the pulse **fires**, not merely that the run finished: drive
reaches 1.000 and the RC charges to 0.9911 by t=5.999e-6, with the fall at 6e-6 as the
parameters say. "0 before, 0 after" is equally true of a source that never moved.

**Gate 8a-2 (non-zero `tr`/`tf` do not move).** Bit-identical on the 17-run baseline, which
contains three `VPulse` circuits with `tr = tf = 1e-8`. The fix substitutes a denominator
only where the branch is unreachable, so anything else means the substitution leaked.
OUTCOME: **PASS.** All 17 baseline runs `max|diff| = 0.000e+00`, and a test checks the substituted
denominator mid-ramp — where a wrong one shows up directly — rather than only at endpoints
where every candidate agrees.

**Gate 8b-1 (`VSin` holds before `td`).** `f(t) == offset + amplitude*sin(phase)` for every
`t < td`, which is SPICE's rule and is what clamping gives.
OUTCOME: **PASS.** `f(t) = 0.2500` at t = 0, 2e-4, 9.99e-4 for `offset=0.25, phase=0, td=1e-3` —
exactly `offset + amplitude*sin(phase)`, SPICE's rule.

**Gate 8b-2 (no growth, at any `theta`).** Declared as a bound rather than a value:
`|f(t)| <= |offset| + |amplitude|` for **all** `t` and all `theta >= 0`. A test that checked
one `theta` would pass against a fix that merely reduced the overshoot.
OUTCOME: **PASS.** Worst `|f|` over `theta` in {0, 1e2, 5e3, 1e5, 1e6} is **1.000000** against an
amplitude of 1.0 — the declared bound. Before the fix, `theta=5000` alone reached
**51.93 V**.

**Gate 8b-3 (`td = 0` is unchanged).** The overwhelmingly common case must be bit-identical;
`td = 0` makes the clamp a no-op, so anything else means the rewrite changed the formula.
OUTCOME: **PASS**, to `rel=1e-12` against the formula written out independently, with `theta=250`
and `phase=30` so the check is not degenerate.

**Gate 8ab-4.** Full suite.
OUTCOME: **PASS. Suite 854 passed, 6 skipped, 3 xfailed, 0 failed**, against 848 plus the six tests
added here.

**AN INTERMEDIATE RUN FAILED TWO EXISTING TESTS, AND THEY WERE RIGHT.**
`test_func.test_sin` and `test_elements.test_vsin` assert the SMOOTH symbolic form of the
source. Applying the clamp unconditionally made `Sin.f` return
`Piecewise((t - td, t > td), (0.0, True))` under the symbolic toolkit, which would have
broken every symbolic consumer — AC, transfer functions, the DDD machinery — for which a
time-domain start delay is meaningless in the first place. **Unlike `test_pulse` in 9(d),
these were not pinning the defect; they were protecting a path the fix had no business
touching.** The clamp is now gated on `toolkit.symbolic`.

**Gate 8-1:** `VPulse()` and `VSin()` with class defaults run to completion.
OUTCOME: **PASS** — see 8(a)/8(b) above.
**Gate 8-2:** `VSin(td=...)` holds at the offset for `t < td`.
OUTCOME: **PASS** — see 8(b) above.
**Gate 8-3:** a `TLine` DC solve gives the same answer before and after a transient, and
`max_step <= TD/2` is enforced per line.
OUTCOME: **PASS on both clauses, and one of the three entry claims was wrong.**

*Entry measurements.* DC gave **v(b) = 0.500000 before** a transient and **0.000000 after**
on a matched line where 0.5 is correct — `TLine.G` selects the DC stamp on
`len(self.history) == 0`, and the history was never cleared. History **accumulated across
runs**: 12 entries after run 1, **73** after run 2.

**REFUTED: "there is no `next_event`".** `TLine` *has* one, inherited from `Circuit`, and it
returns `inf` — so it contributes no breakpoints, which is the same substance stated
wrongly.

**NARROWER THAN RECORDED: the 4x/10x delay error needs `fixed_timestep`.** Under the
adaptive default the controller rescues it — measured **1.01-1.05x TD** at every timestep
from 5e-11 to 5e-9. Under `fixed_timestep` the claim reproduces almost exactly: **2.00x at
dt=1e-9, 4.00x at 2e-9, 8.00x at 5e-9** (recorded as 4x and 10x). *The defect only bites the
configuration nobody checks*, which is why it survived.

*What shipped.* Two hooks on `Circuit`, mirroring the existing `accept_step` idiom — base
no-op, `SubCircuit` propagates, `TLine` overrides:

* `reset_state()`, called at the start of every analysis, so which analysis ran before
  cannot change this one's answer. DC now returns 0.500000 both before and after a
  transient, and history is 12 after both runs.
* `max_timestep()`, returning `TD/2`, with `SubCircuit` taking the **minimum** over children
  — two lines with different `TD` both have to be resolved and the shorter governs. Measured:
  `max dt` pinned at exactly 5.0e-10 whatever timestep is requested, delay 1.01x.

Under `fixed_timestep` the cap **warns and obeys** rather than overriding: the caller has
taken the grid into their own hands, but a silently wrong delay is worse than a warning.

**AND THE POSITION OF `reset_state` COST A TEST TO LEARN.** Placed after the initial
`accept_step(0.0, ...)`, it wiped the history that call had just seeded, so `TLine.G` saw an
empty buffer and stamped the line as a **DC short** — `v(p1)` came out 1.0 where 0.5 is
correct. `test_tline.test_tline_transient_reflection` caught it; verified it passed before
the change, so it was my regression and not a pre-existing one.

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
OUTCOME: **ALL THREE PASS — see "Gates 9-1, 9-2, 9-3" below for the measurements.** (a)
euler h^1, trap h^2, gear h^2, each exactly constant — and it found the 3/4 Gear-2
optimism, stage 4i's fix never having crossed to this backend. (b) passes, but only after
the counter added to state it was itself found defective: see the correction. (c) passes
once the opening step is ramped (9(g)).
**Gate 9-2:** a CPU/JAX agreement test in the suite. The stage-5 cross-check was run by
hand and written into prose; it is not in the suite, so the next divergence is invisible.
OUTCOME: **PASS** — and there *was* a divergence the whole time, the Gear-2 error constant.
Compared against the analytic solution as well as against each other, since agreement alone
is satisfied by two backends wrong in the same way.
**Gate 9-3:** a `VPulse` transient under JAX hits the pulse edges.
OUTCOME: **PASS** — every analytic edge in `(0, tend]` has a time point within 1e-12. Was
blocked until 9(d), which is why it had never been run.

## 9(f) — remove `lte_formula` from both backends

**Started 2026-07-31.** This is 4f's option (C), scheduled here when the staged approach
was approved, and it executes **0.2b's already-measured decision rule** rather than making
a new decision: 0.2b asked what the `J^-1` mapping costs per step, said "under 10%, delete
the JAX `estimate_lte` charge path", and measured **1-3%**. It does not depend on 9(a)'s
`_lte_kernels.py`, so it is takeable now.

**Entry measurement, taken before any code (gate 9f-0).** Declared: the JAX `'classic'`
path's tolerance should prove unreachable, it should never reject, and `dt` should run away
to `dt_max`.

OUTCOME: **The tolerance is confirmed inert; the runaway clause is REFUTED, and the
difference matters.** `lte_error_ratio` computes
`tol = trtol*(lte_rel*max(|q|,1e-12) + lte_abs)` with `lte_abs = 1e-6` — a **voltage**
floor used as a **charge** floor, i.e. one microcoulomb. Against physical node charges it
dominates completely, so the ratio cannot approach 1:

| `\|q\|` | tol | lte | ratio |
|---|---|---|---|
| 1e-6 C | 7.007e-6 | 1e-9 | 1.4e-4 |
| 1e-12 C | 7e-6 | 1e-15 | 1.4e-10 |
| 1e-15 C | 7e-6 | 1e-18 | 1.4e-13 |

But `dt` does **not** run away under `JAXTransient.solve()`, because line 628 sets
`dt_max = timestep`. So an open-loop controller there does not diverge — **it degenerates
to a fixed-step run at the user's timestep**, which on a smooth circuit gives a perfectly
plausible answer. *That is why the defect survived.* The runaway is real only in
`solve_batched`, where `dt_max = tend/10`.

**And the harm is modest on benign circuits, which is worth saying plainly.** On the RC
test, `solve()` at coarse timesteps gives `'classic'` a *smaller* final-point error than
`'ywr'` at two of three step sizes, and identical whole-waveform maxima at two of three —
the error there is dominated by the first step, which `dt_max` stops either formula from
refining. Only `solve_batched` separates them:

| formula | timestep | steps | max dt | max err | RMS err |
|---|---|---|---|---|---|
| ywr | 1e-5 | 23 | 5.0e-4 | 1.457e-2 | 6.486e-3 |
| classic | 1e-5 | 16 | 5.0e-4 | 1.724e-2 | 8.856e-3 |
| ywr | 1e-4 | 19 | 5.0e-4 | 1.397e-2 | 7.052e-3 |
| classic | 1e-4 | 13 | 5.0e-4 | 1.676e-2 | 8.890e-3 |

~20% worse error, and `'classic'`'s result barely responds to the requested timestep —
the open-loop signature. **So the case for removal is that the safety net is disconnected,
not that the answers are visibly wrong.** A reviewer who expects a dramatic waveform will
not find one, and should not go looking.

**Gate 9f-1 (the charge path is gone, not merely unselected).** `lte_error_ratio` and
`estimate_lte` deleted; zero references anywhere outside the plan.
OUTCOME: **PASS, but only after a second sweep.** `estimate_lte` and `lte_error_ratio` are
deleted and both raise `ImportError`. The first grep looked clean and was not: it matched
only `Integrator(lte_formula=...)`, and **five live call sites passed the argument
positionally** — `benchmarks/transient_decisions.py`, `stage0_2a_integrator_choice.py`,
`stage0_2b_lte_solve_cost.py`'s label, `transient.py`'s own `Parameter` help text (which
was advertising the removed knob to users), and four constructions inside `lte_dae.rst`'s
`exec-rst` blocks. Only comments and docstrings mention the name now.

**Gate 9f-2 (the JAX default is untouched).** A JAX transient that previously ran
`lte_formula='ywr'` must produce **bit-identical** step counts and waveforms after the
removal. Declared at bit-identical: anything else means the default path was disturbed.
OUTCOME: **PASS.** Both JAX runs bit-identical, `max|diff| = 0.000e+00`.

**Gate 9f-3 (the CPU is untouched).** 4f-D1 established the parameter is inert end to end,
so removal must be a no-op numerically. Same declaration: bit-identical waveforms and step
counts on a circuit whose opening ramp makes the grid non-uniform.
OUTCOME: **PASS.** All nine CPU runs — three circuits x three integrators, including `rc-pulse`
whose opening ramp makes the grid non-uniform — bit-identical, `max|diff| = 0.000e+00`.
Baseline captured before the first edit and compared after the last.

**Gate 9f-4 (the removal is loud).** `lte_formula=` must raise `TypeError`, not be
silently swallowed. A kwarg accepted and ignored is the "thin advertised feature" 0.1c
warns about, and this whole item exists because one of those hid a broken estimator.
**Reconsider if** anyone is importing pycircuit as a dependency — then this wants a
release with a deprecation shim instead, and the shim must warn rather than accept
silently. Nothing in this tree imports it that way.
OUTCOME: **PASS.** All three integrators raise `TypeError`, and `lte_formula` is absent from
`JAXTransient.solve` and `.solve_batched` signatures.

**Gate 9f-5.** Full suite `-m ""`; doc build verified by content. The tests that exist to
pin "`lte_formula` selects nothing" become unstatable and go with it — that is expected
churn, not a regression, and the invariant they protected is now structural.
OUTCOME: **PASS, and the churn is fully accounted.** Suite **797 passed, 6 skipped, 0 failed**
against 810 before. Every one of the 13 is attributable: `_LTE_CASES` went from six rows
(one per method x formula) to three, across three parametrized tests (-9); the two tests
whose whole subject was "the two selections agree" are deleted (-3); and the JAX
end-to-end test is no longer parametrized over the formula (-1). No test changed its
assertions.

## 9(d) — the breakpoint scan: a live hang, not a latent one

**Entry measurements, 2026-07-31, taken before any code.** All three claims confirmed, and
**one is worse than the plan states.**

*The dict-iteration bug is real, but in ONE entry point.* `JAXTransient.solve` (line 564)
iterates `for elem in self.cir.elements` — a dict, so it yields **string keys**, and
`hasattr('V1', 'next_event')` is False. **0 breakpoints, always.** But `solve_batched`
(line 431) iterates `.items()` and is correct. The plan treats this as one defect; it is
one defect in one of two copies, which is itself the argument for 9(a).

*The infinite loop is real, and `solve_batched` therefore HANGS TODAY.* `Pulse.next_event`
ends with `if tmod == 0: return t` — it returns `t` **itself** at `t = 0`. The enumeration
loop feeds the result back in, so it never advances:

    next_event(0.0)  = 0.0     <- cannot advance
    next_event(1e-9) = 0.0001
    100001 iterations without reaching tend; 1 distinct value appended

Confirmed end to end: `solve_batched` on a `VPulse` circuit **times out**, having printed
its entry message and never returned. So the plan's "fixing either alone converts a silent
wrong answer into a hang" understates it — **the correct copy already hangs**, and only the
buggy iteration in `solve()` was hiding the defect behind a different one. Two bugs were
cancelling.

*And the CPU backend has known about this for some time without fixing it.*
`transient.py:762` reads *"Ensure next_t_break strictly advances time to avoid infinite
dt=0 loops"* and calls `next_event` a second time at a nudged `t`. So the working backend
**works around** a `next_event` that may not advance, rather than requiring one that does.
That is the mechanism, and it is what 9(d) should fix.

*The convention already exists — `Pulse` is the only violator.* Audited every
`TimeFunction`: base returns `inf`; `Sin` returns `td` only when `t < td`; `PWL` does
`if pt > t: return pt`; the exponential source does `if t < td1: return td1`; `AM` and
`SFFM` return `inf`. **Every one is strictly greater than `t` except `Pulse`.**

**Gate 9d-1 (the invariant, not the instance).** `next_event(t) > t` for **every**
`TimeFunction` subclass over a grid of `t` including 0, exact edge times, and exact period
multiples. Declared as a property test over the classes, not a case test on `Pulse` — the
point is the convention, and a new source type should fail this if it breaks it.
OUTCOME: **PASS**, over `Pulse` (periodic and aperiodic), `Sin`, `PWL`, `Exp`, `AM`, `SFFM`.
`Pulse` was the only violator, as the audit predicted.

**Gate 9d-2 (the hang is gone, and cannot come back).** `solve_batched` on a `VPulse`
circuit returns. Declared additionally: the enumeration loop carries its own progress
guard, so a future non-advancing `next_event` degrades to a wrong breakpoint list rather
than a hang. Two independent defences, because this failure mode costs a wall-clock
timeout to diagnose rather than a stack trace.
OUTCOME: **PASS.** `solve_batched` on the `VPulse` circuit returns in 3.6 s where it previously
ran past a 90 s timeout having printed its entry message. **The progress guard earned its
place immediately**: the first version of the fix still stalled, and the guard turned what
would have been another wall-clock timeout into a named `RuntimeWarning` at the exact `t`
where the contract broke. A test breaks the contract deliberately to prove the guard still
holds.

**Gate 9d-3 (`solve` actually finds the edges).** Breakpoint count > 0 for a `VPulse`
circuit under `solve`, and the collected times equal the analytic edge times.
OUTCOME: **PASS: 15 breakpoints, matching the analytic edge times exactly** (`k*per + {td,
td+tr, td+tr+pw, td+tr+pw+tf}` plus `tend`), where it previously collected none, ever.

**Gate 9d-4 (the CPU backend does not move).** The CPU transient's defensive re-call means
a strictly-advancing `next_event` should be invisible there. Declared at **bit-identical**
on the CPU baseline including a `VPulse` circuit. **If it is not**, the difference is
reported with the reason rather than absorbed — a changed breakpoint time is a changed
answer.
OUTCOME: **FAILED AS DECLARED — 6 of 17 runs differ — AND THE GATE WAS THE WRONG BAR.** Every
differing run is pulse-driven and adaptive; `rc-vsin`, `stiff-rlc`, both JAX runs and all
three *fixed-step* pulse runs are bit-identical, and no run changed its step count.
Differences are 4.2e-14 to 5.0e-10 on a 1 V-scale signal, i.e. 6 to 10 orders below the
`reltol` of 1e-4 that generated the grid.

**The reason, measured rather than asserted: the old breakpoint times were wrong and the
new ones are exact.** Walking the same sequence the CPU transient walks, including its own
`minbreak` re-call, against the analytic edges:

| | breakpoints landing EXACTLY on an analytic edge | max error |
|---|---|---|
| before | **7 of 15** | 4.235e-22 |
| after | **15 of 15** | **0** |

So the waveform moved because the step grid now lands on the true discontinuities. A gate
declared at bit-identical cannot pass a change whose entire purpose is to correct those
times — **the declaration was wrong, not the result**, and the honest record is that it was
declared before the drift was understood. What the gate did do correctly is force the
comparison that produced the table above; a looser gate would have waved the difference
through without anyone establishing which version was right.

**Gate 9d-5.** Full suite; doc build verified by content.
OUTCOME: **PASS. Suite 799 passed, 6 skipped, 0 failed**, against 797 before plus the two JAX
tests added here.

## 9(b) + 9(c) — tolerances exist, are settable, and are the CPU's

**Done 2026-07-31.** Entry measurements confirmed every claim: `JAXTransient(cir,
reltol=1e-6)` raised `KeyError: 'parameter reltol not in parameter dictionary'`, and
`nrsolver`/`scaler` were accepted in silence.

**What shipped.** A `parameters` list carrying `reltol`, `iabstol`, `vabstol`,
`lte_vabstol`, `lte_iabstol`, `TRTOL`, `maxiter` — **the CPU's names and defaults**, so
the two backends stop presenting as one analysis while solving to different accuracies
(the kernel's hard-coded `reltol=1e-3, abstol=1e-6` disagreed with `Transient`'s shipped
`1e-4 / 1e-12`). `nrsolver` and `scaler` now raise `NotImplementedError` naming why: the
Newton loop is a traced `while_loop` with a fixed algorithm and no scaling step.

**The absolute LTE tolerance is a per-row VECTOR**, `lte_vabstol` on node rows and
`lte_iabstol` on branch rows, built exactly as `Transient._solve` builds its own. 9(c)
warns in as many words that a scalar re-creates 0.3a's residual-vs-solution defect.

**AND THREADING IT ALONE INTRODUCED THE DEFECT THE PLAN WARNS ABOUT — caught by gate
9-1(c), not by review.** The floor `lte_vabstol = 1e-12` is only safe *because the CPU
ships `relref='sigglobal'`*; the JAX estimator measured `pointlocal` (each unknown against
itself, now). Threading the tight floor into a pointlocal reference is precisely the
combination `stepcontroller.py` records as pathological, and the measurement showed it
immediately — `min dt` collapsing 1000x as `reltol` tightened, for a **worse** answer:

| reltol | steps | min dt | late max err |
|---|---|---|---|
| 1e-3 | 56 | 4.76e-6 | 2.145e-3 |
| 1e-6 | 164 | **5.06e-9** | **2.802e-3** |

I had written in the same commit that `sigglobal` "needs a reduction over the whole vector
inside the traced loop, and that is 9(c)'s remaining work". **That was wrong** — a running
maximum is a scalar carried in `TransientState` and a `jnp.max`, both trivially traceable.
It is now implemented, and updated **only on accept**, since a rejected iterate is not a
signal the circuit ever had. With it, `min dt` holds at 5.5e-8 and step growth falls from
2.93x to 1.60x.

**Gate 9-1(c) (step count and error respond to `reltol`): PASSES ONLY ONCE `dt_max` IS
UNBOUND, and that is the finding.** `JAXTransient.solve` sets `dt_max = timestep`, so on a
run where the requested timestep is already fine enough the controller is subordinate to
the clamp and the whole-waveform error is pinned by the **first** step — `max err` sat at
4.2535e-3 to five figures across four decades of `reltol`, always at index 1. Unbind the
clamp and it responds properly:

| timestep | reltol | steps | max err | max dt |
|---|---|---|---|---|
| 1e-4 | 1e-3 | 53 | 4.792e-3 | 1.000e-4 (clamped) |
| 1e-4 | 1e-5 | 60 | 4.254e-3 | 1.000e-4 (clamped) |
| 1e-5 | 1e-3 | 502 | 1.169e-4 | 1.000e-5 (clamped) |
| 1e-5 | 1e-5 | 506 | **4.930e-5** | 1.000e-5 (clamped) |

**2.4x less error for a 100x tighter tolerance** — the controller works; it is just not
what is choosing the step. `max dt == timestep` on every row: the JAX transient is
effectively a fixed-step run with occasional trimming, and **cannot grow its step at all**.
That is a design property of `dt_max = timestep`, not a bug introduced here, and it is the
same clamp that made 9(f)'s broken estimator look plausible. **Whether `dt_max` should
default to `tend/10` as `solve_batched` uses is a real decision and is NOT taken here.**

**NEW FINDING, not previously recorded: the JAX transient overshoots `tend`.** Measured
`t[-1]` of 5.0559e-3, 5.0286e-3, 5.0702e-3, 5.0351e-3 against `tend = 5e-3` — up to 1.4%
past the requested end, where `Transient` lands on it exactly. The final step is not
clamped. Left open with its measurement rather than fixed inside this item.

**Still open in 9(c):** the Newton loop's convergence test is a scalar norm
(`conv_tol = abstol + reltol*F_norm0`), not the CPU's per-row residual test with `iabstol`
on node rows and `vabstol` on branch rows. `reltol`/`vabstol` are now threaded to it, so it
is settable, but the *flavour* split is done only for the LTE. Recorded so the remaining
asymmetry is not mistaken for finished work.

## 9(e) — the JAX Newton commits non-converged points, silently

**Entry measurement (gate 9e-0), 2026-08-01. CONFIRMED.** `newton_inner_loop`'s
`while_loop` exits on `F_norm <= conv_tol` **or** `iters >= maxiter`, and `time_body` takes
`nr_state.x` either way — so an unconverged iterate is committed as the step's solution and
its LTE is computed from it. On a **linear RC**, where the exact answer is known, squeezing
`maxiter` (settable since 9(b)):

| maxiter | steps | max err | diagnostic |
|---|---|---|---|
| 100 | 54 | 4.2835e-3 | none |
| 3 | 54 | 4.2835e-3 | none |
| 2 | 54 | 4.2835e-3 | none |
| **1** | 57 | **4.9708e-2** | **none — silent** |

**An 11.6x worse answer with nothing reported.** The run completes, lands exactly on
`tend`, and returns a plausible-looking waveform.

**Gate 9e-1 (non-convergence is detected at all).** The statistics must report it when
`maxiter` is squeezed. Declared: `maxiter=1` on this circuit must produce a non-zero count
where `maxiter=100` produces zero.
OUTCOME: **PASS.** `NewtonState` carries a `converged` flag out of the traced loop, and
`JAXTransientStatistics` reports `nonconverged_steps` (rejected and retried) and
`forced_nonconverged_steps` (accepted at the floor anyway).

**Gate 9e-2 (the response is to shrink the step, not to commit it).** SPICE's answer to a
failed Newton is a smaller step, not a committed non-solution. Declared: with `maxiter`
squeezed, the run must either recover by stepping smaller or say plainly that it could
not — **the silent 11.6x error must not survive**.
OUTCOME: **PASS**, after the first attempt traded one failure mode for another. A non-converged
step is now rejected and retried at `0.25x` — a fixed cut, because the LTE ratio that
normally sizes the shrink is computed from a point the circuit never occupied and means
nothing. `maxiter=1` now **raises `NoConvergenceError`** naming the time, the counts and
the three knobs, against silently returning a 4.97e-2 answer before.

**The first version turned the silent wrong answer into a hang** — exactly the trade the
plan flags under 9(d). Rejecting drives `dt` to `dt_min`, and a circuit that cannot
converge at any step size then advances by `dt_min` forever: the linear RC needed ~5e12
steps to reach `tend`, so a warning "after the run" is a warning that never prints. Fixed
by checking **per chunk**, where the loop is bounded by `CHUNK_SIZE`, and raising there.

**Gate 9e-3 (converged runs do not move).** A run whose Newton converges must be
**bit-identical** to before. Declared at bit-identical: this changes the accept predicate,
and anything else means it changed steps that were already fine.
OUTCOME: **PASS. All 17 baseline runs bit-identical**, `max|diff| = 0.000e+00` — 9 CPU, 6 CPU
pulse, 2 JAX.

**This gate caught a real regression in the first implementation.** Reading
`final.F_norm <= conv_tol` as "converged" is wrong: `F_norm` is the residual at the
iterate *before* the update that produced `final.x`, so the loop always does one update
beyond the one its test approved, and hitting the iteration cap reported failure even when
the returned point was good. Measured at `maxiter=2` on the linear RC, which had been
giving exactly the `maxiter=100` answer and started raising. The residual is now
re-evaluated at the returned `x` — one extra residual per step, cheaper than an iteration
since it needs no `G` and no linear solve — so the flag describes the point that is
actually returned.

**Gate 9e-4.** Full suite.
OUTCOME: **PASS. Suite 814 passed, 6 skipped, 0 failed**, against 812 before plus
the two tests added here.

## 9(a) — `_lte_kernels.py`, the shared scalar algebra

**Started 2026-08-01, last of stage 9 and deliberately so:** it is a pure refactor, and the
gates that make it safe (9-1..9-3, 9d, 9e, 9g) did not exist until this session. Doing a
deduplication with no JAX step-control coverage is how the two backends diverged in the
first place.

**The plan's count is off, and in the right direction.** It says *"five copies of the same
three VSS alpha coefficients"*. Measured: **three in source** —
`integrator.py:425` (Gear2 companion current), `integrator.py:510` (Gear2's LTE fallback,
reconstructing `g_n`), `jaxtransient.py:95` (`gear2_step`) — plus two in tests. One source
copy was already removed by 9(f) when `estimate_lte` was deleted. Recorded because a
duplication count is exactly the kind of number that gets quoted without recounting.

**Why this stage earned itself three times over.** Every defect stage 9 found was a fix that
existed on one backend and not the other: 4i's Gear2 error constant (gate 9-1(a)), stage 3's
opening step (9(g)), and 9(d)'s breakpoint scan. None was a new bug; each was a repair that
did not cross the copy.

**Gate 9a-1 (genuinely backend-agnostic).** The module must import neither `numpy` nor
`jax` — the formulas are elementwise arithmetic on `(q, C, g, h)` and need only operators.
Declared: an import check, plus the same kernel called with numpy and JAX inputs agreeing to
float tolerance.
OUTCOME: **PASS.** `_lte_kernels.py` imports nothing at all — asserted structurally, by reading its
own source, so the test fails for the right reason on a machine without JAX rather than by
failing to trace. Verified to trace under `jax.jit`, and the same kernel called with numpy
and JAX inputs agrees to float64 round-off (`bdf2_companion` exactly 0.000e+00;
euler/trapezoidal 3.6e-12 and 7.3e-12 on values of order 1e5, i.e. ~1e-16 relative).

**Gate 9a-2 (nothing moves).** **Bit-identical** on the full baseline — 15 CPU runs and the
JAX runs. Declared at bit-identical rather than "close": a refactor that changes a number is
not a refactor, and this is precisely the gate the rest of stage 9 was built to make
possible.
OUTCOME: **PASS. All 17 baseline runs bit-identical**, `max|diff| = 0.000e+00` — 9 CPU, 6 CPU
pulse, **and both JAX runs**. This is the gate the rest of stage 9 existed to make
possible: a refactor of the shared numerics, checked against a recorded baseline on both
backends at once, which was not expressible at the start of this stage.

**One expression did have to move, and it is recorded rather than absorbed.**
`jaxtransient.trap_step` computed `(2/dt) * (q_n - q_{n-1})` where `integrator.py` computes
`2 * (q_n - q_{n-1}) / dt` — algebraically equal, **two roundings against one**, so the two
backends already disagreed in the last bits. Unifying necessarily moves one; it moves the
JAX one, onto the form with fewer roundings. No production path reaches it today
(`eval_method` is hard-coded to `'gear'` at both call sites), which is why the baseline is
still bit-identical.

**Gate 9a-3 (the duplication is actually gone).** One definition of the VSS alphas in
source, reached by both backends. Declared as a grep, because "shared" is easy to claim and
easy to half-do — a kernel module that one backend quietly does not call is worse than two
honest copies.
OUTCOME: **PASS: one definition in source**, `_lte_kernels.py:42`, reached by both backends.
Asserted by *identity* rather than by value in the test —
`integrator.third_divided_difference is _lte_kernels.third_divided_difference`, and
`gear2_step` returns exactly `bdf2_companion` — because two copies that happen to agree
today is precisely the state this replaced, and a value check passes against it.

**One copy remains in the tests, deliberately.** `test_transient_repairs._gear2_true_lte`
derives the coefficients independently to compute the *exact* truncation error the
estimator is checked against. Pointing it at the kernel would make the estimator its own
reference, which is the one thing a test of numerics must not do.

**Gate 9a-4.** Full suite.
OUTCOME: **PASS. Suite 824 passed, 6 skipped, 0 failed**, against 817 before plus
the seven kernel tests added here. No existing test needed changing, which for a pure
refactor is the expected result rather than a reassuring one.

## 9(g) — the opening step, and reconciling `dt_max`

**Raised 2026-08-01, while explaining the `dt_max` question rather than acting on it.** The
question was mis-framed, and checking it produced a better answer.

**`dt_max = timestep` is NOT a JAX anomaly.** `transient.py:679` sets `max_step = timestep`
on the CPU as well, and measured, the CPU's `max dt / timestep` is **exactly 1.0000**. Both
backends clamp the step to the requested `timestep`, which is a defensible reading of what
`.tran tstep` means. The earlier framing — "the JAX transient is effectively fixed-step, is
that right?" — applies equally to the CPU and was not the real issue.

**The real issue is the OPENING step, and the CPU fixed it in stage 3.**
`Transient._opening_step` records it:

> A run used to open at `timestep`, which is also `max_step` — the largest step the
> controller is ever allowed to take. The step controller accepts the first step
> unevaluated... so that opening step was both the **largest** and the **only unchecked**
> step in the run. Measured: the global error was **1.3212e-01 at reltol 1e-3, 1e-4, 1e-5
> AND 1e-6** — identical to five digits — while the step count went from 24 to 195.

**Gate 9-1(c) measured the identical signature on JAX**: `max err` 4.2535e-3 at reltol
1e-3/1e-4/1e-5/1e-6, identical to five figures, **always at index 1**, while the step count
grew 53 → 85. `jaxtransient.py` opens at `current_dt = timestep` and has no `firststep`.

**So this is the third CPU fix that never crossed to the transcription** — after stage 4i's
Gear2 error constant (found by gate 9-1(a)) and stage 9(d)'s breakpoint scan. That is 9(a)'s
justification stated as measurement for the third time.

**Gate 9g-1 (the opening step is ramped, and settable).** `firststep` on `JAXTransient`
with the CPU's semantics — `None` means `timestep * 1e-3`, a non-positive value raises.
Declared: the same parameter name, default and validation as `Transient`, so the two do not
drift.
OUTCOME: **PASS.** `firststep` on `JAXTransient` with `Transient`'s name, default (`None` means
`timestep*1e-3`) and validation, and `_opening_step` ported with its reasoning. Applied to
both `solve` and `solve_batched`.

**Gate 9g-2 (`reltol` now moves the answer).** THE gate. Declared in advance: the error
must stop being pinned at 4.2535e-3 across four decades of `reltol`, and must fall
monotonically as `reltol` tightens. **If it does not, the diagnosis was wrong** and the
opening step was not what pinned it — say so and revert.
OUTCOME: **PASS, and the diagnosis was right.** The error was pinned at 4.2535e-3 across four
decades, always at index 1. Now:

| reltol | max err | at index |
|---|---|---|
| 1e-3 | 1.1141e-3 | 21 |
| 1e-4 | 1.1079e-3 | 21 |
| 1e-5 | 8.1991e-4 | 27 |
| 1e-6 | **1.8693e-4** | 52 |

Monotone, **5.96x smaller** across the span, and the worst error has moved off the opening
step into the real dynamics. Even the loosest tolerance improved 3.8x, from 4.2535e-3 to
1.1141e-3 — the ramp is not a trade, it is strictly better here.

**Gate 9g-3 (the CPU does not move).** Bit-identical on all 15 CPU baseline runs. This
touches `jaxtransient.py` only, so anything else is a surprise worth chasing.
OUTCOME: **PASS.** All 15 CPU runs bit-identical, `max|diff| = 0.000e+00`. The 2 JAX runs changed
shape (270 -> 305 and 2515 -> 2550 values), which is the intended effect: the ramp costs
~35 extra cheap steps and buys the tolerance response above.

**Gate 9g-4 (`dt_max` reconciled).** `solve_batched` defaults to `dt_max = tend/10` where
`solve` uses `timestep` — two entry points of one class disagreeing by ~50x. Declared:
`solve_batched` takes `timestep`, matching `solve` and the CPU, and the change in its
output is **measured and recorded** rather than assumed harmless, since it changes results
for existing callers.
OUTCOME: **PASS, and the old default was worse than "inconsistent" — it very nearly ignored the
requested timestep.** Same call, old `dt_max` against the new:

| timestep | dt_max | steps | max dt | max err |
|---|---|---|---|---|
| 1e-4 | tend/10 (old) | 29 | 5.000e-4 | 3.6745e-3 |
| 1e-4 | timestep (new) | 61 | 1.000e-4 | **1.1079e-3** |
| 1e-5 | tend/10 (old) | 32 | 5.000e-4 | 3.6680e-3 |
| 1e-5 | timestep (new) | 510 | 1.000e-5 | **1.2148e-5** |

**Asking for a 10x finer timestep used to change the run from 29 steps to 32, and the error
not at all** (3.6745e-3 vs 3.6680e-3) — `timestep` was very nearly inert on this path. It
now costs 2.1x/15.9x the steps and buys 3.3x/302x the accuracy. This does change results for
existing `solve_batched` callers, which is why it was measured rather than assumed.

**AND `solve_batched` WAS ALREADY BROKEN, WITH NOTHING TO SAY SO.** The counters added by
9(c)/9(e) default to scalars, but every field of the batched state must carry a leading
batch axis for `jax.vmap` — so `solve_batched` raised *"rank should be at least 1"* while
the whole suite stayed green. **It had no test coverage at all**; the only callers are
example scripts. It has coverage now.

**Gate 9g-5.** Full suite.
OUTCOME: **PASS. Suite 817 passed, 6 skipped, 0 failed**, against 814 plus the three tests added
here.

**Two existing tests changed premise rather than regressing, and both are worth recording.**
`test_jax_error_responds_to_reltol` compared at `timestep=1e-5` because before the ramp only
a fine timestep made the opening step small; with the ramp, what binds at a fine timestep is
`dt_max` itself, which leaves the tolerance no room (1.2148e-5 at reltol 1e-3 against
1.2141e-5 at 1e-5 — a ratio of 1.00). It now compares at a coarse timestep, over three
decades rather than two, because the two-decade response is 1.36x and too near the bar to
check anything. And `test_gate_9e_...` asserted `maxiter` cannot change a converged answer
*at all* (`rel=1e-12`); that held only while the opening step dominated everything, and is
now a 5% band.

## Gates 9-1, 9-2, 9-3 — the CPU's step-control gates, ported

**Done 2026-07-31.** The plan says of gate 9-1: *"None of these is currently expressible,
which is the asymmetry that let the copied LTE survive."* All three are expressible now
that tolerances are settable (9(b)/(c)), rejections are counted, and breakpoints work
(9(d)). **Porting them found the 3/4 Gear2 optimism for the THIRD time.**

**Gate 9-1(a) — the LTE scales with the right power of `h`.** PASS, and it is the gate
that paid for itself. Synthetic companion-current history, `J = I`, `x` fixed so `etol` is
constant, each method fed a `g` of **its own degree**:

| method | expected | estimate/h^order | observed exponent |
|---|---|---|---|
| euler | h¹ | 5.0000e5 constant | 1.000 |
| trap | h² | 1.6667e5 constant | 2.000 |
| gear | h² | 3.3333e5 constant | 2.000 |

*The first run of this gate got Euler wrong* — one quadratic `g` for all three makes
`g'(0) = 0`, killing Euler's leading term so it reads as second order. **The identical
degeneracy the CPU units test hit under 0.3d, walked into again a day later.** Trap 9 in
the plan is about units; this is its sibling and the harness now takes the degree per
method.

**AND IT FOUND THE REAL DEFECT.** `gear` read **2.5e5 = q'''/4** against the CPU's
**3.3333e5 = q'''/3**. That is YWR's Table I GEAR2 residual, which estimates
`(1/4) h² q'''` against a true `(1/3) h² q'''` — **the solver reported 3/4 of its own
truncation error at every step**, on the only `eval_method` either entry point uses.
The CPU found and fixed this in **stage 4i**; the fix never crossed to `jaxtransient.py`.
So the same defect has now been found three times in two transcriptions, which is the
argument for 9(a) stated as a measurement rather than a preference. Fixed with 4i's form —
`q'''` from a second divided difference of `g`, times the method's own error constant, so
the coefficient is derived rather than transcribed — and pinned against the derived
constant, not a recorded number.

**Gate 9-1(b) — a step is actually rejected. THE COUNTER NOW EXISTS; THE ANSWER IS STILL
NO.** A rejected step advances neither `t` nor `step_idx`, so it leaves no trace in the
output buffers and the gate could not be *stated* before, let alone passed. With
`JAXTransientStatistics.rejected_steps`, the answer across **16 configurations** — RC,
pulse with 1 ns edges, a stiff two-pole circuit, an oversized initial step up to
`timestep = tend`, and a sine at two frequencies, tolerances to 1e-8 — is **zero
rejections, everywhere**. The controller shrinks *predictively* on the accept path and
`dt_max = timestep` caps growth, so the error ratio does not exceed 1 in practice; the
opening step is force-accepted by `first = step_idx < 1` regardless of error. Two
contributing facts: 9(d) *removed* a likely source, since steps now land **on** pulse
edges instead of crossing them, and an extreme tolerance does not force a rejection either
— it collapses `dt` toward `dt_min` and makes the run effectively non-terminating.

**So the reject branch is untested code that no measured circuit reaches.** The test added
asserts the *bookkeeping*, deliberately **not** that rejections are zero: pinning zero
would freeze in a property of the current clamp, and if `dt_max` is ever changed — an open
decision — rejections should start and must not read as a regression.

**CORRECTED 2026-08-01, AND THE ORIGINAL FINDING WAS AN ARTEFACT OF MY OWN COUNTER.**
`do_accept` builds a fresh `TransientState` rather than `_replace`-ing one, so every field
it does not name reverts to the NamedTuple default — and `n_rejected`/`n_nonconverged` are
cumulative counters that it did not name. **The count was zeroed on every accepted step**,
so it could only ever report rejections that happened since the last accept, and the value
read at the end of a run was structurally 0. "Zero rejections across 16 configurations" was
a measurement of that bug, not of the solver.

With the counter carried through `do_accept`, **the reject branch is live**:

| circuit | reltol | accepted | rejected |
|---|---|---|---|
| rc | 1e-4 / 1e-6 / 1e-8 | 60 / 92 / 389 | **0 / 4 / 182** |
| pulse | 1e-4 / 1e-6 / 1e-8 | 81 / 162 / 480 | **0 / 68 / 241** |
| sine | 1e-4 / 1e-6 / 1e-8 | 73 / 384 / 1898 | **34 / 118 / 662** |

**So gate 9-1(b) PASSES**: a step is actually rejected, on every circuit once the tolerance
is tight enough to make the controller work for it. The reject path is not dead code and
the reasoning built on top of the artefact — that the predictive controller plus the
`dt_max` clamp made rejection unreachable — was wrong.

*What the episode is worth.* The counter was added **in order to** state this gate, and it
carried a defect that produced exactly the answer that needed no further explanation. A
plausible zero is the most dangerous reading a new instrument can give, and nothing else in
the suite could have contradicted it. The lesson is trap 13 below: an instrument built to
answer a question has to be shown to be able to answer it *differently*.

**Gate 9-2 — a CPU/JAX agreement test in the suite.** PASS. The stage-5 cross-check had
been run by hand and written into prose, so the next divergence was invisible — **and
there was one the whole time**, the Gear2 constant above. The test compares both backends
against the analytic solution *and* against each other on a shared grid; agreement alone
would be satisfied by two backends wrong in the same way, which is exactly what a copied
transcription produces.

**Gate 9-3 — a `VPulse` under JAX hits the pulse edges.** PASS: every analytic edge in
`(0, tend]` has a time point within 1e-12. Blocked until 9(d), which is why it had never
been run.

**Also fixed here: the `tend` overshoot, in two parts.** `calculate_next_dt` was passed
`state.t` — the time the accepted step *started* — so the breakpoint clamp sized the next
step from the previous position and overshot by about one step. And after that, a residual
one-timestep overshoot remained because `time_cond` used `t < tend` while the breakpoint
filter treated `t` within 1e-12 of `tend` as already reached: the two disagreed, so the
loop took one more full step. Both now use the same epsilon. Measured `t[-1]` = exactly
`tend` on every configuration, against 5.0559e-3 / 5.0286e-3 / 5.0702e-3 / 5.0351e-3 for a
requested 5e-3 before.

**And a chunk-boundary defect found while wiring the counters:** `sig_max` and
`n_rejected` are running totals, but `TransientState` was rebuilt per chunk without them —
so the `sigglobal` reference reset to zero every `CHUNK_SIZE` steps and a long run
silently reverted to a `pointlocal`-like tolerance at each boundary. **Same shape as the
CPU's `_dt_last2` reset**: a per-run quantity re-seeded by a per-call constructor. Pinned
with a test that runs the same transient at `CHUNK_SIZE=7` and `5000` and requires the
same reference and step count.



**One existing test had to change, and it was pinning the bug:** `test_func.py::test_pulse`
asserted `next_event(0) == 0` — precisely the non-advancing return that caused the hang.
Every other assertion in that same block already reads "strictly the next edge after `t`"
(`next_event(td) == td+tr`, `next_event(td+tr+pw+tf) == per`), so the line was the odd one
out rather than a considered contract. It now asserts `td`, and the block gained a walked
invariant check — walked rather than sampled, because the residual defect only appeared
after the accumulated value drifted off the edge grid, which no fixed sample set would have
caught.


**An intermediate run failed 20 tests, all the same cause**, and it is the reason gate
9f-1 got a second sweep: the positional call sites. `TypeError: takes 1 positional
argument but 2 were given` is at least loud, which is the argument for 9f-4's choice.

Doc build: `lte_dae.rst`'s "Selecting the formula" section is replaced by "One estimator,
not a choice", and its build-time table is now per-integrator rather than per-formula —
still generated at build time, so the numbers are re-measured rather than transcribed.
**The doc build exposed the failure mode rule 11 exists for**: `exec-rst` catches an
exception, renders the block's *source* instead, and the build still exits 0. Two blocks
had been silently degraded that way. The directive does emit a WARNING, which is the only
reason it was caught — and the warning count is why it was looked at at all.

---

# STAGE 10 — missing analyses

Scope per decision 0.3c. Ranked by the review's assessment of value:

1. DC sweep (`.dc`) — the most conspicuous absence; no `DCSweep` class exists.

**OUTCOME 10.1, 2026-08-01. `DCSweep` implemented.** SPICE's `.dc`, and the review is right
that it was the most conspicuous absence: there was no way to ask for a transfer curve, an
I-V characteristic or a bias sweep without writing the loop by hand.

**`DC.solve` gained an optional `x0`**, because it had none and a hand-written loop
therefore restarted every point from zeros. `None` is exactly the old behaviour, asserted
bit-identical by test.

**Continuation is the substance, and it took a detour to prove.** Seeding each point with
the previous solution is what makes a sweep across a nonlinearity cheap. The first
measurement said it bought **nothing** — a diode sweep took 232 Newton calls with and
without. That was the circuit being too easy. On a 4-diode chain, 0-5 V in 101 points:

| continuation | time | residual evaluations |
|---|---|---|
| on | 0.114 s | **250** |
| off | 0.378 s | **1052** |

**4.2x less work, same curve to 3.7e-07.** `continuation=False` is kept so the difference
can be measured rather than asserted.

**AND THE DETOUR FOUND A REAL DEFECT.** While checking whether `x0` reached the solver at
all, three DC solves on one circuit gave **15, then 2, then 2** residual evaluations. `None`
and an explicit `zeros` must be identical, so the difference was state: `Diode.limit` stores
`_vlim` on the instance and `G` linearises around it, so **the limiter's remembered junction
voltage outlived its analysis**. That is precisely the defect stage 8(d) found in
`TLine.history` — and `reset_state()`, made general there rather than a `TLine` method, is
the hook that fixes it. Now 15, 15, 15. **It was also masking the continuation measurement**,
which is why the first attempt read as no benefit.

**The sweep restores the swept parameter in a `finally`**, so a sweep leaves the circuit as
it found it; a test forces a mid-sweep failure to prove the restore actually runs. Without
it every later analysis would silently inherit the last swept value — the same shape of
defect again.

**A flaky test of my own was replaced.** 2+.5's `test_construction_is_no_longer_quadratic`
asserted a fitted timing exponent; it passed alone and failed inside the full suite at low
load. A flaky test is worse than none, because the next person learns to re-run it rather
than read it. The mechanism 2+.5 fixed is countable, so it now counts `update_node_map`
rebuilds (which must not grow with circuit size) instead of timing anything.

2. `.tran` output control — `timestep` is currently *both* the initial step and `max_step`,
   and output points are the solver's own non-uniform steps, so every FFT needs hand
   resampling. `nonlinear_leapfrog_sweep.py` does exactly this, and interpolating
   non-uniform data before an FFT is a real accuracy hazard the user is forced to own.

### 10.2 OUTCOME, 2026-08-01 — output control, and two of three problems already gone

**Two of 10.2's three stated problems had been fixed by earlier work**, so what remained was
smaller than the item reads. It says *"`timestep` is currently **both** the initial step and
`max_step`"* — `firststep` (stage 3, and 9(g) for JAX) and `max_step` (decision D2) both
exist and default to `None`. Only output control was left.

**The hazard is real and the plan names the victim correctly.**
`benchmarks/nonlinear_leapfrog_sweep.py` contains

    spec = np.fft.rfft(np.interp(grid, t, v)) / npt

— `np.interp` is **linear**, applied to a solution the integrator computed to second order.
Measured on the solver's own grid (a **1000x** spread in `dt` on an ordinary driven RC),
interpolating a *known* signal so only interpolation error is present:

| signal | linear | quadratic | ratio |
|---|---|---|---|
| fundamental | 1.12e-3 | 8.11e-5 | **13.8x** |
| 5th harmonic | 1.55e-2 | 1.20e-2 | 1.3x |
| decaying exponential | 2.00e-4 | 4.90e-6 | **40.8x** |

**The 5th-harmonic row is the honest one.** Where the adaptive grid barely resolves the
signal, no interpolant recovers what was not sampled — the fix there is a smaller
`max_step`, not a cleverer resample, and the docstring says so.

**What shipped.** `resample_uniform(t, x, npoints=|step=)` — three-point Lagrange, matching
the integrator's order, the same reasoning `TLine._interpolate_history` already gives for
its own lookup — and an `outputstep` Parameter on `Transient`. `None` keeps the solver's own
points, so **nothing existing moves**; measured, the default grid still has its 1000x spread
and `outputstep=1e-6` gives 5000 exactly-uniform points.

**Resampling happens after the run, deliberately.** The adaptive grid is what the error
control is defined on, so the solver goes on choosing its own steps and only the *reported*
points change. Landing the solver exactly on output points would be more accurate still and
is what SPICE does, but it changes the step sequence and therefore every existing result.

**A test asserts exactness for a quadratic**, not merely closeness to a sine: a linear
interpolant is arbitrarily close to a smooth signal on a fine enough grid, so "close to a
sine" would pass for the thing being replaced.

**Also fixed: `test_PSS_nonlinear_C` asserted nothing.** It called `pss.solve(...)` and
checked no property of the result — the same class of defect as `test_sparse_toolkit`
(passed while never exercising the sparse path) and `test_PAC` (skipped). It now asserts
**periodicity**, which is exactly what PSS claims to deliver and needs no external reference
so cannot be fitted: measured 2.2e-05 relative, asserted at 1e-3. A second assertion
requires the tanh knee to be crossed (v spans 18.8 V about a 1 V knee), because without it
the test would also pass on a degenerate solution that never leaves the capacitor's linear
region.

3. `.ic` / `.nodeset` — `uic=True` currently means "start from zeros", not SPICE's
   per-element initial conditions; oscillators and latches are unstartable.
   **PARTLY DONE 2026-08-02: node initial conditions ship as `Transient(..., uic=True,
   ic={node: volts})`.**

   The defect was not an inconvenience. `uic=True` meaning a vector of zeros leaves a class
   of circuit **unsimulable**: an LC tank at zero is *at* an equilibrium and stays there,
   and a latch at zero sits on its metastable point. There was no argument to `solve()` that
   started either. `test_initial_conditions.py` asserts the tank stays dead without `ic`,
   then oscillates at `1/(2π√LC)` with it — frequency measured from zero crossings against
   the analytic value, so it cannot pass on a circuit that merely wobbles.

   Three things raise rather than silently doing nothing: an unknown node name, naming the
   reference node (held at 0 V by construction), and `ic` without `uic=True`.

   **DEFERRED, and for two different reasons.** Element-level ICs — SPICE's `C ... IC=v` and
   `L ... IC=i` — are not implemented. `L`'s is a branch *current* whose unknown already
   exists in the MNA vector, so it needs only a reliable element-to-branch-index mapping;
   that is mechanical but not free, since `SubCircuit.branches` is flattened and does not
   record which element owns each entry. `C`'s is a branch *voltage*, constraining a
   **difference** of two node unknowns rather than either of them, and a set of such
   constraints is a spanning-tree problem rather than an assignment — a floating capacitor
   chain has no unique node-voltage solution without one.

   **`.nodeset` is not implemented either**, and is a genuinely different feature: a *hint*
   that seeds the DC solve and is then released, where `.ic` under UIC is a *starting value*
   that is never released. Stage 5's convergence-aid ladder is the place it belongs.

   **Inductor currents now ship too** (`L(..., ic=)`), via a recorded instance-to-branch
   span rather than a search — `Branch.__eq__` compares node pairs, so two parallel inductors
   produce equal branches and a search returns the first for both.

   **`C ... IC=` REMAINS, and it is the interesting half.** Full write-up in
   `doc/initial_conditions.md`: an inductor's IC is a branch *current* whose unknown already
   exists, so setting it is an assignment; a capacitor's is a branch *voltage*, constraining a
   **difference** of two node unknowns, so a set of them defines a system rather than assigning
   values. It is a spanning-tree problem — root each connected component at ground, propagate
   along a tree, check cycles for consistency — with a real decision to make about components
   that never touch ground, where the node voltages are determined only up to a constant.

   **Reconsider if** a circuit needs a floating capacitor's initial voltage, a starting
   current on an inductor inside a subcircuit, or a DC-convergence hint (`.nodeset`). None of
   the three is expressible by naming node voltages, and none is a workaround away.
4. A SPICE-subset netlist reader — everything else in interop is downstream of getting a
   circuit *in*.
5. Large-signal MOSFET — no CMOS transient is expressible today.
6. Waveform export (raw/PSF/CSV).

**Also, one line:** `Transient` is not exported from `pycircuit/circuit/__init__.py`.
**DONE 2026-08-02 — and it was two, not one.** `PSS` was missing as well; `DC`, `DCSweep`,
`AC` and `Noise` all arrive through the star-imports, so `from pycircuit.circuit import
Transient` — the import every transient example and every doc page uses — raised
`ImportError` while its neighbours worked.

Imported by name rather than with a star: `transient.py` and `shooting.py` both do
`from ...analysis import *` themselves, so a star here would re-export their transitive
imports and make the package namespace depend on what those modules happen to pull in.

`tests/test_package_exports.py` checks all six, in both the `hasattr` and the
`from ... import` form (they can disagree), and asserts the exported object *is* the class
rather than something shadowed. Nothing else in the suite imports the package the way a user
does, which is how a one-line defect survived stages 0-12.

---

# STAGE 11 — `shooting.py`

Blocked on 0.1c. Repair, rewrite against the seams stage 7 creates, or withdraw.

**OUTCOME, 2026-08-01. 0.1c's split recommendation applied: PAC WITHDRAWN, PSS REPAIRED.**
Both of stage 11's prerequisites had landed — stage 5 for limiting, stage 7 for the solver
seams — and each of 0.1c's four defects was re-checked before being acted on. **Two of the
four no longer reproduce**, and the one it ranked *first* turned out to be far worse than
recorded.

**Defect 1 — `method` is dead and PSS is backward-Euler only. CONFIRMED, AND THE COST IS
LARGE.** The Parameter appears only in its own declaration; `solve_timestep` hard-coded
`C/dt`. Measured on a series RLC driven **at resonance**, Q = 20, where the analytic peak is
`Q*va = 20 V`:

| steps/period | Euler | Euler % | trapezoidal | trap % |
|---|---|---|---|---|
| 20 | 2.63 V | **13.2%** | 19.32 V | **96.6%** |
| 50 | 5.61 V | 28.1% | 19.26 V | 96.3% |
| 100 | 8.81 V | 44.1% | 19.24 V | 96.2% |
| 200 | 12.20 V | 61.0% | 19.23 V | 96.2% |

**Euler understates the resonant amplitude by a factor, silently, and is still 39% low at
200 steps per period** — while trapezoidal is accurate and step-independent. Exactly 0.1c's
prediction. `method` now selects, validates its argument, and carries the companion current
through both sweeps; the first step of a sweep falls back to Euler because there is nothing
to average against.

**A CONTRADICTION APPEARED AND WAS CHASED DOWN RATHER THAN AVERAGED OVER.** On the low-Q RLC
of the existing test, Euler matched an AC reference at 0.9886 while trapezoidal read 1.3419
— the opposite verdict. That circuit has `R*C = 1 us` and the test steps at `dt = 1 us`:
**the step equals the time constant.** Refining it:

| dt/RC | Euler | trapezoidal |
|---|---|---|
| 1.000 | 0.9886 | 1.3419 |
| 0.250 | **1.3283** | 1.2121 |
| 0.063 | 1.0452 | 1.0461 |
| 0.016 | **1.0000** | **1.0000** |

Both converge. **Euler's agreement at `dt = RC` was coincidence** — at `dt = RC/4` it is
1.3283, *worse* than trapezoidal. So the resonator table is a statement about damping, not
about luck, and neither method is trustworthy at a step that coarse.

**The default stays `euler`**, so no existing result moves (measured: unchanged to 1.1e-16
relative — not bit-identical only because `(q - qlast)/dt` replaces `q/dt - qlast/dt`, one
fewer rounding). **Whether it should become `trap` is a maintainer's decision**, and the
evidence above is what it should be decided on.

**DEFAULT CHANGED TO `trap`, 2026-08-01, on the maintainer's decision.** The evidence above
is what it was decided on, and one further measurement was taken first — **on the existing
test's own circuit**, so the change is judged where the suite already looks rather than only
on the resonator that motivated it. `test_shooting` drives an RC at `dt/RC = 0.1`, and
against the analytic steady state:

| method | amplitude | relative error |
|---|---|---|
| euler | 1.669512 V | 1.4145e-02 |
| **trap** | **1.693006 V** | **2.7151e-04** |

**52x more accurate**, on a circuit chosen by someone else for a different purpose.

**No test needed changing.** `test_shooting`'s docstring records that `N` was chosen so the
discretisation error is a few per-mille with "a comfortable margin" — it was already
stepping finely enough for the method to matter, which is why it passes either way and why
its error simply falls.

**This is a deliberate behaviour change and existing PSS results will move.** They move
*towards* the analytic answer in every case measured: 13.2% -> 96.6% of a Q=20 resonator's
amplitude, and 1.41e-2 -> 2.72e-4 relative error on the RC above. `method='euler'` remains
available for anyone who needs the old numbers.


**Defect 2 — no limiting, "fails on every circuit it exists for". DOES NOT REPRODUCE.**
Tested a diode driven to 25 V and the common-emitter BJT that stage 5 was written for: all
converge, and the BJT's PSS collector voltage matches its DC operating point to four
decimals (0.7244-0.7275 against 0.7260) with a small swing about it. **Stage 5 overtook
this**: the `_expl()` clamp lives in the *model*, so PSS gets it without calling
`cir.limit()`. The claim was true when written on 2026-07-30 and is not now.

**Defect 3 — a dense `np.linalg.inv` per timestep. CONFIRMED AND FIXED.** `inv(Jf) @ C @
Jshoot` is the solution of `Jf X = C @ Jshoot`; asking for it directly is faster and better
conditioned, since forming an inverse squares the condition number you then multiply
through. **Bit-identical** on both PSS cases. Matrix-free shooting (Telichevesky, Kundert &
White, DAC 1995) needs only products with the monodromy matrix and would be better still;
that is a rewrite, and this is the same computation correctly expressed.

**PAC — WITHDRAWN.** `solve` raises `NotImplementedError` naming the reason: it forms the
whole `(N*M)x(N*M)` operator densely — **419.5 GiB** at N=137, M=1000 — and had never been
validated, its only test being `@unittest.skip("Skip failing test")`. The body is kept,
unreachable, as the starting point for a matrix-free rewrite, with 0.1c's `dtype` slip
(`B = tk.zeros(L.shape)` giving float64 against a complex `L`) corrected so that starting
point is not itself wrong. A test now asserts the withdrawal.

**Also found: `test_PSS_nonlinear_C` asserts nothing at all** — it calls `pss.solve(...)`
and checks no property of the result. It passes whatever PSS returns. Recorded rather than
fixed here, because giving it a real assertion needs a reference to assert against and that
is its own small piece of work.


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

**OUTCOME (2026-08-01): NOT MET, and the condition above was measuring the wrong thing.**
`benchmarks/transient_review/stage12_entry.py`, three stress circuits × two tolerances.

*Rejections* — 1.9 / 0.5 / 2.3 / 0.3 / 3.2 / 0.7 %, worst **3.2%** against the stage's own
10% threshold, and falling as tolerance tightens. By the condition as written, the stage
does not start.

But re-reading the paper rather than this stub's summary of it shows the condition was
aimed at the smaller of Fang's two mechanisms. §4.1 attributes the headline result — 39%
fewer time points, 17% wall clock — to the **lower** bound of eq (15): *"The introduction
of the lower bound γ<sub>min</sub> prevents step sizes from being unnecessarily small."*
Rejection elimination is the secondary effect. The stub half-knew this (see item (c) in
**Work** below, which says the lower bound "is not decoration") and still set the entry
condition on rejections. So the honest entry measurement is **how many accepted steps are
wastefully small**, and that was measured too.

*Wastefully small steps* — **there are essentially none, for a structural reason.**
`IntegralController` predicts `h_next = h · safety · (1/err)^(1/p)`, which is a deadbeat
law: it drives the next step's normalised error to the fixed point `safety^p` exactly. With
`safety = 0.9` and `p` **recovered from the realised step ratios as 2.00** (not assumed —
the benchmark derives it from the running code), that fixed point is 0.81. Measured median
accepted error at `reltol = 1e-6`: **0.8096 / 0.8071 / 0.8050**, IQR **0.0054 / 0.0034 /
0.0076**, with **96.4% / 94.9% / 91.2%** of accepted steps within 5% of the target.

At `reltol = 1e-4` the clustering loosens (61.6 / 55.6 / 24.6%), but the off-target
population sits at `err < 0.01`, not spread between — those are steps clamped by
breakpoints and `max_step`, which are **not LTE-limited and which Fang's method cannot
enlarge either**: a step that ends on a breakpoint ends there whether its size was tested
or solved for.

**So both of Fang's mechanisms are already exhausted on these circuits.** The band's lower
bound exists to move steps up onto the tolerance; pycircuit's steps are already sitting on
it to within half a percent. This is not an argument that the paper is wrong — it is that
TISpice's baseline controller accepted anything under a threshold (§4.1 records the
comparison method redoing a step only above a normalised LTE of **4.63**) without resizing
toward it, and pycircuit's baseline does resize. **The 39% was available there and is not
available here.**

**This is a negative result about the payoff, not about the method**, and it does not by
itself decide the stage — the user has asked for the full algorithm, and the plan below
delivers it. It sets the expectation the gates must be read against: the honest predicted
gain on the measured circuits is **at most the ~3% of steps currently rejected**, and gate
12-2's "at least 20% fewer accepted steps" is, on this evidence, **not reachable** — see
the revised gates.

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

---

## 12.0 What is already there — surveyed 2026-08-01, before planning

The stub says the code "implements **neither** §3.1 nor §3.4" (from 0.1d). That is no
longer accurate, and the difference matters for sizing the work. `_solve_coupled`'s own
comment now reads:

> This is the paper's robust "approximate Newton" (sec. 3.4); it replaces the exact (N+1)
> Schur update, **which is very sensitive to step changes and collapses the step size**.

So a previous attempt at the exact bordered solve exists and was abandoned on evidence.
**That is a recorded negative result about eq (12)/(13) and it must not be re-attempted
blind** — rule 7, a fix reused from a previous failure needs its failure mode checked
first. Whoever starts 12B reproduces the collapse before trying to fix it.

What `_solve_coupled` actually does today: converge the circuit at `h`, evaluate the LTE via
the shared `IntegralController`, and if rejected, shrink and **re-solve from scratch**, up
to `MAX_LTE_ITERS`. That is a rejection loop with the retries hidden inside the time point
rather than exposed as rejected steps. It is *not* §3.4: §3.4's whole point is eq (18),
correcting the solution already computed using the existing factors instead of re-solving.

Measured against the paper, the gap is:

| Fang | present in pycircuit | |
|---|---|---|
| eq (17) Gear step prediction | **yes** — `IntegralController`'s law | |
| eq (15) two-sided band `γ_min τ ≤ ε ≤ γ_max τ` | **no** — one-sided `err ≤ 1` | 12A |
| eq (16) step-change damper `\|Δh\| ≤ η h` | **no** | 12A |
| eq (18) approximate-Newton correction (§3.4) | **no** — full re-solve instead | 12C |
| eq (12) bordered `(N+1)` system | **tried, collapsed** | 12B |
| eq (13) reduced `N`-system + rank-one update | **no** | 12B |
| Fig. 4 two-stage Newton (bordered only when LTE unsatisfied) | **no** | 12B |
| `p = ∂f_ckt/∂h`, `q^T = ∂f_lte/∂v`, `d = ∂f_lte/∂h` | **no** — nothing forms these | 12B |

## The plan — four stages, each gated against the one before

Sequenced cheapest-and-most-likely-to-pay first, so that if the entry measurement's verdict
holds the work stops at a natural boundary instead of mid-way through a bordered solver.
**Every stage is behind the existing `coupled_lte` flag; the default path is untouched
until 12D.**

### 12A — the two-sided band and the step damper (eqs 15, 16). *Least invasive, do first.*

Fang's headline mechanism, and the only part that is nearly free. Add `γ_min`, `γ_max`, `η`
as `Parameter`s (the paper's `ltemin`/`ltemax`) and change the accept test from `err ≤ 1` to
the band, with `|Δh| ≤ η·h` limiting the predicted change.

Touches: `stepcontroller.py` only — the accept branch of `IntegralController.evaluate_step`
and three new parameters. `transient.py` passes them through. **No new derivatives, no
solver changes, no new state.** Roughly 40 lines.

**Gate 12A-1 (the band changes step placement at all).** With `γ_max` at the paper's ratio
to `γ_min`, the distribution of accepted `last_err` must widen measurably from the deadbeat
cluster the entry measurement found. Declared: if it does not move, the band is inert here
and 12B/12C cannot help either — **stop and report that**.

**OUTCOME (2026-08-01): PASS — it moves, but only the lower bound moves it.**
`benchmarks/transient_review/stage12a_band.py`, `reltol=1e-5`, three circuits.

| config | rc-vsin | stiff-rlc | rc-pulse | median err | note |
|---|---|---|---|---|---|
| baseline | 1288 / 19 rej | 490 / 4 | 2823 / 48 | 0.79–0.81 | aims at 0.81 |
| paper 0.7/3.0 | −1.1% / 34 | −1.8% / 25 | −3.4% / 251 | unchanged | 0.81 is *inside* the band |
| γ_max=3 only | −0.5% / **1** | 0.0% / **3** | −0.4% / **33** | unchanged | rejections down, steps flat |
| γ_min=0.95 | −9.6% / 471 | −10.6% / 114 | −11.4% / 1198 | **0.96–0.98** | aim forced up |
| η=0.15 | +0.2% / 17 | +5.3% / 4 | +5.1% / 48 | unchanged | costs steps, changes nothing |

The paper's own `0.7/3.0` is **inert on the aim point** exactly as the entry measurement
predicted — 0.81 lies inside it — and its only effect is *more* rejections, because steps
that drift under 0.7 now get redone. The distribution moves only when `γ_min` is pushed
above the deadbeat fixed point, which is what 12A-2 then has to price.

**Gate 12A-2 (it does not buy steps with accuracy).** Step count against stage 3's analytic
RC at matched final error. Declared: any step-count reduction must come with error no worse
than the standard path's. Given the entry measurement, **the expected honest result is "no
significant change", and that is a pass, not a failure** — it confirms the deadbeat
controller already occupies the band.

**OUTCOME (2026-08-01): the lower bound is a re-parameterisation of the tolerance — and the
more expensive one.** Measured against closed forms, not a fine-mesh run, so the integrator
is never compared with itself: `rc-vsin` has `v_b = A[sin(ωt−φ) + sin(φ)e^{−t/τ}]`, and
`stiff-rlc` satisfies `LC v'' + RC v' + v = 0` with `v(0)=1, v'(0)=0`.

| config on rc-vsin | steps | rejections | max abs error |
|---|---|---|---|
| baseline | 1288 | 19 | 2.014e-03 (1.00×) |
| **γ_min=0.95** | **1164 (−9.6%)** | **471** | **2.208e-03 (1.10×)** |
| **reltol × 1.173** | **1190 (−7.6%)** | **20** | **2.182e-03 (1.08×)** |

Those last two rows are the finding. `γ_min/0.81 = 1.173` was chosen *in advance* as the
tolerance scaling that moves the aim point by the same factor, and the two configurations
land on the same trade — within 2% on step count and within 2% on error — while `γ_min`
pays **24× the rejected steps** to get there. On a controller already sitting on its aim
point, raising `γ_min` does not recover waste, because there is none; it just moves the aim,
which is arithmetically a tolerance change, and doing it by rejection is the costlier way to
express the same thing.

**The one part that pays for itself is `γ_max > 1`**, which is not the mechanism the paper
credits: at zero accuracy cost (1.00×) and no measurable step change it cuts rejections
19→1, 4→3, 48→33. That is the two-threshold controller item 4a-bis already proposed, and it
is cheap — but it is worth at most the ~3% of steps the entry measurement found being
rejected.

**Two defects, both found by these gates, both invisible as exceptions:**

1. **Aiming at `γ_min` put the aim on the rejection edge.** The first implementation clipped
   the target *to* the excluded edge, so every undershoot rejected: **3172 rejections to
   accept 1187 steps**, to save 7.8%. Fixed by aiming at the band's geometric centre. Fang
   does not hit this because there `h` is *solved* to satisfy the band, not predicted and
   then tested — a predict-then-test scheme must aim strictly inside or it rejects on its
   own rounding. Pinned by `test_aim_point_moves_inside_the_band_not_onto_its_edge`.
2. **The damper was throttling the rejection path.** Eq (16) bounds the change between
   *accepted* steps; applied to error-driven shrinking it capped the retreat at 15% per
   retry, so the stiff ringdown exhausted `MAX_REJECT`, force-accepted, and crossed the
   whole ringing transient in **62 steps against the baseline's 490** — reporting an LTE of
   exactly zero, because by then it was integrating a signal that had already decayed.
   **The broken version's step count went DOWN, so a step-count-only gate would have read it
   as the stage's biggest win.** Pinned end-to-end by
   `test_damper_does_not_muzzle_the_rejection_path`, verified to fail against the defect.

**And the mechanism that let (2) hide is fixed, not just the instance:** accuracy was being
measured on `rc-vsin` alone, so a configuration that wrecked the stiff ringdown showed up
only as a step count. `stiff-rlc` now has its own closed form in the gate. Doing that
exposed a second thing worth recording: on that circuit the **total** error is dominated by
the opening step — which is accepted unevaluated by design, so no tolerance governs it — and
saturates at 1.3589e-02 across four decades of `reltol` (1e-5 … 1e-8, unchanged in every
digit). The benchmark now reports the error with and without it, because the all-points
maximum there is blind to everything step control does.

### 12B — the bordered system (eqs 12, 13, 14 and Fig. 4). *The invasive one.*

Make `h` a genuine unknown. Needs three objects nothing currently forms:

- **`q^T = ∂f_lte/∂v`** — tractable; `f_lte` is the divided-difference LTE in
  `_lte_kernels.py`, linear in the charges, so `q^T` is `∂lte/∂q · C`. **But the paper warns
  the controlling node changes:** *"the controlling LTE node can change from iteration to
  iteration"*, so `q^T` is re-formed each iteration, not cached.
- **`d = ∂f_lte/∂h`** — analytic from `gear2_lte`/`trapezoidal_lte`'s explicit `h`
  dependence, plus the implicit dependence through the divided differences. This is what
  the vestigial `analytical_eh` flag was named for.
- **`p = ∂f_ckt/∂h`** — **this is the expensive one, see the invasiveness section.**

Then Fig. 4's structure: solve the regular `N` system; estimate the LTE; **only if the band
is not satisfied** form and solve the bordered system. The `(N+1)` solve is not per
iteration, which is what makes the paper's overhead claim plausible.

**Gate 12B-0 (reproduce the recorded collapse first).** Before writing eq (13), restore the
exact `(N+1)` update and reproduce "collapses the step size" as a number — which circuit,
what `h` sequence. Declared: a described failure mode becomes a measured one, or the note
is wrong and that is the finding. **Do not skip this to save time; it is the cheapest step
in 12B and it decides whether the rest is repair or re-derivation.**
**OUTCOME (2026-08-01): the note was unsupported by the evidence available to it.**
`benchmarks/transient_review/stage12b0_collapse.py` reconstructed the 2026-07 loop from
`git show 9e8fb74^`. `dh` is **identically zero in every run using the estimator of the
time** — `E + TRTOL ~ 0` makes `E_h ~ 0` and `SchurCoupledNewton`'s `if abs(denom) < 1e-20:
dh = 0` guard fires. The code that carries the comment never moved the step at all and
cannot have observed the collapse it records.

**My own first reproduction is WITHDRAWN**: it paired the stage-4 LTE formula with the
pre-stage-4 charge-flavoured tolerance, a combination that never shipped, so the runaway it
measured says nothing about Fang's method. Restoring the band short-circuit verbatim
changed 54 steps to 56 — also not the cause. The real cause was found later and is
different from all of these: pycircuit was solving against the WRONG `f_lte` entirely
(gate 12B-1's outcome, and `doc/fang_stage12_conclusions.md` §3).

**Gate 12B-1 (the derivatives are right).** `p`, `q^T`, `d` each checked against a central
finite difference of the same residual — `p` by perturbing `h` and re-stamping. Declared:
agreement to the finite-difference floor on a nonlinear circuit, not just the RC. **This is
the gate that catches a sign error**, which in a bordered system shows up as a step size
that walks the wrong way rather than as an exception.
**OUTCOME (2026-08-01): PASS, all three, and finding the right `f_lte` first is what made
them cheap.** Eq (6) is a SOLUTION-space quantity, `eps_m = |v_i - v_i,extrapolated|`, not
the charge divided difference this module had been using. With that fixed:

- `q^T` is a signed unit vector on the controlling node — checked column by column against
  a central difference, including that every other column is zero;
- `d` is the extrapolation polynomial's derivative, shared with its value off one
  divided-difference table, checked against a central difference in `h`, sign included;
- `p` is `companion_dh + dudt`, checked by perturbing `h` and re-forming the WHOLE residual
  for all three integrators, plus two further tests that the source term is present and
  not at rounding level, since a `p` missing a term still passes a bare gradient check.

The degenerate case where the solution lands exactly on the extrapolation is pinned too: a
naive `sign(0)` zeroes the whole LTE row of eq (12) and makes the solve singular.

**Gate 12B-2 (eq 13 equals eq 12).** The reduced `N`-system with `Δh` recovered from eq (14)
must give the same `(Δv, Δh)` as a literal dense `(N+1)` factorisation, to solver tolerance.
Declared: measured on a circuit where `d` is small, since that is where the `1/d` division
in both eq (13) and eq (14) is worst conditioned. **`d → 0` is a real case** — it is an LTE
insensitive to the step, i.e. exactly the smooth region where steps grow.
**OUTCOME (2026-08-01): WITHDRAWN, then PARTLY REOPENED 2026-08-02 — eq (12) now works and
is selectable; eq (13)'s reduced form is still not implemented.**

The withdrawal below stands as written, and its reconsider-if turned out to name only one of
two routes. It said to reconsider "if an LTE is adopted whose `d` is not nearly equal and
opposite to `q^T dxh`". The second route is to keep this LTE and **not compute the
denominator by subtraction at all**: `d(eps)/dh = eps w'(h)/w(h)` in closed form, `w'/w`
being a sum of positive reciprocals. Measured against a ground truth from re-solving the
circuit at perturbed `h`: analytic +4.392e6 against a truth of +4.678e6 (ratio 0.939), where
the subtraction gives -9.680e5 — **the wrong sign**.

Eq (12) with that denominator is `coupled_method='bordered'`. Eq (13)'s rank-one reduced
system is still unimplemented and unneeded: the Schur form costs one extra solve against the
same `J`, which the factor/solve split from 7b already makes cheap.

---

*Original outcome, unchanged:* Eq (14)'s denominator `q^T dxh + d` is the solution's sensitivity to the step size
minus the extrapolation's slope; both are approximately `dv/dt`, so their difference is the
truncation error's derivative and is tiny by construction. Measured at `h = 1.6e-7`:
`q^T dxh = +1.818e9`, `d = -1.820e9`, denominator **-2e6**. Three digits lost and the SIGN
decided by the cancellation, so `dh` saturated at the eta limit with an arbitrary sign and
the step drifted down four decades while `err` sat at 0.2 — far BELOW the band that should
have grown it.

Eq (12) computes a small quantity as the difference of two large ones. Sec. 3.4's
approximate Newton takes the step from the error RATIO instead and has no cancellation, so
that is what ships. **Reconsider if** an LTE is adopted whose `d` is not nearly equal and
opposite to `q^T dxh` — the cancellation is a property of this estimator, not of eq (12) in
general.

**Gate 12B-3 (zero discarded time points).** With the bordered path selected, the rejection
counter must be 0 where the standard path records nonzero on the same circuit. This is the
stub's original gate 12-1 and it is the one claim of the method that the entry measurement
does **not** undercut.
**OUTCOME (2026-08-01): PASS.** `benchmarks/transient_review/stage12b_coupled.py`: the
rejection counter is **0** on both circuits at all three tolerances, where the standard path
records 8/19/22 and 4/4/5. Figure 3 has no rejection branch and the implementation now has
none either.

Worth separating: **the eq (6) estimator alone already removes most rejections**, before any
coupling — measured at 0/0/2 and 3/3/5 with `SolutionLTEController` driving the ordinary
adaptive loop. The coupled solve takes the remainder to zero.

### 12C — §3.4's approximate Newton as the default (eqs 17, 18)

Replace the current re-solve-from-scratch retry with eq (18)'s correction
`Δv^{k+1} = Δv^{k+1/2} − J^{-1} p (h^{k+1} − h^k)`, reusing the factors from 7b. The paper:
*"Without the need to modify the matrix solver of the simulator, the approximate Newton
method is straightforward to implement and carries very little overhead."* It needs `p`
from 12B but **no** `q^T`, `d`, or bordered solve.

**Gate 12C-1 (the correction is cheaper than the re-solve and no less accurate).** Declared:
per-time-point Newton iteration count drops, with the final solution matching the
re-solve's to Newton tolerance. If the correction needs a re-solve afterwards anyway to
converge, it has bought nothing — record that.
**OUTCOME (2026-08-01): PASS on accuracy, NOT MEASURED on cost.** Eq (18)'s correction
`dx = dxh * dh` reuses the stage-1 factors and replaced the old re-solve-from-scratch, and
the accuracy is better rather than worse (mean and median waveform error 10-20% below the
standard path at matched tolerance). **The per-iteration cost claim is gate 12-3 and is
still unmeasured** — everything here counts Newton solves, not wall clock.

### 12D — defaults, docs, and the four ignored inputs

Only reached if 12A–12C pass. Covers the stub's gates 12-4 (`fixed_timestep`, breakpoints,
injected controller, `uic` on the coupled path) and 12-5 (suite, and `time_stepping.rst`'s
warning removed only when the code matches the page).

## Invasiveness — the honest assessment

**Contained, except for one thing.**

Low risk, additive, behind `coupled_lte`:

- `stepcontroller.py` — band + damper (12A). Self-contained.
- `_lte_kernels.py` — new `bdf2_alphas_dh`, `gear2_lte_dh` alongside the existing kernels.
  The module imports nothing and that constraint holds for derivatives too.
- `transient.py::_solve_coupled` — the loop is rewritten, but it is one function that the
  default path does not call.
- `linearsolver.py` — eq (13)'s rank-one update wants a `solve`-with-existing-factors entry
  point. 7b already split `factor`/`solve`, so this is a use of the seam, not a change to it.

**The expensive part: `p = ∂f_ckt/∂h` requires a time derivative on every independent
source, and pycircuit's sources do not have one.**

`f_ckt` at `t + h` depends on `h` through two routes. The companion/integrator route is
analytic and cheap — `bdf2_alphas` is a closed form in `h1`, so `∂/∂h` of the Gear-2 stamp
is a few lines. The other route is the **sources**: `u(t + h)`, so `∂f_ckt/∂h ⊃ du/dt`.
`func.py`'s `TimeFunction` defines exactly `f(t)` and `next_event(t)` — **no derivative on
any of the six subclasses** (`Sin`, `Pulse`, `PWL`, `Exp`, `AM`, `SFFM`).

That is not merely six methods to write. **`Pulse` and `PWL` are piecewise linear: `du/dt`
is discontinuous precisely at the breakpoints, which are precisely the time points where
the step size is under the most pressure.** A bordered system whose `p` is evaluated at a
kink is solving for `h` against a derivative that does not exist there. Fang does not
discuss this because a breakpoint-clamped step is not LTE-controlled in the first place —
which is consistent with what the entry measurement found (the `err < 0.01` population at
loose tolerance is breakpoints), and it means the honest design is **`p` from the
integrator terms only, with the source route handled by leaving breakpoint-clamped steps
outside the coupled solve entirely.** That decision belongs in 12B and should be made
before any `du/dt` is written; a finite-difference `p` across a kink would be silently
wrong rather than loudly wrong.

**Blast radius on the default path: none until 12D**, and the two `stepcontroller.py`
changes already made for the entry measurement (`last_err` exposed and cleared on entry) are
the only edits outside the flag so far.

**Size estimate:** 12A ~40 lines and low risk. 12C ~60 lines, contingent on 12B's `p`.
12B is the bulk — new derivative kernels, the reduced solve, Fig. 4's restructuring, and a
recorded prior failure to reproduce first — call it the same order as stage 7b, and the
only stage here with a real chance of not working.

**Reconsider the whole stage if** 12A-1 shows the band inert (the entry measurement predicts
it will), or if 12B-0 reproduces the step-size collapse and the cause is the LTE estimate
rather than the solve. **A negative result remains a good outcome here**, for the reason the
stub already gives: it converts D1's declined "delete `_solve_coupled`" into an
evidence-backed decision.

---

**Gate 12-1 (it is the method, not a rejection loop).** *Superseded by gate 12B-3.*
With the coupled path selected, a
run must contain **zero** discarded time points from LTE rejection — the step size is
solved, not retried. Declared success: the rejection counter is 0 where the standard path
records a nonzero count on the same circuit.
**OUTCOME: see gate 12B-3 — PASS, zero discarded time points on both circuits.**

**Gate 12-2 (the paper's own result, on our circuits).** §4.1 reports **39% fewer time
points and 17% less runtime** on a Class-D amplifier with the LTE bounded between 0.7 and
3.0. Declared success: on the leapfrog and on `test_analysis_transient_stress.py`'s stiff
RLC, **at least 20% fewer accepted steps** than the standard controller at matched accuracy.
**A step reduction bought by a looser effective tolerance does not count** — accuracy is
measured against the analytic reference of stage 3's RC and must not degrade.

**REVISED 2026-08-01, before any code is written.** The entry measurement says 20% is not
reachable on these circuits and explains why: the paper's 39% comes from moving steps up
onto the tolerance, and pycircuit's steps already sit within half a percent of it. Leaving
the target at 20% would mean declaring a failure for a reason known in advance, which is
not what a gate is for. **The declared success is now: the step count does not get
*worse*, and the rejected steps — measured at 3.2% worst — go to zero (gate 12B-3).** If a
20%+ reduction does appear, the first suspicion is a loosened effective tolerance and the
accuracy check above decides it. **The original 20% target is kept in the record rather
than deleted, because the reason it was wrong is the useful part.**
**OUTCOME (2026-08-01): the step count goes UP by 25-28%, and the revised target above is
what it is measured against.** At `reltol=1e-6`, counting rejections: 5136 solves against
4089 on rc-vsin, 1924 against 1504 on stiff-rlc. The declared success after revision was
"the step count does not get worse, and the rejected steps go to zero" — **the second half
passes and the first half fails.**

The reason is §2 of `doc/fang_stage12_conclusions.md` and was predicted before the code was
written: the paper's 39% comes from moving steps up onto the tolerance, and this
controller's steps already sit on it to within half a percent. What was NOT predicted is
that the coupled path is more accurate per unit tolerance — mean and median error 10-20%
lower — so the extra steps are buying accuracy rather than being wasted.

**Gate 12-3 (the overhead claim).** The paper's case rests on "very little overhead". With
stage 7b's factors reused, declared success: the per-step cost of the coupled path is
within **15%** of the standard path's. If it is not, the method is losing in wall time what
it wins in step count, and that is the number that decides whether it ships.
**OUTCOME (2026-08-02): PASSES everywhere, once measured in the right unit.** Declared:
per-step cost within 15% of the standard path.

**Per Newton iteration the coupled path is CHEAPER on every circuit** -- 709 us against 781
on rc-vsin, 595 against 786 on stiff-rlc, 694 against 852 on rc-pulse: 9%, 24% and 18% less.

Its whole cost is doing MORE iterations, and the two reasons separate cleanly:

| | iter/point | total NR iterations | wall clock |
|---|---|---|---|
| rc-vsin  | 2.00 vs 2.01 | 10295 vs 8178 (+26%) | +14% |
| stiff-rlc| 2.17 vs 1.99 | 4241 vs 2983 (+42%) | +8% |
| rc-pulse | 3.74 vs 2.02 | 8310 vs 3030 (+174%) | +124% |

On the smooth circuits the iteration count per point is unchanged, so the extra work is
purely the extra time points. Only the breakpoint circuit needs materially more iterations
per point.

**THE FIRST VERSION OF THIS MEASUREMENT WAS IN THE WRONG UNIT AND IS WITHDRAWN.** It divided
wall clock by TIME POINTS, because `newton_iterations` was flat zero on the coupled path --
`fang_timestep` runs its own loop and never calls `_newton`, where the counter lives. That
gave "2135 us against 2473, the coupled path is 14% cheaper per point", which mixed a
genuine per-iteration saving with a difference in iterations per point and could not tell
them apart. The counter is now wired on both paths.

I also drew the wrong general conclusion from it: that the coupled path "runs several inner
iterations per point" was true of rc-pulse and false of the other two, where the counts are
2.00 against 2.01.

**Gate 12-4 (the four ignored inputs).** Whatever else changes, the coupled path must
honour `fixed_timestep`, breakpoints, an injected step controller, and `uic` — the list
from 0.1d, of which only `uic` is currently fixed.
**OUTCOME (2026-08-01): PARTIAL — breakpoints fixed, and they were worse than "ignored".**
`_solve_coupled` had NO breakpoint handling at all: no `next_event`, no truncation, no order
drop. That matters more on this path than on the standard one, because Figure 3 has no
rejection branch — a coupled step that runs past a pulse edge cannot back up from it. It
also broke the invariant `p` depends on, since `dfdt`'s right-hand limit at a corner is
exact only if a step STARTS on the corner.

It took two changes. Truncating alone is useless: `fang_timestep` solves for `h` and
promptly replaced the truncated step, measured at **0 of 10 pulse edges hit, worst miss
1.24e-7 s — the whole rise time**. A truncated step's size is imposed, so the LTE equation
is dropped for it (`hold_h`) and only the circuit is solved. Both paths now land 10/10 to
2.7e-20 s. Pinned by `test_coupled_breakpoints.py`.

Also fixed: the coupled path recorded only `accepted_steps`, so `breakpoints_hit` and
`order_drops` read as zero on a circuit hitting ten edges, and `tran.statistics` raised
`AttributeError` outright before that.

**`fixed_timestep` is now honoured too** (2026-08-02) — it was never even passed to
`_solve_coupled`. The grid is kept and the LTE equation is dropped on every step, the same
treatment a breakpoint-truncated step gets.

That needed one further distinction: a held step whose error is over the band is normally
reported so the caller can shrink and retry, but on a LOCKED grid shrinking is not an option
the caller has left available, so the step is taken and the accuracy is what was asked for.
Conflating the two broke `test_fixed_timestep_keeps_the_grid_on_the_coupled_path` — the
retry shrank `h` and the uniform grid disappeared.

**All four inputs are now honoured** (2026-08-02). The injected step controller was the
last: it was silently replaced by the path's own, so `tran.step_controller = X` looked
honoured and did nothing.

It is now **used** when it is a `SolutionLTEController` and **refused with the reason**
otherwise, because on this path the step-size law is Fang's and an injected controller's own
accept/predict logic is never consulted — accepting one and using it only for `relref` would
be the same silent no-op in a different costume.

The obvious test, `hasattr(injected, 'lte_gradients')`, is deliberately NOT used: `q^T` and
`d` are implemented and gated but are not called on the shipped path, since sec. 3.4 replaced
the eq (12) branch that used them. Keying on a method nothing calls would pass controllers
that cannot work and fail ones that can.

**Gate 12-5.** Full suite `-m ""`, and the citation in `time_stepping.rst` becomes true —
the page's `.. warning::` recording that the method was documented but absent is removed
only when the code matches it.
**OUTCOME (2026-08-02): PASS.** Suite **934 passed, 6 skipped, 3 xfailed, 0 failed**.

`time_stepping.rst` is rewritten. The page had gone from describing a feature that did not
exist, to describing its absence; it now describes what runs — eq (6) as a solution-space
quantity and why that distinction is load-bearing, Figure 4's two-stage structure, sec. 3.4
in place of eq (12) with the cancellation that forced it, the measured results, and the
limitations. The history is kept in a note rather than deleted, because a
documented-but-absent feature is the same class of defect as a silent wrong answer.

Two companion documents carry what the page cannot: `doc/fang_dac2013_math.md` (the paper's
math, extracted from rendered pages) and `doc/fang_stage12_conclusions.md` (the measurements,
the eleven defects and how each presented, the six negative results, and what is still
undone).

**Reconsider the whole stage if** gate 12-3 fails, or if stage 4 shows the LTE estimate is
still not trustworthy enough to solve against. **A negative result here is a good outcome**:
it would justify deleting `_solve_coupled` after all, which is what 0.1d recommended and D1
declined, and it would do so on evidence rather than on the code-reading that D1 found
unpersuasive.

---

# STAGE 13 — PCNR: the limiting replacement the docs already claim

**Added 2026-08-02, on the maintainer's request to review the paper and fix the
implementation if it is not a proper one.** Reviewed first; the review changed what "fix"
means, so it is recorded before any code.

**Source.** K. V. Aadithya, E. R. Keiter, T. Mei (Sandia National Laboratories),
*"Predictor/Corrector Newton-Raphson (PCNR): A Simple, Flexible, Scalable, Modular, and
Consistent Replacement for Limiting in Circuit Simulation"* —
`/home/andreas/pycircuit_agy/papers/aadithya.pdf`, 8 pages, read from 150/200-dpi renders.

## What the paper actually specifies (Fig. 2)

The solution vector gains a block: `x = [x_MNA ; x_lim]`, where **every limited quantity is
an unknown in its own right**. For the paper's two-parallel-diode example,
`x_lim = [v_D1, v_D2]` and the residual gains

    g_lim = [ v_D1 − (e1 − e2) ; v_D2 − (e1 − e2) ]

The Jacobian becomes four blocks, and the one that matters is `J_lim/lim = I` — the paper
marks it "identity sub-matrix", and it is what makes the elimination cheap. Devices evaluate
their currents at **their own `v_Dk`**, never at `e1 − e2`, which is what makes `g` and its
Jacobian consistent.

Each NR iteration is then two phases:

    (i)   Δx_MNA = −(J_MNA/MNA − J_MNA/lim·J_lim/MNA)⁻¹ (g_MNA − J_MNA/lim·g_lim)
    (ii)  Δx_lim = −(g_lim + J_lim/MNA·Δx_MNA)
    (iii) x_{i+1} = x_i + [Δx_MNA ; Δx_lim]                      ← "Predict"
          x_{i+1,lim} = refine(x_i, x_{i+1})                      ← "Correct"

`refine` is each device limiting **only the variables it owns**, which is why devices cannot
clash. (i) is the Schur complement of the bordered system, so the `Ax=b` solved is the same
size as the original MNA system.

## What is implemented — measured, not read off the docs

**`doc/src/circuit/pcnr_limiting.rst` states that pycircuit implements PCNR. It does not.**
The page makes three specific claims and none of them holds:

| claim on that page | actual |
|---|---|
| "elevates every limited quantity to a formal, explicit unknown" | no such unknown exists anywhere |
| "Predictor Phase … Corrector Phase" | one NR loop with a limiter callback |
| "uses a Schur Complement reduction … the block is strictly the Identity" | the only Schur code is `SchurCoupledNewton`, which is Fang's step-size coupling and unrelated |

What exists is classic SPICE limiting — the thing the paper is a replacement *for* —
implemented **twice, in two incompatible conventions**:

* `elements.Diode.limit` keeps private `_vlim` state and returns `None`, so
  `SubCircuit.limit` writes nothing back. This is the stateful form §2 criticises.
* `Semiconductor.limit` is state-free: it returns a limited vector, driven by a `junctions`
  table of `(anode, cathode, move)`.

**A third finding, WITHDRAWN on checking the record.** `ZenerDiode` declares
`junctions = ()`, so its inherited `limit` is a no-op, and that looked like a limiter which
exists and does nothing. It is a **deliberate decision**, recorded at gate 5(a) and honoured:
*"a forward-junction limiter would step straight through the breakdown knee, which is a
second exponential running the other way."* The device gets the direction-agnostic `_expl`
clamp instead. Withdrawn before it reached any code — but it is the third time in this
review that something reading as a defect turned out to be a decision already taken, which is
an argument for grepping the plan before writing a finding down.

**A hypothesis I tested and had to withdraw.** `SubCircuit.limit` applies limiters in
sequence, each writing back into `x`, so it looked as though two devices sharing a node
would clash with the last one winning — §2's complaint exactly. **Measured on two parallel
BJTs with saturation currents four decades apart, and it does not happen**: both insertion
orders give v(b) = 0.136066066, because the second limiter reads the *already-limited* value
and composes with it rather than overwriting. The structural concern is real; the
consequence I predicted from it is not, on this circuit.

## Gates

**Gate 13-1 (the claim is retired before anything else).** `pcnr_limiting.rst` must stop
asserting an implementation that does not exist, whatever else happens in this stage. A
documented-but-absent feature is the same class of defect as a silent wrong answer — this is
the second one found in this project, after the Fang citation in `_solve_coupled`.
OUTCOME:

**Gate 13-2 (one convention, not two).** `Diode` and `Semiconductor` must limit through the
same mechanism. Declared: the state-free form, because private per-device state cannot cross
a traced backend and because it is what PCNR's `refine` needs.

`ZenerDiode` is explicitly **out of scope**: its empty `junctions` is the recorded decision
above, not a gap.

OUTCOME (2026-08-02): **PARTIALLY MET, and the declared choice is REFUTED BY MEASUREMENT.**

*What was achieved.* There is now one implementation of `pnjlim`, in the new
`pycircuit/circuit/_limiting.py`, which imports nothing from the package so both `elements`
and `semiconductors` can use it without a cycle. Two copies of a numeric formula is how the
two conventions drifted apart in the first place.

*What was refuted.* The gate declared the state-free form — return a limited vector, let
`SubCircuit.limit` write it back, take `G` at the same point `i` is evaluated at. It was
implemented in full and **reverted**: `test_dc_pcnr_diode`, a 1 A source into a diode, stops
converging. Not diverging — **crawling**. After 100 iterations the update is 3.5x over
tolerance and the residual 1e4x over. The limiter is not at fault: traced separately, it
advances the junction 0.136 V per iteration exactly as `pnjlim` should. What changed is that
a consistent `G` taken at a limited point is a much smaller conductance than the residual
implies, so each Newton step buys back less than the limiter gives away.

*And the reason this matters more than the revert.* **Limiting cannot be made consistent by
moving a shared node voltage**, because the node is not the device's to move — everything
else attached to it sees the change. That is precisely the argument for PCNR's first key
idea, that each device gets its own unknown for the quantity it limits. So gate 13-2's
failure is evidence FOR stage 13's remaining gates rather than against them, and the
inconsistency it tried to remove (on a non-autodiff toolkit `Diode.G` linearises around
`_vlim` while `Diode.i` uses the node voltage) stays until 13-3 provides the variable that
can hold it.

**Gate 13-3 (the extra unknowns exist and are consistent).** With `x_lim` in the vector,
`g_lim = 0` must hold at convergence to solver tolerance, and every device must evaluate its
current at its own `v_Dk`. Declared: on the paper's own two-parallel-diode circuit, the
converged `v_D1` and `v_D2` agree with `e1 − e2` and with each other.

OUTCOME (2026-08-02): **PASS.** On the paper's own Fig. 1 — two parallel diodes with
saturation currents six decades apart, fed through 1 kΩ from a 1 V source — PCNR converges in
**8 iterations** with

    e2 = 0.346053885    v_D1 = 0.346053885    v_D2 = 0.346053885

agreeing to nine decimals with each other and with the branch voltage, and matching the
classic solver's operating point exactly. Consistency that changed the answer would be
worthless, so that second check matters as much as the first.

`pycircuit/circuit/pcnr.py` implements it, and **it needed no change to MNA assembly**.
Element stamps are additive, so the "everything except the PCNR devices" part of the residual
and Jacobian is the ordinary assembly minus each participating device's own stamp — measured
before writing the module, and it is why this is a layer rather than a rewrite. A device opts
in by declaring `pcnr_junctions` plus `pcnr_i(v)` and `pcnr_didv(v)`; `Diode` does.

Two properties are asserted beyond the residual, because a residual can be satisfied by a
system that is quietly double-counting: the diode contributes **nothing** to `J_MNA/MNA`
(its whole conductance moves to `J_MNA/lim`, since its current depends on its own unknown
and not on the node voltages), and the two currents divide in the ratio of their `IS`.

**A test of mine was wrong before the code was.** The independence test first drove `refine`
from `vold = 0`, where `pnjlim` takes its `VT·log(vnew/VT)` branch and does not depend on
`IS` at all — so two devices four decades apart limited to the same value and the assertion
was testing something the algorithm does not do. `IS` enters only through
`vc = VT·log(VT/(1.414·IS))`, and only when `vold > 0`. The test now steps to 0.6 V, which is
above `vc` for one device (0.432 V) and below it for the other (0.788 V), so one limits and
one passes through.

**Gate 13-4 (the Schur elimination costs nothing).** The `Ax=b` actually solved must stay
the size of the original MNA system. Declared: matrix dimension unchanged from the
non-PCNR path, and the per-iteration cost within 15% — measured per Newton iteration, not
per time point, which is the unit gate 12-3 had to be corrected to use.

OUTCOME (2026-08-02): **STRUCTURAL HALF PASSES, COST HALF FAILS — and the cause is this
implementation, not the method.** `benchmarks/transient_review/stage13_pcnr.py`.

*Structure: PASS.* The matrix actually factorised is `n × n`, unchanged from the classic
path, asserted on the object `predict` factorises rather than on a timing.

*Cost: FAIL.* On 60 diodes, **+60% to +80% per Newton iteration** across repeated runs. (Only
the ratio is meaningful: the absolute figures moved from 2044 to 4463 µs/iter for the classic
path alone between runs, so the machine load varies more than the effect being measured.)

*Attributed, not guessed.* Two rounds of optimisation were tried and neither closed it:

- **The dense Schur product was the obvious suspect and was not the cause.** Forming
  `J_ml @ J_lm` as an `(n,k)·(k,n)` product is `O(n²k)`, where each column of `J_ml` has two
  nonzeros and each row of `J_lm` has two — so it collapses to `k` rank-one updates touching
  four entries each. Implemented; the ratio did not improve.
- **Nor were the allocations.** Skipping the `(n,k)` and `(k,n)` blocks entirely on the sparse
  path, and caching the per-device parameter dict, also did not close it.

Profiling then gave the answer: classic assembly (`circuit.i` + `circuit.G`) is 1115 µs and
`augmented_system` is 2965 µs, of which the per-device subtract loop is only **14%**. **PCNR
here runs the full ordinary assembly and then does more work on top** — subtract each
participating device's node-voltage stamp, then re-stamp it at the device's own unknown.

That is inherent to building this as a *layer*. The paper's design does not pay it: a device
is asked for its PCNR contribution directly and never stamps the node-voltage form that would
have to be subtracted. Reaching the paper's cost means devices not stamping their nonlinear
part at all — an assembly change, and precisely the one this implementation avoided in order
to be a layer. **Reconsider if** PCNR is ever wanted as a default; the layer is then the wrong
shape and the assembly change is the work.

**Gate 13-5 (it converges where limiting does).** The existing limiting tests
(`test_dc_pcnr.py`, `test_analysis_transient.py`'s diode step) must pass unchanged, and the
six nonlinear stress circuits must still converge. **A method that is more consistent and
converges less often is not an improvement**, and that is the outcome to watch for.

OUTCOME (2026-08-02): **PASS.** Six circuits chosen to make limiting matter — the paper's two
parallel diodes, a 1 A current drive, a 20 V hard forward drive, a four-diode series chain,
six parallel diodes spanning six decades of `IS`, and sixty diodes. **PCNR converges on all
six, and all six answers agree with the classic path** to within 1e-6 relative.

Iteration counts are within one of each other throughout (8/8, 3/4, 8/8, 12/13, 12/13,
13/14), so the extra consistency costs essentially no extra iterations — the cost in 13-4 is
per-iteration work, not more iterations. The existing limiting tests pass unchanged.

**Reconsider the whole stage if** gate 13-4 shows the elimination is not free, or if 13-5
shows convergence degrading — in which case the honest result is to keep classic limiting,
fix gates 13-1 and 13-2, and record why PCNR was not adopted.

**Gate 13-6 (NEW — the Jacobian handed to the step controller is the one that was solved).**
Added 2026-08-02, prompted by the question *"is the ordinary device limiting active when PCNR
is used? I guess it should not be."* It should not be, and it is not — but asking exposed a
**remnant** of it that invalidated gate 13-5's headline transient measurement.

**The defect.** `Diode.G` linearises around `self._vlim`, which only `Diode.limit` writes.
PCNR never calls `limit`, because limiting is the thing PCNR replaces. So `_vlim` keeps
whatever value it was first given and the diode's conductance is frozen there. Instrumented
on a half-wave rectifier:

| path | `limit()` calls | `G()` calls | distinct `_vlim` | `_vlim` range | true junction range | max abs error |
|---|---|---|---|---|---|---|
| limiting | 4112 | 5978 | 3607 | −18.4731 .. 0.7435 | −18.4731 .. 0.7499 | 3.10e-01 |
| **PCNR** | **1** | **2283** | **1** | **0.0000 .. 0.0000** | −18.4712 .. 0.7498 | **1.85e+01** |

`cir.G(x)` therefore carried **no diode conductance at all** for the entire run.

**Why every test passed anyway.** Inside `augmented_system` the stale value is added by
`cir.G(x)` and subtracted again by the per-device loop, so it cancels exactly — the solved
system was always right, which is why gates 13-3 and 13-5 and every DC test were clean and
the converged answers agreed to 1e-6. The corruption escaped through one door only:
`_solve_timestep_pcnr` returned `J = cir.G(x) + Geq` to the **step controller**, which
computes `lte = J^-1 Eg`. That mapped the truncation error through a Jacobian missing the
nonlinearity. **The bug was unreachable from DC and invisible to every assertion that looked
at a voltage.**

**What it invalidates.** Gate 13-5's transient table reported PCNR taking 4.1-6.8x fewer
steps. That number was the defect. Corrected, with the Schur matrix — which `predict` already
factorises, and which at convergence (`v_lim == e_a − e_b`) *is* the Jacobian of the residual
— returned instead:

| circuit | limiting steps | PCNR steps | limiting max err | PCNR max err |
|---|---|---|---|---|
| half wave | 3314 | **3314** | 1.2456e-03 | **1.2456e-03** |
| full wave | 3418 | **3418** | 1.2735e-02 | **1.2735e-02** |
| charge pump | 1095 | **1095** | 4.5138e-03 | **4.5138e-03** |
| clipper | 693 | **693** | 1.3604e-01 | **1.3604e-01** |

Identical in every column, and the accepted-error distributions (median/p10/p90) match to
four figures on all four. Verified the PCNR path really runs rather than silently falling
back: 1766 PCNR timesteps vs 1 classic solve with `pcnr=True`, 0 vs 1767 with `pcnr=False`,
both giving `v(2) max = 9.335271001`.

**And identity is the correct answer.** The limiter and PCNR's `refine` change only the
*iteration path*; the converged solution of a Newton method is whatever the residual says it
is, independent of the route taken. Equal solutions give equal LTE, hence equal steps. **A
measured difference in step count between the two paths is now a defect signature, not a
feature** — which is what the regression test asserts.

OUTCOME: **PASS, and it retracts the stage's headline claim.** Three earlier conclusions fall
with it, all recorded rather than edited away:
  - the 4.1-6.8x step reduction — **withdrawn**, it was this bug;
  - "the DC gates structurally could not see this effect" — **wrong**, and itself derived
    from the buggy measurement; gate 13-4's off-by-default conclusion stands on its original
    evidence;
  - "limiting's LTE estimate is inflated by 25-30% at matched `h`" (0.8091 vs 0.6397 in the
    `1.58e-6 – 6.31e-6` bin) — **withdrawn**, the same artefact seen from another angle.

**Gate 13-7 (NEW — `pcnr=True` is honoured on the coupled path).** Found while tracing the
callers of `cir.G` for gate 13-6's residual hazard. `_solve` dispatches
`if coupled_lte: return self._solve_coupled(...)` **before** it ever looks at `pcnr`, so the
combination silently ran the classic limiter:

| pcnr | coupled | PCNR timesteps | `Diode.limit` calls | v(2) max |
|---|---|---|---|---|
| True | False | 1766 | 1 | 9.335271 |
| True | **True** | **0** | **4869** | 9.334558 |
| False | True | 0 | 4869 | 9.334558 |

Bit-identical to not asking for it. **A parameter that silently does nothing** — the same
class as the orphaned `h_next` and the blind `vabstol` guard.

**The wiring, and why it is small.** The coupled path writes its Newton out by hand, so PCNR
cannot attach the way it does to `solve_timestep`. It does not need to: the Schur-reduced
system is an `n`-sized system whose Newton step **is** `predict`'s `dx_mna`, and
`fang_timestep` already deletes the refnode and solves exactly that shape. Handing it
`(f_eff, J_eff)` makes its existing solve work unchanged — **and its bordered (N+1)
extension too**, since `dxh = -J^-1 p` reuses the same factors. The Schur reduction was
factored into `pcnr.schur_reduce` so the formula has one home; gate 13-6 happened precisely
because a second consumer open-coded it.

- **13-7a (PCNR actually runs).** PASS: 4880 augmented assemblies, `Diode.limit` down to 1
  (the DC operating point, which still uses the classic path).
- **13-7b (agreement with the classic coupled path).** **FAILED AS DECLARED — and the
  declared threshold was the wrong instrument.** Measured 5.58e-6 relative against a declared
  <1e-6. The threshold was imported from the non-coupled case, where the two paths are
  step-for-step identical (gate 13-6). On the coupled path they are not, and *should not be*:
  the step size is solved **inside** the Newton loop from the mid-iteration iterate
  `x_stage1`, which is the limited point under limiting and the full update under PCNR, so
  the grids legitimately diverge. `max|dv|` is then a grid-sampling artefact, and it is not
  even monotonic in `reltol` (1.10e-3, 5.21e-5, 2.64e-4, 6.65e-7 at 1e-4…1e-7) — which is the
  signature of an instrument measuring the wrong thing.

  **Restated and passed on the right instrument:** neither path is systematically less
  accurate, measured against an independent tight reference rather than against each other.
  Both give the *identical* maximum error, to six figures and at the same time point —
  `fig1` 7.98231e-05, `hard-drive` 4.37299e-03. This is rule 8 biting: comparing the
  integrator with itself could not answer the question, and an independent reference did.
- **13-7c (the stress circuits still converge).** PASS: all six converge both ways on the
  coupled path. Two showed a large direct `max|dv|` (`fig1` 4.67e-3, `hard-drive` 9.17e-3
  relative) — **and that was the same flawed instrument**, interpolation across a steep
  waveform, refuted by the reference measurement above.
- **13-7d (suite).** PASS: 973 passed, 6 skipped, 3 xfailed, 0 failed.

Pinned by `test_gate_13_7_pcnr_is_honoured_on_the_coupled_path`, which counts calls rather
than comparing waveforms — the two paths agree to within their own truncation error, so a
waveform comparison **cannot** tell "PCNR ran" from "PCNR was ignored". Verified to fail
against the defect reinstated by injection (`assert 0 > 10`).

**The residual hazard, recorded and not fixed.** `cir.G(x)` returns a *wrong Jacobian* during
a PCNR run, for any caller. One consumer existed and is fixed; nothing prevents the next one.
The clean fix is for devices to stop stamping their nonlinear part into `G` at all — which is
exactly the assembly change gate 13-4 declined in order to keep PCNR a layer. **Reconsider
if** a second consumer of `cir.G` appears on the PCNR path, or if PCNR is ever made default;
either turns this from a latent trap into an active one.

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
OUTCOME: **PASS at the time (744/6/0).** The count has since risen to **840/6/0** as later
stages added tests; the interface clause still holds — 7b added a `linsolver` argument to
all seven `solve_system` signatures at once, so no subclass diverged.

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
OUTCOME: **PASS at the time (744/6/0), and `relref = 'sigglobal'` is the shipped default.**
Now 840/6/0. Worth noting downstream: 9(c) had to port `sigglobal` to the JAX backend as
well, because threading the tight `lte_vabstol` floor into a `pointlocal` reference there
reproduced exactly the pathology this item fixed on the CPU.

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

## 2+.4 — assemble the reduced system directly, never stamping the reference node  **[WITHDRAWN]**

**Raised 2026-08-01 by stage 7a, and scheduled here rather than left in 7a's prose.**
7a made the reference-node *removal* cheap (2.0-4.6x); it did not remove the need for it.
Any approach that copies is O(n^2). Eliminating the copy means the reference row and
column are never assembled, which is an assembly-level change touching every element's
stamp — stage 2's territory, not the linear solver's.

**Why it is worth doing, in numbers 7a measured.** Per Newton iteration, on a numeric RC
ladder: assembly scales **n^1.26**, dense LU **n^1.83**, and the reduction **n^2.53** —
the worst of the three even after 7a. At n=20000 the *reduced* form still costs ~2.5 s per
iteration and copies 3.2 GB per call.

**This resolves an apparent contradiction between stage 2 and stage 7, worth stating
because both are right.** Stage 2 opens with *"at n=137 the solve is 2.1% of runtime and
assembly is 96%, so replacing the solver first would optimise the wrong thing"* — correct,
and 7a's profile reproduces it (LU 8% at n=302). Stage 7 opens with the memory wall at
n=5000. **The orderings differ because the exponents differ**: assembly dominates below
n~2000 and the solver above it, with LU's share running 35% at n=802, 60% at n=5000 and
77% at n=20000. Neither stage is wrong; each is stating the binding constraint at the size
it cares about. **Reconsider if** anyone quotes either percentage without its `n`.

**Gate 2+.4-1 (the reference row and column are never built).** Declared as an allocation
check, not a timing one: the assembled `G`/`C`/`J` must have shape `(n-1, n-1)` at the
point of assembly, so `remove_row_col` is not called on the Newton path at all. A timing
gate would pass on a version that still builds and discards them.
OUTCOME: **NOT RUN — 2+.4 was withdrawn before this gate was reached; see the withdrawal note below for why, and for the reconsider-if that turned out to be the wrong one.**

**Gate 2+.4-2 (bit-identical).** Stage 2's stop condition, inherited: waveform drift
exactly `0.00e+00` and an identical step count on the full baseline, or stop. Note the
ordering of the sum changes when a row is never stamped rather than stamped-and-dropped,
so this may prove *unachievable* — **if so, say that and stop**, rather than relaxing the
bar to fit. A reduced-assembly path that changes the last bits is a different decision and
needs its own justification.
OUTCOME: **NOT RUN — 2+.4 was withdrawn before this gate was reached; see the withdrawal note below for why, and for the reconsider-if that turned out to be the wrong one.**

**Gate 2+.4-3 (the symbolic and sparse toolkits still work).** `remove_row_col` has 20+
call sites and several toolkits reach it. Declared: the reduced-assembly path is opt-in
per analysis, and anything that has not been converted keeps today's behaviour.
OUTCOME: **NOT RUN — 2+.4 was withdrawn before this gate was reached; see the withdrawal note below for why, and for the reconsider-if that turned out to be the wrong one.**

**Reconsider the whole item if** 7b/7c land a sparse solver first: a sparse assembly that
never materialises the dense matrix makes this moot by construction, and doing both is
wasted work. **That is the more likely outcome and this item should be checked against
7b's design before anyone starts it.**

### 2+.4 WITHDRAWN, 2026-08-01 — on the maintainer's decision, and for a reason the item's own reconsider-if got wrong

**The recorded reconsider-if was NOT met, and saying so matters more than the outcome.** It
reads: *"if 7b/7c land a sparse solver first: a sparse ASSEMBLY that never materialises the
dense matrix makes this moot by construction."* What 7b actually shipped is a sparse
**solver** operating on a densely-assembled matrix — `AutoSolver` picks dense / SuperLU /
KLU *after* `G` and `C` have been built and reduced. **The reference row and column are
still stamped and still copied away.** Withdrawing on that clause would have been citing a
condition that never came true.

**The real case is arithmetic, and it is stronger.** After 7a replaced the two `np.delete`
passes with one blocked copy (2.0-4.6x), `remove_row_col`'s share of a transient — measured
at sizes 2+.5 has only now made reachable:

| n | transient | `remove_row_col` | share |
|---|---|---|---|
| 202 | 7.58 s | 0.0198 s | **0.26%** |
| 802 | 32.79 s | 0.2427 s | **0.74%** |
| 1602 | 73.16 s | 0.7280 s | **1.00%** |

The share does grow — it is still the worst-scaling component in isolation — but the growth
is *decelerating* across these points (2.8x per 4x in `n`, then 1.35x per 2x), so even at
n=5000 it is a low single-digit percentage. **Perfect elimination buys about 1%.**

**Against that, the cost is high and was known in advance.** Gate 2+.4-2 says in its own
text that never stamping a row changes the order of the summation, so bit-identity "may
prove **unachievable**". The item would therefore move every waveform in the package for
~1% — and bit-identity is the property that made 7a, 7b, 9(a) and 2+.5 safe to do at all.

**Withdrawn.** The gates stay recorded rather than deleted. **Reconsider if** an assembly
path is ever written that produces the reduced system *directly* — sparse or dense — because
then the copy disappears as a side effect of work being done for another reason, and the
1% is free rather than bought.


---

## 2+.5 — circuit CONSTRUCTION is O(N^2.27): 25 s to build an 800-element ladder

**Raised 2026-08-01 by stage 7b, where an `n=2000` measurement had to be abandoned.**
Written up as an item rather than left in a commit message, per the lesson of 2+.4.

**This is construction, not assembly, and the two must not be conflated.** Assembly —
`cir.G(x)`, `cir.C(x)`, run every Newton iteration — measures **N^1.26** and is fine.
Construction is `cir['R1'] = R(...)`, run once, and it is the one that is quadratic.
Measured on a numeric RC ladder:

| N | elements | construct | per element | vs linear |
|---|---|---|---|---|
| 25 | 50 | 0.049 s | 982 us | 1.00x |
| 50 | 100 | 0.163 s | 1631 us | 1.66x |
| 100 | 200 | 0.707 s | 3536 us | 3.60x |
| 200 | 400 | 4.449 s | 11123 us | 11.33x |
| 400 | 800 | **24.765 s** | 30956 us | **31.53x** |

**Fitted N^2.27.** Per element it is 31x worse at N=400 than at N=25. A 1600-node circuit
could not be built inside a ten-minute budget at all, which is why 7b's `n=2000` row was
killed rather than completed — **the solver work it was meant to measure never ran**.

**Two independent mechanisms, both named by the profile at N=400:**

| cost | calls | time |
|---|---|---|
| `list.index` | 642,400 | 27.6 s |
| `circuit.py:70 Node.__eq__` | **129,042,200** | 17.1 s |
| `update_node_map` | 800 | 3.9 s |

1. **`circuit.py:1144` does `self.nodes.index(node)`** — a linear scan of the node list per
   node per element, each step calling `Node.__eq__`. **`Node` already defines `__hash__`
   (`hash(self.name)`)**, so the nodes are hashable and a dict would make this O(1). The
   129 million `__eq__` calls are that scan; `__eq__` is a `try` / bare `except` around a
   name comparison, so each one also pays exception-handler setup.
2. **`update_node_map()` rebuilds the entire map from scratch**, and is called **once per
   element insertion** — 800 rebuilds for 800 elements. Inserting one element adds one
   entry; it should not re-derive the other 799.

Either fix alone helps; together they should make construction linear. They are listed
separately because they are separately testable and a combined change would not be
attributable.

**Gate 2+.5-1 (node indices are UNCHANGED).** The most important gate, and the reason this
is not a trivial change: `self.nodes.index(node)` produces the **matrix row index**, so any
reordering silently permutes every `G`, `C` and `J` in the package. Declared: for a set of
circuits including one with global nodes and one with subcircuit hierarchy, the full
`elementnodemap` must be **equal element-by-element** to today's. A dict built by insertion
order preserves this in Python 3.7+, but that must be asserted, not assumed.
OUTCOME: **PASS.** Identical `elementnodemap` on four shapes — an 8- and a 40-section ladder, a
mixed circuit with global nodes / inductor / diode / `TLine` / current source, and a nested
subcircuit. The suite test asserts the **invariant** (an element's rows are its nodes'
positions in `cir.nodes`, then its branch rows) rather than a recorded snapshot, and a
separate test pins **first-occurrence-wins**: `list.index` returns the first match while a
dict comprehension keeps the last, so two equal-named nodes would silently renumber matrix
rows.

**Gate 2+.5-2 (bit-identical results).** Stage 2's inherited stop condition: waveform drift
exactly `0.00e+00` and identical step counts on the full baseline, or stop. Construction
feeds indices only — nothing numerical should move, so unlike 2+.4 this bar is expected to
be reachable and a failure means the indices moved.
OUTCOME: **PASS.** All 17 baseline runs `max|diff| = 0.000e+00`.

**Gate 2+.5-3 (the scaling is actually fixed).** Declared *before* implementation: fitted
exponent **<= 1.3** over N = 25..400, and N=800 built in **under 5 s** — a case that is
untestable today at ~100 s. A speedup ratio alone is not the gate; the exponent is, because
a constant-factor win would leave the wall in the same place.
OUTCOME: **PASS on both clauses. N^1.15 against the declared <= 1.30**, and N=800 built in
**0.345 s** against a declared 5 s:

| N | before | after |
|---|---|---|
| 400 | 18.887 s | **0.174 s** (108x) |
| 800 | ~100 s | **0.345 s** |
| 1600 | *unbuildable in 10 min* | **1.25 s** |
| 3200 | — | 4.45 s (6400 elements) |

Two changes, both needed: the dict lookup alone took N^2.23 to N^1.97 (the per-insertion
full rebuild is what remained); deferring the rebuild took it to N^1.15.

**THE ITEM'S OWN RECONSIDER-IF WAS WRONG, and testing it was worth it.** It proposed fixing
`Node.__eq__` alone first — dropping the `try`/bare-`except`. Measured on its own that is a
**2x REGRESSION**: 18.9 s to 40.8 s at N=400, because `getattr(a, 'name', default)` is
*slower* than `try: a.name` when the attribute normally exists. Python's `try` costs almost
nothing unless it fires. Reverted.

**Gate 2+.5-4.** Full suite, and `test_circuit.py` in particular — the node map is what
everything else indexes through.
OUTCOME: **PASS. 865 passed, 6 skipped, 3 xfailed, 0 failed**, against 859 plus the six tests added
here.

**A PERFORMANCE REGRESSION WAS INTRODUCED AND CAUGHT.** Making `elementnodemap` a property
put a Python call in two per-element loops — `limit` (every Newton iteration) and
`accept_step` (every accepted step) — where the stamping paths already hoist. The first
suite run took **16m31s against a ~8m17s baseline**. Hoisting both restored it: a
representative transient measures **1.256 s against 1.287 s pre-change**.

**And the second slow run was NOT the same thing.** It came in at 14m42s with the hoist in
place, which contradicted the transient measurement. Rather than assume load, the two
versions were run back to back on the same box: **52.98 / 54.27 s against 52.80 / 53.86 s**
— identical. Load average was 5.27 from unrelated processes. Trap 2, established by
controlled comparison instead of asserted.

**Note for whoever takes 2+.4:** that item wants assembly to write into a reduced matrix
directly, which also touches `update_node_map`'s index arithmetic. **Do 2+.5 first** — it
is smaller, it is a strict prerequisite for measuring anything at n >= 1000, and 2+.4's
gates are unmeasurable until a circuit that size can be built.

**Reconsider if** profiling shows the `__eq__` cost dominates the `index` cost rather than
being caused by it — then the cheaper first move is `Node.__eq__` alone (drop the bare
`except`, compare `type` then `name`), which is a two-line change and would be worth
measuring before the dict work.

---


## DECISIONS TAKEN, 2026-08-01

### D1 — `LinearSolver.factor()` / 7b-2 factor reuse: **DEFERRED, not rejected**

Reusing Newton's last factorisation for the step controller's `J^-1 Eg` solve would remove
**33.1%** of all linear solves (measured: 228 Newton against 113 controller over 111
steps). It is deferred on arithmetic, not on principle:

- The solve is **~8% of a transient** at n~300, so removing a third of them saves
  **~2.7% of runtime**.
- It costs **bit-identity**, because it substitutes `J(x_k)` for `J(x_conv)` — median
  5.5e-9, max 2.2e-5, exactly zero on linear circuits where `J` depends only on `dt`.
- Bit-identity is what made every safe refactor of this session possible: 7a, 7b and 9(a)
  each passed on `max|diff| = 0.000e+00`, and 9(a) was a five-copy deduplication that
  nobody would have dared without it.

**2.7% is a poor price for that.** The seam exists and returns `None`, with the reasoning
on it.

**Reconsider when 2+.5 lands.** At n=5000, 7a's scaling table puts LU at **60%** of the
step rather than 8%, so the same change becomes **~20%** and the trade flips. It cannot be
evaluated before then, because circuits that size cannot currently be *built* — which is
2+.5.

### D2 — `dt_max`: keep `timestep`, add an explicit `max_step`

**Measured both ways before deciding.** On a run of `5 tau` the clamp is mostly irrelevant:
above `timestep ~ 3e-4` the **error controller** becomes binding and `max dt` stalls at
2.97e-4 however much slack the clamp is given. There the clamp only bites when a finer grid
was explicitly requested, which is the caller getting what they asked for.

On a run of `100 tau` — ~99% quiescent — it is clamped at **every** setting:

| timestep | steps | max dt | final err | |
|---|---|---|---|---|
| 1e-4 | **1027** | 1.000e-4 | 4.4e-16 | clamped |
| 1e-3 | 158 | 1.000e-3 | 0.0 | clamped |
| 1e-2 | 75 | 1.000e-2 | 2.6e-14 | clamped |

**1027 steps across a dead-flat solution whose error is 4e-16.** That is the real cost, and
it is paid by exactly the circuits that idle.

**Decision: `max_step` becomes a parameter defaulting to `None`, meaning `timestep`.**
Changing the *default* would move every waveform in the package for a benefit only
idle-heavy circuits see; making it *reachable* costs nothing and lets a caller ask for
SPICE's own `TMAX` behaviour.

**Gate D2-1 (nothing moves by default).** `max_step=None` must be bit-identical on the full
17-run baseline. Declared at bit-identical because this is a pure plumbing change when the
parameter is not set.
OUTCOME: **PASS.** All 17 baseline runs `max|diff| = 0.000e+00`, and a test asserts `max_step=None`
and `max_step=timestep` give the identical step count and identical `max dt`.

**Gate D2-2 (it actually frees the step).** On the `100 tau` quiescent run, `max_step`
larger than `timestep` must cut the step count substantially while the final error stays at
the same order. Declared: **>= 3x fewer steps at `max_step = tend/50`**, and final error no
worse than 10x. A change that frees the clamp without reducing steps would mean the
controller was never the thing being held back.
OUTCOME: **PASS, comfortably past the declared 3x.** On the `100 tau` quiescent run at
`timestep = 1e-4`:

| max_step | steps | max dt | final err | |
|---|---|---|---|---|
| None (= 1e-4) | 1027 | 1.000e-4 | 4.441e-16 | |
| 1e-3 | 161 | 1.000e-3 | **0.0** | 6.4x fewer |
| 2e-3 (SPICE's `tend/50`) | 116 | 2.000e-3 | **0.0** | **8.9x fewer** |

The error does not degrade — it *improves* to exactly zero, because the settled solution is
flat and fewer steps accumulate less round-off. The test asserts the step-count **ratio**
with the error held, not an absolute count, so it says "the clamp was what held it back"
rather than pinning today's number.

**Gate D2-3 (both backends, one meaning).** `Transient` and `JAXTransient` must take the
same parameter with the same default, since 9(d)/(g) exist because the two drifted. A
`max_step` below `timestep` is a caller error, not a silent clamp — declared to raise.
OUTCOME: **PASS.** `max_step` on `Transient` and `JAXTransient` with the same name and the same
`None` default, asserted by a test that compares the two parameter lists rather than
trusting that both were edited. Wired at all four clamp sites — `_solve`, `_solve_coupled`,
`JAXTransient.solve`, `solve_batched`. `max_step < timestep` raises on both, naming both
values: the step can never exceed `max_step`, so accepting one silently discards the
`timestep` the caller asked for.

**Gate D2-4.** Full suite.
OUTCOME: **PASS. Suite 844 passed, 6 skipped, 0 failed**, one summary line per trap 15, against 840
plus the four tests added here.

**Doc build re-run and verified by content:** `build succeeded, 2 warnings, 0 ERROR`, **0
degraded `exec-rst` blocks**, and both generated tables in `lte_dae.html` hold live numbers
rather than rendered source. This retires the "unverified since `d060b13`" label the resume
block was carrying.
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

## D3, second attempt — 2026-07-31, after 4g(b) and 4i

**Entry condition, checked before any change.** Gate D3-a failed on one property: under
`sigglobal` the trapezoidal error stopped falling monotonically as `reltol` tightened. That
was attributed to the `(-1)^n` mode in the trapezoidal estimator. Re-measured with
`benchmarks/transient_decisions.py --relref` on the current estimators:

| relref | integrator | steps (reltol 1e-3 .. 1e-6) | error |
|---|---|---|---|
| pointlocal | trap-classic | 54 / 56 / 86 / 152 | 1.03e-3 / 8.70e-4 / 2.05e-4 / 4.59e-5 |
| pointlocal | gear2-classic | 53 / 69 / 112 / 204 | 4.02e-3 / 2.01e-3 / 4.44e-4 / 9.76e-5 |
| **sigglobal** | **trap-classic** | 36 / 40 / 56 / 90 | 1.03e-3 / 1.01e-3 / 3.03e-4 / **7.02e-5** |
| **sigglobal** | **trap-ywr** | 36 / 41 / 56 / 89 | 1.03e-3 / 9.94e-4 / 3.03e-4 / **7.03e-5** |
| **sigglobal** | **gear2-classic** | 37 / 40 / 62 / 117 | 3.94e-3 / 2.91e-3 / 6.76e-4 / **1.51e-4** |
| **sigglobal** | euler | 38 / 68 / 190 / 584 | 2.35e-2 / 8.12e-3 / 2.63e-3 / 8.39e-4 |

**Every one of the six configurations is now monotone in both columns.** The blocker is
gone, and it is gone for the stated reason — the two trapezoidal rows are the ones that
used to rise (`2.23e-4 -> 3.14e-4` and `3.82e-5 -> 2.65e-4`) and now fall cleanly.

### A claim in this plan that the re-run does NOT support, and a gate to settle it

The recorded justification for `sigglobal` is "**1.7-2.5x fewer steps**". That comparison is
taken at equal `reltol`, and the table above shows the two modes do not deliver equal
accuracy at equal `reltol` — `sigglobal`'s error at 1e-6 is 1.5x larger for both trapezoidal
and Gear2. **A step count taken at a looser effective tolerance is not a speedup**, it is a
relabelling, and 4i has just produced one instance of exactly that mistake being avoided
only because it was checked. So the claim is not carried forward untested.

**Gate D3-c (NEW — the efficiency is real, not a relabelling).** Compare the steps each
mode needs to reach the **same global error**, not the same `reltol`. Declared success:
`sigglobal` needs at least **1.2x fewer steps at matched accuracy** for every integrator.
**If it does not, `sigglobal` is a change in what `reltol` means and must be described that
way in the parameter documentation**, and the "1.7-2.5x" figure is struck from this plan.
OUTCOME: **PASSED — the efficiency is real — and the recorded 1.7-2.5x is struck as an
overstatement.** Measured by sweeping `reltol` over nine values per mode and interpolating
the steps each needs to reach the *same* global error, inside the error range both modes
cover so neither is extrapolated:

| integrator | matched-accuracy step ratio (pointlocal / sigglobal) | worst |
|---|---|---|
| euler | 1.89x, 1.95x, 2.06x, 2.06x, 1.48x | **1.48x** |
| trap-classic | 1.47x, 1.42x, 1.38x, 1.31x, 1.35x | **1.31x** |
| gear2-classic | 1.44x, 1.47x, 1.54x, 1.60x, 1.55x | **1.44x** |
| gear2-ywr | 1.44x, 1.47x, 1.54x, 1.57x, 1.55x | **1.44x** |

Worst case **1.31x** against the declared 1.2x, so the gate passes and `sigglobal` is a
genuine efficiency gain rather than a relabelling of the tolerance. **But the honest figure
is 1.31-2.06x, not the "1.7-2.5x" recorded under gate D3-a.** That number was taken at equal
`reltol`, and the two modes do not deliver equal accuracy at equal `reltol` — `sigglobal`'s
error is ~1.5x larger there. **This is the fifth magnitude in this plan to shrink on
measurement**, and the first one this work produced itself; the comment in `transient.py`
carries the corrected range.

**Preserved as `benchmarks/transient_decisions.py --matched-accuracy`**, and preserving it
turned up one more thing. The script must **pin `lte_vabstol`**: this gate ran while the
floor was 1e-6, and D3-e then moved the shipped default to 1e-12. The floor is load-bearing
under `pointlocal` and inert under `sigglobal` — that is D3-e's whole finding — so leaving
it at the default penalises the `pointlocal` arm and inflates the very ratio this gate
reports. At the gate's own floor the range is **1.31x .. 2.06x** as recorded; at the shipped
floor it reads **1.70x .. 5.97x**.

That second number is not a better version of the first: **it is `relref` and `lte_vabstol`
together**, and it is reported separately for that reason. It does say that D3's two changes
compound rather than overlap — which was not measured at the time and is worth knowing if
anyone is deciding whether `pointlocal` is still worth offering.

**Gate D3-a (re-run).** Full suite `-m ""`. Expect churn in step counts; every failure
explained individually and either fixed or justified. Declared success: no failure that is
not explained, and specifically **no recurrence of the non-monotone accuracy** that failed
this gate the first time.
OUTCOME: **PASSED. 779 passed, 6 skipped, 0 failed, 670.77 s** (`-m "" --timeout=400`),
against 778 before — the one addition is the new `sigglobal` half of gate 1-5. The
non-monotone accuracy is gone: all six integrator/formula combinations fall monotonically
in both step count and error across reltol 1e-3..1e-6 (the table in the entry condition
above). Runtime 722.88 -> 670.77 s, a 7% fall, which is consistent with `sigglobal` needing
fewer steps but is **not** offered as a measurement of it — this box varies more than that
between runs (trap 2).

Doc build: **build succeeded, 2 warnings, 0 ERROR**, verified by content. The regenerated
tables move as expected under the new default (the stiff RLC at reltol 1e-3 goes 199 -> 60
steps for `gear2-ywr`) and **`ratios over bound` is still 0 on every row**, so 4b's
zero-stability invariant survives the change.

Three test failures during the work, each explained and each a stale assertion rather than
a regression:

1. **`test_gate_1_5_lte_vabstol_moves_the_step_count`** — predicted by gate D3-b. Rewritten;
   see gate D3-d.
2. **`test_lte_formula_ywr`** — asserted `rejections >= 1` at reltol 1e-4 as a proxy for
   "the controller is alive". Under `sigglobal` the tolerance does not collapse early in an
   RC charge where the node voltage is still small, so this circuit needs no rejection until
   1e-5. Measured, `gear2-classic` on that exact circuit:

   | relref | reltol 1e-3 .. 1e-7 | rejections | steps | error |
   |---|---|---|---|---|
   | pointlocal | | 3 / 3 / 2 / 2 / 2 | 31 / 52 / 97 / 192 / 387 | 9.75e-3 .. 1.97e-5 |
   | **sigglobal** | | **0 / 0 / 1 / 3 / 2** | 21 / 29 / 50 / 97 / 195 | 3.08e-2 .. 1.53e-4 |

   **The controller is demonstrably not blind** — step count and error are both monotone
   over five decades. The proxy was widened to "rejects somewhere in the sweep", which a
   blind controller still cannot satisfy and which does not depend on which mode ships.
3. **`test_gate_1_5_newton_and_lte_tolerances_are_separate`** — pinned `lte_vabstol == 1e-6`,
   which gate D3-e moved. Updated, and its docstring now says what actually proves
   separateness (see D3-d).

**Gate D3-d (NEW — gate 1-5 is rewritten, not deleted).** `sigglobal` makes the absolute
LTE floor non-load-bearing on a circuit with a healthy signal, which is the intended
behaviour and was measured under gate D3-b (403 steps at `lte_vabstol` 1e-3, 1e-6 and 1e-9
alike). `test_gate_1_5_lte_vabstol_moves_the_step_count` asserts the opposite and will fail.

**Deleting it is not acceptable**: gate 1-5 exists because `vabstol` used to serve two roles
at once, and the property that matters — *Newton's tolerance must not be a step-control
knob* — is still true and still worth pinning. Declared success: the rewritten test (i)
still fails if `vabstol` reaches the step controller, verified by injection rather than by
assertion alone; (ii) asserts a property that is true under the **shipped default**; and
(iii) the `pointlocal` behaviour it used to assert is kept as its own test, since that mode
still exists and is still selectable.
OUTCOME: **PASSED on all three clauses — and clause (i) caught something that would
otherwise have shipped.**

**The finding.** Clause (i) was written as a formality. Run properly — point the
controller's `abstol` at `self.par.vabstol`, i.e. reintroduce exactly the defect gate 1-5
exists to catch, and run the tests — and under `sigglobal`
**`test_gate_1_5_vabstol_does_not_move_the_step_count` still PASSES.** The guard is blind.

The mechanism is the same one that makes this whole item work: under `sigglobal` the
relative reference cannot degenerate, so `reltol*ref` dominates and *no* absolute term moves
anything — including an injected wrong one. **A guard that cannot see the defect it guards
is worse than no guard, because it reads as coverage.** It now runs under `pointlocal`,
where an absolute floor is load-bearing, and under injection it fails as it should. Verified
both ways: 4 passed clean, 2 failed under injection.

This is the second time in three commits that a test kept passing for a reason unrelated to
the property it names — the first was `test_lte_formula_ywr`'s rejection proxy above. Both
were found by asking "what would make this fail?" rather than by the suite.

Resulting shape, all four passing: `..._newton_and_lte_tolerances_are_separate` pins the
values; `..._lte_vabstol_moves_the_step_count_under_pointlocal` keeps the old property in
the mode where it holds; `..._lte_vabstol_is_not_load_bearing_under_sigglobal` asserts the
default's behaviour, which is what licenses D3-e; and `..._vabstol_does_not_move_the_step_count`
is the guard, now able to see.

**Gate D3-e (NEW — `lte_vabstol` is decided by measurement, not by tidiness).** With
`sigglobal` shipped, D3-b's finding says the absolute floor can return to 1e-12 "at little
cost". *Little* is not a number. Declared: measure step count and global error at
`lte_vabstol` 1e-6 (current) and 1e-12 (the defensible value matching `vabstol`), under
`sigglobal`, on both a healthy-signal circuit and the quiet-node circuit that motivated the
floor in the first place. Change the default **only if** the cost is under 10% of step count
on both; otherwise record the number and leave it, because a floor that is never reached is
harmless and one that collapses the step size is not.
OUTCOME: **PASSED, and the cost is not "little" — it is exactly zero.** `lte_vabstol` is
returned to **1e-12**.

| circuit | relref | 1e-6 | 1e-9 | 1e-12 |
|---|---|---|---|---|
| pulsed RC | **sigglobal** | 403 steps | **403** | **403** |
| quiet node | **sigglobal** | 601 steps | **601** | **601** |
| pulsed RC | pointlocal | 1164 steps | 1260 (+8.2%) | 1263 (**+8.5%**) |
| quiet node | pointlocal | 2463 steps | 2685 (+9.0%) | 2689 (**+9.2%**) |

Bit-identical under `sigglobal` — same step count and same endpoint voltage to six figures
at all three values — against a declared 10% threshold. The `pointlocal` rows are the
control: they show the floor *is* load-bearing there, at 8.5-9.2%, which is what the 1e-6
workaround was buying and why it was needed then.

**What this actually retires.** The 1e-6 value was introduced to stop the timestep
collapsing on the leapfrog's quiet nodes — treating the symptom of `pointlocal`'s
degenerating reference. `sigglobal` removes the cause, so the workaround can go and the
principled value (matching `vabstol`, so the operating point and every step after it use the
same absolute floor) returns at measured zero cost. **The parameter comment states plainly
that selecting `relref='pointlocal'` makes the floor load-bearing again**, so anyone who
does is not left with a value tuned for a different mode.
