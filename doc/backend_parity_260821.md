# CPU/JAX transient feature-parity review — 2026-08-21

**Branch:** `cna-jax-vectorization` (at `3ac9656`)
**Scope:** the *feature surface* of `Transient` vs `JAXTransient` — parameters, entry
points, selectable machinery, and semantics of shared knobs. This complements the
cross-backend conformance harness (`test_backend_conformance.py`), which guards
*behavioral* parity of what both backends already do; this document is about what only
one of them does, and whether that should stay true.
**Method:** mechanical diff of the declared `Parameter` lists and `solve` signatures,
plus source verification of every claimed divergence. New findings are marked; every
gap carries a verdict, an effort size, and — where closure is recommended — the test
that would pin it.

**Goal, as stated:** the backends should be feature-comparable *where possible*. The
verdicts below therefore split three ways: **align** (close the gap), **document**
(deliberate divergence, stated where users look), and **one-sided by design** (the
feature is the backend's reason to exist, or cannot survive a traced loop).

---

## Summary matrix

Parameters and options, from the mechanical diff:

⚠ **READ THE VERDICT COLUMN AS "WHAT WAS DECIDED AT REVIEW TIME", NOT AS STATUS.**
Every row was written before any of it was done, and **P1–P25 have since all been
executed** (Phases A/B/C; the narrative sections below carry each one's outcome). A row
reading "align now (S)" therefore describes a decision that was taken and carried out, not
work outstanding. That is not a hypothetical misreading: on 2026-08-31 this table's P5 and
P13 rows were cited elsewhere as "architecture-review residuals" still open, and both had
been done for weeks. The rows are kept as written because the verdicts are the record;
this note is the correction.

| Feature | CPU `Transient` | `JAXTransient` | Verdict (as decided, 2026-08-21) |
|---|---|---|---|
| `reltol`, `iabstol`, `vabstol`, `lte_vabstol`, `lte_iabstol`, `maxiter`, `timestep_max`, `firststep`, `epar` | ✓ | ✓ | aligned (post-review; `max_step` renamed `timestep_max` and DECOUPLED from `timestep` by owner decision 2026-08-21 — None now means tend/50, SPICE's TMAX, so tolerance knobs are live on gentle circuits where the old timestep-as-cap froze the step count) |
| t=0 point, tend landing, breakpoint order drop, force-accept accounting, per-row Newton, safety-factor law, statistics-on-result | ✓ | ✓ | aligned (Phases 1–3) |
| **P1** `uic` / `minstep` / typo safety | Parameters | `**kwargs`, silently swallowed | **align now (S)** |
| **P2** LTE-ratio knob | `LTERATIO = 7.0` class constant | `TRTOL` Parameter | **align now (S)** — reversed asymmetry |
| **P3** `solve_batched` floor & naming | — | `dt_min=1e-15` vs solve's `1e-18`; arg `irefnode` vs `refnode` | **align now (S)** |
| **P4** `refnode` default | `gnd` | `0` | **align now (S)** |
| **P5** `analysis` parameter | `'tran'`, threaded | inherited, dead; loop hardcodes `'tran'` | **align now (S)** — **DONE**, verified 2026-08-31: `Parameter(name='analysis', default='tran')` declared and threaded as `self.par.analysis` at both `outer_time_loop` call sites, and the prescribed test exists — the dead-knob allowlist entry for `('JAXTransient', 'analysis')` is deleted (`test_backend_parity_phaseA.py`) |
| **P6** integrator selection | `integrator=` Parameter (Euler/Trap/Gear2) | internal `'euler'/'gear'`, not exposed; hardcoded `'gear'` | **align (M)** + a default-order decision |
| **P7** `relref` modes | `pointlocal`/`alllocal`/`sigglobal` | `sigglobal` only, hardwired | align later (M) |
| **P8** LTE band (`lte_gamma_*`, `lte_eta`) | ✓ (standard + coupled) | absent; accept test hardcodes `err <= 1` | align later (M) |
| **P9** `fixed_timestep` | ✓ | absent | align later (M) |
| **P10** `outputstep` uniform-grid resample | ✓ | absent — but the resampler is backend-agnostic | **align now (S)** |
| **P11** `provided_function` extra source | ✓ (`pf(t)`) | absent | align later (M, with a traceability constraint) |
| **P12** initial conditions under `uic` | `ic` dict + element ICs + spanning-tree solve | reads `node.ic` attributes — **dead code, nothing sets them** (new finding) | **align (M)** — delete the dead path, reuse the CPU machinery |
| **P13** `bypass`/`bypasstol` | ✓ | absent | document — meaningless under vmap — **DONE**, verified 2026-08-31: the line is in the `JAXTransient` docstring ("a non-concept here, not a missing feature") |
| **P14** `minbreak` | ✓ | progress guard only, no spacing merge | align now (S) |
| **P15** statistics vocabulary | `force_accepts`, `order_drops`, `breakpoints_hit`, timing | `forced_lte_steps`, `nonconverged_steps`, `signal_max` | partial align (S) |
| **P16** TLine delay interpolation | quadratic (3-point Lagrange) | linear + ring-cap raise; **stamp defect fixed + step cap wired + std/pcnr paths verified to 5e-16 / 4.4e-10** (TLine port); **quadratic upgrade EXECUTED 2026-08-21**, Phase B -- `interp_tlines` carries the CPU's 3-point stencil with a guarded linear fallback (`jaxtransient.py`). This cell read "still open" until the 2026-08-31 audit, four lines above the note recording it done | **DONE** |
| **P17** solver strategy objects (`nrsolver`/`scaler`/`linearsolver`) | ✓, injectable | refused loudly | **document — permanent** (traced loop cannot dispatch) |
| **P18** continuation (gmin/source stepping) | DC path; transient retries dt | dt retries only | roadmap, both sides |
| **P19** `coupled_lte` (Fang), `pcnr` | ✓ | ✓ **PORTED** — 'approx' Fang + PCNR + PCNR-inside-Fang (`bordered` and TLine-under-coupled/pcnr stay CPU-only, refused loudly) | closed (see P19 note) |
| **P20** `solve_batched` per-lane sweeps | — | ✓ — the backend's reason to exist | one-sided by design |
| **P21** batched DC operating point | n/a | missing (`uic=True` required, refused loudly) | roadmap |

---

## Close-now items (Phase A — small, mostly mechanical)

### P1 — `uic`/`minstep` as Parameters; kill the `**kwargs` sink *(new finding)*

`JAXTransient.solve(self, refnode=0, tend, x0, timestep, CHUNK_SIZE, **kwargs)` reads
`kwargs.get('uic', False)` and `kwargs.get('minstep', 1e-18)` and ignores everything
else — so `solve(..., uic_=True)` or `solve(..., minstep=...)` misspelled runs
silently with defaults. That is the dead-knob defect class at the call boundary, in
the backend whose review just spent five findings on it. The CPU side declares both
as Parameters.

**Fix:** declare `uic` and `minstep` as Parameters (same names, same defaults as CPU);
`solve` takes no `**kwargs`, so a typo raises `TypeError`. `CHUNK_SIZE` stays a solve
argument (it is an execution detail, not a model setting) but gets the same treatment
as the rest of the signature: explicit.
**Test:** `solve(uic=True)` works via the Parameter; `solve(uicc=True)` raises;
the dead-knob scan's kwarg half covers regressions.

### P2 — one name for the LTE ratio, on both sides

The asymmetry is *reversed*: JAX declares `TRTOL` as a Parameter while the CPU
hardcodes `LTERATIO = 7.0` as a class constant (recently consolidated to one
definition, but still not settable). A user tuning Spectre's `lteratio` can do it on
one backend only — the backend with fewer knobs elsewhere.

**Fix:** CPU `Transient` gains the same `TRTOL` Parameter (default 7.0) feeding
`self.LTERATIO`'s uses; keep the class attribute as the default's home or retire it
into the Parameter. Same name both sides — `TRTOL`, since it is already public API
on JAX and matches the SPICE literature.
**Test:** conformance-harness variant at `TRTOL=1` on both backends: step counts move
the same direction; parameter-reachability test updated.

### P3 — one `dt_min` floor and one argument vocabulary for `solve_batched`

`solve_batched(..., dt_min=1e-15, ...)` vs `solve`'s `minstep → 1e-18`: the same
physical floor differs by three decades between two entry points of one class, and
the argument is named `dt_min` here, `minstep` there, `minstep` (Parameter) on the
CPU. Also `solve_batched`'s first argument is `irefnode` but accepts a node object
(it calls `get_node_index` on it) — the name lies.

**Fix:** `solve_batched` reads the same `minstep` Parameter (P1) as its floor;
`irefnode` renamed `refnode`. Keep `dt_min`/old name accepted for one release only if
external callers exist — the repo's own callers are the tests.
**Test:** both entry points report the same effective floor on a probe run.

### P4 — `refnode=gnd` default on both

JAX defaults `refnode=0` (works because `get_node_index(0)`… resolves gnd's index by
accident of ordering); CPU defaults `gnd`. Same default object on both removes a
class of silent wrong-reference runs on reordered circuits.
**Fix + test:** one line; assert the pinned row is gnd's on a circuit built with
nodes added before gnd would sort first.

### P5 — `analysis` threaded, not hardcoded

The JAX loop hardcodes `analysis='tran'` at its `circuit.u` call while the inherited
`analysis` Parameter sits dead (currently allowlisted in the dead-knob scan). Set the
declared default to `'tran'` and thread `self.par.analysis` through `outer_time_loop`
— then delete the allowlist entry, which is the test.

### P10 — `outputstep` on the JAX path

`resample_uniform` is a free function on arrays; the CPU path applies it when
`outputstep` is set. Nothing about it is backend-specific — the JAX result is numpy
by the time it exists. Declare the same Parameter, apply the same call before
returning. One test: uniform grid out, quadratic values within tolerance of the
adaptive run.

### P14 — `minbreak` honoured by `collect_breakpoints`

The CPU merges breakpoints closer than `minbreak`; the JAX scan has a progress guard
and a count cap but no spacing merge, so two sources with edges 1e-16 apart produce
two breakpoints and a degenerate truncated step. Apply the same `minbreak`-scaled
merge in `collect_breakpoints` (Parameter added, same default 1e-14).
**Test:** two pulses offset by less than `minbreak` yield one merged breakpoint.

### P15 — statistics vocabulary (partial)

Semantically shared counters should share names: JAX's `forced_lte_steps` **is** the
CPU's `force_accepts` (same definition: accepted with LTE above tolerance) — rename
to `force_accepts`. `rejected_steps`/`accepted_steps` already match. Genuinely
backend-specific fields keep their names and a docstring saying so
(`nonconverged_steps`, `forced_nonconverged_steps`, `signal_max` on JAX;
`newton_iterations`, `order_drops`, `breakpoints_hit`, timing on CPU — the last four
are also candidates for later port, cheap ones first: `breakpoints_hit` and
`order_drops` are one counter each in the traced state).

---

## Align-later items (Phase B — real work, worth scheduling)

### P6 — expose the integrator choice, and decide the shipped default

The JAX loop already implements `'euler'` and `'gear'`; `eval_method` is simply not
reachable from the API (hardcoded at both call sites). Exposing it as
`integrator='gear'|'euler'` (Parameter) is mostly threading. Trapezoidal stays
CPU-only until someone ports the *variable-step* estimator — the uniform-grid trap
branch was deleted for cause; say so in the Parameter doc.

**DONE 2026-09-01, and the framing above was wrong in a way worth keeping.**
"Port the variable-step estimator" implied the mathematics was missing. It was
not: `trapezoidal_lte` and `third_divided_difference` were already in
`_lte_kernels`, already used by the CPU, already plain arithmetic that traces.
What was missing was the WIRING — a much smaller claim, and the gap between the
two is roughly the day of work nobody did because the record said otherwise.

⚠ The one part that genuinely had to be written, and could not be copied: the
estimator must difference the **charge**, not the companion current. The Gear-2
branch beside it differences `g`; for trapezoidal that is forbidden, because the
recursion `iq_n = 2(q_n − q_{n−1})/h − iq_{n−1}` carries an undamped `(−1)^n`
mode and differencing `g` measures that mode (a 1.9x swing from step history
alone). So `ywr_error_ratio` gained a `q_curr` argument. Copying the neighbouring
branch — the obvious move — would have produced a plausible, passing, wrong
estimator; the deletion of the old branch is what forced the question to be asked.

Measured at landing (rc/vsin, reltol 1e-4, ~5x bounds per the conformance
harness's rule): **rc 4.2779e-06, vsin 4.1643e-03**, step counts 61/60 and 92/89
CPU/JAX; Gear-2 measures rc 4.3e-6 / vsin 2.9e-3 for reference. The decisive gate
is the stiff mode: at `h·λ = -1000` both backends ring at **0.9960**, the textbook
`|R_TR(-1000)| = 0.996008`. That is what proves the branch is really trapezoidal
— an estimator that had silently fallen back to Euler would still have passed the
smooth pairings. The ringing is unmitigated on both backends and must stay that
way; damping the JAX side would be a parity break dressed as an improvement.

This also forces the deferred F19(b) decision into the open: **the shipped defaults
differ** (CPU Euler, JAX Gear-2), so identical scripts get different methods by
backend. Recommendation: make **Gear-2 the CPU default too** — the Phase-0 baseline
measured it at half the steps and half the wall-clock of Euler on the same circuit at
the same tolerance, the estimator work of stages 4g/4i was built for it, and the
conformance harness already treats it as the reference configuration. This is a
user-facing waveform change on the CPU path (every default run's step sequence
moves), so it is an owner decision; the fallback is documenting the asymmetry in both
`integrator` Parameter descriptions. Either way the harness pins whichever is chosen.

> **EXECUTED 2026-08-21** (owner decision: Gear-2 on both).  The JAX
> `integrator` Parameter ('gear'|'euler', trap refused with the reason) reaches
> the traced loop; order is measurably live (euler 8.2e-3 vs gear 1.1e-3 on the
> step-response RC).  One finding worth its own line: **the CPU coupled path
> keeps Euler as its own default** — Fang's step law, the acceptance-band
> record, and the 127/127 cross-backend step parity were all derived at
> order 1, and Gear-2 underneath the coupled path was measured to livelock the
> coupled rectifier (h collapsed to 6.3e-12 at t=1.25e-4) and break the parity
> (96 vs 127).  Explicit integrators are honoured on either path; re-deriving
> the coupled estimator at order 2 is the recorded reconsider-if.  Four
> Euler-era method-calibrated records (two QUCS fixed-grid endpoints, the
> Idtmod wrap sample, the BJT storage-time triple) now pin
> `integrator=EulerIntegrator()` explicitly.

### P7 — `relref` modes on JAX

`sigglobal` is carried in traced state (`sig_max`). `pointlocal` is stateless
(`max(|x|, |x_last|)` per row — trivial). `alllocal` needs a per-row running-max
*vector* in `TransientState` — same shape as `x`, same update discipline as
`sig_max`. All three are traceable; the work is threading a mode switch that is
static at trace time (Python string, like `eval_method`). Port the CPU's
unit-group split (`n_nodes`) with it, or `sigglobal` on JAX keeps mixing volts and
amps — check whether it already does (the CPU comment warns about exactly this; the
JAX `ref` is a scalar max over everything, so **it does mix units today** — worth a
line in the port, and until then a note in the code).


> **EXECUTED 2026-08-21.**  All three modes on the traced standard estimator
> (`ref_running` per-row vector in state, accepted-only updates; the mode is a
> trace-static string like `eval_method`), WITH the unit-group split — the
> kernel gate pins that a 1 A branch current no longer deadens millivolt node
> control (ratio restored >100×).  Measured mode ordering on the step-response
> RC at reltol 1e-5: pointlocal 95 steps / 3.4e-4, alllocal 72 / 4.9e-4,
> sigglobal 63 / 6.9e-4.  The COUPLED path deliberately keeps its scalar
> reference — its whole measured record (the TLine kink constants included) is
> derived against it; noted at the site.

### P8 — the LTE acceptance band on JAX

Post-F17 both backends share the prediction law, so the band port is the acceptance
test only: `gamma_min <= err <= gamma_max` with the growth-retry (reject with a
larger dt — the machinery exists: `do_reject` with `retry_dt` can grow) and the eta
damper as a clamp on `factor`. Reuse the CPU's `'auto'` sentinel semantics from F5.
Moderate; do after P6 so the band composes with a chosen integrator.

> **EXECUTED 2026-08-21.**  `_standard_band()` mirrors set_lte_band's 'auto'
> resolution (0, 1, None) and validation; the traced accept takes gamma_max,
> the too-accurate growth-reject carries the CPU's three guards (breakpoint
> landing, at-cap, would-not-actually-grow), and eta clamps the factor
> branchlessly (inf when unset).  Measured: defaults bit-identical to the
> pre-band loop; gamma_min=0.3 fires 10 growth-rejects (63 → 58 accepted);
> eta=0.10 bounds the step ratio at exactly 1.100; invalid bands refused.

### P9 — `fixed_timestep` on JAX

Semantics per the CPU: grid wins, breakpoints don't move it, order still drops when
a step crosses an edge, non-uniform fallback warns. In the traced loop this is
`dt_next = timestep` (skip `calculate_next_dt`) plus the accept test forced true
except on Newton failure. Moderate, mostly tests.

> **EXECUTED 2026-08-21**, and the port found a CPU defect: the fixed-grid
> crossing test `next_t_break <= t + dt` is a float knife-edge under
> accumulated `t += dt` — spied on the pulsed-RC run, the order drop fired at
> the FIRST edge and silently missed every later one.  Toleranced on the CPU
> (`dt·1e-9`), after which the two backends produce the SAME fixed-grid
> waveform to 5.6e-17.  The traced side: grid pinned (breakpoints untouched,
> tend landed exactly), LTE skipped, order drop by a toleranced crossing
> test, non-uniform fallback warned once post-run with the count (a trace
> cannot warn per event), element-cap warning mirrored, coupled+fixed refused
> (the CPU's grid_locked wiring is not ported).  Two CPU behaviors recorded
> as findings, not replicated: after a fixed-grid non-convergence shrink the
> CPU keeps the shrunken step for the rest of the run (its own warn text says
> "for this step"); JAX returns to the grid on the next accept.

### P11 — `provided_function` on JAX, with the honest constraint

The CPU contract (F4) is an extra source term `pf(t)`. In the traced loop, `pf` must
be jax-traceable (pure, jnp-composable) and is baked in at jit time. That constraint
is acceptable — state it in the Parameter doc and raise a clear error if tracing
fails. The DC-seed caveat from F4 applies verbatim (JAX's non-uic seed is a plain
DC solve that won't see `pf`): same warning, same wording.


> **EXECUTED 2026-08-21**, and the port found a latent JAX defect: the traced
> Newton scored its residual test over ALL rows including the reference row,
> which is not an equation of the solved system — an unbalanced pf put the
> full injection there, and the run livelocked at maxiter on every step,
> silently force-accepting at the dt floor (~1e14 steps to finish).  The
> convergence test now scores the reduced rows, as the CPU's always has.  pf
> folds in at every `circuit.u` site, so the coupled FD du/dt pair carries
> pf' automatically; the analytic-dudt branch omits it exactly as the CPU's
> residual_dh does (shared recorded gap).  DC-seed warning mirrored verbatim;
> cross-backend pf waveforms agree within 2% of signal on the gate.

### P12 — real initial conditions on JAX; delete the dead `node.ic` path *(new finding)*

`solve(uic=True)` walks `cir.nodes` reading a `node.ic` attribute **that nothing in
the package or tests ever sets** — dead code posing as a feature, adjacent to the
dead-knob class. Meanwhile the CPU's `_initial_state` (ic dict → element ICs →
spanning-tree capacitor solve) is pure pre-loop Python operating on names and
indices: it does not touch the toolkit until it returns a vector, so it can be
shared as-is. Extract it (module function or mixin), give JAX the same `ic`
Parameter, delete the `node.ic` walk.
**Test:** the CPU's `test_initial_conditions.py` cases, run through `JAXTransient`
— the spanning-tree tests included.

> **EXECUTED 2026-08-21.**  Shared by BINDING rather than extraction: the four
> CPU methods are pure pre-loop Python once `_initial_state` builds its vector
> with numpy instead of `self.toolkit` (a traced array cannot take item
> assignment), so `JAXTransient` assigns them as class attributes unchanged.
> The `ic` Parameter, the ic-without-uic refusal, and the node.ic deletion all
> landed; `solve_batched` seeds every lane from the same resolved vector.
> Gates in `test_jax_initial_conditions.py` (spanning tree, floating-group
> refusal, dead-walk regression included).

### P16 — TLine interpolation order

CPU delay lookup is 3-point quadratic with recorded reasoning ("a first-order
interpolant feeding a second-order method injects error the method never made");
JAX's is linear. Port the quadratic stencil into `interp_tlines` (three gathered
points instead of two; same modulo indexing). Measure the edge-accuracy delta on the
existing TLine tests before/after, per house rules.

> **EXECUTED 2026-08-21** — Phase B complete.  The ring lookup is the CPU's
> 3-point stencil, monotone-limited per channel exactly like the CPU's
> `_interpolate_history` (the phantom-EMF defect class), value and slope
> selected together so `tline_source_dudt` stays d/dt of the source.
> Measured per house rules: parabola recovered exactly (linear cannot), kink
> clamped; matched-line coupled agreement stays bit-close (2.8e-16); the
> rc-load corner gap did NOT move (2.66e-2 before and after) — that gap is
> controller step placement at the τ=66 ps arrival corner, not interpolation
> order, recorded so nobody re-chases it here.

---

## Documented divergences (no code change; one paragraph each where users look)

> **PHASE C EXECUTED 2026-08-21**: both class docstrings now carry the parity
> story (`Transient` — shared vocabulary/defaults, CPU-only-with-cause list,
> P20's why-no-imitation; `JAXTransient` — the aligned-parameter list, P13's
> non-concept line, P17's permanent refusal, P20's purpose statement, the
> CPU-only list).  P17's refusal is stated as PERMANENT in the code and now
> covers `linearsolver` too — it was still being silently swallowed; gated in
> test_backend_parity_phaseA.py.  A stale "Time step is fixed." line in the
> Transient docstring fell to the same edit.

- **P13 `bypass`:** device bypassing skips evaluating quiescent elements; a vmapped
  group evaluates all lanes of all instances in one kernel — there is nothing to
  skip. Not a missing feature; a non-concept on this execution model. One line in
  the JAX class docstring.
- **P17 strategy objects:** `nrsolver`/`scaler`/`linearsolver` are Python objects
  dispatched per iteration; a traced `while_loop` cannot call into them. The loud
  refusal in `__init__` **is** the contract — permanent, and now stated as such
  rather than as "not yet".
- **P19 research paths — verdict AMENDED (owner decision, same day): portable, and
  scheduled.** Feasibility review reversed the original call on three facts: Fang's
  solution-space LTE kernels are trace-ready by the stage-9(a) charter (one bisection
  helper needs a branchless rewrite); `Diode.eval_i_pure` already exists, so junction
  devices already evaluate on the JAX path; and PCNR is architecturally a *better* fit
  for a traced loop than the classic limiting it replaces — limiting mutates per-device
  Python state (`_vlim`), which a `lax.while_loop` cannot host, while PCNR's junction
  unknowns are just more state vector.  Plan: Fang's 'approx' branch first (M/L —
  per-lane coupled (x,h) in `solve_batched` is the payoff), PCNR second (L — framed as
  the JAX backend's nonlinear-robustness strategy, since that path has neither limiting
  nor PCNR today).  The CPU's +60–80 %/iteration PCNR cost figure does NOT transfer to
  vmapped execution and is re-measured, not assumed.  `bordered` and PCNR-inside-Fang
  stay CPU-only until separately justified.  *(PCNR-inside-Fang did not stay
  CPU-only: it landed the same day — see below.  A source comment beside
  `fang_inner_loop` went on asserting it had not, for weeks, and was corrected
  2026-09-01.)*

  **Executed same day** (commits `df28b7c`, `a9a0082`, `e591bca`): branchless
  bisection kernel; Fang 'approx' traced — rc step parity 127/127 CPU/JAX,
  JAX-coupled *more accurate* than CPU-coupled at every tolerance probed, batched
  per-lane coupled (x, h) landing both lanes on their analytic curves; PCNR traced —
  the cold-start junction probe fails outright on plain Newton and completes under
  PCNR (35 steps, 4.367 V), JAX-PCNR tracks CPU-PCNR to 1.1e-2, cost **+29 % wall**
  against the CPU's +60–80 %/iteration (the figure indeed did not transfer).  The
  dead-knob scan caught the Fang draft's unread TLine buffers mid-port, and the traced
  `pnjlim` is bit-exact against the branching original across a 48-case
  warnings-as-errors sweep.  **PCNR-inside-Fang followed the same day**: the cold-start probe that kills plain coupled completes under
  coupled+PCNR (56 steps, 4.367 V), and the rectifier tracks CPU coupled+PCNR to
  9.7e-3 at 7.4× fewer steps.  **The TLine follow-up** then found and fixed the
  standard path's DC-stamp defect (a matched 1 V line had returned 24.5 V since the
  JAX TLine landed, untested), wired the never-applied TD/2 step cap, and enabled
  TLine on the standard (5e-16 vs CPU) and PCNR (4.4e-10) paths.  **coupled+TLine
  landed last**, after a per-iteration trace pinned the livelock: at a from-zero
  kink every signal scales together, so the coupled LTE's dev ∝ ref and h cancels —
  err = 1/(TRTOL·reltol) = 1428.6 measured, h-independent — and the step after that
  extrapolates degree-2 through a pre-kink point (err = 139.2e11·h, in-band only
  below the within-point floor).  The failure point was the SOURCE EDGE, not the
  wavefront as first recorded.  Fix, three composed pieces: the coupled accept path
  zeroes `h_history` on a breakpoint landing (cold-start semantics: band skipped
  one step, degree held at 1 for two), `no_hist` also honours `force_first_order`,
  and `collect_breakpoints` registers source-corner + k·TD wavefront arrivals — the
  arrivals fix that was falsified alone works with the history reset, which is why
  it was reinstated.  Result: the pulsed matched line lands on every tend probed
  (1.5/2.5/8 ns) at 5e-16 vs the CPU standard path; an RC-loaded line at 2.6e-2 vs
  CPU-Gear2 (cross-method).  **CPU coupled+TLine landed the next day**, closing the
  last hole: `TLine.dudt` written (derivative of the same interpolation polynomial
  `u` evaluates, FD-verified to 0.0), the JAX kink discipline ported (step ring
  emptied on breakpoint landings, source corners echoed as corner + k·TD wavefront
  arrivals, and the coupled solve's growth capped at the breakpoint — the entry-h
  truncation test cleared corners that the solved h then straddled), and a defect
  the port EXPOSED fixed at its root: the quadratic history interpolation
  overshoots across a recorded kink (reflected EMF read 1.009 against samples
  bounded by 1.000), a band-blind step accepted a solution against the phantom,
  and the LTE hit an h-independent floor of exactly the pollution.  The
  interpolation is now monotone-limited per channel (quadratic outside the
  bracket envelope falls back to linear), which also makes it indifferent to
  accept_step's pruning.  Measured: matched line 5.6e-16 vs CPU-Gear2 and
  **1.7e-16 vs JAX-coupled** (bit-identical waveforms across backends);
  mismatched RC load completes to the correct Γ=1/3 steady level.  Fang
  coupled + TLine now runs on BOTH backends.
- **P20 `solve_batched`:** JAX-only by design — it is the branch's purpose. The CPU
  gets no imitation API; a Python loop over `Transient` is already expressible and
  honest about its cost.  **Measured 2026-08-21** (doc/batched_sweep_260821.md):
  on the rectifier corner sweep the crossover is ~16 lanes and the payoff 22.5×
  at 512 lanes (54.7 s CPU loop vs 2.4 s warm batched), with the JAX wall nearly
  flat in N — the purpose statement now carries its number.

## Roadmap items (neither backend has it; noting so silence isn't read as parity)

- **P18** continuation at a failed solve — the review's "deliberately
  unscheduled" item, **re-scoped 2026-08-21 (owner discussion)** with the
  industry vocabulary kept straight: **`gmin`** is the junction-parallel
  conductance (across each pn junction, bounding its cutoff conductance —
  the deformed circuit is a physical leaky diode, the linear subnetwork is
  untouched, and the homotopy tracks the physical branch); **`gshunt`** is
  node-to-ground (SPICE3's dynamic diagonal stepping), cruder but the only
  one that rescues a singular G from floating nodes.  Design: a
  junction-`gmin` ladder as the primary continuation, built on the
  machinery that already exists — `pcnr_junctions(cir)` enumerates the
  (node_a, node_b) pairs and PCNR's `_scatter_junction_G` scatters exactly
  the needed 4-entry pattern; F-side term ±gmin·(v_a − v_b) — with a
  `gshunt` rung reserved as the singular-matrix fallback.  Ramp in decades
  to **exactly zero**, and only a gmin-free converged solution may be
  accepted: residue in a committed point is the P22 inconsistency-floor
  trap by its measured signature.  Branch-tracking caveat (bistable
  circuits) recorded.  Order: batched DC first (converts the P21
  junction-lane raise into a rescue; static ramp schedule around the
  existing traced Newton, vmapped per lane for free), CPU DC second, the
  failed-time-point variant last (history rings must never see ladder-era
  iterates).

  > **EXECUTED 2026-08-21**, all three phases, with one scope finding:
  >
  > *Batched DC* — `dc_with_continuation`: plain Newton, then the
  > junction-gmin ladder (decades 1e-2 → 1e-12 → **0**), then the gshunt
  > ladder, `lax.cond`-skipped when plain lands; vmapped per lane.  The
  > P21 junction-lane raise became a rescue: the rectifier lanes that
  > raised yesterday now start at their true per-lane bias (4.366857 /
  > 4.384664 vs CPU-DC-with-limiting 4.366854 / 4.384664) — and the raise
  > survives where it must (a floating node with no DC path is singular at
  > gshunt=0; gated).
  >
  > *CPU DC* — the chain already existed (`GminSteppingNewton` +
  > `SourceSteppingNewton`, final pure solve and all); the P18 deltas were
  > vocabulary (its "Gmin" is **gshunt** by the owner's correction —
  > documented in the class, name kept for compatibility) and the missing
  > PRIMARY stage: `JunctionGminSteppingNewton` (proper `gmin`, rows from
  > `pcnr_junctions`, reduced-system indices) now runs before it.  The
  > unit gate walks the classic junction slam (pure first step 2.5e7,
  > limexp-clamped model) through the ladder to the exact solution with a
  > pure final Jacobian asserted.
  >
  > *Transient failed point* — the standard CPU path's minstep exhaust now
  > attempts ONE chain-wrapped re-solve (junction-gmin → gshunt, no source
  > stepping: scaling u(t) mid-transient would scale the companion
  > history too) before raising; only a pure converged point flows into
  > the accept machinery, counted in `statistics.gmin_rescues`.  SCOPE
  > FINDING, recorded honestly: no legitimate triggering circuit could be
  > fabricated — the CPU's limiting + retry machinery genuinely covers the
  > space (which is why the review left P18 unscheduled), and the one
  > fabrication attempt (neutralizing `cir.limit`) invalidated itself by
  > freezing `Diode._vlim`, which the model's evaluation depends on.  The
  > mechanism is unit-gated; coupled-path and JAX-transient point rescues
  > are recorded reconsider-ifs (the latter would put the ladder inside
  > the traced step).
- **P21** batched DC operating point, so `solve_batched` can start from bias. The
  loud `uic=True` refusal stays the contract until then.

  > **EXECUTED 2026-08-21.**  `dc_operating_point` is the CPU DC class's
  > assembly (`F = i(x) + u(0, 'dc')`, `J = G`) as a traced Newton on the
  > reduced system (convergence scored on reduced rows — the P11 lesson),
  > `jax.vmap`-ed so every lane's bias is solved with ITS OWN swept
  > parameters.  Measured: a 2 V divider with R1 swept {1k, 3k} starts its
  > lanes at exactly 1.0 V / 0.5 V (analytic) and holds them flat.  The
  > honest failure contract replaces the old refusal: no limiting or PCNR
  > at batched DC yet, so a junction lane whose plain Newton fails raises
  > naming the lanes, pointing at `uic=True` or gmin continuation (P18) —
  > gated alongside the bias case.  `solve_batched`'s `uic` now genuinely
  > defaults to the Parameter (False) instead of being effectively
  > mandatory.
- **P22** *(added 2026-08-21, owner request)* — re-derive the coupled (Fang)
  estimator at order 2, so the coupled path can share the shipped Gear-2
  default instead of carving out its own Euler default.  Everything around
  eq (6) was derived and measured at order 1: the acceptance band
  (0.7–3.0), the eta bracket, the sec. 3.4 step-correction law, and the
  `bordered` branch's `w'/w` denominator with its measured record.  Landing
  Gear-2 underneath the untouched machinery was measured to (a) livelock
  the coupled rectifier — NoConvergenceError at t=1.25e-4, h collapsed to
  6.3e-12, the starting evidence for this item — (b) mistune `bordered`'s
  grow-back (1181 points where 'approx' took 350 on the pulsed RC), and
  (c) break the 127/127 cross-backend step parity (96 vs 127).  Not
  uniformly worse — 'approx'+Gear-2 completed the pulsed RC in 2.8× fewer
  steps — which is exactly why the re-derivation is worth doing rather
  than the combination being refused.  Scope: re-derive the band constants
  and step law at degree 2, re-derive `bordered`'s denominator pairing,
  fix whatever the rectifier livelock's per-iteration trace pins (the
  TLine campaign's probe recipe applies verbatim), then flip the coupled
  default on BOTH backends in one commit and retire the
  `_get_integrator(coupled=True)` carve-out.  Explicit
  `integrator=Gear2Integrator()` on the coupled path is honoured today for
  anyone who wants to experiment ahead of it.

  > **EXECUTED 2026-08-21, same day.**  The re-derivation turned out to be one
  > insight, not new constants: the livelock traced to eq (6) being measured
  > on ALGEBRAIC rows — the rectifier's source-current row carries the
  > diode's dq/dt through KCL, so its accepted value holds the OLD grid's
  > derivative convention and the deviation floor (2.5e-6 A vs etol 3.6e-7)
  > is h-independent.  The **state-row mask** (zero row AND column of C ⇒
  > algebraic ⇒ excluded from the coupled LTE via a 1e30 tolerance, one
  > mechanism for the band, the controlling-node argmax, and bordered's
  > gradients) retires the class on BOTH backends; the band constants and
  > step law needed no change.  With it, coupled+Gear-2 completes the
  > rectifier in 259 points vs Euler's 769 at equal accuracy, the coupled
  > default is Gear-2 everywhere, the carve-out (and its then-dead `coupled`
  > parameter — the scan flagged it) is gone, and the old 127/127 parity
  > mystery is solved: it was BOTH backends' algebraic-row artifacts
  > synchronizing order-1 and order-2 runs by accident; masked, CPU and JAX
  > re-align through the shared default instead.  The JAX port surfaced
  > three in-trace livelock hazards, each now fixed and load-bearing: a NaN
  > band error freezing every predicate (→ inf, deterministic shrink), the
  > rejection shrinking fang's within-point-GROWN h (×5 grow vs ×4 shrink
  > never reaches the floor → shrink from the entry h), and the dt-floor
  > death march (~1 s/step on GPU × chunk_size before the raise → the chunk
  > exits on the first forced non-converged accept, which was already an
  > unconditional post-chunk raise).  A pre-mask accident is recorded
  > honestly: the algebraic-row error had been acting as a Newton-rescue
  > governor on cold starts; the mask removed the accident and the early
  > exit replaces it with the deliberate escape.  Bordered-vs-approx smooth
  > agreement re-derived (the old 5% partly rode the artifact; now 11%
  > measured, bounded 20% with both pinned to the analytic).
- **P23** *(added and executed 2026-08-21, owner request)* — the
  Spectre/Mica-style **voltage check**, `max_dv_step` on both backends: the
  largest allowed change of any node voltage in one accepted step, None
  (default) disabling it.  The scenario it exists for is the P22 mask's
  honest gap: a purely resistive/algebraic network (a designer exploring an
  amplifier topology — Rs + controlled sources, no reactances yet) has NO
  error estimator with anything to measure, so the run samples at the step
  cap and nothing reports how coarse that is.  Bounding per-step Δv controls
  output resolution directly and tracks the actual waveform; it is
  h-proportional by construction (|Δv| ≈ h·slew), so it cannot h-cancel the
  way solution-LTE did.  Node rows only; proportional retry (0.9·bound/Δv);
  vetoes compose with the LTE accept on every path (CPU standard, CPU
  coupled, JAX standard, JAX coupled).  Measured on the R+VCCS amplifier at
  1 MHz: default max per-step output change 2.256 V on a 10 V swing (61
  points, blind); with `max_dv_step=0.2` the bound holds at 0.199 on every
  path — CPU 413 points, JAX 412, JAX-coupled 777.  **Amended same day
  (owner corrections)**: (1) the bounds are FACTORS, not volts/amps —
  effective bound = max(factor, 1) · lte_vabstol (resp. lte_iabstol), so
  the check scales with the same abstols the LTE composes with, and the
  clamp-at-1 floor is exactly the solver-noise scale (a first cut
  multiplied Newton's vabstol=1e-12 with volt-scale factors and marched a
  gate at 0.2 µV resolution — 7e8 steps — before the family was
  corrected); (2) `max_di_step` is the current-row sibling — on an
  algebraic network a current waveform is otherwise bounded only through
  Δv times the conductances, with gm as an amplifier.  Worked example:
  factor 2e11 at the default lte_vabstol=1e-12 bounds steps to 0.2 V.
  **`'auto'` (owner request, the scientific setting)**: N points per period
  of a sinusoid ⇔ per-step excursion ≤ 2π·swing/N, so the auto bound is
  (2π/N)·max(declared source swing, running unit-group maximum), N =
  `points_per_period` (default 64 — resolves harmonics to ~20th order with
  Nyquist margin).  The source term comes from the new `signal_scale`
  element hook (VS/VSin/VPulse/IS/ISin/IPulse) and anchors the bound at
  signal BIRTH, where every running-reference scheme h-cancels (documented
  three times over in this campaign); the running term (P7's `ref_running`
  on JAX, accepted-only scalars on the CPU) grows the bound as an
  amplifier's output reveals gain the sources cannot know.  Measured on the
  amplifier: anchor 0.0982 V (va=1, N=64), per-step max 0.847 ≤ 2π/64·10 V
  once the 10 V output emerges, 130/129 points CPU/JAX — no hand-computed
  constant anywhere.  Curvature equidistribution (concentrating points at
  peaks and corners by bounding the interpolation error itself) is the
  recorded reconsider-if beyond this.  The gate's first e2e
  use of VCCS under the JAX toolkit also found and fixed that element's
  in-place mutation (immutable-array crash) — and the first fix attempt
  (numpy intermediate) broke the SYMBOLIC toolkit in turn (sympy gm cannot
  enter a float array; 17 ddd/ua741 tests said so).  Owner call on the root
  cause: stamp construction is the TOOLKIT's decision, not the element's —
  `Toolkit.matrix_from_entries(shape, entries)` now builds each backend's
  matrix from the element's (row, col, value) list (base: nested Python
  lists coerced through the toolkit's own array(); toolkits may override
  with native routes), and VCCS passes entries only.  The documented
  pattern for new elements.  Off by default: default runs are untouched.
- **P24** *(added 2026-08-21, owner request)* — the **stamp-construction
  hygiene sweep**: migrate every element that builds a stamp by mutating
  ``toolkit.zeros(...)`` to ``Toolkit.matrix_from_entries``.  Surveyed at
  entry: **18 zeros-then-mutate sites across 13 ``update()`` methods in
  elements.py** (plus a handful in semiconductors.py), each a latent
  VCCS-class defect — in-place assignment crashes on immutable JAX arrays,
  and any numeric intermediate rejects symbolic values; VCCS broke both
  ways in one day the first time it ran under the JAX toolkit.  One open
  question belongs to the sweep, not ahead of it: several of these
  elements DO run on the JAX path today without crashing, so their
  construction-time arrays evidently are not jnp — establish which toolkit
  each ``update()`` actually receives before rewriting, so the sweep fixes
  the real exposure instead of churning working code.  Mechanical
  otherwise; the stamp tests and the cross-backend suite are the gates,
  and every migrated element should keep its stamp bit-identical on the
  numeric toolkit (assert, don't assume).

  > **EXECUTED 2026-08-21.**  The open question resolved decisively first:
  > NONE of the 12 classes had ever been constructed under the JAX toolkit
  > — every one crashed on immutable-array assignment at the probe (the
  > fail-first record), so the sweep fixed uniform real exposure, not
  > working code.  A mechanical transformer (semicolon-splitting
  > preprocessor + strict block matcher) migrated all 16 mutate sites; the
  > two remaining `toolkit.zeros` uses are pure early-returns.  Gates held:
  > numeric stamps bit-identical 16/16 against pre-sweep captures; all 12
  > classes construct AND evaluate under JAX; the enduring cross-toolkit
  > gate compares numeric-vs-JAX stamps at near-machine tolerance
  > (VCVS_limited's JAX G comes from autodiff of eval_i_pure — same
  > derivative, different rounding) and symbolic construction keeps its
  > symbols.  The gate immediately earned its keep: **ISwitch's
  > `eval_i_pure` returned 0.0 for its control branch's defining equation**
  > (v_cp − v_cm), so on any autodiff toolkit the KVL row vanished from the
  > jacobian and the residual never enforced the equation the manual stamp
  > claims — a latent correctness defect fixed with the sweep.  Two
  > follow-on refinements, both owner-decided: NumericToolkit and
  > JAXToolkit override `matrix_from_entries` to pin float dtype (the base
  > nested-list build let all-integer stamps come out int64, which moved
  > five doctest reprs; the symbolic toolkit keeps the base — exact zeros
  > and live symbols are what it needs), and the cross-toolkit gate is
  > BIT-EXACT for the ten same-construction classes with fp tolerance only
  > for the two whose derivative METHOD legitimately differs per toolkit
  > (VCVS_limited: autodiff vs hand formula; BSource: autodiff vs
  > finite-difference `derivative`).
- **P25** *(added and executed 2026-08-21, owner request)* — **pseudo-transient
  continuation (Ψtc)**, the continuation chain's LAST rung on both backends.
  Embeds F(x) = 0 in dx/dτ = −F(x) and takes backward-Euler pseudo-time
  steps F + (1/δ)(x − x_k) = 0, J + I/δ — with g = 1/δ this is exactly a
  conductance rung with a **moving anchor**, so the adaptive exponent-space
  driver both P18 ladders share is reused verbatim (g marches 1 → 1e-12,
  δ grows 1 s → 1e12 s, halving on failure, escalating to heavier damping
  e_max = +6 if even the first step fails, finishing at exactly g = 0 —
  the P22 zero-rung rule unchanged).  CPU: `PseudoTransientNewton`, wired
  outermost in the DC chain (industry order: gmin → gshunt → source
  stepping → Ψtc); its `rung_solver` parameter keeps deformed solves off
  the failed chain, because SourceSteppingNewton's rungs rebuild F from
  the callback WITHOUT the pseudo term (the option-dropped-in-transit
  class, dodged by construction).  JAX: `ptc_g`/`ptc_anchor` on
  `dc_operating_point` plus a third `lax.cond`-gated ladder in
  `dc_with_continuation`; the traced driver's per-rung reseeding IS the
  pseudo-transient march (anchor = each rung's own seed).  Measured at
  landing: the classic Newton 2-cycle cubic x³ − 2x + 2 (plain Newton
  cycles 0 → 1 → 0 forever) lands on the real root at machine precision
  on both backends, in 3.4× fewer base calls than the zero-anchored
  gshunt ladder (47 vs 158); on tristable x³ − x every seed lands in ITS
  OWN basin (0.9 → +1, −0.9 → −1, 0.05 → 0) where a strong zero-anchored
  deformation commits every seed to 0 (measured on a forced g = 1 start;
  the shipped gshunt ladder's weak 1e-3 start does not exhibit the hazard
  on this problem — recorded, not gated).  The P18 scope finding applies
  verbatim: no legitimate CIRCUIT defeating gmin+gshunt but yielding to
  Ψtc could be fabricated, so the mechanism is unit-gated (CPU:
  `test_nrsolver_variants.py`; JAX: mock reduced system in
  `test_backend_parity_phaseA.py`) and chain-level engagement is
  exercised by the floating-node raise, which now traverses all three
  ladders.  Both reconsider-ifs were CONSIDERED same day (owner request),
  with opposite outcomes:

  > **Ψtc at a failed TRANSIENT point — EXECUTED.**  The P18 rescue chain
  > (extracted to `Transient._rescue_solver` so its topology is testable)
  > is now junction-gmin → gshunt → **Ψtc outermost**, rungs solved by the
  > plain base solver exactly as in the DC wiring.  Ψtc is the one
  > continuation that is mid-transient-safe beyond the conductance
  > ladders: its anchor is the last accepted state (physical continuity)
  > and it scales nothing — where source stepping would scale the
  > companion history, which is why source stepping stays out.  The P18
  > scope finding stands (no legitimate triggering circuit could be
  > fabricated), so ladder behavior is gated at the nrsolver level and
  > the wiring topology in `test_transient_repairs.py`.
  >
  > **SER δ-adaptation — TESTED AND REJECTED**, with the measurement
  > recorded so nobody re-chases it: SER (δ_{k+1} = δ_k·‖F_k‖/‖F_{k+1}‖,
  > one extra pure residual evaluation per rung, march clamped to
  > [0.25, 4] decades) was prototyped against the fixed decade march on
  > the P25 gate problems.  Measured base calls, fixed vs SER: cycle
  > cubic 47 vs 57 (SER worse), junction slam at budget 120: 80 vs 71,
  > junction slam at budget 20: 138 vs 71 (SER 1.9× better), tristable
  > 19 vs 29 (SER worse).  SER wins only where the per-rung iteration
  > budget is artificially tight — the regime the halving-on-failure
  > driver already covers correctly, just paying refinement rungs — and
  > regresses the easy majority while adding a SECOND marching law to a
  > driver whose value is being one shared mechanism.  Reconsider only
  > if a real circuit shows the tight-budget signature (rescue cost
  > dominated by refinement rungs after repeated maxiter failures);
  > the prototype lives in the session record, the numbers here.

  **The JAX transient-point rescue — BUILT 2026-09-01 (sec. 50).**  This
  entry recorded a deferral, "in full so the deferral is a decision, not
  an omission".  It has now been implemented, and the deferral's own
  reasoning is kept below because two thirds of it was right and the
  third that was wrong is the part worth reading.

  *What shipped.*  `transient_point_rescue` in `jaxtransient.py`: the
  CPU chain rung for rung — junction-gmin → gshunt → Ψ-tc, each ladder
  seeded from the last accepted state, each ending at EXACTLY zero, only
  a pure converged solution reported as converged (the P22 no-residue
  rule).  `newton_inner_loop` gained the same three deformation terms
  `dc_operating_point` carries, stamped identically.  It is gated at the
  dt floor — `at_floor AND not converged`, the CPU's "below minstep"
  trigger exactly — and exposed as the `continuation` Parameter,
  **default OFF**, with `statistics.gmin_rescues` named for the CPU's
  counter.

  *Why the CPU rescue could not be ported as-is* (unchanged, and still
  the reason the port looks nothing like the original).  The CPU rescue
  is Python control flow: Newton failure is an exception, caught per
  time point; below `minstep` the loop swaps in a different solver
  object (`_rescue_solver`) and runs the ladder — try/except, dynamic
  rung counts, data-dependent branching, all decided at runtime in the
  interpreter.  The JAX loop is a `lax.while_loop` traced ONCE into a
  single XLA program: no exceptions, no runtime choice of code path.

  *The one that was WRONG — and it was the design.*  The deferral said a
  port would put the ladder "behind a `lax.cond` on the failure flag".
  Built that way it would have been far worse than the entry feared, and
  silently: under `jit` a `lax.cond` executes only the taken branch, so
  it tests clean, but under `vmap` — the batched sweep this backend
  exists for — a batched predicate cannot branch, because XLA has one
  instruction stream for all lanes.  `cond` degenerates to a `select`
  that evaluates BOTH sides.  Measured: a fixed-trip rung driver behind
  such a `cond` cost **59x** against no ladder at all with **ZERO lanes
  failing**.  The fix is to gate the ladder's loop CONDITION instead
  (`_adaptive_ladder_traced(..., active=)`), which makes an unasked
  ladder run zero rungs: **0.9x, i.e. free**.  The chaining follows the
  same rule — "gshunt only if gmin failed" is an AND into the next
  ladder's `active`, not a nested `cond`.  The DC twin
  (`dc_with_continuation`) does use `lax.cond`, which is right THERE
  (small graphs, top level, once per lane) and would have been wrong
  here.

  *Cost (2) was RIGHT, and is now a number.*  The death march is real
  and irreducible: a vmapped `while_loop` runs until EVERY lane's
  condition is false, so one lane in rescue costs what all of them cost
  — measured **55.3x at 1 lane of 8, and 55.3x at 8 of 8, identical**.
  That is why the chain is gated at the dt floor rather than run
  speculatively: at the floor the alternative is force-accepting a point
  that is not a solution, so the trade is the ladder against a wrong
  answer, not against a fast one.

  *Cost (1) was RIGHT at real scale, though a small probe missed it.*
  On a synthetic 20x20 loop the graph doubled while compile time did not
  measurably move, which briefly read as the cost being unsupported.  On
  a real circuit (6-stage RC ladder with a diode, interleaved, min of 3)
  it is plainly there: 136/555/2128 steps cost 2.36x/2.12x/1.74x with
  the chain armed, and the decaying ratio is the tell — the ABSOLUTE gap
  is 1.13/1.16/1.26 s across a 16x growth in steps.  A fixed ~1.15 s of
  compile, plus ~14% per step.  Hence default OFF: on the CPU the chain
  is absent from the interpreter until it runs, here it is compiled into
  every step of every run.  `n_forced_nonconverged > 0` in the run
  report is the signal to turn it on.

  *Cost (3), the scope finding — HELD until the mechanism behind it was
  understood, then fell.*  Three circuits were tried first and none
  triggered anything (a near-singular reverse-diode chain, a floating
  capacitor node, the documented cold-start diode); the armed chain
  reproduced the disarmed run bit for bit in every case.  **The reason is
  the step-size controller**, and it is worth stating plainly because it
  is what four years of "no fabricable trigger" actually measures: a
  failed Newton REJECTS the step and shrinks dt, and with a floor at
  1e-18 there are ~130 halvings between a normal step and the floor.  The
  controller solves the point long before the rescue is reachable.  That
  is a first line of defence, not an absence of failure — and it means a
  trigger must take the controller out of the way.

  Pinning `minstep` at the step size does exactly that, and then the
  trigger is easy: a 5 V, 2 MHz source through 1 k into a diode at ten
  points per period.  The extrapolated PREDICTOR the plain Newton starts
  from overshoots the junction into its exponential, while the last
  ACCEPTED state — which is where the ladders seed — is still good.
  Measured: **the plain Newton raises `NoConvergenceError`; with the
  chain armed the run completes, 8 points rescued, 0 forced.**  The
  rescued waveform agrees with a lowered-floor JAX reference to 6.0e-2 V
  and with the CPU backend to 5.7e-2 V on a 5.633 V swing (~1%, which is
  the coarse step's truncation error, and the two independent references
  corroborate each other).  Both are tests:
  `test_the_chain_completes_a_run_the_plain_newton_cannot` and
  `test_the_rescued_waveform_is_right_not_merely_finite`.

  Note what the demonstration says about the SEED, not the deformation:
  the first thing the chain changes is where Newton starts.  That alone
  carried this case.

  **The finding that forced a second layer — now BUILT.**  The CPU chain
  sits on a solver that LIMITS (`pnjlim`); the first cut of this one sat
  on a plain Newton that does not.  For a junction circuit the dominant
  JAX failure is the diode exponential overflowing — `exp(200)` from one
  bad Newton step — and no conductance in parallel with a junction can
  fix a term that is already infinite.  Measured on the documented cold
  start: the gmin rung at g = 1e-2 diverges exactly as the plain rung
  does.

  So the chain was layered on **PCNR** as well, which is this backend's
  only limiting mechanism and therefore the analogue of the CPU's base
  solver — gmin stepping on top of a limited Newton being the
  combination SPICE has always shipped.  Both PCNR paths take it: the
  pair view (`pcnr_inner_loop`) and the device view
  (`pcnr_vector_timestep`, sec. 49), because a chain that only traced
  through the scalar path would miss every real compact model.

  The three deformations are now stamped by ONE shared helper,
  `_continuation_terms`, used by the plain Newton and both PCNR
  assemblies — a rung has to mean the same thing in all of them, or a
  ladder's final pure rung is not the system its earlier rungs were
  approaching.  For PCNR the ladder carries only `x`: convergence there
  REQUIRES `|g_lim| → 0`, i.e. `v_lim == e_a − e_b`, so a converged
  rung's limited unknowns are recoverable from its `x`, and the rescued
  point's `v_lim` (which the step controller reads) is recomputed from
  the rescued `x` for the same reason.

  *And it pays off exactly where the plain layer did not.*  A series
  diode string driven hard and fast — limiting each junction says
  nothing about the STRING — with the floor pinned: **PCNR alone raises
  `NoConvergenceError`; PCNR plus the chain completes, 12 points
  rescued, 0 forced**, agreeing with a lowered-floor PCNR reference and
  with the CPU's own PCNR to **0.13% of a 9.96 V swing** (both, to
  within 6e-6 V of each other).  Eight of sixteen swept configurations
  showed the same rescue.  Cost on PCNR has the plain path's shape:
  ~0.9–1.0 s fixed, near-constant across a 4x growth in steps.

  *Ψ-tc's rung exponents: the hypothesis above was TESTED AND DISPROVED
  2026-09-01.*  It read: "`e_start=0`, `e_max=+6` are calibrated for DC
  diagonals, and in transient at the dt floor the companion conductance is
  ~1e9, so every rung is negligible against it."  Both halves are wrong.
  The ~1e9 came from a 1e-18 floor; on the circuits where the rescue
  actually fires the diagonal is ~1e-3, so the shipped exponents are too
  LARGE.  And the sensitivity is not to scale at all: with gmin and gshunt
  disabled so Ψ-tc must carry the rescue alone, the shipped grid rescues
  **1 of 4** circuits and a **half-decade** shift rescues **4 of 4**, with
  a finer sweep failing at exactly the offsets ≡ 0 (mod 2) — the ones that
  reproduce the shipped grid of visited `g`, the step being 2 decades.  So
  it is the VISITED RUNGS, not their magnitude, and re-scaling would have
  moved the same window elsewhere.

  Left alone deliberately: gmin lands first on every rescue demonstrated to
  date, so Ψ-tc never executes in production, and tuning an exponent grid
  to four circuits with no mechanism is fitting noise.

  **What the chase did turn up is a real driver defect, now fixed.**
  `_adaptive_ladder_traced` did not search the range it declared: with no
  rung yet landed it escalated toward `e_max`, then refined just below
  `e_start`, halving until it gave up — so the exponents it could ever try
  were `{e_start … e_max}` plus `{e_start-1, e_start-0.5, e_start-0.25}`
  and nothing below, however low `e_end` was set.  Measured with a rung
  converging only for `g <= 10^k`: landed for `k >= -1`, failed for every
  `k <= -2`, against `e_end = -12`.  A descent phase now walks down to
  `e_end` before failing.  Both DC ladders use the same driver, and the
  suite is unchanged — the phase runs exactly where the old code was
  already about to give up, so it can only reach further.  It does NOT
  change the four Ψ-tc results, which is how we know the grid sensitivity
  is a separate, still-unexplained thing.

  *The trigger that reopens what remains,* stated concretely: a JAX run
  whose forced-non-converged exit survives `continuation=True` on the
  PCNR path.  That is the run whose rescue needs Ψ-tc to be doing real
  work, and therefore the one that would justify rescaling its rung
  exponents to the transient's own diagonal.

  **Two P19 defects found 2026-09-01, while scoping the remaining CPU-only
  items.** Neither was a roadmap item; both were found by reading the coupled
  path rather than by a failing test, which is the point worth keeping.

  *(a) The Fang `dh`-derivative ignored the chosen integrator.*  The residual
  is assembled with `compute_integration(..., method=eval_method, ...)`, but
  `p = d(iq)/dh` was hardcoded euler-or-gear2.  So `integrator='euler'` with
  `coupled_lte=True` integrated with Euler and then, on every non-first-order
  step, differentiated a **Gear-2** companion with respect to `h`: the residual
  and its `h`-derivative described different methods, and Fang's step-size
  Newton was solving a slightly wrong equation.  It converged on `x` anyway,
  which is why nothing ever failed.  Measured across the fix (driven RC, Euler
  both sides, reltol 1e-6):

  | | max abs(jax − cpu) | JAX steps |
  |---|---|---|
  | before | 9.1162e-03 V (1.058% of span) | 16,522 |
  | after | 3.8516e-03 V (0.447% of span) | 1,011 |

  The step count is the louder half — a wrong `dh`-derivative made the
  step-size Newton flail, taking **16x the steps to agree half as well**.  No
  test covered `coupled_lte` with a non-default integrator; one does now
  (`test_euler_coupled_uses_an_euler_dh_derivative`).

  *(a-bis) A P9 asymmetry, decided and then FIXED on the CPU.*  With
  `fixed_timestep` the JAX path bypasses the `max_dv`/`max_di` excursion check
  (`dv_ok = ... | fixed_timestep`); the CPU's coupled loop did not, and would
  shrink and reject on it — breaking the grid it was told to keep.  The JAX
  behaviour was kept when `grid_locked` landed, on P9's own principle:
  `fixed_timestep` is the caller stating that the output points are theirs.

  The CPU side was then corrected (same day).  Reproduced first, on a pulsed
  RC with `max_dv_step=1e10` (a 0.01 V bound) and `coupled_lte=True`: the
  caller asked for **21 uniform points** and got **447 spanning dt 9.0e-10 ..
  1.0e-6, with 1722 rejections** — silently, since a waveform is still
  returned.  Guarded now, with an anti-vacuity test (the bound must fire on
  the free-grid run, or the assertion proves nothing).

  ⚠ Only the COUPLED path needed it.  A matching guard was written for the
  standard path and removed again on reading the indentation: that voltage
  check already sits inside `if not fixed_timestep:`, so it is structurally
  unreachable under a fixed grid.  Measured before noticing — no bound from
  1e12 down to 1e6, at 1 V or 5 V, produced a single rejection there.  A
  guard that cannot fire is the dead-but-plausible branch this tree has paid
  for twice.

  *(c) `bordered` measured, and refused rather than ported.*  The scope note
  above says "re-derive `bordered`'s denominator pairing"; that was never
  executed, and the JAX coupled path defaults to Gear-2 — the configuration
  `bordered` is recorded as mistuned in.  Measured on the CPU first
  (`benchmarks/coupled_bordered_gear2.py`, reltol 1e-5), which is where
  experimenting is cheap:

  | circuit | integrator | bordered/approx points | newton |
  |---|---|---|---|
  | pulsed | euler | 1.02x | 1.10x |
  | pulsed | **gear2** | **2.72x** | **2.24x** |
  | smooth | euler | 0.96x | 1.05x |
  | smooth | **gear2** | 0.97x | **1.18x** |

  The 1181-vs-350 signature reproduces: the integrator mistunes it, 1.02x
  under Euler against 2.72x under Gear-2 on the same circuit.  The decisive
  part is not the ratio though — it is that **the trade is gone**.  `bordered`
  buys fewer Newton iterations for a few more time points (5.4% / 1.6% in the
  original record); under Gear-2 it costs more of both, on both circuits.  So
  it is refused on measurement rather than queued as a gap, and the JAX
  CPU-only list is now empty.  Reopening needs a `w'/w` pairing that survives
  Gear-2 on the CPU — noting that slicing the history to the integrator's
  degree is an already-recorded negative result (5x the steps), and the exact
  bordered solve is separately recorded as abandoned and not to be
  re-attempted.

  *(b) coupled + vector PCNR crashed instead of refusing.*  PCNR-inside-Fang
  understands only the PAIR-view meta `(ra, rb, IS, VT, fns)`.  Handed the
  device-view `('vector', devices, epar)` — which is what EVERY real compact
  model produces — it used the literal string `'vector'` as an array index and
  died several frames deep with "JAX does not support string indexing".  Refused at `_pcnr_setup` as a stopgap the
  same day, then **actually wired** rather than left refused: the Schur-reduced
  `(F_eff, J_eff)` is an n-sized system whose Newton step equals predict's
  `dx_mna`, so fang's machinery — eq (18)'s second solve against the same
  factors included — works on it unchanged, and that argument never depended on
  WHICH view produced the blocks.  Five sites: the assembly (through
  `pcnr_vector_blocks`, which shadows participants out of `circuit.i`/`G`, so
  the ordinary i/G must not be called on that path), the CORRECT phase
  (per-device `pcnr_limit_branchless`), the `g_lim` convergence test (PER
  COMPONENT here, `pcnr.lim_converged`'s rule — a 40 V `vds` sharing a device
  with a `vbe` would otherwise loosen the `vbe` component fortyfold), eq (18)'s
  `v_lim` correction, and the seed.  Measured against the CPU on an EKV NMOS:
  standard-path vector PCNR 2.000e-05, coupled 2.030e-05 — as accurate as the
  path already trusted, which is the claim being made.

### The stalled-estimate order drop — built both ways, measured, reverted

2026-09-01.  The CPU drops a Gear-2 step to Euler in two places aimed at a
stalled 2nd-order estimate: `Gear2Integrator.check_order_drop`'s
`h_curr/h_last < 0.1` proxy, and stage 4b's `MAX_REJECT = 3` cap.  Both were
implemented on this backend and both were removed again.
`benchmarks/order_drop_stalled_estimate.py` reproduces it.

*Ported literally, neither fires where it matters.*  The backends give up in
different places: the CPU force-accepts after three rejections at a point,
this loop shrinks to the `dt_min` floor (~130 halvings at 1e-18) and
force-accepts there.  At the floor `dt` equals the previous accepted step, so
the ratio never trips, and the point has been rejected zero times, so the
count never reaches three.  1900 force-accepts, unchanged.

*Ported faithfully — a converged floor step whose LTE failed sent back for one
retry at order 1, which is where this backend actually gives up — it works and
is still wrong.*  Against an error-controlled reference (same run, floor
lowered until nothing is force-accepted):

| reltol | circuit | force-accepts off → on | error off | error on | |
|---|---|---|---|---|---|
| 1e-6 | stiff RLC | 1900 → 516 | 5.142e-01 | 6.949e-01 | 1.35x worse |
| 1e-6 | pulsed RC | 12 → 12 | 7.156e-04 | 1.114e-03 | 1.56x worse |
| 1e-7 | stiff RLC | 1900 → 688 | 5.196e-01 | 6.946e-01 | 1.34x worse |
| 1e-7 | pulsed RC | 55 → 55 | 6.935e-04 | 4.505e-03 | 6.50x worse |
| 1e-8 | stiff RLC | 1900 → 864 | 5.208e-01 | 6.946e-01 | 1.33x worse |
| 1e-8 | pulsed RC | 600 → 480 | 7.035e-04 | 4.498e-03 | 6.39x worse |

**Worse in 6 of 6 while cutting force-accepts by up to 3x.**  The finding:
*fewer force-accepts is not a better answer.*  Dropping to order 1 makes the
ESTIMATE pass while making the actual error larger — the order-1 estimator is
agreeing with a less accurate method about a smaller claim.  The force-accept
count is a label on a step, not a measurement of its error, and optimising the
label made the waveform worse every time.

⚠ *A first version of this entry blamed the difference on the CPU's
`MAX_REJECT` cap producing force-accepts the drop then saves, and named "if
this backend gains a rejection cap" as the trigger.  The measurement above
already refutes that*: the raised-floor runs ARE that regime — 1900
force-accepts in 1908 accepted steps — and the drop lost.  Frequency is not
the variable; the drop was harmful per step it fired on.

What does explain it: `h_curr/h_last < 0.1` is a DETECTOR for "near a
discontinuity", which an Integrator object can only infer from repeated
rejection because it sees nothing else.  This backend already has the real
signal — `force_first_order` drops the order on any step landing on a
declared breakpoint (F11) — so the corners were already at order 1, and the
rule's extra firings landed on ordinary stiff floor steps where order 2 is
better (pulsed RC: ~12 declared edges, 55 firings).

Trigger: a discontinuity with NO declared breakpoint — a kink internal to a
device model, or a piecewise I-V curve whose corner is in the solution but
not in the time-breakpoint list.  There the proxy is the only detector
available.

---

## Suggested order

**Phase A (close now, all S):** P1 → P2 → P3 → P4 → P5 → P10 → P14 → P15 — one
sitting; every item lands with its test; the dead-knob scan and parameter walk guard
them afterward. P1 first: it is the safety item.

> **EXECUTED 2026-08-21**, all eight, gated by
> `test_backend_parity_phaseA.py` plus the dead-knob scan (P5's allowlist
> entry deleted — the scan itself is that item's test).  Two findings made
> during execution: threading `analysis` (P5) reached two inner loops the
> item's text did not name (`newton_inner_loop`, `pcnr_inner_loop`), and the
> item's "set the declared default" half is load-bearing — the inherited
> `Analysis` default is None, so threading without re-declaring turned every
> source's u() to zeros (caught by the rc-charging gate as an all-zero
> waveform, 26 failures); and the
> P2 behavioral assertion is CPU-only by measurement — on every probed
> configuration the JAX gear-2 charge LTE sits so far under tolerance that
> the run is step-cap-limited and NO tolerance knob moves the step count
> (reltol 1e-4 vs 1e-6: identical 209-step runs), an order/headroom artifact
> worth knowing when comparing shipped defaults (P6's territory), so the JAX
> side pins the kernel contract instead (`ywr_error_ratio` scales exactly as
> 1/TRTOL).  `solve_batched`'s floor default moved 1e-15 → 1e-18 (the shared
> `minstep` Parameter), per this item's own text.

**Phase B (schedule):** P6 (+ the default-order decision, which is the owner's) →
P12 (delete dead code, share real machinery) → P8 → P9 → P7 → P11 → P16. P6 before
P8/P9 because both compose with the integrator choice; P12 early because it removes
a false feature.

**Phase C (document):** P13, P17, P19, P20 paragraphs in the two class docstrings —
half an hour, and the parity story becomes checkable: everything is then either
aligned, scheduled, or stated.

The conformance harness is the standing instrument here too: each Phase-B port adds
its configuration to the harness matrix, so "feature-comparable" stays a tested
property rather than a snapshot claim.
