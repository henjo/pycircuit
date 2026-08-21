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

| Feature | CPU `Transient` | `JAXTransient` | Verdict |
|---|---|---|---|
| `reltol`, `iabstol`, `vabstol`, `lte_vabstol`, `lte_iabstol`, `maxiter`, `max_step`, `firststep`, `epar` | ✓ | ✓ | aligned (post-review) |
| t=0 point, tend landing, breakpoint order drop, force-accept accounting, per-row Newton, safety-factor law, statistics-on-result | ✓ | ✓ | aligned (Phases 1–3) |
| **P1** `uic` / `minstep` / typo safety | Parameters | `**kwargs`, silently swallowed | **align now (S)** |
| **P2** LTE-ratio knob | `LTERATIO = 7.0` class constant | `TRTOL` Parameter | **align now (S)** — reversed asymmetry |
| **P3** `solve_batched` floor & naming | — | `dt_min=1e-15` vs solve's `1e-18`; arg `irefnode` vs `refnode` | **align now (S)** |
| **P4** `refnode` default | `gnd` | `0` | **align now (S)** |
| **P5** `analysis` parameter | `'tran'`, threaded | inherited, dead; loop hardcodes `'tran'` | **align now (S)** |
| **P6** integrator selection | `integrator=` Parameter (Euler/Trap/Gear2) | internal `'euler'/'gear'`, not exposed; hardcoded `'gear'` | **align (M)** + a default-order decision |
| **P7** `relref` modes | `pointlocal`/`alllocal`/`sigglobal` | `sigglobal` only, hardwired | align later (M) |
| **P8** LTE band (`lte_gamma_*`, `lte_eta`) | ✓ (standard + coupled) | absent; accept test hardcodes `err <= 1` | align later (M) |
| **P9** `fixed_timestep` | ✓ | absent | align later (M) |
| **P10** `outputstep` uniform-grid resample | ✓ | absent — but the resampler is backend-agnostic | **align now (S)** |
| **P11** `provided_function` extra source | ✓ (`pf(t)`) | absent | align later (M, with a traceability constraint) |
| **P12** initial conditions under `uic` | `ic` dict + element ICs + spanning-tree solve | reads `node.ic` attributes — **dead code, nothing sets them** (new finding) | **align (M)** — delete the dead path, reuse the CPU machinery |
| **P13** `bypass`/`bypasstol` | ✓ | absent | document — meaningless under vmap |
| **P14** `minbreak` | ✓ | progress guard only, no spacing merge | align now (S) |
| **P15** statistics vocabulary | `force_accepts`, `order_drops`, `breakpoints_hit`, timing | `forced_lte_steps`, `nonconverged_steps`, `signal_max` | partial align (S) |
| **P16** TLine delay interpolation | quadratic (3-point Lagrange) | linear + ring-cap raise | align later (M) |
| **P17** solver strategy objects (`nrsolver`/`scaler`/`linearsolver`) | ✓, injectable | refused loudly | **document — permanent** (traced loop cannot dispatch) |
| **P18** continuation (gmin/source stepping) | DC path; transient retries dt | dt retries only | roadmap, both sides |
| **P19** `coupled_lte` (Fang), `pcnr`, `coupled_method` | ✓ — research paths with measured records | absent | **amended: portable — Fang then PCNR scheduled** (see P19 note) |
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

This also forces the deferred F19(b) decision into the open: **the shipped defaults
differ** (CPU Euler, JAX Gear-2), so identical scripts get different methods by
backend. Recommendation: make **Gear-2 the CPU default too** — the Phase-0 baseline
measured it at half the steps and half the wall-clock of Euler on the same circuit at
the same tolerance, the estimator work of stages 4g/4i was built for it, and the
conformance harness already treats it as the reference configuration. This is a
user-facing waveform change on the CPU path (every default run's step sequence
moves), so it is an owner decision; the fallback is documenting the asymmetry in both
`integrator` Parameter descriptions. Either way the harness pins whichever is chosen.

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

### P8 — the LTE acceptance band on JAX

Post-F17 both backends share the prediction law, so the band port is the acceptance
test only: `gamma_min <= err <= gamma_max` with the growth-retry (reject with a
larger dt — the machinery exists: `do_reject` with `retry_dt` can grow) and the eta
damper as a clamp on `factor`. Reuse the CPU's `'auto'` sentinel semantics from F5.
Moderate; do after P6 so the band composes with a chosen integrator.

### P9 — `fixed_timestep` on JAX

Semantics per the CPU: grid wins, breakpoints don't move it, order still drops when
a step crosses an edge, non-uniform fallback warns. In the traced loop this is
`dt_next = timestep` (skip `calculate_next_dt`) plus the accept test forced true
except on Newton failure. Moderate, mostly tests.

### P11 — `provided_function` on JAX, with the honest constraint

The CPU contract (F4) is an extra source term `pf(t)`. In the traced loop, `pf` must
be jax-traceable (pure, jnp-composable) and is baked in at jit time. That constraint
is acceptable — state it in the Parameter doc and raise a clear error if tracing
fails. The DC-seed caveat from F4 applies verbatim (JAX's non-uic seed is a plain
DC solve that won't see `pf`): same warning, same wording.

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

### P16 — TLine interpolation order

CPU delay lookup is 3-point quadratic with recorded reasoning ("a first-order
interpolant feeding a second-order method injects error the method never made");
JAX's is linear. Port the quadratic stencil into `interp_tlines` (three gathered
points instead of two; same modulo indexing). Measure the edge-accuracy delta on the
existing TLine tests before/after, per house rules.

---

## Documented divergences (no code change; one paragraph each where users look)

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
  stay CPU-only until separately justified.
- **P20 `solve_batched`:** JAX-only by design — it is the branch's purpose. The CPU
  gets no imitation API; a Python loop over `Transient` is already expressible and
  honest about its cost.

## Roadmap items (neither backend has it; noting so silence isn't read as parity)

- **P18** transient-level continuation on a hard step (gmin ramp at a failed time
  point) — the review's "deliberately unscheduled" gmin item, unchanged.
- **P21** batched DC operating point, so `solve_batched` can start from bias. The
  loud `uic=True` refusal stays the contract until then.

---

## Suggested order

**Phase A (close now, all S):** P1 → P2 → P3 → P4 → P5 → P10 → P14 → P15 — one
sitting; every item lands with its test; the dead-knob scan and parameter walk guard
them afterward. P1 first: it is the safety item.

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
