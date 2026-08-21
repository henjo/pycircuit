# Phase-0 baseline — 2026-08-21

Recorded per `doc/transient_review_260820.md`, Phase 0 item 6, after the Phase-0 fixes
(F2(a), F6(a), F8, F9(a), F18) landed. These are the yardsticks the Phase-2/3 fixes are
judged against: each behavior-changing fix records its own delta as it lands; suite
expectations are pinned once, at the Phase-3 exit gate.

Machine: ThinkPad P16 Gen 1, i7-12800HX, RTX A1000 (4 GB), venv Python 3.14.4,
jax 0.11.1/CUDA 13, `XLA_PYTHON_CLIENT_PREALLOCATE=false`.

## Performance (benchmarks/phase0_baseline.py — one warm-up, one timed run)

rc-vsin (VSin 1 kHz → R 1k → C 100n), `tend=2e-3`, `timestep=1e-5`:

| configuration | wall | accepted | rejected |
|---|---|---|---|
| CPU `Transient`, Euler (shipped default) | 150.1 ms | 412 | 8 |
| CPU `Transient`, Gear2 | 78.7 ms | 211 | 0 |
| JAX `JAXTransient`, gear (shipped default) | 612.5 ms | 210 | 0 |

Three facts worth naming, all relevant to the review:

- **At matched order the backends agree almost exactly** (211 vs 210 accepted steps) —
  the step-count parity F19(b) says is only visible once the integrator is pinned. The
  shipped defaults compare 412 vs 210, i.e. different methods, which is F19(b)'s point.
- **The JAX path is ~7.8× slower than CPU-Gear2 on this single small circuit.** Not a
  defect — n≈3 matrices cannot feed a GPU, and the JAX path's value proposition is
  `solve_batched` — but until now that trade-off was an impression, not a number. Any
  Phase-3 fix that moves this ratio materially in either direction should say so.
- **CPU-Euler's 8 rejections vs Gear2's 0** on the same circuit is the F17-adjacent
  population to watch when the JAX safety factor lands.

## Coverage (coverage 7.13, `--source=pycircuit/circuit`, the eleven transient-related
test files at `-m ""`; 204 tests, 72 s)

| module | stmts | miss | cover |
|---|---|---|---|
| `_lte_kernels.py` | 107 | 9 | 92% |
| `integrator.py` | 97 | 10 | 90% |
| `transient.py` | 735 | 110 | 85% |
| `jaxtransient.py` | 477 | 73 | 85% |
| `stepcontroller.py` | 174 | 49 | 72% |
| `nrsolver.py` | 241 | 72 | 70% |

Dark territory worth knowing before trusting a green suite (from `--show-missing`):

- **`transient.py:805–857` — `_apply_voltage_ics`, the spanning-tree capacitor-IC solve,
  is essentially untested**, along with most of the element-IC machinery (686–704,
  728–752). The review praised its design from reading; nothing currently guards it.
- `transient.py:536–549` — the `_newton` exception-classification branches
  (SingularMatrix promotion) are unexercised.
- `nrsolver.py` at 70% — the DampedNewton backtracking paths (the F15 territory) and
  `JAXNewtonSolver` are largely dark, consistent with F15's "delete-or-port; grep for
  users first".
- `stepcontroller.py` at 72% — the PI controller and parts of the band logic (the F10
  territory).

  > **ADDRESSED 2026-08-21** (owner request).  Re-measured under the full
  > post-campaign suite first: 78%, with the genuinely dark regions being
  > both controllers' LTE-solve fallback (decision 0.3d's kept-but-loud
  > units-violation warning — measured firing zero times in production,
  > guarded by nothing) and the WHOLE of
  > `SolutionLTEController.evaluate_step`, a selectable standard-path
  > controller whose design record cites measurements no test pinned.
  > `test_stepcontroller_units.py` now gates: the fallback warning + its
  > stand-in verdict on both controllers, IntegralController's
  > lower-band growth-retry with all three suppression guards
  > (h_clamped, at-max_step, damper-would-not-grow), PIController's
  > gamma_min refusal and sentinel pass-through, and the solution
  > controller's full standard-path contract (accept-unevaluated
  > openings, controlling-index record, both band edges).  **After:
  > 99%** — the one uncovered line is the ABC's `pass`, left honestly.

Reproduce: `XLA_PYTHON_CLIENT_PREALLOCATE=false .venv/bin/python -m coverage run
--source=pycircuit/circuit -m pytest -q -m "" <the eleven files>` then
`coverage report --show-missing`.

### Coverage verification pass — 2026-08-21, post-campaign

Re-measured under the FULL suite (1010 passed), so the numbers below are
not directly comparable to the Phase-0 table's eleven-file basis — they
are the campaign-exit state of the same modules:

| module | Phase-0 | verified now |
|---|---|---|
| `stepcontroller.py` | 72% | **99%** (the ABC's `pass` only) |
| `nrsolver.py` | 70% | **90%** |
| `transient.py` | 85% | **96%** |
| `jaxtransient.py` | 85% | **97%** |
| `integrator.py` | 90% | 93% |
| `_lte_kernels.py` | 92% | 96% |
| six-module total | — | **96%** |

Every named dark spot is confirmed addressed: the spanning-tree
capacitor-IC solve and element-IC machinery (P12's binding + the IC
suites), the `_newton` exception-classification branches (gates 6-1/6-2
— one leg remains, see below), the DampedNewton backtracking paths
(F15's fix and gate), and `JAXNewtonSolver` (DELETED per F15; a
tombstone comment marks the site).  `stepcontroller.py`'s territory was
retired separately (see the ADDRESSED note above).

Residual dark territory, characterized so silence isn't read as
clearance — all of it rare-guard code, none of it load-bearing logic:

- `nrsolver.py`: parts of `_structural_singularity`'s diagnosis branches
  (symbolic fall-through, dead-row messages), DampedNewton's
  singular-classify twin of the StandardNewton branch that IS covered,
  and `SchurCoupledNewton`'s linearsolver-failure / tiny-denominator
  fallbacks.
- `transient.py`: scattered single-line guards; the one named leg is
  `_newton`'s raw-`LinAlgError` promotion (713) — rare because the
  solvers wrap factorisation failures as `NoConvergenceError` first.
- `jaxtransient.py`: the coupled path's FD du/dt fallback for circuits
  whose `dudt` raises `NotImplementedError` (1115–1134; the analytic
  branch is what production circuits hit), and `solve_batched`'s TLine
  setup block (2372–2377) — batched TLine has no test, consistent with
  TLine-under-batch never being claimed.

Chasing these to 100% would mean fabricating failures inside guards that
exist precisely because they are hard to reach; the honest disposition
is this record.

## Phase-3 exit gate (same day)

Recorded after the Phase-3 cluster (F19, F11, F17, F6(b)+F16(a), breakpoint cap) landed.

Performance (same benchmark, same machine):

| configuration | Phase-0 | Phase-3 exit | delta |
|---|---|---|---|
| CPU Euler (default) | 150.1 ms / 412 steps | 150.9 ms / 412 | unchanged (untouched) |
| CPU Gear2 | 78.7 ms / 211 | 78.2 ms / 211 | unchanged (untouched) |
| JAX gear (default) | 612.5 ms / 210 | 544.7 ms / 210 | **−11% wall** (F17's rejection cut) |

Conformance harness (matched order, reltol 1e-4):

| pair | Phase-1 landing | Phase-3 exit | bound now |
|---|---|---|---|
| rc | 7.4e-4 | **4.3e-6** (173× better) | 2e-5 |
| vsin | 3.5e-3 | 2.9e-3 | 1.5e-2 |

Step counts at matched order: rc 61/60, vsin **120/120** (was 120/113).
Rejection ratios (rc-vsin): 44.5%→20.1% at reltol 1e-4, 40.1%→**4.0%** at 1e-6.
Pulse-edge rejections: 38→16 of ~190 accepted (F11).
Suite: 1001 passed, 12 skipped, 3 xfailed under `-n 8`.

Per-fix deltas were recorded in each fix's commit message; suite expectations are
pinned in the harness bounds above, per the gate's acceptance rule.

## Execution record — all phases complete (2026-08-21)

All 27 schedule items of `doc/transient_review_260820.md` executed, F1–F19 + R1/R2
hardening + hygiene + F2(b) feature work. 24 commits (`f722372..56a143b`). Suite grew
970 → **1005 passed** (12 skipped, 3 xfailed) under `-n 8`; the dead-knob scan and the
cross-backend conformance harness (bounds tightened at the Phase-3 gate) now guard the
two defect classes the review identified.

Findings made during execution, recorded per the review's own conventions:

- The dead-knob scan found a **fourth** dead knob on landing (`ywr_error_ratio`'s
  `x_last`), beyond the meta-review's three.
- F19's chunk-invariance issue was caught by an existing gate mid-landing: the estimator
  predicate `step_idx < 2` was chunk-local; the run-global history facts replaced it,
  fixing a pre-existing chunk-boundary artifact in the *integration* fallback too.
- F17's before/after: 44.5%/40.1% → 20.1%/**4.0%** rejection ratio at reltol 1e-4/1e-6.
- One hygiene row was **tested and rejected**: the bordered-branch degree slice collapsed
  that method to 5× the steps — the full-history gradient pairs with the full-history
  `w'/w` denominator, not with `err`'s degree. Reverted with the measurement recorded
  at the site.
- Three vacuously-green assertions were re-derived when their counters went live
  (`test_coupled_method`'s pulse zeros, gate 9-1(b)'s reject case, gate 6-3's
  point-count identity).
