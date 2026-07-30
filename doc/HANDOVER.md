# Handover — where the leapfrog/symbolic work stands, 2026-07-30

**Read this first if you are picking the symbolic work back up.** The transient work
(`doc/transient_work_plan.md`) is expected to run first and to span several sessions, so
this document exists to make the leapfrog thread resumable after it, not to be read during
it.

Everything below is committed. A substantial backlog is **unpushed** on
`cna-jax-vectorization` — `git log --oneline origin/cna-jax-vectorization..HEAD | wc -l`
for the current count, deliberately not pasted here because it goes stale on the next
commit. **Pushing is the maintainer's call and has not been done.**

---

## 1. The one-paragraph summary

The goal is full symbolic extraction of a 5th-order leapfrog filter built from five µA741
amplifiers, then approximation to expressions a designer can read. Numerics may *rank*
terms; they must not replace them. Chasing a transient cross-check for that work uncovered
that the benchmark fixture was an **unstable circuit**, then that it was an **uncompensated
one**, and separately that the transient engine had four defects. The fixture is now
repaired and compensated; the transient engine is partly repaired and fully reviewed. What
remains on the symbolic side is finishing the third regeneration of the tables and then the
transient-vs-perturbation comparison that started all of it.

---

## 2. State of the fixture

`leapfrog_5th_order` has changed **twice** today. Any number in any document must be read
against which version produced it.

| version | dim | nnz | Q | tau | passband peak | commit |
|---|---|---|---|---|---|---|
| original (UNSTABLE — 2 RHP poles) | 127 | 536 | — | — | — | before `ff5c6e6` |
| repaired topology | 127 | 536 | 16.76 | 208 us | +8.79 dB | `ff5c6e6` |
| **+ GBW compensation (current)** | **136** | **567** | **5.93** | **73.5 us** | **+0.000 dB** | `95545e5` |

- The instability was a sign error: backward coupling entered the same summing node as the
  forward coupling, so each stage integrated the *sum* of its neighbours where a ladder
  integrates their *difference*. Confirmed on four independent routes.
- The compensation is a 350 ohm resistor in series with each integrating capacitor
  (`_LEAPFROG_RC`), correcting finite-GBW Q enhancement. **`rc=0` recovers the
  uncompensated circuit exactly** — that is how the +9 nodes are attributable.
- `build_leapfrog_network` now builds the **entire** circuit including the source,
  dispatching on toolkit (symbolic -> `VS`; numeric -> `VSin` per tone). Nothing is
  replicated anywhere. That was itself a fix: the IM3 harness had hand-carried a copy of
  the topology and kept the *broken* wiring after the fixture was repaired.

---

## 3. What is finished

- **Transient tolerance fix** (`e37ddad`): the step controller was handed the
  residual-flavoured tolerance vector where it needed the solution-flavoured one, and
  `vabstol` defaulted to 1e-12 V against Spectre's 1 uV. 19x fewer steps. **See the open
  issue in section 5 — this fix has a known side effect.**
- **Gear2 LTE repair** (`doc/transient_repair_plan.md`, stages 1-5, all gates recorded):
  the `'classic'` estimate computed `q''*h^3` where BDF-2 needs `q'''*h^2`. Independently
  re-verified by the later review at est/true = 1.000282 against the 2/9 constant.
- **T2/T3 first redo** (`2d109c5`): every leapfrog experiment re-run against the repaired
  127-unknown fixture, all tables regenerated.
- **Four-lens transient review** (`7d10116`): `doc/transient_review.md` plus evidence in
  `benchmarks/transient_review/`.
- **The transient work plan** (`doc/transient_work_plan.md`) — the next thread.

---

## 4. What is IN FLIGHT and must be finished

### 4.1 T2b — the third regeneration, COMPLETE

The fixture changed again (compensation), so every leapfrog number needed regenerating a
**third** time. **All seven scripts have now run; logs are in `doc/t2b_logs/`.**

| script | result on the 136-unknown fixture |
|---|---|
| `cancellation_leapfrog` | kappa **1.153e+12** (unchanged), 7 913 groups, 14 409 600 terms, log10\|det\| −397.2 |
| `cancellation_blocks` | fails — **pre-existing**, verified identical against the pre-repair fixture |
| `cancellation_parallel` | top kappa 1.194e+12; **blocks bit-identical** (D_k 1.147e+03) |
| `cancellation_compose` | GATE 5 6/6 PASS, GATE 4 2/6 PASS — verdict unchanged |
| `nonlinear_leapfrog` | GATE 14-2 PASS rel **2.60e-13**, GATE 14-3 PASS; speedup **27-29x** (was 61-68x) |
| `transfer_function` | leapfrog tol=1e-3: 7 913 groups -> **134 520 ops = 17.0 ops/group**; uA741 control **bit-identical** |
| `order_convergence` | **`agrees from` IDENTICAL for the third time**; `v_turn` = 5.3621e-02 V |

The `transfer_function` control is worth keeping in view: its uA741 half is untouched by
every fixture change and has now re-measured **734 groups -> 50 377 operations across all
three versions**. Numbers that should not move have not moved, which is what makes the
leapfrog deltas trustworthy rather than merely different.

One recorded claim shifted: **"16 ops/group is a property of the diagram" now reads 17.0**
(it was 16 on the unstable fixture and 16.0 on the repaired one). Still structural, but no
longer identical — record it as a rescale, not a survival.

**THE RESULT WORTH HAVING.** `order_convergence` is the script whose output feeds the
pasted §10.2/10.3 tables, and it is the sharpest test of the campaign's one durable
finding. The `agrees from` column is **identical for the third time**, across three
structurally different circuits:

| amp (V) | unstable 127 | repaired 127 | **compensated 136** |
|---|---|---|---|
| 0.01 | U^5 | U^5 | **U^5** |
| 0.03 | U^5 | U^5 | **U^5** |
| 0.1 | U^5 | U^5 | **U^5** |
| 0.3 | U^7 | U^7 | **U^7** |
| 1 | U^9 | U^9 | **U^9** |
| 3 | not by U^13 | not by U^17 | **not by U^17** |

The scale moved, as it must: `v_turn` 5.3656e-02 -> **5.3621e-02 V** (`g` at `s0_e1`
4.3184e-04 -> 4.3128e-04 S), node voltages down ~19%, and `% of v_turn` now 0/0/1/3/11/32
against 0/0/1/4/13/39. The `kk` sweep holds its shape too: U^7 at 3%, none by U^17 at 34%,
and `kk = 50` still correctly flagged **UNPHYSICAL at 340% of v_turn**.

**So the claim has now survived two independent perturbations of the circuit underneath
it** — which is the strongest form the section-7 finding has taken: it measures the series
against the cubic's own validity limit, and is invariant to what the matrix does.

Run any of these as:
`cd benchmarks && PYTHONPATH=<repo>:<repo>/benchmarks MPLBACKEND=Agg python3 -u <script>.py > log 2>&1`

### 4.2 T3 — the doc rewrite for the 136-unknown fixture, NOT STARTED

The tables in `doc/distortion_ddd_conclusions.md` §10.2/10.3 and
`doc/cancellation_ranking_conclusions.md` currently hold **127-unknown** numbers, correctly
labelled as such. They need the 136-unknown regeneration once 4.1 completes.

**Two traps, both hit before and both recorded in `doc/leapfrog_redo_plan.md`:**
- **A normal `sphinx -b html` does NOT re-run live `exec-rst` blocks** when only the *code*
  changed. Sphinx reuses cached doctrees for unchanged sources, so the built page keeps
  stale numbers while still reporting `build succeeded, 2 warnings`. **Force it:**
  `rm -rf doc/build/doctrees` and `sphinx -E`.
- **Verify in the block's own output format.** The live block prints 3 decimals where
  `order_convergence.py` prints 4, so grepping for `7.1404e-05` returns zero hits — and so
  does the *old* value, which reads exactly like a third kind of staleness. Read the
  rendered table.

### 4.3 T4 — the IM3 transient comparison, BLOCKED

`benchmarks/nonlinear_leapfrog_sweep.py`. Two-tone IM3 (100/110 kHz, product at 90 kHz),
because **HD3 is unmeasurable here by any transient at any cost**: the third harmonic at
300 kHz is attenuated 160 000x by the filter's own stopband while the fundamental loses
only 106x. IM3 lands beside the fundamentals and is 277x larger.

**Blocked on `doc/transient_work_plan.md` stage 4g.** The harness sets
`TrapezoidalIntegrator()` and drives two `VSin` tones — which is exactly the mechanism that
seeds the trapezoidal LTE estimator's parasitic `(-1)^n` mode, making its error O(1/h). The
integrator was chosen on a step-count comparison a contaminated estimator could have
produced. **Gate 0.2a of the transient plan re-verifies this; do not quote a T4 number
before it passes.**

Cost, as last measured (on the *uncompensated* fixture, so re-measure): ~6.8 h per
amplitude, dominated by settling rather than stepping. The compensation cut tau 208 ->
73.5 us, which should reduce that ~2.8x. The larger lever is **seeding `x0` with the linear
two-tone steady state** — the circuit is linear apart from a cubic contributing ~1e-4, so
one AC solve removes nearly all the settling. Not implemented.

---

## 5. Open decisions — these need the maintainer, not an implementer

1. **`vabstol` serves two roles** (`doc/transient_work_plan.md` 0.3a). It is Newton's
   x-tolerance *and* the LTE tolerance. The 1e-12 -> 1e-6 fix in `e37ddad` was reasoned
   about only as a step-control knob, so it **also loosened Newton's node convergence by
   10^6, unmeasured**; and `DC.vabstol` is still 1e-12, so the operating point seeding
   every transient is solved 10^6 times tighter than any step after it. All four review
   lenses flagged it.
2. **Gate 1-3 of `transient_repair_plan.md`** is still open: its 2x threshold was written
   against a metric later shown to be degenerate, and the sharp reading gives 3.94. Needs
   a threshold justified by the 2/9-vs-1/12 error constants, or an explicit "recorded, not
   resolved".
3. **`Gear2`'s default `lte_formula`** (0.3b) — the evidence now runs against the `'ywr'`
   default that was chosen belt-and-braces.
4. **Scope of transient stage 10** (0.3c) — the missing-analyses list is a product
   decision.

---

## 6. Refuted premises — do not re-derive these

Recorded because each cost real time and each was *stated before it was measured*.

- **"The stable circuit will be cheaper to integrate."** False: 2896 steps against 571,
  5x *more*. A smooth exponential blow-up is easy to track; a lightly-damped high-Q
  resonance is not.
- **"kappa is driven by the Q resonance."** False: kappa is **1.153e+12 on both the
  Q=16.76 and Q=5.93 fixtures**, identical to four significant figures. The numerical
  review independently explains why — `log10|det|` is a *unit-scale* artefact, moving +365
  decades as h goes 1e-9 -> 1e-12, which is 3 decades x ~122 ~ `rank(C)`. kappa here is a
  property of the determinant representation, not of the circuit.
- **"Compensation will make the benchmark easier."** False: 6.4x more terms, 32% more
  groups, and the determinant underflows *further* (−397 vs −374).
- **"The transient never ran because the circuit is stiff."** False. The mean step was
  5.03 ns against a 39 ns cap — a ratio that says steps are being *rejected*, not that the
  circuit is stiff.
- **`cancellation_blocks`'s failure is caused by the fixture change.** False — it fails
  identically against the pre-repair fixture. Verified by checking out `ff5c6e6^`'s
  `benchmark_circuits.py` and re-running.

---

## 7. The one durable finding from the redo campaign

**What was refuted measured the determinant's *conditioning*; what survived measured the
*method*.** A topology change alters matrix conditioning and leaves the properties of the
series, and of the diagram representation, alone.

Survived both fixture changes unchanged: the order-convergence pattern (`agrees from` was
identical at all six amplitudes through the first redo), 16.0 ops/group, the symbolic
nonlinear analysis matching its numeric oracle, "readable H(s) over 2 decades only".
Refuted or rescaled: kappa, group counts, term counts.

Read any future claim in these documents against which of those two it depends on.

---

## 8. Environment notes that will save an hour

- `PYTHONPATH=/home/andreas/sources/pycircuit` and `MPLBACKEND=Agg` are **required** — a
  stale root-owned egg in `/usr/local/lib/python3.14/dist-packages` shadows the source.
- **Two venvs, not interchangeable.** The repo's `.venv` has sphinx but no test deps; a
  scratch venv (`python3 -m venv --system-site-packages`, then `pip install pytest pynose
  pytest-timeout`) runs the suite but has no sphinx.
- Suite: `pytest pycircuit -q -p no:cacheprovider -m "" --timeout=400 --timeout-method=signal`.
  **`-m ""` is mandatory** — 17 `slow` tests are deselected by default and they are the only
  ones comparing an analysis against a time-domain reference.
- Doc build: `cd doc && MPLBACKEND=Agg ../.venv/bin/python -m sphinx -b html -d build/doctrees src build/html`.
  Clean baseline is `build succeeded, 2 warnings`.
- **Never pipe a long-running command through `tail`/`grep`/`head`** — output is lost
  entirely if it is killed, and a 0-byte log is indistinguishable from "produced nothing".
  Redirect `> log 2>&1` and read the file. Use `python3 -u`.
- To check whether a background job is alive, match on `comm` (the executable), never on
  the command line — a `bash -c` wrapper contains the full python invocation and matches
  too.

---

## 9. Document map

| document | what it is |
|---|---|
| `doc/transient_work_plan.md` | **the next thread** — staged plan, gates declared, nothing run |
| `doc/transient_review.md` | the four-lens review and its measurements |
| `benchmarks/transient_review/` | evidence probes for the review's performance/solver numbers |
| `doc/transient_repair_plan.md` | the completed Gear2 LTE repair, gates recorded |
| `doc/transient_repair_reasoning.md` | why that repair was scoped as it was |
| `doc/leapfrog_redo_plan.md` | the fixture-change redo campaign; T0/T1/T2/T3 outcomes |
| `doc/distortion_ddd_conclusions.md` | §10.1-10.3, the nonlinear leapfrog results |
| `doc/cancellation_ranking_conclusions.md` | the cancellation/group-ranking results |
| `doc/t2b_logs/` | the five completed T2b runs on the 136-unknown fixture |
