# Fixing the leapfrog topology and redoing every experiment

**Status: written 2026-07-30. BLOCKED on the topology fix. No stage has run. Every
OUTCOME line is deliberately blank and must be filled from a real run.**

Maintainer's instruction, 2026-07-30: *"We have to fix the topology and redo every
experiment."* This plan is the enumeration of "every experiment", because a redo campaign
that works from memory will miss something, and a missed artifact is worse than no redo —
it leaves a stale number that still reads as current.

## Why

`leapfrog_5th_order` appears to have **two right-half-plane poles**
(s = +1.4491e+05, +5.6716e+04; growth time constants 6.9 us and 17.6 us), found as
generalized eigenvalues of the pencil `det(A0 + s*C) = 0`. Corroborated circumstantially:
a two-tone transient reached 6.3e+03 V at the output where 6.7e-05 V was expected, and
`e^(1.4491e5 * 120e-6) = 3.6e7` turns a ~1.8e-4 V startup kick into ~6.4e3 V.

**What this does and does not invalidate.** An unstable circuit still has a well-defined
`H(s)` on the jw axis — the analysis just solves a linear system at `s = jw` — so every
recorded frequency-domain and symbolic number is *arithmetically correct*. What is wrong
is the physical claim. The results describe a circuit that would not operate: it has no
steady state, its "filter response" is the response of a divergent system, and any
statement of the form "a designer could read this expression and understand the filter"
is describing a filter that does not exist. That is why annotating is not sufficient and
the numbers have to be regenerated.

**Anything that only ever needed a matrix is cheap to redo.** Most of these experiments
are symbolic or frequency-domain, so the redo cost is dominated by re-running scripts, not
by new development. `order_convergence.py` in particular keeps `Adrv` and `kk` symbolic, so
one build per order still serves the whole amplitude sweep.

## Blast radius, enumerated

Grepped, not remembered (`grep -rln leapfrog_5th_order --include=*.py`).

### Code that consumes the fixture

| file | what it measures | redo cost |
|---|---|---|
| `pycircuit/circuit/tests/test_benchmark_circuits.py` | fixture invariants | seconds |
| `benchmarks/order_convergence.py` | U^3..U^17 convergence, amplitude + `kk` sweeps, `v_turn` | ~minutes |
| `benchmarks/nonlinear_leapfrog.py` | stage 14: symbolic build vs numeric oracle | ~minutes |
| `benchmarks/cancellation_leapfrog.py` | the 181-group hierarchical result | ~minutes |
| `benchmarks/cancellation_blocks.py` | per-block cancellation | ~minutes |
| `benchmarks/cancellation_compose.py` | composed error | ~minutes |
| `benchmarks/cancellation_parallel.py` | parallel elimination | ~minutes |
| `benchmarks/transfer_function.py` | part 1 re-reads the 181-group verdict in operations | ~minutes |
| `benchmarks/nonlinear_leapfrog_sweep.py` | two-tone IM3 vs transient | **hours** |

### Documents carrying pasted numbers (these do NOT self-update)

- `doc/distortion_ddd_conclusions.md` — sections 10.1-10.3, including the corrected
  amplitude/`v_turn` table and the `kk` sweep table
- `doc/cancellation_ranking_conclusions.md` — leapfrog rows across several sections
- `doc/architecture.md` — P16
- `doc/distortion_ddd_plan.md`, `doc/cancellation_ranking_plan.md`,
  `doc/hierarchical_approximation_plan.md` — stage outcomes referencing the leapfrog
- `doc/ddd_references.md` — if it quotes leapfrog figures

### Documents that DO self-update

- `doc/src/circuit/distortion_ddd.rst:500-525` — a live `exec-rst` block calling
  `bc.leapfrog_5th_order()`. It regenerates on every build, which is the
  "never type a measured number into prose" rule paying for itself. **But** a live block
  that raises renders its own source and the build still says "succeeded", so it must be
  verified by grepping the built HTML for a *computed* value, not for a heading.
- `doc/src/circuit/ddd.rst:1020` — prose referencing the 181-group figure; check whether
  the number is pasted or computed.

## Stages

### Stage T0 — land and verify the topology fix

Blocked on the separate instability investigation. Requirements before anything downstream
runs:

**Gate T0-1 (stability, two independent routes).** All finite poles have Re < 0, confirmed
by two methods that do not share an implementation (e.g. `scipy.linalg.qz` on the pencil,
and a driven transient whose peak stops growing with run length). Report condition numbers
— the original finding rests on a generalized eigenproblem with a singular `C`, and a
2-of-125 result wants corroboration.
OUTCOME:

**Gate T0-2 (it is still the intended filter).** `|H|` across 100 Hz - 1 MHz must look like
a 5th-order lowpass with a cutoff near 16 kHz and a stopband slope approaching
-100 dB/decade. A "fix" that stabilises the circuit by detuning it into something else
fails this gate.
OUTCOME:

**Gate T0-3 (the cause is understood, not just suppressed).** A one-line statement of what
was wrong, with the line numbers. Adding damping until the poles move left does not
qualify: if the cause is not identified, the same defect returns in the next fixture.
OUTCOME:

**Gate T0-4 (record the old poles).** The pre-fix pole locations go in the conclusions
document. Historical numbers in git history need to remain *explainable*, and
"the circuit had RHP poles at +1.4491e+05 and +5.6716e+04" is what explains them.
OUTCOME:

### Stage T1 — decide the fate of the old fixture

**Decision needed from the maintainer, not to be guessed.** Options: replace the fixture
outright; or keep the unstable variant reachable behind a flag for reproducing historical
numbers. Recommendation: **replace outright**, and rely on T0-4 plus git history for
explicability — a fixture with a "deliberately broken" mode invites accidental use, and
nothing in the project needs to reproduce a divergent circuit.
OUTCOME:

### Stage T2 — the cheap experiments (symbolic and frequency-domain)

Re-run, in this order, recording each script's full output to a file:
`test_benchmark_circuits.py`, `cancellation_leapfrog.py`, `cancellation_blocks.py`,
`cancellation_parallel.py`, `cancellation_compose.py`, `transfer_function.py`,
`nonlinear_leapfrog.py`, `order_convergence.py`.

**Gate T2-1 (each script still passes its own internal gates).** Several carry their own
assertions — `nonlinear_leapfrog.py` has GATE 14-2 (symbolic matches numeric oracle to
1e-10) and GATE 14-3 (the harmonic is non-zero). These must still pass. GATE 14-3 matters
especially: the µA741 version of that test once recorded driving the wrong node and reading
an identically-zero harmonic, which produced entirely plausible sizes and timings for an
analysis that did nothing.
OUTCOME:

**Gate T2-2 (qualitative conclusions, held or refuted).** For each, state explicitly
whether the *conclusion* survives, not merely what the new number is. The ones at risk:
- the "agrees from U^5 / U^7 / U^9" order-convergence pattern
- the finding that required order tracks `% of v_turn` rather than drive level
- the 181-group hierarchical result and its "uninterpretable" verdict
- the cancellation `kappa` measurements and the `N_eff` concentration profile
- the "readable H(s) over 2 decades only" limit
A refuted conclusion is a result and gets recorded as one.
OUTCOME:

**Gate T2-3 (v_turn is recomputed, not carried over).** `v_turn = 1/sqrt(3|a|)` with
`a = kk/g`, and `g` is the driving-point conductance at the nonlinear node — which the
topology change moves. Every table row must carry its recomputed `% of v_turn`. A stale
`v_turn` is exactly how a model-caused divergence gets misread as a series-caused one,
which already cost a withdrawn table row once.
OUTCOME:

### Stage T3 — rewrite the pasted tables

Regenerate every table in the documents listed above from the T2 outputs. Do not hand-edit
numbers; paste whole generated blocks, and say in each section that it was regenerated
after the topology fix, with the date.

**Gate T3-1 (no stale numbers survive).** Grep the conclusions documents for the
distinctive old values (e.g. `5.9193e-11`, `4.6643e-06`, `1.8465e-02`, `366%`) and confirm
zero hits outside explicitly-historical passages.
OUTCOME:

**Gate T3-2 (doc build).** "build succeeded, 2 warnings", and `grep -cE 'ERROR'` checked
separately. Then verify the live `distortion_ddd.rst` block actually executed by grepping
the built HTML for a computed digit string that appears nowhere in the source — a failed
block renders its own source, and every heading and identifier is then present whether it
ran or not.
OUTCOME:

### Stage T4 — the transient IM3 comparison

Doubly blocked: needs the topology fix (T0) **and** the transient engine repair
(`transient_repair_plan.md` stages 1-5). Running it before either produces a number that
is growth or integration error rather than distortion.

This is stage 6 of `transient_repair_plan.md`; its gates 6-1 through 6-4 apply unchanged
and are not restated here.
OUTCOME:

## Order and dependencies

```
instability fix ──> T0 ──> T1 ──> T2 ──> T3
                     │
engine repair 1-5 ───┴──────────────────> T4
```

T2 and T3 need only the topology. T4 needs both tracks. The engine repair is already in
flight and is deliberately independent: the engine must be trustworthy before it is used to
judge a circuit, and fixing both at once would leave neither verified.

**Do not `git push`.** 58 commits are already unpushed; that is the maintainer's call.
