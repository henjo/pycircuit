# Fixing the leapfrog topology and redoing every experiment

**Status, 2026-07-30: T0 RUN and passed (all five gates). T1 settled by decision. T2 RUN
TWICE — once against the repaired 127-unknown fixture and again against the COMPENSATED
136-unknown one after `95545e5`; all seven consumers each time, two failures recorded (one
pre-existing). T3 RUN TWICE, the second pass covering the compensated numbers. T4 remains
BLOCKED and unattempted: it needs the
transient engine repair, which is complete (`transient_repair_plan.md` stages 1-5), so it
is now unblocked on both counts but has not been run.**

**The one-line summary of T2, because it is the reusable part:** what was refuted measured
the determinant's *conditioning*; what survived measured the *method*. A topology change
alters matrix conditioning and leaves the properties of the series and of the diagram
representation alone.

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

### Documents that LOOK self-updating but are NOT

- `doc/src/circuit/distortion_ddd.rst:500-560` — a live `exec-rst` block calling
  `bc.leapfrog_5th_order()`. This was listed here as self-updating. **It is not, and
  that was measured on 2026-07-30, not argued.** A build after the topology fix
  succeeded with the clean 2 warnings and 0 ERRORs, the block genuinely executed (the
  HTML holds a rendered table, and it did not fall back to echoing its source) — and
  the table still showed the numbers from the **unstable** fixture: 6.329e-05 /
  5.919e-11 where recomputing gives 7.1404e-05 / 6.2704e-11.

  The cause is sphinx's doctree cache: a source file that has not itself changed is not
  reprocessed, so its live block does not re-run. **A live block is self-updating with
  respect to its own source only, never with respect to the code it calls.** Stage T3
  must therefore force a clean rebuild (`sphinx -E`, or delete `doc/build/doctrees`)
  and then *diff the rendered numbers*, because "build succeeded, 2 warnings" is
  emitted identically whether the block re-ran or not.

  The separate, previously-known failure mode still applies too: a block that raises
  renders its own source while the build still reports success, so verification must
  grep the built HTML for a *computed* value, never for a heading or an identifier.
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
OUTCOME: **PASSED, 2026-07-30, on four routes that agree.** `rank(C) = 125` of 127, so 125
finite poles and 2 infinite, before and after.

| route | before | after |
|---|---|---|
| generalized eig, `scipy.linalg.eig(G, -C)` (LAPACK `ggev`/QZ) | 2 RHP: +1.44910e+05, +5.67162e+04 | 0 RHP, max Re = **-4.81494e+03** |
| shift-invert to a standard eigenproblem, `numpy.linalg.eigvals` of `-(G+sigma*C)^-1 C`, at three shifts sigma = 0, +1e4, -3e5 | same 2 RHP at every shift | 0 RHP at every shift |
| argument principle: winding of `det(G+s*C)` on a D-contour enclosing the RHP out to \|s\|=1e13 | **+2** RHP / 123 LHP | **0** RHP / 125 LHP |
| driven transient (implicit Euler, h = 2 ns, 1 mV at 1 kHz) | envelope grows at +1.45238e+05 /s, 0.2 % from the eigenvalue; peak reaches 3.7e+06 V by 300 us | peak settles to 5.0017e-04 V, i.e. the 0.5 passband gain of a 1 mV drive |

Conditioning: `cond(G) = 5.65e+05`, `cond(G + sigma*C) = 5.5e+05` to `8.7e+05` over the
shifts — no ill-conditioning to explain the result away. The two pre-fix RHP eigenvalues
were **numerically solid**, not artefacts: relative residual `||A(s)x||/(||G||+|s|||C||)`
of 2.8e-16 and 2.0e-16, and they moved by less than 0.2 in the 5th significant figure under
random *relative* entry perturbations of 1e-12, 1e-9 and 1e-6. Likewise the post-fix
`max Re = -4815` is stable to 1e-6 perturbation.

Two recorded traps in the measurement itself:

- An **undriven** transient of a linear circuit started at equilibrium returns exactly
  0.0000e+00 at every run length, because zero is an exact fixed point in floating point.
  It refutes nothing. Use a small nonzero drive or a perturbed initial condition.
- The argument-principle contour must not be a **circle tangent to the imaginary axis at
  the origin**: one parameter step beside the tangent point moves `|s|` by `r*dtheta`, which
  for r = 5e12 jumps clean across every pole and the phase aliases. That version reported
  "0 RHP zeros" for the broken circuit and "56" for the fixed one — both garbage. The
  working contour samples the imaginary-axis segment logarithmically (max phase step
  0.295 rad) and its two halves sum to 125 in both cases, which is the check that catches it.

**Gate T0-2 (it is still the intended filter).** `|H|` across 100 Hz - 1 MHz must look like
a 5th-order lowpass with a cutoff near 16 kHz and a stopband slope approaching
-100 dB/decade. A "fix" that stabilises the circuit by detuning it into something else
fails this gate.
OUTCOME: **PASSED, with one disclosed non-ideality.** Post-fix: DC gain 0.50010 (the
doubly-terminated ladder's 1/2), flat to 0.4 % out to 5 kHz, -3 dB at 27.9 kHz — the
*same* -3 dB frequency as the ideal `L1..L5` ladder reference, and `1/(2*pi*R*C)` = 15.9 kHz.
Stopband slope **-104.0 dB/decade** across 100-200 kHz, against **-101.4** for the ideal
ladder itself; -106.4 over 200-500 kHz. Deviation from the ideal ladder's `|H|` is at most
9.7 % anywhere from 100 Hz to 20 kHz. Nothing was detuned: R and C are untouched at
10 kOhm / 1 nF, and `dim` stays 127 with 536 nonzeros.

The disclosed non-ideality is a **+8.8 dB peak at 25.7 kHz**, from a Q = 16.8 pole pair.
It is not a residual sign error, and this was separated experimentally rather than assumed:

- The ideal all-equal-element ladder *already* has a Q = 5.87 pair at 26.6 kHz. That
  prototype is "every element 10k/1nF", not a designed approximation, so its pole set was
  never well damped.
- Rebuilding the fixed passive network around **nullors** instead of µA741s reproduces the
  ideal ladder's five poles exactly (-7.15225e+04, -5.0e+04 ± 8.66025e+04j,
  -1.42387e+04 ± 1.66615e+05j) with no peaking at all. So the topology is exactly the
  ladder and the peaking is the amplifier's.
- The µA741 model's gain is only ~31 at 16 kHz (GBW ≈ 0.5 MHz), and the excess phase that
  leaves in each integrator raises Q = 5.87 to 16.8. The same enhancement appears in *any*
  correct-sign realisation at this cutoff — including the negative-resistor sign probe,
  which peaks +9.8 dB.

Not "fixed", deliberately: it would take either a faster amplifier — and `ua741()` is the
calibration fixture, which must not move — or a lower cutoff, which is the detuning this
gate forbids. Recorded in the fixture docstring instead.

**Gate T0-3 (the cause is understood, not just suppressed).** A one-line statement of what
was wrong, with the line numbers. Adding damping until the poles move left does not
qualify: if the cause is not identified, the same defect returns in the next fixture.
OUTCOME: **PASSED.** One line: *the backward coupling entered the same summing node as the
forward coupling, so every stage integrated the sum of its neighbours where the ladder
integrates their difference, making each two-integrator loop positive feedback.*

In `pycircuit/circuit/benchmark_circuits.py` as it stood at `3fe5468` (the fixture is
byte-identical there to `e37ddad`; only docs and `nonlinear_leapfrog_sweep.py` moved in
between):

- **line 693** — `cir['rb%d' % k] = R(amp[k+1]['out'], amp[k]['inn'], r=10e3)`, the backward
  coupling, lands on `inn`;
- **line 687** — `cir['rf%d' % k] = R(amp[k-1]['out'], amp[k]['inn'], r=10e3)`, the forward
  coupling, lands on the *same* `inn`;
- **line 669** — `cir['sg%d' % k] = R(amp[k]['inp'], gnd, r=1.0)` grounds every
  non-inverting input through 1 Ohm, so there was no second input available to carry the
  opposite sign.

The amplifier was **not** at fault, and that was checked rather than assumed: driving the
µA741 model at `inn` gives arg H = +174 deg at 1 Hz and +90 deg at 10 kHz, at `inp` -5.8 deg
and -90 deg, so `inn` is genuinely the inverting input. `ua741()` alone has 0 RHP poles.

Two decisive probes:

- Replacing the four `rb` resistors with **negative** 10 kOhm — unrealizable, a pure sign
  probe — removes both RHP poles with no other change.
- The **nullor** rebuild of the *old* passive network still has 2 RHP poles
  (±6.18034e+04, ±1.61803e+05 with a -2e+05), which places the defect in the passive
  network and rules the transistor model out entirely.

Loop gain of a two-integrator loop is `-1/(s*T)**2` when the signs alternate, giving the
lossless `s = ±j/T` an LC resonance needs; with both couplings on one node it is
`+1/(s*T)**2`, and `1 - L = 0` at `s = ±1/T` — one of which is in the RHP. Note the poles
are invariant under relabelling the stage variables, so no sign convention repairs it; the
inversion must be physically in the network. `T = 10e3 * 1e-9`, so `1/T = 1e+05`, the order
of the two poles that were found.

Combinatorially there are only two ways to place the inversion: every forward coupling on
the non-inverting input, or every backward coupling. (Requiring one inversion per loop
*and* the two inputs of each stage on opposite sides forces `f_{k+1} = f_k` for all k, so
no mixed scheme exists.) Both were built; only the backward-coupling one gets the
terminations right, and it is also the smaller diff.

**Gate T0-4 (record the old poles).** The pre-fix pole locations go in the conclusions
document. Historical numbers in git history need to remain *explainable*, and
"the circuit had RHP poles at +1.4491e+05 and +5.6716e+04" is what explains them.
OUTCOME: **Recorded here and in the fixture docstring; NOT yet in
`doc/distortion_ddd_conclusions.md`, which is Stage T3 and has not run.**

Pre-fix (`HEAD = e37ddad`), 125 finite poles:

- 2 in the RHP: **s = +1.44910e+05** and **+5.67162e+04**, both real to within an imaginary
  part of ~1e-10 (growth time constants 6.90 us and 17.6 us).
- 123 in the LHP; slowest stable tau 1.7069e-05 s, fastest 8.6612e-13 s, stiffness 1.97e+07.
- `|s|` spanned 5.6716e+04 to 1.1546e+12.

Post-fix, 125 finite poles, all in the LHP: worst **Re = -4.81494e+03**; slowest tau
2.0769e-04 s, fastest 8.6628e-13 s, stiffness 2.40e+08 (up 12x — the dominant pole pair is
now the lightly damped Q = 16.8 one, which matters for transient step control). Dominant
poles: -6.90499e+04; -4.52474e+04 ± 8.53886e+04j (Q = 1.07); three near -1.0e+05 (the
`cp`/`rb` networks, pole-zero cancelled in `H`); -4.81494e+03 ± 1.61356e+05j (Q = 16.8);
-1.99845e+05.

**Gate T0-5 (added: the hole that hid it is closed).** Not in the original plan, added
because it is the more important half. `test_leapfrog_is_a_fifth_order_lowpass` **passed
for both broken versions of this fixture** — a magnitude response is not a stability test.
The sign error rotates the pole set by roughly 90 degrees, which leaves the DC gain and the
asymptotic slope alone, and the two passband points it samples (10 Hz, 1 kHz) sit two
decades below where the shape differs.
OUTCOME: **PASSED.** `test_leapfrog_has_no_right_half_plane_poles` added to
`pycircuit/circuit/tests/test_benchmark_circuits.py`. It asserts `rank(C) == 125`, exactly
125 finite poles, and `max Re < -1e3` — the margin, not just the sign, so a pole at
Re = -1 cannot pass. Verified in **both** directions: it fails on the stashed old fixture
with `2 poles in the right half plane, worst Re = +1.44910e+05`, and passes on the new one.
The two pre-existing leapfrog tests pass on the old fixture, which is the evidence that
they were never the guard.

### Stage T1 — the fate of the old fixture: SETTLED

**Maintainer's decision, 2026-07-30: replace the fixture outright. Historical tables do
not need to be re-runnable.**

So there is no flag, no `unstable=True` keyword, and no preserved variant. A fixture with a
"deliberately broken" mode invites accidental use, and nothing in the project needs to
reproduce a divergent circuit.

What this obliges instead, since reproducibility is being given up deliberately rather than
by oversight:

- **Gate T0-4 is now load-bearing, not a nicety.** The old pole locations
  (+1.4491e+05, +5.6716e+04) are the only remaining explanation for every historical number,
  so they must be recorded in the conclusions document, not just in a commit message.
- Historical tables get marked as **superseded and not reproducible at this HEAD**, with the
  commit that changed the fixture named. "Superseded" without a pointer is how a stale table
  gets quietly re-cited.
- No stage below may compare a new number against a historical one as though both came from
  the same circuit. They did not. Where a before/after comparison is genuinely wanted, it
  has to be framed as "the unstable fixture gave X" — a statement about a different circuit.

OUTCOME: settled by decision, no measurement required.

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
OUTCOME: **PASSED with two recorded failures, one of them pre-existing.** All seven
scripts re-run sequentially 2026-07-30.

| script | s | result |
|---|---|---|
| `cancellation_leapfrog` | 172 | ran clean |
| `cancellation_blocks` | 11 | **rc=1**, GATE 2 assertion — see below |
| `cancellation_parallel` | 58 | ran; **GATE 1b FAILS at nominal gm** (rel 6.52e-10 vs 1e-10); degraded passes at 5.40e-12 |
| `cancellation_compose` | 607 | ran; GATE 5 PASS 6/6, GATE 4 FAIL 4/6, GATE 6 FAIL 3/6 |
| `nonlinear_leapfrog` | 277 | **GATE 14-2 PASS (rel 6.52e-13), GATE 14-3 PASS** |
| `transfer_function` | 882 | ran clean |
| `order_convergence` | 327 | ran clean, extended to U^17 |

`nonlinear_leapfrog`'s two gates are the ones this stage most needed, and both hold —
the harmonic is 3.903694e-18, not the identically-zero trap. The symbolic build is also
still ~60x faster than the numeric oracle (1.2 s vs 73.4 s at U^3, 2.8 s vs 190.2 s at
U^5, i.e. **61-68x**).

**`cancellation_blocks` is NOT a regression, and checking that mattered.** It fails
`assert not stamps` -- block 1 picks up `_lvl0_*` symbols. It fails **identically against
the pre-repair fixture** (same assertion, same five symbols), verified by checking out
`ff5c6e6^`'s `benchmark_circuits.py` and re-running. It is the known `_suppress`
symbol-laundering defect that `cancellation_parallel.py` was written to avoid -- which is
why the parallel script's equivalent gate passes on the same partition. It was very
nearly reported here as caused by the repair.

`cancellation_parallel`'s GATE 1b failure IS attributable to the repair: precision loss
that tracks the `kappa` blow-up recorded under T2-2, and it appears only at nominal gm
where `kappa` is largest.

**Gate T2-2 (qualitative conclusions, held or refuted).** For each, state explicitly
whether the *conclusion* survives, not merely what the new number is. The ones at risk:
- the "agrees from U^5 / U^7 / U^9" order-convergence pattern
- the finding that required order tracks `% of v_turn` rather than drive level
- the 181-group hierarchical result and its "uninterpretable" verdict
- the cancellation `kappa` measurements and the `N_eff` concentration profile
- the "readable H(s) over 2 decades only" limit
A refuted conclusion is a result and gets recorded as one.
OUTCOME: **DONE. The results split cleanly along one line.**

| conclusion | verdict |
|---|---|
| leapfrog top diagram benign, `kappa = 13.8` | **REFUTED** -- `kappa = 1.153e+12` |
| the 181-group result | **rescaled 33x** -- 5 997 groups, 95 951 ops |
| "uninterpretable" verdict on those groups | **holds** (more so) |
| 16 ops/group is a property of the diagram | **holds** -- 16.0 |
| `N_eff` / concentration profile | superseded with `kappa`; same cause |
| composition not error-controlled | **holds**, refined: 2 of 6 settings do hold |
| symbolic nonlinear analysis matches its oracle | **holds** -- rel 6.52e-13 |
| ~60x faster than numerics | **holds** -- 61-68x |
| required order tracks `% of v_turn`, not drive | **holds** -- thresholds unchanged |
| readable `H(s)` over 2 decades only | **holds** -- 177 ops narrow, wide band fails |

**What was refuted is what measured the DETERMINANT'S CONDITIONING; what survived is
what measured the METHOD.** A topology change alters the conditioning of the matrix and
leaves the properties of the series, and of the diagram representation, alone. That is
the cut to read any future claim in these documents against.

The refutation has a physical cause rather than being a bare number change: the repaired
circuit has a **`Q ~ 16.8` pole pair (+8.8 dB at 25.7 kHz)**, disclosed when the topology
was fixed and not removable without a faster amplifier than `ua741()` or a lower cutoff.
A high-Q resonance is a near-cancellation between stages, so a large `kappa` is that
non-ideality's arithmetic signature. `cancellation_parallel` localises it: every block is
still benign (`D_k` `kappa = 1.147e+03`, cofactors median `2.827e+02`), so the added
cancellation is **entirely in the coupling between amplifiers**.

**Control that makes the deltas trustworthy:** `transfer_function`'s µA741 half is
untouched by the repair and re-measures **bit-identical** -- 734 groups -> 50 377
operations, and the 177-op two-decade `H(s)` at err 7.903e-02. Numbers that should not
have moved did not.

**Gate T2-3 (v_turn is recomputed, not carried over).** `v_turn = 1/sqrt(3|a|)` with
`a = kk/g`, and `g` is the driving-point conductance at the nonlinear node — which the
topology change moves. Every table row must carry its recomputed `% of v_turn`. A stale
`v_turn` is exactly how a model-caused divergence gets misread as a series-caused one,
which already cost a withdrawn table row once.
OUTCOME: **PASSED.** Recomputed, not carried over: `g` at `s0_e1` = **4.3184e-04 S**,
`a = kk/g` = 1.1578e+02, **`v_turn` = 5.3656e-02 V** (was ~5.43e-02 -- the driving-point
conductance moved, as predicted). Every row of the regenerated table carries its own
recomputed `% of v_turn`, and the `kk` sweep still flags `kk = 50` as **419% -- UNPHYSICAL,
past v_turn**, so the guard that cost a withdrawn row is still doing its job.

### Stage T3 — rewrite the pasted tables

Regenerate every table in the documents listed above from the T2 outputs. Do not hand-edit
numbers; paste whole generated blocks, and say in each section that it was regenerated
after the topology fix, with the date.

**Gate T3-1 (no stale numbers survive).** Grep the conclusions documents for the
distinctive old values (e.g. `5.9193e-11`, `4.6643e-06`, `1.8465e-02`, `366%`) and confirm
zero hits outside explicitly-historical passages.
OUTCOME: **PASSED.** `5.9193e-11`, `4.6643e-06`, `1.8465e-02`, `6.3290e-05`,
`5.8261e-07`, `7.6201e-03`, `5.1128e-06` -- all **CLEAN** outside this plan's own
deliberate historical references.

`2 895` and `= 13.8` still appear, in two places and both intentionally: the before/after
comparison tables written into `cancellation_ranking_conclusions.md`, where the old value
is the point. **The sweep did find one document I had not touched** --
`doc/cancellation_ranking_plan.md` carried an un-annotated "the leapfrog's top diagram
barely cancels (`kappa = 13.8`) and needs only 181 groups", now marked WITHDRAWN with the
repaired figures. That is precisely what this gate exists to catch: the rewrite worked
from the sections I remembered, and the grep found the one I did not.

### Second pass, 2026-07-30 — the compensated 136-unknown fixture

**PASSED, and the sweep earned its place a second time — differently.**

The obvious hits (`5 997`, `95 951`, `2 264 448`, `16.0 ops`, `374.3`, `7.1404e-05`,
`5.3656e-02`, `419%`) all resolved to **deliberate** occurrences: the middle column of
the new three-generation comparison tables, and this plan's own first-pass OUTCOME
records. `2.0668e-02` came back clean.

**But assuming that was the whole story would have been wrong.** Checking each hit's
context rather than trusting the classification found **two genuinely stale passages**,
both of them PROSE I wrote during the first pass:

- `cancellation_ranking_conclusions.md` — the verdict table still read *"rescaled 33x —
  5 997 groups / 95 951 ops"* and *"16 ops/group **holds** — 16.0"*. Now 46x /
  7 913 / 134 520, with ops/group recorded as a **rescale** (16 -> 16.0 -> 17.0), not a
  survival.
- `doc/distortion_ddd_conclusions.md` — a paragraph asserting the result *"survived the
  fixture repair"* (singular) with `v_turn` at 5.3656e-02 and node voltages "rose ~12%".
  Now covers both changes.

**The transferable lesson: the stale text is not in the tables, it is in the sentences
that interpret them.** A numeric sweep finds tables reliably; it finds prose only where
the prose happens to quote a number. "Survived the fixture repair" is stale in a way no
grep detects, and it was wrong in a document whose subject is what survived.

Also confirmed as *correctly* untouched: `transient_repair_plan.md`'s copies of
`7.1404e-05` and `6.2704e-11`, which are historical evidence of the stale-doctree trap —
they record what the built page showed at that moment, and updating them would destroy
the evidence.

**Gate T3-2 (doc build).** "build succeeded, 2 warnings", and `grep -cE 'ERROR'` checked
separately. Then verify the live `distortion_ddd.rst` block actually executed by grepping
the built HTML for a computed digit string that appears nowhere in the source — a failed
block renders its own source, and every heading and identifier is then present whether it
ran or not.
OUTCOME: **PASSED, and only because the rebuild was forced.** `rm -rf build/doctrees`
plus `sphinx -E`: **build succeeded, 2 warnings**, `grep -cE 'ERROR'` = **0**.

The live block regenerated against the repaired fixture, verified by computed digit
strings that appear nowhere in the source: `7.140e-05`, `6.270e-11` and `2.067e-02` are
each present exactly once, and the pre-repair `6.329e-05` is **gone**. An ordinary
`sphinx -b html` earlier that same day had left every one of the old values in place
while still reporting "build succeeded, 2 warnings" -- see the trap recorded under
"Documents that LOOK self-updating but are NOT".

**A first-decimal mismatch nearly hid this.** The block prints 3 decimals where
`order_convergence.py` prints 4, so grepping the built page for `7.1404e-05` returned
zero hits -- as did the old `6.3290e-05`. Zero hits for both old and new reads exactly
like "the page is stale in some third way". The check has to be done in the *block's own*
output format, which means reading the rendered table rather than guessing at it.

**One real defect found by reading it rather than grepping it:** the header column
rendered as `%% of v_turn`. The block emits `print("     - %% of v_turn")`, a plain
string with no %-formatting, so the escape survives literally; the neighbouring lines
are %-formatted, which is how it was copied in. Fixed to `percent of v_turn`. Cosmetic,
but it is only visible in the rendered output -- a grep for computed values would never
have found it.

### Second pass, 2026-07-30 — the compensated 136-unknown fixture

**PASSED.** `rm -rf build/doctrees` + `sphinx -E`: **build succeeded, 2 warnings**,
`grep -cE 'ERROR'` = **0**.

The live block regenerated against the compensated fixture, and the check was made in
**both directions** — new values present AND both older generations absent, which is the
form that cannot be satisfied by a stale page:

| value (block's own 3-decimal format) | hits |
|---|---|
| `5.774e-05`, `5.756e-11`, `1.691e-02` (compensated 136) | **1 each** |
| `7.140e-05` (repaired 127) | **0** |
| `6.329e-05` (unstable 127) | **0** |

Checking only for the new values would have been weaker: zero hits for *both* old and new
is a real state (it happened on the first pass, from a decimal-format mismatch) and reads
exactly like staleness. Requiring the old ones to be gone distinguishes them.

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
