# Can the DDD work be reused for symbolic distortion? — feasibility

**Verdict: yes, and the fit is closer than it first appears — but not for the
reason that motivated the question, and the largest available win may not need
DDD at all.**

Written before any code, as the reasoning behind
`doc/distortion_ddd_plan.md`. Everything numeric here was measured, not
argued; the spikes are reproducible from the snippets quoted.

> ## SUPERSEDED, AND NOT REPRODUCIBLE AT THIS HEAD — leapfrog numbers only
>
> **Every number in this document measured on `leapfrog_5th_order` was measured on an
> UNSTABLE circuit, and the fixture has since been replaced.** That is sections 10.1-10.3
> in particular, including the corrected amplitude/`v_turn` table and the `kk` sweep table.
>
> The fixture had **two right-half-plane poles, s = +1.4491e+05 and +5.6716e+04**: the
> backward coupling resistors entered the same summing node as the forward ones, so every
> stage integrated the *sum* of its neighbours where a ladder simulation integrates their
> *difference*, making each two-integrator loop positive feedback. Confirmed on four
> independent routes; see `doc/leapfrog_redo_plan.md`, Gate T0-1.
>
> **What this does and does not invalidate.** `H(s)` on the jw axis is well defined for a
> divergent circuit — the analysis just solves a linear system at `s = jw` — so these
> numbers are *arithmetically correct*. What is void is the physical claim: they describe a
> circuit that has no steady state, so any statement of the form "a designer could read
> this and understand the filter" is about a filter that does not exist.
>
> The topology fix landed in `pycircuit/circuit/benchmark_circuits.py` on 2026-07-30
> (`_build_leapfrog`: the four `rb` resistors moved to the non-inverting input, matched
> `cp0..cp3` capacitors and `rp0` added, `sg0..sg3` removed). By the maintainer's Stage T1
> decision the unstable variant was **replaced outright with no flag**, so the tables below
> cannot be regenerated at this HEAD — the old poles above are their only explanation.
> Regenerating them is Stage T2/T3 of `doc/leapfrog_redo_plan.md` and **has not run.**
> Do not compare any new leapfrog number against a table here as though both came from the
> same circuit. They did not.

## 1. Where this came from

`doc/distortion_mimo_plan.md` §8.3 records an open question. Stage E measured
symbolic expression **size** and found it polynomial in truncation order, with
multi-node costing a small constant factor. It measured nothing about **cost**,
and `sympy.cancel` failed to finish in 900 s at `U^7` on a two-node system.
Those are different claims and only one was established.

The suggestion to look at the DDD work is a good one because the two problems
share a structural property that is unusually specific.

## 2. The structural match

The distortion method's central claim, from the theory page: **every step is a
solve against the same linear operator `Y(s) = G + sC`, evaluated at a
different frequency.** Nothing needs a nonlinear solve or a new matrix.

The DDD work's central objects are the answer to almost exactly that:

| distortion needs | DDD provides |
|---|---|
| `det(G + sC)` at many harmonic frequencies | `s_expand` — determinant split by powers of `s`, **coefficients shared as diagrams**; one construction serves every harmonic |
| one matrix, very many right-hand sides (one per graded `(harmonic, power)` key) | `DDDFamily` — numerators expressed as **cofactors of the original matrix**, so a new right-hand side is a *recombination*, not a rebuild |
| many evaluations of overlapping structure | `eval_roots` — one pass over shared structure, memoised across roots |

The `DDDFamily` docstring states its own motivation as noise analysis, where
"one transfer function per noise source costs little more than a single
transfer function: the transimpedances share almost everything." **A graded
perturbation solve has the same shape**: dozens of right-hand sides against one
matrix, at a handful of distinct frequencies.

## 3. What was measured

Three spikes on the gm-C biquad pencil `[[-g2 + s*C1, -g4], [-g3, s*C2]]`,
symbolic entries, evaluated against `numpy.linalg.solve`.

**Complex `s` works.** This was the first thing that could have killed it: DDD
is built for a symbolic frequency variable, and harmonics need `s = j m w`.

```
m=1  max rel err 2.50e-15
m=2  max rel err 2.13e-16
m=3  max rel err 1.90e-16
```

**`s_expand` reproduces the determinant across harmonics** from one
construction: degree 2, **5 shared vertices**, `det` correct to 5e-15 at
`m = 1` and `m = 3`.

**The right-hand side really is free, but only on the fast path.** This is the
important measurement:

| route | cost per right-hand side |
|---|---|
| `DDDCombination.eval(env)` per solve | **3.05 ms** |
| `eval_roots` memo once per frequency, then recombine cofactors | **3.07 us** |

**A factor of 1000, for the same answer to 2.5e-15.** The DDD memory records
this trap as a 200x `subs` cost fixed by a `_resolve` fast path and keeping
factors unmultiplied; here the same trap is worse. **Any prototype that reaches
for `eval` in a loop will measure DDD as hopeless and be wrong.**

## 4. What DDD does *not* solve — the honest part

**The blow-up in §8.3 is not in the determinant.** After `p` passes of the
graded recurrence, a coefficient is a sum of products of up to `p` cofactor and
determinant factors *at different frequencies*. DDD compacts each factor and
shares them; it does not by itself stop the product tree growing. Anyone
expecting "use DDD and the symbolic cost problem goes away" will be
disappointed.

**DDD is not a numeric accelerator, and must not be sold as one.** Even the
fast path is ~3 us per solve against numpy's sub-microsecond dense solve on a
2x2. The existing numeric path must stay numpy. Nothing in this work should
make the current tests faster, and a claim that it does would be false.

**The 2x2 biquad proves nothing about scale.** Five vertices is trivial.
Whether sharing pays on a matrix where it matters is unmeasured — and that is
the whole question.

## 5. So what would actually help

Ordered by expected value per unit of work, cheapest first. **This ordering is
the main output of this document.**

1. **Stop calling `cancel`, and keep coefficients factored.** The single
   measured fact in §3 is that a memo plus recombination beats substitution by
   1000x. That lesson is *transferable without DDD*: represent each graded
   coefficient as a sum of products of unevaluated factors, and evaluate
   through one memo. **It is entirely possible this alone answers §8.3.**
2. **DDD for the linear solves.** Replace the `solve` callable with a
   DDD-backed one, so every harmonic draws on one shared construction.
3. **Canonical sharing across the whole computation.** DDD's distinctive
   contribution over (1) is that identical minors recurring across harmonics
   *and* across perturbation orders collapse to one vertex. sympy has no
   canonical merge. This is the part that could beat a plain factored
   representation, and the part nobody has evidence for.

**Reconsider-if for the whole exercise:** if stage A of the plan shows that a
plain factored graded ring already makes `U^7` tractable symbolically, then
DDD is an optimisation of an already-solved problem and this work should be
stopped and recorded as such. That outcome is a success, not a failure.

## 6. Why this is a good place to try it

The DDD work's own calibration boundary is a real limitation: for named
op-amps "a published `|DDD|` pins an order of magnitude, not a number", because
formulation and device model differ between papers. Its *exact* checks were
confined to full matrices.

The distortion implementation is a much better oracle. It is gated against
**four published closed forms** across two circuits (eqs. 41, 42, 48, 52), an
**independent ODE integration** (magnitude, phase to 0.0000 deg, and full
waveform to 2e-10), and a **2-D numerical Fourier extraction** for the
two-tone exponential. A DDD-backed path can be required to reproduce all of it
to floating point.

**That is the strongest argument for doing this at all**: not that DDD will be
faster, but that for once there is an oracle good enough to prove a symbolic
representation correct rather than merely plausible.

## 7. Constraints and traps carried over

- **Use `eval_roots`, never `eval`, in any loop.** §3. This is not a
  micro-optimisation; it is the difference between feasible and not.
- **`DDDFamily` defaults to `'min-degree'`, not `'auto'`**, deliberately:
  `'auto'` reorders the matrix and the class addresses cofactors by *original*
  row and column. Do not "improve" this to `'auto'`.
- **The exact-rational blow-up that killed GiNaC cannot arise here** (DDD
  negative result P3): a diagram never multiplies entries symbolically. With
  complex float entries we are further from it still.
- **Keep factors unmultiplied as tuples**, as `DDDCombination` does — a
  premultiplied product is a sympy expression and costs a full substitution.
- **Separate files and a separate label throughout**, per the instruction that
  prompted this: the existing implementation must remain the reference, not a
  thing to be refactored around.

## 8. Which diagram for which regime — the framing that matters

Added after the first draft, prompted by two things: the question of whether a
*determinant* diagram is even the right shape, and the goal of running
distortion analysis on **bigger circuits**. The second inverts part of §4.

The cost has two independent sources, and they dominate in opposite regimes:

| | small circuit, high truncation order | big circuit, low order |
|---|---|---|
| dominant cost | the **product tree** of graded coefficients | the **determinant and cofactors** |
| right tool | a canonical *polynomial* diagram, or just factored evaluation | **DDD, exactly as built** |
| evidence | tree/DAG ratio 1.6 → 4.7 over `U^3`..`U^9`; `cancel` unfinished at `U^7` on 2 nodes | µA741 flat `\|DDD\|` 1040, hierarchy 1040 → 156, noise 11 088 → 26 vertices |

**§4 said "the blow-up is not in the determinant". That is true of the 2x2
biquad and false of a 25x25 op-amp.** On the biquad the determinant has five
vertices and the products dominate; on a real amplifier the determinant is the
millions-of-terms object DDD was built for, and the graded products sit *on
top of* it. For the stated goal — bigger circuits — DDD is not the wrong
shape. It is the necessary substrate, and the polynomial growth rides above it.

The two compose rather than compete: a graded coefficient on a large circuit
is a product of cofactors, each of which DDD represents compactly.

## 9. Other diagram variants, and an honest note on sourcing

**These are directions to evaluate, not results.** The project rule applies:
*only use external data if both the thing measured and the measurement are
understood*. None of the literature below has been re-read for this document
and **no number from it should be quoted until the paper is in hand**. What
follows is a shape-matching argument, which is all it claims to be.

### 9.1 The natural extension of what already exists — a *graded* diagram

`SExpandedDDD` is "one root per power of `s`, sharing subgraphs, memo key
`(rows, cols, power)`". The distortion work grades by `(harmonic, power of the
drive)`. So the direct analogue is **one root per graded key, memo key
`(rows, cols, harmonic, power)`** — the same multiroot construction Shi & Tan
use for `s`, applied to the grading this analysis already carries.

This is the cheapest variant to try because the machinery, the sharing
argument and the trimming of zero coefficients all exist; only the expansion
variable changes. **It is also the one most likely to pay**, because the
grading is exactly where the repeated structure lives.

### 9.2 Diagrams built for polynomials rather than determinants

A determinant is one specific polynomial — a signed sum over permutations —
and DDD's variable order and zero-suppression are tuned to it. The graded
numerator is a *general* multivariate polynomial. Two families were built for
that case:

- **Taylor Expansion Diagrams (TED)** — canonical form for multivariate
  polynomials by Taylor decomposition, `f = f(0) + x f'(0) + x^2 f''(0)/2 +
  ...`, each child a diagram. Developed for datapath verification, where BDDs
  blow up on arithmetic. The decomposition is the same one this analysis
  already performs on its devices, which is at least suggestive.
- **Multiplicative binary moment diagrams (`*BMD`, `K*BMD`)** — moment
  decomposition `f = f_0 + x f_1` with **multiplicative edge weights**. Those
  weights are the canonical version of the trick that gave 1000x in §3:
  keeping a common factor out of the representation instead of distributing
  it.

**Reconsider-if, and it is a real one:** both are canonical forms for
polynomials over a *fixed* variable set. Here each harmonic contributes its
own set of transfer-function values, so the variable count grows with the
truncation order. Whether either stays compact under that is exactly the
question, and it is unmeasured.

### 9.3 Considered and set aside

- **ZDD over factor multisets.** Terms are *multisets* — a factor can appear
  squared — and ZDDs are canonical for *sets*. Encoding multiplicity as
  distinct variables is possible but reintroduces the growth it was meant to
  remove. **Reconsider if** a formulation appears where every factor is
  square-free.
- **ADD / MTBDD.** Terminals carrying values help when many coefficients
  repeat exactly. Here they are distinct rational functions.

### 9.4 The lever for *bigger circuits*, which needs no new diagram at all

`HierarchicalDDD` and `suppression_order` already exist and already work:
1040 → 156 vertices on the µA741, 11 088 → 26 and 86 s → 0.48 s on noise.
Suppressing internal nodes before expanding is a property of **the matrix**,
and the distortion path solves the same matrix at every harmonic.

**So the most promising route to larger circuits may not be a new diagram
variant but the existing hierarchy applied to the existing solve** — cheaper
than anything in §9.2 and already validated on a circuit far larger than
anything the distortion work has run. It goes in the plan ahead of the exotic
options for that reason.

## 10. The circuit axis: fully symbolic nonlinear analysis of a 127-unknown filter

§4 above recorded the circuit axis as the open one — this representation had been
pushed on perturbation *order* and never tried at op-amp scale. These three sections
close it. **They were originally written in
`doc/cancellation_ranking_conclusions.md` §§10.1-29**, which is a document about term
ranking in a linear determinant; they are about perturbation order in a nonlinear
analysis, so they have been moved here. Plan stages: `doc/distortion_ddd_plan.md` §7.
Live table: `doc/src/circuit/distortion_ddd.rst`.

### 10.1 Feasibility: fully symbolic nonlinear analysis of the leapfrog — it works

A different question from the cancellation-ranking work (the cancellation-ranking work (`cancellation_ranking_conclusions.md` §§1-26), asked by the maintainer: nonlinear analysis of the
127-unknown leapfrog, fully symbolic, then evaluated numerically.
`benchmarks/nonlinear_leapfrog.py`.

**Measured, both orders passing every gate:**

| order | graph nodes | symbolic build | evaluate | **numeric oracle** | vs oracle |
|---|---|---|---|---|---|
| `U^3` | 11 945 | 1.4 s | 0.20 s | **87.1 s** | 4.63e-13 |
| `U^5` | 20 992 | 3.7 s | 0.32 s | **238.5 s** | 4.63e-13 |

Third harmonic `8.193583e-26`, non-zero (gate 14-3) and matching an independent
fully-numeric graded solve to **4.6e-13** (gate 14-2, which asked for 1e-10).

**The symbolic route is 60× faster than the numeric one.** At `U^5`, 4.0 s to build and
evaluate the graph against 238.5 s for the numeric path — because the graph is built
once and evaluated by a memoised walk, while the numeric path re-solves a 127×127
matrix for every `(harmonic, power)` key. That inverts the usual expectation that
symbolic costs more, and it is the strongest argument for this representation.

**Growth is mild.** `U^3 → U^5`: nodes 1.76×, build time 2.6×. Not the exponential the
determinant work fought throughout the cancellation-ranking work (the cancellation-ranking work (`cancellation_ranking_conclusions.md` §§1-26).

### Three things stated so the result is not over-read

**The leapfrog fixture is linear.** It is built from `add_small_signal_bjt`, which
stamps only `R`, `VCCS` and `C`, so it contains no nonlinearity of its own.
`graded_response_mimo` *attaches* one, and where it attaches is a modelling choice.
Made explicitly: stage 0's input differential pair (rows 4 and 5), because that is
where an op-amp's distortion originates. A different placement is a different circuit
and would give different numbers.

**"Fully symbolic" here means the expression-graph sense**, not the expanded sense.
`distortion_ddd.Expr` is a straight-line program over the device symbols — 20 992
nodes, hash-consed, never expanded. It evaluates fast and, as the `distortion_ddd`
documentation already records, **cannot rank terms or be read.** That is exactly right
for "symbolic form, numeric evaluation" and exactly wrong for the readability goal
the cancellation-ranking work (the cancellation-ranking work (`cancellation_ranking_conclusions.md` §§1-26) pursued. The two must not be conflated: this stage does *not* deliver a
readable nonlinear expression, and nothing here suggests one is available.

**`U^5` demonstrates scaling, not extra accuracy.** The third harmonic is identical at
both orders to every printed digit, because at a 1 mV drive and `k = 1e-9` the `U^5`
correction to the third harmonic is far below double precision. The `U^5` run is
evidence about *cost growth*; it is not evidence that the answer improved.

### What this says about the session as a whole

the cancellation-ranking work (the cancellation-ranking work (`cancellation_ranking_conclusions.md` §§1-26) asked for a *readable* symbolic expression and got 4 081 operations at best —
20× short. This stage asked for a *fast, exact* symbolic representation of a harder
(nonlinear) problem on a bigger circuit and got it comfortably, 60× faster than
numerics.

**The two goals want opposite representations**, which is the most useful thing to
carry away. An expanded, rankable form can be read and cannot scale; a hash-consed
straight-line program scales beautifully and cannot be read. Every stage of the cancellation-ranking work (the cancellation-ranking work (`cancellation_ranking_conclusions.md` §§1-26) was
spent trying to make the first behave like the second.

### 10.2 Fully symbolic nonlinear analysis of the leapfrog: order vs amplitude

`benchmarks/order_convergence.py`. 127 unknowns, cubic `i = kk*v^3` on stage 0's
input differential pair, drive at **100 kHz**, third harmonic at the output.

**Everything stays symbolic and numbers go in last** — including the drive amplitude
`Adrv` and the nonlinear coefficient `kk`. That is not a constraint but the economy of
the method: **one symbolic build per order serves the whole sweep.**

| order | graph nodes | build | evaluate (6 amplitudes) |
|---|---|---|---|
| `U^3` | 6 415 | 1.8 s | 0.13 s |
| `U^7` | 22 461 | 9.9 s | 0.69 s |
| `U^13` | **54 079** | **51.7 s** | **2.11 s** |

121 s of builds in total, then every amplitude and every `kk` is a graph walk.

### The table

| amp (V) | \|v\| at s0_e1 | HD3 | d(U^5) | d(U^7) | d(U^9) | d(U^11) | d(U^13) | **agrees from** |
|---|---|---|---|---|---|---|---|---|
| 0.01 | 6.33e-05 | 5.92e-11 | 1.6e-06 | 2.9e-12 | 0 | 0 | 0 | **U^5** |
| 0.03 | 1.90e-04 | 5.33e-10 | 1.4e-05 | 2.3e-10 | 4.2e-15 | 0 | 0 | **U^5** |
| 0.1 | 6.33e-04 | 5.92e-09 | 1.6e-04 | 2.9e-08 | 5.7e-12 | 1.2e-15 | 0 | **U^5** |
| 0.3 | 1.90e-03 | 5.32e-08 | 1.4e-03 | 2.3e-06 | 4.1e-09 | 7.7e-12 | 1.5e-14 | **U^7** |
| 1 | 6.31e-03 | 5.83e-07 | 1.6e-02 | 2.9e-04 | 5.8e-06 | 1.2e-07 | 2.6e-09 | **U^9** |
| 3 | 1.85e-02 | 4.66e-06 | 1.7e-01 | 2.7e-02 | 4.7e-03 | 8.8e-04 | 1.7e-04 | **not by U^13** |

`d(U^k) = |H3(U^k) − H3(U^k−2)| / |H3(U^k)|`. The required order climbs monotonically
with amplitude — `U^5` up to 100 mV, `U^7` at 300 mV, `U^9` at 1 V, and at 3 V the
series has not settled by `U^13`, its successive corrections falling only from 1.7e-01
to 1.7e-04. That is the textbook signature of approaching the radius of convergence,
and it is the answer to "up to what amplitude may I trust third-order theory?"

### The frequency mattered as much as the amplitude

A first run at **1 kHz** produced a useless table: every amplitude converged at `U^3`
with HD3 ~ 1e-14. Physically correct, and the reason is instructive — at 1 kHz the
amplifier's loop gain is high, so the input pair is held at a **virtual ground**
(2.5 µV at a 10 mV drive) and the cubic barely acts. **Distortion rises where loop gain
falls.** This µA741 measures +87.9 dB at 10 Hz falling at −20 dB/decade, so at 100 kHz
its open-loop gain is down to roughly 8 dB.

The evidence that the loop gain really has reduced is in the table: at the same 10 mV
drive the input-pair voltage is **6.33e-05 against 2.5e-06 — 25× larger.** The virtual
ground no longer holds, and the higher orders come alive.

An intermediate attempt moved the nonlinearity to the amplifier's *output* to find
signal swing. That works but is the wrong physics: an op-amp's distortion originates in
the input pair, and the right fix was the frequency, not the node.

### The free second sweep, which is the point of substituting last

`kk` is a symbol, so sweeping the nonlinearity's strength needs **no rebuild at all**:

| `kk` | HD3 at 1 V | agrees from |
|---|---|---|
| 0.005 | 5.91e-08 | `U^7` |
| 0.5 | 5.11e-06 | not by `U^13` |
| 50 | 7.62e-03 | not by `U^13` |

Two independent parameter sweeps, one build. Substituting numbers into the circuit
first would have cost a full rebuild per point — 121 s each rather than milliseconds.

### What was dropped, and why

The original request was to compare against **transient** simulation. That was
abandoned on a measured cost wall, not a guess: the leapfrog is stiff — time constants
from ~0.5 ns (the µA741's 27 Ω output degeneration into 20 pF) to 10 µs (10 kΩ × 1 nF
integrators), with 1 Ω input tie-downs — and **2 ms of simulated time did not complete
in 6 minutes**, while the table needed hundreds of times more. A contributing cost:
`BSource`'s Jacobian is a central-difference numerical derivative, two extra evaluations
per Newton iteration per step. The harness for it exists in
`benchmarks/nonlinear_leapfrog_sweep.py` — correct, with a real `BSource` cubic the
transient engine would see, and unaffordable at 1 kHz. Raising the drive frequency, as
was done here for a different reason, would also shrink that integration by ~100×;
that is the obvious way to revive it.

**So this table is self-convergence, not validation against an independent solver.** The
representation itself *was* validated independently — stage 14 matched a fully numeric
graded solve to 4.6e-13 on this circuit — so what is unverified here is the perturbation
series' own truncation, which is exactly what the table measures.

### 10.3 §10.2 corrected: one row was past the cubic's turning point

The maintainer pointed at a finding already recorded in `doc/distortion_mimo_plan.md`
§8.3 and not applied in §10.2:

> `i = g(v + a*v^3)` with `a < 0` turns over at `v_turn = 1/sqrt(3|a|)`, beyond which it
> has negative differential conductance and **is not a physical device**. And the
> *loaded* cell is worse than the isolated device, at 87% of the turning point rather
> than 66%: the same negative `alpha` that compresses the output current also reduces
> the node's loading, so the node voltage **expands** toward the breakdown.

§10.2 reported convergence orders without ever asking where its amplitudes sat relative
to that limit. Measured now: at 100 kHz the driving-point conductance at `s0_e1` is
`g = 4.46e-04 S` (`|Z_ii| = 2.24 kohm`), so with `a = kk/g`:

| case | `a` | `v_turn` | node voltage | % of `v_turn` |
|---|---|---|---|---|
| `kk` = 0.05, 3 V drive | 1.12e+02 | 5.45e-02 | 1.85e-02 | **34%** |
| `kk` = 0.05, 1 V drive | 1.12e+02 | 5.45e-02 | 6.31e-03 | 12% |
| `kk` = 0.5, 1 V drive | 1.12e+03 | 1.72e-02 | ~6.3e-03 | **37%** |
| **`kk` = 50, 1 V drive** | 1.12e+05 | **1.72e-03** | ~6.3e-03 | **~370%** |

**The `kk = 50` row of §10.2's second table is withdrawn.** Its node voltage is several
times `v_turn`, so the element has negative differential conductance there and is not a
device; "not settled by `U^13`" describes a model that has already stopped meaning
anything. A series cannot converge outside the region where its model is physical, and
reporting the non-convergence as a property of the *series* mistakes an artefact for a
result.

**The amplitude table stands, with its scale now stated.** The 3 V row sits at 34% of
`v_turn` — inside the physical region, and consistent with the recorded rule that 1 dB
compression lands at 66% and a loaded node breaks down at 87%. So "approaching the
radius of convergence" was right in direction; what §10.2 lacked was the number that makes
it checkable. Every row should carry its fraction of `v_turn`, and the script now
computes it.

**One thing to state that §10.2 left implicit: the sign.** `kk` here is **positive**, and
with the recorded convention `(A + sC)x = G_h x_in + f_h(x_in) - f(x)` that makes the
cubic **compressive** — confirmed in the data, since scaling the 0.01 V node voltage of
6.33e-05 linearly to 3 V predicts 1.90e-02 against the 1.85e-02 measured, a 2.6%
compression. So this run is not in the expansive regime the recorded finding warns
about; it is the compressive one, and it still loses convergence at a third of
`v_turn`.

**The transferable lesson, which cost a wrong table row:** a truncated polynomial
nonlinearity has a validity limit of its own, independent of the perturbation order, and
an order-convergence study is only meaningful inside it. `v_turn` should be computed and
printed beside any such sweep — it is two lines given the driving-point impedance, and
without it a divergence caused by the *model* is indistinguishable from one caused by
the *series*.

### The corrected table, with the scale on every row

```
amp (V)   |v| s0_e1    %v_turn   HD3          d(U^5)   d(U^7)   d(U^9)   d(U^11)  d(U^13)  agrees from
0.01      6.3290e-05   0%        5.9193e-11   1.6e-06  2.9e-12  0        0        0        U^5
0.03      1.8987e-04   0%        5.3273e-10   1.4e-05  2.3e-10  4.2e-15  0        0        U^5
0.1       6.3287e-04   1%        5.9184e-09   1.6e-04  2.9e-08  5.7e-12  1.2e-15  0        U^5
0.3       1.8981e-03   3%        5.3197e-08   1.4e-03  2.3e-06  4.1e-09  7.7e-12  1.5e-14  U^7
1         6.3080e-03   12%       5.8261e-07   1.6e-02  2.9e-04  5.8e-06  1.2e-07  2.6e-09  U^9
3         1.8465e-02   34%       4.6643e-06   1.7e-01  2.7e-02  4.7e-03  8.8e-04  1.7e-04  not by U^13

kk        HD3          %v_turn   agrees from
0.005     5.9098e-08   4%        U^7
0.5       5.1128e-06   37%       not by U^13
50        7.6201e-03   366%!     UNPHYSICAL -- past v_turn
```

**And the pattern is cleaner stated against `v_turn` than against amplitude.** The
required order tracks *how close the node is to the cubic's turning point*, not the
drive level as such: `U^5` below a few percent, `U^7` at 3-4%, `U^9` at 12%, and no
settling by `U^13` from about 34% onward — which the two independent sweeps agree on
(34% by amplitude, 37% by `kk`). That is one curve in one variable, and §10.2 had it
spread across two.
