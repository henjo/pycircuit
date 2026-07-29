# Perturbation distortion analysis — implementation plan

**Status: nothing implemented. Gates below are declared in advance.**

Companion documents:

- **Reasoning** — `~/pycircuit_claude/distortion_pertubation.md`. Why this
  method, what is in and out of scope, and every rejection with its
  reconsider-if condition. Read that first if you want to argue with the
  approach; this document only says *what will be done and in what order*.
- **Sources** — `~/pycircuit_agy/papers/distortion_pertubation/distrortion_pertybation_reference.md`.
  Per-paper extractions. Every equation and number there was read off a
  600-dpi render, never `pdftotext`, which silently dropped micro prefixes
  twice in this paper set (a factor of 10⁶ each time).

The goal is a new analysis family predicting harmonic and intermodulation
distortion of weakly-nonlinear *driven* circuits, by regular perturbation of
the linearised circuit with harmonic balance closing each order. Primary
reference is **Buonomo & Lo Schiavo, AICSP 77:483–493 (2013)** — hereafter
*the 2013 paper* — with the 2005 TCAS-I paper as the source for the Volterra
equivalence proof.

## 1. Ground rules

1. **The theory document ships in the same commit as its code**, not as a
   documentation phase afterwards. Writing the explanation is what forces
   there to be one, and it becomes an anchor to check the code against.
2. **No measured number is typed into prose.** Generate it at build time so
   it is re-measured on every build. A pasted number is correct exactly once.
3. **Gates are real.** A failed gate stops or redirects the work; it does not
   get quietly reinterpreted as a pass. A refuted premise is a result and
   gets recorded in §5 as one.
4. **Only extract a test where both the circuit and the measurement are
   understood.** Same rule that governed the DDD benchmark survey. Several
   worked examples in these papers publish results only as graph curves;
   those give tolerances in dB, not exact values, and must be labelled
   graph-derived wherever they appear.
5. **Never mix parameter sets across papers.** The 2009 and 2012 papers
   publish the same amplifier with fourteen differing values. A 2009 value
   checked against a 2012 plot yields a small, plausible, entirely spurious
   disagreement — the worst kind, because it looks exactly like an
   implementation bug. The 2013 paper settles it: use the **2009/2013** set.

## 2. Architectural integration — what this forces to change outside itself

Reviewed through a single named lens before any code, per the working
method: *as an architect, what does this force to change outside itself, and
how do I avoid bloating those parts?*

### 2.1 What the existing architecture already gives us

The method's whole computational content is: solve against the **same**
linear operator `G(s) = (A + sC)⁻¹`, evaluated at `s = j·m·ω₀` for
`m = 0,1,2,3` (or at `j(mω₁ + nω₂)` for two tones). There is no nonlinear
solve at any order.

That is *exactly* what `AC` already assembles. So:

- **`circuit.py` needs no change.** The MNA stamp, the `A`/`C` matrices, and
  the reference-node handling are all already there and already correct.
- **The solve is an existing toolkit primitive.** A symbolic `G(s)` held once
  in `s` serves every harmonic — the same lever DDD's H7 used, where
  overriding two primitives moved a whole analysis family without rewriting
  any analysis.
- **`acsolution.py` already defines the result-object contract**
  (`ACSolution`, `NumDenSolution`), and `dddresult.py` shows how a second
  representation slots in beside it.

### 2.2 Where it genuinely presses outside itself

**One place: the element layer must yield second and third derivatives.**
The method needs the Taylor coefficients `b` (quadratic) and `c` (cubic) of
each nonlinearity about the operating point. Elements today expose `i()`/`q()`
and `G()`/`C()` — first derivative only.

### 2.3 The bloat problem, and the fix

The tempting design is to give every nonlinear element a new API — some
`taylor_coefficients()` method. **Reject it.** There are twenty-plus elements;
it would mean touching all of them, it creates a second statement of each
device model that can silently disagree with `i()` (precisely the P9/P10
defect class), and it puts distortion-specific vocabulary into the element
layer permanently.

**The fix: derive the coefficients, do not declare them.** Every nonlinear
element already exposes a *pure* function `eval_i_pure(x, params, epar,
toolkit)`, and as of commit `148bd95` every toolkit can differentiate it —
exactly, via sympy, for the symbolic backends. So

```
b = (1/2)·d²i/dx²|₍x₀₎        c = (1/6)·d³i/dx³|₍x₀₎
```

comes from repeated application of an existing primitive. **No element
changes at all**, no new element API, and the coefficients cannot drift from
the model because they *are* the model.

The one honest addition this needs is a way to ask a toolkit for an
`n`-th derivative rather than only a Jacobian. Keep it minimal — a
`derivative_tensor(func, x, order, ...)` on `Toolkit`, defaulting to repeated
`jacobian`, overridden where a backend can do better. That is one method on
one class, not an API change rippling through the element layer.

### 2.4 Dependency discipline

Nothing here may make any optional backend a hard dependency. The
distortion analysis must import and run with only numpy and sympy present;
JAX stays optional (P1), and GiNaC/symengine stay optional. If a stage wants
a compiled backend for speed, it goes behind the same guarded-import pattern
`_jaxtoolkit` uses, with a test pinning that absence is survivable.

### 2.5 What stays out

The autonomous-oscillator branch of this literature is deliberately excluded
(reasoning doc §2). It shares the perturbation apparatus but requires solving
for amplitude *and* frequency self-consistently, which reintroduces the
nonlinear solve that §2.1's whole argument rests on avoiding.

## 3. The theory document, written alongside

`doc/src/circuit/distortion.rst`, built incrementally — the section covering
a stage lands in that stage's commit, not afterwards. Contents and owning
stage:

| Section | Written by |
|---|---|
| The perturbation setup; why `ε` is fictitious | Stage 1 |
| Where harmonic balance enters (it closes each order; it does *not* find the fundamental) | Stage 1 |
| Order-by-order recurrence, and why every step is a linear solve | Stage 1 |
| Relation to Volterra series, with the cubic coincidence shown live | Stage 1 |
| Cubic coefficients; the multi-nonlinearity matrix form | Stage 2 |
| Exponential nonlinearity via modified Bessel functions | Stage 3 |
| Two-tone bookkeeping: the `(m,n)` index pair | Stage 4 |
| Validity diagnostic and its limits | Stage 5 |

Numbers in that page are generated by `exec-rst` at build time (rule 1.2).
Note the existing hazard: a dead `exec-rst` block used to render its own
source at *info* level, so the build still said "succeeded" — fixed in
`3adfe69`, and `sphinx-build -W` makes it fatal. Verify by grepping the built
HTML for a number the block should have produced.

## 4. Stages and gates

### Prerequisite — element derivatives — **DONE** (`148bd95`)

`semiconductors.py` migrated onto `toolkit.jacobian`; `SymbolicToolkit`
gained exact sympy differentiation; `minimum`/`maximum` backend-parity gaps
closed.

*Gate:* symbolic device Jacobians agree with numeric finite differences.
**Passed** — 3e-10 (BJT, forward-active and saturation), 2.5e-07 (JFET,
saturation and triode). Suite 521 → 536.

### Stage 1 — cubic, single tone, single nonlinearity

Implement the 2013 paper's eqs. 18–20 over `G(s) = (A + sC)⁻¹` with scalar
`b`, `c` from §2.3. Start with `f_h = 0` (no input nonlinearity).

**Gate — the strongest one available, and first deliberately.** The result
must reproduce the **textbook Volterra series** HD2/HD3 for a cubic
nonlinearity, symbolically. The 2005 paper proves exact coincidence there
(its eqs. 43a–b, via Theorems 1 and 2). This checks the implementation
against a *different theory* rather than against the same authors'
arithmetic, and needs no paper-specific numbers at all.

**Decision that must be made in this stage, not later:** carry a *general
frequency index* from the outset. Two-tone (Stage 4) reuses the identical
recurrence with the scalar harmonic index `m` replaced by a pair `(m,n)`. If
Stage 1 hardcodes harmonics of a single `ω₀`, Stage 4 becomes a rewrite
rather than an extension.

### Stage 2 — several nonlinearities; numeric gate

Extend to the matrix form (`b`, `c` as K×K matrices), add the input
nonlinearity `f_h`.

**Gate.** Reproduce the 2013 paper's Fig. 4(a) for the three-stage RNMC
amplifier, using the **2009/2013 parameter set**: HD2 within **±1 dB**, HD3
within **±2 dB**. Those tolerances are set by graph-read uncertainty, not
optimism — there are no printed digits in any of these papers. Independently
reconstructed values to match: HD2 −31.8 dB @ 100 Hz, HD2 minimum −105.3 dB
@ 631 kHz, HD3 peak −65.5 dB @ 631 Hz.

Note the model restriction the paper carries (its eq. 6): cubic polynomials
with **self-terms only**, no cross-products `x_j x_k`. Record it in the
theory page rather than silently inheriting it.

### Stage 3 — exponential nonlinearity

Fourier coefficients via modified Bessel functions —
`F_{m,0} = 2·I_BQ·I_m(X̄₀)` and the derivative forms (2005 eqs. 46–47,
identical to 2009 eq. 20 under an obvious symbol mapping).

**Gate.** The common-emitter BJT of the 2005 paper's Example 1 — every value
printed (`I_BQ = 10 µA`, `V_T = 25 mV`, `β_F = 180`, `R = 250 Ω`,
`C = 37 pF`, `U = 50 mV`, `f₀ = 1 MHz`) — must match that paper's closed-form
eqs. 48a–b, which are *exact expressions* evaluable to arbitrary precision,
not graph reads.

**Why this stage is not optional.** The 2005 paper measures the cost of the
usual shortcut: approximating an exponential device by a cubic — standard
Volterra practice — converges to a materially different answer than true
harmonic balance, roughly **20% off on the second harmonic**. Cubic-only
would be quietly wrong for every BJT and diode circuit.

### Stage 4 — two tones; intermodulation

Generalise the frequency index to `(m,n)` for components at `mω₁ + nω₂`. The
second-order step becomes the 2013 paper's eq. 28 — a two-term enumeration
of the ways `2ω₁ − ω₂` decomposes, which is the 2005 paper's double
convolution evaluated and truncated to terms ∝ input³. Same mathematics,
different presentation; the two papers agree.

**Gate.** Reproduce the 2013 paper's Fig. 4(b)/(c): first-order term
−41.65 dB, second-order −40.92 dB, total −62.73 dB at low `f₁`, and
Fig. 4(c)'s relative-error plateaus (≈0.7 for the complete formula, ≈18.2 for
first-order-only). The second of those is the more useful check — it is a
*ratio*, so it is insensitive to absolute graph-reading error.

### Stage 5 — independent cross-check, and the validity diagnostic

Run pycircuit's own numeric `Transient` on a Stage-2 circuit, FFT the steady
state, compare harmonic amplitudes.

**Gate.** Agreement with the perturbation result inside the method's validity
range, and *visible divergence outside it* — the second half matters as much
as the first, because it demonstrates the diagnostic below actually
discriminates.

**This is the only check in the whole plan that does not depend on Buonomo
being right.** Every gate above verifies that we implemented their algebra
correctly; this one verifies the algebra predicts the circuit.

**Ship the diagnostic with it.** Report `‖M·f(x₀)‖ / ‖G·u‖` alongside every
result. **No paper in this set gives a quantitative bound on "weakly
nonlinear"** — 2009, 2012 and 2013 give none at all, and only 2005
acknowledges the convergence-radius problem exists. Without a reported ratio
the analysis returns a confident wrong number for a circuit outside its
range: a silent failure, and this project's standing preference is to fix the
mechanism that lets an error hide rather than only the error.

## 5. Outcomes

Filled in as stages complete. Negative results are recorded here as results,
not omitted.

| Stage | Outcome |
|---|---|
| Prerequisite — element derivatives | **Gate passed** (`148bd95`): symbolic Jacobians match finite differences to 3e-10 (BJT) / 2.5e-07 (JFET); suite 521 → 536 |
| 1 — cubic, single tone | **Gate passed** (`b8ec918`): HD2 and HD3 agree symbolically with Volterra kernels derived independently, on both a resistive and an RC circuit. Mutation-checked — dropping the second order fails 3 tests, a wrong Fourier factor fails 4, a wrong evaluation frequency fails 1. Suite 536 → 544. **Two findings:** `taylor_coefficients` silently fell back to default device parameters for non-`Semiconductor` elements (fixed; caught by comparing against a symbolic `IS` and getting a number back), and the general frequency index was carried from the start as planned, so `Harmonic((2,-1))` already resolves `2w1-w2` |
| 2 — several nonlinearities | not started |
| 3 — exponential | not started |
| 4 — two tones | not started |
| 5 — cross-check + diagnostic | not started |

## 6. Errata in the sources — implement around these

All found by numerically reconstructing the papers' own circuits, not by
reading them.

- **2013 paper eq. (49) is missing a factor `1/g₁`** — wrong by **+90.1 dB**
  uniformly. As printed it evaluates to A/V where HD3 must be dimensionless,
  and eq. (52) for the same circuit carries the factor. Use `eq.(49)/g₁`.
- **2012 paper eq. (9d) has transposed indices** — must be `M_{q,p}`, not
  `M_{p,q}`. Implemented literally it gives HD3 = +67.5 dB against a plotted
  −72 dB; corrected, −70.8 dB.
- **2012 paper has no eq. (9b)** — labels run 9a, 9c, 9d. The missing item is
  the first-harmonic correction, so the HD denominators use the *uncorrected*
  linear fundamental. Consistent with the 2005/2009/2013 papers, which all
  drop that correction at assembly, so it is a labelling gap and not a
  dropped term.
- Lesser: 2013 eq. (22)'s summation limit should be `k`, not `k+1`; its
  Fig. 6 caption says 10 mV but the plotted curves correspond to 20 mV
  (offsets of exactly +6 dB on HD2 and +12 dB on HD3 — the signature of a
  doubled input).
