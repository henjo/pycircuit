# Multi-node graded perturbation response — implementation plan

**Status: all five stages complete, all gates passed. Gates were declared in
advance; outcomes in section 6.**

**Headline: the multi-node extension works, and going multi-node costs a small
constant factor symbolically — not the blow-up an early mismeasurement
suggested.** Validated against four published closed forms across two
circuits, and against an independent ODE integration for the higher orders the
papers do not print.

Follows `doc/distortion_plan.md` (5 stages, complete) and
`doc/distortion_higher_order_plan.md` (4 stages, complete).

## 1. Why this exists

`graded_response` is **single-node scalar**: it takes `response(s) -> complex`
and one polynomial nonlinearity at one port. Every worked example in the
primary implementation reference — Buonomo & Lo Schiavo, AICSP 2013 — is
**multi-node with several nonlinearities**:

| Example | Nodes | Nonlinearities |
|---|---|---|
| 3-stage RNMC amplifier (eq. 45) | 3 | quadratic + cubic at two nodes |
| gm-C Tow–Thomas biquad (eq. 46) | 2 | four cubics, `g_1c g_2c g_3c g_4c` |

So the feature currently cannot run any published example end to end. That is
the gap this plan closes.

**What is already known to work.** Eliminating node 2 of the biquad by hand
reduces the `g_2c` and `g_3c` nonlinearities *exactly* to the scalar form, and
`graded_response` then reproduces the paper's published eq. (48) to **3e-16**
across five frequencies. So the mathematics and the grading are right; only the
plumbing is scalar. That measurement is also the strongest available evidence
that the sign convention below is correct.

**Sign convention, established by measurement, not assumed:**

```
(A + sC) x = G_h x_in + f_h(x_in) - f(x)
x = M(s) [ G_h x_in + f_h(x_in) - f(x) ],   M(s) = (A + sC)^-1
```

This is the convention under which the scalar reduction reproduced eq. (47)
exactly and eq. (48) to floating point. The opposite sign does not.

## 2. What this is *not* for

**Not** an attempt to obtain a 1 dB compression point. That question prompted
this work and the answer is already recorded: neither published example
compresses by 1 dB anywhere the series converges. The biquad reaches 0.0086 dB
at `X_in = 0.3 V` and its truncations disagree by 0.5 V; 1 dB would need about
3.2 V, roughly 8x beyond the convergence limit. The RNMC amplifier's cubic
coefficients are positive (`g_m2c = +60 mA/V^3`, `g_m3c = +3 mA/V^3`) so it
expands. **Reconsider if** a reference circuit turns up whose compression is
in-phase and whose output is a transconductor current rather than a node
voltage — 1 dB of cubic compression corresponds to a contraction factor of only
~0.44, which is well inside the bound, so the obstacle is the available
circuits, not the method.

**Not** general cross-term nonlinearities in the first instance — but see §4;
the chosen representation gets them for free, so this is a non-restriction.

**Not** two-tone. The existing `intermodulation_response` stays scalar for now.
**Reconsider if** stage C lands cheaply; the graded ring already carries tuple
harmonics, so the extension is mechanical.

## 3. The design

`GradedSpectrum` is unchanged — it stays the scalar per-node object.

New `GradedVector`: one `GradedSpectrum` per node. `__add__`, `scaled` and
`truncated` are componentwise. The only interesting operation is `through`:

```
for each (harmonic m, power p) present across any component:
    gather the length-n vector of coefficients at that key
    solve  M(j m w) . v          <- one linear solve per (m, p)
    scatter back
```

This is exactly where the scalar version multiplied by `response(1j*m*w)`. The
structural claim of the whole method survives intact: **every step is a linear
solve against the same matrix pencil, evaluated at a harmonic frequency.**

The nonlinearity is a **callable** `f(x_vec) -> vec`, written in the graded
ring using the existing `+`, `*` and `scaled`. That is fully general — cross
terms like `x1*x2` cost nothing extra — and needs no new coefficient
representation. The caller truncates; the driver re-truncates each pass.

## 4. What this forces to change outside itself

Reviewed through the architect lens, as the previous two plans were.

- **`GradedSpectrum` unchanged.** New code is additive.
- **`circuit.py` unchanged.** Still a solve against `(A + sC)`.
- **Elements unchanged.** Nothing is declared; coefficients still come from
  `eval_i_pure` where a real circuit is used.
- **`graded_response` unchanged and kept.** The scalar path is not a special
  case to be deleted — it is the regression oracle for stage A.
- **New and local to `distortion.py`:** `GradedVector`,
  `graded_response_mimo`.

The expensive prerequisites — the grading, the consistency filter, the
composition bookkeeping — were paid for by the previous two plans. This is
plumbing on top of them.

## 5. Stages and gates

Every gate below is a **published closed form or an existing verified result**.
None is a self-comparison, and none is fitted.

### Stage A — `GradedVector`, and scalar equivalence

**Gate.** A 1x1 system driven through `graded_response_mimo` must equal
`graded_response` to floating point on every circuit in `test_distortion.py`.
A block-diagonal 2x2 carrying two independent scalar problems must reproduce
both simultaneously. If MIMO does not contain scalar exactly, the plumbing is
wrong and nothing downstream can be trusted.

### Stage B — the biquad, `g_2c` and `g_3c` only

**Gate.** Must match the hand-eliminated scalar reduction, and hence the
`g_3c`/`g_2c` fraction of eq. (48), to ~1e-15 — the accuracy already measured
for the scalar path. This isolates the matrix plumbing from the new
nonlinearities.

### Stage C — the full biquad, all four cubics

Adds `g_4c` (acting on node 2, unreachable scalar-wise because `x2 = L(s) x1`
with `L` frequency-dependent, so `(L x1)^3 != L^3 x1^3` across harmonics) and
`g_1c` (an input nonlinearity, a third-harmonic drive term).

**Gate — the prize.** Must match the **complete** published eq. (48), both
fractions, to floating point. This is a closed form the current code cannot
reach at all, so passing it is new capability and not a regression.

### Stage D — the RNMC amplifier, three nodes

A second, structurally different circuit: three nodes, **quadratic as well as
cubic** terms, nonlinearities at two nodes, and — in the eq. (45) configuration
— external feedback.

**Gate.** Must match published eqs. (41) `HD2` and (42) `HD3` to floating
point. Note the recorded erratum for this paper: **eq. (49) is missing a factor
`1/g_1`** — that is the filter formula, not these, but it is a standing warning
that a printed form may be wrong. If a mismatch appears, dimensional analysis
comes before adjusting the code.

### Stage E — higher order on a multi-node circuit

Raise `max_power` to 5, 7 on the biquad.

**Gate.** The `U^3` truncation must equal the published closed form (the paper
*is* a `U^3` truncation), and higher orders must converge monotonically below
the bound, as they did in the scalar case. Symbolic term count must stay
polynomial in the order. **If a multi-node system makes the term count blow up
super-polynomially, that is a negative result worth recording** — it would mean
the method's suitability for a symbolic tool does not survive going to many
nodes, which is the single most important claim in the theory page.

## 6. Outcomes

_To be filled in as stages complete. Negative results recorded here, not
edited out._

| Stage | Outcome |
|---|---|
| A — `GradedVector`, scalar equivalence | **Gate passed.** 1x1 MIMO equals `graded_response` to **1.7e-15** worst case over 4 drives x 5 truncation orders x 4 harmonics; a block-diagonal 2x2 carries two independent problems to **2.2e-16**. Added `GradedVector` and `graded_response_mimo` |
| B — biquad, `g_2c`/`g_3c` only | **Gate passed**, worst **2.3e-14** against the `g_3c`/`g_2c` fraction of published eq. (48) at five frequencies |
| C — full biquad, all four cubics | **Gate passed**, worst **2.3e-14** against the *complete* published eq. (48). This is capability the scalar path does not have at all. **And the omitted terms were not small:** at 100 kHz the full third harmonic is `1.87e-14` against the scalar-reachable `1.37e-21` — a factor of 1.3e4. A scalar reduction is not a mild approximation off resonance, it misses essentially the whole answer. Guarded by its own test |
| D — RNMC amplifier, 3 nodes | **Gate passed**, worst **9.1e-16** on both `HD2` (eq. 41) and `HD3` (eq. 42) at six frequencies from 100 Hz to 10 MHz. Second circuit, structurally unlike the biquad: three nodes, **quadratic as well as cubic** terms, and a nonlinearity acting on a node other than the one it is injected into. The paper prints only the *feedback* matrix (eq. 45), so the open-loop matrix was reconstructed as its limit `k_in -> 1`, `k_3 -> 0`, `g_03e -> g_03`; that reconstruction has its own gate against published eq. (39), passed at **2.8e-16** |
| E — higher order, multi-node | **Gate passed, both halves.** *Symbolic size:* growth stays polynomial and multi-node costs a **small constant factor**, not an explosion — numerator terms over a common denominator go 1 → 10 → 46 (scalar) against 1 → 18 → 117 (two-node) for `U^3` → `U^5` → `U^7`, a ratio drifting 1.0 → 1.8 → 2.5. Graded keys are exactly 2x throughout, one set per node. *Convergence:* strictly monotone at every drive below the bound, against an independent ODE integration — at `Xin = 0.05` the error falls 1.9e-3 → 9.3e-5 → 2.9e-7 → 6.4e-9 → **2.3e-11** across `U^3`..`U^11` (600 cycles, `rtol=1e-13`; the doc table uses cheaper settings and bottoms out at 3.1e-10 instead — **that last cell measures the integrator, not the method**, and moves with integrator settings while no other cell does). It goes non-monotone only at `Xin = 0.3` (1.0e-1, 1.3e-1, 8.2e-3, 1.2e-2, 1.6e-3), which is the bound, and matches the asymptotic behaviour already documented for the scalar case |

**Two findings from stage E, both about how easily a wrong answer looked
convincing.**

*The multi-node "blow-up" was an artifact of the measurement.* The first run
showed `count_ops` going 52 → 766 → 13396 → 222924 against a scalar 6 → 17 →
31 → 45, which reads as a catastrophic explosion and would have refuted the
theory page's central claim. It was wrong twice over: the scalar baseline was
given `response(s) = 1/Y`, a bare symbol with **no `s` dependence**, so it
could never accumulate denominators; and raw `count_ops` counts redundancy
that putting the expression over a common denominator removes. Both errors
pushed the same way, which is what made the result look solid. With a fair
baseline the factor is ~2.5 and shrinking relative to order.

*Caveat kept deliberately:* `cancel` did not finish in 900 s at `U^7`, so the
counts are terms over a common denominator **without** full cancellation — an
upper bound on the true size. The practical limit for symbolic use may be
simplification *cost* rather than expression size, which is a different claim
from the one the theory page makes and is **not** established here.

**The pencil printed in eq. (46) has right-half-plane poles**
(`Re s = +1.68e6`). Every published result for this circuit is a modulus and
`|q|` is conjugation-invariant, so the sign is invisible in the paper's
figures — and all our stage B/C gates passed with it. It matters anyway: a
time-domain reference needs a stable system. **And the natural
over-generalisation is false** — identical `|q|` does *not* mean identical
magnitudes for the nonlinear result, which differs by 1-2%, because the graded
computation adds complex quantities across harmonics and addition is not
conjugation-invariant. Only products and ratios of `q` are, which is why the
published closed forms survive. Both facts are pinned by tests.

**A finding from stage D worth keeping: the published ratio uses the
*linearised* fundamental.** Dividing by the graded fundamental instead
disagrees with eq. (41) by 5.7e-3 at 100 Hz. That is not an error — the graded
form carries the fundamental's own `U^3` correction and the published form
drops it, exactly as recorded for stage B of the higher-order plan. What
distinguishes the two explanations is the *scaling*: the gap falls by 100x per
decade of drive (measured 99.4, 100.0), the signature of an `O(U^3)` term
against an `O(U)` fundamental. The test asserts the scaling rather than the
value, because only the scaling separates "by design" from "wrong".

The tell that pointed there: **HD2 and HD3 showed identical relative error at
every frequency.** Two different harmonics cannot be wrong by the same factor
unless the fault is in something they share — their denominator.

**Mutation check.** Flipping the sign of the nonlinear feedback term in
`graded_response_mimo` fails 12 of the 13 stage A–C tests, and 6 of the stage D
tests. Measured, not assumed: under that mutation `HD2` is unchanged to
7.96e-14% while `HD3` is wrong by **125%**. So `HD2` cannot detect a sign error
and `HD3` can — `HD3` sums a direct cubic against a second-order contribution
of the *quadratic* coefficient, and only their relative sign moves. This is why
both are asserted and not just the larger. The survivor is
`test_biquad_omitted_nonlinearities_are_not_negligible`, which compares a
ratio and is sign-insensitive by construction — expected, and noted so it is
not mistaken later for a weak test.

## 7. Things already known that bear on this

- **The biquad is a Duffing resonator.** At its centre frequency the cubic
  correction is ~90 degrees out of phase with the fundamental, so it detunes
  rather than compresses. Any amplitude-domain check should be done off
  resonance, or on the harmonics rather than the fundamental.
- **`q(jw) = g_3 g_4 + C_1 C_2 w^2 + g_2 C_2 jw`** with `g_3 g_4 < 0`; centre
  10.6931 MHz, `Q = 20.00`, unity gain at centre. Verified.
- **Paper 10 has four recorded errata**, one severe. Check dimensions before
  trusting a printed formula.
- **A transient cross-check needs a cubic transconductor element**, which
  pycircuit does not currently have. Stages B–D are gated on published closed
  forms instead, which is stronger; a transient would only be needed to go
  beyond what the paper prints. Stage E sidestepped it by integrating the
  circuit's own ODEs directly with `solve_ivp`, which needs no element.

## 8. What is left

Ordered by value, with the gate named where one already exists. Nothing here
is blocking anything shipped; all five stages are complete.

### 8.1 Two tones on the multi-node path — the recommended next step

`intermodulation_response` is still scalar-only, so the multi-node path is
single-tone. This is the natural completion of the work and unusually cheap:

- **The representation already supports it.** `Harmonic` has been a tuple from
  the outset and `GradedSpectrum` keys on `(harmonic, power)` with the harmonic
  index free to be a tuple, so two-tone is a matter of the driver and the
  convolution bookkeeping, not of the data structure.
- **Published gates already exist and are identified**: eq. (43) `IM3` for the
  amplifier and eq. (52) `IM3` for the filter. Both are in the same paper whose
  eqs. (41), (42), (47), (48) the single-tone path already matches to ~1e-15,
  so the gates are of known quality.
- **Note the recorded erratum before starting:** eq. (49) is missing a factor
  `1/g_1`, and eq. (52) — the `IM3` formula — is the one that is *right*. Check
  dimensions first.

**Reconsider if** the double convolution turns out to make term counts grow
faster than the single-tone case did; stage E's measurement should be repeated
for two tones rather than assumed to carry over.

### 8.2 Cross-terms `x_j x_k` are permitted but unverified

`graded_response_mimo` takes the nonlinearity as a callable in the graded ring,
where a cross term costs nothing extra — but **no published example in this
line uses one**, so the capability has never been checked against an external
reference. The docs say so explicitly rather than implying coverage.

**Reconsider if** a reference circuit with a genuine cross term turns up. Until
then the honest position is "the representation permits it", which is not the
same as "it works".

### 8.3 Simplification cost, as distinct from expression size

Stage E measured expression *size* and found it polynomial. It did **not**
measure simplification *cost*: `sympy.cancel` failed to finish in 900 s at
`U^7` on a two-node system. The theory page's claim is about size; the
practical limit for a user may well be cost. **These are different claims and
only one is established.**

Worth measuring properly, because if cost is the binding constraint the useful
fix is a different representation (rational functions kept factored, or
per-harmonic numeric evaluation) rather than anything about the perturbation
method.

### 8.4 A usable 1 dB compression example

The question that prompted the multi-node work, still unanswered, and **not
answerable from either published circuit** — see §2. The route that would work:
a minimal single-node transconductor cell using this paper's cubic model and
`alpha = -0.0535 V^-2`, with the output taken as the **transconductor current**
rather than a node voltage. 1 dB of cubic compression corresponds to a
contraction factor of about 0.44, comfortably inside the bound, so unlike the
published circuits it is reachable.

The topology would be ours rather than a reference's, so the result is a
demonstration of the method, **not** a validation of it. It should be labelled
that way wherever it lands.

### 8.5 A decision, not a task: the eq. (46) instability

The pencil as printed has right-half-plane poles. That is a defect in the
source, alongside the four errata already recorded for the same paper in
`distrortion_pertybation_reference.md`. Nothing in our code depends on
resolving it — the gates pass with the matrix as printed, and stage E
documents the flip. Whether it goes any further than our notes is a call for
the maintainer.

### 8.6 Open elsewhere, recorded here only as pointers

Not part of this plan; listed so this section is a complete answer to "what is
left" rather than a partial one.

- **Two tones on an *exponential* device** — §7 of `doc/distortion_plan.md`.
  Blocked on carrying complex amplitudes through the Bessel path, not on
  mathematics.
- **P15**, `doc/architecture.md` — `test_stress_stiff_rlc_pulse` is roughly
  two-thirds of total suite runtime. Three candidate routes written up, none
  attempted.
- **`supports('autodiff')` conflates two things** — exact differentiation and
  skipping Newton limiting — plus an unfixed `cosh` gap in `_symbolic` that
  leaves `VSwitch.G` unusable symbolically.
