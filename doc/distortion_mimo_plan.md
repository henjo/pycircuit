# Multi-node graded perturbation response — implementation plan

**Status: planned, not started.**

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
| E — higher order, multi-node | **not run.** Depends on nothing from D; the open question it settles — whether term count stays polynomial in the order when there are several nodes — is the one that bears on the theory page's central claim |

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
  beyond what the paper prints.
