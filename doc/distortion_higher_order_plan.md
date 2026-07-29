# Higher-order perturbation by Faà di Bruno — implementation plan

**Status: all four stages complete, all gates passed. Gates were declared in advance; outcomes in section 6.**

**Headline: raising the perturbation order helps, substantially — reversing the conclusion reached earlier from Picard iterations.**

Follows on from `doc/distortion_plan.md`, whose five stages are complete. This
one answers a question that work left genuinely open.

## 1. Why this exists

`harmonic_response` implements the published second-order truncation. Asked
whether raising the perturbation order would improve accuracy, the previous
work substituted **Picard iteration** and reported its behaviour as though it
described the perturbation series. It does not — a Picard iterate is the k-th
truncation *plus an unbalanced fragment of the orders beyond it* — and two
conclusions drawn from that substitution had to be retracted (commit
`0218a29`).

So the question stands, unanswered: **does a genuine third-, fourth- or
fifth-order perturbation truncation improve the circuit results?** Answering it
needs the terms built properly, which is what this plan is for.

**Do not expect a large win, and be willing to record a negative result.** The
scalar evidence is that the series converges monotonically but slowly — errors
6.6e-2, 2.2e-2, 9.1e-3, 4.3e-3, 2.2e-3 — while Picard, which is cheaper to
implement, converged faster on the same problem. It is entirely possible that
the honest outcome is "third order costs much more and buys little", and that
is a result worth having rather than a failure.

> **This prediction was wrong, and is left standing deliberately.** The
> measured outcome (§6) is a large win — roughly 100× at a drive well inside
> the validity range. The scalar problem was a poor guide because its
> convergence is governed by the contraction factor alone, whereas the circuit
> case also benefits from the harmonic structure being resolved more finely at
> higher truncation. Left unedited so the gates can be seen not to have been
> written around a known answer.

## 2. The mathematics, verified before planning around it

For `x = x_0 + ε x_1 + ε² x_2 + …`, the ε^k coefficient of `f(x)` is

```
z_k = sum_{m=1..k}  f^(m)(x_0)/m!  *  sum over ordered (j_1..j_m), sum j_i = k
                                          of  x_{j_1} · x_{j_2} · … · x_{j_m}
```

and the recurrence is `G x_k = M z_{k-1}` (2005 paper eq. 8e).

Checked symbolically against a direct series expansion for a cubic
nonlinearity, orders 1 through 5: **all agree**. (A first check reported
failure at k≥3; that was a broken `Derivative` substitution in the check
itself, not the formula. Recorded because the same trap will recur — verify
with a *concrete* `f`, not with symbolic derivative placeholders.)

## 3. The hard part, and it is not the formula

**Perturbation order is not the same as power of the drive**, and consistency
is required in the latter. This is the lesson `doc/src/circuit/distortion.rst`
records under "What the truncation is actually truncating": the third harmonic
at `U³` receives contributions from *both* first and second perturbation
order, and a truncation that keeps some `U⁴` terms without keeping all of them
can be worse than one that keeps none.

`x_n` carries terms of order `U^(n+1)` **and higher**. So truncating
consistently at `U^N` is not "compute `x_0 … x_{N-1}`" — it is "compute those,
then discard the parts of each that exceed `U^N`". Getting this wrong
reproduces exactly the defect that made the spectral experiment worse than the
published method.

**Proposed representation: grade every quantity by both harmonic index and
power of the drive.** A solution becomes `{harmonic m: {U-power p: coeff}}`,
and products convolve in both indices simultaneously. Truncation is then a
single filter on `p`, applied uniformly, and consistency is structural rather
than something the implementer has to remember.

The structure to expect for a single tone: harmonic `m` first appears at
`U^m`, then receives corrections at `U^(m+2)`, `U^(m+4)`, … — so a `U⁵`
truncation keeps harmonic 1 to `U⁵`, harmonic 2 to `U⁴`, harmonic 3 to `U⁵`,
and introduces harmonics 4 and 5.

## 4. What this forces to change outside itself

Reviewed through the same architect lens as the previous plan.

- **`Toolkit.nth_derivative` already exists** and takes an arbitrary order —
  built in stage 1 of the previous plan, exact under sympy, order-scaled
  central difference otherwise. `f^(m)(x_0)` for any `m` is available today.
- **Elements need no change.** Derivatives continue to come from
  `eval_i_pure`; nothing is declared.
- **`circuit.py` needs no change.** Still `G(jmω)` solves against the same
  matrix.
- **New, and local to `distortion.py`:** composition enumeration (pure
  combinatorics), and the graded arithmetic of §3.

So this is self-contained. That is worth stating because it is the main
argument for attempting it at all: the expensive prerequisite (arbitrary-order
derivatives through the toolkit) was already paid for.

## 5. Stages and gates

### Stage A — the scalar series, symbolically

Implement `z_k` and the recurrence for a scalar problem with no circuit and no
harmonics: `Y x + b x² = u`.

**Gate.** The terms must come out as the known closed forms —
`u/Y`, `−b u²/Y³`, `2b² u³/Y⁵`, `−5b³ u⁴/Y⁷`, `14b⁴ u⁵/Y⁹` — each a **single
monomial**, with coefficients 1, −1, 2, −5, 14 (the Catalan numbers, the
signature of a quadratic fixed point). Any term that expands to a sum means the
composition bookkeeping is wrong.

This gate is exact, symbolic, and needs no circuit, so it isolates the
combinatorics completely.

### Stage B — harmonics, truncated at `U³`

Add the harmonic index and the graded truncation of §3, then truncate at `U³`.

**Gate — a regression, and the strongest available.** The result must equal
`harmonic_response` **exactly** (to floating-point) on every circuit in
`test_distortion.py`. The published second-order form *is* the consistent `U³`
truncation, so a correct general implementation must reproduce it. If it does
not, the graded bookkeeping is wrong, and no higher-order result from it can be
trusted.

### Stage C — truncate at `U⁵`

Raise the truncation and nothing else.

**Gate.** Against the transient cross-check already built
(`test_distortion_vs_transient.py`), at drives *below* the convergence bound
`a < ln 2`: the `U⁵` result must be **at least as accurate** as the `U³` one at
every tested drive, and strictly better at the larger ones. Monotone
improvement is the claim being tested; if it fails, that is the answer and it
gets recorded as such.

Also required: symbolic expression size must stay bounded. The scalar terms
stay single monomials (Stage A); the circuit terms must not blow up
super-polynomially with order. **If they do, the method stops being suitable
for a symbolic tool regardless of accuracy** — which is the reason the whole
family was chosen (see the theory page, "Why this method suits a symbolic
tool").

### Stage D — is it worth it?

Measure cost against benefit and decide.

**Gate — a judgement, stated in advance so it is not rationalised
afterwards.** Third order is worth keeping if it either (a) improves accuracy
by more than 2× at drives where the second-order form is already within its
validity range, or (b) meaningfully widens the usable drive range before the
error exceeds some fixed target. If it does neither, **it gets recorded as a
negative result and not shipped as a default** — the second-order form costs
three linear solves and is within a fraction of a percent wherever it is valid
at all.

## 6. Outcomes

| Stage | Outcome |
|---|---|
| A — scalar series | **Gate passed exactly**, and one order beyond what was asked. Every term is a single monomial with Catalan coefficients: `u/Y`, `−b u²/Y³`, `2b² u³/Y⁵`, `−5b³ u⁴/Y⁷`, `14b⁴ u⁵/Y⁹`, and `−42b⁵ u⁶/Y¹¹`. Added `_compositions`, `composition_coefficient`, `perturbation_terms`. Five tests, including a check of the formula against a direct series expansion and an explicit guard that these terms are *not* the Picard iterates — the confusion this plan exists to correct. Composition counts verified as `2^(k-1)`, the expected total |
| B — `U³`, must reproduce `harmonic_response` | **Gate passed**, with one finding. Harmonics 2 and 3 reproduce `harmonic_response` to floating point (`0` to `2e-16`) at every drive tested. The **fundamental differs by design**: the published form drops its `U³` correction at assembly, the graded form carries it, and the difference is verified to scale as `U²` relative to the fundamental — the signature of an `O(U³)` term, not an error. Added `GradedSpectrum` (two-sided harmonics × drive-power grading) and `graded_response`. **Design note:** substituting-and-re-truncating in the graded ring *terminates exactly* on the consistent truncation rather than converging toward it, because each pass raises the lowest new power by at least one and the filter removes the rest — so the same filter that enforces consistency also removes the unbalanced tail that makes a plain Picard iterate a different object. **Process note:** the first test point had `\|X₂/X₁\| = 0.90`, essentially at the divergence threshold, where the fundamental "correction" came out 5× the fundamental — arithmetically right, operating point meaningless. A test now guards that the test point is weakly nonlinear |
| C — `U⁵` | **Gate passed decisively.** Against the transient cross-check, error falls monotonically with truncation order at every drive tested: at `a=0.221` HD3 goes 2.85% (`U³`) → 0.00% (`U⁵`); at `a=0.442` HD2 goes 8.40% → 0.70% → 0.05% (`U⁷`). Symbolic requirement also met — third-harmonic term count grows 2 → 7 → 16 for `U³` → `U⁵` → `U⁷`, polynomial in the order, against an iterative scheme's doubling |
| D — worth it? | **Yes, on both criteria.** (a) Better than 2× inside the validity range — HD3 improves ~100× at `a=0.221`, where the perturbation ratio is 0.05 and the method is comfortably valid. (b) The usable range widens materially: `U³` is marginal at `a=0.442` (8–12% error) where `U⁷` reaches 0.05%. **This reverses the earlier conclusion**, which was drawn from Picard iterations rather than perturbation order |

## 7. Things already known that bear on this

- **The convergence bound applies regardless.** Above `a = ln 2` the series
  diverges; higher order cannot help there and Stage C should not be tested
  there.
- **Device-model truncation is separately removable.** A quintic fit already
  matches the exact Bessel treatment to two decimals below the bound, so any
  residual error at that point is the perturbation truncation — which is
  exactly what this plan targets.
- **The harmonic cutoff is not the lever.** It saturates at six harmonics;
  raising it changes nothing.
