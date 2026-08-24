---
name: validation-design
description: Design measurements and tests that can actually catch an error — choosing the metric, the window, the tolerance and the coverage when validating an implementation against reference data. Use when writing or reviewing comparison tests, choosing tolerances, deciding a term is negligible, or when a suite is green and the model is still wrong.
---

# Designing a comparison that can fail

A test that cannot fail is not evidence. Most of what follows is about
the ways a comparison quietly loses its power to catch things.

## Choosing the metric

**Spread, not median.** The median ratio across a sweep says whether the
GAIN is right. The spread says whether the BIAS DEPENDENCE is — and that
is what the physics governs and what a missing term actually breaks. Two
terms were switched off here for making the median worse; judged by
spread they were worth 1.80 → 0.38 and the decision was reversed.

**A flat ratio is a well-posed question; a sweeping ratio is a broken
drive.** One number with one cause beats a curve.

**Convert to the unit the error lives in.** In an exponential regime a
current ratio is a voltage: `dV = ln(ours/ref) / (d ln I / dx)`. That
tells you the size in millivolts *and*, from whether it is flat, whether
it is an offset at all.

**Separate observables that respond differently.** Gain terms scale a
curve; electrostatic terms change its slope. Measuring slope and level
separately can distinguish mechanisms that no fit to a single curve can.

## Choosing the tolerance

**One tolerance across all biases of a device, not one per bias.** A bug
that is exact in one regime and wrong in another is *accommodated* by a
per-bias tolerance and *caught* by a shared one. A bug here was exact at
`Vds = 0.05` and 73% high at 1.2; only a shared tolerance sees that.

**Assert the positive; scan for boundaries.** A test asserting "this
breaks above 1e28" fails the day someone stops it breaking. Scan for how
far it gets and assert it has not got worse.

**Record known gaps as measured bands, not prose.** "Up to 1.73 in
saturation, 0.2% in the linear region" is a regression guard. "We do not
model X well" is not.

## Choosing the coverage

**Find the floors and masks.** A comparison helper with `FLOOR = 1e-6`
hid five decades of every transfer curve here — precisely the regime
where a threshold error is visible. Ask what the harness is *excluding*
before trusting what it reports.

**The reference has floors too.** Leakage, noise, terms you do not
model. Measuring below one measures *its* floor. Find it, stay above it,
and write the window down with the reason.

**A zero coefficient does not test a scaling.** Three terms were found
wrong or missing here because the one card they were tested against
happened to zero them. **Two cards is not twice the validation — it is
the difference between exercising a term and assuming it.**

**A bias no sweep visits tests nothing.** A body-bias correction that is
identically 1 at `Vsb = 0` is untested by every sweep at `Vsb = 0`.

**Test both signs, both channel types, both geometries.** Asymmetry is
diagnostic: when one side is right and the other is not, everything
shared is eliminated at a stroke.

**Defaults are not the card.** Defaults leave optional blocks off and
every quantity smaller. If a regime is only reachable with real
parameters, test it with real parameters — as its own test class, so it
cannot be assumed from the other one.

## The failure that survives all of the above

**A model checked only against itself is checked for CONSISTENCY, not
CORRECTNESS.** These are different claims.

A charge model here passed exact conservation, exact source/drain
symmetry, and an unfitted Ward–Dutton partition — and was 24% wrong.
Every one of those checks is a RATIO, and a uniform error in the oxide
capacitance divides out of all of them. No construction property could
have seen it.

The fix is an external quantity. Prefer, in order:

1. a quantity the reference **exports** (`lp_*` scaled parameters,
   operating-point outputs) compared term by term — a current can be
   right for compensating reasons, a parameter cannot;
2. a quantity the reference exports *implicitly* — two outputs plus one
   shared quantity can determine a third the model never prints;
3. a high-precision independent computation of the same equation.

## Golden data

**Regeneration must be additive.** After adding outputs, verify every
existing number reproduces bit-identically before trusting the new ones.

**Record the gap, not just the agreement.** Parameters you do not yet
model belong in the reference file too — a known gap written down is
better than a known gap remembered.

## Know where the floor is

**The reference's own precision is the floor**, and agreeing to it is
success, not a remaining gap. Recorded data carries maybe six
significant figures; a ratio of 1.00003 against it is the last digit,
not a 0.003% error to chase.

**Distinguish a model residual from a harness one before chasing it.**
The last 0.04% here was a **0.15 K** difference between the temperature
the parameters were scaled at and the one the model was evaluated at —
`sid ∝ φt`, so 300.0/300.15 = 0.9995 exactly. A residual that is
IDENTICAL across every geometry and both device types is a constant, and
a constant is far more likely to be a condition mismatch than physics.
Check temperature, bias convention, units and instance parameters first.

## Before publishing a number

- Did I measure in a clean window, or in the reference's noise?
- Is my extraction right? (Reading a `CY` entry at 1 kHz as white noise
  when it is four decades of flicker cost a day here.)
- Does anything else I have measured contradict it? **Two measurements
  that cannot both be true are a finding, not a nuisance** — do not
  reconcile them by inventing a mechanism that makes the contradiction
  sound like physics.
- **A tidy, interpretable answer is not thereby a correct one.** The
  most dangerous result is the one that looks like a clean explanation.
