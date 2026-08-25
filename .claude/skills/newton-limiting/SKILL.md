---
name: newton-limiting
description: Design, debug and validate convergence interventions in a circuit solver — $limit/pnjlim/fetlim/limvds write-backs, PCNR participation, gmin, and deciding whether an intervention is helping at all. Use when a device declares a limiter, when a solve diverges or is slow, and before trusting that a limiter works because its own tests pass.
---

# Limiters, and how to tell whether one is helping

A limiter is a *convergence intervention*: it changes the path Newton
takes and must never change the answer. That makes it unusually easy to
get wrong and unusually easy to believe in — it has no correct output
to check against, only a better or worse trajectory.

## The one test that matters, and the one everybody writes instead

**Measure the intervention against its ABSENCE, on a real solve, at the
scale it is used.** Not the limiter's own arithmetic — the solve.

Everything below passed on a limiter that was making convergence NINE
TIMES WORSE:

- the limited value matched SPICE's `DEVfetlim` exactly, over millions
  of points;
- a step of at most `2·VT` passed through untouched (the no-op-near-
  the-solution property);
- a wild step was compressed to SPICE's own number;
- it compiled and ran on both code generators.

All true. All irrelevant to the question. The defect was in the
*write-back* — where the limited value is put — which no test of the
limiting law can see.

So for every limiter, keep a solve that:

1. **fails, or is slow, without it** — otherwise you are measuring
   nothing;
2. **converges to the SAME point with and without** — a limiter that
   moves the answer is a bug, and this is the assertion that catches
   it;
3. **runs on a plain Newton, bypassing any rescue ladder.** A DC
   analysis with gmin stepping and source ramping will rescue the
   circuit either way and hide the entire effect. Here that meant
   calling `StandardNewton` directly.

Count Jacobian evaluations, not wall-clock.

## ⚠ A bad number with a plausible story attached is how this survives

The circuit-scale test above *existed*. It measured that `fetlim` alone
failed to converge, and recorded that as a **property**, with an
explanation: *SPICE's `fetlim` bounds a step to about a volt, and a
volt is forty thermal voltages to a subthreshold exponential.*

That is true, self-consistent, and not the reason. `fetlim` was failing
because of the write-back bug. Once fixed, `fetlim` alone became the
BEST of the variants — 12 iterations against `limvds`'s 57.

**A helper that measurably makes things worse is a defect report, not a
characteristic.** When a component does the opposite of its purpose,
suspect the component before you explain the measurement. An
explanation that arrives quickly and sounds like physics is exactly how
a wrong number gets closed.

## The write-back is where the bugs are

In a state-free convention the limiter returns a limited copy of the
solution vector, so something has to decide **which node moves** to
make the branch voltage equal `vlim`.

**Never write a bounded quantity as an absolute value derived from an
unbounded one.** `x[a] = x[b] + vlim` looks symmetric and is not:
`vlim` is bounded, `x[b]` is not, so a wild node hands its magnitude to
a sane one. Measured: a drain sitting at a perfectly good 40 V was
dragged to 6.1e8, then 6.1e9, then 5.8e10 — a decade per iteration —
because the source node was wild and elements are limited in instance
order, so the upper device read a node the lower one had not fixed yet.

**Which terminal to move is a RUNTIME question.** Deciding it at
compile time — "always the plus terminal" — cannot know which end is
pinned by a source, which is at a rail, or which one Newton is actually
being wild about. Move the terminal that has **drifted further from the
last accepted point**; the other one is information worth keeping.

**A limiter that did not bite must write NOTHING.** `a - (a - b)` is
not `b` in floating point, so even a no-op write destroys the
"unchanged near the solution" property — and that property is often
load-bearing, because *did limiting fire?* is used as a convergence
signal.

**Two probes may not move the same node**, or the second undoes the
first. When they compete for a shared terminal (a BJT's two junctions
on one base, a FET's gate and drain on one source), give it to the
probe applying the **larger correction**. Beware the tie: those two
junctions see *identical* drift when only the base has moved, so the
tie-break decides real behaviour rather than being cosmetic.

**Keep the result independent of declaration order.** If the choice
depends on the other probes, resolving it in list order makes the
answer depend on that order. Sort by something derived from the data
instead, and assert order-independence by literally reversing the list.

## What a per-probe limiter structurally cannot do

Worth knowing before promising it.

SPICE limits a MOSFET in a specific sequence: `fetlim(vgs)` first, then
`limvds(vds)` — but it recomputes `vds` from the **unlimited** `vgd`
first, so the coupling runs through a third branch that neither probe
names. **No per-probe limiter can reproduce that.** It needs a
device-level limiter receiving all of a device's voltages at once.

Likewise a limiter parameterised on a card value cannot follow a
bias-dependent one: `fetlim` wants SPICE's `von`, recomputed each
iteration through the body effect, and a parameter-only `vto` is
measurably loose — 565 mV at 2 V of body bias in one case. Looser, not
wrong: the step is still bounded. Say which it is.

## Limiting does not fix a singular row

A limiter bounds a *step*. A row that has gone empty is a different
failure and needs a `gmin` anchor or a DC path, not limiting.

But note the interaction: better limiting changes **which circuits
reach** the degenerate region. A stacked pair that reliably raised
`SingularMatrix` began converging after the write-back fix, because the
iterate stopped visiting deep cutoff where the channel conductance
underflows to exactly zero. The hazard was not removed, only made
harder to reach — so if you pin such a hazard in a test, re-check that
the test still reaches it after any convergence change.

## Testing a limiter without pinning its implementation

Four tests had to be rewritten when the write-back changed, and each
had mixed an invariant with an incidental. Separate them:

- **invariant** — each probe's branch ends carrying its own limited
  value; no probe is undone; nothing outside the element is touched;
  the result is finite and bounded; the solution is unchanged.
- **incidental** — *which* node moved, and in what order.

Assert the first. A test that pins the second will fail on every
improvement, and the failure tells you nothing.

One caution when asserting exactness: check the **write**, not the
recomputed difference. `(b + v) - b` is not `v`, and a rounding error
is not what the test is trying to catch. Handle the no-op case
explicitly rather than folding it into one of the write forms.
