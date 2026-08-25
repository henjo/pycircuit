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

### And the cheapest failure of all: the number nobody re-ran

The same table then recorded `both` = 30 against `fet` = 12 as the
open problem a device-level limiter would fix. **`both` was 13**, at
the very commit that wrote 30 — the 30 was the row above, the
before-the-fix column, copied down into the after-table without being
re-measured. A one-iteration difference read as 2.5x for a day and was
used to justify a feature.

**Re-run every number in a before/after table, including the ones you
expect not to have moved.** A limiter's numbers move in whole
directions when the write-back changes, and the row you did not
re-measure is the row you will quote. Where the ratio matters, measure
it over a GRID of operating points and sum: on one point `both` costs
one iteration more than `fet`; over 48 points `both` is the cheapest
variant there is, and the single-point conclusion is simply false.

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

**But that sequence is not why you want one.** Built and measured: the
device-level form can reproduce SPICE's ordering, and SPICE's ordering
is *worse* — 1029 Jacobian evaluations against 927 for an
order-independent reading, over 48 operating points. The real reason
is the write-back, and it shows up as a **refusal, not a count**: a
per-probe write-back gives each probe a terminal of its own, so
SPICE's own four-probe MOSFET declaration — `fetlim(vgs)`,
`limvds(vds)`, `pnjlim(vbs)`, `pnjlim(vbd)`, four probes over four
terminals — has no terminal left for the fourth probe and *cannot be
compiled at all*. Count how many probes a device wants before
promising anything about how well two of them get along.

### If you build one: read the probes in SEQUENCE, not together

The natural device-level implementation — limit every probe from the
same unlimited vector, then solve for node values satisfying all of
them — **over-corrects**, measurably: 1040 against 927 on the same
grid, and the extra iterations come from writes onto source-pinned
nodes.

The reason is worth keeping. Probes share terminals, so one write
often satisfies the next probe outright. Measured: `vgs = 57.6 V` and
`vds = 59.6 V` off a pinned gate and drain with only the middle node
wild. Read in sequence, clamping `vds` brings that node back and `vgs`
no longer bites at all, so nothing but the wild node moves. Read
together, the two clamps disagree about the shared node by 2.5 V, and
the only way to honour both is to move the **gate** — the exact
failure the runtime write-back was built to remove.

So: read in a canonical sequence (largest correction first, ties by row
index — data-derived, so declaration order still cannot matter), each
probe seeing what the earlier ones did; then write the whole device
back at once. The write-back that is worth having is a spanning forest
over the terminals — keep every constraint that does not close a cycle,
drop the smallest one that does, hold the least-drifted node in each
component and derive the rest. That is the per-probe rule ("move the
terminal that drifted further") stated for a tree.

A limiter parameter that depends on the bias -- SPICE's `von`, the
turn-on recomputed each iteration through the body effect; a
self-heating device's saturation current at its *actual* temperature --
is evaluated at the **last accepted iterate**. That is the only
well-defined choice (the proposed iterate is the thing the limiter
exists to distrust) and it is exactly what SPICE does. A rule here used
to refuse such parameters outright -- "the limiter runs BEFORE the
device, so its parameters cannot read the solution" -- which was right
about the order and wrong about the conclusion, and cost one model
565 mV of looseness and another a duplicated saturation current before
it was measured. If a limiter parameter *cannot* follow the bias, say
so; do not accept "looser, not wrong" without pricing it.

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


## A device-level limiter is a capability, not a speed-up

Three measurements, three devices, one direction:

    two-probe FET star, 48-point grid    ties at 34, +1 at the rest
    MOS level 1, 48-point grid           552 vs 554 -- nothing
    Gummel-Poon, six cold starts         never better, once 13x worse
                                         (92 vs 7 Jacobians)

What grouping buys is *expressibility*: four probes over four terminals
cannot compile per-probe at all, because the fourth finds both its
terminals claimed. That is a real gap and it is the reason the
mechanism exists. But **a two-probe device should not declare it**, and
"the grouped form is more principled" is not a measurement. When an
intervention's justification is a capability, do not expect -- or
claim -- an iteration count.
