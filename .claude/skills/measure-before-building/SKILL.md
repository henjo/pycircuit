---
name: measure-before-building
description: Decide whether a performance change is worth building, and take a number you can trust when you do. Use before starting any optimisation, when a plan proposes work justified by speed, when a benchmark result decides something, and whenever a measurement disagrees with intuition.
---

# Measuring before building, and getting a number you can trust

An optimisation is a claim about where time goes. The claim is cheap to
test and expensive to assume. Measured on one campaign: **six proposed
optimisations were refused on measurement** — several of them
"obviously" worthwhile, one estimated at weeks — and every item that was
built was smaller than the item that was planned.

The refusals, with the number that killed each:

| proposed | measured |
|---|---|
| batch the flagship model | the two fast paths are **disjoint**: 22 classes can be compiled, 15 can be batched, **0 both** |
| a second compiler backend for the other path | element arithmetic is 18.2% of the run, 43.3% at scale — **ceiling 1.76x**, against 212x for the backend that exists |
| elide every block a card switches off | 15% of the generated code, exponentials included, removed: **1.00x** |
| a cheaper competitor to the shipped limiter | order-independent *and worse* — 909 → 1878 iterations, 8 new failures |
| memoise a hot failed lookup | 1.4–3.2%, and a recorded decision already argued against it |
| hold every sub-expression automatically | 52% of the holds a real model writes are ≤3 ops |

None of these was lazy thinking. Each had a plausible mechanism, and the
mechanism was usually real — it just wasn't where the time was.

## Before building: ask what it could possibly buy

**Compute the ceiling first.** Take the phase you mean to improve,
assume you make it *infinitely fast*, and see what the run does. If the
ceiling is not worth the work, stop — you have your answer in minutes
and you never needed the design.

Profile **by phase**, not by intuition. Wrap each stage, count and time
it, and check the shares sum to roughly the whole. A phase you did not
think to measure is where the time usually is.

⚠ **Op counts, line counts and call counts are the wrong proxy.** Here a
5-operation and a 130-operation conditional both measured ~7 µs; a
change that removed **15% of a Jacobian's source, exponentials
included,** measured **1.00x**; and dropping a hold made the emitted
source *shorter* while making the compile *not finish*. Measure time, or
measure the specific operation you believe dominates — never size.

## Taking a number you can trust

**Interleave the variants and repeat.** Never A-then-B once. Report the
minimum of several, or the median of pairs; state which.

⚠ **On a busy machine, build both variants in ONE process and interleave
their timings, then report the per-pair ratio.** Sequential runs of
*identical code* gave anything from **1.04x to 1.26x** here; paired, the
same comparison gave a per-pair median of **1.216x** with p10 1.168 and
p90 1.284. Each pair then sees the same load, and the spread tells you
how much to trust it.

⚠ **A profiler inflates what it measures.** Per-call costs read off a
profile put a set of primitives at ~58% of a Jacobian; timed standalone
they were **28%**. Use the profile to find *where*, then time the
suspect directly to learn *how much*.

⚠ **A speedup figure is a property of the whole call, not of the
change.** One memoisation measured 1.16x when it landed and **1.24x an
hour later, unchanged**, because a later fix removed a 425 ns lookup
from the same method and shrank the denominator. Where a figure
attributes a share to one component, write "re-run this" beside it
rather than a number to quote, and prefer publishing the **absolute**
before/after, which does not move when its neighbours do.

**Re-run every number at the commit that publishes it.** One campaign
published `both = 30` when the value was 13 *and had been 13 at that very
commit* — measured mid-fix, three changes later, never re-run. The test
passed, because the wrong number lived only in a docstring.

## On a machine you do not control, sample more -- do not abstain

A busy machine is a reason to take more samples, not a reason to give
up, and the statistic you gate on decides which of those you get.

Measured under four concurrent unrelated jobs, load 17 of 24 cores.
Minimum-of-N estimates of the same ratio, three independent trials:

| samples/side | trial 0 | trial 1 | trial 2 | spread |
|---|---|---|---|---|
| 3 | 1.008 | 0.916 | 1.054 | **14%** |
| 12 | 1.018 | 1.007 | 1.003 | **1.5%** |

At 3 samples the spread is the entire headroom of the bound being
asserted; at 12 it is 1.5% and centred on the published value. Noise
only ever ADDS time, so the minimum converges from above -- the only
question is how many samples it costs to find the floor, and under load
the floor is found late.

⚠ **Gate on whether the ESTIMATE reproduces, not on whether the MACHINE
held still.** A guard here measured the spread of raw timings inside one
measurement and abandoned the run above 1.25x. That statistic reads
1.86x-3.35x in exactly the trials whose estimates agreed to 1.5% -- it
measures the machine, and the question is the number. Gating on it would
have abstained on every good measurement. Take two independent estimates
and require THEM to agree.

⚠ **Retry before abstaining.** A machine busy on average is still quiet
in patches. Bounded retries (four, here) turned an abstain-half-the-time
guard into three clean runs out of three, and the cost is paid only when
the first pair disagrees -- a quiet machine still pays for exactly two.

⚠ **CPU time is NOT a load-immune instrument -- measured, not assumed.**
The obvious fix is to time `process_time` instead of the wall clock, on
the theory that contention shows up as descheduling. Recorded
simultaneously on identical work under that load: drift
`2.371 / 1.149 / 3.353` wall against `2.371 / 1.149 / 3.353` CPU.
Identical to three decimals. Contention there did not deschedule the
process; it slowed it. Check before switching instruments.

⚠ **"A coarse effect needs fewer samples" conflates two things.** A
bite-check resolving a 1.40x injected slowdown was given 3 samples on
that reasoning and abstained immediately, at 1.147x disagreement against
a 1.05x limit. A reproducibility gate constrains the PRECISION OF THE
ESTIMATE, which is +/-14% at 3 samples however large the effect being
measured. A coarse effect needs fewer samples to SEPARATE; it needs
exactly as many to REPRODUCE.

⚠ **Calibrate against the adverse condition when that is the condition.**
Waiting for an idle machine is only available if one is coming. Where
the box carries other people's work for hours, thresholds calibrated on
a quiet machine describe a situation that never occurs -- and a guard
that works only when quiet does not guard.

## Reading the result

⚠ **A ratio that did not move is not a null result.** When the published
figure compares two things that share a code path, a change to the
shared part speeds up *both* and the ratio sits still. A 10%
simulator-wide gain showed here as an unchanged 1.15x. Check the
absolute numbers on both sides before concluding nothing happened.

⚠ **A magnitude that matches a known mechanism is not that mechanism.**
An element method had 959 ns unaccounted for; a *missing* attribute
lookup in that tree costs ~1024 ns; the method visibly contained such a
lookup. Three facts agreeing. Measured before acting: that lookup
**resolves normally, in 39 ns**, and the 959 ns was two entirely
different costs.

⚠ **A confounded comparison reads as a result.** "14.2 s with the new
elements against 13.5 s with the old" was taken to mean device
evaluation was ~5% of the run. Both sides pay the same assembly cost, so
the comparison measures the gap between two implementations, not the
size of the phase they sit in. Profiled by phase: **59%**.

⚠ **When a null is surprising, find the mechanism before believing it.**
Removing 15% of a Jacobian's source bought nothing — because the removed
lines contained **zero** calls to the primitives the cost was made of
(executed primitive counts were *identical to the last call*). That
mechanism turned a disappointing null into the campaign's most useful
finding: a compact model's cost is its **regularisers**, not its physics.

## Refusing well

A refusal is a deliverable. Record:

* **the number that killed it**, so nobody re-derives the case;
* **what would reopen it** — the measurement that would change the
  answer;
* **the half that survives**, if any. "Elide inert blocks" was refused on
  speed while its *correctness* argument stood: a zero parameter does not
  switch a block off, and `0 x inf = NaN`.

And say plainly when a refusal is a judgement rather than a measurement —
a shipped guarantee, a recorded decision someone else made, a cost the
user should weigh. Those belong to whoever owns the trade, with the
numbers in hand.

## When the change ships

Gate it on measurements chosen *before* the work, and name them. Then
expect the gates to be incomplete: a set of four gates here — vendor
agreement, finiteness, library equivalence, cost — all passed while the
change quietly converted eager models to chained ones, forfeiting a 20.7x
batching path. It was caught by a **test model**, not by the library,
because the library's batchable classes never exercised the changed code.

**Ask what the gates cannot see.** Then go looking there.
