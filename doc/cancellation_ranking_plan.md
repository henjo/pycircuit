# Cancellation-aware ranking — implementation plan

Reasoning behind this plan, including what is out of scope and what would bring
each rejected route back: `doc/cancellation_ranking_conclusions.md`. Read that
first; this document only says what will be done, in what order, and what
result makes each stage a success.

Predecessor: `doc/hierarchical_approximation_plan.md`. Its §5 records the
negative result this plan responds to.

**Status: stages 0-5 run. Stage 1's tolerance is met at both operating points
and its term-count budget at one; stage 3 meets its count and accuracy gates
and fails its interpretability gate; stage 5 recovers device parameters but
shows the composition is not error-controlled. See §"Outcomes" — every number
there comes from a logged run of the scripts named in it.**

**The next lever, named three times over and still not built:** a `group_rank`
that accepts a frontier of **several roots** rather than one. It is what stage 5
needs to rank a cofactor combination as one object, and it is the only thing
standing between the current result and an error-controlled hierarchical
approximation.

## Fixture, fixed in advance for every stage

`benchmark_circuits.ua741(symbolic_devices=('q1', 'q2', 'q17'))` — 26 unknowns,
flat diagram 1040 vertices, 2 773 885 product terms. Evaluated at **1 kHz**,
`s = 2πj·10^3`. Two operating points, as in the previous round:

| point | `gm_q1`, `gm_q2`, `gm_q17` |
|---|---|
| nominal | 1.0 mA/V |
| degraded | 0.1 mA/V |

**Corrected after stage 0.** This table first said 0.4 / 2 mA/V, which was taken
from the fixture's own defaults rather than from the previous round, and stage 0
was measured there before the mistake was caught (`κ = 6.66e+03` / `1.18e+05`).
Those points are *not* the ones the 994% baseline was measured at, so stage 0 was
re-run at the points above and only the re-run is reported. The reason it
mattered is not the size of the difference — it is that the baseline to beat
becomes unreproducible if the operating point drifts.

The baseline to beat is the measured failure of `DDD.approximate`: **994%**
relative error at nominal and **17 300%** at degraded, keeping 500 of 2 773 885
terms at `tol = 0.05`.

Everything is reported as `(term count, relative error)` **together**. A count
without an error is the specific mistake the previous round made and wrote down.

## Stage 0 — where does the cancellation live?

One pass over the DAG, no search, no new API. Compute the absolute-value
companion `A[v]` and the exact value `V[v]` for every vertex, and the ratio
`κ[v] = A[v] / |V[v]|`.

Report:

1. `κ[root]`, and the double-precision error it implies (`κ · 2^-53`).
2. The distribution of `κ[v]` by depth — the question of
   `cancellation_ranking_conclusions.md` §4, concentrated versus uniform.
3. The same two numbers recomputed in `mpmath` at the precision `κ[root]`
   demands, whenever `κ[root] · 2^-53 > 10^-6`.

**Gate, declared before running.**

- **PASS** if a cut of the diagram exists on which most groups have `κ[v] < 10`
  — i.e. cancellation is *concentrated*. Then stage 1 is worth building.
- **FAIL** if `κ[v]` stays within an order of magnitude of `κ[root]` at every
  depth. That refutes the premise, and the outcome is a recorded negative
  result about term ranking as such — not an implementation defect. Stop, write
  it down, and reconsider preconditioning (§6 of the reasoning document).
- Anything between is a **partial**: the route survives as a heuristic and
  stage 1 proceeds with its guarantee downgraded in advance.

## Stage 1 — group ranking as a script, not an API

Best-first expansion over groups, exactly as
`cancellation_ranking_conclusions.md` §2 defines them. A script under
`benchmarks/`, not a method on `DDD`, so that a failure costs nothing to throw
away.

- frontier starts as the single group `((), root)`;
- repeatedly split the largest-`|value|` non-terminal group;
- drop any child with `|value| < tol · |V[root]|`, accumulating the dropped sum
  exactly;
- stop when every non-terminal group is below threshold, or at a frontier cap.

**Simplified during implementation, and the simplification is the better
mechanism.** The drop threshold and the dropped-sum accumulator were both
dropped. Because `kept + frontier = det` holds exactly at every step,
`|det − kept| / |det|` is *already* the kept expression's own error, so the
stopping test needs no threshold and no bookkeeping: whatever is still in the
frontier when the error reaches `tol` is simply not kept. A threshold would have
been a tuning parameter standing in for a quantity that is directly available.

**Gate.** At `tol = 0.05`, at **both** operating points:

- relative error `≤ 0.05`, computed against `DDD.eval`, and
- at most **1000** kept groups, and
- the kept expression contains `gm_q1`, `gm_q2`, `gm_q17` or `s` — device
  symbols survive, which is the entire point;
- **FAIL** if the frontier exceeds 200 000 groups without reaching `tol`.

Beating 994% is not sufficient. The gate is `tol`, not "better than before".

## Stage 2 — the retained form, and whether it is readable

Only if stage 1 passes. Two things:

1. **Cancellation-stopped splitting.** Refuse to split a group whose split
   raises `κ` by more than a declared factor; keep it as
   `(Π P) · det(M[rows, cols])`. Requires recording node → `(rows, cols)`,
   which the builder's cache already knows and currently discards.
2. **The size distribution of retained minors.**

**Built differently, deliberately.** The retention rule shipped is *minor
dimension* (`max_minor=k`), not a `κ`-growth factor. Dimension is the quantity
that decides whether the retained determinant is readable, which is what the
gate is about, and it needs no threshold to tune. A `κ`-growth rule remains
worth trying — it would retain exactly the groups that cancel internally rather
than the ones that happen to be small — but it was not needed to answer the
stage's question. `minor_positions` was recovered from the diagram rather than
from the builder's cache, so no construction path changed.

**Gate.** The expression is symbolic in device parameters, and **the retained
minors are at most 6×6** (§3's reconsider-if). Report the distribution
regardless of pass or fail; the distribution is the useful output even when the
gate misses.

## Stage 3 — the leapfrog

Only if stage 2 says the form is readable. Group ranking on the 127-unknown
leapfrog via `HierarchicalDDD`, which is the only representation that reaches it.
(The predecessor plan records 299 vertices for this; that did not reproduce —
see the last section.)

**Gate.** Fewer than ~1000 groups, verified against the numeric solve to a
stated tolerance, **visibly different at a second operating point**, and — the
condition stage A failed — **naming device parameters rather than level
placeholders**.

## Stage 4 — library API, tests, docs

Only after stage 3. `DDD.approximate_groups` (name provisional), tests, and a
section in `doc/src/circuit/ddd.rst` with every number generated at build time.
The explanation ships in the same commit as the code.

## Stage 5 — device parameters back in the answer (declared 2026-07-30)

Stage 3's gate failed on one clause: the leapfrog's groups name level stamps, not
devices. Reasoning for the approach, including why the naive recursion through
111 single-node levels is hopeless: `cancellation_ranking_conclusions.md` §9.

**The change is to the partition, not the ranking.** Suppress each amplifier as
one block of 22 internal nodes — five blocks instead of 111 levels — so that
every block's cofactor family is over device entries and only the top matrix
carries stamps. Then group-rank three kinds of object separately: each block
determinant `D_l`, each stamp's cofactor combination, and the top diagram.

**Gate, declared before running.**

1. The five-block hierarchy builds in under 5 minutes, and its determinant
   agrees with `numpy.linalg.det` of the resolved reduced matrix to better than
   `1e-10` relative. Carried as a log magnitude — stage 3 established that the
   product underflows.
2. **No block's `A_ii` contains a `_lvl` symbol.** Asserted in code, not assumed:
   this is the structural claim the whole stage rests on, and if it is false the
   stage is over.
3. `κ` reported **per block** before any composition. If a block's `κ` is near 1
   there is nothing inside for group ranking to fix, and that is the stage's
   answer (§9's reconsider-if).
4. A composed approximation reaches relative error **≤ 1e-2** against the numeric
   oracle at **both** operating points, and the two are visibly different.
5. The composed expression **names device parameters** (`gm_s*_q*`) — the clause
   stage 3 failed.
6. The nested form's total operation count is reported and is **far below** the
   2 256 398 that plain unfolding costs. FAIL if it is not: an approximation that
   costs as much as the exact answer bought nothing.

**Declared in advance about the output shape:** it will be *nested*, not a flat
sum, because substituting stamps into top-level groups multiplies out. A flat
expression is not a goal of this stage and its absence is not a failure.

## Stage 6 — rank elements, not terms (declared 2026-07-30)

From §15d/§14c of the reasoning document: both literature sweeps converge on
*stop accounting error per term; evaluate the candidate answer*. Hu, Shi, Tai &
Lee (ISCAS 2015) rank circuit **elements** by the exact error their removal
causes, and their output matched textbook hand-derived formulas. Guerra et al.
(DATE 2000) score against a global numerical reference for the same reason.

`κ` cannot enter such a ranking: nothing sums term magnitudes, so the wall this
plan has been climbing is simply not in the path. And it is *simplification
before generation* — fewer symbols means a smaller diagram for everything
downstream.

**Adapted to an MNA matrix rather than a graph.** Their Short/Open pair are GPDD
edge operations. Here "Open" is setting a parameter to zero, which is immediate;
"Short" (parameter to infinity, i.e. node merging) is **not implemented** and is
stated as a limitation rather than skipped silently. So this is the weaker half of
their method, deliberately.

**The composition trap is designed out.** Individual removal errors must not be
added — §11 measured what that does. Selection is greedy with **exact joint
re-evaluation** after every accepted removal, which is the whole point.

**Gate, declared before running.** On `ua741` with all 24 transistor `gm`
symbolic, error measured as the **worst relative error of the transfer function
over a logarithmic frequency sweep** (sampled, and stated as sampled — the
literature is unanimous that a single frequency is not enough):

1. The ranking completes in under 5 minutes.
2. Report how many of the 24 symbols can be dropped at joint error ≤ 1% and
   ≤ 5%, with the joint error measured exactly at every step, never summed.
3. **The decisive one:** report `κ` of the simplified circuit's determinant
   against the full one. **PASS if `κ` falls by ≥ 10×** — that would mean element
   pruning helps the downstream term ranking, which is why it is being tried here
   and not merely as a circuit-reduction tool.
4. Report the surviving symbols **by name**, and report the diagram size and term
   count of the simplified circuit. A designer-readable answer is the deliverable.
5. Report whether the greedy order's joint error is monotone. Non-monotonicity is
   informative, not a failure — but it must be looked at, because it is the
   signature §11 found.

**Declared in advance:** it is entirely possible that `κ` barely moves while the
symbol count drops a lot. That would still be a useful circuit-reduction result
and a *failure* of gate 3, and it must be reported as both.

## Stage 7b — is determinant-side de-cancellation polynomial? (declared 2026-07-30)

§17 measured a path-predicate de-cancellation losing sharing at 1.7-1.8× per
circuit section. Tan & Shi's construction instead removes cancelling terms by
**set operations between sub-diagrams** (`REMAINDER(P, L_y)` inside a modified
coefficient-multiply), which forms no path-dependent state. Before building that
machinery, test *why* it should work, because the reason is checkable in one
function.

**The hypothesis.** The naive predicate carries the whole set of devices used so
far. Most of that is irrelevant: once a device has been used at position `(p,p)`,
its other stamp positions `(p,q)` and `(q,p)` need row or column `p`, which the
minor no longer has — so **exactly one partner position can remain reachable**, and
it stops mattering once the minor drops its row or column. So the state that
actually matters is not "devices used" but "**still-reachable forbidden
positions**", and that set should stay small because a minor only borders a few
devices.

**Gate, declared before running.** Re-count the memo states with the key
canonicalised to the still-reachable forbidden positions only, on RC ladders
`N = 3..10`:

1. **Correctness first:** the canonicalised count must yield the **same number of
   surviving terms** as the naive predicate. If it does not, the canonicalisation
   is wrong and nothing else in this stage means anything.
2. **PASS if the ratio to the plain (un-de-cancelled) state count stays below ~3×
   and is flat or near-flat in `N`**, against the naive key's 2.6→23.4×.
3. **FAIL if the ratio still grows geometrically.** That would mean
   determinant-side de-cancellation is inherently exponential, Tan & Shi's
   `REMAINDER` construction cannot rescue it either, and the topological route is
   the answer rather than one option among three.

Declared in advance: a pass justifies building the diagram machinery; a fail
**settles the architecture question** in favour of GPDD and saves that work.

## Stage 7g — the payoff test (declared 2026-07-30)

Everything up to 7f measured *properties* of the de-cancelled expansion. This
ranks it, which is the experiment the whole plan exists for.

Ranking runs over the de-cancelled state graph directly rather than over a
materialised `DDD`. That is a deliberate scope cut: a per-coefficient diagram would
be roughly 5 million vertices with first-row pivoting, and the question here is
whether the *algorithm* converges, not whether the object can be built. Stated as a
limitation, not skipped.

**Gate, declared before running**, at `tol = 0.05` on the µA741 at 1 kHz:

1. **The control is measured at the same operating point**, not quoted from
   earlier runs at a different one. 7f's near-miss came from a convenient control.
2. **PASS if group ranking on the de-cancelled expansion needs at least 5× fewer
   terms than group ranking on the compact diagram** at the same tolerance and
   point.
3. Report **magnitude** ranking on both as well. Tan & Shi's pipeline is
   magnitude-ranked; if de-cancellation makes plain magnitude ranking converge, that
   reproduces their result and closes the loop on §15b.
4. The kept sum must be verified against the exact de-cancelled value, which is
   itself already verified against `numpy.linalg.det`. A ranking that converges to
   the wrong number is not a result.

**Declared in advance:** the fixed-`s` de-cancelled `κ` is 99, not the 14-17 of the
low-order coefficients, so this tests the weaker of the two available settings. If
it passes here it will do better per coefficient; if it fails here that does not
condemn the per-coefficient route.

## Stage 8 — the concentration diagnostic (declared 2026-07-30)

§20's conclusion: `κ` measures *conditioning* and says nothing about term counts,
which depend on *concentration*. This builds the missing quantity.

**The cheap form.** Alongside `A[v] = Σ|term|`, the sum of squared magnitudes
`S2[v] = Σ|term|²` is computable in the same shaped single pass, because
`|Π|² = Π|·|²`. Their ratio is the standard participation ratio, an **effective
number of terms**:

    N_eff = A[root]² / S2[root]

`N_eff = N` when all terms are equal; `N_eff → 1` when one dominates. It is `O(|V|)`,
exactly like `cancellation`, so it is a genuine peer diagnostic rather than a
disguised enumeration.

**Gate, declared before running.**

1. **Sanity, against cases with known answers:** `N_eff = N` to within a few
   percent for a diagram whose terms are equal in magnitude, and `N_eff ≈ 1` for one
   with a single dominant term.
2. **Validation against the exact quantity** on circuits small enough to enumerate:
   `N_eff` must track the measured number of terms needed to reach 99% of the
   absolute mass — same order of magnitude, monotone in the same direction.
3. **The decisive one — it must predict what `κ` could not.** §20 measured compact
   ranking converging in 870 terms and de-cancelled ranking failing past 72 724.
   **PASS if `N_eff(compact) ≪ N_eff(de-cancelled)`**, i.e. if the diagnostic orders
   the two representations the way their ranking costs actually came out. `κ` ordered
   them backwards (99 against 6 659), so this is the discriminating test.
4. **FAIL if `N_eff` also ranks them backwards**, or if it disagrees with the
   enumerated 99%-mass count by more than an order of magnitude. Either would mean
   the participation ratio is the wrong concentration measure and the right one is
   still missing.

## Stage 9 — what reduces N_eff? (declared 2026-07-30)

Stage 8 made the readability question answerable for the first time. The compact
µA741 has `N_eff ≈ 194`, and group ranking duly needs 870 terms. **A readable
expression means tens of terms, so it means `N_eff` of order ten.** The question is
therefore no longer "how do we approximate better" but **"what transformation
lowers `N_eff` while keeping device symbols?"**

Three levers are available and none has been measured with this diagnostic:

1. **s-expansion.** §13 measured per-coefficient `κ` (worse than the whole
   determinant) and never per-coefficient `N_eff`. Ranking per coefficient is what
   Tan & Shi do, and a coefficient may be far more concentrated than the frequency-
   weighted whole.
2. **Element pruning** (stage 6). It left the compact diagram structurally
   identical and `κ` unchanged, so it was written off — but `N_eff` was never looked
   at, and zeroing devices changes the term *magnitudes* even when it changes no
   vertex.
3. **The `Short` operation**, stage 6's unfinished half, which collapses nodes and
   so genuinely shrinks the matrix.

**Gate, declared before running.** On the µA741 at 1 kHz, both operating points:

1. Report `N_eff`, `κ`, term count and contribution magnitude **per coefficient of
   `s`**, alongside the whole-determinant baseline of ~194.
2. **PASS if some coefficient that materially carries the response has
   `N_eff <= 20`** — an order of magnitude below the whole determinant — **and
   ranking it to 5% then needs at most ~30 terms with device symbols intact.** That
   would be a readable symbolic result and the goal of this entire plan.
3. **PARTIAL if `N_eff` drops usefully but not to ~10**: record the profile, since
   it says which lever to pull next.
4. **FAIL if per-coefficient `N_eff` is no better than the whole determinant's.**
   Then s-expansion is not the lever either, and levers 2 and 3 are what remain.

Declared in advance: the low-order coefficients dominate the response at 1 kHz, so a
small `N_eff` in a coefficient that contributes nothing is not a pass — the
contribution magnitude is reported next to it precisely so that cannot be claimed.

## Stage 10 — how few terms at a tolerance a designer would accept? (declared 2026-07-30)

Every measurement in this plan has used `tol = 0.05`, inherited from the first
round and never questioned. Stage 9 got the dominant coefficient to 155 terms
there. **But a designer reading a hand analysis routinely accepts 20-30%** — the
textbook expressions symbolic analysis is compared against are order-of-magnitude
arguments, not 5% models.

So before reaching for another transformation, measure the trade the existing
machinery already offers: **terms against tolerance, on the coefficient that carries
the response.** This is nearly free and nobody has looked at it.

**Gate, declared before running.** On the µA741's dominant coefficient, both
operating points:

1. Report the full curve — terms and achieved error at `tol` = 0.5, 0.3, 0.2, 0.1,
   0.05, 0.01 — with the device symbols present at each.
2. **PASS if some setting gives ≤ 30 terms at ≤ 20% error with device symbols
   intact**, and **the expression is printed** so the claim "readable" can be
   judged rather than asserted.
3. **FAIL if 20% still costs hundreds of terms**, which would mean the term count is
   insensitive to tolerance and only a change of representation can help.

**Declared in advance, because this gate is the easiest in the plan to fool:** a
small term count at a loose tolerance is only a result if the *expression* is
something a person can read. The printed expression is the evidence, and if it is
twenty products of six factors each then the gate is met on its letter and failed on
its purpose — which has happened once already in this project's history (stage A of
the predecessor plan) and must be called the same way if it happens again.

## Stage 11 — the transfer function, in operations, over a sweep (declared 2026-07-30)

Two things at once, because both are operation counts and stating either honestly
needs the other.

**Part 1 — re-read the existing verdicts with the corrected metric.** §23 found term
count overstates the answer 4-8×. Several "not readable" verdicts were reached on
term counts and never re-measured: the whole determinant's 870-term ranking, and
**the leapfrog's 181 groups** from stage 3. If those collect to a few hundred
operations, the verdicts change.

**Part 2 — the deliverable.** Everything so far approximates `det(A)`, the
*denominator*. A designer wants `H(s) = N(s)/D(s)`. Both are determinants, so
`s_expand` applies to each; the numerator is `A` with the output column replaced by
`b`.

**Two methodological requirements, from lessons already recorded in this plan.**

* **Verify over a frequency sweep, not at one point.** Every measurement in this
  plan bar stage 6 has been at 1 kHz, and the literature is unanimous that
  single-frequency error control guarantees nothing between samples
  (Rodríguez-García et al., DATE 1999). Stated as sampled either way.
* **Choose which coefficients to keep greedily, with exact re-evaluation of `H` over
  the sweep after each drop** — never by a per-coefficient budget. That is §14c's
  global-reference lesson and §11's non-monotonicity, and stage 6 showed exact
  re-evaluation removes it.

**Gate, declared before running.**

1. Part 1 reported: operations for the whole-determinant ranking and for the
   leapfrog's 181 groups, beside their term counts.
2. **PASS if `H(s)` is under 200 operations while holding ≤ 20% error across two
   decades**, with device symbols intact and visibly different at a second operating
   point.
3. **PARTIAL if it holds the error but costs more than 200 operations** — report the
   count, since it is the first honest measure of what a readable transfer function
   would cost.
4. **FAIL if no coefficient subset holds the error across the sweep**, which would
   mean stage 10's single-frequency result does not extend and the reconsider-if
   recorded there was right.

**Declared in advance:** stage 10 succeeded partly because `s^1` carries 97.3% of the
response *at 1 kHz*. Across two decades several coefficients will share it, so the
operation count should be expected to be several times stage 10's eleven, and a
result of 100-200 operations would still be a good outcome rather than a
disappointment.

## Stage 12 — does 177 operations survive a real band? (declared 2026-07-30)

Stage 11's reconsider-if, taken up immediately because it is the one claim in that
result most likely to be wrong. The sweep was 100 Hz-10 kHz, **below the µA741's
compensation pole and its unity-gain frequency**, so the coefficients that shape the
rolloff were never exercised. A symbolic transfer function that holds only below the
dominant pole is not a transfer function a designer would use.

**Gate, declared before running.** Same greedy selection against a global sweep
error, same two operating points, but over bands that actually contain the amplifier's
dynamics — up to 10 MHz, six decades:

1. Report `|H|` across the wide sweep, so the band being claimed is visible and the
   pole locations are evident rather than assumed.
2. Report kept coefficients, operations and sweep error for the narrow band
   (100 Hz-10 kHz) and the wide one side by side, at both operating points.
3. **PASS if the wide band still holds ≤ 20% under 200 operations.**
4. **PARTIAL if the error holds but the operation count grows** — report by how
   much, since that is the real price of covering the band.
5. **FAIL if no coefficient subset holds ≤ 20% across the wide band**, which would
   confine stage 11's result to a narrow band and should be said plainly.

**Declared in advance:** more coefficients matter over six decades than over two, so
the count is *expected* to rise. A result of 300-600 operations holding 20% over six
decades would be a good outcome and should be reported as one, not dressed up as a
pass.

## Stage 13 — the `Short` operation (declared 2026-07-30)

Stage 6's unfinished half, and the last cheap lever. `Open` (parameter to zero) left
the compact diagram structurally identical, because zeroing one addend inside a
composite entry removes no vertex. **`Short` collapses two nodes into one, so the
matrix loses a row and a column** — the only operation measured so far that makes the
determinant genuinely smaller.

Implemented as a node merge on the matrix: `V_a = V_b` means adding column `b` into
column `a` and row `b` into row `a`, then deleting `b`. Chosen greedily with **exact
re-evaluation of `H` over the wide sweep after each merge** — the discipline stage 12
had to relearn.

**Gate, declared before running.**

1. **Calibration first.** The merge must reproduce a hand-computable case exactly: a
   resistive divider whose upper leg is shorted has `H = 1`. A node merge is easy to
   get subtly wrong — double-counting the diagonal, or shifting the output index —
   and a wrong merge would silently produce a smaller matrix for the wrong circuit.
   **If the calibration fails, nothing else in this stage means anything.**
2. Report merges accepted at 20% sweep error over **10 Hz–10 MHz**, the resulting
   matrix dimension, and the reduced circuit's term count and `N_eff`.
3. **PASS if the reduced circuit's `H(s)` comes under 200 operations at ≤ 20% over
   the wide band** — i.e. if `Short` closes stage 12's 36× gap.
4. **PARTIAL if it cuts the operation count materially without reaching 200.**
5. **FAIL if no merge is acceptable**, which would say the µA741 has no two nodes
   that can be identified within 20% and the lever is empty on this circuit.

**Declared in advance:** an amplifier's nodes are mostly at genuinely different
potentials, so few merges should be expected. A result of two or three merges would
be unsurprising and would not by itself close a 36× gap.

## What would make this whole plan fail

- Stage 0 finds uniform cancellation. Most likely single outcome, and the
  cheapest to discover.
- The frontier grows without the error falling — exact at every step, but never
  compact.
- The retained minors are large, so the compact form is compact over objects
  nobody can read. This is stage A's failure recurring in a new shape, and it is
  the reason stage 2 has a size gate rather than only a symbol-presence gate.
- Double precision turns out insufficient everywhere, not just near the root, so
  group values are no longer cheap.

## Outcomes

_From measured runs only. Negative results recorded here, not edited out — the
DDD work kept five and they are among its more useful output._

Scripts: `benchmarks/cancellation_profile.py` (stage 0),
`benchmarks/cancellation_groups.py` (stages 1-2),
`benchmarks/cancellation_leapfrog.py` (stage 3).

| Stage | Outcome |
|---|---|
| 0 — cancellation profile | **PASS, and it explains the previous round quantitatively.** `κ[root] = 9.376e+03` nominal, `1.403e+05` degraded, so magnitude ranking had to capture **99.99947%** / **99.99996%** of the total absolute mass to reach `tol=0.05`. Cancellation is **not uniform**: median `κ` falls from 9.4e3 at the root to 1.0 by depth 27, monotonically. Double precision is ample (`κ·2^-53 = 1.6e-11`), so mpmath was not needed. |
| 1 — group ranking | **Tolerance met at both points; count budget met at one.** `tol=0.05` reached with **734 groups / err 4.11e-02** at nominal (0.2 s) and **2401 groups / err 3.70e-02** at degraded (1.2 s). The gate asked for ≤1000, so degraded **misses on count while converging on accuracy**. Also 337 / 185 groups at `tol=0.2` and 7573 / 8788 at `tol=1e-3`. Device symbols present throughout. |
| baseline check | **The previous round's numbers reproduce exactly.** `DDD.approximate` at 500 terms: err **9.938** (994%) nominal, **1.728e+02** (17 280%) degraded — the recorded figures, re-measured rather than quoted. At *equal* term count the comparison is 11.86 vs 0.041 (nominal, 734 terms) and 13.78 vs 0.037 (degraded, 2401 terms): **289× and 372× lower error for the same expression size.** |
| 2 — retained form readable | **PASS, and the trade is monotone.** At `tol=0.05`, retaining minors instead of expanding them: 2×2 → 516 / 1669 items, 3×3 → 353 / 1007, **6×6 → 182 / 503** (err 3.3e-02 / 8.5e-03). So the 6×6 form is **4.0× and 4.8× smaller** than the fully expanded one, names `gm_q1`, `gm_q2`, `gm_q17`, `s` throughout, and brings the degraded point **inside** stage 1's 1000-item budget that expansion missed. |
| 3 — the leapfrog | **Count and accuracy gates PASS; the interpretability gate FAILS, as predicted.** 127 unknowns, 16 symbols, 536 nonzeros; 111 levels, top diagram 16×16 / 1847 vertices / **374 608 terms**. `κ` of the *top* diagram is only **13.8** nominal and **160** degraded. `tol=1e-3` reached with **181 groups** (0.07 s) and **221** at degraded; `tol=1e-8` with 821 / 1248. Verified against `DDD.eval` (1.3e-15) and against `numpy.linalg.det` of the reduced matrix (4.1e-15). **But every symbol in every group is a level stamp (`_lvl110_*`); not one device parameter appears.** |
| 5 — device parameters back in | **Gate 5 PASSED, gate 4 FAILED, and the failure is the result.** Device parameters are back: the expressions name `gm_s0_q2`, `gm_s0_q17`, `gm_s1_q2`, … — identified transistors in identified stages, which is the clause stage 3 failed. Gate 3 confirmed the motivating diagnosis (`κ = 1.1e3` per block vs `13.8` at the top). **But the composed error is not controlled by the per-piece tolerances.** At degraded `gm` it never reaches 1e-2 and moves *non-monotonically*: 1.47e-2 → **1.13e-1** → 1.48e-2 as every tolerance is tightened 4× then 5×, while the operation count runs 1.6M → 3.2M → 5.5M. Cause: each reduced entry sums 25 heavily-cancelling cofactors, so 25 individually-small residuals combine unpredictably. **Hierarchical symbolic approximation is not compositional.** Verified throughout against `numpy.linalg.slogdet` of the full 127×127 matrix. Details: `cancellation_ranking_conclusions.md` §11. |
| 5 — first attempt | **Gate 2 FAILED, and the failure was narrow.** Sequential five-block suppression does *not* give device-level cofactors: `_suppress` renames **every** nonzero of the reduced matrix, so block 1's `A_ii` came back over `_lvl0_*` stamps. The topology claim it rested on is nonetheless **verified** — 0 entries couple one amplifier's internals to another's — so the blocks are genuinely independent and can be eliminated **in parallel against the original matrix** instead of in sequence. Also measured: the five-block hierarchy is **22 163 vertices / 1 076 448 top terms** against 1 958 / 374 608 for 111 single-node levels, so naming the blocks costs ~11× in representation. Retried with the parallel construction, which worked — see the row above and `cancellation_ranking_conclusions.md` §10. |
| 6 — rank elements, not terms | **Gates 1, 2, 5 pass; gate 3 FAILS and gate 4 fails informatively.** 5 of 24 `gm` symbols drop at 1% joint error, 7 at 5%, in 0.3-0.5 s, greedy order monotone, survivors named, bias devices leaving first exactly as a designer would expect. **But the diagram is bit-for-bit unchanged** — 1040 vertices, 2 773 885 terms — and `κ` moves 1.00× / **0.61×** (worse). Cause: a compact-symbol DDD has one vertex per matrix *entry*, and zeroing a device removes an *addend inside* an entry, so no vertex disappears. **Element pruning and the compact-symbol representation do not compose.** Useful as circuit reduction and design insight; not a route to a smaller symbolic expression. Details: `cancellation_ranking_conclusions.md` §16. |
| 7a — de-cancellation, calibrated | **Both gates PASS on the theory; the naive route is measured dead.** On Song & Shi's ladder (target re-derived here from the matrix *and* from the spanning trees) de-cancellation takes `κ` from **40.61 to exactly 1.000000** with the determinant preserved. But a de-cancelled expansion cannot key its memo on the minor alone, and the sharing loss grows geometrically: **2.6× → 4.3× → 7.3× → 12.8× → 23.4×** for RC ladders N=3..7, ~1.7-1.8× per section. Hopeless at 26×26. Meanwhile the *answer* shrinks (29%→72% of terms removed) and the survivors are the spanning-tree counts 5, 13, 34, 89, 233 — so **the cancellation-free answer is the spanning-tree set**, which strengthens the topological route this plan had put last. Next, if staying determinant-side: Tan & Shi's `CL`/set-operation construction during s-expansion, which forms no path-dependent state. Details: `cancellation_ranking_conclusions.md` §17. |
| 7b — is it polynomial? | **PASS, decisively, and §17's verdict is overturned for the correct construction.** Re-keying the memo on the *still-reachable forbidden labels* instead of the path gives a **flat ~1.2× overhead** (1.20× at N=3, 1.19× at N=10) where the naive predicate ran to **170.5×**. Surviving-term counts identical at every N (asserted). 86% of terms removed at N=10 — Tan & Shi's 70-90%, reproduced. **Gate 7b-1 caught a real bug:** the first version forbade matrix *positions*, but one position carries several devices' addends (`C2*s + 1/R1`), so it silently lost a term — the state must be a set of **labels**, i.e. (device, position) pairs. |
| 7c — does it move `κ`? | **No, and the reason is instructive.** Determinant preserved against `numpy.linalg.det` at every size (1e-15). But `κ` goes 3.08→1.48 at N=4 and **1.11→1.97 (worse) at N=9**. An RC ladder has `κ ≈ 1` already, so it is the wrong test; and the residue is **phase**, not sign — at fixed complex `s`, conductance terms are real and `C·s` terms imaginary. |
| 7d — with s-expansion | **PASS, exactly.** 7c's diagnosis predicted `κ = 1` per coefficient; measured **1.000000000000** on every coefficient of `rc_ladder(4)`, `(6)`, `(8)`. So the pipeline order is explained from our own numbers: s-expansion alone makes `κ` worse (§13), de-cancellation alone does nothing at fixed `s` (7c), **together they give `κ = 1`**. §16's "package" claim is now a positive demonstration rather than an inference from failures. |
| 7e — de-cancellation on a real amplifier | **7e-2 PASS, 7e-1 FAIL, and the FAIL is on my threshold rather than on feasibility.** `ua741(fully_symbolic=True)` added (132 device symbols, 343 addends, only 2 numeric; back-substitution matches the numeric fixture to 2.2e-16; two tests). States 27 577 → **204 679**, an overhead of **7.42×** against ladders' flat 1.20× — the declared reconsider-if was right that ladders are benign. But it is a *constant* factor, 204 679 states is small, and it runs in 14 s. Determinant preserved to **2.9e-14** on an active circuit. `κ` at fixed `s`: **6.66e3 → 99.05 (67×)**. |
| 7f — per coefficient, vs the right control | **The payoff.** Against the like-for-like control (same matrix, same point, s-expanded, *not* de-cancelled) **every one of 24 coefficients improves, by 60× to 394 367×**; the low-order ones land at **`κ ≈ 14-17`** (s^1: 1.85e3 → 17.3). First time in this thread that `κ` is in a range where magnitude ranking is plausible — 99.67% of the mass needed at `tol=0.05`, against 99.99947% at `κ=9.4e3`. `κ = 1` does not hold and was predicted not to: an amplifier is not passive. **Recorded self-correction:** the first reading used the compact *whole determinant* as the control and nearly wrote this up as a failure. |
| 7g — the payoff test | **FAIL, and it corrects 7f's reading.** At `tol=0.05`, control measured at the same point: compact group ranking **870 terms / 3.75e-02, converged**; compact magnitude **1 690 / 3.78e-02, converged**. De-cancelled: group **72 724 terms / 1.40e-01, NOT converged**; magnitude **68 310 / 1.59e-01, NOT converged**. So despite `κ` being 67× better, the de-cancelled expansion needs ≥80× more terms and does not converge. Cause: `κ` measures *conditioning*, term count depends on *concentration*, and de-cancellation spreads the same value from 2.77e6 terms over 1.1e21. **Low `κ` is necessary, not sufficient.** A bug found by gate 7g-4 on the way: a dead-end state (rows remaining, all children pruned) was counted as a completed term, giving an impossible error of 6.7e+55; caught by reasoning from the bound `Σ\|term\| = 99·\|det\|`, fixed, and re-verified against brute-force enumeration. Details: `cancellation_ranking_conclusions.md` §20. |
| 8 — the concentration diagnostic | **PASS on all four gates, and it settles §20.** `DDD.concentration(env)` returns the participation ratio `N_eff = A²/S2` in **one traversal**, the same cost as `cancellation`. Sanity: exactly the term count for equal-magnitude terms, 1.0 for a dominated one, within `[1, N]`, and within an order of magnitude of the *enumerated* 99%-mass count. **The discriminating test:** compact `N_eff = 194` (ranking took 870 terms, converged) against de-cancelled `N_eff = 11 565` (>72 724, did not converge) — correctly ordered and predictive to a factor of ~5, where `κ` ordered them backwards (6 659 against 99). Five tests; doc section with build-time numbers. |
| 9 — what reduces `N_eff`? | **Gate 9-2 FAIL, gate 9-3 PARTIAL — a real 5-12× gain, short of readable.** Ranking per coefficient of `s` instead of the whole determinant: nominal **`s^1`, 155 terms** against 734 (4.7×), carrying **97.3%** of the response; degraded **`s^2`, 206 terms** against 2 401 (11.7×). Device symbols intact throughout. Target was `N_eff ≤ 20` and ≤30 terms; best is 74.5 and 155. **And `N_eff` earned a caveat:** within one circuit its ordering across coefficients is *inverted* (126.8/341.2/709.2 → 950/275/206 terms), so it is a coarse screen for "is this representation hopeless", not a fine predictor. **Bug found and fixed:** `concentration` reported `N_eff = 0` above `s^12` because `Σ|term|²` underflows; rewritten scale-free on the weights, low-order values unchanged, test added. Details: `cancellation_ranking_conclusions.md` §22. |
| 10 — the tolerance curve | **Gate 10-2 FAILS on its threshold; the substance is the best result of the plan.** On the dominant `s^1` coefficient (97.3% of the response), varying `tol`: **5 groups / 11 operations / 26.3% error**, then 43 / 57 / 6.0%, 155 / 91 / 3.2%, 164 / 91 / **0.79%**. The gate wanted ≤30 terms at ≤20% and the table straddles it. **The eleven-operation expression, printed:** `5.9997e-70·gm_q17·(gm_q1 + 4.09e-4)·(−gm_q17 − 1.0023e-2)·(gm_q2 + 4.09e-4)` — four factors, three device symbols, and it reads as a circuit statement. **And the metric was wrong all along:** term count overstates the answer's size 4-8× because sympy collects shared factors — terms grow 33× across the curve while operations grow 8×. Every earlier "not readable" verdict in this plan was reached on term counts. Details: `cancellation_ranking_conclusions.md` §23. |
| 11 — transfer function, in operations, over a sweep | **Part 1 CORRECTS §23; part 2 PARTIAL — a 177-operation `H(s)` at one operating point.** Part 1: the terms-to-operations ratio is a property of the diagram, not a constant — 2.2 ops/group for an s-expanded coefficient but **69 ops/group for the compact whole determinant** (734 groups → **50 377 operations**), and 16 for the leapfrog top (181 → 2 895). So §23's "term count overstates 4-8×" is withdrawn: for the whole determinant the earlier verdicts were too *kind*. Part 2: exact `N/D` verified against the solve over the sweep (3.4e-15), then greedy coefficient choice with exact global re-evaluation. **Degraded: N `[0]`, D `[0,1]`, 177 operations, 7.9% error across two decades, device symbols intact — GATE 11-2 PASS.** Nominal fails: 119 ops at 26.7%, or 439 ops at 6.5%. The kept coefficient sets *differ* between operating points, which is the "different expressions for different symbol values" of the original brief. Details: `cancellation_ranking_conclusions.md` §24. |
| 12 — does 177 ops survive a real band? | **Gate 12-3 FAIL, gate 12-4 PARTIAL — it scopes stage 11's headline.** The µA741 falls at −20 dB/decade from below 100 Hz, so stage 11's window sat entirely on the single-pole rolloff. Over 10 Hz–10 MHz the error *does* hold but costs **6 439 operations against 177 — a 36× increase** (9 228 for 13.8%). My own advance estimate of "300-600 operations" was an order of magnitude optimistic. **The failure was mine and it is the fifth repeat:** the first wide run showed 101% error while keeping 27 of 46 coefficients — impossible if the subset were the problem — because the subset was chosen against the *global* error but each coefficient was then approximated at a **per-coefficient tolerance**, the exact per-piece budget §11 measured as unsound and §14c records being rejected in 2000. Separating the two tolerances fixes it. Details: `cancellation_ranking_conclusions.md` §25. |
| 13 — the `Short` operation | **13-1 PASS, 13-3 FAIL, 13-4 PARTIAL.** Calibrated against a hand-computable divider first (merged node at exactly `I·R2`). Greedy node merging with exact global re-evaluation accepts **8 merges** at 20% — more than the two or three predicted — taking the matrix **26 → 18**, terms **2 773 885 → 84 100 (33×)** and `N_eff` **151 → 41.9**. But the transfer function over the wide band goes only **6 439 → 4 081 operations (1.6×)** at the same 16.6% error. **A 33× smaller term pool buys 1.6× in expression size: the binding constraint is the accuracy demanded across the band, not term supply.** Every remaining lever of this kind attacks supply. **And I defeated my own safeguard:** tightening the coefficient tolerance made the error *worse* (16.6% → 145%) because three of sixteen coefficients hit `max_splits` and returned 58-87% error — each raising the `RuntimeWarning` built for exactly that, which a `simplefilter('ignore')` in the benchmark threw away. Fixed to record and report. Details: `cancellation_ranking_conclusions.md` §26. |
| 4 — library API | **Done.** `DDD.cancellation`, `DDD.subdiagram_values`, `DDD.minor_positions`, `DDD.approximate_groups` in `ddd.py`; twelve tests in `test_ddd.py`; three new subsections of `doc/src/circuit/ddd.rst` with every number generated at build time. `DDD.approximate` now **warns** when it returns without meeting `tol`. |

### Three results worth carrying forward

**1. `κ = Σ|term| / |Σ term|` is a cheap forecast of whether term-ranked
approximation will converge**, and costs one pass over the diagram. Capturing a
`1 − tol/κ` fraction of the absolute mass is *sufficient* for error `tol`; the
µA741's `κ = 9.4e3` at `tol = 0.05` makes that demand 99.99947%, and 500 of 2.77M
terms did leave 994% error.

**Corrected 2026-07-30 after reading the approximation literature:** the
condition is **sufficient, not necessary** — dropped terms may cancel among
themselves, so a large `κ` does not prove convergence impossible. The classical
stopping rule is the exact *signed* partial sum (Fernández et al. 2012,
criterion 8), which `DDD.approximate` already returns as its error; ranking by
magnitude only fixes the *order* of deletion. `κ` is also the textbook condition
number of a summation (Higham 1993), not a new quantity. See
`cancellation_ranking_conclusions.md` §14. The 994% remains a measurement; the
word "impossible" is withdrawn.

**2. The cancellation is concentrated near the root, not uniform.** Median `κ`
runs 9.4e3 at depth 0 down to 1.0 by depth 27. Term ranking fails on the whole
determinant while succeeding on almost every deep subdeterminant of it — which
is the condition under which ranking groups with exact values works. Measured:
734 groups at 4.1% error against 734 terms at 1186%.

**3. Where the readable-output problem actually lives.** The leapfrog's *top*
diagram barely cancels (`κ = 13.8`) and needs only 181 groups. The µA741's does
cancel hard (`κ = 9.4e3`). Since a leapfrog block **is** a µA741, the difficulty
is entirely **inside the blocks** — which is also the only place device
parameters exist. So group ranking is aimed at the right target and applied at
the wrong level.

### The next stage, named rather than left implicit

**Rank inside a suppressed block.** Each stamp symbol `_lvlN_a_b` is a
`DDDCombination` over cofactors of that block's internal matrix, whose entries
*are* device parameters. Group-ranking those diagrams — with the top-level
group's coefficient as the prefix — would put `gm_q17` back into the answer and
close the gap stage A of `doc/hierarchical_approximation_plan.md` opened. The
machinery needed is a `group_rank` that accepts an initial frontier of several
roots rather than a single one; everything else already exists.

**Reconsider-if this is not worth doing:** if a block's own `κ` turns out to be
near 1 after the `/D` normalisation the level walk applies, there is nothing for
group ranking to fix inside the block either, and the compact-but-anonymous form
is the best this representation offers.

### Two findings outside the plan's scope, recorded rather than dropped

**`HierarchicalDDD.eval` returns `0` on the leapfrog, and it is not wrong.**
`log10|det(A)| = -358.6` at nominal and `-372.9` at degraded, so the product of
111 level factors underflows double precision. The value is unrepresentable, not
miscomputed — but a caller cannot tell those apart from a returned `0.0`, and
`ZeroDivisionError` from a relative-error check is how it surfaced here. A
log-magnitude return, or a raise, would be the honest interface.

**The previously recorded hierarchical size does not reproduce.**
`doc/hierarchical_approximation_plan.md` §1 records the leapfrog's hierarchical
form as **299 vertices**; the run here measures **1958** (top diagram 1847) with
the same 16 terminals and the 111-node suppression. The fixtures differ — that
round had 121 symbols against 16 here — but more symbols should not give a
*smaller* diagram, so the discrepancy is **unexplained and left standing rather
than reconciled by assumption**. Re-derive before quoting either figure.
