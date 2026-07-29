# Cancellation-aware ranking — the reasoning

Why this document exists: `doc/hierarchical_approximation_plan.md` §5 records a
negative result that stopped the previous round — **term-ranked approximation
does not converge on a real operational amplifier** — and lists three ways
forward. This is the reasoning behind taking the second one:

> *Attack the cancellation rather than the representation. The failure is not
> hierarchy's; it is magnitude ranking's.*

The plan and its gates are in `doc/cancellation_ranking_plan.md`. This document
says why, and is meant to survive being argued with.

## 1. The measured failure, and what it is actually about

µA741, flat diagram, 1040 vertices standing for 2 773 885 product terms, three
symbolic transconductances, 1 kHz. Keeping the 500 largest terms leaves
**994%** relative error at nominal `gm` and **17 300%** at degraded `gm`.
`tol = 0.05` is never met.

The stated cause was "massive cancellation". That is right but too vague to act
on, so state it as an inequality.

Let `A[v]` be the **absolute-value companion** of the diagram:

    A[ONE]  = 1
    A[ZERO] = 0
    A[v]    = |entry(v)| * A[one(v)] + A[zero(v)]

`A[root]` is the sum of the absolute values of all 2.77M product terms, and —
this is the point — it costs **one pass over 1040 vertices**, because absolute
value distributes over a product. The exact determinant is `V[root]`, computed
by the same shaped pass with signs. Define

    κ[v] = A[v] / |V[v]|                    (the condition number of the sum)

Now: if `S` is any set of terms and the rest are dropped, the error is
`|Σ_{t ∉ S} t| ≤ A_tail(S)`, and magnitude ranking has no information beyond
`A`. So to guarantee relative error `tol`, a magnitude-ranked prefix must
satisfy `A_tail ≤ tol · |V[root]|`, i.e. it must capture

    a fraction  1 − tol/κ  of the total absolute mass.

At `κ = 10^6` and `tol = 0.05` that is 99.999995% of the absolute mass. **This
is not a tuning problem.** No reordering of a magnitude-ranked list fixes it,
no larger `max_terms` fixes it at any tractable count, and the 994% figure is
the expected consequence rather than a bug. The quantity that decides whether
term ranking can ever work on a circuit is `κ`, and nothing in the current code
measures it.

## 2. The object being ranked is wrong

A magnitude-ranked term carries one number: an upper bound on what it can
contribute, `|Π entries|`. That bound is only tight when nothing cancels — the
one case where approximation was never hard.

But a DDD node already knows the exact sum of every term below it. `V[v]` is
what `DDD.eval` computes; it is exact and it is free. So replace the ranked
object:

- a **term** is one 1-path, ranked by `|Π entries|` — an upper bound;
- a **group** is a pair `(P, v)` of a chosen prefix and a node, standing for
  *all* 1-paths through `v` with prefix `P`, and contributing
  `(Π P) · V[v]` — **exactly**.

Splitting a group at its node's two edges gives two groups whose values sum to
the parent's, so a frontier of groups always satisfies

    Σ over frontier of (Π P) · V[v]  =  the determinant, exactly.

Three consequences follow, and they are the whole argument for this route.

**Dropping becomes exact rather than bounded.** Drop a group and the error is
its value, which is known. A million mutually cancelling terms are dropped in
one operation at their true — tiny — cost. Magnitude ranking must keep all
million or lose the answer, because it can only see their individual sizes.

**Ranking becomes by true contribution.** `|(Π P) · V[v]|` is not a bound; it
is what the group contributes. The frontier is ranked by the thing that matters.

**Cancellation becomes a signal instead of a silent failure.** If splitting
`(P, v)` produces children far larger than the parent, the split *created* the
cancellation: `κ` rose across that step. That is locally measurable, and the
available response — stop splitting, keep the group whole — is exactly what a
magnitude ranker cannot do.

## 3. What a retained group *is*, and why it is not stage A's placeholder

Stage A of the previous plan produced a small expression over anonymous
hierarchical placeholders (`_lvl109_16_0`), and the maintainer's verdict was
"this was not a win": a placeholder is a Schur-complement rational function of
a suppressed block, and it names nothing a designer can act on.

A retained DDD group is a different object. A node reached by the builder
corresponds to a minor `(rows, cols)` of the **original** matrix, so the
retained form is

    det ≈ Σ_i  (product of named matrix entries) · det( M[rows_i, cols_i] )

Every symbol in it is a device parameter or `s`. The factored part is a
determinant of named entries, not an opaque reduction. So the output shape sits
between "sum of product terms" (readable, but needs the expansion to converge)
and stage A's placeholders (compact and meaningless), and it is symbolic in the
device parameters throughout — which is the maintainer's actual requirement.

**Reconsider-if:** if the retained minors turn out to be large — say above
6×6 — then "a determinant of named entries" is a formality and the designer is
no better off than with a placeholder. The measurement to make is the *size
distribution of retained minors*, and stage 2's gate is written around it.

## 4. The premise this route rests on, stated so it can be refuted

Everything above is machinery. Whether it works depends on one empirical fact
about real circuit matrices: **where the cancellation lives in the DAG.**

- **Concentrated.** `κ[root]` is large, but there is a shallow cut of the
  diagram on which most groups have `κ[v] = O(1)`. Then a few hundred groups
  with exact values reproduce the determinant, the retained minors are the
  honest carriers of the cancellation, and the route works.
- **Uniform.** `κ[v] ≈ κ[root]` at every depth: every subdiagram cancels as
  badly as the whole. Then no cut helps. Any expansion deep enough to expose
  device symbols has already destroyed the accuracy, and the scheme degenerates
  to either "keep the entire diagram" (exact, unreadable) or "drop everything"
  (wrong). **That would be a real negative result and it would apply to term
  ranking as such, not only to this implementation** — worth recording next to
  the five kept from the DDD work.

This is measurable in one pass, before any ranking code is written. It is
stage 0, and its gate is declared in the plan.

## 5. A precision hazard, up front

`V[v]` computed in double precision carries relative error about `κ[v] · 2^-53`.
So the scheme's "exact" group values are only exact where `κ` is moderate: at
`κ ≥ 10^14` they are noise, and a measurement taken there means nothing.

This project has been bitten by the neighbouring version of this before — H6's
scaling mitigation was borrowed from a *different* failure and did not apply,
and the fix was extended precision (`_suggested_precision`, `roots_of`). The
pattern to reuse exists. But note the difference, because rule 7 of the working
method exists for exactly this reason: **H4/H6 was float *root* conditioning;
this is float *summation* conditioning.** They are not the same failure and the
fix is only coincidentally the same tool.

Therefore stage 0 reports `κ` and the precision it implies *together*, and any
result where `κ · 2^-53 > 10^-6` is redone in `mpmath` before being believed.

**Reconsider-if:** if extended precision is needed at every vertex rather than
near the root, group values stop being free, the per-step cost rises by the
mpmath factor, and the frontier cap in stage 1 must come down accordingly.

## 6. Out of scope, and what would bring each back

**Moment matching / multi-point MOR.** Gives up device symbols for numbers.
The maintainer has already rejected numeric substitution once — "if we replace
all symbols with numbers I can just run SPICE". *Reconsider-if:* the group
route also fails to converge **and** the maintainer accepts numeric poles at
the filter level while keeping symbols inside one amplifier block.

**Direct truncation (DTT).** Passive RLC trees, numeric common poles. Already
retracted once in `doc/ddd_references.md` after being characterised from its
title. *Reconsider-if:* the target circuit becomes a passive ladder, which it
is not.

**Preconditioning the matrix to remove the cancellation.** Row and column
combinations can break the near-zero row sums that KCL builds into an MNA
matrix, which is plausibly where much of `κ` comes from. Excluded now because
it destroys the property section 3 is built on: after a basis change an entry is
a *combination* of device parameters, so the output names combinations rather
than devices. *Reconsider-if:* stage 0 finds cancellation uniform in depth —
in which case no cut exists, preconditioning is the only remaining lever, and
the interpretability cost has to be paid rather than avoided.

**Hierarchical diagrams.** Deliberately last, not first, even though the
leapfrog is the goal. The µA741 flat diagram is the only fixture where the
exact value, a full expansion, and the failing magnitude-ranked baseline all
exist at once, so it is the only place a new ranking can be *falsified*. A
scheme first demonstrated on the leapfrog would have no oracle — which is
precisely how stage B of the previous plan became unreachable.
*Reconsider-if:* stage 1 passes on the µA741; then the leapfrog is the next
stage, and it is stage 3 for that reason.

**Multi-corner ranking.** In scope eventually, not first. `DDD.iter_terms`
already accepts a sequence of environments and ranks by worst case, and the
group frontier extends the same way. Left out of stages 0–2 so that a failure
is attributable to the ranking rather than to the corner set.

## 7. How this could fail even if stage 0 passes

- **The frontier could grow without the error falling.** Group values sum to
  the determinant at every step, so nothing is ever *wrong*, but if each split
  produces two children of comparable size to their parent, the frontier grows
  while nothing becomes droppable. Then the honest output is a heuristic with a
  measured error distribution, not a guarantee — which is what most published
  symbolic approximation actually offers, and should be labelled as such.
- **Retained minors could be too large to read** (§3's reconsider-if).
- **The result could be exact but structureless**: a thousand groups each
  contributing a percent. Readability is the goal, so a term count without an
  error, or an error without a term count, is a misleading report. Both, always
  — this is the mistake §5 of the previous plan records being made once already.

## 8. What this document got right and wrong

Added after the stages ran, so the predictions above can be checked rather than
admired. Measurements: `doc/cancellation_ranking_plan.md` §Outcomes.

**Right, and the central claim.** §1's inequality is the correct account of the
failure. `κ = 9.4e3` on the µA741 at nominal, so `tol = 0.05` demanded
99.99947% of the absolute mass — and the measured error at 500 terms was 994%,
reproducing the previous round exactly. §2's group ranking then met the same
tolerance with 734 terms where magnitude ranking left 1186% at the *same* count.

**Right, and it was the load-bearing assumption.** §4 said the route lives or
dies on whether cancellation is concentrated or uniform. It is concentrated:
median `κ` falls monotonically from 9.4e3 at the root to 1.0 by depth 27.

**Wrong in emphasis.** §5 spent a paragraph on the precision hazard and it never
materialised: `κ·2^-53 = 1.6e-11`, so double precision was ample and mpmath was
not needed. The care was cheap and the reasoning was sound — a `κ` two orders of
magnitude larger would have needed it — but the hazard was over-weighted
relative to the one that actually bit, which was not on the list at all: **on the
leapfrog the determinant itself underflows** (`log10|det| = -358.6`), so
`HierarchicalDDD.eval` returns `0.0` and every relative error is undefined. The
conditioning of the *sum* was analysed; the dynamic range of the *product* was
not.

**Wrong about where the difficulty is.** §6 justified leaving hierarchical
diagrams last on the grounds that the µA741 is the only place a new ranking can
be falsified. That reasoning holds. But the implicit expectation — that the
leapfrog would be the hard case — is refuted: its *top* diagram has `κ = 13.8`
and needs 181 groups. The cancellation is inside the blocks, which is also the
only place device parameters live. So the honest statement is that group ranking
was demonstrated on the circuit that needed it and then applied at the level
that did not.

**Unresolved, and §3's reconsider-if is what to watch.** Retained minors turned
out to be readable at 6×6 and cut the item count 4-5×, so the factored form
earns its place. But 182 items is still not an expression a designer reads at a
glance, and no stage of this plan claimed it would be. The gap between
"converges" and "readable" is now the whole remaining distance.

## 9. Why the obvious recursion is hopeless, and what replaces it

Written before stage 5 runs, so the reasoning can be checked against it.

Stage 3's failure is precise: the leapfrog's groups name `_lvl110_*` stamps, not
devices. The plan's stated next step was "rank inside a suppressed block", and
the naive reading of that does not work.

**Why.** `suppression_order` chose **111 levels, one node each**. Level `k`'s
matrix entries are level `k-1`'s stamp symbols, so a `_lvl110_*` symbol is a
combination of cofactors of a matrix made of `_lvl109_*` symbols, and so on down.
Reaching a device parameter means unwinding all 111 levels — which is exactly the
unfolding the previous round measured at `count_ops = 2 256 398`. Recursing to
device symbols through single-node levels rebuilds the explosion that hierarchy
was introduced to avoid.

**What replaces it: choose blocks the size of the thing you want to name.**
Suppress each amplifier as *one* block of 22 internal nodes — five blocks, not
111 levels. Then the structure is two-deep instead of 111-deep, and the crucial
property is this:

> Amplifier 1's internal nodes do not couple to amplifier 0's internal nodes;
> they couple only through terminals. So block 0's elimination stamps **only**
> into terminal rows and columns, and block 1's `A_ii` is untouched — still
> device entries.

Every block's cofactor family is therefore over **device parameters directly**,
and only the top matrix carries stamp symbols. `det(A) = (∏ D_l)·det(reduced)`
where each `D_l = det(A_ii,l)` is a diagram over device entries, and each stamp
is a cofactor combination of the same. Three kinds of object, all
device-symbolic, all separately group-rankable.

**This is a hierarchy-design choice, not an algorithm change**, and it is the
one thing stage 3 got wrong. `suppression_order` optimises for *diagram size*,
which is the right objective for evaluation and the wrong one for
interpretability: it will happily choose a decomposition whose intermediate
quantities correspond to nothing a designer can name. **Interpretability
constrains the partition, not the ranking.**

**The output will be nested, and that has to be stated up front rather than
discovered.** Substituting the stamp expressions into the top-level groups would
multiply out to something astronomical — 17 factors each a sum of ~100 terms.
So the deliverable is a *factored, few-level* expression: a sum of products of
named sub-expressions, each itself a small sum of device-parameter products.
Whether that counts as readable is a fair question, and the honest answer is that
it is how a designer already reasons about a hierarchical circuit — "the gain is
set by this stage's transconductance through that stage's output conductance" —
whereas a flat sum of a million products is not.

**Reconsider-if:** if a block's own `κ` turns out near 1, group ranking buys
nothing *inside* the block either, and the honest conclusion is that the
compact-but-anonymous form was the best this representation offers. `κ` per block
is therefore the first thing stage 5 measures, before any composition.

## 10. Section 9's claim was refuted, and the refutation is narrow

Measured `benchmarks/cancellation_blocks.py`, first run. **Gate 2 failed:**
block 1's internal matrix came back carrying `_lvl0_*` stamp symbols, so the
five-block partition did *not* give device-level cofactors.

**The mathematics in §9 was right; the claim about the code was wrong.** Checked
separately: the number of matrix entries coupling one amplifier's internal nodes
to another's is **0**, so `A_ii` over the union of internal nodes really is block
diagonal, and the Schur complement really does factorise per amplifier.

What breaks it is `HierarchicalDDD._suppress`: it creates a fresh stamp symbol
for **every nonzero entry of the reduced matrix**, not only for the entries the
elimination modified. So after block 0 is suppressed, amplifier 1's internal
entries have been renamed too — multiplied by `det(A_ii,0)` and hidden behind
`_lvl0_*` — even though the elimination never touched them. Sequential
suppression therefore launders device entries into stamps regardless of
topology.

**So §9's conclusion survives and its route does not.** The fix is not to change
`_suppress` — that would alter every existing hierarchical result — but to stop
suppressing *sequentially*. Because the blocks are provably independent, all five
can be eliminated **in parallel** against the original matrix:

    det(A) = (∏_l D_l) · det( A_tt − Σ_l A_ti,l · adj(A_ii,l)/D_l · A_it,l )

Every `A_ii,l` is taken from `A` itself, so every cofactor is over device
entries by construction, and there is no sequence to launder them. This is a
*construction* in the benchmark script rather than a library change.

**The cost of interpretability, now quantified.** The five-block hierarchy is
**22 163 vertices** against 1 958 for the 111 single-node levels, and its top
diagram stands for 1 076 448 terms against 374 608 — an order of magnitude
larger. That is exactly the trade §9 predicted: `suppression_order` optimises
diagram size, and choosing blocks a designer can name costs about 11× in
representation. Worth paying, and worth stating rather than discovering.

**Kept as a negative result:** *a partition that is mathematically independent is
not automatically independent in the implementation.* The topology argument was
correct and still produced a false prediction, because the renaming step sits
between the topology and the result. Check what the code does to the symbols, not
only what the circuit does to the currents.

## 11. Stage 5's result: the cancellation problem is not compositional

Measured `benchmarks/cancellation_compose.py`. The parallel construction of §10
worked exactly as designed, and it produced a result that is more interesting
than a clean pass.

**What was won.** Gate 5 passes: the expressions name `gm_s0_q2`, `gm_s0_q17`,
`gm_s1_q2`, … — **device parameters of identified transistors in identified
stages**. That is the clause stage 3 failed, and it is a real advance. Gate 3
also confirmed the diagnosis that motivated the stage: `κ = 1.1e3` inside each
block against `13.8` at the top, so the cancellation really is where the devices
are.

**What was lost, and it is the finding.** Gate 4 asked for composed relative
error ≤ 1e-2 at *both* operating points. At nominal `gm` it is met. At degraded
`gm` it is met at **none** of three tolerance settings — and the failure is
**non-monotone**:

| per-piece tol (top/block/cofactor) | composed error, degraded | nested ops |
|---|---|---|
| 5e-2 / 2e-2 / 2e-2 | **1.47e-2** | 1 605 524 |
| 2e-2 / 5e-3 / 5e-3 | **1.13e-1** | 3 180 366 |
| 1e-2 / 1e-3 / 1e-3 | **1.48e-2** | 5 518 235 |

Tightening every piece by 4× made the composed answer **ten times worse**, and
tightening again brought it back. Meanwhile the operation count more than
tripled. **Per-piece tolerances do not control the composed error at all.**

**Why, and it is the same mechanism one level up.** Each reduced entry is a
weighted sum of 25 cofactors, and those cofactors cancel heavily against each
other — that is what `κ = 1.1e3` inside the block *means*. Truncating each
cofactor to a relative error δ leaves 25 residuals that are individually small
but need not cancel the way the values they came from do. So the composed error
is set by how the residuals happen to align, not by δ; and changing δ reshuffles
that alignment. This is precisely the failure §1 describes for magnitude
ranking — an error bounded per-part while the parts cancel — reappearing at the
block interface.

**So the honest statement is: hierarchical symbolic approximation is not
compositional.** Approximating the parts independently and assembling them gives
no guarantee about the whole, exactly when the parts cancel — which is exactly
when approximation is needed. Everything measured here is consistent with that,
and it explains why the published hierarchical-DDD literature reports
uninterpretable results rather than bad ones: the compositional route does not
carry an error bound to hand to a designer.

**The cure is the one this whole document argues for, applied once more.** Rank
the *combination* as a single object, with a frontier seeded by all 25 cofactor
roots at once, so the inter-cofactor cancellation lives **inside** the ranked
object and its contribution is exact. Then per-piece tolerances are not needed:
there is one ranking with one exact error. `group_rank` currently accepts a
single root; the multi-root frontier is the missing machinery, and it is now the
third time the same gap has been the answer.

**Two smaller results, both negative, both kept.**

- **Pruning cofactors by what the top ranking uses buys almost nothing.** The
  kept top groups touch only 42-43 of 96 reduced entries, which sounded
  promising, but those entries still reach **123-125 of the 125 cofactors**. The
  interface is too well connected for that kind of pruning. Implemented,
  measured, retained because it costs nothing — but it is not the lever.
- **Gate 6 was a badly designed gate and is worth flagging as such.** As first
  written it compared the nested form against the exact form's 4.11e+34 terms.
  Nothing could fail that. Restated against the meaningful reference — the
  2 256 398 operations the previous round spent unfolding a *single* µA741 — the
  nested form sits at 1.6M to 5.5M and straddles it. A gate whose bar cannot be
  missed measures nothing, and this one was set by reaching for whatever number
  was to hand.

**And the thing not to lose sight of:** even at the settings where gate 4 passes,
the answer is ~2 000 000 operations. It names devices, it is verified against
`slogdet` to 3e-3, and **no human will read it.** Converging and being readable
remain different problems, and only the first is solved.
