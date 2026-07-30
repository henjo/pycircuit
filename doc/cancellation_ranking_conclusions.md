# Cancellation-aware ranking — the reasoning

> ## READ THIS FIRST
>
> **This document grew by accretion over one session: §§1-13 were written as the
> work proceeded, and §§14-16 correct several of them.** Sections are never
> deleted, because a refuted premise is a result — but that means **an early
> section can state a claim this document no longer holds.** Every superseded
> section now carries a banner pointing forward. If you are here to find out what
> is true, use this table and do not read linearly.
>
> ### Status of every substantive claim
>
> | claim | where | status |
> |---|---|---|
> | `κ = Σ\|term\| / \|Σ term\|` forecasts whether magnitude ranking converges | §1 | **holds**, as a *sufficient* condition only |
> | "a ranking **must** capture `1 − tol/κ` of the absolute mass" | §1 | **CORRECTED** — sufficiency, not necessity (§14a) |
> | "the approximation is impossible as posed" | §1 | **WITHDRAWN** (§14a, §15b) |
> | `κ` is a novel diagnostic | §1 implied | **CORRECTED** — it is Higham's summation condition number (§14b) |
> | term-ranked approximation fails on a real op-amp | §1 | **SCOPED** — true of *our* pipeline; refuted on the same circuit by Tan & Shi 2004 (§15b) |
> | group ranking converges where magnitude ranking does not | §2 | **holds, measured** (734 terms at 4.1% vs 1186%) |
> | a retained group is a named minor, unlike a placeholder | §3 | **holds** |
> | cancellation is concentrated, not uniform | §4 | **holds, measured** |
> | extended precision will be needed | §5 | **did not materialise** (§8) |
> | the leapfrog needs amplifier-sized blocks | §9 | **holds** |
> | sequential block suppression keeps device symbols | §9 | **REFUTED** (§10) |
> | hierarchical approximation is not compositional | §11 | **holds, measured** — and was published in 2000 (§14c) |
> | s-expansion buys cancellation-freedom | §13 implied | **REFUTED** (§13's own measurement, §15c) |
> | de-cancellation caps the `κ` gain at 3-10× | §14e | **WITHDRAWN** — conflates term count with mass (§15c) |
> | de-cancellation is the load-bearing step | §15c | **holds, provable** |
> | element pruning shrinks a compact-symbol diagram | §16 hypothesis | **REFUTED, measured** (§16) |
> | de-cancellation gives `κ = 1` on a calibration circuit | §17 | **holds, exactly** |
> | determinant-side de-cancellation is not a compact representation | §17 | **CORRECTED** — true of a *path*-keyed state, false for a label-keyed one, which costs a flat 1.2× (§18) |
> | de-cancellation alone reduces `κ` | §17 implied | **REFUTED** — no benefit at a fixed frequency, worse above N=6; the residue is phase (§18) |
> | de-cancellation **plus s-expansion** gives `κ = 1` per coefficient | §18 | **holds for PASSIVE networks**, measured to 12 figures; **not** for an amplifier (§19) |
> | the 1.2× state overhead carries to a real amplifier | §18 reconsider-if | **REFUTED** — 7.42× on the µA741, though still a constant factor (§19) |
> | de-cancellation brings `κ` into a workable range on a real amplifier | §19 | `κ` figures **hold**; the inference that ranking would therefore work is **REFUTED** (§20) |
> | low `κ` means a ranking will converge in few terms | §19 inference | **REFUTED, measured** — `κ` is conditioning, term count is *concentration*; de-cancellation dilutes 2.77e6 terms into 1.1e21 (§20) |
> | de-cancellation is the route to a readable expression | §13, §18, §19 | **REFUTED** (§20). The compact diagram with group ranking remains the best route: 870 terms at 3.7% |
>
> ### The three things worth carrying away
>
> 1. **`κ` is cheap and worth reporting**, and it measures the *formulation* at
>    least as much as the circuit (§12). Nothing in symbolic circuit analysis
>    appears to report it.
> 2. **Do not account error per term or per block. Evaluate the candidate
>    answer.** Both literature sweeps converge on this (§14c, §15d), and it is what
>    made stage 6's greedy selection monotone where stage 5's per-piece tolerances
>    were non-monotone.
> 3. **Full symbols, de-cancellation and element pruning are a package** — each
>    one measured alone here, each failing or backfiring alone (§16), and §18 then
>    demonstrated the combination *positively*: de-cancellation with s-expansion
>    gives `κ = 1.000000000000` per coefficient. The remaining blocker is a
>    **fully-symbolic fixture** (the µA741 bakes 159 of 215 addends into numbers),
>    which is plumbing rather than research.
>
> ### The metric caveat, and its own correction
>
> §§16-22 judge readability by **term count**. §23 found that understates the answer
> for an s-expanded coefficient (5 groups → 11 operations) and generalised the claim;
> **§24 withdrew the generalisation.** The ratio is a property of the diagram: ~2
> operations per group for a coefficient diagram, **~69 for the compact whole
> determinant** (734 groups → 50 377 operations). So the earlier verdicts on the
> whole determinant were too *kind*, not too harsh. Report operations, and do not
> assume a ratio.
>
> ### The state of the goal, in one line
>
> **Converging is solved** — group ranking on the compact diagram, shipped and
> tested. **And for one coefficient, so is readability:** the µA741's dominant `s^1`
> coefficient reduces to **11 operations in three device symbols at 26% error**, or
> **91 operations at 0.79%** (§23). That is the original brief — device-symbolic,
> readable, stated error, different at different operating points — met on a
> component of the problem.
>
> **Not met:** the complete transfer function (24 coefficients plus a numerator), and
> the leapfrog, where §11 showed the composition is not error-controlled. Two routes
> are closed by measurement: de-cancellation (§20 — improves `κ` 67× while diluting
> 2.77e6 terms into 1.1e21) and element pruning of a compact-entry diagram (§16).
> Both `κ` (`cancellation`) and `N_eff` (`concentration`) are shipped; `N_eff` is a
> coarse screen only (§22).
>
> ### Where the code is
>
> `DDD.cancellation`, `DDD.approximate_groups`, `DDD.subdiagram_values`,
> `DDD.minor_positions` in `pycircuit/circuit/ddd.py`; user-facing write-up in
> `doc/src/circuit/ddd.rst` (already carries the corrected claims); gates and
> outcomes in `doc/cancellation_ranking_plan.md`; five benchmark scripts
> `benchmarks/cancellation_*.py` and `benchmarks/element_ranking.py`.

Why this document exists: `doc/hierarchical_approximation_plan.md` §5 records a
negative result that stopped the previous round — **term-ranked approximation
does not converge on a real operational amplifier** — and lists three ways
forward. This is the reasoning behind taking the second one:

> *Attack the cancellation rather than the representation. The failure is not
> hierarchy's; it is magnitude ranking's.*

The plan and its gates are in `doc/cancellation_ranking_plan.md`. This document
says why, and is meant to survive being argued with.

## 1. The measured failure, and what it is actually about

> **CORRECTED BY §14a, §14b AND §15b — do not quote this section alone.** The
> inequality below is *sufficient*, not necessary: dropped terms may cancel among
> themselves. "Impossible" is withdrawn. `κ` is Higham's summation condition
> number, not a new quantity. And "does not converge on a real operational
> amplifier" is true of *this pipeline* only — Tan & Shi (TCAD 2004) report the
> opposite on the same µA741.

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

> **DID NOT MATERIALISE — see §8.** `κ·2^-53 = 1.6e-11` on the µA741, so double
> precision was ample and mpmath was never needed. The hazard that actually bit
> was the *dynamic range of the product* (the leapfrog determinant underflows),
> which is not on this list.

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

> **PARTLY REFUTED BY §10.** The topology argument is correct — 0 entries couple
> one amplifier's internals to another's — but the claim that *sequential*
> five-block suppression preserves device symbols is false: `_suppress` renames
> every nonzero of the reduced matrix. The fix is parallel elimination, in §10.

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

> **STILL HOLDS, BUT IT IS NOT NEW — see §14c.** Guerra et al. (DATE 2000)
> rejected per-block error budgets in advance for this reason, and their fix
> (score against a global numerical reference) is the one to adopt.

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

## 12. This was published in 2012, and the paper was already on the shelf

The maintainer asked whether I had read Song & Shi, *Hierarchical Graph
Reduction Approach to Symbolic Circuit Analysis with Data Sharing and
Cancellation-Free Properties* (ASP-DAC 2012). I had not. It is in
`doc/ddd_references.md`, the PDF is local, and this session opened that very
file — for its µA725 schematic — without reading its cancellation content.

**What it says, and it is the diagnosis of everything above.**

1. **§II shows the cancellation is a property of the formulation, not of the
   circuit.** Its example is our own compact-symbol MNA convention. With
   `a = G1`, `b = c = -G1`, `d = G1+G2+G3`, `e = f = -G3`, `g = G3+G4`, the
   determinant reads `adg - aef - bcg`; substituting the composites gives

       G1(G1+G2+G3)(G3+G4) - G1*G3² - G1²(G3+G4)  =  G1G2G3 + G1G2G4 + G1G3G4

   Three terms, **all positive**: `κ = 1` exactly. The cancelling terms are
   manufactured by expanding composite entries that share conductances with
   opposite signs. **So §1's `κ = 9.4e3` characterises MNA-with-composite-symbols
   applied to a µA741, not a µA741.** Everything measured here remains true of
   the representation pycircuit has; the claim that had to be withdrawn is the
   stronger framing that the approximation was "impossible as posed".

2. **§III is titled "Schur Decomposition versus Symbolic Stamp"** — precisely the
   choice §9 made — and states the result of §11 directly: *"Since DDD is not
   cancellation-free, every layer of the hierarchy has the term cancellation
   problem."* §11's "hierarchical symbolic approximation is not compositional" is
   a **rediscovery**. What this work adds is quantification: `κ` as a
   one-pass predictor, and the non-monotone composition
   (1.47e-2 → 1.13e-1 → 1.48e-2) that shows per-piece tolerances are not merely
   loose but *unsound*.

3. **GPDD gets the interpretability property for free.** Its expressions are
   "directly composed of the circuit parameters, rather than intermediate
   symbols" — which is exactly what stages 3–5 chased through `_lvl110_*` and
   `T_a_b` placeholders and never obtained.

**Why it was skipped, which is the part worth keeping.** The reference entry
said: *"Relevant only if we take the graph route; listed for completeness."* The
GPDD entry said its comparable runtime was *"a real argument for staying with
DDD."* Both judgements were made while choosing a data structure, against a
speed-and-memory criterion, and both are defensible **for that question**. They
were then carried forward unexamined into a different question — whether
approximation can converge — where cancellation-freedom is not a nice-to-have but
the whole issue. **The failure was not in reading the papers wrongly; it was in
trusting a filed conclusion after the question changed.** A reference note should
record what a paper *says* separately from what we *decided*, because the
decision expires and the facts do not.

**The cost of the property, measured** (Shi, TCAD 2013, Table VIII): GPDD sizes
10 432 / 17 488 / 197 274 against DDD's 3 579 / 11 506 / 62 794 — about 3× — with
runtimes 0.682 / 0.793 / 6.771 against 0.586 / 2.042 / 10.359, so comparable and
sometimes faster. Three times the diagram for a guarantee that decides whether
dominant-term extraction works at all is cheap.

**What does not change:** the µA725 decline stands. Song & Shi's Tables I–II give
that circuit's module *partition*, not its connectivity, and Fig. 4 still draws
no junction dots.

## 13. A cheaper route than GPDD, and a measurement that locates the target

> **ITS RECOMMENDATION IS SUPERSEDED BY §14e, §15c AND §16.** The measurement in
> this section stands and is important (all cancellation is intra-coefficient,
> `κ` between powers of `s` is 1.03). The conclusion drawn from it — that
> de-cancellation alone is the cheap win — was too kind: de-cancellation is
> necessary but not sufficient, and it only pays together with full symbols and
> element pruning (§16).

Searching the local paper set turned up **Tan & Shi, "Parametric Analog
Behavioral Modeling Based on Cancellation-Free DDDs" (BMAS 2002)** — held here
since before this work started, unread. It matters because it obtains
cancellation-freedom **inside the determinant framework**, with no graph model,
which is the expensive part of the GPDD route.

**Its distinction, which this project had been eliding.**

* **Symbolic cancellation** — terms that cancel *as symbols*, "aris[ing] from the
  use of the MNA formulation and device matching in analog circuits". Their
  example: with `g = k = 1/R3` and `i = j = -1/R3`, the term `agks^0` cancels
  `-ajis^0`. These are detectable from **local 2x2 matrix patterns** (their
  Fig. 5 — e.g. `[[p,-p],[-p,p]]`, which "may come from the rectangular
  appearance of a floating resistor in the nodal admittance formulation").
* **Numerical cancellation** — terms that cancel only once numbers are supplied.
  They explicitly defer this "to the dominant term generation step".

**`DDD.cancellation` measures the numerical kind.** So the two are not
interchangeable, and the question that decides whether de-cancellation is worth
building is how much of `κ` is the symbolic part.

**Their method, and the two costs they state.** De-cancelling *after* building
the s-expanded DDD "may take a huge amount of temporary memory ... and is
typically the most time consuming step in the whole approximation operation".
Better is to remove the terms *during* construction, via a "canceling label list"
`CL(L_x)` per label and a modified coefficient-multiply that never generates them.
Two honest caveats from the paper: "up to 70-90 percent product terms can be
canceling terms for a typical analog circuit" (the prize), and "the
cancellation-free s-expanded DDDs **can be larger** than the normal s-expanded
DDDs" (the price).

### The measurement: where the cancellation actually sits

`benchmarks/cancellation_symbolic_vs_numeric.py`, µA741 at 1 kHz. Reassembly of
the s-expanded form against the complex DDD checks out at 2.3e-15, so the
comparison is sound.

| representation | `κ` (nominal) | `κ` (degraded) |
|---|---|---|
| compact-symbol complex DDD | 9.38e+03 | 1.40e+05 |
| **between powers of `s` only** | **1.028** | **1.116** |
| s-expanded coefficients, best | 1.57e+03 (k=1) | 4.45e+03 (k=2) |
| s-expanded coefficients, worst | 3.09e+06 (k=23) | 3.09e+06 (k=23) |

**Two results, one encouraging and one not.**

1. **Essentially none of the cancellation is between powers of `s`** (`κ = 1.03`).
   All of it lives *within* a single coefficient — which is precisely the
   MNA-and-device-matching kind Tan & Shi remove by construction. So the
   cancellation this project has been fighting is, by their taxonomy, the
   **removable** kind. That is the strongest argument yet for the de-cancellation
   route.
2. **s-expansion alone does not deliver it, and mostly makes `κ` worse.**
   Per-coefficient `κ` runs from 1.5e3 to 3.1e6 against the compact form's
   9.4e3. The dominant coefficient at this frequency (k=1) is 6× better than the
   compact determinant, which is a real but modest gain; the high-order
   coefficients are far worse, though they carry negligible magnitude at 1 kHz.
   **Having `s_expand` is therefore not the same as having the property** — the
   de-cancellation step is the substance, not the s-expansion.

### What this does to the route comparison

| route | cancellation removed | needs | size cost | reuses pycircuit |
|---|---|---|---|---|
| group ranking (built) | none — works around it | nothing | none | yes |
| **de-cancelled sDDD** | **symbolic only** | modified coefficient construction in `_SBuilder` | "can be larger" | **yes — `s_expand` is the hook** |
| GPDD | symbolic and numerical, by construction | a circuit **graph model** | ~3× | no |

**Recommendation, revised from "prototype GPDD".** De-cancellation is the cheaper
experiment and the measurement above says its target is the right one. It also
has a natural gate: apply it and re-measure `κ` per coefficient. If `κ` drops to
O(1), the whole approximation problem changes character and group ranking becomes
unnecessary. If it barely moves, the residue is numerical cancellation, GPDD's
extra guarantee is what is actually needed, and the graph model has to be paid
for. **Either outcome is worth having, and the de-cancellation experiment is the
one that distinguishes them.**

*Reconsider-if:* the `CL` construction turns out to need the full-symbol
(one-vertex-per-device) representation rather than pycircuit's compact-symbol
entries. Their Fig. 5 patterns are stated over matrix entries, and their sDDD
introduces "a unique symbol for each circuit parameter", so this is a live risk
and is the first thing to check when implementing.

## 14. A literature sweep, and it corrects §1's logic

A search of the approximation literature outside Shi's own corpus returned four
things that change what is written above. Sources and quotations:
`doc/ddd_references.md` §F.

### 14a. §1's inequality was stated as a necessity. It is only a sufficiency.

§1 says a magnitude ranking "**must** capture a fraction `1 - tol/κ` of the
total absolute mass". That is wrong as stated. The error of a truncation is
`|Σ dropped|`, and absolute mass only *bounds* it — the dropped terms may cancel
**among themselves**, in which case the true error is far below their mass and a
ranking can converge at large `κ`.

The classical stopping rule makes this explicit. Fernández, Sánchez-López,
Castro-López & Roca-Moreno, *Approximation Techniques in Symbolic Circuit
Analysis* (Bentham, 2012), give three successive criteria; the standard one,
their (8), is

    |Σ_{l=1..P} h_kl(x_o)| / |Σ_{l=1..T} h_kl(x_o)|  <  ε_M

a **signed** partial sum of the deleted terms over the signed total. Their own
note on it: *"Mutually canceling terms do not contribute to (8) because they are
added with their respective signs."* Their (9), which normalises by
`Σ|h_kl|`, is a robustness-to-parameter-variation criterion and **not** a
relative-error bound — satisfying it at 5% permits an actual error up to `0.05κ`.

**`DDD.approximate` already implements criterion (8)**: its returned `err` is
`|Σ kept − exact| / |exact|`. So pycircuit's stopping test was never the problem,
and the field does not rank by absolute mass either — magnitude ranking decides
only the *order* of deletion. Corrected in the `cancellation` docstring and in
`ddd.rst`.

**What survives.** `κ` remains a cheap and useful forecast, and the measured 994%
at 500 terms still stands — but it is a *measurement*, not a consequence of `κ`.
The honest claim is: a large `κ` means the absolute-mass argument gives no
comfort, so convergence has to be checked rather than assumed. The word
"impossible" was too strong and is withdrawn.

### 14b. `κ` is not a new quantity

It is the textbook condition number of a summation, `Σ|x_i| / |Σ x_i|`, around
which the whole compensated-summation literature is organised (Higham, *The
Accuracy of Floating Point Summation*, SIAM J. Sci. Comput. 14(4), 1993; Ogita,
Rump & Oishi 2005). Now cited. The contribution here is the *application* — using
it as an a-priori forecast for symbolic term ranking, and computing it on a
diagram in one pass — not the quantity.

Practical corollary: compensated summation would fix *evaluation* of a term list
at `κ = 9.4e3` (naive summation loses ~4 digits), and does nothing for
*truncation*. Two separate problems; only the first has an off-the-shelf fix.

### 14c. Someone already rejected per-piece error budgets, for our reason

§11 concluded that hierarchical approximation is not compositional. Guerra, Roca,
Fernández & Rodríguez-Vázquez, *A Hierarchical Approach for the Symbolic Analysis
of Large Analog Integrated Circuits* (DATE 2000), considered exactly that design
and rejected it in advance:

> *"a separate application to each block would require an error propagation
> mechanism at this early stage of the analysis process. This necessarily yields
> more conservative results (less reduced circuits) and, consequently, has a
> negative impact on the global performance of the analysis methodology."*

**And their alternative is directly implementable here.** They use hierarchy for
*generation only*, and control error against a **global numerical reference for
the whole circuit** — never as a per-block budget:

> *"The introduced errors are a combination of the error in each
> (trans)admittance (only part of it has been generated) and the contribution of
> such (trans)admittance to the global behavior."*

That is the fix for stage 5's non-monotonicity: stop assigning tolerances per
piece, approximate greedily, and accept or reject each candidate truncation by
its effect on the *global* answer. Their stated limitation is the familiar one:
*"this is guaranteed only at the selected frequency sample."*

The 2012 survey (14a) contains **no hierarchical section at all**, and the one
paper claiming a hierarchical formula-approximation algorithm — Fernández,
Rodríguez-Vázquez, Martín & Huertas, AICSP 3(1), 1993 — could not be obtained.
That is the outstanding gap.

### 14d. The one positive composition result is in model order reduction

Reis & Stykel, *Stability analysis and model order reduction of coupled systems*
(MCMDS 13(5), 2007), reduce subsystems and recouple them through the original
interconnection and *"obtain error bounds for the reduced-order closed-loop
system in terms of the errors in the reduced-order subsystems"* — the shape of
theorem §11 wants. Abstract only; the hypotheses are unverified. Such bounds are
norm-based with a stability side condition, so a per-block *term tolerance* is not
the input they take — which is a plausible explanation for why tightening term
tolerances uniformly made things worse.

### 14e. And a correction to §13's optimism — ITSELF SUPERSEDED BY §15c

**The 3-10× estimate below is withdrawn.** It converts a term-*count* fraction
into a *mass* fraction, and those are different quantities: an exactly-cancelling
pair adds `2|t|` to the mass and `0` to the value, so if the cancelling pairs are
the large terms, removing them collapses `κ` by far more than their count share.
See §15c. The rest of the subsection — the exact/near distinction — stands.

§13 said the µA741's cancellation is "by their taxonomy, the **removable** kind".
That needs qualifying. Every "cancellation-free" result in this literature —
Tan & Shi's patterns, Tan/Guo/Qi's hierarchical de-cancellation, Filaretov's
nullor determinants — concerns **exact** cancellation: term pairs that sum to
identically zero and can be removed structurally. `κ = 9.4e3` is **near**-
cancellation among the terms that survive that removal.

Quantitatively: if 70-90% of terms are exactly cancelling (Tan & Shi's figure),
removing them cuts the absolute mass by that factor while leaving the value
unchanged, so `κ` would fall by roughly 3-10× — from 9.4e3 to order 10**3, not to
order 1. **So de-cancellation alone would not deliver the property, and §13's
route comparison was too kind to it.** Worth doing for the size reduction; not a
solution to the ranking problem.

Also recorded: the sweep found **nobody who treats near-cancellation as such** in
symbolic circuit analysis. That makes `κ` more useful than §14b's deflation
suggests — the quantity is standard, but nothing in this field appears to measure
or report it.

## 15. The Shi/Tan sweep: a published counter-example, and §14e was also wrong

A second sweep, over Guoyong Shi's and Sheldon Tan's corpora. Sources in
`doc/ddd_references.md` §G. Note there are **three** researchers routinely
conflated: Guoyong Shi (SJTU, GPDD), C.-J. Richard Shi (UW, original DDD), and
Sheldon X.-D. Tan (UCR). The cancellation-free-DDD thread is Tan's.

### 15a. Cancellation-free DDD is established, and stays in the determinant model

Three mechanisms, all determinant-side, none needing a graph model:

1. **An index predicate — no numerics at all.** Tan, Qi & Li (DATE 2004),
   Theorem 2: *"For a given product term from a determinant, which consists k
   first-order cofactor ... k ≥ 2, if there are two first-order cofactors that
   share the same row index or column index, then there exists another product
   term which will cancel with this product term."* Root cause named in the same
   paper: *"Term cancellation ... will happen when MNA formulation is used where
   each device admittance may appear more than once."* This is a combinatorial
   test on index sets, implementable against our matrix directly.
2. **Hierarchical Schur, made cancellation-free.** Tan, Guo & Qi (DAC 2004),
   Theorem 1, with partial DDDs plus a complementary DDD per subcircuit,
   producing a **flat** cancellation-free DDD. Their stated cost: *"the new
   construction method will lead to larger DDD size than non-hierarchical method
   in general."*
3. **De-cancellation on the s-expanded multi-root DDD.** Tan & Shi (TCAD 2004).
   Stated cost, and it is a real one for us: *"cancellation-free s-expanded DDDs
   do not satisfy Theorem 1"* — de-cancellation **destroys DDD canonicity**, on
   which the efficient graph operations depend. Their own workaround is to run
   term generation *before* de-cancelling.

**A contradiction in the literature, flagged rather than resolved.** Song & Shi
(ASP-DAC 2012) motivate hierarchical GPDD by asserting that the DDD-based
hierarchical method *"does not have the cancellation-free property"* — while
Tan, Guo & Qi (DAC 2004) prove that hierarchical DDD **is** cancellation-free
once de-cancellation is applied, and Guoyong Shi's own survey lists
"DDD → De-cancellation, Exact" as a peer method citing that paper. Song & Shi do
not cite or rebut it. Per the standing rule about sources disagreeing, this
bounds how much weight either claim can carry.

### 15b. A published counter-example to this project's headline claim

`doc/hierarchical_approximation_plan.md` §5 records, and this document repeats:
*"term-ranked approximation does not converge on a real operational amplifier."*

**Tan & Shi (TCAD 23(6), 2004) report the opposite, on the µA741.** For the `s^1`
denominator coefficient of a two-stage op-amp: *"the first product term amounts
to 86% of the total magnitude of the coefficient and the first two terms amount
to 97%."* Their pipeline differs from ours in three ways at once:

* they **de-cancel** first (*"70-90% terms are canceling terms"*);
* they rank **within each coefficient of `s^k`** on a multi-root DDD, never over
  the whole determinant;
* their symbols are **individual circuit parameters**, not composite matrix
  entries.

So our negative result is **scoped to our pipeline** — compact composite symbols,
no de-cancellation, ranking the whole determinant at one frequency — and is not a
property of the µA741. That is the second time this claim has had to be narrowed;
it should not be restated without the scope attached.

Also relevant: Guoyong Shi's survey defines the term as *"the existence of
identical product terms with opposite term signs"* and judges that
*"numerically speaking, term cancellation might not be a very serious problem"* —
consistent with §14e's point that his "cancellation" and our `κ` are different
quantities.

### 15c. §14e's estimate was unsound, and the correct statement is provable

§14e argued that since 70-90% of *terms* cancel exactly, removing them cuts `κ`
by only 3-10×. **That converts a term-count fraction into a mass fraction, and
those are not the same quantity.** An exactly-cancelling pair contributes
`2|t|` to the absolute mass and exactly `0` to the value, so if the cancelling
pairs are the *large* terms, removing them collapses `A` and hence `κ` by far
more than their count share. The number 3-10× is withdrawn; the effect's size is
not predictable from a term count and has to be measured.

What *is* provable, and worth stating because it settles a question §13 left
open: since `|a + b| <= |a| + |b|`, expanding composite entries into sums of
device symbols can only **increase** the absolute mass while leaving the
determinant unchanged. So

    κ(full symbols, no de-cancellation)  >=  κ(composite entries)

**Expanding to one-symbol-per-device makes `κ` worse, not better, on its own.**
This is consistent with §13's measurement (s-expansion alone raised per-coefficient
`κ` to 1.5e3-3.1e6) and it explains it. The gain in Song & Shi's ladder example —
`adg - aef - bcg` becoming three positive terms — comes entirely from the
**de-cancellation**, not from the change of symbol. Therefore: **de-cancellation
is the load-bearing step, and full symbols without it are strictly worse.** Any
attempt at this must do both.

### 15d. The idea that sidesteps `κ` entirely

Hu, Shi, Tai & Lee, *Topological Symbolic Simplification for Analog Design*
(ISCAS 2015) — open access — **ranks circuit elements rather than product
terms.** Each symbol is tried under two GPDD edge operations, Short and Open; the
reduced transfer function is evaluated **exactly**; the candidate is scored by RMS
relative error on dc gain and phase margin; symbols are sorted by that error and
the worst-scoring `K` removed. A folded-cascode op-amp goes from 123 symbols to
18 in 3.9 s, and the reduced dc-gain expressions **match the textbook
hand-derived formulas**.

**This is immune to `κ` by construction**, because nothing ever sums term
magnitudes — every candidate is judged by an exact evaluation of a reduced
network. It is also the only method found anywhere whose output was checked
against a formula a human wrote. Their exclusivity claim, which would have to be
tested rather than believed: *"The property described so far for element
operation is a unique property owned by GPDD only; DDD does not have an analogous
property."*

Note the convergence with §14c: Guerra et al. score candidates against a global
numerical reference; Hu et al. score element removals against an exact evaluation.
**Both abandon per-term error accounting in favour of evaluating the candidate
answer.** That is the consistent lesson of both sweeps.

### 15e. What to do, revised again

In cost order, cheapest first:

1. **Element-level ranking scored by exact evaluation** (§15d, §14c). Needs no new
   representation: pycircuit can already evaluate a modified circuit exactly. This
   is the highest value per unit of work and it is representation-agnostic — the
   GPDD-exclusivity claim applies to the *edge operations*, not to the idea of
   scoring an element's removal by an exact solve.
2. **Reproduce Tan & Shi's pipeline** (§15b): de-cancel, then rank per `s^k`
   coefficient. Either it fixes our case or it refutes a published result, and
   both outcomes are worth having. Note the canonicity cost in §15a(3).
3. **GPDD**, only if 1 and 2 both fail. And note the cost is not free for our
   circuit class: Kolka et al. record that for VCCS-heavy circuits — i.e.
   transistor circuits — the two-graph intersection rejects an exponentially
   growing fraction of candidates.

*Reconsider-if:* someone measures `κ` for a GPDD expansion of an op-amp. Nobody
in either sweep did. "Cancellation-free" guarantees no identical-opposite pairs;
it does not guarantee `κ ≈ 1` once controlled sources make tree-pair products
signed. **That single number would settle the architecture question**, and it is
currently unknown.

## 16. Stage 6: element ranking works, and does nothing for term ranking

Measured `benchmarks/element_ranking.py`, µA741 with all 24 transistor `gm`
symbolic, worst relative error of the transfer function over 21 points from 10 Hz
to 1 MHz, greedy with exact joint re-evaluation at every step.

| gate | result |
|---|---|
| 1 — under 5 min | **PASS** — 71 s to split the matrix, 0.3-0.5 s to rank |
| 2 — how many symbols drop | **5 of 24** at 1% joint error, **7 of 24** at 5% |
| 3 — `κ` falls by ≥ 10× | **FAIL** — 1.00× at 1%, and **0.61×** at 5%, i.e. `κ` got *worse* |
| 4 — named survivors, smaller diagram | names yes; **diagram identical**: 1040 vertices and 2 773 885 terms in both cases |
| 5 — monotone joint error | **yes**, in both runs |

The order is exactly what a designer would expect — `gm_q22`, `gm_q24`, `gm_q21`,
`gm_q15`, `gm_q13` go first, at costs from 1.7e-08 to 2.0e-03, and they are bias
devices. The method is fast, the answer is readable, and the greedy joint error is
monotone, so the §11 composition trap really is designed out by re-evaluating
exactly.

**But gate 3 failed, and gate 4 failed in a way that explains why.** The diagram
of the "simplified" circuit is *bit-for-bit the same size* as the full one. The
reason is structural and it generalises:

> A compact-symbol DDD has one vertex per **matrix entry**, and an entry is a sum
> like `gm_q22 + g_o + s*c`. Setting one device's parameter to zero removes an
> *addend inside* an entry; the entry stays nonzero. So the matrix's sparsity
> pattern, the diagram, and the term count are all untouched.

**Element pruning and the compact-symbol representation do not compose.** Dropping
seven devices changed the term count by nothing at all and moved `κ` in the wrong
direction. This was declared as a possible outcome before the run — *"it is
entirely possible that `κ` barely moves while the symbol count drops a lot"* — and
it is what happened.

### What this means, and it ties the whole document together

Three ingredients have now each been measured in isolation, and **each one fails
or backfires alone**:

| ingredient | alone | why |
|---|---|---|
| one symbol per device (full symbols) | `κ` **worse** | §15c: `|a+b| <= |a|+|b|`, so expanding entries only adds absolute mass |
| de-cancellation | untested here | removes the mass that expansion added; §15a |
| element pruning | diagram **unchanged** | §16: zeroing an addend inside a composite entry removes no vertex |

They are a **package**, and Tan & Shi's pipeline contains all three for that
reason: per-device symbols make elements individually addressable *and* inflate
the term count, de-cancellation removes the inflation, and element pruning then
actually deletes vertices. Adopting any one of them piecemeal — which is what this
project did three times over — cannot work, and now there is a measurement for
each.

**So the revised order of work in §15e is itself revised.** Element ranking is
worth keeping as a *circuit-reduction and design-insight* tool: 5-7 devices
identified as irrelevant, with an exact error, in under a second. It is **not** a
route to a smaller symbolic expression in this representation, and it should not
be sold as one. The next step that could actually pay is the full-symbol
construction *with* de-cancellation — items 2 of §15e — because that is the
smallest subset of the package that is self-consistent.

*Reconsider-if:* the `Short` operation (parameter to infinity / node merging),
which was not implemented here, does collapse nodes and therefore *does* shrink
the matrix. Hu et al. use Short and Open together, and it is plausible that the
diagram reduction they report comes mostly from Short. That is the cheapest
remaining test of this route and it was left undone.

## 17. Stage 7a: the theory is confirmed exactly, and the naive route is dead

> **THE SECOND HALF OF THIS SECTION IS SUPERSEDED BY §18.** The path-keyed state
> really does blow up, but that is a property of *that* key, not of the
> determinant side: keying on reachable **labels** costs a flat 1.2×. Read §18
> before concluding anything about the route from this section's numbers.

Measured `benchmarks/decancellation_calibration.py`.

### The calibration passes, and it identifies the mechanism

Song & Shi's four-conductance ladder, with the target re-derived here rather than
taken on trust — from the matrix, and independently as the circuit's **spanning
trees**, which agree with the paper term for term (`G1G2G3`, `G1G2G4`, `G1G3G4`).

| form | terms | `Σ\|term\|` | `\|Σ\|` | `κ` |
|---|---|---|---|---|
| compact-symbol DDD (what we build) | 3 | 6.172e-08 | 1.520e-09 | **40.61** |
| full symbols, no de-cancellation | 9 | 6.172e-08 | 1.520e-09 | **40.61** |
| full symbols, **de-cancelled** | **3** | 1.520e-09 | 1.520e-09 | **1.000000** |

Both gates pass: the determinant is preserved exactly, and `κ` is **1 to ten
decimal places**. So the mechanism §12 attributed to the formulation is confirmed
end to end on a case with a published answer.

Two details worth keeping. The compact form and the un-de-cancelled full-symbol
form have *identical* absolute mass — necessarily so, because with all
conductances positive every composite entry's addends share a sign, so expanding
them redistributes the mass without changing it. That is §15c's inequality
holding with equality, and it is a useful consistency check on the code. And the
surviving terms are exactly the spanning trees, which is *why* they are all
positive: Kirchhoff's expansion has no signs to cancel.

The de-cancellation rule that achieves this is simply **each device symbol at most
once along a term** — a floating device stamps `+g` at (i,i),(j,j) and `-g` at
(i,j),(j,i), so a term using it twice pairs with one of equal magnitude and
opposite permutation sign.

### And the naive implementation is not a representation

The hazard predicted before the run — Tan & Shi's *"cancellation-free s-expanded
DDDs do not satisfy Theorem 1"*, i.e. de-cancellation destroys canonicity — is
now quantified. A de-cancelled expansion cannot key its memo on the minor alone,
because what lies below a minor depends on which devices the path already spent:

| RC ladder `N` | dim | memo states, plain | memo states, de-cancelled | ratio | terms kept |
|---|---|---|---|---|---|
| 3 | 4 | 15 | 39 | 2.6× | 5 / 7 |
| 4 | 5 | 26 | 111 | 4.3× | 13 / 23 |
| 5 | 6 | 42 | 305 | 7.3× | 34 / 76 |
| 6 | 7 | 64 | 820 | 12.8× | 89 / 251 |
| 7 | 8 | 93 | 2177 | **23.4×** | 233 / 829 |

**The sharing loss grows geometrically** — about 1.7-1.8× per added section, with
no sign of levelling. On a 26×26 µA741 that is hopeless: the compact diagram is
1040 vertices, and this construction would need a state count larger by a factor
that is itself exponential in the circuit size. **De-cancelling by a path
predicate is correct and useless.**

### The distinction that matters, and where it points

Notice what the last two columns say together: the **answer** gets *smaller* —
de-cancellation removes 29% of terms at `N=3` rising to 72% at `N=7`, heading for
Tan & Shi's reported 70-90% — while the **search for it** blows up. The
de-cancelled answer is compact; filtering a determinant expansion down to it is
not.

And the kept counts are 5, 13, 34, 89, 233: the spanning-tree counts of a ladder.
**The cancellation-free answer *is* the spanning-tree set.** Which reframes the
whole architecture question: enumerating spanning trees directly is the natural
algorithm for that answer, and filtering a determinant expansion to recover them
is doing exponential work to discard exponentially much of it.

**So this measurement strengthens the case for the topological route, which §15e
had put last.** Two-graph / GPDD methods enumerate exactly the surviving objects.
The reason they exist is now visible from our own numbers rather than taken from
their abstracts.

**What to implement instead, if the determinant side is still preferred:** not a
path predicate but Tan & Shi's construction — a *canceling label list* `CL(L_x)`
per symbol, with a modified coefficient-multiply that removes cancelling terms by
**set operations between sub-diagrams** during the s-expansion. That never forms a
path-dependent state, so it does not pay the ratio above. It is a genuinely
different algorithm from the one measured here, and this stage is the evidence for
why the difference matters.

*Reconsider-if:* the ratio column flattens for circuits with a different topology
— ladders maximise the number of distinct device subsets a path can accumulate.
A denser matrix might share better. Untested, and cheap to test with the same
script.

## 18. Stages 7b-7d: the determinant side is viable after all, and why

§17 concluded that determinant-side de-cancellation was dead. **That conclusion
was about the naive implementation, and it is now overturned for the correct one.**
Three measurements, `benchmarks/decancellation_calibration.py`.

### 7b — it is polynomial, once the state is the reachable remainder

§17's path predicate carried every device used so far. Most of that is dead
information: once a device is consumed at `(p,p)`, its partners at `(p,q)` and
`(q,p)` need row or column `p`, which the minor no longer has, so **at most one
partner label stays reachable** — and it stops mattering as soon as the minor drops
its row or column. Re-keying the memo on the *still-reachable forbidden labels*:

| RC ladder `N` | plain states | naive key | label-canonical key | naive | canonical | terms kept |
|---|---|---|---|---|---|---|
| 3 | 15 | 39 | 18 | 2.6× | **1.20×** | 5 / 7 |
| 5 | 42 | 305 | 52 | 7.3× | **1.24×** | 34 / 76 |
| 7 | 93 | 2 177 | 114 | 23.4× | **1.23×** | 233 / 829 |
| 10 | 232 | 39 546 | 277 | **170.5×** | **1.19×** | 4 181 / 29 867 |

**The overhead is flat at about 1.2×** where the naive key runs to 170×, and the
surviving-term count is identical at every `N` (asserted, not assumed).
De-cancellation removes 29% of terms at `N=3` rising to **86% at `N=10`** — Tan &
Shi's reported 70-90%, reproduced independently.

**A bug worth recording, caught by the gate that was there to catch it.** The first
canonicalisation forbade matrix *positions* and silently lost a term. One position
can carry several devices' addends — entry `(2,2)` of an RC ladder is
`C2*s + 1/R1` — so banning the position to keep `R1` out also banned `C2`. **The
state must be a set of labels, a label being one (device, position) pair.** This is
exactly why the literature says "canceling label list per label", a phrase that
read as pedantry until it cost a wrong answer. Gate 7b-1 required the term count to
match the naive predicate and failed loudly; without it the error would have looked
like a *better* result, because it made the state smaller.

### 7c — de-cancellation alone does not reduce `κ`, and the residue is phase

On RC ladders at a fixed complex `s`, with the determinant verified against
`numpy.linalg.det` at every size (2.4e-15 to 1.1e-15, so the rule preserves the
determinant — including across MNA's voltage-source row, which had been a worry):

| `N` | `κ` compact | `κ` de-cancelled | change |
|---|---|---|---|
| 4 | 3.08 | 1.48 | 2.08× better |
| 6 | 1.95 | 2.17 | 0.90× — worse |
| 9 | 1.11 | 1.97 | **0.56× — worse** |

Two things to take from this. First, **an RC ladder is the wrong test**: its
compact `κ` is already 1-3, so there is nothing to fix, and a method that cannot
improve on 1.1 has not been tested. Second, and more useful: the residual `κ ≈ 2`
is **not sign cancellation at all — it is phase.** At a fixed complex `s` a
conductance term is real and a `C·s` term is imaginary, so their moduli cannot sum
to the modulus of their sum, however cancellation-free the term set is.

### 7d — with s-expansion, `κ` is exactly 1

That diagnosis makes a prediction: separate the powers of `s`, and within one
coefficient of a passive network every surviving spanning-tree product is positive,
so `κ` should be **exactly** 1. Measured, on every coefficient of `rc_ladder(4)`,
`(6)` and `(8)`:

```
worst kappa over coefficients: 1.000000000000   -> PASS (== 1)
```

All eight coefficients at `N=8`, twelve significant figures, no exceptions.

**So the order in Tan & Shi's pipeline is not arbitrary, and this project can now
say why from its own measurements:** s-expansion alone makes `κ` worse (§13),
de-cancellation alone does not help at a fixed frequency (7c), and the two
together give `κ = 1` (7d). §16's "the ingredients are a package" is no longer an
inference from three separate failures — it is a positive demonstration.

### What this changes, and the one thing now blocking it

§17's "the naive route is dead" stands; **its implied conclusion that the
determinant side is dead does not.** Determinant-side de-cancellation is
polynomial (1.2×), preserves the determinant, and combined with s-expansion
achieves `κ = 1` exactly. The topological route is no longer the only one that
gets there.

**The blocker is now infrastructure, not theory.** De-cancellation needs
**one symbol per device**, because a numeric device contribution merges into its
entry and its cancellation becomes invisible. No pycircuit fixture provides that:

| fixture | dim | device symbols | addends | numeric-only addends |
|---|---|---|---|---|
| `rc_ladder(6)` | 7 | 11 | 28 | 2 (the source's ±1 incidence) |
| `cauer_lowpass(3)` | 7 | 9 | 28 | 10 |
| `mfb_filter()` | 4 | 5 | 16 | 2 |
| **`ua741`, all `gm` symbolic** | 26 | 24 | 215 | **159** |

The µA741 bakes 159 of its 215 addends into numbers — `rpi`, `ro`, `cpi`, `cmu`
and every resistor — so de-cancellation there would catch only the `gm` pairs, a
small minority. **The next step is therefore a fully-symbolic µA741 fixture**, not
more algorithm work: `add_small_signal_bjt` already accepts each parameter, so it
is a matter of passing symbols instead of floats and threading them into
`BenchSystem.params`.

*Reconsider-if:* the state count on a 26×26 with ~200 device symbols turns out not
to follow the ladder trend. Ladders are a benign topology for this measurement —
every minor borders few devices. A denser matrix could accumulate more live labels
per minor, and the 1.2× could grow. Cheap to check the moment a fully-symbolic
fixture exists, and it should be checked before anything is built on top.

## 19. Stages 7e-7f: on a real amplifier, `κ` finally lands in a usable range

> **ITS CONCLUSION IS REFUTED BY §20.** The `κ` measurements below stand and were
> verified; the inference drawn from them — that `κ ≈ 15` makes ranking plausible —
> does not. Ranking the de-cancelled expansion needs ≥80× *more* terms than the
> compact one and does not converge, because `κ` says nothing about term counts.
> Read §20 before acting on this section.

The prerequisite §18 identified — one symbol per device — is now built:
`ua741(fully_symbolic=True)` gives every transistor's `gm`, `rpi`, `ro`, `cpi`,
`cmu` and every resistor and capacitor its own symbol. **132 device symbols, 343
additive contributions, of which only 2 are numeric** (the source's ±1 incidence,
which is topology rather than a device) against 168 of 168 numeric before.
Substituting the recorded nominals reproduces the numeric fixture entry for entry
to **2.2e-16**, and two tests pin both properties.

### The declared reconsider-if was right: the ladder overhead does not transfer

| | states | terms |
|---|---|---|
| plain full-symbol | 27 577 | 5.03e+21 |
| label-keyed de-cancelled | **204 679** | 1.10e+21 |

**Overhead 7.42×, against the flat 1.20× ladders gave — gate 7e-1 (≤3×) FAILS.**
Ladders are a benign topology, exactly as flagged: each minor borders few devices,
so few labels stay live. An amplifier is denser and keeps more.

But read the absolute numbers before writing the route off: 204 679 states is
small, it took 14 s, and the overhead is a *constant factor*, not the geometric
growth §17 measured for the path-keyed state. 78% of terms are removed —
consistent with the 70-90% the literature reports. Determinant preserved against
`numpy.linalg.det` to **2.9e-14** (gate 7e-2 PASS), which is the important
correctness statement: the rule is sound on an active circuit with controlled
sources, not only on a passive ladder.

### And the payoff, once compared against the right control

At a fixed complex `s`, `κ` goes **6.66e+03 → 99.05**, a 67× improvement.
Per coefficient of `s` it is much better, but only against a **like-for-like
control** — the same matrix at the same operating point, s-expanded but *not*
de-cancelled:

| `s^k` | compact | de-cancelled | improvement |
|---|---|---|---|
| 0 | 7.63e+05 | 1.26e+04 | 60.5× |
| 1 | 1.85e+03 | **1.73e+01** | 107× |
| 2 | 1.50e+03 | **1.42e+01** | 106× |
| 3 | 1.45e+03 | **1.43e+01** | 101× |
| 12 | 3.70e+04 | 3.31e+02 | 112× |
| 23 | 3.09e+06 | **7.84e+00** | 394 367× |

**Every one of the 24 coefficients improves, by 60× to 394 000×**, and the
low-order coefficients — the ones that carry the response at practical
frequencies — land at `κ ≈ 14-17`.

**That is the first time in this whole thread that `κ` has been in a range where
magnitude ranking is plausible.** At `κ = 15` and `tol = 0.05` the sufficient
condition of §1/§14a asks for `1 - 0.05/15 = 99.67%` of the absolute mass, against
99.99947% at `κ = 9.4e3`. The difference between those two numbers is the
difference between a hard approximation and an impossible one.

`κ = 1` does **not** hold, and was predicted not to: an amplifier is not passive,
so its controlled sources make the surviving terms signed rather than all-positive.
The `s^0` (DC) coefficient stays hard at 1.26e4.

### A mistake I made reading my own numbers, recorded because it nearly stuck

The first reading of the per-coefficient table compared the de-cancelled `κ`
against the **compact whole determinant** (6.66e3) and concluded that
de-cancellation had made the `s^0` coefficient *worse*. That is the wrong control:
the whole determinant and a single coefficient are different quantities. Against
the compact *per-coefficient* value at the same operating point, `s^0` improves
60×. **The comparison was nearly written up as a failure because the control was
convenient rather than matched** — the same error as §14a's, in a different guise,
and the reason the script now computes the control itself instead of reaching for a
number already on the page.

### Where this leaves it

The route works and is affordable. What is *not* built is the thing that would
turn it into an answer: a **rankable diagram**. Everything above is a memoised
recursion that computes values and counts, not a data structure
`approximate_groups` can walk. Turning it into one is Tan & Shi's `CL`/`REMAINDER`
construction, and the measurements now say what it should cost — about 7× the
plain full-symbol diagram on an amplifier, 204 679 states for the µA741.

**The payoff test, once that exists:** rank a de-cancelled `s^1` coefficient at
`κ ≈ 17` and see whether it converges in tens of terms rather than thousands. That
is the experiment this entire thread has been trying to reach, and it is now one
implementation away rather than blocked on a premise.

## 20. Stage 7g: the payoff test FAILS, and it corrects §19

The experiment this whole plan was built toward. `benchmarks/decancellation_ranking.py`,
µA741 at 1 kHz, `tol = 0.05`, control measured at the same operating point rather
than quoted:

| | terms | error | converged |
|---|---|---|---|
| compact, group ranking | **870** | 3.75e-02 | **yes** |
| compact, magnitude ranking | **1 690** | 3.78e-02 | **yes** |
| de-cancelled, group ranking | 72 724 | 1.40e-01 | **no** (hit the 400 000-split cap) |
| de-cancelled, magnitude ranking | 68 310 | 1.59e-01 | **no** |

**Gate 7g-2 FAILS, and not narrowly: the de-cancelled expansion needs at least 80×
*more* terms and still does not converge**, despite `κ` being 67× better (99.05
against 6 659).

### Why, and this is the correction to §19

§19 concluded that `κ ≈ 15` per coefficient was "the first time in this thread that
`κ` has been in a range where magnitude ranking is plausible". **That inference was
wrong, and this measurement is what refutes it.**

`κ = Σ|term| / |Σ term|` is a statement about the *aggregate*. The sufficient
condition it supports — capture `1 − tol/κ` of the absolute mass — says how much
mass is needed and **nothing whatever about how many terms that takes.** Those are
independent:

* the compact form packs the determinant into **2.77e+06** terms;
* the de-cancelled full-symbol form spreads the same value over **1.1e+21** — fifteen
  orders of magnitude more, each correspondingly smaller.

So de-cancellation improves the *conditioning* and destroys the *concentration*.
Needing 99.95% of the mass is easy when 870 terms carry it and hopeless when it is
diluted across 10^21. **Low `κ` is necessary, not sufficient.**

That is the third time in this thread a claim about `κ` has had to be narrowed —
first sufficiency stated as necessity (§14a), then scope (§15b), now this — and the
pattern is the finding: **`κ` bounds the error achievable from a given fraction of
the mass; it predicts nothing about term counts, which is what "readable" actually
depends on.**

### What the thread therefore concludes

For *this* representation and this circuit, the compact-symbol diagram with group
ranking remains the best route to a converged approximation: **870 terms at 3.7%**,
and it is shipped and tested. De-cancellation is a real and correct transformation
with a measured, affordable cost, and it is the wrong lever for readability, because
readability is a term-count property and de-cancellation moves term count in the
wrong direction by fifteen orders of magnitude.

**The missing diagnostic, which this thread should have had from the start:** a
*concentration* measure alongside `κ` — how many terms carry (say) 99% of `Σ|term|`.
That is the quantity a ranking's term count depends on, it is as cheap to compute on
a diagram as `κ` is, and had it been measured in §1 the de-cancellation route would
have been assessed differently. It is the obvious next thing to build.

**Scope, stated:** this was measured at fixed `s`, where the de-cancelled `κ` is 99
rather than the 14-17 of the low-order coefficients, and the gate said in advance
that a failure here would not condemn the per-coefficient route. It does not — but
the *dilution* argument applies there too and more strongly, since the coefficient
expansions have astronomically many terms as well. Per-coefficient ranking was not
measured, and that is the one remaining way this conclusion could be wrong.

### A bug the gate caught, worth keeping

The first run reported a magnitude-ranking error of **6.7e+55**. That is
*impossible*: with `Σ|term| = 99·|det|`, no partial sum can be more than about 100
relative units away. Reasoning from that bound rather than accepting the number
found the cause: a state with rows still remaining but **every child pruned** is a
*dead end* contributing zero, and the code was treating it as a completed term and
adding its whole prefix product. Verified after the fix by enumerating a small
circuit both ways — 13 paths, 13 brute-force terms, both rankings agreeing exactly.

**The general lesson: know the bound your number must satisfy before you read it.**
Gate 7g-4 existed because a ranking that converges to the wrong value is not a
result, and it is the only reason this was caught rather than written up as a
property of the method.

## 21. Stage 8: the missing diagnostic, and it predicts what `κ` could not

§20 identified the gap: `κ` measures conditioning, term count depends on
concentration, and only the first was being measured. This builds the second.

**`DDD.concentration(env)`** returns the participation ratio

    N_eff = A[root]**2 / S2[root]     with  S2[v] = Σ |term|**2

computed in **one traversal**, the same cost as `cancellation`, because `|Π|²` is the
product of the squares. It is the *effective number of terms*: equal to the true term
count when all terms have the same magnitude, falling to 1 when one dominates.

### It passes the discriminating test

| | terms | `κ` | `N_eff` | terms a ranking needed |
|---|---|---|---|---|
| compact-symbol | 2.77e+06 | 6 659 | **194** | **870**, converged |
| de-cancelled full-symbol | 1.10e+21 | **99** | **11 565** | >72 724, did not converge |

**`N_eff` orders the two representations the way their ranking costs actually came
out; `κ` orders them backwards.** And it does not merely order them — it predicts
the magnitude: 870 against `N_eff = 194`, and >72 724 against 11 565, both within a
factor of about five.

Gates 8-1 and 8-2 pass too: `N_eff` equals the term count exactly for a diagram of
equal-magnitude terms, falls to 1.0 for a dominated one, stays within `[1, N]`, and
tracks the *enumerated* count of terms needed for 99% of the absolute mass to within
an order of magnitude on a circuit small enough to check. Five tests, including one
that pins the independence directly: a 2×2 matrix can have `κ > 10^8` and
`N_eff = 2`.

### What this means for the thread

The two diagnostics are now both cheap, both shipped, and they say different things:

* **`cancellation`** — can a truncation of this expansion be accurate at all, and
  what fraction of the absolute mass must it capture?
* **`concentration`** — how many terms is that likely to take?

**A ranking is cheap only when both are small.** Reporting one without the other is
what made the de-cancellation route look promising for four stages: its `κ` was 67×
better while its `N_eff` was 60× worse, and only the first was on the page.

**The retrospective judgement, stated plainly:** had `concentration` existed at §1 it
would have been visible immediately that the compact-symbol µA741 has `N_eff = 194`
— i.e. that a few hundred terms should suffice — which is exactly what group ranking
then achieved. The whole excursion through full symbols, de-cancellation and
spanning-tree theory (§§13-20) was chasing a quantity that was never the binding
one. That excursion produced real results, several of them published-and-confirmed,
and it was avoidable.

*Reconsider-if:* `N_eff` is a second moment, so it is blind to distribution shape
beyond the first two moments — a heavy tail with the same variance reports the same
number. If a case appears where `N_eff` is small and ranking is still expensive, the
next refinement is the enumerated mass profile itself (terms to reach 90/99/99.9% of
`A`), which `iter_terms` already provides at the cost of running the ranking.

## 22. Stage 9: per-coefficient ranking helps 5-12×, and `N_eff` has a resolution limit

Stage 8 made the readability question answerable: the compact µA741 has
`N_eff ≈ 151` and group ranking needs 734 terms, so a *readable* answer needs
`N_eff` of order ten. Stage 9 asks which available lever lowers it. Measured
`benchmarks/concentration_profile.py`.

### The result: a real gain, short of the gate

| | whole determinant | best coefficient | gain |
|---|---|---|---|
| nominal | 734 terms, `N_eff` 151 | **`s^1`: 155 terms**, `N_eff` 74.5, carries **97.3%** of the response | **4.7×** |
| degraded | 2 401 terms, `N_eff` 133 | **`s^2`: 206 terms**, `N_eff` 709 | **11.7×** |

Device symbols (`gm_q1`, `gm_q2`, `gm_q17`) survive in every case, and the
contribution share is printed beside each coefficient so a concentrated-but-
irrelevant one cannot be claimed as a win.

**Gate 9-2 FAILS** — the target was `N_eff ≤ 20` and ≤30 terms, and the best is 74.5
and 155. **Gate 9-3, PARTIAL, is what this is:** per-coefficient ranking is a genuine
5-12× reduction with symbols intact, and it is the first lever measured that moves
term count in the right direction at all. 155 terms is still not readable.

### And `N_eff` does not survive the run unqualified

Look at the degraded row by row: `N_eff` of 126.8, 341.2, 709.2 for `s^0`, `s^1`,
`s^2` — and the ranking costs are **950, 275, 206 terms**. That is *inverted*: the
coefficient with the largest `N_eff` was the cheapest to rank.

So the enthusiasm of §21 needs a bound put on it. `N_eff` correctly ordered the two
*representations* in stage 8, where they differed 60-fold (194 against 11 565) and
their ranking costs differed 80-fold. **It does not reliably order quantities whose
`N_eff` differs by only a factor of a few.** It is a coarse screen — good for "is
this representation hopeless?", not for "which of these coefficients should I rank?".

That is exactly the limitation §21 declared in advance — a second moment is blind to
tail shape — arriving sooner than expected. The refinement named there, the
enumerated mass profile, is now the thing to build if this ordering matters.

### A bug in the diagnostic, caught by a bound

The first run reported `N_eff = 0.0` for every coefficient above `s^12`. **Impossible:
`N_eff ≥ 1` by construction.** Cause: accumulating `S2 = Σ|term|²` directly underflows
— a product of 26 squared entries of order `1e-12` is `1e-624`, which is zero in
double precision long before the *ratio* is. Rewritten scale-free, carrying the
inverse participation ratio `r = S2/A²` on the weights:

    w1 = e*A[one]/A[v],  w0 = A[zero]/A[v],   r[v] = w1²*r[one] + w0²*r[zero]

every quantity bounded in `(0, 1]`. The low-order values are unchanged by the fix
(117.3, 74.5, 187.2 before and after), which is the check that it altered only the
underflowing cases. A test now pins `1 ≤ N_eff ≤ term_count` on the high-order
coefficients specifically.

**Second time in three stages that a known bound on a quantity caught a bug its
value alone would not have.** The other was the 6.7e+55 error in §20. Knowing what
range a number must lie in is worth as much as computing it.

## 23. Stage 10: an eleven-operation expression, and the metric was wrong all along

> **ITS GENERALISATION IS WITHDRAWN BY §24.** The eleven-operation measurement stands.
> The claim that term count overstates the answer 4-8× *in general* does not: that
> ratio holds for an s-expanded coefficient diagram and inverts for the compact
> whole-determinant one, where 734 groups are 50 377 operations.

Every stage in this plan used `tol = 0.05`, inherited from the first round and never
questioned. Stage 10 varied it. `benchmarks/tolerance_curve.py`, µA741's dominant
coefficient (`s^1`, which carries 97.3% of the response at 1 kHz, 51 054 561 terms):

| `tol` | groups kept | **operations** | achieved error |
|---|---|---|---|
| 0.5 | 5 | **11** | 26.3% |
| 0.2 | 43 | 57 | 6.0% |
| 0.05 | 155 | 91 | 3.2% |
| 0.01 | 164 | 91 | **0.79%** |

Degraded `gm`: 5 / 11 / 24.9%, then 103 / 67 / 7.5%, 275 / 86 / 2.6%, 577 / 108 /
0.002%.

### The expression

At `tol = 0.5`, five kept groups collect into a single product:

```
5.9996902350538e-70 * gm_q17 * (gm_q1 + 4.09e-4) * (-gm_q17 - 1.0022501e-2)
                             * (gm_q2 + 4.09e-4)
```

**Eleven operations, four factors, three device symbols, 26% error.** And it reads
as a circuit statement: the coefficient goes as the gain stage's `gm_q17`, times each
input transistor's transconductance in parallel with a conductance, times the gain
stage's own loading. That is the shape a designer would write by hand.

**Gate 10-2 FAILS as declared** — 26.3% exceeds the 20% threshold, and the setting
that does meet 20% costs 43 groups, above the 30 allowed. The gate falls in the gap
between two rows of the table. Recording it as a fail because that is what it is;
the 20% figure was my own choice and the substance is not in doubt.

### The metric was wrong, and that is the transferable finding

**Term count is a poor proxy for readability.** Five groups collect to 11
operations; 155 groups collect to 91. Across the whole curve the terms grow 33-fold
(5 → 164) while the operations grow 8-fold (11 → 91), because sympy collects shared
factors that the term count counts separately.

This plan has measured term counts throughout — 870, 734, 155, 72 724 — and
**operations are what a reader actually faces.** By that measure the position is far
better than any previous section suggested: **91 operations for 0.79% error** on a
real amplifier's dominant coefficient, with device parameters intact.

Every earlier "not readable" verdict in this document was reached on term counts. It
was the wrong yardstick, and §§16, 19, 20, 22 should be read with that in mind: their
*measurements* stand, their readability judgements were made against a quantity that
overstates the size of the answer by roughly 4-8×.

### What is and is not achieved

**Achieved:** for one coefficient of one amplifier, a symbolic expression in device
parameters, readable by inspection, with a stated error, and *different at different
operating points* — which is the whole of the maintainer's original brief applied to
a coefficient.

**Not achieved:** the same for the *complete* transfer function, which needs every
coefficient (24 of them) and a numerator; and the leapfrog, where stage 5 showed the
composition is not error-controlled. The honest summary is that the goal is met on a
component of the problem and not on the whole of it.

*Reconsider-if:* the `s^1` coefficient dominating at 1 kHz is what makes one
coefficient nearly sufficient here. At a frequency where two or three coefficients
share the response, the readable object is their combination, and nothing measured
here says that stays at tens of operations.

## 24. Stage 11: a transfer function at 177 operations — and §23 over-generalised

### Part 1 corrects §23, and the correction matters

§23 concluded that "term count overstates the answer's size 4-8×" and that every
earlier readability verdict was therefore reached on the wrong yardstick. **Measured
on the whole determinant, that is false and the error goes the other way:**

| object, at `tol` | groups | **operations** | ratio |
|---|---|---|---|
| `s^1` coefficient, 0.5 | 5 | **11** | 2.2 ops/group |
| whole determinant, 0.5 | 5 | **338** | 68 ops/group |
| whole determinant, 0.05 | 734 | **50 377** | 69 ops/group |
| leapfrog top diagram, 1e-3 | 181 | **2 895** | 16 ops/group |

**So the terms-to-operations ratio is a property of the diagram, not a constant.**
In an s-expanded *coefficient* diagram each vertex payload is a bare number or
symbol, so collection is effective and five groups really do become eleven
operations. In the *compact whole-determinant* diagram each payload is a polynomial
like `gm_q1 + 6e-12*s + 4.09e-4`, so every group is a product of 26 binomials and the
operation count is ~69× the group count.

**§23's generalisation is withdrawn.** Its specific measurement stands — the `s^1`
coefficient really is 11 operations — but "every earlier not-readable verdict was
reached on the wrong yardstick" is wrong. For the whole determinant those verdicts
were, if anything, **too kind**: 734 groups is 50 377 operations. The leapfrog's 181
groups are 2 895 operations over placeholders, so that verdict is unchanged too.

That is the fourth self-correction in this thread and they share a shape:
**generalising from one measurement to a class.** §14a (sufficiency as necessity),
§15b (scope), §22 (`N_eff` as a fine predictor), and now this. The measurement was
right each time; the sentence around it reached further than the measurement did.

### Part 2: the deliverable, and it passes at one operating point

`H(s) = N(s)/D(s)` with the numerator being `A` with the output column replaced by
`b`. Both s-expanded, coefficients chosen greedily with **exact re-evaluation of `H`
over the whole sweep after every drop** — never a per-coefficient budget. Verified
first that the exact `N/D` reproduces the linear solve across the sweep:
**3.4e-15** at nominal, **1.4e-14** at degraded.

Sweep: 13 points, 100 Hz to 10 kHz, two decades, sampled and stated as sampled.

| point | `tol` | kept | **operations** | sweep error | gate |
|---|---|---|---|---|---|
| nominal | 0.2 | N `[0]`, D `[1]` | 119 | 26.7% | fail (error) |
| nominal | 0.05 | N `[0,1]`, D `[0,1,2]` | 439 | 6.5% | fail (size) |
| **degraded** | **0.2** | **N `[0]`, D `[0,1]`** | **177** | **7.9%** | **PASS** |
| degraded | 0.05 | N `[0,1]`, D `[0,1,2]` | 362 | 1.3% | fail (size) |

**Gate 11-2 passes at the degraded operating point and fails at nominal** — 177
operations holding 7.9% relative error across two decades, device symbols intact.
That is the first *transfer function* this project has produced at a readable size
with an error verified over a band rather than at a point.

At nominal the same tolerance keeps a different coefficient set and lands at 26.7%,
which is over the gate; buying 6.5% costs 439 operations. So the honest verdict is
**PARTIAL: it works at one of the two operating points and not the other.**

**Worth noting for the original brief:** the kept coefficient sets *differ* between
operating points — `D [1]` at nominal against `D [0,1]` at degraded. That is
"different expressions for different symbol values", which the maintainer asked for
and which no amount of numeric substitution provides.

### Where this leaves the plan

The gap between 177 and 439 operations is the gap between 8% and 1.3% error. Nothing
here is blocked; what remains is a size/accuracy trade that a designer would choose
rather than a missing capability. The unmeasured levers are unchanged — the `Short`
operation, and the topological route — and neither is needed for what this stage
produced.

*Reconsider-if:* the sweep is 100 Hz-10 kHz. The µA741's compensation pole and the
unity-gain frequency lie above it, so a wider sweep would bring more coefficients
into play and the operation count with them. Widening it is one line and the honest
next check on this result.
