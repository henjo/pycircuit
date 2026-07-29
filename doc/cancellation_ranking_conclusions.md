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
