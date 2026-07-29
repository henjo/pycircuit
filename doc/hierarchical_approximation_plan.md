# Term-ranked approximation on hierarchical diagrams — implementation plan

**Status: planned, not started.**

**Goal, stated by the maintainer:** *full symbolic extraction of a 5th-order
leapfrog filter built from five µA741 amplifiers, with no numeric substitution
of the circuit itself* — then **approximation** to reduce the result to
expressions a designer can read, with **different expressions for different
symbol values**.

The distinction that matters, and that an earlier round of this work blurred:
numerics are used to **rank** terms, not to replace them. Symbols survive into
the answer, including `s`, so poles and zeros remain readable. Substituting
values into the circuit would make the exercise pointless — "if we replace all
symbols with numbers I can just run SPICE".

## 1. What is already measured

All on the circuits themselves, no estimates.

**The target circuit exists and is correct.** `leapfrog_5th_order()` —
127 unknowns, 536 nonzeros, five full µA741s with leapfrog feedback. Verified
as a filter, not merely as a matrix: flat `-6.02 dB` passband (the 0.5 of a
doubly-terminated ladder) and `-100 dB/decade` in the stopband, which is
5 x 20. An earlier version omitted the damped terminating integrators, had no
passband at all, and would have been a convincing-looking fixture for
something that was not a fifth-order filter.

**Flat diagrams do not reach it.** Fully symbolic (121 symbols, `s` symbolic),
`ddd_of_matrix` on the leapfrog was killed at **15 min 16 s and 2.7 GB**,
unfinished.

**Hierarchical diagrams do.** Suppressing the 111 internal nodes and keeping
the 16 terminals — each amplifier exposes only `inp`, `inn`, `out` —
gives **299 vertices in 2.2 s**, verified against `numpy` at three
frequencies: **3.9e-15, 1.4e-15, 6.6e-15**.

**Approximation is exactly the wanted workflow, and works on flat diagrams.**
On the µA741 (1040 vertices, 2 773 885 terms):

| operating point | terms kept | secs | symbols in the result |
|---|---|---|---|
| nominal `gm = 1.0 mA/V` | 260 | 0.15 | `gm_q1, gm_q2, gm_q17, s` |
| degraded `gm = 0.1 mA/V` | **180** | 0.08 | `gm_q1, gm_q2, gm_q17, s` |
| boosted `gm = 5.0 mA/V` | 260 | 0.08 | `gm_q1, gm_q2, gm_q17, s` |

2.77M terms to a few hundred, symbols intact, and **different operating points
genuinely give different expressions** — the degraded point needs 30% fewer
terms because different products fall below tolerance.

**The gap is one capability on one class.**

| | flat `DDD` | `HierarchicalDDD` |
|---|---|---|
| `eval` | yes | yes |
| `size` | yes | yes |
| `to_sympy` | yes (refuses above a cap) | no |
| `term_count` / `iter_terms` | yes | no |
| **`approximate`** | **yes** | **no** |

So: the only representation that reaches the target circuit is the only one
that cannot approximate.

## 2. Why this is not a missing accessor

Worth stating before any stage, because it sets what the stages have to solve.

In a **flat** diagram every root-to-terminal path *is* a product term, its
magnitude is the product of its entries' magnitudes, and ranking is a
branch-and-bound over paths. `DDD.approximate` already does this.

In a **hierarchical** diagram the entries of the reduced matrix at level `k+1`
are themselves *rational functions* of level `k` — Schur complements,
`A22 - A21 A11^-1 A12`. A "term" therefore spans levels, and a magnitude must
be propagated through a division and a subtraction at every level.

Two consequences that are genuine hazards, not difficulties of implementation:

- **Subtraction destroys the ranking premise.** Term-ranked approximation
  assumes small terms may be dropped. A Schur complement subtracts, so the
  *result* can be far smaller than its parts, and a term that is negligible
  against its neighbours can be decisive against the difference. Any bound
  that ignores this is unsound rather than merely loose.
- **Error at a low level is amplified by division.** Dropping a term inside
  `A11` perturbs `A11^-1`, and near a pole of the reduced network that
  perturbation is unbounded.

**Reconsider-if:** if stage B cannot produce a sound bound, the honest outcome
is a *heuristic* ranking with a measured error distribution rather than a
guarantee, clearly labelled as such — which is still useful, and is what most
published symbolic approximation actually offers.

## 3. Stages

Ordered cheapest-first, and the first stage may make the rest unnecessary.

### Stage A — approximate the *reduced* network, treating its entries as atoms

Suppression already produces a 16x16 reduced matrix whose entries are rational
functions of the suppressed nodes. **Build a flat diagram of that**, and
approximate it with each reduced entry treated as an opaque atom.

This needs **no new theory**: the flat machinery applies unchanged, and the
result answers a real design question — *which paths through the
sixteen-terminal network matter* — even though it cannot see inside an
amplifier.

**Gate.** A simplified symbolic expression for the leapfrog's transfer
function at a stated operating point, agreeing with the full hierarchical
evaluation to within the requested tolerance, in under a minute. And two
different operating points must produce **different** term sets, or the
exercise has not demonstrated what it is for.

**If this suffices, stop.** It is the cheapest thing that could work, and the
plan says so in advance so that reaching for the harder stages is a decision
rather than momentum.

### Stage B — magnitude propagation through one Schur level

The core primitive: given magnitudes at level `k`, bound the magnitude of a
term of the reduced matrix at level `k+1`.

**Gate — the flat diagram is the oracle.** On a circuit small enough that both
routes finish (the µA741 is the obvious one), the hierarchical ranking must
select the **same dominant terms as flat ranking**, or a superset, at the same
tolerance. A subset is a failure: it means terms are being dropped that the
flat method keeps.

Record the subtraction hazard explicitly here: construct a case with
near-cancellation in the Schur complement and check the bound still holds. If
it does not, that is the negative result and stage C should not be attempted
on a guarantee it cannot honour.

### Stage C — `approximate` on `HierarchicalDDD`

Only if B produces a sound bound.

**Gate.** On the µA741, where both routes are feasible, hierarchical
approximation must reproduce flat approximation's *value* to the stated
tolerance with a comparable term count. Not identical term sets — the
groupings differ — but the same answer to the same accuracy.

### Stage D — the target

The 5th-order leapfrog, fully symbolic, approximated at a stated operating
point.

**Gate.** A symbolic expression small enough to read — fewer than ~1000 terms
— verified against the numeric solve to a stated tolerance, and **visibly
different** at a second operating point. Report the term count and the error
together; either alone is misleading.

### Stage E — is it worth it?

**Declared in advance.** Stage C is worth keeping over stage A only if it
reaches inside the amplifiers in a way a designer can act on — for example
naming a specific transistor's `gm` as dominant — rather than merely producing
a smaller number of larger terms.

If it does not, **record it and ship stage A**, which is already useful.

## 4. What would make this fail

- **The nested form may not admit meaningful ranking at all.** A single term at
  the top level can stand for an enormous set of flat terms, so "keep the top
  200" may mean something quite different at each level.
- **Near-cancellation**, per section 2. This is the one most likely to turn a
  guarantee into a heuristic.
- **The reduced entries may already be too large to read**, in which case
  stage A produces a compact *structure* over unreadable *atoms* and the
  designer is no better off. Measure the size of a reduced entry before
  assuming otherwise.

## 5. Outcomes

_To be filled in as stages complete. Negative results recorded here, not
edited out — the DDD work kept five and they are among its more useful
output._

| Stage | Outcome |
|---|---|
| A — approximate the reduced network | not started |
| B — magnitude propagation | not started |
| C — hierarchical `approximate` | not started |
| D — the leapfrog | not started |
| E — worth it? | not started |

## 6. Facts carried in

- **`HierarchicalDDD` exposes only `eval` and `size`.** Everything else in
  section 1's table is flat-only.
- **`H.eval()` with no environment does return a symbolic expression** — so
  the hierarchical form is not evaluation-only. On the µA741 with three
  symbolic transconductances it unfolds to `count_ops = 2 256 398` in 22.8 s,
  which is the same order as the flat form's 2.77M terms. **The compression is
  in the representation, not in the answer** — which is precisely why
  approximation, not unfolding, is the thing worth building.
- **`DDD.to_sympy` refuses above a cap** (`DDDSizeError`) rather than
  attempting a 2.77M-term expansion. Any new extraction should fail the same
  way rather than hang.
- **`suppression_order(A, keep=...)`** already chooses the elimination order;
  the leapfrog's terminals are the three per amplifier plus the input.
- **The expression graph in `distortion_ddd.py` is not an alternative here.**
  It is a straight-line program of the elimination's arithmetic: it evaluates
  fast and cannot rank terms, extract coefficients of `s`, or approximate.
  That distinction is documented on the `distortion_ddd` page and is the
  reason this plan exists.
