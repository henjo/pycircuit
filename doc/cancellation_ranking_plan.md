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
