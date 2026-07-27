# DDD for pycircuit — findings, reasoning, conclusions and priorities

Companion to `ddd_references.md`. That file says *what to read and why*; this one
records *what we concluded and why*, so that the implementation plan can be
written from the two together without re-deriving the argument.

Status: research complete, no DDD code written yet. Date: 2026-07-27.

**How to read the claims below.** They are tagged, because the confidence differs
a lot and that matters for planning:

| Tag | Meaning |
|---|---|
| **[LIT]** | Claimed in a paper. Credible, but measured on *their* circuits and implementation. |
| **[OURS]** | Measured by us, in this repo, and reproducible from a doc page or test. |
| **[VERIFY]** | An assumption we are carrying into the plan that must be measured before it is relied on. |

---

## 1. The question

pycircuit's symbolic analysis is fast where the circuit is *mostly numeric with a
symbolic `s`*, and falls over where the circuit is *fully or largely symbolic*.
The roadmap item is whether a Determinant Decision Diagram fixes the second case,
and if so how to add it **without disturbing the existing toolkits** — which is a
hard requirement, since `numeric`, `symbolic` and `symbolic_poly` all have
downstream users.

---

## 2. Starting point — what is actually in the repo

Verified against the tree at commit `5bae9b7`, not from memory:

- **`pycircuit/circuit/ddd.py` (211 lines) — a partial, orphaned implementation.**
  `Node` exists and is reasonable: `cofactor`, `remainder`, `union`, `intersec`,
  `__mul__`, `__sub__`, `eval`, `asdot`. What is missing or broken:
  - `DDD_of_matrix` is **non-functional** — it calls a 3-argument `DDD(...)` while
    `class DDD` is an empty `pass`; it tests `A == sympy.zero(1)` / `sympy.one(1)`
    which do not exist in modern sympy; it calls `A.minorMatrix` (renamed
    `minor_submatrix`); and its docstring example uses invalid `sympy.Matrix`
    syntax. It has clearly never run.
  - **There is no unique-node hash table.** `__eq__` is a recursive structural
    comparison, so nothing is *shared* — which is the entire point of a DDD.
    This is the single most important gap.
  - `asdot` needs `yapgvb` (import already commented out) — a py2-era dead end.
- **`pycircuit/circuit/tests/test_ddd.py` exists** and mostly passes: three
  `Node`-level tests plus one `@unittest.skip`ped test that references a `.top`
  attribute the class never had (it is `.index`). Its second test class is taken
  from *Symbolic Analysis and Reduction of VLSI Circuits* (Tan & He, Springer
  2005) p.198 — worth adding to the reading list as the DDD-era book.
- **Toolkit architecture is ready for this.** `Toolkit` base with a
  `supports(capability)` hook returning `False`;
  `SymbolicPolyToolkit.supports` returns `capability in ('num_den',)`;
  `linearsolver_num_den(self, A, b)` is the extension point (two arguments —
  `GinacToolkit` detects the frequency symbol internally via
  `_ginac.detect_freq`). Singletons `numeric` / `symbolic` / `symbolic_poly`, with
  `ginac_toolkit` / `symengine_toolkit` optional.
- Existing result objects to be a sibling of: `TransferFunction` /
  `CircuitResultACPoly`, `GinacResult`, `SoESolution` (`soe.py`).

**Conclusion 2.1** — the existing `ddd.py` is worth keeping only as the `Node`
vocabulary and its tests. `DDD_of_matrix` should be *rewritten*, not repaired, and
the rewrite must introduce a hash table for sharing.

---

## 3. What we already learned the hard way (prior negative results)

These are ours, and they constrain the plan more than the literature does.

- **[OURS] symengine was not a win.** Its `LUsolve`/`det` are not fraction-free,
  so nested divisions swell; `FFLU` does not cancel (last pivot ~125M chars at
  N=8). End-to-end slower than `symbolic_poly` (N=10: ~1.4 s vs ~0.04 s).
- **[OURS] GiNaC was not a win either, for a different reason.** `normal()` does
  cancel properly, but numeric component values (1e-9) become rationals whose
  products across a determinant explode into ~(1e9)^N integers that stall CLN —
  hanging around dim 16. Fixed only by frequency/impedance scaling to O(1)
  coefficients using **powers of ten** (irrational scale factors reintroduce the
  blow-up via `limit_denominator`).
- **[OURS] For *fully symbolic* circuits, `symbolic_poly` has no advantage at
  all** — the multivariate determinant is inherently exponential (spanning-tree
  terms). Fraction-free elimination only wins with *few generators*.
- **[OURS] Sharing beats inlining, decisively.** `soe.py`'s Gaussian elimination
  as a sequence of assignments stays linear on a fully-symbolic ladder — 73 / 157
  / 241 / 325 ops at N = 4 / 8 / 12 / 16 — while inlining it to one expression
  (`to_ratio`) regrows super-linearly: 213 / 1235 / 3633 / 7983. A single GiNaC
  `compile_ex` does not help, because it has no CSE and re-expands.

**Conclusion 3.1** — the failure mode to design against is not "slow arithmetic",
it is **expression growth**. Two of three previous attempts died on it. Any DDD
work must be justified on *representation size*, and a faster arithmetic backend
is a second-order concern.

**Conclusion 3.2** — we already have one representation that avoids the blow-up
(SoE). So DDD has to clear a real bar, not merely work. See §7 gate.

---

## 4. Findings from the literature

### 4.1 Construction: implement the 2010 algorithm, not the 2000 one

**[LIT]** Shi's ICCAD 2010 "simple implementation" replaces the TCAD-2000
five-step flow and deletes three subsystems we would otherwise have to build:

| TCAD 2000 flow | ICCAD 2010 (LED) |
|---|---|
| Greedy-Labeling symbol-ordering pre-pass | none — on-the-fly *expansion* ordering (min-degree row/col) |
| ZBDD package for vertex-triple sharing | a dict hashing minors by `(row-index tuple, col-index tuple)` |
| separate sign-determination pass over 1-paths | signs fall out during construction |

Mechanism: expand a whole selected row/column at a time into layered queues
(Layered Expansion Diagram), sharing minors through the hash table; then convert
LED→DDD mechanically (existing next-layer pointers become 1-edges; add 0-edges
along each sibling group ending at `0`; bottom queue terminates at `1`).
`|DDD|` = total elements across queues. Property 4: order *within* a sibling
group provably does not change the size.

**Conclusion 4.1** — this is the construction to implement. It is a far better
fit for pure Python than the BDD-package route, and it is the difference between
"port a BDD library" and "write a dict-keyed recursive expansion".

### 4.2 The hard limit, and where the win actually comes from

**[LIT]** Shi, TCAS-II 2010: a row-wise order is *optimal* for a **full** n×n
matrix, and the optimal size is then exactly `n·2^(n-1)`.

**Conclusion 4.2** — DDD is still exponential in the worst case. All of its value
comes from **sparsity and sharing**, so the Phase-0 measurement must be vertex
count on *our* sparse circuit matrices, and `n·2^(n-1)` is the "no sharing
achieved" reference line to plot against. This also settles ordering for the dense
case, so we need not go hunting for a heuristic.

### 4.3 The key structural finding: the coefficient / s-expanded form

This is the most important result of the review, and it arrives twice:

- **[LIT]** *Semi-symbolic* — Pi & Shi, MTDDD (DAC 2000). Terminals are extended
  beyond `{0,1}` to carry arbitrary **numeric values**, so numeric sub-products
  collapse into a single terminal *during construction*. The "coefficient MTDDD"
  is **multi-root**: one root per power of `s`.
- **[LIT]** *General and proven* — Shi & Tan, s-expanded DDD (TCAD 2001).
  "s-expanded DDD" and "multiroot DDD" are the same object. **Theorem 1**: from a
  complex DDD of size `|DDD|`, the s-expanded DDD is built in time and size
  ~`q·d·|DDD|`, where `q` = denominator degree and `d` = max devices per node
  (`d ≤ 2` under MNA compact-symbol form). So **s-expansion is linear in the
  complex DDD**, though the expanded term count is astronomically larger — µA741:
  108 032 complex product terms → ~7.8×10³⁴ expanded terms, held in 99 844
  vertices, built in seconds on a 1998 workstation.

Why this matters more than it first appears: **the s-expanded coefficient form is
the shape pycircuit already speaks.** `TransferFunction.as_num_den` and
`_ginac.poly_coeffs` produce exactly `N(s)/D(s)` coefficient vectors, and
`poles()` / `zeros()` / bode consume them. A multiroot DDD drops into that
interface rather than requiring a new one.

**[LIT]** It also answers a problem we hit independently: MTDDD's numeric
terminals fold numeric values into a single number as they arise, so the exact
huge-rational blow-up that capped the GiNaC backend at ~dim 16 (§3) **never
forms**. Different mechanism, same disease.

**[LIT]** And the noise application is our own `noise_psd` insight, generalised:
noise analysis is a *set* of transfer functions, one per source, sharing most
subexpressions, so represent them all in a **single multiroot DDD** at a cost
comparable to one transfer function.

**Conclusion 4.3** — the target is not "a DDD"; it is a **multiroot s-expanded
DDD**. Plain flat-determinant DDD is a stepping stone to it, not the deliverable.

### 4.4 Approximation depends on s-expansion — an ordering constraint

**[LIT]** TCAD 2001 §I: approximation *requires* the s-expanded form so that the
coefficient of each power of `s` is approximated on an equal footing; otherwise
the result is unreliable. Dominant poles/zeros then come out as ratios of
coefficients of consecutive powers of `s`.

**Conclusion 4.4** — `approximate(tol)` cannot be an independent early feature.
It sits *after* s-expansion. The previous reading list had these as parallel
options; that was wrong.

### 4.5 Hierarchy

**[LIT]** Tan & Shi (TCAD 19(4), April 2000) is the primary paper: suppress each
subcircuit to its terminals via its determinant and cofactors, then Cramer at the
top. On a µA741, three-level two-way hierarchy needs **117 vertices vs 6654 flat
(56×) vs 119 011 sum-of-product terms**, and DDD size grows almost linearly in
circuit size while product terms grow exponentially.

**[LIT]** Xu, Shi & Li (ASP-DAC 2011) is the better *fit*: a subcircuit reduces to
a **"symbolic stamp"** — an n-terminal block that stamps into the parent MNA
matrix exactly as a primitive element does, which is how pycircuit's `SubCircuit`
already composes.

**Conclusion 4.5** — hierarchy is real leverage but it is a *later* stage: it
multiplies whatever the flat engine achieves. Take the symbolic-stamp
formulation when we get there.

### 4.6 A published DDD-vs-SoE head-to-head — directly relevant to us

The TCAD-2000 hierarchical paper benchmarks against **SCAPP, a
sequence-of-expressions analyser** — the same representation as our `soe.py`.

**[LIT]** Their claims: DDD is "much more compact than the sequence-of-expression
representation used in SCAPP", and beats both SCAPP and SPICE on repetitive
evaluation.

Critically, they give us a **comparable metric**: each DDD vertex costs one
addition and one multiplication, so `|DDD|` *is* an operation count — the same
quantity we already report for SoE.

**[VERIFY]** We should not accept "DDD ≪ SoE" as given. Their SCAPP baseline is
hierarchical-suppression SoE on cascaded opamps; ours is Gaussian-elimination SoE
measured on ladders, where **[OURS]** we found linear growth (§3). These may not
be the same contest. This comparison is the core of the Phase-0 gate.

### 4.7 Alternatives considered and rejected

- **GPDD** (Shi, TCAD 2013) — graph-pair reduction, **cancellation-free**, unlike
  DDD which can generate terms that later cancel. **[LIT]** its own conclusion is
  that runtime/memory are only *comparable* to DDD despite generating far more
  terms.
- **PDD / HPDD** (Lasota 2018) — higher-order summative cofactors, cancellation-
  free, built **directly from a netlist**, with a hierarchical variant and native
  pathological elements (nullors, mirrors).

**Conclusion 4.7** — both are rejected for the *same architectural reason*, and it
is decisive: they are **topology/netlist-driven**, while pycircuit is
**MNA-matrix-driven** — elements stamp into `G`/`C`/`u` and there is no
topological circuit-graph model to build from. Adopting either means a second
front end and a second engine, not another strategy on the existing toolkit axis.
DDD/MTDDD consume the matrix we already build. (Noted so this is not
re-litigated. Caveat: the PDD papers are paywalled; assessment is from abstracts.)

### 4.8 Independent corroboration of the SoE result

**[LIT]** Filaretov & Gorshkov (2020) generate the s-expanded form "with every
coefficient being a compact-nested expression", reporting results more compact
than commercial CAS factorisation. That is independent confirmation of
**[OURS]** §3's sharing-beats-inlining finding, and it points at the same
synthesis our roadmap is converging on: **s-expanded coefficient split + nested /
shared coefficient expressions**.

Their tool **CirSym** (freeware, C source, SPICE netlist input) is usable as an
**independent oracle** for cross-checking our compact output and operation counts.

---

## 5. Where DDD fits the architecture

The framing that makes this non-breaking: a "toolkit" today conflates two
independent axes.

- **Axis 1 — stamping backend**: the primitives elements use to build `G`, `C`,
  `u` (`_numeric`, `_symbolic`, `_jaxtoolkit`, `_sparse_numeric`).
- **Axis 2 — solve / representation strategy**: how `A x = b` becomes transfer
  functions, poles and noise.

Every symbolic toolkit is `SomeToolkit(_symbolic)` — they are **identical on
Axis 1 and differ only on Axis 2**. So `symbolic_poly`, `ginac` and a new `ddd`
are interchangeable *representation* strategies over the same sympy-stamped
system. Stamping is untouched; elements and `Circuit` never learn that DDD exists.

DDD adds one genuinely new thing: **the graph's shape depends only on matrix
sparsity/entry identity, while the per-vertex operation `D0 + sign·entry·D1` is
where arithmetic happens.** So DDD composes with a pluggable **arithmetic
backend** — sympy, GiNaC, or plain numerics — and the other symbolic toolkits get
a second life as DDD's arithmetic engines rather than as competing solvers.

Concrete shape (all additive):

1. **`DDDResult`** — sibling of `GinacResult` / `SoESolution`; holds the
   determinant DDD, builds numerator DDDs lazily via Cramer, exposes
   `denominator()`, `numerator(i)`, `tf(i,j)`, `eval(params)`, `poles()`,
   `approximate(tol)`. Never flattens.
2. **A small arithmetic interface** — `add`, `mul`, `zero`, `one`,
   `evaluate(entry, env)` — with sympy / GiNaC / numeric adapters.
3. **`DDDToolkit(SymbolicPolyToolkit)`**, e.g. `DDDToolkit(arithmetic=ginac_toolkit)`.
   It inherits `supports('num_den')`, so `AC(..., toolkit=ddd_toolkit)` works
   through the existing branch. Two integration depths:
   - *conservative*: `linearsolver_num_den` builds DDDs and returns sympy
     `num/den` — correct, zero API change;
   - *full*: add `supports('ddd')` plus a `CircuitResultACPoly` subclass so
     `poles`/`eval`/`approximate` run **on the compressed graph**.

---

## 6. Conclusions in brief

Each conclusion states its reason; §-references point to the evidence.

1. **Rewrite `DDD_of_matrix` using the ICCAD-2010 LED algorithm**, keeping `Node`
   and its tests — because the existing function has never run and its three
   defects are structural rather than typos (§2), and because the 2010 algorithm
   removes the ordering pre-pass, the BDD package and the sign pass that the 2000
   flow needs (§4.1). Sharing via a **minor hash table** is non-negotiable: without
   it there is no compression, which is the only reason to build a DDD at all.
2. **The deliverable is a multiroot s-expanded DDD, not a flat determinant DDD** —
   because the flat form does not match any interface we have, whereas the
   s-expanded coefficient form is exactly what `as_num_den` produces and
   `poles()`/`zeros()` consume, and because Theorem 1 makes the expansion *linear*
   in the flat DDD, so it is cheap to add and irrational to omit (§4.3).
3. **Numeric terminals (MTDDD) are core, not an optimisation** — because
   pycircuit's matrices are mostly numeric in normal use, and because folding
   numerics into terminal values during construction structurally prevents the
   exact-rational blow-up that capped the GiNaC backend at ~dim 16 (§3, §4.3).
4. **Approximation comes after s-expansion, never before** — because approximating
   a non-s-expanded function treats the powers of `s` unequally and gives
   unreliable results (§4.4).
5. **Hierarchy is later** — not because it is weak (it is the largest single
   reported gain, 56× on µA741) but because it *multiplies* whatever the flat
   engine achieves, so its value is unknown until the flat engine is measured
   (§4.5).
6. **Reject GPDD and PDD** — netlist/topology-driven versus our MNA-matrix
   architecture, and GPDD's own results show no performance gain to pay for that
   mismatch (§4.7, with reconsideration conditions in §9).
7. **DDD enters as an Axis-2 representation strategy with a pluggable arithmetic
   backend** — because all symbolic toolkits already differ *only* on that axis
   (§5), so this is the one shape that adds a strategy without touching stamping,
   existing toolkits, or the analyses.
8. **The bar is SoE, not sympy** — because we already own a representation that
   does not blow up (§3), so beating sympy proves nothing. DDD must beat
   `|SoE ops|` on circuits where `symbolic_poly` fails to earn its complexity.

---

## 7. Priorities and staged plan

Each stage ends in a measurable gate. Stop or re-plan if a gate fails.

### P0 — Core construction + measurement (gate: go / no-go for everything else)

Rewrite `DDD_of_matrix` per ICCAD 2010: layered expansion, minor hash table keyed
by `(row-index tuple, col-index tuple)`, min-degree expansion ordering, signs
during construction. Add numeric evaluation over the graph. Standalone module and
tests — **no toolkit wiring yet**.

Measure, on fully-symbolic ladders (N = 4…16), the MFB filter from Example 10, and
a larger opamp-like circuit:

- `|DDD|` (= add/multiply count) vs **our existing SoE op counts** (73/157/241/325);
- `|DDD|` vs flat product-term count, and vs the `n·2^(n-1)` no-sharing line;
- construction time, and numeric evaluation time per sweep point.

**Gate**: `|DDD|` must be competitive with, and ideally well below, SoE ops on
circuits where `symbolic_poly` is known to fail. If it is not, DDD does not earn
its complexity here and SoE + GiNaC remain the answer — record the negative
result (as we did for symengine and GiNaC) and stop.

### P1 — s-expanded / multiroot DDD

Implement the coefficient split: one root per power of `s`, roots sharing
subgraphs. Validate Theorem 1's linearity empirically (`|s-expanded DDD|` vs
`q·d·|DDD|`).

**Gate**: coefficients match `symbolic_poly` / `_ginac.poly_coeffs` on circuits
small enough for both, and size grows linearly in `|DDD|`.

### P2 — `DDDResult` + conservative toolkit integration

`DDDResult` with lazy Cramer numerators; `DDDToolkit.linearsolver_num_den`
returning sympy `num/den`. Existing `AC` path works unchanged; existing tests must
stay green.

### P3 — Semi-symbolic numeric terminals (MTDDD)

Extend terminals to carry numeric values so numeric sub-products collapse during
construction. This is pycircuit's dominant regime; pull it earlier if P0 shows
numeric-heavy matrices dominating.

**Gate**: no exact-rational blow-up at the dimensions where GiNaC stalled (~16+).

### P4 — Noise via a single multiroot DDD

All per-noise-source transfer functions in one multiroot DDD.

**Gate**: agrees with the existing `noise_psd` result, at a cost comparable to a
single transfer function.

### P5 — Approximation / dominant terms

`approximate(tol)` on the s-expanded form; dominant poles/zeros as ratios of
consecutive coefficients. Prune using component *variation ranges*, not a single
nominal value (Yu & Sechen).

### P6 — Hierarchy

Symbolic stamps per `SubCircuit`.

### Cross-cutting

- Live doc page (`exec-rst`, as for `soe_symbolic.rst` and `ginac_native.rst`)
  regenerating the P0/P1 comparison tables at build time.
- Cross-check selected results against **CirSym**.
- Add Tan & He, *Symbolic Analysis and Reduction of VLSI Circuits* (2005) to the
  reading list — `test_ddd.py` already cites it.

---

## 8. Open questions and risks

1. **[VERIFY] Does DDD actually beat SoE for us?** The one published head-to-head
   is on different circuits with a different flavour of SoE. This is the whole
   P0 gate.
2. **Variable/expansion ordering on our matrices.** Min-degree is a heuristic;
   TCAS-II 2010 only settles the dense case. A bad ordering on a real circuit
   could blow up the vertex count — the main technical risk.
3. **Cancellation.** DDD expands a determinant and can generate terms that later
   cancel; GPDD exists precisely because of this. We accept it (§4.7) but should
   watch for it on circuits with dependent sources and nullors.
4. **Python constant factors.** All the cited runtimes are C/C++. A pure-Python
   LED build may be fine (it is dict operations plus small arithmetic) but the
   arithmetic backend choice may matter more than in the papers.
5. **Where numerics enter.** MTDDD's numeric terminals trade exact rationals for
   floating values. Need to decide, per use, whether coefficients stay exact
   (poles from exact polynomials) or go numeric — this is exactly the trade that
   bit us in the GiNaC `_ginac_coeffs` gating bug.

---

## 9. Deliberately not doing — and why

Each of these is a live option we are declining, so the reason is recorded along
with what would make us reconsider.

**Not implementing GPDD or PDD.** Two independent reasons, either sufficient.
*Architecturally* they are netlist/topology-driven while pycircuit is
MNA-matrix-driven, so adopting one means writing a circuit-graph front end *and* a
new arithmetic (HOSC, for PDD) — a second engine beside the existing toolkit axis,
not an addition to it. *On merit*, *[LIT]* GPDD's own conclusion is that its
runtime and memory are only **comparable** to DDD, so that cost buys no measured
speed; its distinctive property is cancellation-freedom, not performance.
→ **Reconsider if** cancellation proves to be a real practical problem for us
(risk §8.3) — then cancellation-freedom becomes worth the front end — or if a
topological model gets built for some unrelated reason.

**Not building a topological/graph front end for the circuit.** pycircuit's
element model is stamping-based: every element knows how to contribute to `G`,
`C` and `u`. A graph model would require expressing topological semantics for
each element *again* — including dependent sources, nullors and the whole
semiconductor set — leaving two representations of the same knowledge to keep in
sync. That is a large, permanent maintenance cost whose only consumers would be
the methods rejected above.
→ **Reconsider if** something else needs it (topology-based structural checks,
say), at which point the cost is shared.

**Not using an off-the-shelf BDD/ZBDD package.** *[LIT]* ICCAD 2010 demonstrates
that sharing by hashing minor indexes needs no package at all, so the dependency
would buy nothing we cannot get from a dict. It would also cost us: a C library
(or its bindings) complicates installation for what is otherwise a pure-Python
library, and the ZBDD route drags back in the two passes — the symbol-ordering
pre-pass and the separate sign pass — that the 2010 algorithm eliminates.

**Not implementing dynamic variable reordering (sifting) in P0.** Ordering is the
classic BDD escape hatch and *[LIT]* finding an optimal order is NP-hard in
general, so this could absorb unbounded effort. We start with the min-degree
expansion ordering the 2010 paper uses and *measure*. Reordering is only worth its
complexity if the measurement shows ordering to be the binding constraint.
→ **Reconsider if** P0 shows vertex counts blowing up on real circuits and traced
to ordering (risk §8.2).

**Not touching Axis 1 (stamping) or `numeric` / `symbolic` / `symbolic_poly`
behaviour.** This is partly an external constraint — those toolkits have
downstream users and forks — but mainly it is what makes the work *cheap to
abandon*. The whole point of the toolkit class hierarchy was to add features
without core surgery; if the P0 gate fails we delete a module instead of unpicking
edits threaded through the analyses. Given that two of the last three backend
attempts ended in negative results (§3), keeping this one independently revertible
is worth more than any convenience gained by relaxing it.

**Not pursuing symengine further; GiNaC kept only as an arithmetic backend.** The
distinction between the two matters, because it is the reason one is dropped and
the other retained. *[OURS]* symengine's failure is **intrinsic to the library**:
its `LUsolve`/`det` are not fraction-free and `FFLU` does not cancel, so
expressions swell no matter how we drive it — there is no input we can supply that
fixes that. GiNaC's failure was **a property of the inputs**: `normal()` does
cancel correctly, and the blow-up came from numeric component values becoming huge
exact rationals, which frequency scaling to O(1) coefficients demonstrably
controls. A library whose failure we can eliminate by conditioning its input is
worth keeping as the per-vertex arithmetic engine (§5); one whose algorithm is
wrong for the job is not.
