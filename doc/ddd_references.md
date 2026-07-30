# DDD (Determinant Decision Diagram) — reference papers

**Companion document: `ddd_conclusions.md`** records what we concluded from these
papers, the reasoning, and the staged priorities. This file is the reading list;
that one is the argument. Use both as input to the implementation plan.

Reading list for the Determinant Decision Diagram work (improvement #4 in
`ginac_fully_symbolic.md`): compact, canonical representation of a symbolic
network determinant for symbolic analysis of *large* / *fully symbolic* circuits,
with fast numeric evaluation and dominant-term (approximate) simplification.

Links marked **(open)** are author-hosted or open repositories and should
download directly. Links marked **(may need access)** are on academia.edu /
ResearchGate / Springer and can prompt for a login. The author PDF directories
at the bottom hold many more of these papers.

Design context this feeds: rewriting the broken `DDD_of_matrix` on top of the
existing `pycircuit/circuit/ddd.py` `Node`, a lazy compressed `DDDResult`
(sibling of `GinacResult` / `SoESolution`), a pluggable *arithmetic backend*
(sympy / GiNaC / numeric), and `approximate(tol)`.

Local PDF copies of the SJTU papers (section A2/F) are in
`~/pycircuit_agy/papers/ddd/` — see the filenames noted per entry.

## A. Core DDD method — the blueprint to implement

**Shi & Tan, "Canonical Symbolic Analysis of Large Analog Circuits with
Determinant Decision Diagrams," IEEE TCAD 19(1), 2000.**
- (open) https://labs.ece.uw.edu/mscad/shi/Papers/TCAD_2000_Jan_1.pdf
- (open, mirror) https://intra.engr.ucr.edu/~stan/papers/tcad00_1.pdf
- Why: *this is DDD.* Defines the diagram, construction from the circuit matrix,
  vertex sharing, signs, cofactors, and how to read off network functions —
  exactly what a correct `DDD_of_matrix` and `DDDResult`
  (`denominator`, `numerator(i)` via cofactors, `tf`) need. Quantifies the
  "vertices ≪ product terms" claim we want to measure in Phase 0.
- **But do not implement its construction flow — implement A2 instead.**

## A2. The construction algorithm we should actually implement

**Shi, "A Simple Implementation of Determinant Decision Diagram," ICCAD 2010,
pp. 70–76.**
- (open, Wayback) https://web.archive.org/web/20240414194754if_/https://aice.sjtu.edu.cn/msda/data/publication/ICCAD2010.pdf
- local: `~/pycircuit_agy/papers/ddd/ICCAD2010.pdf`
- Why: **this is the single most useful paper for us.** It replaces the TCAD-2000
  five-step flow with a much simpler one that is a far better fit for a pure-Python
  implementation. It deletes three things we would otherwise have to build:
  - **no BDD/ZBDD package** — sharing is done by hashing *minors* directly, keyed
    by their (row-index tuple, column-index tuple). A plain dict is the whole
    sharing mechanism. (TCAD-2000 needs ZBDD `Change()`/`Union()` vertex-triple
    sharing.)
  - **no symbol-ordering pre-pass** — the full symbol order is weakened to an
    *expansion order* (which row/column to expand next), chosen on the fly by a
    min-degree heuristic. Order *within* a sibling group provably doesn't affect
    size (Property 4), so no Greedy-Labeling phase.
  - **no separate sign-determination pass** — signs fall out during construction
    because the minor indexes are retained.

  Mechanism: "Layered Expansion Diagram" (LED) — expand a whole selected
  row/column at a time into queues, one queue per layer, minors shared via the
  hash table; then convert LED→DDD mechanically (existing next-layer pointers
  become 1-arrows; add 0-arrows along each sibling group, terminating at `0`;
  bottom-queue elements terminate at `1`). |DDD| = total elements across queues.
  Results: LED+rowwise beats the traditional Greedy-ordered build by >30× in DDD
  size at n=18, and solves the µA725 op-amp by flat expansion in seconds.

**Shi, "Computational Complexity Analysis of Determinant Decision Diagram,"
IEEE TCAS-II 57(10), 2010, pp. 828–832.**
- (open, Wayback) https://web.archive.org/web/20240416073219if_/https://aice.sjtu.edu.cn/msda/data/publication/TCAS2_DDD_2010.pdf
- local: `~/pycircuit_agy/papers/ddd/TCAS2_DDD_2010.pdf`
- Why: the **go/no-go arithmetic for Phase 0**. Proves a row-wise (or column-wise)
  order is *optimal* for a full n×n matrix, and that the optimal DDD size is then
  exactly `n·2^(n-1)`. So DDD is still exponential in the worst case — the win
  comes entirely from *sparsity*, and this gives us the yardstick to measure our
  measured vertex count against. It also settles the ordering question for the
  dense case (use row-wise) so we don't go hunting for a heuristic.

## A3. Semi-symbolic (MTDDD) — the regime pycircuit is *actually* in

**Pi & Shi, "Multi-Terminal Determinant Decision Diagrams: A New Approach to
Semi-Symbolic Analysis of Analog Integrated Circuits," DAC 2000.** (Univ. of
Washington)
- (open) https://www.cs.york.ac.uk/rts/docs/SIGDA-Compendium-1994-2004/papers/2000/dac00/pdffiles/01_5.pdf
- local: `~/pycircuit_agy/papers/ddd/dac2000_york.pdf`
- Why: **the highest-value paper here after A2**, and the one that matches how
  pycircuit is actually used: most parameters numeric, a few symbolic, plus `s`.
  Two ideas, both of which land directly on code we already have:

  1. **Multiple numeric terminals.** A DDD has only the `0` and `1` terminals; an
     MTDDD adds terminals carrying arbitrary nonzero *numeric values*. Numeric
     sub-products are therefore collapsed into a single terminal value **during
     construction** instead of being carried symbolically. Result: far fewer
     symbolic vertices, faster repeated numeric evaluation, and only the
     parameters of interest survive as symbols (so the output stays readable).
     Note this is also a structural answer to the coefficient blow-up that killed
     the GiNaC backend at ~dim 16 (see `ginac_fully_symbolic.md`): there, numeric
     component values became exact rationals whose products exploded into huge
     CLN integers. MTDDD never forms those — numerics are folded into a terminal
     as they arise.
  2. **The "coefficient MTDDD"** — a *multi-root* MTDDD in which each root is the
     coefficient of one power of `s` in the numerator or denominator, i.e. the
     s-expanded rational `H(p,s) = Σ nᵢ(p)sⁱ / Σ dⱼ(p)sʲ`. **This is exactly the
     form `TransferFunction.as_num_den` / `_ginac.poly_coeffs` already produce and
     that `poles()`/`zeros()`/bode consume.** Built from the "complex DDD" (the
     DDD of the MNA matrix, whose entries embed `s`) by a graph algorithm, with
     sub-graph sharing across coefficients. Exact — no approximation.

  Reported result: full semi-symbolic transfer function for the µA741 in under a
  minute — on a 450 MHz Pentium-II, in 2000.

  Caveat worth knowing: this paper is **not** covered by Shi's 2011 survey (E) —
  zero mentions — because the survey follows the UW/SJTU fully-symbolic BDD line.
  Easy to miss.

## A4. s-expanded (multiroot) DDD — the same idea, done rigorously

**Shi & Tan, "Compact Representation and Efficient Generation of s-Expanded
Symbolic Network Functions for Computer-Aided Analog Circuit Design," IEEE TCAD
20(7), July 2001, pp. 813–827.**
- (open, Wayback) https://web.archive.org/web/20240416213605if_/http://labs.ece.uw.edu/mscad/shi/Papers/TCAD_2001_July_1.pdf
- local: `~/pycircuit_agy/papers/ddd/TCAD_2001_Jul_sexpanded.pdf`
- Why: **this is the definitive treatment of the coefficient-form idea**, and the
  companion to A3. Where MTDDD introduced multiple roots for the semi-symbolic
  case, this paper states it generally and proves it. "s-expanded DDD" and
  "multiroot DDD" are the same object: one root per coefficient of a power of `s`
  in the numerator and denominator, with all roots sharing common subgraphs.
  Again, that is exactly `TransferFunction.as_num_den`'s output shape.
  - **Theorem 1** is the reason to care: from a complex DDD of size `|DDD|`, the
    s-expanded DDD is built in time proportional to `q·d·|DDD|` and has no more
    than `q·d·|DDD|` vertices, where `q` = degree of the denominator polynomial
    and `d` = max devices attached to a node (and under MNA "compact symbol"
    form, `d ≤ 2`). So **s-expansion is linear in the complex DDD**, even though
    the number of s-expanded product terms is astronomically larger.
  - The numbers make the point: the µA741 denominator has 108 032 *complex*
    product terms, whose s-expansion is ~7.8×10³⁴ product terms — represented by
    an s-expanded DDD of 99 844 vertices, built in a few seconds on an
    UltraSPARC-I. Under full-symbol representation the term count rises nine
    orders of magnitude while the DDD grows only ~3× (297 115 vertices).
  - **Noise** (§V, directly relevant to us): noise analysis is a *set* of transfer
    functions, one per noise source, which share most subexpressions — so
    represent them all in a **single multiroot DDD**, at cost comparable to one
    transfer function. This is the same insight as `SymbolicPolyToolkit`'s
    shared-denominator `noise_psd`, generalised and scaled.
  - Ordering consequence for our roadmap: §I notes approximation *requires* the
    s-expanded form, so that coefficients of each power of `s` are approximated
    on an equal footing — otherwise the result is unreliable. Likewise dominant
    poles/zeros come out as ratios of coefficients of consecutive powers of `s`.
    **So s-expansion is a prerequisite for section B, not a parallel option.**

## B. Approximate / dominant-term analysis — the "readable formula for a big circuit" payoff

Note: per A4, do this **on the s-expanded form**, not on a raw network function.

**Tan & Shi, DDD-based dominant-term generation** (behavioral-modeling line;
shortest-path & dynamic-programming term generation, linear in the number of DDD
vertices).
- (may need access) Overview: https://link.springer.com/content/pdf/10.1023/A:1015041927036.pdf
- Why: the algorithm behind `DDDResult.approximate(tol)` — extracting the few
  dominant product terms directly on the graph in linear time. Turns a big
  circuit's exact-but-huge function into a short, interpretable expression.
  (The s-expanded representation it consumes is A4.)

**Yu & Sechen, "A unified approach to the approximate symbolic analysis of large
analog integrated circuits," IEEE TCAS-I 43(8), 1996.**
- (may need access) https://www.researchgate.net/publication/3224818
- Why: classic "approximate-during-computation" strategy and, crucially, using
  component *variation ranges* rather than a single nominal value when deciding
  which terms to drop. Informs *how* we prune so the simplification stays valid
  across the design space.

## C. Hierarchical DDD — scaling to very large circuits

**Tan & Shi, "Hierarchical Symbolic Analysis of Analog Integrated Circuits via
Determinant Decision Diagrams," IEEE TCAD 19(4), April 2000, pp. 401–412.**
- (open, Wayback) https://web.archive.org/web/20240416164453if_/http://labs.ece.uw.edu/mscad/shi/Papers/TCAD_2000_Apr_1.pdf
- local: `~/pycircuit_agy/papers/ddd/TCAD_2000_Apr_hierDDD.pdf`
- Why: **the primary hierarchical-DDD paper** (the DAC 2004 entry below is the
  later, shorter follow-up). Method: suppress each subcircuit to its terminals in
  terms of its matrix determinant and cofactors, represent those with DDDs, then
  apply Cramer's rule at the top level. Illustrated on a Cauer low-pass filter.
- **It also contains the head-to-head we most need**, because its baseline SCAPP
  is a *sequence-of-expressions* analyser — i.e. the same representation as our
  `soe.py`. Their claims, on cascaded-opamp circuits: the DDD is "much more
  compact than the sequence-of-expression representation used in SCAPP"; DDD size
  grows almost linearly in circuit size while product terms grow exponentially;
  and DDD beats both SCAPP and SPICE on repetitive evaluation. On a µA741,
  three-level two-way hierarchical DDD needs **117 vertices vs 6654 for flat DDD
  (56×), vs 119 011 sum-of-product terms**.
- Crucially it gives us a **directly comparable metric**: each DDD vertex costs
  one addition and one multiplication, so `|DDD|` is an operation count — exactly
  what we already report for SoE in `soe_symbolic.rst` (73/157/241/325 ops at
  N=4/8/12/16). Phase 0 can therefore put DDD vertices and SoE ops on one axis
  for the same circuits. Treat "DDD ≪ SoE" as a claim to *verify*, not assume:
  their SCAPP baseline is hierarchical-suppression SoE on cascaded opamps, while
  ours is Gaussian-elimination SoE measured on ladders, where we found linear
  growth.

**Tan & Shi, "Hierarchical Approach to Exact Symbolic Analysis of Large Analog
Circuits," DAC 2004.**
- (open) https://intra.engr.ucr.edu/~stan/papers/dac04_1.pdf
- Why: builds DDDs *per subcircuit* (suppress to terminals via
  determinants/cofactors, then Cramer at the top). Maps directly onto pycircuit's
  `SubCircuit` hierarchy; the route to 100+ node circuits without one flat matrix.

**Xu, Shi & Li, "Hierarchical Exact Symbolic Analysis of Large Analog Integrated
Circuits by Symbolic Stamps," ASP-DAC 2011.**
- (open, Wayback) https://web.archive.org/web/20240418011841if_/https://aice.sjtu.edu.cn/msda/data/publication/ASPDAC2011_xuhui.pdf
- local: `~/pycircuit_agy/papers/ddd/ASPDAC2011_xuhui.pdf`
- Why: the newer take on hierarchy, and a notably pycircuit-shaped one — a
  subcircuit is reduced to a **"symbolic stamp"**, i.e. an n-terminal admittance
  block that stamps into the parent MNA matrix exactly like a primitive element
  does. That is precisely how `SubCircuit` already composes, so this maps onto our
  code with less new machinery than the DAC-2004 cofactor-suppression route.

**Song & Shi, "Hierarchical Graph Reduction Approach to Symbolic Circuit Analysis
with Data Sharing and Cancellation-Free Properties," ASP-DAC 2012.**
- (open, Wayback) https://web.archive.org/web/20240415085657if_/https://aice.sjtu.edu.cn/msda/data/publication/songyang_aspdac12.pdf
- local: `~/pycircuit_agy/papers/ddd/songyang_aspdac12.pdf`
- **RE-READ 2026-07-30 and PROMOTED. The old note said "relevant only if we take
  the graph route; listed for completeness" — that was wrong, and it cost a
  session.** This paper contains the diagnosis of the central obstacle in
  `doc/cancellation_ranking_conclusions.md`, and it was on the shelf the whole
  time.
- Why, in three parts:
  1. **§II identifies our cancellation as a formulation artifact.** Its worked
     example is our own compact-symbol MNA convention: with `a=G1, b=c=-G1,
     d=G1+G2+G3, e=f=-G3, g=G3+G4`, `det = adg - aef - bcg`, and expanding the
     composites gives `G1*G2*G3 + G1*G2*G4 + G1*G3*G4` — three terms, **all
     positive**, so the cancellation-free form has `kappa = 1`. Our measured
     `kappa = 9.4e3` on the µA741 therefore describes MNA-plus-composite-symbols,
     not the amplifier.
  2. **It states our stage-5 result five years earlier:** *"since DDD is not
     cancellation-free, every layer of the hierarchy has the term cancellation
     problem."* Our "hierarchical approximation is not compositional" is a
     rediscovery; the quantification (non-monotone 1.47e-2 → 1.13e-1) is what we
     added.
  3. **§III is "Schur Decomposition versus Symbolic Stamp"** — the exact choice
     stage 5 made — and it prefers the stamp as more partition-friendly, then
     recommends hierarchical GPDD *"whenever the property of cancellation-free is
     required"*.
- Also useful: their Tables I–II give the **module partition** of the µA725 (which
  devices in which block) but **not its connectivity**, so the µA725 decline
  recorded below still stands.
- Caveat they state and we confirmed: in plain double-precision *evaluation*
  cancellation is not an accuracy problem (we measured `kappa * 2**-53 = 1.6e-11`).
  It bites in *variational* analysis, sensitivity, and — our addition —
  **approximation**.

## D. Foundational data structure

**Bryant, "Graph-Based Algorithms for Boolean Function Manipulation," IEEE TC
35(8), 1986.**
- (open) https://www.cs.cmu.edu/~bryant/pubdir/ieeetc86.pdf
- Why: DDD is a BDD variant. reduce/apply, the unique-node hash table (how
  sharing is *enforced*), canonical form, and — most important for us —
  **variable ordering**, all transfer directly. The reference for making the
  diagram actually compact and for the ordering heuristic (the main risk).

## D2. The main alternative to DDD — cancellation-free graph methods (GPDD)

DDD expands a *determinant*, so it can generate product terms that later cancel.
Shi's second line of work builds the diagram from the circuit *graph* instead, and
is cancellation-free by construction. Worth knowing before committing to DDD.

**Shi, "Graph-Pair Decision Diagram Construction for Topological Symbolic Circuit
Analysis," IEEE TCAD 32(2), 2013, pp. 275–288.**
- (open, Wayback) https://web.archive.org/web/20240415232304if_/https://aice.sjtu.edu.cn/msda/data/publication/2013_TCAD_GPDD.pdf
- local: `~/pycircuit_agy/papers/ddd/2013_TCAD_GPDD.pdf`
- Why: the mature statement of GPDD — build a pair of graphs from the small-signal
  circuit and reduce them successively into a BDD. **Cancellation-free**, unlike
  DDD, and its terms are composed *directly of circuit parameters rather than
  intermediate symbols*, which is the interpretability property stages 3–5 of the
  approximation work chased and never got.
- **The old note's conclusion is WITHDRAWN.** It read: "GPDD runtime/memory is
  only *comparable* to DDD despite generating many more terms. That is a real
  argument for staying with DDD." The *facts* are right — Table VIII gives
  |GPDD| 10 432 / 17 488 / 197 274 against |DDD| 3 579 / 11 506 / 62 794, so ~3×
  the size, with times 0.682/0.793/6.771 against 0.586/2.042/10.359, i.e.
  comparable and sometimes faster. But the conclusion weighed only speed and
  memory and treated cancellation-freedom as optional. **We then measured
  cancellation to be the binding constraint on approximation**, which inverts the
  trade: ~3× size is cheap for the property that decides whether dominant-term
  extraction converges at all.
- What is still true from the old note, and is the real cost: GPDD needs a
  **circuit graph**, and pycircuit has an MNA matrix and no graph model. That is
  the work item, not the runtime.

**Shi & Chen, "A Graph Reduction Approach to Symbolic Circuit Analysis,"
ASP-DAC 2007.**
- (open) https://cecs.uci.edu/~papers/aspdac07/pdf/p197_2C-2.pdf
- local: `~/pycircuit_agy/papers/ddd/ASPDAC2007_chenweiwei.pdf`
- Why: the shorter, earlier, more readable version of the GPDD idea (graph
  reduction + recursive sign determination). Read this before the 2013 TCAD paper
  if we ever evaluate the graph route.

## E. Context / survey and the semi-symbolic case

**Shi, "A survey on binary decision diagram approaches to symbolic analysis of
analog integrated circuits," Analog Integr. Circ. Sig. Process., 2011.**
- (open, Wayback) https://web.archive.org/web/20240415114254if_/https://aice.sjtu.edu.cn/msda/data/publication/AICSP_survey_2011_gshi.pdf
- (may need access) https://link.springer.com/article/10.1007/s10470-011-9773-8
- local: `~/pycircuit_agy/papers/ddd/AICSP_survey_2011_gshi.pdf`
- Why: a map of the whole BDD-for-symbolic-analysis field by the author of both
  DDD implementations and GPDD; places DDD relative to alternatives before we
  commit. **Start here** — it is the cheapest way to get the whole landscape, and
  the author-hosted PDF above avoids the Springer paywall.

**Shi, Tan & Tlelo-Cuautle, *Advanced Symbolic Analysis for VLSI Systems:
Methods and Applications*, Springer, 2014.** (book)
- (may need access) https://link.springer.com/book/10.1007/978-1-4939-1103-5
- Ch. 4 "Determinant Decision Diagrams": https://link.springer.com/chapter/10.1007/978-1-4939-1103-5_4
- Why: the only *consolidated, post-2000* treatment — the three principals of the
  field writing up DDD, GPDD, approximation and applications in one place, 14
  years after the original papers. If any single item is worth buying/borrowing,
  it is this. Not needed to start Phase 0 (the papers above cover it).

(The semi-symbolic DAC 2000 paper formerly listed here is the MTDDD paper — it has
been corrected and promoted to section A3.)

**Shokouhifar, Yazdanjouei & Weber, "Direct Simplified Symbolic Analysis (DSSA)
Tool," arXiv:2510.15901, Sep 2025.**
- (open) https://arxiv.org/pdf/2510.15901
- Why: the only genuinely *recent* work found on this topic. It is **not** a DDD
  paper — it skips the matrix/graph representation entirely and treats
  simplification as a modelling problem, extracting the dominant transfer-function
  terms with Monte-Carlo sampling driven by a genetic algorithm. Useful as a
  contrast to `approximate(tol)`: term selection by *sampled numeric significance*
  over a parameter distribution, rather than by walking a graph. Low priority for
  Phase 0; potentially interesting later, and it composes with our existing
  numeric-evaluation path.

## D3. Parameter Decision Diagrams (PDD) — a separate, later school

A third representation family, from the Filaretov / Rodanski / Lasota line
(Russian & Polish groups), developed independently of the UW/SJTU BDD work and
running later, to ~2020. Assessed and **not** recommended as an implementation
target for pycircuit — reasons below — but worth knowing about.

**Lasota, "Symbolic analysis of electric networks with higher order summative
cofactors and parameter decision diagrams," Int. J. Circuit Theory Appl. 46,
2018.**
- (paywalled, Wiley) https://doi.org/10.1002/cta.2495
- What: introduces *higher order summative cofactors* (HOSC — "a cofactor of a
  cofactor") as the arithmetic, giving a cancellation-free symbolic technique
  that builds a BDD called a **parameter decision diagram** *directly from the
  netlist*. Rolling up already-analysed parts yields a multilevel **HPDD**
  (hierarchical PDD) that stays cancellation-free at every level and compresses
  via the circuit's self-similarity. Pathological elements (nullors, mirrors)
  fall out naturally.
- Why **not** for us: the same architectural objection as GPDD (D2), and it is
  the decisive one. PDD is *netlist/topology*-driven; pycircuit is
  *MNA-matrix*-driven — elements stamp into `G`/`C`/`u` and there is no
  topological circuit-graph model to build a PDD from. Adopting it means writing
  a parallel front end plus a new arithmetic (HOSC), i.e. a second engine rather
  than another strategy on the existing toolkit axis. DDD/MTDDD, by contrast,
  consume the matrix we already build.
- Worth noting anyway: native handling of nullors is interesting given pycircuit
  *has* a `Nullor` (it carries the MFB example), which currently goes through
  ordinary MNA stamps.

**Filaretov & Gorshkov, "Efficient generation of compact symbolic network
functions in a nested rational form," Int. J. Circuit Theory Appl. 48, 2020.**
- (paywalled, Wiley) https://doi.org/10.1002/cta.2789
- Why this one *does* matter to us, more than PDD itself: it is independent
  confirmation of the central finding from our own SoE experiment — that a
  **nested / shared** representation stays compact where an expanded one does
  not. They generate the s-expanded form "with every coefficient being a
  compact-nested expression," by implicit parameter extraction and factoring by
  grouping, and report expressions more compact than the factorisation
  algorithms of commercial CAS. That is precisely the combination our roadmap is
  converging on: MTDDD's per-power-of-`s` coefficient split (A3) with SoE's
  nested/shared coefficient expressions (`soe.py`, and see the "sharing vs
  inlining" result in `soe_symbolic.rst`).
- **CirSym** — their implementation, freeware, C source at
  https://github.com/k-gorshkov/cirsym, online at http://intersyn.net/en/cirsym.html,
  input is a lightly-modified SPICE netlist. Directly useful as an **independent
  oracle**: feed it the same netlist and compare our nested/compact output and
  operation counts against a published tool.

*Caveat: both Wiley papers are paywalled; the notes above are from abstracts and
secondary sources, not full text. Get the PDFs before implementing anything from
this section.*

## F. How current is this field?

Searched (Jul 2026) for newer DDD work, post-2015 extensions, ML-based symbolic
analysis and recent PhD theses.

The **BDD/DDD line specifically** peaked 1997–2013: its most recent substantive
items are Shi's 2010–2013 SJTU papers (A2, D2) and the 2014 Springer book (E). No
newer thesis or DDD successor surfaced, so implementing against A2 is not chasing
a moving target.

*Symbolic analysis as a whole, however, did not stop there* — the Filaretov /
Lasota parameter-extraction school (D3) runs to 2018–2020 with a working tool, and
the 2025 DSSA paper (B) attacks approximation from a Monte-Carlo/GA direction
outside the BDD world entirely. Neither displaces DDD for our architecture, but
"the field ended in 2013" would be wrong.

## What was actually extracted as tests (surveyed 2026-07-27)

Every paper here was read for reusable benchmarks under one rule: **a published
figure is only usable as a test if the circuit and the measurement are both
understood.** The outcome, so nobody re-runs the survey:

| source | yields | as what |
|---|---|---|
| A2 ICCAD 2010, Table II | full `n × n` matrices, LED and Greedy columns | **exact** — every LED value is `n·2^(n-1)` and ours reproduces it |
| A (TCAD 2000), Table IV | µA741 flat / hierarchical / s-expanded sizes | order-of-magnitude, plus an exact *degree* check |
| A (TCAD 2000), Table II | DDD-vs-SCAPP (= sequence-of-expressions) op counts | ratio check against our `soe.py` |
| A (TCAD 2000) | Cauer low-pass, cascaded op-amp blocks | Tier 2 fixtures |
| C, A3, D2, D3 | — | nothing usable, see below |

Declined, with the reason:

- **µA725** — the tempting one, and the most carefully rejected. Its only
  published schematic (`songyang_aspdac12.pdf`, Fig. 4) draws **no junction
  dots**, so on 26 transistors the connectivity is not recoverable; and the
  papers reporting sizes disagree about the netlist — 34×34 / 1.280604e8 terms
  (A2) against 32×32 / 5.47e7 (D2), while ASP-DAC 2007 says outright "schematic
  omitted". Ambiguous figure, ambiguous target. Reopen only with a dotted
  schematic or a netlist.
- **`miller`, `Cascode`, `bigtst`, `rlctest`, `ccstest`, `vcstest`** — names, not
  specifications.
- **`butter`** — realisation not stated, so its degree cannot be reproduced.
- **`rctreeA` / `rctreeB`** — node counts inferable, but a tree's diagram size
  depends on its *shape*, which is not given.
- **A3 (MTDDD)** — tabulates semi-symbolic µA741 results without saying which
  parameters were kept symbolic, which is the entire content of the measurement.

Two facts worth carrying: published `|DDD|` figures for a *named* circuit pin an
order of magnitude, not a number — the µA741 appears as 23×23 / 6654 / 119 011
terms in A and 25×25 / 13 722 / 4 203 232 terms in A2, and term count is
order-independent. Full matrices have no such ambiguity, which is why they are
the exact check. `benchmarks/paper_extract.py` renders these scans and vector
figures to PNG, which is how the tables above were read at all.

## Author PDF directories (many more open PDFs)

- Sheldon X.-D. Tan (UC Riverside): https://intra.engr.ucr.edu/~stan/papers/
- C.-J. Richard Shi, MSCAD lab (U. Washington): https://labs.ece.uw.edu/mscad/shi/Papers/
- Guoyong Shi, MSDA lab (Shanghai Jiao Tong Univ.) — **host is dead, use the
  Wayback Machine.** The full publication directory was archived in Apr 2024 and
  is the open source of most of the newer papers above. To list it:
  `curl "http://web.archive.org/cdx/search/cdx?url=aice.sjtu.edu.cn/msda/data/publication/*&fl=original,timestamp&collapse=original"`
  then fetch `https://web.archive.org/web/<timestamp>if_/<original>`.

---

Note: the papers already in `~/pycircuit_agy/papers/` (Fang DAC'13, Yao–Wang
ICECS'14) concern *transient / LTE step control*, unrelated to DDD — this is a
separate reading set.

## Added 2026-07-29 — three routes to a *readable* result

Both papers were read at 300 dpi from the local PDFs in
`~/pycircuit_agy/papers/ddd/`. **No number is quoted from either**: they are
recorded for their framing and their stated problem, which is what bears on
our plan. Neither has been reproduced, and the extraction rule applies before
any figure from them is used as a target.

### Qi, Tan, Yu & He — *Wideband Modeling of RF/Analog Circuits via
### Hierarchical Multi-Point Model Order Reduction*

**Why it matters: it states as known the thing we measured independently.**
Its introduction says the sequence-of-expressions methods, the hierarchical
network methods, *and* **"the hierarchical DDD graphs by DDD-based hierarchical
decomposition method are difficult to interpret. The resulting symbolic
expressions are too complicated to gain insights into circuit behavior"**.

That is exactly what we found on the leapfrog: hierarchical reduction gives a
compact diagram (299 vertices for 127 unknowns, fully symbolic) whose entries
are *anonymous generated placeholders* — `_lvl109_16_0` — with `count_ops = 0`.
Compact and uninterpretable. **Note the author: Sheldon X.-D. Tan, of the same
group whose DDD papers our implementation is calibrated against.** The
limitation is theirs to report and they report it.

Two claims worth having:

- **s-domain hierarchical reduction is equivalent to implicit moment matching
  around `s = 0`.** If true, it explains directly why one-point hierarchical
  reduction loses accuracy with frequency — it is a Taylor expansion at the
  origin — and why our numeric-`omega` results are single-frequency by nature.
- **Higher reduction order does not fix that**, because of truncation and
  numerical noise. Their Fig. 2 shows high-order curves *diverging* from exact
  at high frequency rather than converging. **Reconsider-if:** we have not
  reproduced this, and it is the sort of claim that could be an artifact of
  their arithmetic rather than of the method.

Their answer is *not* symbolic interpretation. It is a **pole-residue model**
obtained by expanding at several frequency points, on the argument that poles
near an expansion point are captured more accurately than poles far from it.
So the interpretable output is **poles and residues**, not dominant product
terms.

### Ismail & Friedman — *DTT: Direct Truncation of the Transfer Function — An
### Alternative to Moment Matching for Tree Structured Interconnect*

IEEE TCAD 21(2):131-144, Feb 2003. Reference [10] of the paper above.

**CORRECTION.** An earlier version of this entry described DTT as "a third
route, and the closest to machinery we already have", on the argument that
truncating a rational function directly is the shape `SExpandedDDD` has. That
was written from the **title alone**; the paper had not been read. Having read
the abstract, the characterisation was wrong and the inference does not hold.

What it actually is, from the abstract:

- a method to evaluate **time-domain signals in RLC *trees*** — tree
  structured *passive* interconnect;
- it finds a **low-frequency reduced-order transfer function** by truncating
  the exact one at the nodes of that tree, and determines a **common set of
  poles** characterising every node;
- **numeric** poles, with guaranteed stability below five poles, and
  complexity linear in the number of branches;
- its comparison target is **AWE**, whose problems it names: unstable poles
  even at low order, and numerical breakdown above roughly eight poles.

So it is a numerical interconnect-timing method, an AWE alternative. **It is
not a symbolic-expression technique, and its domain is passive trees.** A
leapfrog filter of five µA741s is neither a tree nor passive, so its
applicability here is *unestablished* rather than close.

**Reconsider-if:** if the goal shifts from symbolic expressions to a compact
*numeric* wideband model — a pole-residue macromodel of the filter — then DTT
and the multi-point MOR paper are both directly relevant, and DTT's stability
argument is the interesting part. For symbolic extraction in the device
parameters they are not.

**The lesson, recorded because it is the rule this project already has.** The
extraction rule says to use external data only when both the thing measured
and the measurement are understood. A title is neither. This entry asserted a
paper's method and its relevance to our machinery from its title, and a
reconsider-if in `hierarchical_approximation_plan.md` was written on that
basis.

**What none of them gives** is what the maintainer asked for in the strongest
form: an expression symbolic *in the device parameters* for a 127-unknown
circuit. Poles and residues are numbers. Term-ranked approximation keeps
symbols but needs a flat diagram. That tension is the real content of
`doc/hierarchical_approximation_plan.md`.

## F. Approximation and term cancellation outside Shi's corpus (added 2026-07-30)

Found while answering "is the cancellation we measured already solved?". Every
entry here was read at least to abstract level; anything not characterised is
marked. Bearing on our work: `doc/cancellation_ranking_conclusions.md` §14.

**Fernández, Sánchez-López, Castro-López & Roca-Moreno, "Approximation Techniques
in Symbolic Circuit Analysis," ch. 7 of *Design of Analog Circuits through
Symbolic Analysis*, Bentham, 2012, pp. 193-226.**
- (open, Wayback) https://web.archive.org/web/2019id_/https://digital.csic.es/bitstream/10261/83058/4/approximation_techniques.pdf
- **The most important item here — it corrected a logical error of ours.** The
  definitive SAG/SBG/SDG survey. Its three simplification criteria, per
  coefficient of `s^k`:
  - (7) ISAAC's, term-vs-largest-term. Their drawback: *"lack of control on the
    accumulated error for each coefficient"*.
  - (8) the standard one, `|Σ deleted| / |Σ all| < ε_M` — a **signed** partial
    sum. *"Mutually canceling terms do not contribute to (8) because they are
    added with their respective signs."* Its limitation, also theirs: *"the
    resulting error at other points may be well beyond ε_M"* — it holds at the
    nominal point only. **`DDD.approximate` already implements exactly this.**
  - (9) `Σ|deleted| / Σ|all| < ε_M`, both sums of moduli — a *robustness* test,
    **not** a relative-error bound; satisfying it at 5% permits real error `0.05κ`.
- So our claim that a ranking "must" capture `1 - tol/κ` of the absolute mass was
  a sufficiency stated as a necessity. Also useful on determinant methods:
  simplification heuristics *"partially palliate the cancellation problem of
  determinant-based approaches"*.
- **Contains no hierarchical section at all** — a meaningful absence.

**Guerra, Roca, Fernández & Rodríguez-Vázquez, "A Hierarchical Approach for the
Symbolic Analysis of Large Analog Integrated Circuits," DATE 2000, pp. 148-152.**
- (open) https://past.date-conference.com/proceedings-archive/2000/DATE00/PDFFILES/01C_2.PDF
- **They considered per-block error budgets and rejected them, for our reason:**
  *"a separate application to each block would require an error propagation
  mechanism at this early stage of the analysis process. This necessarily yields
  more conservative results ... and, consequently, has a negative impact on the
  global performance."* Their alternative — hierarchy for *generation* only, with
  error controlled against a **global numerical reference** — is directly
  implementable for us and is the fix for stage 5's non-monotonicity. Stated
  limitation: *"guaranteed only at the selected frequency sample."*

**Rodríguez-García et al., "An accurate error control mechanism for
simplification before generation algorithms," DATE 1999, pp. 412-416.**
- (open, Wayback) https://web.archive.org/web/2020id_/https://digital.csic.es/bitstream/10261/84968/1/accurate%20error.pdf
- The only genuinely *guaranteed* (non-sampled) error control found, via interval
  extension of the derivatives of the error function — but over **frequency at a
  fixed parameter point**, bounding where the error peaks, not the truncation
  error of a term set. Records the sampling failure honestly: *"the error specs
  are met at the sampling frequencies, but exceeded at intermediate ones."*

**Tan & Shi (C.-J. Richard Shi, U. Washington — NOT Guoyong Shi), "Efficient
Approximation of Symbolic Expressions for Analog Behavioral Modeling and
Analysis," IEEE TCAD 23(4):907-918, 2004.**
- (open) https://escholarship.org/content/qt3pb396zz/qt3pb396zz.pdf
- Same circuit, same compact-entry MNA setup, same magnitude ranking (k-shortest
  path with weight `-log|a_i|`) as ours. Source of the *"70-90% terms are
  canceling terms"* figure. Error control is monitoring, not bounding: *"we
  monitor both magnitudes and phases ... until errors are within the
  user-specified error bounds"*.

**Tan, Guo & Qi, "Hierarchical Approach to Exact Symbolic Analysis of Large Analog
Circuits," DAC 2004, pp. 860-863 (ext. IEEE TCAD 24(8), 2005).**
- (open) https://intra.engr.ucr.edu/~stan/papers/dac04_1.pdf
- *"exact symbolic expressions of a circuit are cancellation-free expressions when
  the circuit is analyzed hierarchically"*, via a *"symbolic decancellation
  process"*. **Scope caveat that matters: this is cancellation-freeness for EXACT
  hierarchical analysis, not for approximation** — it removes exactly-cancelling
  pairs and says nothing about composing truncations.

**Yu & Sechen, "A unified approach to the approximate symbolic analysis of large
analog integrated circuits," IEEE TCAS-I 43(8):656-669, 1996.** Abstract verified.
- Generates common trees of the two-graph *"in the decreasing order of
  magnitude"*. Because two-graph enumeration is cancellation-free, magnitude order
  equals contribution order there — the structural reason this branch never hit
  our wall. Their own limit: *"the limit is imposed mainly by the interpretability
  of the generated symbolic network function."*

**Kolka, Vlk & Horák, "Topology Reduction for Approximate Symbolic Analysis,"
Radioengineering 20(1):252-257, 2011.**
- (open) https://dspace.vut.cz/bitstreams/dedfaf2f-1d6a-40db-88ab-90e6ff1ccf5b/download
- Cleanest statement of why two-graph is cancellation-free: *"If all edges
  represent a unique symbol, there are no two identical tree admittance products
  ... that would cancel each other."* And the honest admission that the guarantee
  is sampled in practice, not proved.
- **The cost of cancellation-freeness, which we must not ignore:** *"for circuits
  with many voltage-controlled current sources (i.e. those circuits composed of
  transistors), the number of terms that are valid for the intersection problem
  but are not spanning trees of the current graph grows exponentially with the
  circuit size."* Transistor circuits are the bad case for the cheap matroid.

**Reis & Stykel, "Stability analysis and model order reduction of coupled
systems," MCMDS 13(5):413-436, 2007.** Abstract only.
- The **only positive composition result found anywhere**: reduce subsystems,
  recouple through the original interconnection, and *"obtain error bounds for the
  reduced-order closed-loop system in terms of the errors in the reduced-order
  subsystems."* Norm-based with a stability side condition, so a per-block *term
  tolerance* is not the input it takes. Hypotheses unverified.

**Higham, "The Accuracy of Floating Point Summation," SIAM J. Sci. Comput.
14(4):783-799, 1993**; **Ogita, Rump & Oishi, "Accurate Sum and Dot Product,"
SIAM J. Sci. Comput. 26(6):1955-1988, 2005.** Abstracts verified.
- Our `κ` is the standard summation condition number `Σ|x_i| / |Σ x_i|`. Fixes
  *evaluation* (compensated summation buys back the ~4 digits lost at
  `κ = 9.4e3`); does nothing for *truncation*. Two separate problems.

**Filaretov, Gorshkov & Kurganov, "A Cancellation-Free Symbolic Sensitivity
Technique Based on Network Determinant Expansion," Advances in Electrical
Engineering 2015:328517.**
- (open) http://downloads.hindawi.com/archive/2015/328517.pdf
- A non-DDD, non-two-graph route to cancellation-freeness via *"high order
  summative cofactors and the generalized parameter extraction method"*. Their
  2018 companion identifies cancelling terms with *"determinants of the circuit
  with singular network elements - norator or nullator"* — a structural
  characterisation testable against our composite MNA entries.

**Not characterised, listed so they are not re-derived from titles:** Daems,
Verhaegen, Wambacq, Gielen & Sansen, TCAS-I 46(5):594-606, 1999 (abstract only;
the field's central error-control comparison, scoped to *flat* analysis);
**Fernández, Rodríguez-Vázquez, Martín & Huertas, AICSP 3(1):43-58, 1993**
(abstract only — **the one paper claiming an algorithm for approximating
*hierarchical* formulae; highest-value gap**); Guerra et al., AICSP 31:131-145,
2002 (abstract only, journal version of the DATE 2000 entry); Sommer, Hennig,
Dröge & Horneber, *Alta Frequenza* 5(6), 1993 (**TITLE ONLY**, not digitised);
Verhaegen & Gielen, AICSP 31:119-130, 2002 (**TITLE ONLY**); Lasota, ICSES 2012
(abstract only — independent claim that a multilevel decomposition stays
cancellation-free *and* that this makes term elimination sound); Hashemian,
IEEE TCAS-I 69(7), 2022 (abstract only — claims the determinant's *magnitude*
factors into an all-passive spanning-tree product with signs isolated in a nullor
part, which would give a `κ ≈ 1` expansion if it applies to VCCS-heavy models).

## G. Cancellation-free DDD and topological simplification (added 2026-07-30)

The determinant-side answers to term cancellation. **The thread is Sheldon X.-D.
Tan's (UCR), not Guoyong Shi's** — three researchers are routinely conflated:
Guoyong Shi (SJTU, GPDD), C.-J. Richard Shi (UW, original DDD), Sheldon X.-D. Tan
(UCR). Bearing on our work: `doc/cancellation_ranking_conclusions.md` §15.

**Tan, Qi & Li, "Hierarchical Modeling and Simulation of Large Analog Circuits,"
DATE 2004.**
- (open) https://intra.ece.ucr.edu/~stan/papers/date04.pdf
- **Theorem 2 gives a purely combinatorial test for cancelling terms** — no
  numerics: *"For a given product term from a determinant, which consists k
  first-order cofactor ... k >= 2, if there are two first-order cofactors that
  share the same row index or column index, then there exists another product
  term which will cancel with this product term."* Root cause: *"Term cancellation
  ... will happen when MNA formulation is used where each device admittance may
  appear more than once."* Directly implementable against our matrix.

**Tan, Guo & Qi, "Hierarchical Approach to Exact Symbolic Analysis of Large Analog
Circuits," DAC 2004, pp. 860-863.**
- (open) https://intra.engr.ucr.edu/~stan/papers/dac04_1.pdf
- Theorem 1 + partial-DDD/complementary-DDD construction gives a **flat**
  cancellation-free DDD from a hierarchical Schur decomposition. µA741 and µA725
  (the latter exactly, for the first time). Their stated cost: *"the new
  construction method will lead to larger DDD size than non-hierarchical method in
  general."*
- **Contradicts Song & Shi (ASP-DAC 2012)**, who motivate hierarchical GPDD by
  asserting hierarchical DDD lacks the cancellation-free property, without citing
  this. Guoyong Shi's own survey lists "DDD -> De-cancellation, Exact" citing this
  paper. Unresolved; treat both claims with the corresponding caution.

**Tan & Shi (C.-J. Richard), "Efficient Approximation of Symbolic Expressions...,"
IEEE TCAD 23(6), 2004** — also listed in §F.
- (open) https://intra.engr.ucr.edu/~stan/papers/tcad04a.pdf
- **Contains a counter-example to our headline claim.** For a two-stage op-amp's
  `s^1` denominator coefficient: *"the first product term amounts to 86% of the
  total magnitude of the coefficient and the first two terms amount to 97%."*
  Their pipeline de-cancels first, ranks **per `s^k` coefficient** on a multi-root
  DDD, and uses **individual circuit parameters** rather than composite entries.
- The price of de-cancellation, stated: *"cancellation-free s-expanded DDDs do not
  satisfy Theorem 1"* — it **destroys DDD canonicity**, which the efficient graph
  operations rely on. They therefore generate terms *before* de-cancelling.

**Hu, Shi, Tai & Lee, "Topological Symbolic Simplification for Analog Design,"
ISCAS 2015.**
- (open) http://hanbinhu.github.io/data/paper/2015_ISCAS_TopoSimp.pdf
- **The idea that sidesteps `κ` completely: rank circuit *elements*, not terms.**
  Each symbol is tried Short and Open on the GPDD, the reduced transfer function
  is evaluated **exactly**, scored by RMS relative error on dc gain and phase
  margin, and the worst-scoring K symbols removed. Folded-cascode op-amp: 123
  symbols to 18 in 3.9 s, and the reduced dc-gain expressions **match the textbook
  hand-derived formulas**. Their claim, to be tested not believed: *"The property
  described so far for element operation is a unique property owned by GPDD only;
  DDD does not have an analogous property."*

**Hao, Tan, Shen & Shi, "Performance bound analysis of analog circuits considering
process variations," DAC 2011, pp. 310-315.** Abstract only.
- *"We show that symbolic de-cancellation is critical for the affine interval
  analysis."* Confirms that when cancellation genuinely must go, this group used
  the **determinant-side** fix rather than switching to a graph model. Their own
  caveat: the Kharitonov bounds *"are conservative given the correlations among
  coefficient intervals"*.

**Guoyong Shi, AICSP survey (already held) — the passage that matters.**
- (open, Wayback) https://web.archive.org/web/2024id_/https://aice.sjtu.edu.cn/msda/data/publication/AICSP_survey_2011_gshi.pdf
- Defines the term narrowly: *"the problem of term cancellation, i.e., the
  existence of identical product terms with opposite term signs."* And judges:
  *"Numerically speaking, term cancellation might not be a very serious problem.
  But some applications do require the elimination of term cancellation, such as
  in the application of variational analysis based on interval algebra."* So his
  "cancellation" is **exact** cancellation, a different quantity from our `κ`.

**Shi & Meng, "Variational Analog IC Design via Symbolic Sensitivity Analysis,"
ISCAS 2009.**
- (open, Wayback) https://web.archive.org/web/2024id_/https://aice.sjtu.edu.cn/msda/data/publication/ISCAS2009_shi_meng.pdf
- The cleanest statement of *which* property buys cancellation-freeness:
  *"matrix-free computation, cancellation-free sum-of-product terms, maintaining
  one-to-one correspondence from the circuit parameters to the simulator symbols,
  and direct circuit-topology based construction."* The third clause — no
  composite symbols — is the load-bearing one. But see §15c: one-symbol-per-device
  *without* de-cancellation provably makes `κ` worse.

**Filaretov & Gorshkov, "Efficient generation of compact symbolic network
functions in a nested rational form," Int. J. Circuit Theory Appl. 48(7), 2020.**
Abstract only. DOI 10.1002/cta.2789
- Cancellation-free **determinant** expansion via generalized parameter extraction,
  output as *"a fully symbolic form of rational expressions or a nested s-expanded
  polynomial"*. A nested form is the other structural answer: never form the
  expanded sum whose terms cancel. Claim: *"CirSym is the only available program
  that provides the exact calculation of the symbolic function of large circuits
  in the s-expanded form with every coefficient being a compact-nested
  expression."*

**Shi, Tan & Tlelo-Cuautle, *Advanced Symbolic Analysis for VLSI Systems*,
Springer, 2014.** DOI 10.1007/978-1-4939-1103-5
- The authoritative monograph; Shi's papers cite it for GPDD implementation
  details. Chapters 6 (Generalized Two-Graph Theory) and 8 (Hierarchical Analysis
  Methods) are the likely home of an explicit hierarchical-cancellation treatment.
  **Only the chapter list and ch. 7's abstract were obtained — highest-value gap.**
- Observation worth recording: **there is no chapter on approximation,
  simplification or term ranking.** For the definitive 2014 treatment of this
  family that absence is itself informative — the approximation work postdates it
  or sits outside it.

**Not characterised, listed so they are not re-derived from titles:** Yang, Ranjan,
Verhaegen, Ding, Vemuri & Gielen, "Efficient symbolic sensitivity analysis ...
using **element-coefficient diagrams**," ASP-DAC 2005 (**TITLE ONLY** — but by its
name a DDD variant whose symbols are circuit *elements* rather than matrix
entries, which per §15c is exactly the missing half of the de-cancellation story;
**highest-value unverified lead**); Verhaegen & Gielen, TCAS-II 49(7), 2002
(**TITLE ONLY** — their rival DDD de-cancellation by vertex duplication is known
to us only through Tan & Shi's description of it); Tan, ICCAD 2003 and IEEE TCAD
24(3), 2005 (**TITLE ONLY**); Filaretov & Korotkov, "Generalized Parameter
Extraction Method," 2003 (**TITLE ONLY**, founding paper of the CirSym method);
Guoyong Shi, "Toward automated reasoning for analog IC design by symbolic
computation - A survey," Integration 60, 2018 (abstract not retrieved; supersedes
the AICSP survey we hold and is worth acquiring).
