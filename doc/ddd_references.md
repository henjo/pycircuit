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
- Why: the hierarchical counterpart on the GPDD side (see D2), keeping the
  cancellation-free property at every level of hierarchy. Relevant only if we take
  the graph route; listed for completeness.

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
  DDD. The paper's own conclusion is the useful part for us: GPDD runtime/memory
  is only *comparable* to DDD despite generating many more terms. That is a real
  argument for staying with DDD (which also fits pycircuit better, since we
  already have an MNA matrix and no graph model of the circuit).

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
