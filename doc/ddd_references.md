# DDD (Determinant Decision Diagram) — reference papers

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

## B. Approximate / dominant-term analysis — the "readable formula for a big circuit" payoff

**Tan & Shi, DDD-based generation of s-expanded / dominant-term symbolic
functions** (behavioral-modeling line; shortest-path & dynamic-programming term
generation, linear in the number of DDD vertices).
- (may need access) Overview: https://link.springer.com/content/pdf/10.1023/A:1015041927036.pdf
- (may need access) s-expanded functions: https://www.academia.edu/77305832/Compact_representation_and_efficient_generation_of_s_expanded_symbolic_network_functions_for_computer_aided_analog_circuit_design
- Why: the algorithm behind `DDDResult.approximate(tol)` — extracting the few
  dominant product terms directly on the graph in linear time. Turns a big
  circuit's exact-but-huge function into a short, interpretable expression.

**Yu & Sechen, "A unified approach to the approximate symbolic analysis of large
analog integrated circuits," IEEE TCAS-I 43(8), 1996.**
- (may need access) https://www.researchgate.net/publication/3224818
- Why: classic "approximate-during-computation" strategy and, crucially, using
  component *variation ranges* rather than a single nominal value when deciding
  which terms to drop. Informs *how* we prune so the simplification stays valid
  across the design space.

## C. Hierarchical DDD — scaling to very large circuits

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

**"A New Approach to Semi-Symbolic Analysis of Analog Integrated Circuits," DAC
2000.**
- (open) https://www.cs.york.ac.uk/rts/docs/SIGDA-Compendium-1994-2004/papers/2000/dac00/pdffiles/01_5.pdf
- Why: *semi-symbolic* = keep some elements numeric, others symbolic — precisely
  pycircuit's common regime (numeric components + symbolic `s`, or a few symbolic
  parameters). Informs which entries to keep symbolic in the DDD.

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

## F. How current is this field?

Searched (Jul 2026) for newer DDD work, post-2015 extensions, ML-based symbolic
analysis and recent PhD theses. Conclusion: **the field peaked 1997–2013 and the
list above is essentially complete.** The most recent substantive items are Shi's
2010–2013 SJTU papers (A2, D2), the 2014 Springer book (E), and the 2025 DSSA
paper (B) which is outside the BDD line altogether. No newer thesis or DDD
successor surfaced. This is stable ground — implement against A2, not against a
moving target.

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
