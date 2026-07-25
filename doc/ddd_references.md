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

## D. Foundational data structure

**Bryant, "Graph-Based Algorithms for Boolean Function Manipulation," IEEE TC
35(8), 1986.**
- (open) https://www.cs.cmu.edu/~bryant/pubdir/ieeetc86.pdf
- Why: DDD is a BDD variant. reduce/apply, the unique-node hash table (how
  sharing is *enforced*), canonical form, and — most important for us —
  **variable ordering**, all transfer directly. The reference for making the
  diagram actually compact and for the ordering heuristic (the main risk).

## E. Context / survey and the semi-symbolic case

**"A survey on binary decision diagram approaches to symbolic analysis of analog
integrated circuits," Analog Integr. Circ. Sig. Process., 2011.**
- (may need access) https://link.springer.com/article/10.1007/s10470-011-9773-8
- Why: a map of the whole BDD-for-symbolic-analysis field; places DDD relative to
  alternatives before we commit.

**"A New Approach to Semi-Symbolic Analysis of Analog Integrated Circuits," DAC
2000.**
- (open) https://www.cs.york.ac.uk/rts/docs/SIGDA-Compendium-1994-2004/papers/2000/dac00/pdffiles/01_5.pdf
- Why: *semi-symbolic* = keep some elements numeric, others symbolic — precisely
  pycircuit's common regime (numeric components + symbolic `s`, or a few symbolic
  parameters). Informs which entries to keep symbolic in the DDD.

## Author PDF directories (many more open PDFs)

- Sheldon X.-D. Tan (UC Riverside): https://intra.engr.ucr.edu/~stan/papers/
- C.-J. Richard Shi, MSCAD lab (U. Washington): https://labs.ece.uw.edu/mscad/shi/Papers/

---

Note: the papers already in `~/pycircuit_agy/papers/` (Fang DAC'13, Yao–Wang
ICECS'14) concern *transient / LTE step control*, unrelated to DDD — this is a
separate reading set.
