# DDD implementation plan

**Inputs:** `ddd_references.md` (what to read) and `ddd_conclusions.md` (what we
concluded and why, including the gates). This document says *what to build, in
what order, and when to stop*. Where the three disagree, `ddd_conclusions.md`
wins on rationale and this file wins on sequencing.

Written 2026-07-27. No DDD code exists yet.

---

## 1. Ground rules

These are carried from `ddd_conclusions.md` and are not re-litigated during
implementation.

1. **A new analysis, not a new tool** (§2.2). Circuits are built with
   `Circuit`/`SubCircuit` and the existing elements, stamped through the existing
   MNA path, and reached as `AC(cir, toolkit=ddd_toolkit)`. No parallel front end,
   no separate netlist, no hand-built matrices as the interface.
2. **`ddd.py` is not a constraint** (§2.1). It is rewritten in place — it has no
   importers other than its own test (verified).
3. **Nothing existing changes behaviour.** Additive only; `numeric`, `symbolic`,
   `symbolic_poly` untouched. The work must stay revertible by deleting modules.
4. **Gates are real.** A failed gate stops or redirects the work; it does not get
   renegotiated after the fact. Test circuits are fixed in Stage B, before results
   exist.
5. **Theory documentation is written with the code, in the same commit** — see §2.
   A stage is not complete without it.

---

## 2. The theory document — written alongside, not afterwards

**Why this is a rule and not a nicety.** A document that explains the algorithm
serves three purposes here, and only the first is conventional:

- it explains the implementation to a reader (including us, later);
- it is an **anchor**: writing the explanation forces the design to be explainable,
  and vague designs fail that test early, while they are still cheap to change;
- it is **executable verification**. Built the way described below, the examples
  and diagrams are produced by *running the real implementation at doc-build time*.
  They cannot drift from the code, because if the code changes incompatibly the
  doc build breaks.

That last property is why this must be concurrent. A theory document written after
the fact describes what was built; one written alongside constrains what gets
built, and then keeps testing it.

### Mechanism

- **Location:** `doc/src/circuit/ddd.rst`, in the circuit toctree.
- **Live examples:** the existing `.. exec-rst::` directive
  (`doc/sphinxext/exec_directive.py`), already used by `soe_symbolic.rst` and
  `ginac_native.rst` — runs Python at build time and inserts its printed reST,
  with a timeout and graceful source-fallback.
- **Diagrams:** `sphinx.ext.graphviz` — **must be added to `extensions` in
  `doc/src/conf.py`** (stock Sphinx extension; graphviz 14.1.2 is installed,
  `/usr/bin/dot`). The implementation provides `to_dot()`, and an `exec-rst` block
  *prints a `.. graphviz::` directive containing DOT generated from a real DDD*.
  So every diagram in the document is a picture of what the code actually built.
  - The old `ddd.py` used `yapgvb` for this. Do not revive it — emit DOT text and
    let Sphinx run `dot`. No Python dependency.
  - **Size guard, mandatory.** A DDD of a real circuit has thousands of vertices
    and must never reach `dot`. `to_dot()` takes a vertex limit and raises above
    it; all rendered examples are deliberately tiny (a 3×3 matrix, a 2–3 section
    ladder). This mirrors the `MAX_COMPILE_CHARS` guard added after the GiNaC
    `compile_ex` OOM incident.
- **Measured numbers in prose are generated**, never typed. Any figure quoted in
  the text comes from an `exec-rst` block, so it cannot go stale.

### Contents, and which stage writes each part

| Section | Written during | Content |
|---|---|---|
| Why symbolic determinants explode | P0 | The problem, in one worked example: `n!` terms, and where sharing appears |
| Laplace expansion → DDD | P0 | The 3×3 example from ICCAD 2010 Fig. 1, **rendered**; 1-paths as product terms |
| Layered expansion (LED) | P0 | Queues layer by layer, minor hash keys, the LED→DDD conversion, **rendered before/after** |
| Sharing, concretely | P0 | The same minor reached by two paths, **rendered**, with the hash key that unified them |
| Signs | P0 | Why they fall out of construction; the `det` cross-check |
| Ordering and complexity | P0 | Min-degree heuristic; the `n·2^(n-1)` identity **plotted: measured vs theory** |
| Cofactors, Cramer, transfer functions | P2 | From determinant to `v(out)/v(in)` |
| s-expansion / multiroot | P1 | One root per power of `s`, **rendered showing shared subgraphs**; measured linearity |
| Numeric terminals | P3 | Semi-symbolic; why the coefficient blow-up cannot form |
| How it compares | P0→ | Live table: DDD vs SoE vs `symbolic_poly`, regenerated at build |

---

## 3. File map

| Path | Status | Contents |
|---|---|---|
| `pycircuit/circuit/ddd.py` | **rewritten** | `DDD` vertex/graph, `ddd_of_matrix` (LED), minor hash table, `eval`, `cofactor`, `to_dot` |
| `pycircuit/circuit/dddresult.py` | new | `DDDResult` — lazy Cramer numerators, `tf`, `poles`, `denominator`, `eval` |
| `pycircuit/circuit/toolkit.py` | edited | `DDDToolkit(SymbolicPolyToolkit)` + `ddd_toolkit` singleton |
| `pycircuit/circuit/symbolic_benchmark.py` | new | Benchmark circuits + runner (size, time, memory, completion) |
| `pycircuit/circuit/tests/test_ddd.py` | rewritten | Core: identities, signs, `det` cross-checks, randomised |
| `pycircuit/circuit/tests/test_dddresult.py` | new | Result object + toolkit integration |
| `doc/src/circuit/ddd.rst` | new | The theory document (§2) |
| `doc/src/conf.py` | edited | add `sphinx.ext.graphviz` |
| `doc/src/circuit/index.rst` | edited | toctree entry |

---

## 4. Stages

Effort figures are rough sizing, not estimates from a decomposition.

### Stage B — Benchmark harness *(~1 day)*

**Goal:** make every later gate executable, and fix the comparison before results
exist.

Tasks:
1. `symbolic_benchmark.py` with the **circuit builders** (all real pycircuit
   circuits): fully-symbolic RC ladder (N = 4, 8, 12, 16, matching
   `soe_symbolic.rst` exactly, same `v_out/v_in`); the Example-10 MFB filter (has
   a `Nullor`); a semi-symbolic opamp-like circuit; and dense symbolic n×n
   matrices for n = 3…8 (matrix-level, for the Tier-1 identity).
2. A runner producing per (circuit, backend): representation size, **raw
   expression size** (the confound), construction time, evaluation time per sweep
   point, **peak memory**, correctness vs a reference, and **completion status** —
   `ok` / `timeout` / `killed`, since several backends hang rather than fail.
3. Every case runs under a timeout and `ulimit -v`.
4. Baselines recorded now, before any DDD exists: `symbolic_poly`, `ginac`, `soe`.

**Gate:** the harness reproduces the published SoE op counts (73/157/241/325) for
the ladder. If it cannot reproduce our own numbers, it cannot referee anything.

**Docs:** the "How we measure" section, including what each metric means and the
symbol-convention note (§5 below).

---

### Stage P0 — Core construction *(largest stage, a few days)*

**Decide first — the symbol convention** (`ddd_conclusions.md` §7). Build the DDD
over *matrix entries* (compact), keeping each entry's symbolic expression as the
vertex payload so component-level results remain available. Record the choice in
the docs beside every number. Validate at P0; do not assume.

Tasks:
1. LED construction: layered expansion, whole row/column at a time; min-degree
   expansion ordering chosen on the fly; **minor hash table keyed by
   `(row-index tuple, col-index tuple)`**; signs determined during construction.
2. LED → DDD conversion (next-layer pointers become 1-edges; 0-edges along sibling
   groups terminating at `0`; bottom queue terminates at `1`).
3. Numeric evaluation over the graph.
4. `to_dot()` with a vertex limit.
5. Drive it from real circuits through the existing MNA path.

Tests:
- **Tier-1 identity:** `|DDD| == n·2^(n-1)` exactly, dense n×n, row-wise, n = 3…8.
  This is the structural check that sharing works.
- **Signs:** numeric evaluation vs `sympy.Matrix.det` on random sparse matrices.
- **Randomised:** dimension 3–8, `det` and Cramer vs `linearsolver_num_den`.
- **Nullor:** MFB result checked against `symbolic_poly`, not merely measured.

**Gate (two criteria — passing *either* continues; both failing stops):**
- **(a) Capability** — a correct flat DDD whose vertex count makes P1 plausible at
  N ≥ 8 on the ladder, where `to_ratio`/`poly_coeffs` currently hang.
- **(b) Size** — `|DDD| ≤ 2×` SoE ops at N = 12, growing no faster across N = 4→16.

**Disambiguate before acting:** the Tier-1 identity must pass first. A poor gate
result with Tier 1 *failing* means our sharing is broken, not that DDD is wrong —
fix the code. Only a gate failure with Tier 1 passing is evidence about the method.

**On failure of both (Tier 1 passing):** write up the negative result in `doc/`
beside the symengine and GiNaC write-ups, leave the module experimental and
unwired, close the roadmap item. **If (b) fails but (a) passes:** continue, but
drop the compactness claim and justify on capability only.

**Docs:** the bulk of the theory document — problem, Laplace→DDD, LED, sharing,
signs, ordering/complexity, all with rendered diagrams and the measured-vs-theory
plot.

---

### Stage P1 — s-expanded / multiroot DDD *(~2 days)*

Tasks: coefficient split, one root per power of `s`, roots sharing subgraphs;
`q`/`d` measurement to test Theorem 1's linearity empirically.

**Gate:**
- coefficients match `symbolic_poly` / `_ginac.poly_coeffs` where both can run;
- size grows linearly in `|DDD|`;
- **pole accuracy**, not merely coefficient agreement — per
  `ddd_conclusions.md` §8.6, root-finding from expanded high-degree coefficients
  is ill-conditioned, and we have already silently dropped poles this way once.

**Docs:** s-expansion section, with a rendered multiroot diagram showing shared
subgraphs, and the measured linearity table.

---

### Stage P2 — `DDDResult` + toolkit integration *(~2 days)*

Tasks: `DDDResult` (lazy Cramer numerators, `tf`, `denominator`, `poles`, `eval`);
`DDDToolkit(SymbolicPolyToolkit)` overriding `linearsolver_num_den` and inheriting
`supports('num_den')`; `ddd_toolkit` singleton.

Conservative integration only — returns sympy `num/den`, so nothing else in the
codebase learns about DDD. The `supports('ddd')` deep path is deferred until there
is a reason for it.

**Gate:** `AC(cir, toolkit=ddd_toolkit)` produces results agreeing with
`symbolic_poly` on every circuit where the latter runs; **the full existing test
suite stays green**.

**Docs:** cofactors → Cramer → transfer functions, and the user-facing API.

---

### Stage P3 — Numeric terminals (semi-symbolic / MTDDD) *(~2 days)*

Terminals carry numeric values; numeric sub-products collapse during construction.

**Gate:** no exact-rational blow-up at dimensions where GiNaC stalled (~16+), and
results still agree with `symbolic_poly` where it runs.

**Docs:** numeric-terminal section — why the blow-up cannot form.

---

### Stage P4 — Noise via a single multiroot DDD *(~2 days)*

All per-noise-source transfer functions in one multiroot DDD.

**Gate:** agrees with existing `noise_psd`, at a cost comparable to a single
transfer function.

---

### Stage P5 — Approximation / dominant terms *(~3 days)*

`approximate(tol)` on the s-expanded form (never before it); dominant poles/zeros
as ratios of consecutive coefficients; prune on component *variation ranges*
rather than a single nominal value.

---

### Stage P6 — Hierarchy *(~3 days)*

Symbolic stamps per `SubCircuit` (ASP-DAC 2011 formulation).

---

### Not scheduled — revisit after P2

**Symbolic sensitivity.** pycircuit has none, and DDD makes it nearly free (the
derivative with respect to a symbol *is* the cofactor). Plausibly the best
value-per-effort item in the programme, but it is a new user-facing feature rather
than part of making DDD work.

---

## 5. Definition of done

**Per stage** — all four, or the stage is not done:
1. code + tests, full suite green;
2. its theory-document section written, **in the same commit as the code**;
3. `make html` builds clean — 0 warnings, 0 errors, no `exec-rst` fallbacks (a
   fallback means a live example failed and the doc is no longer verifying);
4. its gate evaluated and the result recorded — including a negative one.

**Overall** — symbolic poles and zeros of a fully-symbolic ~10–15 node circuit, in
seconds, via `AC(cir, toolkit=ddd_toolkit)`, agreeing with `symbolic_poly`
wherever that can still produce an answer. Today this is impossible at any
non-trivial size.

---

## 6. Sequencing notes

- Stage B strictly first — its baselines must be recorded *before* DDD exists, or
  the comparison is retrospective.
- P0's symbol-convention decision precedes P0 coding; it defines the hash key.
- P1 depends on P0; P2 on P1 (the result object exposes coefficients). P3–P6 are
  independent of each other and can be reordered by value once P2 lands.
- Paper calibration circuits (µA741/µA725, Tier 3) come *after* P0 passes. They
  need a hybrid-π `SubCircuit` — `VCCS` + `R` + `C`, no new element types — and
  should not gate P0, since a µA741 is exactly where a first implementation
  struggles for uninteresting reasons.
- Keep paper PDFs out of the repo; they live in `~/pycircuit_agy/papers/ddd/`.
