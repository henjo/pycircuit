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

## 3. Architectural integration — what changes outside DDD

Ground rule 3 says "additive only". That is achievable, but **not** by the obvious
route, and the places where DDD presses on the rest of pycircuit should be decided
deliberately rather than discovered at P2.

### 3.1 What the existing architecture already gives us (verified)

- **The import graph is clean and hierarchical**, which is the main thing that
  makes this containable. `transferfunction.py` imports only numpy and sympy;
  `toolkit.py` imports numpy, sympy and the `_*` backend modules — **neither knows
  about `circuit` or `analysis`**. `analysis_ss.py` sits on top and imports both.
  - **Rule that follows:** `ddd.py` stays a **leaf** — sympy only, no imports from
    `circuit`, `analysis` or `toolkit`. `dddresult.py` may import
    `transferfunction`. Then `toolkit → dddresult → transferfunction` stays
    acyclic. Violating this is the fastest way to create an import cycle in a
    package that currently has none in this area.
- **Noise needs no external change at all.** `noise_psd` is already a method on the
  `Toolkit` base (with `SymbolicPolyToolkit` overriding it), so **P4 is a pure
  override** — zero lines outside DDD. This extension point already exists and was
  built for exactly this.
- **Stamping needs nothing.** DDD consumes `G`, `C`, `u` as stamped; `_symbolic.py`
  is untouched, and no new backend primitives are required. Axis 1 really is
  unaffected.

### 3.2 Where DDD genuinely presses outside itself

| Site | What is needed | Risk |
|---|---|---|
| `analysis_ss.py` AC path (~line 182) | A route that does **not** flatten to sympy | **High** — see 3.3 |
| `toolkit.py` | `DDDToolkit` + `ddd_toolkit` singleton | Low — purely additive, follows the `ginac_toolkit` pattern |
| `CircuitResultACPoly` | Not reusable as-is — see 3.3 | Medium |
| `doc/src/conf.py` | add `sphinx.ext.graphviz` | Trivial, doc-build only |
| `SubCircuit` (P6 only) | Hierarchical stamps | Deferred; the deepest leak, and a reason to keep P6 last |

### 3.3 The bloat problem, and the fix

**The smell.** `analysis_ss.py` currently does capability-query plus inline
construction:

```python
if self.toolkit.supports('num_den') and not isiterable(ss):
    num, den = self.toolkit.linearsolver_num_den(s*C + G, -u)
    xac = np.array([ni / den for ni in num], dtype=object)
    self.result = CircuitResultACPoly(...)
```

Every new representation adds another capability string, another branch, and
another result class. Two backends is fine; this is precisely how analysis code
rots. DDD would make it three, and P3's semi-symbolic variant four.

`CircuitResultACPoly` cannot simply be reused, either: `tf_i` performs real
arithmetic on the numerator vector (`0 * N`, `s * N`, `A + c0*(D - 1)`), which
presumes sympy expressions. Graphs would need every one of those operations.

**The naive fix does not work.** The instinct is a factory on the toolkit —
`toolkit.make_ac_result(...)`. But result classes live in `analysis_ss.py`, which
imports `toolkit.py`; having `toolkit` import them creates a **cycle**. Worth
stating explicitly because it is the first thing anyone will try.

**Recommended fix — delegation, not inheritance or factories.** Have the toolkit
return a representation-agnostic **solution object** with a uniform interface
(`tf(plus, minus)`, `tf_i(...)`, `poles()`, `eval(params)`), and let
`CircuitResultACPoly` *hold* one and delegate to it:

- `symbolic_poly` returns a `NumDenSolution` wrapping `(num, den, s)` — today's
  behaviour, moved not changed;
- `ddd_toolkit` returns a `DDDSolution` wrapping the graph, which never flattens;
- `analysis_ss.py` keeps **one** path, no capability branch, no `dtype=object`
  array on the graph route;
- the solution classes live beside the toolkits, so nothing imports upward.

Net effect: this **removes** the existing conditional rather than adding to it, and
each future backend costs zero lines in shared analysis code.

**Honest caveat.** This is a change to shared code, which sits against ground rule
3. Reconcile it as follows: it is a **behaviour-preserving refactor**, done as its
own commit **before P2**, with the existing AC/`symbolic_poly` tests as the safety
net (they already cover both branches). If it cannot be made green without
touching behaviour, abandon it and accept one extra branch — the refactor is an
investment against future bloat, not a prerequisite for DDD.

### 3.4 Are the other issues the same shape? Mostly not

**DECIDED (2026-07-27): the 3.3 delegation fix is accepted.** Before generalising
it, note that the remaining issues are *different* problems, and applying the same
remedy to all of them would be wrong.

The distinction that matters: **delegation answers "who decides"; duck typing
answers "what type flows through".** The AC dispatch is a *who decides* problem —
shared code was making the choice and constructing the result — so no amount of
duck typing removes that branch. Conflating the two is the usual mistake here.

| Issue | Same shape? | Right fix |
|---|---|---|
| AC dispatch `supports('num_den')` + inline construction | — this is the one | Delegation (3.3) ✔ accepted |
| `supports()` capability strings generally | **Yes**, same family | Disappears with 3.3 — no separate work |
| `tf_i` arithmetic on `N` (`0*N`, `s*N`, `A + c0*(D-1)`) | **No** — not a branch, an implicit *operator contract* | **Duck typing** — see 3.5 |
| `Toolkit.__getattr__` forwarding to the backend | **No** — this is *already* unbounded duck typing, and it is the cautionary case | Do not extend it; DDD methods explicit on `DDDToolkit` |
| `noise_psd` overridable on `Toolkit` base | **No** — already correct polymorphism | Nothing. Use it as the model |
| `default_toolkit` global | **No** — global mutable state | Out of scope; do not entrench it further |
| `isinstance(self.inputsrc, VS/IS)` (`analysis_ss.py:420,426,466,470`) | **Yes**, same anti-pattern family | **Out of scope.** Noted only as evidence the pattern recurs; fixing it is not DDD's job |

### 3.5 Would duck typing help? Yes, in one place; no in another

**Where it helps — and it is already how the code works.** `tf_i` does
`0 * N`, `s * N`, `A + c0*(D - 1)` on the numerator vector. That is not a type
check, it is an operator protocol, and it *already* serves both numpy arrays
(`numeric`) and sympy expressions (`symbolic`). If the DDD numerator objects
implement `__mul__` and `__add__`, **`tf_i` needs no change at all**.

This is a better fit than it first looks: DDDs form an algebra, and the operations
required are the standard ones — the 2008 `ddd.py` already sketched `__mul__`,
`__sub__` and `union` for exactly this reason. So the protocol DDD must satisfy is
one it naturally has.

*Caveat, and it is a real one:* duck typing here hides **cost**, not correctness. A
DDD will happily accept `s * N` and quietly build a larger graph;
`A + c0*(D - 1)` combines two graphs. Nothing errors — it just gets expensive. So
`DDDSolution.tf_i` should be measured early, and may need its own implementation
rather than inheriting the generic arithmetic. Silent expense is harder to notice
than a silent wrong answer.

**Where it does not help — the solution object itself.** With several toolkits and
multiple solution flavours, a contract that exists only implicitly is a
maintenance hazard, and this codebase has already been bitten by exactly that:
`Toolkit.__getattr__` forwards any unknown attribute to the backend module, which
makes `hasattr()` unreliable, turns typos into confusing errors, and means
anything dropped into `_symbolic.py` silently widens every symbolic toolkit's
surface. That is unbounded duck typing, and it is a argument *against* reflexively
adding more.

**Recommendation — follow the precedent already in this repo.** `integrator.py`
defines its strategy family (`Integrator`, and likewise `StepController`,
`NonLinearSolver`) as an **`ABC` with `@abstractmethod`** and type hints. That is
the codebase's own, most recent answer for a pluggable strategy family, and the
solution objects are the same kind of thing. Concretely:

- a small `ACSolution` ABC with abstract `tf`, `tf_i`, `poles`, `eval`;
- `NumDenSolution` and `DDDSolution` inherit it — free, since we write both and
  there are no third-party implementers;
- benefit over pure duck typing: a missing method fails at construction with a
  clear message, instead of an `AttributeError` surfacing deep inside a frequency
  sweep.

`typing.Protocol` would be the alternative (structural, no inheritance), but it
buys flexibility we do not need and departs from the existing convention. Prefer
consistency with `integrator.py`.

**The cheap guard that makes this safe: one parametrised conformance test.** A
single test parametrised over `[NumDenSolution, DDDSolution]` asserting the whole
contract — including the operator behaviour that stays duck-typed — costs almost
nothing and is what keeps two implementations from drifting. This is the piece
that makes duck typing tolerable at this scale.

### 3.6 Dependency discipline

- **No new runtime dependencies.** `ddd_toolkit` must import with sympy alone.
  Follow the established optional pattern (`ginac_toolkit = None` on `ImportError`)
  only if something optional is ever needed — for pure-Python DDD, nothing is.
- **Graphviz is a documentation-build dependency, not a runtime one.** `to_dot()`
  returns a *string*; nothing in the library shells out to `dot`. Keep it that way,
  and keep `dot` out of `install_requires`.
- **Benchmarks: split by audience.** The circuit *builders* are reusable fixtures
  wanted by both tests and the doc build, so they belong in the package. The
  *harness* (timing, memory, `ulimit`, timeout handling) is developer tooling and
  should not ship inside the library — put it under `doc/` (already on the doc
  build path via `conf.py`) or a top-level `benchmarks/`.
- **Do not add DDD primitives to `_symbolic.py`.** `Toolkit.__getattr__` silently
  forwards unknown attributes to the backend module, so anything placed there
  becomes reachable from *every* symbolic toolkit and quietly widens their surface.
  DDD methods go on `DDDToolkit` explicitly.

## 4. File map

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

## 5. Stages

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
symbol-convention note (see Stage P0).

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
- **(a) Capability** — on the N = 8 ladder, where *[OURS]* `to_ratio` /
  `poly_coeffs` do not complete: construction finishes within the benchmark
  timeout, the graph evaluates numerically to the same transfer function as
  `numeric` (to `1e-10` relative), and `|DDD|` is **under 10 000 vertices** so that
  P1's `q·d·|DDD|` expansion stays tractable. Concrete pass/fail, not a judgement.
- **(b) Size** — `|DDD| ≤ 2×` SoE ops at N = 12, growing no faster across N = 4→16.
  **Caveat on (b): this is not yet an apples-to-apples comparison.** Under the
  compact-symbol convention a DDD vertex carries a whole matrix entry
  (`g₁+g₂+s(c₁+c₂)`) while an SoE operation is over component-level symbols — the
  DDD has fewer, fatter units. Report the per-vertex payload sizes alongside, and
  if the two cannot be reconciled, treat (a) as decisive and (b) as indicative.

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

### Stage P2 — `DDDResult` + toolkit integration *(~3 days)*

**Read this before planning P2 in detail — the obvious integration is a trap.**
`SymbolicPolyToolkit.linearsolver_num_den` returns *sympy expressions*, and the
whole downstream path is eager about it: `analysis_ss.py:184` immediately builds
`xac = np.array([ni / den for ni in num], dtype=object)`, and
`TransferFunction.as_num_den` calls `sympy.fraction(self.canonical())`. So a
`DDDToolkit` that satisfies the existing contract **must flatten its DDDs into
sympy expressions — which is exactly the exponential blow-up the DDD exists to
avoid.** Flattening works only for circuits small enough that `symbolic_poly`
already handles them, i.e. not the motivating case.

Therefore P2 has two deliverables, and the second is **not** optional:

1. **Compatibility path** (`linearsolver_num_den` → sympy). Keep it, because it
   makes every existing analysis work unchanged and is the cheapest correctness
   check against `symbolic_poly`. But label it honestly: it is valid for **small
   circuits only**, and must carry a size guard that refuses (or warns and falls
   back) above a vertex/term threshold rather than hanging — the same discipline as
   `ginac_max_dim` and `MAX_COMPILE_CHARS`.
2. **Graph-preserving path** — `DDDResult` as the primary interface, never
   flattened: `denominator()`, `numerator(i)`, `tf(i,j)` and `eval(params)` all
   operating on the graph, with numeric evaluation done by walking it. This is
   where the value actually is, and deferring it (as an earlier draft of this plan
   did) would leave the project unable to serve its own use case.

Tasks: `DDDResult`; `DDDToolkit(SymbolicPolyToolkit)` with both paths;
`ddd_toolkit` singleton; `supports('ddd')` plus the `CircuitResultACPoly`
subclass needed to reach the graph path from `AC`.

**Integration details to check early** (each has bitten a previous backend):
- `self.toolkit.concatenate` / `array` at `analysis_ss.py:186` operate on the
  numerator vector; confirm they behave for whatever objects the graph path
  returns, or arrange for the graph path to bypass that branch entirely.
- interaction with the deprecated `default_toolkit` global;
- sympy expressions as vertex payloads hash structurally — `Add`/`Mul` are
  canonicalised so `g1+g2` and `g2+g1` collide correctly, but this should have a
  test rather than being assumed.

**Gate:** `AC(cir, toolkit=ddd_toolkit)` agrees with `symbolic_poly` on every
circuit where the latter runs; the graph path returns correct numeric results on a
circuit where the compatibility path cannot complete; **full existing test suite
green**.

**Docs:** cofactors → Cramer → transfer functions, the two paths and why both
exist, and the user-facing API.

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

## 6. Definition of done

**Per stage** — all four, or the stage is not done:
1. code + tests, full suite green;
2. its theory-document section written, **in the same commit as the code**;
3. `make html` builds clean — 0 warnings, 0 errors, no `exec-rst` fallbacks (a
   fallback means a live example failed and the doc is no longer verifying);
4. its gate evaluated and the result recorded — including a negative one.

**Overall.** An earlier draft said "symbolic poles and zeros of a fully-symbolic
~10–15 node circuit". **That target is impossible in principle and has been
corrected.** Two independent reasons:

- *Abel–Ruffini.* A 10–15 node circuit has a denominator of degree ~10–15, and a
  general polynomial of degree ≥5 has no closed-form roots. Verified concretely:
  `sympy.roots` on a general degree-5 polynomial with symbolic coefficients
  returns `{}`. No representation — DDD or otherwise — can produce what does not
  exist.
- *Size.* Even the exact symbolic *coefficients* of such a circuit stand for
  ~10³⁰⁺ product terms (`ddd_conclusions.md` §4.3). They can be represented
  compactly and evaluated, but never printed or read.

The corrected target, which is what the literature actually delivers:

1. **Exact s-expanded coefficients**, held compactly as a multiroot DDD, for a
   fully-symbolic ~10–15 node circuit, in seconds — never flattened, and correct
   when evaluated numerically against `symbolic_poly`/`numeric` where those can
   run.
2. **Dominant pole/zero estimates** as ratios of coefficients of consecutive
   powers of `s` — *[LIT]* TCAD 2001 §I, and the only route to interpretable
   poles for a circuit this size.
3. **Exact numeric poles** via `numpy.roots` once parameters are substituted.
4. **A short, readable symbolic expression** only after approximation prunes to
   dominant terms.

**Consequence for sequencing:** items 2 and 4 are P5 work, so **P5 is on the
critical path to any human-readable result for large circuits** — it is not the
optional finishing touch that §5's ordering implies. P1–P4 produce something
correct, compact and evaluable but not yet *readable*. If the point of the
exercise is designer insight, P5 must be scheduled, not merely listed.

---

## 7. Sequencing notes

- Stage B strictly first — its baselines must be recorded *before* DDD exists, or
  the comparison is retrospective.
- P0's symbol-convention decision precedes P0 coding; it defines the hash key.
- P1 depends on P0; P2 on P1 (the result object exposes coefficients). P3, P4 and
  P6 are independent of each other and can be reordered by value once P2 lands.
  **P5 is not in that set** — per §6 it is the only route to human-readable output
  for large circuits, so it should be scheduled directly after P2 unless the
  semi-symbolic regime (P3) is the more urgent user need.
- **Python-level practicalities to settle in P0, not discover in P4**: recursion
  depth (LED expansion is naturally recursive and will exceed
  `sys.setrecursionlimit` defaults on realistic matrices — write it iteratively or
  raise the limit deliberately), and the cost of sympy expression payloads per
  vertex, which dominates memory (§Stage B).
- Paper calibration circuits (µA741/µA725, Tier 3) come *after* P0 passes. They
  need a hybrid-π `SubCircuit` — `VCCS` + `R` + `C`, no new element types — and
  should not gate P0, since a µA741 is exactly where a first implementation
  struggles for uninteresting reasons.
- Keep paper PDFs out of the repo; they live in `~/pycircuit_agy/papers/ddd/`.
