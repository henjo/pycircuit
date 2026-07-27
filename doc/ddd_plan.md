# DDD implementation plan

**Inputs:** `ddd_references.md` (what to read) and `ddd_conclusions.md` (what we
concluded and why, including the gates). This document says *what to build, in
what order, and when to stop*. Where the three disagree, `ddd_conclusions.md`
wins on rationale and this file wins on sequencing.

Written 2026-07-27; rebuilt the same day once round one finished. **Stages B and
P0–P6 are complete** — see §4 for what each one actually returned, including the
two that contradicted the plan that produced them. §5 is the remaining work, in
priority order.

Sections 1–3 are carried forward unchanged: the ground rules, the
documentation-with-the-code rule and the architectural analysis all still govern,
and the decisions recorded in §3 (solution-object delegation, the ABC-versus-duck-
typing split, dependency discipline) are settled rather than open. Read them as
constraints on round two, not as history.

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
| Approximation, noise, hierarchy | P5, P4, P6 | Written as those stages landed |
| Calibration against the papers | Tier 1/3 | The exact identity, and the µA741 |
| *Ordering* | H1 | The distribution over node numberings, not a single figure |
| *Recursive hierarchy* | H2 | Replaces the "does not reproduce 56×" note with the measured outcome |
| *Reduced noise* | H3 | How internal sources are referred to terminals |
| *Sensitivity* | H4 | New section |

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
| `SubCircuit` (hierarchy) | Hierarchical stamps | Still open: H1 works on matrices, and wiring it to `SubCircuit` remains the deepest leak |

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

## 4. Status — round one is complete

Stages B and P0–P6 are implemented, measured and documented, together with the
solution-object refactor and calibration tiers 1 and 3.

| Stage | Outcome |
|---|---|
| B — benchmark harness | Gate passed: reproduces the published SoE op counts exactly |
| P0 — LED construction | Both gates passed; `\|DDD\|` 63 vs SoE 466 at N=16, where `symbolic_poly` times out |
| P1 — s-expanded multiroot | All three gates passed, including pole *accuracy* |
| P2 — `ddd_toolkit` | Gate passed; graph path answers where expansion cannot (94 vertices for 2.18M terms) |
| P3 — numeric terminals | **Measured negative.** Implemented, correct, off by default: its premise did not hold |
| P4 — noise | Gate passed; whole transimpedance family for ~3× one transfer function |
| P5 — approximation | Delivered, and conditional: prunes only as far as component values differ |
| P6 — hierarchy | Delivered **single-level**; does *not* reach the published 56× |
| Calibration | Tier 1 exact (`n·2^(n-1)`, n = 1…7); Tier 3 µA741 built and compared |
| H1 — node ordering | **Done.** Gate passed: spread 2.7–2.8× → 1.00× (ladder, MFB) and 1.42× (µA741) |
| H2 — sequential suppression | **Done.** Gate passed: µA741 1040 → 156 vertices, 6.7× |
| H3 — noise through the reduction | **Done.** Gate passed: PSD identical, µA741 noise 11 088 → 26 vertices and 86 s → 0.48 s |
| H4 — symbolic sensitivity | **Done.** Gate passed: matches finite differences to ~1e-8; ten parameters cost what one does |
| H5 — Tier 2 calibration | **Done.** Our DDD-vs-SoE ratios 3.3–5.6× land in the published 5.0–6.0× band |
| H6 — conditioning | **Done, by a different route.** Scaling does not help; extended precision takes N=20 from 5e-10 to 1.6e-15 |

The original definition of done is met: exact s-expanded coefficients held
compactly, dominant pole/zero estimates from coefficient ratios, exact numeric
poles after substitution, and a readable expression after approximation.

Two results contradicted the plan that produced them, and both are recorded as
such rather than quietly dropped — P3's premise (the exact-rational blow-up it
was meant to prevent cannot arise in a representation that never multiplies
entries symbolically) and P5's generality (approximation is conditional on
parameter spread). A third, P6's single-level limitation, is what round two
opens with.

---

## 5. Round two — remaining work, in priority order

Ordering is as directed: hierarchy first, then ordering, then the rest.

### H1 — Expansion and node ordering — **DONE**

**Why first.** Rebuilding the µA741 with a different node *numbering* moved
`|DDD|` from 1040 to 2424 — a 2.3× swing from something arbitrary. Every vertex
count reported anywhere in this work carries that much slack, including the
comparisons against SoE and against the published figures.

That matters most for what comes next. H2's gate is "does recursive suppression
beat the flat diagram", and both sides of that comparison are vertex counts. A
marginal result would be uninterpretable — method or numbering? — so the
measurement has to be made trustworthy before it is relied on.

Tasks:
1. Characterise it: permute node order repeatedly and report the **distribution**
   of `|DDD|`, not a single number.
2. Compare policies — the current on-the-fly min-degree, row-wise, and at least
   one stronger heuristic (min-degree with a fill-based tie-break, or a
   Markowitz-style criterion).
3. Decide between ordering once up-front as a permutation and continuing to
   choose per minor.
4. Re-state the standard suite's figures under the chosen policy, with the spread.

**Gate.** A documented policy whose worst case over random node numberings is
within a stated factor of its best, and republished suite numbers carrying that
spread. The factor is the deliverable: it is what every later comparison must be
read against.

*Effort: small to moderate — mostly measurement, and the construction already
takes an ordering strategy.*

### H2 — Recursive hierarchy — **DONE**

**Why second.** `HierarchicalDDD` splits the matrix *once*, and on the µA741 that
loses whichever way it is pointed: suppress a small block and ~21 terminals
remain, so the reduced system is nearly the original; suppress a large block (21
internal, 6 terminals) and its cofactor family costs 4842 vertices against 1072
flat. Tan & Shi's 56× comes from never building a large matrix at all — Table IV
shows leaves of 3, 4, 4 and 3 nodes, suppressed bottom-up into middles of 2 and 3.

Tasks:
1. Generalise `HierarchicalDDD` to take a **partition tree** rather than one index
   set, and evaluate bottom-up: each subcircuit is reduced to its cut-set before
   its parent is built.
2. Terminal/cut-set bookkeeping per level — the piece that does not exist today.
3. Reuse what is already right: `DDDFamily` is the correct primitive for a
   block's cofactors, `eval_roots` already prevents the quadratic re-walk, and the
   `_blk_*` symbol indirection already lets a parent treat a reduced block as an
   opaque payload.
4. Transcribe their Fig. 16 partition of the µA741 so the comparison is
   like-for-like rather than a partition of our choosing.

**Gate.** Total vertices across all levels must beat the flat diagram on the
µA741 — under H1's ordering policy, and by more than H1's measured spread, or the
result says nothing. Reaching their ~117 with a comparable three-level two-way
partition is the target; beating flat by less is a partial pass, and losing to
flat is a negative result to write up as P3 was.

*Effort: the largest remaining item.*

### H3 — Noise through the reduction — **DONE**

**Why this follows H2 immediately.** P4 already builds every transimpedance from
one shared cofactor family, but it does so on the **flat** matrix — and H2 showed
that matrix is several times larger than it needs to be on a real amplifier.
Noise is the analysis with most to gain, because it is the one that needs the
*whole* transimpedance vector: the cost reduction attacks is exactly the cost
noise pays once per source.

**The part that is not just plumbing.** An internal unknown cannot simply be
suppressed and forgotten. Every resistor in a suppressed block *is a noise
source*, and its contribution has to arrive at the terminals. So a reduced block
must stamp **two** things into its parent: the admittance stamp H2 already
produces, and an equivalent noise correlation referred to its terminals. Omitting
the second does not fail loudly — it silently under-reports noise, which is the
kind of error that looks entirely plausible in a plot.

Tasks:
1. Carry each block's noise correlation through its suppression: the
   terminal-referred correlation is ``T · CY_block · Tᴴ``, where ``T`` maps the
   block's internal sources to its terminals. ``T`` is built from the *same
   cofactors* already computed for the admittance stamp, so the sharing that
   makes P4 and H2 cheap should serve this too.
2. Combine the referred correlation with the parent's own ``CY``.
3. Wire it into `DDDToolkit.noise_psd`, so `Noise(cir, toolkit=ddd_toolkit)`
   benefits with no change to that analysis — the same zero-lines-outside-DDD
   property P4 had.

**Gate.** The PSD must be **identical** to the flat result, not merely close, on
the ladder and the MFB — a reduction that loses a noise source would still look
reasonable, so equality is the only safe test. On the µA741 the construction cost
should fall by a factor comparable to H2's 6.7×, rather than staying flat.

**Risk to watch.** Correlation is where reductions usually go wrong: sources
inside a block become *correlated* at the terminals even when they were
independent internally, so the referred correlation is a full matrix rather than
a diagonal. Testing against the flat result on a circuit with several noisy
elements in one block is what catches an incorrectly diagonal assumption.

### H4 — Symbolic sensitivity — **DONE**

The plan's own judgement, deferred pending P0 and now due: *pycircuit has none,
and a DDD makes it nearly free — the derivative of a determinant with respect to
a symbol **is** the cofactor.* The infrastructure exists and is tested;
`DDDFamily` already produces exactly those cofactors, shared.

It is also the only remaining item that gives pycircuit a **user-facing
capability it currently lacks**, rather than improving a number.

Tasks: sensitivity of the determinant to a matrix entry from the cofactor; chain
to component parameters; expose on `DDDSolution` and the AC result.

**Gate.** Agreement with finite differences on the ladder, the MFB and the µA741,
at a cost that is a small multiple of a single solve rather than one solve per
parameter.

### H5 — Tier 2 calibration — **DONE**

The Cauer low-pass filter and the cascaded-opamp series from TCAD 2000, both
constructible from `R`/`C`/`VCCS` today. The cascaded series is the one their
SCAPP comparison uses, so it lets the DDD-versus-SoE question be checked against
their 117 vertices vs 539 operations rather than only against our own SoE.

### H6 — Frequency and impedance scaling — **DONE (premise refuted)**

P1 measured denominator coefficients spanning 10⁹⁷ by N = 20, with pole accuracy
degrading to 5e-10. Scaling toward O(1) coefficients is recorded there as the
mitigation and has never been applied; it is the same fix that rescued the GiNaC
backend from a related blow-up.

**Gate.** Pole accuracy at N = 20 improves by orders of magnitude with no
regression at small N.

### H7 — DC and other analyses

Everything built so far is AC and noise. Nothing rules out DC symbolic solving
through the same machinery and nothing has been thought about. Scoping only.

---

## 6. Definition of done for round two

1. ~~**Vertex counts are trustworthy**~~ — done (H1): every figure carries a known
   ordering policy and a 1.42× worst-case spread.
2. ~~**Hierarchy earns its place**~~ — done (H2): µA741 1040 → 156 vertices, 6.7×,
   well beyond that spread.
3. ~~**Reduction reaches the analysis that needs it most**~~ — done (H3): noise
   through the reduction gives an identical PSD, 426× fewer vertices and 180×
   faster on the µA741.
4. ~~**A designer gains something they did not have**~~ — done (H4): sensitivity
   to every device, agreeing with finite differences, at the cost of two solves
   regardless of how many parameters are asked about.

All four round-two objectives are met. What remains (H5–H7) is calibration and
scoping rather than capability.

Per-stage, the round-one rules stand unchanged: code plus tests with the full
suite green, the theory-document section in the **same commit**, `make html`
clean with no `exec-rst` fallbacks, and the gate evaluated with its result
recorded even when negative.

---

## 7. Sequencing notes

- H1 and H2 are done; the spread H1 measured (1.42× on the µA741) is the
  tolerance every later vertex count is read against.
- H3 depends on H2 and on nothing else — it is the direct application of the
  reduction machinery to the analysis that most wants it.
- H1–H6 are done. Only H7 (DC, scoping) remains.
- H6 refuted its own premise, which is worth remembering when reading the rest of
  this plan: the mitigation it proposed had been borrowed from GiNaC's failure
  (exact-arithmetic growth, which scaling does fix) and applied to a different
  one (floating-point root conditioning, which it does not). Two failures that
  look alike in a summary are not necessarily the same failure.
- H5's outcome is worth carrying forward: the diagram-versus-sequence advantage
  is a *small factor*, a few times, in their measurements and ours alike. Claims
  of orders of magnitude belong to the comparison against the **expanded** form,
  not against another shared representation.
- One thing H3 turned up that applies more widely: the `_resolve` fast path and
  unmultiplied factors sped the reduced solve ~200×, and the same pattern —
  sympy substitution in a hot loop — is worth watching wherever payloads are
  sympy objects. Timings recorded before that change understate what the code
  now does.
- H3 is independent of both and could run in parallel; it touches only new code.
- H4 is cheap and worth doing alongside whichever of H1/H2 is active, since it
  supplies reference points for both.
- H5 is self-contained and only affects the s-expanded path.
- Keep paper PDFs out of the repo; they live in `~/pycircuit_agy/papers/ddd/`.
  `benchmarks/paper_extract.py` renders their figures and tables, which are
  300-dpi scans rather than text.
