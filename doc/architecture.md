# pycircuit architecture

*Orientation for anyone — human or agent — about to change this code. It
describes what the design is, why it is that way, and where it does not hold.
Written 2026-07-27 against branch `cna-jax-vectorization`. Claims here were
checked against the code; where something is inference rather than measurement,
it says so.*

---

## 1. What pycircuit is, in one paragraph

A circuit simulator written as a library rather than a program. You build a
netlist out of Python objects, an *analysis* turns it into a matrix equation and
solves it, and a *result* object answers questions about the solution. The
distinguishing feature is that the same analysis code produces a **numeric**
answer, a **symbolic** expression, or a **shared-graph** representation,
depending on one object you pass in.

## 2. The central idea: the toolkit

This is the thing to understand first, because almost every design decision
follows from it.

An analysis never calls `numpy` or `sympy` directly. It calls **toolkit**
methods — `linearsolver`, `dot`, `zeros`, `det`, `exp`, and physical constants
like `kboltzmann`. Swap the toolkit and the same code produces a different kind
of answer:

```python
AC(cir, toolkit=numeric).solve(1e3)        # a complex number
AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)   # N(s)/D(s)
AC(cir, toolkit=ddd_toolkit).solve(s, complexfreq=True)     # a shared diagram
```

`AC` contains no branch on which of these is in play.

### Two axes, and why the distinction matters

Toolkits vary along two independent axes, and conflating them is the most common
way to get confused here:

- **Axis 1 — arithmetic.** What a number *is*: a float, a `sympy` expression, a
  JAX tracer. This decides how elements stamp their matrices.
- **Axis 2 — solving and representation.** What "solve" *returns*: a vector of
  numbers, a rational function, an unexpanded graph.

`numeric` and `symbolic` differ on Axis 1 only. `symbolic_poly`, `ddd_toolkit`
and `ginac_toolkit` share Axis 1 with `symbolic` (all sympy) and differ on
Axis 2. `JAXToolkit` differs on Axis 1 and adds automatic differentiation.

### How a toolkit is built

A `Toolkit` is a thin class wrapping a **backend module** of primitives:

```python
numeric       = NumericToolkit(_numeric)        # numpy
symbolic      = SymbolicToolkit(_symbolic)      # sympy
symbolic_poly = SymbolicPolyToolkit(_symbolic)  # sympy, fraction-free solve
ddd_toolkit   = DDDToolkit(_symbolic)           # sympy, diagram solve
```

`Toolkit.__getattr__` forwards anything the class does not define to that
module. So `toolkit.zeros(...)` reaches `_numeric.zeros` or `_symbolic.zeros`
without the class knowing either exists.

**This delegation is load-bearing, and it is the main extension mechanism.**
The DDD analysis family was added by overriding exactly two methods —
`det` and `cofactor` — because `TwoPortAnalysis` and `FeedbackLoopAnalysis`
already reached their answers through those. No analysis code changed. When
adding a backend, look first for the primitive the existing analyses already
funnel through; overriding it is usually the whole job.

### The class hierarchy

```
Toolkit                       supports(), jacobian(), derivative(),
├── NumericToolkit               ac_solution(), noise_psd()
│   └── SparseNumericToolkit
├── JAXToolkit                 autodiff; the only Axis-1 outlier
└── SymbolicToolkit            sympy LUsolve
    ├── NumDenMixin + SymbolicToolkit -> SymbolicPolyToolkit
    │       fraction-free N(s)/D(s), flat NumDenSolution
    │       ├── SymengineToolkit   EXPERIMENTAL, measured slower
    │       └── GinacToolkit       EXPERIMENTAL, conditional win
    └── NumDenMixin + SymbolicToolkit -> DDDToolkit
            decision diagrams, its own ac_solution/noise_psd
```

`NumDenMixin` carries only the `N(s)/D(s)` *contract*: implement
`linearsolver_num_den` and get `linearsolver` and `supports('num_den')` for
free. `SymbolicPolyToolkit` and `DDDToolkit` are siblings under it, not
parent/child — they answer `ac_solution`/`noise_psd` too differently (a flat
result vs. a shared diagram) to share a base carrying those. `SymengineToolkit`
and `GinacToolkit` *do* genuinely subclass `SymbolicPolyToolkit`: they override
only `linearsolver_num_den` and keep its `ac_solution`/`noise_psd` unchanged
(formerly problem P3, now fixed this way).

### Optional capabilities

Two mechanisms, deliberately:

- **Class flags** — `symbolic`, `poly`, `jax`, all declared on the base so they
  are always safe to read.
- **`supports(capability)`** — string-keyed, for capabilities an analysis opts
  into. Current keys: `'num_den'`, `'ddd'`, `'autodiff'`.

Use `supports()` for anything new. It is what lets an element ask *"will my
derivatives be computed for me?"* without naming a backend.

## 3. Layers

```
  post/            results, waveforms, plotting, Cadence PSF readers
    ^
  analysis         AC, DC, Noise, Transient, TwoPort, PSS, Volterra, ...
    ^                                        (analysis.py, analysis_ss.py, ...)
  toolkit          arithmetic + solve strategy      <- swap point
    ^
  elements         R, C, L, sources, semiconductors  (elements.py)
    ^
  circuit          topology: nodes, branches, hierarchy  (circuit.py)
```

Dependency weight, from the import graph:

| module | fan-in | lines |
|---|---:|---:|
| `circuit.circuit` | 29 | 1656 |
| `circuit.elements` | 24 | 1947 |
| `circuit.toolkit` | 22 | 604 |
| `circuit.analysis_ss` | 13 | 508 |

## 4. The element contract

An element is a `Circuit` subclass declaring `terminals` and `instparams`, and
implementing as many of these as apply:

| method | meaning | used by |
|---|---|---|
| `G(x, epar)` | conductance matrix ∂i/∂v | all |
| `C(x, epar)` | capacitance matrix ∂q/∂v | AC, transient |
| `u(t, epar, analysis)` | independent sources | all |
| `i(x, epar)` | static currents | DC, transient |
| `q(x, epar)` | charges | transient |
| `CY(x, w, epar)` | noise current correlation | Noise |

Nonlinear elements additionally provide **pure** static forms
`eval_i_pure(x, params, epar, toolkit)` and `eval_q_pure(...)`. These exist so a
differentiating toolkit can produce the Jacobian by autodiff rather than the
element hard-coding it:

```python
def G(self, x, epar=defaultepar):
    return self.toolkit.jacobian(self.eval_i_pure, x,
                                 {'r': self.iparv.r}, epar, self._G)
```

`Toolkit.jacobian` returns the precomputed `fallback` matrix; `JAXToolkit`
differentiates. **An element must never ask which toolkit it has.** Elements
that genuinely need two algorithms — device limiting versus a smooth form —
guard on `self.toolkit.supports('autodiff')`.

### State the element once

The base class already derives one side from the other:

```python
def i(self, x, epar=defaultepar, params_tree=None):
    return self.toolkit.dot(self.G(x), x)     # and q from C likewise
```

So **a linear element states only its matrix.** `i` follows exactly, and the
Jacobian *is* the matrix, so there is nothing to differentiate:

```python
def G(self, x, epar=defaultepar):
    return self._G
```

Twelve linear elements briefly carried an `eval_i_pure` as well — a second
encoding of the same conductance, added so JAX could batch them. It bought
nothing (differentiating a linear function to rediscover a stored constant) and
cost a consistency obligation, so it is gone; the JAX paths measured the same
without it.

**A nonlinear element states only its current**, `eval_i_pure`, and gets its
Jacobian from `toolkit.jacobian`. Four elements — `Diode`, `VCVS_limited`,
`VSwitch`, `ISwitch` — still hand-write a stamp *as well*, because their
non-autodiff path applies device limiting. Those four are the only place the
two-sources-of-truth problem can still arise, and they are exactly the elements
`test_element_jacobians.py` exists to police. Do not add a fifth without adding
a case there.

## 5. Worked trace: what `AC.solve` does

1. `dc_steady_state()` finds the operating point (skipped for symbolic
   toolkits, where `x=None`) and stamps `G`, `C`, `CY`, `u` by walking the
   element tree.
2. The reference node's row and column are removed — the matrix is otherwise
   singular.
3. **The Axis-2 hook:**
   ```python
   solution = self.toolkit.ac_solution(ss*C + G, -u, ss, irefnode)
   ```
   A toolkit that can express the answer compactly returns an `ACSolution`;
   the base returns `None`.
4. If a solution came back, the result wraps it and delegates `tf`, `poles`
   and node voltages to it. If not, the analysis falls back to
   `toolkit.linearsolver(s*C + G, -u)` once per frequency.

Step 3 is why adding a representation does not add a branch to `AC`. Put new
solve strategies behind `ac_solution`, not behind an `if` in the analysis.

## 6. How to extend

**A new element** — subclass `Circuit` in `elements.py`, declare `terminals`
and `instparams`, implement the contract in §4. Provide `eval_i_pure` if it is
nonlinear. Do not import a backend.

**A new toolkit** — subclass the closest existing toolkit, override the *fewest*
primitives that change the answer, and declare capabilities via `supports()`.
Guard optional third-party imports so the singleton becomes `None` when the
dependency is absent, and **do not import the backend module at file top level**
(see P1).

**A new analysis** — subclass `Analysis` or `SSAnalysis`. Reach everything
through `self.toolkit`. If it is a variant of an existing analysis under a
different representation, prefer overriding a toolkit primitive over writing a
new analysis: `dddanalysis.py` adds an entire family without any class defining
`solve`.

---

## 7. Known problems

Ranked by what I would fix first. Each was verified.

### P1 — optional dependencies are guarded at the wrong level *(fixed for JAX; the pattern may recur)*

`toolkit.py` guards its JAX import so `jaxtoolkit` becomes `None` when JAX is
missing. That guard was silently defeated for months by an unconditional
`from . import _jaxtoolkit` at the top of the same file, making JAX a hard
dependency of the whole package. Fixed, and pinned by `test_jax_optional.py`.

**The general rule this yields:** a module that imports an optional third-party
package at top level must itself only ever be imported lazily. `_ginac` and
`symengine` follow this correctly today; check any new backend against it.

### P2 — `__getattr__` delegation gives misleading errors *(fixed)*

A name missing from both class and backend used to raise an `AttributeError`
naming the *backend module*:

```
numeric.sparse -> AttributeError: module 'pycircuit.circuit._numeric' has no attribute 'sparse'
```

which sent a reader to the wrong file. More seriously in principle: a method
removed from a subclass would silently fall through to the backend instead of
failing. **No instance of that second case was found actually happening** —
that risk is a hazard worth naming, not a bug report, and the delegation
itself earns its keep and stays (architecture.md §2/§6: it's the main
extension mechanism). `Toolkit.__getattr__` now catches the backend's
`AttributeError` and re-raises naming the toolkit class first, with the
backend module kept as parenthetical context rather than the headline:
`"'NumericToolkit' toolkit (backend 'pycircuit.circuit._numeric') has no
attribute 'sparse'"`. Pinned by `test_toolkit.py`.

The secondary note is fixed too: `sparse` was set only on
`SparseNumericToolkit`, so unlike `symbolic`/`poly`/`jax` it wasn't safe to
read on an arbitrary toolkit (confirmed nothing actually read it that way
today — a latent hazard, not a live bug, same shape as the first half of
this item). Declared `sparse = False` on the base `Toolkit` alongside the
other three flags.

### P3 — the toolkit hierarchy claims more than it delivers *(fixed)*

`DDDToolkit` extended `SymbolicPolyToolkit` but overrode `supports`,
`ac_solution`, `noise_psd` and `linearsolver_num_den` — nearly everything it
inherited. What it actually reused was one method, `linearsolver`, itself a
two-line wrapper that just calls back into `linearsolver_num_den` — i.e.
`DDDToolkit`'s own override. The inheritance asserted a kinship the code
never used.

Checked first whether anything depended on the exact hierarchy before
changing it: no `isinstance(toolkit, SymbolicPolyToolkit)` checks anywhere in
analyses or tests, and the one place `.poly` is read
(`test_symbolic_poly.py`) checks the `symbolic_poly` singleton specifically —
safe to restructure.

Extracted `NumDenMixin` (`toolkit.py`): implement `linearsolver_num_den` and
get the standard `linearsolver` and `supports('num_den')` for free — nothing
about `ac_solution`/`noise_psd`, since that's exactly the part
`SymbolicPolyToolkit` and `DDDToolkit` answer differently. `DDDToolkit` now
inherits `NumDenMixin` + `SymbolicToolkit` directly, a sibling of
`SymbolicPolyToolkit` rather than a subclass of it. `SymengineToolkit` and
`GinacToolkit` are unchanged: they only ever overrode `linearsolver_num_den`
and genuinely want `SymbolicPolyToolkit`'s `ac_solution`/`noise_psd`
defaults, so they keep subclassing it directly — the fix is precise about
which relationships were real and which weren't. Verified behaviourally
identical: `ddd_toolkit.poly`/`.symbolic`/`.supports(...)` and
`linearsolver`'s output are unchanged before and after, pinned by
`test_toolkit.py`; full DDD/GiNaC/symengine/symbolic_poly suites (247 tests)
pass unmodified.

### P4 — ~1528 lines of maintained dead code *(fixed)*

`circuit_cna.py` (1377 lines) and `elements_cna.py` (151) were imported by
nothing — not the package, tests, docs or benchmarks — yet `circuit_cna.py`
duplicated eight core class names including `Circuit` and `SubCircuit`. Git
showed both receiving the Python-3 port, the NumPy 2.0 fixes and the toolkit
refactor, despite neither actually importing under Python 3 at all (both used
Python-2 implicit relative imports — `from constants import *` /
`from circuit import *`).

CNA turned out to mean Compact Nodal Analysis, a real cited technique
(Tlelo-Cuautle & Sarmiento-Reyes, JART 2003 — downloaded and read, see
`doc/cna_conclusions.md`) for reducing matrix size in nullor-containing
circuits, and the user chose to attempt a real revival rather than a plain
delete. The reasoning, staged plan and outcome are in `doc/cna_conclusions.md`
and `doc/cna_plan.md`: Stage 1's hypothesis — that pycircuit's own
NA-compatible elements would let a single zero-`G`, no-branch `Nullor` reach
the paper's "Compacted PNA" size for free — was tested against the paper's
own worked example and **refuted**: the naive element produces a genuinely
singular matrix, not merely a less-compact one, because the paper's row/column
merge-and-delete rules are load-bearing, not an artifact of their own
methodology. Real compaction would need node-graph surgery at circuit-assembly
time — a materially bigger, more invasive feature than an element class,
scoped separately if pursued further. The negative result is pinned by
`pycircuit/circuit/tests/test_nullor_cna.py` (the wrong candidate stays
test-local, not shipped in `elements.py`, since it gives incorrect answers).
`circuit_cna.py`/`elements_cna.py` are deleted regardless of that outcome —
nothing in the finding depended on repairing them, and 1377 of the 1528 lines
were pure duplicate `Circuit`/`SubCircuit` engine with zero CNA-specific
content in the first place.

### P5 — six symbolic paths and no user-facing story *(fixed)*

`symbolic`, `symbolic_poly`, `soe`, `ginac`, `ddd`, `symengine`. Two are
honestly labelled experimental negative results, which is good practice and
should stay. But the answer to *"which should I use?"* lived in commit
messages and scattered doc pages, not in the API or one table.

Added `doc/src/circuit/symbolic_backends.rst` — a comparison table plus a
short decision-tree paragraph, deliberately additive (doesn't touch or
duplicate the existing per-backend pages, just orients a reader among them).
Left specific benchmark numbers out of it on purpose: the existing pages
already regenerate theirs live via `exec-rst` (rule #6 — never paste a
measured number into prose), so this page states each backend's
representation and "best when" qualitatively and links out for the numbers.

`soe` is the odd one out: a module-level `solve_soe()` used by tests,
benchmarks and docs, reachable from no analysis. It is a research prototype
promoted to package code without an integration path — compare DDD, which got
`ac_solution`, a result object and an analysis family. Investigated *why*
before writing it off as a simple gap: `ACSolution` (the interface every
other representation implements) assumes a shared denominator —
`numerators()` over one `denominator()`, `poles()` as its roots — which is
still true of a DDD graph (mathematically `N(s)/D(s)`, just never expanded)
but not of SoE, whose entire point is that inlining into that shape is the
exponential blow-up it exists to avoid. Giving `soe` a real integration is
therefore a genuine architecture question — extend `ACSolution` to a
denominator-less representation, or give SoE its own analysis family the way
DDD has one — not a wiring task like P3's mixin was. Recorded as the
explanation on the new page and left as an explicit open question rather than
attempted here; the "no user-facing story" complaint this item was about is
resolved regardless of how that question is eventually answered.

### P6 — batched evaluation now lives in the toolkit *(fixed)*

The JAX vectorisation had put its machinery in core `circuit.py`: building the
per-class evaluation groups, and two copies of a vmapped stamping loop. That
made the topology module carry autodiff-shaped code, and hid a bug (below).

It is now three toolkit hooks, with base implementations that do nothing:

| hook | base | `JAXToolkit` |
|---|---|---|
| `evaluation_groups(circuit)` | `{}` | groups elements by class |
| `batched_contributions(...)` | `None` | one vmapped call per class |
| `generate_batched_eval(cls, method)` | — | jit/vmap/jacfwd evaluator |

`circuit.py` offers the opportunity and accumulates the result; it names
nothing JAX-specific. `-96` lines there, and the duplicate
`generate_batched_eval` in `_jaxtoolkit.py` is gone — there is one
implementation now.

**The bug it was hiding:** each stamp asks every group for both a value and a
Jacobian, but a diode has no charge and a capacitor no static current. Any JAX
circuit mixing a nonlinear device with a reactive one raised `AttributeError`
on `G`, `C`, `i` *and* `q` — i.e. every stamp. The missing pure form now
contributes zeros. Pinned by a test that compares the batched stamp against the
per-element one rather than merely checking it returns something.

### P7 — the toolkit reaches into element instances *(fixed)*

The odd coupling, resolved.

`Circuit.__init__` used to end with:

```python
if self.toolkit.supports('autodiff'):
    if hasattr(self, 'eval_i'):
        self.toolkit.generate_eval_i_and_G(self)
    if hasattr(self, 'eval_q'):
        self.toolkit.generate_eval_q_and_C(self)
```

and `JAXToolkit.generate_eval_i_and_G` **installed a method on the element
instance** (`element.eval_i_and_G = ...`, plus a per-instance `_jax_cache_i`),
which `Circuit.G`, `C`, `i` and `q` then discovered with `hasattr`. So the
toolkit mutated element state, and the base class picked the result up by name.

**What was measured (2026-07-27, extended 2026-07-28).** The 2026-07-27 pass
established: no shipped element defines `eval_i` or `eval_q`, so the guard was
false for every element in the package; the only definition anywhere was
`JAXDiode` inside `test_jax_autodiff.py`, a test-local subclass, and that test
called `generate_eval_i_and_G` **by hand** because `JAXDiode.__init__` assigns
its toolkit *after* `super().__init__()`, so the guard never sees an autodiff
toolkit at construction time. The follow-up measurement (2026-07-28) confirmed
this with a `raise` planted in both `generate_*` methods and a full suite run:
**498 passed, 6 skipped, 1 failed**, and the one failure was exactly
`test_jax_autodiff.py`'s manual call. The `generate_eval_q_and_C` probe never
fired **at all** — not automatically, not even by hand — so the charge/
capacitance half of the mechanism had zero exercise anywhere in the codebase,
a stronger form of "dead" than the `eval_i` half.

This was also a *second* pure-form convention: bound `eval_i(self, x, epar)`
here, against static `eval_i_pure(x, params, epar, toolkit)` used by
:meth:`~pycircuit.circuit.toolkit.Toolkit.jacobian` and by batching. Two
conventions for one job was the duplication P9 is about, one level up.

**Decision:** delete, made by the user after the measurement above (deferred
here on 2026-07-27 as a product judgement the code could not make for itself —
vestigial design superseded by `eval_i_pure`/`jacobian`, or a deliberate
escape hatch for an out-of-tree element). The evidence — no in-tree user, a
differing signature convention, and the `eval_q` half unexercised even by
hand — favoured vestigial, and that is the call that was made.

**What changed.** Removed `generate_eval_i_and_G`/`generate_eval_q_and_C` from
`toolkit.py` (~45 lines), the `autodiff`/`hasattr` hook at the end of
`Circuit.__init__`, and the four `hasattr(self, 'eval_i_and_G'/'eval_q_and_C')`
branches in `Circuit.G/C/i/q` — those methods now just do what the base always
did: read the stamped matrix, or derive `i`/`q` from it. `test_jax_autodiff.py`
was rewritten first, per the plan above, to point at the mechanism that was
always the real one: it builds a real `Diode` with `toolkit=jaxtoolkit` and
checks that `Diode.G`/`Diode.i` — which already guard on
`supports('autodiff')` and call `eval_i_pure` / `toolkit.jacobian` — agree with
the numeric-toolkit stamp. Coverage of "JAX autodiff produces a correct
Jacobian for a real element" did not dip; it improved, since the old test only
ever exercised a mechanism nothing else used. Full suite: 499 passed, 6
skipped, 0 failed — same test count as the probe run, confirming nothing else
depended on the removed path.

### P8 — two extension conventions, built opposite ways *(fixed)*

The transient engine had a second pluggability system: ABCs `Integrator`,
`StepController`, `NonLinearSolver`, `Scaler`, three of the four chosen by
**string parameter** through an if/elif factory (`Transient._get_integrator`,
`Analysis._get_nrsolver`, `Analysis._get_scaler`).

| | toolkit layer | strategy layer (before) |
|---|---|---|
| base | plain class, duck-typed | `ABC` |
| selection | pass an object | pass a string |
| capabilities | flags + `supports()` | class identity |

Decided in favour of the toolkit convention (pass an instance), not kept as a
choice between equals: `StepController` — the fourth ABC in the same file —
**already worked this way** (`transient.py` constructs `IntegralController()`
directly, only overriding it if the caller injected one), so this wasn't
picking a side, it was making the other three consistent with a pattern
already established next to them. Concrete reasons beyond consistency: an
`Integrator`/`Scaler` often needs its own sub-configuration
(`EulerIntegrator(lte_formula=...)`, and `Scaler`'s hardcoded
`SinkhornKnoppScaler(max_iter=5)` had no way to be configured differently at
all through the public API) — the string convention has to thread that as a
*second*, separately-named parameter and pair the two back together inside
the factory; passing the object directly carries its own configuration, no
pairing required. `Integrator.check_order_drop()` already returns a (possibly
different) `Integrator` instance internally, so the object-as-unit-of-
exchange story was already how the ABC worked beneath the string entry point.

The string form is **not** kept as compatibility sugar: this project has
never had a version bump past `0.0` (`setup.py`) and has no `pyproject.toml`
or external distribution, so there is no downstream API contract to protect
— keeping a permanent dual-path "accept a string or an instance" convention
exists to protect *other people's* code against a large installed base,
which doesn't apply here. Passing the old string form now raises a clear
`TypeError` immediately (e.g. `"integrator must be an Integrator instance...,
not 'euler'"`), rather than failing confusingly deep inside a Newton loop.

Fixed the concrete wart as part of the same change: `transient.py:188`
compared `self.active_integrator.__class__.__name__` to the literal string
`'EulerIntegrator'`, and referenced `self.par.method`, which no longer
exists. Since the object itself is now what's threaded through,
`_effective_method` is derived directly from the live instance
(`type(self.active_integrator).__name__`) — more accurate than the old
binary euler/not-euler check, and it needed no new API surface to get there.

`Transient.parameters`'s `lte_formula` was removed as a standalone
parameter: it only ever existed to be paired back up with `method` inside the
factory, and is now fully expressible on the `Integrator` instance itself
(`Gear2Integrator(lte_formula='ywr')`). `Analysis.parameters`'s `nrsolver`/
`scaler` keep their names but their default changes from a string to `None`
(resolved to `StandardNewton()`/`NoneScaler()` in the getters, matching
`StepController`'s existing default-construction pattern exactly).

Updated every internal call site using the old string/attribute form (none
were external — none exist to update): `test_breakpoints.py`,
`test_minstep.py`, `test_vss_gear2.py`, `test_bypass.py`,
`test_analysis_transient.py` (three separate patterns: constructor kwarg,
post-construction `tran.par.method = ...` assignment, and the
`lte_formula`-paired-with-`method` helper in `test_lte_formula_ywr`), and
`doc/src/circuit/lte_dae.rst`'s prose and its live `exec-rst` benchmark
table. New `test_strategy_objects.py` pins the convention itself: defaults,
accepting an instance, and the `TypeError` on the old string form. Full
suite: 516 passed, 6 skipped, 0 failed. Doc build succeeded, 2 warnings
(unchanged baseline); the `lte_dae.rst` table verified still rendering live,
not falling back to source.

### P9 — the invariant nothing enforced: `G` vs ∂i/∂x *(fixed, and now tested)*

The promise in §2 is that one element definition serves every toolkit. The place
it failed is worth keeping on record, because the failure was **silent**.

An element may state its conductance matrix twice: as a hand-written `G()`
stamp, and implicitly through `i()`, which a differentiating toolkit turns into
the same matrix by autodiff. Nothing made those agree. When they disagree Newton
gets a wrong Jacobian and often still converges -- more slowly, or to a slightly
different point -- so nothing looks broken.

`VCVS_limited` had three defects at once, all pre-dating the JAX work:

1. `func.Tanh.fprime` returned a hard `0` (P10), so the stamp had no
   input-to-output coupling at all.
2. `i()` computed `v_out - f'(u)*f(u)` and never used the gain `g`, while `G()`
   stamped `g*f'(u)`. Beyond the knee it even carried the opposite sign.
3. `eval_i_pure` inlined `offset + level*tanh(...)` as "func.Tanh", which is not
   what `func.Tanh.f` computes -- so the differentiated and stamped paths were
   limiting differently.

The `G()` stamp is correct and was taken as the specification. The residual with
its derivatives is `v_outn - v_outp - g*f(v_inn - v_inp)`, which reproduces the
stamp to finite-difference precision at every operating point tested; `i()` and
`eval_i_pure` were corrected to that, and now use one definition of `f`.

**`test_element_jacobians.py` pins the invariant** -- it finite-differences
`i()` and compares against `G()` across operating points that straddle the knee.
Verified to fail on all three historical defects when they are reinstated.
Extend `CASES` there when adding a nonlinear element; it is the cheapest guard
this codebase has against a whole class of silent error.

The limiter itself was wrong too. The specified behaviour is

```
vout = level * tanh(g * vin / level)
```

with the gain **inside** the limiter, which is what makes `level`
output-referred: the output clamps at ±`level` whatever the gain, while the
small-signal gain is still exactly `g`. `offset` is input-referred.

`func.Tanh.f` was a bare `tanh((x-offset)/level)`, missing the `level` factor,
so it saturated at ±1 regardless of `level` and had slope `1/level` at the
origin rather than 1 — every `VCVS_limited` gain was off by that factor. It is
now a unit-slope limiter, `level * tanh((x-offset)/level)`, and `VCVS_limited`
passes it the *amplified* input.

`fprime` and `F` were updated with `f`, because the three must stay mutually
consistent: `VCVS_limited` stamps its Jacobian from `fprime` while computing its
residual from `f`, so a change to one alone silently reintroduces exactly the
class of bug above. The stamp also evaluates `f'` at the amplified,
offset-shifted input, since that is what the limiter sees.

Verified across `g` = 1, 3 and 29: gain exactly `g`, saturation exactly
±`level` and independent of `g`, and `G` equal to `∂i/∂x` at every operating
point including past the knee.

### P10 — a 17-year-old shadowed derivative *(fixed)*

`func.Tanh.fprime` was:

```python
def fprime(self, x):
    return 0
    return (1-self.toolkit.power(self.toolkit.tanh(...), 2))/self.level
```

Introduced 2009 (`ad97d98`). Every `VCVS_limited` Jacobian therefore had zero
input-to-output coupling, and Newton converged without it.

The shadowing was **not careless**: the expression below calls `toolkit.power`,
which no backend has ever defined, so simply deleting the `return 0` raises
`AttributeError`. It was a workaround for a missing primitive that was never
finished. The fix uses `**` instead, which numpy, sympy and JAX all implement —
*reach for an operator before adding a primitive to every backend.* All ten
transient stress tests pass with the corrected derivative.

This one also shows why P9's finite-difference test matters: no existing test
could distinguish a correct Jacobian from a zero one, because Newton still
converged.

**Unreachable code is worth grepping for generally** — an AST scan found four
`return`-followed-by-code sites. `elements.py:1840` is another live one: a
transmission-line `u()` returns before its entire transient branch, so a TLine
contributes no transient stimulus. Not investigated further.

### P11 — `_symbolic` was missing `tanh` *(fixed)*

Corrected from the original framing: `power` was never actually a live gap —
nothing in the codebase calls `toolkit.power` (the one call site was inside
the shadowed, unreachable half of `Tanh.fprime`, removed by P10's fix in
favour of `**`), and `_numeric` itself doesn't define `power` either. The only
real asymmetry was `tanh`: `_numeric.py` imports it from `numpy`, `_symbolic.py`
never imported it from `sympy`. One-line fix.

The reason this mattered at all is worth stating precisely, since it's easy
to overstate: `tanh` (hence `func.Tanh`, hence `VCVS_limited`) is nonlinear,
and every toolkit in the `symbolic_poly`/`ddd`/`ginac`/`soe`/`symengine`
family is fundamentally a *linear*-circuit method — `AC.dc_steady_state`
explicitly skips finding an operating point for symbolic toolkits (`x=None`),
so a nonlinear element's `G`/`i` are never meaningfully evaluated through any
of those paths regardless of what primitives `_symbolic` has. The one real,
reachable nonlinear-symbolic path is `SymbolicDC` (`symbolicdc.py`) — it
builds and solves the full symbolic KCL equation system with `sympy.solve`
rather than numeric Newton iteration, and already had a passing test
(`test_nonlinear`) using `Diode`. Fixing `tanh` without checking whether
`SymbolicDC` could actually *use* it would have been declaring victory on an
import, not a capability — the standard this document tries to hold itself
to (§11 of `CLAUDE.md`: verify the output, not the exit code).

It couldn't, not fully: `SymbolicDC.solve()` assumed `sympy.solve` always
returns a dict (`{symbol: value}`). True for linear systems and for the
single-equation nonlinear case (which the code special-cases by hand), but
`sympy.solve` returns a **list of solution tuples** for a genuinely nonlinear
*multi*-equation system — unreachable before `tanh` existed, since nothing
could build one, so this was a latent bug rather than a previously-passing
path that broke. Fixed by handling both shapes. Verified end-to-end (not
just that the import resolves): a `VCVS_limited` circuit now solves
symbolically to `level*tanh(g*Vin/level)`, exactly the corrected limiter
semantics P9 established, pinned by `test_nonlinear_multivariable` in
`test_analysis_symbolicdc.py`.

### P12 — smaller, verified

- `analysis.py:155` — `fsolve(..., toolkit='Numeric', ...)`: a **string** where a
  toolkit object is required, so the default is unusable (`str.linearsolver`).
  Every real caller passes `toolkit=`, so it is a trap rather than a live bug.
- `circuit.default_toolkit` is a process-wide mutable global. Already
  deprecated in favour of `use_toolkit()`, but tests still assign it directly
  and must save/restore.
- `setup_analysis` is an optional hook discovered by `hasattr` and provided by a
  *backend module* rather than a toolkit class — a third capability idiom.
- Comments left mid-thought in production code, e.g. `elements.py` VCVS_limited:
  *"but wait, the plan is to vectorize! ... Let's see if tests pass."*

### P15 — slow transient tests dominate the suite *(route 3 taken; the underlying cost is still open)*

**Route 3 below has been applied**: the transient-heavy tests are now marked
`@pytest.mark.slow` and **deselected by default** by `pytest.ini`. Run them
with `pytest pycircuit -m slow`, or the whole suite with `-m ""`. **CI and any
release check must use `-m ""`.**

**The framing below was measured again on 2026-07-29 and had gone stale.** It
said one test was 70% of the runtime. It is not, and no longer was even before
the marking:

| test | time |
|---|---|
| `test_stress_stiff_rlc_pulse` | 98 s |
| `test_determinant_sensitivity_is_exact[5]` | 46 s |
| `test_prediction_visibly_fails_outside_the_validity_range` | 38 s |
| `test_higher_truncation_improves_accuracy` (2 params) | 61 s |
| `test_operating_point_and_fundamental_agree` | 25 s |
| `test_stress_charge_pump` / `full_wave_bridge` / `lc_tank` | 53 s |

Full suite 502 s over 673 tests. The stress test is **20%**, not 70% — the
suite has roughly doubled in size since the original note, and the cost is now
spread across many transient-based tests rather than concentrated in one. Any
plan built on "fix the one test" would have been aimed at the wrong target.

**What is still open, and it is the substantive part.** Marking tests slow
hides the cost; it does not remove it. The circuits above are genuinely
expensive and several take far longer than the physics requires — the `tanh`
test added the same day integrated 100 drive cycles where the node settles in
under one, a 13x overcharge that a review caught only by accident. **A pass
over the slow set asking "how many cycles does this actually need?" is likely
to recover most of the time without touching what any test checks.**

**And a real risk was accepted by marking them.** These are the *only* tests
that compare the analysis against a time-domain reference. That is exactly the
evidence class that caught the `g_2` sign error which four published-formula
gates could not see (`doc/distortion_mimo_plan.md` section 6). A default that
skips them is a default that skips the strongest checks in the distortion work.
Hence the insistence on `-m ""` in CI.

*Original note follows, kept because its analysis of the stress test is still
accurate.*

`test_analysis_transient_stress.py::test_stress_stiff_rlc_pulse` was measured
at **~266 s** against a full-suite total of ~380 s, at a point when every other
test together accounted for roughly 25 s.

It is not a regression. The timing was measured at HEAD and again at a
`git worktree` checkout of a commit several fixes earlier, identical both
places, when a change was suspected of having caused it. It is a property of
the test as written, and it is why the suite needs `--timeout=300` rather than
the default 120 — see the note in the running-tests recipe.

Where the cost comes from, for anyone picking this up:

- `tend=25e-6` with `timestep=1e-9` asks for 25 000 steps *before* adaptivity
  has any say.
- `_compare_methods` then runs the whole thing **twice** — once with
  `coupled_lte=False` and once with `True` — since comparing the two step
  controllers is the point of the file.
- The circuit is genuinely stiff (`R=1`, `L=1 mH`, `C=1 nF`, with 1 µs pulse
  edges), so the adaptive controller has real work to do at every edge.

**It would be good to reduce this**, and there are three plausible routes,
none of them yet tried:

1. **Start with a larger `timestep`.** It is the *initial* step; the adaptive
   controller is supposed to find its own. If 1 ns is merely a conservative
   opening guess rather than a requirement, raising it may cost nothing.
   Check first whether the step controller actually converges to the same
   trajectory from a coarser start — that is a measurement, not an assumption.
2. **Shorten `tend`.** The assertion depends on landing at 25 µs inside the
   second pulse (20–30 µs is high), so this needs the pulse period changed
   too, not just the end time.
3. **Mark it slow and exclude it by default** (`@pytest.mark.slow`), keeping
   it in full runs. This takes the routine suite from ~6.5 min to about 2 min
   and costs nothing analytically — but it is a change to how the suite
   behaves for everyone, so it is a decision rather than a cleanup.
   **This is the route taken, 2026-07-29, at the maintainer's request.**
   Routes 1 and 2 remain open and are the ones that would actually reduce the
   work rather than defer it.

Route 1 is the cheapest to evaluate and the only one that might be free.
Whatever is done, the guard to keep is that the test's *purpose* — comparing
the adaptive and coupled step controllers on a stiff circuit — survives; a
faster test that no longer stresses the controller would be worse than a slow
one that does.

### P13 — test coverage is shaped inversely to risk *(the doctest gap fixed; the rest stays open)*

Re-measured 2026-07-28 rather than trusted: 511 tests, 180 of them (35%,
was 37% at 489) still `test_ddd.py` — the newest subsystem. The direction the
original count pointed at is real, but "tests in its main file" was a proxy;
`pytest-cov` (worked around a hard crash where its import-tracing hook
collides with `jaxlib`'s native module init — block `jax` at the
`sys.modules` level first, the same trick `test_jax_optional.py` already
relies on) gives an actual number: **`circuit.py` 83% covered (121/699
statements missed), `elements.py` 86% (109/802 missed)** — thinner than
`ddd.py`'s tests-per-line ratio, on the two modules with the highest fan-in
in the codebase (29 and 24 respectively).

**What the gap actually was, found by reading the missing lines rather than
just counting them:** both files carry many doctests documenting their public
API (`get_node`, `add_terminals`, `save_current`, `name_state_vector`,
`Quantity`, most elements' constructors...) but none of them ran under
pytest — gated behind `if __name__ == "__main__": doctest.testmod()`, which
only fires when a module is executed as a script, never on import. Enabling
`--doctest-modules` on just these two files surfaced **20 failures**, and
tracing each one down (not just batch-updating expected output) found:

* **Two real, previously-undetected bugs**, both invisible for the same
  reason. `Quantity.__repr__` (circuit.py) unconditionally raised
  `ValueError("APA")` before its actual logic — introduced 2010, unreachable
  ever since, the same shadowed-dead-code shape as `func.Tanh.fprime`'s
  2009–2026 bug (P10). `Circuit.name_state_vector` had two bugs at once:
  `x[:len(self.nodes)][0]` double-indexed into a scalar instead of slicing
  (only "worked" because its own doctest happened to pass a `(1,1)`-shaped
  array rather than the flat 1-D `x` every other caller in the codebase
  uses), and its branch-naming loop unpacked `enumerate(zip(...))` — a
  2-tuple — into three names. Nothing in the live codebase calls this method,
  so nothing else was affected, but it's public API and the docstring
  advertises a working example.
* **Six stale expected-output strings**, purely numpy-version drift
  (`array([ 0.5,  0. ])` → `array([0.5, 0. ])`; bare floats → `np.float64(...)`
  repr; `array([[ 0.0000e+00, ...`  → `array([[ 0.e+00, ...`) and one
  intentional-but-undocumented API evolution (`Transformer.branches` is a
  list at the instance level, not the tuple the class declares — needed so
  `append_branches`'s `.extend()` works, confirmed before treating it as
  a mere doctest update rather than a behaviour change).
* Nine doctests failing purely on **Python-2-era absolute imports**
  (`from elements import *`, `from dcanalysis import DC`) never updated for
  the package-relative form — the exact same staleness pattern P4 found in
  `circuit_cna.py`/`elements_cna.py`, except here on live, maintained code.

All fixed; all 29 doctests across both files now pass, and
`test_doctests.py` runs them as a permanent part of the ordinary suite (no
special flag needed) so this can't go silently dead again.

**Scope boundary, deliberate.** Checked what `--doctest-modules` finds
project-wide before deciding how far to take this: **15 collection errors**
in largely legacy, peripheral modules — `pycircuit/post/cds`,
`pycircuit/post/jwdb`, `pycircuit/sim/gnucap` (`import sys, StringIO`, a
Python-2-only stdlib module, still unported), a couple of broken example
scripts, and the JAX-dependent modules (an artifact of the jax-blocking
workaround used to measure coverage, not a real gap). That's a separate,
larger cleanup — reviving several already-known-neglected corners of the
tree — not something to fold into closing out this item. *Reconsider if*
those modules are touched for another reason and it becomes cheap to fix
their imports at the same time.

**What's still genuinely open.** The doctest-collection mechanism is fixed
for the two flagged modules, and the two bugs it had been hiding are gone,
but this doesn't make `circuit.py`/`elements.py` well-tested in the
"P13 is large" sense the doc originally meant — 83%/86% line coverage still
leaves real gaps (much of the *rest* of the missing 121+109 lines is
plausible defensive/error-branch code, not audited line-by-line here), and
the test-count *shape* (DDD still 35% of the suite) is unchanged. This item
stays open for anyone who wants to extend it; what's fixed is recorded so it
doesn't get redone by accident.

### P14 — the legacy corners `--doctest-modules` found, deferred out of P13 *(fixed)*

P13's scope check (does turning on doctests project-wide open a bigger can
of worms than the two flagged modules) surfaced **13 real collection
errors** (verified with real `jax` importable, not the coverage
workaround's blocked-jax artifact, which added two more that don't
reproduce otherwise). Resolved bucket by bucket, each on its own merits:

**Expected-missing optional native/vendor dependencies** — correctly failing
without proprietary tools this environment doesn't have, the same shape as
the already-accepted `pycircuit.post.cds` baseline in §8's conventions:
`post/cds/psflibpsf.py` (`libpsf`), `post/cds/psftool.py` (`psf`),
`post/jwdb/jwdbresult.py` + `post/jwdb/test/example.py` (`jwdb`, a
proprietary waveform-database format), `post/jwdb/configure.py`
(`sipconfig`, the PyQt/SIP build tool for `jwdb`'s native extension),
`circuit/xdot.py` (`gobject`, a GTK graph-viewer dependency). Not bugs; left
alone. `post/cds/psftoasc.py`'s `FileNotFoundError` turned out to be a
different, benign shape: a `#!`-shebang CLI script with top-level
side-effecting code (`open(sys.argv[1])`), which only misbehaves because
`--doctest-modules` imports every file as a module — nothing imports it as
a library, so nothing to fix.

**`sim/gnucap/session.py` + `simulation.py` — deleted.** Confirmed fully
orphaned: `sim/gnucap/__init__.py` imports the *generic*
`pycircuit.sim.simulation`, not these; nothing else references them.
Independently broken on top of the `StringIO` import that started this
(an undefined `EngineError`, Python-2-only `types.StringType`) and
superseded by the direct-gnucap-bindings path `test_gnucap.py`/
`test_gnucap_direct.py` actually use (skipped here for lack of a real
gnucap install, not broken). Matches the P4 precedent.

**`sim/tests/simulation_class_example.py` — deleted.** Never valid code: a
design sketch referencing `Circuit`/`Simulation`/`Alter`/`SetVariables`/
`SimGroup`, none defined or imported anywhere.

**`circuit/volterra.py` — fixed, and left honestly incomplete.** The
top-level import (`AC` moved to `analysis_ss.py`) was the smallest of what
running its doctests (not just reading them) turned up: `self.c` where the
base class sets `self.cir`, an undefined `symbolic_linsolve` and dead
`toMatrix` (deleted, never called), `InternalResult` where only
`InternalResultDict` exists (which had a matching bug of its own —
`__len__` read `self.results`, `__init__` sets `self.items`, fixed
separately since the class is used well beyond this file), a stray `self`
parameter on a module-level function, a sympy-version change in how `diff`
applies to a plain array (now elementwise), and `NLVCCS` (a test-only
nonlinear element defined in the same file) hitting the same G-vs-i
mismatch P9 fixed for the shipped elements. `Volterra.solve()`/`.run()`
remain a deliberately unfinished stub — the call to `K()` (the complete,
tested Taylor-series primitive) was commented out in 2008 and never
finished, and investigating why turned up a second, structural reason:
`AC(toolkit=symbolic)` never computes an operating point to linearise a
nonlinear element around (the same gap P11 documented for `SymbolicDC` vs
`AC`), so a circuit with a real nonlinear element fails inside the `AC`
solve itself before `solve()` even lists them. Completing the orchestration
means implementing the Wambacq & Sansen algorithm from a source not in
hand — not attempted; the docstring says so plainly instead of claiming
output the class never produced.

**`circuit/examples/multi_feedback_filter.py` and `regulatedcascode.py` —
fixed and verified, not just made to run.** Both were broken by real API
drift: the 2008-era dict-style `res['out']` never worked against `AC`'s
result object (empty by design), and a symbolic-toolkit `z0` crashed
`NPortS.Y` in three independent ways (`np.sqrt`/`np.real` are ufuncs that
can't dispatch on a sympy object; two of `.Y`'s three `np.linalg.inv` calls
were needless — `Gref`/`Zref` are scalar multiples of the identity, whose
inverse is just the reciprocal; the third, a genuine matrix inverse,
needs `sympy.Matrix.inv()` for object-dtype input). Fixed by adding a
`real` primitive to `_symbolic.py` (sympy has `sqrt` already but not `real`
under that name — it's `re` — the same backend-parity gap P11 found for
`tanh`) and threading a `toolkit` through `NPortS` so `.Y` can call
`self.toolkit.sqrt(self.toolkit.real(...))` instead of hardcoding numpy.
`multi_feedback_filter.py` also needed its own migration to `res.v(node,
gnd)`. Along the way, `symbolicapprox.py` (imported by both examples for
an `approx()` call) turned out to have two independent bugs in one line —
`series(..., point=0, ...)` (sympy renamed the keyword to `x0`), and
`.subs({'t': 1}).removeO()` in the wrong order, which silently zeroed the
result for 18 years (substituting a concrete value into an expression that
still carries a symbolic `O(t**3)` term collapses it before `removeO()`
ever runs) — its own docstring's claimed example output was itself never
verified and is now the checked value. Every fix was checked against real
output, not just a clean exit: the MFB filter's DC gain (`-R1/R3`) and
input impedance share the same characteristic polynomial as its transfer
function, a real cross-check that this is correct and not merely
non-crashing.

New tests: `pycircuit/circuit/tests/test_examples.py` runs both scripts as
subprocesses and checks real output content; `test_doctests.py` (P13)
extended to `volterra.py` and `symbolicapprox.py`;
`pycircuit/post/tests/test_internalresult.py` pins the `InternalResultDict`
fix. Full suite: 521 passed, 6 skipped, 0 failed. Doc build: succeeded, 2
warnings (unchanged baseline).

Aside, found while investigating a suspected regression from the P8 work
(the user correctly asked whether it was implicated, given
`test_stress_stiff_rlc_pulse` is a `Transient` test): measured directly,
via a git worktree at the pre-P8 commit, that this test takes ~266s both
before and after — a pre-existing property of the test, close to or over
the standard 120s per-test timeout independent of any code here, not a
regression. Worth `--timeout=300` (or running it in isolation) rather than
the default when it shows up as a timeout.

## 7.1 State of the list

Fixed during this review, each with a regression test: P1 (JAX made optional
again), P4 (`circuit_cna.py`/`elements_cna.py` deleted; the CNA revival
attempt refuted at Stage 1 and the negative result kept as a regression test
rather than the dead files themselves — see `doc/cna_conclusions.md`,
`doc/cna_plan.md`), P6 (batched evaluation moved into the toolkit, fixing an
`AttributeError` on every stamp for circuits mixing nonlinear and reactive
elements), P7 (the per-instance `eval_i_and_G`/`eval_q_and_C` mechanism
deleted — measured dead, decided vestigial, test rewritten to the real
`eval_i_pure`/`toolkit.jacobian` path), P9 (`G` versus `∂i/∂x`, plus the
`VCVS_limited` residual and limiter semantics), P10 (the shadowed
`Tanh.fprime`), P11 (`_symbolic` missing `tanh` — `power` turned out not to be
a real gap on either side; fixing `tanh` also exposed and fixed a latent
`SymbolicDC.solve()` bug that could only be reached once a multi-equation
nonlinear symbolic system was actually possible to build), P2 (`__getattr__`
now names the toolkit class, not just the backend module; `sparse` moved to
the base `Toolkit` alongside `symbolic`/`poly`/`jax` so it's safe to read on
any toolkit), P13's doctest-collection gap (both real bugs it was hiding —
`Quantity.__repr__`'s shadowed `raise`, `name_state_vector`'s double bug —
found and fixed, all 29 `circuit.py`/`elements.py` doctests now run as part
of the ordinary suite via `test_doctests.py`), P3 (`NumDenMixin` extracted;
`DDDToolkit` is now a sibling of `SymbolicPolyToolkit` under it rather than a
subclass overriding nearly everything inherited; `SymengineToolkit`/
`GinacToolkit` unchanged, since they genuinely wanted what they inherited),
P5 (a new `doc/src/circuit/symbolic_backends.rst` orients a reader among the
six symbolic toolkits; `soe`'s missing integration is explained rather than
fixed — a real architecture question, not a wiring gap), P8 (`Integrator`/
`NonLinearSolver`/`Scaler` now take an instance, matching `StepController`'s
existing convention; the old string form removed outright rather than kept
as sugar, since there is no external API contract to protect), P14 (all
five buckets resolved on their own merits — two deletions, one module fixed
and left honestly incomplete, two example scripts fixed and verified
against real output, plus three supporting library bugs the investigation
turned up along the way: `InternalResultDict.__len__`, `NPortS.Y`'s
symbolic-toolkit incompatibility, and two independent bugs in
`symbolicapprox.approx`).

P13 stays partly open by nature — "large," not a bounded bug — so it isn't
struck off the list, just narrowed: the doctest mechanism and the two bugs
behind it are fixed and won't get redone by accident; the broader
line-coverage gap on `circuit.py`/`elements.py` (83%/86%, measured) and the
test-count shape (DDD still 35% of the suite) remain exactly as open as
before, for whoever wants to extend it next.

**Nothing is open.** Every other item in this list (P1-P14) is now fixed or
explicitly resolved — P13's residual scope is the one exception, by its own
nature (see above), and P5's `soe` integration question is likewise
deliberately left as an open architecture question rather than a queued
task. If this document is being read to decide what to work on next, there
isn't a next item here; it means either picking up one of those two open
questions, or something not yet surveyed at all.

**P15 (slow transient tests dominate the suite) is partly addressed**, added
2026-07-29 while working on the distortion analysis. The transient-heavy tests
are now marked `slow` and deselected by default; **CI must run `-m ""`**,
because those are the only tests comparing the analysis against a time-domain
reference. The underlying cost is untouched and the original "one test is 70%"
framing turned out to be stale when re-measured — see P15 for the current
numbers and for why route 1 (fewer integration cycles) is still the one worth
doing.

## 8. Conventions

**Running things** — the suite needs `PYTHONPATH` set to the repo root, because a
stale root-owned egg on the system path shadows the source:

```bash
PYTHONPATH=$PWD MPLBACKEND=Agg python -m pytest pycircuit -q
```

Docs are built with `python -m sphinx` (the Makefile's bare `sphinx-build` is
not on PATH):

```bash
cd doc && MPLBACKEND=Agg python -m sphinx -b html -d build/doctrees src build/html
```

Two long-standing autodoc warnings about `pycircuit.post.cds` are the expected
baseline; anything beyond them is new.

**Live documentation.** `.. exec-rst::` blocks run Python at build time and
insert the printed reST, so tables and measurements regenerate on every build
rather than going stale. A block that fails renders its source with a warning
and a visible admonition — if you see that admonition in the output, the page is
not proving what it claims. Never paste a measured number into prose.

**Where to look for design reasoning.** `doc/ddd_plan.md` and
`doc/ddd_conclusions.md` record the reasoning and the negative results behind
the DDD subsystem; `doc/ginac_fully_symbolic.md` does the same for the compiled
backends. Negative results are kept deliberately — if something here looks
unfinished, check whether it was measured and rejected before redoing it.
