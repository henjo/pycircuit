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
    └── SymbolicPolyToolkit    fraction-free N(s)/D(s)   <- the num/den contract
        ├── SymengineToolkit   EXPERIMENTAL, measured slower
        ├── GinacToolkit       EXPERIMENTAL, conditional win
        └── DDDToolkit         decision diagrams
```

Read the `SymbolicPolyToolkit` subtree as *"shares the `N(s)/D(s)` contract"*,
not as *"is a kind of polynomial toolkit"* — see problem P3 below.

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

### P2 — `__getattr__` delegation gives misleading errors

A name missing from both class and backend raises an `AttributeError` naming the
*backend module*:

```
numeric.sparse -> AttributeError: module 'pycircuit.circuit._numeric' has no attribute 'sparse'
```

which sends a reader to the wrong file. More seriously in principle: a method
removed from a subclass silently falls through to the backend instead of
failing. **I found no instance of that second case actually happening** — this is
a hazard, not a bug report.

The delegation earns its keep and should stay; the fix is a better
`__getattr__` error message, not removal. Note also that `sparse` is set only on
`SparseNumericToolkit`, so unlike `symbolic`/`poly`/`jax` it is not safe to read
on an arbitrary toolkit. Declare it on the base, or route it through
`supports()`.

### P3 — the toolkit hierarchy claims more than it delivers

`DDDToolkit` extends `SymbolicPolyToolkit` but overrides `supports`,
`ac_solution`, `noise_psd` and `linearsolver_num_den` — nearly everything it
inherits. What it actually reuses is one method, `linearsolver`. The
inheritance asserts a kinship that the code does not use.

A `NumDenMixin` carrying the `N(s)/D(s)` contract would say what is true. Low
urgency, but it misleads every reader who assumes the hierarchy means something.

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

### P5 — six symbolic paths and no user-facing story

`symbolic`, `symbolic_poly`, `soe`, `ginac`, `ddd`, `symengine`. Two are
honestly labelled experimental negative results, which is good practice and
should stay. But the answer to *"which should I use?"* lives in commit messages
and scattered doc pages, not in the API or one table.

`soe` is the odd one out: a module-level `solve_soe()` used by tests, benchmarks
and docs, reachable from no analysis. It is a research prototype promoted to
package code without an integration path — compare DDD, which got
`ac_solution`, a result object and an analysis family.

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

### P8 — two extension conventions, built opposite ways

The transient engine has a second pluggability system: ABCs `Integrator`,
`StepController`, `NonLinearSolver`, `Scaler`, chosen by **string parameter**
through an if/elif factory.

| | toolkit layer | strategy layer |
|---|---|---|
| base | plain class, duck-typed | `ABC` |
| selection | pass an object | pass a string |
| capabilities | flags + `supports()` | class identity |

They are not the same idea twice — one chooses arithmetic, the other an
algorithm inside transient — but there is no stated rule for which a new
extension should use. Concrete wart: `transient.py:188` decides behaviour with
`self.active_integrator.__class__.__name__ == 'EulerIntegrator'`, interrogating
a polymorphic hierarchy by class-name string. An `order` property would do it
properly.

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

### P11 — `_symbolic` is missing primitives `_numeric` has

`_numeric` provides `tanh`; `_symbolic` provides neither `tanh` nor `power`. So
`func.Tanh` — and any element built on it — cannot run symbolically at all,
despite nothing in its design being numeric-specific. The backends have drifted
apart without anything detecting it.

A test that asserts the backend modules expose the same primitive names would
pin this cheaply.

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

### P13 — test coverage is shaped inversely to risk

489 tests. **180 of them (37%) are `test_ddd.py`** — the newest subsystem. The
two most depended-upon modules are the thinnest covered:

| module | fan-in | lines | tests in its main file |
|---|---:|---:|---:|
| `circuit.py` | 29 | 1656 | 20 |
| `elements.py` | 24 | 1947 | 12 |
| `ddd.py` | low | 1972 | 180 |

Newest code best tested, oldest and most-depended-upon least. That is backwards
from where a regression hurts most. (Element behaviour is also exercised
indirectly by analysis tests, so the gap is smaller than the raw counts
suggest — but the direction is real.)

---

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
`Tanh.fprime`).

Open, in the order I would take them:

| # | item | size | blocked on |
|---|---|---|---|
| P11 | `_symbolic` missing `tanh`/`power` | small | nothing |
| P2 | delegation error messages | small | nothing |
| P13 | test coverage inverted against risk | large | nothing |
| P3 | hierarchy claims more than it delivers | medium | nothing |
| P5 | no user-facing backend story | medium | nothing |
| P8 | two extension conventions | medium | a convention decision |

P11 is first among what's left: small and unblocked, and it's the kind of
silent backend drift (§7 P11's own description) worth closing before it grows
another primitive gap.

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
