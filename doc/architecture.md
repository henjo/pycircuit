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

### P4 — ~1528 lines of maintained dead code

`circuit_cna.py` (1377 lines) and `elements_cna.py` (151) are imported by
nothing — not the package, tests, docs or benchmarks — yet `circuit_cna.py`
duplicates eight core class names including `Circuit` and `SubCircuit`. Git
shows both receiving the Python-3 port, the NumPy 2.0 fixes and the toolkit
refactor. Effort is being spent maintaining code nothing runs, and duplicated
core classes invite an edit landing in the wrong file. Establish what CNA was
meant to be, then delete or revive.

### P5 — six symbolic paths and no user-facing story

`symbolic`, `symbolic_poly`, `soe`, `ginac`, `ddd`, `symengine`. Two are
honestly labelled experimental negative results, which is good practice and
should stay. But the answer to *"which should I use?"* lives in commit messages
and scattered doc pages, not in the API or one table.

`soe` is the odd one out: a module-level `solve_soe()` used by tests, benchmarks
and docs, reachable from no analysis. It is a research prototype promoted to
package code without an integration path — compare DDD, which got
`ac_solution`, a result object and an analysis family.

### P6 — two extension conventions, built opposite ways

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

### P7 — where the abstraction actually leaks: `G` is not always ∂i/∂x

The promise in §2 is that one element definition serves every toolkit. It holds
almost everywhere, and the place it fails is worth knowing because the failure is
**silent**.

An element may supply its conductance matrix two ways: hand-written `G()`, or
autodiff of `eval_i_pure`. Nothing checks that these agree. For
`VCVS_limited` they do not — measured on a properly initialised element at
`x = [0.1, 0, 0.3, 0, 0.05]`, `g = 2`, `level = 0.5`:

```
hand-written G, branch row : [ 0.      0.     -1.   1.   0. ]
d i / d x,      branch row : [ 0.8487 -0.8487 -1.   1.   0. ]
```

Two distinct causes, both pre-dating the JAX work (both present in
`py3-support`):

1. `i()` computes `v_out - f'(v_in)·f(v_in)` and never uses the gain `g`, while
   `G()` stamps `g_limit·g`. The `f'·f` product looks wrong for a VCVS on its
   face — you would expect `f(g·v_in)` — but the intended semantics are not
   recorded anywhere, so it is **left alone** pending someone who knows.
2. `Tanh.fprime` returned a hard `0` — see below. **Fixed.**

The general lesson for the architecture: an element with both a manual Jacobian
and a pure form has an invariant nothing enforces. A test that finite-differences
`i()` and compares against `G()` for every nonlinear element would have caught
this, and is the single highest-value test this codebase does not have.

### P8 — a 17-year-old shadowed derivative *(fixed)*

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

This one also shows why P7's proposed test matters: no existing test could
distinguish a correct Jacobian from a zero one, because Newton still converged.

**Unreachable code is worth grepping for generally** — an AST scan found four
`return`-followed-by-code sites. `elements.py:1840` is another live one: a
transmission-line `u()` returns before its entire transient branch, so a TLine
contributes no transient stimulus. Not investigated further.

### P9 — `_symbolic` is missing primitives `_numeric` has

`_numeric` provides `tanh`; `_symbolic` provides neither `tanh` nor `power`. So
`func.Tanh` — and any element built on it — cannot run symbolically at all,
despite nothing in its design being numeric-specific. The backends have drifted
apart without anything detecting it.

A test that asserts the backend modules expose the same primitive names would
pin this cheaply.

### P10 — smaller, verified

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

### P11 — test coverage is shaped inversely to risk

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
