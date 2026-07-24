# Fully symbolic analysis with GiNaC — what actually helps

Context: the numeric-components + symbolic-`s` regime is covered in
`ginac_toolkit_design.md` / `ginac_toolkit_plan.md` (there `symbolic_poly` wins;
GiNaC's `normal()` hits a coefficient-magnitude blow-up). This note is about the
*fully symbolic* regime (all `R`, `C`, ... are symbols — many generators) and
where GiNaC can and cannot help.

## Measured baseline (AC solve of an N-section RC ladder, all symbolic)

| N | `symbolic` (LUsolve) | `symbolic_poly` | `ginac_toolkit` |
|---|----------------------|-----------------|-----------------|
| 4 | 0.04 s / 126 ops     | 0.07 s / 187    | 0.10 s / 187    |
| 6 | 0.10 s / 350 ops     | 0.37 s / 1752   | 1.70 s / 1752   |

Two things to read off:

1. **No coefficient blow-up / no hang** here — that failure mode was specific to
   *numeric* values (`1e-9`). Fully symbolic circuits have no large rationals.
2. But GiNaC still does not win. `GinacToolkit` inherits the fraction-free
   (`det`-based) `linearsolver_num_den`, which **expands the multivariate
   determinant** (1752 ops), while stock `LUsolve` keeps a compact *nested,
   uncanonicalised* form (350 ops). And of `ginac_toolkit`'s 1.7 s, ~1.4 s is
   the sympy↔GiNaC **string bridge** of the large result.

A direct GiNaC check confirms the trap: `GiNaC::matrix::solve()` already
canonicalises (its raw result has the same tree size as `.normal()`, e.g. 2203
nodes at N=6), so there is no free "compact form" to grab. The compact form has
to be produced deliberately.

## The exponential is inherent

A fully symbolic transfer function is a ratio of network determinants, and a
symbolic determinant is a sum over spanning-tree / topological products —
exponentially many terms in general. **No direct solver escapes this**;
`normal()` produces exactly that canonical (large) form. So "use GiNaC" is not by
itself a fix. What helps:

## Improvements, most valuable first

### 1. Go native — kill the string bridge (Option B) — **DONE**

For large symbolic results the `sympify` round-trip dominates (≈1.4 s of 1.7 s
above). Keep the result as a GiNaC `ex` inside a lightweight Python wrapper and
convert to sympy only lazily / for small pieces. Prerequisite for everything
below.

**Implemented.** `_ginac_ext.GinacResult` holds the numerators and shared
denominator as native GiNaC `ex`; `_ginac.solve_native(A, b)` (and
`GinacToolkit.solve_native`) return a light `_ginac.GinacResult` wrapper. The
determinant-sized full result is never round-tripped through sympy — you request
only the compact piece you need, and the cancellation happens *in GiNaC* before
that small piece crosses back:

- `.tf(i, j)` — transfer function `x_i/x_j = num_i/num_j`; the shared
  denominator cancels natively, so the result stays compact,
- `.denominator()` — the network determinant, for poles,
- `.component(i)` / `.numerator(i)` — the full (large) pieces, still available,
- `.eval_tf(i, j, params)` / `.eval_component(i, params)` — numeric sweeps that
  compose the compact `.tf`/`.component` with `eval_sweep` (#2).

Symbol assumptions (`Symbol('s', imaginary=True)`) are re-mapped onto the
returned sympy, as in `linearsolver_num_den`.

Measured on a fully-symbolic RC ladder, extracting a transfer function
`v_{N-1}/v_0` — native `solve_native + tf` vs the eager `linearsolver_num_den`
+ sympy `cancel` (identical result):

| N | eager (sympify full, then cancel) | native (`solve_native`+`tf`) | speed-up |
|---|-----------------------------------|------------------------------|----------|
| 4 | 0.04 s | 0.01 s | 3.3× |
| 6 | 0.47 s | 0.17 s | 2.8× |
| 8 | 10.5 s | 3.8 s  | 2.8× |

The gap widens with size; the compact `tf` string is also ~25× smaller than the
materialised full solution. This unlocks #2 (native manipulation/eval) and gives
#3/#4 a native result object to build on.

*(These numbers regenerate live on every docs build — see the Sphinx page
`doc/src/circuit/ginac_native.rst`; the table above is a point-in-time record.)*

### 2. Use GiNaC where it is actually faster — *after* the solve — **DONE**

The symbolic form is the same size, but GiNaC manipulates and (especially)
*evaluates* it much faster than sympy:

- **Fast numeric evaluation** of a symbolic result (substitute component values,
  sweep frequency). GiNaC's `subs`+`evalf` over its shared DAG is far faster than
  sympy for a big expression — the killer feature for "derive once symbolically,
  evaluate numerically many times" (Bode / Monte-Carlo over a symbolic TF).
- `factor()`, `collect()`, `series()` (Laurent expansion for pole/zero
  approximation) — much faster on large multivariate.
- GiNaC expressions are shared DAGs, so common subexpressions cost memory once.

**Implemented.** Numeric evaluation is `eval_sweep` (native `compile_ex`, see the
`eval_sweep` docs). Native symbolic manipulation, running on the GiNaC `ex` and
returning only the compact structured result:

- `_ginac.poly_coeffs(expr, var)` and `GinacResult.tf_coeffs(i, j, var)` /
  `.denominator_coeffs(var)` — collect a rational into `N(var)/D(var)`
  coefficient vectors (the canonical transfer-function form for poles/zeros or a
  `scipy.signal` handoff). `normal()` + collection run in GiNaC; the common
  factor is cancelled natively.
- `_ginac.series(expr, var, order, point=0)` — truncated Laurent/Taylor series
  for low-/high-frequency and pole/zero approximation.

Extracting the `N(s)/D(s)` coefficients of a fully-symbolic RC-ladder transfer
function, native `tf_coeffs` vs sympy `Poly` (both starting from the native
`tf`): ~4.7× (N=4), ~4× (N=6). The gap grows with expression size; this is the
"manipulate in GiNaC, only the small structured result crosses back" win.

*(Also regenerated live on every docs build — see
`doc/src/circuit/ginac_native.rst`.)*

### 3. Don't force the canonical `N/D` — Sequence of Expressions (SoE)  — **DONE**

**Implemented** as a result object: `pycircuit.circuit.soe.solve_soe(A, b)`
returns an `SoESolution` holding the ordered `(symbol, expr)` assignments and the
node solution. `.eval(params)` resolves the assignments *in order* over a
parameter sweep (vectorised numpy) and `.eval_tf(i, j, params)` gives a swept
transfer function; `.to_ratio(i, j)` inlines to a single sympy expression for the
GiNaC `N(s)/D(s)` tools.

**Key finding:** the SoE's compactness is entirely from *sharing* -- each
intermediate symbol is referenced several times but stored once. The assignment
count stays linear for a ladder (72/156/240/324 ops at N=4/8/12/16), but
*inlining* (`to_ratio`) destroys the sharing and regrows super-linearly
(213/1235/3633/7983 ops). So the scalable path is **sequential evaluation**, not
compiling one expression -- and, notably, a single GiNaC `compile_ex` cannot
help here (it does not CSE, so it would re-expand the shared form). `to_ratio` ->
`poly_coeffs` is therefore a *small-circuit* bridge only: for a fully symbolic
ladder the coefficient extraction hits the inherent multivariate exponential
(the very blow-up SoE avoids). Live demo + check: `doc/src/circuit/soe_symbolic.rst`.

Original prototype notes:

For fully symbolic circuits, poles have no closed form anyway (degree ≥ 5), so
the exploded `N/D` is rarely worth its exponential size. Represent the solution
as a **sequence of expressions**: symbolic Gaussian elimination whose
intermediate results are kept as fresh symbols. This is `O(n^3)` assignments
(and, for a banded/tridiagonal ladder, linear), each small, and directly
evaluable.

`doc/prototypes/soe_symbolic.py` demonstrates it (sympy, to show the
representation). Exploded solve grows exponentially; SoE grows polynomially:

| N | direct (exploded) ops | SoE total ops | SoE assignments |
|---|-----------------------|---------------|-----------------|
| 4  | 49    | 72  | 21 |
| 6  | 141   | 114 | 35 |
| 8  | 382   | 156 | 49 |
| 10 | 1013  | 198 | 63 |
| 12 | 2665  | 240 | 77 |

All correct (numeric evaluation of the SoE matches the direct solve). By N=12
the SoE is ~11× smaller and the gap widens exponentially. In production the SoE
arithmetic and — most importantly — the numeric evaluation of the assignment
list would run through GiNaC (fast compiled arithmetic + DAG sharing), tying #3
to #1/#2.

### 4. The real scalability lever is algorithmic (GiNaC as the backend)

- **DDD — Determinant Decision Diagrams** (Shi & Tan): a compressed, canonical
  representation of the symbolic determinant. pycircuit has an orphaned
  `ddd.py`; reviving it with GiNaC doing the per-node polynomial arithmetic is
  the proper path to *large* symbolic circuits.
- **Approximate symbolic analysis** — keep only dominant terms (the
  `symbolicapprox.py` direction). Exact fully symbolic is exponential by nature.
- **Hierarchical / topological decomposition** — solve subcircuits, exploit
  sparsity / connected components to cut the term count *before* GiNaC
  arithmetic.

## Numeric-coefficient regime: O(1) scaling — **DONE**

The other failure mode was *numeric* components: values like `1e-9` become
rationals whose products across the determinant explode into astronomically
large CLN integers that stall GiNaC around dimension 16 (this is separate from
the fully-symbolic exponential). The fix is **frequency + magnitude scaling** so
GiNaC sees O(1) coefficients: for a numeric AC system `G + s C`, substitute
`s = w0*sigma` and divide by an amplitude `amp`, both chosen as **powers of ten**
so the scaled entries are exact small rationals (irrational geometric-mean
factors reintroduce the blow-up). Solve in `sigma`, then back-substitute
`sigma = s/w0` — a cheap coefficient rescale on the already-compact result.

Implemented in `pycircuit.circuit._ginac` (`detect_freq`, `_frequency_scale`,
auto-applied by `linearsolver_num_den`) and wired through `GinacToolkit`, which
now drops the `ginac_max_dim` guard whenever the system is scalable. Result: the
numeric-AC solve is fast far past the old cliff (e.g. dim 24 in ~0.4 s, machine
precision). `AC(..., toolkit=ginac_toolkit)` therefore handles large numeric
ladders directly, and `TransferFunction.as_num_den` uses GiNaC `poly_coeffs` for
the **symbolic** coefficient collection (poles/zeros) while numeric transfer
functions stay on sympy/`numpy` (a scaled numeric result spans ~1e0…1e-90, where
GiNaC's `limit_denominator` would drop tiny high-order terms).

## Bottom line

GiNaC will not make the exact fully-symbolic transfer function smaller — that is
an inherent exponential. Its value is (a) a native representation to escape the
string bridge, (b) fast manipulation/evaluation of large symbolic results, and
(c) a fast arithmetic backend for compact (SoE) and compressed (DDD) / approximate
methods that tackle the exponential head-on. The current Option-A, det-based
`GinacToolkit` captures none of these; the productive next steps are a
GiNaC-`ex`-backed result object with fast numeric evaluation (#1/#2), then SoE
(#3, prototyped) and DDD (#4).
