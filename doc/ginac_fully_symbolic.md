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

### 1. Go native — kill the string bridge (Option B)

For large symbolic results the `sympify` round-trip dominates (≈1.4 s of 1.7 s
above). Keep the result as a GiNaC `ex` inside a lightweight Python wrapper and
convert to sympy only lazily / for small pieces. Prerequisite for everything
below.

### 2. Use GiNaC where it is actually faster — *after* the solve

The symbolic form is the same size, but GiNaC manipulates and (especially)
*evaluates* it much faster than sympy:

- **Fast numeric evaluation** of a symbolic result (substitute component values,
  sweep frequency). GiNaC's `subs`+`evalf` over its shared DAG is far faster than
  sympy for a big expression — the killer feature for "derive once symbolically,
  evaluate numerically many times" (Bode / Monte-Carlo over a symbolic TF).
- `factor()`, `collect()`, `series()` (Laurent expansion for pole/zero
  approximation) — much faster on large multivariate.
- GiNaC expressions are shared DAGs, so common subexpressions cost memory once.

### 3. Don't force the canonical `N/D` — Sequence of Expressions (SoE)  *(prototyped)*

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

## Bottom line

GiNaC will not make the exact fully-symbolic transfer function smaller — that is
an inherent exponential. Its value is (a) a native representation to escape the
string bridge, (b) fast manipulation/evaluation of large symbolic results, and
(c) a fast arithmetic backend for compact (SoE) and compressed (DDD) / approximate
methods that tackle the exponential head-on. The current Option-A, det-based
`GinacToolkit` captures none of these; the productive next steps are a
GiNaC-`ex`-backed result object with fast numeric evaluation (#1/#2), then SoE
(#3, prototyped) and DDD (#4).
