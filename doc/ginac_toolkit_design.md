# GiNaC / compiled-CAS symbolic toolkit — design & status

Goal: a faster symbolic backend for pycircuit's symbolic circuit analysis, using
a compiled C++ CAS for the linear solve instead of pure sympy. Plugs into the
existing toolkit architecture (`pycircuit/circuit/toolkit.py`) exactly like
`symbolic_poly`, so `numeric`/`symbolic`/`symbolic_poly` stay untouched.

## Decision gate (2026-07) — symengine `LUsolve` does NOT win (corrected)

A first "gate" timed symengine `LUsolve` followed by *numeric* evaluation and
looked ~10× faster than sympy `solve_den`. **That was measuring the wrong
thing** — numeric evaluation, not symbolic `N(s)/D(s)` extraction. Re-run
end-to-end (AC solve + transfer function on numeric-components + symbolic-`s`
ladders), the `SymengineToolkit` is *slower*, and worsens with N (N=10:
symengine ~1.4 s vs sympy `symbolic_poly` ~0.04 s).

Why:

- symengine's `LUsolve`/`det` are **not fraction-free** — the denominator comes
  back as a large *unexpanded nested expression with divisions*
  (`/(1e-9*s+101)…`), exactly the expression swell `symbolic_poly`'s
  fraction-free `DomainMatrix.solve_den` was built to avoid. sympy then has to
  cancel that mess downstream.
- symengine *does* expose fraction-free `FFLU`/`FFLDU`, but they don't cancel to
  a compact determinant here: `FFLU` is instant, yet the last pivot (the
  determinant) **expands to ~125 M characters at N=8** and explodes beyond.

**Conclusion:** the lever is **fraction-free elimination with real polynomial
cancellation**, which sympy's polynomial-domain `solve_den` already does well
(see the live benchmark in `doc/src/circuit/symbolic_poly.rst`: `symbolic_poly`
is ~2–3× faster than stock `symbolic` and ~30× smaller expressions). A compiled
win requires **GiNaC's `normal()`** (multivariate-GCD cancellation), which
symengine lacks and which is not installable in this environment.

`SymengineToolkit` is therefore kept as a **correct-but-experimental plug-in
point** (relabelled in its docstring), not a speed win.

## Environment reality

- **symengine 0.14.1** is installable via pip and present. It is a fast C++ CAS
  with a sympy-compatible Python wrapper and `DenseMatrix` with LU / fraction-free
  LU — the recommended validation vehicle.
- **GiNaC** native lib is *not* installed (no `libginac`/`ginsh`); installing
  `libginac-dev` needs root/apt, and Python bindings (pyginac/swig) are stale. A
  GiNaC C++ extension (pybind11) is the last-increment option once the lib is
  available.

## Recommended path (Option A — hybrid)

Follow the `SymbolicPolyToolkit` pattern: a new `SymengineToolkit`
(later `GinacToolkit`) that **inherits `SymbolicToolkit`** and overrides only the
solve, keeping all matrix assembly in sympy object-arrays:

1. `linearsolver_num_den(A, b) -> (numerators, denominator)`:
   - Bridge sympy→symengine (`symengine.sympify`, `symengine.Matrix`).
   - Fraction-free / LU solve in symengine; denominator = network determinant.
   - Bridge back symengine→sympy (`sympy.sympify(expr)` / `expr._sympy_()`).
2. `noise_psd` and `det` inherit the shared-denominator machinery unchanged
   (they already build on `linearsolver_num_den`).
3. Register a singleton `symengine_tk` (later `ginac`) in `toolkit.py`; opt-in via
   `toolkit=` / `use_toolkit(...)`. Reuses the `num_den` capability + the whole
   `TransferFunction` / `poles()` plumbing.

Validate against `symbolic_poly` on the RC example in the `SymbolicPolyToolkit`
docstring (poles `{-1/(C*R): 1}`) and example2/example5 noise, then benchmark on
larger circuits.

## Known trade-offs (Option A)

The bridge round-trips at the solve boundary only (negligible next to the
solve). Assumptions on symbols (`imaginary=True` for `s`, `positive=True`) don't
survive a naive string round-trip — establish them on both sides. Post-solve
simplification (`cancel`/`poles`) is still sympy; a deeper Option B (native
objects end-to-end) is only worth it if profiling shows assembly or post-solve
simplification dominates. See the earlier `sympy_review.md` notes.

## Status & next step

Done: `SymengineToolkit.linearsolver_num_den` + optional `symengine_toolkit`
singleton + regression test (RC poles vs `symbolic_poly`). It is correct but,
per the corrected gate above, **not faster** — so it stays experimental.

The real next step is a **GiNaC** backend using `normal()` for fraction-free
cancellation, behind the identical interface:

1. Get `libginac-dev` available (needs root/apt in this environment).
2. Minimal pybind11 extension: `solve_numden(matrix, rhs, s) -> (N, D)` that
   does the solve + `normal()` in C++ and returns compact polynomials.
3. Benchmark GiNaC `normal()` vs sympy `solve_den` on the ladder above; only
   then does a `GinacToolkit` replace `symbolic_poly` for large numeric+`s`
   circuits.

Until then, `symbolic_poly` (sympy polynomial-domain, fraction-free) is the
fastest available symbolic backend for the numeric-components + symbolic-`s`
regime.
