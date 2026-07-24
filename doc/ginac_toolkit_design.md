# GiNaC / compiled-CAS symbolic toolkit — design & status

Goal: a faster symbolic backend for pycircuit's symbolic circuit analysis, using
a compiled C++ CAS for the linear solve instead of pure sympy. Plugs into the
existing toolkit architecture (`pycircuit/circuit/toolkit.py`) exactly like
`symbolic_poly`, so `numeric`/`symbolic`/`symbolic_poly` stay untouched.

## Decision gate (2026-07) — PASSED

Benchmark: solve a symbolic-`s` RC/RLC ladder MNA system (numeric component
values + one symbolic generator `s`) — the sweet spot for fraction-free
solving. Compared sympy `DomainMatrix.solve_den` (current `symbolic_poly` path)
vs **symengine** `LUsolve`:

| N (nodes) | sympy solve_den | symengine | correct? |
|-----------|-----------------|-----------|----------|
| 4  | 0.016 s | 0.002 s | yes (machine precision) |
| 6  | 0.022 s | 0.002 s | yes |
| 8  | 0.036 s | 0.003 s | yes |
| 10 | 0.062 s | 0.004 s | yes |

~10–15× faster, correct to machine precision, gap grows with N. A compiled
backend is worth it.

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

## Next step

Implement the minimal `SymengineToolkit.linearsolver_num_den` + a singleton and
one regression test (RC transfer-function poles vs `symbolic_poly`).
