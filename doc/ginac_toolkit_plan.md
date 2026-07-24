# GiNaC toolkit — implementation plan

Concrete plan for a `GinacToolkit` that solves the MNA system fraction-free with
GiNaC's `normal()` (multivariate-GCD cancellation) behind the identical toolkit
interface. Rationale and the negative symengine result are in
`ginac_toolkit_design.md`.

## The one question that gates everything

**Does GiNaC `normal()` actually beat sympy `DomainMatrix.solve_den` on
numeric-components + symbolic-`s` circuits?** symengine failed this (its
fraction-free primitives don't cancel to a compact determinant). GiNaC's
`normal()` is purpose-built for this and *should* win, but it is unproven in
this repo. **Phase 3's benchmark is a hard go/no-go gate** — do not wire
`GinacToolkit` into anything until it demonstrably beats `solve_den` on the
`symbolic_poly.rst` ladder.

## Environment prerequisites (one-time, needs you)

sudo needs a password here, so run these yourself (e.g. via `! <cmd>`):

```
sudo apt install libginac-dev        # 1.8.10-1, pulls libcln too
<venv>/bin/pip install pybind11
pkg-config --modversion ginac        # verify: prints 1.8.10
```

`g++`, `pkg-config`, and Python 3.14 headers are already present; `cmake` is not
needed (build via setuptools).

## Phase 1 — Minimal C++ extension (tiny surface)

One pybind11 module `pycircuit/circuit/_ginac_ext.cpp` exposing a single
function; keep all GiNaC contact in C++.

```cpp
// solve_numden(n, entries, rhs, symbols) -> (numerators[str], denominator[str])
//   entries : n*n row-major GiNaC-parseable strings
//   rhs     : n strings
//   symbols : names to register (e.g. "s", "R", "C")
// Builds a symbol table + GiNaC::parser, assembles a matrix, then:
//   ex det = A.determinant();                 // GiNaC picks Bareiss for symbolic
//   matrix num = A.solve(x, b);               // or Cramer via cofactors
//   for each i: result_num[i] = (num[i]*det).normal();   // <-- the GiNaC win
//   result_den = det.normal();
// Returns the ex's serialized with `ostream << ex` (parseable back).
```

Notes:
- `normal()` is the whole point — it returns a canonical `numer/denom` with the
  GCD cancelled; apply it once to each numerator and the denominator.
- Prefer `A.solve()` with `bareiss`/`gauss` and reuse its determinant if
  exposed; otherwise Cramer (`n+1` determinants) is fine for a first cut since
  only the requested rows are needed.
- Serialize with default `operator<<` (uses `^` for powers, GiNaC function
  names). The Python side maps back.

## Phase 2 — Python bridge + toolkit

- `pycircuit/circuit/_ginac.py`: wraps the extension. `linearsolver_num_den(A, b)`:
  1. sympy → GiNaC strings: `sympy.srepr`-free plain `str()`, then `**`→`^`.
     Collect `A.free_symbols | b.free_symbols` for the symbol table and for the
     assumption re-map (below). Floats: pass as decimals; GiNaC reads them as
     `numeric` (CLN) — acceptable, or rationalize first for exactness.
  2. call `_ginac_ext.solve_numden(...)`.
  3. GiNaC strings → sympy: `sympy.sympify(s.replace('^','**'))`, then
     `xreplace` the plain symbols back to the originals **with assumptions**
     (the exact bug fixed in `SymengineToolkit` — `Symbol('s', imaginary=True)`
     etc. must survive).
- `GinacToolkit(SymbolicPolyToolkit)` in `toolkit.py` overriding only
  `linearsolver_num_den` to call `_ginac.py`; on any parse/ImportError, fall
  back to `super()` (sympy fraction-free). Everything else (noise_psd, det,
  TransferFunction/poles) is inherited unchanged.
- Register an optional `ginac_toolkit` singleton (None if the extension is not
  importable), mirroring `symengine_toolkit`/`jaxtoolkit`.

## Phase 3 — Build integration + the go/no-go benchmark

- `setup.py`: add an **optional** `Extension('pycircuit.circuit._ginac_ext', ...)`
  built with pybind11, guarded by a `pkg-config --exists ginac` probe. If GiNaC
  is absent the extension is skipped and `ginac_toolkit` is `None` — install
  stays green everywhere. (`pybind11.get_include()` for headers;
  `pkg-config --cflags --libs ginac` for flags.)
- **Benchmark gate**: extend the `symbolic_poly.rst` live benchmark to a third
  column, `ginac_toolkit`, on the numeric+`s` ladder, comparing time and
  `count_ops`. GiNaC must be clearly faster (especially as N grows) to proceed.

## Phase 4 — Validation

- Correctness (skip if extension unavailable, like `test_symengine_toolkit.py`):
  RC poles `{-1/(C*R): 1}` == `symbolic_poly`; the `symbolic_poly.rst` noise
  path; example2/example5.
- Robustness: transcendental entries (sympy `EX` domain) → fall back to sympy;
  1×1 and singular systems; imaginary-`s` conjugation in `noise_psd`.

## Phase 5 — Only if the gate passes

Make `GinacToolkit` the recommended backend for large numeric-components +
symbolic-`s` circuits; document in `symbolic_poly.rst` and the design doc.
Otherwise: keep `symbolic_poly` as the default and record GiNaC's result too
(a second negative result would still be worth documenting).

## Risks

1. **normal() may still not beat solve_den** — the whole gate. Mitigation: prove
   it in Phase 3 before any integration.
2. **String round-trip fidelity** — powers (`^`↔`**`), function names, floats,
   and symbol assumptions. Mitigation: the assumption re-map + a targeted
   round-trip test; consider building the GiNaC `ex` from a sympy walk later if
   strings prove lossy.
3. **Build fragility** — pybind11 + pkg-config ginac across environments.
   Mitigation: guard the extension; graceful `None` fallback; never a hard dep.
4. **Float vs exact** — GiNaC `numeric` from doubles vs sympy Rationals; keep
   the `get_exact()` discipline `symbolic_poly` already uses.

## Effort estimate

Phase 1–2: ~a day (small C++ + bridge). Phase 3 build+benchmark: ~half a day.
The C++ surface is intentionally one function, so most risk is in the bridge and
the benchmark outcome, not the GiNaC calls.
