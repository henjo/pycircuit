# The backend spike (roadmap section 14.5), measured 2026-08-26

A measurement, not a build.  `benchmarks/backend_spike.py` reproduces
every row below from the live `PspMosLongChannel` class; nothing here
ships as a backend and `hdl.py` was not touched.

Branch `cna-jax-vectorization`, HEAD `37c9f72`.  Python 3.14.4,
numpy 2.5.2, gcc 15.2.0, clang (Ubuntu), 24 cores.  Installed into the
venv for this spike: numba 0.67.0, llvmlite 0.49.0, cffi 2.1.1
(pycparser 3.0 as a dependency).

**First, a record defect.** `doc/hdl_roadmap_260824.md` at HEAD no
longer contains sections 12 to 15: commit `37c9f72` ("the library's
own header table and the status row catch up") removed 1300 lines,
including section 14 -- the specification this spike runs against --
while the status table still says "see section 13/14".  The text used
here is from `8267714`, the last commit that has it.  Restore before
anything cites section 14 again.

## The answer

| candidate | G us | i us | G speedup | compile s (per class) | worst agreement over the sweep | finiteness | verdict (14.5 words) |
|---|---|---|---|---|---|---|---|
| numpy, as shipped (element path) | 17 230-17 800 | 1 260-1 290 | 1x | 0.6 warm (section 16 cache) | -- | non-finite at 96 of 4096 extreme-grid points | baseline |
| numba, default `NUMBA_OPT=3` | **JIT did not finish in 600 s** | 6.6 (packed) / 6.2 (155 args) | -- | G: > 600 s; i: 9.3 s | i: 5 ulp, 6.5e-16 rel | i: nothing lost | G: **don't** -- 14.7 "JIT cost per class exceeds what a simulation amortises" |
| numba, `NUMBA_OPT=1` | JIT did not finish in 590 s | -- | -- | > 590 s | -- | -- | same |
| numba, `NUMBA_OPT=0` | 622-643 | -- | **28x** | 68 s, +1.8 GB RSS | 2.4e-10 rel (1.28e6 ulp) at G[0,2], bias (-0.05, -0.525, 0, 0) | nothing lost | "between 10x and 50x: report and decide" -- and the recommendation is don't (below) |
| C, gcc -O2, bitwise-faithful (`-fno-builtin-pow`, every `**` is libm `pow`) | **75.5** | **7.7** | **229x** | 22-32 s (G), 1.5 s (i) | **0 ulp, 0.0 abs, every entry of G and i at all 4985 points** | nothing lost, nothing gained | **worth building out (>= 50x, full-sweep agreement)** |
| C, gcc -O2, `x*x` for `**2` (`pow` elsewhere) | **14.0** | 3.4 | **1235x** | 19 s | 2.4e-10 rel at the same entry; i: 5 ulp | nothing lost | worth building out IF 2.4e-10 is accepted as agreement (see "the one ulp") |
| C, gcc -O2, builtins folded (gcc's default) | 9.6 | 3.25 | 1800x | 23 s | 2.4e-10, same origin | nothing lost | as above |
| C, clang -O2, faithful | 79.7 | 7.8 | 216x | 15 s | (same code path as gcc faithful) | -- | -- |
| C, gcc -O1 / -O0, faithful | 78.5 / 186 | 7.9 / 9.8 | 220x / 93x | 20 s / 4.7 s | -- | -- | -- |
| llvmlite | not built | | | | | | skipped: C answers the question, and numba already shows what LLVM's optimiser does with this function (below) |

Per-call numbers for C are from an idle re-timing (min of five 2000-call
runs; the harness's first-pass numbers were taken with five jobs
sharing the machine and are up to 5x worse).  Baseline and numba
numbers are from the harness.  Speedups are against the element path
measured in the same process (17.2-17.8 ms; section 13 quoted 17.3).

The sweep: `Vd, Vg, Vs, Vb` each over `{-1e30, -100, -1, 0, 0.7, 1.2,
100, 1e30}`, all 4096 combinations, plus the 889 bias points of the
twelve `psp_reference.py` sweeps (nmos and pmos alike, fed to this
nmos class) -- 4985 points, 36 entries of G and 6 of i each.  The
numpy path raised nothing and is finite at every reference point; at
96 of the 4096 extreme points it returns NaN in the same 20 entries of
G, always with `Vb = +-1e30` together with `Vs` or `Vd` at `+-1e30`.
That is the model, not a backend, and no candidate changed it: over
the whole sweep no candidate lost a single entry the numpy path had
finite, and none gained one.

## What each candidate needed changed, and why the semantics survive

The property section 14.6 protects: both arms of every conditional
are computed, so the discarded arm's overflow or NaN is a fact about
the arm and not hidden by the selection.  In the emitted Python this
is not a property of `numpy.where` but of the CALL: `where(c, a, b)`
receives `a` and `b` already evaluated.  So a translation keeps the
property exactly when the arms stay call ARGUMENTS.  Both candidates
below do that; neither uses `a if c else b` inline (which would skip
an arm) or `a*c + b*(1-c)` (which turns a losing infinite arm into
NaN -- `_maxc_numpy`'s docstring records why).

Values cannot differ from this: on scalars `where` returns one of its
two arguments, so whichever way the arms are computed the selected one
is the same double.  The only observable a backend could lose is the
floating-point exception flags of the losing arm, which the numpy path
renders as warnings and neither backend reports.

**Common to both (mechanical, on the AST of `_src`, no arithmetic
changed):**

1. `numpy.logical_or.reduce((a, b))` and `numpy.logical_and.reduce(...)`
   (32 + 40 sites in G) become nested binary `numpy.logical_or(a, b)`.
   The tuple was built -- every element evaluated -- before `reduce`;
   the nested form evaluates the same elements in the same order.
   This is sympy's `NumPyPrinter` printing `Or`/`And` past the custom
   printer, the same escape section 8 recorded for `Min`/`Max`.
2. The `return [[...], ...]` of mixed floats and literal `0` becomes
   stores into a preallocated array.  Each cell expression is
   unparsed verbatim.

**numba additionally:**

3. `numpy.where(c, a, b)` becomes `_where(c, a, b)`, an `@njit`
   scalar select.  numba types scalar `numpy.where` as numpy does -- a
   0-d array -- and then has no `abs()` for it, so the emitted source
   does not type as it stands.  The arms remain arguments (property
   kept); numba's `if c` truth test is numpy's (nonzero, and NaN, are
   true).
4. `maxc`, `minc`, `_step`, `_rdiv`, `_recip2`, `_wrapfloor` are
   re-declared as `@njit` twins of their numpy runtime forms with the
   sympy/jax dispatch guards removed -- there is nothing to dispatch on
   inside nopython mode.  Same arithmetic.
5. Parameters either as the 155-argument signature (works; ~0.4 us
   cheaper per call on `i` than packing, so marshalling is not the
   cost) or packed into one array.

**C additionally:**

3. Every primitive is a `static inline` FUNCTION: `_sel(c, a, b)` for
   `where` (arguments evaluated before the call -- C's rule, same as
   Python's), `_npmax`/`_npmin` written as numpy defines them
   (`(a >= b || a != a) ? a : b`, so NaN in either argument wins),
   `_step`, `_rdiv`, `_recip2`, `_sign` (NaN-propagating, zero for
   zero), `_land`/`_lor` on numpy's truth test; comparisons become C
   comparisons; `abs` is `fabs`; `numpy.nan`/`inf` are `NAN`/`INFINITY`.
4. Flags: `-fno-fast-math -ffp-contract=off` (no FMA contraction, no
   `-march`), and `-fno-builtin-pow` for the bitwise variant.
5. Calling convention `(const double *x, const double *p, double *out)`
   through cffi ABI mode (`dlopen`, no cffi compile step); the
   parameter array is cast once, `x` per call.  `out` is reused, as
   any C stamp would be.

Nothing was clamped, reordered, simplified or constant-folded by hand.
Association and precedence follow the Python AST node for node.

## The one ulp, and what it says about "agreement"

The trace at the worst point (`--trace` reproduces it; `tmp/trace.py`
and `tmp/dissect.py` in the job dir) shows every one of the ~740 VALUE
intermediates bitwise identical between numpy and C -- `exp`, `log`,
`sqrt` all agree exactly on this machine -- and the first divergence
inside one derivative line, `_d__v217_x0_2`, at `den**2`.  numpy's
scalar `float64 ** 2` is glibc `pow(x, 2.0)` (measured: equal to
`math.pow` on 20 000 random doubles, different from `x*x` in 13 of
them; `**3`, `**-1.0`, `**0.5` are all `pow`).  glibc's `pow` is
within 0.52 ulp, not correctly rounded; `x*x` is.  gcc folds
`pow(x, 2.0)` to `x*x` at -O1 and above without fast-math, and numba
lowers `x**2` to a multiply -- so both compiled paths were the MORE
accurate side.  That one ulp then went through `-2*den**2*recip2 + 1`
(a cancellation to -1) to 30 ulp in that line and to 2.4e-10 relative
in `G[0,2]` (a 1.7e-14 entry) through the gate-current terms.

Consequences:

- "Bitwise agreement with the numpy path" is achievable in C and costs
  **5.4x** (75.5 vs 14.0 us): 5618 of the 6354 powers in G are `**2`,
  and a real `pow` call is ~11 ns each.  It is not achievable in numba
  without a further rewrite (untested whether `numpy.power(a, 2.0)`
  survives LLVM's folding).
- The 2.4e-10 is not a backend error; it is the CONDITIONING of the
  derivative chain at that bias to a one-ulp perturbation of one
  square.  Any future change of libm, numpy's scalar power, or the
  CPU's exp implementation moves the numpy path by the same amount.
  A test that pins G to 1 ulp against numpy would be pinning glibc's
  0.52-ulp `pow` error, not the model.
- Section 14.6's stated risk ("a different float path can move where
  `safe_*` overflow") did not materialise anywhere in 4985 points; the
  risk that did was the opposite one -- a more accurate path
  disagreeing through cancellation.  Both are the same lesson:
  measure at the extremes, but ALSO trace the first differing
  intermediate, because the worst entry is far from the cause.

## numba, plainly

On `i` numba is fine: 9 s JIT, 6.6 us, 191x.  On `G` it is not:
LLVM's optimiser at numba's default level did not return in 10 minutes
on a 2700-statement straight-line function (the working set reached
1.8 GB), at `-O1` equivalent it did not return in 590 s either, and at
`NUMBA_OPT=0` it took 68 s to produce code that runs at 622 us -- 28x,
inside the "report and decide" band and 45x slower than the C the
system compiler produced in 20 s.  gcc on the same function: 23 s at
-O2, 4.7 s at -O0 with code (186 us) still 3.3x faster than numba's.
Section 14.5's premise for candidate 1 -- "no new code generator at
all, the chained path already emits its ideal input" -- was wrong in
both halves: the source needed three rewrites to type, and the input
is not one LLVM's pass pipeline handles at this size.  Section 14.7
names the outcome: the JIT cost per class exceeds what a simulation
amortises.  Decision: don't.

llvmlite was skipped for the same reason: it is numba's LLVM without
numba's type inference, and the LLVM side is the part that did not
finish.  It would only be worth revisiting with an explicitly reduced
pass pipeline, which is more surface than emitting C.

## What building the winner would take

C via cffi, bitwise-faithful by default with the `x*x` form as a
documented option.  Files and the replacement point:

- **`pycircuit/circuit/_hdl_cbackend.py`** (new): the transpiler.  The
  spike works from the AST of `_src`; the build should instead be a
  second printer beside `_ChainPrinter` -- a `_CChainPrinter` that
  prints the SAME sympy expressions `_chain_compile` already walks --
  because `_src` is downstream of sympy and re-parsing it is an extra
  fragile layer.  The C prelude of the spike (the ten helpers) is the
  kernel's runtime contract in C and belongs beside `_KERNEL_NUMPY`
  as `_KERNEL_C`, one definition per primitive, the way the numpy
  forms are kept beside their sympy classes so they cannot drift.
- **`_chain_compile`** grows a `backend=` choice: after `_chain_prune`
  and the forward-accumulation loop it has `defs`, the Jacobian
  expressions and `args`; the numpy printer emits `_src` as today (and
  it MUST keep doing so -- `explain()`, `_hdl_cache`, `check_jacobians`
  and every test read it), and the C printer emits `G.c` from the same
  lists.  The returned callable keeps `_src` (the numpy text, the
  reference) and gains `_csrc`.  Signature `(x, p, out)`; `_args_of`
  packs the parameter vector into a float array once per
  `update_iparv` instead of per call.
- **`_hdl_cache`**: key the `.so` by the SHA of the C source plus
  compiler version string, store it beside the entry, and on a hit
  `dlopen` instead of `exec`.  The 20-30 s compile is then paid once
  per (class, compiler) like the 41 s sympy compile is today, and
  instantiation stays 0.6 s warm.  Without a C compiler on the path
  the numpy function is used and `_hdl_cache_status` says so.
- **What must not change**: the jax path (traces the numpy source, and
  needs none of this), the eager path, PCNR's chained twins.  Only the
  `funcs['i'/'G'/...]` callables of chained classes gain a C body.
- **Tests, named for the claims**: G and i bitwise against the numpy
  path over the extreme grid and the reference biases (the spike's
  sweep, made a test); both-arms via a model whose losing arm
  overflows; a `pow` sentinel (`x**2` at a value where glibc differs
  from `x*x`) so a future compiler flag cannot silently swap the
  variant; the cache round-trip in a fresh interpreter.
- Threading: cffi releases the GIL for the call, and the function is
  reentrant with no globals, so the vector/batched solvers can call it
  from threads; the numpy path cannot.  Memory: 0.6-1 MB per `.so`.

Effort: the transpiler is ~300 lines in the spike; the printer form,
the cache and the tests are the real work -- a few days, not weeks.

## What in section 14 turned out wrong

1. **14.2 / 13: "3471 dispatched numpy calls, ~16%; the remaining 84%
   is the interpreter".**  The census of the emitted G counts ~21 400
   Python-level calls: `minc` 3833, `sqrt` 3236, `maxc` 2911, `where`
   2642, `_step` 2046, `_rdiv` 1643, `_recip2` 1132, ~2800
   comparisons, `exp` 303.  Section 13 counted only `where`/`min`/
   `max`/`exp`.  cProfile puts the kernel helpers alone (`minc`,
   `maxc`, `_step`, `_rdiv`, `_recip2` and their `isinstance` guards)
   at ~41% of G.  The conclusion -- only a backend that stops
   interpreting reaches parity -- was right; the split it rested on
   was not, and "scalar codegen is worth ~1.2x" was priced on the
   wrong count.
2. **14.5, candidate 1: "no new code generator at all".**  Three
   rewrites to type, and the JIT does not finish.
3. **14.4: the 300-1000x extrapolation "must not be repeated".**  It
   was wrong in reason (14030 operations at 1-5 ns) and right in
   range: 14.0 us for the non-faithful C is 1235x and 9.6 us with the
   compiler's own folding is 1800x.  The faithful form is 229x.  The
   section was right to forbid quoting it and the number it forbade
   was, by accident, about right.
4. **14.5: "compile/JIT time ... against 66.4 s".**  The comparison is
   now against 0.6 s warm (section 16), and the right frame is the
   cache: 20-30 s once per class per compiler is amortised the way
   the sympy compile is.  The numba figure is not amortisable at any
   frame.
5. **14.6: the overflow risk of a different float path.**  Did not
   occur; the accuracy-through-cancellation risk did.  See above.
6. **Not in section 14 but found by it:** the numpy path is non-finite
   at 96 extreme-grid points (`Vb = +-1e30` with `Vs` or `Vd` at
   `+-1e30`, 20 entries of G, always the same 20).  The safety
   primitives hold everywhere else in the grid.  For the model owner.
7. **The document itself**: sections 12-15 are absent at HEAD.

## Reproducing

    .venv/bin/python benchmarks/backend_spike.py                 # C faithful, numba, full sweep
    .venv/bin/python benchmarks/backend_spike.py --no-numba --pow2-mul     # the x*x variant
    .venv/bin/python benchmarks/backend_spike.py --no-numba --cc clang
    NUMBA_OPT=0 .venv/bin/python benchmarks/backend_spike.py --no-c --only G --numba-pack-only
    .venv/bin/python benchmarks/backend_spike.py --trace G -0.05 -0.5250000000000035 0 0

The default numba run on G will not finish; give it `--no-numba` or
`NUMBA_OPT=0` unless the point is to reproduce that.  Per-call timings
are only comparable on an idle machine.
