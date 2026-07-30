# Transient subsystem review — four lenses, 2026-07-30

**Status: review only. Nothing in this document has been applied.** Every diff sketched
here is a proposal. Findings are marked CONFIRMED (the reviewer read the code and
reproduced the behaviour) or SUSPECTED (mechanism identified, not measured end to end).

Reviewed at HEAD `a675455`, before the GBW-compensation commit `95545e5`. Four
independent reviews, each given the same scope but a different lens, and each told which
defects were already fixed so they would not re-report them:

| lens | question it was asked to answer |
|---|---|
| senior Python architect | what does this force to change outside itself, and where can two things not be kept in step? |
| circuit simulation expert | does it converge and give the right answer on circuits real users simulate? |
| numerical analysis expert | is the mathematics right, and is the step-size control theoretically sound? |
| EDA tool expert | what would a user migrating from a commercial simulator miss first? |

Matrix ordering was called out by the maintainer as a specifically missing topic and was
assigned to the EDA lens as top priority, with a supporting angle from the numerical one.

**Why four lenses rather than one pass.** Three of them independently found the same two
defects (the silent DC fallback, and the unevaluated first step) by different routes and
on different circuits. That agreement is worth more than any single report, and it is the
argument for the named-lens discipline: one reviewer asked "is this fast?" would not have
found the first-step defect, and one asked "is this correct?" would not have found that
threaded BLAS makes the solve 19x slower.

---

## 1. Cross-confirmed findings

These were reported by more than one lens, independently.

### 1.1 A failed DC operating point is silently replaced by zeros — 3 of 4 lenses

`transient.py:248-252`, and again at `:443-444` in `_solve_coupled`.

```python
try:
    dc_res = dc.solve()
    x0 = dc_res.x
except Exception:
    x0 = self.toolkit.zeros(n)      # bare except, no warning, no flag
```

Demonstrated: a circuit whose DC fails outright runs a transient to completion and
returns 724 points with a plausible `v(b) = -0.4464` — **no exception, no warning, and
nothing on the result object recording that the operating point was never found.** Every
commercial tool aborts.

Compounding it, both `dcanalysis.py:130-133` and `transient.py:169-172` wrap **`except
Exception`** into `NoConvergenceError`, so a `ZeroDivisionError` from a source model is
reported as a convergence problem and then discarded. That is a chain of three layers
between a programming error and the user.

It is also made much more likely by 1.4 below: any voltage-driven bipolar circuit fails
DC today, so this is the *normal* path for those circuits, not an edge case.

**Effort: half a day. This is the highest-value single fix in the review** and it should
not wait for anything else.

### 1.2 `reltol` does not control transient accuracy — 2 of 4 lenses, different circuits

`stepcontroller.py:38-39` accepts the first step of a run unevaluated (correctly — there
is nothing to difference), and `transient.py:272` opens the run at `dt = timestep`, i.e.
at `max_step`. So every run commits one uncontrolled backward-Euler step at maximum size,
and its error dominates everything after it.

Circuit-simulation lens, RC low-pass + 1 kHz `VSin` against the closed form:

```
reltol=1e-3  steps=270   maxerr=0.009798
reltol=1e-4  steps=343   maxerr=0.009136
reltol=1e-5  steps=626   maxerr=0.009136      7.5x the work,
reltol=1e-6  steps=1240  maxerr=0.009136      bit-identical error
reltol=1e-7  steps=2579  maxerr=0.009136
```

Numerical lens, different circuit, tend=3 ms, max_step=2e-4:

```
Euler          reltol 1e-3 -> 1e-6 : globErr 3.041e-02 both
Trap(classic)  reltol 1e-3 -> 1e-6 : globErr 3.041e-02 both
Gear2(ywr)     reltol        1e-6  : globErr 3.041e-02
```

Identical to four digits across every method and a 1000x tolerance range. Localised to
step index 0 in both reviews. One backward-Euler step of h=2e-5 gives 0.020889 against
the exact 0.011848 — error 0.009041, matching the measured 0.009136.

**The controller itself is fine.** Past the opening transient, error falls 7.48e-4 ->
4.29e-5 for reltol 1e-4 -> 1e-6, a factor 17.4 against the theoretical `tol^(2/3)` = 21.5
for a 2nd-order method. It is only the first step that is uncontrolled.

Measured fix — open at `timestep * 1e-3` (the controller doubles per step, so this costs
~10 steps):

```
              baseline                 ramped first step
Euler   steps=1726 err=9.136e-3  |  steps=1774 err=4.152e-3
Trap    steps=321  err=9.136e-3  |  steps=336  err=6.410e-4   (14x better, +5% steps)
Gear2   steps=343  err=9.136e-3  |  steps=354  err=1.876e-3
and reltol becomes a knob again: 1e-3 -> 4.48e-3 ... 1e-6 -> 1.08e-4
```

The principled version is a Hairer-style initial step from `q'`/`q''` at the operating
point. A `firststep` Parameter is the minimum.

**This matters on a dissipative circuit only because the error decays. On an oscillator,
an integrator, or anything with a marginally-stable mode it is a permanent offset.**

### 1.3 `DC.vabstol` is 1e-12 while `Transient.vabstol` is 1e-6 — 4 of 4 lenses

`dcanalysis.py:48-50` vs `transient.py:91-93`. The operating point that seeds every
transient is solved to a criterion 10^6 times tighter than any subsequent point in the
same run. The 2026-07-30 `vabstol` fix could not propagate: the two `Parameter` lists are
independent literals with no shared source.

**And the same fix had a second, unmeasured effect.** `vabstol` serves *two* roles:
Newton's x-tolerance on node rows (`transient.py:142`) and the LTE tolerance on node rows
(`:294`). The change was reasoned about purely as a step-control knob, so relaxing
1e-12 -> 1e-6 also loosened the Newton solve's node convergence by 10^6 at every time
point. That is a real behaviour change that rode along unmeasured.

**This needs a decision, not a diff:** two names (`vabstol` for Newton, a separate or
derived one for the controller), or one name with the controller's ratio made explicit
rather than coincidental. Whatever is chosen must reach DC.

### 1.4 `jaxtransient.py` is a fork, not an accelerator — 2 of 4 lenses

728 lines with its own integrators (`:35-89`), its own LTE (`:228`), its own step control
(`:337`), its own Newton (`:102`) — and still a dense `jnp.linalg.solve`. It is not a
faster backend for `Transient`; it is a second transcription of the same paper.

The consequence is not hypothetical: the Gear2 LTE defect had to be found and fixed
**twice**, and the open JAX tolerance defect is the same shape. `transient.py` contains a
*third* transcription internally (`_solve_coupled`, `:429-551`) which silently ignores
breakpoints, `uic`, `fixed_timestep`, and any injected step controller.

**Root cause, stated by the architect lens:** there is no seam between "the time-stepping
algorithm" and "the array backend", so every backend gets its own copy of the algorithm.

---

## 2. Performance — measured, prototyped, 10.5x

From the EDA lens, on the leapfrog + a `BSource` cubic (137 unknowns, 708 elements,
J = 784 nnz, 4.18% dense). cProfile over 21 accepted steps, 6.44 s total:

| | cumtime | share |
|---|---|---|
| `circuit.py:1366` `_add_element_subvectors` (i, q, u) | 4.21 s | 65% |
| `circuit.py:1290` `_add_element_submatrices` (G, C) | 1.98 s | 31% |
| **`np.linalg.solve`** | **0.133 s** | **2.1%** |
| step controller | 0.021 s | 0.3% |

Inside assembly the top costs are not arithmetic:

| | tottime | calls |
|---|---|---|
| `numpy.ufunc.at` | 0.917 s | 255 588 |
| **`toolkit.py:53` `Toolkit.__getattr__`** | **0.895 s** | **547 119** |
| builtin `getattr` | 0.893 s | 1 059 088 |
| builtin `hasattr` | 0.317 s | 360 392 |

**~34% of total runtime is attribute-lookup machinery.** `circuit.py:1341` and `:1390`
call `hasattr(self.toolkit, 'add_at')` once per element per stamp; `NumericToolkit` has
no `add_at`, so `Toolkit.__getattr__` fires, misses, and **constructs a formatted error
string and raises** (`toolkit.py:66-68`) 255 000 times. A helpful error message is being
paid for on the hot path.

### 2.1 Threaded BLAS makes the solve 19x slower

`np.linalg.solve`, n=136, same matrix, same machine:

```
OMP_NUM_THREADS=1  ->  0.238 ms
OMP_NUM_THREADS=4  ->  4.462 ms      (18.7x WORSE)
scipy lu_factor/lu_solve, 4 threads -> 0.318 ms
```

Thread-spawn overhead dwarfs 1.7 MFLOP of work. On the full transient, three repeats:
13.19/9.60/11.29 s at 4 threads against 5.48/3.96/4.45 s at 1. **2.3-2.5x for one
environment variable.** Production simulators run single-threaded sparse kernels for
exactly this reason — circuit matrices are far too small and too sparse for threaded
dense BLAS.

### 2.2 A third of all assemblies are redundant

Instrumented over 724 accepted steps: `G=3.17  C=3.17  i=3.17  u=3.17  q=5.22` per step,
where Newton needs ~2.17. `transient.py:218-219` re-assembles everything at the accepted
point immediately after Newton converged there, and `StandardNewton` discards the `(F, J)`
it already has (`nrsolver.py:33`, `:67`). `transient.py:364` and `:409` each recompute
`cir.q(x)` at that same `x` again.

### 2.3 Compound result

| configuration | wall (21 steps) |
|---|---|
| baseline, 4 BLAS threads (as shipped) | **8.80 s** |
| + fast assembly | 5.96 s |
| + no redundant `func(x)` | 4.77 s |
| baseline, 1 BLAS thread | 3.11 s |
| + fast assembly | 1.35 s |
| + no redundant `func(x)` | **0.84 s** |

**10.5x end to end, waveform drift from baseline exactly `0.00e+00`.** Roughly three
days of work, none of it touching numerics.

---

## 3. Matrix ordering and the linear solver

The maintainer flagged this as a missing topic. The finding is sharper than "missing".

### 3.1 It is not missing, it is inapplicable — assembly is dense

`circuit.py:1290` allocates `toolkit.zeros((n,n))` and elements stamp into it;
`_numeric.py:28` is `np.linalg.solve`, i.e. LAPACK `dgesv`, dense LU with partial
pivoting. There is no sparse storage, no symbolic phase, no fill-reducing ordering, no
factor reuse, and no pivot-threshold control. **There is nothing to order.**

Everything downstream is dense-only and would move with it: `analysis.py:63-69`
`remove_row_col` does an O(n^2) `np.delete` copy of J on **every Newton iteration**
purely to drop the ground row and column; `nrsolver.py:61` builds a dense `|J|` temporary
per iteration; `nrsolver.py:153` allocates a dense n x n identity per gmin evaluation.

**The irony worth recording: this repo already contains competent Markowitz ordering with
fill-in tracking, plus minimum-degree and RCM — at `ddd.py:1023`, `:1070`, `:1190-1196`,
`:2063-2086`. It is used only by the symbolic determinant-decision-diagram path. The
numeric solver that runs every transient has none of it.**

### 3.2 What ordering would buy, measured

On the real leapfrog Jacobian (136x136, nnz 567, 3.07% dense):

```
permc_spec=NATURAL          nnz(L+U) =  3145   fill vs A =  5.55
permc_spec=COLAMD           nnz(L+U) =  1590   fill vs A =  2.80
permc_spec=MMD_AT_PLUS_A    nnz(L+U) =  1984   fill vs A =  3.50
dense LU (the default path) nnz(L+U) = 18496   fill vs A = 32.62
```

COLAMD halves natural-ordering fill; the dense path stores 11.6x more than COLAMD. **No
accuracy cost:** relative residual 1.840e-14 for SuperLU/COLAMD against 1.153e-14 dense,
because SuperLU applies threshold partial pivoting on top of the ordering.

Scaling with circuit size (stages of the same fixture):

| stages | n | nnz | dense solve | spsolve | fixed-pattern `splu` | `lu.solve` only | LU fill |
|---|---|---|---|---|---|---|---|
| 1 | 29 | 155 | 0.035 ms | 0.241 ms | 0.114 ms | **0.007 ms** | 2.1x |
| 5 | 137 | 779 | 0.238 ms* | 1.04 ms | 0.437 ms | **0.026 ms** | 4.0x |
| 12 | 326 | 1864 | 15.4 ms | 2.88 ms | 1.34 ms | **0.054 ms** | 7.1x |
| 20 | 542 | 3104 | 39.8 ms | 5.89 ms | 2.73 ms | **0.103 ms** | 11.0x |

\* single-threaded; 2.66 ms with 4 BLAS threads, see 2.1.

**The reason to do the work is the wall, not the 2%:** dense MNA is 200 MB at n=5000 and
3.2 GB at n=20 000, so pycircuit cannot run an ordinary mixed-signal block at all. Note
also that fill grows with size (4x at n=137, 11x at n=542), so ordering matters *more* as
circuits grow — "just use scipy's default" is the right start but not the end state.

### 3.3 A third of linear solves are already redundant

63 calls to `linearsolver` for 21 accepted steps = 42 Newton + 1 DC + **21 from
`stepcontroller.py:60`**, which solves `J^-1 Eg` using the same J the last Newton
iteration just factored. With a factorisation object that survives the iteration this
is a triangular solve: 0.026 ms instead of 0.238 ms.

### 3.4 Recommended path

**Do not start with the solver. Start with the reference node** — `remove_row_col` is
what forces a dense copy and destroys any cached pattern. Every production simulator
assigns ground index -1 at elaboration and never puts it in the matrix.

- **Stage 0 (½ day, no risk).** Build the reduced index map once in
  `circuit.update_node_map()`; stamp straight into the (n-1) system; delete
  `remove_row_col` from the hot path.
- **Stage 1 (1-2 days, scipy only, no new dependency).** A `LinearSolver` strategy object
  in the shape the codebase already uses for `Integrator`/`StepController`/`Scaler`/
  `NonLinearSolver`, with `analyze`/`factor`/`solve` split so the factors survive the
  iteration. `DenseLU` (scipy `lu_factor`/`lu_solve`) below ~200 unknowns, `SuperLU`
  above. Unblocks n ~ 5000.
- **Stage 2 (2-3 days, optional dependency).** scipy's `splu` recomputes COLAMD every
  call — **94% of its cost** (0.437 ms factor+solve against 0.026 ms solve-only). Reusing
  the symbolic phase needs KLU's `klu_analyze` / `klu_factor` / **`klu_refactor`** split.
  That is what Xyce does and it is the single biggest transient lever in production tools.
  Make it optional, discovered the way `_ginac_ext` and `symengine` already are.
- **Rejected: a native Markowitz in pure Python.** The ordering machinery exists in
  `ddd.py`, but a Python LU over 3000 nonzeros will lose to compiled SuperLU.
  *Reconsider if* a Cython/numba build step becomes acceptable anyway for the assembly
  work in section 2.

**Sequencing against section 2:** at n=137 the solve is 2% and assembly is 96%. After the
section-2 fixes the solve becomes ~20% — *that* is when replacing it starts paying.

### 3.5 Dead and broken sparse code

- **`pybsmatrix.py` (340 lines) is unreferenced** — `grep` finds only the file itself.
  It is profile/skyline storage with **no pivoting** (`min_pivot` is set by
  `set_min_pivot` and never read), and `fbsub` (`:248-260`) has its back-substitution
  loop ascending with the division inside the inner loop, so it cannot produce a correct
  solve. Delete it or move it to an attic; leaving it implies a sparse capability that
  does not exist.
- **`_sparse_numeric.py` is slower than dense and its tests do not test it.** `:19-43`
  calls `csc_matrix(A)` on a **dense** array then `spsolve`, paying dense assembly plus
  conversion plus a full analyze+factor per call: 1.04 ms against 0.238 ms dense.
  `SparseNumericToolkit.build_sparse` returns a `coo_matrix` that nothing downstream
  accepts (`np.delete` and `toolkit.array(J, dtype=float)` both fail on it).
  `test_sparse_toolkit.py` passes because it builds the circuit with the **default**
  toolkit and passes `sparse_numeric` only to the analysis, so assembly stays dense and
  the sparse path is never exercised. **Fix the test first, by constructing the circuit
  with the toolkit under test.**

---

## 4. Numerical defects

### 4.1 The Gear2 repair is confirmed correct

Independent re-derivation and measurement of est/true at equal steps as h -> 0:

| estimator | ratio | target constant |
|---|---|---|
| Euler classic | 0.999333 | 1/2 |
| Trapezoidal classic | 1.000305 | 1/12 |
| Trapezoidal ywr | 1.000305 | 1/12 |
| **Gear2 classic (repaired)** | **1.000282** | **2/9** |
| Gear2 ywr | 0.750212 | 3/4 of 2/9 |

Four independently derived constants hit. **Also: the comment at `integrator.py:145-146`
claiming the two trapezoidal formulas "agree to the same 5/6 of the one-step LTE" is
wrong** — both are asymptotically exact.

### 4.2 CRITICAL — the Trapezoidal LTE estimator is contaminated by a `(-1)^n` mode

`integrator.py:108-129`, mechanism at `:103`. Trapezoidal is the only method whose
companion current depends on its own past value: `iq_n = 2*dq/h - iq_{n-1}`, a recursion
with eigenvalue **exactly -1**, so its companion error carries an undamped alternating
component that never decays. Both LTE branches take a **second difference of g**, which
amplifies that mode by exactly 4 and never cancels it.

With a clean seed this is a +-33% oscillation. But `check_order_drop` drops to Euler on
the first step **and at every breakpoint**, and the Euler companion seeds the mode at
O(h) instead of O(h^2), making the relative estimator error **O(1/h)** — worse as the
step shrinks:

```
h=1e-2   est/true:  -10.32   +12.28   -10.31   +12.41   -10.54
h=1e-3             -277.19  +279.90  -278.61  +281.36  -280.11
h=1e-4            -2975.84 +2978.68 -2977.53 +2980.38 -2979.23
```

**Operationally this couples error control to the escape hatch:** reject -> smaller h ->
larger err -> reject -> ... -> `MAX_REJECT` -> force-accept.

**What seeds it in practice is `Sin.next_event` returning quarter-period points
(`func.py:29-35`).** A sinusoid is C-infinity; those are not discontinuities. Treating
them as breakpoints re-arms the order drop four times per period on a smooth waveform.
Measured elsewhere: **120 of 1236 Newton evaluations forced to order 1** on a plain
`VSin`-driven RC — 10% of the run silently first-order, at drive-synchronous positions,
so the error is coherent rather than random.

Three tiers of correction: (a) stop treating smooth-source convenience points as
breakpoints; (b) difference a mode-free quantity — `d_n = (q_n - q_{n-1})/h` is
`(g_n + g_{n-1})/2` and carries no alternating component, needing one more charge in
history and an `h_last2` interface change; (c) replace TR with TR-BDF2, which is L-stable
and whose embedded estimator is genuinely `q'''`-based.

**This bears directly on the IM3 harness**, which sets `TrapezoidalIntegrator()` and
drives two `VSin` tones. The integrator was chosen on a step-count comparison that a
contaminated estimator could have produced; it must be re-verified before any T4 number
is quoted.

### 4.3 `PIController`'s step-size loop is linearly unstable

`stepcontroller.py:100`, `:139`, `:155`. The Gustafsson / Hairer-Wanner gains are
`k_I = 0.3/k` and `k_P = 0.4/k` where k is the order of the error estimate — the `p` that
`compute_lte` returns. The code uses the **numerators undivided**, and the `exponent =
1.0/p` computed at line 139 is **dead code on the accept path**.

```
gains            k   roots              rho
code (0.3,0.4)   2   -1.1165, +0.7165   1.1165  UNSTABLE
code (0.3,0.4)   3   -1.7758, +0.6758   1.7758  UNSTABLE
standard (/k)    2,3 +0.8000, -0.5000   0.8000  stable
```

Simulated with the code's own `min(2.0, max(0.2, .))` clamp, k=3 settles into a
**period-2 limit cycle alternating h = 0.857 / 0.429** — a permanent 2x step oscillation
running against the growth clamp every other step. With standard gains it converges to
h* = 1.000. The only test asserts `len(steps) > 10`, so this is invisible.

Also `self.last_err` is not updated on the rejection path, so after a rejection the PI
term differences against a two-steps-stale error.

### 4.4 `MAX_REJECT` force-accept violates BDF-2's zero-stability bound

`transient.py:380-393`, specifically `next_dt = min(max_step, dt * 10.0)`. Variable-step
BDF-2 is zero-stable only for step ratio `h_n/h_{n-1} < 1 + sqrt(2) = 2.414214`
(Grigorieff). Measured spectral radius of the homogeneous recursion:

```
r        rho         growth over 20 steps
2.4142   1.000000    1            &lt;- the bound, recovered exactly
2.5      1.041667    2.262
3.0      1.285714    152.4
10.0     4.761905    3.59e+13
```

And it fires: on a benign RC + `VSin` at reltol 1e-5, **`Gear2(ywr)` — the shipped
default — is the only configuration that reaches the force-accept path**, 4 times, each
taking a 10x ratio with the parasitic root amplifying 4.76 on that step. Largest accepted
error was 2.55x tolerance. Nothing is warned.

A bounded rejection count is defensible; **growing 10x in response to an error that was
too large is the wrong sign.** Dropping order is what a stalled high-order estimate is
actually asking for.

### 4.5 `Trapezoidal(lte_formula='ywr')` is order-inconsistent on non-uniform grids

`integrator.py:123`. The plain second difference annihilates `g'` only on a uniform grid;
with `h1 != h2` its leading term is `(h1-h2)*q''`, i.e. **O(h) where the truncation error
is O(h^2)**, and its sign flips with the direction of the ratio change:

```
r=h_last/h_curr    0.25     0.5     0.8     1.0    1.25     2.0     4.0
TR ywr           112.09   74.97   30.55    1.00  -35.86 -145.94 -435.79
TR classic         0.75    0.84    0.94    1.00    1.09    1.34    2.01
```

Confirmed O(1/h): the r=0.5 entry is 74.97 at h=1e-3 and 7497 at h=1e-5. The correct
divided-difference generalisation is **algebraically identical to the `classic` branch
five lines below**, so the right fix is to delete the branch and keep `'ywr'` as an alias.
Not the default for Trapezoidal, so this is opt-in — but it is selectable and silent.

### 4.6 Trapezoidal ringing on stiff modes is unmitigated

No TR-BDF2, no damping, no step-ratio guard. In-solver at `h*lambda = -1000`:

```
Backward Euler   |e_n/e_(n-1)| ~ 0.0010
Trapezoidal      |e_n/e_(n-1)| ~ 0.9960     rings at Nyquist, 0.4% loss per step
Gear2            |e_n/e_(n-1)| ~ 0.0972
```

Exactly the predicted `|R_TR(-1000)| = 0.996008`. Textbook and expected of pure TR; the
finding is the absence of mitigation, **and its interaction with 4.2** — the estimator
amplifies this ringing by 4 and then rejects steps because of it.

### 4.7 Backward Euler's LTE carries a `(1+r)/2` bias on variable steps

`integrator.py:83`. `gn - gn_1` approximates `((h1+h2)/2) q''`, not `h1 q''`, because the
`-(h/2)q''` companion offsets at `t_n` and `t_{n-1}` carry different h. Measured against
the prediction:

```
r          0.25     0.5     0.8     1.0    1.25     2.0     4.0
measured 0.6219  0.7458  0.8944  0.9933  1.1169  1.4866  2.4665
(1+r)/2  0.6250  0.7500  0.9000  1.0000  1.1250  1.5000  2.5000
```

Over the reachable range that is a 25% under-estimate after a growth step and up to 3x
over-estimate after a rejection. One-line rescale by `2*h_curr/(h_curr+h_last)`.

### 4.8 Index-2 structures break the exponent, not the estimator

Measured h-exponent of `lte = J^-1 Eg` for an `Eg ~ h^2` residual:

| circuit | slope | cond(J) trend |
|---|---|---|
| RC, RLC (index 1) | **3.00, 2.99** | ~1/h |
| VS \|\| C, IS-L cutset, VS+C-loop (index 2) | **2.00** | ~1/h, one at **~1/h^2** |

**Index-1: the `J^-1` mapping is exactly right** — this is the whole point of mapping the
residual through `J^-1` rather than dividing by `alpha0`. **Index-2: `lte ~ h^p`, one
order lower** (classical index-2 order reduction, correctly *reported*), so the
controller's `1/(order+1)` is the wrong exponent — still contracting, but geometric
rather than deadbeat, and the user silently gets 1st-order local error on those rows. No
index detection exists anywhere.

### 4.9 The Newton KCL residual test becomes vacuous as h shrinks

`nrsolver.py:61,64`. `I_scale = |J| . |x| + |F|` includes `Geq = alpha0*C`, which grows
like 1/h, so the *tolerance* grows like 1/h while physical branch currents do not:

```
h        max geq [S]   conv_f tolerance    (real currents ~1 mA, reltol=1e-4)
1e-05      1.500e-02        1.0e-04
1e-09      1.500e+02        1.5e-02        already 15x the real current
1e-13      1.500e+06        1.5e+02
```

Below h ~ 1.5e-8 s on that circuit `conv_f` is satisfied by any x and Newton converges on
`conv_x` alone — precisely the regime (small steps, stiff transients, nonlinear devices)
where the residual check is the one you want. SPICE sums element branch currents at each
node; the companion conductance is an algorithmic artefact and should not enter the
tolerance.

### 4.10 Smaller numerical items

- **`DampedNewton` can declare convergence on a step 32x larger than the one it tested**
  (`nrsolver.py:97-116`): if backtracking exhausts, `x_next` is the full step but `alpha`
  is left at 0.03125 and the `conv_x` test uses `alpha*xdiff`.
- **`fixed_timestep=True` does not fix the timestep** (`transient.py:340-347` vs
  `:415-416`): `dt` is mutated in place for breakpoint truncation and only restored when
  *not* fixed-step. Measured: expected ~20 steps, got **19 002**, with dt collapsing to
  3.276e-22 s and staying there.
- **`minstep = 1e-18` s is not a defensible floor** — `cond(J) ~ 1/h` (index 1) or
  `1/h^2` (index 2), so at 1e-18 s the solve returns noise and nothing checks.
- **`Gear2.check_order_drop` guards the safe direction** (`integrator.py:158`): it fires
  on *shrink*, which is unconditionally zero-stable, and has no guard on *growth*, which
  is the unstable one — exactly how 4.4 slipped through.
- **`Gear2`'s default `'ywr'` measures worse than the repaired `'classic'`** — bias
  swings over a factor of 16 across step ratio against classic's 4, and end-to-end it
  needed 57 rejections and 4 force-accepts where classic needed 29 and 0. The evidence
  now runs against the comment justifying the default.
- **`TRTOL = 7.0`** is SPICE3's constant, but SPICE3 applies it to a *charge* tolerance
  with `CHGTOL`; here it multiplies a *solution-vector* tolerance. Hard-coded in four
  places and not a Parameter.
- **`analysis.fsolve` returns a non-converged answer silently** (`analysis.py:192-199`);
  both `shooting.py` call sites use the non-`full_output` form.

---

## 5. Circuit-simulation defects

### 5.1 No junction limiting outside `Diode`

`BJT`, `JFET` and `ZenerDiode` define no `limit()`; only `elements.Diode` implements
`pnjlim`. A textbook common-emitter stage with a voltage-driven base fails DC at every
base resistance. Traced:

```
it0  |F|=5      |dx|=5      x=[5.00e+00 8.00e-01  5.00e+00 ...]
it1  |F|=0.292  |dx|=77.2   x=[5.00e+00 8.00e-01 -7.2232e+01 ...]
it2  |F|=inf    |dx|=nan    -> NON-FINITE
```

The BJT's conductances at x=0 are `IS/VT ~ 4e-13`, so the collector is effectively
floating; iteration 1 drives it to -72 V, `exp(73/0.0258)` overflows, and
`NumericToolkit.jacobian`'s central difference turns `inf - inf` into **nan**. Gmin and
source stepping both fail because neither prevents the first overshoot.

Needs three things: per-device `junctions` + a shared `limit()`; an `_expl()` exponential
clamp (the standard SPICE `EXPMAX` treatment) so overflow-to-nan cannot happen; **and a
fix to `SubCircuit.limit`** (`circuit.py:1194-1200`), which does `subx = x[nodemap]` —
fancy indexing, therefore a **copy** — and discards the return value. `Diode` does not
notice because it keeps its limited voltage in `_vlim` and linearises around it, which is
exactly why the copy bug survived.

### 5.2 The transient's `epar` never reaches the devices

`transient.py:211-215` calls `cir.C(x)`, `cir.q(x)`, `cir.i(x)`, `cir.G(x)` with **no
`epar`**, so every device sees the module global. `dcanalysis.py:77` passes it correctly.
Instrumented during `Transient(c, epar=<T=400>, bypass=True, bypasstol=1e-2)`:

```
Transient.epar.T=400.0  bypasstol=0.01
  inner-DC calls      : T=300,  bypasstol=-1.0
  transient-step calls: T=300,  bypasstol=MISSING
```

So devices run at **300 K whatever the user set**, and `bypasstol` is missing from
`defaultepar` entirely, so `Diode` takes its `except AttributeError` branch and uses a
hardcoded 1e-12. **`bypass=False` does not disable bypass and `bypass=True` does not
enable it** — the scheme is not off, it is disconnected. `test_bypass.py` passes because
both settings do the same thing.

The inner `DC(self.cir)` also inherits *nothing* — not `epar`, tolerances, `maxiter`,
`nrsolver`, `scaler`, `toolkit`, nor the `refnode` passed to `solve()`.

### 5.3 Source models

- **`VPulse`/`IPulse` crash on their own class defaults.** `per=0` (SPICE's "one shot")
  reaches `ceil(t/self.per)*self.per` -> `ZeroDivisionError`. `tr=0`/`tf=0` (SPICE
  substitutes TSTEP) divide by zero because `Pulse.f` evaluates **all** branches eagerly
  before `where()` selects. The exception is then caught by a blanket `except`, reported
  as a convergence problem, swallowed by 1.1, and the run dies at `minstep` — four
  masking layers between the bug and the message.
- **`VSin`/`ISin` ignore `td` and produce a *growing* waveform before the source starts.**
  For `t < td` with damping the exponent is positive: measured **2835 V from a 1 V
  source** at `t = td/5`, and 4.71 V into a real RC. SPICE holds SIN at the offset until
  TD.
- **`Sin.next_event` plants a breakpoint every quarter period** on a C-infinity waveform
  — see 4.2 for the numerical consequence.

### 5.4 `TLine`

- **The DC answer depends on whether a transient ran first.** `elements.py:1685` switches
  the stamp on `len(self.history) == 0` rather than on the analysis: `v(b) = 1` before any
  transient, `v(b) = 0` after. History is never reset between runs, so a second `solve()`
  interpolates against the previous run's timeline.
- **No breakpoints and no step cap at TD.** Measured arrival time against a true 1 ns
  delay: 1.02 ns at dt=2e-11, 1.50 ns at dt=5e-10, **4.00 ns at dt=2e-9, 10.0 ns at
  dt=5e-9**. `_interpolate_history` clamps to the newest sample when the step outruns the
  buffer, so the delay is not quantised, it is wrong and silent. Every SPICE caps tmax at
  TD/2 per line.

### 5.5 Device model gaps

- **No charge model on any semiconductor**: `BJT.C(x) == 0`, `JFET.C(x) == 0`,
  `Diode.C(x) == 0`. No `Cje`/`Cjc`/`TF`/`Cj0`. Bipolar and diode circuits have infinite
  bandwidth in transient unless the user adds explicit `C` elements; reverse recovery and
  storage time are not reproducible.
- **`Varactor.eval_q_pure`** clamps with `minimum(v, 0.99*VJ)`, so C -> 0 above the knee
  instead of SPICE's linear extrapolation above `FC*VJ` — a forward-biased varactor loses
  its capacitance instead of gaining it.
- **`MOS_ACM.__init__` calls `super(MOS, self).__init__`** (`mos.py:104`) and `MOS` is
  not in its MRO, so the class cannot be instantiated. **There is no large-signal MOSFET
  anywhere in the repo** — no CMOS transient is expressible.

### 5.6 Convergence-aid ladder — inventory

| aid | state |
|---|---|
| gmin stepping | present, DC only; injects gmin on **branch rows too**, which perturbs KVL; one failed rung aborts the ladder |
| source stepping | present, DC only; ladder `[0, 1e-2, 1e-1, 1.0]` — a **10x jump** from 0.1 to 1.0 on an exponential device, non-adaptive |
| junction limiting | **`Diode` only** — see 5.1 |
| damped Newton | `DampedNewton` exists but is **in no default chain** and is not composed with gmin/source stepping |
| pseudo-transient / dynamic ramping | **absent** |
| `.nodeset` / `.ic` / warm start | **absent**; `uic=True` means "start from zeros", not SPICE's per-element IC |
| device bypass | present but **disconnected** — see 5.2 |
| transient step recovery | present (`dt *= 0.25` on non-convergence); no gmin fallback at a hard time point, no diagnostic |
| Newton predictor | **absent** — `X[-1]` is the initial guess, no extrapolation |

---

## 6. Diagnostics, options, interoperability

### 6.1 A structurally singular matrix reports the wrong thing

A floating node (capacitor with no DC path to ground) produces:

```
NoConvergenceError: Source Stepping failed at lambda=0.0
```

Spectre says `Node 'float' has no DC path to ground`. The diagnosis is destroyed by three
layers of re-wrapping: `nrsolver.py:42` turns the `LinAlgError` into a
`NoConvergenceError` that no longer names the row, `:161-162` re-wraps as "Gmin Stepping
failed", `:197-198` re-wraps that as "Source Stepping failed". **The singular case should
be classified before continuation is attempted** — an LU zero pivot names the row, and
the row maps to a node name via `cir.get_node_name`.

Nothing anywhere reports which node failed to converge, though `xdiff`, `F`, `conv_x` and
`conv_f` are all in scope one `argmax` away, and `Circuit.name_state_vector` already
exists.

### 6.2 No run statistics at all

The result object exposes only `v`, `i`, `keys`, `values`, `table`, `build_waveform`.
There is no iteration count (`solve_system` returns it and `transient.py:158` discards it
with `x_res, _ =`), no rejected-step count, no minimum step reached, no breakpoint count,
no load-vs-solve split, no convergence log. A user whose long run produced a wrong answer
has nothing to look at — and this is also how 1.2 and 4.4 stayed invisible.

### 6.3 `.options` coverage

Present: `reltol`, `iabstol`, `vabstol`, `maxiter`, `minbreak`, `minstep`, `bypass`,
`bypasstol`, `uic`, `integrator`. Absent or hard-coded: `chgtol` (absent entirely),
`trtol` (hard-coded 7.0 in four places), `gmin` (hard-coded ladder), `itl1..itl6` (one
`maxiter` for both the DC OP and per-timestep budgets), `tmax`/`tstart` (conflated with
`timestep`), `pivtol`/`pivrel`, `temp`/`tnom`, and **`MAX_REJECT`**, which silently
accepts an out-of-tolerance step.

Also: `Transient.parameters` contains `analysis`, `epar`, `nrsolver`, `scaler` **three
times each** from `self.parameters = super().parameters + self.parameters` applied at two
inheritance levels — cosmetic now, but it makes any future `.options` introspection
unusable.

### 6.4 Analyses: present, thin, absent

**Solid:** DC (with gmin and source stepping — more than SPICE3 ships), Transient, AC,
Noise (adjoint), TwoPort, Transimpedance, **FeedbackLoopAnalysis** (loop gain/stability —
genuinely better than SPICE3, which has nothing), SymbolicDC, the DDD family, Volterra.

**Present but thin — `shooting.py` (PSS/PAC) does not scale past toy circuits:** fixed
timestep only; its `method` Parameter is **never read** (backward Euler is hardcoded and
the `Integrator` abstraction is not used); no limiter, so any real junction overflows;
`np.linalg.inv` of a dense n x n **per timestep** inside the shooting Newton; and PAC
allocates a dense `(N*M) x (N*M)` complex matrix — **at N=137, M=1000 that is 280 GiB for
`L`, plus 140 GiB for `B`, 420 GiB in all.**

> **CORRECTED 2026-07-30 (stage 0.1c).** This sentence originally read "**~150 TB**",
> which is wrong by about 365x. `(137*1000)^2 = 1.8769e10` elements; at complex128 that is
> 3.0030e11 bytes = 280 GiB. The `150` is `B`'s 1.5015e11 bytes read as 150 **GB** and
> written as 150 **TB** — a unit slip on the smaller of the two arrays. (`B` really is the
> smaller: `shooting.py:195` is `tk.zeros(L.shape)` with **no dtype**, so it is float64
> while `L` is complex.) **The conclusion is unaffected** — 420 GiB is as unallocatable as
> 150 TB — and that is precisely why the error survived: nothing downstream depended on
> the magnitude. Recorded rather than silently fixed. The mechanism that let it hide is
> the one the standing rules already name: the figure was typed into prose instead of
> being generated, and `benchmarks/transient_review/` carried no script for it. It now
> does — `benchmarks/transient_review/stage0_1c_pac_memory.py` reads both dtypes off
> `shooting.py` and prints the table. Regenerate rather than retype if this is ever edited.

**Absent:** DC sweep (`.dc`) — the most conspicuous gap, no `DCSweep` class exists;
`.tran` output/strobe interval (`timestep` is *both* the initial step and `max_step`, and
output points are the solver's own non-uniform steps, so every FFT needs hand
resampling — which `nonlinear_leapfrog_sweep.py` does, and interpolating non-uniform data
before an FFT is a real accuracy hazard); `tstart`; `.ic`/`.nodeset`; transfer function
(`.tf`); transient noise; Monte Carlo/mismatch/corners; numeric sensitivity; harmonic
balance; envelope; reachable `.four`/FFT (a good `fourier_analysis` exists in
`utilities/fourier.py` but is not reachable from a transient result); parameter sweep of
an analysis.

**And `Transient` is not exported from `pycircuit/circuit/__init__.py`** — `from
pycircuit.circuit import *` gives DC and AC but not transient. One line.

### 6.5 Interoperability

No netlist import of any kind. `circuit.py:857` `netlist()` emits pycircuit class names
and Python parameter names, so it round-trips to nothing. No waveform output at all — no
raw, PSF, CSV or Touchstone (`post/cds/psf.py` is a *reader*). The highest-value item by
a wide margin is a **SPICE-subset netlist reader**; everything else is downstream of
being able to get a circuit in.

---

## 7. Prioritised action list

Ranked by value per unit effort. Nothing here has been applied.

| # | item | effort | evidence |
|---|---|---|---|
| **1** | Stop swallowing a failed DC OP (1.1) | **½ day** | CONFIRMED, demonstrated |
| **2** | Single-thread the BLAS around the solve (2.1) | **½ day** | CONFIRMED, **2.4x** |
| **3** | Assembly: hoist `hasattr`, batch the scatter (2) | 1-2 days | CONFIRMED, prototyped **2.0x** |
| **4** | Reuse Newton's final `(F, J)` (2.2) | 1 day | CONFIRMED. **1-4 together: 10.5x, waveform identical** |
| **5** | Ramp the first step (1.2) | 1 day | CONFIRMED by 2 lenses, fix measured |
| **6** | Junction limiting for BJT/JFET + `expl` clamp + the `limit` copy bug (5.1) | 2-3 days | CONFIRMED |
| **7** | Pass `epar` in the transient; `bypasstol` on `defaultepar` (5.2) | ½ day | CONFIRMED |
| **8** | `vabstol` double role, and propagate to DC (1.3) | **decision** | CONFIRMED by 4 lenses |
| **9** | Convergence/singularity diagnostics that name the node (6.1) | 2 days | CONFIRMED |
| **10** | Run statistics on the result object (6.2) | 1-2 days | CONFIRMED |
| **11** | Trapezoidal `(-1)^n` contamination (4.2) | 2-3 days | CONFIRMED — **blocks the IM3 harness** |
| **12** | `PIController` gains (4.3) | ½ day | CONFIRMED |
| **13** | `MAX_REJECT` growth + order drop + warning (4.4) | ½ day | CONFIRMED |
| **14** | Source-model bugs: `VPulse` defaults, `VSin` td, `Sin.next_event` (5.3) | 1-2 days | CONFIRMED |
| **15** | `fixed_timestep` does not fix the timestep (4.10) | 1 hour | CONFIRMED |
| **16** | Reference node out of the matrix (3.4 stage 0) | ½ day | CONFIRMED |
| **17** | `LinearSolver` strategy, scipy only (3.4 stage 1) | 1-2 days | CONFIRMED, unblocks n~5000 |
| **18** | Delete `pybsmatrix.py`; fix `test_sparse_toolkit` to test the sparse path (3.5) | 1 day | CONFIRMED |
| **19** | `TLine` DC/transient stamp and TD step cap (5.4) | 1-2 days | CONFIRMED |
| **20** | KLU with `klu_refactor`, optional dependency (3.4 stage 2) | 2-3 days | CONFIRMED |
| **21** | Consolidate `jaxtransient` behind shared kernels (1.4) | ~1 week | CONFIRMED |
| **22** | DC sweep; `.tran` tstep/tmax/tstart; `.ic`/`.nodeset` (6.4) | ~1 week | CONFIRMED absent |
| **23** | SPICE-subset netlist reader (6.5) | ~1 week | CONFIRMED absent |
| 24 | Large-signal MOSFET; semiconductor charge models (5.5) | large | CONFIRMED absent |

**If only one thread: items 1-4.** Three days, 10.5x faster with a bit-identical
waveform, and it closes the one failure mode that returns confidently wrong answers. The
linear-solver work pays *after* that, when the solve goes from 2% of runtime to ~20%.

---

## 8. What was not reviewed

Collected from all four reports, because the gaps matter as much as the findings.

- **No reviewer ran the test suite.** Every claim rests on standalone scripts.
- `jaxtransient.py` was read structurally by the architect lens and skimmed by the others;
  its chunking/batching machinery (`solve_batched`, `run_chunk`, the TLine ring buffer
  under `vmap`) was not read by anyone.
- `shooting.py` was inventoried, not analysed; `_solve_coupled` was read for consistency
  but its Fang-2013 co-determination loop was not measured.
- `ddd.py`, `distortion*.py`, `volterra.py`, `nport*.py`, `analysis_ss.py` internals,
  `hdl.py`, `xdot.py`, `macromodels.py`, `symbolicapprox.py`, `_ginac.py`, `soe.py`:
  read only where they bear on ordering or the analysis inventory.
- Device models beyond `Diode.limit/i/G`, `BJT`, `JFET` and `Varactor`.
- `pycircuit/post/` beyond an inventory; `sim/gnucap` not read.
- The numerical correctness of `TLine`'s history interpolation.
- The cost of the proposed fixes on circuits that currently converge, and the runtime
  cost of the architectural proposals. Effort estimates are judgement, not measurement.
- `benchmark_circuits.py` and its tests were excluded (uncommitted work in flight).
- **No file outside the scratchpad was modified by any reviewer.**
