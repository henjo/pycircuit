# Integrator architecture — how to make the next one cheap (2026-09-06)

Written right after adding TR-BDF2 and Radau IIA(3) across the whole PSS stack, so
the pain points below are measured, not imagined. The question this answers: **more
integrators are coming — what should change so the next one is tableau-only instead of
a two-day tour through three files?**

This is a design proposal, not a landed change. Nothing here is built yet; the plan at
the end is incremental and test-guarded (the current suite is the regression net).

---

## 0. The one-sentence diagnosis

**TR-BDF2 and Radau are both Runge–Kutta methods, but the code has no Runge–Kutta
abstraction, so each was re-implemented from scratch and wired in by name at ~35 dispatch
sites.** Everything bespoke about them — the stage step, the coupled monodromy, the
forward/adjoint source folds, the Lyapunov injection, the cost transform — is a function
of the Butcher tableau `(A, b, c)` and would be written *once* if the tableau were a
first-class object.

## 1. What adding an integrator costs today (measured)

Adding Radau touched, by category:

**A. The `Integrator` ABC doesn't fit stage methods.** Its abstract surface —
`get_required_history`, `companion_coefficients`, `companion_dT`, `compute_lte`,
`compute_derivatives` — is the *linear-multistep companion* form
(`iq = Σ a_k q_{n-k} + b·iq_{-1}`, `geq = a·h·G`). Euler/trap/Gear implement it. TR-BDF2
and Radau `raise NotImplementedError` on **all five** and are stepped by dedicated code.
So the ABC abstracts the wrong thing for half the methods it now hosts.

**B. Bespoke per-method machinery, copy-pasted.** Each stage method carries a near-parallel
set, differing only by tableau constants:

| concern | TR-BDF2 | Radau |
|---|---|---|
| transient step | `_solve_timestep_trbdf2` | `_solve_timestep_radau` |
| adaptive driver | `_run_trbdf2_adaptive` | `_run_radau_adaptive` |
| error estimate | inline 2(3) | `_radau_error_estimate` |
| factored monodromy | `_traverse_factored_trbdf2` | `_traverse_factored_radau` |
| dense monodromy | `_traverse_trbdf2` | `_traverse_radau` |
| monodromy matvec | `_monodromy_matvec_trbdf2` | `_monodromy_matvec_radau` |
| " transpose | `_monodromy_matvec_transposed_trbdf2` | `_monodromy_matvec_transposed_radau` |
| forward forced replay | `_forced_replay_trbdf2` | `_forced_replay_radau` |
| adjoint forced replay | `_forced_replay_transposed_trbdf2` | `_forced_replay_transposed_radau` |
| sideband fold | `_sideband_forced_trbdf2` | `_sideband_forced_radau` |
| Lyapunov pieces | `_lyapunov_pieces_trbdf2` | `_lyapunov_pieces_radau` |
| Newton residual (shooting) | `func_trbdf2` / `func_autonomous_trbdf2` | `func_radau` / `func_autonomous_radau` |

That is **twelve routines per method**, and TR-BDF2's are literally the two-stage special
case of Radau's coupled versions (its "two-vector fold" is the `A⊗B` coupling with `s=2`
and a lower-triangular `A`).

**C. Name/kind dispatch scattered across three files** (~35 live branches):
`FactoredPeriod.kind ∈ {'solved_history','plain','trbdf2','radau'}` with its own switch in
`matvec`/`matvec_transposed`; `isinstance(...TRBDF2Integrator)` in `Transient.solve_timestep`
and `_solve`; method-string checks in `_integrator_for`, `_companion_reach`,
`factored_period`, the `solve` whitelist, `_want_lte`, `monodromy_twin`'s self-sufficient
set, `_forced_replay`, `_forced_replay_transposed`, `adjoint_sideband_row`,
`_lyapunov_pieces`, and the `ppv` dense-`M` branch. Each is an edit-here-too the next method
must not miss — and a missed one fails *silently self-consistently* (the class of bug that
made an un-limited Radau mixer a pure sinusoid that still "converged").

## 2. The design: make the Butcher tableau a first-class object

### 2a. A `RungeKuttaIntegrator` base

```
class RungeKuttaIntegrator(Integrator):
    A: ndarray          # s×s Butcher matrix
    b: ndarray          # s weights   (== A[-1] when stiffly accurate)
    c: ndarray          # s abscissae
    b_hat: ndarray|None # embedded weights for the error estimate (+ optional
                        #   explicit-f0 weight, radau5's dd form)
    # derived, cached:
    stages, is_stiffly_accurate, is_explicit_first_stage,
    structure ∈ {DIAGONAL, LOWER_TRIANGULAR(DIRK), SDIRK, ESDIRK, FULL},
    Ainv, eig(Ainv)     # for the cost transform
```

The self-starting one-step contract is answered generically:
`get_required_history()→1`, `companion_reach()→1`, `is_self_starting()→True`,
`is_stage_method()→True`; the five LMM abstract methods raise with one shared message.
TR-BDF2 and Radau become **just the tableau** plus their `b_hat`.

### 2b. One generic stage implementation, structure-aware

Every bespoke routine in the table above is written once against `(A, b, c)`:

- **Stage residual / block Jacobian** — `F_i = q(Y_i) − q(x_n) − h Σ_j A_ij K_j`,
  `J[i][j] = δ_ij C(Y_i) + h A_ij G(Y_j)` — is tableau arithmetic.
- **The linear solve is chosen from `structure`**, which is where the DIRK/fully-implicit
  cost difference lives and is the whole point of keeping it generic *and* fast:
  - `LOWER_TRIANGULAR` (DIRK, TR-BDF2) → solve stages **sequentially**, one `m×m` factor
    per distinct diagonal (SDIRK/ESDIRK reuse a single factor — the TR-BDF2 one-LU win,
    now free for any SDIRK).
  - `FULL` (Radau) → the coupled `3m` solve, or the `eig(A⁻¹)` **cost transform** (1 real +
    1 complex solve via `ComplexKLUSolver`) — already built, just parameterised by `eig(Ainv)`.
- **Monodromy, forward/adjoint forced folds, sideband fold, Lyapunov injection** — all read
  the same `A⊗{C,G}` structure and the same per-stage abscissae `t_n + c_k h`. The adjoint
  is the transpose of the forward *by construction* (the dual-consistency gate that caught
  the TR-BDF2 sign flip becomes a property of one code path, not a per-method check).
- **Error estimate** — `b_hat` (or the radau5 `dd` filtered form) is tableau data; the
  filtered `((γ_r/h)C + G)⁻¹` estimate and the `1/(p̂+1)` controller are generic.

### 2c. Dispatch becomes polymorphism

- `FactoredPeriod` holds the **integrator** (or a small `StageMap` capability object built
  from the tableau) instead of a `kind` string. `matvec`/`matvec_transposed`/forced/lyapunov
  ask it. The `'trbdf2'`/`'radau'` arms collapse into one `is_stage_method` arm; `'plain'`
  and `'solved_history'` stay as the LMM arms.
- The scattered `isinstance`/method-string checks become questions the integrator answers:
  `integ.is_self_starting()`, `integ.is_stage_method()`, `integ.companion_reach()`,
  `integ.needs_x0_unknown()`. `monodromy_twin`'s self-sufficient set becomes
  `integ.carries_its_own_monodromy()`.

## 3. What this buys

- **The next integrator is a tableau.** Gauss–Legendre, Lobatto IIIC, ESDIRK4(3)7, a
  higher-order Radau — each is ~20 lines of `(A, b, c, b_hat)` and inherits the step, the
  adaptive control, the full PSS small-signal stack, and the cost transform with **zero new
  stepping/monodromy/noise code**.
- **One place for each concern.** The dual-consistency and limiting invariants (both of which
  cost real bugs this cycle) are guaranteed by construction rather than re-checked per method.
- **The DIRK fast path is recovered generically** from `structure`, so unifying does not cost
  TR-BDF2 its one-LU advantage.

## 4. Migration — incremental, each step green before the next

The current suite (TR-BDF2 + Radau across transient and the PSS stack) is the regression net:
every step below must leave it passing, and each new generic routine is first gated **against
the bespoke one it replaces** (same answer to machine precision) before the bespoke one is
deleted — the §0j discipline (validate the instrument against the shipped implementation).

1. **Introduce `RungeKuttaIntegrator`** with the tableau + derived structure; make
   `TRBDF2Integrator` and `RadauIIA3Integrator` subclasses that supply only `(A,b,c,b_hat)`.
   Add the polymorphic predicates. No behaviour change yet (bespoke code still runs).
2. **Unify the transient step** into `_solve_timestep_rk` (structure-aware solve). Gate vs
   `_solve_timestep_{trbdf2,radau}` on RC + diode mixer, then delete the two.
3. **Unify the shooting monodromy** (`_traverse_factored_rk`, matvec/transpose) behind
   `is_stage_method`; gate vs the pencil and the bespoke transposes; delete the four.
4. **Unify the forced/adjoint/sideband/Lyapunov surfaces**; gate vs dual consistency and the
   gear/AC references; delete the bespoke folds. Collapse `FactoredPeriod.kind`.
5. **Unify the Newton residuals** (`func_rk`, `func_autonomous_rk`) and the adaptive driver.

Steps 1–2 are self-contained and low-risk (transient only); 3–5 touch the PSS stack and are
where the dispatch count actually falls. Each is independently shippable.

## 5. Scope guards (what NOT to fold in)

- The **LMM path stays as it is.** Euler/trap/Gear have a genuinely different structure
  (companion history, the manufactured opener, solved-history pairs); this proposal unifies
  the *stage* methods and leaves the LMM branch a peer, not a subclass.
- **No auto-selection of the cost transform.** It stays opt-in (`_radau_use_transform`) —
  simplified Newton can stall on a strongly nonlinear step, and the choice is the caller's.
  A generic RK base does not change that.
- **Fully-implicit methods need `A` invertible** for the coupled/transform paths (Radau has
  `det A = 1/60`); an explicit-first-stage method (`A` singular) uses the sequential path
  only. The base should expose this so a new tableau selects the right machinery instead of
  failing obscurely.

---

## Implementation log — steps 1–2 landed; a step-3 obstruction, measured (2026-09-06)

**Steps 1–2 are done and green** (280 passed on the transient+shooting files):

- **Step 1** — `RungeKuttaIntegrator` base carrying the Butcher tableau `(A,B,C)` +
  `stage_structure()` (ESDIRK/SDIRK/DIRK/FULL); the four polymorphic predicates on the
  `Integrator` base (`is_stage_method`, `companion_reach`, `carries_own_monodromy`,
  `needs_x0_unknown`); the pure-predicate dispatch sites (`_companion_reach`,
  `monodromy_twin`'s self-sufficient set, the `x0_unknown` forcing, `_want_lte`) rewired to
  ask the method. `_integrator_for` is the one validated method→class map.
- **Step 2** — ONE tableau-driven transient step, `_solve_timestep_rk`, structure-aware:
  DIRK/ESDIRK → `_rk_step_dirk` (stage by stage via `self._newton`); FULL → `_rk_step_coupled`
  (the dense/transform path). One adaptive driver `_run_rk_adaptive`
  (`1/(EMBEDDED_ORDER+1)` exponent). The bespoke `_solve_timestep_trbdf2` and
  `_run_trbdf2_adaptive` deleted (~210 lines). Gate: generic step bit-identical to bespoke on
  the mixer, 4e-14 on van der Pol.

**Step 3 obstruction — ROUTING DIRK THROUGH THE COUPLED MONODROMY IS WRONG ON A DAE, and
the suite proved it.** The first step-3 attempt routed *every* stage method through the
fully-implicit coupled monodromy (the Radau routines, which are tableau-generic). It failed
8 shooting tests with `Singular matrix, diagonal … is exactly zero`. The cause is structural,
not a bug: a lower-triangular tableau with an **explicit first stage** (ESDIRK — TR-BDF2 has
`A[0,0]=0`, `c0=0`) gives a coupled block `J_block[0][0] = C(Y_0) + h·A[0,0]·G = C` alone,
which is **singular on any circuit with an algebraic variable** (a voltage-source branch
current, an inductor cutset — i.e. almost every circuit). Radau (fully implicit, `det A =
1/60 ≠ 0`, no explicit stage) has no such block, which is why it works coupled and TR-BDF2
does not. The transient step never hit this because it *already* solves DIRK sequentially
(`_rk_step_dirk`) — only the shooting monodromy attempt routed DIRK through the coupled solve.
Reverted to the green step-2 state.

**Corrected step 3–5 plan (structure-aware, two families — the honest design):** the shooting
monodromy/folds/Lyapunov must be **structure-aware**, exactly like the transient step:

- **FULL family** (fully-implicit — Radau) — the coupled `sm×sm` solve. The current Radau
  shooting routines already read the tableau, so making them `_rk_full_*` and dispatching FULL
  methods to them makes **any new fully-implicit method (Gauss, Lobatto IIIC, higher Radau)
  free** — low risk, high value.
- **DIRK family** (lower-triangular — TR-BDF2) — a SEQUENTIAL per-stage monodromy that
  eliminates explicit stages (`Y_0 = x_n`, `dY_0/dx_n = I`) and forward-substitutes the
  implicit stages, so no singular block ever forms. TR-BDF2's existing two-stage shooting
  routines ARE this family at `s=2` implicit stages; generalising them to s-stage sequential
  makes a future ESDIRK free. This is the larger, more delicate piece.
- **Dispatch** collapses from method-NAME (`kind=='trbdf2'`/`'radau'`) to STRUCTURE
  (`kind∈{'dirk','full'}` or `integ.stage_structure()`): still two arms per site, but a new
  method reuses the arm its structure selects instead of adding a name — which is the actual
  "next integrator is cheap" win.

⚠ THE LESSON (0j again): the shooting stack's DIRK-vs-fully-implicit split is not incidental
duplication — it is required by the DAE. "One coupled implementation for both" is refuted by
the singular explicit-stage block. The right unification is structure-aware with two generic
families, not one.

---

## Refactor COMPLETE — steps 1–5 landed, ESDIRK4(3)6 proves it (2026-09-06)

All five steps done and green. A new integrator is now tableau-only.

- **Steps 1–2** (transient side): `RungeKuttaIntegrator` base + polymorphic predicates;
  one tableau-driven, structure-aware transient step (`_solve_timestep_rk` → DIRK-sequential
  or FULL-coupled); one adaptive driver. Bespoke TR-BDF2 step/driver deleted.
- **Step 3a** (FULL family): the Radau coupled shooting routines made tableau- and s-generic
  (`_*_full`, `kind='full'`). Any new **fully-implicit** method (Gauss, Lobatto IIIC, higher
  Radau) is now tableau-only.
- **Step 3b** (DIRK family): a new generic **s-stage sequential** shooting family (`_*_dirk`,
  `kind='dirk'`) — monodromy (factored/dense/matvec/transpose), forward/adjoint/sideband
  folds, Lyapunov. TR-BDF2 moved onto it; its bespoke shooting routines (9 methods + 2 nested
  funcs) deleted. Any new **DIRK/ESDIRK** of any stage count is now tableau-only.
- **Step 3c** (dispatch): `FactoredPeriod.kind ∈ {'plain','solved_history','full','dirk'}`;
  every surface routes by STRUCTURE (`is_fully_implicit()`), never by method name.
- **Test vehicle**: `ESDIRK43Integrator` (KenCarp4, s=6, order 4) — added as **only a Butcher
  tableau** — reaches order 4 through the untouched generic transient step AND the generic
  s-stage DIRK shooting family (monodromy ratios →16 vs `exp(μT)`, adjoint exact). That an
  s=6 method works with zero method-specific code is the proof the abstraction is right.

### The obstruction, sharpened (measured, with docs-0d)

Routing every stage method through the FULL coupled monodromy is **wrong on a DAE**, and the
reason is stronger than "the transform needs `det A ≠ 0`". Measured, coupled diagonal block,
n=4:

| | | block(1,1) rank | cond(K) |
|---|---|---|---|
| ODE (M=I) | ESDIRK | 4/4 | 1.11 |
| | Radau | 4/4 | 1.13 |
| DAE (M=C, rank 3/4) | **ESDIRK** | **3/4** | **9.25e15** (singular) |
| | Radau | 4/4 | 3.07e2 |

With `a11 = 0` the diagonal block is `M − h·a11·J = M`, singular exactly when the problem is a
DAE. So an ESDIRK's coupled system is **not formable at all** on a DAE — eliminating the
explicit stages sequentially is not an optimisation, it is the **only formulation that
exists**. The two-family split (FULL coupled ⇔ `det A ≠ 0`; DIRK sequential otherwise) is thus
the **only correct partition**, and the discriminant is exactly `det A`. `cond = 1.11` on the
ODE — DAE-specific, the same place everything else in this stack comes apart.

### PCNR wired for the DIRK stage methods (transient AND shooting)

⚠ Standalone transient is a first-class use of the RK methods, so per-step features matter
there, not only in shooting. **PCNR** (the junction-continuation limiting, `Aadithya et al.`) is
now the per-stage limiting for the DIRK/ESDIRK path (`_rk_stage_pcnr`): each implicit stage's
residual recasts as the DC-flow form `i(Y) + iq_eff + u = 0` (`iq_eff = (q − target)/(h·a_ii)`),
so `pcnr.augmented_system`/`predict`/`refine` apply unchanged, with a `limit(x,x)` at
convergence to sync the devices' `_vlim` for the downstream `K`/`J`. It is NOT a separate
shooting feature — PCNR lives in `solve_timestep`, which shooting's inner transient also drives,
so forwarding a `pcnr` Parameter from `PSS` reaches it there too (the review's PCNR-in-shooting
gap, 1-for-4, is now 2-for-4). Verified: TR-BDF2 / ESDIRK4(3)6 with `pcnr=True` match device
limiting to machine precision (2.8e-17 / 7.5e-17) in transient, and PSS matches to 7e-18 with
identical spectral radius.

### PCNR wired for the FULL coupled (Radau) path — and it is *more accurate* than device limiting

`_rk_step_coupled_pcnr` closes the last PCNR gap: PCNR is now the first-class per-step limiting
for the fully-implicit coupled path too, gated behind `_rk_use_pcnr()` in `_rk_step_coupled`
(alongside the cost-transform branch). The coupled `3m` Newton keeps its exact structure; only
the per-stage residual current and Jacobian change. For each stage `j`, `pcnr.augmented_system`
builds the junction system at the stage's limiting voltages `v_lim[j]` and `pcnr.schur_reduce`
reduces it onto the MNA size:

- **residual current** = `g_mna` (the augmented MNA residual, which STAMPS the junction current
  at `v_lim` via `dev.stamp` — the *physical* current, `== cir.i` once `v_lim` == branch);
- **Jacobian block** = `δ_ij C + h A_ij G_eff`, `G_eff = J_eff` from the Schur reduction;
- the junction is eliminated from the coupled step by the correct phase (`dx_lim_of` on the
  coupled `dx_MNA`, then `refine`) — no `cir.limit` during the solve, one `limit(Y,Y)` sync per
  stage at convergence for the downstream `i`/`G`/estimate.

⚠ **THE ONE TRAP (0j, measured):** the residual current is `g_mna`, **not** the Schur RHS
`f_eff = g_mna − J_ml g_lim`. `f_eff` folds the junction current into the Newton *step*
right-hand side, where it VANISHES as `g_lim → 0`; using it as the residual dropped the junction
current at convergence and converged to a neighbouring wrong root (step-1 node error 8.5e-5).

⚠ **THE FINDING (thesis confirmed) — and a real bug in the DEFAULT coupled path it exposed.**
Validated NOT against device limiting but against an independent limiting-free coupled Newton
(`_vlim := branch` each iterate — a third code path, shared with neither PCNR nor the
coupled-limiting step, so it can rubber-stamp neither). **Radau PCNR reaches that exact
collocation root to 1e-16 at every step** (va = 0.3 / 0.8 / 2.0). The DEFAULT coupled limiting
path did **not** — and the reason was a genuine bug, not a tolerance:

> **The three collocation stages are solved simultaneously but share ONE device `_vlim`.** The
> per-stage step-limit (`cir.limit(Y_trial, Y_prev)`) leaves `_vlim` at the LAST stage; the next
> residual assembly then reads `cir.i(Y[j])`/`cir.G(Y[j])` for stages 0/1 at *another stage's*
> junction voltage (`Diode.i` reads `_vlim`). So the coupled Newton converged **cleanly** — its
> own residual to ~1e-28, limiting never even clamping — to the root of a **wrong** residual:
> node error ~8.5e-5, true stage residual ~3e-16, NOT tightening with `reltol` (the residual
> itself is off), and a spurious 4× step-halve on the hard drive. My first hypothesis (exit on
> the *limited* step) was **wrong** — instrumentation showed raw step == limited step every
> iteration; the culprit is the shared `_vlim`. The sequential DIRK/LMM paths never hit it: each
> stage owns `_vlim` for the duration of its own `self._newton`.

**The fix** (`_rk_step_coupled` and `_rk_step_transformed`): re-sync `_vlim` to each stage
(`cir.limit(Y[j], Y[j])`, zero delta) before reading its `i/q/C/G`. Both coupled paths now reach
the exact collocation root (== PCNR == the limiting-free Newton, bit-for-bit), converge in ~2
Newton iterations instead of 5, and no longer spuriously step-halve. **Measured net win: the
shooting suite dropped from 373 s to 266 s (~29% faster) with 241 pass / 1 skip unchanged.** PCNR
needed no such fix — it carries an explicit per-stage `v_lim`, which is exactly why it was
already correct and is the structurally right way to limit a coupled multi-stage solve (the
user's thesis, now with a concrete bug to show for it). Tests:
`test_pcnr_coupled_radau_solves_the_collocation_exactly` (both paths vs the limiting-free
reference) and `test_radau_cost_transform_matches_the_dense_coupled_solve` (its harmonic check
moved from bin 2, a non-harmonic bin that only rang on the old bug's artifacts, to bin 6, the
real 2nd harmonic). The PCNR-in-shooting scoreboard is now 3-for-4 (LMM, DIRK, FULL all carry
PCNR through `solve_timestep`).

### The same shared-`_vlim` defect was in the MONODROMY — found, measured, fixed

The transient fix above named a *class* of bug, so the monodromy was checked for it rather than
assumed clean. It had it. `_C_at`/`_G_at` — which feed **every** period-map builder (LMM
`_traverse_factored*`, FULL `_traverse_factored_full`, DIRK `_traverse_factored_dirk`, plus forced
replay, sidebands and the Lyapunov/noise paths) — called `cir.C`/`cir.G` with no limiting sync,
while the monodromy evaluates them at `x_n` **and every stage**. The step's solve leaves `_vlim`
at the last stage, so the period map linearised the junction at the wrong voltage.

**Measured against a finite difference of the discrete period map** — a reference this code cannot
influence — on a diode fed through a series resistor:

| | as-is | `_vlim` synced |
|---|---|---|
| Radau, va=15 | 1.65e-03 | 3.0e-09 |
| Radau, va=5 | 7.2e-06 | 1.1e-10 |
| TR-BDF2, va=15 | 9.1e-04 | 4.2e-09 |
| TR-BDF2, va=5 | 5.1e-06 | 3.5e-10 |

⚠ **The error was told from FD noise by a δ-sweep: it stayed FLAT at 1.6533e-03 across four
decades of the FD step (1e-4 … 1e-8).** FD noise makes a V (truncation ∝ δ² down, roundoff ∝ 1/δ
up); a fixed plateau is a real error. It also scaled with how hard the junction was driven and
vanished when it was off — consistent with the mechanism and with nothing else.

⚠ **Two degenerate regimes made this check vacuous first, and both were hit.** A fast RC drove the
whole monodromy to ~0 (multiplier 1e-51: a 0-vs-0 comparison that passes regardless — the
zero-vs-zero gate again); and a bare diode either never conducted (the capacitor shunts the node,
so the multiplier is exactly the diode-off `exp(-T/RC)` and the junction is absent from the answer)
or conducted so hard it shorted the node and the multiplier collapsed to 0. **Name the number the
arithmetic predicts first.** The fix was to bound the junction conductance with a series resistor:
off → `exp(-1)` = 0.368, fully on → `exp(-2)` = 0.135, solution in between. The regression test
asserts the junction is loading, so it cannot silently decay into a vacuous check.

Fix: `_sync_limit_at`, called by `_C_at`/`_G_at`. Regression: shooting 241 pass / 1 skip unchanged
(294 s vs 266 s — ~10% for the extra `limit()` calls, the price of a correct Jacobian). Test:
`test_monodromy_matches_a_finite_difference_of_the_period_map`, **verified to fail (1.653e-03) with
the fix neutered** — a test that passes both ways would prove nothing.

### ⚠ RETRACTION — the "`trap` is off by 8.5e-3" finding does not exist

An earlier revision of this document recorded, as a remaining gap, that `trap`'s monodromy was off
by ~8.5e-3 against the FD reference on a purely linear circuit, and called it a FOURTH independent
pointer at *make `x_0` the unknown*. **That was wrong, and the error was mine: I compared against a
reference I had not shown to be the right map.**

With a MANUFACTURED opening (`x0_unknown=False`, the default for a one-step method) the shooting
unknown is `x_in`, but the monodromy is about the POST-manufacturing state. My forward map started
at `x_in` and skipped the manufacturing step, so it was a **different map** — its periodicity error
was **5.13**. I had checked periodicity for `radau`/`trbdf2` and not for `trap`, and the unchecked
case is the one that bit. Re-measured in the formulation whose map IS a function of the state,
`trap` with `x0_unknown=True` gives periodicity 4.4e-16 and **RELDIFF 4.4e-10** — correct.

The rule this cost: **assert the reference map is periodic at the solution BEFORE reading anything
into a disagreement.** That assertion is now in the test, with this episode named in its docstring.
It also means the `x_0`-as-unknown case gained nothing here; the three prior findings stand alone.

### All four period-map families now FD-verified, and `_G_at` shares PCNR's limiting

`gear` — never checked before — is the multistep **2m PAIR** map (`solved_history`), seeded by
`_install_history(x0, xm1, hs[0], h_prev=hs[-1])`, with `matvec` layout
`v = (v_0, v_{-1}) -> (P_last v, P_prev v)`. Extending the harness to perturb the pair closes it.
Measured against the FD reference, junction active:

| family | method | RELDIFF |
|---|---|---|
| `solved_history` | gear | 4.3e-09 |
| `plain` | trap (`x0_unknown=True`) | 4.4e-10 |
| `full` | radau | 3.0e-09 |
| `dirk` | trbdf2 | 2.4e-09 |

`_G_at` now builds `G` through `pcnr.augmented_system` + `schur_reduce` from an explicitly passed
`v_lim` whenever the circuit has junctions (`_pcnr_junctions()` caches the device scan), so the
transient and the monodromy share ONE limiting. **What that buys is statelessness, not accuracy:**
varying the device's prior `_vlim` before evaluating at a fixed point moves PCNR's `J_eff` by
**0.0** and the `limit(x,x)` route's `G` by up to **15.15** — the old route was right in the
traversal only BY LOCALITY. Monodromy numbers are identical either way, and the suite cost nothing
measurable (269.8 s vs 274.6 s).

### The PCNR coupled path has no ladder — deliberately, with the design recorded

A continuation rung was built for the PCNR coupled solve three times (gshunt, junction-gmin, and a
junction capacitance) and removed each time. The design is kept here so a fourth attempt starts from
the evidence instead of repeating it. The path instead **falls back** to the device-limiting coupled
solve, which does carry a ladder.

**The design, if it is ever needed.** A capacitance across the limited junction, anchored at the last
accepted state: the two-node incidence stamp `JunctionGminSteppingNewton` uses, but carrying
`g (v_j − v_j,n)` instead of `g v_j` — the backward-Euler form of `C = g·h·a_ii` in parallel with the
junction. It rides `pcnr.augmented_system`'s existing `u_extra`/`J_extra` hooks, so it needs no new
plumbing and the Schur reduction carries it by construction. Two rules that cost measurements:

- ⚠ **Across the junction, not on every row.** `g·eye(n)` also anchors a voltage source's
  BRANCH-CURRENT row, where `g(i − i_n)` is a conductance applied to a current unknown. Measured at
  equal strength, the two-node stamp held the reappearing junction gap to **50 V** where the
  whole-diagonal one let it snap back to **359 V**.
- ⚠ **The schedule must start above the circuit's own conductance** (`‖G(x_n)‖∞`). A first attempt
  marched `g ≤ 1 S` against a 10 mΩ (100 S) source and never bit.

**It works as a mechanism**: with the anchor present the strong rungs converge in two iterations with
`|g_lim| = 0`.

⚠⚠ **Why it is not built.** No circuit is known where PCNR fails at a normal iteration budget.
PCNR's one documented failure — the BJT mirror of `test_dc_pcnr.py` from a uniform 20 V start — was
fixed at its source by **limiting the seed** (`pcnr.v_lim_init`; +20 V went `LinAlgError → 8
iterations`), and that docstring already recorded that a continuation could never have fixed it:

> "No ladder around the solve could help, because every rung began by building the same Jacobian at
> the same unlimited seed."

Re-measured: that mirror solves at DC with `pcnr_status='used'`, and driven as a transient (pulsed
0→5 V, rise times to 1 ps, steps to 1 ps, radau and trbdf2, PCNR on and off) every combination
converges with 0 rungs and 0 fallbacks. Transient also suppresses the mode structurally — PCNR fails
on a far-off *initial guess*, and every transient step starts from the last accepted state.

**The trigger to watch for**, if this is revisited: a PCNR stage failure whose `g_lim` stays large
while the MNA state diverges (signature: `|g_lim|` in the hundreds of volts, `ynorm` running to
1e17) on a circuit at a DEFAULT `maxiter`. ⚠ **Do not validate a candidate by starving `maxiter`** —
a ladder must end with a pure solve of the original system, so a starved budget defeats the final
rung whatever the deformation. Two of the three rungs were rejected on evidence from exactly that
broken instrument, and the third was accepted on it before this was understood.

⚠ **And read `pcnr.py`'s own docstrings first.** They answered this question before any of it was
built.

### Remaining gaps (honest)

- ⚠ **`_C_at` keeps `_sync_limit_at` and is NOT pinned by a test.** PCNR re-stamps `i`/`G` at
  `v_lim` but leaves `q` alone (`pcnr.py` treats the algebraic equations; diffusion charge is its
  stated caveat), so the capacitance cannot join. Neutering that sync alone leaves the FD test
  passing — this diode's charge does not read `_vlim` — so it is kept as correct-in-principle
  insurance for a device whose charge does, and should be treated as unverified until one is used.
### The continuation rescue: reaching it at all, and what PCNR does when it needs one

`_solve` arms `_continuation_rescue` at `minstep` as the last resort before giving up. Two things
stood between a bad circuit and that ladder, and the second was much larger than the first.

**(a) The coupled `sm` solve had no ladder.** `_continuation_rescue` is read inside
`Transient._newton`, which the coupled solve does not go through — measured, flag armed over a
40-step run: **TR-BDF2 wrapped the rescue solver 80 times, Radau 0**. ⚠ `self._newton` could NOT be
reused: it is MNA-SIZED (reduces an `n`-vector at `irefnode`, limits a full `n`-vector, per-MNA-row
tolerances and row names), so a `3m` block system cannot be handed to it, and its device-limiting
route would reintroduce the shared-`_vlim` hazard fixed above. So the coupled path carries its own
**gshunt** ladder. ⚠ The shunt enters as a CONDUCTANCE IN THE DEVICE CURRENT (`i + g x`, `G + g I`),
not as `F + g x`: this residual is in CHARGE units with Jacobian `C + h A G`, so `g` must arrive
where `G` does. Exercised: a 6-diode 400 V slam at 5 GHz fails outright at a 3-iteration budget and
converges with **2 gshunt rescues**; an 8-iteration budget fires it 0 times.

**(b) ⚠⚠ THE LADDER WAS UNREACHABLE ON THE DEFAULT PATH ANYWAY.** Adaptive stepping routes every
Runge–Kutta method to `_run_rk_adaptive`, not through `_solve`'s loop, and that driver halved to
`minstep` and then bare-`raise`d — `_continuation_rescue` appeared **0 times in it against 3 times
in `_solve`**. So (a) had been validated through `fixed_timestep=True`, a door users do not come
through, and *no* stage method — Radau, TR-BDF2, ESDIRK, PCNR or not — could reach a continuation in
normal operation. `_run_rk_adaptive` now arms the chain at `minstep` before giving up.

**(c) PCNR reaches the ladder by FALLING BACK, not by carrying one.** A gshunt rung and a
junction-gmin rung were both built for the PCNR coupled solve and **measured not to rescue it**:
instrumenting the failure shows `max|g_lim|` — the junction limiting residual — crawling from 359 V
while the MNA state diverges to 4.6e17. PCNR's bottleneck is the junction limiter's **slew rate**,
which no deformation of the circuit accelerates; and wiring the ladder anyway turned a 0.02 s
failure into a >136 s one (60 rungs at every halving level). So both PCNR stage paths now fall back
to the device-limiting solve that does carry a ladder — the same fallback the LMM step and DC have
always done on a PCNR failure. Measured on the default adaptive path, all four combinations now
attempt a real continuation (`pcnr=False` directly; `pcnr=True` after 3 and 2 fallbacks).

⚠ **The fallback is not free, and the code says so.** The two limitings agree at the root (0.0 /
5e-18 relative on a single junction), so a fallback step is not a different answer — EXCEPT on
PARALLEL junctions on one branch, which is the case PCNR exists for and which per-device limiting
resolves order-dependently. The warning detects duplicate junction pairs and states it. Normal
operation does not fall back at all (0 fallbacks, `pcnr_status='used'`).

⚠ **A live defect this surfaced:** `pcnr_solves`/`pcnr_fallbacks`/`pcnr_status` were initialised
only in `_solve`, which **shooting never calls** — it drives `solve_timestep` directly. The stage
fallback tripped over it at once (`AttributeError`), and the LMM PCNR path carried the identical
latent bug and had simply never been reached that way. Now initialised in `__init__`, with
`_solve`'s per-analysis reset kept. Tests:
`test_a_bad_circuit_reaches_the_continuation_ladder_on_the_default_path` (verified to fail both when
the driver bare-raises and when the PCNR fallback is removed) and
`test_the_continuation_rescue_reaches_the_full_coupled_path`.

### Remaining gaps (honest)

- ⚠ **`_C_at` keeps `_sync_limit_at` and is NOT pinned by a test.** PCNR re-stamps `i`/`G` at
  `v_lim` but leaves `q` alone (`pcnr.py` treats the algebraic equations; diffusion charge is its
  stated caveat), so the capacitance cannot join. Neutering that sync alone leaves the FD test
  passing — this diode's charge does not read `_vlim` — so it is kept as correct-in-principle
  insurance for a device whose charge does, and should be treated as unverified until one is used.
### The continuation rescue now reaches the FULL coupled path

`_solve` arms `_continuation_rescue` once the step has shrunk to `minstep`, as the last resort
before it gives up. That arming used to do **nothing** on a fully-implicit method: the flag is read
inside `Transient._newton`, which the coupled `sm` solve does not go through. Measured with the flag
set over a 40-step run — **TR-BDF2 wrapped the rescue solver 80 times, Radau 0** — and `_solve` then
reported that the "gmin/gshunt/pseudo-transient continuation could not rescue the point" for a
ladder that had never run.

⚠ **`self._newton` could NOT simply be reused, and that refuted the original plan.** It is
MNA-SIZED: it reduces an `n`-vector at `irefnode`, its limiter calls `cir.limit` on a full
`n`-vector, and `abstol`/`xtol`/`row_names` are per-MNA-row — a `3m` block system cannot be handed
to it. Its device-limiting route would also reintroduce the shared-`_vlim` hazard across
simultaneous stages fixed earlier in this document.

So the coupled path carries its **own** gshunt ladder (`_stage_newton(seed, gshunt)` +
`_adaptive_conductance_ladder`). Only the gshunt rung is offered: it is the one deformation that is
structure-free — every node to ground through `g` — so it needs no MNA row map and no per-stage
junction bookkeeping. ⚠ **The shunt enters as a CONDUCTANCE IN THE DEVICE CURRENT** (`i + g x`,
`G + g I`), not as `F + g x` on the residual: this residual is in CHARGE units with Jacobian
`C + h A G`, so `g` has to arrive where `G` does or it is dimensionally wrong.

**Exercised, not merely wired.** On a 6-diode 400 V slam at 5 GHz the coupled Newton cannot converge
in 3 iterations: the step fails outright without the rescue and converges **with 2 gshunt rescues**
with it, while an 8-iteration budget fires the ladder 0 times. Test:
`test_the_continuation_rescue_reaches_the_full_coupled_path`.

⚠ **Still open:** the PCNR variant of the coupled solve has no shunt rung, so arming the flag with
`pcnr=True` on a fully-implicit method still changes nothing. `_honours_continuation_rescue`
reports that honestly and the `_solve` diagnostic follows it. And no circuit has yet been found that
defeats the coupled Newton at a NORMAL iteration budget — 16 parallel diodes behind 10 mΩ under a
500 V slam at 10 GHz, down to two points per period, all converge — so the ladder is insurance whose
natural trigger is still unobserved.

- **The FULL coupled path** uses a hand-rolled Newton (limiting only), not the full nrsolver
  (line-search/continuation-rescue) the DIRK stages get via `self._newton`.
- **The cost transform** stays opt-in (simplified Newton, falls back to dense); the DIRK
  sequential path keeps its efficient one-LU structure.
- **Continuation rescue / breakpoints** still do not reach shooting — armed only in
  `Transient.solve`, which PSS never calls (it drives `solve_timestep` on its own frozen grid).
