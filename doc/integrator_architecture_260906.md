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
