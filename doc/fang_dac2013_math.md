# Fang, *A New Time-Stepping Method for Circuit Simulation* (DAC 2013) — extracted math

**Provenance.** G. Peter Fang, Texas Instruments, DAC'13, May 29 – June 07 2013, Austin TX.
ACM 978-1-4503-2071-9/13/05. Local copy:
`/home/andreas/pycircuit_agy/papers/2463209.2488904.pdf`.

**How this was extracted, because it matters.** The PDF was rendered to PNG at 200 DPI and
the pages were read **as images**. Nothing here comes from `pdftotext` or any other text
extraction: superscripts and table structure do not survive it, and on an earlier paper in
this project that silently produced a fabricated exponent. Equations (1)–(18) below are
transcribed from pages 1–4; the results in §7 are from page 4 col. 2.

**Why this file exists.** Stage 12 built against a summary of this paper and got a method
that collapses the step size, while TI ships the original in a production simulator. This is
the primary source, written down once, so the implementation can be checked against the
paper rather than against a recollection of it.

---

## 1. Notation

| symbol | meaning |
|---|---|
| `v̄ ∈ R^N` | solution vector: **nodal voltages and branch currents** |
| `q̄(v̄) ∈ R^N` | node charges or fluxes |
| `ī(v̄) ∈ R^N` | resistive node currents |
| `ū(t) ∈ R^N` | independent sources |
| `m` | time-point index |
| `h_m` | step size between points `m−1` and `m` |
| `n` | order of the integration method |
| `k` | Newton iteration index |
| `ε_m` | local truncation error estimate at point `m` |
| `τ_m` | LTE tolerance at point `m` |

## 2. The DAE and its discretisation

The system solved (eq 1):

    f̄_ckt(v̄(t)) = d/dt q̄(v̄(t)) + ī(v̄(t)) + ū(t) = 0,     v̄(t₀) = v̄₀,  t > t₀

The differential term is discretised as (eq 2):

    q̇̄_m = (1/h_m) Σ_{i=0..n} α_i q̄_{m−i} + Σ_{i=1..n} β_i q̇̄_{m−i}

For the **trapezoidal rule** the paper states `n = 1, α₀ = 1, α₁ = −2, β₁ = −1`, and writes it
out separately as (eq 3):

    q̇̄_m = (q̄_m − q̄_{m−1}) / (h_m/2) − q̇̄_{m−1}

> **SOURCE INCONSISTENCY — recorded, not silently resolved.** Eq (2) with the stated
> coefficients gives `(q_m − 2q_{m−1})/h − q̇_{m−1}`, which is **not** eq (3); eq (3) is
> `2(q_m − q_{m−1})/h − q̇_{m−1}`. Consistency requires `α₀ = 2, α₁ = −2, β₁ = −1`. Eq (3) is
> the standard trapezoidal rule and is what `_lte_kernels.trapezoidal_companion` implements,
> so **eq (3) is authoritative here and the coefficient list in eq (2) is taken as a typo.**
> Do not "fix" our companion formula to match eq (2)'s literal coefficients.

For Gear / BDF methods `β_i = 0, i = 1..n` (eq 4):

    q̇̄_m = Σ_{i=0..n} α_i(h_m) q̄_{m−i}

Note `α_i` is written as an explicit **function of `h_m`** here — this is what makes
`p = ∂f̄_ckt/∂h_m` analytic. Discretised, eq (1) becomes a nonlinear algebraic system (eq 5):

    f̄_ckt(v̄_m) = 0

## 3. THE LOCAL TRUNCATION ERROR — read this before implementing anything

**Equation (6), verbatim:**

    ε_m = | v_i(t_m) − v_{i,extrapolated} |

> *"a typical industrial circuit simulator calculates the difference between the computed
> solution and the polynomial extrapolation from previous time steps, and considers the
> maximum value a good estimate for the local truncation error"*
>
> *"where `i` is the index of the variable with the maximum difference. In most of cases the
> variable is a node voltage. The corresponding node is often referred to as **controlling
> LTE node**, which may vary from time point to time point."*

**`ε_m` is a quantity in SOLUTION SPACE.** It is a voltage (or branch current) difference at
one node — computed solution minus polynomial extrapolation — maximised over the unknowns.

**This is not what pycircuit computes.** `Integrator.compute_lte` /
`_lte_kernels.third_divided_difference` build a divided difference of the **charge** vector
and scale it by `h²`, then `IntegralController` pushes it through `J⁻¹` to reach solution
units. The two agree in *order* but are different estimators, and three consequences follow
directly:

1. **Conditioning.** Divided differences of charge carry repeated `1/h` factors, so as `h`
   shrinks the estimate is increasingly rounding noise in `q̄` amplified by `1/h`. Fang's
   `ε_m` is a difference of two solution values and has no such amplification. An estimator
   that *grows* as `h` shrinks cannot be used as an equation to solve for `h` — Newton runs
   away from the root. This is the most likely explanation for the stage-12B collapse.
2. **`q̄ᵀ` is trivial.** With eq (6), `q̄ᵀ = ∂f_lte/∂v̄_m` is a **signed unit vector on the
   controlling node**, `±e_i` — nothing to finite-difference. That is why §3.2 can claim
   *"p̄, q̄ᵀ, and d can be computed explicitly with negligible computational costs."* With a
   charge-divided-difference LTE this claim is false, which is a signal that the estimator
   was substituted somewhere.
3. **The controlling node moves.** `ε_m` is a max over unknowns, so `i` changes between
   iterations — hence Fig. 4's two-stage split (below). `q̄ᵀ` must be re-formed each
   iteration, but re-forming it is free.

The classical acceptance test (eq 7), which the new method replaces:

    ε_m < γ · τ_m,        γ > 1

## 4. Classical step prediction (what the new method is compared against)

Gear's 1972 selection (eq 8):

    h_{m+1} = (τ_m / ε_m)^{1/(n+1)} · h_m

Digital-filter controller [3,4] (eq 9):

    h_{m+1} = (τ_m/ε_m)^{1/(4(n+1))} · (τ_{m−1}/ε_{m−1})^{1/(4(n+1))} · (h_m/h_{m−1})^{−1/4} · h_m

**Figure 2 — the classical adaptive flow:**

    Start → Predict a step size → Solve circuit equations (with step size FIXED)
          → Is the LTE condition satisfied?
              No  → Reduce step size ──────────────► back to "Solve circuit equations"
              Yes → Accept current time point
                  → End of transient? Yes → End / No → Advance, back to "Predict"

## 5. The new method

### 5.1 The coupled system (§3.1)

The LTE condition becomes an equation (eq 10):

    f_lte(v̄_m, h_m) = ε_m(v̄_m, h_m) − τ_m = 0

`h_m` is treated as an **unknown**, and solved with the circuit equations as one system
(eq 11):

    { f̄_ckt(v̄_m, h_m) = 0
    { f_lte(v̄_m, h_m) = 0          i.e.  F̄_coupled(v̄_m, h_m) = 0,  F̄_coupled ∈ R^{N+1}

**Figure 3 — top-level flow. Note there is NO rejection branch:**

    Start → Predict a step size
          → Solve circuit equations TOGETHER WITH the LTE condition
            (with step size treated as a variable)
          → Accept current time point
          → End of transient? Yes → End / No → Advance, back to "Predict"

*"There is no backup due to LTE."* The predicted step size is only the **initial guess** for
the Newton iteration.

### 5.2 Figure 4 — the Newton flow, transcribed exactly

    Start
      ↓
    ┌►Evaluate nonlinear devices
    │   ↓
    │ Form the linear system for the circuit equations: Jacobian matrix and RHS
    │   ↓
    │ Solve the linear system and update the solution vector          ← STAGE 1: the N system
    │   ↓
    │ Estimate local truncation error
    │   ↓
    │ Is the condition for the local truncation error satisfied?
    │   ├─ Yes → Is converged or needs to be terminated?
    │   │          ├─ Yes → End
    │   │          └─ No ──────────────────────────────────────┐
    │   └─ No  → Form the combined linear (N+1) system for      │
    │            circuit equations and the LTE conditions       │      ← STAGE 2: the N+1 system
    │              ↓                                            │
    │            Solve the combined linear system and update    │
    │            the solution vector AND the stepsize           │
    │              ↓                                            │
    └──────────────┴────────────────────────────────────────────┘

Three things this flow makes explicit that a summary loses:

- **The (N+1) system is NOT formed every iteration.** Only when the LTE condition fails.
- **The first stage is a plain circuit solve at fixed `h`** — the existing Newton, untouched.
- *"If the LTE condition is not satisfied, **the solution will be discarded** and a LTE
  equation will be formed in the second stage."* (§3.1)

### 5.3 The (N+1) linear system (§3.2)

Equation (12):

    ⎛ J̄  p̄ ⎞ ⎛ Δv̄_m^{k+1} ⎞   ⎛ −f̄_ckt(v̄_m^k, h_m^k) ⎞
    ⎜        ⎟ ⎜             ⎟ = ⎜                        ⎟
    ⎝ q̄ᵀ  d ⎠ ⎝ Δh_m^{k+1} ⎠   ⎝ −f_lte(v̄_m^k, h_m^k) ⎠

with

    p̄  = ∂f̄_ckt/∂h_m  ∈ R^{N×1}
    q̄ᵀ = ∂f_lte/∂v̄_m  ∈ R^{1×N}
    d  = ∂f_lte/∂h_m  ∈ R^{1×1}
    J̄  = ∂f̄_ckt/∂v̄_m ∈ R^{N×N}     (the original circuit Jacobian)

*"p̄, q̄ᵀ, and d can be computed explicitly with negligible computational costs [6]."*

**Two ways to solve it without an (N+1) factorisation:**

(a) *Partial LU (Doolittle-ordered solvers).* In the first stage the top `N` rows of the
`(N+1)×(N+1)` matrix are factorised with the last row left undetermined; the fully factorised
`N×N` sub-matrix still identifies the controlling LTE node. If the LTE condition fails, the
last row is loaded and factorised in the second stage.

(b) *Reduction to an `N` system* — equation (13):

    ( J̄ − (1/d) p̄ q̄ᵀ ) Δv̄_m^{k+1} = −f̄_ckt + (1/d) f_lte p̄

and equation (14):

    Δh_m^{k+1} = −(1/d) ( f_lte + q̄ᵀ Δv̄_m^{k+1} )

*"Rank-one update technique [8] can be readily used to update LU factors in the second
stage."* For iterative solvers, the matrix–vector product is modified per eq (13) and the
first stage's preconditioner is reused.

> Note `1/d` appears in both (13) and (14): `d → 0` is the ill-conditioned case, and it means
> an LTE insensitive to the step size.

### 5.4 Convergence criteria (§3.3) — TWO extra criteria, both required

Beyond the ordinary circuit-equation criteria, the Newton iteration is converged only when
both hold. First the **two-sided LTE band**, equation (15):

    γ_min · τ_m  ≤  ε_m^k(v̄_m^k, h_m^k)  ≤  γ_max · τ_m

*"where γ_max is a coefficient greater than 1, and γ_min is a coefficient between 0 and 1.
The introduction of the lower bound γ_min prevents step sizes from being unnecessarily
small. By adjusting γ_max and γ_min, we can precisely control the LTE allowed at each step."*

Second, the **step-size change criterion**, equation (16):

    | Δh_m^{k+1} |  ≤  η · h_m^k

*"where η is relative tolerance for step size."*

**Both are convergence/termination tests on the Newton iteration.** Eq (15) is what makes the
iteration stop without solving `f_lte = 0` exactly — the step only has to land in the band.

### 5.5 Approximate Newton (§3.4)

*"The coupled nonlinear system sometimes is very sensitive to the change of step size,
especially during a fast transition. Sometimes step size change needs to be limited or damped
to avoid convergence problems. It is unnecessary and computationally expensive to obtain high
accuracy for the step size. A typical value for relative tolerance η is 15%."*

In the second stage, instead of solving the coupled system, use Gear's formula on the step
size from the previous iteration `h_m^k` and the LTE from the first stage `ε_m^{k+1/2}`
(eq 17):

    h_m^{k+1} = ( τ_m / ε_m^{k+1/2} )^{1/(n+1)} · h_m^k

and **correct** the solution already computed in the first stage (eq 18):

    Δv̄_m^{k+1} = Δv̄_m^{k+1/2} − J̄^{−1} p̄ ( h_m^{k+1} − h_m^k )

*"Though polynomial extrapolation or interpolation is sometimes preferred for large step size
changes."*

*"Without the need to modify the matrix solver of the simulator, the approximate Newton
method is straightforward to implement and carries very little overhead."*

**This variant needs `p̄` only** — no `q̄ᵀ`, no `d`, no bordered solve, no rank-one update.

> **There is NO order reduction or backward-Euler fallback anywhere in this paper.** The only
> remedies given for large step changes are eq (16)'s damping and the extrapolation remark
> above. pycircuit's `Gear2Integrator.check_order_drop` (drop to Euler above
> `ZERO_STABILITY_RATIO = 1+√2`) is *our* safeguard, not Fang's, and it is a stability
> measure rather than a convergence one.

## 6. What the paper does NOT specify

Recorded so these are known to be our choices, not the paper's:

- how `τ_m` is formed (absolute/relative mix, and what the relative part references);
- the polynomial order or construction of the extrapolation in eq (6);
- values for `γ_min`, `γ_max` beyond the Class-D example, or for `η` beyond "typically 15%";
- what to do when the coupled solve fails to converge at all;
- any lower bound on `h_m`.

## 7. Reported results (§4, §4.1)

Implemented in TI's in-house simulator **TISpice**, exposing two user options **`ltemin`** and
**`ltemax`** — the `γ_min`/`γ_max` of eq (15). Results generated with the **2nd-order Gear**
method, the simulator's default; the paper notes the method also works for backward Euler and
trapezoidal.

The comparison ("standard") method predicts the step from eq (8) and **recomputes the time
point when the normalised LTE exceeds 4.63**.

Class-D audio amplifier, LTE bounded between **0.7 (`ltemin`)** and **3.0 (`ltemax`)**:

- **39% fewer time points**
- **17% better total wall-clock runtime**
- 35% of time steps altered, about half of them *larger* than predicted
- LTE more evenly distributed along the non-uniform grid (Fig. 5)

*"The maximum allowable LTE for the new method, 3.0 (ltemax) is actually smaller than that
for the standard method (4.63), but the average LTE for the new method could be larger."*

*"The stability condition, output, and some behavioral models can also limit the step size
besides local truncation error. It is possible and acceptable for LTE to go below the lower
bound at certain time points."*

> **On reusing 0.7 / 3.0 as defaults.** These are quoted against a baseline whose redo
> threshold is 4.63, so "1.0" in that normalisation is not "1.0" in
> `IntegralController`'s (where `err ≤ 1` already folds in `TRTOL`). The *ratio* transfers;
> the absolute values do not, until our `τ_m` is put on the same footing as theirs.

## 8. Gap list against pycircuit, as of 2026-08-01

| Fang | pycircuit today |
|---|---|
| eq (6) `ε_m` = ‖computed − extrapolated‖∞ in **solution space** | charge divided difference × `h²`, pushed through `J⁻¹` — **different estimator** |
| eq (8) step prediction | `IntegralController` — present, and deadbeat at `safety^p` |
| eq (15) two-sided band | present since stage 12A (`lte_gamma_min`/`lte_gamma_max`), default inert |
| eq (16) `\|Δh\| ≤ η h` | present as `lte_eta`, applied to accepted steps only |
| Fig. 4 two-stage Newton | **absent** |
| eq (12)/(13)/(14) bordered solve | `SchurCoupledNewton` exists and implements eq (12) via Schur; not wired to a correct `ε_m` |
| eq (17)/(18) approximate Newton | **absent** — `_solve_coupled` re-solves from scratch instead of correcting |
| `p̄ = ∂f̄_ckt/∂h` | **absent**; needs `du/dt` on the six `TimeFunction` subclasses, and `Pulse`/`PWL` have no derivative at their breakpoints |
| `q̄ᵀ`, `d` | **absent**; both become cheap once `ε_m` is eq (6) |

**The first row is the one to fix first.** Every other gap is mechanical; that one changes
what is being solved.
