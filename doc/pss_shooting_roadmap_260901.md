# PSS / shooting — state, records and open items

*Written 2026-09-01. On 2026-09-02 item 6 was re-measured and its verdict reversed,
and item 5's payoff case was closed — both by falsifying a recorded diagnosis.*

**Read this first if you are picking the shooting work back up.** The detailed record
lives in the `PSS` class docstring in `pycircuit/circuit/shooting.py`, which is long
and is the authority; this file is the map to it.

Everything described here is committed and pushed on `cna-jax-vectorization`.
Merging to master is the repo owner's call and has not been done.

## What PSS is now

| | |
|---|---|
| inner solve | `Transient`, driven one step at a time — one integrator definition, limiting, PCNR, breakpoints, order drops |
| methods | `euler`, `trap` (default), `gear`/`gear2` |
| formulations | plain `x_in`; solved history `(x₀, x₋₁)` for Gear-2; free period `(x₀, T)`; composed `(x₀, x₋₁, T)` |
| grid | uniform, or a caller's step **fractions** via `solve(grid=…)` |
| reports | `max_lte`, `total_lte`, `max_lte_seam`, `fundamental_period`, `spectral_radius`, `converged` |

Three convergence criteria, and PSS now says something about all three: the inner
Newton, the shooting Newton, and the **discretisation** (the LTE report).

## The numbers worth keeping

Q=20 series RLC, 100 points/period, analytic peak 20 V, `exp(-π/Q) = 0.854636`:

| method | peak | seam cost | after solved history |
|---|---|---|---|
| euler | 8.815 V | 5.1e-12 V (none) | n/a |
| gear2 | 19.766 V | **1.266e-01 V** | **19.893 V** (2.18×) |
| trap | 19.990 V | 1.3e-11 V (none) | n/a |

Quadrature phase element, true period 1.000e-03: composed autonomous Gear-2 lands on the
free-running period (+332.185 vs +332.184 ppm predicted) and the radius wobble collapses
2.095e-04 → 6.6e-12.

Van der Pol, μ=100, stiffness ratio 5443, true period 162.842412: **20 000 uniform points**
needed; an adaptive grid of **1106** converges **through the analysis** (2026-09-02) at
−73.8 ppm (trap) and −100.6 (gear), against the uniform run's −60.6. The prototype with its
own unknowns still gives the best answer, −47.3. The uniform-grid restriction cost ~18× the
points.

## The three ranked items — all closed 2026-09-02

Each had a gate written down before the work started. Two of the three gates were
answered **no**, and in both cases that is what found the real cause: item 5's exact
Jacobian did not fix it (the manufactured opening step did), and item 6's ceiling was
being measured by an instrument that timed a third of what it was named for.


1. **Item 5's payoff case — DONE 2026-09-02.** Van der Pol now solves through
   `PSS.solve(grid=…)` on its own 1105-step LTE-chosen grid: trapezoidal −73.8 ppm,
   Gear-2 −100.6 ppm (its first convergence on this grid by any driver), against the
   20 000 uniform points a uniform grid needs for −60.6. Pinned by
   `test_the_lte_chosen_grid_solves_van_der_pol_through_the_analysis`.

   ⚠ **The gate's answer was no, and that is what found it.** The gate was "does it
   converge with an exact outer Jacobian?" It does **not** — an exact finite-difference
   Jacobian on the real analysis fails identically. The blocker was the **manufactured
   opening step**: `_traverse` builds `x(0)` with one order-dropped Euler step of
   `hs[0]`, and an LTE grid opens wherever the transient's window started — here
   `h[0] = 1.4845` against a median of `4.62e-04`, 3200× coarser, and the *inner*
   Newton fails on that single step. `_period_grid` now opens a coarse-opening grid on
   its own finest step, gated at 8× so grids that already work are untouched. With that,
   the analytic and finite-difference Jacobians agree to six digits — the plain path's
   ~30% Jacobian was never what stopped this case.

2. **Item 6, matrix-free variational shooting — BUILT 2026-09-02.**
   `PSS.solve(matrix_free=True)`, driven solved-history path only. The gate ("does the
   propagation share pass ~30%?") passed — it crosses near m=220 and reaches 79% at
   m=1002 — and the thing was then built and measured end to end:

   | m | dense (iters) | matrix-free (iters) | speedup |
   |---|---|---|---|
   | 242 | 2.113 s (2) | 1.557 s (2) | 1.36× |
   | 502 | 9.255 s (2) | 6.131 s (3) | 1.51× |
   | 1002 | 52.402 s (2) | 24.636 s (3) | 2.13× |

   Answers agree with the dense path to **1.1e-16**. The traversal-only speedup is
   higher (2.23× at m=502, 3.62× at m=1002) because a `solve` also does setup, replay
   and the DFT; both are true, so say which is being quoted.

   ⚠ **`k` was an assumption and is now measured**: GMRES takes 2/4/7/12 iterations at
   m=40/110/242/502, because `I − M` has almost every eigenvalue within 1% of 1.0 — the
   fast modes decay to nothing over a period, so **k tracks the number of slow modes,
   not m**. That is the property the whole item rests on and it had never been checked.

   ⚠ **What it costs:** `2·N·m²` doubles of stored factorisations (~800 MB at m=1002,
   50 points); one extra Newton iteration at m≥502, from a convergence test that cannot
   be identical (`fsolve` scales by `|J|·|x|`, which no matrix-free method has); and no
   monodromy, so `spectral_radius` is `None`. The autonomous and plain paths are not
   converted and **raise** rather than silently going dense.

3. **Autonomous Gear-2 diagnostics — DONE 2026-09-02.** The composed monodromy is
   2m×2m and mixes physical Floquet multipliers with the discretisation's parasitic
   roots. `spectral_radius` no longer takes a maximum over that mixture:
   `_spectral_report` separates them and reports the maximum over the **physical**
   multipliers, exposing both sets as `floquet_multipliers` / `parasitic_roots` and
   warning when a parasitic root climbs within a decade of the physical spectrum.

   **The discriminator is the eigenvector's block structure, not the eigenvalue.** The
   composed map acts on the pair `(x₀, x₋₁)`: a physical mode's halves are one timestep
   apart on a smooth trajectory (`v₋₁ = e^{−λh} v₀ → v₀`), while a parasitic mode's are
   related by the method's spurious root and stay O(1) apart. That is a prediction, and
   it holds — the physical split falls **0.1281 → 0.0316** when the phase circuit goes
   from 50 to 200 points (factor 4.05 for 4× in h), and the Q=20 RLC's parasitic split
   is **1.9997** against the 2.0 that BDF-2's `v₋₁ = 3v₀` predicts exactly.

   ⚠ **The split is by rank, not by a threshold**, and that was measured into the design:
   a 0.25 threshold returned *no* physical modes on a stiff RC ladder (λh ≈ 40, so every
   mode's halves differ by O(1)) and handed back `spectral_radius = None`. A k-step
   method on m states has exactly m physical multipliers — structural — so the m smallest
   splits are the physical set and no cut has to be chosen. On a stiff circuit the labels
   can still swap, and it does not matter: every mode involved is at the noise floor, so
   the radius is unaffected.

   No method in the tree has a near-unit parasitic root, so this changes no number today.
   The case it exists for is pinned by a synthesised monodromy
   (`test_a_parasitic_root_near_the_unit_circle_is_not_read_as_stability`), because
   waiting for such a method to arrive would mean shipping the separation untested on the
   day it does.

## What is open now

Nothing from the original three. Two candidates, neither started:

- ~~**The phase condition.**~~ **TESTED AND REJECTED 2026-09-02, with numbers.** The
  gate was "a case where the frozen pin demonstrably fails". It could not be built, and
  the attempt overturned the reasoning behind the proposal.

  **The `argmax` compares volts with amperes on purpose.** The phase row removes the
  orbit's tangent in proportion to `|e_k·f̂|`; one step of `|Δx_k|` *is* `|f_k|` up to h,
  so `argmax|Δx_k|` is exactly `argmax |e_k·f̂|` — it maximises the quantity the row
  needs. The "obvious repair" of normalising each coordinate by its own swing picks a row
  **704× worse aligned** (1.4e-03 against 1.0000), measured on a van der Pol carrying a
  VCVS-scaled copy of `v`. Pinned by `test_the_phase_pin_compares_units_on_purpose`.

  **The orthogonality row buys 2–6× conditioning and never decides anything.** On the
  case built to break the pin (seeded at `v`'s turning point, pin sitting 1e-3 into its
  coordinate's range), both rows converge to the same answer at every grid, with bordered
  condition numbers 1.2e3 / 3.0e2 / 8.2e1 (pin) against 2.0e2 / 6.0e1 / 3.7e1
  (orthogonality) at 200/800/3200 points — and the gap *shrinks* as the grid refines. The
  row's alignment at the solution stayed 1.0000 throughout: the tangency the attack aimed
  at never materialised.

  ⚠ **One real failure mode survives, and it is about seeds, not the condition.**
  `phase_pin` is a *value* the orbit must attain, so a far seed can pin one outside its
  range and the system is then inconsistent, not merely hard: on van der Pol at μ=1,
  seeds at 4×/10×/30× the orbit amplitude pin `v` at −5.66/−9.50/−15.46 against an orbit
  range of [−2.01, 2.01], all reported as ordinary non-convergence. The orthogonality
  row's plane through the same far seed misses the orbit too, so it fixes nothing here.
  **The remaining gap is diagnostic** — naming an inconsistent phase plane the way the
  `T = 0` trivial root is named — and is not started.

- **Matrix-free for the autonomous and plain paths — DONE 2026-09-02.** Three of the four
  systems are now converted; only the composed autonomous `(x₀, x₋₁, T)` raises.

  | system | m=242/128 | m=502/308 | m=1002/608 |
  |---|---|---|---|
  | driven solved-history (2m cols) | 1.36× | 1.51× | 2.13× |
  | driven plain (m cols) | 1.10× | 1.34× | 1.63× |
  | autonomous plain (m + border) | 1.02× | 1.11× | 1.34× |

  Propagation share for the plain path is about half the solved-history path's — 42.5%
  vs 63.8% at m=502, 60.7% vs 79.1% at m=1002 — exactly what m columns instead of 2m
  against the same assembly must mean. All three agree with the dense path to ≤2e-16,
  and the autonomous one reproduces the period exactly.

  ⚠ **Found on the way: the plain path had its own inlined copy of `_step_sensitivity`**,
  byte-for-byte, while that method's docstring claimed the systems "share this and not a
  copy". They did not — only the solved-history path called it. Found because a timer
  wrapped around `_step_sensitivity` reported the plain path's propagation share as
  **0.0%**, which is not a small number but a wrong one. De-duplicated.

  ⚠ **And two bugs from `_coeffs` being live state.** The manufacturing step is
  order-dropped to Euler (`b = 0`) while loop steps run the method's own pair. Applying
  the opening's pair to every step put the period column 40–50% out for trap and gear;
  applying the loop's to the opening made `Pq` non-zero where `_traverse` seeds it at
  zero — 100% wrong for trap. **Both times the trajectory itself matched to zero**, and
  `euler` was exact under both, so a one-method test would have passed.


## Traps — things that looked true and were not

Every one of these was settled by measurement, and each cost real time:

- **A converged shooting solve is not evidence of a correct answer.** Three integrators
  reported success while disagreeing by 56%.
- **`k·T` solves the periodicity condition whenever `T` does**, so an autonomous solve
  follows its seed and reports success at 2× or 3× the fundamental. Now detected.
- **`T = 0` is a regular root of every autonomous shooting system.** Any seed below the
  fundamental is drawn to it. Now named.
- **A sub-block of a sensitivity is not a monodromy** — reporting one gave ρ = 1.28 for a
  decaying resonator.
- **Trapezoidal cannot be given an exact Jacobian by reformulation.** Three designs, one
  obstruction: its `iq` has an undamped `(-1)ⁿ` mode, so every coupling through `-G` is
  degenerate or rank-deficient. The plain path's Euler manufacturing step is what makes
  the problem well-posed; the ~30% Jacobian error is what that costs.
- **`doc/transient_review.md` §4.6's ringdown numbers do not transfer to PSS.** A
  periodic steady state has no transient to ring.
- **A test that could not fail, on a circuit that could not show the defect.** The
  matrix-free matvec tests were first written on a linear RC ladder, where every `C`
  along the period is the same matrix — so corrupting the recursion's capacitance ring
  changed nothing and the suite stayed green through the mutation. Seeding both history
  states with zeros hid a second defect the same way. Fixed with a state-dependent `C`
  (a `BSource` `q_func`) and two distinct entering states; both mutations are now caught.
  ⚠ Neither hole was visible from reading the test — only from mutating the code under it.
- **An observation recorded as an exoneration.** The old item-5 diagnosis listed "the
  opening step is large not small" among the obvious causes *eliminated*. It had seen the
  cause and checked it against the wrong worry — that the opening step might be too small
  to open the trajectory — never against its own size. The note that would have solved it
  was already written down, filed under things ruled out.
- **A measurement named for the question it was meant to settle, not for what it timed.**
  Item 6's harness accumulated `toolkit.linearsolver` and was called "the propagation
  share". It was stable, reproducible, and under a third of the propagation, because
  `_step_sensitivity`'s `C @ P` products are the same O(m³) and were counted as
  "assembly". It read 2.2% where the answer was 4.8%, and at m=502 the gap is 21% vs
  64% — the difference between "buys nothing" and a 2.6× ceiling.
- **The recorded blocker for item 5 was stale**, and so was the claim that a variable-step
  trapezoidal estimator had never been written.

⚠ The pattern across all of them: **a number compared against itself, or carried from the
context where it was measured into one where it had not been checked.** Every fix came
from comparing against something independent — the analytic decay, the true period from a
transient, the free-running limit cycle.

## Benchmarks

| file | what it settles |
|---|---|
| `pss_seam_cost.py` | what the cold-start seam costs, per method |
| `pss_stiff_autonomous.py` | whether a stiff oscillator needs Gear-2 (it does not) |
| `pss_lte_grid.py` | the LTE-chosen grid's prize; now a CONTROL, the analysis reaches it |
| `pss_matrix_free_ceiling.py` | item 6's ceiling and its delivered speedups; three failed measurement designs |
