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

## Open items, ranked, each with its gate

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

2. **Item 6, matrix-free variational shooting — GATE PASSED 2026-09-02, ready to build.**
   The gate was "on a quiet box, does the propagation share pass ~30%". It does:
   single-threaded BLAS, it crosses 30% near m=220 and reaches 79% at m=1002.

   | m | propagation | (solve alone) | ceiling k=20 |
   |---|---|---|---|
   | 40 | 4.8% | 2.3% | 1.02× |
   | 110 | 15.0% | 6.2% | 1.14× |
   | 242 | 38.1% | 14.1% | 1.54× |
   | 502 | 63.8% | 21.3% | 2.58× |
   | 1002 | 79.1% | 24.9% | 4.44× |

   ⚠ **The old record was overturned by a fix to the harness, not by the quiet box.**
   It read 2.2% at m=40, ceiling 1.01–1.03×, because it timed `toolkit.linearsolver`
   alone. `_step_sensitivity` also forms `C @ P` — m×m against m×2m, O(m³), the same
   order as the solve — and matrix-free eliminates those too, because it never forms
   `P`. The old number was stable and reproducible; it was answering a different
   question. Its **verdict** at m=40 survives on the corrected figure (4.8%, 1.02×).

   ⚠ Quote the threading condition with the number: threads free, the same sizes read
   4.8 / 7.2 / 12.4 / 18.2%. `k=20` is an assumption, so this is an upper bound on an
   upper bound. `benchmarks/pss_matrix_free_ceiling.py` carries it and the three
   measurement designs that failed on a loaded box.

3. **Autonomous Gear-2 diagnostics.** The composed monodromy is 2m×2m and mixes physical
   multipliers with the discretisation's parasitic roots. Harmless for Gear-2 (its
   parasitic root is 1/3 per step, ~1e-95 over a period) and **not** harmless for a
   method whose root sits nearer the unit circle. A test pins the gap.

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
| `pss_matrix_free_ceiling.py` | item 6's ceiling (gate passed, m>~250), and three failed measurement designs |
