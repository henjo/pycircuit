# PSS / shooting — state, records and open items (2026-09-01)

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
needed; an adaptive grid of **1105** converges in a prototype at −47.3 ppm against the
uniform run's −60.6. The uniform-grid restriction costs ~17× the points.

## Open items, ranked, each with its gate

1. **Item 5's payoff case.** `solve(grid=…)` works and is tested on benign circuits;
   van der Pol still does not solve through it. Diagnosis eliminated the opening step
   size, the plumbing, and the Jacobian-on-non-uniform-grids. What is left: the plain
   path's Jacobian is ~30% off (its flat-history seed, docstring 4b) and a stiff
   relaxation oscillator is the first circuit that cannot tolerate it.
   *Gate:* does it converge with an exact outer Jacobian? `benchmarks/pss_lte_grid.py`
   is the target to hit.
2. **Item 6, matrix-free variational shooting.** Ceiling measured at m=40 as
   **1.01–1.03×** — buys nothing there. **Unmeasured above m=40**: three measurement
   designs failed on a box at load 17–27. *Gate:* on a quiet machine, does the
   propagation share pass ~30%? `benchmarks/pss_matrix_free_ceiling.py` carries the
   harness and the three failed designs.
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
| `pss_lte_grid.py` | the LTE-chosen grid's prize, and the target to hit |
| `pss_matrix_free_ceiling.py` | item 6's ceiling, and three failed measurement designs |
