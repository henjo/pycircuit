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

- **Matrix-free for the autonomous and plain paths — DONE 2026-09-02.** All four
  shooting systems now take `matrix_free=True`; nothing refuses.

  | system | | | |
  |---|---|---|---|
  | driven solved-history (2m cols) | 1.36× (242) | 1.51× (502) | 2.13× (1002) |
  | driven plain (m cols) | 1.10× (242) | 1.34× (502) | 1.63× (1002) |
  | autonomous plain (m + border) | 1.02× (128) | 1.11× (308) | 1.34× (608) |
  | autonomous composed (2m + border) | 1.22× (208) | 1.65× (408) | |

  All agree with the dense path to ≤2e-16, and both autonomous systems reproduce the
  period exactly. The plain path's propagation share is about half the solved-history
  path's — 42.5% vs 63.8% at m=502 — exactly what m columns instead of 2m against the
  same assembly must mean, and why the two 2m systems win earlier and by more.

  The composed bordered matvec is `J[v;s] = [v − Mv − s·(Pt_last, Pt_prev) ; v[k]]` —
  one phase row pinning the `x₀` block only, because time translation slides both states
  along the orbit together so the freedom stays one-dimensional.

  ⚠ **Found on the way: the plain path had its own inlined copy of `_step_sensitivity`**,
  byte-for-byte, while that method's docstring claimed the systems "share this and not a
  copy". Found because a timer wrapped around it reported the plain path's propagation
  share as **0.0%** — not a small number but a wrong one. De-duplicated.

  ⚠ **And two bugs from `_coeffs` being live state.** The manufacturing step is
  order-dropped to Euler (`b = 0`) while loop steps run the method's own pair. Applying
  the opening's pair to every step put the period column 40–50% out for trap and gear;
  applying the loop's to the opening made `Pq` non-zero where `_traverse` seeds it at
  zero — 100% wrong for trap. **Both times the trajectory itself matched to zero**, and
  `euler` was exact under both, so a one-method test would have passed.


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


## The canonical papers, read against the open questions (2026-09-02)

A second review round read Aprille & Trick (Proc IEEE 60(1) 108-114, 1972; IEEE TCT
19(4) 354-360, 1972) and Trick, Colon & Fan (IEEE TCAS 22(5) 391-396, 1975). Four
results, each verified here before being recorded.

**1. `x₀` is the canonical unknown — the manufactured opening step is a deviation, not a
refinement.** [AT-P]'s Step 1: *"for the given initial state x_0^i compute the solution
x^i(t; x_0^i), 0 ≤ t ≤ T"*. The trajectory begins **at** the unknown. So
`solve(x0_unknown=True)` is a **return**, and the design question is settled by
precedent; only the cost is open.

**2. The identity seed is theirs too — which states our frame error exactly.** [AT-P]
gives `Phi = PROD [I − h F(x)]⁻¹`, seeded by `z(0)` at t=0, i.e. the identity **at `x₀`,
which in their formulation is the unknown**. This file **kept A&T's Jacobian and changed
A&T's unknown**. `Px = [I, I]` is correct for the formulation it came from.

**3. The period column has a closed form, and it makes a free permanent check.** [AT-O]
eq.(19): `dH/dT = −f(x(T))`. That is the *continuous* derivative and is **not** what we
compute — the grid scales with T, so the exact discrete derivative carries the `dh/dT`
terms `companion_dT` accumulates, and ours is the more exact object. But the gap must be
**O(h) and must fall**:

| npts | trap | gear |
|---|---|---|
| 100 | 9.51e-02 | 3.25e-02 |
| 200 | 4.73e-02 (2.01×) | 1.61e-02 (2.02×) |
| 400 | 2.36e-02 (2.00×) | 7.99e-03 (2.01×) |
| 800 | 1.18e-02 (2.00×) | 3.99e-03 (2.01×) |

⚠ **Verified to reject the defect it was added for:** with Gear-2's period column back on
the partial `d/dh_n`, the gap **stops falling** and plateaus at exactly 0.5 — 0.486 /
0.495 / 0.498 / 0.499, ratio 1.00. A constant-factor error is invisible to a test that
asks whether a number is small and unmissable to one that asks whether it *falls*. Its
value is its **independence**: it compares against the converged waveform, so an error in
the accumulation cannot hide inside it. Now a permanent test.

**4. The phase rule is canonical; re-selecting it was TRIED AND REJECTED.** [AT-O] Step 3
selects `k = argmax|f_k(x(T))|` — vindicating our rule and our rejection of the Poincaré
upgrade. Their Step 3 also sits **inside** the loop, pinning the iterate's own value, which
looked like a free repair for the far-seed failure mode.

**Built (with a per-iteration hook in `fsolve`) and measured — it regresses the working
case.** Van der Pol at μ=1 from an **on-orbit** seed went from converged to *not*
converged, and far seeds wandered to periods of −52, −1088 and +110 against a true 6.6633.
Reverted.

⚠ **The reason is structural, not a flaw in the attempt.** Pinning the iterate's own value
makes the phase residual `x₀[k] − pin` **identically zero** at every iterate, so the row
carries no information — it constrains the *step* (`dz[k] = 0`) and nothing else, and with
`k` re-chosen each iteration a different coordinate is frozen each time, so the orbit
slides along itself. **A&T do not have this problem because they have no phase equation**:
their unknown vector *substitutes* the period for the pinned coordinate
(`v = [x_01, …, x_0(k−1), T, x_0(k+1), …, x_0n]`), an n×n system where `x₀ₖ` is a
**constant**, not an unknown with a trivially-satisfied equation. The constraint is
structural where ours is algebraic.

So Step 3 is **not portable to a bordered formulation as a drop-in** — taking it means
taking the substitution with it. The frozen pin's failure mode stands, and *that* is its
fix.

⚠ **A correction the reviewer volunteered on index-2:** they had called admitting the
constraint "research-shaped". [TCF] 1975 solves it, by the same authors — *"the method
yields valid results even when the capacitor voltages and inductor currents are not an
independent set"*, with worked capacitor-loop and inductor-cutset examples, and it shoots
on the **independent states** rather than the full MNA vector. **The priority is
unchanged** — index-2 still produces no silent wrong answer here — but if support is ever
wanted, [TCF] §III is the reference and it is not research.

## `x₀` as the unknown — SHIPPED as `solve(x0_unknown=True)` (2026-09-02)

The formulation change three findings point at: the plain path solves for `x_in`, the
pre-image of a manufactured opening step, while its Jacobian is taken with respect to
`x₀`. See `benchmarks/pss_x0_unknown.py`.

⚠ **The obvious version is the FOURTH design to die on the `(−1)ⁿ` mode — and this time
it was predicted, not discovered.** Dropping the manufacturing step removes the L-stable
opener, so trapezoidal's period map becomes `A_trap^K`, which the L-stability theorem
says is singular at even K. Measured on the model problem *before any shooting code was
written* — `sigma_min(I − A^K)`:

| circuit | variant | K=100 (even) | K=101 (odd) | K=200 (even) |
|---|---|---|---|---|
| Q=20 | trap bare | **0.0** | 6.18e-03 | **0.0** |
| Q=20 | euler bare | 6.31e-03 | 6.32e-03 | 7.33e-03 |
| Q=20 | trap **E-open** | 6.04e-03 | 6.10e-03 | 9.84e-03 |
| RC | trap bare | **0.0** | 9.99e-01 | **0.0** |
| RC | trap **E-open** | 9.99e-01 | 9.99e-01 | 9.99e-01 |

Ten lines of linear algebra, not a build. **This is what the standing "measure before
building" warning is for, and it paid.**

**The design that survives** keeps one L-stable step and moves it *inside* the period:
unknown `x₀`, no manufacturing step, first step order-dropped to Euler — which
`_begin_period` already arranges. Period map `A_trap^{K−1} A_euler`, non-singular at
every K.

**What it buys:**
- **Quadratic convergence**, at even and odd point counts alike: 2.9e+00 → 9.1e-08 →
  1.1e-14 in two iterations on the Q=20 resonator.
- **Van der Pol, item 5's payoff case: −47.3 ppm** in four iterations against the shipped
  path's **−73.8 ppm**. That −47.3 is exactly what the throwaway driver in
  `pss_lte_grid.py` reached — and for the same reason: its unknown was `x₀`.
- **Euler is unchanged to the digit** (8.83917 / 12.28401), the control — euler's
  manufacturing step *is* an Euler step, so the two maps coincide.

⚠ **What it costs, and why it is not shipped on a hunch.** Moving the Euler step inside
the period degrades the **orbit**, not just the opening. Against the analytic 20 V:

| npts | shipped (`x_in`) | `x₀` formulation |
|---|---|---|
| 100 | 20.01273 | 19.76939 |
| 200 | 20.02208 | 19.96123 |

18× worse at 100 points, 1.8× at 200 — the gap closes as the single first-order step is
diluted, and the *order* is preserved, but the constant is real.

**So the trade-off is grid-dependent and both directions are measured:** a uniform grid
on a benign circuit favours the default; a non-uniform grid on a stiff one favours `x₀`.
Shipped as an **option**, `solve(x0_unknown=True)`, the same shape as `matrix_free`.
`method='gear'` refuses it — a solved-history formulation already solves for `x₀` and
`x₋₁` directly, so the flag would be a no-op, and a no-op flag the caller believes in is
worse than an error.

⚠ **The −47.3 ppm is NOT the formulation alone, and the first version of this entry said
it was.** Measured on the *same* grid, the `x₀` formulation gives −73.8 subdivided and
−47.3 raw. The gain comes from the formulation making the **opening-step subdivision
unnecessary** — that subdivision exists only to protect a manufactured step taken from an
iterate that may be far from the orbit, and it costs accuracy. With `x₀` the unknown the
first step starts *on* the orbit and the raw grid is solvable. Two changes, one enabling
the other; the confound was visible only by running both grids.

⚠ **Two frame errors surfaced only end to end**, neither predicted by the prototype:
- the **phase pin** was taken one step after the seed (`_x1[phase_k]`), right when `x_in`
  is the unknown and wrong when `x₀` is. With a fine opening step the two nearly agree
  and it hides; on van der Pol's grid (opening step 1.4845 against a median 4.6e-04) it
  pinned a value the orbit need not attain and the solve died.
- the **replay** dropped `X[0]`, correct when it is the pre-image and wrong when it is
  `x(0)` — one column short of `times`.

Both are the same shape as the bugs this campaign keeps finding: a quantity that is right
in one frame carried into another.

## ⚠ The plain-path replay walked a shifted grid (fixed 2026-09-02)

The two traversals pair `(t, h)` differently: `_traverse_solved_history` walks
`times[1:]` with `hs[_j]`, while the plain `_traverse` takes the **manufacturing step**
at `(times[0], hs[0])` first and only then walks `times[1:]` with `hs[_j]` — so the step
after the opening one uses `hs[0]` again, not `hs[1]`. The replay set `walk = times` and
indexed `hs[min(_j, …)]`, pairing `times[k]` with `hs[k]` from k=1 on: **off by one**.

Closure `|x(T) − x(0)|` of the RETURNED waveform, Q=20 resonator at resonance, 200
points, `converged = True` in every row:

| grid | trap (before) | trap (after) | gear (control) |
|---|---|---|---|
| uniform | 5.61e-13 | 5.61e-13 | 4.62e-14 |
| 4:1 | **1.70e-02** | 1.44e-11 | 5.33e-15 |
| 16:1 | **4.88e-03** | 3.55e-14 | 1.78e-14 |

**A uniform grid hides it completely** — every `hs` is the same number — which is why it
survived every uniform test here. Gear closes either way because it takes the
solved-history branch, whose pairing was already right; that made it a clean control.

Now built as explicit `(t, h)` pairs rather than two sequences indexed in parallel,
because the parallel indexing *is* the bug and a pair cannot be misaligned by one.

⚠ **The headline ppm figures were NOT affected, contrary to what was assumed when this
was found** — including by me, in the triage that recommended holding it for a decision.
The period comes from the shooting **solve**, not the replay: van der Pol still reads
**−73.8 ppm** (trap) and **−100.6 ppm** (gear) after the fix. What the shift corrupted is
the returned waveform, `pss.times`, and the LTE reports — and every recorded LTE figure
was taken on a **uniform** grid, so **no recorded number moved** and the full suite is
unmoved. The item was worth fixing on its own terms (a returned waveform that is not the
solution its residual was driven to zero on), not because it invalidated results.

## Matrix-free against a SPARSE baseline — the gate holds (settled 2026-09-02)

The round-one worry: every recorded matrix-free speedup had a dense baseline on *both*
sides, because the path reached past `linearsolver=` to `scipy.linalg` and
`_step_sensitivity` reached past it to `toolkit.linearsolver`. Both now go through the
caller's solver, so the comparison can be run twice. On a quiet box (load 1.7),
single-threaded, three Newton iterations — wall seconds:

| m | dense/dense | dense/mfree | splu/dense | splu/mfree | mf vs dense | mf vs splu |
|---|---|---|---|---|---|---|
| 110 | 0.713 | 0.646 | 0.739 | 0.682 | 1.10× | 1.08× |
| 242 | 2.166 | 1.596 | 2.045 | 1.531 | 1.36× | 1.34× |
| 502 | 9.301 | 6.614 | 7.897 | 5.298 | 1.41× | 1.49× |
| 1002 | 51.083 | 26.049 | 40.804 | **17.674** | 1.96× | **2.31×** |

**The m≈250 gate holds.** Matrix-free gains the same or *more* against SuperLU than
against dense LAPACK at every size — the worry is not borne out. The two are
**complementary, not alternatives**: `splu/mfree` is fastest at every m ≥ 242, and at
m=1002 it is 17.674 s against the 51.083 s this analysis shipped with — **2.89×**, where
a sparse solver alone buys 1.25× and matrix-free alone 1.96×. `splu/dense` never wins a
row, so "just use a sparse solver" is not the better recommendation at any size.

⚠ SuperLU is a small *loss* at m=110 — the usual sparse-bookkeeping crossover. And
`rc_ladder` is tridiagonal, the most favourable structure a sparse solver can be handed,
so this **overstates** the sparse baseline's advantage — which is the direction that
makes the conclusion safe.

## External review, second round — what was verified and what was done

Every finding below was **reproduced here before being acted on**; where my measurement
disagreed with the review, mine is recorded and the disagreement stated.

### Fixed

| # | Finding | Evidence |
|---|---|---|
| 1b | **`_is_autonomous` sampled ~9 points** (`times[1::len//8]`) while calling itself exact. A `VPulse` placed *between* samples: 40% and 20% duty read correctly, **5%, 1% and 0.5% all read AUTONOMOUS** — which routes the solve to the free-period system and **discards the caller's period**. `DEGENERATE_PERIOD_FACTOR` cannot catch it (it tests |T|). Now evaluates every grid point; cost is N `u()` calls once per solve, with early exit. |
| 1c | **Two reference nodes in one analysis.** `__init__` fixes `self.irefnode` (used by the traversal); `solve()` computed its own from `refnode=` (used to reinsert the zero row). Never compared. Now refused — there is no answer to give, since the two disagree about which variable was eliminated. |
| 1e | **Grid ratio > 1+√2 silently demotes Gear-2 to first order.** Measured on a Q=20 resonator *at resonance*, alternating 3:1 grid: **7.998 / 11.429 / 14.550 / 16.853** at N=100/200/400/800 against a uniform **19.915 / 20.010 / 20.022 / 20.024** and an analytic 20 V — 60% low, `converged=True` every time. (The review reported 7.97/11.40/14.52/16.83; agreement to three digits.) `_period_grid` now warns, naming the *ratio* — refining does not help, since a refined 3:1 grid is still 3:1. |
| 3a | **The `(−1)ⁿ` obstruction had the right conclusion and the wrong cause.** Not the `iq` recursion: trapezoidal is A-stable but **not L-stable**, so it maps `null(C)` by exactly −1 per step. `A_trap = (C/h+G/2)⁻¹(C/h−G/2) → −I` on `null(C)`; Euler's → 0. Verified to 1e-9: Q=20 has m=4, rank(C)=2 and **exactly 2** modes at −1 (trap) / 0 (euler); RC has m=3, rank(C)=1 and 2 of each. This makes "the plain path is correct for trapezoidal" **provable**, and adds the corollary that *any* L-stable opening step works, rescuing exactly `m − rank(C)` modes. |

### Documentation corrected in place

- **The "true Newton for both methods" claim was false for the plain path.** The true
  `dF/dx_in` is **singular on every circuit tested** — rank 1/3 (RC), 2/4 (Q20), 1/3
  (nonlinear), σ_min exactly 0 — so an exact Newton for that formulation *does not
  exist*. What `_traverse` returns approximates `dF/dx₀`, a different and well-posed
  derivative, which is why the path works at all and why it converges **linearly**
  (trap's autonomous residual falls at a constant ratio ~0.076).
- **The "~30%" Jacobian error is a single-circuit number and misframed.** Measured:
  1.56 (RC), 0.0066/0.0037 (Q20 euler/trap), 4.73 (nonlinear) — 0.4% to 470%, not 30%,
  and it is a *different operator*, not a perturbation.
- **`max_lte_seam` is `None` on all three shipped methods** — the seam is collected only
  for a companion reaching two charges back (Gear-2 alone), and Gear-2 always takes
  solved-history, which clears the flag. Its `_limits` warning entry is dead too.
- **"driving Transient buys limiting, PCNR, breakpoints and the continuation rescue" is
  1-for-4.** Limiting reaches. `PSS(pcnr=True)` raises `KeyError`; `_rescue_solver` and
  `cir.next_event` appear nowhere in PSS, because they are armed inside `Transient.solve`
  which PSS never calls — the same structural fact as the TLine refusal.
- **`_spectral_report`'s split is an ODE count.** Index-1 MNA with `d = rank(C) < m` has
  `d` physical, `d` parasitic and `2(m−d)` structural zeros, so both arrays are
  mislabelled on real circuits (`parasitic_roots` comes back identically zero).
  ⚠ `spectral_radius` is unaffected — the physical multipliers have the smallest block
  split and are always inside the first `m`. And it cannot be fixed by magnitude:
  Gear-2's parasitic roots are ~1e-95, indistinguishable from a structural zero.

### Index-2: detection is cheap, and a refusal would be WRONG (settled 2026-09-02)

The review proposed detecting `index > 1` (projector test: `NᵀGN` singular for `N` a
null basis of `C`; measured 29.7 ms at m=112, 58.7 ms at m=502) and refusing, naming
`method='gear'` as the workaround. **Both halves of that recommendation fail on wider
evidence.** Four index-2 topologies, all confirmed index-2 by the projector test:

| circuit | trap | euler | gear |
|---|---|---|---|
| CV-loop | F | F | **T** |
| CV-loop2 | raises | F | **F** |
| LI-cutset | **T** | **T** | **T** |
| LI-cutset2 | raises | raises | **F** |

- **`index > 1` is not predictive of failure.** On `LI-cutset` all three methods
  converge, so a refusal keyed on index would reject circuits that work.
- **Gear is not the workaround.** It succeeds on `CV-loop` — the one circuit the review
  tested — and fails on two of the other three. A refusal naming gear would have sent
  users to a method that does not generalise.
- **And the failures are LOUD.** Every index-2 solve that reports `converged=True` is
  *correct*: worst relative error against a settled transient, over every state
  including the algebraic ones, is **6.6e-05** (CV-loop/gear) and 1.3e-05 / 1.7e-04 /
  7.4e-05 (LI-cutset trap/euler/gear). The review's "89% error on the v-source branch
  current" did not reproduce as a *converged* answer here — a non-converged PSS still
  returns a result, and that is where such a number comes from.

So: **no refusal, no warning, no detection code** — an unused detector is a dead knob,
and this tree has a standing guard against those. What went in instead is a test pinning
the property that would matter if it broke: *a converged index-2 solve agrees with a
transient*. A regression to a silently wrong algebraic variable is what it catches.

⚠ One caveat kept from the review, because it is right and independent of the above:
admitting the constraint on the plain path is **not an increment** — it needs consistent
initialisation of the algebraic variables at the period boundary, which the manufactured
flat history destroys. That is the same change as "make `x₀` the unknown".

### Disputed, then retracted by the reviewer — my evidence stands

- The review reported that my Poincaré condition numbers are identical between the pin
  and the orthogonality row (edge 1.000×). **Re-measured: they are not.** Pin
  1.21e3 / 3.02e2 / 8.23e1 against orthogonality 1.95e2 / 5.99e1 / 3.69e1 at
  200/800/3200 points — a 6.2× / 5.0× / 2.2× edge, as recorded. The review's reasoning
  (alignment `|f̂[k]| ≈ 0.9999`, so "the pin row IS the orthogonality row") conflates
  *parallel* with *equal*: `e_k` is a unit vector and `f̂` is dense, and that difference
  is what moves the conditioning. **The decision is unchanged either way** — the
  orthogonality row still buys nothing that decides anything.

  ⚠ **The reviewer re-checked and retracted §7 in full.** Their harness took the code's
  analytic `J` and swapped only the border row, so both readings were the conditioning
  of the *same* operator — and on the plain path that operator is in the wrong frame
  (their own finding), so its defect dominated `cond` in both cases. Identical numbers
  were guaranteed by construction. A defect found by one lens invalidated a measurement
  taken by another; worth remembering as a shape, not just as an episode.

## ⚠ Gear-2's period column was exactly 3/2 too large (fixed 2026-09-02)

Raised by external review, reproduced here before acting. The autonomous Jacobian's
`d/dT` column against central differences:

| method | npts | ratio code/FD | relerr before | relerr after |
|---|---|---|---|---|
| trap | 100→400 | 0.99993 → 1.00000 | 7.1e-3 → 1.8e-3 | unchanged |
| gear | 100→400 | **1.4859 → 1.4972** | **0.49** | **1.5e-9** (FD floor) |

**Cause: a result carried across contexts.** `residual_dh` is Fang's `p` — `d/dh_m` with
the past steps *held fixed*, and `bdf2_alphas_dh`'s own docstring says so ("h2 is a past
step and is held fixed"). Correct for the coupled time-stepping method it was written
for; wrong for shooting, which rebuilds the grid at the current `T` so every step is
`c_k·T` and `h_{n-1}` moves too. Euler and trapezoidal never noticed — their
coefficients depend on `h_n` alone, so the partial *is* the total. **That is exactly why
only Gear-2 was hit, and why a one-method test would have found nothing.**

**The fix is one term and no new derivative.** The `alphas` are homogeneous of degree −1
in the step sizes (they must be — `iq` approximates `dq/dt`), so Euler's theorem gives
`Σⱼ hⱼ ∂aₖ/∂hⱼ = −aₖ` and hence `T·d(iq)/dT = −Σₖ aₖ q_{n−k}` exactly. Verified
numerically for variable-step BDF-2 to 6.4e-08 (the FD floor). Implemented **once** in
`Integrator.companion_dT` from `companion_coefficients` rather than as three more
per-method partials — which is what that file explicitly warns against.

**What it buys:** quadratic convergence restored. Gear on the phase circuit now runs
`1.69e-1 → 1.06e-2 → 3.78e-6 → 3.48e-12 → 2.76e-15`.

⚠ **Trapezoidal's column is still O(h), and that is a different defect.** Its autonomous
route is the *plain* path, whose `Pt` opens at zero although `x₀` is manufactured by a
step of size `c₀·T` and does depend on T. Measured: trap's Newton is visibly **linear**
(`3.91e-3 → 3.14e-4 → 2.66e-5 …`, ratio ~0.076 held constant) where gear's is quadratic.
That is the plain path's seeding, not this — see the review's §3.

## ⚠ Hidden state: PSS silently replaced a TLine with a short (fixed 2026-09-02)

Raised by external review, reproduced here before acting. On a quarter-wave open stub:

| | PSS | Transient |
|---|---|---|
| converged | **True** | — |
| spectral_radius | **0.0** | — |
| warnings | **none** | — |
| amplitude v(b) | **0.999969** (line absent) | **0.244201** (line active) |
| TLine history | **0** | 4 |

**Mechanism.** `TLine.history` is per-element state outside `x`, filled by
`cir.accept_step` — which a forward transient calls at every accepted step and which
PSS never calls, because it drives `solve_timestep` directly. `TLine.G`/`u` branch on
an empty buffer and stamp a **DC short**. This is Kundert's hidden state: the period
map is not a function of `x₀`, so the periodicity condition is imposed on an incomplete
state and Newton converges to something meaningless. `spectral_radius = 0.0` is the
tell — a circuit reporting no state at all. ⚠ **The tree already knew this defect by
name**: `transient.py` carries a comment about the identical failure being found and
fixed *there*.

**Why filling the buffer is not the fix.** Calling `accept_step` per step would populate
the history and make `phi` genuinely history-dependent — so the monodromy would be the
derivative of a neighbouring problem, the exact failure `_begin_period` exists to
prevent. `_begin_period` resets what is *in* `x`; that is the right scope for the
integrator rings and the wrong one for state outside the vector, and no reset of the
rings turns a period map that is not a function of `x₀` into one that is.

**What was done: refuse, and say why.** `Circuit.hidden_state` is declared on the
element (`TLine` sets it), `hidden_state_elements()` finds them by instance path, and
`PSS.solve` raises naming the elements and pointing at `Transient`. ⚠ Overriding
`accept_step` is deliberately *not* the criterion — `_WrapEvents` and the HDL `@cross`
wrapper override it and feed only `next_event`, which PSS ignores because it imposes
its own grid; a diode's `_vlim` is hidden state too but self-erasing. Sniffing the
method would refuse ordinary circuits.

No PSS test used a `TLine`, which is why this was never caught here.

## External review, 2026-09-02 — first round, three findings, all confirmed

A read-only review against the shooting/PSS literature (Telichevesky/Kundert/White
DAC'95, Kundert ICCAD'97, Nastov, ter Maten) raised three; each was verified here
before being acted on.

1. **`nrsolver` / `linearsolver` / `scaler` were accepted and silently discarded.**
   Declared on the base `Analysis`, so `PSS(cir, linearsolver=SuperLUSolver())` was
   valid input — and `_transient()` forwarded the tolerances and not these, so the
   inner `Transient` resolved to `DenseSolver`/`StandardNewton` regardless. **The third
   time this class has taken a parameter it never read** (`method` declared and never
   read; `analysis='PSS'` matching nothing), and invisible the same way each time: the
   constructor validates the name, nothing checks the boundary. Fixed; the test asserts
   the *resolved object*, since a test that the string `linearsolver` appears in the
   source passes with the parameter dropped.

2. **The matrix-free path was hardcoded to dense LAPACK** — and so was the dense
   propagation, one layer further in (`_step_sensitivity` → `toolkit.linearsolver`).
   Both now go through the caller's solver, via a `factor()` hook the `LinearSolver`
   ABC already declared and no subclass implemented; `DenseSolver` and `SuperLUSolver`
   implement it now. ⚠ **The measurement half is OPEN.** Every recorded speedup has a
   dense baseline on both sides, so the m≈250 gate may move. The attempt to re-measure
   read load average 31–43 with a second session benchmarking the same machine, and its
   readings disagreed with earlier ones by 25–30% on the dense/dense pair alone — not
   quotable in either direction. Harness: `benchmarks/pss_matrix_free_sparse.py`.

3. **GMRES's `info` was discarded**, so a Krylov breakdown became the generic outer
   "No convergence" with nothing naming it — in a file whose standard is that failures
   name themselves. Also: scipy's `maxiter` counts *restart cycles*, so the old
   `restart=200, maxiter=400` was a worst case near 80 000 matvecs, each a full replay
   of the period, with no diagnostic while it was being spent. Now surfaced, bounded
   (`KRYLOV_RESTART` × `KRYLOV_MAX_CYCLES`), and exercised by a test that starves the
   solve on purpose.

Two things the review explicitly said were **not** findings and should not be "fixed":
the `conv_f` ordering in `_matrix_free_newton` (it mirrors `analysis.fsolve`; changing
one alone creates the inconsistency), and the absence of a GMRES preconditioner
(correct as-is — `I − M` clusters at 1, which is why k tracks slow modes, and the
Gourary adaptive-preconditioning work is harmonic balance and does not transfer).

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
