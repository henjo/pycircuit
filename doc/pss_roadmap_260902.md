# PSS / shooting — forward roadmap (2026-09-02)

`pss_shooting_roadmap_260901.md` is the **record**: what was built, what was measured,
what was falsified. This file is the **plan**: what is left, what it would cost, and the
gate each item has to pass before anyone starts it.

Nothing here is in progress. Every item is either unstarted or deliberately closed, and
the closed ones are listed because the expensive mistake is re-opening them.

**Standing rule for everything below:** each item names a *gate* — a measurement that
decides whether to build. Four items on this list were killed by their gate for less than
an hour's work each; two were killed after being built. Run the gate first.

---

## A. Capabilities — unbuilt, entry points known

### A1. PAC (periodic AC) — closest to hand

**The operator already exists.** Telichevesky, Kundert & White (DAC 1996, pp.292-297)
reach the iterative form by "reinterpreting the use of `L^-1` … as a preconditioner",
where `L` is the block lower-bidiagonal transient discretisation — **that preconditioner
is our per-step factored solve sequence**. The PAC operator is `I + alpha(f) H`, with
`Hv` one replay of the stored per-step factors, i.e. `_traverse_factored` +
`_monodromy_matvec`, both shipped and tested.

Recycling across the frequency sweep is **two scalars**. Their Theorem 1: the Krylov space
spanned by `{p0, (I + alpha H)p0, …}` is independent of `alpha`, so
`beta = alpha(f_hat)/alpha(f_s)`, `gamma = 1 - beta` converts a matvec at one frequency
into a matvec at another. The first frequency pays; the rest rescale.

⚠ **The withdrawn PAC's 419.5 GiB came from FORMING the operator.** This method never
forms it. That is why the withdrawal does not carry over.

*Reported:* "up to forty times faster than the standard optimized direct methods", on a
thousand-node RF mixer.
*They use GCR, not GMRES; whether that matters here is unexamined.*

**Gate:** none needed to start — the operator is built and tested. The honest first step is
a single-frequency PAC against a dense reference, then the sweep.
**Size:** the largest item here that is not research.

### A2. PPV (perturbation projection vector)

**Unblocked** by `_monodromy_matvec_transposed` (shipped, agrees with dense `M^T` to
1.8e-15, costs 0.75x a forward matvec).

Demir & Roychowdhury (TCAD 22(2) 188-196) Thm II.4: augment the reverse Jacobian and solve
once. `J~_r (x; y) = (0; N)` gives `x` = PPV, `y = 0` — and `y` coming back zero is a free
correctness check. `q = ` sampled `C(t) u_1(t)` with `u_1 = xdot` from the PSS solution,
which we already have.

⚠ **Do NOT build this on `_spectral_report`'s eigen-split.** D&R 2003 is a paper about why
not: "inaccuracies … often corrupt the oscillatory-mode eigenvalue of 1 to the extent that
it cannot be distinguished from other eigenvalues … a potentially large number of candidate
PPVs are typically found and one chosen using heuristics … not entirely reliable." We
independently measured that breakdown — at ~2 points per cycle the parasitic root has the
*smaller* block split and is labelled physical.

⚠ **`J_r` is NOT `J_f^T`.** `J_f = Omega D_C + D_G`, `J_r = D_{C^T} Omega - D_{G^T}`: the
operator differs *and* the sign on `G` differs. No free transpose at the block level.

**Gate:** verify the block-cyclic `Omega` against the BVP operator already written by hand
for the `func_solved_history` cross-check (it matched PSS to 4.0e-13). Confirm `J_r`'s two
sign/transpose differences against p.192 before coding — the algebra was read off a page
image and `pdftotext` drops it.

### A3. pnoise — ⚠ UNBLOCKED 2026-09-02, the paper is on this machine

Needs the **adjoint** formulation: pnoise is many-to-one (hundreds of sources, one output),
so a forward solve costs one solve *per source* and recycling does not help, because the
RHS is what changes.

Shares its machinery with A2 — both bottleneck on the transposed replay, now built.

~~**Blocked:** Okumura et al. 1993 not in the library.~~ **It was acquired 2026-09-02**
and is at `~/docs/02-oscillator-noise-jitter/Okumura-Tanimoto-Itakura-Sugawara-1993-…
(TCAS-I).pdf`. Relayed by the docs session, not read here — the quotes below are theirs.

Its eq. (36) is `T^T X^a = d_j`: **one transposed solve per output**, chosen for exactly
our asymmetry — "it is efficient to use the adjoint method … *because circuits have many
noise sources*". With `M^T` a reverse replay (26d43da), the machinery is in place.

⚠ **Two limits on what the paper gives.** The adjoint is a *brief passage*; it cites [34]
= Vlach & Singhal, an out-of-print textbook that could not be obtained — **and is not
needed**: Trick, Colon & Fan 1975 (in `07-`) argues the classical adjoint-*network*
construction, which is what [34] would describe, is the wrong route here ("no convolution
of the adjoint circuit response with the original circuit response is required"), and
Director & Rohrer 1969 (in `07-`) is the primary source for the classical treatment.
What Okumura 1993 *does* supply is the part we need: LPTV transfer functions **from a
shooting solution**, and cyclostationary sources with aliasing handled by accumulation
under a ratio-test stopping rule.

**Gate: open.** Start from Okumura's §III, not from the DAC'96 deferral.

### A4. Warm start — ⚠ SHIPPED 2026-09-02 as `tstab=`; the *automatic* criterion is what remains

De Luca, Bolcato & Schilders (2019, TCAS-I) frame our exact situation: "none of the works
in the literature addresses the relevant problem of automatically identifying such a proper
initial solution. Usually, heuristics are used." That heuristic is what our
non-convergence warning currently tells the user to do by hand.

**Their Algorithm 1 is `_monodromy_matvec`** — same recursion, same per-step structure, run
during pre-integration to find the first period after which the iteration is inside the
contraction region. "Non-invasive … can be implemented with little effort", ~`M * nnz(A)`
flops.

**Gate — RUN, and it split.** `benchmarks/pss_warm_start.py`.

⚠ **The value is confirmed and it is large.** Seeded near the unstable DC point — the
trivial-root basin the literature describes — van der Pol fails cold and converges after
pre-integration:

| circuit | cold solve | periods needed |
|---|---|---|
| μ = 1 (strongly attracting) | `LinAlgError` | **1** |
| μ = 0.05 (high-Q) | not converged | **~24** |

The ~24 is the `1/mu` envelope time constant, so the count is a property of how strongly
the cycle attracts, not of the seed. A large-amplitude seed is far easier — from 4× and
even 20× the orbit amplitude, one period suffices at μ = 1.

⚠ **But the stopping criterion is the open part, and it is blocked on the paper.**
Carrying a probe through the variational system each period (`u^{k+1} = M_k u^k`) and
watching its direction settle **does not identify the handoff**: the drift is *small*
(1.4e-02) while the solve still fails and *large* (1.2e-01) once it succeeds. It moves the
wrong way. The reason is diagnostic — near the DC point the monodromy is nearly constant,
so the probe settles into its own eigenvector and reports "converged" while the state is
stuck at the trivial root.

⚠ **The obvious alternatives share the defect**, which is why this is not worth guessing
at: "the state stopped changing period to period" and "the shooting residual is small" are
both *also* true at the trivial root — it is a fixed point of the period map, so it passes
every periodicity test. Telling them apart needs something amplitude-like, which is what
`DEGENERATE_PERIOD_FACTOR` does for the period.

⚠ **BOTH PARAGRAPHS ABOVE DIAGNOSE THE WRONG THING. Corrected 2026-09-02** from the
paper itself (relayed by the docs session). Two independent errors, and the second is the
one that actually explains the measurement:

1. **The definition was guessed wrong.** `u^k = x^k − x^{k+1}` (eq. 4) is the **measured**
   shooting error between consecutive pre-integration periods, and the criterion compares
   it against the **Jacobian-propagated prediction**, `u^k = J_φ(x*) u^{k−1}` (eq. 8). It
   asks *"does the error still evolve the way the linear model says it should"* — a test
   of being **inside the linear region**. The guess here — `ũ^k` = the previous iterate —
   made it `||u^k − u^{k−1}||`, a test of whether the error has **settled**. Different
   question, and the one that gave the backwards reading.
2. ⚠ **And the gate was run on a circuit class the paper excludes.** Its title and §I say
   **non-autonomous**; its conclusion offers autonomous circuits only as conditional
   future work ("*may be easily extended* … *could be effective if* the time evolution of
   the period is characterized by a linear region"), conditioned on a variable their
   formulation does not have. **Van der Pol is autonomous**, and its trivial root is a
   genuine fixed point of `φ`. So the diagnosis above — the probe settles and reports
   converged while the state is stuck at the trivial root — is the *predictable* result:
   the criterion detects proximity to **a** fixed point, and on an autonomous circuit the
   DC point **is** one. No variant of it can separate them there, and the paper never
   claimed otherwise.

**Revised gate: a DRIVEN circuit**, where there is no trivial root to be attracted to and
"inside the linear region" is unambiguous. That is the class the paper addresses and a
large share of real PSS use. For the *autonomous* trivial-root basin the right tool is the
probe technique (B5, Bizzarri et al., in `01-`), which pumps energy in so the solve cannot
fall to the DC point — **or, per B5, a device `ic`, which is what that paper's own authors
use on their flagship high-Q oscillator.** A4 and B5 answer different halves of the seeding
problem and neither substitutes for the other, but the cheap half is the `ic`.

⚠ **The measured value above stands** — it was never criterion-dependent.

**What shipped instead, and why first.** `PSS.solve(..., tstab=<seconds>)` runs a plain
`Transient` from the seed and hands its final state to the shooting solve. Directive:
*"We need to add warm start as an option. It is in every commercial tool."* An explicit
option needs no criterion, so it was not blocked on any of the above — and the user's
field experience is that the automatic version is where the trouble is: *"Spectre has an
automatic tstab criterion but it does not work properly on circuits with even moderate Q,
and does not work [on] high Q circuits"* — **opinion, offered as such, not measurement**,
but pointing the same way as the μ = 0.05 row above. The option is therefore the primary
interface; an automatic criterion is an addition on top of it, not a prerequisite.

⚠ **The one limit `tstab` cannot pass**, asserted in the suite rather than left to be
rediscovered: with `x0=None` the seed is the **operating point**, and on an autonomous
circuit that is an *equilibrium* — a transient started exactly there never leaves. The
pre-integration needs somewhere to go: an `x0` off the equilibrium, or a device `ic`.
This is the same reason it does not substitute for B5.

### A5. Envelope-following — last

Linaro et al. (OJCAS 2020) apply EFM to the *variational* problem, with a
composition/"dragging" trick (binary powering of `Phi`) as the headline. Reported: 0.4% as
many cycles.

⚠ **Do not start before the LTE/smoothness question is settled on our terms.** They impose
an LTE tolerance on the variational envelope that governs step size. Our docstring argues
at length that LTE cannot be a controller under shooting, because an adaptive step sequence
makes `phi` a different discrete map for each `x0` and costs the Newton its rate. Their
adaptivity is at the *envelope* level, which may or may not be the same objection. **This
is unresolved, not refuted.**

---

## B. Formulation decisions — measured, awaiting a call

### B1. Make `x0_unknown` the default on non-uniform grids

Shipped as an option. The evidence says it wins exactly there and loses on uniform grids:

| | default | `x0_unknown` |
|---|---|---|
| Q=20 @ 100 pts (analytic 20 V) | 20.01273 | 19.76939 |
| Q=20 @ 200 pts | 20.02208 | 19.96123 |
| van der Pol, LTE grid | −73.8 ppm | **−47.3 ppm** |

**Gate:** a rule that picks correctly without the caller knowing. "Non-uniform" is not
quite it — the gain came from the formulation making the opening-step *subdivision*
unnecessary, so the real predictor is whether the grid opens coarse.

### B2. theta = 1/2 + Ch — the fifth trapezoidal design

Houben (2003): "a pure BDF method should not be used" for autonomous oscillators, and his
theta-method biases theta off 1/2 by `Ch`, making it "virtually second order" while damping
the DAE modes.

⚠ **It is the only design that removes the opening step rather than depending on it**, so
it composes with `x0_unknown` instead of competing. Biasing theta damps `null(C)`, so the
`(-1)^n` obstruction dissolves rather than needing annihilation by a special first step.

⚠ **Only worth it alongside B1** — on its own it replaces a mechanism that already works.
⚠ **Houben gives a two-sided constraint on `C` and no recipe.** Small enough not to damp
the physical oscillation over a period, large enough to kill the DAE modes.

**Gate:** the same falsifier that killed the fourth design — `sigma_min(I - A_theta^K)` at
**even** K, plus the Q=20 peak against 20 V. Ten lines. Run it before writing anything.

### B3. The Aprille & Trick substitution formulation

Their unknown vector *substitutes* the period for the pinned coordinate —
`v = [x_01, …, x_0(k-1), T, x_0(k+1), …, x_0n]` — an n x n system where `x_0k` is a
**constant**, not an unknown with an equation.

**It would subsume two open things at once:** the per-iteration phase re-selection (which
we built and reverted — see C3) and the far-seed failure where a frozen pin names a value
the orbit never attains.

**Gate:** it changes the outer system's shape, so the honest first step is the same
far-seed case that defeated the bordered version: van der Pol at 4x/10x/30x amplitude.

### B4. Index-2 support

[TCF] 1975 §III solves it, by the same authors as the method — capacitor loops and inductor
cutsets, worked examples — by shooting on the **independent states** rather than the full
MNA vector. **Not research.**

⚠ **Priority unchanged and low:** index-2 produces no silent wrong answer here. Every
failure reports `converged=False`; every converged run is correct to <=1.7e-4.

⚠ **And it is not an increment.** It needs consistent initialisation of the algebraic
variables at the period boundary, which the manufactured flat history destroys — so it is
the same change as B1/B3, not a separate one.

### B5. The probe technique for the trivial-root basin

Bizzarri et al., "Probe Based Shooting Method …", **already in the library**. A periodic
voltage source feeds the oscillator until its own current reaches zero; the probe can then
be removed without changing the steady state. Named in our trivial-root warning since
2026-09-02.

⚠ **It widens the basin; it does not remove seed dependence** — the paper says so itself:
"convergence is still not always obtained for high-Q oscillators as long as the initial
estimate is not close enough."

⚠ **PRIORITY DOWN 2026-09-02, and the paper is what lowered it.** Read past its §III (docs
session, `~/docs/pycircuit-probe-analysis.md`) and the authors do **not** use the probe to
get convergence on their own flagship high-Q example. On the Pierce crystal oscillator:
*"it is easy to assign a tentative current to the Ls inductor modeling the crystal and
obtain convergence in a few iterations (**we did this with conventional SH**)"*. That is
exactly the device-`ic` workaround `tstab`'s limit already points at, so **basic
convergence on a crystal oscillator costs one sensible `ic`, not a new element and a solve
mode.** An earlier reading of this entry — that the probe is the only thing that can inject
the energy the equilibrium seed lacks — is **wrong**, and the same paper refutes it.

**What it does buy** is the *sweep*: unstable limit cycles, coexisting solutions, and a
stability verdict — none of which anyone has asked for. It is also not turnkey: probe
placement is circuit-specific, and seeding `omega_p` needs **the eigenvalues of the system
Jacobian at the equilibrium**.

**Worth remembering separately:** its 2×2 power-flow test (`P = dv_R dy_R + dv_I dy_I`)
skips the eigenvalue computation — but it is **one-directional**. Only `P > 0 ⇒ unstable`
is proven; the authors state plainly that the converse and `P < 0 ⇒ stable` "have not been
proven", only tested. A cheap *instability* detector, not a replacement for
`_spectral_report`.

**Gate:** only if unstable cycles or multiple coexisting solutions are actually wanted.
Not for seeding.

---

## C. Closed — do not re-open without new evidence

Each of these cost real time. They are recorded so the next reader spends none.

| # | Item | Why it is closed |
|---|---|---|
| C1 | **Trapezoidal exact-Jacobian reformulation** | Four designs dead on the same `(-1)^n` mode. It is a **theorem**: trapezoidal is A-stable but not L-stable, maps `null(C)` by exactly −1, so any period map `A_trap^K` without an L-stable opener is singular at even K. Verified: `m − rank(C)` modes at −1, exactly. |
| C2 | **Poincaré / orthogonality phase row** | Tested and rejected. 2–6x conditioning edge that shrinks with refinement and never decides an outcome. The `argmax` rule is canonical (A&T Step 3) and compares units *on purpose* — normalising picks a row **704x worse aligned**. |
| C3 | **Per-iteration phase re-selection** | Built and reverted: it **regressed the working case** (on-orbit seed went converged → not). Structural — pinning the iterate's own value makes the residual identically zero, so the row constrains only the step. A&T avoid it by having no phase equation at all (see B3). |
| C4 | **Index-2 detect-and-refuse** | Would reject working circuits: `index > 1` is **not predictive** (all three methods converge on an LI-cutset), and gear is not the workaround (fails on 2 of 4). Failures are loud, not silent. |
| C5 | **Saltation correction for switching** | **Falsified.** The monodromy–FD gap falls at exactly 2.00x per doubling = O(h) discretisation, not the O(1) a missing correction leaves. PSS's monodromy is the derivative of the *discrete* map, which has no undefined instant. |
| C6 | **Reusing the inner Newton's factorisation** | Does not exist to reuse — `StandardNewton` factors and discards. A retained one would be `J(x_k)`, not the converged `J(x_{k+1})`: already measured at median 5.5e-9 and rejected as an approximation. |
| C7 | **Gourary adaptive preconditioning for PAC** | Harmonic balance. His own introduction says the shooting case was already solved by Telichevesky et al., because `A' = I` there — which is our case. Solves a problem we do not have. |
| C8 | **Deflation for the trivial root** | "The condition number of the system's Jacobian matrix grows beyond any bound and convergence stalls." Not the shortcut. |
| C9 | **A GMRES preconditioner** | The per-step factored solves *are* the preconditioner, in both DAC'95 and DAC'96. ⚠ Do not carry this forward as "this operator needs no preconditioning" — that is false, and the clustered spectrum is a *consequence* of the implicit preconditioning. |

---

## D. How these items keep failing — the shapes worth checking for

Sixteen claims were overturned across this campaign. Four shapes account for most:

1. **A quantity right in one frame, carried into another.** The `x_in`/`x_0` Jacobian
   frame; the phase pin taken one step after the seed; the replay's `X[0]`; A&T's Step 3
   ported to a bordered row. **Check what the source formulation was carrying.**
2. **A number compared against itself, or read without its validity flag.** The propagation
   share timed on the wrong function; "89% error" read off `converged=False` runs; a harness
   reusing the operator under test.
3. **Asserting the size when the question is the rate.** A constant-factor error is
   invisible to "is it small" and unmissable to "does it fall". Both the A&T period-column
   check and the saltation falsifier are rate assertions for this reason.
4. **Two things each tested alone.** `x0_unknown` x `matrix_free` crashed on a line written
   in the same commit as the feature.

**And one about measurement itself:** this machine runs more than one agent. Check
`ps -eo pid,pcpu,args --sort=-pcpu` and `uptime` before trusting any wall-clock ratio — a
concurrent run moved readings 25-30% on the *same* configuration.
