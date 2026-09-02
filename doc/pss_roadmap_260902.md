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

⚠ **These are not five independent choices.** A1 → A3 is a dependency chain (A3 consumes
A1's `H_l`), A2 shares A1/A3's transposed replay, and A5 is gated on an unsettled question.
The order that falls out is **A1 first**, because two other items sit on it.

### A1. PAC (periodic AC) — ⚠ BUILT 2026-09-02

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

⚠ **THE WITHDRAWN PAC BODY WAS NEVER WRONG, and this is now RUN rather than reasoned**
(docs session, `~/docs/pycircuit-pac-operator.md`). It is the *un-preconditioned* form of
the right operator. Its `L` is block lower bidiagonal and its `B` is confined to the first
`N` rows and last `N` columns, so `(L + αB)v = -u` is `(I + α L⁻¹B) v = -L⁻¹u`; applying
`L⁻¹` is forward substitution through the timesteps, which is `_monodromy_matvec`'s
recursion, and **`H := L⁻¹B` IS THE MONODROMY**.

Checked against our own code, not assumed from the paper:

| check | result |
|---|---|
| `L` block lower bidiagonal (upper block part) | max `0.0` |
| `B` confined to last `N` cols / first `N` rows | max `0.0` outside, both |
| `L⁻¹Bv` vs forward substitution through the steps | rel `7.4e-19` |
| last block vs a monodromy replay | rel `1.3e-15` |
| `(L+αB)v` vs `L(I + αL⁻¹B)v` | rel `4.3e-17` |

So the **419.5 GiB was entirely the cost of FORMING `L`** — nothing else — and the
preconditioner that removes it is the traversal that already exists. **Keep the operator;
do not re-derive it. Never form `L`.**

⚠ **Three plumbing gaps, none of them in the operator:**
1. `_monodromy_matvec` and `_monodromy_matvec_plain` both do `np.asarray(v, dtype=float)`
   — hard-real. PAC needs complex `v`.
2. `_traverse_factored`'s stored factors are local to the `_build` closures and discarded;
   a PAC hook needs them retained.
3. `pss.Jtvec` / `pss.Cvec` — which the dead body reads — are written only by `_traverse`
   and `_traverse_solved_history`, **not** by either `_traverse_factored*`. After
   `matrix_free=True` they are absent or stale. **The matrix-free PAC must take the
   factors, not `Jtvec`/`Cvec`.**

⚠ **Incidental third source for the order-dropped opening step.** Okumura, Sugawara &
Tanimoto 1990: "if the second-order integration method is used, the first point is
approximated by the backward Euler algorithm, in order to solve the start-up problems."
Same choice, same reason, in the PAC context.

**⚠ BUILT.** `PAC.solve(pss, freqs, recycle=True)`. What is left after the preconditioning
is `m x m`: `(I - α M) y_0 = α w(f)`, with `w` the forced response over one period from a
zero initial state. Nothing forms `L` or `B`.

**The gate was run against a reference PAC cannot influence.** On a linear circuit the LPTV
response collapses to the LTI one, so the `AC` analysis is the answer. Rel error at 700 Hz,
**per doubling of the grid** — the rate is the assertion, not the size:

| method | rate | rel @ 250 pts |
|---|---|---|
| euler | 2.00x (O(h)) | 1.40e-02 |
| gear | **4.00x (O(h²))** | 1.65e-04 |
| trap | 2.00x (O(h)) | 4.13e-03 |

⚠ **Trapezoidal is second order and its PAC is first**, and the first hypothesis was wrong.
Not the `(-1)^n` mode — the monodromy's `null(C)` modes sit at **0**, not −1, because the
Euler opening step annihilates them, exactly as C1's theorem says an L-stable opener does.
Not the trajectory either — trap's waveform converges at ~4.2x per doubling. It is the
**manufacturing step**: it lives outside `_traverse_factored_plain`'s loop, so it is not in
`steps` and the source is never applied there — one step of `u` out of `N`. Falsified as
predicted, with euler as the control:

| | rate | rel @ 250 |
|---|---|---|
| trap, plain | 2.00x | 4.13e-03 |
| trap, `x0_unknown=True` | **4.00x** | 1.09e-04 |
| euler, either | 2.00x | 1.40e-02, identical to five digits |

`PAC` warns and names the two ways out (`x0_unknown=True`, or gear's solved-history path).

**The sweep recycles one Krylov subspace** — `A(α) = I − αM`, so the space is `M`'s and is
α-independent (Thm 1). The RHS is *not* shared, so it minimises the true residual over the
span and extends the basis when a frequency needs it. On RC ladders, 24 frequencies:
72→6 / 168→15 / 302→26 matvecs (**~11.6x**), agreeing to 1.2e-13.

⚠ **Two defects the dead body carried, and only the gate found the second.** Its `L` is
backward-Euler-shaped, which is *not* the discretisation for a two-step method (see the
table in C1's neighbourhood: `-L⁻¹B` has ρ = 0 against our 0.8545/0.8412) — and it read its
source vector **positionally**, `u(0, analysis_name)` into a signature whose second
parameter is `epar`, taking the transient source at `t = 0`, which is zero for every
sinusoid. It would have returned zeros at every frequency, silently.

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

⚠ **THE TRANSCRIPTION IS CONFIRMED TWICE** (2026-09-02), though still worth checking
against the PDF before coding. Traversa & Bonani (IET CDS 5(1):46-51, 2011), in extractable
text, give the direct/adjoint pair as `d/dt[C z] - A z = 0` and `C^T dw/dt + A^T w = 0` —
**including the opposite relative sign on the `A` term**, which is exactly the feature that
makes `J_r != J_f^T` and rules out a free transpose. And Demir (IJCTA 2000) Remark 3.1 gives
`v_i^T(t) C(t) u_j(t) = delta_ij`, which is precisely the meaning assigned to the `q` row.

⚠ **AND THE 2003 IMPROVEMENT IS ONE MOVE, worth understanding before building.** Demir 2000
SELECTED the eigenvector by its inner product against `C(0) xdot(0)` (0.2 against 1e-5,
1e-7, 2e-5). His 2003 paper rejects that heuristic outright — "no guarantee that any of the
candidate eigenvectors will be appreciably more orthonormal than the others". The same
vector then changes role: it becomes the augmented row `q`, so **no selection is performed
at all**. The quantity used to *choose among* candidates becomes the *constraint that makes
the candidate unique*. That is why A2 goes to the augmented solve and not to
`_spectral_report`'s eigenvectors.

⚠ **SCOPE, RECORDED AS A CHOICE RATHER THAN LEFT IMPLICIT.** There are two tiers. Demir's
PPV is **phase noise only**: one augmented-Jacobian solve. Traversa & Bonani's is phase +
**orbital** + correlation, needing *all* direct and adjoint Floquet eigenvectors from a
generalised eigenvalue problem — and they say plainly that Demir's approach "considers only
phase noise … neglecting orbital noise which may in some cases become important".

**The DAE caveat picks the tier for us.** As of 2011 the orbital analyses are "currently
limited to … ordinary differential equations … although this is not the most general case
because the modified nodal analysis of circuits, in general, leads to … DAEs. The
Floquet-based analysis of *phase* noise has been formulated for DAEs [Demir 2000], whereas
the extension … for orbital fluctuations … is currently under development." For an MNA
simulator the cheap tier is the one with a formulation behind it. Both papers assume
**index-1** explicitly, which is this solver's operating assumption too (see B4).

⚠ **Traversa & Bonani does NOT port** — its algorithm is a generalised eigenvalue problem on
**harmonic-balance** Jacobians, and there is no HB path here. It also inherits the objection
D&R raise against eigen-based extraction. Flagged so nobody chases it.

**Gate:** verify the block-cyclic `Omega` against the BVP operator already written by hand
for the `func_solved_history` cross-check (it matched PSS to 4.0e-13). Confirm `J_r`'s two
sign/transpose differences against p.192 before coding — the algebra was read off a page
image and `pdftotext` drops it.

### A3. pnoise — ⚠ BUILT 2026-09-02 for STATIONARY sources; cyclostationary is not

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

⚠ **A3 IS NOT A SIBLING OF A1, IT IS A CONSUMER OF IT. Corrected 2026-09-02** on a full
read of the paper (docs session, `~/docs/pycircuit-pnoise-analysis.md`; the first pass had
read only the abstract and the adjoint paragraph). `H_l(·)`, the Fourier coefficients of
the LPTV transfer function, is an **input** to the noise calculation — and it is exactly
what A1 computes. The paper's own structure is sequential: "*first*, a numerical
calculation method for the time-varying transfer function … *next*, a noise analysis
method is proposed for these circuits."

    PSS  →  H_l  (A1 / PAC)  →  S_alias  (A3 / pnoise), adjoint used INSIDE A3
                                                        because the sources are many

So acquiring the reference did not open A3 for parallel work; it established that **PAC is
a dependency, not a preference**. That is also why DAC'96 is titled "Efficient AC *and
Noise* Analysis" — one recycled-Krylov solve serves both. With A1's operator now confirmed
rather than rewritten, the order is **A1, then A3 on top of it**.

⚠ **AND IT NEEDS NO NEW INPUTS — the paper says so explicitly**, and names the shooting
route specifically: "the specific information required for noise calculation are `h_m`,
`S_m`, and `H_l(·)` … for the Fourier components `H_l(·)`, they are directly obtained by
the harmonic balance method, while they are **obtained via `J_m` and `C_m` matrices as
by-products of the final transient analysis in the shooting method**."

| input | what it is | where it already is |
|---|---|---|
| `h_m` | the timesteps | `hs` |
| `S_m` | per-interval white noise density | `CY(x, w, epar)` across the element library; `analysis_ss.Noise` already consumes it |
| `H_l` | LPTV transfer coefficients | "via `J_m` and `C_m`" — `Jtvec` / `Cvec` ⚠ *but see A1's gap 3: the factored traversals never write them* |

⚠ **The cyclostationary model is one stationary source per timestep**, which is why it
lands on the machinery we have rather than needing new machinery: `c(t) = Σ n_m(t) w_m(t)`
with `n_m` stationary and `w_m` a T-periodic non-overlapping window — and "**the number of
integration time points in one period is used for the number of discrete points**". The
windows *are* the steps, each at its own bias-dependent density. That maps onto a **frozen
grid with no adaptation**, which is what PSS has.

⚠ **A stated validity condition, worth knowing before rather than after:** the
modulated-white model holds because time samples "discretized by more than several tens of
picoseconds can be regarded as uncorrelated". **Below ~tens of ps per step the premise
fails**, and nothing in the formulation will say so.

**Two properties to build the tests around, before any physics:**

1. **The truncation bound is hard, not heuristic** (eq. 32): `N, L <= (w_max - w0)/ws` —
   you cannot alias down from above the grid's own maximum representable frequency. The
   abstract's "accumulated until their contributions become negligible" is a ratio test
   operating *inside* that ceiling, **not instead of one**. An implementation carrying only
   the ratio test is missing the bound that makes it terminate for the right reason.
2. **It reduces to the stationary result at `p = 1`** (eq. 33): one window of width `T`
   collapses the whole cyclostationary machinery to `S_alias(w0) = Σ_l |H_l(w0 - l ws)|²`,
   which "is exactly the same as that derived for a stationary noise". **A free first unit
   test** that exercises the `H_l` path without needing any cyclostationary modelling to be
   correct yet.

⚠ **Not read:** §IV–V (three worked examples) and the construction of the window Fourier
coefficients `R_{m,n}`. Flagged rather than glossed — B5 is the entry where a §III reading
was refuted by §V.

⚠ **THE ADJOINT HALF IS BUILT 2026-09-02.** `PAC.adjoint_transfer_row(pss, freq, output)`
returns every source's contribution to one output in **a single transposed solve**, against
one solve *per source* forward. Measured on RC ladders: agreement **9.6e-16** with the `m`
forward solves, and the speedup grows linearly with `m` — 7.0x / 16.9x / 40.0x at
m = 6 / 14 / 32 — which is the shape the identity predicts, since forward is O(m) solves and
this is O(1). Gear-2 only; the reverse replay is.

⚠ **AND THE SIDEBAND ROWS ARE BUILT TOO** (`PAC.adjoint_sideband_row`). `H_l` is a
functional *distributed over the period*, so its adjoint takes an injection at **every
step** — and the initial state is itself a function of the source through the periodic
boundary condition, which is a **second term**, `α W^T z` with `z = (I − αM)^{-T} g`. The
two are comparable in size (140 against 749 at `l = 0`), so an implementation with only the
first returns a plausible number rather than an obviously broken one.

Verified three ways, each against something the adjoint path cannot influence:

| check | result |
|---|---|
| nonlinear (diode), l = 0 / 1 / −2, vs `m` forward driven solves | **3.1e-16 / 4.1e-16 / 1.5e-15** |
| … and the sidebands are real there | \|H₁\|/\|H₀\| = 0.50, \|H₋₂\|/\|H₀\| = 0.13 |
| linear circuit: every l ≠ 0 | **~10 orders below H₀** |
| `H₀` vs the LTI transfer function, per doubling | **4.06x / 4.03x / 4.02x** = O(h²) |

⚠ **A CONVENTION ERROR THAT ONLY THE LINEAR CASE COULD CATCH.** The first version
decomposed `y(t)` rather than `v(t) = y(t)e^{-jωt}`, the part that is actually T-periodic.
It agreed with a forward reference written the same way to 1e-15 — and reported
\|H₁\| ≈ \|H₀\| on a *linear* circuit, which has nothing to convert with. Self-consistency
is not a check.

**The hard truncation bound is implemented as a refusal**, not a comment: `|l| <= N/2`, the
grid's own Nyquist, per eq. (32). Nothing aliases down from above what the grid represents,
so the remedy is a finer grid and clamping would answer a question nobody asked.

⚠ **AND THE ACCUMULATION IS BUILT** — `PAC.pnoise(pss, freq, output)`:

    S(f) = Σ_l  h_l CY h_l^H ,   h_l = H_l(f - l f0)

Noise entering at `f - l f0` leaves at `f` through sideband `l`; white sources in disjoint
bands are uncorrelated, so the bands add in power. Each `h_l` is one adjoint row — one
transposed solve for *every source in the circuit*.

| gate | result |
|---|---|
| linear divider vs `analysis_ss.Noise` (Okumura's `p = 1`) | ratio **1.000000**, l≠0 terms ~1e-32 |
| …per doubling of the period grid | 7.80x / 7.33x / 6.75x → 3.9e-09 at 400 pts |
| diode mixer: share of output noise from the fold | **62%** |
| …and the folded total under refinement | 1.0001 / 1.0004 / 1.0001 |

⚠ **BOTH STOPPING RULES, AND WHICH ONE FIRED IS PART OF THE ANSWER.** Ending on the ratio
test means the series converged. Ending on the grid's Nyquist means the grid ran out first
and every sideband above it is **missing rather than small** — a lower bound, and it warns.
Measured on the diode: 80 and 160 points end on the bound, 320 and 640 on the ratio test,
totals differing by 0.04%.

⚠ **STATIONARY SOURCES ONLY, AND IT CHECKS RATHER THAN ASSUMES.** A bias-dependent `CY`
makes the sources cyclostationary; the windows' Fourier coefficients then correlate the
sidebands, they stop adding in power, and the cross terms need the `R_{m,n}` construction
from §III-B — **not read, not built**. Summing powers anyway would answer a different
question and look normal doing it, so `pnoise` samples `CY` along the orbit and raises.
Every source in the discrete element library is bias-independent (a resistor's `4kT/R` never
reads `x`), so the refusal is unreachable there; a compact device's `CY` does read `x`.

**Gate: open, and A1 is built** (2026-09-02), so `H_l` is available. Start from Okumura's
§III, not from the DAC'96 deferral. ⚠ The `p = 1` reduction is the first thing to write:
it exercises the `H_l` path against the stationary formula before any cyclostationary
modelling has to be right.

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
