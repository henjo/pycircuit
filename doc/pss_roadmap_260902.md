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

### A2. PPV — ⚠ BUILT 2026-09-02, normalisation MEASURED not transcribed

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

⚠ **BUILT** as `PSS.ppv()` — two bordered solves, both matrix-free, no eigendecomposition
anywhere:

    [ I - M^T   q ] [v]   [0]          [ I - M   q ] [u]   [0]
    [   q^T     0 ] [y] = [1]          [  q^T    0 ] [y] = [1]

`q = C(0) ẋ(0)` is **exact, not differenced** — the DAE gives `dq/dt = -(i(x)+u)`, so
`q = -(i(x₀) + u(0))`, two evaluations at the converged solution. `y` back at **1.4e-11** is
D&R's free correctness check; null residual **4.7e-11**.

**GATED PHYSICALLY, because an identity check cannot settle a scale.** Displace van der Pol
by `eps·delta`, integrate the *full nonlinear* system until the transverse modes die (second
multiplier 8.6e-4 per period), read the surviving displacement along the tangent. Nothing in
that touches the monodromy, the adjoint or the border:

| npts | worst \|1−ratio\| | rel resid | fitted scale |
|---|---|---|---|
| 200 | 5.52e-02 | 1.90e-02 | 1.006266 |
| 400 | 2.45e-02 | 8.94e-03 | 1.003290 |
| 800 | 1.16e-02 | 4.35e-03 | **1.001677** |

The fitted scale converges to **1** — not to some other constant a direction-only test would
have accepted — and the residual falls at O(h).

⚠ **THE NORMALISATION NOW HAS AN INDEPENDENT STATEMENT FROM THE EQUATION, and it agrees.**
Lai, Zhu & Feng (IMS 2009) give the scaling condition as `⟨Ω₂Q, V₁⟩ = 1`, and `Ω₂Q` is
`dq/dt₂ = C·u₁` — so it is literally `V₁ᵀ C u₁ = 1`. Our bordered solve returns what behaves
as `Cᵀv₁`, so `v·ẋ = 1` **is** that condition. The physical experiment and the equation reach
the same place by different routes.

⚠ **AND THE TRANSCRIBED NORMALISATION IS A 7% ERROR HERE.** Demir's Remark 3.1 reads
`v₁ᵀ C(0) u₁(0) = 1`, and bordering with `q` makes `v·q = 1` fall out for free. But the
vector this bordered solve returns behaves as `Cᵀv₁` — it contracts with a **state**
perturbation directly — so the right normalisation is `v·ẋ(0) = 1`. Both statements are true
of different objects; using one where the other belongs is a silent scale error in every
phase-noise number downstream.

⚠ **The obvious repair is also wrong, and was measured so.** Predicting a state jump's shift
as `vᵀCδ` gives residuals **0.36 / 0.40 / 0.42 that GROW under refinement**, with
per-direction ratios scattering from −0.44 to 28.7. `v·δ` converges at O(h). Do not change it
back without re-running that experiment.

⚠ **AND THE GATE ABOVE IS BLIND TO THE MODEL'S OWN VALIDITY BOUNDARY — measured 2026-09-02.**
The phase equation `α' = v₁ᵀ(t+α) b(t)` treats the oscillator's frequency response as
**instantaneous**; the truth is a convolution, and the PPV form assumes the kernel is
`v₁(t)δ(t−τ)`. Real circuits have finite bandwidth, so a **slow node filters the noise of
devices near it**, the PPV cannot see the filtering, and phase noise comes out
**over-estimated**. Lai (Cadence) is explicit that better extraction does not help: "although
the PPV can be extracted correctly, the oscillator noise analysis is still inaccurate: the
phase noise is always over-estimated."

⚠ **He names this branch's own gate as the blind spot:** the phase equation "was verified to
be correct in many previous works … because it was evaluated on SMALL, SIMPLE OSCILLATORS,
and perturbations were applied to OSCILLATOR CORES." That is van der Pol perturbed at its
core — exactly the gate above. It passes whether or not the failure is present, so it
establishes that the extraction and normalisation are right and says **nothing about the
model's range**.

⚠ **THE EXTRACTION DEGRADES TOO, AND SILENTLY — reproduced here on our own code**, matching a
relayed synthetic result to three digits. The border removes the *phase* mode and nothing
else, so a second multiplier approaching 1 takes the conditioning with it:

| τ/T | \|λ₂\| | σ_min(bordered) | null residual |
|---|---|---|---|
| none | 0.000856 | 8.62e-01 | 4.1e-11 |
| 1e2 | 0.990049 | 4.49e-03 | 4.6e-11 |
| 1e4 | 0.999900 | 4.47e-05 | 4.6e-11 |
| 1e6 | 0.999999 | 4.47e-07 | 4.4e-11 |

`σ_min` tracks `T/τ` over six decades **while the residual does not move**. So `ppv()` now
estimates `|λ₂|` by deflated power iteration (recovered to six digits) and **warns**, saying
the result is an upper bound — because no residual can report this.

⚠ **GATED 2026-09-03 — and the result is a NULL, not a falsification.** Monte Carlo on the
full nonlinear circuit (200 realisations, 150 periods, phase from zero-crossing timing, no PPV
in the measurement): control **0.9965**, slow node at τ/T = 10 **0.8016**.

⚠ **BUT τ/T = 10 IS OUTSIDE THE REGIME THE MECHANISM NEEDS, and an earlier version of this
entry read the null as a refutation.** The reported effect bites through ill-conditioning, and
by our own table `σ_min` at τ/T = 10 is ~**4.5e-02** — healthy. The PPV has no large entries
there and nothing is splitting into nearly cancelling components. Lai's own case is a
gated-capacitor tuning bank whose off-caps have RC "larger than 1 second" at 3.15 GHz, i.e.
**τ/T ~ 3e9** — eight orders from what was tested. **A null at 10 is what the mechanism
predicts.**

**The honest record:** not reproduced at τ/T = 10, which is outside the regime where the
mechanism predicts an effect; **untested** at the τ/T ~ 1e9 where it is reported.

**And the fix is still not built for a reason that stands on its own: COST.** The measurement
needs ~15 time constants of settling, so at τ/T = 1e4 that is 150 000 periods per realisation.
That argument justifies the decision; the falsification framing did not, and this ledger
distinguishes them.

⚠ **The 0.80 is in the OPPOSITE direction to the reported effect.** If it survives the ~10%
Monte Carlo uncertainty it is a separate ~20% *under*-prediction at a τ/T where conditioning
is fine — not a weak version of Lai's. At ~2.5σ it is not established either way.

⚠ **It took three attempts, all the same mistake:** a window of 2–4 time constants read the
slow mode's *decay* as diffusion; an impulse test could not resolve a 1e-11 time shift; and one
noise amplitude for both circuits put a **2.5 V jump per step** on an orbit of amplitude 2.
Each time a number was read before the *measurement* was shown to be in the regime it assumes.
That is a §D shape about the instrument rather than the quantity.

⚠ **And the quadratic isochron rung would not fix it either** — two independent
approximations. Linear-isochron is in the perturbation's **amplitude** (flat hyperplanes;
curvature is what quadratic adds); instantaneous-response is in the **dynamics**. Slow nodes
are the second, and noise is small by construction, so quadratic would earn its cost on *large*
perturbations — injection locking, big interferers — not on phase noise.

**The fix, if it is ever wanted, is a frequency-aware PPV, and it is A1's operator.** FW-PPV is this same bordered
system at nonzero `ω_s`: "if we make the small frequency `ω_s = 0`, it is the augmented PPV
extraction equation … the previous PPV extraction methods give us the transfer function at
`ω_s = 0`." **The classical PPV is the DC point of a PAC-like solve.**

⚠ **BUT THE GAIN IS MODEST AND THE RULE IS NOT "DISTANCE FROM DC" — measured here.**
`I − e^{−jωT}M` is singular at **every harmonic of `f₀`**, not only at DC, and only for an
oscillator:

| offset/f₀ | 0 | 0.25 | 0.5 | 0.75 | 1 | 2 | 3 |
|---|---|---|---|---|---|---|---|
| autonomous | **2.8e-11** | 0.51 | 0.65 | 0.51 | **2.8e-11** | **2.8e-11** | **2.8e-11** |
| driven | 0.79 | — | 0.92 | — | 0.79 | 0.79 | — |

`σ_min` falls **linearly with the distance to the nearest harmonic** (2.5e-1 / 2.6e-2 /
2.6e-3 / 2.6e-4 at 0.9 / 0.99 / 0.999 / 0.9999 of the way). So FW-PPV's advantage is bounded
by `min(T/τ, distance-to-nearest-harmonic)` and is **roughly 1–3 orders in the regime of
interest** — an earlier relayed figure of "nine orders" was withdrawn by its author as a
misreading of their own table.

**PAC now refuses an autonomous solve on a harmonic** rather than discovering it when GMRES
fails — it is physics, not conditioning: a perturbation at a harmonic is a perturbation
*along* the orbit, answered with unbounded phase drift.

⚠ **And FW-PPV as formulated is structurally near-DC.** Mei's Lemma 2.2: augmenting by **one**
column restores rank at `s = 0` but *not* at `s = j·i·ω₀`; Lemma 3.2: the full **Toeplitz
block** is needed for rank everywhere. Lai drops the AC columns for GMRES conditioning, which
reinstates the single-column form and with it near-DC-only validity. Fine for phase noise,
not for a general oscillator AC sweep. (That reading is the docs session's inference, checkable
against Lai's (24) and Mei's Lemma 3.2.) Its foundation (Mei & Roychowdhury, *Oscillator-AC*)
is **now in the library**; Armand 1969 and Adler 1946 are not.

**Why bordered and not an eigenvector:** this is the 2003 improvement's entire content, and
it matters most exactly where a PPV is wanted. Multipliers crowd 1 on high-Q oscillators
(four independent witnesses), `_spectral_report`'s split was measured labelling a parasitic
root physical, and a bordered solve never has to tell candidates apart.

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

⚠ **WHAT IT RETURNS IS THE *TIME-AVERAGED* PSD, now said out loud** (2026-09-02, from
Kundert's tutorial — citation, not measurement). **Two mechanisms make output noise
cyclostationary, and only one is about the sources:** bias-dependent sources modulated by the
operating point (refused, below), *and* the **periodic source-to-output transfer function** —
which applies even when every source is stationary, so a circuit whose only noise is constant
resistors still has cyclostationary output noise. The sideband sum handles the second
correctly and returns its time average.

That is right for most uses and **incomplete for two ordinary RF topologies**, both named: a
**nonlinear subsequent stage** ("an oscillator drives a limiter … the same is true when an
oscillator drives a mixer"), and **cascaded stages off a shared reference**, where "the second
mixer is synchronous with, and tracks the variations in, the cyclostationary noise of the
first." The test is whether anything downstream can track the PSD's variation. Also: a scalar
per frequency cannot carry the **correlation between frequencies separated by `k f0`** that
cyclostationary noise has and stationary noise does not — stated in the docstring rather than
left to be discovered.

⚠ **THE RECORDED BLOCKER WAS WRONG, and the real one is different in kind.** This entry said
cyclostationary needs "the `R_{m,n}` construction … **not read, not built**", as though the
window Fourier coefficients were exotic. They are not. The windows are **rectangular and
non-overlapping**, so `W_{m,k}` is the Fourier series of a **boxcar** — a `sinc` times a phase,
closed form. And the per-interval sources are taken **uncorrelated** (justified because
`H(jω,t)` is time-invariant within each interval), so the sum is *incoherent* over intervals
and coherent only over sidebands within one. Nothing is missing.

**The actual barrier is cost**, which is a different decision to make: the source count is
(timepoints per period) × (noisy devices) — a 500-point grid with 50 noisy devices is **25 000
stationary sources** — and the reported analysis ran at **~14× the PSS per frequency point**
"because all aliasing components need to be computed". Unmeasured here.

⚠ **And its authors left the physics open**: "it is further necessary to discuss the
correspondence between the actual physical phenomena of noises and this modeling". The
windowed-stationary decomposition is a numerical construct whose fidelity to a real device is
not settled by its numerical validation — the same tension as synthesised filters versus
physical trap states in A4d.

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

### A4b. The output layer — NEW 2026-09-02, unbuilt

A1–A3 answer "what does the circuit do"; these are "what should the number be". All from
Kundert's 47-page tutorial — **citations, not measurements**, and weighed as literature.

**AM/PM conversion is a change of basis on results PAC already has.** "It is possible, using
a change of basis, to recast these transfer functions in terms of the AM and PM components of
the modulation." Upper and lower sidebands counter-rotate about the carrier: equal magnitudes
parallel to it is pure AM, perpendicular is pure PM, the general case is an ellipse and is
both. **No new solve** — linear algebra on `adjoint_sideband_row`'s output, and a measurement
designers actually ask for. The cheapest real feature on this list.

**PAC's result is indexed by sideband, not one number per frequency** — "for a single output
frequency there may be many transfer functions from a single input". Conversion gain is one
entry; image, LO and supply rejection are the others. `adjoint_sideband_row` already returns
this shape; what is missing is a result object that carries it. Calibration target: measured
mixers agree "to within 0.25 dB".

⚠ **AM/PM NEEDS THE SIDEBAND CORRELATION, and there is a free invariant.** "To find the AM
or PM noise of a carrier, one must perform PNoise analysis [computing] both the noise at the
upper and lower sidebands … along with [their] correlation." And: **"Linear time-invariant
circuits driven [by stationary] noise sources … can be decomposed into AM and PM noise, but
there will always be EQUAL AMOUNTS OF BOTH."** An LTI circuit must give AM = PM exactly — a
free regression test, and the noise counterpart of the change of basis above.

⚠ **THE ANSWER HAS A HARD BOUNDARY**, and the first version of this entry got the
recommendation **backwards**. It said, from Kundert: `S_phi` is valid at all offset
frequencies while `S_v` and `L(f)` hold only for Δf ≫ `f_Δ` — implying `S_phi` is the safe
thing to report everywhere. **Demir 2006 rejects his own derivation of `S_phi` for exactly
that region:** the excess phase is a Wiener process, so the PSD formula's precondition (a
stable LTI system) fails, and "the PSD … has a singularity at the origin … [it] and its total
power has no mathematical or physical meaning". The obvious patch is rejected too — making
the integrator leaky "end[s] up with a qualitatively incorrect phase noise model", because a
leaky integrator's output is stationary and the phase genuinely random-walks.

**The two sources are not in conflict; they answer different questions.** *Above* the corner
(Kundert's `f_Δ` and the Lorentzian corner are the same frequency) `S_phi` and `L` agree and
both are fine — which is where phase noise is normally plotted, and why the loose usage
survives. *Below* it, `L(f)` is the physically correct object and is finite, but a
small-signal analysis **cannot produce it** (the tapering comes from the nonlinear
phase-to-voltage map); `S_phi` **is** computable there and is meaningless. **So below the
corner neither is correctly available from a small-signal computation** — A4e's near-carrier
branch is what supplies `L` there.

⚠ **Net: do not report `S_phi` near the carrier as though it were a spectrum.** The earlier
framing here implied it was the safe choice everywhere; it is the opposite.

Calibration: phase noise on bipolar resonant oscillators predictable "to within 2 dB".

⚠ **AM/PM IS BUILT 2026-09-03** — `PAC.am_pm(pss, freq, output, carrier)`, with
`am_pm_indices` and `carrier_phasor`. No new solve: the sidebands come from
`adjoint_sideband_row` at `±freq`, and the split is one conjugate —
`m_am = a + conj(b)`, `m_pm = a − conj(b)`.

⚠ **THE CONJUGATE IS THE WHOLE THING, and `a ± b` looks equally plausible.** The sidebands
*counter-rotate* about the carrier, so their sum traces an ellipse; pure AM forces
`a = conj(b)`, pure PM `a = −conj(b)`. Drop the conjugate and the split still returns two
numbers, reporting a rotating ellipse as pure AM. Pinned by a test asserting the naive form
does **not** vanish where the correct one does.

⚠ **AND THE TWO SIDEBANDS ARE NOT CONJUGATES OF EACH OTHER** — that holds for an LTI circuit,
and an LPTV analysis exists because it does not. Conjugating one solve instead of taking two
would force `m_pm = 0` or `m_am = 0` depending on which.

| check | result |
|---|---|
| pure AM → PM index; pure PM → AM index | **exactly 0** |
| diode detector (driven), \|m_pm\|/\|m_am\| | **0.005** — no free phase to modulate |
| …across offsets 0.3 / 0.1 / 0.03 f₀ | flat to 1.2× |
| oscillator PM/AM ratio vs offset | **3.90× / 3.91×** per 4× — the `1/ω_m` divergence |

⚠ **OPEN, AND NOT CLAIMED AS WORKING:** on the *autonomous* circuit the sideband rows come
back at ~1e-12 in absolute terms. The PM/AM *ratio* is right and the driven case is healthy
at O(10³), so the decomposition and the driven path are sound — but why the autonomous rows
are that small is **not established**. Do not read an oscillator AM/PM magnitude from this
until it is.

⚠ **A GMRES ROBUSTNESS FIX CAME OUT OF IT.** SciPy reports breakdown (`info=4`) on these
small operators where it has already solved them — the Krylov space is exhausted in a few
steps and the next vector is numerically zero, a *lucky* breakdown. Trusting the flag turned
exact answers into `RuntimeError` at small offsets. All PAC solves now judge by **residual**;
a genuine failure still fails and quotes it, because near a harmonic the operator really is
near-singular and no tolerance helps.

**Still unbuilt:** the AM/PM *noise* split, which needs the sideband **correlation** rather
than a transfer pair, and carries the free LTI invariant (AM = PM exactly).

**Gate:** none needed for AM/PM — it was a basis change on tested output. The others are
interface decisions, not measurements.

### A4c. Time-varying noise statistics — ⚠ BUILT 2026-09-03 for DRIVEN circuits

The covariance route, and it is a PSS problem the existing machinery already solves. Demir,
Liu & Sangiovanni-Vincentelli (TCAD 1996) propagate a covariance alongside the transient:
`K' = EK + KE^T + FF^T`. For a periodic large signal `E` and `F` are T-periodic, so the
cyclostationary covariance is the T-periodic solution — **a shooting problem whose monodromy
is the Kronecker square of the circuit's**, `M ⊗ M`, with multipliers `λᵢλⱼ`.

⚠ **BUILT** — `PAC.covariance(pss, samples=True)` returns `K(t)` over the period.
`(I − M⊗M)vec(K₀) = vec(K₁)`, **one linear solve, no Newton**, because the Lyapunov equation
is linear in `K` — unlike the trajectory it rides on.

⚠ **GATED AGAINST kT/C — exact, famous, and independent of `R`**, so nothing about the
resistor, the drive or the grid should appear in it, and the covariance shares no machinery
with the closed form it is checked against:

| npts | 100 | 200 | 400 | 800 |
|---|---|---|---|---|
| full `CY` | 1.861 | 1.928 | 1.963 | 1.981 |
| **`CY/2`** | 0.931 | 0.964 | 0.982 | **0.991** |

⚠ **A FACTOR OF TWO SETTLED BY MEASUREMENT.** The two candidate conventions differ by exactly
the factor at issue: full `CY` converges to **2**, halved to **1**. So `CY` is a *one-sided*
density and the per-step injection carries `CY/2`. Either could have been argued from the
definitions; the sequence decides it. The assertion is the **rate** — the error halves per
doubling, O(h), which a piecewise-constant approximation to white noise gives, and which a
wrong *constant* would not do.

⚠ **A PRECONDITION THAT PRODUCED A WRONG READING FIRST.** The initial attempt came back at
**0.517** and looked like a factor-of-two bug. It was not: the RC pole sat at 159 kHz while the
grid's Nyquist was 100 kHz, so the *discrete* system genuinely does not carry the noise the
continuous one does. **A kT/C that comes back low is the grid, not the code** — now its own
test, with the under-resolved case reproduced deliberately.

**Cost:** `(2m)²` unknowns, dense, so `O(m⁴)`. Small circuits only.

⚠ **For an OSCILLATOR it does not exist, and the failure is the physics.** `λ₁ = 1` gives
`λ₁² = 1`, so `I − M ⊗ M` is exactly singular — **verified here on our own monodromy:**
σ_min = 3.1e-11 against σ_min(I − M) = 3.8e-11, spectrum matching `{λᵢλⱼ}` to 2.6e-15, and a
driven circuit showing no such obstruction. The covariance does not settle, it **grows** —
variance linear in `t` is a random walk, which is phase diffusion, which is the linewidth.
**So the near-unit multiplier is not an inconvenience here; it is the answer, and it says the
covariance route is the wrong method for an oscillator.**

**The PPV over the period is now available** (`ppv()['samples']`) — `Phi(T,s)ᵀv(T) = v(s)`,
checked against forward-rebuilt propagators at **1.8e-15** per step. That is what Demir's
diffusion constant `c = (1/T)∫ v₁ᵀBBᵀv₁ dt` integrates.

**The diffusion constant is confirmed from its defining paper.** Demir 2002 eq. (65) is
`c_w = (1/T) ∫ v₁ᵀ B_w B_wᵀ v₁ dt`, and Lemma VII.1 gives `σ²(t) = c_w·t + [coloured]` — the
white term exactly linear in `t`, which is the measured growth (1.825 → 1793.4 over 1000
periods, ratios 9.8/10.0/10.0). **The PPV here, that diffusion measurement, and eq. (65) are
three views of one constant.**

⚠ **TWO DIFFERENT FUNCTIONALS SHARE THE PPV, and one is not a square:**

    c_w  = (1/T) ∫ v₁ᵀ B_w B_wᵀ v₁ dt     WHITE    — QUADRATIC
    V_0m = (1/T) ∫ v₁ᵀ B_cm dt            COLOURED — LINEAR

`V_0m` is the zeroth Fourier coefficient of a periodic scalar. **Using the quadratic form for
a coloured source returns a plausible non-zero number from the same PPV**, and nothing that
did not know to look would catch it.

⚠ **RUN HERE, AND IT IS AN EXACT ZERO.** Demir §VIII: on a parallel-RLC oscillator with a
nonlinear current source — van der Pol is that circuit — "the time-average of [the Floquet
vector entry] for the capacitor voltage is 0! … any … colored-noise source connected across
the capacitor has NO contribution to the oscillator spectrum due to phase noise." Measured on
our PPV: **|mean|/rms ~ 1e-11** with rms ~0.40, at μ = 0.5 and 1.0. Zero and non-zero from one
vector, which is exactly the discrimination the coloured functional needs and the quadratic
one destroys. Now a test.

### A4e. Oscillator phase noise — ⚠ CLOSED FORM BUILT 2026-09-03

⚠ **BUILT** — `PAC.diffusion_constant(pss)`, `PAC.lorentzian(...)`,
`PAC.oscillator_spectrum(pss, offsets, output, harmonic)`.

| check | result |
|---|---|
| `c` vs the A2 gate's **Monte Carlo** measurement | 1.5903e-07 vs **1.5959e-07** (0.35%) |
| `∫ S_i df` over the implemented function | **1.000000** for every harmonic and `c` |
| peak, half-width | analytic to 1e-12 |
| far skirt | **4.00× per doubling** = 1/f² |
| harmonic `i` far out | **20·log₁₀(i)** dB, to 0.05 dB |

⚠ **The one input is already tied to a physical measurement.** `c` is the same integral the
A2 gate measured a completely different way — Monte Carlo on the full nonlinear circuit, phase
from zero-crossing timing, no PPV anywhere in it. So the closed form does not rest on a fresh
claim.

⚠ **Total power conservation is the invariant that matters**, and it is asserted by integrating
the *implemented* function rather than re-deriving the algebra: noise spreads the carrier's
power into a line of finite width and creates none. LTV treatments "erroneously predict
infinite noise power … as well as infinite total integrated power", so this is the property
that says the closed form is doing the nonlinear thing.

**White sources only, where the lineshape is EXACT rather than a two-regime limit** — which is
also all `pnoise` supports.

⚠ **NOT A SWEEP, AND NOT pnoise's SHAPE.** Demir 2002 (68)+(69): once the PSS waveform's
Fourier coefficients `X_i`, the single scalar `c_w`, and one `V_0m` per coloured source are
known, "we have an analytical expression that gives us the spectrum at any frequency. The
computation of the spectrum is not performed separately for every frequency of interest."
Closed form in a handful of scalars — no per-frequency solve, no sideband grid, and it
sidesteps the 1/f sweep-grid trap in A4d because there is no sweep to land a point on.

Near carrier a Lorentzian with `c_eff = c_w + Σ|V_0m|² S_Nm(0)`; away from it white gives
`1/f²` and coloured gives `1/f³`. The i-th harmonic's skirt scales `i²` and its corner `i⁴`,
so higher harmonics are `20·log₁₀(i)` dB noisier.

⚠ **(69) IS A TWO-REGIME LIMITING FORM, NOT EXACT** — the transform "does not have a simple
closed form", and there is no exact expression joining the regimes. Whatever is implemented
must say which regime it is reporting. That is the quantitative version of A4b's `f_Δ`
boundary.

⚠ **AND 1/f NEEDS AN EXPOSED CUTOFF.** `K I^a/f` "is not a well-defined spectral density …
It blows up at [f=0]", so Demir *postulates* a cutoff. The near-carrier linewidth depends on
`S_Nm(0)`; **with no cutoff the linewidth diverges**. It is directly visible in the output, so
it must be a parameter and not a buried constant.

**Two assertable invariants:** the carrier PSD is finite and per-harmonic total power is
conserved — analyses "based on linear time-invariant or linear time-varying concepts
erroneously predict infinite noise power … as well as infinite total integrated power".

**Validity (Assumption IV.1):** the coloured theory needs source bandwidth ≪ `f₀`. True of 1/f
and burst noise; **a coloured source with bandwidth comparable to `f₀` is outside it.**

⚠ **JITTER IS WELL-POSED EXACTLY WHERE THE SPECTRUM IS NOT.** "The integral of a stationary
process is itself not necessarily stationary, but it has stationary increments … The
difference operation, in a sense, undoes the nonstationarity." The increment
`α(t+Δt) − α(t)` is the output of a *stable* delayed integrator, so its PSD is legitimate and
"does not explode at f=0". The reported quantity is therefore **σ²(Δt), a function of
accumulation time** — not a bare spectrum and not a single number.

⚠ **The common literature formula is built on the bad object.** Equations relating jitter to
phase noise "appear elsewhere … [where `S_phi`] is defined in the sense of (23) as the
ill-defined PSD"; in others "no distinction is made between the well-defined `L(f)` and the
problematic `S_phi`". **If jitter-from-spectrum is implemented from a typical reference, it
will be the unsound version.** The sound one integrates the macro-source PSD through the
stable delayed-integrator transfer function.

**A cheap sanity check on any implementation:** white gives `σ² ∝ Δt` (jitter ∝ √Δt); 1/f
gives `σ² ∝ Δt²` (jitter **linear** in Δt), because "within small accumulation times, 1/f
noise samples in time are indeed almost fully correlated".

⚠ **TWO LOW-FREQUENCY CUTOFFS AND ONLY ONE IS LEGITIMATE.** A cutoff on the **1/f source
model**, representing finite observation time, is sound — Demir is explicit it is "*not* …
the 3-dB frequency of the Lorentzian" and that "we are NOT trying to fix a problem we created
by doing something ill-defined". A cutoff on `S_phi` to stop it diverging is the illegitimate
one. Because the legitimate cutoff encodes observation time, **a 1/f jitter number without a
stated observation window is incomplete.**

⚠ **EXPOSE `c` AND `c_i` AS FIRST-CLASS OUTPUTS.** Demir 2006 §IX exists *only* as a
workaround for simulators that do not provide them — reverse-engineering the scalars by
curve-fitting `L(f)`, or by switching sources on and off. **The PPV computes them directly.**
They are also the per-source attribution handle: the total jitter decomposes into "contributions
from the individual noise sources, e.g. the on-chip inductor and the transistors" — available
in simulation and not from measurement.

**Naming, on the author's own warning** that the historical terms are "somewhat confusing and
not very precise": call them `ssb_phase_noise_L` and `jitter_variance(delta_t)`, not
`phase_noise` and `jitter`. Three of this branch's defects have been unstated conventions;
this is cheap insurance.

### A4d. Coloured noise — NEW; costs nothing structural for the path we built

⚠ **THE FILTER DOES NOT ENTER THE PSS.** Demir 1996 synthesises 1/f from white sources
through a Lorentzian network because Itô theory admits only white driving noise — "we can not
express a flicker noise source in terms of the standard white Gaussian noise process" — at
"one state variable per decade of frequency". That reads like it would drag 10+ decades of
time constants into the shooting run. **It does not.** The network has no deterministic drive
and is one-directionally coupled, so the augmented `E` is block-triangular and the filter's
deterministic periodic solution is identically zero. Measured: PSS solution identical with and
without the filter to **1.2e-14**, and `M_aug`'s circuit block matching `M` to 2e-15…6e-15 at
every `tau` tried.

⚠ **AND IN A FREQUENCY-DOMAIN pnoise IT DOES NOT EXIST AT ALL** — which is the path A3 took.
The Lorentzian network is an artefact of the SDE formulation specifically. No SDE is formed
here, so nothing needs synthesising: **a coloured source is just a different `S(f)`, folded
like any other.** SpectreRF confirms by omission (no filter, no augmentation in its Pnoise
treatment of flicker), and Kundert states the consequence directly — "S_u(f) is generally pink
or proportional to 1/f. Then S_phi(f) would be proportional to 1/f³ at low frequencies." **A
slope, not a state.**

⚠ **AN EARLIER VERSION OF THIS ENTRY CLAIMED AN ARCHITECTURAL FORK — filter states inside the
MNA versus outside — and there is no fork.** It was recorded as untested and it was wrong.
The reasoning that produced it was sound in isolation (Demir's "large time constants", Kundert
on warm start, A4's own μ=0.05 → ~24 periods) and reached a false conclusion because it never
asked whether those states are *coupled*. Shape 4 from §D: two things each considered alone.

⚠ **AND "COSTS NOTHING STRUCTURAL" NEEDS A DISTINCTION THIS ENTRY DID NOT DRAW.** It is true
of a **synthesised** Lorentzian filter — fictitious, linear, no deterministic drive, no
dependence on the circuit, measured to stay out of the PSS. It is **not** true of a
**physically modelled trap**: real oxide-trap states are bias-modulated and depend on the
circuit solution down to the drain voltage, so they have a deterministic periodic component
and genuinely belong among the circuit variables. Same words, two different objects — the
recurring failure shape on this branch.

**What survives, and only for A4c's route:** the covariance system's spectrum is `{λᵢλⱼ}`, so
slow filter poles enter **squared**. Measured `σ_min(I − M_aug ⊗ M_aug)` tracking `T/tau`
exactly: 8.6e-01 / 1.4e-01 / 1.4e-02 / 1.4e-04 / 1.4e-06 at `tau/T` of 1e0…1e6. A 10-decade
network sits near 1e-8. ⚠ **That ill-conditioning is physical, not numerical:** a process with
a 10⁶-period time constant has not reached periodic steady state within one period, so its
periodic covariance is not a well-posed question. The A4c one-solve trick will not swallow 1/f,
and correctly so.

⚠ **TWO OPPOSITE SWEEP HAZARDS, AND A GRID POLICY HAS TO SATISFY BOTH.** Landing a point
*on* a harmonic gives absurd totals (below); sampling *too coarsely* near the injection
multiples of a driven oscillator hides real structure — a 5 GHz LC oscillator with a 5 MHz
supply perturbation shows spikes at 5 and 10 MHz from up-converted flicker that a 5-point/decade
log sweep misses entirely. Same physics (1/f power at DC, translated by a periodic signal),
opposite sampling failures.

⚠ **THE REAL 1/f ITEM FOR THIS PATH IS A SWEEP-GRID TRAP, NOT AN ARCHITECTURE.** A 1/f source
is singular at DC, and folding puts a copy of that singularity at **every harmonic**.
SpectreRF: "place a cluster of frequencies near each harmonic … but AVOID PUTTING FREQUENCY
POINTS PRECISELY ON THE HARMONICS … you run the risk of generating absurd noise totals because
a very narrow noise peak artificially has its apparent width greatly magnified … and has its
amplitude exaggerated by placing a point precisely at the singularity." Plausible-looking
nonsense, no error raised — the same shape as the convention defects this branch keeps
finding. Not reachable today (every source in the discrete library is white), and it belongs
in `pnoise`'s preconditions before the first coloured source arrives.

⚠ The legacy 1/f-under-modulation rule (use the DC average of the time-varying current) is a
**placeholder its own author disowned** in the same paragraph: "either a theoretical or
experimental derivation of a model for flicker noise associated with a time-varying current is
needed." Do not implement it as if it were the model.

⚠ **AND THE PAPER THAT CLOSES IT SAYS SWITCHING WHITENS 1/f, NOT SCALES IT.** Mahmutoglu &
Demir 2015: under a switched bias the trap state "is perfectly reset in every switching
period … the trap noise turns into white noise for time scales above the switching period …
a WHITENING EFFECT, as opposed to a frequency independent power reduction." **Below
`f_switch` the spectrum is flat, and the corner moves to `f_switch`, not the trap rate.**
Above it the benefit reverses — the switched PSD "falls on top of the PSD for the
non-switching case", so a flat "3 dB better because switched" model is wrong in *both*
directions.

**Consequence for a stationary flicker model at an averaged bias**, which is what would be
implemented first: below `f_switch` the real spectrum is white while the model says 1/f, so
**the error grows without bound as the offset decreases**. Not a constant offset, not
trimmable with a fudge factor. If 1/f is ever added that way, say so explicitly — the same
class of fix as pnoise's two named exceptions.

⚠ **The idealised correction is an upper bound, not a prediction.** Its own §II-F: the
assumption that a trap empties instantly at switch-off "is clearly violated in physical
systems", and with a realistic trap model the whitening is only *partial* — "it is quite
likely that a full trap stays full through the off-state". Truth sits between ideal switching
(fully white) and always-on (full 1/f). Shipping the idealised form trades over-prediction for
under-prediction, which is worth saying rather than discovering.

### A6. Driven oscillators and PLLs — REQUESTED 2026-09-03

Directive: *"we will need driven oscillators for pll analysis"*. Recorded now while the
reading is fresh; nothing built.

⚠ **THE T4-1 FALSIFICATION STANDS AND DOES NOT TRANSFER.** C5 closed saltation correction
after measuring the monodromy–FD gap falling at exactly 2.00x per doubling — O(h)
discretisation, not the O(1) a missing term leaves. **That result holds for pycircuit's switch
models**, whose discrete map has no instant where the field is undefined. A PLL with a
**digital PFD** is the other side of that line: the model is not Lipschitz continuous, and
there saltation is **mandatory and carries the feedback**.

| | field at the switch | saltation |
|---|---|---|
| pycircuit switch models | **defined** | not needed — C5, falsified correctly |
| Verilog PFD / integer divider | **undefined** | **mandatory** |

⚠ **AND THE REASON IS STRUCTURAL, NOT ACCURACY.** Along a *locked* orbit the reference and
divided edges coincide, the PFD never changes state, and "the analog part of the circuit
behaves as if it were not connected to the digital one; the loop … is thus **OPEN**." The
plain variational model then has a null block and **at least two unit multipliers** —
describing an open loop while the circuit runs a closed one. Omitting saltation there does not
give a slightly wrong Jacobian; it gives the Jacobian of a different circuit.

⚠ **A GATE WOULD MISFIRE HERE.** "At least two unit eigenvalues" is exactly the condition for
lacking asymptotic phase, so running that check on a naive locked-PLL variational model flags
it invalid — and the right response is to **add the saltation terms**, not to distrust the
gate. Our own `ppv()` second-multiplier warning would fire for the same reason.

⚠ **THE 1/f OBSERVATION-WINDOW CAVEAT (A4e) DOES NOT APPLY TO A PLL.** Free-running: the
jitter integral diverges and needs a cutoff tied to observation time. In a PLL "the integral
in (61) **converges** even when the ideal flicker noise PSD is used", because the loop's
jitter transfer is **high-pass** — it suppresses exactly the low-frequency content that makes
the free-running integral diverge. So free-running needs the window; a PLL does not, and a
single total rms jitter in ps is a meaningful number there.

**Different output shape too:** free-running variance grows without bound (∝Δt white, ∝Δt²
for 1/f); a PLL's **settles**, and can peak and ring on the way from the locking dynamics.
**Free test:** at short Δt the loop has almost no effect, so PLL jitter must equal
free-running VCO jitter there.

⚠ **THREE JITTER METRICS, NOT INTERCHANGEABLE** — absolute/edge, k-cycle
`τ(t+kT) − τ(t)`, and cycle-to-cycle (the *second* difference). Three numbers from one phase
process. Given three of this branch's defects have been unstated conventions, `jitter` as an
API name without the metric attached is the same trap.

**Noise folding in a PLL is ordinary LPTV aliasing**, present in an all-analog PLL too — not
an artefact of the PFD sampling. No new machinery to explain it.

**Two modelling choices worth copying:** use **zero crossings** rather than edges as switching
events (saltation needs `grad h`, so the surface must be differentiable), and **promote the
reference phase to a state variable**. That last is the fourth appearance of one device in this
line — promote the quantity you want to perturb into an explicit state, or normalise time by
the local frequency, so periodic objects become comparable. Worth recognising up front.

⚠ **AND THE PPV IS AN OPERATING-POINT OBJECT.** For A2's use — noise at a locked PLL's
operating point — one PPV is correct, and sweeping the control voltage means redoing the
analysis at each point, which is ordinary. But PPVs at different control voltages have
**different periods**, so a tabulated PPV cannot be interpolated pointwise; normalise by the
instantaneous free-running frequency first. (The "totally wrong predictions" result in that
literature is about using a fixed PPV as a *transient macromodel* across a capture transient —
a different activity, and not implied by anything in A1–A3.)

⚠ **A BASE-SOLUTION TRAP, AND IT IS NOW GUARDED.** The natural way to model a driven
oscillator — solve the PSS of the bare oscillator, treat the injection as a perturbation — is
**wrong**, because the injection *device* is present even when its *signal* is zero: "in
absence of the injection signal, the injection circuit affects the basic LC oscillator by
**changing the nonlinearity of the feedback loop** … [it] can affect the start-up condition …
**or its oscillation amplitude**, or both." So the free-running orbit of the
circuit-with-the-device is not the orbit of the circuit-without-it, and every Floquet quantity
built on the wrong one inherits the error. **The analysis converges and reports a plausible
number.**

It generalises past dividers to **any driven oscillator whose drive enters through a real
device** — injection-locked VCOs, supply and substrate coupling paths. `PAC` now refuses an
operating point solved on a different circuit object; the existing reference-node check could
not catch it, since two circuits differing by one device have the same reference node and
often the same node count.

**An exact test vector exists for the unlocked driven oscillator** (Armand 1969, on Adler
1946): the whole spectrum from one scalar `K`, sidebands at `w1 + n·Ω` on **one side only**,
magnitudes a geometric ladder with ratio `tan(θ/2)`, phase advancing by exactly `θ` per
sideband, `Σ|A_n|² = E²` **exactly**, and clean degeneration to a single tone at both limits.
A conserved quantity plus two limits catches normalisation, sign and branch errors — the
failure class that has cost this branch the most. Valid within Adler's three conditions, whose
first two are **high Q** and **slow time constants** — the same two triggers as A2's validity
boundary, sixty years earlier.

⚠ **And Adler's equation has limits on what it can PRODUCE, not just where it holds:** phase
only (no envelope), small injection only, first harmonic only. An ILFD locking range is
available in closed form beyond those limits — **inversely proportional to Q**, SPICE-validated
over R = 300–800 Ω — with the authors flagging their own artefact: the first-order expression
is **symmetric about tank resonance** and real locking ranges are not, repaired by substituting
the free-running frequency for `ω₀`.

### A7. Near-carrier oscillator noise — ⚠ BUILT 2026-09-03

⚠ **`Φ(T) − I` IS SINGULAR FOR AN OSCILLATOR, AND ITS NULL VECTOR IS THE PPV.** So a
near-carrier noise computation is ill-conditioned *by construction* — the thing being computed
is the thing that breaks the matrix.

⚠ **THE OBSERVABLE SYMPTOM, worth recording before it is seen:** "the standard time domain
noise analysis yields **flat PSD curves or curves with unexpected slope near the oscillation
frequency**." Flat oscillator noise near the carrier is *that singularity* — not the physics,
the noise models or the source definitions. It points at the right layer immediately.

⚠ **BUILT — and derived here rather than transcribed**, which the ledger argues for. With `u`,
`v` the right and left null vectors of `I − M` (the tangent and the PPV, both of which `ppv()`
already returns):

    [ I − αM   u ] [ w ]   [ b ]
    [   vᵀ     0 ] [ s ] = [ 0 ]

Since `vᵀ(I − αM) = (1 − α)vᵀ` and `vᵀw = 0`, the border variable is **bounded**,
`s = (vᵀb)/(vᵀu)`, with no `1/ε` in it; and since `(I − αM)u = (1 − α)u`, the solution is
`y = w + s·u/(1 − α)` — the vanishing factor in **closed form** rather than inverted
numerically.

| offset/f₀ | σ_min plain | σ_min **bordered** |
|---|---|---|
| 3e-01 | 5.68e-01 | 1.17e-01 |
| 1e-03 | 2.61e-03 | **2.04e-01** |
| 1e-06 | 2.61e-06 | **2.04e-01** |
| 1e-09 | 2.61e-09 | **2.04e-01** |

The plain operator tracks the offset over **nine decades**; the bordered one is flat. The two
agree to **5.7e-12** where the plain solve is still trustworthy, and their disagreement grows
as `1/Δf` — that is the *plain* solve losing digits.

⚠ **THE TEST ASSERTS BOTH HALVES.** Flat conditioning alone would be satisfied by an operator
that had stopped solving the right problem, so agreement with the plain solve *where the plain
solve can be believed* is what says it is still the same equation. And `|y|` growing 10× per
decade closer to the harmonic says the pole was removed from the **conditioning**, not from the
**answer**.

Wired into `adjoint_sideband_row` for autonomous circuits only — a driven circuit has no pole,
and the plain solve is correct and cheaper there. `HARMONIC_GUARD` drops from 1e-6 to
**1e-12**: it used to be the conditioning floor being accepted and now excludes only what has
no finite answer. Transposed variant included, borders swapped.

⚠ **It still refuses an EXACT harmonic**, and should: `1/(1 − α)` is then a division by zero
and the physical response is unbounded.

⚠ **A scope split for the PLL work:** integer-N is tractable by PSS methods; **fractional-N
with a Δ-Σ divider is outside both shooting and HB** — its period is "so large that the
shooting or the harmonic balance methods are inapplicable" — and must go time-domain, where a
measured numerical noise floor in commercial simulators becomes the binding constraint. Worth
saying before fractional-N is promised.

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

0. **A result asserted outside the regime the claim is about.** A4's criterion tested on an
   autonomous circuit when the paper says non-autonomous; the A2 gate's window set to 2–4 time
   constants when the quantity is asymptotic; and the A2 gate's *conclusion* read as a
   falsification when it was run at τ/T = 10 and the mechanism needs τ/T ~ 1e9. **Three
   instances, and the third was in the conclusion rather than the instrument.** Before reading
   a null as evidence, check that the measurement was inside the regime where the effect is
   predicted to exist.

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
