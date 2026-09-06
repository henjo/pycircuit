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

## 0. The organising fact — read this before any item below

⚠⚠⚠ **THE THING DESIGNERS OPTIMISE FOR IS THE THING THAT BREAKS EVERY NUMERICAL METHOD, AND
THIS IS AN IDENTITY RATHER THAN A TRADE-OFF ANYONE CHOSE.**

A designer's whole objective on an LC oscillator is to push `Q` up. Pushing `Q` up *is* pushing
`|λ₂|` toward 1. And `λ₂ → 1` is simultaneously every failure mode in this document. **The
better the circuit, the worse the numerics.**

The link is not an analogy. Wang & Roychowdhury (2017): `Q = log(threshold)/log|λ₂|` — so
"bounded by `Q`" and "bounded by `λ₂`" are **the same sentence**. ⚠ **And we already ship it**:
`PSS.ppv()`'s `info['Q']` is `-1/log(λ₂)`, which is that identity at a `1/e` threshold
(`shooting.py:2558`). A high-Q oscillator is *defined* by slow amplitude decay, and slow
amplitude decay *is* the second multiplier approaching one.

⚠⚠⚠ **AND `Q` IS NOT THE ONLY AXIS — TWO BOUNDS ON THIS SECTION ITSELF, ADDED 2026-09-04.**
Everything below is organised on `Q ↔ λ₂`. That organisation is sound for what it covers and is
**not complete**, and a reader who takes it as complete will look for every problem along one
axis. Both bounds are *cited, not measured here* — read from Traversa & Bonani TCAS-I 2011 p. 5 by
the docs session, whose equations `pdftotext` drops, so they were read as an image.

  **(i) A SECOND, ORTHOGONAL AXIS ON WHICH THE FRAMEWORK ITSELF DEGRADES.** The phase
  decomposition assumes `|α̇(t)| ≪ 1` — *"the decomposition proposed here, **as well as the theory
  in ([1], p. 661)**, is based on the assumption that |α̇(t)| ≪ 1"*, where [1] is **Demir, Mehrotra
  & Roychowdhury**. So the assumption sits under *everything* we compute — the PPV, `c`, the
  Lorentzian — and it is stated by a third party about Demir's paper rather than by Demir.
  ⚠ **IT IS INDEXED BY AM-TO-PM CONVERSION, NOT BY `Q`.** Their controlling group is `v·ε`
  (AM-to-PM coupling × noise amplitude), checked against an exact Fokker–Planck solution: exact at
  `v = 0`, *"still good"* at `v·ε = 0.1265`, *"less favourable"* at `0.3162`, clearly degraded by
  `0.9487`. **A circuit can be low-`Q` and still break this**, which no amount of attention to
  `λ₂` would predict. ⚠ **AND IT IS A BOUNDARY ON OUR METHOD, NOT ON OSCILLATOR NOISE** — Bonnin
  2015, *"Amplitude and Phase Dynamics of Noisy Oscillators"*, gives amplitude and phase SDEs whose
  validity is *"not limited to the weak noise limit"*, with closed forms for the expected angular
  frequency, amplitude and amplitude variance (cited, not verified here). So `|α̇| ≪ 1` bounds the
  Floquet/PPV family; it is not a hard physical wall, and should not be recorded as one.

  **(ii) `λ₂` DOES NOT ORDER THE ORBITAL CONTRIBUTIONS.** From the same paper's HBT results: six
  orders of magnitude separate two Floquet exponents, and the corresponding orbital-noise
  contributions are **not** in that ratio — *"far from the oscillator harmonics, the contribution
  of [the smaller one] is dominant … this clearly shows that **also the eigenvectors may give an
  important contribution** to the orbital noise spectrum, which might also dominate over the
  [exponent] factor"*. ⚠ **So for the FAR-OUT spectrum a `λ₂`-only view is insufficient by
  construction** (A9): eigenvector magnitudes can outrank exponents. This does not weaken any
  near-carrier item below — those really are one edge — but it does mean **`λ₂` is not a
  sufficient statistic for the far skirt.**

⚠ **So the seven places this quantity appears are not seven findings. They are one edge seen
from seven sides**, and they enter through structurally different doors:

| | door | source |
|---|---|---|
| (a) | **numerical distinguishability** — `λ₁ = 1` cannot be separated from `λ₂ ≈ 1`, so eigen-based PPV extraction fails | **Demir & Sangiovanni-Vincentelli 1998, Table 6.4 — a REPORTED failure**, the stated motivation for Demir, Long & Roychowdhury 2000; 2003 is the fuller procedure |
| (b) | **conditioning** — the bordered Jacobian degrades; `σ_min` tracks `T/τ` over six decades | Lai DAC 2006; measured here |
| (c) | **truncation validity** — the single-mode reduction needs `\|exp(η_i)\| ≪ 1` | Demir 1998 (6.72) — ⚠ **the SAME condition as (a), seen from the other side**: the book says (6.72) holds for "most" oscillators and defers the exceptions to its high-Q crowding discussion, which is (a). One condition: an approximation that stops holding, and an algorithm that stops working. ⚠⚠ **AND IT IS THE SAME NUMBER, NOT JUST THE SAME CONDITION**: `η_i` is the exponent times the period, so `exp(η_i)` **is** `λ_i` and (6.72) reads `\|λ₂\| ≪ 1` — the settling `Q` itself. Reading "≪ 1" as 0.05 puts the bound at **`Q_λ < 0.334`**. Our ordinary gates (`Q_λ = 0.14`) sit inside it; **the entire high-Q programme does not** — `\|λ₂\| = 0.730` at `Q_λ = 3.18`, `0.9845` at 64. So the closed-form single-mode variance has **no justification anywhere in the high-Q sweep**. ✅ **RESOLVED BY A PAGE READ (docs session), AND THE ANSWER IS THE OPPOSITE OF THE RECONSTRUCTION THAT WAS REFUSED.** In the book, (6.72) is invoked at exactly one place: to drop the `i = 2..n` terms of the Floquet-basis expansion (6.71) inside the integral (6.69) that gives the `k`-th diagonal entry of `K(t)` — the time-varying variance of a **state variable**. *"When (6.72) is satisfied, the contribution of the terms for i = 2,…,n … will be negligible."* It is **not** invoked for the spectrum, and **not** for the reduction to the phase equation; `c` is *defined* by the phase projection, not approximated by it, and the Lorentzian follows rigorously. My reconstruction — that it attaches to `c` and propagates to the Lorentzian and `oscillator_spectrum` — was **wrong**, and refusing to write it was correct. ⚠ **CHECKED IN THE TREE: NO SHIPPED SURFACE MAKES THE (6.72) APPROXIMATION.** The only state-covariance routines are `covariance` and `oscillator_covariance`, both via `_lyapunov_pieces` on the **full** `I − M⊗M` (`shooting.py:7186`) — an exact bordered Lyapunov solve with no modal truncation. `oscillator_covariance` *keeps* the transverse modes: `K = K_orb + c·t·uuᵀ`, and `orbital_mode_weights` measured the phase-mode weight of `K_orb` at `1.6e-19` with the amplitude mode carrying it. No routine computes a state or node variance from `c` alone; the only `c·T` is `d = cT`, the *phase* variance, which is the rigorous term. **So the honest row is: the book's derivation makes an approximation at high Q that this implementation does not make.** The high-Q programme is outside the *book's* validity for that one equation, and inside ours |
| (d) | **theory validity** — two multipliers at 1 means no asymptotic phase; the PPV is undefined | Demir 2006 |
| (e) | **settling and ringing** — long `tstab`, ringing impulse response, large `M` | SpectreRF; Hull & Meyer; the probe methods |
| (f) | **method-dependence of the value** — backward Euler biases `λ₂` low, so `Q` is method-dependent | measured, docs session |
| (g) | **physical identity** — `Q = log(threshold)/log\|λ₂\|` | Wang & Roychowdhury 2017 |

**(g) is why the other six are not a coincidence:** each of (a)–(f) is a sensitivity to slow
amplitude decay, and (g) says that is precisely what `λ₂` measures.

⚠⚠⚠ **AND (g)'s NAME IS ONLY RIGHT WHILE THE OSCILLATOR'S AMPLITUDE MODE IS THE SLOWEST
NON-UNIT MODE.** A parasitic with `τ_p/T > Q_osc` simply **is** the second multiplier — by
definition, not by error — and then `Q = −1/log|λ₂|` reports **the parasitic's decay time in
periods** under the name `Q`. MEASURED on a `Q = 16` oscillator: at `τ_p/T` = 32 and 100, `λ₂`
is the parasitic and `Q` returns 32 and 100.

⚠ **The number is right and its NAME is wrong**, which is why nothing misbehaves — every
residual stays clean and the value is well converged. **A DCO's gated capacitor sits at
`τ_p/T ≈ 10⁴`, permanently in that regime**, so on exactly the circuits a hierarchical DCO
method exists for, a reported `Q` would be the gated cap's RC in periods. That connects the DCO
thread to the `Q` identity. Read `Q` as *"cycles for the **slowest non-unit mode** to decay by
1/e"*, which is what it computes; it is the *oscillator's* `Q` only when that mode is the
oscillator's. Everywhere else in this record the two have been used interchangeably.

⚠ **AND THE FIELD HAS NOT ABSORBED (g).** A citation-forward search on the anchors:
`Q = log(threshold)/log|λ₂|` has been **cited once, ever**. The shooting-methods literature has
essentially stopped — **9 citations of Telichevesky 1995 since 2021**, most from other fields —
while the PPV line is cited ~150 times, almost entirely by applications. So the identity that
organises this whole document is one we re-derived independently and ship at
`shooting.py:2558`, and it is **not** a result anyone is building on.

⚠ **Read the count carefully, though: a field can stop publishing on something because it is
FINISHED rather than because it is abandoned.** Nine citations since 2021 alone reads as
abandonment; the same nine *plus* Cadence shipping shooting-adjacent PSS acceleration into
**PSpice** in 2023 reads as settled — the vendor datum is the discriminator. And "cited once
ever" could mean nobody noticed *or* that everyone who needed it re-derived it. **We re-derived
it independently**, which is one data point for the second reading. Treat that as a reason to
state it carefully rather than as a reason to doubt it: we verified it numerically here
(`test_the_reported_Q_amplifies_its_own_lambda2_error_by_Q`) and it has now reproduced in three
independent settings.

⚠ **AND IT PREDICTS WHERE THE NEXT ONE WILL BE: anywhere a method assumes the non-oscillatory
modes have died.** That is a search rule, not a summary. It already found one — the
cyclostationary cost blocker's *remedy* is bounded by the same quantity (Hull & Meyer's cheap
construction fails because the impulse response rings), so even the escape route from a cost
problem is door (e).

⚠⚠⚠ **AND (g) IS NOT A SEVENTH DOOR — IT AMPLIFIES ALL THE OTHERS.** Differentiating
`Q = −1/ln λ₂`:

```
(dQ/Q) / (dλ₂/λ₂)  =  −1/ln λ₂  =  Q
```

**The relative error in `Q` is `Q` times the relative error in `λ₂`.** MEASURED end to end in
this solver, not just in the algebra, on van der Pol tuned by the recipe below:

| `Q` | npts | rel err `λ₂` | rel err `Q` | **ratio** |
|---|---|---|---|---|
| 3.18 | 120 | 1.04e-03 | 3.31e-03 | **3.2** |
| 15.92 | 120 | 5.75e-04 | 9.23e-03 | **16.1** |
| 63.66 | 120 | 4.89e-04 | 3.21e-02 | **65.7** |
| 63.66 | 480 | 7.56e-06 | 4.81e-04 | **63.7** |

⚠ **So the exposure is a RESOLUTION REQUIREMENT SCALING WITH `Q`, not a fixed accuracy.** To
report `Q` to 1% needs `λ₂` to 1e-4 relative at `Q = 100` and 1e-5 at `Q = 1000`. **Payable
here** — Gear-2's `λ₂` converges at better than second order (~8× per doubling) — but a method
that *biases* `λ₂` at fixed order has no escape: backward Euler's measured **−5.6e-2** bias
becomes **−85% in `Q` at `Q = 100`**. Now pinned by
`test_the_reported_Q_amplifies_its_own_lambda2_error_by_Q`, which asserts the *sensitivity*
rather than the accuracy.

### The high-Q fixture recipe — `μ = 1/(2πQ)`

`λ₂ ≈ exp(−μT)` with `T ≈ 2π` gives `Q = 1/(2πμ)`. **Verified on this solver to four digits:**
predicted 3.183 / 7.958 / 15.92 / 63.66 against measured 3.182 / 7.959 / 15.92 / 63.67.
Available as `_vdp_at_Q(Q)`.

⚠ **Two costs.** The period seed must be `2π/√(1−μ²/4)` with `reltol = 1e-12`, or the shooting
solve becomes the thing under test. And *transient* settling to 1% takes **≈ 5.4·Q periods** —
door (e) taxing the **test harness** rather than the circuit.

⚠⚠ **A "5.4·Q" CORRECTION WAS PROPOSED, REPRODUCED HERE, AND THEN RETRACTED BY BOTH SESSIONS.
4.6·Q IS RIGHT.** The claim was a transient-growth prefactor `C ≈ 2.1–2.3` in `‖M^k‖ ~ C λ₂^k`,
making 1% settling `Q(ln 100 + ln C)`. Two sessions measured `C` at 2.101 / 2.201 / 2.314 and
both read it as a property of the circuits. **It is a property of the deflation.**

| fixture | `λ₂` | `‖PMP‖` | argmax_k | monotone | `‖PMP‖/λ₂` | `‖P‖` |
|---|---|---|---|---|---|---|
| vdP Q=16 | 0.939426 | 2.1013 | **1** | yes | 2.2368 | **2.236** |
| vdP Q=64 | 0.984505 | 2.2014 | **1** | yes | 2.2360 | **2.236** |
| bulk m=12 | 0.983536 | 2.3143 | **1** | yes | 2.3530 | **3.050** |

⚠ **There is no hump.** `argmax_k = 1` and the sequence is monotone from `k = 1`, so the "peak"
is just `‖PMP‖` — and it equals `λ₂·‖P‖` to four digits. The spectral projector used to remove
the phase mode is **oblique**, and its norm is the whole prefactor: `‖PMP^k‖ = λ₂^k‖P‖`
identically.

**So the measured 1.18 ratio is exactly `ln‖P‖`.** A criterion `‖PAP‖ < 0.01` needs
`λ₂^k < 0.01/‖P‖`, i.e. `k = Q(ln 100 + ln 2.236) = 5.41·Q`. **The physical decay is exactly
`λ₂^k` and 1% settling is `Q·ln 100 = 4.6·Q`.**

⚠ **The non-normality objection to `λ₂` is not small — it is ABSENT.** No transient growth
exists on these monodromies. `λ₂` is the right quantity and the rate is exactly `λ₂^k`.

⚠ **Practice, since a prefactor can be reintroduced by any norm:** measure a **residual**, not a
period count, and if a period count is unavoidable, write down the norm beside it. A basis change
(scaling the current block by 1e3) moves the apparent prefactor 1.81 → 8.37 without touching the
spectrum.

⚠⚠ **AND THEN THE MEASUREMENT SAID OTHERWISE — the gate is NOT badly broken, and the protection
is ACCIDENTAL.** Measured with a *generic* perturbation direction (both components), which is
what the shipped gate actually uses:

| `Q` | `λ₂` | `\|d·x̂\|` | raw error n=1 → n=6 | Aitken(1,2,3) |
|---|---|---|---|---|
| 0.12 | 0.0003 | 0.865 | −0.46% → −0.49% | −0.49% |
| 15.9 | 0.939 | 0.736 | −0.41% → −0.63% | **−1.05%** |
| 63.7 | 0.984 | 0.734 | −0.73% → −0.88% | **−1.87%** |

⚠ **A random direction in 2-D is ~71% tangential**, so the phase signal dominates and the
transverse contamination enters as a small additive term. The gate is accurate to <1% at
`Q = 64` **because the fixture has two states**, not because the gate is sound.

⚠⚠ **THAT PROTECTION SCALES AS `1/√m` AND VANISHES ON A REAL CIRCUIT.** A random unit vector's
tangential fraction is `~1/√m`: 0.71 at `m = 2`, 0.32 at `m = 10`, 0.10 at `m = 100`. So on any
circuit with more than a handful of states a random kick is **mostly transverse**, the
contamination dominates, and the three-period premise does bite. **The gate's validity is a
property of the fixture's DIMENSION** — a fourth variation on the fixture theme, and the first
one that is about size rather than structure.

⚠⚠ **AND AITKEN IS NOT A FREE IMPROVEMENT.** Applied where the geometric mode dominates the
residual it removes ~90% of the error (measured elsewhere on a transverse kick); applied where
it does not, it **amplifies** the remaining systematic — 0.63% → 1.05% and 0.88% → 1.87% above.
So the repair is **conditional on the contamination being the dominant residual**, which must be
established before applying it rather than assumed.

⚠ **The transverse-kick configuration that motivated all this is itself degenerate three times
over** and should not be used: near-zero signal (a transverse kick barely shifts phase — that is
what high `Q` means), a prediction that is a **5× cancellation** (`v·d̂ = +0.01010 − 0.00820 =
+0.00190`), and the direction where `C = diag(1,−1)`'s sign structure is most active. **The safe
kick is constructed from `ẋ_s(0)` and is neither a coordinate axis nor its exact orthogonal
complement** — the convenient choice and the clever choice are both degenerate.

⚠ **The stated premise is still false, and that stands separately from the measured error.**
`test_the_ppv_predicts_a_phase_shift_the_oscillator_actually_has` integrates **3 periods** and
reads the surviving tangential displacement, on the stated premise that *"the transverse
components have died (the second multiplier is 8.6e-4 per period)"*. At `Q = 16`, `λ₂³ = 0.83`
— **nothing has died**. The gate's *method*, not just its coverage, is `λ₂`-bound.

⚠ **One gate that does NOT degrade, and why that is not reassurance.** `oscillator_covariance`'s
`d/T` against `diffusion_constant` *improves* with Q — 1.0054 → 1.0013 → 1.000016 → 0.99967 at
`Q` = 0.14 → 15.9. But both routes ride the **same monodromy**, so a shared `λ₂` error cancels:
this is the shared-instrument shape (§D 0b) wearing a reassuring number. Its agreement says
nothing about `λ₂` accuracy.

### 0b. Which gates the high-Q fixture actually re-tests — measured 2026-09-03

The exposure recorded above ("every gate that passes on van der Pol at `μ = 1` was tested at
`|λ₂| = 8.6e-4`") is real but **narrower than it reads**, and the difference is worth stating
before anyone spends time re-running the suite at high `Q`.

⚠⚠ **MOST OF THE OSCILLATOR-NOISE GATES ARE INSENSITIVE TO `Q` BY CONSTRUCTION.** Run on the
`m = 12`, `Q = 60` fixture and on van der Pol at `μ = 1`, they agree **to every printed digit**:

| gate | `Q = 0.12` | `Q = 60` |
|---|---|---|
| `∫S df`, harmonics 1/2/3 | 0.99993570 / 0.99998154 / 0.99998720 | **identical** |
| far skirt per doubling | 4.000000 | **identical** |
| harmonic 3 vs 1 | 9.5424 dB | **identical** |
| `phase_psd` vs `lorentzian` | 1.000e-06 | **identical** |

They are self-consistency checks on the **closed-form** `lorentzian`, which knows the circuit
only through the scalar `c`. Re-running them at high `Q` establishes nothing — the fixture never
enters them.

⚠ **The `Q`-sensitive gates are the ones that touch the CIRCUIT, and those are done:** Ritz
(err 3.5e-12 at `k = m`), `λ₂`/`Q` convergence (the amplification ratio 60 reappearing),
`oscillator_covariance` (`d/T` vs `c` = 0.999603), A4d's zero pattern (`Γ/c` = 4.9e-27), the PPV
physical gate (**breaks**, 24%, and cannot be repaired), and the frequency-shift gate (**1.000006
at Q = 75**).

⚠⚠ **WHAT REMAINS EXPOSED, precisely:** `diffusion_constant`'s **absolute scale** at high `Q`.
It is anchored by a nonlinear Monte Carlo at `μ = 1` only. At `Q = 60` the only cross-check is
`oscillator_covariance`'s `d/T`, which **shares the monodromy** — so it rules out a defect in
either route and not a defect common to both. The frequency-shift gate anchors `⟨v⟩` (the *DC*
functional) at `Q = 75`; it does not anchor `c` (the *quadratic* one).

✅✅ **CLOSED 2026-09-04 — see §0e.** `CY` is exact, `v₀` is anchored **pointwise** at `Q = 8`
and `Q = 30`, and the quadrature converges at **third order** (converged to `~1e-4` at the shipped
480 points). Three measurements, no shared instrument, no Monte Carlo. The paragraph below is kept
as written because the route to closing it is the useful part.

⚠⚠ **AND THE NARROWED CHECK FOUND A DEFECT — see §0d.** `CY` is exact, but `c` comes back
**exactly 0.0** for an oscillator whose only noise is its series tank loss, because the PPV's
entry for a purely algebraic node row is structurally zero. Three independent confirmations and a
validated correction are in §0d, and it is **FIXED** — `c` series/parallel now agrees to 0.999973
at `Q = 8`, with every existing PPV number bit-for-bit unchanged.

⚠ **NARROWED 2026-09-04 by §0c below.** The state-localised frequency-shift probe anchors `v₀`
**pointwise** at `Q = 8` and `Q = 30`, to between 0.015% and 0.13%, against an instrument that
shares no code with the adjoint replay. So the PPV is no longer the suspect: what is still
unanchored at high `Q` is `diffusion_constant`'s **assembly** — the quadrature and the `B CY Bᵀ`
contraction — which is a far smaller surface than "the absolute scale" and needs no Monte Carlo.

**So the honest statement is one sentence, not the blanket one:** every *shape* result holds at
any `Q` because it never sees `Q`; every *circuit* result has been re-run; and one *scale* — `c`
— rests at high `Q` on a check that shares an instrument. A Monte Carlo at `Q = 60` would close
it, at 4.6·Q ≈ 280 periods per realisation.

**Consequence for planning:** the failures below are not a collection of unrelated sharp edges
to be patched one at a time. Any item whose gate passes on van der Pol at `μ = 1` has been
tested at `|λ₂| = 8.5e-4` — six orders from where a real LC oscillator sits — and has therefore
not been tested at all in the regime that matters.

---

### 0c. The windowed PPV probe — ⚠ **the time-localised form is IMPOSSIBLE**, the state-localised form works — 2026-09-04

§0b leaves one thing exposed: `diffusion_constant`'s **absolute scale** at high `Q`. The
frequency-shift gate anchors `⟨v⟩`, the *DC* functional, at `Q = 75`; `c` is the *quadratic* one,
and a quadratic functional is only anchored if the PPV is anchored **pointwise**. So the natural
next instrument is a windowed version of the same gate: inject over a window rather than over the
whole period, and compare `ΔT/A` against `∫ w(t) v₀(t) dt`.

⚠⚠ **THE OBVIOUS CONSTRUCTION — A SOURCE LOCALISED IN TIME — CANNOT WORK, AND NOT FOR A
NUMERICAL REASON.** A window in `t` has a period of its own, so the circuit stops being
autonomous. Measured on the `Q = 1` fixture with a `C¹` periodic raised-cosine bump at
`A = 1e-6`:

    unforced:  converged=True   autonomous=True    T0 = 6.298385479
               |λ| = [1, 0.4482, 0]      min|1 − λ| = 1.108e-11
    bumped:    converged=False  autonomous=FALSE   period = 6.298385479  (Δ = 0.0 EXACTLY)

Two independent obstructions, and either alone is fatal:

  * **the period stops being an unknown.** PSS infers a driven problem and solves at the period
    it was handed, so `ΔT ≡ 0` *by construction* — the quantity the probe exists to measure is
    defined away. The `Δ = 0.0` above is exact, not small.
  * **the Newton has no Jacobian.** A driven solve's Jacobian is `M − I`, and the *unforced* `M`
    carries `λ₁ = 1` to `1.1e-11` — the phase mode. At `reltol = 1e-12` the solve cannot
    converge, and at `A = 1e-6` the forced `M` is still that singular.

Physically this is **injection locking**: a perturbation periodic at the free-running period pulls
the *phase* and leaves the *frequency* alone. The locked phase satisfies `∫ v₀(t+φ)w(t) dt = 0`,
which pins the **zeros** of the PPV correlation and says nothing about its scale — so even the
converged version of this experiment would not do the job it was built for.

**WHAT WORKS INSTEAD: LOCALISE IN THE STATE, NOT IN TIME.** A source `i = A·g(v)` is autonomous —
the period stays an unknown — and because the orbit passes through each `v` at known phases it is
localised in phase all the same. Built as a second `BSource` on the same node; note
`terminals = ('inp','inn','outp','outn')`, so `v_ctrl = +v`.

⚠ **The control comes free and is the reason to trust the rest:** `g ≡ 1` is the DC gate, and it
reproduces it to `1.1e-4` (`Q = 1`) and `2e-5` (`Q = 8`), fixing the sign convention at **−1**.

⚠⚠ **AND THE DISAGREEMENT THAT LOOKED LIKE A DEFECT IS MY QUADRATURE, WHICH ONLY A REFINEMENT
STUDY COULD SAY.** The bumps missed by 4–8% at `Q = 1` and by up to 32% at `Q = 8`, with
`lin(A/2) ≈ 1.000` ruling out nonlinearity. Refining 16× at `Q = 8`:

    npts     predicted        measured        ratio
     240   +6.582886e-03   -4.061242e-03   -0.61694
     480   +5.341411e-03   -4.075562e-03   -0.76301
     960   +4.717965e-03   -4.076759e-03   -0.86409
    1920   +4.406292e-03   -4.076873e-03   -0.92524
    3840   +4.250592e-03   -4.076954e-03   -0.95915

**The measured column moves 0.4% and the predicted column moves 55%.** The PSS is converged and
the prediction integral is first order — rate 0.87. Same shape at `Q = 30`, so it is not a
high-`Q` effect.

⚠ **TWO PLAUSIBLE CAUSES, BOTH FALSIFIED BY THEIR OWN DIAGNOSTIC**, which is worth recording
because each was convincing enough to have been written up as fact:

  * **a fixed fractional sample association.** `samples` has 479 rows for 480 times, so a sample
    indexes an *interval*; pairing it with the left and the right endpoint **brackets** the
    measurement, so the truth is interior. But the `θ` that would land exactly on `−1` DRIFTS —
    0.67166, 0.67635, 0.68596, 0.70517, 0.74358 — so there is no fixed offset to correct.
  * **an oscillatory `O(h)` component in `v₀`**, which would cancel in a whole-period integral
    and not against a narrow weight. The zig-zag measure decays 4.016×, 4.008×, 4.004×, 4.002×
    per doubling — **exactly `O(h²)`**, which is ordinary curvature. There is no saw.

**§D SHAPE 0g — A CONTROL THAT IS A WHOLE-PERIOD INTEGRAL IS BLIND TO ANY ERROR THAT INTEGRATES
TO ZERO.** The DC gate agrees to `1e-4` while the localised functional is 32% out, and both read
the same `v₀`. An `O(h)` *time shift* of the PPV is invisible to `∫ v₀ dt` — exactly — and is
first order in every localised functional. So passing the DC gate was never evidence that the
pointwise PPV was accurate, and the roadmap's own sentence "the frequency-shift gate anchors
`⟨v⟩` … it does not anchor `c`" was more literally true than intended.

**THE INSTRUMENT IS USABLE ONCE BOTH ERRORS ARE EXTRAPOLATED AWAY.** The prediction's error is a
clean `O(h)` and the measurement's a clean `O(A)`; Richardson in `h` alone leaves a *stable*
0.44%, which is exactly the size `lin(A/2) = 1.00218` predicts. Extrapolating both —
`pred₀ = 2·pred(h/2) − pred(h)`, `meas₀ = 2·meas(A/2) − meas(A)`:

    Q     v*        pred_0          meas_0          ratio      grids
     8    -1.0733   +4.094979e-03   -4.093816e-03   -0.99972   480/960
     8    -0.3577   +3.576679e-03   -3.576144e-03   -0.99985
     8    +1.0733   -2.616944e-03   +2.613608e-03   -0.99873
    30    -1.0733   +1.092002e-03   -1.091794e-03   -0.99981   960/1920
    30    -0.3578   +9.535700e-04   -9.537214e-04   -1.00016
    30    +1.0733   -6.973740e-04   +6.971334e-04   -0.99965

⚠ **Six independent localised functionals of `v₀`, at three phases and two `Q`, agreeing to
between 0.015% and 0.13%** — against a measurement that shares no code with the adjoint replay.
That is a *pointwise* anchor, which is what `c` needs and what `∫ v₀ dt` could never give. The
`Q = 30` column is the one that matters: it is the regime §0b says is untested, and the agreement
does not degrade there.

⚠⚠ **WHAT THIS DOES AND DOES NOT CLOSE.** It anchors `v₀` **pointwise** at `Q = 30` against an
independent instrument. It does **not** measure `c`, which is a *quadratic* functional of `v₀`
against `CY`. So the exposure in §0b narrows rather than closes: the PPV is no longer the
suspect, and what remains unanchored at high `Q` is `diffusion_constant`'s **assembly** — the
quadrature and the `B CY Bᵀ` contraction. That is a much smaller surface than "the absolute
scale", and it is checkable without a Monte Carlo.

**⚠ THE COST IS THE POINT, THOUGH: four PSS solves per functional** (two grids × two amplitudes),
and the grid has to be fine enough that the `O(h)` term dominates the `O(h²)` one. At `Q = 8`
that is 480/960; at `Q = 30` it is 960/1920. This is a *gate*, not something to run in the suite
at every `Q`.

---

### 0d. The assembly check found a DEFECT — ⚠ **noise on an ALGEBRAIC row is dropped** — 2026-09-04

§0c narrowed the last open scale question to `diffusion_constant`'s **assembly**: with `v₀`
anchored pointwise, all that was left in `c = (1/T) ∫ vᵀ (CY/2) v dt` was `CY` and the
quadrature. Checking them found something else.

**`CY` is exact.** For a fixture whose only noisy element is one resistor, `_cy_reduced` returns
`4kT/r` at `T = 300 K` to every printed digit, on exactly one entry — the reduced row of the node
the resistor is on — and zero elsewhere, at `Q = 8` and `Q = 30`. That piece is clean and needs no
further work.

⚠⚠ **AND THEN `c` CAME BACK AS EXACTLY 0.0 FOR A MANIFESTLY NOISY OSCILLATOR.** The PPV's entry
for a purely **algebraic** node row is structurally zero, and that is precisely where a series
loss resistor's noise current lands. The two facts meet in the contraction and the answer is
silently zero — no warning, no exception.

**MEASURED, THREE INDEPENDENT WAYS**, on a van der Pol tank whose only noisy element is its loss,
drawn two equivalent ways: `series` puts it in the inductor branch (node `x` has no capacitance,
so its KCL row is algebraic); `parallel` puts `Rp = L/(C·Rs)` across the capacitor (a
differential row). The two are matched to `3e-6` in amplitude and `1e-5` in `Q`:

    Q     c series (shipped)   c parallel (shipped)   d/T series (Lyapunov)   corrected series
     8    0.000000000e+00      5.147901730e-24        5.147855204e-24         5.147858663e-24
    30    0.000000000e+00      1.372766117e-24        1.372650093e-24         1.372755969e-24

  * **the equivalent circuit** — the same physics on a differential row gives a nonzero `c`;
  * **`oscillator_covariance`** — a SHIPPED function reaching `CY` through the Lyapunov
    recursion rather than the PPV gets it right, `d/T` matching the corrected value to `7e-7`;
  * **a direct phase-sensitivity measurement** — injecting a DC current at the algebraic node
    shifts the period by `+1.937090e-06` per amp while the PPV predicts `0`.

⚠ **THE MECHANISM, AND WHY THE PPV IS NOT SIMPLY WRONG.** The differential rows agree with the
measurement to `1.5e-5`. An algebraic row's perturbation reaches the dynamics through the
**constraint**: eliminating `v_x = r(i_L + b)` puts `−r·b` into the inductor's row, so the true
sensitivity to that row is proportional to the BRANCH row's. Measured: `|r·∫v_branch dt| =
1.937064e-06` against `1.937090e-06`. The vector carries zero where that belongs, so the generic
fill-in is `v_A = (G[A,Z]ᵀ)^{-1} G[D,Z]ᵀ v_D` — ⚠ **whose SIGN this fixture cannot settle**, for
the reason shape 0h below records; the divider does.

⚠ **SCOPE, MEASURED RATHER THAN ASSUMED — `pnoise` IS NOT AFFECTED.** On a driven linear circuit
whose only noise is a series R at an algebraic node, `pnoise` agrees with `analysis_ss.Noise` to
**1.000000**, and so does the parallel form. The drop is specific to the **PPV path** —
`diffusion_constant` and whatever else contracts `CY` against `ppv()` — not to the adjoint
transfer machinery in general.

⚠ **WHY NO GATE CAUGHT IT.** `_vdp_at_Q`, the high-`Q` fixture, takes its noise from
`IS('v', gnd, noisePSD=…)` — a differential row. `_lc_osc(rs=…)` does have the series resistor,
but it also adds an explicit `IS` at node `v` at `psd = 1e-6` against the resistor's `4kT/0.2 ≈
8.3e-20`, so the dropped term is **14 orders down** and invisible. The `d/T` vs `c` gate that
would have caught it runs on the fixture that cannot show it. **No shipped test is numerically
wrong; the gap was in fixture placement, not in any assertion.**

**PINNED BY THREE TESTS** — `test_the_ppv_carries_the_slaved_sensitivity_on_an_algebraic_row`
(the divider's topology-fixed ratio, and every row on one convention),
`test_the_fill_is_skipped_entirely_without_an_algebraic_row` (the bit-for-bit guarantee, asserted
on the pattern), and `test_diffusion_constant_sees_noise_on_an_algebraic_row` (series against
parallel, cross-checked against the Lyapunov route).

⚠⚠ **FRAMING CORRECTED IN §0h — THIS AND §0g ARE ONE DEFECT.** The algebraic entries are zero
because `ppv()` returns `Cᵀ v_1` and `Cᵀ` annihilates the algebraic-state columns; the fill below
is an empirical reconstruction of half of `v_1`. It stands and it is gated, but read §0h for the
mechanism before building on it.

✅ **FIXED 2026-09-04 in `ppv()` itself**, by `_algebraic_adjoint_pattern` and
`_algebraic_adjoint_fill`. Results: `c` series/parallel **0.999973** at `Q = 8` and **0.999984**
at `Q = 30`; the PPV route and the Lyapunov route now agree; and the DC-injection probe reads
**+0.9999849** (differential row), **+0.9999982** and **+0.9998744** (the two algebraic rows).

⚠ **THE NORMALISATION QUESTION DISSOLVED RATHER THAN BEING PAID.** The worry was that filling
`v_A` moves `v·ẋ`. It does — by a measured 1.58e-05 — but it **should not be allowed to**, and
the reason is structural: `v·ẋ = 1` is about a **state** perturbation, and a state perturbation
of a DAE lies **on** the constraint manifold, its algebraic components determined by its
differential ones rather than free. The algebraic entries of `v` answer a different question — the
sensitivity to a perturbation of an **equation row**, which is exactly what a noise current
injected into an algebraic KCL row is. So the normalisation line is untouched and **every PPV
number on every circuit is bit-for-bit what it was**, which is the guarantee the fix was gated on.
`test_the_fill_is_skipped_entirely_without_an_algebraic_row` asserts it on the PATTERN, because an
empty pattern means `ppv()` runs exactly the lines it ran before.

⚠⚠ **AND FILLING IN THE WRONG PLACE WAS MEASURED WORSE, WHICH IS THE PART WORTH KEEPING.** The
first implementation filled `v` before the normalisation — which also made it the **replay's
seed**. The DC probe went from 0.9999849 to **0.9992364**: the algebraic components are SLAVED,
so propagating them through the step map corrupts the differential ones. There is nothing to
propagate; they are a pointwise function of the differential entries, and the fill belongs
**after** the replay, applied to each sample at its own operating point.

⚠⚠⚠ **§D SHAPE 0h — A FIXTURE WHERE THE CANDIDATE ANSWERS ARE NUMERICALLY DEGENERATE.** The
sign of the fill could NOT be settled on the single-series-resistor tank: there
`|∫v₀| = 1.9371e-06` and `|r∫v_branch| = 1.9374e-06` agree to `1.5e-4`, so a prediction matching
in magnitude matches for **either sign**, and the agreement I first recorded as a confirmation was
no evidence at all. **A RESISTIVE DIVIDER breaks the degeneracy:** two algebraic nodes fold into
the branch row with `(r1 + r2)` and `r2`, so the topology fixes the RATIO of their entries —
measured **10.000000** against a chosen 10.000000. The magnitude is then structural and only the
overall sign is measured, fixed by requiring algebraic and differential rows to share one
convention. **Before trusting a sign, check that the fixture can express the wrong one.**

⚠ **INDEX > 1 IS REFUSED, NOT GUESSED.** When the algebraic equations and algebraic states differ
in count, or the block `G[A, Z]` is singular, the entries are left at zero and a `RuntimeWarning`
says the noise entering those rows is under-counted. That is §B4's territory and a guess there
would be worse than a known zero.

---

### 0e. The quadrature — ✅ **§0b's EXPOSURE IS CLOSED** — 2026-09-04

The half of the assembly check that §0d's defect interrupted. `asm2` died on a `ZeroDivisionError`
at its first row **because `c` was zero** — the defect was hiding inside the measurement meant to
audit it — so `c`'s convergence order was never taken. With the fix in, it is.

⚠ **§0c MADE A FALSIFIABLE PREDICTION HERE AND IT HELD.** The PPV's pointwise error behaves like
an `O(h)` TIME SHIFT, and a shift is *exactly* invisible to a whole-period integral — so `c`, a
whole-period integral, had to converge at **second order or better** even though every localised
functional of the same `v₀` converges at first. A first-order result would have meant the shipped
value was wrong by percents and would have falsified that reading. Measured:

    fixture      Q     orders (per doubling)        c at npts = 960
    parallel      8    3.00, 3.57, 1.97             5.147901730e-24
    parallel     30    2.98, 2.93, 3.00             1.372766117e-24
    series        8    3.14, 5.11, (floor)          5.147858663e-24
    series       30    2.99, 2.97, 3.16             1.372755969e-24

**Third order**, at both `Q` and in both topologies. ⚠ **AND IT FLOORS**, which is worth stating
rather than hiding: past `npts ≈ 960` the differences reach `~1e-28` absolute (`~2e-5` relative)
and the order estimate goes to `−1.00` and `nan` at `Q = 8`. That is the solver's own noise, not
the rule — so refining past ~960 points buys nothing, and the *reported* order beyond that point
is meaningless.

⚠⚠ **WHAT THIS CLOSES.** §0b has carried "`diffusion_constant`'s **absolute scale** at high `Q`"
as the one open scale question since 2026-09-03, wanting a 280-period Monte Carlo. Three
measurements now stand in its place, none sharing an instrument with the others:

  * **`CY` is exact** — `4kT/r` at 300 K to every printed digit, right row, zero elsewhere (§0d);
  * **`v₀` is anchored POINTWISE** at `Q = 8` and `Q = 30` to 0.015–0.13%, against a
    frequency-shift instrument that shares no code with the adjoint replay (§0c);
  * **the quadrature is third order** and converged to `~1e-4` relative at the shipped 480
    points, per the table above.

⚠ **AND ONE AXIS WAS STILL UNSWEPT, which §0g found the next day: every fixture here uses
`C = 1 F`, so none of this says anything about the capacitance scale — and `diffusion_constant` is
wrong by `C²`.** What closes below is the `Q`-dependence and the quadrature, which is what §0b
asked about; it is not a statement about dimensional correctness.

**So the absolute scale is anchored end to end, at high `Q`, with no Monte Carlo** — which is what
§0b asked for and priced at 280 periods per realisation. The series and parallel routes agree to
`0.9999916` at `npts = 960`.

---

### 0f. Reading from the docs session, 2026-09-04 — ⚠ **RELAYED, NOT VERIFIED HERE**

Four findings arrived from the docs session against ~60 new papers. Recorded because two bear on
shipped code, with the standing caveat this file already uses for that channel: **the quotations
are relayed and have not been checked against the sources from here.**

⚠⚠ **1. KRYLOV MAKES THE DEFLATION WORSE — AND IT CORROBORATES A NUMBER WE MEASURED.** Mei &
Roychowdhury 2007 TCAD are quoted saying partial Floquet decomposition via **Krylov subspace
methods**, "while substantially faster, EXACERBATE THE IMPERFECT-CANCELLATION ISSUE". That is not
a new claim to us — it is the **published explanation for the tolerance floor B6 measured**:
`matrix_free=True` at `reltol = 1e-12` exhausts the inner GMRES restarts at `Q ≥ 16` with exactly
126 matvecs, while `1e-9` converges in 28–35, flat in `Q`. Two sessions had a fixture disagreement
over whether that was an artefact; on this reading the two fixtures sit on either side of a
**known boundary**, which is a better explanation than either fixture being wrong. ⚠ It does not
change the shipped default (`matrix_free` is off) and it does not need a code change; it needs to
be in B6 so the floor is not re-litigated as a bug.

⚠ **2. A NAMED LIMIT ON A7.** The same paper is quoted proving the singularity as its Lemma 2.1 —
`J(s)` loses rank by one at **every** `s = j k ω₀`, i.e. at DC and every harmonic, which is what
both sessions derived independently — and then saying deflation only **AMELIORATES** it. A7's
`b79b458` carries the harmonic pole analytically, which is deflation, so it inherits a documented
residual weakness. Their remedy is a different formulation (GeMPDE with augmenting phase
conditions); two other published routes are named, least squares with no phase condition and the
probe. **No change proposed** — recorded so the weakness is known rather than discovered.

⚠⚠ **3. C2 DOES NOT COVER LSOAC, AND CHECKING SAID SO MORE SHARPLY THAN THE QUESTION DID.** The
docs session asked whether `96a06ac`'s rejection might have been of a different method. Checked:
`96a06ac` rejects the **Poincaré / orthogonality phase ROW** against the frozen-coordinate pin,
recorded as **C2**. Mei & Roychowdhury 2006 DATE (LSOAC) is quoted as removing phase conditions
**entirely** — minimum-norm least squares on the underdetermined system, motivated by "the use of
phase condition equations can cause various numerical artifacts". ⚠ **So C2 is correctly scoped
and stays closed**; the finding is not that a rejection needs reopening but that **LSOAC was never
considered at all** — `grep` finds no mention of least squares or LSOAC anywhere in this file. It
is now B10, open, so the two are not confused again.

⚠ **4. AN EXTERNAL PPV ORACLE, ON A FIXTURE THAT ALREADY EXISTS — THE MOST ACTIONABLE ITEM.**
Ghanta, Li & Roychowdhury 2004 ASP-DAC are cited as giving **analytical** PPV expressions for
generic LC oscillators, with **symmetry** in the negative-resistance mechanism giving "particularly
simple forms". `_lc_osc(a=…)` sweeps exactly that symmetry — `a` breaks the nonlinearity's
half-wave symmetry — so the fixture for it is already built and already used. This is the standing
gap §0b names: the PPV physical gate cannot verify the PPV at high `Q`, and every check we have is
internal or shares the monodromy. **An analytical PPV would be the first fully external oracle.**
Blocked only on the expressions themselves, which are not in the relayed message.

Two smaller ones, no action implied: the idealised 3-stage ring's `{1, φ⁻⁶, φ⁻¹²}` is published
(Srivastava & Roychowdhury 2007 TCAS, golden mean "central to our exact analytical phase model") —
**cite rather than re-derive**; and the `τ_p/T` semantic failure recorded in `ppv`'s `Q` note has a
published counterpart *with a remedy* (Lai & Roychowdhury 2006 DAC, DCO gated capacitors, hierarchical HB).

⚠ **The pattern the docs session named is worth keeping:** four of the first ten papers are primary
sources for results one session or the other derived from scratch. The improvement is that the
reading is now arriving **before** the implementation commits rather than after.


---

### 0g. The external oracle — ✅ **BUILT**, and it found a SECOND PPV defect — 2026-09-04

The docs session supplied Ghanta, Li & Roychowdhury 2004 ASP-DAC **Lemma 5.2**, for an LC
oscillator with an ODD-symmetric `i-v` and a sinusoidal steady state:

    c = (N²/2) (L/C) / A²

This is the first check on `c` that shares **nothing** with our monodromy — §0b's standing
complaint that "the PPV physical gate cannot verify the PPV at high `Q`".

⚠⚠ **THE SWEEP IS THE TEST, NOT THE CONSTANT — AND THAT IS WHERE IT BROKE.** The docs session
measured ratio 0.50003 constant over 40× in `Q` and read the residual factor 2 as a
one-sided/two-sided convention. That reading is right, but `Q` barely moves `A`, `L` or `C`, so
their sweep pins a **scale factor** — and a scale factor is exactly what a PSD convention looks
like. ⚠ **A factor of two against an external reference is also the precise shape that already bit
this campaign once** (full `CY` vs `CY/2`, which only `kT/C` could see), so it was settled rather
than labelled. Sweeping `L` and `C` **independently** tests the functional form:

    C     L     c (PPV route)   d/T (Lyapunov)   c/(d/T)      d/T ÷ Lemma 5.2
    0.1   1     6.250857e-09    6.248771e-07     0.010003     0.499840
    1     1     6.250857e-08    6.248771e-08     1.000334     0.499840
    10    1     6.250860e-07    6.248764e-09     100.033536   0.499839
    1     0.1   6.250854e-09    6.248769e-09     1.000334     0.499840
    1     10    6.250862e-07    6.248767e-07     1.000335     0.499839

⚠⚠⚠ **`oscillator_covariance` MATCHES THE ORACLE AT 0.49984 ACROSS 100× IN `C` AND 100× IN `L`,
AND `diffusion_constant` IS WRONG BY EXACTLY `C²`.** Varying `L` is fine in both; varying `C` is
not. That asymmetry is the whole diagnosis, and it is the same shape as §0d — two shipped routes
to one number, the Lyapunov one right and the PPV one wrong.

⚠ **CONFIRMED INDEPENDENTLY BY §0c's INSTRUMENT.** The state-localised bump probe measures a
LINEAR functional of `v₀` against a real period shift. Over the same `C` sweep the measured shift
is essentially constant (−1.53e-03, −1.52e-03, −1.49e-03) while the prediction `∫v₀·g dt` scales
as `C`: ratios **9.337, 0.9287, 0.0912** — `1/C` over two decades. So the error is in `v₀`
itself, not in the quadratic assembly.

**THE MECHANISM IS IN `ppv()`'s OWN DOCSTRING.** It records that "the vector this bordered solve
returns behaves as `Cᵀ v₁`" and that it "is contracted with a state perturbation directly". But
`CY` is an **equation-row** covariance — a current injected into a KCL row — and Demir's
`c = (1/T) ∫ v₁ᵀ B Bᵀ v₁ dt` uses `v₁`, not `Cᵀ v₁`. An impulse `b` in the equation rows produces
a state jump `C⁻¹ b`, so the sensitivity to `b` is `C⁻ᵀ v`, and contracting `v` instead
over-counts by `C` — squared, in a quadratic functional. **Exactly the measured `C²`.**

⚠⚠ **WHY NOTHING CAUGHT IT: EVERY FIXTURE IN THIS CAMPAIGN USES `C = 1 F`.** `_vdp_at_Q`,
`_lc_osc`, `_loss_osc`, the high-`Q` recipe — all of them. At `C = 1` the factor is 1 and the
`kT/C` anchor, the Monte Carlo, §0c's pointwise gate and §0e's quadrature study all pass while
saying nothing about it. **A dimensionless fixture cannot test a dimensional error.** New §D
shape 0i.

**NOT FIXED HERE, and the reason is specific rather than caution.** The correction is `C⁻ᵀ v`
where `C` is **singular** — that is the whole index-1 structure, and it is the same algebraic-row
territory where §0d's sign came out wrong when derived from scratch. The docs session has named
the reference that settles it: **Demir 2000, "Floquet Theory and Non-Linear Perturbation Analysis
for Oscillators with Differential-Algebraic Equations", IJCTA 28:163–185** — §2.2, §3.2, §3.4 and
**Remark 3.1** for the orthogonality/biorthogonality conditions, with Traversa & Bonani 2011's
`Vᵀ C U = [[I,0],[0,0]]` as the block form. Requested; the fix waits on it.

**PINNED BY TWO TESTS:** `test_the_lyapunov_route_matches_an_analytic_external_oracle` (passing —
the oracle itself, asserted as a CONSTANT ratio across `L` and `C`, not as a value) and
`test_diffusion_constant_should_not_depend_on_the_capacitance_scale` (**strict xfail**).

✅⚠⚠ **THE FACTOR OF TWO IS NOW CITED, AND THE ORACLE IS EXACT — 2026-09-04.** Winkler
(Oberwolfach Report 18/2006 p.1160, relayed) states Nyquist as `I_th = √(2kT/R) ξ(t)`, i.e.
`2kT/R` **TWO-SIDED** — which is exactly the `cy/2` that `diffusion_constant` contracts, since
`_cy_reduced` returns the ONE-SIDED `4kT/R` (measured against the analytic value to every printed
digit in §0d). So Ghanta's `N²` is two-sided and the gate had been feeding it a **one-sided**
`noisePSD`, double-counting by exactly the factor observed. Halving it:

    0.49984 × 2 = 0.99968        with no free parameter

**The gate now asserts ONE**, with the conversion named as `n_sq_two_sided` and the citation
beside it rather than absorbed.

⚠⚠ **AND THE OLD FORM OF THIS GATE IS ITS OWN LESSON.** "Assert the ratio is CONSTANT and equal to
0.5" **passed**, across four decades of `L/C` and two decades of `Q`, while carrying a factor of
two — because a constant wrong factor is exactly what a convention mismatch looks like, and
asserting constancy tests everything about the functional form and nothing about the scale. This
campaign had already lost time to one factor of two that only `kT/C` could see. **A gate that
tolerates an unexplained constant preserves it forever while reporting success.**

⚠ **ALSO FILED, next to the Monte-Carlo measurement it explains:** Sickenberger & Winkler (PAMM
2007) simulate noisy oscillators directly as SDAEs with **stochastic analogues of BDF2 and the
trapezoidal rule** — the same two integrators shipped here — and give the error bound
`O(h² + εh + ε²h^{1/2})` for small noise `ε`. That is the theory for a shape measured empirically
and recorded without explanation: a **two-sided optimum in `npts`**, where the Monte-Carlo error
fell and then GREW again with refinement. Deterministic-dominated at usable step sizes, with
`ε²h^{1/2}` as the floor. ⚠ The pairing only works because the unexplained measurement was written
down in the form it came out rather than smoothed.

⚠⚠ **CORRECTED 2026-09-04, SAME DAY — TWO WAYS, AND THE ATTRIBUTION ABOVE IS AT BEST
INCOMPLETE.**

**(i) THE EXPONENT THAT APPLIES TO OUR ESTIMATOR IS THE MEAN ONE, NOT THE MEAN-SQUARE ONE.**
Sickenberger, Weinmüller & Winkler (ASC Report 17/2007, "Local Error Estimates … Part II — SDEs
and SDAEs with Small Noise") give a **half-order gap**:

    ‖E(L_i | F_{t_{i-2}})‖_L2 = O(h^{γ+1})     consistency in the MEAN
    ‖L_i‖_L2                  = O(h^{γ+1/2})   consistency in the MEAN-SQUARE

We fit the **slope of `Var(θ)`**, which is an EXPECTATION, so it converges at the **mean** rate —
the deterministic order — not the mean-square rate that gives strong order 1/2. **The `ε²h^{1/2}`
floor applies to PATHS, not to the fitted slope.** ⚠ Anyone sizing a burn-in or a step against the
strong exponent would be using the wrong one.

**(ii) AND THE TWO-SIDED OPTIMUM PROBABLY HAS A BETTER EXPLANATION THAN THE ERROR BOUND.**
Römisch, Sickenberger & Winkler ("Simultaneous Step-size and Path Control") is not about the local
error at all — it is about tuning the **number of Monte-Carlo paths**: *"our aim in tuning the
number of paths is to balance the LOCAL ERROR and the SAMPLING ERROR"*, with a per-time-point
tolerance and a path count that varies over time. **Our Monte Carlo runs a FIXED path count and a
fixed burn-in.** Refining `npts` while the sampling error dominates buys nothing; adding paths
while the local error dominates buys nothing either. **The two-sided optimum in `npts` is one half
of that trade with the other half held fixed** — which explains the shape better than the error
bound alone, and means the bound above should not be quoted as "the" explanation.

⚠ **Recorded, not acted on.** Both are relayed and unverified here, and neither changes a shipped
number; they change what a future Monte-Carlo campaign should be designed against.

⚠ Also from the same design: *"we concentrate on two-step schemes, since the higher numerical
effort for higher deterministic order pays off ONLY IF THE NOISE IS VERY SMALL"* — which bears on
any future IRK adoption (B8's companion), and their stated motivation is ours: small noise in
circuit simulation, "where especially the BDF and the trapezoidal rule have proven valuable in the
deterministic case".

⚠ **AND THE ORACLE'S OWN PRECONDITION IS ASSERTED**, because it is easy to lose: the lemma needs a
SINUSOIDAL orbit and an ODD-symmetric nonlinearity. `rms/peak` is checked at 0.70785 against
0.70711. ⚠ `_lc_osc` is the WRONG fixture for it — its `i_func` coefficient is 1.0, so the orbit
is strongly non-sinusoidal, and its `a` breaks the odd symmetry the lemma requires; the docs
session measured the ratio drifting 0.646 / 0.768 / 1.029 there. `a` sweeps the right parameter on
the wrong operating point.


---

### 0h. ⚠⚠ **§0d AND §0g ARE ONE DEFECT** — the primary source, and a correction to my own framing — 2026-09-04

The docs session supplied Demir 2000 (IJCTA 28:163–185) verbatim, and it collapses two findings
into one. **Relayed, not verified here.**

**THE CONTRACTION CONVENTION, BOTH FORMS ON ONE PAGE.** For a STATE initial condition, eq (41)
carries `v_iᵀ(0) C(0) x(0)` — **with `C`**. For an EQUATION-ROW input `b`, eq (42) carries
`v_iᵀ(s) b(s)` and the phase equation (44) reads `dα/dt = v_1ᵀ(t+α) B(x_s) b(t)` — **bare `v_1`,
no `C`**. Demir's DAE form is `q(x) + g(x) + B(x)b(t) = 0`, so `B(x)b(t)` **is** an equation-row
input: our `CY` case exactly. So `c = (1/T)∫ v_1ᵀ B Bᵀ v_1 dt` takes `v_1`, and §0g's `C²`
diagnosis is confirmed by the primary source rather than by measurement alone.

⚠⚠⚠ **AND THAT MAKES §0d THE SAME DEFECT, WHICH I DID NOT SEE.** `ppv()` returns `Cᵀ v_1`.
`(Cᵀ v_1)_i` is the `i`-th COLUMN of `C` dotted with `v_1` — and the algebraic-state columns of
`C` are **zero by definition** (measured on the series fixture: "zero COLUMNS of C: [1]"). So:

  * on **differential** rows `Cᵀ` multiplies by the capacitance → a factor `C`, squared in `c`
    → §0g's `C²`;
  * on **algebraic** rows `Cᵀ` **annihilates** → exactly `0.0` → §0d.

**One cause, two symptoms.** §0d's record calls the algebraic entries "slaved to the differential
ones and left at zero", which is true but is not the mechanism: they are zero *because `Cᵀ`
annihilates that column*. The fill built there is an empirical reconstruction of half of `v_1`,
and it stands — it was gated against three independent references — but its framing was wrong and
is corrected here rather than in place, so the sequence of understanding stays readable.

⚠⚠ **DO NOT INVERT `C` — THE PRACTICAL POINT.** `C` is singular; that is the whole index-1
structure. Demir's adjoint, eq (24), is

    Cᵀ(t) (d/dt) y  −  Gᵀ(t) y  =  0

with `Cᵀ` a **multiplier on `dy/dt`, never inverted**, and the paper flags the trap in its own
sentence: "the time derivative operates on `y` only, NOT on the product `Cᵀ(t)y`, in contrast with
Equation (19)". ⚠ **That asymmetry is the sign I got wrong deriving §0d from scratch.** The fix is
to OBTAIN `v_1` — integrate (24) backwards, or re-pose the border so the returned object is `v_1`
— not to undo a `Cᵀ` through a singular matrix.

✅✅ **THE FILL IS IDENTIFIED, NOT MERELY VALIDATED — AND THE SIGN IS NOW DERIVED.** The open
question above ("is the fill `v_1`'s algebraic entries, or `Γ` in disguise?") is **answered**, by
a test that goes through none of the three outcome references. Write eq (24) componentwise: row
`i` is `(col i of C)ᵀ ẏ = (col i of G)ᵀ y`, and for an ALGEBRAIC state the column of `C` is zero,
so the left side vanishes and the row degenerates to a **pointwise constraint with no time
derivative in it at all**:

    (col i of G)ᵀ v_1(t) = 0        for every algebraic i, at every t

MEASURED on the series fixture, at nine times over the period, scaled by `max|v|`:

    v_1 = C⁻ᵀ v on the differential rows   ->  0.0000e+00   EXACTLY
    the returned vector read as v_1        ->  1.9870e+00
    the same with the fill's sign flipped  ->  1.9870e+00

**Machine zero, and both alternatives are O(1) out** — so it discriminates rather than merely
tolerates. Three things fall out at once:

  * **the fill IS `v_1`'s algebraic entries.** `Γ` is not in the phase equation at all (the docs
    session checked: every occurrence of `Γ` is in §3.3, and §4's perturbation analysis contains
    none; `Γ` lives in the ORBITAL deviation `z`, which decays when the perturbation is removed).
    ⚠ **The caveat stands for `oscillator_covariance`**, which computes the FULL state covariance
    including the orbital part — `Γ` can contribute there, and a series-loss tank injects exactly
    where `Γ` is nonzero. Clean for `c`; not assumed clean for the covariance.
    ⚠⚠ **NO SOURCE IN THE LIBRARY CONSTRUCTS `Γ`** — the docs session checked Traversa & Bonani
    2011 (no hits), Demir 2006 and Demir & Sangiovanni-Vincentelli 1998 (hits are "instantaneous
    frequency" and "instantaneous spectral density", different things). Eq (40) fixes its null
    space and rank and leaves its action on the complement free. **So "we cannot bound that term"
    is the correct position against this library, and it is written down rather than left
    implicit.** ⚠ **BUT IT IS MEASURABLE FROM ITS DEFINING ROLE, WHICH IS THE NEXT CHEAP ITEM:**
    `Γ(t)b(t)` IS the instantaneous non-propagating part of the response, so applying an impulsive
    `b = e_j` and reading the state jump at `t⁺` gives `Γ(t)e_j`; `n` applications give the matrix
    column by column and `T`-periodicity gives every other `t`. It comes with two acceptance gates
    from eq (40) — `rank Γ(t) = n − m`, and `Γ(t) C(t) u_i(t) = 0` for `i = 1..m` — so a measured
    `Γ` that fails either is a bad measurement rather than a bad theory. **That is the same
    identity-not-agreement shape that worked above**, and it would turn the covariance caveat into
    a bound or a defect without needing a source.
  * **THE EMPIRICALLY-FLIPPED SIGN WAS DERIVED AFTER ALL.** `C = diag(1, 0, −L)` on this fixture:
    the INDUCTOR BRANCH ROW CARRIES `−L`, so reading `v = Cᵀv_1` as `v_1` negates that row — and
    the term in the fill is dominated by the branch. The derivation and the measurement were
    describing **different vectors**, and both were right. That is a better outcome than either
    "the derivation was wrong" or "the measurement was noisy", and it is why §D shape 0h wants the
    fixture to be able to express the wrong answer.
  * **`v = Cᵀ v_1` IS CONFIRMED TO MACHINE PRECISION**, which is §0h's whole claim, established
    here independently of the `C`-sweep that suggested it.

⚠ **AND IT HANDS US THE FIX WITHOUT INVERTING A SINGULAR MATRIX.** `C⁻ᵀ` is applied only to the
DIFFERENTIAL block, which is invertible by construction; the algebraic entries come from the fill.
So the recipe is complete and gated: `v_1 = C[D,D]⁻ᵀ v_D` on differential rows, the fill on
algebraic ones, with the constraint above as the acceptance test. ⚠ **Existing gates would not
move:** every fixture puts its noise on a row with `C = 1 F`, so `c` is unchanged there — which is
also, exactly, why the defect survived.

**AND ONE SENTENCE SETTLES A THING WE MEASURED.** §3.5, verbatim: "`v_i(0)` are NOT the
eigenvectors of the transposed monodromy matrix `Φ(T,0)ᵀ`." `v_1` is an eigenvector of the
ADJOINT system's monodromy `Ψ(T,0)`, and §3.4 gives
`Ψ(t,s) = Vᵀ(t) D(s−t) Uᵀ(s) Cᵀ(s) ≠ Φᵀ(s,t)` — "NOT simply given by `Φᵀ(s,t)` ... as it was the
case for ODEs". The docs session's measurement — the left eigenvector of `M` aligns with
`C(0)ᵀv_1` to 1.000000000000 and with `v_1` to 0.9657 — **is that sentence**. Also: the
`Φ = U D V C` factorisation is eq (37), **Demir 2000**, not Traversa & Bonani 2011.

**Consequence for the fix:** it is now one change, not two, and it has a derivation rather than a
fitted sign. It is still not applied — re-posing the border is a change to the most-gated function
in this campaign, and §0d's fill will need re-deriving (or retiring) inside it rather than beside
it.


---

### 0j. The `C²` fix APPLIED, and Γ MEASURED — 2026-09-04

Andreas authorised both. **`c` now agrees with the Lyapunov route at every `C`:** the ratio is
1.000334 / 1.000334 / 1.000335 across `C` = 0.1, 1, 10 and 1.000334 / 1.000335 across `L` = 0.1,
10, where it had been 0.010003 / 1.000334 / 100.033536.

⚠⚠ **THE FIX IS TWO NAMED OBJECTS, NOT A CONVERSION — AND THAT IS A DESIGN DECISION, NOT THE
SOURCE'S.** Demir carries ONE vector and TWO contraction rules: eq (41), a state initial
condition, with `C`; eq (42) and the phase equation (44), an equation-row input, bare. This
implementation stores TWO vectors instead: `info['samples']` stays `Cᵀ v_1`, and
`info['samples_eq']` is the new `v_1`, produced by `_equation_row_ppv`. `diffusion_constant` and
`colour_projection` now read the latter.

  * **Why not convert:** `ppv()`'s return value is what a STATE perturbation contracts with, which
    is what its docstring promises AND what an in-repo measurement already established —
    predicting a state jump as `vᵀCδ` gives residuals of 0.36/0.40/0.42 that GROW with refinement
    while `v·δ` converges at `O(h)`. Converting would have silently changed a shipped, measured
    contract. Every existing gate keeps its anchor.
  * ⚠ **The cost of the choice, stated because it is real:** a caller can now pick the wrong
    array. Demir's one-vector-two-rules formulation has no such failure mode. Mitigated by naming
    and by a pointed comment at each contraction site; not eliminated.
  * ⚠ **`C` IS NEVER INVERTED.** The solve is on `C[D, NZ]ᵀ` — differential equations against
    non-algebraic states — which is square and invertible by construction. The algebraic entries
    come from the eq (24) constraint, which now carries the DERIVED sign because it acts on `v_1`.

⚠ **AND THE "SKIP WHEN THERE ARE NO ALGEBRAIC ROWS" SHORTCUT HAD TO GO**, which is the whole `C²`
lesson in one line: that shortcut was right for the FILL and wrong for the CONVERSION. A plain
ODE circuit with no algebraic row at all still needs `C⁻ᵀ` whenever its capacitance is not 1 F.

---

**Γ MEASURED — both eq (40) gates pass.** No source in the collection constructs it: the docs
session ran a full pass over **339 PDFs** for Drazin, spectral/index-1 projector, projector chain,
consistent initialisation, perturbation index, algebraic jump, jump condition, non-propagating,
impulse response of a DAE, and algebraic variables + covariance — **one hit, a false positive**
(Demir & Sangiovanni-Vincentelli 1998 p.75, "instantaneous jumps ... at the times of the
instantaneous carrier crossings", which is Poisson shot noise). ⚠ **Worth recording as a scoping
fact: DAE projector theory is not circuit literature, so it is not in this collection, and
obtaining Γ properly would be an ACQUISITION rather than another search.**

So it was measured off the solver's own matrices. Eq (42) has `Γ(t)b(t)` proportional to `b` at
the SAME INSTANT, so Γ is the `h → 0` limit of the step operator:

    Γ(t) = lim_{h→0} (C(t)/h + G(t))⁻¹

Converges first order (differences 9.0e-04, 9.0e-05, 9.0e-06, 9.0e-07, 9.0e-08 per decade) to a
single nonzero entry, `Γ[x,x] = 0.003978874`, which is `G[A,Z]⁻¹` — **exactly the series
resistance `r`**. Gates: rank → 1 = `n − m`; `|Γ C u|` → 1.788e-08, falling linearly in `h`.

⚠⚠ **AND A DIMENSIONAL SLIP OF MINE, RECORDED BECAUSE IT IS SHAPE 0i AGAIN.** I first compared
`Γ (CY/2) Γᵀ = 3.294507e-23` (which is exactly `2kTr`, so it looked right) against `K_orb[x,x]`
and reported a ratio of 76. **The comparison is invalid:** Γ is in ohms and `CY` in A²/Hz, so
`Γ CY Γᵀ` is a **PSD in V²/Hz**, while `K_orb` is a **variance in V²**. Two quantities that are
numerically comparable and dimensionally not — the same failure shape as `C = 1 F` hiding a
missing capacitance.

**⚠ THE COVARIANCE CAVEAT IS NOT CLOSED, AND IT NOW HAS A SHARPER OPEN QUESTION.** My hypothesis
was that an algebraic node, having no state, has no bandwidth limit, so its variance is not finite
and there is nothing for `K_orb` to carry. **Testing it FALSIFIED it:** adding a parasitic
capacitor at node `x` makes the node differential (`algebraic rows []`) but leaves `K_orb[x,x]`
*unchanged* — 2.508914e-21, 2.510968e-21, 2.511174e-21, 2.511194e-21 for `C_par` = 1e-3 … 1e-6,
against a `kT/C_par` of 4.14e-18 … 4.14e-15. **Flat over three decades, and a factor ~1e6 BELOW
`kT/C_par`.** So the parasitic node's thermal equilibrium does not appear in `K_orb` at all.

✅✅ **RESOLVED, AND THERE IS NO THIRD DEFECT.** The docs session supplied the definition —
Demir eq (68): mode 1 is the phase and modes `2..m` are where `z` lives, with **no**
"orbit-coupled only" clause — so a parasitic mode *should* be in `K_orb`, and their proposed
explanation was that the node lacked a thermal bath. ⚠ **That did not fit: node `x` already has
`Rs` across it, noisy at `4kT/r`.** The actual explanation is arithmetic:
`rs·C_par = 4e-9 s` against a timestep of `0.013 s` — **the parasitic mode is SIX ORDERS faster
than the grid**, and a mode the discretisation cannot represent cannot reach its equilibrium.

Tested by hanging a weakly-coupled RC branch off the tank with `τ = R_par·C_par` **chosen** rather
than inherited:

    tau/h = 152.79  ->  0.995110      <- kT/C to 0.5%
    tau/h =  15.28  ->  0.953586
    tau/h =   1.53  ->  0.688971
    tau/h =   0.15  ->  0.214472
    tau/h =   0.02  ->  0.029182

⚠⚠ **AND THE DISCRIMINATOR IS THAT IT DEPENDS ON `τ/h` AND NOT ON `C`:** at fixed `τ/h` the ratio
is **0.953586 identically across three decades of `C_par`**. A missing-term defect would scale
with something; a resolution limit is a pure function of `τ/h`, and that is what it is. So
`oscillator_covariance` reaches the `kT/C` external anchor — the same one that settled the `CY/2`
convention — whenever the mode is resolved, and degrades monotonically toward `τ/h` when it is
not. **Pinned by `test_the_orbital_covariance_reaches_kTC_when_the_mode_is_RESOLVED`**, which
asserts all three: the anchor, the `C`-independence, and that an unresolved mode does NOT reach it.

⚠ **TWO WRONG EXPLANATIONS DIED BEFORE THE RIGHT ONE, and that is the record worth keeping.**
Mine — an algebraic node has no state, hence no bandwidth limit, hence no finite variance — was
falsified by adding `C_par`. The docs session's — the node has no thermal bath — did not fit,
because node `x` already carries `Rs` at `4kT/r`. **Both would have closed the question**, and both
were argued as fitting the numbers. The arithmetic that decided it (`4e-9 s` against `0.013 s`) was
checkable before either was written.

⚠ **What this does NOT license.** It is a statement about a mode the grid resolves. A real
circuit's parasitics sit far below any PSS timestep — the original `rs·C_par` case is the typical
one, not the exotic one — so `K_orb` will routinely omit their thermal equilibrium, correctly and
silently. That is a property of sampling a periodic steady state, not a bug, but a caller reading
`K_orb` as "the" state covariance should know it contains only what the grid can see.

And Γ itself: measured, passing both its gates, and **absent from the phase equation**, so `c` is
unaffected — which was the question that mattered for §0g.

✅⚠⚠ **THE ACQUISITION LANDED AND IT GIVES A PRECONDITION, NOT Γ — AND OUR FIXTURE VIOLATES IT.**
Winkler 2004 (JCAM 163:435–463) **Definition 2**, relayed verbatim: an SDAE is **index 1** when
"the noise sources do not appear in the constraints", i.e. `im G ⊆ im A` — in our notation

    im B  ⊆  im C          the noise input must lie in the image of the capacitance matrix

and otherwise it is an SDAE **WITH DIRECT NOISE**, outside the class the theory covers.
⚠ **This is checkable on any netlist from matrices we already build**: `CY = B Bᵀ`, so
`im B = im CY`, and the test is the residual of projecting `CY`'s columns onto `im C`. Measured:

    fixture                  rank C   max|resid|/|CY|   verdict
    parallel loss             2/2       0.000e+00       index-1 SDAE
    _vdp_at_Q (IS at node v)  2/2       0.000e+00       index-1 SDAE
    series loss               2/3       1.000e+00       *** DIRECT NOISE ***
    series + parasitic C      3/3       0.000e+00       index-1 SDAE

⚠⚠⚠ **AND IT COMPLETES §0j.** The series-loss tank — the fixture that produced §0d's exact zero
and §0j's missing `kT/C` — is an SDAE **with direct noise**. White noise applied to a variable
determined by a CONSTRAINT rather than an integrator is filtered by nothing, so **that node has no
finite variance for `K_orb` to report** — not a missing term and not a resolution artefact. Adding
a parasitic capacitor moves the circuit back INTO the class (row 4), which is exactly the case
§0j's `τ/h` study measured and where `kT/C` duly appears. **Two answers, two regimes, and the
fixtures separate them cleanly:** the `τ/h` result explains the DIFFERENTIAL parasitic case; this
explains the ALGEBRAIC case it started from.

⚠ **WHAT IT DOES NOT OVERTURN.** `diffusion_constant` on the series fixture still agrees with the
equivalent parallel circuit to 0.999973 and with the Lyapunov route — so the PHASE diffusion is
well defined even where the state covariance is not, and §0d's fill stands. Winkler's condition
bears on the covariance, not on `c`. **Proposed gate, not yet built:** assert `im B ⊆ im C` and
WARN on the covariance paths, naming the class rather than a symptom — cheap, and the same shape
as the refusals already in `_cy_reduced`.

⚠ **AND THE PROJECTOR MACHINERY IS NOW AVAILABLE BUT STILL DOES NOT CONSTRUCT Γ.** Lamour, März &
Winkler 1998 (JMAA 217:372–394) eq (2.5) gives `X(t) = P_can(t) U(t) P(0)` — the structural home
of Demir's `Φ = U D V C` — but **homogeneous only**: no variation-of-constants, no input response.
Winkler's constructive route for the algebraic part is `x = Px + Qx = u + v̂(u,t)` with a
pseudo-inverse `A⁻ = D(I−R)`, `A⁻A = P`. So the measured Γ stands as the empirical object.

✅✅ **AND IT IS NOW CLOSED AS *CORRECT BEHAVIOUR, EXPLAINED* — THE DIVERGENCE IS THE ANSWER, AND
ITS EXPONENT IS A DIAGNOSTIC. 2026-09-04.** Alabert & Ferrante (arXiv math/0507159v2, 2006), relayed:
in Kronecker canonical form each nilpotent block of size `q` returns its algebraic variables as
`⟨v_j, φ⟩ = Σ_{k=j..q} ⟨c_k, φ^(k−j)⟩` — **the `j`-th algebraic variable carries the `(k−j)`-th
DERIVATIVE of its input.** Their Theorem 4.4: the law is absolutely continuous, but only against a
test function of order `r`, the nilpotency index. **The solution exists and is well behaved as a
GENERALISED PROCESS, with no pointwise value at all.**

⚠⚠ **SO THE QUANTITY A COVARIANCE ROUTINE ASKS FOR DOES NOT EXIST**, and no refinement produces
it. §0j's measurement was not failing to converge on a number; **there is no number.**

**AND THE PREDICTED SIGNATURE IS MEASURED HERE.** Discretising at step `h` tests white noise
against a bump of width `h`, so `Var ∝ 1/h` when the noise reaches an algebraic row directly and
`1/h³` one nilpotent level up. On the series-loss tank:

    npts   K_orb[x,x]        ratio
     120   6.231932e-22
     240   1.252673e-21      2.0101
     480   2.511196e-21      2.0047
     960   5.028027e-21      2.0022

**Exactly 2× per doubling — `1/h`, not `1/h³`** — so the noise reaches an algebraic row **directly**,
which is what a resistor injecting into node `x`'s KCL should do. ⚠ **An unbounded number becomes a
diagnostic that reports the structure of the netlist.**

⚠ **MY FIRST "CONTROL" WAS NOT A CONTROL, and its failure is §0j's result rather than a
contradiction.** Adding a parasitic capacitor at `x` gave ratios 2.0101 / 2.0047 / 2.0022 —
**identical to four digits** — because `τ = rs·C_par = 8e-8 s` against `h = 6.5e-3 s` is
`τ/h ≈ 1e-5`: the capacitor is INVISIBLE at that step, so the node is still effectively algebraic.
The genuine differential control is §0j's own — a resolvable `τ/h ≈ 153` converges to `kT/C` at
0.995. **The two together bracket it: algebraic diverges at `1/h`, resolved-differential
converges.**

⚠⚠ **THE GAP, FLAGGED AS THE RELAY FLAGGED IT.** Alabert & Ferrante assume **CONSTANT**
coefficients. A shooting linearisation is periodically TIME-VARYING, so this settles the LTI
special case and ours follows **by analogy**. Römisch & Winkler's `im B ⊆ im C` is the
time-varying condition and its classification is still the unanswered "future work". The analogy
is strong enough to explain the measurement and to predict the exponent — which it did, to three
digits — and it is **not a theorem about our system.**

✅⚠⚠ **THE COVARIANCE ITEM WAS CLOSED AS OPEN — AND IT IS OPEN IN THE LITERATURE, NOT IN OUR
READING.** Römisch & Winkler, "Stochastic DAEs in Circuit Simulation" (ISNM 146:303–318,
Birkhäuser 2003), relayed verbatim, gives the condition a **circuit-topological** form:

    im G(x,t) ⊆ im A   ⟺   THERE ARE ALWAYS CAPACITANCES IN PARALLEL TO A NOISE SOURCE

and calls it *"quite restrictive in the actual noise modelling"* — i.e. real noise models routinely
violate it. Our series tank-loss resistor has no capacitance across it, so it violates the
condition **structurally, not by an accident of the fixture**.

⚠⚠⚠ **AND THE NEXT SENTENCE SETTLES WHY THE ITEM STAYS OPEN.** Verbatim: *"one can also handle
many situations where this condition is violated. Often noisy constraints are only needed for the
determination of algebraic solution components that DO NOT INTERACT WITH THE DYNAMICAL ONES.
FUTURE WORK SHOULD BE DIRECTED TO A CLASSIFICATION OF SUCH SITUATIONS."* **As of 2003 the
classification of when a violated constraint is benign is EXPLICITLY OPEN WORK, by the authors who
defined the condition.** So "unbounded, measurement pending" was the correct entry and it is not a
gap in our reading — there is no result to have found. The practical read is that our case is very
likely the benign one (a series-loss resistor determines a branch current that does not feed back
into the tank dynamics), and **"very likely" is the strongest statement the literature supports** —
which is exactly why measuring Γ was the right move and why the `τ/h` result stands on its own.

⚠⚠ **BUT THE TOPOLOGICAL PHRASING IS LOOSE, AND TAKEN LITERALLY IT OVER-FLAGS — MEASURED BEFORE
IMPLEMENTING IT.** "Capacitances in parallel to a noise source" reads as a capacitor across the
SAME node pair. Tested against the authoritative matrix test:

    R from x to gnd, no capacitance at x            bad=True    residual 1.000e+00
    R from x to gnd, capacitor AT x   (strict)      bad=False   residual 0.000e+00
    R from x to y, capacitors to gnd on BOTH        bad=False   residual 0.000e+00   ← no cap ACROSS it

**The third circuit has no capacitor in parallel with the resistor and satisfies the condition
anyway**, because grounded capacitors on both terminals already put the injection direction
`e_x − e_y` inside `im C`. So a strict topological test would report a violation the matrices deny.
⚠ **This is the C-only-loop trap a second time in one day** — a relayed phrasing that, taken
literally, disagrees with a direct computation — and it was caught the same way, by measuring
before building. **The matrix test stays authoritative; what the topological form buys is
LOCALISATION, not a second verdict.**

⚠ Also there: §5 is transient noise simulation of a **ring-oscillator** model with drift-implicit
Euler, with trapezoidal and Milstein variants discussed — the nearest thing in the collection to a
reference implementation of the path-wise route.

⚠ **THE ONE STANDING OPEN ITEM FROM THIS ARC IS AN ACQUISITION, NOT A SEARCH.** Getting Γ
*properly* — a construction rather than the `h → 0` measurement above — needs März, Lamour or
Tischendorf, or a DAE-numerics text. None is among the 339 papers, because DAE projector theory is
not circuit literature. Recorded so nobody spends another pass searching what is already known not
to contain it.


---

### 0k. ⚠⚠ **AN INCONSISTENT OPENING STEP MAKES TRAPEZOIDAL RETURN EXACTLY 2×** — 2026-09-04

⚠⚠⚠ **CORRECTED WITHIN THE HOUR, AND THE CORRECTION IS THE IMPORTANT PART. An earlier version of
this entry — and the commit that shipped it — said "`converged` reports `False`, so this is not
silent". THAT WAS MEASURED ON ODD-STEP GRIDS ONLY AND IS WRONG.** `timestep = per/N` gives `N`
points and `N−1` STEPS, and on an **even** number of steps the 2× answer is **converged AND
periodic to 1e-13**:

    method  pts  steps  parity  converged   |v1|/exact   |v(T)−v(0)|/V
    trap    100    99   odd      False        2.000672      2.0e+00
    trap    101   100   even     TRUE         2.000658      2.7e-13
    trap    400   399   odd      False        2.000041      2.0e+00
    trap    401   400   even     TRUE         2.000041      1.3e-12

**The flag does not merely fail to discriminate — on half the grids it AFFIRMS the wrong answer**,
and a periodicity residual and a refinement study pass alongside it. Whether a caller is warned
depends on the **parity of the point count**, which nobody would think to vary.

⚠⚠ **THE MECHANISM IS AN INCONSISTENT INITIAL VALUE, NOT A BAD MODE** — diagnosed by the docs
session, verified here. The samples are exactly

✅ **AND A 1975 SOURCE SAYS THE 2× IS DISCRETISATION-ONLY** (Trick, Colon & Fan, TCAS 1975,
sensitivities w.r.t. initial conditions — i.e. the monodromy — read firsthand by the docs
session): on a degenerate network *"these additional dependencies are higher order dependencies
which result in derivatives of distributions which **do not affect the 0+ initial conditions**"*,
with their Appendix II deriving the same for Newton applied to oscillators. So nothing about
index-2 structure itself produces the 2×; it is entirely the discretisation inventing an
inconsistent `v(0)`. A stronger statement than "index 2 is where the opener is inconsistent", and
the one to make.

    v_n = V ( cos(ω t_n) − (−1)^n )

so the **smooth part is right** (second order) and a unit ripple rides on it, doubling
peak-to-peak. The index-2 constraint FIXES `v(0) = ωLI`, but the plain path **manufactures** `x(0)`
from the entering state with one order-dropped step and starts the algebraic variable at **zero** —
an error of exactly `−V`. What each method does with that seed is its stability function at the
algebraic limit `|sh| → ∞`:

    trap    (2+sh)/(2−sh) → −1    marginally stable: carried FOREVER at constant amplitude
    euler   1/(1−sh)      →  0    L-stable: killed in one step
    gear-2                →  1/3  killed in a few

That predicts the whole table, **including why the ripple is exactly `V`** rather than an arbitrary
null-space coefficient: it is pinned by `v[0] = 0`. ⚠ And it makes **Euler's `False` HONEST** — it
damps the seed within the period, so its endpoints genuinely differ by that one seed and
`|v(T)−v(0)|/V = 1.0` exactly. Not a false alarm.

✅✅ **AND `x0_unknown=True` FIXES IT AT EVERY PARITY — the remedy is a SHIPPED OPTION.** Measured:

    trap + x0_unknown   1.001343  1.001316  1.000083  1.000082   converged, residual 0.0, v(0)/V ≈ 1
    euler + x0_unknown  0.999832  0.999342  0.999990  0.999959   converged, residual 0.0

Because it makes `x(0)` a genuine unknown instead of manufacturing it — **which is exactly why
Gear-2 was never affected**: its solved-history path already solves for `x(0)`, and it starts
consistent (`v(0)/V = 1.00134`).

⚠ **SO THE CLAIM NARROWS PROPERLY.** Not "trapezoidal is unusable on index-2" but **"the
manufactured opening step is inconsistent on index-2, and `x0_unknown=True` removes it"**. That is
a second and stronger reason for **B1** (`x0_unknown` as the default), which until now rested only
on non-uniform grids.

✅✅ **THE DEFAULT CHANGE IS MADE, AND IT IS CONDITIONAL — 2026-09-04.** `x0_unknown` now defaults
to `None`, meaning *decide from the topology*: switched ON for a netlist the criterion **proves** is
index 2, left OFF otherwise.

    L-I cutset  DEFAULT            1.000082   v(0)/V 0.99996   warns, naming the cutset
    L-I cutset  x0_unknown=False   2.000041   v(0)/V 0.00000   honoured, silent
    index-1 RC  DEFAULT            untouched, `_open_at_x0` False, no warning

⚠⚠ **IT IS NOT A NEW GLOBAL DEFAULT, AND THE REASON IS MEASURED IN `x0_unknown`'s OWN DOCSTRING.**
Trapezoidal still needs an L-stable opener, so switching it on moves the Euler step INSIDE the
period, degrading the ORBIT rather than just the opening — on a `Q = 20` resonator against its
analytic 20 V peak, `x0_unknown` gives **19.76939 against the default's 20.01273** at 100 points.
**Turning it on everywhere would trade a real defect on a few circuits for a real regression on
most.** B1's unconditional form is still open and still needs its own case.

**Three refusals, each deliberate:** an explicit `True`/`False` is honoured untouched; a two-step
method is left alone (its solved-history path already solves for `x(0)`); and ⚠⚠ **a PROVISIONAL
verdict does NOT trigger it — a REFUSAL ON THE THEORY, not caution.** Estevez Schwarz &
Tischendorf's closing page gives up BOTH halves of the criterion for controlled sources: *"if
arbitrary controlling elements for the controlled sources are considered then THE INDEX OF THE
NETWORK EQUATIONS MAY DEPEND ON THE PARAMETERS"*, and *"if controlled sources are allowed to form
a part of L-I cutsets or C-V loops then IT IS POSSIBLE TO BE CONFRONTED WITH HIGHER INDEX (> 2)
PROBLEMS"*. **So `provisional` is not a lower-confidence index-2 verdict — it is not an index-2
verdict at all**, and neither the premise ("the criterion proves index 2") nor the remedy's
justification (`x0_unknown` fixes an inconsistent opening step *on an index-2 algebraic row*)
survives. If the true index is 3 the remedy is not known to apply and would mask a worse problem
while reporting a fix. Same for a structurally singular netlist. It **warns** when it fires,
naming the loop or cutset.

⚠ **AND `topological_index`'s OWN DOCSTRING CLAIMED TOO MUCH — "the DAE index from the netlist
alone", shipped this morning, is FALSE in exactly the controlled-source case.** Corrected to "from
the netlist, WITHIN A STATED CLASS", with both quotes. **The criterion is decidable from topology
only inside its class; outside it, the index can depend on element VALUES.**

⚠ **AND IT BROKE TWO TESTS THAT WERE RIGHT TO BREAK.** `_resolve_x0_unknown` runs BEFORE `solve`
validates its arguments, so a bad `method` reached `_solves_history()` and came back as
`KeyError: 'bogus'` instead of the `ValueError('method must be …')` the caller is owed. **A
defaulting helper has no business changing which exception an invalid call raises**; it is now
best-effort and cannot raise. Pinned by an assertion in the new test.

✅✅ **AND THE CHAIN IS NOW CLOSED AT SOURCE LEVEL** — the one link neither session had read. In
`_traverse`:

    if open_at_x0:                                    # x0_unknown=True
        x = copy(x_in);  x0 = copy(x_in)              # the caller's unknown IS x_0
    else:
        x = self.solve_timestep(x_in, times[0], hs[0])   # ← MANUFACTURED
        x0 = copy(x)

On the default path `x(0)` is one order-dropped Euler step off the entering state. **On this
cutset that step cannot produce the right answer, and the arithmetic says exactly why.** The branch
row gives `v1 = L (i_L(0) − i_L(−h)) / h`; the seed has `i_L(−h) = 0`, and the source contributes
`i_s(0) = I sin(0) = 0`, **so both terms vanish and `v1(0) = 0`** while the exact value is
`ωLI = V`. **The seed error is exactly `−V`, which is why the ripple amplitude is exactly `V`
rather than an arbitrary null-space coefficient.** Measurement, stability argument and source now
agree, and nothing in the diagnosis rests on a relayed claim.

⚠ **AND IT SHOWS WHY THE FIXTURE IS NOT SPECIAL.** Any index-2 netlist whose algebraic variable is
fixed by a *derivative* of the input has this: the manufactured step differences a seed that has no
history, so it returns zero where the constraint demands a nonzero value. `topological_index`
detects exactly that class.

**PINNED** by `test_the_manufactured_opening_step_is_INCONSISTENT_on_an_L_I_cutset`, which asserts
the EVEN-step case is converged-and-wrong (the dangerous one), that the odd-step case fails loudly,
that `x0_unknown` fixes both parities with a consistent `v(0)`, and that gear was never affected.

---

### 0l. Gourary ECCTD 2007 — ⚠ **WE NEVER ASSESSED IT, AND A GENERALISATION IS WHY** — 2026-09-04

Andreas asked why we had not implemented Gourary, Rusakov, Ulyanov, Zharov, Gullapalli & Mulvaney,
*"A numerical technique for time domain noise analysis of oscillators"*, ECCTD 2007. **The answer
is that we never assessed it.**

⚠⚠ **THE RECORD DISMISSES THE AUTHOR WHOLESALE.** C7 rejects *adaptive preconditioning for PAC* as
harmonic balance solving a problem we do not have — sound for C7 — and the B-note explains *"why
all five Gourary papers are about HB"*. **This one is TIME-DOMAIN and says so in its title.** A
blanket claim about a body of work, never measured, used to close a line of inquiry: the same
shape as the MPE and commercial-stepping provenance errors, and it cost us a paper that speaks
directly to A7.

**What the paper addresses**, verbatim: *"numerical difficulties arise because
`J(0) = Φ(T) − I` is singular for oscillator circuits. The singularity of `J(0)` can lead to
ill-conditioning of `J(Δω)` … and subsequently to **distorted PSD curves**."* Its fix uses
`J(Δω) = J(0) + I(1 − e^{jΔωT})` and the LEFT null vector `u` of `J(0)` to substitute the exact
row `uᵀJ(Δω) = (1 − e^{jΔωT}) uᵀ`.

⚠⚠⚠ **RETRACTED — THE ANALYSIS BELOW IS WRONG, AND SO IS EVERY CONCLUSION I DREW FROM IT.** The
docs session obtained the paper and found they had measured a **RECONSTRUCTION** of the method
built from my own one-line description, not the method. The reconstruction substituted the row
`(1 − e^{jΔωT}) uᵀ`, which vanishes at `Δω = 0` — hence "zero row", "no regularisation", and the
"division of labour" conclusion. **The paper DIVIDES THAT FACTOR OUT.** The substituted row is the
CONSTANT `uᵀ`, and the vanishing factor moves to the right-hand side as an explicit scalar:

    uᵀ γ = (uᵀ ρ) / (1 − e^{jΔωT})

`uᵀ` does not vanish. Re-measured on the same `_vdp_at_Q(15.92, 400)` monodromy:

    Δω      cond(J)        cond(Gourary as published)   min row norm
    1e-01   8.588e+00           8.264e+00                  0.962
    1e-03   7.954e+02           8.039e+01                  0.623
    1e-05   7.954e+04           8.079e+01                  0.623
    1e-09   7.954e+08           8.079e+01                  0.623
    0       4.318e+12           8.079e+01                  0.623

⚠⚠ **THE CONDITION NUMBER IS FLAT AT ~80.8 THROUGH `Δω = 0`, against `4.3e12` untransformed —
seven orders at `Δω = 1e-9`. It IS a genuine regularisation and it WORKS AT `Δω = 0`.** The
abstract's claim is exactly right, and my "it hands back a zero row there, so it could never
replace the bordered solve" was wrong.

**THE CORRECTED COMPARISON.** Both are regularisations built on the null vector and they are much
closer than recorded. Bordering augments to `(n+1)×(n+1)`, adding an unknown and a normalisation.
Gourary stays `n×n`, replaces one equation with `uᵀ`, and puts the vanishing factor on the RHS in
closed form. For `Δω ≠ 0` the swap is an exact rescaling of an equation the system already
implies; **at `Δω = 0` it is the continuous extension of that equation, supplying precisely the
information the singular system was missing.** The RHS divergence is not a defect — it is the
phase term's genuine `1/Δω²` PSD, isolated into one closed-form scalar instead of being extracted
from an ill-conditioned solve.

⚠ **WHAT SURVIVES OF THE OLD ANALYSIS:** the identity itself, and that accuracy is capped by `u`'s
own residual (the `1.158e-12` floor, unchanged across every `Δω`). Nothing else.

⚠⚠ **AND THE APPLICABILITY CONCLUSION IS UNAFFECTED, BECAUSE IT NEVER RESTED ON THIS.** Our
near-carrier path holds `S·Δf²` to seven digits down to `Δf/f0 = 1e-9` — measured directly, not
inferred from any characterisation of Gourary's method. **We do not appear to need it.** That a
conclusion survives the retraction of a premise it did not use is worth noticing rather than
assuming.

✅⚠⚠ **ALL NINETEEN NOW READ BY THE DOCS SESSION. THE REGULARISATION IS A FAMILY OF THREE, AND ONE
OF THEM IS THE PAC CASE.** Cited, not measured here:

  * *"A Numerical Technique for Time Domain Noise Analysis of Oscillators"* — the shooting/PPV case,
    the one Andreas asked about;
  * ⚠⚠ *"New numerical technique for cyclostationary noise analysis of oscillators"* — **THE SAME
    TRANSFORMATION APPLIED TO THE PAC FREQUENCY-CONVERSION SYSTEM.** That is **our** harmonic
    singularity, the one `_deflated_solve` handles;
  * *"A New Simulation Technique for Periodic Small-Signal Analysis"* — Krylov specialised to PAC
    under FREQUENCY SWEEPING, i.e. the multi-RHS recycling idea aimed straight at a PAC sweep, which
    is where B6 left matvec counts open.

**Of these, only the cyclostationary one looks capable of changing a design decision, and it would
need measuring against our bordered solve before anyone believes it.** Our own near-carrier
measurement (seven digits at `Δf/f0 = 1e-9`) says we are not currently hurting.

⚠⚠ **AND ONE CLAIM AGAINST OUR RECORD IS DECLINED, BECAUSE IT DESCRIBES A DIFFERENT DOCUMENT.**
The docs session reports that "our document says higher order is unattractive because the
variational recursion must carry extra history per level", and that Obreshkov single-step order-4
voids it. **This roadmap says the opposite** (§B-note): *"the cost is not of parallelism but of
parallelising a method with HISTORY, so a ONE-STEP high-order method pays none of it"*, citing
Chebyshev-IRK for exactly that property. Our recorded objection to adopting one is **specific and
different**: every adjoint path refuses anything but `gear`, so it costs the whole A1–A4d surface —
an ARCHITECTURE cost, which is B8. **Accepting the correction would have replaced a correct
objection with a wrong one.** Second time today a relayed correction did not apply (after MPE).

⚠ **THE PAPER IS STILL USEFUL THOUGH, AND FOR B8 SPECIFICALLY:** *"The Periodic Steady-State
Analysis Based on Single-Step High Order Integration Methods"* supplies the SENSITIVITY-MATRIX
formulas shooting needs, for **charge-oriented** equations, with single-step Obreshkov formulas of
orders 1–4. **B8's whole cost is deriving a reverse recursion for a one-step companion** — this is
material for exactly that, and it should be read before that work starts rather than after.

⚠ **A defect worth knowing if a phase macromodel is ever built here:** Gourary et al. state that
Floquet phase macromodel "variations prevent the determination of the DC solution of the
macromodel differential equation, which in turn prevents its application to standard small-signal,
stability, and noise analyses near the DC operating point"; their smoothed form removes the
oscillatory terms and, for an LC oscillator under sinusoidal excitation, **reduces to ADLER'S
EQUATION** — a closed-form check of the same species as §0g's oracle.

⚠ **AND THE GENERALISATION IS WORSE THAN RECORDED: there are NINETEEN Gourary-authored papers in
the library and only about SIX are HB.** The rest include time-domain oscillator noise,
cyclostationary noise, PSS by single-step high-order (Obreshkov) integration with sensitivity-matrix
formulas, periodic small-signal analysis with Krylov under frequency sweeping, oscillator phase and
frequency transfer functions, stability for large analog circuits, and a PLL/jitter line.

---

**SUPERSEDED ANALYSIS, kept because the failure shape is the useful part:**

⚠⚠⚠ ~~AND IT IS NOT A COMPETITOR TO OUR BORDERED SOLVE — MEASURED BY THE DOCS SESSION ON OUR OWN
`_vdp_at_Q(15.92, 400)` MONODROMY.~~ Because `uᵀΦ(T) = uᵀ`, the substituted row is a **linear
combination of `J`'s own rows**, so the substitution is left-multiplication by an invertible `E`
and `E J x = E f` is **exactly the same system** — same solution set, same dimension, same
singularities. It **selects nothing**. Bordering adds a row AND a column, changes the dimension to
`n+1`, and is nonsingular at `Δω = 0` because it picks one solution from a one-parameter family.

     Δω        ‖uᵀJ − (1−e)uᵀ‖   rel        cond(J)     cond(Gourary)   min row norm
     1e-01        1.1579e-12    1.87e-12   8.588e+00     8.117e+00      6.181e-01
     1e-05        1.1580e-12    1.84e-08   7.954e+04     7.576e+04      6.284e-05
     1e-09        1.1580e-12    1.84e-04   7.954e+08     7.576e+08      6.284e-09
     0            1.1580e-12    1.00e+00   4.318e+12     3.069e+17      0.000e+00

**(1)** the identity is exact — the residual is `1.158e-12` at every `Δω`, which is `u`'s own null
residual and nothing else; **(2)** `cond(Gourary)` tracks `cond(J)` (7576 against 7954, ~5%) and
grows as `1/Δω` just the same — **NO REGULARISATION**; **(3)** ⚠ **at `Δω = 0` the substituted row
is identically zero** and the condition number goes from `4.3e12` to `3.1e17` — **worse**. It makes
the singularity structural and visible rather than curing it.

⚠ **WHAT IT ACTUALLY BUYS (inference from that table, not directly measured):** forming
`uᵀJ = uᵀM − e uᵀ` is a difference of two `O(1)` quantities giving an `O(Δω)` result, so it costs
about `log10(1/Δω)` digits — nine of them at `Δω = 1e-9`. Substituting the closed form recovers
that row to full precision **given `u`**, and is capped by `u`'s own accuracy.

**~~SO THE DIVISION OF LABOUR IS CLEAN~~ (VOID — it works at `Δω = 0`):**

  * **Gourary is for `Δω ≠ 0`**, where `J` is nonsingular in exact arithmetic and the only problem
    is CONDITIONING — cheap, keeps dimension `n`, needs no normalisation;
  * **bordering is for `Δω = 0`**, where the singularity is genuine and a solution must be chosen.

**Gourary CANNOT replace the bordered solve at `Δω = 0` — it hands back a zero row there.** Both
use the left null vector; one supplies an exact ROW, the other closes a RANK DEFICIENCY.

⚠⚠⚠ **RETRACTED WITHIN THE HOUR — I MEASURED THE WRONG QUANTITY. THE NEAR-CARRIER PATH HAS NO
BREAKDOWN AT ALL.** Kundert §3.5 (eq. 15) says the small-signal analysis predicts noise **RISING**
as `Δf → 0` from the carrier. My sweep showed it **FALLING**, and I did not notice that the sign
contradicted the physics. Re-measured at `f = f0 + Δf`, which is the regime Gourary's "distorted
PSD curves at small offset" and Kundert's window are both about:

    Δf/f0     pnoise at f0+Δf   S·Δf²             vs c·f0²
    1e-4      1.251561e-01      3.168626e-09      2.001392
    1e-5      1.248996e+01      3.162131e-09      1.997289
    1e-6      1.249021e+03      3.162194e-09      1.997329
    1e-7      1.249026e+05      3.162208e-09      1.997338
    1e-8      1.249027e+07      3.162209e-09      1.997339
    1e-9      1.249027e+09      3.162209e-09      1.997339

**`S·Δf²` is constant to SEVEN DIGITS over five decades, down to `Δf/f0 = 1e-9`. There is no
breakdown.** What I swept before was `f = f0·10⁻ᵏ` — **BASEBAND**, a different quantity, at signal
levels of `~1e-18` against `~1e+3` near the carrier.

⚠⚠ **SO THE BASEBAND FLOOR IS REAL AND IRRELEVANT.** The invariance study below (reltol, npts, Q,
noisePSD) is a correct measurement of a quantity nobody reads, and its identification with
Gourary's problem — and with "the edge of the practically relevant range" — was **wrong**. Phase
noise is quoted at offsets FROM the carrier, which is the column above.

⚠ **THE SHAPE IS §D 0j AGAIN, AND THIS IS THE FIFTH TIME TODAY:** a clean power law breaking at a
sharp point is convincing enough that I did not ask whether the probe was pointed at the quantity
in question. **`pnoise(pss, f, 0)` at `f ≪ f0` is not the near-carrier noise.** The instrument was
sound and aimed elsewhere.

**WHAT SURVIVES:** Gourary was never assessed (that stands, and the reason stands); it addresses a
real numerical mechanism; and **we do not appear to suffer from it** — the deflated route carries
the pole analytically and holds seven digits at `Δf/f0 = 1e-9`. ✅⚠⚠ **AND THE "LOOSE END" IS CLOSED — IT WAS NEVER A FACTOR OF TWO.** I recorded the ratio
`1.997339` as "the same one-sided/two-sided family §0g already caught once". **It is `A²/2`, the
CARRIER POWER**, and van der Pol's amplitude is 2, so `A²/2 = 1.999` — numerically
indistinguishable from a PSD convention. **§D shape 0i, in my own record, one entry after naming
it.**

`pnoise` returns OUTPUT VOLTAGE noise (V²/Hz); `c f0²/Δf²` is a PHASE PSD (rad²/Hz); the
conversion between them is the carrier power. Measured over a 36× range of `A²/2` — **the ratio
MOVES, which a convention factor could not do:**

    A        A²/2       pnoise/(c f0²/Δf²)   ÷(A²/2)     phase_psd/L
    0.9998   0.49975    0.499334             0.999164    1.00000000
    1.9995   1.99901    1.997338             0.999164    1.00000000
    3.9990   7.99604    7.989351             0.999164    1.00000000
    5.9985   17.99108   17.976041            0.999164    1.00000000

⚠ **`phase_psd` EQUALS KUNDERT eq (15) EXACTLY** — `1.00000000` at every amplitude — so there was
never a convention error to find. And the residual `0.999164` is **discretisation**, converging to
1 as the grid refines: `0.995796 / 0.999164 / 0.999870 / 1.000012` at `npts = 120/240/480/960`.
**Nothing is left unexplained**, which is the standard §0g's own lesson demands.

**PINNED** by `test_pnoise_is_phase_psd_times_the_CARRIER_POWER_not_a_psd_convention`, which
asserts the equality with Kundert at two amplitudes AND that the ratio **scales as `A²`** — a test
that would pass on a constant factor is exactly the test that preserved the last one.

~~APPLICABILITY MEASURED 2026-09-04 — WE DO HAVE A SMALL-OFFSET BREAKDOWN~~ (superseded): `phase_psd` is a CLOSED FORM built from `c` and never touches `J(Δω)`, so A7's
analytic pole-carrying is not the path at issue. The path that IS is **`pnoise` on an autonomous
circuit** → `_deflated_solve`. Swept on `_vdp_at_Q`-style fixture, `f0 = 0.159114371`:

    f/f0      pnoise S(f)       S(f)·f²
    1e-2      1.000041e-10      2.531841e-16
    1e-3      9.998426e-13      2.531340e-20
    1e-4      9.998406e-15      2.531335e-24
    1e-5      9.998478e-17      2.531353e-28     ← clean power law, mantissa 2.5313 to 5 digits
    1e-6      1.071931e-18      2.713848e-32     ← breaks
    1e-7      7.204088e-18      1.823887e-33
    1e-8      7.143311e-16      1.808499e-33     ← RISING
    1e-9      7.215805e-14      1.826853e-33

`S ∝ f²` holds to **five digits** down to `f/f0 ≈ 1e-5`, then breaks, and below that the values
**rise** — the signature of round-off dominating rather than underflow (1e-18 is nowhere near
denormal).

⚠⚠ **AND IT LANDS AT THE EDGE OF THE PRACTICALLY RELEVANT RANGE, WHICH IS WHY IT MATTERS.** Phase
noise is quoted at 1 kHz–10 MHz offsets from a GHz carrier, i.e. `1e-6` to `1e-2` relative. **The
breakdown begins exactly where the useful range ends.** Not academic.

✅⚠⚠ **CAUSATION DISCRIMINATED 2026-09-04 — EVERY ALTERNATIVE IS FALSIFIED BY MEASUREMENT.** The
breakpoint (where `S/f²` leaves its plateau by >0.1%, sampled at 1/3-decade resolution) sits at
**`e = 5.6667` — `Δω/ω₀ ≈ 2.15e-06` — INVARIANT ACROSS EVERY AXIS:**

    reltol     1e-6 → 1e-12   (6 decades)     e = 5.6667 throughout
    noisePSD   1e-2 → 1e-14  (12 decades)     e = 5.6667 throughout
    npts       120 → 480      (16× in h²)     e = 5.6667 throughout
    Q          4 → 32                          e = 5.6667 throughout

⚠ **THE `reltol = 1e-6` ROW IS THE ONE THAT SETTLES IT.** A hardcoded solver floor
(`tol = max(KRYLOV_FACTOR·reltol, 1e-14)`) would be reltol-invariant ONLY once
`KRYLOV_FACTOR·reltol` drops below `1e-14`; at `reltol = 1e-6` the tolerance is genuinely
reltol-driven, and **the breakpoint still does not move.** So it is not the GMRES tolerance.

**What each invariance kills:** `reltol` → not the solver tolerance; `npts` → not the monodromy's
discretisation error; `Q` → not the deflation conditioning (`λ₂ → 1`); `noisePSD` over twelve
decades → not an absolute round-off floor, the limit is RELATIVE.

**What survives is a floor invariant to everything except the offset itself — the signature of
FLOATING-POINT CANCELLATION AT MACHINE EPSILON**, which is Gourary's mechanism and **which no
tolerance can fix**. `log10(1/Δω)` digits lost at `Δω/ω₀ ≈ 2e-6` is ≈ 5.7 of ~16, which is the
right order for a 1e-3 criterion.

⚠⚠ **AND MY OWN FIRST FRAMING OF THE DISCRIMINATOR WAS WRONG, which is why the first run looked
like a double falsification.** I predicted a cancellation floor would MOVE with `M`'s accuracy
(`reltol`, `npts`). **It does not — the cancellation is in the SUBTRACTION, not in `M`**: even with
`M` exact to machine precision, forming `uᵀM − e^{jΔωT}uᵀ` in IEEE arithmetic loses the digits.
The corrected prediction is total invariance, and that is what four sweeps show.

⚠ **STILL NOT PROVEN, and the gap is narrow but real:** the specific subtraction has not been
instrumented. What is established is that **every alternative mechanism I could construct is
falsified by a measured invariance**, which is a much stronger position than the inference this
entry previously recorded — and weaker than a proof.

⚠ **THE ORIGINAL QUESTION IS ANSWERED EITHER WAY.** Gourary was not "considered and rejected"; it
was never assessed, and the thing it fixes turns out to be something we measurably have. If A7's recorded weakness is
ill-conditioning of `J(Δω)` at small offset, Gourary addresses it and our bordered solve does not.
If it is `J(0)` singular, ours addresses it and Gourary does not. **Whether our PAC sweep actually
loses those digits at small offset — A7 already carries the harmonic pole analytically — has not
been measured here.** That is one experiment, not a reading.


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

⚠⚠ **AND THAT ONE GUARD CLOSES A SECOND, UNRELATED HAZARD: THE ITÔ/STRATONOVICH CHOICE.**
Demir ch.2: the Itô SDE `dX = f dt + G dW` and the Stratonovich one agree *"as long as
`G(t,x) = G(t)` is independent of `x`"*; otherwise they are **two distinct Markov processes**
differing *"in the systematic (drift) behavior but not in the fluctuational (diffusion)
behavior"*. `CY = GGᵀ`, so a state-dependent `CY` **is** a state-dependent `G` — and the
cyclostationarity guard refuses it. **The shipped code never faces the interpretation choice**,
so it need not lean on Demir's small-noise resolution (*"the noise signals are small compared
with the deterministic signals"*). Now pinned on all four noise paths.

⚠ **The fixture point is sharper than it looks: for ADDITIVE noise the two interpretations are
IDENTICAL**, so a suite whose sources are all state-independent could not detect an
interpretation error *even in principle*. Ours are all additive. **What makes that safe is not
the fixtures — it is that the code path does not exist**, which is a stronger position than an
untested one and worth distinguishing. ⚠ The tell if it ever arrives: the drift shifts by
`½G∂ₓG` and the diffusion does not, so look for a discrepancy in a **mean** but not a variance.

⚠⚠⚠ **AND THE REFUSAL'S SCOPE IS FAR WIDER THAN "AN UNUSUAL CASE" — IT IS A BLANKET REFUSAL OF
MOS pnoise.** There is **no physically correct MOS noise model whose `CY` is state-independent**:
thermal channel noise is `4kT·γ·g_d0` with `g_d0` bias-dependent, flicker goes as `I_D^AF`, gate
shot noise as `2qI_G`, and Mahmutoglu & Demir (2015) are explicit that trap capture and emission
rates *"depend on the voltages across the MOSFET which can considerably vary with time during
large-signal operation"* — state-dependent **twice over**, since the Wiener process is also
*"modulated with a state (N_t) dependent term"*, which *"in fact makes the equation nonlinear."*

⚠ So the guard is **not a filter that would pass some MOS models and refuse others**. The answer
to "will a real device pass it" is already determined, and it is **no**. Currently invisible only
because `CY ≡ 0` means nothing reaches it. **Same shape as the `TLine` refusal: a scope broader
than its author intended, masked by a feature that does not exist yet.** The refusal message now
says so.

⚠⚠ **THAT REORDERS THIS ROADMAP: the cyclostationary path is not an enhancement for MOS pnoise,
it is the PRECONDITION.** And the literature's answer to a bias-dependent `CY` is not to refuse
it — Hull & Meyer's construction *is* the standard treatment of exactly this case, their worked
example being **shot noise modulated by the collector current**. One source per device at the
cycle-averaged current, modulation in the `H_l` that A1 already computes, under a **checkable
condition** where we currently have a blanket refusal — and it fails in the familiar direction,
degrading as `λ₂ → 1` when the impulse response rings.

⚠ **SECOND-ORDER CONSEQUENCE, WORTH KNOWING BEFORE THE MODEL LANDS.** The same check is what
keeps Itô/Stratonovich out of reach. Relaxing it for MOS makes the interpretations diverge, and
**Demir's small-noise escape may not carry for trap noise** — a trap occupancy is a two-state
Markov chain, not a small perturbation of a large signal, and the same paper calls the state
dependence *nonlinear*. The tell stays: a discrepancy in a **mean** but not a variance.

⚠⚠ **CORRECTED 2026-09-03: `PspMosLongChannel` HAS NOISE. The earlier "identically zero" here
was measured on a DEFAULT-CONSTRUCTED ELEMENT.** The model declares
`white_noise(mult·n_sid) + flicker_noise(mult·n_sfl, ef)` at `compact.py:834`; what is zero is
`fnt` and `nfa`, because *"an element built without a card is noiseless"* (`compact.py:1049`).
With `fnt = 1` the element grows a noise branch (`n`: 4 → 5) and `CY` is nonzero, **exactly
white**, and **bias-dependent**.

⚠ **§D shape 0c, and the fixture was a constructor call.** Every sweep read `CY = 0` and the
conclusion drawn was that the feature did not exist. The model's own docstring lists the noise
under *"Since built, and no longer absent"* and warns two paragraphs later that *"a stale gap
note is worse than none: it is trusted like a measurement and it is not one."* **The note was
current; the reader was not.**

⚠⚠ **SO THE BIAS-DEPENDENT REFUSAL IS REACHABLE FROM A REAL DEVICE TODAY**, and
`modulated=True` is the route past it. MOS pnoise is no longer waiting on a device model.

⚠⚠⚠ **AND THE NOISE HAS AN EXTERNAL ANCHOR — THERMODYNAMICS.** At `Vds = 0` a MOSFET is in
equilibrium, so the fluctuation-dissipation theorem fixes `S_id = 4kT·g_ds` with **no model
freedom**. Measured at `fnt = 1`, `T = 300 K`:

| `Vg` | `g_ds` | `CY[0,0]` | `4kT·g_ds` | ratio |
|---|---|---|---|---|
| 0.40 | 5.221e-05 | 8.895e-25 | 8.649e-25 | 1.0283 |
| 0.80 | 2.180e-04 | 3.691e-24 | 3.612e-24 | 1.0218 |
| 1.20 | 3.392e-04 | 5.711e-24 | 5.620e-24 | 1.0161 |
| 1.50 | 3.919e-04 | 6.579e-24 | 6.492e-24 | 1.0134 |

**Satisfied to 1.3–2.8%**, and the residual is **structural rather than scatter** — it has a
sign and falls monotonically with `Vg`. ⚠ **Deliberately not tuned:** `fnt` is an exact linear
scale (0.509293 / 1.018586 / 2.037172 at `fnt` = 0.5/1/2), so `fnt = 0.98175` would make the
ratio read 1.000000 and would be **fitting a physical constant to a discrepancy we do not
understand**. This is the MOS analogue of `kT/C`.

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

⚠⚠⚠ **THE COST BLOCKER IS LIFTED — HULL & MEYER 1993 NEED ONE SOURCE PER DEVICE, NOT `p` PER
DEVICE.** This entry recorded the barrier as cost: Okumura's construction wants one independent
stationary source per timestep interval per device (~25,000 sources, ~14× PSS per frequency
point). **Hull & Meyer carry the modulation in the RESPONSE instead of in the SOURCES:**

> *"cyclostationary noise sources, such as shot noise, may be modeled as **modulated stationary
> noise sources**. The impulse response that is calculated **includes the effect of this
> modulation**. … the hypothetical stationary noise source has spectral density
> `S_i(ω) = 2q·Ī_c` where `Ī_c = (1/T)∫₀ᵀ I_c(u) du`."*

One stationary source per device, its density using the **cycle-averaged** current, with the
modulation carried by the impulse response `h(u)` we compute anyway. **Same physics, `p` times
cheaper**, and `H_l` is already built.

⚠ **Its condition is stated and checkable:** *"valid when the impulse response duration is much
less than the time it takes for the mixer circuit to significantly change its state … **none of
the large-signal state variables may change significantly over the decay time of the impulse
response**."*

⚠ **Note the regime — it is the OPPOSITE of high-Q, and they say so:** *"high-Q filters should be
avoided, since they cause the impulse response to ring, and thus require a very large value of
M."* So the cheap construction is valid for fast-settling circuits and degrades exactly where
`λ₂ → 1`. **The same boundary as everything else in this record, from a fourth direction** — and
it means the two constructions are complementary rather than redundant: Okumura's expensive one
is what a high-Q circuit needs.

Validated against silicon: a 1 GHz monolithic mixer, predicted vs measured NF over LO power
−12…+6 dBm, *"in good agreement"*, ~**1 dB** error at low LO powers.

⚠ **GENEALOGY, TWO CORRECTIONS TO THIS RECORD'S IMPLICIT CHRONOLOGY.** **Held & Kerr 1978** is
the conversion matrix's source, and its headline is *empirical*: *"correlation of downconverted
components of the time-varying shot noise is shown to explain the **'anomalous' noise observed
in millimeter-wave mixers**"*, validated at 87 and 115 GHz. **Sideband correlation — the fact
the whole cyclostationary line rests on — entered the literature as the explanation of an
experimental anomaly, measured before it was formalised.** And **Rizzoli 1994** had a general
**autonomous** HB noise analysis *"[overcoming] the limitations of the traditional
frequency-conversion approach"* **six years before** Demir/Mehrotra/Roychowdhury — a different
object (frequency-conversion, not Floquet), but this record's implicit story in which oscillator
noise analysis begins with the PPV line is too narrow.

⚠ **CITATION HAZARD ON THE SSB RULE, because Hull & Meyer state it most loosely and theirs is
the reference a mixer designer reaches for:**

| source | condition | verdict |
|---|---|---|
| Ström & Signell 1977 | LPF, all `T_k = 0` above `f_s/2`, **and input band-limited** | *"approximately"* |
| RLF98 | two-sided BPF, `BW < ω₀/2` | **NOT stationary** (`i = 0, ±2` survive) |
| RLF98 | one-sided BPF, `BW < ω₀` | stationary |
| **Hull & Meyer 1993** | **any filter, `BW < f_LO`** | *"stationary"* |

Taken literally the last is too loose — RLF98 Result 1 leaves `i = 0, ±2` for a two-sided filter
just under `ω₀`. **But they are not wrong for their configuration:** an IF filter sits at `f_IF`
far from any half-multiple of `f_LO`, so the `±2` terms have no support overlap and vanish
anyway. **Sufficient in their case, not the general rule.** Our support-overlap test is the
general form; a docstring line should say why the textbook mixer statement looks weaker.

⚠ **AND WHEN IT IS BUILT, THE FIRST GATE IS ALREADY CHOSEN, AND IT IS EXTERNAL.**
Roychowdhury, Long & Feldmann (1998) Fig. 1: stationary noise → mixer(×cos ω₀t) → bandpass at
`f₀`, `BW ≪ f₀` → mixer(×cos ω₀t). A stationary-only path shifts and scales by ¼ twice and
returns **¼** of the input power. The truth is `o(t) = i(t)·cos²(ω₀t)`, so
`⟨o²⟩ = ⟨i²⟩·E[cos⁴] =` **⅜** — *"50% more than that predicted by the previous naïve
analysis."* **1.5× = 1.76 dB, a number that separates**, not a shape agreement. And it satisfies
the discipline rule the `CY/2` error taught: `E[cos⁴] = 0.375` is a two-line time-domain
identity using **none** of the HPSD machinery it gates, so its answer cannot be influenced by
the implementation under test. It is this problem's `kT/C`. (A 4M-sample Monte Carlo of the
cascade gives 0.374918.) A result of ¼ names exactly what was dropped: the cyclostationary
components' contribution to the *stationary* output, which the HPSD route reaches by 1 → 3 → 5
nonzero HPSDs, the extra ¼ landing on the lobe at zero.

⚠ **THE "FILTERING MAKES IT STATIONARY" SHORTCUT IS THREE RULES, AND THE MIDDLE ONE IS A TRAP.**
The support-overlap test on `S_xx,i(ω) = H(−ω)H(ω + i·ω₀)S_nn,i(ω)` — checkable from the
transfer function alone, before any noise analysis runs:

| filter | bandwidth | result |
|---|---|---|
| low-pass | `< ω₀/2` | ⚠ **"approximately" stationary, and TWO conditions — see below** |
| two-sided bandpass | `< ω₀/2` | `i = 0, ±2` survive — **NOT stationary** |
| one-sided (SSB) bandpass | `< ω₀` | **stationary** (RLF98 Result 2) |

A narrow two-sided bandpass is the natural thing to reach for and it does **not** license the
cheap path: **the one-sidedness does the work, not the narrowness** — and the one-sided case
tolerates *twice* the bandwidth. ⚠ Do not implement the shortcut from the one-line version.

⚠⚠ **AND THE LOW-PASS ROW ABOVE WAS ITSELF A ONE-LINE VERSION, CORRECTED 2026-09-03 FROM THE
PRIMARY.** Ström & Signell 1977 p.538 Example 2 states three things the citation compresses
away, and the second one bites:

1. the condition is on **`T_k(f)` for EVERY harmonic `k`**, not on the filter's own passband —
   every harmonic transfer function must vanish above `f_s/2`;
2. ⚠⚠ **the INPUT must ALSO be band-limited to `f_s/2`** — *"although {T_k(f)} are assumed
   ideal, we will obtain FOLDING OF THE INPUT unless also `R_u(f) = 0, |f| > f_s/2`, i.e. the
   input has to be lowpass filtered before transmission"*, restated as *"the necessity of band
   limiting the input u(t) is evident also in this relation"*. **A white noise source violates
   this by construction**, and white is the default assumption everywhere else in this stack;
3. the conclusion is *"**approximately** weakly stationary"*. The paper's alternative route is
   *"regarding the sampling time as random over one sampling interval (or equivalently
   **averaging** the mean and spectral density over one period)"* — a time-averaging argument,
   which is a **weaker and different claim** than genuine stationarity.

⚠ So a precondition test built on the low-pass row must check **the input's band limit as well
as the filter's**, and must not promise stationarity where the paper says "approximately".
Whether (2) binds in practice depends on whether the circuit ahead of the filter has already
band-limited the noise — which is a per-circuit question, not a property of the filter. Same
shape as the two-sided/one-sided trap: a condition asserted without being checked.

#### A4b-note. The invariance any AM/PM split must have — Kärtner 1990 §3.2, and it is a TEST

Under a linear change of state variables `x' = Ax` the Floquet basis transforms as `u' = Au`,
`v'ᵀ = vᵀA⁻¹`, and Kärtner obtains *"the same equation for the time shift θ(t) … therefore the
separation in amplitude and phase is **independent of the co-ordinate system used** … this is by
no means a trivial result, since there are **arbitrarily many other definitions of amplitude and
phase which seem to be more illustrative but do not have this invariance**, and therefore a
change of co-ordinates also **transforms a part of phase noise into amplitude noise** and vice
versa."*

⚠ **THIS IS THE TEST `PAC.am_pm` HAS TO PASS, AND IT DOES NOT HAVE ONE.** Kundert's framing —
*"AM/PM is a change of basis"* on PAC's per-sideband transfer functions — makes Kärtner's §3.2
the property that change of basis must have, and his warning is that most plausible-looking
definitions **do not have it**. A controlled transformation: apply `A`, recompute, require the
phase process **unchanged** and the amplitude process to transform as `dX' = A·dX`. **External
in the strong sense** — a with/without difference whose answer the implementation cannot
influence. Not built; the netlist-level route to a coordinate change is the open part (a
different `refnode` is one, and is the cheapest thing to try first).

⚠ **AND IT IS WHY LEESON FAILS, WHICH JOINS IT TO DEMIR'S COROLLARY 6.1.** Kärtner: Leeson
carries an extra `f⁻ᵅ` term *"due to the amplitude noise, since in Leeson's derivation no
distinction is made between amplitude and phase noise … it is clear that the so-defined
single-sideband phase noise **is no longer independent of the state variable to be measured**."*
**That is Demir's Corollary 6.1 used as a criterion.** A decomposition lacking Kärtner's
invariance yields a "phase noise" that depends on where you probe — exactly the property Demir
proves the correct definition has. ⚠ **The two results are one fact seen from two sides**, which
is what Demir meant by "the same characterization by a completely different derivation" — and
reading both *confirms* it rather than taking his word.

⚠ **AND LEESON IS MORE CAREFUL THAN HIS REPUTATION — he stated the Demir 2006 boundary in
1966.** His own conditions: the RF spectrum equals the two-sided phase spectrum *"subject to the
limitations that Δθ ≪ 1 (small total modulation index) and that AM ≪ FM components"*. **`S_φ`
and the RF spectrum coincide only for small modulation index** — which is exactly the near-
carrier boundary `phase_psd` refuses at. He *conditioned* it; **what Demir added is that the
condition NECESSARILY fails there.**

He also calls it *"a **heuristic** derivation, presented **without formal proof**"* of *"a linear
feedback oscillator"*, and handles upconversion by **inflating the noise figure empirically** —
9 dB *"taken high to account for nonlinear mixing of noise at third harmonic and higher
frequencies which is mixed into the pass band by second harmonic periodic parameter
variations"*, concluding *"the excellent fit of the data implies that this degradation of
effective noise figure may well be an adequate description."*

⚠ **So Kärtner's critique is right and narrower than it sounds: Leeson KNEW about upconversion
and absorbed it into a fitted parameter. Kärtner's advance is that his coupling coefficient —
our `Γ` — PREDICTS what Leeson FITS.** That is the cleanest statement of what A4d buys over the
textbook formula, and it is worth having in those terms rather than as "Leeson is wrong".

What Kärtner captures that Leeson cannot: *"the feedback of the oscillation onto the noise
sources, which results in **multiplicative noise**"*, and *"the mixing and upconversion of noise
due to the non-linearities"* through a coefficient that *"determines **how much of this
low-frequency noise is upconverted** to f₀"* — where in Leeson *"only the noise figure for
small-signal operation of the active element enters, which can hardly describe the discussed
effects."* That upconversion coefficient is our `Γ`.

⚠ **And his numerical route is SHOOTING, in 1990** — *"the so-called shooting methods, which are
based on the fact that the computation of the limit cycle can be formulated as a boundary value
problem."* So the "completely different derivation" differs in its **analysis** (Langevin plus
perturbation methods rather than Floquet DAE machinery) and lands on the same numerical object
computed by the same means. **His `v₁` is what our bordered solve returns.**

#### A4c-note. Demir Ch.6 is the primary for the covariance split — and one open question

Theorem 6.1 / eq. (6.79)–(6.81) give the variance as a **node-dependent prefactor** `ẋ_{s,k}²`
times a **node-independent integral** `(1/T)∫v₁ᵀFFᵀv₁`, hence *"a linear ramp envelope"* with
slope `α`. **That is our `d`, in closed form, as the adjoint quadratic form** — so our two
routes are Demir's derivation and its numerical dual. The split is standard, not improvised.

✅ **Corollary 6.1 HOLDS EXACTLY — and the earlier "110% spread" recorded here was MY
TRANSCRIPTION ERROR, now corrected.** Demir (6.18) defines `S = max_t ẋ_s(t)`, the **maximum
SLEW RATE**; I read it as `max x_s`, the **amplitude**. With `β = 2πf_c/S` and `α = [max ẋ_s²]·c`
the node-dependence cancels identically:

```
β²α  =  (2πf_c)²/(max ẋ_k)² · (max ẋ_k)² · c  =  (2πf_c)² c
```

**MEASURED: rel. deviation 2.22e-16 and 0.00e+00 across the nodes of both fixtures.**

⚠ **The same failure shape as the `v·q = 1` PPV normalisation** — a symbol transcribed as the
plausible quantity rather than the defined one, self-consistent afterwards. Refusing to call the
corollary false on one measurement was right, and *for the stated reason*: the conditions were
unseen, and the unseen thing was in my own reading, not the paper's.

⚠ **And a non-trivial confirmation falls out.** Demir notes that taking the high-to-low
transition *"yields exactly the same results"*. On the **asymmetric** oscillator `max ẋ = 1.767`
and `max(−ẋ) = 2.606` — **47.5% apart** — and the invariant holds for either, because `α`'s
square and `β`'s reciprocal cancel whichever extremum is chosen. That is his mini-invariance
surviving a case where the two branches are nowhere near equal.

⚠⚠ **BUT THE PRECONDITION IS NOT NEAR-SINUSOIDALITY — IT IS THAT `Γ` BE WELL-DEFINED, PER
CIRCUIT VARIABLE.** §6.2.3's examples are a *ring*, a *relaxation* and a *harmonic* oscillator,
so waveform shape is not the issue. What fails is a **triangle wave**: *"ẋ_s(t) for a
triangle-wave is a periodic piecewise constant function, and hence x_s(t) does not have
well-defined low-to-high transition times."* His heuristic is checkable — *"the periodic
waveforms obtained as their time derivatives look like themselves"*, i.e. differentiate and see
whether isolated peaks survive.

⚠⚠⚠ **AND THAT IS A HAZARD INSIDE ANY TWO-NODE TEST BUILT ON THIS.** Demir says only that
every practical oscillator has **a** circuit variable with well-defined `Γ` — *not every*
variable. Point such a test at a node carrying a **ramp** — an integrator output, a relaxation
oscillator's timing capacitor — and `ẋ` is piecewise constant, `Γ` is ill-defined, and **the
test fails with no bug present.** A 0c/0d hazard *inside the gate*. Cheap guard: check that
`ẋ_k` attains its maximum on isolated points. Measured on our fixtures, the max is attained over
**1.9%–3.7%** of the period — isolated, so both are admissible.

⚠ **Its validity boundary is `λ₂` again.** The step to (6.74) drops every Floquet mode but the
first under `|exp(η_i)| ≪ 1` — *"satisfied for 'most' oscillator circuits"* — and Remark 6.2:
*"a second eigenvalue that has a magnitude close to 1 suggests that the oscillator circuit is
close to being unstable, which is usually the case for **high-Q oscillators**."* That is
**`Q ↔ λ₂` written as intuition in 1998**, nineteen years before it became an equality.

⚠ **A disambiguation to keep next to the ISF note below.** (6.79) implies the **node-voltage
variance** peaks at that node's transitions, via the `ẋ_{s,k}²` prefactor. **The PPV does not.**
Both true of the same circuit at once — and *"peaks at transitions"* is precisely the phrase
that would get mis-transferred between them.

⚠ **Kärtner is settled better than the secondary sources put it.** Demir §6.4: his methodology
*"arrives at exactly the same phase noise characterization … even though his definition of phase
noise, and his derivation … is completely different than ours."* Not "Kärtner came first" —
**two completely different derivations, eight years apart, producing the identical formula**,
certified by the later author about his own result.

#### B-note. Why shooting is robust, and what parallel shooting costs

Kundert Ch.8: shooting converges *"if the state-transition function is near linear … it is quite
often the case (usually by design) that the state-transition function is linear even when the
overall circuit behavior is not … **numerical integration is a natural continuation method where
time is the continuation parameter**. This **hiding of the nonlinear behavior** gives shooting
methods a considerable advantage."*

⚠ **That is the mechanism behind two things this record had without an explanation** —
"shooting-Newton needs no preconditioner", and why all five Gourary papers are about HB. The
outer Newton never sees the nonlinearity because time-stepping already continued through it.

⚠⚠⚠ **AND A SECOND, WORSE COST — MEASURED: PARALLEL SHOOTING BIASES `λ₂`, WHICH `Q` THEN
AMPLIFIES.** Kundert ch.4: *"at each step a high-order integration method needs the **history**
of the solution over several past time-steps. **This history cannot extend beyond a shooting
interval boundary.** Thus it is necessary to build up to higher order integration methods by
**taking several steps of a low order method at the beginning of each interval**."*

Low-order steps are what biases `λ₂`, so parallel shooting **injects that bias at every
subinterval boundary**, into the one quantity §0 says everything is bounded by. Measured
(trapezoidal throughout, backward Euler at each boundary, 400 pts/period):

| `Q` | subintervals | low-order fraction | rel err `λ₂` | `Q` error | `Q ×` λ₂err |
|---|---|---|---|---|---|
| 3.18 | 0 (pure trap) | 0% | +2.9e-05 | +0.0% | +0.0% |
| 3.18 | 20 | 5% | −1.34e-03 | −0.4% | −0.4% |
| 3.18 | 100 | 25% | −7.01e-03 | −2.2% | −2.2% |
| 15.92 | 20 | 5% | −2.22e-03 | −3.4% | −3.5% |
| 15.92 | 100 | 25% | −1.11e-02 | **−15.1%** | −17.7% |

Three things: the bias is **negative and linear** in the low-order fraction (≈ −0.028 × fraction
at `Q = 3.18`); **the amplification law holds through it** (predicted vs measured agree at −0.4%
and −2.2%, departing at −15.1% vs −17.7% where the first-order expansion gives out); and ⚠ **the
`λ₂` bias itself grows with `Q`** — −7.0e-3 at `Q = 3.18` against −1.11e-2 at `Q = 15.92` for the
same 100 subintervals — so the `Q` error grows **superlinearly**: 6.9× for a 5× in `Q`.

⚠⚠ **AND THIS IS A GENUINE ASYMMETRY BETWEEN THE TWO BRANCHES, NEW TO THIS RECORD: HB DOES NOT
HAVE THE PROBLEM.** Harmonic balance has **no timestep history to break at an interval
boundary**, so the forced low-order restarts do not exist there and no `λ₂` bias is injected. So
the GPU-parallelism argument is *cleaner on the HB side than on the shooting side* — which is the
first thing in this record that favours HB over shooting on anything but sparsity.

⚠⚠ **THE BIAS IS AVOIDABLE BY CHANGING THE INTEGRATOR, NOT ONLY BY ABANDONING THE PARALLELISM
— and we already hold the note, argued for a different reason.** The cost is not of *parallelism*
but of *parallelising a method with history*, so a **one-step** high-order method pays none of it:
there is no multistep history to break at a boundary. `shooting.py:689` already records Wambacq,
Vandersteen, Phillips, Roychowdhury, Eberle, Yang, Long & Demir arguing for one-step
Chebyshev-IRK precisely because *"each step is independent of the ones before and after"* — from
**stability and step adaptivity**. The parallel-shooting bias is a **second, independent
motivation for the same property, which that paper does not state.**

⚠⚠ **BUT THE ESCAPE COLLIDES WITH OUR OWN ARCHITECTURE, AND THAT IS THE LIVE CONSTRAINT.** We
ship `euler`, `trap` and `gear` (BDF-2). Trapezoidal *is* one-step and would pay no bias — but
**every adjoint path refuses anything but `gear`**: `factored_period().matvec_transposed`, `ppv`,
`_forced_replay_transposed`, `covariance`, `oscillator_covariance` and `_lyapunov_pieces` all
require the solved-history factors. So "switch the integrator to dodge the parallel-shooting
bias" would cost **the entire A1–A4d surface**. Adopting IRK means rebuilding the transposed
replay for a one-step companion, which is a real project and not a flag change.

⚠ **So parallel shooting has THREE costs, not one: memory, the lost nonlinearity-hiding, and a
systematic `λ₂` bias amplified by `Q`.** The parallel-in-time and GPU-shooting papers filed here
as "parallelisation, not correctness" were **mis-filed** — they are parallelisation *at a
correctness price that scales with the thing you care about most*. Not extrapolated past
`Q = 16`; the trend is superlinear and the direction unambiguous.

⚠⚠ **AND THE COST KUNDERT STATES DIRECTLY:** *"the advantage that shooting methods
enjoy by hiding nonlinear behavior from the outer loop is **often lost with parallel shooting
methods**."* The parallel-in-time and GPU-shooting papers are filed here as "parallelisation,
not correctness". **From the method's own author, the trade is not free: multiple/parallel
shooting forfeits the property that makes shooting converge on strongly nonlinear circuits.**
That bounds a whole class of speedups before anyone reaches for one.

Also: shooting *"cannot handle distributed devices"* and lumping them *"considerably increase[s]
the cost"* — the 1990 position our `TLine` refusal reproduces, with Yang & Phillips 2002 as the
dated escape. And eq. (8.1) shows **HB and finite-difference are one method in two bases**:
*"though both methods give the same answer, the matrices in the finite-difference method are
denser."* The basis choice is a sparsity choice.

#### A2-note. The closed-form ISF is a different object from the PPV, and they peak in
#### different places — recorded 2026-09-03 from primaries

Hajimiri & Lee's closed-form ISF is **eq. (36)**, `Γ_i(x) = f'_i / Σ_j f'_j²` — *not* eq. (31),
which is the projection step, and which is the number the secondary literature repeats. Nothing
in this tree cited it; recorded so nothing starts to.

⚠ **THE TWO PAPERS MAKE DIRECTLY OPPOSED CLAIMS ABOUT WHERE THE SENSITIVITY PEAKS**, and this is
sharper than "they differ in shape and magnitude":

| | peaks when |
|---|---|
| Hajimiri, on his eq. (36) | *"maximum during transitions … waveforms with larger slope show a smaller peak"* — **this** node's transitions |
| Srivastava, on the exact PPV | *"the PPV's discontinuities … take place when the oscillator's response is smooth … a node is most sensitive to noise when **the next node in the ring** experiences rapid transitions"* |

**We compute the exact PPV**, so the second is ours. A designer reasoning from the closed form
would place a noise-critical device at exactly the wrong node in a ring.

⚠ **And the closed form assumes HOW the perturbation enters** — eqs. (34)–(36) are derived for
*capacitive node perturbations*, `Δq_i/C_i`. That is the same assumption Andreani identifies as
failing for linear-region transistors: two independent routes to one caveat. Our bordered solve
makes no such assumption.

#### A2-note-2. When the PPV framework itself diverges — and it is not "non-stationary"

Vanassche, Gielen & Sansen (ICCAD 2002) exist to locate the split between Demir/Mehrotra/
Roychowdhury and Hajimiri & Lee. The models differ only in whether the shift is inside the
argument: exact `θ' = ε·Γ(t+θ)·n(t)` versus approximate `θ' = ε·Γ(t)·n(t)`.

> *"for `n(t)` a **stationary** (noise) source, equations (1) and (2) will, up to 0-th order in
> ε, **predict the same output phase noise**. On the other hand, when `n(t)` is no longer
> stationary, results diverge."*

⚠ **The operational form is better than "non-stationary", and it is the sentence to keep:**
*"note that at first, near `t = 0`, the predicted phases are the same. However, **when θ becomes
too large** [they diverge]."* So the failure condition is not the source's stationarity as such
— **it is that `θ` grows large.** A stationary source makes `θ` *diffuse*; a non-stationary one
makes it grow *secularly*, which is what carries it out of range. Same secular growth PPV-HB
splits off as `(Δf/f₀)·t`, and the same unbounded drift Kundert describes: three descriptions,
one mechanism.

⚠ **CONSEQUENCE, AND IT IS A CLEAN BOUNDARY: nothing in our noise path depends on the
distinction**, because every source we support is stationary and the two frameworks then agree
to 0th order. **It becomes load-bearing the moment anything drives the oscillator** — injection
locking, a PLL in lock, coupled oscillators. That is exactly A6, so A6 must use the exact form
and cannot inherit the linearisation the stationary path is allowed.

### A4. Warm start — ✅ **CLOSED 2026-09-06: the automatic criterion is BUILT**

⚠ **This heading's "what remains" is now done.** `PSS.find_initial_solution` implements De Luca,
Bolcato & Schilders Algorithm 2 (commit `a7887cb`). Two things this record had wrong:

- **The paper was NOT unacquired.** It is at
  `~/docs/07-shooting-methods/DeLuca-Bolcato-Schilders-2019-...pdf`. Verify a blocker before
  repeating it.
- **The criterion is not "Algorithm 1 is `_monodromy_matvec`" alone.** Alg. 1 is only the
  Jacobian-vector product. The CRITERION (eqs. 11–13, 16) compares TWO sequences — the LINEAR
  prediction `u_{k+1} = J_phi(x_khat) u_k` against the ACTUAL `utilde_{k+1} = x_{k+1} -
  phi(x_{k+1})` — accepted componentwise for `n_iter` CONSECUTIVE iterations. The guess this
  project had made ("carry a probe and watch it settle") was the wrong SHAPE, not a wrong
  constant: a settled probe only says the Jacobian stopped changing, which is equally true at an
  equilibrium.

Measured: the paper's own RLC (Q=100) gives `khat=0, periods=7` — forced, because an affine phi
makes the linear generator exact; the same tank with a diode refuses the first iterate and finds
the region at 8 periods; and shooting that does NOT converge cold converges from the returned
iterate.

⚠⚠ **NON-AUTONOMOUS ONLY, which is the paper's scope.** For an autonomous oscillator the
equilibrium IS a fixed point of the period map and the map is linear around it, so the criterion
certifies the TRIVIAL ROOT — the van der Pol case in `benchmarks/pss_warm_start.py` is NOT solved
by this and must not be handed to it. That case remains open and belongs with B5.

--- original heading kept below for the record ---

### A4 (original). Warm start — SHIPPED 2026-09-02 as `tstab=`; the *automatic* criterion is what remains

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

⚠ **THE AM/PM NOISE SPLIT IS BUILT 2026-09-06** — `PAC.am_pm_noise(pss, freq, output, carrier)`,
returning `(S_am, S_pm, bands)`. It is NOT `|m_am|^2` from `am_pm`: that is the transfer pair for
a deterministic input, and noise asks whether the two sidebands are CORRELATED.

**The band bookkeeping is the derivation.** `adjoint_sideband_row(pss, g, output, l)` is the
coefficient at output `g + l f0` for a unit source at `g`, so a REAL noise band whose positive
component is at `g = freq + p f0` reaches the UPPER output at `+g` through `l = carrier - p` and
the LOWER at `-g` through `l = carrier + p` (a real process has `N(-g) = conj(N(g))`, and that
shared realisation IS the correlation). Contributions from one band combine coherently, different
`p` sum in power. ⚠ `am_pm` is exactly the `p = 0` term of this sum.

⚠ **THE GATE IS AN IDENTITY, NOT A TOLERANCE.** `pnoise` at the upper sideband folds precisely
those bands and at the lower precisely their negatives, so
`S_am + S_pm == pnoise(carrier*f0 + freq) + pnoise(carrier*f0 - freq)` exactly (no cross term
survives `|a+c|^2 + |a-c|^2`). Measured residual **9.0e-3 -> 3.3e-11** as the sideband count goes
4 -> 64 — it converges away, which a wrong pairing would not.

⚠⚠ **AND THE IDENTITY ALONE IS NOT A SUFFICIENT GATE**: returning HALF the total in each of AM and
PM satisfies it exactly while computing nothing. The test therefore also requires the split to be
NON-DEGENERATE (measured `S_pm/S_am = 1.42`), since equal AM and PM is precisely the
uncorrelated-sideband answer. Both that neuter and a wrong band pairing are verified to fail.

⚠ **The "free LTI invariant (AM = PM exactly)" recorded here is NOT what gates it.** It is the
statement about narrowband noise through a time-INVARIANT system; a circuit with a periodic
operating point has no such limit to take on this path, so the identity above is the gate that was
actually available. Test:
`test_am_pm_noise_splits_the_sideband_pair_and_obeys_its_identity`.

⚠ The autonomous caveat above applies unchanged: on a free-running oscillator the sideband rows
come back at ~1e-12 for reasons not established, so do not read an oscillator AM/PM MAGNITUDE from
this. The driven case is what it is built and gated for.

**Gate:** none needed for AM/PM — it was a basis change on tested output. The others are
interface decisions, not measurements.

### A4c. Time-varying noise statistics — ⚠ BUILT 2026-09-03, DRIVEN **and AUTONOMOUS**

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


#### A4c-osc. The oscillator's covariance — **BUILT 2026-09-03**, `PAC.oscillator_covariance`

`covariance` refuses an oscillator, correctly: `λ₁ = 1` ⇒ `λ₁² = 1`, so `I − M⊗M` is exactly
singular (**σ_min 2.3e-11 against a next singular value of 0.997** — a cleanly *one-dimensional*
null space) and no periodic covariance exists. The answer is not a number, it is a **split**:

```
K(t₀ + nT) = K_orb + n·d·u·uᵀ

[ I − M⊗M    u⊗u ] [ vec(K_orb) ]   [ vec(K₁) ]
[ (v⊗v)ᵀ      0  ] [     d      ] = [    0    ]
```

The border is the pair `ppv()` already computes — `u⊗u` right null, `v⊗v` left null — so this
deflates exactly as A7's PAC solve does. **Bordered σ_min comes back to 5.4e-02: nine orders.**

⚠ **`d` has a closed form and never needed the Kronecker.** Left-multiplying by `(v⊗v)ᵀ` kills
the singular block: `d = (vᵀK₁v)/(v·u)²`, an O(n²) contraction behind an O(n⁴) solve. Both are
computed; they agree to **4e-15**, which makes a disagreement diagnostic (the border pair is
wrong) rather than a precision question.

⚠ **THE GATE IS A PREDICTION, NOT A RESIDUAL.** A bordered system can always be solved. Running
the real Lyapunov recursion forward **forty periods (9,600 steps) from `K = 0`**, touching
nothing the bordered solve produced:

| periods | 1 | 2 | 5 | 10 | 20 | 40 |
|---|---|---|---|---|---|---|
| rel. err | 2.6e-07 | 1.1e-10 | 2.2e-10 | 5.4e-10 | 1.2e-09 | 2.5e-09 |
| walk / total trace | 0.953 | 0.976 | 0.990 | 0.995 | 0.998 | 0.999 |

The first period is the loosest **because of physics** — starting from `K = 0` leaves a
transient in the bounded part, which decays with `|λ₂| = 8.5e-4` and is gone by period two. And
`P(T) − P(0) = d·u uᵀ` to **3.2e-15**.

⚠ **`d/T` IS AN INDEPENDENT ROUTE TO THE DIFFUSION CONSTANT, AND THE ANCHORS ARE INDEPENDENT
TOO** — which is precisely the property that was missing when a 2× error survived a 0.9965
agreement. A phase deviation `α` displaces the state by `α·u`, so the walk is `Var(α)·u uᵀ =
c·t·u uᵀ` and `d = c·T`. `c` is a quadratic form in the **adjoint**-replayed PPV; `d` is a
**forward** Lyapunov recursion closed by a bordered Kronecker solve. `covariance`'s injection is
anchored to **kT/C**; `diffusion_constant` to a **nonlinear Monte Carlo** on zero crossings.

Asserted as a *convergence*, not a tolerance — both carry an O(h) piecewise-constant white-noise
approximation, so what must hold is that the gap halves:

| npts | 120 | 240 | 480 |
|---|---|---|---|
| `(d/T)/c − 1` | 1.87% | 1.03% | 0.54% |

⚠ **`(v·u) = 0.663`, NOT 1, AND ASSUMING OTHERWISE IS THE 2.31×.** `ppv()` normalises on the
first block (`v[:m]·ẋ = 1`) — right for a perturbation entering the first block, which is where
an injected current lands. The full *pair* contraction is a different number, so `(v·u)⁻² = 2.27`.
That exact slip was chased as a code defect. `d` is written with the pair product spelled out.

⚠ **`d` ALONE IS NOT AN INVARIANT.** `u → s·u` sends `d → d/s²`; only `d·u uᵀ` — returned as
`info['growth']` — is a property of the circuit. `u` is pinned by `C u = q`, the same condition
`ppv()` scales the tangent with, which is what gives `d` its `d = cT` reading.

⚠ **And the split is the answer, not a numerical device.** Demir 2002: an oscillator's noise is
*stationary*, not cyclostationary, because "noisy autonomous systems cannot provide a perfect
time reference". `K_orb` is the orbital/amplitude noise a designer can read; `n·d·u uᵀ` is the
walk along the orbit that no periodic object can hold. Reporting only their sum at some finite
time is what makes an oscillator covariance look divergent and useless.

`_lyapunov_pieces` is now **shared** by both routes, so the `CY/2` convention, the `b = 0`
restriction and the `C` ring cannot drift apart the way `diffusion_constant` and `covariance`
once did over exactly that factor of two.

### A4e. Oscillator phase noise — ⚠ CLOSED FORM BUILT 2026-09-03

⚠ **RESOLVED 2026-09-03 — and it WAS a 2× error in shipped code.** `diffusion_constant` used
full `CY` where `covariance` used `CY/2`. Its Monte Carlo validation (0.9965) was against an
injection carrying the *same* hot convention, so it confirmed the bug rather than catching it.
**kT/C**, external to both, settled it: that injection gives **1.92× kT/C** over ten runs.
Fixed; `c` = 7.9516e-08 against a corrected MC at 7.7083e-08, **ratio 1.0316**.

⚠ **The 2.31× was NOT a code defect** — it was the diagnostic contracting a full *pair*
deviation with a `v` normalised on the first block only (`1/(v·u_pair) = 1.508`). Every shipped
path contracts the first block, which is where an injected current lands, so the shipped
normalisation was right. Corrected, the variance ratio goes 2.13 → **1.07**.

⚠ **A SHARED DISCRETISATION IS SHAPE 0b WEARING A RATE.** A convergence rules out a *defect*;
only independent anchors rule out a *shared* one. Two implementations carrying the same O(h)
white-noise approximation would converge together to the wrong answer — the O(h) approximation
is then the instrument, and `kT/C` and the zero-crossing Monte Carlo are what stand outside it.

⚠ **The lesson, recorded in §D:** *a measurement built on the assumption under test cannot
test it.* Two of this session's gates agreed with the code because they shared its convention.

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

### A4d. Coloured noise — ⚠ **BUILT 2026-09-03**; it did cost nothing structural

✅⚠ **2026-09-04 (late): THE FOLD FINALLY HAS A SOURCE TO RUN ON, AND IT IS GATED.** Every library
`CY(x, w)` ignored `w` — no coloured element existed — so A4d's fold had only ever been exercised
on a white source treated as if coloured. `IS` now takes `noiseTau` (Lorentzian,
`psd/(1+(wτ)²)` — exactly white noise through an RC, hence realisable in-netlist) and `noiseFc`
(flicker corner `1 + 2πf_c/|w|`, not realisable by a finite filter; returns the white value at
`w = 0` rather than infinity). **Gate:** the coloured element against the SAME noise realised as
a white `IS` → RC → linear `BSource` into the node — two representations of one physics —

    driven RLC (linear):   coloured vs filtered agree to 1e-5; both track the LTI closed form
                           |Z(f)|²·psd/(1+(2πfτ)²) to 0.2–1.0 % (200-point discretisation of the
                           sideband transfer, identical in both)
    van der Pol, τ=0.3T:   coloured vs filtered agree to 1.3e-4 at Δf/f₀ = 1e-2, 1e-3, 1e-4;
                           PSS periods identical to ten digits (the filter stays out — reconfirmed)

✅ **`pnoise` was correct for colour from the start**: it folds `CY` at the **source-side**
frequency `f − l·f₀` per sideband (A3's design), so a coloured source needs nothing there.

⚠⚠ **`phase_psd`'s `c + Γ(f)` IS NOT THE COLOURED FOLD, and the reason it passed is a fixture
accident.** `diffusion_constant` reads `CY` at ONE frequency (`ω₀`); `Γ` is the `l = 0` term
(square of the mean) at `Δf`. The correct coloured phase diffusion is **per harmonic**:
`c_res(Δf) = Σ_l |V_l|²·CY(Δf − l f₀)/2` — which is what `pnoise` does. Measured on the
Lorentzian-coloured van der Pol at `Δf/f₀ = 1e-3, 1e-4`: `pnoise`-derived, `phase_psd`, `c`, and
`c_res` all agree to ≤ 2e-3 — **because van der Pol's PPV is dominated by `l = ±1`, so the
source-side frequency is ≈ `f₀`, the one frequency `c` reads.** A non-sinusoidal oscillator (strong
`|l| ≥ 2`) or a flicker source (the `l = 0` term at `Δf` dominates) separates them, and `c + Γ`
double-counts `l = 0` (at `ω₀` and at `Δf`). White: `Γ` is 22 orders down, harmless. ✅ **The
build: the harmonic-resolved fold, exact, reducing to `c` for white and containing `Γ` as its
`l = 0` term; `phase_psd` uses it alone; `diffusion_constant` refuses colour, as its docstring
already claims but `_cy_reduced` never enforced.** `_lyapunov_pieces` and `orbital_correlation`
now refuse a coloured source with the reason (they fold `CY` at one frequency; A4d's "plausible
number, not an error" shape). ⚠ Orbital colour — eq (22) with `CY(ω)` — is NOT built; the
collapse to `Ṽᵀ CY Ṽ*` is lost and each term needs a frequency integral.

⚠ **AND A RE-READING OF THIS AFTERNOON'S CORNER CHECK.** At `Δf/f₀ = 1e-2`, `pnoise` exceeds every
phase-only route by 16 % — the same excess the corner check showed (7.51e-8 against 6.25e-8) and
filed as "outside Kundert's window". `1e-2·f₀ ≪ f₀` is inside the window. The better reading is
**the orbital term**: on a Q = 8 van der Pol the amplitude mode's Lorentzian has half-width
`|μ₂|/2π ≈ 0.02 f₀`, peaking exactly there, and `pnoise` is the full sideband transfer — phase
AND orbital — while `phase_psd`, `c` and `c_res` are phase only. **That is the far-out floor A9
exists for, seen at 1e-2.** Candidate, with `S_yy(ω)` as its test; not asserted.

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

#### A4d-built. `colour_projection`, `coloured_diffusion`, `phase_psd`

```
Gamma(f) = <v>ᵀ (CY(2πf)/2) <v>          the SQUARE OF THE MEAN
c        = (1/T)∫ vᵀ (CY/2) v dt          the MEAN OF THE SQUARE
S_φ,i(f) = i² f₀² (c + Gamma(f)) / f²
```

Same vector, same matrix, mean and square exchanged — and **substituting one for the other
returns a plausible non-zero number, not an error.** Measured separation on van der Pol:
`c = 7.95e-08` against `Γ = 1.9e-29`, **22 orders**. They are not approximations of each other.

⚠ **THE GATE IS A PATTERN OF ZEROS, WHICH THE QUADRATIC FUNCTIONAL CANNOT FAKE** because it is
large in every row:

| `a` (even term) | `Rs` (tank loss) | `Γ/c` | `c` |
|---|---|---|---|
| 0.00 | 0.00 | 2.4e-22 | 7.95e-08 |
| 0.00 | 0.20 | 9.7e-23 | 1.01e-07 |
| 0.25 | 0.00 | 4.9e-23 | 8.22e-08 |
| **0.25** | **0.05** | **2.1e-04** | 8.65e-08 |
| **0.25** | **0.20** | **4.1e-03** | 1.20e-07 |

⚠ **TWO INDEPENDENT MECHANISMS FORCE THE ZERO, AND ONLY ONE IS THE ONE DESIGNERS KNOW.**
A *symmetric* waveform gives `⟨v⟩ = 0` — Hajimiri & Lee, why symmetry is the first thing a VCO
designer reaches for. A *lossless LC tank* gives `⟨v⟩[0] = 0` **structurally**, whatever the
waveform does: `v` behaves as `Cᵀv₁` and `dv/dt = Gᵀv₁`, whose inductor row is exactly `v[0]`,
so periodicity of `v[1]` forces `∫v[0] dt = 0`. **A property of the topology, not the orbit.**

⚠ **Which makes van der Pol useless as a POSITIVE fixture** — it reports zero at *every*
asymmetry — **and perfect as a negative one.** A gate built only on van der Pol would have
passed an implementation that returns zero always. §D shape 0c, caught before it shipped this
time rather than after.

⚠ **NO SECOND CONVENTION IS INTRODUCED.** `phase_psd` is checked against `lorentzian`'s far
skirt, already gated by power conservation to 1.000000 — and **the residual is fully accounted
for**, which is stronger than it being small: the Lorentzian carries an `f_h²` term its skirt
drops, so the disagreement must be exactly `(i²·corner/f)²`, quartic in the harmonic.

| harmonic | 1 | 2 | 3 |
|---|---|---|---|
| max \|1−ratio\| | 1.0e-06 | 1.6e-05 | 8.1e-05 |
| ratio | 1 | **16** | **81** = 1 : 2⁴ : 3⁴ |

⚠ **A SLOPE, NOT A STATE — confirmed.** With `CY ∝ 1/f` the skirt steepens to `1/f³` below the
flicker corner and returns to `1/f²` above it, the corner landing where `Γ(f) = c`:

| offset (Hz) | 1e-6 | 1e-5 | 1e-4 | 6.3e-4 | 1.6e-3 | 1e-2 | 6.3e-2 |
|---|---|---|---|---|---|---|---|
| slope | −2.998 | −2.979 | −2.824 | **−2.428** | −2.230 | −2.045 | −2.008 |

`Γ = c` at **6.31e-04 Hz**, the middle of the sweep. A *white* source on the same circuit gives
−2.000000 everywhere, so the steepening is the frequency dependence and nothing else. **The PSS
never sees a filter and no state was added** — Demir's Lorentzian synthesis network is an
artefact of the SDE formulation, and no SDE is formed on this path.

⚠ **`Γ ≤ c` always**, by Cauchy-Schwarz on the weighted mean, and both functionals share a
quadrature *deliberately* so the bound holds exactly at the discrete level — an assertion
rather than an expectation.

⚠ **REFUSED BELOW THE LORENTZIAN CORNER**, a validity boundary rather than a conditioning one:
there the excess phase is a Wiener process with a singular spectrum, and the finite value the
real lineshape attains comes from the **nonlinear** phase-to-voltage map, which
`oscillator_spectrum` carries and this does not.

⚠⚠ **A SECOND FLOOR, ADDED 2026-09-03, AND IT CAUGHT A LIVE DEFECT IN THE TEST ABOVE.** The
normalised lineshape integrates to 1, and one box of width `df` each side is a lower bound on
that integral, so

```
2 · df · S_φ(df)  ≤  1
```

is **necessary** for the linearised skirt. Vanassche, Gielen & Sansen (2003) §6 derive the same
statement for a `1/f` input and reduce it to `df_c ≥ ε·f₀·√(2·f_1f)`; the form asserted here
needs **no assumption about the source's colour** and reproduces their worked example exactly —
**100.000 Hz** against their *"≥ 100 Hz"*.

⚠ **The Lorentzian corner does not catch it, by 306×.** It is built from `c` alone and knows
nothing about a `Γ(f)` that grows as the offset falls. The first version of the `1/f³` test
swept to 1e-6 Hz where **`2f·S_φ = 3.10`** — three times the carrier's total power — and every
assertion in it passed. Two independent floors, and for a coloured source the *new* one binds.

| | value on the flicker fixture |
|---|---|
| Lorentzian corner (white only) | 8.2e-09 Hz |
| **power bound (colour-aware)** | **2.5e-06 Hz** |

⚠⚠ **AND A THIRD FLOOR EXISTS THAT IS NOT A MODEL PROPERTY AT ALL — Kärtner 1990.** He frames
the `1/f` divergence as a dilemma with two exits. *Every* other source we hold takes the first
(postulate a cutoff): Demir 2002 calls it a postulate, Demir 2006 replaces it with a continuum,
Vanassche 2003 shows the chosen value is observable through `|ln γ|`. **Kärtner takes the
second — account for the finite measurement time** — after noting that *"other measurements
show that the 1/f characteristic is conserved down to the microhertz"*, i.e. the cutoff may not
exist. That gives `f > 1/(2T_obs)`.

| floor | depends on | blind to |
|---|---|---|
| Lorentzian corner | `c` alone | colour |
| `2Δf·S_φ ≤ 1` | the spectrum's own values | — (but needs monotonicity) |
| **`f > 1/(2T_obs)`** | **the measurement, not the model** | — inescapable |

⚠ **The third explains why the other two exist**: they are *model-side substitutes for a
measurement-side limit*. A measurement of duration `T_obs` cannot resolve structure below
`1/(2T_obs)` whatever the model says. **If `phase_psd` ever takes an observation window, this is
the floor that belongs to it — and it is the one a user will recognise.** Not built.

⚠⚠ **AND IT INTERLOCKS WITH A DISCRIMINATING EXPONENT THAT IS THE RIGHT NEXT GATE FOR A4d.**
Kärtner: *"in the case α = 1 the phase fluctuations are growing proportionally to **τ²**, in
comparison with the white noise case where the phase fluctuations are proportional to **τ**
… 1/f noise with its infinite correlation time leads to a **quasi-deterministic motion** of the
phase."* The phase does not diffuse, it **drifts**.

**DERIVED HERE, NOT YET MEASURED** — flagged as such deliberately. From `Var[α(t+τ)−α(t)] =
4∫₀^∞ S_α(f)·sin²(πfτ) df`:

| source | `S_α` | `Var(τ)` |
|---|---|---|
| white | `∝ f⁻²` | `∝ τ` — Wiener, and we measured 1.825 → 1793.4 over 1 → 1000 periods, **dead linear** |
| flicker | `∝ f⁻³` | `∝ τ²` |

⚠ **And the `f⁻³` integral is log-divergent at the origin** — `sin²(πfτ) ~ (πfτ)²` makes the
integrand `~1/f` — **which is exactly why the third floor is needed to make the second law
finite.** Kärtner says the same: an *"additional logarithmic part at exactly α = 1 … but this
logarithm is essentially constant in the range τ ≪ 2T"*. **The τ² law and the `1/(2T_obs)` floor
are one construction, not two findings.**

⚠ This is a gate whose wrong answer differs in **slope on a log-log plot**, not in magnitude —
**the shape-vs-scale problem inverted, for once in our favour**, after a `1/f³` slope of −2.998
failed to notice a skirt carrying 3.10× unit power.

⚠⚠ **AND THE BOUND ITSELF HAD AN UNSTATED PRECONDITION — §D shape 0e, one hour after shape 0e
was written up.** The box argument is `2Δf·S(Δf) ≤ ∫_{−Δf}^{+Δf}S ≤ 1`, and the **first**
inequality needs `S(f) ≥ S(Δf)` for every `|f| ≤ Δf`: the spectrum must not dip below its edge
value further in. True of a monotone skirt, of the flattened near-carrier shape, and even with
a spur (which *adds* power inside rather than creating a dip).

⚠ **FALSE FOR A LOCKED PLL**, whose phase-noise transfer function is **high-pass** — suppressed
at DC, rising to the free-running level beyond the loop bandwidth, so it dips below its edge
value everywhere inside. The bound is not thereby shown to be *violated* there (total power is
still 1); it is **no longer derived**, and a floor that is not derived cannot be used as one.
Now **checked on the shape of the returned spectrum**, so it will catch the PLL case when A6
lands without needing to know about loops. *A correction is not self-certifying, and neither is
a generalisation.*

⚠ **A side-finding from building the falsifier:** `diffusion_constant` samples `CY` at the
single frequency `f₀`. That is exactly right for a white source and **meaningless for a coloured
one** — a rising source made `c = 53.3` and drowned every offset below `f₀`, so the guard never
fired and the test passed for the wrong reason. Not a defect today (no shipped source is
coloured) but it is the next thing to get wrong.

⚠ **It is a LOWER BOUND on the breakdown, not the breakdown.** On Vanassche's own example the
observed flattening is at ~300 Hz, **3× the bound**. So it refuses what is definitely invalid
and admits a band that is already suspect — deliberately, since refusing at 3× would be fitting
a threshold to one example.

⚠ **AND THE NEAR-CARRIER FLATTENING NEVER GOES AWAY.** Below the bound the real spectrum
flattens to a DC level with a steep edge, and the position of that edge is set by `γ`, the
low-frequency cutoff of the input `1/f` PSD — the postulate Demir 2002 makes and Demir 2006
replaces with a continuum. Vanassche: *"the DC-level keeps dropping in a manner that becomes
proportional to `|ln γ|` while the corner frequency **increases** like `|ln γ|`"*. Over four
decades of `γ` (5 kHz → 0.5 Hz) `|ln γ|` moves only 8.52 → 0.69. **So pushing the cutoff
arbitrarily low does not remove the flattening; it moves it logarithmically** — the quantitative
form of "`1/f` has no valid stationary model without a cutoff". Note the direction is
counter-intuitive: *as `γ` falls the corner RISES* while the DC level drops.

⚠ **And this is why the paper exists:** Demir's two asymptotes are, in its §1's words,
*"claimed to hold for both white and colored input noise"* — asserted for coloured rather than
derived for it. Our far-skirt 4.00×-per-doubling measurement is untouched (their own scope:
*"for large frequency offsets, the phase noise spectrum assumes the well-known 1/f³–1/f²
characteristic"*); this is strictly near-carrier.

**Still open from this entry:** the sweep-grid trap (a cluster near each harmonic, never *on*
one), the bias-modulated physical trap (which genuinely does belong among the circuit
variables, unlike a synthesised filter), and Mahmutoglu & Demir's whitening under switched bias
— below `f_switch` the real spectrum is flat while a stationary flicker model says `1/f`, so
**the error grows without bound as the offset decreases**. `_Flicker` is a test-only element;
no shipped source is coloured yet.

#### A4d-fold. The phase diffusion is folded PER HARMONIC — ✅ **BUILT 2026-09-04**; `c + Γ` is retired

`phase_psd` computed `S_φ = i²f₀²(c + Γ(f))/f²`, with `c` reading `CY` at ONE frequency (`2π/T`) as if
it held at every harmonic and `Γ` the `l = 0` term at the offset. The correct object is the fold per
harmonic, each read at its own **source-side** frequency — what `pnoise` has done since A3:

    c(f) = Σ_l V_lᴴ (CY(2π|f − l f₀|)/2) V_l,      V_l = Fourier coefficients of the equation-row PPV

`PAC.coloured_diffusion_resolved` is that sum; `phase_psd` uses it. Pinned three ways:

* **Parseval at round-off.** With `CY` constant the sum IS `c`, and with the same step-weighted
  quadrature the discrete identity holds to **4e-16** (asserted at 1e-12). That equality pins the
  transform's normalisation — the one thing a Fourier fold gets wrong silently (by `N`, `T`, `2π`).
* **`pnoise/P_carrier`** on a Lorentzian-coloured van der Pol (`τ = 0.3T`): 1.8e-3 at `Δf/f₀ = 1e-3`,
  2.6e-4 at `1e-4` — the SAME residual as the white control row, so it is the sideband
  discretisation, not the colour.
* **The double count, on the one fixture that can see it.** `Γ` is exactly the `l = 0` term, so
  `c + Γ − c_res == Γ` to round-off. On van der Pol that is `1e-22·c` and proves nothing. On
  `_lc_osc(a=0.25, rs=0.2)` `Γ/c = 4e-3`, so the retired form was **0.4% high for a WHITE source**
  there; the test asserts `Γ/c > 1e-3` first so the identity is not proved on a zero (§D 0m).

⚠ **WHY NO FIXTURE COULD SHOW IT, measured not argued.** Two "asymmetric" cores were built to expose
the double count: `0.25μ(u² − 2)` (`|V₀|/|V₁| = 2.8e-4`) and `0.3u²` (`0.77` by norm, period 6.28 →
6.73 — genuinely bias-sensitive). `Γ` was `1e-27` on BOTH. The `0.77` was the **inductor-current row**:
an inductor to ground shorts the tank node at DC, so a `u²` term moves the PPV's DC component into
`i_L`'s row, which an `IS` at the node never contracts. A flicker source at the tank node of ANY LC
core cannot upconvert through `l = 0`, whatever the core does — the same structural identity
`test_coloured_upconversion_needs_asymmetry_AND_loss` records, met again from the other side. A third
attempt (series `R` after `L`, `0.3u²`) collapsed to `T = 0` with symmetry 1.000 and ratio exactly
2.0 — the constant-PPV limit, a vacuous fixture. §D 0b, three times in one afternoon.

⚠ **THE CORNER IS THE WHITE LORENTZIAN'S, and the first build got that wrong.** Taking the folded
value nearest the carrier as `c` for `f_h = πi²f₀²c` put a `1/f` source's corner ABOVE the offsets
and in front of the power bound — the floor that actually binds for colour — and the existing flicker
test caught it on the first run. `_white_diffusion_at(pss, w)` is the split-out white functional at
one frequency; `phase_psd` reads it at the carrier for the corner as it always did, and
`diffusion_constant` is now that call behind a refusal of colour — the refusal its docstring claimed
since it was written. `oscillator_spectrum` refuses through it (its Lorentzian is exact for white
only, by its own docstring).

⚠ **WHAT IS STILL OPEN HERE (not built):** the orbital term `S_yy(ω)` for a coloured source
(`orbital_correlation` refuses colour). The "14% pnoise deficit" recorded here on 2026-09-04 was
**`c`, not `pnoise`** — the Gear pair's first block is a first-order PPV; see §0l below. ✅ **The vdp
EXCESS is the amplitude mode, pinned by its POSITION (2026-09-05).** `E(Δf) = pnoise/(P_c f₀² c/Δf²) − 1`
is a STEP (a Lorentzian AM sideband over the `1/Δf²` PM part), not a peak; its half-rise corner moves as
`Q_λ^−1.03` over Q = 4…32. A step alone misfits (rms 0.05–0.13, corner 1.4× `f₀/(2πQ_λ)`); a step PLUS a
term linear in `Δf/f₀` fits to rms 0.003 with the corner at 1.12 / 1.06 / 1.04 / 1.02× the prediction
(→ 1 with Q), `E_∞ = 1.01` at every Q (AM = PM far out), and a **Q-independent linear coefficient 1.74 /
1.81 / 1.83 / 1.84** — a `1/Δf` piece of the spectrum. ⚠ **Its attribution to T&B's phase-orbital
correlation term is ARGUED AGAINST by the scaling, the review session's own constraint:** that term's
coefficient goes as `D_lhj = 1/O(μ_l) ∝ Q_λ`, and in `E` it is divided by the phase part `∝ c`; `c` on van
der Pol is exactly flat in Q (6.2537 / 6.2515 / 6.2514 / 6.2504e-8 over Q = 4…32, slope −0.000), so a
correlation-term coefficient would scale as `Q_λ` and the measured one does not. And the term persists
to `Δf = 0.3 f₀`, fifteen corners out, where a phase×amplitude cross term would fall as `1/Δf²`. ✅ **It
is the TANK'S OWN FIRST-ORDER ASYMMETRY, found by the review session's parity test and then derived.**
On the LOWER sideband the coefficient flips sign (upper +1.809 / +1.842, lower −2.197 / −2.162 at Q = 8 /
32): odd in `Δf`, so an asymmetry of the response, which a correction to the even phase-only reference
cannot produce. Its odd part is 2.003 / 2.002, and 2 is closed form: `|Z|² ∝ ω²/(ω² − ω₀²)² =
(1/4k²)(1 + k + …)` at `ω = ω₀(1 + k)` — the upper sideband is `(1 + Δf/f₀)` stronger than the symmetric
leading term, the lower weaker — and with the far-out total twice the phase part (`E_∞ = 1`) the
coefficient in `E` is 2. Q-independent because the tank's asymmetry does not know `μ`. **And the even
remainder is not a term** — the review session's window sweep: with the fit window's lower edge moved
from `Δf/f₀ = 0.001` to 0.15 the odd part stays 2.002 / 2.018 / 2.010 / 2.004 / 1.998 (Q = 8) and
2.003…1.994 (Q = 32), and model-free `(E₊ − E₋)/2k` reads 1.988–1.998 at the outer points, while the
"even" coefficient drifts −0.19 → −0.32 and drags `E_∞` and the corner with it: absorbed step curvature
(a step and a line are not orthogonal on a finite window) plus the tank's own even `−k²/4`, and an even
term linear in `|Δf|` would have been non-analytic at `Δf = 0` anyway. **The excess is fully explained:
the amplitude mode's Lorentzian step plus the resonator's `(1 ± Δf/f₀)` detuning asymmetry showing
through a symmetric reference. Nothing open here.** Test:
`test_the_pnoise_excess_over_phase_only_is_the_amplitude_mode`. (An instrument note: the reference
`c` here is the pair-consistent PPV's; before §0l it was 3e-4 high on this fixture, invisible at this
level.)

#### Hardening pass (measurement-discipline follow-through) — ✅ 2026-09-05

Four guards adopted after the self-confirming-measurement review, all internal (no external-tool data):

* **tnom units guard.** `Behavioural.__init__` warns once when a model's `tnom` looks like a Kelvin
  temperature passed as Celsius (`> 200 C`) — the reversed hazard the Celsius switch introduced, where
  `tnom=300` meaning Kelvin now silently means 300 C = 573 K and moves the drain current ~2x. The
  Kelvin defaults (273, 300, 300.15) all trip it; a real card (≤ ~200 C) does not.
* **Mutation checks.** Three tests inject a defect and assert the gate fires:
  `test_the_kTC_gate_rejects_an_unscaled_CY` (the `CY` vs `CY/2` factor a Monte Carlo once confirmed
  rather than caught) and `test_the_sideband_gate_rejects_the_endpoint_and_the_unconjugated_fold` (the
  two `PAC.solve` reporting defects, on a converting circuit so they bite). A gate that survives its
  own defect is not a gate.
* **Cache-disabled standing test.** `test_cache_disabled_matches_the_cached_compile` compiles the
  constant-folding models with the cache off and asserts bit-identity with the cached compile — the
  stale-constant divergence is now a test, not luck. (The cache key already includes the physical
  constants since 7f5cf33; a rough timing shows a cache HIT costs ~0.44s against a ~0.68s cold compile,
  so a two-level key-on-generated-source scheme is not obviously worth it — A3 is the closure.)
* **Precondition asserts (B1), adopted as convention, and the FULL AUDIT done (2026-09-05).** Every
  gate should assert its precondition (this circuit converts, the variance varies, the fold has a
  phase). The whole of `test_analysis_shooting.py` was audited -- 212 tests: a mechanical triage
  flagged 13 with an assertion but no obvious guard signal and 1 with no assertion at all; hand-reading
  all 14 found the file thoroughly guarded already, in forms a text search does not see -- two-sided
  brackets (`trap > 0.9 Q` with `euler < 0.5 Q`), exact structural equalities (`companion_reach ==`,
  algebraic-row counts), closed-form identities (the Lorentzian integrates to 1; the `20 log10 i` law),
  two-branch detectors with a magnitude floor (`not bad_p and bad_s and res_s > 0.5`), explicit notes
  (`algebraic > 0`, `abs(a-b) > 0.5` for the conjugate), and `assert_array_equal` (the no-assert one).
  The single marginal gate -- `test_gamma_never_exceeds_c_at_the_same_density`, an inequality two
  near-zeros would satisfy -- got a `c > 1e-8` precondition. Plus the three periodic gates from the
  blind-to review (§ committed 68a8d8e). The peer's read of the periodic third generalised: the file
  practises the habit widely.

The general lesson, next to §D: a measurement that shares an assumption with the thing it measures
confirms it at full precision — a cached artefact keyed on source but not constants, a Monte Carlo in
the convention under test, a gate on a circuit whose `v(t)` cannot express the defect. The remedy is a
check that reads the value back from where the assumption did not reach.

#### §0l. The Gear pair's FIRST BLOCK is a first-order PPV — ⚠⚠ **FOUND AND FIXED 2026-09-05**; `c` was 16.6% high on a non-isochronous oscillator and every fixture before it was blind

**How it was found.** The "14% deficit" of `pnoise` under `c` on `vdp + 0.3u²` (bias-sensitive:
period 6.28 → 6.73, `c` 100× van der Pol's). ⚠ **This core is now load-bearing in several gates and has
an amplitude VALIDITY RANGE:** its `0.3u²` term overwhelms the `μ = 0.02` cubic at large excursions, so
under a strong drive or a far-off seed the pre-roll leaves the basin (measured in B16-preroll: `|u|`
1.8 → 92 in three periods under a 0.3 A drive). On its own limit cycle (amplitude 2.2) it is fine; do
not use it as a driven or large-signal fixture. Partitioned by KNOBS before theory:

    grid 400/800/1600:  pn/c 0.855 / 0.924 / 0.961    c 6.264 / 5.800 / 5.581e-6    pnoise 5.355 / 5.361 / 5.363e-6
    fundamental share of P_carrier 0.9555 (reference-bug hypothesis: out)
    Q 4/8/16/32 at 400 pts: deficit 8 / 15 / 26 / 47%  (∝ Q_λ)
    Euler: c and pnoise agree to <1% at every grid  (Gear-specific)

**The reference has no shooting code in it.** The fixture is an explicit ODE; DOP853 at `1e-12`
gives the orbit (period 6.730654; the shooting periods extrapolate to it at second order), the exact
monodromy by the variational equations, its left null vector normalised `v·f = 1`, `v(t)` by
transport: **`c_true = 5.3703e-06`**. `pnoise` at 400 points: 5.355e-06. `c`: 16.6 / 8.0 / 3.9 / 1.9
/ 1.0% high at 400…6400 — clean first order. ⚠ Two kick instruments disagreed with the adjoint first
and BOTH were mine (§D 0j, twice): one subtracted a record-dependent mean before finding zero
crossings (the kicked record spans a non-integer number of periods); the other started the transient
off the orbit, where a non-isochronous oscillator's start-up relaxation shifts the phase by ~0.15 T.
Rebuilt on an exact Poincaré section, the kick matches the adjoint to `1e-4` at 11 of 16 phases (the
rest give exactly one period per kick — an event-count artefact, recognisable on sight).

**Mechanism, then measured.** Gear-2's adjoint state is the pair `(w1, w2) = (∂φ/∂x_k, ∂φ/∂x_{k−1})`.
`w1` alone answers a perturbation of `x_k` with `x_{k−1}` HELD — an inconsistent history, resolved
through the parasitic root `1/3`. Pair biorthogonality makes `w1 ⟂ (u₂,k − u₂,k−1/3)`, i.e. orthogonal
to the amplitude direction ROTATED by `O(h)`; `v·ẋ = 1` then amplifies the rotation by `|v||ẋ|`,
which is the near-cancellation a non-isochronous oscillator has (`|v||ẋ| ≈ 12` here, `≈ 1` on van
der Pol). The physical functional is

    v(t_k) = w1 + Φ(t_{k−1}, t_k)ᵀ w2 = w1 + (C_{k−1} + h G)ᵀ z,   z = −α₂ t_k   (w2 = Cᵀz exactly)

differential rows of `z` only (algebraic multipliers are `O(1/h)`, the decomposition is non-unique
there, and with them in a DC probe flipped sign). **The invariant `v(t)·ẋ(t) = 1` ALONG THE ORBIT is
the test**: first block std `2.7e-2` (mean 1.12), consistent `8e-5`. `c`: `1.8e-3` at 400, `8e-5` at
1600 — second order. Van der Pol unchanged to `1e-4`: its two rows are in quadrature, so the rotation
averaged out of `⟨v²⟩` — every fixture before this one shared the claim's assumption (§D 0b).

⚠ **THE DC CONTENT — a decision taken by measurement, and an open mechanism.** The consistent object
is second order pointwise, but its `O(h²)` pointwise errors do not cancel in the mean: an absolute
floor of `~1e-5 |v|`. The raw block's orbit integral matches a same-grid DC-injection probe to `1e-5`
on the divider fixture (node row, true mean `4e-6 |v|`: below the floor, so the consistent object had
the wrong SIGN at 480 points) — but on the bias fixture's INDUCTOR row (DC voltage in series with L,
`dT/dV = 16.20` by a second-order re-solve) the raw block reads 17.49 / 16.83 / 16.51 at
400/800/1600, first order, 8% off, while the consistent object's invariant pins its mean in every row
to `~1e-5 |v|`. Stitching the raw mean into the consistent object was TRIED and broke the invariant by
±0.3. So `samples` is one object, second order everywhere; the raw pair is kept as
`info['samples_pair']`, and the three structural DC gates read it through `_raw_pair_integrals`, with
the consistent object held within its measured floor beside them. **The inductor-row 8% was the same
defect** — the review session asked for exactly this partition, and with the consistent object the
inductor-row integral reads 16.176 / 16.187 / 16.190 against the re-solved 16.193 / 16.200 / 16.202
(1.0e-3 → 7e-4, both second order). **What remains unexplained is only why the RAW block's integral is
exact to 3e-11 ABSOLUTE on the divider's node rows** where the true mean is `4e-6 |v|` (the consistent
object is `2.7e-6` off there) — recorded as a row-specific identity, not a rule. ⚠ **And the review
session's follow-up partition (a node row with a REAL mean) gave a third answer, 2026-09-05:** on the
asymmetric core with series tank loss (`_lc_osc(0.25, 0.2)` topology, node-`v` mean 0.216 = 44% of the
rms), raw and consistent AGREE to `2e-4` and BOTH sit below the re-solved `dT/di` by `6.9e-4 / 3.7e-4 /
1.9e-4` at 240/480/960 points — **first order, common to both objects** — while the re-solve itself
converges at second order (0.216029 → 0.215976). So (i) the three structural DC gates do not separate
the two objects on any row where the mean is resolvable — they pin the identities they name, not the
object; (ii) there is an `O(h)` DC error in `samples_eq` on this DAE fixture that the ODE fixtures do
not show (the bias fixture's inductor row flattened at `7e-4` with no trend). ✅ **CLOSED the same
day, and the candidate was wrong.** The review session asked for a LINEAR-DAE cell to separate
state-dependence from structure; that cell is empty here (`ppv` is autonomous-only and a linear
autonomous DAE has no limit cycle) — but the fixture's algebraic block of `G` is CONSTANT by inspection,
which killed the state-dependence candidate without a run. The exact reference (the DAE reduced to an
ODE, `x = R·i_L` eliminated; scipy adjoint: `c_true = 1.204953e-07`, `⟨v_v⟩ = 3.137167e-02`) then showed
EVERYTHING on that fixture was first order — `c/c_true` 0.9940 / 0.9971 / 0.9986, the invariant
drifting at `1.1e-3` — not a DC property. The tangent scale was ruled out (an exact differential-block
solve for `ẋ` gives the same `v(0)·ẋ` to `1e-5`). The mechanism: the consistent propagation treated the
algebraic state as a free coordinate and then zeroed it; on a DAE it is SLAVED, and its coupling into
the differential propagation is `O(h)`. The right object is the Schur complement `G[D,NZ] − G[D,Z]
G[A,Z]⁻¹ G[A,NZ]` — the same elimination the algebraic fill performs for the rows. With it `c` is `2.7e-4
/ 7e-5 / 2e-5` from exact, the mean `8.6e-4 / 2e-4 / 5e-5`, the drift `4.3e-4 / 1.1e-4 / 2.7e-5`: second
order on all three. Test: `test_the_pair_consistent_ppv_is_second_order_on_a_DAE_too`. ⚠ **Its
boundary is sharp and named:** `G[A,Z]` nonsingular IS the index-1 condition, so the complement does
not exist at index ≥ 2 (L-I cutset, C-V loop); the review session showed from the pencil that where it
exists it is the exact reduced generator (`eig(−G_red, C[D,NZ])` = the finite generalised eigenvalues
of `(C, G)` to `1e-12`) and that it is undefined on `li_plus_rc` / `cv_plus_rc`. The propagation now
warns once with the reason and falls back to the full `G` (first order there, as the fill already is),
pinned on autonomous index-2 oscillators that `PSS` solves. **The fallback is PRICED on the fixtures
that can see the dropped term** — the non-isochronous core (rows not in quadrature) with each index-2
topology added, so the ODE and its exact `c_true = 5.3703e-06` are unchanged:

    index-1 consistent object            c/c_true  0.99825 / 0.99957 / 0.99992   (400/800/1600)
    L-I cutset (tank inductor split)               0.99822 / 0.99956 / 0.99991
    C-V loop (tank cap split to a bias rail)       0.99853 / 0.99964 / 0.99993

second order on all three, the fallback within `3e-4` of the index-1 object and under its own
discretisation error at every grid. **The guard is the whole answer at index 2, both topologies, and
the März projector chain comes off the roadmap.** (An earlier attempt on plain van der Pol with its
inductor split read 8e-5 and was BLIND — an inert added structure with rows in quadrature can prove
harmlessness and nothing else; the review session named that limitation and both halves of the pricing
fixture.) ⚠ The C-V loop cell first read EMPTY — `NoConvergenceError` at every grid — because the seed
held the bias node at 0 V against a 1 V source; seeded at the source value it converges everywhere with
the index-1 periods to all digits. A fourth partition outcome: "the cell reads empty because the
instrument was mis-set", worth one seed check before recording an empty cell. The review session's
pencil identity stands as the closed-form reason the complement is exact where it exists, and
Weierstrass giving the finite spectrum untouched at index 2 (`ordqz` on `li_plus_rc` / `cv_plus_rc`)
remains the linear reference for driven circuits. The fill's per-sample warning is now one per `ppv`
call. ✅ **And the raw block's
`3e-11` "exactness" on the divider's node rows is CLOSED (2026-09-05): it was a units mix.** The raw
block is the consistent object plus an `O(h)` rotation, and on every row with a resolvable mean the
rotation's DC content follows a SCALE law: `mean(raw) = 1.5 s · mean(consistent)`, `s` the
pair-consistency scale read off the stored second blocks (`w2` unscaled against `w2/s`), exactly 2/3 for
an isochronous pair. Verified to five digits: bias core inductor row 1.081098 vs `1.5 s = 1.08110`
(`s = 0.7207` — the 8%); series-loss tank 0.99858 / 0.99842 vs 0.99842; divider inductor row 1.000143
vs 1.00014. On the divider `s − 2/3 = 9e-5`, so the raw block's RELATIVE DC error is 1.4e-4, and on a
node-`v` mean of 3e-7 that is 4e-11 ABSOLUTE — recorded earlier as an exactness because one number was
absolute and the other relative to `|v|` (the review session's shape 0i, on my side this time). The
consistent object's DC error is ADDITIVE (~1e-6 |v|, the `h Gᵀz` mean) and the raw block's is
MULTIPLICATIVE (`1.5 s − 1`), which is why the raw pair remains the better DC estimator on a tiny-mean
row and the three structural DC gates keep reading it — now for a stated reason. Test:
`test_the_raw_pair_dc_is_the_consistent_dc_times_1p5_s`. **Nothing in §0l is open.** Also structural: the consistent object is ZERO on the algebraic
COLUMNS, as `Cᵀv₁` is (the `hGᵀz` term would leave `4e-3` there; the full suite's Demir-(24) gate
caught it, the targeted subset had not).

**Also fixed on the way, from the peer sessions:** `_cy_reduced`'s third probe was the ZERO VECTOR (on
the orbit by accident; a DC-clocked LTI RC was refused as cyclostationary — a reference-simulator cross-check
suite, 2026-09-05), now the stored mid-period state, with a test; the "7% error" note on the `v·q`
normalisation was `|v·q| − 1`, not the error (`v·q = −1.0696` against `v·ẋ = 1`: a factor 2.07 and
the opposite sign — review audit). The exact vdp `c` at `Q = 8` is `6.250850e-08` (scipy adjoint);
the consistent object is `+2.9e-5` from it, the swept-noise path `−1.4e-4`, and that bound now says so.

✅ **The reference cross-check items are CLOSED (2026-09-05).** `CY` is now evaluated per step at the step's own
state on both Lyapunov paths (`PAC._cy_at`, no cyclostationarity check — the covariance routes model a
modulated source exactly; `pnoise`'s stationary sum still refuses it, as it must), and
`_refuse_coloured` asks its colour question at one state and two frequencies instead of through
`_cy_reduced`, which had refused every modulated source before the covariance could reach it. On the
suite's sample-and-hold to the parameter, covariance over the period in units of `kT/C`:

    points   200      400      800      1600
    hold     0.9920   0.9987   0.9998   1.0000     second order — a reference simulator's sampled pnoise: 0.99999
    track    0.7398   0.8467   0.9157   0.9556     the covariance's known O(h/τ), identical to the control

The held variance is `kT/C` to `1e-4`; the tracking phase sits on the switch-held-closed control at
every grid, so that floor is the separate item it always was. `modulated=True` now says it is for gentle
modulation and fails as a factor, with the suite's table (1.000 / 4.33 / 13.2 / 15.7 / 16.0 over
`goff/gon = 1…1e-6`). Tests: `test_a_switched_capacitor_holds_kTC_with_per_step_CY`; the refusal test
no longer lists `oscillator_covariance`. **Verified against a reference simulator's sampled (time-domain) noise
at matched instants, the whole PROFILE and not only the held number: 0.99878 track, 0.99915 edge, 0.99999
hold** — the transition is what only a per-step `CY` can produce. Two things it corrected in the prose:
the held variance converges at BETTER than second order (6.3× / 7.7× / 13.3× per doubling, the last
against the reference's own floor), and **the tracked variance is 0.957 `kT/C`, not `kT/C`** — a sinusoidal
clock holds the switch at full `gon` only instantaneously, so the capacitor is never in equilibrium
with `Ron`; both tools agree independently, which is worth more than a round number.

**The suite's three further findings, 2026-09-05, all acted on:**

* **`PAC.solve` reported its sidebands wrongly — two defects, neither in the solve (FIXED).** (a) The
  DFT ran over `fp.times`, `[0, T]` INCLUSIVE, so the endpoint repeated the first sample and `dt =
  T/(N−1)` put the sidebands at `f₀(N−1)/N`: 109 500 / 89 500 Hz for 110 000 / 90 000 at N = 200
  (109 875 at 800), costing an ORDER (O(h) for O(h²), 68× at 800 points). (b) `|sb + f|` folded a negative
  sideband frequency to positive with the coefficient untouched; the physical response there is the
  CONJUGATE (uncorrected, `l = −1` was 166% off and did not converge). After the fix the reported
  coefficients equal `adjoint_sideband_row · u_ac` to `1e-15` at every grid with `l = −1` the conjugate —
  the adjoint never calls `freq_analysis`, so it is the proof, and it exonerates the solve. ⚠ Both were
  invisible to every earlier PAC gate: their `v(t)` is CONSTANT over the period (5.6e-16 variation) and a
  constant's DFT is exact for any window — §D 0b, a fixture that could not express the defect. The gate
  now runs on the switched capacitor.
* **A frequency-dependent `CY` could not be swept (FIXED, two layers).** The generated `CY` came out
  RAGGED under an array frequency (a flicker entry array-valued, thermal entries scalar — at `kf = 0`
  too, the term being emitted unconditionally), and on top of that the small-signal analysis handed `CY`
  the whole sweep while the assembly takes scalar entries — which blocked a handwritten coloured `IS`
  just the same. The generated `CY` now broadcasts; `Noise` evaluates `CY` per frequency;
  `dc_steady_state` builds its representative at the first frequency. A coloured `IS` and a level-1 MOS
  sweep to exactly the per-frequency values; three compact models return `(n, n, nf)`. The LIST form
  `freqs=[…]` fails on ANY circuit with a `TypeError` (arrays work) — a pre-existing API quirk, recorded,
  not changed.
* **`oscillator_spectrum`'s `S_v` is exactly 0.5000× a one-sided PSD** (`|X₁|² = A²/4` against the
  carrier power `A²/2`; a reference simulator, four decades). `L_dBc` unaffected. Wording fixed, scale kept: a return
  value callers may already divide by `|X₁|²`.
* **From the same doc, confirmations worth keeping:** a reference simulator's own PPV agrees in scale to 6e-5 and in
  shape to 2.5e-3 with `ppv()` — the `v·ẋ(0) = 1` normalisation settled from outside the project, after
  the pair-consistent contraction landed; `diffusion_constant` reproduces a reference simulator's swept pnoise to four
  digits and `L_dBc` to 0.001 dB over three decades. **High Q is the LIMITING, not the tank loss:** `λ₂ =
  exp(−3bA²T/4C)` has no `g_l` in it, and at `λ₂ = 0.53 / 0.9 / 0.99` the PPV machinery does not degrade
  (multiplier to five digits against the describing function) — what shrinks is the validity window of
  a phase-only spectrum, collapsing onto one curve in `f/f_amp` with `f_amp = −ln λ₂/(2πT)` (−0.03 dB at
  `f_amp/10`, −0.4 dB at `f_amp/3`). `null_residual` stays flat at 1e-9 across it and tells a caller
  nothing; `σ_min` is not in `info`. And a converged autonomous PSS needs a grid that grows with `λ₂`: at
  0.99, 400 and 800 points burn their budget and fail where 1600 converges faster than either, `reltol`
  does nothing, and `pss.converged` is the only thing separating a stall from a solution — the failure
  message sends users to the one knob that provably does not help. **Open, not built.**
* ✅ **`tnom` and Boltzmann's constant, FIXED 2026-09-05 on the user's instruction.** `tnom` on all six HDL
  library models now defaults to the ambient (`float(defaultepar.T)` = 300 K) instead of 300.15 K, so a
  default-constructed device is no longer temperature-scaled by 4.9e-4; `kboltzmann` is the SI-2019 exact
  1.380649e-23 (it was 1.38e-23, 4.7e-4 low — the constant an external reference had to carry as a
  parameter on both sides of every noise test). What moved in the tests, each with its reason in place:
  five compile-record digests (the explain text carries the default), the library3 reference helpers'
  own 300.15 assumption, the 59.5 mV/decade literal (59.53 with the exact `k`), the diode
  series-resistance gate's `rtol` (the junction's log shift is 3.7e-9 of the drop now), and the EKV
  accumulation/seam cards, whose ~1e-25 F entries moved clear of their quantisation and now resolve at
  noise level under an explicit 1e-24 floor (the round-off branch keeps its coverage in library4/5) —
  a reminder that a classification boundary measured in units of a physical constant moves with it.
  ⚠ **And the HDL compile cache did not know about constants:** `vt()` and every noise density
  constant-fold `k` and `q` at compile time, and the cache (13 511 objects) was keyed on source and
  library versions only — so after the change a model read one thermal voltage while the tree's constant
  said another, and only the limiting gate, which reads `VT` back from the compiled spec, noticed (it
  passed with the cache disabled). The physical constants are now in the key. The cross-check analysis names
  the general form, and it belongs next to §D's Monte-Carlo lesson: a cached artefact keyed on source
  but not on the physical constants is a measurement built on the convention under test — it inherits
  the assumption and confirms it — and the gate that caught it did so for the same reason `kT/C` caught
  the `CY/2` factor: it read the value back from somewhere the assumption had not propagated to. Further re-pins with their
  reasons: twelve adoption digests, the limiting gate to two ulp (the C backend and the Python fold now
  round `kT/q` differently), a PCNR tail-node literal, a batched-DC bias pair, a limiter write-back to an
  ulp, and four iteration-count pins that sat on 19 and read 20 or 24 against a budget of 200. ⚠⚠ **And
  three convergence-basin knife edges, exposed by a 4.7e-4 change in the thermal voltage and recorded
  as OPEN solver items, not repaired:** (i) `DC()`'s ladder now hits a SINGULAR Jacobian ('tail' in no
  equation) on the unlimited differential pair at `vin = 0.3` and `1.0` from its default start (the test
  that used it as a reference now uses PCNR); (ii) PCNR stops converging at two of the 48 cascode grid
  points on the 4-terminal level-1 model (`vdd = 20`, `vg2 = 2`, `vg1 = 1.2` and `2.0`; budget 200 and 800
  alike) where the limited plain Newton converges in 45 and 24 — pinned as failures in the test so a fix is
  noticed; (iii) the circuit-level forest at `vin = +1.0` flipped from failing to converging in 18. None of
  these is the constant's fault; each is a basin boundary that a 5e-4 perturbation crosses.
  ⚠ **Still open, two decisions:** `qelectron = 1.602e-19` is the same class of number (1.1e-4 from the
  exact 1.602176634e-19) and was NOT changed — with it exact the slope reads 59.526 mV/decade; and `tnom`
  stays KELVIN where every SPICE card is Celsius (a transcribed `TNOM=27` returns `1.9e92` A silently).

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

⚠⚠ **THE RATIONALE BELOW IS RIGHT ABOUT SALTATION AND WRONG ABOUT WHY — corrected 2026-09-03
from Andreas's domain knowledge plus a measurement.** It rests on the locked orbit sitting at
*exactly zero* phase error, where the PFD never switches. That is the **dead zone**: the interval
in edge *timing* (not voltage, and only picoseconds wide) in which the PFD cannot resolve the
difference, the pulses are too narrow to switch the charge pump fully, `K_d` collapses to zero
and the loop gain with it.

⚠ **But real designs deliberately keep the operating point OUT of it** — a fixed DC current
into the CP output, or an offset between the switched up and down sources, moves the static
phase error away from zero. So along the locked orbit of a *shipped* PLL the PFD **does** switch,
every cycle, with finite pulse widths, and the loop is **not** open in the variational sense.

⚠ **MEASURED, and the pathology does not reproduce.** Two van der Pol oscillators, 1% detuned,
coupled through a narrow zero-gain notch — `|λ| = [1.000129, 0.983971, 0.166, 0.138]`, **one**
unit multiplier, indistinguishable from the same pair under a linear coupling
(`[1, 0.981052, 0.166, 0.136]`). The extra unit multiplier needs the orbit to *live inside* the
notch, which is the unmitigated-dead-zone case — a broken design, not the normal one.

⚠⚠⚠ **AND THE REPLACEMENT REASON IS FALSIFIED TOO — MEASURED. SALTATION IS NOT NEEDED FOR A
PFD AT ALL.** I proposed that saltation is still mandatory because the PFD switches twice per
cycle where the field is discontinuous. C5 measured a switched *conductance* and found the
monodromy-vs-finite-difference gap falling at exactly **2.00× per doubling** — O(h), not the
O(1) a missing term leaves. A PFD is a discontinuous **injection**, not a conductance, so C5's
result had to be re-run with that element and everything else held:

| element | 200 | 400 | 800 | rate |
|---|---|---|---|---|
| `VSwitch` (conductance) | 9.74e-04 | 4.86e-04 | 2.43e-04 | 2.01×, 2.00× |
| sign source (**injection**) | 7.54e-05 | 3.76e-05 | 1.88e-05 | 2.01×, 2.00× |

Both toggle twice per period with `|M|` = 0.653 and 0.990, so neither is the numerical-zero trap
C5's own docstring records falling into. **C5's reason turns out to be general**: each step uses
its own converged `Jf` and `C`, which already describe whichever side of the switch that step is
on, and a *discrete* map has no instant at which the field is undefined.

⚠ **So A6's saltation claim was wrong twice over — the stated reason and my replacement for it.**
The stated reason described an unmitigated dead zone; the replacement described a
continuous-time construct a discrete map does not need. Neither survived a measurement, and the
second was falsified by the instrument that had already settled the first.

⚠⚠ **THE DIVIDER WAS THE OPEN QUESTION AND IT IS NOW MEASURED — saltation is not needed there
either, and the real problem is something else entirely.**

⚠ First, a correction to this record: **`Idtmod` folds the STATE, not only the output map.**
`I.idt_node` stays inside one modulus and is periodic to 2.6e-15. So it *is* the state-reset
object the question needed.

**Off a grid point the flow map is differentiable and the gap falls at O(h)** — wrap placed
off-grid with `ic = 0.31`:

| npts | rel err | rate | FD noise floor |
|---|---|---|---|
| 250 | 2.754e-06 | — | 7.7e-08 |
| 500 | 1.349e-06 | **2.04×** | 1.3e-07 |
| 1000 | 6.738e-07 | **2.00×** | 2.7e-07 |
| 2000 | 5.551e-07 | 1.21× | **5.5e-07** ← floor |

Exactly first order, same as the switched conductance and the discontinuous injection. The last
row's 1.21× is the FD instrument's **own noise meeting the signal**, which is why the
eps-stability is measured alongside rather than assumed away.

⚠⚠⚠ **BUT ON A GRID POINT THE MAP IS GENUINELY DISCONTINUOUS, AND THE PSS CONVERGES ANYWAY.**

| ε | npts=600 (off) | npts=1200 (**on** a grid point) |
|---|---|---|
| 1e-10 | 1.732085 | 8.202463e+07 |
| 1e-08 | 1.732026 | 8.202453e+05 |
| 1e-06 | 1.732025 | 8.201463e+03 |
| 1e-05 | 1.732025 | 8.192475e+02 |

Constant across six decades on the left; on the right `|Δφ|` is a **constant ≈8.2e-3 independent
of `ε`** — a perturbation of *any* size gives the same finite jump, because an infinitesimal
change flips which step the reset lands in and a whole modulus propagates.

⚠ **SO A6'S REAL PROBLEM IS EVENT LOCALISATION, NOT SALTATION.** `shooting.py` contains **no
reference to `next_event`** — the consumer is `transient.py` alone. So the transient breaks its
steps at events and the PSS traversal does not.

✅ **HALF OF THAT IS BUILT 2026-09-06 as `PSS.event_grid(period, npts=... | grid=...)`**, which
returns step fractions with every event in the period landed ON a grid point. Measured on an RC
driven by a `VPulse`, against a 4000-point reference, a 40-step uniform grid versus the same grid
with its 3-4 event times landed:

    edge offset    uniform      + events     gain
    td = 0         6.787e-03    8.032e-04    8.4x
    td = 0.0125T   5.720e-03    2.099e-04     27x
    td = 0.0092T   2.483e-03    7.910e-04    3.1x

⚠ **THE SNAP IS WHAT KEEPS IT SAFE:** an event near an existing point MOVES that point onto it
rather than inserting a second beside it, so no sliver is created (smallest resulting step: 8% of a
uniform one). Inserting unconditionally is how a merge acquires arbitrarily small steps — the same
lesson B7c's separation rule encodes.

⚠ **It is a HELPER in the `lte_grid`/`refine_grid` idiom, not a default change.** Making event
breaking automatic would alter step selection for every PSS solve in the suite; that is a separate
decision and has not been taken.

⚠⚠ **AND ONLY THE TIME-DRIVEN HALF IS SOLVED.** `next_event(t)` is parameterised by time, so a
source's edges can be walked out once and placed — which is why this needed no change to the six
traversal loops, and is identical to shortening steps in flight for those events. A STATE-DEPENDENT
reset cannot be: `Idtmod.next_event` is a linear prediction from the last accepted point and returns
`inf` before a traversal starts, and its wrap time MOVES as the Newton iterates. **That half remains
exactly as described below — it needs the event time to become an unknown the Newton solves for**,
and `event_grid` will not help it.
Test: `test_event_grid_lands_the_period_on_its_event_times`.

⚠⚠ **BUT "AND THE PSS REPORTS CONVERGENCE THERE" WAS WRONG, and the correction matters more
than the claim.** It warns, loudly, three times over on this fixture:

* *"Local truncation error reaches **4.58e+05 times tolerance** accumulated over the period"* —
  the LTE check already catches the under-resolved wrap;
* *"the returned waveform **IS STILL A FULL RESULT** — it is the last iterate, not a periodic
  steady state — so a reader who does not check `converged` gets an array that looks like an
  answer and is not"*;
* and a third that caught a defect in **my own fixture** — see below.

So the existing machinery does guard this. The gap is that the grid cannot *break* at the wrap,
not that the failure is silent.

⚠⚠⚠ **AND THE THIRD WARNING FOUND A FIXTURE DEFECT: the divider fixture is solved as AUTONOMOUS
at 2× the fundamental.** *"this autonomous solve returned a period that is a MULTIPLE of the
fundamental … the fundamental is about 5.0e-4 s and the returned 1e-3 s traverses it about 2.0
times."* `VS` is DC-only, so there is no periodic drive and PSS infers autonomy — the same trap
as the MOS amp at `va = 0`.

**What that does and does not invalidate.** The monodromy-vs-finite-difference **rate** (2.00×)
is an internal consistency of one traversal — same grid, same period on both sides — so it
stands. What it does invalidate is the separate comparison of the PSS *solution* against a
settled adaptive transient: those were a `k·T` orbit and a `T` orbit, i.e. different objects,
and the 1e-2 "error" measured there is not evidence of anything. ⚠ **A fixture running at
4.58e+05× its LTE tolerance is not one to draw quantitative conclusions from without saying
so**, and the earlier entry did not say so.

⚠ **What the dead zone IS still worth to this roadmap:** an unmitigated loop, or one whose
offset is too small, genuinely does sit where the linearisation describes a disconnected
circuit. That is a real failure mode to be able to *detect* — and `ppv()`'s
`PPV_SECOND_MULTIPLIER_WARN` would fire on it, which is the right behaviour arrived at for an
unrelated reason.

⚠ **The literature's original claim, kept for the record.** Along a *locked* orbit the reference
and divided edges coincide, the PFD never changes state, and "the analog part of the circuit
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

⚠ **A NAMED RESIDUAL WEAKNESS — see §0f item 2.** Carrying the harmonic pole analytically is
DEFLATION, and Mei & Roychowdhury's Lemma 2.1 is relayed as proving the singularity at every
harmonic and then saying deflation only **ameliorates** it. Known, not fixed, no change proposed.


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

### B6. Floquet multipliers from the GMRES basis — ⚠ **HALF BUILT, HALF CLOSED 2026-09-03**

⚠⚠ **THE TOLERANCE FLOOR HAS A PUBLISHED EXPLANATION — see §0f item 1.** Mei & Roychowdhury 2007
are relayed as saying Krylov partial decomposition **exacerbates** the imperfect cancellation that
deflation leaves. That is the floor measured below (`reltol = 1e-12` exhausting the inner restarts
at `Q ≥ 16`, 126 matvecs, against 28–35 flat in `Q` at `1e-9`) — **not a bug, a known boundary**.
Do not re-litigate it as one.


García, Romero & Acha (IEEE Trans. Power Systems 37(1), 2022) determine periodic-orbit
stability *"by computing the Floquet multipliers using **Ritz values and the Hessenberg
matrix**"* of the GMRES that already solved the Newton correction. **We build that Hessenberg
matrix every shooting-Newton solve.** Ritz values `θᵢ` of `I − M` map back as `λᵢ = 1 − θᵢ`, so
the multipliers are free: **no eigendecomposition, no monodromy formed** — which also sidesteps
the eigenvector-selection problem Demir & Roychowdhury 2003 retire the monodromy method over.
(They add a Givens-rotation QR of `H` worth up to 50% more, and report up to **8×** against the
standard Poincaré-map method on a 118-node system.)

⚠ **THE MAPPING IS VERIFIED ON OUR REAL MONODROMIES** — machine precision at `k` well below the
full dimension:

| fixture | true spectrum | k=2 err(λ₂) | k=3 | k=4 |
|---|---|---|---|---|
| Q=0.12 | {1, 2.9e-4, 0, 0} | 3.88e-04 | 1.3e-14 | 7.9e-16 |
| Q=15.9 | {1, 0.939, 0, 0} | 6.69e-02 | 3.7e-15 | 3.7e-15 |
| Q=63.7 | {1, 0.984, 0, 0} | 1.71e-02 | 3.2e-14 | 1.1e-15 |
| slow τ/T=1e4 (n=6) | {1, 0.9999, 8.6e-4, 0,0,0} | 1.16e-04 | 1.04e-04 | 2.2e-15 |

⚠⚠ **BUT THE CLAIMED `Q`-IMPROVEMENT DOES NOT REPRODUCE HERE, AND THE REASON IS SHAPE 0c.** The
docs session measured `λ₂`'s Ritz error *improving* with `Q` on a 40-dimensional synthetic
monodromy with a damped **bulk** of 38 eigenvalues in `U(0, 0.35)`. Their mechanism is
**separation**: at low `Q`, `λ₂` is buried in the bulk and Arnoldi cannot pick it out; at high
`Q` it is pulled away and becomes isolated.

**Our fixtures have no bulk.** The spectrum is `{1, λ₂, 0, 0}`, so `λ₂` is *always* isolated and
the mechanism cannot operate. At `k = 2` we measure the **opposite** direction — and that has its
own explanation: with no bulk to hide in, a high-`Q` `λ₂ ≈ 1` instead competes with **`λ₁ = 1`**
for the same Krylov direction. **That is door (a), numerical distinguishability, appearing inside
the Arnoldi basis.**

⚠ **So the synthesis is more informative than either measurement: whether the Ritz route improves
or degrades with `Q` depends on WHAT `λ₂` is competing with.** Against a damped bulk it improves
(separation grows); against `λ₁ = 1` it degrades (distinguishability shrinks). A real circuit has
both — a bulk of fast modes *and* `λ₁ = 1` — so **the gate must be run on a circuit with a
genuine damped bulk before this is adopted**, and neither of our fixtures qualifies. Not their
result refuted; **untestable on what we have**, which is exactly §D 0c.

⚠⚠⚠ **THE GATE WAS RUN — on a synthetic `{1, λ₂, U(0, bulk_max)^38}` carrying BOTH a bulk and
`λ₁`. Both effects are real, they cross over, and THEY ARE NOT COMPARABLE IN SIZE.** Median Ritz
error in `λ₂` over 12 draws at `k = 8`, `bulk_max = 0.35`:

| `Q` | 0.14 | 1 | 2 | 8 | 32 | 64 | 256 |
|---|---|---|---|---|---|---|---|
| err(`λ₂`) | 8.9e-03 | 2.9e-02 | 2.1e-04 | **2.5e-05** | 9.0e-05 | 2.1e-04 | 2.8e-04 |

**Non-monotonic, sweet spot at `Q ≈ 4–16`.**

⚠ **The lower edge is a CLIFF and it is predictable:** the error drops immediately above
`Q = −1/ln(bulk_max)` — verified at 0.53, 0.95 and 2.80 for `bulk_max` = 0.15, 0.35, 0.70.
**`λ₂` must be slower than the fastest of the fast modes before Arnoldi can separate it.**

⚠⚠ **AND THE TWO EFFECTS ARE VERY UNEQUAL, WHICH SETTLES IT:** the bulk cliff is **3–4 orders
over a factor ~2 in `Q`**; the `λ₁` degradation is **1–1.5 orders over four decades**. So the
Ritz route is governed by the **bulk** question, not the `λ₁` collision. My concern was real and
is **second-order**. The honest statement is neither "improves with `Q`" nor "degrades": it
**improves sharply up to the bulk edge, then degrades gently forever**.

⚠⚠⚠ **AND THE FINDING NEITHER OF US ANTICIPATED — THE BULK WIDTH DOMINATES `Q` ENTIRELY.** At
`bulk_max = 0.15` the error is ~1e-7 across the whole usable range; at 0.70 it never gets below
1e-2. **Five orders, from the PARASITICS rather than from the oscillator.** A circuit with fast,
well-separated parasitic modes gets an essentially free `λ₂`; one whose fast modes are sluggish
does not, **whatever its `Q`**. That inverts the usual reading of §0: here the thing that decides
is not the designed quantity but the incidental one.

**Gate remaining:** confirm on a real circuit with ≥10 states — the synthetic settles the
mechanism, not the applicability. Cheap; the mapping is already verified above.

#### B6-note. Matrix-free shooting at high `Q` has a TOLERANCE FLOOR — bisected 2026-09-04

Two sessions measured `matrix_free=True` on `_vdp_at_Q` and got opposite answers — "cost
saturates with `Q` at 28–35 matvecs" against "fails above `Q ≈ 3`". **Both are right; the
variable is `reltol`, and the circuit was identical.**

| `Q` | `reltol = 1e-9` | `reltol = 1e-12` |
|---|---|---|
| 3.18 | 28, converged | 28–39, converged |
| 16 | 28–35, converged | **126, FAILED** |
| 64 | 28–35, converged | **126–127, FAILED** |

⚠ **The failure is not an outer-iteration budget** — `maxiterations` 150 and 600 both give
exactly 126 matvecs. The inner GMRES exhausts its own restarts.

⚠⚠ **So the honest statement is neither "saturates" nor "fails": matrix-free shooting at high
`Q` has an ATTAINABLE-TOLERANCE FLOOR set by `λ₂`.** Below it the cost is flat in `Q`; demand
more accuracy and it stops converging. That is door (b) — conditioning — quantified, and it is
the first measurement in this record of what `λ₂ → 1` costs the *solver* rather than the
*answer*.

#### B6 outcome: the Ritz values are IN, reusing the Newton's basis is OUT

✅ **The valuable half shipped** (`fef3d60`): `ppv()` now takes `λ₂` from **Arnoldi Ritz values**
instead of a deflated power iteration — machine precision even at exact degeneracy, at `k`
matvecs instead of 30. That was the accuracy win, and it is independent of where the basis
comes from.

⚠ **The remaining half — reusing the Newton's Hessenberg matrix — is closed, on measurement:**

1. **The default Newton has no GMRES at all.** `solve(matrix_free=False)` is the default and
   factors the Jacobian directly. Instrumenting the `m = 12` bulk fixture counted **zero** inner
   GMRES matvecs. García's route presupposes a Krylov shooting-Newton; ours is dense.
2. **Where there is one, the saving is 12 matvecs** — 27% of `ppv()` but only **~1.5%** of a
   PSS-plus-`ppv` workflow (Arnoldi 0.077 s against PSS 4.17 s + `ppv` 0.29 s), and *shrinking*
   on the large circuits where `matrix_free` is actually worth using, because the PSS dominates
   more there.
3. **scipy's `gmres` does not expose `H`**, so realising it means replacing the core Newton's
   inner solve — the highest-risk change available — for that 1.5%.
4. ⚠ **And nothing depends on the choice.** `factored_period()` re-traverses at the *converged*
   solution, so `λ₂`, `Q` and the PPV are **bit-identical** under `matrix_free=True` and
   `False`. Now pinned by `test_the_ppv_is_invariant_to_the_newtons_inner_solver`, because it is
   a fragile property: the class docstring records that `Jtvec`/`Cvec` are written by *neither*
   factored traversal, so an analysis reading them after a matrix-free solve would rebuild an
   operator for a different trajectory, silently.

⚠⚠ **AND THE STRUCTURAL REASON, WHICH IS STRONGER THAN ANY OF THE FOUR ABOVE AND WAS FOUND
LAST.** Even on the matrix-free path, *which operator GMRES solves depends on the case*:

| case | matrix-free operator | Ritz values give |
|---|---|---|
| **driven** | `v ↦ v − α·M v` — plain **`I − M`** | `λ = 1 − θ`, the multipliers directly ✓ |
| **autonomous** | `[v; s] ↦ [(I−M)v − s·dφ/dT ; v_k]` — **bordered** | eigenvalues of a bordered operator, **not** `1 − λᵢ` |

**The two halves do not meet.** Where the operator is clean (driven) we do not need the
multipliers — `ppv()` refuses driven circuits, and `Q`, `λ₂` and the PPV are all oscillator
quantities. Where we *need* them (autonomous) the operator is bordered, **and the bordering
exists precisely to remove the `λ = 1` singularity** — i.e. to destroy the spectrum one would be
trying to read.

That is a structural mismatch rather than a cost argument, and it survives every improvement to
the cost side. ⚠ An earlier statement of this omitted the driven exception and claimed the
bordering blanket-wide; the driven operator really is plain `I − M`.

⚠ **Undoing the bordering analytically** to recover `I − M`'s spectrum from the bordered `H` is
possible in principle, more work than the 12 matvecs it saves, and error-prone in the one place
where a wrong `λ₂` is amplified by `Q`.

⚠ **What would re-open it:** a hand-written Arnoldi-GMRES for the Newton would give `H` for
free *and* retire `_gmres_checked`. ✅ **That half is now BUILT** (`0c9f0fe`): `_arnoldi_gmres`
keeps `H` and judges by residual, and is swapped into `_gmres_checked` and `ppv()`'s two
bordered solves. So "scipy does not expose `H`" is no longer a reason — the remaining ones,
especially the structural mismatch above, are.

### B7. Adaptive time stepping in the inner transient — ⚠ **GATE RUN 2026-09-04, IT FAILS; SPLIT IN TWO**

`PSS` builds a **fixed** grid with `_period_grid(T, npts, fracs)` and traverses it; the
`Transient` it is built on is adaptive and breaks its steps at events. So the shooting solve
throws away step control it already owns.

⚠ **THIS IS THE SAME ITEM AS A6's REAL DEFECT, ARRIVING FROM THE OTHER SIDE.** `shooting.py`
contains **no reference to `next_event`** — the only consumer is `transient.py`. Measured
consequence on a wrapping `Idtmod`: `LTE 4.58e+05 ×` tolerance, and a monodromy that is
*discontinuous* when the reset lands on a grid point (`|Δφ|` constant in `ε` over four decades).

⚠ **The hard part is not the stepping, it is that shooting needs a REPRODUCIBLE grid.** The
monodromy is a product of per-step maps, and `factored_period()` re-traverses at the converged
solution: if the grid moves between traversals the map is not the one the Newton converged on.
So an adaptive scheme has to be *frozen* after the first pass, or made a function of the state
only.

---

⚠⚠ **THE GATE WAS RUN 2026-09-04 AND IT FAILS. The entry above is left as written because what it
got wrong is the useful part.**

**FIRST, MOST OF THE MECHANISM ALREADY EXISTS AND THE ENTRY DID NOT SAY SO.** `PSS.solve(grid=…)`
already takes non-uniform step FRACTIONS and freezes them, and `benchmarks/pss_lte_grid.py`
already derives such a grid from an adaptive `Transient` — measured on van der Pol at `μ = 100` to
converge on 1105 steps where 1105 uniform steps do not, and to beat a 20000-point uniform grid
(−47.3 ppm against −60.6). So "shooting throws away step control it already owns" is only half
true: the *consumption* side is shipped, and what is missing is the *derivation* being automatic.

⚠⚠ **AND ON THE WRAPPING FIXTURE THE DERIVED GRID DOES NOT HELP — MEASURED:**

    grid                    steps   max LTE (× tol)   at t/T
    uniform                   500      4.83e+05       0.348697
    uniform                  1428      1.67e+05       0.346181
    transient-derived        1429      2.64e+05       0.348172

**The derived grid is WORSE than the uniform grid of the same count.** So the premise — that
adaptive stepping fixes this fixture — is falsified, and B7 cannot be justified on it.

⚠ **THE LTE PEAK IS AT THE RESET, ON EVERY GRID.** With `ic = 0.31` and the integral advancing
2.0 per period the wraps are at `t/T = 0.345` and `0.845`; all three peaks sit at ≈ 0.348. The
number being reported is the **discontinuity**, and no step size makes a discontinuity's local
truncation error small — which is C5's finding ("a discrete map has no undefined instant")
arriving from the LTE side.

⚠⚠⚠ **AND THE STRUCTURAL REASON A FROZEN GRID CANNOT CARRY AN EVENT.** `Idtmod.next_event` is a
LINEAR PREDICTION FROM THE LAST ACCEPTED POINT, returning `inf` before the first step — measured
on a fresh instance: `inf` at every `t`. It is meaningful only *during* a traversal. Measured
during one, the transient lands NEAR the resets but never on them (gaps 0.3%–23% of the local
step), and its own docstring says that is by design: it "only needs to BRACKET the corner".

**So a grid derived from a past traversal cannot represent an event whose time depends on the
state — because the event time MOVES as the Newton iterates.** Freezing the grid and localising
events are in direct tension. That is the reproducibility problem the entry named, now with a
mechanism instead of a worry, and it means the two halves of B7 are **not one item**:

  * **B7a — automatic grid derivation for STIFF SMOOTH problems.** ✅ **BUILT 2026-09-04** as
    `PSS.lte_grid(period, ...)`, returning `(fracs, seed)` ready for `solve(grid=…, x0=…)`.
    The `μ = 100` gate was re-run before promoting and still holds: **1105 LTE-chosen steps
    converge at −47.3 ppm where 1105 UNIFORM steps do not converge at all**, for `trap` and
    `gear` alike. Step spread on that grid is 16438×. Pinned by
    `test_lte_grid_derives_a_frozen_nonuniform_grid_that_solve_accepts`, which asserts the
    fractions sum to one period, are genuinely non-uniform, and are not WORSE than a uniform
    grid of the same count against a fine reference.
  * **B7b — event localisation under shooting.** NOT solved by a frozen adaptive grid, and the
    measurements above say so. It needs the event time to be an unknown the Newton solves for, or
    a formulation where the reset is a state-dependent map rather than a grid feature. That is a
    different and much larger item, and it is A6's defect proper.

**Recommendation: build B7a, and re-file B7b against A6 rather than here.**

---

⚠⚠ **B7c — MONOTONE GRID REFINEMENT ACROSS SHOOTING ITERATIONS. ANDREAS'S PROPOSAL, MEASURED
2026-09-06. THE TWO-STAGE FORM WORKS; THE FROM-SCRATCH FORM COSTS 3x FOR NOTHING; AND THE CASE
THAT WOULD JUSTIFY IT OVER B7a COULD NOT BE BUILT.**

The proposal: do not freeze the grid at the first iteration. Let each shooting iteration ADD the
points it needs to hold LTE at tolerance, keeping the previous iteration's points, so the step
controller keeps working and the grid only grows. ⚠ With the constraint that a new point may not
be placed too close to an existing one -- that creates arbitrarily small steps and its own
numerical trouble. All numbers below: van der Pol at `mu = 100`, `gear`, `reltol = 1e-7`, one
fixture, one seed.

**1. The naive form -- union whole grids -- is WORSE THAN USELESS.** Merging each iterate's
independently-derived grid gives **5.21x** the points AND degrades accuracy (+619 ppm against the
solution grid's -3.8 ppm), because the merged spacing becomes wildly uneven. More points in the
wrong distribution is worse than fewer well-placed ones.

**2. ⚠ AND A WARMUP DOES NOT FIX IT — the mechanism is not what I predicted.** I expected a bad
early iterate to pin points, which a warmup removes. Measured with post-warmup iterates that are
SHRINKING perturbations (3e-2 -> 0) of the settled point, the union still bloats **3.6x - 5.2x**.
The per-iterate grids are nearly the same SIZE (1238, 1209, 1185, 1169, 1137) but their points
barely coincide: **a tiny perturbation of `x_0` shifts every step boundary slightly, so the union
is close to a SUM.** Adaptive step POSITIONS are not stable under small state perturbations, and
that is warmup-proof.

**3. The TARGETED form -- add only where the current grid is coarser than the controller asks --
preserves accuracy** (-2.3 to -4.1 ppm) but still costs **3.0x - 3.2x**. `gamma`, the tolerated
coarseness before subdividing, is a weak knob: 1.5 -> 5.0 buys only 3.24x -> 2.79x and starts
costing accuracy at 5.

⚠ **3b. THE SEPARATION RADIUS IS A STRONGER KNOB THAN `gamma`, AND ITS SETTING IS `delta = 1`.**
Scaled to the CANDIDATE's own intended local step -- "do not add a point if one already sits
within a full step of the one you wanted" -- swept on the post-warmup union:

    delta   points  growth   period err
     0.0     5925    5.21x    +619.4      (no rule at all)
     0.3     4112    3.62x     -28.2
     0.5     3627    3.19x     +20.6
     0.7     3361    2.96x      -5.9
     1.0     3115    2.74x      -9.3      <- the knee
     1.5     1708    1.50x    +339.6      <- rejects points the solution needs

**`delta = 1` nearly HALVES the bloat at no accuracy cost** (5.21x -> 2.74x, -9.3 ppm against the
solution grid's -3.8). ⚠ And the knee is sharp: at 1.5 the count falls to 1.50x -- close to the
two-stage form's 1.24x -- but accuracy collapses by 40x, because the rule has started starving the
solution's own resolution. **So the count cannot be bought by tightening the rule; the two-stage
STRUCTURE is what makes 1.24x possible.**

**4. ⚠⚠ THE TWO-STAGE FORM IS THE ONE THAT WORKS (Andreas's follow-up): solve on a fixed grid
FIRST, then refine.** Because the iterates are then already at the solution, the refinement
criterion stops firing and the grid REACHES A FIXED POINT:

    stage 0 (decimated fixed grid)   569 pts   +424.8 ppm
    after one refinement            1407 pts    +18.2 ppm
    stage 2                         1411 pts    +18.2 ppm   (+4 points)
    stage 3                         1415 pts    +18.2 ppm   (+4 points)
    (the solution's own grid        1137 pts     -3.8 ppm)

**1.24x the solution's grid, and a 23x accuracy recovery in one pass.** That is the whole
difference from the from-scratch form, and it is exactly the "start after a warmup" insight taken
to its conclusion.

**5. It cannot improve a grid that is already good.** From a B7a-quality grid: 1137 pts at
-3.81 ppm -> 1153 pts at -3.98 ppm. It adds ~16 points and drifts marginally worse. So this is a
REPAIR MECHANISM FOR AN UNDER-RESOLVED GRID, not a replacement for `lte_grid`.

⚠ **THE REMAINING GAP AND WHY IT IS STRUCTURAL.** Refined-from-decimated lands at +18.2 ppm where
the solution's own grid gets -3.8, because the scheme can only SUBDIVIDE -- never move or remove a
point -- so it inherits the starting grid's placement. Closing that needs COARSENING, which breaks
the monotonicity that buys reproducibility for `factored_period()`. A real trade, not an oversight.

⚠⚠ **THE CASE THAT WOULD JUSTIFY THIS OVER B7a COULD NOT BE CONSTRUCTED, AND THE FAILURE IS
INFORMATIVE.** B7a derives a grid once from a settled run and freezes it; the scheme above only
wins where the SOLUTION's stiff regions are not where the WARMUP's were. Attempted: a high-Q tank
(Q = 100, envelope 32 periods) with a diode clamp, warmed up for only 3 periods so the diode is
still OFF (0.449 V) while the steady state conducts each cycle (0.549 V). Measured against a
4259-point reference, the warmup-derived grid (250 pts, rel err **1.540e-3**) is barely worse than
the settled-derived one (260 pts, **1.482e-3**) -- a 4% difference. The diode CLAMPS the tank, so
the two states stay qualitatively similar and the controller places points similarly.

⚠ **And the obstacle looks structural, not a failure of imagination.** What makes a warmup-derived
grid wrong is a SLOW MODE (the warmup has not settled) -- and a slow mode is precisely what drives
`|lambda_2| -> 1`, i.e. what makes the PSS solve itself hard (section 0's organising fact). The
condition that motivates the feature and the condition that breaks the solver are the same
condition.

✅ **BUILT 2026-09-06 as `PSS.refine_grid(grid, x0, period=...)`** on exactly that basis -- a repair
path, not a replacement. Pinned by
`test_refine_grid_repairs_an_under_resolved_grid_and_reaches_a_fixed_point`, which asserts the
recovery, the fixed point, AND that the union rule's `delta = 1` is inert here (see below).
⚠ **THE SEPARATION CONSTANT IS NOT PORTABLE BETWEEN THE TWO RULES AND THIS WAS HIT FOR REAL:**
the method was first written with `delta = 1` from the union sweep and silently did nothing --
4 points a stage, +424.8 -> +423.5 ppm. In the SUBDIVISION rule the inserted points already sit
about one WANTED step apart, so a full step of clearance refuses them all. `0.25` is the measured
value for this rule.

**Recommendation: do NOT adopt this in place of `lte_grid`.** It is worth building only as an
explicit repair path -- "I have a grid I suspect is too coarse, improve it" -- where its measured
behaviour (converges, 1.24x, 23x recovery) is exactly right. ⚠ And if it is built, use the TARGETED
rule with the separation radius scaled to the CANDIDATE's own intended step AT `delta = 1` (see 3b):
scaling it to the existing grid's local gap was tried instead and admitted only 33 of 1158 points on
a coarse grid, producing grids that did not converge at all.


---

⚠⚠⚠ **A DESIGN HYPOTHESIS FROM ANDREAS, AND IT GOES TO THE CENTRAL POINT.** ⚠ **PROVENANCE
CORRECTED 2026-09-04 — AN EARLIER VERSION OF THIS PARAGRAPH PRESENTED IT AS ESTABLISHED
COMMERCIAL PRACTICE ("from practice with commercial SPICE PSS engines: they do not use a fixed
grid at all"). ANDREAS HAS SINCE SAID PLAINLY THAT IT IS A GUESS — the commercial simulator's source is not
visible, so what it does inside is not knowable from outside.** The hypothesis is that the
stepping is controlled entirely by the inner transient with **the last step placed on the period
boundary**; it is a plausible design from an experienced user, and it is **NOT** an appeal to
authority about any shipping tool.

⚠⚠ **THE HYPOTHESIS'S PROVENANCE AND ITS TECHNICAL CONTENT ARE SEPARATE, AND THE SECOND SURVIVES
INTACT** — everything below was checked against THIS codebase and holds regardless of what any
commercial tool does:

  * **The shooting Newton never needed a reproducible grid.** `_traverse` and
    `_traverse_solved_history` return the endpoint AND the sensitivity `P = dx/d(x_0, x_{-1})`
    **from one walk**. So `φ` and `M` already come from the *same* traversal and are consistent
    with each other on whatever steps that traversal took — which is exactly the commercial
    arrangement. The reproducibility requirement recorded above comes only from
    `factored_period()` re-traversing SEPARATELY for the stored factors, and a deterministic step
    controller is "a function of the state only" — the escape clause the original entry already
    allowed and I did not follow up.
  * **What the frozen grid actually buys is the PERIOD COLUMN, and only under one construction.**
    Fractions exist so that `h_i = fr_i · T` and hence `dh_i/dT = h_i/T`, which the autonomous
    period column differentiates through (see the comment at `_traverse`: "Every step scales
    together (`h = T/(N-1)`), so `dh/dT = …`"). That is the whole reason the grid must be a
    fraction list rather than a step list.

⚠⚠ **AND ANDREAS'S "LAST STEP ON THE PERIODIC BORDER" IS PRECISELY THE CONSTRUCTION THAT REMOVES
IT.** Let the transient step freely and truncate the final step to land on `T`. Then
`dh_i/dT = 0` for every interior step and `dh_N/dT = 1`, so the period column comes from the LAST
STEP ALONE and tends to `ẋ(T) = f(x(T))` — Aprille & Trick's exact continuous statement, with no
assumption about how steps scale. ⚠ It is also arguably the *better* derivative: spreading `δT`
proportionally over every step perturbs the trajectory in a way the physical period change does
not.

**So B7 should be rebuilt on the closing-step formulation rather than on frozen fractions**, and the
"direct tension between freezing and event localisation" recorded above dissolves — because
nothing needs freezing. What survives from the gate above is still true and still useful:

  * the derived-grid experiment's LTE numbers (a **frozen** grid derived from a past traversal
    does not help, and is worse than uniform at equal count) — that measures the frozen approach,
    which is the one being replaced;
  * the LTE peak sitting at the reset on every grid, so the number reported there is the
    discontinuity;
  * `next_event` being a prediction from the last accepted point, meaningful only DURING a
    traversal — which under the commercial scheme is exactly where it now WOULD be consulted,
    so this stops being an obstruction and becomes the mechanism.

⚠⚠ **STANDING VERDICT AFTER THE RETRACTION AND THE PROVENANCE CORRECTION (2026-09-04):** the
closing-step design is a **hypothesis with no demonstrated advantage**. Its claimed advantage
(gates 1 and 4) was my own finite-difference artefact and is withdrawn; its appeal to commercial
practice is a guess and is withdrawn too. What remains is a built, tested, opt-in implementation
and gate 2's determinism result. **That is a reasonable place to stop until there is a reference
accurate enough to rank the two period columns** — and building further on either the retracted
measurement or the guessed provenance would be building on nothing.

**ORIGINAL PLAN — B7c, superseding B7a and B7b:** drive the period with the inner transient's own
step control; force the closing step onto `T`; take the period column from that step alone. Gates,
in order: (1) the autonomous period column against the proportional-scaling one on a smooth
fixture, where both are valid and must agree; (2) `M` from `factored_period()` against `M` from
the solving traversal, which the deterministic controller should make reproducible; (3) the
wrapping fixture's LTE, which is the case the frozen grid could not fix; (4) van der Pol at
`μ = 100`, which must not regress.

---

⚠⚠⚠ **RETRACTED 2026-09-04, SAME DAY, BY BUILDING IT. GATES 1 AND 4 BELOW ARE WRONG — THE ERROR
WAS IN MY FINITE-DIFFERENCE INSTRUMENT, NOT IN THE SHIPPED COLUMN.** Read them for the method and
not for the conclusion.

Implementing the closing-step column made the analytic derivative available on BOTH conventions,
and comparing all four against the same reference says:

    npts   ana_prop     ana_close    fd_prop      fd_close     (relative to ẋ)
     120   1.359e-05    2.270e-07    8.395e-03    8.794e-08
     240   3.430e-06    1.357e-08    4.183e-03    3.004e-07
     480   8.579e-07    7.813e-10    2.087e-03    3.132e-07

**`ana_prop` — the SHIPPED column — converges at ~`O(h²)` and is small.** The `O(h)` sequence
gates 1 and 4 reported is `fd_prop`, MY finite-difference construction, which is `O(h)` wrong on
its own. Rebuilding the grid at `T + δT` also moves the MANUFACTURED OPENING STEP, and that
artefact is what was being measured. Gate 4's "4.2% relative on the stiff grid" comes from the
same instrument and falls with it.

⚠⚠ **§D SHAPE 0j — I COMPARED TWO INSTRUMENTS AND BLAMED THE SUBJECT.** Both columns were taken
by finite difference "off the SAME `_traverse`, so the comparison is between the two STEP
CONVENTIONS and nothing else" — which was true and beside the point, because one of the two FDs
was wrong. **The analytic column the code already computes was available the whole time via
`want_dT=True`, and checking either FD against it would have caught this in one line.** When a
measurement says a shipped implementation is wrong, validate the instrument against that
implementation before believing it.

⚠ **WHAT SURVIVES.** The closing-step column is BUILT and verified against its own finite
difference to `5.6e-07`, flat across grids — that is FD truncation noise, so the implementation is
correct. It is behind `_period_column = 'closing'`, **default unchanged**. But its JUSTIFICATION
is gone: `ana_close` looks better than `ana_prop` above (7.8e-10 against 8.6e-7 at 480 points),
and **that ranking is NOT established** — the reference is a centred difference of the discrete
waveform whose own accuracy is nowhere near `1e-9`, so agreement at that level means the two share
a construction rather than that one is nearer the truth. **A reference accurate enough to rank
them is the missing piece, and B7c should not be built further until there is one.**

---

✅⚠ **GATE 1 RUN 2026-09-04 — IT PASSES, AND MORE STRONGLY THAN IT ASKED.** ⚠ **(RETRACTED — see
above.)** The gate was "both
constructions are valid on a smooth orbit, so they must agree". They do not merely agree: **they
differ by `O(h)`, and the shipped one is the one that is wrong.**

Both period columns taken by finite difference off the SAME `_traverse`, so the comparison is
between the two STEP CONVENTIONS and nothing else. `ẋ` is a centred difference of the converged
waveform about `t = 0` (at a periodic solution `ẋ(T) = ẋ(0)`):

    npts   |prop − ẋ|      |last − ẋ|      |prop − last|
     120   1.500319e-02    1.571766e-07    1.500321e-02
     240   7.480561e-03    5.371494e-07    7.480558e-03
     480   3.733940e-03    5.602383e-07    3.733939e-03
     960   1.865252e-03    5.617951e-07    1.865252e-03

⚠⚠ **THE PROPORTIONAL COLUMN IS FIRST ORDER.** `|prop − ẋ|` halves cleanly per doubling — 1.50e-2,
7.48e-3, 3.73e-3, 1.87e-3 — while `|last − ẋ|` sits at ~5.6e-07 and stops improving, which is the
REFERENCE's floor and not the column's. And `|prop − last|` equals `|prop − ẋ|` to six digits, so
the discrepancy is entirely in `prop`.

**So the shipped autonomous period column carries an `O(h)` error that the closing-step
construction does not.** That is a stronger reason to build B7c than the adaptive stepping was:
it is a correctness improvement to the autonomous Newton's Jacobian, independent of whether the
grid ever adapts.

⚠ **WHAT THE MEASUREMENT DOES AND DOES NOT ESTABLISH.** The reference is a centred difference of
the DISCRETE solution, so what is pinned is that the two conventions disagree at `O(h)` and that
the closing-step one agrees with the solution's own derivative. That is the right target: the
Newton solves the DISCRETE equations, so its period column should be the discrete map's
derivative, which is exactly what the closing-step convention computes exactly and the
proportional one approximates. It is NOT a comparison against continuous truth, and it should not
be quoted as one.

⚠ **AND MY FIRST TWO ATTEMPTS AT THIS MEASUREMENT WERE WRONG**, both in the fixture rather than
the finding. `ppv()` is gear-only so the first run died on `method='trap'`; then I modified
`hs[-1]`, which is not the closing step — the walk is
`for _j, t in enumerate(times[1:]): dt = hs[_j]`, so the closing step is `hs[len(times) − 2]` and
it reaches `times[-1]`. With the wrong entry modified the last-step column came back ≈ 0 and
briefly looked like a falsification of the whole construction.

---

✅ **GATE 2 RUN 2026-09-04 — PASSES ON BOTH HALVES, and it removes the last stated obstruction.**

    (a) factored_period() replay vs the solving traversal, on the STORED grid
        trap  (plain)           max|ΔM| = 3.442e-15   relative 3.4e-15
        gear  (solved_history)  max|ΔM| = 8.091e-15   relative 5.4e-15

    (b) adaptive grid determinism, re-derived from the SAME state
        2292 steps both runs    max|t1 − t2| = 0.000e+00   BIT-IDENTICAL
        (4908 steps at reltol 1e-10, so it is a function of tolerance too)

**(a)** the replay is faithful to round-off — not bit-identical, because the summation order
differs, but at machine precision. **(b)** the adaptive controller is **bit-deterministic** given
the same state and tolerance. So "a deterministic step controller is a function of the state only"
is no longer an assumption: at the converged `x_0`, `factored_period()` would regenerate exactly
the grid the solve used.

⚠⚠ **BUT THERE IS A CONSEQUENCE GATE 3 MUST MEASURE, AND IT IS NOT A DEFECT — IT IS WHAT AN
ADAPTIVE PSS MEANS.** `factored_period()` re-traverses at the CONVERGED `x_0`, while the Newton's
last iterate traversed from `x_0^(k)`. With a fixed grid those give the same map. With an adaptive
one they need not: **the converged solution satisfies periodicity on the grid THAT traversal
chose**, so a re-traversal from a slightly different state can pick a slightly different grid and
the residual reappears at the size of the local error. That is the accepted behaviour of every
commercial adaptive PSS — the solution is accurate to the transient tolerance and re-integration
moves it by that much — but it must be MEASURED here rather than assumed bounded, because this
codebase's convergence test does not currently expect it. **Gate 3 should report the residual at
`x_0*` under a re-derived grid, alongside the wrapping fixture's LTE.**

---

⚠✅ **GATE 3 RUN 2026-09-04 — BOTH PARTS ANSWERED, AND THEY POINT OPPOSITE WAYS.**

**(A) THE PHASE HYPOTHESIS WAS RIGHT AND IT DOES NOT RESCUE THE WRAPPING FIXTURE.** The failed B7
gate derived its grid from a free-running transient window whose phase was arbitrary — finest
steps at relative 0.0 and 0.5 while the resets are at 0.345 and 0.845. Deriving at the CONVERGED
`x_0` fixes that by construction: the six finest steps land at `t/T = 0.3443`, on the first reset,
at a step ratio of **417199×**. And the LTE barely moves:

    uniform 1428      max_lte = 4.1377e+04   at t/T = 0.346181
    derived at x_0*   max_lte = 3.8970e+04   at t/T = 0.845054   (1743 steps)

**6%, and THE PEAK RELOCATED TO THE OTHER RESET.** Refine one reset to a step ratio of four
hundred thousand and the peak simply appears at the second one, at the same magnitude. ⚠ **That is
the signature of a discontinuity, not of under-resolution**, and it closes the question the
earlier gate left ambiguous: the wrapping fixture's LTE is irreducible by stepping, however the
steps are chosen. B7c must NOT be justified on it, and A6's item should say so.

**(B) RE-DERIVING THE GRID REOPENS THE RESIDUAL, AND THAT BOUNDS WHAT "CONVERGED" CAN MEAN.**

    own grid    ( 480 steps)   ||x_end − x(0)|| = 2.6990e-15   rel 1.5e-15
    re-derived  (4909 steps)   ||x_end − x(0)|| = 1.6009e-04   rel 8.9e-05
    solve reltol = 1e-12

Machine zero on the grid the solve used, `8.9e-05` on a finer one. ⚠ **The right reading is that
this measures the SOLVE'S DISCRETISATION ERROR, not a defect** — the 480-step trapezoidal answer
is accurate to ~1e-4, and the finer grid exposes it. But the consequence for B7c is concrete:
**with an adaptive grid the achievable residual is bounded by the difference between grids, not by
`reltol`**, so the Newton must stop at that floor rather than chase 1e-12. That is exactly how
commercial adaptive PSS behaves, and it is now measured here rather than assumed.

⚠ **A number NOT to quote from this run:** `||phi_own − phi_rederived|| = 2.3304e-02` looks like
the grid dependence and is mostly the OPENING STEP — the two grids manufacture `x(0)` with
different first steps. The residual above is the meaningful figure.

⚠ **AND THE FIXTURE TRAP THAT PRODUCED A FALSE 2.3e-02 RESIDUAL FIRST TIME.**
`_period_state[1]` is the MANUFACTURED `x(0)`, not the solved unknown `x_in`. Feeding it back to
`_traverse` as an input manufactures a SECOND opening step, so the "residual" came out at
`h·|ẋ| = (2π/480)·1.79 = 0.0234` — which is what was measured, to three digits. A converged solve
showing a 1e-2 residual should have been read as a fixture error immediately; the arithmetic
identifying it took one line.

**GATE 3 VERDICT: B7c's justification is gate 1 (the `O(h)` period column) and gate 2 (determinism),
NOT the wrapping fixture.**

---

✅✅ **GATE 4 RUN 2026-09-04 — PASSES, AND GATE 1's RESULT IS AMPLIFIED ON THE STIFF CASE.**
⚠⚠ **(RETRACTED — the 4.2% is my finite-difference instrument, not the convention. See the
retraction at the head of gate 1.)** The
regression risk was the opposite regime from gate 1: van der Pol at `μ = 100` on its LTE-chosen
grid, **step ratio 16438×**, where the proportional convention scales every UNEQUAL step while the
closing-step one moves only the last. If the closing-step column were ever going to be the worse
choice, it would be here.

    uniform,  1105 steps                    NoConvergenceError   ← benchmark reproduced
    LTE grid, 1105 steps   T = 162.830391669
        |prop − ẋ| = 8.0500e-02    relative 4.238e-02
        |last − ẋ| = 1.7595e-03    relative 9.263e-04

**The two conventions differ by 4.2% RELATIVE on the stiff grid, and the closing-step column is
46× closer to the solution's own derivative.** On the smooth uniform grid of gate 1 the gap was
`O(h)` and shrinking; here it is a percentage-level error in a Jacobian column. That is the
expected direction — at a 16438× ratio, scaling every step is a very different perturbation from
extending the last — and it means the proportional convention is worst exactly where the
non-uniform grids B7c exists to enable put it.

⚠ **WHAT RESTS ON WHAT, because the reference is weaker here than in gate 1.** `ẋ` is a centred
difference about `t = 0`, and on a grid with a 16438× ratio the spacing either side of `t = 0` is
NOT symmetric, so the reference's own accuracy is not established. What gate 4 establishes on its
own is that **the two conventions differ by 4.2%**; WHICH of them is right rests on gate 1, where
the grid was uniform, the reference sound, and `prop` converged to `last` at a clean `O(h)`. Quote
the 46× as indicative, and the 4.2% divergence as the measured fact.

⚠ **And the regression baseline is intact:** uniform at 1105 steps still fails to converge while
the LTE grid at 1105 succeeds, which is the benchmark's headline and the thing gate 4 was there to
protect.

**ALL FOUR GATES ARE NOW RUN.** B7c is justified by gates 1 and 4 (the period column is `O(h)` on
smooth grids and 4.2% wrong on stiff ones) and made safe by gate 2 (bit-deterministic
regeneration). Gate 3 removed the wrapping fixture from its justification and bounded what
"converged" can mean under an adaptive grid. **Build order: the closing-step period column first —
it is a self-contained correctness fix that needs none of the adaptive machinery.**

### B8. All integration methods in PAC, pnoise and the adjoint paths — ⚠ **BUILT 2026-09-04**

✅✅ **THE PLAIN TRANSPOSED REPLAY SHIPPED** as `_monodromy_matvec_transposed_plain`, so
`matvec_transposed` now dispatches on kind instead of refusing. The reverse recursion is derived
separately for the two one-step companions: euler (`b = 0`) has one term, trapezoidal (`b = -1`)
keeps the pair `(P, Pq)` — a *different* pair from Gear-2's `(P_n, P_{n-1})`, which is why the
solved-history replay could not be reused — and collapses to **one transposed solve per step**
through a shared bracket, the same cost as gear.

Gated against the FORWARD replay transposed densely, with **gear as a harness control**:

    driven  euler/120  1.16e-13   trap/120  3.90e-13   gear/120  1.46e-13
    driven  euler/240  5.40e-13   trap/240  1.45e-12   gear/240  7.92e-13
    vdp                trap/120  2.89e-15   gear/120  1.04e-15

⚠ **TWO TEST BUGS FOUND ON THE WAY IN, BOTH OF WHICH PASSED FIRST.** See §D shape 0m. The
driven fixture was `Q ~ 1e-3`, so the monodromy decayed to numerical zero and the check compared
zero against zero, reporting **1.9e-46 and printing OK**. And the autonomous assertion demanded
`|λ| ≈ 1` to 5e-3 and failed at 0.992007 — *the test* was wrong, not the code: `M` itself carries
that deficit at 120 points on a `Q = 8` cycle, and measured rather than loosened it converges at
`O(h²)` (7.99e-3 → 1.87e-3, ratio 4.3, gear exact).

The original entry follows.


⚠ **The shooting SOLVE already supports every integrator that exists.** `integrator.py` defines
exactly three — `EulerIntegrator`, `TrapezoidalIntegrator`, `Gear2Integrator` — and
`PSS.solve` accepts `'euler'`, `'trap'`/`'trapezoidal'`, `'gear'`/`'gear2'`. So the request is
already satisfied there, and the gap is one level up.

**The gap is the ADJOINT surface, which is Gear-2 only.** `factored_period().matvec_transposed`,
`ppv`, `_forced_replay_transposed`, `covariance`, `oscillator_covariance` and `_lyapunov_pieces`
all require the solved-history factors and refuse otherwise. So **PAC, pnoise, the PPV and every
noise result are `method='gear'` only** — `euler` and `trap` reach the PSS and stop there.

⚠ **What it costs to lift:** the transposed replay is *derived* for a two-step companion with
`b = 0` — the docstring states the pair map and its transpose explicitly, and raises for
`len(alphas) < 3` and for `b != 0`. A one-step method needs **its own reverse recursion**, not a
special case of this one. That is the whole of the work, and it is the same recursion an IRK
adoption would need (see the parallel-shooting entry), so the two should be costed together.

⚠ **And there is a reason to want it beyond completeness:** `trap` is one-step, so a `trap`
adjoint would pay none of the parallel-shooting low-order-restart bias, and `x0_unknown=True`
exists precisely for the one-step path. Today that path cannot produce a PPV.

### B12. Gourary's regularisation against our bordered solve — ⚠ **MEASURED 2026-09-04; IT IS THE SAME METHOD, AND ADOPTING IT WOULD BE A REGRESSION**

The one item in the Gourary corpus that could have changed a design decision, so it was measured
rather than filed. *M. M. Gourary et al., "A numerical technique for time domain noise analysis of
oscillators", ECCTD 2007, 1002–1005.*

**The method.** Both schemes exploit the *same* identity, `vᵀ(I − αM) = (1 − α)vᵀ`, where `v` is
the PPV and `α = exp(−2πjfT)` vanishes to 1 at every harmonic. Gourary **substitutes the constant
row `vᵀ`** for one equation and moves the vanishing factor to that RHS entry. We **border**, and
carry the pole analytically as `y = w + s·u/(1 − α)`.

**Near-carrier sweep (`f₀ + Δf`, van der Pol at `Q = 15.92`, 400 points):**

    Δf/f₀     cond(plain)   cond(bordered)  cond(gourary)   |y_bord − y_gour|/|y|
    1e-01      8.59e+00       1.57e+01        9.35e+00          3.22e-13
    1e-03      7.95e+02       8.25e+01        7.57e+01          4.16e-13
    1e-06      7.95e+05       8.29e+01        7.61e+01          4.15e-13
    1e-09      7.95e+08       8.29e+01        7.61e+01          4.15e-13
    1e-12      7.82e+11       8.29e+01        7.61e+01          4.15e-13

⚠⚠ **THE TWO ANSWERS AGREE TO 4.15e-13 AT EVERY OFFSET, FLAT.** Gourary is not a better
regularisation; it is the *same* regularisation with the pole parked somewhere else. Its 9 %
edge on `cond` (76.1 against 82.9) is noise.

⚠⚠ **AND IT FAILS AT AN EXACT HARMONIC EXACTLY AS WE DO.** Its matrix stays at cond 76.1 at
`α = 1` — which is what makes the published claim *look* stronger — but its RHS entry is
`(vᵀb)/(1 − α)`, a division by zero. The flat condition number describes the matrix, not the
solve. Ours divides by zero in `s·u/(1 − α)`. Same pole, same place, different clothing.

**Two ways it is WORSE, neither visible in a condition number:**

⚠ **1. IT HAS A FREE PARAMETER AND WE DO NOT.** The method replaces *one* equation, and the paper
fixes which by fiat ("the equation corresponding to the output node"). Measured over every
possible row, against the bordered answer:

    n = 4    best 4.15e-13   worst 7.46e-10
    n = 10   best 1.19e-12   worst 5.04e-02
    n = 16   best 6.85e-13   worst 5.48e-02

**Ten orders of magnitude on a choice the paper makes by convention.** `argmax|v|` happens to land
near the best row, but that is a heuristic nobody derived and it is not the paper's rule.

⚠⚠ **2. IT IS SENSITIVE TO `v`; WE ARE IMMUNE TO IT.** This is the finding worth keeping, and it
is a property of *our* code that had not been stated. Measured **through the shipped
`_deflated_solve`**, not through a reimplementation of it:

    perturbation   border ROW (v)    border COLUMN (u)    gourary (v)
      1e-08          9.2e-15             9.0e-09           3.2e-08
      1e-04          2.7e-14             1.0e-04           4.2e-04
      1e-02          3.1e-15             7.8e-03           3.3e-02

For `α ≠ 1` the system `(I − αM)y = b` is **nonsingular**, so `y` is already determined by `b`;
the border row only picks a well-conditioned route to it, and *any* `v` not orthogonal to the null
direction gives the same answer. Gourary **replaces an equation** with `v`, so an error in `v`
corrupts the system itself and passes straight through.

✅ **THE OPERATIONAL CONSEQUENCE, WHICH IS THE REAL DELIVERABLE:** PAC's accuracy is capped by the
**orbit tangent `u`**, which enters the reconstruction, and **not by the PPV `v`**. Tightening
`ppv()`'s tolerance to improve a PAC result optimises the wrong vector. Pinned by
`test_the_deflated_solve_is_capped_by_the_TANGENT_not_by_the_PPV`.

⚠ **PRIORS, STATED BEFORE THE RUN AND SCORED AFTER.** Mine: "Gourary wins on conditioning and that
does not make it better, because they exploit the same identity" — **held**. The docs session's:
"Gourary holds accuracy closer to the carrier, both hit a null-vector floor at the same Δf" —
**falsified in both halves**: they are identical near the carrier, and only *one* vector caps
either method, `u` for ours and `v` for Gourary. The shared-floor hypothesis was the reasonable
guess and the measurement did not support it.

**VERDICT: do not adopt.** The bordered solve stays. This is the first item from that corpus
measured against what we already have, and what it establishes is that the published method is
our method with an extra free parameter and a vector sensitivity we do not have.

### B13. Krylov cost at high Q above m=12 — ✅ **MEASURED 2026-09-04; THE WORRY IS FALSIFIED AND THE REAL DRIVER IS SLOW NODES**

The binding gap for the matrix-free question: every existing benchmark varies `m` at **one** `Q`,
and the concern was the *pairing* — at high `Q` the multipliers crowd the unit circle, so `I − M`
crowds zero, and GMRES was expected to need `O(m)` iterations exactly where large `m` makes
matrix-free worth having.

**Fixture.** A van der Pol at a target `Q` plus an RC ladder whose first `nslow` sections have time
constants **straddling the period**. ⚠ The straddle is the whole fixture: a ladder that decays
inside one step raises `m` without adding modes near the unit circle, which is what makes `nslow`
and `m` separable. Bordered `(I − M)w = b`, GMRES at `rtol = 1e-10`.

**Fixed `m = 32`, sweeping only `nslow`:**

    nslow      0    4    8   14   22   30
    Q =   8    4    8   11   16   23   29
    Q = 256    5    8   11   17   23   30

⚠⚠ **A 32× CHANGE IN `Q` MOVES THE ITERATION COUNT BY AT MOST ONE.** Iterations track `nslow`
(roughly `1 + nslow`) and ignore `Q` entirely, and ignore `m` — which is held fixed here.

**Why:** Krylov iteration count is set by the number of **distinct eigenvalue clusters**, a
spectral-spread property. `Q` controls *conditioning*, which is a different thing, and the
oscillator contributes **one** mode near 1 however high `Q` goes. Slow nodes each contribute
another.

✅ **THE OPERATIONAL RULE INVERTS.** Matrix-free is **safe on a high-Q oscillator** — the case
everyone worried about — and degrades on a circuit with **many slow nodes**, at any `Q`. A
designer's high-Q tank costs nothing; a bias network with a dozen long time constants costs
linearly. Pinned by `test_krylov_cost_ignores_Q_and_tracks_the_SLOW_NODE_COUNT`.

⚠ **`|λ| > 0.9` IS A POOR PROXY** and was nearly reported as the driver: it counted 1/2/2/3/4/5
across that sweep while iterations went 4/8/11/16/23/29. A threshold count of near-unit
multipliers is not the number of distinct clusters, and only the latter predicts.

⚠ **THE FIRST FIXTURE MEASURED THE WRONG OBJECT AND IS KEPT AS THE CONTROL.** Its ladder ran `τ`
from 1e-2 down to 1e-9 against a period of `2π`, so every added mode decayed inside one step:
**4 iterations at every `(Q, m)` across `Q` = 8…1024 and `m` = 8…32**. Uniform, with nothing
predicting uniformity. That run is the `nslow = 0` column — at `m = 32` iterations go **4 → 33**
purely by making the ladder slow, which is the cleanest separation in the whole measurement.

⚠ **BOUNDARY: `Q = 1024` WITH A FULLY SLOW LADDER DOES NOT CONVERGE** at 400 points. That is the
shooting solve failing, not the Krylov solve, and it is a different limit from the one this entry
closes.

⚠ **A CAVEAT ON SCOPE, since it bounds what this licenses:** iterations reach 29–33 at `m = 32`
(width 64), i.e. approaching `n/2`, where a direct solve wins regardless. The result says `Q` is
free, not that the matrix-free route is unconditionally cheap.

⚠⚠ **AND A SECOND SCOPE LIMIT, FROM A SOURCE READ AFTER THIS WAS MEASURED:** this entry measured
GMRES for the bordered **SOLVE** `(I − M)w = b`. It says **nothing** about Krylov **eigenvalue
extraction**, which is a different problem with the opposite property. Garcia, Romero & Acha 2022
(the Ritz-route source, read firsthand by the docs session): Arnoldi resolves the largest-magnitude
eigenvalues first, the Ritz route works on `I − M` with `λ = 1 − θ`, so the physical `λ₂ → 1` is
`θ₂ → 0` — **the smallest, resolved LAST** — with separation `1/θ₂ ≈ Q_λ` (3.7 / 16.4 / 64.5 /
128.5 at `Q_λ` = 3.18 / 15.9 / 64 / 128). The difficulty of the *eigen* problem scales with the
quantity being measured, exactly where the *solve* was shown free. ✅ This also explains why
`ppv()`'s Ritz selection — "the second-smallest Ritz value of `I − M`" — is the paper's own
construction and not the patch it was filed as. ⚠ Two independent structural arguments (this and
Mei & Roychowdhury's imperfect-cancellation one) now point the same way on Krylov *extraction* at
high Q; neither is a measurement. The paper validates on power networks, not RF, so it reports no
evidence about the high-Q regime either way.

### B14. Krylov recycling across the PAC sweep — ✅ **ALREADY BUILT; MEASURED 2026-09-04**

Raised as an open question ("does multi-RHS recycling cut the PAC-sweep matvec counts?"). It is
not unbuilt: `_solve_subspace` has shared one Krylov basis across the sweep since PAC landed, and
`recycle=True` is the default. What was never measured is **what it buys**.

The identity it rests on is Telichevesky's Theorem 1: `A(α) = I − αM`, so
`span{r, Ar, A²r, …} = span{r, Mr, M²r, …}` for **every** `α`. The basis is frequency-independent;
each frequency then costs a small dense least-squares over it.

**Measured against `recycle=False`, driven RLC with an RC ladder:**

    m    K    mv recycle   mv each   ratio    t recyc   t each    max rel Δy
    4    4         5          12      2.4×     0.149     0.106      3.7e-13
    4   16         5          48      9.6×     0.275     0.423      3.7e-13
    4   64         5         192     38.4×     0.787     1.702      5.1e-13
   18    4         8          24      3.0×     0.270     0.179      3.7e-14
   18   64        10         384     38.4×     0.974     2.862      3.3e-13

⚠⚠ **MATVECS ARE FLAT IN `K`; THE WALL CLOCK IS NOT — 38.4× against 2.9×.** That gap is the real
finding. The Krylov solve has **already been removed** from the sweep's cost, and what remains is
the **one forced replay per frequency** outside it, which recycling cannot touch. Further work on
this sweep must target the replays, not the linear solve. Filing "add Krylov recycling" as an
optimisation would have been work with a 1.0× ceiling.

⚠ **CITATION ORDER, corrected 2026-09-04 from the docs session's audit:** the recycling was relayed
this afternoon as an unexplored route sourced to Gourary's multi-RHS paper. It is older and
closer than that — **Telichevesky, Kundert & White, DAC 1996** already carries it, and this
record's own reading log has the mechanism: *"β = α(f_new)/α(f_old), γ = 1 − β converts a matvec
between frequencies — a nearly exact fit for a shooting sweep."* DAC'96 is the primary; Gourary
DATE 2003 is the paper that *names* the structural property (`A′ = I`) the recycling exploits.

⚠ **AND IT IS A LOSS ON SHORT SWEEPS:** at `K = 4` recycling costs *more* wall clock than solving
each (0.149 against 0.106) despite fewer matvecs — the dense least-squares over the shared basis
dominates. The win needs roughly `K ≥ 8`. Not an argument against the default, but the reason the
ratio must be read on matvecs against sweep length rather than on one timing.

Pinned by `test_pac_sweep_recycling_makes_matvecs_INDEPENDENT_of_sweep_length`, which asserts the
*contrast* (unrecycled scales with `K`, recycled does not) rather than absolute counts.

### B15. Is the orbit tangent losing accuracy in the bordered solve? — ✅ **MEASURED 2026-09-04: NO. THE B12 FOLLOW-UP IS CLOSED WITH NOTHING TO GAIN**

B12 measured that PAC's accuracy is capped **linearly by the orbit tangent `u`** and is immune to
the PPV `v`. Gourary (ECCTD 2007 §III.A) states the null vector needs no special algorithm because
`u = dx/dt|_{t=T}`. That suggested an obvious improvement: get `u` from the DAE directly instead
of from the bordered GMRES solve in `ppv()`, and lift the cap. **It was recorded as the
highest-value follow-up available.**

**REFERENCE THAT SHARES NO MACHINERY WITH THE SOLVE:** spectral (FFT) differentiation of the
converged periodic waveform. Exponentially accurate for a smooth periodic function on a uniform
grid, and it touches neither the monodromy, nor the border, nor GMRES.

    npts     |xdot| shipped      rel gap vs FFT      ratio
     200      1.999199347         3.360e-04            —
     400      1.999904917         8.313e-05          4.04
     800      1.999991761         2.067e-05          4.02
    1600      2.000002554         5.154e-06          4.01

⚠⚠ **EXACTLY `O(h²)`, WITH NO PLATEAU. The bordered Krylov solve is NOT a floor** — the tangent is
already as accurate as the trajectory it is taken from permits. (`|xdot| → 2.000003`, van der Pol's
own amplitude, is the free sanity check.)

✅ **SO THE PROPOSED IMPROVEMENT BUYS NOTHING, and that is the deliverable.** `dx/dt` taken
directly would differentiate **the same discrete trajectory** and inherit the same `O(h²)`. The
B12 cap is real, but it is the *trajectory's* accuracy showing through `u`, not the solve's.
**To improve PAC's accuracy, refine the grid or raise the order — not the tangent extraction.**

⚠ **THE PREDICTION WAS STATED FIRST AND HELD:** *"if the gap plateaus, the bordered solve is the
floor and Gourary's direct route would lift the cap; if it falls with refinement, there is nothing
to gain."* It fell at the textbook rate. A plateau was the falsifiable outcome and it did not
occur.

#### B15-obreshkov. If Obreshkov/Gourary is ever built — ⚠ DE-RISKED BY THE REVIEW SESSION 2026-09-04, recorded so it is not re-derived

From the peer session, verified by exact rational arithmetic and an order check on all five methods:
Gourary's (10) coefficients in closed form are `a_i = (−1)^i (l+m−i)!/(l+m)! · m!/(i!(m−i)!)`,
`b_i = (l+m−i)!/(l+m)! · l!/(i!(l−i)!)`; `l=0,m=1` gives backward Euler and `l=m=1` trapezoid.
**Table 1's fifth row is a typo** (`l=0, m=2` prints `a₁ = 1`; Taylor forces `−1`, and the printed
value does not converge at any `h`). **Every `a₁` is negative**, so the Jacobian block
`C − a₁hG − a₂h²G′v̇` is a positive addition (`C + hG/2` at order 4) — anyone carrying `C + hG/2`
from the trapezoidal companion will "correct" the source's minus and be wrong twice. (17), (18),
(19), (22) all carry the minus form; an earlier claim that the paper was sign-inconsistent between
(17) and (22) is retracted by its author. The sensitivity right-hand side `p₁` carries BOTH the
incoming state and the incoming DERIVATIVE sensitivity (`−b₂hG(vₙ)·h dv̇ₙ/dv₀`); dropping the second
corrupts the monodromy while Newton still converges. And (22) drops the `G′`, `C′` terms of (20)/(21)
deliberately — an inexact Newton, fine, and an APPROXIMATE monodromy, not fine when `λ₂` is the
quantity under study.

### B16. The manufactured opener DOES cap λ₂'s order on an INDEX-1 circuit — ⚠⚠ **CONFIRMED 2026-09-04, AFTER TWO WRONG REFUTATIONS OF MY OWN**

⚠⚠⚠ **THE SHIPPED `x0_unknown` DEFAULT IS TOO NARROW FOR MONODROMY ACCURACY.** It keys off
topology (on when the criterion proves index 2). This circuit is **index 1** — a plain series RLC,
no C-V loop, no L-I cutset — and on the manufactured opener its `λ₂` does not converge at the
method's order.

**Series loop `vs → R → L → C → gnd`, `Q = 100`, `trap`, reference `|(2+sh)/(2−sh)|^N`:**

    N       discrete ref     default      rel err    ratio    x0_unknown   rel err    ratio
    199     0.969080011709  0.969035749  4.568e-05    —       0.968597471  4.979e-04    —
    399     0.969074313510  0.969082426  8.371e-06   5.46     0.968954201  1.239e-04   4.02
    799     0.969072896947  0.969084458  1.193e-05   0.70     0.969042938  3.091e-05   4.01
    1599    0.969072543820  0.969080193  7.894e-06   1.51     0.969065063  7.720e-06   4.00
    3199    0.969072455665  0.969076747  4.428e-06   1.78     0.969070587  1.929e-06   4.00

✅ **`x0_unknown` GIVES TEXTBOOK `h²` — FOUR CONSECUTIVE RATIOS OF 4.00. The default is erratic
and below first order.** Independently obtained by two sessions, matching to every digit.

✅ **AND THE CROSSOVER IS REAL, so this is a TRADE and not a winner.** The default is 10× better
at `N = 199`, level at `N = 1599`, and beaten at `N = 3199`. Coarse grids favour the manufactured
opener; fine grids favour `x0_unknown`, which is **the only one that converges at the method's
order**. This reconciles with B1's opposing `Q = 20` resonator measurement rather than
contradicting it — they sit on opposite sides of the crossover.

⚠⚠ **A REFRAME FROM THE READING RECORD, WITH ITS TENSION KEPT (docs session, 2026-09-04):**
Aprille & Trick's driven-case paper, as this record filed it months ago, has *"the unknown is
`x(0)` itself; **no manufactured opening step**. The identity seed is canonical *because* `x₀` is
the unknown."* So `x0_unknown` is not a repair for an index-2 defect — it is the **return to the
source formulation**, and the manufactured opener is the **departure**; B16 measures what the
departure costs everywhere, not only on the topologies the criterion detects. ⚠ **But the tension
is real and this record already carries the other half:** the manufacturing step *"is not
scaffolding to be removed; it is what makes trapezoidal's shooting problem well-posed"* — three
reformulations without it excited the `(−1)ⁿ` companion mode. A&T's formulation had no multistep
companion current to seed. Both statements stand; their conflict is the opener problem itself,
and neither the `x0_unknown` default (30 % amplitude cost) nor the pre-roll (no convergence)
resolved it today.

⚠ **THE CONSEQUENCE FOR B1:** the topology-keyed default is defensible for *waveform accuracy on
coarse grids* and **wrong for `λ₂` accuracy on fine ones**. If a caller wants monodromy accuracy,
the **crossover** is the number that should drive the choice, not the index.

⚠⚠⚠ **I REFUTED THIS TWICE AND WAS WRONG BOTH TIMES. BOTH WERE INSTRUMENT ERRORS, AND NEITHER WAS
ARITHMETIC.**

  1. **Wrong topology.** `R('a','b')`, `L('b',gnd)`, `C('c',gnd)` put L and C **both to ground, in
     parallel**, and I applied the *series* `Q = ω₀L/R`. Signature: a flat `3.09e-02` with ratio
     **exactly 1.00** across five doublings, identical for both flags. Caught by the reflex — a
     result too clean to be real.
  2. ⚠ **The `npts`-versus-steps off-by-one, FOR THE THIRD TIME TODAY.** I set `h = T/npts` and
     raised it to the power `N = npts − 1`, so the reference covered `T·N/(N+1)` — **one step
     short of a period**. That injected an `O(h)` error into the *reference*, larger than the
     effect under test: it flattened both columns, hid the h² convergence, and manufactured the
     "default is 2× better at every N" and "`x0_unknown` at N equals default at N/2" relationships
     I reported. **Both were artefacts of the reference, not properties of the solver.**

✅ **THE FIX THAT MAKES IT UNREPEATABLE: take `N` and `h` FROM THE RETURNED WAVEFORM, never from
the requested `npts`** — `t = pss.waveform[0]; N = len(t) − 1; h = per/N` — so the exponent and
the step are consistent by construction. The script now also asserts `h·N == T` and that `h`
matches the waveform's actual first step before using either.

#### B16-review. Two corrections from the review session and a measured order gap at INDEX 1 — ⚠⚠ 2026-09-04, DECISION FOR THE USER

* **März's norm is `C¹_N`, not `C¹`.** The "close in `C¹`" quote is his index-1 passage (p. 271);
  the index-2 result is Theorem 4.5 (pp. 284–285), stated in `‖x‖ = ‖x‖_∞ + ‖(Px)′‖_∞` — derivative
  closeness on the **P-projected (differential) component only**. The pre-roll test as first stated
  (`‖ẋ_manufactured(0) − ẋ_orbit(0)‖` vs `k`) measures the wrong object; **project first**. And März
  explicitly permits an INCONSISTENT `x₀` ("no need for `x₀` to satisfy the second equation … but
  also the hidden constraint"), so "the pre-image is off the constraint manifold" is not a mechanism
  he licenses. Sharper prediction: `‖P(ẋ − ẋ_orbit)‖` grows with `k` while the value distance need
  not; if it is flat in `k` the `C¹_N` story is dead and constraint violation cannot replace it.
* **"Circuit MNA is Form-B-shaped" is the reviewer's inference, not Bereza's**; the thesis never
  mentions MNA. Bereza also does not flag index 2 as a known gap — he says only that his guarantees
  cover two index-1 forms. Nothing on either side predicts a state/monodromy gap FROM THE INDEX.
* **The `li_plus_rc` "refutation branch" is STRUCK — it was confounded.** That topology auto-enables
  `x0_unknown`, so trap's monodromy there was repaired before the comparison started. Held at
  `x0_unknown=False` it gives `5.053e-04, 1.253e-04, order 1.01` — bit-identical to plain RC.
* ⚠⚠ **The gap is measured on an ORDINARY INDEX-1 RC, on the default method**, against the exact
  pencil monodromy and the analytic RC transfer function (`VSin − R(1e4) − C(1e-6)`, trapezoidal):

        x0_unknown        state order   monodromy order   state err@100   mono err@100
        False (default)    2.01/2.00      1.01/1.00         3.36e-04        5.05e-04
        True               2.01/2.00      2.01/2.00         2.15e-02        5.01e-07

  Second-order state, FIRST-order monodromy; trap's monodromy is 10× WORSE than backward Euler's
  (`5.05e-05`). `x0_unknown=True` restores order 2 and a thousandfold at `N = 100` — at a **64× state
  cost**, which is why the unconditional default failed 13 tests (B1) and why "enable regardless of
  topology" is already measured as not free. Gear is immune (solved history, no opener) and refuses
  the flag with the reason. **The index has nothing to do with it; the manufactured opener does.**

**The decision, the user's:** for a one-step method, either enable `x0_unknown` whenever a
monodromy is REPORTED (`floquet_modes`, `ppv`, `info['Q']`, the Floquet paths) and keep the state
path as it is, or state plainly that `λ₂` from a manufactured opener is first-order. Today the number
looks like a second-order method's output and is not one.

#### B16-decision. ✅ **DECIDED AND BUILT 2026-09-05: the oscillator monodromy is Gear-2's, whatever method solved the state**

The user's instruction was "select the most accurate". Measured on the fixture that shows the gap
(`vdp + 0.3u²`, exact `Q_λ = 5.9083`, `c_true = 5.3703e-06`), the choice the question offered does not
exist — trapezoidal's own monodromy is unusable with EITHER opener:

    trap, default opener    Q_λ  11.1 / 28.4 / 63.9   at 400/800/1600   c/c_true 0.936 / 0.968 / 0.984
    trap, x0_unknown=True   Q_λ  3086 / 12228 / 48699                    c/c_true 1.022 / 1.011 / 1.006
    gear                    Q_λ  5.9094 / 5.9086 / 5.9084               c/c_true 0.998 / 0.9996 / 0.9999
    (trap's period and state are second order both ways: 1.2e-5 / 2.9e-6 / 6.7e-7)

The default opener's second multiplier DIVERGES toward 1 under refinement and `x0_unknown` puts a
spurious multiplier AT 1 — the one-step companion's parasitic mode, not the amplitude mode. So "most
accurate" is per quantity: **the state keeps the method asked for; every monodromy-derived quantity
comes from Gear-2 on the same orbit.** `PSS.monodromy_twin()` re-converges a Gear-2 `PSS` on the same
grid from the converged state (cached, a few warm iterations); `factored_period`, `ppv` and
`floquet_modes` delegate to it for an autonomous circuit under a one-step method; `pss.monodromy =
'native'` keeps the method's own factorisation for the gates that measure the plain path. Under trap
the bias fixture now reads `Q_λ = 5.9094` and `c` to `1.7e-3` at 400 points with its period untouched;
the trap oscillator covariance, which used to refuse with its reason, runs through the twin and matches
a direct Gear-2 solve to `1e-6`. An orbit too poor to seed the twin (Euler at 400 points: period 5% off,
amplitude 55% off) gets a `RuntimeError` with the reason. Driven circuits are unaffected: their plain
path is what the reference cross-check validated to six digits, and B16's first-order `λ₂` there is
recorded above as a known property of the manufactured opener, not repaired.

#### B16-preroll. A COMPUTABLE criterion for when a pre-roll may hand over to shooting — relayed 2026-09-05 at the user's request; ✅ **RUN the same day, results below**

De Luca, Bolcato & Schilders, "Proper Initial Solution to Start Periodic Steady-State-Based Methods",
IEEE TCAS-I 2019 (doi:10.1109/TCSI.2018.2874570; on disk under 07-shooting-methods). Instead of testing
closeness to the unknown `x*`, test the LINEARITY of the shooting error `u_k = x_k − φ(x_k)` during the
pre-integration: freeze `J_φ` at a candidate `k̂`, predict `u_{k+1} = J_φ(x_k̂) u_k` (their eq. 12),
measure `ũ_{k+1} = x_{k+1} − φ(x_{k+1})` (13), accept when `|u − ũ|_j ≤ ε_rel |u_k̂,j| + ε_abs` for all
`j` (16) on `n_iter` consecutive periods; a failure re-freezes `J_φ`. Their Algorithm 1 for `J_φ u` is
`FactoredPeriod.matvec` verbatim, so on this side it is the existing matvec pointed at the pre-roll,
running alongside the integration. Settings: four preliminary periods before checking, `n_iter = 7`,
`ε_rel = 1e-2`, `ε_abs = 1e-3`; found `k̂ = 4` on an RLC, an LNA at `n = 21` and an industrial LNA at
`n = 607`, matching manual tuning. Two caveats of theirs: the LTE degrades `J_φ` and the detector
(constrain `h`; their Fig. 3 is the noise floor arriving earlier with `h` unconstrained), and a
finite-difference `J_φ u` must reuse the same time points and method. ⚠ **Scope: the paper is
NON-AUTONOMOUS only** — known `T`, and `J_φ` without a unit multiplier. For an oscillator `u_k` converges
onto the phase direction rather than to zero. The review session's suggested adaptation, UNTESTED: apply
(16) to the phase-projected error `Π u_k`, `Π = I − u v₁ᵀ/(v₁ᵀu)` (the oblique projector A9 built), so
the amplitude part's linearity is what is detected. If `Π u_k` fails (16) deep into a converged
pre-roll the adaptation is wrong and the criterion is a driven-circuit tool — check the setup before
recording that (§D 0w). Sits beside März (the theorem, in the `C¹_N` norm, needing the distance to `x*`)
as the detector that never mentions `x*`, and replaces the guessed number of pre-roll periods Kundert
describes and the paper opens by criticising.

**RUN 2026-09-05, with pycircuit's own one-period map** (`_traverse_factored` gives both `φ(x)` and the
factored steps, so `J_φ u` is `FactoredPeriod.matvec` on the pre-roll state — the paper's Algorithm 1
verbatim; pair state under Gear-2, 400 points, their settings `ε_rel = 1e-2`, `ε_abs = 1e-3`, `n_iter =
7`, four preliminary periods; Algorithm 2's re-freeze on failure). Ground truth per period: does a
shooting Newton from that state converge in ≤ 6 iterations.

    fixture                                      detector fires at k     Newton converges from k
    DRIVEN, the paper's setting (stable lossy       raw 15  orth 15  obl 15          3
      nonlinear tank, 5% off resonance, seed 3 V)
    AUTONOMOUS, period guess 7e-6 off, seed 0.1 V   raw 81  orth 10  obl 66         48
    AUTONOMOUS, period guess 2.5% off               none in 90 periods              marginal (yes/no alternating)

Three findings. (i) **On its own ground the criterion works as advertised and is CONSERVATIVE**: it never
fired early, and fired twelve periods after the Newton could already have converged — a sufficient
condition, not a sharp one (their `k̂ = 4` on an RLC is under their `h` and tolerances, not a property
of the test). (ii) **The orthogonal-tangent adaptation is UNSAFE for an oscillator**: it fired at `k = 10`
with the amplitude at 0.15 V, thirty-eight periods before a Newton could converge — because the error
recursion IS linear there, around the unstable equilibrium the pre-roll is leaving. Linearity of `u_k`
detects linearity of whatever the local dynamics are, and the paper's implicit premise (a unique
attractor being approached) is what an oscillator started near its equilibrium violates. (iii) The
oblique projector built from the FROZEN Jacobian's own near-unit multiplier pair — which needs no `x*`
and no converged PPV, and which is only armed when `J_φ` actually has a multiplier within 0.3 of 1 —
fires at 66 against a truth of 48: late, but never false; the arming condition is what saves it. With
the period guess 2.5% off nothing fires: `u_k` is dominated by the phase slip `(T_guess − T*) ẋ`, and the
frozen Jacobian never carries a near-unit multiplier away from the orbit. So for oscillators the
usable form is "`J_φ` has a near-unit multiplier AND the obliquely projected error recursion is linear",
and it costs a dense eigen-pair of the frozen map — fine at `2m = 4`, a Krylov job at size. ⚠ **An
instrument failure on the way, caught by 0w before it was recorded:** the first driven control was the
bias-sensitive CORE under a 0.3 A drive, and its pre-roll blew up (`|u|` 1.8 → 92 in three periods); the
traversal was checked against an independent transient (`φ(x*) − x*` at `1e-16`) and cleared, and the
cause was the fixture — that core's `0.3u²` term overwhelms its `μ = 0.02` cubic at large excursions —
not the instrument. Not built into the tree: the criterion is a pre-roll policy, and the pre-roll
opener itself (the user's two-steps-back idea) is still the open B16 build.

### B9. Outer damped Newton — ✅ **ALREADY BUILT**, recorded so it is not re-requested

Requested 2026-09-04; it is in. All three `fsolve` calls pass `line_search=True`, and
`shooting.py` carries the rationale: *"the outer Newton is damped, which it was not … a departure
from standard practice rather than a neutral choice: Brachtendorf et al. describe 'shooting,
finite difference, or harmonic balance techniques in conjunction with a DAMPED NEWTON METHOD' as
what is widely employed for limit cycles. The full step is still tried first and kept whenever it
improves the residual, so a solve that was converging is unchanged; the halving only runs where
the undamped iteration would have moved uphill."*

⚠ **What is NOT damped is the INNER solve**, and that is where the high-`Q` failure lives: see
the B6-note tolerance floor, where `matrix_free=True` at `reltol = 1e-12` exhausts the inner
GMRES restarts at `Q ≥ 16` with exactly 126 matvecs regardless of the outer budget.

### B10. LSOAC — least squares with NO phase condition — NEW 2026-09-04, unbuilt

Mei & Roychowdhury 2006 DATE, relayed by the docs session: resolve the phase ambiguity by taking
**minimum-norm** solutions of the underdetermined system — a particular solution, then subtract the
null-space component — rather than adding a phase equation at all.

⚠ **THIS IS NOT WHAT C2 REJECTED, AND THE DISTINCTION IS THE ENTRY.** C2 tested one phase ROW
against another (orthogonality against the frozen-coordinate pin) and found the pin canonical and
the alternative 704× worse aligned when "fixed". LSOAC removes the row. Its stated motivation is
that phase conditions "cause various numerical artifacts" and that good ones "are not easy" to
choose — which is a claim about the whole family C2 lives in, not about which member wins.

⚠ **AND IT CONNECTS TO B3.** C3 already records that "A&T avoid it by having no phase equation at
all (see B3)", so this codebase has met the no-phase-condition idea once before, from Aprille &
Trick. B10 and B3 should be read together and probably costed together.

**Gate before building:** does the minimum-norm solve reproduce the shipped `λ₂` and `Q` on the
`m = 12` bulk fixture, and does it survive `λ₂ → 1` better than the pin — measured on the same
`Q`-sweep that produced B6's tolerance floor? If it is merely equivalent, it is not worth the
second formulation.

### B11. The index is DECIDABLE FROM TOPOLOGY — ✅ **BUILT 2026-09-04**, a diagnostic, not a refusal

Estevez Schwarz & Tischendorf (IJCTA 28(2):131–162, 2000), relayed by the docs session:

    the index of the DAE is 2 IF AND ONLY IF the network contains a C-V loop or an L-I cutset;
    otherwise the index is 1

for nonlinear time-independent networks without controlled sources, "assuming the positive
definiteness of the Jacobians of the element-characterizing functions", and extended in their
ref [26] to RLCTG networks (independent sources, resistive/capacitive/inductive subnetworks,
ideal transformers, gyrators).

⚠⚠ **THIS DOES NOT REOPEN C4, AND SAYING WHY IS THE POINT.** C4 closed *index-2
detect-and-refuse* because `index > 1` **is not predictive**: all three methods converge on an
LI-cutset and Gear-2 fails on 2 of 4. That finding is untouched. What the criterion changes is
that **the INDEX never needed measuring — it is decidable**, while **which INTEGRATOR converges
still does**. Knowing the index exactly and still not knowing which method to use is a *sharper*
result than not knowing either, and it retires the sloppier reading ("no method generalises across
index-2 topologies; measure per circuit") that conflated the two.

⚠ **WHAT IS WORTH BUILDING IS THE DIAGNOSTIC, BECAUSE THE CRITERION IS LOCAL.** The authors'
stated design goal is exactly this codebase's refusal-message problem: topological criteria "that
can be checked very fast", based on "LOCAL assumptions, i.e. we want to provide the opportunity to
LOCALIZE critical element modellings", motivated by circuits of ~1e7 elements where "it is often
difficult to find the circuit configurations that lead to numerical difficulties". So a failing
solve could say **WHICH elements form the offending loop or cutset** instead of "index > 1". That
is graph work on the netlist — no matrices, no solve — and it is a strictly better error message
than anything currently in the tree.

⚠ **AND ONE CASE WE MAY NOT COVER.** *"C-only loops have to be added to the class of C-V loops
since the currents through C-only loops belong to the network variables whereas these currents are
excluded in MNA formulations."* **A pure capacitor loop with no voltage source in it is index 2
and does not look like it.** The three shipped fixtures are named in the theory — CV-loop and
V-across-C are both C-V loops (the second the degenerate case), LI-cutset is an L-I cutset — but a
C-only loop is a fourth case and is not among them.

**Gate before building:** does the topological test agree with the measured index on all three
existing fixtures, and does it flag a newly built C-only loop that no current check catches?
Cheap, and it needs no solver.

---

✅ **BUILT 2026-09-04 — `topological_index(cir)` and `noise_enters_constraints(C, CY)`.** Both are
netlist/matrix-level checks with no solver in them. The index criterion is gated against a DIRECT
computation on the MNA matrices (`C`'s null basis `N`, then the rank of `NᵀGN`) over **seven**
topologies, and agrees on all of them.

✅ **THE C-ONLY-LOOP DISAGREEMENT IS RESOLVED, AND THE PAPER AGREES WITH THE MEASUREMENT.** The
docs session went back to the source: the sentence they relayed is from p.144, comparing against
their reference [10], and its second half says which formulation it is about — *"C-only loops have
to be added ... SINCE THE CURRENTS THROUGH C-ONLY LOOPS BELONG TO THE NETWORK VARIABLES WHEREAS
THESE CURRENTS ARE EXCLUDED IN MNA FORMULATIONS"*. **"In this case" is [10]'s formulation, which
carries capacitor currents as unknowns; MNA excludes them, which is exactly why MNA does not need
them added.** And the paper states it independently a page after Theorem 2.2: *"loops containing
only capacitances are EXCLUDED under point 4, whereas cutsets containing only inductances are
INCLUDED under point 3."*

⚠⚠ **THE CRITERION IS NOT SYMMETRIC, AND THAT ASYMMETRY IS NOW THE THING TO GET RIGHT.** A C-V
loop requires a voltage source; an L-I cutset does NOT require a current source — **L-only cutsets
count**. This implementation removes all `L` and `I` branches when testing connectivity, so
L-only cutsets are detected; the C-V search closes only on a voltage source, so C-only loops are
not. **Both halves match the paper.** So there is no split to record after all — the measurement
and the theory agree, and what disagreed was a relayed causal clause with its sign inverted.

⚠ **The grounded capacitor ring remains the cleanest demonstration:** `det C ≠ 0`, so it is not
index 2 and not even a DAE.

⚠⚠ **AND THE ORIGINAL DISAGREEMENT WAS STILL WORTH HAVING.** The quote said C-only loops
must be counted as C-V loops; a first version did, and disagreed with the measurement on **three
separate C-only topologies** — a grounded ring, a floating triangle, and the triangle with every
node resistively grounded, all measuring **index 1**. The arithmetic is checkable by hand: a
grounded ring has `det C = c1c2 + c1c3 + c2c3 ≠ 0`, so it is not even a DAE; a floating triangle
has `C` singular (its Laplacian) but `NᵀGN = (1/R)/3 ≠ 0`, so the constraint is uniquely solvable.
**A C-only loop makes `C` singular WITHOUT making the index 2** — index 2 needs a VOLTAGE SOURCE
fixing the loop. ✅ **CONFIRMED AGAINST THE PAPER, above:** the quote describes a formulation whose
variables differ from ours, exactly as the measurement implied.

⚠ **A UNION-FIND SUBTLETY THAT WOULD HAVE SHIPPED SILENTLY.** The first fix discarded any loop
found to contain no voltage source — wrong on a netlist carrying BOTH a C-only loop and a C-V
loop, since union-find reports only the FIRST closing edge and the C-only one can close first,
hiding the real one. Restructured to union capacitors first and close only on a voltage source, so
the loop is guaranteed to contain one. **Gated by a fixture carrying both**, which correctly names
`vs, c4` and not the `a/b/cc` triangle.

⚠⚠ **V-ONLY LOOPS AND I-ONLY CUTSETS ARE A DIFFERENT CATEGORY — ANDREAS ASKED, AND THEY ARE NOT
INDEX 2.** A loop of voltage sources over-determines KVL and a cutset of current sources
over-determines KCL: the MNA system is **structurally singular** and has no solution at all,
barring an exact cancellation. Calling that "index 2" would send a reader hunting a solver problem
when the netlist is the error. Reported as `info['ill_posed']` with the offending elements, and **`index` comes
back `None`** — the DAE index presumes a solvable system, so there is no honest value to give.

⚠ **ASKING WHETHER THE CHECK WAS THERE IMPROVED IT TWICE, WHICH IS WORTH RECORDING BECAUSE THE
FIRST VERSION *DID* PASS ITS TEST.** Demonstrating it rather than asserting it exposed two
defects the gate had not:

  * **the V loop named only its CLOSING source**, not the loop. On three sources in a ring that
    points at one and leaves the reader to find the other two — the opposite of the localisation
    the criterion exists for. Now shares the C-V loop's path reconstruction and names all three.
  * **and it reported `index = 2` with a "C-V loop" containing no capacitor.** A pure-V loop
    closes on a source, so the C-V search claimed it. Both symptoms point a reader at the SOLVER
    when the NETLIST is the error, which is precisely what this split exists to prevent.

    3 sources in a loop    index=None  ill_posed=True   V LOOP: v3, v2, v1
    I-only cutset          index=None  ill_posed=True   I CUTSET: i1, i2
    well-posed RC          index=1     ill_posed=False

⚠ **AND ONE FAILURE WAS IN MY REFERENCE, NOT THE CRITERION.** The first full run had `cv_loop`
disagreeing — criterion 2, reference 1 — because the reference guarded with `s2.max() > 0` and
`NᵀGN` for that fixture is **identically zero**, which is the MOST singular case rather than the
least. Shape 0j again, one day later: **the instrument was wrong, not the subject.**

### A10. Crystal oscillators (Q ≥ 10⁴) — ✅ **MEASURED 2026-09-04. THE BINDING LIMIT IS FREQUENCY ACCURACY, NOT Q**

Asked what it takes to support a crystal: `Q > 10⁴` for the motional arm, loaded lower by the
circuit but still far above anything previously tested here.

**Three things were candidates. Two are fine and the third is the whole answer.**

✅ **KRYLOV COST — FINE, and already settled by B13**: a 32× change in `Q` moves the GMRES count by
at most one. Iterations track the SLOW-NODE count, not `Q`. Loading the crystal only helps.

✅ **CONVERGENCE AND THE PPV AT `Q = 10⁴` — FINE.** Both methods converge, `ppv()` works, and
`|λ₂| = 0.99991` matches the analytic `exp(−μT) = 0.99990` to five digits.

⚠⚠ **AND `Q = 10⁵` IS NOT A CEILING EITHER — AN EARLIER READING HERE WAS WRONG.** A single
400-point run failed to converge at `Q = 10⁵` and was reported as a Q limit. It is a **GRID**
limit: at 1600 points both methods converge, and at 6400 `trap` reaches **0.0803 ppm**. Nothing
in the range tested shows a `Q` ceiling.

    Q = 1e5   400 pts    1600 pts            6400 pts
    gear      False      True  5.1468 ppm    True  0.3214 ppm
    trap      False      True  1.2867 ppm    True  0.0803 ppm

❌ **WHAT ACTUALLY BINDS: THE INTEGRATOR'S WARPING ERROR, AND IT IS INDEPENDENT OF `Q`.** Period
error at `Q = 10⁴`:

    method   400 pts    800 pts   1600 pts   3200 pts    rate
    gear     82.652     20.613     5.1468     1.2859    4.01, 4.00, 4.00
    trap     20.665      5.1533    1.2867     0.32146   4.01, 4.01, 4.00

⚠ **EXACTLY `O(h²)`, and `trap` is EXACTLY 4× BETTER THAN `gear` AT EVERY GRID** (82.652/20.665,
20.613/5.1533, 5.1468/1.2867, 1.2859/0.32146 — all 4.00). So `trap` ≡ `gear` at **half the
points**, which is a clean constant-factor statement rather than a trend. ⚠ And the error is
**constant across four decades of `Q`** (82.5 / 82.7 / 82.7 ppm at `Q` = 10², 10³, 10⁴), which is
what identifies it as the integrator rather than the physics.

**COST OF ppb, WHICH IS THE SPEC A CRYSTAL IS WRITTEN TO.** Second order, so from `trap`'s
0.0803 ppm at 6400 points, 1 ppb needs `√80.3 ≈ 9×` more — **≈ 57,000 points per period**; `gear`
needs **≈ 115,000**. Expensive but not absurd for a single period, and **refinement does get
there** — an earlier reading here implied it could not.

⚠ **THIS IS THE "WARPING ERROR" THE LITERATURE NAMES, and it is the one place our recorded
objection to higher-order methods is weak.** Brachtendorf-adjacent: Brambilla & Storti-Gajani,
TCAS-I 50:904 (2003) (*cited, not verified here*) characterise integration-induced `λ₂` bias as
*"equivalent to a perturbation of the eigenvalues of the linearized ordinary differential
problem"*, usually negligible — *"nevertheless an exception … is found when simulating
**high-quality factor circuits** where even very small warping errors can lead to qualitatively
wrong solutions"* — and conclude that *"higher order linear multistep methods, while characterized
by weaker stability properties, introduce **less** of a warping error and are **well suited** to
the simulation of high-quality factor circuits."*

⚠ **Our objection to higher order is a COST argument** (carrying variational history), **not a
numerical one**, and Gourary's Obreshkov single-step orders 1–4 carry no history at all. So the
objection has a published way round it. **Nothing here recommends building that** — it records
that the argument we were leaning on is weaker than it read.

⚠⚠ **A SECOND, STRUCTURAL ARGUMENT FOR A SINGLE-STEP HIGH-ORDER METHOD, from Kundert, White &
Sangiovanni-Vincentelli 1990 §4.2.5 on PARALLEL shooting** (read firsthand by the docs session):
*"At each step a high-order integration method needs the history of the solution over several past
time-steps. **This history cannot extend beyond a shooting interval boundary.** … When an interval
contains only a few time-points, high order methods lose their advantage because of the large
percentage of time-steps taken with the low order methods."* A single-step method carries no
history, so **full order is available on the first step of every subinterval** — the penalty
does not arise. That is an incompatibility, not an error-constant comparison, and it is the
stronger of the two arguments. ✅ The same section says parallel shooting *"increases the region
of convergence … as the number of subintervals increases"* (unstable modes cannot grow far over a
short interval) — the high-Q seeding basin attacked by shortening the interval rather than by
homotopy. **The two compose**: parallel shooting enlarges the basin but wastes high order; a
single-step high-order method restores the order without history. If the opener is ever fixed
(B16) and Obreshkov revisited, the parallel case is where it pays most. ✅ **And the recorded
"forced low-order restarts bias `λ₂`, −15.1 %" emulation is the one reconstruction in this
document the source CONFIRMS** — the passage above is its mechanism, stated in the foundational
text. Acquisition pointers, absent from the corpus: `skelboe80`, `smith87` (extrapolation
shooting); `keller68`, `keller76`, `stoer80` (multiple shooting).

⚠ **AND `info['Q']` DOES NOT REPORT A CRYSTAL'S DATASHEET `Q`** — see A9: it is the settling
rate `1/(2π(a − G))`, measured bit-identical across an 8× sweep of the component-set tank `Q`. A
designer reading `info['Q']` and expecting the motional `Q` gets an unrelated number.

**Not tested:** a real motional-arm + `C₀` model with a sustaining amplifier; startup; and
loaded-vs-unloaded `Q`. The fixture is a van der Pol at `μ = 1/(2πQ)`, which isolates `Q` and is
second order — see A9 on why second order is where a resonator `Q` is unambiguous.

### A9. Orbital (AM) noise and the far-out floor — ⚠ **THE PUBLISHED ANSWER IS IN OUR OWN LIBRARY**, 2026-09-04

⚠⚠ **THIS ITEM WAS SCOPED WRONG TWICE IN ONE DAY, BY TWO SESSIONS INDEPENDENTLY, AND THE
CORRECTION IS THE ENTRY.** It was reported as "nobody computes the far-out floor / the AM half is
unbuilt anywhere". The correct statement is: **we do not compute it, and a fifteen-year-old
principled method sits in our own `~/docs`.**

**F. L. Traversa and F. Bonani, "Including orbital fluctuations in the noise spectrum of
autonomous circuits", IEEE TCAS-I, 2011** — in `02-oscillator-noise-jitter/` under exactly that
title. ⚠ **Cited, not verified here**: nobody in this repo has read the paper. Arrived
independently from two sessions on the same afternoon, which is the only reason it is recorded
this firmly.

**The decomposition is the gap, verbatim:**

    x(t) = x_s(t + a(t)) + y(t)        a = phase,  y = orbital (amplitude)

and the two spectra are **summed**, not crossfaded:

> *"although the dominant component near to the oscillator frequency f₀ (and to its harmonics
> k f₀) is phase noise, **orbital fluctuations become the stronger contribution at large offset
> frequencies**"*

with *"evidence of its relevance for **high-Q oscillators**"* — **our regime**, and directly
relevant to the crystal question (A10), not an aside.

✅⚠ **AND IT IS PROBABLY A SMALLER BUILD THAN "PHASE + ORBITAL + CORRELATION" IMPLIES.** On the
paper's own worked example: *"the correlation spectrum is **negligible**, while orbital noise
becomes the dominant term for frequencies away from the harmonics"*. So the usable form looks
like **two terms added**, with the phase–orbital cross term dropped — materially less than the
three-term object this gap was described as needing.

⚠ **THE APPARENT 1 GHz / 5 GHz CONTRADICTION RESOLVED, AND THE RESOLUTION STRENGTHENS THE CLAIM.**
Two sessions described the worked example differently and I recorded both as UNKNOWN. Read from
the PDFs, **both were right about different papers**: the **IJMWT companion** is an InGaP/GaAs HBT
Gummel-Poon at **5 GHz, HB with 30 harmonics**; the **TCAS-I theory paper** is a Colpitts with the
same device model at **1 GHz, HB with 300 harmonics**. ⚠ **Quote the frequency and the harmonic
count together — the pair is what identifies which paper is meant.**

✅ **AND THE NEGLIGIBLE-CORRELATION RESULT IS REPORTED INDEPENDENTLY IN BOTH, ON TWO DIFFERENT
CIRCUITS** — TCAS-I: *"the correlation between phase and orbital noise (99), on the other hand, is
negligible"*; IJMWT: *"We found that the correlation spectrum is negligible, while orbital noise
becomes the dominant term for frequencies away from the harmonics."* Two circuits, not one.

⚠⚠ **STILL A DEFAULT TO TRY FIRST, NOT A LICENCE TO DROP THE TERM.** Same two authors on related
designs, so two circuits is not two independent confirmations.

✅⚠ **BUT THE CROSS TERM HAS A KNOWN SIGN, WHICH IS A BETTER STATEMENT THAN "NEGLIGIBLE"** (cited,
not measured here; read from p. 5 as an image): *"the approximate full normalized spectrum is
**lower** than the phase noise contribution, thus showing that **the correlation between the phase
and orbital deviations can decrease the total noise**. This effect is not present for `v = 0`,
since in this case the correlation spectrum is zero."*

So the cross term is **identically zero when there is no AM-to-PM coupling**, and when present it
**reduces** the total. ✅ **Dropping it therefore OVERESTIMATES noise — conservative for a design
margin, and wrong in a KNOWN DIRECTION against measurement.** That converts a two-term
implementation from "safe on one circuit" to "safe to ship with a documented bias", which is a
materially stronger position to build from.

⚠ **AND THERE ARE TWO Traversa & Bonani 2011 PAPERS**, which is how the confusion is most likely
to have arisen. The THEORY is in the TCAS-I paper; the MOTIVATION is in a companion — *Int. J.
Microwave and Wireless Technologies* 3(1):11–18, 2011, which carries the title used above. The
companion's conclusion is the one that matters for us: orbital noise *"becomes more significant
for **high-Q oscillators**, since its magnitude is … an increasing function of the Q factor"*.

**The ingredients are things we already compute or nearly do:** the PPV **is** the adjoint Floquet
eigenvector for the zero exponent; `ppv()` already returns `|λ₂| = exp(T·μ₂)`, the next Floquet
exponent; and the PSS waveform's Fourier coefficients are already consumed by
`oscillator_spectrum`.

⚠⚠ **"COST IS A FEW MORE FLOQUET PAIRS" WAS RELAYED WITHOUT CHECKING AND IS WRONG — corrected
2026-09-04.** Traversa & Bonani, IET CDS 2011, state the requirement flatly: the orbital and
correlation terms *"require the availability of **ALL** the direct and adjoint Floquet
eigenvectors"*. ✅ **What rescues it is their TCAD 2013 paper, not an assumption**: a chosen number
of exponents and both eigenvector sets, for the linearisation of **index-1 DAEs** around a limit
cycle — our formulation — with the error *"proved to tend to zero along with the ratio between the
norms of the neglected and retained rows"*. So truncation is legitimate **with that computable
bound, and only through that method**. `floquet_modes` now defaults to every non-null mode and
says so; a caller who truncates owes the norm ratio as the gate.

⚠ **THE BOUND, READ FROM TCAD 2013 (relayed; cited not verified here).** Their construction keeps
the first `m` rows of a matrix `R(t)` and neglects the last `n − m`, with small parameter
`ε = |smallest diagonal element KEPT| / |largest NEGLECTED|`, and *"both δμ_k and δũ_k tend to
zero **linearly** with ε → 0 … the error induced on the FEs and eigenvectors by the approximation
of R is at worst of the same order as the system approximation itself."* First order in `ε`, same
constant for exponents and eigenvectors.

⚠⚠ **THE DIRECTION IS THE OPPOSITE OF THE INTUITIVE ONE.** *"The larger the absolute value of the
neglected FEs, the smaller the error induced on the calculated FEs."* You drop the **fast** modes,
not the small ones, and accuracy improves the faster the dropped ones are. The step-2 residual
explanation — that it is the DAE's slaved algebraic directions — is exactly consistent: those are
the infinitely-fast end, which is the safe end to drop. That explanation was right; this is its
quantitative form.

⚠⚠ **CONSEQUENCE FOR THE API: `m` IS AN OUTPUT, NOT AN INPUT.** The paper does not take the
retained count as a parameter; it tests the condition per time sample and *"proceeds by changing
the m value to guarantee that the condition is met"*. So the principled form of `nmodes` is
*"the ε you will tolerate"*, with the routine returning however many modes that needs — possibly
differing between time points on one orbit. `nmodes=None` returning everything is the right safe
default; the ε-driven form is what makes truncation **legitimate** rather than merely permitted.

❌ **THE OPEN HALF, NOT SKIPPED:** `R(t)` is **not our monodromy and not `K_orb`**. It is a
specific object reached after *"a normalization procedure should be applied"* to the linearised
system, because *"the procedure is based on the estimation of the rank of matrix C(t)"*. The
ratio is between diagonal elements of **that** matrix in **their** normalisation. Transferring the
bound to `orbital_mode_weights`'s reconstruction needs the mapping established first, and it has
not been. The hand-set `1e-2` in the step-2 test stays until it is — replacing it with an `ε`
computed on the wrong matrix would be §D 0q again.

✅ **THE MAPPING IS NOW ESTABLISHED (relayed from the paper read as an image; NOT run).** `R(t)`
is the R factor of a **QR factorisation with column pivoting of the capacitance matrix** `C(t)`,
`C(t) = Q(t)R(t)Eᵀ` (their eq 3), after **row normalisation**. The recipe, with the three details
that are not optional:

  1. ⚠⚠ **NORMALISE ROWS FIRST OR THE DIAGONAL COMPARISON IS MEANINGLESS.** Each row of the
     linearised system is divided by the row norm of `S(t) = [dC/dt − A(t), (1/T)C(t)]` (eq 7),
     `‖S_j‖ = √((1/T) max_k ∫₀ᵀ S²_{jk} dt)` (eq 8), *"guaranteeing that all the vectors
     corresponding to the n discretized equations are versors"*. MNA rows carry wildly different
     physical scales; without this the diagonals compare a current row against a voltage row.
  2. **`E` is constant, from the FIRST sample.** *"The same E matrix is used for all time steps."*
     Re-pivoting per sample breaks their eq (4).
  3. The diagonal test proxies a row test because pivoted QR puts each row's largest element on
     the diagonal.

⚠⚠⚠ **EQUATION (9) IS THE REVERSE OF THE PAPER'S OWN PROSE, AND THIS WOULD HAVE BEEN IMPLEMENTED
BACKWARDS FROM THE PROSE.** The prose says *"the ratio between the smallest diagonal element kept
and the largest neglected"*. Equation (9), read as an image:

    max_{t∈]0,T]} |R_{m+1,m+1}(t)|   ≪   min_{t∈]0,T]} |R_{m,m}(t)|

So the **neglected** diagonal is the **small** one and the kept one is the large one —
`ε = max_t|R_{m+1,m+1}| / min_t|R_{m,m}|`, **neglected in the numerator, and a worst case over the
whole period on both sides**, not per-timepoint. A formula transcribed from the prose has the
grouping inverted. (§D 0m's pre-commitment rule, applied to a citation: name what the equation
says before accepting what the sentence about it says.)

**The identically-zero tail of `R(t)` is the nullspace of `C(t)`** — *"the last n − ρ rows of R(t)
are identically zero … correspondence between the nullspace of C(t) and the infinite FEs"*. That
is a **first, exact** reduction, separate from and prior to the approximate truncation `ε`
governs.

⚠ **SCOPE: THE METHOD ASSUMES INDEX 1** — *"according to the assumption of index-1 DAE, the rank ρ
of C(t) is time independent."* It does not cover the L-I cutset circuits of §0k. They also flag
that constant rank *"does not necessarily imply that the entire nullspace of C(t) is time
independent"*, which their method handles and their ref [16] does not.

⚠⚠ **A MAPPING SUBTLETY OF OURS THAT STILL BLOCKS BUILDING IT, found while recording this.** The
"null modes" `floquet_modes` drops on the **gear** path are **not** `C(t)`'s nullspace. Van der
Pol's `C` is 2×2 and full rank (capacitor voltage and inductor current are both differential), yet
the gear monodromy is 4-wide with two exact zeros — those are the **solved-history pair's**
artefacts, not infinite Floquet exponents of the DAE. On the plain (`trap`) path `n = m` and van
der Pol has no null modes at all. So the paper's *"structural zeros = nullspace of C"*
identification maps onto the **plain** monodromy directly and onto the **gear** one only after the
pair structure is accounted for. ✅ **Precondition before any `ε` is trusted: on a fixture with a
genuine algebraic node (a resistive-only node, so `C` really is rank-deficient), check that the
identically-zero tail of `R(t)` has exactly the rank the netlist predicts, on the plain path
first.** The peer session that supplied the recipe recommends the same gate and has not run it.

⚠ **THEIR OWN WARNING, AND IT CONSTRAINS OUR OUTPUT LAYER:** orbital contributions are
**asymmetric about the carrier**, so a symmetric single-sideband report **cannot carry them** —
see A4b. And *"the identification of the oscillator classes mostly impacted by this effect is not
an easy task"*: high-Q is a candidate, but eigenvector magnitudes matter as much as exponents.
**Do not present it as strictly better than what we ship without measuring it.**

⚠⚠ **A FREE PRIOR WE ALREADY RETURN — AND IT IS WEAKER THAN THIS ENTRY FIRST CLAIMED.** The
measurement below stands (`info['Q']` *is* the resonator `Q`), but the **use** proposed for it does
not follow: the same paper reports that `λ₂` does **not order the orbital contributions** —
eigenvector magnitudes can outrank exponents by more than the exponents differ (see §0 bound (ii)).
⚠ **So a `Q` in hand is NOT a sufficient prior on the far-out spectrum**, which is exactly what it
was proposed for. It remains a correct statement about settling and about the near-carrier items;
it is not the one-number-two-jobs shortcut it looked like an hour ago.

✅ **THE IDENTIFICATION ITSELF, MEASURED.** If orbital
noise magnitude increases with `Q`, then a `Q` already in hand is also a prior on **how badly a
phase-only spectrum reads FAR from the carrier** — one number doing two jobs, no new machinery.
`ppv()` does return `info['Q']`, built from the second Floquet multiplier as
`log(threshold)/log|λ₂|` — a SETTLING count in periods, where the `Q` in "high-Q oscillator" is
the resonator's. I predicted these were **different objects** and recorded the identification as
unverified.

⚠⚠⚠ **AN EARLIER VERSION OF THIS ENTRY SAID "MEASURED — THEY AGREE" AND THAT OVERCLAIMED. THE
FIXTURE CANNOT TEST THE IDENTIFICATION, BECAUSE BOTH SIDES COME FROM THE SAME `μ`.** For van der
Pol the amplitude envelope decays as `exp(−μT)` with `T = 2π`, so setting `μ = 1/(2πQ_target)`
gives `|λ₂| = exp(−1/Q_target)` and therefore `Q = −1/ln|λ₂| = Q_target` **identically**. Pure van
der Pol's linear part is undamped — it has **no independent resonator `Q`** to compare against.

✅ **WHAT THE RUN DOES ESTABLISH, which is worth having:** `ppv()` computes the settling count it
claims to, to 0.1 %, and the ratio landing at 1.00 rather than `ln(20) = 3.00` pins its threshold
convention as `1/e` **from data** rather than from reading the source. ❌ **What it does NOT
establish is "settling `Q` = resonator `Q`"** — §D shape 2, a number compared against itself.

⚠⚠⚠ **AND THE FIXTURE THAT CAN TEST IT SAYS `info['Q']` IS *NOT* THE UNLOADED RESONATOR `Q`.
MEASURED 2026-09-04.** Parallel `LC` tank, `L = C = 1`, with the two knobs **separated**: a linear
conductance `G` sets the **unloaded tank** `Q = 1/G`, and a nonlinear negative conductance `a`
supplies the loss back, so the amplitude relaxation goes as `(a − G)`.

**Hold the settling rate fixed at `a − G = 0.02` and sweep the tank `Q` across 8×:**

    G          0.010     0.020     0.040     0.080
    Q_res     100.0      50.0      25.0      12.5
    info['Q']   7.958163  7.958163  7.958163  7.958163

⚠ **BIT-IDENTICAL TO SEVEN DIGITS.** It does not move with the resonator `Q` at all. And holding
`Q_res = 50` while varying `a` moves it freely — 7.958 → 2.652 → 1.134. Across every row
`Q·(a − G) = 0.1592 = 1/(2π)`, so **`info['Q'] = 1/(2π(a − G))` exactly**: the settling rate, and
nothing else.

✅ **THE PRACTICAL STATEMENT, and it matters for A10:** a crystal's datasheet `Q` is the **unloaded
motional** `Q`, and **`info['Q']` will not report it**. Do not read one and expect the other.

⚠⚠ **AND IT RESOLVES A CONFLICT WITH A PEER RESULT THAT LOOKS LIKE A CONTRADICTION AND IS NOT.**
A peer session measured a **PASSIVE** resonator and found `Q_d / Q_λ → π` to six digits, concluding
the two quantities are *proportional with constant π* and that building an LC fixture "would
reproduce π — I would spend the time elsewhere". **Both results are correct; they are different
systems**, and the distinction is the whole point:

| system | amplitude decay per period | settling `Q` |
|---|---|---|
| **passive** resonator | `exp(−ω₀T/2Q_d) = exp(−π/Q_d)` | `Q_d/π` — **proportional to the tank `Q`** |
| **self-sustained** oscillator | set by the NET `(a − G)`, the active device having cancelled the loss | `1/(2π(a−G))` — **independent of the tank `Q`** |

Their π against the measurement above: at `Q_res` = 50 / 100 / 12.5 it predicts `info['Q']` =
15.9 / 31.8 / 3.98, and the measured value is **7.958 in all three**. ⚠ **π does not fit once the
circuit oscillates** — the prediction varies 8× where the measurement is flat.

⚠⚠ **AND THE OSCILLATOR CASE IS THE ONE THAT MATTERS, because every circuit PSS analyses is
self-sustained.** So the π conversion must NOT be applied to compare `info['Q']` against a
datasheet `Q`: it is right for a passive tank and wrong for the oscillator built around it. Taking
the peer's advice to skip the fixture would have shipped that conversion.

✅⚠ **AND WANG & ROYCHOWDHURY SAY THE SAME THING, IN TERMS — SO THERE IS NO CLAIM LEFT TO
ADJUDICATE.** An earlier version of this entry recorded that they *"argue it IS the energy `Q`"*
and left the question open pending their definition. That attribution was wrong. From their p. 2
(relayed verbatim by a peer session that went back to the paper; *cited, not verified here*):

> *"While the frequency- and energy-based definitions of `Q` for second-order linear resonators
> are equivalent, **different `Q` definitions for oscillators are not**."*
> *"We emphasize that the proposed `Q` factor formulation is indeed **different** and more
> suitable for characterizing amplitude stability."*
> *"the widely-used energy-based `Q` formulation in (2) for **linear** systems is in fact just a
> special case of our definition (4), with the amplitude-stable state being the zero state."*

⚠ **The containment runs the OTHER WAY from how it was first recorded:** the classical energy `Q`
is a *special case of theirs for linear systems*, and **for oscillators they state the definitions
diverge deliberately.** So the measurement above is not in tension with the paper — it is what the
paper predicts.

**THE THRESHOLD ALGEBRA, exact, for the passive case only:** `Q_th = −Q_d·ln(th)/π`, so the ratio
is a pure function of the reporting threshold —

    threshold  1/e     = 0.367879   →   Q_d/Q = 3.141593    ← ours
    threshold  0.05                 →   Q_d/Q = 1.048689    ← Wang & Roychowdhury's
    threshold  e^(−π)  = 0.043214   →   Q_d/Q = 1.000000    ← would make them identical

⚠ **W&R's 5 % lands within 4.9 % of the threshold that would make the two coincide, which is why
their agreement READS as an identity and is not one.**

✅⚠ **THREE CASES, NOT TWO — AND THE THIRD RETURNS TO π FOR A REASON.** The ratio is a property of
the **core's compression**, not of the formulation. All three measured or derived-then-measured
here:

| core | `Q_λ` | `Q_d/Q_λ` |
|---|---|---|
| **passive** resonator | `Q_d/π` | **π** |
| **soft-compressing** (cubic) | `1/(2π(g − α))` | `2π(g − α)/α` — a **design quantity** |
| **amplitude-clamping** (hard limiter) | `1/(πα)` | **π**, drive-independent |

**The clamping case was relayed as an assertion and is now derived and measured.** A hard
limiter's describing function is `N(A) = 4a/(πA)`; at equilibrium `4a/(πA) = α`, so
`dN/dA = 4a/(πA²) = α/A` and the amplitude rate is `(A/2)(α/A) = α/2` — **independent of the
drive `a`**, because the equilibrium condition pins the slope. Over one period `exp(−πα)`, giving
`Q_λ = 1/(πα)` and ratio π.

    G      a      k·A     info['Q']    pred 1/(πG)   Q_d/Q_λ
    0.02   0.02    25.4   15.935730     15.915494    3.13760
    0.02   0.04    50.9   15.920465     15.915494    3.14061
    0.02   0.08   101.9   15.916669     15.915494    3.14136
    0.01   0.02    50.9   31.841200     31.830989    3.14059
    0.04   0.04    25.4    7.966809      7.957747    3.13802
    0.08   0.08    25.5    3.981092      3.978874    3.13984

⚠ **THE APPROACH TO π IS MONOTONE IN CLAMP HARDNESS** — 3.1376 → 3.1406 → 3.14136 as `k·A` goes
25 → 51 → 102 — which confirms the *mechanism* and not merely the number, since the derivation
assumes a hard limiter. And `info['Q']` tracks `1/(πG)` while **ignoring `a`**: the exact opposite
dependence from the cubic core, where it tracked `1/(2π(a − G))` and ignored `G`.

⚠⚠ **SO "MULTIPLY BY π" IS RIGHT FOR TWO OF THE THREE AND WRONG FOR THE ONE IN BETWEEN**, which is
the soft-compressing core most transistor oscillators actually are. A constant that holds at both
ends and fails in the middle is the worst possible shape for a rule of thumb.

⚠⚠⚠ **AND THE SOFT-COMPRESSING SWEEP CONTAINS A COINCIDENCE POINT WHERE THE WRONG RULE IS EXACTLY
RIGHT.** `2π(a − G)/G = π` when `a = 1.5G`, which at `a − G = 0.02` is `G = 0.040` — **the third
row of the measured table above**. ⚠⚠ **AND THE COINCIDENCE IS ALGEBRAICALLY EXACT, NOT APPROXIMATE:** at `a = 1.5G` both
expressions are literally `25/π` — `Q_res/π = (1/0.04)/π` and `1/(2π(a−G)) = 1/(0.04π)`. The
*measured* 3.14142 against `π = 3.14159` differs only by the grid's own error. ⚠ **A SINGLE-POINT
CHECK AT THAT ROW WOULD NOT HAVE RETURNED A PLAUSIBLE WRONG ANSWER — IT WOULD HAVE RETURNED AN
EXACT ONE**, and the row was in the sweep by accident rather than design. The general form:
**a law that matches at one point has been tested at one point.**
**The sweep is what carries this result, not the fixture** — spot-checking a proposed constant at
one operating point cannot distinguish "constant" from "passes through that value here". (Spotted
by a peer session reading the table, not by the session that produced it.)

⚠⚠ **TWO CONSTANTS, NOT INTERCHANGEABLE.** `ln(20) = 2.9957` converts between the `1/e` and 5 %
**reporting conventions** — a choice. `π` converts a **passive** resonator's settling `Q` to its
energy `Q` — a per-period decay. ⚠ **Neither converts an oscillator's `info['Q']` into anything
about its tank**, because the measurement above shows no such information is present.

✅✅ **CLOSED 2026-09-04 BY READING THE PAPER — AND THE ANSWER IS "IT DEPENDS ON THE REGIME",
WHICH IS WHY THE TAXONOMY ABOVE IS THE THING THAT SETTLES IT.**

Read from the IJMWT companion (`02-oscillator-noise-jitter/`, extracted with `pdftotext -layout`;
these passages are in the body text, not the dropped-font equations):

> *"This circuit was chosen since it can be shown to be equivalent [6] to a parallel RLC circuit
> with a **Q factor corresponding to the Q coefficient** in Fig. 6: this allows for a simple
> modulation of the Q value to study its impact on orbital noise."*
> Fig. 7 caption: *"Floquet exponent μ₂ … as a function of the **Q factor of the equivalent RLC
> circuit**."*

⚠ **SO THEIR `Q` IS THE TANK `Q`** — the component-set RLC quality factor, not a settling count.
On the bare reading that would make `info['Q']` the wrong object and the link void.

✅⚠ **BUT THEIR CIRCUIT IS IN THE CLAMPING REGIME, WHERE THE TWO ARE PROPORTIONAL.** Their inverter
is *"approximated by the input–output relation `v_out = tanh(2a·v_in)`, where `a` is a parameter
representing the slope"*, run at **`a = 23`** — slope 46, a hard clamp. That is **row three of the
taxonomy above**, measured here at ratio π with `Q_λ = 1/(πα) ∝ Q_d`. And they say exactly this:
*"a high-Q oscillator is characterized by at least a second Floquet exponent near to zero, which
in turn should result into a larger amplitude noise component"* — which is `Q_λ ∝ Q_d`, the
clamping relation, not a general one.

✅✅⚠ **AND THERE IS A DIRECT ROUTE THAT NEEDS NO TANK `Q` AT ALL, WHICH IS BETTER THAN THE
REGIME ARGUMENT AND SUPERSEDES IT.** T&B's *knob* is the tank `Q`; their *stated mechanism* is the
multiplier:

> *"a high-Q oscillator is characterized by at least a **second Floquet exponent near to zero**,
> which in turn should result into a larger amplitude noise component."*

Since `λ₂ = exp(μ₂T)`, **`info['Q'] = −1/ln|λ₂| = −1/(μ₂T)`** — so a large settling `Q` **IS**
`μ₂` near zero. Not a proxy for their mechanism: **the same quantity as their mechanism**, up to
the period. Verified against the measured table: at `Q_target = 8`, `|λ₂| = 0.8825206810` gives
`−1/ln|λ₂| = 8.0014` against `info['Q'] = 8.001725`.

✅ **So `info['Q']` is a legitimate and DIRECT prior on orbital noise, in EVERY regime, and no tank
`Q` is needed at any point.** The regime taxonomy governs the *tank-`Q`* route only.

**THE DIRECTION TABLE — this is what the fixture actually bought:**

| direction | verdict |
|---|---|
| settling `Q` → *"`μ₂` is near zero, expect larger orbital noise"* | ✅ **SOUND** — their own mechanism, and an identity |
| settling `Q` → tank `Q` | ❌ **VOID** — measured bit-identical across an 8× sweep |
| tank `Q` → orbital noise | ⚠ holds in Tow-Thomas, **not in general** — needs the clamping regime |

⚠⚠ **THE FIXTURE DID NOT KILL THE LINK. IT KILLED THE WRONG DIRECTION OF IT** — which is the one
that had been written down, and the one that would have been used silently.

⚠ **ONE THING TO KEEP MARKED:** T&B *show* `Q → μ₂` (Fig. 7) and `Q → noise` (Figs 8–10, at 1, 10
and 100 Hz offsets). That `μ₂ → noise` holds *independently of `Q`* is an **inference from their
data plus their stated mechanism**, not something they isolate. Their fixture *"has two state
variables, therefore only two Floquet exponents"*, so it is genuinely the second-order case their
own caveat names.

### A9 step 3 — `C_lhj` assembled and GATED three ways; a defect in step 1 found on the way — 2026-09-04

**The gate.** Eq (23): `Σ_{l≥2,h,j} C_lhj = R∞_yy(0)`, the stationary transverse covariance. The
reference is `oscillator_covariance`'s Lyapunov solve, which shares nothing with the modal sum.
Prototyped in the scratchpad first; nothing shipped until the identity held.

⚠ **THE REFERENCE HAD TO BE THE CYCLE-MEAN, NOT `K_orb(0)`.** Lemma 3.5's `R∞_yy` "depends on τ
only" — it is the *stationary part*, i.e. the cycle average. `K_orb(0)` on van der Pol is all on
the voltage (at `t = 0` the orbit sits at `[2, 0]`, so the amplitude direction is pure-v) while
the modal sum is isotropic — which is exactly a rotating radial direction averaged over a cycle.
Comparing to `t = 0` was the wrong time reference. The per-sample transverse part is
`P(t_j) − (t_j/T)·growth_samples[j]`; with the along-orbit growth removed the transverse trace is
**flat around the cycle** (6.275–6.293e-6), and its cycle mean is isotropic at exactly half the
radial variance. That is the reference.

⚠⚠ **THE GATE FAILED BY 1.75×, AND LOCALISING IT FOUND A REAL DEFECT IN `floquet_modes`.** A
**third route** — `R∞_yy(0)` from its definition, a 1-D Lyapunov integral along the single orbital
mode with **no Fourier sum** — matched eq (22)'s sum to **3.5e-4**, isolating the discrepancy to a
*shared input*. It was the adjoint's scale: `v_k` was biorthonormalised on the width-`n` **pair**
vectors, then sliced to the width-`m` state block, leaving `q(0)ᵀp(0) = 1.324143` there — constant
around the cycle to four digits (the invariant is preserved by the flow, which is itself a check)
— and `q` enters the covariance quadratically: `1/c₀² × 1.7535 = 1.0001`. **The periodicity gate
`p(T) = p(0)` could not see it: periodicity is scale-free.** Fixed by renormalising `q` on the
state block (a no-op on the plain path, where `n = m`); the test now asserts `q(t)ᵀp(t) = 1`
around the cycle. ⚠ I named the number before reading it — `c₀ ≈ 0.755` — and got the
**direction** wrong (`c₀ = 1.324`); the mechanism was right. The pre-commitment still did its job
once the sign of the guess was corrected.

✅ **AFTER THE FIX, THREE ROUTES AGREE IN MAGNITUDE TO < 1e-3:**

    eq (22) modal sum  vs  cycle-mean Lyapunov   ratio 0.99976
    definition integral vs  cycle-mean Lyapunov   ratio 1.00008
    definition integral vs  eq (22) modal sum     ratio 1.00031, residual 3.5e-4

which validates the transcription of eq (22), the **two-sided `CY/2`** convention (consistent with
the `kT/C`-calibrated Monte Carlo injection `Var(i) = CY/(2h)` already in the record), and the
renormalisation, together.

✅ **SHIPPED as `PAC.orbital_correlation(pss, H)` → `(R, C_lhj)`**, gated by
`test_orbital_correlation_is_gated_three_ways`, which builds the definition-route integral inside
the test as the independent check. **Stationary white sources only — and that is a scoping
DECISION, not a gap: coloured noise is deferred by Andreas (2026-09-04, "we will add coloured
noise later").** With `CY` constant the `Λ̃` products collapse to `Ṽᵀ CY Ṽ*` and `B` is never
needed; a coloured source breaks that collapse and needs `B(t)` per source, which is the shape of
that later work.

❌ **OPEN, NOT TUNED: a 3 % SHAPE residual** between the two modal routes (which agree with each
other to 3.5e-4) and the Lyapunov reference — the scalar-fit residual is 2.98e-2 before and after
the fix, so it is not a factor. Attributed, **not proven**, to the reference being the **pair**
covariance on gear sliced to its state block. The clean comparison needs the **plain-path**
Lyapunov solve — `_lyapunov_pieces` is still gear-only, the one adjoint surface B8's wiring did
not reach. That is the next precondition, not a tolerance to widen. ⚠ **One named thing to RULE OUT for
that residual, offered by the docs session as a search term and explicitly NOT as a hypothesis:**
Traversa & Bonani 2013's *"Floquet eigenvalue split set"* — replicas of the true exponents shifted
in imaginary part by integer multiples of the fundamental. If a pair formulation sliced to its
state block admitted a replica alongside the physical mode, it would produce a shape error at the
correct magnitude, which is the signature seen. Nobody has checked whether gear's pair map does
this; it is recorded so the thing has a name, not because anyone believes it yet.

✅⚠ **THE PLAIN-PATH LYAPUNOV SOLVE IS NOW WIRED (2026-09-04, `_lyapunov_pieces_plain`), AND IT
FALSIFIES THE ATTRIBUTION ABOVE.** On **euler-plain** — `n = m`, no pair, nothing to slice —
`orbital_correlation` against the cycle-mean transverse Lyapunov covariance gives magnitude
**0.99993** and a shape residual of **2.4–2.7 %, flat between 1600 and 3200 points**. So the
residual is neither the pair slice nor discretisation. It stays open, bounded at 5 % in the
test with a floor asserting it has not silently vanished. **Candidate, UNTESTED, that fits its
size and the theory:** eq (22)'s `l ≥ 2` sum is the *pure orbital* term; Theorem 4.1's total
adds the **phase–orbital correlation** `S_corr` (eqs 18/20, `D_lhj`), whose `τ = 0` value the
Lyapunov covariance contains and this sum excludes — reported by the authors as small, with a
sign. The test that settles it is assembling `D_lhj` and adding its `τ = 0` contribution.

✅✅ **CLOSED THE SAME EVENING: THE RESIDUAL WAS THE REFERENCE, NOT THE SUM.** Two things arrived
together. The docs session **disproved the correlation candidate from the equation** before it
was built: at `τ = 0` eq (18a)'s brace is `{exp[0] − exp[0]} = 0` identically, for every `l, h, j`,
and (23) states `R_yy(0) = Σ C_lhj` alone — there is no correlation term in a `τ = 0` covariance.
And the remaining suspect was the reference: subtracting only the **secular** growth
`(t/T)·d·uuᵀ` leaves the phase direction's **bounded within-period variance**, which the `l ≥ 2`
sum correctly excludes. Demir's `y` is the **oblique projection** `v₁ᵀy = 0`. Prediction named
before running: projecting the samples with `Π = I − uvᵀ/(vᵀu)` drops the residual to quadrature
(~1e-3). Measured, euler-plain, `n = m`:

    Q_λ            4.44      9.97     19.94     39.87
    growth-subtracted   5.37e-2   2.39e-2   1.20e-2   5.99e-3     ← falls as 1/Q_λ
    oblique-projected   5.9e-4    6.0e-4    3.1e-4    1.5e-4     ← quadrature, improving

✅ **Three routes now agree to ~1e-4**, and the peer's Q-sweep discriminator — physics grows with
`Q`, numerics is flat — returned a *third* outcome, **1/Q**, which is the same fact: orbital
variance `∝ Q_λ` (`C_lhj ∼ 1/μ₂`) against a constant phase-bounded part. ⚠ `K_orb` is
**bounded, not transverse**: read it as the bounded part, and project it if `R_yy(0)` is wanted
(docstring updated). The tests now compare against the projected reference at 3e-3 and keep the
growth-subtracted 2.4 % as the documented signature.

⚠ **THE 1/Q LAW IS AN INVARIANT, checked by the docs session rather than taken:** `r · Q_λ` =
0.2384 / 0.2383 / 0.2393 / 0.2388 — **constant to 0.4 % of its mean over a 9× range in `Q`**, fitted
exponent −0.9988. The exponent is safe; **the constant is vdp-specific until a second oscillator
family is run**, and the number should not be read as meaning anything before that. (The projected
column fits −0.645, which is a floor, not a law — correctly not quoted as one.) The projector
algebra was also checked independently: `‖Πu‖ = 2e-16`, `‖vᵀΠ‖ = 7e-16`, `‖Π² − Π‖ = 2e-15`.

⚠ **AND A CORRECTION TO HOW THE SWEEP WAS CREDITED, which is §D 0v.** It was recorded as "returning
a third outcome". That is the test's error, not a bonus: it offered *grows with Q → physics, flat →
numerical* on a **ratio** whose numerator and denominator scale differently — a dimensionless
residual can *fall* while the physics it measures *grows*. Read literally, "falling" is neither
outcome and the honest conclusion would have been "the test does not apply"; the right answer was
reached only because the projection had already found it from the inside.

**What the wiring is.** Euler (`b = 0`): step state `x` alone, `n = m`, consumers untouched.
Trapezoidal (`b = −1`): step state is the pair `(x, iq)`, the maps are `2m × 2m`, the noise
reaches `iq` through `a₀C_kK` as well as `x` through `K`. ⚠⚠ **The un-reset trap pair is
SINGULAR, measured** — carrying `iq` across the boundary puts a marginal `(−1)ⁿ` mode in
`I − M⊗M` (`LinAlgError` on a driven RLC): the obstruction this document records for **every**
formulation that keeps the companion across a period, met for the fourth time, from the
covariance side. The period map re-seeds `iq` (`M = (ΠA)·diag(I, 0)`), as the solve does. **Tie
gate:** the product of the per-step maps reproduces `fp.matvec` to 0 (euler), 1e-12 (trap, on
`(x, 0)` read out on `x`), 0 (gear). **`kT/C` gate:** `covariance` under euler and trap converges
on 1.0 exactly as gear does (0.969 / 0.969 / 0.955 at 1600 points). ❌ `oscillator_covariance`
is **refused for trap-plain with its reason**: it borders with `ppv()`'s width-`m` vectors and the
trap pair map is `2m × 2m`; the pair's own null vectors would be needed, with a normalisation
this record was burned on twice today. Euler-plain and gear both work.

⚠ **A FIXTURE CATCH, PREDICTED BEFORE IT WAS READ.** The first `kT/C` gate read **2.0** under
every method. The prediction — "if `R` already carries thermal noise, an explicit
`IS(noisePSD = 4kT/R)` beside it double-counts, and gear will read 2.0 on the same fixture too"
— held exactly: with `IS` all three → 2.0 (1.939 / 1.939 / 1.911), with `R` alone all three →
1.0. The code was right; the fixture was mine.

⚠ **AN EARLIER VERSION OF THIS ENTRY CONCLUDED "the link holds in the clamping regime and is void
in the soft-compressing one".** That is correct about the tank-`Q` route and **generalised too
far** — it treated the only route it had found as the only route there is.

**The table, read for what it can support:** van der Pol at `μ = 1/(2πQ_target)`, 400 points,
gear:

    Q_target      2.0        8.0       32.0      128.0      512.0
    info['Q']   1.998794   8.001725  32.017759 128.222454 515.346075
    ratio        0.9994     1.0002    1.0006     1.0017     1.0065

⚠ **The ratio DRIFTS MONOTONICALLY UPWARD** — 0.65 % high at `Q = 512`. ⚠ **AND THE OBVIOUS
READING OF THAT DRIFT IS ALSO WRONG.** It is not the identification degrading; it is the
**formula's conditioning**. `Q = −1/ln|λ₂|` has `dQ/Q = Q·ε` for a fixed *absolute* error `ε` in
`|λ₂|`, and `|λ₂| → 1` as `Q` rises. Inverting the table for the implied `ε`:

    Q_target        2         8        32       128       512
    implied ε   1.83e-04  2.38e-05  1.68e-05  1.35e-05  1.27e-05

**Essentially constant from `Q = 32` upward**, while the relative error in `Q` grows linearly in
`Q`. So the quantity to distrust at crystal `Q` is **`|λ₂|`**, and `Q` merely inherits its error
multiplied by `Q`. (Diagnosis relayed by a peer session; the inversion is arithmetic on the table
above.)

✅ **CONFIRMED BY REFINEMENT AT FIXED `Q = 512`** — the test that separates "conditioning" from
"breakdown". The drift **falls with `h`**, so it is discretisation in `|λ₂|` being magnified and
the formula is sound:

    npts        400        800       1600       3200
    rel err   6.54e-03   8.27e-04   1.08e-04   1.48e-05
    implied ε 1.28e-05   1.62e-06   2.11e-07   2.88e-08

⚠ The rate is ≈ 7.9/7.7/7.3 per doubling rather than the 4 a second-order method would give on a
generic functional — recorded as measured and **not explained here**; `|λ₂|` may converge faster
than the trajectory does. Do not build on the rate without accounting for it.

**The supporting paper is also in the library** (cited, not verified here): T. Wang and
J. Roychowdhury, *"Rigorous Q Factor Formulation and Characterization for Nonlinear Oscillators"*,
arXiv:1710.02015, `09-phase-macromodels-and-prc/`. It defines `Q` from exactly the object `ppv()`
computes and argues it **is** the energy `Q` rather than a correlate.

⚠ **TWO CAVEATS, EITHER OF WHICH PRODUCES A PLAUSIBLE WRONG ANSWER.** (a) **The threshold
convention differs by ≈ 3**: `ppv()` reports cycles to `1/e`, Wang & Roychowdhury to 5 %, and
`ln(20) = 3.00`. The measurement above lands at ratio ≈ 1 rather than ≈ 3, which independently
confirms `ppv()` uses the `1/e` convention — compare across conventions without that factor and a
threefold disagreement looks like a bug. (b) **Traversa & Bonani's own claim is narrower than it
sounds**: *"at least for the second-order oscillator considered in this study, an increasing
function of the Q factor"*. Second order is where a resonator `Q` is unambiguous — and the fixture
above is second order, so the measurement confirms it exactly where the theory says it holds and
**not one step further**.

✅ **AN INDEPENDENT CHECK ON SOMETHING WE ALREADY SHIP.** Kundert, *Introduction to RF
Simulation* §3.5, gives the swept small-signal result's validity window as `f_Δ ≪ Δf ≪ f₀`,
states it is in error below it, and shows **how to obtain `f_Δ` from the swept result itself**.
`PAC.phase_psd` predicts the corner in closed form — `f_h = π i² f₀² c`, `shooting.py:7419` — and
**refuses** below it.

⚠⚠ **AN EARLIER VERSION OF THIS PARAGRAPH SAID THAT PREDICTION "HAS NEVER BEEN CHECKED AGAINST
ANYTHING". THAT WAS WRONG — checked 2026-09-04 by reading the tests rather than assuming.** Three
gates already stand:

  * `test_the_lorentzian_conserves_the_carrier_power_exactly` — `∫S_i df = 1` to 1e-6, integrated
    **numerically over the implemented function** (not re-derived) for 3 values of `c` × 3
    harmonics. This is Kundert's own eq (19) invariant.
  * `test_the_lineshape_is_lorentzian_where_it_should_be` — peak `1/(π²i²f₀²c)`, the `1/f²` skirt
    asserted **as a rate**, and the half-width `π i² f₀² c`. ⚠ That half-width assertion is
    algebra against the same algebra, so it is self-consistency rather than an external check.
  * `c` itself is gated against a **Monte Carlo**, with a `kT/C` reference recorded as having
    caught a factor-of-two that a same-assumption Monte Carlo had confirmed.

✅ **AND THE CLOSED FORM AGREES WITH KUNDERT'S, READ FROM THE SOURCE** (`07-shooting-methods/`,
*Introduction to RF Simulation and its Application*, p. 11): *"The corner frequency fΔ is known as
the **linewidth** of the oscillator and is given by **fΔ = cπfo²**"* — identical to ours at the
fundamental, and our `i²` is the `i`-th harmonic's scaling.

⚠ **WHAT KUNDERT'S "GET IT FROM THE SWEEP" ACTUALLY IS, which is not what it was relayed as.** He
does not give a separate corner formula. He says the small-signal analysis *"does not show the
roll off"* but that *"it is possible to use (15) to determine fΔ"* — i.e. **fit `c` from the far
`1/Δf²` skirt of the swept result and apply the same formula.** So the genuinely independent
check is not on the corner formula but on **`c` itself, by two code paths**: the swept noise path
against `diffusion_constant()`. ✅ **THAT was the unrun item, and it is now RUN — the two paths agree to seven figures.**

    Δf/f₀      c from the swept skirt        c from diffusion_constant()
    1e-2       7.510951726e-08               ← outside Kundert's window
    3e-3       6.381291626e-08
    1e-3       6.263389371e-08
    3e-4       6.251214799e-08
    1e-4       6.250576334e-08          ←→   6.250576786e-08

The routes share the PSS and the noise sources and **nothing else** — path A contracts the PPV
against `CY`; path B propagates sidebands through the adjoint rows and normalises to carrier
power. Agreement at `Δf/f₀ = 1e-4` is **7e-8 relative**. ✅ And `P_carrier = 2.000009` against the
analytic `A²/2 = 2`, which independently re-confirms the carrier-power normalisation from §0i.

⚠⚠ **THE ESTIMATOR HAD TO BE THE LIMIT, NOT AN AVERAGE — AND THE FIRST VERSION GOT THAT WRONG.**
Taking the *median* across all five offsets reported the two paths agreeing to 0.2 %, which is not
a measurement of their disagreement but of **how many invalid offsets were included**. Kundert
states the window as `fΔ ≪ Δf ≪ f₀`; the outermost point is 20 % high because it is outside it.
**An average over a validity boundary reports the boundary, not the quantity.**

Pinned by `test_c_agrees_between_the_ppv_form_and_the_swept_noise_path`, which asserts both the
tight agreement inside the window and that it *improves* as the window is entered — the second
being what distinguishes a window effect from a constant offset.

⚠ **WE ARE NOT EXPOSED TO §10 OF THE CITATION MAP.** It warns that a crossfade weight,
`√(1.5 − |f|/W)` band edges and a Lorentzian cap are implementation choices with **no literature**
behind them. Checked: **`shooting.py` contains no crossfade, blend or band-edge construct at
all.** What it has is the corner *refusal* above — a switch with a stated validity window, which
is exactly what Demir §X.C and Kundert (15) support. Nothing to disown here.

⚠⚠ **THE FAILURE SHAPE THAT HID IT — §D 0p.** The paper was in the library, had been **opened**,
and was **cited for a different claim it does not own** (the `Φ = U·D·V·C` state-transition
factorisation, which belongs to Demir 2000 eq (37)). Its *actual subject* was then reported as an
open gap hours later. **Reading a paper for one claim and filing it under that claim is how a
library loses a result it already contains.** The tell available in advance: a citation whose
title does not match the claim it is attached to.

⚠ **METADATA, so nobody trusts the wrong field.** The Gourary regularisation cluster is three
papers — **DATE 2003** (periodic small-signal), **ECCTD 2007 pp. 1002–1005** (time-domain
oscillator noise, measured and rejected as B12), and **37th EuMC Munich, October 2007**
(cyclostationary). Their filenames carry no author or venue, against the library convention, and
the EuMC paper's embedded PDF metadata is **placeholder junk dated 1999** — do not read the date
off the file.

### A8. Sampled / edge-jitter noise (`noisetype=timedomain`) — NEW 2026-09-04, unbuilt

⚠ **FROM A REAL CIRCUIT, NOT FROM THE LITERATURE.** Andreas put the case up: a free-running
oscillator followed by several inverter buffers, output taken at the **last inverter**, the
inverters fed from a **noisy LDO**.

⚠⚠ **AND HIS CORRECTION IS WHAT MAKES THIS A SMALL ITEM.** The first reading of it reached for
envelope-following and MPDE on time-scale grounds. That was wrong: the LDO is characterised
separately and injected as a **supply-node noise source**, and inverter delays are a few percent
of the period, so this is **one ordinary autonomous PSS over osc + buffers at a few hundred
points**. No multirate anything. It moved from *"needs machinery we do not have"* to **"needs an
analysis we do not have, on a netlist we can already solve."**

**THE TWO JITTER MECHANISMS ARE PHYSICALLY DIFFERENT AND ONLY ONE IS COVERED:**

| mechanism | source | status |
|---|---|---|
| **Accumulating** (random walk) | the oscillator's own noise; supply *pushing* on the core | ✅ **covered exactly** — this is `c`, which already reads as **jitter per second** |
| **Additive** (white) | the buffer chain's own thermal noise; the LDO modulating each inverter's **delay**, displacing each edge independently | ❌ **not computed anywhere** |

The second is the number a clock designer actually wants at the last inverter. It does **not**
accumulate and it does **not** appear in `c`.

**THE RIGHT OBJECT IS A SAMPLED / TIME-DOMAIN NOISE ANALYSIS** — noise evaluated at the
**threshold crossings** rather than averaged over the cycle. A commercial RF simulator exposes this as
`noisetype=timedomain`.

⚠ **AND THE SPECTRAL `pnoise` DOES NOT SUBSTITUTE, BY OUR OWN DOCSTRING.** `PAC.pnoise` records
the boundary in terms — *"an oscillator drives a limiter … the same is true when an oscillator
drives a mixer"* — with the test being **whether anything downstream can track the PSD's
variation over the cycle**. A switching inverter samples at its crossing instant, so it **can**.
That is the same fact that makes the sampled analysis the *correct* object rather than a
convenience: **the physical circuit samples, so the analysis must sample.**

✅ **IT UNIFIES WITH SOMETHING ALREADY OPEN, AND THAT IS THE STRONGEST ARGUMENT FOR SCOPING IT.**
The far-out noise floor — where our PM-only spectrum keeps falling at 20 dB/decade while a real
oscillator **flattens** — is set by exactly this additive buffer noise. **Both gaps are one
missing object, not two.** Worth knowing before either is scoped separately.

**REFERENCE, CITED NOT VERIFIED HERE** (relayed by the docs session; nobody in this repo has read
the paper): Demir, Liu & Sangiovanni-Vincentelli, *"Time-Domain Non-Monte Carlo Noise Simulation
for Nonlinear Dynamic Circuits with Arbitrary Excitations"*, **TCAD 15:493 (1996)** — Demir's own,
from *before* the PPV theory. Built on SDE theory; reported to return *"the noise variances and
covariances of circuit variables **as a function of time**"* and *"noise correlations between
circuit variables at **different time points**"* — precisely what a time-averaged PSD throws away
and a sampling stage consumes. Non-Monte-Carlo, and reported to need **no steady state**
(*"any nonlinear dynamic circuit with any kind of excitation which can be simulated by the
transient analysis routine"*), so it would cover the buffer chain and the oscillator in one
analysis. ⚠ If this is ever scoped, that is the starting point rather than a fresh derivation.

⚠⚠ **A GATE WARNING THAT MUST BE HONOURED BEFORE ANY OF IT IS BELIEVED.** Brambilla et al.,
*"Effects of numerical noise floor on the accuracy of time domain noise analysis in circuit
simulators"* (cited, not verified here): time-domain noise analyses implemented by *extending*
linear multistep formulas, or by *introducing sampled versions of noise generators*, are
*"often affected by a **relevant numerical noise floor hiding the effects of noise sources**"*.
**The floor we would be chasing can be manufactured by the method chasing it.** So any gate for
this item needs a **control that measures the floor with the sources OFF** — section D shape 0m
in its sharpest form, since a spurious floor is exactly a clean number with no mechanism
predicting it.

⚠ **AND WE ALREADY CARRY A SECOND WARNING ABOUT THIS CLASS OF ANALYSIS**, in `PAC.pnoise`'s
docstring: Gourary et al. name the symptom as *"the standard time domain noise analysis yields
FLAT PSD CURVES OR CURVES WITH UNEXPECTED SLOPE NEAR THE OSCILLATION FREQUENCY"* — the
`Φ(T) − I` singularity showing through. The same note records that a published removal exists
*in a time-domain form written for shooting*, using the PPV as the null vector, and is not built.
**Two independent warnings that a naive time-domain noise analysis produces plausible-looking
wrong curves.** Neither is a reason not to build it; both are reasons the gate comes first.

**ABSENCE VERIFIED HERE** (2026-09-04, not taken on report): nothing in `pycircuit/` implements
this. The only occurrences of "jitter" outside tests are `elements_hdl.py`'s `idtmod` docstrings
and `shooting.py`'s two uses, both of which are `c` described as jitter per second. There is no
`noisetype`, no sampled-noise generator, and no crossing-time statistic.

⚠⚠ **MEASURED 2026-09-04 ON THE ACTUAL CHAIN — `pnoise` RUNS, RETURNS PLAUSIBLE NUMBERS, AND THE
ADDITIVE MECHANISM IS ABSENT FROM THEM.** Asked whether *unseparated* total output noise is
available for this circuit, since that is the simpler-sounding request. It is available and it is
the wrong number, which is worse than it being unavailable.

Built the topology — van der Pol core, three `tanh` transconductance buffers into RC loads, each
buffer carrying its own noise source, output at the last one. `m = 5`, converged, **`autonomous =
True`** (the buffers do NOT make it driven — the oscillator makes the whole netlist autonomous).

`PAC.pnoise(pss, f0(1+Δf/f0), output)` at the last buffer:

    Δf/f₀      1e-1       1e-2       1e-3       1e-4       1e-5
    S        2.218e-05  1.472e-03  1.250e-01  1.250e+01  1.250e+03

⚠ **AND THE CONTROL — the buffers' own sources ALONE, oscillator source off:**

    Δf/f₀      1e-1       1e-2       1e-3       1e-4       1e-5
    S        1.578e-43  1.516e-41  1.519e-39  1.517e-37  1.517e-35
    ratio         96.0      100.2       99.9      100.0

⚠⚠ **THE BUFFERS' THERMAL NOISE COMES BACK AS A `1/Δf²` PHASE TAIL AND NEVER FLATTENS.** `pnoise`
represents it only through its **orbit-perturbing** effect on the oscillator's phase. The
**additive** mechanism — each inverter's delay modulated independently, displacing edges *without
accumulating* — is not in the number at all. There is no flat floor because the object that
produces one does not exist. (The ~38-decade gap between the two columns is an artifact of the
source magnitudes chosen for the fixture; the **shape** is the finding, and shape is
magnitude-independent.)

⚠ **AND THERE IS NO AUTONOMOUS GUARD.** `pnoise` carries its oscillator refusal in the DOCSTRING
ONLY — *"an oscillator is not this function's problem at all"* — and nothing in the body enforces
it. It computes, on an autonomous circuit, and returns a smooth Lorentzian at every offset. The
near-carrier `1/Δf²` is real physics, which is exactly what makes the result unfalsifiable by
inspection: **a designer reading this number has no signal that the mechanism they care about is
missing.**

⚠ Three things stack up for this one topology, and they are independent: (1) no guard; (2)
`pnoise`'s own docstring names *"an oscillator drives a limiter"* as its incompleteness case, and
an inverter samples at its threshold crossing so it CAN track the PSD's variation, which is the
stated test; (3) the far-out floor never appears. ⚠⚠ **AN EARLIER VERSION OF THIS PARAGRAPH RECOMMENDED ADDING A GUARD. RETRACTED.**
`adjoint_sideband_row` — `pnoise`'s own solve path — **explicitly branches on `pss.autonomous`**
and routes to `_deflated_solve`, so oscillator support there is deliberate and documented, and
the `1/Δf²` above is *correct* phase noise. The docstring line quoted above is about the
**state-covariance** route (`I − M⊗M` singular), which `pnoise` does not take. What stands is
that the answer is **incomplete**, which is a property of the whole PM-only framework (A9) and
not of this function — so the fix is A9, not an annotation here.

**Status: RECORDED, NOT REQUESTED.** No cost estimate has been made, and the item is written down
because the use case is **ordinary rather than exotic** and arrived from a circuit somebody
actually wants to build.

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

### B1. Make `x0_unknown` the default — ⚠⚠ **TRIED UNCONDITIONALLY 2026-09-04 ON B16's EVIDENCE, FAILED THE SUITE, REVERTED**

⚠⚠⚠ **THE UNCONDITIONAL DEFAULT WAS BUILT, MEASURED, AND REVERTED WITHOUT COMMITTING. THE RECORD
OF WHY IS THE DELIVERABLE.** On B16's evidence (the manufactured opener caps `λ₂` at a floor on an
index-1 circuit; `x0_unknown` converges at `h²`) the recommendation was: unconditional
`x0_unknown=True` for `trap`/`euler`. It was implemented, the topology test inverted, and both
motivating fixtures re-measured under the new default — they behaved exactly as predicted
(resonator 19.76939 = explicit `True`; RLC `λ₂` ratios 4.02/4.01/4.00/4.00). **Then the full suite
returned 13 failures**, and they were not thirteen tests pinning the old default:

  * **5 — the known coarse-grid cost, on tests that assert against it**: 19.756 V against 20 V, and
    the `euler < gear < trap` damping order broken because the in-period Euler step makes trap
    *more* damped than gear.
  * **4 — the Newton got better and the tolerance tests lost their lever**: residual 1.8e-15 in 3
    iterations regardless of `reltol`. The old path was a contraction; the new one is a true
    Newton. A premise loss, not a defect.
  * ⚠⚠ **4 — regressions the two fixtures did not predict, and one is disqualifying on its face:
    `trap` recovers only 69.6 % of the analytic amplitude (13.92 V against 20)** in
    `test_backward_euler_damps_the_limit_cycle_and_trapezoidal_does_not`. That is not a 1.2 %
    coarse-grid cost. Also a 4.7e-4 orbit-radius error on an autonomous `trap` solve (bound 1e-4),
    and the `Idtmod` grid-aligned reset losing its discontinuity signature (10⁴× spread where it
    should be constant) — the event interaction B7 already warned the in-period step has.

**VERDICT: the trade is not "a small visible cost for a silent floor". On grids real tests use,
`x0_unknown` costs up to 30 % of the waveform amplitude.** B16's finding stands — the manufactured
opener does cap `λ₂`'s order — but this remedy is the wrong one to ship as a default. Reverted;
the topology-keyed default remains.

⚠ **THE SHAPE, FOR §D: A DEFAULT CHANGED ON THE TWO FIXTURES THAT MOTIVATED IT, NOT ON THE ONES
THAT WOULD FALSIFY IT.** Both re-measurements passed *because they were the same measurements*.
The suite is the falsifier, and it was run — which is the only reason this is a record rather
than a regression.

✅ **THE LIVE ALTERNATIVE, FROM ANDREAS: START THE TRAVERSAL TWO STEPS EARLY.** Prime the multistep
method with its low-order opening steps *before* `t = 0` (run the period two steps longer), so the
period map `[0, T]` is taken entirely with full-order steps and the opener's order-drop is
excluded from it — **without moving an Euler step inside the period**, which is what costs the
amplitude. It targets B16's floor by a different mechanism from `x0_unknown` and may keep the
coarse-grid waveform. ⚠ It does **not** obviously address §0k's index-2 *inconsistency* (an
algebraic variable seeded at a forbidden value is carried forever by trap wherever the period
starts), so the two remedies may be complementary rather than alternatives.

⚠⚠⚠ **BUILT AND MEASURED 2026-09-04 — IT NEVER CONVERGES FOR `k ≥ 2`. REVERTED.** Implemented as
an opt-in `preroll=k` on the plain path: `k` order-dropped Euler priming steps at
`t = −(k−1)h … 0`, then the period, with the sensitivity seeded at `x₀` exactly as today.

  * ✅ **`k = 1` reproduces the default TO EVERY DIGIT on both fixtures** — resonator peak
    `20.0127317714`, B16 `λ₂` errors `4.568e-05 / 8.371e-06 / 1.193e-05 / 7.894e-06 / 4.428e-06` —
    which is the built-in gate that says the plumbing is right.
  * ❌ **`k = 2` and `k = 3` never converge**: resonator at 100, 101 and 200 points (so not the
    `(−1)ᴺ` parity obstruction), RLC at every `N` tried. The non-converged waveforms sit at
    20.44 / 21.04 V — wrong in a consistent direction, not noise.

⚠ **THE MECHANISM IS NOT ESTABLISHED, and a wrong one was nearly recorded.** I first attributed
the failure to the last priming step being a *trap* step and seeding `Pq` with trap's
coefficients. Then a "second variant" with all-Euler priming returned **bit-identical** numbers —
because `solve_timestep`'s `iq_last` already defaults to `None`, so every priming step was an
Euler restart in *both* runs and the "variant" was the same code. The identical number is what
caught it. So the priming steps were Euler all along, the trap-seed explanation is **wrong**, and
what stands is only the record's own statement: the plain path's fixed-point iteration depends on
`x₀` being **exactly one** Euler step from the unknown, and the reason is not yet understood here.

**VERDICT:** the "start earlier" family is closed as a remedy for B16's floor. The floor remains
unremedied except by `x0_unknown`, which costs the waveform (see above). ✅ **What would actually
move this is understanding why `k = 1` converges as a contraction and `k = 2` does not** — the
frame error the record names, with the true `dF/dx_in` singular. That is an analysis item, not a
build.

⚠ **NARROWED 2026-09-04 (docs session, measured):** the obvious mechanism — `dx(0)/dx_in = A_eu^k`
becoming ill-conditioned with `k` — is **refuted**. On the series RLC, cond over the non-null
directions is 1.60 / 2.04 / 2.67 at `k` = 1/2/3 (100 pts) and 1.43 / 1.47 / 1.52 (400 pts), with
**rank 2 at every `k`**. The algebraic directions are annihilated by ONE Euler step and no more
annihilated by two, so the singularity of `dF/dx_in` is *identical* for `k = 1` and `k = 2`.
Whatever separates them is **not** rank and **not** the conditioning of the pre-image map.

⚠ **AND THE CORPUS DOES NOT REACH IT** — looked for, not assumed. Estévez Schwarz & Tischendorf
2000 characterise which `x(0)` are *consistent* and defer initialisation to three Humboldt
technical reports we do not hold. The closest title in existence to this question is
**März & Rodríguez Santiesteban, "Analyzing the stability behaviour of DAE solutions and their
approximations", TR 99-2, Humboldt-Universität Berlin, 1999** — an acquisition, not a search, if
this becomes worth chasing.

⚠⚠ **RETRACTED BY THE DOCS SESSION THE SAME EVENING: THE CORPUS DOES REACH IT.** We hold
**R. März, "On linear differential-algebraic equations and linearizations", Applied Numerical
Mathematics 18 (1995) 267–292** — filed under a mangled DOI with no author or title, which is why
a filename-keyed audit missed it. Its abstract: linearizations of **nonlinear index-2** systems,
*"the local convergence of the **Newton–Kantorovich method (quasilinearization)** result
immediately … this applies also to fully implicit index-1 systems whose leading nullspace is
allowed to vary with all its arguments."* Newton–Kantorovich on a boundary value problem is the
framework shooting lives in. ⚠ **AND IT NAMES A NORM:** convergence *"with any initial guess
`x₀` being close enough to `x_*` in **C¹**"* — the guess's **derivative** must be close, not only
its value. ✅ **A HYPOTHESIS WITH ITS TEST ATTACHED, marked as such:** if the `C¹` ball is what
separates `k = 1` from `k ≥ 2`, then the discriminating measurement is not `x(0)` but
`‖ẋ_manufactured(0) − ẋ_orbit(0)‖` as a function of `k` — it should degrade where convergence
does. If it is flat while convergence fails, the `C¹` story is wrong. It also fits the shape the
conditioning measurement could not: settling 2–5 % off in a consistent direction is what falling
outside a local-convergence ball looks like, which is a failure of the *initial guess*, not of
the Jacobian. **Unrun.** ⚠ The general lesson: ~30 files in the library carry no author or title
(DOIs, arXiv ids, `selting1997.pdf`, `2763.pdf`); an index keyed on filenames is blind to exactly
those, and one of them answered a question declared unanswerable. **The index is not the
thing.** Two more surfaced from the same blind spot, unassessed: Selting & Zheng 1997 (stability
of self-excited oscillating circuits, J. Comp. Appl. Math. 82) and a 2020 Russian paper on
adaptive stepping for oscillatory circuits (DOI 10.31114/2078-7707-2020-3-28-34).

⚠⚠ **AN INDEX-2 GAP IN THE SENSITIVITY SYSTEM — HYPOTHESIS WITH ITS TEST, UNTESTED (docs session,
from Bereza, *Identification of Non-Linear DAEs*, KTH licentiate 2024, filed as `FULLTEXT01.pdf`).**
Propositions 3.1/3.2: for an **index-1** DAE of two common forms (circuit MNA is the second), the
forward sensitivity system `F_ẋ ẋ_θ + F_x x_θ + F_θ = 0` **is also index 1**, and concatenations
stay index 1. That is the licence — assumed everywhere in this tree, justified nowhere until now —
to propagate monodromy columns with the same integrator, order and step control as the state.
**The guarantee covers index 1 only.** On the index-2 configurations (L-I cutset, C-V loop) *nothing
guarantees the sensitivity system inherits the nominal index*. ✅ **Prediction:** if the sensitivity
index exceeds the nominal one there, the **monodromy rows lose convergence order FASTER than the
state rows** on the same fixture. **Refutation:** they degrade together, or the monodromy is no
worse — then the index-2 defect lives in the state solve alone. Every measurement so far (§0k, B16)
covers the **state** only. If real, it locates the 2× in sensitivity propagation — a different fix
from re-seeding. ⚠ The fixture needs *both* a state reference and a monodromy reference on an
index-2 circuit; not built.

✅✅ **THE MONODROMY REFERENCE EXISTS AND NEEDS NO INTEGRATION (docs session, 2026-09-04 late).**
For a linear constant-coefficient circuit the variational system is the pencil `(C, G)`:
`eig(M) = {exp(T·μᵢ): μᵢ a finite generalised eigenvalue of −Gv = μCv} ∪ {exactly 0 on the
nilpotent block}` — by the Weierstrass–Kronecker form, **valid at index 2 identically to index 1**
(a higher index only enlarges the nilpotent block). Computed with `scipy.linalg.ordqz(−G, C)`,
which never forms `C⁻¹`. A genuinely external instrument.

⚠⚠ **AND IT SAYS THE L-I CUTSET ALONE HAS NO MONODROMY AT ALL** — `n = 2`, `rank(C) = 1`, zero
finite eigenvalues, `exact |eig(M)| = [0, 0]`: every mode is algebraic. The state has dynamics
because the source drives it, but the homogeneous variational system has no non-trivial solution.
All three methods return `[0, 0]` and "agree" trivially. **The fixture had nothing to reference**,
which is why one could not be built for it. The fixture that has both — index 2 *plus* a
differential mode — is `li_plus_rc` (`I` in series with `L`; an independent `V → R → C` with
`τ = T` so `μT = −1`): `exact |eig(M)| = [e⁻¹, 0, 0, 0, 0]`, and against it all three methods hit
textbook order — gear 2.02/2.00, trap 2.01/2.00, euler 1.00/1.00. ✅ **That is the REFUTATION branch
of the Bereza prediction on this fixture** — with the stated caveat that here the index-2
constraint and the differential mode sit in separate branches; a fixture where the same branch
carries both is the harder test and is not built.

⚠ **AN UNEXPLAINED OBSERVATION, flagged so nobody builds on it:** on the plain index-1 RC control
(`V → R → C`, same values) **trapezoidal converges at order 1**, agreeing with euler to four digits
(0.36974 / 0.36973 against exact 0.36788) — while the SAME branch inside `li_plus_rc` gives trap a
clean order 2. Not a top-k truncation artefact (3×3, clean spectrum). The trap/euler agreement
points at the opening step dominating, consistent with B16 — **but untested, and recorded as an
open observation about the tree, not a result.**

⚠ **Two instrument errors in this attempt, both mine, both from the record's own list:** a
`pkill -f "pre_b[.]py"` that killed its own launcher because the *other* lines of the same command
mentioned `pre_b.py` (the self-matching trap in this campaign's notes, third time); and the
non-variant above. Neither reached a conclusion; both cost a run.


Shipped as an option. The evidence says it wins exactly there and loses on uniform grids:

| | default | `x0_unknown` |
|---|---|---|
| Q=20 @ 100 pts (analytic 20 V) | 20.01273 | 19.76939 |
| Q=20 @ 200 pts | 20.02208 | 19.96123 |
| van der Pol, LTE grid | −73.8 ppm | **−47.3 ppm** |

**Gate:** a rule that picks correctly without the caller knowing. "Non-uniform" is not
quite it — the gain came from the formulation making the opening-step *subdivision*
unnecessary, so the real predictor is whether the grid opens coarse.

---

⚠⚠ **THE GATE WAS RUN 2026-09-04. "OPENS COARSE" IS FALSIFIED AS A PREDICTOR, AND THE
RECOMMENDATION IS TO LEAVE THE DEFAULT ALONE.**

**Driven**, Q=20 resonator, grid graded geometrically so it opens coarse while staying inside
the `1 + √2` zero-stability bound (absolute error against the analytic 20 V):

    h0/mean    1.00    1.39    2.01    2.55    3.49    4.59
    default   0.0127  0.0148  0.0075  0.0073  0.0417  0.1247
    x0_unk    0.2306  0.4376  0.9156  1.4459  2.5582  4.0442

⚠ `x0_unknown` loses at **every** grading and loses **monotonically more** as the grid opens
coarser — the OPPOSITE of the predicted direction.

**Autonomous**, van der Pol on its own `lte_grid`, period error in ppm against a fine reference:

    μ        1        4       10
    default  1.4    101.0      0.5
    x0_unk  21.9    109.5    108.8

⚠ It loses there too, at every `μ` tried. (`gear` returns the documented "nothing to change"
refusal, since its solved-history path already solves for `x₀`.)

⚠⚠ **AND THE ONE RECORDED WIN DOES NOT SUPPORT THE CHANGE, BECAUSE IT IS NOT ABOUT THIS FLAG.**
Re-running `benchmarks/pss_lte_grid.py` reproduces **−47.3 ppm**, but that is the **LTE GRID**
beating uniform — where 1105 uniform steps *do not converge at all* — measured through the
benchmark's own hand-rolled Newton. The `−73.8` it is contrasted against is the same grid with
the opening-step SUBDIVISION, which is a different knob. So the table above entangles a grid
change and a formulation change, and the formulation half is unsupported.

✅ **RECOMMENDATION: DO NOT make `x0_unknown` a grid-driven default.** Keep the conditional
**topology**-driven default (index-2 → on), which rests on the L-I cutset measurement in §0k and
is the only regime where the flag is measured to win. ⚠ A caller with a specific reason can still
pass it explicitly — that is what an option is for.

⚠ **ONE HARNESS TRAP WORTH THE RECORD.** The first version of the driven sweep put the entire
opening ratio into step one, giving `h₀/h₁` of 10 and 100 — far outside the `1 + √2` bound that
this file's own `test_a_grid_that_outruns_zero_stability_says_so` documents. It was therefore
measuring the integrator's silent **demotion to Euler**, not the formulation, and one grid
returned `nan`. Geometric grading spreads the same overall opening across every step and stays
inside the bound by construction. **A grid is an instrument, and an invalid grid measures the
integrator instead of the thing under test.**

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
| C2 | **Poincaré / orthogonality phase row** | Tested and rejected. 2–6x conditioning edge that shrinks with refinement and never decides an outcome. The `argmax` rule is canonical (A&T Step 3) and compares units *on purpose* — normalising picks a row **704x worse aligned**. ⚠ **SCOPE: this is one phase ROW against another. It does NOT cover removing the phase condition entirely** — that is LSOAC, which was never considered here and is now **B10**. |
| C3 | **Per-iteration phase re-selection** | Built and reverted: it **regressed the working case** (on-orbit seed went converged → not). Structural — pinning the iterate's own value makes the residual identically zero, so the row constrains only the step. A&T avoid it by having no phase equation at all (see B3). |
| C4 | **Index-2 detect-and-refuse** | Would reject working circuits: `index > 1` is **not predictive** (all three methods converge on an LI-cutset), and gear is not the workaround (fails on 2 of 4). Failures are loud, not silent. ⚠ **STILL CLOSED after B11:** the index turns out to be *decidable* from topology, which changes nothing here — deciding it was never the obstacle, and it still does not predict convergence. B11 proposes a DIAGNOSTIC that names the offending loop or cutset, not a refusal. |
| C5 | **Saltation correction for switching** | **Falsified.** The monodromy–FD gap falls at exactly 2.00x per doubling = O(h) discretisation, not the O(1) a missing correction leaves. PSS's monodromy is the derivative of the *discrete* map, which has no undefined instant. |
| C6 | **Reusing the inner Newton's factorisation** | Does not exist to reuse — `StandardNewton` factors and discards. A retained one would be `J(x_k)`, not the converged `J(x_{k+1})`: already measured at median 5.5e-9 and rejected as an approximation. |
| C7 | **Gourary adaptive preconditioning for PAC** | Harmonic balance. His own introduction says the shooting case was already solved by Telichevesky et al., because `A' = I` there — which is our case. Solves a problem we do not have. |
| C8 | **Deflation for the trivial root** | "The condition number of the system's Jacobian matrix grows beyond any bound and convergence stalls." Not the shortcut. |
| C9 | **A GMRES preconditioner** | The per-step factored solves *are* the preconditioner, in both DAC'95 and DAC'96. ⚠ Do not carry this forward as "this operator needs no preconditioning" — that is false, and the clustered spectrum is a *consequence* of the implicit preconditioning. |

---

## D. How these items keep failing — the shapes worth checking for

Sixteen claims were overturned across this campaign. Four shapes account for most:

0c. **A fixture whose blind spots were never recorded next to its coverage.** "Verified to
   1e-15" says nothing about which errors the check is *capable* of seeing. Van der Pol's
   `C = diag(1, −1)`, so along the capacitor node `v·δ` and `vᵀCδ` are numerically
   **identical** — the 7% transcription error is invisible to any single-direction probe
   there, and two gates in this campaign probed exactly there. The shipped physical gate
   survives only because it uses four *random* directions. ⚠ **A unit reactance makes `C` the
   identity up to a sign; run the fixture at `c = 2` and the same blind direction separates
   the two formulations by exactly the capacitance.** Now asserted on both fixtures, so the
   blindness is a measured property rather than a hypothesis. (Independently found by the
   docs session on Srivastava's ring oscillator, where `τ = 1` makes `C = I` outright and its
   PPV oracle passes a missing-`C` implementation to 2e-15. Same shape, two circuits.)

   ⚠ **AND THE BLINDNESS IS COMPUTABLE, SO THE PROBE CAN BE DESIGNED RATHER THAN DRAWN.**
   `v·δ = vᵀCδ` exactly when `vᵀ(C − I)δ = 0`, so the blind set is a **hyperplane with a known
   normal** and `(C − I)ᵀv` is the direction of maximum sensitivity — one matrix-vector product
   from quantities already in hand. Van der Pol's `C = diag(1,−1)` gives `(0, −2v₁)`, and `e₀`
   is orthogonal to it *exactly*, which is why the capacitor ratio is 1.0000 rather than merely
   small. Along the designed probe the two hypotheses predict **opposite signs**, so no
   tolerance is needed to separate them; the circuit picks one (measured −5.4215e-01 against
   `v·δ` = −5.4131e-01 and `vᵀCδ` = +5.4131e-01, converging 0.9984 → 0.9992 at O(h)).

   ⚠ **The structural rule generalises past this circuit: in an MNA-shaped `C` the ALGEBRAIC
   rows are the best probes.** `C` is *zero* there, so `C − I = −I` and the discrepancy is
   maximal — the rows a DAE solver already treats specially discriminate best, and the
   capacitor nodes everyone reaches for first are the blind ones.

   ⚠ **A designed probe's optimality is about ONE error**, so the random-direction gate stays:
   it is not aimed at a hypothesis. Recording what a check can see applies to the check built
   from this rule too.

0d. **A blind spot with NO good probe, because the discriminating quantity is identically
   zero.** The dropped-`C` error has a blind *hyperplane* — bad directions alongside good
   ones. The **transposed**-`C` error is separated by `C − Cᵀ`, so on a symmetric `C` no
   direction distinguishes the implementations at all. ⚠ **MEASURED, and the code it certifies
   is correct only by inspection:** `_monodromy_matvec_transposed` does `cs0[j].T @ t`, and its
   recorded gate ("agreement 1.8e-15" against the dense `Mᵀ`) is passed *identically* by an
   implementation with no transpose in it.

   | fixture | shipped | **with the `.T` dropped** |
   |---|---|---|
   | symmetric `C` (vdP, RC — every fixture we had) | 3.994e-16 | **3.994e-16 — blind** |
   | non-symmetric `C` | 3.994e-16 | **4.667e-01 — caught** |

   ⚠ **AND THE CLASS IS REAL IN THIS TREE.** `compact.PspMosLongChannel`'s `C` is
   non-symmetric — `Cgd = −4.31 fF` against `Cdg = −0.13 fF`, a factor of **33**. That is
   Ward-Dutton charge partition, the model being *right*. So the transpose matters the moment a
   PSS runs on a MOS circuit, and nothing in the suite would have changed state to say so.
   Now covered by `_TransCap`, a fixture that exists for no other reason.

   ⚠ **QUOTED AGAINST `Cox`, AND THE FIRST ATTEMPT USED THE WRONG DENOMINATOR.** McAndrew's
   figure is a nonreciprocity over `Cox`; the 0.44 first written here was a ratio to `max|C|`,
   and 33 is a spread between two entries — three different normalisations, one of which was
   compared to the other two. Measured properly: `|C_ij − C_ji| = 4.97 fF` against
   `Cox = 15.35 fF` → **0.324 = 32× McAndrew's 0.01**, not 44×. `Cox` anchored *outside* the
   `C()` code by the model's own geometry (`ε_ox W L / tox` = 15.70 fF, **2.2%** away).

   ⚠⚠ **AND FIXING THE DENOMINATOR DID NOT FIX THE COMPARISON — THE CONDITIONS WERE ALSO
   MISMATCHED, WHICH IS A DIFFERENT DEFECT WITH A DIFFERENT REPAIR.** McAndrew's bound is
   stated at **`VDS = 0`**, under his eq. (2), derived there. The 4.97 fF is a worst case over
   a box *containing saturation*, where Ward-Dutton partition makes `Cgd`/`Cdg` asymmetric **by
   design**. At his condition:

   | | nonreciprocity / `Cox` | × his 0.01 |
   |---|---|---|
   | `Vds = 0`, worst over `Vg` | 0.844% | **0.84 — under the bound** |
   | `Vds = 1.2`, `Vg = 1.2` | 31.7% | 31.7 |

   monotone in `Vds`: 0.47, 1.63, 3.77, 9.47, 20.9, 30.6, 31.7 at `Vds` = 0 … 1.2.

   ⚠ **So his 1% is VERIFIED, not contradicted, and the "floor for the ideal case, not a
   typical value" written here one commit ago was wrong.** A real compact model gated to 1.3e-6
   against a compiled PSP103 satisfies his bound at his stated condition to within a factor of
   **1.2**. The 32× is a statement about **saturation** — what this gate cares about, and not a
   claim about the paper.

   ⚠ **Useful corollary: the nonreciprocity is essentially CREATED by `Vds` (68× growth), so a
   DC-biased fixture is a weak transpose gate and a swinging one is a strong one.**

   ⚠ *The normalisation is where every error in this campaign has lived* — and this one had a
   **conditions** defect underneath it that survived the normalisation fix entirely. See shape
   0e.

   ⚠ **The fixture's sensitivity is LINEAR in the asymmetry, not thresholded** — 0.6667×asym
   across four decades, so at McAndrew's own 1% the error is still 1e-2 against a 4e-16 floor.
   That is a different structure from a probe-direction failure (a *direction* misaligning, with
   a distribution over draws) and is why "weak asymmetry defeats it" does not transfer: here the
   *quantity* shrinks, deterministically. The discriminator dies only at exactly zero.

   ⚠ **The earlier "max|C − Cᵀ| = 0 across the element library" was scoped to the DISCRETE
   library and used as though it covered the tree.** The compact models were never in that
   sweep. Shape 0 again, on my own measurement.

0e. **A comparison whose CONDITIONS do not match, which survives fixing the normalisation.**
   Distinct repair: a normalisation defect is fixed by finding the comparable *denominator*; a
   conditions defect is still there afterwards, because you can make a number comparable in its
   denominator and still have measured it in a regime the bound was never asserted for. ⚠ The
   check is not "what is this divided by" but **"under what conditions was the thing I am
   comparing against asserted"**. Instance: quoting a worst-over-bias-box nonreciprocity against
   McAndrew's at-`VDS = 0` bound — and doing it *in the message correcting the denominator*, so
   the correction itself carried an over-claim.

0f. **"Bounded" and "already bounded" are different claims** — an asymptotic property used as
   if it held from `t = 0`. ⚠ **Instance, mine, 2026-09-04:** designing a high-`Q` Monte Carlo
   to *skip* settling by fitting the **slope** of `Var(θ)` vs `n`, reasoning from
   `oscillator_covariance`'s own split that the transverse part is **bounded** (`K_orb`) while
   only the phase random-walks — so the un-decayed transverse would fall into the *intercept*.
   The argument is sound and the conclusion is false: `K_orb` is reached only after
   `~1/(1−λ₂²) ≈ 30` periods, and **the approach to it is monotone growth**, so during the ramp
   the transverse contribution is not a constant offset but **a second slope**. Measured:
   `Var/n` climbing 3.9e-11 → 3.1e-10 → 3.3e-9 over 40 periods, linear-fit residual **16% of the
   range**. ⚠ **The settling is what establishes the split; it cannot be used to avoid the
   settling.** The tell is free and was ignored: a slope fitted through a curve has a residual,
   and the residual was reported in the same table as the answer.

0b. **A measurement that shares the assumption under test.** `diffusion_constant`'s Monte
   Carlo used the same one-sided-as-two-sided injection as the code, so it agreed to 0.9965
   while both were 2× wrong; only `kT/C`, external to both, could see it. **Ask what the
   measurement assumes before trusting what it confirms.**

0j. **Comparing two instruments and blaming the subject.** B7c's gates 1 and 4 took BOTH period
   columns by finite difference "off the same `_traverse`, so the comparison is between the two
   conventions and nothing else" — true, and beside the point, because one of the two finite
   differences was itself `O(h)` wrong. That produced a confident, twice-committed claim that a
   shipped Jacobian column was `O(h)` (later 4.2%) wrong; **it is not, it converges at `O(h²)`**.
   ⚠ **The analytic column the code already computed was one keyword away (`want_dT=True`) the
   whole time.** When a measurement says a shipped implementation is wrong, VALIDATE THE
   INSTRUMENT AGAINST THAT IMPLEMENTATION before believing it — the implementation is the cheaper
   thing to check, and it is the one with a track record.

0i. **A fixture in which the WRONG QUANTITY IS NUMERICALLY INDISTINGUISHABLE FROM THE RIGHT ONE.**
   Four instances in two days, and every one passed every gate that existed at the time: unit
   component values (`C = 1 F`, `L = 1 H`) hid a **dimensional** error — a missing `C` in the `CY`
   contraction was exactly 1, invisible to `kT/C`, to the Monte Carlo, to §0c's pointwise gate and
   to §0e's quadrature study alike; unit amplitude hides a **scale** error; a time vector read as a
   node voltage hid a **units** error (the docs session got `A = 6.66`, the period, and a plausible
   ratio of 7.10); and a single series resistor hid a **sign**, because `|∫v₀|` and `|r∫v_branch|`
   agreed to 1.5e-4 there. **Sweep the units, not just the regime.** ⚠⚠ **AND THE RULE IS
   CONSTRUCTIVE, NOT A POST-HOC CHECK.** "Could the answer have come out otherwise?" is applied
   after the fact and is easy to answer wrongly; **BUILD THE FIXTURE SO IT CAN EXPRESS THE WRONG
   ANSWER** is a design applied before. §0j's `K_orb` case is the worked example: sweeping `C` at
   FIXED `τ/h` killed two competing explanations at once, because a missing term must scale with
   something and a resolution limit cannot — a discrimination no amount of agreement between
   references could have produced. (Sharpened with the docs session, 2026-09-04; shape 0h is the
   sign case of it.)

0h. **A fixture where the candidate answers are numerically degenerate.** The single
   series-resistor tank makes `|∫v₀|` and `|r∫v_branch|` agree to 1.5e-4, so a prediction that
   matched in magnitude matched for EITHER SIGN — and the match was recorded as a confirmation
   before a divider fixture, whose two algebraic nodes must stand in a topology-fixed ratio of
   10.000000, showed the sign was inverted. **Before trusting a sign, check the fixture can
   express the wrong one.**

0g. **A control that is a whole-period integral, guarding a pointwise quantity.** The DC
   frequency-shift gate agrees to `1e-4` while a *localised* functional of the same `v₀` is 32%
   out. An `O(h)` error that integrates to zero over the period — a time shift is exactly one —
   is **invisible** to `∫ v₀ dt` and first order in every windowed functional. **A control is
   only evidence for the functional it actually is;** passing an integral gate says nothing
   about the integrand.

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

0m. ⚠⚠ **A DEGENERATE FIXTURE, WHICH AGREES WITH ANYTHING AND REPORTS OK.** B8's driven gate
   used `R = 1k` with `L = 1mH`, `C = 1nF` — `Q ~ 1e-3`, so the monodromy decayed to numerical
   zero over the period. The gate compared the zero matrix against the zero matrix, divided by
   `1e-30`, and printed **relative agreement of 1.9e-46, verdict OK**. A clean measurement of
   nothing. ⚠ **THE TELL WAS THAT THE NUMBER WAS TOO GOOD** — but "suspect a clean residual" is
   the wrong rule, because `1e-15` on a float64 algebraic identity is clean *and correct*, and
   that rule would have thrown away the same gate's real result. The right form is a
   **pre-commitment**: *name the number the arithmetic predicts before accepting the number it
   produced.* Nothing predicts `1e-46`. Every gate that divides by a scale must now assert that
   scale is non-degenerate.

0n. ⚠ **A FAILING ASSERTION WIDENED INSTEAD OF EXPLAINED.** B8's autonomous test demanded the
   transposed map's dominant multiplier be `1.0` to `5e-3` and failed at `0.992007`. The
   instinct is to loosen the bound. Measured instead: **`M` itself carries that deficit** — 120
   points on a `Q = 8` limit cycle is a coarse grid — and it converges at `O(h²)` (7.99e-3 →
   1.87e-3, ratio 4.3, gear exact to 0). The assertion is now against the *forward* map's
   spectrum plus that rate, which is a stronger test than the one intended. **When an assertion
   fails, find the true value before widening the bound**; a loosened tolerance would have
   hidden a second-order effect behind a green light. This is 0m pointed at the failing
   direction rather than the passing one.

0o. ⚠ **CREDIT MISASSIGNED *AWAY* FROM ONESELF, WHICH DESTROYS THE RECORD JUST AS EFFECTIVELY.**
   Told a peer "I have not reproduced your parity numbers and will not represent them as mine"
   — while §0k of this document already carried the full table, measured here, hours earlier,
   including the correction to my own earlier wrong claim. Guarding hard against claiming
   another session's work produced the mirror error. **Check the record before disclaiming, not
   only before claiming.**

0p. ⚠⚠ **A PAPER READ FOR ONE CLAIM AND FILED UNDER THAT CLAIM.** Traversa & Bonani 2011 was
   in the library, had been opened, and was cited for the `Φ = U·D·V·C` factorisation — a claim
   it does not own (that is Demir 2000 eq (37)). Its **actual subject**, the orbital half of the
   oscillator spectrum, was then reported as an open gap with "no published answer" hours later,
   by two sessions independently. **A library loses a result it already contains this way, and
   the loss is invisible** — the paper is present, indexed and cited, so every search for it
   succeeds while every search for its subject fails. ⚠ The tell available in advance: **a
   citation whose title does not match the claim it is attached to.** Check the title against
   the claim, not just the claim against the source.

0q. ⚠⚠ **A RIGHT FORMULA ON THE WRONG OBJECT — AND CHECKING THE ARITHMETIC CATCHES NONE OF IT.**
   Across a four-session exchange on one question (does `ppv()`'s `Q` equal a designer's `Q`),
   **every single reversal was a SCOPE error, not an arithmetic one.** `π` was exact — for a
   passive resonator, applied to an oscillator. Wang & Roychowdhury's containment was real — for
   linear systems, quoted about oscillators. A van der Pol fixture measured its own `μ` back. The
   describing-function slope was right — for a hard limiter, generalised to a cubic. **Four
   times, right formula, wrong object.** ⚠ Recomputing any of them would have confirmed them.
   The check that works is naming the object a derivation is *about* and asking whether the thing
   in front of you is that object — and the cheap version is to **sweep rather than spot-check**,
   because a wrong constant can pass through the right value at one operating point (see A9's
   coincidence at `a = 1.5G`). ⚠ **Corollary for advice: "the algebra says it will reproduce X,
   don't bother measuring" is this shape at its most expensive**, because it spends someone
   else's time and removes the check that would have caught it.

0r. ⚠⚠⚠ **THE `npts`-VERSUS-STEPS OFF-BY-ONE — THREE TIMES IN ONE DAY, IN THREE DIFFERENT
   PLACES.** `timestep = per/npts` yields `npts` POINTS and `npts − 1` STEPS, spaced
   `per/(npts−1)`. It appeared as the parity confusion in §0k, as the "converged=False is not
   silent" error, and in B16 as a REFERENCE built with `h = T/npts` raised to the power
   `npts − 1`. ⚠ **The B16 instance is the instructive one because the error was in the
   INSTRUMENT, not the subject**: an `O(h)` term one step of period wide, larger than the effect
   under test, which flattened the convergence being measured and manufactured two plausible
   relationships that were then reported as findings. ✅ **The fix is structural, not vigilance:
   take the step count and the step size FROM THE RETURNED WAVEFORM, never from the requested
   `npts`**, so the exponent and the step cannot disagree — and assert `h·N == T` before using
   either.

0s. ⚠⚠ **A RELAYED CLAIM HAS NO AUTHOR PRESENT TO FEEL UNEASY ABOUT IT.** Three times today a
   sentence arrived from another session, correct-sounding and sourced, and was wrong in scope:
   "a few more Floquet pairs, not a sweep" (the requirement is ALL eigenvectors — IET CDS 2011);
   "PAC currently refuses the Gourary transformation" (it implements the same identity by
   bordering — B12); and "the closed-form variance has no justification anywhere in the high-Q
   programme" (true of the book's (6.72) derivation, false of the tree — §0 row (c)). The first
   reached shipped code as `nmodes=2`. ⚠ **Relayed claims and derived claims fail the same way
   (§D 0q) but are FOUND differently.** A derived claim carries its author's unease about the
   step they were unsure of; a relayed one arrives finished, with the unease stripped off in
   transit, and repetition by a second session reads as corroboration when it is only
   propagation. The two that were caught were caught by the same act: **refusing to write the
   relayed sentence into the tree until the source page was read** — and in the (6.72) case the
   page put the condition somewhere neither the relay nor the reconstruction had considered.
   ✅ The rule: a claim that arrives from another session is a POINTER to a source, not a
   finding; it is recorded with "cited, not verified here" until someone reads the page, and it
   is not written into code or a docstring caveat before that. ⚠ **And the read-vs-cited audit
   must report its NULLS**: the docs session's pass found three load-bearing relay-only citations
   AND one (Ngoya 1995) that came back verbatim faithful to a row written without reading him. An
   audit that reports only its hits is not an audit; the clean result is what licenses trusting
   the rows it did not flag.
   is not written into code or a docstring caveat before that.

0t. ⚠⚠ **A PERIODICITY GATE IS A NECESSARY CONDITION ONLY — it is blind to ANY error that
   returns to its starting value.** Two failures today, from different mechanisms, passed
   periodicity while being wrong everywhere in between. §0k: trapezoidal returned **exactly 2×**
   on an L-I cutset with `converged=True` and a periodicity residual of **1e-13**, because the
   ripple `v_n = V[cos(ωt_n) − (−1)ⁿ]` closes on itself over an even step count. A9 step 3:
   `floquet_modes` carried a state-block scale of `q(0)ᵀp(0) = 1.324`, **constant around the
   cycle**, invariant under the flow, and `p(T) = p(0)` held to 3e-15 — periodicity is
   scale-free. ✅ `p(T) = p(0)` constrains the endpoints of a trajectory with many degrees of
   freedom; it can never be sufficient. **Every periodicity assertion needs a companion that
   fixes the scale or the interior** — `qᵀp = 1` around the cycle in one case, a closed-form
   amplitude in the other. Both sessions reached for the same repair independently after being
   burned, which is some evidence it is the right one. (Named by the docs session, from the
   "scale-free" observation.)

0u. ✅ **BUILD THE THIRD ROUTE SO THAT AGREEMENT BETWEEN THE FIRST TWO CANNOT BE
   SELF-CONFIRMING.** The shared-instrument trap (§D 0r, the B16 reference) used constructively:
   in A9 step 3 the modal sum and a definition-route integral agreed to 3.5e-4 while BOTH
   disagreed with the Lyapunov reference by 1.75×. Two formulas agreeing is not evidence when
   they share an input — but two formulas agreeing while a third disagrees ISOLATES the fault to
   the shared input rather than to either formula, which is what found the scale defect. The
   third route was built for that purpose, not for extra confidence. **When two routes agree,
   ask what they share; when a third disagrees, the shared thing is the suspect.**

0v. ⚠⚠ **REFUTATION CONDITIONS MUST PARTITION THE OUTCOMES.** A9's Q-sweep discriminator offered
   *grows with Q → physics, flat → numerical* for a **ratio** whose numerator and denominator
   scale differently. The measurement did neither — it **fell** as 1/Q — and read literally the
   test did not apply. The right conclusion came from a different route (the projection), not
   from the test. Naming falsifiers in advance (§D 0m) is worth something only if the named
   outcomes exhaust the space; a two-outcome test on a three-outcome quantity licenses nothing,
   and a "third outcome" is the test's failure, not a bonus. (Named by the docs session, of its
   own test.)

0w. ⚠⚠ **A NULL RECORDED AS A RESULT.** The measurement does not run — a non-convergence, a crash, an
empty result set — and its SILENCE is written down as a value. The C-V-loop cell of the index-2 pricing
(§0l, 2026-09-05) read `NoConvergenceError` at every grid, and the standing advice was "record rather
than retry"; the cause was the seed (a bias node held at 0 V against a 1 V source), and correctly seeded
the cell gave the answer that closed the question. Following the advice literally would have filed a
false negative in the direction that closes an open question. The instrument declining to report is not
the quantity being absent, and the two look identical in a table. Repair: name the null as its own
branch when the outcomes are named, and attach to it what must be checked before it counts — seeds,
operating point, whether the code path is reachable on that fixture at all. "Record rather than retry"
governs MEASURED VALUES; it must never be applied to FAILURES TO MEASURE. Sits beside 0b, not inside it:
0b is an instrument that runs and lies (repair: an external reference), 0w one that does not run and
whose silence is read as a value (repair: a setup check). Named by the review session from this
session's report.

**And one about measurement itself:** this machine runs more than one agent. Check
`ps -eo pid,pcpu,args --sort=-pcpu` and `uptime` before trusting any wall-clock ratio — a
concurrent run moved readings 25-30% on the *same* configuration.

---

## TR-BDF2 (a two-stage DIRK alongside the LMM tree), 2026-09-05

The integrator tree was three linear multistep methods (Euler, trapezoidal, Gear-2). TR-BDF2
is the first STAGE method: a trapezoid sub-step to an internal stage at `t_n + gamma h`
(`gamma = 2 - sqrt(2)`), then a BDF2-shaped stage to `x_{n+1}`, two Newton solves sharing one
factorisation (`a22 == a33`). L-stable, order 2, self-starting.

T1. **Fixed-step transient first** (committed 1318db9). `TRBDF2Integrator` is NOT an LMM and
   the ABC's companion-recursion contract does not fit it: `compute_derivatives`,
   `companion_coefficients`, `companion_dT`, `compute_lte` all raise with reasons rather than
   return a plausible-but-wrong LMM companion. The step lives in `Transient._solve_timestep_trbdf2`,
   dispatched by an `isinstance` check; `_solve` refuses a non-fixed grid (the embedded 2(3)
   estimator is item T3, not built). Verified order 2 on the analytic RC step.

T2. **The shooting monodromy is its own map, `m x m`, no opener, no pair** (this increment).
   ⚠ THE ONE-STEP PLAIN PATH IS THE WRONG TEMPLATE. `_step_sensitivity`'s recursion
   `S = sum_k a_k C_{n-k} P_{n-k} + b Pq` differentiates a LINEAR-MULTISTEP update; a two-stage
   DIRK's per-step Jacobian is a COMPOSITION of two implicit solves, not a companion sum.
   Differentiating the two stage residuals w.r.t. the entering `x_n`:

       [C1 + (g h/2) G1] dY1 = [Cn - (g h/2) Gn] dxn
       [C2 + a33 h G2]   dY2 = A1 C1 dY1 + A0 Cn dxn

   so the step stores `(lu1, B1, lu2, C1, Cn, A1, A0)` and the matvec is two back-substitutions.
   THREE linearisation points per step (`x_n`, the internal stage `Y1`, `x_{n+1}`), each with its
   own stage coefficient -- which is why `_G_at` had to exist: recovering a physical `G` from the
   accepted step's stored `Geq = a h G` divides out one coefficient and mislabels the other two.

   ⚠ WHY THIS EARNS ITS KEEP. Trapezoidal is a second-order method whose SHOOTING monodromy is
   FIRST-order on a limit cycle: its opening manufacturing step is order-dropped to Euler and that
   seam sits inside the period map (see `_traverse_factored_plain`, `monodromy_twin`, the B16
   decision). TR-BDF2 is self-starting -- every step, first included, is the full two-stage method
   reading only `x_n` -- so there is no opener seam and the monodromy stays second-order without a
   Gear-2 twin. Measured against the pencil `exp(mu T)` on a source-free RC network: eigenvalue
   error falls at ratio ~4.0 per grid doubling (100/200/400 pts), and the full monodromy matrix
   matches `expm(A T)` in the 2-norm at the same ratio. A first-order map would halve, not quarter.

   ⚠ DELIVERED AS `factored_period_trbdf2(x0, T, npts)`, NOT `method='trbdf2'` in `solve`. The
   dense shooting Newton propagates its Jacobian through `companion_coefficients`, which the DIRK
   does not have; wiring TR-BDF2 into that Newton is a separate, larger change. The monodromy --
   the Floquet spectrum, which is the whole reason to carry a DIRK here -- stands on its own and
   goes first. Verified against the pencil in scratch BEFORE the formal test (the "wrong integrator
   gives plausible wrong numbers" trap).

T3. **Adaptive step control from the embedded 2(3) estimate** (this increment). The estimator
   coefficients were DERIVED, not quoted (the corpus has nothing on TR-BDF2): Taylor-matching the
   combination `est_raw = h (c0 f_n + c1 f_gamma + c2 f_{n+1})` to the order-2 solution's leading
   local truncation error gives `c0 = (1-sqrt2)/3`, `c1 = 1/3`, `c2 = -(2-sqrt2)/3` and the
   principal LTE coefficient `(4-3sqrt2)/6`. Verified against the analytic LTE of `y' = a y`:
   `est/true -> 1` as `h -> 0`.

   ⚠ THE RAW ESTIMATE IS A TRAP ON STIFF MODES, and the scalar `C=1` check HID a second trap. For
   `a h -> -inf` the raw `est_raw` grows like `|a h|` while the true error is damped to zero by
   L-stability -- a naive `||est_raw||` forces the controller to crawl through exactly the stiff
   transient the method exists to step over. H&S filter it through the stage matrix once. And the
   filter is `(C + a33 h G) Est = est_raw`, NOT `= C est_raw`: `est_raw` is in CHARGE units
   (`h dq/dt`), the mass matrix `C` maps state to charge, so the STATE error is `C^-1 est_raw` and
   the stiff replacement is `(C + a33 h G)^-1`. The extra `C` multiply (which the scalar `C=1`
   test could not see) left the estimate in charge units and 7 orders too small on a real circuit
   -- caught by checking the stored estimate against the true one-step error on an RC network,
   where `est/true -> 1` only after the `C` was removed.

   Delivered as a dedicated driver `Transient._run_trbdf2_adaptive`, reached for the non-fixed
   grid; the LMM controller path is untouched. Error per step, order 2:
   `err = rms(Est_i/(reltol|x_i| + atol_i))`, accept at `err <= 1`,
   `dt_next = dt clip(0.9 err^(-1/3), K, 1/K)` with `K = 1/2`. No `_push_history` (a DIRK reads no
   charge rings; only `accept_step` for stateful elements). `compute_lte` STILL raises -- it is the
   LMM divided-difference interface, which a stage method does not fit; the estimate is a stage
   combination computed in `_solve_timestep_trbdf2` and consumed by the driver. Verified: tighter
   reltol spends more steps and lands closer to analytic; and the adaptive waveform on a driven
   2-node RC matches a tight DOP853 reference to 8e-6 relative.

T4. **`method='trbdf2'` in the shooting Newton** (this increment). The dense Newton needs the
   monodromy, not `companion_coefficients`, so TR-BDF2 gets its own dense traversal
   `_traverse_trbdf2` (the same per-step two-stage map as the factored one, accumulated into a
   full `m x m` `P`, plus a period column `Pt` for the autonomous case). Both were
   FINITE-DIFFERENCE checked before use (roadmap 0j: the dT column has been got wrong here twice):
   `M` vs FD to 7e-12, `Pt` vs central FD to its truncation floor. Wiring:
   `_companion_reach` returns 1 for the DIRK (never solves-history); `x0_unknown` is forced True
   (self-starting, no manufacturing step); `func_trbdf2` (driven, `F = x0 - phi(x0)`,
   `J = I - M`) and `func_autonomous_trbdf2` (free period, with the `Pt` column and a phase row);
   `factored_period()` routes to the two-stage map; the replay skips the LMM seam/interior LTE
   (which does not apply to a self-starting method); and `monodromy_twin` EXEMPTS trbdf2 -- the
   twin exists to give a first-order one-step LMM a second-order monodromy, and the DIRK's native
   one is already second-order, so a Gear-2 twin would be pure cost and would hide its spectrum.
   Verified: driven RC matches the AC steady state and reports `rho = exp(-T/tau)` exactly;
   autonomous van der Pol converges to the LMM period (to O(h^2)) with a unit multiplier. Matrix-
   free trbdf2 shooting refuses (the monodromy is a dense two-stage product).

   ⚠ AN EXTERNAL CROSS-CHECK (docs session, from Bank 1985 + Hosea & Shampine + Kennedy &
   Carpenter, read firsthand) landed while this was built and CONFIRMED the constants
   independently: `gamma = 2 - sqrt(2)` is BOTH the one-LU condition AND the truncation-error
   optimum (Bank eq. 38, `C(gamma) = (-3g^2+4g-2)/(12(2-g))` minimised there at `~ -0.0404` =
   this file's derived estimator coefficient `(4-3sqrt2)/6`); and there is a factor-of-two
   convention trap in `gamma` between Bank/Hosea (the abscissa, our value) and Kennedy & Carpenter
   (its half). Both now noted in the integrator docstring. Index-2: HLR Theorem 5.9's no-order-
   reduction guarantee does NOT apply to TR-BDF2 (its `A` is singular -- the ESDIRK zero first
   row), but the docs session MEASURED full order 2 on a Hessenberg index-2 problem anyway, so no
   gate is written either way.

T5. **The adjoint transpose, and the whole AUTONOMOUS phase-noise stack** (this increment).
   `_monodromy_matvec_transposed_trbdf2` gives `M^T v` -- the transpose of the two-stage product,
   replayed in reverse step order (`M_j^T w = A1 B1^T K1^T C1^T z + A0 Cn^T z`, `z = K2^T w`) --
   plus `collect=True` returning width-`m` per-step `states` (`v(t_j) = Phi(T,t_j)^T v(T)`, no
   pair). Verified: `M^T` equals the dense `(M)^T` to 1.3e-15, complex splits into two real
   replays.

   ⚠ THAT UNLOCKS MORE THAN THE EIGENVECTOR. `ppv` reads only `states` (it discards `ts`), and
   its per-period PPV needs no pair reconstruction for a width-`m` map -- so `ppv`,
   `diffusion_constant`, and `oscillator_spectrum` ALL run over TR-BDF2, matching Gear-2: PPV
   vector to 3.5e-4, `Q` 0.1417 vs 0.1416, `c = 8.045e-08` vs `8.042e-08` on the noisy van der
   Pol, lineshape to 0.002 dBc across three decades. No Gear-2 twin -- the DIRK's native
   second-order monodromy carries it.

   ⚠ ONE `ppv` FIX WAS NEEDED. `ppv`'s `Q`/second-multiplier came from a matrix-free Arnoldi on
   `I - M`; on the small width-`m` DIRK map its clustered near-null spectrum left the unit root at
   `1 - 1.5e-6`, past the `1e-6` deflation, so it reported the ORBIT TANGENT as the second
   multiplier (`Q ~ 6e5`). The DIRK map is dense and small, so its exact eigenvalues are cheap:
   for `kind='trbdf2'` `ppv` now eigen-solves the densified `M` directly (deflates cleanly,
   `lam2 = 8.59e-4`). Gear/solved-history keep the Arnoldi byte-for-byte.

T6. **Still open: the DRIVEN forced surfaces.** `covariance`, `pnoise`, and PAC's adjoint sideband
   REFUSE loudly for TR-BDF2 (clear `NotImplementedError`, not a tuple-shape crash and not a
   plausible wrong number). A source injected into a two-stage step enters BOTH stages, so the
   per-step forced response and its noise covariance `Q_j` are two-stage quantities the LMM
   single-companion replay (`_forced_replay`, `_lyapunov_pieces`) does not represent -- a matching
   two-stage forward forced replay is the remaining derivation. Autonomous noise (the oscillator
   case, T5) is unaffected and complete.

   ⚠ THE `gamma` CONSTANT: one quadratic `Q(gamma) = gamma^2 - 4gamma + 2`, reached from FOUR
   different requirements (verified symbolically here; the docs session's first "three independent
   legs" framing was corrected -- they are the SAME quadratic, each `Q` times a factor NONVANISHING
   on (0,1), which is what makes `Q=0` the ONLY solution not merely one): one-LU
   `a33-gamma/2 = Q/(2(2-gamma))` (a one-char sign fix vs the first ad159a2 docstring); Bank eq.38
   `dC/dgamma = Q/(4(2-gamma)^2)`; Rosenbrock 1963 `d^2-2d+1/2 = Q/4` (d=gamma/2); and the radius of
   absolute monotonicity `R(A,b)=2(2-gamma)/(1+(1-gamma)^2)` (Bonaventura & Della Rocca 2015),
   `dR/dgamma = 2Q/(1+(1-gamma)^2)^2`, MAX `1+sqrt2` at the root -- verified here from the
   Kraaijevanger 1991 definition (validated on CN=2, BE=inf first); TR-BDF2 tolerates a ~21% larger
   step than trapezoid before monotonicity/positivity/TVD can fail. TR-BDF2 is L-stable for EVERY
   gamma (`R(inf)=0` identically), so L-stability alone selects nothing. All in the integrator
   docstring.

## TR-BDF2 as the default Floquet twin, 2026-09-05

T7. **`monodromy_twin` defaults to TR-BDF2, Gear-2 selectable** (committed bf6f891). The twin that
   supplies a one-step LMM's (trap/euler) Floquet/PPV/phase-noise quantities was always Gear-2;
   it now defaults to TR-BDF2. Measured against exact references: on a source-free damped LC
   (`exp(A T)`) and van der Pol (Abel's `exp(mu integral(1-v^2))`), TR-BDF2's lambda2/Q error is
   1-2 orders smaller than Gear-2's in the hard regime (coarse grid, high Q) -- 22% vs 0.27% at
   Q=100/N=50 -- and on the asymmetric `vdp+0.3u^2` fixture the advantage carries to the PPV-based
   `c` (1.4e-4 vs 1.7e-3) against the scipy-adjoint exact. Regime-dependent, not a fixed factor;
   at a fine grid or low Q both reach the O(h^2) floor. `pss.monodromy in {'trbdf2','gear','native'}`;
   twin-build factored into `_solve_twin(method)`, cached per method.

   ⚠⚠ IMPROVING ROBUSTNESS DEGRADED SAFETY -- the guard is the fix, and it is the review's most
   repeated finding arriving from the OPPOSITE direction. From a poor seed (euler at 400 pts, orbit
   55% off) the Gear-2 twin failed LOUDLY ('did not converge') but the more robust TR-BDF2 twin
   CONVERGED to a SPURIOUS limit cycle and reported Q=1.97 vs exact 5.91 with no error. `_solve_twin`
   now checks the twin's converged orbit against its seed (period + entering state); two convergent
   methods on one limit cycle agree to O(h^p) (measured 6e-6/1e-4 good vs 2.17/0.91 spurious), and
   the 0.25 gate -- set from the universal good-case agreement, not the fixture-dependent failure
   size -- refuses a too-poor seed. Keeping Gear-2 gives a fail-loud cross-check on the fail-silent
   default. Definitive test is refinement (spurious orbit fails h/2); this is the conservative
   stand-in.

T8. **The noise-injection surfaces fall back to Gear-2 (`_lyapunov_host`), and a real TR-BDF2 Q_j is
   QUEUED, not guessed.** covariance/pnoise/PAC-sideband need the per-step injection `Q_j`; a source
   in a two-stage step enters BOTH stages, so `Q_j` is a two-stage quantity. Until it is built, a
   TR-BDF2 Floquet source routes these surfaces to a Gear-2 twin of the same orbit (correct,
   validated). ⚠ MEASURED: the naive two-stage injection is WRONG BY ~27% (an O(1) bias, not
   asymptotic) on the scalar RC vs kT/C -- `Q = W1^2 S gamma h + W2^2 S (1-gamma)h` converges to
   1.268 kT/C, not kT/C. Building the stochastic DIRK stage weights by analogy misses the BDF2
   stage's previous-Wiener-increment coupling (Denk/Sickenberger/Winkler; the weight is
   `kappa^2/(2kappa+1)` = the deterministic parasitic root = ZERO_STABILITY_RATIO at
   kappa=1+sqrt2). So it is NOT good to inject yet (the user's own condition: "queue it as soon as
   you think it good").

   THE CORRECT ROUTE (peer research, firsthand): for ADDITIVE linearised noise the Lévy areas
   vanish and `Q_n = integral_0^h Phi(h,s) (CY/2) Phi(h,s)^T ds` is a DETERMINISTIC integral with
   an exact VAN LOAN oracle (exponentiate `[[-A, D],[0, A^T]] h`). ⚠ But MNA is a DAE (E singular);
   the nilpotent block differentiates white noise, so Van Loan must be applied AFTER projecting to
   the differential subspace (the "inherent regular SDE"; Demir 1996 propagates covariance exactly
   there -- "nodes connected to a capacitor"). Winkler's noise-free-constraint `im A_N subset im A_C`
   (a capacitive path in parallel with every noise source) is what keeps the covariance an ordinary
   process at all. Plan: project -> Van Loan on the differential subspace (2nd order) -> use with
   trbdf2's discrete `A_n`; validate against Van Loan-exact (matrix, per step) + kT/C + Monte Carlo.
   Scalar-RC Van Loan already verified 2nd order (0.9938/0.9987/0.9996 at h=0.5/0.2/0.1).

T8a. **When the trbdf2 Lyapunov IS built, one more test — and it is invisible to the obvious three**
   (docs session, MEASURED). If a noise source enters an ALGEBRAIC constraint (violating Winkler's
   `im A_N subset im A_C` -- e.g. a thermal-noise source between two nodes with no capacitance on
   either), the covariance's ALGEBRAIC entries diverge as EXACTLY 1/h (the nilpotent block
   differentiates white noise; discretised white noise has variance S/h). ⚠ The DIFFERENTIAL entries
   (the `P_xx` a kT/C test reads) stay perfectly healthy -- measured P_xx -> 0.5 correctly while
   P_yy doubles every step-halving. So per-step Van Loan (run on the projected system), kT/C, AND
   Monte Carlo are ALL blind to it (they look at differential entries; the MC shares the defect).
   The test: a fixture that violates `im A_N subset im A_C` must REFUSE (a cheap incidence-matrix
   rank check, pre-integration), or the runtime detector is refinement again (run h and h/2, assert
   no covariance entry doubles). ⚠ This makes the differential-subspace PROJECTION load-bearing, not
   a formality: Demir's capacitive-node-only propagation never forms an algebraic-node covariance
   and is structurally immune -- the projected route inherits that. This fixture is what catches a
   later "simplification" to the full MNA state; nothing else in the stack would.
   (Refinement is now the general instrument for THREE failure shapes in this arc: opener divergence,
   spurious twin orbit, and noise-on-a-constraint.)

## TR-BDF2 Lyapunov covariance built; pnoise sideband fold queued, 2026-09-05

T9. **The TR-BDF2 per-step injection is BUILT (`covariance`/`oscillator_covariance` native)** --
   the deterministic DAE-projected Van Loan integral, no stochastic stage weights (additive noise
   -> Levy areas vanish). `_vanloan_step_injection`: split the reduced `(C,G)` into differential
   (capacitive) and algebraic rows, Schur-complement the algebraic ones out, route the algebraic-row
   noise into the differential rows through the same elimination (`R_proj`), Van Loan on the
   differential subspace (`expm([[-A,D],[0,A^T]] h)`), embed back. `_lyapunov_pieces_trbdf2` pairs it
   with the dense two-stage step map `A_n`. VALIDATED: converges to kT/C at SECOND order on R||C
   (ODE) and VS-R-C (DAE) -- 3.8e-4 at 100 pts vs Gear-2's 6.9e-2 (~180x, and 2nd order vs gear's
   1st). ⚠ The assertion is the RATE (~4x/doubling), NOT machine precision -- a machine-zero kT/C
   would mean a method-consistent `Q=P(1-A^2)` fudge that corrupts the transient covariance (peer
   trap). `_lyapunov_host` now returns the monodromy twin (trbdf2 by default; gear selectable);
   Winkler's `im A_N subset im A_C` handled structurally by the projection (immune to the 1/h
   algebraic divergence -- it never forms an algebraic-node covariance).

T10. **pnoise native over TR-BDF2 is QUEUED, not shipped -- the sideband FOLD is the blocker, and it
   was measured.** The forward forced replay and its CHAINED two-stage transpose were built and are
   dual-consistent to 1.8e-16 (`<lam, J du> == <J^T lam, du>`, exact on the rectangular operator).
   BUT `adjoint_sideband_row`'s forced part is `-sum_j phase[j] ts[j]` -- ONE source-injection time
   per step (the endpoint) -- while a two-stage step injects the source at THREE abscissae
   (`t_n`, `t_n+gamma h`, `t_{n+1}`) with three phases, which that fold cannot represent. Shipped it
   naively: measured 1.25e11 relative error on the linear divider (99 spurious sidebands) vs the
   AC-noise reference -- so it was reverted rather than shipped. pnoise over a TR-BDF2 Floquet source
   now falls back to a GEAR-2 twin (`_adjoint_host`), correct (rel err 2.65e-8, matches gear). The
   remaining piece: extend `adjoint_sideband_row`'s forced fold to inject at the three stage
   abscissae with their phases (the per-step collect must carry the stage-2 solve `z` AND the
   stage-1 feed `A1 p`), gated on the dual-consistency test against `m` forward driven solves. The
   forced replays are the validated building blocks; only the fold's one-time-per-step assumption
   needs lifting.

## pnoise native over TR-BDF2 -- the two-stage sideband fold, 2026-09-05

T11. **pnoise is now NATIVE over TR-BDF2** (the "advanced substantial extension", requested).
   The blocker was `adjoint_sideband_row`'s forced fold `-sum_j phase[j] ts[j]` -- one source-
   injection time per step -- while a two-stage step injects at THREE abscissae. The fix is the
   TWO-VECTOR fold (`_sideband_forced_trbdf2` + `_forced_replay_transposed_trbdf2`): the reverse
   pass reads the source coupling through BOTH stages per step --

       t3 = K2^-T lam ;  p = K1^-T (C1^T t3) ;  t2 = A1 p
       forced -= a33 h e^{jw t_{n+1}} t3
       forced -= (gamma h/2)(e^{jw t_n} + e^{jw(t_n+gamma h)}) t2
       lam <- A1 B1^T p + A0 Cn^T t3          # monodromy transpose
       lam <- lam + e^{-j(l w0 + w) t_n}/N d   # output injection, AFTER the update

   ⚠ THE INJECT ORDER IS LOAD-BEARING: the output injection is added AFTER the step's costate
   update (so the output at step n couples to sources at steps < n -- causality), and the source
   coupling is read BEFORE it. With the inject added first the fold was 4.8-52% wrong; with the
   correct order it matches `m` forward driven solves on the diode mixer to MACHINE PRECISION
   (1.5e-16 .. 1.1e-15 over l = 0/1/-2/3). The peer independently validated the same two-vector
   fold on an LPTV system (4.49e-16; endpoint-only 42% wrong), and the (1+e^{jw gamma h})/2 shape
   is the two TR abscissae combining into one per-sideband constant.

   ⚠ TWO GATES, TWO LEVELS (peer): the CHAINED TRANSPOSE is exact at the STEP level
   (`<lam, J du> == <J^T lam, du>`, dual-consistent to 1.8e-16) while the FOLD was still 42% wrong
   -- passing the step-level test says nothing about the fold, so both are kept. End-to-end:
   trbdf2 pnoise reduces to the AC-noise analysis on the linear divider (rel err 1.2e-9, better
   than gear's 2.7e-8, stops on the ratio test) and folds on the diode mixer to match gear within
   0.01% (79 sidebands). `_adjoint_host` now returns the monodromy twin (no gear fallback);
   `adjoint_transfer_row` also works over TR-BDF2 as a bonus (single-instant, the W^T alone).

   Still refusing over TR-BDF2: the FORWARD `_forced_replay` (PAC.solve forward) -- its forward
   fold has the same three-abscissa structure and is not yet validated; out of scope for pnoise,
   which needs only the reverse path.

---

## Radau IIA(3) across the PSS stack — order 5, fully implicit, 2026-09-06

Radau IIA 3-stage added the way TR-BDF2 was, one validated commit at a time. Where
TR-BDF2 earns its place by the failure modes it lacks (order 2, but no opener seam, no
`(-1)^n` mode, no step-ratio limit), Radau earns its place by **coverage and accuracy**:
order 5, stage order 3, L-stable, stiffly accurate, and `det A = 1/60 != 0` — the one
candidate the H&W VI.2 Thm 2.3 index-1 DAE convergence result actually covers (a singular
`A`, any explicit-first-stage DIRK included, is not).

**The tableau (gate R1).** `c = [2/5 - √6/10, 2/5 + √6/10, 1]` (the two interior Radau
points and the stiffly-accurate endpoint), `A` the 3×3 collocation matrix, `b` = last row
of `A`. Constants DERIVED from √6 in the class, not transcribed, and self-checked against
the reference relations: row sums == c, det A == 1/60, B(1..5) (order 5), C(1..3) (stage
order 3), R(∞) == 1 - bᵀA⁻¹1 == 0 (L-stable), stiff accuracy (b == A[-1], c₃ == 1), and
eig(A⁻¹) == the stored cost-transform eigenvalues (one real γ_r = 3.6378, one pair
α±iβ = 2.6811 ± 3.0504i). All to machine precision.

**The step is FULLY IMPLICIT — one coupled 3n solve, no DIRK shortcut.** The three stages
couple into ` J_block[i][j] = δ_ij C(Y_i) + h A_ij G(Y_j)` (= `I₃⊗C + h A⊗G` in the linear
case); `x_{n+1} == Y₃` by stiff accuracy. `Transient._solve_timestep_radau` solves it with
a dense coupled Newton. The `A⁻¹`-eigenbasis transform (1 real + 1 complex LU, needing the
complex `klu_z_*` binding) is the documented EFFICIENCY follow-up — the coupled real solve
is correct and, for the small circuits here, cheap. Gate: order 5 on the analytic RC step
through the real loop — ratios 30.6/31.1/31.5 → 32 per doubling.

⚠ **LIMITING IS LOAD-BEARING, and the pnoise cross-check is what found it.** The hand-rolled
coupled Newton first shipped WITHOUT junction limiting (TR-BDF2 gets it free — each stage
runs through `self._newton`). Without it the diode Newton overshoots the exponential and
settles on a spurious near-linear solution: a diode mixer came out a PURE SINUSOID
(H0,H2,H3 ~1e-16 against the reference's 0.39/0.18/0.046), and it reported `converged`.
Found ONLY because pnoise was validated against gear (a reference the Radau path cannot
influence) rather than self-consistency — the sideband rows for l≠0 were ~1e-12 and even
agreed with a brute-force forward of the SAME (wrong) discrete model. Fixed by
`cir.limit(Y_trial, Y_prev)` per stage per iteration; Radau transient then matches gear's
rectified spectrum to the digit. The lesson is §0j again: a self-consistent instrument
(adjoint == transpose of my own forward) proves nothing about physical correctness.

**The factored shooting monodromy (gate R2).** `_traverse_factored_radau` stores one dense
`3m×3m` factor + `Cn` per step; the matvec stacks `[Cn v; Cn v; Cn v]`, solves, and reads
the third `m`-block. Self-starting, so the map is `m×m` and order-5 round the whole period
— no opener, no solved-history pair. `_monodromy_matvec_transposed_radau` gives
`Cnᵀ(p₁+p₂+p₃)` with `p = J⁻ᵀ[0;0;w]`; `ts[j]` keeps the full `3m` coupled solve for the
fold. Gate (same source-free RC as the TR-BDF2 pencil test): eigenvalues vs `exp(μT)` —
1.5e-9 at 25 pts, ratios 31.3/31.6 → 32, and the adjoint is the exact transpose (4.6e-16).

**`method='radau'` in the dense shooting Newton (gate R3).** `_traverse_radau` propagates
the dense monodromy `P` and, for the free-period system, the period column `Pt = dx/dT`
from d/dT of the stage residuals (`h_j = frac_j T`). `func_radau`/`func_autonomous_radau`
mirror the TR-BDF2 pair. `_companion_reach`, the method whitelist, `x0_unknown` forcing,
`_want_lte`, and `monodromy_twin`'s self-sufficient set all learn `'radau'` (a self-starting
method reads its own native map, no twin). Gates: dense `P` vs FD 1.4e-7 and `Pt` vs FD
2.4e-5 on van der Pol (the dT column checked before use — §0j, wrong twice before); driven
RC PSS fundamental matches AC to <1e-3 and spectral radius to exp(-T/τ) at 1e-6; autonomous
van der Pol converges to gear's period to <1e-3 with a unit multiplier.

**The whole small-signal stack rides on the factored map, without a twin.**
- ppv/diffusion/oscillator_spectrum: the dense width-m stage map gets its spectrum by m
  matvecs + eigvals (the Arnoldi mis-resolves the unit-root cluster), and the coupled forced
  adjoint `_forced_replay_transposed_radau` (W^T xa) carries the source through all three
  stages — `acc += -h Σ_k exp(jw t_{n,k}) Σ_i A_ik p_i`, verified EXACT by dual consistency
  `<xa, W u> == <Wᵀxa, u>` to 7.7e-16. Gate: c and the lineshape match gear on van der Pol.
- covariance/oscillator_covariance: `_lyapunov_pieces_radau` reuses the SAME
  `_vanloan_step_injection` (the exact continuous per-step Van Loan integral is
  method-independent); only the transition A_n changes. On an LTI RC the injection is exact,
  so the covariance error is the transition's order — O(h²) for TR-BDF2, O(h⁵) for Radau:
  measured ~30x per doubling, ~1e-9 at 100 pts, >20x closer than gear.
- pnoise: `_sideband_forced_radau` is the coupled three-vector fold (no `A⊗B` shortcut),
  the injected sibling of the forced adjoint; injection AFTER the costate update (causality).
  Gate: linear divider == AC-noise to <1e-6 (stops on the ratio test); diode mixer folds >5
  sidebands and lands on gear to <5e-3.

**The one deferred parity piece: adaptive step control.** The embedded 5(3) estimate needs
the radau5 `dd` weights (a published constant with a free parameter fixed for L-stability);
the FILTER is known (the transform's real factor) but reconstructing the weights from memory
risks shipping a wrong instrument, which §0j forbids. DEFERRED, not guessed:
`_solve_timestep_radau` raises under the (never-set) want-est guard, and `_solve` refuses the
adaptive grid for Radau (fixed_timestep=True, or TR-BDF2). It would be built the way
TR-BDF2's 2(3) estimate was — derive, then validate the estimate/true-LTE ratio → 1 and the
adaptive-reltol gate — before it is trusted. Also still deferred (shared with TR-BDF2): the
FORWARD `_forced_replay` (PAC.solve forward), out of scope for pnoise's reverse path.

---

## Radau IIA(3): the adaptive 5(3) estimator and the cost transform, 2026-09-06

The two pieces deferred from the first Radau arc, now built and validated.

**The adaptive 5(3) estimator (`_radau_error_estimate` + `_run_radau_adaptive`).**
Hairer & Wanner's radau5 estimator: a lower-order (order 3) embedded solution differs
from the order-5 step by a combination of the three stage increments ``Z_i = Y_i - x_n``
plus a fictitious explicit stage ``f(x_n)``:

    F1  = (dd1 Z1 + dd2 Z2 + dd3 Z3)/h
    rhs = C(x_n) F1 + f0,     f0 = -(i(x_n) + u(t_n))
    est = ((gamma_r/h) C(x_n) + G(x_n))^{-1} rhs            (STATE units)

with ``dd1 = -(13+7√6)/3``, ``dd2 = (-13+7√6)/3``, ``dd3 = -1/3`` and ``gamma_r`` the real
eigenvalue of ``A^{-1}``.  ⚠ THE FILTER (the transform's real factor) is what makes it
stiff-robust: an unfiltered embedded difference grows like ``|λh|`` on a stiff mode while
the true error is L-damped to zero; the real-factor inverse maps it back to a bounded state
error.  The controller uses the ``1/(3+1) = 1/4`` step exponent (the peer's ``1/4``).

⚠ THE dd WEIGHTS WERE VALIDATED, NOT TRUSTED (0j -- the whole reason this was deferred
rather than guessed last round).  Two references the estimator cannot influence: on a smooth
RC (VSin, consistent IC) the single-step estimate at the RC node falls as ``h^4`` (ratios
6.4/13.2/15.1 → 16) -- a wrong lower-order construction would give 8 (``h^3``); and on a
diode rectifier the accepted-step count rises monotonically as ``reltol`` tightens
(160/478/1697/4895 over 1e-2..1e-8), so the estimate is steering.

**The cost transform (`ComplexKLUSolver` + `_solve_timestep_radau_transformed`).**  Opt-in
via ``_radau_use_transform``; the dense coupled ``3n`` solve stays the default and the
correctness reference.  Multiplying the coupled system by ``(A^{-1} ⊗ I)/h`` and
diagonalising ``A^{-1} = V diag(λ) V^{-1}`` decouples it into ``(λ_k C/h + G) w_k = rhs_k``:
the real eigenvalue → one REAL ``m×m`` solve, the complex pair → one COMPLEX solve plus its
free conjugate.  An ``O((3m)^3)`` dense solve becomes two sparse ones.

- ``ComplexKLUSolver`` binds the ``klu_z_*`` family (``klu_analyze`` is pattern-only and
  shared with the real path; KLU's packed-complex layout is exactly numpy ``complex128``, so
  values marshal by ``.view(float64)``).  Validated to 4e-16 vs a dense complex solve, with
  the analyze-once/refactor-many reuse KLUSolver uses.
- The transform step is SIMPLIFIED Newton (Jacobian frozen at ``x_n``, full per-stage
  residual, junction limiting), reusing the two factorisations every iteration.  ⚠ On a
  strongly nonlinear step the frozen Jacobian can stall; then ``_solve_timestep_radau`` FALLS
  BACK to the dense full-Newton solve, so the answer is never wrong, only occasionally slower.

Gates: transform == dense to 4.7e-17 on a linear RC ladder (0 fallbacks) and to 2.3e-12 on
the diode mixer (2 of 160 steps fall back at the diode-switching instants, harmonics
preserved).  Speedup grows with size: 1.18x at m=202, 1.92x at m=402 (load 0.74, so
wall-clock is trustworthy) -- the dense ``3m`` solve is cubic, the transform two sparse
solves, so the gap widens with ``m``.

**Remaining deferrals (shared with TR-BDF2):** the FORWARD ``_forced_replay`` / ``PAC.solve``
forward, out of scope for pnoise's reverse path.  Radau's parity with TR-BDF2 is otherwise
complete, two orders higher, and now fast on a large sparse circuit.
