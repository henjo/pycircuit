# Periodic small-signal and noise analyses — reference (2026-09-03)

What is built on top of `PSS`, how to call it, what each capability was gated against, and
— at more length, because it is the part that cost the most to learn — **where each one
stops being valid and what it refuses to do**.

`doc/pss_roadmap_260902.md` is the plan. This is the built surface.

---

## The chain

Everything here is one object away from the shooting solve:

    PSS.solve()
      └── PSS.factored_period()      one period, kept FACTORED  (lazy, cached)
            ├── fp.matvec(v)              M v          real or complex
            ├── fp.matvec_transposed(v)   Mᵀ v         Gear-2 only
            └── fp.samples / fp.steps     per-step factorisations
                  ├── PAC.solve()                    sideband spectrum, swept
                  ├── PAC.adjoint_transfer_row()     every source → one output
                  ├── PAC.adjoint_sideband_row()     …at one sideband
                  ├── PAC.pnoise()                   time-averaged output noise
                  ├── PAC.am_pm()                    AM/PM split of a modulation
                  └── PSS.ppv()                      perturbation projection vector

**Nothing forms an operator.** The withdrawn PAC allocated the whole `(N m) × (N m)`
system — 419.5 GiB at `N = 137`, `m = 1000`. The per-step factorisations `PSS` already makes
*are* the preconditioner, and what remains is `m × m`.

⚠ `factored_period()` re-traverses at the **converged** solution rather than reusing the
Newton's factors. The last `build` call inside the Newton is at the last *trial* iterate;
reusing it would give an operator for a trajectory near the solution instead of at it.

---

## PAC — periodic AC

```python
pac = PAC(cir)
res = pac.solve(pss, freqs, recycle=True)      # sideband spectrum vs frequency
```

After preconditioning the system is `(I − αM) y₀ = α w(f)`, with `w` the forced response
over one period from a zero initial state.

**Gated** against `analysis_ss.AC` on a linear circuit, where the LPTV response collapses to
the LTI one. Per doubling of the period grid, at 700 Hz — **the rate is the assertion, not
the size**:

| method | rate | rel @ 250 pts |
|---|---|---|
| gear | **4.00× (O(h²))** | 1.65e-04 |
| euler | 2.00× (O(h)) | 1.40e-02 |
| trap | 2.00× (O(h)) | 4.13e-03 |

⚠ **Trapezoidal is second order and its PAC is first.** Not the `(-1)ⁿ` mode — those sit at
0, annihilated by the Euler opening step — and not the trajectory, whose waveform converges
at ~4.2×. It is the **manufacturing step**, which lives outside the traversal loop and so
carries no small-signal source: one step of `u` out of `N`. `x0_unknown=True` restores
4.00×; euler is unmoved either way, which is the control. **`PAC` warns.**

`recycle=True` shares one Krylov subspace across the sweep — `A(α) = I − αM`, so the space is
`M`'s and α-independent. Measured on RC ladders, 24 frequencies: 72→6, 168→15, 302→26
matvecs (~11.6×), agreeing with the per-frequency solve to 1.2e-13.

---

## The adjoint rows — many sources, one output

```python
row  = pac.adjoint_transfer_row(pss, freq, output)          # response at t = 0
rows = pac.adjoint_sideband_row(pss, freq, output, [0,1,-1])  # H_l coefficients
```

`d^T y₀ = α ((I − αM)^-T d)^T W u`, so **one transposed solve** replaces one forward solve
*per source* — the asymmetry pnoise is shaped by, and unhelped by recycling because the
right-hand side is what changes.

| m | adjoint | forward (m solves) | speedup | agreement |
|---|---|---|---|---|
| 6 | 0.013 s | 0.090 s | 7.0× | 9.9e-16 |
| 14 | 0.020 s | 0.346 s | 16.9× | 6.9e-16 |
| 32 | 0.032 s | 1.299 s | 40.0× | 9.3e-16 |

The speedup growing *linearly* with `m` is the claim; forward is O(m) solves, this is O(1).
It needs no reverse integrator — the stored factorisations already solve transposed.

⚠ **A sideband is a functional distributed over the period**, so its adjoint takes an
injection at every step *and* picks up a second term for the initial state's own source
dependence. The two are comparable (140 against 749 at `l = 0`), so an implementation with
only the first returns a plausible number.

⚠ **The sideband convention is the DFT of `v(t) = y(t)e^{−jωt}`, not of `y`.** Decomposing
`y` is self-consistent with a matching forward reference to 1e-15 and still wrong — it smears
a Dirichlet kernel across every `l`. Only a **linear** circuit, which converts nothing,
catches it.

⚠ **`|l| ≤ N/2` is a refusal**, not a tolerance: nothing aliases down from above the grid's
own Nyquist. Use a finer grid.

---

## pnoise — time-averaged output noise

```python
S, sidebands = pac.pnoise(pss, freq, output)
```

`S(f) = Σ_l h_l CY h_lᴴ` with `h_l = H_l(f − l f₀)`; white sources in disjoint bands are
uncorrelated, so the bands add in power.

**Gated** against `analysis_ss.Noise` — a different analysis, an `(sC + G)` solve, no
monodromy or period in it. On a linear circuit the fold collapses to the stationary formula
(Okumura's `p = 1`): ratio **1.000000**, every `l ≠ 0` term ~1e-32. On a diode mixer the fold
carries **62%** of the output noise, which is why the linear gate is not sufficient alone.

⚠ **"Time-averaged" is the specification, not a hedge.** Two mechanisms make output noise
cyclostationary and only one is about the sources: bias-dependent sources (refused, below),
**and the periodic source-to-output transfer function**, which applies even when every source
is stationary. The sum handles the second and returns its time average — right for most uses,
**incomplete for two ordinary RF topologies**: a nonlinear subsequent stage, and cascaded
stages off a shared reference. A scalar per frequency also cannot carry the correlation at
`k f₀` that cyclostationary noise has.

⚠ **`maxsidebands` is an accuracy knob here and a reporting knob in `PAC.solve`** — the
opposite of how it reads. A driven signal lives at the frequencies it is driven at; noise
lives at all of them. Capping it always makes `S` a **lower bound**.

⚠ **Which stopping rule fired is part of the answer.** The ratio test means the series
converged; the grid's Nyquist means the grid ran out first and sidebands above it are
**missing rather than small**. It warns.

⚠ **Stationary sources only, and it checks.** A bias-dependent `CY` makes the sources
cyclostationary; the windows then correlate the sidebands, they stop adding in power, and the
cross terms need a construction that is not built. It samples `CY` along the orbit and
**raises**. Unreachable in the discrete element library — every `CY` there ignores `x`.

---

## PPV — perturbation projection vector

```python
v, info = pss.ppv()        # info: samples, xdot, q, residuals, second_multiplier
```

Two bordered solves, both matrix-free, **no eigendecomposition** — which is the whole content
of the 2003 improvement, and matters most where a PPV is wanted, since multipliers crowd 1 on
high-Q oscillators and a bordered solve never has to tell candidates apart.

`q = C(0)ẋ(0)` is **exact, not differenced**: the DAE gives `q = −(i(x₀) + u(0))`.
`info['border_residual']` ~1.4e-11 is the free correctness check — a nonzero `y` means the
border absorbed a residual belonging to the null space.

**Gated physically**, because an identity check cannot settle a scale: displace van der Pol,
integrate the full nonlinear system until the transverse modes die, read the surviving
displacement along the tangent. Fitted scale → **1.0017** at 800 points, residual O(h).

⚠ **The normalisation is `v·ẋ(0) = 1`, not the transcribed `v·q = 1`** — a 7% error. The
bordered solve returns what behaves as `Cᵀv₁`, so it contracts with a state perturbation
directly. Independently confirmed from the equation: `⟨Ω₂Q, V₁⟩ = 1` with `Ω₂Q = C u₁`.

⚠ **The obvious repair is also wrong.** Predicting a state jump's shift as `vᵀCδ` gives
residuals 0.36 / 0.40 / 0.42 that **grow** under refinement, with per-direction ratios from
−0.44 to 28.7. Do not change it back without re-running that experiment.

⚠ **VALIDITY BOUNDARY — SLOW NODES, and the gate above cannot see it.** The phase equation
treats the frequency response as **instantaneous**; a slow node *filters* the noise of devices
near it and the PPV cannot see the filtering. The literature names this test's own regime as
the blind spot: the equation "was verified … on small, simple oscillators, [with]
perturbations applied to oscillator cores." The extraction degrades too, and silently:

| τ/T | \|λ₂\| | σ_min(bordered) | null residual |
|---|---|---|---|
| none | 0.000856 | 8.62e-01 | 4.1e-11 |
| 1e2 | 0.990049 | 4.49e-03 | 4.6e-11 |
| 1e6 | 0.999999 | 4.47e-07 | 4.4e-11 |

`σ_min` tracks `T/τ` over six decades **while the residual does not move**. So `ppv()`
estimates `|λ₂|` by deflated power iteration and **warns**, saying the result is an upper
bound — no residual can report this.

### Q — one number that subsumes four diagnostics

`info['Q']` comes free with the PPV, from the `|λ₂|` the deflated power iteration already
produced. An amplitude perturbation decays to `|λ₂|` of its size each cycle, so the cycles
needed to fall below a threshold **is** the oscillator's Q — reported for a `1/e` threshold,
`Q = −1/log|λ₂|`.

The usual definitions do not apply to an autonomous circuit: `f_r/Δf` presumes a Bode plot of
a BIBO-stable linear system, stored/dissipated presumes damping a self-sustaining oscillator
does not have, and it is **not** the Q of the resonator inside it.

⚠ **"High Q", "a second multiplier near 1", "slow amplitude restoration" and "a long time
constant" are four vocabularies for one condition.** That is why the same circuits defeat the
phase row, the eigen-split, the probe's continuation and the PPV's instantaneous-response
assumption — not four coincidences. This is the cheapest early warning available: it costs
nothing beyond what `ppv()` already does.

Validated by construction — a slow node built with a known `τ/T` has multiplier `exp(−T/τ)`,
so Q must come back as `τ/T` itself: **10.00** and **99.99** for 10 and 100, against 0.14 for
plain van der Pol.

---

## AM/PM

```python
m_am, m_pm = pac.am_pm(pss, freq, output, carrier=1)
```

`m_am = a + conj(b)`, `m_pm = a − conj(b)`, with `a`, `b` the upper and lower sidebands over
the carrier phasor. **No new solve** — a change of basis on `adjoint_sideband_row`.

⚠ **The conjugate is the whole thing**, and `a ± b` looks equally plausible: the sidebands
counter-rotate, so their sum traces an ellipse, and without the conjugate a rotating ellipse
reports as pure AM. ⚠ **The two sidebands are not conjugates of each other** — that is LTI,
and LPTV analysis exists because it is not — so two solves at `±freq`, not one and a
conjugate.

Verified: pure AM → PM index exactly 0 and vice versa; a driven diode detector gives 0.005
PM/AM (no free phase to modulate), flat across offsets; an oscillator's PM/AM ratio grows
3.90×/3.91× per 4× offset reduction — the `1/ω_m` divergence.

⚠ **Open:** on an *autonomous* circuit the sideband rows come back at ~1e-12 in absolute
terms. The ratio is right and the driven case is healthy at O(10³) — but the magnitude is
**not established**. Do not read an oscillator AM/PM magnitude from this yet.

---

## What every one of these refuses

Thirty-three refusals in `shooting.py`. The ones a caller will meet:

| refusal | why it is not a tolerance |
|---|---|
| PAC on a **harmonic** of an oscillator | `I − e^{−jωT}M` becomes `I − M`; a perturbation along the orbit gives unbounded phase drift. Measured 2.8e-11 at every harmonic, linear in the distance, and never below 0.78 on a driven circuit. |
| An operating point from a **different circuit** | The injection *device* is present with its *signal* at zero, so the bare oscillator's orbit is the wrong base solution. The reference-node check cannot catch it. |
| `ppv()` on a **driven** circuit | Its phase is the source's; there is no unit multiplier to project onto. |
| `ppv()` / transposed replay on a **one-step method** | The reverse recursion is derived for a two-step companion. |
| `pnoise` with **bias-dependent `CY`** | Cyclostationary sources stop adding in power — a different model, not a harder sum. |
| AM/PM at a **harmonic with no carrier** | Relative to the signal, not against zero: dividing by a 1e-16 phasor reports an enormous modulation of nothing. |
| Sidebands **above the grid's Nyquist** | Nothing aliases down from above what the grid represents. |

**GMRES is judged by residual, not by its status flag.** SciPy reports breakdown on these
small operators where it has already solved them — the Krylov space is exhausted and the next
vector is numerically zero. Trusting the flag turned exact answers into errors at small
offsets, which is where phase noise lives.

---

## Method requirements

| capability | method |
|---|---|
| `PAC.solve`, `pnoise` | any; **gear** for the method's own order |
| adjoint rows, `pnoise`, `ppv`, `am_pm` | **gear only** — they need the transposed replay |
| PAC accuracy on the plain path | `x0_unknown=True`, or gear |

---

## Open, and stated as open

- **Oscillator AM/PM magnitude** — ratio right, absolute rows ~1e-12, cause unestablished.
- **A2's slow-node validity boundary** — detected and warned, not fixed. The fix is a
  frequency-aware PPV, which is this same bordered system at `ω_s ≠ 0`; measured, the offset
  puts a floor under the conditioning at ~`ω_s/f₀`, worth ~1–3 orders in the regime of
  interest.
- **Cyclostationary sources** — refused, not supported.
- **Trapezoidal over long runs** — it has no numerical *damping*, which is a measured result
  and stands; but a constant-bias error growing as `t²` eventually drags the amplitude error
  up with it. "No damping" is not "no amplitude error over many periods", and it bears on
  warm-start length.
- **Oscillator phase noise / jitter output** — a different shape entirely (closed form in a
  few scalars, no sweep). Not built. ⚠ And `S_phi` must not be reported near the carrier.
