# Handover

> **⚠ THIS FILE COVERS TWO THREADS.**  The CURRENT state is immediately below.
> Everything from "Handover — the leapfrog/symbolic work" onward is the OLDER
> thread (2026-07-30) and is kept because it is still resumable, not because it
> is current.

---

## CURRENT STATE — 2026-09-17

**Branch `cna-jax-vectorization`.**
⚠ **The commit/push state is deliberately NOT written here.**  A line that
records it is changed by the commit that records it, so it is stale the moment
it lands — this file carried a wrong count twice in one day (first six days
behind, then one commit behind in the other direction, after a push).  Ask git,
which is the only statement of it that stays true:

    git rev-list --count origin/cna-jax-vectorization..HEAD   # 0 = nothing unpushed
    git log --oneline -1

Merging to master remains the repo owner's call and has NOT been done.
⚠ Push only when asked, and only this branch.

Suite: **3235 passed, 6 skipped, 3 xfailed** — **27 min as ONE run**.

✅ **RUN IT DETACHED, IN ONE PIECE** (measured 2026-09-18, 1631 s):

    setsid nohup .venv/bin/python -m pytest pycircuit -q -n 4 -p no:randomly \
        > suite.log 2>&1 < /dev/null &

then poll the log rather than holding the foreground.

⚠⚠ **THIS REPLACES THE FIVE-WAY CHUNKED RECIPE.** That recipe existed only
because six background runs in one day were killed as "low memory" with ~20 GB
free and PSI zero, and because `test_B16` alone runs **785 s** against this
harness's 600 s cap. `setsid nohup` puts the run outside the harness's process
group and it simply survives — `test_B16` included. It is also FASTER than
chunking (~36 min) and, unlike a split, it cannot silently omit a test.
(Credit: the private suite's session, which measured the same kills at 18–19 GB
available against a 260 MB peak and found detaching, not retrying, was the fix.)

⚠ **WHEN YOU POLL, EXIT ON THE PID DISAPPEARING AS WELL AS ON A RESULT LINE.**
A kill and "still running" look identical to a grep that only matches success;
a watcher that greps only for `passed` hangs forever on a crash.
⚠ **Use `kill -0 <pid>`, NOT `pgrep -f`** — `pgrep` matches its own shell, and
it has produced a false "process gone" and a hung `until` loop repeatedly here.
Take the PID from `ps`, not from a `pgrep` list: the list can hold transients.
⚠ Run ONE suite at a time; `-n 4`–`6` is what survives, `-n 10` OOMs.

**Fallback, if a detached run is ever killed after all:** the old five-way
chunking — everything except `test_analysis_shooting.py`, plus that file split
four ways by node id, `-n 4`. ⚠ RE-COLLECT the node ids before splitting
(`--collect-only -q`): a stale split once produced four green chunks that
silently omitted the very test being committed.

**The active thread is PSS / shooting / DAE integration**, and it is now TWO
files (split 2026-09-11, because the plan was buried in 14.7k lines of log):

  - `doc/pss_roadmap_260902.md` — **the PLAN**: the organising fact, §A
    capabilities, §B decisions awaiting a call, §C closed, §D the failure
    shapes.  It ends with a 128-entry INDEX of the log.
  - `doc/pss_log_260902.md` — **the RECORD**: what was built and measured, in
    order.  ⚠ Code docstrings reference sections BY NAME, and some of those
    names are here rather than in the plan; the plan's index resolves them.

Almost everything in §A and §B is now marked BUILT or CLOSED.  What is
genuinely open:

| item | status |
|---|---|
| A4b output layer (AM/PM as a change of basis) | unbuilt; citations, not measurements |
| A8 sampled / edge-jitter noise | OPENED 2026-09-16, first number GATED in the suite, unbuilt as a feature. ⚠ "point `sampled_noise` at the crossings" is WITHDRAWN — it refuses autonomous PSS by design. Oscillator share = `c` (n-scaling reproduces it to 1.0005). ✅ The ADDITIVE number needs no new machinery: the chain is DRIVEN given the waveform — 3 tanh buffers, gm·R=0.5, PSD 4e-21, f0=1 MHz → **σ_t 8.27e-11 s**, band [1e3, 5e5], test `test_the_additive_edge_jitter_is_method_independent_and_nothing_manufactures_its_floor`. ⚠⚠ THREE CLAIMS IN THE FIRST PASS WERE WITHDRAWN, all mine: σ_t 5.712532e-11 was GRID-BOUND (pole at 80·f0 → the grid's Nyquist set the noise bandwidth, +14.7/+8.7/+4.9 % per doubling; fix = pole below f0); the "242× sources-OFF floor" was RESISTOR THERMAL NOISE (`R.CY` = 4kT/r), predicted exactly as 1+PSD/(4kT/R) = 242.312 vs 242.311695 — the TRUE control is `R(noisy=False)`, which gives **exactly 0.0**; and "σ_t·slew constant" is VACUOUS since σ_t ≡ √var/\|slew\| (PSD×4 → 1.99690 is likewise arithmetic). Gate is now CROSS-FAMILY (radau/gear/trap within 0.17 %). ✅ PLUS AN ABSOLUTE ANCHOR (the peer named the gap: cross-family agreement cannot show a number is RIGHT): a ONE-STAGE LINEAR fixture has a closed form — `var = (S_inj+4kT/R)·R/(4C)`, band `1−fmin/(f0/2)`, fold captured `(2/π)arctan(f_N/f_c)` — matched to **3.7e-04** at npts 400 and 800, slew to 1.1e-05, σ_t 3.614654e-11 vs analytic 3.618558e-11; the Nyquist deficit was PREDICTED (2.533e-03 → 1.266e-03, halving) and measured 2.178e-03 → 8.976e-04, test `test_the_edge_jitter_of_a_linear_stage_matches_its_closed_form_and_the_grid_truncation_it_names`. ✅ **AND σ_t IS NO LONGER A DEFINITION — the Monte Carlo is BUILT and it agrees**: a noisy transient's crossings scatter by `σ_t(MC)/closed form = 1.01099 ± 0.01156` pooled over 9 seeds × 3 grids and 3528 crossings (0.95σ from unity; factor-2 excluded at 26–35σ; seed spread 0.03467 vs `1/√(2N)` = 0.03571, so the error model checks itself), test `test_the_edge_jitter_is_the_scatter_of_a_noisy_transients_crossings_not_just_a_variance`. ⚠⚠ Phase 1 measured NO jitter on purpose — an MC sharing the analysis's noise convention cannot validate it (that pair once agreed to 0.9965 with both 2× wrong), so `σ_i² = S/(2h)` was pinned on `kT/C`: predicted 0.993789, measured 1.005742 ± 0.014891, `S/h` rejected at ~66σ. ⚠ Two withdrawals: a harness whose controls failed their own criteria (its startup-contaminated midpoint was the crossing threshold, worth ~half an apparent deficit), and a "1.2 % deficit" that was 0.95σ chased by a refinement sweep whose per-grid error equalled the effect. ⚠ LDO delay modulation + ρ_k → **A11**. Full record in `doc/pss_log_260902.md` "A8 opened" |
| A11 non-stationary delay modulation + ρ_k | ✅ **ρ_k HALF BUILT 2026-09-16** as `PAC.jitter_metrics` (σ_t, ρ_k, k-cycle `√(2(R_0−R_k))/\|s\|`, cycle-to-cycle `√(6R_0−8R_1+2R_2)/\|s\|`), test `test_the_across_period_correlation_is_the_cosine_transform_of_the_sample_series`. ⚠⚠ MY OWN ENTRY WAS WRONG: ρ_k is NOT Demir's — `sampled_noise` returns the PSD of the SAMPLE SERIES, so ρ_k is its cosine transform, `R_k = ∫S cos(2πf kT)df`. Validated three ways on a fixture built so the answer is exact (τ=RC=T → ρ_k = e^{−k}): transform **1.0005 at every lag**, MC over 1176 crossings within 1σ. ⚠ A8's τ=0.2T fixture could NOT have tested it (e⁻⁵ = 0.0067 sits under the MC's 0.03 floor). ⚠ `dc_rectangle` (the [0,fmin) rectangle) makes ρ_k fmin-independent over 100× but assumes a FLAT PSD below fmin — exact for white, wrong for 1/f — so it is OFF by default. ⚠ Guard is on the LINEARISATION (σ_t ≥ T/2 refused), not on slope: I had compared a difference (−8.86e-07) against a slope, the real slope being −176.7 = 3.58e-04 of steepest, which is also the smallest ratio a 400-point grid can produce. ⚠⚠ "non-stationary delay modulation" is **WITHDRAWN** (Andreas concurred 2026-09-17) — A8's entry already had it: the LDO is a **supply-node noise source** in an ordinary autonomous PSS, "no multirate anything". The remaining half is **ADDITIVE EDGE JITTER ON AN AUTONOMOUS CIRCUIT**: `sampled_noise` refuses autonomous PSS, so the non-accumulating part at the last buffer has no route (`c` is only the walk). ⚠ RECONNAISSANCE 2026-09-17, NOT a result: object is `σ_t² = eᵀΠP(t_j)Πᵀe/slew²` (Π = I − u_jv_jᵀ/(v_jᵀu_j), `u_j` from the rank-one `growth_samples[j]`, all `[:m,:m]` out of pair space) — NOT raw `P` (15.8 % different) and NOT `P−(t/T)G` (superseded). Fixture: vdP + 3 tanh buffers, **τ_buf = 0.5** (at 0.1 the FROZEN free-period Jacobian is singular, T/τ=66.6; `reselect` finds the same orbit to 4.44e-16 — and the library's "seed below the fundamental" message is WRONG for that failure). Analysis converges: 2A = 1.023870/1.055112/1.062936e-06 at npts 240/480/960 → ≈1.066e-06; `c_from_growth/c` = 0.998992. ✅ **BUILT 2026-09-17 as `PAC.oscillator_edge_jitter`** (σ_t, A, c, slew, k_cycle, projection_share), test `test_the_additive_edge_jitter_of_an_oscillator_is_the_projected_bounded_covariance` (10.7 s). ✅ **VALIDATED: `MC/analysis = 1.0066 ± 0.0102` (0.64σ) over 124 seed-runs on three grids.** ⚠⚠ The 4σ "defect" committed in `7abd2ea` (1.0571 ± 0.0141, "under-predicts by 5.7 %, grid-independent") was **MY HARNESS**: `σ_t = √var/slew` goes as `1/slew²`, and the MC crosses on an EULER orbit while the analysis divides by the PSS's GEAR slope. Euler's slope is low by 2.717/1.361/0.686 % at npts 240/480/960 (halving, first order) → predicted 1.0566/1.0278/1.0139, the 960 pinned before its seeds ran; dividing raw 1.0493/1.0513/1.0165 by those gives 0.9931/1.0229/1.0025. ⚠ **"Grid-independent" was asserted from TWO POINTS 0.99σ APART — two noisy points cannot tell flat from halving.** When validating against a transient, match the integrator or divide by the MC orbit's own slope. The nine-seed "1.032 ± 0.033, consistent at 0.99σ" is SUPERSEDED (centre barely moved, error bar halved). ✅ Centre shift not tail (skew 0.050, median 1.0648, 67.5 % of draws >1.0); ✅ GRID-INDEPENDENT (240: 1.0493 ± 0.0162; 480: 1.0824 ± 0.0291; 0.99σ apart) so it is the ANALYSIS's ~6 %, not MC discretisation. **STANDS**: factor-2 excluded >15σ, so the object is right and the structure is gated; the test asserts structure, NOT the ratio. ⚠ `ρ_k` excluded BY SIGN (it drives k-lag variance down; the discrepancy is up). Open candidates: projection removing the wrong component, `P` at the nearest sample vs slope at the exact instant, MC threshold from a noisy `0.5(max+min)`. ⚠⚠ TWO EARLIER NUMBERS WITHDRAWN, both from one slope error: "the analysis side converges" (committed in `39acd11`) was MY SLOPE ERROR shrinking with the grid — with the slope taken at the REQUESTED INSTANT, 2A is flat (1.055708/1.055112/1.062936e-06); and the 15.8 % projection share is a property of the SOURCE MIX (0.1611/0.0173/0.0017/0.0002 as tank noise falls 1e-6→1e-9), not of the method. ⚠⚠ The "UNEXPLAINED gear-240 anomaly" committed in `3e0ff8e` is **WITHDRAWN**: at 8 seeds the spreads are Euler 0.0804 / gear 0.1369 (ratio 1.7, not 5.4) and an 8-sample sd is itself ±0.021/±0.037, so it is **1.35σ — never an anomaly**. Controls: zero miscounted crossings either way; per-lag sd/mean flat in k. **Fourth over-read of a 3-sample spread in three days.** ✅ 16 seeds give **1.03 ± 0.03**, reproducing the committed figure. ⚠⚠ A straddling difference is NOT a slew estimator — it read −0.72/−0.08/**+2.03 %** across refinement while interpolation gives a flat 1.5275; the 1.486101 used earlier was 2.8 % wrong, inflating every `var/slew²` by 5.6 %. Andreas's to scope: a long MC campaign plus an API, not a research problem |
| B10 LSOAC (least squares, no phase condition) | ⛔ **MEASURED AND REJECTED 2026-09-17** — minimum-norm is strictly worse. At offsets 1e-1/1e-2 it returns the SAME vector as the bordered solve (1.4e-12 / 1.6e-11), because for `α ≠ 1` the system is NONSINGULAR and there is **no ambiguity to resolve** — our border fixes CONDITIONING, not arbitrariness, while LSOAC's motivation is arbitrariness. Below 1e-6 its error grows linearly in 1/offset (1.6e-4 at 1e-9) since it inherits the plain operator, whose σ_min tracks the offset while the bordered one is flat at 2.858e-02. The revival condition (minimum-norm holding at 1e-9) did not occur. ⚠ Null on ONE fixture (vdP Q=15.92, n=4); structural argument not separately demonstrated. No test added — it would gate a rejected method. Record in `doc/pss_log_260902.md` "B10 closed" |
| A6 driven oscillators and PLLs | ⚠ **"nothing built" IS LONG STALE.** ✅ **The ANALYSIS half is DONE and gated (2026-09-18)**, on a minimal loop (`VcoHdl` + multiplier PD + RC; no PFD, which would add realism rather than capability and open a `hidden_state` question). **Step 2 — closing the loop PINS the marginal phase mode**: a free-running integrator's phase row has `dx_end/dx_0 = 1.000000`, and feedback moves it to `exp(∓π·kvco·K·T)`, matched to 4e-09…4e-05 over two decades of gain on two independent knobs, with K = 0 correctly failing to converge. **Step 4 — the loop SHAPES the VCO's phase noise** against a closed form, `sf(1+2f_c/f_RC)/(4π²(f_m²+f_c²))`, to ~1e-4 over three decades with NOTHING fitted. **Step 3 — the LOCK RANGE as a LOCUS**, `−ln|λ|/(π·kvco·K·T) = s + s²(f_c/f_RC)`, 6e-07…8e-06, lock lost between 99 and 100.5 Hz against `Δf_max = 100 Hz`. **So `f_c` is pinned THREE independent ways — a Floquet multiplier, a noise corner and a detuning edge — and step 4's second-pole term, found while chasing an offset in the NOISE gate, predicts the TILT of the LOCK-RANGE curve.** ⚠⚠ Two traps the tests encode: **convergence does NOT select stability** (the natural seed lands on the SADDLE and reports success — read lock from the MULTIPLIER), and the branch is chosen by the integrator's **accumulator**, not the dependent `ph` node. ⚠ Earlier A6 work also stands: injection-locking gates, `event_grid` + `break_events`, `_fold_periodic`, and the trivial-root flag fix; saltation for a PFD/divider is **FALSIFIED twice** and state-dependent resets are **not** an event-localisation item — do not redo either. **REMAINING**: the PFD or a phase-domain detector, and realistic divider ratios, which wait on A5. See §A6 of the plan |
| A5 envelope-following | marked LAST — ⚠ **but that ordering predates A6, and the reason to do it has CHANGED**: it is now the only thing gating realistic PLL divider ratios. ⚠⚠ **SPIKED 2026-09-19 AND THE ANSWER IS NOT "BUILD IT YET"**: plain PSS reaches **N ≈ 32 in ~37 s** with cost ≈ N^1.5 over N = 4…32, but **the binding constraint is CONVERGENCE FRAGILITY, not cost** — each seeding strategy fails at exactly one N and the failure MOVES with the seeding (fixed seed fails at 32, succeeds at 64; continuation succeeds at 32, fails at 64). A basin problem, which envelope-following does not obviously fix. Attack the SEEDING first (the probe technique, B5, already exists and widens the basin), then re-make A5's case. Record in `doc/pss_log_260902.md` "The A5 spike" 🔬 **SEEDING EXPERIMENT 2026-09-19:** with the analytic locked seed, N=32 FAILS fixed and with `tstab`=2T but converges in 27.6 s with 8T; N=64 converges all three ways (297/153/290 s); the spike's continuation seed failed at 64. Three seedings, three failure sets, every converged |λ| = 0.999371 → a FRAGILE NEWTON BASIN, not a cost wall; envelope-following would not address it. Candidates (not built): a seed already periodic in the fast variable (B5 probe / a VCO phase condition), or a globalised Newton on the period map. ⚠⚠ **BOTH OF THE ABOVE ARE WITHDRAWN (2026-09-19, same day): the spike's "ends at N≈32, ~N^1.5" AND the "fragile Newton basin".** The solves were CONVERGED at iteration 2 (|F| 7e-13); the STEP test then failed on ONE unknown — the VCO output sampled at its zero crossing, tolerance `vabstol`=1e-12 — against rounding noise of 4e-12…1e-09 (the VCO's frequency node is 3.2e7, one ulp = 3.7e-09), so pass/fail was luck. **With a tolerance above the floor (`vabstol`=1e-9): 5.7 / 12.9 / 16.3 s at N = 32 / 64 / 128** (N=32 was 142 s and FAILED; 128 was never reached) — about LINEAR in N. ⚠ My first repair DECLARED SUCCESS on the stall signature and the suite refused it twice: a SINGULAR `I−M` has the same signature and no unique solution. Now COUNTED, NOT ACTED ON (`fsolve(floor_detect=True)` → `pss.step_floor`): the iteration is bit-identical and a failed solve says why, naming both causes. Plain PSS is not a wall here; A5 stays last for a measured reason. Free-period and matrix-free solves covered the same way on 2026-09-20 (counted, never acted on; unit-level tests). ⚠ `DCSweep` could not be asked for a tolerance at all until 2026-09-20 (peer finding after the default change): it now declares and forwards DC's Newton parameters. |
| ~~`coupled_method='bordered'` does not grow its step back at the default Newton tolerance~~ (found 2026-09-19 by the `vabstol` 1e-12 → 1e-6 default change) | ✅ **FIXED the same day — the `qᵀdv₀` term was COUNTED TWICE.** Eq (12) is `dh = −(f + qᵀdv₀)/denom` with the LTE residual `f` taken at the iterate BEFORE the Newton update; this code takes `err` at `x_stage1 = x + dx0`, where the update is already in it, and kept the term. `denom = err·w′/w` is tiny wherever the error is, so the spurious term decided the SIGN of `dh`: per time point the step grew 15 % on iteration 1 and SHRANK 15 % on iteration 2 (0.9775 on 8821 of 8828 points). At `vabstol`=1e-12 the loop kept iterating until `dx0` (and the term) had decayed, which is why the old default hid it. Driven RC, steps: approx 93/100, bordered 93/**8828 → 100** at 1e-12/1e-6, same error 2.5e-4. Prediction stated before the run and held. Test: a step controller's step COUNT must not depend on the Newton tolerance (fails on the pre-fix code). ⚠ HONEST CONSEQUENCE: that term was what distinguished `bordered` from `approx`; without it the two agree to every printed digit on a smooth circuit — what remains is a Newton step on the LTE equation vs the error-ratio inversion. ⚠ NOT REVISITED: the 'TESTED AND REJECTED' degree-slicing note in the same branch was tuned WITH the defect in place. |
| Unimprovable Newton step, free-period | ✅ FIXED 2026-09-16 (peer report; their "wrong Jacobian" headline WITHDRAWN — the mismatch is the documented deliberate approximation, `dF/dx_in` being singular). What survived: on a weakly limited tank (λ₂ = 0.99) the step is UPHILL against the true derivative, `fsolve` commits it anyway and said NOTHING. Now counted (`infodict['ls_unimproved']`) and named by `_diagnose_lmm_free_period_stall`; **observational only** — a case that climbs ~90 iterations then converges by 200 must not be refused. Message's "iterations do not help" narrowed (measured for gear's solved-history stall only). ⚠ `x0_unknown` and `tstab` do NOT fix that fixture; `radau` converges there in 25 iterations. Test `test_an_unimprovable_step_is_counted_and_named_instead_of_committed_in_silence` |
| B3 Aprille & Trick substitution | ✅ BUILT 2026-09-16 as `PSS.solve(phase_rule='reselect')`, **OPT-IN** (default stays `'frozen'`): the substitution alone is the bordered solve (iterates ≤ 1e-9); A&T's per-iterate re-selection is the gain — van der Pol 4× seeds 0/6 → 6/6, 10× 0/6 → 2/6, on-orbit unchanged; C3's 09-02 rejection was a pin-frame error. ⚠ Costs (found by the full suite, 6 failures on a reselect default): the grid-aligned `Idtmod` wrap stops converging, and it lands on a different PHASE of the same orbit (slow-node PPV mode content at 0.1 f0 1.64e-6 vs 2.45e-6). `_solve_twin` 399-vs-400 grid floor fixed on the way. ✅ **THAT PHASE-SPECIFIC PIN IS FIXED 2026-09-17 (`51e2d54`)** and is no longer a reason to avoid `reselect`: `mode_content` is read from the PPV AT t = 0, so it tracks the orbit phase the solve lands on, and the phase is set by the seed — a fresh sweep on ONE orbit moves the absolute plateau 17 % (2.2563…2.6410e-6) while **the T/τ RATIO is invariant to 0.0 %** (10.349…10.353), so the ratio carries the gate and the absolute became the INTERVAL [1.3e-6, 3.0e-6] (accepts all twelve recorded phases incl. reselect's 1.6437e-6 with ≥1.14× margin; still rejects a wrong power of T/τ by 7.8×). ⚠ A first rewrite as `2.2e-6 ±25 %` put its lower edge 0.08 % above the reselect reading — a knife edge, caught before it shipped: **for a phase- or seed-dependent quantity pin an INTERVAL against the measured spread, never centre±percent** |
| B4 index-2 support via independent states | awaiting a call; **priority low** — no silent wrong answer |
| ~~B10 LSOAC (minimum-norm, no phase row)~~ | ⚠ **DUPLICATE ROW, superseded — see the B10 row above**: this is the same item, and it was MEASURED AND REJECTED on 2026-09-17. Kept rather than deleted so the duplication is visible instead of silently repaired |
| ~~`orbital_spectrum` amplitude~~ | ✅ 2026-09-14 VALIDATED against pnoise on SYMMETRIC orbits (0.1 %, C=1/C=4/Q=50); on ASYMMETRIC ones `S_ph + S_orb` over-states up to 3.2× — see E6 |
| E6 Floquet decomposition over-states asymmetric orbits above f_amp | ✅ EXPLAINED 2026-09-15: PHASE half = DC PPV (fixed, frequency-aware `oscillator_spectrum`); ORBITAL half = the orbital mode's PM projection at the output, cancelled in pnoise by the phase–orbital correlation (AM share 0.307 vs flat factor 0.317 at a = 0.30; total closes ≤ 0.5 %; residual O(h)); harmonic guard ✅ BUILT. ✅ FULL CORRELATION BUILT: `PAC.modal_spectrum` (phase + orbital + full-harmonic correlation, one modal transfer) closes on pnoise (6e-4 symmetric; ≤1.4 % at 3–10 f_amp asymmetric, O(h)) |
| ~~`oscillator_covariance`'s orbital samples (the time-domain Lyapunov route) are FIRST order in the step~~ (recorded 2026-09-20, withdrawn the same day) | ✅ **WITHDRAWN — it was a ONE-NODE SHIFT in the consumers' phase-vector list, not the route.** `vs = [v0] + samples` prepends the anchor and pairs node j's covariance with node j−1's phase vector: first order in the step. It sat in three tests' route C and in shipped `oscillator_edge_jitter`; all four unshifted (`samples[j]` IS node j). Unshifted, the transverse cycle mean self-converges at 3.3e-4 / 5.1e-4 / 1.8e-4 (N=200/400/800) with the shipped `ppv` samples — the old `ppv` was fine; a `_ppv_propagate` patch that made its samples the mode's Cᵀq changed nothing there (cos 1.000000 to the mode) and cost the unit-C `v·xdot` gate, so it was REVERTED. The hostile three-way gate now pins the two routes agreeing at SECOND order: relAC 3.07e-3 / 7.90e-4 / 2.12e-4, ratios 3.9 / 3.7. Chain and rule in the log. |
| ~~`modal_spectrum`'s parts sit ~0.2 % low on a 3:1 non-uniform grid, no clean order~~ (found 2026-09-19 while converting the period harmonics) | ✅ **EXPLAINED 2026-09-20 — it was GEAR's, on any grid.** Split by ingredient: `c`, `T`, exponents and every coefficient of `p` converge; every coefficient of the adjoint mode `q` sits ~1e-3 low (q enters the parts quadratically). Gear's adjoint modes are FIRST-order consistent — the invariant `qᵀCp` drifts along the orbit on a UNIFORM grid too (9.9e-3 / 5.0e-3 / 2.5e-3 at N=200/400/800, halving); vs a uniform N=3200 reference gear's mode coefficients err 3e-2 at N=200 (rates 2.2–3) where trap's err 7e-5 (rates ~4): ~100× less accurate. The TOTAL closes on pnoise at second order (stationary in the modes); the PARTS are not. **Radau on the 3:1 grid reproduces the uniform grid to 1e-10 and keeps the phase multiplier at 1+1e-11, so `modal_spectrum` RUNS there** — ⚠ corrects b6e874a's "the phase multiplier leaves 1 at O(h²) on a non-uniform grid": true of gear/trap, not radau (the default). Refused on the way: a ±1 index shift of `q` (worse by 2.7 %), a step-size power law (α≈0.01). **BUILT 2026-09-20 (Andreas: "Build it"): gear's adjoint modes are second order on a uniform grid.** `collect` returns both the per-step transposed solve `ts[k]` (= the adjoint at node k+1) and the pair; `floquet_modes` had taken the PAIR's first block `w1 = Jfᵀt` through pinv(Cᵀ) — `a0·q(t+2h/3)`, a fractional-step stagger. Now `q_j = a0_j·ts[j−1]·exp(μt_j)` for `solved_history` only (one-step kinds keep their state-block path: radau exact, trap 2nd order). Invariant spread 9.3e-4/2.3e-4/5.6e-5 (quartering; was 9.9e-3 halving); gear parts vs trap at N=400 within 1e-3; gear's closure on pnoise 1.0004 (was 1.0145). **On a NON-uniform grid (Andreas: "Do it", 2026-09-20): `_continuous_adjoint` integrates `CᵀdQ/dt = GᵀQ` backwards with BDF2 on the REVERSE grid's own step pair, gated on STRUCTURE (`solved_history` AND non-uniform grid, never on order) — invariant quarters there too (1.8e-3/4.6e-4/1.1e-4 vs the transpose's 2.2e-2/1.1e-2/5.9e-3); two adjoints in the tree on purpose (the exact transpose stays for PPV/noise/sidebands); dense 2m×2m, ceiling 200 unknowns. Found on the way: `orbital_correlation` silently returned N²-growing garbage on a non-uniform gear grid (the phase mode at 1−5e-5 swept into the orbital sum). **Then (Andreas: gear as a first-class choice on non-uniform grids): the phase mode is now identified by its eigenvector's alignment with the orbit tangent (`PAC._phase_mode_split`), not a 1e-6 window — gear and trap RUN on non-uniform grids with a warning stating the measured departure; refused only when a second multiplier is comparably close to the circle. Gear's modal total converges to radau's on the 3:1 grid at second order (7.6e-2 / 1.76e-2 / 4.2e-3), the PPV's diffusion constant too (2.4e-2 / 5.1e-3 / 1.1e-3). Item 3 (matrix-free continuous adjoint above the 200-unknown ceiling): PROTOTYPED and validated — Arnoldi on the backward pair map with Ritz-residual certification equals the dense eigenvector to cos 1.00000000 at a basis of 3 / 16 / 24 for 2m = 4 / 32 / 124 (residuals 1e-16 … 1e-52); BUILT after 8c9a5a9: dense up to `CONTINUOUS_ADJOINT_DENSE_M`=8 unknowns, Arnoldi (Ritz-certified, per-node `lu_factor`) above, refusal by name if a mode never certifies; the 200-unknown refusal is gone. All three items done; gear runs on non-uniform grids at any size.** Test pins radau exact, gear/trap refused on 3:1, gear's drift QUARTERING on uniform. |
| ~~`break_events` is OFF for gear by default (`companion_reach()==1`), and a multistep restart after a landed edge may drop order~~ (item 1 of gear-first-class, 2026-09-20) | ✅ **MEASURED AND FIXED 2026-09-20.** No order is lost by gear across a landed edge (second order on the ladder, 3.9x/doubling) and no restart is needed. The gear-off default rested on ONE step count: on a ladder gear+events beats uniform at 5 of 6 N — now ON for every method. **The real defect was `event_grid`'s, shared by every method:** it snapped the second of two events ONTO the first, collapsing a tr=0 pulse's 1e-18 ramp to one node on the post-jump side → EVERY method first order at a true jump (gear 1.3e-3 / trap 9.5e-4 / radau 3.0e-4 at 1600, halving). Fixed (never snap an event onto an event): radau exact, trap 2nd order, gear 2.3e-3 → 3.4e-6 — BDF2 with ω→∞ over a consistent tiny step IS the trapezoidal rule. The ratio warning now counts only REPEATED up-steps (`RATIO_ISOLATION`). Gear's remaining cost on hard corners is a CONSTANT (30× trap's at the one step after a corner), recorded in `_resolve_break_events`. Open: free-period on non-uniform grids, index-2, the `method` row. |
| ~~Gear's free-period solve on a non-uniform grid unmeasured~~ (item 2 of gear-first-class, 2026-09-20) | ✅ **MEASURED 2026-09-20, nothing to fix for gear:** second order on uniform, smooth and 2:1-alternating grids (the 2:1 constant = uniform's), FIRST order beyond the zero-stability bound (3:1: ratios 1.84–1.98, warned). **Found instead: the DEFAULT method's period column was wrong wherever an autonomous circuit has a DC source** — `_traverse_full`/`_traverse_dirk` built the stage derivative as `−i(Y)` ("no source term"); an autonomous circuit's `u` is CONSTANT, not zero, so source-pinned rows read `−u/T` and radau FAILED the free-period solve on the tree's own phase fixture at every N ("seed below the fundamental" — mislabelled; the stage matrix went singular after a 2.8e6 Newton step). Fixed (`_k_at` = `−(i + u)`): FD agreement 1e-10, radau order 5 on the phase fixture. The suite's radau free-period test had a source-free fixture where the two expressions coincide. Open: index-2, the `method` row. |
| ~~Gear at index 2 unmeasured (`_ppv_propagate`'s algebraic rows "can be first order"; the continuous adjoint on a singular C)~~ (item 3 of gear-first-class, 2026-09-20) | ✅ **MEASURED AND FIXED 2026-09-20.** Forward: order 2 in BOTH subspaces on uniform and smooth grids (no BDF index-2 reduction). Autonomous index-2 fixture (van der Pol + DC source in a C–V loop): period and mode invariant second order on both grids; **`c` was FIRST order on the non-uniform grid only**, from TWO stacked causes — the index-2 fallback's PPV samples (fixed: the continuous adjoint's phase mode as `Cᵀq`, normalised `q₀ᵀC₀ẋ₀ = 1`, under the fallback's own condition) and PAC's LEFT-RECTANGLE period weights at four sites (fixed: `PAC._period_weights`, periodic trapezoid, bit-identical on uniform grids). Smooth-grid `c` −2.4e-3 / −1.2e-3 → +4.0e-5 / +3.6e-5 at 200 / 400. **Recorded, not changed: the one-step kinds' `factored_period_full/_dirk` and the TR-BDF2 twin replay on a uniform `linspace` grid whatever the solve used** — radau/trap's PPV, modes and noise on a "non-uniform grid" are uniform-grid replays (why radau read as exact on 3:1); gear is now the only method whose oscillator surfaces are on the caller's grid. Open: a noise source on an index-2 constraint gives `c = 0` silently (index-2 noise inputs); item 4. |
| ~~The one-step kinds' factored periods and the TR-BDF2 twin replay on a uniform `linspace` grid whatever grid the solve used~~ (found in item 3, 2026-09-20) | ✅ **FIXED 2026-09-20 (Andreas: "start the honour the caller's grid").** `factored_period()` hands the solved fractions to `factored_period_full/_dirk/_glm` (`grid=`; `_replay_grid`), so radau/trbdf2/esdirk43 and the twin replay on the grid the solve was on; uniform solves bit-identical, bare-count calls uniform as before. **It was live on the default method:** with events landed (7735444) `carrier_phasor` under radau was a plain mean over non-uniform nodes — 5.4 % / 2.8 % / 0.70 % off the exact harmonic at N = 100 / 200 / 400 (first order) → 3.8e-4 / 9.4e-5 / 2.3e-5 (second order). Earlier "radau exact on the 3:1 grid" readings were of uniform replays; re-measure on genuinely non-uniform grids when convenient. |
| ~~Radau's PPV / modes / noise on genuinely non-uniform grids unmeasured (earlier 1e-10 "grid-independence" was of uniform replays)~~ (2026-09-21) | ✅ **MEASURED 2026-09-21.** On alternating grids radau is order 5 for the period and `c` (3:1 at 1e-10 .. 1e-12), multipliers at 1 to 2e-11, invariant 1e-15 — the old claims survive on a real replay. **On a SMOOTHLY varying grid every method's `c` is capped at second order by the period quadrature** (radau +1.3e-5 / 3.3e-6 / 8.4e-7 at 200 / 400 / 800, the trapezoid's O(h²)); the alternating grids escape because two interleaved uniform sums are spectral. Open: a higher-order non-uniform period quadrature (lifts `c` for every method on derived grids). |
| ~~A higher-order non-uniform period quadrature (the trapezoid caps every method's `c` at second order on smoothly varying grids)~~ (found in (c), 2026-09-21) | ✅ **BUILT 2026-09-21.** `periodic_spline_weights` — the periodic cubic spline's integral as a weight vector (sparse cyclic-tridiagonal, O(n)) — used by `_period_quadrature`, `PAC._period_weights` and `fpss` on non-uniform EVENT-FREE grids; radau's `c` on a 1 + 0.5 sin grid 1.3e-5 → 2.4e-10 at N = 200 (reference floor from 400). Landed edges keep the trapezoid (a spline rings through a kink — pinned); uniform grids untouched (weights equal the trapezoid to 1.7e-18, path not taken). The `fpss`/`carrier_phasor` 1e-12 identity caught a split rule on the way — one rule, one vector, every consumer. |
| ~~Gear on ADAPTIVE grids (`lte_grid` / `refine_grid`) unmeasured~~ (2026-09-21) | ✅ **MEASURED 2026-09-21, two `lte_grid` defects FIXED, one silent error now WARNED.** `lte_grid`'s window ended at the `tend`-landing step (truncated to a tenth) → a 10.5× growth across the seam → two Euler steps → gear's transposed replay REFUSED the grid gear had produced; now the window ends at the last natural step, the cut moves only when the seam is outside the bound (by the least — an unconditional rotation broke the μ = 4 test's basin), `_period_grid`'s first-step subdivision is a doubling ramp (one `fr.min()` step before a coarse one was a 138× growth), and Transient's shrink-to-Euler heuristic is off in PSS-driven transients (`Gear2Integrator.shrink_guard`). The period passed must be within ~1 % (fractions of it; 16 % off fails every method's per-step Newton). At μ = 10 gear's period on its 195-point adaptive grid is second order and equals a 780-point uniform grid; **its `c` was 52 % high with `converged=True` and silence** — the period map's unit multiplier 0.10 off the circle, which `ppv` (solved at exactly 1) cannot see: **`ppv` now warns from `spectral_radius`** (`PPV_UNIT_MULTIPLIER_WARN`); `c` tracks the departure at 5–12×. **trbdf2 uses gear's own grid 8× better** (−66 vs −519 ppm). At μ ≥ 20 a reltol 1e-5 adaptive grid is too coarse for any PSS (the μ = 100 benchmark used 1e-7). Open: the 'closing' period-column convention (B7c) for the free period on adaptive grids; `lte_grid` measuring its own period. |
| ~~The 'closing' period-column convention (B7c: 46× closer on adaptive grids) exists on the plain path only; default "unchanged pending the rest of B7c"~~ (2026-09-21) | ✅ **BUILT 2026-09-21 on every kind** (solved-history ×3, dirk, full, glm), `period_column` Parameter, 'auto' = closing on a caller's grid for an autonomous run. It is the adaptive-grid fix: from a period seed 16 % low, where 'proportional' fails every method's per-step Newton (the fine regions slide off the edges), 'auto' converges — closing as the basin device, then ONE proportional polish on the caller's fractions at the solved period (the closing step had absorbed the whole correction: an Euler drop for gear, −1.1 % / −97 % for trbdf2 otherwise). Lands on the reference-seeded answers (gear −1464 ppm, trbdf2 −195). Bit-identical where the seed was good. |
| ~~Sampled-noise tail closure (peer: held variance short by 3.0 × (2/π)·fc/F; the 3.0 to be derived, never fitted)~~ (2026-09-21) | ✅ **DERIVED AND BUILT 2026-09-21.** The fold's kernel has coefficient ONE (symmetric |n| ≤ L covers |ν| < (L+½)f0); the peer's 3.0 is the DISCRETISATION of the covered top sidebands at ω·h ≈ 3 rad per step, held steady by a fixed M/npts — at fixed M = 100 the coefficient falls 3.03 → 1.05 as npts goes 204 → 3200. `tail=True` adds the 1/ν² extrapolation of the two outermost covered sidebands (exact for one pole, (fc/F)² in general), default off; `SAMPLED_RESOLUTION_WARN` names a top sideband above ω·h = 1. Measured: radau at 100 sidebands, 204 points: 0.9479 without (the pure tail, its covered bands accurate even at omega h = 3) -> 0.9981 with; gear at fixed M = 100: 0.847 -> 0.874 / 0.926 / 0.969 / 0.989 / 0.995 at 204 / 400 / 800 / 1600 / 3200 points -- the tail removed, the remainder the discretisation error falling with h.  ⚠ With the covered edge INSIDE the spectrum's corner (L = 6 at 204 points, F = 0.8 fc) the closure reads 0.71: the 1/nu^2 form needs F well above fc, so the warning's remedy is points first, the tail second. |
| ~~A noise source on an index-2 constraint gives `c = 0` silently~~ (found in item 3, 2026-09-20) | ✅ **WARNED 2026-09-21.** `_white_diffusion_at` names it once when `G[A,Z]` is singular at the orbit point (every kind) and `CY` has power on an algebraic row; silent for differential-row sources. |
| ~~`lte_grid` takes the period as an input and does not measure it (the grid is fractions of it; a hint 16 % off fails every method's per-step Newton)~~ (2026-09-21) | ✅ **BUILT 2026-09-21.** `lte_grid` measures the period from the run's own recurrences (rising crossings of the fastest state over the last 8 hints, interpolated; 19.1003 vs 19.0986 from a 16 %-low hint), cuts the window there, leaves it in `pss.lte_period`, warns when the hint is > 1 % off. Seeded from it: gear −17 ppm, trbdf2 +10, no second pass. The −1436 ppm recorded for gear was mostly the window's length. |
| ~~`lte_grid`/`refine_grid` step with the Transient default (`Gear2Integrator()`) whatever the PSS's method~~ (2026-09-21) | ✅ **BUILT 2026-09-21.** Both adaptive runs pass `integrator=self._integrator_for(self.par.method)`. Each method on its own grid (μ = 10): gear 195 pts −1461 ppm, second order (−361 / −92 under splitting; a window reading −22 sat on a SIGN CHANGE, −22 → +7 → +3 — the '70×' was a cancellation, and trbdf2 spreads on the same two grids too, +10 / −180: window quality, not gear), trbdf2 286 pts −8.7, radau 110 pts −0.32. The `method` table carries the guidance: trbdf2 or radau on a relaxation oscillator's adaptive grid; reltol tighter than 1e-5 above μ ~ 10; seed `solve` from `pss.lte_period`. |
| ~~PSP tests with `pytest.approx` defaults~~ | ✅ 2026-09-15 tightened, each mutation-checked; exposed 4 wrong-thing assertions (fixed) and 2 real residuals now pinned (p-ch overlap 2e-3; drain density vs ign split 2e-7). Outside PSP: ✅ 2026-09-15 `test_spicecard.py` (29, `rel=1e-12, abs=0.0`; mutation: 6/6 tiny values caught, 4/6 missed before) and `test_chained_first_class.py` (2, 1e-14 A) tightened. ⚠ `abs=` alone makes `approx` EXACT |
| ~~PSP p-channel overlap capacitance, 2e-3 on the short device~~ | ✅ NOT the model: the test's reference omits PSP's gate–bulk overlap `cgbol`; the test now subtracts `cgbov` (8.8e-7 / 3.4e-6) |
| ~~PSP sampler: flat ~+0.1 % in the sampled FLICKER part (peer external comparison, 2026-09-15)~~ | ✅ **EXPLAINED AND FIXED 2026-09-19 — it was never numerical: the coloured folds took `sqrt(PSD)` and were SIGN-BLIND.** A MOSFET's 1/f current follows sgn(Vds); Vds crosses zero twice per cycle while the sampler's switch conducts, at every clock amplitude; a coloured source is correlated across the period, so `R(t,t') = k(t)k(t')R_n` keeps the sign product and `CY = k²S` has lost it. That is why the residual survived every grid, sideband, tolerance and frequency-axis knob (A, B, C, step 0 — all refuted or unbound over 2026-09-17/18; history in `doc/pss_log_260902.md`). **Mechanism confirmed by the peer on a route independent of our fold** (periodic steady state of the linearised orbit under signed vs unsigned modulation: predicted 1.0009 against 1.0011 measured at the fixture, and the reference's converged NOTCH at clock amplitude 0.20). **Built:** `Element.noise_amplitudes` (hdl.py — the signed injection columns of the coloured sources; emitted ONLY where the model states a sign, i.e. never for a source whose power follows x under a constant scale factor), `PAC._signed_amplitudes` (used only if `W Wᴴ` rebuilds the fitted component), consumed by pnoise's cyclostationary fold and `sampled_noise`; PSP now writes `sgn * flicker_noise(n_sfl)` (CY unchanged). Gates: reference-free identity 1.000000000 where sqrt(PSD) reads 811×; PSP sweep blind/signed = 1.0007 / 1.0069 / 1.0602 / 1.42 / 3.58 / **401** / 5.16 / 2.18 against the peer's old-fold/reference 1.0011 / 1.0073 / 1.0582 / 1.39 / 3.24 / **567** / 5.87 / 2.27. ✅ **CLOSED BY THE PEER'S RUN ON `986eac8`**: ours/reference = 1.00024 / 1.00017 / 1.00013 / 1.00018 / 1.00049 / **1.036 at the notch (was 567)** / 0.99984 / 1.00018 over clock amplitude 0.75…0.10; all three pre-stated predictions held; their notch gate PASSES (and XFAILed the two earlier commits, each on the branch written for it). ⚠ `accf5a4` alone is NOT the fix: it silently fell back to the old fold at 0.375 and 0.25 (bit-identical to the old commit) — `986eac8` repairs that and makes any such fallback WARN. ⚠ Residual at 0.25 is ~2× the others (beside the notch; plausible, not explained). ✅ The six SPICE-style `kf·|I|^af` sites (diode, GP base, EKV, MOS1/3, Statz) converted the same day to `(I/|I|)·flicker_noise(...)`, CY unchanged, gated per device at element level, and **MOS level 1 at CIRCUIT level against an exact answer** (af=2 ⇒ amplitude = κ·ids(t) ⇒ a constant 1/f source × the SENSED drain current is the same physics): signed/exact = 1.00003…1.00011, sign-blind/exact = 13…58 over the clock amplitude; liveness asserted (current reverses, flicker 44 % of the total). The other five share the construction but are not each measured on a circuit. ⚠ Also found: `_uniform_exponent` let a 1e-12-of-scale entry vote and silently dropped one amplitude of the sweep into the per-band route (a bit-identical pair in a sweep is a FALLBACK, not agreement). |
| ~~PSP flicker keeps a (0.2 mV)² floor at Vds = 0 (peer)~~ | ✅ FIXED 2026-09-15: core on `vdsa = sqrt(Vds²+(2e-4)²)` made Sfl ∝ vdsa²; flicker now × `sgn²` → exactly 0 at the origin, ∝ Vds² (+16 % at 0.5 mV before). Gate-tunnelling, junction and avalanche SHOT noise then ADDED the same day (`module:1886-1906, 1951-1954`; identities 2q|I| against the element's own currents); the edge transistor stays out (SWEDGE = 0, no edge path) |
| ~~PSP `rg` uses drawn W·L~~ | ✅ fixed: `geometry()` returns `Lf`/`Wf`, `rg` divides by them (PMOS was +9.4 % / +55 %) |
| ~~`NPortS.Z` float cast; `.Z`/`.Y` fail on sweeps~~ | ✅ fixed (5103× on an RC two-port at 1e7 Hz before); swept S converted per frequency |
| ~~NPortY/Z/A `.S`/`.CS` renormalise to 50 Ω (dead `z0` property argument)~~ | ✅ fixed 2026-09-14 (Andreas confirmed the peer's spec): conversions carry the source's `z0`, `.S`/`.CS` at the stored `z0`, direct construction 50 Ω, explicit `to_s(z0)` on every class (1.18 off in S at z0 = 1 on HEAD) |
| Time-sampled noise with coloured sources (peer request) | ✅ BUILT 2026-09-15 first cut: `PAC.sampled_noise` / `sampled_variance` (sample-series PSD on (0, f0/2], one seeded adjoint per (t0, f) covers all sidebands; white held variance = `covariance` to 1e-6); stage methods (radau, trbdf2) native since the same day (CY at stage states; held 1.000000 radau); per-frequency time-average identity pinned (mean over t0 = fold of pnoise over every band to Nyquist); GLM refused (tested since the cleanup batch); Sepke et al. 2009 eq. (33) reset-RC case validated (white and 1/f, 1e-5 at 1000 points, 1.2e-3 at 400; stationary and one-step-shifted references fail) |
| ~~`pnoise(cyclostationary=True)` + 1/f ON a clock harmonic (peer)~~ | ✅ FIXED 2026-09-15: absurd finite value (6.3e-2 vs 9.2e-15, band 1.5e-11 Hz off DC by 1/T rounding) or ZeroDivisionError at exact k/T → refuses naming the harmonic; also the coloured ratio stop was dead (NaN at negative band w), fixed |
| ~~`covariance` under radau/trbdf2 at a switching edge~~ | ✅ FIXED 2026-09-15: injection at the STAGES (`_stage_injection`, CY at stage states) — radau held 0.124 → 1.3e-7 off at 400 points, trbdf2 0.130 → 2.5e-3 (second order); ESDIRK43 (non-positive weights) since the same day takes the Van Loan injection at the stage states with positive trapezoid weights over the abscissae — held 0.124 → 3.7e-4 off at 400 points, constant-operating-point behaviour unchanged |
| ~~`pnoise(cyclostationary=True)` non-additive for independent sources~~ | ✅ FIXED 2026-09-15: joint sqrt of summed CY → one root per element × {white, power-law}; was +7.3 % of the total (switch + 1/f at one node); MOS test pin 0.319 → 0.3055 |
| ~~Stale `nport.py` doctests (4)~~ | ✅ fixed 2026-09-15: `symbols('a b c d')`, today's Matrix/NumPy reprs, a rounded float view for the object-dtype `NPortS.A`; all 4 pass (suite still does not collect doctests) |
| ~~NPort swept CS "not renormalised per frequency"~~ | ✅ WITHDRAWN 2026-09-15: an untested negative claim — the pinning test passed on the unmodified code (swept `to_s(50)` = scalar per frequency to 1e-12); no change made, test kept |
| E1 GLM PCNR stage path | ✅ BUILT 2026-09-16: `pcnr=True` under a GLM method had been doing DEVICE limiting while reporting `pcnr_status='used'` (3 solves — the startup's Radau substeps — against 39987 `Diode.limit` calls). Each GLM stage now goes through `_rk_stage_pcnr` with the DIRK path's per-stage fallback; device evals 0.34–0.38× of the limiting run. ⚠ Agreement with device limiting is a TOLERANCE statement on every sequential stage path (reltol 1e-9, 1001-step fixed grid: glm3 6.2e-9, esdirk43 1.6e-8, trbdf2 8.4e-11; coupled radau 0.0), and GLM4 floors at ~1e-8 on its own (badly scaled tableau). PCNR scoreboard 4-for-4 |
| E2 GLM on the JAX backend | unbuilt — ⚠⚠ and the premise was WRONG (checked 2026-09-16): NO stage or multivalue method is on that backend. `jaxtransient.py` mentions radau/trbdf2/esdirk zero times and refuses anything but `'gear'`/`'euler'`/`'trap'` in three places; the parity ledger's P6 row agrees. GLM would be the FIRST multi-stage method there, so this is a large item, not the last family to arrive. Docstring's "CPU-only: nothing" corrected. ✅ SPIKE DONE 2026-09-16 (`benchmarks/stage_method_batched.py`, owner's choice of the three options): TR-BDF2 vs gear-2 on identical machinery, 128 lanes, CUDA — gear-2 1.7–1.8× cheaper per step, TR-BDF2's error constant ~10× smaller, so **TR-BDF2 ~1.4–1.6× faster at equal accuracy**; march validated against the CPU on an identical grid to 1.9e-12. ⚠ RECOMMENDATION: do NOT port GLM on this (4-unknown fixture, fixed-step, no LTE/rejection/breakpoints/rescue); if any stage method goes there it is TR-BDF2. ✅ **THE LARGER-`m` MEASUREMENT IS DONE 2026-09-17 (`3118a04`) AND BOTH HALVES SURVIVE**: at m = 4 / 13 / 28 (grids 400/800/1600 vs a 12800-step reference, every row `nonconv 0`) gear-2 stays cheaper per step (1.97 / 1.88 / 1.77, and 1.44–1.97 over all nine cells, never approaching 1) and **TR-BDF2 still wins at equal accuracy: 1.48 / 1.57 / 1.66**. The error-constant ratio is the stable part (8.55–9.32 everywhere), so "gear needs ~3× the steps" is a property of the METHODS, not the fixture size. ⚠ NO TREND IN m IS CLAIMED — the per-step column is monotone at 1600 steps but not at 800, and every cell is a single timing whose grid-to-grid scatter equals the m-to-m variation. ⚠⚠ Two fixture defects of mine had to be fixed first, both caught by a guard rather than by an implausible number: the chain's far sections were DEAD at the original drive (nodes at 5e-11 down to 1e-161 V, so a first "m=28" was nine live sections plus sixteen unknowns doing linear-algebra work only — the tell was error columns IDENTICAL TO EVERY DIGIT at two sizes), and **`params_tree` is keyed by CLASS, not instance name**, replacing the whole group's params, so it needs one column per element. So the recommendation's PREMISE changed (it is no longer "a four-unknown artefact") while the recommendation itself did not: 28 is still not the "hundreds" a batched backend exists for, and fixed-step is still fixed-step |
| ~~E3 `trap` non-monotone near Q = 100~~ | ✅ RESOLVED 2026-09-14 — a DEFECT (twin's `c` over trap's own period), fixed; withdraws the "trap cannot be estimated" half of the radau-default argument |
| E4 Gear-2/trap +3.0 % on a coarse fixed grid | ✅ EXPLAINED 2026-09-15 (decision unchanged, not carved out): two knee solves — no t = 0 predictor node (start reset by design) gives a linear seed 7.5 VT past the knee, then the clamp edge seeds below the root; a scratch t = 0 node → +0.7 %. Consistent-start t = 0 node = owner's call |
| ~~E5 ~2 % pnoise residual~~ | ✅ CLOSED 2026-09-14: NOT SIGNIFICANT — at 36/36 seeds the MC slow/core ratio moved −1.5 %/−1.7 %; residual +1.07 ± 1.04 % demod (1.0σ), −0.23 ± 1.01 % crossing (0.2σ). It was core-fixture seed scatter (16 seeds); the "constant" level fit is consistent with no residual at all |

⚠ **E1–E5 USED TO LIVE ONLY IN THIS FILE.**  They were real, repeatedly
recommended, and absent from the roadmap — so a reader of the plan could not
see them.  They are now **§E of the plan**, which is the authoritative list;
this table is a summary of it, not the original.  A list that exists only in a
handover is a list that dies with the handover.

⛔ Wright §3.11 is PARKED by owner decision (2026-09-10) and IS recorded in the
plan; `Ag` standing `gmin` is ⛔ DECIDED and in §A.  Neither is missing.

### What a fresh session most needs to know

⚠⚠ **`cir.G(x)` IS NOT A PURE FUNCTION OF `x`.**  `Diode` linearises around a
stored `_vlim`.  DC, AC and transient are safe because a CONVERGED solve leaves
the limiting state consistent with its own answer — but any NEW site that
evaluates `G` outside a converged solve inherits the bug.  The defence is
per-call-site by design (declaring `hidden_state` would make PSS refuse every
diode circuit).  `test_no_analysis_answer_depends_on_a_stale_limiting_state` is
the net.

⚠⚠ **A NEGATIVE CLAIM CANNOT BE MADE PROVISIONALLY.**  A test asserting an
ABSENCE passes vacuously on anything a classifier does not recognise.  Presence
claims survive that ignorance; absence claims do not — so after widening the
range of any quantity, the audit is FINITE: find the negatives.  Both live
instances are now guarded (`topological_index`'s index-0 rung reads the absence
off the MNA dimension; `Circuit.hidden_state` is gated by an allow-list with a
reason per entry).

⚠ **A diagnostic that changes the simulation is a defect** — twice now
(`branch_check` left `_vlim` at a speculative solve's value; the algebraic-block
probe needed a re-sync).  Both bracket their writes with snapshot/restore.

⚠ `copy.deepcopy` of a circuit FAILS (`cannot pickle 'module' object`, the
toolkit reference).  Do not retry it.

The measurement discipline this campaign accumulated lives in
`doc/pss_roadmap_260902.md` §D ("How these items keep failing — the shapes
worth checking for"), and it is worth reading before starting anything new.

---

# Handover — the leapfrog/symbolic work, 2026-07-30

**Read this first if you are picking the symbolic work back up.** The transient work
(`doc/transient_work_plan.md`) is expected to run first and to span several sessions, so
this document exists to make the leapfrog thread resumable after it, not to be read during
it.

Everything below is committed. ⚠ **The "unpushed backlog" this paragraph used to describe
was gone as of 2026-09-01**, when `cna-jax-vectorization` was pushed to its tip. ⚠⚠ **THAT IS NO
LONGER TRUE — 187 commits were unpushed as of 2026-09-07**, so do not read the sentence above as a
current statement about the remote.  ⚠ AND NEITHER IS THAT ONE: as of **2026-09-11 the branch is
pushed to `adf0fbc` with nothing unpushed** — see CURRENT STATE at the top of this file.  The
lesson the two stale sentences make is the point: **no count written in a document is a statement
about the remote**, and each correction has itself gone stale within days. Merging to master remains the repo owner's call and has NOT been done. Use
`git log --oneline origin/cna-jax-vectorization..HEAD | wc -l` to check the current state
rather than trusting any count written here.

---

## 1. The one-paragraph summary

The goal is full symbolic extraction of a 5th-order leapfrog filter built from five µA741
amplifiers, then approximation to expressions a designer can read. Numerics may *rank*
terms; they must not replace them. Chasing a transient cross-check for that work uncovered
that the benchmark fixture was an **unstable circuit**, then that it was an **uncompensated
one**, and separately that the transient engine had four defects. The fixture is now
repaired and compensated; the transient engine is partly repaired and fully reviewed. The third
regeneration of the tables is **done**; what remains on the symbolic side is **only T4**,
the transient-vs-perturbation comparison that started all of it.

⚠ **T4 IS RUN AND ANSWERED AT BOTH PROBES** (2026-08-31, see 4.3): the transient and
the perturbation series agree to **+0.01%** at `s0_e1` and **−0.48%** at the output —
the comparison the whole thread was built to make, from two solvers sharing no code on
the measurement side. It was unblocked by gate 0.2a (which refuted its own premise),
made affordable by running single-threaded (14–20x) and by the `x0` seeding, and
resolved by getting BOTH the tolerance and the settle right at once.

---

## 2. State of the fixture

`leapfrog_5th_order` has changed **twice** today. Any number in any document must be read
against which version produced it.

| version | dim | nnz | Q | tau | passband peak | commit |
|---|---|---|---|---|---|---|
| original (UNSTABLE — 2 RHP poles) | 127 | 536 | — | — | — | before `ff5c6e6` |
| repaired topology | 127 | 536 | 16.76 | 208 us | +8.79 dB | `ff5c6e6` |
| **+ GBW compensation (current)** | **136** | **567** | **5.93** | **73.5 us** | **+0.000 dB** | `95545e5` |

- The instability was a sign error: backward coupling entered the same summing node as the
  forward coupling, so each stage integrated the *sum* of its neighbours where a ladder
  integrates their *difference*. Confirmed on four independent routes.
- The compensation is a 350 ohm resistor in series with each integrating capacitor
  (`_LEAPFROG_RC`), correcting finite-GBW Q enhancement. **`rc=0` recovers the
  uncompensated circuit exactly** — that is how the +9 nodes are attributable.
- `build_leapfrog_network` now builds the **entire** circuit including the source,
  dispatching on toolkit (symbolic -> `VS`; numeric -> `VSin` per tone). Nothing is
  replicated anywhere. That was itself a fix: the IM3 harness had hand-carried a copy of
  the topology and kept the *broken* wiring after the fixture was repaired.

---

## 3. What is finished

- **Transient tolerance fix** (`e37ddad`): the step controller was handed the
  residual-flavoured tolerance vector where it needed the solution-flavoured one, and
  `vabstol` defaulted to 1e-12 V against a commercial simulator's 1 uV. 19x fewer steps. **See the open
  issue in section 5 — this fix has a known side effect.**
- **Gear2 LTE repair** (`doc/transient_repair_plan.md`, stages 1-5, all gates recorded):
  the `'classic'` estimate computed `q''*h^3` where BDF-2 needs `q'''*h^2`. Independently
  re-verified by the later review at est/true = 1.000282 against the 2/9 constant.
- **T2/T3 first redo** (`2d109c5`): every leapfrog experiment re-run against the repaired
  127-unknown fixture, all tables regenerated.
- **Four-lens transient review** (`7d10116`): `doc/transient_review.md` plus evidence in
  `benchmarks/transient_review/`.
- **The transient work plan** (`doc/transient_work_plan.md`) — the next thread.

---

## 4. What is IN FLIGHT and must be finished

### 4.1 T2b — the third regeneration, COMPLETE

The fixture changed again (compensation), so every leapfrog number needed regenerating a
**third** time. **All seven scripts have now run; logs are in `doc/t2b_logs/`.**

| script | result on the 136-unknown fixture |
|---|---|
| `cancellation_leapfrog` | kappa **1.153e+12** (unchanged), 7 913 groups, 14 409 600 terms, log10\|det\| −397.2 |
| `cancellation_blocks` | fails — **pre-existing**, verified identical against the pre-repair fixture |
| `cancellation_parallel` | top kappa 1.194e+12; **blocks bit-identical** (D_k 1.147e+03) |
| `cancellation_compose` | GATE 5 6/6 PASS, GATE 4 2/6 PASS — verdict unchanged |
| `nonlinear_leapfrog` | GATE 14-2 PASS rel **2.60e-13**, GATE 14-3 PASS; speedup **27-29x** (was 61-68x) |
| `transfer_function` | leapfrog tol=1e-3: 7 913 groups -> **134 520 ops = 17.0 ops/group**; uA741 control **bit-identical** |
| `order_convergence` | **`agrees from` IDENTICAL for the third time**; `v_turn` = 5.3621e-02 V |

The `transfer_function` control is worth keeping in view: its uA741 half is untouched by
every fixture change and has now re-measured **734 groups -> 50 377 operations across all
three versions**. Numbers that should not move have not moved, which is what makes the
leapfrog deltas trustworthy rather than merely different.

One recorded claim shifted: **"16 ops/group is a property of the diagram" now reads 17.0**
(it was 16 on the unstable fixture and 16.0 on the repaired one). Still structural, but no
longer identical — record it as a rescale, not a survival.

**THE RESULT WORTH HAVING.** `order_convergence` is the script whose output feeds the
pasted §10.2/10.3 tables, and it is the sharpest test of the campaign's one durable
finding. The `agrees from` column is **identical for the third time**, across three
structurally different circuits:

| amp (V) | unstable 127 | repaired 127 | **compensated 136** |
|---|---|---|---|
| 0.01 | U^5 | U^5 | **U^5** |
| 0.03 | U^5 | U^5 | **U^5** |
| 0.1 | U^5 | U^5 | **U^5** |
| 0.3 | U^7 | U^7 | **U^7** |
| 1 | U^9 | U^9 | **U^9** |
| 3 | not by U^13 | not by U^17 | **not by U^17** |

The scale moved, as it must: `v_turn` 5.3656e-02 -> **5.3621e-02 V** (`g` at `s0_e1`
4.3184e-04 -> 4.3128e-04 S), node voltages down ~19%, and `% of v_turn` now 0/0/1/3/11/32
against 0/0/1/4/13/39. The `kk` sweep holds its shape too: U^7 at 3%, none by U^17 at 34%,
and `kk = 50` still correctly flagged **UNPHYSICAL at 340% of v_turn**.

**So the claim has now survived two independent perturbations of the circuit underneath
it** — which is the strongest form the section-7 finding has taken: it measures the series
against the cubic's own validity limit, and is invariant to what the matrix does.

Run any of these as:
`cd benchmarks && PYTHONPATH=<repo>:<repo>/benchmarks MPLBACKEND=Agg python3 -u <script>.py > log 2>&1`

### 4.2 T3 — the doc rewrite for the 136-unknown fixture, COMPLETE (`83c8c22`)

All tables regenerated; gate T3-1 (stale sweep) and T3-2 (forced clean rebuild, checked in
both directions) passed. `benchmarks/order_convergence.py` also had a hardcoded
`127 unknowns` in the header of the very table that gets pasted — now `% system.dim`.

**The two traps below are kept because they will apply to the NEXT rewrite too**, not
because this one still needs doing.

**Two traps, both hit before and both recorded in `doc/leapfrog_redo_plan.md`:**
- **A normal `sphinx -b html` does NOT re-run live `exec-rst` blocks** when only the *code*
  changed. Sphinx reuses cached doctrees for unchanged sources, so the built page keeps
  stale numbers while still reporting `build succeeded, 2 warnings`. **Force it:**
  `rm -rf doc/build/doctrees` and `sphinx -E`.
- **Verify in the block's own output format.** The live block prints 3 decimals where
  `order_convergence.py` prints 4, so grepping for `7.1404e-05` returns zero hits — and so
  does the *old* value, which reads exactly like a third kind of staleness. Read the
  rendered table.

### 4.3 T4 — the IM3 transient comparison: ANSWERED at both probes (2026-08-31)

`benchmarks/nonlinear_leapfrog_sweep.py`. Two-tone IM3 (100/110 kHz, product at 90 kHz),
because **HD3 is unmeasurable here by any transient at any cost**: the third harmonic at
300 kHz is attenuated 160 000x by the filter's own stopband while the fundamental loses
only 106x. IM3 lands beside the fundamentals and is 277x larger.

⚠ **UNBLOCKED, verified 2026-08-31 — this paragraph stood for a month after its own
blocker cleared.** It read: *blocked on stage 4g; gate 0.2a re-verifies this; do not quote
a T4 number before it passes.* Gate 0.2a **passed**, and its premise was **refuted**: the
integrator choice is not contaminated, the 10x stands, and stage 4 never had to re-measure
it — the plan records that as the largest single saving stage 0 produced. Stage 4g landed
the same day (4g(a) `45f4fe0`, 4g(b) as 4i `1122c31`). T4 is the oldest open thread on the
repo and nothing is holding it back.

The original concern, kept because it explains the harness design: it sets
`TrapezoidalIntegrator()` and drives two `VSin` tones, which is the mechanism that seeds
the trapezoidal LTE estimator's parasitic `(-1)^n` mode. That mode was real and was fixed;
it just never reached this harness.

Cost, as last measured (on the *uncompensated* fixture, so re-measure): ~6.8 h per
amplitude, dominated by settling rather than stepping. The compensation cut tau 208 ->
73.5 us, which should reduce that ~2.8x. The larger lever is **seeding `x0` with the linear
two-tone steady state** — the circuit is linear apart from a cubic contributing ~1e-4, so
one AC solve removes nearly all the settling. Not implemented.

⚠ **The `x0` seeding is BUILT, 2026-08-31** — `linear_steady_state_x0` in the harness.
It computes the linear circuit's state at t=0 by AC superposition (`x_ss(0) = Im{X}`,
the conventions of `func.Sin` and `VS`'s AC stimulus agreeing because they share the
`phase` parameter). Verified exact: `e(t=0) = 0` on the real circuit, and thereafter the
seeded run deviates from the AC-predicted trajectory only by the integrator's own `h^2`
error — 2.455e-05 at the harness tolerance against 2.670e-07 at 1000x tighter, a 92x fall
for a 9.6x step reduction, which is the trapezoidal rule's order and nothing else. On an
RC with a closed-form answer it reproduced `v_out(0)` to all printed digits.

**T4 is now affordable, and mostly for a reason that has nothing to do with seeding.**
numpy was threading a 139-unknown dense solve across every core: single-threaded is
**14-20x faster** (`OMP_NUM_THREADS=1`, three interleaved pairs, and 14.2x with the cubic
live). At 0.73 s per simulated microsecond a full 5-tau run is ~14 min per amplitude
against the "hours" recorded above. The seeding then cuts the settle on top of that: the
nonlinear-node IM3 is within **0.004%** of its 5-tau value by 2 tau, so most of the
1.04 ms settle is no longer needed.

### T4 RUN, 2026-08-31 — it passes at the nonlinear node

**The two independent paths agree.** Perturbation series against transient, amplitude
1.0, seeded, 2 tau settle. The series converges to six digits by U^9 (U^13:
`IM3/f @out = 1.816614e-04`, `IM3/f @nl = 2.844758e-03`), so it is a fixed oracle with
no integrator in it. The transient against it, over a tolerance ladder:

| tolerance | IM3/f @out | vs pert | IM3/f @nl | vs pert | secs |
|---|---|---|---|---|---|
| harness (vabstol 1e-9, reltol 1e-6) | 7.157515e-04 | **+294.00%** | 2.811437e-03 | −1.17% | 234 |
| 30x (1e-11, 3e-8) | 1.760965e-04 | −3.06% | 2.840985e-03 | −0.13% | 677 |
| 300x (1e-13, 3e-9) | 1.696682e-04 | **−6.60%** | 2.843814e-03 | **−0.03%** | 1573 |

**At the nl node this is T4 answered.** Monotone convergence onto the series —
−1.17% → −0.13% → −0.03% — a different solver, a different representation and a
different failure mode agreeing to 3 parts in 10^4. That is the property the
transient-vs-perturbation gate exists to establish, and the amplifier's own
`s0_e1` is where the nonlinearity acts, so it is the node that tests the nonlinear
machinery rather than five stages of linear filtering.

⚠ **RESOLVED AT BOTH PROBES the same day — the block that stood here said the output
was "not answered" and that the remaining work was a tolerance ladder. It was half
right.**

**Converged configuration:** seeded, `reltol 3e-9` / `vabstol 1e-13` (300x the harness
defaults), settle **5 tau**, quadratic resampling. Against the series (U^13:
`out 1.816614e-04`, `nl 2.844758e-03`):

| probe | transient | vs perturbation |
|---|---|---|
| `s0_e1` (where the nonlinearity is) | 2.844924e-03 | **+0.01%** |
| output | 1.807923e-04 | **-0.48%** |

Two independent paths — a different solver, a different representation, a different
failure mode — agreeing to better than half a percent at both probes. **T4 is answered**,
at **1087 s per amplitude** single-threaded against the ~6.8 h recorded above.

⚠ **IT TOOK TWO VARIABLES, AND MOVING ONE AT A TIME READ AS NON-CONVERGENCE.** The
tolerance ladder at a fixed 2 tau settle gave +294% -> -3.06% -> **-6.60%** — it crossed
the series value and kept going, which looked like a real disagreement, or a floor that
tightening could not reach. It was neither. Holding tolerance at 300x and moving the
SETTLE instead gave **-6.60% -> -0.74%**, a factor of nine. Both knobs were wrong at
once, and **a ladder in one variable while another is wrong is not a convergence
study** — it is a slice through a surface, and it can point away from the answer.

The seed is what put the settle under suspicion, and the distinction held exactly:
`linear_steady_state_x0` removes the LINEAR settling (`e(t=0) = 0`), and what remained
was the cubic's own approach to periodic steady state — which reaches the output through
five filter stages, so the nl probe was settled at 2 tau while the output was not.

The last 0.26 percentage points came from resampling at the integrator's order
(`resample_uniform`, `c34881c`) rather than `np.interp`: -0.74% -> -0.48%, matching the
~0.2% a synthetic test predicted for this grid's step spread.

---

## 5. Open decisions — these need the maintainer, not an implementer

1. **`vabstol` serves two roles** (`doc/transient_work_plan.md` 0.3a). It is Newton's
   x-tolerance *and* the LTE tolerance. The 1e-12 -> 1e-6 fix in `e37ddad` was reasoned
   about only as a step-control knob, so it **also loosened Newton's node convergence by
   10^6, unmeasured**; and `DC.vabstol` is still 1e-12, so the operating point seeding
   every transient is solved 10^6 times tighter than any step after it. All four review
   lenses flagged it.
2. **Gate 1-3 of `transient_repair_plan.md`** is still open: its 2x threshold was written
   against a metric later shown to be degenerate, and the sharp reading gives 3.94. Needs
   a threshold justified by the 2/9-vs-1/12 error constants, or an explicit "recorded, not
   resolved".
3. **`Gear2`'s default `lte_formula`** (0.3b) — the evidence now runs against the `'ywr'`
   default that was chosen belt-and-braces.
4. **Scope of transient stage 10** (0.3c) — the missing-analyses list is a product
   decision.

---

## 6. Refuted premises — do not re-derive these

Recorded because each cost real time and each was *stated before it was measured*.

- **"The stable circuit will be cheaper to integrate."** False: 2896 steps against 571,
  5x *more*. A smooth exponential blow-up is easy to track; a lightly-damped high-Q
  resonance is not.
- **"kappa is driven by the Q resonance."** False: kappa is **1.153e+12 on both the
  Q=16.76 and Q=5.93 fixtures**, identical to four significant figures. The numerical
  review independently explains why — `log10|det|` is a *unit-scale* artefact, moving +365
  decades as h goes 1e-9 -> 1e-12, which is 3 decades x ~122 ~ `rank(C)`. kappa here is a
  property of the determinant representation, not of the circuit.
- **"Compensation will make the benchmark easier."** False: 6.4x more terms, 32% more
  groups, and the determinant underflows *further* (−397 vs −374).
- **"The transient never ran because the circuit is stiff."** False. The mean step was
  5.03 ns against a 39 ns cap — a ratio that says steps are being *rejected*, not that the
  circuit is stiff.
- **`cancellation_blocks`'s failure is caused by the fixture change.** False — it fails
  identically against the pre-repair fixture. Verified by checking out `ff5c6e6^`'s
  `benchmark_circuits.py` and re-running.

---

## 7. The one durable finding from the redo campaign

**What was refuted measured the determinant's *conditioning*; what survived measured the
*method*.** A topology change alters matrix conditioning and leaves the properties of the
series, and of the diagram representation, alone.

Survived both fixture changes unchanged: the order-convergence pattern (`agrees from` was
identical at all six amplitudes through the first redo), 16.0 ops/group, the symbolic
nonlinear analysis matching its numeric oracle, "readable H(s) over 2 decades only".
Refuted or rescaled: kappa, group counts, term counts.

Read any future claim in these documents against which of those two it depends on.

---

## 8. Environment notes that will save an hour

- `PYTHONPATH=/home/andreas/sources/pycircuit` and `MPLBACKEND=Agg` are **required** — a
  stale root-owned egg in `/usr/local/lib/python3.14/dist-packages` shadows the source.
- **Two venvs, not interchangeable.** The repo's `.venv` has sphinx but no test deps; a
  scratch venv (`python3 -m venv --system-site-packages`, then `pip install pytest pynose
  pytest-timeout`) runs the suite but has no sphinx.
- Suite: `pytest pycircuit -q -p no:cacheprovider -m "" --timeout=400 --timeout-method=signal`.
  **`-m ""` is mandatory** — 17 `slow` tests are deselected by default and they are the only
  ones comparing an analysis against a time-domain reference.
- Doc build: `cd doc && MPLBACKEND=Agg ../.venv/bin/python -m sphinx -b html -d build/doctrees src build/html`.
  Clean baseline is `build succeeded, 2 warnings`.
- **Never pipe a long-running command through `tail`/`grep`/`head`** — output is lost
  entirely if it is killed, and a 0-byte log is indistinguishable from "produced nothing".
  Redirect `> log 2>&1` and read the file. Use `python3 -u`.
- To check whether a background job is alive, match on `comm` (the executable), never on
  the command line — a `bash -c` wrapper contains the full python invocation and matches
  too.

---

## 9. Document map

| document | what it is |
|---|---|
| `doc/transient_work_plan.md` | **the next thread** — staged plan, gates declared, nothing run |
| `doc/transient_review.md` | the four-lens review and its measurements |
| `benchmarks/transient_review/` | evidence probes for the review's performance/solver numbers |
| `doc/transient_repair_plan.md` | the completed Gear2 LTE repair, gates recorded |
| `doc/transient_repair_reasoning.md` | why that repair was scoped as it was |
| `doc/leapfrog_redo_plan.md` | the fixture-change redo campaign; T0/T1/T2/T3 outcomes |
| `doc/distortion_ddd_conclusions.md` | §10.1-10.3, the nonlinear leapfrog results |
| `doc/cancellation_ranking_conclusions.md` | the cancellation/group-ranking results |
| `doc/t2b_logs/` | the five completed T2b runs on the 136-unknown fixture |
