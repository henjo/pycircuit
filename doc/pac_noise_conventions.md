# PAC noise surfaces: the conventions

What each public noise, jitter and phase-noise method of `PAC`
(`pycircuit/circuit/shooting/`) returns, on what scale, and with which
arguments -- and where the surfaces disagree with each other.  Written
2026-09-28 from a read of the code (the docstrings where they agree with
it); a renaming is not implied by anything here.

**Changed 2026-09-28** (Andreas, aligning with a commercial simulator's
pnoise): `oscillator_spectrum`'s `S_v` is ONE-SIDED (it was 0.5x; `L_dBc`
unchanged), `am_pm_noise` returns PER-SIDEBAND densities (they were the
pair's totals, 2x), and `output` takes a node NAME too.  Then the modal,
orbital and correlation spectra ONE-SIDED as well (they were 0.5x): a real
output has `S(-f) = S(f)`, so the one-sided PSD is twice the two-sided one
on BOTH sidebands, which are NOT symmetric about the carrier (measured:
the lower 2.45x the upper on an asymmetric van der Pol, as `pnoise`) and
keep their own values.  The tables below are the new conventions; the
list at the end marks what these closed.

**Changed 2026-09-29** (Andreas: match a commercial simulator where it has
a convention, else IEEE one-sided): the swept frequency follows that
simulator's `sweeptype` rule (item 4); the band arguments are named for
their band (item 6); the covariances are `m x m` whatever the method
(item 7); one jitter vocabulary (item 8); flag pairs, a missing carrier and
a swapped argument refused (items 9-11).  The tables are the code as of
then.

## Shared conventions

- **`CY` is one-sided.** An element's noise density is the one-sided PSD
  (a resistor's `4kT/R`), the scale of `analysis_ss.Noise`'s `Svnout`. The
  variance and diffusion routes use `CY/2` internally (a two-sided
  intensity), which is not visible in their results.
- **`output` is a node NAME, a REDUCED-state index, or a weight vector.**
  A name (a string, or a `Node`) is resolved to its reduced index
  (`_numerics.output_index`; the reference node is refused).  The reduced
  state is the circuit's unknown vector with the reference node's row
  removed.  A vector is taken as weights over the reduced unknowns (a
  differential output).  Names at the top level of the circuit.
- **A swept frequency follows a commercial RF simulator's `sweeptype`
  rule** (`_numerics.sweep_kind`, since 2026-09-29): `'absolute'` is the
  frequency itself, `'relative'` an offset from a harmonic of the PSS
  fundamental, and `None` -- that simulator's 'unspecified', the default --
  relative on an AUTONOMOUS PSS and absolute on a driven one (an
  oscillator's period is an output of the solve, so a frequency near a
  harmonic is only writable as an offset).  `pnoise` and `PAC.solve` take
  `sweeptype` and `relharmnum` (the harmonic, default 1); `am_pm_noise`
  and `am_pm` take `sweeptype` with their `carrier` as the harmonic.  The
  spectra take OFFSETS and the sampled family SERIES frequencies, by name.
- **A band argument is named for its band:** `colour_fmin` /
  `colour_fmax`, the band a COLOURED source is integrated over (white
  sources always over every frequency); `series_fmin` / `series_fmax`, the
  band of a SAMPLED series, which cuts white noise too; `offset_fmin` /
  `offset_fmax`, the lineshape's phase-offset band.
- **Oscillator surfaces refuse a driven circuit, and the driven-circuit
  surfaces refuse an oscillator**, each by name.
- **A coloured source rooted from its PSD is WARNED where that PSD touches
  zero along the orbit** (`_noise_components.warn_sign_blind`, one policy
  since 2026-09-29).  An element that states no signed amplitudes
  (`Element.noise_amplitudes`) is factored by `sqrt(PSD)`: the `|m|`
  process, exact while its modulation keeps its sign and wrong in either
  direction where it changes sign.  Every surface that roots a coloured
  PSD says so -- the spectra and the modal family, and since 2026-09-29
  the sample family and the covariance family, which were silent.

## The scales, side by side

| scale | methods |
|---|---|
| one-sided PSD of the output, unit^2/Hz | `pnoise`, `sampled_noise`, `oscillator_spectrum` (`S_v`: the Lorentzian times the carrier's one-sided power `2 |X|^2 = A^2/2`), `orbital_spectrum`, `modal_spectrum`, `correlation_spectrum`, `am_pm_noise` per sideband: `2 (S_am + S_pm) = pnoise(k f0 + f) + pnoise(k f0 - f)` (and `analysis_ss.Noise`) |
| `L(f)`, the single-sideband phase noise per Hz relative to the carrier (a commercial simulator's phase noise; the IEEE one-sided `S_phi` is `2 L`) | `phase_psd`, `lorentzian` (integrates to 1 over `(-inf, inf)`) |
| dBc/Hz, `10 log10(S_v / (2 |X|^2))` = `L(f)` | `oscillator_spectrum` (`L_dBc`) |
| variance, unit^2 | `sampled_variance`, `covariance`, `oscillator_covariance` (`K_orb`), `orbital_correlation` |
| seconds (jitter), s^2 (growth) | `jitter_metrics`, `event_jitter`, `oscillator_edge_jitter`; `oscillator_covariance`'s `d` |
| diffusion constant `c`, seconds | `diffusion_constant`, `frequency_aware_diffusion`, `coloured_diffusion`, `coloured_diffusion_resolved` |

So near a carrier `pnoise ~ S_v ~ S_pm` (per sideband, where AM is
small), and the modal total is `pnoise` at every offset (to its truncation).

## Driven circuits

| method | arguments, in order | returns | frequency | key defaults |
|---|---|---|---|---|
| `pnoise` | `pss, freq, output, ratio_tol=None, maxsidebands=None, modulated=False, cyclostationary=False, sweeptype=None, relharmnum=None` | `(S, sidebands_used)`; also sets `self.alias_stop`, `self.sidebands_used` | `freq`, scalar, the OUTPUT frequency by the sweep rule: itself when absolute (the default on a driven PSS), `relharmnum f0 + freq` when relative (the default on an oscillator) | `maxsidebands` -> `N//2` (an explicit count above it raises; clamped silently until 2026-09-29); `ratio_tol` -> 1e-9, stops after two quiet pairs |
| `am_pm_noise` | `pss, freq, output, harmonic=1, maxsidebands=None, modulated=False, sweeptype=None` | `(S_am, S_pm, bands)`, PER SIDEBAND (half the pair's total; it was the total until 2026-09-28) | `freq` by the sweep rule with `harmonic` the reference: an OFFSET from `harmonic f0` when relative (the default on an oscillator), the upper sideband's own frequency when absolute (the default on a driven PSS) | `maxsidebands` -> every pair within Nyquist, `N//2 - harmonic` (an explicit count above it raises) |
| `band_spread` | `pss, output, band, points=9, harmonic=1, quantity='pnoise', **kw` | `(spread, info)` | `band` in units of f0: an OFFSET from `harmonic f0` for every quantity (`**kw` takes no `sweeptype`) | `points=9` |
| `sampled_noise` | `pss, output, times, freqs, maxsidebands=None, tail=False` | array `(len(times), len(freqs))` | `freqs`: SERIES frequency, `0 < f <= f0/2` | `maxsidebands` -> `N//2 - 1` (raises above) |
| `sampled_variance` | `pss, output, times, series_fmin, series_fmax, points_per_decade=40, maxsidebands=None, tail=False` | array `(len(times),)` | `series_fmin`, `series_fmax`: the SERIES band, required; the band cuts white noise too | `ppd=40` |
| `jitter_metrics` | `pss, output, time, series_fmin, series_fmax, kmax=8, maxsidebands=None, nfreq=601, dc_rectangle=False, points_per_decade=40` | dict: `sigma_t, rho, k_cycle` (exact), `cycle_to_cycle, slew` (`edge_slope`), `R, instant` | as `sampled_variance` | `kmax=8` (>= 2) |
| `covariance` | `pss, samples=False, colour_fmin=None, colour_fmax=None, points_per_decade=40, pair=False` | `(K0, info)`: `info['samples']` at `info['times']` with `samples=True`, `info['K_coloured']` (+ `'coloured_samples'`) with a band; `m x m` whatever the method (a gear/trap pair map's `2m x 2m` with `pair=True`) | `colour_fmin`, `colour_fmax`: the SOURCE band of the COLOURED part only; `colour_fmax` -> the grid's Nyquist | `colour_fmin` required when a source is coloured |
| `event_jitter` | `pss, colour_fmin=None, colour_fmax=None, points_per_decade=40` | dict: `sigma_t` (s), `cov_fraction`, `fractions`, `nodes` | as `covariance` | refuses a solve with no landed events |

## Oscillators

| method | arguments, in order | returns | frequency | key defaults |
|---|---|---|---|---|
| `oscillator_spectrum` | `pss, offsets, output, harmonic=1, frequency_aware=True, offset_fmin=None, offset_fmax=None, all_orders=None` | `(S_v, L_dBc)`; sets `self.lineshape_info` | OFFSETS from `harmonic f0`, any sign, EVEN in them | `offset_fmin` needed only with colour; `offset_fmax` -> f0/2; `all_orders` with `frequency_aware=False` refused |
| `phase_psd` | `pss, offsets, harmonic=1, frequency_aware=True` | array, `L(f)` per Hz (the IEEE `S_phi` is `2 L`) | OFFSETS, `> 0` | refuses offsets at or below the Lorentzian corner |
| `lorentzian` (static) | `offsets, c, f0, harmonic=1` | array, 1/Hz relative to the harmonic's power | OFFSETS, any sign, even in them | harmonic 0 returns zeros |
| `modal_spectrum` | `pss, offsets, output, harmonic=1, maxharmonics=None, maxsidebands=None` | dict: `phase, orbital, correlation, total` | OFFSETS, negative = the lower sideband (not the upper's) | `maxharmonics` -> 32 capped by the grid (an explicit count above the grid raises), `maxsidebands` -> `2 maxharmonics` |
| `orbital_spectrum`, `correlation_spectrum` | `pss, offsets, output, harmonic=1, maxharmonics=None` (+ `maxsidebands` for correlation) | array | as `modal_spectrum` | as `modal_spectrum` |
| `orbital_correlation` | `pss, maxharmonics=None` | `(R, C)`: `R_yy(0)` `m x m`, `C[(l, h, j)]` | -- | refuses colour |
| `oscillator_covariance` | `pss, samples=False, colour_fmin=None, colour_fmax=None, points_per_decade=40, pair=False` | `(K_orb, info)`: the growth `info['d']`, the bounded part per node `info['samples']`; `K_orb` `m x m` whatever the method (the pair with `pair=True`) | as `covariance` (colour enters the transverse part only) | |
| `oscillator_edge_jitter` | `pss, output, time, kmax=8` | dict: `sigma_t, A` (the ONE-sided projection, the k-lag law's intercept / 2), `c, slew, k_cycle` (EXACT at every k, as `jitter_metrics`), `instant, d, projection_share` | -- | `kmax=8`; `output` an index or a name |
| `orbital_mode_weights` | `pss, nmodes=None` | `(cw, info)`: `info['modes']`, `info['K']` | -- | |
| `diffusion_constant` | `pss` | `c`, s | -- | refuses colour |
| `frequency_aware_diffusion` | `pss, offset` | `c(f)`, s | `offset` scalar, abs taken | refuses colour |
| `colour_projection` | `pss` | `(vbar, info)`: `rms, symmetry, samples, times` | -- | |
| `coloured_diffusion` | `pss, offsets` | `Gamma(f)`, the `l = 0` term ONLY | OFFSETS, any sign, even in them | do not add it to `c` (counts `l = 0` twice); use `coloured_diffusion_resolved` |
| `coloured_diffusion_resolved` | `pss, offsets, harmonics=None, frequency_aware=True` | `c(f)` over all harmonics | OFFSETS, any sign, even in them | `harmonics` -> every `l` with > 1e-14 of the PPV energy |

## Transfers (not noise)

| method | arguments | returns |
|---|---|---|
| `solve` | `pss, freqs, refnode=gnd, recycle=True, sweeptype=None, relharmnum=None` | a `CircuitResult` over the absolute OUTPUT frequencies; `freqs` the SOURCE frequencies by the sweep rule, as `pnoise`'s `freq` |
| `am_pm` | `pss, freq, output, harmonic=1, sweeptype=None` | `(m_am, m_pm)` per unit source; `freq` by the sweep rule, as `am_pm_noise`'s; refuses a missing carrier |
| `am_pm_indices` | `a, b` | `(a + conj(b), a - conj(b))` |
| `carrier_phasor` | `pss, output, harmonic=1` | the complex Fourier coefficient (`A/2` for `A cos`) |

## Where the surfaces disagree

**Can give a wrong number:**

1. ~~**Two sidedness scales under one unit.**~~ CLOSED 2026-09-28: every
   spectrum is one-sided (the orbital/modal/correlation spectra and `S_v`
   were 0.5x); `phase_psd` stays the two-sided `S_phi = L(f)` (item 2);
   only `oscillator_spectrum` returns dBc.
2. ~~**`phase_psd` is the two-sided `S_phi` (= `L(f)`),** not the IEEE
   one-sided `S_phi = 2 L(f)`; its docstring says only "rad^2/Hz".~~
   CLOSED 2026-09-29 (Andreas: match a commercial simulator, else IEEE):
   that simulator reports an oscillator's phase noise as `L(f)` alone, and
   `phase_psd` returns exactly that -- documented so; no number moved.
3. ~~**`am_pm_noise` returns pair totals** over both sidebands, not a
   per-sideband density; on an oscillator `S_pm = 4 S_v`.~~ CLOSED
   2026-09-28: per sideband, and on an oscillator `S_pm = S_v`.
4. ~~**`freq` is ABSOLUTE in `pnoise` and an OFFSET in `am_pm_noise` /
   `am_pm`** -- the same name in the same position.~~ CLOSED 2026-09-29
   (Andreas: check what a commercial simulator supports and match; PAC.solve
   too): that simulator's `sweeptype` rule, `absolute` / `relative` /
   unspecified (relative on an autonomous PSS, absolute on a driven one),
   with `relharmnum` the harmonic of a relative sweep -- on `pnoise` and
   `PAC.solve`; `am_pm_noise` and `am_pm` take `sweeptype` with `carrier`
   as the harmonic.  ⚠ The DEFAULT moved for `pnoise` / `PAC.solve` on an
   oscillator (now an offset from `f0`) and for `am_pm_noise` / `am_pm` on
   a driven circuit (now the upper sideband's frequency); the 44 test
   calls that relied on the old reading pass it explicitly.
5. ~~**`band_spread`'s `band` is an offset for three quantities and
   absolute for 'pnoise'**; its `**kw` is not forwarded to
   `oscillator_spectrum`; an unknown `quantity` falls through to
   'pnoise'.~~ CLOSED 2026-09-29: an offset from `harmonic f0` for every
   quantity (the pnoise it samples is the relative sweep), forwarded, an
   unknown name raises.
6. ~~**`fmin`/`fmax` name three different bands:** the coloured SOURCE band
   (`covariance`, `event_jitter`, `oscillator_covariance`; ignored when all
   sources are white), the SERIES band that cuts white noise too
   (`sampled_variance`, `jitter_metrics`), and the phase-offset band for
   colour (`oscillator_spectrum`).~~ CLOSED 2026-09-29 (Andreas: rename
   per meaning): `colour_fmin/colour_fmax`, `series_fmin/series_fmax`,
   `offset_fmin/offset_fmax`; the messages say the new names.
7. ~~**`covariance` is `m x m` on radau/euler/GLM and `2m x 2m` on gear/trap**
   (the pair state), so traces, eigenvalues and shapes depend on the method.~~ CLOSED 2026-09-29: `m x m` whatever the method,
   the pair on request (`pair=True`), `oscillator_covariance` too.
8. ~~**The jitter dictionaries disagree:** `sigma_t` vs `sigma`
   (`event_jitter`), `k_cycle` exact in `jitter_metrics` but a large-`k`
   upper bound in `oscillator_edge_jitter`, `slew` from a central
   difference vs a quadratic fit, `kmax` 4 vs 8.~~ CLOSED 2026-09-29:
   `sigma_t` everywhere, the bound named `k_cycle_bound`, one slope
   estimator (`edge_slope`, the quadratic fit), `kmax` 8.  Then (#17 B0/B4,
   the same day) that bound was found to drop a transverse-phase cross term
   (a committed Monte Carlo, 8.9 sigma): `oscillator_edge_jitter` returns the
   EXACT `k_cycle` and the one-sided `A` -- one meaning of `k_cycle` again.
9. ~~**Flag combinations dropped silently:**~~ CLOSED 2026-09-29 for the
   two that lost a flag (both raise now; `frequency_aware` defaults to
   True, as it behaved): `pnoise` with both `modulated`
   and `cyclostationary` (cyclostationary wins); `oscillator_spectrum` with
   `frequency_aware=False` drops `all_orders`; `all_orders=None` means
   different things for white and coloured sources; `frequency_aware`
   defaults to None in `oscillator_spectrum` and True elsewhere.
10. ~~**A missing carrier:**~~ CLOSED 2026-09-29: `am_pm_noise` refuses
    it too (`lorentzian`'s zeros at harmonic 0 are the right answer: DC does
    not broaden).  Was: `am_pm` refuses it, `am_pm_noise` proceeds
    unrotated, `lorentzian` returns zeros for harmonic 0, the oscillator
    spectra refuse.
11. ~~**Argument order:**~~ CLOSED 2026-09-29 as a silent wrong number:
    `output_index` validates an index (an integer row) and a weight vector
    (exactly the reduced width), so a swapped argument raises.  The orders
    themselves still differ: `(pss, freq/offsets, output)` for `pnoise`,
    `am_pm_noise` and the spectra, but `(pss, output, ...)` for the sampled
    family, `band_spread` and `oscillator_edge_jitter`. Because `output`
    accepts a vector, a swapped array can be taken as weights silently.

**Naming only:**

12. ~~The truncation knob has many names: `maxsidebands` (`N//2`, clamped, in
    `pnoise`; `N//2 - 1`, raising, in the sampled family), `H`,
    `sidebands` (a count in `modal_spectrum`, a list in
    `adjoint_sideband_row`), `harmonics`, `kmax`, `nmodes`.~~ CLOSED
    2026-09-29 (Andreas: "Rename + raise"): a sideband COUNT is
    `maxsidebands` everywhere (`modal_spectrum`, `correlation_spectrum`
    said `sidebands`), the PPV / mode harmonic count `maxharmonics` (was
    `H`); `sidebands` stays a SELECTION of indices (`adjoint_sideband_row`,
    `mixer_response`), `harmonics` a selection of `l`, and `kmax` /
    `nmodes` are result lengths, not truncations.  An EXPLICIT count above
    what the grid resolves raises everywhere (`pnoise`, `am_pm_noise` and
    the modal harmonic count clamped silently); the defaults are as many as
    the grid resolves.
13. ~~The harmonic is `harmonic` in the spectra and `carrier` in `am_pm`,
    `am_pm_noise`, `carrier_phasor`.~~ CLOSED 2026-09-29: `harmonic`
    everywhere.
14. ~~Frequency names: `freq` (scalar), `freqs` (a series frequency or an
    offset), `offsets`, `offset`, `band` (units of f0), `time` vs `times`;
    offset signs differ (`phase_psd` refuses <= 0, `oscillator_spectrum`
    accepts any, modal/orbital read negative as the lower sideband).~~
    CLOSED 2026-09-29: an offset is `offsets` everywhere (the coloured
    diffusion said `freqs`); `freqs` is left to input and series
    frequencies, `time` / `times` to one instant or many.  Each offset
    surface states its sign convention, measured at +-o on an asymmetric
    orbit: `phase_psd` refuses <= 0; `oscillator_spectrum`, `lorentzian`
    and the coloured diffusion are EVEN in the offset;
    `frequency_aware_diffusion` takes its magnitude; the modal family reads
    a negative offset as the LOWER sideband (1.02-1.06x the upper there,
    and the correlation term changes sign).
15. ~~Returns: tuples, dicts, arrays and floats; `covariance` changes arity
    with `samples` while `oscillator_covariance` puts samples in `info`;
    diagnostics go to instance attributes (`alias_stop`, `sidebands_used`,
    `sampled_instants`, `lineshape_info`, ...).~~ CLOSED 2026-09-29 for the
    covariance family (Andreas: "Redesign #15 too"): `(value, info)`, the
    house pattern of `colour_projection`, `band_spread`, `PSS.ppv` --
    `covariance -> (K0, info)` whatever `samples` is,
    `oscillator_covariance -> (K_orb, info)` with `d` in `info['d']` and
    `'orbital_samples'` renamed `'samples'`, `orbital_mode_weights ->
    (cw, info)`; the jitter surfaces keep their metric dicts (as
    `jitter_metrics`).  The instance-attribute diagnostics belong to other
    surfaces (`pnoise`, `solve`, the lineshape, the sampled family) and are
    left as they are.
16. ~~The same refusal raises different exceptions: a missing lower band
    edge is a `NotImplementedError` in `covariance` (`colour_fmin`) and
    `oscillator_spectrum` (`offset_fmin`) but a `TypeError` in the sampled
    family (`series_fmin`, a required argument).~~ CLOSED 2026-09-29: a
    `TypeError` everywhere -- a required argument missing, as Python's own
    in the sampled family (`covariance`, `event_jitter`,
    `oscillator_covariance`, `oscillator_spectrum`).  A surface that cannot
    take colour at all still raises `NotImplementedError`.
17. ~~Several oscillator surfaces cannot take a coloured source at all
    (`oscillator_edge_jitter`, `orbital_mode_weights` have no
    `colour_fmin` to pass; `orbital_*` refuse colour; `modal_spectrum` accepts it).~~
    CLOSED 2026-09-29 (plan approved): `oscillator_edge_jitter` and
    `orbital_mode_weights` take `colour_fmin` / `colour_fmax` /
    `points_per_decade`.  The edge jitter's coloured PHASE is the colour
    fold's increment at the requested instant (a coloured source has memory:
    the stationary structure function was 33 % off at a van der Pol edge),
    gated on the edge time against the source's exact white realisation
    (1.4e-3); the mode weights' coloured part is built in the map's own
    space from the bordered solution, so its phase row vanishes (1.1e-12).
    On the way (B0/B4): the WHITE k-cycle law dropped a transverse-phase
    cross term -- `k_cycle` is now the exact law (a committed Monte Carlo,
    8.9 sigma) and `A` the one-sided projection.  `orbital_correlation` /
    `orbital_spectrum` still refuse colour (`modal_spectrum` covers it).
18. ~~Stale text: `colour_projection`'s docstring promises "the ratio
    |mean|/rms" under a key that is `symmetry`; three comments still
    describe a Gear-2 fallback for TR-BDF2 that the accuracy host no longer
    uses (it uses the monodromy twin).~~ CLOSED 2026-09-29: the key named;
    the three comments (and `_event_closure`'s host sentence) say what the
    hosts do.
