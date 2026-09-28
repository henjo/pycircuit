# PAC noise surfaces: the conventions

What each public noise, jitter and phase-noise method of `PAC`
(`pycircuit/circuit/shooting/`) returns, on what scale, and with which
arguments -- and where the surfaces disagree with each other.  Written
2026-09-28 from a read of the code (the docstrings where they agree with
it); a renaming is not implied by anything here.

**Changed 2026-09-28** (Andreas, aligning with a commercial simulator's
pnoise): `oscillator_spectrum`'s `S_v` is ONE-SIDED (it was 0.5x; `L_dBc`
unchanged), `am_pm_noise` returns PER-SIDEBAND densities (they were the
pair's totals, 2x), and `output` takes a node NAME too.  The tables below
are the new conventions; the list at the end marks what these closed.

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
- **Frequency arguments are one of three kinds** -- an ABSOLUTE output
  frequency, an OFFSET from a carrier harmonic, or a SERIES frequency of a
  sampled sequence -- and the argument names do not tell which (see the
  inconsistencies below).
- **Oscillator surfaces refuse a driven circuit, and the driven-circuit
  surfaces refuse an oscillator**, each by name.

## The scales, side by side

| scale | methods |
|---|---|
| one-sided PSD of the output, unit^2/Hz | `pnoise`, `sampled_noise`, `oscillator_spectrum` (`S_v`: the Lorentzian times the carrier's one-sided power `2 |X|^2 = A^2/2`), `am_pm_noise` per sideband: `2 (S_am + S_pm) = pnoise(k f0 + f) + pnoise(k f0 - f)` (and `analysis_ss.Noise`) |
| **0.5 x** a one-sided PSD (the carrier PHASOR's square `|X|^2 = A^2/4`) | `orbital_spectrum`, `modal_spectrum`, `correlation_spectrum` |
| **two-sided** `S_phi`, equal to `L(f)` in linear units (the IEEE one-sided `S_phi` is `2 L`) | `phase_psd`, `lorentzian` (integrates to 1 over `(-inf, inf)`) |
| dBc/Hz, `10 log10(S_v / (2 |X|^2))` = `L(f)` | `oscillator_spectrum` (`L_dBc`) |
| variance, unit^2 | `sampled_variance`, `covariance`, `oscillator_covariance` (`K_orb`), `orbital_correlation` |
| seconds (jitter), s^2 (growth) | `jitter_metrics`, `event_jitter`, `oscillator_edge_jitter`; `oscillator_covariance`'s `d` |
| diffusion constant `c`, seconds | `diffusion_constant`, `frequency_aware_diffusion`, `coloured_diffusion`, `coloured_diffusion_resolved` |

So near a carrier `pnoise ~ S_v ~ S_pm` (per sideband, where AM is
small), and the modal total is `~ S_v / 2`.

## Driven circuits

| method | arguments, in order | returns | frequency | key defaults |
|---|---|---|---|---|
| `pnoise` | `pss, freq, output, ratio_tol=None, maxsidebands=None, modulated=False, cyclostationary=False` | `(S, sidebands_used)`; also sets `self.alias_stop`, `self.sidebands_used` | `freq`: ABSOLUTE output frequency, scalar | `maxsidebands` -> `N//2` (clamped); `ratio_tol` -> 1e-9, stops after two quiet pairs |
| `am_pm_noise` | `pss, freq, output, carrier=1, maxsidebands=None, modulated=False` | `(S_am, S_pm, bands)`, PER SIDEBAND (half the pair's total; it was the total until 2026-09-28) | `freq`: OFFSET from `carrier f0` | `maxsidebands` -> every pair within Nyquist, `N//2 - carrier` (fixed 2026-09-28: it raised) |
| `band_spread` | `pss, output, band, points=9, harmonic=1, quantity='pnoise', **kw` | `(spread, info)` | `band` in units of f0: an OFFSET from `harmonic f0` for 'S_pm', 'S_am', 'oscillator_spectrum', but ABSOLUTE for 'pnoise' | `points=9` |
| `sampled_noise` | `pss, output, times, freqs, maxsidebands=None, tail=False` | array `(len(times), len(freqs))` | `freqs`: SERIES frequency, `0 < f <= f0/2` | `maxsidebands` -> `N//2 - 1` (raises above) |
| `sampled_variance` | `pss, output, times, fmin, fmax, points_per_decade=40, maxsidebands=None, tail=False` | array `(len(times),)` | `fmin`, `fmax`: the SERIES band, required; the band cuts white noise too | `ppd=40` |
| `jitter_metrics` | `pss, output, time, fmin, fmax, kmax=4, maxsidebands=None, nfreq=601, dc_rectangle=False, points_per_decade=40` | dict: `sigma_t, rho, k_cycle, cycle_to_cycle, slew, R, instant` | as `sampled_variance` | `kmax=4` (>= 2) |
| `covariance` | `pss, samples=False, fmin=None, fmax=None, points_per_decade=40` | `K0`, or `(K0, [K at each node])` with `samples=True`; `n x n` with `n = m`, or `2m` on gear/trap pair maps | `fmin`, `fmax`: the SOURCE band of the COLOURED part only; `fmax` -> the grid's Nyquist | `fmin` required when a source is coloured |
| `event_jitter` | `pss, fmin=None, fmax=None, points_per_decade=40` | dict: `sigma` (s), `cov_fraction`, `fractions`, `nodes` | as `covariance` | refuses a solve with no landed events |

## Oscillators

| method | arguments, in order | returns | frequency | key defaults |
|---|---|---|---|---|
| `oscillator_spectrum` | `pss, offsets, output, harmonic=1, frequency_aware=None, fmin=None, fmax=None, all_orders=None` | `(S_v, L_dBc)`; sets `self.lineshape_info` | OFFSETS from `harmonic f0`, any sign | `frequency_aware` -> True; `fmin` needed only with colour; `fmax` -> f0/2 |
| `phase_psd` | `pss, offsets, harmonic=1, frequency_aware=True` | array, rad^2/Hz (two-sided, = `L(f)`) | OFFSETS, `> 0` | refuses offsets at or below the Lorentzian corner |
| `lorentzian` (static) | `offsets, c, f0, harmonic=1` | array, 1/Hz relative to the harmonic's power | OFFSETS, any sign | harmonic 0 returns zeros |
| `modal_spectrum` | `pss, offsets, output, harmonic=1, H=None, sidebands=None` | dict: `phase, orbital, correlation, total` | OFFSETS, negative = the lower sideband | `H` -> 32 (capped by the grid), `sidebands` -> `2H` |
| `orbital_spectrum`, `correlation_spectrum` | `pss, offsets, output, harmonic=1, H=None` (+ `sidebands` for correlation) | array | as `modal_spectrum` | `H` -> 32 |
| `orbital_correlation` | `pss, H=None` | `(R, C)`: `R_yy(0)` `m x m`, `C[(l, h, j)]` | -- | refuses colour |
| `oscillator_covariance` | `pss, samples=False, fmin=None, fmax=None, points_per_decade=40` | `(K_orb, d, info)`; `K_orb` in pair space on gear/trap | as `covariance` (colour enters the transverse part only) | |
| `oscillator_edge_jitter` | `pss, output, time, kmax=8` | dict: `sigma_t, A, c, slew, k_cycle, instant, d, projection_share` | -- | `kmax=8`; `output` an index or a name |
| `orbital_mode_weights` | `pss, nmodes=None` | `(cw, modes, K_orb)` | -- | |
| `diffusion_constant` | `pss` | `c`, s | -- | refuses colour |
| `frequency_aware_diffusion` | `pss, offset` | `c(f)`, s | `offset` scalar, abs taken | refuses colour |
| `colour_projection` | `pss` | `(vbar, info)`: `rms, symmetry, samples, times` | -- | |
| `coloured_diffusion` | `pss, freqs` | `Gamma(f)`, the `l = 0` term ONLY | `freqs`: OFFSETS | do not add it to `c` (counts `l = 0` twice); use `coloured_diffusion_resolved` |
| `coloured_diffusion_resolved` | `pss, freqs, harmonics=None, frequency_aware=True` | `c(f)` over all harmonics | OFFSETS | `harmonics` -> every `l` with > 1e-14 of the PPV energy |

## Transfers (not noise)

| method | arguments | returns |
|---|---|---|
| `am_pm` | `pss, freq, output, carrier=1` | `(m_am, m_pm)` per unit source; `freq` an OFFSET; refuses a missing carrier |
| `am_pm_indices` | `a, b` | `(a + conj(b), a - conj(b))` |
| `carrier_phasor` | `pss, output, carrier=1` | the complex Fourier coefficient (`A/2` for `A cos`) |

## Where the surfaces disagree

**Can give a wrong number:**

1. **Two sidedness scales under one unit.** `pnoise`, `am_pm_noise`,
   `sampled_noise` and (since 2026-09-28) `oscillator_spectrum`'s `S_v` are
   one-sided; the orbital/modal/correlation spectra are 0.5x one-sided --
   HALF `S_v`. Normalising the latter by `A^2/2` instead of `|X|^2 = A^2/4`
   reads 3 dB low; only `oscillator_spectrum` returns dBc.
2. **`phase_psd` is the two-sided `S_phi` (= `L(f)`),** not the IEEE
   one-sided `S_phi = 2 L(f)`; its docstring says only "rad^2/Hz".
3. ~~**`am_pm_noise` returns pair totals** over both sidebands, not a
   per-sideband density; on an oscillator `S_pm = 4 S_v`.~~ CLOSED
   2026-09-28: per sideband, and on an oscillator `S_pm = S_v`.
4. **`freq` is ABSOLUTE in `pnoise` and an OFFSET in `am_pm_noise` / `am_pm`**
   -- the same name in the same position.
5. **`band_spread`'s `band` is an offset for three quantities and absolute
   for 'pnoise'**, and its `**kw` is not forwarded to
   `oscillator_spectrum` although the docstring says so; an unknown
   `quantity` silently falls through to 'pnoise'.
6. **`fmin`/`fmax` name three different bands:** the coloured SOURCE band
   (`covariance`, `event_jitter`, `oscillator_covariance`; ignored when all
   sources are white), the SERIES band that cuts white noise too
   (`sampled_variance`, `jitter_metrics`), and the phase-offset band for
   colour (`oscillator_spectrum`).
7. **`covariance` is `m x m` on radau/euler/GLM and `2m x 2m` on gear/trap**
   (the pair state), so traces, eigenvalues and shapes depend on the method.
8. **The jitter dictionaries disagree:** `sigma_t` vs `sigma`
   (`event_jitter`), `k_cycle` exact in `jitter_metrics` but a large-`k`
   upper bound in `oscillator_edge_jitter`, `slew` from a central
   difference vs a quadratic fit, `kmax` 4 vs 8.
9. **Flag combinations dropped silently:** `pnoise` with both `modulated`
   and `cyclostationary` (cyclostationary wins); `oscillator_spectrum` with
   `frequency_aware=False` drops `all_orders`; `all_orders=None` means
   different things for white and coloured sources; `frequency_aware`
   defaults to None in `oscillator_spectrum` and True elsewhere.
10. **A missing carrier:** `am_pm` refuses it, `am_pm_noise` proceeds
    unrotated, `lorentzian` returns zeros for harmonic 0, the oscillator
    spectra refuse.
11. **Argument order:** `(pss, freq/offsets, output)` for `pnoise`,
    `am_pm_noise` and the spectra, but `(pss, output, ...)` for the sampled
    family, `band_spread` and `oscillator_edge_jitter`. Because `output`
    accepts a vector, a swapped array can be taken as weights silently.

**Naming only:**

12. The truncation knob has many names: `maxsidebands` (`N//2`, clamped, in
    `pnoise`; `N//2 - 1`, raising, in the sampled family), `H`,
    `sidebands` (a count in `modal_spectrum`, a list in
    `adjoint_sideband_row`), `harmonics`, `kmax`, `nmodes`.
13. The harmonic is `harmonic` in the spectra and `carrier` in `am_pm`,
    `am_pm_noise`, `carrier_phasor`.
14. Frequency names: `freq` (scalar), `freqs` (a series frequency or an
    offset), `offsets`, `offset`, `band` (units of f0), `time` vs `times`;
    offset signs differ (`phase_psd` refuses <= 0, `oscillator_spectrum`
    accepts any, modal/orbital read negative as the lower sideband).
15. Returns: tuples, dicts, arrays and floats; `covariance` changes arity
    with `samples` while `oscillator_covariance` puts samples in `info`;
    diagnostics go to instance attributes (`alias_stop`, `sidebands_used`,
    `sampled_instants`, `lineshape_info`, ...).
16. The same refusal raises different exceptions: a missing `fmin` is a
    `NotImplementedError` in `covariance` and `oscillator_spectrum` but a
    `TypeError` in the sampled family (a required argument).
17. Several oscillator surfaces cannot take a coloured source at all
    (`oscillator_edge_jitter`, `orbital_mode_weights` have no `fmin` to
    pass; `orbital_*` refuse colour; `modal_spectrum` accepts it).
18. Stale text: `colour_projection`'s docstring promises "the ratio
    |mean|/rms" under a key that is `symmetry`; three comments still
    describe a Gear-2 fallback for TR-BDF2 that the accuracy host no longer
    uses (it uses the monodromy twin).
