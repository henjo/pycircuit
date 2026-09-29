"""Oscillator phase noise: the diffusion constants (white, coloured, frequency-
aware), phase_psd and the lineshape.
"""
import numpy as np
import warnings
from ._noise_components import psd_sqrt
from ._numerics import insert_ref, output_index


class _PhaseNoise(object):
    """Oscillator phase noise: the diffusion constants (white, coloured,
    frequency-aware), phase_psd and the lineshape.  A theme of `PAC` (see
    `pac.py`)."""

    def diffusion_constant(self, pss):
        """`c` — the phase diffusion constant, in seconds.

        `c = (1/T) ∫ v₁ᵀ(t) B(t) Bᵀ(t) v₁(t) dt` with `B Bᵀ = CY`, so this
        is the time-average of a QUADRATIC form in the PPV.  It is the one
        scalar the whole free-running phase-noise spectrum is built from,
        and it reads, for a designer, as JITTER PER SECOND.

        ⚠ QUADRATIC FOR WHITE SOURCES, LINEAR FOR COLOURED ONES, and the
        two are different functionals of the same vector: a coloured
        source contributes `V_0m = (1/T) ∫ v₁ᵀ B_cm dt`, with no square.
        Using this one for a coloured source returns a plausible non-zero
        number from the same PPV.  Only WHITE sources are supported here.

        A MODULATED white source (its `CY` follows the orbit: a MOS
        channel's `4kT gamma g_d0(x)`, a shot noise `2qI(x)`) is Demir's
        own case, `B = B(x(t))`: `CY` is read at each PPV sample's state.

        ⚠ `CY/2`, as in `covariance`: settled against `kT/C`, which is
        external to both (an injection of `Var(i) = CY/h` per step
        reproduces 1.92x `kT/C`).  A Monte Carlo built on the convention
        under test cannot test it.

        ⚠ `ppv()` normalises on the FIRST BLOCK (`v[:m] . xdot = 1`), which
        is right for a perturbation entering the first block -- an injected
        current, and what every shipped path does -- but wrong for
        contracting against a full PAIR deviation, where the factor is
        `1/(v . u_pair)`.  The sign is a convention: a later zero crossing
        means DELAYED, while projecting onto the tangent makes positive
        mean ADVANCED.

        History: `doc/shooting_history.md`, `PAC.diffusion_constant`.
        """
        self._check_circuit(pss)
        self._refuse_coloured(
            pss, 'diffusion_constant',
            'a coloured source has no diffusion constant; phase_psd and '
            'coloured_diffusion_resolved give its phase noise.')
        self._refuse_driven(pss, 'diffusion_constant')
        return self._white_diffusion_at(pss, 2.0 * np.pi / float(pss.period))

    def _white_diffusion_at(self, pss, w, cy=None):
        """`(1/T) integral v_1^T (CY(w)/2) v_1 dt` with `CY` FROZEN at `w`.

        The white functional at one frequency, with no refusal: it is `c`
        when the source is white, and for a coloured source it is the
        value `diffusion_constant` refuses.  `phase_psd`
        reads it at the carrier for the Lorentzian CORNER, which is a
        white-noise construct whatever the source's colour; the spectrum
        itself comes from `coloured_diffusion_resolved`.

        History: `doc/shooting_history.md`, `PAC._white_diffusion_at`.
        """
        v, info = pss.ppv()
        ## `lambda_2` is computed here anyway; `oscillator_spectrum` needs it to
        ## report its own validity limit and a second `ppv()` would be a full
        ## extra solve.  Recorded, not returned, so this method's signature is
        ## unchanged -- and read ONLY immediately after a call, which is how
        ## `oscillator_spectrum` uses it.
        self._last_second_multiplier = (
            info.get('second_multiplier'),
            info.get('second_multiplier_certified'))
        m = pss.cir.n - 1
        ## ⚠ `samples_eq`, NOT `samples`.  `CY` is an EQUATION-ROW
        ## covariance and `samples` is `C^T v_1`; contracting that would make
        ## `c` wrong by `C^2` on the differential rows and exactly zero on
        ## the algebraic ones.  See `_equation_row_ppv`.
        S = np.asarray(info['samples_eq'])[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        ## the samples' own orbit, not `pss.period` -- see `ppv()`'s 'period'
        T = float(info['period'])
        h = self._period_weights(tms, S.shape[0], T, pss)
        ## `cy` given: the functional of THAT source matrix (the white part
        ## of a coloured circuit, `_modal_colour`), one matrix or one per
        ## PPV sample `(N, m, m)`
        ## ⚠ A MODULATED source: `CY` at each sample's OWN state -- Demir's
        ## `B(x(t))` -- where `_cy_reduced` refuses.
        ## The states are the PPV's orbit (`_ppv_states`: a twin's where one
        ## serves the PPV).
        ## History: `doc/shooting_history.md`, `PAC._white_diffusion_at`.
        if cy is None:
            try:
                cy = self._cy_reduced(pss, float(w))
            except NotImplementedError:
                cy = self._noise_components(
                    pss, self._ppv_states(pss)).cy_at_states(float(w))
        cy = np.asarray(cy)
        if cy.ndim == 3 and cy.shape[0] != S.shape[0]:
            raise ValueError(
                'PAC: %d source samples against %d PPV samples -- they must '
                'be taken at the same orbit states.' % (cy.shape[0], S.shape[0]))
        ## ⚠ A NOISE SOURCE ON AN INDEX-2 CONSTRAINT GIVES c = 0, SILENTLY: a
        ## voltage noise in series with a DC source inside a capacitor loop
        ## perturbs an algebraic constraint -- a DIFFERENTIATED input, whose
        ## response is a charge jump the PPV projection cannot represent --
        ## and its share of `c` is exactly 0 for every method.  Named here
        ## once: the PPV's algebraic fallback fired (index >= 2) and `CY` has
        ## power on an algebraic row.  (The gate is the index-2 condition
        ## itself -- `G[A,Z]` singular at the orbit point, the test the PPV's
        ## algebraic fallback makes -- computed here so every kind is
        ## covered.)
        try:
            x0r = np.asarray(pss._period_state[1], dtype=float).ravel()
            irn = pss.irefnode
            x0f = insert_ref(x0r, irn)
            _arows, _acols = pss._algebraic_adjoint_pattern(x0f)
            _idx2 = False
            if _arows and _acols and len(_arows) == len(_acols):
                _Gz = np.asarray(pss._G_at(x0r), dtype=float)[np.ix_(
                    np.asarray(_arows, dtype=int), np.asarray(_acols, dtype=int))]
                _sv = np.linalg.svd(_Gz, compute_uv=False)
                _idx2 = (float(_sv[-1]) <= 1e-12 * max(float(_sv[0]), 1e-300))
            elif _arows:
                _idx2 = True
        except Exception:                                      # noqa: BLE001
            _idx2 = False
        if _idx2:
            try:
                _dcy = np.abs(np.real(np.diagonal(cy, axis1=-2, axis2=-1)))
                if _dcy.ndim == 2:
                    _dcy = _dcy.max(axis=0)
                _on_alg = [r for r in _arows if _dcy[r] > 0.0]
                if _on_alg and float(np.max(_dcy)) > 0.0:
                    warnings.warn(
                        'PAC.diffusion_constant: this circuit is index >= 2 and '
                        'a noise source sits on an algebraic row (%s). A '
                        'perturbation of an index-2 constraint is a '
                        'differentiated input -- a charge jump -- which the '
                        'PPV projection cannot represent, and its share of c '
                        'is 0 here whatever the source (measured: exactly 0 '
                        'for every method). Only the differential rows\' '
                        'sources are counted.' % (_on_alg,),
                        RuntimeWarning, stacklevel=3)
            except Exception:                                  # noqa: BLE001
                pass
        ## ⚠ `cy/2`, THE SAME ONE-SIDED-TO-TWO-SIDED CONVERSION `covariance`
        ## USES.  `CY` is a one-sided density (a resistor's `4kT/R`); an
        ## injection of `Var(i) = CY/h` per step reproduces `1.92x kT/C`, so
        ## that convention carries TWICE the physical noise power.
        if cy.ndim == 3:
            quad = np.einsum('ij,ijk,ik->i', S, 0.5 * np.real(cy), S)
        else:
            quad = np.einsum('ij,jk,ik->i', S, 0.5 * np.real(cy), S)
        return float((quad * h).sum() / T)

    def colour_projection(self, pss):
        """`<v_1>` — the PPV's TIME AVERAGE, which is a different functional.

        Returns `(vbar, info)`.  `vbar` is `(1/T) integral v_1(t) dt` over
        the orbit; `info` carries the per-entry `rms`, the ratio
        `|mean|/rms` under `'symmetry'` -- the number that says whether a
        coloured source at that node can upconvert at all -- and the PPV
        `samples` at `times`.

        ⚠ COLOURED SOURCES CONTRACT THE SQUARE OF THE MEAN; WHITE ONES
        CONTRACT THE MEAN OF THE SQUARE.  `diffusion_constant` computes
        `(1/T) integral v^T (CY/2) v dt`.  A coloured source's low-frequency
        power cannot be modulated away, so what survives is
        `V_0m = (1/T) integral v_1^T B_cm dt` -- LINEAR, no square -- and the
        contraction is `vbar^T (CY/2) vbar`.  Same vector, same matrix,
        the mean and the square exchanged.

        ⚠ AND USING THE QUADRATIC ONE FOR A COLOURED SOURCE RETURNS A
        PLAUSIBLE NUMBER, NOT AN ERROR.  It is never zero where the white
        answer is not, so nothing downstream would look wrong -- and the
        two are not close approximations of each other (22 orders apart on
        van der Pol), so they cannot be substituted.

        ⚠ TWO INDEPENDENT MECHANISMS FORCE `vbar` TO ZERO, AND ONLY ONE OF
        THEM IS THE ONE DESIGNERS KNOW.  NEITHER ASYMMETRY ALONE NOR LOSS
        ALONE UPCONVERTS (on an LC oscillator with an even nonlinearity
        term and a series tank resistance, `Gamma/c` ~ 1e-22 with either
        alone and 1e-4 .. 1e-3 with both, while `c` barely moves):

        * A SYMMETRIC waveform gives `vbar = 0` -- Hajimiri & Lee, and the
          reason symmetry is the first thing a VCO designer reaches for.
        * A LOSSLESS LC TANK gives `vbar[0] = 0` STRUCTURALLY, whatever the
          waveform does.  `v` behaves as `C^T v_1` and `dv/dt = G^T v_1`,
          whose inductor row is exactly `v[0]`; periodicity of `v[1]` then
          forces `integral v[0] dt = 0`.  ⚠ THIS IS A PROPERTY OF THE
          TOPOLOGY, NOT OF THE ORBIT, and it is why van der Pol reports
          zero at every asymmetry -- it makes van der Pol useless as a
          POSITIVE fixture and perfect as a negative one.

        ⚠ `Gamma <= c` ALWAYS, at the same `CY`, by Cauchy-Schwarz on the
        weighted mean -- with equality only if `v` is constant over the
        orbit.  Both use the same quadrature here so the bound holds
        exactly at the discrete level, which makes it an assertion rather
        than an expectation.

        History: `doc/shooting_history.md`, `PAC.colour_projection`.
        """
        self._check_circuit(pss)
        if not getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.colour_projection: the PPV time-average is the kernel '
                "of a FREE-RUNNING oscillator's coloured-noise upconversion. "
                "A driven circuit's phase is its source's.")
        _v, info = pss.ppv()
        m = pss.cir.n - 1
        ## ⚠ the EQUATION-ROW adjoint, for the same reason
        ## `diffusion_constant` uses it: a coloured source is an
        ## equation-row input too.
        S = np.asarray(info['samples_eq'])[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        T = float(info['period'])
        h = self._period_weights(tms, S.shape[0], T, pss)
        ## ⚠ THE SAME QUADRATURE `diffusion_constant` USES, deliberately:
        ## it is what makes `Gamma <= c` exact rather than approximate.
        vbar = (S * h[:, None]).sum(0) / T
        rms = np.sqrt((S ** 2 * h[:, None]).sum(0) / T)
        with np.errstate(divide='ignore', invalid='ignore'):
            sym = np.where(rms > 0, np.abs(vbar) / rms, 0.0)
        return vbar, {'rms': rms, 'symmetry': sym,
                      'samples': S, 'times': tms}

    def coloured_diffusion(self, pss, offsets):
        """`Gamma(f) = vbar^T (CY(2 pi f)/2) vbar` — the coloured analogue of `c`.

        Returns an array over `offsets` (from the harmonic; `freqs`
        until 2026-09-29), of either sign and even in them.  `CY` is evaluated at each offset,
        so a source whose density varies with frequency -- which is what
        "coloured" means -- is folded in exactly as a white one is.

        ⚠ NO FILTER, NO EXTRA STATE, NO SDE.  Demir 1996 synthesises 1/f
        from white sources through a Lorentzian network at "one state
        variable per decade", because Ito theory admits only white driving
        noise.  That is an artefact of the SDE formulation.  This path
        never forms an SDE, so a coloured source is just a different
        `S(f)` -- a SLOPE, NOT A STATE.  A commercial RF simulator confirms by omission:
        no filter and no augmentation in its treatment of flicker.

        ⚠ THE `CY/2` IS THE SAME ONE-SIDED-TO-TWO-SIDED CONVERSION THE
        REST OF THIS CLASS USES, and it is shared rather than repeated so
        the functions cannot drift apart over that factor.

        A source that FOLLOWS THE ORBIT: the
        l = 0 term of the modulated fold, `(s(f)/2) |<v_1^T G>|^2` per
        component, `G` its columns (signed where the element states them,
        `colour_groups`) -- `vbar^T (CY/2) vbar` when `G` does not move, and
        never above `coloured_diffusion_resolved` (Jensen).  ⚠ A modulated
        WHITE part enters through the root of its PSD, a convention: the
        sign of a white source is unobservable (`g xi` and `|g| xi` are one
        process), so its share of this DC term is not a property of the
        noise -- only `c` is.  (Its stationary realisation through a
        signed multiplier reads a different share.)

        History: `doc/shooting_history.md`, `PAC.coloured_diffusion`.
        """
        vbar, _ = self.colour_projection(pss)
        if self._modulated_present(pss):
            ## ⚠ A source that follows the orbit: the l = 0 term
            ## of the modulated fold, `(s(f)/2) |<v_1^T G>|^2` per component
            ## (`_colour_fold`) -- `vbar^T (CY/2) vbar` when `G` does not move
            ## History: `doc/shooting_history.md`, `PAC.coloured_diffusion`.
            fr = np.atleast_1d(np.asarray(offsets, dtype=float))
            pos = np.abs(fr[fr != 0.0])
            return self._colour_fold(
                pss, float(pos.min()) if pos.size else None, None,
                'coloured_diffusion').dc(fr)
        out = []
        for f in np.atleast_1d(np.asarray(offsets, dtype=float)):
            cy = np.real(self._cy_reduced(pss, 2.0 * np.pi * float(f)))
            out.append(float(vbar @ (0.5 * cy) @ vbar))
        return np.asarray(out)

    def coloured_diffusion_resolved(self, pss, offsets, harmonics=None,
                                    frequency_aware=True):
        """`c(f) = sum_l V_l^H (CY(2 pi |f - l f_0|)/2) V_l` — the fold PER HARMONIC.

        `V_l` are the Fourier coefficients of the equation-row PPV `v_1(t)`
        (the rows `diffusion_constant` contracts), so a source's density is
        read at the SOURCE-SIDE frequency `f - l f_0` for each harmonic it
        folds through -- which is what `pnoise` has done from the start and
        what a coloured source requires.  Returns an array over `offsets` (`freqs`
        until 2026-09-29), of either sign and even in them (`V_{-l}` is
        `conj(V_l)`; measured equal at +-o on an asymmetric orbit).

        ⚠ `c + Gamma(f)` IS NOT THIS: `c` reads `CY` at ONE frequency
        (`2 pi / T`) as if it held at every harmonic, and `Gamma` is exactly
        the `l = 0` term of this sum, so `c + Gamma` counts `l = 0` twice.
        Neither shows on van der Pol, whose PPV at the tank node averages to
        zero (the inductor shorts the node at DC, so no core can bias it).

        EXACT FOR WHITE, BY PARSEVAL: with `CY` constant the sum is
        `(1/T) integral v_1^T (CY/2) v_1 dt = c`, and the discrete version
        with the grid's step weights reproduces `diffusion_constant` to
        round-off -- that equality pins the transform's normalisation, and
        it is asserted.  For a DC-centred colour (Lorentzian, flicker) and
        `f << f_0` the `l != 0` terms read `CY(l f_0)` to `O(f/f_0)`, so
        the sum differs from `c + Gamma` only where `V_0` is not small.

        `harmonics` caps `|l|`; by default every harmonic carrying more
        than 1e-14 of the PPV's energy is kept, which is all of them that
        can move the sum at double precision.

        ⚠ `frequency_aware` (default True): `V_l` from the
        frequency-aware PPV at each offset (`PSS.frequency_aware_ppv`, one
        bordered solve per offset), with the same bands `f - l f0`.  The DC
        PPV assumes a source moves the phase instantly, and behind a slow
        path it over-states by the path's filter (pnoise PM: DC 0.729 /
        0.026 at 1e-3 / 1e-2 f0 behind a tau = 100 T node, this 1.0001 /
        1.0003, white, Lorentzian or 1/f alike).  For a white source this is
        `frequency_aware_diffusion(f)` exactly, not `c` -- `False` keeps the
        DC fold, which reproduces `diffusion_constant` (the Parseval
        identity below).

        A source whose level follows the orbit takes
        `_coloured_diffusion_modulated`: the harmonics of the PRODUCT
        `v_1^T G(t)` per component, and the white parts as Demir's `c`.

        History: `doc/shooting_history.md`, `PAC.coloured_diffusion_resolved`.
        """
        self._check_circuit(pss)
        self._refuse_driven(pss, 'coloured_diffusion_resolved')
        if self._modulated_present(pss):
            return self._coloured_diffusion_modulated(pss, offsets, harmonics,
                                                      frequency_aware)
        m = pss.cir.n - 1
        v0, info = pss.ppv()
        S = np.asarray(info['samples_eq'], dtype=float)[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        n = S.shape[0]
        ## the samples' orbit: both the 1/T AND the harmonic frequency `w0`
        ## must come from it, or Parseval leaks (1.5e-09 under a trap twin)
        T = float(info['period'])
        ## ⚠ THE SAME QUADRATURE `diffusion_constant` USES: one sample per
        ## step, weighted by that step, so that Parseval closes exactly.
        ## ⚠ Sample `j` is at `t_j` (the PPV pairs it with the orbit's column
        ## `j`).  One step late is identical on a uniform grid (a common
        ## phase), but on a smoothly varying one it caps the fold's order.
        ## History: `doc/shooting_history.md`, `PAC.coloured_diffusion_resolved`.
        t = tms[:n]
        h = self._period_weights(t, n, T, pss)
        w0 = 2.0 * np.pi / T
        L = n // 2 if harmonics is None else int(harmonics)
        ls = np.arange(-L, L + 1) if harmonics is not None else np.arange(-L, L)
        E = np.exp(-1j * np.outer(ls, w0 * t)) * h[None, :]          # (nl, n)
        if frequency_aware:
            ## ⚠ THE FREQUENCY-AWARE PPV (the default): at each
            ## offset the coefficients of `frequency_aware_ppv(f)`'s periodic
            ## samples, with the SAME bands `f - l f0` -- measured, not read:
            ## a source peaked at f0 + f against its white-through-resonator
            ## realisation reads 0.9999 this way and 1.083 with `f + l f0`.
            ## Behind a slow node the DC PPV over-states by the path's filter
            ## (pnoise PM / phase_psd: DC 0.729 / 0.026 at 1e-3 / 1e-2 f0,
            ## this 1.0001 / 1.0003 for a Lorentzian or 1/f source).
            out = []
            for f in np.atleast_1d(np.asarray(offsets, dtype=float)):
                Sf = S if f == 0.0 else self._fa_samples(pss, f)
                Vf = (E @ Sf) / T
                en = np.sum(np.abs(Vf) ** 2, axis=1)
                tot = 0.0
                for l, vl in zip(ls[en > 1e-14 * en.sum()],
                                 Vf[en > 1e-14 * en.sum()]):
                    cy = np.real(self._cy_reduced(
                        pss, 2.0 * np.pi * abs(float(f) - l / T)))
                    tot += float(np.real(np.conj(vl) @ (0.5 * cy) @ vl))
                out.append(tot)
            return np.asarray(out)
        V = (E @ S) / T                                               # (nl, m)
        energy = np.sum(np.abs(V) ** 2, axis=1)
        keep = energy > 1e-14 * energy.sum()
        ls, V = ls[keep], V[keep]
        out = []
        for f in np.atleast_1d(np.asarray(offsets, dtype=float)):
            tot = 0.0
            for l, vl in zip(ls, V):
                cy = np.real(self._cy_reduced(pss, 2.0 * np.pi * abs(float(f) - l / T)))
                tot += float(np.real(np.conj(vl) @ (0.5 * cy) @ vl))
            out.append(tot)
        return np.asarray(out)

    def _fa_samples(self, pss, f):
        """`frequency_aware_ppv(f)`'s equation-row samples, `(n, m)` complex:
        the periodic envelope on the DC PPV's times, its Fourier coefficient
        `l` the phase transfer of the source band at `f - l f0` (the DC fold's
        convention; see `coloured_diffusion_resolved`).  Cached per offset."""
        cache = self.__dict__.setdefault('_fa_cache', {})
        fp = pss.factored_period()
        key = (id(fp), float(f))
        hit = cache.get(key)
        ## ⚠ THE ENTRY HOLDS ITS FACTORED PERIOD AND IS MATCHED BY IDENTITY,
        ## as `_transverse_cache` is: keyed on the id alone, a re-solve frees
        ## the old period, a new one can be born at the same address, and the
        ## OLD grid's samples come back.
        ## History: `doc/shooting_history.md`, `PAC._fa_samples`.
        if hit is None or hit[0] is not fp:
            m = pss.cir.n - 1
            _v, fi = pss.frequency_aware_ppv(float(f))
            hit = cache[key] = (fp, np.asarray(fi['samples_eq'])[:, :m])
        return hit[1]

    def _coloured_diffusion_modulated(self, pss, freqs, harmonics=None,
                                      frequency_aware=False):
        """`coloured_diffusion_resolved` for sources that follow the orbit.
        The phase moves as `v_1(t)^T G(t) xi(t)` for each
        component, `xi` a unit process of power `s(nu)` and `G` its columns
        at the PPV's own states, so

            c(f) = c_white + sum_groups sum_l (s(|f - l f0|)/2) |V_l|^2

        with `V_l` the Fourier coefficients of the PRODUCT `v_1^T G` -- a
        modulation's harmonics shift the source-side frequency the way
        the PPV's do.  The white parts give `c_white`, Demir's `c` with
        `CY` at each sample's state, flat in `f` (a white process modulated
        is still white).  The same period weights as the stationary form,
        with the samples at their own times.

        History: `doc/shooting_history.md`, `PAC._coloured_diffusion_modulated`."""
        fr = np.atleast_1d(np.asarray(freqs, dtype=float))
        pos = np.abs(fr[fr != 0.0])
        fold = self._colour_fold(pss, float(pos.min()) if pos.size else None,
                                 harmonics, 'coloured_diffusion_resolved')
        if frequency_aware:
            return fold.resolved_fa(fr)
        return fold.coloured(fr, start=fold.c_white)

    def _colour_fold(self, pss, flo, harmonics, what):
        """The pieces of `_coloured_diffusion_modulated`, built once: returns
        an object with `c_white` (Demir's `c` of the WHITE parts at the PPV's
        states) and `coloured(freqs, start=0.0)`, the coloured part of
        `c(f)` (vectorised over `freqs`, added to `start`).  A STATIONARY
        source is the case whose columns do not move.  `flo`: the lowest
        frequency the caller will ask (for the per-band classification)."""
        m = pss.cir.n - 1
        _v0, info = pss.ppv()
        S = np.asarray(info['samples_eq'], dtype=float)[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        n = S.shape[0]
        T = float(info['period'])
        t = tms[:n]
        h = self._period_weights(t, n, T, pss)
        f0 = 1.0 / T
        w0 = 2.0 * np.pi * f0
        L = n // 2 if harmonics is None else int(harmonics)
        ls = np.arange(-L, L + 1) if harmonics is not None else np.arange(-L, L)
        E = np.exp(-1j * np.outer(ls, w0 * t)) * h[None, :] / T      # (nl, n)
        states = self._ppv_states(pss)
        nc = self._noise_components(pss, states)
        model = nc.model(1e-3 * f0, f0)
        if model is None:
            raise NotImplementedError(
                'PAC.%s: this circuit\'s CY is not the sum of its elements\' '
                'and, as a whole, not thermal-plus-power-law, so its white and '
                'coloured parts cannot be separated.' % what)
        c_white = float(self._white_diffusion_at(
            pss, w0, cy=np.real(np.asarray(model.white))))
        wlo = 2.0 * np.pi * float(flo) if flo else 1e-3 * w0
        fixed, band = [], []
        groups_all = nc.colour_groups(model, wlo, f0, L, what)
        for kind, G, s in groups_all:
            if kind == 'fixed':
                V = E @ np.einsum('jm,jmr->jr', S, G)
                pw = np.sum(np.abs(V) ** 2, axis=1)
                keep = pw > 1e-14 * pw.sum()
                fixed.append((pw[keep], ls[keep], s, G))
                continue
            ## per band: the harmonics that carry weight at the carrier's
            ## band, then each read at its own frequency
            Vr = E @ np.einsum('jm,jmr->jr', S, G(w0))
            pr = np.sum(np.abs(Vr) ** 2, axis=1)
            band.append((G, np.where(pr > 1e-14 * pr.sum())[0]))

        class _Fold:
            pass
        fold = _Fold()
        fold.c_white = c_white
        fold.has_colour = bool(fixed or band)
        ## the l = 0 row: `<v_1^T G>`, the PPV's time average through the
        ## columns -- `coloured_diffusion`'s `vbar^T (CY/2) vbar` when `G`
        ## does not move.  A WHITE part enters through the root of its PSD
        ## (its share of the DC term is a convention: the sign of a white
        ## source is unobservable, and only `c` is physical for it).
        E0 = E[int(np.flatnonzero(ls == 0)[0])]
        white_dc = 0.0
        for _key, A in model.white_parts:
            v0 = E0 @ np.einsum('jm,jmr->jr', S,
                                psd_sqrt(np.real(np.asarray(A))))
            white_dc += 0.5 * float(np.sum(np.abs(v0) ** 2))
        fixed_dc = [(float(np.sum(np.abs(E0 @ np.einsum('jm,jmr->jr', S, G))
                                  ** 2)), s)
                    for kind, G, s in groups_all if kind == 'fixed']
        band_G = [G for kind, G, _s in groups_all if kind == 'band']

        def dc(freqs):
            fr = np.atleast_1d(np.asarray(freqs, dtype=float))
            out = np.full(fr.shape, white_dc)
            for p0, s in fixed_dc:
                out += 0.5 * p0 * np.asarray(s(2.0 * np.pi * np.abs(fr)),
                                             dtype=float)
            for G in band_G:
                for i, f in enumerate(fr):
                    v0 = E0 @ np.einsum('jm,jmr->jr', S,
                                        G(2.0 * np.pi * abs(float(f))))
                    out[i] += 0.5 * float(np.sum(np.abs(v0) ** 2))
            return out
        fold.dc = dc

        def coloured(freqs, start=0.0):
            fr = np.atleast_1d(np.asarray(freqs, dtype=float))
            out = np.full(fr.shape, float(start))
            for pw, lk, s, _G in fixed:
                for i, f in enumerate(fr):
                    nu = 2.0 * np.pi * np.abs(f - lk * f0)
                    out[i] += 0.5 * float(np.sum(np.asarray(s(nu)) * pw))
            for G, idx in band:
                for i, f in enumerate(fr):
                    for k in idx:
                        Vk = E[k] @ np.einsum(
                            'jm,jmr->jr', S,
                            G(2.0 * np.pi * abs(f - ls[k] * f0)))
                        out[i] += 0.5 * float(np.sum(np.abs(Vk) ** 2))
            return out
        fold.coloured = coloured
        white_A = 0.5 * np.real(np.asarray(model.white))

        def resolved_fa(freqs, parts=False):
            ## the whole `c(f)` with the FREQUENCY-AWARE PPV's samples per
            ## offset (`_fa_samples`): the white parts Demir's functional on
            ## them, each coloured group through its columns as above.
            ## `parts`: `(white, coloured)` instead of the sum (the lineshape
            ## cuts the colour at `fmin` and not the white)
            fr = np.atleast_1d(np.asarray(freqs, dtype=float))
            out = np.zeros(fr.shape)
            wpart = np.zeros(fr.shape)
            cpart = np.zeros(fr.shape)
            for i, f in enumerate(fr):
                Sf = S if f == 0.0 else self._fa_samples(pss, f)
                out[i] = wpart[i] = float(np.sum(h * np.real(np.einsum(
                    'jm,jmk,jk->j', np.conj(Sf), white_A, Sf)))) / T
                for _pw, _lk, s, G in fixed:
                    V = E @ np.einsum('jm,jmr->jr', Sf, G)
                    p = np.sum(np.abs(V) ** 2, axis=1)
                    k = p > 1e-14 * p.sum()
                    nu = 2.0 * np.pi * np.abs(f - ls[k] * f0)
                    term = 0.5 * float(np.sum(np.asarray(s(nu)) * p[k]))
                    out[i] += term
                    cpart[i] += term
                for G, idx in band:
                    for k in idx:
                        Vk = E[k] @ np.einsum(
                            'jm,jmr->jr', Sf,
                            G(2.0 * np.pi * abs(f - ls[k] * f0)))
                        term = 0.5 * float(np.sum(np.abs(Vk) ** 2))
                        out[i] += term
                        cpart[i] += term
            return (wpart, cpart) if parts else out
        fold.resolved_fa = resolved_fa
        return fold

    def phase_psd(self, pss, offsets, harmonic=1, frequency_aware=True):
        """`L(f)`, the SINGLE-SIDEBAND phase noise at `offsets` from harmonic
        `i`, per Hz relative to the carrier (`10 log10` of it is dBc/Hz) --
        white AND coloured.  `offsets > 0`: a zero or negative offset is
        refused (`L` diverges at zero offset, and that is physical).

            L_i(f) = i^2 f_0^2 c(f) / f^2,   c(f) = sum_l V_l^H (CY(f - l f_0)/2) V_l

        ⚠ THE CONVENTION IS A COMMERCIAL SIMULATOR'S PHASE NOISE, and not
        the IEEE one-sided `S_phi(f) = 2 L(f)` (rad^2/Hz).  Its pnoise
        reports an oscillator's phase noise as `L(f)` alone, the Lorentzian
        `c f0^2 / ((pi c f0^2)^2 + f^2)` whose skirt this is; it has no
        rad^2/Hz output.  (Until 2026-09-29 this docstring called the same
        number "the two-sided S_phi in rad^2/Hz".)

        `c(f)` is `coloured_diffusion_resolved`: the phase diffusion with
        each harmonic's colour read at its own source-side frequency.  For
        a white source it is `c` exactly (`c + Gamma(f)` would count the
        `l = 0` term twice).

        ⚠ THE CONVENTION IS PINNED BY `oscillator_spectrum`, NOT ARGUED.
        `lorentzian`'s far skirt is `i^2 f_0^2 c / f^2` exactly, and that
        object was gated by power conservation to 1.000000.  So this
        expression is the same quantity its tail already reports, with the
        coloured term added -- no second convention is introduced.

        ⚠ THE TWO TERMS ADD BECAUSE THE SOURCES ARE INDEPENDENT, and with
        `CY ~ 1/f` the coloured term gives `S_phi ~ 1/f^3` -- Kundert's
        "S_u(f) is generally pink ... then S_phi(f) would be proportional
        to 1/f^3 at low frequencies".

        ⚠ AND THIS IS THE LINEARISED PHASE MODEL, WHICH IS EXACT ENOUGH
        ONLY BECAUSE THE SOURCES ARE STATIONARY.  Vanassche, Gielen &
        Sansen (ICCAD 2002) locate the split between the exact phase
        equation `theta' = eps Gamma(t + theta) n(t)` and the approximate
        `theta' = eps Gamma(t) n(t)`: for a STATIONARY source the two
        "will, up to 0-th order in eps, predict the same output phase
        noise", and they diverge otherwise.  Their operational form is
        better than "non-stationary" -- "at first, near t = 0, the
        predicted phases are the same. However, when THETA BECOMES TOO
        LARGE [they diverge]".  A stationary source makes `theta` DIFFUSE;
        a driven one makes it grow SECULARLY, which is what carries it out
        of range.  ⚠ So this is sound for free-running noise and must NOT
        be reused for injection locking, a PLL in lock, or coupled
        oscillators -- there the shift has to stay inside the argument.
        ⚠ A source MODULATED BY THE OSCILLATOR'S OWN STATE is
        the stationary case in this sense: its level `B(x(t + theta))`
        moves with the phase, so `v^T B` is one periodic function of
        `t + theta` driven by a stationary process -- Demir's own form.  A
        modulation by an EXTERNAL clock would not be; only a driven circuit
        has one, and those are refused.

        ⚠ REFUSED BELOW THE LORENTZIAN CORNER, and this is a validity
        boundary rather than a conditioning one.  There the excess phase is
        a Wiener process whose spectrum is singular at the origin; the
        finite value the real lineshape attains comes from the NONLINEAR
        phase-to-voltage map, which `oscillator_spectrum` carries and this
        does not.  Reporting `S_phi` near the carrier is the mistake this
        object invites, so it raises instead.

        ⚠ `frequency_aware` (default True): `c(f)` from the
        frequency-aware PPV (`coloured_diffusion_resolved`), as
        `oscillator_spectrum` takes it for white sources -- the DC PPV
        over-states a source behind a slow node by the path's filter.  The
        Lorentzian corner and the power bound's probe stay on the DC PPV:
        the corner is a white-noise construct at the carrier, and the probe
        over-states, which keeps the bound it serves conservative.

        History: `doc/shooting_history.md`, `PAC.phase_psd`.
        """
        self._check_circuit(pss)
        f0 = 1.0 / float(pss.period)
        i = int(harmonic)
        if i < 1:
            raise ValueError('PAC.phase_psd: harmonic must be >= 1.')
        offs = np.atleast_1d(np.asarray(offsets, dtype=float))
        if np.any(offs <= 0.0):
            raise ValueError(
                'PAC.phase_psd: offsets must be positive; S_phi diverges '
                'at zero offset and that divergence is physical.')
        cres = self.coloured_diffusion_resolved(
            pss, offs, frequency_aware=frequency_aware)
        ## ⚠ THE CORNER IS THE WHITE LORENTZIAN'S, read at the carrier.  For
        ## a coloured source `f_h = pi i^2 f0^2 c` is not a lineshape
        ## parameter at all -- there is no Lorentzian -- and taking the
        ## folded value nearest the carrier instead would put a 1/f source's
        ## corner ABOVE the offsets, in front of the power bound below, which
        ## is the floor that actually binds for colour.
        c = self._white_diffusion_at(pss, 2.0 * np.pi * f0)
        ## The i-th harmonic's Lorentzian half-width.  `S_i(f) =
        ## i^2 f0^2 c / (pi^2 i^4 f0^4 c^2 + f^2)` is a Lorentzian in `f`
        ## whose denominator is `f_h^2 + f^2`, so `f_h = pi i^2 f0^2 c`.
        corner = np.pi * (i ** 2) * (f0 ** 2) * c
        if offs.min() <= corner:
            raise ValueError(
                'PAC.phase_psd: offset %.6g Hz is at or below the '
                'Lorentzian corner %.6g Hz for harmonic %d, where S_phi is '
                'not the right object -- the excess phase is a Wiener '
                'process and its spectrum is singular at the origin. The '
                'finite value the LINESHAPE attains there comes from the '
                'nonlinear phase-to-voltage map: use oscillator_spectrum().'
                % (float(offs.min()), corner, i))
        sphi = (i ** 2) * (f0 ** 2) * cres / offs ** 2

        ## ⚠ POWER CONSERVATION AS A SECOND, INDEPENDENT FLOOR -- and for a
        ## COLOURED source it is the binding one, by orders.  The
        ## normalised lineshape integrates to 1, and the integral over one
        ## box of width `df` on each side is a lower bound on it, so
        ##
        ##     2 df S_phi(df) <= 1
        ##
        ## is NECESSARY for the linearised skirt to be consistent with
        ## unit power.  Vanassche, Gielen & Sansen (2003) derive the same
        ## statement for a 1/f input as `df_c >= eps f0 sqrt(2 f_1f)`; the
        ## form here needs no assumption about the source's colour and
        ## reproduces their worked example exactly.
        ##
        ## ⚠ THE LORENTZIAN CORNER ABOVE DOES NOT CATCH THIS.  It is built
        ## from `c` alone, so it knows nothing about a `Gamma(f)` that
        ## grows as the offset falls (on this class's own flicker fixture
        ## the power bound bites 306x above the Lorentzian corner).
        ##
        ## ⚠ IT IS A LOWER BOUND ON THE BREAKDOWN, NOT THE BREAKDOWN.
        ## Passing it is not a guarantee (on Vanassche's own example the
        ## observed flattening sits at 3x the bound): this refuses what is
        ## definitely invalid and admits a band that is already suspect --
        ## deliberately, because refusing at 3x would be fitting a
        ## threshold to one example.
        ## ⚠ AND THE DERIVATION HAS A PRECONDITION THE BOUND DOES NOT
        ## STATE, so it is checked rather than assumed.  The box argument
        ## is `2 df S(df) <= integral_{-df}^{+df} S <= 1`, and the FIRST
        ## inequality needs `S(f) >= S(df)` for every `|f| <= df` -- the
        ## spectrum must not dip below its edge value anywhere further in.
        ## True of a monotone skirt, of the flattened near-carrier shape,
        ## and even with a spur, which ADDS power inside.
        ##
        ## ⚠ FALSE FOR A LOCKED PLL, whose phase-noise transfer function
        ## is HIGH-PASS: the spectrum dips below its edge value everywhere
        ## inside the loop bandwidth.  The bound is then NO LONGER DERIVED,
        ## and a floor that is not derived cannot be used as one.  This
        ## method refuses a driven circuit today; the check is there for the
        ## driven-oscillator work.
        probe = np.unique(np.concatenate((
            offs, np.logspace(np.log10(offs.min() / 1e3),
                              np.log10(offs.max()), 32))))
        ## the probe on the DC PPV: it over-states behind a slow node, so the
        ## bound it serves stays conservative, and it costs no solves
        sprobe = ((i ** 2) * (f0 ** 2)
                  * self.coloured_diffusion_resolved(
                      pss, probe, frequency_aware=False) / probe ** 2)
        if np.any(np.diff(sprobe) > 1e-12 * np.abs(sprobe[:-1])):
            k = int(np.argmax(np.diff(sprobe) > 0)) + 1
            raise ValueError(
                'PAC.phase_psd: the spectrum RISES with offset near '
                '%.6g Hz, so it dips below its edge value further in and '
                'the power bound below is no longer derived -- its box '
                'argument needs S(f) >= S(df) for every |f| <= df. That '
                'happens for a high-pass-shaped spectrum such as a locked '
                'loop, and for a source whose density grows faster than '
                'f^2. The bound may still hold; it is not established '
                'here, so it is refused rather than applied.'
                % float(probe[k]))

        power = 2.0 * offs * sphi
        bad = power >= 1.0
        if np.any(bad):
            k = int(np.argmax(bad))
            raise ValueError(
                'PAC.phase_psd: at offset %.6g Hz the linearised skirt '
                'already carries %.3f times the TOTAL power of the '
                'carrier (2 f S_phi >= 1), so it has broken down there -- '
                'a normalised spectrum integrates to 1. This bound is '
                'independent of the Lorentzian corner (%.6g Hz here) and '
                'for a coloured source it binds far earlier, because '
                'Gamma(f) grows as the offset falls. Sweep above it, or '
                'use oscillator_spectrum() for the lineshape. Note the '
                'TRUE breakdown is higher still: this is a lower bound.'
                % (float(offs[k]), float(power[k]), corner))
        return sphi

    @staticmethod
    def lorentzian(offsets, c, f0, harmonic=1):
        """The `i`-th harmonic's normalised lineshape at `offsets` from it.

            S_i(f) = i² f₀² c / (π² i⁴ f₀⁴ c² + f²)

        `offsets` of either sign; even in them.

        ⚠ EXACT FOR WHITE SOURCES, not a limiting form.  With coloured
        sources the transform "does not have a simple closed form" and only
        two-regime approximations exist — which is why the diffusion
        constants that feed it (`diffusion_constant`,
        `frequency_aware_diffusion`) refuse a coloured source
        (`_refuse_coloured`).

        ⚠ AND ITS TOTAL POWER IS EXACTLY 1.  `∫ a/(b²+f²) df = aπ/b`, and
        here `a = i² f₀² c`, `b = π i⁴ f₀⁴ c² ^ ½`… concretely `b = π i²
        f₀² c`, so the integral is exactly one.  **The carrier's power is
        redistributed, never created or destroyed** — which is the
        invariant that separates this from LTV small-signal treatments,
        which "erroneously predict infinite noise power [at the carrier] as
        well as infinite total integrated power".  It is asserted in the
        suite.

        The half-width is `π i² f₀² c` and the peak `1/(π² i² f₀² c)`, so a
        higher harmonic has a skirt scaling as `i²` and a corner as `i⁴` —
        `20 log₁₀(i)` dB noisier far out.
        """
        i = int(harmonic)
        if i == 0:
            return np.zeros_like(np.asarray(offsets, dtype=float))
        f = np.asarray(offsets, dtype=float)
        a = (i * i) * f0 * f0 * c
        b = np.pi * (i * i) * f0 * f0 * c
        return a / (b * b + f * f)

    def oscillator_spectrum(self, pss, offsets, output, harmonic=1,
                            frequency_aware=True, offset_fmin=None,
                            offset_fmax=None, all_orders=None):
        """Free-running output spectrum at `offsets` from harmonic `harmonic`.

        `offsets` of either sign, and the spectrum is EVEN in them (measured
        equal at +-o on an asymmetric orbit): the lineshape is symmetric
        about the harmonic.  (The modal family reads a negative offset as
        the LOWER sideband, which is not the upper's.)

        ⚠⚠ THIS DOES NOT GO THROUGH `pnoise`'s SIDEBAND FOLD, AND IT CANNOT.
        The fold is a FREQUENCY-CONVERSION computation, complete for a
        driven circuit and structurally incomplete for an AUTONOMOUS
        oscillator: what it omits is exactly the near-carrier phase-noise
        skirt this method returns.  Rizzoli, Mastri & Masotti (IEEE MTT
        42-807, 1994): "frequency-conversion techniques alone are not
        sufficient ... for general autonomous circuits", because the
        noise-induced FREQUENCY MODULATION OF THE CARRIER at low offsets is
        not a frequency-conversion effect.  Their Section III names the two
        stacks -- CONVERSION noise, rising as 1/f for f -> 0, and MODULATION
        noise, "a jitter of the oscillatory steady state", rising as 1/f^3
        -- which DECOUPLE exactly at the steady state, and are "usually
        nearly equal" in an intermediate offset band (a cross-stack
        agreement test not built here).  Diagnostic value: a FLAT PSD near
        the carrier is neither slope -- it is the Phi(T) - I singularity.

        So the Floquet/PPV stack (`ppv`, `diffusion_constant`, this method)
        and the sideband fold (`pnoise`) ARE NOT TWO IMPLEMENTATIONS OF ONE
        QUANTITY, and unifying them is not a simplification waiting to be
        made.  ⚠ THE HAZARD IS THAT THE WRONG ONE STILL RETURNS A NUMBER:
        oscillator phase noise from the fold alone is a spectrum missing the
        dominant contribution near the carrier.  That is the completeness
        argument for the split; the efficiency argument (Floquet is cheaper)
        is the weaker one.

        Returns `(S_v, L_dBc)`.  `S_v` is the ONE-SIDED PSD of the output
        voltage, as `pnoise`'s: the Lorentzian lineshape times the carrier's
        one-sided power `2 |X_1|^2 = A^2/2` (`X_1 = A/2` the carrier phasor).
        Against a reference simulator at every offset over four decades.
        `L_dBc` is `S_v` over that carrier power, in dBc/Hz.  `output`: a
        reduced index, a weight vector, or a node name.
        History: `doc/shooting_history.md`, `PAC.oscillator_spectrum`.

        ⚠ NO SWEEP AND NO PER-FREQUENCY SOLVE.  Once the PSS waveform's
        Fourier coefficients and the scalar `c` are known, "we have an
        analytical expression that gives us the spectrum at any frequency.
        The computation of the spectrum is not performed separately for
        every frequency of interest."  Which also means it never meets the
        near-carrier singularity that a swept small-signal computation
        would, and never meets the 1/f sweep-grid trap — there is no sweep
        to place a point on.

        ⚠⚠ SCOPE: A SOURCE BEHIND A SLOW NODE.  The DC-PPV Lorentzian, for
        a noise source that reaches the core through a slow path (RC leg,
        tau >> T), holds only BELOW the source's corner `T/(2 pi tau)`;
        above it the true skirt is scaled by the PPV-harmonic-weighted
        filter `sum_k |G_k|^2 F_k(f) / sum_k |G_k|^2 F_k(0)` (G_k the PPV
        entry's Fourier coefficients at the source node, F_k the path's
        transfer at k f0 + f) -- 1/1000 at 0.1 f0 on a one-RC-leg fixture.
        `c` is still right, and a Monte Carlo of `c` cannot see it.
        `frequency_aware=True` (the default) replaces `c` by `c(f)` from the
        frequency-aware PPV (`frequency_aware_diffusion`), which matches
        pnoise through the corner (test ..._behind_a_slow_node_...) and to
        <= 2 % above f_amp on an orbit with AM-to-PM coupling.
        `frequency_aware=False` is the closed form, one `c` for every offset,
        and costs no solve.

        ⚠ AND IT IS THE ONLY ROUTE THAT IS VALID BELOW THE CORNER.  A
        small-signal analysis cannot produce `L(f)` there however well
        conditioned it is: the excess phase is a Wiener process, its
        spectrum has a singularity at the origin and no physical meaning,
        and the finite value `L` attains comes from the NONLINEAR
        phase-to-voltage map — which is what this closed form carries.
        Reporting `S_phi` near the carrier instead is the mistake that
        object invites.

        ⚠ A COLOURED SOURCE: with a 1/f source the phase is
        no longer a Wiener process and the line is not a Lorentzian.  The
        lineshape is then the Fourier transform of `exp(-D(tau)/2)`, with
        `D` the excess phase's structure function built from `phase_psd`'s
        `c(f)` (`_lineshape`): the WHITE part of the sources gives the
        Lorentzian in closed form, the coloured part `c(f) - c_w` is
        integrated over `[offset_fmin, offset_fmax]`.  `offset_fmin` is REQUIRED -- a 1/f^3
        phase has no stationary lineshape without a low cutoff; it plays
        the part of the observation time -- and `offset_fmax` defaults to f0/2,
        the phase model's reach.  A white-only circuit is untouched.
        `frequency_aware` (default True) takes `c_fa(nu)`
        from the frequency-aware PPV: to first order in the change (the
        skirt's `i^2 f0^2 (c_fa - c_dc) / f^2` per offset) while that order's
        estimated error is below `FA_FIRST_ORDER_TOL`, and to ALL orders
        above it (the correction's own structure function inside `D`,
        ~25 solves; `_fa_lineshape`, `self.lineshape_info` says which);
        `False` is the DC lineshape.  Each offset takes the transform or the
        far skirt -- `phase_psd`'s linear skirt plus its SECOND-order
        correction (`_lineshape.second_order_skirt`) -- whichever estimates
        the smaller error (`_lineshape.handover`).  Against an mpmath
        reference (`benchmarks/lineshape_reference.py`) the worst is 3.4e-7
        on a real oscillator's line and 4.5e-6 on a very broad one.

        `all_orders`: the frequency-aware correction to ALL
        orders -- `c_fa(nu)` inside the structure function `D`, not per
        offset.  A white source's default is the Lorentzian with `c(f)` per
        offset, first order in the same sense (the core keeps its DC weight,
        `exp(-D_corr(inf)/2)` off; ~8 % with a slow corner 10 linewidths
        out); `all_orders=True` builds the full line instead, `offset_fmax`
        (default f0/2) bounding the band, past which the white part is held
        at its corrected level.  It is OFF by default for a white source
        (~25 bordered solves).  For a
        coloured source None is the estimate's choice (above), True forces
        all orders and False the first.  `self.lineshape_info` says which
        ran.

        History: `doc/shooting_history.md`, `PAC.oscillator_spectrum`.
        """
        output = output_index(pss, output)
        fmin, fmax = offset_fmin, offset_fmax     # (the names the internals use)
        if all_orders and frequency_aware is False:
            raise ValueError(
                'PAC.oscillator_spectrum: all_orders=True is the frequency-'
                'aware correction to all orders, and frequency_aware=False '
                'turns the frequency-aware PPV off; the pair was taken as '
                'frequency_aware=False silently until 2026-09-29.')
        ## (the lineshape is built on the carrier PHASOR's square, `S_v / 2`;
        ## the one-sided PSD doubles it, exactly, and leaves `L_dBc` as is)
        if self._coloured_present(pss):
            Sv, L = self._coloured_spectrum(
                pss, offsets, output, harmonic, fmin, fmax,
                frequency_aware=(frequency_aware is None
                                 or bool(frequency_aware)),
                all_orders=all_orders)
            return 2.0 * np.asarray(Sv), L
        if frequency_aware is None:
            frequency_aware = True
        c = self.diffusion_constant(pss)
        f0 = 1.0 / float(pss.period)
        self._warn_above_amplitude_pole(offsets, f0)
        X = self._carrier_line(pss, output, harmonic)
        self.lineshape_info = {'frequency_aware': None}
        if frequency_aware and all_orders:
            Sv = abs(X) ** 2 * self._white_all_orders(pss, offsets, harmonic,
                                                      c, f0, fmax)
        elif frequency_aware:
            self.lineshape_info = {'frequency_aware': 'first order',
                                   'estimate': None}
            ## `c(f)` per offset -- one bordered adjoint solve each, cached on
            ## `|offset|` so a symmetric sweep pays once per magnitude.  `c(0)`
            ## is `c` exactly, so the near-carrier lineshape is unchanged.
            off = np.asarray(offsets, dtype=float)
            _cache = {}
            flat = []
            for o in np.atleast_1d(off).ravel():
                key = abs(float(o))
                if key not in _cache:
                    _cache[key] = (c if key == 0.0 else
                                   self.frequency_aware_diffusion(pss, key))
                flat.append(float(self.lorentzian(np.array([o]), _cache[key],
                                                  f0, harmonic)[0]))
            Sv = abs(X) ** 2 * np.asarray(flat).reshape(np.shape(off))
        else:
            Sv = abs(X) ** 2 * self.lorentzian(offsets, c, f0, harmonic)
        with np.errstate(divide='ignore'):
            L = 10.0 * np.log10(np.maximum(Sv / max(abs(X) ** 2, 1e-300),
                                           1e-300))
        return 2.0 * Sv, L

    def _white_all_orders(self, pss, offsets, harmonic, c, f0, fmax):
        """The WHITE line with the frequency-aware PPV to all orders
        (`oscillator_spectrum(all_orders=True)`): the Lorentzian's
        `exp(-a |tau|)` becomes `exp(-D/2)` with `D` built from `c_fa(nu) = c
        (1 + rho(nu))`, `rho` from `frequency_aware_diffusion` (one bordered
        solve a frequency, `c(0) = c` exactly).  The machinery is the
        coloured line's (`_fa_orders`) with no colour: `rho` as a Chebyshev
        series from where it dies away up to `fmax`, the change in `D` by
        signed quadrature, the white part held at its corrected level past
        `fmax`, and each offset from the transform or the linear skirt by
        estimated error.  Normalised, as `lorentzian`.

        History: `doc/shooting_history.md`, `PAC._white_all_orders`."""
        i = int(harmonic)
        fmax = 0.5 * f0 if fmax is None else float(fmax)
        if not (0.0 < fmax <= 0.5 * f0 * (1.0 + 1e-12)):
            raise ValueError(
                'PAC.oscillator_spectrum: need 0 < offset_fmax <= f0/2 '
                '(%.6g Hz); got offset_fmax=%r.' % (0.5 * f0, fmax))
        a = 2.0 * np.pi ** 2 * i * i * f0 * f0 * c
        pref = 4.0 * i * i * f0 * f0
        rho_at = lambda nu: np.array(
            [self.frequency_aware_diffusion(pss, nu) / c - 1.0, 0.0])
        off = np.asarray(offsets, dtype=float)
        vals, errs, _shapes = self._fa_orders(
            rho_at, None, None, off, None, i, f0, c, a, pref, None, fmax, True)
        worst = max(errs) if errs else 0.0
        self._warn_lineshape('all-orders', worst, errs, off)
        return np.asarray(vals).reshape(np.shape(off))

    def _warn_lineshape(self, which, worst, errs, off):
        """Warn when the `which` lineshape's estimated error passes
        `LINESHAPE_WARN`, at the offset where it is worst.  (`stacklevel`
        4: the user's call, through `oscillator_spectrum` and its path.)"""
        if worst > self.LINESHAPE_WARN:
            warnings.warn(
                'PAC.oscillator_spectrum: the %s lineshape carries an '
                'estimated relative error of %.1e at offset %.6g Hz (neither '
                'the transform nor the linear skirt is resolved better '
                'there).' % (which, worst, float(np.atleast_1d(off).ravel()[
                    int(np.argmax(errs))])), RuntimeWarning, stacklevel=4)

    def _carrier_line(self, pss, output, harmonic):
        """The carrier phasor `X` at `harmonic`, refusing an output with no
        line there."""
        X = self.carrier_phasor(pss, output, harmonic)
        ## ⚠⚠ NO CARRIER, NO LINE -- AND THE ANSWER WOULD BE A PLAUSIBLE ZERO.
        ## This is a LINE-SHAPE model: it broadens the carrier's own harmonic.
        ## Where the output has no component at `harmonic` (a half-wave
        ## symmetric orbit's even harmonics, or `harmonic = 0`, which the
        ## Lorentzian returns as zeros by construction) there is nothing to
        ## broaden, and the true density is BROADBAND noise this method does
        ## not represent (~0 against 5.6e-6 V^2/Hz at 2 f0 on a symmetric
        ## van der Pol).  Refused, as `am_pm` refuses the same case.  ⚠ NOT
        ## caught: away from the fundamental the model also misses where a
        ## line DOES exist (asymmetric orbit, 2 f0: 0.40 of a Monte Carlo)
        ## -- use `pnoise` away from the fundamental.
        _scale = float(np.max(np.abs(self._output_waveform_row(pss, output))))
        if int(harmonic) == 0 or abs(X) <= 1e-9 * max(_scale, 1e-300):
            raise ValueError(
                'PAC.oscillator_spectrum: the output carries no component at '
                'harmonic %d (|X| = %.3e against a signal scale of %.3e), so '
                'there is no line for this line-shape model to broaden; the '
                'noise there is broadband and this method would return ~0 '
                '(measured: ~0 against 5.6e-6 V^2/Hz on a symmetric van der '
                'Pol at 2 f0). Use PAC.pnoise at that frequency.'
                % (int(harmonic), abs(X), _scale))
        return X

    #: the frequency-aware coloured lineshape stays FIRST ORDER in the
    #: correction while that order's estimated error is below this, and is
    #: taken to all orders above it (`_fa_lineshape`): tighter bought ~1e-6
    #: for 4x the time on an LC with no slow node
    #: History: `doc/shooting_history.md`, `PAC.FA_FIRST_ORDER_TOL`.
    FA_FIRST_ORDER_TOL = 1e-4
    #: the frequency-aware correction `rho = c_fa/c_dc - 1` is probed a decade
    #: at a time down from `fmax` until it falls below this
    FA_RHO_FLOOR = 1e-9
    #: how the all-orders lineshape represents `rho`: 'rational'
    #: (`_lineshape.RationalRho`, ~25 solves, falling back to the Chebyshev
    #: series when a guard refuses it) or 'chebyshev' (`LogChebyshev`, ~65)
    FA_RHO_FIT = 'rational'

    #: the coloured lineshape warns when its estimated error at an offset
    #: exceeds this, relative; the true worst is 3.4e-7 / 4.5e-6 on a real /
    #: a very broad line (`_lineshape.handover`)
    #: History: `doc/shooting_history.md`, `PAC.LINESHAPE_WARN`.
    LINESHAPE_WARN = 1e-3

    def _coloured_spectrum(self, pss, offsets, output, harmonic, fmin, fmax,
                           frequency_aware=True, all_orders=None):
        """`oscillator_spectrum` for a coloured source -- see there and
        `_lineshape`."""
        from . import _lineshape
        self._check_circuit(pss)
        self._refuse_driven(pss, 'oscillator_spectrum')
        f0 = 1.0 / float(pss.period)
        i = int(harmonic)
        ## (a TypeError, a required argument missing, as `_coloured_prepare`
        ## and the sampled family; NotImplementedError until 2026-09-29)
        if fmin is None:
            raise TypeError(
                'PAC.oscillator_spectrum: a noise source in this circuit is '
                'COLOURED, and a 1/f^3 phase has no stationary lineshape '
                'without a low cutoff -- pass offset_fmin (the reciprocal of the '
                'observation time; the coloured part is integrated over '
                '[offset_fmin, offset_fmax], offset_fmax defaulting to f0/2 '
                '= %.6g Hz).'
                % (0.5 * f0))
        fmin = float(fmin)
        fmax = 0.5 * f0 if fmax is None else float(fmax)
        if not (0.0 < fmin < fmax <= 0.5 * f0 * (1.0 + 1e-12)):
            raise ValueError(
                'PAC.oscillator_spectrum: need 0 < offset_fmin < offset_fmax '
                '<= f0/2 (%.6g Hz); got offset_fmin=%r, offset_fmax=%r.'
                % (0.5 * f0, fmin, fmax))
        self._warn_above_amplitude_pole(offsets, f0)
        X = self._carrier_line(pss, output, i)
        ## the white part's `c` and the coloured part of `c(f)` from one fold
        ## (`_colour_fold`: every coloured component through its own
        ## columns, a stationary one the case whose columns do not move)
        fold = self._colour_fold(pss, fmin, None, 'oscillator_spectrum')
        c_w = fold.c_white
        cfun = lambda nus: np.maximum(fold.coloured(nus), 1e-300)
        pc, converged = _lineshape.refine(cfun, fmin, fmax)
        if not converged:
            warnings.warn(
                'PAC.oscillator_spectrum: the coloured c(f) was not resolved '
                'to %.0e at every node within %d nodes; the lineshape is '
                'accurate to about the largest midpoint mismatch.'
                % (_lineshape.NODE_TOL, _lineshape.MAX_NODES),
                RuntimeWarning, stacklevel=3)
        a = 2.0 * np.pi ** 2 * i * i * f0 * f0 * c_w
        pref = 4.0 * i * i * f0 * f0
        ## ⚠ THE FAR SKIRT IS A SKIRT, AND THE TRANSFORM CANNOT SAY SO: there
        ## the lineshape is ~1e-7 of the integrand's scale, and the transform
        ## finds it by cancellation.  So each offset takes the transform or
        ## the second-order skirt, whichever carries the smaller estimated
        ## error (`_lineshape.handover`: the transform's
        ## estimate is its grid-phase move, QUADPACK's bound only where it
        ## failed; the skirt's is its split move plus its correction
        ## squared).
        ## ⚠ THE TWO GRIDS ARE ONE DENSITY AT TWO PHASES: the returned grid
        ## (2 x TAU_PER_DECADE) against the same density with its interior
        ## nodes half a step along, so the estimate reads the returned
        ## value's own error (1.2 .. 2.1x it).  Against HALF the density it
        ## would measure the COARSE grid's.
        ## History: `doc/shooting_history.md`, `PAC._coloured_spectrum`.
        off = np.asarray(offsets, dtype=float)

        def dc():
            ## the DC line; built only where it is used: the all-orders path
            ## replaces it (and it costs a quarter of that path's time)
            ## History: `doc/shooting_history.md`, `PAC._coloured_spectrum`.
            shapes = [_lineshape.ColouredLineshape(
                a, pc, pref, per_decade=2 * _lineshape.TAU_PER_DECADE, shift=sh)
                for sh in (True, False)]
            ## the second-order skirt's `c(nu)`: white plus the colour on its
            ## band, tabulated once from below fmin (the table holds its end
            ## values, and below fmin `c` is the white level)
            flat = np.atleast_1d(off).ravel()
            top = max(1e3 * float(np.max(np.abs(flat))) if flat.size else 0.0,
                      10.0 * fmax)
            ctab = _lineshape.ClampedTable(
                lambda v: c_w + np.where((v >= fmin) & (v <= fmax),
                                         pc(np.clip(v, fmin, fmax)), 0.0),
                1e-3 * fmin, top, per_decade=400)

            def skirt(o):
                ao = abs(float(o))
                cf = c_w + (float(fold.coloured([ao])[0]) if fmin <= ao <= fmax
                            else 0.0)
                return tuple(_lineshape.second_order_skirt(
                    ao, ctab, i * i * f0 * f0, 1e-3 * fmin, fmax, c_at_f=cf,
                    split=sp)
                    for sp in (_lineshape.SKIRT_SPLIT,
                               _lineshape.SKIRT_SPLIT_ALT)) + (
                    i * i * f0 * f0 * cf / (ao * ao),)
            vals, errs = _lineshape.handover(shapes, off, skirt)
            return vals, errs, shapes
        self.lineshape_info = {'frequency_aware': None}
        if frequency_aware:
            vals, errs, shapes = self._fa_lineshape(
                pss, fold, pc, off, dc, i, f0, c_w, a, pref, fmin, fmax,
                mode=all_orders)
        else:
            vals, errs, shapes = dc()
        S = np.asarray(vals).reshape(np.shape(off))
        worst = max(errs) if errs else 0.0
        self._warn_lineshape('coloured', worst, errs, off)
        shape = shapes[1]
        if c_w == 0.0 and shape.line_weight > 1e-12:
            warnings.warn(
                'PAC.oscillator_spectrum: no WHITE source broadens the line, '
                'so a coherent carrier of weight %.3e remains (exp(-D/2) at '
                'infinite lag, set by offset_fmin); it is not in the returned '
                'density.' % shape.line_weight, RuntimeWarning, stacklevel=3)
        Sv = abs(X) ** 2 * S
        with np.errstate(divide='ignore'):
            L = 10.0 * np.log10(np.maximum(S, 1e-300))
        return Sv, L

    def _fa_lineshape(self, pss, fold, pc, off, dc, i, f0, c_w, a, pref,
                      fmin, fmax, mode=None):
        """The frequency-aware PPV in the coloured lineshape (the default).
        The fold's frequency-aware samples (one bordered solve a
        frequency) give the WHITE and the COLOURED parts of `c_fa(nu)`
        separately: `c_w (1 + rho_w)` at every `nu`, and `c_c (1 + rho_c)` on
        `[fmin, fmax]` only, where the colour lives.  ⚠ The white part's
        filtering is physical below `fmin`, which is the colour's
        observation-time cutoff and not the white's: one ratio cut at `fmin`
        dropped `rho = -3.7e-5` there, ~1e-4 at the core (measured).

        FIRST ORDER while it is enough.  Each offset gains `i^2 f0^2 Delta
        / f^2`, the change of the linear skirt: exact past the core, where
        the linear skirt IS `c_fa`.  The carrier is left alone.  Its error
        is ESTIMATED as the larger of two terms:
          * `|D_corr(inf)| / 2`, the core's weight, which the first order
            leaves at the DC one.  `D_corr(inf) = 4 i^2 f0^2 int Delta /
            nu^2`, by the trapezoid in `ln nu` over one probe a decade,
            down from `fmax` until both ratios fall below `FA_RHO_FLOOR`;
            within 6 % of the full value, measured;
          * `|Delta / c_dc| linear_error(f)` per offset, the core's spread
            acting on the change.

        TO ALL ORDERS above `FA_FIRST_ORDER_TOL`.  Both ratios as one
        rational over the probed span (`_lineshape.RationalRho`, ~25
        solves, `FA_RHO_FIT`; the Chebyshev series in `ln nu` if a guard
        refuses it), and each part's `Delta` in
        `D` by signed quadrature on its own band (`_lineshape.SignedTable`,
        `correction_structure`); the transform as before.
        ⚠ Measured with the slow corner 10 linewidths out: the first order
        was 8 .. 13 % low through the core.  Against a closed form
        (`_lineshape`'s test) it went NEGATIVE, -4.55 of the true value,
        where this path stays within 4.8e-7.  On the slow-node fixture the
        corner is 7e6 linewidths out, the estimate ~1e-7, and the first
        order stands.

        `self.lineshape_info` says which, and why.  `mode`: None the
        estimate's choice, True all orders, False first order.

        History: `doc/shooting_history.md`, `PAC._fa_lineshape`."""
        colour = lambda nus: np.asarray(
            fold.coloured(np.atleast_1d(np.asarray(nus, dtype=float))), dtype=float)

        def rho(nu):
            ## `(rho_w, rho_c)` at `nu`: each part of `c_fa` over its DC one
            wf, cf = fold.resolved_fa([nu], parts=True)
            cc = float(colour(nu)[0])
            return np.array([
                float(wf[0]) / c_w - 1.0 if c_w > 0.0 else 0.0,
                float(cf[0]) / cc - 1.0 if cc > 0.0 else 0.0])
        return self._fa_orders(rho, colour, pc, off, dc, i, f0, c_w, a, pref,
                               fmin, fmax, mode)

    def _fa_orders(self, rho_at, colour, pc, off, dc, i, f0, c_w, a, pref,
                   fmin, fmax, mode):
        """`_fa_lineshape`'s choice and both its paths, for a coloured line
        (`colour`, the DC coloured `c` on `[fmin, fmax]`) or a white one
        (`colour` None, `pc` None, `dc` None: `_white_all_orders`).
        `rho_at(nu)`: `(rho_w, rho_c)`, one bordered solve; `dc()` the DC
        line `(vals, errs, shapes)`, which the first order corrects."""
        from . import _lineshape
        flat = np.atleast_1d(off).ravel()
        if colour is None:
            colour = lambda nus: np.zeros(np.shape(np.atleast_1d(nus)))
            fmin = fmax
        rhos = {}

        def rho(nu):
            nu = float(nu)
            if nu not in rhos:
                rhos[nu] = np.asarray(rho_at(nu), dtype=float)
            return rhos[nu]

        def delta(nu):
            ## the change of `c` at `nu`: the colour's part on its band only
            r = rho(nu)
            d = c_w * r[0]
            if nu >= fmin:
                d += float(colour(nu)[0]) * r[1]
            return d
        if mode is False:
            vals, errs, shapes = dc()
            for k, o in enumerate(flat):
                if o != 0.0:
                    ao = abs(float(o))
                    vals[k] += i * i * f0 * f0 * delta(ao) / ao ** 2
            self.lineshape_info = {'frequency_aware': 'first order',
                                   'estimate': None, 'solves': len(rhos)}
            return vals, errs, shapes
        ## one probe a decade down from fmax, to where both parts have died
        ## away (the white one below fmin too); bounded at 20 decades below
        probes, nu, floored = [], fmax, False
        for _k in range(int(np.ceil(np.log10(fmax / fmin))) + 20):
            probes.append(nu)
            r = rho(nu)
            if (abs(r[0]) < self.FA_RHO_FLOOR
                    and (nu < fmin or abs(r[1]) < self.FA_RHO_FLOOR)):
                floored = True
                break
            nu /= 10.0
        pn = np.asarray(sorted(probes))
        ## (with the white part's tail past fmax, held at its level there)
        cinf = pref * (np.log(10.0) * float(sum(delta(v) / v for v in pn))
                       + c_w * rho(fmax)[0] / fmax)
        est = 0.5 * abs(cinf)
        for o in flat:
            if o != 0.0:
                ao = abs(float(o))
                cdc = c_w + (float(colour(ao)[0]) if ao >= fmin else 0.0)
                est = max(est, abs(delta(ao)) / cdc * min(
                    _lineshape.linear_error(o, i, f0, c_w, pc), 1.0))
        info = {'estimate': est, 'cinf_estimate': cinf,
                'probes': len(probes)}
        if mode is None and est <= self.FA_FIRST_ORDER_TOL:
            vals, errs, shapes = dc()
            for k, o in enumerate(flat):
                if o != 0.0:
                    ao = abs(float(o))
                    vals[k] += i * i * f0 * f0 * delta(ao) / ao ** 2
            info.update(frequency_aware='first order', solves=len(rhos))
            self.lineshape_info = info
            return vals, errs, shapes
        ## to all orders
        nu_lo = float(pn[0])
        if not floored:
            warnings.warn(
                'PAC.oscillator_spectrum: the frequency-aware correction had '
                'not died away (below %.0e) at %.6g Hz, the lowest probe; the '
                'part below is left out.' % (self.FA_RHO_FLOOR, nu_lo),
                RuntimeWarning, stacklevel=4)
        ## `rho` as ONE rational from ~25 adaptive solves; any guard that
        ## refuses it (verification, a pole on the band, the budget) falls
        ## back to the Chebyshev series, which reuses every solve made
        cheb, rho_fit = None, 'chebyshev'
        if self.FA_RHO_FIT == 'rational':
            fit = _lineshape.RationalRho(rho, nu_lo, fmax)
            if fit.converged:
                cheb, rho_fit = fit, 'rational'
            else:
                rho_fit = 'chebyshev (the rational fit: %s)' % fit.reason
        if cheb is None:
            ## ⚠ relative to `1 + rho` as the rational fit is (see
            ## `RationalRho`): the smallest the probes saw, floored
            w_min = max(min(float(np.min(np.abs(1.0 + rho(v)))) for v in pn),
                        _lineshape.RationalRho.REL_FLOOR)
            cheb = _lineshape.LogChebyshev(rho, nu_lo, fmax,
                                           tol=_lineshape.CHEB_TOL * w_min)
            if not cheb.converged:
                warnings.warn(
                    'PAC.oscillator_spectrum: the frequency-aware correction '
                    'was resolved to %.1e only (Chebyshev degree %d); the '
                    'lineshape is no better than that.'
                    % (cheb.err, cheb.c.shape[0] - 1), RuntimeWarning,
                    stacklevel=4)
        tabs = []
        if c_w > 0.0:
            tabs.append(_lineshape.SignedTable(
                lambda v: c_w * cheb(v)[0], nu_lo, fmax))
            ## past fmax the white part stays at its frequency-aware level:
            ## returning to `c_w` there is a jump whose ringing swamped `D`
            ## behind a slow node (`_lineshape.ConstantTail`)
            tabs.append(_lineshape.ConstantTail(c_w * rho(fmax)[0], fmax))
        clo = max(nu_lo, fmin)
        if clo < fmax:
            tabs.append(_lineshape.SignedTable(
                lambda v: colour(v) * cheb(v)[1], clo, fmax))
        ## (one density at two phases, as `_coloured_spectrum`'s)
        shapes = [_lineshape.ColouredLineshape(
            a, pc, pref, per_decade=2 * _lineshape.TAU_PER_DECADE, corr=tabs,
            shift=sh) for sh in (True, False)]
        ## the second-order skirt's `c_fa(nu)`: the fit below fmax, the white
        ## part held at its level there above it (as `D` holds it),
        ## tabulated once; the value AT each offset is its exact solve
        rw_top = float(rho(fmax)[0])
        top = max(1e3 * float(np.max(np.abs(flat))) if flat.size else 0.0,
                  10.0 * fmax)
        t_lo = 1e-3 * min(nu_lo, fmin)

        def c_fa(v):
            v = np.atleast_1d(np.asarray(v, dtype=float))
            r = np.asarray(cheb(np.clip(v, nu_lo, fmax)), dtype=float).reshape(2, -1)
            r = np.where(v[None, :] >= nu_lo, r, 0.0)
            out = c_w * (1.0 + np.where(v <= fmax, r[0], rw_top))
            band = (v >= fmin) & (v <= fmax)
            col = np.where(band, colour(np.clip(v, fmin, fmax)), 0.0)
            return out + col * (1.0 + r[1])
        ctab = _lineshape.ClampedTable(c_fa, t_lo, top, per_decade=400)

        def skirt(o):
            ao = abs(float(o))
            cf = (c_w + (float(colour(ao)[0]) if fmin <= ao <= fmax else 0.0)
                  + delta(ao))
            return tuple(_lineshape.second_order_skirt(
                ao, ctab, i * i * f0 * f0, t_lo, fmax, c_at_f=cf, split=sp)
                for sp in (_lineshape.SKIRT_SPLIT, _lineshape.SKIRT_SPLIT_ALT)
                ) + (i * i * f0 * f0 * cf / (ao * ao),)
        vals, errs = _lineshape.handover(shapes, flat, skirt)
        info.update(frequency_aware='all orders', solves=len(rhos),
                    rho_fit=rho_fit, rho_err=cheb.err, cinf=shapes[1].cinf)
        self.lineshape_info = info
        return vals, errs, shapes

    def frequency_aware_diffusion(self, pss, offset):
        """`c(f)` — the phase diffusion constant seen at modulation offset `f`.

        `offset` a scalar of either sign; its magnitude is used.

        `c(f) = (1/T) integral v_f^H (CY/2) v_f dt` with `v_f` the
        frequency-aware PPV (`PSS.frequency_aware_ppv`, Lai 2008 eq. 23) on
        the equation rows, integrated by the SAME quadrature as
        `diffusion_constant` over the orbit's full grid -- so `c(0)` is `c`
        exactly, not approximately.

        ⚠⚠ WHY IT EXISTS.  The Lorentzian from `c` uses the DC PPV at every
        offset: a noise current is assumed to move the phase instantly.
        Wherever part of that response goes THROUGH a slow mode -- the
        amplitude mode on an orbit with AM-to-PM coupling, or a slow node in
        the source's path -- it is filtered above that mode's corner, and
        the DC-PPV Lorentzian over-states.  With `c(f)` the Lorentzian
        matches pnoise's PM content to ~2 % over 0.3-10 f_amp on van der Pol
        with an asymmetric orbit (the DC PPV: up to 3.3x high) and to 0.4 %
        behind a slow node (the DC PPV: 0.73x and 0.027x); on a symmetric
        orbit `c(f)/c` stays within 1e-3.

        ⚠ WHITE sources only, like `diffusion_constant`; a MODULATED one
        is read at each sample's state, as there.  Cost:
        one bordered adjoint GMRES per offset (0.25-0.7 s on these fixtures);
        the solve can fail to converge (Lai's own warning about eq. 23), and
        then this raises rather than returning the DC value silently.

        History: `doc/shooting_history.md`, `PAC.frequency_aware_diffusion`.
        """
        self._check_circuit(pss)
        self._refuse_coloured(
            pss, 'frequency_aware_diffusion',
            'phase_psd and coloured_diffusion_resolved give a coloured '
            "source's phase noise.")
        self._refuse_driven(pss, 'frequency_aware_diffusion')
        off = abs(float(offset))
        if off == 0.0:
            return self.diffusion_constant(pss)
        _v, fi = pss.frequency_aware_ppv(off)
        base = fi['ppv']
        m = pss.cir.n - 1
        S = np.asarray(fi['samples_eq'])[:, :m]
        tms = np.asarray(base['times'], dtype=float)
        n = min(len(tms) - 1, S.shape[0])
        T = float(base['period'])
        ## ⚠ `_period_weights`, as `diffusion_constant` integrates: the LEFT
        ## RECTANGLE (`np.diff(times)`) is first order on a smoothly
        ## non-uniform grid, and `c(f)` would jump at f = 0 instead of
        ## reaching `c`
        ## History: `doc/shooting_history.md`, `PAC.frequency_aware_diffusion`.
        h = self._period_weights(tms, n, T, pss)
        w0 = 2.0 * np.pi / float(pss.period)
        try:
            cy = 0.5 * np.real(self._cy_reduced(pss, w0))
        except NotImplementedError:
            cy = None
        if cy is not None:
            quad = np.real(np.einsum('ij,jk,ik->i', np.conj(S[:n]), cy, S[:n]))
        else:
            ## ⚠ A MODULATED source: `CY` at each sample's own
            ## state, as `diffusion_constant` reads it -- the PPV's orbit,
            ## `_ppv_states`.  ⚠ Not `pss.waveform`: trap and
            ## euler serve `factored_period()` -- and so this PPV -- from
            ## their monodromy TWIN, and the solve's own orbit differs from
            ## the twin's by the discretisation (-1.3e-4 in c(0+)/c on trap).
            ## History: `doc/shooting_history.md`, `PAC.frequency_aware_diffusion`.
            cys = 0.5 * np.real(self._noise_components(
                pss, self._ppv_states(pss)[:n]).cy_at_states(w0))
            quad = np.real(np.einsum('ij,ijk,ik->i', np.conj(S[:n]), cys, S[:n]))
        return float((quad * h[:n]).sum() / T)

    def _warn_above_amplitude_pole(self, offsets, f0):
        """⚠ THE PHASE-ONLY SPECTRUM IS A LOWER BOUND ABOVE `f_amp`.

        `oscillator_spectrum` returns the PHASE contribution only.  A real
        oscillator also carries AMPLITUDE noise, which is suppressed near the
        carrier because the limit cycle restores the amplitude -- but only at
        the amplitude-relaxation rate.  Above the pole where that restoring
        action runs out, amplitude noise stops decaying within a period and
        adds to the total, so this method UNDER-reports (a relayed
        measurement against a commercial simulator's total noise, ~3 dB low
        well above f_amp).

        ⚠⚠ AND THE VALID REGION SHRINKS AS `1/Q`.  With
        `f_amp = -ln(lam2)/(2 pi T)` and `Q = -1/ln(lam2)`,

            f_amp = f0 / (2 pi Q)

        so the better the oscillator, the narrower the band in which its
        phase-only spectrum is the whole answer; at `lam2 = 0.999` it has
        collapsed below ~253 Hz.

        ⚠ THIS IS THE OPPOSITE SIGN FROM THE ERROR `PSS.ppv` ALREADY WARNS
        ABOUT.  That one says the instantaneous phase equation misses slow
        nodes which FILTER device noise, so phase noise is OVER-estimated.
        This one is a second, independent mechanism in which the phase-only
        answer is UNDER-estimated.  Both are live and they are not the same
        effect.

        History: `doc/shooting_history.md`, `PAC._warn_above_amplitude_pole`.
        """
        lam2, certified = getattr(self, '_last_second_multiplier',
                                  (None, None))
        if lam2 is None:
            return
        lam2 = float(lam2)
        ## `lam2 <= 0` is a real or overdamped mode with no relaxation pole to
        ## speak of, and `lam2 >= 1` is not a decaying mode at all -- in both
        ## cases there is no `f_amp` and inventing one would be worse than
        ## silence.
        if not (0.0 < lam2 < 1.0):
            return
        ## `f_amp = -ln(lam2)/(2 pi T)` and `T = 1/f0`.
        f_amp = -np.log(lam2) * float(f0) / (2.0 * np.pi)
        off = np.atleast_1d(np.asarray(offsets, dtype=float))
        worst = float(np.max(np.abs(off))) if off.size else 0.0
        if worst < f_amp:
            return
        warnings.warn(
            'PAC.oscillator_spectrum: this is a PHASE-ONLY spectrum and %g Hz '
            'is above the amplitude-relaxation pole f_amp = %.4g Hz '
            '(lambda_2 = %.6f, f_amp = f0/(2*pi*Q)). Above f_amp the '
            'amplitude noise no longer decays within a period and adds to the '
            'total, so on a half-wave-symmetric orbit the value returned here '
            'is a LOWER BOUND: measured excess of a commercial simulator over '
            'the phase-only prediction is -0.54 dB at 1 kHz and -2.90 dB at '
            '10 kHz for lambda_2 = 0.99. On an ASYMMETRIC orbit it can instead '
            'OVER-state the total (pnoise/phase-only = 0.61 at half-wave '
            'asymmetry 0.10 with frequency_aware=False, confirmed by Monte '
            'Carlo; the default frequency_aware=True corrects the phase part '
            'to ~2 %%) -- use PAC.pnoise for '
            'the total above f_amp. '
            '%sThe valid band scales as 1/Q, so it NARROWS as the oscillator '
            'improves.'
            % (worst, f_amp, lam2,
               ('' if certified is not False else
                'lambda_2 itself is NOT certified here (see '
                "info['second_multiplier_certified']), so f_amp is uncertain "
                'too. ')),
            RuntimeWarning, stacklevel=3)
