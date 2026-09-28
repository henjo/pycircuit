"""Driven-circuit noise: pnoise (stationary, modulated and cyclostationary
folds), its AM/PM split and the band spread.
"""
import numpy as np
import warnings


class _DrivenNoise(object):
    """Driven-circuit noise: pnoise (stationary, modulated and
    cyclostationary folds), its AM/PM split and the band spread.  A theme of
    `PAC` (see `pac.py`)."""

    ## How small a sideband's contribution must be, relative to the running
    ## total, before the accumulation stops.  Okumura et al.: powers are
    ## "accumulated until their contributions become negligible".
    ALIAS_RATIO_TOL = 1e-9

    ## ⚠ THE FOLD BELOW IS FOR DRIVEN CIRCUITS.  It is a frequency-conversion
    ## computation and is complete for one; for an AUTONOMOUS oscillator it is
    ## structurally incomplete -- the near-carrier phase-noise skirt is not a
    ## conversion effect (Rizzoli, Mastri & Masotti, MTT 42-807, 1994).  Free-
    ## running phase noise goes through the Floquet/PPV stack instead; see
    ## `oscillator_spectrum` for why the two cannot be unified.
    def pnoise(self, pss, freq, output, ratio_tol=None, maxsidebands=None,
               modulated=False, cyclostationary=False):
        """TIME-AVERAGED output noise PSD at `freq`, sidebands folded in.

            S(f) = sum_l  h_l CY h_l^H ,   h_l = H_l(f - l f0)

        Noise entering at `f - l f0` leaves at `f` through sideband `l`, and
        white sources in disjoint bands are uncorrelated, so the bands add
        in POWER.  Each `h_l` is one adjoint row -- one transposed solve for
        every source in the circuit, which is the whole reason this is
        affordable.

        Returns `(S, sidebands_used)`.  `S` is the one-sided
        **time-averaged** PSD at the output, in the same units as
        `analysis_ss.Noise`'s `Svnout`.  Gated against `analysis_ss.Noise`
        on a linear circuit, where the sidebands vanish and this reduces
        to the stationary answer (Okumura's `p = 1` case).

        The sources' `CY`, three ways:

          * default: STATIONARY SOURCES ONLY, AND IT CHECKS.  A
            bias-dependent `CY` raises (`_cy_reduced`) rather than return a
            number that is quietly the wrong model.
          * `modulated=True`: Hull & Meyer's route -- one stationary source
            at the CYCLE-AVERAGED `CY`, the modulation carried by `H_l`.  It
            keeps the power and drops the correlation between sidebands.
          * `cyclostationary=True`: the construction for a bias-dependent
            `CY`.  A source whose PSD follows the orbit is white noise
            modulated by `B(t) = sqrt(CY(x(t)))`; its band at `f - p f0`
            reaches the output through every modulation harmonic,
            COHERENTLY over `k` and incoherently over `p`, through the SAME
            rows as the stationary fold.  Summing the bands turns the square
            root into the PSD's own harmonics `P_j` (the DFT of `CY(x(t))`,
            no matrix square root, no window count `p`):

                S(f) = sum_{l,l'} a_l P_{l'-l} a_{l'}^H,

            exact on the grid; constant `CY` collapses it to the stationary
            sum.  A coloured source is folded band by band (see
            `_cyclostationary_fold`).  Like the stationary fold it is a
            LOWER bound at a sideband cap.  Coherence is the whole content:
            it differs from `modulated=True` wherever the modulated noise
            crosses a periodically varying transfer, and equals it (only
            `P_0` survives) through a time-invariant one
            (`test_..._cyclostationary_...`).  Cost: white = the cycle
            average's; coloured (any frequency-dependent `CY`, a negligible
            flicker coefficient included) ~6x.
            ⚠ FLICKER: a PSD cannot carry the modulation's SIGN (for a
            coloured source `m xi` and `|m| xi` are different processes), so
            this fold is the `|m|` one: exact when the modulation is
            sign-definite, different physics when it changes sign --
            Okumura's "cannot be modeled as a cyclostationary process by
            using this method, because it has very long time constants".

        ⚠ "TIME-AVERAGED" IS NOT A HEDGE, IT IS THE SPECIFICATION.  Output
        noise is cyclostationary through bias-dependent sources AND through
        the PERIODIC SOURCE-TO-OUTPUT TRANSFER, which applies even when
        every source is stationary.  The time average is sufficient unless
        something downstream tracks the PSD's variation -- a NONLINEAR
        SUBSEQUENT STAGE, or CASCADED STAGES OFF A SHARED REFERENCE
        (Kundert, *Introduction to RF Simulation*).  One number per output
        frequency cannot carry the correlation between frequencies `k f0`
        apart.  (A commercial RF simulator's PNoise computes the same time
        average.)

        ⚠ `maxsidebands` IS AN ACCURACY KNOB HERE AND A REPORTING KNOB IN
        `PAC.solve`.  Noise lives at every frequency, so capping sidebands
        drops power that belonged in the total: `S` becomes a LOWER bound,
        never a cheaper estimate of the same number.

        ⚠ TWO STOPPING RULES, AND THE BOUND IS NOT THE RATIO TEST.  The
        accumulation stops when a sideband pair adds less than `ratio_tol`
        of the running total -- and it can never pass `|l| <= N/2`, the
        grid's own Nyquist (Okumura eq. 32).  `alias_stop` says which
        fired; ending on the bound warns.

        ⚠ HARMONICS: folding puts a copy of a 1/f source's DC singularity
        at every harmonic.  A `freq` on a harmonic where the folded `CY` is
        non-finite or frequency-dependent raises; one just beside it warns.
        Cluster frequencies NEAR each harmonic, never ON it.

        ⚠ AN OSCILLATOR is not this function's problem: its output noise
        is STATIONARY (Demir 2002; `I - M kron M` is singular,
        `test_no_periodic_covariance_exists_for_an_oscillator`), and its
        phase noise is `oscillator_spectrum`'s.  Near a harmonic the
        operator is singular with the PPV as its null vector, which a plain
        solve shows as "flat PSD curves or curves with unexpected slope
        near the oscillation frequency" (Gourary et al.).  The rows use the
        deflated solve (`_deflated_solve`, Gourary et al. eq. 27/28,
        bordered with BOTH null vectors); `PAC.deflated` says which route
        ran.

        History: `doc/shooting_history.md`, `PAC.pnoise`.
        """
        self._check_circuit(pss)
        ## pnoise folds sidebands through the ADJOINT (adjoint_sideband_row ->
        ## _forced_replay_transposed), whose two-stage chained transpose is
        ## not built for TR-BDF2, so it falls back to a Gear-2 twin -- see
        ## `_adjoint_host`.  (covariance/oscillator_covariance use the built
        ## TR-BDF2 Lyapunov injection via `_lyapunov_host`.)
        pss = pss._adjoint_host()
        fp = pss.factored_period()
        m = pss.cir.n - 1
        N = len(fp.steps)
        T = float(fp.T)
        f0 = 1.0 / T
        tol = self.ALIAS_RATIO_TOL if ratio_tol is None else float(ratio_tol)
        lmax = N // 2 if maxsidebands is None else min(int(maxsidebands),
                                                       N // 2)

        w = 2.0 * np.pi * float(freq)
        ## ⚠ `modulated=True` IS HULL & MEYER'S ROUTE, NOT A TOLERANCE
        ## RELAXATION.  Off, a bias-dependent `CY` raises, because the
        ## stationary sum would be the wrong model.  On, the source is
        ## replaced by ONE stationary source at the CYCLE-AVERAGED bias and
        ## the modulation is carried by `H_l` -- which is the standard
        ## treatment of exactly this case, and the only route to MOS
        ## pnoise, since no physically correct MOS noise model has a
        ## state-independent `CY`.
        colour = None
        if cyclostationary:
            ## the stop rule and the harmonic probes below run on the
            ## cycle-averaged power (the modulation's B_0 B_0^H); the fold
            ## itself is the convolution after the rows are gathered.  The
            ## colour model is fitted ONCE here and serves both, per element
            ## so independent sources ADD in the coloured fold (see
            ## `_cy_components_model`); its call is the summed model the
            ## stop rule reads
            colour = self._cy_components_model(pss, float(freq), f0)
            if colour is None:
                cyfn = self._cy_cycle_averaged
            else:
                ## the cycle average with `_cy_cycle_averaged`'s weights
                ## History: `doc/shooting_history.md`, `PAC.pnoise`.
                fp_ = pss.factored_period()
                tms_ = np.asarray(fp_.times, dtype=float)
                def cyfn(pss_, w_, _m=colour, _t=tms_, _T=float(fp_.T)):
                    Cs = _m(w_)
                    ns = min(len(_t) - 1, Cs.shape[0])
                    _h = self._period_weights(_t, ns, _T, pss_)
                    return np.einsum('k,kij->ij', _h[:ns], Cs[:ns]) / float(_h[:ns].sum())
        else:
            cyfn = (self._cy_cycle_averaged if modulated else self._cy_reduced)
        cy = cyfn(pss, w)

        ## ⚠⚠ ON A HARMONIC, A SIDEBAND FOLDS THE SOURCES TO DC -- AND
        ## SOME DEVICE MODELS ARE NOT DEFINED THERE.  Sideband `l`
        ## evaluates `CY` at `f - l f0`, so `f = k f0` evaluates it at
        ## ZERO: a `1/f` term is infinite there, and a flicker term with its
        ## coefficient set to ZERO is `0/0` = `nan` (`PspMosLongChannel`,
        ## `fnt=1, nfa=0`) -- a caller who sets `nfa = 0` believing flicker
        ## is off still gets `nan`.  White sources fold to DC harmlessly, so
        ## this checks the SOURCES at the frequency actually used rather
        ## than refusing a harmonic on principle.
        f0_ = 1.0 / float(pss.period)
        lscan = max(1, int(maxsidebands or 8))
        offs = np.abs(float(freq) - np.arange(-lscan, lscan + 1) * f0_)
        near = float(np.min(offs))
        if near <= self.HARMONIC_GUARD * f0_:
            probe = cyfn(pss, 2.0 * np.pi * near)
            if not np.all(np.isfinite(np.asarray(probe))):
                raise ValueError(
                    'PAC.pnoise: %.12g Hz sits on a harmonic of %.12g Hz, '
                    'so a sideband folds the noise sources to DC -- and at '
                    'DC this circuit\'s CY is not finite. A 1/f term is '
                    'infinite there; a flicker term whose COEFFICIENT IS '
                    'ZERO is 0/0 and gives nan, so disabling flicker does '
                    'not avoid this. Offset from the harmonic: a commercial RF simulator\'s '
                    'own advice is to cluster frequencies NEAR each '
                    'harmonic and never place one ON it.'
                    % (float(freq), f0_))
            ## ⚠⚠ A FINITE PROBE IS NOT A SAFE ONE.  `1/T` rounds, so at
            ## `f = f0` the folded band sits ~1e-11 Hz from DC, not ON it: a
            ## 1/f source there is finite and enormous (6.3e-2 V^2/Hz against
            ## 9.2e-15 at 0.1 % either side).  So a frequency-DEPENDENT CY at
            ## the folded band refuses too; a white one stays allowed.
            probe_ref = cyfn(pss, 2.0 * np.pi * max(abs(float(freq)) * 2.0, f0_))
            if not np.allclose(np.asarray(probe), np.asarray(probe_ref),
                               rtol=1e-9, atol=0.0):
                raise ValueError(
                    'PAC.pnoise: %.12g Hz sits on harmonic %d of %.12g Hz, so '
                    'a sideband folds the noise sources to %.3g Hz -- DC up to '
                    'rounding -- and this circuit\'s CY is frequency-dependent '
                    'there: a 1/f source is read next to its singularity and '
                    'the fold returns a finite, absurd number (measured '
                    '6.3e-2 V^2/Hz against 9.2e-15 at 0.1 %% either side). '
                    'Offset from the harmonic, or use PAC.sampled_variance, '
                    'whose explicit fmin keeps every band off DC.'
                    % (float(freq), int(round(abs(float(freq)) / f0_)), f0_,
                       near))

        ## ⚠ THE STEEP REGION BESIDE A HARMONIC IS A SWEEP HAZARD RATHER
        ## THAN A WRONG NUMBER, so it warns instead of raising: the VALUE is
        ## right (2 % above the plateau at `f0 + 0.01` Hz with a real flicker
        ## source), but a grid that lands there by accident integrates a
        ## spike it never resolved.
        elif near < 1e-6 * f0_ and float(freq) > 0.0:
            cy_hi = cyfn(pss, 2.0 * np.pi * max(float(freq) * 2.0, f0_))
            if not np.allclose(cy, cy_hi, rtol=1e-9, atol=0.0):
                warnings.warn(
                    'PAC.pnoise: %.12g Hz is %.3g Hz from a harmonic of '
                    '%.12g Hz and a source has a frequency-dependent CY, '
                    'so the folded density varies steeply here. The VALUE '
                    'is correct; a swept grid landing this close will '
                    'misrepresent the integrated total. Cluster near each '
                    'harmonic deliberately rather than by accident.'
                    % (float(freq), near, f0_), RuntimeWarning, stacklevel=2)

        total = 0.0
        used = []
        quiet = 0
        self.alias_stop = 'bound'
        rows = {}
        for l in range(0, lmax + 1):
            step = 0.0
            for sl in ((0,) if l == 0 else (l, -l)):
                fin = float(freq) - sl * f0
                h = self.adjoint_sideband_row(pss, fin, output, sl)[0]
                rows[sl] = np.asarray(h, dtype=complex)
                step += float(np.real(h @ cyfn(
                    pss, 2.0 * np.pi * fin) @ np.conj(h)))
                used.append(sl)
            total += step
            if total > 0 and abs(step) < tol * abs(total):
                ## ⚠ TWO QUIET PAIRS, NOT ONE.  A single sideband can come
                ## back near zero by symmetry while its neighbours do not,
                ## and stopping there would truncate a series that had not
                ## converged.
                quiet += 1
                if quiet >= 2:
                    self.alias_stop = 'ratio'
                    break
            else:
                quiet = 0
        self.sidebands_used = used
        if cyclostationary:
            total = self._cyclostationary_fold(pss, float(freq), rows, model=colour)
        ## ⚠ WHICH RULE STOPPED IT IS PART OF THE ANSWER.  Ending on the
        ## ratio test means the series converged; ending on the Nyquist
        ## bound means the grid ran out before the series did, and the
        ## number is a LOWER bound on the folded noise -- every sideband
        ## above the grid's own maximum frequency is missing, not small.
        ## A strongly switching circuit does this readily.
        if self.alias_stop == 'bound' and lmax > 0:
            warnings.warn(
                'PAC.pnoise: the sideband accumulation stopped at the '
                "grid's Nyquist (|l| = %d at %d points per period), not "
                'because the contributions became negligible. Sidebands '
                'above the grid\'s maximum frequency are MISSING rather '
                'than small, so this is a lower bound on the folded noise. '
                'Re-solve the PSS on a finer period grid and compare.'
                % (lmax, N),
                RuntimeWarning, stacklevel=2)
        return total, used

    def _cy_harmonics(self, pss, w):
        """`P_j`: the Fourier coefficient matrices of `CY(x(t), w)` over the
        orbit, `(N, n, n)` indexed like `numpy.fft.fftfreq`.  `P_0` is the
        cycle average; `P_j` with `j != 0` carry the modulation and vanish
        for a bias-independent source.  No square root: the fold uses the
        PSD's own harmonics (`a P a^H`), exact on the grid at any window.
        (A sqrt-modulation route is not: the root of a PSD that crosses
        zero has a kink and a slowly decaying harmonic tail, which the
        sideband window truncates.)

        History: `doc/shooting_history.md`, `PAC._cy_harmonics`."""
        return self._period_dft(pss, self._cy_at_states(pss, w))

    @staticmethod
    def _sqrt_harmonics_of(Cs, dft=None):
        """The DFT of the symmetric square root of the sampled `CY` (see
        `_cy_sqrt_harmonics`), for a `(N, n, n)` array already in hand."""
        Bs = []
        for cyk in Cs:
            cyk = 0.5 * (cyk + cyk.conj().T)
            lam, U = np.linalg.eigh(cyk)
            lam = np.clip(np.real(lam), 0.0, None)
            Bs.append((U * np.sqrt(lam)[None, :]) @ U.conj().T)
        Bs = np.asarray(Bs, dtype=complex)
        return (np.fft.fft(Bs, axis=0) / Bs.shape[0]) if dft is None else dft(Bs)

    def _cy_sqrt_harmonics(self, pss, w):
        """`B_k`: the DFT of the symmetric square root of `CY(x(t), w)` over
        the orbit, `(N, n, n)`, for the band-resolved (coloured) fold."""
        return self._sqrt_harmonics_of(self._cy_at_states(pss, w),
                                       lambda Bs: self._period_dft(pss, Bs))

    def _cyclostationary_fold(self, pss, freq, rows, model=None):
        """`S(f) = sum_{l,l'} a_l Q_{l,l'} a_{l'}^H` over the gathered
        sideband rows (`rows[l]` = the row for a source at `f - l f0`,
        output at `f`).  WHITE source: `Q_{l,l'} = P_{l'-l}`, the DFT of
        `CY(x(t))` itself -- no square root, exact on the grid.  COLOURED
        source (`CY` depends on `w`; detected by comparing two bands): the
        white band `p = l + k` shared by rows `l` and `l'` carries its OWN
        `CY`, so `Q_{l,l'} = sum_k B_k^{(l+k)} B_{k+l-l'}^{(l+k) H}` with
        `B^{(p)}` the sqrt-DFT at the band's frequency `|f - p f0|`, summed
        over ALL `N` modulation harmonics `k` (which is what makes the
        square root exact here; a window on `k` is not).  ⚠ On a flicker
        source `||P_0||` differs 24x across the bands the fold sums, so
        "the band of l" is not an approximation to use; the band-resolved
        form is pinned against the stationary fold of a stationary FLICKER
        source through the same multiplier.

        The colour model (`_cy_components_model`, fitted once in `pnoise`
        and shared with the stop rule) and a vectorised pair sum keep the
        coloured call near the white one, exact to 1e-11 against the
        per-band evaluation, which remains the fallback for a colour the
        model does not fit.

        History: `doc/shooting_history.md`, `PAC._cyclostationary_fold`."""
        f0 = 1.0 / float(pss.period)
        ls = sorted(rows)
        f = float(freq)
        ## coloured or white?  two bands, same test the stationary path uses
        w_a = 2.0 * np.pi * abs(f - ls[0] * f0)
        w_b = 2.0 * np.pi * max(abs(f) * 2.0, f0)
        Pa = self._cy_harmonics(pss, w_a)
        Pb = self._cy_harmonics(pss, w_b)
        coloured = not np.allclose(Pa, Pb, rtol=1e-9, atol=0.0)
        total = 0.0
        if not coloured:
            P = Pa
            N = P.shape[0]
            for l in ls:
                for lp in ls:
                    total += complex(rows[l] @ P[(lp - l) % N] @ np.conj(rows[lp]))
            return float(np.real(total))
        ## band-resolved: every white band the modulation harmonics reach.
        ## The cost is the circuit's CY, not the algebra: every colour in
        ## the library is thermal plus flicker in 1/f^ef, so the model fixes
        ## each entry's shape from three evaluations per sample (A + B
        ## (w1/w)^ef, ef by a root find), verifies it at further
        ## frequencies, and gives every band with no further circuit calls;
        ## a colour not of that shape gets the full evaluation.
        Nn = Pa.shape[0]
        if model is None:
            model = self._cy_components_model(pss, f, f0)
        ## ⚠ A SPECIFICATION LIMIT, NOT AN IMPLEMENTATION ONE: a coloured
        ## source under a modulation that CHANGES SIGN is not representable
        ## by any fold built from a PSD -- R(t,t') = m(t) m(t') R_c(t-t')
        ## keeps the sign product and CY cannot carry it -- so for a source
        ## that states no signed amplitudes this fold computes the |m|
        ## process.  A device's own 1/f current is a conductance fluctuation
        ## TIMES the current and follows its sign.  Where the element states
        ## its signed amplitudes (`Element.noise_amplitudes`) the folds use
        ## them and nothing below applies.  (White sources are untouched:
        ## uncorrelated across the period, no sign product survives.)
        ## The sign is invisible here; its NECESSARY condition is a PSD that
        ## touches zero along the orbit with a KINK in its square root, so
        ## that is warned on.  The touch threshold: a zero crossing SAMPLED
        ## on an N-point grid bottoms out near (pi/N)^2 of the maximum, a
        ## sign-definite PSD with a ten-fold swing sits at 1e-2 -- so 1e-2;
        ## a heuristic, and a warning for that reason.
        ## ⚠ NOT WHEN EVERY COLOURED COMPONENT CARRIES ITS SIGN: then nothing
        ## below takes a square root of a PSD and there is nothing to warn of
        _signed = getattr(model, 'amplitude', None) or {}
        self._warn_signed_unused(model, 'PAC.pnoise(cyclostationary=True)')
        _all_signed = (getattr(model, 'flicker', None) is not None
                       and all(k_ in _signed and (
                               self._uniform_exponent(B_, E_) is not None
                               or self._exponent_columns(B_, E_, _signed[k_]) is not None)
                               for k_, B_, E_ in model.flicker)
                       and all(self._perband_mode(pss, k_, [2.0 * np.pi * f0], None)
                               is not None for k_ in model.perband))
        Cs0 = np.asarray([np.abs(np.diag(np.fft.ifft(Pa, axis=0)[k])) for k in range(Nn)])
        dmax = Cs0.max(axis=0)
        touches = (dmax > 0) & (Cs0.min(axis=0) <= 1e-2 * dmax)
        ## ⚠ THE ORDER OF THE ZERO: a LINEAR sign crossing gives sqrt(PSD) a
        ## first-derivative KINK, a sign-definite quadratic touch a smooth
        ## minimum.  The circular second difference of sqrt(PSD) divided by
        ## h/T and by the maximum is a DERIVATIVE JUMP: grid-independent at a
        ## kink (~4 pi for a sinusoidal slope) and falling as h/T where
        ## smooth, so 3 separates them down to ~50 points per period (a raw
        ## threshold would encode the grid).  ⚠ STILL NECESSARY, NOT
        ## SUFFICIENT, AND THE SENSITIVITY RUNS INVERSE TO THE EFFECT: a
        ## crossing flatter than linear (an LO shaped v |v|^(p-1), p > 1)
        ## stays O(1) wrong while the indicator falls by orders.  A quiet
        ## warning is not evidence of a small discrepancy.
        kinked = np.zeros_like(touches)
        hT = 1.0 / float(Nn)
        for jj in np.where(touches)[0]:
            sq = np.sqrt(Cs0[:, jj])
            d2 = np.abs(sq - 0.5 * (np.roll(sq, 1) + np.roll(sq, -1)))
            kinked[jj] = bool(d2.max() / (sq.max() * hT) > 3.0)
        if bool(np.any(kinked)) and not _all_signed:
            warnings.warn(
                'PAC.pnoise(cyclostationary=True): a COLOURED source whose PSD '
                'touches zero along the orbit -- if its modulation changes sign '
                '(a switching gain), no PSD-specified model can represent the '
                'coloured process (Okumura eq. 23 in concrete form), and this '
                'fold computes the |m| one (its square root has a first-derivative '
                'kink at the zero, the signature of a LINEAR sign crossing; a '
                'necessary condition -- a shallow crossing shows no kink and '
                'errs MORE): measured 0.56x and 1.33x of the signed '
                'physics at two offsets on a flicker source through a '
                'zero-crossing gain -- EITHER direction, the sign of the '
                'discrepancy is set by the offset, not the mechanism -- and '
                'exact for a sign-definite one. Only the element knows the sign.',
                RuntimeWarning, stacklevel=3)
        ks = np.fft.fftfreq(Nn, d=1.0 / Nn).astype(int)
        pmin = min(ls) + int(ks.min()); pmax = max(ls) + int(ks.max())
        wband = lambda p: 2.0 * np.pi * abs(f - p * f0)
        ## ⚠⚠ ONE SQUARE ROOT PER INDEPENDENT COMPONENT, NOT OF THE SUM: a
        ## joint `sqrt(CY)` makes independent sources with different
        ## modulations NON-ADDITIVE -- see `_cy_components_model`.
        ## The white parts stay in the exact P-form (linear in `CY`, so
        ## additive already); each coloured part gets its own root, scaled
        ## per band when its exponent is uniform (`sqrt(c B) = sqrt(c)
        ## sqrt(B)`, no per-band eigendecomposition).
        ## the period's Fourier coefficients: the index DFT on a uniform grid,
        ## the trapezoid-weighted sum on a non-uniform one -- see `_period_dft`
        _dft = lambda B: self._period_dft(pss, B)
        if model is None:
            ## no component model (the whole `CY` is not thermal-plus-power-
            ## law and not the sum of its elements'): one root per band
            groups = [lambda p: self._cy_sqrt_harmonics(pss, wband(p))]
        else:
            Pw = _dft(model.white)
            for l in ls:
                for lp in ls:
                    total += complex(rows[l] @ Pw[(lp - l) % Nn] @ np.conj(rows[lp]))
            groups = []
            for _key, Bc, EF in model.flicker:
                ef = self._uniform_exponent(Bc, EF)
                if ef is not None:
                    ## the element's SIGNED amplitudes where it states them
                    ## (any factor with `W W^dagger = B` serves the pair sum;
                    ## only this one knows the sign), else the PSD's root
                    _W = _signed.get(_key)
                    groups.append(lambda p, SB=(_dft(_W) if _W is not None else
                                                self._sqrt_harmonics_of(Bc, _dft)), ef=ef:
                                  (model.w1 / wband(p)) ** (0.5 * ef) * SB)
                else:
                    ## signed columns grouped by their own exponents, where
                    ## they carry one each (`_exponent_columns`)
                    _split = (self._exponent_columns(Bc, EF, _signed[_key])
                              if _key in _signed else None)
                    if _split is not None:
                        for _Wg, _efg in _split:
                            _SBg = _dft(_Wg)
                            groups.append(lambda p, SB=_SBg, ef=_efg:
                                          (model.w1 / wband(p)) ** (0.5 * ef) * SB)
                        continue
                    groups.append(lambda p, Bc=Bc, EF=EF: self._sqrt_harmonics_of(
                        Bc * (model.w1 / wband(p)) ** EF, _dft))
            ## (the ONE element, `_one_element_cy`; its SIGNED amplitudes
            ## where it states them, `_perband_amplitudes`, else its root)
            ## History: `doc/shooting_history.md`, `PAC._cyclostationary_fold`.
            def _perband_at(p, key, mode):
                if mode is not None:
                    return _dft(self._perband_amplitudes(pss, key, wband(p),
                                                         None, mode))
                return self._sqrt_harmonics_of(
                    self._one_element_cy(pss, key, wband(p), None), _dft)
            for key in model.perband:
                ## (its columns where its amplitudes fit its CY at f0)
                mode = self._perband_mode(pss, key, [2.0 * np.pi * f0], None)
                groups.append(lambda p, key=key, mode=mode: _perband_at(p, key, mode))
        for sqrt_at in groups:
            cache = {}
            ## every band the sum reaches, stacked once: BB[pi, k] =
            ## B_k^{(p)} with pi = p - pmin.  The (l, l') pair sum is then two
            ## fancy indexings and one einsum instead of N small products in
            ## Python.
            def _B(p, cache=cache, sqrt_at=sqrt_at):
                key = round(abs(f - p * f0) / f0, 12)
                if key not in cache:
                    cache[key] = sqrt_at(p)
                return cache[key]
            ## (a band whose signed amplitudes did not rebuild its CY takes
            ## the root, with another column count: zero columns pad it,
            ## and add nothing to the pair sum)
            Bs = [_B(p) for p in range(pmin, pmax + 1)]
            r = max(b.shape[-1] for b in Bs)
            BB = np.asarray([b if b.shape[-1] == r else np.concatenate(
                (b, np.zeros(b.shape[:-1] + (r - b.shape[-1],), dtype=complex)),
                axis=-1) for b in Bs], dtype=complex)
            total = self._band_resolved_pairs(rows, ls, ks, Nn, pmin, BB, total)
        return float(np.real(total))

    @staticmethod
    def _band_resolved_pairs(rows, ls, ks, Nn, pmin, BB, total):
        for l in ls:
            for lp in ls:
                ## (B B^H)_j = sum_k B_k B_{k-j}^H: the partner index is
                ## k + l - l', NOT k + l' - l -- the mirror is invisible to
                ## the constant-modulation reduction (only k = 0 there).
                ## ⚠ NO CIRCULAR WRAP HERE: a partner beyond N/2 would be
                ## paired with the wrong BAND (each band carries its own
                ## weight), harmless in the white P-form and wrong here.
                ## History: `doc/shooting_history.md`,
                ## `PAC._band_resolved_pairs`.
                kp = ks + l - lp
                ok = np.abs(kp) <= Nn // 2
                kk, kk2 = ks[ok], kp[ok]
                pi = (l + kk) - pmin
                X = BB[pi, kk % Nn]
                Y = BB[pi, kk2 % Nn]
                Q = np.einsum('kij,klj->il', X, Y.conj())
                total += complex(rows[l] @ Q @ np.conj(rows[lp]))
        return float(np.real(total))

    def band_spread(self, pss, output, band, points=9, harmonic=1,
                    quantity='pnoise', **kw):
        """How much `S(r)·r²` VARIES across a band — the number that says
        whether a band mean and a point value are the same measurement.

        Returns `(spread, info)` with `spread = max/min` of `S(r)·r²` over
        `points` offsets spanning `band = (r_lo, r_hi)` in units of `f0`,
        and `info` carrying the samples, the band MEAN, the value at the
        band's midpoint, and their ratio.

        ⚠⚠ WHY THIS EXISTS.  Far above the AM corner both AM and PM fall as
        `1/r²`, so `S·r²` is flat and a band mean IS a point value, which
        makes the distinction invisible.  It is NOT general: a source behind
        a slow RC node has an in-band spectrum that is not `1/r²` at all
        (its `k = 0` term is filtered at the RC corner while the `k >= 1`
        terms are not, and their mix moves across the band), so a band mean
        and a point value can differ by ~4 % -- larger than most of the
        agreements this file asserts.  A comparison that takes a band mean
        on one side and a point value on the other is then measuring the
        convention, not the physics.

        ⚠ So: call this before comparing a measured band-averaged number
        against a computed point value, or vice versa.  A spread near 1
        licenses the shortcut; anything else says put both sides on the same
        footing.  `quantity` selects the surface (`'pnoise'`, `'S_pm'`,
        `'S_am'`, `'oscillator_spectrum'`); `**kw` is forwarded to it.

        History: `doc/shooting_history.md`, `PAC.band_spread`.
        """
        import numpy as _np
        f0 = 1.0 / float(pss.period)
        rs = _np.linspace(float(band[0]), float(band[1]), int(points))
        vals = []
        for r in rs:
            f = float(r) * f0
            if quantity == 'oscillator_spectrum':
                Sv, _i = self.oscillator_spectrum(pss, _np.array([f]), output,
                                                  harmonic=harmonic)
                v = float(_np.real(Sv[0]))
            elif quantity in ('S_pm', 'S_am'):
                am, pm, _b = self.am_pm_noise(pss, f, output, carrier=harmonic,
                                              **kw)
                v = float(_np.real(pm if quantity == 'S_pm' else am))
            else:
                v = float(_np.real(self.pnoise(pss, f, output, **kw)[0]))
            vals.append(v * float(r) ** 2)
        vals = _np.asarray(vals, dtype=float)
        lo = float(_np.min(_np.abs(vals)))
        spread = float(_np.max(_np.abs(vals)) / lo) if lo > 0.0 else _np.inf
        mean = float(_np.mean(vals))
        mid = float(_np.interp(0.5 * (rs[0] + rs[-1]), rs, vals))
        return spread, {'offsets': rs, 'values': vals, 'band_mean': mean,
                        'midpoint': mid,
                        'mean_over_point': (mean / mid) if mid != 0.0 else _np.inf}

    def am_pm_noise(self, pss, freq, output, carrier=1, maxsidebands=None,
                    modulated=False):
        """Output NOISE split into its AM and PM parts at `freq` from `carrier`.

        Returns `(S_am, S_pm, bands_used)`.  The two add to the noise in the
        pair of sidebands they decompose -- see the identity below -- and are in
        the same units as :meth:`pnoise`.

        ⚠ THIS NEEDS THE SIDEBAND *CORRELATION*, WHICH IS WHY IT IS NOT
        `|m_am|^2` FROM :meth:`am_pm`.  That method is the TRANSFER pair for a
        deterministic input; noise asks a different question, because whether
        the upper and lower sidebands are CORRELATED is exactly what decides the
        split.  Uncorrelated sidebands carry equal AM and PM -- the classical
        result for narrowband noise through an LTI system -- and it is the
        periodic operating point that correlates them.

        THE BAND BOOKKEEPING, which is the whole of the derivation and the one
        place a sign error would produce a plausible wrong answer.
        `adjoint_sideband_row(pss, g, output, l)` is the coefficient at output
        `g + l f0` for a unit source at `g`.  The two output sidebands sit at
        `carrier*f0 ± freq`, so a REAL noise band whose positive-frequency
        component is at `g = freq + p f0` reaches

            the UPPER output at `+g` through sideband `l = carrier - p`,
            the LOWER output at `-g` through sideband `l = carrier + p`,

        the second because a real process has `N(-g) = conj(N(g))` -- and that
        shared realisation IS the correlation.  Both contributions come from ONE
        band, so they are combined coherently; different `p` are different
        bands and are summed in power.  :meth:`am_pm` is exactly the `p = 0`
        term of this sum.

        The split per band is the same conjugate one :meth:`am_pm_indices`
        makes, `a + conj(b)` and `a - conj(b)` -- ⚠ the CONJUGATE, not `a ± b`:
        the sidebands counter-rotate about the carrier, and dropping it reports
        a rotating ellipse as pure AM.

        ⚠ THE GATE IS AN IDENTITY, NOT A TOLERANCE.  `pnoise` at the upper
        sideband folds precisely the bands `g = freq + p f0`, and at the lower
        precisely their negatives, so with the factor of one half below

            S_am + S_pm  ==  pnoise(carrier*f0 + freq) + pnoise(carrier*f0 - freq)

        exactly, because `|a+c|^2 + |a-c|^2 = 2|a|^2 + 2|c|^2` leaves no cross
        term.  A pairing error breaks it, which is what the test asserts.

        ⚠ ON A FREE-RUNNING OSCILLATOR this split sits on the SAME absolute
        scale as `pnoise` (the identity to 1e-12) and as the externally
        certified `oscillator_spectrum` (`S_pm = 4 S_v` at every offset: the
        PM content of the pair IS the Lorentzian, 2 S_v per sideband), with
        `S_am` rising from ~0 below the AM corner `f0/(2 pi Q_lambda)` to
        `S_pm` above it -- the ratio is an exact Lorentzian
        `u^2/(u_c^2 + u^2)` in `u = offset/f0`, `u_c = 1/(2 pi Q_lambda)`,
        `Q_lambda = -1/ln|lambda_2|` -- so the pair total is 4 S_v there and
        8 S_v far out.  "~1e-12 rows" from `am_pm` on a half-wave-symmetric
        fixture are a symmetry zero, see `am_pm`.  Oscillator magnitudes
        from this are trustworthy.

        ⚠ WHAT A MEASUREMENT MUST BE TO BE COMPARED WITH `S_pm`: `S_pm` is
        PM BY QUADRATURE OF THE FUNDAMENTAL'S SIDEBANDS.  A "phase" read by
        a one-period demodulation of the fundamental leaks the other
        harmonics' sidebands through its boxcar (sinc(pi(1 - r)) ~ 0.1 in
        amplitude for the second harmonic's), and a phase read from zero
        crossings converts EVERY harmonic's sidebands; both drift (6 %)
        against this quantity on a harmonic-rich asymmetric orbit.  So
        compare `S_pm` with the fundamental's sideband PM (a spectrum
        analyser's sidebands around f0, or the forward-tone gate
        `test_pnoise_oscillator_pm_matches_a_forward_tone_
        transient_with_no_adjoint`), never with a demodulated or
        crossing-time phase, and BAND WITH BAND (see `band_spread`).

        History: `doc/shooting_history.md`, `PAC.am_pm_noise`.
        """
        self._check_circuit(pss)
        pss = pss._adjoint_host()
        fp = pss.factored_period()
        N = len(fp.steps)
        f0 = 1.0 / float(fp.T)
        k = int(carrier)
        ## ⚠ BOTH SIDEBANDS OF A PAIR, `k - p` AND `k + p`, WITHIN THE GRID'S
        ## NYQUIST (`|l| <= N//2`, `adjoint_sideband_row`): so `|p|` up to
        ## `N//2 - |k|`.  The default used to be `N//2` itself, and every
        ## call with `carrier >= 1` and no `maxsidebands` raised (2026-09-28).
        cap = N // 2 - abs(k)
        if cap < 0:
            raise ValueError(
                'PAC.am_pm_noise: carrier %d is above the grid\'s Nyquist '
                '(%d harmonics at %d points per period) -- use a finer period '
                'grid.' % (k, N // 2, N))
        lmax = cap if maxsidebands is None else min(int(maxsidebands), cap)
        cyfn = (self._cy_cycle_averaged if modulated else self._cy_reduced)
        ## ⚠⚠ THE SPLIT IS TAKEN IN THE CARRIER'S FRAME, NOT THE TIME
        ## ORIGIN'S.  AM is the envelope component ALONG the carrier phasor,
        ## so `a + conj(b)` is right only for a cosine-phased carrier;
        ## unrotated, the split depends on where t = 0 sits, and the
        ## identity cannot see it (`|a_r|`, `|b_r|` are `|a|`, `|b|`).
        ## `am_pm` divides by the COMPLEX carrier phasor instead.  With no
        ## carrier at this harmonic the phase is undefined and the split is
        ## left unrotated, as `am_pm` refuses the same case.
        _C = self.carrier_phasor(pss, output, k)
        _scale = float(np.max(np.abs(self._output_waveform_row(pss, output))))
        _rot = (np.exp(-1j * np.angle(_C))
                if abs(_C) > 1e-9 * max(_scale, 1e-300) else 1.0)
        S_am = 0.0
        S_pm = 0.0
        bands = []
        for p in range(-lmax, lmax + 1):
            g = float(freq) + p * f0
            a = self.adjoint_sideband_row(pss, g, output, k - p)[0]
            b = self.adjoint_sideband_row(pss, -g, output, k + p)[0]
            cy = cyfn(pss, 2.0 * np.pi * g)
            a_r = a * _rot
            b_r = np.conj(b) * np.conj(_rot)
            m_am = a_r + b_r
            m_pm = a_r - b_r
            S_am += 0.5 * float(np.real(m_am @ cy @ np.conj(m_am)))
            S_pm += 0.5 * float(np.real(m_pm @ cy @ np.conj(m_pm)))
            bands.append(p)
        return S_am, S_pm, bands
