"""Driven-circuit noise: pnoise (stationary, modulated and cyclostationary
folds), its AM/PM split and the band spread.
"""
import numpy as np
from ._noise_components import warn_sign_blind
from ._numerics import output_index, sweep_frequency, sweep_offset
from pycircuit.circuit.simwarnings import AccuracyWarning, warn


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
               modulated=False, cyclostationary=False, sweeptype=None,
               relharmnum=None):
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

        ⚠ THE SWEEP (`sweeptype`, `relharmnum`), a commercial RF simulator's
        rule: `'absolute'` reads `freq` as the output frequency itself,
        `'relative'` as the offset `relharmnum * f0 + freq` (`relharmnum`
        default 1), and `None` is relative on an AUTONOMOUS PSS and absolute
        on a driven one (`_numerics.sweep_kind`).

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
        fired: 'ratio', 'bound' (the Nyquist) or 'cap' (an explicit
        `maxsidebands` below it); ending on either of the last two warns.

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
        bordered with BOTH null vectors) on an oscillator, the plain one on
        a driven circuit (`_sideband_family`).

        History: `doc/shooting_history.md`, `PAC.pnoise`.
        """
        output = output_index(pss, output)
        self._check_circuit(pss)
        if modulated and cyclostationary:
            raise ValueError(
                'PAC.pnoise: modulated=True and cyclostationary=True are two '
                'models of the same bias-dependent source -- one stationary '
                'source at the cycle-averaged bias, or the cyclostationary '
                'construction; choose one.')
        ## pnoise folds sidebands through the ADJOINT (`_sideband_family`: one
        ## reverse pass per output frequency, `_reverse_points`) on the
        ## monodromy twin -- see `_adjoint_host`: the stage methods' fold is
        ## built, so there is no Gear-2 fallback; a Gear-2 twin only if
        ## `monodromy='gear'` asks for one.
        pss = pss._adjoint_host()
        ## ⚠ THE SWEEP ON THE HOST: a relative frequency is an offset from
        ## the carrier of the operator the rows are solved on -- the twin's
        ## on a trap/euler oscillator, whose period differs from the run's by
        ## O(h^2): read off the run's period, an offset below that gap would
        ## land on the wrong side of the pole.
        freq = sweep_frequency(pss, freq, sweeptype, relharmnum, 'pnoise')
        fp = pss.factored_period()
        m = pss.cir.n - 1
        N = len(fp.steps)
        T = float(fp.T)
        f0 = 1.0 / T
        tol = self.ALIAS_RATIO_TOL if ratio_tol is None else float(ratio_tol)
        ## ⚠ A COUNT ABOVE THE GRID'S NYQUIST RAISES (as the sampled family
        ## does)
        if maxsidebands is not None and int(maxsidebands) > N // 2:
            raise ValueError(
                f'PAC.pnoise: maxsidebands={maxsidebands} is above the grid\'s '
                f'Nyquist ({N // 2} sidebands at {N} points per period) -- '
                'use a finer period grid, or leave it None (the default, '
                'every sideband the grid resolves).')
        lmax = N // 2 if maxsidebands is None else int(maxsidebands)

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
            ## `NoiseComponents.model`); its call is the summed model the
            ## stop rule reads
            colour = self._noise_components(pss).model(float(freq), f0)
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
        ## ⚠ EVERY SIDEBAND THE FOLD USES (`lmax`): a scan shorter than the
        ## fold lets a harmonic above it read a 1/f source next to DC
        ## unrefused.
        lscan = max(1, lmax)
        offs = np.abs(float(freq) - np.arange(-lscan, lscan + 1) * f0_)
        self._dc_fold_guard(pss, cyfn, float(freq), float(np.min(offs)), f0_,
                            cy, 'pnoise')

        total = 0.0
        used = []
        quiet = 0
        self.alias_stop = 'bound'
        rows = {}
        ## every row here has OUTPUT frequency `freq`: one family, one
        ## adjoint solve for all the sidebands (`_sideband_family`)
        fam = self._sideband_family(pss, float(freq), output)
        for l in range(0, lmax + 1):
            step = 0.0
            for sl in ((0,) if l == 0 else (l, -l)):
                fin = float(freq) - sl * f0
                h = fam.row(sl, fin)
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
        ## ⚠ A CAP THE CALLER SET IS NAMED AS THAT CAP, not as the grid's
        ## Nyquist.
        if (self.alias_stop == 'bound' and maxsidebands is not None
                and lmax < N // 2):
            self.alias_stop = 'cap'
        ## (an explicit `maxsidebands=0` asks for the unfolded term alone,
        ## and stays quiet, as it always did)
        if self.alias_stop == 'cap' and lmax > 0:
            warn(
                'PAC.pnoise: the sideband accumulation stopped at '
                f'maxsidebands={lmax} (the grid resolves {N // 2} at {N} '
                'points per period), not because the contributions became '
                'negligible. Sidebands above the cap are MISSING rather '
                'than small, so this is a lower bound on the folded noise. '
                'Raise maxsidebands (or leave it None) and compare.', AccuracyWarning)
        elif self.alias_stop == 'bound' and lmax > 0:
            warn(
                'PAC.pnoise: the sideband accumulation stopped at the '
                "grid's Nyquist (|l| = %d at %d points per period), not "
                'because the contributions became negligible. Sidebands '
                'above the grid\'s maximum frequency are MISSING rather '
                'than small, so this is a lower bound on the folded noise. '
                'Re-solve the PSS on a finer period grid and compare.'
                % (lmax, N), AccuracyWarning)
        return total, used

    def _dc_fold_guard(self, pss, cyfn, freq, near, f0_, cy, what):
        """Refuse (or warn about) a fold that reads the sources next to DC.

        `freq` the output frequency, `near` the smallest |source frequency|
        the fold evaluates `cyfn` at, `cy` the sources at `freq` (None: read
        here when needed).  Used by `pnoise` and `am_pm_noise`.  See `pnoise`
        for the cases.

        History: `doc/shooting_history.md`, `PAC._dc_fold_guard`."""
        if near <= self.HARMONIC_GUARD * f0_:
            probe = cyfn(pss, 2.0 * np.pi * near)
            if not np.all(np.isfinite(np.asarray(probe))):
                raise ValueError(
                    'PAC.%s: %.12g Hz sits on a harmonic of %.12g Hz, '
                    'so a sideband folds the noise sources to DC -- and at '
                    'DC this circuit\'s CY is not finite. A 1/f term is '
                    'infinite there; a flicker term whose COEFFICIENT IS '
                    'ZERO is 0/0 and gives nan, so disabling flicker does '
                    'not avoid this. Offset from the harmonic: a commercial RF simulator\'s '
                    'own advice is to cluster frequencies NEAR each '
                    'harmonic and never place one ON it.'
                    % (what, float(freq), f0_))
            ## ⚠⚠ A FINITE PROBE IS NOT A SAFE ONE.  `1/T` rounds, so at
            ## `f = f0` the folded band sits ~1e-11 Hz from DC, not ON it: a
            ## 1/f source there is finite and enormous (6.3e-2 V^2/Hz against
            ## 9.2e-15 at 0.1 % either side).  So a frequency-DEPENDENT CY at
            ## the folded band refuses too; a white one stays allowed.
            probe_ref = cyfn(pss, 2.0 * np.pi * max(abs(float(freq)) * 2.0, f0_))
            if not np.allclose(np.asarray(probe), np.asarray(probe_ref),
                               rtol=1e-9, atol=0.0):
                raise ValueError(
                    'PAC.%s: %.12g Hz sits on harmonic %d of %.12g Hz, so '
                    'a sideband folds the noise sources to %.3g Hz -- DC up to '
                    'rounding -- and this circuit\'s CY is frequency-dependent '
                    'there: a 1/f source is read next to its singularity and '
                    'the fold returns a finite, absurd number (measured '
                    '6.3e-2 V^2/Hz against 9.2e-15 at 0.1 %% either side). '
                    'Offset from the harmonic, or use PAC.sampled_variance, '
                    'whose explicit fmin keeps every band off DC.'
                    % (what, float(freq), int(round(abs(float(freq)) / f0_)),
                       f0_, near))

        ## ⚠ THE STEEP REGION BESIDE A HARMONIC IS A SWEEP HAZARD RATHER
        ## THAN A WRONG NUMBER, so it warns instead of raising: the VALUE is
        ## right (2 % above the plateau at `f0 + 0.01` Hz with a real flicker
        ## source), but a grid that lands there by accident integrates a
        ## spike it never resolved.
        elif near < 1e-6 * f0_ and float(freq) > 0.0:
            if cy is None:
                cy = cyfn(pss, 2.0 * np.pi * float(freq))
            cy_hi = cyfn(pss, 2.0 * np.pi * max(float(freq) * 2.0, f0_))
            if not np.allclose(cy, cy_hi, rtol=1e-9, atol=0.0):
                warn(
                    'PAC.%s: %.12g Hz is %.3g Hz from a harmonic of '
                    '%.12g Hz and a source has a frequency-dependent CY, '
                    'so the folded density varies steeply here. The VALUE '
                    'is correct; a swept grid landing this close will '
                    'misrepresent the integrated total. Cluster near each '
                    'harmonic deliberately rather than by accident.'
                    % (what, float(freq), near, f0_), AccuracyWarning)

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
        return self._period_dft(pss, self._noise_components(pss).cy_at_states(w))

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
        return self._sqrt_harmonics_of(self._noise_components(pss).cy_at_states(w),
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

        The colour model (`NoiseComponents.model`, fitted once in `pnoise`
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
        nc = self._noise_components(pss)
        if model is None:
            model = nc.model(f, f0)
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
        ## The sign is invisible here: a rooted component whose PSD touches
        ## zero is warned by the one verdict every surface gives
        ## (`NoiseComponents.sign_blind`, below; until 2026-10-01 this fold
        ## had its own test -- a first-derivative KINK of the root of the
        ## circuit's TOTAL `CY` -- which missed a crossing flatter than
        ## linear (2.98x the signed physics, silent) and any touch a white
        ## source on the same node filled in (1.67x, silent)).
        ks = np.fft.fftfreq(Nn, d=1.0 / Nn).astype(int)
        pmin = min(ls) + int(ks.min()); pmax = max(ls) + int(ks.max())
        wband = lambda p: 2.0 * np.pi * abs(f - p * f0)
        ## ⚠⚠ ONE SQUARE ROOT PER INDEPENDENT COMPONENT, NOT OF THE SUM: a
        ## joint `sqrt(CY)` makes independent sources with different
        ## modulations NON-ADDITIVE -- see `NoiseComponents.model`.
        ## The white parts stay in the exact P-form (linear in `CY`, so
        ## additive already); each coloured part gets its own root, scaled
        ## per band when its exponent is uniform (`sqrt(c B) = sqrt(c)
        ## sqrt(B)`, no per-band eigendecomposition).
        ## the period's Fourier coefficients: the index DFT on a uniform grid,
        ## the trapezoid-weighted sum on a non-uniform one -- see `_period_dft`
        _dft = lambda B: self._period_dft(pss, B)
        rooted = []
        if model is None:
            ## no component model (the whole `CY` is not thermal-plus-power-
            ## law and not the sum of its elements'): one root per band
            groups = [lambda p: self._cy_sqrt_harmonics(pss, wband(p))]
            rooted.append(nc.JOINT_KEY)
        else:
            Pw = _dft(model.white)
            for l in ls:
                for lp in ls:
                    total += complex(rows[l] @ Pw[(lp - l) % Nn] @ np.conj(rows[lp]))
            ## the components as EVERY surface groups them
            ## (`colour_components`, review O2, 2026-10-01), each band read
            ## EXACTLY (no classification shortcut): a uniform power law's
            ## fixed columns scaled per band (`sqrt(c B) = sqrt(c) sqrt(B)`,
            ## no per-band eigendecomposition) -- the element's SIGNED
            ## amplitudes where it states them (any factor with `W W^dagger
            ## = B` serves the pair sum; only this one knows the sign), else
            ## the PSD's root -- and the rest rooted (or its signed columns)
            ## per band.  The sign-blind verdict is given there.
            ## History: `doc/shooting_history.md`, `PAC._cyclostationary_fold`.
            _ws = [wband(p) for p in range(pmin, pmax + 1)]
            comps = nc.colour_components(
                model, max(min(w_ for w_ in _ws if w_ > 0.0), 2e-3 * np.pi * f0),
                max(_ws), f0, 'pnoise(cyclostationary=True)', shortcuts=False)
            ## a power law's ONE root and its per-band weight `(w1/w)^EF`
            ## (`_band_resolved_pairs`' `scale2`): the root is not copied per
            ## band (it was, stacked, until 2026-10-01 -- the review's M4)
            ws_all = np.array([wband(p) for p in range(pmin, pmax + 1)])
            for _key, _W, ef, _Bc in comps.fixed:
                SB = (_dft(_W) if _Bc is None
                      else self._sqrt_harmonics_of(_Bc, _dft))
                with np.errstate(divide='ignore'):
                    scale2 = (model.w1 / ws_all) ** float(ef)
                total = self._band_resolved_pairs(rows, ls, ks, Nn, pmin, SB,
                                                  total, scale2=scale2)
            groups = []
            for band in comps.bands:
                if band.signed:
                    groups.append(lambda p, root=band.exact_root:
                                  _dft(root(wband(p))))
                else:
                    groups.append(lambda p, psd=band.exact_psd:
                                  self._sqrt_harmonics_of(psd(wband(p)), _dft))
        if model is None:
            blind = nc.sign_blind(None, rooted)
            if blind:
                warn_sign_blind('pnoise(cyclostationary=True)', blind)
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
    def _band_resolved_pairs(rows, ls, ks, Nn, pmin, BB, total, scale2=None):
        """The pair sum over one component: `BB[pi, k]` its root per band
        (`pi = p - pmin`), or, with `scale2`, `BB[k]` one root for every band
        and `scale2[pi]` the band's weight on `B B^H` (a power law's
        `(w1/w)^EF`)."""
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
                if scale2 is None:
                    X = BB[pi, kk % Nn]
                    Y = BB[pi, kk2 % Nn]
                    Q = np.einsum('kij,klj->il', X, Y.conj())
                else:
                    Q = np.einsum('k,kij,klj->il', scale2[pi], BB[kk % Nn],
                                  BB[kk2 % Nn].conj())
                total += complex(rows[l] @ Q @ np.conj(rows[lp]))
        return float(np.real(total))

    def band_spread(self, pss, output, band, points=9, harmonic=1,
                    quantity='pnoise', **kw):
        """How much `S(r)·r²` VARIES across a band — the number that says
        whether a band mean and a point value are the same measurement.

        Returns `(spread, info)` with `spread = max/min` of `S(r)·r²` over
        `points` offsets spanning `band = (r_lo, r_hi)` in units of `f0`,
        and `info` carrying the samples, the band MEAN, the value at the
        band's midpoint, and their ratio.  The offsets are from harmonic
        `harmonic` for EVERY quantity; the sweep is fixed here, so `**kw`
        takes no `sweeptype` / `relharmnum`.

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
        output = output_index(pss, output)
        if quantity not in ('pnoise', 'S_pm', 'S_am', 'oscillator_spectrum'):
            raise ValueError(
                "PAC.band_spread: quantity must be 'pnoise', 'S_pm', 'S_am' "
                f"or 'oscillator_spectrum', not {quantity!r}")
        for k in ('sweeptype', 'relharmnum'):
            if k in kw:
                raise ValueError(
                    f'PAC.band_spread: {k} is not a knob here -- the band is '
                    'an offset from `harmonic` for every quantity.')
        f0 = 1.0 / float(pss.period)
        rs = np.linspace(float(band[0]), float(band[1]), int(points))
        vals = []
        for r in rs:
            f = float(r) * f0
            if quantity == 'oscillator_spectrum':
                Sv, _i = self.oscillator_spectrum(pss, np.array([f]), output,
                                                  harmonic=harmonic, **kw)
                v = float(np.real(Sv[0]))
            elif quantity in ('S_pm', 'S_am'):
                am, pm, _b = self.am_pm_noise(pss, f, output, harmonic=harmonic,
                                              sweeptype='relative', **kw)
                v = float(np.real(pm if quantity == 'S_pm' else am))
            else:
                v = float(np.real(self.pnoise(
                    pss, f, output, sweeptype='relative',
                    relharmnum=int(harmonic), **kw)[0]))
            vals.append(v * float(r) ** 2)
        vals = np.asarray(vals, dtype=float)
        lo = float(np.min(np.abs(vals)))
        spread = float(np.max(np.abs(vals)) / lo) if lo > 0.0 else np.inf
        mean = float(np.mean(vals))
        mid = float(np.interp(0.5 * (rs[0] + rs[-1]), rs, vals))
        return spread, {'offsets': rs, 'values': vals, 'band_mean': mean,
                        'midpoint': mid,
                        'mean_over_point': (mean / mid) if mid != 0.0 else np.inf}

    def am_pm_noise(self, pss, freq, output, harmonic=1, maxsidebands=None,
                    modulated=False, sweeptype=None):
        """Output NOISE split into its AM and PM parts at `freq` from `harmonic`.

        Returns `(S_am, S_pm, bands_used)`: the AM and PM noise PER SIDEBAND,
        one-sided densities in the units of :meth:`pnoise` -- half the pair
        of sidebands they decompose, see the identity below.  `output`: a
        reduced index, a weight vector, or a node name.
        History: `doc/shooting_history.md`, `PAC.am_pm_noise`.

        ⚠ THE SWEEP (`sweeptype`), a commercial RF simulator's rule with
        `harmonic` as the reference harmonic: `'relative'` reads `freq` as
        the offset from `harmonic*f0`, `'absolute'` as the UPPER output
        sideband's own frequency (the offset is `freq - harmonic*f0`), and
        `None` is relative on an AUTONOMOUS PSS and absolute on a driven one
        (`_numerics.sweep_kind`).  Below, `freq` is the offset.

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
        `harmonic*f0 ± freq`, so a REAL noise band whose positive-frequency
        component is at `g = freq + p f0` reaches

            the UPPER output at `+g` through sideband `l = harmonic - p`,
            the LOWER output at `-g` through sideband `l = harmonic + p`,

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
        precisely their negatives, so

            2 (S_am + S_pm)  ==  pnoise(harmonic*f0 + freq) + pnoise(harmonic*f0 - freq)

        exactly, because `|a+c|^2 + |a-c|^2 = 2|a|^2 + 2|c|^2` leaves no cross
        term.  A pairing error breaks it, which is what the test asserts.

        ⚠ ON A FREE-RUNNING OSCILLATOR this split sits on the SAME absolute
        scale as `pnoise` (the identity to 1e-12) and as the externally
        certified `oscillator_spectrum` (`S_pm = S_v` at every offset: the
        PM content of a sideband IS the Lorentzian), with
        `S_am` rising from ~0 below the AM corner `f0/(2 pi Q_lambda)` to
        `S_pm` above it -- the ratio is an exact Lorentzian
        `u^2/(u_c^2 + u^2)` in `u = offset/f0`, `u_c = 1/(2 pi Q_lambda)`,
        `Q_lambda = -1/ln|lambda_2|` -- so a sideband's total is S_v there and
        2 S_v far out.  "~1e-12 rows" from `am_pm` on a half-wave-symmetric
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
        output = output_index(pss, output)
        self._check_circuit(pss)
        pss = pss._adjoint_host()
        ## (the offset from the HOST's carrier, as in `pnoise`)
        freq = sweep_offset(pss, freq, sweeptype, harmonic, 'am_pm_noise')
        fp = pss.factored_period()
        N = len(fp.steps)
        f0 = 1.0 / float(fp.T)
        k = int(harmonic)
        ## ⚠ BOTH SIDEBANDS OF A PAIR, `k - p` AND `k + p`, WITHIN THE GRID'S
        ## NYQUIST (`|l| <= N//2`, `adjoint_sideband_row`): so `|p|` up to
        ## `N//2 - |k|`.
        cap = N // 2 - abs(k)
        if cap < 0:
            raise ValueError(
                'PAC.am_pm_noise: harmonic %d is above the grid\'s Nyquist '
                '(%d harmonics at %d points per period) -- use a finer period '
                'grid.' % (k, N // 2, N))
        ## (an explicit count above that raises, as in `pnoise`)
        if maxsidebands is not None and int(maxsidebands) > cap:
            raise ValueError(
                f'PAC.am_pm_noise: maxsidebands={maxsidebands} is above what '
                f'the grid resolves at harmonic {k} ({cap}: both sidebands of '
                f'a pair within the Nyquist, {N // 2} at {N} points per '
                'period) -- use a finer period grid, or leave it None.')
        lmax = cap if maxsidebands is None else int(maxsidebands)
        cyfn = (self._cy_cycle_averaged if modulated else self._cy_reduced)
        ## the sources are read at `freq + p f0`: next to DC when the output
        ## sits on a harmonic (`pnoise`'s guard)
        gs = float(freq) + np.arange(-lmax, lmax + 1) * f0
        self._dc_fold_guard(pss, cyfn, k * f0 + float(freq),
                            float(np.min(np.abs(gs))), f0, None, 'am_pm_noise')
        ## ⚠⚠ THE SPLIT IS TAKEN IN THE CARRIER'S FRAME, NOT THE TIME
        ## ORIGIN'S.  AM is the envelope component ALONG the carrier phasor,
        ## so `a + conj(b)` is right only for a cosine-phased carrier;
        ## unrotated, the split depends on where t = 0 sits, and the
        ## identity cannot see it (`|a_r|`, `|b_r|` are `|a|`, `|b|`).
        ## `am_pm` divides by the COMPLEX carrier phasor instead.  With no
        ## carrier at this harmonic the phase is undefined, and the split is
        ## REFUSED, as `am_pm` refuses it: unrotated, the split would move
        ## with where t = 0 sits.
        _C = self.carrier_phasor(pss, output, k)
        _scale = float(np.max(np.abs(self._output_waveform_row(pss, output))))
        if abs(_C) <= 1e-9 * max(_scale, 1e-300):
            raise ValueError(
                'PAC.am_pm_noise: the output carries no component at '
                f'harmonic {k} (|C| = {abs(_C):.3e} against a signal scale '
                f'of {_scale:.3e}), so there is no carrier to split the noise '
                'against; `pnoise` at `harmonic*f0 +- freq` is the total.')
        _rot = np.exp(-1j * np.angle(_C))
        S_am = 0.0
        S_pm = 0.0
        bands = []
        ## the `a` rows share OUTPUT frequency `k f0 + freq`, the `b` rows
        ## `k f0 - freq`: two families, two adjoint solves in all
        fam_a = self._sideband_family(pss, k * f0 + float(freq), output)
        fam_b = self._sideband_family(pss, k * f0 - float(freq), output)
        for p in range(-lmax, lmax + 1):
            g = float(freq) + p * f0
            a = fam_a.row(k - p, g)
            b = fam_b.row(k + p, -g)
            cy = cyfn(pss, 2.0 * np.pi * g)
            a_r = a * _rot
            b_r = np.conj(b) * np.conj(_rot)
            m_am = a_r + b_r
            m_pm = a_r - b_r
            S_am += 0.5 * float(np.real(m_am @ cy @ np.conj(m_am)))
            S_pm += 0.5 * float(np.real(m_pm @ cy @ np.conj(m_pm)))
            bands.append(p)
        ## the pair's totals, halved (exactly): per sideband
        return 0.5 * S_am, 0.5 * S_pm, bands
