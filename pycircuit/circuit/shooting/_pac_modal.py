"""The Floquet-mode spectra: the orbital correlation and spectrum, the modal
and correlation spectra, and the phase-mode split.
"""
import numpy as np
import warnings
from ._noise_components import orbit_states
from ._numerics import _output_row, insert_ref, integer_arg, output_index
from pycircuit.circuit._limiting import devices_at
from pycircuit.circuit.simwarnings import (
    AccuracyWarning,
    CostWarning,
    ModelWarning,
    warn,
)


class _ModalSpectra(object):
    """The Floquet-mode spectra: the orbital correlation and spectrum, the
    modal and correlation spectra, and the phase-mode split.  A theme of
    `PAC` (see `pac.py`)."""

    ORBITAL_HARMONICS = 32

    ## Half-wave asymmetry above which `orbital_correlation` warns of its O(h)
    ## residual and `orbital_spectrum` of its over-statement; 0.02 is a decade
    ## inside the smallest asymmetry at which the error was visible.
    ORBITAL_ASYMMETRY_LIMIT = 0.02

    def _harmonic_count(self, maxharmonics, N, what):
        """The PPV / Floquet-mode harmonics kept, `|h| <= H`: `maxharmonics`,
        or by default `ORBITAL_HARMONICS` capped by what the `N`-point grid
        resolves (`N//2 - 1`).  ⚠ An EXPLICIT count above that raises; it
        was capped without a word until 2026-09-29 (as `H`)."""
        cap = N // 2 - 1
        if maxharmonics is None:
            return min(self.ORBITAL_HARMONICS, cap)
        maxharmonics = integer_arg(maxharmonics, 'maxharmonics', what,
                                   minimum=0)
        if maxharmonics > cap:
            raise ValueError(
                f'PAC.{what}: maxharmonics={maxharmonics} is above what the '
                f'grid resolves ({cap} at {N} points per period) -- use a '
                'finer period grid, or leave it None (the default, '
                f'{self.ORBITAL_HARMONICS} capped by the grid).')
        return maxharmonics

    def _sideband_count(self, maxsidebands, H, N, what):
        """The input sidebands the modal sum folds, `|m| <= M`: `2 H` by
        default, an explicit `maxsidebands` at most `2 (N//2 - 1)` -- the
        most the default can reach -- above which it raises (it was taken
        as given until 2026-10-01)."""
        if maxsidebands is None:
            return 2 * H
        M = integer_arg(maxsidebands, 'maxsidebands', what, minimum=0)
        cap = 2 * (N // 2 - 1)
        if M > cap:
            raise ValueError(
                f'PAC.{what}: maxsidebands={M} is above what the grid '
                f'resolves ({cap} at {N} points per period: twice the '
                'harmonics it holds) -- use a finer period grid, or leave it '
                'None (twice maxharmonics).')
        return M

    def _orbit_asymmetry(self, pss):
        """Half-wave asymmetry of the orbit, in [0, ~1].

        `max|x(t) + x(t + T/2)| / max|x|` on the first state row -- zero for a
        half-wave symmetric orbit (van der Pol), growing as the orbit
        distorts.  Cheap: the waveform is already stored.
        """
        W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                      pss.irefnode, axis=0)
        if W.size == 0 or W.shape[1] < 4:
            return 0.0
        row = W[0]
        half = len(row) // 2
        den = float(np.max(np.abs(row)))
        if den <= 0.0:
            return 0.0
        return float(np.max(np.abs(row[:half] + row[half:2 * half]))) / den

    def orbital_correlation(self, pss, maxharmonics=None):
        """`R_yy(0)` and the `C_lhj` of Traversa & Bonani eq (22) — A9 step 3.

        Returns `(R, C)`.  `R` is the STATIONARY transverse (orbital)
        state covariance, `m x m` real symmetric — eq (23),
        `R = Σ_{l≥2,h,j} C_lhj`.  `C` maps `(l, h, j)` to the `m x m`
        complex coefficient, over every non-null orbital mode `l ≥ 2` and
        harmonics `|h|, |j|, |j'| ≤ H`, `H = maxharmonics`, default
        `ORBITAL_HARMONICS` (capped by the grid; an explicit count above
        it raises); van der Pol converges by `H = 4`, and a
        strongly non-sinusoidal orbit needs more — check by raising it.

        ⚠ STATIONARY WHITE SOURCES ONLY, and that is what makes it
        computable without `B`.  Eq (22) needs the Fourier coefficients
        of `v_l(t)^T B(t)`; with `CY = B B^T` constant those products
        collapse to `V~_{l'k}^T CY V~*_{lk'}`, so the noise enters only
        through the reduced `CY` that `_cy_reduced` already refuses to
        hand over when it is bias-dependent.

        ⚠ `CY/2`, NOT `CY`.  The library's `CY` is one-sided; eq (22)
        integrates `B B^T` as a two-sided intensity -- the `kT/C`-calibrated
        Monte Carlo injection `Var(i) = CY/(2h)`, confirmed here by three
        routes agreeing.

        ⚠⚠ GATED THREE WAYS, because a modal sum transcribed from an image
        of an equation is exactly the object to distrust.  (i) This sum
        against `R_yy(0)` evaluated from its DEFINITION as a 1-D Lyapunov
        integral along the orbital mode, no Fourier machinery; (ii) both
        against the CYCLE-MEAN transverse part of `oscillator_covariance`'s
        samples, which shares no machinery with either.

        ⚠ THE REFERENCE IS THE CYCLE MEAN, NOT `K_orb(0)`.  Lemma 3.5's
        `R∞_yy` depends on `τ` only — the stationary part.  At `t = 0`
        van der Pol's amplitude direction is pure-v while this is
        isotropic, which is a rotating radial direction averaged over a
        cycle, not a disagreement.

        ⚠ AND THE REFERENCE IS OBLIQUELY PROJECTED.  Subtracting only the
        SECULAR growth `(t/T) d u u^T` from the Lyapunov samples leaves the
        phase direction's BOUNDED within-period variance, which eq (22)'s
        `l >= 2` sum correctly excludes.  Demir's `y` is defined by the
        OBLIQUE projection `v_1^T y = 0`: project the samples with
        `Pi = I - u v^T/(v^T u)`.  Unprojected, the reference is off by a
        residual falling as 1/Q.
        (The phase-orbital CORRELATION's tau = 0 value is not a missing
        term: eq (18a)'s brace is {1 - 1} = 0 there, and eq (23) states
        R_yy(0) = sum C_lhj alone.)

        History: `doc/shooting_history.md`, `PAC.orbital_correlation`.
        """
        return self._orbital_correlation(pss, maxharmonics,
                                         'orbital_correlation')

    def _orbital_correlation(self, pss, maxharmonics, what):
        """`orbital_correlation`, its refusals naming `PAC.<what>`, the
        method the user called."""
        self._check_circuit(pss)
        self._refuse_driven(pss, what)
        self._refuse_coloured(
            pss, what,
            'eq (22) is a white-noise residue sum; modal_spectrum takes a '
            'stationary coloured source per sideband, and '
            'oscillator_covariance(colour_fmin=...) its transverse covariance.')
        modes = pss.floquet_modes()
        m = self.cir.n - 1
        Tp = float(pss.period)
        w0 = 2.0 * np.pi / Tp
        CY2 = 0.5 * np.real(np.asarray(self._cy_reduced(pss, 0.0)))
        ## the phase mode by its tangent alignment (see `_phase_mode_split`),
        ## and NEVER in the orbital sum: swept in with a near-zero exponent it
        ## blows up as 1/|mu|^2
        _kph, orb = self._phase_mode_split(pss, modes, 'PAC.' + what)
        if not orb:
            raise ValueError(
                f'PAC.{what}: no orbital mode -- every non-null '
                'multiplier is the phase mode.')

        def fcoef(P):
            ## `_period_dft`: the index DFT on a uniform grid, unchanged, and
            ## the trapezoid-weighted sum at the true times on a non-uniform one
            X = np.asarray(P)[:, :-1]
            N = X.shape[1]
            return self._period_dft(pss, X.T).T, N

        U, V, N = {}, {}, None
        for k in orb:
            U[k], N = fcoef(modes[k]['p'])
            V[k], _ = fcoef(modes[k]['q'])
        H = self._harmonic_count(maxharmonics, N, what)
        hs = np.arange(-H, H + 1)
        idx = lambda k: k % N

        C = {}
        R = np.zeros((m, m), dtype=complex)
        for l in orb:
            mul = modes[l]['mu']
            for lp in orb:
                mulp = modes[lp]['mu']
                for j in hs:
                    Ulj = U[l][:, idx(j)]
                    ## the Lambda products for every (h, j') at once
                    Vl_hj = V[l][:, idx(hs - j)]              # m x nh  (h - j)
                    for jp in hs:
                        res = 1.0 / (1j * (j - jp) * w0 - mulp - np.conj(mul))
                        outer = res * np.outer(U[lp][:, idx(jp)], np.conj(Ulj))
                        Vlp_hjp = V[lp][:, idx(hs - jp)]      # m x nh  (h - j')
                        sc = np.einsum('ih,ik,kh->h', Vlp_hjp, CY2, np.conj(Vl_hj))
                        for hi, h in enumerate(hs):
                            term = sc[hi] * outer
                            key = (l, int(h), int(j))
                            C[key] = C.get(key, 0.0) + term
                            R += term
        return np.real(R), C

    def orbital_spectrum(self, pss, offsets, output, harmonic=1,
                         maxharmonics=None):
        """`S_yy` — the ORBITAL (amplitude) noise spectrum. A9 step 4.

        Returns `S` at `harmonic*f0 + offsets` (a negative offset is the
        LOWER sideband, which is not the upper's), a ONE-SIDED PSD in V^2/Hz,
        the scale of `oscillator_spectrum`'s `S_v` and of `pnoise` (since
        2026-09-28; it was half that), so **the two are summed** — which is
        what Traversa & Bonani (TCAS-I 2011) say to do:

            x(t) = x_s(t + a(t)) + y(t)      a = phase, y = orbital

        with the phase--orbital CROSS term dropped.  ⚠ That is a documented
        approximation with a KNOWN SIGN, not an oversight: the paper reports
        the correlation spectrum negligible on two circuits, and that when
        present it *decreases* the total.  **Dropping it therefore OVER-states
        noise** — conservative for a design margin, wrong in a known
        direction.  It is identically zero with no AM-to-PM coupling.

        ⚠⚠ "NEGLIGIBLE" IS A PROPERTY OF THOSE CIRCUITS, and the
        over-statement can be large.  Against pnoise (the total linear
        sideband noise, confirmed by a Monte Carlo of the SDE), on van der
        Pol with an `a u^2` asymmetry, the sum is right to 0.1 % on a
        SYMMETRIC orbit and over-states by up to 3.2x (5 dB) at half-wave
        asymmetry 0.1, with no grid dependence.

        ⚠⚠ THE DROPPED CROSS TERM IS THE CAUSE -- WITH EVERY HARMONIC KEPT.
        `S_corr` from eq (92) keeps only the PPV's DC harmonic AT THE NOISE
        SOURCE's row (~1e-8 of the total where a tank inductor shorts that
        node at DC); the full-harmonic correlation is -1.1 to -2.4x this
        spectrum at `a = 0.30`, and `modal_spectrum`'s three terms sum to
        pnoise.  The over-statement is in the decomposition's
        frequency-independent terms above f_amp: the phase half is the
        Lorentzian's frequency-independent PPV, and the orbital half
        over-states by a factor FLAT in offset -- the orbital mode's AM
        share at the output, `sin^2 arg(U_{l,1}/U_{0,1})`.  Traversa &
        Bonani's own Figs 1-2 show the same limit.  A warning fires above
        `ORBITAL_ASYMMETRY_LIMIT`.

        **Lemma 3.5**: the orbital spectrum is a sum of Lorentzians centred at
        `j*w0 + Im(mu_l)` with half-width `|Re(mu_l)| + (1/2) h^2 w0^2 c`,
        weighted by the `C_lhj` of eq (22).  Every input already exists:
        `orbital_correlation` returns `C_lhj` (gated three ways), the
        exponents come from `floquet_modes`, and `c` from
        `diffusion_constant`.

        ⚠ UNIT CONVERSION, DONE ONCE HERE.  Lemma 3.5's widths are ANGULAR.
        `(1/2) h^2 w0^2 c` rad/s is `pi h^2 f0^2 c` Hz -- exactly the
        half-width `lorentzian` already uses for the phase line -- and
        `|Re(mu_l)|` rad/s is `|Re(mu_l)|/(2 pi)` Hz.  The two half-widths
        ADD, so an orbital mode's line is the phase line broadened by the
        mode's own relaxation rate.

        ⚠⚠ AND THAT IS WHY IT MATTERS AT LARGE OFFSET, WHICH IS THE WHOLE
        POINT OF THE ITEM.  The phase line's width is `pi h^2 f0^2 c`, which
        for a good oscillator is tiny, so its skirt has fallen as `1/f^2` long
        before the orbital line -- width `|Re(mu_2)|/(2 pi)`, i.e. the
        AMPLITUDE RELAXATION RATE -- has even started to roll off.  The
        crossover therefore sits near

            f_amp = -ln(lam2) f0 / (2 pi) = f0 / (2 pi Q)

        the same pole `oscillator_spectrum` warns about from the other side.
        ⚠ Those two arrived independently -- one from a commercial
        simulator's excess over our phase-only answer, one from this paper's
        modal sum -- and they must land in the same place.  That is the gate
        (`test_the_orbital_spectrum_is_a_lorentzian_of_half_width_f_amp`: the
        `h = 0` line's half-width IS `f_amp`), and it is the check that can
        actually fail.

        ⚠ `output` follows `oscillator_spectrum`: an integer indexes the
        REDUCED state (the reference row already removed), an array is a
        weight vector over it.

        ⚠ STATIONARY WHITE SOURCES ONLY -- inherited from
        `orbital_correlation`, which needs `CY` constant for eq (22)'s
        products to collapse.

        A NOISELESS circuit has no lines at all and returns zeros (it was
        refused as "no orbital line" until 2026-10-01).

        History: `doc/shooting_history.md`, `PAC.orbital_spectrum`.
        """
        output = output_index(pss, output)
        self._check_circuit(pss)
        self._refuse_driven(pss, 'orbital_spectrum')
        harmonic = integer_arg(harmonic, 'harmonic', 'orbital_spectrum',
                               minimum=1, why='the line is a carrier\'s')
        try:
            _asym = self._orbit_asymmetry(pss)
        except Exception:
            _asym = 0.0
        if _asym > self.ORBITAL_ASYMMETRY_LIMIT:
            ## ⚠ NOT the grid residual the deleted `_warn_if_orbit_is_asymmetric`
            ## named (`doc/shooting_history.md`):
            ## the sum this spectrum is meant for over-states the TOTAL on an
            ## asymmetric orbit, and no refinement changes it -- see the
            ## docstring.
            warn(
                'PAC.orbital_spectrum: this orbit has half-wave asymmetry '
                '%.3f. On an asymmetric orbit S_ph + S_orb over-states the '
                'total sideband noise above f_amp: x1.03 / x1.45 / x3.2 at '
                'asymmetry 0.033 / 0.067 / 0.100 (van der Pol, C=4, Q=8, '
                '10 f_amp), confirmed by Monte Carlo, which agrees with '
                'pnoise to 1 %%. Refining the grid does not change it. The '
                'cause is the phase-orbital correlation this sum drops (with '
                'every harmonic kept it is -1.1 to -2.4x the orbital term '
                'there): use PAC.modal_spectrum for a phase/orbital/'
                'correlation split that sums to the total, or PAC.pnoise for '
                'the total.'
                % (_asym,), ModelWarning)
        _R, C = self._orbital_correlation(pss, maxharmonics, 'orbital_spectrum')
        modes = pss.floquet_modes()
        c = float(self._diffusion_constant(pss, 'orbital_spectrum'))
        f0 = 1.0 / float(pss.period)
        m = pss.cir.n - 1

        row = _output_row(output, m)

        ## ⚠⚠ NO ORBITAL LINE AT THIS HARMONIC -- refuse rather than return the
        ## tails of the others.  The line weight at `j f0` is
        ## `W_j = sum_{l,h} Re(row C_lhj row)`; where it is zero (a symmetric
        ## orbit's even harmonics: the modes' own Fourier content vanishes)
        ## what this would return is the neighbouring lines' Lorentzian tails
        ## (3.2x LOW at 2 f0 on van der Pol against a Monte Carlo).  ⚠ NOT
        ## caught: DC, where a small line can exist and the model reads ~100x
        ## HIGH (the tank inductor shorts the node, which Lorentzian tails do
        ## not know), and 2 f0 on an asymmetric orbit -- away from the
        ## fundamental use `pnoise`.
        _W = {}
        for (_l, _h, _j), _cl in C.items():
            _W[_j] = _W.get(_j, 0.0) + float(np.real(row @ _cl @ row))
        _Wtot = sum(abs(v_) for v_ in _W.values())
        f = float(harmonic) * f0 + np.atleast_1d(
            np.asarray(offsets, dtype=float))
        if _Wtot == 0.0:
            ## no source reaches the output: a density of zero, not a refusal
            return np.zeros_like(f, dtype=float)
        if abs(_W.get(harmonic, 0.0)) <= 1e-9 * _Wtot:
            raise ValueError(
                'PAC.orbital_spectrum: no orbital line at harmonic %d for this '
                'output (weight %.3e of %.3e), so the value here would be the '
                'tails of other lines, not the noise at that frequency '
                '(measured 3.2x low at 2 f0 on a symmetric van der Pol). Use '
                'PAC.pnoise there.' % (harmonic, _W.get(harmonic, 0.0),
                                       _Wtot))

        S = np.zeros_like(f, dtype=float)
        for (l, h, j), Clhj in C.items():
            ## The weight is the output's own share of this term.  It is real
            ## for the total (`R` is real symmetric); an individual `(l,h,j)`
            ## can carry a small imaginary part that cancels against its
            ## conjugate partner, so take the real part per term rather than
            ## asserting each is real.
            w = float(np.real(row @ Clhj @ row))
            if w == 0.0:
                continue
            mul = modes[l]['mu']
            ## Hz, both terms -- see the unit note above.
            gam = abs(float(np.real(mul))) / (2.0 * np.pi) \
                + np.pi * float(h) ** 2 * f0 ** 2 * c
            fc = float(j) * f0 + float(np.imag(mul)) / (2.0 * np.pi)
            if gam <= 0.0:
                continue
            ## Normalised Lorentzian: integrates to 1 over all `f`, so the
            ## total power is `sum(w) = row^T R row` by construction.
            S = S + w * (gam / np.pi) / ((f - fc) ** 2 + gam ** 2)
        ## (`R` is the two-sided covariance: doubled, exactly, one-sided)
        return 2.0 * S

    def modal_spectrum(self, pss, offsets, output, harmonic=1,
                       maxharmonics=None, maxsidebands=None):
        """Phase, orbital AND phase-orbital CORRELATION spectra from ONE modal
        transfer, which sum to the total.

        Returns a dict of arrays at `harmonic*f0 + offsets` (a negative offset
        is the lower sideband), each a ONE-SIDED PSD at its own absolute
        frequency, the scale of `pnoise`, `oscillator_spectrum`'s `S_v` and
        `orbital_spectrum` (since 2026-09-28; they were half that).  The
        sidebands are NOT symmetric about the carrier (the lower 2.45x the
        upper on an asymmetric van der Pol, as `pnoise` there): each keeps
        its own value.  `output`: a reduced index, a weight vector, or a
        node name.  History: `doc/shooting_history.md`, `PAC.modal_spectrum`.

            'phase', 'orbital', 'correlation', 'total'
            total = phase + orbital + correlation

        Every Floquet mode `l` -- the phase mode (`mu = 0`) and each orbital
        mode -- carries noise from input sideband `m` to the output at `w`:

            T_m^l(w) = sum_j (d . U_{l,j}) V_{l,m-j}^T / (i(w - j w0) - mu_l + a_j)

        with `U`, `V` the Fourier coefficients of `p_l`, `q_l` (the
        conventions of `orbital_correlation`) and `a_j = j^2 w0^2 c / 2` the
        phase-diffusion rate of output harmonic `j`.  With `T = T^0 + sum_l
        T^l` the output is `sum_m T_m (CY/2) T_m^H`; `phase`, `orbital` and
        `correlation` are its phase-phase, orbital-orbital and
        `2 Re(phase-orbital)` blocks.

        ⚠⚠ WHY IT EXISTS.  `oscillator_spectrum(frequency_aware=False) +
        orbital_spectrum` over-states an asymmetric orbit's total by up to
        3.2x (Monte-Carlo-confirmed).  The missing piece IS the correlation
        -- but with every harmonic kept, not in the form Traversa & Bonani
        keep: their eq (92) retains only its DC harmonic, ~1e-8 of the total
        on van der Pol (the tank inductor shorts the source node at DC).  It
        removes the DC-PPV phase excess above f_amp AND the orbital mode's PM
        projection at the output -- the orbital line's AM share
        `sin^2 arg(U_{l,1}/U_{0,1})` is the flat factor `orbital_spectrum`
        over-states by.  `total / (pnoise/2)` is 1.0006 on a symmetric van
        der Pol; on an asymmetric one the excess is the modes' O(h) grid
        error (it halves with the grid).
        ⚠ Because the correlation cancels most of the other two, a few
        percent of error in any part is AMPLIFIED in the total -- which is
        why the three are computed together here rather than the
        correlation being offered as an add-on to `oscillator_spectrum +
        orbital_spectrum`: those line-shape spectra keep only the resonant
        term of each line (2.5 % short at 10 f_amp even on a symmetric
        orbit), which is harmless alone and not under cancellation.  The
        modal sum also reproduces pnoise's upper/lower sideband asymmetry,
        which the two-term sum cannot.

        Near the carrier `phase` IS the library Lorentzian (symmetric orbit)
        and `correlation` is ~1e-6 of it.  Above f_amp `total` agrees with
        `pnoise`; within the linewidth pnoise has no meaning and this is the
        route.

        ⚠ Free-running oscillators and the dense `floquet_modes` only
        (inherited).  `harmonic >= 1`: harmonic 0 was never measured.  The
        PPV / mode harmonics kept, `H = maxharmonics`, default to
        `ORBITAL_HARMONICS` capped by the grid (an explicit count above it
        raises); the sidebands summed, `maxsidebands`, to `2 H`.  (`H` and
        `sidebands` until 2026-09-29.)  `output` follows `orbital_spectrum`.

        ⚠ A source whose LEVEL FOLLOWS THE ORBIT (a MOS
        channel's thermal and flicker noise, a shot noise): its sidebands
        are correlated, so they no longer add in power.  WHITE parts enter
        in the P-form ``sum_{m,m'} T_m (P_{m'-m}/2) T_{m'}^H`` with `P_k`
        the harmonics of `CY(x(t))` -- exact, no square root (a root of
        `(k V)^2` is `|k V|`, whose kink the sideband window would cut);
        `c` is Demir's `B(x(t))` form.  Each COLOURED component is a unit
        process through its own columns `G(t)` -- its signed amplitudes
        where the element states them, else the root of its PSD (warned
        where that PSD touches zero) -- and its rows are the harmonics of
        `q_l^T G`, band `p` weighted by the colour at ``|w - p w0|``
        (`_modal_modulated`, `NoiseComponents.colour_groups`).  Gated against the same
        physics built as a stationary source times the modulating voltage:
        every part to ~1e-13.

        ⚠ A COLOURED SOURCE: input sideband `m` carries the
        source at ``w - m w0``, and reads `CY` THERE (as `pnoise` does)
        rather than one `CY` for every sideband; the phase-diffusion widths
        `a_j` take `c` from the WHITE part of the sources alone (zero for a
        pure 1/f or Lorentzian source).  The phase part is then the
        linearised skirt, valid only where `phase_psd` is: offsets it refuses
        are refused here (`_modal_colour`).

        History: `doc/shooting_history.md`, `PAC.modal_spectrum`.
        """
        return self._modal_spectrum(pss, offsets, output, harmonic,
                                    maxharmonics, maxsidebands,
                                    'modal_spectrum')

    def _modal_spectrum(self, pss, offsets, output, harmonic, maxharmonics,
                        maxsidebands, what):
        """`modal_spectrum`, its refusals naming `PAC.<what>`, the method the
        user called."""
        output = output_index(pss, output)
        self._check_circuit(pss)
        coloured = self._coloured_present(pss)
        self._refuse_driven(pss, what)
        harmonic = integer_arg(harmonic, 'harmonic', what,
                               minimum=1, why='DC was never measured against '
                               'pnoise; use PAC.pnoise there')
        offs = np.atleast_1d(np.asarray(offsets, dtype=float))
        ## ⚠ A MODULATED source (its `CY` follows the orbit) takes its own
        ## sum, `_modal_modulated`; a stationary circuit runs the code
        ## below as before
        modulated = self._modulated_present(pss)
        if modulated:
            ## the band reach, for the per-band sources' classification
            N_ = len(self._ppv_states(pss))
            H_ = self._harmonic_count(maxharmonics, N_, what)
            M_ = self._sideband_count(maxsidebands, H_, N_, what)
            c_mod, P2, groups = self._modal_modulated(
                pss, offs, harmonic, coloured, M_ + harmonic,
                what)
        elif coloured:
            c_col, cy_at = self._modal_colour(pss, offs, harmonic,
                                              what)
        modes = pss.floquet_modes()
        ## the phase mode by its tangent alignment, on any grid -- see
        ## `_phase_mode_split` (a 1e-6 window on |lam| - 1 refuses every gear
        ## solve on a non-uniform grid)
        _kph, orb = self._phase_mode_split(pss, modes, 'PAC.' + what)
        ph = [_kph]
        m = pss.cir.n - 1
        row = _output_row(output, m)
        if modulated:
            c = c_mod
        else:
            c = c_col if coloured else float(
                self._diffusion_constant(pss, what))
        w0 = 2.0 * np.pi / float(pss.period)
        CY2 = (None if coloured or modulated
               else 0.5 * np.real(np.asarray(self._cy_reduced(pss, 0.0))))
        N = np.asarray(modes[ph[0]]['p']).shape[1] - 1
        H = self._harmonic_count(maxharmonics, N, what)
        M = self._sideband_count(maxsidebands, H, N, what)
        js = np.arange(-H, H + 1)
        ms = np.arange(-M, M + 1)
        a_j = 0.5 * js.astype(float) ** 2 * w0 ** 2 * c

        ## per mode: the output's share of each harmonic of p_l, and q_l's
        ## Fourier coefficients (m x N).  ⚠ The phase mode's exponent is set
        ## to 0 exactly: its multiplier is 1 to rounding, and a 1e-16 real part
        ## would put a spurious pole width on the Lorentzian.
        coef = []
        ## modulated: per mode, q_l's samples and, per FIXED coloured group,
        ## the Fourier coefficients of `q_l^T G(t)` -- the modulation enters
        ## as a product in time, i.e. the convolution
        ## `R_p = sum_k T_{p-k} G_k` over every harmonic the grid holds
        extra = []
        for l in ph + orb:
            ## ⚠ `_period_dft`, not an index DFT, which is 8-13 % off and does
            ## not converge on a 3:1 grid (see `PSS._period_quadrature`)
            Ul = self._period_dft(pss, np.asarray(modes[l]['p'])[:, :-1].T).T
            Vl = self._period_dft(pss, np.asarray(modes[l]['q'])[:, :-1].T).T
            coef.append((l, row @ Ul[:, js % N], Vl,
                         0.0 if l == ph[0] else complex(modes[l]['mu'])))
            if modulated:
                ql = np.asarray(modes[l]['q'])[:, :-1]
                extra.append((ql, [
                    self._period_dft(pss, np.einsum('mn,nmr->nr', ql, G)).T
                    for kind, G, _s in groups if kind == 'fixed']))

        def transfer(w, entry):
            ## (2M+1) x m; looped over j so memory stays m x (2M+1) per mode
            _l, u, Vl, mul = entry
            g = u / (1j * (w - js * w0) - mul + a_j)
            ## columns: the source rows `Vl` carries (`m`; `r` for a group)
            ## ⚠ NO HARMONIC PAST THE GRID'S NYQUIST: the coefficient `m - j`
            ## with `|m - j| > N/2` is not on the grid -- the index wrapped
            ## `% N` onto a LOW harmonic until 2026-10-01, 51 % high at
            ## N = 64 on an asymmetric van der Pol (every grid where the
            ## default harmonic count meets the Nyquist, N <= 66); dropped,
            ## as `quadP` drops it
            T = np.zeros((ms.size, Vl.shape[0]), dtype=complex)
            for ji, j in enumerate(js):
                if g[ji] != 0.0:
                    on = (np.abs(ms - j) <= N // 2)[:, None]
                    T += g[ji] * (Vl[:, (ms - j) % N].T * on)
            return T

        def quad(A, B):
            if CY2.ndim == 3:
                return complex(np.einsum('mi,mik,mk->', A, CY2, np.conj(B)))
            return complex(np.einsum('mi,ik,mk->', A, CY2, np.conj(B)))

        def quadP(A, B):
            ## a modulated WHITE source: `sum_{m,m'} A_m (P_{m'-m}/2) B_{m'}^H`
            ## with `P_k` the harmonics of its `CY(x(t))` -- exact, no root
            ## (a root of `(k V)^2` is `|k V|`, whose kink spreads it over
            ## harmonics the sideband window then cuts).  ⚠ No circular wrap
            ## past N/2: that harmonic is not on the grid.
            S_ = A.shape[0]
            tot = 0.0j
            for d in range(-(S_ - 1), S_):
                if abs(d) > N // 2:
                    continue
                lo, hi = max(0, -d), min(S_, S_ - d)
                tot += complex(np.einsum('mi,ik,mk->', A[lo:hi], P2[d % N],
                                         np.conj(B[lo + d:hi + d])))
            return tot

        def coloured_rows(w, gi, kind, G):
            ## a modulated COLOURED group: the rows `R_p` of the unit
            ## process in band `p` (at `|w - p w0|`), phase mode and the
            ## orbital modes' sum, `(2M+1) x r` each
            if kind == 'fixed':
                Rp = transfer(w, coef[0][:2] + (extra[0][1][gi],) + coef[0][3:])
                Ro = np.zeros_like(Rp)
                for k in range(1, len(coef)):
                    Ro += transfer(w, coef[k][:2] + (extra[k][1][gi],)
                                   + coef[k][3:])
                return Rp, Ro
            ## the root per band: `q_l^T G(t; nu_p)`, row `p` alone
            rows = []
            for k in range(len(coef)):
                _l, u, _Vl, mul = coef[k]
                g = u / (1j * (w - js * w0) - mul + a_j)
                R = []
                for p in ms:
                    Wp = self._period_dft(pss, np.einsum(
                        'mn,nmr->nr', extra[k][0], G(abs(w - float(p) * w0)))).T
                    ## (no harmonic past the Nyquist, as in `transfer`)
                    R.append(Wp[:, (p - js) % N] @ (g * (np.abs(p - js) <= N // 2)))
                rows.append(np.asarray(R))
            return rows[0], sum(rows[1:], np.zeros_like(rows[0]))

        res = {k: np.zeros(offs.shape, dtype=float)
               for k in ('phase', 'orbital', 'correlation', 'total')}
        for i, o in enumerate(offs.ravel()):
            w = float(harmonic) * w0 + 2.0 * np.pi * float(o)
            if coloured and not modulated:
                ## input sideband `m` is the source at `w - m w0`
                CY2 = np.array([0.5 * np.real(np.asarray(
                    cy_at(abs(w - float(mm) * w0)))) for mm in ms])
            Tp = transfer(w, coef[0])
            To = np.zeros_like(Tp)
            for entry in coef[1:]:
                To += transfer(w, entry)
            if modulated:
                sp = float(np.real(quadP(Tp, Tp)))
                so = float(np.real(quadP(To, To)))
                sc = 2.0 * float(np.real(quadP(Tp, To)))
                nu = np.abs(w - ms * w0)
                gi = 0
                for kind, G, s in groups:
                    Rp, Ro = coloured_rows(w, gi, kind, G)
                    if kind == 'fixed':
                        gi += 1
                        wt = 0.5 * np.asarray(s(nu), dtype=float)
                    else:
                        wt = 0.5 * np.ones(ms.size)
                    sp += float(np.sum(wt * np.sum(np.abs(Rp) ** 2, axis=1)))
                    so += float(np.sum(wt * np.sum(np.abs(Ro) ** 2, axis=1)))
                    sc += 2.0 * float(np.sum(wt * np.real(
                        np.sum(Rp * np.conj(Ro), axis=1))))
            else:
                sp = float(np.real(quad(Tp, Tp)))
                so = float(np.real(quad(To, To)))
                sc = 2.0 * float(np.real(quad(Tp, To)))
            ix = np.unravel_index(i, offs.shape)
            res['phase'][ix] = sp
            res['orbital'][ix] = so
            res['correlation'][ix] = sc
            res['total'][ix] = sp + so + sc
        ## (a real output: `S(-f) = S(f)`, so the one-sided PSD is twice the
        ## two-sided one at every frequency, both sidebands -- exact)
        return {k: 2.0 * v for k, v in res.items()}

    def _phase_psd_gate(self, pss, offs, harmonic, what):
        """With a COLOURED source the modal spectra's phase part is the
        linearised skirt: refuse the offsets where `phase_psd` would (its
        corner and power bound), naming the caller.  Returns `|offs|`."""
        ao = np.abs(np.asarray(offs, dtype=float)).ravel()
        try:
            self.phase_psd(pss, np.unique(ao), harmonic=int(harmonic),
                           frequency_aware=False)
        except ValueError as e:
            raise ValueError(
                'PAC.%s: with a COLOURED source the phase part is the '
                'linearised skirt, valid only where phase_psd is -- %s'
                % (what, e)) from None
        return ao

    def _modal_colour(self, pss, offs, harmonic, what):
        """For `modal_spectrum` on a coloured circuit: `(c_white, cy_at)`.
        STATIONARY sources (a modulated one takes `_modal_modulated`;
        `_cy_reduced` still guards this sum).  Refuses the offsets where the
        linearised phase has broken down (`phase_psd`'s corner and power
        bound); `c_white` is the diffusion of the sources' WHITE part alone
        (`NoiseComponents.model`); `cy_at(w)` the reduced `CY` at `w`."""
        f0 = 1.0 / float(pss.period)
        w0 = 2.0 * np.pi * f0
        self._cy_reduced(pss, w0)
        self._phase_psd_gate(pss, offs, harmonic, what)
        x0r = np.asarray(pss._period_state[1], dtype=float).ravel()
        irn = pss.irefnode
        x0f = insert_ref(x0r, irn)
        with warnings.catch_warnings():
            ## (`model`'s one cost note, by its category -- a text filter
            ## until 2026-10-01)
            warnings.simplefilter('ignore', CostWarning)
            model = self._noise_components(pss, [x0f]).model(1e-3 * f0, f0)
        if model is None:
            raise NotImplementedError(
                'PAC.%s: this circuit\'s CY is not the sum of its elements\' '
                'and, as a whole, not thermal-plus-power-law, so the white part '
                'of its sources cannot be separated for the phase-diffusion '
                'widths.' % what)
        c_white = float(self._white_diffusion_at(
            pss, w0, cy=np.real(np.asarray(model.white[0]))))
        return c_white, (lambda w: self._cy_at(pss, w, x0r))

    def _modal_modulated(self, pss, offs, harmonic, coloured, L, what):
        """For `modal_spectrum` with a MODULATED source: `(c_white, P2,
        groups)` -- `c_white` Demir's `c` of the WHITE parts at their own
        states, `P2` half the harmonics of the white `CY(x(t))` (the
        P-form), `groups` the coloured components (`colour_groups`).  A
        coloured circuit is first held to `phase_psd`'s validity, as
        `_modal_colour` does."""
        f0 = 1.0 / float(pss.period)
        w0 = 2.0 * np.pi * f0
        states = self._ppv_states(pss)
        if not coloured:
            white = np.real(self._noise_components(pss, states).cy_at_states(w0))
            groups = []
        else:
            ao = self._phase_psd_gate(pss, offs, harmonic, what)
            nc = self._noise_components(pss, states)
            model = nc.model(1e-3 * f0, f0)
            if model is None:
                raise NotImplementedError(
                    'PAC.%s: this circuit\'s CY is not the sum of its '
                    'elements\' and, as a whole, not thermal-plus-power-law, '
                    'so its white and coloured parts cannot be separated.  '
                    'PAC.pnoise(cyclostationary=True) gives the total.' % what)
            white = np.real(np.asarray(model.white))
            groups = nc.colour_groups(model, 2.0 * np.pi * float(np.min(ao)),
                                      f0, L, what)
        c_white = float(self._white_diffusion_at(pss, w0, cy=white))
        return c_white, 0.5 * self._period_dft(pss, white), groups

    def correlation_spectrum(self, pss, offsets, output, harmonic=1,
                             maxharmonics=None, maxsidebands=None):
        """The FULL-harmonic phase-orbital correlation spectrum:
        `modal_spectrum(...)['correlation']`, at `harmonic*f0 + offsets` (a
        negative offset is the LOWER sideband; the term can change sign
        between the two).

        ⚠ It sums to the total with `modal_spectrum`'s own `phase` and
        `orbital` -- NOT with `oscillator_spectrum + orbital_spectrum`, whose
        line-shape approximations are a few percent off under the cancellation
        this term produces (see `modal_spectrum`).  Negative where AM-to-PM
        coupling exists: -1.1 to -2.4x the orbital term on van der Pol at
        half-wave asymmetry 0.10.
        """
        return self._modal_spectrum(pss, offsets, output, harmonic,
                                    maxharmonics, maxsidebands,
                                    'correlation_spectrum')['correlation']

    #: the phase mode may sit this far off the unit circle before it is refused
    PHASE_MODE_MAX_DEPARTURE = 1e-3

    def _phase_mode_split(self, pss, modes, where):
        """`(phase_index, orbital_indices)` -- the phase mode identified by WHAT
        DEFINES IT, its eigenvector's alignment with the orbit tangent, not by
        a window on `|lam| - 1`.

        On a uniform grid the phase multiplier is 1 to rounding.  On a grid
        whose step varies, a multistep or trapezoidal solve loses time-
        translation symmetry and the multiplier leaves the circle at O(h^2)
        (radau keeps it at 1 to ~1e-11), so a window on `|lam| - 1` would
        refuse gear there.  The right eigenvector of the phase mode is the
        tangent `xdot(0)` (`C xdot = -i(x)` for the autonomous circuit),
        which no orbital mode shares, so alignment picks it on any grid;
        its exponent is then forced to 0 exactly, as the consumers already
        do.  The departure is WARNED with its size when it exceeds rounding,
        and the split is REFUSED when a second multiplier lies within ten
        times that departure of the circle with any alignment -- the case a
        window ever protected against.

        History: `doc/shooting_history.md`, `PAC._phase_mode_split`.
        """
        m = pss.cir.n - 1
        irn = pss.irefnode
        ## ⚠ `waveform` is FULL width (the reference row is in it): a zero
        ## inserted a second time reads every unknown past the reference one
        ## slot late
        ## History: `doc/shooting_history.md`, `PAC._phase_mode_split`.
        xf = orbit_states(pss, [np.asarray(pss.waveform[1],
                                           dtype=float)[:, 0]])[0]
        xr = np.delete(xf, irn)
        ## ⚠ `C xdot = -(i + u)`: the SOURCE too -- a DC source on a row with
        ## capacitance is part of the tangent (without it a 1 A current on
        ## the tank node read `xdot` 56 % off, the phase mode's alignment
        ## 0.83 under the 0.9 bar, and `orbital_correlation` REFUSED the
        ## oscillator; until 2026-10-01, the review's D1).  As `ppv()` forms it.
        with devices_at(pss.cir, xf, pss.epar):
            i_red = np.delete(np.asarray(pss.cir.i(xf, pss.epar), dtype=float).ravel()
                              + np.asarray(pss.cir.u(0.0, epar=pss.epar,
                                                     analysis=pss.par.analysis),
                                           dtype=float).ravel(), irn)
        C0 = np.asarray(pss._C_at(xr), dtype=float)
        try:
            xdot0 = np.linalg.solve(C0, -i_red)
        except np.linalg.LinAlgError:
            xdot0 = np.linalg.lstsq(C0, -i_red, rcond=None)[0]
        nx = float(np.linalg.norm(xdot0))
        dep, cos = [], []
        for md in modes:
            u = np.asarray(md['u0'])[:m]
            dep.append(abs(abs(complex(md['lam'])) - 1.0))
            cos.append(abs(complex(np.vdot(u, xdot0))) / max(float(np.linalg.norm(u)) * nx, 1e-300))
        cand = [k for k in range(len(modes)) if dep[k] <= self.PHASE_MODE_MAX_DEPARTURE and cos[k] > 0.9]
        if not cand:
            raise ValueError(
                '%s: no Floquet mode is both within %.0e of the unit circle and '
                'aligned with the orbit tangent (|lam|-1: %s; alignment: %s).  Is '
                'this an autonomous oscillator solved at its own period?'
                % (where, self.PHASE_MODE_MAX_DEPARTURE,
                   ', '.join('%.1e' % d for d in dep), ', '.join('%.2f' % c for c in cos)))
        k = max(cand, key=lambda j: cos[j])
        rival = [j for j in range(len(modes)) if j != k and dep[j] <= max(10.0 * dep[k], 1e-6)]
        if rival:
            raise ValueError(
                '%s: a second Floquet multiplier (|lam|-1 = %s) lies as close to '
                'the unit circle as the phase mode (%.1e), so the phase mode '
                'cannot be told from an orbital one.  A uniform grid or '
                'method=\'radau\' puts the phase multiplier at 1 to rounding.'
                % (where, ', '.join('%.1e' % dep[j] for j in rival), dep[k]))
        if dep[k] > 1e-6:
            warn(
                '%s: the phase mode sits %.1e off the unit circle (alignment '
                'with the orbit tangent %.4f); its exponent is forced to 0.  '
                'On a non-uniform grid a multistep or trapezoidal solve leaves '
                'it there at O(h^2); the modal parts are then the method\'s '
                'order (measured second order for gear on a 3:1 grid).'
                % (where, dep[k], cos[k]), AccuracyWarning)
        return k, [j for j in range(len(modes)) if j != k]
