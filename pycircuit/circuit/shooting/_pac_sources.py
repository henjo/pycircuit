"""The noise sources as PAC reads them: the component model's factory
(`_noise_components`; the model itself is `_noise_components.py`), `CY`
at one state and cycle-averaged, the stationarity and colour guards, and
the period quadrature.
"""
import numpy as np
from pycircuit.circuit.analysis import remove_row_col
from ._noise_components import NoiseComponents, orbit_states
from ._numerics import insert_ref, periodic_spline_weights


class _NoiseSources(object):
    """The noise sources as PAC reads them: the component model's factory,
    `CY` at one state and cycle-averaged, the guards, and the period
    quadrature.  A theme of `PAC` (see `pac.py`)."""

    def _period_dft(self, pss, S):
        """Fourier coefficients over the period of samples `S` `(N, ...)` taken
        at `fp.times[:N]`, laid out like `numpy.fft.fftfreq`.  Uniform grid:
        the index DFT, unchanged.  Non-uniform: the weighted sum at the TRUE
        times with `PSS._period_quadrature`'s weights -- O(N^2), paid only by
        a caller who chose that grid."""
        S = np.asarray(S, dtype=complex)
        fp = pss.factored_period()
        wq = pss._period_quadrature(fp)
        if wq is None or S.shape[0] != len(wq):
            return np.fft.fft(S, axis=0) / S.shape[0]
        N = S.shape[0]
        tms = np.asarray(fp.times, dtype=float)
        ks = np.fft.fftfreq(N, d=1.0 / N)
        E = np.exp(-2j * np.pi * np.outer(ks, (tms[:N] - tms[0]) / float(tms[N] - tms[0]))) * wq[None, :]
        return np.tensordot(E, S, axes=(1, 0))

    def _noise_components(self, pss, states=None):
        """The noise sources of `pss`'s circuit at `states` (None: the
        orbit samples), built per call: `_noise_components.NoiseComponents`,
        the component model, the per-element shares and roots.  A test
        hands its own subclass in here."""
        return NoiseComponents(pss, states)

    def _cy_cycle_averaged(self, pss, w):
        """`CY` time-averaged over the orbit — Hull & Meyer's construction.

        ⚠⚠ VALID FOR GENTLE MODULATION ONLY, AND IT FAILS AS A FACTOR, NOT A
        PERCENTAGE: on a series switch + shunt capacitor it reads 16x a
        reference simulator at `goff/gon` = 1e-6 (1.000 unmodulated).  The
        averaged source injects `4kT <g>` for the WHOLE period, including
        the hold phase, where Hull & Meyer's own condition fails by six
        orders.  So this route is for a mixer's `gm` or a bias-dependent
        shot noise -- not for a switch.  A switch's noise is reachable
        exactly: `covariance` and `oscillator_covariance` evaluate `CY` at
        every step and need no averaging.

        This is what `_cy_reduced` refuses, done instead of refused, after
        Hull & Meyer (1993): *"cyclostationary noise sources, such as shot
        noise, may be modeled as MODULATED STATIONARY NOISE SOURCES.  The
        impulse response that is calculated INCLUDES THE EFFECT OF THIS
        MODULATION"*, with the stationary source at the cycle-averaged
        value.  The modulation is carried by the RESPONSE, `H_l`, rather
        than by the sources: ONE source per device.

        ⚠ Its condition: *"NONE OF THE LARGE-SIGNAL STATE VARIABLES MAY
        CHANGE SIGNIFICANTLY OVER THE DECAY TIME OF THE IMPULSE
        RESPONSE."*  A ringing (high-Q) impulse response breaks it, so it
        degrades as `lambda_2 -> 1`; `info` reports `|lambda_2|` so the
        caller can see which regime they are in.

        ⚠ SAMPLED ON THE ORBIT, NOT AT THE OPERATING POINT: `CY` is
        evaluated at every stored state and averaged with the step weights
        -- the same quadrature `diffusion_constant` uses, so the two remain
        comparable.

        History: `doc/shooting_history.md`, `PAC._cy_cycle_averaged`.
        """
        irn = pss.irefnode
        fp = pss.factored_period()
        tms = np.asarray(fp.times, dtype=float)
        T = float(fp.T)
        xs = np.asarray(pss.waveform[1], dtype=float)
        nsamp = min(len(tms) - 1, xs.shape[1])
        hs = self._period_weights(tms, nsamp, T, pss)
        acc = None
        for k in range(nsamp):
            ## ⚠ `waveform` is FULL width: `xs[:m]` then a second zero reads
            ## every unknown past the reference one slot late
            ## History: `doc/shooting_history.md`, `PAC._cy_cycle_averaged`.
            xf = orbit_states(pss, [xs[:, k]])[0]
            cyk = np.asarray(pss.cir.CY(xf, w, epar=pss.epar), dtype=complex)
            (cyk,) = remove_row_col((cyk,), irn, pss.toolkit)
            cyk = np.asarray(cyk, dtype=complex) * hs[k]
            acc = cyk if acc is None else acc + cyk
        if acc is None:
            raise NotImplementedError(
                'PAC: the orbit has no stored samples to average CY over.')
        return acc / float(hs[:nsamp].sum())

    def _cy_reduced(self, pss, w):
        """`CY` with the reference node removed, refusing a moving one.

        ⚠ THE CHECK IS THE POINT.  A bias-dependent `CY` makes the sources
        CYCLOSTATIONARY, and then sidebands stop adding in power -- the
        stationary sum this class computes would be the wrong model, not
        merely an inaccurate one, and nothing downstream would say so.
        Sampled at three states on the converged orbit rather than argued
        from the element types, because a compact model's `CY` reads `x`
        and the discrete library's does not.

        ⚠⚠ IN EFFECT THIS REFUSES MOS pnoise: no physically correct MOS
        noise model has a state-independent `CY` (thermal `4kT gamma g_d0`
        with `g_d0` bias-dependent, flicker `I_D^AF`, gate shot `2qI_G`,
        trap rates reading the terminal voltages -- Mahmutoglu & Demir
        2015).  `PspMosLongChannel` is noiseless only by default (`fnt =
        0`).  The routes past the refusal are `cyclostationary=True` (exact
        for thermal and shot noise; flicker is the `|m|` fold, see
        `pnoise`) and `modulated=True` (Hull & Meyer's cycle average, see
        `_cy_cycle_averaged`).

        ⚠ This same check keeps the Ito/Stratonovich choice out of reach
        (`CY = GG^T`, so a state-dependent `CY` is a state-dependent `G`).
        Past it the two interpretations diverge, and Demir's escape -- "the
        noise signals are small compared with the deterministic signals" --
        may NOT carry for trap noise (a two-state Markov chain, not a small
        perturbation).  The tell would be a discrepancy in a MEAN but not in
        a variance.  (The surfaces that read `CY` state by state --
        `covariance`, `pnoise(cyclostationary=True)`, `diffusion_constant`,
        `modal_spectrum` -- report DIFFUSION, where the two agree, and no
        mean.)

        History: `doc/shooting_history.md`, `PAC._cy_reduced`.
        """
        irn = pss.irefnode
        fp = pss.factored_period()
        ## ⚠ THREE STATES ON THE ORBIT: the stored state half a period in,
        ## `x_last` and `x_prev`.  Never the zero vector -- it is on the
        ## orbit only by accident, and a switch model reading `goff` at
        ## v(ck) = 0 would refuse a linear time-invariant RC held by a DC
        ## clock as cyclostationary.
        _W = np.asarray(pss.waveform[1], dtype=float)
        _mid = np.delete(_W[:, _W.shape[1] // 2], irn, axis=0)
        states = [np.asarray(_mid, dtype=float).ravel()[:pss.cir.n - 1],
                  np.asarray(fp.x_last, dtype=float).ravel(),
                  np.asarray(fp.x_prev, dtype=float).ravel()[:pss.cir.n - 1]]
        mats = []
        for xr in states:
            xf = insert_ref(xr, irn)
            cy = np.asarray(pss.cir.CY(xf, w, epar=pss.epar), dtype=complex)
            (cy,) = remove_row_col((cy,), irn, pss.toolkit)
            mats.append(np.asarray(cy, dtype=complex))
        scale = max(float(np.max(np.abs(mats[0]))), 1e-300)
        for other in mats[1:]:
            drift = float(np.max(np.abs(other - mats[0]))) / scale
            if drift > 1e-9:
                raise NotImplementedError(
                    'PAC: this circuit has a BIAS-DEPENDENT CY (varies by '
                    '%.3g over the orbit), so its noise sources are '
                    'cyclostationary. The sidebands are then correlated '
                    'through the window Fourier coefficients and no longer '
                    'add in power, so the stationary sum here would be the '
                    'wrong model rather than an imprecise one. '
                    '⚠ THIS IS THE NORMAL CASE FOR ANY COMPACT MOS MODEL, '
                    'not an exotic one: thermal channel noise is '
                    '4kT.gamma.gd0 with gd0 bias-dependent, flicker goes '
                    'as I_D^AF, gate shot noise as 2qI_G, and trap capture '
                    'and emission rates read the terminal voltages. So '
                    'this is in effect a refusal of MOS pnoise, and the '
                    'route out is the CYCLOSTATIONARY construction rather '
                    'than a different device model -- BUILT: pass '
                    'cyclostationary=True (the PSD\'s own harmonics, exact; '
                    'modulated=True is the cycle average, which drops the '
                    'sideband correlation and read 0.53 of the truth on a '
                    'driven multiplier). Hull & Meyer (1993) '
                    'make it affordable -- ONE stationary source per '
                    'device at the CYCLE-AVERAGED current, with the '
                    'modulation carried by the impulse response, valid '
                    'while no large-signal state variable changes much '
                    'over the impulse response decay time. ⚠ THAT ROUTE '
                    'IS BUILT: pass modulated=True to use it. It is a '
                    'MODEL CHOICE with the validity condition above, not '
                    'a tolerance relaxation, which is why it is opt-in '
                    'and why this refusal is the default.'
                    % drift)
        return mats[0]

    ## How close to a harmonic of `f0` counts as "on" it, as a fraction of
    ## `f0` (`HARMONIC_GUARD`, below).  The deflated solve's conditioning is
    ## FLAT down to 1e-9 of `f0`, so the guard excludes only what has no
    ## finite answer: at an EXACT harmonic `1/(1 - alpha)` is a division by
    ## zero and the response is genuinely unbounded.
    ## History: `doc/shooting_history.md`, `PAC.HARMONIC_GUARD`.
    def _cy_at(self, pss, w, xr):
        """`CY` at ONE reduced state `xr`, reference row/column removed, with
        NO cyclostationarity check -- for the routes that evaluate the
        source at every step and therefore model a modulated source
        exactly (`covariance`, `oscillator_covariance`).  The stationary
        sum in `pnoise` cannot, which is why `_cy_reduced` refuses there.
        """
        irn = pss.irefnode
        xr = np.asarray(xr, dtype=float).ravel()[:pss.cir.n - 1]
        xf = insert_ref(xr, irn)
        cy = np.asarray(pss.cir.CY(xf, w, epar=pss.epar), dtype=complex)
        (cy,) = remove_row_col((cy,), irn, pss.toolkit)
        return np.asarray(cy, dtype=complex)

    def _refuse_coloured(self, pss, what, instead=None):
        """Refuse a coloured source where the machinery assumes WHITE.
        (`covariance` and `event_jitter` take a band instead and do not
        come here with one -- `_coloured_prepare`.)

        ⚠ THE TRAP IS THAT NOTHING ELSE WOULD OBJECT. The Lyapunov
        recursion, `diffusion_constant` and eq (22)'s collapse all read
        `CY` at ONE frequency and treat it as the noise intensity at every
        frequency; a coloured source folded that way returns a plausible
        number, not an error (A4d names exactly this shape). Detected by
        evaluating the reduced `CY` at two frequencies -- colour is
        frequency dependence, bias dependence is what `_cy_reduced`
        refuses separately.
        """
        if self._coloured_present(pss):
            raise NotImplementedError(
                'PAC.%s: a noise source in this circuit is COLOURED (its CY '
                'differs between w0 and 10 w0), and this routine assumes '
                'white sources -- it would fold CY at one frequency as if '
                'it held at every frequency and return a plausible wrong '
                'number. Use the frequency-resolved surfaces (pnoise, '
                'sampled_variance, phase_psd/coloured_diffusion on an '
                'oscillator, where 1/f noise makes the phase growth '
                'non-diffusive), covariance/event_jitter with a band '
                '(fmin, fmax) on a driven circuit, or the white-through-filter '
                'form of the source.%s' % (what, (' Here: ' + instead)
                                            if instead else ''))

    def _coloured_present(self, pss):
        """Whether any noise source of the circuit is COLOURED -- see
        `_refuse_coloured`, which asks this."""
        w1 = 2.0 * np.pi / float(pss.period)
        ## ⚠ Colour is asked at fixed state, two frequencies: it is separable
        ## from the bias question, and asking it through `_cy_reduced` would
        ## refuse every MODULATED source the covariance routes (which
        ## evaluate `CY` per step) handle exactly.
        ## ⚠⚠ PER ENTRY, AND AT MORE THAN ONE STATE: each entry is judged
        ## against ITS OWN magnitude (a drain's flicker colour sits orders
        ## below a gate resistor's white 4kT/rg in the same matrix), at
        ## `x_last` and at states spread over the orbit (a switch OFF at
        ## t = 0 is white there and coloured elsewhere).  Exact zeros and
        ## white entries give identical values at both frequencies, so
        ## neither can fire.
        ## `_cy_at` takes a REDUCED state: the reference row is removed here.
        ## History: `doc/shooting_history.md`, `PAC._refuse_coloured`.
        _xl = np.asarray(pss.factored_period().x_last, dtype=float).ravel()
        _states = [_xl]
        _wf = getattr(pss, 'waveform', None)
        if _wf is not None:
            _W = np.delete(np.asarray(_wf[1], dtype=float), pss.irefnode,
                           axis=0)
            for _k in sorted(set(np.linspace(0, _W.shape[1] - 1,
                                             8).astype(int))):
                _states.append(_W[:, _k])
        for _xr in _states:
            c1 = self._cy_at(pss, w1, _xr)
            c2 = self._cy_at(pss, 10.0 * w1, _xr)
            den = np.maximum(np.abs(c1), np.abs(c2))
            if np.any(np.abs(c1 - c2) > 1e-9 * den):
                return True
        return False

    @classmethod
    def _power_law_weights(cls, nus, ef, richardson=True):
        """Weights `q_i` with ``int nu^-ef g(nu) dnu = sum_i q_i g(nu_i)``
        EXACT for `g` linear in ``ln nu`` between the points -- the power law
        integrated analytically, the response interpolated.  On ``[u_i,
        u_i + h]`` in ``u = ln nu``, with ``lam = 1 - ef``, ``z = lam h``:
        ``nu_i^lam h (phi1(z) - phi2(z))`` to the left point and ``nu_i^lam h
        phi2(z)`` to the right, ``phi1 = (e^z - 1)/z``, ``phi2 = int_0^1 t
        e^{z t} dt``, by series near ``z = 0``.  `ef = 1` is the trapezoid in
        `ln nu`; `ef = 0` is a density sampled at the points.

        ⚠ AND RICHARDSON ON TOP, on an even number of intervals: ``(4 Q_h -
        Q_2h) / 3`` with ``Q_2h`` on every other point.  The product rule is
        exact for a flat response but carries ``h^2 g''`` where the response
        BENDS, and at `ef != 1` that term does not telescope to the band
        ends (it does at `ef = 1`): measured on an RC under 1/f^0.8, +6.85e-5
        / +1.71e-5 / +4.28e-6 / +1.07e-6 at 20 / 40 / 80 / 160 per decade --
        4.0x per halving, a clean `h^2`, which the combination removes.  The
        weights stay positive (Simpson's pattern at `ef = 1`), and it stays
        exact wherever the product rule is.

        History: `doc/shooting_history.md`, `PAC._coloured_covariance`."""
        nus = np.asarray(nus, dtype=float)
        if richardson and nus.size >= 3 and (nus.size - 1) % 2 == 0:
            q1 = cls._power_law_weights(nus, ef, richardson=False)
            q2 = np.zeros(nus.size)
            q2[::2] = cls._power_law_weights(nus[::2], ef, richardson=False)
            return (4.0 * q1 - q2) / 3.0
        wl, wr = cls._interval_weights(nus[:-1], nus[1:], ef)
        q = np.zeros(nus.size)
        q[:-1] += wl
        q[1:] += wr
        return q

    @staticmethod
    def _interval_weights(lo, hi, ef):
        """The product rule on intervals ``[lo, hi]`` (arrays): the weights
        of the left and right values -- see `_power_law_weights`."""
        lo = np.asarray(lo, dtype=float)
        h = np.log(np.asarray(hi, dtype=float) / lo)
        lam = 1.0 - float(ef)
        z = lam * h
        small = np.abs(z) < 1e-3
        zs = np.where(small, 1.0, z)
        phi1 = np.where(small, 1.0 + z / 2.0 + z ** 2 / 6.0 + z ** 3 / 24.0
                        + z ** 4 / 120.0, np.expm1(zs) / zs)
        phi2 = np.where(small, 0.5 + z / 3.0 + z ** 2 / 8.0 + z ** 3 / 30.0
                        + z ** 4 / 144.0, (np.exp(zs) * (zs - 1.0) + 1.0) / zs ** 2)
        base = lo ** lam * h
        return base * (phi1 - phi2), base * phi2

    @staticmethod
    def _loglog_integral(S, fs):
        """``int S df`` over the grid `fs` (last axis of `S`), `S` a power
        law between neighbouring points: on `[f1, f2]`, ``S = S1
        (f/f1)^p`` through both ends, integrated exactly --
        ``S1 f1 ln(r) (e^z - 1)/z``, ``z = ln(S2 f2 / (S1 f1))``, stable at
        ``z -> 0`` (1/f).  An interval with an end that is not positive
        (no power law passes through it) takes the linear trapezoid.

        History: `doc/shooting_history.md`, `PAC.sampled_variance`."""
        S = np.asarray(S, dtype=float)
        f1, f2 = fs[:-1], fs[1:]
        S1, S2 = S[..., :-1], S[..., 1:]
        L = np.log(f2 / f1)
        ok = (S1 > 0.0) & (S2 > 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            z = np.log(np.where(ok, S2, 1.0) / np.where(ok, S1, 1.0)) + L
            phi = np.where(np.abs(z) < 1e-8, 1.0 + 0.5 * z,
                           np.expm1(z) / np.where(z == 0.0, 1.0, z))
        pw = S1 * f1 * L * phi
        lin = 0.5 * (S1 + S2) * (f2 - f1)
        return np.sum(np.where(ok, pw, lin), axis=-1)

    def _refuse_driven(self, pss, what):
        if not getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.%s: phase diffusion is a property of a '
                "FREE-RUNNING oscillator. A driven circuit's phase is its "
                "source's, and its noise is pnoise's problem, not this one."
                % what)

    @staticmethod
    def _period_weights(tms, nsamp, T, pss=None):
        """Periodic TRAPEZOID weights for samples at `tms[0..nsamp-1]` over
        a period `T`: `w_j = (g_j + g_{j-1}) / 2` with `g_j` the gap to the
        next sample and the last gap closing the period.

        ⚠ NOT `h = diff(times)`: that is the LEFT RECTANGLE rule, spectrally
        accurate on a uniform periodic grid but FIRST order on a NON-UNIFORM
        one (its error is (1/2) integral h'(t) y(t) dt, not zero).  On a
        uniform grid `0.5 h + 0.5 h == h` exactly.  ⚠ AND THE TRAPEZOID IS
        ITSELF A SECOND-ORDER CAP on a smoothly varying grid: with `pss`
        given and a non-uniform grid, these are the periodic cubic-spline
        weights of `periodic_spline_weights`, broken at the nodes of landed
        events; a uniform grid is unchanged.

        History: `doc/shooting_history.md`, `PAC._period_weights`."""
        tms = np.asarray(tms, dtype=float).ravel()
        n = int(nsamp)
        t = tms[:n]
        g = np.empty(n)
        g[:-1] = t[1:] - t[:-1]
        g[-1] = float(T) + t[0] - t[-1]
        if (pss is not None and n >= 4
                and float(np.max(g)) / float(np.min(g)) - 1.0 > pss.UNIFORM_GRID_TOL):
            return periodic_spline_weights(t, T, pss._event_nodes(t, T))
        w = 0.5 * g
        w[1:] += 0.5 * g[:-1]
        w[0] += 0.5 * g[-1]
        return w

    def _ppv_states(self, pss):
        """The full-width orbit states the PPV samples and the Floquet modes
        live on, one per sample: the monodromy twin's orbit where a twin
        serves them (trap, euler), else the solve's own."""
        tw = pss.monodromy_twin() if hasattr(pss, 'monodromy_twin') else pss
        return orbit_states(tw)

    def _modulated_present(self, pss):
        """Whether a noise source of the circuit follows the orbit -- the
        bias dependence `_cy_reduced` refuses, asked without refusing."""
        try:
            self._cy_reduced(pss, 2.0 * np.pi / float(pss.period))
        except NotImplementedError:
            return True
        return False
