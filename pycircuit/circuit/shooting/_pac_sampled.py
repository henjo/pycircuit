"""Sampled noise: the noise at sampling instants, its variance and the jitter
metrics.
"""
import numpy as np
import warnings
from ._noise_components import cached_root, psd_sqrt
from ._numerics import _output_weights
from .events import EventColumns


class _SampledNoise(object):
    """Sampled noise: the noise at sampling instants, its variance and the
    jitter metrics.  A theme of `PAC` (see `pac.py`)."""

    def sampled_noise(self, pss, output, times, freqs, maxsidebands=None,
                      tail=False):
        """The one-sided PSD of the SAMPLE SERIES `y(t0 + kT)` -- DRIVEN
        circuits, white AND coloured sources.

        Returns `S` of shape `(len(times), len(freqs))` in `output`'s units
        squared per Hz, for `0 < f <= f0/2`.  `sum` over the band of `S` is
        the variance at `t0` that a sampler sees; see `sampled_variance`.
        The instants actually used (the nearest period-grid points) are left
        in `self.sampled_instants`.  ⚠ When comparing at a round instant,
        pass the grid's own times (`pss.factored_period().times`): the grid
        need not have the step count the `timestep` suggests (T/1000 can give
        999 steps), and with a 1/f source the held variance moves measurably
        between neighbouring grid points.

        ⚠ WHAT THIS IS, AND WHAT IT IS NOT.  The noise at a sampling instant
        folds every source band `nu = f + n f0` onto the series frequency
        `f`: `S(f; t0) = sum_n G_n(f; t0) CY(|nu_n|) G_n^H` with `G_n` the
        response AT `t0` to a source at `nu_n`.  It is NOT the time-averaged
        output PSD (`pnoise`); only the INTEGRATED statistics connect -- the
        cycle average of the variance at `t0` is the integral of the
        time-averaged PSD.  ⚠ A band `[fmin, fmax]` on `f` excludes
        `+-fmin` around EVERY clock harmonic, not only DC.

        ⚠⚠ ONE TRANSPOSED SOLVE PER `(t0, f)` COVERS EVERY SIDEBAND.  For a
        driven circuit the periodic adjoint at `nu = f + n f0` differs from
        the one at `f` only through phases (`exp(-j nu T)` is the same), so
        the solve seeded at `t0` is shared and each `n` is a weighted sum of
        the per-step sensitivities the reverse pass already collects.  The
        sideband count costs nothing; the grid's Nyquist (`N/2 - 1`) bounds
        it, as in `adjoint_sideband_row`.

        ⚠ SOURCE MODEL, A NAMED CHOICE: MODULATED-STATIONARY, per band --
        `sqrt(CY(x(t), nu))`, the convention of `pnoise(cyclostationary=True)`
        -- with ONE square root per independent component (per leaf element,
        white part and power-law part separately; `NoiseComponents.model`).
        (A joint root of the summed `CY` makes independent sources
        non-additive.)  ⚠ Mahmutoglu & Demir (TCAS-I 62(4), 2015) show
        that a SWITCHED MOSFET's trap (1/f) noise is whitened below the
        switching frequency and that a modulated-stationary 1/f model
        over-predicts it; the physical fix needs trap states in the device
        model, which is outside this analysis.  An agreement with another
        tool on this quantity is agreement on the convention.  A flicker
        spectrum is also singular at DC: keep `f` away from 0 (the band
        integral takes an explicit `fmin`).

        Gated against `covariance` (white, held and tracking, at the same
        instant), `adjoint_transfer_row` (the seeded row at `t0 = 0`) and,
        on an LTI circuit, the fold of `pnoise` over the same sidebands
        (white and 1/f).

        ⚠ STAGE METHODS (radau, trbdf2) RUN NATIVELY: the source enters
        every stage, so the sensitivities are read at the stage abscissae
        and `CY` at the STAGE states (one re-traversal of the orbit,
        cached).  GLM period maps are refused.

        ⚠ TIME AVERAGE, PER FREQUENCY: the mean over `t0` of this PSD is
        the fold of the time-averaged PSD, `sum_k pnoise(|f + k f0|)`, over
        EVERY output band to the grid's Nyquist (a fold cut at |k| <= 10
        leaves 0.4 %).

        ⚠ TWO THINGS LIMIT A HELD VARIANCE, AND THEY ARE NOT THE SAME
        THING.  The fold stops at `|n| <= maxsidebands`, so the source
        spectrum beyond `F = (L + 1/2) f0` is dropped: for a Lorentzian that
        is the tail `1 - (2/pi) atan(F/fc)` ~ `(2/pi) fc/F`, first order in
        the sideband count.  `tail=True` adds, per instant and series
        frequency, the `1/nu^2` extrapolation of the two OUTERMOST covered
        sidebands over the uncovered half-lines -- `A = nu_L^2 dens_L`,
        `A / (f0 (F + f0/2))` each side -- which is exact for a single pole
        and correct to `(fc/F)^2` in general; nothing is fitted.  AND the
        covered sidebands are computed on the grid: at `omega h ~ 3` per
        step (the top bands of a fold to the grid's Nyquist) gear's
        discrete transfer is far from the continuous one.  A run whose top
        sideband sits above `omega h = 1` is WARNED
        (`SAMPLED_RESOLUTION_WARN`): the remedy is more points; `tail=True`
        then closes what lies beyond the covered edge, which needs that
        edge well above the spectrum's corner (inside the corner the
        `1/nu^2` form is wrong).

        History: `doc/shooting_history.md`, `PAC.sampled_noise`.
        """
        return self._sampled_series(pss, output, times, freqs, maxsidebands,
                                    tail=tail)

    def sampled_variance(self, pss, output, times, fmin, fmax,
                         points_per_decade=40, maxsidebands=None, tail=False):
        """The variance at sampling instants over the SERIES band
        `[fmin, fmax]`, `0 < fmin < fmax <= f0/2`: the integral of
        `sampled_noise` on a log grid of `points_per_decade`, nothing added
        below `fmin`.  Returns an array over `times`.

        ⚠ POWER LAW BETWEEN THE POINTS (`_loglog_integral`): the density is
        interpolated linearly in log-log and each interval integrated
        exactly, so a white band and any pure power law are exact and the
        error is second order where the density BENDS (+6.6e-5 on a 1/f
        switched sampler at 40 per decade).

        ⚠ `fmin` AND `fmax` ARE REQUIRED.  With a 1/f source the integral
        grows as `ln(fmax/fmin)` and has no limit at `fmin -> 0`; with white
        sources only, the band removes `fmin/(f0/2)` of the full variance
        because the series PSD is flat.  Nothing is extrapolated into
        `[0, fmin]`.

        History: `doc/shooting_history.md`, `PAC.sampled_variance`.
        """
        f0 = 1.0 / float(pss.factored_period().T)
        fmin, fmax = float(fmin), float(fmax)
        if not (0.0 < fmin < fmax <= 0.5 * f0 * (1.0 + 1e-12)):
            raise ValueError(
                'PAC.sampled_variance: need 0 < fmin < fmax <= f0/2 = %.6g Hz '
                '(the sample series\' Nyquist); got fmin = %.6g, fmax = %.6g. '
                'fmin has no default: a 1/f source makes the variance grow '
                'as ln(fmax/fmin) without limit.' % (0.5 * f0, fmin, fmax))
        nf = max(2, int(np.ceil(points_per_decade * np.log10(fmax / fmin))) + 1)
        fs = np.logspace(np.log10(fmin), np.log10(fmax), nf)
        S = self._sampled_series(pss, output, times, fs, maxsidebands,
                                 tail=tail)
        return self._loglog_integral(np.asarray(S, dtype=float), fs)

    def jitter_metrics(self, pss, output, time, fmin, fmax, kmax=4,
                       maxsidebands=None, nfreq=601, dc_rectangle=False,
                       points_per_decade=40):
        """Edge jitter at ONE instant: `sigma_t`, the across-period
        correlation `rho_k`, and the three metrics that are functions of it.
        DRIVEN circuits (an oscillator is refused by `_sampled_series`; its
        accumulating share is `diffusion_constant`'s `c`).

        `sampled_noise` returns the one-sided PSD of the SAMPLE SERIES
        `y(t0 + kT)`, so that series' own autocovariance is its cosine
        transform -- no new machinery, and none of Demir 1996:

            R_k   = int_fmin^fmax S(f; t0) cos(2 pi f k T) df
            rho_k = R_k / R_0

        ⚠ TWO GRIDS.  ``R_k = int S df + int S (cos - 1) df``:
        the first takes each grid where it is the finer -- a power law
        between points (`_loglog_integral`) on a log grid of
        `points_per_decade` joined with the linear one below their
        crossover (the 1/f low end), the plain trapezoid on the linear grid
        above it (the top of the band); the second vanishes as f^2 at low f
        and oscillates at high f, so it keeps the LINEAR grid of `nfreq`
        points.
        (The linear grid alone spans the whole 1/f low end with one
        interval.)  ``R_0 - R_k`` -- the k-cycle and cycle-to-cycle metrics
        -- is the second integral alone.

        A crossing is displaced by `delta_y(t0)/slew`, so with `s` the slope
        at `t0` the three metrics A8 names follow directly:

            absolute / edge   sigma_t  = sqrt(R_0)/|s|
            k-cycle           sigma_k  = sqrt(2 (R_0 - R_k))/|s|
            cycle-to-cycle    sigma_cc = sqrt(6 R_0 - 8 R_1 + 2 R_2)/|s|

        the last from coefficients `[1, -2, 1]` on the second difference.
        ⚠ The check on that algebra: for an UNCORRELATED series they must
        collapse to `sqrt(2) sigma_t` and `sqrt(6) sigma_t`, and they do.

        Gated on a one-stage linear fixture built so the answer is known
        exactly (`tau = RC = T`, LTI noise path, so `rho_k = e^-k`), and
        against a Monte Carlo of noisy crossings.

        ⚠ `dc_rectangle` EXTRAPOLATES, WHICH THE REST OF THIS FAMILY REFUSES
        TO DO.  `S` is known only on `[fmin, fmax]`; adding `S(fmin)*fmin`
        assumes the series PSD is FLAT below `fmin`.  That is exact for white
        sources and WRONG for `1/f`, where the integral has no limit as
        `fmin -> 0` (see `sampled_variance`).  Hence off by default.
        Without it, `rho_k` drifts with `fmin` as truncation should, more at
        larger `k`.  ⚠ With a coloured source, LOWER `fmin` -- do not reach
        for the rectangle.

        ⚠ `slew` IS A FINITE DIFFERENCE ON THE PSS GRID, central about the
        instant actually used, and it converges at FIRST order.  It is
        returned so a caller can check it rather than trust it; every metric
        here is inversely proportional to it.

        ⚠ THE INSTANT IS THE CALLER'S, deliberately.  This does not hunt for a
        threshold crossing: a threshold inferred from a simulated record can
        be biased by startup, which moves the crossing off the steepest point.
        Pass the instant you mean; `instant` in the result is the grid point
        actually used.

        Returns a dict: `sigma_t`, `rho` (k = 1..kmax), `k_cycle`
        (k = 1..kmax), `cycle_to_cycle`, `slew`, `R` (k = 0..kmax), `instant`.

        History: `doc/shooting_history.md`, `PAC.jitter_metrics`.
        """
        from scipy.integrate import trapezoid
        fp = pss.factored_period()
        T = float(fp.T)
        f0 = 1.0 / T
        fmin, fmax = float(fmin), float(fmax)
        if not (0.0 < fmin < fmax <= 0.5 * f0 * (1.0 + 1e-12)):
            raise ValueError(
                'PAC.jitter_metrics: need 0 < fmin < fmax <= f0/2 = %.6g Hz; '
                'got fmin = %.6g, fmax = %.6g. As in sampled_variance, fmin '
                'has no default -- a 1/f source makes R_0 grow as '
                'ln(fmax/fmin) without limit.' % (0.5 * f0, fmin, fmax))
        kmax = int(kmax)
        if kmax < 2:
            raise ValueError(
                'PAC.jitter_metrics: kmax >= 2, because cycle-to-cycle jitter '
                'is a SECOND difference and needs R_2; got %d.' % kmax)

        fs = np.linspace(fmin, fmax, int(nfreq))
        S = np.asarray(self._sampled_series(pss, output, [time], fs,
                                            maxsidebands), dtype=float)[0]
        nlog = max(2, int(np.ceil(points_per_decade
                                  * np.log10(fmax / fmin))) + 1)
        fg = np.logspace(np.log10(fmin), np.log10(fmax), nlog)
        Sg = np.asarray(self._sampled_series(pss, output, [time], fg,
                                             maxsidebands), dtype=float)[0]
        rect = float(S[0]) * fmin if dc_rectangle else 0.0
        ## `int S df`, each grid where it is the finer: below the crossover
        ## `f*` (log spacing = linear spacing) a power law between the points
        ## of BOTH grids, which resolves a 1/f low end; above it the plain
        ## trapezoid on the linear grid, spectrally accurate on a smooth
        ## spectrum whose slope vanishes at f0/2.  (The log grid alone moved
        ## R_0 by +3.5e-5 on the white A11 fixture -- its top intervals are
        ## 6 % wide -- and power-law interpolation on the union by -1.9e-5,
        ## each 0.1 .. 0.2 % on rho_4.)
        fstar = (fs[1] - fs[0]) / (fg[1] / fg[0] - 1.0)
        js = min(int(np.searchsorted(fs, fstar)), len(fs) - 1)
        low = fg < fs[js]
        fu = np.concatenate((fs[:js + 1], fg[low]))
        order = np.argsort(fu, kind='stable')
        fu, Su = fu[order], np.concatenate((S[:js + 1], Sg[low]))[order]
        keep = np.concatenate(([True], np.diff(fu) > 0.0))
        R0 = (float(self._loglog_integral(Su[keep], fu[keep]))
              + float(trapezoid(S[js:], fs[js:])) + rect)
        R = np.array([R0 + float(trapezoid(
            S * (np.cos(2.0 * np.pi * fs * k * T) - 1.0), fs))
            for k in range(kmax + 1)])
        if not R[0] > 0.0:
            raise ValueError(
                'PAC.jitter_metrics: R_0 = %.6g is not positive, so there is '
                'no jitter to report -- are any sources noisy?' % R[0])

        t0 = float(np.asarray(self.sampled_instants, dtype=float).ravel()[0])
        row = np.asarray(self._output_waveform_row(pss, output), dtype=float)
        times = np.asarray(fp.times, dtype=float)[:len(row)]
        j = int(np.argmin(np.abs(times - t0)))
        jm, jp = max(j - 1, 0), min(j + 1, len(row) - 1)
        slew = float((row[jp] - row[jm]) / (times[jp] - times[jm]))
        ## ⚠ A slope threshold RELATIVE to the steepest slope only describes
        ## the grid (a grid point never lands exactly on a peak) and moves
        ## with N.  What does not move with the grid is the LINEARISATION
        ## this whole family rests on: a crossing is displaced by
        ## delta_y/slew only while that displacement stays small against the
        ## period.
        if not abs(slew) > 0.0:
            raise ValueError(
                'PAC.jitter_metrics: the slope at t = %.6g is exactly zero, '
                'so delta_y/slew is undefined. Pass an instant on an edge.'
                % t0)
        sigma_t = float(np.sqrt(R[0]) / abs(slew))
        if not sigma_t < 0.5 * T:
            steepest = float(np.max(np.abs(np.gradient(row, times))))
            raise ValueError(
                'PAC.jitter_metrics: at t = %.6g the implied displacement is '
                'sigma_t = %.3g s, %.3g of the period -- the first-order '
                'picture (a crossing moved by delta_y/slew) does not hold '
                'there, so every metric here would be meaningless. The slope '
                'is %.3g, %.2e of the waveform\'s steepest. Pass an instant '
                'on an edge, or check the source levels.'
                % (t0, sigma_t, sigma_t / T, abs(slew),
                   abs(slew) / max(steepest, 1e-300)))

        return {
            'sigma_t': sigma_t,
            'rho': R[1:] / R[0],
            'k_cycle': np.sqrt(np.maximum(2.0 * (R[0] - R[1:]), 0.0)) / abs(slew),
            'cycle_to_cycle': float(
                np.sqrt(max(6.0 * R[0] - 8.0 * R[1] + 2.0 * R[2], 0.0)) / abs(slew)),
            'slew': slew,
            'R': R,
            'instant': t0,
        }

    def _sampled_series(self, pss, output, times, freqs, maxsidebands,
                        tail=False):
        import scipy.sparse.linalg as spla
        self._check_circuit(pss)
        if getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.sampled_noise: an OSCILLATOR has no sampling instant '
                'fixed to its own phase -- its phase diffuses, so there is no '
                'cyclostationary sample series (see covariance()). Use '
                'oscillator_spectrum or modal_spectrum.')
        fp = pss._state_map()
        ## a Nordsieck GLM's map on the state passes as a stage map: its
        ## sources enter every stage, the output rows and the startup that
        ## opens a step (`_GLMStateStep`), so the pass collects per
        ## injection point, as a stage method's does
        stage = fp.is_stage or fp.is_glm
        T = float(fp.T)
        f0 = 1.0 / T
        tms = np.asarray(fp.times, dtype=float)
        N = len(fp.steps)
        n = fp.width
        m = pss.cir.n - 1
        fr = np.atleast_1d(np.asarray(freqs, dtype=float))
        if np.any(fr <= 0.0) or np.any(fr > 0.5 * f0 * (1.0 + 1e-12)):
            raise ValueError(
                'PAC.sampled_noise: series frequencies must lie in (0, f0/2] '
                '= (0, %.6g Hz]; the sample series has nothing above its '
                'Nyquist and a 1/f source is singular at 0.' % (0.5 * f0))
        lmax = N // 2 - 1
        L = lmax if maxsidebands is None else int(maxsidebands)
        if L > lmax:
            raise ValueError(
                "PAC.sampled_noise: maxsidebands = %d is above the grid's "
                'Nyquist (%d at %d points per period). Nothing aliases down '
                'from above what the grid represents -- use a finer period '
                'grid.' % (L, lmax, N))
        d = _output_weights(output, m)
        ts = np.atleast_1d(np.asarray(times, dtype=float))
        grid = tms[:N]
        k0s = []
        for t in ts:
            tt = float(t) % T
            dist = np.abs(grid - tt)
            dist = np.minimum(dist, T - dist)
            k0s.append(int(np.argmin(dist)))
        self.sampled_instants = grid[k0s]

        ## History: `doc/shooting_history.md`, `PAC._sampled_series`.
        ## components, rolled to the INJECTION index: step j's source enters
        ## at t_{j+1} (the reverse pass's own pairing)
        ## ⚠ WHERE THE SOURCE ENTERS: an LMM step's source enters at
        ## `t_{j+1}`, so `CY` is sampled at `x(t_{j+1})`; a stage method's
        ## enters at every stage abscissa `t_j + c_k h`, so at the STAGE
        ## states (the end-of-step state for every stage is not the limit
        ## under refinement).
        if stage:
            tinj = self._stage_times(pss, fp)
            states = self._stage_states(pss, fp)
        else:
            tinj = tms[1:N + 1]
            xs = np.asarray(pss.waveform[1], dtype=float)
            states = [xs[:, (j + 1) % N] for j in range(N)]
        nc = self._noise_components(pss, states)
        model = nc.model(float(np.min(fr)), f0)
        white, scaled, perband = [], [], []
        if model is None:
            ## the elements do not sum to the circuit's CY (warned) and the
            ## whole is not thermal-plus-power-law: one root of the whole
            ## circuit's CY per band
            perband.append(cached_root(nc.cy_at_states))
        else:
            white = [psd_sqrt(A) for _key, A in model.white_parts]
            ## the coloured components as the modal spectra and the folds
            ## take them (`colour_groups`): a fixed column set with its power
            ## weight, or a root per band frequency
            groups = nc.colour_groups(
                model, 2.0 * np.pi * float(np.min(fr)), f0, L,
                'sampled_noise', warn_touch=False)
            scaled = [(G, s) for kind, G, s in groups if kind == 'fixed']
            perband = [G for kind, G, _s in groups if kind == 'band']

        tol = max(self.KRYLOV_FACTOR * pss.par.reltol, 1e-14)
        ns = np.arange(-L, L + 1)
        ## ⚠ THE TOP SIDEBANDS ARE COMPUTED ON THE GRID.  At `omega h > 1`
        ## per step their discrete transfer is not the continuous one, and
        ## the held variance comes out short by a factor that looks like a
        ## kernel constant.
        _hmax = float(np.max(np.diff(tms))) if len(tms) > 1 else T / max(N, 1)
        _wh = 2.0 * np.pi * (L + 0.5) * f0 * _hmax
        if _wh > pss.SAMPLED_RESOLUTION_WARN and not getattr(self, '_sampled_res_warned', False):
            self._sampled_res_warned = True
            warnings.warn(
                'PAC.sampled_noise: the top sideband (|n| = %d, %.3g Hz) sits '
                'at omega h = %.2f per step on this grid; a two-step method\'s '
                'discrete transfer there is far from the continuous one, and '
                'the held variance came out low by 3x the pure tail at '
                'omega h = 3 (1.05x at 0.2) on a switched capacitor. Use '
                'more points per period; tail=True closes the spectrum '
                'beyond the covered edge, but only once that edge is well '
                'above the spectrum\'s corner (measured: 0.71 x kT/C with the '
                'edge inside the corner, 0.998 with it 12x beyond).'
                % (L, (L + 0.5) * f0, _wh), RuntimeWarning, stacklevel=3)
        S = np.zeros((len(ts), len(fr)))
        ## ⚠ ON A STAGED SOLVE THE SAMPLE'S ADJOINT IS BORDERED -- the dual
        ## of the bordered forward solve, as `adjoint_sideband_row`'s: the
        ## operator is the TOTAL map's transpose, the sample is read at
        ## FIXED time (its costate carries `dtheta/dx_0^T Pk_fixed[k0]^T
        ## d`), and the source's own motion of the crossings enters as a
        ## third reverse pass carrying `-zeta_k W_k` at the event nodes,
        ## `zeta = Gt^-T (g_theta + a P_theta^T z)`.  Unbordered, the jitter
        ## sampler's held node has no path to the threshold's noise at all.
        _ev = EventColumns.of(pss, n)
        if _ev is not None:
            _dth = np.asarray(_ev.dth, dtype=float)
            _MtT = _ev.total_matvec(fp.matvec_transposed, transposed=True)
            _Pkf, _t_, _x_ = self._fixed_time_event_columns(pss)
        for ti, k0 in enumerate(k0s):
            for fi, f in enumerate(fr):
                alpha = np.exp(-2j * np.pi * f * T)
                inject = np.zeros((N, m), dtype=complex)
                inject[k0] = np.exp(-2j * np.pi * f * tms[k0]) * d
                if _ev is None:
                    A_ = spla.LinearOperator(
                        (n, n), dtype=complex,
                        matvec=lambda v, a=alpha: np.asarray(v) - a * fp.matvec_transposed(v))
                else:
                    A_ = spla.LinearOperator(
                        (n, n), dtype=complex,
                        matvec=lambda v, a=alpha: np.asarray(v) - a * _MtT(v))
                if stage:
                    seedv = np.exp(-2j * np.pi * f * tms[k0]) * d
                    g, cA = self._stage_pass(pss, fp, np.zeros(m), (k0, seedv))
                    if _ev is not None:
                        g_theta = np.exp(-2j * np.pi * f * tms[k0]) * (_Pkf[k0].T @ d)
                        g = np.asarray(g) + _dth.T @ g_theta
                    z = self._gmres_checked(A_, g, tol, 'the sampled adjoint solve')
                    _l, cZ = self._stage_pass(pss, fp, z)
                    Sv = -(cA + alpha * cZ)                          # N s x m
                    if _ev is not None:
                        zeta = _ev.collapsed_zeta(g_theta, alpha, z)
                        _l2, cE = self._stage_pass(pss, fp, np.zeros(m), _ev.injection_dict(zeta))
                        Sv = Sv - cE
                else:
                    g, t_inj, _st = fp.matvec_transposed(
                        np.zeros(n, dtype=complex), collect=True, inject=inject)
                    if _ev is not None:
                        g_theta = np.exp(-2j * np.pi * f * tms[k0]) * (_Pkf[k0].T @ d)
                        g = np.asarray(g) + _dth.T @ g_theta
                    z = self._gmres_checked(A_, g, tol, 'the sampled adjoint solve')
                    _e, t_z, _st = fp.matvec_transposed(z, collect=True)
                    Sv = -(np.asarray(t_inj) + alpha * np.asarray(t_z))  # N x m
                    if _ev is not None:
                        zeta = _ev.collapsed_zeta(g_theta, alpha, z)
                        _g2, t_ev, _st2 = fp.matvec_transposed(
                            np.zeros(n, dtype=complex), collect=True,
                            inject=_ev.injection(zeta, N, m))
                        Sv = Sv - np.asarray(t_ev)
                nu = f + ns * f0
                E = (np.exp(2j * np.pi * nu[:, None] * tinj[None, :])
                     * np.exp(-2j * np.pi * ns * f0 * tms[k0])[:, None])
                dens = 0.0
                ## per-sideband densities, kept for the tail closure
                _pb = np.zeros(len(ns))
                for SA in white:
                    R = E @ np.einsum('ji,jik->jk', Sv, SA)
                    _pb += np.sum(np.abs(R) ** 2, axis=1)
                for SB, weight in scaled:
                    R = E @ np.einsum('ji,jik->jk', Sv, SB)
                    c = weight(2.0 * np.pi * np.abs(nu))
                    _pb += c * np.sum(np.abs(R) ** 2, axis=1)
                ## (each per-band component hands back its ROOT, cached per
                ## frequency: `perband_root_sampler`)
                for comp in perband:
                    for bi, nb in enumerate(nu):
                        R = E[bi] @ np.einsum(
                            'ji,jik->jk', Sv, comp(2.0 * np.pi * abs(nb)))
                        _pb[bi] += float(np.sum(np.abs(R) ** 2))
                dens = float(np.sum(_pb))
                if tail and len(ns) >= 3:
                    ## the 1/nu^2 extrapolation of each outermost covered
                    ## sideband over its uncovered half-line: int_F^inf A/nu^2
                    ## per f0 of series bandwidth = A / (f0 F), F the edge
                    ## half a band beyond the last centre.  Coefficient 1,
                    ## derived; exact for one pole.
                    for _edge in (0, -1):
                        _nu_e = float(nu[_edge])
                        _A = _nu_e ** 2 * float(_pb[_edge])
                        _F = abs(_nu_e) + 0.5 * f0
                        dens += _A / (f0 * _F)
                S[ti, fi] = dens
        return S
