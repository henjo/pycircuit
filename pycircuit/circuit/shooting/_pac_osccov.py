"""An oscillator's covariance: the transverse responses, the orbital lines, the
edge jitter and the mode weights.
"""
import numpy as np
from ._numerics import edge_slope, output_index
import warnings
from .events import EventColumns


class _OscillatorCovariance(object):
    """An oscillator's covariance: the transverse responses, the orbital
    lines, the edge jitter and the mode weights.  A theme of `PAC` (see
    `pac.py`)."""

    def _node_projectors(self, pss):
        """``Pi_j = I - xdot_j s_j^T / (s_j . xdot_j)`` at every node of an
        oscillator's orbit, `(N + 1, m, m)`: the OBLIQUE projection onto the
        transverse (orbital) deviation, Demir's `v_1^T y = 0` --
        `xdot_j` the DAE's own rate at the node (`_orbit_rate`), `s_j` the
        PPV in state coordinates (`ppv()['samples']`, `C^T v_1`).  The
        normalisation ``s_j . xdot_j`` is 1 on the exact orbit; it is
        divided by, not assumed."""
        m = pss.cir.n - 1
        _v, info = pss.ppv()
        S = np.asarray(info['samples'], dtype=float)[:, :m]
        xd = np.asarray(self._orbit_rate(pss, []), dtype=float)[:, :m]
        ## (the PPV samples cover nodes 0..N-1; node N is node 0 a period on)
        I = np.eye(m)
        return np.array([I - np.outer(xd[j], S[j % S.shape[0]])
                         / float(S[j % S.shape[0]] @ xd[j])
                         for j in range(xd.shape[0])])

    def _orbital_lines(self, pss, fmin, fmax):
        """``(centre, half-width)`` in Hz of every orbital line an
        oscillator's transverse response has in ``[fmin, fmax]``: each
        non-phase Floquet mode ``mu`` at ``k f0 +- Im(mu)/2 pi``, every `k`,
        half-width ``|Re mu| / 2 pi`` -- for `_coloured_covariance` to
        resolve before it adapts (a van der Pol at Q = 8 has lines 0.02 f0
        wide at EVERY harmonic)."""
        modes = pss.floquet_modes()
        _kph, orb = self._phase_mode_split(pss, modes,
                                           'PAC.oscillator_covariance')
        f0 = 1.0 / float(pss.period)
        out = []
        for l in orb:
            mu = complex(modes[l]['mu'])
            hw = abs(mu.real) / (2.0 * np.pi)
            fi = mu.imag / (2.0 * np.pi)
            ## ⚠ ONLY THE HARMONICS THE MODE COUPLES: an input at `nu` drives
            ## the mode through harmonic `k` of its adjoint `q_l`, and a
            ## harmonic with no energy there carries no line.  Seeding every
            ## harmonic to Nyquist (200 at 400 points) cost thousands of
            ## solves on lines of no weight.
            ## ⚠ AND ONLY ISOLATED LINES: half-width below f0/4.  A strongly
            ## damped mode (a relaxation oscillator's: multipliers 0.02, 0,
            ## lines 0.62 and 3.4 f0 wide) OVERLAPS its neighbours into a
            ## smooth response, which the adaptive rule handles; seeding its
            ## 140 lines cost more than 20 minutes.
            if hw >= 0.25 * f0:
                continue
            Vl = self._period_dft(pss, np.asarray(modes[l]['q'])[:, :-1].T).T
            Nh = Vl.shape[1]
            en = np.sum(np.abs(Vl) ** 2, axis=0)
            emax = float(en.max()) if en.size else 0.0
            for k in range(Nh // 2 + 1):
                if emax <= 0.0 or max(en[k], en[-k % Nh]) < 1e-10 * emax:
                    continue
                for sgn in ((1.0, -1.0) if fi != 0.0 else (1.0,)):
                    c0 = k * f0 + sgn * fi
                    if float(fmin) < c0 < float(fmax) and hw > 0.0:
                        out.append((c0, hw))
        return out

    def _transverse_responses(self, pss, fp, freqs, u_ac, u_points=None):
        """`_forced_responses` for an OSCILLATOR's transverse deviation: per
        frequency the steady response at every node projected with
        `_node_projectors`, and no crossings (`[None]`).

        ⚠ REPLAYED FROM THE BOUNDED PART.  The deflated solve's answer is
        ``y = w + s u / (1 - alpha)``; replaying `y` and projecting would
        cancel the `1/offset` term against itself and leave its
        discretisation error times `1/offset` -- which, against a 1/f
        density, diverges as `fmin` falls.  Replaying `w` never forms it:
        the pole's part propagates along the tangent, which the projection
        removes.

        ⚠ ON A STAGED OSCILLATOR the map is the TOTAL one
        (`EventColumns.total_matrix`), the source moves the crossings itself
        (``dtheta_f``, `forced_shift`: its ``P_theta dtheta_f`` enters the
        right-hand side), and the node responses are at FIXED times: the
        crossings' motion ``dtheta = dth w + dtheta_f`` enters through
        `_fixed_time_event_columns`, as `_forced_responses` does it.  The
        pole's part, dropped with `y`, moves the crossings along the orbit
        with it, and its fixed-time response is again along ``xdot_j``.

        History: `doc/shooting_history.md`, `PAC._transverse_responses`."""
        T = float(fp.T)
        m = self.cir.n - 1
        N = len(fp.steps)
        ## ⚠ ONCE PER MAP, NOT PER FREQUENCY: the projectors, the border pair
        ## and -- below `FLOQUET_DENSE_LIMIT` -- the monodromy itself, which
        ## does not depend on the frequency.  Per frequency the bordered
        ## system is then a small dense solve (its `w` exact: no GMRES, no
        ## refinement).  Recomputing `ppv()` inside every deflated solve was
        ## 65 of 145 s on a 100-point van der Pol.  Cleared by the caller.
        cache = getattr(self, '_transverse_cache', None)
        if cache is None or cache[0] is not fp:
            nw = fp.width
            _v, info = pss.ppv()
            vb = np.asarray(_v, dtype=float).ravel()
            ub = np.asarray(info['tangent_pair'], dtype=float).ravel()
            Md = (np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                                   for e in np.eye(nw)])
                  if nw <= pss.FLOQUET_DENSE_LIMIT else None)
            staged = None
            _evd = EventColumns.of(pss, nw)
            if _evd is not None:
                ## the total map; the crossings' state sensitivity; the event
                ## columns at fixed time
                if Md is not None:
                    Md = np.asarray(_evd.total_matrix(Md), dtype=float)
                staged = (_evd, np.asarray(_evd['P_end'], dtype=complex),
                          np.asarray(pss._event_columns.dth, dtype=float),
                          self._fixed_time_event_columns(pss)[0])
            cache = (fp, self._node_projectors(pss), vb, ub, Md, staged)
            self._transverse_cache = cache
        _fp, Pi, vb, ub, Md, staged = cache
        tol = max(pss.par.reltol * self.KRYLOV_FACTOR, 1e-14)
        out = []
        for f in freqs:
            w, _ = pss._forced_replay(fp, f, u_ac, u_points=u_points)
            a = np.exp(-2j * np.pi * f * T)
            rhs = a * np.asarray(w, dtype=complex)
            dth_f = None
            if staged is not None:
                ## the source's own motion of the crossings
                _evd, Pth, _dthx, _Pkf = staged
                _e0, f_steps = pss._forced_replay(
                    fp, f, u_ac, y0=np.zeros(fp.width, dtype=complex),
                    collect=True, u_points=u_points)
                f_nodes = [np.zeros(m, dtype=complex)] + [
                    np.asarray(v_, dtype=complex)[:m] for v_ in f_steps]
                dth_f = _evd.forced_shift(f_nodes)
                rhs = rhs + a * (Pth @ dth_f)
            if Md is not None:
                nw = Md.shape[0]
                B = np.zeros((nw + 1, nw + 1), dtype=complex)
                B[:nw, :nw] = np.eye(nw) - a * Md
                B[:nw, nw] = ub
                B[nw, :nw] = vb
                wb = np.linalg.solve(B, np.concatenate((rhs, [0.0])))[:nw]
            else:
                _y, wb = self._deflated_solve(pss, a, rhs, tol=tol, parts=True)
            _e, ysteps = pss._forced_replay(fp, f, u_ac, y0=wb, collect=True,
                                            u_points=u_points)
            Y = np.array([np.asarray(wb)[:m]]
                         + [np.asarray(v)[:m] for v in ysteps])[:N + 1]
            if staged is not None:
                ## the crossings' motion from the bounded part, at fixed time
                _evd, Pth, _dthx, _Pkf = staged
                dth = _dthx @ np.asarray(wb)[:_dthx.shape[1]] + dth_f
                Y = Y + np.tensordot(_Pkf[:len(Y)], dth, axes=(2, 0))
            out.append(np.einsum('jab,jb->ja', Pi[:len(Y)], Y))
        return out, [None] * len(freqs)

    ## ⚠ ON AN ASYMMETRIC ORBIT this route is validated by a direct SDE
    ## simulation of the variational system (the calibrated `Var(i) =
    ## CY/(2h)` injection, phase projected out every step), which shares no
    ## Lyapunov solve and no modal sum: 0.02 % at `a = 0.30` on van der Pol
    ## + `a u^2`.  ⚠ Asymmetry changes the MODE SHAPES, so the noise
    ## projected onto the orbital direction grows and the transverse
    ## variance RISES despite the faster relaxation (`|lam2|` falls) -- the
    ## physical argument that it should shrink is wrong.
    ## History: `doc/shooting_history.md`, `PAC.oscillator_covariance`.
    def oscillator_covariance(self, pss, samples=False, colour_fmin=None,
                              colour_fmax=None, points_per_decade=40,
                              pair=False):
        """The state covariance of a FREE-RUNNING oscillator, split in two.

        Returns `(K_orb, d, info)`.  `K_orb` is the BOUNDED periodic
        (orbital) part of the covariance at `t = 0`; `d` is the growth per
        period along the orbit tangent (see the split below).

        ⚠ "BOUNDED" IS NOT "TRANSVERSE".  `K_orb` has the SECULAR growth
        removed and still contains the phase direction's bounded
        within-period variance.  Demir's orbital deviation `y` is the
        OBLIQUE projection `v_1^T y = 0`, so the transverse covariance is
        `Pi K_orb Pi^T` with `Pi = I - u v^T/(v^T u)` -- which is what
        `orbital_correlation`'s eq (23) sum equals, and what `K_orb` itself
        does NOT equal (2-6 %, falling as 1/Q).  Read `K_orb` as the
        bounded part; project it if you want `R_yy(0)`.  See
        `orbital_correlation`.

            K(t_0 + n T) = K_orb + n d u u^T

        exactly, for every integer `n`, with `u` the pair-space tangent
        scaled so its first block is `xdot(0)`.

        ⚠ WITH `samples=True` THE SPLIT MOVES WITH THE ORBIT, AND THE
        OBVIOUS READING IS WRONG.  `info['orbital_samples'][j]` is `P(t_j)`,
        the solution started from `K_orb` at `t = 0`, and it satisfies

            K(t_j + n T) = P(t_j) + n d u_j u_j^T,   u_j = Phi(t_j, 0) u

        so `P` is periodic UP TO the growth -- `P(T) = P(0) + d u u^T`, not
        `P(T) = P(0)`.  The walk is along the orbit and the orbit turns, so
        the growth DIRECTION is the propagated tangent rather than a fixed
        `u`.

        ⚠ THIS IS THE OBJECT `covariance` REFUSES TO RETURN: there is no
        periodic solution.  `lambda_1 = 1` gives `lambda_1^2 = 1`, so
        `I - M kron M` is exactly singular, with a cleanly ONE-DIMENSIONAL
        null space spanned by `u kron u` and left null `v kron v`.  So it
        borders exactly as the PPV and the deflated PAC solve do, and the
        border is the pair the rest of this class already computes.

        ⚠ THE SPLIT IS NOT A NUMERICAL DEVICE, IT IS THE ANSWER.  Demir
        2002: an oscillator's noise is STATIONARY, not cyclostationary,
        because "noisy autonomous systems cannot provide a perfect time
        reference".  `K_orb` is the part a designer can read as an
        amplitude/orbital noise -- it settles, it is periodic, it is
        finite.  `n d u u^T` is the random walk ALONG the orbit, which
        never settles and which no periodic object can hold.

            [ I - M kron M    u kron u ] [ vec(K_orb) ]   [ vec(K_1) ]
            [ (v kron v)^T        0    ] [     d      ] = [     0    ]

        ⚠ AND `d` HAS A CLOSED FORM THAT NEEDS NO KRONECKER AT ALL.
        Left-multiplying the first row by `(v kron v)^T` kills the singular
        block, leaving

            d = (v^T K_1 v) / (v . u)^2

        an `O(n^2)` contraction against the `n^4` solve.  Both are computed
        and `info['d_residual']` is their relative difference; they are the
        same quantity by construction, so a disagreement means the border
        pair is wrong rather than that one route is less accurate.

        ⚠ `(v . u)` IS NOT 1 AND ASSUMING IT IS COSTS A FACTOR OF 2.3.
        `ppv()` normalises on the FIRST BLOCK, `v[:m] . xdot = 1` -- the
        normalisation an injected current sees, and what every other
        shipped path does.  The FULL PAIR contraction is a different number
        (0.663 on van der Pol), which is why `d` is written with the pair
        inner product spelled out.

        ⚠ `d` ALONE IS MEANINGLESS WITHOUT PINNING `u`'s SCALE.  Rescaling
        `u -> s u` sends `d -> d / s^2`, so only the PRODUCT `d u u^T` --
        returned as `info['growth']` -- is an invariant of the circuit.
        `u` is pinned here by `C u = q`, the same condition `ppv()` uses to
        scale the tangent, which makes its first block exactly `xdot(0)`
        and gives `d` its physical reading below.

        ⚠ WHICH MAKES `d / T` A COMPLETELY INDEPENDENT ROUTE TO THE
        DIFFUSION CONSTANT, and that is this method's real gate.  A phase
        deviation `alpha` displaces the state by `alpha u`, so the growing
        covariance is `Var(alpha) u u^T = c t u u^T`, giving `d = c T`.
        The two computations share only the `CY/2` convention: `c` is a
        quadratic form in the ADJOINT-replayed PPV, while `d` comes from a
        FORWARD Lyapunov recursion closed by a bordered Kronecker solve.
        Their anchors are independent too (`covariance`'s injection to
        `kT/C`, `diffusion_constant` to a nonlinear Monte Carlo reading
        phase from zero crossings), so `info['c_from_growth']` against
        `diffusion_constant` closes a loop between two separately anchored
        quantities.

        ⚠ COST: the bordered solve has `(2m)^2 + 1` unknowns and is dense,
        so it is `O(m^4)` like `covariance`.  Small circuits only.  The
        closed form for `d` is cheap; pass `samples=False` and read
        `info['c_from_growth']` if the orbital part is not wanted.

        `info['K_transverse']` is the TRANSVERSE covariance at `t = 0`,
        ``Pi K_orb Pi^T`` in the node space (`_node_projectors`), and with
        `samples=True` `info['transverse_samples']` at every node.

        ⚠ A COLOURED SOURCE needs the band `colour_fmin` / `colour_fmax` /
        `points_per_decade`, as `covariance`.  Its phase does not diffuse --
        1/f FM grows faster than linearly -- so it has no bounded-plus-
        random-walk split, and it enters the TRANSVERSE covariance alone:
        the white sources go through the recursion above (`K_orb`, `d` and
        `c_from_growth` are theirs, and a warning says so), the coloured
        ones through the deflated forced responses over the band, projected
        at every node (`_transverse_responses`), as
        `info['K_coloured']` / `info['coloured_samples']`, and
        `info['K_transverse']` / `info['transverse_samples']` hold both.
        The coloured PHASE is `phase_psd`'s.  On a staged oscillator the
        responses are the total map's at fixed time (`_transverse_responses`).

        History: `doc/shooting_history.md`, `PAC.oscillator_covariance`.
        """
        self._check_circuit(pss)
        if not getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.oscillator_covariance: this splits a covariance that '
                'GROWS into a bounded part plus a random walk along the '
                'orbit. A driven circuit has neither -- its covariance '
                'settles, and I - M kron M is nonsingular. Use '
                'covariance().')
        ## the source-injection surfaces use the gear twin when the Floquet
        ## source is TR-BDF2 (its two-stage Q_j is not built).  Swapped BEFORE
        ## both the Lyapunov pieces and `ppv` below, so the bordering keeps
        ## them on one host -- see `_lyapunov_host`.
        pss = pss._lyapunov_host()
        fmin, fmax = colour_fmin, colour_fmax     # (the names the internals use)
        col = self._coloured_prepare(pss, fmin, fmax, points_per_decade,
                                     'oscillator_covariance')
        As, Qs, K1, M, m, n = self._lyapunov_pieces(
            pss, 'oscillator_covariance',
            white=None if col is None else col['white'])
        ## a staged oscillator closes on the TOTAL map with the crossings'
        ## noise-driven motion in the injection -- the same `_event_closure`
        ## as `covariance`, whose `u`, `v` below are the total map's already
        _bordered = self._event_closure(pss, As, Qs, M, m, n)
        if _bordered is not None:
            M, K1, _samples_unused, _pieces_unused = _bordered

        v, pinfo = pss.ppv()
        v = np.asarray(v, dtype=float).ravel()
        u = np.asarray(pinfo['tangent_pair'], dtype=float).ravel()
        xdot = np.asarray(pinfo['xdot'], dtype=float).ravel()
        if n == 2 * m and v.shape[0] == m:
            ## ⚠ THE TRAPEZOIDAL PLAIN PAIR `(x, iq)`: its map re-seeds `iq`
            ## at every period start, so its
            ## last `m` columns are zero and its null vectors follow from the
            ## state map's -- the left one `[v; 0]`, the right one the tangent
            ## with the `iq` block it carries, `M[:, :m] u`.  The result is as
            ## good as trap's own map: first order on a limit cycle (d -2.5 %
            ## at 200 points, -1.1 % at 800 on van der Pol, against a radau
            ## reference; the trbdf2 twin -4.8e-6).
            ## History: `doc/shooting_history.md`, `PAC.oscillator_covariance`.
            u = np.asarray(M, dtype=float)[:, :m] @ u[:m]
            v = np.concatenate((v, np.zeros(m)))
        ## rescale the bordered solve's DIRECTION onto the tangent `ppv`
        ## already scaled by `C u = q`; least squares so this is stable
        ## even where `u[:m]` is small, and exact where it is not.
        uu = float(u[:m] @ u[:m])
        if uu == 0.0:
            raise ValueError(
                'PAC.oscillator_covariance: the tangent has no first '
                'block, so its scale cannot be pinned to xdot(0).')
        u = u * (float(u[:m] @ xdot) / uu)

        vu = float(v @ u)
        if vu == 0.0:
            raise ValueError(
                'PAC.oscillator_covariance: the left and right null '
                'directions are orthogonal in the pair space, so the '
                'bordered system is singular. That should not happen on a '
                'converged limit cycle.')
        d_closed = float(v @ K1 @ v) / (vu * vu)

        S = np.eye(n * n) - np.kron(M, M)
        uk = np.kron(u, u)
        vk = np.kron(v, v)
        B = np.zeros((n * n + 1, n * n + 1))
        B[:n * n, :n * n] = S
        B[:n * n, n * n] = uk
        B[n * n, :n * n] = vk
        rhs = np.concatenate((K1.reshape(-1), [0.0]))
        z = np.linalg.solve(B, rhs)
        K_orb = z[:n * n].reshape(n, n)
        K_orb = 0.5 * (K_orb + K_orb.T)
        d = float(z[n * n])

        scale = max(abs(d), abs(d_closed), 1e-300)
        info = {'d_closed_form': d_closed,
                'd_residual': abs(d - d_closed) / scale,
                'growth': d * np.outer(u, u),
                'c_from_growth': d / float(pss.period),
                'tangent_pair': u,
                'ppv_pair': v,
                'pair_inner': vu,
                'sigma_min': float(np.linalg.svd(S, compute_uv=False)[-1]),
                'sigma_min_bordered':
                    float(np.linalg.svd(B, compute_uv=False)[-1]),
                'null_residual': float(np.linalg.norm(S @ uk))
                                 / max(float(np.linalg.norm(uk)), 1e-300),
                'ppv': pinfo}
        if samples:
            ## ⚠ THE GROWTH DIRECTION MOVES WITH THE ORBIT.  The invariant
            ## is `K(t_j + nT) = K_orb(t_j) + n d u_j u_j^T` with `u_j` the
            ## FORWARD-propagated tangent, not `u` held fixed -- the walk
            ## is along the orbit, and the orbit turns.
            orb, grw, K, uj = [K_orb], [d * np.outer(u, u)], K_orb, u
            ## a GLM's native steps are `(x, P)` wide: pad, read the `x`
            ## block (`_lyap_walk`)
            na = As[0].shape[0] if As else n
            if na != n:
                K = np.pad(K_orb, ((0, na - n), (0, na - n)))
                uj = np.pad(u, (0, na - n))
            for A, Q in zip(As, Qs):
                K = A @ K @ A.T + Q
                uj = A @ uj
                orb.append((0.5 * (K + K.T))[:n, :n])
                grw.append(d * np.outer(uj[:n], uj[:n]))
            info['orbital_samples'] = orb
            info['growth_samples'] = grw
            info['times'] = np.asarray(pss.factored_period().times,
                                       dtype=float)
        ## the TRANSVERSE covariance in the node space; the coloured sources'
        ## enter here alone (see the docstring)
        Pi = self._node_projectors(pss)
        info['K_transverse'] = Pi[0] @ K_orb[:m, :m] @ Pi[0].T
        if samples:
            info['transverse_samples'] = [
                Pi[j] @ Kj[:m, :m] @ Pi[j].T
                for j, Kj in enumerate(info['orbital_samples'][:len(Pi)])]
        if col is not None:
            try:
                Kc, _none = self._coloured_covariance(
                    pss, col, m, m, responses=self._transverse_responses,
                    lines=self._orbital_lines(pss, col['fmin'], col['fmax']),
                    all_nodes=bool(samples))
            finally:
                self._transverse_cache = None
            info['K_coloured'] = Kc[0]
            info['K_transverse'] = info['K_transverse'] + Kc[0]
            if samples:
                info['coloured_samples'] = list(Kc)
                info['transverse_samples'] = [
                    a + b for a, b in zip(info['transverse_samples'], Kc)]
            warnings.warn(
                'PAC.oscillator_covariance: this circuit has a COLOURED source. '
                'K_orb, d and c_from_growth are the WHITE sources\' alone: a '
                'coloured source\'s phase does not diffuse (1/f FM grows '
                'faster than linearly) and has no growth rate -- its phase '
                'noise is phase_psd\'s. Its transverse covariance is in '
                "info['K_coloured'], and info['K_transverse'] holds both.",
                RuntimeWarning, stacklevel=2)
        ## ⚠ `m x m` WHATEVER THE METHOD (2026-09-29), as `covariance`: a
        ## two-step method's pair space is its own, `pair=True` keeps it
        ## (`orbital_mode_weights` reads the Floquet modes there)
        if not pair:
            K_orb = K_orb[:m, :m]
            info['growth'] = info['growth'][:m, :m]
            for key in ('orbital_samples', 'growth_samples'):
                if key in info:
                    info[key] = [K[:m, :m] for K in info[key]]
        return K_orb, d, info

    def oscillator_edge_jitter(self, pss, output, time, kmax=8):
        """The ADDITIVE (non-accumulating) edge jitter of a FREE-RUNNING
        oscillator -- the number a clock designer wants at the last buffer,
        and the one `c` does not contain.

        `sampled_noise` refuses an autonomous PSS by design: a diffusing phase
        has no sampling instant fixed to it, so there is no cyclostationary
        sample series and `jitter_metrics` cannot be used here.  But the split
        `oscillator_covariance` already returns says exactly what to do --

            K(t_j + n T) = P(t_j) + n d u_j u_j^T

        `n d u_j u_j^T` is the random walk ALONG the orbit (that is `c`);
        `P(t_j)` is bounded and does not accumulate.  A crossing is displaced
        by `delta_y/slew`, so with `s` the slope at `t_j`

            sigma_t^2 = e^T Pi P(t_j) Pi^T e / s^2,  Pi = I - u_j v_j^T/(v_j^T u_j)

        and the k-lag law that a designer actually measures is

            Var(tau_{n+k} - tau_n) = c k T + 2 sigma_t^2 (1 - rho_k)

        ⚠⚠ `P` ITSELF IS THE WRONG OBJECT, and by a margin that hides easily.
        "Bounded is not transverse": `P` keeps the phase direction's bounded
        within-period variance, which the ORBITAL deviation excludes (Demir's
        `v_1^T y = 0`).  `projection_share` in the result is how much it would
        cost you here.  ⚠ IT IS A PROPERTY OF THE SOURCE MIX, NOT OF THIS
        METHOD: the phase direction is exactly what TANK noise drives, so a
        fixture whose oscillator is quiet shows the projection doing nothing
        (0.02 % with the tank 1000x quieter than the buffers) and teaches the
        wrong lesson; one with a noisy tank shows 16 %.  It also falls as
        1/Q, so a high-Q fixture hides it too.
        ⚠ `P - (t/T) G` is the SUPERSEDED prescription and is also wrong; see
        `orbital_correlation`, which gates the projection three ways.

        ⚠ `u_j` COMES BACK FROM `growth_samples[j]`, WHICH IS A RANK-ONE
        MATRIX `d u_j u_j^T`, not a vector -- its leading eigenpair gives
        `sqrt(d) u_j`, and `Pi` is invariant to that scale (and to `v`'s), so
        taking `v` from `ppv()` in a separate call is safe here.  Everything
        is sliced `[:m, :m]` out of PAIR space.

        ⚠ THE SLOPE IS A LOCAL QUADRATIC FIT AT THE REQUESTED INSTANT, NOT A
        TWO-POINT DIFFERENCE.  A threshold crossing sits at a different
        fraction of a step on every grid, so a straddling difference moves
        with the grid (it changes sign under refinement, and was 2.8 % low
        at 240 points); every quantity here goes as `1/s^2`.

        Gated against a Monte Carlo of a noisy transient with no PSS, no
        adjoint and no Lyapunov solve in it (van der Pol tank driving three
        tanh buffers): `MC / analysis = 1.0066 +/- 0.0102`, with the same
        orbit on both sides.
        ⚠ When validating this against a transient, run the Monte Carlo on the
        SAME integrator as the PSS, or divide by the slope of the orbit the
        Monte Carlo actually runs on -- otherwise the mismatch enters squared
        (an Euler orbit's slope against gear's: 5.7 % at 240 points).

        ⚠ `k_cycle_bound` IS THE LARGE-`k` FORM, with `rho_k` taken to zero,
        and named for what it is: the orbital part's across-period
        correlation is not computed here, so at small `k` the true k-cycle
        jitter is LOWER (the `(1 - rho_k)` factor is below 1).  (Until
        2026-09-29 it was `k_cycle`, the key under which `jitter_metrics`
        returns the EXACT value.)  For a DRIVEN circuit use `jitter_metrics`,
        which computes `rho_k` from the sample series.

        ⚠ THE INSTANT IS THE CALLER'S.  This does not hunt for a crossing: a
        threshold taken from a simulated record can be biased by startup, and
        that moves the instant off the steepest point.  `instant` in the
        result is the grid point used.

        Returns a dict: `sigma_t`, `A` (= sigma_t^2), `c`, `slew`,
        `k_cycle_bound` (k = 1..kmax), `instant`, `d`, `projection_share`.

        History: `doc/shooting_history.md`, `PAC.oscillator_edge_jitter`.
        """
        output = output_index(pss, output)
        import warnings as _warnings
        self._check_circuit(pss)
        ## ONE ORBIT: the covariance's host (a GLM's or trap's twin) supplies
        ## the period, the factored period and the PPV as well -- read off
        ## the run itself they come from another discretisation
        ## History: `doc/shooting_history.md`, `PAC.oscillator_edge_jitter`.
        pss = pss._lyapunov_host()
        K_orb, d, info = self.oscillator_covariance(pss, samples=True)
        m = self.cir.n - 1
        T = float(pss.period)
        fp = pss.factored_period()
        times = np.asarray(fp.times, dtype=float)
        row = np.asarray(self._output_waveform_row(pss, output), dtype=float)
        nt = int(min(len(times), len(row)))
        if nt < 5:
            raise ValueError(
                'PAC.oscillator_edge_jitter: the period grid has %d points; '
                'the slope needs at least 5.' % nt)
        times, row = times[:nt], row[:nt]

        ## ⚠ THE SLOPE IS TAKEN AT THE REQUESTED INSTANT, NOT AT THE SNAPPED
        ## GRID POINT, and the difference is the whole correction.  `time` is
        ## typically a threshold crossing, which sits at a different FRACTION
        ## of a step on every grid; differentiating at the nearest sample
        ## instead reproduces the straddling value, and every quantity here
        ## goes as 1/s^2.  (`edge_slope`, shared with `jitter_metrics`.)
        slew, j = edge_slope(times, row, float(time))
        if not abs(slew) > 0.0:
            raise ValueError(
                'PAC.oscillator_edge_jitter: the slope at t = %.6g is exactly '
                'zero, so delta_y/slew is undefined. Pass an instant on an '
                'edge.' % times[j])

        Ps = [np.asarray(P, dtype=float)[:m, :m] for P in info['orbital_samples']]
        G = [np.asarray(g, dtype=float)[:m, :m] for g in info['growth_samples']]
        with _warnings.catch_warnings():
            _warnings.simplefilter('ignore')
            v0, pinfo = pss.ppv()
        ## ⚠ `samples[j]` IS node j: prepending `v0` pairs node j's covariance
        ## with the phase vector of node j - 1 -- a one-node shift, first
        ## order in the step.
        vs = [np.asarray(sv, dtype=float)[:m] for sv in pinfo['samples']]
        jj = int(min(j, len(Ps) - 1, len(vs) - 1))

        w, U = np.linalg.eigh(G[jj])
        if not float(w.max()) > 0.0:
            raise ValueError(
                'PAC.oscillator_edge_jitter: the growth term has collapsed at '
                'this instant (largest eigenvalue %.3g), so the orbit tangent '
                'cannot be recovered from it and the phase direction cannot '
                'be projected out. Is any source noisy?' % float(w.max()))
        uj = U[:, int(np.argmax(w))] * np.sqrt(float(w.max()))
        den = float(vs[jj] @ uj)
        if den == 0.0:
            raise ValueError(
                'PAC.oscillator_edge_jitter: the left and right null '
                'directions are orthogonal at this instant, so the oblique '
                'projection is undefined.')
        Pi = np.eye(m) - np.outer(uj, vs[jj]) / den
        e = np.zeros(m)
        e[int(output)] = 1.0
        var_prj = float(e @ (Pi @ Ps[jj] @ Pi.T) @ e)
        var_raw = float(e @ Ps[jj] @ e)
        if not var_prj > 0.0:
            raise ValueError(
                'PAC.oscillator_edge_jitter: the projected variance is %.3g, '
                'not positive -- there is no additive jitter to report.'
                % var_prj)

        A = var_prj / (slew * slew)
        sigma_t = float(np.sqrt(A))
        if not sigma_t < 0.5 * T:
            raise ValueError(
                'PAC.oscillator_edge_jitter: the implied displacement is '
                'sigma_t = %.3g s, %.3g of the period -- the first-order '
                'picture (a crossing moved by delta_y/slew) does not hold '
                'there. Pass an instant on an edge, or check the source '
                'levels.' % (sigma_t, sigma_t / T))

        c = float(info['c_from_growth'])
        ks = np.arange(1, int(kmax) + 1)
        return {
            'sigma_t': sigma_t,
            'A': A,
            'c': c,
            'slew': slew,
            'k_cycle_bound': np.sqrt(c * ks * T + 2.0 * A),
            'instant': float(times[j]),
            'd': float(d),
            'projection_share': float(var_raw / var_prj - 1.0),
        }

    def orbital_mode_weights(self, pss, nmodes=None):
        """`K_orb` resolved onto the Floquet modes — A9's second step.

        ⚠⚠⚠ READ THIS FIRST: THE BASIS OMITS THE ANNIHILATED MODES, AND WHAT
        THEY CARRY IS A FLOOR NOTHING BELOW CAN GO UNDER.  `floquet_modes`
        returns the NON-NULL directions, so `sum cw[k,k'] u_k u_k'^H`
        reproduces only the part of `K_orb` that lives on them, and how much
        that is depends entirely on WHERE THE NOISE ENTERS.  The annihilated
        modes enter the stationary covariance only through the `j = 0`
        term, which is not small when the noise is injected there -- and
        that is where device noise actually is: every resistor in a bias or
        tuning network.  On `_osc_with_ladder` (`nslow = 4`) the
        reconstruction residual is 0.18 % injected at the oscillator node
        and 99.96 % injected at a FAST ladder node.

        ⚠ So a modal orbital spectrum built on this basis is complete only
        for noise that enters the slow subspace; the suite's own gate on
        this (`rel < 1e-2`) holds because its fixture injects at the
        oscillator node.

        ⚠ AND THE RESIDUAL IS A DETECTOR, NOT A TRUNCATION BOUND: it
        catches a DROPPED NON-NULL MODE, but it SATURATES at the floor above,
        so it cannot certify a truncation below whatever the null modes
        carry, however many modes are kept.

        Returns `(cw, modes, K_orb)` with `cw[k, k'] = v_k† K_orb v_k'`,
        the weight of each pair of Floquet directions in the bounded
        (orbital) part of the state covariance.

        ⚠ **THIS IS THE BRIDGE BETWEEN THE TWO ROUTES WE ALREADY OWN.**
        `oscillator_covariance` gets `K_orb` from a bordered Kronecker
        solve; `floquet_modes` gets the eigen-directions from the
        monodromy. Traversa & Bonani's eq (22) sums over exactly these
        mode pairs, so resolving the covariance we already trust onto the
        modes is the step that connects them — and, unlike the spectrum
        itself, it has an **exact identity** to check against:

            Σ_{k,k'} cw[k,k'] · u_k u_{k'}†  =  K_orb

        because `(u, v)` are biorthonormal. A wrong pairing, a wrong
        normalisation, or a dropped mode all break that reconstruction
        while leaving every individual eigenvector residual clean.

        ⚠ **THE PHASE MODE IS INCLUDED AND ITS WEIGHT SHOULD BE SMALL, NOT
        ZERO.** `K_orb` is the part of the covariance that stays bounded,
        with the along-orbit growth `n·d·uuᵀ` already removed — so the
        `k = k' = 1` entry is what the split left behind rather than a
        quantity that must vanish. Reading it as an error is a
        misinterpretation of `oscillator_covariance`'s own contract.

        ⚠ **NOT THE SPECTRUM.** `S_yy` additionally needs the Fourier
        coefficients of the periodic parts (`floquet_modes` returns them
        as `p`/`q`) and the resolvent `1/(i(j−j')ω₀ − μ_l' − μ_l*)` of
        eq (22): that is `orbital_correlation` and `orbital_spectrum`.

        History: `doc/shooting_history.md`, `PAC.orbital_mode_weights`.
        """
        ## ONE ORBIT: the covariance's host (a GLM's or trap's twin) reads
        ## the modes too, or they would come from another discretisation
        pss = pss._lyapunov_host()
        K_orb, _d, _info = self.oscillator_covariance(pss, pair=True)
        K = np.asarray(K_orb, dtype=float)
        n = K.shape[0]
        modes = pss.floquet_modes(pss, nmodes=(n if nmodes is None
                                               else int(nmodes)))
        V = np.column_stack([m['v0'] for m in modes])
        cw = V.conj().T @ K @ V
        return cw, modes, K
