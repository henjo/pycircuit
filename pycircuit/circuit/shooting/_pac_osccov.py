"""An oscillator's covariance: the transverse responses, the orbital lines, the
edge jitter and the mode weights.
"""
import numpy as np
from ._factored import dense_map
from ._numerics import _output_row, edge_slope, output_index
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

    def _transverse_responses(self, pss, fp, freqs, u_ac, u_points=None,
                              map_node0=False):
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

        `map_node0`: the bordered solution `wb` itself -- the MAP's state at
        node 0, map-width (the pair on gear), `v_pair' wb = 0` exactly by the
        border -- instead of the projected node responses; no replay
        (`orbital_mode_weights`).

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
            Md = (dense_map(fp, nw)
                  if nw <= pss.FLOQUET_DENSE_LIMIT else None)
            staged = None
            _evd = EventColumns.of(pss, nw)
            if _evd is not None:
                ## the total map; the crossings' state sensitivity; the event
                ## columns at fixed time
                if Md is not None:
                    Md = np.asarray(_evd.total_matrix(Md), dtype=float)
                staged = (_evd, np.asarray(_evd['P_end'], dtype=complex),
                          np.asarray(_evd.dth, dtype=float),
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
            if map_node0:
                out.append(np.asarray(wb, dtype=complex)[None, :])
                continue
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
                              pair=False, _coloured=True):
        """The state covariance of a FREE-RUNNING oscillator, split in two.

        Returns `(K_orb, info)`.  `K_orb` is the BOUNDED periodic
        (orbital) part of the covariance at `t = 0`; `info['d']` is the
        growth per period along the orbit tangent (see the split below).

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
        OBVIOUS READING IS WRONG.  `info['samples'][j]` is `P(t_j)`,
        the solution started from `K_orb` at `t = 0`, and it satisfies

            K(t_j + n T) = P(t_j) + n d u_j u_j^T,   u_j = Phi(t_j, 0) u

        so `P` is periodic UP TO the growth -- `P(T) = P(0) + d u u^T`, not
        `P(T) = P(0)`.  The walk is along the orbit and the orbit turns, so
        the growth DIRECTION is the propagated tangent rather than a fixed
        `u`.  On a STAGED solve both are at FIXED times -- the event
        closure's samples, and ``u_j = R_j u`` with ``R_j`` the total
        fixed-time map to node j (`_event_closure`): the plain walk would
        miss the crossings' motion.
        `info['growth_samples'][j]` is ``d u_j u_j^T`` and
        `info['tangent_samples'][j]` is ``u_j`` itself -- its first block the
        orbit's own rate at node j.

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

        The rest of `info`: `d_closed_form` and `d_residual` (above),
        `growth` (``d u u^T``), `c_from_growth` (``d / T``), `tangent_pair`
        (`u`), `ppv_pair` (`v`), `pair_inner` (``v . u``), `sigma_min` and
        `sigma_min_bordered` (the smallest singular values of ``I - M kron M``
        and of the bordered system), `null_residual` (``||(I - M kron M)
        (u kron u)||`` relative), `ppv` (`ppv()`'s info); with `samples=True`
        also `samples` (``P(t_j)``), `growth_samples`, `tangent_samples` and
        `times`.

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
        ## the host the Lyapunov surfaces read: the monodromy twin (a
        ## gear/trbdf2 run is its own; TR-BDF2's injection is built).  Swapped
        ## BEFORE both the Lyapunov pieces and `ppv` below, so the bordering
        ## keeps them on one host -- see `_lyapunov_host`.
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
        closure_samples = closure_pieces = None
        if _bordered is not None:
            M, K1, closure_samples, closure_pieces = _bordered

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

        self._check_kron(n, 'oscillator_covariance')
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
            if closure_samples is not None:
                ## ⚠ ON A STAGED SOLVE, THE CLOSURE'S OWN SAMPLES, AT FIXED
                ## TIMES (`_event_closure`): the plain walk misses the
                ## crossings' noise-driven motion and closed on `K_orb +
                ## growth` only to 3.5e-7 of its size on the comparator
                ## oscillator (5.8e-5 with `event_window_steps=2`), against
                ## 1.6e-15 for the closure's.  The growth direction
                ## at node j is the TOTAL fixed-time map's image of `u`.
                orb = [np.asarray(P, dtype=float)[:n, :n]
                       for P in closure_samples(K_orb)]
                grw, tng = [], []
                for j in range(len(orb)):
                    uj = closure_pieces['map_to'](j) @ u
                    grw.append(d * np.outer(uj[:n], uj[:n]))
                    tng.append(uj[:n])
            else:
                orb, grw, K, uj = [K_orb], [d * np.outer(u, u)], K_orb, u
                tng = [u]
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
                    tng.append(uj[:n])
            info['samples'] = orb
            info['growth_samples'] = grw
            info['tangent_samples'] = tng
            info['times'] = np.asarray(pss.factored_period().times,
                                       dtype=float)
        ## the TRANSVERSE covariance in the node space; the coloured sources'
        ## enter here alone (see the docstring)
        Pi = self._node_projectors(pss)
        info['K_transverse'] = Pi[0] @ K_orb[:m, :m] @ Pi[0].T
        if samples:
            info['transverse_samples'] = [
                Pi[j] @ Kj[:m, :m] @ Pi[j].T
                for j, Kj in enumerate(info['samples'][:len(Pi)])]
        ## (`_coloured=False`: the white part alone, for a caller that builds
        ## the coloured part itself -- `oscillator_edge_jitter`)
        if col is not None and _coloured:
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
        ## ⚠ `m x m` WHATEVER THE METHOD, as `covariance`: a
        ## two-step method's pair space is its own, `pair=True` keeps it
        ## (`orbital_mode_weights` reads the Floquet modes there)
        if not pair:
            K_orb = K_orb[:m, :m]
            info['growth'] = info['growth'][:m, :m]
            for key in ('samples', 'growth_samples'):
                if key in info:
                    info[key] = [K[:m, :m] for K in info[key]]
            if 'tangent_samples' in info:
                info['tangent_samples'] = [t[:m] for t in info['tangent_samples']]
        info['d'] = d
        return K_orb, info

    def _edge_coloured_law(self, pss, col, output, nodes, kmax,
                           projectors=None):
        """The EXACT coloured k-lag law of an oscillator's output sampled at
        grid nodes: ``V[r, k-1] = Var(y_{n+k} - y_n)``, `y_n` the output at
        node ``nodes[r]`` of period `n`, from the COLOURED components
        `_coloured_prepare` hands out (`col`; its white remainders are the
        Lyapunov path's) -- in the output's units squared (the caller
        divides by the node's rate squared).

        The ONE-SIDED PSD of the sample series, folded into ``f in (0,
        f0/2]``, is ``S(f) = sum_n CY(|f + n f0|) |G_n(f)|^2`` (the sampled
        series' convention, `_sampled_series`), and the k-lag kernel is the
        same at every alias, ``sin^2(pi (f + n f0) k T) = sin^2(pi f k T)``:

            V_k = int S(f) 4 sin^2(pi f k T) df.

        `G_n` is the sample's transfer from a source at ``f + n f0`` -- the
        FULL response, the phase pole included: one transposed solve per `f`
        gives every alias at once, the pole carried as ``1/(1 - alpha)`` by
        the bordered solve (`_deflated_solve`'s, and its rule: the plain
        operator's own answer wherever ``|1 - alpha| >=
        DEFLATION_REFINE_MIN``).  The pole's ``1/f^2`` meets the kernel's
        ``f^2``: the integrand is finite as ``f -> 0`` for a flat density,
        ``1/f`` for a 1/f one (whose band starts at `colour_fmin`).  Every
        alias is masked to the source band ``[fmin, fmax]``; the hole ``(0,
        fmin)`` the fold cuts around every harmonic (``n != 0``) is closed
        by a rectangle, its integrand being finite there.  Quadrature:
        `_lineshape.increment_nodes` (panels ``1/(2 kmax T)`` wide, log
        panels below).

        ⚠ THE TRANSVERSE PART'S MEMORY AND ITS CROSS TERM WITH THE PHASE
        ARE IN IT, with the phase: the terms the colour fold's increment
        plus ``2 A_col`` left out.

        ⚠ LINEAR IN THE SEED, SO THE REVERSE PASSES ARE PER NODE, NOT PER
        FREQUENCY: the pass seeded at node j is ``exp(-2 pi i f t_j)``
        times the unit-seeded one, the pass over the solved costate `z` a
        sum of `n` unit passes (which also give the dense transposed map),
        and on a staged solve the crossings' pass (``zeta = Gt^-T (g_theta
        + alpha P_end^T z)``, `EventColumns.collapsed_zeta`) a sum of one
        per crossing.  Per frequency only the ``(n + 1)`` solve and one
        alias product remain.  Built densely: a map up to
        `FLOQUET_DENSE_LIMIT` wide.

        Returns ``{'V': (len(nodes), kmax), 'alpha': (len(nodes),), 'gint':
        (len(nodes),), 'settled': bool}``.  THE LARGE-k INTERCEPT: near
        ``f = 0`` the folded PSD is a random walk's, ``F ~ alpha / (4
        sin^2(pi f T))`` (`F` is even in `f`, so the next term is flat), and
        with ``G = F - alpha / (4 sin^2(pi f T))`` the Fejer integral gives

            V_k = alpha k / (2 T) + 2 int G df - 2 int G cos(2 pi f k T) df

        -- a slope `alpha / (2T)` per period and the intercept ``2 int G``
        (`gint`), the cosine term vanishing as k grows.  ``alpha = 4 sin^2(pi
        f T) F`` at the lowest frequency; `settled` says it is flat there (a
        decade up it agrees to 1 %) and no component is a power law -- a
        1/f source's phase grows faster than linearly and has no intercept.

        With `projectors` (``Pi_j`` per node, `_node_projectors`'), also the
        TRANSVERSE part: the output's ``e' Pi_j y`` sampled -- a second seed,
        ``Pi_j^T e``, whose transfer has no pole -- its VARIANCE (the integral
        of its folded PSD: `transverse`, the same quantity
        `oscillator_covariance`'s coloured samples integrate over the source
        band), and the PHASE part's k-lag law, the full transfer less the
        transverse one (linear in the seed; `Vphase`)."""
        from ._lineshape import increment_nodes
        fp = pss._state_map()
        stage = fp.is_stage or fp.is_glm
        T = float(fp.T)
        f0 = 1.0 / T
        tms = np.asarray(fp.times, dtype=float)
        N = len(fp.steps)
        n = fp.width
        m = pss.cir.n - 1
        if n > pss.FLOQUET_DENSE_LIMIT:
            raise NotImplementedError(
                f'PAC.oscillator_edge_jitter: the period map is {n} wide, '
                f'above FLOQUET_DENSE_LIMIT = {pss.FLOQUET_DENSE_LIMIT}; the '
                'exact coloured k-lag law is built on the dense map.')
        L = N // 2 - 1
        ns = np.arange(-L, L + 1)
        fmin, fmax = float(col['fmin']), float(col['fmax'])
        d = _output_row(output, m)
        tinj = self._stage_times(pss, fp) if stage else tms[1:N + 1]

        def couplings(lam0=None, seed=None, extra=None):
            ## one reverse pass: its couplings, one row per injection point
            if stage:
                lam = np.zeros(m, dtype=complex) if lam0 is None else lam0
                return np.asarray(self._stage_pass(
                    pss, fp, lam, seed if seed is not None else extra)[1])
            lam = np.zeros(n, dtype=complex) if lam0 is None else lam0
            if seed is not None:
                inject = np.zeros((N, m), dtype=complex)
                inject[seed[0]] = seed[1]
            else:
                inject = extra
            return np.asarray(fp.matvec_transposed(
                np.asarray(lam, dtype=complex), collect=True, inject=inject)[1])

        ## the dense transposed map, and the couplings of a unit costate
        eye = np.eye(n)
        MT = np.column_stack([np.real(np.asarray(fp.matvec_transposed(
            eye[i].astype(complex)))) for i in range(n)])
        Tz = np.array([couplings(eye[i].astype(complex)) for i in range(n)])
        _ev = EventColumns.of(pss, n)
        if _ev is not None:
            P = np.asarray(_ev['P_end'], dtype=float)
            D = np.asarray(_ev.dth, dtype=float)
            MT = MT + D.T @ P.T                      # the TOTAL map's
            Z1 = np.linalg.inv(np.asarray(_ev['Gt'], dtype=float).T)
            Z2 = Z1 @ P.T
            Pkf = self._fixed_time_event_columns(pss)[0]
            K = Z1.shape[0]
            Tev = np.array([couplings(extra=(
                _ev.injection_dict(np.eye(K)[k]) if stage
                else _ev.injection(np.eye(K)[k], N, m))) for k in range(K)])
        ## the border: the null vectors of the (total) map, as `_deflated_solve`
        _v, pinfo = pss.ppv()
        vb = np.asarray(_v, dtype=float).ravel()
        ub = np.asarray(pinfo['tangent_pair'], dtype=float).ravel()
        bcol = vb / max(float(np.linalg.norm(vb)), 1e-300)
        brow = ub / max(float(np.linalg.norm(ub)), 1e-300)

        def solve_z(alpha, g):
            one = 1.0 - alpha
            if abs(one) >= self.DEFLATION_REFINE_MIN and not fp.is_glm:
                return np.linalg.solve(eye - alpha * MT, g)
            Bm = np.zeros((n + 1, n + 1), dtype=complex)
            Bm[:n, :n] = eye - alpha * MT
            Bm[:n, n] = bcol
            Bm[n, :n] = brow
            sol = np.linalg.solve(Bm, np.concatenate((g, [0.0])))
            return sol[:n] + sol[n] * bcol / one

        ## the coloured components as groups: FIXED columns per point with a
        ## weight per band frequency, or a STATIONARY density on a support
        ## per band frequency, or columns per point per band frequency
        w1 = col['w1']
        fixed, stationary, moving = [], [], []
        for _key, W, ef in col['comps']:
            fixed.append((np.asarray(W, dtype=complex),
                          lambda nu, ef=ef: (w1 / (2.0 * np.pi * nu)) ** ef))
        cache = {}

        def cy_at(key, xref, nu):
            k_ = (key, None if xref is None else id(xref), float(nu))
            if k_ not in cache:
                cache[k_] = np.asarray(self._noise_components(
                    pss, [col['state0'] if xref is None else xref]).one_element_cy(
                        key, 2.0 * np.pi * float(nu))[0], dtype=complex)
            return cache[k_]
        for key, W0, xref, pq, cref in col.get('separable', ()):
            fixed.append((np.asarray(W0, dtype=complex),
                          lambda nu, key=key, xref=xref, pq=pq, cref=cref:
                          np.array([float(np.real(cy_at(key, xref, x)[pq] / cref))
                                    for x in np.atleast_1d(nu)])))
        for key, supp in col.get('perband', ()):
            stationary.append((key, np.asarray(supp)))
        for _key, root_at in col.get('nonseparable', ()):
            moving.append(root_at)

        stacked = {}

        def alias_power(Sp, B, nu, keep, fi):
            ## one sampled transfer's power per alias, over every group
            pb = np.zeros(len(ns))
            for W, wt in fixed:
                R = B[keep] @ np.einsum('ji,jik->jk', Sp, W)
                pb[keep] += wt(nu[keep]) * np.sum(np.abs(R) ** 2, axis=1)
            if stationary or moving:
                Y = B[keep] @ Sp                                 # (kept, m)
                idx = np.nonzero(keep)[0]
                for key, supp in stationary:
                    ## the density at every kept alias, stacked once per
                    ## folded frequency (both nodes, every seed read it)
                    if (key, fi) not in stacked:
                        stacked[(key, fi)] = np.array(
                            [cy_at(key, None, nu[bi])[np.ix_(supp, supp)]
                             for bi in idx])
                    Ys = Y[:, supp]
                    pb[idx] += np.real(np.einsum('ai,aij,aj->a', Ys,
                                                 stacked[(key, fi)], Ys.conj()))
                for root_at in moving:
                    for bi in idx:
                        Wn = np.asarray(root_at(nu[bi]), dtype=complex)
                        R = B[bi] @ np.einsum('ji,jik->jk', Sp, Wn)
                        pb[bi] += float(np.sum(np.abs(R) ** 2))
            return pb

        def seeded(k0, dv):
            ## the reverse pass seeded at node k0 with the functional `dv`:
            ## its final costate (with the crossings' term on a staged
            ## solve), its couplings, and its theta-sensitivity
            dv = np.asarray(dv, dtype=complex)
            cA = couplings(seed=(k0, dv))
            if stage:
                g_ = np.asarray(self._stage_pass(
                    pss, fp, np.zeros(m, dtype=complex), (k0, dv))[0])
            else:
                inject = np.zeros((N, m), dtype=complex)
                inject[k0] = dv
                g_ = np.asarray(fp.matvec_transposed(
                    np.zeros(n, dtype=complex), collect=True, inject=inject)[0])
            gth = None
            if _ev is not None:
                gth = Pkf[k0].T @ dv
                g_ = g_ + D.T @ gth
            return g_, cA, gth

        ## ⚠ THE ORBITAL LINES, FOLDED INTO (0, f0/2]: a lightly damped
        ## mode's line is narrower than a panel and falls between its
        ## points -- a Q = 1000 resonator's read 29 % low (`increment_nodes`'
        ## `lines`); every harmonic of a mode folds to one place
        folded = set()
        for c0, hw in self._orbital_lines(pss, fmin, fmax):
            x = float(c0) % f0
            folded.add((round(min(x, f0 - x) / f0, 12) * f0,
                        round(float(hw) / f0, 12) * f0))
        fq, wq = increment_nodes(fmin, 0.5 * f0, T, int(kmax),
                                 lines=sorted(folded))
        ks = np.arange(1, int(kmax) + 1)
        ker = 4.0 * np.sin(np.pi * fq[None, :] * ks[:, None] * T) ** 2
        V = np.zeros((len(nodes), len(ks)))
        Vph = np.zeros((len(nodes), len(ks)))
        var_t = np.zeros(len(nodes))
        alphas, gints, settled = [], [], not col['comps']
        s2f = 4.0 * np.sin(np.pi * fq * T) ** 2
        i10 = min(int(np.searchsorted(fq, 10.0 * fq[0])), len(fq) - 1)
        for r_, j in enumerate(nodes):
            k0 = int(j) % N
            ## the output, and with `projectors` its transverse part
            dvs = [d] if projectors is None else [
                d, np.asarray(projectors[r_], dtype=float).T @ d]
            seeds = [seeded(k0, dv) for dv in dvs]
            B = np.exp(2j * np.pi * ns[:, None] * f0
                       * (tinj[None, :] - tms[k0]))              # (2L+1, K)
            nparts = 1 if projectors is None else 3              # full, transverse, phase
            per_f = np.zeros((nparts, len(fq)))
            hole = np.zeros(nparts)
            for fi, f in enumerate(fq):
                alpha = np.exp(-2j * np.pi * f * T)
                ph = np.exp(-2j * np.pi * f * tms[k0])
                at_f = np.exp(2j * np.pi * f * tinj)[:, None]
                Sps = []
                for g0, cA0, gth0 in seeds:
                    z = solve_z(alpha, ph * g0)
                    Sv = -(ph * cA0 + alpha * np.tensordot(z, Tz, axes=1))
                    if _ev is not None:
                        zeta = Z1 @ (ph * gth0) + alpha * (Z2 @ z)
                        Sv = Sv - np.tensordot(zeta, Tev, axes=1)
                    Sps.append(at_f * Sv)
                if projectors is not None:
                    Sps.append(Sps[0] - Sps[1])                  # the phase part
                nu = np.abs(f + ns * f0)
                keep = (nu >= fmin) & (nu <= fmax * (1.0 + 1e-12))
                for q, Sp in enumerate(Sps):
                    pb = alias_power(Sp, B, nu, keep, fi)
                    per_f[q, fi] = float(np.sum(pb))
                    if fi == 0:
                        ## (the hole below fmin around every harmonic n != 0)
                        hole[q] = float(np.sum(pb[ns != 0]))
            V[r_] = (ker * (wq * per_f[0])[None, :]).sum(axis=1) + fmin * ker[:, 0] * hole[0]
            if projectors is not None:
                var_t[r_] = float(np.sum(wq * per_f[1])) + fmin * hole[1]
                Vph[r_] = ((ker * (wq * per_f[2])[None, :]).sum(axis=1)
                           + fmin * ker[:, 0] * hole[2])
            ## the random walk's slope and the intercept's remainder
            al = float(s2f[0] * per_f[0, 0])
            a10 = float(s2f[i10] * per_f[0, i10])
            settled = settled and (abs(a10 - al) <= 1e-2 * abs(al) if al != 0.0
                                   else a10 == 0.0)
            G = per_f[0] - al / s2f
            alphas.append(al)
            gints.append(float(np.sum(wq * G) + fmin * G[0]))
        return {'V': V, 'alpha': np.asarray(alphas), 'gint': np.asarray(gints),
                'settled': bool(settled), 'transverse': var_t, 'Vphase': Vph}

    def oscillator_edge_jitter(self, pss, output, time, kmax=8,
                               colour_fmin=None, colour_fmax=None,
                               points_per_decade=40, intercept='white'):
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
        by `delta_y/slew`, so with `s` the slope at `t_j` and `M_j` the period
        map from `t_j` (`M_j u_j = u_j`), the k-lag law a designer measures is
        EXACTLY (linear noise)

            Var(tau_{n+k} - tau_n) = e^T (2 P_j + k G_j - M_j^k P_j - P_j M_j^k^T) e / s^2

        (`k_cycle`, `G_j = d u_j u_j^T`), and as `k` grows `M_j^k -> u_j v_j^T /
        (v_j^T u_j)`, so it tends to `c k T + 2 sigma_t^2` with

            sigma_t^2 = e^T Pi P(t_j) e / s^2,  Pi = I - u_j v_j^T/(v_j^T u_j)

        -- the ONE-sided projection, which keeps the cross term
        ``X = e^T Pi P (I - Pi)^T e / s^2``, the correlation between the
        transverse and the phase deviation at the edge (X/A = -0.150 on the
        A11 fixture; the two-sided ``Pi P Pi^T`` drops it).  A committed
        Monte Carlo (`benchmarks/oscillator_edge_jitter_probe.py`: noisy
        radau transients at the PSS's step, 5760 crossings) agrees with the
        exact law at k = 1..8 (1.2-1.7 sigma).

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

        `u_j` is `oscillator_covariance`'s `tangent_samples[j]`, the
        tangent carried to node j; `Pi` is invariant to its scale (and to
        `v`'s), so taking `v` from `ppv()` in a separate call is safe here.
        Everything is sliced `[:m, :m]` out of PAIR space.

        ⚠ THE LAW IS SELF-CONSISTENT PER NODE AND LINEAR IN THE INSTANT
        BETWEEN NODES.  At node j every piece is node j's -- `P_j`, `G_j`,
        `M_j` and the slope ``s_j = e . u_j``, the orbit's own rate there, so
        ``e^T G_j e / s_j^2 = c T`` exactly -- and the law at the requested
        instant is the linear blend of the two nodes around it: 1.5182 /
        1.5181 / 1.5180e-6 at 240 / 480 / 960 points (k = 1, the A11 chain).
        `slew` in the result is the local quadratic fit AT the instant
        (`edge_slope`), reported: a two-point difference moves with the grid.

        ⚠ When validating this against a transient, run the Monte Carlo on the
        SAME integrator as the PSS, or divide by the slope of the orbit the
        Monte Carlo actually runs on -- otherwise the mismatch enters squared
        (an Euler orbit's slope against gear's: 5.7 % at 240 points).

        `k_cycle` is exact at every `k` -- the same key and meaning as
        `jitter_metrics`' for a DRIVEN circuit.  On a STAGED solve (state
        events) the period map from node `j`
        carries the crossings' motion: ``M_j^k = R_j M_tot^{k-1} S_j`` from
        the event closure (`_event_closure`'s `period_map_from`), `P_j` and
        `G_j` its fixed-time samples.  ⚠ Gated by its construction only: on
        the one staged oscillator here the switch's transition sits inside
        the landed window, and the event term moves `k_cycle^2` by 2e-7.  A
        Nordsieck GLM's native-width map is refused.

        ⚠ THE INSTANT IS THE CALLER'S.  This does not hunt for a crossing: a
        threshold taken from a simulated record can be biased by startup, and
        that moves the instant off the steepest point.  `instant` in the
        result is the instant the law describes (`time` modulo the period),
        `nodes` the two grid nodes around it and `th` its fraction between
        them.

        ⚠ A COLOURED SOURCE needs its band, as in `oscillator_covariance`:
        `colour_fmin` (and `colour_fmax`, default the grid's Nyquist).  Its
        part of `k_cycle` is EXACT: the output sampled once a period at the
        edge's nodes, its one-sided PSD folded into (0, f0/2] and integrated
        against the k-lag kernel (`_edge_coloured_law`) -- the phase with its
        memory (the increment depends on WHERE in the period the edge sits),
        the transverse part with its correlation across periods, and their
        cross term, together; per node over the node's own rate, blended as
        the white law is.  `coloured_variance` in the result is that part
        (s^2); `coloured_phase_variance` is its PHASE part alone -- the full
        transfer less the transverse one, the oblique split at the edge's
        nodes (reported, not added) -- and `coloured_transverse_variance` its
        transverse variance, both from the same folded solve.

        ⚠ WHAT `A` AND `sigma_t` MEAN WITH A COLOURED SOURCE is `intercept`'s
        choice.  A coloured source's TRANSVERSE variance is, for a 1/f
        source, a slow wander rather than additive jitter (2A = 7.2e-6 s^2
        against k_cycle_1^2 = 1.0e-7 s^2 on an orbit-modulated 1/f van der
        Pol, where that wander cancels in the increments), so by default it
        is not in them:

          'white' (default)  the WHITE sources' exact intercept; a coloured
                             source enters `k_cycle` only (`A = 0`,
                             `sigma_t = 0` with no white source).
          'exact'            the whole law's large-k intercept: the white
                             part's plus the coloured part's, from the folded
                             PSD's random-walk remainder (`_edge_coloured_law`)
                             -- nan, warned, where a coloured source has none
                             (a power law: its phase grows faster than
                             linearly in k).  A Lorentzian's equals its exact
                             white realisation's.

        `coloured_transverse_variance` (s^2) is the coloured transverse
        variance at the edge, and `c_coloured` the coloured sources' linear
        diffusion (the slope of `k_cycle^2` per `k T`; nan for a power law),
        both reported in either case.  Without a coloured source the two
        choices are the same.

        Returns a dict: `sigma_t`, `A` (= sigma_t^2), `c`, `slew`,
        `k_cycle` (k = 1..kmax), `instant`, `nodes`, `th`, `d`,
        `projection_share` (the
        white part's; None without one), `coloured_variance` (the coloured
        part of `k_cycle^2` per k, s^2; zeros when white),
        `coloured_transverse_variance` (s^2), `c_coloured` (s; 0 when white),
        `coloured_phase_variance` (the coloured phase part's k-lag
        variance per k, s^2; zeros when white), `band` (the colour band, or
        None).

        History: `doc/shooting_history.md`, `PAC.oscillator_edge_jitter`.
        """
        output = output_index(pss, output)
        if intercept not in ('white', 'exact'):
            raise ValueError(
                "PAC.oscillator_edge_jitter: intercept must be 'white' (the "
                "white sources' exact intercept) or 'exact' (the whole law's, "
                f"nan for a power-law source); got {intercept!r}.")
        self._check_circuit(pss)
        ## ONE ORBIT: the covariance's host (a GLM's or trap's twin) supplies
        ## the period, the factored period and the PPV as well -- read off
        ## the run itself they come from another discretisation
        ## History: `doc/shooting_history.md`, `PAC.oscillator_edge_jitter`.
        pss = pss._lyapunov_host()
        ## (the covariance in the MAP's own space: the exact law contracts
        ## it with the period map; the node block is `[:m, :m]` of it)
        with warnings.catch_warnings():
            ## (its "K_orb, d and c_from_growth are the WHITE sources' alone"
            ## is what this method handles below)
            warnings.filterwarnings('ignore', message='PAC.oscillator_covariance: '
                                     'this circuit has a COLOURED source')
            ## (the WHITE part alone: the coloured one -- its transverse
            ## variance and phase included -- is `_edge_coloured_law`'s, one
            ## folded solve; the covariance's coloured samples at every node
            ## would cost 29 s of a 67 s call)
            K_orb, info = self.oscillator_covariance(
                pss, samples=True, pair=True, colour_fmin=colour_fmin,
                colour_fmax=colour_fmax, points_per_decade=points_per_decade,
                _coloured=False)
        coloured = bool(self._coloured_present(pss))
        d = info['d']
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

        ## ⚠ THE LAW AT THE REQUESTED INSTANT, from the two nodes around it,
        ## each SELF-CONSISTENT: `P_j`, `G_j`, `M_j` and the slope `e . u_j`
        ## (the orbit's own rate at the node) at ONE node, then linear in the
        ## instant.  The slope at the instant (`edge_slope`, a local
        ## quadratic fit) is reported as `slew`.
        tc = float(time) % T
        slew, _jn = edge_slope(times, row, tc)
        if not abs(slew) > 0.0:
            raise ValueError(
                'PAC.oscillator_edge_jitter: the slope at t = %.6g is exactly '
                'zero, so delta_y/slew is undefined. Pass an instant on an '
                'edge.' % tc)
        N = len(info['samples']) - 1
        tn = np.asarray(fp.times, dtype=float)[:N + 1]
        ## node N is node 0 a period on (the law is the same there: the
        ## growth `d u u^T` cancels in it), so an instant in the last step
        ## blends N-1 and N
        a = int(np.clip(np.searchsorted(tn, tc, side='right') - 1, 0, N - 1))
        b = a + 1
        th = float(np.clip((tc - tn[a]) / (tn[b] - tn[a]), 0.0, 1.0))

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _v0, pinfo = pss.ppv()
        ## (the PPV samples cover nodes 0..N-1; node N is node 0)
        vs = [np.asarray(sv, dtype=float)[:m] for sv in pinfo['samples']]
        ## the output as a weight vector (a node, or a differential output)
        e = _output_row(output, m)
        white = float(d) > 0.0
        if not white and not coloured:
            raise ValueError(
                'PAC.oscillator_edge_jitter: the growth per period is zero and '
                'no source is coloured -- is any source noisy?')

        c = float(info['c_from_growth'])
        col = None
        if coloured:
            col = self._coloured_prepare(pss, colour_fmin, colour_fmax,
                                         points_per_decade,
                                         'oscillator_edge_jitter')
        As, Qs, _K1, M, _m, n = self._lyapunov_pieces(
            pss, 'oscillator_edge_jitter',
            white=None if col is None else col['white'])
        na = np.asarray(As[0]).shape[0] if As else n
        width = np.asarray(info['samples'][0]).shape[0]
        if na != width:
            raise NotImplementedError(
                f'PAC.oscillator_edge_jitter: the period map is {na} wide and '
                f'the covariance {width} -- a Nordsieck GLM read on its '
                'native map; solve with the default monodromy (its radau '
                'twin).')
        ## ⚠ ON A STAGED SOLVE the period map from node `j` carries the
        ## crossings' motion, ``M_j^k = R_j M_tot^{k-1} S_j`` (`_event_closure`;
        ## no inverse -- the step maps of a DAE are singular), and the
        ## samples are its fixed-time ones; otherwise the step maps' product
        with warnings.catch_warnings():
            ## (a host it cannot border, `oscillator_covariance` warned above)
            warnings.filterwarnings('ignore', message='PAC.covariance: the '
                                     'solve is staged on its state events')
            staged = self._event_closure(pss, As, Qs, M, _m, n)
        ef = _output_row(e, na)

        def at_node(jn):
            """Node `jn`'s law over its own slope: `(Var_k for k = 1..kmax,
            the projected, raw and coloured transverse variance)`, s^2."""
            uj = np.asarray(info['tangent_samples'][jn], dtype=float)[:m]
            s2 = float(e @ uj) ** 2
            if not s2 > 0.0:
                raise ValueError(
                    'PAC.oscillator_edge_jitter: the orbit is flat in this '
                    f'output at node {jn} (t = {tn[jn]:.6g}), next to the '
                    'instant, so delta_y/slew is undefined there. Pass an '
                    'instant on an edge.')
            prj = raw = 0.0
            if white:
                vj = vs[jn % len(vs)]
                den = float(vj @ uj)
                if den == 0.0:
                    raise ValueError(
                        'PAC.oscillator_edge_jitter: the left and right null '
                        f'directions are orthogonal at node {jn}, so the '
                        'oblique projection is undefined.')
                ## the ONE-sided projection: the exact law's large-k
                ## intercept / 2 (the two-sided `Pi P Pi^T` would drop the
                ## cross term)
                Pm = np.asarray(info['samples'][jn], dtype=float)[:m, :m]
                Pi = np.eye(m) - np.outer(uj, vj) / den
                prj = float(e @ (Pi @ Pm) @ e)
                raw = float(e @ Pm @ e)
            ## (a coloured source's transverse part: `_edge_coloured_law`'s)
            var_col = 0.0
            ## ⚠ THE EXACT k-LAG LAW, in the covariance's own space (a pair
            ## map's `(x_n, x_{n-1})` on gear)
            Pf = np.asarray(info['samples'][jn], dtype=float)
            Gf = np.asarray(info['growth_samples'][jn], dtype=float)
            if staged is None:
                Mj = np.eye(na)
                for i in list(range(jn, len(As))) + list(range(jn)):
                    Mj = np.asarray(As[i], dtype=float) @ Mj
            else:
                M_tot = staged[0]
                Rj, Sj = staged[3]['period_map_from'](jn)
            kc, Mk, Mp = [], np.eye(na), np.eye(na)
            for k in range(1, int(kmax) + 1):
                if staged is None:
                    Mk = Mj @ Mk
                else:
                    Mk = Rj @ Mp @ Sj
                    Mp = M_tot @ Mp
                kc.append(float(ef @ (2.0 * Pf + k * Gf - Mk @ Pf - Pf @ Mk.T) @ ef))
            return np.asarray(kc) / s2, prj / s2, raw / s2, var_col / s2

        la, lb = at_node(a), at_node(b)
        kc, A_prj, A_raw, A_col = ((1.0 - th) * x + th * y for x, y in zip(la, lb))
        if not coloured and A_prj == 0.0:
            raise ValueError(
                'PAC.oscillator_edge_jitter: the projected variance is zero '
                '-- there is no additive jitter to report.')

        inc = np.zeros(len(kc))
        cvar = np.zeros(len(kc))
        band = None
        A_cx = c_col = 0.0
        if coloured:
            band = (float(col['fmin']), float(col['fmax']))
            ## ⚠ THE EXACT COLOURED k-LAG LAW at both nodes, each over its
            ## own rate (`_edge_coloured_law`), the transverse part's memory
            ## and its cross term with the phase included
            Pi = self._node_projectors(pss)
            law = self._edge_coloured_law(pss, col, output, (a, b),
                                          int(kmax), projectors=(Pi[a], Pi[b]))
            sa, sb = (float(e @ np.asarray(info['tangent_samples'][jn],
                                           dtype=float)[:m]) ** 2
                      for jn in (a, b))
            Vn = law['V']
            cvar = (1.0 - th) * Vn[0] / sa + th * Vn[1] / sb
            kc = kc + cvar
            ## the coloured part's own intercept and slope, where they exist
            if law['settled']:
                A_cx = (1.0 - th) * law['gint'][0] / sa + th * law['gint'][1] / sb
                c_col = ((1.0 - th) * law['alpha'][0] / sa
                         + th * law['alpha'][1] / sb) / (2.0 * T * T)
            else:
                A_cx = c_col = float('nan')
            ## the coloured TRANSVERSE variance and the PHASE part's k-lag
            ## variance, reported: the same solve's
            A_col = (1.0 - th) * law['transverse'][0] / sa + th * law['transverse'][1] / sb
            inc = (1.0 - th) * law['Vphase'][0] / sa + th * law['Vphase'][1] / sb

        A = A_prj if intercept == 'white' else A_prj + A_cx
        if np.isnan(A):
            warnings.warn(
                'PAC.oscillator_edge_jitter: intercept=\'exact\' and a '
                'coloured source here is a POWER LAW, whose phase grows faster '
                'than linearly in k: the k-cycle law has no large-k intercept, '
                'so A and sigma_t are nan; k_cycle is exact.',
                RuntimeWarning, stacklevel=2)
            sigma_t = float('nan')
        elif A == 0.0:
            ## (no white source, intercept='white': none of the jitter is
            ## additive white jitter; the coloured part is in `k_cycle`)
            sigma_t = 0.0
        elif A > 0.0:
            sigma_t = float(np.sqrt(A))
            if not sigma_t < 0.5 * T:
                raise ValueError(
                    'PAC.oscillator_edge_jitter: the implied displacement is '
                    'sigma_t = %.3g s, %.3g of the period -- the first-order '
                    'picture (a crossing moved by delta_y/slew) does not hold '
                    'there. Pass an instant on an edge, or check the source '
                    'levels.' % (sigma_t, sigma_t / T))
        else:
            ## ⚠ A NEGATIVE INTERCEPT IS NOT AN ERROR: `A` is the exact
            ## law's large-k intercept / 2, and where the transverse and the
            ## phase deviation at the edge are ANTI-correlated (a slow state
            ## that carries the noise the phase later takes up: white noise
            ## through an RC into the tank reads X/A = -1.94) the intercept is
            ## negative while every `Var_k` is positive -- the walk dominates.
            ## There is then no additive variance; `k_cycle` is exact anyway.
            warnings.warn(
                'PAC.oscillator_edge_jitter: the k-cycle law\'s intercept is '
                'NEGATIVE here (A = %.3g s^2): the transverse and the phase '
                'deviation at this edge are anti-correlated, so there is no '
                'additive variance (sigma_t is nan); k_cycle is exact.' % A,
                RuntimeWarning, stacklevel=2)
            sigma_t = float('nan')

        return {
            'sigma_t': sigma_t,
            'A': A,
            'c': c,
            'slew': slew,
            'k_cycle': np.sqrt(np.clip(kc, 0.0, None)),
            'instant': tc,
            'nodes': (a, b),
            'th': th,
            'd': float(d),
            'projection_share': (float(A_raw / A_prj - 1.0) if white
                                 else None),
            'coloured_variance': cvar,
            'coloured_transverse_variance': float(A_col),
            'c_coloured': float(c_col),
            'coloured_phase_variance': inc,
            'band': band,
        }

    def orbital_mode_weights(self, pss, nmodes=None, colour_fmin=None,
                             colour_fmax=None, points_per_decade=40):
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

        Returns `(cw, info)` with `cw[k, k'] = v_k† K_orb v_k'`, the weight
        of each pair of Floquet directions in the bounded (orbital) part of
        the state covariance, `info['modes']` those directions and
        `info['K']` that covariance (pair space on a pair map).

        ⚠ A COLOURED SOURCE needs its band, as in
        `oscillator_covariance` (`colour_fmin`, `colour_fmax`,
        `points_per_decade`).  Its part of the bounded covariance is built in
        the MAP's own space from the bordered solution at node 0
        (`_transverse_responses(map_node0=True)`): the border makes
        `v_pair' wb = 0` exactly, so the phase row of the coloured weights is
        zero by construction (stacking node-projected responses into a pair
        leaks it at first order).  `info['K_coloured']`,
        `info['cw_coloured']` hold it; `cw` and `info['K']` the total.  A
        staged solve is refused for the coloured part (the map-space
        crossing terms are not built).

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
        with warnings.catch_warnings():
            ## (its "K_orb ... the WHITE sources' alone" is handled below)
            warnings.filterwarnings('ignore', message='PAC.oscillator_covariance: '
                                    'this circuit has a COLOURED source')
            K_orb, _info = self.oscillator_covariance(
                pss, pair=True, colour_fmin=colour_fmin,
                colour_fmax=colour_fmax, points_per_decade=points_per_decade)
        K = np.asarray(K_orb, dtype=float)
        n = K.shape[0]
        modes = pss.floquet_modes(pss, nmodes=(n if nmodes is None
                                               else int(nmodes)))
        V = np.column_stack([m['v0'] for m in modes])
        cw = V.conj().T @ K @ V
        info = {'modes': modes, 'K': K}
        if 'K_coloured' in _info:
            if EventColumns.of(pss) is not None:
                raise NotImplementedError(
                    'PAC.orbital_mode_weights: a COLOURED source on a STAGED '
                    'solve -- the map-space crossing terms are not built.')
            col = self._coloured_prepare(pss, colour_fmin, colour_fmax,
                                         points_per_decade,
                                         'orbital_mode_weights')
            try:
                Kc, _none = self._coloured_covariance(
                    pss, col, self.cir.n - 1, n,
                    responses=lambda *a, **k: self._transverse_responses(
                        *a, map_node0=True, **k),
                    lines=self._orbital_lines(pss, col['fmin'], col['fmax']),
                    all_nodes=False, map_node0=True)
            finally:
                self._transverse_cache = None
            Kc = np.asarray(Kc[0], dtype=float)
            cw_c = V.conj().T @ Kc @ V
            info.update(K_coloured=Kc, cw_coloured=cw_c, K=K + Kc)
            cw = cw + cw_c
        return cw, info
