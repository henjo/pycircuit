"""The periodic covariance by a Lyapunov walk over the period: the per-step
pieces, the coloured band integral, event jitter.
"""
import numpy as np
import warnings
from ._noise_components import (exponent_columns, psd_sqrt,
                               uniform_exponent, warn_signed_unused)


class _LyapunovCovariance(object):
    """The periodic covariance by a Lyapunov walk over the period: the per-
    step pieces, the coloured band integral, event jitter.  A theme of `PAC`
    (see `pac.py`)."""

    def _lyap_cy(self, pss, w, xr, white=None):
        """`CY` as the Lyapunov pieces read it: `_cy_at`, or for a COLOURED
        covariance the WHITE part `A(x)` of the component model alone
        (`white`, `_coloured_prepare`'s) -- the coloured part is added in
        the frequency domain (`_coloured_covariance`), and `CY` at `w0`
        would count it again."""
        if white is None:
            return self._cy_at(pss, w, xr)
        return white(xr)

    def _lyapunov_pieces(self, pss, what, white=None):
        """The per-step maps, injections and one-period accumulation.

        Returns `(As, Qs, K1, M, m, n)`: the step maps `A_j`, the noise
        injections `Q_j`, the covariance `K1` reached after one period
        starting from zero, the monodromy `M`, and the two widths.  `white`:
        the white part of a COLOURED source's `CY` (`_coloured_prepare`),
        read in place of `CY`; without it a coloured source is refused.

        ⚠ SHARED BY THE DRIVEN AND AUTONOMOUS ROUTES ON PURPOSE.  The two
        differ only in what they do with `I - M kron M`: `covariance`
        inverts it, `oscillator_covariance` borders it because it is
        singular there.  Everything upstream -- the `CY/2` convention, the
        `b = 0` restriction, the `C` ring the forward recursion sees -- is
        one implementation, so the pair cannot drift apart over exactly this
        factor of two.

        History: `doc/shooting_history.md`, `PAC._lyapunov_pieces`.
        """
        if white is None:
            ## (a coloured `covariance` / `event_jitter` hands in the WHITE
            ## part for the pieces and adds the coloured one itself)
            self._refuse_coloured(pss, what)
        fp = pss.factored_period()
        if fp.is_glm:
            ## reached with monodromy='native' only (`_lyapunov_host`)
            return self._lyapunov_pieces_glm(pss, fp.state_map(), what, white)
        if fp.is_stage:
            ## the stage method's per-step map + its stage injection (or the
            ## SAME exact Van Loan integral) -- see `_lyapunov_pieces_stage`
            return self._lyapunov_pieces_stage(pss, fp, what, white)
        if not fp.is_pair:
            return self._lyapunov_pieces_plain(pss, fp, what, white)
        m = pss.cir.n - 1
        n = fp.width
        hs = np.diff(np.asarray(fp.times, dtype=float))
        ## ⚠ `CY` PER STEP, AT THE STEP'S OWN STATE -- not one `CY` for the
        ## period: a MODULATED source (a switch's `4kT g(t)`, a MOS channel's
        ## `4kT gamma gd0(t)`) is then inside the formulation, and needs no
        ## cyclostationarity refusal here.  Evaluated at the state the
        ## step's companion was factored at (the implicit step's own
        ## solution); the colour refusal still applies -- colour is a
        ## different axis.
        w0 = 2.0 * np.pi / float(fp.T)
        _W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                       pss.irefnode, axis=0)
        cys = [np.real(self._lyap_cy(pss, w0,
                                     _W[:, min(k + 1, _W.shape[1] - 1)],
                                     white))
               for k in range(len(fp.steps))]

        ## the C ring as the forward recursion sees it -- see the replays
        cs0, cs1, ring = [], [], list(fp.opening)
        for _lu, C_new, _a, _b in fp.steps:
            cs0.append(ring[0])
            cs1.append(ring[1])
            ring = [C_new, ring[0]]

        def step_map(k):
            lu, _Cn, alphas, b = fp.steps[k]
            if b:
                raise NotImplementedError(
                    'PAC.%s: derived for a b = 0 companion (Gear-2).' % what)
            ## a one-step companion (gear's Euler backstop past the
            ## zero-stability bound on an event grid) has no third alpha: the
            ## pair map holds with it zero -- see
            ## `_monodromy_matvec_transposed`
            a2 = float(alphas[2]) if len(alphas) > 2 else 0.0
            A = np.zeros((n, n))
            for j in range(n):
                p0 = np.zeros(m)
                p1 = np.zeros(m)
                (p0 if j < m else p1)[j if j < m else j - m] = 1.0
                A[:m, j] = -lu.solve(alphas[1] * (cs0[k] @ p0)
                                     + a2 * (cs1[k] @ p1))
                A[m:, j] = p0
            return A

        As, Qs = [], []
        for k, (lu, _Cn, _a, _b) in enumerate(fp.steps):
            ## Q = Jf^-1 (CY / 2h) Jf^-T, symmetrised against round-off
            half = cys[k] / (2.0 * hs[k])
            left = np.column_stack([lu.solve(half[:, j]) for j in range(m)])
            Q1 = np.column_stack([lu.solve(left[j, :]) for j in range(m)]).T
            Q = np.zeros((n, n))
            Q[:m, :m] = 0.5 * (Q1 + Q1.T)
            Qs.append(Q)
            As.append(step_map(k))

        K = np.zeros((n, n))
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
        M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
        return As, Qs, K, M, m, n

    def _lyapunov_pieces_plain(self, pss, fp, what, white=None):
        """`_lyapunov_pieces` for the PLAIN path — the one-step companions.

        ⚠ THE PER-STEP STATE DEPENDS ON THE METHOD, and that is the whole
        content of this routine.  `_monodromy_matvec_plain` writes every
        one-step companion as

            S    = a1 C_{k-1} x_{k-1} + b iq_{k-1}
            x_k  = -K S,           K = Jf_k^-1
            iq_k = a0 C_k x_k + S

        **Euler** (`b = 0`): `iq` never re-enters, the state is `x` alone,
        `A_k = -a1 K C_{k-1}` is `m x m`, and nothing downstream changes --
        `n = m`, `M = fp.matvec`, and `ppv()`'s width-`m` vectors border
        it directly.

        **Trapezoidal** (`b = -1`): `iq` DOES re-enter, so the per-step
        state is the PAIR `(x, iq)`, `A_k` is `2m x 2m`, and the noise --
        which enters the KCL rows and reaches `x_k` through `K` -- reaches
        `iq_k` through `a0 C_k K` as well:

            G_k = [ K ; a0 C_k K ],     Q_k = G_k (CY/2h_k) G_k^T

        The period map on that pair RE-SEEDS `iq` at zero at the boundary,
        as the shooting solve does (the manufactured opener, B16): it is
        the product of the `A_k` applied to `(x, 0)`.  Its `x -> x` block
        IS `fp.matvec`, and that tie is the gate.  Carrying `iq` across the
        boundary instead makes `I - M kron M` singular (see the comment at
        the end).

        ⚠ `oscillator_covariance` on TRAP-PLAIN borders `I - M kron M` with
        the PAIR's null vectors, `2m` wide where `ppv()`'s are `m`: since the
        map re-seeds `iq`, they are `[v; 0]` and `M[:, :m] u` (see
        `oscillator_covariance`).  It reaches this only under
        `monodromy='native'` -- a trap oscillator otherwise reads a twin --
        and is as accurate as trap's own map, first order on a limit cycle.

        History: `doc/shooting_history.md`, `PAC._lyapunov_pieces_plain`.
        """
        m = pss.cir.n - 1
        hs = np.diff(np.asarray(fp.times, dtype=float))
        ## ⚠ `CY` PER STEP, AT THE STEP'S OWN STATE -- as in
        ## `_lyapunov_pieces`: a modulated source is inside the formulation;
        ## the colour refusal still applies.
        w0 = 2.0 * np.pi / float(fp.T)
        _W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                       pss.irefnode, axis=0)
        cys = [np.real(self._lyap_cy(pss, w0,
                                     _W[:, min(k + 1, _W.shape[1] - 1)],
                                     white))
               for k in range(len(fp.steps))]
        C_open = np.asarray(fp.opening[0], dtype=float)
        prevC = [C_open] + [np.asarray(st[1], dtype=float)
                            for st in fp.steps[:-1]]
        bs = {bool(st[3]) for st in fp.steps}
        if len(bs) != 1:
            raise NotImplementedError(
                'PAC.%s: the plain period mixes b = 0 and b != 0 steps, '
                'which have different per-step states.' % what)
        pair = bs.pop()
        n = 2 * m if pair else m
        As, Qs = [], []
        for k, (lu, C_new, alphas, b) in enumerate(fp.steps):
            Ck = np.asarray(C_new, dtype=float)
            Cp = prevC[k]
            A = np.zeros((n, n))
            for j in range(n):
                e = np.zeros(n)
                e[j] = 1.0
                p0, p1 = e[:m], (e[m:] if pair else None)
                S = alphas[1] * (Cp @ p0)
                if pair:
                    S = S + b * p1
                x = -np.asarray(lu.solve(S), dtype=float)
                A[:m, j] = x
                if pair:
                    A[m:, j] = alphas[0] * (Ck @ x) + S
            ## noise: K (CY/2h) K^T on the state block, built the same way
            ## as the solved-history route so the two cannot drift apart
            half = cys[k] / (2.0 * hs[k])
            left = np.column_stack([lu.solve(half[:, j]) for j in range(m)])
            Q1 = np.column_stack([lu.solve(left[j, :]) for j in range(m)]).T
            Q1 = 0.5 * (Q1 + Q1.T)
            Q = np.zeros((n, n))
            Q[:m, :m] = Q1
            if pair:
                Bk = alphas[0] * Ck
                Q[:m, m:] = Q1 @ Bk.T
                Q[m:, :m] = Bk @ Q1
                Q[m:, m:] = Bk @ Q1 @ Bk.T
                Q = 0.5 * (Q + Q.T)
            As.append(A)
            Qs.append(Q)
        K = np.zeros((n, n))
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
        if pair:
            ## ⚠ THE PERIOD MAP RE-SEEDS THE COMPANION, AND THAT IS
            ## LOAD-BEARING.  The plain product of the A_k carries `iq`
            ## across the boundary, and its `I - M kron M` is SINGULAR:
            ## trapezoidal maps an algebraic row's companion by exactly -1
            ## per step, so the un-reset pair carries a marginal mode (the
            ## `(-1)^n` obstruction of every formulation that keeps `iq`
            ## across a period).  The shooting solve is well-posed because
            ## the manufactured opener re-seeds `iq` at zero; the
            ## covariance's period map must do the same.  With `iq` zeroed
            ## at the start, the x->x block of the product IS `fp.matvec`,
            ## and the map on the pair is the product applied to `(x, 0)`.
            Mp = np.eye(n)
            for A in As:
                Mp = A @ Mp
            M = np.zeros((n, n))
            M[:, :m] = Mp[:, :m]
        else:
            M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                                 for e in np.eye(n)])
        return As, Qs, K, M, m, n

    def _vanloan_step_injection(self, Cr, Gr, CYr, h):
        """The per-step process-noise covariance `Q_n` for TR-BDF2, by the
        DAE-projected VAN LOAN integral.

        For ADDITIVE (linearised) noise the injection is the DETERMINISTIC
        integral `Q = integral_0^h Phi(h,s) D Phi(h,s)^T ds` -- the Levy
        areas vanish, so there are no stochastic stage weights to derive
        (Roemisch & Winkler; a naive two-stage scheme is 27 % biased on
        kT/C).  Van Loan evaluates it exactly: the
        upper-right block of `expm([[-A, D],[0, A^T]] h)` premultiplied by
        the flow.

        ⚠ BUT MNA IS A DAE (`C` singular), and the nilpotent block
        DIFFERENTIATES white noise -- discretised white noise has variance
        `S/h`, so a covariance formed on an algebraic row diverges as `1/h`
        (measured).  So Van Loan is applied on the DIFFERENTIAL SUBSPACE
        only (the capacitive nodes -- Demir 1996 propagates exactly there),
        after eliminating the algebraic variables by their Schur complement.
        The algebraic noise is routed to the differential rows through the
        same elimination (`R_proj`), so a source with a capacitive path
        (Winkler's `im A_N subset im A_C`) is handled; a source on a bare
        constraint has no differential image and is dropped rather than
        divergently amplified -- the projection is structurally immune to
        the `1/h` blow-up.

        ⚠ Against kT/C the stationary error is the METHOD's O(h^2), NOT
        machine zero: a machine-zero kT/C would mean a method-consistent
        `Q = P(1-A^2)` fudge that corrupts the transient covariance.

        History: `doc/shooting_history.md`, `PAC._vanloan_step_injection`.
        """
        import scipy.linalg as sla
        Cr = np.asarray(Cr, dtype=float)
        Gr = np.asarray(Gr, dtype=float)
        CYr = np.asarray(np.real(CYr), dtype=float)
        m = Cr.shape[0]
        d = [i for i in range(m)
             if np.any(np.abs(Cr[i, :]) > 0) or np.any(np.abs(Cr[:, i]) > 0)]
        a = [i for i in range(m) if i not in d]
        if not d:
            raise NotImplementedError(
                'PAC: this circuit has no capacitive (differential) node, so '
                'there is no covariance to propagate -- every state is '
                'algebraic and a white source on it is differentiated by the '
                'DAE. Add the capacitance that shunts the noise, or ask for a '
                'quantity that does not need a covariance.')
        di = np.ix_(d, d)
        Emb = np.zeros((m, len(d)))
        for k, i in enumerate(d):
            Emb[i, k] = 1.0
        if a:
            Gaa = Gr[np.ix_(a, a)]
            Gai = np.linalg.inv(Gaa)
            Gad = Gr[np.ix_(a, d)]
            Sc = Gr[di] - Gr[np.ix_(d, a)] @ Gai @ Gad
            ## R_proj = [I_d, -G_da G_aa^-1] routes the algebraic-row noise
            ## into the differential rows through the same elimination
            Rproj = np.zeros((len(d), m))
            for k, i in enumerate(d):
                Rproj[k, i] = 1.0
            Rproj[:, a] = -Gr[np.ix_(d, a)] @ Gai
            CYred = Rproj @ CYr @ Rproj.T
            ## the algebraic variables are slaved to the differential ones
            Emb[np.ix_(a, range(len(d)))] = -Gai @ Gad
        else:
            Sc = Gr[di]
            CYred = CYr[di]
        Cdd = Cr[di]
        Cinv = np.linalg.inv(Cdd)
        Ared = -Cinv @ Sc
        ## CY is a ONE-SIDED density; CY/2 is the two-sided intensity, the
        ## same convention `_lyapunov_pieces` and `diffusion_constant` use
        Dred = Cinv @ (0.5 * CYred) @ Cinv.T
        Dred = 0.5 * (Dred + Dred.T)
        md = len(d)
        Z = np.zeros((md, md))
        E = sla.expm(np.block([[-Ared, Dred], [Z, Ared.T]]) * float(h))
        Phi = E[md:, md:].T
        Qd = Phi @ E[:md, md:]
        Qd = 0.5 * (Qd + Qd.T)
        return Emb @ Qd @ Emb.T

    def _lyapunov_pieces_glm(self, pss, sm, what, white=None):
        """`_lyapunov_pieces` on a Nordsieck GLM's OWN map (`monodromy=
        'native'`; the default hands a GLM run's covariance to a radau twin,
        `_lyapunov_host`).

        The per-step state is the map on the state's ``(x, P)``
        (`_GLMStateStep`), width ``(r+1) m``, stacked `x` first, so `A_j` is
        that step on it, dense.  The closure stays on `x`: step 0 opens with
        a startup, which reads `x` alone, so the covariance's Nordsieck
        block at the period start never matters, `M` is the map on the state
        and `K1` the `x` block of the accumulation -- the Kronecker system
        stays ``m^2``.  The consumers' walks pad `K0` to the step width
        (`_lyap_walk`).

        ⚠ ONE SHARED SAMPLE PER STEP, AND SO FIRST ORDER.  A white source
        reaches a GLM step's state through its effective weights ``w = l^T
        B`` (GLM3 0.359, -0.0167, 0.067, 0.591; GLM4 -26 .. +166), which no
        set of independent per-stage samples with positive variances
        carries.  So the source is one constant over the step, entering
        every stage, output row and opening startup substage as a
        transient's source does (`T_j`, the step's response to it), with
        variance ``CYbar / 2h``: `CY` at the stage states averaged with
        positive weights over the abscissae (`_abscissa_weights`,
        normalised).

        History: `doc/shooting_history.md`, `PAC._lyapunov_pieces_glm`.
        """
        m = pss.cir.n - 1
        irn = pss.irefnode
        steps = sm.step_objects()
        r = len(sm.steps[0].Qin)
        na = (r + 1) * m
        w0 = 2.0 * np.pi / float(sm.T)

        def carry(Z):
            return ([Z[(k + 1) * m:(k + 2) * m] for k in range(r)], Z[:m])

        def stack(c):
            return np.vstack([np.asarray(c[1])] + [np.asarray(b) for b in c[0]])

        E = np.eye(na)
        I = np.eye(m)
        Z0 = carry(np.zeros((na, m)))
        As, Qs = [], []
        for st in steps:
            As.append(stack(st.solve(carry(E))))
            sub = (None if st.startup is None
                   else [[I] * 3 for _ in st.startup.lus])
            Tj = stack(st.solve(Z0, ([I] * st.s, sub)))
            wts = self._abscissa_weights(st.rec.c)
            wts = wts / float(np.sum(wts))
            CY = np.zeros((m, m))
            for i, y in enumerate(st.rec.Ys):
                if wts[i] > 0.0:
                    CY += wts[i] * np.real(np.asarray(self._lyap_cy(
                        pss, w0, np.delete(np.asarray(y, dtype=float), irn),
                        white),
                        dtype=complex))
            Q = Tj @ (CY / (2.0 * st.rec.h)) @ Tj.T
            Qs.append(0.5 * (Q + Q.T))
        K = np.zeros((na, na))
        for A_j, Q_j in zip(As, Qs):
            K = A_j @ K @ A_j.T + Q_j
        M = np.column_stack([np.asarray(sm.matvec(e), dtype=float)
                             for e in np.eye(m)])
        return As, Qs, K[:m, :m], M, m, m

    @staticmethod
    def _lyap_walk(As, Qs, K0):
        """The per-node covariances from `K0` at node 0: ``K_{j+1} = A_j K_j
        A_j^T + Q_j``.  A step wider than `K0` (a GLM's ``(x, P)``,
        `_lyapunov_pieces_glm`) starts from `K0` padded -- its startup reads
        `x` alone -- and each sample is the `x` block."""
        n = K0.shape[0]
        na = As[0].shape[0] if As else n
        K = K0 if na == n else np.pad(K0, ((0, na - n), (0, na - n)))
        seq = [K0]
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
            seq.append((0.5 * (K + K.T))[:n, :n])
        return seq

    def _lyapunov_pieces_stage(self, pss, fp, what, white=None):
        """`_lyapunov_pieces` for a Runge-Kutta stage method's Floquet source
        (Radau IIA, TR-BDF2, ESDIRK).

        The per-step transition `A_n` is the stage step map (dense, `m x m`,
        via `_monodromy_matvec_stage` one step at a time) and the per-step
        injection `Q_n` is the stage injection (`_stage_injection`) -- the
        source enters every STAGE; the end-of-step DAE-projected Van Loan
        integral (`_vanloan_step_injection`, first order across a switch
        edge) stays the fallback.  State width `m`, so `n = m`.

        ⚠ THE VAN LOAN INJECTION IS EXACT, THE METHOD SETS ONLY THE
        PROPAGATION.  Under it the covariance still converges to the
        stationary target (kT/C on an RC) at the injection's O(h^2), not at
        Radau's O(h^5): Van Loan already integrates the step exactly, so
        refining the grid gains on the recursion's discretisation of a
        continuous Lyapunov flow, which the higher-order transition does not
        change.

        History: `doc/shooting_history.md`, `PAC._lyapunov_pieces_stage`.
        """
        if white is None:
            self._refuse_coloured(pss, what)
        m = pss.cir.n - 1
        n = m
        hs = np.diff(np.asarray(fp.times, dtype=float))
        w0 = 2.0 * np.pi / float(fp.T)
        _W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                       pss.irefnode, axis=0)
        As, Qs = [], []
        for k, step in enumerate(fp.steps):
            xk = _W[:, min(k + 1, _W.shape[1] - 1)]
            Cn = np.asarray(pss._C_at(xk), dtype=float)
            Gn = np.asarray(pss._G_at(xk), dtype=float)
            CYn = self._lyap_cy(pss, w0, xk, white)
            A_k = np.column_stack([
                np.asarray(pss._monodromy_matvec_stage([step], e), dtype=float)
                for e in np.eye(m)])
            As.append(A_k)
            Q_k = self._stage_injection(pss, fp, k, w0, white)
            Qs.append(Q_k if Q_k is not None
                      else self._vanloan_step_injection(Cn, Gn, CYn, hs[k]))
        K = np.zeros((n, n))
        for A_k, Q_k in zip(As, Qs):
            K = A_k @ K @ A_k.T + Q_k
        M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                             for e in np.eye(n)])
        return As, Qs, K, M, m, n

    def _stage_injection(self, pss, fp, k, w, white=None):
        """The per-step process-noise covariance `Q_k` of a STAGE method with
        the source entering EVERY stage, or None when the period is not a
        stage method's (then the caller keeps its end-of-step Van Loan).

            Q_k = sum_i T_i (CY(Y_i) / (2 h b_i)) T_i^T,
            T_i = d x_{k+1} / d u_i   through the method's own stage solve

        with `CY` at the STAGE states `Y_i`.  White noise over the step is
        the method's quadrature `h sum_i b_i u_i` with independent stage
        samples of variance `CY/(2 h b_i)`, so the increment's variance is
        `h CY/2` -- the diffusion -- and each sample reaches the step's end
        through the stage equations exactly as a stage source does.

        ⚠⚠ The Van Loan injection freezes `C`, `G`, `CY` at the END of the
        step, so across a switch-off edge it integrates the injection with
        the OFF conductance and is first order (a switched capacitor's held
        variance 12 % low at 400 points under radau); the stage injection
        takes the method's own order there.
        ⚠ Needs stiff accuracy (`x_{k+1} = Y_s`) and positive weights:
        radau and trbdf2 qualify.  ⚠⚠ A TABLEAU WITH A NON-POSITIVE WEIGHT
        (ESDIRK43: `b = 0.158, 0, 0.187, 0.681, -0.275, 0.25`) takes the Van
        Loan injection at the STAGE states instead, averaged with POSITIVE
        trapezoid weights over the stage abscissae in time
        (`_abscissa_weights`).  (Equal-variance stage samples, the other
        tableau-independent candidate, were measured and are worse; see
        history.)

        History: `doc/shooting_history.md`, `PAC._stage_injection`.
        """
        m = pss.cir.n - 1
        tms = np.asarray(fp.times, dtype=float)
        h = float(tms[k + 1] - tms[k])
        if not fp.is_stage:
            return None
        st = fp.steps[k]
        bvec, cvec = st.b, st.c
        s = st.s
        states = self._stage_states(pss, fp)
        irn = pss.irefnode
        if np.any(bvec <= 0.0):
            wts = self._abscissa_weights(cvec)
            Q = np.zeros((m, m))
            for i in range(s):
                if wts[i] <= 0.0:
                    continue
                yi = np.delete(np.asarray(states[k * s + i], dtype=float), irn)
                CYi = np.real(np.asarray(self._lyap_cy(pss, w, yi, white),
                                         dtype=complex))
                Q += wts[i] * self._vanloan_step_injection(
                    np.asarray(pss._C_at(yi), dtype=float),
                    np.asarray(pss._G_at(yi), dtype=float), CYi, h)
            return 0.5 * (Q + Q.T)
        Q = np.zeros((m, m))
        for i in range(s):
            yi = np.delete(np.asarray(states[k * s + i], dtype=float), irn)
            CYi = np.real(np.asarray(self._lyap_cy(pss, w, yi, white),
                                         dtype=complex))
            Ti = st.source_response(i)
            Q += Ti @ (CYi / (2.0 * h * bvec[i])) @ Ti.T
        return 0.5 * (Q + Q.T)

    @staticmethod
    def _abscissa_weights(cvec):
        """Positive trapezoid weights over stage abscissae `c` in [0, 1]
        (summing to one), shared equally among stages with the same `c`."""
        c = np.asarray(cvec, dtype=float)
        uniq = np.unique(np.round(c, 12))
        wu = np.zeros(len(uniq))
        for i, u in enumerate(uniq):
            lo = uniq[i - 1] if i > 0 else u
            hi = uniq[i + 1] if i < len(uniq) - 1 else u
            wu[i] = 0.5 * (hi - lo)
            if i == 0:
                wu[i] += 0.5 * u
            if i == len(uniq) - 1:
                wu[i] += 0.5 * (1.0 - u)
        wts = np.zeros(len(c))
        for i, u in enumerate(uniq):
            same = [k for k in range(len(c)) if abs(c[k] - u) < 1e-12]
            for k in same:
                wts[k] = wu[i] / len(same)
        return wts

    def _event_closure(self, pss, As, Qs, M, m, n):
        """The BORDERED Lyapunov closure on a staged solve, or `None` when
        the solve is not staged.

        On a staged solve the state's linearised period map is not `M`
        alone: the per-step noise `w_j` moves the landed crossings,
        ``dtheta = -Gt^-1 (G dx_0 + sum_j d_j w_j)`` with ``d_j[k] = W_k
        P_{nd_k <- j+1}`` (the event row's sensitivity to the noise of
        step `j`), and the state at the period carries ``P_theta dtheta``.
        Stationarity then closes on the TOTAL monodromy `M + P_theta
        dtheta/dx_0` with the injection ``Q_tot = Cov(u - P_theta Gt^-1
        v)``: ``[I, -P_theta Gt^-1] Cov([u; v]) [.]^T`` from ``Cov(u) =
        K_1`` (the plain forward recursion), ``Cov(u, v) = E = Z_N`` with
        ``Z_{j+1} = A_j Z_j + Q_j d_j^T`` and ``Cov(v) = D = sum_j d_j Q_j
        d_j^T``.

        ⚠ THE UNBORDERED CLOSURE ON A STAGED SOLVE IS NOT MERELY
        INCOMPLETE, IT IS WRONG BY O(1): the landed window step's OWN
        linearisation carries the threshold noise through the switch with
        a gain the three collocation points invent (a comparator-jitter
        sampler reads 3.8x its analytic held variance).  With the crossing
        conditions pinned at both window edges the bordered system cancels
        that internal sensitivity, which is why the moving events must be
        unknowns of the noise problem too.

        SAMPLES ARE AT FIXED TIMES.  The grid's nodes move with the
        events (`_event_remap`: the steps between two crossings scale
        together), so the covariance of "node j" would include the node's
        own motion along the orbit -- ``xdot_j tau_j^T dtheta``, `tau_j =
        sum_{i<j} dh_i/dtheta` -- an artefact the size of the physics (a
        NOISELESS ramp source would read the threshold's kT/C).  The
        sample is the state at the node's unperturbed time: ``Pk_j^fixed
        = Pk_j - xdot_j tau_j^T``, ``R_j = P_{j<-0} + Pk_j^fixed dth``,
        ``Cov_j = R_j K_0 R_j^T + K_j^fwd - Z_j Gt^-T Pk_j^T - Pk_j Gt^-1
        Z_j^T + Pk_j Gt^-1 D Gt^-T Pk_j^T`` (`Pk_j` fixed throughout); at
        `j = 0` this is `K_0` and at `j = N` it closes back to `K_0`.
        ``dtheta`` depends on the noise of EVERY step, the future ones
        included -- a crossing later in the period moves the node's time
        now -- so the sample recursion is not causal step by step; the
        period-level objects are exact.

        Returns ``(M_tot, Q_tot, samples)`` with ``samples(K0)`` the list
        of per-node covariances.  Built for the one-step hosts whose
        per-step maps are the state maps (`n == m`: radau; trbdf2 borrows
        its gear twin, which has no columns -- warned, unbordered).

        History: `doc/shooting_history.md`, `PAC._event_closure`."""
        ev = getattr(pss, '_event_columns', None)
        if ev is None:
            return None
        N = len(As)
        Pk_nodes = np.asarray(ev['Pk_nodes'], dtype=float)
        P_end = np.asarray(ev['P_end'], dtype=float)
        pair = (n == 2 * m and P_end.shape[0] == 2 * m)
        if (n != m and not pair) or Pk_nodes.shape[0] != N + 1:
            warnings.warn(
                'PAC.covariance: the solve is staged on its state events, '
                'but this Floquet host (%s, %d steps for %d event-column '
                'nodes) is not the one the event columns were built on -- '
                'the closure runs UNBORDERED and its answer through the '
                'switching instants is not to be trusted. Solve with '
                'method=\'radau\' for the bordered closure.'
                % (getattr(pss.par, 'method', '?'), N, Pk_nodes.shape[0] - 1),
                RuntimeWarning, stacklevel=3)
            return None
        nodes = [int(j) for j in ev['nodes']]
        ## gear's PAIR form: the state is (x_j, x_{j-1}), the event row acts
        ## on the first block, the per-node column of node j is the pair
        ## (Pk_j, Pk_{j-1}) and the map to node j the pair of `P_nodes` rows;
        ## the samples come out as pair covariances, as the plain gear path
        ## returns them
        ## the recursion runs at the STEP's width -- `n`, or a GLM's native
        ## `(x, P)` (`_lyapunov_pieces_glm`), read back at `n`
        na = As[0].shape[0] if As else n
        W = [np.pad(np.asarray(w, dtype=float).ravel(), (0, na - m)) for w in ev['W']]
        K = len(nodes)
        P_nodes = np.asarray(ev['P_nodes'], dtype=float)
        if pair:
            Pk_prev = np.concatenate((Pk_nodes[:1] * 0.0, Pk_nodes[:-1]), axis=0)
            Pk_nodes = np.concatenate((Pk_nodes, Pk_prev), axis=1)          # (N+1, 2m, K)
            P_prev = np.concatenate((np.zeros((1,) + P_nodes.shape[1:]), P_nodes[:-1]), axis=0)
            P_prev[0] = np.hstack((np.zeros((m, m)), np.eye(m)))          # node -1 is the pair's second block
            P_nodes = np.concatenate((P_nodes, P_prev), axis=1)            # (N+1, 2m, 2m)
        Gt = np.asarray(ev['Gt'], dtype=float)
        dth = np.asarray(pss._event_sensitivity, dtype=float)
        ## d_j[k] = W_k A_{nd_k - 1} ... A_{j+1}: the event row's response to
        ## the noise landing at node j+1 (zero once the crossing is past)
        d = np.zeros((N, K, na))
        for k, nd in enumerate(nodes):
            r = W[k].copy()
            for j in range(nd - 1, -1, -1):
                d[j, k] = r
                r = r @ As[j]
        Z = np.zeros((N + 1, na, K))
        Kf = np.zeros((N + 1, na, na))
        D = np.zeros((K, K))
        for j in range(N):
            Z[j + 1] = As[j] @ Z[j] + Qs[j] @ d[j].T
            Kf[j + 1] = As[j] @ Kf[j] @ As[j].T + Qs[j]
            D = D + d[j] @ Qs[j] @ d[j].T
        Z, Kf = Z[:, :n], Kf[:, :n, :n]
        Gi = np.linalg.inv(Gt)
        E = Z[N]
        M_tot = M + P_end @ dth
        Q_tot = (Kf[N] - E @ Gi.T @ P_end.T - P_end @ Gi @ E.T
                 + P_end @ Gi @ D @ Gi.T @ P_end.T)
        Q_tot = 0.5 * (Q_tot + Q_tot.T)
        Pk_fixed, _tau, _xdot = self._fixed_time_event_columns(pss)
        if pair:
            Pkf_prev = np.concatenate((Pk_fixed[:1] * 0.0, Pk_fixed[:-1]), axis=0)
            Pk_fixed = np.concatenate((Pk_fixed, Pkf_prev), axis=1)

        def samples(K0):
            seq = []
            for j in range(N + 1):
                Pkf = Pk_fixed[j]
                Rj = P_nodes[j] + Pkf @ dth
                Cj = (Rj @ K0 @ Rj.T + Kf[j]
                      - Z[j] @ Gi.T @ Pkf.T - Pkf @ Gi @ Z[j].T
                      + Pkf @ Gi @ D @ Gi.T @ Pkf.T)
                seq.append(0.5 * (Cj + Cj.T))
            return seq
        pieces = {'dth': dth, 'Gi': Gi, 'D': D, 'nodes': nodes, 'E': E}
        return M_tot, Q_tot, samples, pieces

    def event_jitter(self, pss, fmin=None, fmax=None, points_per_decade=40):
        """The noise-driven JITTER of every landed crossing of a staged,
        driven solve: ``sigma`` in seconds per crossing, and the crossings'
        covariance in fractions of the period.

        The bordered Lyapunov closure (`_event_closure`) already carries
        it: the crossings move as ``dtheta = (dtheta/dx_0) dx_0 - Gt^-1
        sum_j d_j w_j`` -- the stationary state at the period start
        (covariance `K_0`, from the previous periods' noise) and this
        period's per-step injections, independent of each other -- so
        ``Cov(dtheta) = dth K_0 dth^T + Gt^-1 D Gt^-T`` with ``D = sum_j
        d_j Q_j d_j^T``.  On the comparator-jitter sampler
        (`_jitter_sampler`: a sawtooth of slope s_1 crossing a threshold
        node with kT/C_n of noise) the turn-off crossing's sigma is the
        analytic ``sqrt(kT/C_n) / s_1``.

        Returns ``{'sigma_t': (K,) s, 'cov_fraction': (K, K), 'fractions':
        (K,) the crossings' positions, 'nodes': (K,) their grid nodes}``.
        A COLOURED source needs the band `fmin` / `fmax` /
        `points_per_decade`, as `covariance`: the crossings' coloured motion
        is read off the same bordered forced responses (the coloured part
        of a 1/f threshold's sigma^2 is `Var_band(v_th) / s_1^2` to 1e-6).
        An oscillator's crossings diffuse without bound with its phase;
        that is `oscillator_covariance`'s object, and this refuses one.
        Every source of the circuit is in it together; a per-source
        split is the per-source `Q_j`, which `_lyapunov_pieces` does not
        keep.

        History: `doc/shooting_history.md`, `PAC.event_jitter`."""
        self._check_circuit(pss)
        if getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.event_jitter: an OSCILLATOR\'s crossings diffuse with its '
                'phase and have no stationary jitter; use '
                'oscillator_covariance() for the growth and the bounded '
                'orbital part.')
        if getattr(pss, '_event_columns', None) is None:
            raise ValueError(
                'PAC.event_jitter: the solve has no landed state events -- '
                'solve with state_events=True on a circuit that declares '
                'them (a VSwitch).')
        pss = pss._lyapunov_host()
        col = self._coloured_prepare(pss, fmin, fmax, points_per_decade,
                                     'event_jitter')
        As, Qs, K1, M, m, n = self._lyapunov_pieces(
            pss, 'covariance', white=None if col is None else col['white'])
        bordered = self._event_closure(pss, As, Qs, M, m, n)
        if bordered is None:
            raise ValueError(
                'PAC.event_jitter: this Floquet host carries no event '
                'columns (see the warning above); solve with method=\'radau\'.')
        M_tot, Q_tot, _samples, pieces = bordered
        S = np.eye(n * n) - np.kron(M_tot, M_tot)
        K0 = np.linalg.solve(S, Q_tot.reshape(-1)).reshape(n, n)
        K0 = 0.5 * (K0 + K0.T)
        dth, Gi, D = pieces['dth'], pieces['Gi'], pieces['D']
        cov = dth @ K0 @ dth.T + Gi @ D @ Gi.T
        if col is not None:
            ## the crossings' coloured motion, from the same bordered forced
            ## responses (`_forced_responses`' shifts)
            _Kc, Dc = self._coloured_covariance(pss, col, m, n, all_nodes=False)
            cov = cov + Dc
        cov = 0.5 * (cov + cov.T)
        T = float(pss.period)
        ## (`sigma_t`, the key the other jitter surfaces use; `sigma`
        ## until 2026-09-29)
        return {'sigma_t': np.sqrt(np.clip(np.diag(cov), 0.0, None)) * T,
                'cov_fraction': cov,
                'fractions': np.asarray(pss._state_event_fracs, dtype=float).copy(),
                'nodes': list(pieces['nodes'])}

    def _injection_points(self, pss, fp):
        """`(counts, states)`: how many injection points each step has and
        their states (full width), in `injection_times` order -- where a
        MODULATED source is evaluated (`_coloured_covariance`).  A stage
        method's and a GLM's are the stage points (`_stage_states`); a
        multistep step's source enters at its END, so its one point is the
        next node (as `_sampled_series` reads it)."""
        N = len(fp.steps)
        if fp.is_glm:
            return ([len(st.injection_times(0.0)) for st in fp.step_objects()],
                    self._stage_states(pss, fp))
        if fp.is_stage:
            return [st.s for st in fp.steps], self._stage_states(pss, fp)
        ## (the column the Lyapunov pieces read `CY` at)
        xs = np.asarray(pss.waveform[1], dtype=float)
        return [1] * N, [xs[:, min(j + 1, xs.shape[1] - 1)] for j in range(N)]

    def _coloured_prepare(self, pss, fmin, fmax, points_per_decade, what):
        """None on a circuit whose sources are all white.  Otherwise the
        coloured components and the band, and under `'white'` the WHITE part
        of each source, `xr -> A(x)`, for the Lyapunov pieces to read in
        place of `CY` (`_lyap_cy`).  Refuses what cannot be integrated: no `fmin` (a
        1/f variance grows as ``ln(fmax/fmin)`` without limit), a circuit
        whose `CY` is not the sum of its elements'.  A colour that is not a
        power law is taken three ways: STATIONARY (its own `CY(nu)`),
        modulated SEPARABLE (a level per point, a shape per frequency), or
        modulated with a moving shape (its density per point per frequency,
        warned)."""
        if not self._coloured_present(pss):
            return None
        fp = pss._state_map()
        T = float(fp.T)
        N = len(fp.steps)
        f0 = 1.0 / T
        fnyq = 0.5 * N / T
        if fmin is None:
            raise NotImplementedError(
                'PAC.%s: a noise source in this circuit is COLOURED (a 1/f '
                'source), and its variance grows as ln(fmax/fmin) without '
                'limit -- pass fmin (and fmax, default the grid\'s Nyquist, '
                '%.6g Hz): the coloured part is integrated over [fmin, fmax] '
                'in the frequency domain, the white part as for white '
                'sources.' % (what, fnyq))
        fmin = float(fmin)
        fmax = fnyq if fmax is None else float(fmax)
        if not (0.0 < fmin < fmax <= fnyq * (1.0 + 1e-12)):
            raise ValueError(
                'PAC.%s: need 0 < fmin < fmax <= the grid\'s Nyquist '
                '(N/2T = %.6g Hz); got fmin = %.6g, fmax = %.6g.'
                % (what, fnyq, fmin, fmax))
        counts, states = self._injection_points(pss, fp)
        nc = self._noise_components(pss, states)
        with warnings.catch_warnings():
            ## (the per-band FOLD's caveat -- one root per element -- does not
            ## apply here: a per-band component enters through its `CY`
            ## itself, below, with no square root)
            warnings.filterwarnings(
                'ignore', message='PAC: the noise of .* is not '
                'thermal-plus-power-law')
            model = nc.model(fmin, f0)
        if model is None:
            raise NotImplementedError(
                'PAC.%s: this circuit\'s CY is not the sum of its elements\' '
                '(see the warning above), and as a whole it is not thermal-'
                'plus-power-law, so its coloured part cannot be separated from '
                'the white one.' % what)
        ## ⚠ A COLOUR THAT IS NOT A POWER LAW (a Lorentzian `IS(noiseTau)`)
        ## has no density to factor out, but a STATIONARY one -- the same
        ## `CY(w)` at every point of the orbit -- needs none: its response is
        ## linear in the source, so per band frequency `K += sum_kl
        ## CY_kl(nu) Re[y_k y_l^H]` over unit sources `e_k` on its support
        ## (`_coloured_covariance`).  A MODULATED one: below.
        ## History: `doc/shooting_history.md`, `PAC._coloured_prepare`.
        perband, separable, nonseparable = [], [], []
        wref = 2.0 * np.pi * f0
        ## the band's TOP as well as its middle (`pi fmax`), so a shape that
        ## departs only above fmax/2 does not read as separable.  A probe
        ## more can only send a source to the exact path.
        ## History: `doc/shooting_history.md`, `PAC._coloured_prepare`.
        wt = sorted({2.0 * np.pi * fmin, wref, 20.0 * np.pi * f0,
                     np.pi * fmax, 2.0 * np.pi * fmax})
        white_split = []
        for key in model.perband:
            ## classified from the element's signed amplitudes where it
            ## states them, else from its `CY` (`perband_classify`)
            kind, Cs, Ws, mode = nc.perband_classify(key, wt)
            if mode == 'white':
                ## ⚠ ITS WHITE REMAINDER JOINS THE WHITE PART: the band
                ## integral covers [fmin, fmax] only, and a white source
                ## belongs to the Lyapunov path over every frequency (band-
                ## limited, it read -1.4 % / -1.8 % against a separate white
                ## source).  Its signed columns alone are read per band
                ## frequency -- right for any shape of theirs.
                white_split.append(key)
                root = nc.perband_root(key, 'signed')
                nonseparable.append((key, lambda nu, root=root: root(
                    2.0 * np.pi * nu)))
                continue
            if kind == 'stationary':
                C0 = np.asarray(Cs[wt.index(wref)][0], dtype=complex)
                supp = np.nonzero(np.any(np.abs(C0) > 0.0, axis=1))[0]
                if supp.size:
                    perband.append((key, supp))
                continue
            ## ⚠ MODULATED AND NOT A POWER LAW.  SEPARABLE -- a level that
            ## follows the state under a fixed spectral shape, ``C(x, w) =
            ## C(x, w_ref) s(w)``, the usual burst / G-R noise -- replays one
            ## amplitude per point and weights each band frequency by `s`;
            ## otherwise the element is read at EVERY point for EVERY band
            ## frequency (the quasi-static model `pnoise` and
            ## `sampled_variance` use, per band).  Either takes the element's
            ## SIGNED amplitudes where it states them, else the root of its
            ## PSD, warned: the |m| process, and one column per point.
            if Ws is None:
                warnings.warn(
                    'PAC.%s: the noise of %s is coloured, not a power law, and '
                    'modulated by the orbit; it states no signed amplitudes, so '
                    'it enters as the square root of its PSD -- the |m| '
                    'process, SIGN-BLIND where the modulation changes sign, and '
                    'independent sources inside the element merged into one '
                    '(Element.noise_amplitudes states both).'
                    % (what, '.'.join(key)), RuntimeWarning, stacklevel=3)
            if kind == 'separable':
                Cref = Cs[wt.index(wref)]
                W0 = psd_sqrt(Cref) if Ws is None else Ws[wt.index(wref)]
                jr, pi_, qi = np.unravel_index(int(np.argmax(np.abs(Cref))),
                                               Cref.shape)
                separable.append((key, W0, states[jr], (pi_, qi),
                                  complex(Cref[jr, pi_, qi])))
            else:
                root = nc.perband_root(key, mode)
                nonseparable.append((key, lambda nu, root=root: root(
                    2.0 * np.pi * nu)))
                warnings.warn(
                    'PAC.%s: the noise of %s is coloured, not a power law, '
                    'and its spectral SHAPE changes along the orbit, so it is '
                    'read at every point for every band frequency -- the '
                    'quasi-static model pnoise and sampled_variance use per '
                    'band, and costly here.'
                    % (what, '.'.join(key)), RuntimeWarning, stacklevel=3)
        warn_signed_unused(model, 'PAC.%s' % what)
        amp = getattr(model, 'amplitude', None) or {}
        comps = []
        for key, B, EF in model.flicker:
            ef = uniform_exponent(B, EF)
            split = (exponent_columns(B, EF, amp[key])
                     if ef is None and key in amp else None)
            if split is not None:
                ## the element's signed columns, grouped by their own
                ## exponents: a uniform power law each (`exponent_columns`)
                comps.extend((key, np.asarray(Wg, dtype=complex), float(efg))
                             for Wg, efg in split)
                continue
            if ef is None:
                ## ⚠ EXPONENTS THAT DIFFER BETWEEN ENTRIES.  Entries in
                ## DISJOINT index blocks, one exponent
                ## each (sources of different slope on branches that share no
                ## node), are independent components: split exactly.
                parts = self._split_by_exponent(B, EF)
                if parts is not None:
                    for Bg, efg in parts:
                        comps.append((key, np.asarray(psd_sqrt(Bg),
                                                      dtype=complex), efg))
                    continue
                ## otherwise no one amplitude to replay: the moving-shape way,
                ## the density `B (w1/w)^EF` rooted at every point per band
                ## frequency (its white part is already in `white`: the
                ## element's own `CY` would count it twice)
                nonseparable.append((key, lambda nu, B=B, EF=EF, w1=model.w1:
                                     psd_sqrt(
                                         B * (w1 / (2.0 * np.pi * nu)) ** EF)))
                warnings.warn(
                    'PAC.%s: the coloured noise of %s carries different '
                    'power-law exponents in different entries, so it has no '
                    'one amplitude to replay; its density is rooted at every '
                    'point per band frequency -- costlier, and SIGN-BLIND (a '
                    'square root per point).' % (what, '.'.join(key)),
                    RuntimeWarning, stacklevel=3)
                continue
            ## ⚠ THE SIGN: the element's stated amplitudes where it has them
            ## (`W W^H = B` with the sign of the modulation); `sqrt(B)` is
            ## the sign-blind |m| process (`warn_signed_unused` said so)
            W = amp.get(key)
            W = np.asarray(W if W is not None else psd_sqrt(B), dtype=complex)
            comps.append((key, W, float(ef)))
        ## the white part of each source, at the states the pieces read:
        ## the injection points from the batch model, any other state (a
        ## step end under the Van Loan fallback) fitted on demand
        irn = pss.irefnode
        m = pss.cir.n - 1

        def remainder(sts):
            ## the white remainders ``C - W W^H`` of the split per-band
            ## elements at `sts` (white: read at one frequency)
            out = 0.0
            for key in white_split:
                at = self._noise_components(pss, sts)
                C_ = at.one_element_cy(key, wref)
                W_ = at.one_element_amplitudes(key, wref)
                out = out + (C_ - np.einsum('kis,kjs->kij', W_, W_.conj()))
            return out
        rem = remainder(states) if white_split else None
        cache = {}
        for k, (x, A) in enumerate(zip(states, model.white)):
            xr = np.delete(np.asarray(x, dtype=float), irn)
            cache[xr.tobytes()] = np.asarray(A, dtype=complex) + (
                rem[k] if rem is not None else 0.0)

        def white(xr):
            xr = np.asarray(xr, dtype=float).ravel()[:m]
            key = xr.tobytes()
            if key not in cache:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    mdl = self._noise_components(pss, [xr]).model(fmin, f0)
                if mdl is None:
                    raise NotImplementedError(
                        'PAC.%s: the circuit\'s CY, taken as a whole, is not '
                        'thermal-plus-power-law at a state the white part is '
                        'read at, so the white part cannot be separated '
                        'there.' % what)
                cache[key] = np.asarray(mdl.white[0], dtype=complex) + (
                    remainder([xr])[0] if white_split else 0.0)
            return cache[key]
        return {'fp': fp, 'counts': counts, 'comps': comps, 'w1': model.w1,
                'white': white,
                'perband': perband, 'separable': separable,
                'nonseparable': nonseparable, 'state0': states[0],
                'states': states,
                'fmin': fmin, 'fmax': fmax, 'ppd': int(points_per_decade)}

    ## The band integral's error target and its bounds (see
    ## `_coloured_covariance`): relative to the integral, estimated per pair
    ## of intervals as `|Q_h - Q_2h|`.
    COLOURED_REFINE_TOL = 1e-6
    COLOURED_REFINE_PASSES = 16
    COLOURED_REFINE_MAXPTS = 20000

    def _coloured_covariance(self, pss, col, m, n, responses=None, lines=(),
                             all_nodes=True):
        """The COLOURED sources' covariance at every node and the crossings'
        -- the frequency-domain half of a coloured `covariance` /
        `event_jitter`.  Returns `(K (N + 1, n, n), Cov(dtheta) or None)`.

        A component is ``u(t) = W(x(t)) zeta(t)``, `W` its amplitudes at the
        injection points, `zeta` independent unit processes whose one-sided
        density ``(w1/w)^EF`` makes ``W W^H (w1/w)^EF`` the component's `CY`.
        Per column of `W` and per input frequency `nu` of a log grid over
        `[fmin, fmax]`, the steady response `y(nu, t_j)` to the MODULATED
        source ``W e^{j 2 pi nu t}`` (`_forced_responses`: bordered and at
        fixed time on a staged solve), and

            K(t_j) = int_fmin^fmax (w1 / 2 pi nu)^EF Re[y y^H] dnu

        -- the two-sided density `CY/2` over +-nu, `-nu` the conjugate.  A
        log grid, `points_per_decade` as `sampled_variance`'s, the power law
        integrated EXACTLY on each interval and the response linear in `ln
        nu` (`_power_law_weights`).  A stationary colour that is
        not a power law enters through its `CY(nu)` and unit sources on its
        support.  Gear's PAIR covariance carries `(x_j,
        x_{j-1})`: the previous node's response, node -1 being node N - 1
        a period back (``e^{-j 2 pi nu T}``).

        ⚠ A SLOPE, NOT A STATE: no shaping filter, no fitted Lorentzian
        ladder -- the exact power law over a hard band, as `sampled_variance`
        and `coloured_diffusion` take it.

        `responses` replaces `_forced_responses` (same signature): an
        oscillator's TRANSVERSE responses (`_transverse_responses`).

        ⚠ THE GRID IS ADAPTIVE.  A response with a narrow line -- a high-Q
        tank in a driven circuit, an oscillator's orbital modes at EVERY
        harmonic -- is under-resolved by a fixed log grid (a Q = 20 tank
        reads -10.5 % at 40 per decade).  Adaptive Simpson on the
        log axis: a panel (two equal intervals, the rule above) is compared
        with itself on five points, ``|S_2 - S_1| / 15``, against its share
        (by `ln`-width) of `COLOURED_REFINE_TOL`, and split where it fails.
        ⚠ Not ``|Q_h - Q_2h|``: that is the PLAIN rule's second-order error,
        and driving it below the target refines smooth regions the fourth-
        order rule already has.
        `lines` -- ``(centre, half-width)`` in Hz, an oscillator's orbital
        lines -- are resolved first, so a line narrower than the starting
        grid cannot be missed.  `all_nodes=False` keeps node 0 only (a
        caller that reads `K[0]` or the crossings).

        History: `doc/shooting_history.md`, `PAC._coloured_covariance`."""
        fp = col['fp']
        if n != m and not fp.is_pair:
            raise NotImplementedError(
                'PAC.covariance: the plain trapezoidal map\'s covariance is '
                'on the pair (x, iq), whose second block is a companion '
                'current and not a node -- the coloured part is built on the '
                'nodes. Use gear, radau or trbdf2 for a coloured covariance.')
        T = float(fp.T)
        N = len(fp.steps)
        fmin, fmax = col['fmin'], col['fmax']
        ## the power law integrated EXACTLY between the grid points, the
        ## response linear in ln(nu) (`_power_law_weights`): exact for a pure
        ## power law under a flat response, second order in the grid ratio
        ## where the response bends
        nn = max(2, int(np.ceil(col['ppd'] * np.log10(fmax / fmin)))) + 1
        nn += (nn - 1) % 2            # an even number of intervals: Richardson
        grid = np.geomspace(fmin, fmax, nn)
        offs = np.concatenate(([0], np.cumsum(col['counts'])))
        zero = np.zeros(m, dtype=complex)
        resp = self._forced_responses if responses is None else responses
        nk = N + 1 if all_nodes else 1

        def node_responses(y, nu):
            y = np.asarray(y, dtype=complex)[:N + 1]
            if n != m:
                prev = np.vstack((y[N - 1:N] * np.exp(-2j * np.pi * nu * T),
                                  y[:N]))
                y = np.hstack((y, prev))
            return y[:nk]

        ## the integrand's TERMS: per term a map from band frequencies to
        ## per-frequency contributions `(G, Dd)` (node covariances, crossings'),
        ## its power-law exponent, and the evaluated frequencies
        terms = []

        def reduced(W):
            ## ⚠ THE COMPONENT'S OWN RANK, not `m` columns: ``W_k U`` with `U`
            ## an orthonormal basis of the stacked rows is the same process
            ## (``U^T zeta`` are independent unit processes) in `rank`
            ## columns -- a two-terminal source's symmetric root has two
            ## proportional columns, three resistors' seven (one term each)
            W = np.asarray(W, dtype=complex)
            _u, sv_, vh = np.linalg.svd(W.reshape(-1, W.shape[2]),
                                        full_matrices=False)
            r = int(np.sum(sv_ > 1e-12 * max(float(sv_[0]), 1e-300))) \
                if sv_.size else 0
            return np.einsum('kms,sr->kmr', W, vh[:r].conj().T)

        def column_ev(u_points, shape=None):
            ## ONE unit column through the response: per band frequency the
            ## node covariance `y y^H` and the event term `dd dd^H`, scaled by
            ## `shape(nu)` for a separable source (evaluated first, as it was
            ## when each kind of column had its own copy of this)
            def ev(batch):
                ys, dths = resp(pss, fp, batch, zero, u_points=u_points)
                out = []
                for i, nu in enumerate(batch):
                    sc = shape(nu) if shape is not None else None
                    y = node_responses(ys[i], nu)
                    Dd = None
                    if dths[i] is not None:
                        dd = np.asarray(dths[i], dtype=complex)
                        Dd = np.real(np.outer(dd, dd.conj()))
                        if sc is not None:
                            Dd = sc * Dd
                    G = np.real(np.einsum('ji,jk->jik', y, y.conj()))
                    out.append((G if sc is None else sc * G, Dd))
                return out
            return ev

        for _key, W, ef in col['comps']:
            W = reduced(W)
            for s_ in range(W.shape[2]):
                Wc = W[:, :, s_]
                if not np.any(Wc):
                    continue
                u_points = [Wc[offs[j]:offs[j + 1]] for j in range(N)]
                terms.append((float(ef), column_ev(u_points), {}))
        ## a STATIONARY colour that is not a power law: unit sources on its
        ## support, weighted per band frequency by its own `CY(nu)`
        for key, supp in col.get('perband', ()):
            def ev(batch, key=key, supp=supp):
                per_k = []
                for k in supp:
                    e = np.zeros(m, dtype=complex)
                    e[k] = 1.0
                    per_k.append(resp(pss, fp, batch, e))
                out = []
                for i, nu in enumerate(batch):
                    ## (the ONE element: `element_cy_samples` evaluated
                    ## every element to keep this one, per band frequency)
                    cy = np.asarray(self._noise_components(
                        pss, [col['state0']]).one_element_cy(
                            key, 2.0 * np.pi * nu)[0],
                        dtype=complex)[np.ix_(supp, supp)]
                    Y = np.array([node_responses(pk[0][i], nu) for pk in per_k])
                    Dd = None
                    if per_k[0][1][i] is not None:
                        dd = np.array([np.asarray(pk[1][i], dtype=complex)
                                       for pk in per_k])
                        Dd = np.real(np.einsum('ka,kl,lb->ab', dd, cy, dd.conj()))
                    out.append((np.real(np.einsum('kja,kl,ljb->jab', Y, cy,
                                                  Y.conj())), Dd))
                return out
            terms.append((0.0, ev, {}))

        ## MODULATED, SEPARABLE: one amplitude per point, the spectral shape
        ## `s(nu)` (the element at one point, its dominant entry) per band
        ## frequency
        for key, W0, xref, pq, cref in col.get('separable', ()):
            Wr = reduced(W0)

            def shape(nu, key=key, xref=xref, pq=pq, cref=cref):
                c = self._noise_components(pss, [xref]).one_element_cy(
                    key, 2.0 * np.pi * nu)[0]
                return float(np.real(c[pq] / cref))
            for s_ in range(Wr.shape[2]):
                u_points = [Wr[offs[j]:offs[j + 1], :, s_] for j in range(N)]
                terms.append((0.0, column_ev(u_points, shape), {}))
        ## MODULATED, THE SHAPE MOVING: per band frequency the columns at
        ## every point, one replay per column (`root_at(nu)`: the element's
        ## signed amplitudes or the root of its density, or the root of a
        ## power law whose exponent differs between entries)
        for _key, root_at in col.get('nonseparable', ()):
            def ev(batch, root_at=root_at):
                out = []
                for nu in batch:
                    Wn = reduced(root_at(nu))
                    G = np.zeros((nk, n, n))
                    Dd = None
                    for s_ in range(Wn.shape[2]):
                        u_points = [Wn[offs[j]:offs[j + 1], :, s_]
                                    for j in range(N)]
                        Gs, Di = column_ev(u_points)(np.array([nu]))[0]
                        G += Gs
                        if Di is not None:
                            Dd = Di if Dd is None else Dd + Di
                    out.append((G, Dd))
                return out
            terms.append((0.0, ev, {}))

        def evaluate(g):
            for _ef, ev, store in terms:
                new = [float(x) for x in g if float(x) not in store]
                if new:
                    for x, r in zip(new, ev(np.asarray(new))):
                        store[x] = r

        fac = lambda ef: (col['w1'] / (2.0 * np.pi)) ** ef   # noqa: E731

        def summary(r):
            G, Dd = r
            v = float(np.sum(np.trace(G, axis1=1, axis2=2)))
            return v + (float(np.trace(Dd)) if Dd is not None else 0.0)

        ## PANELS `(a, c)` (their log midpoint implied): the known lines
        ## first -- a panel holding a line centre is split until it is
        ## narrower (in ln nu) than half the line's own width
        mid = lambda a, c: float(np.sqrt(a * c))              # noqa: E731
        panels = [(float(grid[k]), float(grid[k + 2]))
                  for k in range(0, grid.size - 1, 2)]
        for _p in range(self.COLOURED_REFINE_PASSES):
            split = False
            nxt = []
            for a, c in panels:
                if any(a <= c0 <= c and np.log(c / a) > 0.5 * hw / c0
                       for c0, hw in lines):
                    b = mid(a, c)
                    nxt += [(a, b), (b, c)]
                    split = True
                else:
                    nxt.append((a, c))
            panels = nxt
            if not split or 2 * len(panels) > self.COLOURED_REFINE_MAXPTS:
                break
        evaluate(sorted({x for a, c in panels for x in (a, mid(a, c), c)}))

        def panel_value(pts):
            ## the rule on `pts` (3 or 5 points), all terms, as a scalar
            v = 0.0
            for ef, _ev, store in terms:
                sv = np.array([summary(store[float(x)]) for x in pts])
                v += fac(ef) * float(self._power_law_weights(pts, ef) @ sv)
            return v

        total = sum(panel_value(np.array([a, mid(a, c), c])) for a, c in panels)
        span = float(np.log(fmax / fmin))
        budget = self.COLOURED_REFINE_TOL * abs(total)
        accepted, active = [], panels
        for _p in range(self.COLOURED_REFINE_PASSES + 1):
            if not active:
                break
            five = [(a, mid(a, mid(a, c)), mid(a, c), mid(mid(a, c), c), c)
                    for a, c in active]
            if _p == self.COLOURED_REFINE_PASSES or \
                    2 * len(set().union(*[set(f) for f in five])) \
                    > self.COLOURED_REFINE_MAXPTS:
                accepted += five
                warnings.warn(
                    'PAC: the coloured band integral stopped with %d panels '
                    'short of its target error %.0e: a response line narrower '
                    'than the refinement can resolve, or a band too wide -- '
                    'raise points_per_decade.' % (len(active),
                                                  self.COLOURED_REFINE_TOL),
                    RuntimeWarning, stacklevel=3)
                break
            evaluate(sorted({x for f in five for x in f}))
            nxt = []
            for (a, c), f in zip(active, five):
                s1 = panel_value(np.array([f[0], f[2], f[4]]))
                s2 = panel_value(np.array(f))
                share = budget * float(np.log(c / a)) / span
                if abs(s2 - s1) / 15.0 <= share:
                    accepted.append(f)
                else:
                    nxt += [(f[0], f[2]), (f[2], f[4])]
            active = nxt
        ## the rule on every accepted panel (its five points: two pairs),
        ## summed per point
        K = np.zeros((nk, n, n))
        D = None
        for ef, _ev, store in terms:
            wts = {}
            for f in accepted:
                q = self._power_law_weights(np.array(f), ef) * fac(ef)
                for x, qx in zip(f, q):
                    wts[float(x)] = wts.get(float(x), 0.0) + qx
            for x, qx in wts.items():
                G, Dd = store[x]
                K += qx * G
                if Dd is not None:
                    D = qx * Dd if D is None else D + qx * Dd
        K = 0.5 * (K + np.swapaxes(K, 1, 2))
        return K, D

    def covariance(self, pss, samples=False, fmin=None, fmax=None,
                   points_per_decade=40, pair=False):
        """The periodic (cyclostationary) state covariance — DRIVEN circuits.

        ⚠ A GRID CHOSEN FOR `kT/C` IS NOT A GRID FOR THE PROFILE.  The
        injection is piecewise constant, so the covariance converges to the
        exact continuous answer at FIRST order in both phases of a switched
        circuit (gated against a closed-form time-varying reference).  On a
        `kT/C` switched capacitor the HELD value converges faster only
        because the exact profile is a constant there; the TRACKING phase
        sits at the O(h/tau) floor (4 % at 800 points).  ⚠ On a fixture
        whose noise is tied to its own conductance (`g V` and
        `white_noise(4 kT g)` with the SAME `g`), fluctuation-dissipation
        makes `V(t) = kT/C` exact at every instant for ANY `g(t)`: a
        tracking value below `kT/C` there is the discretisation floor, not
        physics, and two tools agreeing on it agree about a SHARED
        artefact.  A real device need not balance (a PSP switch's
        `sid/(4kT g)` runs 1.09 to 3.17), so its tracking limit need not be
        `kT/C`.

        Returns `K0`, the covariance of the circuit state at `t = 0`,
        `m x m`; with `samples=True`, `(K0, [K_j])`, the covariance at every
        step, which is the time-varying statistic this exists to produce.
        ⚠ `m x m` WHATEVER THE METHOD (2026-09-29): a two-step method's map
        carries `(x_n, x_{n-1})`, and its covariance was returned on that
        pair, `2m x 2m`, the shape depending on the integrator; `pair=True`
        returns it so.

        ⚠ A COLOURED SOURCE NEEDS A BAND.  With a 1/f source
        (a `flicker_noise`, a MOS channel's flicker) pass `fmin` -- and
        `fmax`, default the grid's Nyquist `N/2T` -- because a 1/f variance
        grows as `ln(fmax/fmin)` without limit.  The WHITE part of every
        source then goes through the recursion below and the COLOURED part
        is integrated over the band in the frequency domain, per input
        frequency the forced response to the modulated source
        (`_coloured_covariance`; `points_per_decade` as
        `sampled_variance`'s).  Measured: against the closed form on an RC
        to each method's transfer error (radau 4e-9, gear 5.5e-5 at 100
        points, second order); against `sampled_variance`'s adjoint route
        on a switched sampler to that route's own quadrature (5e-6).
        Refused: trap's plain map (its covariance is on the (x, iq) pair).

        The noise covariance obeys a Lyapunov recursion alongside the
        trajectory, `K_{j+1} = A_j K_j A_jᵀ + Q_j`, so over one period
        `K_N = M K_0 Mᵀ + K_1`.  Periodicity closes it:

            (I - M ⊗ M) vec(K_0) = vec(K_1)

        ⚠ ONE LINEAR SOLVE, NO NEWTON.  The Lyapunov equation is LINEAR in
        `K`, so shooting on it is exact in a single step — unlike the
        trajectory it rides on.  The monodromy of the covariance system is
        the KRONECKER SQUARE of the circuit's, so its multipliers are the
        pairwise products `lambda_i lambda_j`.

        ⚠ AND THAT IS WHY IT REFUSES AN OSCILLATOR.  There `lambda_1 = 1`
        gives `lambda_1^2 = 1`, so `I - M ⊗ M` is exactly as singular as
        `I - M`, and the covariance does not settle, it GROWS.  Variance
        linear in `t` is a random walk, which is phase diffusion, which is
        the linewidth.  Demir 2002: an oscillator's output noise is
        STATIONARY, not cyclostationary, because "noisy autonomous systems
        cannot provide a perfect time reference".  `oscillator_covariance`
        and `oscillator_spectrum` are the routes there.

        ⚠ `CY/2` IS THE ONE-SIDED-TO-TWO-SIDED CONVERSION AND IT IS NOT
        COSMETIC.  `CY` is a one-sided density (a resistor's `4kT/R`), so
        the per-step injection is `Q_j = Jf_j^-1 (CY_j / 2h_j) Jf_j^-T`.
        Against `kT/C` on an RC the full `CY` converges to 2 and the halved
        one to 1, at first order (a piecewise-constant approximation to
        white noise).

        ⚠ AND THE GRID MUST RESOLVE THE NOISE BANDWIDTH, which is a real
        precondition rather than an accuracy note: with the RC pole above
        the grid's Nyquist the discrete system does not carry the noise the
        continuous one does.  A `kT/C` that comes back low is the grid, not
        the code.

        ⚠ COST: the solve has `(2m)^2` unknowns and is dense here, so it is
        `O(m^4)`.  Small circuits only until that is replaced.

        ⚠ ON A STAGED SOLVE (`state_events=True`) THE CLOSURE IS BORDERED:
        the noise moves the landed crossings, and the plain closure on such
        a solve is wrong by O(1), not merely incomplete -- see
        `_event_closure`.  Samples are the covariance at FIXED times.

        History: `doc/shooting_history.md`, `PAC.covariance`.
"""
        self._check_circuit(pss)
        if getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.covariance: an OSCILLATOR has no periodic covariance. '
                'Its unit multiplier squares to one, so I - M kron M is '
                'singular and the covariance grows without bound rather '
                'than settling -- that growth IS the phase diffusion, and '
                'its output noise is stationary rather than '
                'cyclostationary. Use oscillator_covariance() for the '
                'split into a bounded orbital part and that growth, or '
                'oscillator_spectrum() for the lineshape it produces.')
        ## the source-injection surfaces use the gear twin when the Floquet
        ## source is TR-BDF2 (its two-stage Q_j is not built) -- see
        ## `_lyapunov_host`
        pss = pss._lyapunov_host()
        ## a COLOURED source: the white part through the recursion below,
        ## the coloured part in the frequency domain over `[fmin, fmax]`
        ## (`_coloured_covariance`)
        col = self._coloured_prepare(pss, fmin, fmax, points_per_decade,
                                     'covariance')
        As, Qs, K1, M, m, n = self._lyapunov_pieces(
            pss, 'covariance', white=None if col is None else col['white'])
        ## a staged solve closes on the TOTAL monodromy with the events'
        ## noise-driven motion in the injection -- see `_event_closure`
        bordered = self._event_closure(pss, As, Qs, M, m, n)
        if bordered is not None:
            M, K1, _samples, _pieces = bordered
        S = np.eye(n * n) - np.kron(M, M)
        K0 = np.linalg.solve(S, K1.reshape(-1)).reshape(n, n)
        K0 = 0.5 * (K0 + K0.T)
        seq = None
        if samples:
            seq = (_samples(K0) if bordered is not None
                   else self._lyap_walk(As, Qs, K0))
        if col is not None:
            Kc, _dth = self._coloured_covariance(pss, col, m, n,
                                                 all_nodes=bool(samples))
            K0 = K0 + Kc[0]
            if seq is not None:
                seq = [a + b for a, b in zip(seq, Kc)]
        if not pair:
            K0 = K0[:m, :m]
            if seq is not None:
                seq = [K[:m, :m] for K in seq]
        if not samples:
            return K0
        return K0, seq
