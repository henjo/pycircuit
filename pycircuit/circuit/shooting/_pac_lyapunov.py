"""The periodic covariance by a Lyapunov walk over the period: the per-step
pieces, the coloured band integral, event jitter.
"""
import numpy as np
import warnings
from ._factored import dense_map
from ._steps import dense_c
from .events import EventColumns
from pycircuit.circuit.simwarnings import AccuracyWarning, CostWarning, warn


class _LyapunovCovariance(object):
    """The periodic covariance by a Lyapunov walk over the period: the per-
    step pieces, the coloured band integral, event jitter.  A theme of `PAC`
    (see `pac.py`)."""

    ## The dense Lyapunov solves form ``I - M (x) M``, `n^2 x n^2`: 8 n^4 bytes
    ## before the solve's own copies (n = 200, a gear pair of m = 100: 12.8
    ## GB).  Above this width they refuse by name instead of being killed by
    ## the operating system (the review's X4, 2026-10-01); raise it on the
    ## PAC to go further.
    LYAPUNOV_DENSE_LIMIT = 110

    def _check_kron(self, n, what):
        if n > self.LYAPUNOV_DENSE_LIMIT:
            raise MemoryError(
                f'PAC.{what}: the dense Lyapunov solve forms an n^2 x n^2 '
                f'operator (n = {n}: {8.0 * n ** 4 / 1e9:.1f} GB before the '
                f'solve\'s copies), above LYAPUNOV_DENSE_LIMIT = '
                f'{self.LYAPUNOV_DENSE_LIMIT}; set pac.LYAPUNOV_DENSE_LIMIT '
                'higher if the memory is there.')

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
        prevC = [C_open] + [dense_c(st[1]) for st in fp.steps[:-1]]
        bs = {bool(st[3]) for st in fp.steps}
        if len(bs) != 1:
            raise NotImplementedError(
                'PAC.%s: the plain period mixes b = 0 and b != 0 steps, '
                'which have different per-step states.' % what)
        pair = bs.pop()
        n = 2 * m if pair else m
        As, Qs = [], []
        for k, (lu, C_new, alphas, b) in enumerate(fp.steps):
            Ck = dense_c(C_new)
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
            M = dense_map(fp, n)
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
            A_k = np.column_stack([
                np.asarray(pss._monodromy_matvec_stage([step], e), dtype=float)
                for e in np.eye(m)])
            As.append(A_k)
            Q_k = self._stage_injection(pss, fp, k, w0, white)
            if Q_k is None:
                ## the end-of-step Van Loan fallback alone reads the
                ## devices at the node (until 2026-09-30 every step did,
                ## ~8.7 s of a compact-model call, for a fallback radau
                ## and trbdf2 never take)
                xk = _W[:, min(k + 1, _W.shape[1] - 1)]
                Cn = np.asarray(pss._C_at(xk), dtype=float)
                Gn = np.asarray(pss._G_at(xk), dtype=float)
                CYn = self._lyap_cy(pss, w0, xk, white)
                Q_k = self._vanloan_step_injection(Cn, Gn, CYn, hs[k])
            Qs.append(Q_k)
        K = np.zeros((n, n))
        for A_k, Q_k in zip(As, Qs):
            K = A_k @ K @ A_k.T + Q_k
        M = dense_map(fp, n)
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

        Returns ``(M_tot, Q_tot, samples, pieces)`` with ``samples(K0)``
        the list of per-node covariances and ``pieces`` the closure's
        parts (`dth`, `Gi`, `D`, `nodes`) `event_jitter` reads, and
        `map_to(j)` / `period_map_from(j)`: the TOTAL fixed-time maps
        ``R_j`` (node 0 -> t_j) and ``(R_j, S_j)`` with ``S_j`` (t_j -> T)
        -- the period map from node j is ``M_j = R_j S_j``, its powers
        ``R_j M_tot^{k-1} S_j`` (no inverse: the step maps of a DAE are
        singular).  ``S_j`` re-borders the crossings AFTER node j on the
        state at node j at its FIXED time: ``G^(j)[k] = W_k P_{nd_k<-j}``,
        ``Gt^(j)[k, l] = W_k (Pk_{nd_k}[:, l] - P_{nd_k<-j} Pk_j^fixed[:, l])``,
        ``S_j = P_{N<-j} - P_end^(j) Gt^(j)^-1 G^(j)``; ``S_0 = M_tot``,
        ``S_N = I``, and an event AT node j counts as past (in ``R_j``).
        The crossings BEFORE node j are HELD: the steps between two
        crossings scale together, so the grid after node j still moves
        with the last one -- in the continuum that coupling cancels
        (fixed-time states propagate by the plain maps), on the grid it
        leaves ``S_j R_j - M_tot`` nonzero inside the events' span: 4.4e-7
        / 1.3e-7 of `|M_tot|` at 200 / 400 radau points on the comparator
        oscillator, 1e-14 with the held term added back, 1e-13 outside
        the span.

        Built for a host whose per-step maps are the state maps (`n ==
        m`) and for gear's PAIR form (`n == 2m`, below); any other host,
        or event columns built on another grid, runs UNBORDERED, warned.

        History: `doc/shooting_history.md`, `PAC._event_closure`."""
        ev = EventColumns.of(pss)
        if ev is None:
            return None
        N = len(As)
        Pk_nodes = np.asarray(ev['Pk_nodes'], dtype=float)
        P_end = np.asarray(ev['P_end'], dtype=float)
        pair = (n == 2 * m and P_end.shape[0] == 2 * m)
        if (n != m and not pair) or Pk_nodes.shape[0] != N + 1:
            warn(
                'PAC.covariance: the solve is staged on its state events, '
                'but this Floquet host (%s, %d steps for %d event-column '
                'nodes) is not the one the event columns were built on -- '
                'the closure runs UNBORDERED and its answer through the '
                'switching instants is not to be trusted. Solve with '
                'method=\'radau\' for the bordered closure.'
                % (getattr(pss.par, 'method', '?'), N, Pk_nodes.shape[0] - 1), AccuracyWarning)
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
        def map_to(j):
            return P_nodes[j] + Pk_fixed[j] @ dth

        def period_map_from(j):
            Rj = map_to(j)
            Pkf = Pk_fixed[j]
            ## the plain maps from node j: P_{i <- j}, i = j .. N
            Pf = [np.eye(n)]
            for i in range(j, N):
                Pf.append(np.asarray(As[i], dtype=float)[:n, :n] @ Pf[-1])
            F = [k for k, nd in enumerate(nodes) if nd > j]
            if not F:
                return Rj, Pf[-1]
            Gj = np.array([W[k][:n] @ Pf[nodes[k] - j] for k in F])
            Gtj = np.array([[W[k][:n] @ (Pk_nodes[nodes[k]][:, l]
                                         - Pf[nodes[k] - j] @ Pkf[:, l])
                             for l in F] for k in F])
            Pendj = P_end[:, F] - Pf[-1] @ Pkf[:, F]
            return Rj, Pf[-1] - Pendj @ np.linalg.solve(Gtj, Gj)

        pieces = {'dth': dth, 'Gi': Gi, 'D': D, 'nodes': nodes, 'E': E,
                  'map_to': map_to, 'period_map_from': period_map_from}
        return M_tot, Q_tot, samples, pieces

    def event_jitter(self, pss, colour_fmin=None, colour_fmax=None,
                     points_per_decade=40):
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
        A COLOURED source needs the band `colour_fmin` / `colour_fmax` /
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
        if EventColumns.of(pss) is None:
            raise ValueError(
                'PAC.event_jitter: the solve has no landed state events -- '
                'solve with state_events=True on a circuit that declares '
                'them (a VSwitch).')
        pss = pss._lyapunov_host()
        fmin, fmax = colour_fmin, colour_fmax     # (the names the internals use)
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
        self._check_kron(n, 'event_jitter')
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
        ## ⚠ A TypeError, a REQUIRED ARGUMENT MISSING, as the sampled family
        ## raises for its `series_fmin` (Python's own); NotImplementedError
        ## until 2026-09-29, though nothing here is unimplemented
        if fmin is None:
            raise TypeError(
                'PAC.%s: a noise source in this circuit is COLOURED (a 1/f '
                'source), and its variance grows as '
                'ln(colour_fmax/colour_fmin) without limit -- pass '
                'colour_fmin (and colour_fmax, default the grid\'s Nyquist, '
                '%.6g Hz): the coloured part is integrated over '
                '[colour_fmin, colour_fmax] in the frequency domain, the '
                'white part as for white sources.' % (what, fnyq))
        fmin = float(fmin)
        fmax = fnyq if fmax is None else float(fmax)
        if not (0.0 < fmin < fmax <= fnyq * (1.0 + 1e-12)):
            raise ValueError(
                'PAC.%s: need 0 < colour_fmin < colour_fmax <= the grid\'s '
                'Nyquist (N/2T = %.6g Hz); got colour_fmin = %.6g, '
                'colour_fmax = %.6g.'
                % (what, fnyq, fmin, fmax))
        counts, states = self._injection_points(pss, fp)
        nc = self._noise_components(pss, states)
        with warnings.catch_warnings():
            ## (the per-band FOLD's caveat -- one root per element -- does not
            ## apply here: a per-band component enters through its `CY`
            ## itself, below, with no square root)
            ## (`model`'s one cost note, by its category -- a text filter
            ## until 2026-10-01)
            warnings.simplefilter('ignore', CostWarning)
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
        ## the components as EVERY surface groups them (`colour_components`,
        ## review O2, 2026-10-01): a power law through fixed columns (`comps`),
        ## the rest per band frequency -- a STATIONARY per-band element
        ## through its own `CY` (no root), a SEPARABLE one replaying one
        ## amplitude per point weighted by its level, otherwise the element
        ## read at EVERY point for EVERY band frequency (the quasi-static
        ## model `pnoise` and `sampled_variance` use, per band).  The
        ## sign-blind verdict and the signed-amplitude notes are given there;
        ## an element's white remainder is already in `model.white`.
        ## History: `doc/shooting_history.md`, `PAC._coloured_prepare`.
        perband, separable, nonseparable = [], [], []
        groups = nc.colour_components(model, 2.0 * np.pi * fmin,
                                      2.0 * np.pi * fmax, f0, what)
        comps = [(key, W, ef) for key, W, ef, _B in groups.fixed]
        for band in groups.bands:
            if band.kind == 'stationary':
                if band.supp.size:
                    perband.append((band.key, band.supp))
            elif band.kind == 'separable':
                separable.append((band.key, band.W0, band.xref, band.pq,
                                  band.cref))
            else:
                nonseparable.append((band.key, lambda nu, root=band.exact_root:
                                     root(2.0 * np.pi * nu)))
                if band.kind == 'moving':
                    warn(
                        'PAC.%s: the noise of %s is coloured, not a power law, '
                        'and its spectral SHAPE changes along the orbit, so it '
                        'is read at every point for every band frequency -- '
                        'the quasi-static model pnoise and sampled_variance '
                        'use per band, and costly here.'
                        % (what, '.'.join(band.key)), CostWarning)
        ## the white part of each source, at the states the pieces read:
        ## the injection points from the batch model, any other state (a
        ## step end under the Van Loan fallback) fitted on demand
        irn = pss.irefnode
        m = pss.cir.n - 1

        cache = {}
        for x, A in zip(states, model.white):
            xr = np.delete(np.asarray(x, dtype=float), irn)
            cache[xr.tobytes()] = np.asarray(A, dtype=complex)

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
                cache[key] = np.asarray(mdl.white[0], dtype=complex)
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
                             all_nodes=True, map_node0=False):
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
        caller that reads `K[0]` or the crossings).  `map_node0`: the
        responses are the MAP's state at node 0, `n` wide, taken as they
        are (`orbital_mode_weights`: no pair stacking of node responses).

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
            if map_node0:
                return np.asarray(y, dtype=complex)[:1]
            y = np.asarray(y, dtype=complex)[:N + 1]
            if n != m:
                prev = np.vstack((y[N - 1:N] * np.exp(-2j * np.pi * nu * T),
                                  y[:N]))
                y = np.hstack((y, prev))
            return y[:nk]

        ## the integrand's TERMS: per term a map from band frequencies to
        ## per-frequency contributions `(Y, Cw, sc, Dd)` -- the node
        ## covariances as their FACTOR, ``G = sc Re[sum_kl Y_k Cw_kl Y_l^H]``
        ## per node (`Cw` None: the identity), and the crossings' `Dd` --
        ## its power-law exponent, and the evaluated frequencies.  ⚠ The
        ## factor, not `G`: `(N + 1, n, n)` per term per band frequency was
        ## the integral's memory (the review's M4, until 2026-10-01); `G` is
        ## formed once per frequency where it is summed (`node_cov`).
        terms = []

        def node_cov(r):
            ## ``sc Re[sum_kl Y_k Cw_kl Y_l^H]`` at every node
            Y, Cw, sc, _Dd = r
            if Cw is None:
                G = np.real(np.einsum('kja,kjb->jab', Y, Y.conj()))
            else:
                G = np.real(np.einsum('kja,kl,ljb->jab', Y, Cw, Y.conj()))
            return G if sc is None else sc * G

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
                    ## (a COPY: `y` is a view into the batch's responses,
                    ## which it would otherwise keep alive -- every node, every
                    ## frequency of the batch, the full width)
                    out.append((y[None].copy(), None, sc, Dd))
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
                    out.append((Y, cy, None, Dd))
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
                    Ys = []
                    Dd = None
                    for s_ in range(Wn.shape[2]):
                        u_points = [Wn[offs[j]:offs[j + 1], :, s_]
                                    for j in range(N)]
                        Yi, _c, _s, Di = column_ev(u_points)(np.array([nu]))[0]
                        Ys.append(Yi[0])
                        if Di is not None:
                            Dd = Di if Dd is None else Dd + Di
                    Y = (np.asarray(Ys, dtype=complex) if Ys
                         else np.zeros((0, nk, n), dtype=complex))
                    out.append((Y, None, None, Dd))
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
            ## the trace of `G` summed over the nodes, from the factor
            Y, Cw, sc, Dd = r
            gram = np.einsum('kja,lja->kl', Y, Y.conj())
            v = float(np.real(np.trace(gram) if Cw is None
                              else np.sum(Cw * gram)))
            v = v if sc is None else sc * v
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
                warn(
                    'PAC: the coloured band integral stopped with %d panels '
                    'short of its target error %.0e: a response line narrower '
                    'than the refinement can resolve, or a band too wide -- '
                    'raise points_per_decade.' % (len(active),
                                                  self.COLOURED_REFINE_TOL), AccuracyWarning)
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
                Dd = store[x][3]
                K += qx * node_cov(store[x])
                if Dd is not None:
                    D = qx * Dd if D is None else D + qx * Dd
        K = 0.5 * (K + np.swapaxes(K, 1, 2))
        return K, D

    def covariance(self, pss, samples=False, colour_fmin=None,
                   colour_fmax=None, points_per_decade=40, pair=False):
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

        Returns `(K0, info)`: `K0` the covariance of the circuit state at
        `t = 0`, `m x m`; with `samples=True` `info['samples']`, the
        covariance at every node -- the time-varying statistic this exists
        to produce -- at `info['times']`; with a colour band the coloured
        part alone in `info['K_coloured']` (and `info['coloured_samples']`).
        ⚠ ONE SHAPE (2026-09-29), the family's `(value, info)`: it returned
        `K0`, or `(K0, [K_j])` with `samples=True`.
        ⚠ `m x m` WHATEVER THE METHOD (2026-09-29): a two-step method's map
        carries `(x_n, x_{n-1})`, and its covariance was returned on that
        pair, `2m x 2m`, the shape depending on the integrator; `pair=True`
        returns it so.

        ⚠ A COLOURED SOURCE NEEDS A BAND.  With a 1/f source (a
        `flicker_noise`, a MOS channel's flicker) pass `colour_fmin` -- and
        `colour_fmax`, default the grid's Nyquist `N/2T` -- because a 1/f
        variance grows as `ln(colour_fmax/colour_fmin)` without limit.  The
        WHITE part of every source then goes through the recursion below and the COLOURED part
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
        ## the host the Lyapunov surfaces read: the monodromy twin (a
        ## gear/trbdf2 run is its own; TR-BDF2's injection is built) -- see
        ## `_lyapunov_host`
        pss = pss._lyapunov_host()
        ## a COLOURED source: the white part through the recursion below,
        ## the coloured part in the frequency domain over `[fmin, fmax]`
        ## (`_coloured_covariance`)
        fmin, fmax = colour_fmin, colour_fmax     # (the names the internals use)
        col = self._coloured_prepare(pss, fmin, fmax, points_per_decade,
                                     'covariance')
        As, Qs, K1, M, m, n = self._lyapunov_pieces(
            pss, 'covariance', white=None if col is None else col['white'])
        ## a staged solve closes on the TOTAL monodromy with the events'
        ## noise-driven motion in the injection -- see `_event_closure`
        bordered = self._event_closure(pss, As, Qs, M, m, n)
        if bordered is not None:
            M, K1, _samples, _pieces = bordered
        self._check_kron(n, 'covariance')
        S = np.eye(n * n) - np.kron(M, M)
        K0 = np.linalg.solve(S, K1.reshape(-1)).reshape(n, n)
        K0 = 0.5 * (K0 + K0.T)
        seq = None
        if samples:
            seq = (_samples(K0) if bordered is not None
                   else self._lyap_walk(As, Qs, K0))
        Kc = None
        if col is not None:
            Kc, _dth = self._coloured_covariance(pss, col, m, n,
                                                 all_nodes=bool(samples))
            K0 = K0 + Kc[0]
            if seq is not None:
                seq = [a + b for a, b in zip(seq, Kc)]
        cut = (lambda K: K) if pair else (lambda K: K[:m, :m])
        K0 = cut(K0)
        info = {}
        if seq is not None:
            info['samples'] = [cut(K) for K in seq]
            info['times'] = np.asarray(pss.factored_period().times,
                                       dtype=float)
        if Kc is not None:
            info['K_coloured'] = cut(Kc[0])
            if samples:
                info['coloured_samples'] = [cut(K) for K in Kc]
        return K0, info
