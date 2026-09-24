"""One step of a period map as an object: `_StageStep` (Runge-Kutta stage
methods), `_LMMStep` (linear multistep, plain and gear's pair), `_GLMStep`
(Nordsieck GLMs), and the recursions they share.
"""
import numpy as np
from ._numerics import _complex_solve
from ._numerics import _complex_solve_transposed
from ._numerics import _lu_solve_split


class _StageStep(object):
    """One step of a Runge-Kutta STAGE method's period map, kept factored --
    an entry of a 'full' or 'dirk' `FactoredPeriod.steps`, and what
    `PSS._walk_stage` solves every column of a step with.

    The stage system is ``J Z = B`` with ``J[i][j] = delta_ij C(Y_i) + h
    A_ij G(Y_j)`` -- the stage residuals ``F_i = q(Y_i) - q(x_n) - h sum_j
    A_ij K_j`` linearised, ``K_j = -(i(Y_j) + u)`` -- and the step's result is
    the LAST stage by stiff accuracy (``x_{n+1} == Y_s``).  The tableau's
    structure picks the factorisation:

    * FULLY IMPLICIT (Radau IIA): the coupled `s m x s m` block, `lu`.  The
      stages are solved together, so the step map is not a product of
      per-stage solves.
    * LOWER TRIANGULAR (TR-BDF2, ESDIRK): `J` is block lower-triangular, so
      one factor per stage, ``Kf[i] = LU(C(Y_i) + h a_ii G(Y_i))`` -- None
      for the EXPLICIT first stage, which carries its input straight through
      (``Y_0 = x_n``, so ``D_0 = carry``; `C` may be singular, so it is never
      solved for).  The coupled block would be singular there on a DAE
      (``J[0][0] = C``), which is why the two are not one factorisation.

    ⚠ ONE COPY OF EVERY REPLAY.  The two structures differ only inside
    `solve` and `adjoint`; everything else -- the forward and transposed
    mat-vecs, the forced replay and its transpose, the sideband fold, the
    stage pass, the stage-source response and the Lyapunov pieces -- reads
    the stage costates those return and the step's OWN tableau (`A`, `b`,
    `c`; never `par.method`'s, which is wrong for a `factored_period_stage`
    built with another `method=`).

    History: `doc/shooting_history.md`, `_StageStep`.
    """

    __slots__ = ('Cn', 'Gs', 'h', 'A', 'b', 'c', 'lu', 'Kf', 'm', 's')

    def __init__(self, Cn, Gs, h, A, b, c, lu=None, Kf=None):
        self.Cn, self.Gs, self.h = Cn, Gs, h
        self.A, self.b, self.c = A, b, c
        self.lu, self.Kf = lu, Kf
        self.m = Cn.shape[0]
        self.s = A.shape[0]

    def solve(self, carry, forcing=None):
        """The last stage of ``J Z = [C_n carry + sum(forcing[i])]_i`` --
        the step map on `carry` (a vector or an `m x k` block, real or
        complex) plus a per-stage forcing, a tuple of terms per stage added
        in order: ``(dh/dT) S_i`` for a period column, ``w S_i - h U_i`` for
        an event column, a source's `sources`."""
        s, m = self.s, self.m
        base = self.Cn @ carry
        if self.lu is not None:
            blocks = []
            for i in range(s):
                r = base
                if forcing is not None:
                    for term in forcing[i]:
                        r = r + term
                blocks.append(r)
            Z = _lu_solve_split(self.lu, np.vstack(blocks) if np.ndim(carry) == 2
                                else np.concatenate(blocks))
            return Z[(s - 1) * m:s * m]
        A, h, Gs = self.A, self.h, self.Gs
        D = [None] * s
        for i in range(s):
            if self.Kf[i] is None:
                D[i] = carry          # the explicit first stage: D_0 = carry
                continue
            r = base
            if forcing is not None:
                for term in forcing[i]:
                    r = r + term
            r = r - h * sum(A[i, jj] * (Gs[jj] @ D[jj]) for jj in range(i))
            D[i] = np.asarray(_complex_solve(self.Kf[i], r))
        return D[s - 1]

    def adjoint(self, w):
        """The transposed step: ``(M_j^T w, r)`` with `r` the per-stage
        costates, ``r = J^{-T} [0; ..; 0; w]`` and ``M_j^T w = C_n^T sum_i
        r_i``.  Coupled: one transposed solve of the block.  Lower
        triangular: reverse mode on the forward recursion, stages in
        REVERSE -- ``r_i = K_i^{-T} dbar_i``, ``wbar += C_n^T r_i``, ``dbar_j
        -= h A_ij G_j^T r_i`` for ``j < i`` -- and the explicit first stage
        (no costate, ``r_0 = None``) passes ``dbar_0`` straight to `wbar`."""
        s, m = self.s, self.m
        if self.lu is not None:
            p = _lu_solve_split(
                self.lu, np.concatenate([np.zeros(m)] * (s - 1) + [w]), trans=1)
            r = [p[i * m:(i + 1) * m] for i in range(s)]
            return self.Cn.T @ sum(r), r
        A, h, Gs = self.A, self.h, self.Gs
        dbar = [np.zeros_like(w) for _ in range(s)]
        dbar[s - 1] = w
        wbar = np.zeros_like(w)
        r = [None] * s
        for i in range(s - 1, -1, -1):
            if self.Kf[i] is None:
                wbar = wbar + dbar[i]
                continue
            rb = _complex_solve_transposed(self.Kf[i], dbar[i])
            if rb is None:
                raise NotImplementedError(
                    'PSS: this linear solver cannot solve transposed, so a '
                    'stage method\'s adjoint cannot be replayed. Use '
                    'DenseSolver or SuperLUSolver.')
            r[i] = rb
            wbar = wbar + self.Cn.T @ rb
            for jj in range(i):
                dbar[jj] = dbar[jj] - h * A[i, jj] * (Gs[jj].T @ rb)
        return wbar, r

    def reach(self, r, k):
        """``sum_i A_ik r_i`` -- what the stage costates `r` read of a source
        at abscissa `k` (every stage residual carries it as ``-h A_ik u_k``);
        None when no stage reaches it.  A lower-triangular tableau's sum
        starts at the diagonal (``A_ik = 0`` above it)."""
        lo = 0 if self.lu is not None else k
        cp = sum(self.A[i, k] * r[i] for i in range(lo, self.s)
                 if r[i] is not None)
        return None if np.isscalar(cp) else cp

    def sources(self, u, jw, ts, _te):
        """A source ``u e^{jw t}`` on the step from `ts`, as the stage
        residuals carry it: ``-h sum_k A_ik u e^{jw (ts + c_k h)}`` per stage
        `i`, one-term forcing tuples for `solve`."""
        s, A, c, h = self.s, self.A, self.c, self.h
        return [(-h * sum(A[i, k] * u * np.exp(jw * (ts + c[k] * h))
                          for k in range(s if self.lu is not None else i + 1)),)
                for i in range(s)]

    def source_adjoint(self, acc, r, jw, ts, _te):
        """`acc` less the step's source coupling to the costates `r` --
        ``acc - h sum_k e^{jw (ts + c_k h)} sum_i A_ik r_i``, the transpose of
        `sources`."""
        for k in range(self.s):
            cp = self.reach(r, k)
            if cp is not None:
                acc = acc - self.h * np.exp(jw * (ts + self.c[k] * self.h)) * cp
        return acc

    def source_response(self, i_src):
        """``d x_{n+1} / d u`` (`m x m`) for a unit source entering stage
        `i_src` alone -- the stage residuals carry it as ``-h A_{i,i_src} u``."""
        m = self.m
        return self.solve(np.zeros((m, m)),
                          [(-self.h * self.A[i, i_src] * np.eye(m),)
                           for i in range(self.s)])


def _butcher(integ):
    """A Runge-Kutta integrator's `(A, b, c)` as float arrays."""
    return tuple(np.array(v, dtype=float) for v in integ.butcher())


def _lmm_recursion(Px, Cs, Pq, C_new, alphas, b, solve, source=None,
                   forcing=()):
    """One step of the linear-multistep sensitivity recursion (see
    `PSS._step_sensitivity`, which documents it):

        S    = sum_{k>=1} a_k C_{n-k} P_{n-k} + b Pq + sum(forcing)
        P_n  = -Jf_n^-1 (S + source)
        Pq_n = a_0 C_n P_n + S

    `solve` takes the `Jf` solve.  ⚠ `forcing` and `source` are different
    things: a `forcing` term is part of the step's COMPANION derivative (the
    period column's `dr/dh`, an event column's step-size and source-time
    terms), so it enters `Pq` too; a small-signal `source` is an injected
    current, not a charge, so it enters the solve alone.  The forcing terms
    are added one at a time, in order.  Returns `(P_n, Pq_n)`."""
    S = b * Pq if b else np.zeros_like(Px[0])
    for k in range(1, len(alphas)):
        S = S + alphas[k] * (Cs[k - 1] @ Px[k - 1])
    for f in forcing:
        S = S + f
    S_solve = S if source is None else S + source
    Px_new = -solve(S_solve)
    Pq_new = alphas[0] * (C_new @ Px_new) + S
    return Px_new, Pq_new


class _LMMStep(object):
    """One step of a linear multistep method's period map -- the plain map
    (euler, trap, theta) and gear's solved-history pair alike -- with the
    `_StageStep` interface, built from the stored record `(lu, C_new, alphas,
    b)`, the capacitance ring the forward pass saw (`C1 = C_{n-1}`, `C2 =
    C_{n-2}`) and the step's end time.

    The per-step state is ``(P_n, P_{n-1}, Pq_n)``: the recursion
    `_lmm_recursion` reads ``a_k C_{n-k} P_{n-k}`` for ``k <= 2`` and ``b
    Pq``.  Which of it the MAP exposes -- ``P_n`` for the plain map, the pair
    ``(P_n, P_{n-1})`` for gear's -- is the `FactoredPeriod`'s business, and
    so is the opening; the step algebra is one.

    ⚠ ONE TRANSPOSE FOR EVERY COMPANION, the plain map's and gear's pair
    alike: with ``(w1, w2, w3)`` the adjoints of ``(P_n, P_{n-1},
    Pq_n)``,

        Sbar  = w3 - Jf^-T (w1 + a_0 C_n^T w3)
        (P_{n-1}, P_{n-2}, Pq_{n-1})bar = (a_1 C_{n-1}^T Sbar + w2,
                                           a_2 C_{n-2}^T Sbar,  b Sbar)

    -- trap's shared bracket (``b != 0``, one-step) and gear's ``(-a_1 C^T t
    + w2, -a_2 C^T t)`` (``b = 0``, ``Sbar = -t``) are its two special cases,
    operation for operation.

    History: `doc/shooting_history.md`, `_LMMStep`."""

    __slots__ = ('lu', 'C_new', 'alphas', 'b', 'C1', 'C2', 't_end')

    def __init__(self, record, C1, C2, t_end):
        self.lu, self.C_new, self.alphas, self.b = record
        self.C1, self.C2, self.t_end = C1, C2, t_end

    def solve(self, carry, forcing=None):
        """The step on the state `(P_n, P_{n-1}, Pq_n)`; `forcing` is a
        source vector (see `sources`)."""
        P1, P2, Pq = carry
        Px_new, Pq_new = _lmm_recursion(
            [P1, P2], [self.C1, self.C2], Pq, self.C_new, self.alphas,
            self.b, lambda S, _l=self.lu: _complex_solve(_l, S), forcing)
        return (Px_new, P1, Pq_new)

    def adjoint(self, w):
        """The transposed step on `(w1, w2, w3)` (see the class note), and
        the transposed solve ``t = Jf^-T (...)`` -- the step's costate, which
        a source entering the solve reads."""
        w1, w2, w3 = w
        a, b = self.alphas, self.b
        rhs = w1 + a[0] * (np.asarray(self.C_new).T @ w3) if b else w1
        t = _complex_solve_transposed(self.lu, rhs)
        if t is None:
            raise NotImplementedError(
                'PSS: this linear solver cannot solve transposed, so the '
                'monodromy transpose cannot be replayed. Use DenseSolver or '
                'SuperLUSolver.')
        Sbar = (w3 - t) if b else -t
        p1 = a[1] * (np.asarray(self.C1).T @ Sbar) + w2
        p2 = (a[2] * (np.asarray(self.C2).T @ Sbar) if len(a) > 2
              else np.zeros_like(w1))
        pq = b * Sbar if b else np.zeros_like(w1)
        return (p1, p2, pq), t

    def sources(self, u, jw, _ts, te):
        """A source ``u e^{jw t}``: it enters the step's solve (not the
        companion -- an injected current is not a charge) at the step's
        END, ``t_{n+1}``."""
        return u * np.exp(jw * float(te))

    def source_adjoint(self, acc, t, jw, _ts, te):
        """`acc` less the source's coupling to the transposed solve `t`: the
        forward source is ``P_n = -Jf^-1 (S + u e^{jw t_{n+1}})``."""
        return acc - np.exp(jw * float(te)) * np.asarray(t)


class _GLMStep(object):
    """One step of a Nordsieck GLM's period map with the `_StageStep`
    interface, on the MULTIVALUE state (`r` blocks of width `m`), from the
    record `(Kfacs, Gs, h, A, U, B, V, Ks)` of `_glm_period_blocks`."""

    __slots__ = ('rec',)

    def __init__(self, rec):
        self.rec = rec

    def solve(self, P, forcing=None):
        if forcing is not None:
            raise NotImplementedError('_GLMStep: no forcing')
        return _glm_step(self.rec, P)[0]

    def adjoint(self, W):
        """The reverse-mode adjoint of one step (see
        `PSS._monodromy_matvec_transposed`): per step, with `W` the adjoint
        of the output vector,

            Dbar_i  = -h sum_k B_ki G_i^T W_k
            Pbar_j  =  sum_k V_kj W_k
            for i = s-1 .. 0:   rbar_i = K_i^{-T} Dbar_i
                                Pbar_j += U_ij rbar_i          (all j)
                                Dbar_j += -h A_ij G_j^T rbar_i (j < i)

        Returns `(Pbar, rbars)`."""
        Kfacs, Gs, h, A, U, B, V, _Ks = self.rec
        s = len(Kfacs)
        r = len(W)
        Dbar = [-h * sum(B[k, i] * (Gs[i].T @ W[k]) for k in range(r))
                for i in range(s)]
        Pbar = [sum(V[k, jj] * W[k] for k in range(r)) for jj in range(r)]
        rbars = [None] * s
        for i in range(s - 1, -1, -1):
            rb = Kfacs[i].solve_transposed(Dbar[i])
            if rb is None:
                raise NotImplementedError(
                    'PSS: this linear solver cannot solve transposed, so '
                    'the GLM monodromy transpose cannot be replayed. Use '
                    'DenseSolver or SuperLUSolver.')
            rbars[i] = rb
            for jj in range(r):
                Pbar[jj] = Pbar[jj] + U[i, jj] * rb
            for jj in range(i):
                Dbar[jj] = Dbar[jj] - h * A[i, jj] * (Gs[jj].T @ rb)
        return Pbar, rbars

    def sources(self, *_args):
        raise NotImplementedError(
            'PAC: the forced (small-signal) response is not built on a '
            "Nordsieck GLM's own map (monodromy='native') -- its sources "
            'would enter every stage and the Nordsieck output rows.  A GLM '
            "run reads its small-signal response from a twin: leave "
            "monodromy at 'trbdf2' or 'radau'.")

    source_adjoint = sources


class _GLMStartup(object):
    """A Nordsieck GLM's STARTUP, linearised (`PSS._glm_startup_linearisation`):
    how the starting vector ``Q_k`` (k = 0..p) moves with the unknown `x_0`
    -- and, on an autonomous circuit, with the period.

    The startup takes p Radau IIA(3) substeps of ``h_s = h/p`` from `x_0`
    and interpolates their charges, ``Q_k = p^k k! sum_j (V^-1)_kj q(x_j)``
    (``Q_0 = q(x_0)``), so

        dQ_k = sum_j W_kj C(x_j) dx_j,      W_kj = p^k k! (V^-1)_kj

    with ``dx_j`` carried through each substep's converged stage system:
    ``dY = J^-1 [C(x_{j-1}) dx_{j-1}]_i``, ``J[i][l] = delta_il C(Y_i) + h_s
    A_il G(Y_l)``, ``dx_j = dY_3`` (stiff accuracy).  The period enters
    through ``h_s``: ``dY/dT`` gains ``J^-1 (h_s/T) [sum_l A_il K_l]_i``.

    ⚠ THIS IS WHAT MAKES THE GLM'S MAP ON `x` EXACT.  Seeding the Nordsieck
    recursion with ``[C(x_0), 0, ...]`` -- the startup's higher components
    held fixed -- made the shooting Jacobian approximate and the factored
    map not the Newton's Jacobian (matrix-free, glm2 and glm3 diverged on a
    driven RLC)."""
    __slots__ = ('Cx', 'lus', 'fT', 'W', 'hs', 'm')

    def __init__(self, Cx, lus, fT, W, hs):
        self.Cx, self.lus, self.fT, self.W, self.hs = Cx, lus, fT, W, hs
        self.m = Cx[0].shape[0]

    def _substeps(self, d0, T=None):
        """``dx_j``, j = 0..p, from ``dx_0 = d0`` (a vector or a block);
        with `T`, the period's own forcing added on every substep."""
        from scipy.linalg import lu_solve
        m = self.m
        ds = [d0]
        for j, lu in enumerate(self.lus):
            blk = self.Cx[j] @ ds[-1]
            rhs = (np.concatenate([blk] * 3) if blk.ndim == 1
                   else np.vstack([blk] * 3))
            if T is not None:
                rhs = rhs + (self.hs / float(T)) * self.fT[j]
            ds.append(lu_solve(lu, rhs)[2 * m:3 * m])
        return ds

    def _combine(self, ds):
        return [sum(self.W[k, j] * (self.Cx[j] @ ds[j])
                    for j in range(len(ds)))
                for k in range(self.W.shape[0])]

    def matrix(self):
        """`dQ_k/dx_0`, `r` blocks of ``m x m``."""
        return self._combine(self._substeps(np.eye(self.m)))

    def matvec(self, v):
        """`dQ_k/dx_0 v`, `r` vectors."""
        return self._combine(self._substeps(np.asarray(v, dtype=float).ravel()))

    def dT(self, T):
        """`dQ_k/dT` with `x_0` held, `r` vectors (the substeps scale with
        the period)."""
        return self._combine(self._substeps(np.zeros(self.m), T=T))


def _glm_step(rec, P, fT=None):
    """One step of the GLM sensitivity recursion (`_glm_propagate`'s body):
    the stage solves on the multivalue input `P` (`r` blocks), and the
    output vector.  `fT` (not None) adds the step's period-column forcing
    ``fT sum_j A_ij K_j`` / ``fT sum_i B_ki K_i``.  Returns `(P_out, D)`."""
    Kfacs, Gs, h, A, U, B, V, Ks = rec
    s = len(Kfacs)
    r = len(P)
    D = [None] * s
    for i in range(s):
        rhs = sum(U[i, j] * P[j] for j in range(r))
        if fT is not None:
            rhs = rhs + fT * sum(A[i, j] * Ks[j] for j in range(i + 1))
        for j in range(i):
            rhs = rhs - h * A[i, j] * (Gs[j] @ D[j])
        D[i] = Kfacs[i].solve(rhs)
    P = [sum(V[k, j] * P[j] for j in range(r))
         - h * sum(B[k, i] * (Gs[i] @ D[i]) for i in range(s))
         + (fT * sum(B[k, i] * Ks[i] for i in range(s))
            if fT is not None else 0.0)
         for k in range(r)]
    return P, D
