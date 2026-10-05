"""The coupled Radau IIA(3) stages: the full and the transformed Newton.  A
theme of `Transient` (see `transient.py`).
"""

import os

import numpy as np

from pycircuit.circuit import (
    _evalhint,
    _paths,
    _tran_core,
    _tran_radau_c,
    _tran_radau_tc,
)
from pycircuit.circuit import pcnr as _pcnr
from pycircuit.circuit._limiting import (
    limit_sync,
    limiter_snapshot,
    state_restore,
    stateful_limiters,
)
from pycircuit.circuit.analysis import (
    NoConvergenceError,
    remove_row_col,
)
from pycircuit.circuit.simwarnings import (
    UsageWarning,
    warn,
)

#: the transform's frozen factors made once a step (`_radau_frozen`, speed
#: round 10, B3.1); env `PYCIRCUIT_RADAU_FROZEN=0` makes them every iteration
RADAU_FROZEN = os.environ.get('PYCIRCUIT_RADAU_FROZEN', '1') != '0'
#: the stages' readers' passes in one core call a stage (`_stage_end_passes`,
#: speed round 10, B3.2); env `PYCIRCUIT_STAGE_FUSE=0` leaves each reader its
#: own evaluation
STAGE_FUSE = os.environ.get('PYCIRCUIT_STAGE_FUSE', '1') != '0'
#: numpy's solve as imported: a caller's stand-in takes the per-iteration path
_NP_SOLVE = np.linalg.solve
_FROZEN = {}


def _transform_rhs(F0, F1, F2, p):
    """The transform's right-hand sides ``-(P[k,0] F0 + P[k,1] F1 + P[k,2]
    F2)``, k = 0, 1 (`p` the six entries).  ONE PLACE for the per-iteration
    and the frozen solve: a warning their arithmetic gives comes from one
    line, as before the frozen solve (numpy's warnings are shown once a
    line)."""
    p00, p01, p02, p10, p11, p12 = p
    rhs0 = -(p00 * F0 + p01 * F1 + p02 * F2)
    rhs1 = -(p10 * F0 + p11 * F1 + p12 * F2)
    return rhs0, rhs1


def _transform_back(w0, w1, v):
    """``dY_i = V[i,0].real w0 + 2 Re(V[i,1] w1)``, concatenated (`v` the
    pairs ``(V[i,0].real, V[i,1])``): one place for both solves, as
    `_transform_rhs`."""
    dY = [a * w0 + 2.0 * np.real(b * w1) for a, b in v]
    return np.concatenate(dY)


class _FrozenTransform:
    """One step's transform, frozen (`_RadauStages._radau_frozen`): its `P`
    entries, `V`'s, the real factor with its kept LU (or the analysis
    solver to call) and the complex factor marshalled for its solver (or
    the per-call complex solve)."""

    __slots__ = ('cf', 'cs', 'ls', 'lu', 'm', 'p', 'prep', 'rf', 'tk', 'v', 'zs')

    def solve(self, R3):
        """`_RadauStages._radau_transform_solve(R3, Cr, Gr, h)` against the
        frozen factors: its operations, in its order.  The real solve is the
        kept LU's for a finite right-hand side (numpy's solve, bit for bit),
        the analysis solver's otherwise."""
        m = self.m
        R3 = np.asarray(R3, dtype=float)
        rhs0, rhs1 = _transform_rhs(R3[0:m], R3[m:2 * m], R3[2 * m:3 * m], self.p)
        b0 = np.real(rhs0)
        if self.lu is not None and np.isfinite(b0).all():
            w0 = self.lu.solve(b0)
        else:
            w0 = np.asarray(self.ls.solve(self.rf, b0, self.tk), dtype=float)
        if self.prep is not None:
            w1 = np.asarray(self.zs.solve_prepared(self.prep, rhs1), dtype=complex)
        else:
            w1 = np.asarray(self.cs(self.cf, rhs1), dtype=complex)
        return _transform_back(w0, w1, self.v)


def _frozen_chain():
    """The transform's pieces as defined (`_paths.genuine`), read once; None
    while one is not."""
    if 'chain' not in _FROZEN:
        from pycircuit.circuit import _numeric
        from pycircuit.circuit.analysis import Analysis
        from pycircuit.circuit.linearsolver import (
            AutoSolver,
            ComplexKLUSolver,
            DenseSolver,
        )
        from pycircuit.circuit.toolkit import NumericToolkit
        R, L = 'pycircuit.circuit._tran_radau', 'pycircuit.circuit.linearsolver'
        own = ((_RadauStages._radau_transform_solve, '_RadauStages._radau_transform_solve', R),
               (_RadauStages._radau_complex_solve, '_RadauStages._radau_complex_solve', R),
               (_RadauStages._radau_transform_matrices,
                '_RadauStages._radau_transform_matrices', R),
               (Analysis._get_linearsolver, 'Analysis._get_linearsolver',
                'pycircuit.circuit.analysis'),
               (AutoSolver.solve, 'AutoSolver.solve', L),
               (AutoSolver._select, 'AutoSolver._select', L),
               (DenseSolver.solve, 'DenseSolver.solve', L),
               (_numeric.linearsolver, 'linearsolver', 'pycircuit.circuit._numeric'),
               (ComplexKLUSolver.solve, 'ComplexKLUSolver.solve', L),
               (ComplexKLUSolver.prepare, 'ComplexKLUSolver.prepare', L),
               (ComplexKLUSolver.solve_prepared, 'ComplexKLUSolver.solve_prepared', L))
        if not all(_paths.genuine(*o) for o in own):
            return None
        _FROZEN['chain'] = tuple(o for o, _q, _m in own)
        _FROZEN['types'] = (ComplexKLUSolver, AutoSolver, DenseSolver, NumericToolkit, _numeric)
    return _FROZEN['chain']


class _RadauStages:
    """The coupled Radau IIA(3) stages: the full and the transformed Newton.  A
    theme of `Transient` (see `transient.py`)."""

    def _coupled_stage_context(self, x0, t, provided_function):
        """What every Radau IIA(3) coupled step shares before its Newton:
        the tableau, the step, the stage times, the source, the reduced-vector
        map (the reference row dropped) and the entering charge."""
        from types import SimpleNamespace
        integ = self.base_integrator
        tk = self.toolkit
        iref = self.irefnode
        h = self._dt
        tn = t - h
        cvec = np.array(integ.C, dtype=float)
        arr = lambda v: tk.array(v, dtype=float)

        def red(v):
            return np.concatenate((np.asarray(v)[:iref],
                                   np.asarray(v)[iref + 1:]))
        _rec = self._memo_get(x0)
        if _rec is None or 'q' not in _rec:
            ## (the evaluate core's passes where it serves: speed round 9,
            ## stage 8 -- `_tran_core.passes`)
            P = _tran_core.passes(self, x0, 'q')
            qn = P['q'] if P is not None else arr(self.cir.q(x0, self.epar))
        else:
            qn = _rec['q']
        return SimpleNamespace(
            Amat=np.array(integ.A, dtype=float), h=h, tn=tn, iref=iref,
            arr=arr, red=red, src=self._stage_source(provided_function),
            qn=qn, m=red(qn).shape[0],
            tstage=[tn + cvec[i] * h for i in range(3)])

    def _coupled_stage_system(self, ctx, qi, Ki, Ci=None, Gi=None):
        """The coupled stage residual ``R_i = q(Y_i) - q(x_n) - h sum_j A_ij
        K_j`` (reduced) and, with `Ci`/`Gi`, its block Jacobian ``J[i][j] =
        delta_ij C_i + h A_ij G_j`` -- from per-stage lists, whatever
        linearisation produced them (limited device evaluations, PCNR's
        Schur-reduced conductance, ...)."""
        m, h, A, red = ctx.m, ctx.h, ctx.Amat, ctx.red
        R = np.empty(3 * m)
        Jbig = None if Ci is None else np.zeros((3 * m, 3 * m))
        for i in range(3):
            Fi = qi[i] - ctx.qn - h * sum(A[i, j] * Ki[j] for j in range(3))
            R[i * m:(i + 1) * m] = red(Fi)
            if Jbig is not None:
                for j in range(3):
                    if i == j:
                        blk = Ci[i] + h * A[i, j] * Gi[j]
                    else:
                        blk = h * A[i, j] * Gi[j]
                    (blk_r,) = remove_row_col((blk,), ctx.iref, self.toolkit)
                    Jbig[i * m:(i + 1) * m, j * m:(j + 1) * m] = \
                        np.asarray(blk_r)
        return R, Jbig

    def _rk_step_coupled_pcnr(self, x0, t, provided_function=None):
        """The coupled Radau IIA(3) step with PCNR as the limiting, IN EVERY
        STAGE, instead of per-device ``cir.limit``.

        Structurally identical to :meth:`_rk_step_coupled` -- the SAME coupled
        ``3m`` Newton on ``F_i = q(Y_i) - q(x_n) - h sum_j A_ij K_j`` -- but the
        stage residual current and Jacobian are the JUNCTION-CONTINUATION
        effective ones.  For each stage ``j`` an augmented system is built at the
        device limiting voltages ``v_lim[j]`` and Schur-reduced onto the MNA
        size, giving ``i_eff_j``/``G_eff_j`` (:func:`pcnr.schur_reduce`, with
        ``u_extra=0`` so it is the plain effective ``i``/``G``, the companion
        living in ``q`` and ``A`` as in the non-PCNR coupled step).  These
        replace ``cir.i``/``cir.G`` in the identical assembly, so
        ``K_j = -(i_eff_j + u(t_j))`` and ``J[i][j] = delta_ij C + h A_ij
        G_eff_j``.  The coupled solve delivers each stage's ``dx_MNA``; the
        limiting voltages advance by the paper's correct phase (``dx_lim`` from
        :func:`pcnr.dx_lim_of` on the coupled ``dx_MNA``, then
        :func:`pcnr.refine`), never by ``cir.limit``.

        At convergence ``v_lim`` equals the junction voltage, so
        ``i_eff``/``G_eff`` equal ``cir.i``/``cir.G`` and the accepted point is
        bit-identical to what device limiting reaches on a single-junction stage
        -- while PARALLEL junctions on one branch, which per-device limiting
        cannot resolve (they fight over the shared branch voltage, order-
        dependently), are limited JOINTLY here.  Same return and side-effect
        contract as :meth:`_rk_step_coupled`.
        """
        junctions = _pcnr.pcnr_devices(self.cir)
        ctx = self._coupled_stage_context(x0, t, provided_function)
        iref = ctx.iref
        arr, red, src, tstage, m = (ctx.arr, ctx.red, ctx.src, ctx.tstage,
                                    ctx.m)
        epar = self.epar
        tk = self.toolkit
        xn = x0

        ## Stage values and their PCNR limiting voltages.  `v_lim[j]` is the
        ## per-stage state that stands in for each device's internal `_vlim`;
        ## seeded (and limited) from the stage guess, exactly like the DC and
        ## single-stage paths.
        ## STAGE PREDICTOR -- and here it seeds the LIMITING too: `v_lim_init`
        ## reads the stage guess, and limiting the seed is what fixed PCNR's
        ## one documented failure (`pcnr.py`), so a seed nearer the answer is
        ## the same medicine.
        Y = [self._pred_or(np.array(xn, dtype=float), tstage[j])
             for j in range(3)]
        v_lim = [_pcnr.v_lim_init(junctions, Y[j]) for j in range(3)]
        reltol = self.par.reltol
        abstol = float(self.par.vabstol)
        maxit = int(self.par.maxiter)
        converged = False
        for _ in range(maxit):
            qi, Ki, Ci, Geff, glim = [], [], [], [], []
            for j in range(3):
                qi.append(arr(self.cir.q(Y[j], epar)))
                Ci.append(arr(self.cir.C(Y[j], epar)))
                ## EFFECTIVE i/G at the stage's limiting voltages: the augmented
                ## junction system Schur-reduced onto the MNA size.  With
                ## `u_extra=0`/`J_extra=0` these are the plain effective `i`/`G`
                ## (the companion stays in q and A, not folded in as a DC-flow
                ## source the way the sequential DIRK stage PCNR does it).
                g_mna, g_lim, J_mm, _J_ml, _J_lm, didv = _pcnr.augmented_system(
                    self.cir, Y[j], v_lim[j], junctions, epar,
                    u_extra=0.0, dense_blocks=False, J_extra=0.0)
                ## ⚠ THE RESIDUAL CURRENT IS `g_mna`, NOT the Schur RHS
                ## `f_eff`.  `augmented_system` STAMPS the junction current at
                ## `v_lim` into `g_mna` (via `dev.stamp`), so `g_mna` is the
                ## PHYSICAL MNA current with junctions at `v_lim` -- exactly
                ## `cir.i` once `v_lim` == the branch voltage.  `f_eff =
                ## g_mna - J_ml g_lim` folds the junction current into the
                ## Newton-STEP right-hand side instead, where it VANISHES as
                ## `g_lim -> 0`; using it as the residual drops the junction
                ## current at convergence and converges the coupled Newton to a
                ## neighbouring, wrong root.  Only `G_eff = J_eff` (the
                ## Schur-reduced Jacobian) is taken here; the junction is
                ## eliminated from the step by the correct phase (`dx_lim_of` +
                ## `refine`) below.
                ## History: `doc/transient_history.md`, `Transient._rk_step_coupled_pcnr`.
                _f_eff, G_eff = _pcnr.schur_reduce(
                    g_mna, g_lim, J_mm, junctions=junctions, didv=didv)
                Ki.append(-(np.asarray(g_mna, dtype=float) + src(tstage[j])))
                Geff.append(np.asarray(G_eff, dtype=float))
                glim.append(g_lim)
            ## residual blocks (reduced) and dense 3m x 3m Jacobian -- identical
            ## to _rk_step_coupled with i_eff/G_eff in place of cir.i/cir.G.
            R, Jbig = self._coupled_stage_system(ctx, qi, Ki, Ci, Geff)
            dY = np.linalg.solve(Jbig, -R)
            scale = 0.0
            lim_ok = True
            for i in range(3):
                di = dY[i * m:(i + 1) * m]
                di_full = np.asarray(tk.insert(di, iref, 0.0))
                ## The correct phase: dx_lim from THIS stage's g_lim and the
                ## coupled dx_MNA, then each device limits only its own probes
                ## (refine).  No `cir.limit` anywhere -- PCNR IS the limiting.
                dx_lim = _pcnr.dx_lim_of(junctions, glim[i], di_full)
                Y_prev = Y[i]
                v_new = _pcnr.refine(junctions, v_lim[i], v_lim[i] + dx_lim,
                                     epar, x_old=Y_prev)
                if not _pcnr.lim_converged(glim[i], v_new, reltol,
                                           self.par.vabstol):
                    lim_ok = False
                Y[i] = Y_prev + di_full
                v_lim[i] = v_new
                scale = max(scale, np.max(np.abs(di)))
            ynorm = max(np.max(np.abs(red(Y[i]))) for i in range(3))
            if lim_ok and scale <= reltol * ynorm + abstol:
                converged = True
                break
        if not converged:
            ## ⚠ THERE IS DELIBERATELY NO CONTINUATION LADDER HERE.  The step
            ## instead FALLS BACK to the device-limiting coupled solve (see the
            ## caller, `_rk_step_coupled`), which carries one.
            ##
            ## No circuit is known where PCNR fails at a normal iteration
            ## budget, and a transient suppresses its failure mode structurally
            ## -- PCNR fails on a far-off INITIAL GUESS, and every step here
            ## starts from the last accepted state.  The design a ladder would
            ## take, and the failure signature that would justify building it,
            ## are in the history.  ⚠ Do NOT validate a candidate by starving
            ## `maxiter`: a ladder must end with a PURE solve of the original
            ## system, so a starved budget defeats the final rung whatever the
            ## deformation.
            ## History: `doc/transient_history.md`, `Transient._rk_step_coupled_pcnr`.
            raise NoConvergenceError(
                'Radau IIA(3) coupled PCNR stage Newton did not converge')

        ## SYNC each device's internal `_vlim` to the converged solution: PCNR
        ## never called `cir.limit`, so the epilogue's `i`/`G` at the step
        ## end (`_iq`, the returned J) would otherwise linearise at a stale
        ## `_vlim`.  The step end is the last stage, and `limit_sync` puts
        ## the state exactly there (see _rk_stage_pcnr); the chain of
        ## `limit(Y[j], Y[j])` it replaces reached it only while consecutive
        ## stages were close (its first syncs landed up to 25 mV short).
        limit_sync(self.cir, Y[-1], epar)
        return self._finish_radau(ctx, x0, t, provided_function, Y)

    def _stage_limiter_states(self, Y, lims, epar):
        """Each coupled stage's own limiting state, `S[j]` (None per stage
        without stateful limiters): the step's entering state synced to the
        stage's seed `Y[j]`.  The full and the simplified (transformed)
        coupled Newton's one set-up (written out in each until 2026-10-01,
        the review's O13); why each stage owns one is at
        `_coupled_stage_solver`."""
        if not lims:
            return [None] * 3
        S0 = limiter_snapshot(lims)
        S = []
        for j in range(3):
            state_restore(S0)
            self.cir.limit(Y[j], Y[j], epar)
            S.append(limiter_snapshot(lims))
        return S

    def _stages_converged(self, Y, S, lims, scale, reltol, abstol, red):
        """The coupled Newton's convergence test: the largest stage update
        `scale` within `reltol` of the largest reduced stage (`red`) plus
        `abstol`, AND every stateful limiter at rest (`_limiters_at_rest`):
        a small step is not a converged one while a limiter still holds a
        stage short of its node voltage."""
        ynorm = max(np.max(np.abs(red(Y[i]))) for i in range(3))
        return scale <= reltol * ynorm + abstol and (
            not lims or self._limiters_at_rest(Y, S, lims))

    def _limiters_at_rest(self, Y, S, lims):
        """Whether every stage's device current, read at the limiting state
        its Newton solved with (`S[j]`), is the device's EXACT current at
        `Y[j]` (`limit_sync`) -- within `reltol` and `iabstol`.

        ⚠ A SMALL STEP IS NOT CONVERGENCE WHILE A LIMITER IS CLAMPING.  A
        `Diode` reads its current as the tangent at `_vlim`; clamped far
        below the node voltage, that tangent carries almost nothing, and
        the Newton's next step is tiny because the diode is effectively
        OFF.  Measured with per-stage states and the step test alone: a
        diode driven through 1 ohm came back at 6.18 V on a 50 us step (the
        whole source voltage), "converged" after two iterations.  The
        ordinary Newton is saved by its residual test (`conv_f`); the
        coupled one tests the step alone, so it asks the limiters.  Each
        stage's state is put back."""
        epar = self.epar
        reltol = float(self.par.reltol)
        iabstol = float(self.par.iabstol)
        at_rest = True
        for j in range(3):
            state_restore(S[j])
            i_state = np.asarray(self.cir.i(Y[j], epar), dtype=float)
            limit_sync(self.cir, Y[j], epar, lims)
            i_exact = np.asarray(self.cir.i(Y[j], epar), dtype=float)
            state_restore(S[j])
            if np.any(np.abs(i_state - i_exact) > reltol * np.abs(i_exact) + iabstol):
                at_rest = False
        return at_rest

    def _coupled_stage_solver(self, x0, t, provided_function):
        """The dense coupled step's Newton and residual for the step
        entering at `x0`: ``(ctx, stage_newton, block_residual, seed0)``.
        `_rk_step_coupled` solves with them; the PCNR and transform paths
        build them only once the branch screen fires, to confirm a second
        root of the SAME step equation (its roots do not depend on which
        Newton found the first)."""
        ctx = self._coupled_stage_context(x0, t, provided_function)
        Amat, h, iref = ctx.Amat, ctx.h, ctx.iref
        arr, red, src, qn, tstage, m = (ctx.arr, ctx.red, ctx.src, ctx.qn,
                                        ctx.tstage, ctx.m)
        epar = self.epar
        tk = self.toolkit
        xn = x0

        ## Stage values, initialised at the previous solution.  Solved in the
        ## reference-removed space of dimension m = n-1 per stage; the coupled
        ## system is 3m.  Newton to the transient tolerances.
        reltol = self.par.reltol
        abstol = float(self.par.vabstol)
        maxit = int(self.par.maxiter)
        lims = stateful_limiters(self.cir)
        ## (a bypassed device returns its last evaluation for a NEARBY point,
        ## so its values depend on history: no memo then)
        _nobypass = float(getattr(epar, 'bypasstol', -1.0) or -1.0) < 0.0

        def _stage_newton(seed, gshunt=0.0, damped=False):
            """The coupled `3m` Newton, optionally with a node-to-ground shunt.

            ⚠ THE SHUNT ENTERS AS A CONDUCTANCE IN THE DEVICE CURRENT, not as
            `F + g x` on the residual.  This residual is in CHARGE units
            (`q - q_n - h sum A K`) and its Jacobian is `C + h A G`, so adding a
            conductance straight to either is dimensionally wrong -- `g` has to
            arrive where `G` does.  Adding it to `i` and to `G` is also what
            makes the deformed problem a real circuit (every node shunted to
            ground by `g`), so the ladder tracks a physical branch.
            """
            if not gshunt and not damped:
                ## (the undamped Newton in one C call where it serves, the
                ## memo left as below: `_tran_radau_c`, speed round 9, B2;
                ## None -- declined or handed back, nothing left behind)
                Yc = _tran_radau_c.solve(self, ctx, seed, src, provided_function, lims,
                                         _nobypass, reltol, abstol, maxit)
                if Yc is not None:
                    return Yc
            Y = [np.array(y, dtype=float) for y in seed]
            ## ⚠ EACH STAGE OWNS ITS LIMITING STATE.  The three stages are
            ## solved SIMULTANEOUSLY, and a stateful limiter (`Diode`) keeps
            ## ONE `_vlim`, which its `i` / `G` read.  Shared, each stage was
            ## limited and re-synced against ANOTHER stage's state, and the
            ## three states could settle into a cycle with none of them on its
            ## stage: a diode driven to 0.85 V came out 1.05e-4 V off the
            ## exact solution at 2.5 us steps, and 0.87 V off (1.67 V where
            ## 0.80 V is right) on the first 25 us step.  So each stage keeps
            ## its own (`S[j]`), swapped in before it is read or limited --
            ## the sequential DIRK's Newton per stage, which met the exact
            ## solution to 5.6e-16 on the same circuit.  Each starts from the
            ## step's entering state, synced to its seed.
            S = self._stage_limiter_states(Y, lims, epar)

            def assemble(Y, S):
                """Residual and block Jacobian at the stage vector `Y`, each
                stage read at its own limiting state `S[j]`."""
                qi, Ki, Ci, Gi = [], [], [], []
                for j in range(3):
                    if S[j] is not None:
                        state_restore(S[j])
                    ## (the stage's four passes in one evaluate-core call where
                    ## it serves -- `_tran_core.passes`, speed round 9, stage
                    ## 8; else one evaluation session: `_evalhint`)
                    P = _tran_core.passes(self, Y[j], 'qiGC')
                    if P is not None:
                        qi.append(P['q'])
                        i_j, G_j = P['i'], P['G']
                        Ci.append(P['C'])
                    else:
                        with _evalhint.evaluating('q', 'i', 'G', 'C'):
                            qi.append(arr(self.cir.q(Y[j], epar)))
                            i_j = arr(self.cir.i(Y[j], epar))
                            G_j = arr(self.cir.G(Y[j], epar))
                            Ci.append(arr(self.cir.C(Y[j], epar)))
                    ## (the evaluations at this stage, for the step's
                    ## readers -- `_memo_get`; never under a shunt or a
                    ## stateful limiter)
                    if not gshunt and not lims and _nobypass:
                        self._memo_put(Y[j], {'q': qi[-1], 'i': i_j,
                                              'C': Ci[-1], 'G': G_j})
                    if gshunt:
                        i_j = i_j + gshunt * np.asarray(Y[j], dtype=float)
                        G_j = G_j + gshunt * np.eye(np.asarray(G_j).shape[0])
                    Ki.append(-(i_j + src(tstage[j])))
                    Gi.append(G_j)
                ## residual blocks (reduced) and dense 3m x 3m Jacobian
                return self._coupled_stage_system(ctx, qi, Ki, Ci, Gi)

            ## ⚠ BACKTRACKING LINE SEARCH, AS A RETRY (owner decision).  A
            ## whole class of circuits cannot be integrated undamped at ANY
            ## grid: a nonlinearity fed by a branch current through a
            ## capacitor sees the companion difference quotient, whose stage
            ## sensitivity is tau/h, and the undamped basin is 0.94 h/(k tau)
            ## (READING-LOG 2.156) -- refining the grid SHRINKS it in
            ## proportion.  The rule and floor are `DampedNewton`'s (Armijo
            ## 1e-4, alpha >= 0.05).
            ## ⚠⚠ WHY A RETRY AND NOT A BLANKET SEARCH: a blanket line search
            ## CHANGES CONVERGED ANSWERS.  At h = 100 tau on a tanh charge
            ## circuit the coupled stage system has more than one root (its
            ## a_23 < 0 makes the coupled equations non-monotone where the
            ## scalar map is monotone); the undamped-plus-limited path finds
            ## the physical one, v in [0, 1], and the damped path a spurious
            ## one at v = -0.010.  So `damped=False` (the first attempt) is the
            ## undamped Newton -- the full step always, the undamped
            ## convergence test on it -- and the search is tried only after it
            ## has failed, before the shunt ladder.
            ## History: `doc/transient_history.md`, `Transient._coupled_stage_solver`.
            R, Jbig = assemble(Y, S)
            converged = False
            for _ in range(maxit):
                dY = np.linalg.solve(Jbig, -R)
                Rnorm = float(np.sum(np.abs(R)))
                alpha = 1.0
                while True:
                    Y_trial, S_t = [], []
                    scale = 0.0
                    for i in range(3):
                        if S[i] is not None:
                            state_restore(S[i])
                        di = alpha * dY[i * m:(i + 1) * m]
                        ## LIMITING IS LOAD-BEARING ON A NONLINEAR JUNCTION.
                        ## Without it the coupled Newton on a diode overshoots
                        ## the exponential and settles on a spurious
                        ## near-linear solution (a mixer produces a pure
                        ## sinusoid with no harmonics -- measured).
                        Y_prev = Y[i]
                        Y_new = self.cir.limit(Y_prev + tk.insert(di, iref, 0.0),
                                               Y_prev, epar)
                        Y_trial.append(Y_new)
                        S_t.append(limiter_snapshot(lims) if lims else None)
                        step_i = red(np.asarray(Y_new) - np.asarray(Y_prev))
                        scale = max(scale, np.max(np.abs(step_i)))
                    R_t, J_t = assemble(Y_trial, S_t)
                    ## ⚠ THE UNDAMPED CONVERGENCE TEST, ON THE FULL STEP,
                    ## BEFORE THE LINE SEARCH SEES IT.  Testing the size of a
                    ## DAMPED step lets a step at the floor stop the iteration
                    ## short of the root; requiring a full step to pass Armijo
                    ## hangs at the roundoff floor, where both residual norms
                    ## are noise and the test fails by chance.  A converged
                    ## full step is accepted here exactly as the undamped
                    ## Newton accepts it; the search engages only when the full
                    ## step is neither small nor residual-reducing.
                    ## ⚠ AND ONLY WITH EVERY LIMITER AT REST (`_limiters_at_rest`):
                    ## a small step is not a converged one while a stateful
                    ## limiter still holds a stage short of its node voltage.
                    ## History: `doc/transient_history.md`, `Transient._coupled_stage_solver`.
                    if alpha == 1.0 and self._stages_converged(
                            Y_trial, S_t, lims, scale, reltol, abstol, red):
                        converged = True
                        break
                    if not damped:
                        break                      # the undamped Newton: the full step, always
                    if float(np.sum(np.abs(R_t))) <= Rnorm * (1.0 - 1e-4 * alpha):
                        break                      # Armijo satisfied
                    if alpha * 0.5 <= 0.05:
                        break                      # floor: keep the smallest tried
                    alpha *= 0.5
                Y, R, Jbig, S = Y_trial, R_t, J_t, S_t
                if converged:
                    break
            if not converged:
                raise NoConvergenceError(
                    'Radau IIA(3) coupled stage Newton did not converge'
                    + ('' if not gshunt else ' at gshunt=%g S' % gshunt))
            ## the step end (the last stage) is what the epilogue reads
            if S[2] is not None:
                state_restore(S[2])
            return Y

        def _block_residual(Ylist):
            """`max |F_i|` for the coupled system, assembled from the SAME
            formula the solve uses: `F_i = q(Y_i) - q(x_n) - h sum_j A_ij K_j`
            with `K_j = -(i(Y_j) + u(t_j))`.

            ⚠ This exists because a FIXED-POINT test was not enough: if the
            solve hands its seed back, re-solving from that seed hands it back
            again and the fixed-point test passes vacuously.  A residual is a
            measurement; a fixed point of a broken solve is not.  So each
            stage's devices are read EXACTLY at it (`limit_sync`), not at
            whatever limiting state the last solve left, which is put back.
            """
            Ks = []
            snap = limiter_snapshot(lims) if lims else None
            try:
                for j, Yj in enumerate(Ylist):
                    if lims:
                        limit_sync(self.cir, Yj, epar, lims)
                    Ks.append(-(arr(self.cir.i(Yj, epar)) + src(tstage[j])))
            finally:
                if snap is not None:
                    state_restore(snap)
            worst = 0.0
            for i in range(3):
                Fi = arr(self.cir.q(Ylist[i], epar)) - qn \
                    - h * sum(Amat[i, j] * Ks[j] for j in range(3))
                worst = max(worst, float(np.max(np.abs(np.asarray(Fi)))))
            return worst

        seed0 = [self._pred_or(np.array(xn, dtype=float), tstage[i])
                 for i in range(3)]
        return ctx, _stage_newton, _block_residual, seed0

    def _rk_step_coupled(self, x0, t, provided_function=None):
        """One Radau IIA(3) step: the three collocation stages solved as ONE
        coupled ``3n`` Newton system.  Fully implicit -- no explicit first
        stage, no per-stage one-LU shortcut -- and stiffly accurate, so the
        step IS the last stage (``x_{n+1} = Y_3``, ``c_3 == 1``).  Self-starting
        (the only past state is ``x0``), so no history ring is read.

        The stage system in the charge formulation ``q'(x) = -(i(x) + u(t))`` is

            F_i(Y) = q(Y_i) - q(x_n) - h sum_j A_ij K_j = 0,   K_j = -(i(Y_j)+u(t_j))

        with block Jacobian ``J[i][j] = delta_ij C(Y_i) + h A_ij G(Y_j)`` --
        exactly ``(I3 (x) C + h A (x) G)`` in the linear case.  This dense
        coupled solve is the DEFAULT and the correctness reference.

        ⚠ THE COST TRANSFORM is the fast path, opt-in via the
        ``radau_transform`` Parameter.  It block-diagonalises the coupled
        system through ``eig(A^{-1})`` into one REAL and one COMPLEX `m x m`
        solve (see :meth:`_rk_step_transformed`) -- an
        ``O((3m)^3)`` dense solve becomes two sparse ones.  It is SIMPLIFIED
        Newton (one Jacobian per step), so on a strongly nonlinear step it can
        fail to converge; then this falls back to the dense full-Newton path
        below, so the answer is never wrong, only occasionally slower.

        When ``_rk_want_est`` is set (by the stepping loop's `_StageSteps` on
        an adaptive run) it also
        leaves the filtered embedded 5(3) error estimate in ``_rk_est`` (see
        :meth:`_radau_error_estimate`); the fixed-step path does not set the flag
        and pays nothing for it.  Returns ``(x, None, J, None)`` like
        `solve_timestep`, with ``J`` the last-stage operator, and leaves
        ``_iq``/``_q_cache`` set so the history push after the step is consistent.
        """
        ## ⚠ PCNR, WHERE ASKED FOR AND APPLICABLE, TAKES PRECEDENCE over the
        ## cost transform: it is a robustness the caller asked for, it has
        ## no transform variant, and run first the transform answered every
        ## step it converged on with device limiting, so `pcnr=True` did
        ## nothing and said nothing.
        use_pcnr = self._rk_use_pcnr()
        transform = self._newton_option(self.par.radau_transform,
                                        'radau_transform')
        ## (warned only where the transform was ASKED for: 'auto' gives way
        ## to PCNR silently)
        if self.par.radau_transform is True and use_pcnr and \
                not getattr(self, '_transform_pcnr_warned', False):
            self._transform_pcnr_warned = True
            warn(
                'Transient: radau_transform=True is not combined with '
                'pcnr=True -- PCNR has no transform variant, and it takes '
                'precedence: each step is solved by the dense coupled PCNR '
                'Newton.', UsageWarning)
        if transform and use_pcnr:
            _paths.COUNTS['radau.transform:pcnr'] += 1
        if transform and not use_pcnr:
            try:
                return self._rk_step_transformed(
                    x0, t, provided_function)
            except NoConvergenceError:
                ## simplified Newton stalled on this (nonlinear) step -- fall
                ## through to the dense full-Newton solve, which is the
                ## correctness reference and always converges here.
                _paths.COUNTS['radau.transform:fallback'] += 1
                self._radau_transform_fallbacks = getattr(
                    self, '_radau_transform_fallbacks', 0) + 1
        if use_pcnr:
            ## PCNR is the first-class limiting here too: the coupled solve keeps
            ## its structure but limits every junction, IN EVERY STAGE, by the
            ## joint continuation instead of per-device `cir.limit` -- which is
            ## the case device limiting cannot handle (parallel junctions on one
            ## branch fight over the shared voltage).  See _rk_step_coupled_pcnr.
            ##
            ## ⚠ AND IT FALLS BACK, exactly as the LMM step does (and DC before
            ## it): a PCNR failure on ONE step drops to the device-limiting
            ## coupled solve below rather than ending the transient.  That is
            ## what answers "what happens on a circuit bad enough to need the
            ## ladder?" -- the PCNR Newton has no continuation ladder of its own
            ## (a gshunt one and a junction-gmin one were both built and MEASURED
            ## not to rescue it: its bottleneck is the junction limiter's slew,
            ## which no deformation of the circuit accelerates), but the path it
            ## falls back to HAS one.  So the rescue is reached by falling back
            ## to the limiting that carries it, not by duplicating a ladder that
            ## does not work here.
            ##
            ## ⚠ THE TWO LIMITINGS AGREE AT THE ROOT (measured 0.0 / 5e-18 rel on
            ## a single junction), so a fallback step is not a different answer
            ## -- EXCEPT where PCNR was load-bearing: PARALLEL junctions on one
            ## branch, which per-device limiting resolves order-dependently.
            ## The warning says so, because that is the one case where a silent
            ## fallback would hand back a subtly different orbit.
            def _say(exc):
                _pairs = [(ra, rb) for _i, _e, ra, rb
                          in _pcnr.pcnr_junctions(self.cir)]
                _parallel = len(_pairs) != len(set(_pairs))
                self._pcnr_failed(
                    'coupled PCNR', t, exc,
                    ' -- ⚠ THIS CIRCUIT HAS PARALLEL JUNCTIONS ON ONE BRANCH, '
                    'which is the case PCNR exists for; the fallback resolves '
                    'them order-dependently' if _parallel else '')
            out, ok = self._pcnr_attempt(
                lambda: self._rk_step_coupled_pcnr(x0, t, provided_function),
                _say)
            if ok:
                return out
        ## (a new step for the device memo: this step's stages, and the
        ## previous step's, whose last is this step's start)
        self._memo_step()
        ctx, _stage_newton, _block_residual, seed0 = self._coupled_stage_solver(
            x0, t, provided_function)
        from pycircuit.circuit.nrsolver import _adaptive_conductance_ladder
        try:
            Y = _stage_newton(seed0)
        except NoConvergenceError:
            ## ⚠ THE CONTINUATION RESCUE REACHES THE COUPLED PATH TOO.  It is
            ## armed by `_solve` once the step has shrunk to `minstep`, and
            ## `_continuation_rescue` is read inside `self._newton`, which the
            ## coupled solve does not go through -- so it is read here too.
            ##
            ## Only the gshunt rung is offered: it is the one deformation that
            ## is structure-free (every node to ground through `g`, entering
            ## exactly where `G` does), so it needs no MNA row map and no
            ## per-stage junction bookkeeping.  Rungs march in exponent space
            ## and ONLY A PURE SOLVE IS RETURNED -- see
            ## `_adaptive_conductance_ladder`.
            ##
            ## ⚠ THE LINE SEARCH IS THE LAST RESORT, AFTER THE LADDER (owner
            ## decision).  As the FIRST retry it pre-empts the ladder and can
            ## beat it to a SPURIOUS root (a tanh charge circuit at
            ## h = 100 tau); the ladder tracks a physical branch and the
            ## search does not.  Where no ladder is available -- the PSS path
            ## never arms one, it drives `solve_timestep` on its own grid --
            ## there is nothing to pre-empt, so the search is the one retry
            ## there, opt-in through `_damped_last_resort`, which PSS sets; a
            ## plain adaptive transient keeps its step sequence unchanged.
            ## History: `doc/transient_history.md`, `Transient._rk_step_coupled`.
            if not getattr(self, '_continuation_rescue', False):
                if getattr(self, '_damped_last_resort', False):
                    Y = _stage_newton(seed0, damped=True)
                else:
                    raise
            else:
                try:
                    Y = _adaptive_conductance_ladder(
                        lambda seed, g: _stage_newton(seed, g),
                        lambda seed: _stage_newton(seed, 0.0),
                        seed0, label='coupled-stage gshunt stepping')
                except NoConvergenceError:
                    Y = _stage_newton(seed0, damped=True)
                self.statistics.gmin_rescues += 1
        return self._finish_radau(ctx, x0, t, provided_function, Y,
                                  solver=(_stage_newton, seed0, _block_residual))

    def _finish_radau(self, ctx, x0, t, provided_function, Y, solver=None):
        """The end every Radau IIA(3) coupled step shares, once its stages `Y`
        are solved and the devices sit at the last stage's limiting state: the
        branch check, the step end, the embedded error estimate.

        `solver` is the dense path's own ``(stage_newton, seed0, residual)``.
        The paths with their own Newton (PCNR, the transform) pass None: the
        branch check's confirmation then re-solves the same step equation with
        the dense coupled solver, built only if the screen fires."""
        ## (the passes the step's readers make at its stages, one core call
        ## a stage: `_stage_end_passes`, speed round 10, B3.2 -- on the paths
        ## with a Newton of their own; the dense Newton has recorded all
        ## four passes at every assembly's stages, the converged ones too)
        if solver is None:
            self._stage_end_passes(Y)

        ## BRANCH DETECTION on the COUPLED path.  ⚠ This path does NOT go
        ## through `self._newton`, so the check wired there does not reach
        ## the fully-implicit method -- which is the PSS default.  Same screen, same
        ## confirmation, but the solve being re-run is the coupled `3m`
        ## `_stage_newton` rather than a single-stage residual, so it needs its
        ## own call: a stage of the block can be at a rank drop while the
        ## others are not, and it is the BLOCK that has to be re-solved.
        ## History: `doc/transient_history.md`, `Transient._rk_step_coupled`.
        if solver is None:
            self._branch_after_coupled(
                None, None, Y, None,
                build=lambda: self._coupled_stage_solver(x0, t,
                                                         provided_function))
        else:
            self._branch_after_coupled(solver[0], solver[1], Y, solver[2])

        ## x_{n+1} == the last stage (stiff accuracy); the stage values are
        ## kept for the shooting monodromy, which needs all three
        xnp1, J = self._finish_stage_step(t, ctx.tstage, Y, ctx.Amat[2, 2],
                                          ctx.h, ctx.src)

        ## THE EMBEDDED 5(3) ERROR ESTIMATE (Hairer & Wanner Vol II, IV.8, the
        ## radau5 estimator), gated so the fixed-step path pays nothing.  See
        ## :meth:`_radau_error_estimate` for the construction and its gates.
        if getattr(self, '_rk_want_est', False):
            self._rk_est = self._radau_error_estimate(x0, Y, ctx.tn, ctx.h,
                                                      ctx.src, ctx.arr)
        return xnp1, None, J, None

    def _stage_end_passes(self, Y):
        """THE STAGES' READERS' PASSES, ONE CORE CALL A STAGE (speed round
        10, B3.2).  At the step's end ``x_{n+1} = Y[-1]`` the branch screen
        reads `C` and `_finish_stage_step` then evaluates `q`, `i` and `G`:
        two evaluations of one state, the device kernels run twice (on the
        PSP stage `C` 199 k instructions and `qiG` ~290 k; `qiCG` in one
        call 308 k: the kernels share their statements).  Here the passes
        the memo lacks there are made in one call (`_tran_core.passes`, bit
        for bit the circuit's own) and recorded (`_memo_put`), and each
        reader finds its own in the memo.  At the other stages the screen
        reads `C` and -- where the shooting walk factors the step at once
        (`_stage_G_read`, set around the step by `_walk_stage`) -- the
        shooting reads `G` (`_InnerTransient._G_at`): `CG` in one call
        there.  Only where the memo records (`_memo_ok`: a rolling memo, no
        stateful limiter, no bypass -- a recorded value is the one
        re-evaluating gives); where the core does not serve, each reader
        evaluates as before (counted)."""
        if not STAGE_FUSE:
            return _paths.no('radau.fuse:off')
        if not self._memo_ok():
            return _paths.no('radau.fuse:memo')
        get, put = self._memo_get, self._memo_put
        x = Y[-1]
        rec = get(x)
        need = ''.join(k for k in 'qiCG' if rec is None or k not in rec)
        if need:
            vals = _tran_core.passes(self, x, need)
            if vals is None:
                return _paths.no('radau.fuse:core')
            put(x, vals)
        if self.__dict__.get('_stage_G_read'):
            for y in Y[:-1]:
                rec = get(y)
                if rec is None or ('C' not in rec and 'G' not in rec):
                    vals = _tran_core.passes(self, y, 'CG')
                    if vals is None:
                        return _paths.no('radau.fuse:core')
                    put(y, vals)
        _paths.COUNTS['radau.fuse:served'] += 1

    def _radau_error_estimate(self, xn, Y, tn, h, src, arr):
        """The filtered embedded 5(3) error estimate for one Radau IIA(3) step
        (Hairer & Wanner Vol II, IV.8 -- the radau5 estimator), in STATE units.

        A lower-order (order 3) embedded solution ``yhat`` differs from the
        order-5 step by a combination of the three stage INCREMENTS
        ``Z_i = Y_i - x_n`` plus a fictitious explicit stage ``f(x_n)``:

            F1  = (dd1 Z1 + dd2 Z2 + dd3 Z3) / h            (state/time)
            rhs = C(x_n) F1 + f0,   f0 = -(i(x_n) + u(t_n)) (charge-rate)
            est = ((gamma_r/h) C(x_n) + G(x_n))^{-1} rhs     (state)

        with the radau5 constants ``dd1 = -(13+7 sqrt6)/3``,
        ``dd2 = (-13+7 sqrt6)/3``, ``dd3 = -1/3`` and ``gamma_r`` the real
        eigenvalue of ``A^{-1}`` (``RadauIIA3Integrator.GAMMA_REAL``).

        ⚠ THE FILTER IS WHAT MAKES IT STIFF-ROBUST, and it is the SAME real
        factor the cost transform uses.  An unfiltered embedded difference
        grows like ``|lambda h|`` on a stiff mode while the true error is
        damped to zero by L-stability, forcing the controller to crawl through
        exactly the transient the method exists to step over; the real-factor
        inverse maps it back to a bounded state error.  Solved in the
        reference-removed space (the full factor is singular on that row).

        ⚠ THE dd WEIGHTS WERE VALIDATED, NOT TRUSTED (roadmap 0j).  On the
        smooth RC the estimate scales as ``h^4`` (the order-3 embedded), and
        the adaptive controller it drives spends more steps and lands closer
        as ``reltol`` tightens -- both measured before this shipped.
        """
        import math
        integ = self.base_integrator
        gamma_r = integ.GAMMA_REAL
        s6 = math.sqrt(6.0)
        dd1 = -(13.0 + 7.0 * s6) / 3.0
        dd2 = (-13.0 + 7.0 * s6) / 3.0
        dd3 = -1.0 / 3.0
        epar = self.epar
        iref = self.irefnode
        tk = self.toolkit
        Y1, Y2, Y3 = Y
        xn = np.asarray(xn, dtype=float)
        Z1 = np.asarray(Y1, dtype=float) - xn
        Z2 = np.asarray(Y2, dtype=float) - xn
        Z3 = np.asarray(Y3, dtype=float) - xn
        F1 = (dd1 * Z1 + dd2 * Z2 + dd3 * Z3) / h
        ## ⚠ THE DEVICES AT `x_n`, EXACTLY: after the step a stateful
        ## limiter sits at the step END, and `i(x_n)` read there is the
        ## tangent at `x_{n+1}` -- so it is synced to `x_n` for these reads
        ## and put back (`limit_sync`).
        lims = stateful_limiters(self.cir)
        snap = limiter_snapshot(lims) if lims else None
        try:
            if lims:
                limit_sync(self.cir, xn, epar, lims)
            P = _tran_core.passes(self, xn, 'CGi')
            if P is not None:
                Cn, Gn = P['C'], P['G']
                f0 = -(P['i'] + src(tn))
            else:
                Cn = arr(self.cir.C(xn, epar))
                Gn = arr(self.cir.G(xn, epar))
                f0 = -(arr(self.cir.i(xn, epar)) + src(tn))
        finally:
            if snap is not None:
                state_restore(snap)
        rhs = np.asarray(Cn @ F1) + np.asarray(f0)
        real_factor = (gamma_r / h) * np.asarray(Cn) + np.asarray(Gn)
        (Rf,) = remove_row_col((real_factor,), iref, tk)
        rhs_r = tk.concatenate((rhs[:iref], rhs[iref + 1:]))
        est_r = self.toolkit.linearsolver(Rf, rhs_r)
        return self.toolkit.insert(est_r, iref, 0.0)

    def _radau_transform_matrices(self):
        """``(lam, V, Tinv)`` for the Radau cost transform, cached.

        ``A^{-1} = V diag(lam) V^{-1}`` with the eigenvalues ORDERED as
        ``[gamma_r (real), alpha + i beta, alpha - i beta]`` -- the real one
        first, then the pair with positive imaginary part, then its conjugate.
        Since ``A`` is real, ``V[:,0]`` is real and ``V[:,1] = conj(V[:,2])``,
        so the conjugate stage is free and the transform advances with one real
        and one complex ``m x m`` solve.
        """
        cached = getattr(self, '_radau_Tmats', None)
        if cached is not None:
            return cached
        A = np.array(self.base_integrator.A, dtype=float)
        Ainv = np.linalg.inv(A)
        lam, V = np.linalg.eig(Ainv)
        order = np.argsort(np.abs(lam.imag))
        ir = order[0]
        rest = [i for i in range(3) if i != ir]
        ip = rest[0] if lam[rest[0]].imag > 0 else rest[1]
        ic = rest[1] if ip == rest[0] else rest[0]
        idx = [ir, ip, ic]
        lam = lam[idx]
        V = V[:, idx]
        Tinv = np.linalg.inv(V)
        cached = (lam, V, Tinv)
        self._radau_Tmats = cached
        return cached

    def _radau_complex_solve(self, A, b):
        """Solve the complex system ``A x = b`` for the transform's complex
        stage, through :class:`ComplexKLUSolver` (cached) with a dense fallback
        when libklu is absent."""
        zs = getattr(self, '_radau_zsolver', 'unset')
        if zs == 'unset':
            try:
                from pycircuit.circuit.linearsolver import ComplexKLUSolver
                zs = ComplexKLUSolver()
            except ImportError:
                zs = None
            self._radau_zsolver = zs
        if zs is not None:
            return zs.solve(A, b)
        return np.linalg.solve(np.asarray(A, dtype=complex), np.asarray(b))

    def _radau_transform_solve(self, R3, Cr, Gr, h):
        """One transformed coupled solve: ``dY = (I3 (x) C + h A (x) G)^{-1}
        (-R3)`` via the ``A^{-1}`` eigenbasis, returning the reduced ``3m``
        update.

        Multiplying the coupled system by ``(A^{-1} (x) I)/h`` and
        diagonalising ``A^{-1} = V Lam V^{-1}`` decouples it into
        ``(lam_k C/h + G) w_k = rhs_k`` with
        ``rhs = -(Lam V^{-1} (x) I)/h R3``.  The real ``k=0`` block is a real
        solve; the ``k=1`` block a complex solve; ``k=2`` is its conjugate
        (free).  Then ``dY_i = V[i,0] w0 + 2 Re(V[i,1] w1)`` (real)."""
        lam, V, Tinv = self._radau_transform_matrices()
        m = Cr.shape[0]
        R3 = np.asarray(R3, dtype=float)
        F = (R3[0:m], R3[m:2 * m], R3[2 * m:3 * m])
        P = (np.diag(lam) @ Tinv) / h
        rhs0, rhs1 = _transform_rhs(F[0], F[1], F[2], (P[0, 0], P[0, 1], P[0, 2],
                                                       P[1, 0], P[1, 1], P[1, 2]))
        Cr = np.asarray(Cr, dtype=float)
        Gr = np.asarray(Gr, dtype=float)
        real_factor = (lam[0].real / h) * Cr + Gr
        w0 = np.asarray(self._get_linearsolver().solve(
            real_factor, np.real(rhs0), self.toolkit), dtype=float)
        comp_factor = (lam[1] / h) * Cr + Gr
        w1 = np.asarray(self._radau_complex_solve(comp_factor, rhs1),
                        dtype=complex)
        return _transform_back(w0, w1, [(V[i, 0].real, V[i, 1]) for i in range(3)])

    def _radau_frozen(self, Cr, Gr, h):
        """THE TRANSFORM'S FROZEN FACTORS, ONCE A STEP (speed round 10, B3.1):
        a `_FrozenTransform` whose `solve(R3)` is `_radau_transform_solve(R3,
        Cr, Gr, h)` -- or None (declined, counted), and each iteration calls
        that.

        The transform's Newton is SIMPLIFIED: `C` and `G` are frozen at `x_n`
        for the whole step, so `P`, the real and complex factors and their
        factorisations are the same at every iteration -- and
        `_radau_transform_solve` made them again each time: `P`, both
        factors, a dense LU of the real one (`numpy.linalg.solve`), a scipy
        CSC copy of the complex one (382 k instructions of a 578 k complex
        solve on the PSP stage) and its KLU refactor.  Here they are made
        once: the real factor's LU kept (`linearsolver.NumpyLU`, bit for bit
        numpy's solve) where the analysis solver's choice is the dense one,
        the complex factor marshalled once (`ComplexKLUSolver.prepare`) and
        refactored once (`solve_prepared`: a refactor of the same values is
        not repeated; the residual check and its fallback run on every
        solve).  Each iteration forms the right-hand sides and the update in
        `_radau_transform_solve`'s operations and order.  Made under
        `numpy.errstate(all='raise')`: an operation that would warn leaves
        the step to the per-iteration solve, which warns as before."""
        if not RADAU_FROZEN:
            return _paths.no('radau.frozen:off')
        chain = _FROZEN.get('chain') or _frozen_chain()
        if chain is None:
            return _paths.no('radau.frozen:patched')
        T = type(self)
        d = self.__dict__
        if ((T._radau_transform_solve, T._radau_complex_solve, T._radau_transform_matrices,
             T._get_linearsolver) != chain[:4]
                or '_radau_transform_solve' in d or '_radau_complex_solve' in d
                or '_radau_transform_matrices' in d or '_get_linearsolver' in d):
            return _paths.no('radau.frozen:patched')
        ZS, AS, DS, NTK, NUM = _FROZEN['types']
        lam, V, Tinv = self._radau_transform_matrices()
        Cr = np.asarray(Cr, dtype=float)
        Gr = np.asarray(Gr, dtype=float)
        try:
            with np.errstate(all='raise'):
                P = (np.diag(lam) @ Tinv) / h
                real_factor = (lam[0].real / h) * Cr + Gr
                comp_factor = (lam[1] / h) * Cr + Gr
        except FloatingPointError:
            return _paths.no('radau.frozen:fp')
        fz = _FrozenTransform()
        fz.m = Cr.shape[0]
        fz.p = (P[0, 0], P[0, 1], P[0, 2], P[1, 0], P[1, 1], P[1, 2])
        fz.v = tuple((V[i, 0].real, V[i, 1]) for i in range(3))
        fz.rf, fz.cf, fz.tk = real_factor, comp_factor, self.toolkit
        ls = fz.ls = self._get_linearsolver()
        fz.lu = None
        ## the real solve is numpy's where the solver is the dense one --
        ## `AutoSolver` once it chose it, or a `DenseSolver` -- on the numeric
        ## toolkit, each piece as defined
        if type(ls) is AS:
            dense = ls._choice
            ok = (AS.solve is chain[4] and AS._select is chain[5]
                  and 'solve' not in ls.__dict__ and '_select' not in ls.__dict__)
        else:
            dense, ok = ls, True
        if (ok and type(dense) is DS and DS.solve is chain[6] and 'solve' not in dense.__dict__
                and type(self.toolkit) is NTK and self.toolkit.linearsolver is chain[7]
                and NUM.linearsolver is chain[7] and np.linalg.solve is _NP_SOLVE
                and np.isfinite(real_factor).all()):
            from pycircuit.circuit.linearsolver import NumpyLU
            fz.lu = NumpyLU.make(real_factor)
        zs = getattr(self, '_radau_zsolver', 'unset')
        if zs == 'unset':
            ## (made as `_radau_complex_solve` makes it, at its first call)
            try:
                zs = ZS()
            except ImportError:
                zs = None
            self._radau_zsolver = zs
        fz.zs, fz.prep, fz.cs = zs, None, self._radau_complex_solve
        if (type(zs) is ZS and (ZS.solve, ZS.prepare, ZS.solve_prepared) == chain[8:11]
                and not ('solve' in zs.__dict__ or 'prepare' in zs.__dict__
                         or 'solve_prepared' in zs.__dict__)):
            fz.prep = zs.prepare(comp_factor)
        _paths.COUNTS['radau.frozen:served'] += 1
        return fz

    def _rk_step_transformed(self, x0, t, provided_function=None):
        """One Radau IIA(3) step by the COST TRANSFORM -- simplified Newton
        (one Jacobian, evaluated at ``x_n``) with each iteration's coupled
        solve done through the ``A^{-1}`` eigenbasis: one real and one complex
        ``m x m`` solve (:meth:`_radau_transform_solve`) instead of a dense
        ``3m`` factorisation.  The residual is the FULL per-stage residual, so
        the root is the same as the dense path; only the Jacobian is frozen.

        Raises ``NoConvergenceError`` if the frozen Jacobian does not carry the
        iteration to tolerance in ``maxiter`` steps -- the caller
        (:meth:`_rk_step_coupled`) then falls back to the dense
        full-Newton solve.  Sets the same downstream state (``_rk_Y``,
        ``_iq``, ``_q_cache``, ...) so every consumer is identical to the dense
        path."""
        ctx = self._coupled_stage_context(x0, t, provided_function)
        h, iref = ctx.h, ctx.iref
        arr, src, tstage = ctx.arr, ctx.src, ctx.tstage
        epar = self.epar
        tk = self.toolkit
        xn = x0

        ## this step's device memo (the dense path rolls it too, after the
        ## transform's dispatch): the frozen `C` and `G` at `x_n` are the
        ## previous step's final stage's, recorded by `_finish_stage_step`
        self._memo_step()
        _rec = self._memo_get(xn) or {}
        need = [k for k in ('C', 'G') if k not in _rec]
        ## the FROZEN Jacobian pieces, at x_n, reduced -- the whole point:
        ## factored implicitly once per step and reused every iteration
        P = _tran_core.passes(self, xn, ''.join(need)) if need else None
        if P is not None:
            Cn = _rec['C'] if 'C' in _rec else P['C']
            Gn = _rec['G'] if 'G' in _rec else P['G']
        else:
            with _evalhint.evaluating(*need):
                Cn = _rec['C'] if 'C' in _rec else arr(self.cir.C(xn, epar))
                Gn = _rec['G'] if 'G' in _rec else arr(self.cir.G(xn, epar))
        if need and self._memo_ok():
            self._memo_put(xn, {'C': Cn, 'G': Gn})
        (Cr,) = remove_row_col((Cn,), iref, tk)
        (Gr,) = remove_row_col((Gn,), iref, tk)
        Cr = np.asarray(Cr, dtype=float)
        Gr = np.asarray(Gr, dtype=float)
        ## (the transform's factors, once a step: `_radau_frozen`, speed round
        ## 10, B3.1; None -- each iteration makes them)
        fz = self._radau_frozen(Cr, Gr, h)

        ## STAGE PREDICTOR.  ⚠ This path is SIMPLIFIED Newton (one Jacobian per
        ## step), so it is the one that benefits most from starting near the
        ## root and the one that most needs the clamp: it has no fresh Jacobian
        ## to recover from a wild guess, only the dense fallback.
        Y = [self._pred_or(np.array(xn, dtype=float), tstage[i])
             for i in range(3)]
        reltol = self.par.reltol
        abstol = float(self.par.vabstol)
        maxit = int(self.par.maxiter)
        ## ⚠ EACH STAGE OWNS ITS LIMITING STATE, and the step converges only
        ## with every limiter at rest -- as the dense path, and for its
        ## reasons (`_coupled_stage_solver`): shared, the transform landed
        ## 1.09e-4 V off the exact solution on a diode driven to 0.85 V.
        lims = stateful_limiters(self.cir)
        S = self._stage_limiter_states(Y, lims, epar)
        ## (the loop in one C call where it serves -- its iterates, the
        ## source memo and the two factorisations' state as the loop leaves
        ## them: `_tran_radau_tc`, speed round 10, B3.3; None -- declined or
        ## handed back, and the loop runs from the same seed)
        _nobypass = float(getattr(epar, 'bypasstol', -1.0) or -1.0) < 0.0
        Yc = _tran_radau_tc.solve(self, ctx, fz, Y, src, provided_function, lims,
                                  _nobypass, reltol, abstol, maxit)
        if Yc is not None:
            Y = Yc
        elif not self._transform_loop(ctx, fz, Y, S, lims, Cr, Gr, reltol, abstol, maxit):
            raise NoConvergenceError(
                'Radau IIA(3) transform (simplified Newton) did not converge')
        ## the step end (the last stage) is what the epilogue reads
        if S[2] is not None:
            state_restore(S[2])
        return self._finish_radau(ctx, x0, t, provided_function, Y)

    def _transform_loop(self, ctx, fz, Y, S, lims, Cr, Gr, reltol, abstol, maxit):
        """The transform's simplified Newton from the stages `Y` (updated in
        place, as their limiter states `S`): True once converged within
        `maxit` iterations.  Per iteration each stage's `q` and `K` at its
        own limiting state, the residual, the transform's solve (the step's
        frozen factors, `fz`; None: made each iteration), and per stage the
        step, the limiter and its snapshot; the test after the update.
        `_tran_radau_tc` runs it in C where it serves."""
        h, iref = ctx.h, ctx.iref
        arr, red, src, tstage, m = (ctx.arr, ctx.red, ctx.src, ctx.tstage,
                                    ctx.m)
        epar = self.epar
        tk = self.toolkit
        for _ in range(maxit):
            ## each stage's q and K read ONCE, at its own limiting state
            qi_all, Ki_all = [], []
            for j in range(3):
                if S[j] is not None:
                    state_restore(S[j])
                P = _tran_core.passes(self, Y[j], 'qi')
                if P is not None:
                    qi_all.append(P['q'])
                    Ki_all.append(-(P['i'] + src(tstage[j])))
                else:
                    with _evalhint.evaluating('q', 'i'):
                        qi_all.append(arr(self.cir.q(Y[j], epar)))
                        Ki_all.append(-(arr(self.cir.i(Y[j], epar))
                                        + src(tstage[j])))
            R, _J = self._coupled_stage_system(ctx, qi_all, Ki_all)
            dY = fz.solve(R) if fz is not None else self._radau_transform_solve(R, Cr, Gr, h)
            scale = 0.0
            for i in range(3):
                if S[i] is not None:
                    state_restore(S[i])
                di = dY[i * m:(i + 1) * m]
                Y_prev = Y[i]
                Y_trial = Y_prev + tk.insert(di, iref, 0.0)
                Y_new = self.cir.limit(Y_trial, Y_prev, epar)
                Y[i] = Y_new
                if lims:
                    S[i] = limiter_snapshot(lims)
                step_i = red(np.asarray(Y_new) - np.asarray(Y_prev))
                scale = max(scale, np.max(np.abs(step_i)))
            if self._stages_converged(Y, S, lims, scale, reltol, abstol, red):
                return True
        return False
