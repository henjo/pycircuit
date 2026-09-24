"""States defined up to a modulus (idtmod): the state fold of the closing
residual and the jump across an output wrap.
"""
import numpy as np
import warnings


class _PeriodicStates(object):
    """States defined up to a modulus (idtmod): the state fold of the
    closing residual and the jump across an output wrap.  A theme of `PSS`
    (see `pss.py`)."""

    ## No fold until a solve collects one -- `_fold_periodic` is then the
    ## identity, which is exactly right for every circuit without a folding
    ## state and for any caller that reaches a residual outside `solve`.
    _periodic_fold = []

    def _collect_periodic_fold(self):
        """`[(reduced_row, modulus, offset)]` for every state defined up to
        `n*modulus`.

        `Circuit.periodic_states()` reports GLOBAL rows; the shooting residual
        lives in reduced coordinates, so each is mapped through the same
        `irefnode` deletion the stamps use.  A declared row that IS the
        reference node has no reduced coordinate and is dropped.

        The state row's own fold needs only the modulus -- a DIFFERENCE of
        two states is defined up to `n*modulus` whatever window each sits in.
        The offset is for `_wrap_jump`: the declared window's edges
        (``offset + k*modulus``) are where the element's OUTPUT wraps (the
        `periodic_states` contract), so it says on which branch of the
        output a state sits.
        """
        cir = self.cir
        if not hasattr(cir, 'periodic_states'):
            return []
        try:
            declared = cir.periodic_states()
        except Exception:
            ## A late-bound modulus that cannot be resolved is not a reason to
            ## refuse the solve -- it is a reason not to fold.
            return []
        iref = self.irefnode
        out = []
        for row, modulus, offset in declared or []:
            row = int(row)
            m = float(modulus)
            if row == iref or not np.isfinite(m) or m <= 0.0:
                continue
            out.append((row if row < iref else row - 1, m, float(offset)))
        return out

    def _fold_periodic(self, F):
        """Fold a shooting residual's periodic rows into `[-m/2, m/2)`.

        ⚠ THIS IS WHAT MAKES THE PERIOD MAP'S FIXED-POINT PROBLEM WELL POSED
        FOR A FOLDING STATE, and it belongs in the RESIDUAL rather than in the
        traversal.  `Idtmod`'s state is defined only up to `n*modulus` -- it
        says so itself through `periodic_states()`, and the transient engine
        already uses that declaration to keep the state bounded by exact gauge
        translations.  Shooting did not: it asked for `x_0 - phi(x_0) == 0`
        literally, which on a folding row demands the SAME REPRESENTATIVE, not
        the same state.  An orbit that closes after advancing exactly one
        modulus -- the normal case for a phase -- then has NO root at all, and
        near the fold the raw difference jumps by a whole modulus while the
        state moves infinitesimally.

        Measured (see `test_a_state_fold_breaks_the_period_map_at_the_ENDPOINT_
        not_on_the_grid`): that jump is grid-INDEPENDENT -- 1.414214e+09 at
        seven different grids, with the wrap on a node and off it alike -- so
        it is not an event-localisation defect and no refinement of the time
        grid can reach it.  It is the output map, and the output map is what
        this folds.

        ⚠ THE JACOBIAN IS DELIBERATELY NOT TOUCHED.  `d/dx0` of
        `wrap(x0 - phi(x0))` equals `d/dx0 (x0 - phi(x0))` almost everywhere --
        the wrap has unit slope between its jumps -- so `D - alpha*Mx` is
        already the right derivative of the folded residual.  The fold moves
        the residual onto the branch the Jacobian was always describing; that
        is the whole reason this is a local change and not surgery on the six
        traversal loops.
        """
        rows = self._periodic_fold
        if not rows:
            return F
        F = np.asarray(F, dtype=float).copy()
        ## every STATE of the unknown folds on its own rows -- gear's pair
        ## `(x_0, x_{-1})` carries two (its autonomous and event residuals did
        ## not fold at all before 2026-09-23; the driven one folded each half)
        w = self.cir.n - 1
        for off in range(0, F.shape[0], w):
            for r, m, _o in rows:
                if off + r < F.shape[0]:
                    F[off + r] -= m * np.round(F[off + r] / m)
        return F

    def _close_periodic(self, z0, z_end, tms):
        """The closing residual ``z_0 - z_end``, folded: the orbit's jump
        across an output wrap removed (`_wrap_jump`), then every periodic
        state row folded (`_fold_periodic`).  `tms` is the period's grid;
        block `b` of the unknown (gear's pair has two) ends at
        ``tms[-1 - b]``."""
        z0 = np.asarray(z0, dtype=float)
        z_end = np.asarray(z_end, dtype=float)
        return self._fold_periodic(z0 - z_end - self._wrap_jump(z0, z_end, tms))

    def _wrap_jump(self, z0, z_end, tms):
        """The orbit's jump across an output wrap that lies between `z_end`
        and `z_0`, per block of the unknown; zero when none does.

        ⚠ THE STATE FOLD ALONE LEAVES THE OUTPUTS DISCONTINUOUS.  An
        idtmod's OUTPUT is `wrap(state)`, and every algebraic quantity it
        feeds (the phase node, the branch current into a load, a phase
        detector's node) jumps with it.  When the orbit wraps exactly at
        ``t = 0`` -- a free-running VCO pinned at a zero crossing of
        `sin(2 pi phase)`, or a reference phase that starts at 0 -- `z_0`
        sits just after the wrap and `z_end` just before it, the state row
        folds to zero, and those rows still differ by the whole jump.  At
        rounding level the Newton iterates straddle the wrap and the jump
        flips in and out of the residual: radau on the free-running
        `VcoHdl` failed there (2026-09-23), with every other row at 1e-14.

        ⚠ WHICH SIDE A STATE IS ON IS DECIDED BY THE STATE, NEVER BY THE
        RESIDUAL.  The declared window's edges are where the output wraps,
        so the fractional window positions `f_0`, `f_end` of the two states
        say whether the state fold's nearest-representative path crosses an
        edge: ``k = -round(f_0 - f_end)``, which is -1, 0 or +1.  Folding an
        output row by the modulus instead (round the residual) accepts a
        start point one modulus off its own constraint -- measured: radau,
        gear and trap all "converged" with `ph(0) = 1.1` against a phase of
        0.1 -- and cannot reach a loaded output at all (1 mA on a load
        resistor's current, 2 V on a gain-2 detector's node).

        The jump itself is the algebraic part re-solved on each side of the
        edge with the differential part held: ``x = x_end' + N y`` with
        ``N = ker C``, solving ``Z^T (i(x) + u(t)) = 0`` with ``Z = ker
        C^T``, once with the straddling states just past their edges on
        `z_0`'s side and once on `z_end`'s.  The difference ``N (y_+ - y_-)``
        is the jump, exact for a nonlinear load as well; it lies in `ker C`,
        so the states' own rows never see it.  The Jacobian is untouched:
        `k` is constant between straddles.  ⚠ A capacitor ON a wrapped node
        makes the jump an impulse (index 2); the algebraic block is then
        singular, no jump is subtracted, and the solve warns once.
        """
        z0 = np.asarray(z0, dtype=float)
        z_end = np.asarray(z_end, dtype=float)
        jump = np.zeros_like(z_end)
        rows = self._periodic_fold
        if not rows:
            return jump
        w = self.cir.n - 1
        for b, off in enumerate(range(0, z_end.shape[0], w)):
            xe = z_end[off:off + w]
            sides = []
            for r, m, o in rows:
                f0 = (z0[off + r] - o) / m
                fe = (xe[r] - o) / m
                k = -np.round((f0 - np.floor(f0)) - (fe - np.floor(fe)))
                if k:
                    sides.append((r, m, o, int(k)))
            if sides:
                jump[off:off + w] = self._jump_across(
                    xe, float(tms[len(tms) - 1 - b]), sides)
        return jump

    def _jump_across(self, x, t, sides):
        """The algebraic jump at `x` (time `t`) when each state in `sides`
        (``(row, modulus, offset, k)``) crosses the output edge it sits
        next to, from its own side to the other.  See `_wrap_jump`."""
        plus, minus = x.copy(), x.copy()
        for r, m, o, k in sides:
            e = o + m * np.round((x[r] - o) / m)
            d = 64.0 * np.finfo(float).eps * max(abs(e), m)
            plus[r], minus[r] = e + k * d, e - k * d
        ## ker C and ker C^T from an EQUILIBRATED C: a femtofarad column next
        ## to an idtmod state's unit charge must not read as rank-deficient
        C = np.asarray(self._C_at(x), dtype=float)
        cs, rs = np.max(np.abs(C), axis=0), np.max(np.abs(C), axis=1)
        dc = np.where(cs > 0.0, 1.0 / np.where(cs > 0.0, cs, 1.0), 1.0)
        dr = np.where(rs > 0.0, 1.0 / np.where(rs > 0.0, rs, 1.0), 1.0)
        U, S, Vt = np.linalg.svd(dr[:, None] * C * dc[None, :])
        rank = int(np.sum(S > 1e-10 * (S[0] if S.size and S[0] > 0 else 1.0)))
        N = dc[:, None] * Vt[rank:].T
        Z = dr[:, None] * U[:, rank:]
        if N.shape[1] == 0:
            return np.zeros_like(x)

        def algebraic(xs):
            y = np.zeros(N.shape[1])
            for _it in range(30):
                xk = xs + N @ y
                g = Z.T @ np.asarray(self._k_at(xk, t), dtype=float)
                dy = np.linalg.solve(Z.T @ np.asarray(self._G_at(xk),
                                                      dtype=float) @ N, g)
                y = y + dy
                if np.max(np.abs(N @ dy)) <= 1e-13 * max(1.0, np.max(np.abs(xk))):
                    return y
            raise np.linalg.LinAlgError('the algebraic re-solve did not converge')

        try:
            return N @ (algebraic(plus) - algebraic(minus))
        except np.linalg.LinAlgError:
            if not getattr(self, '_wrap_jump_warned', False):
                self._wrap_jump_warned = True
                warnings.warn(
                    'PSS: the orbit starts on an idtmod output wrap, and the '
                    'jump across it could not be resolved (the algebraic '
                    'block is singular there -- a capacitor on a wrapped '
                    'node?); the solve may not converge.', RuntimeWarning,
                    stacklevel=2)
            return np.zeros_like(x)
