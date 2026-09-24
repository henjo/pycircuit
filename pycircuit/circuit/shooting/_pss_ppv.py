"""The adjoint side of the orbit: the PPV, the frequency-aware PPV, the Floquet
modes and the continuous adjoint.
"""
import numpy as np
import warnings
from pycircuit.circuit.analysis import remove_row_col
from ._numerics import _arnoldi_gmres
from .events import EventColumns


class _PPVFloquet(object):
    """The adjoint side of the orbit: the PPV, the frequency-aware PPV, the
    Floquet modes and the continuous adjoint.  A theme of `PSS` (see
    `pss.py`)."""

    ## How many deflated power iterations estimate the second multiplier.
    ## Convergence is at |lambda_3|/|lambda_2|, which is fast in the case
    ## that matters -- a lone slow node leaves everything below it tiny.
    ## ⚠ RETIRED 2026-09-03, kept as a name so the history reads.  The
    ## deflated power iteration this sized converged at
    ## `|lambda_3|/|lambda_2|` and lost three digits at a ratio of 1.065;
    ## Arnoldi replaced it at machine precision and fewer matvecs.
    PPV_DEFLATION_ITERS = 30
    ## Arnoldi basis size for the second-multiplier estimate.  Exact at
    ## `k = n`; below that it is a truncation and `lam2` is a lower bound
    ## ONLY for a normal `M` -- see `ppv`, where a circuit monodromy was
    ## measured to break the bound in both directions.  This is now the
    ## STARTING basis: `_ritz_second_multiplier` grows it until the pair's
    ## own Ritz residual certifies it.
    PPV_RITZ_BASIS = 12

    ## ⚠ THE GATE ON A TRUNCATED `lam2`: the SELECTED PAIR's Ritz residual,
    ## `|h_{k+1,k}| |y_i[last]| / ||y_i||`.  It needs no extra matvec -- both
    ## factors are already in the `H` this class forms -- and it is the one
    ## quantity that separates a converged pair from a leaked one.  The
    ## GMRES-style residual cannot: `_arnoldi_gmres`'s own note says a
    ## drifted basis "gives multipliers that are wrong in a way the residual
    ## cannot see", which is true of the SOLVE residual and false of the
    ## EIGENPAIR one.
    ##
    ## MEASURED (`_osc_with_ladder`, k = 12): <= 3.1e-07 at every `nslow`
    ## the truncated path gets right, and 2.8e-04 / 3.5e-04 / 2.1e-03 at the
    ## three it gets wrong -- and 1.5e-16 once `k` is large enough to be
    ## exact.  A peer session's independent sweep puts the two populations
    ## thirteen decades apart at the median, ⚠ TOUCHING at ~1e-5 (right 90th
    ## percentile 1.27e-05 against wrong 10th percentile 1.03e-05).  So the
    ## robust band is BELOW ~1e-6, and that -- not a magic 1e-8 -- is what
    ## this is set to.
    PPV_RITZ_RESIDUAL_TOL = 1e-6

    ## ⚠ A COST CEILING, NOT A CORRECTNESS THRESHOLD.  `k` doubles until the
    ## residual certifies or this is reached; overrunning it produces a
    ## WARNING and an uncertified number, never a silently wrong one, which
    ## is what makes an arbitrary-ish constant safe here.
    ##
    ## The size is set by what `k` has to reach: measured, `k` tracks the
    ## SLOW-MODE COUNT and not `n` (this tree's ladder cannot separate the
    ## two -- it sets `nslow = nladder` -- but a peer's synthetic can, and
    ## reports `k_min` flat under a doubling of `n` at fixed `nslow`).  The
    ## largest published case is Lai's 64-gated-capacitor DCO, so 128 leaves
    ## 2x headroom on it while costing 1/6 of that circuit's `n = 813`.
    PPV_RITZ_MAX_BASIS = 128
    ## Above this, the bordered extraction is losing digits AND the phase
    ## equation's instantaneous-response assumption is in doubt.
    PPV_SECOND_MULTIPLIER_WARN = 0.9

    def _ritz_second_multiplier(self, fp, kk):
        """`(lam2, residual)` from a `kk`-dimensional Arnoldi on `I - M`.

        Garcia, Romero & Acha (IEEE Trans. Power Systems 37(1), 2022): the
        Ritz values of `I - M` map back as `lambda = 1 - theta`.  Returns the
        selected pair's own RITZ RESIDUAL alongside it,
        `|h_{k+1,k}| |y_i[last]| / ||y_i||`, which costs nothing -- both
        factors are already in `H` and its eigenvectors.

        ⚠ `eig`, NOT `eigvals`.  The eigenVECTOR's last component is half the
        residual, so asking only for the values is what made this estimate
        uncheckable for as long as it was.

        ⚠ EVERY LOCAL IS UNDERSCORED ON PURPOSE, inherited from when this was
        inline in `ppv`: the first version used `q` for the Arnoldi start
        vector, silently overwriting `C(0) xdot(0)` -- returned as
        `info['q']` and consumed by two tests -- and the suite caught it as a
        shape mismatch three frames away.

        `residual` is `inf` when no Ritz value survives the deflation, so a
        caller that gates on it cannot read "nothing found" as "certified".
        """
        _n = fp.width
        kk = int(min(_n, kk))
        _rng = np.random.default_rng(12345)
        _q0 = _rng.standard_normal(_n)
        _q0 = _q0 / np.linalg.norm(_q0)
        _Qb = [_q0]
        _H = np.zeros((kk + 1, kk))
        for _j in range(kk):
            _wj = _Qb[_j] - np.asarray(fp.matvec(_Qb[_j]))
            for _i in range(_j + 1):
                _H[_i, _j] = float(_Qb[_i] @ _wj)
                _wj = _wj - _H[_i, _j] * _Qb[_i]
            _H[_j + 1, _j] = float(np.linalg.norm(_wj))
            if _H[_j + 1, _j] < 1e-13:
                ## ⚠ AN INVARIANT SUBSPACE: the basis closed on itself, so
                ## every Ritz pair in it is EXACT and the residual is zero by
                ## construction rather than by convergence.
                kk = _j + 1
                _theta, _Y = np.linalg.eig(_H[:kk, :kk])
                _lams = 1.0 - _theta
                _mask = np.abs(_lams - 1.0) > 1e-6
                if not _mask.any():
                    return 0.0, float('inf')
                return float(max(np.max(np.real(_lams)[_mask]), 0.0)), 0.0
            _Qb.append(_wj / _H[_j + 1, _j])
        _theta, _Y = np.linalg.eig(_H[:kk, :kk])
        _lams = 1.0 - _theta
        ## drop the unit root the border already accounts for
        _mask = np.abs(_lams - 1.0) > 1e-6
        if not _mask.any():
            return 0.0, float('inf')
        _idx = int(np.where(_mask)[0][int(np.argmax(np.real(_lams)[_mask]))])
        _lam2 = float(max(np.real(_lams)[_idx], 0.0))
        _res = (abs(_H[kk, kk - 1]) * abs(_Y[kk - 1, _idx])
                / max(float(np.linalg.norm(_Y[:, _idx])), 1e-300))
        return _lam2, float(_res)

    def _algebraic_adjoint_pattern(self, xf):
        """`(rows, cols)` — the ALGEBRAIC equations and the algebraic states.

        `rows` are the equations with no time-derivative (a zero row of
        `C`); `cols` are the state variables that appear under no
        derivative anywhere (a zero column).  For an index-1 DAE the two
        have the same count and the block between them is invertible,
        which is what makes the fill below well posed.

        Returns `([], [])` for a plain ODE, and the caller then does no
        per-sample work at all -- which is why this costs nothing on every
        fixture that has no algebraic row.
        """
        m = self.cir.n - 1
        Cr, = remove_row_col((np.asarray(self.cir.C(xf), dtype=float),),
                             self.irefnode, self.toolkit)
        Cr = np.asarray(Cr, dtype=float)
        rows = [i for i in range(m) if not np.any(Cr[i, :])]
        cols = [j for j in range(m) if not np.any(Cr[:, j])]
        return rows, cols

    def _equation_row_ppv(self, vblock, xf, rows, cols):
        """`v_1` — the adjoint contracted with an EQUATION-ROW input.

        ⚠⚠ THERE ARE TWO ADJOINT VECTORS AND CONFLATING THEM WAS A DEFECT.
        Demir 2000 puts both conventions on one page: a STATE initial
        condition contracts as `v_1^T(0) C(0) x(0)` (eq 41, WITH `C`), while
        an EQUATION-ROW input contracts as `v_1^T(s) b(s)` (eq 42, and the
        phase equation 44, BARE).  `ppv()` returns `C^T v_1`, which is the
        right object for a state perturbation and is documented as such.
        `CY` is an equation-row covariance -- a current injected into a KCL
        row -- so `diffusion_constant` and `colour_projection` need `v_1`,
        and this produces it.

        ⚠ ONE DEFECT, TWO SYMPTOMS, both measured.  `(C^T v_1)_j` is COLUMN
        `j` of `C` dotted with `v_1`, so on DIFFERENTIAL states `C^T`
        multiplies by the capacitance -- `diffusion_constant` was wrong by
        `C^2`, ratios 0.010003 / 1.000334 / 100.033536 over a 100x sweep --
        and on ALGEBRAIC states `C^T` ANNIHILATES, so `c` came back EXACTLY
        0.0 for an oscillator whose only noise was its series tank loss.
        Neither was visible because every fixture here uses `C = 1 F`.

        ⚠⚠ AND `C` IS NEVER INVERTED.  The solve is on `C[D, NZ]^T`, the
        block between the DIFFERENTIAL equations and the NON-algebraic
        states, which is square and invertible by construction; the
        singular `C` as a whole is not touched.  The algebraic entries come
        from the constraint below instead.
        """
        m = self.cir.n - 1
        Cr, = remove_row_col((np.asarray(self.cir.C(xf), dtype=float),),
                             self.irefnode, self.toolkit)
        Cr = np.asarray(Cr, dtype=float)
        diff = [i for i in range(m) if i not in rows]
        nz = [j for j in range(m) if j not in cols]
        out = np.zeros(m, dtype=float)
        blkC = Cr[np.ix_(diff, nz)].T
        try:
            out[diff] = np.linalg.solve(blkC, np.asarray(vblock)[nz])
        except np.linalg.LinAlgError:
            warnings.warn(
                'PSS.ppv: the differential block C[D, NZ] is singular, so '
                'the equation-row adjoint cannot be recovered; falling back '
                'to the state-perturbation vector, which is wrong by a '
                'factor of the capacitance in any CY contraction.',
                RuntimeWarning, stacklevel=2)
            return np.array(vblock, dtype=float, copy=True)
        if not rows:
            return out
        return self._algebraic_adjoint_fill(out, xf, rows, cols)

    def _algebraic_adjoint_fill(self, vblock, xf, rows, cols):
        """Fill the adjoint's ALGEBRAIC entries, which are SLAVED, not free.

        ⚠ THE REPLAY LEAVES THEM AT ZERO AND ZERO IS NOT THEIR VALUE.  The
        PPV entry for a row IS the phase sensitivity to a perturbation
        entering that row, and an algebraic row's perturbation reaches the
        dynamics through the CONSTRAINT rather than through its own row.
        On a tank with series loss, eliminating `v_x = r (i_L + b)` puts
        `-r b` into the inductor's equation, so the sensitivity to node
        `x` is `r` times the branch row's -- and a noise current landing
        there was being contracted against a structural zero, which made
        `diffusion_constant` return EXACTLY 0.0 for an oscillator whose
        only noise was its tank loss.  Measured against three independent
        references; see the roadmap's section 0d.

        The adjoint equation's ALGEBRAIC-STATE columns are what determines
        them.  For a column `j` with no `C` entry the equation carries no
        derivative, so it reads `sum_i G_ij v_i = 0`, and splitting `i`
        into algebraic and differential rows gives

            v_A  =  (G[A, Z]^T)^-1 G[D, Z]^T v_D

        ⚠ THE MAGNITUDE IS STRUCTURAL AND THE SIGN IS MEASURED, and saying
        which is which is the point.  On a RESISTIVE DIVIDER between the
        inductor and ground the two algebraic nodes fold into the branch
        row with coefficients `(r1 + r2)` and `r2`, so their entries must
        stand in a ratio the topology fixes -- measured `10.000000` against
        a chosen `10.000000`, on a fixture built so the single-resistor
        degeneracy cannot hide a mistake.  ⚠ THE SINGLE SERIES RESISTOR
        CANNOT SETTLE THIS: there `|integral v_0|` and `|r integral
        v_branch|` agree to 1.5e-4, so BOTH SIGNS FIT and an agreement
        there is no evidence at all.

        The overall sign is then fixed by requiring algebraic and
        differential rows to share ONE convention -- `dT/dA = +integral
        v_j` -- and measured on the divider: with it, the DC-injection
        probe reads +0.9999849 (node v, differential), +0.9999982 and
        +0.9998744 (the two algebraic nodes); with the sign the naive
        derivation gives, the last two come back NEGATIVE.

        Returns `vblock` unchanged, with a warning, when the structure is
        not index-1: `len(rows) != len(cols)` or a singular block.  That
        case is section B4's, and guessing at it would be worse than
        leaving a known zero.
        """
        if not rows:
            return vblock
        m = self.cir.n - 1
        if len(rows) != len(cols):
            warnings.warn(
                'PSS.ppv: %d algebraic equations against %d algebraic '
                'states, so the adjoint\'s algebraic entries are not '
                'determined by a square solve -- this is an index > 1 '
                'structure (roadmap B4). They are left at zero, and noise '
                'entering those rows will be UNDER-COUNTED.'
                % (len(rows), len(cols)), RuntimeWarning, stacklevel=2)
            return vblock
        Gr, = remove_row_col((np.asarray(self.cir.G(xf), dtype=float),),
                             self.irefnode, self.toolkit)
        Gr = np.asarray(Gr, dtype=float)
        diff = [i for i in range(m) if i not in rows]
        blk = Gr[np.ix_(rows, cols)].T
        ## ⚠ THE SIGN IS NOW THE DERIVED ONE, because this acts on `v_1`.
        ## It was FLIPPED while this fill acted on `C^T v_1`: on that
        ## fixture `C = diag(1, 0, -L)`, the INDUCTOR BRANCH ROW CARRIES
        ## `-L`, and the term below is dominated by that branch -- so
        ## reading `C^T v_1` as `v_1` negated it.  The derivation and the
        ## measurement were describing DIFFERENT VECTORS and both were
        ## right.  Pinned by the eq (24) constraint, which returns
        ## 0.0000e+00 exactly here and 1.9870e+00 for either alternative.
        rhs = -(Gr[np.ix_(diff, cols)].T @ np.asarray(vblock)[diff])
        try:
            va = np.linalg.solve(blk, rhs)
        except np.linalg.LinAlgError:
            warnings.warn(
                'PSS.ppv: the algebraic block G[A, Z] is singular, so the '
                'adjoint\'s algebraic entries cannot be recovered; they '
                'are left at zero and noise entering those rows will be '
                'UNDER-COUNTED.', RuntimeWarning, stacklevel=2)
            return vblock
        out = np.array(vblock, dtype=float, copy=True)
        out[rows] = va
        return out

    def _ppv_propagate(self, fp, v, m, xdot, alg_rows, alg_cols, inject=None):
        """The pair-consistent SECOND-ORDER propagation of an anchor vector
        `v` over the period (the block `ppv()` applies to its null vector,
        lifted 2026-09-08 so `frequency_aware_ppv` can run it on a COMPLEX
        anchor).  Returns `(states, states_pair, ts, Xf)`: `states` the
        per-step state-space samples (`C^T v`, second order, rescaled by
        `v . xdot = 1` at the first sample), `states_pair` the raw
        pair-space replay, `ts` the transposed per-step states, `Xf` the
        orbit.  Only the `solved_history` (LMM) kind carries the
        correction; the others return the replay as it is."""
        _dt = complex if np.iscomplexobj(v) else float
        _alg_rows, _alg_cols = alg_rows, alg_cols
        self._ppv_alg_fallback = False
        _end, _ts, states = fp.matvec_transposed(v, collect=True, inject=inject)
        states_pair = [np.array(st, dtype=_dt, copy=True) for st in states]
        _Xf = np.asarray(self.waveform[1], dtype=float)
        ## ⚠⚠ THE PAIR'S FIRST BLOCK IS NOT THE PPV, AND THE ERROR IS FIRST
        ## ORDER AND GROWS WITH Q.  For Gear-2 the adjoint state is the pair
        ## `(w1, w2) = (dphi/dx_k, dphi/dx_{k-1})`, and `w1` alone is the
        ## response to a perturbation of `x_k` WITH `x_{k-1}` HELD -- an
        ## inconsistent history, which the two-step method resolves through
        ## its parasitic root.  A physical state perturbation moves both:
        ## `dx_{k-1} = Phi(t_{k-1}, t_k) dx_k`, so the phase functional is
        ##
        ##     v(t_k) = w1 + Phi(t_{k-1}, t_k)^T w2,   Phi ~ I - h J + O(h^2)
        ##
        ## and with `w2 = C_{k-1}^T z` (`z = -a2 t_k`, exact by the
        ## recursion) that is `w1 + (C_{k-1} + h G)^T z` -- no inverse of
        ## `C`, so it holds for a DAE.  Equivalently `w1` is orthogonal to
        ## the amplitude eigenvector's first block, which is the true
        ## amplitude direction ROTATED by `O(h)`; `v . xdot = 1` then
        ## amplifies that rotation by `|v||xdot|`, the near-cancellation a
        ## non-isochronous oscillator has (its PPV grows with `Q_lambda`).
        ## MEASURED against the exact continuous adjoint (DOP853 at 1e-12,
        ## no shooting code in the reference) on `vdp + 0.3 u^2`, whose
        ## `c` is 100x van der Pol's: the first block gave `c` 16.6 / 8.0 /
        ## 3.9 / 1.9 / 1.0% high at 400..6400 points -- clean first order
        ## -- and violated `v(t) . xdot(t) = 1` along the orbit by 12%
        ## (std 2.7e-2).  This contraction holds the invariant to 8e-5 and
        ## gives `c` to 1.8e-3 at 400 and 8e-5 at 1600, second order.  On
        ## van der Pol both agree to 1e-4: the two rows are in quadrature
        ## there, so the rotation averaged out of `<v^2>` -- the fixture
        ## shared the claim's assumption (failure shape 0b), and `pnoise`,
        ## which contracts in PAIR space, was right all along and 14%
        ## below `c` on the fixture that could see it.
        ## The seed's scale is `w1(0) . xdot = 1`; the consistent object
        ## is renormalised ONCE by its own `v(0) . xdot`, which is why the
        ## per-step invariant is the test and not the definition.
        if (fp.is_pair and len(states) > 0
                and len(states[0]) == 2 * m):
            _cs1, _ring = [], list(fp.opening)
            for _lu, _Cn, _al, _b in fp.steps:
                _cs1.append(_ring[1])
                _ring = [_Cn, _ring[0]]
            _hs = np.diff(np.asarray(fp.times, dtype=float))
            _vphys = []
            for _j, st in enumerate(states):
                _lu, _Cn, _al, _b = fp.steps[_j]
                ## (a one-step companion -- gear's Euler backstop past the
                ## zero-stability bound on an event grid -- has no third alpha)
                _z = -(float(_al[2]) if len(_al) > 2 else 0.0) * np.asarray(_ts[_j], dtype=_dt)
                ## ⚠ DIFFERENTIAL ROWS ONLY.  `w2 = C^T z` does not see the
                ## algebraic rows of `z` (their rows of `C` are zero), so
                ## the decomposition is non-unique there, and those
                ## multipliers are O(1/h): `h G^T z` would carry an O(1)
                ## component along the constraint normal into the
                ## differential entries.  The consistent propagation
                ## `C_D dx_{k-1} = (C_D + h G_D) dx_k` involves only the
                ## differential equations, which is the choice that makes
                ## it unique.  Measured: with the algebraic rows in, a DC
                ## injection probe on a series-loss tank flipped sign.
                if _alg_rows:
                    _z[np.asarray(_alg_rows, dtype=int)] = 0.0
                _xj = _Xf[:, _j if _j < _Xf.shape[1] else -1]
                _Gj, = remove_row_col(
                    (np.asarray(self.cir.G(_xj), dtype=float),),
                    self.irefnode, self.toolkit)
                _Gj = np.asarray(_Gj, dtype=float)
                _Cj = np.asarray(_cs1[_j], dtype=float)
                if _alg_rows:
                    ## ⚠ ON A DAE THE ALGEBRAIC STATE IS SLAVED, AND ITS
                    ## COUPLING INTO THE DIFFERENTIAL PROPAGATION IS O(h).
                    ## A consistent perturbation propagates as
                    ## `C_D dx_{k-1} = (C_D + h G_red) dx_k` on the
                    ## differential states, with `G_red` the Schur
                    ## complement `G[D,NZ] - G[D,Z] G[A,Z]^-1 G[A,NZ]`.
                    ## With the full `G` instead, the series-loss tank had
                    ## `c` 0.6 / 0.3 / 0.15% low at 240/480/960 points
                    ## (first order) and the invariant drifting at 1.1e-3;
                    ## with the complement `c` is 2.7e-4 / 7e-5 / 2e-5 from
                    ## the exact reduced-ODE value and the drift 4e-4 /
                    ## 1.1e-4 / 2.7e-5 -- second order (found through the
                    ## review session's linear-DAE partition, 2026-09-05).
                    _A = np.asarray(_alg_rows, dtype=int)
                    _Zc = np.asarray(_alg_cols, dtype=int)
                    _D = np.array([i for i in range(m) if i not in _alg_rows],
                                  dtype=int)
                    _NZ = np.array([j for j in range(m) if j not in _alg_cols],
                                   dtype=int)
                    ## ⚠ AND `G[A,Z]` NONSINGULAR *IS* THE INDEX-1 CONDITION.
                    ## At index >= 2 (an L-I cutset, a C-V loop) it is
                    ## singular by definition and the algebraic variables
                    ## come from a differentiation, not a solve; the
                    ## complement does not exist.  Same shape as the fill
                    ## above: warn once with the reason and fall back to
                    ## the full-`G` propagation, which is then first order.
                    ## (Boundary named by the review session from the
                    ## pencil: `eig(-G_red, C[D,NZ])` equals the finite
                    ## generalised eigenvalues of `(C, G)` to 1e-12 on the
                    ## series-loss tank, and the reduction is undefined on
                    ## `li_plus_rc` and `cv_plus_rc`.)
                    try:
                        if len(_alg_rows) != len(_alg_cols):
                            raise np.linalg.LinAlgError('not square')
                        _Gred = (_Gj[np.ix_(_D, _NZ)]
                                 - _Gj[np.ix_(_D, _Zc)] @ np.linalg.solve(
                                     _Gj[np.ix_(_A, _Zc)], _Gj[np.ix_(_A, _NZ)]))
                    except np.linalg.LinAlgError:
                        self._ppv_alg_fallback = True
                        if _j == 0:
                            warnings.warn(
                                'PSS.ppv: the algebraic block G[A,Z] is '
                                'singular (index > 1: an L-I cutset or a '
                                'C-V loop), so the pair-consistent '
                                'propagation cannot eliminate the algebraic '
                                'state and falls back to the full G -- the '
                                'PPV samples can then be FIRST order in the '
                                'step, as they are for the algebraic fill. '
                                'Priced on the fixture that can see the '
                                'dropped term (a non-isochronous core with '
                                'its inductor split, and the same core with '
                                'a C-V loop through a bias rail): within 3e-4 '
                                'of the index-1 object and second order on '
                                'both, so at index 2 this fallback is the '
                                'whole answer.',
                                RuntimeWarning, stacklevel=2)
                        _Gred = _Gj[np.ix_(_D, _NZ)]
                    _corr = np.zeros(m, dtype=_dt)
                    _corr[_NZ] = (_Cj[np.ix_(_D, _NZ)]
                                  + _hs[_j] * _Gred).T @ _z[_D]
                    _vp_j = st[:m] + _corr
                else:
                    _vp_j = st[:m] + (_Cj + _hs[_j] * _Gj).T @ _z
                ## ⚠ AND ZERO ON THE ALGEBRAIC COLUMNS, as `C^T v_1` is:
                ## the state functional contracts a perturbation ON the
                ## constraint manifold, whose algebraic components are
                ## slaved, and `h G^T z` would otherwise leave 4e-3 there
                ## (caught by the full suite's Demir-(24) gate).  The
                ## equation-row conversion never reads these entries.
                if _alg_cols:
                    _vp_j[np.asarray(_alg_cols, dtype=int)] = 0.0
                _vphys.append(_vp_j)
            _scale = (_vphys[0] @ xdot) if _dt is complex else float(_vphys[0] @ xdot)
            if _scale == 0.0:
                raise ValueError(
                    'PSS.ppv: the pair-consistent adjoint is orthogonal to '
                    'the orbit tangent at t = 0.')
            ## ⚠ AND ITS DC CONTENT IS THE CONSISTENT OBJECT'S TOO -- taking
            ## the mean from the raw block was TRIED AND MEASURED WRONG.
            ## The raw block's orbit integral reproduces a same-grid
            ## DC-injection probe to 1e-5 on the divider fixture (node row,
            ## true mean 4e-6 |v|), where the consistent object's O(h^2)
            ## pointwise errors leave an absolute floor of ~1e-5 |v| --
            ## the wrong sign at 480 points.  But on the bias-sensitive
            ## fixture's INDUCTOR row (a DC voltage in series with L, true
            ## dT/dV = 16.20 by a second-order re-solve) the raw block
            ## reads 17.49 / 16.83 / 16.51 at 400/800/1600 -- first order,
            ## 8% off -- while the consistent object holds `v . xdot = 1`
            ## to 8e-5 along the orbit, which pins its mean in EVERY row to
            ## ~1e-5 |v|.  Stitching the raw mean in broke that invariant
            ## by +-0.3.  So the raw block's DC exactness is row- or
            ## fixture-specific (mechanism open, recorded in the roadmap),
            ## and `samples` is one object, second order everywhere, with a
            ## ~1e-5 |v| absolute floor on its mean.  The raw pair is kept
            ## as `samples_pair` for the structural gates that live on its
            ## discrete identities.
            states = [np.concatenate((vp / _scale, st[m:] / _scale))
                      for vp, st in zip(_vphys, states)]
        return states, states_pair, _ts, _Xf

    def ppv(self, tol=None):
        """The perturbation projection vector at `t = 0` (Demir & Roychowdhury).

        Returns `(v, info)`.  `v` is the pair-space left null vector of
        `I - M`, normalised so that `v . xdot(0) = 1` -- see the note on the
        normalisation below, which was MEASURED rather than transcribed.
        The phase shift caused by a state perturbation `delta` at `t = 0` is
        then `v[:m] . delta`.  `info` carries both border residuals, the
        null residual, `q` and the scaled tangent.

        ⚠ AN AUGMENTED SOLVE, NOT AN EIGENVECTOR -- AND IT IS THE FIX FOR A
        NAMED FAILURE OUR OWN FIXTURES SIT INSIDE.  Demir &
        Sangiovanni-Vincentelli, 1998 (the book, read firsthand by the docs
        session 2026-09-04), report the eigenvector route BREAKING on a
        high-Q circuit, with a table of the crowded eigenvalues (their
        Table 6.4): "Since this circuit is a high-Q one, Phi(T,0) has
        eigenvalues with magnitudes close to 1 other than the one which is
        supposed to be equal to 1 ... Because of numerical errors, we can
        not identify the eigenvalue that is supposed to be equal to 1 ...
        so it is not feasible to identify the correct" one.  Not
        ill-conditioned there -- INFEASIBLE.  That reported failure is the
        stated motivation for the single-solve method two years later, so
        the lineage is firsthand end to end: 1998 selection fails at high
        Q; 2000 the single linear solve; 2001 "particularly useful for
        high-Q oscillators"; 2003 the fuller procedure.  ⚠ It also changes
        what the second-multiplier warning below MEANS: not "this result is
        degrading" but "you are in the regime this method was invented to
        escape".  And the same book's eq (6.72), `|exp(eta_i)| << 1`, is
        the closed-form variance's validity condition -- the book says it
        "is satisfied for 'most' oscillator circuits" and defers the rest
        to the high-Q discussion above, so that condition and the crowding
        are ONE condition seen from two sides.

        ⚠ ATTRIBUTION CORRECTED
        2026-09-04: the idea ORIGINATES in Demir, Long & Roychowdhury,
        ICCAD 2000 ("Computing Phase Noise Eigenfunctions Directly from
        Steady-State Jacobian Matrices" -- "a single linear solution of the
        oscillator's ... steady-state Jacobian matrix ... dispenses with the
        need to select the correct one eigenfunction"), with the 2001
        companion carrying it to HB/shooting matrices and noting the
        selection heuristic is worst "for high-Q oscillators".  The 2003
        paper is the fuller procedure and the source of the quote below,
        not the origin.  Demir's IJCTA 2000 method SELECTED the
        right eigenvector by its inner product against `C(0) xdot(0)` --
        measured 0.2 against 1e-5, 1e-7, 2e-5 on a Colpitts.  His 2003
        paper rejects that: "no guarantee that any of the candidate
        eigenvectors will be appreciably more orthonormal than the others,
        leading to a potential breakdown."  The same vector then changes
        role -- it becomes the BORDER, so the candidate is unique and
        nothing is selected:

            [ I - M^T   q ] [v]     [0]
            [   q^T     0 ] [y]  =  [1]

        This matters here specifically.  `_spectral_report`'s eigenvector
        split was measured labelling a parasitic root physical at ~2 points
        per cycle, and multipliers cluster near 1 on exactly the high-Q
        oscillators a PPV is wanted for -- four independent sources say so.
        A bordered solve does not care how close the other multipliers are;
        it never has to tell them apart.

        ⚠ AND THAT IS THE METHOD'S STATED DESIGN DRIVER, not a lucky
        property of it.  Demir, Long & Roychowdhury (ICCAD 2000), who
        introduced the single-solve route: it is "especially advantageous
        for HIGH-Q OSCILLATORS, MONODROMY MATRICES OF WHICH OFTEN HAVE MANY
        EIGENVALUES CLOSE TO 1 THAT ARE NUMERICALLY INDISTINGUISHABLE from
        the oscillatory [unit eigenvalue]", and "a key advantage is that it
        DISPENSES WITH THE NEED TO SELECT THE CORRECT ONE-EIGENFUNCTION
        from amongst a potentially large set of choices".  So the hardest
        case in this codebase is the case the method was aimed at.

        ⚠ WHICH DOES NOT RETRACT THE SLOW-NODE BOUNDARY BELOW, and keeping
        the two apart is the point.  Near-degenerate multipliers make
        EIGENANALYSIS ILL-POSED -- there is no fact of the matter about
        which eigenfunction is the PPV -- while they make this bordered
        solve merely ILL-CONDITIONED, which is measured above and warned
        about.  Ill-conditioned beats ill-posed, and neither is the same as
        the PHASE EQUATION's own limit, which is about the response being
        treated as instantaneous and is not an extraction question at all.
        Three separate things that a "high-Q oscillators are hard" summary
        would blur into one.

        ⚠ THE QUADRATIC RUNG BELOW WOULD NOT FIX THE SLOW-NODE BOUNDARY,
        which is the obvious hope and is wrong.  TWO INDEPENDENT
        approximations are in play: the LINEAR ISOCHRON one is in the
        perturbation's AMPLITUDE -- it treats isochrons as flat
        hyperplanes, and the quadratic rung adds their curvature -- while
        the INSTANTANEOUS-RESPONSE one is in the DYNAMICS, ignoring the
        bandwidth between injection point and core.  Slow nodes are the
        second.  Noise is small by construction, so the linear term
        dominates it by definition; the quadratic rung would earn its cost
        on LARGE perturbations -- injection locking, big supply or
        substrate interferers -- not on phase noise.

        ⚠ AND THE PPV IS ONE RUNG ON A LADDER, worth knowing before it is
        mistaken for exact.  Suvak & Demir (TCAD 2011) place it: an EXACT
        phase equation exists and is "practically unusable"; the PPV
        equation is its LINEAR isochron approximation; a QUADRATIC one is
        more accurate.  Isochrons are the geometric form of asymptotic
        phase, so an oscillator without asymptotic phase is one whose
        isochrons do not exist -- the same fact as the Floquet condition,
        seen in the geometry.  Computing exact isochrons is exponential, so
        the only live question is which local approximation is affordable.

        `c` -- the diffusion constant this vector feeds -- has the
        designer-facing reading "JITTER PER SECOND".

        ⚠ THREE NAMES FOR THIS OBJECT, AND ONE NEAR-MISS THAT IS NOT IT.
        The PPV, Kaertner's adjoint LPTV impulse response, the PRC of
        mathematical biology, and Hajimiri's NUMERICAL ISF are the same
        thing.  His CLOSED-FORM ISF is NOT: it is the normalised tangent,
        and the difference is not a scale factor but a SIGN -- for a noise
        impulse at one point in the cycle "the closed-form ISF predicts a
        POSITIVE phase change, whereas in fact the correct phase change is
        in the NEGATIVE direction and of a different magnitude".  It also
        does not scale with the perturbation, where the PPV does.

        ⚠ AND `xdot` IS NOT A CHEAP SUBSTITUTE FOR IT -- "time-shifts and
        amplitudes are both different ... the two waveforms scale in
        OPPOSITE DIRECTIONS with respect to the RC time constant".  Nothing
        here offers it as one: `xdot` appears only as the NORMALISATION
        (`v . xdot = 1`, the PPV's defining property), as the border `q =
        C(0) xdot(0)`, and in the record above of the SELECTION heuristic
        that was rejected.  Stated because the substitution is a documented
        point of common confusion, and the failure would be a sign error
        rather than a visible one.

        ⚠ `y` COMING BACK ZERO IS A FREE CORRECTNESS CHECK, and it is not
        decoration.  With a zero first block on the right-hand side,
        `(I - M^T) v + y q = 0` forces `y q = 0`, so a nonzero `y` means the
        border absorbed a residual the null space should have taken -- the
        computed `v` is not in the null space.  Measured on van der Pol:
        1.4e-11.

        ⚠ ITS VALIDITY BOUNDARY IS SLOW NODES, AND THE VAN DER POL GATE
        CANNOT SEE IT.  The phase equation `alpha' = v_1^T(t+alpha) b(t)`
        treats the oscillator's frequency response as INSTANTANEOUS; the
        truth is a convolution, and the PPV form is what you get by
        assuming the kernel is `v_1(t) delta(t - tau)`.  Real circuits have
        finite bandwidth, so a slow node FILTERS the noise of devices near
        it, the PPV cannot see the filtering, and phase noise is
        OVER-ESTIMATED.  Lai (Cadence) is explicit that better extraction
        does not help: "although the PPV can be extracted correctly, the
        oscillator noise analysis is still inaccurate: the phase noise is
        always over-estimated."

        ⚠ AND HE NAMES THIS TEST'S OWN REGIME AS THE BLIND SPOT: "the
        phase equation was verified to be correct in many previous works
        ... because it was evaluated on SMALL, SIMPLE OSCILLATORS, and
        perturbations were applied to OSCILLATOR CORES.  Since oscillator
        cores have very wide bandwidth, ignoring the dynamics may not
        compromise the macromodelling accuracy very much."  That is
        `test_the_ppv_predicts_a_phase_shift_the_oscillator_actually_has`
        exactly -- van der Pol, perturbed at its core.  It passes whether
        or not this failure mode is present, so it establishes that the
        extraction and normalisation are right and says NOTHING about the
        model's range.  The fix is a frequency-aware PPV, which is this
        same bordered system at nonzero `w_s` -- the classical PPV is its
        DC point, and `PAC` already solves at nonzero frequency.  Not
        built.  VERIFIED at the source (docs session, 2026-09-08): Lai
        2008 eq. (24) at `w_s = 0` "is the augmented PPV extraction
        equation (6) and (7)", verbatim; two scope limits: (24) is a
        NEAR-DC approximation of (23) (the AC Toeplitz columns dropped,
        "if we are only interested in ... w_s close to DC"), and (23)
        "is very difficult to solve using iterative solvers (such as
        GMRES)" because the border degrades the block preconditioner.
        ⚠ AND `_vdp_with_slow_node` CANNOT SHOW THE EFFECT AT ANY tau/T:
        its slow node couples through `Rs = 1e6` against a tank impedance
        of 1, so its PPV entry is 7.3e-6 of the core's, flat over
        tau/T = 1e2..1e6 while lambda_2 moves four decades -- it tests
        conditioning, not the filtering of noise that REACHES the phase.
        The fixture owed is a slow node IN the phase path (small Rs, large
        Cs), with tau/T and coupling as separate knobs.  BUILT AND MEASURED
        the same day: with an asymmetric core AND tank loss (G_0 != 0 needs
        both) the slow node's PPV entry is DC-dominated (|G0|/|G1| = 30) and
        the Lorentzian over-states a source behind it by 1000x at 0.1 f0,
        predicted to four digits from this PPV's harmonics and the RC
        filter -- so `c` from this PPV is right and the SHAPE above
        T/(2 pi tau) is what the frequency-aware PPV corrects; see
        `oscillator_spectrum`.  A2's "gated at tau/T = 10" was a Monte
        Carlo of `c`, which cannot see it.

        ⚠ `q` IS EXACT, NOT DIFFERENCED.  `q = C(0) xdot(0)` looks like it
        needs the orbit's tangent, and differencing the waveform for it
        would be O(h) at best.  The DAE gives it directly: `dq/dt + i(x) +
        u(t) = 0` and `dq/dt = C xdot`, so `q = -(i(x_0) + u(0))` -- two
        evaluations at the converged solution, no derivative anywhere.
        
        ⚠ ON A STAGED SOLVE (`state_events=True`) THE PPV IS BORDERED
        (2026-09-22, events phase B): the null vector is that of the TOTAL
        monodromy `M + P_theta dtheta/dx_0`, and the samples carry the
        crossings' motion as costate injections `-zeta_k W_k` at the event
        nodes, `zeta = Gt^-T P_theta^T v` -- one reverse pass, and it IS
        the phase gradient at fixed time (the tail-restricted Newton plus
        the node's own motion agrees with it to 1e-4).  VERIFIED against
        the exact piecewise-linear PPV of a comparator relaxation
        oscillator (linear flows joined by saltation matrices): every
        sample within 0.4 % before, between and after the crossings,
        where the fixed-grid PPV of the same solve is 140 % off before
        the first crossing and its `M^T v - v` residual is 1.5.  ⚠ THE
        TRANSIENT FD THAT WAS TO BE THE INSTRUMENT read `c` at +4.4e-8
        for the exact -7.4e-9 s/V, converged in its own step, period to
        4e-6, along the flow to 1e-4: its phase was read at the `c`
        waveform's crossing of ``(max + min) / 2`` OF THE PERTURBED RECORD
        -- a level the perturbation itself moved.  Read at the
        comparator's own crossing it agrees with the exact value to 0.4 %.
        An instrument's reference level must not come from the record it
        measures.
"""
        _tw = self.monodromy_twin()
        if _tw is not self:
            return _tw.ppv(tol)
        import scipy.sparse.linalg as spla
        fp = self.factored_period()
        ## ⚠ NO LONGER GEAR-ONLY (B8). The refusal that stood here said the
        ## transposed replay was "implemented for the solved-history map
        ## only"; since `_monodromy_matvec_transposed_plain` shipped that
        ## sentence is false, and every call below goes through
        ## `fp.matvec_transposed`/`fp.matvec` so the map's kind is the
        ## dispatcher's business rather than this method's.
        if not self.autonomous:
            raise ValueError(
                'PSS.ppv: a perturbation projection vector describes the '
                'phase of a FREE-RUNNING oscillator. This circuit is '
                'driven, so its phase is the source\'s and there is no '
                'unit Floquet multiplier to project onto.')

        m = self.cir.n - 1
        n = fp.width
        irn = self.irefnode
        x0r = np.asarray(self._period_state[1], dtype=float).ravel()
        x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))
        qf = -(np.asarray(self.cir.i(x0f)).ravel()
               + np.asarray(self.cir.u(0.0,
                                       analysis=self.par.analysis)).ravel())
        q = np.delete(np.asarray(qf, dtype=float), irn)
        qp = np.concatenate((q, np.zeros(n - m)))
        nq = float(np.linalg.norm(qp))
        if nq == 0.0:
            raise ValueError(
                'PSS.ppv: C(0) xdot(0) is zero, so the orbit has no tangent '
                'at t=0 and there is nothing to normalise against. That '
                'should not happen on a converged limit cycle.')

        ## ⚠ ON A STAGED SOLVE THE MONODROMY IS THE TOTAL ONE (2026-09-22,
        ## events phase B): a perturbation moves the crossings, `M_tot = M
        ## + P_theta dtheta/dx_0`, and the PPV is its left null vector.
        ## Its samples along the orbit carry the same correction as costate
        ## injections at the event nodes -- `-zeta_k w_k`, `zeta = Gt^-T
        ## P_theta^T v` -- which the reverse pass carries to every earlier
        ## node: the saltation matrix's transpose, derived from the bordered
        ## system rather than guessed.  A comparator oscillator's PPV jumps
        ## at its switching instants, and this is where the jump comes from.
        _ev = EventColumns.of(self, n)
        if _ev is None:
            _MtT = lambda v_: np.asarray(fp.matvec_transposed(v_))
            _Mt = lambda u_: np.asarray(fp.matvec(u_))
        else:
            _MtT = _ev.total_matvec(fp.matvec_transposed, transposed=True)
            _Mt = _ev.total_matvec(fp.matvec)
        def _mv(z):
            z = np.asarray(z)
            v_, y_ = z[:n], z[n]
            top = v_ - _MtT(v_) + y_ * qp
            return np.concatenate((top, [float(qp @ v_)]))

        rtol = max(self.par.reltol * 1e-2 if tol is None else tol, 1e-14)
        A = spla.LinearOperator((n + 1, n + 1), matvec=_mv, dtype=float)
        rhs = np.zeros(n + 1)
        rhs[n] = 1.0
        ## ⚠ JUDGED BY ITS RESIDUAL, NOT BY A STATUS CODE.  SciPy returns
        ## `info = 4` on a HAPPY breakdown -- the Krylov space exhausted
        ## because the answer is exact -- and these bordered operators are
        ## small enough to hit that routinely.  `_arnoldi_gmres` detects
        ## the breakdown where it happens and returns the exact answer.
        z, relres, _Ha, _ka = _arnoldi_gmres(
            _mv, rhs, rtol=rtol, maxiter=min(n + 1, 200))
        if relres > max(1e3 * rtol, 1e-8):
            raise RuntimeError(
                'PSS.ppv: the augmented solve did not converge (relative '
                'residual %.3e). The border is `q` itself; if the orbit '
                'tangent is nearly orthogonal to the null direction the '
                'bordering is poor.' % relres)
        v, y = z[:n], float(z[n])
        resid = float(np.linalg.norm(v - _MtT(v)))

        ## ⚠ THE SCALE NEEDS THE TANGENT, so the RIGHT null vector is solved
        ## for too -- by the same bordering, not by an eigendecomposition,
        ## for the same reason: on a high-Q oscillator the other multipliers
        ## crowd 1 and no selection among candidates is reliable.
        def _mvf(zz):
            zz = np.asarray(zz)
            u_, yy = zz[:n], zz[n]
            top = u_ - _Mt(u_) + yy * qp
            return np.concatenate((top, [float(qp @ u_)]))

        Af = spla.LinearOperator((n + 1, n + 1), matvec=_mvf, dtype=float)
        zf, relresf, _Hf, _kf = _arnoldi_gmres(
            _mvf, rhs, rtol=rtol, maxiter=min(n + 1, 200))
        if relresf > max(1e3 * rtol, 1e-8):
            raise RuntimeError(
                'PSS.ppv: the tangent solve did not converge (relative '
                'residual %.3e).' % relresf)
        u, yf = zf[:n], float(zf[n])

        ## `u` is the tangent's DIRECTION; its scale comes from `C u = q`,
        ## which is the definition of `q` read backwards.  Least squares
        ## because `C` is singular for a DAE and only its range is
        ## determined.
        x0red = np.asarray(self.cir.C(x0f))
        Cm = np.delete(np.delete(x0red, irn, 0), irn, 1)
        Cu = np.asarray(Cm, dtype=float) @ u[:m]
        denom = float(Cu @ Cu)
        if denom == 0.0:
            raise ValueError(
                'PSS.ppv: C(0) annihilates the orbit tangent, so its scale '
                'is not determined by C u = q.')
        xdot = u[:m] * (float(q @ Cu) / denom)

        ## ⚠ NORMALISED BY `v . xdot = 1`, WHICH IS NOT WHAT `v . q = 1`
        ## GIVES, and the difference is not cosmetic: on `_vdp_ppv(400)`
        ## `v . xdot = 1.0` against `v . q = -1.0696` -- a factor 2.07 AND
        ## the opposite sign (an earlier version of this note said "7%",
        ## which was the |v . q| - 1 residual and not the error; corrected
        ## by the review session's audit, 2026-09-04).  The defining
        ## property is that displacing the state
        ## ALONG the orbit by `eps xdot` advances the phase by `eps`, so
        ## `v . xdot = 1` is the normalisation a state perturbation sees.
        ## Demir's Remark 3.1 reads `v_1^T C u_1 = 1`; the vector this
        ## bordered solve returns behaves as `C^T v_1` -- it is contracted
        ## with a state perturbation directly -- so the two statements agree
        ## about different objects.
        ##
        ## ⚠ AND THAT IS NOT A QUIRK OF THIS FORMULATION, WHICH THIS
        ## COMMENT USED TO IMPLY.  The conserved pairing propagates to
        ## `M^T (C(0)^T v_1) = C(0)^T v_1`, so the left eigenvector of the
        ## STATE-SPACE monodromy simply IS `C(0)^T v_1` -- for ANY `C`,
        ## symmetric or not, and whatever the augmentation.  MEASURED on a
        ## limit cycle with a constant NON-SYMMETRIC `C`: alignment with
        ## `C(0)^T v_1` is 1.000000000000 against 0.9657 for `v_1` itself,
        ## and bordering with `xdot(0)` reproduces Demir's normalisation
        ## exactly while bordering with `C(0) xdot(0)` gives 0.805.
        ## (Derived and measured by the docs session, 2026-09-04.)  ⚠ TREATING THEM AS THE SAME OBJECT WAS
        ## MEASURED WRONG: predicting a state jump's phase shift as
        ## `v^T C delta` gives residuals of 0.36/0.40/0.42 that GROW with
        ## refinement and per-direction ratios scattering from -0.44 to
        ## 28.7, while `v . delta` converges at O(h).
        ## ⚠ THE ALGEBRAIC ENTRIES ARE FILLED *AFTER* THIS, AND THAT IS A
        ## DECISION RATHER THAN AN ORDERING ACCIDENT.  Filling first was
        ## tried and MEASURED WORSE: the DC-injection probe went from
        ## 0.9999849 to 0.9992364, because `v` here is also the REPLAY'S
        ## SEED and the algebraic components are SLAVED -- propagating them
        ## through the step map corrupts the differential ones.
        ##
        ## ⚠ AND THE NORMALISATION SHOULD NOT SEE THEM EITHER.  `v . xdot`
        ## is about a STATE perturbation, and a state perturbation of a DAE
        ## lies ON the constraint manifold: its algebraic components are
        ## determined by its differential ones, not free.  The algebraic
        ## entries of `v` answer a different question -- the sensitivity to
        ## a perturbation of an EQUATION ROW, which is what a noise current
        ## injected into an algebraic KCL row is.  So this line is
        ## unchanged, and every PPV number on every circuit is
        ## bit-for-bit what it was.
        ## ⚠ THE UNIT MULTIPLIER OFF THE CIRCLE IS A SILENT ERROR IN `c`
        ## (2026-09-21).  This solves the bordered system AT lambda = 1
        ## whatever the discrete period map's own unit multiplier is, so
        ## nothing here can see that multiplier sit at 1.10 -- and on a
        ## relaxation van der Pol (mu = 10) on its own 195-point `lte_grid`,
        ## gear's did: the PPV then gave a diffusion constant 52 % high (14.7 %
        ## at 2x, 4.0 % under trbdf2 at |rho - 1| = 3.4e-3), tracking the
        ## departure at 5-12x, with `converged = True` and no message.  The
        ## solve already records `spectral_radius`; an autonomous run whose
        ## unit multiplier is more than `PPV_UNIT_MULTIPLIER_WARN` off the
        ## circle is told so here, once, with the size and the remedy.
        _rho = getattr(self, 'spectral_radius', None)
        if (getattr(self, 'autonomous', False) and _rho is not None
                and np.isfinite(_rho)
                and abs(float(_rho) - 1.0) > self.PPV_UNIT_MULTIPLIER_WARN):
            warnings.warn(
                'PSS.ppv: the discrete period map\'s unit multiplier sits '
                '%.2e off the unit circle (spectral_radius %.6f). The PPV is '
                'solved at exactly 1, so it cannot see this, and the '
                'diffusion constant built from it is uncertain by about '
                '5-12x that departure (measured on a relaxation oscillator: '
                '1.0e-1 -> c 52 %% high). Refine the grid, or use a method '
                'whose multiplier stays on the circle at this step count '
                '(radau, or trbdf2 on the same grid).'
                % (abs(float(_rho) - 1.0), float(_rho)), RuntimeWarning,
                stacklevel=2)
        _alg_rows, _alg_cols = self._algebraic_adjoint_pattern(x0f)
        vx = float(v[:m] @ xdot)
        if vx == 0.0:
            raise ValueError(
                'PSS.ppv: the null vector is orthogonal to the orbit '
                'tangent, so no normalisation makes it a phase projector.')
        v = v / vx
        ## ⚠ THE PPV OVER THE PERIOD, not just at `t = 0`, because that is
        ## what an oscillator noise calculation needs: Demir's diffusion
        ## constant is `c = (1/T) integral v_1^T(t) B(t) B^T(t) v_1(t) dt`,
        ## an integral over the orbit.  `Phi(T,s)^T v(T) = v(s)`, and the
        ## reverse replay computes exactly that sequence on its way to the
        ## answer -- it was being discarded.
        _inject = self._event_costate_injection(fp, v, n)
        states, states_pair, _ts, _Xf = self._ppv_propagate(fp, v, m, xdot, _alg_rows, _alg_cols,
                                                          inject=_inject)
        ## ⚠ INDEX >= 2 ON A NON-UNIFORM SOLVED-HISTORY GRID (2026-09-20, item
        ## 3 of "gear as a first-class choice on non-uniform grids"): the
        ## fallback above keeps only the differential block of the
        ## pair-consistent correction, and the coupling it drops -- a
        ## DERIVATIVE term at index 2, with no Schur complement to stand in
        ## for it -- cancels between steps on a uniform grid and not on a
        ## non-uniform one.  Measured on a van der Pol with a DC source inside
        ## a capacitor loop, smooth grid, `c` against radau: -5.0e-3 / -2.4e-3
        ## / -1.2e-3 / -5.8e-4 at N = 100..800 (FIRST order) where the uniform
        ## grid gives +1.2e-3 / +3.3e-4 / +8.7e-5 / +2.2e-5 and trap on the
        ## same smooth grid is second order.  So on that grid the samples come
        ## from the continuous adjoint's phase mode instead -- the object
        ## `floquet_modes` already uses there, whose invariant quarters on
        ## this fixture -- as the STATE-SPACE PPV `v_j = C_j^T q_j` (the
        ## continuous `q` IS Demir's equation-row `v_1`, and `samples` carry
        ## `C^T v_1` -- see `_equation_row_ppv`; `q` put there bare lands
        ## 16.7x off with the right shape), scaled so `q_0^T C_0 xdot_0 = 1`,
        ## which is `ppv()`'s own `v . xdot = 1`.  Measured with a periodic
        ## trapezoid over the actual grid: +4.0e-5 / +3.6e-5 / +1.3e-5 at
        ## N = 200 / 400 / 800, the reference's floor.  ⚠ The second half of
        ## that measurement was PAC's period weights (`_period_weights`):
        ## a left-rectangle rule on a non-uniform grid is first order by
        ## itself.  The anchor `v` and the pair block are untouched; the
        ## uniform grid never comes here, and the one-step kinds cannot:
        ## their `factored_period_full` / `_dirk` replays (now
        ## `factored_period_stage`) were built on a uniform `linspace` grid
        ## whatever the solve used, which is WHY
        ## radau read as "exact on the 3:1 grid" -- its adjoint surfaces
        ## never saw that grid.
        if (getattr(self, '_ppv_alg_fallback', False)
                and getattr(fp, 'is_pair', False)
                and self._period_quadrature(fp) is not None):
            _tms = np.asarray(self.waveform[0], dtype=float)
            _qc = np.real(self._continuous_adjoint(fp, 1.0 + 0.0j, 0.0, _tms)[0])
            _Wr = np.delete(np.asarray(self.waveform[1], dtype=float),
                            self.irefnode, axis=0)
            _C0 = np.asarray(self._C_at(_Wr[:, 0]), dtype=float)
            _sN = float(_qc[:, 0] @ (_C0 @ np.asarray(xdot, dtype=float)))
            if _sN != 0.0 and _qc.shape[1] >= len(states):
                _new = []
                for _sj, st in enumerate(states):
                    _Cj = np.asarray(self._C_at(_Wr[:, _sj]), dtype=float)
                    _vj = (_Cj.T @ _qc[:, _sj]) / _sN
                    _new.append(np.concatenate((_vj, np.asarray(st)[m:])))
                states = _new
        ## ⚠ AND FILL EVERY SAMPLE TOO, at ITS OWN operating point, because
        ## `G` is state-dependent and the algebraic entries are a pointwise
        ## function of the differential ones.  Done here rather than by
        ## seeding the replay: these components are SLAVED, so there is
        ## nothing to propagate, and post-processing leaves the validated
        ## step map untouched.  The pair's SECOND block is the history term
        ## and is deliberately not filled -- `v(t)` is the first block.
        ## ⚠⚠ THE EQUATION-ROW ADJOINT IS A SECOND OBJECT, NOT A CORRECTION
        ## TO THE FIRST.  `states` and `v` stay exactly what they were --
        ## `C^T v_1`, the vector a STATE perturbation contracts with, which
        ## is what this method documents and what every existing gate
        ## measures.  Demir gives both conventions on one page (eq 41 with
        ## `C`, eq 42 and the phase equation 44 bare), so naming both is
        ## the fix; converting one into the other would have silently
        ## changed what `ppv()` returns.
        ## ⚠ ONE WARNING PER CALL, NOT ONE PER SAMPLE: the fill warns when
        ## `G[A,Z]` is singular, and at index 2 it is singular at every
        ## sample -- 240 identical warnings for one call, which trains a
        ## reader to filter this module's warnings and miss a real one.
        with warnings.catch_warnings(record=True) as _caught:
            warnings.simplefilter('always')
            _eq = [self._equation_row_ppv(
                       st[:m], _Xf[:, _sj if _sj < _Xf.shape[1] else -1],
                       _alg_rows, _alg_cols)
                   for _sj, st in enumerate(states)]
        _seen = set()
        for _w in _caught:
            _key = (str(_w.message), _w.category)
            if _key not in _seen:
                _seen.add(_key)
                warnings.warn(str(_w.message), _w.category, stacklevel=2)
        _v_eq = self._equation_row_ppv(v[:m], x0f, _alg_rows, _alg_cols)
        ## ⚠ A SECOND MULTIPLIER NEAR 1 BREAKS THIS SILENTLY, and none of
        ## the residuals above can see it.  The border removes the PHASE
        ## mode's singularity and does nothing about any OTHER root
        ## approaching the unit circle -- which a slow node puts there.
        ## MEASURED on van der Pol with one weakly coupled RC node:
        ##
        ##     tau/T    |lambda_2|    sigma_min(bordered)   null residual
        ##     none      0.000856          8.62e-01            4.1e-11
        ##     1e2       0.990049          4.49e-03            4.6e-11
        ##     1e4       0.999900          4.47e-05            4.6e-11
        ##     1e6       0.999999          4.47e-07            4.4e-11
        ##
        ## `sigma_min` tracks `T/tau` over six decades while the residual
        ## does not move at all: GMRES converges, the answer looks clean,
        ## and the conditioning has lost six digits.  So this estimates
        ## `|lambda_2|` explicitly rather than trusting a small residual.
        ##
        ## ⚠ AND THE ACCURACY COST WAS GATED, WITH A NEGATIVE RESULT worth
        ## recording so nobody re-derives a fix from the warning alone.
        ## Monte Carlo on the FULL NONLINEAR circuit -- 200 realisations,
        ## 150 periods, phase read from zero-crossing timing, so no PPV
        ## appears anywhere in the measurement:
        ##
        ##     core injection (control)   c_ppv/c_meas = 0.9965
        ##     slow node, tau/T = 10      c_ppv/c_meas = 0.8016
        ##
        ## Within 20%, about 2 sigma at this sample count, and in the
        ## UNDER-predicting direction.
        ##
        ## ⚠ BUT THAT IS NOT A FALSIFICATION, AND THIS DOCSTRING SAID IT
        ## WAS.  `tau/T = 10` is OUTSIDE the regime the reported mechanism
        ## needs: it bites through ill-conditioning, and by the table above
        ## `sigma_min` at `tau/T = 10` is ~4.5e-02 -- healthy.  The PPV has
        ## no large entries there and nothing is splitting into two nearly
        ## cancelling components.  Lai's own case is a gated-capacitor
        ## tuning bank (226 MOSFETs, 3.15 GHz) whose off-caps have RC
        ## exceeding ~1 s, i.e. `tau/T ~ 3e9` -- eight orders from what was
        ## tested.  A null result at 10 is what the mechanism PREDICTS, not
        ## evidence against it.
        ##
        ## ⚠⚠ PROVENANCE, AND THE CHAIN IS NOW FULLY TRACED -- A UNIT WAS
        ## MANUFACTURED IN TWO STEPS.  This used to render "larger than 1
        ## second" AS A QUOTATION.  The primary source IS on disk, at
        ## `~/docs/09-phase-macromodels-and-prc/Lai-2008-Frequency-Aware
        ## PPV ... (Cadence).pdf`, and p.4 reads, verbatim:
        ##
        ##     "Since the RC time constants of the "off" gated capacitors is
        ##      very large (LARGER THAN 1), it is safe to assume that these
        ##      gates have very small contribution to the total phase noise
        ##      when offset frequency is reasonably large."
        ##
        ## **NO UNIT.**  Our own reading of the paper
        ## (`~/docs/pycircuit-frequency-aware-ppv.md`) paraphrased it as
        ## "their RC constants exceed 1 s" -- ADDING the unit, and unmarked,
        ## beside that file's properly marked quotations.  This comment then
        ## promoted the paraphrase to a QUOTATION, carrying the added unit
        ## with it.  Two steps, each small, and the result was a quoted unit
        ## the source does not contain.
        ##
        ## ⚠ Seconds remains the natural reading (the `tau/T ~ 3e9` above
        ## follows from it and nothing downstream moves), but it is OURS and
        ## is marked as such.  ⚠⚠ AND THE FIRST VERSION OF THIS CORRECTION
        ## SAID THE PDF WAS "NOT ON DISK AT ALL" -- it is, in a
        ## SUBDIRECTORY, and the search that missed it looked only at the
        ## top level of `~/docs`.  Search a library recursively before
        ## reporting a source missing.
        ##
        ## ⚠ SO THE HONEST RECORD IS: not reproduced at `tau/T = 10`, which
        ## is outside the regime where the mechanism predicts an effect;
        ## UNTESTED at the `tau/T ~ 1e9` where it is reported.  And the
        ## reason the fix is still not built is COST, not falsification:
        ## the measurement needs ~15 time constants of settling, so at
        ## `tau/T = 1e4` that is 150 000 periods per realisation.  That
        ## argument stands on its own; the falsification framing does not,
        ## and this codebase's ledger distinguishes them.
        ##
        ## ⚠ AND THE 0.80 IS IN THE OPPOSITE DIRECTION TO THE REPORTED
        ## EFFECT.  If it survives the ~10% Monte Carlo uncertainty at 200
        ## realisations it is a separate ~20% UNDER-prediction at a `tau/T`
        ## where the conditioning is fine -- not a weak version of Lai's.
        ## At ~2.5 sigma it is not established either way, and it is
        ## recorded rather than resolved.
        ##
        ## ⚠ Larger `tau/T` is untested and the cost is why: the
        ## measurement needs ~15 time constants of settling.
        ## ⚠ AND IT TOOK THREE ATTEMPTS.  A window of 2-4 time constants
        ## read the slow mode's DECAY as diffusion; an impulse test could
        ## not resolve a 1e-11 time shift; and one noise amplitude for both
        ## circuits put a 2.5 V jump per step on an orbit of amplitude 2,
        ## because the slow node's capacitance is 6.7e-5 F against the
        ## core's 1.0.  Each time the number was read before the
        ## MEASUREMENT was shown to be in the regime it assumes.
        ##
        ## ⚠ ARNOLDI RITZ VALUES, NOT A DEFLATED POWER ITERATION -- and the
        ## replacement is BOTH more accurate and cheaper, which is rare
        ## enough to state plainly.  Power iteration converges at
        ## `|lambda_3|/|lambda_2|`, so it fails exactly where a parasitic
        ## multiplier crowds the oscillatory one.  Arnoldi does not care
        ## about that ratio.  MEASURED on van der Pol at `Q = 16` with one
        ## parasitic RC swept through it:
        ##
        ##     tau_p/T   lam2/lam3   POWER err   RITZ err
        ##       1        2.554      3.53e-14    0
        ##       4        1.206      1.52e-06    1.11e-15
        ##       8        1.065      1.41e-03    0
        ##      16        1.000      1.10e-05    2.22e-16
        ##      32        1.032      4.61e-03    2.22e-16
        ##     100        1.054      2.48e-03    2.22e-16
        ##
        ## Machine precision everywhere INCLUDING at exact degeneracy,
        ## against a power iteration losing three digits at a ratio of
        ## 1.065 -- and at `k` matvecs rather than
        ## `PPV_DEFLATION_ITERS = 30`.
        ##
        ## ⚠ THE ROUTE IS GARCIA, ROMERO & ACHA (IEEE Trans. Power
        ## Systems 37(1), 2022): Ritz values of `I - M` map back as
        ## `lambda = 1 - theta`.  They take `H` from the GMRES that
        ## already solved the Newton correction; here it is a small
        ## dedicated Arnoldi, because this call has no GMRES of its own.
        ##
        ## ⚠ EXACT AT `k = n` AND A TRUNCATION OTHERWISE.  The cap keeps a
        ## large circuit from paying `n` matvecs for a diagnostic.
        ##
        ## ⚠ A TRUNCATED `lam2` IS A LOWER BOUND, so the near-unit warning
        ## below can only UNDER-fire -- and that is Cauchy interlacing,
        ## not an accident: `theta_j >= lambda_j(A)` for a Rayleigh-Ritz
        ## projection, hence `1 - theta_2 <= lam2`.  MEASURED on a
        ## synthetic 40x40 with a verified-normal `M`
        ## (`||M^H M - M M^H||/||M||^2 = 8.9e-16`): a lower bound in
        ## 100.0% of 200 draws at `k` = 3, 5, 8 and 12.
        ##
        ## ⚠⚠ AND THIS SURVIVED A ROUND TRIP THROUGH A FALSE REFUTATION,
        ## which is why it is written out.  An intermediate version of
        ## this comment said "over-estimate, not a lower bound", on a
        ## measurement whose selection rule took the largest `|1 - theta|`
        ## after discarding `|lam - 1| < 1e-8`.  On a truncated basis that
        ## picks the UNCONVERGED UNIT-MODE Ritz value -- below 1 but above
        ## `lam2` -- so it measured its own filter and attributed the
        ## result to the phase mode contaminating `lam2`.  Selecting
        ## `theta_2` as the second-smallest Ritz value of `I - M`, which
        ## is what interlacing is about, restores the bound.
        ##
        ## ⚠ THE BOUND IS FOR A NORMAL `M`.  At forced eigenvector
        ## conditioning `cond(V) = 1e6` it fails in 100% of draws at
        ## `k = 20` -- but by exactly `1 - lam2`, i.e. the unit mode being
        ## SELECTED as `lam2`, a selection failure rather than a Ritz one:
        ## at that conditioning `|lam_1 - 1|` is 1.5e-8 to 4.3e-8 and no
        ## value-based rule separates it.  A fixed absolute tolerance is
        ## the weak point; deflating the phase mode explicitly with `q`,
        ## which this method already has, would sidestep it.
        ##
        ## ⚠⚠⚠ AND THE ESCAPE CLAUSE IS LOAD-BEARING: A CIRCUIT MONODROMY IS
        ## NOT NORMAL, AND THE BOUND FAILS ON ONE.  MEASURED 2026-09-07 on
        ## THIS FILE'S OWN `_osc_with_ladder(Q, 14, nslow)` against the dense
        ## spectrum of the SAME operator (`n` matvecs, the route the
        ## `dirk`/`full` branch below already takes), scored in the GAP
        ## because `Q ~ 1/(1 - lam2)`::
        ##
        ##     nslow   dense lam2     k=12 Arnoldi   gap ratio   Q dense/Arnoldi
        ##      <=11   0.995706203    0.995706197      1.000       232 / 232
        ##        12   0.996324417    1.000114048     -0.031       271 / inf
        ##        13   0.996818781    0.942674586     18.020       313 / 16.9
        ##        14   0.997220139    0.999318472      0.245       359 / 1467
        ##
        ## **THE ERROR IS NOT ONE-SIGNED**, so "can only UNDER-fire" does not
        ## hold here: 13 under-estimates (19x low in `Q`), 12 and 14
        ## OVER-estimate, and 12 returns `lam2 > 1` -- a spurious UNSTABLE
        ## multiplier, which `Q` reports as `inf`.  The two failures are
        ## different: at 13 the Arnoldi never resolves `0.99682` and selects
        ## the next TRUE eigenvalue down (`0.9427`); at 14 it selects a
        ## SPURIOUS Ritz value at `0.99932` that is no eigenvalue at all.
        ##
        ## ⚠⚠ SO "NOT LIVE ON A CIRCUIT MONODROMY" (below) WAS MEASURED ON THE
        ## WRONG AXIS.  It was checked against EIGENVECTOR CONDITIONING; the
        ## trigger here is the number of distinct near-unit CLUSTERS, which
        ## `_osc_with_ladder` varies BY CONSTRUCTION and which the fixture's
        ## own test already reports reaching ~29 at `nslow = 14`.  It is not
        ## Q-specific either: the same `nslow` fails at Q = 8/16/256.
        ##
        ## ⚠⚠ IT IS A SIZING PROBLEM, AND RAISING THIS CONSTANT IS THE WRONG
        ## FIX -- MEASURED.  `k = 16` is exact on the fixture above and FAILS
        ## on a longer ladder, because the required `k` grows with `n`::
        ##
        ##     ladder/nslow   n    k=12 ratio   k=16 ratio   k=20 ratio
        ##        14 / 14     32      0.245        1.000        1.000
        ##        20 / 20     44      0.277        0.279        1.000
        ##        26 / 26     56      2.313        0.410        0.265
        ##
        ## `k ~ n/2` and rising, against a DENSE route that costs `n` and needs
        ## no threshold at all.
        ##
        ## ⚠⚠ BUT `k ~ n/2` IS AN ARTEFACT OF THIS FIXTURE, AND THE FIXTURE
        ## CANNOT SEE IT.  `_osc_with_ladder` sets `nslow = nladder`, so `n`
        ## and the slow-mode count move together here and no measurement on
        ## it can separate "k tracks n" from "k tracks nslow".  A peer
        ## session's synthetic CAN separate them and reports `k_min` rising
        ## with `nslow` and FLAT under a doubling of `n` at fixed `nslow`
        ## (24->24, 24->16, 48->48, 48->48).  If that transfers, the rule is
        ## **cost tracks the SLOW-MODE COUNT, not the system size** -- a
        ## large fast circuit is cheap and a small one with a big tuning
        ## bank is not, which also says Lai's 813-equation oscillator is
        ## expensive because of the BANK and not the 813.  Recorded with
        ## that provenance: measured on a synthetic, consistent with
        ## everything measured here, and NOT separable on this fixture.
        ##
        ## ⚠⚠ BUT THE DENSE ROUTE IS OUT ON THE CIRCUITS THAT MOTIVATE THIS.
        ## `FLOQUET_DENSE_LIMIT = 400`, and the published cases are LARGER:
        ## Lai's 64-gated-capacitor DCO is "about 200 transistors, and the
        ## system size is more than 500 ... We have trouble to apply direct
        ## harmonic balance in this case due to memory issue" (DAC 2006
        ## p.1021, verified on disk), and [L08]'s tuning oscillator is 813.
        ## So dense is the right default only in the `n <= 400` band this
        ## class already draws, and the RITZ-RESIDUAL gate is what the large
        ## end needs.  ⚠ And a 64-element bank with any realistic fraction
        ## off is an order of magnitude past the >= 3 clusters that break
        ## `k = 12` -- i.e. the extension AT ITS CURRENT BASIS SIZE would
        ## fail on exactly the circuits it exists for.
        ## THE DIAGNOSTIC THAT SEPARATES THEM CLEANLY IS THE PER-PAIR RITZ
        ## RESIDUAL `|h_{k+1,k}| |y_i[last]|`, free from `H`: 1.0e-02 at
        ## k=8, 2.1e-03 at k=12 (both wrong), 1.5e-16 at k=16 (right), and
        ## <=3.1e-07 at every `nslow` the shipped path gets right.  Neither
        ## is built -- see the roadmap; `lam2` and `Q` are REPORTED
        ## DIAGNOSTICS with no non-test consumer, so nothing computes wrong,
        ## but a caller reading `info['Q']` on a bias network with many long
        ## time constants can be off by 4x to 19x, silently.
        ##
        ## ⚠ ONE FAILURE MODE CHECKED AND NOT LIVE HERE -- ⚠⚠ SUPERSEDED BY
        ## THE MEASUREMENT ABOVE, KEPT BECAUSE IT RECORDS WHAT WAS TESTED.  At
        ## `cond(V) >= 1e4` the `|lam - 1|` filter itself fails: the phase
        ## mode stops being resolved to the tolerance, survives the
        ## discard, and is selected as `lam2`, sending `Q` to infinity.
        ## MEASURED on what was then this class's stiffest realistic fixture
        ## -- a Q=60 oscillator with a 10-mode damped bulk, `m = 12` --
        ## `cond(V) = 92` and `|lam_1 - 1| = 3.0e-13`, seven orders inside
        ## the 1e-6 filter.  That axis is still clean; the CLUSTER-COUNT axis
        ## is not.  Sorted by real part, not magnitude, because an
        ## amplitude mode is real and positive while a complex pair of
        ## larger modulus would be an oscillation about the orbit.
        vu = float(v[:m] @ u[:m] + v[m:] @ u[m:])
        lam2 = 0.0
        ## `_resid`/`_certified` describe how `lam2` was obtained; the
        ## degenerate `n < 2` fall-through never enters either branch, and
        ## `lam2 = 0` there is exact rather than estimated.
        _resid, _certified = 0.0, True
        kk = int(min(n, self.PPV_RITZ_BASIS))
        ## ⚠⚠ DENSE WHENEVER IT IS AFFORDABLE, AND THAT IS NOW THE DEFAULT
        ## RATHER THAN A STAGE-METHOD CARVE-OUT.  Forming `M` by `n` matvecs
        ## and taking its exact spectrum has no threshold, no basis size and
        ## no selection ambiguity; the truncated Arnoldi below has all three.
        ##
        ## It used to run only for `dirk`/`full`, on the argument quoted
        ## below -- and that argument was never stage-specific.  MEASURED
        ## 2026-09-07 on `_osc_with_ladder(16, 14, nslow)` (`gear`, so the
        ## Arnoldi path), against this same dense spectrum, scored in the GAP
        ## because `Q ~ 1/(1 - lam2)`::
        ##
        ##     nslow   dense lam2     k=12 Arnoldi   gap ratio   Q dense/Arn
        ##      <=11   0.995706203    0.995706197      1.000       232 / 232
        ##        12   0.996324417    1.000114048     -0.031       271 / inf
        ##        13   0.996818781    0.942674586     18.020       313 / 16.9
        ##        14   0.997220139    0.999318472      0.245       359 / 1467
        ##
        ## Not one-signed, so the Cauchy lower bound recorded above does not
        ## hold on a circuit monodromy (it is stated for a NORMAL `M`, and
        ## this is not one); and at `nslow = 12` it reports `lam2 > 1`, a
        ## spurious UNSTABLE multiplier, which `Q` turns into `inf`.
        ##
        ## ⚠ RAISING `PPV_RITZ_BASIS` IS NOT THE FIX AND WAS MEASURED NOT TO
        ## BE: `k = 16` is exact on that fixture and fails on a longer
        ## ladder (20/20 -> 0.279, 26/26 -> 0.410), because the basis has to
        ## grow with the problem.  A constant cannot.
        ##
        ## ⚠ `FLOQUET_DENSE_LIMIT` is the same cap `floquet_modes` applies to
        ## the same assembly, so the two agree about what "affordable" means.
        ## `dirk`/`full` keep the dense route ABOVE it as well: there it is
        ## expensive, but the alternative is not slower, it is WRONG, and
        ## those paths have never had the truncated one.
        _dense_ok = (fp.is_stage
                     or n <= self.FLOQUET_DENSE_LIMIT)
        if _dense_ok:
            ## ⚠ THE STAGE MAP IS DENSE AND WIDTH `m`, so its exact spectrum
            ## is cheap -- and the Arnoldi below resolves it BADLY here.
            ## `I - M` has `M`'s annihilated modes clustered at eigenvalue 1
            ## and the physical unit root also at 1 after `1 - theta`;
            ## measured, the Arnoldi left the unit root at `1 - 1.5e-6`, past
            ## the `1e-6` deflation, so it reported the ORBIT TANGENT as the
            ## second multiplier and `Q ~ 6e5`.  Forming `M` by `m` matvecs
            ## and taking its eigenvalues directly gives the unit root to
            ## machine precision (it deflates cleanly) and the true second
            ## multiplier -- 8.59e-4 on van der Pol, matching Gear-2's
            ## 8.58e-4.
            _Md = np.column_stack([np.asarray(fp.matvec(_e), dtype=float)
                                   for _e in np.eye(n)])
            _lams = np.linalg.eigvals(_Md)
            _keep = np.real(_lams)[np.abs(_lams - 1.0) > 1e-6]
            if _keep.size == _lams.size:
                ## ⚠ NO MULTIPLIER IN THE WINDOW: for a MULTISTEP or trapezoidal
                ## solve the phase multiplier is 1 only as far as the
                ## discretisation is time-translation invariant, and a
                ## non-uniform `grid=` breaks that at O(h^2) (measured
                ## 1 - 5.1e-05 gear / 1 + 4.7e-05 trap at N = 400, 3:1, rate
                ## 4).  ⚠ NOT radau: its collocation solve keeps the
                ## multiplier at 1 + 1.2e-11 on the same grid (2026-09-20).  Without
                ## this the phase multiplier itself was reported as the
                ## SECOND one -- silently, f_amp 600x too small.  Drop the
                ## one nearest 1 instead; a uniform grid never gets here.
                _keep = np.real(np.delete(_lams, int(np.argmin(np.abs(_lams - 1.0)))))
            if _keep.size:
                lam2 = float(max(np.max(_keep), 0.0))
            ## the spectrum is exact, so there is nothing to certify against
            _resid, _certified = 0.0, True
        elif kk >= 2:
            ## ⚠⚠ THE TRUNCATED PATH, NOW GATED ON THE PAIR'S OWN RITZ
            ## RESIDUAL AND GROWN UNTIL IT CERTIFIES.  This branch runs only
            ## where the dense spectrum is unaffordable -- which is exactly
            ## where the truncation is least trustworthy, since a big circuit
            ## is the one likely to carry the many slow nodes that break the
            ## selection.  A fixed basis cannot work here: `k` has to track
            ## the slow-mode count, so `PPV_RITZ_BASIS = 16` was measured to
            ## be exact on one ladder and wrong on a longer one.  Doubling
            ## until the residual certifies is the same rule at every size.
            ##
            ## ⚠ THE LOOP TERMINATES ON THREE THINGS and only one of them is
            ## a threshold: the residual certifying, the basis reaching `n`
            ## (where the Arnoldi IS the spectrum), or the cost ceiling
            ## `PPV_RITZ_MAX_BASIS` -- which produces a WARNING and an
            ## uncertified number, never a silently wrong one.
            _budget = int(min(n, self.PPV_RITZ_MAX_BASIS))
            _resid = float('inf')
            while True:
                lam2, _resid = self._ritz_second_multiplier(fp, kk)
                if _resid <= self.PPV_RITZ_RESIDUAL_TOL or kk >= _budget:
                    break
                kk = int(min(2 * kk, _budget))
            _certified = _resid <= self.PPV_RITZ_RESIDUAL_TOL
            if not _certified:
                ## ⚠ AND IT WARNS ONLY WHEN IT FAILS.  An unconditional
                ## warning on every truncated call is noise a caller learns
                ## to ignore, which is worse than none: the point of the gate
                ## is that silence now MEANS something.
                warnings.warn(
                    'PSS.ppv: `second_multiplier` (%.6f) is NOT CERTIFIED. '
                    'n = %d exceeds FLOQUET_DENSE_LIMIT = %d, so it comes '
                    'from a truncated Arnoldi, and the selected pair\'s Ritz '
                    'residual is %.2e against a tolerance of %.0e after '
                    'growing the basis to %d (ceiling %d). Measured on a '
                    'ladder oscillator, an uncertified value is wrong by '
                    '4x-19x, in BOTH directions, and can exceed 1. Treat '
                    '`second_multiplier` and `Q` as indicative; '
                    "`info['second_multiplier_certified']` says which. "
                    'Raising PPV_RITZ_BASIS is measured NOT to be the fix -- '
                    'the basis has to grow with the slow-mode count, which '
                    'is what the loop above does; raise PPV_RITZ_MAX_BASIS '
                    'if the cost is acceptable.'
                    % (lam2, n, self.FLOQUET_DENSE_LIMIT, _resid,
                       self.PPV_RITZ_RESIDUAL_TOL, kk, _budget),
                    RuntimeWarning, stacklevel=2)
        if lam2 > self.PPV_SECOND_MULTIPLIER_WARN:
            warnings.warn(
                'PSS.ppv: a SECOND Floquet multiplier sits at %.6f, near '
                'the unit circle. The bordered extraction removes only the '
                'phase mode, so its conditioning degrades as that root '
                'approaches 1 -- measured losing six digits over six '
                'decades of time constant while every residual stayed at '
                '1e-11. ⚠ AND THE PHASE EQUATION ITSELF IS THE DEEPER '
                'ISSUE: it treats the frequency response as instantaneous, '
                'so slow nodes that FILTER a device\'s noise are not seen '
                'and phase noise is OVER-ESTIMATED. Neither a smaller '
                'tolerance nor a better extraction fixes that; it needs a '
                'frequency-aware PPV. Treat this result as an upper bound.'
                % lam2, RuntimeWarning, stacklevel=2)
        ## ⚠ ONE NUMBER THAT SUBSUMES FOUR DIAGNOSTICS.  An amplitude
        ## perturbation decays to `|lambda_2|` of its size each cycle, so
        ## the cycles needed to fall below a threshold IS the oscillator's
        ## Q: `Q = log(threshold)/log|lambda_2|` (Wang & Roychowdhury).
        ## The usual definitions do not apply to an autonomous circuit --
        ## `f_r/df` presumes a Bode plot of a BIBO-stable linear system,
        ## and stored/dissipated presumes damping a self-sustaining
        ## oscillator does not have.  Nor is it the resonator's Q.
        ##
        ## ⚠ AND IT IS THE SAME CONDITION AS EVERY FAILURE THIS CLASS
        ## WARNS ABOUT.  "High Q", "a second multiplier near 1", "slow
        ## amplitude restoration" and "a long time constant" are four
        ## vocabularies for one thing -- which is why the same circuits
        ## defeat the phase row, the eigen-split, the probe's continuation
        ## and the PPV's instantaneous-response assumption.  Not four
        ## coincidences.  It costs nothing here: the Arnoldi above already
        ## produced `|lambda_2|`.
        ##
        ## ⚠⚠ AND ITS NAME IS ONLY RIGHT WHILE THE OSCILLATOR'S AMPLITUDE
        ## MODE IS THE SLOWEST NON-UNIT MODE.  A parasitic with
        ## `tau_p/T > Q_osc` simply IS the second multiplier -- by
        ## definition, not by error -- and then this reports THE
        ## PARASITIC'S DECAY TIME IN PERIODS under the name `Q`.  MEASURED:
        ## at `tau_p/T` = 32 and 100 on a `Q = 16` oscillator, `lam2` is
        ## the parasitic and `Q` returns 32 and 100.
        ##
        ## ⚠ THE NUMBER IS RIGHT AND ITS NAME IS WRONG, which is why
        ## nothing misbehaves: every residual stays clean and the value is
        ## well converged.  A DCO's gated capacitor sits at
        ## `tau_p/T ~ 1e4`, i.e. permanently in that regime, so on exactly
        ## the circuits a hierarchical DCO method exists for, a reported
        ## `Q` would be the gated cap's RC in periods.  Read `Q` as
        ## "cycles for the SLOWEST NON-UNIT MODE to decay by 1/e", which is
        ## what it computes; it is the oscillator's Q only when that mode
        ## is the oscillator's.
        ##
        ## Reported for a `1/e` threshold, so `Q` is cycles-to-1/e.
        ##
        ## ⚠⚠ AND THIS LINE IS WANG & ROYCHOWDHURY'S IDENTITY
        ## `Q = log(threshold)/log|lambda_2|`, which does double duty and
        ## was shipped before either use was noticed.  It is what makes
        ## "bounded by Q" and "bounded by lambda_2" the SAME SENTENCE --
        ## the organising fact of this whole area, since a designer's
        ## objective (raise Q) IS the numerics' failure mode (lambda_2 ->
        ## 1).  It is also an ERROR AMPLIFIER:
        ##
        ##     (dQ/Q) / (dlambda_2/lambda_2)  =  -1/ln(lambda_2)  =  Q
        ##
        ## ⚠ SO THE RELATIVE ERROR IN `Q` IS `Q` TIMES THE RELATIVE ERROR
        ## IN `lambda_2`, and a caller reading `Q` at high Q is reading a
        ## quantity far less accurate than the multiplier behind it.
        ## MEASURED end to end on van der Pol tuned by `mu = 1/(2 pi Q)`,
        ## against the finest grid:
        ##
        ##     Q      npts   rel err lam2   rel err Q   ratio
        ##      3.18   120    1.04e-03      3.31e-03      3.2
        ##     15.92   120    5.75e-04      9.23e-03     16.1
        ##     63.66   120    4.89e-04      3.21e-02     65.7
        ##     63.66   480    7.56e-06      4.81e-04     63.7
        ##
        ## ⚠ THAT IS A RESOLUTION REQUIREMENT SCALING WITH `Q`, NOT A
        ## FIXED ACCURACY: 120 points/period gives `Q` to 0.3% at Q = 3
        ## and only 3.2% at Q = 64.  Payable here because Gear-2's
        ## `lambda_2` converges at better than second order (~8x per
        ## doubling); a method that BIASES `lambda_2` at fixed order has
        ## no such escape, and backward Euler's 5.6e-2 bias would become
        ## 85% in `Q` at Q = 100.
        Q = (-1.0 / np.log(lam2) if 0.0 < lam2 < 1.0 else float('inf'))
        info = {'border_residual': y,
                'tangent_border_residual': yf,
                'Q': Q,
                'null_residual': resid / max(float(np.linalg.norm(v)), 1e-300),
                ## ⚠ MULTIPLY `null_residual` BY THIS TO GET THE RELATIVE
                ## ERROR IN `v` THE RESIDUAL CANNOT EXCLUDE.  `null_residual`
                ## is `||v - M^T v|| / ||v||`, so an error component along the
                ## `lam2` left-eigendirection enters it scaled by `1 - lam2`
                ## and is nearly INVISIBLE exactly when `lam2 -> 1`.
                ## MEASURED on `_vdp_with_slow_node`, injecting a 1% error
                ## into a converged `v` (floor 4.6e-11):
                ##
                ##     lam2        r(random dir)   r(lam2 dir)   0.01*(1-lam2)
                ##     0.000856      1.65e-02       1.003e-02      9.99e-03
                ##     0.990049      1.65e-02       9.950e-05      9.95e-05
                ##     0.999900      1.65e-02       1.000e-06      1.00e-06
                ##     0.999999      1.65e-02       1.000e-08      1.00e-08
                ##
                ## Exact to every digit printed.  A RANDOM error is caught
                ## nine orders above the floor, so `null_residual` is a real
                ## gate and this module's assertions on it can fail -- but it
                ## loses sensitivity in the ONE direction that matters as the
                ## circuit gets better, which is the opposite of the
                ## reassurance a flat residual gives.
                ##
                ## ⚠⚠ THIS IS WHY A FLAT `null_residual` IS NOT EVIDENCE OF
                ## ACCURACY.  A residual that does not move while `lam2`
                ## sweeps toward 1 is not reporting that the answer stayed
                ## good; the bordered system is well-conditioned BY
                ## CONSTRUCTION, and the quantity it fails to see is
                ## precisely the one that grows.  Read the two numbers
                ## together or neither.
                'null_residual_amplification': (
                    1.0 / max(1.0 - lam2, np.finfo(float).eps)),
                'second_multiplier': lam2,
                ## Which route produced it, so a caller can tell an exact
                ## spectrum from a truncated estimate without re-deriving
                ## the rule.  See the branch above.
                'second_multiplier_route': ('dense' if _dense_ok
                                            else 'arnoldi'),
                ## The selected pair's Ritz residual and whether it cleared
                ## `PPV_RITZ_RESIDUAL_TOL`.  Exact on the dense route (0.0,
                ## True).  A caller that reads `Q` should read this too.
                'second_multiplier_residual': float(_resid),
                'second_multiplier_certified': bool(_certified),
                'q': q, 'xdot': xdot, 'tangent_pair': u,
                'samples': np.asarray(states),
                ## ⚠ `samples_eq` IS THE ONE TO CONTRACT `CY` AGAINST.
                ## `samples` is `C^T v_1` (a state perturbation's
                ## sensitivity); this is `v_1` (an equation-row input's).
                ## A noise current injected into a KCL row is the latter.
                'samples_pair': np.asarray(states_pair),
                'monodromy_method': getattr(self.par, 'method', '?'),
                'samples_eq': np.asarray(_eq),
                'v_eq': _v_eq,
                'times': np.asarray(fp.times, dtype=float),
                ## ⚠ THE PERIOD OF THE ORBIT THESE SAMPLES LIVE ON, which is
                ## not the caller's `period` when this came from a twin (trap
                ## and euler read a TR-BDF2 twin whose period differs by
                ## O(h^2)).  A quadrature over `times` divides by THIS; mixing
                ## it with the host's period was the E3 "sign change".
                'period': float(fp.T)}
        return v, info

    def frequency_aware_ppv(self, offset, tol=None):
        """The PPV at a nonzero modulation frequency (Lai 2008, eq. 23).

        ⚠ EQ. 23, NOT 24 (docs session, 2026-09-09): Lai's eq. 24 drops the
        AC columns of the Toeplitz block and is justified only "if we are
        only interested in the transfer functions when w_s is close to
        DC"; this shooting form has no such truncation -- `I - exp(-j w_s T)
        M^T` is the exact sampled LPTV adjoint at ANY offset, the monodromy
        already carrying the full time variation -- so it is eq. 23 for
        what the object IS, exact for the discretised system, and eq. 24
        only for the DC-reduction sentence pinned below.

        The classical PPV is the left null vector of `I - M^T`, bordered by
        `q = C(0) xdot(0)`; it is the phase response to a perturbation that
        is SLOW against every other Floquet mode.  This is the SAME
        bordered system at `alpha = exp(-j w_s T)`,

            [[I - alpha M^T,  q], [q^T, 0]] [v; y] = [0; 1],

        whose solution `v(w_s)` is the phase sensitivity to a perturbation
        modulated at `w_s`: at `w_s = 0` it IS `ppv()` (pinned), and away
        from it the AMPLITUDE mode admixes with weight
        `(1 - alpha)/(1 - alpha mu_2)` -- zero at DC, rising ten-fold per
        decade, cornering where `2 pi f_s T = 1 - mu_2` and flat above
        (docs session, 2026-09-08, on the slow-node fixture; the corner
        tracks the slow multiplier over two decades).  Verified at the
        source: Lai's eq. (24) at `w_s = 0` "is the augmented PPV
        extraction equation", verbatim; his construction is harmonic
        balance, this is the shooting basis, and the two agree on what the
        object is.

        Returns `(v, info)`: `v` the pair-space anchor vector (complex);
        `info['samples_pair']` the T-periodic envelope over the period,
        `lambda_k exp(+j w_s t_k)` with `lambda` the transposed replay of
        `v` -- the sideband rows' own convention, so its Fourier
        coefficient at harmonic `k` is the phase transfer of a source band
        at `k f0 + f_s`; `info['samples']` the same in state space
        (`C^T v`), SECOND order in the step through `ppv()`'s own
        pair-consistent propagation run on the complex anchor
        (`_ppv_propagate`, lifted 2026-09-08; at w_s = 0 it is `ppv()`'s
        `samples`); `info['admixture']` the norm fraction of `v(w_s)`
        orthogonal to the DC PPV; `info['corner']` the predicted corner
        `|1 - mu_2| / (2 pi T)` in Hz; `info['alpha']`; `info['ppv']` the
        DC object's info.

        ⚠ WHAT IT IS FOR.  A source that reaches the phase through a slow
        path (an RC leg, tau >> T) is filtered at its own corner, and the
        DC PPV cannot see that; the harmonic sum built from THIS object's
        coefficients carries the filter inside `c_k(w_s)` with no explicit
        model of the path (A2, roadmap).  MEASURED (2026-09-08): the ratio
        `sum_k |c_k(w_s)|^2 / sum_k |c_k(0)|^2` reproduces `pnoise`'s
        `S_pm/(4 S_v)` for a source behind the slow node within 0.6 % to
        r = 1e-2 and 2 % at 5e-2 (a = 0.4, loss 0.2, tau/T = 100), where the
        DC sum with the filter by hand was 4 % off; at r = 0.1 the two
        differ by -6 %, unchanged at twice the grid -- a gap between PM by
        sideband quadrature and phase-mode projection, both 1e-3 of DC
        there -- and put to a nonlinear Monte Carlo that instantiates
        neither construction (2026-09-09; 4 seeds x 10 000 periods per
        point, phase read two ways).  ⚠ A first reading at one asymmetry
        (a = 0.25) assigned each construction to ONE phase definition
        crosswise (crossing phase to S_pm at 1.027, demodulated phase to
        this sum at 0.977, each +-1.8 %); the asymmetry sweep a = 0.25 /
        0.12 / 0.05 / 0 REPLACED it.  Double ratios (Monte Carlo slow/core
        over predicted slow/core, each estimator calibrated on the core):

            a      crossing vs S_pm / this sum   demod vs S_pm / this sum
            0.25   1.032 / 1.066  (+-1.8 %)      0.954 / 0.985
            0.12   1.093 / 1.121  (+-3.5 %)      0.974 / 0.999
            0.05   1.122 / 1.149  (+-4 %)        1.004 / 1.027
            0.00   1.134 / 1.162  (+-4 %)        1.027 / 1.052

        ⚠⚠ AT 16 SEEDS PER POINT (Andreas, same day; +-1.5-1.7 %):

            a      crossing vs S_pm / this sum   demod vs S_pm / this sum
            0.25   1.038 / 1.072                 0.957 / 0.989
            0.12   1.085 / 1.113                 0.979 / 1.004
            0.05   1.117 / 1.143                 1.011 / 1.035
            0.00   1.129 / 1.156                 1.036 / 1.061

        The statistic is the SLOPE in a, not any one point (the docs
        session's framing): demod vs this sum -0.27 +- 0.09 per unit a
        (3.2 sigma), vs S_pm 3.6 sigma; crossing 4.2 / 4.7 sigma.  So the
        ratio is NOT constant in a at ~3 sigma for the demodulated phase
        and above 4 for the crossing: NEITHER construction describes
        EITHER measured phase across the range.  The two constructions
        track each other to 1 % over the sweep while both estimators --
        a point sample at a crossing and an average over a period -- drift
        together, in the same direction, against both.  The a = 0.25
        agreement of the demodulated phase with this sum is where its
        curve crosses the sum, not a match.  What the crossing carries
        beyond that: two thirds of its excess is waveform content beyond
        0.5 f0 from the carrier (an instantaneous crossing aliases the
        additive noise a one-period demodulation cannot see), the rest
        the demod's own boxcar loss (1/sinc^2 = 1.045 at the band centre)
        plus a common amplitude-to-crossing gain of ~0.6 (band-limited
        coherence).  The common drift of BOTH estimators against BOTH
        linear constructions as a -> 0 is the open object; candidate, a
        second-order amplitude-to-phase conversion the linear theory
        cannot contain -- REFUTED the same evening (PSD/4: the slow
        fixture scales linearly, 0.993 +- 0.013; grid doubling moves the
        constructions < 0.7 %; band conventions identical).  RESOLVED by a
        forward tone-transient route on the MC's own discretised system
        (no adjoint, no sideband assembly): pnoise's S_pm agrees with the
        forward LPTV PM sidebands to 1 % at both asymmetries, and the
        estimators' a-dependent double ratios are REPRODUCED by that
        deterministic linear route (demod 0.938 -> 1.003 against the MC's
        0.957 -> 1.036).  The drift is the ESTIMATORS: a one-period
        fundamental demodulation leaks the other harmonics' sidebands
        through its boxcar (sinc(pi(1-r)) ~ 0.1 for the second harmonic's,
        which is ~a), zero crossings convert every harmonic's; both read a
        given sideband PM with a fixture-dependent gain.  This object and
        S_pm are PM by quadrature of the FUNDAMENTAL'S sidebands; compare
        them with that, not with a demodulated or crossing phase.  ⚠ AND
        COMPARE BAND WITH BAND: a source behind the slow node has an
        in-band spectrum that is not 1/r^2 (its slow/core ratio swings
        1.16 -> 0.87 across 0.08-0.15 f0 at a = 0.25), so a band mean and a
        point value differ by ~4 % there; a 2 % constant between the noisy
        Monte Carlo and the deterministic route remains after that, within
        the excursion-amplitude bound, unresolved.  ⚠ The
        premise "a -> 0 makes the definitions coincide" was wrong: the
        asymmetry removes even harmonics only, van der Pol's third stays
        at 9.7 % of the fundamental, and the construction gap GROWS as
        a -> 0 (core 1.010 -> 1.055); the sinusoidal limit is mu -> 0.  The
        instrument hypothesis (spectrum analyser <-> this sum, time-
        interval analyser <-> S_pm) is refuted in its crossing half.  The slow multiplier's own coefficient
        (`mode_content[0]`) corners at 1.6e-3 f0 for tau/T = 100 with a
        plateau of 2.45e-6 (the docs session's 2.29e-6), scaling as
        T/tau.  ⚠ DO NOT GATE ON `|v|`: with
        `q^T v = 1` the `1/(1 - alpha)` pole cancels between numerator
        and denominator and the norm is frequency-flat by construction; a
        one-percent orthogonal admixture moves it by 5e-5.  The change is a
        DIRECTION -- read `admixture`, or the per-harmonic coefficients.
        """
        import scipy.sparse.linalg as spla
        v0, info0 = self.ppv(tol)
        fp = self.factored_period()
        m = self.cir.n - 1
        n = fp.width
        irn = self.irefnode
        T = float(fp.T)
        alpha = np.exp(-2j * np.pi * float(offset) * T)
        q = np.asarray(info0['q'], dtype=float)
        qp = np.concatenate((q, np.zeros(n - m))).astype(complex)

        def _mv(z):
            z = np.asarray(z, dtype=complex)
            v_, y_ = z[:n], z[n]
            top = v_ - alpha * np.asarray(fp.matvec_transposed(v_), dtype=complex) + y_ * qp
            return np.concatenate((top, [complex(qp @ v_)]))
        rtol = max(self.par.reltol * 1e-2 if tol is None else tol, 1e-14)
        rhs = np.zeros(n + 1, dtype=complex)
        rhs[n] = 1.0
        z, relres, _H, _k = _arnoldi_gmres(_mv, rhs, rtol=rtol, maxiter=min(n + 1, 200))
        if relres > max(1e3 * rtol, 1e-8):
            ## Lai's own warning about this object: eq. 23 "is very
            ## difficult to solve using iterative solvers (such as GMRES),
            ## because the extra columns and rows from the Toeplitz block
            ## degrade the block diagonal preconditioner" -- a property of
            ## the formulation, not a defect of the fixture
            raise RuntimeError(
                'PSS.frequency_aware_ppv: the bordered solve at offset %g did '
                'not converge (relative residual %.3e) -- a known property of '
                'the bordered LPTV adjoint (Lai 2008), not a fixture defect.'
                % (float(offset), relres))
        v = np.asarray(z[:n], dtype=complex)
        ## `ppv()`'s own normalisation, `v . xdot(0) = 1` -- the border fixes
        ## `q^T v = 1` only, and the two differ by a factor AND a sign on
        ## the van der Pol (measured 2.07 there; -0.94 on the slow-node
        ## fixture): at alpha = 1 this is what makes the object `ppv()`.
        xdot = np.asarray(info0['xdot'], dtype=float)
        vx = complex(v[:m] @ xdot)
        if vx == 0.0:
            raise ValueError('PSS.frequency_aware_ppv: v(w_s) is orthogonal to the orbit tangent.')
        v = v / vx
        ## the envelope over the period, in the sideband rows' convention:
        ## `ppv()`'s own second-order propagation on the complex anchor
        ## (lifted into `_ppv_propagate`), then the per-step phase
        irn = self.irefnode
        x0r = np.asarray(self._period_state[1], dtype=float).ravel()
        x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))
        _alg_rows, _alg_cols = self._algebraic_adjoint_pattern(x0f)
        states, states_pair, _ts, _Xf = self._ppv_propagate(fp, v, m, xdot, _alg_rows, _alg_cols)
        st = np.asarray(states_pair, dtype=complex)
        tms = np.asarray(fp.times, dtype=float)[:st.shape[0]]
        phase = np.exp(2j * np.pi * float(offset) * tms)
        samples_pair = st * phase[:, None]
        samples = np.asarray(states, dtype=complex) * phase[:, None]
        ## the admixture: what of v(w_s) is NOT along the DC PPV.  ⚠ This is
        ## dominated by whichever mode has the largest weight, which on a
        ## core with a fast amplitude mode is THAT one (weight ~ 2 pi r,
        ## cornering at r ~ (1 - lambda_2)/2pi ~ 0.16 for the van der Pol),
        ## so a slow node's own admixture -- 1e-6 to 1e-2 of it -- is
        ## invisible in the norm.  `mode_content` reads each mode's own
        ## coefficient (dense route below FLOQUET_DENSE_LIMIT).
        v0c = np.asarray(v0, dtype=complex)
        v0n = v0c / np.linalg.norm(v0c)
        proj = np.vdot(v0n, v) * v0n
        admixture = float(np.linalg.norm(v - proj) / np.linalg.norm(v))
        mode_content = None
        multipliers = None
        if n <= self.FLOQUET_DENSE_LIMIT:
            Md = np.column_stack([np.asarray(fp.matvec(e), dtype=float) for e in np.eye(n)])
            mu, P = np.linalg.eig(Md.T)
            order = np.argsort(-np.abs(mu))
            mu, P = mu[order], P[:, order]
            a = np.linalg.solve(P, v)
            mode_content = np.abs(a[1:]) / abs(a[0])
            multipliers = mu
        mu2 = float(info0.get('second_multiplier', 0.0))
        info = {'alpha': alpha, 'offset': float(offset), 'admixture': admixture,
                'mode_content': mode_content, 'multipliers': multipliers,
                'corner': abs(1.0 - mu2) / (2.0 * np.pi * T), 'second_multiplier': mu2,
                'samples_pair': samples_pair, 'samples': samples, 'times': tms,
                'residual': float(relres), 'ppv': info0}
        ## ⚠ THE EQUATION-ROW SAMPLES, which is what `CY` contracts against
        ## (see `ppv()`'s `samples_eq`).  `_equation_row_ppv` casts to float,
        ## so the complex samples go through it as real and imaginary parts:
        ## the map is linear (measured to round-off) and commutes with the
        ## per-step phase factor.  At offset 0 these equal `ppv()`'s
        ## `samples_eq` to 7e-14.  ⚠ `times` above is truncated to the sample
        ## count, ONE ENTRY SHORT of the orbit's grid -- a quadrature over it
        ## drops the last step (0.4 % of `c` on an asymmetric orbit, measured);
        ## integrate over `info['ppv']['times']` and `['period']` instead.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            info['samples_eq'] = np.array([
                self._equation_row_ppv(
                    np.real(st[:m]), _Xf[:, _j if _j < _Xf.shape[1] else -1],
                    _alg_rows, _alg_cols)
                + 1j * self._equation_row_ppv(
                    np.imag(st[:m]), _Xf[:, _j if _j < _Xf.shape[1] else -1],
                    _alg_rows, _alg_cols)
                for _j, st in enumerate(samples)])
        info['period'] = T
        return v, info

    FLOQUET_DENSE_LIMIT = 400
    ## below this a multiplier is an annihilated algebraic
    ## direction, not a mode -- see `floquet_modes`
    FLOQUET_NULL_TOL = 1e-12

    def floquet_modes(self, pss_unused=None, nmodes=None, fp=None):
        """⚠ THE MODES' ACCURACY IS THE METHOD'S -- AND GEAR'S ADJOINT MODES
        WERE FIRST ORDER UNTIL 2026-09-20, from a second-order method, on
        every grid.  The reconstruction of `q` from the discrete adjoint took
        the pair's first block through pinv(C^T), a fraction of a step off
        the node (see the note at the reconstruction below).  The invariant
        `q^T C p`, constant along the orbit for the true adjoint, measured on
        the asymmetric van der Pol at N = 200 / 400 / 800:

            gear, uniform, before   9.9e-03  5.0e-03  2.5e-03   halving: FIRST order
            gear, uniform, now      9.3e-04  2.3e-04  5.6e-05   quartering: SECOND
            gear, 3:1 grid, now     2.2e-02  1.1e-02  5.9e-03   halving (was 3.4e-02):
                                    a variable-step multistep adjoint is first
                                    order and no rescaling makes it more

        On a uniform grid at N = 400 gear's modal-spectrum PARTS now agree
        with trap's (second order) to 4e-4 / 9e-4 / 5e-4 and the total closes
        on pnoise at 1.0004 (it read 1.0145 before).  ⚠ A mode's k = 1
        Fourier coefficient still differs from a uniform N = 3200 reference
        by 3e-2 at N = 200 -- for `p` AND `q` alike, so it is the mode's
        phase across N, not the adjoint; the invariant and the spectra are
        phase-insensitive and are the gates.

        Radau: on a 3:1 non-uniform grid its modes reproduce the uniform
        grid's to 1e-10 and its phase multiplier stays at 1 + 1e-11; gear and
        trap leave it at O(h^2) there and the modal spectra refuse.  For the
        modes on a non-uniform grid use radau (the default).

        The Floquet pairs `(λ_l, μ_l, p_l(t), q_l(t))` — A9's prerequisite.

        Returns a list of dicts, one per mode, ordered by `|λ|` descending.
        `nmodes=None` returns EVERY non-null mode, and that default is the
        requirement rather than a convenience:

        ⚠⚠ ALL OF THEM ARE REQUIRED, BY THE SOURCE. Traversa & Bonani, IET
        CDS 2011: "The calculation of orbital fluctuations and of the
        phase-orbital correlation within Floquet-based noise analysis of
        autonomous systems requires the availability of ALL the direct and
        adjoint Floquet eigenvectors associated with the noiseless limit
        cycle."  (Cited, not verified here; relayed from the paper.)  An
        earlier cost estimate for A9 -- "a few more Floquet pairs" -- was
        relayed without checking it against that sentence, and is wrong.

        ⚠ TRUNCATION IS LEGITIMATE ONLY WITH A BOUND. Traversa & Bonani,
        TCAD 2013, compute a CHOSEN number of exponents and both
        eigenvector sets for the linearisation of index-1 DAEs around a
        limit cycle -- this formulation -- with the error "proved to tend
        to zero along with the ratio between the norms of the NEGLECTED
        AND RETAINED ROWS".  So passing `nmodes` is allowed, but a caller
        who does owes that ratio as the gate; this routine does not
        compute it.  The dense route below returns everything anyway, so
        at the sizes it serves the question does not arise.


            lam    the Floquet MULTIPLIER, eigenvalue of the monodromy
            mu     the Floquet EXPONENT, `log(λ)/T` (complex)
            u0,v0  right and left eigenvectors at `t = 0`, biorthonormal
                   (`v_k† u_l = δ_kl`)
            p      `p_l(t_j) = Φ(t_j,0) u_l(0) · exp(−μ_l t_j)` — the
                   T-PERIODIC part, sampled on the PSS grid
            q      the adjoint counterpart from the reverse replay
            times  the grid `p` and `q` are sampled on

        ⚠⚠ **WHY THIS EXISTS: `S_yy` NEEDS THE EIGENVECTORS OVER THE
        PERIOD, NOT JUST THE EXPONENTS.** Traversa & Bonani (TCAS-I 2011)
        Lemma 3.5 makes the orbital spectrum a sum of Lorentzians centred
        at `jω₀ + Im{μ_l}` with half-width `|Re{μ_l}| + ½h²ω₀²c`, weighted
        by `C_lhj` (their eq 22) — and `C_lhj` is built from the FOURIER
        COEFFICIENTS of `u_l(t)` and of `v_l(t)ᵀ B(t)`. Their own text is
        explicit that the exponents alone do not order the result: *"a
        major role in the C and D coefficients is also played by the
        Floquet eigenvectors, which could determine large orbital
        fluctuations contributions even when the Floquet exponents are not
        near to zero."* So `|λ₂|` — the only mode information this class
        used to expose — is not sufficient, by the source's own statement.

        ⚠ THE PERIODIC PART IS THE OUTPUT, NOT `Φ(t,0)u(0)`. Floquet's
        theorem says the solution is `p_l(t)exp(μ_l t)` with `p_l`
        T-periodic; the raw propagated vector is not periodic and its
        Fourier series is not the one eq (22) wants. Dividing out
        `exp(μ_l t)` is what makes `p_l(T) = p_l(0)` — which is also the
        gate below, and the only check here that needs no reference.

        ⚠ DENSE, AND REFUSED ABOVE `FLOQUET_DENSE_LIMIT`. The monodromy is
        assembled column by column (`n` matvecs) and diagonalised. That is
        honest for the sizes this is useful at and wrong to hide at larger
        ones: an Arnoldi route would return Ritz VECTORS rather than only
        the Ritz values `ppv()` currently keeps, and is the extension.

        ⚠⚠ AND THAT EXTENSION CONVERGES TO THE PHYSICAL MODE *LAST*, WORSE
        AS Q RISES -- a structural fact, not a measurement (Garcia, Romero
        & Acha 2022, read firsthand by the docs session).  Arnoldi resolves
        the LARGEST-magnitude eigenvalues first; the Ritz route works on
        `A = I - M` and recovers `lam = 1 - theta`, so the physical
        `lam_2 -> 1` maps to `theta_2 -> 0`, the SMALLEST, while the fast
        parasitic modes (`lam ~ 0`) sit at `theta ~ 1` and are resolved
        first.  The separation to resolve is `1/theta_2 ~ Q_lambda`: 3.7,
        16.4, 64.5, 128.5 at `Q_lambda` = 3.18, 15.9, 64, 128.  The
        difficulty scales with the very quantity being measured.  ⚠ This
        is a DIFFERENT Krylov problem from B13's, which measured GMRES
        iterations for the bordered SOLVE `(I - M) w = b` and found them
        independent of Q -- solving a system and extracting its smallest
        eigenvalue are not the same question, and B13 says nothing about
        the second.  The paper is a sound source for the method and was
        validated on power networks, not RF oscillators, so it reports no
        evidence either way about the high-Q regime; and it states that a
        truncated run "cannot compute ALL the Floquet multipliers" -- which
        is the requirement eq (22) carries (IET CDS 2011, above).
                """
        _tw = self.monodromy_twin()
        if _tw is not self:
            ## ⚠ the twin gets `None`, not `pss_unused`: it must read ITS OWN
            ## factored period, which is the whole reason a twin exists.  An
            ## explicit `fp` from the caller still wins.
            return _tw.floquet_modes(None, nmodes, fp)
        ## ⚠⚠ `pss_unused` IS IGNORED, AS ITS NAME SAYS -- and it used to be
        ## DEREFERENCED here, so the documented default call `floquet_modes()`
        ## raised `AttributeError: 'NoneType' object has no attribute
        ## 'factored_period'` on every method (measured on radau, trap and a
        ## GLM alike).  Every call site inside this file passes `self`, so
        ## reading `self` changes no existing answer and makes the no-argument
        ## call work; the parameter stays in the signature because callers
        ## pass it positionally.
        fp = self.factored_period() if fp is None else fp
        n = fp.width
        T = float(fp.T)
        if n > self.FLOQUET_DENSE_LIMIT:
            raise NotImplementedError(
                'PSS.floquet_modes: the monodromy is assembled densely and '
                'this one is %d wide, past the %d limit. The extension is '
                'an Arnoldi that keeps its Ritz VECTORS -- ppv() already '
                'builds the basis and discards them.'
                % (n, self.FLOQUET_DENSE_LIMIT))

        M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                             for e in np.eye(n)])
        ## on a staged solve the map is the TOTAL one (2026-09-22): the
        ## crossings move with the state, `M + P_theta dtheta/dx_0`
        _ev_fm = EventColumns.of(self, n)
        if _ev_fm is not None:
            M = _ev_fm.total_matrix(M)
        lam, U = np.linalg.eig(M)
        lam_l, V = np.linalg.eig(M.T)

        order = np.argsort(-np.abs(lam))
        lam, U = lam[order], U[:, order]
        ## pair each right eigenvalue with its left partner by value
        pair = []
        used = set()
        for k in range(len(lam)):
            d = np.abs(lam_l - lam[k])
            for j in np.argsort(d):
                if j not in used:
                    used.add(int(j))
                    pair.append(int(j))
                    break
        V = V[:, pair]

        times = np.asarray(fp.times, dtype=float)
        out = []
        ## ⚠ NULL MODES ARE DROPPED, NOT RETURNED WITH A BAD RESIDUAL. A
        ## DAE's monodromy has exact zeros (the algebraic directions the
        ## step map annihilates); their "eigenvectors" are arbitrary, the
        ## exponent `log(0)` does not exist, and the periodic part comes
        ## back as noise -- measured, residual 0.56 and periodicity 8.9
        ## against 1e-15 for the physical pair. Returning them invites a
        ## caller to average over a mode that means nothing.
        keep = [k for k in range(n) if abs(lam[k]) > self.FLOQUET_NULL_TOL]
        if nmodes is not None and int(nmodes) < len(keep):
            ## ⚠⚠ TRUNCATING BY MULTIPLIER MAGNITUDE IS REFUTED BY THE SOURCE'S
            ## OWN WORKED EXAMPLE.  Traversa & Bonani TCAS-I 2011 Sec. V, on
            ## their Colpitts: "six orders of magnitude separate mu_2 and mu_3,
            ## while the corresponding contribution to orbital noise are not in
            ## the same ratio.  Rather, far from the oscillator harmonics, the
            ## contribution of mu_3 is dominant with respect to mu_2".  The
            ## ordering INVERTS.  Four statements agree: eq (8)'s sum over
            ## k = 2..n (structural), p.4 (asserted), Sec. V (measured on a
            ## real circuit), this repo's concentration sweep (m/n = 0.97).
            ## The caller who truncates owes the dropped weight as a gate;
            ## this says so at the call rather than only in the docstring.
            warnings.warn(
                'PSS.floquet_modes: nmodes=%d keeps %d of %d non-null modes, '
                'selected by |lambda|. Orbital-noise weight does NOT follow '
                'multiplier magnitude -- Traversa & Bonani (TCAS-I 2011, '
                'Sec. V) show the contribution ordering INVERTING across six '
                'orders in mu on their Colpitts, and this repo measured no '
                'concentration (m/n = 0.97). A covariance or spectrum built '
                'from a truncated set is missing weight you have not bounded; '
                'pass nmodes=None for all modes.'
                % (int(nmodes), int(nmodes), len(keep)),
                RuntimeWarning, stacklevel=2)
        ## ⚠ `None` means ALL non-null modes -- the default since the IET CDS
        ## 2011 correction -- and it used to fall into `int(None)` here because
        ## the only test passed a number. A default nobody exercises is not a
        ## default.
        for k in (keep if nmodes is None else keep[:int(nmodes)]):
            lk = complex(lam[k])
            uk, vk = U[:, k].astype(complex), V[:, k].astype(complex)
            nrm = complex(np.vdot(vk, uk))
            if abs(nrm) < 1e-30:
                raise ValueError(
                    'PSS.floquet_modes: mode %d has left and right '
                    'eigenvectors orthogonal to each other (v.u = %.3e), so '
                    'it cannot be biorthonormalised. That happens at a '
                    'defective eigenvalue -- two multipliers have collided.'
                    % (k, abs(nrm)))
            vk = vk / np.conj(nrm)                     ## v_k† u_k = 1
            muk = np.log(lk) / T

            ## ⚠ EVERYTHING BELOW IS THE WIDTH-`m` STATE BLOCK, NOT THE
            ## WIDTH-`n` MAP INPUT. Under a solved-history map `n = 2m`
            ## and the second block is the history term, not a second
            ## state; the replays collect `m`-wide states either way, and
            ## `u_l(t)` in eq (22) is a state-space function. Mixing the
            ## two is a shape error that surfaces three frames away.
            m = self.cir.n - 1

            ## forward: Phi(t_j,0) u_k(0), by an UNFORCED driven replay
            zero = np.zeros(m)
            _end, fwd = self._forced_replay(fp, 0.0, zero, y0=uk, collect=True)
            traj = ([np.asarray(uk, dtype=complex)[:m]]
                    + [np.asarray(z, dtype=complex).ravel()[:m] for z in fwd])
            tt = times[:len(traj)]
            traj = traj[:len(tt)]
            p = np.column_stack([traj[j] * np.exp(-muk * tt[j])
                                 for j in range(len(traj))])

            ## adjoint: Phi(T,s_j)^T v_k(T) -- B8 made this available under
            ## every integrator, not only the solved-history one
            ## ⚠⚠ THE ADJOINT IS THE PER-STEP TRANSPOSED SOLVE `t`, NOT THE
            ## PAIR'S FIRST BLOCK -- AND IT BELONGS TO THE NEXT NODE
            ## (2026-09-20, measured).  `collect` hands back both: `ts[k]`,
            ## the solve `Jf_k^-T w1` made while replaying step k backwards,
            ## and `states[k]`, the pair (w1; w2) it leaves behind.  This
            ## took the pair's first block and mapped it through pinv(C^T);
            ## since `w1 = Jf^T t = (a0 C + G)^T t` and the adjoint equation
            ## `C^T dq/dt = G^T q` turns the `G^T t` part into a time
            ## derivative, that `q` was `a0 * q(t + 2h/3)` -- staggered by a
            ## fraction of a step, so the invariant `q^T C p` drifted along
            ## the orbit by 1e-2 at N = 200 and HALVED per doubling: FIRST
            ## order, from a second-order method, on every grid.  The solve
            ## `t` itself obeys the BDF2-discretised adjoint recursion, and
            ## the replay computes it for step k from the pair at node k+1,
            ## so `ts[k]` is the adjoint at node k + 1.  Scored on the
            ## invariant's spread at N = 200 / 400 / 800, uniform grid:
            ##
            ##     pair block, pinv(C^T)  (this, before)   9.9e-03 5.0e-03 2.5e-03   x2 per doubling
            ##     ts[k] at node k        (one node off)   1.5e-02 7.6e-03 3.8e-03   x2
            ##     ts[k] at node k + 1    (this, now)      9.3e-04 2.3e-04 5.6e-05   x4  SECOND ORDER
            ##
            ## ⚠ THE SCALE.  `t = Jf^-T w1` carries the step through `a0 ~ 1/h`,
            ## invisible on a uniform grid (absorbed by `c0` below) and a
            ## factor-3 modulation on a 3:1 one, so `a0` of the node's own step
            ## is put back.  ⚠ THAT IS FIRST ORDER ON A NON-UNIFORM GRID and
            ## cannot be more: the discrete adjoint of a variable-step
            ## multistep method draws `a1` and `a2` from LATER steps, so its
            ## recursion is the continuous adjoint's only to O(h) once the step
            ## changes (Sandu's inconsistency).  Four scalings were measured
            ## on the 3:1 grid and all halve per doubling; this one is the
            ## best of them at 2x the previous accuracy.  A non-uniform gear
            ## grid is refused by the modal spectra anyway (its phase
            ## multiplier leaves the unit circle); radau is exact there.
            ## Node 0 is node N by periodicity of the periodic part.
            ## ⚠ GEAR ONLY (`solved_history`).  A one-step kind's `ts` is
            ## NESTED -- per-stage solves per step -- and its state-block
            ## adjoint through pinv(C^T) was measured exact (radau, 1e-10 on
            ## a 3:1 grid) and second order (trap); those keep their path.
            _e2, _tsolves, _st = fp.matvec_transposed(
                vk, collect=True, inject=self._event_costate_injection(fp, vk, n))
            _gear_pair = getattr(fp, 'is_pair', False)
            ## ⚠ GEAR ON A NON-UNIFORM GRID: THE TRANSPOSE IS FIRST ORDER AND
            ## NO RESCALING LIFTS IT, so the adjoint is integrated SEPARATELY
            ## there -- gated on STRUCTURE (a multistep pair on a grid whose
            ## step changes), never on order.  See `_continuous_adjoint`.
            _nonuniform = self._period_quadrature(fp) is not None
            if _gear_pair and _nonuniform:
                q, ts2 = self._continuous_adjoint(fp, lk, muk, times)
            elif _gear_pair:
                _tsolves = [np.asarray(z, dtype=complex).ravel()[:m]
                            for z in _tsolves]
                _a0 = [float(np.asarray(_step[2][0])) for _step in fp.steps]
                _nq = min(len(_tsolves), len(times) - 1)
                ts2 = times[:_nq + 1]
                qtraj = [None] * (_nq + 1)
                for _k in range(_nq):
                    qtraj[_k + 1] = (_a0[_k] * _tsolves[_k]
                                     * np.exp(muk * ts2[_k + 1]))
                qtraj[0] = qtraj[_nq]
                q = np.column_stack(qtraj)
            else:
                qtraj = ([np.asarray(z, dtype=complex).ravel()[:m] for z in _st]
                         + [np.asarray(vk, dtype=complex)[:m]])
                ts2 = times[:len(qtraj)]
                qtraj = qtraj[:len(ts2)]
                q = np.column_stack([qtraj[j] * np.exp(muk * ts2[j])
                                     for j in range(len(qtraj))])

            ## ⚠⚠⚠ THE REPLAYED VECTOR IS `C^T q`, NOT `q`.  The conserved
            ## bilinear form of the variational DAE is `w^T C delta`, so over
            ## a period `M_a^T C M = C`, which makes the LEFT eigenvector of
            ## the state monodromy `C^T w(0)` -- the adjoint mode in the
            ## "left-eigenvector coordinates", one factor of `C^T` away from
            ## the state-space adjoint `q` that eq (22) and every covariance
            ## here need.  The transposed replay propagates that object, so
            ## every sample of `q` above is `C(t)^T q_true(t)`.
            ##
            ## ⚠ INVISIBLE ON EVERY FIXTURE THIS REPO HAD, for a geometric
            ## reason: van der Pol's reduced `C` is `diag(1, -1)`, and at
            ## `t = 0` the orbit sits at `[2, 0]` where the adjoint is nearly
            ## axis-aligned, so `C^T q` and `q` point the same way up to sign
            ## (`|cos| = 0.9972`).  On an ASYMMETRIC orbit the seed is off-axis
            ## and the two separate -- measured `|cos(v_k, q_true)| = 0.5738`
            ## at `a = 0.30` on van der Pol + `a u^2` -- while the two adjoints
            ## there are nearly PARALLEL (`|cos(q_2, q_1)| = 0.997`), so the
            ## wrong vector is mostly phase adjoint.  Result: the orbital
            ## covariance was 81x LOW against a Monte Carlo (0.0123 of the
            ## truth), and `|cos(C^-T v_k, q_true)| = 1.0000` at both
            ## asymmetries.  Applying `C^-T` here takes it to 1.06 of the
            ## Monte Carlo at `a = 0.30` and 1.0004 at `a = 0`.
            ##
            ## ⚠ THIS ALSO DISSOLVES THE "q IS STORED IN REVERSE TIME"
            ## finding recorded the same day: with the right vector,
            ## `q(t)^T C p(t)` is conserved at the SAME index (4.2e-04) and
            ## NOT the reversed one (2.0).  `diag(1, -1)` flips one
            ## component, which on a half-wave symmetric orbit is exactly the
            ## relation between `q(t)` and `q(T - t)` -- a sign flip read as a
            ## time reversal.
            ##
            ## ⚠ Per sample, because `C` may depend on the state.  `pinv`
            ## rather than `inv` so a singular reduced `C` (an index-2 MNA,
            ## algebraic rows) does not raise; the algebraic components of `q`
            ## are then the minimum-norm choice, which is a SCOPE LIMIT and
            ## not a solution -- recorded, not hidden.
            ## The pinv(C^T) map belongs to the STATE-BLOCK adjoint of the
            ## one-step kinds; gear's transposed solve is already the adjoint
            ## of the DAE variable (its invariant is `q^T C p`, see above).
            if not _gear_pair:
                _Wq = np.delete(np.asarray(self.waveform[1], dtype=float),
                                self.irefnode, axis=0)
                _nw = _Wq.shape[1]
                for _j in range(q.shape[1]):
                    _Cj = np.asarray(self._C_at(_Wq[:, min(_j, _nw - 1)]),
                                     dtype=float)
                    q[:, _j] = np.linalg.pinv(_Cj.T) @ q[:, _j]

            ## ⚠⚠ RENORMALISE ON THE STATE BLOCK. `v_k` was biorthonormalised
            ## against `u_k` at the map's FULL width `n`; under a
            ## solved-history map that is the pair `[x_n; x_{n-1}]`, and the
            ## width-`m` state block then carries `q(0)^T p(0) = c0 != 1`.
            ## MEASURED on van der Pol under gear: c0 = 1.324143, constant
            ## around the cycle to four digits -- and the orbital covariance
            ## assembled from these parts came out too large by EXACTLY
            ## c0^2 = 1.7535 against two independent routes, because `q`
            ## enters it quadratically. The periodicity gate p(T) = p(0)
            ## cannot see this: periodicity is scale-free. On the plain path
            ## n = m and c0 = 1, so this is a no-op there. The adjoint takes
            ## the scale (the right vector is the physical direction).
            ## ⚠⚠ THE INNER PRODUCT IS `C`-WEIGHTED, AND THE UNWEIGHTED ONE
            ## WAS WRONG BY A FACTOR OF `C` -- INVISIBLE ON EVERY FIXTURE
            ## THIS REPO HAD.  The conserved bilinear form of the variational
            ## DAE is `q(t)^T C(t) p(t)`, not `q(t)^T p(t)`: differentiating
            ## `G p + d(C p)/dt = 0` against the adjoint gives
            ## `d/dt [q^T C p] = 0`, so `q^T C p` is the invariant and the
            ## biorthonormality that eq (22) assumes is `q_k^T C p_l = d_kl`.
            ##
            ## ⚠ ON A UNIT-REACTANCE FIXTURE THE TWO ARE THE SAME NUMBER,
            ## which is exactly why this survived: van der Pol with
            ## `c = L = 1` gives `q^T C p = 0.9992` against `q^T p = 1`.
            ## Sweep the capacitance at fixed `w0` and the two separate --
            ## MEASURED `q^T C p` = 0.2495 / 0.9992 / 3.9982 at
            ## `C` = 0.25 / 1 / 4, i.e. exactly `C`, while `q^T p` stayed
            ## pinned at 1.000000.
            ##
            ## ⚠⚠ AND `q` ENTERS THE COVARIANCE QUADRATICALLY, so the orbital
            ## covariance came out too large by exactly `C^2`.  Measured
            ## against the independent Lyapunov reference before the fix:
            ## ratio 0.0624 / 1.0001 / 16.043 at those same `C` -- right ONLY
            ## at `C = 1`, which is the only place A9's three-way gate ever
            ## ran.  §D 0c, on the very circuit that produced that entry: a
            ## unit reactance makes `C` the identity and the two inner
            ## products indistinguishable.
            ##
            ## The pair-slicing correction this block was written for is
            ## subsumed: normalising on `q^T C p` fixes the slice scale and
            ## the weighting in one step.
            _x0r = np.delete(np.asarray(self.waveform[1], dtype=float)[:, 0],
                             self.irefnode)
            _Cm = np.asarray(self._C_at(_x0r), dtype=float)
            c0 = complex(np.vdot(q[:, 0], _Cm @ p[:, 0]))
            if abs(c0) < 1e-30:
                raise ValueError(
                    'PSS.floquet_modes: mode %d has q(0)^T C p(0) = %.3e on '
                    'the state block, so it cannot be biorthonormalised '
                    'there.' % (k, abs(c0)))
            q = q / np.conj(c0)
            out.append({'lam': lk, 'mu': muk, 'u0': uk, 'v0': vk,
                        'p': p, 'q': q, 'times': tt, 'c0': c0,
                        'residual': float(np.linalg.norm(M @ uk - lk * uk)
                                          / max(abs(lk), 1e-300))})
        return out

    def _continuous_adjoint(self, fp, lam, mu, times):
        """A mode's adjoint `q(t_j)` on a NON-UNIFORM gear grid, by integrating
        the continuous adjoint equation SEPARATELY (2026-09-20).

        The exact transpose of gear's two-step recursion draws `a1` and `a2`
        from LATER steps, so on a grid whose step changes it is a consistent
        scheme for the adjoint equation only to first order, and no per-node
        rescaling lifts it (four measured, all halving per doubling).  Here
        `C^T dq/dt = G^T q` is integrated backwards with BDF2 whose
        coefficients belong to the REVERSE grid's own step pair --
        `companion_coefficients(h_{n+1}, h_{n+2})`:

            (a0' C_n + G_n)^T q_n = -C_n^T (a1' q_{n+1} + a2' q_{n+2})

        The state is the pair (q_{n+1}, q_{n+2}); the backward map over one
        period is built from 2m basis propagations and the mode's eigenvector
        matched to the forward multiplier `lam` (they agree to O(h^2):
        1.2e-3 / 3.1e-4 / 7.7e-5 at N = 200 / 400 / 800).  MEASURED on the
        asymmetric van der Pol, 3:1 grid, invariant `q^T C p` spread:

            transpose, a0-scaled   2.2e-02  1.1e-02  5.9e-03   halving
            this                   1.8e-03  4.6e-04  1.1e-04   QUARTERING

        and `orbital_correlation`'s R against radau's (exact there) 2.3e-02 /
        7.2e-03 / 2.0e-03, monotone at ~x3.5, onto the floor `p` sets.

        ⚠ TWO ADJOINTS IN THE TREE, ON PURPOSE.  This one is NOT the transpose
        of the discrete map: the PPV, the adjoint noise folds and the sideband
        rows keep the exact transpose (pinned at 1e-15) because their
        identities need it.  This serves the MODE SHAPES only, where the
        continuous adjoint is the object wanted.  On a uniform grid the two
        coincide (`a0' = a0`) and the transpose is used; one-step kinds never
        come here.  Dense 2m x 2m up to `CONTINUOUS_ADJOINT_DENSE_M`
        unknowns; above that, Arnoldi on the backward map with Ritz-residual
        certification (see the note at the branch).  ⚠ Biorthonormality against the OTHER modes is not
        exact here (it is for the transpose) -- each mode is normalised on its
        own `q(0)^T C p(0)`, as before.
        """
        N = len(fp.steps)
        m = self.cir.n - 1
        tms = np.asarray(times, dtype=float)
        h = np.diff(tms[:N + 1])
        W = np.delete(np.asarray(self.waveform[1], dtype=float),
                      self.irefnode, axis=0)
        Cs = [np.asarray(self._C_at(W[:, j]), dtype=float) for j in range(N + 1)]
        Gs = [np.asarray(self._G_at(W[:, j]), dtype=float) for j in range(N + 1)]
        integ = self._integrator_for(self.par.method)
        import scipy.linalg as _sla
        coef, lus = [], []
        for n in range(N):
            a = [float(x) for x in
                 integ.companion_coefficients(h[n % N], h[(n + 1) % N])[0]]
            coef.append(a)
            ## one factorisation per node, reused by every propagation
            lus.append(_sla.lu_factor((a[0] * Cs[n] + Gs[n]).T))

        def propagate(qN, qN1):
            q = [None] * (N + 2)
            q[N], q[N + 1] = qN, qN1
            for n in range(N - 1, -1, -1):
                _a0, a1, a2 = coef[n]
                q[n] = _sla.lu_solve(lus[n], -Cs[n].T @ (a1 * q[n + 1] + a2 * q[n + 2]))
            return q

        if m <= self.CONTINUOUS_ADJOINT_DENSE_M:
            ## small: the dense 2m x 2m map, exact
            Mp = np.zeros((2 * m, 2 * m), dtype=complex)
            for c in range(2 * m):
                e = np.zeros(2 * m, dtype=complex)
                e[c] = 1.0
                qq = propagate(e[:m], e[m:])
                Mp[:, c] = np.concatenate((qq[0], qq[1]))
            lams, vecs = np.linalg.eig(Mp)
            k = int(np.argmin(np.abs(lams - lam)))
            w = vecs[:, k]
        else:
            ## ⚠ MATRIX-FREE ABOVE THAT (2026-09-20, item 3 of "gear as a
            ## first-class choice on non-uniform grids").  The wanted modes --
            ## the phase mode and the slow orbital ones -- are the DOMINANT
            ## eigenvalues of the backward map, so a plain Arnoldi on it with
            ## the Ritz-residual certification `_ritz_second_multiplier` uses
            ## (|h_{k+1,k}| |y_last| / ||y||) reaches them at a basis far below
            ## 2m: MEASURED equal to the dense eigenvector to cos 1.00000000 at
            ## a basis of 3 / 16 / 24 for 2m = 4 / 32 / 124, residuals 1e-16 ..
            ## 1e-52, on the hostile fixture and the ladder oscillator
            ## re-solved on a 3:1 grid.  Grown from `PPV_RITZ_BASIS` toward
            ## `PPV_RITZ_MAX_BASIS` until the matched pair certifies; a pair
            ## that never certifies is refused, not returned.
            n2 = 2 * m
            kmax = int(min(n2, self.PPV_RITZ_MAX_BASIS))
            kk = int(min(n2, self.PPV_RITZ_BASIS))
            _rng = np.random.default_rng(12345)
            q0 = _rng.standard_normal(n2).astype(complex)
            Qb = [q0 / np.linalg.norm(q0)]
            H = np.zeros((kmax + 1, kmax), dtype=complex)
            kdone = 0
            while True:
                while kdone < kk:
                    qq = propagate(Qb[kdone][:m], Qb[kdone][m:])
                    wv = np.concatenate((qq[0], qq[1]))
                    for i in range(kdone + 1):
                        H[i, kdone] = np.vdot(Qb[i], wv)
                        wv = wv - H[i, kdone] * Qb[i]
                    for i in range(kdone + 1):
                        cc = np.vdot(Qb[i], wv)
                        H[i, kdone] += cc
                        wv = wv - cc * Qb[i]
                    H[kdone + 1, kdone] = np.linalg.norm(wv)
                    kdone += 1
                    if H[kdone, kdone - 1] < 1e-14:
                        break
                    Qb.append(wv / H[kdone, kdone - 1])
                th, Y = np.linalg.eig(H[:kdone, :kdone])
                i = int(np.argmin(np.abs(th - lam)))
                resid = (abs(H[kdone, kdone - 1]) * abs(Y[kdone - 1, i])
                         / max(float(np.linalg.norm(Y[:, i])), 1e-300)
                         if kdone < n2 else 0.0)
                if resid <= self.PPV_RITZ_RESIDUAL_TOL or kdone >= kmax:
                    break
                kk = int(min(2 * kk, kmax))
            if resid > self.PPV_RITZ_RESIDUAL_TOL:
                raise ValueError(
                    'PSS.floquet_modes: the continuous adjoint\'s Arnoldi did '
                    'not certify the mode at |lam| = %.6f (Ritz residual %.1e '
                    'at a basis of %d of %d). Use a uniform grid, or '
                    'method=\'radau\'.' % (abs(lam), resid, kdone, n2))
            w = np.zeros(n2, dtype=complex)
            for j in range(kdone):
                w += Y[j, i] * Qb[j]
        qq = propagate(w[:m], w[m:])
        ts2 = tms[:N + 1]
        q = np.column_stack([qq[j] * np.exp(mu * ts2[j]) for j in range(N + 1)])
        return q, ts2

    #: below this many unknowns `_continuous_adjoint` forms the dense 2m x 2m
    #: backward map (exact); above it, Arnoldi on the map (Ritz-certified)
    CONTINUOUS_ADJOINT_DENSE_M = 8

    #: `sampled_noise` warns when its top sideband sits above this many
    #: radians per grid step -- the regime where the covered bands are
    #: discretisation-limited (measured: 3x the pure tail low at 3.1)
    SAMPLED_RESOLUTION_WARN = 1.0

    #: `ppv` warns on an autonomous run whose period map's unit multiplier
    #: is further than this from 1 -- `c` is then uncertain by 5-12x that
    #: departure (measured); see the note in `ppv`
    PPV_UNIT_MULTIPLIER_WARN = 1e-2
