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

    ## STARTING Arnoldi basis size for the second-multiplier estimate:
    ## `_ritz_second_multiplier` is grown from it until the selected pair's
    ## own Ritz residual certifies.  Exact at `k = n`; below that a
    ## truncation, and on a circuit monodromy (not normal) a truncated `lam2`
    ## errs in both directions -- see `ppv`.
    ## History: `doc/shooting_history.md`, `_PPVFloquet.PPV_RITZ_BASIS`.
    PPV_RITZ_BASIS = 12

    ## ⚠ THE GATE ON A TRUNCATED `lam2`: the SELECTED PAIR's Ritz residual,
    ## `|h_{k+1,k}| |y_i[last]| / ||y_i||`, free from the `H` already formed.
    ## It separates a converged pair from a leaked one, which the GMRES
    ## (solve) residual cannot.  Measured, the two populations touch at
    ## ~1e-5, so the robust band is below ~1e-6.
    ## History: `doc/shooting_history.md`, `_PPVFloquet.PPV_RITZ_RESIDUAL_TOL`.
    PPV_RITZ_RESIDUAL_TOL = 1e-6

    ## ⚠ A COST CEILING, NOT A CORRECTNESS THRESHOLD: `k` doubles until the
    ## residual certifies or this is reached, and overrunning it gives a
    ## WARNING and an uncertified number, never a silently wrong one.  `k`
    ## tracks the slow-mode count rather than `n` (measured on a synthetic);
    ## 128 leaves 2x headroom on the largest published case, Lai's
    ## 64-gated-capacitor DCO (`n = 813`).
    ## History: `doc/shooting_history.md`, `_PPVFloquet.PPV_RITZ_MAX_BASIS`.
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

        ⚠ `eig`, NOT `eigvals`: the eigenVECTOR's last component is half the
        residual.

        Every local is underscored, a convention inherited from when this was
        inline in `ppv`.

        `residual` is `inf` when no Ritz value survives the deflation, so a
        caller that gates on it cannot read "nothing found" as "certified".

        History: `doc/shooting_history.md`, `_PPVFloquet._ritz_second_multiplier`.
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

        ⚠⚠ THERE ARE TWO ADJOINT VECTORS; DO NOT CONFLATE THEM.  Demir 2000
        puts both conventions on one page: a STATE initial condition
        contracts as `v_1^T(0) C(0) x(0)` (eq 41, WITH `C`), while an
        EQUATION-ROW input contracts as `v_1^T(s) b(s)` (eq 42, and the phase
        equation 44, BARE).  `ppv()` returns `C^T v_1`, which is the right
        object for a state perturbation.  `CY` is an equation-row covariance
        -- a current injected into a KCL row -- so `diffusion_constant` and
        `colour_projection` need `v_1`, and this produces it.  Contracted with
        `C^T v_1` instead, `diffusion_constant` is wrong by `C^2` on
        differential states and annihilated on algebraic ones -- invisible
        on any fixture with `C = 1 F`.

        ⚠⚠ AND `C` IS NEVER INVERTED.  The solve is on `C[D, NZ]^T`, the
        block between the DIFFERENTIAL equations and the NON-algebraic
        states, which is square and invertible by construction; the
        singular `C` as a whole is not touched.  The algebraic entries come
        from the constraint (`_algebraic_adjoint_fill`) instead.

        History: `doc/shooting_history.md`, `_PPVFloquet._equation_row_ppv`.
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
        dynamics through the CONSTRAINT rather than through its own row (on
        a tank with series loss, eliminating `v_x = r (i_L + b)` puts `-r b`
        into the inductor's equation).  Left at zero, noise landing there is
        contracted against a structural zero.

        The adjoint equation's ALGEBRAIC-STATE columns are what determines
        them.  For a column `j` with no `C` entry the equation carries no
        derivative, so it reads `sum_i G_ij v_i = 0`, and splitting `i`
        into algebraic and differential rows gives

            v_A  =  (G[A, Z]^T)^-1 G[D, Z]^T v_D

        ⚠ THE MAGNITUDE IS STRUCTURAL AND THE SIGN IS MEASURED.  On a
        resistive divider the algebraic entries stand in the ratio the
        topology fixes; the sign follows from requiring algebraic and
        differential rows to share ONE convention, `dT/dA = +integral v_j`,
        measured on that divider with a DC-injection probe.  ⚠ A SINGLE
        SERIES RESISTOR CANNOT SETTLE THE SIGN: both signs fit there.

        Returns `vblock` unchanged, with a warning, when the structure is
        not index-1: `len(rows) != len(cols)` or a singular block.  That
        case is section B4's, and guessing at it would be worse than
        leaving a known zero.

        History: `doc/shooting_history.md`, `_PPVFloquet._algebraic_adjoint_fill`.
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
        ## ⚠ THE DERIVED SIGN, because this acts on `v_1`: applied to
        ## `C^T v_1` it comes out negated wherever the term is dominated by an
        ## inductor branch row (which carries `-L` in `C`).  Pinned by the
        ## eq (24) constraint.
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
        ## complex for a Floquet mode's adjoint (`_floquet_mode`)
        out = np.array(vblock, dtype=np.result_type(np.asarray(vblock), float),
                       copy=True)
        out[rows] = va
        return out

    def _ppv_propagate(self, fp, v, m, xdot, alg_rows, alg_cols, inject=None):
        """The pair-consistent SECOND-ORDER propagation of an anchor vector
        `v` over the period -- the block `ppv()` applies to its null vector,
        shared with `frequency_aware_ppv`, which runs it on a COMPLEX anchor.
        Returns `(states, states_pair, ts, Xf)`: `states` the per-step
        state-space samples (`C^T v`, second order, rescaled by
        `v . xdot = 1` at the first sample), `states_pair` the raw
        pair-space replay, `ts` the transposed per-step states, `Xf` the
        orbit.  Only the `solved_history` (LMM) kind carries the
        correction; the others return the replay as it is.

        History: `doc/shooting_history.md`, `_PPVFloquet._ppv_propagate`.
        """
        _dt = complex if np.iscomplexobj(v) else float
        _alg_rows, _alg_cols = alg_rows, alg_cols
        self._ppv_alg_fallback = False
        _end, _ts, states = fp.matvec_transposed(v, collect=True, inject=inject)
        states_pair = [np.array(st, dtype=_dt, copy=True) for st in states]
        _Xf = np.asarray(self.waveform[1], dtype=float)
        ## ⚠⚠ THE PAIR'S FIRST BLOCK IS NOT THE PPV: its error is first order
        ## and grows with Q.  For Gear-2 the adjoint state is the pair
        ## `(w1, w2) = (dphi/dx_k, dphi/dx_{k-1})`, and `w1` alone is the
        ## response to a perturbation of `x_k` WITH `x_{k-1}` HELD -- an
        ## inconsistent history.  A physical state perturbation moves both:
        ## `dx_{k-1} = Phi(t_{k-1}, t_k) dx_k`, so the phase functional is
        ##
        ##     v(t_k) = w1 + Phi(t_{k-1}, t_k)^T w2,   Phi ~ I - h J + O(h^2)
        ##
        ## and with `w2 = C_{k-1}^T z` (`z = -a2 t_k`, exact by the
        ## recursion) that is `w1 + (C_{k-1} + h G)^T z` -- no inverse of
        ## `C`, so it holds for a DAE.  This is second order and holds
        ## `v(t) . xdot(t) = 1` along the orbit.  ⚠ Van der Pol cannot tell
        ## the two apart (its rows are in quadrature); a non-isochronous core
        ## such as `vdp + 0.3 u^2` can.
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
                ## ⚠ DIFFERENTIAL ROWS ONLY.  `w2 = C^T z` does not see `z`'s
                ## algebraic rows, where the decomposition is non-unique and
                ## the multipliers are O(1/h), so `h G^T z` would carry an
                ## O(1) component along the constraint normal.  The
                ## consistent propagation `C_D dx_{k-1} = (C_D + h G_D) dx_k`
                ## involves only the differential equations and is unique.
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
                    ## complement `G[D,NZ] - G[D,Z] G[A,Z]^-1 G[A,NZ]`;
                    ## the full `G` instead is first order.
                    _A = np.asarray(_alg_rows, dtype=int)
                    _Zc = np.asarray(_alg_cols, dtype=int)
                    _D = np.array([i for i in range(m) if i not in _alg_rows],
                                  dtype=int)
                    _NZ = np.array([j for j in range(m) if j not in _alg_cols],
                                   dtype=int)
                    ## ⚠ AND `G[A,Z]` NONSINGULAR *IS* THE INDEX-1 CONDITION.
                    ## At index >= 2 (an L-I cutset, a C-V loop) it is
                    ## singular by definition and the complement does not
                    ## exist (`li_plus_rc`, `cv_plus_rc`): warn once with the
                    ## reason and fall back to the full-`G` propagation (its
                    ## order: see the warning, and the non-uniform-grid note
                    ## in `ppv`).
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
                ## slaved, and `h G^T z` would otherwise leave a residue
                ## there (the Demir-(24) gate).  The equation-row conversion
                ## never reads these entries.
                if _alg_cols:
                    _vp_j[np.asarray(_alg_cols, dtype=int)] = 0.0
                _vphys.append(_vp_j)
            _scale = (_vphys[0] @ xdot) if _dt is complex else float(_vphys[0] @ xdot)
            if _scale == 0.0:
                raise ValueError(
                    'PSS.ppv: the pair-consistent adjoint is orthogonal to '
                    'the orbit tangent at t = 0.')
            ## ⚠ AND ITS DC CONTENT IS THE CONSISTENT OBJECT'S TOO: the mean
            ## taken from the raw block is first order and 8 % off on an
            ## inductor row (exact on some other rows, mechanism open), and
            ## stitching it in breaks `v . xdot = 1`.  So `samples` is one
            ## object, second order everywhere, with a ~1e-5 |v| absolute
            ## floor on its mean.  The raw pair is kept as `samples_pair` for
            ## the structural gates that live on its discrete identities.
            states = [np.concatenate((vp / _scale, st[m:] / _scale))
                      for vp, st in zip(_vphys, states)]
        return states, states_pair, _ts, _Xf

    def _state_map(self):
        """The period map on the STATE that every state-space consumer reads
        -- `ppv`, `floquet_modes`, `PAC.solve`, the adjoint rows, the
        deflated solve, `sampled_noise`: the factored period, except that a
        Nordsieck GLM's is taken on the state (`_GLMPeriod.state_map`) --
        its own acts on the Nordsieck vector, whose null vector's first
        block holds the higher components fixed.  ⚠ Until 2026-09-24 a
        native GLM's `ppv()` returned that Nordsieck object, `r*m` wide; until
        2026-09-25 PAC and the noise folds read a GLM run from a radau twin
        (named `_ppv_map` then, `ppv`'s alone)."""
        fp = self.factored_period()
        return fp.state_map() if fp.is_glm else fp

    def ppv(self, tol=None):
        """The perturbation projection vector at `t = 0` (Demir & Roychowdhury).

        Returns `(v, info)`.  `v` is the pair-space left null vector of
        `I - M`, normalised so that `v . xdot(0) = 1` -- see the note on the
        normalisation in the body.  The phase shift caused by a state
        perturbation `delta` at `t = 0` is then `v[:m] . delta`.  `info`
        carries both border residuals, the null residual, `q` and the scaled
        tangent.  `c` -- the diffusion constant this vector feeds -- has the
        designer-facing reading "JITTER PER SECOND".

        ⚠ AN AUGMENTED SOLVE, NOT AN EIGENVECTOR.  On a high-Q oscillator the
        monodromy has several multipliers numerically indistinguishable from
        the unit one, and SELECTING the unit eigenvector is infeasible (Demir
        & Sangiovanni-Vincentelli 1998, Table 6.4); selecting it by its inner
        product against `C(0) xdot(0)` (Demir, IJCTA 2000) risks "a potential
        breakdown" (Demir 2003).  The single linear solve (Demir, Long &
        Roychowdhury, ICCAD 2000; the fuller procedure in Demir 2003) makes
        that vector the BORDER instead, so the candidate is unique and
        nothing is selected:

            [ I - M^T   q ] [v]     [0]
            [   q^T     0 ] [y]  =  [1]

        A bordered solve never has to tell the other multipliers apart:
        near-degenerate multipliers make EIGENANALYSIS ILL-POSED but this
        solve merely ILL-CONDITIONED (measured and warned about below).
        Neither is the PHASE EQUATION's own limit (slow nodes, below), which
        is not an extraction question at all.  So the second-multiplier
        warning means "you are in the regime this method was invented to
        escape", not "this result is degrading".

        ⚠ `y` COMING BACK ZERO IS A FREE CORRECTNESS CHECK.  With a zero first
        block on the right-hand side, `(I - M^T) v + y q = 0` forces
        `y q = 0`, so a nonzero `y` means the border absorbed a residual the
        null space should have taken -- the computed `v` is not in the null
        space.

        ⚠ `q` IS EXACT, NOT DIFFERENCED.  `q = C(0) xdot(0)`, and the DAE
        gives it directly: `dq/dt + i(x) + u(t) = 0` and `dq/dt = C xdot`, so
        `q = -(i(x_0) + u(0))` -- two evaluations at the converged solution,
        no derivative anywhere.

        ⚠ NAMES.  The PPV is Kaertner's adjoint LPTV impulse response, the PRC
        of mathematical biology and Hajimiri's NUMERICAL ISF.  His
        CLOSED-FORM ISF is NOT: it is the normalised tangent and can get the
        SIGN of a phase change wrong.  Nor is `xdot` a substitute (the two
        scale oppositely with the RC time constant): here `xdot` is only the
        normalisation and, through `q`, the border.

        ⚠ THE PPV IS ONE RUNG ON A LADDER (Suvak & Demir, TCAD 2011): the
        LINEAR-isochron approximation of an exact but unusable phase
        equation.  A QUADRATIC rung adds isochron curvature and would earn
        its cost on LARGE perturbations (injection locking, big interferers),
        not on phase noise -- and it would NOT fix the slow-node boundary
        below, which is in the DYNAMICS, not the amplitude.

        ⚠ ITS VALIDITY BOUNDARY IS SLOW NODES, AND THE VAN DER POL GATE
        CANNOT SEE IT.  The phase equation `alpha' = v_1^T(t+alpha) b(t)`
        treats the oscillator's frequency response as INSTANTANEOUS, so a
        slow node FILTERS the noise of devices near it, the PPV cannot see
        the filtering, and phase noise is OVER-ESTIMATED; better extraction
        does not help (Lai 2008).
        `test_the_ppv_predicts_a_phase_shift_the_oscillator_actually_has`
        (van der Pol, perturbed at its core) passes whether or not this
        failure mode is present: it establishes the extraction and
        normalisation, not the model's range.  `c` from this PPV is right;
        the spectral SHAPE above `T/(2 pi tau)` is what the frequency-aware
        PPV (`frequency_aware_ppv`, this bordered system at nonzero `w_s`)
        corrects -- see `oscillator_spectrum`.  A Monte Carlo of `c` cannot
        see the effect, and `_vdp_with_slow_node` tests conditioning only
        (its slow node barely reaches the phase).

        ⚠ ON A STAGED SOLVE (`state_events=True`) THE PPV IS BORDERED: the
        null vector is that of the TOTAL monodromy `M + P_theta dtheta/dx_0`,
        and the samples carry the crossings' motion as costate injections
        `-zeta_k W_k` at the event nodes, `zeta = Gt^-T P_theta^T v` -- one
        reverse pass, and it IS the phase gradient at fixed time.  Verified
        against the exact piecewise-linear PPV of a comparator relaxation
        oscillator (every sample within 0.4 %).  ⚠ A transient-FD check must
        not take its phase reference level from the record it measures:
        read the phase at the comparator's own crossing.

        History: `doc/shooting_history.md`, `_PPVFloquet.ppv`.
        """
        _tw = self.monodromy_twin()
        if _tw is not self:
            return _tw.ppv(tol)
        import scipy.sparse.linalg as spla
        fp = self._state_map()
        ## Every call below goes through `fp.matvec_transposed`/`fp.matvec`,
        ## so the map's kind is the dispatcher's business rather than this
        ## method's (B8).
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

        ## ⚠ ON A STAGED SOLVE THE MONODROMY IS THE TOTAL ONE: a perturbation
        ## moves the crossings, `M_tot = M + P_theta dtheta/dx_0`, and the PPV
        ## is its left null vector.  Its samples along the orbit carry the
        ## same correction as costate injections at the event nodes --
        ## `-zeta_k w_k`, `zeta = Gt^-T P_theta^T v` -- which the reverse pass
        ## carries to every earlier node: the saltation matrix's transpose.
        ## A comparator oscillator's PPV jumps at its switching instants, and
        ## this is where the jump comes from.
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
        ## GIVES: on `_vdp_ppv(400)` the two differ by a factor 2.07 AND a
        ## sign.  Displacing the state ALONG the orbit by `eps xdot` advances
        ## the phase by `eps`, so `v . xdot = 1` is the normalisation a state
        ## perturbation sees.  Demir's Remark 3.1 reads `v_1^T C u_1 = 1`; the
        ## vector this solve returns is `C(0)^T v_1` -- the left eigenvector of
        ## the STATE-SPACE monodromy for ANY `C`, since the conserved pairing
        ## gives `M^T (C(0)^T v_1) = C(0)^T v_1` -- so the two statements agree
        ## about different objects.  ⚠ Predicting a state jump's phase shift
        ## as `v^T C delta` is wrong (it does not converge); `v . delta` is
        ## right.
        ## ⚠ THE ALGEBRAIC ENTRIES ARE FILLED *AFTER* THIS, ON PURPOSE: `v` is
        ## also the REPLAY'S SEED, and the algebraic components are SLAVED --
        ## propagating them through the step map corrupts the differential
        ## ones.  Nor does the normalisation see them: a state perturbation of
        ## a DAE lies ON the constraint manifold, while the algebraic entries
        ## answer the equation-row question (a noise current injected into an
        ## algebraic KCL row).
        ## ⚠ THE UNIT MULTIPLIER OFF THE CIRCLE IS A SILENT ERROR IN `c`.  This
        ## solves the bordered system AT lambda = 1 whatever the discrete
        ## period map's own unit multiplier is, so nothing here can see that
        ## multiplier leave the circle (gear's does on a coarse `lte_grid`),
        ## and the solve still reports `converged = True`.  The solve already
        ## records `spectral_radius`; an autonomous run whose unit multiplier
        ## is more than `PPV_UNIT_MULTIPLIER_WARN` off the circle is told so
        ## here, once, with the size and the remedy.
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
        ## ⚠ THE PPV OVER THE PERIOD, not just at `t = 0`: Demir's diffusion
        ## constant `c = (1/T) integral v_1^T(t) B(t) B^T(t) v_1(t) dt` is an
        ## integral over the orbit, and the reverse replay computes
        ## `Phi(T,s)^T v(T) = v(s)` on its way to the answer.
        _inject = self._event_costate_injection(fp, v, n)
        states, states_pair, _ts, _Xf = self._ppv_propagate(fp, v, m, xdot, _alg_rows, _alg_cols,
                                                          inject=_inject)
        ## ⚠ INDEX >= 2 ON A NON-UNIFORM SOLVED-HISTORY GRID: the fallback in
        ## `_ppv_propagate` keeps only the differential block of the
        ## pair-consistent correction, and the coupling it drops -- a
        ## DERIVATIVE term at index 2, with no Schur complement to stand in
        ## for it -- cancels between steps on a uniform grid and not on a
        ## non-uniform one, so the samples would be first order.  On that grid
        ## they come from the continuous adjoint's phase mode instead (the
        ## object `floquet_modes` uses there), as the STATE-SPACE PPV
        ## `v_j = C_j^T q_j` -- the continuous `q` IS Demir's equation-row
        ## `v_1`, and `samples` carry `C^T v_1` (see `_equation_row_ppv`), so
        ## `q` put there bare is wrong -- scaled so `q_0^T C_0 xdot_0 = 1`,
        ## which is `ppv()`'s own `v . xdot = 1`.  A quadrature of them needs
        ## a periodic trapezoid (PAC's `_period_weights`): a left-rectangle
        ## rule on a non-uniform grid is first order by itself.  The anchor
        ## `v` and the pair block are untouched; a uniform grid never comes
        ## here, and neither do the one-step kinds.
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
        ## nothing to propagate.  The pair's SECOND block is the history term
        ## and is deliberately not filled -- `v(t)` is the first block.
        ## ⚠⚠ THE EQUATION-ROW ADJOINT IS A SECOND OBJECT, NOT A CORRECTION
        ## TO THE FIRST.  `states` and `v` stay `C^T v_1`, the vector a STATE
        ## perturbation contracts with, which is what this method documents
        ## and what every existing gate measures; converting one into the
        ## other would silently change what `ppv()` returns.
        ## ⚠ ONE WARNING PER CALL, NOT ONE PER SAMPLE: at index 2 the fill
        ## warns at every sample.
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
        ## Measured, `sigma_min(bordered)` tracks `T/tau` over six decades
        ## while the null residual does not move: GMRES converges, the answer
        ## looks clean, and the conditioning has lost six digits.  So this
        ## estimates `|lambda_2|` explicitly rather than trusting a small
        ## residual.  Its accuracy cost in `c` is not established: a Monte
        ## Carlo at `tau/T = 10` read 0.80 of the prediction (~2.5 sigma, in
        ## the direction opposite to the reported effect), and larger `tau/T`
        ## (Lai's case is ~1e9) is untested because the measurement needs ~15
        ## time constants of settling.
        ##
        ## `lam2` is exact on the dense route below.  Above
        ## `FLOQUET_DENSE_LIMIT` it comes from a truncated Arnoldi on `I - M`
        ## (`_ritz_second_multiplier`), and a truncated `lam2` is a lower
        ## bound only for a NORMAL `M` (Cauchy interlacing).  ⚠ A circuit
        ## monodromy is not normal: with several near-unit clusters (slow
        ## modes) the error is NOT one-signed -- `Q` 19x low, 4x high, or
        ## `lam2 > 1` (a spurious unstable multiplier) -- and not Q-specific.
        ## ⚠ Raising `PPV_RITZ_BASIS` is not the fix: the basis has to grow
        ## with the slow-mode count, so it doubles until the selected pair's
        ## Ritz residual certifies (`PPV_RITZ_RESIDUAL_TOL`), up to
        ## `PPV_RITZ_MAX_BASIS`, and `info` reports whether it did.
        ## ⚠ The fixed `|lam - 1| > 1e-6` unit-mode filter is the weak point
        ## at eigenvector conditioning `cond(V) >= 1e4`: the phase mode
        ## survives the discard and is selected as `lam2` (`Q -> inf`);
        ## deflating it explicitly with `q` would sidestep that.
        ## Sorted by real part, not magnitude, because an amplitude mode is
        ## real and positive while a complex pair of larger modulus would be
        ## an oscillation about the orbit.
        vu = float(v[:m] @ u[:m] + v[m:] @ u[m:])
        lam2 = 0.0
        ## `_resid`/`_certified` describe how `lam2` was obtained; the
        ## degenerate `n < 2` fall-through never enters either branch, and
        ## `lam2 = 0` there is exact rather than estimated.
        _resid, _certified = 0.0, True
        kk = int(min(n, self.PPV_RITZ_BASIS))
        ## ⚠⚠ DENSE WHENEVER IT IS AFFORDABLE.  Forming `M` by `n` matvecs and
        ## taking its exact spectrum has no threshold, no basis size and no
        ## selection ambiguity; the truncated Arnoldi below has all three.
        ## `FLOQUET_DENSE_LIMIT` is the same cap `floquet_modes` applies to
        ## the same assembly, so the two agree about what "affordable" means.
        ## Stage maps keep the dense route ABOVE it as well: there it is
        ## expensive, but the alternative is not slower, it is WRONG.
        _dense_ok = (fp.is_stage
                     or n <= self.FLOQUET_DENSE_LIMIT)
        if _dense_ok:
            ## ⚠ THE STAGE MAP IS DENSE AND WIDTH `m`, so its exact spectrum
            ## is cheap -- and the Arnoldi below resolves it BADLY: `I - M`
            ## has `M`'s annihilated modes clustered at 1 along with the
            ## physical unit root, which an Arnoldi leaves past the `1e-6`
            ## deflation, reporting the ORBIT TANGENT as the second
            ## multiplier.  The exact eigenvalues deflate the unit root to
            ## machine precision and give the true second multiplier.
            _Md = np.column_stack([np.asarray(fp.matvec(_e), dtype=float)
                                   for _e in np.eye(n)])
            _lams = np.linalg.eigvals(_Md)
            _keep = np.real(_lams)[np.abs(_lams - 1.0) > 1e-6]
            if _keep.size == _lams.size:
                ## ⚠ NO MULTIPLIER IN THE WINDOW: for a MULTISTEP or
                ## trapezoidal solve the phase multiplier is 1 only as far as
                ## the discretisation is time-translation invariant, and a
                ## non-uniform `grid=` breaks that at O(h^2) (radau's
                ## collocation keeps it at 1 to ~1e-11).  Without this the
                ## phase multiplier itself is reported as the SECOND one,
                ## silently.  Drop the one nearest 1 instead; a uniform grid
                ## never gets here.
                _keep = np.real(np.delete(_lams, int(np.argmin(np.abs(_lams - 1.0)))))
            if _keep.size:
                lam2 = float(max(np.max(_keep), 0.0))
            ## the spectrum is exact, so there is nothing to certify against
            _resid, _certified = 0.0, True
        elif kk >= 2:
            ## ⚠⚠ THE TRUNCATED PATH, GATED ON THE PAIR'S OWN RITZ RESIDUAL
            ## AND GROWN UNTIL IT CERTIFIES.  It runs only where the dense
            ## spectrum is unaffordable -- exactly where the truncation is
            ## least trustworthy, since a big circuit is the one likely to
            ## carry the many slow nodes that break the selection.  `k` has to
            ## track the slow-mode count, so no fixed basis works; doubling
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
        ## `Q = log(threshold)/log|lambda_2|` (Wang & Roychowdhury): an
        ## amplitude perturbation decays to `|lambda_2|` of its size each
        ## cycle, so the cycles needed to fall below a threshold IS the
        ## oscillator's Q.  Reported for a `1/e` threshold, so `Q` is
        ## cycles-to-1/e.  (`f_r/df` and stored/dissipated do not apply to an
        ## autonomous circuit, and this is not the resonator's Q.)  "High Q",
        ## "a second multiplier near 1", "slow amplitude restoration" and "a
        ## long time constant" are one condition -- the one behind every
        ## failure this class warns about.
        ##
        ## ⚠⚠ ITS NAME IS ONLY RIGHT WHILE THE OSCILLATOR'S AMPLITUDE MODE IS
        ## THE SLOWEST NON-UNIT MODE.  A parasitic with `tau_p/T > Q_osc` IS
        ## the second multiplier, and `Q` then reports THE PARASITIC'S DECAY
        ## TIME IN PERIODS (a DCO's gated capacitors sit permanently there).
        ## Read `Q` as "cycles for the SLOWEST NON-UNIT MODE to decay by 1/e".
        ##
        ## ⚠ `Q` AMPLIFIES THE ERROR IN `lambda_2`:
        ##
        ##     (dQ/Q) / (dlambda_2/lambda_2)  =  -1/ln(lambda_2)  =  Q
        ##
        ## so the resolution needed for a given accuracy in `Q` scales with
        ## `Q`, and a method that BIASES `lambda_2` at fixed order has no
        ## escape from it.
        Q = (-1.0 / np.log(lam2) if 0.0 < lam2 < 1.0 else float('inf'))
        info = {'border_residual': y,
                'tangent_border_residual': yf,
                'Q': Q,
                'null_residual': resid / max(float(np.linalg.norm(v)), 1e-300),
                ## ⚠ MULTIPLY `null_residual` BY THIS TO GET THE RELATIVE
                ## ERROR IN `v` THE RESIDUAL CANNOT EXCLUDE.  `null_residual`
                ## is `||v - M^T v|| / ||v||`, so an error component along the
                ## `lam2` left-eigendirection enters it scaled by `1 - lam2`
                ## and is nearly INVISIBLE exactly when `lam2 -> 1` (a random
                ## error is still caught).  A FLAT `null_residual` IS NOT
                ## EVIDENCE OF ACCURACY: read the two numbers together or
                ## neither.
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
                ## O(h^2)).  A quadrature over `times` divides by THIS.
                'period': float(fp.T)}
        return v, info

    def frequency_aware_ppv(self, offset, tol=None):
        """The PPV at a nonzero modulation frequency (Lai 2008, eq. 23).

        The classical PPV is the left null vector of `I - M^T`, bordered by
        `q = C(0) xdot(0)`; it is the phase response to a perturbation that
        is SLOW against every other Floquet mode.  This is the SAME
        bordered system at `alpha = exp(-j w_s T)`,

            [[I - alpha M^T,  q], [q^T, 0]] [v; y] = [0; 1],

        whose solution `v(w_s)` is the phase sensitivity to a perturbation
        modulated at `w_s`: at `w_s = 0` it IS `ppv()` (pinned), and away
        from it the AMPLITUDE mode admixes with weight
        `(1 - alpha)/(1 - alpha mu_2)` -- zero at DC, rising ten-fold per
        decade, cornering where `2 pi f_s T = 1 - mu_2` and flat above.
        Eq. 23, not 24: `I - exp(-j w_s T) M^T` is the exact sampled LPTV
        adjoint at ANY offset (the monodromy carries the full time
        variation), where Lai's eq. 24 drops the AC columns of the Toeplitz
        block and holds only near DC; eq. (24) at `w_s = 0` "is the augmented
        PPV extraction equation".  His construction is harmonic balance,
        this is the shooting basis.

        Returns `(v, info)`: `v` the pair-space anchor vector (complex);
        `info['samples_pair']` the T-periodic envelope over the period,
        `lambda_k exp(+j w_s t_k)` with `lambda` the transposed replay of
        `v` -- the sideband rows' own convention, so its Fourier
        coefficient at harmonic `k` is the phase transfer of a source band
        at `k f0 + f_s`; `info['samples']` the same in state space
        (`C^T v`), SECOND order in the step through `ppv()`'s own
        pair-consistent propagation run on the complex anchor
        (`_ppv_propagate`; at w_s = 0 it is `ppv()`'s `samples`);
        `info['admixture']` the norm fraction of `v(w_s)` orthogonal to the
        DC PPV; `info['corner']` the predicted corner `|1 - mu_2| / (2 pi T)`
        in Hz; `info['alpha']`; `info['ppv']` the DC object's info.

        ⚠ WHAT IT IS FOR.  A source that reaches the phase through a slow
        path (an RC leg, tau >> T) is filtered at its own corner, and the
        DC PPV cannot see that; the harmonic sum built from THIS object's
        coefficients carries the filter inside `c_k(w_s)` with no explicit
        model of the path (A2, roadmap).  The ratio
        `sum_k |c_k(w_s)|^2 / sum_k |c_k(0)|^2` reproduces `pnoise`'s
        `S_pm/(4 S_v)` for a source behind a slow node to within 2 % up to
        r = 5e-2; at r = 0.1 the two differ by -6 %, a gap between PM by
        sideband quadrature and phase-mode projection.

        ⚠ COMPARE LIKE WITH LIKE.  This object and `S_pm` are PM by
        quadrature of the FUNDAMENTAL'S sidebands; compare them with that,
        not with a demodulated or zero-crossing phase -- both of those
        estimators read a given sideband PM with a fixture-dependent gain
        (they pick up the other harmonics' sidebands).  And compare band
        with band: behind a slow node the in-band spectrum is not 1/r^2.
        The slow multiplier's own coefficient is `mode_content[0]` (its
        plateau scales as T/tau).

        ⚠ DO NOT GATE ON `|v|`: with `q^T v = 1` the `1/(1 - alpha)` pole
        cancels between numerator and denominator and the norm is
        frequency-flat by construction.  The change is a DIRECTION -- read
        `admixture`, or the per-harmonic coefficients.

        History: `doc/shooting_history.md`, `_PPVFloquet.frequency_aware_ppv`.
        """
        import scipy.sparse.linalg as spla
        v0, info0 = self.ppv(tol)
        fp = self._state_map()
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
        ## the map is linear and commutes with the per-step phase factor.
        ## ⚠ `times` above is truncated to the sample count, ONE ENTRY SHORT
        ## of the orbit's grid -- a quadrature over it drops the last step;
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

    ## the widest monodromy `floquet_modes` assembles densely (every mode --
    ## what the modal spectra need); above it only the DOMINANT modes, from a
    ## Ritz-certified Arnoldi (`_floquet_modes_ritz`), and a spectrum is
    ## `pnoise`'s.  A requirement of the spectra, not a budget: their weight
    ## is spread over m/n = 0.97 of the modes.
    FLOQUET_DENSE_LIMIT = 400
    ## below this a multiplier is an annihilated algebraic
    ## direction, not a mode -- see `floquet_modes`
    FLOQUET_NULL_TOL = 1e-12

    def floquet_modes(self, pss_unused=None, nmodes=None, fp=None):
        """⚠ THE MODES' ACCURACY IS THE METHOD'S.  Gear's adjoint modes are
        second order on a uniform grid.  On a non-uniform grid the transpose
        of a variable-step multistep method is first order and no rescaling
        makes it more, so gear's modes there come from `_continuous_adjoint`;
        gear and trap also leave the phase multiplier O(h^2) off 1 there and
        the modal spectra refuse.  Radau's modes on a 3:1 grid reproduce the
        uniform grid's to 1e-10, its phase multiplier at 1 + 1e-11: for the
        modes on a non-uniform grid use radau (the default).  Gate on the
        invariant `q^T C p` and the spectra, which are phase-insensitive; a
        mode's Fourier coefficients also carry its phase across N.

        The Floquet pairs `(λ_l, μ_l, p_l(t), q_l(t))` — A9's prerequisite.

        Returns a list of dicts, one per mode, ordered by `|λ|` descending.
        `nmodes=None` returns EVERY non-null mode, and that default is the
        requirement rather than a convenience:

        ⚠⚠ ALL OF THEM ARE REQUIRED, BY THE SOURCE (Traversa & Bonani, IET
        CDS 2011): orbital fluctuations and the phase-orbital correlation
        need "ALL the direct and adjoint Floquet eigenvectors associated with
        the noiseless limit cycle."  (Cited, not verified here.)

        ⚠ TRUNCATION IS LEGITIMATE ONLY WITH A BOUND (Traversa & Bonani, TCAD
        2013: the error tends to zero with the ratio between the norms of the
        NEGLECTED AND RETAINED ROWS).  A caller who passes `nmodes` owes that
        ratio as the gate; this routine does not compute it.

            lam    the Floquet MULTIPLIER, eigenvalue of the monodromy
            mu     the Floquet EXPONENT, `log(λ)/T` (complex)
            u0,v0  right and left eigenvectors at `t = 0`, biorthonormal
                   (`v_k† u_l = δ_kl`)
            p      `p_l(t_j) = Φ(t_j,0) u_l(0) · exp(−μ_l t_j)` — the
                   T-PERIODIC part, sampled on the PSS grid
            q      the adjoint counterpart from the reverse replay
            times  the grid `p` and `q` are sampled on

        WHY THIS EXISTS: the orbital spectrum `S_yy` (Traversa & Bonani,
        TCAS-I 2011, Lemma 3.5) is a sum of Lorentzians centred at
        `jω₀ + Im{μ_l}` with half-width `|Re{μ_l}| + ½h²ω₀²c`, weighted by
        `C_lhj` (their eq 22), which is built from the FOURIER COEFFICIENTS
        of `u_l(t)` and of `v_l(t)ᵀ B(t)`.  The exponents alone (`|λ₂|`
        included) do not order the result: the eigenvectors can make a mode
        whose exponent is far from zero dominate.

        ⚠ THE PERIODIC PART IS THE OUTPUT, NOT `Φ(t,0)u(0)`: Floquet's
        solution is `p_l(t)exp(μ_l t)` with `p_l` T-periodic, and eq (22)
        wants `p_l`'s Fourier series.  `p_l(T) = p_l(0)` is the gate that
        needs no reference.

        ⚠ DENSE UP TO `FLOQUET_DENSE_LIMIT` (`n` matvecs, then `eig`: every
        mode); ABOVE IT, THE DOMINANT `nmodes` ONLY (`_floquet_modes_ritz`,
        2026-09-25): a Ritz-certified Arnoldi on the map and on its
        transpose, run on `M` itself -- whose outer spectrum, the slow modes,
        converges first -- not on `I - M`, where the physical `lam_2 -> 1`
        is the smallest `theta` and resolves LAST (Garcia, Romero & Acha
        2022).  Each mode dict then carries `certified` and `ritz_residual`;
        an uncertified mode is warned and flagged, never silent.
        `nmodes=None` is refused there, and so are the modal spectra built
        on it: eq (22) needs ALL the modes, whose weight this repo measured
        spread over m/n = 0.97 of them -- above the limit a spectrum is
        `pnoise`'s.

        History: `doc/shooting_history.md`, `_PPVFloquet.floquet_modes`.
        """
        _tw = self.monodromy_twin()
        if _tw is not self:
            ## ⚠ the twin gets `None`, not `pss_unused`: it must read ITS OWN
            ## factored period, which is the whole reason a twin exists.  An
            ## explicit `fp` from the caller still wins.
            return _tw.floquet_modes(None, nmodes, fp)
        ## `pss_unused` IS IGNORED, AS ITS NAME SAYS (the modes are read from
        ## `self`); the parameter stays in the signature because callers pass
        ## it positionally.
        fp = self._state_map() if fp is None else fp
        if getattr(fp, 'is_glm', False) and hasattr(fp, 'state_map'):
            ## a GLM's Nordsieck map handed in: its modes are read on the
            ## state (the Nordsieck eigenvectors are not state-space modes)
            fp = fp.state_map()
        n = fp.width
        T = float(fp.T)
        if n > self.FLOQUET_DENSE_LIMIT:
            return self._floquet_modes_ritz(fp, n, T, nmodes)

        M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                             for e in np.eye(n)])
        ## on a staged solve the map is the TOTAL one: the crossings move
        ## with the state, `M + P_theta dtheta/dx_0`
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
        ## back as noise.  Returning them invites a caller to average over a
        ## mode that means nothing.
        keep = [k for k in range(n) if abs(lam[k]) > self.FLOQUET_NULL_TOL]
        if nmodes is not None and int(nmodes) < len(keep):
            ## ⚠⚠ TRUNCATING BY MULTIPLIER MAGNITUDE IS REFUTED BY THE SOURCE'S
            ## OWN WORKED EXAMPLE (Traversa & Bonani TCAS-I 2011 Sec. V: on
            ## their Colpitts the contribution ordering INVERTS across six
            ## orders in mu).  The caller who truncates owes the dropped
            ## weight as a gate; this says so at the call rather than only in
            ## the docstring.
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
        ## `None` means ALL non-null modes -- the default.
        for k in (keep if nmodes is None else keep[:int(nmodes)]):
            out.append(self._floquet_mode(
                fp, k, complex(lam[k]), U[:, k].astype(complex),
                V[:, k].astype(complex), T, times, n,
                lambda uk, lk: float(np.linalg.norm(M @ uk - lk * uk)
                                     / max(abs(lk), 1e-300))))
        return out

    def _floquet_modes_ritz(self, fp, n, T, nmodes):
        """`floquet_modes` above `FLOQUET_DENSE_LIMIT`: the DOMINANT `nmodes`
        modes (by `|lambda|`), from a Ritz-certified Arnoldi on the map for
        the right vectors and on its transpose for the left ones, paired by
        value as the dense path pairs them; each then goes through the same
        `_floquet_mode`.  On a staged solve both are the TOTAL map's.

        ⚠ DOMINANT MODES ONLY, BY REQUIREMENT, NOT BY BUDGET.  The modal
        spectra need EVERY mode -- orbital weight is spread over m/n = 0.97
        of them (measured on two circuits) and does not follow `|lambda|`
        (Traversa & Bonani, TCAS-I 2011, sec. V) -- so `nmodes=None` is
        refused here, and a spectrum above the limit is `pnoise`'s.  What
        this serves is stability and mode inspection: the multipliers,
        exponents and shapes of the slow modes.  (Built 2026-09-25; until
        then refused outright.)"""
        if nmodes is None:
            raise NotImplementedError(
                'PSS.floquet_modes: this monodromy is %d wide, past '
                'FLOQUET_DENSE_LIMIT = %d, and ALL its modes are the dense '
                "route's. Pass nmodes=k for the k dominant modes (Ritz-"
                'certified). The modal spectra (orbital_correlation, '
                'orbital_spectrum, modal_spectrum) need every mode -- the '
                'orbital weight does not concentrate (m/n = 0.97 measured) -- '
                'so above the limit use PAC.pnoise for the spectrum.'
                % (n, self.FLOQUET_DENSE_LIMIT))
        k = int(nmodes)
        kmax = int(min(n, self.PPV_RITZ_MAX_BASIS)) // 2
        if not 1 <= k <= kmax:
            raise NotImplementedError(
                'PSS.floquet_modes: nmodes=%d on a %d-wide monodromy (past '
                'FLOQUET_DENSE_LIMIT = %d): an Arnoldi of at most %d vectors '
                '(PPV_RITZ_MAX_BASIS) certifies at most %d dominant modes. '
                'The modal spectra need every mode; above the limit use '
                'PAC.pnoise for the spectrum.'
                % (k, n, self.FLOQUET_DENSE_LIMIT, int(min(n, self.PPV_RITZ_MAX_BASIS)), kmax))
        mv, mvT = fp.matvec, fp.matvec_transposed
        _ev = EventColumns.of(self, n)
        if _ev is not None:
            ## the TOTAL map, as the dense path's `total_matrix`
            mv = _ev.total_matvec(fp.matvec)
            mvT = _ev.total_matvec(fp.matvec_transposed, transposed=True)
        lam, U, res_r, kk_r = self._ritz_modes(mv, n, k)
        ## the left set with a margin, so every right mode finds its partner
        lam_l, V, res_l, kk_l = self._ritz_modes(mvT, n, k + 2)
        certified = [bool(r <= self.PPV_RITZ_RESIDUAL_TOL) for r in res_r]
        if not all(certified):
            warnings.warn(
                'PSS.floquet_modes: %d of %d dominant modes did not certify '
                'within an Arnoldi of %d vectors (PPV_RITZ_MAX_BASIS; Ritz '
                'residuals %s against PPV_RITZ_RESIDUAL_TOL = %.0e). They are '
                "returned flagged 'certified': False -- read their multipliers "
                'and shapes as estimates.'
                % (certified.count(False), len(certified), kk_r,
                   ', '.join('%.1e' % r for r in res_r), self.PPV_RITZ_RESIDUAL_TOL),
                RuntimeWarning, stacklevel=3)
        warnings.warn(
            'PSS.floquet_modes: nmodes=%d returns the %d DOMINANT modes (by '
            '|lambda|) of a %d-wide monodromy. Orbital-noise weight does NOT '
            'follow multiplier magnitude -- Traversa & Bonani (TCAS-I 2011, '
            'Sec. V) show the contribution ordering INVERTING across six '
            'orders in mu, and this repo measured no concentration (m/n = '
            '0.97). These modes are for stability and inspection; a spectrum '
            'above the limit is PAC.pnoise\'s.' % (k, len(lam), n),
            RuntimeWarning, stacklevel=3)
        times = np.asarray(fp.times, dtype=float)
        tol = max(self.PPV_RITZ_RESIDUAL_TOL * 10.0,
                  10.0 * max(max(res_r), max(res_l)))
        used = set()
        out = []
        for i in range(len(lam)):
            d = np.abs(lam_l - lam[i])
            j = next((int(jj) for jj in np.argsort(d) if int(jj) not in used), None)
            if j is None or d[j] > tol * max(1.0, abs(lam[i])):
                raise ValueError(
                    'PSS.floquet_modes: the dominant multiplier %s found no '
                    'partner among the transposed map\'s Ritz values (nearest '
                    '%s) within %.1e -- the left and right Arnoldi runs did '
                    'not converge to the same eigenvalue. Ask for fewer modes.'
                    % (complex(lam[i]), None if j is None else complex(lam_l[j]), tol))
            used.add(j)
            md = self._floquet_mode(
                fp, i, complex(lam[i]), U[:, i].astype(complex),
                V[:, j].astype(complex), T, times, n,
                lambda uk, lk: float(np.linalg.norm(np.asarray(mv(uk)) - lk * uk)
                                     / max(abs(lk), 1e-300)))
            md['certified'] = certified[i]
            md['ritz_residual'] = float(res_r[i])
            out.append(md)
        return out

    def _ritz_modes(self, mv, n, k):
        """`(lams, X, residuals, kk)`: the `k` dominant eigenpairs (by
        `|lambda|`, a conjugate pair completed, null multipliers dropped) of
        the operator `mv` on `R^n`, from an Arnoldi grown from
        `PPV_RITZ_BASIS` by doubling up to `PPV_RITZ_MAX_BASIS` until every
        selected pair's Ritz residual ``|h_{kk+1,kk}| |y_last| / ||y||``
        (``||M x - theta x||`` for the unit Ritz vector) is within
        `PPV_RITZ_RESIDUAL_TOL` -- the gate `ppv`'s truncated `lam2` uses
        (`_ritz_second_multiplier`).

        ⚠ ON THE MAP ITSELF, NOT ON `I - M`: Arnoldi resolves the OUTER
        spectrum first, and the dominant modes are the outer ones; on `I -
        M` they are the smallest `theta` and resolve last (Garcia, Romero &
        Acha 2022).  The basis is EXTENDED, not restarted, and fully
        reorthogonalised (twice), from `ppv`'s seeded start.  An invariant
        subspace (the basis closes on itself) makes every pair exact."""
        budget = int(min(n, self.PPV_RITZ_MAX_BASIS))
        kk = int(min(n, max(self.PPV_RITZ_BASIS, 2 * k + 2), budget))
        rng = np.random.default_rng(12345)
        q0 = rng.standard_normal(n)
        Q = np.zeros((n, budget + 1))
        Q[:, 0] = q0 / np.linalg.norm(q0)
        H = np.zeros((budget + 1, budget))
        j = 0
        closed = False
        while True:
            while j < kk and not closed:
                w = np.asarray(mv(Q[:, j]), dtype=float).ravel()
                for _pass in range(2):
                    hcol = Q[:, :j + 1].T @ w
                    w = w - Q[:, :j + 1] @ hcol
                    H[:j + 1, j] += hcol
                H[j + 1, j] = float(np.linalg.norm(w))
                if H[j + 1, j] < 1e-13 * max(1.0, float(np.max(np.abs(H[:j + 2, :j + 1])))):
                    closed = True
                    kk = j + 1
                    break
                Q[:, j + 1] = w / H[j + 1, j]
                j += 1
            theta, Y = np.linalg.eig(H[:kk, :kk])
            order = [int(i) for i in np.argsort(-np.abs(theta))
                     if abs(theta[i]) > self.FLOQUET_NULL_TOL]
            sel = order[:k]
            ## a conjugate pair is completed, never split
            for i in list(sel):
                if abs(theta[i].imag) > 1e-12 * max(1.0, abs(theta[i])):
                    c = min(order, key=lambda jj: abs(theta[jj] - np.conj(theta[i])))
                    if c not in sel:
                        sel.append(c)
            ynorm = np.linalg.norm(Y[:, sel], axis=0)
            res = (np.zeros(len(sel)) if closed else
                   abs(H[kk, kk - 1]) * np.abs(Y[kk - 1, sel]) / np.maximum(ynorm, 1e-300))
            if closed or np.all(res <= self.PPV_RITZ_RESIDUAL_TOL) or kk >= budget:
                X = Q[:, :kk] @ (Y[:, sel] / ynorm[None, :])
                return theta[sel], X, [float(r) for r in res], kk
            kk = int(min(2 * kk, budget))

    def _floquet_mode(self, fp, k, lk, uk, vk, T, times, n, residual):
        """One Floquet mode's dict from its multiplier `lk` and its right and
        left eigenvectors `uk`, `vk` at `t = 0` (the map's width `n`):
        biorthonormalised, its periodic parts `p` (a forward replay) and `q`
        (the transposed replay) sampled on `times`.  `residual(uk, lk)` is
        its relative eigen-residual -- the dense matrix's, or the map's
        mat-vec above `FLOQUET_DENSE_LIMIT`; `k` names the mode in a refusal.
        (The loop body of `floquet_modes` until 2026-09-25, moved verbatim.)"""
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
        if hasattr(fp, 'forward_states'):
            ## a GLM on the state (`_GLMStateMap`): its own forward pass
            _end, fwd = fp.forward_states(uk)
        else:
            _end, fwd = self._forced_replay(fp, 0.0, zero, y0=uk,
                                            collect=True)
        traj = ([np.asarray(uk, dtype=complex)[:m]]
                + [np.asarray(z, dtype=complex).ravel()[:m] for z in fwd])
        tt = times[:len(traj)]
        traj = traj[:len(tt)]
        p = np.column_stack([traj[j] * np.exp(-muk * tt[j])
                             for j in range(len(traj))])

        ## adjoint: Phi(T,s_j)^T v_k(T), available under every integrator
        ## (B8).
        ## ⚠⚠ UNDER GEAR (`solved_history`) THE ADJOINT IS THE PER-STEP
        ## TRANSPOSED SOLVE `t`, NOT THE PAIR'S FIRST BLOCK, AND IT BELONGS
        ## TO THE NEXT NODE.  `collect` returns `ts[k]`, the solve
        ## `Jf_k^-T w1` made while replaying step k backwards from the pair
        ## at node k + 1 -- so the adjoint at node k + 1, second order --
        ## and `states[k]`, the pair (w1; w2).  With
        ## `w1 = (a0 C + G)^T t`, the first block through pinv(C^T) is
        ## `a0 * q(t + 2h/3)`: staggered by a fraction of a step, FIRST
        ## order.  `t` carries the step through `a0 ~ 1/h`, so the node's
        ## own `a0` is put back (on a uniform grid `c0` below absorbs it).
        ## On a non-uniform grid the transpose of a variable-step multistep
        ## method is the continuous adjoint's only to O(h) whatever the
        ## scaling (Sandu's inconsistency), so that case is integrated
        ## separately (below).  Node 0 is node N by periodicity of the
        ## periodic part.
        ## ⚠ GEAR ONLY.  A one-step kind's `ts` is NESTED (per-stage solves
        ## per step); its state-block adjoint through pinv(C^T) is exact
        ## (radau) and second order (trap), and keeps its path.
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
        ## the state monodromy `C^T w(0)` -- one factor of `C^T` away from
        ## the state-space adjoint `q` that eq (22) and every covariance
        ## here need.  The transposed replay propagates that object, so
        ## every sample of `q` above is `C(t)^T q_true(t)`.  ⚠ A symmetric
        ## orbit (van der Pol, reduced `C = diag(1, -1)`) hides the
        ## difference, and can make its sign flip read as a time
        ## reversal; an asymmetric orbit shows it.
        ##
        ## ⚠ Per sample, because `C` may depend on the state.  `pinv`
        ## rather than `inv` so a singular reduced `C` (an index-2 MNA,
        ## algebraic rows) does not raise; the algebraic components of `q`
        ## are then its minimum-norm 0, and are filled from the constraint
        ## below.
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
            ## ⚠⚠ AND THE ALGEBRAIC ENTRIES ARE SLAVED, NOT ZERO
            ## (2026-09-26).  `pinv` leaves them at its minimum-norm 0, so a
            ## source on an algebraic node reached no mode: `modal_spectrum`
            ## read EXACTLY 0 for one on radau, silently (gear's transposed
            ## solve carries them).  An algebraic state's column of the
            ## adjoint equation holds no derivative, whatever the mode's
            ## exponent, so the PPV's constraint fill applies as it stands
            ## (`_algebraic_adjoint_fill`), at each sample's own state.
            _irn = self.irefnode
            _xf0 = np.insert(_Wq[:, 0], _irn, 0.0)
            _arows, _acols = self._algebraic_adjoint_pattern(_xf0)
            if _arows:
                with warnings.catch_warnings(record=True) as _caught:
                    warnings.simplefilter('always')
                    for _j in range(q.shape[1]):
                        q[:, _j] = self._algebraic_adjoint_fill(
                            q[:, _j], np.insert(_Wq[:, min(_j, _nw - 1)],
                                                _irn, 0.0),
                            _arows, _acols)
                for _msg in sorted({str(_w.message) for _w in _caught}):
                    warnings.warn(_msg.replace('PSS.ppv', 'PSS.floquet_modes'),
                                  RuntimeWarning, stacklevel=3)

        ## ⚠⚠ RENORMALISE ON THE STATE BLOCK, WITH THE `C`-WEIGHTED INNER
        ## PRODUCT.  `v_k` was biorthonormalised against `u_k` at the map's
        ## FULL width `n`; under a solved-history map that is the pair
        ## `[x_n; x_{n-1}]`, and the width-`m` state block then carries
        ## `q(0)^T p(0) = c0 != 1` (on the plain path `n = m`).  And the
        ## conserved bilinear form of the variational DAE is
        ## `q(t)^T C(t) p(t)`, not `q(t)^T p(t)`: differentiating
        ## `G p + d(C p)/dt = 0` against the adjoint gives
        ## `d/dt [q^T C p] = 0`, so the biorthonormality eq (22) assumes is
        ## `q_k^T C p_l = d_kl`.  Normalising on `q^T C p` fixes the slice
        ## scale and the weighting in one step.  ⚠ `q` enters the
        ## covariance QUADRATICALLY, so either error squares; the
        ## periodicity gate p(T) = p(0) is scale-free, and a unit-reactance
        ## fixture (`C` the identity) cannot tell the two inner products
        ## apart.  The adjoint takes the scale (the right vector is the
        ## physical direction).
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
        return {'lam': lk, 'mu': muk, 'u0': uk, 'v0': vk,
                'p': p, 'q': q, 'times': tt, 'c0': c0,
                'residual': residual(uk, lk)}

    def _continuous_adjoint(self, fp, lam, mu, times):
        """A mode's adjoint `q(t_j)` on a NON-UNIFORM gear grid, by integrating
        the continuous adjoint equation SEPARATELY.

        The exact transpose of gear's two-step recursion draws `a1` and `a2`
        from LATER steps, so on a grid whose step changes it is a consistent
        scheme for the adjoint equation only to first order, and no per-node
        rescaling lifts it.  Here `C^T dq/dt = G^T q` is integrated backwards
        with BDF2 whose coefficients belong to the REVERSE grid's own step
        pair -- `companion_coefficients(h_{n+1}, h_{n+2})`:

            (a0' C_n + G_n)^T q_n = -C_n^T (a1' q_{n+1} + a2' q_{n+2})

        The state is the pair (q_{n+1}, q_{n+2}); the backward map over one
        period is built from 2m basis propagations and the mode's eigenvector
        matched to the forward multiplier `lam` (they agree to O(h^2)).  The
        invariant `q^T C p` is then second order on a 3:1 grid.

        ⚠ TWO ADJOINTS IN THE TREE, ON PURPOSE.  This one is NOT the transpose
        of the discrete map: the PPV's anchor, the adjoint noise folds and the
        sideband rows keep the exact transpose (pinned at 1e-15) because
        their identities need it.  This serves the MODE SHAPES
        (`floquet_modes`) and, at index 2 on a non-uniform solved-history
        grid, `ppv`'s samples (see the note there) -- where the continuous
        adjoint is the object wanted.  On a uniform grid the two
        coincide (`a0' = a0`) and the transpose is used; one-step kinds never
        come here.  Dense 2m x 2m up to `CONTINUOUS_ADJOINT_DENSE_M`
        unknowns; above that, Arnoldi on the backward map with Ritz-residual
        certification (see the note at the branch).  ⚠ Biorthonormality
        against the OTHER modes is not exact here (it is for the transpose)
        -- each mode is normalised on its own `q(0)^T C p(0)`.

        History: `doc/shooting_history.md`, `_PPVFloquet._continuous_adjoint`.
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
        coef, lus = [], []
        for n in range(N):
            a = [float(x) for x in
                 integ.companion_coefficients(h[n % N], h[(n + 1) % N])[0]]
            coef.append(a)
            ## one factorisation per node, reused by every propagation --
            ## by the caller's solver (`_factorise`; `lu_factor` under the
            ## default `DenseSolver`, as it was)
            lus.append(self._factorise((a[0] * Cs[n] + Gs[n]).T))

        def propagate(qN, qN1):
            q = [None] * (N + 2)
            q[N], q[N + 1] = qN, qN1
            for n in range(N - 1, -1, -1):
                _a0, a1, a2 = coef[n]
                q[n] = lus[n].solve(-Cs[n].T @ (a1 * q[n + 1] + a2 * q[n + 2]))
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
            ## ⚠ MATRIX-FREE ABOVE THAT.  The wanted modes -- the phase mode
            ## and the slow orbital ones -- are the DOMINANT eigenvalues of the
            ## backward map, so a plain Arnoldi on it with the Ritz-residual
            ## certification `_ritz_second_multiplier` uses
            ## (|h_{k+1,k}| |y_last| / ||y||) reaches them at a basis far below
            ## 2m.  Grown from `PPV_RITZ_BASIS` toward `PPV_RITZ_MAX_BASIS`
            ## until the matched pair certifies; a pair that never certifies
            ## is refused, not returned.
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
